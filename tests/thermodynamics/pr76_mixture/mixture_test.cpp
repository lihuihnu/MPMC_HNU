#include <mpmc/thermodynamics/pr76_mixture.hpp>
#include <mpmc/ad/differentiate.hpp>
#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <source_location>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

double mixture_plain_header(const mpmc::thermodynamics::Pr76Mixture<double>& mixture);

namespace {
namespace th = mpmc::thermodynamics;
namespace ad = mpmc::ad;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}
template <typename Error, typename Function>
void expect_error(Function&& function) {
    bool caught = false;
    try { function(); } catch (const Error&) { caught = true; }
    require(caught, "expected exception was not thrown");
}
template <std::floating_point T>
void near(T actual, long double expected,
          std::source_location where = std::source_location::current()) {
    // New sums/products and derivative transformations allow 512 rounding units.
    // Relative only: a spurious zero must not pass a small nonzero derivative test.
    const long double error = std::abs(static_cast<long double>(actual) - expected);
    const long double limit = 512.0L * std::numeric_limits<T>::epsilon() * std::abs(expected);
    if (!std::isfinite(actual) || !std::isfinite(expected) || error > limit) {
        std::ostringstream message;
        message << std::setprecision(24) << "mixture mismatch: actual=" << actual
                << ", expected=" << expected << ", absolute error=" << error;
        require(false, message.str(), where);
    }
}

// Explicitly manufactured parameter records, not properties of real substances.
struct Spec { double tc, pc, omega; };
constexpr std::array<Spec, 5> specs{{{400, 4e6, -0.125}, {500, 3e6, 0.25},
                                    {600, 5e6, 0.5}, {350, 2.5e6, 0.125},
                                    {450, 6e6, 0.375}}};
constexpr double interactions[5][5] = {
    {0, .125, -.0625, .03125, .25}, {.125, 0, .0625, -.125, .015625},
    {-.0625, .0625, 0, .25, -.03125}, {.03125, -.125, .25, 0, .125},
    {.25, .015625, -.03125, .125, 0}};
std::string id(std::size_t i) { return "test:mix:" + std::to_string(i); }
th::Provenance source(std::string locator) {
    return {th::SourceKind::synthetic_test, "test://pr76-mixture", "v1", std::move(locator),
            "Artificial formula/derivative fixture, not fluid data", "Constructed in tests",
            "No third-party property data copied"};
}
th::SourcedScalar datum(double value, th::Unit unit) {
    return {value, unit, source("parameter"), "specified molar SI unit", "identity"};
}
struct Fixture {
    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    Fixture() {
        input.model_id = th::pr76_profile;
        input.dataset_id = "manufactured-mixture";
        input.revision = "v1";
        input.applicability = {std::nullopt, std::nullopt, source("unknown bounds")};
        for (std::size_t i = 0; i < specs.size(); ++i) {
            catalog.push_back({id(i), "Artificial", th::ComponentKind::pure, source(id(i)), {}});
            input.pure.push_back({id(i), datum(specs[i].tc, th::Unit::kelvin),
                                  datum(specs[i].pc, th::Unit::pascal),
                                  datum(specs[i].omega, th::Unit::dimensionless)});
            for (std::size_t j = 0; j < i; ++j) {
                // Intentionally reverse the ID order relative to requested snapshots.
                input.binary.push_back({id(i), id(j),
                                        datum(interactions[i][j], th::Unit::dimensionless)});
            }
        }
    }
    th::PrParameterSet select(const std::vector<std::size_t>& indices) const {
        std::vector<std::string> order;
        for (const auto i : indices) { order.push_back(id(i)); }
        return th::PrParameterSet::create(catalog, order, input,
                                        th::DataPolicy::allow_synthetic_tests);
    }
    void set_spec(std::size_t i, const Spec& spec) {
        input.pure[i].critical_temperature->value = spec.tc;
        input.pure[i].critical_pressure->value = spec.pc;
        input.pure[i].acentric_factor->value = spec.omega;
    }
};
template <typename T>
th::Pr76Mixture<T> kernel(const std::vector<std::size_t>& order = {0, 1, 2}) {
    return th::Pr76Mixture<T>::from_parameters(Fixture{}.select(order));
}

struct PureReference { long double a, b, da; };
PureReference pure_reference(long double temperature, Spec spec) {
    // Original literals and raw fixture data; no production functions or coefficients.
    const long double r = 8.31446261815324L, tc = spec.tc, pc = spec.pc, w = spec.omega;
    const long double ac = 0.45724L * r * r * tc * tc / pc;
    const long double k = 0.37464L + 1.54226L*w - 0.26992L*w*w;
    const long double q = 1 + k*(1 - std::sqrt(temperature/tc));
    return {ac*q*q, 0.07780L*r*tc/pc, -ac*k*q/std::sqrt(temperature*tc)};
}
struct Reference {
    long double a = 0, b = 0, da = 0;
    std::vector<long double> ax, bx;
};
template <typename T>
Reference reference(T temperature, const std::vector<T>& x,
                    const std::vector<std::size_t>& order) {
    const std::size_t n = order.size();
    Reference result;
    result.ax.resize(n); result.bx.resize(n);
    std::vector<PureReference> pure;
    for (const auto i : order) { pure.push_back(pure_reference(temperature, specs[i])); }
    // Full double sum as printed, deliberately distinct from production's triangle.
    for (std::size_t i = 0; i < n; ++i) {
        result.b += static_cast<long double>(x[i]) * pure[i].b;
        result.bx[i] = pure[i].b;
        for (std::size_t j = 0; j < n; ++j) {
            const long double aij = (1 - static_cast<long double>(interactions[order[i]][order[j]]))
                                    * std::sqrt(pure[i].a) * std::sqrt(pure[j].a);
            const long double weight = static_cast<long double>(x[i]) * x[j];
            result.a += weight * aij;
            result.da += weight * aij * (pure[i].da/pure[i].a + pure[j].da/pure[j].a) / 2;
            result.ax[i] += 2 * static_cast<long double>(x[j]) * aij;
        }
    }
    return result;
}

template <typename T>
void printed_binary() {
    Fixture fixture;
    fixture.set_spec(0, {400, 4e6, 0});
    fixture.set_spec(1, {400, 1e6, 0}); // a1=4*a0, b1=4*b0, same alpha.
    const auto mix = th::Pr76Mixture<T>::from_parameters(fixture.select({0, 1}));
    th::Pr76MixtureWorkspace<T> workspace;
    const std::array<T, 2> x{T{0.25}, T{0.75}};
    const auto value = mix.evaluate_full(T{100}, x, workspace);
    // x=(1/4,3/4), kij=1/8: a_mix=(95/32)*a0, b_mix=(13/4)*b0.
    const auto expected = pure_reference(100, {400, 4e6, 0});
    near(value.a, expected.a * 95 / 32);
    near(value.b, expected.b * 13 / 4);
    // Independent rational mixture multipliers test the cross term's factor of two.
    using Number = ad::Dual<T, 1>;
    th::Pr76MixtureWorkspace<Number> ad_workspace;
    const std::array<Number, 2> ax{Number{x[0]}, Number{x[1]}};
    const auto d = mix.evaluate_full(Number::variable(T{100}, 0), ax, ad_workspace);
    near(d.a.derivative(0), expected.da * 95 / 32);
    require(d.b.derivative(0) == T{0}, "b must not depend on T at fixed x");
}

template <typename T>
void full_gradients() {
    const auto mix = kernel<T>();
    const std::vector<T> x{T{0.25}, T{0.375}, T{0.375}};
    const auto ref = reference(T{300}, x, {0, 1, 2});
    using Number = ad::Dual<T, 4>;
    th::Pr76MixtureWorkspace<Number> work;
    std::vector<Number> ax;
    for (std::size_t i = 0; i < x.size(); ++i) {
        ax.push_back(Number::variable(x[i], i + 1));
    }
    const auto result = mix.evaluate_full(Number::variable(T{300}, 0), ax, work);
    near(result.a.value(), ref.a); near(result.b.value(), ref.b);
    near(result.a.derivative(0), ref.da);
    require(result.b.derivative(0) == T{0}, "full b temperature column");
    for (std::size_t i = 0; i < x.size(); ++i) {
        near(result.a.derivative(i + 1), ref.ax[i]);
        near(result.b.derivative(i + 1), ref.bx[i]);
    }
}

template <typename T>
void reduced_gradients() {
    const auto mix = kernel<T>();
    const auto ref = reference(T{300}, std::vector<T>{T{0.25}, T{0.375}, T{0.375}}, {0, 1, 2});
    const auto result = ad::value_and_jacobian([&](const auto& p) {
        using Number = typename std::remove_cvref_t<decltype(p)>::value_type;
        th::Pr76MixtureWorkspace<Number> work;
        const auto value = mix.evaluate_reduced(p[0], std::span<const Number>{p}.subspan(1), work);
        return std::array{value.a, value.b};
    }, std::array<T, 3>{T{300}, T{0.25}, T{0.375}});
    near(result.values[0], ref.a); near(result.values[1], ref.b);
    near(result.jacobian[0][0], ref.da);
    require(result.jacobian[1][0] == T{0}, "reduced b temperature column");
    for (std::size_t j = 0; j < 2; ++j) {
        near(result.jacobian[0][j+1], ref.ax[j] - ref.ax[2]);
        near(result.jacobian[1][j+1], ref.bx[j] - ref.bx[2]);
    }
}

template <typename T>
void chain_rule() {
    const auto mix = kernel<T>();
    const auto ref = reference(T{300}, std::vector<T>{T{0.25}, T{0.375}, T{0.375}}, {0, 1, 2});
    using Number = ad::Dual<T, 3>;
    const Number t{T{300}, {T{2}, T{-1}, T{0}}};
    const std::array<Number, 2> x{
        Number{T{0.25}, {T{1}, T{3}, T{0}}}, Number{T{0.375}, {T{-2}, T{1}, T{0}}}};
    th::Pr76MixtureWorkspace<Number> work;
    const auto result = mix.evaluate_reduced(t, x, work);
    const long double ga = ref.ax[0]-ref.ax[2], ha = ref.ax[1]-ref.ax[2];
    const long double gb = ref.bx[0]-ref.bx[2], hb = ref.bx[1]-ref.bx[2];
    near(result.a.derivative(0), 2*ref.da + ga - 2*ha);
    near(result.a.derivative(1), -ref.da + 3*ga + ha);
    near(result.b.derivative(0), gb - 2*hb);
    near(result.b.derivative(1), 3*gb + hb);
    require(result.a.derivative(2) == T{0} && result.b.derivative(2) == T{0}, "zero lane");
    const std::array<Number, 2> constants{Number{T{0.25}}, Number{T{0.375}}};
    const auto zero = mix.evaluate_reduced(Number{T{300}}, constants, work);
    require(zero.a.derivatives() == typename Number::Gradient{}, "zero seeds not preserved");
}

template <typename T>
void zero_fractions() {
    const auto mix = kernel<T>();
    const std::vector<T> x{T{1}, T{0}, T{0}};
    const auto ref = reference(T{300}, x, {0, 1, 2});
    using Number = ad::Dual<T, 3>;
    const std::array<Number, 3> ax{Number::variable(x[0], 0), Number::variable(x[1], 1),
                                    Number::variable(x[2], 2)};
    th::Pr76MixtureWorkspace<Number> work;
    const auto result = mix.evaluate_full(Number{T{300}}, ax, work);
    near(result.a.value(), ref.a);
    for (std::size_t i = 0; i < 3; ++i) {
        near(result.a.derivative(i), ref.ax[i]); near(result.b.derivative(i), ref.bx[i]);
    }
    require(result.a.derivative(1) != T{0}, "zero-fraction insertion derivative lost");
    const std::array<Number, 2> reduced{ax[0], ax[1]};
    const auto constrained = mix.evaluate_reduced(Number{T{300}}, reduced, work);
    near(constrained.a.derivative(1), ref.ax[1]-ref.ax[2]);
    // This boundary result is the polynomial continuation/tangent derivative,
    // not a claim that every signed perturbation remains in the physical simplex.
}

template <typename T>
void permutations() {
    std::vector<std::size_t> order{0, 1, 2, 3};
    const std::array<T, 4> by_id{T{0.125}, T{0.25}, T{0.125}, T{0.5}};
    const auto original = reference(T{300}, std::vector<T>{by_id.begin(), by_id.end()}, order);
    std::size_t count = 0;
    do {
        const auto mix = kernel<T>(order);
        std::vector<T> x;
        for (const auto i : order) { x.push_back(by_id[i]); }
        const auto ref = reference(T{300}, x, order);
        using Number = ad::Dual<T, 5>;
        std::vector<Number> ax;
        for (std::size_t i = 0; i < 4; ++i) { ax.push_back(Number::variable(x[i], i+1)); }
        th::Pr76MixtureWorkspace<Number> work;
        const auto full = mix.evaluate_full(Number::variable(T{300}, 0), ax, work);
        near(full.a.value(), original.a); near(full.b.value(), original.b);
        near(full.a.derivative(0), original.da);
        const auto reduced = mix.evaluate_reduced(Number::variable(T{300}, 0),
                                                  std::span<const Number>{ax}.first(3), work);
        near(reduced.a.value(), original.a); near(reduced.b.value(), original.b);
        for (std::size_t i = 0; i < 4; ++i) {
            near(full.a.derivative(i+1), original.ax[order[i]]);
            near(full.b.derivative(i+1), original.bx[order[i]]);
            require(mix.parameters().components().at(i).id == id(order[i]), "model ID order");
            if (i < 3) {
                near(reduced.a.derivative(i+1), ref.ax[i]-ref.ax[3]);
                near(reduced.b.derivative(i+1), ref.bx[i]-ref.bx[3]);
            }
        }
        ++count;
    } while (std::next_permutation(order.begin(), order.end()));
    require(count == 24, "not all four-component permutations checked");
}

template <typename T>
void runtime_shapes() {
    using Number = ad::Dual<T, 2>;
    ad::RuntimeJacobianWorkspace<T, 2> jacobian_work;
    th::Pr76MixtureWorkspace<Number> mixing_work;
    // Same workspaces, changing n; all components share ONE temperature.
    for (const std::size_t n : std::array<std::size_t, 5>{1, 5, 2, 4, 1}) {
        std::vector<std::size_t> order(n);
        std::iota(order.begin(), order.end(), std::size_t{0});
        const auto mix = kernel<T>(order);
        std::vector<T> inputs(n, T{0.125}); inputs[0] = T{300};
        std::vector<T> x(n, T{0.125}); x.back() = T{1}-T{0.125}*static_cast<T>(n-1);
        const auto ref = reference(T{300}, x, order);
        std::size_t calls = 0;
        const auto result = ad::value_and_jacobian_runtime<2>(
            [&](std::span<const Number> p, std::span<Number> y) {
                ++calls;
                const auto v = mix.evaluate_reduced(p[0], p.subspan(1), mixing_work);
                y[0] = v.a; y[1] = v.b;
            }, inputs, 2, jacobian_work);
        require(calls == (n+1)/2 && result.input_count == n && result.output_count == 2,
                "runtime shape or block count");
        near(result.values[0], ref.a); near(result.values[1], ref.b);
        near(result.jacobian[0], ref.da);
        require(result.jacobian[n] == T{0}, "runtime db/dT");
        for (std::size_t j = 1; j < n; ++j) {
            near(result.jacobian[j], ref.ax[j-1]-ref.ax.back());
            near(result.jacobian[n+j], ref.bx[j-1]-ref.bx.back());
        }
    }
    require(mixing_work.capacity() >= 5, "workspace should retain capacity after shrink");
}

template <typename T>
void normalization() {
    Fixture fixture;
    fixture.set_spec(1, specs[0]); fixture.input.binary[0].kij->value = 0;
    const auto mix = th::Pr76Mixture<T>::from_parameters(fixture.select({0, 1}));
    th::Pr76MixtureWorkspace<T> work;
    const T eps = std::numeric_limits<T>::epsilon();
    const std::array<T, 2> close{T{0.5}, T{0.5}+T{16}*eps};
    const auto value = mix.evaluate_full(T{300}, close, work);
    const auto pure = th::Pr76Pure<T>::from_parameters(mix.parameters(), 0).evaluate(T{300});
    // If the kernel silently normalized x, both differences would disappear.
    const T sum = close[0] + close[1];
    near(value.a, static_cast<long double>(pure.a)*sum*sum);
    near(value.b, static_cast<long double>(pure.b)*sum);
    require(value.b > pure.b, "accepted near-normalized input was silently normalized");
    const std::array<T, 2> bad{T{0.5}, T{0.5}+T{128}*eps};
    expect_error<std::domain_error>([&] { (void)mix.evaluate_full(T{300}, bad, work); });
    for (T v : std::array<T, 6>{T{-0.01}, T{1.01},
            std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::infinity(),
            -std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::min()}) {
        const std::array<T, 2> invalid{v, T{0.5}};
        expect_error<std::domain_error>([&] { (void)mix.evaluate_full(T{300}, invalid, work); });
    }
    expect_error<std::invalid_argument>([&] { (void)mix.evaluate_full(T{300}, {}, work); });
    expect_error<std::invalid_argument>([&] { (void)mix.evaluate_reduced(T{300}, close, work); });
    const auto three = kernel<T>();
    const std::array<T, 2> overshoot{T{0.5}, T{0.5}+eps};
    expect_error<std::domain_error>([&] { (void)three.evaluate_reduced(T{300}, overshoot, work); });
    const std::array<T, 2> zeros{T{0}, T{0}};
    expect_error<std::domain_error>([&] { (void)mix.evaluate_full(T{300}, zeros, work); });
}

template <typename T>
void temperature_and_seeds() {
    Fixture fixture;
    fixture.input.applicability.temperature_k = th::ClosedInterval{200, 600};
    const auto mix = th::Pr76Mixture<T>::from_parameters(fixture.select({0, 1}));
    th::Pr76MixtureWorkspace<T> work;
    const std::array<T, 1> x{T{0.5}};
    (void)mix.evaluate_reduced(T{200}, x, work); (void)mix.evaluate_reduced(T{600}, x, work);
    for (T t : std::array<T, 6>{T{0}, T{-1}, T{199}, T{601},
            std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::infinity()}) {
        expect_error<std::domain_error>([&] { (void)mix.evaluate_reduced(t, x, work); });
    }
    using Number = ad::Dual<T, 1>;
    th::Pr76MixtureWorkspace<Number> ad_work;
    const T inf = std::numeric_limits<T>::infinity();
    const std::array<Number, 1> bad{Number{T{0.5}, {inf}}};
    const std::array<Number, 1> good{Number{T{0.5}}};
    expect_error<std::domain_error>([&] {
        (void)mix.evaluate_reduced(Number{T{300}}, bad, ad_work);
    });
    expect_error<std::domain_error>([&] {
        (void)mix.evaluate_reduced(Number{T{300}, {inf}}, good, ad_work);
    });
}

template <typename T>
void pair_semantics() {
    for (const double kij : std::array{0.0, -0.25, 1.0, 3.0}) {
        Fixture f; f.set_spec(1, specs[0]); f.input.binary[0].kij->value = kij;
        const auto mix = th::Pr76Mixture<T>::from_parameters(f.select({0, 1}));
        th::Pr76MixtureWorkspace<T> work;
        const auto result = mix.evaluate_reduced(T{300}, std::array<T, 1>{T{0.5}}, work);
        const auto pure = pure_reference(300, specs[0]);
        near(result.a, pure.a*(1-static_cast<long double>(kij)/2));
        near(result.b, pure.b);
        // In particular kij=3 yields NEGATIVE a_mix: no clipping or validity claim.
    }
    Fixture missing; missing.input.binary.erase(missing.input.binary.begin());
    expect_error<th::ContractError>([&] {
        (void)th::Pr76Mixture<T>::from_parameters(missing.select({0, 1}));
    });
    Fixture f;
    expect_error<std::length_error>([&] {
        (void)th::Pr76Mixture<T>::from_parameters(f.select({0, 1}), {1, 10, 25});
    });
    if constexpr (std::numeric_limits<T>::max_exponent <
                  std::numeric_limits<double>::max_exponent) {
        f.input.binary[0].kij->value = std::numeric_limits<double>::max();
        expect_error<std::range_error>([&] {
            (void)th::Pr76Mixture<T>::from_parameters(f.select({0, 1}));
        });
    }
}

template <typename T>
void principal_roots() {
    Fixture f; f.set_spec(0, {100, 4e6, .5}); f.set_spec(1, {1600, 3e6, 0});
    const auto mix = th::Pr76Mixture<T>::from_parameters(f.select({0, 1}));
    using Number = ad::Dual<T, 1>;
    th::Pr76MixtureWorkspace<Number> work;
    const std::array<Number, 1> x{Number{T{0.5}}};
    const auto result = mix.evaluate_reduced(Number::variable(T{1600}, 0), x, work);
    const auto a = pure_reference(1600, {100, 4e6, .5});
    const auto b = pure_reference(1600, {1600, 3e6, 0});
    const long double cross = .875L*std::sqrt(a.a)*std::sqrt(b.a);
    near(result.a.value(), (a.a+b.a+2*cross)/4);
    near(result.a.derivative(0), (a.da+b.da+cross*(a.da/a.a+b.da/b.a))/4);
    // First q is negative; using signed q_i*q_j instead of principal roots is wrong.
}

void zero_attraction_float() {
    Fixture f; f.set_spec(0, {400, 4e6, .125});
    // This binary32 temperature makes the documented PR76 bracket round to zero
    // with the existing scalar evaluation order; not an empirical fluid state.
    constexpr float temperature = 3081.476806640625F;
    const auto pure = th::Pr76Pure<float>::from_parameters(f.select({0}), 0);
    require(pure.evaluate(temperature).a == 0.0F, "zero-bracket boundary fixture changed");
    th::Pr76MixtureWorkspace<float> work;
    const auto single = th::Pr76Mixture<float>::from_parameters(f.select({0}));
    require(single.evaluate_reduced(temperature, {}, work).a == 0.0F, "pure reduction at ai=0");
    const auto multi = th::Pr76Mixture<float>::from_parameters(f.select({0, 1}));
    expect_error<std::domain_error>([&] {
        (void)multi.evaluate_reduced(temperature, std::array<float, 1>{0.0F}, work);
    });
}

template <typename T>
void ownership_and_recovery() {
    auto mix = kernel<T>(); // Factory's source parameter set has already died.
    th::Pr76MixtureWorkspace<T> first_work, second_work;
    const std::array<T, 2> x{T{0.25}, T{0.375}};
    const auto first = mix.evaluate_reduced(T{300}, x, first_work);
    const auto copy = mix;
    expect_error<std::domain_error>([&] { (void)mix.evaluate_reduced(T{0}, x, first_work); });
    const auto recovered = copy.evaluate_reduced(T{300}, x, first_work);
    const auto independent = mix.evaluate_reduced(T{300}, x, second_work);
    require(first.a == recovered.a && first.b == recovered.b && independent.a == first.a,
            "failed call or independent workspace changed result");
    require(copy.parameters().dataset_id() == "manufactured-mixture", "source lifetime");
    const auto source_kind = copy.parameters().binary_records()[0].kij->source.kind;
    require(source_kind == th::SourceKind::synthetic_test, "binary provenance lost");
    const auto single = kernel<T>({2});
    const auto reduced = single.evaluate_reduced(T{300}, {}, first_work);
    near(reduced.a, pure_reference(300, specs[2]).a);
    near(first.a, reference(T{300}, std::vector<T>{T{0.25}, T{0.375}, T{0.375}}, {0, 1, 2}).a);
}

template <typename T>
void numerical_scale() {
    const std::array<double, 2> pressures = std::same_as<T, float>
        ? std::array{1e-20, 1e30} : std::array{1e-180, 1e180};
    for (const double pc : pressures) {
        Fixture f; f.set_spec(0, {400, pc, 0}); f.set_spec(1, {400, pc, 0});
        const auto mix = th::Pr76Mixture<T>::from_parameters(f.select({0, 1}));
        using Number = ad::Dual<T, 1>;
        th::Pr76MixtureWorkspace<Number> work;
        const auto value = mix.evaluate_reduced(Number::variable(T{300}, 0),
                                                std::array<Number, 1>{Number{T{0.5}}}, work);
        const auto ref = pure_reference(300, {400, pc, 0});
        near(value.a.value(), ref.a*0.9375L);
        near(value.a.derivative(0), ref.da*0.9375L);
        near(value.b.value(), ref.b);
        // Product ai*aj could overflow/underflow although the mixing result does not.
    }
}

template <typename Number>
concept AcceptsMixture = requires(const th::Pr76Mixture<double>& m, const Number& t,
                                  th::Pr76MixtureWorkspace<Number>& w) {
    m.evaluate_full(t, std::span<const Number>{}, w);
};
static_assert(AcceptsMixture<double> && AcceptsMixture<ad::Dual<double, 4>>);
static_assert(!AcceptsMixture<float> && !AcceptsMixture<int>);
static_assert(!std::is_default_constructible_v<th::Pr76Mixture<double>>);
static_assert(!std::is_copy_assignable_v<th::Pr76Mixture<double>>);
static_assert(!std::is_copy_constructible_v<th::Pr76MixtureWorkspace<double>>);
template <typename Kernel>
concept TemporaryParameters = requires(Kernel&& k) { std::move(k).parameters(); };
static_assert(!TemporaryParameters<th::Pr76Mixture<double>>);

void headers() {
    const auto mix = kernel<double>({0});
    near(mixture_plain_header(mix), pure_reference(300, specs[0]).a);
}

template <typename Function>
void all_precisions(Function function) {
    function.template operator()<float>();
    function.template operator()<double>();
    function.template operator()<long double>();
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 2) { std::cerr << "Expected a named mixture test\n"; return 2; }
    const std::string_view name = argv[1];
    try {
#define MIX_CASE(case_name) \
        if (name == #case_name) { all_precisions([]<typename T> { case_name<T>(); }); } else
        MIX_CASE(printed_binary)
        MIX_CASE(full_gradients)
        MIX_CASE(reduced_gradients)
        MIX_CASE(chain_rule)
        MIX_CASE(zero_fractions)
        MIX_CASE(permutations)
        MIX_CASE(runtime_shapes)
        MIX_CASE(normalization)
        MIX_CASE(temperature_and_seeds)
        MIX_CASE(pair_semantics)
        MIX_CASE(principal_roots)
        MIX_CASE(ownership_and_recovery)
        MIX_CASE(numerical_scale)
        if (name == "zero_attraction_float") { zero_attraction_float(); }
        else if (name == "headers") { headers(); }
        else { std::cerr << "Unknown test\n"; return 2; }
#undef MIX_CASE
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
