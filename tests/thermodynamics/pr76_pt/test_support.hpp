#ifndef MPMC_TEST_PR76_PT_SUPPORT_HPP
#define MPMC_TEST_PR76_PT_SUPPORT_HPP

#include <mpmc/thermodynamics/pr76_phase.hpp>
#include <mpmc/ad/math.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <numbers>
#include <optional>
#include <source_location>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace test {
namespace th = mpmc::thermodynamics;
namespace ad = mpmc::ad;
inline void require(bool condition, std::string_view message,
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
    require(caught, "expected error not reported");
}
template <std::floating_point T>
void near(T actual, long double expected,
          std::source_location where = std::source_location::current()) {
    // Root solving, logs, mixing and analytic derivative crosschecks: relative
    // 4096*epsilon; no absolute floor that could hide a spurious zero result.
    const long double error = std::abs(static_cast<long double>(actual) - expected);
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        error > 4096.0L * std::numeric_limits<T>::epsilon() * std::abs(expected)) {
        std::ostringstream out;
        out << std::setprecision(24) << "actual=" << actual << " expected=" << expected
            << " absolute error=" << error;
        require(false, out.str(), where);
    }
}

struct Spec { double tc, pc, omega; };
inline constexpr std::array<Spec, 4> specs{{{400, 4e6, -.125}, {500, 3e6, .25},
                                           {600, 5e6, .5}, {350, 2.5e6, .125}}};
inline constexpr double kij[4][4] = {{0,.125,-.0625,.03125}, {.125,0,.0625,-.125},
                                     {-.0625,.0625,0,.25}, {.03125,-.125,.25,0}};
inline std::string id(std::size_t i) { return "test:pt:" + std::to_string(i); }
inline th::Provenance source(std::string locator) {
    return {th::SourceKind::synthetic_test, "test://pr76-pt", "v1", std::move(locator),
            "Manufactured formula tests, NOT real fluids", "Constructed in tests",
            "No third-party parameter data copied"};
}
inline th::SourcedScalar datum(double value, th::Unit unit) {
    return {value, unit, source("parameter"), "specified molar SI", "identity"};
}
struct Fixture {
    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    Fixture() {
        input.model_id = th::pr76_profile;
        input.dataset_id = "manufactured-pt";
        input.revision = "v1";
        input.applicability = {std::nullopt, std::nullopt, source("unknown bounds")};
        for (std::size_t i = 0; i < specs.size(); ++i) {
            catalog.push_back({id(i), "Artificial", th::ComponentKind::pure, source(id(i)), {}});
            input.pure.push_back({id(i), datum(specs[i].tc, th::Unit::kelvin),
                                  datum(specs[i].pc, th::Unit::pascal),
                                  datum(specs[i].omega, th::Unit::dimensionless)});
            for (std::size_t j = 0; j < i; ++j) {
                input.binary.push_back({id(i), id(j), datum(kij[i][j], th::Unit::dimensionless)});
            }
        }
    }
    th::PrParameterSet select(const std::vector<std::size_t>& order) const {
        std::vector<std::string> names;
        for (const auto i : order) { names.push_back(id(i)); }
        return th::PrParameterSet::create(catalog, names, input, th::DataPolicy::allow_synthetic_tests);
    }
};
template <typename T>
th::Pr76Phase<T> kernel(const std::vector<std::size_t>& order = {0,1,2}) {
    return th::Pr76Phase<T>::from_parameters(Fixture{}.select(order));
}

// Independent reference uses the ORIGINAL Z polynomial and Cardano/trigonometry,
// not the production shifted-domain interval solver. Inputs here stay well clear
// of discriminant cancellation; dedicated manufactured tests cover degeneracy.
inline std::vector<long double> reference_roots(long double a, long double b) {
    const long double c2 = b-1, c1 = a-3*b*b-2*b, c0 = b*b*b+b*b-a*b;
    const long double q = (3*c1-c2*c2)/9;
    const long double r = (9*c2*c1-27*c0-2*c2*c2*c2)/54;
    const long double d = q*q*q+r*r;
    std::vector<long double> values;
    if (d >= 0) {
        values.push_back(std::cbrt(r+std::sqrt(d))+std::cbrt(r-std::sqrt(d))-c2/3);
    } else {
        const long double theta = std::acos(r/std::sqrt(-q*q*q));
        for (int k = 0; k < 3; ++k) {
            values.push_back(2*std::sqrt(-q)*std::cos((theta+2*k*std::numbers::pi_v<long double>)/3)-c2/3);
        }
    }
    std::erase_if(values, [b](long double z) { return z <= b; });
    std::sort(values.begin(), values.end());
    return values;
}
struct PureReference { long double a, b, da; };
inline PureReference pure_reference(long double temperature, Spec spec) {
    const long double r = 8.31446261815324L, tc = spec.tc, pc = spec.pc, omega = spec.omega;
    const long double ac = .45724L*r*r*tc*tc/pc;
    const long double kappa = .37464L+1.54226L*omega-.26992L*omega*omega;
    const long double q = 1+kappa*(1-std::sqrt(temperature/tc));
    return {ac*q*q, .07780L*r*tc/pc, -ac*kappa*q/std::sqrt(temperature*tc)};
}
struct ReferencePhase {
    long double z{};
    std::vector<long double> ln_phi, dz;
    std::vector<std::vector<long double>> dln; // [component][p,T,w_0,...].
};
template <std::floating_point T>
std::vector<ReferencePhase> reference(T pressure, T temperature, const std::vector<T>& w,
                                     const std::vector<std::size_t>& order) {
    const std::size_t n = w.size(), columns = n+2;
    const long double p = pressure, t = temperature, r = 8.31446261815324L;
    std::vector<PureReference> pure;
    for (const auto i : order) { pure.push_back(pure_reference(t, specs[i])); }
    std::vector<std::vector<long double>> aij(n, std::vector<long double>(n));
    std::vector<long double> s(n), st(n);
    long double a = 0, b = 0, at = 0;
    for (std::size_t i = 0; i < n; ++i) {
        b += static_cast<long double>(w[i])*pure[i].b;
        for (std::size_t j = 0; j < n; ++j) {
            aij[i][j] = (1-static_cast<long double>(kij[order[i]][order[j]]))*
                         std::sqrt(pure[i].a)*std::sqrt(pure[j].a);
            const long double derivative = aij[i][j]*(pure[i].da/pure[i].a+pure[j].da/pure[j].a)/2;
            s[i] += static_cast<long double>(w[j])*aij[i][j];
            st[i] += static_cast<long double>(w[j])*derivative;
            a += static_cast<long double>(w[i])*w[j]*aij[i][j];
            at += static_cast<long double>(w[i])*w[j]*derivative;
        }
    }
    const long double A = a*p/(r*r*t*t), B = b*p/(r*t), sqrt2 = std::sqrt(2.0L);
    const auto zs = reference_roots(A,B);
    std::vector<ReferencePhase> phases;
    for (const auto z : zs) {
        ReferencePhase result;
        result.z = z;
        result.ln_phi.resize(n);
        result.dz.resize(columns);
        result.dln.assign(n, std::vector<long double>(columns));
        const long double g = A/(2*sqrt2*B);
        const long double plus = z+(1+sqrt2)*B, minus = z+(1-sqrt2)*B;
        const long double L = std::log(plus/minus);
        for (std::size_t i = 0; i < n; ++i) {
            const long double ratio = pure[i].b/b;
            result.ln_phi[i] = ratio*(z-1)-std::log(z-B)-g*(2*s[i]/a-ratio)*L;
        }
        for (std::size_t column = 0; column < columns; ++column) {
            const long double da = column == 1 ? at : (column >= 2 ? 2*s[column-2] : 0);
            const long double db = column >= 2 ? pure[column-2].b : 0;
            const long double dA = column == 0 ? A/p : (column == 1 ? A*(at/a-2/t) : A*da/a);
            const long double dB = column == 0 ? B/p : (column == 1 ? -B/t : B*db/b);
            // Differentiated ORIGINAL cubic, independent of the production IFT expression.
            const long double fz = 3*z*z+2*(B-1)*z+A-3*B*B-2*B;
            const long double fa = z-B, fb = z*z-(6*B+2)*z+3*B*B+2*B-A;
            const long double dz = -(fa*dA+fb*dB)/fz;
            result.dz[column] = dz;
            const long double dg = (dA*B-A*dB)/(2*sqrt2*B*B);
            const long double dL = (dz+(1+sqrt2)*dB)/plus-(dz+(1-sqrt2)*dB)/minus;
            for (std::size_t i = 0; i < n; ++i) {
                const long double ds = column == 1 ? st[i] : (column >= 2 ? aij[i][column-2] : 0);
                const long double ratio = pure[i].b/b, dratio = -ratio*db/b;
                const long double e = 2*s[i]/a-ratio;
                const long double de = 2*(ds*a-s[i]*da)/(a*a)-dratio;
                result.dln[i][column] = dratio*(z-1)+ratio*dz-(dz-dB)/(z-B)-
                                       (dg*e+g*de)*L-g*e*dL;
            }
        }
        phases.push_back(std::move(result));
    }
    return phases;
}
} // namespace test
#endif
