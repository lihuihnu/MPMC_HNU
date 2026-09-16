#include <mpmc/thermodynamics/pr76_phase_ad.hpp>
#include "test_support.hpp"

#include <array>
#include <cstddef>
#include <exception>
#include <iostream>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace {
using namespace test;

template <std::floating_point T>
void check_direct_ad() {
    constexpr std::size_t width = 2;
    th::Pr76PhaseAdWorkspace<T, width> workspace;
    auto model = kernel<T>();
    const T pressure = T{1000000};
    const T temperature = T{450};
    const std::vector<T> x{T{1}/T{4}, T{1}/T{2}, T{1}/T{4}};
    const auto ref = reference(pressure, temperature, x, {0, 1, 2});

    for (std::size_t root = 0; root < ref.size(); ++root) {
        const auto full = workspace.value_and_jacobian_full(
            model, pressure, temperature, std::span<const T>{x}, root);
        require(full.input_count == 5 && full.output_count == 4,
                "direct full AD shape/order contract");
        near(full.values[0], ref[root].z);
        for (std::size_t i = 0; i < x.size(); ++i) {
            near(full.values[i + 1], ref[root].ln_phi[i]);
        }
        for (std::size_t column = 0; column < full.input_count; ++column) {
            near(full.jacobian[column], ref[root].dz[column]);
            for (std::size_t i = 0; i < x.size(); ++i) {
                near(full.jacobian[(i + 1) * full.input_count + column],
                     ref[root].dln[i][column]);
            }
        }

        const std::array<T, 2> independent{x[0], x[1]};
        const auto reduced = workspace.value_and_jacobian_reduced(
            model, pressure, temperature, std::span<const T>{independent}, root);
        require(reduced.input_count == 4 && reduced.output_count == 4,
                "direct reduced AD shape/order contract");
        for (std::size_t row = 0; row < reduced.output_count; ++row) {
            near(reduced.values[row], full.values[row]);
        }
        const std::size_t dependent_column = 4; // full x_2 column.
        for (std::size_t column = 0; column < reduced.input_count; ++column) {
            const long double expected_z = column < 2
                ? ref[root].dz[column]
                : ref[root].dz[column] - ref[root].dz[dependent_column];
            near(reduced.jacobian[column], expected_z);
            for (std::size_t i = 0; i < x.size(); ++i) {
                const long double expected_phi = column < 2
                    ? ref[root].dln[i][column]
                    : ref[root].dln[i][column] - ref[root].dln[i][dependent_column];
                near(reduced.jacobian[(i + 1) * reduced.input_count + column], expected_phi);
            }
        }
    }

    // The same adapter workspace must support a runtime component-count change,
    // including the n=1 reduced coordinate system with no composition inputs.
    auto pure_model = kernel<T>({0});
    const std::vector<T> pure_x{T{1}};
    const auto pure_ref = reference(T{1000000}, T{350}, pure_x, {0});
    th::Pr76PhaseWorkspace<T> plain;
    const auto roots = pure_model.roots_reduced(T{1000000}, T{350}, {}, plain);
    require(roots.status == th::Pr76RootStatus::success && roots.count == pure_ref.size(),
            "pure direct-AD root fixture");
    const std::size_t root = roots.count - 1;
    const auto pure = workspace.value_and_jacobian_reduced(
        pure_model, T{1000000}, T{350}, std::span<const T>{}, root);
    require(pure.input_count == 2 && pure.output_count == 2,
            "pure reduced direct AD shape");
    near(pure.values[0], pure_ref[root].z);
    near(pure.values[1], pure_ref[root].ln_phi[0]);
    for (std::size_t column = 0; column < 2; ++column) {
        near(pure.jacobian[column], pure_ref[root].dz[column]);
        near(pure.jacobian[2 + column], pure_ref[root].dln[0][column]);
    }

    // Resource limits are checked by the adapter before growing its input scratch.
    ad::RuntimeJacobianLimits limits;
    limits.max_inputs = 4;
    expect_error<std::length_error>([&] {
        (void)workspace.value_and_jacobian_full(
            model, pressure, temperature, std::span<const T>{x}, 0, {}, limits);
    });
}

} // namespace

int main() {
    try {
        check_direct_ad<float>();
        check_direct_ad<double>();
        check_direct_ad<long double>();
        std::cout << "[PASS] direct_ad\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
