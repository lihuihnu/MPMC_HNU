#include <mpmc/flash/pr76_split.hpp>
#include <mpmc/flash/sw92_profile_c_phase_set.hpp>
#include <mpmc/physics/component_inventory.hpp>
#include <mpmc/physics/pr76_thermodynamic_closure.hpp>
#include <mpmc/physics/sw92_thermodynamic_closure.hpp>

#include "../../flash/sw92_phase_assigned_pt/physical_sample6.hpp"
#include "../../thermodynamics/sw92/test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool component_inventory_header();

namespace {
namespace fl = mpmc::flash;
namespace ph = mpmc::physics;
namespace th = mpmc::thermodynamics;
namespace sample6 = sw92_profile_c_sample6;
using Vec = std::vector<double>;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

void near(double actual, double expected, double relative = 2e-10,
          double absolute = 2e-10,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute + relative * std::max(std::abs(actual), std::abs(expected))) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "numeric mismatch", where);
    }
}

void near_fd(double actual, double expected,
             double relative = 7e-3, double absolute = 2e-7,
             std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute + relative * std::max(std::abs(actual), std::abs(expected))) {
        std::cerr << "analytic=" << actual << " fresh_fd=" << expected << '\n';
        require(false, "0D fresh-resolve Jacobian mismatch", where);
    }
}

th::SourcedScalar sourced_scalar(double value, th::Unit unit,
                                  const th::Provenance& source,
                                  std::string original_unit,
                                  std::string conversion) {
    return {value, unit, source, std::move(original_unit), std::move(conversion)};
}

th::Pr76Phase<double> pr76_binary_model() {
    const th::Provenance source{
        th::SourceKind::literature,
        "https://academicweb.nd.edu/~markst/zm97a.pdf",
        "revised October 1997",
        "Section 4.3 and Table 4",
        "PR numerical regression, not experimental validation",
        "Read author-hosted PDF pages",
        "Limited attributed factual parameters; no article text or code reproduced"
    };
    const std::array<std::string, 2> ids{"nitrogen", "ethane"};
    const std::array<std::array<double, 3>, 2> specs{{
        {126.2, 3.39e6, 0.04},
        {305.4, 4.88e6, 0.098}
    }};
    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "Hua1997-nitrogen-ethane";
    input.revision = "split-regression-v1";
    input.applicability = {std::nullopt, std::nullopt, source};
    for (std::size_t i = 0; i < ids.size(); ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, source, {}});
        input.pure.push_back({
            ids[i],
            sourced_scalar(specs[i][0], th::Unit::kelvin, source, "K", "identity"),
            sourced_scalar(specs[i][1], th::Unit::pascal, source, "bar",
                           "bar * 100000 -> Pa"),
            sourced_scalar(specs[i][2], th::Unit::dimensionless, source,
                           "dimensionless", "identity")
        });
    }
    input.binary.push_back({
        ids[0], ids[1],
        sourced_scalar(0.08, th::Unit::dimensionless, source,
                       "dimensionless", "identity")
    });
    const std::vector<std::string> order{ids.begin(), ids.end()};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input));
}

th::Sw92Phase<double> sw92_binary_model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::co2, reverse));
}

void require_inventory_basics(const ph::PtComponentInventorySnapshot& inventory,
                              std::size_t components,
                              std::size_t phases,
                              bool derivative_expected = true) {
    require(inventory.inventory_available(), inventory.diagnostic);
    require(inventory.component_count == components &&
                inventory.phase_count == phases &&
                inventory.primal &&
                inventory.primal->component_molar_density_mol_per_m3.size() == components,
            "inventory primal shape changed");
    double sum = 0.0;
    for (std::size_t i = 0; i < components; ++i) {
        const double value = inventory.primal->component_molar_density_mol_per_m3[i];
        require(std::isfinite(value) && value >= 0.0,
                "inventory component density is invalid");
        near(value,
             inventory.primal->total_molar_density_mol_per_m3 * inventory.feed[i],
             2e-8, 2e-8);
        sum += value;
    }
    near(sum, inventory.primal->total_molar_density_mol_per_m3, 2e-8, 2e-8);
    if (derivative_expected) {
        require(inventory.linearization_available() && inventory.linearization,
                inventory.diagnostic);
        require(inventory.linearization->input_count == components + 1U,
                "inventory reduced-feed input count changed");
        for (std::size_t column = 0;
             column < inventory.linearization->input_count; ++column) {
            double derivative_sum = 0.0;
            for (std::size_t i = 0; i < components; ++i) {
                derivative_sum +=
                    inventory.linearization->d_component_molar_density(i, column);
            }
            near(derivative_sum,
                 inventory.linearization->d_total_molar_density(column),
                 2e-8, 2e-8);
        }
    }
}

ph::PtPhaseSetThermodynamicClosureSnapshot synthetic_phase_set(bool swapped) {
    ph::PtPhaseSetThermodynamicClosureSnapshot source;
    source.primal_status = ph::ThermodynamicClosurePrimalStatus::valid;
    source.linearization_status =
        ph::ThermodynamicClosureLinearizationStatus::available;
    source.linearization_reason =
        ph::ThermodynamicClosureLinearizationReason::none;
    source.pressure_pa = 1.0e6;
    source.temperature_k = 300.0;
    source.feed = {0.4, 0.6};
    source.component_ids = {"A", "B"};
    source.thermodynamic_model = "synthetic";
    source.dataset_id = "inventory-algebra";
    source.revision = "v1";
    source.primal.emplace();
    ph::ThermodynamicPhaseState a{0.25, {0.8, 0.2}, 1.0, 10000.0};
    ph::ThermodynamicPhaseState b{
        0.75, {0.2666666666666666667, 0.7333333333333333333},
        1.0, 2000.0};
    source.primal->phases = swapped
        ? std::vector<ph::ThermodynamicPhaseState>{b, a}
        : std::vector<ph::ThermodynamicPhaseState>{a, b};

    ph::PtPhaseSetThermodynamicLinearization linearization;
    linearization.component_count = 2U;
    linearization.phase_count = 2U;
    linearization.input_count = 3U;
    const double dx_b_pressure =
        -1.0e-8 * (0.8 - 0.2666666666666666667) / 0.75;
    const std::array<double,3> db_a{1.0e-8, 0.0, 0.0};
    const std::array<double,3> db_b{-1.0e-8, 0.0, 0.0};
    const std::array<double,3> dx_a0{0.0, 0.0, 1.0};
    const std::array<double,3> dx_a1{0.0, 0.0, -1.0};
    const std::array<double,3> dx_b0{dx_b_pressure, 0.0, 1.0};
    const std::array<double,3> dx_b1{-dx_b_pressure, 0.0, -1.0};
    const std::array<double,3> dc_a{0.001, -1.0, 0.0};
    const std::array<double,3> dc_b{0.0002, -0.2, 0.0};

    const auto append_phase = [&](const std::array<double,3>& db,
                                  const std::array<double,3>& dx0,
                                  const std::array<double,3>& dx1,
                                  const std::array<double,3>& dc) {
        linearization.phase_fraction_jacobian.insert(
            linearization.phase_fraction_jacobian.end(), db.begin(), db.end());
        linearization.composition_jacobian.insert(
            linearization.composition_jacobian.end(), dx0.begin(), dx0.end());
        linearization.composition_jacobian.insert(
            linearization.composition_jacobian.end(), dx1.begin(), dx1.end());
        linearization.compressibility_jacobian.insert(
            linearization.compressibility_jacobian.end(), 3U, 0.0);
        linearization.molar_density_jacobian.insert(
            linearization.molar_density_jacobian.end(), dc.begin(), dc.end());
    };
    if (swapped) {
        append_phase(db_b, dx_b0, dx_b1, dc_b);
        append_phase(db_a, dx_a0, dx_a1, dc_a);
    } else {
        append_phase(db_a, dx_a0, dx_a1, dc_a);
        append_phase(db_b, dx_b0, dx_b1, dc_b);
    }
    linearization.equilibrium_jacobian_rcond = 1.0;
    source.linearization = std::move(linearization);
    return source;
}

void synthetic_phase_permutation() {
    const auto normal = ph::build_pt_component_inventory(synthetic_phase_set(false));
    const auto swapped = ph::build_pt_component_inventory(synthetic_phase_set(true));
    require_inventory_basics(normal, 2U, 2U);
    require_inventory_basics(swapped, 2U, 2U);
    near(normal.primal->total_molar_density_mol_per_m3,
         swapped.primal->total_molar_density_mol_per_m3, 0.0, 0.0);
    require(normal.primal->component_molar_density_mol_per_m3 ==
                swapped.primal->component_molar_density_mol_per_m3,
            "phase-slot permutation changed component inventory");
    require(normal.linearization->total_molar_density_gradient ==
                swapped.linearization->total_molar_density_gradient &&
                normal.linearization->component_molar_density_jacobian ==
                    swapped.linearization->component_molar_density_jacobian,
            "phase-slot permutation changed inventory Jacobian");
}

ph::PtComponentInventorySnapshot solve_pr76_inventory(double pressure_pa,
                                                      double z0) {
    const auto model = pr76_binary_model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto split = fl::solve_pr76_pt_vle(
        pressure_pa, 270.0, Vec{z0, 1.0 - z0}, evaluator);
    const auto closure = ph::build_pr76_pt_vle_closure(split, evaluator);
    return ph::build_pt_component_inventory(closure);
}

ph::PtComponentInventorySnapshot solve_sw92_inventory(double pressure_pa,
                                                      double z0) {
    const auto model = sw92_binary_model();
    const auto source = fl::solve_sw92_profile_c_pt_phase_set(
        pressure_pa, 340.0, Vec{z0, 1.0 - z0}, model, 0.0);
    const auto closure =
        ph::build_sw92_profile_c_thermodynamic_closure(source, model);
    return ph::build_pt_component_inventory(closure.closure);
}

void pr76_fixed_vle() {
    const auto inventory = solve_pr76_inventory(7.6e6, 0.30);
    require_inventory_basics(inventory, 2U, 2U);
    require(inventory.source_closure_convention ==
                ph::ThermodynamicClosureSnapshot::convention,
            "PR76 inventory lost fixed-VLE source convention");
}

void sw92_one_and_three_phase() {
    const auto binary = sw92_binary_model();
    const auto dry_source = fl::solve_sw92_profile_c_pt_phase_set(
        3.0e6, 340.0, Vec{0.995, 0.005}, binary, 0.0);
    const auto dry_closure =
        ph::build_sw92_profile_c_thermodynamic_closure(dry_source, binary);
    const auto dry = ph::build_pt_component_inventory(dry_closure.closure);
    require_inventory_basics(dry, 2U, 1U);
    near(dry.primal->total_molar_density_mol_per_m3,
         dry_closure.closure.primal->phases.front().molar_density_mol_per_m3,
         2e-12, 2e-10);

    const auto normal_model = sample6::model(false);
    const auto reverse_model = sample6::model(true);
    const auto normal_source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(false), normal_model, 0.0);
    const auto reverse_source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(true), reverse_model, 0.0);
    const auto normal_closure = ph::build_sw92_profile_c_thermodynamic_closure(
        normal_source, normal_model);
    const auto reverse_closure = ph::build_sw92_profile_c_thermodynamic_closure(
        reverse_source, reverse_model);
    const auto normal = ph::build_pt_component_inventory(normal_closure.closure);
    const auto reverse = ph::build_pt_component_inventory(reverse_closure.closure);
    require_inventory_basics(normal, 8U, 3U);
    require_inventory_basics(reverse, 8U, 3U);
    near(normal.primal->total_molar_density_mol_per_m3,
         reverse.primal->total_molar_density_mol_per_m3, 3e-8, 3e-8);
    for (std::size_t i = 0; i < 8U; ++i) {
        near(normal.primal->component_molar_density_mol_per_m3[i],
             reverse.primal->component_molar_density_mol_per_m3[7U - i],
             4e-8, 4e-8);
        for (std::size_t column = 0; column < 2U; ++column) {
            near(normal.linearization->d_component_molar_density(i, column),
                 reverse.linearization->d_component_molar_density(
                     7U - i, column),
                 8e-6, 5e-7);
        }
    }
}

void unavailable_policy() {
    const auto sw_model = sample6::model();
    const auto sw_source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(), sw_model, 0.0);
    ph::Sw92ProfileCClosureOptions sw_options;
    sw_options.sensitivity.minimum_derivative_phase_fraction = 0.04;
    const auto sw_closure = ph::build_sw92_profile_c_thermodynamic_closure(
        sw_source, sw_model, sw_options);
    const auto sw_inventory = ph::build_pt_component_inventory(sw_closure.closure);
    require_inventory_basics(sw_inventory, 8U, 3U, false);
    require(sw_inventory.linearization_status ==
                ph::ThermodynamicClosureLinearizationStatus::unavailable &&
                sw_inventory.linearization_reason ==
                    ph::ThermodynamicClosureLinearizationReason::phase_boundary &&
                !sw_inventory.linearization_available(),
            "SW derivative-only failure did not preserve inventory primal");

    const auto pr_model = pr76_binary_model();
    fl::Pr76VleEvaluator evaluator(pr_model);
    const auto split = fl::solve_pr76_pt_vle(
        7.6e6, 270.0, Vec{0.30, 0.70}, evaluator);
    ph::Pr76PtVleClosureOptions pr_options;
    pr_options.sensitivity.minimum_derivative_phase_fraction = 0.45;
    const auto pr_closure = ph::build_pr76_pt_vle_closure(
        split, evaluator, pr_options);
    const auto pr_inventory = ph::build_pt_component_inventory(pr_closure);
    require_inventory_basics(pr_inventory, 2U, 2U, false);
    require(pr_inventory.linearization_reason ==
                ph::ThermodynamicClosureLinearizationReason::phase_boundary &&
                !pr_inventory.linearization_available(),
            "PR derivative-only failure did not preserve inventory primal");

    auto inconsistent = sw_closure.closure;
    inconsistent.feed[0] += 1.0e-4;
    inconsistent.feed[1] -= 1.0e-4;
    const auto rejected = ph::build_pt_component_inventory(inconsistent);
    require(!rejected.inventory_available() &&
                rejected.primal_status ==
                    ph::ThermodynamicClosurePrimalStatus::indeterminate,
            "phase/feed inconsistency was accepted as component inventory");

    bool caught = false;
    try {
        ph::PtComponentInventoryOptions invalid;
        invalid.maximum_primal_identity_error = 0.0;
        (void)ph::build_pt_component_inventory(sw_closure.closure, invalid);
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require(caught, "invalid inventory integrity options were not rejected");
}

struct ZeroDResidual {
    std::vector<double> values;
    std::vector<double> jacobian; // N residuals x N variables (p,z_0..z_{N-2}).
};

ZeroDResidual zero_d_fixed_volume_residual(
    const ph::PtComponentInventorySnapshot& inventory,
    double fluid_volume_m3, std::span<const double> target_moles) {
    require(inventory.inventory_available() && inventory.linearization_available() &&
                inventory.primal && inventory.linearization,
            "0D residual requires inventory primal and local Jacobian");
    const std::size_t n = inventory.component_count;
    require(target_moles.size() == n && std::isfinite(fluid_volume_m3) &&
                fluid_volume_m3 > 0.0,
            "0D residual fixture shape/volume invalid");
    ZeroDResidual result;
    result.values.assign(n, 0.0);
    result.jacobian.assign(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        result.values[i] = fluid_volume_m3 *
            inventory.primal->component_molar_density_mol_per_m3[i] -
            target_moles[i];
        for (std::size_t column = 0; column < n; ++column) {
            const std::size_t q_column = column == 0U ? 0U : column + 1U;
            result.jacobian[i * n + column] = fluid_volume_m3 *
                inventory.linearization->d_component_molar_density(i, q_column);
        }
    }
    return result;
}

template <class Solve>
void check_zero_d_binary(Solve&& solve, double pressure_pa, double z0,
                         double pressure_step, double feed_step) {
    constexpr double volume = 0.125;
    const auto base = solve(pressure_pa, z0);
    require_inventory_basics(base, 2U, 2U);
    std::array<double,2> target{
        volume * base.primal->component_molar_density_mol_per_m3[0],
        volume * base.primal->component_molar_density_mol_per_m3[1]};
    const auto residual = zero_d_fixed_volume_residual(base, volume, target);
    near(residual.values[0], 0.0, 0.0, 2e-12);
    near(residual.values[1], 0.0, 0.0, 2e-12);

    const std::array<double,2> step{pressure_step, feed_step};
    for (std::size_t column = 0; column < 2U; ++column) {
        double pp = pressure_pa;
        double pm = pressure_pa;
        double zp = z0;
        double zm = z0;
        if (column == 0U) {
            pp += step[column];
            pm -= step[column];
        } else {
            zp += step[column];
            zm -= step[column];
        }
        const auto plus = solve(pp, zp);
        const auto minus = solve(pm, zm);
        require(plus.inventory_available() && minus.inventory_available(),
                "fresh 0D perturbation left accepted inventory domain");
        for (std::size_t row = 0; row < 2U; ++row) {
            const double finite_difference = volume *
                (plus.primal->component_molar_density_mol_per_m3[row] -
                 minus.primal->component_molar_density_mol_per_m3[row]) /
                (2.0 * step[column]);
            near_fd(residual.jacobian[row * 2U + column], finite_difference);
        }
    }
}

void zero_d_pr76() {
    check_zero_d_binary(
        [](double p, double z) { return solve_pr76_inventory(p, z); },
        7.6e6, 0.30, 500.0, 1.0e-5);
}

void zero_d_sw92() {
    check_zero_d_binary(
        [](double p, double z) { return solve_sw92_inventory(p, z); },
        3.0e6, 0.70, 500.0, 1.0e-5);
}

void headers() {
    require(component_inventory_header(),
            "component inventory public-header probe failed");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"synthetic_phase_permutation", synthetic_phase_permutation},
    {"pr76_fixed_vle", pr76_fixed_vle},
    {"sw92_one_and_three_phase", sw92_one_and_three_phase},
    {"unavailable_policy", unavailable_policy},
    {"zero_d_pr76", zero_d_pr76},
    {"zero_d_sw92", zero_d_sw92},
    {"headers", headers}
};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument("one test name required");
        }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout << "[PASS] " << name << '\n';
                return 0;
            }
        }
        throw std::invalid_argument("unknown test");
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
