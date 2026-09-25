#include <mpmc/flow_discretization/tpfa_phase_darcy_flux.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool tpfa_phase_darcy_flux_header();

namespace {

namespace disc = mpmc::discretization;
namespace flow = mpmc::flow;
namespace fd = mpmc::flow_discretization;
namespace mesh = mpmc::mesh;

void require(
    bool condition,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(
            std::string{where.file_name()} +
            ":" +
            std::to_string(where.line()) +
            ": " +
            std::string{message});
    }
}

void near(
    double actual,
    double expected,
    double relative = 2.0e-13,
    double absolute = 2.0e-13,
    std::source_location where =
        std::source_location::current()) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute +
                relative *
                    std::max(
                        std::abs(actual),
                        std::abs(expected))) {
        std::cerr
            << "actual=" << actual
            << " expected=" << expected
            << '\n';
        require(false, "numeric mismatch", where);
    }
}

template <class Function>
void expect_invalid(
    Function&& function,
    std::string_view fragment) {
    try {
        function();
    } catch (const std::invalid_argument& error) {
        require(
            std::string_view{error.what()}.find(fragment) !=
                std::string_view::npos,
            "invalid_argument diagnostic changed");
        return;
    }
    throw std::runtime_error(
        "expected std::invalid_argument");
}

disc::CombinedTransmissibilityAdmissibility3D
direct_admissibility() {
    return {
        disc::
            CombinedTransmissibilityAdmissibilityDisposition3D::
                direct_normal_projection_k_orthogonal_candidate,
        {
            disc::
                TransmissibilityGeometryDisposition3D::
                    direct_normal_projection_allowed,
            0.0,
            0.05},
        {
            disc::
                KOrthogonalityDisposition3D::
                    k_orthogonal_within_policy,
            std::nullopt,
            std::nullopt,
            {std::nullopt, std::nullopt},
            0.05}};
}

disc::TpfaInternalFaceTransmissibilityEntry3D
materialized_entry(double tf_m3 = 2.0e-12) {
    return {
        mesh::LocalIndex{
            mesh::LocalIndex::value_type{7}},
        disc::
            TpfaInternalFaceTransmissibilityDisposition3D::
                materialized,
        direct_admissibility(),
        disc::TpfaStaticFaceTransmissibility3D{
            disc::
                TpfaStaticFaceTransmissibilityDisposition3D::
                    positive_harmonic_combination,
            1.5,
            tf_m3}};
}

flow::NaturalVariableStateIdentity3P
identity(bool owner) {
    const auto pivot =
        flow::NaturalVariableCompositionPivot3P::
            from_dependent_components(
                3U,
                owner
                    ? std::array<std::size_t, 3>{1U, 0U, 2U}
                    : std::array<std::size_t, 3>{0U, 2U, 1U});

    return {
        flow::NaturalVariableLayout3P{pivot},
        {"A", "B", "C"},
        owner ? 1.0e6 : 0.99e6,
        owner ? 350.0 : 351.0,
        owner
            ? std::array<double, 3>{0.20, 0.30, 0.50}
            : std::array<double, 3>{0.25, 0.35, 0.40},
        owner
            ? std::array<std::vector<double>, 3>{
                  std::vector<double>{0.10, 0.70, 0.20},
                  std::vector<double>{0.60, 0.20, 0.20},
                  std::vector<double>{0.20, 0.30, 0.50}}
            : std::array<std::vector<double>, 3>{
                  std::vector<double>{0.55, 0.25, 0.20},
                  std::vector<double>{0.25, 0.25, 0.50},
                  std::vector<double>{0.20, 0.60, 0.20}}};
}

std::array<std::vector<double>, 3>
gradient(
    const flow::NaturalVariableLayout3P& layout,
    double scale) {
    std::array<std::vector<double>, 3> result;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        result[phase].resize(layout.unknown_count());
        for (std::size_t column = 0U;
             column < layout.unknown_count();
             ++column) {
            result[phase][column] =
                scale *
                static_cast<double>(
                    1U +
                    phase * 10U +
                    column);
        }
    }
    return result;
}

flow::TwoCellPhasePotentialUpwindLinearization3P
potential_payload() {
    const auto owner_identity =
        identity(true);
    const auto neighbour_identity =
        identity(false);

    const auto owner_dphi =
        gradient(owner_identity.layout, 10.0);
    const auto neighbour_dphi =
        gradient(neighbour_identity.layout, 20.0);
    const auto owner_dlambda =
        gradient(owner_identity.layout, 1.0e-3);
    const auto neighbour_dlambda =
        gradient(neighbour_identity.layout, 2.0e-3);

    flow::TwoCellPhasePotentialUpwindLinearization3P result{
        flow::FacePhaseDensityPolicy3P::
            arithmetic_mean_owner_neighbour,
        {0.0, 0.0, -9.81},
        {0.0, 0.0, 1.0},
        owner_identity,
        neighbour_identity,
        {}};

    result.phase[0] = {
        900.0,
        std::vector<double>(
            owner_identity.layout.unknown_count(),
            0.0),
        std::vector<double>(
            neighbour_identity.layout.unknown_count(),
            0.0),
        -9.81,
        -8829.0,
        10000.0,
        owner_dphi[0],
        neighbour_dphi[0],
        flow::UpwindCellSelection3P::
            neighbour_positive_phase_potential,
        4.0,
        std::vector<double>(
            owner_identity.layout.unknown_count(),
            0.0),
        neighbour_dlambda[0]};

    result.phase[1] = {
        500.0,
        std::vector<double>(
            owner_identity.layout.unknown_count(),
            0.0),
        std::vector<double>(
            neighbour_identity.layout.unknown_count(),
            0.0),
        -9.81,
        -4905.0,
        -5000.0,
        owner_dphi[1],
        neighbour_dphi[1],
        flow::UpwindCellSelection3P::
            owner_negative_phase_potential,
        2.0,
        owner_dlambda[1],
        std::vector<double>(
            neighbour_identity.layout.unknown_count(),
            0.0)};

    result.phase[2] = {
        750.0,
        std::vector<double>(
            owner_identity.layout.unknown_count(),
            0.0),
        std::vector<double>(
            neighbour_identity.layout.unknown_count(),
            0.0),
        -9.81,
        -7357.5,
        0.0,
        owner_dphi[2],
        neighbour_dphi[2],
        flow::UpwindCellSelection3P::
            owner_exact_zero_tie,
        3.0,
        owner_dlambda[2],
        std::vector<double>(
            neighbour_identity.layout.unknown_count(),
            0.0)};

    return result;
}

void materialized_phase_flux() {
    constexpr double tf = 2.0e-12;
    const auto entry =
        materialized_entry(tf);
    const auto potential =
        potential_payload();

    const auto result =
        fd::
            build_materialized_tpfa_internal_face_phase_darcy_flux(
                entry,
                potential);

    require(
        result.face ==
                mesh::LocalIndex{
                    mesh::LocalIndex::value_type{7}} &&
            result.owner_state_identity.component_ids ==
                potential.owner_state_identity.component_ids &&
            result.neighbour_state_identity.component_ids ==
                potential.neighbour_state_identity.component_ids,
        "materialized phase flux lost face/state identity");
    near(
        result.static_transmissibility_m3,
        tf);

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto& source =
            potential.phase[phase];
        const auto& flux =
            result.phase[phase];
        const double expected =
            -tf *
            source.upwind_mobility_per_pa_s *
            source.phase_potential_difference_pa;

        near(
            flux.volumetric_flux_m3_per_s,
            expected,
            2.0e-13,
            2.0e-20);
        require(
            flux.upwind_selection ==
                source.upwind_selection,
            "phase flux lost upwind selection");

        for (std::size_t column = 0U;
             column <
             source.owner_phase_potential_gradient.size();
             ++column) {
            const double expected_gradient =
                -tf *
                (source.upwind_mobility_per_pa_s *
                     source.owner_phase_potential_gradient[
                         column] +
                 source.phase_potential_difference_pa *
                     source.owner_upwind_mobility_gradient[
                         column]);
            near(
                flux.owner_flux_gradient[column],
                expected_gradient,
                2.0e-13,
                2.0e-20);
        }

        for (std::size_t column = 0U;
             column <
             source.neighbour_phase_potential_gradient.size();
             ++column) {
            const double expected_gradient =
                -tf *
                (source.upwind_mobility_per_pa_s *
                     source.neighbour_phase_potential_gradient[
                         column] +
                 source.phase_potential_difference_pa *
                     source.neighbour_upwind_mobility_gradient[
                         column]);
            near(
                flux.neighbour_flux_gradient[column],
                expected_gradient,
                2.0e-13,
                2.0e-20);
        }
    }

    require(
        result.phase[0].volumetric_flux_m3_per_s < 0.0 &&
            result.phase[1].volumetric_flux_m3_per_s > 0.0 &&
            result.phase[2].volumetric_flux_m3_per_s == 0.0,
        "phase volumetric flux sign convention changed");
}

flow::TwoCellPhasePotentialUpwindLinearization3P
perturb_potential(
    flow::TwoCellPhasePotentialUpwindLinearization3P value,
    bool owner,
    std::size_t column,
    double delta) {
    for (auto& phase : value.phase) {
        if (owner) {
            phase.phase_potential_difference_pa +=
                phase.owner_phase_potential_gradient[column] *
                delta;
            phase.upwind_mobility_per_pa_s +=
                phase.owner_upwind_mobility_gradient[column] *
                delta;
        } else {
            phase.phase_potential_difference_pa +=
                phase.neighbour_phase_potential_gradient[column] *
                delta;
            phase.upwind_mobility_per_pa_s +=
                phase.neighbour_upwind_mobility_gradient[column] *
                delta;
        }
    }
    return value;
}

void flux_jacobian_fresh_perturbation() {
    constexpr double tf = 2.0e-12;
    constexpr double step = 1.0e-6;
    const auto entry =
        materialized_entry(tf);
    const auto base_potential =
        potential_payload();
    const auto base =
        fd::
            build_materialized_tpfa_internal_face_phase_darcy_flux(
                entry,
                base_potential);

    for (std::size_t column = 0U;
         column <
         base.owner_state_identity.layout.unknown_count();
         ++column) {
        const auto plus =
            fd::
                build_materialized_tpfa_internal_face_phase_darcy_flux(
                    entry,
                    perturb_potential(
                        base_potential,
                        true,
                        column,
                        step));
        const auto minus =
            fd::
                build_materialized_tpfa_internal_face_phase_darcy_flux(
                    entry,
                    perturb_potential(
                        base_potential,
                        true,
                        column,
                        -step));

        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            const double finite_difference =
                (plus.phase[phase]
                     .volumetric_flux_m3_per_s -
                 minus.phase[phase]
                     .volumetric_flux_m3_per_s) /
                (2.0 * step);
            near(
                base.phase[phase]
                    .owner_flux_gradient[column],
                finite_difference,
                2.0e-7,
                2.0e-13);
        }
    }

    for (std::size_t column = 0U;
         column <
         base.neighbour_state_identity.layout.unknown_count();
         ++column) {
        const auto plus =
            fd::
                build_materialized_tpfa_internal_face_phase_darcy_flux(
                    entry,
                    perturb_potential(
                        base_potential,
                        false,
                        column,
                        step));
        const auto minus =
            fd::
                build_materialized_tpfa_internal_face_phase_darcy_flux(
                    entry,
                    perturb_potential(
                        base_potential,
                        false,
                        column,
                        -step));

        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            const double finite_difference =
                (plus.phase[phase]
                     .volumetric_flux_m3_per_s -
                 minus.phase[phase]
                     .volumetric_flux_m3_per_s) /
                (2.0 * step);
            near(
                base.phase[phase]
                    .neighbour_flux_gradient[column],
                finite_difference,
                2.0e-7,
                2.0e-13);
        }
    }
}

void exact_zero_tie_branch() {
    constexpr double tf = 2.0e-12;
    const auto potential =
        potential_payload();
    const auto result =
        fd::
            build_materialized_tpfa_internal_face_phase_darcy_flux(
                materialized_entry(tf),
                potential);

    const auto& source = potential.phase[2];
    const auto& flux = result.phase[2];

    require(
        source.upwind_selection ==
                flow::UpwindCellSelection3P::
                    owner_exact_zero_tie &&
            flux.volumetric_flux_m3_per_s == 0.0,
        "exact-zero phase-potential tie semantics changed");

    for (std::size_t column = 0U;
         column <
         source.owner_phase_potential_gradient.size();
         ++column) {
        near(
            flux.owner_flux_gradient[column],
            -tf *
                source.upwind_mobility_per_pa_s *
                source.owner_phase_potential_gradient[
                    column],
            2.0e-13,
            2.0e-20);
    }
    for (std::size_t column = 0U;
         column <
         source.neighbour_phase_potential_gradient.size();
         ++column) {
        near(
            flux.neighbour_flux_gradient[column],
            -tf *
                source.upwind_mobility_per_pa_s *
                source.neighbour_phase_potential_gradient[
                    column],
            2.0e-13,
            2.0e-20);
    }
}

void gate_rejection() {
    const auto potential =
        potential_payload();

    auto blocked =
        materialized_entry();
    blocked.disposition =
        disc::
            TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_geometry_non_orthogonal;
    expect_invalid(
        [&] {
            (void)fd::
                build_materialized_tpfa_internal_face_phase_darcy_flux(
                    blocked,
                    potential);
        },
        "materialized");

    auto wrong_admissibility =
        materialized_entry();
    wrong_admissibility.admissibility.disposition =
        disc::
            CombinedTransmissibilityAdmissibilityDisposition3D::
                requires_k_non_orthogonal_treatment;
    expect_invalid(
        [&] {
            (void)fd::
                build_materialized_tpfa_internal_face_phase_darcy_flux(
                    wrong_admissibility,
                    potential);
        },
        "materialized");

    auto zero_tf =
        materialized_entry();
    zero_tf.static_transmissibility =
        disc::TpfaStaticFaceTransmissibility3D{
            disc::
                TpfaStaticFaceTransmissibilityDisposition3D::
                    zero_due_to_one_half,
            1.5,
            0.0};
    expect_invalid(
        [&] {
            (void)fd::
                build_materialized_tpfa_internal_face_phase_darcy_flux(
                    zero_tf,
                    potential);
        },
        "positive static transmissibility");

    auto missing_tf =
        materialized_entry();
    missing_tf.static_transmissibility.reset();
    expect_invalid(
        [&] {
            (void)fd::
                build_materialized_tpfa_internal_face_phase_darcy_flux(
                    missing_tf,
                    potential);
        },
        "materialized");

    auto malformed =
        potential;
    malformed.phase[0]
        .owner_phase_potential_gradient
        .pop_back();
    expect_invalid(
        [&] {
            (void)fd::
                build_materialized_tpfa_internal_face_phase_darcy_flux(
                    materialized_entry(),
                    malformed);
        },
        "malformed");
}

void headers() {
    require(
        tpfa_phase_darcy_flux_header(),
        "flow-discretization phase-flux header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"materialized_phase_flux", materialized_phase_flux},
    {"flux_jacobian_fresh_perturbation", flux_jacobian_fresh_perturbation},
    {"exact_zero_tie_branch", exact_zero_tie_branch},
    {"gate_rejection", gate_rejection},
    {"headers", headers}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "one test name required");
        }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout
                    << "[PASS] "
                    << name
                    << '\n';
                return 0;
            }
        }
        throw std::invalid_argument(
            "unknown test");
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';
        return 1;
    }
}
