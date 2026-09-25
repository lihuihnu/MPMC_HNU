#include <mpmc/flow_discretization_petsc/sw92_transactional_phase_transition_restart.hpp>

#include "../../flash/sw92_phase_assigned_pt/physical_sample6.hpp"

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace {

namespace dp = mpmc::discretization_petsc;
namespace fdp = mpmc::flow_discretization_petsc;
namespace flow = mpmc::flow;
namespace fl = mpmc::flash;
namespace mesh = mpmc::mesh;
namespace sample6 = sw92_profile_c_sample6;
namespace th = mpmc::thermodynamics;

void require_sample6_transaction(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

th::Provenance synthetic_mass_source(
    std::string locator) {
    return {
        th::SourceKind::synthetic_test,
        "MPMC_HNU manufactured SW92 Sample-6 transactional fixture",
        "v1",
        std::move(locator),
        "Nonphysical positive molar mass used only to exercise selected-phase mass density and PETSc transaction mechanics.",
        "Constructed inside the regression from the published Sample-6 thermodynamic snapshot.",
        "Test-only generated fixture"};
}

th::Sw92ParameterSet
sample6_transaction_parameters() {
    const auto source =
        sample6::parameters();

    std::vector<th::Component> catalog{
        source.components().items().begin(),
        source.components().items().end()};
    std::vector<std::string> order;
    order.reserve(catalog.size());
    for (std::size_t i = 0U;
         i < catalog.size();
         ++i) {
        order.push_back(catalog[i].id);
        catalog[i].molar_mass =
            th::SourcedScalar{
                0.020 +
                    0.005 *
                        static_cast<double>(i),
                th::Unit::kilogram_per_mole,
                synthetic_mass_source(
                    catalog[i].id +
                    std::string{
                        " manufactured molar mass"}),
                "manufactured kg/mol",
                "identity"};
    }

    th::Sw92ParameterInput input;
    input.model_id =
        std::string{
            th::sw92_corrected_profile};
    input.dataset_id =
        source.dataset_id();
    input.revision =
        source.revision();
    input.applicability =
        source.applicability();
    input.pure.assign(
        source.pure_records().begin(),
        source.pure_records().end());
    input.water_binary.assign(
        source.water_binary_records().begin(),
        source.water_binary_records().end());
    input.nonwater_binary.assign(
        source.nonwater_binary_records().begin(),
        source.nonwater_binary_records().end());

    return th::Sw92ParameterSet::create(
        catalog,
        order,
        input,
        th::DataPolicy::allow_synthetic_tests);
}

struct ManufacturedSample6TransportCaloricProvider {
    template <typename Number>
    [[nodiscard]]
    flow::SelectedPhaseTransportCaloricValues<Number>
    operator()(
        const th::Sw92SelectedPhase<double>&,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number>,
        const Number&,
        const Number&,
        const Number& mass_density_kg_per_m3) const {
        // Test-only transport/caloric closure.  The production selected-phase
        // closure computes u=h-p/rho, so this manufactured h gives
        // u(T)=u0+cp*(T-350 K): topology-independent yet temperature-sensitive.
        // The nonzero du/dT keeps the energy row independent from the component
        // accumulation rows while still isolating the transaction from missing
        // eight-component source-complete caloric data.
        constexpr double u0_j_per_kg = 2.5e5;
        constexpr double cp_j_per_kg_k = 1.0e3;
        const Number internal_energy =
            Number{u0_j_per_kg} +
            Number{cp_j_per_kg_k} *
                (temperature_k - Number{350.0});
        return {
            Number{1.0e-5},
            internal_energy +
                pressure_pa /
                    mass_density_kg_per_m3};
    }
};

flow::SelectedPhasePropertyProvenance
manufactured_sample6_provenance() {
    return {
        {
            "SW92 selected molar density x manufactured test molar mass",
            "sample6-transaction-manufactured",
            "v1"},
        {
            "manufactured constant viscosity",
            "sample6-transaction-manufactured",
            "v1"},
        {
            "manufactured h=u0+cp*(T-350K)+p/rho",
            "sample6-transaction-manufactured",
            "v1"},
        {
            "production identity u=h-p/rho",
            "sample6-transaction-manufactured",
            "v1"}};
}

mesh::PartitionSnapshot
sample6_partition() {
    mesh::Topology topology{
        {
            {},
            {},
            {},
            {mesh::GlobalEntityId{
                UINT64_C(9201)}}},
        {}};
    mesh::EntityOwnerRanks owners;
    owners.cells = {
        mesh::PartitionRank{0U}};
    return mesh::PartitionSnapshot::create(
        topology,
        mesh::PartitionRank{0U},
        1U,
        std::move(owners));
}

dp::ParallelOwnedConnectionSchedule3D
sample6_schedule() {
    return {
        mesh::PartitionRank{0U},
        1U,
        {},
        {}};
}

dp::PetscMpiAijSymbolicPreallocation3D
sample6_bridge() {
    return {
        mesh::PartitionRank{0U},
        1U,
        0,
        1,
        1,
        {mesh::LocalIndex{0U}},
        {mesh::GlobalEntityId{
            UINT64_C(9201)}},
        {0},
        {1},
        {0},
        {0}};
}

dp::OwnedCellStructuralColumnPatternSnapshot3D
sample6_pattern() {
    return {
        mesh::PartitionRank{0U},
        1U,
        1U,
        0,
        1,
        1,
        {{
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{
                UINT64_C(9201)},
            0,
            0U,
            1U,
            0U,
            0U}},
        {0},
        {}};
}

PetscErrorCode sample6_zero_rock(
    const flow::NaturalVariableLayoutDescriptor& layout,
    std::span<const double>,
    void*,
    std::optional<
        fdp::Sw92RockThermalStorageLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->emplace(
        fdp::Sw92RockThermalStorageLinearization3D{
            0.0,
            std::vector<double>(
                layout.unknown_count(),
                0.0),
            {
                "zero rock storage",
                "sample6-transaction-manufactured",
                "v1"}});
    *status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    return PETSC_SUCCESS;
}

PetscErrorCode sample6_two_phase_kr(
    const flow::NaturalVariableCellState2P& state,
    void*,
    std::optional<
        fdp::Sw92TwoPhaseRelativePermeabilityLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    const std::size_t q =
        state.layout().unknown_count();
    output->emplace(
        fdp::Sw92TwoPhaseRelativePermeabilityLinearization3D{
            {1.0, 1.0},
            {
                std::vector<double>(q, 0.0),
                std::vector<double>(q, 0.0)}});
    *status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    return PETSC_SUCCESS;
}

struct Sample6LinearThreePhaseKr {
    template <typename Number>
    [[nodiscard]]
    flow::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const flow::ThreePhaseSaturationState3P<Number>& state)
        const {
        return {state.saturation};
    }
};

PetscErrorCode sample6_three_phase_saturation(
    const flow::NaturalVariableCellState3P& state,
    void*,
    std::optional<
        flow::
            ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P>*
                output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    try {
        const auto primal =
            flow::
                evaluate_three_phase_saturation_constitutive(
                    state,
                    Sample6LinearThreePhaseKr{},
                    flow::NoCapillaryPressure3P{});
        flow::
            ThreePhaseSaturationCoordinateDerivatives3P
                derivatives{};
        derivatives.relative_permeability[0] =
            {1.0, 0.0};
        derivatives.relative_permeability[1] =
            {0.0, 1.0};
        derivatives.relative_permeability[2] =
            {-1.0, -1.0};
        output->emplace(
            flow::
                make_saturation_constitutive_natural_variable_linearization(
                    state,
                    primal,
                    derivatives));
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    } catch (const std::exception&) {
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                domain_error;
    }
    return PETSC_SUCCESS;
}

PetscErrorCode sample6_no_absent_provider(
    const flow::
        AbsentPhaseThermodynamicCoordinateExtension&,
    void* user_context,
    std::optional<
        flow::
            AbsentPhasePotentialExtensionLinearization>*,
    fdp::NaturalVariableSnesEvaluationStatus3D*) {
    return user_context == nullptr
        ? PETSC_ERR_ARG_NULL
        : PETSC_ERR_SUP;
}

flow::FrozenPhysicalPhaseIdentity
sample6_identity(
    std::string key) {
    return {
        "SW92/Profile-C/sample6-transaction-test",
        std::move(key)};
}

PetscErrorCode sample6_target_identity(
    mesh::GlobalEntityId,
    const fdp::Sw92AuthoritativeTargetMaterialization3D& target,
    const flow::FrozenActivePhaseIdentityMap& source,
    void*,
    std::optional<
        flow::FrozenActivePhaseIdentityMap>* output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    const auto aq =
        sample6_identity("aqueous");
    const auto na =
        sample6_identity("nonaqueous");
    const auto extra =
        sample6_identity("extra");

    if (source.phase_count() == 2U &&
        target.phases.size() == 3U &&
        target.phases[0].metadata.thermodynamic_family ==
            th::SwPhaseFamily::aqueous &&
        target.phases[1].metadata.thermodynamic_family ==
            th::SwPhaseFamily::nonaqueous &&
        target.phases[2].metadata.thermodynamic_family ==
            th::SwPhaseFamily::nonaqueous &&
        source.find(aq).has_value() &&
        source.find(na).has_value()) {
        output->emplace(
            std::vector<
                flow::FrozenPhysicalPhaseIdentity>{
                aq,
                na,
                extra});
        return PETSC_SUCCESS;
    }

    if (source.phase_count() == 3U &&
        target.phases.size() == 2U &&
        target.phases[0].metadata.thermodynamic_family ==
            th::SwPhaseFamily::aqueous &&
        target.phases[1].metadata.thermodynamic_family ==
            th::SwPhaseFamily::nonaqueous &&
        source.find(aq).has_value() &&
        source.find(na).has_value()) {
        output->emplace(
            std::vector<
                flow::FrozenPhysicalPhaseIdentity>{
                aq,
                na});
        return PETSC_SUCCESS;
    }

    return PETSC_ERR_ARG_INCOMP;
}

std::size_t dependent_largest(
    std::span<const double> composition) {
    const auto found =
        std::max_element(
            composition.begin(),
            composition.end());
    require_sample6_transaction(
        found != composition.end() &&
            std::isfinite(*found) &&
            *found > 0.0,
        "invalid Sample-6 composition");
    return static_cast<std::size_t>(
        std::distance(
            composition.begin(),
            found));
}

std::vector<double> sample6_edge_feed(
    std::span<const double> first,
    std::span<const double> second,
    double second_fraction) {
    require_sample6_transaction(
        first.size() == second.size() &&
            second_fraction > 0.0 &&
            second_fraction < 1.0,
        "invalid Sample-6 edge feed");
    std::vector<double> result(
        first.size(),
        0.0);
    for (std::size_t i = 0U;
         i < result.size();
         ++i) {
        result[i] =
            (1.0 - second_fraction) *
                first[i] +
            second_fraction *
                second[i];
    }
    return result;
}

template <typename Closure>
std::vector<double> two_phase_q(
    const flow::NaturalVariableLayout2P& layout,
    double pressure_pa,
    double temperature_k,
    const std::array<std::vector<double>, 2>& compositions,
    const std::array<double, 2>& mole_fractions,
    const Closure& closure) {
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        pressure_pa;
    q[layout.temperature_unknown_index()] =
        temperature_k;

    std::array<double, 2> weight{};
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        const auto properties =
            closure.evaluate(
                phase,
                pressure_pa,
                temperature_k,
                std::span<const double>{
                    compositions[phase]});
        weight[phase] =
            mole_fractions[phase] /
            properties
                .molar_density_mol_per_m3;
    }
    q[layout.independent_saturation_unknown_index()] =
        weight[0] /
        (weight[0] + weight[1]);

    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        for (std::size_t component = 0U;
             component <
                 compositions[phase].size();
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        phase,
                        component);
            if (column.has_value()) {
                q[*column] =
                    compositions[phase][component];
            }
        }
    }
    return q;
}

template <typename Closure>
std::vector<double> three_phase_q(
    const flow::NaturalVariableLayout3P& layout,
    double pressure_pa,
    double temperature_k,
    const std::array<std::vector<double>, 3>& compositions,
    const std::array<double, 3>& mole_fractions,
    const Closure& closure) {
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        pressure_pa;
    q[layout.temperature_unknown_index()] =
        temperature_k;

    std::array<double, 3> weight{};
    double sum = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto properties =
            closure.evaluate(
                phase,
                pressure_pa,
                temperature_k,
                std::span<const double>{
                    compositions[phase]});
        weight[phase] =
            mole_fractions[phase] /
            properties
                .molar_density_mol_per_m3;
        sum += weight[phase];
    }
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        const auto column =
            layout
                .independent_saturation_unknown_index(
                    static_cast<
                        flow::PhaseSlot3>(
                            phase));
        require_sample6_transaction(
            column.has_value(),
            "missing Sample-6 3P saturation column");
        q[*column] =
            weight[phase] / sum;
    }

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<
                flow::PhaseSlot3>(
                    phase);
        for (std::size_t component = 0U;
             component <
                 compositions[phase].size();
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (column.has_value()) {
                q[*column] =
                    compositions[phase][component];
            }
        }
    }
    return q;
}

struct Sample6TargetLinearDiagnostics {
    std::size_t rank{};
    std::size_t component_count{};
    std::size_t phase_count{};
    double minimum_pivot{};
    double residual_l2{};
    double row_scaled_residual_l2{};
    double component_residual_l2{};
    double energy_residual_abs{};
    double fugacity_residual_l2{};
    double component_row_scaled_l2{};
    double energy_row_scaled_abs{};
    double fugacity_row_scaled_l2{};
    std::vector<double> raw_rows;
    std::vector<double> row_scaled_rows;
};

struct Sample6TargetStateConsistencyDiagnostics {
    bool checked_three_phase{};
    double sidecar_q_max_relative_difference{};
    double evaluator_sidecar_max_relative_difference{};
    double component_history_max_relative_difference{};
    double energy_history_relative_difference{};
    double fugacity_max_abs{};
};

Sample6TargetLinearDiagnostics
sample6_target_linear_diagnostics(
    fdp::PhaseTransitionRebuiltNaturalVariableSystem3D&
        system) {
    const auto& numbering =
        system.numbering();
    const PetscInt start =
        numbering.petsc_owned_scalar_start();
    const PetscInt end =
        numbering.petsc_owned_scalar_end();
    require_sample6_transaction(
        start >= 0 &&
            end > start &&
            numbering.rank_count() == 1U,
        "Sample-6 target rank diagnostic requires PETSC_COMM_SELF");

    const std::size_t n =
        static_cast<std::size_t>(
            end - start);
    std::vector<PetscInt> indices(n);
    for (std::size_t i = 0U; i < n; ++i) {
        indices[i] =
            start +
            static_cast<PetscInt>(i);
    }

    Vec residual = nullptr;
    Mat jacobian = nullptr;
    require_sample6_transaction(
        VecDuplicate(
            system.initial_state(),
            &residual) ==
                PETSC_SUCCESS &&
            MatDuplicate(
                system.jacobian_structure(),
                MAT_DO_NOT_COPY_VALUES,
                &jacobian) ==
                PETSC_SUCCESS,
        "failed to allocate Sample-6 target rank diagnostic objects");

    auto cleanup =
        [&]() {
            if (jacobian != nullptr) {
                (void)MatDestroy(&jacobian);
            }
            if (residual != nullptr) {
                (void)VecDestroy(&residual);
            }
        };

    try {
        require_sample6_transaction(
            VecSet(residual, 0.0) ==
                PETSC_SUCCESS &&
            MatZeroEntries(jacobian) ==
                PETSC_SUCCESS,
            "failed to clear Sample-6 target diagnostic objects");

        auto evaluator =
            system.snes_evaluator();
        fdp::NaturalVariableSnesEvaluationStatus3D
            status =
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success;
        require_sample6_transaction(
            evaluator.function(
                system.initial_state(),
                residual,
                evaluator.user_context,
                &status) ==
                PETSC_SUCCESS &&
            status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            VecAssemblyBegin(residual) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(residual) ==
                PETSC_SUCCESS,
            "failed to evaluate Sample-6 target residual diagnostic");

        status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
        require_sample6_transaction(
            evaluator.jacobian(
                system.initial_state(),
                jacobian,
                evaluator.user_context,
                &status) ==
                PETSC_SUCCESS &&
            status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            MatAssemblyBegin(
                jacobian,
                MAT_FINAL_ASSEMBLY) ==
                PETSC_SUCCESS &&
            MatAssemblyEnd(
                jacobian,
                MAT_FINAL_ASSEMBLY) ==
                PETSC_SUCCESS,
            "failed to evaluate Sample-6 target Jacobian diagnostic");

        PetscReal residual_norm = 0.0;
        require_sample6_transaction(
            VecNorm(
                residual,
                NORM_2,
                &residual_norm) ==
                PETSC_SUCCESS,
            "failed to norm Sample-6 target residual");

        std::vector<PetscScalar>
            matrix_values(n * n);
        std::vector<PetscScalar>
            residual_values(n);
        require_sample6_transaction(
            MatGetValues(
                jacobian,
                static_cast<PetscInt>(n),
                indices.data(),
                static_cast<PetscInt>(n),
                indices.data(),
                matrix_values.data()) ==
                PETSC_SUCCESS &&
            VecGetValues(
                residual,
                static_cast<PetscInt>(n),
                indices.data(),
                residual_values.data()) ==
                PETSC_SUCCESS,
            "failed to read Sample-6 target dense diagnostic");

        const std::size_t component_count =
            numbering.component_count();
        const std::size_t phase_count =
            numbering
                .cell(mesh::LocalIndex{0U})
                .phase_count;
        require_sample6_transaction(
            n ==
                phase_count *
                    component_count +
                    1U,
            "Sample-6 target equation/cardinality diagnostic mismatch");

        std::vector<double> a(n * n);
        std::vector<double> raw_r(n);
        std::vector<double> scaled_r(n);
        for (std::size_t row = 0U;
             row < n;
             ++row) {
            double maximum = 0.0;
            for (std::size_t column = 0U;
                 column < n;
                 ++column) {
                const double value =
                    static_cast<double>(
                        PetscRealPart(
                            matrix_values[
                                row * n +
                                column]));
                require_sample6_transaction(
                    std::isfinite(value),
                    "non-finite Sample-6 target Jacobian diagnostic");
                a[row * n + column] =
                    value;
                maximum =
                    std::max(
                        maximum,
                        std::abs(value));
            }
            require_sample6_transaction(
                maximum > 0.0 &&
                    std::isfinite(maximum),
                "zero Sample-6 target Jacobian row");
            for (std::size_t column = 0U;
                 column < n;
                 ++column) {
                a[row * n + column] /=
                    maximum;
            }
            raw_r[row] =
                static_cast<double>(
                    PetscRealPart(
                        residual_values[row]));
            scaled_r[row] =
                raw_r[row] /
                maximum;
        }

        auto block_l2 =
            [](std::span<const double> values) {
                double sum2 = 0.0;
                for (const double value : values) {
                    sum2 += value * value;
                }
                return std::sqrt(sum2);
            };

        const auto component_raw =
            std::span<const double>{
                raw_r.data(),
                component_count};
        const auto fugacity_raw =
            std::span<const double>{
                raw_r.data() +
                    component_count +
                    1U,
                n -
                    component_count -
                    1U};
        const auto component_scaled =
            std::span<const double>{
                scaled_r.data(),
                component_count};
        const auto fugacity_scaled =
            std::span<const double>{
                scaled_r.data() +
                    component_count +
                    1U,
                n -
                    component_count -
                    1U};

        double scaled_norm2 = 0.0;
        for (double value : scaled_r) {
            scaled_norm2 +=
                value * value;
        }

        std::size_t rank = 0U;
        double minimum_pivot =
            std::numeric_limits<double>::
                infinity();
        constexpr double pivot_tolerance =
            1.0e-11;
        for (std::size_t column = 0U;
             column < n &&
             rank < n;
             ++column) {
            std::size_t pivot_row =
                rank;
            double pivot_magnitude = 0.0;
            for (std::size_t row = rank;
                 row < n;
                 ++row) {
                const double magnitude =
                    std::abs(
                        a[row * n + column]);
                if (magnitude >
                    pivot_magnitude) {
                    pivot_magnitude =
                        magnitude;
                    pivot_row = row;
                }
            }
            if (!(pivot_magnitude >
                  pivot_tolerance)) {
                continue;
            }
            if (pivot_row != rank) {
                for (std::size_t j = 0U;
                     j < n;
                     ++j) {
                    std::swap(
                        a[rank * n + j],
                        a[pivot_row * n + j]);
                }
            }
            const double pivot =
                a[rank * n + column];
            minimum_pivot =
                std::min(
                    minimum_pivot,
                    std::abs(pivot));
            for (std::size_t row =
                     rank + 1U;
                 row < n;
                 ++row) {
                const double factor =
                    a[row * n + column] /
                    pivot;
                if (factor == 0.0) {
                    continue;
                }
                for (std::size_t j = column;
                     j < n;
                     ++j) {
                    a[row * n + j] -=
                        factor *
                        a[rank * n + j];
                }
            }
            ++rank;
        }

        Sample6TargetLinearDiagnostics result;
        result.rank = rank;
        result.component_count =
            component_count;
        result.phase_count =
            phase_count;
        result.minimum_pivot =
            rank == 0U
                ? 0.0
                : minimum_pivot;
        result.residual_l2 =
            static_cast<double>(
                residual_norm);
        result.row_scaled_residual_l2 =
            std::sqrt(
                scaled_norm2);
        result.component_residual_l2 =
            block_l2(component_raw);
        result.energy_residual_abs =
            std::abs(
                raw_r[component_count]);
        result.fugacity_residual_l2 =
            block_l2(fugacity_raw);
        result.component_row_scaled_l2 =
            block_l2(component_scaled);
        result.energy_row_scaled_abs =
            std::abs(
                scaled_r[component_count]);
        result.fugacity_row_scaled_l2 =
            block_l2(fugacity_scaled);
        result.raw_rows =
            std::move(raw_r);
        result.row_scaled_rows =
            std::move(scaled_r);

        cleanup();
        return result;
    } catch (...) {
        cleanup();
        throw;
    }
}

template <typename Provider>
Sample6TargetStateConsistencyDiagnostics
sample6_target_state_consistency_diagnostics(
    fdp::PhaseTransitionRebuiltNaturalVariableSystem3D&
        system,
    fdp::Sw92TransactionalPhaseTransitionRebuildContext3D<
        Provider>& context) {
    Sample6TargetStateConsistencyDiagnostics
        result;

    const auto& numbering =
        system.numbering();
    const auto& record =
        numbering.cell(mesh::LocalIndex{0U});
    const auto* target =
        context.scanner_context == nullptr
        ? nullptr
        : context.scanner_context
              ->find_target(
                  record.cell_global);
    require_sample6_transaction(
        target != nullptr,
        "Sample-6 rebuilt target lost authoritative sidecar");

    std::vector<PetscInt>
        indices(record.scalar_count);
    std::vector<PetscScalar>
        petsc_values(record.scalar_count);
    for (std::size_t slot = 0U;
         slot < record.scalar_count;
         ++slot) {
        indices[slot] =
            record.petsc_global_scalar_start +
            static_cast<PetscInt>(slot);
    }
    require_sample6_transaction(
        VecGetValues(
            system.initial_state(),
            static_cast<PetscInt>(
                indices.size()),
            indices.data(),
            petsc_values.data()) ==
            PETSC_SUCCESS,
        "failed to read Sample-6 rebuilt initial q");

    std::vector<double> q(
        record.scalar_count);
    for (std::size_t slot = 0U;
         slot < record.scalar_count;
         ++slot) {
        q[slot] =
            static_cast<double>(
                PetscRealPart(
                    petsc_values[slot]));
    }

    const auto authoritative_q =
        target->projection
            .natural_variables();
    require_sample6_transaction(
        q.size() ==
            authoritative_q.size(),
        "Sample-6 rebuilt q/sidecar shape mismatch");

    auto relative_difference =
        [](double first,
           double second) {
            const double scale =
                std::max(
                    {1.0,
                     std::abs(first),
                     std::abs(second)});
            return std::abs(
                       first - second) /
                scale;
        };

    for (std::size_t column = 0U;
         column < q.size();
         ++column) {
        result.sidecar_q_max_relative_difference =
            std::max(
                result
                    .sidecar_q_max_relative_difference,
                relative_difference(
                    q[column],
                    authoritative_q[column]));
    }

    if (record.phase_count != 3U) {
        return result;
    }

    std::optional<
        fdp::MixedCardinalityPhysicalCurrentCellLinearization3D>
        current;
    double porosity = 0.0;
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    require_sample6_transaction(
        system
                .evaluate_current_cell_for_phase_transition(
                    record.cell,
                    q,
                    &current,
                    &porosity,
                    &status) ==
                PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            current.has_value(),
        "Sample-6 rebuilt evaluator cannot reproduce target q");

    const auto* three =
        std::get_if<
            fdp::FixedThreePhaseCurrentCellLinearization3D>(
                &*current);
    require_sample6_transaction(
        three != nullptr &&
            target->phases.size() == 3U &&
            target->projection
                    .target_saturations()
                    .size() == 3U,
        "Sample-6 rebuilt evaluator/sidecar is not authoritative 3P");

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<
                flow::PhaseSlot3>(
                    phase);
        result.evaluator_sidecar_max_relative_difference =
            std::max(
                result
                    .evaluator_sidecar_max_relative_difference,
                relative_difference(
                    three->state
                        .phase_saturation(slot),
                    target->projection
                        .target_saturations()[
                            phase]));
        result.evaluator_sidecar_max_relative_difference =
            std::max(
                result
                    .evaluator_sidecar_max_relative_difference,
                relative_difference(
                    three->molar_density
                        .molar_density_mol_per_m3[
                            phase],
                    target->phases[phase]
                        .molar_density_mol_per_m3));

        const auto composition =
            three->state
                .phase_composition(slot);
        require_sample6_transaction(
            composition.size() ==
                target->phases[phase]
                    .composition.size(),
            "Sample-6 rebuilt evaluator/sidecar composition shape mismatch");
        for (std::size_t component = 0U;
             component <
                 composition.size();
             ++component) {
            result.evaluator_sidecar_max_relative_difference =
                std::max(
                    result
                        .evaluator_sidecar_max_relative_difference,
                    relative_difference(
                        composition[component],
                        target->phases[phase]
                            .composition[
                                component]));
        }
    }

    require_sample6_transaction(
        context.baseline_cells.size() == 1U &&
            context.baseline_cells[0]
                .previous_component_accumulation
                .has_value() &&
            context.baseline_cells[0]
                .previous_energy_accumulation
                .has_value(),
        "Sample-6 rebuilt consistency diagnostic lacks frozen history");

    const auto current_component =
        flow::build_pore_volume_component_accumulation(
            three->state,
            porosity);
    const auto& previous_component =
        *context.baseline_cells[0]
             .previous_component_accumulation;
    require_sample6_transaction(
        current_component.component_ids ==
                previous_component.component_ids &&
            current_component
                    .component_accumulation_mol_per_bulk_m3
                    .size() ==
                previous_component
                    .component_accumulation_mol_per_bulk_m3
                    .size(),
        "Sample-6 rebuilt component history identity mismatch");
    for (std::size_t component = 0U;
         component <
             current_component
                 .component_accumulation_mol_per_bulk_m3
                 .size();
         ++component) {
        result.component_history_max_relative_difference =
            std::max(
                result
                    .component_history_max_relative_difference,
                relative_difference(
                    current_component
                        .component_accumulation_mol_per_bulk_m3[
                            component],
                    previous_component
                        .component_accumulation_mol_per_bulk_m3[
                            component]));
    }

    const auto current_energy =
        flow::build_pore_volume_energy_accumulation_snapshot(
            three->state,
            porosity,
            three->transport,
            three->caloric,
            three->rock);
    const auto& previous_energy =
        *context.baseline_cells[0]
             .previous_energy_accumulation;
    result.energy_history_relative_difference =
        relative_difference(
            current_energy
                .total_internal_energy_j_per_bulk_m3,
            previous_energy
                .total_internal_energy_j_per_bulk_m3);

    for (const double value :
         three->fugacity.values()) {
        result.fugacity_max_abs =
            std::max(
                result.fugacity_max_abs,
                std::abs(value));
    }
    result.checked_three_phase = true;
    return result;
}

template <typename Provider>
PetscErrorCode sample6_preflight_rebuild(
    const fdp::PhaseTransitionRebuiltNaturalVariableSystem3D&
        current_system,
    Vec converged_state,
    const fdp::VariableCardinalityNaturalVariableSnesSolveReport3D&
        solve_report,
    std::span<
        const fdp::PostSnesPhaseTransitionProposal3D>
        local_owned_proposals,
    std::span<
        const fdp::AcceptedPhaseTransitionSummary3D>
        accepted_global_batch,
    void* raw_context,
    std::unique_ptr<
        fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>*
            rebuilt_system) {
    const PetscErrorCode rebuild_error =
        fdp::
            rebuild_sw92_transactional_phase_transition_system_3d<
                Provider>(
                current_system,
                converged_state,
                solve_report,
                local_owned_proposals,
                accepted_global_batch,
                raw_context,
                rebuilt_system);
    if (rebuild_error != PETSC_SUCCESS ||
        rebuilt_system == nullptr ||
        *rebuilt_system == nullptr) {
        return rebuild_error != PETSC_SUCCESS
            ? rebuild_error
            : PETSC_ERR_PLIB;
    }

    Vec probe_state = nullptr;
    std::optional<
        fdp::VariableCardinalityNaturalVariableSnesSolveReport3D>
        probe_report;
    std::optional<
        fdp::NaturalVariableSnesFailureDiagnostics3D>
        diagnostics;
    const auto linear_diagnostics =
        sample6_target_linear_diagnostics(
            **rebuilt_system);
    auto* rebuild_context =
        static_cast<
            fdp::Sw92TransactionalPhaseTransitionRebuildContext3D<
                Provider>*>(
                    raw_context);
    const auto state_diagnostics =
        sample6_target_state_consistency_diagnostics(
            **rebuilt_system,
            *rebuild_context);
    const PetscErrorCode solve_error =
        (*rebuilt_system)
            ->solve(
                &probe_state,
                &probe_report,
                &diagnostics);
    if (probe_state != nullptr) {
        (void)VecDestroy(&probe_state);
    }
    if (solve_error != PETSC_SUCCESS) {
        std::cerr
            << "[Sample-6 target preflight] error="
            << static_cast<int>(solve_error)
            << " dense_rank="
            << linear_diagnostics.rank
            << " dense_size="
            << (*rebuilt_system)
                   ->numbering()
                   .petsc_local_owned_scalar_count()
            << " min_scaled_pivot="
            << linear_diagnostics.minimum_pivot
            << " raw_residual_l2="
            << linear_diagnostics.residual_l2
            << " row_scaled_residual_l2="
            << linear_diagnostics.row_scaled_residual_l2
            << " component_raw_l2="
            << linear_diagnostics.component_residual_l2
            << " energy_raw_abs="
            << linear_diagnostics.energy_residual_abs
            << " fugacity_raw_l2="
            << linear_diagnostics.fugacity_residual_l2
            << " component_scaled_l2="
            << linear_diagnostics.component_row_scaled_l2
            << " energy_scaled_abs="
            << linear_diagnostics.energy_row_scaled_abs
            << " fugacity_scaled_l2="
            << linear_diagnostics.fugacity_row_scaled_l2
            << " sidecar_q_max_rel_diff="
            << state_diagnostics
                   .sidecar_q_max_relative_difference
            << " evaluator_sidecar_max_rel_diff="
            << state_diagnostics
                   .evaluator_sidecar_max_relative_difference
            << " component_history_max_rel_diff="
            << state_diagnostics
                   .component_history_max_relative_difference
            << " energy_history_rel_diff="
            << state_diagnostics
                   .energy_history_relative_difference
            << " fugacity_max_abs="
            << state_diagnostics.fugacity_max_abs;
        if (diagnostics.has_value()) {
            std::cerr
                << " snes_reason="
                << static_cast<int>(
                       diagnostics->snes_reason)
                << " ksp_reason="
                << static_cast<int>(
                       diagnostics->ksp_reason)
                << " pc_reason="
                << diagnostics->pc_failed_reason
                << " sub_ksp_reason="
                << static_cast<int>(
                       diagnostics->asm_sub_ksp_reason)
                << " sub_pc_reason="
                << diagnostics->asm_sub_pc_failed_reason
                << " nonlinear_iterations="
                << static_cast<long long>(
                       diagnostics->nonlinear_iterations)
                << " function_evaluations="
                << static_cast<long long>(
                       diagnostics->function_evaluations)
                << " jacobian_evaluations="
                << static_cast<long long>(
                       diagnostics->jacobian_evaluations)
                << " function_domain_errors="
                << static_cast<long long>(
                       diagnostics->function_domain_errors)
                << " jacobian_domain_errors="
                << static_cast<long long>(
                       diagnostics->jacobian_domain_errors)
                << " line_search_prechecks="
                << static_cast<long long>(
                       diagnostics->line_search_prechecks)
                << " line_search_direction_changes="
                << static_cast<long long>(
                       diagnostics->line_search_direction_changes)
                << " function_l2_norm="
                << diagnostics->function_l2_norm;
        }
        std::cerr << " raw_rows=[";
        for (std::size_t row = 0U;
             row <
                 linear_diagnostics.raw_rows.size();
             ++row) {
            if (row != 0U) {
                std::cerr << ',';
            }
            std::cerr <<
                linear_diagnostics.raw_rows[row];
        }
        std::cerr << "] scaled_rows=[";
        for (std::size_t row = 0U;
             row <
                 linear_diagnostics
                     .row_scaled_rows.size();
             ++row) {
            if (row != 0U) {
                std::cerr << ',';
            }
            std::cerr <<
                linear_diagnostics
                    .row_scaled_rows[row];
        }
        std::cerr << "]\n";
        return solve_error;
    }
    return PETSC_SUCCESS;
}

template <typename Provider>
void run_controller(
    std::unique_ptr<
        fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
        initial_system,
    fdp::PostSnesSw92ProfileCPhaseTransitionScannerContext3D*
        scanner,
    fdp::Sw92TransactionalPhaseTransitionRebuildContext3D<
        Provider>* rebuild,
    std::size_t source_count,
    std::size_t target_count) {
    std::unique_ptr<
        fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
        final_system;
    Vec final_state = nullptr;
    std::optional<
        fdp::PostSnesPhaseTransitionControllerReport3D>
        report;

    const PetscErrorCode error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_SELF,
                std::move(initial_system),
                {
                    &fdp::
                        scan_post_snes_sw92_profile_c_phase_transitions_3d,
                    scanner,
                    &sample6_preflight_rebuild<
                        Provider>,
                    rebuild},
                {.max_transition_restarts = 4U},
                &final_system,
                &final_state,
                &report);

    std::string lifecycle_diagnostic =
        "Sample-6 SW92 " +
        std::to_string(source_count) +
        "P->" +
        std::to_string(target_count) +
        "P lifecycle failed: error=" +
        std::to_string(static_cast<int>(error)) +
        " final_system=" +
        std::to_string(final_system != nullptr ? 1 : 0) +
        " final_state=" +
        std::to_string(final_state != nullptr ? 1 : 0) +
        " report=" +
        std::to_string(report.has_value() ? 1 : 0) +
        " rebuild_runtimes=" +
        std::to_string(rebuild->runtimes.size()) +
        " scanner_targets=" +
        std::to_string(scanner->local_owned_targets.size());
    if (report.has_value()) {
        lifecycle_diagnostic +=
            " outcome=" +
            std::to_string(
                static_cast<int>(
                    report->outcome)) +
            " restarts=" +
            std::to_string(
                report->transition_restarts) +
            " generations=" +
            std::to_string(
                report->generations.size());
        if (!report->generations.empty()) {
            lifecycle_diagnostic +=
                " gen0_batch=" +
                std::to_string(
                    report->generations[0]
                        .accepted_transition_batch
                        .size()) +
                " gen0_snes_reason=" +
                std::to_string(
                    static_cast<int>(
                        report->generations[0]
                            .nonlinear_solve
                            .converged_reason));
        }
    }
    if (final_system != nullptr) {
        lifecycle_diagnostic +=
            " final_phase_count=" +
            std::to_string(
                final_system->numbering()
                    .cell(mesh::LocalIndex{0U})
                    .phase_count) +
            " final_dt=" +
            std::to_string(
                final_system
                    ->time_step_seconds());
    }

    require_sample6_transaction(
        error == PETSC_SUCCESS &&
            final_system != nullptr &&
            final_state != nullptr &&
            report.has_value() &&
            report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            report->transition_restarts == 1U &&
            report->generations.size() == 2U &&
            report->generations[0]
                    .accepted_transition_batch
                    .size() == 1U &&
            report->generations[0]
                    .accepted_transition_batch[0]
                    .source_phase_count ==
                source_count &&
            report->generations[0]
                    .accepted_transition_batch[0]
                    .target_phase_count ==
                target_count &&
            report->generations[1]
                    .accepted_transition_batch
                    .empty() &&
            final_system->numbering()
                    .cell(mesh::LocalIndex{0U})
                    .phase_count ==
                target_count &&
            final_system->time_step_seconds() == 1.0,
        lifecycle_diagnostic);

    bool history_matches = false;
    require_sample6_transaction(
        final_system
                ->accepted_history_matches_state(
                    final_state,
                    &history_matches) ==
                PETSC_SUCCESS &&
            history_matches,
        "Sample-6 SW92 restart changed frozen component/energy history");

    (void)VecDestroy(&final_state);
}

void sample6_two_to_three() {
    auto parameters =
        sample6_transaction_parameters();
    auto model =
        th::Sw92Phase<double>::
            from_parameters(parameters);

    constexpr double pressure_pa = 1.0e7;
    constexpr double temperature_k = 350.0;
    const auto feed =
        sample6::feed();
    const auto owned =
        fl::solve_sw92_phase_assigned_pt_boundary_aware(
            pressure_pa,
            temperature_k,
            feed,
            model,
            0.0);
    require_sample6_transaction(
        owned.base.c1.has_value() &&
            owned.base.c1->candidate() != nullptr,
        "Sample-6 C1 2P source is unavailable");
    const auto& c1 =
        *owned.base.c1->candidate();

    std::array<std::vector<double>, 2>
        compositions{
            c1.aqueous_phase.composition,
            c1.nonaqueous_phase.composition};
    std::array<double, 2>
        fractions{
            c1.aqueous_phase.mole_phase_fraction,
            c1.nonaqueous_phase.mole_phase_fraction};
    std::vector<
        th::Sw92SelectedPhase<double>>
        selections{
            {
                0.0,
                th::SwPhaseFamily::aqueous,
                c1.aqueous_phase.activity.branch,
                {}},
            {
                0.0,
                th::SwPhaseFamily::nonaqueous,
                c1.nonaqueous_phase.activity.branch,
                {}}};

    ManufacturedSample6TransportCaloricProvider
        provider;
    auto closure =
        flow::make_sw92_selected_phase_property_closure(
            model,
            selections,
            provider,
            manufactured_sample6_provenance());
    using Closure = decltype(closure);

    const flow::NaturalVariableLayout2P layout{
        flow::NaturalVariableCompositionPivot2P::
            from_dependent_components(
                feed.size(),
                {
                    dependent_largest(
                        compositions[0]),
                    dependent_largest(
                        compositions[1])})};
    const auto q =
        two_phase_q(
            layout,
            pressure_pa,
            temperature_k,
            compositions,
            fractions,
            closure);

    fdp::Sw92TwoPhaseProductionCellEvaluatorContext3D<
        Closure>
        source_context{
            &closure,
            {&sample6_two_phase_kr, nullptr},
            {&sample6_zero_rock, nullptr},
            {}};
    std::optional<
        fdp::TwoPhaseCurrentCellLinearization3D>
        source_current;
    fdp::NaturalVariableSnesEvaluationStatus3D
        source_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    require_sample6_transaction(
        fdp::
            evaluate_sw92_two_phase_production_cell_3d<
                Closure>(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{UINT64_C(9201)},
                q,
                layout,
                closure.component_ids(),
                &source_context,
                &source_current,
                &source_status) ==
                PETSC_SUCCESS &&
            source_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            source_current.has_value(),
        "Sample-6 2P source evaluator failed");

    constexpr double porosity = 0.20;
    const auto previous_component =
        flow::build_two_phase_component_accumulation(
            source_current->state,
            porosity);
    const auto previous_energy =
        flow::build_two_phase_energy_accumulation_snapshot(
            source_current->state,
            porosity,
            source_current->transport,
            source_current->caloric,
            source_current->rock);

    auto partition =
        sample6_partition();
    auto schedule =
        sample6_schedule();
    auto bridge =
        sample6_bridge();
    auto pattern =
        sample6_pattern();

    fdp::CellScopedMixedCardinalityEvaluatorDispatcher3D
        initial_dispatcher{
            {},
            {{
                mesh::GlobalEntityId{UINT64_C(9201)},
                {
                    &fdp::
                        evaluate_sw92_two_phase_production_cell_3d<
                            Closure>,
                    &source_context}}},
            {}};

    fdp::FrozenPhaseTransitionRebuildCell3D
        cell;
    cell.cell = mesh::LocalIndex{0U};
    cell.cell_global =
        mesh::GlobalEntityId{UINT64_C(9201)};
    cell.bulk_volume_m3 = 1.0;
    cell.porosity = porosity;
    cell.component_ids.assign(
        closure.component_ids().begin(),
        closure.component_ids().end());
    cell.target_layout = layout.descriptor();
    cell.target_natural_variables = q;
    cell.target_active_phases =
        flow::FrozenActivePhaseIdentityMap{
            {
                sample6_identity("aqueous"),
                sample6_identity("nonaqueous")}};
    cell.previous_component_accumulation =
        previous_component;
    cell.previous_energy_accumulation =
        previous_energy;
    cell.transition_evidence_profile =
        "SW92/sample6/transaction-source-2P/v1";

    int absent_token = 1;
    std::unique_ptr<
        fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
        initial_system;
    require_sample6_transaction(
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_SELF,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                {cell},
                {},
                initial_dispatcher.bindings(),
                {},
                &sample6_no_absent_provider,
                &absent_token,
                &initial_system) ==
                PETSC_SUCCESS &&
            initial_system != nullptr,
        "failed to build Sample-6 2P source system");

    fdp::PostSnesSw92ProfileCPhaseTransitionScannerContext3D
        scanner{&model, {}, {}};

    using Provider =
        ManufacturedSample6TransportCaloricProvider;
    fdp::Sw92TransactionalPhaseTransitionRebuildContext3D<
        Provider>
        rebuild;
    rebuild.comm = PETSC_COMM_SELF;
    rebuild.schedule = &schedule;
    rebuild.partition = &partition;
    rebuild.cell_bridge = &bridge;
    rebuild.cell_pattern = &pattern;
    rebuild.model = &model;
    rebuild.provider = provider;
    rebuild.provenance =
        manufactured_sample6_provenance();
    rebuild.scanner_context = &scanner;
    rebuild.target_identity = {
        &sample6_target_identity,
        nullptr};
    rebuild.initial_evaluators =
        initial_dispatcher.bindings();
    rebuild.two_phase_relative_permeability = {
        &sample6_two_phase_kr,
        nullptr};
    rebuild.three_phase_saturation = {
        &sample6_three_phase_saturation,
        nullptr};
    rebuild.rock_storage = {
        &sample6_zero_rock,
        nullptr};
    rebuild.project_target_to_conservation_storage =
        true;
    rebuild.baseline_cells.push_back(
        {
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(9201)},
            1.0,
            porosity,
            cell.component_ids,
            previous_component,
            previous_energy});

    run_controller(
        std::move(initial_system),
        &scanner,
        &rebuild,
        2U,
        3U);
}

void sample6_three_to_two() {
    auto parameters =
        sample6_transaction_parameters();
    auto model =
        th::Sw92Phase<double>::
            from_parameters(parameters);

    constexpr double pressure_pa = 1.0e7;
    constexpr double temperature_k = 350.0;
    const auto feed =
        sample6_edge_feed(
            sample6::w(),
            sample6::h0(),
            0.45);
    const auto authoritative =
        fl::solve_sw92_profile_c_pt_phase_set(
            pressure_pa,
            temperature_k,
            feed,
            model,
            0.0);
    require_sample6_transaction(
        authoritative.solution.status ==
                fl::PtPhaseSetStatus::accepted &&
            authoritative
                .accepted_phase_set_published() &&
            authoritative.solution
                    .accepted_phase_count() == 2U &&
            authoritative.phase_metadata.size() == 2U,
        "Sample-6 W+H target is not authoritative 2P");
    const auto* accepted =
        authoritative.solution
            .accepted_phase_set();
    require_sample6_transaction(
        accepted != nullptr &&
            accepted->phases.size() == 2U,
        "Sample-6 W+H accepted phases missing");

    const auto& aq =
        accepted->phases[0];
    const auto& na =
        accepted->phases[1];
    std::array<std::vector<double>, 3>
        compositions{
            aq.composition,
            na.composition,
            na.composition};
    std::array<double, 3>
        fractions{
            aq.mole_phase_fraction,
            0.5 * na.mole_phase_fraction,
            0.5 * na.mole_phase_fraction};
    std::vector<
        th::Sw92SelectedPhase<double>>
        selections{
            {
                0.0,
                th::SwPhaseFamily::aqueous,
                aq.activity.branch,
                {}},
            {
                0.0,
                th::SwPhaseFamily::nonaqueous,
                na.activity.branch,
                {}},
            {
                0.0,
                th::SwPhaseFamily::nonaqueous,
                na.activity.branch,
                {}}};

    ManufacturedSample6TransportCaloricProvider
        provider;
    auto closure =
        flow::make_sw92_selected_phase_property_closure(
            model,
            selections,
            provider,
            manufactured_sample6_provenance());
    using Closure = decltype(closure);

    const flow::NaturalVariableLayout3P layout{
        flow::NaturalVariableCompositionPivot3P::
            from_dependent_components(
                feed.size(),
                {
                    dependent_largest(
                        compositions[0]),
                    dependent_largest(
                        compositions[1]),
                    dependent_largest(
                        compositions[2])})};
    const auto q =
        three_phase_q(
            layout,
            pressure_pa,
            temperature_k,
            compositions,
            fractions,
            closure);

    fdp::Sw92ThreePhaseProductionCellEvaluatorContext3D<
        Closure>
        source_context{
            &closure,
            {&sample6_three_phase_saturation, nullptr},
            {&sample6_zero_rock, nullptr},
            {}};
    std::optional<
        fdp::FixedThreePhaseCurrentCellLinearization3D>
        source_current;
    fdp::NaturalVariableSnesEvaluationStatus3D
        source_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    require_sample6_transaction(
        fdp::
            evaluate_sw92_three_phase_production_cell_3d<
                Closure>(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{UINT64_C(9201)},
                q,
                layout,
                closure.component_ids(),
                &source_context,
                &source_current,
                &source_status) ==
                PETSC_SUCCESS &&
            source_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            source_current.has_value(),
        "Sample-6 stale 3P source evaluator failed");

    constexpr double porosity = 0.20;
    const auto previous_component =
        flow::build_pore_volume_component_accumulation(
            source_current->state,
            porosity);
    const auto previous_energy =
        flow::build_pore_volume_energy_accumulation_snapshot(
            source_current->state,
            porosity,
            source_current->transport,
            source_current->caloric,
            source_current->rock);

    auto partition =
        sample6_partition();
    auto schedule =
        sample6_schedule();
    auto bridge =
        sample6_bridge();
    auto pattern =
        sample6_pattern();

    fdp::CellScopedMixedCardinalityEvaluatorDispatcher3D
        initial_dispatcher{
            {},
            {},
            {{
                mesh::GlobalEntityId{UINT64_C(9201)},
                {
                    &fdp::
                        evaluate_sw92_three_phase_production_cell_3d<
                            Closure>,
                    &source_context}}}};

    fdp::FrozenPhaseTransitionRebuildCell3D
        cell;
    cell.cell = mesh::LocalIndex{0U};
    cell.cell_global =
        mesh::GlobalEntityId{UINT64_C(9201)};
    cell.bulk_volume_m3 = 1.0;
    cell.porosity = porosity;
    cell.component_ids.assign(
        closure.component_ids().begin(),
        closure.component_ids().end());
    cell.target_layout =
        flow::NaturalVariableLayoutDescriptor{
            layout};
    cell.target_natural_variables = q;
    cell.target_active_phases =
        flow::FrozenActivePhaseIdentityMap{
            {
                sample6_identity("aqueous"),
                sample6_identity("nonaqueous"),
                sample6_identity("extra")}};
    cell.previous_component_accumulation =
        previous_component;
    cell.previous_energy_accumulation =
        previous_energy;
    cell.transition_evidence_profile =
        "SW92/sample6/transaction-source-3P/v1";

    int absent_token = 1;
    std::unique_ptr<
        fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
        initial_system;
    require_sample6_transaction(
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_SELF,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                {cell},
                {},
                initial_dispatcher.bindings(),
                {},
                &sample6_no_absent_provider,
                &absent_token,
                &initial_system) ==
                PETSC_SUCCESS &&
            initial_system != nullptr,
        "failed to build Sample-6 stale 3P source system");

    fdp::PostSnesSw92ProfileCPhaseTransitionScannerContext3D
        scanner{&model, {}, {}};

    using Provider =
        ManufacturedSample6TransportCaloricProvider;
    fdp::Sw92TransactionalPhaseTransitionRebuildContext3D<
        Provider>
        rebuild;
    rebuild.comm = PETSC_COMM_SELF;
    rebuild.schedule = &schedule;
    rebuild.partition = &partition;
    rebuild.cell_bridge = &bridge;
    rebuild.cell_pattern = &pattern;
    rebuild.model = &model;
    rebuild.provider = provider;
    rebuild.provenance =
        manufactured_sample6_provenance();
    rebuild.scanner_context = &scanner;
    rebuild.target_identity = {
        &sample6_target_identity,
        nullptr};
    rebuild.initial_evaluators =
        initial_dispatcher.bindings();
    rebuild.two_phase_relative_permeability = {
        &sample6_two_phase_kr,
        nullptr};
    rebuild.three_phase_saturation = {
        &sample6_three_phase_saturation,
        nullptr};
    rebuild.rock_storage = {
        &sample6_zero_rock,
        nullptr};
    rebuild.project_target_to_conservation_storage =
        true;
    rebuild.baseline_cells.push_back(
        {
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(9201)},
            1.0,
            porosity,
            cell.component_ids,
            previous_component,
            previous_energy});

    run_controller(
        std::move(initial_system),
        &scanner,
        &rebuild,
        3U,
        2U);
}

} // namespace

void sw92_transactional_phase_transition_sample6_test() {
    sample6_two_to_three();
    sample6_three_to_two();
}
