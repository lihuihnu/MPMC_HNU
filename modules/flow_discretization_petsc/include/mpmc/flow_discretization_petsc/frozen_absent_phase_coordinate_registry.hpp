#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_FROZEN_ABSENT_PHASE_COORDINATE_REGISTRY_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_FROZEN_ABSENT_PHASE_COORDINATE_REGISTRY_HPP

#include <mpmc/flow_discretization_petsc/thermodynamic_cross_cardinality_tpfa_adapter.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    frozen_absent_phase_coordinate_registry_convention =
        "flow_discretization_petsc/frozen-absent-phase-coordinate-registry/v1";

struct FrozenAbsentPhaseCoordinateEntry3D {
    mpmc::flow::FrozenPhysicalPhaseIdentity identity;
    mpmc::flow::AbsentPhaseThermodynamicCoordinateExtension
        reference_coordinates;
    std::string selected_branch_provenance;
};

struct FrozenAbsentPhaseCoordinateCell3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::flow::FrozenActivePhaseIdentityMap
        active_phases{
            {
                mpmc::flow::FrozenPhysicalPhaseIdentity{
                    "invalid",
                    "invalid"}}};
    std::vector<
        FrozenAbsentPhaseCoordinateEntry3D>
        absent_phases;
};

struct FrozenAbsentPhaseFaceEndpoint3D {
    mpmc::mesh::GlobalEntityId face_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::mesh::GlobalEntityId owner_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    mpmc::mesh::GlobalEntityId neighbour_cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
};

namespace frozen_absent_phase_registry_detail {

[[nodiscard]] inline bool
near_roundoff(
    double first,
    double second) {
    if (!std::isfinite(first) ||
        !std::isfinite(second)) {
        return false;
    }
    const double scale =
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second)});
    return std::abs(first - second) <=
        8192.0 *
            std::numeric_limits<double>::
                epsilon() *
            scale;
}

[[nodiscard]] inline std::vector<double>
natural_variables_from_state_identity(
    const mpmc::flow::
        NaturalVariableStateIdentity3P&
            state) {
    const auto& layout =
        state.layout;
    const std::size_t q =
        layout.unknown_count();
    const std::size_t n =
        layout.component_count();
    const std::size_t p =
        layout.phase_count();
    if (state.component_ids.size() != n ||
        p == 0U ||
        p >
            mpmc::flow::
                fixed_three_phase_count ||
        !std::isfinite(
            state.reference_pressure_pa) ||
        !(state.reference_pressure_pa > 0.0) ||
        !std::isfinite(
            state.temperature_k) ||
        !(state.temperature_k > 0.0)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: malformed natural-variable state identity in frozen coordinate registry");
    }

    std::vector<double> result(
        q,
        0.0);
    result[
        layout.pressure_unknown_index()] =
        state.reference_pressure_pa;
    result[
        layout.temperature_unknown_index()] =
        state.temperature_k;

    double saturation_sum = 0.0;
    for (std::size_t phase = 0U;
         phase < p;
         ++phase) {
        const double saturation =
            state.saturation[phase];
        if (!std::isfinite(saturation) ||
            !(saturation > 0.0)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: frozen registry host saturation is invalid");
        }
        saturation_sum +=
            saturation;
        const auto slot =
            static_cast<
                mpmc::flow::PhaseSlot3>(
                    phase);
        if (const auto column =
                layout
                    .independent_saturation_unknown_index(
                        slot)) {
            result[*column] =
                saturation;
        }

        if (state.phase_composition[phase]
                .size() != n) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: frozen registry host composition shape mismatch");
        }
        double composition_sum = 0.0;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double value =
                state.phase_composition[
                    phase][component];
            if (!std::isfinite(value) ||
                !(value > 0.0)) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: frozen registry host composition left positive support");
            }
            composition_sum += value;
            if (const auto column =
                    layout
                        .independent_composition_unknown_index(
                            slot,
                            component)) {
                result[*column] =
                    value;
            }
        }
        const double composition_tolerance =
            8192.0 *
            std::numeric_limits<double>::
                epsilon() *
            static_cast<double>(n);
        if (!std::isfinite(
                composition_sum) ||
            std::abs(
                composition_sum - 1.0) >
                composition_tolerance) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: frozen registry host composition is not normalized");
        }
    }

    const double saturation_tolerance =
        8192.0 *
        std::numeric_limits<double>::
            epsilon() *
        static_cast<double>(p);
    if (!std::isfinite(
            saturation_sum) ||
        std::abs(
            saturation_sum - 1.0) >
            saturation_tolerance) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: frozen registry host saturations are not normalized");
    }

    return result;
}

[[nodiscard]] inline bool
same_layout(
    const mpmc::flow::
        NaturalVariableLayoutDescriptor& first,
    const mpmc::flow::
        NaturalVariableLayoutDescriptor& second) {
    if (first.component_count() !=
            second.component_count() ||
        first.phase_count() !=
            second.phase_count() ||
        first.unknown_count() !=
            second.unknown_count()) {
        return false;
    }
    for (std::size_t phase = 0U;
         phase < first.phase_count();
         ++phase) {
        const auto slot =
            static_cast<
                mpmc::flow::PhaseSlot3>(
                    phase);
        if (first
                .dependent_composition_component(
                    slot) !=
            second
                .dependent_composition_component(
                    slot)) {
            return false;
        }
    }
    return true;
}

} // namespace frozen_absent_phase_registry_detail

class FrozenAbsentPhaseCoordinateRegistry3D {
public:
    FrozenAbsentPhaseCoordinateRegistry3D(
        std::vector<
            FrozenAbsentPhaseCoordinateCell3D>
            cells,
        std::vector<
            FrozenAbsentPhaseFaceEndpoint3D>
            faces)
        : cells_(
              std::move(cells)),
          faces_(
              std::move(faces)) {
        validate();
    }

    [[nodiscard]] std::span<
        const FrozenAbsentPhaseCoordinateCell3D>
    cells() const noexcept {
        return cells_;
    }

    [[nodiscard]] std::span<
        const FrozenAbsentPhaseFaceEndpoint3D>
    faces() const noexcept {
        return faces_;
    }

    [[nodiscard]]
    const FrozenAbsentPhaseCoordinateCell3D&
    cell(
        mpmc::mesh::GlobalEntityId
            cell_global) const {
        const auto found =
            std::lower_bound(
                cells_.begin(),
                cells_.end(),
                cell_global.value(),
                [](const auto& entry,
                   std::uint64_t value) {
                    return entry
                               .cell_global
                               .value() <
                        value;
                });
        if (found ==
                cells_.end() ||
            found->cell_global !=
                cell_global) {
            throw std::out_of_range(
                "mpmc::flow_discretization_petsc: frozen coordinate registry does not contain requested cell");
        }
        return *found;
    }

    [[nodiscard]]
    mpmc::mesh::GlobalEntityId
    face_cell(
        mpmc::mesh::GlobalEntityId
            face_global,
        CrossCardinalityAbsentPhaseSide3D
            side) const {
        const auto found =
            std::lower_bound(
                faces_.begin(),
                faces_.end(),
                face_global.value(),
                [](const auto& entry,
                   std::uint64_t value) {
                    return entry
                               .face_global
                               .value() <
                        value;
                });
        if (found ==
                faces_.end() ||
            found->face_global !=
                face_global) {
            throw std::out_of_range(
                "mpmc::flow_discretization_petsc: frozen coordinate registry does not contain requested authoritative face");
        }
        return side ==
                CrossCardinalityAbsentPhaseSide3D::
                    owner
            ? found->owner_cell_global
            : found->neighbour_cell_global;
    }

    [[nodiscard]]
    const FrozenAbsentPhaseCoordinateEntry3D&
    reference_entry(
        mpmc::mesh::GlobalEntityId
            cell_global,
        const mpmc::flow::
            FrozenPhysicalPhaseIdentity&
                identity) const {
        const auto& cell_entry =
            cell(cell_global);
        const auto found =
            std::lower_bound(
                cell_entry.absent_phases.begin(),
                cell_entry.absent_phases.end(),
                identity,
                [](const auto& entry,
                   const auto& value) {
                    return mpmc::flow::
                        frozen_phase_identity_less(
                            entry.identity,
                            value);
                });
        if (found ==
                cell_entry.absent_phases.end() ||
            found->identity !=
                identity) {
            throw std::out_of_range(
                "mpmc::flow_discretization_petsc: frozen coordinate registry lacks requested absent physical phase");
        }
        return *found;
    }

    [[nodiscard]]
    mpmc::flow::
        AbsentPhaseThermodynamicCoordinateExtension
    resolve_affine(
        mpmc::mesh::GlobalEntityId
            cell_global,
        const mpmc::flow::
            FrozenPhysicalPhaseIdentity&
                identity,
        const mpmc::flow::
            NaturalVariableStateIdentity3P&
                current_host_state) const {
        using namespace
            frozen_absent_phase_registry_detail;

        const auto& entry =
            reference_entry(
                cell_global,
                identity);
        const auto& reference =
            entry.reference_coordinates;
        if (!same_layout(
                reference
                    .host_state_identity
                    .layout,
                current_host_state.layout) ||
            reference
                    .host_state_identity
                    .component_ids !=
                current_host_state
                    .component_ids ||
            reference.identity !=
                identity) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: current host state is incompatible with frozen absent-phase coordinate chart");
        }

        const auto q_reference =
            natural_variables_from_state_identity(
                reference
                    .host_state_identity);
        const auto q_current =
            natural_variables_from_state_identity(
                current_host_state);
        if (q_reference.size() !=
            q_current.size()) {
            throw std::logic_error(
                "mpmc::flow_discretization_petsc: frozen/current coordinate widths differ");
        }

        mpmc::flow::
            AbsentPhaseThermodynamicCoordinateExtension
            result = reference;
        result.host_state_identity =
            current_host_state;

        for (std::size_t column = 0U;
             column < q_current.size();
             ++column) {
            const double delta =
                q_current[column] -
                q_reference[column];
            result.phase_pressure_pa +=
                result.phase_pressure_gradient[
                    column] *
                delta;
            for (std::size_t component = 0U;
                 component <
                    result.component_count();
                 ++component) {
                result
                    .hypothetical_composition[
                        component] +=
                    result.d_composition(
                        component,
                        column) *
                    delta;
            }
        }

        result.provenance =
            reference.provenance +
            "|frozen-branch=" +
            entry.selected_branch_provenance;
        result.validate();
        return result;
    }

private:
    void validate() {
        std::sort(
            cells_.begin(),
            cells_.end(),
            [](const auto& first,
               const auto& second) {
                return first.cell_global <
                    second.cell_global;
            });
        std::sort(
            faces_.begin(),
            faces_.end(),
            [](const auto& first,
               const auto& second) {
                return first.face_global <
                    second.face_global;
            });

        for (std::size_t index = 0U;
             index < cells_.size();
             ++index) {
            auto& cell_entry =
                cells_[index];
            if (index > 0U &&
                cells_[index - 1U]
                        .cell_global ==
                    cell_entry.cell_global) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: duplicate cell in frozen absent-phase coordinate registry");
            }
            std::sort(
                cell_entry.absent_phases.begin(),
                cell_entry.absent_phases.end(),
                [](const auto& first,
                   const auto& second) {
                    return mpmc::flow::
                        frozen_phase_identity_less(
                            first.identity,
                            second.identity);
                });

            for (std::size_t absent = 0U;
                 absent <
                    cell_entry
                        .absent_phases
                        .size();
                 ++absent) {
                const auto& entry =
                    cell_entry
                        .absent_phases[absent];
                entry.reference_coordinates
                    .validate();
                if (entry.identity !=
                        entry.reference_coordinates
                            .identity ||
                    entry.selected_branch_provenance
                        .empty() ||
                    cell_entry
                        .active_phases
                        .find(entry.identity)
                        .has_value() ||
                    entry.reference_coordinates
                            .host_state_identity
                            .layout
                            .phase_count() !=
                        cell_entry
                            .active_phases
                            .phase_count() ||
                    (absent > 0U &&
                     cell_entry
                             .absent_phases[
                                 absent - 1U]
                             .identity ==
                         entry.identity)) {
                    throw std::invalid_argument(
                        "mpmc::flow_discretization_petsc: malformed frozen absent-phase coordinate entry");
                }
            }
        }

        for (std::size_t index = 0U;
             index < faces_.size();
             ++index) {
            const auto& face =
                faces_[index];
            if (index > 0U &&
                faces_[index - 1U]
                        .face_global ==
                    face.face_global) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: duplicate face in frozen absent-phase coordinate registry");
            }
            (void)cell(
                face.owner_cell_global);
            (void)cell(
                face.neighbour_cell_global);
        }
    }

    std::vector<
        FrozenAbsentPhaseCoordinateCell3D>
        cells_;
    std::vector<
        FrozenAbsentPhaseFaceEndpoint3D>
        faces_;
};

[[nodiscard]] inline PetscErrorCode
resolve_frozen_absent_phase_thermodynamic_coordinates_3d(
    const MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D&
        face_input,
    const mpmc::flow::
        CrossCardinalityFacePhaseIdentityBinding&
            phase_binding,
    CrossCardinalityAbsentPhaseSide3D
        absent_side,
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        owner,
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        neighbour,
    void* raw_context,
    std::optional<
        mpmc::flow::
            AbsentPhaseThermodynamicCoordinateExtension>*
                output,
    NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        NaturalVariableSnesEvaluationStatus3D::
            success;

    try {
        auto* registry =
            static_cast<
                FrozenAbsentPhaseCoordinateRegistry3D*>(
                    raw_context);
        const auto cell_global =
            registry->face_cell(
                face_input.face_global,
                absent_side);
        const auto& current_cell =
            absent_side ==
                    CrossCardinalityAbsentPhaseSide3D::
                        owner
                ? owner
                : neighbour;
        const auto current_identity =
            std::visit(
                [](const auto& typed) {
                    return typed.transport
                        .state_identity;
                },
                current_cell);
        output->emplace(
            registry->resolve_affine(
                cell_global,
                phase_binding.identity,
                current_identity));
        return PETSC_SUCCESS;
    } catch (const std::domain_error&) {
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (const std::range_error&) {
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        return PETSC_ERR_ARG_INCOMP;
    } catch (const std::out_of_range&) {
        return PETSC_ERR_ARG_INCOMP;
    } catch (...) {
        return PETSC_ERR_LIB;
    }
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_FROZEN_ABSENT_PHASE_COORDINATE_REGISTRY_HPP
