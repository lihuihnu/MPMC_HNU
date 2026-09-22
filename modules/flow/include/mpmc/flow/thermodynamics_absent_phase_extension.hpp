#ifndef MPMC_FLOW_THERMODYNAMICS_ABSENT_PHASE_EXTENSION_HPP
#define MPMC_FLOW_THERMODYNAMICS_ABSENT_PHASE_EXTENSION_HPP

#include <mpmc/ad/dual.hpp>
#include <mpmc/flow/cross_cardinality_phase_identity.hpp>
#include <mpmc/flow/phase_identity_continuation.hpp>
#include <mpmc/thermodynamics/selected_phase_density.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace mpmc::flow {

inline constexpr std::string_view
    thermodynamics_absent_phase_extension_convention =
        "flow/thermodynamics/absent-phase-potential-extension/v1";

struct AbsentPhaseThermodynamicCoordinateExtension {
    FrozenPhysicalPhaseIdentity identity;
    NaturalVariableStateIdentity3P host_state_identity;

    double phase_pressure_pa{};
    std::vector<double> phase_pressure_gradient;

    std::vector<double> hypothetical_composition;
    std::vector<double> hypothetical_composition_jacobian;

    // Optional opaque selected-branch key. Generic coordinate validation does
    // not interpret it; cell-scoped providers may require it.
    std::string selected_branch_provenance;
    std::string provenance;

    [[nodiscard]] std::size_t input_count() const noexcept {
        return host_state_identity.layout.unknown_count();
    }

    [[nodiscard]] std::size_t component_count() const noexcept {
        return host_state_identity.layout.component_count();
    }

    [[nodiscard]] double d_composition(
        std::size_t component,
        std::size_t column) const {
        const std::size_t q = input_count();
        if (component >= component_count() || column >= q) {
            throw std::out_of_range(
                "mpmc::flow: absent-phase hypothetical composition Jacobian index out of range");
        }
        return hypothetical_composition_jacobian.at(
            component * q + column);
    }

    void validate() const {
        const std::size_t q = input_count();
        const std::size_t n = component_count();
        if (identity.provenance_scope.empty() ||
            identity.opaque_phase_key.empty() ||
            q == 0U || n == 0U ||
            host_state_identity.component_ids.size() != n ||
            !std::isfinite(phase_pressure_pa) ||
            !(phase_pressure_pa > 0.0) ||
            phase_pressure_gradient.size() != q ||
            hypothetical_composition.size() != n ||
            hypothetical_composition_jacobian.size() != n * q ||
            provenance.empty()) {
            throw std::invalid_argument(
                "mpmc::flow: malformed absent-phase thermodynamic coordinate extension");
        }

        double composition_sum = 0.0;
        for (double value : hypothetical_composition) {
            if (!std::isfinite(value) || !(value > 0.0)) {
                throw std::invalid_argument(
                    "mpmc::flow: hypothetical absent-phase composition must remain on strict positive support");
            }
            composition_sum += value;
        }
        const double sum_tolerance =
            4096.0 * std::numeric_limits<double>::epsilon() *
            static_cast<double>(n);
        if (!std::isfinite(composition_sum) ||
            std::abs(composition_sum - 1.0) > sum_tolerance) {
            throw std::invalid_argument(
                "mpmc::flow: hypothetical absent-phase composition is not normalized");
        }

        for (double value : phase_pressure_gradient) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "mpmc::flow: absent-phase pressure gradient contains non-finite derivative");
            }
        }

        for (std::size_t column = 0U; column < q; ++column) {
            double tangent_sum = 0.0;
            double tangent_scale = 0.0;
            for (std::size_t component = 0U; component < n; ++component) {
                const double value =
                    hypothetical_composition_jacobian[
                        component * q + column];
                if (!std::isfinite(value)) {
                    throw std::invalid_argument(
                        "mpmc::flow: absent-phase composition Jacobian contains non-finite derivative");
                }
                tangent_sum += value;
                tangent_scale += std::abs(value);
            }
            const double tolerance =
                4096.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, tangent_scale);
            if (!std::isfinite(tangent_sum) ||
                std::abs(tangent_sum) > tolerance) {
                throw std::invalid_argument(
                    "mpmc::flow: absent-phase composition derivative must remain tangent to the simplex");
            }
        }
    }
};

template <typename Selection>
struct FrozenSelectedPhaseBranchBinding {
    FrozenPhysicalPhaseIdentity identity;
    Selection selection;
    std::string provenance;
};

template <typename Selection>
class FrozenSelectedPhaseBranchRegistry {
public:
    FrozenSelectedPhaseBranchRegistry(
        std::vector<FrozenSelectedPhaseBranchBinding<Selection>> bindings,
        std::string transition_evidence_profile)
        : bindings_(std::move(bindings)),
          transition_evidence_profile_(
              std::move(transition_evidence_profile)) {
        validate();
    }

    [[nodiscard]] const FrozenSelectedPhaseBranchBinding<Selection>&
    binding(const FrozenPhysicalPhaseIdentity& identity) const {
        const auto found =
            std::find_if(
                bindings_.begin(),
                bindings_.end(),
                [&](const auto& entry) {
                    return entry.identity == identity;
                });
        if (found == bindings_.end()) {
            throw std::out_of_range(
                "mpmc::flow: physical phase identity has no frozen thermodynamic branch binding");
        }
        return *found;
    }

    [[nodiscard]] std::span<
        const FrozenSelectedPhaseBranchBinding<Selection>>
    bindings() const noexcept {
        return bindings_;
    }

    [[nodiscard]] std::string_view
    transition_evidence_profile() const noexcept {
        return transition_evidence_profile_;
    }

private:
    void validate() const {
        if (bindings_.empty() ||
            transition_evidence_profile_.empty()) {
            throw std::invalid_argument(
                "mpmc::flow: frozen selected-phase branch registry requires bindings and transition provenance");
        }
        for (std::size_t index = 0U; index < bindings_.size(); ++index) {
            const auto& entry = bindings_[index];
            if (entry.identity.provenance_scope.empty() ||
                entry.identity.opaque_phase_key.empty() ||
                entry.provenance.empty()) {
                throw std::invalid_argument(
                    "mpmc::flow: selected-phase branch binding lacks explicit identity/provenance");
            }
            for (std::size_t previous = 0U;
                 previous < index;
                 ++previous) {
                if (bindings_[previous].identity ==
                    entry.identity) {
                    throw std::invalid_argument(
                        "mpmc::flow: duplicate physical phase identity in selected-phase branch registry");
                }
            }
        }
    }

    std::vector<FrozenSelectedPhaseBranchBinding<Selection>>
        bindings_;
    std::string transition_evidence_profile_;
};

template <typename Selection>
[[nodiscard]] inline FrozenSelectedPhaseBranchRegistry<Selection>
make_transition_selected_phase_branch_registry(
    const PhaseIdentityContinuationSnapshot& continuation,
    std::vector<FrozenSelectedPhaseBranchBinding<Selection>> bindings) {
    const auto universe = continuation.phase_universe();
    if (bindings.size() != universe.size()) {
        throw std::invalid_argument(
            "mpmc::flow: transition selected-phase branch registry must cover the complete source/target physical phase universe");
    }
    for (const auto& identity : universe) {
        const bool present =
            std::any_of(
                bindings.begin(),
                bindings.end(),
                [&](const auto& entry) {
                    return entry.identity == identity;
                });
        if (!present) {
            throw std::invalid_argument(
                "mpmc::flow: transition selected-phase branch registry is missing an explicit physical phase binding");
        }
    }
    for (const auto& entry : bindings) {
        const bool belongs =
            std::find(
                universe.begin(),
                universe.end(),
                entry.identity) !=
            universe.end();
        if (!belongs) {
            throw std::invalid_argument(
                "mpmc::flow: selected thermodynamic branch identity is not part of the transition phase universe");
        }
    }
    return {
        std::move(bindings),
        std::string{continuation.evidence_profile()}};
}

namespace absent_phase_thermodynamics_detail {

template <typename Parameters>
[[nodiscard]] inline std::vector<double>
molar_masses_kg_per_mol(
    const Parameters& parameters,
    std::span<const std::string> component_ids) {
    const auto& components = parameters.components();
    if (components.size() != component_ids.size()) {
        throw std::invalid_argument(
            "mpmc::flow: thermodynamic model component count differs from absent-phase extension state");
    }
    std::vector<double> result;
    result.reserve(component_ids.size());
    for (std::size_t index = 0U;
         index < component_ids.size();
         ++index) {
        const auto& component =
            components.at(index);
        if (component.id != component_ids[index] ||
            !component.molar_mass.has_value() ||
            component.molar_mass->unit !=
                thermodynamics::Unit::kilogram_per_mole ||
            !std::isfinite(component.molar_mass->value) ||
            !(component.molar_mass->value > 0.0)) {
            throw std::invalid_argument(
                "mpmc::flow: selected-phase mass density requires explicit ordered component molar masses");
        }
        result.push_back(
            component.molar_mass->value);
    }
    return result;
}

template <std::size_t K, typename DensityEvaluator>
[[nodiscard]] inline AbsentPhasePotentialExtensionLinearization
evaluate_blocked_mass_density(
    const AbsentPhaseThermodynamicCoordinateExtension& coordinates,
    std::span<const double> molar_masses,
    std::string provider_provenance,
    DensityEvaluator&& evaluate_molar_density) {
    using D = mpmc::ad::Dual<double, K>;
    coordinates.validate();

    const std::size_t q = coordinates.input_count();
    const std::size_t n = coordinates.component_count();
    if (molar_masses.size() != n) {
        throw std::invalid_argument(
            "mpmc::flow: molar-mass vector shape mismatch in absent-phase density provider");
    }

    std::vector<double> mass_density_gradient(q, 0.0);
    double mass_density_value =
        std::numeric_limits<double>::quiet_NaN();

    const std::size_t block_count =
        (q + K - 1U) / K;
    for (std::size_t block = 0U;
         block < block_count;
         ++block) {
        typename D::Gradient pressure_seed{};
        typename D::Gradient temperature_seed{};
        std::vector<typename D::Gradient>
            composition_seed(n);

        for (std::size_t lane = 0U;
             lane < K;
             ++lane) {
            const std::size_t column =
                block * K + lane;
            if (column >= q) {
                continue;
            }
            pressure_seed[lane] =
                coordinates
                    .phase_pressure_gradient[column];
            if (column ==
                coordinates.host_state_identity
                    .layout.temperature_unknown_index()) {
                temperature_seed[lane] = 1.0;
            }
            for (std::size_t component = 0U;
                 component < n;
                 ++component) {
                composition_seed[component][lane] =
                    coordinates.d_composition(
                        component,
                        column);
            }
        }

        const D pressure{
            coordinates.phase_pressure_pa,
            pressure_seed};
        const D temperature{
            coordinates.host_state_identity.temperature_k,
            temperature_seed};
        std::vector<D> composition;
        composition.reserve(n);
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            composition.emplace_back(
                coordinates.hypothetical_composition[component],
                composition_seed[component]);
        }

        const D molar_density =
            evaluate_molar_density(
                pressure,
                temperature,
                std::span<const D>{composition});

        D mixture_molar_mass{0.0};
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            mixture_molar_mass +=
                composition[component] *
                molar_masses[component];
        }
        const D mass_density =
            molar_density *
            mixture_molar_mass;

        if (!std::isfinite(mass_density.value()) ||
            !(mass_density.value() > 0.0)) {
            throw std::range_error(
                "mpmc::flow: selected thermodynamic branch produced invalid absent-phase mass density");
        }
        if (block == 0U) {
            mass_density_value =
                mass_density.value();
        } else {
            const double tolerance =
                4096.0 *
                std::numeric_limits<double>::epsilon() *
                std::max(
                    {1.0,
                     std::abs(mass_density_value),
                     std::abs(mass_density.value())});
            if (std::abs(
                    mass_density.value() -
                    mass_density_value) >
                tolerance) {
                throw std::runtime_error(
                    "mpmc::flow: selected-phase density primal changed across AD blocks");
            }
        }

        for (std::size_t lane = 0U;
             lane < K;
             ++lane) {
            const std::size_t column =
                block * K + lane;
            if (column >= q) {
                continue;
            }
            const double derivative =
                mass_density.derivative(lane);
            if (!std::isfinite(derivative)) {
                throw std::range_error(
                    "mpmc::flow: selected thermodynamic branch produced non-finite absent-phase density derivative");
            }
            mass_density_gradient[column] =
                derivative;
        }
    }

    if (provider_provenance.empty()) {
        throw std::invalid_argument(
            "mpmc::flow: absent-phase thermodynamic provider provenance must not be empty");
    }
    return {
        coordinates.identity,
        coordinates.phase_pressure_pa,
        mass_density_value,
        q,
        coordinates.phase_pressure_gradient,
        std::move(mass_density_gradient),
        std::move(provider_provenance)};
}

inline std::string provider_provenance(
    std::string_view backend,
    std::string_view transition,
    std::string_view branch,
    std::string_view coordinates) {
    if (backend.empty() ||
        transition.empty() ||
        branch.empty() ||
        coordinates.empty()) {
        throw std::invalid_argument(
            "mpmc::flow: incomplete absent-phase provider provenance");
    }
    return std::string{backend} +
        "|transition=" + std::string{transition} +
        "|branch=" + std::string{branch} +
        "|coordinates=" + std::string{coordinates};
}

} // namespace absent_phase_thermodynamics_detail

template <std::floating_point T = double>
class Pr76AbsentPhasePotentialExtensionProvider {
public:
    using Selection = thermodynamics::Pr76SelectedPhase;

    Pr76AbsentPhasePotentialExtensionProvider(
        const thermodynamics::Pr76Phase<T>& model,
        FrozenSelectedPhaseBranchRegistry<Selection> registry)
        : model_(&model),
          registry_(std::move(registry)) {}

    [[nodiscard]] AbsentPhasePotentialExtensionLinearization
    evaluate(
        const AbsentPhaseThermodynamicCoordinateExtension&
            coordinates) const {
        const auto& binding =
            registry_.binding(coordinates.identity);
        const auto masses =
            absent_phase_thermodynamics_detail::
                molar_masses_kg_per_mol(
                    model_->parameters(),
                    coordinates
                        .host_state_identity
                        .component_ids);

        return absent_phase_thermodynamics_detail::
            evaluate_blocked_mass_density<4U>(
                coordinates,
                masses,
                absent_phase_thermodynamics_detail::
                    provider_provenance(
                        thermodynamics::pr76_pt_convention,
                        registry_
                            .transition_evidence_profile(),
                        binding.provenance,
                        coordinates.provenance),
                [&](const auto& pressure,
                    const auto& temperature,
                    auto composition) {
                    using Number =
                        std::remove_cvref_t<
                            decltype(pressure)>;
                    thermodynamics::
                        Pr76PhaseWorkspace<Number>
                        workspace;
                    return thermodynamics::
                        evaluate_selected_phase_molar_density(
                            *model_,
                            pressure,
                            temperature,
                            composition,
                            binding.selection,
                            workspace)
                            .molar_density_mol_per_m3;
                });
    }

private:
    const thermodynamics::Pr76Phase<T>* model_;
    FrozenSelectedPhaseBranchRegistry<Selection>
        registry_;
};

template <std::floating_point T = double>
class Pr76CellScopedAbsentPhasePotentialExtensionProvider {
public:
    using Selection =
        thermodynamics::Pr76SelectedPhase;

    struct Binding {
        FrozenPhysicalPhaseIdentity identity;
        Selection selection;
        std::string transition_evidence_profile;
        std::string branch_provenance;
    };

    Pr76CellScopedAbsentPhasePotentialExtensionProvider(
        const thermodynamics::Pr76Phase<T>& model,
        std::vector<Binding> bindings)
        : model_(&model),
          bindings_(std::move(bindings)) {
        std::sort(
            bindings_.begin(),
            bindings_.end(),
            [](const auto& first,
               const auto& second) {
                return first.branch_provenance <
                    second.branch_provenance;
            });
        for (std::size_t index = 0U;
             index < bindings_.size();
             ++index) {
            const auto& binding =
                bindings_[index];
            if (binding.identity.provenance_scope.empty() ||
                binding.identity.opaque_phase_key.empty() ||
                binding.transition_evidence_profile.empty() ||
                binding.branch_provenance.empty() ||
                (index > 0U &&
                 bindings_[index - 1U]
                         .branch_provenance ==
                     binding.branch_provenance)) {
                throw std::invalid_argument(
                    "mpmc::flow: invalid cell-scoped PR76 absent-phase branch binding");
            }
        }
    }

    [[nodiscard]] std::size_t
    binding_count() const noexcept {
        return bindings_.size();
    }

    [[nodiscard]] std::size_t
    evaluation_count() const noexcept {
        return evaluation_count_;
    }

    [[nodiscard]]
    AbsentPhasePotentialExtensionLinearization
    evaluate(
        const AbsentPhaseThermodynamicCoordinateExtension&
            coordinates) const {
        if (coordinates
                .selected_branch_provenance
                .empty()) {
            throw std::invalid_argument(
                "mpmc::flow: cell-scoped PR76 absent-phase coordinates lack selected-branch provenance");
        }

        const auto found =
            std::lower_bound(
                bindings_.begin(),
                bindings_.end(),
                coordinates
                    .selected_branch_provenance,
                [](const auto& binding,
                   const auto& provenance) {
                    return binding
                               .branch_provenance <
                        provenance;
                });
        if (found == bindings_.end() ||
            found->branch_provenance !=
                coordinates
                    .selected_branch_provenance ||
            found->identity !=
                coordinates.identity) {
            throw std::out_of_range(
                "mpmc::flow: cell-scoped PR76 absent-phase branch binding not found");
        }

        ++evaluation_count_;
        const auto masses =
            absent_phase_thermodynamics_detail::
                molar_masses_kg_per_mol(
                    model_->parameters(),
                    coordinates
                        .host_state_identity
                        .component_ids);

        return absent_phase_thermodynamics_detail::
            evaluate_blocked_mass_density<4U>(
                coordinates,
                masses,
                absent_phase_thermodynamics_detail::
                    provider_provenance(
                        thermodynamics::
                            pr76_pt_convention,
                        found
                            ->transition_evidence_profile,
                        found->branch_provenance,
                        coordinates.provenance),
                [&](const auto& pressure,
                    const auto& temperature,
                    auto composition) {
                    using Number =
                        std::remove_cvref_t<
                            decltype(pressure)>;
                    thermodynamics::
                        Pr76PhaseWorkspace<Number>
                        workspace;
                    return thermodynamics::
                        evaluate_selected_phase_molar_density(
                            *model_,
                            pressure,
                            temperature,
                            composition,
                            found->selection,
                            workspace)
                            .molar_density_mol_per_m3;
                });
    }

private:
    const thermodynamics::Pr76Phase<T>* model_;
    std::vector<Binding> bindings_;
    mutable std::size_t evaluation_count_{};
};

template <std::floating_point T = double>
class Sw92AbsentPhasePotentialExtensionProvider {
public:
    using Selection =
        thermodynamics::Sw92SelectedPhase<T>;

    Sw92AbsentPhasePotentialExtensionProvider(
        const thermodynamics::Sw92Phase<T>& model,
        FrozenSelectedPhaseBranchRegistry<Selection> registry)
        : model_(&model),
          registry_(std::move(registry)) {}

    [[nodiscard]] AbsentPhasePotentialExtensionLinearization
    evaluate(
        const AbsentPhaseThermodynamicCoordinateExtension&
            coordinates) const {
        const auto& binding =
            registry_.binding(coordinates.identity);
        const auto masses =
            absent_phase_thermodynamics_detail::
                molar_masses_kg_per_mol(
                    model_->parameters(),
                    coordinates
                        .host_state_identity
                        .component_ids);

        return absent_phase_thermodynamics_detail::
            evaluate_blocked_mass_density<4U>(
                coordinates,
                masses,
                absent_phase_thermodynamics_detail::
                    provider_provenance(
                        thermodynamics::sw92_pt_convention,
                        registry_
                            .transition_evidence_profile(),
                        binding.provenance,
                        coordinates.provenance),
                [&](const auto& pressure,
                    const auto& temperature,
                    auto composition) {
                    using Number =
                        std::remove_cvref_t<
                            decltype(pressure)>;
                    thermodynamics::
                        Sw92PhaseWorkspace<Number>
                        workspace;
                    return thermodynamics::
                        evaluate_selected_phase_molar_density(
                            *model_,
                            pressure,
                            temperature,
                            composition,
                            binding.selection,
                            workspace)
                            .molar_density_mol_per_m3;
                });
    }

private:
    const thermodynamics::Sw92Phase<T>* model_;
    FrozenSelectedPhaseBranchRegistry<Selection>
        registry_;
};

class CpaAbsentPhasePotentialExtensionProvider {
public:
    using Selection =
        thermodynamics::CpaSelectedPhase;

    CpaAbsentPhasePotentialExtensionProvider(
        const thermodynamics::CpaPtPhase& model,
        FrozenSelectedPhaseBranchRegistry<Selection> registry)
        : model_(&model),
          registry_(std::move(registry)) {}

    [[nodiscard]] AbsentPhasePotentialExtensionLinearization
    evaluate(
        const AbsentPhaseThermodynamicCoordinateExtension&
            coordinates) const {
        const auto& binding =
            registry_.binding(coordinates.identity);
        const auto masses =
            absent_phase_thermodynamics_detail::
                molar_masses_kg_per_mol(
                    model_->parameters(),
                    coordinates
                        .host_state_identity
                        .component_ids);

        return absent_phase_thermodynamics_detail::
            evaluate_blocked_mass_density<4U>(
                coordinates,
                masses,
                absent_phase_thermodynamics_detail::
                    provider_provenance(
                        thermodynamics::cpa_pt_convention,
                        registry_
                            .transition_evidence_profile(),
                        binding.provenance,
                        coordinates.provenance),
                [&](const auto& pressure,
                    const auto& temperature,
                    auto composition) {
                    return thermodynamics::
                        evaluate_selected_phase_molar_density(
                            *model_,
                            pressure,
                            temperature,
                            composition,
                            binding.selection)
                            .molar_density_mol_per_m3;
                });
    }

private:
    const thermodynamics::CpaPtPhase* model_;
    FrozenSelectedPhaseBranchRegistry<Selection>
        registry_;
};

} // namespace mpmc::flow

#endif // MPMC_FLOW_THERMODYNAMICS_ABSENT_PHASE_EXTENSION_HPP
