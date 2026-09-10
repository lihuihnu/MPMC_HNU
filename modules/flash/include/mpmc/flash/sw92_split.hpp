#ifndef MPMC_FLASH_SW92_SPLIT_HPP
#define MPMC_FLASH_SW92_SPLIT_HPP

#include <mpmc/flash/pt_split.hpp>
#include <mpmc/flash/sw92_stability.hpp>

#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flash {

/// Gate-1 equilibrium-algorithm identity. This labels the internally consistent
/// fixed-family primitive only; it is neither of the cross-family profiles.
inline constexpr std::string_view sw92_family_vle_algorithm =
    "SW92-equilibrium/fixed-family-vle-primitive/v1";

/// One-family SW92 VLE provider for the generic PT split driver.
///
/// One instance owns exactly one thermodynamic snapshot, NaCl molality and
/// phase family. The three-argument call is the all-admissible-root stability
/// provider. The four-argument call requests the lowest/highest-Z mechanically
/// admissible candidate root inside that SAME family. The requested role is a
/// numerical split role, not a proof of physical aqueous/non-aqueous identity.
class Sw92FamilyVleEvaluator {
public:
    Sw92FamilyVleEvaluator(
        const thermodynamics::Sw92Phase<double>& model,
        double nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily family,
        thermodynamics::Sw92RootOptions root_options = {})
        : stability_(model, nacl_molality_mol_per_kg_water, family, root_options) {}

    Sw92FamilyVleEvaluator(const Sw92FamilyVleEvaluator&) = delete;
    Sw92FamilyVleEvaluator& operator=(const Sw92FamilyVleEvaluator&) = delete;
    Sw92FamilyVleEvaluator(Sw92FamilyVleEvaluator&&) = delete;
    Sw92FamilyVleEvaluator& operator=(Sw92FamilyVleEvaluator&&) = delete;

    [[nodiscard]] const thermodynamics::Sw92Phase<double>& model() const & noexcept {
        return stability_.model();
    }
    const thermodynamics::Sw92Phase<double>& model() const && = delete;
    [[nodiscard]] double nacl_molality_mol_per_kg_water() const noexcept {
        return stability_.nacl_molality_mol_per_kg_water();
    }
    [[nodiscard]] thermodynamics::SwPhaseFamily family() const noexcept {
        return stability_.family();
    }
    [[nodiscard]] thermodynamics::Sw92RootOptions root_options() const noexcept {
        return stability_.root_options();
    }

    [[nodiscard]] StabilityPhase operator()(
        double pressure_pa, double temperature_k, std::span<const double> composition) {
        return stability_(pressure_pa, temperature_k, composition);
    }

    [[nodiscard]] PtSplitPhase operator()(
        double pressure_pa, double temperature_k, std::span<const double> composition,
        PtPhaseRole role) {
        namespace th = thermodynamics;
        if (role != PtPhaseRole::liquid_candidate && role != PtPhaseRole::vapor_candidate) {
            throw std::invalid_argument("SW92 VLE: unknown candidate role");
        }
        try {
            const auto roots = model().roots(
                pressure_pa, temperature_k, composition,
                nacl_molality_mol_per_kg_water(), family(), workspace_, root_options());
            switch (roots.status) {
            case th::Sw92RootStatus::near_multiple:
                throw StabilityPropertyError(StabilityPropertyIssue::root_topology,
                    "SW92 VLE: unresolved root topology");
            case th::Sw92RootStatus::iteration_limit:
                throw StabilityPropertyError(StabilityPropertyIssue::root_iteration_limit,
                    "SW92 VLE: root iteration limit");
            case th::Sw92RootStatus::unrepresentable:
                throw StabilityPropertyError(StabilityPropertyIssue::root_range,
                    "SW92 VLE: nonrepresentable root");
            case th::Sw92RootStatus::success:
                break;
            }

            std::optional<std::size_t> selected;
            for (std::size_t root_index = 0; root_index < roots.count; ++root_index) {
                if (roots.roots[root_index].slope_sign <= 0) { continue; }
                selected = root_index;
                if (role == PtPhaseRole::liquid_candidate) { break; }
            }
            if (!selected) {
                throw StabilityPropertyError(StabilityPropertyIssue::no_admissible_branch,
                    "SW92 VLE: no mechanically admissible root in selected family");
            }
            if (!roots.roots[*selected].derivative_valid) {
                throw StabilityPropertyError(StabilityPropertyIssue::ill_conditioned_root,
                    "SW92 VLE: ill-conditioned candidate root");
            }

            auto values = model().evaluate(
                pressure_pa, temperature_k, composition,
                nacl_molality_mol_per_kg_water(), family(), *selected,
                workspace_, root_options());
            return {{std::move(values.ln_phi), *selected, true}, values.z};
        } catch (const th::Sw92PhaseError& error) {
            StabilityPropertyIssue issue = StabilityPropertyIssue::root_range;
            switch (error.code()) {
            case th::Sw92PhaseErrorCode::near_multiple:
                issue = StabilityPropertyIssue::root_topology;
                break;
            case th::Sw92PhaseErrorCode::iteration_limit:
                issue = StabilityPropertyIssue::root_iteration_limit;
                break;
            case th::Sw92PhaseErrorCode::unrepresentable_root:
                break;
            }
            throw StabilityPropertyError(issue, error.what());
        } catch (const std::range_error& error) {
            throw StabilityPropertyError(StabilityPropertyIssue::root_range, error.what());
        }
    }

private:
    Sw92FamilyStabilityEvaluator stability_;
    thermodynamics::Sw92PhaseWorkspace<double> workspace_;
};

/// Gate-1 one-family result. This is an internally consistent equilibrium result
/// for one fixed SW92 family. It is NOT the dual-model observable profile and it
/// is NOT an AQ/NA asymmetric joint-equilibrium result.
struct Sw92FamilyPtSplitResult {
    PtSplitResult solution;
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::sw92_corrected_profile};
    std::string phase_convention{thermodynamics::sw92_pt_convention};
    std::string equilibrium_algorithm{sw92_family_vle_algorithm};
    double nacl_molality_mol_per_kg_water{};
    thermodynamics::SwPhaseFamily family{thermodynamics::SwPhaseFamily::aqueous};
    thermodynamics::Sw92RootOptions root_options;
};

/// Full one-family baseline:
/// fixed-family feed stability -> material-balanced split -> fixed-family final
/// common-tangent review. Both provider roles are served by the SAME family,
/// model snapshot, molality and component reference states.
[[nodiscard]] inline Sw92FamilyPtSplitResult solve_sw92_pt_family_vle(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    Sw92FamilyVleEvaluator& evaluator, PtSplitOptions options = {},
    std::span<const std::vector<double>> initial_starts = {},
    std::span<const std::vector<double>> final_starts = {}) {
    if (feed.size() != evaluator.model().size()) {
        throw std::invalid_argument("SW92 VLE: feed does not match ordered model snapshot");
    }
    Sw92FamilyPtSplitResult result;
    const auto& parameters = evaluator.model().parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    result.nacl_molality_mol_per_kg_water = evaluator.nacl_molality_mol_per_kg_water();
    result.family = evaluator.family();
    result.root_options = evaluator.root_options();
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }
    result.solution = solve_pt_vle(
        pressure_pa, temperature_k, feed, evaluator, evaluator,
        options, initial_starts, final_starts);
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_SPLIT_HPP
