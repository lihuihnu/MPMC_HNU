#ifndef MPMC_PT_PROCESS_PARAMETER_SNAPSHOT_SUPPLIER_HPP
#define MPMC_PT_PROCESS_PARAMETER_SNAPSHOT_SUPPLIER_HPP

#include <mpmc/pt_process/composition_root.hpp>
#include <mpmc/thermodynamics/cpa_parameters.hpp>
#include <mpmc/thermodynamics/pr_parameters.hpp>
#include <mpmc/thermodynamics/sw92_parameters.hpp>

#include <iosfwd>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::pt_process {

inline constexpr std::string_view pt_parameter_snapshot_bundle_convention =
    "MPMC/PT/parameter-snapshot-bundle/v1";
inline constexpr std::string_view repository_curated_pt_bundle_id =
    "MPMC/PT/repository-curated-literature-snapshots/v1";
inline constexpr std::string_view repository_curated_pt_bundle_revision = "r1";

// One cited source used by a compiled, immutable parameter snapshot. This is
// operational discovery metadata; complete per-datum provenance remains in
// each validated thermodynamic parameter set.
struct PtParameterSnapshotSource {
    std::string reference;
    std::string revision;
    std::string scope;

    [[nodiscard]] bool structurally_valid() const noexcept;
};

struct PtParameterSnapshotDescriptor {
    std::string configured_backend_id;
    std::string snapshot_location;
    std::string model_profile;
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string applicability_scope;
    std::vector<PtParameterSnapshotSource> sources;

    [[nodiscard]] bool structurally_valid() const noexcept;
};

// The supplier returns both owned backend graphs and their public audit
// manifest. Construction validates source records through the established
// thermodynamic contracts, then cross-checks the manifest against backend
// capabilities. It never executes a flash solve.
struct PtParameterSnapshotBundle {
    std::string convention;
    std::string bundle_id;
    std::string revision;
    std::string location;
    std::vector<PtParameterSnapshotDescriptor> snapshots;
    std::vector<OwnedConfiguredPtBackend> backends;

    [[nodiscard]] bool structurally_valid() const noexcept;
};

// Model-aware data supplier kept separate from the process host and adapter.
// Tests and non-service library consumers may audit the exact owning snapshots
// before any backend graph is constructed.
struct RepositoryCuratedPtParameterSnapshotsV1 {
    thermodynamics::PrParameterSet pr76;
    thermodynamics::Sw92ParameterSet sw92;
    thermodynamics::CpaParameterSet cpa;
};

[[nodiscard]] RepositoryCuratedPtParameterSnapshotsV1
load_repository_curated_pt_parameters_v1();

// Repository-curated, literature-sourced snapshots with deliberately narrow
// applicability. This is not a general component database and does not infer,
// fit, or fill any missing parameter.
[[nodiscard]] PtParameterSnapshotBundle
load_repository_curated_pt_parameter_snapshots_v1();

// Emits identifiers, component order, applicability, and citations only. No
// numerical model parameters, request data, or credentials are serialized.
void write_pt_parameter_snapshot_manifest_json(
    std::ostream& output, const PtParameterSnapshotBundle& bundle);

} // namespace mpmc::pt_process

#endif // MPMC_PT_PROCESS_PARAMETER_SNAPSHOT_SUPPLIER_HPP
