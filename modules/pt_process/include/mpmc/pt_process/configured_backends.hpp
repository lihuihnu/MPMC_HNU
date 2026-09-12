#ifndef MPMC_PT_PROCESS_CONFIGURED_BACKENDS_HPP
#define MPMC_PT_PROCESS_CONFIGURED_BACKENDS_HPP

#include <mpmc/flash/cpa_pt_flash_backend.hpp>
#include <mpmc/flash/pr76_pt_flash_backend.hpp>
#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>
#include <mpmc/pt_process/composition_root.hpp>
#include <mpmc/thermodynamics/cpa_parameters.hpp>
#include <mpmc/thermodynamics/pr_parameters.hpp>
#include <mpmc/thermodynamics/sw92_parameters.hpp>

#include <string>

namespace mpmc::pt_process {

// These factories only assemble already-validated parameter snapshots with the
// established model/evaluator/backend constructors. They contain no model
// equations and never call solve(). The returned aliasing pointer retains the
// entire referenced object graph.
[[nodiscard]] OwnedConfiguredPtBackend make_pr76_configured_backend(
    std::string configured_backend_id,
    const thermodynamics::PrParameterSet& parameters,
    thermodynamics::Pr76RootOptions root_options = {},
    flash::Pr76PtFlashBackendOptions backend_options = {});

[[nodiscard]] OwnedConfiguredPtBackend make_sw92_configured_backend(
    std::string configured_backend_id,
    const thermodynamics::Sw92ParameterSet& parameters,
    flash::Sw92ProfileCPtFlashBackendOptions backend_options = {});

[[nodiscard]] OwnedConfiguredPtBackend make_cpa_configured_backend(
    std::string configured_backend_id,
    const thermodynamics::CpaParameterSet& parameters,
    thermodynamics::CpaPtOptions phase_options =
        flash::cpa_pt_vle_default_phase_options(),
    flash::CpaPtFlashBackendOptions backend_options = {});

} // namespace mpmc::pt_process

#endif // MPMC_PT_PROCESS_CONFIGURED_BACKENDS_HPP
