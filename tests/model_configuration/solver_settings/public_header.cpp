#include <mpmc/model_configuration/pt_solver_settings.hpp>

#if defined(MPMC_THERMODYNAMICS_COMPONENTS_HPP) || defined(MPMC_FLASH_PT_STABILITY_HPP) || \
    defined(MPMC_FLASH_PR76_PT_FLASH_BACKEND_HPP)
#error "Public solver settings must not import scientific implementations"
#endif

bool solver_public_header_compiles();
bool solver_public_header_compiles() {
    const mpmc::model_configuration::PtSolverSettings empty;
    return !empty.initial_stability.automatic_multistart.has_value() &&
           !empty.two_phase.max_property_evaluations.has_value();
}
