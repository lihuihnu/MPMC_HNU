#include <mpmc/model_configuration/model_definition.hpp>

// The public DTO must compile without importing scientific or transport headers.
#if defined(MPMC_THERMODYNAMICS_COMPONENTS_HPP) || defined(MPMC_FLASH_PT_FLASH_BACKEND_HPP)
#error "Public model definition leaked a backend dependency"
#endif

bool public_header_compiles();
bool public_header_compiles() {
    mpmc::model_configuration::ThermodynamicModelDefinition value;
    return value.family == mpmc::model_configuration::ThermodynamicModelFamily::unspecified &&
           value.version.empty() && std::holds_alternative<std::monostate>(value.parameters);
}
