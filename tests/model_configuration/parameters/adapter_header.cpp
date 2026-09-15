#include <mpmc/model_configuration/pr76_parameters.hpp>

bool adapter_header_compiles();
bool adapter_header_compiles() {
    return mpmc::model_configuration::ModelConfigurationLimits{}.max_components == 256;
}
