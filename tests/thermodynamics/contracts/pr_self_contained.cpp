#include <mpmc/thermodynamics/pr_parameters.hpp>
#include <mpmc/thermodynamics/pr_parameters.hpp>

#include <string_view>

std::string_view pr_profile_from_separate_translation_unit() {
    return mpmc::thermodynamics::PrParameterSet::model_id();
}
