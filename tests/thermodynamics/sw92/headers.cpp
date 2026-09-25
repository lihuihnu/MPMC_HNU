#include <mpmc/thermodynamics/sw92_types.hpp>
#include <mpmc/thermodynamics/sw92_parameters.hpp>
#include <mpmc/thermodynamics/sw92_correlations.hpp>
#include <mpmc/thermodynamics/sw92_pure.hpp>
#include <mpmc/thermodynamics/sw92_mixture.hpp>
#include <mpmc/thermodynamics/sw92_phase.hpp>
#include <mpmc/thermodynamics/selected_phase_fugacity.hpp>

bool sw92_headers() {
    return mpmc::thermodynamics::Sw92ParameterSet::model_id() ==
           mpmc::thermodynamics::sw92_corrected_profile;
}
