#include <mpmc/physics/sw92_thermodynamic_closure.hpp>

#include <string_view>

bool sw92_thermodynamic_closure_header() {
    return mpmc::physics::PtPhaseSetThermodynamicClosureSnapshot::convention ==
               std::string_view{"PT/phase-set/thermodynamic-closure-primal-v1"} &&
           mpmc::physics::sw92_profile_c_closure_convention ==
               std::string_view{
                   "SW92/phase-assigned/Profile-C/physics-thermodynamic-closure-primal/v1"};
}
