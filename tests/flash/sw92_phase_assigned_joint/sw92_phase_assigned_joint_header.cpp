#include <mpmc/flash/sw92_phase_assigned_joint.hpp>

bool sw92_phase_assigned_joint_header() {
    return !mpmc::flash::sw92_phase_assigned_aq_na_joint_profile.empty() &&
           !mpmc::flash::sw92_phase_assigned_aq_na_joint_primitive.empty();
}
