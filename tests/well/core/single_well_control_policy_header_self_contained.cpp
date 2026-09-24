#include <mpmc/well/single_well_control_policy.hpp>

#include <string_view>

static_assert(
    mpmc::well::
        single_well_control_policy_convention ==
    std::string_view{
        "well/single-well-control-policy/fixed-total-molar-rate-minimum-bhp/v1"});

static_assert(
    mpmc::well::
        AcceptedSingleWellControlState{
            mpmc::well::SingleWellControlMode::
                fixed_total_molar_rate,
            1.0}
        .valid());

int main() {
    return 0;
}
