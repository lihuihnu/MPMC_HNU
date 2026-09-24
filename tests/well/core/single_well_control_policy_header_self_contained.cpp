#include <mpmc/well/single_well_control_policy.hpp>

#include <string_view>

static_assert(
    mpmc::well::
        single_well_control_policy_convention ==
    std::string_view{
        "well/single-well-control-policy/fixed-total-molar-rate-minimum-bhp/v1"});

int main() {
    const mpmc::well::
        AcceptedSingleWellControlState
            state{
                mpmc::well::
                    SingleWellControlMode::
                        fixed_total_molar_rate,
                1.0};
    return state.valid()
        ? 0
        : 1;
}
