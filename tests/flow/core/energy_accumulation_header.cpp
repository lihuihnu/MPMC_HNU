#include <mpmc/flow/energy_accumulation.hpp>

#include <string_view>

static_assert(
    mpmc::flow::
        BackwardEulerEnergyAccumulationResidual3P::
            convention ==
    mpmc::flow::
        backward_euler_energy_accumulation_convention);

bool energy_accumulation_header() {
    return true;
}
