#include <mpmc/thermodynamics/cpa_pt_phase.hpp>
#include <mpmc/thermodynamics/cpa_pt_phase_ad.hpp>
#include <mpmc/thermodynamics/selected_phase_fugacity.hpp>

int main() {
    return mpmc::thermodynamics::cpa_pt_convention.size() > 0U ? 0 : 1;
}
