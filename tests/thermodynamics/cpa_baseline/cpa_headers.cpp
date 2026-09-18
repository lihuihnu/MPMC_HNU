#include <mpmc/thermodynamics/cpa_parameters.hpp>
#include <mpmc/thermodynamics/cpa_association.hpp>
#include <mpmc/thermodynamics/cpa_phase.hpp>

int main() {
    return mpmc::thermodynamics::cpa_profile.size() > 0U ? 0 : 1;
}
