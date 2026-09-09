#include <mpmc/thermodynamics/pr76_phase.hpp>
#include <mpmc/thermodynamics/pr76_phase.hpp>
#include <array>

double pt_plain_header(const mpmc::thermodynamics::Pr76Phase<double>& model) {
    mpmc::thermodynamics::Pr76PhaseWorkspace<double> workspace;
    const std::array<double,3> w{.25,.5,.25};
    return model.evaluate_full(1e6,450.0,w,2,workspace).ln_phi[0];
}
