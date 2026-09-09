#include <mpmc/thermodynamics/pr76_mixture.hpp>
#include <mpmc/thermodynamics/pr76_mixture.hpp>

// The first/only public include must suffice for plain evaluation without AD.
double mixture_plain_header(const mpmc::thermodynamics::Pr76Mixture<double>& mixture) {
    mpmc::thermodynamics::Pr76MixtureWorkspace<double> workspace;
    return mixture.evaluate_reduced(300.0, {}, workspace).a;
}
