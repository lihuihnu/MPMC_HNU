#include <mpmc/thermodynamics/pr76_pure.hpp>
#include <mpmc/thermodynamics/pr76_pure.hpp>

// No AD header is included here. The ordinary scalar path remains usable alone.
static_assert(mpmc::thermodynamics::Pr76Pure<double>::gas_constant() > 8.0);
long double pr76_header_gas_constant() {
    return mpmc::thermodynamics::Pr76Pure<long double>::gas_constant();
}

// Instantiate a plain numerical evaluation in this separate translation unit.
// Its caller, like any public user, must supply a valid prepared object.
double pr76_plain_evaluation(const mpmc::thermodynamics::Pr76Pure<double>& pure, double t) {
    return pure.evaluate(t).a;
}
