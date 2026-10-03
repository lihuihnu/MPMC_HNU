#include <mpmc/ad/dual.hpp>
#include <mpmc/flow/natural_variable_cell_state.hpp>
#include <mpmc/thermodynamics/pr76_pure.hpp>

// Exercise usage requirements of mpmc::flow_thermodynamics in isolation.
// These synthetic values are not a physical flow or EOS validation case.
int main() {
    const mpmc::flow::NaturalVariableLayoutDescriptor layout{2, 1, {1}};
    using Number = mpmc::ad::Dual<double, 1>;
    const mpmc::thermodynamics::Pr76PureValues<Number> coefficients{
        Number::variable(2.0, 0), Number{1.0}};
    return layout.unknown_count() == 3 && coefficients.a.derivative(0) == 1.0
               ? 0
               : 1;
}
