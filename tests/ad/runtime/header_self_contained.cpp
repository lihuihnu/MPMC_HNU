#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>

mpmc::ad::RuntimeValueAndJacobian<double> runtime_from_separate_translation_unit() {
    using Number = mpmc::ad::Dual<double, 2>;
    mpmc::ad::RuntimeJacobianWorkspace<double, 2> workspace;
    return mpmc::ad::value_and_jacobian_runtime<2>(
        [](std::span<const Number> input, std::span<Number> output) {
            for (std::size_t i = 0; i < input.size(); ++i) {
                output[i] = input[i];
            }
        }, std::vector<double>(5, 2.0), 5, workspace);
}
