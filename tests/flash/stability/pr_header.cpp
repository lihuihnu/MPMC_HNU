#include <mpmc/flash/pr76_stability.hpp>
#include <mpmc/flash/pr76_stability.hpp>

bool stability_pr_header(const mpmc::thermodynamics::Pr76Phase<double>& model) {
    mpmc::flash::Pr76StabilityEvaluator evaluator(model);
    return evaluator(1e5,700,std::vector<double>{1}).ln_phi.size()==1;
}
