#include <mpmc/pt_process/configured_backends.hpp>

#include <memory>
#include <utility>

namespace mpmc::pt_process {
namespace {

struct Pr76BackendGraph {
    thermodynamics::Pr76Phase<double> model;
    flash::Pr76VleEvaluator evaluator;
    flash::Pr76PtFlashBackend backend;

    Pr76BackendGraph(
        const thermodynamics::PrParameterSet& parameters,
        thermodynamics::Pr76RootOptions root_options,
        flash::Pr76PtFlashBackendOptions backend_options)
        : model(thermodynamics::Pr76Phase<double>::from_parameters(parameters)),
          evaluator(model, root_options),
          backend(evaluator, std::move(backend_options)) {}
};

struct Sw92BackendGraph {
    thermodynamics::Sw92Phase<double> model;
    flash::Sw92ProfileCPtFlashBackend backend;

    Sw92BackendGraph(
        const thermodynamics::Sw92ParameterSet& parameters,
        flash::Sw92ProfileCPtFlashBackendOptions backend_options)
        : model(thermodynamics::Sw92Phase<double>::from_parameters(parameters)),
          backend(model, std::move(backend_options)) {}
};

struct CpaBackendGraph {
    thermodynamics::CpaPtPhase model;
    flash::CpaVleEvaluator evaluator;
    flash::CpaPtFlashBackend backend;

    CpaBackendGraph(
        const thermodynamics::CpaParameterSet& parameters,
        thermodynamics::CpaPtOptions phase_options,
        flash::CpaPtFlashBackendOptions backend_options)
        : model(thermodynamics::CpaPtPhase::from_parameters(parameters)),
          evaluator(model, std::move(phase_options)),
          backend(evaluator, std::move(backend_options)) {}
};

} // namespace

OwnedConfiguredPtBackend make_pr76_configured_backend(
    std::string configured_backend_id,
    const thermodynamics::PrParameterSet& parameters,
    thermodynamics::Pr76RootOptions root_options,
    flash::Pr76PtFlashBackendOptions backend_options) {
    auto graph = std::make_shared<Pr76BackendGraph>(
        parameters, root_options, std::move(backend_options));
    auto& backend = graph->backend;
    return retain_configured_pt_backend(
        std::move(configured_backend_id), std::move(graph), backend);
}

OwnedConfiguredPtBackend make_sw92_configured_backend(
    std::string configured_backend_id,
    const thermodynamics::Sw92ParameterSet& parameters,
    flash::Sw92ProfileCPtFlashBackendOptions backend_options) {
    auto graph = std::make_shared<Sw92BackendGraph>(
        parameters, std::move(backend_options));
    auto& backend = graph->backend;
    return retain_configured_pt_backend(
        std::move(configured_backend_id), std::move(graph), backend);
}

OwnedConfiguredPtBackend make_cpa_configured_backend(
    std::string configured_backend_id,
    const thermodynamics::CpaParameterSet& parameters,
    thermodynamics::CpaPtOptions phase_options,
    flash::CpaPtFlashBackendOptions backend_options) {
    auto graph = std::make_shared<CpaBackendGraph>(
        parameters, std::move(phase_options), std::move(backend_options));
    auto& backend = graph->backend;
    return retain_configured_pt_backend(
        std::move(configured_backend_id), std::move(graph), backend);
}

} // namespace mpmc::pt_process
