#include <mpmc/runtime/pt_service.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <new>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool pt_service_contract_headers();
bool pt_service_headers();

namespace {
namespace fl = mpmc::flash;
namespace rt = mpmc::runtime;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(
            std::string(where.file_name()) + ":" +
            std::to_string(where.line()) + ": " + std::string(message));
    }
}

template <typename Error, typename Function>
Error expect_error(Function&& function) {
    try {
        function();
    } catch (const Error& error) {
        return error;
    }
    throw std::runtime_error("expected exception missing");
}

fl::PtFlashBackendCapability make_capability(
    std::string backend_id = "test/backend/v1") {
    fl::PtFlashBackendCapability result;
    result.backend_id = std::move(backend_id);
    result.model_profile = "test/model/v1";
    result.algorithm_profile = "test/algorithm/v1";
    result.publication_profile = "test/publication/v1";
    result.configuration_profile = "test/configuration/v1";
    result.dataset_id = "test/dataset";
    result.revision = "test-revision";
    result.component_ids = {"methane", "water"};
    result.supported_phase_counts = {1U, 2U, 3U};
    result.scalar_settings = {{"test_setting", 2.5, "dimensionless"}};
    result.transition_capability.edges = {
        {1U, 2U, fl::PtPhaseTransitionSupport::fresh_target_resolve, true},
        {2U, 1U, fl::PtPhaseTransitionSupport::detection_only, true},
        {2U, 3U, fl::PtPhaseTransitionSupport::fresh_target_resolve, true},
        {3U, 2U, fl::PtPhaseTransitionSupport::fresh_target_resolve, true}};
    result.performs_initial_stability_search = true;
    result.performs_final_phase_set_review = true;
    result.performs_boundary_neighbor_resolve = true;
    result.global_stability_proven = false;
    result.phase_metadata_namespace = "test/phase-metadata/v1";
    return result;
}

fl::PtCandidatePhase phase(double fraction, const std::vector<double>& composition,
                           std::size_t branch) {
    fl::PtCandidatePhase result;
    result.mole_phase_fraction = fraction;
    result.composition = composition;
    result.activity.ln_phi.assign(composition.size(), 0.1 * static_cast<double>(branch + 1U));
    result.activity.branch = branch;
    result.activity.smooth = true;
    result.compressibility_factor = 0.8 + 0.05 * static_cast<double>(branch);
    return result;
}

enum class BackendMode {
    accepted,
    indeterminate,
    unstable,
    throw_invalid_argument,
    throw_runtime_error,
    throw_logic_error,
    throw_bad_alloc,
    throw_nonstandard,
    wrong_state
};

class ScriptedBackend final : public fl::PtFlashBackend {
public:
    explicit ScriptedBackend(std::string id = "test/backend/v1")
        : capability_(make_capability(std::move(id))) {}

    [[nodiscard]] const fl::PtFlashBackendCapability& capability()
        const noexcept override {
        return capability_;
    }

    [[nodiscard]] fl::PtFlashBackendResult solve(
        const fl::PtFlashRequest& request) override {
        ++call_count;
        last_request = request;
        switch (mode) {
        case BackendMode::throw_invalid_argument:
            throw std::invalid_argument("provider-specific request rejection");
        case BackendMode::throw_runtime_error:
            throw std::runtime_error("provider execution failed");
        case BackendMode::throw_logic_error:
            throw std::logic_error("provider contract failed");
        case BackendMode::throw_bad_alloc:
            throw std::bad_alloc();
        case BackendMode::throw_nonstandard:
            throw 7;
        case BackendMode::accepted:
        case BackendMode::indeterminate:
        case BackendMode::unstable:
        case BackendMode::wrong_state:
            break;
        }

        fl::PtFlashBackendResult result;
        result.capability = capability_;
        result.solution.capability.maximum_phase_count = 3U;
        result.solution.pressure_pa = request.pressure_pa;
        result.solution.temperature_k = request.temperature_k;
        result.solution.feed = request.feed;
        result.solution.global_stability_proven = false;
        result.provider_result_convention = "test/provider-result/v1";
        result.morphology_resolved = false;

        if (mode == BackendMode::wrong_state) {
            result.solution.pressure_pa += 1.0;
        }

        if (mode == BackendMode::accepted) {
            result.solution.status = fl::PtPhaseSetStatus::accepted;
            result.solution.diagnostic = "accepted by scripted boundary fixture";
        } else if (mode == BackendMode::unstable) {
            result.solution.status = fl::PtPhaseSetStatus::phase_set_unstable;
            result.solution.diagnostic = "final phase-set review found instability";
        } else {
            result.solution.status = fl::PtPhaseSetStatus::indeterminate;
            result.solution.diagnostic = "finite search did not close";
        }

        fl::PtCandidatePhaseSet candidates;
        if (phase_count == 1U) {
            candidates.phases.push_back(phase(1.0, request.feed, 0U));
        } else if (phase_count == 2U) {
            candidates.phases.push_back(phase(0.4, request.feed, 0U));
            candidates.phases.push_back(phase(0.6, request.feed, 1U));
        } else {
            candidates.phases.push_back(phase(0.2, request.feed, 0U));
            candidates.phases.push_back(phase(0.3, request.feed, 1U));
            candidates.phases.push_back(phase(0.5, request.feed, 2U));
        }
        result.solution.candidate_phase_set = std::move(candidates);
        for (std::size_t i = 0; i < phase_count; ++i) {
            result.phase_metadata.push_back(
                {"opaque-role-" + std::to_string(i), "opaque-family"});
        }
        return result;
    }

    BackendMode mode{BackendMode::accepted};
    std::size_t phase_count{2U};
    std::size_t call_count{};
    fl::PtFlashRequest last_request;
    fl::PtFlashBackendCapability capability_;
};

rt::PtServiceRequest valid_request() {
    return {
        "test.primary",
        5.0e6,
        325.0,
        {{"water", 0.6}, {"methane", 0.4}}};
}

rt::PtService make_service(ScriptedBackend& backend) {
    const std::array<rt::PtServiceBackendRegistration, 1> backends{{
        {"test.primary", &backend}}};
    return rt::PtService(backends);
}

void require_error_code(const rt::PtServiceResponse& response,
                        rt::PtServiceErrorCode code) {
    require(response.structurally_valid(), "error response is structurally invalid");
    require(response.outcome == rt::PtServiceOutcome::error &&
                response.error.has_value() && response.error->code == code &&
                !response.result.has_value(),
            "unexpected service error code or shape");
}

void discovery_inventory() {
    ScriptedBackend first;
    // The implementation backend ID may repeat: configured_backend_id is the
    // service selector for distinct parameter/configuration instances.
    ScriptedBackend second;
    const std::array<rt::PtServiceBackendRegistration, 2> backends{{
        {"test.primary", &first},
        {"test.secondary", &second}}};
    const rt::PtService service(backends);

    const auto discovered = service.discover_capabilities();
    require(discovered.size() == 2U, "configured backend discovery count changed");
    require(discovered[0].structurally_valid() &&
                discovered[0].configured_backend_id == "test.primary" &&
                discovered[0].capability.backend_id == "test/backend/v1" &&
                discovered[0].component_inventory.components.size() == 2U &&
                discovered[0].component_inventory.components[0].component_id ==
                    "methane" &&
                discovered[0].component_inventory.components[0].feed_index == 0U &&
                discovered[0].component_inventory.components[1].component_id ==
                    "water" &&
                discovered[0].component_inventory.components[1].feed_index == 1U,
            "runtime component inventory changed");
    require(discovered[0].capability.scalar_settings.size() == 1U &&
                discovered[0].capability.supported_phase_counts ==
                    std::vector<std::size_t>({1U, 2U, 3U}),
            "capability discovery lost configured provenance or phase counts");
    require(discovered[1].configured_backend_id == "test.secondary" &&
                discovered[1].capability.backend_id ==
                    discovered[0].capability.backend_id &&
                service.find_capability("test.secondary") == &discovered[1] &&
                service.find_capability("missing") == nullptr,
            "configured-instance capability lookup mismatch");
}

void request_dispatch_provenance() {
    ScriptedBackend backend;
    auto service = make_service(backend);
    const auto response = service.solve(valid_request());

    require(response.structurally_valid() &&
                response.outcome == rt::PtServiceOutcome::accepted &&
                response.accepted_phase_count() == 2U && response.result.has_value(),
            "accepted service response invalid");
    require(backend.call_count == 1U &&
                backend.last_request.feed == std::vector<double>({0.4, 0.6}),
            "ID-keyed request was not mapped to discovered backend order");

    const auto& result = *response.result;
    require(result.feed[0].component_id == "methane" &&
                result.feed[0].mole_fraction == 0.4 &&
                result.feed[1].component_id == "water" &&
                result.feed[1].mole_fraction == 0.6,
            "result feed is not in canonical runtime inventory order");
    require(result.provenance.backend.capability.backend_id ==
                "test/backend/v1" &&
                result.provenance.backend.configured_backend_id ==
                    "test.primary" &&
                result.provenance.backend.capability.model_profile ==
                    "test/model/v1" &&
                result.provenance.backend.capability.dataset_id ==
                    "test/dataset" &&
                result.provenance.backend.capability.revision ==
                    "test-revision" &&
                result.provenance.backend_result_convention ==
                    fl::PtFlashBackendResult::convention &&
                result.provenance.phase_set_convention ==
                    fl::PtPhaseSetResult::convention &&
                result.provenance.phase_transition_convention ==
                    fl::PtPhaseTransitionReport::convention &&
                result.provenance.provider_result_convention ==
                    "test/provider-result/v1",
            "result provenance chain changed");
    require(result.phases[0].components[0].component_id == "methane" &&
                result.phases[0].components[0].ln_fugacity_coefficient == 0.1 &&
                result.phases[0].provider_metadata.has_value() &&
                result.phases[0].provider_metadata->role_id == "opaque-role-0" &&
                !result.morphology_resolved &&
                !result.global_stability_proven,
            "phase payload or scientific limitations changed");

    auto roundoff_request = valid_request();
    roundoff_request.feed[0].mole_fraction =
        std::nextafter(0.6, 1.0);
    const auto roundoff_response = service.solve(roundoff_request);
    require(roundoff_response.structurally_valid() &&
                roundoff_response.outcome == rt::PtServiceOutcome::accepted &&
                backend.last_request.feed[1] ==
                    roundoff_request.feed[0].mole_fraction &&
                roundoff_response.result->feed[1].mole_fraction ==
                    roundoff_request.feed[0].mole_fraction,
            "service normalized a feed that was inside the roundoff contract");
}

void variable_phase_count() {
    ScriptedBackend backend;
    auto service = make_service(backend);
    for (const std::size_t phase_count : {1U, 2U, 3U}) {
        backend.phase_count = phase_count;
        const auto response = service.solve(valid_request());
        require(response.structurally_valid() &&
                    response.outcome == rt::PtServiceOutcome::accepted &&
                    response.accepted_phase_count() == phase_count &&
                    response.result->phases.size() == phase_count,
                "service encoded phase count in a fixed result shape");
        for (std::size_t i = 0; i < phase_count; ++i) {
            require(response.result->phases[i].phase_index == i,
                    "phase diagnostic index changed");
        }
    }
}

void scientific_outcomes() {
    ScriptedBackend backend;
    ScriptedBackend alternate;
    const std::array<rt::PtServiceBackendRegistration, 2> backends{{
        {"test.primary", &backend},
        {"test.alternate", &alternate}}};
    rt::PtService service(backends);

    backend.mode = BackendMode::indeterminate;
    auto response = service.solve(valid_request());
    require(response.structurally_valid() &&
                response.outcome == rt::PtServiceOutcome::indeterminate &&
                response.result.has_value() && response.result->phases.empty() &&
                !response.error.has_value() &&
                response.result->diagnostic == "finite search did not close" &&
                response.accepted_phase_count() == 0U &&
                backend.call_count == 1U && alternate.call_count == 0U,
            "indeterminate was treated as success phases or service failure");

    backend.mode = BackendMode::unstable;
    response = service.solve(valid_request());
    require(response.structurally_valid() &&
                response.outcome == rt::PtServiceOutcome::phase_set_unstable &&
                response.result.has_value() && response.result->phases.empty() &&
                !response.error.has_value() &&
                response.result->diagnostic ==
                    "final phase-set review found instability" &&
                backend.call_count == 2U && alternate.call_count == 0U,
            "unstable candidate was promoted or collapsed into indeterminate");
}

void request_errors() {
    ScriptedBackend backend;
    auto service = make_service(backend);
    auto request = valid_request();

    request.configured_backend_id.clear();
    require_error_code(service.solve(request),
                       rt::PtServiceErrorCode::invalid_configured_backend_id);
    request = valid_request();
    request.configured_backend_id = "missing";
    require_error_code(service.solve(request),
                       rt::PtServiceErrorCode::configured_backend_not_found);
    request = valid_request();
    request.pressure_pa = 0.0;
    require_error_code(service.solve(request),
                       rt::PtServiceErrorCode::invalid_pressure);
    request = valid_request();
    request.temperature_k = std::numeric_limits<double>::quiet_NaN();
    require_error_code(service.solve(request),
                       rt::PtServiceErrorCode::invalid_temperature);
    request = valid_request();
    request.feed.pop_back();
    require_error_code(service.solve(request),
                       rt::PtServiceErrorCode::component_count_mismatch);
    request = valid_request();
    request.feed[0].component_id.clear();
    require_error_code(service.solve(request),
                       rt::PtServiceErrorCode::invalid_component_id);
    request = valid_request();
    request.feed[0].component_id = " water";
    require_error_code(service.solve(request),
                       rt::PtServiceErrorCode::invalid_component_id);
    request = valid_request();
    request.feed[0].component_id = "methane";
    request.feed[1].component_id = "methane";
    require_error_code(service.solve(request),
                       rt::PtServiceErrorCode::duplicate_component);
    request = valid_request();
    request.feed[0].component_id = "carbon-dioxide";
    require_error_code(service.solve(request),
                       rt::PtServiceErrorCode::unknown_component);
    request = valid_request();
    request.feed[0].mole_fraction = -0.1;
    require_error_code(service.solve(request),
                       rt::PtServiceErrorCode::invalid_mole_fraction);
    request = valid_request();
    request.feed[0].mole_fraction = 0.5;
    request.feed[1].mole_fraction = 0.4;
    require_error_code(service.solve(request),
                       rt::PtServiceErrorCode::composition_not_normalized);

    require(backend.call_count == 0U,
            "invalid service request reached the scientific backend");
}

void backend_errors() {
    ScriptedBackend backend;
    auto service = make_service(backend);

    backend.mode = BackendMode::throw_invalid_argument;
    require_error_code(service.solve(valid_request()),
                       rt::PtServiceErrorCode::backend_rejected_request);
    backend.mode = BackendMode::throw_runtime_error;
    require_error_code(service.solve(valid_request()),
                       rt::PtServiceErrorCode::backend_execution_failure);
    backend.mode = BackendMode::throw_logic_error;
    require_error_code(service.solve(valid_request()),
                       rt::PtServiceErrorCode::backend_contract_violation);
    backend.mode = BackendMode::throw_nonstandard;
    require_error_code(service.solve(valid_request()),
                       rt::PtServiceErrorCode::backend_execution_failure);
    backend.mode = BackendMode::throw_bad_alloc;
    (void)expect_error<std::bad_alloc>([&] {
        (void)service.solve(valid_request());
    });
    backend.mode = BackendMode::wrong_state;
    require_error_code(service.solve(valid_request()),
                       rt::PtServiceErrorCode::backend_contract_violation);

    backend.mode = BackendMode::accepted;
    backend.capability_.revision = "mutated-after-discovery";
    require_error_code(service.solve(valid_request()),
                       rt::PtServiceErrorCode::backend_contract_violation);
}

void contract_guards() {
    ScriptedBackend backend;
    auto service = make_service(backend);

    backend.mode = BackendMode::indeterminate;
    auto invalid_outcome = service.solve(valid_request());
    require(invalid_outcome.structurally_valid(),
            "indeterminate fixture response invalid");
    invalid_outcome.outcome = static_cast<rt::PtServiceOutcome>(999);
    require(!invalid_outcome.structurally_valid(),
            "unknown service outcome passed structural validation");

    auto invalid_error = valid_request();
    invalid_error.pressure_pa = 0.0;
    auto error_response = service.solve(invalid_error);
    require(error_response.structurally_valid() &&
                error_response.error.has_value(),
            "error fixture response invalid");
    error_response.error->code = static_cast<rt::PtServiceErrorCode>(999);
    require(!error_response.structurally_valid(),
            "unknown service error code passed structural validation");

    auto invalid_descriptor = service.discover_capabilities().front();
    invalid_descriptor.configured_backend_id = "contains space";
    require(!invalid_descriptor.structurally_valid(),
            "invalid configured backend ID passed descriptor validation");
}

void configuration_guards() {
    ScriptedBackend backend;

    const std::array<rt::PtServiceBackendRegistration, 1> invalid_id{{
        {"invalid id", &backend}}};
    const auto invalid_id_error =
        expect_error<rt::PtServiceConfigurationError>([&] {
            (void)rt::PtService(invalid_id);
        });
    require(invalid_id_error.code() ==
                rt::PtServiceConfigurationErrorCode::invalid_registration_id,
            "invalid configured backend ID error changed");

    const std::array<rt::PtServiceBackendRegistration, 1> null_backend{{
        {"test.null", nullptr}}};
    const auto& null_error = expect_error<rt::PtServiceConfigurationError>([&] {
        (void)rt::PtService(null_backend);
    });
    require(null_error.code() == rt::PtServiceConfigurationErrorCode::null_backend,
            "null backend configuration error changed");

    ScriptedBackend second("test/backend/second/v1");
    const std::array<rt::PtServiceBackendRegistration, 2> duplicate{{
        {"test.duplicate", &backend},
        {"test.duplicate", &second}}};
    const auto& duplicate_error =
        expect_error<rt::PtServiceConfigurationError>([&] {
            (void)rt::PtService(duplicate);
        });
    require(duplicate_error.code() ==
                rt::PtServiceConfigurationErrorCode::duplicate_configured_backend_id,
            "duplicate backend ID configuration error changed");

    ScriptedBackend malformed;
    malformed.capability_.supported_phase_counts = {0U};
    const std::array<rt::PtServiceBackendRegistration, 1> invalid{{
        {"test.invalid", &malformed}}};
    const auto& invalid_error =
        expect_error<rt::PtServiceConfigurationError>([&] {
            (void)rt::PtService(invalid);
        });
    require(invalid_error.code() ==
                rt::PtServiceConfigurationErrorCode::invalid_backend_capability,
            "invalid capability configuration error changed");

    const std::array<rt::PtServiceBackendRegistration, 1> one{{
        {"test.primary", &backend}}};
    rt::PtServiceLimits limits;
    limits.max_components_per_backend = 1U;
    const auto& component_error =
        expect_error<rt::PtServiceConfigurationError>([&] {
            (void)rt::PtService(one, limits);
        });
    require(component_error.code() ==
                rt::PtServiceConfigurationErrorCode::component_limit,
            "component quota configuration error changed");
}

void headers() {
    require(pt_service_contract_headers(),
            "PT service contract public header probe failed");
    require(pt_service_headers(), "PT service public header probe failed");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"discovery_inventory", discovery_inventory},
    {"request_dispatch_provenance", request_dispatch_provenance},
    {"variable_phase_count", variable_phase_count},
    {"scientific_outcomes", scientific_outcomes},
    {"request_errors", request_errors},
    {"backend_errors", backend_errors},
    {"contract_guards", contract_guards},
    {"configuration_guards", configuration_guards},
    {"headers", headers}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one test name required"); }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout << "[PASS] " << name << '\n';
                return 0;
            }
        }
        throw std::invalid_argument("unknown test");
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
