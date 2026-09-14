#include <jni.h>

#include <mpmc/flash/pt_phase_transition.hpp>
#include <mpmc/pt_process/parameter_snapshot_supplier.hpp>
#include <mpmc/runtime/pt_service.hpp>

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <memory>
#include <optional>
#include <ostream>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace fl = mpmc::flash;
namespace pp = mpmc::pt_process;
namespace rt = mpmc::runtime;

inline constexpr std::string_view kAndroidBridgeConvention =
    "MPMC/PT/android-product-bridge/v1";

class AndroidPtRuntime final {
public:
    AndroidPtRuntime()
        : bundle_(pp::load_repository_curated_pt_parameter_snapshots_v1()) {
        if (!bundle_.structurally_valid()) {
            throw std::runtime_error(
                "Android PT product bridge: repository snapshot bundle is invalid");
        }
        registrations_.reserve(bundle_.backends.size());
        for (auto& owned : bundle_.backends) {
            if (!owned.backend) {
                throw std::runtime_error(
                    "Android PT product bridge: null configured backend");
            }
            registrations_.push_back(
                {owned.configured_backend_id, owned.backend.get()});
        }
        service_ = std::make_unique<rt::PtService>(
            std::span<const rt::PtServiceBackendRegistration>(
                registrations_.data(), registrations_.size()));
    }

    [[nodiscard]] rt::PtService& service() noexcept { return *service_; }

private:
    pp::PtParameterSnapshotBundle bundle_;
    std::vector<rt::PtServiceBackendRegistration> registrations_;
    std::unique_ptr<rt::PtService> service_;
};

[[nodiscard]] AndroidPtRuntime& android_runtime() {
    static AndroidPtRuntime runtime;
    return runtime;
}

void write_json_string(std::ostream& output, std::string_view value) {
    output << '"';
    for (const char raw : value) {
        const auto character = static_cast<unsigned char>(raw);
        switch (character) {
        case '"':
            output << "\\\"";
            break;
        case '\\':
            output << "\\\\";
            break;
        case '\b':
            output << "\\b";
            break;
        case '\f':
            output << "\\f";
            break;
        case '\n':
            output << "\\n";
            break;
        case '\r':
            output << "\\r";
            break;
        case '\t':
            output << "\\t";
            break;
        default:
            if (character < 0x20U) {
                const auto flags = output.flags();
                const auto fill = output.fill();
                output << "\\u" << std::hex << std::setw(4)
                       << std::setfill('0')
                       << static_cast<unsigned int>(character);
                output.flags(flags);
                output.fill(fill);
            } else {
                output << raw;
            }
        }
    }
    output << '"';
}

void write_json_number(std::ostream& output, double value) {
    if (!std::isfinite(value)) {
        throw std::runtime_error(
            "Android PT product bridge: non-finite number cannot be serialized");
    }
    output << std::setprecision(17) << value;
}

[[nodiscard]] std::string_view transition_support(
    fl::PtPhaseTransitionSupport value) {
    switch (value) {
    case fl::PtPhaseTransitionSupport::detection_only:
        return "detection_only";
    case fl::PtPhaseTransitionSupport::fresh_target_resolve:
        return "fresh_target_resolve";
    }
    throw std::runtime_error("Android PT product bridge: unknown transition support");
}

[[nodiscard]] std::string_view transition_trigger(
    fl::PtPhaseTransitionTrigger value) {
    switch (value) {
    case fl::PtPhaseTransitionTrigger::initial_stability_witness:
        return "initial_stability_witness";
    case fl::PtPhaseTransitionTrigger::final_phase_set_instability:
        return "final_phase_set_instability";
    case fl::PtPhaseTransitionTrigger::phase_disappearance:
        return "phase_disappearance";
    case fl::PtPhaseTransitionTrigger::provider_topology_witness:
        return "provider_topology_witness";
    case fl::PtPhaseTransitionTrigger::provider_boundary_route:
        return "provider_boundary_route";
    }
    throw std::runtime_error("Android PT product bridge: unknown transition trigger");
}

[[nodiscard]] std::string_view transition_resolution(
    fl::PtPhaseTransitionResolution value) {
    switch (value) {
    case fl::PtPhaseTransitionResolution::accepted_target:
        return "accepted_target";
    case fl::PtPhaseTransitionResolution::target_resolve_required:
        return "target_resolve_required";
    case fl::PtPhaseTransitionResolution::target_resolve_failed:
        return "target_resolve_failed";
    case fl::PtPhaseTransitionResolution::candidate_not_accepted:
        return "candidate_not_accepted";
    case fl::PtPhaseTransitionResolution::broader_topology_required:
        return "broader_topology_required";
    case fl::PtPhaseTransitionResolution::indeterminate:
        return "indeterminate";
    }
    throw std::runtime_error(
        "Android PT product bridge: unknown transition resolution");
}

[[nodiscard]] std::string_view service_outcome(rt::PtServiceOutcome value) {
    switch (value) {
    case rt::PtServiceOutcome::accepted:
        return "accepted";
    case rt::PtServiceOutcome::phase_set_unstable:
        return "phase_set_unstable";
    case rt::PtServiceOutcome::indeterminate:
        return "indeterminate";
    case rt::PtServiceOutcome::error:
        break;
    }
    throw std::runtime_error(
        "Android PT product bridge: error outcome has no computation label");
}

[[nodiscard]] std::string_view service_error_code(rt::PtServiceErrorCode value) {
    switch (value) {
    case rt::PtServiceErrorCode::invalid_configured_backend_id:
        return "invalid_configured_backend_id";
    case rt::PtServiceErrorCode::configured_backend_not_found:
        return "configured_backend_not_found";
    case rt::PtServiceErrorCode::invalid_pressure:
        return "invalid_pressure";
    case rt::PtServiceErrorCode::invalid_temperature:
        return "invalid_temperature";
    case rt::PtServiceErrorCode::component_count_mismatch:
        return "component_count_mismatch";
    case rt::PtServiceErrorCode::invalid_component_id:
        return "invalid_component_id";
    case rt::PtServiceErrorCode::duplicate_component:
        return "duplicate_component";
    case rt::PtServiceErrorCode::unknown_component:
        return "unknown_component";
    case rt::PtServiceErrorCode::invalid_mole_fraction:
        return "invalid_mole_fraction";
    case rt::PtServiceErrorCode::composition_not_normalized:
        return "composition_not_normalized";
    case rt::PtServiceErrorCode::backend_rejected_request:
        return "backend_rejected_request";
    case rt::PtServiceErrorCode::backend_contract_violation:
        return "backend_contract_violation";
    case rt::PtServiceErrorCode::backend_execution_failure:
        return "backend_execution_failure";
    }
    throw std::runtime_error("Android PT product bridge: unknown service error code");
}

void write_transition_capability(
    std::ostream& output, const fl::PtPhaseTransitionCapability& value) {
    output << "{\"convention\":";
    write_json_string(output, fl::PtPhaseTransitionCapability::convention);
    output << ",\"edges\":[";
    for (std::size_t index = 0; index < value.edges.size(); ++index) {
        if (index != 0U) { output << ','; }
        const auto& edge = value.edges[index];
        output << "{\"sourcePhaseCount\":" << edge.source_phase_count
               << ",\"targetPhaseCount\":" << edge.target_phase_count
               << ",\"support\":";
        write_json_string(output, transition_support(edge.support));
        output << ",\"requiresFreshTargetSolve\":"
               << (edge.requires_fresh_target_solve ? "true" : "false")
               << '}';
    }
    output << "]}";
}

void write_backend_descriptor(
    std::ostream& output, const rt::PtServiceBackendDescriptor& descriptor) {
    if (!descriptor.structurally_valid()) {
        throw std::runtime_error(
            "Android PT product bridge: invalid capability descriptor");
    }
    const auto& capability = descriptor.capability;
    output << "{\"convention\":";
    write_json_string(output, rt::PtServiceBackendDescriptor::convention);
    output << ",\"configuredBackendId\":";
    write_json_string(output, descriptor.configured_backend_id);
    output << ",\"capability\":{\"convention\":";
    write_json_string(output, fl::PtFlashBackendCapability::convention);
    output << ",\"backendId\":";
    write_json_string(output, capability.backend_id);
    output << ",\"modelProfile\":";
    write_json_string(output, capability.model_profile);
    output << ",\"algorithmProfile\":";
    write_json_string(output, capability.algorithm_profile);
    output << ",\"publicationProfile\":";
    write_json_string(output, capability.publication_profile);
    output << ",\"configurationProfile\":";
    write_json_string(output, capability.configuration_profile);
    output << ",\"datasetId\":";
    write_json_string(output, capability.dataset_id);
    output << ",\"revision\":";
    write_json_string(output, capability.revision);
    output << ",\"componentIds\":[";
    for (std::size_t index = 0; index < capability.component_ids.size(); ++index) {
        if (index != 0U) { output << ','; }
        write_json_string(output, capability.component_ids[index]);
    }
    output << "],\"supportedPhaseCounts\":[";
    for (std::size_t index = 0; index < capability.supported_phase_counts.size();
         ++index) {
        if (index != 0U) { output << ','; }
        output << capability.supported_phase_counts[index];
    }
    output << "],\"scalarSettings\":[";
    for (std::size_t index = 0; index < capability.scalar_settings.size(); ++index) {
        if (index != 0U) { output << ','; }
        const auto& setting = capability.scalar_settings[index];
        output << "{\"id\":";
        write_json_string(output, setting.id);
        output << ",\"value\":";
        write_json_number(output, setting.value);
        output << ",\"unit\":";
        write_json_string(output, setting.unit);
        output << '}';
    }
    output << "],\"transitionCapability\":";
    write_transition_capability(output, capability.transition_capability);
    output << ",\"performsInitialStabilitySearch\":"
           << (capability.performs_initial_stability_search ? "true" : "false")
           << ",\"performsFinalPhaseSetReview\":"
           << (capability.performs_final_phase_set_review ? "true" : "false")
           << ",\"performsBoundaryNeighborResolve\":"
           << (capability.performs_boundary_neighbor_resolve ? "true" : "false")
           << ",\"globalStabilityProven\":"
           << (capability.global_stability_proven ? "true" : "false")
           << ",\"phaseMetadataNamespace\":";
    write_json_string(output, capability.phase_metadata_namespace);
    output << "},\"componentInventory\":{\"convention\":";
    write_json_string(output, rt::PtRuntimeComponentInventory::convention);
    output << ",\"components\":[";
    for (std::size_t index = 0;
         index < descriptor.component_inventory.components.size(); ++index) {
        if (index != 0U) { output << ','; }
        const auto& component = descriptor.component_inventory.components[index];
        output << "{\"componentId\":";
        write_json_string(output, component.component_id);
        output << ",\"feedIndex\":" << component.feed_index << '}';
    }
    output << "]}}";
}

void write_transition_report(
    std::ostream& output, const fl::PtPhaseTransitionReport& report) {
    output << "{\"convention\":";
    write_json_string(output, fl::PtPhaseTransitionReport::convention);
    output << ",\"evidence\":[";
    for (std::size_t index = 0; index < report.evidence.size(); ++index) {
        if (index != 0U) { output << ','; }
        const auto& item = report.evidence[index];
        output << "{\"sourcePhaseCount\":" << item.source_phase_count;
        if (item.target_phase_count.has_value()) {
            output << ",\"targetPhaseCount\":" << *item.target_phase_count;
        }
        output << ",\"trigger\":";
        write_json_string(output, transition_trigger(item.trigger));
        output << ",\"resolution\":";
        write_json_string(output, transition_resolution(item.resolution));
        output << ",\"freshTargetSolveAttempted\":"
               << (item.fresh_target_solve_attempted ? "true" : "false")
               << ",\"targetTopologyClosed\":"
               << (item.target_topology_closed ? "true" : "false")
               << ",\"providerEvidenceProfile\":";
        write_json_string(output, item.provider_evidence_profile);
        output << ",\"diagnostic\":";
        write_json_string(output, item.diagnostic);
        output << '}';
    }
    output << "]}";
}

void write_result(
    std::ostream& output, rt::PtServiceOutcome outcome,
    const rt::PtServiceComputationResult& result) {
    if (!result.structurally_valid()) {
        throw std::runtime_error(
            "Android PT product bridge: invalid service result cannot be serialized");
    }
    output << "{\"outcome\":";
    write_json_string(output, service_outcome(outcome));
    output << ",\"provenance\":{\"backend\":";
    write_backend_descriptor(output, result.provenance.backend);
    output << ",\"backendResultConvention\":";
    write_json_string(output, result.provenance.backend_result_convention);
    output << ",\"phaseSetConvention\":";
    write_json_string(output, result.provenance.phase_set_convention);
    output << ",\"phaseTransitionConvention\":";
    write_json_string(output, result.provenance.phase_transition_convention);
    output << ",\"providerResultConvention\":";
    write_json_string(output, result.provenance.provider_result_convention);
    output << "},\"pressurePa\":";
    write_json_number(output, result.pressure_pa);
    output << ",\"temperatureK\":";
    write_json_number(output, result.temperature_k);
    output << ",\"feed\":[";
    for (std::size_t index = 0; index < result.feed.size(); ++index) {
        if (index != 0U) { output << ','; }
        output << "{\"componentId\":";
        write_json_string(output, result.feed[index].component_id);
        output << ",\"moleFraction\":";
        write_json_number(output, result.feed[index].mole_fraction);
        output << '}';
    }
    output << "],\"phases\":[";
    for (std::size_t phase_index = 0; phase_index < result.phases.size();
         ++phase_index) {
        if (phase_index != 0U) { output << ','; }
        const auto& phase = result.phases[phase_index];
        output << "{\"phaseIndex\":" << phase.phase_index
               << ",\"molePhaseFraction\":";
        write_json_number(output, phase.mole_phase_fraction);
        output << ",\"components\":[";
        for (std::size_t component_index = 0;
             component_index < phase.components.size(); ++component_index) {
            if (component_index != 0U) { output << ','; }
            const auto& component = phase.components[component_index];
            output << "{\"componentId\":";
            write_json_string(output, component.component_id);
            output << ",\"moleFraction\":";
            write_json_number(output, component.mole_fraction);
            output << ",\"lnFugacityCoefficient\":";
            write_json_number(output, component.ln_fugacity_coefficient);
            output << '}';
        }
        output << "],\"providerBranch\":";
        write_json_string(output, std::to_string(phase.provider_branch));
        output << ",\"providerBranchSmooth\":"
               << (phase.provider_branch_smooth ? "true" : "false");
        if (phase.compressibility_factor.has_value()) {
            output << ",\"compressibilityFactor\":";
            write_json_number(output, *phase.compressibility_factor);
        }
        if (phase.provider_metadata.has_value()) {
            output << ",\"providerMetadata\":{\"roleId\":";
            write_json_string(output, phase.provider_metadata->role_id);
            output << ",\"familyId\":";
            write_json_string(output, phase.provider_metadata->family_id);
            output << '}';
        }
        output << '}';
    }
    output << "],\"transitionReport\":";
    write_transition_report(output, result.transition_report);
    output << ",\"globalStabilityProven\":"
           << (result.global_stability_proven ? "true" : "false")
           << ",\"morphologyResolved\":"
           << (result.morphology_resolved ? "true" : "false")
           << ",\"diagnostic\":";
    write_json_string(output, result.diagnostic);
    output << '}';
}

[[nodiscard]] std::string discovery_json() {
    const auto descriptors = android_runtime().service().discover_capabilities();
    std::ostringstream output;
    output << "{\"bridgeConvention\":";
    write_json_string(output, kAndroidBridgeConvention);
    output << ",\"serviceBoundaryConvention\":";
    write_json_string(output, rt::pt_service_boundary_convention);
    output << ",\"capabilityConvention\":";
    write_json_string(output, rt::pt_service_capability_convention);
    output << ",\"backends\":[";
    for (std::size_t index = 0; index < descriptors.size(); ++index) {
        if (index != 0U) { output << ','; }
        write_backend_descriptor(output, descriptors[index]);
    }
    output << "]}";
    return output.str();
}

[[nodiscard]] std::string solve_json(
    std::string configured_backend_id, double pressure_pa, double temperature_k,
    std::vector<std::string> component_ids, std::vector<double> fractions) {
    if (component_ids.size() != fractions.size()) {
        throw std::invalid_argument(
            "Android PT product bridge: component/value array length mismatch");
    }
    rt::PtServiceRequest request;
    request.configured_backend_id = std::move(configured_backend_id);
    request.pressure_pa = pressure_pa;
    request.temperature_k = temperature_k;
    request.feed.reserve(component_ids.size());
    for (std::size_t index = 0; index < component_ids.size(); ++index) {
        request.feed.push_back(
            {std::move(component_ids[index]), fractions[index]});
    }

    const auto response = android_runtime().service().solve(request);
    if (!response.structurally_valid()) {
        throw std::runtime_error(
            "Android PT product bridge: PtService returned invalid response");
    }

    std::ostringstream output;
    output << "{\"bridgeConvention\":";
    write_json_string(output, kAndroidBridgeConvention);
    output << ",\"response\":";
    if (response.outcome == rt::PtServiceOutcome::error) {
        if (!response.error.has_value()) {
            throw std::runtime_error(
                "Android PT product bridge: missing service error payload");
        }
        output << "{\"kind\":\"service_error\",\"error\":{\"code\":";
        write_json_string(output, service_error_code(response.error->code));
        output << ",\"field\":";
        write_json_string(output, response.error->field);
        output << ",\"diagnostic\":";
        write_json_string(output, response.error->diagnostic);
        output << "}}";
    } else {
        if (!response.result.has_value()) {
            throw std::runtime_error(
                "Android PT product bridge: missing computation result payload");
        }
        output << "{\"kind\":\"result\",\"result\":";
        write_result(output, response.outcome, *response.result);
        output << '}';
    }
    output << '}';
    return output.str();
}

[[nodiscard]] std::string java_utf8(JNIEnv* environment, jstring value) {
    if (value == nullptr) {
        throw std::invalid_argument("Android PT product bridge: null Java string");
    }
    const char* raw = environment->GetStringUTFChars(value, nullptr);
    if (raw == nullptr) {
        throw std::runtime_error(
            "Android PT product bridge: unable to access Java UTF-8 string");
    }
    std::string result(raw);
    environment->ReleaseStringUTFChars(value, raw);
    return result;
}

[[nodiscard]] std::vector<std::string> java_strings(
    JNIEnv* environment, jobjectArray values) {
    if (values == nullptr) {
        throw std::invalid_argument("Android PT product bridge: null component array");
    }
    const jsize count = environment->GetArrayLength(values);
    std::vector<std::string> result;
    result.reserve(static_cast<std::size_t>(count));
    for (jsize index = 0; index < count; ++index) {
        auto* value = static_cast<jstring>(environment->GetObjectArrayElement(
            values, index));
        try {
            result.push_back(java_utf8(environment, value));
        } catch (...) {
            if (value != nullptr) { environment->DeleteLocalRef(value); }
            throw;
        }
        environment->DeleteLocalRef(value);
    }
    return result;
}

[[nodiscard]] std::vector<double> java_doubles(
    JNIEnv* environment, jdoubleArray values) {
    if (values == nullptr) {
        throw std::invalid_argument("Android PT product bridge: null fraction array");
    }
    const jsize count = environment->GetArrayLength(values);
    std::vector<double> result(static_cast<std::size_t>(count));
    if (count > 0) {
        environment->GetDoubleArrayRegion(values, 0, count, result.data());
        if (environment->ExceptionCheck()) {
            environment->ExceptionClear();
            throw std::runtime_error(
                "Android PT product bridge: unable to read Java fraction array");
        }
    }
    return result;
}

[[nodiscard]] jstring java_string(JNIEnv* environment, const std::string& value) {
    return environment->NewStringUTF(value.c_str());
}

[[nodiscard]] std::string bridge_failure(std::string_view message) {
    std::ostringstream output;
    output << "{\"bridgeConvention\":";
    write_json_string(output, kAndroidBridgeConvention);
    output << ",\"bridgeError\":";
    write_json_string(output, message);
    output << '}';
    return output.str();
}

} // namespace

extern "C" JNIEXPORT jstring JNICALL
Java_org_mpmc_ptandroid_NativeBridge_discoverJson(
    JNIEnv* environment, jclass) noexcept {
    try {
        return java_string(environment, discovery_json());
    } catch (const std::exception& error) {
        return java_string(environment, bridge_failure(error.what()));
    } catch (...) {
        return java_string(environment, bridge_failure("unknown native discovery failure"));
    }
}

extern "C" JNIEXPORT jstring JNICALL
Java_org_mpmc_ptandroid_NativeBridge_solveJson(
    JNIEnv* environment, jclass, jstring configured_backend_id,
    jdouble pressure_pa, jdouble temperature_k, jobjectArray component_ids,
    jdoubleArray fractions) noexcept {
    try {
        return java_string(
            environment,
            solve_json(java_utf8(environment, configured_backend_id),
                       static_cast<double>(pressure_pa),
                       static_cast<double>(temperature_k),
                       java_strings(environment, component_ids),
                       java_doubles(environment, fractions)));
    } catch (const std::exception& error) {
        return java_string(environment, bridge_failure(error.what()));
    } catch (...) {
        return java_string(environment, bridge_failure("unknown native solve failure"));
    }
}
