#include <jni.h>

#include <mpmc/pt_process/parameter_snapshot_supplier.hpp>
#include <mpmc/runtime/pt_service.hpp>

#include <array>
#include <cstddef>
#include <exception>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace rt = mpmc::runtime;
namespace pp = mpmc::pt_process;

constexpr std::string_view kPr76Id =
    "pr76.methane-ethane-propane.literature-r1";
constexpr std::string_view kSw92Id =
    "sw92.carbon-dioxide-water.freshwater.literature-r1";
constexpr std::string_view kCpaId =
    "cpa.methanol-water-333.15k.cr1.literature-r1";

[[nodiscard]] std::size_t solve_accepted(
    rt::PtService& service, std::string configured_backend_id,
    double pressure_pa, double temperature_k,
    std::vector<rt::PtServiceCompositionEntry> feed) {
    rt::PtServiceRequest request;
    request.configured_backend_id = std::move(configured_backend_id);
    request.pressure_pa = pressure_pa;
    request.temperature_k = temperature_k;
    request.feed = std::move(feed);

    const auto response = service.solve(request);
    if (!response.structurally_valid()) {
        throw std::runtime_error("PtService returned a structurally invalid response");
    }
    if (response.outcome != rt::PtServiceOutcome::accepted ||
        !response.result.has_value() || response.accepted_phase_count() == 0U) {
        throw std::runtime_error(
            "repository-curated Android smoke did not publish an accepted phase set");
    }
    if (response.result->provenance.backend.configured_backend_id !=
        request.configured_backend_id) {
        throw std::runtime_error("PtService provenance selected the wrong backend");
    }
    return response.accepted_phase_count();
}

[[nodiscard]] std::string run_repository_curated_smoke() {
    auto bundle = pp::load_repository_curated_pt_parameter_snapshots_v1();
    if (!bundle.structurally_valid()) {
        throw std::runtime_error(
            "repository-curated PT snapshot bundle is structurally invalid");
    }

    std::vector<rt::PtServiceBackendRegistration> registrations;
    registrations.reserve(bundle.backends.size());
    for (auto& owned : bundle.backends) {
        if (!owned.backend) {
            throw std::runtime_error("repository-curated backend ownership is null");
        }
        registrations.push_back(
            {owned.configured_backend_id, owned.backend.get()});
    }

    rt::PtService service(std::span<const rt::PtServiceBackendRegistration>(
        registrations.data(), registrations.size()));
    const auto descriptors = service.discover_capabilities();
    if (descriptors.size() != 3U) {
        throw std::runtime_error("Android PtService did not discover three backends");
    }

    std::set<std::string, std::less<>> discovered;
    for (const auto& descriptor : descriptors) {
        if (!descriptor.structurally_valid()) {
            throw std::runtime_error("Android PtService discovered an invalid capability");
        }
        discovered.insert(descriptor.configured_backend_id);
    }
    const std::array<std::string_view, 3> expected{kPr76Id, kSw92Id, kCpaId};
    for (const auto id : expected) {
        if (!discovered.contains(id)) {
            throw std::runtime_error("Android PtService discovery omitted a configured backend");
        }
    }

    const auto pr76_phases = solve_accepted(
        service, std::string(kPr76Id), 1.0e6, 350.0,
        {{"methane", 0.8}, {"ethane", 0.1}, {"propane", 0.1}});
    const auto sw92_phases = solve_accepted(
        service, std::string(kSw92Id), 3.0e6, 340.0,
        {{"carbon-dioxide", 0.7}, {"water", 0.3}});
    const auto cpa_phases = solve_accepted(
        service, std::string(kCpaId), 48852.0, 333.15,
        {{"METHANOL", 0.5}, {"WATER", 0.5}});

    std::ostringstream output;
    output << "ANDROID_JNI_PT_SERVICE_OK"
           << " backends=" << descriptors.size()
           << " pr76_phases=" << pr76_phases
           << " sw92_phases=" << sw92_phases
           << " cpa_phases=" << cpa_phases;
    return output.str();
}

[[nodiscard]] jstring java_string(JNIEnv* environment, const std::string& value) {
    return environment->NewStringUTF(value.c_str());
}

} // namespace

extern "C" JNIEXPORT jstring JNICALL
Java_org_mpmc_ptandroid_NativeBridge_runSmoke(JNIEnv* environment, jclass) noexcept {
    try {
        return java_string(environment, run_repository_curated_smoke());
    } catch (const std::exception& error) {
        return java_string(
            environment,
            std::string("ANDROID_JNI_PT_SERVICE_FAIL native=") + error.what());
    } catch (...) {
        return java_string(environment,
                           "ANDROID_JNI_PT_SERVICE_FAIL native=unknown");
    }
}
