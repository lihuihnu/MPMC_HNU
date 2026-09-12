#include <mpmc/runtime/pt_service.hpp>

#include <array>
#include <string_view>

bool pt_service_headers() {
    const std::array<mpmc::runtime::PtServiceBackendRegistration, 0> backends{};
    const mpmc::runtime::PtService service(backends);
    return mpmc::runtime::pt_service_boundary_convention ==
               std::string_view{"PT/service-boundary/v1"} &&
           service.discover_capabilities().empty() &&
           service.find_capability("missing") == nullptr;
}
