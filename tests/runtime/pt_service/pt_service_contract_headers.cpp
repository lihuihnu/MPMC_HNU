#include <mpmc/runtime/pt_service_contract.hpp>

#include <string_view>

bool pt_service_contract_headers() {
    mpmc::runtime::PtRuntimeComponentInventory inventory;
    mpmc::runtime::PtServiceBackendDescriptor descriptor;
    mpmc::runtime::PtServiceRequest request;
    mpmc::runtime::PtServiceResponse response;
    return mpmc::runtime::PtRuntimeComponentInventory::convention ==
               std::string_view{"PT/service-component-inventory/v1"} &&
           mpmc::runtime::PtServiceBackendDescriptor::convention ==
               std::string_view{"PT/service-capability-discovery/v1"} &&
           mpmc::runtime::PtServiceRequest::convention ==
               std::string_view{"PT/service-request/v1"} &&
           mpmc::runtime::PtServiceResponse::convention ==
               std::string_view{"PT/service-result/v1"} &&
           !inventory.structurally_valid() &&
           !descriptor.structurally_valid() &&
           request.configured_backend_id.empty() &&
           !response.structurally_valid();
}
