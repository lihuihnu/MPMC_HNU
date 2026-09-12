#include <mpmc/pt_process/composition_root.hpp>
#include <mpmc/pt_process/configured_backends.hpp>
#include <mpmc/pt_process/json_line_observer.hpp>
#include <mpmc/pt_process/parameter_snapshot_supplier.hpp>
#include <mpmc/pt_process/process_host.hpp>

#include <sstream>
#include <string_view>

bool pt_process_headers() {
    std::ostringstream output;
    mpmc::pt_process::PtJsonLineObserver observer(output);
    (void)observer;
    mpmc::pt_process::PtProcessHostOptions options;
    return !options.structurally_valid() &&
           !mpmc::pt_process::PtProcessTlsIdentity{}.structurally_valid() &&
           mpmc::pt_process::pt_process_max_pem_file_bytes == 1024U * 1024U &&
           mpmc::pt_process::pt_parameter_snapshot_bundle_convention ==
               std::string_view{"MPMC/PT/parameter-snapshot-bundle/v1"};
}
