#include <mpmc/flow/phase_identity_continuation.hpp>
#include <mpmc/flow/thermodynamics_absent_phase_extension.hpp>

bool thermodynamics_absent_phase_extension_header() {
    return !mpmc::flow::
        phase_identity_continuation_convention.empty() &&
        !mpmc::flow::
        thermodynamics_absent_phase_extension_convention.empty();
}
