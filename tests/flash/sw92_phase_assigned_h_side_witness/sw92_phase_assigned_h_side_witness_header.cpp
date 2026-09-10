#include <mpmc/flash/sw92_phase_assigned_h_side_witness.hpp>

bool sw92_phase_assigned_h_side_witness_header() {
    mpmc::flash::Sw92PhaseAssignedHSideWitnessOptions options;
    return options.log_composition_separation > 0.0;
}
