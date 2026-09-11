#include <mpmc/flash/sw92_profile_c_phase_set.hpp>

#include "../sw92_phase_assigned_pt/physical_sample6.hpp"

#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
namespace fl = mpmc::flash;
namespace sample6 = sw92_profile_c_sample6;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

void require_rejected(
    const fl::Sw92PhaseAssignedBoundaryAwareResult& source,
    std::string_view message) {
    const auto result = fl::project_sw92_profile_c_pt_phase_set(source);
    require(result.solution.status == fl::PtPhaseSetStatus::indeterminate &&
                !result.accepted_phase_set_published() &&
                result.solution.accepted_phase_set() == nullptr &&
                !result.solution.candidate_phase_set,
            message);
}

void convention_provenance_guards() {
    const auto source = fl::solve_sw92_phase_assigned_pt_boundary_aware(
        1.0e7, 350.0, sample6::feed(), sample6::model(), 0.0);
    require(source.status == fl::Sw92PhaseAssignedPtStatus::w_h0_h1_locally_closed &&
                source.base.c1 && source.base.c2a1 && source.base.c2b1 &&
                source.base.c2b2,
            "Sample-6 convention-provenance fixture lost its complete C2 chain");

    auto no_w_adapter = source;
    no_w_adapter.base.no_w.adapter_convention = "tampered/no-w-adapter";
    require_rejected(no_w_adapter,
                     "publication accepted a tampered no-W adapter convention");

    auto family_algorithm = source;
    family_algorithm.base.no_w.hydrocarbon_flash.equilibrium_algorithm =
        "tampered/fixed-family-algorithm";
    require_rejected(family_algorithm,
                     "publication accepted a tampered fixed-family algorithm identity");

    auto c1_primitive = source;
    c1_primitive.base.c1->primitive_convention = "tampered/c1-primitive";
    require_rejected(c1_primitive,
                     "publication accepted a tampered C1 primitive convention");

    auto c2a1_witness = source;
    c2a1_witness.base.c2a1->witness_convention = "tampered/c2a1-witness";
    require_rejected(c2a1_witness,
                     "publication accepted a tampered C2a1 witness convention");

    auto c2b1_primitive = source;
    c2b1_primitive.base.c2b1->primitive_convention = "tampered/c2b1-primitive";
    require_rejected(c2b1_primitive,
                     "publication accepted a tampered C2b1 primitive convention");

    auto c2b2_review = source;
    c2b2_review.base.c2b2->review_convention = "tampered/c2b2-review";
    require_rejected(c2b2_review,
                     "publication accepted a tampered C2b2 review convention");
}

} // namespace

int main() {
    try {
        convention_provenance_guards();
        std::cout << "[PASS] convention_provenance_guards\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
