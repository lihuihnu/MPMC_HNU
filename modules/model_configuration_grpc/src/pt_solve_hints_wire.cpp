#include <mpmc/model_configuration_grpc/model_wire.hpp>

#include <string>
#include <utility>
#include <vector>

namespace mpmc::model_configuration_grpc {
namespace {
namespace mc = ::mpmc::model_configuration;
namespace wire = ::mpmc::model_configuration::v1;

std::vector<double> decode_composition(const wire::PtCompositionStart& input) {
    return {input.composition().begin(), input.composition().end()};
}
void encode_composition(const std::vector<double>& input, wire::PtCompositionStart& output) {
    for (double value : input) { output.add_composition(value); }
}
[[noreturn]] void invalid_shape(std::string field, const char* reason) {
    throw mc::ModelConfigurationError(
        mc::ModelConfigurationErrorCode::invalid_value, std::move(field), reason);
}
} // namespace

mc::PtSolveHints decode_solve_hints(const wire::PtSolveHints& input) {
    mc::PtSolveHints output;
    // Missing proto presence maps to the public DTO's empty version. The public
    // model then applies the same unsupported-version rule as direct C++.
    output.version = input.version();
    output.initial_stability_starts.reserve(
        static_cast<std::size_t>(input.initial_stability_starts_size()));
    for (const auto& start : input.initial_stability_starts()) {
        output.initial_stability_starts.push_back(decode_composition(start));
    }
    output.final_two_phase_stability_starts.reserve(
        static_cast<std::size_t>(input.final_two_phase_stability_starts_size()));
    for (const auto& start : input.final_two_phase_stability_starts()) {
        output.final_two_phase_stability_starts.push_back(decode_composition(start));
    }
    output.three_phase_continuation_starts.reserve(
        static_cast<std::size_t>(input.three_phase_continuation_starts_size()));
    for (int index = 0; index < input.three_phase_continuation_starts_size(); ++index) {
        const auto& source = input.three_phase_continuation_starts(index);
        const std::string prefix =
            "hints.three_phase_continuation_starts[" + std::to_string(index) + "]";
        if (source.compositions_size() != 3) {
            invalid_shape(prefix + ".compositions",
                          "exactly three phase compositions required");
        }
        if (source.phase_fraction_seed_size() != 2) {
            invalid_shape(prefix + ".phase_fraction_seed",
                          "exactly two independent phase fractions required");
        }
        mc::PtThreePhaseContinuationHint hint;
        for (std::size_t phase = 0; phase < hint.compositions.size(); ++phase) {
            hint.compositions[phase] = decode_composition(
                source.compositions(static_cast<int>(phase)));
        }
        hint.phase_fraction_seed = {
            source.phase_fraction_seed(0), source.phase_fraction_seed(1)};
        output.three_phase_continuation_starts.push_back(std::move(hint));
    }
    return output;
}

void encode_solve_hints(const mc::PtSolveHints& input, wire::PtSolveHints& output) {
    output.Clear();
    output.set_version(input.version);
    for (const auto& start : input.initial_stability_starts) {
        encode_composition(start, *output.add_initial_stability_starts());
    }
    for (const auto& start : input.final_two_phase_stability_starts) {
        encode_composition(start, *output.add_final_two_phase_stability_starts());
    }
    for (const auto& hint : input.three_phase_continuation_starts) {
        auto* target = output.add_three_phase_continuation_starts();
        for (const auto& composition : hint.compositions) {
            encode_composition(composition, *target->add_compositions());
        }
        target->add_phase_fraction_seed(hint.phase_fraction_seed[0]);
        target->add_phase_fraction_seed(hint.phase_fraction_seed[1]);
    }
}

} // namespace mpmc::model_configuration_grpc
