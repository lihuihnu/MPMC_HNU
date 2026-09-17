#include <jni.h>

#include <mpmc/model_configuration/pr76_executable_model.hpp>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
namespace mc = mpmc::model_configuration;
namespace fl = mpmc::flash;

constexpr std::size_t kMaximumRecords = 200000U;

struct WireError final : std::runtime_error {
    int code;
    std::string validation;
    std::string field;
    WireError(int value, std::string message, std::string validation_code = {},
              std::string location = {})
        : std::runtime_error(std::move(message)), code(value),
          validation(std::move(validation_code)), field(std::move(location)) {}
};

std::string java_utf8(JNIEnv* environment, jstring value) {
    if (value == nullptr) throw WireError(3, "null string", "wire.missing_field");
    const char* raw = environment->GetStringUTFChars(value, nullptr);
    if (raw == nullptr) throw std::bad_alloc();
    std::string result(raw);
    environment->ReleaseStringUTFChars(value, raw);
    if (result.find('\0') != std::string::npos) throw WireError(3, "embedded NUL", "wire.unknown_field");
    return result;
}

std::vector<std::string> java_strings(JNIEnv* environment, jobjectArray values) {
    if (values == nullptr) throw WireError(3, "null records", "wire.missing_field");
    const jsize count = environment->GetArrayLength(values);
    if (count < 0 || static_cast<std::size_t>(count) > kMaximumRecords * 3U) {
        throw WireError(8, "record limit", "wire.feed_limit");
    }
    std::vector<std::string> result;
    result.reserve(static_cast<std::size_t>(count));
    for (jsize index = 0; index < count; ++index) {
        auto* item = static_cast<jstring>(environment->GetObjectArrayElement(values, index));
        result.push_back(java_utf8(environment, item));
        environment->DeleteLocalRef(item);
    }
    return result;
}

jobjectArray java_string_array(JNIEnv* environment, const std::vector<std::string>& values) {
    jclass string_class = environment->FindClass("java/lang/String");
    if (string_class == nullptr) return nullptr;
    jobjectArray result = environment->NewObjectArray(
        static_cast<jsize>(values.size()), string_class, nullptr);
    if (result == nullptr) return nullptr;
    for (std::size_t index = 0; index < values.size(); ++index) {
        jstring value = environment->NewStringUTF(values[index].c_str());
        if (value == nullptr) return nullptr;
        environment->SetObjectArrayElement(result, static_cast<jsize>(index), value);
        environment->DeleteLocalRef(value);
    }
    return result;
}

struct FieldValue {
    std::string kind;
    std::string value;
    bool used{false};
};

class InputRecords {
public:
    explicit InputRecords(std::vector<std::string> records) {
        if (records.size() % 3U != 0U || records.size() / 3U > kMaximumRecords) {
            throw WireError(3, "malformed records", "wire.unknown_field");
        }
        for (std::size_t index = 0; index < records.size(); index += 3U) {
            const auto& path = records[index];
            const auto& kind = records[index + 1U];
            if (path.empty() || (kind != "string" && kind != "number" &&
                                 kind != "boolean" && kind != "array")) {
                throw WireError(3, "invalid record", "wire.unknown_field");
            }
            if (!values_.emplace(path, FieldValue{kind, records[index + 2U], false}).second) {
                throw WireError(3, "duplicate record", "wire.unknown_field");
            }
        }
    }

    bool has(std::string_view path) const {
        return values_.find(std::string(path)) != values_.end();
    }

    std::string string(std::string_view path) { return take(path, "string"); }
    std::optional<std::string> optional_string(std::string_view path) {
        return has(path) ? std::optional<std::string>{string(path)} : std::nullopt;
    }
    double number(std::string_view path) {
        const std::string raw = take(path, "number");
        char* end = nullptr;
        const double value = std::strtod(raw.c_str(), &end);
        if (end == raw.c_str() || *end != '\0' || !std::isfinite(value)) {
            throw WireError(3, "invalid finite number", "wire.unknown_field");
        }
        return value;
    }
    std::optional<double> optional_number(std::string_view path) {
        return has(path) ? std::optional<double>{number(path)} : std::nullopt;
    }
    bool boolean(std::string_view path) {
        const std::string raw = take(path, "boolean");
        if (raw == "true") return true;
        if (raw == "false") return false;
        throw WireError(3, "invalid boolean", "wire.unknown_field");
    }
    bool optional_boolean(std::string_view path, bool fallback = false) {
        return has(path) ? boolean(path) : fallback;
    }
    std::size_t array_size(std::string_view path) {
        const std::string raw = take(path, "array");
        std::uint64_t value = 0;
        const auto parsed = std::from_chars(raw.data(), raw.data() + raw.size(), value);
        if (parsed.ec != std::errc{} || parsed.ptr != raw.data() + raw.size() ||
            value > kMaximumRecords) {
            throw WireError(8, "invalid array size", "wire.feed_limit");
        }
        return static_cast<std::size_t>(value);
    }
    std::int64_t integer(std::string_view path) {
        const auto iterator = values_.find(std::string(path));
        if (iterator == values_.end()) throw WireError(3, "missing integer", "wire.missing_field");
        if (iterator->second.kind != "string" && iterator->second.kind != "number") {
            throw WireError(3, "integer kind", "wire.unknown_field");
        }
        iterator->second.used = true;
        const std::string& raw = iterator->second.value;
        std::int64_t value = 0;
        const auto parsed = std::from_chars(raw.data(), raw.data() + raw.size(), value);
        if (parsed.ec != std::errc{} || parsed.ptr != raw.data() + raw.size()) {
            throw WireError(3, "invalid integer", "wire.unknown_field");
        }
        return value;
    }
    std::optional<std::int64_t> optional_integer(std::string_view path) {
        return has(path) ? std::optional<std::int64_t>{integer(path)} : std::nullopt;
    }

    void finish() const {
        const auto unused = std::find_if(values_.begin(), values_.end(),
                                         [](const auto& entry) { return !entry.second.used; });
        if (unused != values_.end()) throw WireError(3, "unknown input field", "wire.unknown_field");
    }

private:
    std::string take(std::string_view path, std::string_view kind) {
        const auto iterator = values_.find(std::string(path));
        if (iterator == values_.end()) throw WireError(3, "missing input field", "wire.missing_field");
        if (iterator->second.kind != kind) throw WireError(3, "input type mismatch", "wire.unknown_field");
        iterator->second.used = true;
        return iterator->second.value;
    }
    std::unordered_map<std::string, FieldValue> values_;
};

class RecordWriter {
public:
    void text(const std::string& path, std::string_view value) { put(path, "string", std::string(value)); }
    void number(const std::string& path, double value) {
        if (!std::isfinite(value)) {
            text(path, std::isnan(value) ? "NaN" : value > 0 ? "Infinity" : "-Infinity");
            return;
        }
        char buffer[128]{};
        const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value,
                                          std::chars_format::general,
                                          std::numeric_limits<double>::max_digits10);
        if (result.ec != std::errc{}) throw std::runtime_error("numeric serialization failed");
        put(path, "number", std::string(buffer, result.ptr));
    }
    void boolean(const std::string& path, bool value) { put(path, "boolean", value ? "true" : "false"); }
    void array(const std::string& path, std::size_t size) { put(path, "array", std::to_string(size)); }
    void uint64_string(const std::string& path, std::uint64_t value) { text(path, std::to_string(value)); }
    void int64_string(const std::string& path, std::int64_t value) { text(path, std::to_string(value)); }
    std::vector<std::string> success() && {
        std::vector<std::string> result;
        result.reserve(values_.size() + 1U);
        result.push_back("ok");
        result.insert(result.end(), std::make_move_iterator(values_.begin()),
                      std::make_move_iterator(values_.end()));
        return result;
    }
private:
    void put(std::string path, std::string kind, std::string value) {
        values_.push_back(std::move(path));
        values_.push_back(std::move(kind));
        values_.push_back(std::move(value));
    }
    std::vector<std::string> values_;
};

std::string indexed(const std::string& path, std::size_t index) {
    return path + "[" + std::to_string(index) + "]";
}
std::string child(const std::string& path, std::string_view name) {
    return path.empty() ? std::string(name) : path + "." + std::string(name);
}

mc::SourceKind source_kind(const std::string& value) {
    if (value == "SOURCE_KIND_LITERATURE") return mc::SourceKind::literature;
    if (value == "SOURCE_KIND_DATABASE") return mc::SourceKind::database;
    if (value == "SOURCE_KIND_USER_SUPPLIED") return mc::SourceKind::user_supplied;
    if (value == "SOURCE_KIND_ASSUMPTION") return mc::SourceKind::assumption;
    if (value == "SOURCE_KIND_SYNTHETIC_TEST") return mc::SourceKind::synthetic_test;
    return mc::SourceKind::unspecified;
}
std::string_view source_kind(mc::SourceKind value) {
    switch (value) {
    case mc::SourceKind::literature: return "SOURCE_KIND_LITERATURE";
    case mc::SourceKind::database: return "SOURCE_KIND_DATABASE";
    case mc::SourceKind::user_supplied: return "SOURCE_KIND_USER_SUPPLIED";
    case mc::SourceKind::assumption: return "SOURCE_KIND_ASSUMPTION";
    case mc::SourceKind::synthetic_test: return "SOURCE_KIND_SYNTHETIC_TEST";
    default: return "SOURCE_KIND_UNSPECIFIED";
    }
}
mc::ComponentKind component_kind(const std::string& value) {
    if (value == "COMPONENT_KIND_PURE") return mc::ComponentKind::pure;
    if (value == "COMPONENT_KIND_PSEUDO") return mc::ComponentKind::pseudo;
    return mc::ComponentKind::unspecified;
}
std::string_view component_kind(mc::ComponentKind value) {
    switch (value) {
    case mc::ComponentKind::pure: return "COMPONENT_KIND_PURE";
    case mc::ComponentKind::pseudo: return "COMPONENT_KIND_PSEUDO";
    default: return "COMPONENT_KIND_UNSPECIFIED";
    }
}

mc::ModelProvenance parse_provenance(InputRecords& input, const std::string& path) {
    mc::ModelProvenance result;
    if (const auto value = input.optional_string(child(path, "kind"))) result.kind = source_kind(*value);
    result.reference = input.optional_string(child(path, "reference")).value_or("");
    result.revision = input.optional_string(child(path, "revision")).value_or("");
    result.locator = input.optional_string(child(path, "locator")).value_or("");
    result.note = input.optional_string(child(path, "note")).value_or("");
    result.acquisition = input.optional_string(child(path, "acquisition")).value_or("");
    result.usage_terms = input.optional_string(child(path, "usageTerms")).value_or("");
    return result;
}

std::optional<mc::ModelScalar> parse_scalar(InputRecords& input, const std::string& path) {
    if (!input.has(child(path, "value"))) return std::nullopt;
    mc::ModelScalar result;
    result.value = input.number(child(path, "value"));
    result.provenance = parse_provenance(input, child(path, "provenance"));
    result.original_unit = input.optional_string(child(path, "originalUnit")).value_or("");
    result.conversion = input.optional_string(child(path, "conversion")).value_or("");
    return result;
}

mc::ThermodynamicModelDefinition parse_definition(InputRecords& input) {
    mc::ThermodynamicModelDefinition result;
    result.version = input.string("definition.version");
    const std::string family = input.string("definition.family");
    result.family = family == "MODEL_FAMILY_PR76"
        ? mc::ThermodynamicModelFamily::peng_robinson_1976
        : mc::ThermodynamicModelFamily::unspecified;
    result.display_name = input.string("definition.displayName");
    result.dataset_id = input.string("definition.datasetId");
    result.revision = input.string("definition.revision");
    result.provenance = parse_provenance(input, "definition.provenance");

    const std::size_t component_count = input.array_size("definition.components");
    result.components.reserve(component_count);
    for (std::size_t index = 0; index < component_count; ++index) {
        const std::string path = indexed("definition.components", index);
        mc::ComponentDefinition component;
        component.component_id = input.string(child(path, "componentId"));
        component.display_name = input.string(child(path, "displayName"));
        component.kind = component_kind(input.string(child(path, "kind")));
        component.provenance = parse_provenance(input, child(path, "provenance"));
        component.molar_mass_kg_per_mol = parse_scalar(input, child(path, "molarMassKgPerMol"));
        result.components.push_back(std::move(component));
    }

    const std::string applicability = "definition.applicability";
    result.applicability.temperature_lower_k = input.optional_number(child(applicability, "temperatureLowerK"));
    result.applicability.temperature_upper_k = input.optional_number(child(applicability, "temperatureUpperK"));
    result.applicability.pressure_lower_pa = input.optional_number(child(applicability, "pressureLowerPa"));
    result.applicability.pressure_upper_pa = input.optional_number(child(applicability, "pressureUpperPa"));
    result.applicability.provenance = parse_provenance(input, child(applicability, "provenance"));
    result.applicability.temperature_lower_exclusive = input.optional_boolean(child(applicability, "temperatureLowerExclusive"));
    result.applicability.temperature_upper_exclusive = input.optional_boolean(child(applicability, "temperatureUpperExclusive"));
    result.applicability.pressure_lower_exclusive = input.optional_boolean(child(applicability, "pressureLowerExclusive"));
    result.applicability.pressure_upper_exclusive = input.optional_boolean(child(applicability, "pressureUpperExclusive"));

    mc::Pr76ParameterDefinition pr;
    const std::size_t pure_count = input.array_size("definition.parameters.pr76.pure");
    pr.pure.reserve(pure_count);
    for (std::size_t index = 0; index < pure_count; ++index) {
        const std::string path = indexed("definition.parameters.pr76.pure", index);
        mc::Pr76PureParameters pure;
        pure.component_id = input.string(child(path, "componentId"));
        pure.critical_temperature_k = parse_scalar(input, child(path, "criticalTemperatureK"));
        pure.critical_pressure_pa = parse_scalar(input, child(path, "criticalPressurePa"));
        pure.acentric_factor = parse_scalar(input, child(path, "acentricFactor"));
        pr.pure.push_back(std::move(pure));
    }
    const std::size_t binary_count = input.array_size("definition.parameters.pr76.binary");
    pr.binary.reserve(binary_count);
    for (std::size_t index = 0; index < binary_count; ++index) {
        const std::string path = indexed("definition.parameters.pr76.binary", index);
        mc::Pr76BinaryInteraction pair;
        pair.first_component_id = input.string(child(path, "firstComponentId"));
        pair.second_component_id = input.string(child(path, "secondComponentId"));
        pair.kij = parse_scalar(input, child(path, "kij"));
        pr.binary.push_back(std::move(pair));
    }
    result.parameters = std::move(pr);
    return result;
}

template <typename Settings>
void parse_stability(InputRecords& input, const std::string& path, Settings& settings) {
    settings.tpd_tolerance = input.optional_number(child(path, "tpdTolerance"));
    settings.stationarity_tolerance = input.optional_number(child(path, "stationarityTolerance"));
    settings.line_search_armijo_coefficient = input.optional_number(child(path, "lineSearchArmijoCoefficient"));
    settings.max_log_composition_step = input.optional_number(child(path, "maxLogCompositionStep"));
    settings.max_iterations = input.optional_integer(child(path, "maxIterations"));
    settings.max_backtracks = input.optional_integer(child(path, "maxBacktracks"));
    settings.max_property_evaluations = input.optional_integer(child(path, "maxPropertyEvaluations"));
    if (input.has(child(path, "automaticMultistart"))) settings.automatic_multistart = input.boolean(child(path, "automaticMultistart"));
    settings.max_starts = input.optional_integer(child(path, "maxStarts"));
}

mc::PtSolverSettings parse_settings(InputRecords& input) {
    if (input.has("presetId")) {
        return mc::resolve_pt_solver_preset(input.string("presetId"));
    }
    mc::PtSolverSettings settings;
    settings.version = input.string("settings.version");
    const std::string kind = input.string("settings.kind");
    settings.kind = kind == "SOLVER_SETTINGS_KIND_CUSTOM"
        ? mc::PtSolverSettingsKind::custom
        : kind == "SOLVER_SETTINGS_KIND_PRESET"
            ? mc::PtSolverSettingsKind::preset : mc::PtSolverSettingsKind::unspecified;
    settings.preset_id = input.optional_string("settings.presetId").value_or("");
    settings.eos_root.max_iterations = input.optional_integer("settings.eosRoot.maxIterations");
    parse_stability(input, "settings.initialStability", settings.initial_stability);
    auto& two = settings.two_phase;
    two.fugacity_equilibrium_tolerance = input.optional_number("settings.twoPhase.fugacityEquilibriumTolerance");
    two.absolute_mass_balance_tolerance = input.optional_number("settings.twoPhase.absoluteMassBalanceTolerance");
    two.relative_mass_balance_tolerance = input.optional_number("settings.twoPhase.relativeMassBalanceTolerance");
    two.minimum_phase_fraction = input.optional_number("settings.twoPhase.minimumPhaseFraction");
    two.minimum_log_composition_separation = input.optional_number("settings.twoPhase.minimumLogCompositionSeparation");
    two.minimum_relative_z_separation = input.optional_number("settings.twoPhase.minimumRelativeZSeparation");
    two.max_log_equilibrium_ratio_step = input.optional_number("settings.twoPhase.maxLogEquilibriumRatioStep");
    two.residual_progress_coefficient = input.optional_number("settings.twoPhase.residualProgressCoefficient");
    two.gibbs_progress_coefficient = input.optional_number("settings.twoPhase.gibbsProgressCoefficient");
    two.max_iterations = input.optional_integer("settings.twoPhase.maxIterations");
    two.max_backtracks = input.optional_integer("settings.twoPhase.maxBacktracks");
    two.max_property_evaluations = input.optional_integer("settings.twoPhase.maxPropertyEvaluations");
    two.rachford_rice_max_iterations = input.optional_integer("settings.twoPhase.rachfordRiceMaxIterations");
    two.max_split_attempts = input.optional_integer("settings.twoPhase.maxSplitAttempts");
    parse_stability(input, "settings.finalTwoPhaseStability", settings.final_two_phase_stability);
    auto& three = settings.three_phase;
    three.chemical_potential_tolerance = input.optional_number("settings.threePhase.chemicalPotentialTolerance");
    three.absolute_mass_balance_tolerance = input.optional_number("settings.threePhase.absoluteMassBalanceTolerance");
    three.relative_mass_balance_tolerance = input.optional_number("settings.threePhase.relativeMassBalanceTolerance");
    three.generalized_rr_balance_tolerance = input.optional_number("settings.threePhase.generalizedRrBalanceTolerance");
    three.minimum_phase_fraction = input.optional_number("settings.threePhase.minimumPhaseFraction");
    three.minimum_log_composition_separation = input.optional_number("settings.threePhase.minimumLogCompositionSeparation");
    three.max_log_step = input.optional_number("settings.threePhase.maxLogStep");
    three.residual_progress_coefficient = input.optional_number("settings.threePhase.residualProgressCoefficient");
    three.max_iterations = input.optional_integer("settings.threePhase.maxIterations");
    three.max_line_search_backtracks = input.optional_integer("settings.threePhase.maxLineSearchBacktracks");
    three.max_balance_iterations = input.optional_integer("settings.threePhase.maxBalanceIterations");
    three.max_balance_backtracks = input.optional_integer("settings.threePhase.maxBalanceBacktracks");
    three.max_property_evaluations = input.optional_integer("settings.threePhase.maxPropertyEvaluations");
    three.max_three_phase_attempts = input.optional_integer("settings.threePhase.maxThreePhaseAttempts");
    three.new_phase_seed_fraction = input.optional_number("settings.threePhase.newPhaseSeedFraction");
    parse_stability(input, "settings.finalThreePhaseStability", settings.final_three_phase_stability);
    return settings;
}

void write_provenance(RecordWriter& output, const std::string& path,
                      const mc::ModelProvenance& value) {
    output.text(child(path, "kind"), source_kind(value.kind));
    output.text(child(path, "reference"), value.reference);
    output.text(child(path, "revision"), value.revision);
    output.text(child(path, "locator"), value.locator);
    output.text(child(path, "note"), value.note);
    output.text(child(path, "acquisition"), value.acquisition);
    output.text(child(path, "usageTerms"), value.usage_terms);
}
void write_scalar(RecordWriter& output, const std::string& path,
                  const std::optional<mc::ModelScalar>& value) {
    if (!value) return;
    output.number(child(path, "value"), value->value);
    write_provenance(output, child(path, "provenance"), value->provenance);
    output.text(child(path, "originalUnit"), value->original_unit);
    output.text(child(path, "conversion"), value->conversion);
}
void write_definition(RecordWriter& output, const std::string& path,
                      const mc::ThermodynamicModelDefinition& value) {
    output.text(child(path, "version"), value.version);
    output.text(child(path, "family"), "MODEL_FAMILY_PR76");
    output.text(child(path, "displayName"), value.display_name);
    output.text(child(path, "datasetId"), value.dataset_id);
    output.text(child(path, "revision"), value.revision);
    write_provenance(output, child(path, "provenance"), value.provenance);
    output.array(child(path, "components"), value.components.size());
    for (std::size_t index = 0; index < value.components.size(); ++index) {
        const std::string item = indexed(child(path, "components"), index);
        const auto& component = value.components[index];
        output.text(child(item, "componentId"), component.component_id);
        output.text(child(item, "displayName"), component.display_name);
        output.text(child(item, "kind"), component_kind(component.kind));
        write_provenance(output, child(item, "provenance"), component.provenance);
        write_scalar(output, child(item, "molarMassKgPerMol"), component.molar_mass_kg_per_mol);
    }
    const auto& applicability = value.applicability;
    const std::string ap = child(path, "applicability");
    if (applicability.temperature_lower_k) output.number(child(ap, "temperatureLowerK"), *applicability.temperature_lower_k);
    if (applicability.temperature_upper_k) output.number(child(ap, "temperatureUpperK"), *applicability.temperature_upper_k);
    if (applicability.pressure_lower_pa) output.number(child(ap, "pressureLowerPa"), *applicability.pressure_lower_pa);
    if (applicability.pressure_upper_pa) output.number(child(ap, "pressureUpperPa"), *applicability.pressure_upper_pa);
    write_provenance(output, child(ap, "provenance"), applicability.provenance);
    output.boolean(child(ap, "temperatureLowerExclusive"), applicability.temperature_lower_exclusive);
    output.boolean(child(ap, "temperatureUpperExclusive"), applicability.temperature_upper_exclusive);
    output.boolean(child(ap, "pressureLowerExclusive"), applicability.pressure_lower_exclusive);
    output.boolean(child(ap, "pressureUpperExclusive"), applicability.pressure_upper_exclusive);
    const auto& pr = std::get<mc::Pr76ParameterDefinition>(value.parameters);
    const std::string pr_path = child(path, "parameters.pr76");
    output.array(child(pr_path, "pure"), pr.pure.size());
    for (std::size_t index = 0; index < pr.pure.size(); ++index) {
        const std::string item = indexed(child(pr_path, "pure"), index);
        output.text(child(item, "componentId"), pr.pure[index].component_id);
        write_scalar(output, child(item, "criticalTemperatureK"), pr.pure[index].critical_temperature_k);
        write_scalar(output, child(item, "criticalPressurePa"), pr.pure[index].critical_pressure_pa);
        write_scalar(output, child(item, "acentricFactor"), pr.pure[index].acentric_factor);
    }
    output.array(child(pr_path, "binary"), pr.binary.size());
    for (std::size_t index = 0; index < pr.binary.size(); ++index) {
        const std::string item = indexed(child(pr_path, "binary"), index);
        output.text(child(item, "firstComponentId"), pr.binary[index].first_component_id);
        output.text(child(item, "secondComponentId"), pr.binary[index].second_component_id);
        write_scalar(output, child(item, "kij"), pr.binary[index].kij);
    }
}

template <typename Settings>
void write_stability(RecordWriter& output, const std::string& path, const Settings& settings) {
    if (settings.tpd_tolerance) output.number(child(path, "tpdTolerance"), *settings.tpd_tolerance);
    if (settings.stationarity_tolerance) output.number(child(path, "stationarityTolerance"), *settings.stationarity_tolerance);
    if (settings.line_search_armijo_coefficient) output.number(child(path, "lineSearchArmijoCoefficient"), *settings.line_search_armijo_coefficient);
    if (settings.max_log_composition_step) output.number(child(path, "maxLogCompositionStep"), *settings.max_log_composition_step);
    if (settings.max_iterations) output.int64_string(child(path, "maxIterations"), *settings.max_iterations);
    if (settings.max_backtracks) output.int64_string(child(path, "maxBacktracks"), *settings.max_backtracks);
    if (settings.max_property_evaluations) output.int64_string(child(path, "maxPropertyEvaluations"), *settings.max_property_evaluations);
    if (settings.automatic_multistart) output.boolean(child(path, "automaticMultistart"), *settings.automatic_multistart);
    if (settings.max_starts) output.int64_string(child(path, "maxStarts"), *settings.max_starts);
}

void write_settings(RecordWriter& output, const std::string& path,
                    const mc::PtSolverSettings& settings) {
    output.text(child(path, "version"), settings.version);
    output.text(child(path, "kind"), settings.kind == mc::PtSolverSettingsKind::preset
        ? "SOLVER_SETTINGS_KIND_PRESET" : "SOLVER_SETTINGS_KIND_CUSTOM");
    if (!settings.preset_id.empty()) output.text(child(path, "presetId"), settings.preset_id);
    if (settings.eos_root.max_iterations) output.int64_string(child(path, "eosRoot.maxIterations"), *settings.eos_root.max_iterations);
    write_stability(output, child(path, "initialStability"), settings.initial_stability);
    const auto& two = settings.two_phase;
    if (two.fugacity_equilibrium_tolerance) output.number(child(path, "twoPhase.fugacityEquilibriumTolerance"), *two.fugacity_equilibrium_tolerance);
    if (two.absolute_mass_balance_tolerance) output.number(child(path, "twoPhase.absoluteMassBalanceTolerance"), *two.absolute_mass_balance_tolerance);
    if (two.relative_mass_balance_tolerance) output.number(child(path, "twoPhase.relativeMassBalanceTolerance"), *two.relative_mass_balance_tolerance);
    if (two.minimum_phase_fraction) output.number(child(path, "twoPhase.minimumPhaseFraction"), *two.minimum_phase_fraction);
    if (two.minimum_log_composition_separation) output.number(child(path, "twoPhase.minimumLogCompositionSeparation"), *two.minimum_log_composition_separation);
    if (two.minimum_relative_z_separation) output.number(child(path, "twoPhase.minimumRelativeZSeparation"), *two.minimum_relative_z_separation);
    if (two.max_log_equilibrium_ratio_step) output.number(child(path, "twoPhase.maxLogEquilibriumRatioStep"), *two.max_log_equilibrium_ratio_step);
    if (two.residual_progress_coefficient) output.number(child(path, "twoPhase.residualProgressCoefficient"), *two.residual_progress_coefficient);
    if (two.gibbs_progress_coefficient) output.number(child(path, "twoPhase.gibbsProgressCoefficient"), *two.gibbs_progress_coefficient);
    if (two.max_iterations) output.int64_string(child(path, "twoPhase.maxIterations"), *two.max_iterations);
    if (two.max_backtracks) output.int64_string(child(path, "twoPhase.maxBacktracks"), *two.max_backtracks);
    if (two.max_property_evaluations) output.int64_string(child(path, "twoPhase.maxPropertyEvaluations"), *two.max_property_evaluations);
    if (two.rachford_rice_max_iterations) output.int64_string(child(path, "twoPhase.rachfordRiceMaxIterations"), *two.rachford_rice_max_iterations);
    if (two.max_split_attempts) output.int64_string(child(path, "twoPhase.maxSplitAttempts"), *two.max_split_attempts);
    write_stability(output, child(path, "finalTwoPhaseStability"), settings.final_two_phase_stability);
    const auto& three = settings.three_phase;
    if (three.chemical_potential_tolerance) output.number(child(path, "threePhase.chemicalPotentialTolerance"), *three.chemical_potential_tolerance);
    if (three.absolute_mass_balance_tolerance) output.number(child(path, "threePhase.absoluteMassBalanceTolerance"), *three.absolute_mass_balance_tolerance);
    if (three.relative_mass_balance_tolerance) output.number(child(path, "threePhase.relativeMassBalanceTolerance"), *three.relative_mass_balance_tolerance);
    if (three.generalized_rr_balance_tolerance) output.number(child(path, "threePhase.generalizedRrBalanceTolerance"), *three.generalized_rr_balance_tolerance);
    if (three.minimum_phase_fraction) output.number(child(path, "threePhase.minimumPhaseFraction"), *three.minimum_phase_fraction);
    if (three.minimum_log_composition_separation) output.number(child(path, "threePhase.minimumLogCompositionSeparation"), *three.minimum_log_composition_separation);
    if (three.max_log_step) output.number(child(path, "threePhase.maxLogStep"), *three.max_log_step);
    if (three.residual_progress_coefficient) output.number(child(path, "threePhase.residualProgressCoefficient"), *three.residual_progress_coefficient);
    if (three.max_iterations) output.int64_string(child(path, "threePhase.maxIterations"), *three.max_iterations);
    if (three.max_line_search_backtracks) output.int64_string(child(path, "threePhase.maxLineSearchBacktracks"), *three.max_line_search_backtracks);
    if (three.max_balance_iterations) output.int64_string(child(path, "threePhase.maxBalanceIterations"), *three.max_balance_iterations);
    if (three.max_balance_backtracks) output.int64_string(child(path, "threePhase.maxBalanceBacktracks"), *three.max_balance_backtracks);
    if (three.max_property_evaluations) output.int64_string(child(path, "threePhase.maxPropertyEvaluations"), *three.max_property_evaluations);
    if (three.max_three_phase_attempts) output.int64_string(child(path, "threePhase.maxThreePhaseAttempts"), *three.max_three_phase_attempts);
    if (three.new_phase_seed_fraction) output.number(child(path, "threePhase.newPhaseSeedFraction"), *three.new_phase_seed_fraction);
    write_stability(output, child(path, "finalThreePhaseStability"), settings.final_three_phase_stability);
}

std::string_view transition_support(fl::PtPhaseTransitionSupport value) {
    return value == fl::PtPhaseTransitionSupport::detection_only
        ? "PT_PHASE_TRANSITION_SUPPORT_DETECTION_ONLY"
        : "PT_PHASE_TRANSITION_SUPPORT_FRESH_TARGET_RESOLVE";
}
std::string_view transition_trigger(fl::PtPhaseTransitionTrigger value) {
    switch (value) {
    case fl::PtPhaseTransitionTrigger::initial_stability_witness: return "PT_PHASE_TRANSITION_TRIGGER_INITIAL_STABILITY_WITNESS";
    case fl::PtPhaseTransitionTrigger::final_phase_set_instability: return "PT_PHASE_TRANSITION_TRIGGER_FINAL_PHASE_SET_INSTABILITY";
    case fl::PtPhaseTransitionTrigger::phase_disappearance: return "PT_PHASE_TRANSITION_TRIGGER_PHASE_DISAPPEARANCE";
    case fl::PtPhaseTransitionTrigger::provider_topology_witness: return "PT_PHASE_TRANSITION_TRIGGER_PROVIDER_TOPOLOGY_WITNESS";
    case fl::PtPhaseTransitionTrigger::provider_boundary_route: return "PT_PHASE_TRANSITION_TRIGGER_PROVIDER_BOUNDARY_ROUTE";
    }
    return "PT_PHASE_TRANSITION_TRIGGER_UNSPECIFIED";
}
std::string_view transition_resolution(fl::PtPhaseTransitionResolution value) {
    switch (value) {
    case fl::PtPhaseTransitionResolution::accepted_target: return "PT_PHASE_TRANSITION_RESOLUTION_ACCEPTED_TARGET";
    case fl::PtPhaseTransitionResolution::target_resolve_required: return "PT_PHASE_TRANSITION_RESOLUTION_TARGET_RESOLVE_REQUIRED";
    case fl::PtPhaseTransitionResolution::target_resolve_failed: return "PT_PHASE_TRANSITION_RESOLUTION_TARGET_RESOLVE_FAILED";
    case fl::PtPhaseTransitionResolution::candidate_not_accepted: return "PT_PHASE_TRANSITION_RESOLUTION_CANDIDATE_NOT_ACCEPTED";
    case fl::PtPhaseTransitionResolution::broader_topology_required: return "PT_PHASE_TRANSITION_RESOLUTION_BROADER_TOPOLOGY_REQUIRED";
    case fl::PtPhaseTransitionResolution::indeterminate: return "PT_PHASE_TRANSITION_RESOLUTION_INDETERMINATE";
    }
    return "PT_PHASE_TRANSITION_RESOLUTION_UNSPECIFIED";
}

void write_transition_capability(RecordWriter& output, const std::string& path,
                                 const fl::PtPhaseTransitionCapability& value) {
    output.text(child(path, "convention"), fl::PtPhaseTransitionCapability::convention);
    output.array(child(path, "edges"), value.edges.size());
    for (std::size_t index = 0; index < value.edges.size(); ++index) {
        const auto& edge = value.edges[index];
        const std::string item = indexed(child(path, "edges"), index);
        output.number(child(item, "sourcePhaseCount"), static_cast<double>(edge.source_phase_count));
        output.number(child(item, "targetPhaseCount"), static_cast<double>(edge.target_phase_count));
        output.text(child(item, "support"), transition_support(edge.support));
        output.boolean(child(item, "requiresFreshTargetSolve"), edge.requires_fresh_target_solve);
    }
}

void write_capability(RecordWriter& output, const std::string& path,
                      const fl::PtFlashBackendCapability& value) {
    output.text(child(path, "convention"), fl::PtFlashBackendCapability::convention);
    output.text(child(path, "backendId"), value.backend_id);
    output.text(child(path, "modelProfile"), value.model_profile);
    output.text(child(path, "algorithmProfile"), value.algorithm_profile);
    output.text(child(path, "publicationProfile"), value.publication_profile);
    output.text(child(path, "configurationProfile"), value.configuration_profile);
    output.text(child(path, "datasetId"), value.dataset_id);
    output.text(child(path, "revision"), value.revision);
    output.array(child(path, "componentIds"), value.component_ids.size());
    for (std::size_t index = 0; index < value.component_ids.size(); ++index) output.text(indexed(child(path, "componentIds"), index), value.component_ids[index]);
    output.array(child(path, "supportedPhaseCounts"), value.supported_phase_counts.size());
    for (std::size_t index = 0; index < value.supported_phase_counts.size(); ++index) output.number(indexed(child(path, "supportedPhaseCounts"), index), static_cast<double>(value.supported_phase_counts[index]));
    output.array(child(path, "scalarSettings"), value.scalar_settings.size());
    for (std::size_t index = 0; index < value.scalar_settings.size(); ++index) {
        const std::string item = indexed(child(path, "scalarSettings"), index);
        output.text(child(item, "id"), value.scalar_settings[index].id);
        output.number(child(item, "value"), value.scalar_settings[index].value);
        output.text(child(item, "unit"), value.scalar_settings[index].unit);
    }
    write_transition_capability(output, child(path, "transitionCapability"), value.transition_capability);
    output.boolean(child(path, "performsInitialStabilitySearch"), value.performs_initial_stability_search);
    output.boolean(child(path, "performsFinalPhaseSetReview"), value.performs_final_phase_set_review);
    output.boolean(child(path, "performsBoundaryNeighborResolve"), value.performs_boundary_neighbor_resolve);
    output.boolean(child(path, "globalStabilityProven"), value.global_stability_proven);
    output.text(child(path, "phaseMetadataNamespace"), value.phase_metadata_namespace);
}

void write_parameter_limits(RecordWriter& output, const mc::ModelConfigurationLimits& value) {
    output.uint64_string("parameterLimits.maxComponents", value.max_components);
    output.uint64_string("parameterLimits.maxPairRecords", value.max_pair_records);
    output.uint64_string("parameterLimits.maxMatrixEntries", value.max_matrix_entries);
    output.uint64_string("parameterLimits.maxIdentifierBytes", value.max_identifier_bytes);
    output.uint64_string("parameterLimits.maxDisplayNameBytes", value.max_display_name_bytes);
    output.uint64_string("parameterLimits.maxProvenanceFieldBytes", value.max_provenance_field_bytes);
    output.uint64_string("parameterLimits.maxTotalTextBytes", value.max_total_text_bytes);
}
void write_solver_limits(RecordWriter& output, const mc::PtSolverSafetyLimits& value) {
    output.uint64_string("solverLimits.maxRootIterations", value.max_root_iterations);
    output.uint64_string("solverLimits.maxIterations", value.max_iterations);
    output.uint64_string("solverLimits.maxBacktracks", value.max_backtracks);
    output.uint64_string("solverLimits.maxPropertyEvaluations", value.max_property_evaluations);
    output.uint64_string("solverLimits.maxStabilityStarts", value.max_stability_starts);
    output.uint64_string("solverLimits.maxSplitAttempts", value.max_split_attempts);
    output.uint64_string("solverLimits.maxThreePhaseAttempts", value.max_three_phase_attempts);
    output.uint64_string("solverLimits.maxComponents", value.max_components);
    output.uint64_string("solverLimits.maxStabilityStartEntries", value.max_stability_start_entries);
    output.uint64_string("solverLimits.maxThreePhaseStarts", value.max_three_phase_starts);
    output.uint64_string("solverLimits.maxThreePhaseStartEntries", value.max_three_phase_start_entries);
}

void write_snapshot(RecordWriter& output, const mc::Pr76ExecutableModel& model) {
    write_definition(output, "definition", model.parameter_snapshot().definition());
    write_settings(output, "settings", model.solver_configuration().settings());
    write_capability(output, "capability", model.capability());
    write_parameter_limits(output, model.parameter_limits());
    write_solver_limits(output, model.solver_configuration().safety_limits());
}

std::string_view result_outcome(fl::PtPhaseSetStatus value) {
    switch (value) {
    case fl::PtPhaseSetStatus::accepted: return "PT_COMPUTATION_OUTCOME_ACCEPTED";
    case fl::PtPhaseSetStatus::phase_set_unstable: return "PT_COMPUTATION_OUTCOME_PHASE_SET_UNSTABLE";
    case fl::PtPhaseSetStatus::indeterminate: return "PT_COMPUTATION_OUTCOME_INDETERMINATE";
    }
    return "PT_COMPUTATION_OUTCOME_UNSPECIFIED";
}

void write_transition_report(RecordWriter& output, const std::string& path,
                             const fl::PtPhaseTransitionReport& report) {
    output.text(child(path, "convention"), fl::PtPhaseTransitionReport::convention);
    output.array(child(path, "evidence"), report.evidence.size());
    for (std::size_t index = 0; index < report.evidence.size(); ++index) {
        const auto& evidence = report.evidence[index];
        const std::string item = indexed(child(path, "evidence"), index);
        output.number(child(item, "sourcePhaseCount"), static_cast<double>(evidence.source_phase_count));
        if (evidence.target_phase_count) output.number(child(item, "targetPhaseCount"), static_cast<double>(*evidence.target_phase_count));
        output.text(child(item, "trigger"), transition_trigger(evidence.trigger));
        output.text(child(item, "resolution"), transition_resolution(evidence.resolution));
        output.boolean(child(item, "freshTargetSolveAttempted"), evidence.fresh_target_solve_attempted);
        output.boolean(child(item, "targetTopologyClosed"), evidence.target_topology_closed);
        output.text(child(item, "providerEvidenceProfile"), evidence.provider_evidence_profile);
        output.text(child(item, "diagnostic"), evidence.diagnostic);
    }
}

void write_result(RecordWriter& output, const fl::PtFlashBackendResult& result) {
    output.text("backendResultConvention", fl::PtFlashBackendResult::convention);
    output.text("phaseSetConvention", fl::PtPhaseSetResult::convention);
    write_capability(output, "capability", result.capability);
    output.text("outcome", result_outcome(result.solution.status));
    output.number("maximumPhaseCount", static_cast<double>(result.capability.maximum_phase_count()));
    output.number("pressurePa", result.solution.pressure_pa);
    output.number("temperatureK", result.solution.temperature_k);
    output.array("feed", result.solution.feed.size());
    for (std::size_t index = 0; index < result.solution.feed.size(); ++index) output.number(indexed("feed", index), result.solution.feed[index]);
    if (result.solution.candidate_phase_set) {
        output.array("candidatePhaseSet.phases", result.solution.candidate_phase_set->phases.size());
        for (std::size_t index = 0; index < result.solution.candidate_phase_set->phases.size(); ++index) {
            const auto& phase = result.solution.candidate_phase_set->phases[index];
            const std::string item = indexed("candidatePhaseSet.phases", index);
            output.number(child(item, "molePhaseFraction"), phase.mole_phase_fraction);
            output.array(child(item, "composition"), phase.composition.size());
            for (std::size_t component = 0; component < phase.composition.size(); ++component) output.number(indexed(child(item, "composition"), component), phase.composition[component]);
            output.array(child(item, "lnFugacityCoefficient"), phase.activity.ln_phi.size());
            for (std::size_t component = 0; component < phase.activity.ln_phi.size(); ++component) output.number(indexed(child(item, "lnFugacityCoefficient"), component), phase.activity.ln_phi[component]);
            output.uint64_string(child(item, "providerBranch"), phase.activity.branch);
            output.boolean(child(item, "providerBranchSmooth"), phase.activity.smooth);
            if (phase.compressibility_factor) output.number(child(item, "compressibilityFactor"), *phase.compressibility_factor);
        }
    }
    output.boolean("globalStabilityProven", result.solution.global_stability_proven);
    output.text("diagnostic", result.solution.diagnostic);
    write_transition_report(output, "transitionReport", result.transition_report);
    output.text("providerResultConvention", result.provider_result_convention);
    output.array("phaseMetadata", result.phase_metadata.size());
    for (std::size_t index = 0; index < result.phase_metadata.size(); ++index) {
        const std::string item = indexed("phaseMetadata", index);
        output.text(child(item, "roleId"), result.phase_metadata[index].role_id);
        output.text(child(item, "familyId"), result.phase_metadata[index].family_id);
    }
    output.boolean("morphologyResolved", result.morphology_resolved);
}

fl::PtFlashRequest parse_solve(InputRecords& input) {
    fl::PtFlashRequest request;
    request.pressure_pa = input.number("pressurePa");
    request.temperature_k = input.number("temperatureK");
    const std::size_t count = input.array_size("feed");
    request.feed.reserve(count);
    for (std::size_t index = 0; index < count; ++index) request.feed.push_back(input.number(indexed("feed", index)));
    input.finish();
    return request;
}

std::string config_code(mc::ModelConfigurationErrorCode code) {
    using C = mc::ModelConfigurationErrorCode;
    switch (code) {
    case C::unsupported_version: return "unsupported_version";
    case C::unsupported_family: return "unsupported_family";
    case C::missing_field: return "missing_field";
    case C::invalid_identifier: return "invalid_identifier";
    case C::duplicate_identifier: return "duplicate_identifier";
    case C::unknown_component: return "unknown_component";
    case C::invalid_value: return "invalid_value";
    case C::invalid_unit: return "invalid_unit";
    case C::invalid_source: return "invalid_source";
    case C::invalid_range: return "invalid_range";
    case C::duplicate_parameter: return "duplicate_parameter";
    case C::missing_parameter: return "missing_parameter";
    case C::invalid_pair: return "invalid_pair";
    case C::resource_limit: return "resource_limit";
    case C::unsupported_preset: return "unsupported_preset";
    case C::invalid_settings: return "invalid_settings";
    }
    return "invalid_settings";
}
int config_status(mc::ModelConfigurationErrorCode code) {
    using C = mc::ModelConfigurationErrorCode;
    if (code == C::unsupported_version || code == C::unsupported_family || code == C::unsupported_preset) return 12;
    if (code == C::resource_limit) return 8;
    return 3;
}
std::vector<std::string> failure(int code, std::string reason,
                                 std::string validation = {}, std::string field = {}) {
    if (validation.empty()) return {"error", std::to_string(code), std::move(reason)};
    return {"error", std::to_string(code), std::move(reason), std::move(validation), std::move(field)};
}
std::vector<std::string> config_failure(const mc::ModelConfigurationError& error) {
    return failure(config_status(error.code()), "workbench.invalid_input",
                   "configuration." + config_code(error.code()), error.field());
}
std::vector<std::string> request_failure(const std::exception& error) {
    const auto* located = dynamic_cast<const mc::Pr76SolveRequestError*>(&error);
    return failure(3, "workbench.invalid_input", "request.rejected",
                   located == nullptr ? std::string{} : located->field());
}

struct WorkbenchState {
    std::mutex mutex;
    std::shared_ptr<mc::Pr76ExecutableModel> model;
    std::string active_request;
    std::string published_by;
    bool cancelled{false};
};
WorkbenchState& workbench() {
    static WorkbenchState state;
    return state;
}

void begin_request(const std::string& id, bool retire_current) {
    auto& state = workbench();
    std::lock_guard lock(state.mutex);
    if (!state.active_request.empty()) throw WireError(8, "busy");
    state.published_by.clear();
    state.active_request = id;
    state.cancelled = false;
    if (retire_current) state.model.reset();
}
void abandon_request(const std::string& id, bool clear_model) {
    auto& state = workbench();
    std::lock_guard lock(state.mutex);
    if (state.active_request == id) {
        state.active_request.clear();
        state.cancelled = false;
    }
    if (clear_model) state.model.reset();
}
bool publish_model(const std::string& id, std::shared_ptr<mc::Pr76ExecutableModel> model) {
    auto& state = workbench();
    std::lock_guard lock(state.mutex);
    if (state.active_request != id || state.cancelled) {
        if (state.active_request == id) state.active_request.clear();
        state.cancelled = false;
        state.model.reset();
        return false;
    }
    state.model = std::move(model);
    state.published_by = id;
    state.active_request.clear();
    return true;
}
std::shared_ptr<mc::Pr76ExecutableModel> current_model(const std::string& id) {
    auto& state = workbench();
    std::lock_guard lock(state.mutex);
    if (state.active_request != id || state.cancelled) return {};
    return state.model;
}
bool finish_solve(const std::string& id) {
    auto& state = workbench();
    std::lock_guard lock(state.mutex);
    if (state.active_request != id || state.cancelled) {
        state.model.reset();
        if (state.active_request == id) state.active_request.clear();
        state.cancelled = false;
        return false;
    }
    state.active_request.clear();
    return true;
}

std::vector<std::string> apply_records(const std::string& id,
                                       std::vector<std::string> records) {
    try {
        InputRecords input(std::move(records));
        const auto definition = parse_definition(input);
        const auto settings = parse_settings(input);
        input.finish();
        begin_request(id, true);
        try {
            auto next = std::make_shared<mc::Pr76ExecutableModel>(definition, settings);
            RecordWriter output;
            write_snapshot(output, *next);
            if (!publish_model(id, next)) return failure(1, "workbench.document_changed");
            return std::move(output).success();
        } catch (const mc::ModelConfigurationError& error) {
            abandon_request(id, true);
            return config_failure(error);
        } catch (...) {
            abandon_request(id, true);
            throw;
        }
    } catch (const WireError& error) {
        return failure(error.code, error.code == 8 ? "workbench.busy" : "workbench.invalid_input",
                       error.validation, error.field);
    } catch (...) {
        return failure(13, "workbench.failed");
    }
}

std::vector<std::string> solve_records(const std::string& id,
                                       std::vector<std::string> records) {
    try {
        InputRecords input(std::move(records));
        const auto request = parse_solve(input);
        begin_request(id, false);
        auto model = current_model(id);
        if (!model) {
            abandon_request(id, false);
            return failure(5, "workbench.no_model");
        }
        try {
            const auto result = model->solve(request);
            RecordWriter output;
            write_result(output, result);
            if (!finish_solve(id)) return failure(1, "workbench.document_changed");
            return std::move(output).success();
        } catch (const mc::Pr76ModelBusyError&) {
            abandon_request(id, false);
            return failure(8, "workbench.busy");
        } catch (const std::invalid_argument& error) {
            abandon_request(id, false);
            return request_failure(error);
        } catch (const std::domain_error& error) {
            abandon_request(id, false);
            return request_failure(error);
        } catch (const std::length_error& error) {
            abandon_request(id, false);
            return request_failure(error);
        } catch (...) {
            abandon_request(id, true);
            throw;
        }
    } catch (const WireError& error) {
        return failure(error.code, error.code == 8 ? "workbench.busy" : "workbench.invalid_input",
                       error.validation, error.field);
    } catch (...) {
        return failure(13, "workbench.failed");
    }
}

std::vector<std::string> release_records(const std::string& id) {
    try {
        begin_request(id, true);
        abandon_request(id, false);
        return {"ok-null"};
    } catch (const WireError& error) {
        return failure(error.code, error.code == 8 ? "workbench.busy" : "workbench.invalid_input");
    } catch (...) {
        return failure(13, "workbench.failed");
    }
}

void cancel_request(const std::string& id) {
    auto& state = workbench();
    std::lock_guard lock(state.mutex);
    if (state.active_request == id) {
        state.cancelled = true;
        state.model.reset();
    } else if (state.published_by == id) {
        state.model.reset();
        state.published_by.clear();
    }
}

} // namespace

extern "C" JNIEXPORT jobjectArray JNICALL
Java_org_mpmc_ptandroid_NativeBridge_modelApplyRecords(
    JNIEnv* environment, jclass, jstring request_id, jobjectArray records) noexcept {
    try {
        return java_string_array(environment,
            apply_records(java_utf8(environment, request_id), java_strings(environment, records)));
    } catch (...) {
        return java_string_array(environment, failure(13, "workbench.failed"));
    }
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_org_mpmc_ptandroid_NativeBridge_modelSolveRecords(
    JNIEnv* environment, jclass, jstring request_id, jobjectArray records) noexcept {
    try {
        return java_string_array(environment,
            solve_records(java_utf8(environment, request_id), java_strings(environment, records)));
    } catch (...) {
        return java_string_array(environment, failure(13, "workbench.failed"));
    }
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_org_mpmc_ptandroid_NativeBridge_modelReleaseRecords(
    JNIEnv* environment, jclass, jstring request_id) noexcept {
    try {
        return java_string_array(environment, release_records(java_utf8(environment, request_id)));
    } catch (...) {
        return java_string_array(environment, failure(13, "workbench.failed"));
    }
}

extern "C" JNIEXPORT void JNICALL
Java_org_mpmc_ptandroid_NativeBridge_modelCancel(
    JNIEnv* environment, jclass, jstring request_id) noexcept {
    try { cancel_request(java_utf8(environment, request_id)); }
    catch (...) { /* Best-effort cancellation never throws across JNI. */ }
}
