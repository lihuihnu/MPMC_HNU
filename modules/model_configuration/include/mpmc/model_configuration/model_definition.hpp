#ifndef MPMC_MODEL_CONFIGURATION_MODEL_DEFINITION_HPP
#define MPMC_MODEL_CONFIGURATION_MODEL_DEFINITION_HPP

#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace mpmc::model_configuration {

// Parameter-domain version only. This is not a solve request or the future
// complete model-creation protocol (resolved solver settings are a later gate).
inline constexpr std::string_view model_parameter_definition_v1 =
    "thermodynamic-model/parameter-definition/v1";

enum class ThermodynamicModelFamily {
    unspecified, peng_robinson_1976, soreide_whitson_1992, cubic_plus_association
};
enum class ComponentKind { unspecified, pure, pseudo };
enum class SourceKind {
    unspecified, literature, database, user_supplied, assumption, synthetic_test
};
// Host policy, never a client-controlled field in the definition.
enum class ModelDataPolicy { ordinary, allow_synthetic_tests };

struct ModelProvenance {
    SourceKind kind = SourceKind::unspecified;
    std::string reference;
    std::string revision;
    std::string locator;
    std::string note;
    std::string acquisition;
    std::string usage_terms;
};

// The containing field specifies the canonical SI unit. Original-unit metadata
// records history only; neither conversion nor unit inference happens here.
struct ModelScalar {
    double value = std::numeric_limits<double>::quiet_NaN();
    ModelProvenance provenance;
    std::string original_unit;
    std::string conversion;
};

struct ComponentDefinition {
    std::string component_id; // Stable case-sensitive key; not an array index.
    std::string display_name;
    ComponentKind kind = ComponentKind::unspecified;
    ModelProvenance provenance; // Includes pseudo-component characterization.
    std::optional<ModelScalar> molar_mass_kg_per_mol;
};

struct ModelApplicability {
    std::optional<double> temperature_lower_k;
    std::optional<double> temperature_upper_k;
    std::optional<double> pressure_lower_pa;
    std::optional<double> pressure_upper_pa;
    ModelProvenance provenance;
    // Absence means UNKNOWN. The first PR adapter requires both endpoints of
    // each declared interval; it rejects a one-sided range rather than inventing
    // an endpoint unsupported by the existing thermodynamics contract.
};

struct Pr76PureParameters {
    std::string component_id;
    std::optional<ModelScalar> critical_temperature_k;
    std::optional<ModelScalar> critical_pressure_pa;
    std::optional<ModelScalar> acentric_factor;
};
struct Pr76BinaryInteraction {
    std::string first_component_id;
    std::string second_component_id;
    std::optional<ModelScalar> kij;
};
struct Pr76ParameterDefinition {
    std::vector<Pr76PureParameters> pure;
    // Exactly one explicit value for each unordered pair. No self-pair records.
    std::vector<Pr76BinaryInteraction> binary;
};

// Owning, unvalidated draft. Only PR76 parameters can currently be prepared.
// SW92/CPA families are reserved; their existing configured solves are separate.
// Future families get their own variant arm, not PR-only pair semantics.
struct ThermodynamicModelDefinition {
    std::string version; // Required exact version; no implicit current version.
    ThermodynamicModelFamily family = ThermodynamicModelFamily::unspecified;
    std::string display_name;
    std::string dataset_id;
    std::string revision;
    ModelProvenance provenance;
    std::vector<ComponentDefinition> components; // Explicit canonical model order.
    ModelApplicability applicability;
    std::variant<std::monostate, Pr76ParameterDefinition> parameters;
};

// Server-owned limits. Bound shape and text before any mapping/copy or dense
// matrix construction. Transports must separately bound serialized bytes BEFORE
// parsing; max_total_text_bytes is not a serialized-payload or allocator quota.
struct ModelConfigurationLimits {
    std::size_t max_components{256};
    std::size_t max_pair_records{32640};
    std::size_t max_matrix_entries{65536};
    std::size_t max_identifier_bytes{128};
    std::size_t max_display_name_bytes{512};
    std::size_t max_provenance_field_bytes{4096};
    std::size_t max_total_text_bytes{1024U * 1024U};
};

enum class ModelConfigurationErrorCode {
    unsupported_version, unsupported_family, missing_field, invalid_identifier,
    duplicate_identifier, unknown_component, invalid_value, invalid_unit,
    invalid_source, invalid_range, duplicate_parameter, missing_parameter,
    invalid_pair, resource_limit
};
class ModelConfigurationError : public std::invalid_argument {
public:
    ModelConfigurationError(ModelConfigurationErrorCode code, std::string field,
                            std::string_view reason)
        : std::invalid_argument(field + ": " + std::string(reason)), code_(code),
          field_(std::move(field)) {}
    [[nodiscard]] ModelConfigurationErrorCode code() const noexcept { return code_; }
    [[nodiscard]] const std::string& field() const & noexcept { return field_; }
    const std::string& field() const && = delete;
private:
    ModelConfigurationErrorCode code_;
    std::string field_;
};

} // namespace mpmc::model_configuration
#endif // MPMC_MODEL_CONFIGURATION_MODEL_DEFINITION_HPP
