#ifndef MPMC_MODEL_CONFIGURATION_PR76_PARAMETERS_HPP
#define MPMC_MODEL_CONFIGURATION_PR76_PARAMETERS_HPP

#include <mpmc/model_configuration/model_definition.hpp>
#include <mpmc/thermodynamics/pr_parameters.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

namespace mpmc::model_configuration {
namespace detail {

inline void model_require(bool condition, ModelConfigurationErrorCode code,
                          const std::string& field, std::string_view reason) {
    if (!condition) { throw ModelConfigurationError(code, field, reason); }
}

class ModelTextBudget {
public:
    explicit ModelTextBudget(const ModelConfigurationLimits& limits)
        : limits_(limits), remaining_(limits.max_total_text_bytes) {}
    void text(const std::string& value, std::size_t maximum, const std::string& field) {
        model_require(value.size() <= maximum && value.size() <= remaining_,
                      ModelConfigurationErrorCode::resource_limit, field,
                      "configuration text quota exceeded");
        remaining_ -= value.size();
    }
    void id(const std::string& value, const std::string& field) {
        text(value, limits_.max_identifier_bytes, field);
    }
    void label(const std::string& value, const std::string& field) {
        text(value, limits_.max_display_name_bytes, field);
    }
    void provenance(const ModelProvenance& source, const std::string& field) {
        const auto add = [&](const std::string& value, const char* name) {
            text(value, limits_.max_provenance_field_bytes, field + "." + name);
        };
        add(source.reference, "reference");
        add(source.revision, "revision");
        add(source.locator, "locator");
        add(source.note, "note");
        add(source.acquisition, "acquisition");
        add(source.usage_terms, "usage_terms");
    }
    void scalar(const std::optional<ModelScalar>& value, const std::string& field) {
        if (!value) { return; }
        provenance(value->provenance, field + ".provenance");
        text(value->original_unit, limits_.max_provenance_field_bytes, field + ".original_unit");
        text(value->conversion, limits_.max_provenance_field_bytes, field + ".conversion");
    }
private:
    const ModelConfigurationLimits& limits_;
    std::size_t remaining_;
};

inline void bound_model_definition(const ThermodynamicModelDefinition& definition,
                                   const ModelConfigurationLimits& limits) {
    const auto* pr = std::get_if<Pr76ParameterDefinition>(&definition.parameters);
    model_require(definition.version == model_parameter_definition_v1,
                  ModelConfigurationErrorCode::unsupported_version, "version",
                  "unsupported parameter-definition version");
    model_require(definition.family == ThermodynamicModelFamily::peng_robinson_1976,
                  ModelConfigurationErrorCode::unsupported_family, "family",
                  "custom parameters are currently supported only for Peng-Robinson 1976");
    model_require(pr != nullptr, ModelConfigurationErrorCode::missing_field,
                  "parameters", "Peng-Robinson 1976 parameters are required");
    const auto n = definition.components.size();
    model_require(n != 0, ModelConfigurationErrorCode::missing_field,
                  "components", "at least one component is required");
    model_require(n <= limits.max_components && pr->pure.size() <= limits.max_components &&
                      pr->binary.size() <= limits.max_pair_records &&
                      n <= limits.max_matrix_entries / n,
                  ModelConfigurationErrorCode::resource_limit, "parameters",
                  "component, pure-record, pair or matrix quota exceeded");
    ModelTextBudget budget(limits);
    budget.id(definition.version, "version");
    budget.label(definition.display_name, "display_name");
    budget.id(definition.dataset_id, "dataset_id");
    budget.id(definition.revision, "revision");
    budget.provenance(definition.provenance, "provenance");
    budget.provenance(definition.applicability.provenance, "applicability.provenance");
    for (std::size_t i = 0; i < n; ++i) {
        const auto& component = definition.components[i];
        const auto field = "components[" + std::to_string(i) + "]";
        budget.id(component.component_id, field + ".component_id");
        budget.label(component.display_name, field + ".display_name");
        budget.provenance(component.provenance, field + ".provenance");
        budget.scalar(component.molar_mass_kg_per_mol, field + ".molar_mass_kg_per_mol");
    }
    for (std::size_t i = 0; i < pr->pure.size(); ++i) {
        const auto& record = pr->pure[i];
        const auto field = "parameters.pure[" + std::to_string(i) + "]";
        budget.id(record.component_id, field + ".component_id");
        budget.scalar(record.critical_temperature_k, field + ".critical_temperature_k");
        budget.scalar(record.critical_pressure_pa, field + ".critical_pressure_pa");
        budget.scalar(record.acentric_factor, field + ".acentric_factor");
    }
    for (std::size_t i = 0; i < pr->binary.size(); ++i) {
        const auto& record = pr->binary[i];
        const auto field = "parameters.binary[" + std::to_string(i) + "]";
        budget.id(record.first_component_id, field + ".first_component_id");
        budget.id(record.second_component_id, field + ".second_component_id");
        budget.scalar(record.kij, field + ".kij");
    }
}

inline thermodynamics::SourceKind source_kind(SourceKind kind) {
    switch (kind) {
    case SourceKind::literature: return thermodynamics::SourceKind::literature;
    case SourceKind::database: return thermodynamics::SourceKind::database;
    case SourceKind::user_supplied: return thermodynamics::SourceKind::user_supplied;
    case SourceKind::assumption: return thermodynamics::SourceKind::assumption;
    case SourceKind::synthetic_test: return thermodynamics::SourceKind::synthetic_test;
    default: return thermodynamics::SourceKind::unspecified;
    }
}
inline thermodynamics::Provenance source(const ModelProvenance& value) {
    return {source_kind(value.kind), value.reference, value.revision, value.locator,
            value.note, value.acquisition, value.usage_terms};
}
inline std::optional<thermodynamics::SourcedScalar> scalar(
    const std::optional<ModelScalar>& value, thermodynamics::Unit unit) {
    if (!value) { return std::nullopt; }
    return thermodynamics::SourcedScalar{
        value->value, unit, source(value->provenance), value->original_unit, value->conversion};
}

inline void validate_applicability_endpoint(
    const std::optional<double>& value, bool exclusive,
    const std::string& value_field, const std::string& exclusive_field) {
    model_require(!exclusive || value.has_value(),
                  ModelConfigurationErrorCode::invalid_range, exclusive_field,
                  "exclusive marker requires a declared endpoint");
    if (!value) { return; }
    model_require(std::isfinite(*value) && *value > 0.0,
                  ModelConfigurationErrorCode::invalid_range, value_field,
                  "applicability endpoint must be finite and positive");
}

inline void validate_applicability_axis(
    const std::optional<double>& lower, bool lower_exclusive,
    const std::optional<double>& upper, bool upper_exclusive,
    const std::string& lower_field, const std::string& lower_exclusive_field,
    const std::string& upper_field, const std::string& upper_exclusive_field,
    const std::string& axis_field) {
    validate_applicability_endpoint(lower, lower_exclusive, lower_field,
                                    lower_exclusive_field);
    validate_applicability_endpoint(upper, upper_exclusive, upper_field,
                                    upper_exclusive_field);
    if (!lower || !upper) { return; }
    model_require(*lower <= *upper, ModelConfigurationErrorCode::invalid_range, axis_field,
                  "applicability lower endpoint exceeds upper endpoint");
    model_require(*lower < *upper || (!lower_exclusive && !upper_exclusive),
                  ModelConfigurationErrorCode::invalid_range, axis_field,
                  "exclusive equal endpoints declare an empty applicability interval");
}

inline void validate_applicability(const ModelApplicability& bounds) {
    validate_applicability_axis(
        bounds.temperature_lower_k, bounds.temperature_lower_exclusive,
        bounds.temperature_upper_k, bounds.temperature_upper_exclusive,
        "applicability.temperature_lower_k", "applicability.temperature_lower_exclusive",
        "applicability.temperature_upper_k", "applicability.temperature_upper_exclusive",
        "applicability.temperature_k");
    validate_applicability_axis(
        bounds.pressure_lower_pa, bounds.pressure_lower_exclusive,
        bounds.pressure_upper_pa, bounds.pressure_upper_exclusive,
        "applicability.pressure_lower_pa", "applicability.pressure_lower_exclusive",
        "applicability.pressure_upper_pa", "applicability.pressure_upper_exclusive",
        "applicability.pressure_pa");
}

// Native thermodynamics currently stores complete closed intervals. A complete
// public interval maps to that conservative numeric envelope; exclusive endpoints
// are enforced by the executable public guard. One-sided public bounds remain
// absent natively instead of fabricating their unknown opposite endpoint.
inline std::optional<thermodynamics::ClosedInterval> interval(
    std::optional<double> lower, std::optional<double> upper) {
    if (!lower || !upper) { return std::nullopt; }
    return thermodynamics::ClosedInterval{*lower, *upper};
}

inline ModelConfigurationErrorCode public_code(thermodynamics::ContractErrorCode code) {
    using C = thermodynamics::ContractErrorCode;
    using P = ModelConfigurationErrorCode;
    switch (code) {
    case C::missing_field: return P::missing_field;
    case C::invalid_identifier: return P::invalid_identifier;
    case C::duplicate_identifier: return P::duplicate_identifier;
    case C::unknown_component: return P::unknown_component;
    case C::invalid_value: return P::invalid_value;
    case C::invalid_unit: return P::invalid_unit;
    case C::invalid_source: return P::invalid_source;
    case C::invalid_range: return P::invalid_range;
    case C::unsupported_model: return P::unsupported_family;
    case C::duplicate_parameter: return P::duplicate_parameter;
    case C::missing_parameter: return P::missing_parameter;
    case C::invalid_pair: return P::invalid_pair;
    case C::size_limit: return P::resource_limit;
    }
    throw std::logic_error("unmapped thermodynamics contract error");
}

inline std::string public_field(std::string field) {
    if (field.starts_with("pure[") || field.starts_with("binary[")) {
        field.insert(0, "parameters.");
    } else if (field == "binary") {
        field = "parameters.binary";
    } else if (field.starts_with("catalog[") || field.starts_with("component[")) {
        field.replace(0, field.find('['), "components");
    } else if (field.starts_with("order")) {
        field.replace(0, 5, "components");
    }
    const std::size_t start = field.rfind(']');
    const std::size_t suffix = start == std::string::npos ? 0 : start + 1;
    const auto replace = [&](std::string_view old_name, std::string_view new_name) {
        const auto pos = field.find(old_name, suffix);
        if (pos != std::string::npos) { field.replace(pos, old_name.size(), new_name); }
    };
    replace(".critical_temperature", ".critical_temperature_k");
    replace(".critical_pressure", ".critical_pressure_pa");
    replace(".molar_mass", ".molar_mass_kg_per_mol");
    replace(".definition", ".provenance");
    replace(".declaration", ".provenance");
    replace(".source", ".provenance");
    if (field == "component.id") { field = "components.component_id"; }
    return field;
}
} // namespace detail

class Pr76ModelParameters {
public:
    Pr76ModelParameters(const Pr76ModelParameters&) = default;
    Pr76ModelParameters(Pr76ModelParameters&&) noexcept = default;
    Pr76ModelParameters& operator=(const Pr76ModelParameters&) = delete;
    Pr76ModelParameters& operator=(Pr76ModelParameters&&) = delete;

    [[nodiscard]] static Pr76ModelParameters create(
        const ThermodynamicModelDefinition& definition,
        ModelConfigurationLimits limits = {},
        ModelDataPolicy policy = ModelDataPolicy::ordinary) {
        detail::bound_model_definition(definition, limits);
        detail::model_require(policy == ModelDataPolicy::ordinary ||
                                  policy == ModelDataPolicy::allow_synthetic_tests,
                              ModelConfigurationErrorCode::invalid_source, "data_policy",
                              "unknown host data policy");
        const auto data_policy = policy == ModelDataPolicy::ordinary
            ? thermodynamics::DataPolicy::ordinary : thermodynamics::DataPolicy::allow_synthetic_tests;
        try {
            thermodynamics::detail::require_text(definition.display_name, "display_name");
            thermodynamics::detail::validate_source(
                detail::source(definition.provenance), "provenance", data_policy);
            detail::validate_applicability(definition.applicability);
            thermodynamics::PrParameterInput input;
            input.model_id = thermodynamics::pr76_profile;
            input.dataset_id = definition.dataset_id;
            input.revision = definition.revision;
            const auto& bounds = definition.applicability;
            input.applicability = {
                detail::interval(bounds.temperature_lower_k, bounds.temperature_upper_k),
                detail::interval(bounds.pressure_lower_pa, bounds.pressure_upper_pa),
                detail::source(bounds.provenance)};
            std::vector<thermodynamics::Component> catalog;
            std::vector<std::string> order;
            catalog.reserve(definition.components.size());
            order.reserve(definition.components.size());
            for (const auto& component : definition.components) {
                auto kind = thermodynamics::ComponentKind::unspecified;
                if (component.kind == ComponentKind::pure) { kind = thermodynamics::ComponentKind::pure; }
                if (component.kind == ComponentKind::pseudo) { kind = thermodynamics::ComponentKind::pseudo; }
                catalog.push_back({component.component_id, component.display_name, kind,
                                   detail::source(component.provenance),
                                   detail::scalar(component.molar_mass_kg_per_mol,
                                                  thermodynamics::Unit::kilogram_per_mole)});
                order.push_back(component.component_id);
            }
            const auto& pr = std::get<Pr76ParameterDefinition>(definition.parameters);
            input.pure.reserve(pr.pure.size());
            input.binary.reserve(pr.binary.size());
            for (const auto& record : pr.pure) {
                input.pure.push_back({record.component_id,
                    detail::scalar(record.critical_temperature_k, thermodynamics::Unit::kelvin),
                    detail::scalar(record.critical_pressure_pa, thermodynamics::Unit::pascal),
                    detail::scalar(record.acentric_factor, thermodynamics::Unit::dimensionless)});
            }
            for (const auto& record : pr.binary) {
                input.binary.push_back({record.first_component_id, record.second_component_id,
                    detail::scalar(record.kij, thermodynamics::Unit::dimensionless)});
            }
            auto parameters = thermodynamics::PrParameterSet::create(catalog, order, input,
                data_policy, {limits.max_components, limits.max_pair_records, limits.max_matrix_entries});
            return Pr76ModelParameters(definition, std::move(parameters));
        } catch (const thermodynamics::ContractError& error) {
            throw ModelConfigurationError(detail::public_code(error.code()),
                                          detail::public_field(error.field()), error.what());
        }
    }

    [[nodiscard]] const ThermodynamicModelDefinition& definition() const & noexcept { return definition_; }
    const ThermodynamicModelDefinition& definition() const && = delete;
    [[nodiscard]] const thermodynamics::PrParameterSet& parameters() const & noexcept { return parameters_; }
    const thermodynamics::PrParameterSet& parameters() const && = delete;
private:
    Pr76ModelParameters(ThermodynamicModelDefinition definition,
                        thermodynamics::PrParameterSet parameters)
        : definition_(std::move(definition)), parameters_(std::move(parameters)) {}
    ThermodynamicModelDefinition definition_;
    thermodynamics::PrParameterSet parameters_;
};

} // namespace mpmc::model_configuration
#endif // MPMC_MODEL_CONFIGURATION_PR76_PARAMETERS_HPP
