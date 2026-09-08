#ifndef MPMC_THERMODYNAMICS_COMPONENTS_HPP
#define MPMC_THERMODYNAMICS_COMPONENTS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::thermodynamics {

/// Stable diagnostics for input contracts; allocation exceptions propagate unchanged.
enum class ContractErrorCode {
    missing_field, invalid_identifier, duplicate_identifier, unknown_component,
    invalid_value, invalid_unit, invalid_source, invalid_range, unsupported_model,
    duplicate_parameter, missing_parameter, invalid_pair, size_limit
};

class ContractError : public std::invalid_argument {
public:
    ContractError(ContractErrorCode code, std::string field, std::string_view reason)
        : std::invalid_argument(field + ": " + std::string(reason)), code_(code),
          field_(std::move(field)) {}
    [[nodiscard]] ContractErrorCode code() const noexcept { return code_; }
    [[nodiscard]] const std::string& field() const & noexcept { return field_; }
    const std::string& field() const && = delete;
private:
    ContractErrorCode code_;
    std::string field_;
};

/// These are input records, NOT validated scientific data. Builders validate them.
enum class SourceKind {
    unspecified, literature, database, user_supplied, assumption, synthetic_test
};
enum class DataPolicy { ordinary, allow_synthetic_tests };
enum class Unit { unspecified, kelvin, pascal, dimensionless, kilogram_per_mole };
enum class ComponentKind { unspecified, pure, pseudo };

struct Provenance {
    SourceKind kind = SourceKind::unspecified;
    std::string reference; // DOI, database/file identity, or explicit user/test record identity.
    std::string revision;  // Publication/data version or content digest; never inferred.
    std::string locator;   // Table/row/section/record key locating this particular datum.
    std::string note;      // Required for assumptions and synthetic data; explains their purpose.
    std::string acquisition; // How the record was obtained (not a claim of verified authenticity).
    std::string usage_terms; // License/permission statement; redistribution needs separate review.
};

struct SourcedScalar {
    double value = std::numeric_limits<double>::quiet_NaN(); // Missing is never silently zero.
    Unit unit = Unit::unspecified; // Stored value's unit, checked against the field's contract.
    Provenance source;
    std::string original_unit;
    std::string conversion; // Explicit conversion record, e.g. "identity" or "MPa * 1e6 -> Pa".
};

struct Component {
    std::string id;           // Stable, case-sensitive ASCII key, NOT a display name or position.
    std::string display_name; // Non-unique label; UTF-8 allowed. No chemical identity inference.
    ComponentKind kind = ComponentKind::unspecified;
    Provenance definition;   // Pseudocomponents must identify their characterization definition.
    std::optional<SourcedScalar> molar_mass; // kg/mol; not required by the PR76 parameter contract.
};

/// Logical limits, not a memory quota. Services must also bound bytes and input parsing.
struct ContractLimits {
    std::size_t max_components = std::numeric_limits<std::size_t>::max();
    std::size_t max_pair_records = std::numeric_limits<std::size_t>::max();
    std::size_t max_matrix_entries = std::numeric_limits<std::size_t>::max();
};

namespace detail {
inline void require(bool condition, ContractErrorCode code, const std::string& field,
                    std::string_view reason) {
    if (!condition) {
        throw ContractError(code, field, reason);
    }
}

inline void require_text(const std::string& text, const std::string& field) {
    const bool nonblank = std::any_of(text.begin(), text.end(),
                                      [](unsigned char c) { return c > 32; });
    require(nonblank && text.find('\0') == std::string::npos,
            ContractErrorCode::missing_field, field, "nonblank text without NUL is required");
}

inline void validate_id(const std::string& id, const std::string& field) {
    require(!id.empty() && std::all_of(id.begin(), id.end(), [](unsigned char c) {
        return c >= 33 && c <= 126;
    }), ContractErrorCode::invalid_identifier, field, "expected a nonempty visible ASCII key");
}

inline void validate_source(const Provenance& source, const std::string& field, DataPolicy policy) {
    require(policy == DataPolicy::ordinary || policy == DataPolicy::allow_synthetic_tests,
            ContractErrorCode::invalid_source, field, "unknown data policy");
    switch (source.kind) {
    case SourceKind::literature:
    case SourceKind::database:
    case SourceKind::user_supplied:
        break;
    case SourceKind::assumption:
        require_text(source.note, field + ".note");
        break;
    case SourceKind::synthetic_test:
        require(policy == DataPolicy::allow_synthetic_tests, ContractErrorCode::invalid_source,
                field, "synthetic test data require explicit opt-in");
        require_text(source.note, field + ".note");
        break;
    default:
        throw ContractError(ContractErrorCode::invalid_source, field, "unknown source kind");
    }
    require_text(source.reference, field + ".reference");
    require_text(source.revision, field + ".revision");
    require_text(source.locator, field + ".locator");
    require_text(source.acquisition, field + ".acquisition");
    require_text(source.usage_terms, field + ".usage_terms");
}

inline void validate_scalar(const SourcedScalar& datum, Unit unit, bool positive,
                            const std::string& field, DataPolicy policy) {
    require(datum.unit == unit, ContractErrorCode::invalid_unit, field, "wrong stored unit");
    require(std::isfinite(datum.value) && (!positive || datum.value > 0),
            ContractErrorCode::invalid_value, field, "expected finite value in the field's domain");
    validate_source(datum.source, field + ".source", policy);
    require_text(datum.original_unit, field + ".original_unit");
    require_text(datum.conversion, field + ".conversion");
    // Metadata is checked for completeness, not truth; no unit conversion is performed here.
}

inline void validate_component(const Component& component, DataPolicy policy) {
    validate_id(component.id, "component.id");
    const std::string field = "component[" + component.id + "]";
    require_text(component.display_name, field + ".display_name");
    require(component.kind == ComponentKind::pure || component.kind == ComponentKind::pseudo,
            ContractErrorCode::invalid_value, field + ".kind", "unknown component kind");
    validate_source(component.definition, field + ".definition", policy);
    if (component.molar_mass) {
        validate_scalar(*component.molar_mass, Unit::kilogram_per_mole, true,
                        field + ".molar_mass", policy);
    }
}
} // namespace detail

/// An owning snapshot of a requested order. No mutating element access or assignment.
/// Rebuild after adding, removing, replacing, or reordering IDs; old snapshots survive.
/// Different labels may share text, but IDs in a catalog/order must be unique.
class OrderedComponents {
public:
    OrderedComponents(const OrderedComponents&) = default;
    OrderedComponents(OrderedComponents&&) noexcept = default;
    OrderedComponents& operator=(const OrderedComponents&) = delete;
    OrderedComponents& operator=(OrderedComponents&&) = delete;

    [[nodiscard]] static OrderedComponents select(
        std::span<const Component> catalog, std::span<const std::string> order,
        DataPolicy policy = DataPolicy::ordinary, ContractLimits limits = {}) {
        detail::require(!order.empty(), ContractErrorCode::missing_field, "order", "empty order");
        detail::require(catalog.size() <= limits.max_components &&
                            order.size() <= limits.max_components,
                        ContractErrorCode::size_limit, "components", "component limit exceeded");
        std::map<std::string, const Component*, std::less<>> available;
        for (const auto& component : catalog) {
            detail::validate_component(component, policy);
            detail::require(available.emplace(component.id, &component).second,
                            ContractErrorCode::duplicate_identifier,
                            "catalog[" + component.id + "]",
                            "duplicate component ID");
        }
        OrderedComponents result;
        result.items_.reserve(order.size());
        for (const auto& id : order) {
            detail::validate_id(id, "order.id");
            const auto found = available.find(id);
            detail::require(found != available.end(), ContractErrorCode::unknown_component,
                            "order[" + id + "]", "ID not in catalog");
            detail::require(result.indices_.emplace(id, result.items_.size()).second,
                            ContractErrorCode::duplicate_identifier, "order[" + id + "]",
                            "duplicate selection");
            result.items_.push_back(*found->second);
        }
        return result;
    }

    [[nodiscard]] std::size_t size() const noexcept { return items_.size(); }
    [[nodiscard]] std::span<const Component> items() const & noexcept { return items_; }
    std::span<const Component> items() const && = delete;
    [[nodiscard]] const Component& at(std::size_t index) const & { return items_.at(index); }
    const Component& at(std::size_t) const && = delete;
    [[nodiscard]] bool contains(std::string_view id) const { return indices_.contains(id); }
    [[nodiscard]] std::size_t index_of(std::string_view id) const {
        const auto found = indices_.find(id);
        detail::require(found != indices_.end(), ContractErrorCode::unknown_component,
                        "order[" + std::string(id) + "]", "ID not in this snapshot");
        return found->second;
    }
private:
    OrderedComponents() = default;
    std::vector<Component> items_;
    std::map<std::string, std::size_t, std::less<>> indices_;
};

/// Optional bounds are UNKNOWN, never shorthand for an unbounded validated range.
struct ClosedInterval { double lower; double upper; };
enum class RangeAssessment { inside_declared_bounds, outside_declared_bounds, unknown };
struct Applicability {
    std::optional<ClosedInterval> temperature_k;
    std::optional<ClosedInterval> pressure_pa;
    Provenance declaration; // Scope/evidence for the WHOLE dataset, including omissions.

    void validate(DataPolicy policy = DataPolicy::ordinary) const {
        detail::validate_source(declaration, "applicability.declaration", policy);
        validate_bounds();
    }
    /// Only tests declared T/p bounds; this is NOT proof of physical model accuracy.
    [[nodiscard]] RangeAssessment assess(double temperature, double pressure) const {
        validate_bounds();
        detail::require(std::isfinite(temperature) && temperature > 0 &&
                            std::isfinite(pressure) && pressure > 0,
                        ContractErrorCode::invalid_value, "state", "expected finite positive K/Pa");
        const auto outside = [](const auto& bounds, double value) {
            return bounds && (value < bounds->lower || value > bounds->upper);
        };
        if (outside(temperature_k, temperature) || outside(pressure_pa, pressure)) {
            return RangeAssessment::outside_declared_bounds;
        }
        return temperature_k && pressure_pa ? RangeAssessment::inside_declared_bounds
                                            : RangeAssessment::unknown;
    }
private:
    void validate_bounds() const {
        for (const auto* interval : {&temperature_k, &pressure_pa}) {
            if (*interval) {
                const auto& bounds = **interval;
                detail::require(std::isfinite(bounds.lower) && std::isfinite(bounds.upper) &&
                                    bounds.lower > 0 && bounds.lower <= bounds.upper,
                                ContractErrorCode::invalid_range, "applicability",
                                "bounds must be finite, positive and ordered in K/Pa");
            }
        }
    }
};

} // namespace mpmc::thermodynamics
#endif // MPMC_THERMODYNAMICS_COMPONENTS_HPP
