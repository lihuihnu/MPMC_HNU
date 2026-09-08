#ifndef MPMC_THERMODYNAMICS_PR_PARAMETERS_HPP
#define MPMC_THERMODYNAMICS_PR_PARAMETERS_HPP

#include <mpmc/thermodynamics/components.hpp>

#include <cstddef>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::thermodynamics {

/// This profile fixes PR76 alpha, classical quadratic a / linear b mixing,
/// constant symmetric kij and NO volume translation. It is not SW, CPA or PR78.
/// Only the data contract is implemented; no EOS coefficients/roots are calculated.
inline constexpr std::string_view pr76_profile = "PR76/classical-vdw/constant-kij/no-translation";

struct PrPureRecord {
    std::string component_id;
    std::optional<SourcedScalar> critical_temperature; // K, >0; may be characterized.
    std::optional<SourcedScalar> critical_pressure;    // Pa, > 0.
    std::optional<SourcedScalar> acentric_factor;      // Finite, dimensionless; negatives allowed.
};
struct PrBinaryRecord {
    std::string first_id;
    std::string second_id;
    std::optional<SourcedScalar> kij; // Finite, dimensionless; no universal bound.
};
struct PrParameterInput {
    std::string model_id; // Mandatory exact profile match; no default or silent fallback.
    std::string dataset_id;
    std::string revision; // Dataset revision, distinct from model identity.
    Applicability applicability;
    std::vector<PrPureRecord> pure;
    std::vector<PrBinaryRecord> binary;
};

/// Complete, owning PR snapshot. Numeric arrays and records use the requested ID order.
/// Building a new snapshot never mutates the catalog, input records, or prior snapshots.
/// Common components/provenance can be reused by SW/CPA; PR-only pair rules cannot.
class PrParameterSet {
public:
    PrParameterSet(const PrParameterSet&) = default;
    PrParameterSet(PrParameterSet&&) noexcept = default;
    PrParameterSet& operator=(const PrParameterSet&) = delete;
    PrParameterSet& operator=(PrParameterSet&&) = delete;

    [[nodiscard]] static PrParameterSet create(
        std::span<const Component> catalog, std::span<const std::string> order,
        const PrParameterInput& input, DataPolicy policy = DataPolicy::ordinary,
        ContractLimits limits = {}) {
        detail::require(input.model_id == pr76_profile, ContractErrorCode::unsupported_model,
                        "model_id", "requires the exact PR76 profile; no model substitution");
        detail::require_text(input.dataset_id, "dataset_id");
        detail::require_text(input.revision, "revision");
        input.applicability.validate(policy);
        detail::require(!order.empty(), ContractErrorCode::missing_field, "order", "empty order");
        const std::size_t n = order.size();
        detail::require(n <= limits.max_components && catalog.size() <= limits.max_components &&
                            input.binary.size() <= limits.max_pair_records &&
                            n <= limits.max_matrix_entries / n &&
                            n <= std::vector<double>{}.max_size() / n,
                        ContractErrorCode::size_limit, "parameters", "shape or quota exceeded");
        auto components = OrderedComponents::select(catalog, order, policy, limits);
        std::set<std::string, std::less<>> known_ids;
        for (const auto& component : catalog) {
            known_ids.insert(component.id);
        }
        std::map<std::string, const PrPureRecord*, std::less<>> pure;
        for (const auto& record : input.pure) {
            const std::string field = "pure[" + record.component_id + "]";
            detail::require(known_ids.contains(record.component_id),
                            ContractErrorCode::unknown_component, field, "ID not in catalog");
            detail::require(pure.emplace(record.component_id, &record).second,
                            ContractErrorCode::duplicate_parameter, field, "duplicate pure record");
            const auto check = [&](const auto& datum, Unit unit, bool positive, const char* name) {
                detail::require(datum.has_value(), ContractErrorCode::missing_parameter,
                                field + "." + name, "required parameter is absent");
                detail::validate_scalar(*datum, unit, positive, field + "." + name, policy);
            };
            check(record.critical_temperature, Unit::kelvin, true, "critical_temperature");
            check(record.critical_pressure, Unit::pascal, true, "critical_pressure");
            check(record.acentric_factor, Unit::dimensionless, false, "acentric_factor");
        }
        // Canonical unordered keys belong ONLY to this symmetric constant-kij profile.
        // SW phase-dependent correlations and CPA site pairs need different records.
        using Pair = std::pair<std::string, std::string>;
        std::map<Pair, const PrBinaryRecord*> pairs;
        std::size_t selected_pair_count = 0;
        for (const auto& record : input.binary) {
            const std::string field = "binary[" + record.first_id + "," + record.second_id + "]";
            detail::require(known_ids.contains(record.first_id) &&
                                known_ids.contains(record.second_id),
                            ContractErrorCode::unknown_component, field, "pair ID not in catalog");
            detail::require(record.first_id != record.second_id, ContractErrorCode::invalid_pair,
                            field, "self-pair is structural kij=0, not a fitted record");
            const Pair key = record.first_id < record.second_id
                                 ? Pair{record.first_id, record.second_id}
                                 : Pair{record.second_id, record.first_id};
            detail::require(pairs.emplace(key, &record).second,
                            ContractErrorCode::duplicate_parameter,
                            field, "duplicate unordered pair (including reversed pairs)");
            detail::require(record.kij.has_value(), ContractErrorCode::missing_parameter,
                            field + ".kij", "required interaction parameter is absent");
            detail::validate_scalar(*record.kij, Unit::dimensionless, false,
                                    field + ".kij", policy);
            if (components.contains(record.first_id) && components.contains(record.second_id)) {
                ++selected_pair_count;
            }
        }

        // Uniqueness plus this count proves selected-pair completeness BEFORE dense allocation.
        detail::require(selected_pair_count == n * (n - 1) / 2,
                        ContractErrorCode::missing_parameter, "binary",
                        "selected component pairs are incomplete; no implicit zero is allowed");
        PrParameterSet result(std::move(components), input);
        result.pure_.reserve(n);
        result.tc_.reserve(n);
        result.pc_.reserve(n);
        result.omega_.reserve(n);
        for (const auto& id : order) {
            const auto found = pure.find(id);
            detail::require(found != pure.end(), ContractErrorCode::missing_parameter,
                            "pure[" + id + "]", "selected component lacks PR parameters");
            const auto& record = *found->second;
            result.pure_.push_back(record);
            result.tc_.push_back(record.critical_temperature->value);
            result.pc_.push_back(record.critical_pressure->value);
            result.omega_.push_back(record.acentric_factor->value);
        }
        result.kij_.assign(n * n, 0.0); // Only the diagonal is an implicit, structural zero.
        result.binary_.reserve(n * (n - 1) / 2); // n*n was checked before multiplication.
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                const Pair key = order[i] < order[j] ? Pair{order[i], order[j]}
                                                     : Pair{order[j], order[i]};
                const auto found = pairs.find(key);
                detail::require(found != pairs.end(), ContractErrorCode::missing_parameter,
                                "binary[" + order[i] + "," + order[j] + "]",
                                "no implicit zero or inherited interaction for a missing pair");
                PrBinaryRecord record = *found->second;
                record.first_id = order[i]; // Reorient metadata consistently with the snapshot.
                record.second_id = order[j];
                result.kij_[i * n + j] = result.kij_[j * n + i] = record.kij->value;
                result.binary_.push_back(std::move(record));
            }
        }
        return result; // A missing/invalid record can never publish a partial snapshot.
    }

    [[nodiscard]] static constexpr std::string_view model_id() noexcept { return pr76_profile; }
    [[nodiscard]] const std::string& dataset_id() const & noexcept { return dataset_id_; }
    const std::string& dataset_id() const && = delete;
    [[nodiscard]] const std::string& revision() const & noexcept { return revision_; }
    const std::string& revision() const && = delete;
    [[nodiscard]] const OrderedComponents& components() const & noexcept { return components_; }
    const OrderedComponents& components() const && = delete;
    [[nodiscard]] const Applicability& applicability() const & noexcept { return applicability_; }
    const Applicability& applicability() const && = delete;
    [[nodiscard]] std::span<const PrPureRecord> pure_records() const & noexcept { return pure_; }
    std::span<const PrPureRecord> pure_records() const && = delete;
    /// Records ordered by upper-triangle traversal i=0..n-1, j=i+1..n-1, with provenance.
    [[nodiscard]] std::span<const PrBinaryRecord> binary_records() const & noexcept {
        return binary_;
    }
    std::span<const PrBinaryRecord> binary_records() const && = delete;
    [[nodiscard]] std::span<const double> critical_temperatures_k() const & noexcept { return tc_; }
    std::span<const double> critical_temperatures_k() const && = delete;
    [[nodiscard]] std::span<const double> critical_pressures_pa() const & noexcept { return pc_; }
    std::span<const double> critical_pressures_pa() const && = delete;
    [[nodiscard]] std::span<const double> acentric_factors() const & noexcept { return omega_; }
    std::span<const double> acentric_factors() const && = delete;
    /// Contiguous row-major matrix, n*n entries. Values have not been clipped or rescaled.
    [[nodiscard]] std::span<const double> kij_matrix() const & noexcept { return kij_; }
    std::span<const double> kij_matrix() const && = delete;
    [[nodiscard]] double kij(std::size_t i, std::size_t j) const {
        if (i >= components_.size() || j >= components_.size()) {
            throw std::out_of_range("PrParameterSet::kij: index outside component order");
        }
        return kij_.at(i * components_.size() + j);
    }
private:
    PrParameterSet(OrderedComponents components, const PrParameterInput& input)
        : components_(std::move(components)), dataset_id_(input.dataset_id),
          revision_(input.revision),
          applicability_(input.applicability) {}
    OrderedComponents components_;
    std::string dataset_id_;
    std::string revision_;
    Applicability applicability_;
    std::vector<PrPureRecord> pure_;
    std::vector<PrBinaryRecord> binary_;
    std::vector<double> tc_, pc_, omega_, kij_;
};

} // namespace mpmc::thermodynamics
#endif // MPMC_THERMODYNAMICS_PR_PARAMETERS_HPP
