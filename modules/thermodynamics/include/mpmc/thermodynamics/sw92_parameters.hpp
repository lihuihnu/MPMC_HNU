#ifndef MPMC_THERMODYNAMICS_SW92_PARAMETERS_HPP
#define MPMC_THERMODYNAMICS_SW92_PARAMETERS_HPP

#include <mpmc/thermodynamics/sw92_types.hpp>

#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::thermodynamics {

/// Owning, ordered SW92 data snapshot. It is intentionally separate from
/// PrParameterSet because AQ/NA interactions have different semantics.
class Sw92ParameterSet {
public:
    Sw92ParameterSet(const Sw92ParameterSet&) = default;
    Sw92ParameterSet(Sw92ParameterSet&&) noexcept = default;
    Sw92ParameterSet& operator=(const Sw92ParameterSet&) = delete;
    Sw92ParameterSet& operator=(Sw92ParameterSet&&) = delete;

    [[nodiscard]] static Sw92ParameterSet create(
        std::span<const Component> catalog, std::span<const std::string> order,
        const Sw92ParameterInput& input, DataPolicy policy = DataPolicy::ordinary,
        ContractLimits limits = {}) {
        detail::require(input.model_id == sw92_corrected_profile,
                        ContractErrorCode::unsupported_model, "model_id",
                        "requires exact SW92/corrected-original profile");
        detail::require_text(input.dataset_id, "dataset_id");
        detail::require_text(input.revision, "revision");
        input.applicability.validate(policy);
        detail::require(!order.empty(), ContractErrorCode::missing_field,
                        "order", "empty order");
        const std::size_t n = order.size();
        detail::require(n <= limits.max_components &&
                            catalog.size() <= limits.max_components &&
                            input.water_binary.size() <= limits.max_pair_records &&
                            input.nonwater_binary.size() <=
                                limits.max_pair_records - input.water_binary.size() &&
                            n <= limits.max_matrix_entries / n &&
                            n <= std::vector<double>{}.max_size() / n,
                        ContractErrorCode::size_limit, "parameters",
                        "shape or quota exceeded");

        auto components = OrderedComponents::select(catalog, order, policy, limits);
        std::set<std::string, std::less<>> known;
        for (const auto& c : catalog) known.insert(c.id);

        std::map<std::string, const Sw92PureRecord*, std::less<>> pure;
        for (const auto& r : input.pure) validate_pure(r, known, pure, policy);

        Sw92ParameterSet out(std::move(components), input);
        out.pure_.reserve(n); out.species_.reserve(n);
        out.tc_.reserve(n); out.pc_.reserve(n); out.omega_.reserve(n);
        std::size_t water_count = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const auto f = pure.find(order[i]);
            detail::require(f != pure.end(), ContractErrorCode::missing_parameter,
                            "pure[" + order[i] + "]",
                            "selected component lacks SW92 pure parameters");
            const auto& r = *f->second;
            out.pure_.push_back(r); out.species_.push_back(r.species);
            out.tc_.push_back(r.critical_temperature->value);
            out.pc_.push_back(r.critical_pressure->value);
            out.omega_.push_back(r.acentric_factor->value);
            if (r.species == Sw92Species::water) {
                out.water_index_ = i; ++water_count;
            }
        }
        detail::require(water_count == 1, ContractErrorCode::invalid_value,
                        "water", "selected snapshot requires exactly one water component");

        std::map<std::string, const Sw92WaterBinaryRecord*, std::less<>> water;
        for (const auto& r : input.water_binary) {
            const std::string field = "water_binary[" + r.component_id + "]";
            detail::require(known.contains(r.component_id), ContractErrorCode::unknown_component,
                            field, "ID not in catalog");
            const auto f = pure.find(r.component_id);
            detail::require(f != pure.end(), ContractErrorCode::missing_parameter,
                            field, "water-pair component lacks pure record");
            detail::require(f->second->species != Sw92Species::water,
                            ContractErrorCode::invalid_pair, field,
                            "water self-pair is structural kij=0");
            detail::require(water.emplace(r.component_id, &r).second,
                            ContractErrorCode::duplicate_parameter, field,
                            "duplicate water/non-water pair");
            validate_water_rule(r, f->second->species, field, policy);
        }
        out.water_rules_.assign(n, Sw92NonAqueousWaterRule::unspecified);
        out.water_na_constant_.assign(n, std::numeric_limits<double>::quiet_NaN());
        for (std::size_t i = 0; i < n; ++i) {
            if (i == out.water_index_) continue;
            const auto f = water.find(order[i]);
            detail::require(f != water.end(), ContractErrorCode::missing_parameter,
                            "water_binary[" + order[i] + "]",
                            "selected non-water component lacks water-pair data");
            out.water_binary_.push_back(*f->second);
            out.water_rules_[i] = f->second->nonaqueous_rule;
            if (f->second->nonaqueous_kij)
                out.water_na_constant_[i] = f->second->nonaqueous_kij->value;
        }

        using Pair = std::pair<std::string, std::string>;
        std::map<Pair, const Sw92NonWaterBinaryRecord*> pairs;
        std::size_t selected_pairs = 0;
        for (const auto& r : input.nonwater_binary) {
            const std::string field = "nonwater_binary[" + r.first_id + "," + r.second_id + "]";
            detail::require(known.contains(r.first_id) && known.contains(r.second_id),
                            ContractErrorCode::unknown_component, field, "ID not in catalog");
            detail::require(r.first_id != r.second_id, ContractErrorCode::invalid_pair,
                            field, "self-pair is structural kij=0");
            const auto a = pure.find(r.first_id), b = pure.find(r.second_id);
            detail::require(a != pure.end() && b != pure.end(),
                            ContractErrorCode::missing_parameter, field,
                            "pair component lacks pure record");
            detail::require(a->second->species != Sw92Species::water &&
                                b->second->species != Sw92Species::water,
                            ContractErrorCode::invalid_pair, field,
                            "water pairs belong in water_binary");
            const Pair key = r.first_id < r.second_id ? Pair{r.first_id, r.second_id}
                                                      : Pair{r.second_id, r.first_id};
            detail::require(pairs.emplace(key, &r).second,
                            ContractErrorCode::duplicate_parameter, field,
                            "duplicate unordered pair");
            validate_kij(r.aqueous_kij, field + ".aqueous_kij", policy);
            validate_kij(r.nonaqueous_kij, field + ".nonaqueous_kij", policy);
            if (out.components_.contains(r.first_id) && out.components_.contains(r.second_id))
                ++selected_pairs;
        }
        const std::size_t nw = n - 1;
        const std::size_t required = nw < 2 ? 0 : nw * (nw - 1) / 2;
        detail::require(selected_pairs == required, ContractErrorCode::missing_parameter,
                        "nonwater_binary", "selected non-water pairs are incomplete");
        out.aq_nonwater_.assign(n * n, 0.0);
        out.na_nonwater_.assign(n * n, 0.0);
        out.nonwater_binary_.reserve(selected_pairs);
        for (std::size_t i = 0; i < n; ++i) {
            if (i == out.water_index_) continue;
            for (std::size_t j = i + 1; j < n; ++j) {
                if (j == out.water_index_) continue;
                const Pair key = order[i] < order[j] ? Pair{order[i], order[j]}
                                                     : Pair{order[j], order[i]};
                const auto f = pairs.find(key);
                detail::require(f != pairs.end(), ContractErrorCode::missing_parameter,
                                "nonwater_binary[" + order[i] + "," + order[j] + "]",
                                "missing selected non-water pair");
                auto r = *f->second; r.first_id = order[i]; r.second_id = order[j];
                const double aq = r.aqueous_kij->value, na = r.nonaqueous_kij->value;
                out.aq_nonwater_[i*n+j] = out.aq_nonwater_[j*n+i] = aq;
                out.na_nonwater_[i*n+j] = out.na_nonwater_[j*n+i] = na;
                out.nonwater_binary_.push_back(std::move(r));
            }
        }
        return out;
    }

    [[nodiscard]] static constexpr std::string_view model_id() noexcept { return sw92_corrected_profile; }
    [[nodiscard]] const std::string& dataset_id() const & noexcept { return dataset_id_; }
    const std::string& dataset_id() const && = delete;
    [[nodiscard]] const std::string& revision() const & noexcept { return revision_; }
    const std::string& revision() const && = delete;
    [[nodiscard]] const OrderedComponents& components() const & noexcept { return components_; }
    const OrderedComponents& components() const && = delete;
    [[nodiscard]] const Sw92Applicability& applicability() const & noexcept { return applicability_; }
    const Sw92Applicability& applicability() const && = delete;
    [[nodiscard]] std::size_t water_index() const noexcept { return water_index_; }
    [[nodiscard]] Sw92Species species(std::size_t i) const { return species_.at(i); }
    [[nodiscard]] std::span<const Sw92PureRecord> pure_records() const & noexcept { return pure_; }
    std::span<const Sw92PureRecord> pure_records() const && = delete;
    [[nodiscard]] std::span<const double> critical_temperatures_k() const & noexcept { return tc_; }
    [[nodiscard]] std::span<const double> critical_pressures_pa() const & noexcept { return pc_; }
    [[nodiscard]] std::span<const double> acentric_factors() const & noexcept { return omega_; }
    [[nodiscard]] std::span<const Sw92WaterBinaryRecord> water_binary_records() const & noexcept { return water_binary_; }
    [[nodiscard]] std::span<const Sw92NonWaterBinaryRecord> nonwater_binary_records() const & noexcept { return nonwater_binary_; }

    [[nodiscard]] Sw92NonAqueousWaterRule water_nonaqueous_rule(std::size_t i) const {
        check_nonwater(i); return water_rules_.at(i);
    }
    [[nodiscard]] double water_nonaqueous_constant_kij(std::size_t i) const {
        check_nonwater(i);
        if (water_rules_.at(i) != Sw92NonAqueousWaterRule::constant)
            throw std::invalid_argument("Sw92ParameterSet: water pair is correlated, not constant");
        return water_na_constant_.at(i);
    }
    [[nodiscard]] double nonwater_kij(SwPhaseFamily family, std::size_t i, std::size_t j) const {
        const std::size_t n = components_.size();
        if (i >= n || j >= n) throw std::out_of_range("Sw92ParameterSet::nonwater_kij: index");
        if (i == water_index_ || j == water_index_)
            throw std::invalid_argument("Sw92ParameterSet::nonwater_kij: water pair");
        if (i == j) return 0.0;
        switch (family) {
        case SwPhaseFamily::aqueous: return aq_nonwater_.at(i*n+j);
        case SwPhaseFamily::nonaqueous: return na_nonwater_.at(i*n+j);
        }
        throw std::invalid_argument("Sw92ParameterSet: unknown phase family");
    }

private:
    Sw92ParameterSet(OrderedComponents c, const Sw92ParameterInput& in)
        : components_(std::move(c)), dataset_id_(in.dataset_id), revision_(in.revision),
          applicability_(in.applicability) {}

    static bool valid_species(Sw92Species s) noexcept {
        switch (s) {
        case Sw92Species::water: case Sw92Species::hydrocarbon:
        case Sw92Species::nitrogen: case Sw92Species::carbon_dioxide:
        case Sw92Species::hydrogen_sulfide: return true;
        case Sw92Species::unspecified: return false;
        }
        return false;
    }
    static void validate_kij(const std::optional<SourcedScalar>& v,
                             const std::string& f, DataPolicy p) {
        detail::require(v.has_value(), ContractErrorCode::missing_parameter, f,
                        "required interaction parameter is absent");
        detail::validate_scalar(*v, Unit::dimensionless, false, f, p);
    }
    static void validate_pure(const Sw92PureRecord& r,
                              const std::set<std::string, std::less<>>& known,
                              std::map<std::string, const Sw92PureRecord*, std::less<>>& out,
                              DataPolicy p) {
        const std::string f = "pure[" + r.component_id + "]";
        detail::require(known.contains(r.component_id), ContractErrorCode::unknown_component, f, "ID not in catalog");
        detail::require(out.emplace(r.component_id, &r).second, ContractErrorCode::duplicate_parameter, f, "duplicate pure record");
        detail::require(valid_species(r.species), ContractErrorCode::invalid_value, f+".species", "unknown SW92 species routing");
        const auto check = [&](const auto& v, Unit u, bool positive, const char* name) {
            detail::require(v.has_value(), ContractErrorCode::missing_parameter, f+"."+name, "required parameter is absent");
            detail::validate_scalar(*v, u, positive, f+"."+name, p);
        };
        check(r.critical_temperature, Unit::kelvin, true, "critical_temperature");
        check(r.critical_pressure, Unit::pascal, true, "critical_pressure");
        check(r.acentric_factor, Unit::dimensionless, false, "acentric_factor");
        if (r.species == Sw92Species::hydrocarbon)
            detail::require(r.acentric_factor->value > 0.0, ContractErrorCode::invalid_value,
                            f+".acentric_factor", "SW92 Eq.(12) requires omega > 0");
    }
    static void validate_water_rule(const Sw92WaterBinaryRecord& r, Sw92Species s,
                                    const std::string& f, DataPolicy p) {
        if (s == Sw92Species::hydrogen_sulfide) {
            detail::require(r.nonaqueous_rule == Sw92NonAqueousWaterRule::hydrogen_sulfide_eq17,
                            ContractErrorCode::invalid_value, f+".nonaqueous_rule", "H2S/water must use Eq.(17)");
            detail::require(!r.nonaqueous_kij, ContractErrorCode::invalid_value,
                            f+".nonaqueous_kij", "Eq.(17) pair must not also provide a constant");
            return;
        }
        detail::require(s != Sw92Species::water && s != Sw92Species::unspecified,
                        ContractErrorCode::invalid_pair, f, "invalid water-pair species");
        detail::require(r.nonaqueous_rule == Sw92NonAqueousWaterRule::constant,
                        ContractErrorCode::invalid_value, f+".nonaqueous_rule", "explicit NA constant required");
        validate_kij(r.nonaqueous_kij, f+".nonaqueous_kij", p);
    }
    void check_nonwater(std::size_t i) const {
        if (i >= components_.size()) throw std::out_of_range("Sw92ParameterSet: component index");
        if (i == water_index_) throw std::invalid_argument("Sw92ParameterSet: water self-pair");
    }

    OrderedComponents components_;
    std::string dataset_id_, revision_;
    Sw92Applicability applicability_;
    std::vector<Sw92PureRecord> pure_;
    std::vector<Sw92Species> species_;
    std::vector<Sw92WaterBinaryRecord> water_binary_;
    std::vector<Sw92NonWaterBinaryRecord> nonwater_binary_;
    std::vector<double> tc_, pc_, omega_, water_na_constant_, aq_nonwater_, na_nonwater_;
    std::vector<Sw92NonAqueousWaterRule> water_rules_;
    std::size_t water_index_{};
};

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_SW92_PARAMETERS_HPP
