#ifndef MPMC_THERMODYNAMICS_SW92_TYPES_HPP
#define MPMC_THERMODYNAMICS_SW92_TYPES_HPP

#include <mpmc/thermodynamics/components.hpp>

#include <cmath>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::thermodynamics {

/// Soreide-Whitson (1992), using the authors' corrections to Eq.(13) and
/// Table 2. NaCl is an external molality state, not an explicit EOS component.
inline constexpr std::string_view sw92_corrected_profile =
    "SW92/corrected-original/PR76-base/NaCl-molality";

enum class SwPhaseFamily { aqueous, nonaqueous };

/// Correlation routing is explicit data; component names are never parsed.
enum class Sw92Species {
    unspecified,
    water,
    hydrocarbon,
    nitrogen,
    carbon_dioxide,
    hydrogen_sulfide
};

enum class Sw92NonAqueousWaterRule {
    unspecified,
    constant,
    hydrogen_sulfide_eq17
};

struct Sw92Applicability {
    Applicability state;
    // Unit: mol NaCl / kg H2O. Missing means unknown, not unbounded validation.
    std::optional<ClosedInterval> nacl_molality_mol_per_kg_water;

    void validate(DataPolicy policy = DataPolicy::ordinary) const {
        state.validate(policy);
        if (!nacl_molality_mol_per_kg_water) return;
        const auto& b = *nacl_molality_mol_per_kg_water;
        detail::require(std::isfinite(b.lower) && std::isfinite(b.upper) &&
                            b.lower >= 0.0 && b.lower <= b.upper,
                        ContractErrorCode::invalid_range,
                        "applicability.nacl_molality_mol_per_kg_water",
                        "bounds must be finite, nonnegative and ordered");
    }
};

struct Sw92PureRecord {
    std::string component_id;
    Sw92Species species = Sw92Species::unspecified;
    std::optional<SourcedScalar> critical_temperature; // K, >0.
    std::optional<SourcedScalar> critical_pressure;    // Pa, >0.
    std::optional<SourcedScalar> acentric_factor;      // finite, dimensionless.
};

/// The component_id is the non-water side. AQ behavior is fixed by species and
/// corrected SW92 equations. NA is Eq.(17) for H2S or an explicit sourced value.
struct Sw92WaterBinaryRecord {
    std::string component_id;
    Sw92NonAqueousWaterRule nonaqueous_rule =
        Sw92NonAqueousWaterRule::unspecified;
    std::optional<SourcedScalar> nonaqueous_kij;
};

/// SW92 has no universal non-water/non-water BIP correlation. Both family
/// values are explicit; missing values are never silently zero or shared.
struct Sw92NonWaterBinaryRecord {
    std::string first_id;
    std::string second_id;
    std::optional<SourcedScalar> aqueous_kij;
    std::optional<SourcedScalar> nonaqueous_kij;
};

struct Sw92ParameterInput {
    std::string model_id;
    std::string dataset_id;
    std::string revision;
    Sw92Applicability applicability;
    std::vector<Sw92PureRecord> pure;
    std::vector<Sw92WaterBinaryRecord> water_binary;
    std::vector<Sw92NonWaterBinaryRecord> nonwater_binary;
};

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_SW92_TYPES_HPP
