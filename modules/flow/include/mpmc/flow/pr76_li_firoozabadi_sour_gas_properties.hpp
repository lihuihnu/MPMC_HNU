#ifndef MPMC_FLOW_PR76_LI_FIROOZABADI_SOUR_GAS_PROPERTIES_HPP
#define MPMC_FLOW_PR76_LI_FIROOZABADI_SOUR_GAS_PROPERTIES_HPP

#include <mpmc/ad/math.hpp>
#include <mpmc/flow/pr76_selected_phase_property_closure.hpp>
#include <mpmc/thermodynamics/pr76_mixture.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow {

inline constexpr std::string_view
    pr76_li_firoozabadi_sour_gas_property_convention =
        "flow/pr76/li-firoozabadi-2012-sour-gas/NIST-lowT-LBC/v1";

namespace pr76_li_firoozabadi_sour_gas_property_detail {

inline constexpr std::string_view required_pr76_dataset_id =
    "Li-Firoozabadi-2012-acid-gas-PR76";

inline constexpr std::array<std::string_view, 6> component_ids{
    "carbon-dioxide",
    "nitrogen",
    "hydrogen-sulfide",
    "methane",
    "ethane",
    "propane"};

/// NIST Chemistry WebBook SRD 69 molar masses [kg/mol].
inline constexpr std::array<double, 6>
    nist_molar_mass_kg_per_mol{
        0.0440095,
        0.0280134,
        0.034081,
        0.0160425,
        0.0300690,
        0.0440956};

/// Critical molar volumes [m3/mol] used only by the already-adopted
/// Lohrenz-Bray-Clark dense-fluid viscosity correlation.
///
/// CO2: 0.0919 L/mol, NIST WebBook, Li & Kiran (1988).
/// N2: 1 / (11.18 mol/L), NIST WebBook critical density (Jacobsen et al., 1986).
/// H2S: 1 / (10.2 mol/L), NIST WebBook critical density.
/// CH4: 0.09860 L/mol, Ambrose & Tsonopoulos (1995), NIST WebBook.
/// C2H6: 0.147 L/mol, NIST WebBook average critical volume.
/// C3H8: 0.200 L/mol, Ambrose & Tsonopoulos (1995), NIST WebBook.
inline constexpr std::array<double, 6>
    critical_molar_volume_m3_per_mol{
        0.0919e-3,
        (1.0 / 11.18) * 1.0e-3,
        (1.0 / 10.2) * 1.0e-3,
        0.09860e-3,
        0.147e-3,
        0.200e-3};

/// Common low-temperature ideal-gas Cp grid [K].
///
/// Values are NIST-JANAF / NIST Chemistry WebBook recommended ideal-gas
/// heat capacities.  The common interval deliberately stops at 298.15 K;
/// this provider must not extrapolate the Li-Firoozabadi 178.8 K benchmark
/// outside its sourced low-temperature caloric data.
inline constexpr std::array<double, 3> nist_cp_temperature_k{
    100.0,
    200.0,
    298.15};

inline constexpr std::array<std::array<double, 3>, 6>
    nist_ideal_gas_cp_j_per_mol_k{{
        {29.208, 32.359, 37.129}, // CO2, NIST-JANAF C-095
        {29.10366858, 29.10731520, 29.12379083308665}, // N2, NIST Shomate 100-500 K
        {33.259, 33.380, 34.192}, // H2S, NIST-JANAF H-080
        {33.28, 33.51, 35.69},    // CH4, NIST recommended
        {35.70, 42.30, 52.49},    // C2H6, NIST recommended
        {41.30, 56.07, 73.60}     // C3H8, NIST recommended
    }};

template <typename Number>
[[nodiscard]] inline double primal_value(
    const Number& value) {
    if constexpr (std::floating_point<Number>) {
        return static_cast<double>(value);
    } else {
        return primal_value(value.value());
    }
}

template <typename Number>
inline void require_finite(
    const Number& value,
    const char* name,
    bool positive = false) {
    const double primal = primal_value(value);
    if (!std::isfinite(primal) ||
        (positive && !(primal > 0.0))) {
        throw std::range_error(
            std::string{
                "mpmc::flow::PR76 Li-Firoozabadi sour-gas property: "} +
            name +
            (positive
                 ? " must be finite and strictly positive"
                 : " must be finite"));
    }
}

template <typename Number>
inline void validate_composition(
    std::span<const Number> composition) {
    if (composition.size() != component_ids.size()) {
        throw std::invalid_argument(
            "mpmc::flow::PR76 Li-Firoozabadi sour-gas property requires exactly six components");
    }

    long double sum = 0.0L;
    for (const auto& fraction : composition) {
        const double primal = primal_value(fraction);
        if (!std::isfinite(primal) ||
            !(primal > 0.0) ||
            !(primal < 1.0)) {
            throw std::invalid_argument(
                "mpmc::flow::PR76 Li-Firoozabadi sour-gas property requires positive-support composition");
        }
        sum += static_cast<long double>(primal);
    }

    const long double tolerance =
        4096.0L *
        std::numeric_limits<double>::epsilon();
    if (std::abs(sum - 1.0L) > tolerance) {
        throw std::invalid_argument(
            "mpmc::flow::PR76 Li-Firoozabadi sour-gas composition must sum to one");
    }
}

/// Piecewise-linear integral from the explicit common caloric reference
/// h_i^ig(298.15 K)=0 down to T in [100, 298.15] K.
template <typename Number>
[[nodiscard]] inline Number integrated_low_temperature_nist_cp(
    const Number& temperature_k,
    std::size_t component) {
    const double primal = primal_value(temperature_k);
    const double lower = nist_cp_temperature_k.front();
    const double upper = nist_cp_temperature_k.back();

    if (!std::isfinite(primal) ||
        primal < lower ||
        primal > upper) {
        throw std::out_of_range(
            "mpmc::flow::Li-Firoozabadi sour-gas NIST ideal-gas enthalpy is validated only on [100,298.15] K");
    }
    if (component >= nist_ideal_gas_cp_j_per_mol_k.size()) {
        throw std::out_of_range(
            "mpmc::flow::Li-Firoozabadi sour-gas NIST Cp component index out of range");
    }

    Number magnitude{0.0};
    const auto& cp =
        nist_ideal_gas_cp_j_per_mol_k[component];

    for (std::size_t segment = 0U;
         segment + 1U < nist_cp_temperature_k.size();
         ++segment) {
        const double t0 = nist_cp_temperature_k[segment];
        const double t1 = nist_cp_temperature_k[segment + 1U];
        if (primal >= t1) {
            continue;
        }

        const double slope =
            (cp[segment + 1U] - cp[segment]) /
            (t1 - t0);

        if (primal <= t0) {
            const double width = t1 - t0;
            magnitude +=
                cp[segment] * width +
                0.5 * slope * width * width;
        } else {
            const Number delta =
                Number{t1} - temperature_k;
            magnitude +=
                cp[segment + 1U] * delta -
                0.5 * slope * delta * delta;
        }
    }

    const Number result = -magnitude;
    require_finite(
        result,
        "ideal-gas sensible molar enthalpy");
    return result;
}

template <typename Number>
[[nodiscard]] inline Number lbc_polynomial(
    const Number& reduced_density) {
    const Number square =
        reduced_density * reduced_density;
    const Number cube =
        square * reduced_density;
    const Number fourth =
        square * square;
    return Number{0.1023} +
        Number{0.023364} * reduced_density +
        Number{0.058533} * square -
        Number{0.040758} * cube +
        Number{0.0093324} * fourth;
}

} // namespace pr76_li_firoozabadi_sour_gas_property_detail

/// Source-complete transport/caloric provider for the repository's
/// Li-Firoozabadi (2012) six-component PR76 sour-gas benchmark.
///
/// This class introduces no new constitutive model.  It specializes the same
/// model chain already used by Pr76MethaneEthanePropanePropertyProvider:
/// - selected PR76 density/fugacity branch;
/// - Stiel-Thodos dilute viscosity;
/// - Herning-Zipperer dilute-mixture blending;
/// - Kay pseudo-critical mixing;
/// - Lohrenz-Bray-Clark dense-fluid correction;
/// - piecewise-linear integration of sourced ideal-gas Cp;
/// - standard PR76 residual-enthalpy departure expression.
///
/// The caloric reference is h_i^ig(298.15 K)=0 for every component.  This is a
/// nonreactive-flow reference convention, not a heat-of-formation model.
class Pr76LiFiroozabadiSourGasPropertyProvider {
public:
    static constexpr std::string_view convention =
        pr76_li_firoozabadi_sour_gas_property_convention;

    explicit Pr76LiFiroozabadiSourGasPropertyProvider(
        const thermodynamics::Pr76Phase<double>& model)
        : model_(&model),
          mixture_(
              thermodynamics::Pr76Mixture<double>::
                  from_parameters(model.parameters())) {
        using namespace
            pr76_li_firoozabadi_sour_gas_property_detail;

        const auto& parameters = model.parameters();
        if (parameters.dataset_id() !=
                required_pr76_dataset_id ||
            parameters.components().size() !=
                component_ids.size()) {
            throw std::invalid_argument(
                "mpmc::flow::PR76 Li-Firoozabadi sour-gas provider only supports the repository Li-Firoozabadi-2012 acid-gas dataset");
        }

        for (std::size_t component = 0U;
             component < component_ids.size();
             ++component) {
            const auto& record =
                parameters.components().at(component);
            if (record.id != component_ids[component]) {
                throw std::invalid_argument(
                    "mpmc::flow::PR76 Li-Firoozabadi sour-gas component order mismatch");
            }
            if (!record.molar_mass ||
                record.molar_mass->unit !=
                    thermodynamics::Unit::kilogram_per_mole) {
                throw std::invalid_argument(
                    "mpmc::flow::PR76 Li-Firoozabadi sour-gas provider requires sourced NIST molar masses");
            }

            const double expected =
                nist_molar_mass_kg_per_mol[component];
            const double actual =
                record.molar_mass->value;
            const double scale =
                std::max(1.0, std::abs(expected));
            if (!std::isfinite(actual) ||
                std::abs(actual - expected) >
                    64.0 *
                        std::numeric_limits<double>::epsilon() *
                        scale) {
                throw std::invalid_argument(
                    "mpmc::flow::PR76 Li-Firoozabadi sour-gas molar-mass snapshot does not match the NIST SRD 69 provider record");
            }
        }

        std::copy(
            parameters.critical_temperatures_k().begin(),
            parameters.critical_temperatures_k().end(),
            critical_temperature_k_.begin());
        std::copy(
            parameters.critical_pressures_pa().begin(),
            parameters.critical_pressures_pa().end(),
            critical_pressure_pa_.begin());
    }

    [[nodiscard]] static SelectedPhasePropertyProvenance
    provenance() {
        return {
            {
                "PR76 selected molar density x NIST SRD 69 molar mass",
                "Li-Firoozabadi-2012-acid-gas-PR76__NIST-SRD69-MW",
                "v1"},
            {
                "Lohrenz-Bray-Clark + Stiel-Thodos + Herning-Zipperer",
                "LBC-1964__StielThodos-1961__HerningZipperer-1936__NIST-SRD69-critical-volume",
                "six-component-low-temperature-v1"},
            {
                "NIST-JANAF + NIST Chemistry WebBook low-temperature ideal-gas Cp + PR76 residual enthalpy",
                "NIST-JANAF-C095-H080__NIST-WebBook-N2-CH4-C2H6-C3H8-lowT-Cp__PR76-departure",
                "hig-298.15K-zero__100-298.15K__v1"},
            {
                "thermodynamic identity u=h-p/rho_mass",
                "PR76-NIST-SRD69-LBC",
                "v1"}};
    }

    template <typename Number>
    [[nodiscard]] Number ideal_gas_molar_enthalpy_j_per_mol(
        const Number& temperature_k,
        std::span<const Number> composition) const {
        using namespace
            pr76_li_firoozabadi_sour_gas_property_detail;
        validate_composition(composition);

        Number result{0.0};
        for (std::size_t component = 0U;
             component < component_ids.size();
             ++component) {
            result +=
                composition[component] *
                integrated_low_temperature_nist_cp(
                    temperature_k,
                    component);
        }
        require_finite(
            result,
            "mixture ideal-gas molar enthalpy");
        return result;
    }

    template <typename Number>
    [[nodiscard]] Number lbc_dynamic_viscosity_pa_s(
        const Number& temperature_k,
        std::span<const Number> composition,
        const Number& mixture_molar_mass_kg_per_mol,
        const Number& molar_density_mol_per_m3) const {
        using namespace
            pr76_li_firoozabadi_sour_gas_property_detail;
        validate_composition(composition);
        require_finite(
            temperature_k,
            "temperature",
            true);
        require_finite(
            mixture_molar_mass_kg_per_mol,
            "mixture molar mass",
            true);
        require_finite(
            molar_density_mol_per_m3,
            "molar density",
            true);

        using std::pow;
        using std::sqrt;

        constexpr double pa_per_psia =
            6894.757293168;
        constexpr double kelvin_to_rankine =
            1.8;

        Number dilute_numerator{0.0};
        Number dilute_denominator{0.0};
        Number pseudo_tc_rankine{0.0};
        Number pseudo_pc_psia{0.0};
        Number pseudo_vc_m3_per_mol{0.0};

        for (std::size_t component = 0U;
             component < component_ids.size();
             ++component) {
            const double tc =
                critical_temperature_k_[component];
            const double pc =
                critical_pressure_pa_[component];
            const double mw_g_per_mol =
                nist_molar_mass_kg_per_mol[component] *
                1000.0;
            const double tc_rankine =
                tc * kelvin_to_rankine;
            const double pc_psia =
                pc / pa_per_psia;
            const double xi =
                5.4402 *
                std::pow(
                    tc_rankine,
                    1.0 / 6.0) /
                (std::sqrt(mw_g_per_mol) *
                 std::pow(
                     pc_psia,
                     2.0 / 3.0));
            if (!std::isfinite(xi) ||
                !(xi > 0.0)) {
                throw std::range_error(
                    "mpmc::flow::invalid Li-Firoozabadi LBC pure-component viscosity parameter");
            }

            const Number tr =
                temperature_k / tc;
            const double tr_primal =
                primal_value(tr);
            Number pure_viscosity_cp{0.0};
            if (tr_primal <= 1.5) {
                pure_viscosity_cp =
                    34.0e-5 *
                    pow(tr, 0.94) /
                    xi;
            } else {
                const Number bracket =
                    4.58 * tr - 1.67;
                require_finite(
                    bracket,
                    "Stiel-Thodos high-Tr bracket",
                    true);
                pure_viscosity_cp =
                    17.78e-5 *
                    pow(
                        bracket,
                        0.625) /
                    xi;
            }
            require_finite(
                pure_viscosity_cp,
                "Stiel-Thodos dilute viscosity",
                true);

            const double sqrt_mw =
                std::sqrt(mw_g_per_mol);
            dilute_numerator +=
                composition[component] *
                pure_viscosity_cp *
                sqrt_mw;
            dilute_denominator +=
                composition[component] *
                sqrt_mw;

            pseudo_tc_rankine +=
                composition[component] *
                tc_rankine;
            pseudo_pc_psia +=
                composition[component] *
                pc_psia;
            pseudo_vc_m3_per_mol +=
                composition[component] *
                critical_molar_volume_m3_per_mol[
                    component];
        }

        const Number dilute_cp =
            dilute_numerator /
            dilute_denominator;

        const Number mixture_mw_lbm_per_lbmol =
            mixture_molar_mass_kg_per_mol *
            1000.0;
        const Number xi_m =
            5.4402 *
            pow(
                pseudo_tc_rankine,
                1.0 / 6.0) /
            (sqrt(mixture_mw_lbm_per_lbmol) *
             pow(
                 pseudo_pc_psia,
                 2.0 / 3.0));
        require_finite(
            xi_m,
            "LBC mixture viscosity parameter",
            true);

        const Number reduced_density =
            molar_density_mol_per_m3 *
            pseudo_vc_m3_per_mol;
        require_finite(
            reduced_density,
            "LBC reduced density",
            true);

        const Number polynomial =
            lbc_polynomial(reduced_density);
        const Number square =
            polynomial * polynomial;
        const Number fourth =
            square * square;
        const Number dense_increment_cp =
            (fourth - 1.0e-4) /
            xi_m;
        const Number viscosity_cp =
            dilute_cp +
            dense_increment_cp;
        require_finite(
            viscosity_cp,
            "LBC phase viscosity",
            true);

        return viscosity_cp * 1.0e-3;
    }

    template <typename Number>
    [[nodiscard]] Number pr76_residual_molar_enthalpy_j_per_mol(
        const thermodynamics::Pr76SelectedPhase& selection,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition) const {
        using namespace
            pr76_li_firoozabadi_sour_gas_property_detail;
        validate_composition(composition);
        require_finite(
            pressure_pa,
            "pressure",
            true);
        require_finite(
            temperature_k,
            "temperature",
            true);

        using TemperatureDual =
            mpmc::ad::Dual<Number, 1U>;
        TemperatureDual nested_temperature =
            TemperatureDual::variable(
                temperature_k,
                0U);
        std::vector<TemperatureDual>
            nested_composition;
        nested_composition.reserve(
            composition.size());
        for (const auto& fraction : composition) {
            nested_composition.emplace_back(
                fraction);
        }

        thermodynamics::
            Pr76MixtureWorkspace<TemperatureDual>
                mixture_workspace;
        const auto nested_mixed =
            mixture_.evaluate_full(
                nested_temperature,
                std::span<const TemperatureDual>{
                    nested_composition},
                mixture_workspace);

        const Number a =
            nested_mixed.a.value();
        const Number b =
            nested_mixed.b.value();
        const Number da_dtemperature =
            nested_mixed.a.derivative(0U);

        thermodynamics::
            Pr76PhaseWorkspace<Number>
                phase_workspace;
        const auto phase =
            model_->evaluate_full(
                pressure_pa,
                temperature_k,
                composition,
                selection.root_index,
                phase_workspace,
                selection.root_options);

        const Number rt =
            thermodynamics::Pr76Pure<double>::
                gas_constant() *
            temperature_k;
        const Number B =
            b * pressure_pa / rt;

        constexpr double sqrt2 =
            1.41421356237309504880168872420969808;
        const Number numerator =
            phase.z +
            (1.0 + sqrt2) * B;
        const Number denominator =
            phase.z +
            (1.0 - sqrt2) * B;
        require_finite(
            numerator,
            "PR enthalpy logarithm numerator",
            true);
        require_finite(
            denominator,
            "PR enthalpy logarithm denominator",
            true);

        using std::log;
        const Number logarithm =
            log(numerator / denominator);
        const Number departure =
            rt * (phase.z - 1.0) +
            (temperature_k *
                 da_dtemperature -
             a) /
                (2.0 * sqrt2 * b) *
                logarithm;
        require_finite(
            departure,
            "PR residual molar enthalpy");
        return departure;
    }

    template <typename Number>
    [[nodiscard]]
    SelectedPhaseTransportCaloricValues<Number>
    operator()(
        const thermodynamics::Pr76SelectedPhase& selection,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition,
        const Number& mixture_molar_mass_kg_per_mol,
        const Number& molar_density_mol_per_m3,
        const Number&) const {
        const Number viscosity =
            lbc_dynamic_viscosity_pa_s(
                temperature_k,
                composition,
                mixture_molar_mass_kg_per_mol,
                molar_density_mol_per_m3);

        const Number ideal_molar =
            ideal_gas_molar_enthalpy_j_per_mol(
                temperature_k,
                composition);
        const Number residual_molar =
            pr76_residual_molar_enthalpy_j_per_mol(
                selection,
                pressure_pa,
                temperature_k,
                composition);
        const Number specific_enthalpy =
            (ideal_molar + residual_molar) /
            mixture_molar_mass_kg_per_mol;

        pr76_li_firoozabadi_sour_gas_property_detail::
            require_finite(
                specific_enthalpy,
                "total specific enthalpy");

        return {
            viscosity,
            specific_enthalpy};
    }

private:
    const thermodynamics::Pr76Phase<double>* model_;
    thermodynamics::Pr76Mixture<double> mixture_;
    std::array<double, 6> critical_temperature_k_{};
    std::array<double, 6> critical_pressure_pa_{};
};

[[nodiscard]] inline auto
make_pr76_li_firoozabadi_sour_gas_property_closure(
    const thermodynamics::Pr76Phase<double>& model,
    std::vector<thermodynamics::Pr76SelectedPhase>
        selections) {
    return make_pr76_selected_phase_property_closure(
        model,
        std::move(selections),
        Pr76LiFiroozabadiSourGasPropertyProvider{
            model},
        Pr76LiFiroozabadiSourGasPropertyProvider::
            provenance());
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_PR76_LI_FIROOZABADI_SOUR_GAS_PROPERTIES_HPP
