#ifndef MPMC_FLOW_SW92_CO2_WATER_PROPERTIES_HPP
#define MPMC_FLOW_SW92_CO2_WATER_PROPERTIES_HPP

#include <mpmc/ad/math.hpp>
#include <mpmc/flow/sw92_selected_phase_property_closure.hpp>
#include <mpmc/thermodynamics/sw92_mixture.hpp>

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
    sw92_co2_water_property_convention =
        "flow/sw92/co2-water/zero-salinity/Chung1988-NIST-SW92-departure/v1";

namespace sw92_co2_water_property_detail {

inline constexpr std::string_view
    required_dataset_id =
        "SW92-Table3-Table5-corrected";
inline constexpr std::string_view
    required_revision =
        "user-PDF-sha256-cb5b1d50";

inline constexpr std::array<std::string_view, 2>
    canonical_component_ids{
        "carbon-dioxide",
        "water"};

/// NIST Chemistry WebBook SRD 69 molecular weights [kg/mol].
inline constexpr std::array<double, 2>
    molar_mass_kg_per_mol{
        0.0440095,
        0.0180153};

/// Chung transport critical molar volumes [cm3/mol].
/// CO2: NIST SRD 69, Li & Kiran (1988), Vc=0.0919 L/mol.
/// H2O: IAPWS critical density 322 kg/m3 combined with the NIST molar mass.
inline constexpr std::array<double, 2>
    critical_molar_volume_cm3_per_mol{
        91.9,
        55.94813664596274};

/// Chung polar parameters.
/// CO2 has zero permanent dipole and no association correction.
/// H2O values follow the Chung/Poling parameterization used in the public
/// IDAES Chung reference example: dipole=1.8546 D, kappa=0.076.
inline constexpr std::array<double, 2>
    dipole_debye{
        0.0,
        1.8546};
inline constexpr std::array<double, 2>
    association_factor{
        0.0,
        0.076};

/// Repository SW92 Table-3 values used by this provider.
/// They are checked, not substituted into the thermodynamic model.
inline constexpr std::array<double, 2>
    sw92_critical_temperature_k{
        304.2,
        647.3};
inline constexpr std::array<double, 2>
    sw92_critical_pressure_pa{
        73.8e5,
        221.2e5};
inline constexpr std::array<double, 2>
    sw92_acentric_factor{
        0.2273,
        0.3434};

struct ShomateSensibleEnthalpy {
    double a;
    double b;
    double c;
    double d;
    double e;
    double f;
    double h;
};

/// NIST SRD 69 / Chase (1998) gas-phase Shomate coefficients.
/// H-H_298.15 is in kJ/mol.  The common sourced interval is deliberately
/// restricted to 500--1200 K: CO2 is tabulated from 298--1200 K while H2O's
/// first gas-phase Shomate interval begins at 500 K.
inline constexpr std::array<
    ShomateSensibleEnthalpy,
    2>
    nist_shomate{{
        {
            24.99735,
            55.18696,
            -33.69137,
            7.948387,
            -0.136638,
            -403.6075,
            -393.5224},
        {
            30.09200,
            6.832514,
            6.793435,
            -2.534480,
            0.082139,
            -250.8810,
            -241.8264}}};

inline constexpr double
    minimum_temperature_k = 500.0;
inline constexpr double
    maximum_temperature_k = 1200.0;

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
                "mpmc::flow::SW92 CO2/H2O property: "} +
            name +
            (positive
                 ? " must be finite and strictly positive"
                 : " must be finite"));
    }
}

inline void require_source_value(
    double actual,
    double expected,
    double absolute_tolerance,
    const char* name) {
    if (!std::isfinite(actual) ||
        std::abs(actual - expected) >
            absolute_tolerance) {
        throw std::invalid_argument(
            std::string{
                "mpmc::flow::SW92 CO2/H2O property: sourced "} +
            name +
            " does not match the curated provider snapshot");
    }
}

template <typename Number>
inline void validate_composition(
    std::span<const Number> composition) {
    if (composition.size() != 2U) {
        throw std::invalid_argument(
            "mpmc::flow::SW92 CO2/H2O property requires exactly two components");
    }
    long double sum = 0.0L;
    for (const auto& fraction : composition) {
        const double primal = primal_value(fraction);
        if (!std::isfinite(primal) ||
            !(primal > 0.0) ||
            !(primal < 1.0)) {
            throw std::domain_error(
                "mpmc::flow::SW92 CO2/H2O property requires strict-positive interior mole fractions");
        }
        sum += static_cast<long double>(primal);
    }
    if (std::abs(sum - 1.0L) > 2.0e-12L) {
        throw std::domain_error(
            "mpmc::flow::SW92 CO2/H2O property requires normalized mole fractions");
    }
}

template <typename Number>
inline void validate_temperature(
    const Number& temperature_k) {
    require_finite(
        temperature_k,
        "temperature [K]",
        true);
    const double primal =
        primal_value(temperature_k);
    if (primal < minimum_temperature_k ||
        primal > maximum_temperature_k) {
        throw std::domain_error(
            "mpmc::flow::SW92 CO2/H2O NIST caloric provider is validated only on [500,1200] K");
    }
}

inline void validate_zero_salinity(
    const thermodynamics::Sw92SelectedPhase<double>&
        selection) {
    if (selection.nacl_molality_mol_per_kg_water !=
        0.0) {
        throw std::domain_error(
            "mpmc::flow::SW92 CO2/H2O property provider supports only zero NaCl molality; no brine transport/caloric model is installed");
    }
}

template <typename Number>
[[nodiscard]] inline Number
shomate_sensible_molar_enthalpy_j_per_mol(
    const Number& temperature_k,
    const ShomateSensibleEnthalpy& c) {
    const Number t =
        temperature_k / 1000.0;
    const Number t2 = t * t;
    const Number t3 = t2 * t;
    const Number t4 = t3 * t;
    const Number h_kj_per_mol =
        c.a * t +
        0.5 * c.b * t2 +
        (c.c / 3.0) * t3 +
        0.25 * c.d * t4 -
        c.e / t +
        c.f -
        c.h;
    require_finite(
        h_kj_per_mol,
        "NIST Shomate sensible molar enthalpy");
    return h_kj_per_mol * 1000.0;
}

template <typename Number>
[[nodiscard]] inline Number
neufeld_collision_integral(
    const Number& reduced_temperature) {
    require_finite(
        reduced_temperature,
        "Chung reduced temperature",
        true);
    const double primal =
        primal_value(reduced_temperature);
    if (primal < 0.3 || primal > 100.0) {
        throw std::domain_error(
            "mpmc::flow::SW92 CO2/H2O Chung collision integral requires 0.3 <= T* <= 100");
    }
    using std::exp;
    using std::pow;
    const Number omega =
        1.16145 /
            pow(reduced_temperature, 0.14874) +
        0.52487 /
            exp(0.77320 * reduced_temperature) +
        2.16178 /
            exp(2.43787 * reduced_temperature);
    require_finite(
        omega,
        "Chung collision integral",
        true);
    return omega;
}

} // namespace sw92_co2_water_property_detail

/// Source-backed transport/caloric provider for the repository-curated
/// zero-salinity SW92 CO2/H2O binary.
///
/// Viscosity:
///   Chung, Ajlan, Lee & Starling (1988), DOI 10.1021/ie00076a024,
///   high-pressure polar-mixture corresponding-states correlation.  The
///   implementation follows the Poling 5th-edition equation numbering
///   9-5.25..44 and 9-6.18..23, including the Neufeld collision integral.
///   The Chung xi/zeta binary transport parameters are unity.
///
/// Caloric:
///   NIST SRD 69 / Chase (1998) ideal-gas Shomate H-H_298.15 on the common
///   500--1200 K interval plus the Peng-Robinson departure enthalpy evaluated
///   with the *SW92 family-specific* a(T,x), b(x), water alpha and BIPs.
///
/// The enthalpy reference is intentionally sensible H-H_298.15; it is not a
/// heat-of-formation convention.  NaCl molality is required to be exactly zero:
/// this provider does not extrapolate pure-water transport/caloric parameters
/// into brine.
class Sw92Co2WaterPropertyProvider {
public:
    explicit Sw92Co2WaterPropertyProvider(
        const thermodynamics::Sw92Phase<double>&
            model)
        : model_(&model),
          mixture_(
              thermodynamics::Sw92Mixture<double>::
                  from_parameters(
                      model.parameters())) {
        using namespace
            sw92_co2_water_property_detail;
        const auto& parameters =
            model.parameters();
        if (parameters.dataset_id() !=
                required_dataset_id ||
            parameters.revision() !=
                required_revision ||
            parameters.components().size() != 2U) {
            throw std::invalid_argument(
                "mpmc::flow::SW92 CO2/H2O property provider only supports the repository-curated corrected-original binary snapshot");
        }

        std::array<bool, 2> seen{false, false};
        for (std::size_t model_index = 0U;
             model_index < 2U;
             ++model_index) {
            const auto& component =
                parameters.components().at(
                    model_index);
            std::size_t canonical = 2U;
            for (std::size_t candidate = 0U;
                 candidate < 2U;
                 ++candidate) {
                if (component.id ==
                    canonical_component_ids[
                        candidate]) {
                    canonical = candidate;
                    break;
                }
            }
            if (canonical >= 2U ||
                seen[canonical]) {
                throw std::invalid_argument(
                    "mpmc::flow::SW92 CO2/H2O property provider requires exactly carbon-dioxide and water");
            }
            seen[canonical] = true;
            model_to_canonical_[model_index] =
                canonical;

            if (!component.molar_mass ||
                component.molar_mass->unit !=
                    thermodynamics::Unit::
                        kilogram_per_mole) {
                throw std::invalid_argument(
                    "mpmc::flow::SW92 CO2/H2O property provider requires sourced kg/mol molar masses");
            }
            require_source_value(
                component.molar_mass->value,
                molar_mass_kg_per_mol[
                    canonical],
                5.0e-10,
                "molar mass");
            require_source_value(
                parameters
                    .critical_temperatures_k()[
                        model_index],
                sw92_critical_temperature_k[
                    canonical],
                5.0e-10,
                "critical temperature");
            require_source_value(
                parameters
                    .critical_pressures_pa()[
                        model_index],
                sw92_critical_pressure_pa[
                    canonical],
                5.0e-4,
                "critical pressure");
            require_source_value(
                parameters.acentric_factors()[
                    model_index],
                sw92_acentric_factor[
                    canonical],
                5.0e-12,
                "acentric factor");
        }
    }

    [[nodiscard]] static
    SelectedPhasePropertyProvenance
    provenance() {
        return {
            {
                "SW92 selected molar density x NIST SRD 69 molar mass",
                "SW92-Table3-Table5-corrected__NIST-SRD69-CO2-H2O-MW",
                "zero-salinity-v1"},
            {
                "Chung-Ajlan-Lee-Starling 1988 high-pressure polar-mixture viscosity",
                "doi:10.1021/ie00076a024__NIST-CO2-Vc__IAPWS-H2O-rhoc__Chung-H2O-polar",
                "zero-salinity-v1"},
            {
                "NIST SRD 69 Shomate sensible enthalpy + SW92 family-specific PR departure",
                "NIST-SRD69-Chase1998-CO2-H2O-Shomate__SW92-corrected-original",
                "500-1200K-zero-salinity-v1"},
            {
                "thermodynamic identity u=h-p/rho_mass",
                "SW92-Chung-NIST-zero-salinity",
                "v1"}};
    }

    template <typename Number>
    [[nodiscard]] Number
    ideal_gas_molar_enthalpy_j_per_mol(
        const Number& temperature_k,
        std::span<const Number> composition)
        const {
        using namespace
            sw92_co2_water_property_detail;
        validate_temperature(temperature_k);
        validate_composition(composition);
        Number result{0.0};
        for (std::size_t model_index = 0U;
             model_index < 2U;
             ++model_index) {
            const std::size_t canonical =
                model_to_canonical_[model_index];
            result +=
                composition[model_index] *
                shomate_sensible_molar_enthalpy_j_per_mol(
                    temperature_k,
                    nist_shomate[canonical]);
        }
        require_finite(
            result,
            "mixture ideal-gas sensible molar enthalpy");
        return result;
    }

    template <typename Number>
    [[nodiscard]] Number
    chung_dynamic_viscosity_pa_s(
        const Number& temperature_k,
        std::span<const Number> composition,
        const Number& molar_density_mol_per_m3)
        const {
        using namespace
            sw92_co2_water_property_detail;
        validate_temperature(temperature_k);
        validate_composition(composition);
        require_finite(
            molar_density_mol_per_m3,
            "molar density [mol/m3]",
            true);

        std::array<double, 2> sigma{};
        std::array<double, 2> epsilon_over_k{};
        std::array<double, 2> molecular_weight_g_per_mol{};
        std::array<double, 2> omega{};
        std::array<double, 2> dipole{};
        std::array<double, 2> kappa{};

        for (std::size_t model_index = 0U;
             model_index < 2U;
             ++model_index) {
            const std::size_t canonical =
                model_to_canonical_[model_index];
            const double vc =
                critical_molar_volume_cm3_per_mol[
                    canonical];
            sigma[model_index] =
                0.809 * std::cbrt(vc);
            epsilon_over_k[model_index] =
                model_->parameters()
                    .critical_temperatures_k()[
                        model_index] /
                1.2593;
            molecular_weight_g_per_mol[
                model_index] =
                molar_mass_kg_per_mol[
                    canonical] *
                1000.0;
            omega[model_index] =
                model_->parameters()
                    .acentric_factors()[
                        model_index];
            dipole[model_index] =
                dipole_debye[canonical];
            kappa[model_index] =
                association_factor[canonical];
        }

        Number sigma_mix_cubed{0.0};
        Number epsilon_numerator{0.0};
        Number molecular_weight_numerator{0.0};
        Number omega_numerator{0.0};
        Number dipole_sum{0.0};
        Number kappa_mix{0.0};

        for (std::size_t i = 0U;
             i < 2U;
             ++i) {
            for (std::size_t j = 0U;
                 j < 2U;
                 ++j) {
                const Number xij =
                    composition[i] *
                    composition[j];
                const double sigma_ij =
                    std::sqrt(
                        sigma[i] * sigma[j]);
                const double sigma_ij2 =
                    sigma_ij * sigma_ij;
                const double sigma_ij3 =
                    sigma_ij2 * sigma_ij;
                const double epsilon_ij =
                    std::sqrt(
                        epsilon_over_k[i] *
                        epsilon_over_k[j]);
                const double mw_ij =
                    2.0 *
                    molecular_weight_g_per_mol[i] *
                    molecular_weight_g_per_mol[j] /
                    (molecular_weight_g_per_mol[i] +
                     molecular_weight_g_per_mol[j]);
                const double omega_ij =
                    0.5 *
                    (omega[i] + omega[j]);
                const double kappa_ij =
                    std::sqrt(
                        kappa[i] * kappa[j]);

                sigma_mix_cubed +=
                    xij * sigma_ij3;
                epsilon_numerator +=
                    xij * epsilon_ij *
                    sigma_ij3;
                molecular_weight_numerator +=
                    xij * epsilon_ij *
                    sigma_ij2 *
                    std::sqrt(mw_ij);
                omega_numerator +=
                    xij * omega_ij *
                    sigma_ij3;
                dipole_sum +=
                    xij *
                    ((dipole[i] * dipole[i]) *
                     (dipole[j] * dipole[j]) /
                     sigma_ij3);
                kappa_mix +=
                    xij * kappa_ij;
            }
        }

        require_finite(
            sigma_mix_cubed,
            "Chung sigma_mix^3",
            true);
        using std::cbrt;
        using std::exp;
        using std::expm1;
        using std::pow;
        using std::sqrt;
        const Number sigma_mix =
            cbrt(sigma_mix_cubed);
        const Number epsilon_mix =
            epsilon_numerator /
            sigma_mix_cubed;
        const Number mw_scale =
            molecular_weight_numerator /
            (epsilon_mix *
             sigma_mix * sigma_mix);
        const Number molecular_weight_mix =
            mw_scale * mw_scale;
        const Number omega_mix =
            omega_numerator /
            sigma_mix_cubed;
        const Number dipole_mix_fourth =
            sigma_mix_cubed *
            dipole_sum;
        require_finite(
            dipole_mix_fourth,
            "Chung dipole mixture fourth power",
            true);
        const Number dipole_mix =
            sqrt(sqrt(dipole_mix_fourth));
        const Number critical_volume_mix =
            (sigma_mix / 0.809) *
            (sigma_mix / 0.809) *
            (sigma_mix / 0.809);
        const Number critical_temperature_mix =
            1.2593 * epsilon_mix;
        const Number reduced_dipole =
            131.3 * dipole_mix /
            sqrt(
                critical_volume_mix *
                critical_temperature_mix);
        const Number reduced_temperature =
            temperature_k /
            epsilon_mix;
        const Number density_mol_per_cm3 =
            molar_density_mol_per_m3 *
            1.0e-6;
        const Number y =
            density_mol_per_cm3 *
            critical_volume_mix /
            6.0;

        require_finite(
            molecular_weight_mix,
            "Chung mixture molecular weight",
            true);
        require_finite(
            critical_volume_mix,
            "Chung mixture critical volume",
            true);
        require_finite(
            critical_temperature_mix,
            "Chung mixture critical temperature",
            true);
        require_finite(
            reduced_dipole,
            "Chung mixture reduced dipole");
        require_finite(
            y,
            "Chung reduced density",
            true);
        if (!(primal_value(y) < 1.0)) {
            throw std::domain_error(
                "mpmc::flow::SW92 CO2/H2O Chung viscosity requires reduced density y < 1");
        }

        const Number reduced_dipole2 =
            reduced_dipole * reduced_dipole;
        const Number reduced_dipole4 =
            reduced_dipole2 *
            reduced_dipole2;

        constexpr std::array<double, 10> a{
            6.324,
            1.210e-3,
            5.283,
            6.623,
            19.745,
            -1.900,
            24.275,
            0.7972,
            -0.2382,
            0.06863};
        constexpr std::array<double, 10> b{
            50.412,
            -1.154e-3,
            254.209,
            38.096,
            7.630,
            -12.537,
            3.450,
            1.117,
            0.06770,
            0.3479};
        constexpr std::array<double, 10> c{
            -51.680,
            -6.257e-3,
            -168.48,
            -8.464,
            -14.354,
            4.958,
            -11.291,
            0.01235,
            -0.8163,
            0.5926};
        constexpr std::array<double, 10> d{
            1189.0,
            0.03728,
            3898.0,
            31.42,
            31.53,
            -18.15,
            69.35,
            -4.117,
            4.025,
            -0.727};

        std::array<Number, 10> e{};
        for (std::size_t i = 0U;
             i < e.size();
             ++i) {
            e[i] =
                a[i] +
                b[i] * omega_mix +
                c[i] * reduced_dipole4 +
                d[i] * kappa_mix;
        }

        const Number one_minus_y =
            1.0 - y;
        require_finite(
            one_minus_y,
            "Chung 1-y",
            true);
        const Number g1 =
            (1.0 - 0.5 * y) /
            (one_minus_y *
             one_minus_y *
             one_minus_y);
        const Number exprel =
            -expm1(-e[3] * y) / y;
        const Number g2_denominator =
            e[0] * e[3] +
            e[1] +
            e[2];
        require_finite(
            g2_denominator,
            "Chung G2 denominator");
        if (std::abs(
                primal_value(
                    g2_denominator)) <
            64.0 *
                std::numeric_limits<double>::
                    epsilon()) {
            throw std::range_error(
                "mpmc::flow::SW92 CO2/H2O Chung G2 denominator is singular");
        }
        const Number g2 =
            (e[0] * exprel +
             e[1] * g1 * exp(e[4] * y) +
             e[2] * g1) /
            g2_denominator;
        require_finite(
            g2,
            "Chung G2",
            true);

        const Number eta2 =
            e[6] * y * y * g2 *
            exp(
                e[7] +
                e[8] / reduced_temperature +
                e[9] /
                    (reduced_temperature *
                     reduced_temperature));
        const Number collision_integral =
            neufeld_collision_integral(
                reduced_temperature);
        const Number fc =
            1.0 -
            0.2756 * omega_mix +
            0.059035 * reduced_dipole4 +
            kappa_mix;
        const Number eta_star =
            sqrt(reduced_temperature) /
                collision_integral *
                (fc *
                     (1.0 / g2 +
                      e[5] * y)) +
            eta2;
        const Number viscosity_micropoise =
            eta_star *
            36.344 *
            sqrt(
                molecular_weight_mix *
                critical_temperature_mix) /
            pow(
                critical_volume_mix,
                2.0 / 3.0);
        const Number viscosity_pa_s =
            viscosity_micropoise *
            1.0e-7;
        require_finite(
            viscosity_pa_s,
            "Chung high-pressure mixture viscosity [Pa s]",
            true);
        return viscosity_pa_s;
    }

    template <typename Number>
    [[nodiscard]] Number
    sw92_attraction_temperature_derivative(
        const thermodynamics::Sw92SelectedPhase<double>&
            selection,
        const Number& temperature_k,
        std::span<const Number> composition)
        const {
        using namespace
            sw92_co2_water_property_detail;
        validate_zero_salinity(selection);
        validate_temperature(temperature_k);
        validate_composition(composition);

        const auto& parameters =
            model_->parameters();
        constexpr double gas_constant =
            8.31446261815324;

        std::array<Number, 2> pure_a{};
        std::array<Number, 2> pure_da_dt{};

        using std::exp;
        using std::sqrt;

        for (std::size_t model_index = 0U;
             model_index < 2U;
             ++model_index) {
            const std::size_t canonical =
                model_to_canonical_[model_index];
            const double tc =
                parameters
                    .critical_temperatures_k()[
                        model_index];
            const double pc =
                parameters
                    .critical_pressures_pa()[
                        model_index];
            const double rtc =
                gas_constant * tc;
            const double ac =
                0.45724 * rtc * rtc / pc;

            if (canonical == 0U) {
                const double omega =
                    parameters.acentric_factors()[
                        model_index];
                const double kappa =
                    0.37464 +
                    omega *
                        (1.54226 -
                         0.26992 * omega);
                const Number root_temperature =
                    sqrt(temperature_k);
                const double root_tc =
                    std::sqrt(tc);
                const Number q =
                    1.0 +
                    kappa *
                        (1.0 -
                         root_temperature /
                             root_tc);
                pure_a[model_index] =
                    ac * q * q;
                pure_da_dt[model_index] =
                    -ac * q * kappa /
                    (root_temperature *
                     root_tc);
            } else {
                const Number tr =
                    temperature_k / tc;
                const Number tr2 = tr * tr;
                const Number tr3 = tr2 * tr;
                const Number tr4 = tr2 * tr2;
                const Number q =
                    1.0 +
                    0.4530 * (1.0 - tr) +
                    0.0034 *
                        (1.0 / tr3 - 1.0);
                const Number dq_dt =
                    (-0.4530 -
                     0.0102 / tr4) /
                    tc;
                pure_a[model_index] =
                    ac * q * q;
                pure_da_dt[model_index] =
                    2.0 * ac * q * dq_dt;
            }
            require_finite(
                pure_a[model_index],
                "SW92 pure attraction",
                true);
            require_finite(
                pure_da_dt[model_index],
                "SW92 pure attraction temperature derivative");
        }

        const std::size_t co2_index =
            model_to_canonical_[0] == 0U
                ? 0U
                : 1U;
        const std::size_t water_index =
            1U - co2_index;

        Number kij{0.0};
        Number dkij_dt{0.0};
        if (selection.family ==
            thermodynamics::SwPhaseFamily::
                nonaqueous) {
            kij =
                parameters
                    .water_nonaqueous_constant_kij(
                        co2_index);
        } else if (
            selection.family ==
            thermodynamics::SwPhaseFamily::
                aqueous) {
            const double co2_tc =
                parameters
                    .critical_temperatures_k()[
                        co2_index];
            const Number tr =
                temperature_k / co2_tc;
            const Number exponential =
                exp(-6.7222 * tr);
            kij =
                -0.31092 +
                0.23580 * tr -
                21.2566 * exponential;
            dkij_dt =
                (0.23580 +
                 21.2566 * 6.7222 *
                     exponential) /
                co2_tc;
        } else {
            throw std::invalid_argument(
                "mpmc::flow::SW92 CO2/H2O property provider received an invalid phase family");
        }
        require_finite(
            kij,
            "SW92 CO2/water kij");
        require_finite(
            dkij_dt,
            "SW92 CO2/water dkij/dT");

        const Number cross_attraction =
            sqrt(
                pure_a[co2_index] *
                pure_a[water_index]);
        const Number cross_da_dt =
            0.5 * cross_attraction *
            (pure_da_dt[co2_index] /
                 pure_a[co2_index] +
             pure_da_dt[water_index] /
                 pure_a[water_index]);

        const Number x_co2 =
            composition[co2_index];
        const Number x_water =
            composition[water_index];
        const Number derivative =
            x_co2 * x_co2 *
                pure_da_dt[co2_index] +
            x_water * x_water *
                pure_da_dt[water_index] +
            2.0 * x_co2 * x_water *
                (cross_da_dt *
                     (1.0 - kij) -
                 cross_attraction *
                     dkij_dt);
        require_finite(
            derivative,
            "SW92 mixture attraction temperature derivative");
        return derivative;
    }

    template <typename Number>
    [[nodiscard]] Number
    sw92_residual_molar_enthalpy_j_per_mol(
        const thermodynamics::Sw92SelectedPhase<double>&
            selection,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition)
        const {
        using namespace
            sw92_co2_water_property_detail;
        validate_zero_salinity(selection);
        validate_temperature(temperature_k);
        validate_composition(composition);
        require_finite(
            pressure_pa,
            "pressure [Pa]",
            true);

        thermodynamics::
            Sw92MixtureWorkspace<Number>
                mixture_workspace;
        const auto mixed =
            mixture_.evaluate_full(
                temperature_k,
                composition,
                selection
                    .nacl_molality_mol_per_kg_water,
                selection.family,
                mixture_workspace);

        const Number a = mixed.a;
        const Number b = mixed.b;
        const Number da_dtemperature =
            sw92_attraction_temperature_derivative(
                selection,
                temperature_k,
                composition);

        thermodynamics::
            Sw92PhaseWorkspace<Number>
                phase_workspace;
        const auto phase =
            model_->evaluate_full(
                pressure_pa,
                temperature_k,
                composition,
                selection
                    .nacl_molality_mol_per_kg_water,
                selection.family,
                selection.root_index,
                phase_workspace,
                selection.root_options);

        const Number rt =
            thermodynamics::Sw92Pure<double>::
                gas_constant() *
            temperature_k;
        const Number capital_b =
            b * pressure_pa / rt;

        constexpr double sqrt2 =
            1.41421356237309504880168872420969808;
        const Number numerator =
            phase.z +
            (1.0 + sqrt2) * capital_b;
        const Number denominator =
            phase.z +
            (1.0 - sqrt2) * capital_b;
        require_finite(
            numerator,
            "SW92 departure logarithm numerator",
            true);
        require_finite(
            denominator,
            "SW92 departure logarithm denominator",
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
            "SW92 residual molar enthalpy");
        return departure;
    }

    template <typename Number>
    [[nodiscard]]
    SelectedPhaseTransportCaloricValues<Number>
    operator()(
        const thermodynamics::Sw92SelectedPhase<double>&
            selection,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number> composition,
        const Number&
            mixture_molar_mass_kg_per_mol,
        const Number&
            molar_density_mol_per_m3,
        const Number&) const {
        using namespace
            sw92_co2_water_property_detail;
        validate_zero_salinity(selection);
        validate_temperature(temperature_k);
        validate_composition(composition);
        require_finite(
            mixture_molar_mass_kg_per_mol,
            "mixture molar mass [kg/mol]",
            true);

        Number expected_molar_mass{0.0};
        for (std::size_t model_index = 0U;
             model_index < 2U;
             ++model_index) {
            expected_molar_mass +=
                composition[model_index] *
                molar_mass_kg_per_mol[
                    model_to_canonical_[
                        model_index]];
        }
        const double expected_primal =
            primal_value(
                expected_molar_mass);
        const double provided_primal =
            primal_value(
                mixture_molar_mass_kg_per_mol);
        if (std::abs(
                expected_primal -
                provided_primal) >
            2.0e-10 *
                std::max(
                    {1.0,
                     std::abs(expected_primal),
                     std::abs(provided_primal)})) {
            throw std::invalid_argument(
                "mpmc::flow::SW92 CO2/H2O provider received a molar mass inconsistent with the sourced component snapshot");
        }

        const Number viscosity =
            chung_dynamic_viscosity_pa_s(
                temperature_k,
                composition,
                molar_density_mol_per_m3);
        const Number ideal =
            ideal_gas_molar_enthalpy_j_per_mol(
                temperature_k,
                composition);
        const Number residual =
            sw92_residual_molar_enthalpy_j_per_mol(
                selection,
                pressure_pa,
                temperature_k,
                composition);
        const Number specific_enthalpy =
            (ideal + residual) /
            mixture_molar_mass_kg_per_mol;
        require_finite(
            specific_enthalpy,
            "total specific enthalpy [J/kg]");
        return {
            viscosity,
            specific_enthalpy};
    }

private:
    const thermodynamics::Sw92Phase<double>*
        model_;
    thermodynamics::Sw92Mixture<double>
        mixture_;
    std::array<std::size_t, 2>
        model_to_canonical_{};
};

template <typename SelectionRange>
[[nodiscard]] inline auto
make_sw92_co2_water_selected_phase_property_closure(
    const thermodynamics::Sw92Phase<double>& model,
    SelectionRange&& selections) {
    std::vector<
        thermodynamics::Sw92SelectedPhase<double>>
        owned{
            std::forward<SelectionRange>(
                selections)};
    for (const auto& selection : owned) {
        sw92_co2_water_property_detail::
            validate_zero_salinity(selection);
    }
    Sw92Co2WaterPropertyProvider provider{
        model};
    return make_sw92_selected_phase_property_closure(
        model,
        std::move(owned),
        std::move(provider),
        Sw92Co2WaterPropertyProvider::
            provenance());
}

} // namespace mpmc::flow

#endif // MPMC_FLOW_SW92_CO2_WATER_PROPERTIES_HPP
