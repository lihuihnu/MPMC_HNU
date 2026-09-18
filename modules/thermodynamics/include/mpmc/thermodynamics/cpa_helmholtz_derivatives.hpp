#ifndef MPMC_THERMODYNAMICS_CPA_HELMHOLTZ_DERIVATIVES_HPP
#define MPMC_THERMODYNAMICS_CPA_HELMHOLTZ_DERIVATIVES_HPP

#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/thermodynamics/cpa_residual_helmholtz.hpp>

#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace mpmc::thermodynamics {

struct CpaHelmholtzPressure {
    double pressure_physical_pa{};
    double pressure_association_pa{};
    double pressure_pa{};
};

struct CpaHelmholtzChemicalPotentials {
    std::vector<double> cubic_over_rt;
    std::vector<double> association_over_rt;
    std::vector<double> total_over_rt;
};

// Production first-derivative adapter for the canonical residual Helmholtz
// kernel. The association state is solved once by the caller at the primal
// (T,V,n) state and is held stationary while differentiating the algebraic
// cubic and Q terms.
[[nodiscard]] inline CpaHelmholtzPressure cpa_helmholtz_pressure(
    double temperature_k,
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const CpaParameterSet& parameters,
    const CpaAssociationResult& primal_association) {
    const double total_moles =
        cpa_detail::validate_cpa_composition(composition, parameters.size());
    if (!std::isfinite(temperature_k) || !(temperature_k > 0.0) ||
        !std::isfinite(molar_density_mol_per_m3) ||
        !(molar_density_mol_per_m3 > 0.0)) {
        throw std::domain_error(
            "CPA Helmholtz pressure: positive finite T and molar density required");
    }

    using Number = mpmc::ad::Dual<double, 1>;
    const Number temperature{temperature_k};
    const Number volume = Number::variable(
        total_moles / molar_density_mol_per_m3, 0U);

    std::vector<Number> mole_numbers;
    mole_numbers.reserve(composition.size());
    for (const double value : composition) {
        mole_numbers.emplace_back(value);
    }

    const auto mole_span = std::span<const Number>{mole_numbers};
    const Number cubic = cpa_cubic_residual_helmholtz_reduced(
        temperature, volume, mole_span, parameters);
    const Number association = cpa_association_q_reduced(
        temperature, volume, mole_span, parameters, primal_association);

    const double rt = cpa_gas_constant_j_per_mol_k * temperature_k;
    const double ideal_pressure =
        total_moles * rt / (total_moles / molar_density_mol_per_m3);

    CpaHelmholtzPressure result;
    result.pressure_physical_pa =
        ideal_pressure - rt * cubic.derivative(0U);
    result.pressure_association_pa =
        -rt * association.derivative(0U);
    result.pressure_pa =
        result.pressure_physical_pa + result.pressure_association_pa;

    if (!std::isfinite(result.pressure_physical_pa) ||
        !std::isfinite(result.pressure_association_pa) ||
        !std::isfinite(result.pressure_pa)) {
        throw std::range_error(
            "CPA Helmholtz pressure: nonrepresentable derivative state");
    }
    return result;
}

[[nodiscard]] inline CpaHelmholtzChemicalPotentials
cpa_helmholtz_residual_chemical_potentials(
    double temperature_k,
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const CpaParameterSet& parameters,
    const CpaAssociationResult& primal_association) {
    const double total_moles =
        cpa_detail::validate_cpa_composition(composition, parameters.size());
    if (!std::isfinite(temperature_k) || !(temperature_k > 0.0) ||
        !std::isfinite(molar_density_mol_per_m3) ||
        !(molar_density_mol_per_m3 > 0.0)) {
        throw std::domain_error(
            "CPA Helmholtz chemical potential: positive finite T and molar density required");
    }

    using Number = mpmc::ad::Dual<double, 4>;
    mpmc::ad::RuntimeJacobianWorkspace<double, 4> workspace;
    const double volume_m3 =
        total_moles / molar_density_mol_per_m3;

    const auto callback = [&](std::span<const Number> mole_numbers,
                              std::span<Number> outputs) {
        const Number temperature{temperature_k};
        const Number volume{volume_m3};
        outputs[0] = cpa_cubic_residual_helmholtz_reduced(
            temperature, volume, mole_numbers, parameters);
        outputs[1] = cpa_association_q_reduced(
            temperature, volume, mole_numbers,
            parameters, primal_association);
    };

    const auto result = mpmc::ad::value_and_jacobian_runtime<4>(
        callback,
        composition,
        2U,
        workspace,
        {1024U, 2U, 2048U});

    if (result.input_count != parameters.size() ||
        result.output_count != 2U) {
        throw std::runtime_error(
            "CPA Helmholtz chemical potential: unexpected AD result shape");
    }

    CpaHelmholtzChemicalPotentials derived;
    derived.cubic_over_rt.resize(parameters.size());
    derived.association_over_rt.resize(parameters.size());
    derived.total_over_rt.resize(parameters.size());

    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        const double cubic =
            result.jacobian[i];
        const double association =
            result.jacobian[result.input_count + i];
        const double total = cubic + association;
        if (!std::isfinite(cubic) ||
            !std::isfinite(association) ||
            !std::isfinite(total)) {
            throw std::range_error(
                "CPA Helmholtz chemical potential: nonrepresentable derivative state");
        }
        derived.cubic_over_rt[i] = cubic;
        derived.association_over_rt[i] = association;
        derived.total_over_rt[i] = total;
    }

    return derived;
}

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_CPA_HELMHOLTZ_DERIVATIVES_HPP
