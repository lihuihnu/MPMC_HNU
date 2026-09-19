#ifndef MPMC_THERMODYNAMICS_CPA_PT_PHASE_AD_HPP
#define MPMC_THERMODYNAMICS_CPA_PT_PHASE_AD_HPP

#include <mpmc/ad/dual.hpp>
#include <mpmc/ad/math.hpp>
#include <mpmc/thermodynamics/cpa_pt_phase.hpp>
#include <mpmc/thermodynamics/cpa_residual_helmholtz.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::thermodynamics {
namespace cpa_detail {

struct CpaStateDirection {
    double temperature{};
    double volume{};
    std::vector<double> mole_numbers;
};

template <class Scalar>
[[nodiscard]] std::vector<Scalar>
cpa_association_stationarity_residuals(
    const Scalar& temperature_k,
    const Scalar& volume_m3,
    std::span<const Scalar> mole_numbers,
    const CpaParameterSet& parameters,
    const CpaAssociationResult& association,
    std::span<const Scalar> site_fractions) {
    if (mole_numbers.size() != parameters.size() ||
        site_fractions.size() != association.sites.size()) {
        throw std::invalid_argument(
            "CPA association derivative: state dimension mismatch");
    }
    if (association.sites.empty()) { return {}; }

    const Scalar extensive_b =
        cpa_extensive_b(mole_numbers, parameters);
    const Scalar radial_distribution =
        Scalar{1.0} /
        (Scalar{1.0} -
         Scalar{0.475} * extensive_b / volume_m3);
    const Scalar rt =
        Scalar{cpa_gas_constant_j_per_mol_k} *
        temperature_k;

    std::vector<Scalar> residuals(
        association.sites.size(),
        Scalar{});
    for (std::size_t first_index = 0U;
         first_index < association.sites.size();
         ++first_index) {
        const auto& first =
            association.sites[first_index];
        const Scalar first_x =
            site_fractions[first_index];
        Scalar residual =
            Scalar{1.0} / first_x -
            Scalar{1.0};

        for (std::size_t second_index = 0U;
             second_index < association.sites.size();
             ++second_index) {
            const auto& second =
                association.sites[second_index];
            const auto* pair =
                parameters.association_pair(
                    first.component_index,
                    first.site_id,
                    second.component_index,
                    second.site_id);
            if (pair == nullptr ||
                pair->beta_dimensionless == 0.0) {
                continue;
            }
            const double b_ij =
                0.5 *
                (parameters
                     .pure(first.component_index)
                     .b_m3_per_mol +
                 parameters
                     .pure(second.component_index)
                     .b_m3_per_mol);
            const Scalar delta =
                radial_distribution *
                cpa_scalar_expm1(
                    Scalar{
                        pair->epsilon_j_per_mol} /
                    rt) *
                Scalar{
                    b_ij *
                    pair->beta_dimensionless};

            residual -=
                mole_numbers[
                    second.component_index] *
                Scalar{
                    static_cast<double>(
                        second.multiplicity)} *
                site_fractions[second_index] *
                delta /
                volume_m3;
        }
        residuals[first_index] =
            std::move(residual);
    }
    return residuals;
}

class CpaAssociationSensitivitySystem {
public:
    CpaAssociationSensitivitySystem(
        double temperature_k,
        double volume_m3,
        std::span<const double> mole_numbers,
        const CpaParameterSet& parameters,
        const CpaAssociationResult& association)
        : temperature_k_(temperature_k),
          volume_m3_(volume_m3),
          mole_numbers_(
              mole_numbers.begin(),
              mole_numbers.end()),
          parameters_(&parameters),
          association_(&association) {
        if (!association.converged()) {
            throw std::invalid_argument(
                "CPA association derivative: converged association state required");
        }
        if (!std::isfinite(temperature_k_) ||
            !(temperature_k_ > 0.0) ||
            !std::isfinite(volume_m3_) ||
            !(volume_m3_ > 0.0) ||
            mole_numbers_.size() !=
                parameters.size()) {
            throw std::invalid_argument(
                "CPA association derivative: invalid primal state");
        }
        if (association.sites.empty()) {
            return;
        }
        build_and_factorize();
    }

    [[nodiscard]] std::size_t
    site_count() const noexcept {
        return association_->sites.size();
    }

    [[nodiscard]] std::vector<double>
    solve_direction(
        const CpaStateDirection& direction) const {
        if (direction.mole_numbers.size() !=
            mole_numbers_.size()) {
            throw std::invalid_argument(
                "CPA association derivative: direction dimension mismatch");
        }
        if (site_count() == 0U) {
            return {};
        }

        using D =
            mpmc::ad::Dual<double, 1U>;
        const D temperature{
            temperature_k_,
            D::Gradient{
                direction.temperature}};
        const D volume{
            volume_m3_,
            D::Gradient{
                direction.volume}};

        std::vector<D> moles;
        moles.reserve(mole_numbers_.size());
        for (std::size_t i = 0U;
             i < mole_numbers_.size();
             ++i) {
            moles.emplace_back(
                mole_numbers_[i],
                D::Gradient{
                    direction.mole_numbers[i]});
        }

        std::vector<D> sites;
        sites.reserve(site_count());
        for (const auto& site :
             association_->sites) {
            sites.emplace_back(
                site.unbonded_fraction);
        }

        const auto residuals =
            cpa_association_stationarity_residuals(
                temperature,
                volume,
                std::span<const D>{moles},
                *parameters_,
                *association_,
                std::span<const D>{sites});

        std::vector<double> rhs(
            site_count(),
            0.0);
        for (std::size_t i = 0U;
             i < site_count();
             ++i) {
            const double derivative =
                residuals[i].derivative(0U);
            if (!std::isfinite(derivative)) {
                throw std::range_error(
                    "CPA association derivative: non-finite stationarity derivative");
            }
            rhs[i] = -derivative;
        }
        return solve_factored(
            std::move(rhs));
    }

private:
    void build_and_factorize() {
        const std::size_t n =
            site_count();
        lu_.assign(n * n, 0.0);
        pivot_rows_.assign(n, 0U);

        const double total_moles =
            [&] {
                double total = 0.0;
                for (double value :
                     mole_numbers_) {
                    total += value;
                }
                return total;
            }();
        if (!std::isfinite(total_moles) ||
            !(total_moles > 0.0)) {
            throw std::invalid_argument(
                "CPA association derivative: non-positive total moles");
        }

        const double extensive_b =
            cpa_extensive_b(
                std::span<const double>{
                    mole_numbers_},
                *parameters_);
        const double denominator =
            1.0 -
            0.475 * extensive_b /
                volume_m3_;
        if (!std::isfinite(denominator) ||
            !(denominator > 0.0)) {
            throw std::domain_error(
                "CPA association derivative: radial-distribution singularity");
        }
        const double g =
            1.0 / denominator;

        for (std::size_t a = 0U;
             a < n;
             ++a) {
            const auto& first =
                association_->sites[a];
            const double first_x =
                first.unbonded_fraction;
            if (!std::isfinite(first_x) ||
                !(first_x > 0.0)) {
                throw std::invalid_argument(
                    "CPA association derivative: invalid primal site fraction");
            }
            for (std::size_t b = 0U;
                 b < n;
                 ++b) {
                const auto& second =
                    association_->sites[b];
                double value = 0.0;
                const auto* pair =
                    parameters_->association_pair(
                        first.component_index,
                        first.site_id,
                        second.component_index,
                        second.site_id);
                if (pair != nullptr &&
                    pair->beta_dimensionless !=
                        0.0) {
                    const double b_ij =
                        0.5 *
                        (parameters_
                             ->pure(
                                 first.component_index)
                             .b_m3_per_mol +
                         parameters_
                             ->pure(
                                 second.component_index)
                             .b_m3_per_mol);
                    const double delta =
                        cpa_association_strength(
                            temperature_k_,
                            g,
                            b_ij,
                            *pair);
                    value -=
                        mole_numbers_[
                            second.component_index] *
                        static_cast<double>(
                            second.multiplicity) *
                        delta /
                        volume_m3_;
                }
                if (a == b) {
                    value -=
                        1.0 /
                        (first_x * first_x);
                }
                if (!std::isfinite(value)) {
                    throw std::range_error(
                        "CPA association derivative: non-finite site Jacobian");
                }
                lu_[a * n + b] = value;
            }
        }

        double matrix_scale = 0.0;
        for (double value : lu_) {
            matrix_scale =
                std::max(
                    matrix_scale,
                    std::abs(value));
        }
        const double pivot_guard =
            4096.0 *
            std::numeric_limits<double>::
                epsilon() *
            std::max(1.0, matrix_scale);

        for (std::size_t k = 0U;
             k < n;
             ++k) {
            std::size_t pivot = k;
            double pivot_abs =
                std::abs(
                    lu_[k * n + k]);
            for (std::size_t row =
                     k + 1U;
                 row < n;
                 ++row) {
                const double candidate =
                    std::abs(
                        lu_[row * n + k]);
                if (candidate >
                    pivot_abs) {
                    pivot = row;
                    pivot_abs = candidate;
                }
            }
            if (!std::isfinite(pivot_abs) ||
                !(pivot_abs > pivot_guard)) {
                throw std::runtime_error(
                    "CPA association derivative: singular or ill-conditioned site Jacobian");
            }
            pivot_rows_[k] = pivot;
            if (pivot != k) {
                for (std::size_t column =
                         0U;
                     column < n;
                     ++column) {
                    std::swap(
                        lu_[k * n + column],
                        lu_[pivot * n +
                            column]);
                }
            }

            const double diagonal =
                lu_[k * n + k];
            for (std::size_t row =
                     k + 1U;
                 row < n;
                 ++row) {
                double& factor =
                    lu_[row * n + k];
                factor /= diagonal;
                for (std::size_t column =
                         k + 1U;
                     column < n;
                     ++column) {
                    lu_[row * n + column] -=
                        factor *
                        lu_[k * n + column];
                }
            }
        }
    }

    [[nodiscard]] std::vector<double>
    solve_factored(
        std::vector<double> rhs) const {
        const std::size_t n =
            site_count();
        for (std::size_t k = 0U;
             k < n;
             ++k) {
            if (pivot_rows_[k] != k) {
                std::swap(
                    rhs[k],
                    rhs[pivot_rows_[k]]);
            }
        }
        for (std::size_t row = 0U;
             row < n;
             ++row) {
            for (std::size_t column = 0U;
                 column < row;
                 ++column) {
                rhs[row] -=
                    lu_[row * n + column] *
                    rhs[column];
            }
        }
        for (std::size_t reverse = n;
             reverse > 0U;
             --reverse) {
            const std::size_t row =
                reverse - 1U;
            for (std::size_t column =
                     row + 1U;
                 column < n;
                 ++column) {
                rhs[row] -=
                    lu_[row * n + column] *
                    rhs[column];
            }
            rhs[row] /=
                lu_[row * n + row];
            if (!std::isfinite(rhs[row])) {
                throw std::range_error(
                    "CPA association derivative: non-finite site sensitivity");
            }
        }
        return rhs;
    }

    double temperature_k_{};
    double volume_m3_{};
    std::vector<double> mole_numbers_;
    const CpaParameterSet* parameters_{};
    const CpaAssociationResult* association_{};
    std::vector<double> lu_;
    std::vector<std::size_t> pivot_rows_;
};

struct CpaSecondDirectional {
    double value{};
    double first{};
    double second{};
    double cross{};

    constexpr CpaSecondDirectional() =
        default;
    constexpr CpaSecondDirectional(
        double primal)
        : value(primal) {}
    constexpr CpaSecondDirectional(
        double primal,
        double d_first,
        double d_second,
        double d_cross = 0.0)
        : value(primal),
          first(d_first),
          second(d_second),
          cross(d_cross) {}
};

inline CpaSecondDirectional&
operator+=(
    CpaSecondDirectional& a,
    const CpaSecondDirectional& b) {
    a.value += b.value;
    a.first += b.first;
    a.second += b.second;
    a.cross += b.cross;
    return a;
}

inline CpaSecondDirectional&
operator-=(
    CpaSecondDirectional& a,
    const CpaSecondDirectional& b) {
    a.value -= b.value;
    a.first -= b.first;
    a.second -= b.second;
    a.cross -= b.cross;
    return a;
}

[[nodiscard]] inline CpaSecondDirectional
operator+(
    const CpaSecondDirectional& a,
    const CpaSecondDirectional& b) {
    return {
        a.value + b.value,
        a.first + b.first,
        a.second + b.second,
        a.cross + b.cross};
}

[[nodiscard]] inline CpaSecondDirectional
operator-(
    const CpaSecondDirectional& a,
    const CpaSecondDirectional& b) {
    return {
        a.value - b.value,
        a.first - b.first,
        a.second - b.second,
        a.cross - b.cross};
}

[[nodiscard]] inline CpaSecondDirectional
operator-(
    const CpaSecondDirectional& a) {
    return {
        -a.value,
        -a.first,
        -a.second,
        -a.cross};
}

[[nodiscard]] inline CpaSecondDirectional
operator*(
    const CpaSecondDirectional& a,
    const CpaSecondDirectional& b) {
    return {
        a.value * b.value,
        a.first * b.value +
            a.value * b.first,
        a.second * b.value +
            a.value * b.second,
        a.cross * b.value +
            a.first * b.second +
            a.second * b.first +
            a.value * b.cross};
}

[[nodiscard]] inline CpaSecondDirectional
reciprocal(
    const CpaSecondDirectional& a) {
    const double inverse =
        1.0 / a.value;
    const double first_derivative =
        -inverse * inverse;
    const double second_derivative =
        2.0 * inverse * inverse *
        inverse;
    return {
        inverse,
        first_derivative * a.first,
        first_derivative * a.second,
        second_derivative *
                a.first *
                a.second +
            first_derivative *
                a.cross};
}

[[nodiscard]] inline CpaSecondDirectional
operator/(
    const CpaSecondDirectional& a,
    const CpaSecondDirectional& b) {
    return a * reciprocal(b);
}

[[nodiscard]] inline CpaSecondDirectional
unary_second(
    const CpaSecondDirectional& a,
    double value,
    double slope,
    double curvature) {
    return {
        value,
        slope * a.first,
        slope * a.second,
        curvature *
                a.first *
                a.second +
            slope * a.cross};
}

[[nodiscard]] inline CpaSecondDirectional
sqrt(
    const CpaSecondDirectional& a) {
    const double value =
        std::sqrt(a.value);
    const double slope =
        0.5 / value;
    const double curvature =
        -0.25 /
        (a.value * value);
    return unary_second(
        a,
        value,
        slope,
        curvature);
}

[[nodiscard]] inline CpaSecondDirectional
log(
    const CpaSecondDirectional& a) {
    const double slope =
        1.0 / a.value;
    return unary_second(
        a,
        std::log(a.value),
        slope,
        -slope * slope);
}

[[nodiscard]] inline CpaSecondDirectional
log1p(
    const CpaSecondDirectional& a) {
    const double denominator =
        1.0 + a.value;
    const double slope =
        1.0 / denominator;
    return unary_second(
        a,
        std::log1p(a.value),
        slope,
        -slope * slope);
}

[[nodiscard]] inline CpaSecondDirectional
expm1(
    const CpaSecondDirectional& a) {
    const double exponential =
        std::exp(a.value);
    return unary_second(
        a,
        std::expm1(a.value),
        exponential,
        exponential);
}

[[nodiscard]] inline double
cpa_residual_helmholtz_cross_derivative(
    double temperature_k,
    double volume_m3,
    std::span<const double> mole_numbers,
    const CpaParameterSet& parameters,
    const CpaAssociationResult& association,
    const CpaStateDirection& first,
    const CpaStateDirection& second,
    std::span<const double> first_site_sensitivity,
    std::span<const double> second_site_sensitivity) {
    if (first.mole_numbers.size() !=
            mole_numbers.size() ||
        second.mole_numbers.size() !=
            mole_numbers.size() ||
        first_site_sensitivity.size() !=
            association.sites.size() ||
        second_site_sensitivity.size() !=
            association.sites.size()) {
        throw std::invalid_argument(
            "CPA second derivative: direction dimension mismatch");
    }

    using H = CpaSecondDirectional;
    const H temperature{
        temperature_k,
        first.temperature,
        second.temperature};
    const H volume{
        volume_m3,
        first.volume,
        second.volume};

    std::vector<H> moles;
    moles.reserve(mole_numbers.size());
    for (std::size_t i = 0U;
         i < mole_numbers.size();
         ++i) {
        moles.emplace_back(
            mole_numbers[i],
            first.mole_numbers[i],
            second.mole_numbers[i]);
    }

    std::vector<H> sites;
    sites.reserve(
        association.sites.size());
    for (std::size_t i = 0U;
         i < association.sites.size();
         ++i) {
        sites.emplace_back(
            association.sites[i]
                .unbonded_fraction,
            first_site_sensitivity[i],
            second_site_sensitivity[i]);
    }

    const H cubic =
        cpa_cubic_residual_helmholtz_reduced(
            temperature,
            volume,
            std::span<const H>{moles},
            parameters);
    const H association_term =
        cpa_association_q_reduced_with_site_fractions(
            temperature,
            volume,
            std::span<const H>{moles},
            parameters,
            association,
            std::span<const H>{sites});
    const double result =
        cubic.cross +
        association_term.cross;
    if (!std::isfinite(result)) {
        throw std::range_error(
            "CPA second derivative: non-finite Helmholtz cross derivative");
    }
    return result;
}

[[nodiscard]] inline double
cpa_pressure_directional_derivative(
    double temperature_k,
    double volume_m3,
    std::span<const double> mole_numbers,
    const CpaParameterSet& parameters,
    const CpaPhaseState& phase_state,
    const CpaAssociationSensitivitySystem&
        association_system,
    const CpaStateDirection& direction,
    std::span<const double> volume_site_sensitivity,
    std::span<const double> direction_site_sensitivity) {
    double total_moles = 0.0;
    double total_mole_direction = 0.0;
    for (std::size_t i = 0U;
         i < mole_numbers.size();
         ++i) {
        total_moles +=
            mole_numbers[i];
        total_mole_direction +=
            direction.mole_numbers[i];
    }

    const double rt =
        cpa_gas_constant_j_per_mol_k *
        temperature_k;
    const double ideal_pressure =
        total_moles * rt /
        volume_m3;
    const double f_volume =
        (ideal_pressure -
         phase_state.pressure_pa) /
        rt;

    CpaStateDirection volume_direction;
    volume_direction.volume = 1.0;
    volume_direction.mole_numbers.assign(
        mole_numbers.size(),
        0.0);

    const double f_volume_direction =
        cpa_residual_helmholtz_cross_derivative(
            temperature_k,
            volume_m3,
            mole_numbers,
            parameters,
            phase_state.association,
            volume_direction,
            direction,
            volume_site_sensitivity,
            direction_site_sensitivity);

    const double ideal_direction =
        cpa_gas_constant_j_per_mol_k *
        ((total_mole_direction *
              temperature_k +
          total_moles *
              direction.temperature) /
             volume_m3 -
         total_moles *
             temperature_k *
             direction.volume /
             (volume_m3 *
              volume_m3));

    const double derivative =
        ideal_direction -
        cpa_gas_constant_j_per_mol_k *
            (direction.temperature *
                 f_volume +
             temperature_k *
                 f_volume_direction);
    if (!std::isfinite(derivative)) {
        throw std::range_error(
            "CPA PT derivative: non-finite pressure derivative");
    }
    (void)association_system;
    return derivative;
}

} // namespace cpa_detail

template <std::size_t K>
[[nodiscard]] inline std::vector<
    mpmc::ad::Dual<double, K>>
evaluate_cpa_selected_pt_ln_phi_first_order(
    const CpaPtPhase& model,
    const mpmc::ad::Dual<double, K>&
        pressure_pa,
    const mpmc::ad::Dual<double, K>&
        temperature_k,
    std::span<const mpmc::ad::Dual<double, K>>
        composition,
    std::size_t root_index,
    CpaPtOptions options = {}) {
    using D =
        mpmc::ad::Dual<double, K>;

    std::vector<double> primal_composition;
    primal_composition.reserve(
        composition.size());
    for (const auto& value :
         composition) {
        primal_composition.push_back(
            value.value());
    }
    (void)cpa_detail::
        validate_cpa_composition(
            primal_composition,
            model.size());

    for (std::size_t lane = 0U;
         lane < K;
         ++lane) {
        double sum = 0.0;
        double scale = 0.0;
        for (const auto& value :
             composition) {
            const double derivative =
                value.derivative(lane);
            sum += derivative;
            scale +=
                std::abs(derivative);
        }
        const double tolerance =
            4096.0 *
            std::numeric_limits<double>::
                epsilon() *
            std::max(1.0, scale);
        if (!std::isfinite(sum) ||
            std::abs(sum) >
                tolerance) {
            throw std::invalid_argument(
                "CPA selected phase derivative: composition direction must remain tangent to the mole-fraction simplex");
        }
    }

    const double p =
        pressure_pa.value();
    const double t =
        temperature_k.value();
    const auto roots =
        model.roots(
            p,
            t,
            primal_composition,
            options);
    if (roots.status !=
        CpaPtRootStatus::success) {
        if (roots.status ==
            CpaPtRootStatus::
                near_multiple) {
            throw std::domain_error(
                "CPA selected phase derivative: near-multiple/tangent PT root is not differentiable");
        }
        throw std::runtime_error(
            "CPA selected phase derivative: PT root set unavailable: " +
            roots.diagnostic);
    }
    if (root_index >=
        roots.roots.size()) {
        throw std::out_of_range(
            "CPA selected phase derivative: selected root index out of range");
    }

    const auto& root =
        roots.roots[root_index];
    const double rho =
        root.molar_density_mol_per_m3;
    const double volume =
        1.0 / rho;

    const auto phase_state =
        evaluate_cpa_phase_at_density(
            t,
            rho,
            primal_composition,
            model.parameters(),
            options.phase);

    const auto mole_span =
        std::span<const double>{
            primal_composition};
    cpa_detail::
        CpaAssociationSensitivitySystem
            association_system{
                t,
                volume,
                mole_span,
                model.parameters(),
                phase_state.association};

    cpa_detail::CpaStateDirection
        volume_direction;
    volume_direction.volume = 1.0;
    volume_direction.mole_numbers.assign(
        primal_composition.size(),
        0.0);
    const auto x_volume =
        association_system.solve_direction(
            volume_direction);

    const double pressure_volume =
        cpa_detail::
            cpa_pressure_directional_derivative(
                t,
                volume,
                mole_span,
                model.parameters(),
                phase_state,
                association_system,
                volume_direction,
                x_volume,
                x_volume);

    const double pressure_scale =
        std::max(
            1.0,
            std::abs(p / volume));
    const double pressure_guard =
        4096.0 *
        std::numeric_limits<double>::
            epsilon() *
        pressure_scale;
    if (!std::isfinite(
            pressure_volume) ||
        std::abs(pressure_volume) <=
            pressure_guard) {
        throw std::runtime_error(
            "CPA selected phase derivative: ill-conditioned density root");
    }

    std::vector<
        cpa_detail::CpaStateDirection>
        component_directions;
    std::vector<std::vector<double>>
        x_components;
    component_directions.reserve(
        primal_composition.size());
    x_components.reserve(
        primal_composition.size());
    for (std::size_t component = 0U;
         component <
         primal_composition.size();
         ++component) {
        cpa_detail::CpaStateDirection
            direction;
        direction.mole_numbers.assign(
            primal_composition.size(),
            0.0);
        direction.mole_numbers[
            component] = 1.0;
        x_components.push_back(
            association_system
                .solve_direction(
                    direction));
        component_directions.push_back(
            std::move(direction));
    }

    std::vector<
        typename D::Gradient>
        gradients(
            primal_composition.size());

    for (std::size_t lane = 0U;
         lane < K;
         ++lane) {
        cpa_detail::CpaStateDirection
            no_volume;
        no_volume.temperature =
            temperature_k.derivative(
                lane);
        no_volume.mole_numbers.resize(
            primal_composition.size());
        for (std::size_t component =
                 0U;
             component <
             primal_composition.size();
             ++component) {
            no_volume.mole_numbers[
                component] =
                composition[component]
                    .derivative(lane);
        }
        const auto x_no_volume =
            association_system
                .solve_direction(
                    no_volume);
        const double pressure_no_volume =
            cpa_detail::
                cpa_pressure_directional_derivative(
                    t,
                    volume,
                    mole_span,
                    model.parameters(),
                    phase_state,
                    association_system,
                    no_volume,
                    x_volume,
                    x_no_volume);

        const double d_volume =
            (pressure_pa.derivative(
                 lane) -
             pressure_no_volume) /
            pressure_volume;

        auto total_direction =
            no_volume;
        total_direction.volume =
            d_volume;

        std::vector<double>
            x_total(
                x_no_volume.size(),
                0.0);
        for (std::size_t site = 0U;
             site < x_total.size();
             ++site) {
            x_total[site] =
                x_no_volume[site] +
                d_volume *
                    x_volume[site];
        }

        const double d_ln_z =
            pressure_pa.derivative(
                lane) /
                p +
            d_volume /
                volume -
            temperature_k.derivative(
                lane) /
                t;

        for (std::size_t component =
                 0U;
             component <
             primal_composition.size();
             ++component) {
            const double d_mu_res =
                cpa_detail::
                    cpa_residual_helmholtz_cross_derivative(
                        t,
                        volume,
                        mole_span,
                        model.parameters(),
                        phase_state.association,
                        component_directions[
                            component],
                        total_direction,
                        x_components[
                            component],
                        x_total);
            gradients[component][lane] =
                d_mu_res -
                d_ln_z;
            if (!std::isfinite(
                    gradients[component]
                             [lane])) {
                throw std::range_error(
                    "CPA selected phase derivative: non-finite ln(phi) derivative");
            }
        }
    }

    std::vector<D> result;
    result.reserve(
        primal_composition.size());
    for (std::size_t component = 0U;
         component <
         primal_composition.size();
         ++component) {
        result.emplace_back(
            root.ln_phi[component],
            gradients[component]);
    }
    return result;
}

} // namespace mpmc::thermodynamics

#endif // MPMC_THERMODYNAMICS_CPA_PT_PHASE_AD_HPP
