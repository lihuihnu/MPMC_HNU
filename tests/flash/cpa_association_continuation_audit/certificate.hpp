#ifndef MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_HPP
#define MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_HPP

#include <mpmc/thermodynamics/cpa_association.hpp>

#include "contract.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace cpa_association_observable_certificate {

namespace th = mpmc::thermodynamics;
namespace contract = cpa_association_continuation_contract;

using Real = long double;
using Matrix = std::vector<Real>;

struct Result {
    bool enclosure_certified{};
    bool pressure_scale_certified{};
    bool ln_phi_scale_certified{};
    double fresh_fixed_point_residual_inf{};
    double max_site_error_bound{};
    double pressure_error_bound_pa{};
    double max_ln_phi_error_bound{};
    double krawczyk_max_inclusion_ratio{};
    std::vector<double> site_error_bounds;
    std::vector<double> ln_phi_error_bounds;

    [[nodiscard]] bool observable_certified() const noexcept {
        return enclosure_certified &&
               pressure_scale_certified &&
               ln_phi_scale_certified;
    }
};

inline std::size_t index(std::size_t n, std::size_t row, std::size_t col) {
    return row * n + col;
}

inline bool finite(Real value) {
    return std::isfinite(value);
}

inline bool invert(Matrix matrix, std::size_t n, Matrix& inverse) {
    inverse.assign(n * n, 0.0L);
    for (std::size_t i = 0U; i < n; ++i) {
        inverse[index(n, i, i)] = 1.0L;
    }

    constexpr Real pivot_floor = 1024.0L *
        std::numeric_limits<Real>::epsilon();
    for (std::size_t column = 0U; column < n; ++column) {
        std::size_t pivot_row = column;
        Real pivot_magnitude = std::abs(matrix[index(n, column, column)]);
        for (std::size_t row = column + 1U; row < n; ++row) {
            const Real candidate = std::abs(matrix[index(n, row, column)]);
            if (candidate > pivot_magnitude) {
                pivot_row = row;
                pivot_magnitude = candidate;
            }
        }
        if (!finite(pivot_magnitude) || !(pivot_magnitude > pivot_floor)) {
            return false;
        }
        if (pivot_row != column) {
            for (std::size_t j = 0U; j < n; ++j) {
                std::swap(matrix[index(n, column, j)],
                          matrix[index(n, pivot_row, j)]);
                std::swap(inverse[index(n, column, j)],
                          inverse[index(n, pivot_row, j)]);
            }
        }

        const Real pivot = matrix[index(n, column, column)];
        for (std::size_t j = 0U; j < n; ++j) {
            matrix[index(n, column, j)] /= pivot;
            inverse[index(n, column, j)] /= pivot;
        }
        for (std::size_t row = 0U; row < n; ++row) {
            if (row == column) { continue; }
            const Real factor = matrix[index(n, row, column)];
            if (factor == 0.0L) { continue; }
            for (std::size_t j = 0U; j < n; ++j) {
                matrix[index(n, row, j)] -=
                    factor * matrix[index(n, column, j)];
                inverse[index(n, row, j)] -=
                    factor * inverse[index(n, column, j)];
            }
        }
    }
    for (const Real value : inverse) {
        if (!finite(value)) { return false; }
    }
    return true;
}

inline std::vector<Real> fixed_point_map(
    std::span<const Real> candidate,
    const Matrix& a_matrix,
    std::size_t n) {
    std::vector<Real> mapped(n, 0.0L);
    for (std::size_t a = 0U; a < n; ++a) {
        Real denominator = 1.0L;
        for (std::size_t b = 0U; b < n; ++b) {
            denominator += a_matrix[index(n, a, b)] * candidate[b];
        }
        if (!finite(denominator) || !(denominator > 0.0L)) {
            return {};
        }
        mapped[a] = 1.0L / denominator;
    }
    return mapped;
}

inline Result certify(
    double temperature_k,
    double molar_density_mol_per_m3,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaAssociationResult& association) {
    Result result;
    if (!association.converged()) { return result; }
    if (association.status == th::CpaAssociationStatus::no_associating_sites) {
        result.enclosure_certified = true;
        result.pressure_scale_certified = true;
        result.ln_phi_scale_certified = true;
        result.ln_phi_error_bounds.assign(parameters.size(), 0.0);
        return result;
    }

    const std::size_t n = association.sites.size();
    if (n == 0U || composition.size() != parameters.size()) { return result; }

    std::vector<Real> candidate(n, 0.0L);
    std::vector<Real> site_weight(n, 0.0L);
    Matrix a_matrix(n * n, 0.0L);
    for (std::size_t a = 0U; a < n; ++a) {
        const auto& first = association.sites[a];
        if (first.component_index >= parameters.size() ||
            !(first.unbonded_fraction > 0.0) ||
            !std::isfinite(first.unbonded_fraction)) {
            return result;
        }
        candidate[a] = static_cast<Real>(first.unbonded_fraction);
        site_weight[a] =
            static_cast<Real>(composition[first.component_index]) *
            static_cast<Real>(first.multiplicity);
    }

    const Real temperature = static_cast<Real>(temperature_k);
    const Real rho = static_cast<Real>(molar_density_mol_per_m3);
    const Real gas_constant = static_cast<Real>(th::cpa_gas_constant_j_per_mol_k);
    if (!finite(temperature) || !(temperature > 0.0L) ||
        !finite(rho) || !(rho > 0.0L)) {
        return result;
    }

    for (std::size_t a = 0U; a < n; ++a) {
        const auto& first = association.sites[a];
        for (std::size_t b = 0U; b < n; ++b) {
            const auto& second = association.sites[b];
            const auto* pair = parameters.association_pair(
                first.component_index, first.site_id,
                second.component_index, second.site_id);
            if (pair == nullptr || pair->beta_dimensionless == 0.0 ||
                composition[second.component_index] == 0.0) {
                continue;
            }
            const Real b_ij = 0.5L *
                (static_cast<Real>(parameters.pure(first.component_index).b_m3_per_mol) +
                 static_cast<Real>(parameters.pure(second.component_index).b_m3_per_mol));
            const Real exponent = static_cast<Real>(pair->epsilon_j_per_mol) /
                (gas_constant * temperature);
            if (!finite(exponent) || exponent > 700.0L) { return result; }
            const Real delta =
                static_cast<Real>(association.radial_distribution) *
                std::expm1(exponent) * b_ij *
                static_cast<Real>(pair->beta_dimensionless);
            if (!finite(delta) || delta < 0.0L) { return result; }
            a_matrix[index(n, a, b)] = rho * site_weight[b] * delta;
        }
    }

    const auto mapped = fixed_point_map(candidate, a_matrix, n);
    if (mapped.size() != n) { return result; }

    std::vector<Real> h(n, 0.0L);
    Real residual_inf = 0.0L;
    for (std::size_t i = 0U; i < n; ++i) {
        h[i] = candidate[i] - mapped[i];
        residual_inf = std::max(residual_inf, std::abs(h[i]));
    }
    result.fresh_fixed_point_residual_inf = static_cast<double>(residual_inf);

    Matrix jacobian(n * n, 0.0L);
    for (std::size_t a = 0U; a < n; ++a) {
        Real denominator = 1.0L;
        for (std::size_t b = 0U; b < n; ++b) {
            denominator += a_matrix[index(n, a, b)] * candidate[b];
        }
        if (!finite(denominator) || !(denominator > 0.0L)) { return result; }
        const Real denominator_squared = denominator * denominator;
        for (std::size_t b = 0U; b < n; ++b) {
            jacobian[index(n, a, b)] =
                (a == b ? 1.0L : 0.0L) +
                a_matrix[index(n, a, b)] / denominator_squared;
        }
    }

    Matrix inverse;
    if (!invert(jacobian, n, inverse)) { return result; }

    std::vector<Real> correction(n, 0.0L);
    for (std::size_t i = 0U; i < n; ++i) {
        Real value = 0.0L;
        for (std::size_t j = 0U; j < n; ++j) {
            value -= inverse[index(n, i, j)] * h[j];
        }
        correction[i] = value;
    }

    const Real epsilon = std::numeric_limits<Real>::epsilon();
    std::vector<Real> radius(n, 0.0L);
    for (std::size_t i = 0U; i < n; ++i) {
        const Real floor = 256.0L * epsilon *
            std::max(1.0L, std::abs(candidate[i]));
        radius[i] = std::max(4.0L * std::abs(correction[i]), floor);
    }

    Real final_ratio = std::numeric_limits<Real>::infinity();
    constexpr std::size_t max_enclosure_iterations = 48U;
    for (std::size_t iteration = 0U;
         iteration < max_enclosure_iterations; ++iteration) {
        std::vector<Real> lower(n, 0.0L);
        std::vector<Real> upper(n, 0.0L);
        bool positive_box = true;
        for (std::size_t i = 0U; i < n; ++i) {
            lower[i] = candidate[i] - radius[i];
            upper[i] = candidate[i] + radius[i];
            if (!(lower[i] > 0.0L) || !finite(upper[i])) {
                positive_box = false;
            }
        }
        if (!positive_box) { return result; }

        Matrix midpoint(n * n, 0.0L);
        Matrix jacobian_radius(n * n, 0.0L);
        for (std::size_t a = 0U; a < n; ++a) {
            Real denominator_min = 1.0L;
            Real denominator_max = 1.0L;
            for (std::size_t b = 0U; b < n; ++b) {
                const Real coefficient = a_matrix[index(n, a, b)];
                denominator_min += coefficient * lower[b];
                denominator_max += coefficient * upper[b];
            }
            if (!(denominator_min > 0.0L) ||
                !(denominator_max >= denominator_min)) {
                return result;
            }
            const Real inverse_min_squared =
                1.0L / (denominator_min * denominator_min);
            const Real inverse_max_squared =
                1.0L / (denominator_max * denominator_max);
            for (std::size_t b = 0U; b < n; ++b) {
                const Real diagonal = a == b ? 1.0L : 0.0L;
                const Real low = diagonal +
                    a_matrix[index(n, a, b)] * inverse_max_squared;
                const Real high = diagonal +
                    a_matrix[index(n, a, b)] * inverse_min_squared;
                midpoint[index(n, a, b)] = 0.5L * (low + high);
                jacobian_radius[index(n, a, b)] = 0.5L * (high - low);
            }
        }

        Matrix multiplier(n * n, 0.0L);
        for (std::size_t i = 0U; i < n; ++i) {
            for (std::size_t j = 0U; j < n; ++j) {
                Real product = 0.0L;
                Real uncertainty = 0.0L;
                for (std::size_t k = 0U; k < n; ++k) {
                    product += inverse[index(n, i, k)] *
                        midpoint[index(n, k, j)];
                    uncertainty += std::abs(inverse[index(n, i, k)]) *
                        jacobian_radius[index(n, k, j)];
                }
                const Real center = (i == j ? 1.0L : 0.0L) - product;
                multiplier[index(n, i, j)] = std::abs(center) + uncertainty;
            }
        }

        std::vector<Real> required(n, 0.0L);
        bool included = true;
        final_ratio = 0.0L;
        for (std::size_t i = 0U; i < n; ++i) {
            Real value = std::abs(correction[i]);
            for (std::size_t j = 0U; j < n; ++j) {
                value += multiplier[index(n, i, j)] * radius[j];
            }
            const Real roundoff_pad = 1024.0L * epsilon *
                std::max(1.0L, std::abs(candidate[i]) + radius[i]);
            value += roundoff_pad;
            required[i] = value;
            final_ratio = std::max(final_ratio, value / radius[i]);
            if (!(value < radius[i])) { included = false; }
        }

        if (included) {
            result.enclosure_certified = true;
            break;
        }

        for (std::size_t i = 0U; i < n; ++i) {
            radius[i] = std::max(1.5L * radius[i], 1.10L * required[i]);
        }
    }

    result.krawczyk_max_inclusion_ratio = static_cast<double>(final_ratio);
    if (!result.enclosure_certified) { return result; }

    result.site_error_bounds.resize(n, 0.0);
    Real max_site_error = 0.0L;
    Real weighted_radius = 0.0L;
    for (std::size_t i = 0U; i < n; ++i) {
        result.site_error_bounds[i] = static_cast<double>(radius[i]);
        max_site_error = std::max(max_site_error, radius[i]);
        weighted_radius += site_weight[i] * radius[i];
    }
    result.max_site_error_bound = static_cast<double>(max_site_error);

    const Real pressure_coefficient = 0.5L * gas_constant * temperature * rho *
        (1.0L + static_cast<Real>(association.rho_dln_g_drho));
    const Real pressure_bound = pressure_coefficient * weighted_radius;
    if (!finite(pressure_bound) || pressure_bound < 0.0L) { return Result{}; }
    result.pressure_error_bound_pa = static_cast<double>(pressure_bound);
    result.pressure_scale_certified =
        pressure_bound <= static_cast<Real>(contract::pressure_error_guard_pa);

    result.ln_phi_error_bounds.assign(parameters.size(), 0.0);
    Real max_ln_phi_bound = 0.0L;
    for (std::size_t component = 0U;
         component < parameters.size(); ++component) {
        Real log_bound = 0.0L;
        for (std::size_t site = 0U; site < n; ++site) {
            const auto& site_state = association.sites[site];
            if (site_state.component_index != component) { continue; }
            const Real lower = candidate[site] - radius[site];
            const Real upper = candidate[site] + radius[site];
            if (!(lower > 0.0L) || !(upper > 0.0L)) { return Result{}; }
            const Real downward = std::log(candidate[site] / lower);
            const Real upward = std::log(upper / candidate[site]);
            log_bound += static_cast<Real>(site_state.multiplicity) *
                std::max(downward, upward);
        }
        const Real sensitivity = (1.9L / 8.0L) * rho *
            static_cast<Real>(parameters.pure(component).b_m3_per_mol) *
            static_cast<Real>(association.radial_distribution);
        const Real bound = log_bound + sensitivity * weighted_radius;
        if (!finite(bound) || bound < 0.0L) { return Result{}; }
        result.ln_phi_error_bounds[component] = static_cast<double>(bound);
        max_ln_phi_bound = std::max(max_ln_phi_bound, bound);
    }
    result.max_ln_phi_error_bound = static_cast<double>(max_ln_phi_bound);
    result.ln_phi_scale_certified =
        max_ln_phi_bound <= static_cast<Real>(contract::ln_phi_error_guard);
    return result;
}

} // namespace cpa_association_observable_certificate

#endif // MPMC_TEST_CPA_ASSOCIATION_OBSERVABLE_CERTIFICATE_HPP
