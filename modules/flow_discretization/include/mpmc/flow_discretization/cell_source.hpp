#ifndef MPMC_FLOW_DISCRETIZATION_CELL_SOURCE_HPP
#define MPMC_FLOW_DISCRETIZATION_CELL_SOURCE_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::flow_discretization {

inline constexpr std::string_view
    cell_source_convention =
        "flow_discretization/cell-source-linearization/v1";

/// Model-neutral linearization of one external source attached to one finite
/// volume cell.
///
/// Positive component/energy rates inject material/energy INTO the control
/// volume. With residual convention accumulation + outward flux - injection,
/// source residual contributions are -q_i/V_b and -Q_E/V_b.
struct CellSourceLinearization3D {
    static constexpr std::string_view convention =
        cell_source_convention;

    std::string provenance;
    std::vector<std::string> component_ids;
    std::size_t input_count{};
    std::vector<double>
        component_molar_rate_mol_per_s;
    std::vector<double>
        component_molar_rate_jacobian_mol_per_s;
    double energy_rate_w{};
    std::vector<double>
        energy_rate_gradient_w;

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_ids.size();
    }

    [[nodiscard]] double
    d_component_rate(
        std::size_t component,
        std::size_t column) const {
        if (component >= component_count() ||
            column >= input_count ||
            input_count == 0U ||
            component >
                (std::numeric_limits<std::size_t>::max() -
                 column) /
                    input_count) {
            throw std::out_of_range(
                "mpmc::flow_discretization: cell source Jacobian index out of range");
        }
        return component_molar_rate_jacobian_mol_per_s.at(
            component * input_count + column);
    }

    [[nodiscard]] double
    d_energy_rate(
        std::size_t column) const {
        return energy_rate_gradient_w.at(column);
    }
};

struct NormalizedCellSourceContribution3D {
    static constexpr std::string_view convention =
        cell_source_convention;

    std::string provenance;
    std::vector<std::string> component_ids;
    std::size_t input_count{};
    double bulk_volume_m3{};
    std::vector<double>
        component_residual_mol_per_bulk_m3_s;
    std::vector<double>
        component_jacobian_mol_per_bulk_m3_s;
    double energy_residual_w_per_bulk_m3{};
    std::vector<double>
        energy_gradient_w_per_bulk_m3;

    [[nodiscard]] double
    d_component_residual(
        std::size_t component,
        std::size_t column) const {
        return component_jacobian_mol_per_bulk_m3_s.at(
            component * input_count + column);
    }

    [[nodiscard]] double
    d_energy_residual(
        std::size_t column) const {
        return energy_gradient_w_per_bulk_m3.at(column);
    }
};

inline void validate_cell_source_linearization(
    const CellSourceLinearization3D& source,
    std::span<const std::string> expected_component_ids,
    std::size_t expected_input_count) {
    const std::size_t n =
        expected_component_ids.size();
    if (source.provenance.empty() ||
        n == 0U ||
        expected_input_count == 0U ||
        source.component_ids.size() != n ||
        !std::equal(
            source.component_ids.begin(),
            source.component_ids.end(),
            expected_component_ids.begin(),
            expected_component_ids.end()) ||
        source.input_count != expected_input_count ||
        n >
            std::numeric_limits<std::size_t>::max() /
                expected_input_count ||
        source.component_molar_rate_mol_per_s.size() != n ||
        source.component_molar_rate_jacobian_mol_per_s.size() !=
            n * expected_input_count ||
        source.energy_rate_gradient_w.size() !=
            expected_input_count ||
        !std::isfinite(source.energy_rate_w)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: malformed cell source linearization");
    }

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        if (source.component_ids[component].empty()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: cell source component identity must be nonempty");
        }
        for (std::size_t previous = 0U;
             previous < component;
             ++previous) {
            if (source.component_ids[previous] ==
                source.component_ids[component]) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: cell source component identity must be unique and ordered");
            }
        }
        if (!std::isfinite(
                source.component_molar_rate_mol_per_s[
                    component])) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: cell source component rate must be finite");
        }
    }

    for (double value :
         source.component_molar_rate_jacobian_mol_per_s) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: cell source component derivative must be finite");
        }
    }
    for (double value :
         source.energy_rate_gradient_w) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: cell source energy derivative must be finite");
        }
    }
}

[[nodiscard]] inline
NormalizedCellSourceContribution3D
normalize_cell_source_by_bulk_volume(
    const CellSourceLinearization3D& source,
    double bulk_volume_m3) {
    validate_cell_source_linearization(
        source,
        source.component_ids,
        source.input_count);
    if (!std::isfinite(bulk_volume_m3) ||
        !(bulk_volume_m3 > 0.0)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: cell source bulk volume must be finite and strictly positive [m^3]");
    }

    const double inverse_volume =
        1.0 / bulk_volume_m3;
    NormalizedCellSourceContribution3D result;
    result.provenance = source.provenance;
    result.component_ids = source.component_ids;
    result.input_count = source.input_count;
    result.bulk_volume_m3 = bulk_volume_m3;
    result.component_residual_mol_per_bulk_m3_s.resize(
        source.component_count());
    result.component_jacobian_mol_per_bulk_m3_s.resize(
        source.component_molar_rate_jacobian_mol_per_s.size());
    result.energy_gradient_w_per_bulk_m3.resize(
        source.input_count);

    for (std::size_t component = 0U;
         component < source.component_count();
         ++component) {
        result.component_residual_mol_per_bulk_m3_s[
            component] =
            -source.component_molar_rate_mol_per_s[
                 component] *
            inverse_volume;
    }
    for (std::size_t index = 0U;
         index <
         source.component_molar_rate_jacobian_mol_per_s.size();
         ++index) {
        result.component_jacobian_mol_per_bulk_m3_s[
            index] =
            -source.component_molar_rate_jacobian_mol_per_s[
                 index] *
            inverse_volume;
    }
    result.energy_residual_w_per_bulk_m3 =
        -source.energy_rate_w * inverse_volume;
    for (std::size_t column = 0U;
         column < source.input_count;
         ++column) {
        result.energy_gradient_w_per_bulk_m3[
            column] =
            -source.energy_rate_gradient_w[column] *
            inverse_volume;
    }
    return result;
}

} // namespace mpmc::flow_discretization

#endif // MPMC_FLOW_DISCRETIZATION_CELL_SOURCE_HPP
