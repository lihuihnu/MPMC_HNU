#ifndef MPMC_FLOW_FUGACITY_EQUILIBRIUM_LINEARIZATION_HPP
#define MPMC_FLOW_FUGACITY_EQUILIBRIUM_LINEARIZATION_HPP

#include <mpmc/flow/fugacity_equilibrium_residual.hpp>

#include <cmath>
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
    fugacity_equilibrium_linearization_convention =
        "flow/fugacity-equilibrium-residual-linearization/fixed-three-phase/v1";

/// Frozen local 2*Nc fugacity-equilibrium residual and its already-computed
/// Jacobian with respect to one frozen NaturalVariableLayout3P chart.
///
/// This type does not evaluate an EOS and does not choose a phase/root/family.
/// The caller supplies derivatives produced by the existing scalar-generic
/// fugacity residual path (for example through AD). The carrier validates only
/// identity, shape and finiteness before global assembly consumes it.
class FugacityEquilibriumResidualLinearization3P {
public:
    static constexpr std::string_view convention =
        fugacity_equilibrium_linearization_convention;

    FugacityEquilibriumResidualLinearization3P(
        NaturalVariableLayout3P layout,
        std::vector<std::string> component_ids,
        FugacityEquilibriumResidual3P<double> residual,
        std::size_t input_count,
        std::vector<double> jacobian)
        : layout_(std::move(layout)),
          component_ids_(std::move(component_ids)),
          residual_(std::move(residual)),
          input_count_(input_count),
          jacobian_(std::move(jacobian)) {
        validate();
    }

    [[nodiscard]] const NaturalVariableLayout3P&
    layout() const noexcept {
        return layout_;
    }

    [[nodiscard]] std::span<const std::string>
    component_ids() const noexcept {
        return component_ids_;
    }

    [[nodiscard]] std::size_t
    component_count() const noexcept {
        return component_ids_.size();
    }

    [[nodiscard]] std::size_t
    input_count() const noexcept {
        return input_count_;
    }

    [[nodiscard]] std::size_t
    residual_count() const noexcept {
        return residual_.residual_count();
    }

    [[nodiscard]] std::span<const double>
    values() const noexcept {
        return residual_.values();
    }

    [[nodiscard]] std::span<const double>
    jacobian() const noexcept {
        return jacobian_;
    }

    [[nodiscard]] std::size_t
    local_row_index(
        PhaseSlot3 non_reference_phase,
        std::size_t component) const {
        return residual_.local_row_index(
            non_reference_phase,
            component);
    }

    [[nodiscard]] std::size_t
    equation_index(
        PhaseSlot3 non_reference_phase,
        std::size_t component) const {
        return layout_.
            fugacity_equilibrium_equation_index(
                non_reference_phase,
                component);
    }

    [[nodiscard]] double residual(
        PhaseSlot3 non_reference_phase,
        std::size_t component) const {
        return residual_.residual(
            non_reference_phase,
            component);
    }

    [[nodiscard]] double d_residual(
        PhaseSlot3 non_reference_phase,
        std::size_t component,
        std::size_t column) const {
        if (column >= input_count_) {
            throw std::out_of_range(
                "mpmc::flow: fugacity-equilibrium Jacobian column out of range");
        }
        const std::size_t row =
            local_row_index(
                non_reference_phase,
                component);
        if (input_count_ == 0U ||
            row >
                (std::numeric_limits<std::size_t>::max() -
                 column) /
                    input_count_) {
            throw std::out_of_range(
                "mpmc::flow: fugacity-equilibrium Jacobian index out of range");
        }
        return jacobian_.at(
            row * input_count_ +
            column);
    }

private:
    void validate() const {
        const std::size_t n =
            layout_.component_count();
        const std::size_t q =
            layout_.unknown_count();

        if (n < 2U ||
            component_ids_.size() != n ||
            residual_.component_count() != n ||
            input_count_ != q ||
            residual_.residual_count() !=
                2U * n ||
            q == 0U ||
            residual_.residual_count() >
                std::numeric_limits<std::size_t>::max() /
                    q ||
            jacobian_.size() !=
                residual_.residual_count() * q) {
            throw std::invalid_argument(
                "mpmc::flow: fugacity-equilibrium linearization shape/layout mismatch");
        }

        for (std::size_t component = 0U;
             component < n;
             ++component) {
            if (component_ids_[component].empty()) {
                throw std::invalid_argument(
                    "mpmc::flow: fugacity-equilibrium component identity must be nonempty");
            }
            for (std::size_t previous = 0U;
                 previous < component;
                 ++previous) {
                if (component_ids_[previous] ==
                    component_ids_[component]) {
                    throw std::invalid_argument(
                        "mpmc::flow: fugacity-equilibrium component identity must be unique and ordered");
                }
            }
        }

        for (double value : residual_.values()) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "mpmc::flow: fugacity-equilibrium residual contains non-finite value");
            }
        }
        for (double derivative : jacobian_) {
            if (!std::isfinite(derivative)) {
                throw std::invalid_argument(
                    "mpmc::flow: fugacity-equilibrium Jacobian contains non-finite derivative");
            }
        }

        for (std::size_t phase = 1U;
             phase < fixed_three_phase_count;
             ++phase) {
            const auto slot =
                static_cast<PhaseSlot3>(
                    phase);
            for (std::size_t component = 0U;
                 component < n;
                 ++component) {
                const std::size_t local =
                    residual_.local_row_index(
                        slot,
                        component);
                const std::size_t equation =
                    layout_.
                        fugacity_equilibrium_equation_index(
                            slot,
                            component);
                if (equation !=
                    n + 1U + local) {
                    throw std::invalid_argument(
                        "mpmc::flow: fugacity-equilibrium local/global row ordering mismatch");
                }
            }
        }
    }

    NaturalVariableLayout3P layout_;
    std::vector<std::string> component_ids_;
    FugacityEquilibriumResidual3P<double>
        residual_;
    std::size_t input_count_{};
    std::vector<double> jacobian_;
};

} // namespace mpmc::flow

#endif // MPMC_FLOW_FUGACITY_EQUILIBRIUM_LINEARIZATION_HPP
