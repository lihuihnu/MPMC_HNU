#include <mpmc/flow/detail/composition_coordinates.hpp>

#include <mpmc/flow/single_phase_properties.hpp>
#include <mpmc/flow/two_phase_properties.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace mpmc::flow;

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <class Operation>
void expect_out_of_range(Operation operation, const char* message) {
    try {
        operation();
    } catch (const std::out_of_range& error) {
        require(std::string{error.what()} == message,
                "coordinate failure must preserve its diagnostic");
        return;
    }
    throw std::runtime_error("coordinate failure must remain out_of_range");
}

// Differentiate the actual state reconstruction, using exactly representable
// binary inputs. This oracle neither inverts a column index nor calls the
// shared derivative helper. Every pivot and every unknown column is covered.
template <class IndependentColumn, class Derivative>
void check_reconstruction(std::size_t components, std::size_t phases,
                          const std::vector<std::size_t>& pivots,
                          IndependentColumn independent_column,
                          Derivative derivative) {
    constexpr double step = 1.0 / 1024.0;
    for (std::size_t column = 0U; column < phases * components + 1U; ++column) {
        for (std::size_t phase = 0U; phase < phases; ++phase) {
            std::vector<double> plus(components - 1U, 1.0 / 32.0);
            auto minus = plus;
            std::size_t rank = 0U;
            for (std::size_t component = 0U; component < components; ++component) {
                if (component == pivots[phase]) {
                    continue;
                }
                if (independent_column(phase, component) == column) {
                    plus[rank] += step;
                    minus[rank] -= step;
                }
                ++rank;
            }
            const auto positive = natural_variable_detail::reconstruct_positive_composition(
                plus, components - 1U, phase, pivots[phase]);
            const auto negative = natural_variable_detail::reconstruct_positive_composition(
                minus, components - 1U, phase, pivots[phase]);
            double sum = 0.0;
            for (std::size_t component = 0U; component < components; ++component) {
                const double expected = (positive[component] - negative[component]) / (2.0 * step);
                const double actual = derivative(phase, component, column);
                require(actual == expected, "composition derivative must match state reconstruction");
                sum += actual;
            }
            require(sum == 0.0, "a chart column must preserve composition normalization");
        }
    }
}

} // namespace

void composition_coordinates_contract() {
    for (std::size_t components = 2U; components <= 12U; ++components) {
        for (std::size_t pivot = 0U; pivot < components; ++pivot) {
            const NaturalVariableLayout1P one{
                NaturalVariableCompositionPivot1P::from_dependent_component(components, pivot)};
            check_reconstruction(components, 1U, {pivot},
                [&](std::size_t, std::size_t component) {
                    return one.independent_composition_unknown_index(component);
                },
                [&](std::size_t, std::size_t component, std::size_t column) {
                    return single_phase_detail::d_composition(one, component, column);
                });

            const std::array<std::size_t, 3> pivots{
                pivot, (pivot + 1U) % components, (pivot + 2U) % components};
            const NaturalVariableLayout2P two{
                NaturalVariableCompositionPivot2P::from_dependent_components(
                    components, {pivots[0], pivots[1]})};
            check_reconstruction(components, 2U, {pivots[0], pivots[1]},
                [&](std::size_t phase, std::size_t component) {
                    return two.independent_composition_unknown_index(phase, component);
                },
                [&](std::size_t phase, std::size_t component, std::size_t column) {
                    return two_phase_detail::d_composition(two, phase, component, column);
                });

            const NaturalVariableLayout3P three{
                NaturalVariableCompositionPivot3P::from_dependent_components(components, pivots)};
            check_reconstruction(components, 3U, {pivots.begin(), pivots.end()},
                [&](std::size_t phase, std::size_t component) {
                    return three.independent_composition_unknown_index(
                        static_cast<PhaseSlot3>(phase), component);
                },
                [&](std::size_t phase, std::size_t component, std::size_t column) {
                    return composition_coordinate_detail::derivative(three, phase, component, column);
                });

            // The descriptor path is consumed by well source terms and mixed charts.
            for (std::size_t phases = 1U; phases <= 3U; ++phases) {
                std::vector<std::size_t> active_pivots;
                for (std::size_t phase = 0U; phase < phases; ++phase) {
                    active_pivots.push_back(pivots[phase]);
                }
                const NaturalVariableLayoutDescriptor descriptor{components, phases, active_pivots};
                check_reconstruction(components, phases, active_pivots,
                    [&](std::size_t phase, std::size_t component) {
                        return descriptor.independent_composition_unknown_index(
                            static_cast<PhaseSlot3>(phase), component);
                    },
                    [&](std::size_t phase, std::size_t component, std::size_t column) {
                        return composition_coordinate_detail::derivative(descriptor, phase, component, column);
                    });
            }
        }
    }

    const NaturalVariableLayout1P one{3U};
    const NaturalVariableLayout2P two{3U};
    const NaturalVariableLayout3P three{3U};
    for (const auto column : {one.unknown_count(), std::numeric_limits<std::size_t>::max()}) {
        require(single_phase_detail::d_composition(one, 0U, column) == 0.0,
                "1P out-of-block columns retain zero derivatives");
    }
    for (const auto column : {two.unknown_count(), std::numeric_limits<std::size_t>::max()}) {
        require(two_phase_detail::d_composition(two, 1U, 0U, column) == 0.0,
                "2P out-of-block columns retain zero derivatives");
    }
    expect_out_of_range([&] { (void)two_phase_detail::d_composition(two, 2U, 0U, 0U); },
        "mpmc::flow: two-phase composition derivative index out of range");
    expect_out_of_range([&] { (void)two_phase_detail::d_composition(two, 0U, 3U, 0U); },
        "mpmc::flow: two-phase composition derivative index out of range");
    expect_out_of_range([&] {
        (void)composition_coordinate_detail::derivative(three, 0U, 0U, three.unknown_count());
    }, "mpmc::flow::NaturalVariableLayout3P: unknown index out of range");
    const auto descriptor = two.descriptor();
    expect_out_of_range([&] {
        (void)composition_coordinate_detail::derivative(descriptor, 0U, 0U, descriptor.unknown_count());
    }, "mpmc::flow::NaturalVariableLayoutDescriptor: unknown index out of range");
}
