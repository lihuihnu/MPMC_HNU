#ifndef MPMC_FLOW_SINGLE_PHASE_NATURAL_VARIABLE_HPP
#define MPMC_FLOW_SINGLE_PHASE_NATURAL_VARIABLE_HPP

#include <mpmc/flow/detail/validation.hpp>
#include <mpmc/flow/component_accumulation_time.hpp>
#include <mpmc/flow/energy_accumulation.hpp>
#include <mpmc/flow/phase_potential_upwind.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Compatibility entry point. New code can include only the responsibility it uses.
#include <mpmc/flow/single_phase_cell_state.hpp>
#include <mpmc/flow/single_phase_properties.hpp>
#include <mpmc/flow/single_phase_accumulation.hpp>
#include <mpmc/flow/single_phase_transport.hpp>

#endif // MPMC_FLOW_SINGLE_PHASE_NATURAL_VARIABLE_HPP
