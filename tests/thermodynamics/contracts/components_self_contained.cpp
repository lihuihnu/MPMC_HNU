#include <mpmc/thermodynamics/components.hpp>
#include <mpmc/thermodynamics/components.hpp>

#include <type_traits>

bool components_header_compiles() {
    using Components = mpmc::thermodynamics::OrderedComponents;
    static_assert(std::is_move_constructible_v<Components>);
    static_assert(!std::is_copy_assignable_v<Components>);
    return true;
}
