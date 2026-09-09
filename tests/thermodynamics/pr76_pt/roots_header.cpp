#include <mpmc/thermodynamics/pr76_roots.hpp>
#include <mpmc/thermodynamics/pr76_roots.hpp>

bool roots_plain_header() {
    const auto result=mpmc::thermodynamics::pr76_roots(51.0/128,1.0/16);
    return result.status==mpmc::thermodynamics::Pr76RootStatus::success && result.count==3;
}
