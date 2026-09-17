#include "certificate_roundoff_auto.hpp"

#define certify(...) certify_with_full_roundoff_auto(request.pressure_pa, __VA_ARGS__)
#define main mpmc_cpa_association_certificate_base_main
#include "certificate_audit.cpp"
#undef main
#undef certify

int main() {
    return mpmc_cpa_association_certificate_base_main();
}
