#include <mpmc/flash/detail/sw92_profile_c_sensitivity_local.hpp>

#include "../sw92_phase_assigned_pt/physical_sample6.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>

int main() {
    namespace fl = mpmc::flash;
    namespace sample6 = sw92_profile_c_sample6;

    const auto model = sample6::model();
    const auto source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(), model, 0.0);
    if (!source.accepted_phase_set_published()) {
        std::cerr << "Sample-6 publication unavailable\n";
        return 1;
    }
    const auto local = fl::detail::sw92_profile_c_sensitivity_local_system(
        source, model);
    const std::size_t chemical_rows =
        (local.phase_count - 1U) * local.component_count;
    double chemical_max = 0.0;
    double material_max = 0.0;
    std::size_t max_row = 0U;
    double global_max = 0.0;
    for (std::size_t row = 0; row < local.residual_count; ++row) {
        const double value = std::abs(local.values.values[row]);
        if (value > global_max) {
            global_max = value;
            max_row = row;
        }
        if (row < chemical_rows) chemical_max = std::max(chemical_max, value);
        else material_max = std::max(material_max, value);
    }
    std::cout << std::setprecision(17)
              << "chemical_max=" << chemical_max
              << " material_max=" << material_max
              << " global_max=" << global_max
              << " max_row=" << max_row << '\n';
    for (std::size_t row = 0; row < local.residual_count; ++row) {
        std::cout << "row[" << row << "]="
                  << local.values.values[row] << '\n';
    }
    return 0;
}
