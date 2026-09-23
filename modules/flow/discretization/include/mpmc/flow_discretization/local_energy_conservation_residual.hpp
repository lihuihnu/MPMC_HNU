#ifndef MPMC_FLOW_DISCRETIZATION_LOCAL_ENERGY_CONSERVATION_RESIDUAL_HPP
#define MPMC_FLOW_DISCRETIZATION_LOCAL_ENERGY_CONSERVATION_RESIDUAL_HPP

#include <mpmc/flow/energy_accumulation.hpp>
#include <mpmc/flow_discretization/energy_face_flux.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace mpmc::flow_discretization {

inline constexpr std::string_view
    local_energy_conservation_residual_convention =
        "flow_discretization/local-energy-conservation/multi-internal-face/v1";

enum class EnergyIncidentFaceLocalSide3D : std::uint8_t {
    owner,
    neighbour
};

struct LocalCellNormalizedEnergyFaceBinding3D {
    EnergyIncidentFaceLocalSide3D local_side{
        EnergyIncidentFaceLocalSide3D::owner};
    const NormalizedEnergyFaceContributionLinearization3D*
        contribution{};
};

struct LocalEnergyConservationNeighbourBlock3D {
    mpmc::mesh::LocalIndex face{
        mpmc::mesh::LocalIndex::value_type{0}};
    EnergyIncidentFaceLocalSide3D local_side{
        EnergyIncidentFaceLocalSide3D::owner};
    mpmc::flow::NaturalVariableStateIdentity3P
        neighbour_state_identity;
    std::vector<double> gradient;

    [[nodiscard]] double d_residual(
        std::size_t column) const {
        return gradient.at(column);
    }
};

/// One cell's energy equation restricted to accumulation + internal-face
/// advection/conduction.
///
/// R_E,c = R_E,c^acc + sum_f R_E,c^(face,V)
///
/// Units: W / bulk-m^3.
struct LocalEnergyConservationResidualLinearization3D {
    static constexpr std::string_view convention =
        local_energy_conservation_residual_convention;

    mpmc::flow::NaturalVariableStateIdentity3P
        cell_state_identity;
    double cell_bulk_volume_m3{};
    double residual_w_per_bulk_m3{};
    std::vector<double> local_gradient;
    std::vector<LocalEnergyConservationNeighbourBlock3D>
        neighbour_blocks;

    [[nodiscard]] double d_local(
        std::size_t column) const {
        return local_gradient.at(column);
    }
};

namespace local_energy_conservation_detail {

[[nodiscard]] inline bool near_roundoff(
    double first,
    double second,
    double scale = 0.0) {
    if (!std::isfinite(first) ||
        !std::isfinite(second) ||
        !std::isfinite(scale) ||
        scale < 0.0) {
        return false;
    }
    const double reference =
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second),
             scale});
    return std::abs(first - second) <=
        16384.0 *
            std::numeric_limits<double>::epsilon() *
            reference;
}

inline void validate_face(
    const NormalizedEnergyFaceContributionLinearization3D&
        face) {
    const std::size_t owner_q =
        face.owner_state_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        face.neighbour_state_identity.layout.unknown_count();

    if (owner_q == 0U ||
        neighbour_q == 0U ||
        !std::isfinite(
            face.bulk_volume.owner_bulk_volume_m3) ||
        !(face.bulk_volume.owner_bulk_volume_m3 > 0.0) ||
        !std::isfinite(
            face.bulk_volume.neighbour_bulk_volume_m3) ||
        !(face.bulk_volume.neighbour_bulk_volume_m3 > 0.0) ||
        !std::isfinite(
            face.owner_contribution_w_per_bulk_m3) ||
        !std::isfinite(
            face.neighbour_contribution_w_per_bulk_m3) ||
        face.owner_row_owner_column_gradient.size() !=
            owner_q ||
        face.neighbour_row_owner_column_gradient.size() !=
            owner_q ||
        face.owner_row_neighbour_column_gradient.size() !=
            neighbour_q ||
        face.neighbour_row_neighbour_column_gradient.size() !=
            neighbour_q) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: malformed normalized energy face contribution");
    }

    const double weighted =
        face.bulk_volume.owner_bulk_volume_m3 *
            face.owner_contribution_w_per_bulk_m3 +
        face.bulk_volume.neighbour_bulk_volume_m3 *
            face.neighbour_contribution_w_per_bulk_m3;
    const double scale =
        face.bulk_volume.owner_bulk_volume_m3 *
            std::abs(
                face.owner_contribution_w_per_bulk_m3) +
        face.bulk_volume.neighbour_bulk_volume_m3 *
            std::abs(
                face.neighbour_contribution_w_per_bulk_m3);
    if (!near_roundoff(
            weighted,
            0.0,
            scale)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: normalized energy face violates volume-weighted conservation");
    }

    const auto check_gradient =
        [&](const std::vector<double>& owner_row,
            const std::vector<double>& neighbour_row,
            std::size_t q) {
            for (std::size_t column = 0U;
                 column < q;
                 ++column) {
                if (!std::isfinite(
                        owner_row[column]) ||
                    !std::isfinite(
                        neighbour_row[column])) {
                    throw std::invalid_argument(
                        "mpmc::flow_discretization: normalized energy face contains non-finite derivative");
                }
                const double sum =
                    face.bulk_volume.owner_bulk_volume_m3 *
                        owner_row[column] +
                    face.bulk_volume.neighbour_bulk_volume_m3 *
                        neighbour_row[column];
                const double derivative_scale =
                    face.bulk_volume.owner_bulk_volume_m3 *
                        std::abs(owner_row[column]) +
                    face.bulk_volume.neighbour_bulk_volume_m3 *
                        std::abs(neighbour_row[column]);
                if (!near_roundoff(
                        sum,
                        0.0,
                        derivative_scale)) {
                    throw std::invalid_argument(
                        "mpmc::flow_discretization: normalized energy face Jacobian violates volume-weighted conservation");
                }
            }
        };

    check_gradient(
        face.owner_row_owner_column_gradient,
        face.neighbour_row_owner_column_gradient,
        owner_q);
    check_gradient(
        face.owner_row_neighbour_column_gradient,
        face.neighbour_row_neighbour_column_gradient,
        neighbour_q);
}

} // namespace local_energy_conservation_detail

[[nodiscard]] inline
LocalEnergyConservationResidualLinearization3D
build_local_energy_conservation_residual(
    const mpmc::flow::
        BackwardEulerEnergyAccumulationResidual3P&
            accumulation,
    double cell_bulk_volume_m3,
    std::span<
        const LocalCellNormalizedEnergyFaceBinding3D>
        incident_faces) {
    using mpmc::flow::energy_accumulation_detail::
        same_state_identity;
    using namespace local_energy_conservation_detail;

    const auto& cell_identity =
        accumulation.current_state_identity;
    const std::size_t local_q =
        cell_identity.layout.unknown_count();

    if (local_q == 0U ||
        accumulation.input_count !=
            local_q ||
        accumulation.gradient.size() !=
            local_q ||
        !std::isfinite(
            accumulation.residual_w_per_bulk_m3) ||
        !std::isfinite(
            cell_bulk_volume_m3) ||
        !(cell_bulk_volume_m3 > 0.0)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization: malformed energy accumulation/local bulk-volume payload");
    }
    for (double derivative :
         accumulation.gradient) {
        if (!std::isfinite(derivative)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: energy accumulation gradient contains non-finite derivative");
        }
    }

    LocalEnergyConservationResidualLinearization3D
        result{
            cell_identity,
            cell_bulk_volume_m3,
            accumulation.residual_w_per_bulk_m3,
            accumulation.gradient,
            {}};
    result.neighbour_blocks.reserve(
        incident_faces.size());

    for (const auto& binding :
         incident_faces) {
        if (binding.contribution ==
            nullptr) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: incident normalized energy face binding must not be null");
        }
        const auto& face =
            *binding.contribution;
        validate_face(face);

        for (const auto& existing :
             result.neighbour_blocks) {
            if (existing.face ==
                face.face) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization: duplicate incident energy face");
            }
        }

        bool local_is_owner = false;
        switch (binding.local_side) {
        case EnergyIncidentFaceLocalSide3D::owner:
            local_is_owner = true;
            break;
        case EnergyIncidentFaceLocalSide3D::neighbour:
            local_is_owner = false;
            break;
        default:
            throw std::invalid_argument(
                "mpmc::flow_discretization: unsupported incident energy-face local side");
        }

        const auto& selected_identity =
            local_is_owner
                ? face.owner_state_identity
                : face.neighbour_state_identity;
        const auto& neighbour_identity =
            local_is_owner
                ? face.neighbour_state_identity
                : face.owner_state_identity;
        const double face_local_volume =
            local_is_owner
                ? face.bulk_volume.owner_bulk_volume_m3
                : face.bulk_volume.neighbour_bulk_volume_m3;

        if (!same_state_identity(
                selected_identity,
                cell_identity) ||
            !near_roundoff(
                face_local_volume,
                cell_bulk_volume_m3,
                std::max(
                    face_local_volume,
                    cell_bulk_volume_m3))) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: incident energy face does not match local state/chart/bulk volume");
        }

        const double local_value =
            local_is_owner
                ? face.owner_contribution_w_per_bulk_m3
                : face.neighbour_contribution_w_per_bulk_m3;
        const auto& local_gradient =
            local_is_owner
                ? face.owner_row_owner_column_gradient
                : face.neighbour_row_neighbour_column_gradient;
        const auto& neighbour_gradient =
            local_is_owner
                ? face.owner_row_neighbour_column_gradient
                : face.neighbour_row_owner_column_gradient;

        if (local_gradient.size() !=
                local_q ||
            neighbour_gradient.size() !=
                neighbour_identity.layout.unknown_count()) {
            throw std::invalid_argument(
                "mpmc::flow_discretization: selected energy face Jacobian shape mismatch");
        }

        result.residual_w_per_bulk_m3 +=
            local_value;
        if (!std::isfinite(
                result.residual_w_per_bulk_m3)) {
            throw std::range_error(
                "mpmc::flow_discretization: local energy residual became non-finite");
        }

        for (std::size_t column = 0U;
             column < local_q;
             ++column) {
            result.local_gradient[column] +=
                local_gradient[column];
            if (!std::isfinite(
                    result.local_gradient[column])) {
                throw std::range_error(
                    "mpmc::flow_discretization: local energy diagonal Jacobian became non-finite");
            }
        }

        result.neighbour_blocks.push_back(
            LocalEnergyConservationNeighbourBlock3D{
                face.face,
                binding.local_side,
                neighbour_identity,
                neighbour_gradient});
    }

    return result;
}

} // namespace mpmc::flow_discretization

#endif // MPMC_FLOW_DISCRETIZATION_LOCAL_ENERGY_CONSERVATION_RESIDUAL_HPP
