#include <mpmc/discretization/transmissibility_admissibility_3d.hpp>
#include <mpmc/discretization/tpfa_half_connection_3d.hpp>
#include <mpmc/discretization/tpfa_half_transmissibility_3d.hpp>
#include <mpmc/discretization/tpfa_static_face_transmissibility_3d.hpp>
#include <mpmc/discretization/tpfa_internal_face_transmissibility_snapshot_3d.hpp>
#include <mpmc/mesh/active_corner_point.hpp>
#include <mpmc/mesh/cell_face_geometric_operator_3d.hpp>
#include <mpmc/mesh/grdecl.hpp>
#include <mpmc/mesh/permeability_tensor_3d.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {
namespace mesh = mpmc::mesh;
namespace discretization = mpmc::discretization;

void require(
    bool condition,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(
            std::string(where.file_name()) + ":" +
            std::to_string(where.line()) + ": " +
            std::string(message));
    }
}

void require_close(
    double actual,
    double expected,
    double tolerance,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    if (!std::isfinite(actual) ||
        std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(
            std::string(where.file_name()) + ":" +
            std::to_string(where.line()) + ": " +
            std::string(message));
    }
}

template <typename Exception, typename Function>
void expect_throw(Function&& function) {
    bool caught = false;
    try {
        std::forward<Function>(function)();
    } catch (const Exception&) {
        caught = true;
    }
    require(caught, "expected exception was not thrown");
}

std::string grdecl_two_cell_all_active_fixture() {
    return R"grdecl(-- minimal 2x1x1 Cartesian corner-point deck
SPECGRID
  2 1 1 1 F /
COORD
  0 0 0   0 0 1
  1 0 0   1 0 1
  2 0 0   2 0 1
  0 1 0   0 1 1
  1 1 0   1 1 1
  2 1 0   2 1 1 /
ZCORN
  8*0 8*1 /
ACTNUM
  1 1 /
PORO
  0.20 0.35 /
PERMX
  100 200 /
PERMY
  50 75 /
PERMZ
  1.0D+1 2.0D+1 /
)grdecl";
}

std::string grdecl_vertical_two_cell_fixture() {
    return R"grdecl(SPECGRID
  1 1 2 1 F /
COORD
  0 0 0   0 0 2
  1 0 0   1 0 2
  0 1 0   0 1 2
  1 1 0   1 1 2 /
ZCORN
  4*0 8*1 4*2 /
ACTNUM
  2*1 /
PORO
  0.20 0.30 /
PERMX
  100 110 /
PERMY
  90 95 /
PERMZ
  10 12 /
)grdecl";
}

std::string grdecl_skewed_two_cell_fixture() {
    return R"grdecl(-- valid skewed 2x1x1 corner-point deck
SPECGRID
  2 1 1 1 F /
COORD
  0.0 0 0   0.0 0 1
  1.2 0 0   1.2 0 1
  2.0 0 0   2.0 0 1
  0.0 1 0   0.0 1 1
  0.8 1 0   0.8 1 1
  2.0 1 0   2.0 1 1 /
ZCORN
  8*0 8*1 /
ACTNUM
  2*1 /
PORO
  0.20 0.35 /
PERMX
  100 200 /
PERMY
  50 75 /
PERMZ
  10 20 /
)grdecl";
}

const mesh::DenseFieldSnapshot&
active_corner_point_field(
    const mesh::ActiveCornerPointGrid& processed,
    std::string_view id) {
    const auto found =
        std::find_if(
            processed.cell_fields.begin(),
            processed.cell_fields.end(),
            [id](const auto& field) {
                return field.metadata().id == id;
            });
    require(
        found != processed.cell_fields.end(),
        "processed active field ID missing");
    return *found;
}

mesh::LocalIndex only_shared_face(
    const mesh::Topology& topology) {
    const auto& face_cells =
        topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    std::optional<mesh::LocalIndex> found;
    for (std::size_t face = 0U;
         face < topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        if (face_cells.adjacent(local).size() != 2U) {
            continue;
        }
        require(
            !found.has_value(),
            "expected exactly one shared face");
        found = local;
    }
    require(found.has_value(), "expected one shared face");
    return *found;
}

mesh::CellCartesianDiagonalPermeability3D
processed_diagonal_permeability(
    const mesh::ActiveCornerPointGrid& processed) {
    return mesh::make_cell_cartesian_diagonal_permeability_3d(
        processed.topology,
        active_corner_point_field(
            processed, "PERMX"),
        active_corner_point_field(
            processed, "PERMY"),
        active_corner_point_field(
            processed, "PERMZ"));
}

void tpfa_half_connection_3d() {
    const auto horizontal =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_two_cell_all_active_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto horizontal_geometry =
        mesh::make_cell_face_geometric_operator_3d(
            horizontal.topology,
            horizontal.vertex_coordinates_m,
            horizontal.face_geometry);
    const auto horizontal_permeability =
        processed_diagonal_permeability(
            horizontal);
    const auto horizontal_interface =
        only_shared_face(
            horizontal.topology);

    const auto owner =
        discretization::make_owner_tpfa_half_connection_coefficient_3d(
            horizontal_geometry,
            horizontal_permeability,
            horizontal_interface);
    require(
        owner.projection ==
            discretization::TpfaHalfConnectionProjection3D::
                positive_projection,
        "horizontal owner half connection positive projection");
    require_close(
        owner.normal_permeability_displacement_m3,
        50.0e-15,
        1.0e-28,
        "horizontal owner nTKd");
    require_close(
        owner.squared_distance_m2,
        0.25,
        1.0e-14,
        "horizontal owner dTd");
    require_close(
        owner.coefficient_m,
        200.0e-15,
        1.0e-28,
        "horizontal owner one-sided coefficient");

    const auto neighbour =
        discretization::make_neighbour_tpfa_half_connection_coefficient_3d(
            horizontal_geometry,
            horizontal_permeability,
            horizontal_interface);
    require(
        neighbour.projection ==
            discretization::TpfaHalfConnectionProjection3D::
                positive_projection,
        "horizontal neighbour half connection positive projection");
    require_close(
        neighbour.normal_permeability_displacement_m3,
        100.0e-15,
        1.0e-28,
        "horizontal neighbour nTKd");
    require_close(
        neighbour.squared_distance_m2,
        0.25,
        1.0e-14,
        "horizontal neighbour dTd");
    require_close(
        neighbour.coefficient_m,
        400.0e-15,
        1.0e-28,
        "horizontal neighbour one-sided coefficient");

    const auto vertical =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_vertical_two_cell_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto vertical_geometry =
        mesh::make_cell_face_geometric_operator_3d(
            vertical.topology,
            vertical.vertex_coordinates_m,
            vertical.face_geometry);
    const auto vertical_permeability =
        processed_diagonal_permeability(
            vertical);
    const auto vertical_interface =
        only_shared_face(
            vertical.topology);

    const auto vertical_owner =
        discretization::make_owner_tpfa_half_connection_coefficient_3d(
            vertical_geometry,
            vertical_permeability,
            vertical_interface);
    const auto vertical_neighbour =
        discretization::make_neighbour_tpfa_half_connection_coefficient_3d(
            vertical_geometry,
            vertical_permeability,
            vertical_interface);
    require_close(
        vertical_owner.coefficient_m,
        20.0e-15,
        1.0e-28,
        "vertical owner one-sided coefficient");
    require_close(
        vertical_neighbour.coefficient_m,
        24.0e-15,
        1.0e-28,
        "vertical neighbour one-sided coefficient");

    const auto skewed =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_skewed_two_cell_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto skewed_geometry =
        mesh::make_cell_face_geometric_operator_3d(
            skewed.topology,
            skewed.vertex_coordinates_m,
            skewed.face_geometry);
    const auto skewed_permeability =
        processed_diagonal_permeability(
            skewed);
    const auto skewed_interface =
        only_shared_face(
            skewed.topology);
    const double sqrt29 =
        std::sqrt(29.0);

    const auto skewed_owner =
        discretization::make_owner_tpfa_half_connection_coefficient_3d(
            skewed_geometry,
            skewed_permeability,
            skewed_interface);
    const auto skewed_neighbour =
        discretization::make_neighbour_tpfa_half_connection_coefficient_3d(
            skewed_geometry,
            skewed_permeability,
            skewed_interface);
    require_close(
        skewed_owner.normal_permeability_displacement_m3,
        250.0e-15 / sqrt29,
        1.0e-28,
        "skewed owner nTKd");
    require_close(
        skewed_owner.coefficient_m,
        1000.0e-15 / sqrt29,
        1.0e-28,
        "skewed owner one-sided coefficient");
    require_close(
        skewed_neighbour.normal_permeability_displacement_m3,
        500.0e-15 / sqrt29,
        1.0e-28,
        "skewed neighbour nTKd");
    require_close(
        skewed_neighbour.coefficient_m,
        2000.0e-15 / sqrt29,
        1.0e-28,
        "skewed neighbour one-sided coefficient");

    const auto strict_admissibility =
        discretization::classify_internal_face_transmissibility_admissibility(
            skewed_geometry,
            skewed_permeability,
            skewed_interface,
            discretization::TransmissibilityGeometryAdmissibilityPolicy3D{
                0.0},
            discretization::KOrthogonalityAdmissibilityPolicy3D{
                0.0});
    require(
        strict_admissibility.disposition ==
            discretization::CombinedTransmissibilityAdmissibilityDisposition3D::
                requires_geometry_and_k_non_orthogonal_treatment,
        "skewed one-sided coefficients do not bypass strict admissibility");
}

void tpfa_half_connection_3d_invalid() {
    const auto positive_tensor =
        mesh::CartesianDiagonalPermeabilityTensor3D{
            100.0e-15,
            50.0e-15,
            10.0e-15};
    const auto unit_x =
        mesh::UnitVector3D{
            1.0, 0.0, 0.0};
    const auto half_x =
        mesh::Displacement3D{
            0.5, 0.0, 0.0};

    const auto zero_projection =
        discretization::make_tpfa_half_connection_coefficient_3d(
            mesh::CartesianDiagonalPermeabilityTensor3D{
                0.0,
                50.0e-15,
                10.0e-15},
            unit_x,
            half_x);
    require(
        zero_projection.projection ==
            discretization::TpfaHalfConnectionProjection3D::
                zero_projection &&
            zero_projection
                    .normal_permeability_displacement_m3 ==
                0.0 &&
            zero_projection.coefficient_m == 0.0,
        "zero normal permeability projection remains explicit zero coefficient");
    require_close(
        zero_projection.squared_distance_m2,
        0.25,
        0.0,
        "zero projection keeps positive dTd");

    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_tpfa_half_connection_coefficient_3d(
                mesh::CartesianDiagonalPermeabilityTensor3D{
                    -1.0,
                    1.0,
                    1.0},
                unit_x,
                half_x);
        });
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_tpfa_half_connection_coefficient_3d(
                positive_tensor,
                mesh::UnitVector3D{
                    2.0, 0.0, 0.0},
                half_x);
        });
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_tpfa_half_connection_coefficient_3d(
                positive_tensor,
                unit_x,
                mesh::Displacement3D{
                    0.0, 0.0, 0.0});
        });
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_tpfa_half_connection_coefficient_3d(
                positive_tensor,
                unit_x,
                mesh::Displacement3D{
                    std::numeric_limits<double>::
                        quiet_NaN(),
                    0.0,
                    0.0});
        });
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_tpfa_half_connection_coefficient_3d(
                positive_tensor,
                mesh::UnitVector3D{
                    -1.0, 0.0, 0.0},
                half_x);
        });

    const auto processed =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_two_cell_all_active_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto geometry =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);
    const auto permeability =
        processed_diagonal_permeability(
            processed);

    const auto& face_cells =
        processed.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    bool checked_boundary = false;
    for (std::size_t face = 0U;
         face < processed.topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        if (face_cells.adjacent(local).size() != 1U) {
            continue;
        }
        const auto boundary_owner =
            discretization::make_owner_tpfa_half_connection_coefficient_3d(
                geometry,
                permeability,
                local);
        require(
            boundary_owner.squared_distance_m2 > 0.0,
            "boundary owner half connection remains independently defined");
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::make_neighbour_tpfa_half_connection_coefficient_3d(
                    geometry,
                    permeability,
                    local);
            });
        checked_boundary = true;
        break;
    }
    require(
        checked_boundary,
        "TPFA half-connection invalid test found boundary face");

    const auto mismatched_permeability =
        mesh::CellCartesianDiagonalPermeability3D{
            {mesh::CartesianDiagonalPermeabilityTensor3D{
                100.0e-15,
                50.0e-15,
                10.0e-15}}};
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_owner_tpfa_half_connection_coefficient_3d(
                geometry,
                mismatched_permeability,
                only_shared_face(
                    processed.topology));
        });
}


void tpfa_half_transmissibility_3d() {
    const auto horizontal =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_two_cell_all_active_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto horizontal_geometry =
        mesh::make_cell_face_geometric_operator_3d(
            horizontal.topology,
            horizontal.vertex_coordinates_m,
            horizontal.face_geometry);
    const auto horizontal_permeability =
        processed_diagonal_permeability(
            horizontal);
    const auto horizontal_interface =
        only_shared_face(
            horizontal.topology);

    const auto horizontal_owner =
        discretization::make_owner_area_scaled_tpfa_half_transmissibility_3d(
            horizontal_geometry,
            horizontal_permeability,
            horizontal_interface);
    const auto horizontal_neighbour =
        discretization::make_neighbour_area_scaled_tpfa_half_transmissibility_3d(
            horizontal_geometry,
            horizontal_permeability,
            horizontal_interface);

    require_close(
        horizontal_owner.face_area_m2,
        1.0,
        1.0e-14,
        "horizontal owner TPFA half area");
    require_close(
        horizontal_owner
            .half_connection.coefficient_m,
        200.0e-15,
        1.0e-28,
        "horizontal owner half coefficient retained");
    require_close(
        horizontal_owner
            .half_transmissibility_m3,
        200.0e-15,
        1.0e-28,
        "horizontal owner area-scaled half transmissibility");
    require_close(
        horizontal_neighbour
            .half_transmissibility_m3,
        400.0e-15,
        1.0e-28,
        "horizontal neighbour area-scaled half transmissibility");

    const auto vertical =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_vertical_two_cell_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto vertical_geometry =
        mesh::make_cell_face_geometric_operator_3d(
            vertical.topology,
            vertical.vertex_coordinates_m,
            vertical.face_geometry);
    const auto vertical_permeability =
        processed_diagonal_permeability(
            vertical);
    const auto vertical_interface =
        only_shared_face(
            vertical.topology);
    const auto vertical_owner =
        discretization::make_owner_area_scaled_tpfa_half_transmissibility_3d(
            vertical_geometry,
            vertical_permeability,
            vertical_interface);
    const auto vertical_neighbour =
        discretization::make_neighbour_area_scaled_tpfa_half_transmissibility_3d(
            vertical_geometry,
            vertical_permeability,
            vertical_interface);
    require_close(
        vertical_owner.half_transmissibility_m3,
        20.0e-15,
        1.0e-28,
        "vertical owner area-scaled half transmissibility");
    require_close(
        vertical_neighbour.half_transmissibility_m3,
        24.0e-15,
        1.0e-28,
        "vertical neighbour area-scaled half transmissibility");

    const auto skewed =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_skewed_two_cell_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto skewed_geometry =
        mesh::make_cell_face_geometric_operator_3d(
            skewed.topology,
            skewed.vertex_coordinates_m,
            skewed.face_geometry);
    const auto skewed_permeability =
        processed_diagonal_permeability(
            skewed);
    const auto skewed_interface =
        only_shared_face(
            skewed.topology);

    const auto skewed_owner =
        discretization::make_owner_area_scaled_tpfa_half_transmissibility_3d(
            skewed_geometry,
            skewed_permeability,
            skewed_interface);
    const auto skewed_neighbour =
        discretization::make_neighbour_area_scaled_tpfa_half_transmissibility_3d(
            skewed_geometry,
            skewed_permeability,
            skewed_interface);
    require_close(
        skewed_owner.face_area_m2,
        std::sqrt(29.0) / 5.0,
        1.0e-14,
        "skewed shared face area retained");
    require_close(
        skewed_owner
            .half_connection.coefficient_m,
        1000.0e-15 / std::sqrt(29.0),
        1.0e-28,
        "skewed owner coefficient retained before area scaling");
    require_close(
        skewed_owner.half_transmissibility_m3,
        200.0e-15,
        1.0e-28,
        "skewed owner area-scaled half transmissibility");
    require_close(
        skewed_neighbour.half_transmissibility_m3,
        400.0e-15,
        1.0e-28,
        "skewed neighbour area-scaled half transmissibility");

    const auto strict_admissibility =
        discretization::classify_internal_face_transmissibility_admissibility(
            skewed_geometry,
            skewed_permeability,
            skewed_interface,
            discretization::TransmissibilityGeometryAdmissibilityPolicy3D{
                0.0},
            discretization::KOrthogonalityAdmissibilityPolicy3D{
                0.0});
    require(
        strict_admissibility.disposition ==
            discretization::CombinedTransmissibilityAdmissibilityDisposition3D::
                requires_geometry_and_k_non_orthogonal_treatment,
        "area scaling does not bypass strict geometry/K admissibility");
}

void tpfa_half_transmissibility_3d_invalid() {
    const auto zero_half =
        discretization::TpfaHalfConnectionCoefficient3D{
            discretization::TpfaHalfConnectionProjection3D::
                zero_projection,
            0.0,
            0.25,
            0.0};
    const auto zero_scaled =
        discretization::make_area_scaled_tpfa_half_transmissibility_3d(
            2.0,
            zero_half);
    require(
        zero_scaled.half_connection.projection ==
            discretization::TpfaHalfConnectionProjection3D::
                zero_projection &&
            zero_scaled.face_area_m2 == 2.0 &&
            zero_scaled.half_transmissibility_m3 ==
                0.0,
        "zero half projection survives positive area scaling");

    const auto positive_half =
        discretization::TpfaHalfConnectionCoefficient3D{
            discretization::TpfaHalfConnectionProjection3D::
                positive_projection,
            50.0e-15,
            0.25,
            200.0e-15};

    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_area_scaled_tpfa_half_transmissibility_3d(
                0.0,
                positive_half);
        });
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_area_scaled_tpfa_half_transmissibility_3d(
                -1.0,
                positive_half);
        });
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_area_scaled_tpfa_half_transmissibility_3d(
                std::numeric_limits<double>::
                    quiet_NaN(),
                positive_half);
        });

    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_area_scaled_tpfa_half_transmissibility_3d(
                1.0,
                discretization::TpfaHalfConnectionCoefficient3D{
                    discretization::TpfaHalfConnectionProjection3D::
                        positive_projection,
                    50.0e-15,
                    0.25,
                    0.0});
        });
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_area_scaled_tpfa_half_transmissibility_3d(
                1.0,
                discretization::TpfaHalfConnectionCoefficient3D{
                    discretization::TpfaHalfConnectionProjection3D::
                        zero_projection,
                    50.0e-15,
                    0.25,
                    200.0e-15});
        });
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_area_scaled_tpfa_half_transmissibility_3d(
                1.0,
                discretization::TpfaHalfConnectionCoefficient3D{
                    discretization::TpfaHalfConnectionProjection3D::
                        positive_projection,
                    50.0e-15,
                    0.5,
                    200.0e-15});
        });
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_area_scaled_tpfa_half_transmissibility_3d(
                1.0,
                discretization::TpfaHalfConnectionCoefficient3D{
                    static_cast<
                        discretization::TpfaHalfConnectionProjection3D>(
                            255U),
                    50.0e-15,
                    0.25,
                    200.0e-15});
        });

    const auto processed =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_two_cell_all_active_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto geometry =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);
    const auto permeability =
        processed_diagonal_permeability(
            processed);

    const auto& face_cells =
        processed.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    bool checked_boundary = false;
    for (std::size_t face = 0U;
         face < processed.topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        if (face_cells.adjacent(local).size() != 1U) {
            continue;
        }
        const auto owner_half =
            discretization::make_owner_area_scaled_tpfa_half_transmissibility_3d(
                geometry,
                permeability,
                local);
        require(
            owner_half.face_area_m2 > 0.0 &&
                owner_half.half_transmissibility_m3 >=
                    0.0,
            "boundary owner area-scaled half transmissibility is independently defined");
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::make_neighbour_area_scaled_tpfa_half_transmissibility_3d(
                    geometry,
                    permeability,
                    local);
            });
        checked_boundary = true;
        break;
    }
    require(
        checked_boundary,
        "area-scaled TPFA half transmissibility invalid test found boundary face");
}


discretization::TpfaAreaScaledHalfTransmissibility3D
synthetic_area_scaled_tpfa_half(
    double face_area_m2,
    double half_transmissibility_m3,
    discretization::TpfaHalfConnectionProjection3D projection) {
    if (projection ==
        discretization::TpfaHalfConnectionProjection3D::
            zero_projection) {
        return discretization::TpfaAreaScaledHalfTransmissibility3D{
            discretization::TpfaHalfConnectionCoefficient3D{
                projection,
                0.0,
                0.25,
                0.0},
            face_area_m2,
            half_transmissibility_m3};
    }

    const double coefficient_m =
        half_transmissibility_m3 /
        face_area_m2;
    return discretization::TpfaAreaScaledHalfTransmissibility3D{
        discretization::TpfaHalfConnectionCoefficient3D{
            projection,
            0.25 * coefficient_m,
            0.25,
            coefficient_m},
        face_area_m2,
        half_transmissibility_m3};
}

void tpfa_static_face_transmissibility_3d() {
    const auto horizontal =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_two_cell_all_active_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto horizontal_geometry =
        mesh::make_cell_face_geometric_operator_3d(
            horizontal.topology,
            horizontal.vertex_coordinates_m,
            horizontal.face_geometry);
    const auto horizontal_permeability =
        processed_diagonal_permeability(
            horizontal);
    const auto horizontal_interface =
        only_shared_face(
            horizontal.topology);

    const auto horizontal_owner =
        discretization::make_owner_area_scaled_tpfa_half_transmissibility_3d(
            horizontal_geometry,
            horizontal_permeability,
            horizontal_interface);
    const auto horizontal_neighbour =
        discretization::make_neighbour_area_scaled_tpfa_half_transmissibility_3d(
            horizontal_geometry,
            horizontal_permeability,
            horizontal_interface);
    const auto horizontal_face =
        discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
            horizontal_owner,
            horizontal_neighbour);
    const auto horizontal_swapped =
        discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
            horizontal_neighbour,
            horizontal_owner);
    const auto horizontal_wrapper =
        discretization::make_internal_face_static_tpfa_transmissibility_3d(
            horizontal_geometry,
            horizontal_permeability,
            horizontal_interface);

    require(
        horizontal_face.disposition ==
            discretization::TpfaStaticFaceTransmissibilityDisposition3D::
                positive_harmonic_combination,
        "horizontal internal face uses positive harmonic combination");
    require_close(
        horizontal_face.face_area_m2,
        1.0,
        1.0e-14,
        "horizontal combined face area");
    require_close(
        horizontal_face.face_transmissibility_m3,
        400.0e-15 / 3.0,
        1.0e-28,
        "horizontal harmonic face transmissibility");
    require(
        horizontal_swapped.disposition ==
            horizontal_face.disposition,
        "owner/neighbour swap preserves harmonic disposition");
    require_close(
        horizontal_swapped.face_area_m2,
        horizontal_face.face_area_m2,
        0.0,
        "owner/neighbour swap preserves face area");
    require_close(
        horizontal_swapped.face_transmissibility_m3,
        horizontal_face.face_transmissibility_m3,
        0.0,
        "owner/neighbour swap preserves face transmissibility");
    require_close(
        horizontal_wrapper.face_transmissibility_m3,
        horizontal_face.face_transmissibility_m3,
        0.0,
        "internal-face wrapper matches explicit half combination");

    const auto vertical =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_vertical_two_cell_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto vertical_geometry =
        mesh::make_cell_face_geometric_operator_3d(
            vertical.topology,
            vertical.vertex_coordinates_m,
            vertical.face_geometry);
    const auto vertical_permeability =
        processed_diagonal_permeability(
            vertical);
    const auto vertical_interface =
        only_shared_face(
            vertical.topology);
    const auto vertical_face =
        discretization::make_internal_face_static_tpfa_transmissibility_3d(
            vertical_geometry,
            vertical_permeability,
            vertical_interface);
    require_close(
        vertical_face.face_transmissibility_m3,
        120.0e-15 / 11.0,
        1.0e-28,
        "vertical harmonic face transmissibility");

    const auto skewed =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_skewed_two_cell_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto skewed_geometry =
        mesh::make_cell_face_geometric_operator_3d(
            skewed.topology,
            skewed.vertex_coordinates_m,
            skewed.face_geometry);
    const auto skewed_permeability =
        processed_diagonal_permeability(
            skewed);
    const auto skewed_interface =
        only_shared_face(
            skewed.topology);
    const auto skewed_face =
        discretization::make_internal_face_static_tpfa_transmissibility_3d(
            skewed_geometry,
            skewed_permeability,
            skewed_interface);
    require_close(
        skewed_face.face_area_m2,
        std::sqrt(29.0) / 5.0,
        1.0e-14,
        "skewed combined face area");
    require_close(
        skewed_face.face_transmissibility_m3,
        400.0e-15 / 3.0,
        1.0e-28,
        "skewed harmonic static face transmissibility");

    const auto strict_admissibility =
        discretization::classify_internal_face_transmissibility_admissibility(
            skewed_geometry,
            skewed_permeability,
            skewed_interface,
            discretization::TransmissibilityGeometryAdmissibilityPolicy3D{
                0.0},
            discretization::KOrthogonalityAdmissibilityPolicy3D{
                0.0});
    require(
        strict_admissibility.disposition ==
            discretization::CombinedTransmissibilityAdmissibilityDisposition3D::
                requires_geometry_and_k_non_orthogonal_treatment,
        "static harmonic combination does not bypass strict geometry/K admissibility");

    const auto huge_owner =
        synthetic_area_scaled_tpfa_half(
            1.0,
            1.0e300,
            discretization::TpfaHalfConnectionProjection3D::
                positive_projection);
    const auto huge_neighbour =
        synthetic_area_scaled_tpfa_half(
            1.0,
            5.0e299,
            discretization::TpfaHalfConnectionProjection3D::
                positive_projection);
    const auto huge_combined =
        discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
            huge_owner,
            huge_neighbour);
    require(
        std::isfinite(
            huge_combined.face_transmissibility_m3),
        "stable harmonic form avoids positive-half product overflow");
    require_close(
        huge_combined.face_transmissibility_m3,
        1.0e300 / 3.0,
        1.0e286,
        "large finite harmonic combination");
}

void tpfa_static_face_transmissibility_3d_invalid() {
    const auto positive_owner =
        synthetic_area_scaled_tpfa_half(
            1.0,
            2.0,
            discretization::TpfaHalfConnectionProjection3D::
                positive_projection);
    const auto positive_neighbour =
        synthetic_area_scaled_tpfa_half(
            1.0,
            4.0,
            discretization::TpfaHalfConnectionProjection3D::
                positive_projection);
    const auto zero =
        synthetic_area_scaled_tpfa_half(
            1.0,
            0.0,
            discretization::TpfaHalfConnectionProjection3D::
                zero_projection);

    const auto one_zero =
        discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
            zero,
            positive_neighbour);
    const auto one_zero_swapped =
        discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
            positive_neighbour,
            zero);
    require(
        one_zero.disposition ==
            discretization::TpfaStaticFaceTransmissibilityDisposition3D::
                zero_due_to_one_half &&
            one_zero.face_transmissibility_m3 ==
                0.0,
        "one zero half blocks static face transmissibility");
    require(
        one_zero_swapped.disposition ==
            one_zero.disposition &&
            one_zero_swapped.face_transmissibility_m3 ==
                0.0,
        "one-zero semantics are owner/neighbour symmetric");

    const auto both_zero =
        discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
            zero,
            zero);
    require(
        both_zero.disposition ==
            discretization::TpfaStaticFaceTransmissibilityDisposition3D::
                zero_due_to_both_halves &&
            both_zero.face_transmissibility_m3 ==
                0.0,
        "both zero halves produce explicit both-zero state");

    expect_throw<std::invalid_argument>(
        [&] {
            auto negative =
                positive_owner;
            negative.half_transmissibility_m3 =
                -1.0;
            (void)discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
                negative,
                positive_neighbour);
        });
    expect_throw<std::invalid_argument>(
        [&] {
            auto not_finite =
                positive_owner;
            not_finite.half_transmissibility_m3 =
                std::numeric_limits<double>::
                    quiet_NaN();
            (void)discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
                not_finite,
                positive_neighbour);
        });
    expect_throw<std::invalid_argument>(
        [&] {
            auto not_finite =
                positive_owner;
            not_finite.half_transmissibility_m3 =
                std::numeric_limits<double>::
                    infinity();
            (void)discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
                not_finite,
                positive_neighbour);
        });
    expect_throw<std::invalid_argument>(
        [&] {
            auto positive_state_zero_value =
                positive_owner;
            positive_state_zero_value
                .half_transmissibility_m3 = 0.0;
            (void)discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
                positive_state_zero_value,
                positive_neighbour);
        });
    expect_throw<std::invalid_argument>(
        [&] {
            auto zero_state_positive_value =
                zero;
            zero_state_positive_value
                .half_transmissibility_m3 = 1.0;
            (void)discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
                zero_state_positive_value,
                positive_neighbour);
        });
    expect_throw<std::invalid_argument>(
        [&] {
            const auto different_area =
                synthetic_area_scaled_tpfa_half(
                    2.0,
                    4.0,
                    discretization::TpfaHalfConnectionProjection3D::
                        positive_projection);
            (void)discretization::combine_internal_face_tpfa_half_transmissibilities_3d(
                positive_owner,
                different_area);
        });

    const auto processed =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_two_cell_all_active_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto geometry =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);
    const auto permeability =
        processed_diagonal_permeability(
            processed);
    const auto& face_cells =
        processed.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);

    bool checked_boundary = false;
    for (std::size_t face = 0U;
         face < processed.topology.entity_count(
             mesh::EntityKind::face);
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        if (face_cells.adjacent(local).size() !=
            1U) {
            continue;
        }
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::make_internal_face_static_tpfa_transmissibility_3d(
                    geometry,
                    permeability,
                    local);
            });
        checked_boundary = true;
        break;
    }
    require(
        checked_boundary,
        "static TPFA invalid test found boundary face");
}


void tpfa_internal_face_transmissibility_snapshot_3d() {
    const auto horizontal =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_two_cell_all_active_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto horizontal_geometry =
        mesh::make_cell_face_geometric_operator_3d(
            horizontal.topology,
            horizontal.vertex_coordinates_m,
            horizontal.face_geometry);
    const auto horizontal_permeability =
        processed_diagonal_permeability(
            horizontal);
    const auto horizontal_interface =
        only_shared_face(
            horizontal.topology);

    const auto strict_geometry_policy =
        discretization::TransmissibilityGeometryAdmissibilityPolicy3D{
            0.0};
    const auto strict_k_policy =
        discretization::KOrthogonalityAdmissibilityPolicy3D{
            0.0};

    const auto horizontal_snapshot =
        discretization::make_admissibility_gated_internal_face_transmissibility_snapshot_3d(
            horizontal_geometry,
            horizontal_permeability,
            strict_geometry_policy,
            strict_k_policy);

    require(
        horizontal_snapshot.total_face_count() ==
                horizontal_geometry.face_count() &&
            horizontal_snapshot.internal_face_count() ==
                1U &&
            horizontal_snapshot.materialized_face_count() ==
                1U &&
            horizontal_snapshot.blocked_face_count() ==
                0U,
        "strict orthogonal snapshot materializes exactly one internal face");
    require_close(
        horizontal_snapshot
            .geometry_policy()
            .max_direct_normal_projection_angle_rad,
        0.0,
        0.0,
        "snapshot retains strict geometry policy");
    require_close(
        horizontal_snapshot
            .k_policy()
            .max_half_face_co_normal_angle_rad,
        0.0,
        0.0,
        "snapshot retains strict K policy");

    const auto& horizontal_entry =
        horizontal_snapshot.entry(
            horizontal_interface);
    require(
        horizontal_entry.disposition ==
                discretization::TpfaInternalFaceTransmissibilityDisposition3D::
                    materialized &&
            horizontal_entry.admissibility.disposition ==
                discretization::CombinedTransmissibilityAdmissibilityDisposition3D::
                    direct_normal_projection_k_orthogonal_candidate &&
            horizontal_entry
                .static_transmissibility
                .has_value(),
        "strict orthogonal internal face is materialized");
    require_close(
        horizontal_entry
            .static_transmissibility
            ->face_transmissibility_m3,
        400.0e-15 / 3.0,
        1.0e-28,
        "gated snapshot stores static face transmissibility");

    const auto& face_cells =
        horizontal.topology.relation(
            mesh::EntityKind::face,
            mesh::EntityKind::cell);
    bool checked_boundary = false;
    for (std::size_t face = 0U;
         face < horizontal_geometry.face_count();
         ++face) {
        const auto local =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        face)};
        if (face_cells.adjacent(local).size() !=
            1U) {
            continue;
        }
        require(
            !horizontal_snapshot
                .contains_internal_face(local),
            "boundary face is absent from internal-face snapshot");
        expect_throw<std::invalid_argument>(
            [&] {
                (void)horizontal_snapshot.entry(
                    local);
            });
        checked_boundary = true;
        break;
    }
    require(
        checked_boundary,
        "gated snapshot test found boundary face");

    const auto skewed =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_skewed_two_cell_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto skewed_geometry =
        mesh::make_cell_face_geometric_operator_3d(
            skewed.topology,
            skewed.vertex_coordinates_m,
            skewed.face_geometry);
    const auto skewed_permeability =
        processed_diagonal_permeability(
            skewed);
    const auto skewed_interface =
        only_shared_face(
            skewed.topology);
    const double geometry_angle =
        std::acos(
            5.0 / std::sqrt(29.0));
    const double owner_k_angle =
        std::atan(0.2);

    const auto blocked_both =
        discretization::make_admissibility_gated_internal_face_transmissibility_snapshot_3d(
            skewed_geometry,
            skewed_permeability,
            strict_geometry_policy,
            strict_k_policy);
    const auto& blocked_both_entry =
        blocked_both.entry(
            skewed_interface);
    require(
        blocked_both.materialized_face_count() ==
                0U &&
            blocked_both.blocked_face_count() ==
                1U &&
            blocked_both_entry.disposition ==
                discretization::TpfaInternalFaceTransmissibilityDisposition3D::
                    blocked_geometry_and_k_non_orthogonal &&
            !blocked_both_entry
                 .static_transmissibility
                 .has_value(),
        "strict skewed face is blocked without materializing T_f");

    const auto blocked_k =
        discretization::make_admissibility_gated_internal_face_transmissibility_snapshot_3d(
            skewed_geometry,
            skewed_permeability,
            discretization::TransmissibilityGeometryAdmissibilityPolicy3D{
                geometry_angle + 1.0e-12},
            strict_k_policy);
    require(
        blocked_k.entry(
            skewed_interface)
                .disposition ==
            discretization::TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_k_non_orthogonal &&
            !blocked_k.entry(
                skewed_interface)
                 .static_transmissibility
                 .has_value(),
        "relaxed geometry policy leaves K-only block");

    const auto blocked_geometry =
        discretization::make_admissibility_gated_internal_face_transmissibility_snapshot_3d(
            skewed_geometry,
            skewed_permeability,
            strict_geometry_policy,
            discretization::KOrthogonalityAdmissibilityPolicy3D{
                owner_k_angle + 1.0e-12});
    require(
        blocked_geometry.entry(
            skewed_interface)
                .disposition ==
            discretization::TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_geometry_non_orthogonal &&
            !blocked_geometry.entry(
                skewed_interface)
                 .static_transmissibility
                 .has_value(),
        "relaxed K policy leaves geometry-only block");

    const auto relaxed_both =
        discretization::make_admissibility_gated_internal_face_transmissibility_snapshot_3d(
            skewed_geometry,
            skewed_permeability,
            discretization::TransmissibilityGeometryAdmissibilityPolicy3D{
                geometry_angle + 1.0e-12},
            discretization::KOrthogonalityAdmissibilityPolicy3D{
                owner_k_angle + 1.0e-12});
    const auto& relaxed_entry =
        relaxed_both.entry(
            skewed_interface);
    require(
        relaxed_entry.disposition ==
                discretization::TpfaInternalFaceTransmissibilityDisposition3D::
                    materialized &&
            relaxed_entry.static_transmissibility.has_value(),
        "explicitly relaxed geometry/K policies materialize skewed face");
    require_close(
        relaxed_entry
            .static_transmissibility
            ->face_transmissibility_m3,
        400.0e-15 / 3.0,
        1.0e-28,
        "relaxed gated snapshot stores skewed static transmissibility");

    const auto zero_permeability =
        mesh::CellCartesianDiagonalPermeability3D{
            {
                mesh::CartesianDiagonalPermeabilityTensor3D{
                    0.0, 0.0, 0.0},
                mesh::CartesianDiagonalPermeabilityTensor3D{
                    0.0, 0.0, 0.0},
            }};
    const auto degenerate_snapshot =
        discretization::make_admissibility_gated_internal_face_transmissibility_snapshot_3d(
            horizontal_geometry,
            zero_permeability,
            strict_geometry_policy,
            strict_k_policy);
    const auto& degenerate_entry =
        degenerate_snapshot.entry(
            horizontal_interface);
    require(
        degenerate_entry.disposition ==
                discretization::TpfaInternalFaceTransmissibilityDisposition3D::
                    blocked_degenerate_permeability_direction &&
            !degenerate_entry
                 .static_transmissibility
                 .has_value(),
        "degenerate permeability blocks transmissibility materialization");
}

void tpfa_internal_face_transmissibility_snapshot_3d_invalid() {
    const auto processed =
        mesh::process_active_corner_point_grid(
            mesh::import_grdecl(
                grdecl_two_cell_all_active_fixture(),
                mesh::GrdeclImportOptions{
                    1.0, 1.0e-15}));
    const auto geometry =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);
    const auto permeability =
        processed_diagonal_permeability(
            processed);
    const auto interface =
        only_shared_face(
            processed.topology);

    const auto valid =
        discretization::make_admissibility_gated_internal_face_transmissibility_snapshot_3d(
            geometry,
            permeability,
            discretization::TransmissibilityGeometryAdmissibilityPolicy3D{
                0.0},
            discretization::KOrthogonalityAdmissibilityPolicy3D{
                0.0});
    const auto valid_entry =
        valid.entry(interface);

    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::TpfaInternalFaceTransmissibilitySnapshot3D{
                geometry,
                discretization::TransmissibilityGeometryAdmissibilityPolicy3D{
                    -1.0},
                discretization::KOrthogonalityAdmissibilityPolicy3D{
                    0.0},
                {}};
        });
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::TpfaInternalFaceTransmissibilitySnapshot3D{
                geometry,
                discretization::TransmissibilityGeometryAdmissibilityPolicy3D{
                    0.0},
                discretization::KOrthogonalityAdmissibilityPolicy3D{
                    std::numeric_limits<double>::
                        quiet_NaN()},
                {}};
        });

    {
        auto duplicate =
            std::vector<
                discretization::TpfaInternalFaceTransmissibilityEntry3D>{
                valid_entry,
                valid_entry};
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::TpfaInternalFaceTransmissibilitySnapshot3D{
                geometry,
                    valid.geometry_policy(),
                    valid.k_policy(),
                    std::move(duplicate)};
            });
    }

    {
        auto missing =
            valid_entry;
        missing.static_transmissibility =
            std::nullopt;
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::TpfaInternalFaceTransmissibilitySnapshot3D{
                geometry,
                    valid.geometry_policy(),
                    valid.k_policy(),
                    {missing}};
            });
    }

    {
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::TpfaInternalFaceTransmissibilitySnapshot3D{
                    geometry,
                    valid.geometry_policy(),
                    valid.k_policy(),
                    {}};
            });
    }

    {
        auto mismatched_policy =
            valid_entry;
        mismatched_policy
            .admissibility.geometry
            .max_direct_normal_projection_angle_rad =
            1.0e-6;
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::TpfaInternalFaceTransmissibilitySnapshot3D{
                    geometry,
                    valid.geometry_policy(),
                    valid.k_policy(),
                    {mismatched_policy}};
            });
    }

    {
        auto non_finite =
            valid_entry;
        non_finite
            .static_transmissibility
            ->face_transmissibility_m3 =
            std::numeric_limits<double>::
                quiet_NaN();
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::TpfaInternalFaceTransmissibilitySnapshot3D{
                geometry,
                    valid.geometry_policy(),
                    valid.k_policy(),
                    {non_finite}};
            });
    }

    {
        auto blocked =
            valid_entry;
        blocked.disposition =
            discretization::TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_geometry_non_orthogonal;
        blocked.admissibility.disposition =
            discretization::CombinedTransmissibilityAdmissibilityDisposition3D::
                requires_geometry_non_orthogonal_treatment;
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::TpfaInternalFaceTransmissibilitySnapshot3D{
                geometry,
                    valid.geometry_policy(),
                    valid.k_policy(),
                    {blocked}};
            });
    }

    const auto mismatched_permeability =
        mesh::CellCartesianDiagonalPermeability3D{
            {mesh::CartesianDiagonalPermeabilityTensor3D{
                100.0e-15,
                50.0e-15,
                10.0e-15}}};
    expect_throw<std::invalid_argument>(
        [&] {
            (void)discretization::make_admissibility_gated_internal_face_transmissibility_snapshot_3d(
                geometry,
                mismatched_permeability,
                valid.geometry_policy(),
                valid.k_policy());
        });

    expect_throw<std::out_of_range>(
        [&] {
            (void)valid.contains_internal_face(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            geometry.face_count())});
        });
}


} // namespace

int main(int argc, char** argv) {
    try {
        require(
            argc == 2,
            "provide one named discretization TPFA test");
        const std::string_view name{argv[1]};
        if (name == "tpfa_half_connection_3d") {
            tpfa_half_connection_3d();
        } else if (name == "tpfa_half_connection_3d_invalid") {
            tpfa_half_connection_3d_invalid();
        } else if (name == "tpfa_half_transmissibility_3d") {
            tpfa_half_transmissibility_3d();
        } else if (name == "tpfa_half_transmissibility_3d_invalid") {
            tpfa_half_transmissibility_3d_invalid();
        } else if (name == "tpfa_static_face_transmissibility_3d") {
            tpfa_static_face_transmissibility_3d();
        } else if (name == "tpfa_static_face_transmissibility_3d_invalid") {
            tpfa_static_face_transmissibility_3d_invalid();
        } else if (name == "tpfa_internal_face_transmissibility_snapshot_3d") {
            tpfa_internal_face_transmissibility_snapshot_3d();
        } else if (name == "tpfa_internal_face_transmissibility_snapshot_3d_invalid") {
            tpfa_internal_face_transmissibility_snapshot_3d_invalid();
        } else {
            throw std::invalid_argument(
                "unknown discretization TPFA test");
        }
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
