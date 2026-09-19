#include <mpmc/discretization/transmissibility_admissibility_3d.hpp>
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
        active_corner_point_field(processed, "PERMX"),
        active_corner_point_field(processed, "PERMY"),
        active_corner_point_field(processed, "PERMZ"));
}

mesh::ActiveCornerPointGrid
processed_fixture(const std::string& deck) {
    return mesh::process_active_corner_point_grid(
        mesh::import_grdecl(
            deck,
            mesh::GrdeclImportOptions{
                1.0, 1.0e-15}));
}

void geometry_admissibility_orthogonal() {
    const auto processed =
        processed_fixture(
            grdecl_two_cell_all_active_fixture());
    const auto geometry =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);
    const auto interface =
        only_shared_face(processed.topology);

    const auto strict =
        discretization::
            classify_internal_face_transmissibility_geometry(
                geometry,
                interface,
                discretization::
                    TransmissibilityGeometryAdmissibilityPolicy3D{
                        0.0});
    require(
        strict.disposition ==
            discretization::
                TransmissibilityGeometryDisposition3D::
                    direct_normal_projection_allowed,
        "orthogonal face allows direct projection under strict policy");
    require_close(
        strict.non_orthogonality_angle_rad,
        0.0,
        1.0e-14,
        "orthogonal admissibility angle");
    require_close(
        strict.max_direct_normal_projection_angle_rad,
        0.0,
        0.0,
        "orthogonal strict policy retained");

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
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::
                    classify_internal_face_transmissibility_geometry(
                        geometry,
                        local,
                        discretization::
                            TransmissibilityGeometryAdmissibilityPolicy3D{
                                0.0});
            });
        checked_boundary = true;
        break;
    }
    require(
        checked_boundary,
        "orthogonal admissibility fixture exposes boundary face");
}

void geometry_admissibility_skewed() {
    const auto processed =
        processed_fixture(
            grdecl_skewed_two_cell_fixture());
    const auto geometry =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);
    const auto interface =
        only_shared_face(processed.topology);
    const double expected_angle =
        std::acos(5.0 / std::sqrt(29.0));

    const auto strict =
        discretization::
            classify_internal_face_transmissibility_geometry(
                geometry,
                interface,
                discretization::
                    TransmissibilityGeometryAdmissibilityPolicy3D{
                        0.0});
    require(
        strict.disposition ==
            discretization::
                TransmissibilityGeometryDisposition3D::
                    requires_non_orthogonal_treatment,
        "skewed face requires non-orthogonal treatment under strict policy");
    require_close(
        strict.non_orthogonality_angle_rad,
        expected_angle,
        1.0e-14,
        "skewed geometry classifier reports actual angle");

    const auto relaxed_policy =
        discretization::
            TransmissibilityGeometryAdmissibilityPolicy3D{
                expected_angle + 1.0e-12};
    const auto relaxed =
        discretization::
            classify_internal_face_transmissibility_geometry(
                geometry,
                interface,
                relaxed_policy);
    require(
        relaxed.disposition ==
            discretization::
                TransmissibilityGeometryDisposition3D::
                    direct_normal_projection_allowed,
        "explicit relaxed geometry policy allows direct projection");
    require_close(
        relaxed.max_direct_normal_projection_angle_rad,
        relaxed_policy
            .max_direct_normal_projection_angle_rad,
        0.0,
        "relaxed geometry policy retained");

    for (const double invalid_policy : {
             -1.0e-6,
             std::numeric_limits<double>::
                 quiet_NaN(),
             0.5 * std::acos(-1.0)}) {
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::
                    classify_internal_face_transmissibility_geometry(
                        geometry,
                        interface,
                        discretization::
                            TransmissibilityGeometryAdmissibilityPolicy3D{
                                invalid_policy});
            });
    }
}

void require_strict_axis_aligned(
    const mesh::ActiveCornerPointGrid& processed) {
    const auto geometry =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);
    const auto permeability =
        processed_diagonal_permeability(processed);
    const auto interface =
        only_shared_face(processed.topology);

    const auto k_result =
        discretization::
            classify_internal_face_k_orthogonality(
                geometry,
                permeability,
                interface,
                discretization::
                    KOrthogonalityAdmissibilityPolicy3D{
                        0.0});
    require(
        k_result.disposition ==
            discretization::
                KOrthogonalityDisposition3D::
                    k_orthogonal_within_policy,
        "axis-aligned fixture is K-orthogonal under strict policy");
    require(
        k_result.owner_half_face.has_value() &&
            k_result.neighbour_half_face.has_value(),
        "axis-aligned half-face K directions exist");
    require_close(
        k_result.owner_half_face->angle_rad,
        0.0,
        1.0e-14,
        "axis-aligned owner half-face K angle");
    require_close(
        k_result.neighbour_half_face->angle_rad,
        0.0,
        1.0e-14,
        "axis-aligned neighbour half-face K angle");
    require(
        k_result.center_line_diagnostic.owner.has_value() &&
            k_result.center_line_diagnostic.neighbour.has_value(),
        "axis-aligned K*d_cc diagnostics exist");
    require_close(
        k_result.center_line_diagnostic.owner->angle_rad,
        0.0,
        1.0e-14,
        "axis-aligned owner K*d_cc diagnostic angle");
    require_close(
        k_result.center_line_diagnostic.neighbour->angle_rad,
        0.0,
        1.0e-14,
        "axis-aligned neighbour K*d_cc diagnostic angle");

    const auto combined =
        discretization::
            classify_internal_face_transmissibility_admissibility(
                geometry,
                permeability,
                interface,
                discretization::
                    TransmissibilityGeometryAdmissibilityPolicy3D{
                        0.0},
                discretization::
                    KOrthogonalityAdmissibilityPolicy3D{
                        0.0});
    require(
        combined.disposition ==
            discretization::
                CombinedTransmissibilityAdmissibilityDisposition3D::
                    direct_normal_projection_k_orthogonal_candidate,
        "axis-aligned geometry and K combine to direct candidate");
}

void k_orthogonality_axis_aligned() {
    require_strict_axis_aligned(
        processed_fixture(
            grdecl_two_cell_all_active_fixture()));
    require_strict_axis_aligned(
        processed_fixture(
            grdecl_vertical_two_cell_fixture()));
}

void k_orthogonality_skewed() {
    const auto processed =
        processed_fixture(
            grdecl_skewed_two_cell_fixture());
    const auto geometry =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);
    const auto permeability =
        processed_diagonal_permeability(processed);
    const auto interface =
        only_shared_face(processed.topology);

    const double geometry_angle =
        std::acos(5.0 / std::sqrt(29.0));
    const double owner_k_angle =
        std::atan(0.2);
    const double neighbour_k_angle =
        std::atan(0.15);

    const auto strict_k =
        discretization::
            classify_internal_face_k_orthogonality(
                geometry,
                permeability,
                interface,
                discretization::
                    KOrthogonalityAdmissibilityPolicy3D{
                        0.0});
    require(
        strict_k.disposition ==
            discretization::
                KOrthogonalityDisposition3D::
                    requires_k_non_orthogonal_treatment,
        "skewed anisotropic fixture fails strict K-orthogonality");
    require(
        strict_k.owner_half_face.has_value() &&
            strict_k.neighbour_half_face.has_value(),
        "skewed half-face K directions exist");
    require_close(
        strict_k.owner_half_face->angle_rad,
        owner_k_angle,
        1.0e-14,
        "skewed owner half-face co-normal angle");
    require_close(
        strict_k.neighbour_half_face->angle_rad,
        neighbour_k_angle,
        1.0e-14,
        "skewed neighbour half-face co-normal angle");
    require(
        strict_k.center_line_diagnostic.owner.has_value() &&
            strict_k.center_line_diagnostic.neighbour.has_value(),
        "skewed K*d_cc diagnostics exist");
    require_close(
        strict_k.center_line_diagnostic.owner->angle_rad,
        geometry_angle,
        1.0e-14,
        "skewed owner K*d_cc diagnostic angle");
    require_close(
        strict_k.center_line_diagnostic.neighbour->angle_rad,
        geometry_angle,
        1.0e-14,
        "skewed neighbour K*d_cc diagnostic angle");
    require(
        std::abs(
            strict_k.center_line_diagnostic.owner->angle_rad -
            strict_k.owner_half_face->angle_rad) > 0.1,
        "K*d_cc diagnostic must not alias half-face K-orthogonality");

    const auto strict_both =
        discretization::
            classify_internal_face_transmissibility_admissibility(
                geometry,
                permeability,
                interface,
                discretization::
                    TransmissibilityGeometryAdmissibilityPolicy3D{
                        0.0},
                discretization::
                    KOrthogonalityAdmissibilityPolicy3D{
                        0.0});
    require(
        strict_both.disposition ==
            discretization::
                CombinedTransmissibilityAdmissibilityDisposition3D::
                    requires_geometry_and_k_non_orthogonal_treatment,
        "strict skewed fixture requires geometry and K treatment");

    const auto geometry_relaxed =
        discretization::
            classify_internal_face_transmissibility_admissibility(
                geometry,
                permeability,
                interface,
                discretization::
                    TransmissibilityGeometryAdmissibilityPolicy3D{
                        geometry_angle + 1.0e-12},
                discretization::
                    KOrthogonalityAdmissibilityPolicy3D{
                        0.0});
    require(
        geometry_relaxed.disposition ==
            discretization::
                CombinedTransmissibilityAdmissibilityDisposition3D::
                    requires_k_non_orthogonal_treatment,
        "relaxed geometry isolates K non-orthogonality");

    const auto k_relaxed =
        discretization::
            classify_internal_face_transmissibility_admissibility(
                geometry,
                permeability,
                interface,
                discretization::
                    TransmissibilityGeometryAdmissibilityPolicy3D{
                        0.0},
                discretization::
                    KOrthogonalityAdmissibilityPolicy3D{
                        owner_k_angle + 1.0e-12});
    require(
        k_relaxed.disposition ==
            discretization::
                CombinedTransmissibilityAdmissibilityDisposition3D::
                    requires_geometry_non_orthogonal_treatment,
        "relaxed K policy isolates geometry non-orthogonality");

    const auto both_relaxed =
        discretization::
            classify_internal_face_transmissibility_admissibility(
                geometry,
                permeability,
                interface,
                discretization::
                    TransmissibilityGeometryAdmissibilityPolicy3D{
                        geometry_angle + 1.0e-12},
                discretization::
                    KOrthogonalityAdmissibilityPolicy3D{
                        owner_k_angle + 1.0e-12});
    require(
        both_relaxed.disposition ==
            discretization::
                CombinedTransmissibilityAdmissibilityDisposition3D::
                    direct_normal_projection_k_orthogonal_candidate,
        "explicit relaxed geometry and K policies produce candidate");
}

void k_orthogonality_invalid() {
    const auto processed =
        processed_fixture(
            grdecl_two_cell_all_active_fixture());
    const auto geometry =
        mesh::make_cell_face_geometric_operator_3d(
            processed.topology,
            processed.vertex_coordinates_m,
            processed.face_geometry);
    const auto permeability =
        processed_diagonal_permeability(processed);
    const auto interface =
        only_shared_face(processed.topology);

    for (const double invalid_policy : {
             -1.0e-6,
             std::numeric_limits<double>::
                 quiet_NaN(),
             0.5 * std::acos(-1.0)}) {
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::
                    classify_internal_face_k_orthogonality(
                        geometry,
                        permeability,
                        interface,
                        discretization::
                            KOrthogonalityAdmissibilityPolicy3D{
                                invalid_policy});
            });
    }

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
        expect_throw<std::invalid_argument>(
            [&] {
                (void)discretization::
                    classify_internal_face_k_orthogonality(
                        geometry,
                        permeability,
                        local,
                        discretization::
                            KOrthogonalityAdmissibilityPolicy3D{
                                0.0});
            });
        checked_boundary = true;
        break;
    }
    require(
        checked_boundary,
        "K-orthogonality invalid fixture exposes boundary face");

    const mesh::CellCartesianDiagonalPermeability3D
        degenerate{{
            mesh::CartesianDiagonalPermeabilityTensor3D{
                0.0, 0.0, 0.0},
            mesh::CartesianDiagonalPermeabilityTensor3D{
                0.0, 0.0, 0.0}}};

    const auto degenerate_k =
        discretization::
            classify_internal_face_k_orthogonality(
                geometry,
                degenerate,
                interface,
                discretization::
                    KOrthogonalityAdmissibilityPolicy3D{
                        0.0});
    require(
        degenerate_k.disposition ==
            discretization::
                KOrthogonalityDisposition3D::
                    degenerate_permeability_direction,
        "zero permeability produces explicit degenerate K direction");

    const auto combined =
        discretization::
            classify_internal_face_transmissibility_admissibility(
                geometry,
                degenerate,
                interface,
                discretization::
                    TransmissibilityGeometryAdmissibilityPolicy3D{
                        0.0},
                discretization::
                    KOrthogonalityAdmissibilityPolicy3D{
                        0.0});
    require(
        combined.disposition ==
            discretization::
                CombinedTransmissibilityAdmissibilityDisposition3D::
                    degenerate_permeability_direction,
        "combined admissibility preserves degenerate permeability state");
}

} // namespace

int main(int argc, char** argv) {
    try {
        require(
            argc == 2,
            "provide one named discretization admissibility test");
        const std::string_view name{argv[1]};
        if (name == "geometry_admissibility_orthogonal") {
            geometry_admissibility_orthogonal();
        } else if (name == "geometry_admissibility_skewed") {
            geometry_admissibility_skewed();
        } else if (name == "k_orthogonality_axis_aligned") {
            k_orthogonality_axis_aligned();
        } else if (name == "k_orthogonality_skewed") {
            k_orthogonality_skewed();
        } else if (name == "k_orthogonality_invalid") {
            k_orthogonality_invalid();
        } else {
            throw std::invalid_argument(
                "unknown discretization admissibility test");
        }
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
