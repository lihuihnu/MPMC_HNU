#include <mpmc/mesh/linear_cell_mesh_2d.hpp>

void linear_cell_mesh_2d_header() {
    using namespace mpmc::mesh;
    const std::array cells{LinearCell2D{GlobalEntityId{9}, LinearCellType2D::triangle,
        {LocalIndex{0}, LocalIndex{1}, LocalIndex{2}}}};
    const auto mesh = make_linear_mesh_2d(
        {GlobalEntityId{1}, GlobalEntityId{2}, GlobalEntityId{3}},
        {{0,0}, {1,0}, {0,1}}, cells);
    if (mesh.geometry.cell_area_m2(LocalIndex{0}) != 0.5) {
        throw std::runtime_error("standalone linear 2D header probe");
    }
}
