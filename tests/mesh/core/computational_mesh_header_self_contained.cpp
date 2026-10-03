#include <mpmc/mesh/computational_mesh.hpp>

// Compiled with only the public header and exercised by the dedicated test.
void computational_mesh_header_probe() {
    using namespace mpmc::mesh;
    const std::vector<LinearCell2D> cells{{GlobalEntityId{9},LinearCellType2D::triangle,
                                        {LocalIndex{0},LocalIndex{1},LocalIndex{2}}}};
    const auto grid = make_linear_mesh_2d({GlobalEntityId{3},GlobalEntityId{4},GlobalEntityId{5}},
                                        {{0,0},{1,0},{0,1}},cells);
    const auto document = MeshExchangeDocument::create(MeshExchangeFormat::vtu_ascii,2,
        grid.topology,{{0,0,0},{1,0,0},{0,1,0}},grid.face_boundary,{}, {},std::nullopt);
    if (prepare_linear_mesh_2d(document).geometry.cell_area_m2(LocalIndex{0}) != .5) {
        throw std::runtime_error("computational mesh public header probe");
    }
}
