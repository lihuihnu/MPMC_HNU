#ifndef MPMC_MESH_FACE_MESH_HDF5_HPP
#define MPMC_MESH_FACE_MESH_HDF5_HPP
#include <mpmc/mesh/face_mesh.hpp>
#include <filesystem>
namespace mpmc::mesh {
/// Optional mpmc::mesh_hdf5 adapter; core headers expose no HDF5 types.
/// Schema v1, strict typed datasets, no external links or implicit conversions.
/// Existing output files are refused. No proprietary MATLAB runtime required.
void write_face_mesh_hdf5(const std::filesystem::path& file,const FaceMesh& mesh);
[[nodiscard]] FaceMesh read_face_mesh_hdf5(const std::filesystem::path& file,
                                        std::uint64_t maximum_bytes=8ULL*1024*1024*1024);
}
#endif
