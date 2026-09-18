#include <mpmc/mesh_petsc/adapter.hpp>

#include <type_traits>

static_assert(std::is_same_v<
              decltype(&mpmc::mesh_petsc::create_section_mapping),
              PetscErrorCode (*)(MPI_Comm,
                                 const mpmc::mesh::DofLayout&,
                                 const mpmc::mesh::DofNumberingSnapshot&,
                                 PetscSection*,
                                 std::vector<PetscInt>*)>);
