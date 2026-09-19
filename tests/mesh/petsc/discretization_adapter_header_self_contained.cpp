#include <mpmc/discretization_petsc/adapter.hpp>

#include <type_traits>

static_assert(std::is_same_v<
              decltype(&mpmc::discretization_petsc::
                           create_empty_petsc_mpiaij_symbolic_matrix_3d),
              PetscErrorCode (*)(
                  MPI_Comm,
                  const mpmc::discretization_petsc::
                      PetscMpiAijSymbolicPreallocation3D&,
                  Mat*)>);
