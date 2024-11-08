// Driver code for the first testcase for the obtuse angle fix

#include <iostream>
#include <cmath>

// OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Local include files
#include "parameters.h"
#include "singular_navier_stokes_elements.h"
#include "singular_sector_problem.h"
#include "my_element.h"

using namespace oomph;

// Start of the main function
int main(int argc, char** argv)
{
#ifdef OOMPH_HAS_MPI
  // Setup mpi but don't make a copy of mpi_comm_world because
  // mumps wants to work with the real thing.
  bool make_copy_of_mpi_comm_world = false;
  MPI_Helpers::init(argc, argv, make_copy_of_mpi_comm_world);
#endif

  oomph_info << "obtuse_corner_flow" << std::endl;

  // Problem parameters
  Parameters parameters;

  // Create problem
  SingularSectorProblem<SingularNavierStokesElement<MyElement>> problem(140,
                                                                        15);

  // Steady problem
  problem.steady_newton_solve();
  problem.doc_solution();

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

  // Return 0 to tell everyone that the program finished successfully
  return 0;
}
