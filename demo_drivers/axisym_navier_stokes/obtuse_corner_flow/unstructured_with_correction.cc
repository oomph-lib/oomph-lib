// Driver code for the first testcase for the obtuse angle fix

#include <iostream>
#include <fstream>

// OOMPH-LIB include files
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Local include files
#include "singular_unstructured_axisym_sector_problem.h"
#include "my_element.h"
#include "utility_functions.h"

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

  oomph_info << "unstructured_with_correction_regions" << std::endl;

  // Create problem
  Params params = create_parameters_from_file("parameters.dat");
  SingularUnstructuredAxisymSectorProblem<SingularAxisymNavierStokesElement<
    ProjectableAxisymmetricTaylorHoodElement<MyElement>>>
    problem(params);
  problem.setup();

  // Steady problem
  problem.steady_newton_solve(1);
  problem.doc_solution();

  // Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

  // Return 0 to tell everyone that the program finished successfully
  return 0;
}
