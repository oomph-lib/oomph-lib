// Driver code for the first testcase for the obtuse angle fix

#include <iostream>
#include <cmath>

// OOMPH-LIB include files
#include "generic.h"

// Local include files
#include "region_axisym_sector_problem.h"
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

  oomph_info << "structured_no_correction_region" << std::endl;

  // Create problem
  Params params = create_parameters_from_file("parameters.dat");
  RegionAxisymSectorProblem<MyElement> problem(params);
  problem.setup();

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
