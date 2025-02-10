// Driver code for the first testcase for the obtuse angle fix

#include <iostream>
#include <fstream>

// OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Local include files
#include "singular_region_axisym_sector_problem.h"
#include "my_element.h"
#include "utility_functions.h"

using namespace oomph;

// Start of the main function
int main(int argc, char** argv)
{
  oomph_info << "structured_with_correction_regions" << std::endl;

  // Create problem
  Params params = create_parameters_from_file("parameters.dat");
  SingularRegionAxisymSectorProblem<
    SingularAxisymNavierStokesElement<MyElement>>
    problem(params);
  problem.setup();

  // Steady problem
  problem.steady_newton_solve();
  problem.doc_solution();

  // Return 0 to tell everyone that the program finished successfully
  return 0;
}
