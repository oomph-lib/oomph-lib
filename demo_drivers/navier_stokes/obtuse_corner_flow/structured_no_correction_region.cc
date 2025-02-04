// Driver code for the first testcase for the obtuse angle fix

#include <iostream>
#include <cmath>

// OOMPH-LIB include files
#include "generic.h"
#include "navier_stokes.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Local include files
#include "region_sector_problem.h"
#include "my_element.h"

using namespace oomph;

// Start of the main function
int main(int argc, char** argv)
{
  oomph_info << "structured_no_correction" << std::endl;

  // Create problem
  RegionSectorProblem<MyElement> problem;
  problem.setup();

  // Steady problem
  problem.steady_newton_solve();
  problem.doc_solution();

  // Return 0 to tell everyone that the program finished successfully
  return 0;
}
