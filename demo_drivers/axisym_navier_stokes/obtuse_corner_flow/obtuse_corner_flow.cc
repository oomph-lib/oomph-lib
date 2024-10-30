// Driver code for the first testcase for the obtuse angle fix

#include <iostream>
#include <cmath>

#include "generic.h"
#include "navier_stokes.h"
#include "axisym_navier_stokes.h"

#include "axisym_navier_stokes.h"
#include "axisym_navier_stokes_with_singularity.h"
#include "my_navier_stokes_elements_with_singularity.h"

#include "axisym_sector_problem.h"
#include "parameter_values.h"

#include "my_navier_stokes_elements_with_singularity.h"

#include "parse_arguments.h"

using namespace oomph;

// Start of the main function
int main(int argc, char** argv)
{
  std::cout << "Test case 1" << std::endl;
  Arguments passed_arguments(argc, argv);
  bool debug = passed_arguments.has_debug();

  // Set parameter values
  parameters::x_centre[0] = 2.0;

  // Create problem
  AxisymSectorProblem<AxisymNavierStokesElementWithSingularity<
    Hijacked<ProjectableAxisymmetricTaylorHoodElement<
      PseudoSolidNodeUpdateElement<AxisymmetricTTaylorHoodElement,
                                   TPVDElement<2, 3>>>>>>
    problem;

  if (debug)
  {
    problem.debug_residuals();
    problem.debug_jacobian();
  }
  else
  {
    // Unsteady problem
    const double dt = 2e-2;
    const double ft = 5 * dt;
    const unsigned nt = std::ceil(ft / dt);

    for (unsigned it = 0; it < nt; it++)
    {
      std::cout << "unsteady_newton_solve" << std::endl;
      problem.unsteady_newton_solve(dt);
      problem.doc_solution();
    }

    // Steady problem
    problem.make_steady();
    problem.steady_newton_solve();
    problem.doc_solution();
  }

  return 0;
}
