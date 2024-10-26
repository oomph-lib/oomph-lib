#ifndef RUN_TESTS_HEADER
#define RUN_TESTS_HEADER

#include <iostream>

#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Local include files
#include "singular_axisym_dynamic_cap_problem.h"

#include "singular_axisym_navier_stokes_elements.h"
#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"

#include "parameters.h"

#include "utility_functions.h"

namespace oomph
{
  typedef SingularAxisymNavierStokesElement<
    HijackedProjectableAxisymmetricTTaylorHoodPVDElement>
    BASE_ELEMENT;
  typedef BDF<2> TIMESTEPPER;
  typedef SingularAxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER>
    AXISYM_PROBLEM;

  std::shared_ptr<AXISYM_PROBLEM> createBaseProblem(
    const double& contact_angle = 90.0 * oomph::MathematicalConstants::Pi /
                                  180.0);

  enum
  {
    upper,
    outer,
    lower,
    inner,
  };

} // namespace oomph

#endif
