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
#include "axisym_dynamic_cap_problem.h"

#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"

#include "parameters.h"

#include "utility_functions.h"

namespace oomph
{
  typedef HijackedProjectableAxisymmetricTTaylorHoodPVDElement BASE_ELEMENT;
  typedef BDF<2> TIMESTEPPER;
  typedef AxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER> AXISYM_PROBLEM;

  std::shared_ptr<AXISYM_PROBLEM> createBaseProblem();

  enum
  {
    upper,
    outer,
    lower,
    inner,
  };

} // namespace oomph

#endif
