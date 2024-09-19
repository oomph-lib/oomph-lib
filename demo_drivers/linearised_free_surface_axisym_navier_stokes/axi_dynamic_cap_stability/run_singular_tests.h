#ifndef RUN_SINGULAR_TESTS_HEADER
#define RUN_SINGULAR_TESTS_HEADER

#include <iostream>

#include "generic.h"
#include "axisym_navier_stokes.h"
#include "singular_axisym_navier_stokes_elements.h"
#include "bo_height_control_singular_axisym_dynamic_cap_problem.h"
#include "ca_height_control_singular_axisym_dynamic_cap_problem.h"
#include "linearised_axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Local include files
#include "axisym_dynamic_cap_problem.h"
#include "singular_axisym_dynamic_cap_problem.h"
#include "perturbed_linear_stability_cap_problem.h"

#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"
#include "overlaying_my_linear_element.h"
#include "overlaying_Tlinear_axisym_ns_pvd_elements.h"

#include "parameters.h"

#include "utility_functions.h"

namespace oomph
{
  typedef BDF<2> TIMESTEPPER;

  typedef HijackedProjectableAxisymmetricTTaylorHoodPVDElement BASE_ELEMENT;
  typedef OverlayingMyLinearElement<BASE_ELEMENT> PERTURBED_BASE_ELEMENT;
  typedef AxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER> AXISYM_PROBLEM;
  typedef PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                             PERTURBED_BASE_ELEMENT,
                                             TIMESTEPPER>
    PERTURBED_BASE_PROBLEM;

  typedef SingularAxisymNavierStokesElement<
    HijackedProjectableAxisymmetricTTaylorHoodPVDElement>
    SINGULAR_ELEMENT;
  typedef OverlayingMyLinearElement<SINGULAR_ELEMENT> PERTURBED_ELEMENT;
  typedef SingularAxisymDynamicCapProblem<SINGULAR_ELEMENT, TIMESTEPPER>
    SINGULAR_PROBLEM;
  typedef PerturbedLinearStabilityCapProblem<SINGULAR_ELEMENT,
                                             PERTURBED_ELEMENT,
                                             TIMESTEPPER>
    PERTURBED_PROBLEM;

  std::shared_ptr<AXISYM_PROBLEM> createBaseProblem();
  std::shared_ptr<PERTURBED_BASE_PROBLEM> createLinearBaseProblem(
    std::shared_ptr<AXISYM_PROBLEM> base_problem_pt,
    const unsigned& azimuthal_mode_number = 0);

  std::shared_ptr<SINGULAR_PROBLEM> createSingularProblem();
  std::shared_ptr<PERTURBED_PROBLEM> createLinearProblem(
    std::shared_ptr<SINGULAR_PROBLEM> base_problem_pt,
    const unsigned& azimuthal_mode_number = 0);

  enum
  {
    upper,
    outer,
    lower,
    inner,
  };

} // namespace oomph

#endif
