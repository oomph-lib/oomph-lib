#ifndef RUN_TESTS_HEADER
#define RUN_TESTS_HEADER

#include <iostream>

#include "generic.h"
#include "axisym_navier_stokes.h"
#include "linearised_axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Local include files
#include "singular_axisym_dynamic_cap_problem.h"
#include "perturbed_linear_stability_cap_problem.h"

#include "projectable_axisymmetric_Ttaylor_hood_elements.h"
#include "singular_axisym_navier_stokes_elements.h"
#include "overlaying_my_linear_element.h"
#include "overlaying_Tlinear_axisym_ns_pvd_elements.h"

#include "parameters.h"

#include "utility_functions.h"

namespace oomph
{
  typedef SingularAxisymNavierStokesElement<
    ProjectableAxisymmetricTTaylorHoodPVDElement>
    BASE_ELEMENT;
  // typedef OverlayingTLinearisedAxisymNSPVDElement PERTURBED_ELEMENT;
  typedef OverlayingMyLinearElement<BASE_ELEMENT> PERTURBED_ELEMENT;
  typedef BDF<2> TIMESTEPPER;
  typedef SingularAxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER>
    AXISYM_PROBLEM;
  typedef PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                             PERTURBED_ELEMENT,
                                             TIMESTEPPER>
    PERTURBED_PROBLEM;

  std::shared_ptr<AXISYM_PROBLEM> createBaseProblem();
  std::shared_ptr<PERTURBED_PROBLEM> createLinearProblem(
    std::shared_ptr<AXISYM_PROBLEM> base_problem_pt,
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
