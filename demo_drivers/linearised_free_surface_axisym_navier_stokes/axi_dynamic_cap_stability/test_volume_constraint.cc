#include <iostream>

#include "algorithm"

// OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"


#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"
#include "overlaying_linearised_elastic_axisym_fluid_interface_element.h"
#include "overlaying_my_linear_element.h"
#include "axisym_dynamic_cap_problem.h"
#include "perturbed_linear_stability_cap_problem.h"
#include "parameters.h"

using namespace oomph;

int main()
{
  Parameters parameters;
  read_parameters_from_file("default_parameters.dat", parameters);

  typedef HijackedProjectableAxisymmetricTTaylorHoodPVDElement BASE_ELEMENT;
  typedef BDF<2> TIMESTEPPER;
  AxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER> base_problem(
    Global_Physical_Parameters::Equilibrium_contact_angle, false);

  base_problem.set_contact_angle(
    Global_Physical_Parameters::Equilibrium_contact_angle);
  base_problem.set_bond_number(Global_Physical_Parameters::Bo);
  base_problem.set_capillary_number(Global_Physical_Parameters::Ca);
  base_problem.set_reynolds_number(Global_Physical_Parameters::Re);
  base_problem.set_max_adapt(parameters.max_adapt);
  base_problem.set_directory(parameters.dir_name);

  base_problem.steady_newton_solve();
  base_problem.doc_solution();

  typedef OverlayingMyLinearElement<BASE_ELEMENT> PERTURBED_ELEMENT;
  PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                     PERTURBED_ELEMENT,
                                     TIMESTEPPER>
    perturbed_problem(base_problem.bulk_mesh_pt(),
                      base_problem.free_surface_mesh_pt(),
                      base_problem.slip_surface_mesh_pt(),
                      parameters.azimuthal_mode_number);

  perturbed_problem.set_directory(parameters.dir_name);

  perturbed_problem.set_contact_angle(
    Global_Physical_Parameters::Equilibrium_contact_angle);
  perturbed_problem.set_bond_number(Global_Physical_Parameters::Bo);
  perturbed_problem.set_capillary_number(Global_Physical_Parameters::Ca);
  perturbed_problem.set_reynolds_number(Global_Physical_Parameters::Re);

  perturbed_problem.set_initial_condition();
  perturbed_problem.isolate_volume();
  perturbed_problem.doc_full_eigenmatrices();
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  if (abs(perturbed_problem.get_volume() + 0.00875)/0.00875 < 1e-2)
  {
    std::cout << "." << std::endl;
  }
  else
  {
    std::cout << "F" << std::endl;
  }

  Slip_Parameters::wall_velocity = 0.02;
  base_problem.steady_newton_solve();
  base_problem.doc_solution();

  PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                     PERTURBED_ELEMENT,
                                     TIMESTEPPER>
    perturbed_problem2(base_problem.bulk_mesh_pt(),
                       base_problem.free_surface_mesh_pt(),
                      base_problem.slip_surface_mesh_pt(),
                       parameters.azimuthal_mode_number);

  perturbed_problem2.set_directory(parameters.dir_name);

  perturbed_problem2.set_contact_angle(
    Global_Physical_Parameters::Equilibrium_contact_angle);
  perturbed_problem2.set_bond_number(Global_Physical_Parameters::Bo);
  perturbed_problem2.set_capillary_number(Global_Physical_Parameters::Ca);
  perturbed_problem2.set_reynolds_number(Global_Physical_Parameters::Re);

  perturbed_problem2.set_initial_condition();
  perturbed_problem2.isolate_volume();
  perturbed_problem2.steady_newton_solve();
  perturbed_problem2.doc_solution();

  if (abs(perturbed_problem2.get_volume() + 0.00875)/0.00875 < 1e-2)
  {
    std::cout << "." << std::endl;
  }
  else
  {
    std::cout << "F" << std::endl;
  }

  return 0;
}
