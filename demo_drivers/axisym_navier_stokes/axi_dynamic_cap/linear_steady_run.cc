// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Driver code for an axisymmetric free-surface hydrostatics problem.
// The system consists of a layer of fluid
// in a domain of height 1 and radius 0.5.
// The program solves for the interface position as the contact angle
// at the wall, alpha, decreases from pi/2. The resulting shapes should all be
// spherical shells and the pressure jump across the interface should be
// 2 cos(alpha)/0.5 = 4 cos(alpha)/Ca.

#include "algorithm"

// OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"

// The mesh
#include "meshes/triangle_mesh.h"


#include "linearised_axisym_navier_stokes.h"

#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"
#include "linearised_axisymmetric_fluid_interface_elements.h"
#include "decomposed_linear_elasticity_elements.h"

#include "linearised_elastic_axisym_fluid_interface_element.h"
#include "overlaying_linearised_elastic_axisym_fluid_interface_element.h"
#include "overlaying_my_linear_element.h"
#include "overlaying_Tlinear_axisym_ns_pvd_elements.h"

//#include "axisym_linear_stability_cap_problem.h"
#include "axisym_dynamic_cap_problem.h"
#include "perturbed_linear_stability_cap_problem.h"

#include "parameters.h"

using namespace std;
using namespace oomph;


//===start_of_main=======================================================
/// Main driver: Build problem and initiate parameter study
//======================================================================
int main(int argc, char** argv)
{
  // Check number of arguments
  int number_of_arguments = argc - 1;
  if (number_of_arguments == 0 || number_of_arguments > 1)
  {
    cout << "Wrong number of arguments." << std::endl;
    return 1;
  }

#ifdef OOMPH_HAS_MPI
  // Setup mpi but don't make a copy of mpi_comm_world because
  // mumps wants to work with the real thing.
  bool make_copy_of_mpi_comm_world = false;
  MPI_Helpers::init(argc, argv, make_copy_of_mpi_comm_world);
#endif

  // Problem parameters
  Parameters parameters;
  read_parameters_from_file(argv[1], parameters);

  // Construct the base problem
  bool has_restart = false;
  if (parameters.restart_filename != "")
  {
    cout << "restarting" << endl;
    has_restart = true;
  }
  typedef HijackedProjectableAxisymmetricTTaylorHoodPVDElement BASE_ELEMENT;
  typedef BDF<2> TIMESTEPPER;
  AxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER> base_problem(
    Global_Physical_Parameters::Equilibrium_contact_angle, has_restart);

  // Load in restart file
  if (parameters.restart_filename != "")
  {
    try
    {
      ifstream restart_filestream;
      restart_filestream.open(parameters.restart_filename);
      bool is_unsteady_restart = false;
      base_problem.read(restart_filestream, is_unsteady_restart);
      restart_filestream.close();
    }
    catch (exception& e)
    {
      cout << "Restart filename can't be set, or opened, or read." << endl;
      cout << "File: " << parameters.restart_filename << endl;
      return 1;
    }
  }

  base_problem.set_contact_angle(
    Global_Physical_Parameters::Equilibrium_contact_angle);
  base_problem.set_bond_number(Global_Physical_Parameters::Bo);
  base_problem.set_capillary_number(Global_Physical_Parameters::Ca);
  base_problem.set_reynolds_number(Global_Physical_Parameters::Re);

  // Set maximum number of mesh adaptations per solve
  base_problem.set_max_adapt(parameters.max_adapt);

  // Set output directory
  base_problem.set_directory(parameters.dir_name);

  // Setup trace file
  base_problem.open_trace_files(true);

  ofstream parameters_filestream(
    (parameters.dir_name + "/parameters.dat").c_str());
  parameters.doc(parameters_filestream);
  parameters_filestream.close();

  // Solve for the steady state adapting if needed by the Z2 error estimator
  base_problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);

  // Output result
  base_problem.create_restart_file();
  base_problem.doc_solution();

  // Close the trace files
  base_problem.close_trace_files();

  // Create the linear stability problem
  typedef OverlayingMyLinearElement PERTURBED_ELEMENT;
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

  perturbed_problem.assign_initial_values_impulsive();

  // Document the solution before the solve for testing
  perturbed_problem.doc_solution();

  // Solve steady problem
  perturbed_problem.make_steady();
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();
  DoubleVector residuals;
  DenseMatrix<double> actual_jacobian;
  perturbed_problem.get_fd_jacobian(residuals, actual_jacobian);

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

} // end_of_main
