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
#include "singular_axisym_navier_stokes_elements.h"

#include "linearised_elastic_axisym_fluid_interface_element.h"
#include "overlaying_linearised_elastic_axisym_fluid_interface_element.h"
#include "overlaying_my_linear_element.h"

//#include "axisym_linear_stability_cap_problem.h"
#include "singular_axisym_dynamic_cap_problem.h"
#include "perturbed_linear_stability_cap_problem.h"

#include "parameters.h"

using namespace std;
using namespace oomph;

void parse_arguments(const int& argc,
                     char** const& argv,
                     string& parameters_filename,
                     double& starting_step,
                     double*& continuation_param_pt)
{
  CommandLineArgs::setup(argc, argv);

  CommandLineArgs::specify_command_line_flag(
    "--parameters", &parameters_filename, "Required: parameter filename");

  double bo_initial_step = 0.0;
  CommandLineArgs::specify_command_line_flag(
    "--Bo",
    &bo_initial_step,
    "Optional: Continue in Bond number with initial step");

  double ca_initial_step = 0.0;
  CommandLineArgs::specify_command_line_flag(
    "--wall_velocity",
    &ca_initial_step,
    "Optional: Continue in Capillary number with initial step");

  double angle_initial_step = 0.0;
  CommandLineArgs::specify_command_line_flag(
    "--angle",
    &angle_initial_step,
    "Optional: Continue in angle with initial step");

  double slip_length_initial_step = 0.0;
  CommandLineArgs::specify_command_line_flag(
    "--slip_length",
    &slip_length_initial_step,
    "Optional: Continue in slip length with initial step");

  CommandLineArgs::specify_command_line_flag("--arc",
                                             "Optional: Use arc continuation");

  CommandLineArgs::specify_command_line_flag(
    "--height_control", "Optional: Use height control continuation");


  // Parse and assign command line arguments
  bool has_unrecognised_arg = false;
  CommandLineArgs::parse_and_assign(argc, argv, has_unrecognised_arg);
  if (has_unrecognised_arg)
  {
    throw std::invalid_argument("Unrecognised args.");
  }
  if (parameters_filename == "")
  {
    throw std::invalid_argument("Parameter file not set.");
  }

  // Set continuation variables
  try
  {
    if (bo_initial_step != 0)
    {
      starting_step = bo_initial_step;
      continuation_param_pt = &Global_Physical_Parameters::Bo;
    }
    else if (angle_initial_step != 0)
    {
      starting_step = angle_initial_step;
      continuation_param_pt =
        &Global_Physical_Parameters::Equilibrium_contact_angle;
    }
    else if (ca_initial_step != 0)
    {
      starting_step = ca_initial_step;
      continuation_param_pt = &Slip_Parameters::wall_velocity;
    }
    else if (slip_length_initial_step != 0)
    {
      starting_step = slip_length_initial_step;
      continuation_param_pt = &Slip_Parameters::slip_length;
    }
    else
    {
      throw std::invalid_argument("No continuation parameter specified.");
    }
  }
  catch (exception& e)
  {
    throw std::invalid_argument("Continuation parameter can't be set.");
  }
}

//===start_of_main=======================================================
/// Main driver: Build problem and initiate parameter study
//======================================================================
int main(int argc, char** argv)
{
#ifdef OOMPH_HAS_MPI
  // Setup mpi but don't make a copy of mpi_comm_world because
  // mumps wants to work with the real thing.
  bool make_copy_of_mpi_comm_world = false;
  MPI_Helpers::init(argc, argv, make_copy_of_mpi_comm_world);
#endif

  string parameters_filename = "";
  double starting_step = 0.0;
  double* continuation_param_pt = 0;
  parse_arguments(
    argc, argv, parameters_filename, starting_step, continuation_param_pt);

  // Problem parameters
  Parameters parameters;
  read_parameters_from_file(parameters_filename, parameters);

  // Construct the base problem
  bool has_restart = false;
  if (parameters.restart_filename != "")
  {
    cout << "restarting" << endl;
    has_restart = true;
  }
  typedef SingularAxisymNavierStokesElement<
    HijackedProjectableAxisymmetricTTaylorHoodPVDElement>
    BASE_ELEMENT;
  typedef BDF<2> TIMESTEPPER;
  SingularAxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER> base_problem(
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
    catch (OomphLibError& e)
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

  // Output the initial condition
  base_problem.create_restart_file();
  base_problem.doc_solution();

  //====================================================================

  // Solve steady problem
  base_problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);

  // Output result
  base_problem.create_restart_file();
  base_problem.doc_solution();

  // Solve the eigenvalue problem
  Vector<std::complex<double>> eigenvalues(1, 0.0);

  //====================================================================

  // Create the linear stability problem
  typedef OverlayingMyLinearElement<BASE_ELEMENT> PERTURBED_ELEMENT;
  PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                     PERTURBED_ELEMENT,
                                     TIMESTEPPER>
    initial_perturbed_problem(base_problem.bulk_mesh_pt(),
                              base_problem.free_surface_mesh_pt(),
                              base_problem.slip_surface_mesh_pt(),
                              parameters.azimuthal_mode_number);

  initial_perturbed_problem.set_directory(parameters.dir_name);
  initial_perturbed_problem.set_contact_angle(
    Global_Physical_Parameters::Equilibrium_contact_angle);
  initial_perturbed_problem.set_bond_number(Global_Physical_Parameters::Bo);
  initial_perturbed_problem.set_capillary_number(
    Global_Physical_Parameters::Ca);
  initial_perturbed_problem.set_reynolds_number(Global_Physical_Parameters::Re);

  // Document the solution before the solve for testing
  initial_perturbed_problem.doc_solution();
  initial_perturbed_problem.steady_newton_solve();

  eigenvalues =
    initial_perturbed_problem.solve_n_most_unstable_eigensolutions(1);

  double residual = real(eigenvalues[0]);

  const double tolerance = 1e-5;
  bool has_converged = false;
  if (abs(residual) < tolerance)
  {
    has_converged = true;
  }

  //====================================================================

  unsigned n_iterations = 0;
  const unsigned max_n_iterations = floor(abs(parameters.ft / parameters.dt));
  // const double step_tolerance = 1e-3;
  double step_param = starting_step;
  double old_param;
  double current_param = *continuation_param_pt;
  double old_residual;
  while (!has_converged && n_iterations < max_n_iterations)
  {
    TerminateHelper::setup();
    cout << "-------------" << endl;
    cout << "Start of loop" << endl;
    cout << "-------------" << endl;
    cout << endl;
    cout << n_iterations << ", " << current_param << "," << residual << ","
         << current_param + step_param << endl;
    old_param = current_param;
    current_param += step_param;
    *continuation_param_pt = current_param;


    // Solve steady problem
    bool has_base_state = false;
    unsigned base_state_iterations = 0;
    const unsigned max_base_state_iterations = 8;
    while (!has_base_state && base_state_iterations < max_base_state_iterations)
    {
      DoubleVector dofs;
      base_problem.get_dofs(dofs);
      try
      {
        base_problem.set_contact_angle(
          Global_Physical_Parameters::Equilibrium_contact_angle);
        base_problem.set_bond_number(Global_Physical_Parameters::Bo);
        base_problem.set_capillary_number(Global_Physical_Parameters::Ca);
        base_problem.set_reynolds_number(Global_Physical_Parameters::Re);

        // Solve steady problem
        int exit_flag = base_problem.steady_newton_solve_adapt_if_needed(
          parameters.max_adapt);

        // Output result
        base_problem.create_restart_file();
        base_problem.doc_solution();

        if (exit_flag >= 0)
        {
          has_base_state = true;
        }
      }
      catch (OomphLibError& e)
      {
        cout << "Caught exception" << endl;
        cout << "Resetting problem" << endl;
        base_problem.set_dofs(dofs);
        current_param -= step_param;
        step_param /= 3.0;
        current_param += step_param;
        *continuation_param_pt = current_param;
        cout << "Reducing step size.";
        cout << "Number of attempts: " << base_state_iterations;
        cout << ", Step size: " << step_param;
        cout << ", Target wall velocity: " << current_param << endl;
      }
      base_state_iterations++;
    }

    if (!has_base_state)
    {
      cout << "WARNING: Base state not found." << endl;
      break;
    }

    //====================================================================

    // Create the linear stability problem
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

    // Document the solution before the solve for testing
    perturbed_problem.doc_solution();

    perturbed_problem.steady_newton_solve();

    eigenvalues = perturbed_problem.solve_n_most_unstable_eigensolutions(1);

    old_residual = residual;
    residual = real(eigenvalues[0]);

    if (abs(residual) < tolerance)
    {
      has_converged = true;
      perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
    }
    else
    {
      double deriv = (residual - old_residual) / (current_param - old_param);
      const double relaxation = 1.0;
      step_param = -relaxation * residual / (deriv);
      if (n_iterations == (max_n_iterations - 1))
      {
        perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
      }
    }


    n_iterations++;
  }

  if (!has_converged)
  {
    cout << "WARNING: Critical Bond number not converged." << endl;
  }

  // Close the trace files
  base_problem.close_trace_files();

  TerminateHelper::clean_up_memory();

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

} // end_of_main
