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

//#include "axisym_linear_stability_cap_problem.h"
#include "axisym_dynamic_cap_problem.h"
#include "perturbed_linear_stability_cap_problem.h"

#include "parameters.h"

using namespace std;
using namespace oomph;

typedef HijackedProjectableAxisymmetricTTaylorHoodPVDElement BASE_ELEMENT;
typedef OverlayingMyLinearElement PERTURBED_ELEMENT;
typedef BDF<2> TIMESTEPPER;

void solve_for_base_state(
  AxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER>& base_problem,
  Parameters& parameters)
{
  base_problem.set_contact_angle(
    Global_Physical_Parameters::Equilibrium_contact_angle);
  base_problem.set_bond_number(Global_Physical_Parameters::Bo);
  base_problem.set_capillary_number(Global_Physical_Parameters::Ca);
  base_problem.set_reynolds_number(Global_Physical_Parameters::Re);

  base_problem.create_restart_file();
  base_problem.doc_solution();

  DoubleVector dofs;
  base_problem.get_dofs(dofs);
  try
  {
    //====================================================================
    // Solve steady problem
    cout << "Solving base problem." << endl;
    base_problem.steady_newton_solve_adapt_if_needed(
       parameters.max_adapt);

    // Output result
    base_problem.create_restart_file();
    base_problem.doc_solution();
  }
  catch (exception& e)
  {
    base_problem.set_dofs(dofs);
  }
}

bool solve_steady_problem(
  AxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER>& base_problem,
  const Parameters& parameters,
  double& current_param,
  double& step_param,
  double*& param_pt)
{
  cout << "solve_steady_problem" << endl;
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
    catch (exception& e)
    {
      cout << "Caught exception" << endl;
      cout << "Resetting problem" << endl;
      base_problem.set_dofs(dofs);
      current_param -= step_param;
      step_param /= 3.0;
      current_param += step_param;
      *param_pt = current_param;
      cout << "Halving step size.";
      cout << "Number of attempts: " << base_state_iterations;
      cout << ", Step size: " << step_param;
      cout << ", Target wall velocity: " << current_param << endl;
    }
    base_state_iterations++;
  }
  return has_base_state;
}

double get_most_unstable_eigenvalue(
  AxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER>& base_problem,
  const Parameters& parameters)
{
  cout << "get_most_unstable_eigenvalue" << endl;

  PerturbedLinearStabilityCapProblem<BASE_ELEMENT, PERTURBED_ELEMENT, TIMESTEPPER>
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

  Vector<std::complex<double>> eigenvalues(8, 0.0);
  eigenvalues =
    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(8);

  return real(eigenvalues[0]);
}

bool newton_iterate_to_critical_parameter(
  AxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER>& base_problem,
  const Parameters& parameters,
  double*& param_pt)
{
  //====================================================================
  // Solve the eigenvalue problem
  double residual = get_most_unstable_eigenvalue(base_problem, parameters);

  const double tolerance = 1e-3;
  bool has_converged = false;
  if (abs(residual) < tolerance)
  {
    has_converged = true;
  }


  unsigned n_iterations = 0;
  const unsigned max_n_iterations = 10;
  // const double step_tolerance = 1e-3;
  double step_param = 1e-3;
  double old_param;
  double current_param = *param_pt;
  double old_residual;
  while (!has_converged && n_iterations < max_n_iterations)
  {
    cout << "-------------" << endl;
    cout << "Start of loop" << endl;
    cout << "-------------" << endl;
    cout << endl;
    cout << n_iterations << ", " << current_param << "," << residual << ","
         << current_param + step_param << endl;
    old_param = current_param;
    current_param += step_param;
    *param_pt = current_param;
    cout << current_param << endl;
    cout << *param_pt << endl;
    cout << Global_Physical_Parameters::Bo << endl;

    // Solve steady problem
    bool has_base_state = solve_steady_problem(
      base_problem, parameters, current_param, step_param, param_pt);

    if (!has_base_state)
    {
      cout << "WARNING: Base state not found." << endl;
      break;
    }

    old_residual = residual;
    residual = get_most_unstable_eigenvalue(base_problem, parameters);

    if (abs(residual) < tolerance)
    {
      has_converged = true;
    }
    else
    {
      double deriv = (residual - old_residual) / (current_param - old_param);
      const double relaxation = 1.0;
      step_param = -relaxation * residual / (deriv);
    }

    n_iterations++;
  }

  return has_converged;
}


//===start_of_main=======================================================
/// Main driver: Build problem and initiate parameter study
//======================================================================
int main(int argc, char** argv)
{
  cout << "neutral_stability_continuation_run" << endl;
  oomph_info.stream_pt() = new ofstream("output.dat");

  // Check number of arguments
  int number_of_arguments = argc - 1;
  if (number_of_arguments != 3)
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
  read_parameters_from_file(argv[3], parameters);

  // Construct the base problem
  bool has_restart = false;
  if (parameters.restart_filename != "")
  {
    cout << "restarting" << endl;
    has_restart = true;
  }

  //====================================================================
  cout << "Creating base problem" << endl;
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

  // Output the initial condition
  base_problem.create_restart_file();
  base_problem.doc_solution();

  //====================================================================

  string critical_parameter_string = "";
  double* critical_param_pt = 0;
  try
  {
    critical_parameter_string = argv[1];
    if (critical_parameter_string == "--Bo")
    {
      critical_param_pt = &Global_Physical_Parameters::Bo;
    }
    else if (critical_parameter_string == "--angle")
    {
      critical_param_pt =
        &Global_Physical_Parameters::Equilibrium_contact_angle;
    }
    else if (critical_parameter_string == "--wall_velocity")
    {
      critical_param_pt = &Slip_Parameters::wall_velocity;
    }
    else
    {
      throw std::invalid_argument("First argument incorrect.");
    }
  }
  catch (exception& e)
  {
    cout << "Continuation parameter can't be set." << endl;
    cout << "Argument in: " << critical_parameter_string << endl;
    return 1;
  }

  string continuation_parameter_string = "";
  double* continuation_param_pt = 0;
  try
  {
    continuation_parameter_string = argv[2];
    if (continuation_parameter_string == "--Bo")
    {
      continuation_param_pt = &Global_Physical_Parameters::Bo;
    }
    else if (continuation_parameter_string == "--angle")
    {
      continuation_param_pt =
        &Global_Physical_Parameters::Equilibrium_contact_angle;
    }
    else if (continuation_parameter_string == "--wall_velocity")
    {
      continuation_param_pt = &Slip_Parameters::wall_velocity;
    }
    else
    {
      throw std::invalid_argument("First argument incorrect.");
    }
  }
  catch (exception& e)
  {
    cout << "Continuation parameter can't be set." << endl;
    cout << "Argument in: " << continuation_parameter_string << endl;
    return 1;
  }

  const double target_value = 30 * oomph::MathematicalConstants::Pi / 180.0;
  double step = -1.0 * oomph::MathematicalConstants::Pi / 180.0;
  *continuation_param_pt -= step;
  unsigned outer_iterations = 0;
  while (*continuation_param_pt >= target_value && outer_iterations < 2)
  {
    outer_iterations++;
    cout << "value: " << *continuation_param_pt << " target: " << target_value
         << endl;
    step = -1.0 * oomph::MathematicalConstants::Pi / 180.0;
    bool successful_step = false;
    unsigned iterations = 0;
    while (!successful_step && iterations < 10)
    {
      iterations++;
      cout << "step: " << step << endl;
      *continuation_param_pt += step;
      try
      {
        solve_for_base_state(base_problem, parameters);

        bool has_converged = newton_iterate_to_critical_parameter(
          base_problem, parameters, critical_param_pt);

        if (!has_converged)
        {
          *continuation_param_pt -= step;
          step /= 2.0;
        }
        else
        {
          successful_step = true;
        }
      }
      catch (exception& e)
      {
        cout << "halving step" << endl;
        *continuation_param_pt -= step;
        step /= 2.0;
      }
    }
  }

  // Close the trace files
  base_problem.close_trace_files();

  TerminateHelper::clean_up_memory();

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

} // end_of_main
