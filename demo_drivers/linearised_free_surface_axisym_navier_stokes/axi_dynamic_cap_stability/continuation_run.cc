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

// OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Local include files
#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"
#include "axisym_dynamic_cap_problem.h"
#include "parameters.h"

using namespace std;
using namespace oomph;

int main(int argc, char** argv)
{
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

  // Construct the problem
  bool has_restart = false;
  if (parameters.restart_filename != "")
  {
    cout << "restarting" << endl;
    has_restart = true;
  }
  AxisymDynamicCapProblem<HijackedProjectableAxisymmetricTTaylorHoodPVDElement,
                        BDF<2>>
    problem(Global_Physical_Parameters::Equilibrium_contact_angle, has_restart);

  // Load in restart file
  if (parameters.restart_filename != "")
  {
    try
    {
      ifstream restart_filestream;
      restart_filestream.open(parameters.restart_filename);
      bool is_unsteady_restart = false;
      problem.read(restart_filestream, is_unsteady_restart);
      restart_filestream.close();
    }
    catch (exception& e)
    {
      cout << "Restart filename can't be set, or opened, or read." << endl;
      cout << "File: " << parameters.restart_filename << endl;
      return 1;
    }
  }

  problem.set_contact_angle(
    Global_Physical_Parameters::Equilibrium_contact_angle);
  problem.set_bond_number(Global_Physical_Parameters::Bo);
  problem.set_capillary_number(Global_Physical_Parameters::Ca);
  problem.set_reynolds_number(Global_Physical_Parameters::Re);

  // Set maximum number of mesh adaptations per solve
  problem.set_max_adapt(parameters.max_adapt);

  // Set output directory
  problem.set_directory(parameters.dir_name);

  // Setup trace file
  problem.open_trace_files(true);

  ofstream parameters_filestream(
    (parameters.dir_name + "/parameters.dat").c_str());
  parameters.doc(parameters_filestream);
  parameters_filestream.close();

  // Document initial condition
  problem.create_restart_file();
  problem.doc_solution();

  string continuation_parameter_string = "";
  double* continuation_param_pt = 0;
  try
  {
    continuation_parameter_string = argv[1];
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
    else if (continuation_parameter_string == "--slip_length")
    {
      continuation_param_pt = &Slip_Parameters::slip_length;
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

  double starting_step = 0.0;
  try
  {
    starting_step = atof(argv[2]);
  }
  catch (exception& e)
  {
    cout << "Target step size can't be set." << endl;
    cout << "Argument in: " << argv[2] << endl;
    return 1;
  }

  double step = starting_step;
  *continuation_param_pt -= step;
  const unsigned max_iterations = 20;
  unsigned iterations = 0;
  while (iterations < max_iterations)
  {
    // Output info
    cout << "Iteration: " << iterations << ", param: " << *continuation_param_pt
         << ", step: " << step << endl;

    // Update iteration number
    iterations++;

    // Continue parameter
    *continuation_param_pt += step;
    problem.set_contact_angle(
      Global_Physical_Parameters::Equilibrium_contact_angle);
    problem.set_bond_number(Global_Physical_Parameters::Bo);
    problem.set_capillary_number(Global_Physical_Parameters::Ca);
    problem.set_reynolds_number(Global_Physical_Parameters::Re);

    // Store current state of problem
    DoubleVector dofs;
    problem.get_dofs(dofs);
    try
    {
      // Solve for the steady state adapting if needed by the Z2 error
      // estimator
      problem.steady_newton_solve_adapt_if_needed(
        parameters.max_adapt);

      // Create the restart file - needed before the doc solution
      problem.create_restart_file();

      // Document the solution
      problem.doc_solution();
    }
    // If error in solving for the steady state
    catch (OomphLibException& e)
    {
      // Restore state of the problem
      problem.set_dofs(dofs);

      // Undo step
      *continuation_param_pt -= step;

      // Reduce step size
      step /= 2.0;
    }
  }

  // Close the trace files
  problem.close_trace_files();

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

  return 0;
}
