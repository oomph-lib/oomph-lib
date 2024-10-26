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

// Solve for the steady state of a axisymmetric fluid surface in a capillary,
// under the force of gravity and a moving wall.
// Reads from a parameter file, see the default parameter file, and outputs
// bulk and surface data and restart files.

// OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Local include files
#include "singular_axisym_navier_stokes_elements.h"
#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"
#include "singular_axisym_dynamic_cap_problem.h"
#include "parameters.h"

using namespace std;
using namespace oomph;

int main(int argc, char** argv)
{
#ifdef OOMPH_HAS_MPI
  // Setup mpi but don't make a copy of mpi_comm_world because
  // mumps wants to work with the real thing.
  bool make_copy_of_mpi_comm_world = false;
  MPI_Helpers::init(argc, argv, make_copy_of_mpi_comm_world);
#endif

  // Check number of arguments
  int number_of_arguments = argc - 1;
  if (number_of_arguments == 0 || number_of_arguments > 1)
  {
    std::cout << "Wrong number of arguments." << std::endl;
    return 1;
  }

  // Problem parameters
  Parameters parameters;
  read_parameters_from_file(argv[1], parameters);

  // Check if we are restarting
  bool has_restart = false;
  if (parameters.restart_filename != "")
  {
    std::cout << "restarting" << std::endl;
    has_restart = true;
  }

  // Construct the problem
  SingularAxisymDynamicCapProblem<
    SingularAxisymNavierStokesElement<
      HijackedProjectableAxisymmetricTTaylorHoodPVDElement>,
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
      std::cout << "Restart filename can't be set, or opened, or read."
                << std::endl;
      std::cout << "File: " << parameters.restart_filename << std::endl;
      return 1;
    }
  }

  // Set the problem parameters
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

  // Save a copy of the parameters
  ofstream parameters_filestream(
    (parameters.dir_name + "/parameters.dat").c_str());
  parameters.doc(parameters_filestream);
  parameters_filestream.close();

  // Document initial condition
  problem.create_restart_file();
  problem.doc_solution();
  problem.max_newton_iterations() = 40;

  // If the final time is zero (or less) then we are doing a steady solve,
  if (parameters.ft <= 0)
  {
    // Solve for the steady state adapting if needed by the Z2 error estimator
    problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);

    // Create the restart file - needed before the doc solution
    problem.create_restart_file();

    // Document the solution
    problem.doc_solution();
  }
  // ...otherwise, we are doing an unsteady run
  else
  {
    // If the contact angle is acute, then timestep
    if (Global_Physical_Parameters::Equilibrium_contact_angle <=
        90.0 * MathematicalConstants::Pi / 180.0)
    {
      // Timestep until the desired final time
      problem.timestep(parameters.dt, parameters.ft);
    }
    // otherwise, throw a warning as we haven't implemented this yet
    else
    {
      throw(OomphLibWarning(
        "Timestepping is not implemented for obtuse contact angles yet.",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION));
    }
  }

  // Close the trace files
  problem.close_trace_files();

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

  // Return 0 to tell everyone that the program finished successfully
  return 0;
}
