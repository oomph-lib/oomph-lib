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

// extern template class AxisymDynamicCapProblem<
//   HijackedProjectableAxisymmetricTTaylorHoodPVDElement,
//   BDF<2>>;

int main(int argc, char** argv)
{
  // Check number of arguments
  int number_of_arguments = argc - 1;
  if (number_of_arguments != 1 && number_of_arguments != 2)
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

  // Set output directory
  problem.set_directory(parameters.dir_name);

  // Solve for the steady state adapting if needed by the Z2 error estimator
  problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);
  problem.create_restart_file();
  problem.doc_solution();

  // problem.check_first_eigensolution_solves_the_problem();

  // If we have more than one argument (this is a flag for debugging the
  // jacobian),
  if (number_of_arguments == 2)
  {
    // then debug the jacobian and mass matrix by comparing to the full FD
    // version
    problem.debug_jacobian();
    problem.debug_mass_matrix();

    std::ofstream output_stream(parameters.dir_name + "/dofs.txt");
    problem.describe_dofs(output_stream);
    output_stream.close();
  }
  else
  {
    // otherwise, just document the Jacobian
    DoubleVector residuals;
    DenseDoubleMatrix jacobian;

    problem.get_jacobian(residuals, jacobian);

    std::ofstream output_stream(parameters.dir_name + "/jac.dat");
    output_stream.precision(16);
    jacobian.output(output_stream);
    output_stream.close();

    CRDoubleMatrix jacobianCR;
    CRDoubleMatrix mass_matrix;
    problem.get_eigenproblem_matrices(mass_matrix, jacobianCR);

    output_stream.open(parameters.dir_name + "/jacEig.dat");
    jacobianCR.sparse_indexed_output(output_stream, 16, true);
    output_stream.close();

    output_stream.open(parameters.dir_name + "/mass.dat");
    mass_matrix.sparse_indexed_output(output_stream, 16, true);
    output_stream.close();

    output_stream.open(parameters.dir_name + "/dofs.txt");
    problem.describe_dofs(output_stream);
    output_stream.close();
  }

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

  return 0;
}
