#include <iostream>
#include <algorithm>

// OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"
#include "linearised_axisym_navier_stokes.h"

// Parameter
#include "parameters.h"

// Elements
#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"
#include "overlaying_my_linear_element.h"

// Problem
#include "axisym_dynamic_cap_problem.h"
#include "perturbed_linear_stability_cap_problem.h"


// Namespaces
using namespace std;
using namespace oomph;

// Define types
typedef HijackedProjectableAxisymmetricTTaylorHoodPVDElement BASE_ELEMENT;
typedef OverlayingMyLinearElement<BASE_ELEMENT> PERTURBED_ELEMENT;
typedef BDF<2> TIMESTEPPER;
typedef AxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER> BASE_PROBLEM;
typedef PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                           PERTURBED_ELEMENT,
                                           TIMESTEPPER>
  PERTURBED_PROBLEM;

// Parse input arguments
bool parse_arguments(int argc,
                     char** argv,
                     string& continuation_parameter_string,
                     double& parameter_step,
                     string& neutral_stability_parameter_string,
                     double& step_stability,
                     string& filename)
{
  // Check number of arguments
  int number_of_arguments = argc - 1;
  if (number_of_arguments != 5)
  {
    cout << "Wrong number of arguments." << std::endl;
    return 1;
  }
  continuation_parameter_string = argv[1];
  parameter_step = stod(argv[2]);
  neutral_stability_parameter_string = argv[3];
  step_stability = stod(argv[4]);
  filename = argv[5];
  return 0;
}

// Create the base problem
BASE_PROBLEM* create_base_problem(Parameters& parameters)
{
  bool has_restart = false;
  if (parameters.restart_filename != "")
  {
    cout << "restarting" << endl;
    has_restart = true;
  }

  BASE_PROBLEM* base_problem_pt = new BASE_PROBLEM(
    Global_Physical_Parameters::Equilibrium_contact_angle, has_restart);

  // Load in restart file
  if (parameters.restart_filename != "")
  {
    try
    {
      ifstream restart_filestream;
      restart_filestream.open(parameters.restart_filename);
      bool is_unsteady_restart = false;
      base_problem_pt->read(restart_filestream, is_unsteady_restart);
      restart_filestream.close();
    }
    catch (exception& e)
    {
      cout << "Restart filename can't be set, or opened, or read." << endl;
      cout << "File: " << parameters.restart_filename << endl;
      throw(e);
    }
  }

  base_problem_pt->set_contact_angle(
    Global_Physical_Parameters::Equilibrium_contact_angle);
  base_problem_pt->set_bond_number(Global_Physical_Parameters::Bo);
  base_problem_pt->set_capillary_number(Global_Physical_Parameters::Ca);
  base_problem_pt->set_reynolds_number(Global_Physical_Parameters::Re);

  // Set maximum number of mesh adaptations per solve
  base_problem_pt->set_max_adapt(parameters.max_adapt);

  // Set output directory
  base_problem_pt->set_directory(parameters.dir_name);

  // Setup trace file
  base_problem_pt->open_trace_files(true);

  ofstream parameters_filestream(
    (parameters.dir_name + "/parameters.dat").c_str());
  parameters.doc(parameters_filestream);
  parameters_filestream.close();

  return base_problem_pt;
}


PERTURBED_PROBLEM* create_perturbed_problem(
  BASE_PROBLEM* const& base_problem_pt, Parameters const& parameters)
{
  PERTURBED_PROBLEM* perturbed_problem_pt =
    new PERTURBED_PROBLEM(base_problem_pt->bulk_mesh_pt(),
                          base_problem_pt->free_surface_mesh_pt(),
                            base_problem_pt->slip_surface_mesh_pt(),
                          parameters.azimuthal_mode_number);

  perturbed_problem_pt->set_directory(parameters.dir_name);
  perturbed_problem_pt->set_contact_angle(
    Global_Physical_Parameters::Equilibrium_contact_angle);
  perturbed_problem_pt->set_bond_number(Global_Physical_Parameters::Bo);
  perturbed_problem_pt->set_capillary_number(Global_Physical_Parameters::Ca);
  perturbed_problem_pt->set_reynolds_number(Global_Physical_Parameters::Re);

  return perturbed_problem_pt;
}

void solve_the_steady_problem(BASE_PROBLEM* const& base_problem_pt,
                              Parameters const& parameters,
                              double* const& stability_pt,
                              double& step_stability)
{
  bool has_base_state = false;
  unsigned base_state_iterations = 0;
  const unsigned max_base_state_iterations = 8;
  while (!has_base_state && base_state_iterations < max_base_state_iterations)
  {
    DoubleVector dofs;
    base_problem_pt->get_dofs(dofs);
    try
    {
      base_problem_pt->set_contact_angle(
        Global_Physical_Parameters::Equilibrium_contact_angle);
      base_problem_pt->set_bond_number(Global_Physical_Parameters::Bo);
      base_problem_pt->set_capillary_number(Global_Physical_Parameters::Ca);
      base_problem_pt->set_reynolds_number(Global_Physical_Parameters::Re);

      // Solve steady problem
      int exit_flag = base_problem_pt->steady_newton_solve_adapt_if_needed(
        parameters.max_adapt);

      if (exit_flag >= 0)
      {
        has_base_state = true;
      }
    }
    catch (exception& e)
    {
      cout << "Caught exception" << endl;
      cout << "Resetting problem" << endl;
      base_problem_pt->set_dofs(dofs);
      *stability_pt -= step_stability;
      step_stability /= 3.0;
      *stability_pt += step_stability;
      cout << "Reducing step size.";
      cout << "Number of attempts: " << base_state_iterations;
      cout << ", Step size: " << step_stability;
      cout << ", Target wall velocity: " << *stability_pt << endl;
    }
    base_state_iterations++;
  }

  if (!has_base_state)
  {
    cout << "WARNING: Base state not found." << endl;
    throw("Base state not found");
  }
}

double get_the_eigenvalue(BASE_PROBLEM* const& base_problem_pt,
                          Parameters const& parameters)
{
  // Solve the eigenvalue problem
  Vector<std::complex<double>> eigenvalues(1, 0.0);

  // Create the linear stability problem
  PERTURBED_PROBLEM* perturbed_problem_pt =
    create_perturbed_problem(base_problem_pt, parameters);

  perturbed_problem_pt->set_doc_number(base_problem_pt->get_doc_number());
  const unsigned number_of_eigenvalues = 1;
  perturbed_problem_pt->steady_newton_solve();
  eigenvalues = perturbed_problem_pt->solve_n_most_unstable_eigensolutions(
    number_of_eigenvalues);

  delete perturbed_problem_pt;

  return real(eigenvalues[0]);
}

void newton_iterate_to_neutral_stability_point(
  BASE_PROBLEM* const& base_problem_pt,
  Parameters const& parameters,
  double*& stability_pt,
  double step_stability)
{
  unsigned n_iterations = 0;
  const unsigned max_n_iterations = 10;
  // const double step_tolerance = 1e-3;
  double old_param;
  double old_residual;
  double residual = get_the_eigenvalue(base_problem_pt, parameters);
  const double tolerance = 1e-5;
  bool has_converged = false;

  if (abs(residual) < tolerance)
  {
    has_converged = true;
  }

  while (!has_converged && n_iterations < max_n_iterations)
  {
    cout << "-------------" << endl;
    cout << "Start of loop" << endl;
    cout << "-------------" << endl;
    cout << endl;
    cout << n_iterations << ", " << *stability_pt << "," << residual << ","
         << *stability_pt + step_stability << endl;

    old_param = *stability_pt;
    *stability_pt += step_stability;

    // Solve steady problem
    solve_the_steady_problem(
      base_problem_pt, parameters, stability_pt, step_stability);


    //====================================================================

    old_residual = residual;
    residual = get_the_eigenvalue(base_problem_pt, parameters);

    //====================================================================

    if (abs(residual) < tolerance)
    {
      has_converged = true;
    }
    else
    {
      double deriv = (residual - old_residual) / (*stability_pt - old_param);
      const double relaxation = 1.0;
      step_stability = -relaxation * residual / (deriv);
    }

    n_iterations++;
  }

  base_problem_pt->create_restart_file();
  base_problem_pt->doc_solution();

  if (!has_converged)
  {
    cout << "WARNING: Critical Bond number not converged." << endl;
    throw("Critical Bond number not converged");
  }
}


void step_continuation_parameter(BASE_PROBLEM* const& base_problem_pt,
                                 Parameters const& parameters,
                                 double*& param_pt,
                                 double& step_param,
                                 double*& stability_pt,
                                 double step_stability)
{
  *param_pt += step_param;
  bool has_completed_a_step = false;
  unsigned n = 0;
  while (!has_completed_a_step && n < 10)
  {
    try
    {
      solve_the_steady_problem(
        base_problem_pt, parameters, param_pt, step_param);

      newton_iterate_to_neutral_stability_point(
        base_problem_pt, parameters, stability_pt, step_stability);
      has_completed_a_step = true;
    }
    catch (exception& e)
    {
      *param_pt -= step_param;
      step_param /= 2.0;
      *param_pt += step_param;
    }
    n++;
  }
  if (!has_completed_a_step)
  {
    throw("Did not complete continuation step.");
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

  string continuation_parameter_string;
  double step_param;
  string neutral_stability_parameter_string;
  double step_stability;
  string filename;
  if (parse_arguments(argc,
                      argv,
                      continuation_parameter_string,
                      step_param,
                      neutral_stability_parameter_string,
                      step_stability,
                      filename))
  {
    return 1;
  }

  //====================================================================

  double* param_pt = 0;
  try
  {
    if (continuation_parameter_string == "--Bo")
    {
      param_pt = &Global_Physical_Parameters::Bo;
    }
    else if (continuation_parameter_string == "--angle")
    {
      param_pt = &Global_Physical_Parameters::Equilibrium_contact_angle;
    }
    else if (continuation_parameter_string == "--slip_length")
    {
      param_pt = &Slip_Parameters::slip_length;
    }
    else if (continuation_parameter_string == "--wall_velocity")
    {
      param_pt = &Slip_Parameters::wall_velocity;
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

  double* stability_pt = 0;
  try
  {
    if (neutral_stability_parameter_string == "--Bo")
    {
      stability_pt = &Global_Physical_Parameters::Bo;
    }
    else if (neutral_stability_parameter_string == "--angle")
    {
      stability_pt = &Global_Physical_Parameters::Equilibrium_contact_angle;
    }
    else if (neutral_stability_parameter_string == "--slip_length")
    {
      stability_pt = &Slip_Parameters::slip_length;
    }
    else if (neutral_stability_parameter_string == "--wall_velocity")
    {
      stability_pt = &Slip_Parameters::wall_velocity;
    }
    else
    {
      throw std::invalid_argument("Third argument incorrect.");
    }
  }
  catch (exception& e)
  {
    cout << "Stability parameter can't be set." << endl;
    cout << "Argument in: " << neutral_stability_parameter_string << endl;
    return 1;
  }

  //====================================================================

  // Problem parameters
  Parameters parameters;
  read_parameters_from_file(filename, parameters);

  // Construct the base problem
  BASE_PROBLEM* base_problem_pt = create_base_problem(parameters);
  int exit_flag =
    base_problem_pt->steady_newton_solve_adapt_if_needed(parameters.max_adapt);
  if (exit_flag < 0)
  {
    oomph_info << "Couldn't solve for initial state." << endl;
    throw("Couldn't solve for initial state.");
  }

  //====================================================================

  newton_iterate_to_neutral_stability_point(
    base_problem_pt, parameters, stability_pt, step_stability);

  //====================================================================

  const unsigned n_steps = 30;
  for (unsigned n = 0; n < n_steps; n++)
  {
    step_continuation_parameter(base_problem_pt,
                                parameters,
                                param_pt,
                                step_param,
                                stability_pt,
                                step_stability);
  }

  //====================================================================

  // Close the trace files
  base_problem_pt->close_trace_files();

  TerminateHelper::clean_up_memory();

  delete base_problem_pt;

  // Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

} // end_of_main
