#include <iostream>
#include <algorithm>

// OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Parameter
#include "parameters.h"

// Elements
#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"

// Problem
#include "axisym_dynamic_cap_problem.h"


// Namespaces
using namespace std;
using namespace oomph;

// Define types
typedef HijackedProjectableAxisymmetricTTaylorHoodPVDElement BASE_ELEMENT;
typedef BDF<2> TIMESTEPPER;
typedef AxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER> BASE_PROBLEM;

// Parse input arguments
bool parse_arguments(int argc,
                     char** argv,
                     std::string& continuation_parameter_string,
                     double& parameter_step,
                     std::string& neutral_stability_parameter_string,
                     double& step_stability,
                     std::string& filename)
{
  // Check number of arguments
  int number_of_arguments = argc - 1;
  if (number_of_arguments != 5)
  {
    std::cout << "Wrong number of arguments." << std::endl;
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
    std::cout << "restarting" << std::endl;
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
      std::cout << "Restart filename can't be set, or opened, or read." << std::endl;
      std::cout << "File: " << parameters.restart_filename << std::endl;
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

  // base_problem_pt->set_analytic_dparameter(&Slip_Parameters::wall_velocity);

  // Setup trace file
  base_problem_pt->open_trace_files(true);

  ofstream parameters_filestream(
    (parameters.dir_name + "/parameters.dat").c_str());
  parameters.doc(parameters_filestream);
  parameters_filestream.close();

  return base_problem_pt;
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
      std::cout << "Caught exception" << std::endl;
      std::cout << "Resetting problem" << std::endl;
      base_problem_pt->set_dofs(dofs);
      *stability_pt -= step_stability;
      step_stability /= 3.0;
      *stability_pt += step_stability;
      std::cout << "Reducing step size.";
      std::cout << "Number of attempts: " << base_state_iterations;
      std::cout << ", Step size: " << step_stability;
      std::cout << ", Target wall velocity: " << *stability_pt << std::endl;
    }
    base_state_iterations++;
  }

  if (!has_base_state)
  {
    std::cout << "WARNING: Base state not found." << std::endl;
    throw("Base state not found");
  }
}

void arc_length_continue_to_the_limit_point(
  BASE_PROBLEM* const& base_problem_pt,
  Parameters const& parameters,
  double* const& param_pt,
  double*& stability_pt,
  double step_stability)
{
  unsigned n_iterations = 0;
  const unsigned max_n_iterations = 30;
  // const double step_tolerance = 1e-3;
  double old_param;
  bool has_converged = false;
  bool old_direction = false;
  bool direction = false;

  double ds = -0.01;

  while (!has_converged && n_iterations < max_n_iterations)
  {
    std::cout << "-------------" << std::endl;
    std::cout << "Start of loop" << std::endl;
    std::cout << "-------------" << std::endl;
    std::cout << std::endl;
    std::cout << n_iterations << ", " << *stability_pt << "," << direction << ","
         << ds << std::endl;

    old_param = *stability_pt;

    //====================================================================

    // Arc length step forward
    DoubleVector backup_dofs;
    base_problem_pt->get_dofs(backup_dofs);
    base_problem_pt->arc_length_step_solve(stability_pt, ds);

    //====================================================================

    old_direction = direction;
    direction = old_param > *stability_pt;

    // If we are going in the wrong direction on the first step,
    if ((!direction) && (n_iterations == 0))
    {
      // then reverse the ds.
      ds = -ds;
    }
    // If we have swapped direction after the first step,
    if ((old_direction != direction) && (n_iterations > 0))
    {
      // then we have converged
      base_problem_pt->set_dofs(backup_dofs);
      has_converged = true;
    }
    // If the step size has become too small
    // if (std::fabs(*stability_pt - old_param) < 1e-4)
    // {
    //   //, then assume we have converged too.
    //   oomph_info << "Step size is small, assuming convergence." << std::endl;
    //   has_converged = true;
    // }

    base_problem_pt->create_restart_file();
    base_problem_pt->doc_solution();

    n_iterations++;
  }

  base_problem_pt->create_restart_file();
  base_problem_pt->doc_solution();
  // Document
  ofstream my_trace("my_trace.dat", std::ios_base::app);
  my_trace << base_problem_pt->get_doc_number() << " " << *param_pt << " "
           << *stability_pt << " " << old_param << " "
           << (*stability_pt + old_param) / 2.0 << std::endl;
  my_trace.close();

  if (!has_converged)
  {
    std::cout << "WARNING: Critical Bond number not converged." << std::endl;
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

      arc_length_continue_to_the_limit_point(
        base_problem_pt, parameters, param_pt, stability_pt, step_stability);
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

  std::string continuation_parameter_string;
  double step_param;
  std::string neutral_stability_parameter_string;
  double step_stability;
  std::string filename;
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
    std::cout << "Continuation parameter can't be set." << std::endl;
    std::cout << "Argument in: " << continuation_parameter_string << std::endl;
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
    std::cout << "Stability parameter can't be set." << std::endl;
    std::cout << "Argument in: " << neutral_stability_parameter_string << std::endl;
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
    oomph_info << "Couldn't solve for initial state." << std::endl;
    throw("Couldn't solve for initial state.");
  }

  //====================================================================

  arc_length_continue_to_the_limit_point(
    base_problem_pt, parameters, param_pt, stability_pt, step_stability);

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
