#include <iostream>

#include "generic.h"
#include "ode.h"

using namespace std;
using namespace oomph;

namespace sine_problem
{
  /// ODEFunctor class. Used to create an ODE object with the appropriate set of
  /// properties and functions. Implements the specific sinusoidal ODE to solve
  class ODEFunctor : public SolutionFunctorBase
  {
  public:
    /// Constructor
    ODEFunctor() {}

    /// Destructor
    virtual ~ODEFunctor() {}

    /// Exact or approximate solution. Used for initialisation and/or error
    /// checking
    Vector<double> operator()(const double& t, const Vector<double>& x) const
    {
      Vector<double> output(1);

      //output[0] = sin(4 * atan(1) * t);
      output[0] = exp(- t);

      return output;
    }

    /// Derivative function. Specifies the ODE that we are solving
    Vector<double> derivative(const double& t,
                              const Vector<double>& x,
                              const Vector<double>& u) const
    {
      Vector<double> output(1);

      //output[0] = 4 * atan(1) * cos(4 * atan(1) * t);
      //output[0] = -exp(- t);
      output[0] = -u[0];

      return output;
    }
  };
} // namespace sine_problem

/// ODEProblem class. Specific class to solve the given ODEFunctor.
class ODEProblem : public Problem
{
public:
  /// Constructor
  ODEProblem(const bool& Doc_error) : doc_error(Doc_error)
  {
    /// Create and add time stepper
    bool is_adaptive = true;
    TimeStepper* time_stepper_pt = new BDF<1>(is_adaptive);
    this->add_time_stepper_pt(time_stepper_pt);

    /// Create ODE element
    sine_problem::ODEFunctor* ode_functor_pt = new sine_problem::ODEFunctor;
    ODEElement* ode_element_pt =
      new ODEElement(time_stepper_pt, ode_functor_pt);

    /// Create mesh
    this->mesh_pt() = new Mesh;
    this->mesh_pt()->add_element_pt(ode_element_pt);

    /// Setup equation numbering
    this->assign_eqn_numbers();
    cout << "Number of equations: " << this->ndof() << endl;
  }

  /// Destructor
  ~ODEProblem(){};

  /// Set the initial conditions including previous history
  void set_initial_condition()
  {
    unsigned element_index = 0;
    unsigned data_index = 0;
    unsigned value_index = 0;

    unsigned n_tstorage = this->mesh_pt()
                            ->element_pt(element_index)
                            ->internal_data_pt(data_index)
                            ->time_stepper_pt()
                            ->ntstorage();
    for (unsigned time_level = 0; time_level < n_tstorage; time_level++)
    {
      double t = this->mesh_pt()
                   ->element_pt(element_index)
                   ->internal_data_pt(data_index)
                   ->time_stepper_pt()
                   ->time_pt()
                   ->time(time_level);
      Vector<double> u0 =
        dynamic_cast<ODEElement*>(this->mesh_pt()->element_pt(element_index))
          ->exact_solution(t);

      this->mesh_pt()
        ->element_pt(element_index)
        ->internal_data_pt(data_index)
        ->set_value(time_level, value_index, u0[0]);
    }

    assign_initial_values_impulsive();
  }

  /// Iterate forward in time using adaptive steps until t_final
  void iterate_timestepper(const double& t_step,
                           const double& t_tol,
                           const double& t_final,
                           DocInfo& doc_info)
  {
    double current_t_step = t_step;
    while (this->time_stepper_pt()->time() < t_final)
    {
      current_t_step =
        this->adaptive_unsteady_newton_solve(current_t_step, t_tol);
      this->doc_solution(doc_info);
    }
  }

  /// Must include this function for the adaptive scheme to work
  double global_temporal_error_norm()
  {
    return pow(compute_current_error(), 2.0);
  }

  double compute_current_error()
  {
    double actual =
      this->mesh_pt()->element_pt(0)->internal_data_pt(0)->value(0);

    double t = this->time();

    Vector<double> u = dynamic_cast<ODEElement*>(this->mesh_pt()->element_pt(0))
                         ->exact_solution(t);

    double expected = u[0];

    printf("error: %g\n", actual - expected);

    return actual - expected;
  }
  /// Document the solution
  void doc_solution(DocInfo& doc_info)
  {
    /// Create and open the file stream
    ofstream output_stream;
    string filename = doc_info.directory() + "/soln.dat";
    output_stream.open(filename, ofstream::out | ofstream::app);

    /// Write the current time
    unsigned element_index = 0;
    unsigned data_index = 0;
    double t = this->mesh_pt()
                 ->element_pt(element_index)
                 ->internal_data_pt(data_index)
                 ->time_stepper_pt()
                 ->time();
    output_stream << "t: " << t << ", ";

    /// Write the current value
    unsigned value_index = 0;
    double actual = this->mesh_pt()
                      ->element_pt(element_index)
                      ->internal_data_pt(data_index)
                      ->value(value_index);
    output_stream << "actual: " << actual << ", ";

    /// If we are documenting the error then write the exact value and
    /// difference
    if (this->doc_error)
    {
      Vector<double> u =
        dynamic_cast<ODEElement*>(this->mesh_pt()->element_pt(element_index))
          ->exact_solution(t);
      double expected = u[0];
      output_stream << "expected: " << expected << ", ";

      double error = actual - expected;
      output_stream << "error: " << error << ", ";
    }

    /// End line and close file stream
    output_stream << endl;
    output_stream.close();
  }

private:
  /// Document error flag
  bool doc_error;
};

/// Main function to handle the input arguments, setup the ODE problem and then
/// call the iterator to solve the problem forward in time.
int main(int argc, char* argv[])
{
  /// Parse the command line arguments
  CommandLineArgs::setup(argc, argv);
  /// Command line variables and default values
  string output_directory = "RESLT";
  double initial_t_step = 1e-3;
  double t_tol = 1e-3;
  const string error_flag_string = "--with-error";
  bool has_unrecognised_arg = false;

  CommandLineArgs::specify_command_line_flag(
    "-o",
    &output_directory,
    "Optional: Output directory without the trailing slash (e.g. RESLT )");
  CommandLineArgs::specify_command_line_flag(
    error_flag_string,
    "Optional: Use the exact solution to compute the error.");
  CommandLineArgs::specify_command_line_flag(
    "-s",
    &initial_t_step,
    "Optional: Adjust the initial time step value. Default: 0.1");
  CommandLineArgs::specify_command_line_flag(
    "--tol",
    &t_tol,
    "Optional: Adjust the adaptive timestep time step tolerance. Default: "
    "1e-3");

  CommandLineArgs::parse_and_assign(argc, argv, has_unrecognised_arg);
  /// Parse bool arguments
  bool doc_error =
    CommandLineArgs::command_line_flag_has_been_set(error_flag_string);

  /// Create a DocInfo object
  DocInfo doc_info(output_directory);

  /// Setup the ODE problem
  ODEProblem my_ode_problem(doc_error);

  /// Set the initial condition
  my_ode_problem.set_initial_condition();

  /// Document the initial condition
  my_ode_problem.doc_solution(doc_info);

  /// Iterate forward in time
  const double t_final = 2.0;
  my_ode_problem.iterate_timestepper(initial_t_step, t_tol, t_final, doc_info);

  return 0;
}
