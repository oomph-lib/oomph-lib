#ifndef OOMPH_CONVERGENCE_DATA_H
#define OOMPH_CONVERGENCE_DATA_H


/*


TODO: This, along with all output from the newton solve could be better
implemented using the observer pattern. That way we only need a single hook
in the newton solve at ouput points and all the "junk" that makes the
algorithm itself harder to understand would go elsewhere. 

Could be a problem using this for general output because debugging could
get misleading (output to stdout would not be immediate). So maybe not that
part...


*/

#include "linear_solver.h"
#include "iterative_linear_solver.h"

using namespace oomph;

namespace oomph
{

 // ============================================================
 /// Struct to store convergence data for a single Newton step.
 // ============================================================
 struct NewtonStepConvergenceData
 {
 public:
  double Max_residual;
  unsigned Iterative_solver_iterations;
  double Iterative_solver_time;
  double Preconditioner_setup_time;
  double Jacobian_setup_time;
 };

 // =================================================================
 /// Class to store convergence data for a complete Problem (i.e. possibly
 /// multiple Newton solves each consistenting of multiple Newton steps).
 // =================================================================
 class ConvergenceData
 {
 public:
  /// Move on to a new Newton solve
  void new_newton_solve()
  {
   Vector<NewtonStepConvergenceData> new_data_vector;
   Solve_data.push_back(new_data_vector);
   Dt_list.push_back(-1.0);
  }

  /// Add all the data for a single Newton step.
  void add_step(const double& max_residual,
                const LinearSolver* const lin_solver_pt,
                const Time* const time_pt)
  {
   // Store dt for this Newton solve. In time-adaptive methods we want to
   // store the final choice of dt so we need to overwrite the value
   // everytime we take a Newton step (which is why this doesn't go in
   // "new_newton_solve()").
   if(time_pt != 0) Dt_list.back() = time_pt->dt();

   // Make a new data structure.
   NewtonStepConvergenceData conv;

   // Fill in the general data.
   conv.Max_residual = max_residual;
   conv.Jacobian_setup_time = lin_solver_pt->jacobian_setup_time();
   conv.Iterative_solver_time = lin_solver_pt->linear_solver_solution_time();

   // If we are using an iterative solver then store iteration data.
   const IterativeLinearSolver* its_pt
    = dynamic_cast<const IterativeLinearSolver*>(lin_solver_pt);
   if(its_pt == 0)
    {
     conv.Iterative_solver_iterations = 0;
     conv.Preconditioner_setup_time = 0;
    }
   else
    {
     conv.Iterative_solver_iterations = its_pt->iterations();
     conv.Preconditioner_setup_time = its_pt->preconditioner_setup_time();
    }

   // Store the data as a new Newton step in the current Newton solve.
   Solve_data.back().push_back(conv);
  }

  /// Get the final residual of Newton solve number isolve.
  double final_residual(const unsigned &isolve) const
  {return Solve_data[isolve].back().Max_residual;}

  /// Get the number of Newton steps taken in Newton solve number isolve.
  unsigned n_newton_steps(const unsigned &isolve) const
  {return Solve_data[isolve].size();}

  /// \short Get the average number of linear solver iterations per Newton
  /// step over the whole problem so far.
  double average_linear_solver_iterations()
  {
   unsigned count = 0;
   unsigned sum = 0;
   for(unsigned i=0; i < Solve_data.size(); i++)
    {
     for(unsigned j=0; j < Solve_data[i].size(); j++)
      {
       count++;
       sum += Solve_data[i][j].Iterative_solver_iterations;
      }
    }
   return double(sum)/count;
  }

  /// Get the average number of Newton steps per Newton solve.
  double average_newton_iterations()
  {
   unsigned sum = 0;
   for(unsigned i=0; i < Solve_data.size(); i++)
    { sum += n_newton_steps(i); }
   return double(sum)/double(Solve_data.size());
  }

  // Output functions
  // ============================================================

  /// Output headers to a stream
  void write_headers(std::ostream &out) const
  {
   out << "dt N_Newton_steps Final_residual N_Iterative_solver_steps"
       << std::endl;
  }

  // Output headers to a file
  void write_headers(const std::string &filename) const
  {
   std::ofstream out(filename.c_str(),std::ios::app);
   write_headers(out);
   out.close();
  }

  /// Output a single Newton step of data
  void output(const unsigned &i, std::ostream &out) const
  {
   out << Dt_list[i] << " " << n_newton_steps(i) << " " << final_residual(i);
   for(unsigned j=0; j < Solve_data[i].size(); j++)
    {
     out << " " << Solve_data[i][j].Iterative_solver_iterations;
    }
   out << std::endl;
  }

  /// Output all Newton steps to a stream.
  void output(std::ostream &out) const
  {
   write_headers(out);
   for(unsigned i=0; i < Solve_data.size(); i++) output(i,out);
  }

  /// Output only most recent Newton step
  void output_this_newton_step(std::ostream &out) const
  {output(Solve_data.size(),out);}

  /// Output only most recent Newton step to a file
  void output_this_newton_step(const std::string &filename) const
  {
   std::ofstream out(filename.c_str(),std::ios::app);
   output_this_newton_step(out);
   out.close();
  }

 private:

  /// \short Storage for all the convergence data (this vector is 2D because we
  /// store data for each Newton step within each Newton solve).
  Vector< Vector<NewtonStepConvergenceData> > Solve_data;

  /// Storage for simulation time at each step
  Vector<double> Dt_list;
 };

} // End of oomph namespace

#endif
