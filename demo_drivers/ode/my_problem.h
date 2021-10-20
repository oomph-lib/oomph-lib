//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================

#ifndef OOMPH_MY_GENERIC_PROBLEM_H
#define OOMPH_MY_GENERIC_PROBLEM_H

#include "generic.h"

#include <ctime>
#include <ostream>
#include <string>

namespace oomph
{

 /// Given a preconditioner:
 /// 1) if it's a the right type of preconditioner return it
 /// 2) otherwise if its a SoM preconditioner containing the right type
 ///    of preconditioner then return a pointer to the underlying
 ///    preconditioner.
 /// 3) otherwise return null
//??ds compatability wrapper, remove?
 template<class T>
  T smart_cast_preconditioner(Preconditioner* prec_pt)
 {
  T bp_pt = dynamic_cast<T> (prec_pt);
  if(bp_pt != 0)
  {
   return bp_pt;
  }
  else
  {
   return 0;
  }
 }

 using namespace StringConversion;

 class ElementalFunction;

 inline std::string real_date_time()
 {
  time_t rawtime;
  struct tm * timeinfo;
  char buffer [80];

  // Get time
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  // Write into cstring
  strftime(buffer, 80, "%Y-%m-%d-%H-%M-%S", timeinfo);

  // Return as "real" string
  return std::string(buffer);
 }

 struct SolverParameters
 {
  SolverParameters()
   {
    linear_solver_pt = 0;
    mass_matrix_solver_for_explicit_timestepper_pt = 0;
   }

  LinearSolver* linear_solver_pt;

  // Newton options
  double newton_solver_tolerance;
  unsigned max_newton_iterations;
  double max_residuals;
  bool shut_up_in_newton_solve;
  bool always_take_one_newton_step;

  // Linear optimisations
  bool jacobian_reuse_is_enabled;
  bool jacobian_has_been_computed;
  bool problem_is_nonlinear;

  // Explicit solver
  LinearSolver* mass_matrix_solver_for_explicit_timestepper_pt;
  bool mass_matrix_reuse_is_enabled;
  bool mass_matrix_has_been_computed;
  bool discontinuous_element_formulation;
 };


 class MyDocInfo : public DocInfo
 {
 public:
  /// Default constructor
 MyDocInfo() : DocInfo(), output_jacobian("never") {}

  /// Copy dump of args into args_str.
  void copy_args_string()
  {
   std::ostringstream stream;
   CommandLineArgs::doc_all_flags(stream);
   args_str.assign(stream.str());
  }

  std::string output_jacobian;
  std::string args_str;
 };


 class MyProblem : public Problem
 {
 public:
  /// Default constructor
 MyProblem() :
  Trace_filename("trace"),
   Info_filename("info"),
   Trace_seperator("; "),
   Dummy_doc_data(-1)
   {
    Dim = 0;
    Output_precision = 8;
    Error_norm_limit = -1.0;
    Solution_norm_limit = -1.0;

    // Throw a real error (not just a warning) if the output directory
    // does not exist.
    Doc_info.enable_error_if_directory_does_not_exist();

    Disable_mass_matrix_solver_optimisations = false;

    // By default output to trace file every step
    Always_write_trace = true;

    Dump = false;
    Output_ltes = false;
    Output_predictor_values = false;
    Want_doc_exact = false;
    Output_initial_condition = false;
    Should_doc_boundaries = false;

    N_steps_taken = 0;
    Total_step_time= 0;
   }

  /// Virtual destructor. Policy decision: my problem classes won't call
  /// delete on anything. That's up to the driver code or
  /// whatever. Ideally we should just start using c++11 smart pointers
  /// and not have to worry so much about memory management!
  virtual ~MyProblem() {}

  /// Do the newton solve or explicit step (different ones depending flags
  /// set).
  double smart_time_step(double dt, const double& tol)
  {
   double step_time_start = TimingHelpers::timer();

   // The Newton step itself, adaptive if requested.
   if(explicit_flag())
   {
    explicit_timestep(dt);
   }
   else if(tol != 0.0)
   {
    dt = adaptive_unsteady_newton_solve(dt, tol);
   }
   else
   {
    unsteady_newton_solve(dt);
   }

   double step_time_stop = TimingHelpers::timer();
   Total_step_time = step_time_stop - step_time_start;


   oomph_info << "Time for step " << Total_step_time
	      << std::endl;

   N_steps_taken++;

   return dt;
  }

  /// Use to specify a condition for time stepping to halt early. By
  /// default never halt early.
  virtual bool finished() const
  {
   return false;
  }

  void get_solver_parameters(SolverParameters& sp)
  {
   sp.linear_solver_pt = linear_solver_pt();

   // Newton options
   sp.newton_solver_tolerance = newton_solver_tolerance();
   sp.max_newton_iterations = max_newton_iterations();
   sp.max_residuals = max_residuals();
   sp.shut_up_in_newton_solve = Shut_up_in_newton_solve;
   sp.always_take_one_newton_step = Always_take_one_newton_step;

   // Linear optimisations
   sp.jacobian_reuse_is_enabled = jacobian_reuse_is_enabled();
   sp.jacobian_has_been_computed = Jacobian_has_been_computed;
   sp.problem_is_nonlinear = Problem_is_nonlinear;

   // Explicit solver
   sp.mass_matrix_solver_for_explicit_timestepper_pt = mass_matrix_solver_for_explicit_timestepper_pt();
   sp.mass_matrix_reuse_is_enabled = mass_matrix_reuse_is_enabled();
   sp.mass_matrix_has_been_computed = Mass_matrix_has_been_computed;
   sp.discontinuous_element_formulation = Discontinuous_element_formulation;
  }

  void set_solver_parameters(SolverParameters& sp)
  {
   linear_solver_pt() = sp.linear_solver_pt;

   // Newton options
   newton_solver_tolerance() = sp.newton_solver_tolerance;
   max_newton_iterations() = sp.max_newton_iterations;
   max_residuals() = sp.max_residuals;
   Shut_up_in_newton_solve = sp.shut_up_in_newton_solve;
   Always_take_one_newton_step = sp.always_take_one_newton_step;

   // Linear optimisations
   Jacobian_reuse_is_enabled = sp.jacobian_reuse_is_enabled;
   Jacobian_has_been_computed = sp.jacobian_has_been_computed;
   Problem_is_nonlinear = sp.problem_is_nonlinear;

   // Explicit solver
   mass_matrix_solver_for_explicit_timestepper_pt() = sp.mass_matrix_solver_for_explicit_timestepper_pt;
   Mass_matrix_reuse_is_enabled = sp.mass_matrix_reuse_is_enabled;
   Mass_matrix_has_been_computed = sp.mass_matrix_has_been_computed;
   Discontinuous_element_formulation = sp.discontinuous_element_formulation;
  }

  bool explicit_flag()
  {
   return (explicit_time_stepper_pt() != 0)
    && (time_stepper_pt()->is_steady());
  }

  bool is_steady()
  {
   return (explicit_time_stepper_pt() == 0)
    && (time_stepper_pt()->is_steady());
  }

  virtual void actions_before_newton_step()
  {
   // Output Jacobian and residuals if requested
   if(to_lower(Doc_info.output_jacobian) == "always")
   {
    // label with doc_info number and the newton step number
    std::string label = to_string(Doc_info.number())
     + "_"
     + to_string(nnewton_step_this_solve() + 1);

    dump_current_mm_or_jacobian_residuals(label);
   }
  }

  unsigned nnewton_step_this_solve() const
  {
   return Jacobian_setup_times.size();
  }

  virtual void actions_before_explicit_stage()
  {MyProblem::actions_before_newton_step();}

  virtual void actions_after_explicit_stage()
  {
   Jacobian_setup_times.push_back
    (this->mass_matrix_solver_for_explicit_timestepper_pt()->jacobian_setup_time());
   Solver_times.push_back
    (this->mass_matrix_solver_for_explicit_timestepper_pt()->linear_solver_solution_time());

   // No non-linear residuals to store

   const IterativeLinearSolver* its_pt
    = dynamic_cast<const IterativeLinearSolver*>(this->mass_matrix_solver_for_explicit_timestepper_pt());
   if(its_pt != 0)
   {
    Solver_iterations.push_back(its_pt->iterations());
    Preconditioner_setup_times.push_back(its_pt->preconditioner_setup_time());
   }
   else
   {
    // Fill in dummy data
    Solver_iterations.push_back(Dummy_doc_data);
    Preconditioner_setup_times.push_back(Dummy_doc_data);
   }

   // Not quite the same as actions after newton step because we are
   // interested in what happened in the explicit solver instead of the
   // main solver.
  }

  virtual void actions_before_time_integration() {}

  virtual void actions_before_explicit_timestep()
  {
   MyProblem::actions_before_newton_solve();
   check_norm_limits();
  }

  virtual void actions_after_explicit_timestep()
  {
   MyProblem::actions_after_newton_solve();
  }

  virtual void actions_after_implicit_timestep() {}

  virtual void actions_before_implicit_timestep()
  {
   check_norm_limits();
  }


  virtual void actions_after_newton_step()
  {
   Jacobian_setup_times.push_back
    (this->linear_solver_pt()->jacobian_setup_time());
   Solver_times.push_back
    (this->linear_solver_pt()->linear_solver_solution_time());

   const IterativeLinearSolver* its_pt
    = dynamic_cast<const IterativeLinearSolver*>(this->linear_solver_pt());
   if(its_pt != 0)
   {
    Solver_iterations.push_back(its_pt->iterations());
    Preconditioner_setup_times.push_back(its_pt->preconditioner_setup_time());
   }
   else
   {
    // Fill in dummy data
    Solver_iterations.push_back(Dummy_doc_data);
    Preconditioner_setup_times.push_back(Dummy_doc_data);
   }
  }

  virtual void actions_before_newton_solve()
  {
   // Clean up times vectors
   Jacobian_setup_times.clear();
   Solver_times.clear();
   Solver_iterations.clear();
   Preconditioner_setup_times.clear();
  }



  /// Pin all dofs with index in node one of indices in all nodes. Uses a
  /// different magic pinning number so that it can be easily undone
  /// using undo_segregated_pinning(). Does not handle: global data,
  /// element data, nodes with varying nvalue.
  void segregated_pin_indices(const Vector<unsigned>& indices)
  {
   for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
   {
    Mesh* mesh_pt = this->mesh_pt(msh);
    for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
    {
     Node* nd_pt = mesh_pt->node_pt(nd);
     for(unsigned j=0; j<indices.size(); j++)
     {
      if(!nd_pt->is_pinned(indices[j]))
      {
       nd_pt->eqn_number(indices[j])
	= Data::Is_segregated_solve_pinned;
      }
     }
    }
   }
   // Output information
   oomph_info << "Segregated solve, without indices:\n";

   // Loop over the entries of indices and output
   unsigned indices_length=indices.size();
   if (indices_length==0)
   {
    oomph_info << Trace_seperator << "[]";
   }
   else
   {
    oomph_info << Trace_seperator << "[" << indices[0];
    if (indices_length>1)
    {
     for (unsigned i=1;i<indices_length;i++)
     {
      oomph_info << ", " << indices[i];
     }
    }
    oomph_info << "]";
   }
  
   // Output the number of equations
   oomph_info << "Number of equations: " << assign_eqn_numbers()  << std::endl;
  }

  /// Remove pinning set up by segregated_pin_indices.
  void undo_segregated_pinning()
  {
   for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
   {
    Mesh* mesh_pt = this->mesh_pt(msh);
    for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
    {
     Node* nd_pt = mesh_pt->node_pt(nd);
     for(unsigned j=0; j<nd_pt->nvalue(); j++)
     {
      if(nd_pt->eqn_number(j) == Data::Is_segregated_solve_pinned)
      {
       nd_pt->eqn_number(j) = Data::Is_unclassified;
      }
     }
    }
   }
   oomph_info << "un-segregated n eqn " << assign_eqn_numbers() << std::endl;
  }


  /// Check that nothing is currently pinned for a segregated solve.
  void check_not_segregated(const char* function) const
  {
#ifdef PARANOID
   for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
   {
    Mesh* mesh_pt = this->mesh_pt(msh);
    for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
    {
     Node* nd_pt = mesh_pt->node_pt(nd);
     for(unsigned j=0; j<nd_pt->nvalue(); j++)
     {
      if(nd_pt->eqn_number(j) == Data::Is_segregated_solve_pinned)
      {
       throw OomphLibError("Some dofs already segregated",
			   OOMPH_EXCEPTION_LOCATION,
			   function);
      }
     }
    }
   }
#endif
  }


  void check_norm_limits()
  {
   // If a limit has been set
   if(Error_norm_limit != -1.0)
   {
    double error_norm = get_error_norm();

    if((error_norm != Dummy_doc_data)
       && (error_norm > Error_norm_limit))
    {
     std::string err = "Error norm " + to_string(error_norm);
     err += " exceeds the limit " + to_string(Error_norm_limit);
     throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
			 OOMPH_CURRENT_FUNCTION);
    }
   }

   // Same for solution norm
   if(Solution_norm_limit != -1.0)
   {
    double solution_norm = get_solution_norm();

    if((solution_norm != Dummy_doc_data)
       && (solution_norm > Solution_norm_limit))
    {
     std::string err = "Solution norm " + to_string(solution_norm);
     err += " exceeds the limit " + to_string(Solution_norm_limit);
     throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
			 OOMPH_CURRENT_FUNCTION);
    }
   }
  }

  /// ??ds
  double min_element_size();

  ///  Write some general data about the previous time step to a
  /// trace file. Extend by overloading write_additional_trace_data(...).
  void write_trace(const unsigned& t_hist=0);


  ///  Overload to write problem specific data into trace
  /// file. Note: don't add any line endings, seperate data with the
  /// Trace_seperator string (BEFORE each data entry).
  virtual void write_additional_trace_data(const unsigned& t_hist,
					   std::ofstream& trace_file) const {}

  ///  Overload to write problem specific headers into trace
  /// file. Note: don't add any line endings, seperate headers with the
  /// Trace_seperator string (BEFORE each data entry).
  virtual void write_additional_trace_headers(std::ofstream& trace_file)
   const {}

  /// Overload to write any problem specific data
  virtual void output_solution(std::ofstream& soln_file) const {}
  virtual void final_doc_additional() const {}
  virtual void initial_doc_additional() const {}


  ///  Outputs to be done at the start of a run (i.e. outputing
  /// basic info on command line args etc, writing trace file headers,
  /// outputting initial conditions).
  void initial_doc();

  ///  Outputs to be done at the end of a run (e.g. closing tags
  /// for XML files).
  void final_doc()
  {
   // Output Jacobian if requested
   if((to_lower(Doc_info.output_jacobian) == "at_end")
      || (to_lower(Doc_info.output_jacobian) == "always"))
   {
    dump_current_mm_or_jacobian_residuals("at_end");
   }

   // Write end of .pvd XML file
   std::ofstream pvd_file((Doc_info.directory() + "/" + "soln.pvd").c_str(),
			  std::ios::app);
   pvd_file << "</Collection>" << std::endl
	    << "</VTKFile>" << std::endl;
   pvd_file.close();

   // Write end of exact.pvd XML file
   if(doc_exact())
   {
    std::ofstream pvd_file((Doc_info.directory() + "/" + "exact.pvd").c_str(),
			   std::ios::app);
    pvd_file << "</Collection>" << std::endl
	     << "</VTKFile>" << std::endl;
    pvd_file.close();
   }

   // Write out anything requested from an inheriting class.
   final_doc_additional();
  }

  ///  General output function: output to trace file. Maybe output
  /// Jacobian depending on Doc_info.output_jacobian. Maybe output full
  /// solution depending on should_doc_this_step(..) function. Maybe
  /// output ltes depending on value of Output_ltes. Extend by
  /// overloading output_solution(...). Optional prefix for special
  /// outputs (special outputs don't go in pvd file, yet? too messy...)
  void doc_solution(const unsigned& t_hist=0,
		    const std::string& prefix = "" );


  /// Standard output function: loop over all elements in all meshes and
  /// output.
  virtual void output_solution(const unsigned& t, std::ostream& outstream,
			       const unsigned& npoints=2) const
  {
   const unsigned n_msh = nsub_mesh();
   for(unsigned msh=0; msh<n_msh; msh++)
   {
    Mesh* msh_pt = mesh_pt(msh);

    const unsigned n_ele = msh_pt->nelement();
    for(unsigned ele=0; ele<n_ele; ele++)
    {
     FiniteElement* ele_pt = msh_pt->finite_element_pt(ele);
     ele_pt->output(t, outstream, npoints);
    }
   }
  }

  /// Output nodes with boundary number to csv file
  virtual void doc_boundaries(const std::string& boundary_file_basename) const;

  /// output_solution(...) with default output time step = 0 = current
  /// time.
  void output_solution(std::ostream& outstream,
		       const unsigned& npoints=2) const
  {
   output_solution(0, outstream, npoints);
  }

  /// Standard output function: loop over all elements in all meshes and
  /// output exact solution.
  virtual void output_exact_solution(const unsigned& t_hist,
				     std::ostream& outstream,
				     const unsigned& npoints=2) const
  {
   const double time = time_pt()->time();

   const unsigned n_msh = nsub_mesh();
   for(unsigned msh=0; msh<n_msh; msh++)
   {
    Mesh* msh_pt = mesh_pt(msh);

    const unsigned n_ele = msh_pt->nelement();
    for(unsigned ele=0; ele<n_ele; ele++)
    {
     FiniteElement* ele_pt = msh_pt->finite_element_pt(ele);
     ele_pt->output_fct(outstream, npoints, time,
			*Exact_solution_pt);
    }
   }
  }

  ///  Error norm calculator
  virtual double get_error_norm(const unsigned& t_hist=0) const
  {
   if(Exact_solution_pt != 0)
   {
    // ExactFunctionDiffSquared f;
    // f.Exact_pt = Exact_solution_pt;
    // return std::sqrt(integrate_over_problem(&f));

    // Nodal rms difference
    const double t = time_pt()->time(t_hist);

    double diffsq = 0;

    const unsigned n_node = mesh_pt()->nnode();
    for(unsigned nd=0; nd<n_node; nd++)
    {
     Node* nd_pt = mesh_pt()->node_pt(nd);
     Vector<double> values(nd_pt->nvalue(), 0.0), x(dim(), 0.0);
     nd_pt->position(t_hist, x);
     nd_pt->value(t_hist, values);

     Vector<double> exact = exact_solution(t, x);

     const unsigned ni = values.size();
     for(unsigned i=0; i<ni; i++)
     {
      diffsq += std::pow(values[i] - exact[i], 2);
     }
    }

    return std::sqrt(diffsq);
   }
   else
   {
    return Dummy_doc_data;
   }
  }

  ///  Dummy solution norm calculator (overload in derived classes).
  virtual double get_solution_norm(const unsigned& t_hist=0) const
  {
   DoubleVector dofs;
   get_dofs(t_hist, dofs);
   return dofs.norm();
  }

  ///  should the previous step be doc'ed? Check if we went past an
  /// entry of Doc_times in the last step. If no Doc_times have been set
  /// then always output.
  virtual bool should_doc_this_step(const double &dt, const double &time) const
  {
   // I'm sure there should be a more efficient way to do this if we
   // know that Doc_times is ordered, but it shouldn't really matter I
   // guess--Jacobian calculation and solve times will always be far
   // far larger than this brute force search.

   // If no Doc_times have been set then always output.
   if(Doc_times.empty()) return true;

   // Loop over entries of Doc_times and check if they are in the
   // range (t - dt, t].
   for(unsigned j=0; j<Doc_times.size(); j++)
   {
    if(( time >= Doc_times[j]) && ((time - dt) < Doc_times[j]))
    {
     return true;
    }
   }
   return false;
  }

  ///  Assign a vector of times to output the full solution at.
  void set_doc_times(Vector<double> &doc_times)
  {
   Doc_times = doc_times;
  }


  /// Get an lte error norm using the same norm and values as the
  /// adaptive time stepper used.
  double lte_norm()
  {
   if(time_stepper_pt()->adaptive_flag())
   {
    // Just return the error that would be used for adaptivity.
    return global_temporal_error_norm();
   }
   else
   {
    return Dummy_doc_data;
   }
  }

  bool doc_exact() const
  {
   return Want_doc_exact && (Exact_solution_pt != 0);
  }

  virtual Vector<double> trace_values(const unsigned& t_hist=0) const
  {
   unsigned nele = mesh_pt()->nelement();
   unsigned e = nele/2;

   Vector<double> values;
   if(dynamic_cast<FiniteElement*>(mesh_pt()->element_pt(e)))
   {
    // Just use an element somewhere in the middle...
    Node* nd_pt = mesh_pt()->finite_element_pt(e)->node_pt(0);
    values.assign(nd_pt->nvalue(), 0.0);
    nd_pt->value(t_hist, values);
   }
   else
   {
    // Not finite elements so no idea what to use
    values.assign(1, Dummy_doc_data);
   }

   return values;
  }


  void dump_current_mm_or_jacobian_residuals(const std::string& label);


  IterativeLinearSolver* iterative_linear_solver_pt() const
  {
   return dynamic_cast<IterativeLinearSolver*>
    (this->linear_solver_pt());
  }


  ///  Perform set up of problem.
  virtual void build(Vector<Mesh*>& bulk_mesh_pt);

  ///  Get problem dimension (nodal dimension).
  const unsigned dim() const {return this->Dim;}

  virtual std::string problem_name() const {return "unknown";}

  /// Set all history values/dts to be the same as the present values/dt.
  virtual void set_up_impulsive_initial_condition();

  /// Assign initial conditions from function pointer
  virtual void my_set_initial_condition(const SolutionFunctorBase& ic);

  /// Hook to be overloaded with any calculations needed after setting of
  /// initial conditions.
  virtual void actions_after_set_initial_condition();

  ///  Integrate a function given by func_pt over every element
  /// in every bulk mesh in this problem.
  virtual double integrate_over_problem(const ElementalFunction* func_pt,
					const Integral* quadrature_pt=0) const;


  virtual void dump(std::ofstream& dump_file) const
  {
   // Set very high precision to avoid any issues
   dump_file.precision(14);

   dump_file << Doc_info.number() << " # Doc_info.number()" << std::endl;
   dump_file << N_steps_taken << " # N_steps_taken" << std::endl;
   Problem::dump(dump_file);
  }

  virtual void read(std::ifstream& restart_file)
  {
   // buffer
   std::string input_string;

   // Read in Doc_info number. Read line up to termination sign then
   // ignore.
   getline(restart_file, input_string, '#');
   restart_file.ignore(80,'\n');
   Doc_info.number() = std::atoi(input_string.c_str());

   // Read in number of steps taken. Read line up to termination sign
   // then ignore.
   getline(restart_file, input_string, '#');
   restart_file.ignore(80,'\n');
   N_steps_taken = std::atoi(input_string.c_str());

   // Let base class handle the rest
   Problem::read(restart_file);

   // Decrement doc info number so that it is correct after initial doc
   // Doc_info.number()--;
  }


  /// Output lte values of each nodal value along with spatial position
  /// for plotting with paraview. To plot load the csv file, use "Table
  /// To Points", add a 3D view, set the ltes to visible and color by the
  /// lte of interest. Useful to set "Representation" to "Points" and
  /// increase point size so that things are easily visible. Not
  /// implemented for nodes with varying numbers of values (you might
  /// need something fancier than a csv file for this).
  void output_ltes(const unsigned& t_hist, std::ostream& output) const
  {
   // Output position labels
   for(unsigned j=0; j<Dim; j++)
   {
    output << "x" << j << ", ";
   }

   // Output labels for ltes, assuming that all nodes have the same
   // number of values...
   for(unsigned j=0; j<mesh_pt()->node_pt(0)->nvalue(); j++)
   {
    output << "error" << j << ", ";
   }

   output << std::endl;

   // Output actual positions and ltes
   for(unsigned i=0, ni=mesh_pt()->nnode(); i<ni; i++)
   {
    Node* nd_pt = mesh_pt()->node_pt(i);

    // Output position of node
    for(unsigned j=0; j<Dim; j++)
    {
     output << nd_pt->x(j) << ", ";
    }

    // Output ltes of node
    for(unsigned j=0; j<nd_pt->nvalue(); j++)
    {
     // Get timestepper's error estimate for this direction of m
     // at this point.
     double error = nd_pt->time_stepper_pt()->
      temporal_error_in_value(nd_pt, j);

     output << error << ", ";
    }

    output << std::endl;
   }
  }


  MyDocInfo Doc_info;
  unsigned Output_precision;

  unsigned N_steps_taken;
  double Total_step_time;

  std::string Trace_filename;
  std::string Info_filename;

  double Error_norm_limit;
  double Solution_norm_limit;

  /// Option to turn off optimisation of the linear solves needed for
  /// explicit timestepping (for debugging purposes).
  bool Disable_mass_matrix_solver_optimisations;

  /// Should we output to trace file every step?
  bool Always_write_trace;

  /// Should we try to output exact solution?
  bool Want_doc_exact;

  /// Should we output the locations of the boundary nodes?
  bool Should_doc_boundaries;

  /// Should we dump ready for a restart?
  bool Dump;

  /// Should we output the local truncation error at each node as well?
  bool Output_ltes;

  /// Should we output the predicted values too?
  bool Output_predictor_values;

  /// Should we write output for initial conditions as though they are time steps as well?
  bool Output_initial_condition;

  /// Function pointer for exact solution
  SolutionFunctorBase* Exact_solution_pt;

  /// Get exact solution
  Vector<double> exact_solution(const double& t, const Vector<double>& x) const
  {
#ifdef PARANOID
   if(Exact_solution_pt == 0)
   {
    std::string err = "Exact_solution_pt is null!";
    throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
			OOMPH_CURRENT_FUNCTION);
   }
#endif
   return (*Exact_solution_pt)(t, x);
  }

 protected:

  ///  String to insert between fields in trace file. Use "; " by
  /// default (so that "," or " " can be used for lists if needed).
  const std::string Trace_seperator;
  double Dummy_doc_data;
  unsigned Dim;

 private:

  Vector<double> Jacobian_setup_times;
  Vector<double> Solver_times;
  Vector<double> Solver_iterations;
  Vector<double> Preconditioner_setup_times;

  /// Times at which we want to output the full solution.
  Vector<double> Doc_times;

  /// Inaccessible copy constructor
  MyProblem(const MyProblem &dummy)
  {BrokenCopy::broken_copy("MyProblem");}

  /// Inaccessible assignment operator
  void operator=(const MyProblem &dummy)
   {BrokenCopy::broken_assign("MyProblem");}
 };

 ///  Integrate a function given by func_pt over every element
 /// in every bulk mesh in this problem.
 double MyProblem::integrate_over_problem(const ElementalFunction* func_pt,
					  const Integral* quadrature_pt) const
 {
  throw OomphLibError("Not implemented (yet?).", OOMPH_CURRENT_FUNCTION, 
		      OOMPH_EXCEPTION_LOCATION);

 }


 void MyProblem::dump_current_mm_or_jacobian_residuals(const std::string& label)
 {
  throw OomphLibError("Not implemented (yet?).", OOMPH_CURRENT_FUNCTION, 
		      OOMPH_EXCEPTION_LOCATION);
 }

 ///  Perform set up of problem.
 void MyProblem::build(Vector<Mesh*>& bulk_mesh_pt)
 {
  // Copy the first mesh's first timestepper to the problem

  FiniteElement* fele_pt = dynamic_cast<FiniteElement*>
   (bulk_mesh_pt[0]->element_pt(0));

  // Finite element mesh: grab ts from node
  if(fele_pt != 0)
  {
   TimeStepper* ts_pt = bulk_mesh_pt[0]->node_pt(0)->time_stepper_pt();
   this->add_time_stepper_pt(ts_pt);

   // ??ds assumed any timesteppers hiding somewhere else are added elsewhere

#ifdef PARANOID
   for(unsigned j=0; j<bulk_mesh_pt.size(); j++)
   {
    if(bulk_mesh_pt[j]->node_pt(0)->time_stepper_pt()
       != ts_pt)
    {
     std::string err = "Multiple timesteppers, you need to do somedhing more fancy here";
     throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
			 OOMPH_CURRENT_FUNCTION);
    }
   }
#endif
  }

  // Non finite element mesh: grab ts from internal data
  else
  {
   TimeStepper* ts_pt = bulk_mesh_pt[0]->element_pt(0)->
    internal_data_pt(0)->time_stepper_pt();
   this->add_time_stepper_pt(ts_pt);

   // ??ds again assumed any timesteppers hiding somewhere else are added elsewhere

#ifdef PARANOID
   for(unsigned j=0; j<bulk_mesh_pt.size(); j++)
   {
    if(bulk_mesh_pt[j]->element_pt(0)->
       internal_data_pt(0)->time_stepper_pt() != ts_pt)
    {
     std::string err = "Multiple timesteppers? you need to do somedhing more fancy here";
     throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
			 OOMPH_CURRENT_FUNCTION);
    }
   }
#endif
  }


  // Push all the meshes into the problem's sub mesh list
  for(unsigned j=0; j<bulk_mesh_pt.size(); j++)
  {
   add_sub_mesh(bulk_mesh_pt[j]);
  }

  // If we have an iterative solver with a block preconditioner then
  // add all the meshes to the block preconditioner as well.
  IterativeLinearSolver* its_pt = iterative_linear_solver_pt();
  if(its_pt != 0)
  {
// RRR add comments to why this is incorrect:
// We cannot call set_nmesh and set_mesh, this is handled by the derived
// classes.
   // Try to get a block preconditioner from the preconditioner
   BlockPreconditioner<CRDoubleMatrix>* bp_pt
    = smart_cast_preconditioner<BlockPreconditioner<CRDoubleMatrix>*>
    (its_pt->preconditioner_pt());

   if(bp_pt != 0)
   {
#ifdef PARANOID
    {
     std::string err = "IS THIS EVER USED?";
     throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
			 OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // This part of the code is never executed in the driver code.
    // Commented out since it no longer make sense. The functions
    // set_nmesh() and set_mesh() are made protected. Each derived
    // block preconditioner has to call set_nmesh() and set_mesh(), 
    // since the preconditioner knows which mesh goes where, the 
    // writer of the driver code does not necessarily know.

//            // Set up meshes
//            bp_pt->set_nmesh(nsub_mesh());
//            for(unsigned i=0; i< nsub_mesh(); i++)
//              {
//                bp_pt->set_mesh(i, mesh_pt(i));
//              }
   }
  }

  // Get the problem dimension
  if(fele_pt != 0)
  {
   Dim = fele_pt->nodal_dimension();
  }
  else
  {
   // Presumably if a "bulk" mesh contains non-finite elements
   // then this is not a pde, so no dimension as such.
   Dim = 0;
  }

  if(!Disable_mass_matrix_solver_optimisations)
  {
   // Set the solver for explicit timesteps (mass matrix) to CG with a
   // diagonal predconditioner.
   IterativeLinearSolver* expl_solver_pt = new CG<CRDoubleMatrix>;
   expl_solver_pt->preconditioner_pt() =
    new MatrixBasedLumpedPreconditioner<CRDoubleMatrix>;

   // If it takes more than 100 iterations then something has almost
   // certainly gone wrong!
   expl_solver_pt->max_iter() = 100;
   expl_solver_pt->enable_error_after_max_iter();
   mass_matrix_solver_for_explicit_timestepper_pt() = expl_solver_pt;

   // expl_solver_pt->enable_doc_convergence_history();

   // Store + re-use the mass matrix used in explicit steps (since we
   // are almost certainly not going to do spatially adaptivity
   // anytime soon this is safe).
   this->enable_mass_matrix_reuse();
  }


  // If we requested exact solution output then check we have a solution
  if(Want_doc_exact && Exact_solution_pt == 0)
  {
   std::string warning = "Requested doc'ing exact solution, but we don't have an exact solution function pointer.";
   OomphLibWarning(warning, OOMPH_CURRENT_FUNCTION,
		   OOMPH_EXCEPTION_LOCATION);

  }
 }


 /// Hook to be overloaded with any calculations needed after setting of
 /// initial conditions.
 void MyProblem::actions_after_set_initial_condition()
 {
  // If using TR calculate initial derivative with these initial conditions
  TR* tr_pt = dynamic_cast<TR*>(time_stepper_pt());
  if(tr_pt != 0)
  {
   tr_pt->setup_initial_derivative(this);
  }
 }

 void MyProblem::initial_doc()
 {
#ifdef PARANOID
  if(*(Doc_info.directory().end()-1) == '/')
  {
   std::string error_msg = "Don't put a / on the end of results dir";
   throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
		       OOMPH_EXCEPTION_LOCATION);
  }
#endif

  const std::string& dir = Doc_info.directory();

  // Output Jacobian if requested
  if((to_lower(Doc_info.output_jacobian) == "at_start")
     || (to_lower(Doc_info.output_jacobian) == "always"))
  {
   dump_current_mm_or_jacobian_residuals("at_start");
  }

  // pvd files
  // ============================================================
  // Write start of .pvd XML file
  std::ofstream pvd_file((dir + "/" + "soln.pvd").c_str(),
			 std::ios::out);
  pvd_file << "<?xml version=\"1.0\"?>" << std::endl
	   << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"
	   << std::endl
	   << "<Collection>" << std::endl;
  pvd_file.close();

  // Write start of exact.pvd XML file
  if(doc_exact())
  {
   std::ofstream pvd_file((dir + "/" + "exact.pvd").c_str(),
			  std::ios::out);
   pvd_file << "<?xml version=\"1.0\"?>" << std::endl
	    << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"
	    << std::endl
	    << "<Collection>" << std::endl;
   pvd_file.close();
  }

  // Trace file
  // ============================================================

  // Clear (by overwriting) and write headers
  std::ofstream trace_file((dir + "/" + Trace_filename).c_str());
  trace_file
   << "DocInfo_numbers"
   << Trace_seperator << "times"
   << Trace_seperator << "dts"
   << Trace_seperator << "error_norms"

   << Trace_seperator << "n_newton_iters"
   << Trace_seperator << "n_solver_iters"

   << Trace_seperator << "solver_times"
   << Trace_seperator << "jacobian_setup_times"
   << Trace_seperator << "preconditioner_setup_times"

   << Trace_seperator << "LTE_norms"
   << Trace_seperator << "trace_values"

   << Trace_seperator << "unix_timestamp"
   << Trace_seperator << "newton_residuals"
   << Trace_seperator << "solution_norms"
   << Trace_seperator << "total_step_time"


   // Reserved slots in case I think of more things to add later
   << Trace_seperator << "dummy"
   << Trace_seperator << "dummy"
   << Trace_seperator << "dummy"
   << Trace_seperator << "dummy"
   << Trace_seperator << "dummy"
   << Trace_seperator << "dummy";

  // Other headers depending on the derived class
  write_additional_trace_headers(trace_file);

  // Finish the line and close
  trace_file << std::endl;
  trace_file.close();


  // Info file
  // ============================================================
  std::ofstream info_file((dir + "/" + Info_filename).c_str());
  info_file
   << "real_time " << real_date_time() << std::endl
   << "unix_time " << std::time(0) << std::endl
   << "driver_name " << problem_name() << std::endl
   << "initial_nnode " << mesh_pt()->nnode() << std::endl
   << "initial_nelement " << mesh_pt()->nelement() << std::endl
   << "initial_nsub_mesh " << nsub_mesh() << std::endl;

  info_file << Doc_info.args_str;
  info_file.close();

  // If requested then output history values before t=0.
  if(Output_initial_condition)
  {
#ifdef PARANOID
   if(ntime_stepper() != 1)
   {
    std::string err = "More/less that 1 time stepper, not sure what to output.";
    throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
			OOMPH_EXCEPTION_LOCATION);
   }
#endif

   // Output info for each history value in order.
   const unsigned nval = time_stepper_pt()->nprev_values();
   for(unsigned it=nval-1; it>0; it--)
   {
    doc_solution(it);
   }
  }

  // Output boundary numbers for each boundary node
  // ============================================================
  if(Should_doc_boundaries)
  {
   this->doc_boundaries(dir + "/nodes_on_boundary");
  }


  // Write initial solution and anything else problem specific
  // (e.g. more trace file headers)
  initial_doc_additional();

  // Doc the initial condition
  this->doc_solution();
 }

 void MyProblem::doc_solution(const unsigned& t_hist,
			      const std::string& prefix)
 {
  bool doc_this_step = true;
  if(!is_steady())
  {
   doc_this_step = should_doc_this_step(time_pt()->dt(t_hist),
					time_pt()->time(t_hist));
  }

  const std::string dir = Doc_info.directory();
  const std::string num = to_string(Doc_info.number());

  if(Always_write_trace || doc_this_step)
  {
   // Always output trace file data
   write_trace(t_hist);
  }

  // Output full set of data if requested for this timestep
  if(doc_this_step)
  {

   // Solution itself
   std::ofstream soln_file((dir + "/" + prefix + "soln" + num + ".dat").c_str(),
			   std::ios::out);
   soln_file.precision(Output_precision);
   output_solution(t_hist, soln_file);
   soln_file.close();

   // Exact solution if available and requested
   if(doc_exact())
   {
    std::ofstream exact_file((dir + "/" + prefix + "exact" + num + ".dat").c_str(),
			     std::ios::out);
    exact_file.precision(Output_precision);
    output_exact_solution(t_hist, exact_file);
    exact_file.close();
   }

   // If not a steady state problem then write time information
   if(!is_steady() && prefix == "")
   {
    // Write the simulation time and filename to the solution pvd
    // file
    std::ofstream pvd_file((dir + "/" + "soln.pvd").c_str(),
			   std::ios::app);
    pvd_file.precision(Output_precision);

    pvd_file << "<DataSet timestep=\"" << time_pt()->time(t_hist)
	     << "\" group=\"\" part=\"0\" file=\"" << "soln"
	     << num << ".vtu"
	     << "\"/>" << std::endl;

    pvd_file.close();


    // Write the simulation time and filename to the exact solution
    // pvd file
    if(doc_exact())
    {
     std::ofstream exact_pvd_file((dir + "/" + "exact.pvd").c_str(),
				  std::ios::app);
     exact_pvd_file.precision(Output_precision);

     exact_pvd_file << "<DataSet timestep=\"" << time()
		    << "\" group=\"\" part=\"0\" file=\"" <<  "exact"
		    << num << ".vtu"
		    << "\"/>" << std::endl;

     exact_pvd_file.close();
    }
   }


   // Maybe dump the restart data
   if(Dump)
   {
    std::ofstream dump_file((dir + "/" + "dump" + num + ".dat").c_str(),
			    std::ios::out);
    this->dump(dump_file);
    dump_file.close();
   }

   // Maybe dump ltes
   if(Output_ltes)
   {
    std::ofstream ltefile((dir + "/ltes" + num + ".csv").c_str(),
			  std::ios::out);
    output_ltes(t_hist, ltefile);
    ltefile.close();
   }

   // Maybe write out predicted values, check we only have one time
   // stepper first though.
#ifdef PARANOID
   if(Output_predictor_values && ntime_stepper() != 1)
   {
    std::string err = "Can only output predictor values for a single time stepper";
    throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
			OOMPH_EXCEPTION_LOCATION);

    // Otherwise we would have multiple "predictor_time" values
    // below.
   }
#endif

   if(t_hist != 0)
   {
    std::string err = "Can't output history of predicted value: they aren't stored.";
    OomphLibWarning(err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
   }
   else if(Output_predictor_values && time_stepper_pt()->adaptive_flag())
   {
    const unsigned predictor_time =
     time_stepper_pt()->predictor_storage_index();

    std::ofstream pred_file((dir + "/" + "predsoln" + num + ".dat").c_str(),
			    std::ios::out);
    pred_file.precision(Output_precision);
    output_solution(predictor_time, pred_file, 2);
    pred_file.close();
   }


   Doc_info.number()++;
  }
 }

 double MyProblem::min_element_size()
 {
  // Check that we have finite elements ??ds this will still go wrong
  // if there are only some none finite elements in the mesh...

  //??ds what happens with face elements?
  FiniteElement* test_pt = dynamic_cast<FiniteElement*>
   (mesh_pt(0)->element_pt(0));

  if(test_pt != 0)
  {
   double min_size = mesh_pt(0)->finite_element_pt(0)->compute_physical_size();

   // Loop over all meshes in problem
   for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
   {
    Mesh* mesh_pt = this->mesh_pt(msh);
    for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
    {
     FiniteElement* ele_pt = mesh_pt->finite_element_pt(ele);
     double new_size = ele_pt->compute_physical_size();
     if(new_size < min_size)
     {
      min_size = new_size;
     }
    }
   }

   return min_size;
  }
  // If it's not a finite element then we can't get a size so return
  // a dummy value.
  else
  {
   return Dummy_doc_data;
  }
 }

 void MyProblem::write_trace(const unsigned& t_hist)
 {
  std::ofstream trace_file((Doc_info.directory() + "/" + Trace_filename).c_str(),
			   std::ios::app);
  trace_file.precision(Output_precision);

  double time = Dummy_doc_data, dt = Dummy_doc_data;
  if(!is_steady())
  {
   time = this->time_pt()->time(t_hist);
   dt = this->time_pt()->dt(t_hist);
  }

  // Some values are impossible to get for history values, print dummy
  // values instead.
  double nnewton_iter_taken = Dummy_doc_data,
   lte_norm = Dummy_doc_data,
   real_time = Dummy_doc_data,
   total_step_time = Dummy_doc_data;

  Vector<double> solver_iterations(1, Dummy_doc_data),
   solver_times(1, Dummy_doc_data),
   jacobian_setup_times(1, Dummy_doc_data),
   preconditioner_setup_times(1, Dummy_doc_data),
   max_res(1, Dummy_doc_data);

  if(t_hist == 0)
  {
   nnewton_iter_taken = Nnewton_iter_taken;
   solver_iterations = Solver_iterations;
   solver_times = Solver_times;
   jacobian_setup_times = Jacobian_setup_times;
   preconditioner_setup_times = Preconditioner_setup_times;
   real_time = std::time(0);
   total_step_time = Total_step_time;
   max_res = Max_res;

   if(!is_steady())
   {
    lte_norm = this->lte_norm();
   }
  }

  // Write out data that can be done for every problem
  trace_file
   << Doc_info.number()
   << Trace_seperator << time
   << Trace_seperator << dt
   << Trace_seperator << get_error_norm(t_hist)

   << Trace_seperator << nnewton_iter_taken;
     
  // Loop over the entries of solver_iterations and output
  unsigned solver_iterations_length=solver_iterations.size();
  if (solver_iterations_length==0)
  {
   trace_file << Trace_seperator << "[]";
  }
  else
  {
   trace_file << Trace_seperator << "[" << solver_iterations[0];
   if (solver_iterations_length>1)
   {
    for (unsigned i=1;i<solver_iterations_length;i++)
    {
     trace_file << ", " << solver_iterations[i];
    }
   }
   trace_file << "]";
  }
    
  // Loop over the entries of solver_times and output
  unsigned solver_times_length=solver_times.size();
  if (solver_times_length==0)
  {
   trace_file << Trace_seperator << "[]";
  }
  else
  {
   trace_file << Trace_seperator << "[" << solver_times[0];
   if (solver_times_length>1)
   {
    for (unsigned i=1;i<solver_times_length;i++)
    {
     trace_file << ", " << solver_times[i];
    }
   }
   trace_file << "]";
  }
    
  // Loop over the entries of jacobian_setup_times and output
  unsigned jacobian_setup_times_length=jacobian_setup_times.size();
  if (jacobian_setup_times_length==0)
  {
   trace_file << Trace_seperator << "[]";
  }
  else
  {
   trace_file << Trace_seperator << "[" << jacobian_setup_times[0];
   if (jacobian_setup_times_length>1)
   {
    for (unsigned i=1;i<jacobian_setup_times_length;i++)
    {
     trace_file << ", " << jacobian_setup_times[i];
    }
   }
   trace_file << "]";
  }
    
  // Loop over the entries of preconditioner_setup_times and output
  unsigned preconditioner_setup_times_length=preconditioner_setup_times.size();
  if (preconditioner_setup_times_length==0)
  {
   trace_file << Trace_seperator << "[]";
  }
  else
  {
   trace_file << Trace_seperator << "[" << preconditioner_setup_times[0];
   if (preconditioner_setup_times_length>1)
   {
    for (unsigned i=1;i<preconditioner_setup_times_length;i++)
    {
     trace_file << ", " << preconditioner_setup_times[i];
    }
   }
   trace_file << "]";
  }

  trace_file << Trace_seperator << lte_norm;

  // Create a temporary vector
  Vector<double> output_vector=trace_values(t_hist);

  // Loop over the entries of output_vector and output
  unsigned output_vector_length=output_vector.size();
  if (output_vector_length==0)
  {
   trace_file << Trace_seperator << "[]";
  }
  else
  {
   trace_file << Trace_seperator << "[" << output_vector[0];
   if (output_vector_length>1)
   {
    for (unsigned i=1;i<output_vector_length;i++)
    {
     trace_file << ", " << output_vector[i];
    }
   }
   trace_file << "]";
  }
  
  trace_file << Trace_seperator << real_time;

  // Loop over the entries of output_vector and output
  unsigned max_res_length=max_res.size();
  if (max_res_length==0)
  {
   trace_file << Trace_seperator << "[]";
  }
  else
  {
   trace_file << Trace_seperator << "[" << max_res[0];
   if (max_res_length>1)
   {
    for (unsigned i=1;i<max_res_length;i++)
    {
     trace_file << ", " << max_res[i];
    }
   }
   trace_file << "]";
  }
    
  trace_file << Trace_seperator << get_solution_norm(t_hist)
	     << Trace_seperator << total_step_time
     
   // Reserved slots in case I think of more things to add later
	     << Trace_seperator << Dummy_doc_data
	     << Trace_seperator << Dummy_doc_data
	     << Trace_seperator << Dummy_doc_data
	     << Trace_seperator << Dummy_doc_data
	     << Trace_seperator << Dummy_doc_data
	     << Trace_seperator << Dummy_doc_data;


  // Add problem specific data
  write_additional_trace_data(t_hist, trace_file);

  // Finish off this line
  trace_file << std::endl;
  trace_file.close();
 }

 void MyProblem::doc_boundaries(const std::string& boundary_file_basename) const
 {
  const unsigned n_boundary = mesh_pt()->nboundary();
  for(unsigned b=0; b<n_boundary; b++)
  {
   std::ofstream boundary_file((boundary_file_basename +
				to_string(b) + ".csv").c_str(),
			       std::ios::out);

   // write headers
   boundary_file << "x,y,z,b" << std::endl;

   const unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
   for(unsigned nd=0; nd<n_boundary_node; nd++)
   {
    Node* node_pt = mesh_pt()->boundary_node_pt(b, nd);
    Vector<double> position = node_pt->position();

    // Write out this node
    for(unsigned i=0; i<dim(); i++)
    {
     boundary_file << position[i] << ",";
    }
    for(unsigned i=dim(); i<3; i++)
    {
     boundary_file << 0.0 << ",";
    }
    boundary_file << b << std::endl;
   }

   boundary_file.close();
  }
 }

 void MyProblem::set_up_impulsive_initial_condition()
 {

#ifdef PARANOID
  if(nglobal_data() != 0)
  {
   std::string err = "Problem has global data which cannot be set from function pt.";
   throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
		       OOMPH_CURRENT_FUNCTION);
  }
#endif
  unsigned nprev_steps=this->time_stepper_pt()->nprev_values();
  for(unsigned t=0; t< nprev_steps; t++)
  {
   // Loop over all nodes in all meshes in problem and set values.
   for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
   {
    Mesh* mesh_pt = this->mesh_pt(msh);

    for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
    {
     Node* nd_pt = mesh_pt->node_pt(nd);
     for(unsigned j=0, nj=nd_pt->nvalue(); j<nj; j++)
     {
      nd_pt->set_value(t, j, nd_pt->value(0, j));
     }
    }

#ifdef PARANOID
    for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
    {
     FiniteElement* ele_pt = mesh_pt->finite_element_pt(ele);
     if(ele_pt->ninternal_data() != 0)
     {
      std::string err =
       "Element with non-nodal data, cannot set via function...";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
			  OOMPH_CURRENT_FUNCTION);
     }
    }
#endif

   }
  }

  actions_after_set_initial_condition();
 }

 void MyProblem::my_set_initial_condition(const SolutionFunctorBase& ic)
 {
#ifdef PARANOID
  // Can't set global data from a function of space... have to overload this
  // if you have any
  if(nglobal_data() != 0)
  {
   std::string err = "Problem has global data which cannot be set from function pt.";
   throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
		       OOMPH_CURRENT_FUNCTION);
  }
#endif

  // Loop over current & previous timesteps (so we need t<nprev_values+1).
  const int nprev_values = this->time_stepper_pt()->nprev_values();
  for(int tindex=0; tindex<nprev_values+1; tindex++)
  {
   double time = time_pt()->time(tindex);

   // Loop over all nodes in all meshes in problem and set values.
   const unsigned nmsh = nsub_mesh();
   for(unsigned msh=0; msh<nmsh; msh++)
   {
    Mesh* mesh_pt = this->mesh_pt(msh);

    for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
    {
     Node* nd_pt = mesh_pt->node_pt(nd);

     // Get the position at present time
     const unsigned dim = nd_pt->ndim();
     Vector<double> x(dim);
     nd_pt->position(x);

     // Set position at tindex time to be the same as at present
     // (impulsive positions).
     for(unsigned j=0; j<dim; j++)
     {
      nd_pt->x(tindex, j) = x[j];
     }

     // Get the values
     Vector<double> values = ic(time, x);

#ifdef PARANOID
     if(values.size() != nd_pt->nvalue())
     {
      std::string err = "Wrong number of values in initial condition.";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
			  OOMPH_CURRENT_FUNCTION);
     }
#endif
     // Copy into dofs
     for(unsigned j=0, nj=values.size(); j<nj; j++)
     {
      nd_pt->set_value(tindex, j, values[j]);
     }
    }

#ifdef PARANOID
    // Can't set internal data like this so check that we have none.
    for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
    {
     FiniteElement* ele_pt = mesh_pt->finite_element_pt(ele);
     if(ele_pt->ninternal_data() != 0)
     {
      std::string err =
       "Element with non-nodal data, cannot set via function...";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
			  OOMPH_CURRENT_FUNCTION);
     }
    }
#endif
   }
  }

  actions_after_set_initial_condition();
 }

} // End of oomph namespace

#endif
