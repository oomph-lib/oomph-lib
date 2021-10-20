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
#ifndef OOMPH_SEGREGATED_FSI_SOLVER
#define OOMPH_SEGREGATED_FSI_SOLVER


#include "../generic/problem.h"
#include "../generic/geom_objects.h"
#include "../generic/mesh.h"

namespace oomph
{
  //===============================================================
  /// Object that collates convergence data of Picard iteration
  //===============================================================
  class PicardConvergenceData
  {
  public:
    /// Constructor initialises all data
    PicardConvergenceData()
      : Niter(0),
        CPU_total(0.0),
        Essential_cpu_total(0.0),
        CPU_for_global_residual(0.0),
        Tol_achieved(0.0),
        Has_converged(false)
    {
    }

    /// Empty destructor
    ~PicardConvergenceData(){};

    /// Number of iterations performed
    unsigned& niter()
    {
      return Niter;
    }

    /// Total CPU time for segregated solve
    double& cpu_total()
    {
      return CPU_total;
    }

    ///  Total essential CPU time for segregated solve
    /// (excluding any actions that merely doc the progress
    /// of the iteration, etc.)
    double& essential_cpu_total()
    {
      return Essential_cpu_total;
    }

    ///  CPU time for computation of global residual vectors
    /// Note: This time is contained in Total_CPU and is
    /// only used if convergence is based on the residual
    /// of the fully coupled system.
    double& cpu_for_global_residual()
    {
      return CPU_for_global_residual;
    }

    /// Final tolerance achieved by the iteration
    double& tol_achieved()
    {
      return Tol_achieved;
    }

    /// Flag to indicate if the solver has converged
    bool has_converged() const
    {
      return Has_converged;
    }

    /// Set the flag to indicate that the solver has converged
    void set_solver_converged()
    {
      Has_converged = true;
    }

    /// Set the flag to indicate that the solver has not converged
    void set_solver_not_converged()
    {
      Has_converged = false;
    }

  private:
    /// Number of iterations performed
    unsigned Niter;

    /// Total CPU time for segregated solve
    double CPU_total;

    ///  Total essential CPU time for segregated solve
    /// (excluding any actions that merely doc the progress
    /// of the iteration, etc.)
    double Essential_cpu_total;

    ///   CPU time for computation of global residual vectors
    /// Note: This time is contained in Total_CPU and is
    /// only used if convergence is based on the residual
    /// of the fully coupled system
    double CPU_for_global_residual;

    /// Final tolerance achieved by the iteration
    double Tol_achieved;

    /// Flag to indicate if the solver has converged
    bool Has_converged;
  };


  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// A class to handle errors in the Segregated solver
  //=======================================================================
  class SegregatedSolverError
  {
  public:
    /// Default constructor, does nothing
    SegregatedSolverError(const bool& ran_out_of_iterations = false)
    {
      Ran_out_of_iterations = ran_out_of_iterations;
    }

    /// Boolean indicating if the error occured because
    /// we ran out of iterations
    bool Ran_out_of_iterations;
  };


  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////


  //===============================================================
  /// Base class for problems that can be solved by segregated
  /// FSI solver
  //===============================================================
  class SegregatableFSIProblem : public virtual Problem
  {
  protected:
    ///  This function is called once at the start of each
    /// segregated solve.
    virtual void actions_before_segregated_solve() {}

    ///  This function is called once at the end of each
    /// segregated solve.
    virtual void actions_after_segregated_solve() {}

    ///  This function is to be filled with actions that take place
    /// before the check for convergence of the entire segregated solve
    virtual void actions_before_segregated_convergence_check() {}

  public:
    ///  Constructor. Set default values for solver parameters:
    /// - Don't use pointwise Aitken extrapolation but if it's used at
    ///   all, it's used immediately.
    /// - No under-relaxation at all (neither classical nor Irons&Tuck)
    /// - Base convergence on max. residual of coupled system of eqns
    /// - Convergence tolerance = tolerance for Newton solver
    ///   defined in Problem base class
    /// - Max. 50 Picard iterations
    SegregatableFSIProblem()
    {
      // Use pointwise Aitken extrapolation?
      Use_pointwise_aitken = false;

      // Default: No under-relaxation
      Omega_relax = 1.0;

      // Don't use of Irons and Tuck's extrapolation for solid values
      Use_irons_and_tuck_extrapolation = false;

      // Start using pointwise Aitken immediately
      Pointwise_aitken_start = 0;

      // By default we don't recheck convergence
      Recheck_convergence_after_pointwise_aitken = false;

      // Default solve type is full solve
      Solve_type = Full_solve;

      // Convergence criterion
      Convergence_criterion = Assess_convergence_based_on_max_global_residual;

      // Convergence tolerance (as in global Newton solver)
      Convergence_tolerance = Problem::Newton_solver_tolerance;

      // Doc max. global residual during iteration?
      Doc_max_global_residual = false;

      // Max. number of Picard iterations
      Max_picard = 50;

      // Pointer to Mesh containing only fluid elements -- the elements in this
      // Mesh will be excluded from the assembly process when
      // the solid problem is solved
      Fluid_mesh_pt = 0;

      // Pointer to Mesh containing only solid elements -- the elements in this
      // mesh will be excluded from the assembly process when
      // the fluid problem is solved
      Solid_mesh_pt = 0;

      // Initialise timer that allows doc of iteration/cpu time
      T_ref = clock();
      T_spent_on_actual_solve = 0.0;

      ///  boolean flag to indicate if timer has been halted
      Timer_has_been_halted = false;
    }

    /// Empty virtual destructor
    virtual ~SegregatableFSIProblem() {}

    ///  Identify the fluid and solid Data. This is a pure virtual
    /// function that MUST be implemented for every specific problem that
    /// is to be solved by the segregated solver.
    /// The two mesh pointers identify meshes that contain
    /// elements and nodes used during the fluid
    /// or solid solves respectively. Elements that feature during
    /// both phases of the segretated solution must be included in both.
    /// These pointers may be set to NULL. In this case, all elements
    /// in the global mesh (set up during the monolithic discretisation
    /// of the problem) contribute to the global Jacobian matrix
    /// and the residual vector, even if their contributions only contain
    /// zero entries. This can be costly, though the code will
    /// still compute the correct results.
    virtual void identify_fluid_and_solid_dofs(Vector<Data*>& fluid_data_pt,
                                               Vector<Data*>& solid_data_pt,
                                               Mesh*& fluid_mesh_pt,
                                               Mesh*& solid_mesh_pt) = 0;

    ///  Setup the segregated solver: Backup the pinned status of
    /// the fluid and solid dofs and allocate the internal storage
    /// based on the input provided by identify_fluid_and_solid_dofs(...)
    /// In addition, reset storage associated with convergence acceleration
    /// techniques.
    /// If the problem and degrees of freedom has not changed between
    /// calls to the solver then it is wasteful to call
    /// identify_fluid_and_solid_dofs(...) again and again. If the optional
    /// boolean flag is set to false then the storage for convergence
    /// acceleration techniques is reset, but the fluid and solid dofs
    /// are not altered.
    void setup_segregated_solver(
      const bool& full_setup_of_fluid_and_solid_dofs = true);

    ///  Segregated solver. Peform a segregated step from
    /// the present state of the system.
    /// Returns PicardConvergenceData object that contains the vital
    /// stats of the iteration
    PicardConvergenceData segregated_solve();

    ///  Steady version of segregated solver. Makes all
    /// timesteppers steady before solving.
    /// Returns PicardConvergenceData object that contains the
    /// vital stats of the iteration.
    PicardConvergenceData steady_segregated_solve();


    ///  Unsteady segregated solver, advance time by dt and solve
    /// by the segregated solver. The time values are always shifted by
    /// this function.
    /// Returns PicardConvergenceData object that contains the
    /// vital stats of the iteration.
    PicardConvergenceData unsteady_segregated_solve(const double& dt);


    ///  Unsteady segregated solver. Advance time by dt and solve
    /// the system by a segregated method. The boolean flag is used to
    /// control whether the time values should be shifted. If it is true
    /// the current data values will be shifted (stored as previous
    /// timesteps) before the solution step.
    /// Returns PicardConvergenceData object that contains the
    /// vital stats of the iteration.
    PicardConvergenceData unsteady_segregated_solve(const double& dt,
                                                    const bool& shift_values);


    ///  Assess convergence based on max. residual of coupled system of
    /// eqns. The argument specifies the convergence tolerance.
    void assess_convergence_based_on_max_global_residual(const double& tol)
    {
      Convergence_criterion = Assess_convergence_based_on_max_global_residual;
      Convergence_tolerance = tol;
    }

    ///  Assess convergence based on max. residuals of coupled
    /// system of eqns. This interface has no argument
    /// and the default convergence tolerance
    /// for the Newton solver, Problem::Newton_solver_tolerance is used.
    void assess_convergence_based_on_max_global_residual()
    {
      assess_convergence_based_on_max_global_residual(
        Problem::Newton_solver_tolerance);
    }

    ///  Assess convergence based on max. absolute change of solid
    /// dofs. The argument specifies the convergence tolerance.
    void assess_convergence_based_on_absolute_solid_change(const double& tol)
    {
      Convergence_criterion = Assess_convergence_based_on_absolute_solid_change;
      Convergence_tolerance = tol;
    }

    ///  Assess convergence based on max. absolute change of solid
    /// dofs. This interface has no argument and the default
    /// convergence tolerance
    /// for the Newton solver, Problem::Newton_solver_tolerance is used.
    void assess_convergence_based_on_absolute_solid_change()
    {
      assess_convergence_based_on_absolute_solid_change(
        Problem::Newton_solver_tolerance);
    }

    ///  Assess convergence based on max. relative change of solid
    /// dofs. The argument specifies the convergence tolerance.
    void assess_convergence_based_on_relative_solid_change(const double& tol)
    {
      Convergence_criterion = Assess_convergence_based_on_relative_solid_change;
      Convergence_tolerance = tol;
    }

    ///  Assess convergence based on max. relative change of solid
    /// dofs. This interface has no argument and the default
    /// convergence tolerance
    /// for the Newton solver, Problem::Newton_solver_tolerance is used.
    void assess_convergence_based_on_relative_solid_change()
    {
      assess_convergence_based_on_relative_solid_change(
        Problem::Newton_solver_tolerance);
    }


    ///  Use pointwise Aitken extrapolation. The argument is used to
    /// specify the Picard iteration after which pointwise Aitken extrapolation
    /// is to be used for the first time.
    void enable_pointwise_aitken(const unsigned& pointwise_aitken_start)
    {
      Pointwise_aitken_start = pointwise_aitken_start;
      Use_pointwise_aitken = true;
    }

    ///  Use pointwise Aitken extrapolation. This interface has
    /// no argument and the current value of Pointwise_aitken_start will
    /// be used. The default is zero, extrapolation starts immediately
    void enable_pointwise_aitken()
    {
      Use_pointwise_aitken = true;
    }

    ///  Disable the use of pointwise Aitken extrapolation
    void disable_pointwise_aitken()
    {
      Use_pointwise_aitken = false;
    }

    /// Use under-relaxation and (optionally) specify under-relaxation
    /// parameter. Default: omega=1.0 (i.e. no actual under-relaxation;
    /// Other extreme: omega=0.0 (freeze wall shape). Under-relaxation
    /// parameter can also be computed dynamically by setting
    /// use_irons_and_tuck_extrapolation()
    void enable_under_relaxation(const double& omega = 1.0)
    {
      Omega_relax = omega;
    }

    /// Use Irons and Tuck extrapolation for solid dofs
    void enable_irons_and_tuck_extrapolation()
    {
      Use_irons_and_tuck_extrapolation = true;
    }

    /// Do not use Irons and Tuck extrapolation for solid dofs
    void disable_irons_and_tuck_extrapolation()
    {
      Use_irons_and_tuck_extrapolation = false;
    }

    /// Enumerated flags for convergence criteria
    enum convergence_criteria
    {
      Assess_convergence_based_on_absolute_solid_change,
      Assess_convergence_based_on_relative_solid_change,
      Assess_convergence_based_on_max_global_residual
    };


    /// Enumerated flags to indicate which solve is taking place
    enum solve_type
    {
      Full_solve,
      Fluid_solve,
      Solid_solve
    };

    ///  Get rms of change in the solid dofs; the max. change of the
    /// solid dofs and the rms norm of the solid dofs themselves.
    /// Change is computed relative to the reference values stored when
    /// store_solid_dofs() was last called.
    void get_solid_change(double& rms_change,
                          double& max_change,
                          double& rms_norm);

    ///  Store the current solid values as reference values for
    /// future convergence check. Also add another entry to pointwise
    /// Aitken history if required.
    void store_solid_dofs();

    ///  Reset timer
    void reset_timer()
    {
      T_spent_on_actual_solve = 0.0;
      T_ref = clock();
      Timer_has_been_halted = false;
    }


    ///  (Re-)start timer (e.g. after completing non-essential
    /// parts of the code such as documentation of the iteration's
    /// progress)
    void restart_timer()
    {
      T_ref = clock();
      Timer_has_been_halted = false;
    }


    ///  Halt timer (e.g. before performing non-essential
    /// parts of the code such as documentation of the iteration's
    /// progress)
    void halt_timer()
    {
      if (!Timer_has_been_halted)
      {
        T_spent_on_actual_solve += double(clock() - T_ref) / CLOCKS_PER_SEC;
        Timer_has_been_halted = true;
      }
    }


    ///  Total elapsed time since start of solve
    double t_spent_on_actual_solve()
    {
      halt_timer();
      double time = T_spent_on_actual_solve;
      restart_timer();
      return time;
    }


  protected:
    /// Rebuild global mesh for monolithic discretisation
    void rebuild_monolithic_mesh();

    ///  Number of Aitken histories available (int because after
    /// extrapolation it's re-initialised to -1 to force the computation
    /// of three new genuine iterates).
    int Pointwise_aitken_counter;

    /// Use pointwise Aitken extrapolation?
    bool Use_pointwise_aitken;

    ///  Start pointwise Aitken extrpolation after specified number
    /// of Picard iterations
    unsigned Pointwise_aitken_start;

    /// Solve that is taking place (enumerated flag)
    int Solve_type;

    /// Convergence tolerance for Picard iteration
    double Convergence_tolerance;

    /// Max. number of Picard iterations
    unsigned Max_picard;

    /// Doc maximum global residual during iteration? (default: false)
    bool Doc_max_global_residual;

    /// Restore pinned status of fluid dofs
    void restore_fluid_dofs();

    /// Pin solid dofs
    void pin_solid_dofs();

    /// Restore pinned status of solid dofs
    void restore_solid_dofs();


    /// Do pointwise Aitken extrapolation for solid
    void pointwise_aitken_extrapolate();

    ///  Vector storing the Data objects associated with the fluid
    /// problem: Tyically the nodal and internal data of the elements in the
    /// fluid bulk mesh
    Vector<Data*> Fluid_data_pt;

    ///  Vector of vectors that store the pinned status of
    /// fluid Data values
    Vector<std::vector<bool>> Fluid_value_is_pinned;

    ///  Vector storing the Data objects associated with the solid
    /// problem: Typically the positional data of solid nodes and
    /// any quantities associated with displacement control, say.
    Vector<Data*> Solid_data_pt;

    ///  Vector of vectors that store the pinned status of
    /// solid Data values
    Vector<std::vector<bool>> Solid_value_is_pinned;

    ///  Vector storing the previous solid values -- used for
    /// convergence check
    Vector<double> Previous_solid_value;

    ///  Mesh containing only fluid elements -- the elements in this
    /// Mesh will be excluded from the assembly process when
    /// the solid problem is solved
    Mesh* Fluid_mesh_pt;

    ///  Mesh containing only solid elements -- the elements in this
    /// mesh will be excluded from the assembly process when
    /// the fluid problem is solved
    Mesh* Solid_mesh_pt;

    ///  Backup for the pointers to the submeshes in the original problem
    Vector<Mesh*> Orig_sub_mesh_pt;

    /// Vector of changes in Irons and Tuck under-relaxation
    Vector<double> Del_irons_and_tuck;

    /// Irons and Tuck relaxation factor
    double R_irons_and_tuck;

    ///  Vector of Vectors containing up to three previous
    /// iterates for the solid dofs; used for pointwise Aitken extrapolation
    Vector<Vector<double>> Pointwise_aitken_solid_value;

    /// Have we just done a pointwise Aitken step
    bool Recheck_convergence_after_pointwise_aitken;

  private:
    /// Extrapolate solid data and update fluid mesh during unsteady run
    void extrapolate_solid_data();

    ///  Under-relax the most recently computed solid variables, either
    /// by classical relaxation or by Irons & Tuck
    void under_relax_solid();

    ///  Only include fluid elements in the Problem's mesh. This is
    /// called before the segregated fluid solve. The fluid elements are
    /// identified by the user via the fluid_mesh_pt argument
    /// in the pure virtual function identify_fluid_and_solid_dofs(...).
    /// If a NULL pointer is returned by this function (i.e. if the user
    /// hasn't bothered to identify the fluids elements in a submesh, then
    /// no stripping of non-fluid elements is performed. The result
    /// of the computation will be correct but
    /// it won't be very efficient.
    void use_only_fluid_elements();

    ///  Only include solid elements in the Problem's mesh. This is
    /// called before the segregated solid solve. The solid elements are
    /// identified by the user via the solid_mesh_pt argument
    /// in the pure virtual function identify_fluid_and_solid_dofs(...).
    /// If a NULL pointer is returned by this function (i.e. if the user
    /// hasn't bothered to identify the solid elements in a submesh, then
    /// no stripping of non-solid elements is performed. The result
    /// of the computation will be correct but
    /// it won't be very efficient.
    void use_only_solid_elements();

    /// Pin fluid dofs
    void pin_fluid_dofs();

    ///  Under-relaxation parameter. (1.0: no under-relaxation;
    /// 0.0: Freeze wall shape)
    double Omega_relax;

    ///  Boolean flag to indicate use of Irons and Tuck's extrapolation
    /// for solid values
    bool Use_irons_and_tuck_extrapolation;

    /// Convergence criterion (enumerated flag)
    int Convergence_criterion;

    ///  Reference time for segregated solve. Can be
    /// re-initialised whenever total elapsed time has been stored
    /// (before entering non-essential doc sections of the code)
    clock_t T_ref;

    ///  Total elapsed time since start of solve, can be
    /// accumulated by adding bits of time spent in relevant parts of
    /// code (bypassing sections that only document the progress)
    double T_spent_on_actual_solve;

    ///  boolean flag to indicate if timer has been halted
    bool Timer_has_been_halted;
  };

} // namespace oomph


#endif
