// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_ELASTIC_PROBLEMS_HEADER
#define OOMPH_ELASTIC_PROBLEMS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// oomph-lib headers
#include "geom_objects.h"
#include "timesteppers.h"
#include "problem.h"
#include "frontal_solver.h"
#include "mesh.h"
#include "mumps_solver.h"

namespace oomph
{
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Dummy mesh that can be created and deleted in
  /// SolidICProblem
  //======================================================================
  class DummyMesh : public Mesh
  {
  public:
    // Empty constructor
    DummyMesh(){};
  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// IC problem for an elastic body discretised on a given (sub)-mesh.
  /// We switch the elements' residuals and Jacobians to the system of
  /// equations that forces the wall shape to become that of
  /// a specified "initial condition object".
  ///
  ///
  /// Boundary conditions for all data items associated with the mesh's
  /// elements' are temporarily overwritten so that all positional dofs (and
  /// only those!) become free while all other data (nodal data and the
  /// element's internal and external data) are either pinned (for nodal and
  /// internal data) or completely removed (for pointers to external data). The
  /// complete removal of the (pointers to the) external data is necessary
  /// because in FSI problems, positional data of certain elements
  /// can feature as external data of others. Hence pinning them
  /// in their incarnation as external data would also pin them
  /// in their incarnation as positional data.
  ///
  ///
  /// All data and its pinned-status are restored at the end of the
  /// call to setup_ic().
  //======================================================================
  class SolidICProblem : public Problem
  {
   
  public:

   /// Constructor. Initialise pointer
    /// to IC object to NULL. Create a dummy mesh that can be deleted
    /// when static problem finally goes out of scope at end of
    /// program execution.
   SolidICProblem() : IC_pt(0)
    {

#ifdef OOMPH_HAS_MPI
     Mumps_solver_pt=new MumpsSolver;
#endif
     SuperLU_solver_pt=new SuperLUSolver;
     
     // Create dummy mesh
     mesh_pt() = new DummyMesh;
     
     // Default value for checking of consistent assignment of Newmark IC
     Max_residual_after_consistent_newton_ic = 1.0e-12;
    }
   
   /// Destructor
   ~SolidICProblem() 
    {     
#ifdef OOMPH_HAS_MPI
     delete Mumps_solver_pt;
#endif
     delete SuperLU_solver_pt;
    }

   
   /// Broken copy constructor
   SolidICProblem(const SolidICProblem&) = delete;

    /// Broken assignment operator
    void operator=(const SolidICProblem&) = delete;

    /// Update after solve (empty)
    void actions_after_newton_solve() {}

    /// Update the problem specs before solve.  (empty)
    void actions_before_newton_solve() {}

    /// Force the elastic structure that is discretised on the specified
    /// mesh to deform in the shape of the initial condition object
    /// (evaluated at the time specified)
    void set_static_initial_condition(Problem* problem_pt,
                                      Mesh* mesh_pt,
                                      SolidInitialCondition* ic_pt,
                                      const double& time);

    /// Force the elastic structure that is discretised on the specified
    /// mesh to deform in the shape of the initial condition object (wrapper)
    void set_static_initial_condition(Problem* problem_pt,
                                      Mesh* mesh_pt,
                                      SolidInitialCondition* ic_pt)
    {
      double time;
      set_static_initial_condition(problem_pt, mesh_pt, ic_pt, time);
    }

    /// Setup initial condition for time-integration
    /// with Newmark's method. History values are assigned to that
    /// the velocity and accelerations determined by the Newmark
    /// scheme are exact at the initial time.
    template<class TIMESTEPPER>
    void set_newmark_initial_condition_directly(Problem* problem_pt,
                                                Mesh* mesh_pt,
                                                TIMESTEPPER* timestepper_pt,
                                                SolidInitialCondition* ic_pt,
                                                const double& dt);


    /// Setup initial condition for time-integration
    /// with Newmark's method. Past displacements and velocities are assigned
    /// directly (consistent with the idea that a second order ODE
    /// can take ICs up to 1st order, while the history value for
    /// the previous acceleration is determined by the condition that
    /// the weak equation is satisfied at the initial time.)
    /// The multiplier function needs to specify the factor that
    /// multiplies the inertia terms -- typically this is a
    /// constant, given by the ratio \f$ \Lambda^2 \f$ of the
    /// problem's intrinsic timescale to the time used to non-dimensionalise
    /// the equations. If the function (pointer) is not specified
    /// the multiplier is assumed to be equal to 1.0 -- appropriate
    /// for a non-dimensionalisation based on the problem's intrinsic timescale.
    template<class TIMESTEPPER>
    void set_newmark_initial_condition_consistently(
      Problem* problem_pt,
      Mesh* mesh_pt,
      TIMESTEPPER* timestepper_pt,
      SolidInitialCondition* ic_pt,
      const double& dt,
      SolidFiniteElement::MultiplierFctPt multiplier_fct_pt = 0);


    /// Max. tolerated residual after application of consistent
    /// Newmark IC. Used to check if we have specified the correct
    /// timescale ratio (non-dim density).
    double& max_residual_after_consistent_newton_ic()
    {
      return Max_residual_after_consistent_newton_ic;
    }


  private:

   
    /// Backup original state of all data associated with mesh
    void backup_original_state();

    /// Reset original state of all data associated with mesh
    void reset_original_state();

    /// Change pinned status of all data associated with mesh
    /// so that the IC problem can be solved.
    void setup_problem();

    /// Pointer to initial condition object
    SolidInitialCondition* IC_pt;

    /// Vector to store pinned status of all data
    Vector<int> Backup_pinned;

    /// Vector of Vectors  to store pointers to exernal data in the elements
    Vector<Vector<Data*>> Backup_ext_data;

    /// Max. tolerated residual after application of consistent
    /// Newmark IC. Used to check if we have specified the correct
    /// timescale ratio (non-dim density).
    double Max_residual_after_consistent_newton_ic;


#ifdef OOMPH_HAS_MPI
   /// Pointer to mumps solver
   MumpsSolver* Mumps_solver_pt;
#endif
   
   /// Pointer to mumps solver
   SuperLUSolver* SuperLU_solver_pt;

  };


  //======================================================================
  /// Setup initial condition for time-integration
  /// with Newmark's method. History values are assigned to that
  /// the velocity and accelerations determined by the Newmark
  /// scheme are exact at the initial time.
  //======================================================================
  template<class TIMESTEPPER>
  void SolidICProblem::set_newmark_initial_condition_directly(
    Problem* problem_pt,
    Mesh* wall_mesh_pt,
    TIMESTEPPER* timestepper_pt,
    SolidInitialCondition* ic_pt,
    const double& dt)
  {
#ifdef PARANOID
    if (timestepper_pt->type() != "Newmark")
    {
      std::ostringstream error_message;
      error_message
        << "SolidICProblem::set_newmark_initial_condition_directly()\n"
        << "can only be called for Newmark type timestepper whereas\n "
        << "you've called it for " << timestepper_pt->type() << std::endl;

      throw OomphLibError(
        error_message.str(),
        "SolidICProblem::set_newmark_initial_condition_directly()",
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Set value of dt
    timestepper_pt->time_pt()->dt() = dt;

    // Set the weights
    timestepper_pt->set_weights();

    // Delete dummy mesh
    delete mesh_pt();

    // Set pointer to mesh
    mesh_pt() = wall_mesh_pt;

    // Set pointer to initial condition object
    IC_pt = ic_pt;

    // Backup the pinned status of all dofs and remove external data
    // of all elements
    backup_original_state();

    // Now alter the pinned status so that the IC problem for the
    // positional variables can be solved; setup equation numbering
    // scheme
    setup_problem();


    // Choose the right linear solver
#ifdef OOMPH_HAS_MPI
    if (MPI_Helpers::mpi_has_been_initialised())
     {
      linear_solver_pt()=Mumps_solver_pt;
     }
    else
     {
      linear_solver_pt()=SuperLU_solver_pt;
     }
#else
    linear_solver_pt()=SuperLU_solver_pt;
#endif
    
    // Store times at which we need to assign ic:
    double current_time = timestepper_pt->time_pt()->time();
    double previous_time = timestepper_pt->time_pt()->time(1);

    // Stage 1: Set values and time derivs at current time
    //----------------------------------------------------

    // [Note: this acts on time everywhere!]
    IC_pt->geom_object_pt()->time_stepper_pt()->time_pt()->time() =
      current_time;

    // Loop over time-derivatives
    for (unsigned t_deriv = 0; t_deriv <= 2; t_deriv++)
    {
      // Set flag to ensure that the t_deriv-th time derivative
      // of the prescribed solution gets stored in displacements
      IC_pt->ic_time_deriv() = t_deriv;

      // Solve the problem for initial shape
      newton_solve();

      // Loop over all the nodes
      unsigned n_node = mesh_pt()->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        // Assign current time derivative in its temporary storage
        // position
        timestepper_pt->assign_initial_data_values_stage1(
          t_deriv,
          dynamic_cast<SolidNode*>(mesh_pt()->node_pt(n))
            ->variable_position_pt());
      }
    }


    // Stage 2: Now get position at previous time and adjust previous
    //---------------------------------------------------------------
    // values of veloc and accel in Newmark scheme so that current
    //---------------------------------------------------------------
    // veloc and accel (specified in step 1) are represented exactly.
    //---------------------------------------------------------------

    // [Note: this acts on time everywhere and needs to be reset!]
    IC_pt->geom_object_pt()->time_stepper_pt()->time_pt()->time() =
      previous_time;

    // Set flag to ensure that the t_deriv-th time derivative
    // of the prescribed solution gets stored in displacements
    IC_pt->ic_time_deriv() = 0;

    // Solve the problem for initial shape
    newton_solve();

    // Loop over all the nodes and make the final adjustments
    unsigned n_node = mesh_pt()->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      timestepper_pt->assign_initial_data_values_stage2(
        dynamic_cast<SolidNode*>(mesh_pt()->node_pt(n))
          ->variable_position_pt());
    }

    // Reset time
    IC_pt->geom_object_pt()->time_stepper_pt()->time_pt()->time() =
      current_time;

    // Reset the pinned status and re-attach the external data to the elements
    reset_original_state();

    // Set pointer to dummy mesh so there's something that can be deleted
    // when static problem finally goes out of scope.
    mesh_pt() = new DummyMesh;

    // We have temporarily over-written equation numbers -- need
    // to reset them now
    oomph_info << "Number of equations in big problem: "
               << problem_pt->assign_eqn_numbers() << std::endl;
  }


  //======================================================================
  /// Setup initial condition for time-integration
  /// with Newmark's method. Past displacements and velocities are assigned
  /// directly (consistent with the idea that a second order ODE
  /// can take ICs up to 1st order, while the history value for
  /// the previous acceleration is determined by the condition that
  /// that the weak equation is satisfied at the initial time.
  /// The multiplier function needs to specify the factor that
  /// multiplies the inertia terms -- typically this is a
  /// constant, given by the ratio \f$ \Lambda^2 \f$ of the
  /// problem's intrinsic timescale to the time used to non-dimensionalise
  /// the equations. If the function (pointer) is not specified,
  /// the multiplier is assumed to be equal to 1.0 -- appropriate
  /// for a non-dimensionalisation based on the problem's intrinsic timescale.
  //======================================================================
  template<class TIMESTEPPER>
  void SolidICProblem::set_newmark_initial_condition_consistently(
    Problem* problem_pt,
    Mesh* wall_mesh_pt,
    TIMESTEPPER* timestepper_pt,
    SolidInitialCondition* ic_pt,
    const double& dt,
    SolidFiniteElement::MultiplierFctPt multiplier_fct_pt)
  {
#ifdef PARANOID
    if (timestepper_pt->type() != "Newmark")
    {
      std::ostringstream error_message;
      error_message
        << "SolidICProblem::set_newmark_initial_condition_consistently()\n"
        << "can only be called for Newmark type timestepper whereas\n "
        << "you've called it for " << timestepper_pt->type() << std::endl;

      throw OomphLibError(
        error_message.str(),
        "SolidICProblem::set_newmark_initial_condition_consistently()",
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Set value of dt
    timestepper_pt->time_pt()->dt() = dt;

    // Set the weights
    timestepper_pt->set_weights();

    // Delete dummy mesh
    delete mesh_pt();

    // Set pointer to mesh
    mesh_pt() = wall_mesh_pt;

    // Set pointer to initial condition object
    IC_pt = ic_pt;

    // Backup the pinned status of all dofs and remove external data
    // of all elements
    backup_original_state();

    // Now alter the pinned status so that the IC problem for the
    // positional variables can be solved; setup equation numbering
    // scheme
    setup_problem();

    // Choose the right linear solver
#ifdef OOMPH_HAS_MPI
    if (MPI_Helpers::mpi_has_been_initialised())
     {
      linear_solver_pt()=Mumps_solver_pt;
     }
    else
     {
      linear_solver_pt()=SuperLU_solver_pt;
     }
#else
    linear_solver_pt()=SuperLU_solver_pt;
#endif

    
    // Number of history values
    unsigned ntstorage =
      IC_pt->geom_object_pt()->time_stepper_pt()->ntstorage();

    // Set values at previous time
    //----------------------------

    // Loop over number of previous times stored
    unsigned nprevtime =
      IC_pt->geom_object_pt()->time_stepper_pt()->nprev_values();

    // Backup previous times:
    Vector<double> prev_time(nprevtime + 1);
    for (unsigned i = 0; i <= nprevtime; i++)
    {
      prev_time[i] =
        IC_pt->geom_object_pt()->time_stepper_pt()->time_pt()->time(i);
    }

    // Loop over previous times & set values themselves
    //-------------------------------------------------
    for (unsigned i = 1; i <= nprevtime; i++)
    {
      // Set time for geometric object that specifies initial condition
      // [Note: this acts on time everywhere!]
      IC_pt->geom_object_pt()->time_stepper_pt()->time_pt()->time() =
        prev_time[i];

      // Set flag to ensure that the t_deriv-th time derivative
      // of the prescribed solution gets stored in displacements
      IC_pt->ic_time_deriv() = 0;

      // Solve the problem for initial shape: After this solve
      // The nodes's current positions represent the position at
      // previous time level i.
      newton_solve();

      // Loop over all the nodes
      unsigned n_node = mesh_pt()->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        // Get the variable position data
        Data* position_data_pt = dynamic_cast<SolidNode*>(mesh_pt()->node_pt(n))
                                   ->variable_position_pt();
        // Get number of values
        unsigned nval = position_data_pt->nvalue();

        // Assign values at previous times into their corresponding
        // slots
        for (unsigned ival = 0; ival < nval; ival++)
        {
          position_data_pt->set_value(
            i, ival, position_data_pt->value(0, ival));
        }
      }
    }

    // Set veloc (1st time deriv) at previous time and store in appropriate
    //---------------------------------------------------------------------
    // history value. Also assign zero value for 2nd time deriv. at
    //-------------------------------------------------------------
    // previous time
    //--------------

    // Set time for geometric object that specifies initial condition
    // to previous time [Note: this acts on time everywhere!]
    IC_pt->geom_object_pt()->time_stepper_pt()->time_pt()->time() =
      prev_time[1];

    // Set flag to ensure that the t_deriv-th time derivative
    // of the prescribed solution gets stored in displacements
    IC_pt->ic_time_deriv() = 1;

    // Solve the problem for initial shape: After this solve
    // The nodes's current positions represent the time derivatives of at
    // the positons at previous time level.
    newton_solve();

    // Loop over all the nodes
    unsigned n_node = mesh_pt()->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      // Get the position data
      Data* position_data_pt =
        dynamic_cast<SolidNode*>(mesh_pt()->node_pt(n))->variable_position_pt();

      // Get number of values
      unsigned nval = position_data_pt->nvalue();

      // Assign values at previous times into their corresponding
      // slots (last but one history value of Newmark scheme)
      for (unsigned ival = 0; ival < nval; ival++)
      {
        position_data_pt->set_value(
          ntstorage - 2, ival, position_data_pt->value(0, ival));
        position_data_pt->set_value(ntstorage - 1, ival, 0.0);
      }
    }


    // Set values at current time
    //---------------------------

    // Reset time to current value
    // [Note: this acts on time everywhere!]
    IC_pt->geom_object_pt()->time_stepper_pt()->time_pt()->time() =
      prev_time[0];

    // Set flag to ensure that the t_deriv-th time derivative
    // of the prescribed solution gets stored in displacements
    IC_pt->ic_time_deriv() = 0;

    // Solve the problem for initial shape
    newton_solve();


    // Now solve for the correction to the Newmark accelerations
    //----------------------------------------------------------
    // at previous time:
    //------------------

    // Loop over the elements
    unsigned Nelement = mesh_pt()->nelement();
    for (unsigned i = 0; i < Nelement; i++)
    {
      // Cast to proper element type
      SolidFiniteElement* elem_pt =
        dynamic_cast<SolidFiniteElement*>(mesh_pt()->element_pt(i));

      // Switch system to the one that determines the Newmark accelerations
      // by setting the Jacobian to the mass matrix
      elem_pt->enable_solve_for_consistent_newmark_accel();

      // Set pointer to multiplier function
      elem_pt->multiplier_fct_pt() = multiplier_fct_pt;

      // Switch off pointer to initial condition object
      elem_pt->solid_ic_pt() = 0;
    }

    // Correction vector
    DoubleVector correction;

    /// Pointer to member type solver
    typedef void (LinearSolver::*SolverMemPtr)(Problem* const& problem,
                                               DoubleVector& result);
    SolverMemPtr solver_mem_pt = &LinearSolver::solve;

    // Now do the linear solve
    LinearSolver* lin_solver_pt=0;
#ifdef OOMPH_HAS_MPI
    if (MPI_Helpers::mpi_has_been_initialised())
     {
      lin_solver_pt=Mumps_solver_pt;
     }
    else
     {
      lin_solver_pt=SuperLU_solver_pt;
     }
#else
    lin_solver_pt=SuperLU_solver_pt;
#endif

    
    (lin_solver_pt->*solver_mem_pt)(this, correction);

    // Update discrete 2nd deriv at previous time so that it's consistent
    // with PDE at current time

    // Loop over all the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // Get the pointer to the position data
      Data* position_data_pt =
        dynamic_cast<SolidNode*>(mesh_pt()->node_pt(n))->variable_position_pt();

      // Get number of values
      unsigned nval = position_data_pt->nvalue();

      // Assign values for the history value that corresponds to the
      // previous accel in Newmark scheme so that the PDE is satsified
      // at current time
      for (unsigned ival = 0; ival < nval; ival++)
      {
        // Global equation number
        int ieqn = position_data_pt->eqn_number(ival);

#ifdef PARANOID
        if (ieqn < 0)
        {
          throw OomphLibError(
            "No positional dofs should be pinned at this stage!",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Update the value
        *(position_data_pt->value_pt(ntstorage - 1, ival)) -= correction[ieqn];
      }
    }


#ifdef PARANOID
    // Check the residual
    DoubleVector residuals;
    get_residuals(residuals);
    double res_max = residuals.max();
    oomph_info
      << "Max. residual after assigning consistent initial conditions: "
      << res_max << std::endl;
    if (res_max > Max_residual_after_consistent_newton_ic)
    {
      std::ostringstream error_message;
      error_message << "Residual is  bigger than allowed! [Current tolerance: "
                    << Max_residual_after_consistent_newton_ic << "]\n\n";
      error_message << "This is probably because you've not specified the "
                    << "correct multiplier \n(the product of growth factor "
                    << "and timescale ratio [the non-dim density]). \nPlease "
                    << "check the Solid Mechanics Theory Tutorial for "
                    << "details. \n\n"
                    << "If you're sure that the residual is OK, overwrite "
                    << "the default tolerance using\n";
      error_message
        << "SolidICProblem::max_residual_after_consistent_newton_ic()"
        << std::endl
        << "or recompile without the PARANOID flag." << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Reset the pinned status and re-attach the external data to the elements
    reset_original_state();

    // Set pointer to dummy mesh so there's something that can be deleted
    // when static problem finally goes out of scope.
    mesh_pt() = new DummyMesh;

    // We have temporarily over-written equation numbers -- need
    // to reset them now
    oomph_info << "Number of equations in big problem: "
               << problem_pt->assign_eqn_numbers() << std::endl;
  }

} // namespace oomph

#endif
