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
// Generic oomph-lib stuff
#include "generic.h"
#include "navier_stokes.h"

// Oscillating cylinder mesh
#include "meshes/rectangle_with_moving_cylinder_mesh.h"

// General libraries
#include <limits>

using namespace oomph;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// The specialisation of the PMLLayerElement and FaceGeometry element
// has to be inside the oomph namespace
namespace oomph
{
  //======start_of_MyRefineableQTaylorHoodElement=============================
  /// Overloaded element that allows projection of use as PML element
  //==========================================================================
  class MyRefineableQTaylorHoodElement :
    public virtual PMLElementBase<2>,
    public virtual RefineableQTaylorHoodElement<2>
  {
  public:
    /// Constructor
    MyRefineableQTaylorHoodElement() {}

    ///  Pure virtual function in which we specify the
    /// values to be pinned (and set to zero) on the outer edge of
    /// the "pml" layer. None since we're not using this as pml functionality
    void values_to_be_pinned_on_outer_pml_boundary(Vector<unsigned>& values_to_pin)
    {}
  };

  //======start_of_PMLLayerElement============================================
  /// Policy class defining the elements to be used in the PML layers. Same!
  //==========================================================================
  template<>
  class PMLLayerElement<MyRefineableQTaylorHoodElement> :
    public virtual MyRefineableQTaylorHoodElement
  {
  public:
    ///  Constructor: Call the constructor for the
    /// appropriate QElement
    PMLLayerElement() : MyRefineableQTaylorHoodElement()
    {}
  };

  //======start_of_FaceGeometry===============================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //==========================================================================
  template<>
  class FaceGeometry<MyRefineableQTaylorHoodElement> :
    public virtual QElement<1,3>
  {
  public:
    ///  Constructor: Call the constructor for the 1D quadratic element
    FaceGeometry() : QElement<1,3>() {}
  };
} // End of namespace oomph

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//======start_of_OscillatingCylinder========================================
///  Oscillating cylinder class
//==========================================================================
class OscillatingCylinder : public GeomObject
{
public:

  ///  Constructor: Pass radius, amplitude of the motion, the excitation
  /// wavelength and pointer to time object.
  OscillatingCylinder(double* radius_pt,
                      double* amplitude_pt,
                      Time* time_pt,
                      const bool& use_transition_phase=true) :
    GeomObject(1,2),
    Radius_pt(radius_pt),
    Amplitude_pt(amplitude_pt),
    Time_pt(time_pt),
    Use_transition_phase(use_transition_phase)
  {}

  /// Destructor: Empty
  virtual ~OscillatingCylinder() {}

  /// Access function for the amplitude (lvalue)
  double& amplitude()
  {
    // Return the value of the amplitude
    return *Amplitude_pt;
  } // End of amplitude

  ///  Current position vector to material point at Lagrangian
  /// coordinate xi (steady version)
  void position(const Vector<double>& xi,
                Vector<double>& r) const
  {
    // X-coordinate
    r[0]=(*Radius_pt)*cos(xi[0]);

    // Y-coordinate
    r[1]=(*Radius_pt)*sin(xi[0]);
  } // End of position

  ///  Current position vector to material point at Lagrangian
  /// coordinate xi (unsteady version). Implementation includes a
  /// transition phase where the cylinder oscillates to a smaller
  /// amplitude than the target value. Used to ensure that the solution
  /// isn't drastically different to that at the next time step. This
  /// can be disabled by setting Use_transition_phase to false.
  void position(const unsigned& t,
                const Vector<double>& xi,
                Vector<double>& r) const
  {
    // Calculate the coordinate before translation
    position(xi,r);

    // Get current time
    double time=Time_pt->time(t);

    // Scaling factor
    double arg=2.0*MathematicalConstants::Pi;

    // Calculate the translation
    double translation=(*Amplitude_pt)*sin(arg*time);

    // If the user wishes to use a transition phase
    if (Use_transition_phase)
    {
      // Time for transition (equal to two periods)
      double t_transition=2.0;

      // If we're not past the transition point
      if (time<t_transition)
      {
        // Scale the translation
        translation*=(time/t_transition)*exp(1.0-(time/t_transition));
      }
    } // if (Use_transition_phase)

    // Update the y-coordinate
    r[1]+=translation;
  } // End of position

private:

  /// Radius of the cylinder
  double* Radius_pt;

  /// Non-dimensionalised amplitude of the cylinder motion
  double* Amplitude_pt;

  /// Pointer to the current time in the problem
  Time* Time_pt;

  ///  Boolean variable to indicate whether or not to use a transition
  /// phase at the start of the simulation. Defaults to true.
  bool Use_transition_phase;
}; // End of OscillatingCylinder class

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//======start_of_GlobalParameters_namespace=================================
/// Global parameters
//==========================================================================
namespace GlobalParameters
{
  ///------------------------------------------------CYLINDER MOTION------
  ///  Amplitude of the cylinder motion used by Williamson &
  /// Roshko (1988). Probably best to fix the wavelength and vary this
  /// to get different wake modes. Also the easiest way to get different
  /// wake patterns to shown in Leontini et al. (2006). SIDE NOTE: as
  /// D=1 in our mesh this is also the dimensionless amplitude.
  double Amplitude=0.50;

  /// The ratio T_e/T_s
  double Period_ratio=1.0;
  ///------------------------------------------------CYLINDER MOTION------

  ///---------------------------------------NAVIER-STOKES PARAMETERS------
  /// The current Reynolds number
  double Re=0.0;

  /// Reynolds number for unsteady run
  double Re_target=20.0;

  /// The default Strouhal number (overloaded with input value if given)
  double St=1.0;

  /// The Womersley number
  double ReSt=Re*St;

  ///  Function to calculate the appropriate Strouhal number for
  /// the simulation.
  /// Some noteworthy Re-St values:
  ///                  (1)  Re=100 --> St=0.1643
  ///                  (2)  Re=200 --> St=0.198
  /// NOTE: The Re-St values for 46<Re<180 can be found in:
  ///         Williamson, C.H.K, (1988)."Defining a universal and
  ///         continuous Strouhal–Reynolds number relationship for
  ///         the laminar vortex shedding of a circular cylinder".
  double calculate_strouhal_number(const double& re)
  {
    // The min. Reynolds number
    double min_re=46.0;

    // The max. Reynolds number
    double max_re=180.0;

    // The first coefficient of the Re-St polynomial
    double a=-3.3265;

    // The second coefficient of the Re-St polynomial
    double b=0.1816;

    // The third coefficient of the Re-St polynomial
    double c=0.00016;

    // If we're above the maximum Reynolds number
    if (re>max_re)
    {
      // Throw an error
      throw OomphLibError("Don't know what to do for this Reynolds number!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // If we're below the minimum Reynolds number
    else if (re<min_re)
    {
      // Just return the Strouhal value at the minimum Reynolds number
      return a/min_re+b+c*min_re;
    }
    // Otherwise, use the relationship in the Williamson paper
    else
    {
      // Return the Strouhal value at this Reynolds number
      return a/re+b+c*re;
    }
  } // End of calculate_strouhal_number

  ///  Update physical parameters. This updates:
  ///               (1) The Reynolds number, and;
  ///               (2) The Strouhal number,
  /// and should ALWAYS be called after either the Reynolds number or
  /// Period_ratio value changes.
  void update_simulation_parameters()
  {
    // Update the Strouhal number ourselves to ensure "lock-in"
    St=calculate_strouhal_number(Re)/Period_ratio;

    // Calculate the Womersley number
    ReSt=Re*St;
  } // End of update_simulation_parameters

  /// Document the value of the Reynolds number and Womersley number
  void doc_navier_stokes_parameters()
  {
    // Tell the user
    oomph_info << ANSIEscapeCode::Red
               << "Solving for (Re,ReSt): "
               << ANSIEscapeCode::Reset
               << "("
               << GlobalParameters::Re << "," << GlobalParameters::ReSt
               << ")" << std::endl;
  } // End of doc_navier_stokes_parameters
  ///---------------------------------------NAVIER-STOKES PARAMETERS------

  ///------------------------------------TIME-INTEGRATION PARAMETERS------
  /// Number of periods for unsteady run
  unsigned N_period_unsteady=1;

  /// Number of timesteps per period for unsteady run
  unsigned N_step_per_period_unsteady=100;
  ///------------------------------------TIME-INTEGRATION PARAMETERS------

  ///----------------------------------------------DOMAIN PROPERTIES------
  /// Pointer to the cylinder
  OscillatingCylinder* Cylinder_pt=0;

  /// Height of domain
  double Height=20.0;

  /// X-coordinate of upstream end of domain
  double X_left=-10.0;

  /// X-coordinate of downstream end of domain
  double X_right=40.0;

  /// Side-length of the square box in the mesh surrounding the cylinder
  double Length_of_central_box=10.0;

  /// Radius of the cylinder
  double Radius=0.5;

  /// The radius of the annular region surrounding the cylinder
  /// NOTE: The annular rings are used to resolve the boundary layers so they
  /// should not be made too large (hence the use of std::min).
  double Annular_region_radius=
    std::min(Radius+1.0,Radius+0.5*((0.5*Length_of_central_box)-Radius));

  /// Number of uniform refinements before any solve
  unsigned N_uniform_refinement_before_solve=2;
  ///----------------------------------------------DOMAIN PROPERTIES------

  ///------------------------------------------DOCUMENTATION HELPERS------
  /// The number of plot points in each direction
  unsigned N_plot_point=2;

  /// Doc info object
  DocInfo Doc_info;

  /// Document the maximum deformation inside the central box
  void doc_maximum_central_box_deformation()
  {
    // Calculate the distance from the edge of the annular ring to the
    // box boundary. NOTE: We check from the annular ring because the
    // region between the cylinder and annular ring is made rigid so
    // no compression occurs there.
    double compression_region_width=(Length_of_central_box/2.0-
                                     Annular_region_radius);

    // Calculate the current deformation of the inner box
    double deformation_ratio=((compression_region_width-Amplitude)/
                              compression_region_width);

    // If the deformation is too large for the mesh
    if (deformation_ratio<0.0)
    {
      // Used to create an error message
      std::ostringstream error_message_stream;

      // Create an error message
      error_message_stream << "The cylinder amplitude exceeds the size of "
                           << "the central box! Make the box larger!"
                           << std::endl;

      // Throw an error to the user
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // If the deformation is large then warn the user
    else if (deformation_ratio<0.5)
    {
      // Used to create a warning message
      std::ostringstream warning_message_stream;

      // Create a warning message
      warning_message_stream << "Maximal mesh compression results in elements "
                             << "being reduced\nto "
                             << deformation_ratio*100.0
                             << "% of their original width. It is therefore\n"
                             << "recommended that the central box be "
                             << "made larger." << std::endl;

      // Throw a warning to the user
      OomphLibWarning(warning_message_stream.str(),
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      // Document the deformation (special ASCII characters used to output
      // the text in bold red)
      oomph_info << "\033[1;31m"
                 << "\nMaximum element compression ratio inside central box: "
                 << "\033[0m"
                 << deformation_ratio*100.0 << "%\n" << std::endl;
    }
  } // End of doc_maximum_central_box_deformation
  ///------------------------------------------DOCUMENTATION HELPERS------

  ///--------------------------------------------------MISCELLANEOUS------
  /// Function to round a double to the nearest integral value
  double round(const double& d)
  {
    // Round it
    return std::floor(d+0.5);
  } // End of round
  ///--------------------------------------------------MISCELLANEOUS------
} // End of GlobalParameters

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//======start_of_FlowAroundCylinderProblem_class============================
/// Flow around a cylinder in rectangular domain
//==========================================================================
template<class ELEMENT>
class FlowAroundCylinderProblem : public Problem
{
public:

  /// Constructor:
  FlowAroundCylinderProblem();

  /// Destructor: delete all dynamically allocated data
  ~FlowAroundCylinderProblem();

  /// Apply boundary conditions
  void apply_boundary_conditions();

  /// Set the chosen solver and preconditioner
  void set_up_solver_and_preconditioner();

  /// Timestep the problem
  void unsteady_simulation();

  /// Update the problem specs before solve
  void actions_before_newton_solve() {}

  /// Update the problem specs after solve (empty)
  void actions_after_newton_solve() {}

  /// Actions before adapt: empty
  void actions_before_adapt() {}

  /// After adaptation: Unpin pressure and pin redundant pressure dofs.
  void actions_after_adapt();

  ///  Actions to complete before an implicit timestep:
  /// Provide a perturbation to the solution
  void actions_before_implicit_timestep();

  /// Doc the solution
  void doc_solution(const bool& in_unsteady=false);

  ///  Output the velocity components at a point on the centerline
  /// to observe the oscillating behaviour without having to output all
  /// the data (i.e. V1, V2 etc.):
  void doc_trace_node_velocities_and_pressure();

  /// Dump the entire problem
  void dump_problem(char* filename) const;

  /// Read the entire problem
  void read_problem(char* restart_file);

private:

  /// Pointer to the mesh with oscillating cylinder
  RefineableQuadMeshWithMovingCylinder<ELEMENT>* Bulk_mesh_pt;

  /// Oomph-lib iterative linear solver
  IterativeLinearSolver* Solver_pt;

  /// LSC Preconditioner for the linear solve
  NavierStokesSchurComplementPreconditioner* Prec_pt;

  /// Inexact solver for F block
  Preconditioner* F_matrix_preconditioner_pt;

  // Enumeration of the boundaries of the space-time mesh
  enum
  {
    Lower_wall_boundary_id=0,
    Outflow_boundary_id=1,
    Upper_wall_boundary_id=2,
    Inflow_boundary_id=3,
    Cylinder_surface_boundary_id=4
  };
};


//=====start_of_FlowAroundCylinderProblem===================================
/// Constructor
//==========================================================================
template<class ELEMENT>
FlowAroundCylinderProblem<ELEMENT>::FlowAroundCylinderProblem() :
  Bulk_mesh_pt(0),
  Solver_pt(0),
  Prec_pt(0),
  F_matrix_preconditioner_pt(0)
{
  // Use BDF2
  add_time_stepper_pt(new BDF<2>);

  // Indicate that we're solving a steady problem
  time_stepper_pt(0)->make_steady();

  //-----------------
  // Set up the mesh:
  //-----------------
  // Create a OscillatingCylinder object to define the cylinder motion
  GlobalParameters::Cylinder_pt=
    new OscillatingCylinder(&GlobalParameters::Radius,
                            &GlobalParameters::Amplitude,
                            time_pt());

  // Make a new mesh and assign its pointer
  Bulk_mesh_pt=new RefineableQuadMeshWithMovingCylinder<ELEMENT>(
    GlobalParameters::Cylinder_pt,
    GlobalParameters::Annular_region_radius,
    GlobalParameters::Length_of_central_box,
    GlobalParameters::X_left,
    GlobalParameters::X_right,
    GlobalParameters::Height,
    time_stepper_pt());

  // Now assign its pointer
  Problem::mesh_pt()=Bulk_mesh_pt;

  // Set the error estimator: incase we want to use adaptive mesh refinement
  Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

  // No limit on the number of refinements
  Bulk_mesh_pt->max_refinement_level()=100;

  // Set the maximum error bound
  Bulk_mesh_pt->max_permitted_error()=1.0e-04;

  // Set the minimum error bound
  Bulk_mesh_pt->min_permitted_error()=1.0e-08;

  //-----------------------------
  // Set the boundary conditions:
  //-----------------------------
  // Apply the boundary conditions
  apply_boundary_conditions();

  //-----------------------------------------------------------------
  // Complete the build of all elements so they are fully functional:
  //-----------------------------------------------------------------
  // Make it pseudo-traction-free
  NavierStokesEquations<2>::Gamma[0]=0.0;
  NavierStokesEquations<2>::Gamma[1]=0.0;

  // Pin redundant pressure dofs
  RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());

  // How many elements are there in the mesh?
  unsigned n_element=Bulk_mesh_pt->nelement();

  // Loop over the elements in the mesh
  for (unsigned e=0; e<n_element; e++)
  {
    // Upcast the e-th element in the mesh
    ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // Set the Reynolds number pointer
    el_pt->re_pt()=&GlobalParameters::Re;

    // Set the Womersley number (same as Re since St=1) pointer
    el_pt->re_st_pt()=&GlobalParameters::ReSt;
  }

  // Attach the boundary conditions to the mesh
  oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
} // End of FlowAroundCylinderProblem


//=====start_of_~FlowAroundCylinderProblem==================================
/// The destructor. Delete all dynamically allocated data.
//==========================================================================
template<class ELEMENT>
FlowAroundCylinderProblem<ELEMENT>::~FlowAroundCylinderProblem()
{
  // If we set up a preconditioner for the momentum block
  if (F_matrix_preconditioner_pt!=0)
  {
    // Kill it
    delete F_matrix_preconditioner_pt;

    // Make it a null pointer
    F_matrix_preconditioner_pt=0;
  }

  // If we're not using the default linear solver
  if (Prec_pt!=0)
  {
    // Kill the preconditioner
    delete Prec_pt;

    // Make it a null pointer
    Prec_pt=0;
  }

  // If we're not using the default linear solver
  if (Solver_pt!=0)
  {
    // Kill the linear solver
    delete Solver_pt;

    // Make it a null pointer
    Solver_pt=0;

    // Make the Problem's linear solver pointer null too
    linear_solver_pt()=0;
  }

  // Delete the error estimator
  delete Bulk_mesh_pt->spatial_error_estimator_pt();

  // Make it a null pointer
  Bulk_mesh_pt->spatial_error_estimator_pt()=0;

  // Delete the mesh
  delete Bulk_mesh_pt;

  // Make it a null pointer
  Bulk_mesh_pt=0;
} // End of ~FlowAroundCylinderProblem


//=====start_of_apply_boundary_conditions===================================
/// Apply boundary conditions (BX = boundary X):
///   -- Pin both velocities on the cylinder (B4) and the inlet (B3).
///   -- Zero vertical velocity on the upper (B2) and lower (B0) boundaries.
///   -- Do nothing (non-stress divergence form) on the outlet (B1):
//==========================================================================
template<class ELEMENT>
void FlowAroundCylinderProblem<ELEMENT>::apply_boundary_conditions()
{
  // Get the number of boundaries in the mesh
  unsigned n_boundary=Bulk_mesh_pt->nboundary();

  // Loop over the boundaries
  for (unsigned b=0; b<n_boundary; b++)
  {
    // Find the number of nodes on the b-th boundary
    unsigned n_boundary_node=Bulk_mesh_pt->nboundary_node(b);

    // Inlet and side walls (tow-tank conditions)
    if ((b==Lower_wall_boundary_id)||
        (b==Upper_wall_boundary_id)||
        (b==Inflow_boundary_id)||
        (b==Cylinder_surface_boundary_id))
    {
      // Loop over the nodes
      for (unsigned n=0; n<n_boundary_node; n++)
      {
        // Pointer to the boundary node
        Node* boundary_node_pt=Bulk_mesh_pt->boundary_node_pt(b,n);

        // Pin the horizontal velocity
        boundary_node_pt->pin(0);

        // Pin the vertical velocity
        boundary_node_pt->pin(1);

        // Cylinder surface: enforce a no-slip condition
        // NOTE: the no-slip condition will be updated before each timestep in
        // response to the cylinder motion using the FSI helpers.
        if (b==Cylinder_surface_boundary_id)
        {
          // Set horizontal velocity to zero
          boundary_node_pt->set_value(0,0.0);

          // Set the vertical velocity to zero too
          boundary_node_pt->set_value(1,0.0);
        }
        // Inflow and sidewalls: enforce tow-tank boundary conditions
        else
        {
          // Tangential inflow
          boundary_node_pt->set_value(0,1.0);

          // No vertical velocity
          boundary_node_pt->set_value(1,0.0);
        }
      } // for (unsigned b=0;b<n_boundary_node;b++)
    } // if ((b==Lower_wall_boundary_id)||...
  } // for (unsigned b=0; b<n_boundary; b++)
} // End of apply_boundary_conditions


//=====start_of_set_up_solver_and_preconditioner============================
/// Set up the linear solver and preconditioner
//==========================================================================
template<class ELEMENT>
void FlowAroundCylinderProblem<ELEMENT>::set_up_solver_and_preconditioner()
{
  // Use GMRES
  Solver_pt=new GMRES<CRDoubleMatrix>;

  // Set linear solver
  linear_solver_pt()=Solver_pt;

  // Create and assign the Navier-Stokes preconditioner
  Prec_pt=new NavierStokesSchurComplementPreconditioner(this);

  // Pass the Navier-Stokes mesh to the preconditioner
  Prec_pt->set_navier_stokes_mesh(this->Bulk_mesh_pt);

  // Pass the preconditioner to the solver
  Solver_pt->preconditioner_pt()=Prec_pt;

  // By default, the LSC Preconditioner uses SuperLU as an exact
  // preconditioner (i.e. a solver) for the momentum and Schur complement
  // blocks. Can overwrite this by passing pointers to other preconditioners
  // that perform the (approximate) solves of these blocks:

  // Create internal preconditioners used on momentum block
  F_matrix_preconditioner_pt=new BlockDiagonalPreconditioner<CRDoubleMatrix>;

  // Assign the subsidiary preconditioner to the linear solver preconditioner
  Prec_pt->set_f_preconditioner(F_matrix_preconditioner_pt);
} // End of set_up_solver_and_preconditioner


//=====start_of_unsteady_simulation=========================================
/// Timestep the problem
//==========================================================================
template<class ELEMENT>
void FlowAroundCylinderProblem<ELEMENT>::unsteady_simulation()
{
  // The FSI_Functions namespace enforces the no-slip condition on the
  // cylinder boundary as the nodes on this boundary move with time. This
  // requires a handle to the Strouhal number value, so set it here
  FSI_functions::Strouhal_for_no_slip=GlobalParameters::St;

  // Make the problem unsteady again
  time_stepper_pt()->undo_make_steady();

  // Set the period (working on the non-dimensionalised time)
  double period=1.0;

  // Timestep
  double dt=period/double(GlobalParameters::N_step_per_period_unsteady);

  // Initialise the timestep value
  initialise_dt(dt);

  // Doc the number of timesteps per period and timestep size
  oomph_info << "\nTimestepping with "
             << GlobalParameters::N_step_per_period_unsteady
             << " steps per period, corresponding to a timestep of "
             << dt << std::endl;

  // Make the problem unsteady again
  time_stepper_pt()->undo_make_steady();

  // Assign history values for an impulsive start
  assign_initial_values_impulsive();

  // Reset the documentation counter so the file indexing (re)starts from zero
  GlobalParameters::Doc_info.number()=0;

  // Set up the solver and preconditioner
  set_up_solver_and_preconditioner();

  // Make sure that the Reynolds number is exactly equal to our target value
  GlobalParameters::Re=GlobalParameters::Re_target;

  // Document the Reynolds number and Womersley number
  GlobalParameters::doc_navier_stokes_parameters();

  // Initialise the number of timesteps
  unsigned n_timestep=
    GlobalParameters::N_period_unsteady*GlobalParameters::round(period/dt);

  // Tell the user
  oomph_info << "\nTotal number of timesteps: " << n_timestep << std::endl;

  // Target global temporal error
  double epsilon_t=1.0e-04;

  // Max. number of spatial adaptations
  unsigned max_adapt=0;

  // Inform the user
  oomph_info << "\nRunning until time: " << GlobalParameters::N_period_unsteady
             << std::endl;

  // Loop over the timesteps
  for (unsigned t=1; t<=n_timestep; t++)
  {
    // Run the bastard...
    doubly_adaptive_unsteady_newton_solve(dt,epsilon_t,max_adapt,false);

    // Tell the user
    oomph_info << "Obtained solution at timestep " << t
               << " out of " << n_timestep << "." << std::endl;
  }
} // End of unsteady_simulation


//=====start_of_actions_after_adapt=========================================
/// Actions after adaptation. Upin pressure dofs, pin all redundant pressure
/// dofs, apply boundary conditions and pin the smoothed vorticity in any
/// newly created elements.
//==========================================================================
template<class ELEMENT>
void FlowAroundCylinderProblem<ELEMENT>::actions_after_adapt()
{
  // Unpin all pressure dofs
  RefineableNavierStokesEquations<2>::
  unpin_all_pressure_dofs(Bulk_mesh_pt->element_pt());

  // Pin redundant pressure dofs
  RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());

  // Call the boundary conditions function
  apply_boundary_conditions();
} // End of actions_after_adapt


//=====start_of_actions_before_implicit_timestep============================
/// Actions before implicit timestep: Update the nodal positions and update
/// the no-slip condition on the cylinder boundary
//==========================================================================
template<class ELEMENT>
void FlowAroundCylinderProblem<ELEMENT>::actions_before_implicit_timestep()
{
  // Update the domain shape
  Bulk_mesh_pt->node_update();

  // ID of the cylinder boundary
  unsigned b=4;

  // Find the number of nodes on the cylinder boundary
  unsigned num_nod=Bulk_mesh_pt->nboundary_node(b);

  // Cylinder boundary: No slip; this implies that the velocity needs to be
  // updated in response to wall motion
  for (unsigned n=0; n<num_nod; n++)
  {
    // Which node are we dealing with?
    Node* node_pt=Bulk_mesh_pt->boundary_node_pt(b,n);

    // Apply no slip
    FSI_functions::apply_no_slip_on_moving_wall(node_pt);
  }
} // End of actions_before_implicit_timestep()


//=====start_of_doc_solution================================================
/// Document the solution
//==========================================================================
template<class ELEMENT>
void FlowAroundCylinderProblem<ELEMENT>::doc_solution(const bool& in_unsteady)
{
  // Tell the user
  oomph_info << "Documentation step: " << GlobalParameters::Doc_info.number()
             << std::endl;

  // Only output on root processor
  if (this->communicator_pt()->my_rank()==0)
  {
    // Create an output stream
    std::ofstream some_file;

    // Allocate space for the filename
    char filename[10000];

    // Number of plot points
    unsigned n_plot=GlobalParameters::N_plot_point;

    // If we're in the unsteady simulation
    if (in_unsteady)
    {
      // Output solution
      sprintf(filename,"%s/unsteady%i.dat",
              GlobalParameters::Doc_info.directory().c_str(),
              GlobalParameters::Doc_info.number());
    }
    else
    {
      // Output solution
      sprintf(filename,"%s/soln%i.dat",
              GlobalParameters::Doc_info.directory().c_str(),
              GlobalParameters::Doc_info.number());
    }

    // Open a file with the formed filename
    some_file.open(filename);

    // Output solution
    Bulk_mesh_pt->output(some_file,n_plot);

    // Close the file
    some_file.close();
  } // if (this->communicator_pt()->my_rank()==0)

  // Increment the documentation number
  GlobalParameters::Doc_info.number()++;
} // End of doc_solution

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//==========================================================================
/// Driver
//==========================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);

  // Switch off output modifier
  oomph_info.output_modifier_pt()=&default_output_modifier;

  // Switch off oomph_info output for all processors but rank 0
  if (MPI_Helpers::communicator_pt()->my_rank()!=0)
  {
    oomph_info.stream_pt()=&oomph_nullstream;
    OomphLibWarning::set_stream_pt(&oomph_nullstream);
    OomphLibError::set_stream_pt(&oomph_nullstream);
  }
  else
  {
    oomph_info << "\n\n====================================================\n"
               << "Number of processors: "
               << MPI_Helpers::communicator_pt()->nproc()
               << "\n====================================================\n"
               << std::endl;
  }
#endif

  // Start the timer
  double timer_s=TimingHelpers::timer();

  // Generic update function that calls all the appropriate update functions
  GlobalParameters::update_simulation_parameters();

  // Output directory
  GlobalParameters::Doc_info.set_directory("RESLT");

  // Build the problem
  FlowAroundCylinderProblem<MyRefineableQTaylorHoodElement> problem;

  // Refine the problem a couple of times
  for (unsigned i=0; i<GlobalParameters::N_uniform_refinement_before_solve; i++)
  {
    // Refine uniformly
    problem.refine_uniformly();
  }

  // Number of steps to increment the Reynolds number to the target value
  unsigned n_re_step=1;

  // Calculate the increment required to get to the Reynolds number
  double re_increment=GlobalParameters::Re_target/n_re_step;

  //-----------------------------
  // Crank up the Reynolds number
  //-----------------------------
  // Loop over the Reynolds increments
  for (unsigned i=0; i<n_re_step+1; i++)
  {
    // Calculate the Reynolds number
    GlobalParameters::Re=double(i)*re_increment;

    // Update the Womersley number
    GlobalParameters::update_simulation_parameters();

    // Output the current N.St. parameter values
    GlobalParameters::doc_navier_stokes_parameters();

    // Solve the problem
    problem.newton_solve();
  }

  // Document the solution
  problem.doc_solution();

  //------------------------------
  // Start the unsteady simulation
  //------------------------------
  // Run the simulation
  problem.unsteady_simulation();

  // Document the solution
  problem.doc_solution(true);

  // Tell the user we're done and how long the whole thing took
  oomph_info << "\nUnsteady simulation complete!"
             << "\nTotal time for full simulation [sec]: "
             << TimingHelpers::timer()-timer_s << "\n" << std::endl;
}
