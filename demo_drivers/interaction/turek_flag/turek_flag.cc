//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
//Generic stuff
#include "generic.h"

// Solid mechanics
#include "solid.h"

// Navier Stokes
#include "navier_stokes.h"

// Multiphysics for FSI preconditioner
#include "multi_physics.h"

// The meshes
#include "meshes/cylinder_with_flag_mesh.h"
#include "meshes/rectangular_quadmesh.h"


using namespace std;

using namespace oomph;


//=====start_of_global_parameters=================================
/// Global variables
//================================================================
namespace Global_Parameters
{

 /// Default case ID
 string Case_ID="FSI1";

 /// Reynolds number (default assignment for FSI1 test case)
 double Re=20.0; 

 /// Strouhal number (default assignment for FSI1 test case)
 double St=0.5;

 /// Product of Reynolds and Strouhal numbers (default
 /// assignment for FSI1 test case)
 double ReSt=10.0;

 /// FSI parameter (default assignment for FSI1 test case)
 double Q=1.429e-6;
 
 /// Density ratio (solid to fluid; default assignment for FSI1
 /// test case)
 double Density_ratio=1.0; 

 /// Height of flag
 double H=0.2;
 
 /// x position of centre of cylinder
 double Centre_x=2.0;

 /// y position of centre of cylinder
 double Centre_y=2.0;

 /// Radius of cylinder
 double Radius=0.5;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Timescale ratio for solid (dependent parameter
 /// assigned in set_parameters())
 double Lambda_sq=0.0;

 /// Timestep
 double Dt=0.1;

 /// Ignore fluid (default assignment for FSI1 test case)
 bool Ignore_fluid_loading=false;

 /// Elastic modulus
 double E=1.0;

 /// Poisson's ratio
 double Nu=0.4; 

 /// Non-dim gravity (default assignment for FSI1 test case)
 double Gravity=0.0; 

 /// Non-dimensional gravity as body force
 void gravity(const double& time, 
              const Vector<double> &xi,
              Vector<double> &b)
 {
  b[0]=0.0;
  b[1]=-Gravity;
 }
 
 /// Period for ramping up in flux
 double Ramp_period=2.0;
 
 /// Min. flux 
 double Min_flux=0.0; 
 
 /// Max. flux
 double Max_flux=1.0;

 /// Flux increases between Min_flux and Max_flux over 
 /// period Ramp_period
 double flux(const double& t)
 {    
  if (t<Ramp_period)
   {
    return Min_flux+(Max_flux-Min_flux)*
     0.5*(1.0-cos(MathematicalConstants::Pi*t/Ramp_period));
   }
  else
   {
    return Max_flux;
   }
 } // end of specification of ramped influx
 

 /// Set parameters for the various test cases
 void set_parameters(const string& case_id)
 {

  // Remember which case we're dealing with
  Case_ID=case_id;

  // Setup independent parameters depending on test case
  if (case_id=="FSI1")
   {
    // Reynolds number based on diameter of cylinder
    Re=20.0;

    // Strouhal number based on timescale of one second
    St=0.5;

    // Womersley number
    ReSt=Re*St;

    // FSI parameter
    Q=1.429e-6;
    
    // Timestep -- aiming for about 40 steps per period
    Dt=0.1;

    // Density ratio
    Density_ratio=1.0;

    // Gravity
    Gravity=0.0;
    
    // Max. flux
    Max_flux=1.0;

    // Ignore fluid
    Ignore_fluid_loading=false;
 
    // Compute dependent parameters
    
    // Timescale ratio for solid 
    Lambda_sq=Re*Q*Density_ratio*St*St;
  }
  else if (case_id=="FSI2")
   {
    // Reynolds number based on diameter of cylinder
    Re=100.0;

    // Strouhal number based on timescale of one second
    St=0.1;

    // Womersley number
    ReSt=Re*St;

    // FSI parameter
    Q=7.143e-6;
    
    // Timestep -- aiming for about 40 steps per period
    Dt=0.01;

    // Density ratio
    Density_ratio=10.0;

    // Gravity
    Gravity=0.0;
    
    // Max. flux
    Max_flux=1.0;

    // Ignore fluid
    Ignore_fluid_loading=false;

    // Compute dependent parameters
    
    // Timescale ratio for solid 
    Lambda_sq=Re*Q*Density_ratio*St*St;
   }
  else if (case_id=="FSI3")
   {
    // Reynolds number based on diameter of cylinder
    Re=200.0;

    // Strouhal number based on timescale of one second
    St=0.05;

    // Womersley number
    ReSt=Re*St;

    // FSI parameter
    Q=3.571e-6;

    // Timestep -- aiming for about 40 steps per period
    Dt=0.005;

    // Density ratio
    Density_ratio=1.0;

    // Gravity
    Gravity=0.0;
    
    // Max. flux
    Max_flux=1.0;

    // Ignore fluid
    Ignore_fluid_loading=false;

    // Compute dependent parameters
    
    // Timescale ratio for solid 
    Lambda_sq=Re*Q*Density_ratio*St*St;
   } 
  else if (case_id=="CSM1")
   {
    // Reynolds number based on diameter of cylinder
    Re=0.0;

    // Strouhal number based on timescale of one second
    // (irrelevant here)
    St=0.0;

    // Womersley number
    ReSt=Re*St;

    // FSI parameter
    Q=0.0;

    // Timestep -- irrelevant here
    Dt=0.005;

    // Density ratio (switch off wall inertia)
    Density_ratio=0.0;

    // Gravity
    Gravity=1.429e-4;

    // Max. flux
    Max_flux=0.0;
    
    // Ignore fluid
    Ignore_fluid_loading=true;

    // Compute dependent parameters
    
    // Timescale ratio for solid 
    Lambda_sq=0.0;
   }
  else if (case_id=="CSM3")
   {
    // Reynolds number based on diameter of cylinder
    Re=0.0;

    // Strouhal number based on timescale of one second
    // (irrelevant here)
    St=0.0;

    // Womersley number
    ReSt=Re*St;

    // Timestep -- aiming for about 40 steps per period
    Dt=0.025;

    // FSI parameter
    Q=0.0;

    // Density ratio (not used here)
    Density_ratio=0.0;

    // Gravity
    Gravity=1.429e-4;

    // Max. flux
    Max_flux=0.0;
    
    // Ignore fluid
    Ignore_fluid_loading=true;

    // Compute dependent parameters
    
    // Set timescale ratio for solid based on the
    // observed period of oscillation of about 1 sec)
    Lambda_sq=7.143e-6;
   }
  else
   {
    std::cout << "Wrong case_id: " << case_id << std::endl;
    exit(0);
   }
 
  
  // Ramp period (20 timesteps)
  Ramp_period=Dt*20.0;

  // "Big G" Linear constitutive equations:
  Constitutive_law_pt = new GeneralisedHookean(&Nu,&E);
  
  // Doc
  oomph_info << std::endl;
  oomph_info << "-------------------------------------------" 
             << std::endl;
  oomph_info << "Case: " << case_id << std::endl;
  oomph_info << "Re            = " << Re << std::endl;
  oomph_info << "St            = " << St << std::endl;
  oomph_info << "ReSt          = " << ReSt << std::endl;
  oomph_info << "Q             = " << Q << std::endl;
  oomph_info << "Dt            = " << Dt << std::endl;
  oomph_info << "Ramp_period   = " << Ramp_period << std::endl;
  oomph_info << "Max_flux      = " << Max_flux << std::endl;
  oomph_info << "Density_ratio = " << Density_ratio << std::endl;
  oomph_info << "Lambda_sq     = " << Lambda_sq << std::endl;
  oomph_info << "Gravity       = " << Gravity << std::endl;
  oomph_info << "Ignore fluid  = " << Ignore_fluid_loading<< std::endl;
  oomph_info << "-------------------------------------------"
             << std::endl << std::endl;
 }

}// end_of_namespace



/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//====start_of_problem_class=========================================== 
/// Problem class
//====================================================================== 
template< class FLUID_ELEMENT,class SOLID_ELEMENT >
class TurekProblem : public Problem
{

public:

 /// Constructor: Pass length and height of domain
 TurekProblem(const double &length, const double &height);
 
 /// Access function for the fluid mesh 
 RefineableAlgebraicCylinderWithFlagMesh<FLUID_ELEMENT>* fluid_mesh_pt() 
  { return Fluid_mesh_pt;}

 /// Access function for the solid mesh
 ElasticRefineableRectangularQuadMesh<SOLID_ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

 /// Access function for the i-th mesh of FSI traction elements
 SolidMesh*& traction_mesh_pt(const unsigned& i)
  {return Traction_mesh_pt[i];} 
 
 /// Actions after adapt: Re-setup the fsi lookup scheme
 void actions_after_adapt();

 /// Doc the solution
 void doc_solution(DocInfo& doc_info, ofstream& trace_file);

 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve(){}

 /// Update the (dependent) fluid node positions following the
 /// update of the solid variables before performing Newton convergence
 /// check
 void actions_before_newton_convergence_check();

 /// Update the time-dependent influx
 void actions_before_implicit_timestep();

private:

 /// Create FSI traction elements
 void create_fsi_traction_elements();

 /// Pointer to solid mesh
 ElasticRefineableRectangularQuadMesh<SOLID_ELEMENT>* Solid_mesh_pt;

 /// Pointer to fluid mesh
 RefineableAlgebraicCylinderWithFlagMesh<FLUID_ELEMENT>* Fluid_mesh_pt;

 /// Vector of pointers to mesh of FSI traction elements
 Vector<SolidMesh*> Traction_mesh_pt;

 /// Combined mesh of traction elements -- only used for documentation
 SolidMesh* Combined_traction_mesh_pt;

 /// Overall height of domain
 double Domain_height;

 /// Overall length of domain
 double Domain_length;

 /// Pointer to solid control node
 Node* Solid_control_node_pt;

 /// Pointer to fluid control node
 Node* Fluid_control_node_pt;
 
};// end_of_problem_class


//=====start_of_constructor============================================= 
/// Constructor: Pass length and height of domain
//====================================================================== 
template< class FLUID_ELEMENT,class SOLID_ELEMENT >
TurekProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
TurekProblem(const double &length,
             const double &height) :  Domain_height(height),
                                      Domain_length(length)
 
{
 // Increase max. number of iterations in Newton solver to
 // accomodate possible poor initial guesses
 Max_newton_iterations=20;
 Max_residuals=1.0e4;

 // Build solid mesh
 //------------------

 // # of elements in x-direction
 unsigned n_x=20;

 // # of elements in y-direction
 unsigned n_y=2;

 // Domain length in y-direction 
 double l_y=Global_Parameters::H;

 // Create the flag timestepper (consistent with BDF<2> for fluid)
 Newmark<2>* flag_time_stepper_pt=new Newmark<2>;
 add_time_stepper_pt(flag_time_stepper_pt); 

 /// Left point on centreline of flag so that the top and bottom
 /// vertices merge with the cylinder.
 Vector<double> origin(2);
 origin[0]=Global_Parameters::Centre_x+
  Global_Parameters::Radius*
  sqrt(1.0-Global_Parameters::H*Global_Parameters::H/
       (4.0*Global_Parameters::Radius*Global_Parameters::Radius));
 origin[1]=Global_Parameters::Centre_y-0.5*l_y;

 // Set length of flag so that endpoint actually stretches all the
 // way to x=6:
 double l_x=6.0-origin[0];

 //Now create the mesh
 solid_mesh_pt() = new ElasticRefineableRectangularQuadMesh<SOLID_ELEMENT>(
  n_x,n_y,l_x,l_y,origin,flag_time_stepper_pt);
 
 // Set error estimator for the solid mesh
 Z2ErrorEstimator* solid_error_estimator_pt=new Z2ErrorEstimator;
 solid_mesh_pt()->spatial_error_estimator_pt()=solid_error_estimator_pt;


 // Element that contains the control point
 FiniteElement* el_pt=solid_mesh_pt()->finite_element_pt(n_x*n_y/2-1);

 // How many nodes does it have?
 unsigned nnod=el_pt->nnode();

 // Get the control node
 Solid_control_node_pt=el_pt->node_pt(nnod-1);

 std::cout << "Coordinates of solid control point " 
           << Solid_control_node_pt->x(0) << " " 
           << Solid_control_node_pt->x(1) << " " << std::endl;

 // Refine the mesh uniformly
 solid_mesh_pt()->refine_uniformly();

 //Do not allow the solid mesh to be refined again
 solid_mesh_pt()->disable_adaptation();


 // Build mesh of solid traction elements that apply the fluid
 //------------------------------------------------------------
 // traction to the solid elements
 //-------------------------------

 // Create storage for Meshes of FSI traction elements at the bottom
 // top and left boundaries of the flag
 Traction_mesh_pt.resize(3);
 
 // Now construct the traction element meshes
 Traction_mesh_pt[0]=new SolidMesh;
 Traction_mesh_pt[1]=new SolidMesh;
 Traction_mesh_pt[2]=new SolidMesh;

 // Build the FSI traction elements
 create_fsi_traction_elements();

 // Loop over traction elements, pass the FSI parameter and tell them 
 // the boundary number in the bulk solid mesh -- this is required so 
 // they can get access to the boundary coordinates!
 for (unsigned bound=0;bound<3;bound++)
  {
   unsigned n_face_element = Traction_mesh_pt[bound]->nelement();
   for(unsigned e=0;e<n_face_element;e++)
    {
     //Cast the element pointer and specify boundary number
     FSISolidTractionElement<SOLID_ELEMENT,2>* elem_pt=
     dynamic_cast<FSISolidTractionElement<SOLID_ELEMENT,2>*>
      (Traction_mesh_pt[bound]->element_pt(e));

     // Specify boundary number
     elem_pt->set_boundary_number_in_bulk_mesh(bound);

     // Function that specifies the load ratios
     elem_pt->q_pt() = &Global_Parameters::Q;
    
    }
  } // build of FSISolidTractionElements is complete


 // Turn the three meshes of FSI traction elements into compound
 // geometric objects (one Lagrangian, two Eulerian coordinates)
 // that determine the boundary of the fluid mesh
 MeshAsGeomObject*
  bottom_flag_pt=
  new MeshAsGeomObject
  (Traction_mesh_pt[0]);
 
 MeshAsGeomObject* tip_flag_pt=
  new MeshAsGeomObject
  (Traction_mesh_pt[1]);
 
 MeshAsGeomObject* top_flag_pt=
  new MeshAsGeomObject
  (Traction_mesh_pt[2]);


 // Build fluid mesh
 //-----------------

 //Create a new Circle object as the central cylinder
 Circle* cylinder_pt = new Circle(Global_Parameters::Centre_x,
                                  Global_Parameters::Centre_y,
                                  Global_Parameters::Radius);

 // Allocate the fluid timestepper
 BDF<2>* fluid_time_stepper_pt=new BDF<2>;
 add_time_stepper_pt(fluid_time_stepper_pt);
 
 // Build fluid mesh
 Fluid_mesh_pt=
  new RefineableAlgebraicCylinderWithFlagMesh<FLUID_ELEMENT>
  (cylinder_pt, 
   top_flag_pt,
   bottom_flag_pt,
   tip_flag_pt,
   length, height, 
   l_x,Global_Parameters::H,
   Global_Parameters::Centre_x,
   Global_Parameters::Centre_y,
   Global_Parameters::Radius,
   fluid_time_stepper_pt);
 

 // I happen to have found out by inspection that
 // node 5 in the hand-coded fluid mesh is at the 
 // upstream tip of the cylinder
 Fluid_control_node_pt=Fluid_mesh_pt->node_pt(5);

 // Set error estimator for the fluid mesh
 Z2ErrorEstimator* fluid_error_estimator_pt=new Z2ErrorEstimator;
 fluid_mesh_pt()->spatial_error_estimator_pt()=fluid_error_estimator_pt;

 // Refine uniformly
 Fluid_mesh_pt->refine_uniformly();


 // Build combined global mesh
 //---------------------------

 // Add Solid mesh the problem's collection of submeshes
 add_sub_mesh(solid_mesh_pt());

 // Add traction sub-meshes
 for (unsigned i=0;i<3;i++)
  {
   add_sub_mesh(traction_mesh_pt(i));
  }

 // Add fluid mesh
 add_sub_mesh(fluid_mesh_pt());
 
 // Build combined "global" mesh
 build_global_mesh();
 


 // Apply solid boundary conditons
 //-------------------------------
 
 //Solid mesh: Pin the left boundary (boundary 3) in both directions
 unsigned n_side = mesh_pt()->nboundary_node(3);
 
 //Loop over the nodes
 for(unsigned i=0;i<n_side;i++)
  {
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(0);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(1);
  }
 
 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());
 
 // Apply fluid boundary conditions
 //--------------------------------
 
 //Fluid mesh: Horizontal, traction-free outflow; pinned elsewhere
 unsigned num_bound = fluid_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= fluid_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Parallel, axially traction free outflow at downstream end
     if (ibound != 1)
      {
       fluid_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
       fluid_mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
     else
      {
       fluid_mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
    }
  }//end_of_pin
 
 // Pin redundant pressure dofs in fluid mesh
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(fluid_mesh_pt()->element_pt());


 // Apply boundary conditions for fluid
 //-------------------------------------

 // Impose parabolic flow along boundary 3
 // Current flow rate
 double t=0.0;
 double ampl=Global_Parameters::flux(t);
 unsigned ibound=3; 
 unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   double ycoord = Fluid_mesh_pt->boundary_node_pt(ibound,inod)->x(1); 
   double uy = ampl*6.0*ycoord/Domain_height*(1.0-ycoord/Domain_height);
   Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,uy);
   Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);    
  }
 

 // Complete build of solid elements
 //---------------------------------

 //Pass problem parameters to solid elements
 unsigned  n_element =solid_mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   SOLID_ELEMENT *el_pt = dynamic_cast<SOLID_ELEMENT*>(
    solid_mesh_pt()->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Parameters::Constitutive_law_pt;
   
   //Set the body force
   el_pt->body_force_fct_pt() = Global_Parameters::gravity;

   // Timescale ratio for solid
   el_pt->lambda_sq_pt() = &Global_Parameters::Lambda_sq;
  }
 
 

 // Complete build of fluid elements
 //---------------------------------

 // Set physical parameters in the fluid mesh
 unsigned nelem=fluid_mesh_pt()->nelement();
 for (unsigned e=0;e<nelem;e++)
  {
   // Upcast from GeneralisedElement to the present element
   FLUID_ELEMENT* el_pt = dynamic_cast<FLUID_ELEMENT*>
    (fluid_mesh_pt()->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Parameters::Re;   
   
   //Set the Womersley number
   el_pt->re_st_pt() = &Global_Parameters::ReSt;
   
  }//end_of_loop



 // Setup FSI
 //----------
 
 // Pass Strouhal number to the helper function that automatically applies
 // the no-slip condition
 FSI_functions::Strouhal_for_no_slip=Global_Parameters::St;

 // The velocity of the fluid nodes on the wall (fluid mesh boundary 5,6,7)
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 if (!Global_Parameters::Ignore_fluid_loading)
  {
   for(unsigned ibound=5;ibound<8;ibound++ )
    { 
     unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {   
       Fluid_mesh_pt->boundary_node_pt(ibound, inod)->
        set_auxiliary_node_update_fct_pt(
         FSI_functions::apply_no_slip_on_moving_wall);
      }
    } // done automatic application of no-slip
  
   // Work out which fluid dofs affect the residuals of the wall elements:
   // We pass the boundary between the fluid and solid meshes and 
   // pointers to the meshes. The interaction boundary are boundaries 5,6,7
   // of the 2D fluid mesh.
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,5,Fluid_mesh_pt,Traction_mesh_pt[0]);  
   
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,6,Fluid_mesh_pt,Traction_mesh_pt[2]);  
   
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,7,Fluid_mesh_pt,Traction_mesh_pt[1]); 
  } 


 // Build solver/preconditioner
 //----------------------------

 // Build iterative linear solver
 GMRES<CRDoubleMatrix>* iterative_linear_solver_pt =
  new GMRES<CRDoubleMatrix>;
 
 // Set maximum number of iterations
 iterative_linear_solver_pt->max_iter() = 100;
 
 // Set tolerance
 iterative_linear_solver_pt->tolerance() = 1.0e-10;
 
 // Assign solver
 linear_solver_pt()=iterative_linear_solver_pt;

 // Build preconditioner
 FSIPreconditioner* prec_pt=new FSIPreconditioner(this);

 // Set Navier Stokes mesh:
 prec_pt->set_navier_stokes_mesh(Fluid_mesh_pt);

 // Set solid mesh:
 prec_pt->set_wall_mesh(solid_mesh_pt());

 // Retain fluid onto solid terms
 prec_pt->use_block_triangular_version_with_fluid_on_solid();

 // Set preconditioner
 iterative_linear_solver_pt->preconditioner_pt()= prec_pt;

 // By default, the LSC Preconditioner uses SuperLU as
 // an exact preconditioner (i.e. a solver) for the
 // momentum and Schur complement blocks.
 // Can overwrite this by passing pointers to
 // other preconditioners that perform the (approximate)
 // solves of these blocks.

#ifdef OOMPH_HAS_HYPRE
//Only use HYPRE if we don't have MPI enabled
#ifndef OOMPH_HAS_MPI

 // Create internal preconditioners used on Schur block
 Preconditioner* P_matrix_preconditioner_pt = new HyprePreconditioner;

 HyprePreconditioner* P_hypre_solver_pt =
  static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt);
 
 // Set defaults parameters for use as preconditioner on Poisson-type
 // problem
 Hypre_default_settings::
  set_defaults_for_2D_poisson_problem(P_hypre_solver_pt);

 // Use Hypre for the Schur complement block
 prec_pt->navier_stokes_preconditioner_pt()->
  set_p_preconditioner(P_matrix_preconditioner_pt);

 // Shut up
 P_hypre_solver_pt->disable_doc_time();

#endif
#endif

 // Assign equation numbers
 cout << assign_eqn_numbers() << std::endl; 


}//end_of_constructor



//====start_of_actions_before_newton_convergence_check===================
/// Update the (dependent) fluid node positions following the
/// update of the solid variables
//=======================================================================
template <class FLUID_ELEMENT,class SOLID_ELEMENT> 
void TurekProblem<FLUID_ELEMENT,SOLID_ELEMENT>
::actions_before_newton_convergence_check ()
{
 fluid_mesh_pt()->node_update();
}

 

//===== start_of_actions_before_implicit_timestep=========================
/// Actions before implicit timestep: Update inflow profile
//========================================================================
template <class FLUID_ELEMENT,class SOLID_ELEMENT> 
void TurekProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
actions_before_implicit_timestep()
{
 // Current time
 double t=time_pt()->time();
 
 // Amplitude of flow
 double ampl=Global_Parameters::flux(t);
 
 // Update parabolic flow along boundary 3
 unsigned ibound=3; 
 unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   double ycoord = Fluid_mesh_pt->boundary_node_pt(ibound,inod)->x(1); 
   double uy = ampl*6.0*ycoord/Domain_height*(1.0-ycoord/Domain_height);
   Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,uy);
   Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);    
  }

} //end_of_actions_before_implicit_timestep


//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Re-setup FSI
//========================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT >
void TurekProblem<FLUID_ELEMENT,SOLID_ELEMENT>::actions_after_adapt()
{
 // Unpin all pressure dofs
 RefineableNavierStokesEquations<2>::
  unpin_all_pressure_dofs(fluid_mesh_pt()->element_pt());
 
 // Pin redundant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(fluid_mesh_pt()->element_pt());

 // Unpin all solid pressure dofs
 PVDEquationsBase<2>::
  unpin_all_solid_pressure_dofs(solid_mesh_pt()->element_pt());
 
 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());


 // The velocity of the fluid nodes on the wall (fluid mesh boundary 5,6,7)
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 if (!Global_Parameters::Ignore_fluid_loading)
  {
   for(unsigned ibound=5;ibound<8;ibound++ )
    { 
     unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {   
       Fluid_mesh_pt->boundary_node_pt(ibound, inod)->
        set_auxiliary_node_update_fct_pt(
         FSI_functions::apply_no_slip_on_moving_wall);
      }
    }
 

   // Re-setup the fluid load information for fsi solid traction elements
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,5,Fluid_mesh_pt,Traction_mesh_pt[0]); 
   
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,6,Fluid_mesh_pt,Traction_mesh_pt[2]); 
   
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,7,Fluid_mesh_pt,Traction_mesh_pt[1]); 
  }

 
}// end of actions_after_adapt


 
//============start_of_create_traction_elements==============================
/// Create FSI traction elements 
//=======================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT >
void TurekProblem<FLUID_ELEMENT,SOLID_ELEMENT>::create_fsi_traction_elements()
{

 // Container to collect all nodes in the traction meshes
 std::set<SolidNode*> all_nodes;

 // Traction elements are located on boundaries 0-2:
 for (unsigned b=0;b<3;b++)
  {
   // How many bulk elements are adjacent to boundary b?
   unsigned n_element = solid_mesh_pt()->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     SOLID_ELEMENT* bulk_elem_pt = dynamic_cast<SOLID_ELEMENT*>(
      solid_mesh_pt()->boundary_element_pt(b,e));
     
     //What is the index of the face of the element e along boundary b
     int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);
          
     // Create new element and add to mesh
     Traction_mesh_pt[b]->add_element_pt(
      new FSISolidTractionElement<SOLID_ELEMENT,2>(bulk_elem_pt,face_index));
     
    } //end of loop over bulk elements adjacent to boundary b

   // Identify unique nodes
   unsigned nnod=solid_mesh_pt()->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     all_nodes.insert(solid_mesh_pt()->boundary_node_pt(b,j));
    }
  }

 // Build combined mesh of fsi traction elements
 Combined_traction_mesh_pt=new SolidMesh(Traction_mesh_pt);
 
 // Stick nodes into combined traction mesh
 for (std::set<SolidNode*>::iterator it=all_nodes.begin();
      it!=all_nodes.end();it++)
  {
   Combined_traction_mesh_pt->add_node_pt(*it);
  }

} // end of create_traction_elements




//=====start_of_doc_solution========================================
/// Doc the solution
//==================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT >
void TurekProblem<FLUID_ELEMENT,SOLID_ELEMENT>::doc_solution(
 DocInfo& doc_info, ofstream& trace_file)
{
 
//  FSI_functions::doc_fsi<AlgebraicNode>(Fluid_mesh_pt,
//                                        Combined_traction_mesh_pt,
//                                        doc_info);

//  pause("done");

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned n_plot = 5; 

 // Output solid solution
 sprintf(filename,"%s/solid_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 solid_mesh_pt()->output(some_file,n_plot);
 some_file.close();
 
 // Output fluid solution
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 fluid_mesh_pt()->output(some_file,n_plot);
 some_file.close();


//Output the traction
 sprintf(filename,"%s/traction%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
// Loop over the traction meshes
 for(unsigned i=0;i<3;i++)
  {
   // Loop over the element in traction_mesh_pt
   unsigned n_element = Traction_mesh_pt[i]->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     FSISolidTractionElement<SOLID_ELEMENT,2>* el_pt = 
      dynamic_cast<FSISolidTractionElement<SOLID_ELEMENT,2>* > (
       Traction_mesh_pt[i]->element_pt(e) );
     
     el_pt->output(some_file,5);
    }
  }
 some_file.close();  


 // Write trace (we're only using Taylor Hood elements so we know that
 // the pressure is the third value at the fluid control node...
 trace_file << time_pt()->time() << " " 
            << Solid_control_node_pt->x(0) << " " 
            << Solid_control_node_pt->x(1) << " " 
            << Fluid_control_node_pt->value(2) << " "
            << Global_Parameters::flux(time_pt()->time()) << " " 
            << std::endl;
 
 cout << "Doced solution for step " 
      << doc_info.number() 
      << std::endl << std::endl << std::endl;

}//end_of_doc_solution



//=======start_of_main==================================================
/// Driver 
//======================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Get case id as string
 string case_id="FSI1";
 if (CommandLineArgs::Argc==1)
  {
   oomph_info << "No command line arguments; running self-test FSI1" 
              << std::endl;
  }
 else if (CommandLineArgs::Argc==2)
  {
   case_id=CommandLineArgs::Argv[1];
  }
 else
  {
   oomph_info << "Wrong number of command line arguments" << std::endl;
   oomph_info << "Enter none (for default) or one (namely the case id" 
              << std::endl;
   oomph_info << "which should be one of: FSI1, FSI2, FSI3, CSM1" 
              << std::endl;
  } 
 std::cout << "Running case " << case_id << std::endl;

 // Setup parameters for case identified by command line
 // argument
 Global_Parameters::set_parameters(case_id);

 // Prepare output
 DocInfo doc_info;
 ofstream trace_file; 
 doc_info.set_directory("RESLT");
 trace_file.open("RESLT/trace.dat");
 
 // Length and height of domain
 double length=25.0;
 double height=4.1;

 //Set up the problem
 TurekProblem<AlgebraicElement<RefineableQTaylorHoodElement<2> >,
  RefineableQPVDElement<2,3> > problem(length, height);

 // Default number of timesteps
 unsigned nstep=4000;
 if (Global_Parameters::Case_ID=="FSI1")
  {
   std::cout << "Reducing number of steps for FSI1 " << std::endl;
   nstep=400;
  }

 if (CommandLineArgs::Argc==1)
  {
   std::cout << "Reducing number of steps for validation " << std::endl;
   nstep=2;
  }

 //Timestep: 
 double dt=Global_Parameters::Dt;

 // Initialise timestep 
 problem.initialise_dt(dt);

 // Impulsive start
 problem.assign_initial_values_impulsive(dt);

 // Doc the initial condition
 problem.doc_solution(doc_info,trace_file);
 doc_info.number()++; 
 
 // Don't re-set the initial conditions when adapting during first
 // timestep
 bool first = false;
 
 // Max number of adaptation for time-stepping
 unsigned max_adapt=1;
 
 for(unsigned i=0;i<nstep;i++)
  { 
  // Solve the problem
  problem.unsteady_newton_solve(dt,max_adapt,first); 
  
  // Output the solution
  problem.doc_solution(doc_info,trace_file);
  
  // Step number
  doc_info.number()++;
 }
 
 trace_file.close(); 

}//end of main








