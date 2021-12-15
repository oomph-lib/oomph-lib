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

 /// Actions before distribute: Remove traction elements
 void actions_before_distribute();

 /// Actions after distribute: Add traction elements, 
 /// re-setup the fsi lookup scheme
 void actions_after_distribute();

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

 /// Delete FSI traction elements
 void delete_fsi_traction_elements();

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

 /// Backed-up x coordinate of fluid control node
 double Fluid_control_x0;

 /// Backed-up y coordinate of fluid control node
 double Fluid_control_x1;

 /// Boolean indicating if the current processor contains the
 /// fluid control node
 bool I_have_the_fluid_control_node;
 
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

 // By default all processors contain the fluid control node
 I_have_the_fluid_control_node=true;
 
 
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

 // Complete build of solid elements - needs to happen before refinement
 //---------------------------------------------------------------------

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
  } // end build of solid elements

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
  new MeshAsGeomObject(Traction_mesh_pt[0]);
 
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

 std::cout << "Coordinates of fluid control point " 
           << Fluid_control_node_pt->x(0) << " " 
           << Fluid_control_node_pt->x(1) << " " << std::endl;

 // Back it up so we can check if the node still exists
 // once the problem has been distributed
 Fluid_control_x0=Fluid_control_node_pt->x(0);
 Fluid_control_x1=Fluid_control_node_pt->x(1);

 // Set error estimator for the fluid mesh
 Z2ErrorEstimator* fluid_error_estimator_pt=new Z2ErrorEstimator;
 fluid_mesh_pt()->spatial_error_estimator_pt()=fluid_error_estimator_pt;

 // Refine uniformly
 Fluid_mesh_pt->refine_uniformly();


 // Build combined global mesh
 //---------------------------

 // Add Fluid mesh. Note: It's important that the fluid mesh is
 // added before the solid mesh!
 add_sub_mesh(fluid_mesh_pt());
 
 // Now add the solid mesh to the problem's collection of submeshes
 add_sub_mesh(solid_mesh_pt());

 // Add traction sub-meshes
 for (unsigned i=0;i<3;i++)
  {
   add_sub_mesh(traction_mesh_pt(i));
  }

 // Build combined "global" mesh
 build_global_mesh();
 


 // Apply solid boundary conditions
 //--------------------------------
 
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

 // If the solid is to be loaded by the fluid, then set up the interaction
 // and specify the velocity of the fluid nodes based on the wall motion
 if (!Global_Parameters::Ignore_fluid_loading)
  {

#ifdef OLD_FSI

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

#else

   // Package fsi solid traction meshes and boundary IDs in 
   // fluid mesh
   Vector<unsigned> fluid_fsi_boundary_id(3);
   Vector<Mesh*> traction_mesh_pt(3);
   fluid_fsi_boundary_id[0]=5;
   traction_mesh_pt[0]=Traction_mesh_pt[0];
   fluid_fsi_boundary_id[1]=6;
   traction_mesh_pt[1]=Traction_mesh_pt[2];
   fluid_fsi_boundary_id[2]=7;
   traction_mesh_pt[2]=Traction_mesh_pt[1];
   
   // Vector based FSI setup
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,fluid_fsi_boundary_id,Fluid_mesh_pt,
     traction_mesh_pt);

#endif


   // The velocity of the fluid nodes on the wall (fluid mesh boundary 5,6,7)
   // is set by the wall motion -- hence the no-slip condition must be
   // re-applied whenever a node update is performed for these nodes. 
   // Such tasks may be performed automatically by the auxiliary node update 
   // function specified by a function pointer:
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
  } // end of FSI setup

 // Use SuperLU_dist as the solver
 linear_solver_pt() = new SuperLUSolver;
 static_cast<SuperLUSolver*>(linear_solver_pt())
  ->set_solver_type(SuperLUSolver::Distributed);
 static_cast<SuperLUSolver*>(linear_solver_pt())
  ->use_distributed_solve_in_superlu_dist();

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
///  Actions after adapt: Re-setup traction elements and FSI
//========================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT >
void TurekProblem<FLUID_ELEMENT,SOLID_ELEMENT>::actions_after_adapt()
{
 // If the solid mesh has been allowed to refine, then we need to
 // delete the previous solid traction elements and create new ones,
 // as the traction meshes themselves are not refineable
 delete_fsi_traction_elements();

 // (Re-)build the FSI traction elements
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
 // that determine particular boundaries of the fluid mesh
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

 // Tell the fluid mesh about the new "refined" MeshAsGeomObjects
 delete fluid_mesh_pt()->bottom_flag_pt();
 fluid_mesh_pt()->set_bottom_flag_pt(bottom_flag_pt);
 delete fluid_mesh_pt()->top_flag_pt();
 fluid_mesh_pt()->set_top_flag_pt(top_flag_pt);
 delete fluid_mesh_pt()->tip_flag_pt();
 fluid_mesh_pt()->set_tip_flag_pt(tip_flag_pt);

 // Call update_node_update for all the fluid mesh nodes, as the compound
 // geometric objects representing the interaction boundaries have changed
 unsigned n_fluid_node=fluid_mesh_pt()->nnode();
 for (unsigned n=0;n<n_fluid_node;n++)
  {
   // Get the (algebraic) node
   AlgebraicNode* alg_nod_pt=dynamic_cast<AlgebraicNode*>
    (fluid_mesh_pt()->node_pt(n));

   // Call update_node_update for this node
   fluid_mesh_pt()->update_node_update(alg_nod_pt);
  }

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

 // Now rebuild the global mesh
 rebuild_global_mesh();

 // If the solid is to be loaded by the fluid, then set up the interaction
 // and specify the velocity of the fluid nodes based on the wall motion
 if (!Global_Parameters::Ignore_fluid_loading)
  {

#ifdef OLD_FSI

   // Re-setup the fluid load information for fsi solid traction elements
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,5,Fluid_mesh_pt,Traction_mesh_pt[0]); 

   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,6,Fluid_mesh_pt,Traction_mesh_pt[2]); 

   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,7,Fluid_mesh_pt,Traction_mesh_pt[1]); 

#else

   // Package fsi solid traction meshes and boundary IDs in 
   // fluid mesh
   Vector<unsigned> fluid_fsi_boundary_id(3);
   Vector<Mesh*> traction_mesh_pt(3);
   fluid_fsi_boundary_id[0]=5;
   traction_mesh_pt[0]=Traction_mesh_pt[0];
   fluid_fsi_boundary_id[1]=6;
   traction_mesh_pt[1]=Traction_mesh_pt[2];
   fluid_fsi_boundary_id[2]=7;
   traction_mesh_pt[2]=Traction_mesh_pt[1];
   
   // Vector based FSI setup
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,fluid_fsi_boundary_id,Fluid_mesh_pt,
     traction_mesh_pt);

#endif

   // The velocity of the fluid nodes on the wall (fluid mesh boundary 5,6,7)
   // is set by the wall motion -- hence the no-slip condition must be
   // re-applied whenever a node update is performed for these nodes. 
   // Such tasks may be performed automatically by the auxiliary node update 
   // function specified by a function pointer:
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

  }

}// end of actions_after_adapt

//==================start_of_actions_before_distribute====================
/// Actions before distribute: Make sure that the bulk solid elements 
/// attached to the FSISolidTractionElements are kept as halo elements.
/// Unlike in most other parallel codes we DON'T delete the 
/// FSISolidTractionElements here, though, because they need to 
/// be around while the fluid mesh is adapted.
//========================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT >
void TurekProblem<FLUID_ELEMENT,SOLID_ELEMENT>::actions_before_distribute()
{
 // The bulk elements attached to the traction elements need to be kept
 // as halo elements

 // There are 3 traction meshes
 for (unsigned b=0;b<3;b++)
  {
   // Loop over elements in traction meshes
   unsigned n_element=Traction_mesh_pt[b]->nelement();
   for (unsigned e=0;e<n_element;e++)
    {
     FSISolidTractionElement<SOLID_ELEMENT,2>* traction_elem_pt=
      dynamic_cast<FSISolidTractionElement<SOLID_ELEMENT,2>* >
      (Traction_mesh_pt[b]->element_pt(e));

     // Get the bulk element (which is a SOLID_ELEMENT)
     SOLID_ELEMENT* solid_elem_pt = dynamic_cast<SOLID_ELEMENT*>
      (traction_elem_pt->bulk_element_pt());

     // Require bulk to be kept as a (possible) halo element
     // Note: The traction element itself will "become" a halo element 
     // when it is recreated after the distribution has taken place
     solid_elem_pt->set_must_be_kept_as_halo();
    }
  } // end of loop over meshes of fsi traction elements
 

 // Flush all the submeshes out but keep the meshes of FSISolidTractionElements
 // alive (i.e. don't delete them)
 flush_sub_meshes();

 // Add the fluid mesh and the solid mesh back again
 // Remember that it's important that the fluid mesh is
 // added before the solid mesh!
 add_sub_mesh(fluid_mesh_pt());
 add_sub_mesh(solid_mesh_pt());

 // Rebuild global mesh
 rebuild_global_mesh();

} // end of actions before distribute


//==================start_of_actions_after_distribute=====================
///  Actions after distribute: Re-setup FSI
//========================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT >
void TurekProblem<FLUID_ELEMENT,SOLID_ELEMENT>::actions_after_distribute()
{
 // The solid mesh has now been distributed, so it now has halo elements
 // on certain processors. The traction elements attached to these new
 // halo elements need to be halo themselves, so we need to delete the
 // old ones and re-attach new ones. Recall that FaceElements attached
 // to bulk halo elements become halos themselves.
 delete_fsi_traction_elements();

 // (Re-)Build the FSI traction elements
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
 // that determine particular boundaries of the fluid mesh
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


 // Delete the old MeshAsGeomObjects and tell the fluid mesh 
 // about the new ones.
 delete fluid_mesh_pt()->bottom_flag_pt();
 fluid_mesh_pt()->set_bottom_flag_pt(bottom_flag_pt);
 delete fluid_mesh_pt()->top_flag_pt();
 fluid_mesh_pt()->set_top_flag_pt(top_flag_pt);
 delete fluid_mesh_pt()->tip_flag_pt();
 fluid_mesh_pt()->set_tip_flag_pt(tip_flag_pt);

 // Call update_node_update for all the fluid mesh nodes, as the
 // geometric objects representing the fluid mesh boundaries have changed
 unsigned n_fluid_node=fluid_mesh_pt()->nnode();
 for (unsigned n=0;n<n_fluid_node;n++)
  {
   // Get the (algebraic) node
   AlgebraicNode* alg_nod_pt=dynamic_cast<AlgebraicNode*>
    (fluid_mesh_pt()->node_pt(n));

   // Call update_node_update for this node
   fluid_mesh_pt()->update_node_update(alg_nod_pt);
  }

 // Add the traction meshes back to the problem
 for (unsigned i=0;i<3;i++)
  {
   add_sub_mesh(traction_mesh_pt(i));
  }

 // Rebuild global mesh
 rebuild_global_mesh();

 // If the solid is to be loaded by the fluid, then set up the interaction
 // and specify the velocity of the fluid nodes based on the wall motion
 if (!Global_Parameters::Ignore_fluid_loading)
  {

#ifdef OLD_FSI

   // Re-setup the fluid load information for fsi solid traction elements
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,5,Fluid_mesh_pt,Traction_mesh_pt[0]); 

   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,6,Fluid_mesh_pt,Traction_mesh_pt[2]); 

   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,7,Fluid_mesh_pt,Traction_mesh_pt[1]); 

#else

   // Package fsi solid traction meshes and boundary IDs in 
   // fluid mesh
   Vector<unsigned> fluid_fsi_boundary_id(3);
   Vector<Mesh*> traction_mesh_pt(3);
   fluid_fsi_boundary_id[0]=5;
   traction_mesh_pt[0]=Traction_mesh_pt[0];
   fluid_fsi_boundary_id[1]=6;
   traction_mesh_pt[1]=Traction_mesh_pt[2];
   fluid_fsi_boundary_id[2]=7;
   traction_mesh_pt[2]=Traction_mesh_pt[1];
   
   // Vector based FSI setup
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,fluid_fsi_boundary_id,Fluid_mesh_pt,
     traction_mesh_pt);

#endif

   // The velocity of the fluid nodes on the wall (fluid mesh boundary 5,6,7)
   // is set by the wall motion -- hence the no-slip condition must be
   // re-applied whenever a node update is performed for these nodes. 
   // Such tasks may be performed automatically by the auxiliary node update 
   // function specified by a function pointer:
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

  } // end of (re-)assignment of the auxiliary node update fct

 // Re-set control nodes

 // Loop over fluid nodes
 unsigned n_fluid_nod=fluid_mesh_pt()->nnode();
 for (unsigned j=0;j<n_fluid_nod;j++)
  {
   if ((fluid_mesh_pt()->node_pt(j)->x(0)==Fluid_control_x0) &&
       (fluid_mesh_pt()->node_pt(j)->x(1)==Fluid_control_x1))
    {
     // The node still exists on this process...
     I_have_the_fluid_control_node=true;
     break;
    }

   // if the end has been reached then we've lost the control node
   if (j==n_fluid_nod-1)
    {
     I_have_the_fluid_control_node=false;
    }
  }
 
} // end of actions after distribute

 
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


//============start_of_delete_traction_elements==========================
/// Delete FSI traction elements 
//=======================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT >
void TurekProblem<FLUID_ELEMENT,SOLID_ELEMENT>::delete_fsi_traction_elements()
{
 // There are 3 traction meshes
 for (unsigned b=0;b<3;b++)
  {
   unsigned n_element=Traction_mesh_pt[b]->nelement();
   for (unsigned e=0;e<n_element;e++)
    {
     // Kill the element
     delete Traction_mesh_pt[b]->element_pt(e);
    }

   // Wipe the mesh
   Traction_mesh_pt[b]->flush_element_and_node_storage();
  }
} // end of delete traction elements


//=====start_of_doc_solution========================================
/// Doc the solution
//==================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT >
void TurekProblem<FLUID_ELEMENT,SOLID_ELEMENT>::doc_solution(
 DocInfo& doc_info, ofstream& trace_file)
{
 // Current process
 int my_rank=this->communicator_pt()->my_rank();

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned n_plot = 5; 

 // Output solid solution
 sprintf(filename,"%s/solid_soln%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),my_rank);
 some_file.open(filename);
 solid_mesh_pt()->output(some_file,n_plot);
 some_file.close();
 
 // Output fluid solution
 sprintf(filename,"%s/soln%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),my_rank);
 some_file.open(filename);
 fluid_mesh_pt()->output(some_file,n_plot);
 some_file.close();


//Output the traction
 sprintf(filename,"%s/traction%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),my_rank);
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

     // Only output halo elements
     if (!el_pt->is_halo())
      {
       el_pt->output(some_file,5);
      }
    }
  }
 some_file.close();  


 // Write trace if we still have the fluid control node
 // (the solid control node
 if (I_have_the_fluid_control_node)
  {
   trace_file << time_pt()->time() << " " 
              << Solid_control_node_pt->x(0) << " " 
              << Solid_control_node_pt->x(1) << " " 
              << Fluid_control_node_pt->value(2) << " "
              << Global_Parameters::flux(time_pt()->time()) << " " 
              << std::endl;
  }
 
 oomph_info << "Doced solution for step " 
            << doc_info.number()
            << std::endl << std::endl << std::endl;

}//end_of_doc_solution



//=======start_of_main==================================================
/// Driver 
//======================================================================
int main(int argc, char* argv[])
{

#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc,argv);
#endif

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

 // Length and height of domain
 double length=25.0;
 double height=4.1;

 //Set up the problem
 TurekProblem<AlgebraicElement<RefineableQTaylorHoodElement<2> >,
  RefineableQPVDElement<2,3> > problem(length, height);

 // Prepare output
 DocInfo doc_info;
 ofstream trace_file; 
 char filename[100];
 doc_info.set_directory("RESLT_TUREK");
 sprintf(filename,"%s/trace_on_proc%i.dat",doc_info.directory().c_str(),
         problem.communicator_pt()->my_rank());
 trace_file.open(filename);
 
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

#ifdef OOMPH_HAS_MPI
 // Distribute the problem
 bool report_stats=true;

 // Are there command-line arguments?
 if (CommandLineArgs::Argc==1)
  {
   // No arguments, so it's a validation run
   std::ifstream input_file;

   // All meshes are partitioned
   unsigned n_partition=problem.mesh_pt()->nelement();
   Vector<unsigned> in_element_partition(n_partition,0);

   // Get partition from file
   sprintf(filename,"turek_flag_partition.dat");
   input_file.open(filename);
   std::string input_string;
   for (unsigned e=0;e<n_partition;e++)
    {
     getline(input_file,input_string,'\n');
     in_element_partition[e]=atoi(input_string.c_str());
    }
  
   // Now distribute the (still uniformly refined) problem
   problem.distribute(in_element_partition,report_stats);
  }
 else
  {
   // There were command-line arguments, so distribute without
   // reference to any partition vector
   problem.distribute(report_stats);
  }


#endif
 
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

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

}//end of main








