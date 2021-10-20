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
// Driver to demonstrate the interaction of a fluid flow with a rigid body.
// The driver solves a classic problem of the motion of an rigid ellipse in
// a shear flow. The problem imposes a smooth, but rapid, start-up from
// no flow to a shear imposed on the boundaries of a rectangular domain
// and the ellipse soon settles into approximate Jeffery orbits for small
// enough Reynolds number.

//Generic routines
#include "generic.h"

// The equations
#include "navier_stokes.h"
#include "solid.h"
#include "constitutive.h"
#include "rigid_body.h"

// The mesh
#include "meshes/triangle_mesh.h"

//Thin wrapper to "custom" TaylorHoodElements that overload output functions
#include "my_taylor_hood_elements.h"

#include <algorithm>

using namespace std;
using namespace oomph;

//==start_of_namespace==============================
/// Namespace for Problem Parameters
//==================================================
 namespace Problem_Parameter
 {    
  /// Reynolds number
  double Re=1.0;

  /// Strouhal number
  double St = 1.0;
 
  ///  Density ratio (Solid density / Fluid density)
  double Density_ratio = 1.0;

  /// Initial axis of the elliptical solid in x-direction
  double A = 0.25;
  
  /// Initial axis of the elliptical solid in y-direction
  /// (N.B. 2B = 1 is the reference length scale)
  double B = 0.5;

  /// Pseudo-solid (mesh) Poisson ratio
  double Nu=0.3;

  ///  Pseudo-solid (mesh) "density" 
  /// Set to zero because we don't want inertia in the node update!
  double Lambda_sq=0.0;

  /// Constitutive law used to determine the mesh deformation
  ConstitutiveLaw *Constitutive_law_pt=   
   new GeneralisedHookean(&Problem_Parameter::Nu);

 } // end_of_namespace


//=======================================================================
///Exact solution for the rotation of an ellipse in unbounded shear flow
///as computed by Jeffery (1922)
//=======================================================================
namespace Jeffery_Solution
{
 ///Null function
 double null(const double &t) {return 0.0;}

 ///Angular position as a function of time t
 double angle(const double &t)
 {
  const double a = Problem_Parameter::A;
  const double b = Problem_Parameter::B;
  
  return atan((b/a)*tan((a*b*t)/(b*b + a*a)));
 }

 ///Angular velocity as function of time t
 double velocity(const double &t)
 {
  const double a = Problem_Parameter::A;
  const double b = Problem_Parameter::B;

  //Get the angle
  double chi = angle(t);

  //Now return the velocity
  return (a*a*sin(chi)*sin(chi) + b*b*cos(chi)*cos(chi))/(a*a + b*b);
 }

 ///Angular acceleration as a function of time t (should always be zero)
 double acceleration(const double &t)
 {
  const double a = Problem_Parameter::A;
  const double b = Problem_Parameter::A;

  //Get the angle and velocity
  double chi = angle(t);
  double chi_dot = velocity(t);

  //Now return the acceleration
  return 2.0*(a*a - b*b)*(sin(chi)*cos(chi))*chi_dot/(a*a + b*b);
 }
}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//===================start_of_general_ellipse=================================
///  A geometric object for an ellipse with initial centre of mass at
/// (centre_x, centre_y) with axis in the x direction given by 2a
/// and in the y-direction given by 2b. The boundary of the ellipse is
/// parametrised by its angle.
//============================================================================
class GeneralEllipse : public GeomObject
{
private:
 
 //Storage for the centre of mass and semi-major and semi-minor axes
 double Centre_x, Centre_y, A, B;

public:
 
 ///  Simple Constructor that transfers appropriate geometric 
 /// parameters into internal data
 GeneralEllipse(const double &centre_x, const double &centre_y,
                const double &a, const double &b)
  : GeomObject(1,2), Centre_x(centre_x), Centre_y(centre_y), A(a), B(b)
  {}

 /// Empty Destructor
 ~GeneralEllipse() {}

 ///Return the position of the ellipse boundary as a function of 
 ///the angle xi[0]
 void position(const Vector<double> &xi, Vector<double> &r) const
  {
   r[0] = Centre_x + A*cos(xi[0]);
   r[1] = Centre_y + B*sin(xi[0]);
  }

 //Return the position which is always fixed
 void position(const unsigned &t,
               const Vector<double> &xi, Vector<double> &r) const
  {
   return position(xi,r);
  }

};
//end_of_general_ellipse


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Unstructured Navier-Stokes ALE Problem for a rigid ellipse 
/// immersed within a viscous fluid
//====================================================================
template<class ELEMENT>
class UnstructuredImmersedEllipseProblem : public Problem
{

public:

 /// Constructor
 UnstructuredImmersedEllipseProblem();
 
 /// Destructor
 ~UnstructuredImmersedEllipseProblem();

 /// Reset the boundary conditions when timestepping
 void actions_before_implicit_timestep()
  {
   this->set_boundary_velocity();
  }
 
 /// Wipe the meshes of Lagrange multiplier and drag elements
 void actions_before_adapt();
 
 /// Rebuild the meshes of Lagrange multiplier and drag elements
 void actions_after_adapt();
 
 ///  Re-apply the no slip condition (imposed indirectly via dependent
 /// velocities)
 void actions_before_newton_convergence_check()
  {
   // Update mesh -- this applies the auxiliary node update function
   Fluid_mesh_pt->node_update();
  }
 
 ///  Set boundary condition, assign auxiliary node update fct.
 /// Complete the build of all elements, attach power elements that allow
 /// computation of drag vector
 void complete_problem_setup();

 ///Set the boundary velocity
 void set_boundary_velocity();

 /// Function that solves a simplified problem to ensure that 
 ///the positions of the boundary nodes are initially consistent with
 ///the lagrange multiplier formulation
 void solve_for_consistent_nodal_positions();

 /// Doc the solution
 void doc_solution(const bool& project=false);
  
 /// Output the exact solution
 void output_exact_solution(std::ofstream &output_file);

private:
 
 ///  Create elements that enforce prescribed boundary motion
 /// for the pseudo-solid fluid mesh by Lagrange multipliers
 void create_lagrange_multiplier_elements();

 ///  Delete elements that impose the prescribed boundary displacement
 /// and wipe the associated mesh
 void delete_lagrange_multiplier_elements();

 ///  Create elements that calculate the drag and torque on
 /// the boundaries
 void create_drag_elements();
 
 ///  Delete elements that calculate the drag and torque on the 
 /// boundaries
 void delete_drag_elements();

 ///Pin the degrees of freedom associated with the solid bodies
 void pin_rigid_body();

 ///Unpin the degrees of freedom associated with the solid bodies
 void unpin_rigid_body();
 
 /// Pointers to mesh of Lagrange multiplier elements
 SolidMesh* Lagrange_multiplier_mesh_pt;

 /// Pointer to Fluid_mesh
 RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;
 
 /// Triangle mesh polygon for outer boundary 
 TriangleMeshPolygon* Outer_boundary_polygon_pt;

 /// Mesh of drag elements
 Vector<Mesh*> Drag_mesh_pt;

 /// Mesh of the generalised elements for the rigid bodies
 Mesh* Rigid_body_mesh_pt;

 /// Storage for the geom object
 Vector<GeomObject*> Rigid_body_pt;
 
 /// Internal DocInfo object
 DocInfo Doc_info;
 
 /// File to document the norm of the solution (for validation purposes)
 ofstream Norm_file;
 
 /// File to document the motion of the centre of gravity
 ofstream Cog_file;
 
 /// File to document the exact motion of the centre of gravity
 ofstream Cog_exact_file;

}; // end_of_problem_class


//==start_constructor=====================================================
/// Constructor: Open output files, construct time steppers, build 
///  fluid mesh, immersed rigid body and combine to form the problem
//========================================================================
template<class ELEMENT>
UnstructuredImmersedEllipseProblem<ELEMENT>::
UnstructuredImmersedEllipseProblem()
{ 
 // Output directory
 this->Doc_info.set_directory("RESLT");

 // Open norm file
 this->Norm_file.open("RESLT/norm.dat");

 // Open file to trace the centre of gravity
 this->Cog_file.open("RESLT/cog_trace.dat");

 // Open file to document the exact motion of the centre of gravity
 this->Cog_exact_file.open("RESLT/cog_exact_trace.dat");

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 this->add_time_stepper_pt(new BDF<2>);

 // Allocate a timestepper for the rigid body
 this->add_time_stepper_pt(new Newmark<2>);

 // Define the boundaries: Polyline with 4 different
 // boundaries for the outer boundary and 1 internal elliptical hole
 
 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separate polyline segments
 //---------------------------------
 Vector<TriangleMeshCurveSection*> boundary_segment_pt(4);

 //Set the length of the channel
 double half_length = 5.0;
 double half_height = 2.5;
 
 // Initialize boundary segment
 Vector<Vector<double> > bound_seg(2);
 for(unsigned i=0;i<2;i++) {bound_seg[i].resize(2);}

 // First boundary segment
 bound_seg[0][0]=-half_length;
 bound_seg[0][1]=-half_height;
 bound_seg[1][0]=-half_length;
 bound_seg[1][1]=half_height;
 
 // Specify 1st boundary id
 unsigned bound_id = 0;

 // Build the 1st boundary segment
 boundary_segment_pt[0] = new TriangleMeshPolyLine(bound_seg,bound_id);
 
 // Second boundary segment
 bound_seg[0][0]=-half_length;
 bound_seg[0][1]=half_height;
 bound_seg[1][0]=half_length;
 bound_seg[1][1]=half_height;

 // Specify 2nd boundary id
 bound_id = 1;

 // Build the 2nd boundary segment
 boundary_segment_pt[1] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Third boundary segment
 bound_seg[0][0]=half_length;
 bound_seg[0][1]=half_height;
 bound_seg[1][0]=half_length;
 bound_seg[1][1]=-half_height;

 // Specify 3rd boundary id
 bound_id = 2;

 // Build the 3rd boundary segment
 boundary_segment_pt[2] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Fourth boundary segment
 bound_seg[0][0]=half_length;
 bound_seg[0][1]=-half_height;
 bound_seg[1][0]=-half_length;
 bound_seg[1][1]=-half_height;

 // Specify 4th boundary id
 bound_id = 3;

 // Build the 4th boundary segment
 boundary_segment_pt[3] = new TriangleMeshPolyLine(bound_seg,bound_id);
  
 // Create the triangle mesh polygon for outer boundary using boundary segment
 Outer_boundary_polygon_pt = new TriangleMeshPolygon(boundary_segment_pt);

 // Now build the moving rigid body
 //-------------------------------------

 // We have one rigid body
 Rigid_body_pt.resize(1);
 Vector<TriangleMeshClosedCurve*> hole_pt(1);

 // Build Rigid Body
 //-----------------
 double x_center = 0.0;
 double y_center = 0.0;
 double A = Problem_Parameter::A;
 double B = Problem_Parameter::B;
 GeomObject* temp_hole_pt = new GeneralEllipse(x_center,y_center,A,B);
 Rigid_body_pt[0] = new ImmersedRigidBodyElement(temp_hole_pt,
                                                 this->time_stepper_pt(1));

 // Build the two parts of the curvilinear boundary from the rigid body
 Vector<TriangleMeshCurveSection*> curvilinear_boundary_pt(2);

 //First section (boundary 4)
 double zeta_start=0.0;
 double zeta_end=MathematicalConstants::Pi;
 unsigned nsegment=8; 
 unsigned boundary_id=4; 
 curvilinear_boundary_pt[0]=new TriangleMeshCurviLine(
  Rigid_body_pt[0],zeta_start,zeta_end,nsegment,boundary_id);

 //Second section (boundary 5)
 zeta_start=MathematicalConstants::Pi;
 zeta_end=2.0*MathematicalConstants::Pi;
 nsegment=8; 
 boundary_id=5; 
 curvilinear_boundary_pt[1]=new TriangleMeshCurviLine(
  Rigid_body_pt[0],zeta_start,zeta_end, 
  nsegment,boundary_id);
  
 // Combine to form a hole in the fluid mesh
 Vector<double> hole_coords(2);
 hole_coords[0]=0.0;
 hole_coords[1]=0.0;
 Vector<TriangleMeshClosedCurve*> curvilinear_hole_pt(1);
 hole_pt[0]=
  new TriangleMeshClosedCurve(
   curvilinear_boundary_pt,hole_coords);
 
 // Now build the mesh, based on the boundaries specified by
 //---------------------------------------------------------
 // polygons just created
 //----------------------

 TriangleMeshClosedCurve* closed_curve_pt=Outer_boundary_polygon_pt;

 double uniform_element_area=1.0;

 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters(
   closed_curve_pt);

 // Define the holes on the domain
 triangle_mesh_parameters.internal_closed_curve_pt() =
   hole_pt;

 // Define the maximum element area
 triangle_mesh_parameters.element_area() =
   uniform_element_area;

 // Create the mesh
 Fluid_mesh_pt =
   new RefineableSolidTriangleMesh<ELEMENT>(
     triangle_mesh_parameters, this->time_stepper_pt());

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Fluid_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 Fluid_mesh_pt->max_permitted_error()=0.005;
 Fluid_mesh_pt->min_permitted_error()=0.001; 
 Fluid_mesh_pt->max_element_size()=1.0;
 Fluid_mesh_pt->min_element_size()=0.001; 

 // Use coarser mesh during validation
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   Fluid_mesh_pt->min_element_size()=0.01; 
  }

 // Set boundary condition, assign auxiliary node update fct,
 // complete the build of all elements, attach power elements that allow
 // computation of drag vector
 complete_problem_setup();
  
 //Set the parameters of the rigid body elements
 ImmersedRigidBodyElement* rigid_element1_pt = 
  dynamic_cast<ImmersedRigidBodyElement*>(Rigid_body_pt[0]);
 rigid_element1_pt->initial_centre_of_mass(0) = x_center;
 rigid_element1_pt->initial_centre_of_mass(1) = y_center; 
 rigid_element1_pt->mass_shape() = MathematicalConstants::Pi*A*B;
 rigid_element1_pt->moment_of_inertia_shape() = 
 0.25*MathematicalConstants::Pi*A*B*(A*A + B*B);
 rigid_element1_pt->re_pt() = &Problem_Parameter::Re;
 rigid_element1_pt->st_pt() = &Problem_Parameter::St;
 rigid_element1_pt->density_ratio_pt() = &Problem_Parameter::Density_ratio;

 //Pin the position of the centre of mass
 rigid_element1_pt->pin_centre_of_mass_coordinate(0);
 rigid_element1_pt->pin_centre_of_mass_coordinate(1);

  // Create the mesh for the rigid bodies
 Rigid_body_mesh_pt = new Mesh;
 Rigid_body_mesh_pt->add_element_pt(rigid_element1_pt);

 // Create the drag mesh for the rigid bodies
 Drag_mesh_pt.resize(1);
 for(unsigned m=0;m<1;m++) {Drag_mesh_pt[m] = new Mesh;}
 this->create_drag_elements();

 //Add the drag mesh to the appropriate rigid bodies
 rigid_element1_pt->set_drag_mesh(Drag_mesh_pt[0]);



 // Create Lagrange multiplier mesh for boundary motion
 //----------------------------------------------------
 // Construct the mesh of elements that enforce prescribed boundary motion
 // of pseudo-solid fluid mesh by Lagrange multipliers
 Lagrange_multiplier_mesh_pt=new SolidMesh;
 create_lagrange_multiplier_elements();


 // Combine meshes
 //---------------
 
 // Add Fluid_mesh_pt sub meshes
 this->add_sub_mesh(Fluid_mesh_pt);

 // Add Lagrange_multiplier sub meshes
 this->add_sub_mesh(this->Lagrange_multiplier_mesh_pt);

 this->add_sub_mesh(this->Rigid_body_mesh_pt);
 
 // Build global mesh
 this->build_global_mesh();
    
 // Setup equation numbering scheme
 cout <<"Number of equations: " << this->assign_eqn_numbers() << std::endl;
 
} // end_of_constructor


//========================================================================
/// Destructor that cleans up memory and closes files
//========================================================================
template<class ELEMENT>
UnstructuredImmersedEllipseProblem<ELEMENT>::
~UnstructuredImmersedEllipseProblem()
{
 // Delete Fluid timestepper
 delete this->time_stepper_pt(0);
 
 //Delete solid timestepper
 delete this->time_stepper_pt(1);
 
 // Kill data associated with outer boundary
 unsigned n=Outer_boundary_polygon_pt->npolyline();
 for (unsigned j=0;j<n;j++)
  {
   delete Outer_boundary_polygon_pt->polyline_pt(j);
  }
 delete Outer_boundary_polygon_pt;
 
 //Flush the drag mesh from the rigid body
 dynamic_cast<ImmersedRigidBodyElement*>(Rigid_body_pt[0])->
  flush_drag_mesh();

 // Flush Lagrange multiplier mesh
 delete_lagrange_multiplier_elements();
 delete Lagrange_multiplier_mesh_pt;
 
 // Delete error estimator
 delete Fluid_mesh_pt->spatial_error_estimator_pt();
 
 // Delete fluid mesh
 delete Fluid_mesh_pt;
 
 // Kill const eqn
 delete Problem_Parameter::Constitutive_law_pt;
 
 // Close norm and trace files
 this->Cog_exact_file.close();
 this->Cog_file.close();
 this->Norm_file.close();
}


//==========start_actions_before_adapt==================================
/// Actions before adapt: Wipe the mesh of Lagrange multiplier elements
//=======================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::actions_before_adapt()
{
 //Flush the drag mesh from the rigid body
 dynamic_cast<ImmersedRigidBodyElement*>(Rigid_body_pt[0])->
  flush_drag_mesh();

 //Kill the drag element
 delete_drag_elements();

 // Kill the  elements and wipe surface mesh
 delete_lagrange_multiplier_elements();

 // Rebuild the Problem's global mesh from its various sub-meshes
 this->rebuild_global_mesh();
} // end of actions_before_adapt


//============start_actions_after_adapt===================================== 
/// Actions after adapt: Rebuild the mesh of Lagrange multiplier elements
//==========================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::actions_after_adapt()
{
 //Reset the lagrangian coordinates for the solid mechanics
 //an updated lagrangian approach
 Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

 // Create the elements that impose the displacement constraint 
 create_lagrange_multiplier_elements();
 
 // Create the drag elements anew
 create_drag_elements();
 
 // Add the drag elements to the rigid body
 dynamic_cast<ImmersedRigidBodyElement*>(Rigid_body_pt[0])->
  set_drag_mesh(this->Drag_mesh_pt[0]);
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 this->rebuild_global_mesh();
 
 // Setup the problem again -- remember that fluid mesh has been
 // completely rebuilt and its element's don't have any
 // pointers to Re etc. yet
 complete_problem_setup();
 
 // Output solution after adaptation/projection
 bool doc_projection=true;
 doc_solution(doc_projection);
}// end of actions_after_adapt


//============start_complete_problem_setup=================================
///  Set boundary condition, assign auxiliary node update fct.
/// Complete the build of all elements, attach power elements that allow
/// computation of drag vector
//=========================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::complete_problem_setup()
{   
 // Set the boundary conditions for fluid problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned nbound=Fluid_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   unsigned num_nod=Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Cache pointer to node
     Node* const nod_pt=Fluid_mesh_pt->boundary_node_pt(ibound,inod);
     
     //Pin x-velocity unless on inlet (0) and outlet (2) boundaries
     //of the external rectangular box
     if((ibound!=0) && (ibound!=2)) {nod_pt->pin(0);}
     //Pin the y-velocity on all boundaries
     nod_pt->pin(1);
     
     // Pin pseudo-solid positions apart from on the 
     // ellipse boundary that is allowed to move
     // Cache cast pointer to solid node
     SolidNode* const solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);
          
     //Pin the solid positions on all external boundaries
     if(ibound < 4)
      {
       solid_node_pt->pin_position(0);
       solid_node_pt->pin_position(1);
      }
     // Unpin the position of all the nodes on hole boundaries:
     // since they will be moved using Lagrange Multiplier
     else
      {
       solid_node_pt->unpin_position(0);
       solid_node_pt->unpin_position(1);
       
       // Assign auxiliary node update fct, which determines the
       // velocity on the moving boundary using the position history
       // values
       // A more accurate version may be obtained by using velocity
       // based on the actual position of the geometric object,
       // but this introduces additional dependencies between the
       // Data of the rigid body and the fluid elements.
       nod_pt->set_auxiliary_node_update_fct_pt(
        FSI_functions::apply_no_slip_on_moving_wall); 
      }
    } //End of loop over boundary nodes
  } // End loop over boundaries
 
 // Complete the build of all elements so they are fully functional
 unsigned n_element = Fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
   
   // Set the Reynolds number
   el_pt->re_pt() = &Problem_Parameter::Re;
   
   // Set the Womersley number (same as Re since St=1)
   el_pt->re_st_pt() = &Problem_Parameter::Re;
   
   // Set the constitutive law for pseudo-elastic mesh deformation
   el_pt->constitutive_law_pt()=Problem_Parameter::Constitutive_law_pt;
   
   // Set the "density" for pseudo-elastic mesh deformation
   el_pt->lambda_sq_pt()=&Problem_Parameter::Lambda_sq;
  }
 
 // Re-apply Dirichlet boundary conditions for current and history values
 // (projection ignores boundary conditions!)
 this->set_boundary_velocity();

} //end_of_complete_problem_setup


//=================start_set_boundary_velocity===========================
///Set the boundary velocity for current and history values
//=======================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::
set_boundary_velocity()
{
 //Loop over top and bottom walls, inlet and outlet
 for(unsigned ibound=0;ibound<4;ibound++)
  {
   //If we are on the inlet or outlet only zero the 
   //y velocity, leave x alone
   if((ibound==0) || (ibound==2))
    {
     unsigned num_nod=this->Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* const nod_pt=this->Fluid_mesh_pt->boundary_node_pt(ibound,inod);
       
       // Get number of previous (history) values
       unsigned n_prev=nod_pt->time_stepper_pt()->nprev_values();
       
       //Zero all current and previous y-velocities
       for(unsigned t=0;t<=n_prev;t++)
        {
         nod_pt->set_value(t,1,0.0);
        }
      }
    }
   //Otherwise on the top and bottom walls set a smooth x-velocity
   //and zero y-velocity
   else
    {
     unsigned num_nod=this->Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* nod_pt=this->Fluid_mesh_pt->boundary_node_pt(ibound,inod);
       
       // Get number of previous (history) values
       unsigned n_prev=nod_pt->time_stepper_pt()->nprev_values();
       
       //Now set the boundary velocity
       double y = nod_pt->x(1);
       //Get the previous time
       for(unsigned t=0;t<=n_prev;t++)
        {
         //Get the time
         double time_ = this->time_pt()->time(t);
         
         //Get the velocity ramp
         //Initially zero (nothing at all is going on)
         double ramp = 0.0;
         double delta = 5.0;
         
         double e1 = exp(-delta);
         double a1 = 1.0 - (1.0 + delta + 0.5*delta*delta)*e1;
         double b1 = (3.0 + 3.0*delta + delta*delta)*e1 - 3.0;
         double c1 = 3.0 - (3.0 + 2.0*delta + 0.5*delta*delta)*e1;
         //Smooth start
         if((time_ >= 0.0) & (time_ <= 1.0)) 
          { 
           ramp = a1*time_*time_*time_
            + b1*time_*time_
              + c1*time_;
          }
         //Coupled to exponential levelling
         else if (time_ > 1.0)
          {
           ramp = 1.0 - exp(-delta*time_);
          }
         
         nod_pt->set_value(t,0,-y*ramp);
        }
       
       // Zero all current and previous veloc values
       // for the v-velocity
       for (unsigned t=0;t<=n_prev;t++)
        {
         nod_pt->set_value(t,1,0.0);
        }
      }
    }
  }
}


//====================start_of_pin_rigid_body=====================
///Pin the degrees of freedom associated with the solid bodies
//================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::pin_rigid_body()
{
 unsigned n_rigid_body = Rigid_body_pt.size();
 for(unsigned i=0;i<n_rigid_body;++i)
  {
   unsigned n_geom_data = Rigid_body_pt[i]->ngeom_data();
   for(unsigned j=0;j<n_geom_data;j++)
    {
     Rigid_body_pt[i]->geom_data_pt(j)->pin_all();
    }
  }
}


//=================start_unpin_rigid_body==========================
///Unpin the degrees of freedom associated with the solid bodies
//=================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::unpin_rigid_body()
{
 unsigned n_rigid_body = Rigid_body_pt.size();
 for(unsigned i=0;i<n_rigid_body;++i)
  {
   unsigned n_geom_data = Rigid_body_pt[i]->ngeom_data();
   for(unsigned j=0;j<n_geom_data;j++)
    {
     Rigid_body_pt[i]->geom_data_pt(j)->unpin_all();
    }
  }
}


//==========start_solve_for_consistent_nodal_positions================
///Assemble and solve a simplified problem that ensures that the 
///positions of the boundary nodes are consistent with the weak 
///imposition of the displacement boundary conditions on the surface
///of the ellipse.
//===================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::
solve_for_consistent_nodal_positions()
{
 //First pin all degrees of freedom in the rigid body
 this->pin_rigid_body();

 //Must reassign equation numbrs
 this->assign_eqn_numbers();

 //Do a steady solve to map the nodes to the boundary of the ellipse
 this->steady_newton_solve();

 //Now unpin the rigid body...
 this->unpin_rigid_body();

 //...and then repin the position of the centre of mass
 ImmersedRigidBodyElement* rigid_element1_pt = 
  dynamic_cast<ImmersedRigidBodyElement*>(Rigid_body_pt[0]);
 rigid_element1_pt->pin_centre_of_mass_coordinate(0);
 rigid_element1_pt->pin_centre_of_mass_coordinate(1);

 //and then reassign equation numbers
 this->assign_eqn_numbers();

} //end_solve_for_consistent_nodal_positions




//============start_of_create_lagrange_multiplier_elements===============
/// Create elements that impose the prescribed boundary displacement
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::
create_lagrange_multiplier_elements()
{ 
 // The idea is to apply Lagrange multipliers to the boundaries in 
 // the mesh that have associated geometric objects

 //Find the number of boundaries
 unsigned n_boundary = Fluid_mesh_pt->nboundary();

 // Loop over the boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   //Get the geometric object associated with the boundary
   GeomObject* boundary_geom_obj_pt = 
    Fluid_mesh_pt->boundary_geom_object_pt(b);

   //Only bother to do anything if there is a geometric object
   if(boundary_geom_obj_pt!=0)
    {
     // How many bulk fluid elements are adjacent to boundary b?
     unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
     
     // Loop over the bulk fluid elements adjacent to boundary b?
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk fluid element that is 
       // adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Fluid_mesh_pt->boundary_element_pt(b,e));
       
       //Find the index of the face of element e along boundary b
       int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
       
       // Create new element. Note that we use different Lagrange
       // multiplier fields for each distinct boundary (here indicated
       // by b.
       ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>* el_pt =
        new ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(
         bulk_elem_pt,face_index,b);   
       
       // Add it to the mesh
       Lagrange_multiplier_mesh_pt->add_element_pt(el_pt);
       
       // Set the GeomObject that defines the boundary shape and set
       // which bulk boundary we are attached to (needed to extract
       // the boundary coordinate from the bulk nodes)
       el_pt->set_boundary_shape_geom_object_pt(
        boundary_geom_obj_pt,b);
       
       // Loop over the nodes to pin Lagrange multiplier
       unsigned nnod=el_pt->nnode();
       for(unsigned j=0;j<nnod;j++)
        {
         Node* nod_pt = el_pt->node_pt(j);
         
         // How many nodal values were used by the "bulk" element
         // that originally created this node?
         unsigned n_bulk_value=el_pt->nbulk_value(j);
         
         // Pin two of the four Lagrange multipliers at vertices
         // This is not totally robust, but will work in this application
         unsigned nval=nod_pt->nvalue();
         if (nval==7)
          {
           for (unsigned i=0;i<2;i++) 
            { 
             // Pin lagrangian multipliers
             nod_pt->pin(n_bulk_value+2+i);
            }
          }
        }
      } // end loop over the element
    } //End of  case if there is a geometric object
  } //End of loop over boundaries
}
// end of create_lagrange_multiplier_elements


//===============start_delete_lagrange_multiplier_elements==================
///  Delete elements that impose the prescribed boundary displacement
/// and wipe the associated mesh
//==========================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::
delete_lagrange_multiplier_elements()
{
 // How many surface elements are in the surface mesh
 unsigned n_element = Lagrange_multiplier_mesh_pt->nelement();
 
 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete Lagrange_multiplier_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 Lagrange_multiplier_mesh_pt->flush_element_and_node_storage();
 
} // end of delete_lagrange_multiplier_elements




//============start_of_create_drag_elements===============
/// Create elements that calculate the drag and torque on
/// the obstacles in the fluid mesh
//=======================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::create_drag_elements()
{ 
 //The idea is only to attach drag elements to those
 //boundaries associated with each particular rigid body
 //The loop is slightly inefficient, but should work in general

 // Get the number of rigid bodies
 unsigned n_rigid = Rigid_body_pt.size();

 // Get the number of boundaries in the mesh
 unsigned n_boundary = Fluid_mesh_pt->nboundary();

 //Loop over the rigid bodies
 for(unsigned r=0;r<n_rigid;r++)
  {
   //Get the rigid boundary geometric object
   ImmersedRigidBodyElement* rigid_el_pt = 
    dynamic_cast<ImmersedRigidBodyElement*>(this->Rigid_body_pt[r]);
   
   // Loop over all boundaries
   for(unsigned b=0;b<n_boundary;b++)
    {
     //Does the boundary correspond to the current rigid body
     if(dynamic_cast<ImmersedRigidBodyElement*>
        (Fluid_mesh_pt->boundary_geom_object_pt(b)) == rigid_el_pt)
      {
       // How many bulk fluid elements are adjacent to boundary b?
       unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
       
       // Loop over the bulk fluid elements adjacent to boundary b?
       for(unsigned e=0;e<n_element;e++)
        {
         // Get pointer to the bulk fluid element that is 
         // adjacent to boundary b
         ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(b,e));
       
         //Find the index of the face of element e along boundary b
         int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
         
         // Create new element. Note that we use different Lagrange
         // multiplier fields for each distinct boundary (here indicated
         // by b.
         NavierStokesSurfaceDragTorqueElement<ELEMENT>* el_pt =
          new NavierStokesSurfaceDragTorqueElement<ELEMENT>(
           bulk_elem_pt,face_index);   
         
         // Add it to the mesh
         Drag_mesh_pt[r]->add_element_pt(el_pt);
         
         //Set the original centre of rotation 
         //from the rigid body in the drag element
         for(unsigned i=0;i<2;i++) 
          {el_pt->centre_of_rotation(i) = 
            rigid_el_pt->initial_centre_of_mass(i);}
         
         //Set the data to the translation and rotation data
         //as well because these data will enter the Jacobian terms
         //in the DragTorqueElement
         el_pt->set_translation_and_rotation(rigid_el_pt->geom_data_pt(0));
        } // end loop over the element
      }
    } //End of loop over boundaries
  } // end loop over the rigid bodies
}
// end of create_drag_elements


//=======================================================================
///  Delete elements that calculate the drag and torque on the 
/// boundaries
//=======================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::delete_drag_elements()
{
 unsigned n_bodies = Drag_mesh_pt.size();
 for(unsigned n=0;n<n_bodies;n++)
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Drag_mesh_pt[n]->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Drag_mesh_pt[n]->element_pt(e);
    }
   
   // Wipe the mesh
   Drag_mesh_pt[n]->flush_element_and_node_storage();
  }
} // end of delete_drag_elements




//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::doc_solution(
 const bool& project)
{ 

 oomph_info << "Docing step: " << this->Doc_info.number()
            << std::endl;

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 
 // Output solution and projection files
 if(!project)
  {
   sprintf(filename,"%s/soln%i.dat",
           this->Doc_info.directory().c_str(),
           this->Doc_info.number());
  }
 else
  {
   sprintf(filename,"%s/proj%i.dat",
           this->Doc_info.directory().c_str(),
           this->Doc_info.number()-1);
  }

 // Assemble square of L2 norm 
 double square_of_l2_norm=0.0;
 unsigned nel=Fluid_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   square_of_l2_norm+=
    dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(e))->
    square_of_l2_norm();
  }

 std::cout << "Checking " << sqrt(square_of_l2_norm) << "\n";
 this->Norm_file << sqrt(square_of_l2_norm) << "\n";
 
 some_file.open(filename);
 some_file << dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(0))
  ->variable_identifier();
 this->Fluid_mesh_pt->output(some_file,npts);   
 some_file.close();

 // No trace file writing after projection
 if(project) return;

 //Output the boundary only if not projecting
 sprintf(filename,"%s/int%i.dat",
         this->Doc_info.directory().c_str(),
         this->Doc_info.number());
 some_file.open(filename);
 this->Lagrange_multiplier_mesh_pt->output(some_file);
 some_file.close();

 // Increment the doc_info number
 this->Doc_info.number()++;

 //Output the motion of the centre of gravity
 dynamic_cast<ImmersedRigidBodyElement*>(this->Rigid_body_pt[0])->
  output_centre_of_gravity(this->Cog_file);

 //Output the exact solution
 this->output_exact_solution(this->Cog_exact_file);

}


//=====================================================================  
/// Output the exact solution
//=====================================================================  
template<class ELEMENT>
void UnstructuredImmersedEllipseProblem<ELEMENT>::
output_exact_solution(std::ofstream &output_file)
  {
   //Get the current time
   double time = this->time();
   //Output the exact solution
   output_file << time << " " << Jeffery_Solution::angle(time)
               << " " << Jeffery_Solution::velocity(time)
               << " " << Jeffery_Solution::acceleration(time) << std::endl;
  }



//==========start_of_main======================================
/// Driver code for immersed ellipse problem
//============================================================
int main(int argc, char **argv)
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
  
 // Create problem in initial configuration
 UnstructuredImmersedEllipseProblem<
 ProjectableTaylorHoodElement<MyTaylorHoodElement> > problem;  

 //Initially ensure that the nodal positions are consistent with 
 //their weak imposition
 problem.solve_for_consistent_nodal_positions();
 
 // Initialise timestepper
 double dt=0.05;
 problem.initialise_dt(dt);

 // Perform impulsive start
 problem.assign_initial_values_impulsive();

 // Output initial conditions
 problem.doc_solution();

 // Solve problem a few times on given mesh
 unsigned nstep=3;
 for (unsigned i=0;i<nstep;i++)
  {
   // Solve the problem
   problem.unsteady_newton_solve(dt);    
   problem.doc_solution();
  }

 // Now do a couple of adaptations
 unsigned ncycle=200;
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   ncycle=1;
   oomph_info << "Only doing one cycle during validation\n";
  }

 for (unsigned j=0;j<ncycle;j++)
  {       
   // Adapt the problem
   problem.adapt();

   //Solve problem a few times
   for (unsigned i=0;i<nstep;i++)
    {     
     // Solve the problem
     problem.unsteady_newton_solve(dt);
     problem.doc_solution();
    }
  }

} //end of main
