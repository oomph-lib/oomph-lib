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
#include <iostream>
#include <cstring>

// Generic oomph-lib includes
#include "generic.h"
#include "navier_stokes.h"
#include "beam.h"

// The wall mesh
#include "meshes/one_d_lagrangian_mesh.h"

//Include the fluid mesh
#include "meshes/collapsible_channel_mesh.h"

using namespace std;

using namespace oomph;


//==========start_of_BL_Squash =========================================
/// Namespace to define the mapping [0,1] -> [0,1] that re-distributes
/// nodal points across the channel width.
//======================================================================
namespace BL_Squash
{
 
 /// Boundary layer width
 double Delta=0.1;

 /// Fraction of points in boundary layer
 double Fract_in_BL=0.5;

 /// Mapping [0,1] -> [0,1] that re-distributes
 /// nodal points across the channel width
 double squash_fct(const double& s)
 {
  // Default return
  double y=s;
  if (s<0.5*Fract_in_BL)
   {
    y=Delta*2.0*s/Fract_in_BL;
   }
  else if (s>1.0-0.5*Fract_in_BL)
   {
    y=2.0*Delta/Fract_in_BL*s+1.0-2.0*Delta/Fract_in_BL;
   }
  else
   {
    y=(1.0-2.0*Delta)/(1.0-Fract_in_BL)*s+
      (Delta-0.5*Fract_in_BL)/(1.0-Fract_in_BL);
   }

  return y;
 }
}// end of BL_Squash






//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//====start_of_underformed_wall============================================
/// Undeformed wall is a steady, straight 1D line in 2D space 
///  \f[ x = X_0 + \zeta \f]
///  \f[ y = H \f]
//=========================================================================
class UndeformedWall : public GeomObject
{

public:

 /// Constructor: arguments are the starting point and the height
 /// above y=0.
 UndeformedWall(const double& x0, const double& h): GeomObject(1,2)
  {
   X0=x0;
   H=h;
  }
 

 /// Position vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position Vector
   r[0] = zeta[0]+X0;
   r[1] = H;
  }


 /// Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep. Calls steady version.
 void position(const unsigned& t, const Vector<double>& zeta,
               Vector<double>& r) const
  {
   // Use the steady version
   position(zeta,r);

  } // end of position


 /// Posn vector and its  1st & 2nd derivatives
 /// w.r.t. to coordinates:
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i). 
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ = 
 /// ddrdzeta(alpha,beta,i). Evaluated at current time.
 void d2position(const Vector<double>& zeta,
                 Vector<double>& r,
                 DenseMatrix<double> &drdzeta,
                 RankThreeTensor<double> &ddrdzeta) const
  {
   // Position vector
   r[0] = zeta[0]+X0;
   r[1] = H;

   // Tangent vector
   drdzeta(0,0)=1.0;
   drdzeta(0,1)=0.0;

   // Derivative of tangent vector
   ddrdzeta(0,0,0)=0.0;
   ddrdzeta(0,0,1)=0.0;

  } // end of d2position

 private :

 /// x position of the undeformed beam's left end. 
 double X0;

 /// Height of the undeformed wall above y=0.
 double H;

}; //end_of_undeformed_wall


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//====start_of_control_flags============================
/// Namespace for phyical parameters
//======================================================
namespace Control_Flags
{
 /// Normal run or validation run?
 bool Validation_run=false;

 /// Steady run or unsteady run?
 bool Steady_run=false;
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//====start_of_physical_parameters=====================
/// Namespace for phyical parameters
//======================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=50.0;

 /// Womersley = Reynolds times Strouhal
 double ReSt=50.0;
 
 /// Default pressure on the left boundary
 double P_up=0.0;
 
 /// Prescribed y coordinate of control point
 double Prescribed_y=1.0;

 /// Pointer to pressure load
 Data* Pext_data_pt=0;

 /// Access function to value of external pressure
 double &external_pressure() {return *Pext_data_pt->value_pt(0);}

 /// Traction applied on the fluid at the left (inflow) boundary
 void prescribed_traction(const double& t,
                          const Vector<double>& x,
                          const Vector<double>& n,
                          Vector<double>& traction)
 {
  traction.resize(2);
  traction[0]=P_up;
  traction[1]=0.0;

 } //end traction

 /// Non-dimensional wall thickness. As in Jensen & Heil (2003) paper.
 double H=1.0e-2;
 
 /// 2nd Piola Kirchhoff pre-stress. As in Jensen & Heil (2003) paper.
 double Sigma0=1.0e3;

 /// Load function: Apply a constant external pressure to the wall.
 /// Note:  This is the load without the fluid contribution!
 /// Fluid load gets added on by FSIWallElement.
 void load(const Vector<double>& xi, const Vector<double>& x,
           const Vector<double>& N, Vector<double>& load)
 { 
  double p_ext = Pext_data_pt->value(0);
  for(unsigned i=0;i<2;i++) 
   {
    load[i] = -p_ext*N[i];
   }
 } //end of load


 /// Fluid structure interaction parameter: Ratio of stresses used for
 /// non-dimensionalisation of fluid to solid stresses. 
 double Q=1.0e-5;


} // end of namespace




//====start_of_problem_class==========================================
///Problem class
//====================================================================
template <class ELEMENT>
class FSICollapsibleChannelProblem : public Problem
{

 public :

/// Constructor: The arguments are the number of elements and
/// the lengths of the domain.
 FSICollapsibleChannelProblem(
  const unsigned& nup, 
  const unsigned& ncollapsible,
  const unsigned& ndown,
  const unsigned& ny,
  const double& lup,
  const double& lcollapsible, 
  const double& ldown,
  const double& ly,
  const bool& wall_jacobian_ignores_fluid_shear_stress_data,
  const bool& wall_jacobian_ignores_geometric_data,
  const bool& fluid_jacobian_ignores_geometric_data);
 
 /// Destructor (empty)
 ~FSICollapsibleChannelProblem(){}


 /// Access function for the specific bulk (fluid) mesh
 AlgebraicCollapsibleChannelMesh<ELEMENT>* bulk_mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<
    AlgebraicCollapsibleChannelMesh<ELEMENT>*>
    (Bulk_mesh_pt);
  }



 /// Access function for the wall mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* wall_mesh_pt() 
  {
   return Wall_mesh_pt;

  } // end of access to wall mesh


 /// Update the problem specs before solve (empty) 
 void actions_before_newton_solve(){}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}
  
 /// Update before checking Newton convergence: Update the
 /// nodal positions in the fluid mesh in response to possible 
 /// changes in the wall shape
 void actions_before_newton_convergence_check()
  {
   Bulk_mesh_pt->node_update();
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info,ofstream& trace_file);
 
 /// Apply initial conditions
 void set_initial_condition();

private : 

 /// Create the prescribed traction elements on boundary b
 void create_traction_elements(const unsigned &b, 
                               Mesh* const &bulk_mesh_pt,
                               Mesh* const &traction_mesh_pt);
 
 ///Number of elements in the x direction in the upstream part of the channel
 unsigned Nup;

 /// Number of elements in the x direction in the collapsible part of 
 /// the channel
 unsigned Ncollapsible;

 ///Number of elements in the x direction in the downstream part of the channel
 unsigned Ndown;

 ///Number of elements across the channel
 unsigned Ny;

 ///x-length in the upstream part of the channel
 double Lup;

 ///x-length in the collapsible part of the channel
 double Lcollapsible;

 ///x-length in the downstream part of the channel
 double Ldown;

 ///Transverse length
 double Ly;

 /// Pointer to the "bulk" mesh
 AlgebraicCollapsibleChannelMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh that applies the traction at the
 /// inflow
 Mesh* Applied_fluid_traction_mesh_pt; 

 /// Pointer to mesh containing the displacement control element (only!)
 Mesh* Displacement_control_mesh_pt;
 
 /// Pointer to the "wall" mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* Wall_mesh_pt; 

 ///Pointer to the left control node
 Node* Left_node_pt;
 
 ///Pointer to right control node
 Node* Right_node_pt;
 
 /// Pointer to control node on the wall
 Node* Wall_node_pt;
 
 /// Flag controlling whether geometric data is ignored when 
 /// the fluid elements calculate their Jacobian contribution
 bool Fluid_jacobian_ignores_geometric_data;
 
 /// Flag controlling whether geometric data is ignored when 
 /// the wall elements calculate their Jacobian contribution
 bool Wall_jacobian_ignores_geometric_data;

 /// Flag controlling whether fluid shear stress is ignored when 
 /// wall elements calculate their Jacobian contribition
 bool Wall_jacobian_ignores_fluid_shear_stress_data;

};//end of problem class




//=====start_of_constructor======================================
/// Constructor for the collapsible channel problem
//===============================================================
template <class ELEMENT>
FSICollapsibleChannelProblem<ELEMENT>::FSICollapsibleChannelProblem(
 const unsigned& nup, 
 const unsigned& ncollapsible,
 const unsigned& ndown,
 const unsigned& ny,
 const double& lup,
 const double& lcollapsible, 
 const double& ldown,
 const double& ly,
 const bool& wall_jacobian_ignores_fluid_shear_stress_data,
 const bool& wall_jacobian_ignores_geometric_data,
 const bool& fluid_jacobian_ignores_geometric_data)
{
 // Store problem parameters
 Nup=nup;
 Ncollapsible=ncollapsible;
 Ndown=ndown;
 Ny=ny;
 Lup=lup;
 Lcollapsible=lcollapsible;
 Ldown=ldown;
 Ly=ly;
 Wall_jacobian_ignores_fluid_shear_stress_data = 
  wall_jacobian_ignores_fluid_shear_stress_data;
 Wall_jacobian_ignores_geometric_data=wall_jacobian_ignores_geometric_data;
 Fluid_jacobian_ignores_geometric_data=fluid_jacobian_ignores_geometric_data; 

 // Overwrite maximum allowed residual to accomodate bad initial guesses
 Problem::Max_residuals=100000.0;

 // Allow a few extra Newton iterations
 Max_newton_iterations=20;

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 add_time_stepper_pt(new BDF<2>);

 // Geometric object that represents the undeformed wall: 
 // A straight line at height y=ly; starting at x=lup.
 UndeformedWall* undeformed_wall_pt=new UndeformedWall(lup,ly);

 //Create the "wall" mesh with FSI Hermite elements
 Wall_mesh_pt = new OneDLagrangianMesh<FSIHermiteBeamElement>
  //(2*Ncollapsible+5,Lcollapsible,undeformed_wall_pt);
  (Ncollapsible,Lcollapsible,undeformed_wall_pt);
 

 // Build a geometric object (one Lagrangian, two Eulerian coordinates)
 // from the wall mesh
 MeshAsGeomObject* wall_geom_object_pt=
  new MeshAsGeomObject(Wall_mesh_pt); 


 //Build bulk (fluid) mesh
 Bulk_mesh_pt = 
  new AlgebraicCollapsibleChannelMesh<ELEMENT>
  (nup, ncollapsible, ndown, ny,
   lup, lcollapsible, ldown, ly,
   wall_geom_object_pt,
   &BL_Squash::squash_fct,
   time_stepper_pt());


 // Create "surface mesh" that will contain only the prescribed-traction 
 // elements. The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Applied_fluid_traction_mesh_pt = new Mesh;
 
 // Create prescribed-traction elements from all elements that are 
 // adjacent to boundary 5 (left boundary), but add them to a separate mesh
 create_traction_elements(5,Bulk_mesh_pt,Applied_fluid_traction_mesh_pt);
  

 // Complete build of fluid mesh
 //----------------------------- 
 
 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 unsigned n_element=Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Set flag which controls whether geometric data is ignored
   // when calculating the fluid Jacobian
   if(Fluid_jacobian_ignores_geometric_data)
    {
     el_pt->enable_bypass_fill_in_jacobian_from_geometric_data();
    }
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
   
  } // end loop over elements



 // Apply boundary conditions for fluid
 //------------------------------------

 //Pin the velocity on the boundaries
 //x and y-velocities pinned along boundary 0 (bottom boundary) :
 unsigned ibound=0; 
 unsigned num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   for(unsigned i=0;i<2;i++)
    {
     bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(i);
    }
  }
  
 //x and y-velocities pinned along boundaries 2, 3, 4 (top boundaries) :
 for(ibound=2;ibound<5;ibound++)
  { 
   num_nod= bulk_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     for(unsigned i=0;i<2;i++)
      {
       bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(i);
      }
    }
  }

 //y-velocity pinned along boundary 1 (right boundary):
 ibound=1; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(1);
  }


 //y-velocity pinned along boundary 5 (left boundary):
 ibound=5; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(1);
  }
//end of pin_velocity

 // Complete build of applied traction elements
 //--------------------------------------------

 // Loop over the traction elements to pass pointer to prescribed 
 // traction function
 unsigned n_el=Applied_fluid_traction_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   // Upcast from GeneralisedElement to NavierStokes traction element
   NavierStokesTractionElement<ELEMENT> *el_pt = 
    dynamic_cast< NavierStokesTractionElement<ELEMENT>*>(
     Applied_fluid_traction_mesh_pt->element_pt(e));
    
   // Set the pointer to the prescribed traction function
   el_pt->traction_fct_pt() = &Global_Physical_Variables::prescribed_traction;
  }


 //-------------------------------------------------------------
 // Set up displacement control mesh if required for steady run
 //-------------------------------------------------------------
 Displacement_control_mesh_pt=0;

 if (Control_Flags::Steady_run)
  {
   // Choose element in which displacement control is applied as
   // halfway along wall
   unsigned nel_ctrl=Ncollapsible/2;
   
   // Deal with even/odd number of wall elements 
   Vector<double> s_displ_control(1);
   if (n_element%2==1)
    {
     s_displ_control[0]=0.0;
    }
   else
    {
     s_displ_control[0]=-1.0;
    }
   
   // Controlled element
   FSIHermiteBeamElement* controlled_element_pt=
    dynamic_cast<FSIHermiteBeamElement*>(Wall_mesh_pt->element_pt(nel_ctrl));

   // Fix the displacement in the y (1) direction...
   unsigned controlled_direction=1;

   // Output control point
   Vector<double> xi(1);
   Vector<double> x(2);
   controlled_element_pt->interpolated_xi(s_displ_control,xi);
   controlled_element_pt->interpolated_x(s_displ_control,x);

   oomph_info << "Displacement control for element: " 
              << nel_ctrl << std::endl;
   oomph_info << "Displacement control applied at xi = ("
              << xi[0] << ")" << std::endl;
   oomph_info << "Corresponding to                x  = ("
              << x[0] << ", " << x[1] << ")" << std::endl;

   // Build displacement control element
   DisplacementControlElement* displ_control_el_pt=
    new DisplacementControlElement(controlled_element_pt,
                                   s_displ_control,
                                   controlled_direction,
                                   &Global_Physical_Variables::Prescribed_y);

   // set Pext_data_pt to point at displacement control load
   Global_Physical_Variables::Pext_data_pt=
    displ_control_el_pt->displacement_control_load_pt();

   // Create dedicated mesh
   Displacement_control_mesh_pt=new Mesh;

   // Add the displacement-control element to the mesh
   Displacement_control_mesh_pt->add_element_pt(displ_control_el_pt);
  }
 else
  {
   // Create Data object whose one-and-only value contains the
   // (in principle) adjustable load
   Global_Physical_Variables::Pext_data_pt=new Data(1);

   // Pin the external pressure because it isn't actually adjustable.
   Global_Physical_Variables::Pext_data_pt->pin(0);

   // Set the value of external pressure
   Global_Physical_Variables::external_pressure() = 0.0;
  }




 // Complete build of wall elements
 //--------------------------------
  
 //Loop over the elements to set physical parameters etc.
 n_element = wall_mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast to the specific element type
   FSIHermiteBeamElement *elem_pt = 
    dynamic_cast<FSIHermiteBeamElement*>(wall_mesh_pt()->element_pt(e));

   // Do elememts take account of fluid shear stress when calculating Jacobian?
   if(Wall_jacobian_ignores_fluid_shear_stress_data)
    {
     elem_pt->disable_shear_stress_in_jacobian();
    }
   else
    {
     elem_pt->enable_shear_stress_in_jacobian();
    }
    
   // Set physical parameters for each element:
   elem_pt->sigma0_pt() = &Global_Physical_Variables::Sigma0;
   elem_pt->h_pt() = &Global_Physical_Variables::H;
    
   // Set the load vector for each element
   elem_pt->load_vector_fct_pt() = &Global_Physical_Variables::load;

   // Function that specifies the load ratios
   elem_pt->q_pt() = &Global_Physical_Variables::Q;

   // Set the undeformed shape for each element
   elem_pt->undeformed_beam_pt() = undeformed_wall_pt;


   // The normal on the wall elements as computed by the FSIHermiteElements
   // points away from the fluid rather than into the fluid (as assumed
   // by default)
   elem_pt->set_normal_pointing_out_of_fluid();

   // Displacement control for steady problem? If so, the load on *all* 
   // elements is affected by an unknown -- the external pressure, stored
   // as the one-and-only value in a Data object: Add it to the
   // elements' external Data.
   if (Control_Flags::Steady_run)
    {
     //The external pressure is external data for all elements
     elem_pt->add_external_data(Global_Physical_Variables::Pext_data_pt);
    }

  } // end of loop over elements



 // Boundary conditions for wall mesh
 //----------------------------------

 // Set the boundary conditions: Each end of the beam is fixed in space
 // Loop over the boundaries (ends of the beam)
 for(unsigned b=0;b<2;b++)
  {
   // Pin displacements in both x and y directions
   wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(0); 
   wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(1);
  }
  
 //Choose control nodes
 //---------------------
  
 // Left boundary
 ibound=5; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 unsigned control_nod=num_nod/2;
 Left_node_pt= bulk_mesh_pt()->boundary_node_pt(ibound, control_nod);
  
 // Right boundary
 ibound=1; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 control_nod=num_nod/2;
 Right_node_pt= bulk_mesh_pt()->boundary_node_pt(ibound, control_nod);
  
 
 // Set the pointer to the control node on the wall
 Wall_node_pt=wall_mesh_pt()->node_pt(Ncollapsible/2);

 // Add the sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Applied_fluid_traction_mesh_pt);
 add_sub_mesh(Wall_mesh_pt);
 if (Displacement_control_mesh_pt != 0)
  {
   add_sub_mesh(Displacement_control_mesh_pt);
  }

 // Combine all submeshes into a single Mesh
 build_global_mesh();

 // Setup FSI
 //----------

 // The velocity of the fluid nodes on the wall (fluid mesh boundary 3)
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 ibound=3; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   static_cast<AlgebraicNode*>(
    bulk_mesh_pt()->boundary_node_pt(ibound, inod))->
    set_auxiliary_node_update_fct_pt(
     FSI_functions::apply_no_slip_on_moving_wall);
  }
  
  
 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary is boundary 3 of the 
 // 2D fluid mesh.
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,3,Bulk_mesh_pt,Wall_mesh_pt);


 // Setup equation numbering scheme
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
  

}//end of constructor




//====start_of_doc_solution===================================================
/// Doc the solution
//============================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>:: doc_solution(DocInfo& doc_info, 
                                                          ofstream& trace_file)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output fluid solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 bulk_mesh_pt()->output(some_file,npts);
 some_file.close();

 // Document the wall shape
 sprintf(filename,"%s/beam%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 wall_mesh_pt()->output(some_file,npts);
 some_file.close();
   

 // Write trace file 
 trace_file << time_pt()->time() << " "
            << Wall_node_pt->x(1) << " "
            << Left_node_pt->value(0) << " "
            << Right_node_pt->value(0) << " "
            << Global_Physical_Variables::external_pressure()  << " " 
            << std::endl; 

} // end_of_doc_solution




//=====start_of_create_traction_elements======================================
/// Create the traction elements
//============================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::create_traction_elements(
 const unsigned &b, Mesh* const &bulk_mesh_pt, Mesh* const &traction_mesh_pt)
{

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
    (bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-traction element
   NavierStokesTractionElement<ELEMENT>* flux_element_pt = 
    new  NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);
   
   //Add the prescribed-traction element to the surface mesh
   traction_mesh_pt->add_element_pt(flux_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_traction_elements




//====start_of_apply_initial_condition========================================
/// Apply initial conditions
//============================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::set_initial_condition()
{ 
 // Check that timestepper is from the BDF family
 if (time_stepper_pt()->type()!="BDF")
  {
   std::ostringstream error_stream;
   error_stream << "Timestepper has to be from the BDF family!\n"
                << "You have specified a timestepper from the "
                << time_stepper_pt()->type() << " family" << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Update the mesh
 bulk_mesh_pt()->node_update();
 
 // Loop over the nodes to set initial guess everywhere
 unsigned num_nod = bulk_mesh_pt()->nnode();
 for (unsigned n=0;n<num_nod;n++)
  {
   // Get nodal coordinates
   Vector<double> x(2);
   x[0]=bulk_mesh_pt()->node_pt(n)->x(0);
   x[1]=bulk_mesh_pt()->node_pt(n)->x(1);
   
   // Assign initial condition: Steady Poiseuille flow
   bulk_mesh_pt()->node_pt(n)->set_value(0,6.0*(x[1]/Ly)*(1.0-(x[1]/Ly)));
   bulk_mesh_pt()->node_pt(n)->set_value(1,0.0);
  } 

 // Assign initial values for an impulsive start
 bulk_mesh_pt()->assign_initial_values_impulsive();

} // end of set_initial_condition





//============start_of_main====================================================
/// Driver code for a collapsible channel problem with FSI.
///
/// Run with the -help flag for more information on command line arguments to
/// run with an approximation to the Jacobian, to run steady calculations and
/// to perform validation runs with coarse resolution and small number of 
/// timesteps. 
///
/// Disclaimer: this code is intended to demonstrate how to ask for various
/// approximations to be made to the Jacobian matrix which can result in reduced
/// calculations per Newton step (and to test this).However these approximations
/// to the Jacobian may increase the overall number of Newton steps, or even 
/// cause the Newton solve to diverge. 
/// Therefore for any particular FSI problem some experiementation is required
/// to determine whether approximating the Jacobian is beneficial.
/// Note: In tests using this code, approximating the Jacobain increased
/// the overall time per Newton solve, however for a 3D equivalent to this 
/// problem significant speed ups were found.
//=============================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Parse command line
 int arg_index = 0;
 bool print_help = false;
 bool approx_jacobian = false;

 while (arg_index < argc)
  {
   if ( strcmp(argv[arg_index], "-validation_run") == 0 )
    {
     arg_index++;
     Control_Flags::Validation_run=true;
    }
   else if ( strcmp(argv[arg_index], "-steady_run") == 0 )
    {
     arg_index++;
     Control_Flags::Steady_run=true;
    }
   else if ( strcmp(argv[arg_index], "-approx_jacobian") == 0 )
    {
     arg_index++;
     approx_jacobian=true;
    }
   else if ( strcmp(argv[arg_index], "-help") == 0 )
    {
     print_help = true;
     break;
    }
   else
    {
     arg_index++;
    }
  }
 
 if (print_help)
  {
   oomph_info << "\n\nOption flags:\n\n";
   oomph_info << "-validation run   Perform validation run with coarse mesh\n"
              << "                  and fewer Newton solves\n\n";
   oomph_info << "-steady_run       Solve series of steady problems with\n"
              << "                  different wall positions, otherwise\n"
              << "                  solve unsteady problem\n\n";
   oomph_info << "-approx_jacobian  Use Jacobian approximation\n\n";
   return (0);
  }

 // Set flags controlling approximations to the Jacobian
 bool wall_jacobian_ignores_fluid_shear_stress_data = false;
 bool wall_jacobian_ignores_geometric_data = false;
 bool fluid_jacobian_ignores_geometric_data = false;
 if (approx_jacobian)
  {
   wall_jacobian_ignores_fluid_shear_stress_data = true;
   wall_jacobian_ignores_geometric_data = true;
   
   if (Control_Flags::Steady_run)
    {
     // This can speed up the calculation of the Jacobian but causes the
     // Newton method to fail for unsteady problems
     fluid_jacobian_ignores_geometric_data = true;
    }
  }
 
 // Reduction in resolution for validation run?
 unsigned coarsening_factor=1;
 if (Control_Flags::Validation_run)
  {
   coarsening_factor=4;
  }

 // Number of elements in the domain
 unsigned nup=20/coarsening_factor;
 unsigned ncollapsible=40/coarsening_factor;
 unsigned ndown=40/coarsening_factor;
 unsigned ny=16/coarsening_factor;

 // Length of the domain
 double lup=5.0;
 double lcollapsible=10.0;
 double ldown=10.0;
 double ly=1.0;

 // Build the problem with QTaylorHoodElements
 FSICollapsibleChannelProblem<AlgebraicElement<QTaylorHoodElement<2> > > 
  problem(nup,
          ncollapsible,
          ndown,
          ny, 
          lup,
          lcollapsible,
          ldown,
          ly,
          wall_jacobian_ignores_fluid_shear_stress_data,
          wall_jacobian_ignores_geometric_data,
          fluid_jacobian_ignores_geometric_data);

 // Set external pressure (on the wall stiffness scale). 
 Global_Physical_Variables::external_pressure() = 1.0e-1;
 
 // Pressure on the left boundary: This is consistent with steady
 // Poiseuille flow
 Global_Physical_Variables::P_up=12.0*(lup+lcollapsible+ldown);

 // Set initial wall displacement
 Global_Physical_Variables::Prescribed_y = 1.0;

 // Set Re and ReSt - it appears small values are required in this 2D problem, 
 // with the Jacobian approximations otherwise the Newton method converges, this
 // is not the experience found with the 3D equivalent problem.
 Global_Physical_Variables::Re=1;
 Global_Physical_Variables::ReSt=Global_Physical_Variables::Re;

 // Timestep. Note: Preliminary runs indicate that the period of
 // the oscillation is about 1 so this gives us 40 steps per period.
 double dt=1.0/40.0; 
 
 // Initial time for the simulation
 double t_min=0.0;
 
 // Maximum time for simulation
 double t_max=3.5; 
 
 // Initialise timestep 
 problem.time_pt()->time()=t_min;
 problem.initialise_dt(dt);

 // Apply initial condition
 problem.set_initial_condition();
 
 //Set output directory
 DocInfo doc_info;
 if (Control_Flags::Validation_run)
  {
   if (approx_jacobian)
    {
     oomph_info << "Using approximate Jacobian\n";
     if (Control_Flags::Steady_run)
      doc_info.set_directory("RESLT_APPROX_STEADY");
     else
      doc_info.set_directory("RESLT_APPROX_UNSTEADY");
    }
   else
    {
     oomph_info << "Using exact Jacobian\n";
     if (Control_Flags::Steady_run)
      doc_info.set_directory("RESLT_EXACT_STEADY");
     else
      doc_info.set_directory("RESLT_EXACT_UNSTEADY");
    }
  }
 else
  {
   doc_info.set_directory("RESLT");
  }
 
 // Open a trace file 
 ofstream trace_file;
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);
   
 // Output the initial solution
 problem.doc_solution(doc_info, trace_file);
 
 // Increment step number
 doc_info.number()++;
 
 if (Control_Flags::Steady_run)
  {
   // Number of solves (reduced for validation)
   unsigned nstep = 20;
   if (Control_Flags::Validation_run)
    {
     nstep=3;
    }

   // Loop over different wall displacement
   for (unsigned istep=0;istep<nstep;istep++)
    {
     oomph_info << "\nSteady Newton solve " << istep
                << " for wall displacement of " 
                << Global_Physical_Variables::Prescribed_y << "\n"
                << "====================================================\n";
     
     // Solve the problem
     problem.newton_solve();

     cout << "External pressure="
          << Global_Physical_Variables::external_pressure()
          << "\n";
     
     // Output the solution
     problem.doc_solution(doc_info, trace_file);
     
     // Step number
     doc_info.number()++;

     // Change wall position
     Global_Physical_Variables::Prescribed_y -= 0.01;
    }
  }
 else
  {
   // Find number of timesteps (reduced for validation)
   unsigned nstep = unsigned((t_max-t_min)/dt);
   if (Control_Flags::Validation_run)
    {
     nstep=10;
    }
   
   // Timestepping loop
   for (unsigned istep=0;istep<nstep;istep++)
    {
     oomph_info << "\nNewton solve " << istep << "\n"
                << "================\n";
     
     // Solve the problem
     problem.unsteady_newton_solve(dt);
     
     // Output the solution
     problem.doc_solution(doc_info, trace_file);
     
     // Step number
     doc_info.number()++;
    }
   
  }

 // Close trace file.
 trace_file.close();

}//end of main
  
