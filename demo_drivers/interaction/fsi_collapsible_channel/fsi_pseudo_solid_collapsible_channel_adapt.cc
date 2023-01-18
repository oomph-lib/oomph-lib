//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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

// Generic oomph-lib includes
#include "generic.h"
#include "navier_stokes.h"
#include "beam.h"
#include "solid.h"
#include "constitutive.h"

// The wall mesh
#include "meshes/one_d_lagrangian_mesh.h"

//Include the fluid mesh
#include "meshes/collapsible_channel_mesh.h"

using namespace std;

using namespace oomph;


//======================================================================
/// Upgrade mesh to solid mesh
//======================================================================
template <class ELEMENT>
class ElasticRefineableCollapsibleChannelMesh :
 public virtual RefineableCollapsibleChannelMesh<ELEMENT>,
 public virtual SolidMesh
{


public:

 /// Constructor: Build mesh and copy Eulerian coords to Lagrangian
 /// ones so that the initial configuration is the stress-free one.
 ElasticRefineableCollapsibleChannelMesh<ELEMENT>(const unsigned& nup, 
                                                  const unsigned& ncollapsible,
                                                  const unsigned& ndown, 
                                                  const unsigned& ny,
                                                  const double& lup, 
                                                  const double& lcollapsible, 
                                                  const double& ldown, 
                                                  const double& ly,
                                                  GeomObject* wall_pt,
                                                  TimeStepper* time_stepper_pt=
                                                  &Mesh::Default_TimeStepper) :
  CollapsibleChannelMesh<ELEMENT>(nup, ncollapsible, ndown, ny,
                                  lup, lcollapsible, ldown, ly,
                                  wall_pt,
                                  time_stepper_pt),
  RefineableCollapsibleChannelMesh<ELEMENT>(nup, ncollapsible, ndown, ny,
                                            lup, lcollapsible, ldown, ly,
                                            wall_pt,
                                            time_stepper_pt)

  {
   /// Make the current configuration the undeformed one by
   /// setting the nodal Lagrangian coordinates to their current
   /// Eulerian ones
   set_lagrangian_nodal_coordinates();
  }

};




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



/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////


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
 virtual void d2position(const Vector<double>& zeta,
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

 /// Pseudo-solid mass density
 double Lambda_sq=0.0;

 /// Pseudo-solid Poisson ratio
 double Nu=0.1;

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

 /// External pressure
 double P_ext=0.0;

 /// Load function: Apply a constant external pressure to the wall.
 /// Note:  This is the load without the fluid contribution!
 /// Fluid load gets added on by FSIWallElement.
 void load(const Vector<double>& xi, const Vector<double>& x,
           const Vector<double>& N, Vector<double>& load)
 { 
  for(unsigned i=0;i<2;i++) 
   {
    load[i] = -P_ext*N[i];
   }
 } //end of load


 /// Fluid structure interaction parameter: Ratio of stresses used for
 /// non-dimensionalisation of fluid to solid stresses. 
 double Q=1.0e-5;


} // end of namespace




//====start_of_problem_class==========================================
/// Problem class
//====================================================================
template <class ELEMENT>
class FSICollapsibleChannelProblem : public Problem
{

 public :

/// Constructor: The arguments are the number of elements and
/// the lengths of the domain.
 FSICollapsibleChannelProblem(const unsigned& nup, 
                       const unsigned& ncollapsible,
                       const unsigned& ndown,
                       const unsigned& ny,
                       const double& lup,
                       const double& lcollapsible, 
                       const double& ldown,
                       const double& ly);
 
 /// Destructor (empty)
 ~FSICollapsibleChannelProblem(){}

 /// Access function for the specific bulk (fluid) mesh
 ElasticRefineableCollapsibleChannelMesh<ELEMENT>* bulk_mesh_pt()
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<
    ElasticRefineableCollapsibleChannelMesh<ELEMENT>*>
    (Bulk_mesh_pt);
  }


 /// Access function for the wall mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* wall_mesh_pt() 
  {
   return Wall_mesh_pt;

  } // end of access to wall mesh


 /// Actions before adapt: Wipe the mesh of prescribed traction elements
 void actions_before_adapt();
 
 /// Actions after adapt: Rebuild the mesh of prescribed traction elements
 /// and reset FSI
 void actions_after_adapt();

 /// Update the problem specs before solve (empty) 
 void actions_before_newton_solve() {}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}
  
 /// Update no slip before Newton convergence check
 void actions_before_newton_convergence_check()
  {
   // Moving wall: No slip; this implies that the velocity needs
   // to be updated in response to wall motion
   unsigned ibound=3;
   unsigned num_nod=bulk_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Which node are we dealing with?
     Node* node_pt=bulk_mesh_pt()->boundary_node_pt(ibound,inod);
     
     // Apply no slip
     FSI_functions::apply_no_slip_on_moving_wall(node_pt);
    }
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
 
 /// Delete prescribed traction elements from the surface mesh
 void delete_traction_elements(Mesh* const &traction_mesh_pt);

 /// Create elements that enforce prescribed boundary motion
 /// by Lagrange multiplilers
 void create_lagrange_multiplier_elements();

 /// Delete elements that enforce prescribed boundary motion
 /// by Lagrange multiplilers
 void delete_lagrange_multiplier_elements();

 /// Number of elements in the x direction in the upstream part of the channel
 unsigned Nup;

 /// Number of elements in the x direction in the collapsible part of 
 /// the channel
 unsigned Ncollapsible;

 /// Number of elements in the x direction in the downstream part of the channel
 unsigned Ndown;

 /// Number of elements across the channel
 unsigned Ny;

 /// x-length in the upstream part of the channel
 double Lup;

 /// x-length in the collapsible part of the channel
 double Lcollapsible;

 /// x-length in the downstream part of the channel
 double Ldown;

 /// Transverse length
 double Ly;
 
 /// Pointer to the "bulk" mesh 
 ElasticRefineableCollapsibleChannelMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh that applies the traction at the
 /// inflow
 Mesh* Applied_fluid_traction_mesh_pt; 
 
 /// Pointers to meshes of Lagrange multiplier elements
 SolidMesh* Lagrange_multiplier_mesh_pt;

 /// Pointer to the "wall" mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* Wall_mesh_pt; 

 /// Geometric object incarnation of the wall mesh
 MeshAsGeomObject* Wall_geom_object_pt;

 /// Pointer to the left control node
 Node* Left_node_pt;
 
 /// Pointer to right control node
 Node* Right_node_pt;
 
 /// Pointer to control node on the wall
 Node* Wall_node_pt;

 /// Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt;

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
 const double& ly)
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


 // Overwrite maximum allowed residual to accomodate bad initial guesses
 Problem::Max_residuals=1000.0;

 // Allocate the timestepper for the Navier-Stokes equations
 BDF<2>* fluid_time_stepper_pt=new BDF<2>;

 // Add the fluid timestepper to the Problem's collection of timesteppers.
 add_time_stepper_pt(fluid_time_stepper_pt);

 // Create a dummy Steady timestepper that stores two history values
 Steady<2>* wall_time_stepper_pt = new Steady<2>;

 // Add the wall timestepper to the Problem's collection of timesteppers.
 add_time_stepper_pt(wall_time_stepper_pt);

 // Geometric object that represents the undeformed wall: 
 // A straight line at height y=ly; starting at x=lup.
 UndeformedWall* undeformed_wall_pt=new UndeformedWall(lup,ly);

 //Create the "wall" mesh with FSI Hermite beam elements, passing the
 //dummy wall timestepper to the constructor
 Wall_mesh_pt = new OneDLagrangianMesh<FSIHermiteBeamElement>
  (Ncollapsible,Lcollapsible,undeformed_wall_pt,wall_time_stepper_pt);
 
 // Build a geometric object (one Lagrangian, two Eulerian coordinates)
 // from the wall mesh
 Wall_geom_object_pt=
  new MeshAsGeomObject(Wall_mesh_pt); 



 //Build bulk (fluid) mesh
 Bulk_mesh_pt = 
  new ElasticRefineableCollapsibleChannelMesh<ELEMENT>
  (nup, ncollapsible, ndown, ny,
   lup, lcollapsible, ldown, ly,
   Wall_geom_object_pt,
   fluid_time_stepper_pt);

 // Set a non-trivial boundary-layer-squash function...
 Bulk_mesh_pt->bl_squash_fct_pt() = &BL_Squash::squash_fct; 

 // ... and update the nodal positions accordingly
 bool update_all_solid_nodes=true;
 Bulk_mesh_pt->node_update(update_all_solid_nodes);
 Bulk_mesh_pt->set_lagrangian_nodal_coordinates();

 if (CommandLineArgs::Argc>1)
  {
   // One round of uniform refinement
   Bulk_mesh_pt->refine_uniformly();
  }

 // Create "surface mesh" that will contain only the prescribed-traction 
 // elements. The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Applied_fluid_traction_mesh_pt = new Mesh;
 
 // Create prescribed-traction elements from all elements that are 
 // adjacent to boundary 5 (left boundary), but add them to a separate mesh.
 create_traction_elements(5,Bulk_mesh_pt,Applied_fluid_traction_mesh_pt);
 
 // Construct the mesh of elements that enforce prescribed boundary motion
 // by Lagrange multipliers
 Lagrange_multiplier_mesh_pt=new SolidMesh;
 create_lagrange_multiplier_elements();

 // Add the sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Applied_fluid_traction_mesh_pt);
 add_sub_mesh(Wall_mesh_pt);
 add_sub_mesh(Lagrange_multiplier_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();
   
 //Set errror estimator 
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 bulk_mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
 //Set the constitutive law
 Constitutive_law_pt = new GeneralisedHookean(&Global_Physical_Variables::Nu);


 // Complete build of fluid mesh
 //----------------------------- 
 
 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 unsigned n_element=Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {

   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
   
   //Set the constitutive law
   el_pt->constitutive_law_pt() = Constitutive_law_pt;

   // Density of pseudo-solid
   el_pt->lambda_sq_pt()=&Global_Physical_Variables::Lambda_sq;

  } // end loop over elements


 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());


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

  }//end of pin_velocity

 
 //Pin the nodal position on all boundaries apart from 3 (the moving wall)
 for (unsigned ibound=0;ibound<6;ibound++)
  {
   if (ibound!=3)
    {
     unsigned num_nod= bulk_mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       for(unsigned i=0;i<2;i++)
        {
         dynamic_cast<SolidNode*>(bulk_mesh_pt()->
                                  boundary_node_pt(ibound, inod))
          ->pin_position(i);
        }
      }
    }
  }

 // Now pin all nodal positions in the rigid bits
 unsigned nnod=bulk_mesh_pt()->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   SolidNode* nod_pt=dynamic_cast<SolidNode*>(bulk_mesh_pt()->node_pt(j));
   if ((nod_pt->x(0)<=lup)||(nod_pt->x(0)>=(lup+lcollapsible)))
    {
     for(unsigned i=0;i<2;i++)
      {
       nod_pt->pin_position(i);
      }
    }
  }
 

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




 // Complete build of wall elements
 //--------------------------------
  
 //Loop over the elements to set physical parameters etc.
 n_element = wall_mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast to the specific element type
   FSIHermiteBeamElement *elem_pt = 
    dynamic_cast<FSIHermiteBeamElement*>(wall_mesh_pt()->element_pt(e));
    
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
 num_nod= wall_mesh_pt()->nnode();
 Wall_node_pt=wall_mesh_pt()->node_pt(Ncollapsible/2);



 // Setup FSI
 //----------

 // The velocity of the fluid nodes on the wall (fluid mesh boundary 3)
 // is set by the wall motion -- hence the no-slip condition needs to be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 ibound=3; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   bulk_mesh_pt()->boundary_node_pt(ibound, inod)->
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
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
  

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
   
 // Loop over all elements do dump out previous solutions
 // (get the number of previous timesteps available from the wall
 // timestepper)
 unsigned nsteps=time_stepper_pt(1)->nprev_values();
 for (unsigned t=0;t<=nsteps;t++)
  {     
   sprintf(filename,"%s/wall%i-%i.dat",doc_info.directory().c_str(),
           doc_info.number(),t);
   some_file.open(filename);
   unsigned n_elem=wall_mesh_pt()->nelement();
   for (unsigned ielem=0;ielem<n_elem;ielem++)
    {
     dynamic_cast<FSIHermiteBeamElement*>(wall_mesh_pt()->element_pt(ielem))->
      output(t,some_file,npts);
    }
   some_file.close();
  } // end of output of previous solutions

   

 // Write trace file 
 trace_file << time_pt()->time() << " "
            << Wall_node_pt->x(1) << " "
            << Left_node_pt->value(0) << " "
            << Right_node_pt->value(0) << " "
            << Global_Physical_Variables::P_ext  << " " 
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
   
   //What is the index of the face of element e along boundary
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-traction element
   NavierStokesTractionElement<ELEMENT>* flux_element_pt = 
    new  NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);
   
   //Add the prescribed-traction element to the surface mesh
   traction_mesh_pt->add_element_pt(flux_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_traction_elements



//============start_of_delete_traction_elements==============================
/// Delete traction elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::
delete_traction_elements(Mesh* const &surface_mesh_pt)
{
 // How many surface elements are in the surface mesh
 unsigned n_element = surface_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete surface_mesh_pt->element_pt(e);
  }

 // Wipe the mesh
 surface_mesh_pt->flush_element_and_node_storage();

} // end of delete_traction_elements


//============start_of_create_lagrange_multiplier_elements===============
/// Create elements that impose the prescribed boundary displacement
//=======================================================================
template<class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::
create_lagrange_multiplier_elements()
{
 // Lagrange multiplier elements are located on boundary 3:
 unsigned b=3;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt()->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt()->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt()->face_index_at_boundary(b,e);
      
   // Create new element and add to mesh
   Lagrange_multiplier_mesh_pt->add_element_pt(
    new ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(
     bulk_elem_pt,face_index));   
  }  

 
 // Loop over the elements in the Lagrange multiplier element mesh
 // for elements on the top boundary (boundary 3)
 n_element=Lagrange_multiplier_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a Lagrange multiplier element
   ImposeDisplacementByLagrangeMultiplierElement<ELEMENT> *el_pt = 
    dynamic_cast<ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>*>
    (Lagrange_multiplier_mesh_pt->element_pt(i));

   // Set the GeomObject that defines the boundary shape and set
   // which bulk boundary we are attached to(needed to extract
   // the boundary coordinate from the bulk nodes)
   el_pt->set_boundary_shape_geom_object_pt(Wall_geom_object_pt,b);

   // Loop over the nodes 
   unsigned nnod=el_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = el_pt->node_pt(j);
     
     // Is the node also on boundary 2 or 4?
     if ((nod_pt->is_on_boundary(2))||(nod_pt->is_on_boundary(4)))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       unsigned n_bulk_value=el_pt->nbulk_value(j);

       // The remaining ones are Lagrange multipliers and we pin them.
       unsigned nval=nod_pt->nvalue();
       for (unsigned j=n_bulk_value;j<nval;j++)
        {
         nod_pt->pin(j);
        }
      }
    }
  }
  
} // end of create_lagrange_multiplier_elements




//====start_of_delete_lagrange_multiplier_elements=======================
/// Delete elements that impose the prescribed boundary displacement
/// and wipe the associated mesh
//=======================================================================
template<class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::
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
 wall_mesh_pt()->assign_initial_values_impulsive();

} // end of set_initial_condition





//=========start_of_actions_before_adapt==================================
/// Actions before adapt: Wipe the mesh of prescribed traction elements
//========================================================================
template<class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::actions_before_adapt()
{
 // Kill the traction elements and wipe surface mesh
 delete_traction_elements(Applied_fluid_traction_mesh_pt);
 
 // Kill the  elements and wipe surface mesh
 delete_lagrange_multiplier_elements();
 
 // Rebuild the global mesh. 
 rebuild_global_mesh();

} // end of actions_before_adapt



//==========start_of_actions_after_adapt==================================
/// Actions after adapt: Rebuild the mesh of prescribed traction elements
//========================================================================
template<class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::actions_after_adapt()
{
 // Create prescribed-flux elements from all elements that are 
 // adjacent to boundary 5 and add them to surface mesh
 create_traction_elements(5,Bulk_mesh_pt,Applied_fluid_traction_mesh_pt);

 // Create the elements that impose the displacement constraint 
 // and attach them to the bulk elements that are
 // adjacent to the moving boundary
 create_lagrange_multiplier_elements();

 // Rebuild the global mesh
 rebuild_global_mesh();
 
 // Unpin all pressure dofs
 RefineableNavierStokesEquations<2>::
  unpin_all_pressure_dofs(Bulk_mesh_pt->element_pt());
 
 // Pin redundant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());
   
 // Loop over the traction elements to pass pointer to prescribed 
 // traction function
 unsigned n_element=Applied_fluid_traction_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to NavierStokesTractionElement element
   NavierStokesTractionElement<ELEMENT> *el_pt = 
    dynamic_cast<NavierStokesTractionElement<ELEMENT>*>(
     Applied_fluid_traction_mesh_pt->element_pt(e));
   
   // Set the pointer to the prescribed traction function
   el_pt->traction_fct_pt() = &Global_Physical_Variables::prescribed_traction;
  }

 // The functions used to update the no slip boundary conditions 
 // must be set on any new nodes that have been created during the 
 // mesh adaptation process. 
 // There is no mechanism by which auxiliary update functions 
 // are copied to newly created nodes.
 // (because, unlike boundary conditions, they don't occur exclusively 
 // at boundaries)
 
 // The velocity of the fluid nodes on the wall (fluid mesh boundary 3)
 // is set by the wall motion -- hence the no-slip condition needs to be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 unsigned ibound=3; 
 unsigned num_nod=bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   bulk_mesh_pt()->boundary_node_pt(ibound, inod)->
    set_auxiliary_node_update_fct_pt(
     FSI_functions::apply_no_slip_on_moving_wall);
  }
 
 // (Re-)setup fsi: Work out which fluid dofs affect wall elements
 // the correspondance between wall dofs and fluid elements is handled
 // during the remeshing, but the "reverse" association must be done
 // separately. We need to set up the interaction every time because the fluid
 // element adjacent to a given solid element's integration point may have 
 // changed.We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary is boundary 3 of 
 // the Fluid mesh.
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,3,Bulk_mesh_pt,Wall_mesh_pt);

} // end of actions_after_adapt



//============start_of_main====================================================
/// Driver code for a collapsible channel problem with FSI.
/// Presence of command line arguments indicates validation run with 
/// coarse resolution and small number of timesteps.
//=============================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // Reduction in resolution for validation run?
 unsigned coarsening_factor=4; 
 if (CommandLineArgs::Argc>1)
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
 
 // Set external pressure (on the wall stiffness scale). 
 Global_Physical_Variables::P_ext = 1.0e-1;
 
 // Pressure on the left boundary: This is consistent with steady
 // Poiseuille flow
 Global_Physical_Variables::P_up=12.0*(lup+lcollapsible+ldown);


#ifdef TAYLOR_HOOD

 // Build the problem with QtaylorHoodElements
 FSICollapsibleChannelProblem
  <RefineablePseudoSolidNodeUpdateElement<
  RefineableQTaylorHoodElement<2>,
  RefineableQPVDElement<2,3> > >
  problem(nup, ncollapsible, ndown, ny, 
          lup, lcollapsible, ldown, ly);

#else

 // Build the problem with QCrouzeixRaviartElements
 FSICollapsibleChannelProblem
  <RefineablePseudoSolidNodeUpdateElement<
  RefineableQCrouzeixRaviartElement<2>,
  RefineableQPVDElement<2,3> > >
  problem(nup, ncollapsible, ndown, ny, 
          lup, lcollapsible, ldown, ly);

#endif

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
 doc_info.set_directory("RESLT");
 
 // Open a trace file 
 ofstream trace_file;
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);

 // Output the initial condition
 problem.doc_solution(doc_info, trace_file);

 // Increment step number
 doc_info.number()++;

 // Find number of timesteps (reduced for validation)
 unsigned nstep = unsigned((t_max-t_min)/dt);
 if (CommandLineArgs::Argc>1)
  {
   nstep=3;
  }


 // Set targets for spatial adaptivity
 problem.bulk_mesh_pt()->max_permitted_error()=1.0e-3;
 problem.bulk_mesh_pt()->min_permitted_error()=1.0e-5;

 // Overwrite for validation run
 if (CommandLineArgs::Argc>1)
  {
   // Set targets for spatial adaptivity
   problem.bulk_mesh_pt()->max_permitted_error()=0.5e-2;
   problem.bulk_mesh_pt()->min_permitted_error()=0.5e-4;
  }

 // When performing the first timestep, we can adapt the mesh as many times
 // as we want because the initial condition can be re-set
 unsigned max_adapt=3; 
 bool first=true;

 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   // Solve the problem
   problem.unsteady_newton_solve(dt,max_adapt,first);
   
   // Outpt the solution
   problem.doc_solution(doc_info, trace_file);
   
   // Step number
   doc_info.number()++;

   // We've done the first step
   first=false;
   max_adapt=1; 
  }

 // Close trace file.
 trace_file.close();

}//end of main
  
