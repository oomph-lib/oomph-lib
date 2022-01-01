//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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

// The wall mesh
#include "meshes/one_d_lagrangian_mesh.h"

//Include the fluid mesh
#include "meshes/fsi_driven_cavity_mesh.h"


using namespace std;
using namespace oomph;



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







//====start_of_physical_parameters=====================
/// Namespace for phyical parameters
//======================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=100.0;

 /// Womersley = Reynolds times Strouhal
 double ReSt=100.0;
 
 /// Non-dimensional wall thickness. 
 double H=0.002;
 
 /// Fluid structure interaction parameter: Ratio of stresses used for
 /// non-dimensionalisation of fluid to solid stresses. 
 double Q=4.0e-5;

 /// Timescale ratio for solid (dependent parameter
 /// assigned in set_parameters())
 double Lambda_sq=2.0;

} // end of namespace




//====start_of_problem_class==========================================
/// Problem class
//====================================================================
template <class ELEMENT>
class FSIDrivenCavityProblem : public virtual Problem
{

 public :

 /// Constructor: The arguments are the number of elements,
 /// the lengths of the domain, the fractional height of the gap
 /// next to the moving lid and the period of the lid's oscillation 
 FSIDrivenCavityProblem(const unsigned& nx, 
                        const unsigned& ny,
                        const double& lx,
                        const double& ly,
                        const double& gap_fraction,
                        const double& period);
 
 /// Destructor
 ~FSIDrivenCavityProblem()
  { 
   // Mesh gets killed in general problem destructor
  }


 /// Access function for the specific bulk (fluid) mesh
 AlgebraicFSIDrivenCavityMesh<ELEMENT>* bulk_mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<
    AlgebraicFSIDrivenCavityMesh<ELEMENT>*>
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
  
 /// Update the velocity boundary condition on the moving lid
 void actions_before_implicit_timestep()
  {
   // Oscillating lid
   unsigned ibound=0;
   unsigned num_nod=bulk_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Which node are we dealing with?
     Node* node_pt=bulk_mesh_pt()->boundary_node_pt(ibound,inod);
     
     // Set velocity
     double veloc=1.0-cos(2.0*MathematicalConstants::Pi*time_pt()->time()/T);
     
     // Apply no slip
     node_pt->set_value(0,veloc);
    }
  }


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

protected : 
 
 /// Number of elements in the x direction 
 unsigned Nx;

 /// Number of elements in the y direction 
 unsigned Ny;

 /// Width of domain
 double Lx;

 /// Height of domain
 double Ly;

 /// Period of oscillation
 double T;

 /// Pointer to the "bulk" mesh
 AlgebraicFSIDrivenCavityMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "wall" mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* Wall_mesh_pt; 

 /// Pointer to the left control node
 Node* Left_node_pt;
 
 /// Pointer to right control node
 Node* Right_node_pt;
 
 /// Pointer to control node on the wall
 Node* Wall_node_pt;

 /// Pointer to geometric object (one Lagrangian, two Eulerian 
 /// coordinates) that will be built from the wall mesh
 MeshAsGeomObject* Wall_geom_object_pt;

};//end of problem class




//=====start_of_constructor======================================
/// Constructor for the collapsible channel problem
//===============================================================
template <class ELEMENT>
FSIDrivenCavityProblem<ELEMENT>::FSIDrivenCavityProblem(
 const unsigned& nx, 
 const unsigned& ny,
 const double& lx,
 const double& ly,
 const double& gap_fraction,
 const double& period)
{

 // Store problem parameters
 Nx=nx;
 Ny=ny;
 Lx=lx;
 Ly=ly;

 // Period of lid oscillation
 T=period;


 // Allow for crappy initial guess
 Max_newton_iterations = 20;
 Max_residuals = 1.0e8;

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 BDF<2>* fluid_time_stepper_pt=new BDF<2>;
 add_time_stepper_pt(fluid_time_stepper_pt);

 // Create the timestepper for the solid
 Newmark<2>* solid_time_stepper_pt=new Newmark<2>;
 //Steady<2>* solid_time_stepper_pt=new Steady<2>;
 add_time_stepper_pt(solid_time_stepper_pt);

 // Geometric object that represents the undeformed wall: 
 // A straight line at height y=ly; starting at x=0.
 UndeformedWall* undeformed_wall_pt=new UndeformedWall(0.0,ly);

 //Create the "wall" mesh with FSI Hermite elements
 Wall_mesh_pt = new OneDLagrangianMesh<FSIHermiteBeamElement>
  (nx,lx,undeformed_wall_pt,solid_time_stepper_pt);


 // Build a geometric object (one Lagrangian, two Eulerian coordinates)
 // from the wall mesh
 Wall_geom_object_pt=
  new MeshAsGeomObject(Wall_mesh_pt); 

 //Build bulk (fluid) mesh
 Bulk_mesh_pt =  new AlgebraicFSIDrivenCavityMesh<ELEMENT>
  (nx, ny, lx, ly, gap_fraction, Wall_geom_object_pt,fluid_time_stepper_pt);

 
 // Add the sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Wall_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();
   

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
   
  } // end loop over elements



 // Apply boundary conditions for fluid
 //------------------------------------

 // Pin the velocity on all boundaries apart from 1 and 5
 // (the gaps above the driven lid)
 for (unsigned ibound=0;ibound<6;ibound++)
  {
   if ((ibound!=1)&&(ibound!=5))
    {
     unsigned num_nod= bulk_mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       for(unsigned i=0;i<2;i++)
        {
         bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(i);
        }
      }
    }
  }//end of pin_velocity


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
   elem_pt->h_pt() = &Global_Physical_Variables::H;
    
   // Function that specifies the load ratios
   elem_pt->q_pt() = &Global_Physical_Variables::Q;

   // Set the undeformed shape for each element
   elem_pt->undeformed_beam_pt() = undeformed_wall_pt;

   // The normal on the wall elements as computed by the FSIHermiteElements
   // points away from the fluid rather than into the fluid (as assumed
   // by default)
   elem_pt->set_normal_pointing_out_of_fluid();

   // Timescale ratio for solid
   elem_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;

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
 unsigned ibound=5; 
 unsigned num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 unsigned control_nod=num_nod/2;
 Left_node_pt= bulk_mesh_pt()->boundary_node_pt(ibound, control_nod);
  
 // Right boundary
 ibound=1; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 control_nod=num_nod/2;
 Right_node_pt= bulk_mesh_pt()->boundary_node_pt(ibound, control_nod);
  
 
 // Set the pointer to the control node on the wall
 num_nod= wall_mesh_pt()->nnode();
 Wall_node_pt=wall_mesh_pt()->node_pt(nx/2);




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
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
  

}//end of constructor




//====start_of_doc_solution===================================================
/// Doc the solution
//============================================================================
template <class ELEMENT>
void FSIDrivenCavityProblem<ELEMENT>:: doc_solution(DocInfo& doc_info, 
                                                    ofstream& trace_file)
{ 

 // Doc fsi
 if (CommandLineArgs::Argc>1)
  {
   FSI_functions::doc_fsi<AlgebraicNode>(Bulk_mesh_pt,Wall_mesh_pt,doc_info);
  }

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
            << std::endl; 

} // end_of_doc_solution




//====start_of_apply_initial_condition========================================
/// Apply initial conditions
//============================================================================
template <class ELEMENT>
void FSIDrivenCavityProblem<ELEMENT>::set_initial_condition()
{ 
 // Impulsive start for wall
 wall_mesh_pt()->assign_initial_values_impulsive();

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
   
   // Assign initial condition: Zero flow
   bulk_mesh_pt()->node_pt(n)->set_value(0,0.0);
   bulk_mesh_pt()->node_pt(n)->set_value(1,0.0);
  } 

 // Assign initial values for an impulsive start
 bulk_mesh_pt()->assign_initial_values_impulsive();


} // end of set_initial_condition






