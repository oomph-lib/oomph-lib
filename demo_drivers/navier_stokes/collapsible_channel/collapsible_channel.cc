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
#include <iostream>

// Generic oomph-lib includes
#include "generic.h"

// Generic oomph-lib includes
#include "navier_stokes.h"

//Include the mesh
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





//===============start_of_oscillating_wall=================================
/// Straight, horizontal channel wall at \f$ y=H \f$ deforms into an 
/// oscillating parabola. The amplitude of the oscillation 
/// \f$ A \f$ and its period is \f$ T \f$.
/// The position vector to a point on the wall, parametrised by
/// the Lagrangian coordinate \f$ \zeta \in [0,L]\f$, is therefore given by
/// \f[ {\bf R}(\zeta,t) = 
///   \left(
///   \begin{array}{c}
///   L_{up} + \zeta  \\   1 
///   \end{array}
///   \right)
///   + A
///   \left(
///   \begin{array}{l}
///   - B \sin\left(\frac{2\pi}{L_{collapsible}}\zeta\right) \\ \left(
///   \frac{2}{L_{collapsible}}\right)^2 \zeta \ (L-\zeta)
///   \end{array}
///   \right)   
///   \ \sin\left(\frac{2\pi t}{T}\right)
///   \ {\cal R}(t)
///  \f]
/// The parameter \f$ B \f$ is zero by default. If it is set to a nonzero
/// value, the material particles on the wall also perform some 
/// horizontal motion. The "ramp" function 
/// \f[ 
/// {\cal R}(t) =  \left\{ 
/// \begin{array}{ll}  
/// \frac{1}{2}(1-\cos(\pi t/T)) & \mbox{for $t<T$} \\ 1 & \mbox{for $t \ge T$}
/// \end{array} 
/// \right.
/// \f]
/// provides a "smooth" startup of the oscillation. 
//=========================================================================
class OscillatingWall : public GeomObject
{

public:

 /// Constructor : It's a 2D object, parametrised by 
 /// one Lagrangian coordinate. Arguments: height at ends, x-coordinate of 
 /// left end, length, amplitude of deflection, period of oscillation, and
 /// pointer to time object
 OscillatingWall(const double& h, const double& x_left, const double& l, 
                 const double& a, const double& period, Time* time_pt) : 
  GeomObject(1,2), H(h), Length(l), X_left(x_left), A(a), B(0.0), T(period), 
  Time_pt(time_pt) 
  {}
 
 /// Destructor:  Empty
 ~OscillatingWall(){}

/// Access function to the amplitude
 double& amplitude(){return A;}

/// Access function to the period
 double& period(){return T;}

 /// Position vector at Lagrangian coordinate zeta 
 /// at time level t.
 void position(const unsigned& t, const Vector<double>&zeta, 
               Vector<double>& r) const
  {
   using namespace MathematicalConstants;

   // Smoothly ramp up the oscillation during the first period
   double ramp=1.0;
   if (Time_pt->time(t)<T)
    {
     ramp=0.5*(1.0-cos(Pi*Time_pt->time(t)/T));
    }
   
   // Position vector
   r[0] = zeta[0]+X_left 
    -B*A*sin(2.0*3.14159*zeta[0]/Length)*
    sin(2.0*Pi*(Time_pt->time(t))/T)*ramp;

   r[1] = H+A*((Length-zeta[0])*zeta[0])/pow(0.5*Length,2)*
    sin(2.0*Pi*(Time_pt->time(t))/T)*ramp;

  } // end of "unsteady" version


 /// "Current" position vector at Lagrangian coordinate zeta 
 void position(const Vector<double>&zeta, Vector<double>& r) const
  {
   position (0, zeta, r);
  }
 
 /// Number of geometric Data in GeomObject: None.
 unsigned ngeom_data() const {return 0;}              

private:
 
 /// Height at ends
 double H;

 /// Length
 double Length;

 /// x-coordinate of left end
 double X_left;

 /// Amplitude of oscillation
 double A;

 /// Relative amplitude of horizontal wall motion
 double B;

 /// Period of the oscillations
 double T;

 /// Pointer to the global time object
 Time* Time_pt;

}; // end of oscillating wall




//====start_of_Global_Physical_Variables================
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

 /// Traction required at the left boundary
 void prescribed_traction(const double& t,
                          const Vector<double>& x,
                          const Vector<double>& n,
                          Vector<double>& traction)
 {
  traction.resize(2);
  traction[0]=P_up;
  traction[1]=0.0;
 } 

} // end of Global_Physical_Variables



//=======start_of_problem_class=======================================
/// Problem class
//====================================================================
template <class ELEMENT>
class CollapsibleChannelProblem : public Problem
{

 public :

 /// Constructor : the arguments are the number of elements,
 /// the length of the domain and the amplitude and period of 
 /// the oscillations
 CollapsibleChannelProblem(const unsigned& nup, 
                           const unsigned& ncollapsible,
                           const unsigned& ndown,
                           const unsigned& ny,
                           const double& lup,
                           const double& lcollapsible, 
                           const double& ldown,
                           const double& ly,
                           const double& amplitude,
                           const double& period);
 
 /// Empty destructor
 ~CollapsibleChannelProblem() {} 
 
 /// Access function for the specific mesh
 RefineableCollapsibleChannelMesh<ELEMENT>* bulk_mesh_pt() 
  {

   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<RefineableCollapsibleChannelMesh<ELEMENT>*>
    (Bulk_mesh_pt);

  } // end of access to bulk mesh
 
 /// Update the problem specs before solve (empty) 
 void actions_before_newton_solve(){}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}
 
 /// Update the velocity boundary condition on the moving wall
 void actions_before_implicit_timestep();

 /// Actions before adapt: Wipe the mesh of prescribed traction elements
 void actions_before_adapt();
 
 /// Actions after adapt: Rebuild the mesh of prescribed traction elements
 void actions_after_adapt();

 /// Apply initial conditions
 void set_initial_condition();

 /// Doc the solution
 void doc_solution(DocInfo& doc_info, ofstream& trace_file);
 
private : 

 /// Create the prescribed traction elements on boundary b
 /// of the bulk mesh and stick them into the surface mesh.
 void create_traction_elements(const unsigned &b, 
                               Mesh* const &bulk_mesh_pt,
                               Mesh* const &surface_mesh_pt);
 
 /// Delete prescribed traction elements from the surface mesh
 void delete_traction_elements(Mesh* const &surface_mesh_pt);

 /// Number of elements in the x direction in the upstream part of the channel
 unsigned Nup;

 /// Number of elements in the x direction in the "collapsible" 
 /// part of the channel
 unsigned Ncollapsible;

 /// Number of elements in the x direction in the downstream part of the channel
 unsigned Ndown;

 /// Number of elements across the channel
 unsigned Ny;

 /// x-length in the upstream part of the channel
 double Lup;

 /// x-length in the "collapsible" part of the channel
 double Lcollapsible;

 /// x-length in the downstream part of the channel
 double Ldown;

 /// Transverse length
 double Ly;

 /// Pointer to the geometric object that parametrises the "collapsible" wall
 OscillatingWall* Wall_pt;
 
 /// Pointer to the "bulk" mesh
 RefineableCollapsibleChannelMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Pointer to the "surface" mesh that contains the applied traction
 /// elements
 Mesh* Surface_mesh_pt; 
 
 /// Pointer to the left control node
 Node* Left_node_pt;
 
 /// Pointer to right control node
 Node* Right_node_pt;
  
}; // end of problem class




//===start_of_constructor=======================================
/// Constructor for the collapsible channel problem
//===============================================================
template <class ELEMENT>
CollapsibleChannelProblem<ELEMENT>::CollapsibleChannelProblem(
 const unsigned& nup, 
 const unsigned& ncollapsible,
 const unsigned& ndown,
 const unsigned& ny,
 const double& lup,
 const double& lcollapsible, 
 const double& ldown,
 const double& ly,
 const double& amplitude,
 const double& period)
{
 // Number of elements
 Nup=nup;
 Ncollapsible=ncollapsible;
 Ndown=ndown;
 Ny=ny;

 // Lengths of domain
 Lup=lup;
 Lcollapsible=lcollapsible;
 Ldown=ldown;
 Ly=ly;

 // Overwrite maximum allowed residual to accomodate possibly
 // poor initial guess for solution
 Problem::Max_residuals=10000;

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 add_time_stepper_pt(new BDF<2>);

 // Parameters for wall object
 double height=ly;
 double length=lcollapsible;
 double x_left=lup;

 //Create the geometric object that represents the wall
 Wall_pt=new OscillatingWall(height, x_left, length, amplitude, period,
                             time_pt());

 //Build mesh
 Bulk_mesh_pt = new RefineableCollapsibleChannelMesh<ELEMENT>(
  nup, ncollapsible, ndown, ny,
  lup, lcollapsible, ldown, ly,
  Wall_pt,
  time_stepper_pt());


 // Enable boundary layer squash function?
#ifdef USE_BL_SQUASH_FCT

 // Set a non-trivial boundary-layer-squash function...
 Bulk_mesh_pt->bl_squash_fct_pt() = &BL_Squash::squash_fct; 

 // ... and update the nodal positions accordingly
 Bulk_mesh_pt->node_update();

#endif
 // end of boundary layer squash function


 // Create "surface mesh" that will contain only the prescribed-traction 
 // elements at the inflow. The default constructor just creates the mesh 
 // without giving it any elements, nodes, etc.
 Surface_mesh_pt = new Mesh;
 
 // Create prescribed-traction elements from all elements that are 
 // adjacent to boundary 5 (inflow boundary), and add them to the surface mesh.
 create_traction_elements(5,Bulk_mesh_pt,Surface_mesh_pt);

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes added so far into a single Mesh
 build_global_mesh();
   
 //Set errror estimator  for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 dynamic_cast<RefineableCollapsibleChannelMesh<ELEMENT>*>
  (Bulk_mesh_pt)->spatial_error_estimator_pt()=error_estimator_pt;
 
 
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
   
  } // end loop over bulk elements



  // Loop over the traction elements to pass pointer to prescribed 
  // traction function 
  unsigned n_el=Surface_mesh_pt->nelement();
  for(unsigned e=0;e<n_el;e++)
   {
    // Upcast from GeneralisedElement to NavierStokes traction element
    NavierStokesTractionElement<ELEMENT> *el_pt = 
     dynamic_cast< NavierStokesTractionElement<ELEMENT>*>(
      Surface_mesh_pt->element_pt(e));
    
    // Set the pointer to the prescribed traction function
    el_pt->traction_fct_pt() = &Global_Physical_Variables::prescribed_traction;

   }  // end loop over applied traction elements

  
  

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
 
 
  //x and y-velocities pinned along boundary 2, 3, 4 (top boundaries) :
 for(unsigned ib=2;ib<5;ib++)
  { 
   num_nod= bulk_mesh_pt()->nboundary_node(ib);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     for(unsigned i=0;i<2;i++)
      {
       bulk_mesh_pt()->boundary_node_pt(ib, inod)->pin(i);
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
   }// end of pin_velocity



  //Select control nodes "half way" up the inflow/outflow cross-sections
  //--------------------------------------------------------------------
  
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
  
  // Setup equation numbering scheme
  cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
  
} //end of constructor




//====start_of_doc_solution===================================================
/// Doc the solution
//============================================================================
template <class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>:: doc_solution(DocInfo& doc_info, 
                                                       ofstream& trace_file)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 bulk_mesh_pt()->output(some_file,npts);
 some_file.close();


 // Get the position of the midpoint on the geometric object
 Vector<double> zeta(1);
 zeta[0]=0.5*Lcollapsible;
 Vector<double> wall_point(2);
 Wall_pt->position(zeta,wall_point);
 
 // Write trace file
 trace_file << time_pt()->time() << " " 
            << wall_point[1]  << " "
            << Left_node_pt->value(0) << " "
            << Right_node_pt->value(0) << " "
            << std::endl;

 // Output wall shape
 sprintf(filename,"%s/wall%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned nplot=100; 
 for (unsigned i=0;i<nplot;i++)
  {
   zeta[0]=double(i)/double(nplot-1)*Lcollapsible;
   Wall_pt->position(zeta,wall_point);
   some_file << wall_point[0]  << " "
             << wall_point[1]  << std::endl;
  }
 some_file.close();

} // end_of_doc_solution




//==========start_of_create_traction_elements==================================
/// Create the traction elements
//============================================================================
template <class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::create_traction_elements(
 const unsigned &b, Mesh* const &bulk_mesh_pt, Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
    (bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e that lies along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-traction element
   NavierStokesTractionElement<ELEMENT>* traction_element_pt = 
    new  NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);
   
   //Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(traction_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_traction_elements


//============start_of_delete_traction_elements==============================
/// Delete traction elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::
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



//=======start_of_apply_initial_conditions===================================
/// Apply initial conditions: Impulsive start from steady Poiseuille flow
//============================================================================
template <class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::set_initial_condition()
{ 

 // Check that timestepper is from the BDF family
 if (time_stepper_pt()->type()!="BDF")
  {
   std::ostringstream error_stream;
   error_stream 
    << "Timestepper has to be from the BDF family!\n"
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



//=====start_of_actions_before_implicit_timestep==============================
/// Execute the actions before timestep: Update the velocity
/// boundary condition on the moving wall
//============================================================================
template <class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::actions_before_implicit_timestep()
{
 // Update the domain shape
 bulk_mesh_pt()->node_update();
 
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
} //end of actions_before_implicit_timestep




//=========start_of_actions_before_adapt==================================
/// Actions before adapt: Wipe the mesh of prescribed traction elements
//========================================================================
template<class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::actions_before_adapt()
{
 // Kill the traction elements and wipe surface mesh
 delete_traction_elements(Surface_mesh_pt);
 
 // Rebuild the global mesh. 
 rebuild_global_mesh();

} // end of actions_before_adapt


//==========start_of_actions_after_adapt==================================
/// Actions after adapt: Rebuild the mesh of prescribed traction elements
//========================================================================
template<class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::actions_after_adapt()
{
 // Create prescribed-flux elements from all elements that are 
 // adjacent to boundary 5 and add them to surface mesh
 create_traction_elements(5,Bulk_mesh_pt,Surface_mesh_pt);

 // Rebuild the global mesh
 rebuild_global_mesh();
 
 // Loop over the traction elements to pass pointer to prescribed traction function
 unsigned n_element=Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to NavierStokesTractionElement element
   NavierStokesTractionElement<ELEMENT> *el_pt = 
    dynamic_cast<NavierStokesTractionElement<ELEMENT>*>(
     Surface_mesh_pt->element_pt(e));
   
   // Set the pointer to the prescribed traction function
   el_pt->traction_fct_pt() = &Global_Physical_Variables::prescribed_traction;
  }
} // end of actions_after_adapt



//=======start_of_driver_code==================================================
/// Driver code for an unsteady adaptive collapsible channel problem
/// with prescribed wall motion. Presence of command line arguments
/// indicates validation run with coarse resolution and small number of
/// timesteps.
//=============================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // Reduction in resolution for validation run?
 unsigned coarsening_factor=1;
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

 // Initial amplitude of the wall deformation
 double amplitude=1.0e-2; // ADJUST 
  
 // Period of oscillation
 double period=0.45;

 // Pressure/applied traction on the left boundary: This is consistent with 
 // steady Poiseuille flow
 Global_Physical_Variables::P_up=12.0*(lup+lcollapsible+ldown);
 

 //Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 
 // Open a trace file 
 ofstream trace_file;
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);

 // Build the problem with Crouzeix Raviart Elements
 CollapsibleChannelProblem<RefineableQCrouzeixRaviartElement<2> > 
  problem(nup, ncollapsible, ndown, ny, 
          lup, lcollapsible, ldown, ly, 
          amplitude,period);
 

 // Number of timesteps per period
 unsigned nsteps_per_period=40;

 // Number of periods
 unsigned nperiod=3; 

 // Number of timesteps (reduced for validation)
 unsigned nstep=nsteps_per_period*nperiod;
 if (CommandLineArgs::Argc>1)
  {
   nstep=3;
  }
 
 //Timestep: 
 double dt=period/double(nsteps_per_period);

 // Start time
 double t_min=0.0;

 // Initialise timestep and set initial conditions
 problem.time_pt()->time()=t_min;
 problem.initialise_dt(dt);
 problem.set_initial_condition();
 
  // Output the initial solution
 problem.doc_solution(doc_info, trace_file);
 
 // Step number
 doc_info.number()++;


 // Set targets for spatial adaptivity
 problem.bulk_mesh_pt()->max_permitted_error()=1.0e-3;
 problem.bulk_mesh_pt()->min_permitted_error()=1.0e-5;

 // Overwrite with reduced targets for validation run to force
 // some refinement during the first few timesteps
 if (CommandLineArgs::Argc>1)
  {
   problem.bulk_mesh_pt()->max_permitted_error()=1.0e-4;
   problem.bulk_mesh_pt()->min_permitted_error()=1.0e-6;
  }


 // First timestep: We may re-assign the initial condition
 // following any mesh adaptation.
 bool first=true;

 // Max. number of adaptations during first timestep
 unsigned max_adapt=10;

 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   // Solve the problem
   problem.unsteady_newton_solve(dt, max_adapt, first);
    

   // Outpt the solution
   problem.doc_solution(doc_info, trace_file);
   
   // Step number
   doc_info.number()++;

   // We've done one step: Don't re-assign the initial conditions
   // and limit the number of adaptive mesh refinements to one
   // per timestep.
   first=false;
   max_adapt=1;
  }

 trace_file.close();
 
} //end of driver code

  

//  {
//   unsigned n=100;
//   for (unsigned i=0;i<n;i++)
//    {
//     double s=double(i)/double(n-1);
//     cout << s << " " <<  BL_Squash::squash_fct(s) << std::endl;
//    }
//   pause("done");
//  }
