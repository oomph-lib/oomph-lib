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
#include <cstring>

// Generic oomph-lib includes
#include "generic.h"

// Generic oomph-lib includes
#include "navier_stokes.h"

//Include the mesh
#include "meshes/collapsible_channel_mesh.h"

//Include Womersley stuff
#include "womersley.h"

//Include FSI preconditioner
#include "multi_physics.h"

//#include "flux_control_elements.h"
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
  {
   Remain_steady_at_maximum_amplitude=false;
  }
 
 /// Destructor:  Empty
 ~OscillatingWall(){}

 /// Access function to the amplitude
 double& amplitude(){return A;}

 /// Access function to the period
 double& period(){return T;}
 
 /// Set the flag so that the wall remains steady 
 /// once it reaches maximum amplitude
 void enable_remain_steady_at_maximum_amplitude()
  {Remain_steady_at_maximum_amplitude = true;}

 /// Set the flag so that the wall continues to oscillate after it
 /// reaches maximum amplitude
 void disable_remain_steady_at_maximum_amplitude()
  {Remain_steady_at_maximum_amplitude = false;}

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

   if (Remain_steady_at_maximum_amplitude)
    {
     // Position vector
     r[0] = zeta[0]+X_left 
      +B*A*sin(2.0*3.14159*zeta[0]/Length)*ramp;
     
     r[1] = H-A*((Length-zeta[0])*zeta[0])/pow(0.5*Length,2)*ramp;
    }
   else
    {
     // Position vector
     r[0] = zeta[0]+X_left 
      -B*A*sin(2.0*3.14159*zeta[0]/Length)*
      sin(2.0*Pi*(Time_pt->time(t))/T)*ramp;
     
     r[1] = H+A*((Length-zeta[0])*zeta[0])/pow(0.5*Length,2)*
      sin(2.0*Pi*(Time_pt->time(t))/T)*ramp;
    }

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

 /// Flag which if true causes the wall to move to it's maximum amplitude
 /// position over one period and then remains steady
 bool Remain_steady_at_maximum_amplitude;

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
 
 /// Default pressure on the left boundary when controlling
 /// flow by a pressure drop
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


 /// Pointer to Data holding downstream pressure load 
 Data* Pout_data_pt;

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
                           const double& limpedance,
                           const double& amplitude,
                           const double& period,
                           const unsigned& outflow);
 
 /// Empty destructor
 ~CollapsibleChannelProblem() {} 
 
 /// Access function for the fluid bulk mesh
 RefineableCollapsibleChannelMesh<ELEMENT>* bulk_mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<RefineableCollapsibleChannelMesh<ELEMENT>*>
    (Bulk_mesh_pt);

  } // end of access to bulk mesh

 /// Access function to the pointer to the master impedance mesh which
 /// is used when using impedance boundary conditions with flux control elements
 Mesh* &outflow_impedance_master_mesh_pt()
  {
   return Outflow_impedance_master_mesh_pt;
  }
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info, ofstream& trace_file);

 /// Run unsteady problem
 void unsteady_run(string directory_for_data, double nstep, const bool& validation_run=false);

 /// Refine any element whoes centre node lies between x coorodinates
 /// x_min and x_max
 void refine_elements_based_on_x_coord(const double& x_min, const double& x_max);


 /// Function to switch off wall oscillations
 void switch_off_wall_oscillations() 
  {
   Wall_pt->enable_remain_steady_at_maximum_amplitude();
  }

 /// Create a mesh for the NavierStokesSchurComplementPreconditioner
 Mesh* create_mesh_for_navier_stokes_preconditioner();

protected:


 /// Apply Poiseuille flow on outlet
 void set_poiseuille_outflow();

 /// Apply parallel flow on outlet
 void set_parallel_outflow();
 
 /// Update the problem specs before solve
 void actions_before_newton_solve();

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

 
private : 

 
 /// Delete flux control elements from outflow boundary
 void delete_outflow_flux_control_elements();

 /// Attach flux control elements to outflow boundary 
 void setup_outflow_flux_control_elements();

 /// Create the prescribed traction elements on the inflow
 void setup_inflow_traction_elements();
 
 /// Delete prescribed traction elements from the surface mesh
 void delete_inflow_traction_elements();

 // Delete impedance elements for outflow boundary
 void delete_outflow_impedance_elements();

 // Set up impedance elements for outflow boundary
 void setup_outflow_impedance_elements();

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
 /// elements at the inflow
 Mesh* Inflow_traction_mesh_pt; 

 /// Pointer to mesh containing the net flux control element (only!)
 Mesh* Outflow_flux_control_master_mesh_pt;

 /// Pointer to mesh of flux control elements at outflow
 Mesh* Outflow_flux_control_sub_mesh_pt;

 /// Pointer to impendance tube
 WomersleyOutflowImpedanceTube<QWomersleyElement<1,3>,1>* 
 Womersley_impedance_tube_pt;

 /// Pointer to the impedance mesh
 Mesh* Outflow_impedance_mesh_pt;

 /// Pointer to the master impedance mesh when using flux control elements
 Mesh* Outflow_impedance_master_mesh_pt;
 
 /// Pointer to the left control node
 Node* Left_node_pt;
 
 /// Pointer to right control node
 Node* Right_node_pt;

 /// Length of optional impedance tube
 double L_impedance;

 /// Prescribed volume flux: set to -1 for value consistent with Re
 double Prescribed_volume_flux;

 /// Period of oscillation
 double Period;

 /// Amplitude of osciallation
 double Amplitude;

 /// Outflow boundary type:
 /// 0: prescribed flux boundary conditions
 /// 1: impedance boundary conditions
 /// 2: prescribed outflow
 /// Otherwise we have pressure driven flow with zero pressure at outlet
 unsigned Outflow_type;

 /// Id of the fluid bulk mesh
 unsigned Bulk_mesh_id;

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
 const double& limpedance,
 const double& amplitude,
 const double& period,
 const unsigned& outflow)
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

 // Other input data
 Period=period;
 Amplitude=amplitude;

 // Outflow boundary type
 // 0: prescribed flux boundary conditions
 // 1: impedance boundary conditions
 // 2: impedance boundary conditions using flux control elements 
 // 3: prescribed outflow
 // Otherwise we have pressure driven flow with zero pressure at outlet
 Outflow_type = outflow;

 /// Default length of optional impedance tube
 L_impedance=limpedance;
 
 // Overwrite maximum allowed residual to accomodate possibly
 // poor initial guess for solution
 Problem::Max_residuals=100000000;

 // Point optional mesh pointers and other stuff to null 
 Inflow_traction_mesh_pt=0;
 Outflow_flux_control_master_mesh_pt=0;
 Outflow_flux_control_sub_mesh_pt=0;
 Outflow_impedance_mesh_pt=0;
 Outflow_impedance_master_mesh_pt=0;
 Womersley_impedance_tube_pt=0;

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

 // Add the bulk sub mesh to the problem
 Bulk_mesh_id=add_sub_mesh(Bulk_mesh_pt);

 // Create "surface mesh" that will contain only the prescribed-traction 
 // elements at the inflow. The default constructor just creates the mesh 
 // without giving it any elements, nodes, etc.
 Inflow_traction_mesh_pt = new Mesh;
 
 // Create prescribed-traction elements 
 setup_inflow_traction_elements();
 
 // Add the new sub mesh to the problem
 add_sub_mesh(Inflow_traction_mesh_pt);

 if (Outflow_type==0)
  {
   oomph_info << "Flux control flow\n";
   
   // Create the flux control meshes
   
   // Construct Outflow_flux_control_master_mesh_pt
   Outflow_flux_control_master_mesh_pt = new Mesh;
   
   // Construct Outflow_flux_control_sub_mesh_pt
   Outflow_flux_control_sub_mesh_pt = new Mesh;
   
   // Set up elements for sub meshes
   setup_outflow_flux_control_elements();
   
   // Add new sub meshes
   add_sub_mesh(Outflow_flux_control_master_mesh_pt);
   add_sub_mesh(Outflow_flux_control_sub_mesh_pt);
  }
 else if (Outflow_type==1)
  {
   oomph_info << "Pressure driven flow with impedance tube of length "
              << L_impedance << "\n";
   
   if (L_impedance>0.0)
    {
     // Construct the mesh
     Outflow_impedance_mesh_pt = new Mesh;
     
     // Set up elements and impedance tube
     setup_outflow_impedance_elements();
     
     // Add new sub meshes
     add_sub_mesh(Outflow_impedance_mesh_pt);
    }
  }
 else if (Outflow_type==2)
  {
   oomph_info << "Pressure driven flow with impedance tube of length "
              << L_impedance << " using flux control elements\n";
   
   if (L_impedance>0.0)
    {
     // Construct Outflow_flux_control_sub_mesh_pt
     Outflow_flux_control_sub_mesh_pt = new Mesh;
     
     // Construct Mesh
     Outflow_impedance_master_mesh_pt = new Mesh;
     
     // Set up elements and impedance tube
     setup_outflow_impedance_elements();
     
     // Add new sub meshes
     add_sub_mesh(Outflow_flux_control_sub_mesh_pt);
     add_sub_mesh(Outflow_impedance_master_mesh_pt);
    }
  }
 else if (Outflow_type==3)
  {
   oomph_info << "Prescribed outflow\n";
  }
 else
  {
   oomph_info << "Pressure driven flow\n";
  }
 
 // Combine all submeshes added so far into a single Mesh
 build_global_mesh();
   
 //Set errror estimator for bulk mesh
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

 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());

 // Loop over the traction elements to pass pointer to prescribed 
 // traction function 
 if (Inflow_traction_mesh_pt!=0)
  {
   unsigned n_el=Inflow_traction_mesh_pt->nelement();
   for(unsigned e=0;e<n_el;e++)
    {
     // Upcast from GeneralisedElement to NavierStokes traction element
     NavierStokesTractionElement<ELEMENT> *el_pt = 
      dynamic_cast< NavierStokesTractionElement<ELEMENT>*>(
       Inflow_traction_mesh_pt->element_pt(e));
     
     // Set the pointer to the prescribed traction function
     el_pt->traction_fct_pt() = 
      &Global_Physical_Variables::prescribed_traction;
     
    }  // end loop over applied traction elements
  }
  
  

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
 
 //y-velocity pinned along boundary 1 (right boundary) and also x-velocity for 
 //prescribed outflow
 ibound=1; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(1);
   if (Outflow_type==3)
    {
     bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(0);
    }
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
  oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} //end of constructor



//============================================================================
/// Update the inflow and outflow boundary conditioner before solve
//============================================================================
template <class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::actions_before_newton_solve()
{
 Bulk_mesh_pt->node_update();

 // For prescribed outflow set Poiseuille flow on the outlet
 if (Outflow_type==3)
  {
   set_poiseuille_outflow();
  }
 // Otherwise parallel flow on the outlet
 else
  {
   set_parallel_outflow();
  }
}


//============================================================================
/// Refine any element whoes centre node lies between x coorodinates
/// x_min and x_max
//============================================================================
template <class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::refine_elements_based_on_x_coord(
 const double& x_min,
 const double& x_max)
  {
   // Vector to store pointers to elements to be refined
   Vector<RefineableElement*> element_pt;
   
   // Get elements we want to refine
   unsigned n_el = Bulk_mesh_pt->nelement();
   
   for (unsigned e=0; e<n_el; e++)
    {
     // Get element
     RefineableElement* el_pt = 
      dynamic_cast<RefineableElement*>(Bulk_mesh_pt->element_pt(e));
     
     // Get x-coord of node at centre
     double x_centre=el_pt->node_pt(4)->x(0);

     // Add to Vector if to be refined
     if ((x_centre>=x_min) && (x_centre<=x_max))
      {
       element_pt.push_back(el_pt);
      }     
    }
    
   // Refine the elements
   oomph_info << "Number of elements to refine between x=" << x_min
              << " & " << x_max << " is " 
              << element_pt.size() << "\n";
   refine_selected_elements(Bulk_mesh_id, element_pt);
  }

//====start_of_doc_solution===================================================
/// Doc the solution
//============================================================================
template <class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::doc_solution(DocInfo& doc_info, 
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




//===========================================================================
/// Create the traction elements
//============================================================================
template <class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::setup_inflow_traction_elements()
{
 // inflow is boundary 5
 unsigned b=5;
 
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
    (Bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e that lies along boundary b
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-traction element
   NavierStokesTractionElement<ELEMENT>* traction_element_pt = 
    new  NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);
   
   //Add the prescribed-flux element to the surface mesh
   Inflow_traction_mesh_pt->add_element_pt(traction_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} 


//============start_of_delete_traction_elements==============================
/// Delete traction elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::delete_inflow_traction_elements()
{
 // How many surface elements are in the surface mesh
 unsigned n_element = Inflow_traction_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill element
   delete Inflow_traction_mesh_pt->element_pt(e);
  }

 // Wipe the mesh
 Inflow_traction_mesh_pt->flush_element_and_node_storage();

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
 
 // Shift timesteps in impedance tube problem
 if (Womersley_impedance_tube_pt!=0)
  {
   Womersley_impedance_tube_pt->shift_time_values(
    time_stepper_pt()->time_pt()->dt());
  }

} //end of actions_before_implicit_timestep




//=========start_of_actions_before_adapt==================================
/// Actions before adapt: Wipe the mesh of prescribed traction elements
//========================================================================
template<class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::actions_before_adapt()
{
 // Kill the traction elements and wipe surface mesh
 if (Inflow_traction_mesh_pt!=0)
  {
   delete_inflow_traction_elements();
  }
 
 // Delete the flux control elements if required
 if (Outflow_flux_control_sub_mesh_pt!=0)
  {
   if (Outflow_type==0)
    {
     delete_outflow_flux_control_elements();
    }
   else if (Outflow_type==2)
    {
     delete_outflow_impedance_elements();
    }
  }
 
 // Delete impedance elements if required
 if (Outflow_impedance_mesh_pt!=0)
  {
   delete_outflow_impedance_elements();
  }
 
 // Rebuild the global mesh. 
 rebuild_global_mesh();
 
} // end of actions_before_adapt


//==========start_of_actions_after_adapt==================================
/// Actions after adapt: Rebuild the mesh of prescribed traction elements
//========================================================================
template<class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::actions_after_adapt()
{
 // Create prescribed-flux elements on inflow
 if (Inflow_traction_mesh_pt!=0)
  {
   setup_inflow_traction_elements();
  }

 // Attach flux control elements to their mesh if required
 if (Outflow_flux_control_sub_mesh_pt!=0)
  {
   if (Outflow_type==0)
    {
     setup_outflow_flux_control_elements();
    }
   else if (Outflow_type==2)
    {
     setup_outflow_impedance_elements();
    }
  }
 
 // Attach impedance elements to their mesh
 if (Outflow_impedance_mesh_pt!=0)
  {
   setup_outflow_impedance_elements();
  }
 
 // Rebuild the global mesh
 rebuild_global_mesh();
 
 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());
 
 // Loop over the traction elements to pass pointer to prescribed 
 // traction function
 if (Inflow_traction_mesh_pt!=0)
  {
   unsigned n_element=Inflow_traction_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to NavierStokesTractionElement element
     NavierStokesTractionElement<ELEMENT> *el_pt = 
      dynamic_cast<NavierStokesTractionElement<ELEMENT>*>(
       Inflow_traction_mesh_pt->element_pt(e));
   
     // Set the pointer to the prescribed traction function
     el_pt->traction_fct_pt() = 
      &Global_Physical_Variables::prescribed_traction;
    }
  }
} // end of actions_after_adapt


//========================================================================
/// Create a mesh for the NavierStokesSchurComplementPreconditioner
//========================================================================
template<class ELEMENT>
Mesh* CollapsibleChannelProblem<ELEMENT>::
create_mesh_for_navier_stokes_preconditioner()
{
 //Vector to hold the meshes
 Vector<Mesh*> meshes;
 
 //Add the bulk mesh
 meshes.push_back(Bulk_mesh_pt);
 
 // Add the flux control master mesh if it exists
 if (Outflow_flux_control_master_mesh_pt!=0)
  {
   meshes.push_back(Outflow_flux_control_master_mesh_pt);
  }

 // Build "combined" mesh from vector of submeshes
 Mesh* mesh_pt = new Mesh(meshes);

 return mesh_pt;
}


//========================================================================
/// Delete flux control elements from outflow boundary
//========================================================================
template<class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::delete_outflow_flux_control_elements()
{
 // Delete the master flux control element
 delete Outflow_flux_control_master_mesh_pt->element_pt(0);
 
 // Wipe the mesh
 Outflow_flux_control_master_mesh_pt->flush_element_and_node_storage();
 
 // Loop over all flux control elements at outflow
 unsigned n_bound_el=Outflow_flux_control_sub_mesh_pt->nelement();
 for (unsigned e=0;e<n_bound_el;e++)
  { 
   //Delete the element
   delete Outflow_flux_control_sub_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 Outflow_flux_control_sub_mesh_pt->flush_element_and_node_storage();
}


//========================================================================
/// Attach flux control elements to outflow boundary 
//========================================================================
template<class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::setup_outflow_flux_control_elements()
{
 // Attach elements to apply traction on boundary 1
 unsigned b = 1;
 
 // Loop over all elements on this boundary
 unsigned n_bound_el=Bulk_mesh_pt->nboundary_element(b);
 for (unsigned e=0;e<n_bound_el;e++)
  { 
   // Get pointer to bulk element
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
      
   // Get the index of the face of element e along boundary b
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
      
   // Build the corresponding flux control element
   NavierStokesFluxControlElement<ELEMENT>* flux_element_pt = new 
    NavierStokesFluxControlElement<ELEMENT>(bulk_elem_pt, face_index);

   //Add the new element to its mesh
   Outflow_flux_control_sub_mesh_pt->add_element_pt(flux_element_pt);     
  }
  
    
 // Build master element
 NetFluxControlElement* flux_control_el_pt = new NetFluxControlElement(
  Outflow_flux_control_sub_mesh_pt,
  &Prescribed_volume_flux);
 
 // Add NetFluxControlElement to its mesh
 Outflow_flux_control_master_mesh_pt->
  add_element_pt(flux_control_el_pt);
 
 // Get pointer to the outflow pressure data
 Global_Physical_Variables::Pout_data_pt = 
  flux_control_el_pt->pressure_data_pt();
 
 // Loop over the elements in the sub mesh and pass 
 // Global_Physical_Variables::Pout_data_pt to the element
 unsigned n_el = Outflow_flux_control_sub_mesh_pt->nelement();
 for (unsigned e=0; e< n_el; e++)
  {
   // Get pointer to the element
   GeneralisedElement* el_pt = 
    Outflow_flux_control_sub_mesh_pt->element_pt(e);
   
   // Dynamic cast
   dynamic_cast<NavierStokesFluxControlElement<ELEMENT>* >(el_pt)->   
    add_pressure_data(Global_Physical_Variables::Pout_data_pt);
  }

}


//========================================================================
/// Delete impedance elements from outflow boundary
//========================================================================
template<class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::delete_outflow_impedance_elements()
{
 if (Outflow_type==1)
  {
   // Loop over all elements at outflow
   unsigned n_bound_el=Outflow_impedance_mesh_pt->nelement();
   for (unsigned e=0;e<n_bound_el;e++)
    { 
     //Delete the element
     delete Outflow_impedance_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Outflow_impedance_mesh_pt->flush_element_and_node_storage();

   // Delete the impedance tube
   delete Womersley_impedance_tube_pt;
   Womersley_impedance_tube_pt=0; 
  }
 else 
  {
   // Delete the elements in the master mesh
   delete Outflow_impedance_master_mesh_pt->element_pt(0);
   delete Outflow_impedance_master_mesh_pt->element_pt(1);
 
   // Wipe the mesh
   Outflow_impedance_master_mesh_pt->flush_element_and_node_storage();
 
   // Loop over all flux control elements at outflow
   unsigned n_bound_el=Outflow_flux_control_sub_mesh_pt->nelement();
   for (unsigned e=0;e<n_bound_el;e++)
    { 
     //Delete the element
     delete Outflow_flux_control_sub_mesh_pt->element_pt(e);
    }
 
   // Wipe the mesh
   Outflow_flux_control_sub_mesh_pt->flush_element_and_node_storage();

   // Delete the impedance tube
   delete Womersley_impedance_tube_pt;
   Womersley_impedance_tube_pt=0; 
  }
  
}

//========================================================================
/// Attach impedance elements to outflow boundary 
//========================================================================
template<class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::setup_outflow_impedance_elements()
{
 // Check we are at time t=0, otherwise this function will not work 
 // since an ImpedanceTube is contstructed so it's been at steady state
 // for all previous time.
 if (time_stepper_pt()->time()>0.0)
  {
   std::ostringstream error_stream;
   error_stream 
    << "Impedance is set up assuming problem has been at rest for all previous "
    << "times, \ndid you mean to set up the impedance elements at time="
    << time_stepper_pt()->time() << " ?\n";
   
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Elements apply traction on boundary 1
 unsigned b = 1;
 
 // Loop over all elements on this boundary
 unsigned n_bound_el=Bulk_mesh_pt->nboundary_element(b);

 if (Outflow_type==1)
  {
   for (unsigned e=0;e<n_bound_el;e++)
    { 
     // Get pointer to bulk element
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));
     
     // Get the index of the face of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
   
     // Build the corresponding prescribed traction element
     NavierStokesImpedanceTractionElement
      <ELEMENT,QWomersleyElement<1,3>,1>* traction_element_pt =
      new NavierStokesImpedanceTractionElement
      <ELEMENT,QWomersleyElement<1,3>,1>(bulk_elem_pt,face_index);
        
     //Add the prescribed traction element to the mesh
     Outflow_impedance_mesh_pt->add_element_pt(traction_element_pt);
    }
 
   // fixed coord in bulk mesh which is constant at the outflow
   unsigned fixed_coordinate=0;
 
   // nodal index of velocity component which defines the outflow
   unsigned v_index=0;
 
   // Build Womersley impedance tube and let it provide the
   // outflow traction to the elements in the Outflow_traction_mesh_pt.
   Womersley_impedance_tube_pt=
    new WomersleyOutflowImpedanceTube<QWomersleyElement<1,3>,1>
    (L_impedance,
     Outflow_impedance_mesh_pt,
     fixed_coordinate,
     v_index);    
  }
 else
  {
   // Loop over all elements on this boundary
   unsigned n_bound_el=Bulk_mesh_pt->nboundary_element(b);
   for (unsigned e=0;e<n_bound_el;e++)
    { 
     // Get pointer to bulk element
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));
     
     // Get the index of the face of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
     
     // Build the corresponding flux control element
     NavierStokesFluxControlElement<ELEMENT>* flux_element_pt = new 
      NavierStokesFluxControlElement<ELEMENT>(bulk_elem_pt, face_index);
     
     //Add the new element to its mesh
     Outflow_flux_control_sub_mesh_pt->add_element_pt(flux_element_pt);     
    }
   
   // fixed coord in bulk mesh which is constant at the outflow
   unsigned fixed_coordinate=0;
 
   // nodal index of velocity component which defines the outflow
   unsigned v_index=0;
   
   // Build Womersley impedance tube and let it provide the
   // outflow traction to the elements in the Outflow_traction_mesh_pt.
   Womersley_impedance_tube_pt=
    new WomersleyOutflowImpedanceTube<QWomersleyElement<1,3>,1>(
     L_impedance,
     Outflow_flux_control_sub_mesh_pt,
     fixed_coordinate,
     v_index);    
   
   // Build the pressure control element and add to the master mesh
   NavierStokesWomersleyPressureControlElement* pressure_control_element
   = new NavierStokesWomersleyPressureControlElement(
    Womersley_impedance_tube_pt);
   
   Outflow_impedance_master_mesh_pt->add_element_pt(pressure_control_element);
   
   // Build the flux control master element and add to the master mesh
   NetFluxControlElementForWomersleyPressureControl* flux_control_element
   = new NetFluxControlElementForWomersleyPressureControl
      (Outflow_flux_control_sub_mesh_pt, pressure_control_element);
   
   Outflow_impedance_master_mesh_pt->add_element_pt(flux_control_element);
  }
}


//============================================================================
/// Apply Poiseuille flow on outlet
//============================================================================
template <class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::set_poiseuille_outflow()
{ 
 // Outflow is boundary 1
 unsigned b = 1;
 
 // Loop over all elements on this boundary
 unsigned n_bound_el=Bulk_mesh_pt->nboundary_element(b);
 for (unsigned e=0;e<n_bound_el;e++)
  { 
   // Get nodal coordinates
   Vector<double> x(2);
   x[0]=bulk_mesh_pt()->node_pt(e)->x(0);
   x[1]=bulk_mesh_pt()->node_pt(e)->x(1);
   
   // Assign Poiseuille flow
   bulk_mesh_pt()->node_pt(e)->set_value(0,6.0*(x[1]/Ly)*(1.0-(x[1]/Ly)));
   bulk_mesh_pt()->node_pt(e)->set_value(1,0.0);
  } 

}

//============================================================================
/// Apply parallel flow on outlet
//============================================================================
template <class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::set_parallel_outflow()
{ 
 // Outflow is boundary 1
 unsigned b = 1;
 
 // Loop over all elements on this boundary
 unsigned n_bound_el=Bulk_mesh_pt->nboundary_element(b);
 for (unsigned e=0;e<n_bound_el;e++)
  { 
   // Assign parallel flow
   bulk_mesh_pt()->node_pt(e)->set_value(1,0.0);
  } 
}

//============================================================================
/// Run unsteady problem
//============================================================================
template <class ELEMENT>
void CollapsibleChannelProblem<ELEMENT>::
unsteady_run(string directory_for_data, double nstep, const bool& validation_run)
{ 
 // Set volume flux consistent with Re
 Prescribed_volume_flux = 1.0;

 // Pressure/applied traction on the left boundary: This is consistent with 
 // steady Poiseuille flow
 Global_Physical_Variables::P_up=
  12.0*(Lup+Lcollapsible+Ldown+L_impedance);

 DocInfo doc_info;
 doc_info.set_directory(directory_for_data);

 // Open a trace file 
 ofstream trace_file;
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);
 
 // Number of timesteps per period
 unsigned nsteps_per_period=40;
 
 //Timestep: 
 double dt=Period/double(nsteps_per_period);

 if (validation_run)
  {
   nstep=10;
   dt=0.25;
  }

 oomph_info << "Timestep=" << dt << "\n";

 // Start time
 double t_min=0.0;

 // Initialise timestep and set initial conditions
 time_pt()->time()=t_min;
 initialise_dt(dt);
 set_initial_condition();
 
 // Output the initial solution
 doc_solution(doc_info, trace_file);
 
 // Step number
 doc_info.number()++;
  
 // Set targets for spatial adaptivity
 bulk_mesh_pt()->max_permitted_error()=1.0e-3;
 bulk_mesh_pt()->min_permitted_error()=1.0e-5;

 // Setup impedance tube so its initial conditions 
 // correspond to an impulsive start from the current flow rate
 if (Womersley_impedance_tube_pt!=0)
  {
   double q_initial=Womersley_impedance_tube_pt->
    total_volume_flux_into_impedance_tube();
   oomph_info << "q_initial=" << q_initial << "\n";
   Womersley_impedance_tube_pt->setup(&Global_Physical_Variables::ReSt,
                                      dt,
                                      q_initial,
                                      new BDF<2>); 
  }

 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   oomph_info << "\n\nNewton solve: " << istep << "\n";
   oomph_info << "-----------------\n";
  
   // Solve the problem
   unsteady_newton_solve(dt);
   
   // Outpt the solution
   doc_solution(doc_info, trace_file);
   
   // Step number
   doc_info.number()++;
  }

 trace_file.close();

}

//=======start_of_driver_code==================================================
/// Driver code for an unsteady adaptive collapsible channel problem
/// with prescribed wall motion. Presence of command line arguments
/// indicates validation run with coarse resolution and small number of
/// timesteps.
//=============================================================================
int main(int argc, char *argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 //Uncomment for parallel stuff
 /*#ifdef OOMPH_HAS_MPI
 // Initialise MPI
 MPI_Helpers::init(argc,argv);

 // Swtich off output modifier
 oomph_info.output_modifier_pt() = &default_output_modifier;
 
 // switch off oomph_info output for all processors but rank 0
 if (MPI_Helpers::My_rank!=0)
  {
   oomph_info.stream_pt() = &oomph_nullstream;
   OomphLibWarning::set_stream_pt(&oomph_nullstream);
   OomphLibError::set_stream_pt(&oomph_nullstream);
  }
 else
  {
   oomph_info << "\n\n=====================================================\n";
   oomph_info << "Number of processors: " 
              << MPI_Helpers::Nproc << "\n";
  }
  #endif */

 // Set default values
 string directory_for_data = "RESLT";
 unsigned outflow = 0;
 Global_Physical_Variables::Re = 10;
 unsigned solver=1;
 unsigned p_solver=2;
 unsigned f_solver=2;
 unsigned refinement=0;
 unsigned nsteps=200;
 double amplitude=1.0e-2;
 double period=0.45;
 bool validation_run=false;

 //Set Trilinos default value
#ifdef OOMPH_HAS_TRILINOS 
 unsigned f_ml_settings=1;
#endif
 
 // Set Hypre default values
#ifdef OOMPH_HAS_HYPRE
 unsigned f_bamg_smoother=1;
 double f_bamg_damping=1.0;
 double f_bamg_strength=0.25;
#endif


 
 // Parse command line
 int arg_index = 0;
 unsigned print_help = 0;

 while (arg_index < argc)
  {
   if ( strcmp(argv[arg_index], "-dir") == 0 )
    {
     arg_index++;
     string directory_number = argv[arg_index++];
     directory_for_data += directory_number;
    }
   else if ( strcmp(argv[arg_index], "-outflow") == 0 )
    {
     arg_index++;
     outflow = atoi(argv[arg_index++]);
    }
   else if ( strcmp(argv[arg_index], "-Re") == 0 )
    {
     arg_index++;
     Global_Physical_Variables::Re = atof(argv[arg_index++]);
    }
   else if ( strcmp(argv[arg_index], "-solver") == 0 )
    {
     arg_index++;
     solver = atoi(argv[arg_index++]);
    }
   else if ( strcmp(argv[arg_index], "-p_solver") == 0 )
    {
     arg_index++;
     p_solver = atoi(argv[arg_index++]);
    }
   else if ( strcmp(argv[arg_index], "-f_solver") == 0 )
    {
     arg_index++;
     f_solver = atoi(argv[arg_index++]);
    }
#ifdef OOMPH_HAS_HYPRE
   else if ( strcmp(argv[arg_index], "-f_bamg_smoother") == 0 )
    {
     arg_index++;
     f_bamg_smoother = atoi(argv[arg_index++]);
    }
   else if ( strcmp(argv[arg_index], "-f_bamg_damping") == 0 )
    {
     arg_index++;
     f_bamg_damping = atof(argv[arg_index++]);
    }  
   else if ( strcmp(argv[arg_index], "-f_bamg_strength") == 0 )
    {
     arg_index++;
     f_bamg_strength = atof(argv[arg_index++]);
    }
#endif
#ifdef OOMPH_HAS_TRILINOS
   else if ( strcmp(argv[arg_index], "-f_ml_settings") == 0 )
    {
     arg_index++;
     f_ml_settings = atoi(argv[arg_index++]);
    }
#endif
   else if ( strcmp(argv[arg_index], "-refine") == 0 )
    {
     arg_index++;
     refinement = atoi(argv[arg_index++]);
    }
   else if ( strcmp(argv[arg_index], "-amplitude") == 0 )
    {
     arg_index++;
     amplitude = atof(argv[arg_index++]);
    }
   else if ( strcmp(argv[arg_index], "-nsteps") == 0 )
    {
     arg_index++;
     nsteps = atoi(argv[arg_index++]);
    }
   else if ( strcmp(argv[arg_index], "-period") == 0 )
    {
     arg_index++;
     period = atof(argv[arg_index++]);
    }
   else if ( strcmp(argv[arg_index], "-validation_run") == 0 )
    {
     arg_index++;
     validation_run=true;
    }
   else if ( strcmp(argv[arg_index], "-help") == 0 )
    {
     print_help = 1;
     break;
    }
   else
    {
     arg_index++;
    }
  }
 
 if (print_help)
  {
   if (MPI_Helpers::communicator_pt()->my_rank()==0)
    {
     oomph_info << "\n\nOption flags:\n";
     oomph_info << "-dir <n>                  Data saved to /RESLTn\n";
     oomph_info << "-outflow <ID>             Outflow type:\n"
                << "                          ID=0 flux control\n"
                << "                          ID=1 pressure driven with impedance"
                << " with LSC preconditioner\n"
                << "                          ID=2 pressure driven with impedance"
                << " with block triangluar and LSC preconditioner\n"
                << "                          ID=3 prescribed outflow\n"
               
                << "                          otherwise pressure driven flow\n"; 
     oomph_info << "-Re <value>               Sets Re=value\n";
     oomph_info << "-nsteps <value>           Number of timesteps\n";
     oomph_info << "-solver <ID>              Linear solver type:\n"
                << "                          ID=0 GMRES\n"
                << "                          otherwise SuperLU is used\n";
     oomph_info << "-p_solver <ID>            P matrix solve type:\n"
                << "                          ID=0 BoomerAMG\n"
                << "                          ID=1 ML SA-AMG\n"
                << "                          otherwise SuperLU is used\n";
     oomph_info << "-f_solver <ID>            F matrix solve type\n"
                << "                          ID=0 BoomerAMG\n"
                << "                          ID=1 ML SA-AMG\n"
                << "                          otherwise SuperLU is used\n";
     oomph_info << "-f_bamg_smoother <ID>     BoomerAMG smoother for F matrix\n"
                << "                          ID=0 damped Jacobi\n"
                << "                          ID=1 Gauss-Seidel\n";
     oomph_info << "-f_bamg_damping <value>   BoomerAMG Jacobi damping for F matrix\n";
     oomph_info << "-f_bamg_strength <value>  BoomerAMG strength paramater F matrix\n";
     oomph_info << "-f_ml_settings <ID>       Settings for ML on F matrix\n"
                << "                          ID=0 NSSA\n"
                << "                          otherwise SA\n";
     oomph_info << "-refine <n>               Refines mesh n times\n";
     oomph_info << "-period <value>           Set period to value\n";
     oomph_info << "-amplitude <value>        Set amplitude to value\n";
     oomph_info << "-validation_run           Generate validation data for 1D Womersley\n";
    }
//#ifdef OOMPH_HAS_MPI 
//   // finalize MPI
//  MPI_Helpers::finalize();
//#endif
   return (0);
  }

 // Set validation values
 if (validation_run)
  {
   if (outflow==0)
    {
     directory_for_data = "RESLT_flux_control";
    }
   else if (outflow==1)
    {
     directory_for_data = "RESLT_impedance_tube";
    }
   else if (outflow==2)
    {
     directory_for_data = "RESLT_impedance_tube_with_flux_control";
    }
   Global_Physical_Variables::Re = 100;
   Global_Physical_Variables::ReSt = 100;
   solver=0;
   refinement=0;
   amplitude=0.5;
   period=10.0;
  }
 
 // Let St=1
 Global_Physical_Variables::ReSt = Global_Physical_Variables::Re;

 // Output settings
 oomph_info << "Storing data in /"
            << directory_for_data << "\n";
 oomph_info << "Outflow type: "
            << outflow << "\n";
 oomph_info << "Re="
            <<  Global_Physical_Variables::Re
            << "\n";
 if (solver==0)
  {
   oomph_info << "Using GMRES solver\n";
   if (p_solver==0)
    {
#ifdef OOMPH_HAS_HYPRE
     oomph_info << "Using BoomerAMG on P matrix\n";
#else
     oomph_info 
       << "Warning: Hypre not available. Using SuperLU on P matrix\n";
#endif
    }
   else if (p_solver==1)
    {
#ifdef OOMPH_HAS_TRILINOS
     oomph_info << "Using ML SA-AMG on P matrix\n";
#else
     oomph_info 
       << "Warning: Trilinos not available. Using SuperLU on P matrix\n";
#endif
    }
   else
    {
     oomph_info << "Using SuperLU on P matrix\n";
    }
   if (f_solver==0)
    {
#ifdef OOMPH_HAS_HYPRE
     oomph_info << "Using BoomerAMG on F matrix\n"
                << "smoother: " << f_bamg_smoother
                << "\ndamping: " << f_bamg_damping
                << "\nstrength: " << f_bamg_strength 
                << "\n";
#else
     oomph_info 
       << "Warning: Hypre not available. Using SuperLU on F matrix\n";
#endif
    }
   else if (f_solver==1)
    {
#ifdef OOMPH_HAS_TRILINOS
     oomph_info << "Using ML SA-AMG on F matrix ";
     if (f_ml_settings==0)
      {
       oomph_info << "with NSSA settings\n";
      }
     else
      {
       oomph_info << "with SA settings\n";
      }
#else
     oomph_info 
       << "Warning: Trilinos not available. Using SuperLU on F matrix\n";
#endif
    }
   else
    {
     oomph_info << "Using SuperLU on F matrix\n";
    }
  }
 else
  {
   oomph_info << "Using SuperLU solver\n";
  }
 
 oomph_info << "Refining mesh " << refinement << " times\n";

 // Number of elements in the domain
 unsigned nup=20;
 unsigned ncollapsible=8;
 unsigned ndown=40;
 unsigned ny=16;

 // Length of the domain
 double lup=5.0;
 double lcollapsible=1.0;
 double ldown=10.0;
 double ly=1.0;
 double limpedance=0.0;
 
 // Modify values if using impedance tube
 if ((outflow==1) || (outflow==2))
  {
   limpedance = 5.0;
   ldown-=limpedance;
   ndown=20;
  }

 // Modify values for validation run
 if (validation_run)
  {
   limpedance = 8.0;
   ldown=2.0;
   nup=10;
   ncollapsible=4;
   ny=8;
   ndown=4;
  }
 
 oomph_info << "Amplitude=" << amplitude
            << "\nPeriod=" << period << "\n";
 
 // Solver stuff
 GMRES<CRDoubleMatrix>* iterative_solver_pt=0;
 NavierStokesSchurComplementPreconditioner* ns_preconditioner_pt=0;

 // Don't build the problem with Crouzeix Raviart Elements if using the 
 // SchurComplement preconditioner!
 CollapsibleChannelProblem<RefineableQTaylorHoodElement<2> > 
  problem(nup, ncollapsible, ndown, ny, 
          lup, lcollapsible, ldown, ly, limpedance,
          amplitude, period, outflow);

 //problem.switch_off_wall_oscillations();

 // Refine problem
 for (unsigned i=0;i<refinement;i++)
  {
   problem.refine_uniformly();
  }

 if (solver==0)
  {
   // Set up the solver and pass to the problem
   iterative_solver_pt = new GMRES<CRDoubleMatrix>;
   //  iterative_solver_pt->set_preconditioner_RHS();
   iterative_solver_pt->max_iter() = 200;
   iterative_solver_pt->tolerance() = 1.0e-8;
   problem.linear_solver_pt() = iterative_solver_pt;

   // Set up the preconditioner
   NavierStokesSchurComplementPreconditioner* ns_preconditioner_pt = 0;
   FSIPreconditioner* fsi_preconditioner_pt = 0;

   if (outflow==2)
    {
     // Construct the FSI preconditioner
     fsi_preconditioner_pt = new FSIPreconditioner(&problem);
     fsi_preconditioner_pt->enable_doc_time();

     // Get a pointer to the preconditioner
     ns_preconditioner_pt = 
      dynamic_cast<NavierStokesSchurComplementPreconditioner*>
      (fsi_preconditioner_pt->navier_stokes_preconditioner_pt());
     
     // Set solid mesh with "true" to tolerate multiple element types.
     fsi_preconditioner_pt->
      set_wall_mesh(problem.outflow_impedance_master_mesh_pt(),true);

     // Set fluid mesh
     fsi_preconditioner_pt->
      set_navier_stokes_mesh
      (problem.create_mesh_for_navier_stokes_preconditioner());

     // Pass the FSI preconditioner to the solver
     iterative_solver_pt->preconditioner_pt() = fsi_preconditioner_pt;
    }
   else
    {
     // Construct the preconditioner and pass in the problem pointer.
     ns_preconditioner_pt = new NavierStokesSchurComplementPreconditioner
      (&problem);

     // Setup the fluid mesh
     ns_preconditioner_pt->set_navier_stokes_mesh
      (problem.create_mesh_for_navier_stokes_preconditioner());
     
     // Pass the preconditioner to the solver
     iterative_solver_pt->preconditioner_pt() = ns_preconditioner_pt;
    }

   // Settings for the LCS preconditioner
   ns_preconditioner_pt->set_navier_stokes_mesh(
    problem.create_mesh_for_navier_stokes_preconditioner());
   ns_preconditioner_pt->enable_doc_time();
   
   // Set up the P sub block preconditioner
#ifdef OOMPH_HAS_HYPRE
   if (p_solver==0)
    {
#ifdef OOMPH_HAS_HYPRE
     Preconditioner* p_preconditioner_pt = new HyprePreconditioner;
     HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(p_preconditioner_pt);
     Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif
    }
#endif

#ifdef OOMPH_HAS_TRILINOS
   if (p_solver==1)
    {
#ifdef OOMPH_HAS_TRILINOS
     Preconditioner* p_preconditioner_pt = new TrilinosMLPreconditioner;
     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif
    }
#endif
    
   // Set up the F sub block precondioner
#ifdef OOMPH_HAS_HYPRE
   if (f_solver==0)
    {
#ifdef OOMPH_HAS_HYPRE
     Preconditioner* f_preconditioner_pt = new HyprePreconditioner;
     HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(f_preconditioner_pt);
     hypre_preconditioner_pt->use_BoomerAMG();
     hypre_preconditioner_pt->set_amg_iterations(1);
     hypre_preconditioner_pt->amg_using_simple_smoothing();
     hypre_preconditioner_pt->amg_simple_smoother()=f_bamg_smoother;
     hypre_preconditioner_pt->amg_damping()=f_bamg_damping;
     hypre_preconditioner_pt->amg_strength()=f_bamg_strength;
     ns_preconditioner_pt->set_f_preconditioner(f_preconditioner_pt);
#endif 
   }
#endif

#ifdef OOMPH_HAS_TRILINOS    
   if (f_solver==1)
    {
#ifdef OOMPH_HAS_TRILINOS
     Preconditioner* f_preconditioner_pt = new TrilinosMLPreconditioner;
     TrilinosMLPreconditioner* trilinos_preconditioner_pt = 
      static_cast<TrilinosMLPreconditioner*>(f_preconditioner_pt);
     if (f_ml_settings==1)
      {
       trilinos_preconditioner_pt->set_NSSA_default_values();
      }
     else
      {
       trilinos_preconditioner_pt->set_SA_default_values();
      }
     ns_preconditioner_pt->set_f_preconditioner(f_preconditioner_pt);
#endif
    }
#endif
  }
 
 // run the problem
 problem.unsteady_run(directory_for_data, nsteps, validation_run);

 // delete stuff - remembering P and F matrix preconditioners are 
 // deleted in the Schur complement preconditioner
 delete iterative_solver_pt;
 iterative_solver_pt=0;
 delete ns_preconditioner_pt;
 ns_preconditioner_pt=0;
 

//#ifdef OOMPH_HAS_MPI 
// // finalize MPI
// MPI_Helpers::finalize();
//#endif
 
} //end of driver code

  
