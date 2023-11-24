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
//Driver for a multi-physics problem that couples the free-surface
//Navier--Stokes equations to the surface transport equations for interfacial
//surfactant motion as well as the advection-diffusion equations in the bulk
//describing the transport of soluble surfactant and miceles. 

//At the moment the assumption that bulk surfactant is only present in
//the lower layer is handled by pinning all values in the upper layer
//and hijacking the values along the interface to prevent the upper
//layer making any contribution at all. If the problem is to be generalised
//to allow surfactant transport in both layers then an additional field
//will have to be added with a new variable that is pinned (and hijacked)
//in the lower field. Wherever possible I've tried to stick to the notation
//used in Kalogriou and Blyth

//N.B. This all requires more careful validation
//Simple checks that pass are that mass of surfactant is conserved, for
//a suitable small timestep; and that the steady solution matches the
//analytic expression.

//Oomph-lib headers, 
//We require the generic header
#include "generic.h"
//Our custom coupling of advection-diffusion and Navier--Stokes
#include "double_buoyant_navier_stokes_elements.h"
//The fluid interface elements
#include "fluid_interface.h"
//The surfactant transport equations
#include "soluble_surfactant_transport_equations.h"

// The mesh is our standard rectangular quadmesh
#include "meshes/two_layer_spine_mesh.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;

namespace oomph
{

namespace Control_Parameters
{
  bool Periodic_BCs = true;
}

  
//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{
  /// Geometry
  //----------------
  double L = 28.0;

  /// Fluid property Ratios
  //----------------------------
  
  //Density ratio:
  /// Ratio of density in upper fluid to density in lower
  /// fluid. Reynolds number etc. is based on density in lower fluid.
  double R = 1.0;

  //Vecosity ratio;
  /// Ratio of viscosity in upper fluid to viscosity in lower
  /// fluid. Reynolds number etc. is based on viscosity in lower fluid.
  double M = 0.5;

  //Dimensionless position of interface (relative to a total domain height of 1)
  double H0 = 0.2;

  //Dimensionless size of perturbation to the interface
  double Ha = H0/5.0;
  
  /// Hydrodynamic Parameters
  //----------------------------
  
  /// Reynolds number
  double Re = 0.0; 

  /// We do not need the Weber number,
  /// because we have non-dimensionalised pressure
  /// on the viscous scale and so multiplying the normal stress balanced by the
  /// Reynolds number gives a term in the Capillary number only (Ca Re = We).
  
  /// Capillary number (of which the results are independent
  /// for a pinned surface)
  double Ca = 0.001;
  
  /// In our non-dimensionalisation, we have a
  /// Reynolds number/ Froude number squared in
  /// (a material parameter that doesn't involve the velocity scale).
  /// 
  double ReInvFr = 0.0; // (Fr = 1)

  /// Surfactant Parameters
  //--------------------------
  
  /// Marangoni number
  double Ma = 10.0;

  /// Surface Elasticity number (Capillary number x Marangoni number)
  double Beta_s =Ca*Ma;

  /// Surface Peclet number
  double Pe_s = 10.0;
  
  /// Bulk Peclet number
  double Pe_b = 10.0;
  
  /// Micelle Pelect number
  double Pe_m = 10.0;

  /// Solubility Parameters
  //-------------------------
 
  /// Biot number
 double Biot = 0.1; 
  
  /// The ratio of adsorption-desorption times
  double K_b = 3.0; // 5.0; 

  // ratio of equilibrium concentrations
  double Beta_b = 1.0;

  // Reaction rate between bulk and micelle 
 double K_m = 1.0;

 /// Power of the concentration in bulk -> micelle flux expression
 double N = 10.0;
  
 /// The imposed pressure gradient
 double Delta_P = 1.0; 
  
 /// Timescales for transport equations (identically one from our
 /// non-dimensionalisation)
 Vector<double> Tau(2,1.0);
 
 /// Diffusivity  (should be 1/Pe_b, 1/Pe_m), which
 /// will be set in the main code
 Vector<double> D(2,1.0);
 
 /// Gravity vector, will be set in the main code
 Vector<double> Direction_of_gravity(2);
 
 /// This next set of functions is only used if we do NOT have
 /// periodic conditions

 
 /// Function that prescribes the hydrostatic pressure field at the outlet
 /// Let's fix things so that the pressure at the top of the channel is zero.
 void hydrostatic_pressure_outlet_upper(const double& time,
                                        const Vector<double> &x, 
                                        const Vector<double> &n, 
                                        Vector<double> &traction)
 {
  traction[0] = ReInvFr*R*Direction_of_gravity[1]*(1.0 - x[1]);
  traction[1] = 0.0;
 }

 /// Function that prescribes hydrostatic pressure field at the inlet
 void hydrostatic_pressure_inlet_upper(const double& time, const Vector<double> &x, 
				       const Vector<double> &n,
				       Vector<double> &traction)
 {
   traction[0] = Delta_P + -ReInvFr*R*Direction_of_gravity[1]*(1.0 - x[1]);
   traction[1] = 0.0;
 }

 
 /// Function that prescribes the hydrostatic pressure field at the outlet
 /// Must match pressure in lower fluid --- This may be tricky if the
 /// interface is not pinned
 /// (i.e. we'll need to read out the interfacial position
 /// on the boundary). For now assume it's at H0.
  void hydrostatic_pressure_outlet_lower(const double& time,
                                         const Vector<double> &x, 
                                  const Vector<double> &n, 
                                  Vector<double> &traction)
 {
  traction[0] = ReInvFr*Direction_of_gravity[1]*(R*(1.0 - H0) + H0 - x[1]);
  traction[1] = 0.0;
 }

 /// Function that prescribes hydrostatic pressure field at the inlet
 void hydrostatic_pressure_inlet_lower(const double& time, const Vector<double> &x, 
				       const Vector<double> &n,
				       Vector<double> &traction)
 {
  traction[0] = Delta_P +
   -ReInvFr*Direction_of_gravity[1]*(R*(1.0 - H0) + H0 - x[1]);
   traction[1] = 0.0;
 }

  //end of traction functions

 //Set specificied angle if required
  double Inlet_Angle = 2.0*atan(1.0);

  
 /// Direction of the wall normal vector (at the inlet)
 Vector<double> Wall_normal;

 /// Function that specifies the wall unit normal at the inlet
 void wall_unit_normal_inlet_fct(const Vector<double> &x, 
                                 Vector<double> &normal)
 {
  normal=Wall_normal;
 }

 /// Function that specified the wall unit normal at the outlet
 void wall_unit_normal_outlet_fct(const Vector<double> &x, 
                                 Vector<double> &normal)
 {
  //Set the normal
  normal = Wall_normal;
  //and flip the sign
  unsigned n_dim = normal.size();
  for(unsigned i=0;i<n_dim;++i) {normal[i] *= -1.0;}
 }

  
} // end_of_namespace


} //end of oomph namespace

/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////

//====== start_of_problem_class=======================================
/// 2D surfactant transport problem on rectangular domain, discretised 
/// with spine elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT,class INTERFACE_ELEMENT> 
class SurfactantProblem : public Problem
{

public:

 /// Constructor. The boolean indicates whether the free surface
 //should be pinned or not in the first instance
 SurfactantProblem(const bool &pin=true);

 /// Destructor. Empty
 ~SurfactantProblem() {}

 /// Release the free surface so that it can move
 void unpin_surface()
  {
   //Only bother if the surface is pinned
   if(Surface_pinned)
    {
     Surface_pinned = false;
     
     //Unpin the heights of all the spines in the middle
     unsigned n_spine = Bulk_mesh_pt->nspine();
     for(unsigned n=0;n<n_spine;n++)
      {
       Bulk_mesh_pt->spine_pt(n)->spine_height_pt()->unpin(0);
      }

     // Loop over the elements to re-enable ALE
     unsigned n_element = Bulk_mesh_pt->nelement();
     for(unsigned i=0;i<n_element;i++)
      {
       // Upcast from GeneralsedElement to the present element
       ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

       el_pt->enable_ALE();

       //For testing this is where you can pin the pressures
       //unsigned n_internal = el_pt->ninternal_data();
       //for(unsigned i=0;i<n_internal;++i)
       // {
       //  el_pt->internal_data_pt(i)->pin_all();
       // }
      }

     //Pin the velocities 
     /*unsigned n_node = Bulk_mesh_pt->nnode();
     for(unsigned n=0;n<n_node;++n)
      {
       Bulk_mesh_pt->node_pt(n)->pin(0);
       Bulk_mesh_pt->node_pt(n)->pin(1);
       }*/

     
     
     //Unpin the interfacial surfactant concentrations
     unsigned n_interface = Surface_mesh_pt->nelement();
     for(unsigned i=0;i<n_interface;i++)
      {
       FiniteElement *el_pt = Surface_mesh_pt->finite_element_pt(i);
       
       //Need to unpin the values of concentration and micelle
       unsigned n_el_node = el_pt->nnode();
       for(unsigned n=0;n<n_el_node;++n)
        {
         el_pt->node_pt(n)->unpin(4);
        }
      }

     
     //Now unpin the bulk concentrations in the lower region
     const unsigned n_lower = Bulk_mesh_pt->nlower();
      // Loop over bulk elements in lower fluid
      for(unsigned e=0;e<n_lower;e++)
	{
	  // Upcast from GeneralisedElement to the present element
	  ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
						  lower_layer_element_pt(e));
	  unsigned n_node = el_pt->nnode();
	  for(unsigned n=0;n<n_node;++n)
	    {
	      el_pt->node_pt(n)->unpin(2); //Bulk Concentration
              el_pt->node_pt(n)->unpin(3); //Micelle Concentration
	    }
	   }
      
     //Reassign the equation number
     std::cout << "Surface unpinned to give " 
               << assign_eqn_numbers() << " equation numbers\n";
    }
  }


 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Remember to update the nodes if the surface is not pinned
 void actions_before_newton_convergence_check()
  {if(!Surface_pinned) {Bulk_mesh_pt->node_update();}}

 /// Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep() 
  {set_boundary_conditions(time_pt()->time());}

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure


 /// UnFix pressure in element e at pressure dof pdof and set to pvalue
 void unfix_pressure(const unsigned &e, const unsigned &pdof)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    unfix_pressure(pdof);
  } // end_of_unfix_pressure


//Apply a prescribed deforamtion to the interface
void deform_interface(const double &epsilon,
		      const unsigned &n_periods)
{
  // Determine number of spines in mesh
  const unsigned n_spine = Bulk_mesh_pt->nspine();
  
  // Loop over spines in mesh
  for(unsigned i=0;i<n_spine;i++)
    {
      // Determine x coordinate of spine
      double x_value = Bulk_mesh_pt->boundary_node_pt(0,i)->x(0);
      
      // Set spine height
      Bulk_mesh_pt->spine_pt(i)->height() = 
	Global_Physical_Variables::H0 +
	epsilon*(cos(n_periods*
		     MathematicalConstants::Pi*(x_value - Global_Physical_Variables::L)/Global_Physical_Variables::L));
    }
  
  // Update nodes in bulk mesh
  Bulk_mesh_pt->node_update();
} // End of deform_free_surface

  
 /// Doc the solution.
 void doc_solution(std::ofstream &trace);

 /// Set the boundary conditions
 void set_boundary_conditions(const double &time);

 //Return the global error norm to be used in adaptive timestepping
 double global_temporal_error_norm();
  
 /// Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 TwoLayerSpineMesh<ELEMENT>* Bulk_mesh_pt;

 /// Storage to the mesh of Surface interface elements
 Mesh* Surface_mesh_pt;

 /// Storage to point elements if needed for non-periodic domains
 Mesh* Point_mesh_pt;

 /// Storage for any traction elements applied to the inlet
 Mesh* Inlet_traction_mesh_pt;

 /// Storage for any traction elements applied to the outlet
 Mesh* Outlet_traction_mesh_pt;

  void interface_min_max(double &min, double &max)
  {
    //Loop over the interface and add each elemental contribution
    const unsigned n_interface  = Surface_mesh_pt->nelement();
    if(n_interface > 0)
      {
	//Set the initial values to the first nodal position
	max = Surface_mesh_pt->finite_element_pt(0)->node_pt(0)->x(1);
	min = max;

	//Loop over all elements and find the max
	for(unsigned i=0;i<n_interface;i++)
	  {
	    // Upcast from GeneralsedElement to the present element
	    INTERFACE_ELEMENT *el_pt =
	      dynamic_cast<INTERFACE_ELEMENT*>(
					       Surface_mesh_pt->element_pt(i));

	    const unsigned n_node = el_pt->nnode();
	    if(n_node != 3)
	      {
		std::cout << "Can't do max and min for elements that are not quadratic\n";
		return;
	      }


	    //Read out the y values from the nodes
	    Vector<double> y(3);
	    for(unsigned n=0;n<3;++n)
	      {
		y[n] = el_pt->node_pt(n)->x(1);
	      }

	    double local_max = y[0];
	    double local_min = y[0];
	    //Find maximum and minimum nodes
	    for(unsigned n=1;n<3;++n)
	      {
		if(y[n] > local_max) {local_max = y[n];}
		if(y[n] < local_min) {local_min = y[n];}
	      }

	    //If we have a linear case then we are done
	    //Check that it's not a degenerate (linear) case
	    if(std::abs(y[0] - 2*y[1] + y[2]) > 1.0e-10)
	      {
		//Calculate extreme value of the local coordinate based on known
		//quadratic basis functions (This shoudl really be inside the element class)
		Vector<double> extreme_s(1,0.5*(y[0] - y[2])/(y[0] - 2.0*y[1] + y[2]));

		//Find the extreme height if the local coordinate is within the
		//rane of the element
		if(std::abs(extreme_s[0]) <= 1.0)
		  {
		    double extreme_h = el_pt->interpolated_x(extreme_s,1);
		    //Check whether the extreme value is greater than any of the nodes.
		    if(extreme_h > local_max) {local_max = extreme_h;}
		    if(extreme_h < local_min) {local_min = extreme_h;}
		  }
		  }

	    //Now check whether local max and min are global
	    if(local_max > max) {max = local_max;}
	    if(local_min < min) {min = local_min;}
	  }
      }
  }
    
    
 //Return the l2 norm of height difference between the interface
 //and its undeformed value
 double l2_norm_of_height(const double &h0)
  {
   double norm = 0.0;
   //Loop over the interface and add each elemental contribution
   const unsigned n_interface  = Surface_mesh_pt->nelement();
   for(unsigned i=0;i<n_interface;i++)
    {
     // Upcast from GeneralsedElement to the present element
     INTERFACE_ELEMENT *el_pt =
      dynamic_cast<INTERFACE_ELEMENT*>(
       Surface_mesh_pt->element_pt(i));
     
     norm += el_pt->l2_norm_of_height(h0);
    }
   return norm;
  }


  /// Return the total concentrations of the surfactant
  /// integrated over the bulk or surface accordingly
 void compute_integrated_concentrations(double &surface,
                                        double &bulk,
                                        double &micelle)
  {
   //Initialise to zero
   surface = 0.0;
   //Loop over the interface and add each elemental contribution
   const unsigned n_interface  = Surface_mesh_pt->nelement();
   for(unsigned i=0;i<n_interface;i++)
    {
     // Upcast from GeneralsedElement to the present element
     INTERFACE_ELEMENT *el_pt =
      dynamic_cast<INTERFACE_ELEMENT*>(
       Surface_mesh_pt->element_pt(i));
     
     surface += el_pt->integrated_C();
    }
   
   //Initialise to zero
   bulk = 0.0; micelle = 0.0;

   // Loop over bulk elements in lower fluid and add each elemental
   // contribution
   const unsigned n_lower = Bulk_mesh_pt->nlower();
   for(unsigned e=0;e<n_lower;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                             lower_layer_element_pt(e));
     double int_M=0.0, int_C=0.0;
     el_pt->integrated_C_and_M(int_C,int_M);
     bulk += int_C;
     micelle += int_M;
    }
  }


 
private:
 
 /// DocInfo object
 DocInfo Doc_info;

 /// Boolean to indicate whether the surface is pinned
 bool Surface_pinned;
  
}; // end of problem class

//===========start_of_constructor=========================================
/// Constructor for convection problem
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
SurfactantProblem<ELEMENT,INTERFACE_ELEMENT>::
SurfactantProblem(const bool &pin) : Surface_pinned(pin)
{
 //Allocate an (adaptive) timestepper
  add_time_stepper_pt(new BDF<2>(true));

 // Set output directory
 Doc_info.set_directory("RESLT");
 
 // # of elements in x-direction
 unsigned n_x=20;//100;

 // # of elements in y-direction
 unsigned n_y1=4;//8;

 unsigned n_y2=4;//8;

 //Domain length in x direction
 double lx = 2.0*Global_Physical_Variables::L;
 
 // Domain length in y-direction
 double h1=Global_Physical_Variables::H0;
 double h2=1.0 - h1;
 
 // Build a standard rectangular quadmesh
 Bulk_mesh_pt = 
   new TwoLayerSpineMesh<ELEMENT>(n_x,n_y1,n_y2,lx,h1,h2,
                                  Control_Parameters::Periodic_BCs,
                                  time_stepper_pt());


 //Create the surface mesh that will contain the interface elements
 //First create storage, but with no elements or nodes
 Surface_mesh_pt = new Mesh;
 //Make point elements at the end to compensate for
 //the unbalanced line tension terms 
 //if we DON'T have periodic boundaryc conditions
 if(!Control_Parameters::Periodic_BCs)
   {
     Point_mesh_pt = new Mesh;
   }
 
 const unsigned n_node_1d = Bulk_mesh_pt->finite_element_pt(0)->nnode_1d();
 //Loop over the horizontal elements
 for(unsigned i=0;i<n_x;i++)
  {
   //Construct a new 1D line element on the face on which the local
   //coordinate 1 is fixed at its max. value (1) --- This is face 2
   //This relies on knowing the ordering of the elements in the mesh
   //This is the upper face
   INTERFACE_ELEMENT *interface_element_pt =
    new INTERFACE_ELEMENT(
     Bulk_mesh_pt->finite_element_pt(n_x*(n_y1-1)+i),2);
    
   //Push it back onto the stack
   this->Surface_mesh_pt->add_element_pt(interface_element_pt); 
     

   //If it's not periodic, we need to make point elements to balance the surface tension
   if(!Control_Parameters::Periodic_BCs)
     {
       //Make a point element from left-hand side of the 
       //first surface element (note that this relies on knowledge of 
       //the element order within the mesh)
       if(i==0)
	 {
	   FluidInterfaceBoundingElement* point_element_pt =
	     interface_element_pt->make_bounding_element(-1);
	   //Add element to the point mesh
	   Point_mesh_pt->add_element_pt(point_element_pt);
	   //Set the capillary number
	   point_element_pt->ca_pt() =
	     &Global_Physical_Variables::Ca;
	   //Set the wall normal
	   point_element_pt->wall_unit_normal_fct_pt() = 
	     &Global_Physical_Variables::wall_unit_normal_inlet_fct;
	   //Set the contact angle
	   /*point_element_pt->set_contact_angle(
	     &Global_Physical_Variables::Inlet_Angle);*/
	 }
       
       //Make another point element from the right-hand side of the 
       //last surface element (note that this relies on knowledge of 
       //the element order within the mesh)
       if(i==n_x-1)
	 {
	   FluidInterfaceBoundingElement* point_element_pt =
	     interface_element_pt->make_bounding_element(1);
	   //Add element to the mesh
	   Point_mesh_pt->add_element_pt(point_element_pt);
	   //Set the capillary number
	   point_element_pt->ca_pt() =
	     &Global_Physical_Variables::Ca;
	   // Set the function that specifies the wall normal
	   point_element_pt->wall_unit_normal_fct_pt() = 
	     &Global_Physical_Variables::wall_unit_normal_outlet_fct;
	 }
     }
   
   //While we're here let's hijack the concentrations from the upper side
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->finite_element_pt(n_x*n_y1 + i));
   //Hijack all nodes of all elements
   unsigned n_max_node = n_node_1d;
   for(unsigned n=0;n<n_max_node;++n)
    {
     (void)el_pt->hijack_nodal_value(n,2,false);
     (void)el_pt->hijack_nodal_value(n,3,false);
    }
  }


 //If it's not periodic we also need to add traction
 //boundary conditions to impose the pressure drop
 if(!Control_Parameters::Periodic_BCs)
   {
     //Make the inlet and outlet traction meshes
     Inlet_traction_mesh_pt = new Mesh;
     Outlet_traction_mesh_pt = new Mesh;
     
     //Inlet boundary conditions (boundary 4 and 5)
     //Upper inlet
     {
       unsigned b=4;
       //Find the number of elements adjacent to mesh boundary
       unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(b);
       //Loop over these elements and create the traction elements
       for(unsigned e=0;e<n_boundary_element;e++)
	 {
	   SpineElement<NavierStokesTractionElement<ELEMENT> >*surface_element_pt =
	     new SpineElement<NavierStokesTractionElement<ELEMENT> >
	     (Bulk_mesh_pt->boundary_element_pt(b,e),
	      Bulk_mesh_pt->face_index_at_boundary(b,e));
	   //Add the elements to the mesh
	   Inlet_traction_mesh_pt->add_element_pt(surface_element_pt);
	   //Set the traction function
	   surface_element_pt->traction_fct_pt() = 
	     &Global_Physical_Variables::hydrostatic_pressure_inlet_upper;
	 }
     }
     
     //Lower inlet
     {
       unsigned b=5;
       //Find the number of elements adjacent to mesh boundary
       unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(b);
       //Loop over these elements and create the traction elements
       for(unsigned e=0;e<n_boundary_element;e++)
	 {
	   SpineElement<NavierStokesTractionElement<ELEMENT> >*surface_element_pt =
	     new SpineElement<NavierStokesTractionElement<ELEMENT> >
	     (Bulk_mesh_pt->boundary_element_pt(b,e),
	      Bulk_mesh_pt->face_index_at_boundary(b,e));
	   //Add the elements to the mesh
	   Inlet_traction_mesh_pt->add_element_pt(surface_element_pt);
	   //Set the traction function
	   surface_element_pt->traction_fct_pt() = 
	     &Global_Physical_Variables::hydrostatic_pressure_inlet_lower;
	 }
     }
     
     //Outlet boundary conditions (boundary 1 and 2)
     //Lower outlet
     {
       unsigned b=1;
       //Find the number of elements adjacent to mesh boundary
       unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(b);
       //Loop over these elements and create the traction elements
       for(unsigned e=0;e<n_boundary_element;e++)
	 {
	   SpineElement<NavierStokesTractionElement<ELEMENT> >*surface_element_pt =
	     new SpineElement<NavierStokesTractionElement<ELEMENT> >
	     (Bulk_mesh_pt->boundary_element_pt(b,e),
	      Bulk_mesh_pt->face_index_at_boundary(b,e));
	   //Add the elements to the mesh
	   Outlet_traction_mesh_pt->add_element_pt(surface_element_pt);
	   //Set the traction function
	   surface_element_pt->traction_fct_pt() = 
	     &Global_Physical_Variables::hydrostatic_pressure_outlet_lower;
	 }
     }
     
     //Upper outlet
     {
       unsigned b=2;
       //Find the number of elements adjacent to mesh boundary
       unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(b);
       //Loop over these elements and create the traction elements
       for(unsigned e=0;e<n_boundary_element;e++)
	 {
	   SpineElement<NavierStokesTractionElement<ELEMENT> >*surface_element_pt =
	     new SpineElement<NavierStokesTractionElement<ELEMENT> >
	     (Bulk_mesh_pt->boundary_element_pt(b,e),
	      Bulk_mesh_pt->face_index_at_boundary(b,e));
	   //Add the elements to the mesh
	   Outlet_traction_mesh_pt->add_element_pt(surface_element_pt);
	   //Set the traction function
	   surface_element_pt->traction_fct_pt() = 
	     &Global_Physical_Variables::hydrostatic_pressure_outlet_upper;
	 }
     }
   }


 // Add the two sub-meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);
 if(!Control_Parameters::Periodic_BCs)
   {
     add_sub_mesh(Point_mesh_pt);
     add_sub_mesh(Inlet_traction_mesh_pt);
     add_sub_mesh(Outlet_traction_mesh_pt);
   }
 
 // Combine all sub-meshes into a single mesh
 build_global_mesh();

 
 //Pin the heights of all the spines if the surface is pinned
 if(Surface_pinned)
  {
   unsigned n_spine = Bulk_mesh_pt->nspine();
   for(unsigned n=0;n<n_spine;n++)
    {
     Bulk_mesh_pt->spine_pt(n)->spine_height_pt()->pin(0);
    }
  }

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here
 
 //Loop over the boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //If we are on the side-walls, the concentrations
   //satisfy natural boundary conditions, so we only pin the
   //v-velocity for now if not periodic
    if(!Control_Parameters::Periodic_BCs)
      {
	if((ibound==1) || (ibound==2) || (ibound==4) || (ibound==5))
	  {
	    //Loop over the number of nodes on the boundary
	    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
	    for(unsigned inod=0;inod<num_nod;inod++)
	      {
		Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
		/*if((ibound==4) || (ibound==5))
		  {
		  Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
		  }*/
	      }
	  }
      }

   //If we on the top or bottom wall, velocity is pinned
   if((ibound==3) || (ibound==0))
    {
      //Loop over the number of nodes on the boundary
      unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
      for(unsigned inod=0;inod<num_nod;inod++)
	{
	  Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
	  Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
	}
    }
  }

 
 //Pin the zero-th pressure dof in element 0 and set its value to
 //zero:
 fix_pressure(0,0,0.0);

 // Complete the build of all elements so they are fully functional 

 // Determine number of bulk elements in lower/upper fluids
 const unsigned n_lower = Bulk_mesh_pt->nlower();
 const unsigned n_upper = Bulk_mesh_pt->nupper();

 // Loop over bulk elements in lower fluid
 for(unsigned e=0;e<n_lower;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                           lower_layer_element_pt(e));

   // Set the diffusivities number
   el_pt->diff_pt() = &Global_Physical_Variables::D;

   // Set the timescales
   el_pt->tau_pt() =&Global_Physical_Variables::Tau;

   // Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::Re;

   // Set the product of the Reynolds number and the inverse of the
   // Froude number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

   // Set the km parameter
   el_pt->km_pt() = &Global_Physical_Variables::K_m;

   // Set the N parameter
   el_pt->n_pt() = &Global_Physical_Variables::N;

   
   // Set the direction of gravity
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;

  } // End of loop over bulk elements in lower fluid

 // Loop over bulk elements in upper fluid 
 for(unsigned e=0;e<n_upper;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                           upper_layer_element_pt(e));

   // Set the diffusivities number
   el_pt->diff_pt() = &Global_Physical_Variables::D;

   // Set the timescales
   el_pt->tau_pt() =&Global_Physical_Variables::Tau;
   
   // Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::Re;

   // Set the product of the Reynolds number and the inverse of the
   // Froude number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

   // Set the viscosity ratio
   el_pt->viscosity_ratio_pt() = &Global_Physical_Variables::M;

   // Set the density ratio
   el_pt->density_ratio_pt() = &Global_Physical_Variables::R;

   // Set the km parameter
   el_pt->km_pt() = &Global_Physical_Variables::K_m;

   // Set the N parameter
   el_pt->n_pt() = &Global_Physical_Variables::N;
   
   // Set the direction of gravity
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;
   
   //Need to pin the values of concentration and micelle up here
   unsigned n_el_node = el_pt->nnode();
   for(unsigned n=0;n<n_el_node;++n)
     {
       //el_pt->node_pt(n)->set_value(2,0.0);
       //el_pt->node_pt(n)->set_value(3,0.0);
       el_pt->node_pt(n)->pin(2);
       el_pt->node_pt(n)->pin(3);
     }
   
  } // End of loop over bulk elements in upper fluid


  // Loop over the interface elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor.
 unsigned n_interface  = Surface_mesh_pt->nelement();
 for(unsigned i=0;i<n_interface;i++)
  {
   // Upcast from GeneralsedElement to the present element
   INTERFACE_ELEMENT *el_pt = dynamic_cast<INTERFACE_ELEMENT*>(
    Surface_mesh_pt->element_pt(i));
   
   // Set the Biot number
   el_pt->bi_pt() = &Global_Physical_Variables::Biot;

   // Set the Marangoni number
   el_pt->ma_pt() =&Global_Physical_Variables::Ma;

   // Set the Ca number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;

   // Set the surface elasticity number
   el_pt->beta_pt() = &Global_Physical_Variables::Beta_s;

   // Set the surface peclect number
   el_pt->peclet_s_pt() = &Global_Physical_Variables::Pe_s;

   // Set the surface peclect number multiplied by strouhal number
   el_pt->peclet_strouhal_s_pt() = &Global_Physical_Variables::Pe_s;

   // Set the reaction ratio
   el_pt->k_pt() = &Global_Physical_Variables::K_b;

   //el_pt->pe_b_pt() = &Global_Physical_Variables::Pe_b;

   el_pt->beta_b_pt() = &Global_Physical_Variables::Beta_b;
   
   //Need to unpin the values of concentration and micelle
   unsigned n_el_node = el_pt->nnode();
   for(unsigned n=0;n<n_el_node;++n)
     {
       el_pt->node_pt(n)->unpin(2); //unpin
       el_pt->node_pt(n)->unpin(3); //unpin
       //el_pt->node_pt(n)->pin(4);
     }

  }

 //Pin all bulk concentrations ... to save
 //time and to avoid having to conserve global mass
 unsigned n_node = Bulk_mesh_pt->nnode();
 for(unsigned n=0;n<n_node;++n)
   {
     Bulk_mesh_pt->node_pt(n)->pin(2);
     Bulk_mesh_pt->node_pt(n)->pin(3);
     }

 
 //Pin one surface concentration at the surface, if only
 //solving for insoluble surfactant
 //Surface_mesh_pt->finite_element_pt(0)->node_pt(1)->pin(4);
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor



//===========start_of_set_boundary_conditions================
/// Set the boundary conditions as a function of continuous 
/// time
//===========================================================
template<class ELEMENT,class INTERFACE_ELEMENT>
void SurfactantProblem<ELEMENT,INTERFACE_ELEMENT>::set_boundary_conditions(
 const double &time)
{
 //Set initial temperature profile
 if(time <= 0.0)
  {
    const double Gamma_init = 0.75;
    const double C_init = Gamma_init/(Global_Physical_Variables::K_b*(1.0-Gamma_init));
    const double M_init = pow(C_init,Global_Physical_Variables::N);
    unsigned n_lower = Bulk_mesh_pt->nlower();//nnode();
    for(unsigned e=0;e<n_lower;e++)
      {
	FiniteElement* el_pt =
	  Bulk_mesh_pt->lower_layer_element_pt(e);
	unsigned n_node = el_pt->nnode();
	for(unsigned n=0;n<n_node;n++)
	  {
	    Node* nod_pt = el_pt->node_pt(n);
	    //And in uniformly distributed surfactant
	    //Be careful about upper and lower layers
	    //If these are not set 
	    nod_pt->set_value(2,C_init); 
	    nod_pt->set_value(3,M_init);
	    
	    //Set the velocity
	    /*double y = nod_pt->x(1);
	      nod_pt->set_value(0,0.5*y);
	      nod_pt->set_value(1,0.0);
	      std::cout << nod_pt->x(0) << " "
	      << nod_pt->x(1) << " " <<
		 nod_pt->eqn_number(0) << "\n";*/
	  }
      }
  
   //Set the initial surface concentration to be one
   unsigned n_surface_element = Surface_mesh_pt->nelement();
   for(unsigned e=0;e<n_surface_element;++e)
     {
       unsigned n_el_node = Surface_mesh_pt->finite_element_pt(e)->nnode();
       for(unsigned n=0;n<n_el_node;++n)
	 {
	   Surface_mesh_pt->finite_element_pt(e)->node_pt(n)->set_value(4,Gamma_init);
	 }
     }
  }

 // Loop over the boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

     if(!Control_Parameters::Periodic_BCs)
       {
	 //If we are on the side walls we only set the v-velocity.
	 if((ibound==1) || (ibound==2) || (ibound==4) || (ibound==5))
	   {
	     nod_pt->set_value(1,0.0); 
	   }
       }

     //If we are on the top boundary, do not set the velocities
     //(yet)
     if(ibound==3)
       {
	 nod_pt->set_value(0,1.0); nod_pt->set_value(1,0.0);
       }

     //If we are on the bottom boundary
     if(ibound==0)
       {
	 nod_pt->set_value(0,0.0); nod_pt->set_value(1,0.0);
       }
    }
  }
} // end_of_set_boundary_conditions

//===============start_doc_solution=======================================
/// Doc the solution
//========================================================================
template<class ELEMENT,class INTERFACE_ELEMENT>
void SurfactantProblem<ELEMENT,INTERFACE_ELEMENT>::doc_solution(
 ofstream &trace)
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   Bulk_mesh_pt->finite_element_pt(e)->output(some_file,npts);
  }
 some_file.close();

 //Output the interface
 sprintf(filename,"%s/int%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);

 unsigned n_interface  = Surface_mesh_pt->nelement();
 for(unsigned i=0;i<n_interface;i++)
  {
   Surface_mesh_pt->finite_element_pt(i)->output(some_file,npts);
  }
 some_file.close();

 unsigned n_surface_element = Surface_mesh_pt->nelement();
 
 Node* monitor_node_pt = Surface_mesh_pt->finite_element_pt(n_surface_element/2)->node_pt(1);

 //Let's get the mases
 double surface=0.0, bulk=0.0, micelle=0.0;
 this->compute_integrated_concentrations(surface,bulk,micelle);

 //Compute the max and min
 double max = 0.0; double min =0.0;
 this->interface_min_max(min,max);
 
 trace << time_pt()->time() << " " 
       << monitor_node_pt->x(1) << " " 
       << monitor_node_pt->value(2) << " "
       << monitor_node_pt->value(3) << " "
       << monitor_node_pt->value(4) << " "
       << std::sqrt(this->l2_norm_of_height(Global_Physical_Variables::H0)/(2.0*Global_Physical_Variables::L)) << " "
       << surface << " " << bulk << " " << micelle << " "
       << surface + (bulk + micelle)/(Global_Physical_Variables::Beta_b) << " "
       << min << " " << max << " "
       << std::endl;


 Doc_info.number()++;
} // end of doc



//Specify the global error norm
template<class ELEMENT,class INTERFACE_ELEMENT>
double SurfactantProblem<ELEMENT,INTERFACE_ELEMENT>::global_temporal_error_norm()
{
 //Temp
 double global_error = 0.0;
   
 //Find out how many nodes there are in the problem
 unsigned long Nnode = Bulk_mesh_pt->nnode();

 //Loop over the nodes and calculate the errors in the positions
 for(unsigned long i=0;i<Nnode;i++)
  {
   //Find number of dimensions of the node
   unsigned Ndim = Bulk_mesh_pt->node_pt(i)->ndim();
   //Set the position error to zero
   double node_position_error = 0.0;
   //Loop over the dimensions (only need j=1)
   //for(unsigned j=0;j<Ndim;j++)
   unsigned j=1;
   {
     //Get position error
     double error = 
      Bulk_mesh_pt->node_pt(i)->position_time_stepper_pt()->
      temporal_error_in_position(Bulk_mesh_pt->node_pt(i),j);

     //Add the square of the individual error to the position error
     node_position_error += error*error;
    }
    
   //Divide the position error by the number of dimensions
   node_position_error /= Ndim;
   //Now add to the global error
   global_error += node_position_error;
  }
   
   //Now the global error must be divided by the number of nodes
 global_error /= Nnode;

 //Return the square root of the errr
 return sqrt(global_error);
}
 

//=======start_of_main================================================
/// Driver code for 2D Boussinesq convection problem
//====================================================================
int main(int argc, char **argv)
{
 ofstream trace("RESLT/trace.dat");

 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;

 //Set the diffusivities (inverse peclect numbers)
 Global_Physical_Variables::D[0] = 1.0/Global_Physical_Variables::Pe_b;
 Global_Physical_Variables::D[1] = 1.0/Global_Physical_Variables::Pe_m;

 Global_Physical_Variables::Wall_normal.resize(2);
 Global_Physical_Variables::Wall_normal[0] = -1.0;
 Global_Physical_Variables::Wall_normal[1] = 0.0;
 
 //Construct the problem
SurfactantProblem<SpineElement<Hijacked<DoubleBuoyantQCrouzeixRaviartElement<2> > >,
		   SpineLineSolubleSurfactantTransportInterfaceElement<DoubleBuoyantQCrouzeixRaviartElement<2> > > 
                   problem;

 
 // Apply the boundary condition at time zero
 problem.set_boundary_conditions(0.0);
 
 //Perform a single steady Newton solve
 problem.steady_newton_solve();

 //Document the solution
 problem.doc_solution(trace);

 //Now release the interface for real fun
 problem.unpin_surface();

 //Set the timestep
 double dt = 0.1;

 //Initialise the value of the timestep and set initial values 
 //of previous time levels assuming an impulsive start.
 problem.deform_interface(Global_Physical_Variables::Ha,1);
 problem.doc_solution(trace);
 problem.assign_initial_values_impulsive(dt);

 //Set the number of timesteps to our default value
 unsigned n_steps = 1000;

 //If we have a command line argument, perform fewer steps 
 //(used for self-test runs)
 if(argc > 1) {n_steps = 5;}

 //Perform n_steps timesteps
 for(unsigned i=0;i<n_steps;++i)
  {
   // problem.unsteady_newton_solve(dt);
   double dt_next = problem.adaptive_unsteady_newton_solve(dt,1.0e-5);
   dt = dt_next;
   //Limit timestep
   //if(dt > 1.0) {dt = 1.0;}
   problem.doc_solution(trace);
  }

} // end of main









