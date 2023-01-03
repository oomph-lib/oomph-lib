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
// Driver for a non-axisymmetric two fluid relaxing interface problem,
// where the mesh is deformed using a spine-based node-update strategy.

// Generic oomph-lib header
#include "generic.h"

// Linearised axisymmetric Navier-Stokes headers
#include "linearised_axisym_navier_stokes_elements.h"
#include "linearised_axisym_navier_stokes_elements.cc"
#include "multi_domain_linearised_axisym_navier_stokes_elements.h"

// Interface headers
#include "fluid_interface.h"
#include "linearised_axisymmetric_fluid_interface_elements.h"
#include "linearised_axisymmetric_fluid_interface_elements.cc"

// Perturbed spine headers
#include "perturbed_spines.h"
#include "perturbed_spines.cc"

// Bessel function headers
#include "oomph_crbond_bessel.h"

// The mesh
#include "meshes/two_layer_spine_mesh.h"
#include "two_layer_perturbed_spine_mesh.template.h"
#include "two_layer_perturbed_spine_mesh.template.cc"

// Catch floating point exceptions
#include <fenv.h> 

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;



//==start_of_namespace_for_physical_parameters===========================
/// Namespace for physical parameters
//=======================================================================
namespace GlobalPhysicalVariables
{

 /// Reynolds number
 double Re = 50.0;

 /// Strouhal number
 double St = 1.0;

 /// Womersley number (Reynolds x Strouhal, computed automatically)
 double ReSt;
 
 /// Product of Reynolds number and inverse of Froude number
 double ReInvFr = 50.0; // (Fr = 1)

 /// Ratio of viscosity in upper fluid to viscosity in lower
 /// fluid. Reynolds number etc. is based on viscosity in lower fluid.
 double Viscosity_Ratio = 0.1;

 /// Ratio of density in upper fluid to density in lower
 /// fluid. Reynolds number etc. is based on density in lower fluid.
 double Density_Ratio = 1.0;

 /// Capillary number
 double Ca = 1.0;

 /// Direction of gravity
 Vector<double> G(3);

 /// Vector of azimuthal mode numbers to investigate
 Vector<int> Vector_of_azimuthal_mode_numbers;

} // End of namespace



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//=======================================================================
/// Function-type-object to perform comparison of elements
//=======================================================================
class ElementCmp
{
 
public:
 
 /// Comparison. Are the values identical or not?
 bool operator()(GeneralisedElement* const &x,
                 GeneralisedElement* const &y) const
  {
   FiniteElement* cast_x = dynamic_cast<FiniteElement*>(x);
   FiniteElement* cast_y = dynamic_cast<FiniteElement*>(y);
   
   // Orders elements vertically, starting from the bottom of each spine
   if((cast_x==0) || (cast_y==0)) { return 0; }
   else
    {
     return ((cast_x->node_pt(0)->x(0) + 1.0e-8*cast_x->node_pt(0)->x(1))
             < (cast_y->node_pt(0)->x(0) + 1.0e-8*cast_y->node_pt(0)->x(1)));
    }
  }
};



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//=====start_of_base_state_problem_class=================================
/// Base state problem class for viscous two-layer rotating
/// cylinder problem
//=======================================================================
template<class BASE_ELEMENT,class TIMESTEPPER>
class BaseStateProblem : public Problem
{

public:

 /// Constructor: Pass the width of the domain in the r direction,
 /// and the heights of both bottom (fluid 1) and top (fluid 2) layers.
 /// Also pass the number of elements in both horizontal regions in the
 /// r direction and the number of elements in all three vertical regions
 /// in the z direction, along with the fractions which determine the
 /// spacings of those regions.
 BaseStateProblem(const unsigned &n_r, const unsigned &n_z1,
                  const unsigned &n_z2, const double &h1, const double &h2);

 /// Destructor (empty)
 ~BaseStateProblem() {}

 /// Set initial condition (incl previous timesteps)
 void set_initial_condition();

 /// Set the boundary conditions
 void set_boundary_conditions();

 /// Access function for the specific timestepper
 TIMESTEPPER* time_stepper_pt()
  {
   return dynamic_cast<TIMESTEPPER*>(Problem::time_stepper_pt());
  }

 /// Doc the solution
 void doc_solution(DocInfo* &doc_info_pt);

 /// Create interface elements at boundary between upper and lower layers
 void create_interface_elements();

 /// Access function for bulk mesh
 TwoLayerSpineMesh<BASE_ELEMENT>* bulk_mesh_pt() { return Bulk_mesh_pt; }

 /// Access function for surface mesh
 Mesh* surface_mesh_pt() { return Surface_mesh_pt; }

private:

 /// Spine heights/lengths are unknowns in the problem so their values get
 /// corrected during each Newton step. However, changing their value does
 /// not automatically change the nodal positions, so we need to update all
 /// of them here.
 void actions_before_newton_convergence_check()
  {
   Bulk_mesh_pt->node_update();
  }

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep() {}

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   dynamic_cast<BASE_ELEMENT*>(mesh_pt()->element_pt(e))
    ->fix_pressure(pdof,pvalue);
  }
 
 /// Pointer to the (specific) "bulk" mesh
 TwoLayerSpineMesh<BASE_ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

 /// Index at which the i-th velocity component is stored
 Vector<unsigned> U_nodal_index;

}; // End of base_state_problem class



//========start_of_base_state_constructor================================
/// Constructor for base state Tuckerman counter-rotating lids problem
//=======================================================================
template<class BASE_ELEMENT,class TIMESTEPPER>
BaseStateProblem<BASE_ELEMENT,TIMESTEPPER>::
BaseStateProblem(const unsigned &n_r, const unsigned &n_z1,
                 const unsigned &n_z2, const double &h1, const double &h2)
{
 // Always take one newton step even if the initial residuals are
 // below the required tolerance
 Problem::Always_take_one_newton_step = true;

 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new TIMESTEPPER);
 
 // Build and assign "bulk" mesh
 Bulk_mesh_pt = new TwoLayerSpineMesh<BASE_ELEMENT>
  (n_r,n_z1,n_z2,1.0,h1,h2,time_stepper_pt());

 // Create "surface mesh" that will contain only the interface elements.
 // The constructor just creates the mesh without giving it any elements,
 // nodes, etc.
 Surface_mesh_pt = new Mesh;

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();

 // -------------------------------------------------------------------
 // Get information from elements about the order of nodal data storage
 // -------------------------------------------------------------------

 // Get a pointer to the first element in the mesh -- note that we
 // are assuming that the indices will be the same in each element
 BASE_ELEMENT* el_pt = dynamic_cast<BASE_ELEMENT*>
  (Bulk_mesh_pt->element_pt(0));

 // Determine indices at which velocities are stored
 this->U_nodal_index.resize(3);
 for(unsigned i=0;i<3;i++)
  {
   U_nodal_index[i] = el_pt->u_index_axi_nst(i);
  }

 // --------------------------------------------
 // Set the boundary conditions for this problem
 // --------------------------------------------
 
 // Determine number of nodes in the mesh
 const unsigned n_node = Bulk_mesh_pt->nnode();

 // Pin all azimuthal velocities throughout the bulk of the domain
 for(unsigned n=0;n<n_node;n++)
  {
   for(unsigned i=0;i<3;i++)
    {
     Bulk_mesh_pt->node_pt(n)->pin(U_nodal_index[0]); 
    }
  }

 // ------------------------------------
 // Pin all velocity dofs in the problem
 // ------------------------------------
  
 // Loop over all nodes
 for(unsigned n=0;n<n_node;n++)
  {
   // Loop over the velocity components and pin the value
   for(unsigned i=0;i<3;i++) { Bulk_mesh_pt->node_pt(n)->pin(i); }
  }

 // ----------------------
 // Prescribe the pressure
 // ----------------------

 // Determine number of bulk elements in lower and upper fluid
 const unsigned n_lower = Bulk_mesh_pt->nlower();
 const unsigned n_upper = Bulk_mesh_pt->nupper();

 // Loop over elements in the lower layer
 for(unsigned e=0;e<n_lower;e++)
  {
   // Upcast from GeneralisedElement to the present element
   BASE_ELEMENT *el_pt = dynamic_cast<BASE_ELEMENT*>
    (Bulk_mesh_pt->lower_layer_element_pt(e));

   // In these 9-node elements, the 4-th node is always positioned at s1=0,s2=0
   // where s1, s2 are the local coordinates (and the "first" node is the 0-th)
   const double eulerian_z_pos_middle_node = el_pt->node_pt(4)->x(1);

   // Determine the value of the pressure at this node
   const double p_val_at_middle_node =
    GlobalPhysicalVariables::G[1]*
    GlobalPhysicalVariables::ReInvFr*eulerian_z_pos_middle_node;

   // Specify the pressure analytically
   el_pt->fix_pressure(0,p_val_at_middle_node);
   el_pt->fix_pressure(1,0.0);
   el_pt->fix_pressure(2,(GlobalPhysicalVariables::G[1]*
                          GlobalPhysicalVariables::ReInvFr/(n_z1+n_z2)));
  }

 // Loop over elements in the upper layer
 for(unsigned e=0;e<n_upper;e++)
  {
   // Upcast from GeneralisedElement to the present element
   BASE_ELEMENT *el_pt = dynamic_cast<BASE_ELEMENT*>
    (Bulk_mesh_pt->upper_layer_element_pt(e));

   // In these elements, the 4-th node is always positioned at s1=0,s2=0
   // where s1, s2 are the local coordinates (and the "first" node is the 0-th)
   const double eulerian_z_pos_middle_node = el_pt->node_pt(4)->x(1);

   // Determine the value of the pressure at this node
   const double p_val_at_middle_node =
    GlobalPhysicalVariables::G[1]*
    GlobalPhysicalVariables::Density_Ratio*
    GlobalPhysicalVariables::ReInvFr*eulerian_z_pos_middle_node;

   // Specify the pressure analytically
   el_pt->fix_pressure(0,p_val_at_middle_node
                       + (GlobalPhysicalVariables::G[1]*
                          GlobalPhysicalVariables::ReInvFr*
                          (1.0-GlobalPhysicalVariables::Density_Ratio)));
   el_pt->fix_pressure(1,0.0);
   el_pt->fix_pressure(2,(GlobalPhysicalVariables::G[1]*
                          GlobalPhysicalVariables::Density_Ratio*
                          GlobalPhysicalVariables::ReInvFr/(n_z1+n_z2)));
  }

 // ---------------------
 // Pin all spine heights
 // ---------------------

 // Determine the number of spines in the mesh
 const unsigned n_spine = Bulk_mesh_pt->nspine();

 // Loop over all spines
 for(unsigned i=0;i<n_spine;i++)
  {
   // Pin the spine height
   Bulk_mesh_pt->spine_pt(i)->spine_height_pt()->pin(0);

   // Set the value to the height of the lower fluid layer (h1)
   Bulk_mesh_pt->spine_pt(i)->spine_height_pt()->set_value(0,h1);
  }

 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------

 // Loop over bulk elements in lower fluid
 for(unsigned e=0;e<n_lower;e++)
  {
   // Upcast from GeneralisedElement to the present element
   BASE_ELEMENT* el_pt = dynamic_cast<BASE_ELEMENT*>
    (Bulk_mesh_pt->lower_layer_element_pt(e));
   
   // Set the Reynolds number
   el_pt->re_pt() = &GlobalPhysicalVariables::Re;
   
   // Set the Womersley number
   el_pt->re_st_pt() = &GlobalPhysicalVariables::ReSt;
   
   // Set the product of the Reynolds number and inverse Froude number
   el_pt->re_invfr_pt() = &GlobalPhysicalVariables::ReInvFr;
   
   // Set the direction of gravity
   el_pt->g_pt() = &GlobalPhysicalVariables::G;

  } // End of loop over bulk elements in lower fluid

 // Loop over bulk elements in upper fluid
 for(unsigned e=0;e<n_upper;e++)
  {
   // Upcast from GeneralisedElement to the present element
   BASE_ELEMENT* el_pt = dynamic_cast<BASE_ELEMENT*>
    (Bulk_mesh_pt->upper_layer_element_pt(e));

   // Set the Reynolds number
   el_pt->re_pt() = &GlobalPhysicalVariables::Re;
   
   // Set the Womersley number
   el_pt->re_st_pt() = &GlobalPhysicalVariables::ReSt;
   
   // Set the product of the Reynolds number and inverse Froude number
   el_pt->re_invfr_pt() = &GlobalPhysicalVariables::ReInvFr;
   
   // Set the direction of gravity
   el_pt->g_pt() = &GlobalPhysicalVariables::G;

   // Set the viscosity ratio
   el_pt->viscosity_ratio_pt() = &GlobalPhysicalVariables::Viscosity_Ratio;

   // Set the density ratio
   el_pt->density_ratio_pt() = &GlobalPhysicalVariables::Density_Ratio;

  } // End of loop over bulk elements in upper fluid

 // Set the pressure in the first element at 'node' 0 to 0.0
 fix_pressure(0,0,0.0);

 // Pin all the spine heights
 for(unsigned s=0;s<n_spine;s++)
  {
   Bulk_mesh_pt->spine_pt(s)->spine_height_pt()->pin(0);
  }

 // Set up equation numbering scheme
 std::cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

} // End of base state constructor



//======start_of_base_state_set_initial_condition========================
/// Set the initial conditions to be zero everywhere for base state
//=======================================================================
template<class BASE_ELEMENT,class TIMESTEPPER>
void BaseStateProblem<BASE_ELEMENT,TIMESTEPPER>::
set_initial_condition()
{
 // Determine number of nodes in mesh
 const unsigned n_node = mesh_pt()->nnode();
 
 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Loop over the three velocity components
   for(unsigned i=0;i<3;i++)
    {
     // Set velocity component i of node n to zero
     mesh_pt()->node_pt(n)->set_value(U_nodal_index[i],0.0);
    }
  }

 // Initialise the previous velocity values and positions for timestepping
 // corresponding to an impulsive start
 assign_initial_values_impulsive();

} // End of set_initial_condition for base state



//==start_of_base_state_doc_solution=====================================
/// Document the base state solution
//=======================================================================
template<class BASE_ELEMENT,class TIMESTEPPER>
void BaseStateProblem<BASE_ELEMENT,TIMESTEPPER>::
doc_solution(DocInfo* &doc_info_pt)
{
 ofstream some_file;
 char filename[256];
 
 // Set number of plot points (in each coordinate direction)
 const unsigned npts_bulk = 2;
 
 // Open solution output file
 sprintf(filename,"%s/base_soln%i.dat",
         doc_info_pt->directory().c_str(),
         doc_info_pt->number());
 some_file.open(filename);
 
 // Output solution to file
 Bulk_mesh_pt->output(some_file,npts_bulk);
 
 // Close solution output file
 some_file.close();
 
} // End of doc_solution for base state



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//==start_of_perturbed_state_problem_class===============================
/// Perturbed state problem class for Tuckerman counter-rotating lids
/// problem
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT,class TIMESTEPPER>
class PerturbedStateProblem : public Problem
{
 
public:
 
 /// Constructor: Pass the width of the domain in the r direction,
 /// and the heights of both bottom (fluid 1) and top (fluid 2) layers.
 /// Also pass the number of elements in both horizontal regions in the
 /// r direction and the number of elements in all three vertical regions
 /// in the z direction, along with the fractions which determine the
 /// spacings of those regions. Also pass a pointer to the base state mesh.
 PerturbedStateProblem(const unsigned &n_r,
                       const unsigned &n_z1, const unsigned &n_z2,
                       const double &h1, const double &h2,
                       SpineMesh* external_bulk_mesh_pt,
                       int azimuthal_mode_number);
 
 /// Destructor (empty)
 ~PerturbedStateProblem() {}

 /// Set initial conditions (velocities given random values in range [0,1])
 void set_initial_condition();

 /// Set up the (homogeneous) boundary conditions
 void set_boundary_conditions();
 
 /// Access function for the base state mesh
 TwoLayerSpineMesh<BASE_ELEMENT>* base_state_bulk_mesh_pt()
  {
   return dynamic_cast<TwoLayerSpineMesh<BASE_ELEMENT>*>
    (Base_state_bulk_mesh_pt);
  }

 /// Access function for the specific timestepper
 TIMESTEPPER* time_stepper_pt()
  {
   return dynamic_cast<TIMESTEPPER*>(Problem::time_stepper_pt());
  }

 /// Doc the solution
 void doc_solution(DocInfo* &doc_info_pt,const bool& output_soln=true);

 /// Create and initialise a trace file
 void create_trace_file(DocInfo* &doc_info_pt)
  {
   // Open trace file
   char filename[256];
   sprintf(filename,"%s/perturbed_trace_k%i.dat",
           doc_info_pt->directory().c_str(),
           Azimuthal_mode_number);
   Trace_file.open(filename);

   // Initialise
   Trace_file << "time, height (cosine), height (sine), growth_rate, d_kinetic_energy_dt, kinetic_energy" 
              << std::endl;
  }
 
 /// Clear and close trace file
 void close_trace_file() { Trace_file.clear(); Trace_file.close(); }
 
 /// Access function for trace file
 ofstream& trace_file() { return Trace_file; }

 /// Create interface elements at boundary between upper and lower layers
 void create_interface_elements();

 /// Access function for bulk mesh
 TwoLayerPerturbedSpineMesh<PERTURBED_ELEMENT>* bulk_mesh_pt()
  {
   return Bulk_mesh_pt;
  }

 /// Access function for surface mesh
 Mesh* surface_mesh_pt()
  {
   return Surface_mesh_pt;
  }

private:

 /// Update mesh to have nodal positions equal to that of the
 /// base state mesh
 void actions_before_newton_convergence_check()
  {
   Bulk_mesh_pt->node_update();
  }

 /// Update before solve (empty)
 void actions_before_newton_solve() {}
 
 /// Update after solve (empty)
 void actions_after_newton_solve() {}
 
 /// Actions before timestep (empty)
 void actions_before_implicit_timestep()
  {
   // Shift perturbed spine height history values
   // -------------------------------------------

   const unsigned n_perturbed_spine = Bulk_mesh_pt->nspine();
   for(unsigned i=0;i<n_perturbed_spine;i++)
    {
     Bulk_mesh_pt->perturbed_spine_pt(i)->height_pt()->time_stepper_pt()->
      shift_time_values(Bulk_mesh_pt->perturbed_spine_pt(i)->height_pt());
    }
  }
 
 /// Update the problem specs after solve (empty)
 void actions_after_implicit_timestep() {}

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned& e, const unsigned& pdof, 
                   const double& pvalue)
  {
   dynamic_cast<PERTURBED_ELEMENT*>(mesh_pt()->element_pt(e))
    ->fix_pressure(pdof,pvalue);
  }

 /// Pointer to the base state mesh
 SpineMesh* Base_state_bulk_mesh_pt;

 /// Pointer to the (specific) "bulk" mesh
 TwoLayerPerturbedSpineMesh<PERTURBED_ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

 /// Trace file
 ofstream Trace_file;

 /// Azimuthal mode number
 int Azimuthal_mode_number;

 /// Index at which the i-th velocity component is stored
 Vector<unsigned> U_nodal_index;

 /// Index at which the i-th component of the perturbation
 /// to the nodal coordinate is stored.
 Vector<unsigned> Xhat_nodal_index;

}; // End of perturbed_state_problem class



//==start_of_perturbed_state_constructor=================================
/// Constructor for perturbed state Tuckerman counter-rotating lids
/// problem
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT,class TIMESTEPPER>
PerturbedStateProblem
<BASE_ELEMENT,PERTURBED_ELEMENT,TIMESTEPPER>::
PerturbedStateProblem(const unsigned &n_r, const unsigned &n_z1,
                      const unsigned &n_z2, const double &h1, const double &h2,
                      SpineMesh* external_bulk_mesh_pt,
                      int azimuthal_mode_number)
 : Base_state_bulk_mesh_pt(external_bulk_mesh_pt),
   Azimuthal_mode_number(azimuthal_mode_number)
{
 // Always take one newton step even if the initial residuals are
 // below the required tolerance
 Problem::Always_take_one_newton_step = true;

 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new TIMESTEPPER);
 
 // Build and assign mesh
 Bulk_mesh_pt = new TwoLayerPerturbedSpineMesh<PERTURBED_ELEMENT>
  (n_r,n_z1,n_z2,1.0,h1,h2,Base_state_bulk_mesh_pt,time_stepper_pt());

 // Create "surface mesh" that will contain only the interface elements.
 // The constructor just creates the mesh without giving it any elements,
 // nodes, etc.
 Surface_mesh_pt = new Mesh;

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();

 // -------------------------------------------------------------------
 // Get information from elements about the order of nodal data storage
 // -------------------------------------------------------------------

 // Get a pointer to the first element in the mesh -- note that we
 // are assuming that the indices will be the same in each element
 PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>
  (Bulk_mesh_pt->element_pt(0));

 // Determine indices at which the perturbations to the nodal positions
 // are stored
 this->Xhat_nodal_index.resize(4);
 for(unsigned i=0;i<4;i++)
  {
   Xhat_nodal_index[i] = el_pt->xhat_index_lin_axi_nst(i);
  }

 // Determine indices at which velocities are stored
 this->U_nodal_index.resize(6);
 for(unsigned i=0;i<6;i++)
  {
   U_nodal_index[i] = el_pt->u_index_lin_axi_nst(i);
  }

 // Pass nodal indices of cosine (i=2) and sine (i=3) components of the
 // perturbation to the nodal z-position to the bulk mesh. This
 // information is used in the node_update function
 Bulk_mesh_pt->set_perturbation_to_nodal_positions_indices
  (Xhat_nodal_index[2],Xhat_nodal_index[3]);

 // --------------------------------------------
 // Set the boundary conditions for this problem
 // --------------------------------------------
 
 // All nodes are free by default -- just pin the ones that have
 // Dirichlet conditions here
 
 // Determine number of mesh boundaries
 const unsigned n_boundary = mesh_pt()->nboundary();

 // Loop over mesh boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine number of nodes on boundary b
   const unsigned n_node = mesh_pt()->nboundary_node(b);

   // Loop over nodes on boundary b
   for(unsigned n=0;n<n_node;n++)
    {
     // On top and bottom solid boundaries, pin all velocity components
     if(b==0 || b==2)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[0]);
       mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[1]);
       mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[2]);
       mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[3]);
       mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[4]);
       mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[5]);
      }
     // On outer wall, pin only radial and azimuthal velocities
     // (slippery outer wall)
     else if(b==1)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[0]);
       mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[1]);
       mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[4]);
       mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[5]);
      }
     // On symmetry boundary, boundary conditions differ depending on
     // the azimuthal mode number
     else if(b==3)
      {
       // If k=0, pin only radial and azimuthal velocity components
       if(Azimuthal_mode_number==0)
        {
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[0]);
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[1]);
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[4]);
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[5]);
        }
       // If k is ODD, pin axial velocity only
       else if(Azimuthal_mode_number%2)
        {
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[2]);
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[3]);
        }
       // If k is EVEN (and non-zero), pin all velocity components
       else
        {
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[0]);
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[1]);
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[2]);
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[3]);
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[4]);
         mesh_pt()->boundary_node_pt(b,n)->pin(U_nodal_index[5]);
        }
      }
    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries
 
 // Pin all perturbations to the nodal positions (since these are
 // dependent variables and therefore not dofs)
 const unsigned n_node = Bulk_mesh_pt->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   for(unsigned i=0;i<4;i++)
    {
     Bulk_mesh_pt->node_pt(n)->pin(Xhat_nodal_index[i]);
    }
  }

 // ---------------------------------------------------
 // If k==0, pin all dofs corresponding to "sine parts"
 // ---------------------------------------------------
 if(Azimuthal_mode_number==0)
  {
   // Pin all the sine components of the velocity values
   for(unsigned n=0;n<n_node;n++)
    {
     Bulk_mesh_pt->node_pt(n)->pin(U_nodal_index[1]);
     Bulk_mesh_pt->node_pt(n)->pin(U_nodal_index[3]);
     Bulk_mesh_pt->node_pt(n)->pin(U_nodal_index[5]);
    }

   // Pin all the sine components of the pressure dofs
   const unsigned n_element = Bulk_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Loop over the three pressure dofs (sine parts)
     for(unsigned i=0;i<3;i++)
      {
       // Pin the values and set them equal to zero
       dynamic_cast<PERTURBED_ELEMENT*>(mesh_pt()->element_pt(e))->
        fix_sine_component_of_pressure(i,0.0);
      }
    }

   // Pin all the sine components of the perturbed spine heights and set
   // them to zero
   const unsigned n_perturbed_spine = Bulk_mesh_pt->nspine();
   for(unsigned i=0;i<n_perturbed_spine;i++)
    {
     Bulk_mesh_pt->perturbed_spine_pt(i)->height_pt()->pin(1);
     Bulk_mesh_pt->perturbed_spine_pt(i)->height(1) = 0.0;
    }
  } // End of if k==0

 // -------------------------------------------------------------------
 // If k>0, pin interface height on the symmetry axis (both components)
 // -------------------------------------------------------------------
 if(Azimuthal_mode_number>0)
  {
   Bulk_mesh_pt->perturbed_spine_pt(0)->height_pt()->pin(0);
   Bulk_mesh_pt->perturbed_spine_pt(0)->height_pt()->pin(1);
   Bulk_mesh_pt->perturbed_spine_pt(0)->height(0) = 0.0;
   Bulk_mesh_pt->perturbed_spine_pt(0)->height(1) = 0.0;
  }

 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------

 // Determine number of bulk elements in lower and upper fluid
 const unsigned n_lower = Bulk_mesh_pt->nlower();
 const unsigned n_upper = Bulk_mesh_pt->nupper();

 // Loop over bulk elements in lower fluid
 for(unsigned e=0;e<n_lower;e++)
  {
   // Upcast from GeneralisedElement to the present element
   PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>
    (Bulk_mesh_pt->lower_layer_element_pt(e));
   
   // Set the Reynolds number
   el_pt->re_pt() = &GlobalPhysicalVariables::Re;
   
   // Set the Womersley number
   el_pt->re_st_pt() = &GlobalPhysicalVariables::ReSt;
   
   // Set the product of the Reynolds number and inverse Froude number
   el_pt->re_invfr_pt() = &GlobalPhysicalVariables::ReInvFr;
   
   // Set the direction of gravity
   el_pt->g_pt() = &GlobalPhysicalVariables::G;

   // Set the azimuthal mode number
   el_pt->azimuthal_mode_number_pt() = &Azimuthal_mode_number;

   // These terms in the Jacobian matrix are provided analytically, so
   // by filling them in from geometric data we would be "double counting"
   el_pt->enable_bypass_fill_in_jacobian_from_geometric_data();

  } // End of loop over bulk elements in lower fluid

 // Loop over bulk elements in upper fluid
 for(unsigned e=0;e<n_upper;e++)
  {
   // Upcast from GeneralisedElement to the present element
   PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>
    (Bulk_mesh_pt->upper_layer_element_pt(e));

   // Set the Reynolds number
   el_pt->re_pt() = &GlobalPhysicalVariables::Re;
   
   // Set the Womersley number
   el_pt->re_st_pt() = &GlobalPhysicalVariables::ReSt;
   
   // Set the product of the Reynolds number and inverse Froude number
   el_pt->re_invfr_pt() = &GlobalPhysicalVariables::ReInvFr;
   
   // Set the direction of gravity
   el_pt->g_pt() = &GlobalPhysicalVariables::G;

   // Set the viscosity ratio
   el_pt->viscosity_ratio_pt() = &GlobalPhysicalVariables::Viscosity_Ratio;

   // Set the density ratio
   el_pt->density_ratio_pt() = &GlobalPhysicalVariables::Density_Ratio;

   // Set the azimuthal mode number
   el_pt->azimuthal_mode_number_pt() = &Azimuthal_mode_number;

   // These terms in the Jacobian matrix are provided analytically, so
   // by filling them in from geometric data we would be "double counting"
   el_pt->enable_bypass_fill_in_jacobian_from_geometric_data();

  } // End of loop over bulk elements in upper fluid

 // Fix pressure ONLY when looking for axisymmetric modes
 if(Azimuthal_mode_number==0) { fix_pressure(0,0,0.0); }
 
 // ------------------------------------------------------------
 // Set up interaction between base and perturbed state problems
 // ------------------------------------------------------------

 // First interaction (base state velocities)
 Multi_domain_functions::setup_multi_domain_interaction<BASE_ELEMENT>
  (this,mesh_pt(),external_bulk_mesh_pt,0);

 // Second interaction (base state velocity derivatives w.r.t.
 // global spatial coordinates)
 Multi_domain_functions::setup_multi_domain_interaction<BASE_ELEMENT>
  (this,mesh_pt(),external_bulk_mesh_pt,1);

 // Third interaction (base state pressure)
 Multi_domain_functions::setup_multi_domain_interaction<BASE_ELEMENT>
  (this,mesh_pt(),external_bulk_mesh_pt,2);

 // Fourth interaction (base state velocity derivatives w.r.t. time)
 Multi_domain_functions::setup_multi_domain_interaction<BASE_ELEMENT>
  (this,mesh_pt(),external_bulk_mesh_pt,3);

 // Fourth interaction (base state velocity derivatives w.r.t. time)
 Multi_domain_functions::setup_multi_domain_interaction<BASE_ELEMENT>
  (this,mesh_pt(),external_bulk_mesh_pt,3);

 // Fifth interaction (base state velocity derivatives w.r.t.
 // local spatial coordinates)
 Multi_domain_functions::setup_multi_domain_interaction<BASE_ELEMENT>
  (this,mesh_pt(),external_bulk_mesh_pt,4);

 // Set up equation numbering scheme
 std::cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

} // End of perturbed state constructor



//======start_of_perturbed_state_set_initial_condition===================
/// Perturb the interface
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT,class TIMESTEPPER>
void PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT,TIMESTEPPER>::
set_initial_condition()
{
 // Determine number of nodes in mesh
 const unsigned n_node = mesh_pt()->nnode();
 
 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Loop over the six velocity components
   for(unsigned i=0;i<6;i++)
    {
     // Set velocity component i of node n to zero
     mesh_pt()->node_pt(n)->set_value(U_nodal_index[i],0.0);
    }
  }

 // Set value of epsilon
 const double epsilon = 0.01;

 // Determine number of perturbed_spines in mesh
 const unsigned n_perturbed_spine = Bulk_mesh_pt->nspine();
 
 // Determine number of history values used by timestepper (assume
 // same timestepper is used for all perturbed spine heights)
 const unsigned n_history_vals = Bulk_mesh_pt->perturbed_spine_pt(0)->
  height_pt()->time_stepper_pt()->ntstorage();

 // Initialise height
 double height = 0.0;
 
 // Initialise Bessel functions (only need the first!)
 double j0, j1, y0, y1, j0p, j1p, y0p, y1p;

 // Loop over perturbed_spines in mesh
 for(unsigned i=0;i<n_perturbed_spine;i++)
  {
   // Determine radial coordinate of perturbed_spine
   const double r = mesh_pt()->boundary_node_pt(0,i)->x(0);  
   
   // Compute Bessel functions
   const double k_bessel = 3.8317;
   CRBond_Bessel::bessjy01a(k_bessel*r,j0,j1,y0,y1,j0p,j1p,y0p,y1p);

   // Compute profile
   // In the axisymmetric case this has to be volume conserving
   if(Azimuthal_mode_number==0)
    {
     height = epsilon*j0;
    }
   // In the non-axisymmetric case this has to equal zero at r=0, but
   // the volume conservation is automatically handled by the theta-
   // dependence.
   else { height = epsilon*0.5*(1.0 - cos(2.0*MathematicalConstants::Pi*r)); }

   // Set perturbed_spine height (copy into all history values)
   for(unsigned t=0;t<n_history_vals;t++)
    {
     Bulk_mesh_pt->perturbed_spine_pt(i)->height(t,0) = height;
     Bulk_mesh_pt->perturbed_spine_pt(i)->height(t,1) = height;
    }
  } // End of loop over perturbed spines  

 // Update nodes in bulk mesh
 Bulk_mesh_pt->node_update();

 // Initialise the previous velocity values for timestepping
 // corresponding to an impulsive start
 assign_initial_values_impulsive();

 // Set the boundary conditions
 set_boundary_conditions();

} // End of set_initial_condition for perturbed state



//============start_of_perturbed_state_create_interface_elements=========
/// Create Interface Elements at the free surface between the two fluids
/// in the mesh pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT,class TIMESTEPPER>
void PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT,TIMESTEPPER>::
create_interface_elements()
{
 // Determine the number of horizontal elements
 const unsigned n_r = this->Bulk_mesh_pt->nx();

 // Determine number of bulk elements in the lower layer
 const unsigned n_z1 = this->Bulk_mesh_pt->ny1();

 // Loop over the horizontal elements
 for(unsigned e=0;e<n_r;e++)
  {
   // Construct a new 1D line element on the face on which the local
   // coordinate 1 is fixed at its max. value (1) -- Face 2
   FiniteElement* interface_element_element_pt =
    new PerturbedSpineLinearisedAxisymmetricFluidInterfaceElement
    <PERTURBED_ELEMENT>
    (this->Bulk_mesh_pt->finite_element_pt(n_r*(n_z1-1)+e),2);

   // Add the interface element to the surface mesh
   this->Surface_mesh_pt->add_element_pt(interface_element_element_pt);
  }

 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------

 // Determine number of 1D interface elements in mesh
 const unsigned n_interface_element = this->Surface_mesh_pt->nelement();

 // Loop over interface elements
 for(unsigned e=0;e<n_interface_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   PerturbedSpineLinearisedAxisymmetricFluidInterfaceElement
    <PERTURBED_ELEMENT>* el_pt = 
    dynamic_cast<PerturbedSpineLinearisedAxisymmetricFluidInterfaceElement
    <PERTURBED_ELEMENT>*>(this->Surface_mesh_pt->element_pt(e));

   // Set the Strouhal number
   el_pt->st_pt() = &GlobalPhysicalVariables::St;
   
   // Set the Capillary number
   el_pt->ca_pt() = &GlobalPhysicalVariables::Ca;
   
   // Set the azimuthal mode number
   el_pt->azimuthal_mode_number_pt() = &Azimuthal_mode_number;

   // These terms in the Jacobian matrix are provided analytically, so
   // by filling them in from geometric data we would be "double counting"
   el_pt->enable_bypass_fill_in_jacobian_from_geometric_data();
  }

} // End of create_interface_elements



//======start_of_perturbed_state_set_boundary_conditions=================
/// Set the (homogeneous) boundary conditions for the perturbed state
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT,class TIMESTEPPER>
void PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT,TIMESTEPPER>::
set_boundary_conditions()
{
 // Determine number of mesh boundaries
 const unsigned n_boundary = mesh_pt()->nboundary();
 
 // Loop over mesh boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine number of nodes on boundary b
   const unsigned n_node = mesh_pt()->nboundary_node(b);
   
   // Loop over nodes on boundary b
   for(unsigned n=0;n<n_node;n++)
    {
     // For the top and bottom lids set all components to zero
     if(b==0 || b==2)
      {
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[0],0.0);
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[1],0.0);
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[2],0.0);
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[3],0.0);
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[4],0.0);
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[5],0.0);
      }
     // For the outer (slippery) wall set the radial and azimuthal
     // components only to zero
     else if(b==1)
      {
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[0],0.0);
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[1],0.0);
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[4],0.0);
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[5],0.0);
      }
     // Symmetry boundary
     else if(b==3)
      {
       // If k=0, set only radial and azimuthal velocity components to zero
       if(Azimuthal_mode_number==0)
        {
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[0],0.0);
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[1],0.0);
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[4],0.0);
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[5],0.0);
        }
       // If k is ODD, set only axial velocity component to zero
       else if(Azimuthal_mode_number%2)
        {
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[2],0.0);
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[3],0.0);
        }
       // If k is EVEN (and non-zero), set all velocity components to zero
       else
        {
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[0],0.0);
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[1],0.0);
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[2],0.0);
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[3],0.0);
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[4],0.0);
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,U_nodal_index[5],0.0);
        }
      }
    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries
} // End of set_boundary_conditions for perturbed state



//==start_of_perturbed_state_doc_solution================================
/// Document the perturbed state solution
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT,class TIMESTEPPER>
void PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT,TIMESTEPPER>::
doc_solution(DocInfo* &doc_info_pt,const bool& output_soln)
{
 // Decide which spine to use for documenting heights in trace file
 // ---------------------------------------------------------------

 // Initialise "trace spine"
 unsigned trace_spine = 0;

 // Determine the number of spines in the mesh
 const unsigned n_spine = Bulk_mesh_pt->nspine();

 // Make decision based on azimuthal mode number
 if(Azimuthal_mode_number==0) { trace_spine = 0; }
 else { trace_spine = (n_spine-1)/2; }

 // Compute kinetic energy of perturbation
 // --------------------------------------

 // Initialise kinetic energy and its derivative w.r.t. time
 double kinetic_energy = 0.0;
 double d_kinetic_energy_dt = 0.0;

 // Loop over bulk elements in lower layer
 const unsigned n_lower = Bulk_mesh_pt->nlower();
 for(unsigned e=0;e<n_lower;e++)
  {
   // Upcast from GeneralisedElement to the present element
   PERTURBED_ELEMENT *el_pt = dynamic_cast<PERTURBED_ELEMENT*>
    (Bulk_mesh_pt->lower_layer_element_pt(e));

   // Initialise element's contributions to kinetic energy and its deriv
   double el_kinetic_energy = 0.0;
   double el_d_kinetic_energy_dt = 0.0;

   // Compute element's contributions to kinetic energy and its deriv
   el_pt->dkin_energy_dt(d_kinetic_energy_dt,el_kinetic_energy);

   // Add contribution to "global" kinetic energy and its deriv
   // (Don't forget to multiply through by the element's density ratio
   // in both cases)
   kinetic_energy += (el_pt->density_ratio())*el_kinetic_energy;
   d_kinetic_energy_dt += (el_pt->density_ratio())*el_d_kinetic_energy_dt;
  }

 // Loop over bulk elements in upper layer
 const unsigned n_upper = Bulk_mesh_pt->nupper();
 for(unsigned e=0;e<n_upper;e++)
  {
   // Upcast from GeneralisedElement to the present element
   PERTURBED_ELEMENT *el_pt = dynamic_cast<PERTURBED_ELEMENT*>
    (Bulk_mesh_pt->upper_layer_element_pt(e));

   // Initialise element's contributions to kinetic energy and its deriv
   double el_kinetic_energy = 0.0;
   double el_d_kinetic_energy_dt = 0.0;

   // Compute element's contributions to kinetic energy and its deriv
   el_pt->dkin_energy_dt(d_kinetic_energy_dt,el_kinetic_energy);

   // Add contribution to "global" kinetic energy and its deriv
   // (Don't forget to multiply through by the element's density ratio
   // in both cases)
   kinetic_energy += (el_pt->density_ratio())*el_kinetic_energy;
   d_kinetic_energy_dt += (el_pt->density_ratio())*el_d_kinetic_energy_dt;
  }

 // Compute growth rate by considering relative growth
 // (see A. Hazel and R. Hewitt's torus linear stab paper 2011, page 21)
 double growth_rate = 0.0;
 if(kinetic_energy!=0.0)
  {
   growth_rate = (0.5*d_kinetic_energy_dt)/kinetic_energy;
  }

 // Document in trace file
 Trace_file << time_pt()->time() << " " 
            << Bulk_mesh_pt->perturbed_spine_pt(trace_spine)->height(0) 
            << " "
            << Bulk_mesh_pt->perturbed_spine_pt(trace_spine)->height(1) 
            << " "
            << growth_rate << " "
            << d_kinetic_energy_dt << " "
            << kinetic_energy << std::endl;

 // If desired, output solution to file
 if(output_soln)
  {
   ofstream some_file;
   char filename[256];
   
   // Set number of plot points (in each coordinate direction)
   const unsigned npts_bulk = 5;
   const unsigned npts_surface = 5;
   
   // Open solution output file
   sprintf(filename,"%s/perturbed_k%i_soln%i.dat",
           doc_info_pt->directory().c_str(),
           Azimuthal_mode_number,
           doc_info_pt->number());
   some_file.open(filename);
   
   // Output solution to file
   Bulk_mesh_pt->output(some_file,npts_bulk);
   
   // Close solution output file
   some_file.close();

   // Open interface solution output file
   sprintf(filename,
           "%s/perturbation_to_interface_k%i_soln%i.dat",
           doc_info_pt->directory().c_str(),
           Azimuthal_mode_number,
           doc_info_pt->number());
   some_file.open(filename);
   
   // Output solution to file
   const unsigned n_interface_element = Surface_mesh_pt->nelement();
   for(unsigned e=0;e<n_interface_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     PerturbedSpineLinearisedAxisymmetricFluidInterfaceElement
      <PERTURBED_ELEMENT>* el_pt = 
      dynamic_cast<PerturbedSpineLinearisedAxisymmetricFluidInterfaceElement
      <PERTURBED_ELEMENT>*>(Surface_mesh_pt->element_pt(e));
     
     // Output solution to file
     el_pt->output_perturbation_to_interface(some_file,npts_surface);
    }
 
   // Close solution output file
   some_file.close();

   // Open interface solution output file
   sprintf(filename,
           "%s/combined_interface_position_k%i_soln%i.dat",
           doc_info_pt->directory().c_str(),
           Azimuthal_mode_number,
           doc_info_pt->number());
   some_file.open(filename);
   
   // Output solution to file
   for(unsigned e=0;e<n_interface_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     PerturbedSpineLinearisedAxisymmetricFluidInterfaceElement
      <PERTURBED_ELEMENT>* el_pt = 
      dynamic_cast<PerturbedSpineLinearisedAxisymmetricFluidInterfaceElement
      <PERTURBED_ELEMENT>*>(Surface_mesh_pt->element_pt(e));
     
     // Output solution to file
     el_pt->output_interface_position(some_file,npts_surface);
    }
 
   // Close solution output file
   some_file.close();

  }

} // End of doc_solution for perturbed state



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//==start_of_main========================================================
/// Driver for two-layer rotating cylinder problem, where we directly
/// simulate the instability by kicking off the base and perturbed state
/// problems together and watching to see whether the perturbation
/// grows or decays.
//=======================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Set duration of timestep
 const double dt = 0.05;

 // Set number of elements in radial (r) direction
 const unsigned base_n_r = 10;
 const unsigned perturbed_n_r = 10;
 
 // Set number of elements in axial (z) direction in lower fluid
 const unsigned base_n_z1 = 10;
 const unsigned perturbed_n_z1 = 10;

 // Set number of elements in axial (z) direction in upper fluid
 const unsigned base_n_z2 = 10;
 const unsigned perturbed_n_z2 = 10;

 // Height of lower fluid layer
 const double h1 = 1.0;
   
 // Height of upper fluid layer
 const double h2 = 1.0;

 // Compute the Womersley number
 GlobalPhysicalVariables::ReSt =
  GlobalPhysicalVariables::Re*GlobalPhysicalVariables::St;

 // Set direction of gravity (vertically downwards)
 GlobalPhysicalVariables::G[0] = 0.0;
 GlobalPhysicalVariables::G[1] = -1.0;
 GlobalPhysicalVariables::G[2] = 0.0;
 
 // Set azimuthal mode numbers to investigate
 for(int i=0;i<3;i++)
  {
   GlobalPhysicalVariables::Vector_of_azimuthal_mode_numbers.push_back(i);
  }

 // Set maximum time to run simulation for
 double t_max = 8.0;

 // If we are doing a validation run, set maximum time to be lower
 if(CommandLineArgs::Argc>1) { t_max = 0.1; }

 // Build base state problem
 BaseStateProblem<SpineElement<AxisymmetricQCrouzeixRaviartElement >,BDF<2> >
  base_problem(base_n_r,base_n_z1,base_n_z2,h1,h2);

 // Determine the size of the azimuthal mode number vector
 const unsigned n_perturbed_problems = 
  GlobalPhysicalVariables::Vector_of_azimuthal_mode_numbers.size();

 // Initialise vector of pointers to perturbation problems
 Vector<PerturbedStateProblem<SpineElement
  <AxisymmetricQCrouzeixRaviartElement >,
  PerturbedSpineElement
  <LinearisedAxisymmetricQCrouzeixRaviartMultiDomainElement >,
  BDF<2> >* > perturbed_problem_pt(n_perturbed_problems);

 // Build perturbed state problems
 for(unsigned i=0;i<n_perturbed_problems;i++)
  {
   perturbed_problem_pt[i] = new
    PerturbedStateProblem<SpineElement<AxisymmetricQCrouzeixRaviartElement >,
    PerturbedSpineElement
    <LinearisedAxisymmetricQCrouzeixRaviartMultiDomainElement >,
    BDF<2> > (perturbed_n_r,perturbed_n_z1,perturbed_n_z2,h1,h2,
              base_problem.bulk_mesh_pt(),
              GlobalPhysicalVariables::Vector_of_azimuthal_mode_numbers[i]);
  }
 
 // Create interface elements (this has to be done after
 // interaction is set up in perturbed state problem constructor)
 for(unsigned i=0;i<n_perturbed_problems;i++)
  {
   perturbed_problem_pt[i]->create_interface_elements();
  }
 
 // Rebuild global mesh from the sub-meshes
 base_problem.rebuild_global_mesh();
 for(unsigned i=0;i<n_perturbed_problems;i++)
  {
   perturbed_problem_pt[i]->rebuild_global_mesh();
  }

 // Assign equation numbers
 base_problem.assign_eqn_numbers();
 for(unsigned i=0;i<n_perturbed_problems;i++)
  {
   perturbed_problem_pt[i]->assign_eqn_numbers();
  }

 // Sort the elements in the base state problem (for frontal solver)
 std::sort(base_problem.mesh_pt()->element_pt().begin(),
           base_problem.mesh_pt()->element_pt().end(),
           ElementCmp());

 // Sort the elements in the perturbed state problem (for frontal solver)
 for(unsigned i=0;i<n_perturbed_problems;i++)
  {
   std::sort(perturbed_problem_pt[i]->mesh_pt()->element_pt().begin(),
             perturbed_problem_pt[i]->mesh_pt()->element_pt().end(),
             ElementCmp());
  }

 // Doc number of elements in both problems
 cout << "Number of elements:" << endl;
 cout << " - base, bulk         = "
      << base_problem.bulk_mesh_pt()->nelement() << endl;
 cout << " - base, surface      = "
      << base_problem.surface_mesh_pt()->nelement() << endl;
 for(unsigned i=0;i<n_perturbed_problems;i++)
  {
   cout << " - perturbed (" << i+1 << " of " 
        << n_perturbed_problems << "), bulk    = "
        << perturbed_problem_pt[i]->bulk_mesh_pt()->nelement() << endl;
   cout << " - perturbed (" << i+1 << " of " 
        << n_perturbed_problems << "), surface = "
        << perturbed_problem_pt[i]->surface_mesh_pt()->nelement() << endl;
  }

 // Initialise doc_info object
 DocInfo* doc_info_pt = new DocInfo();

 // Set output directory for solutions
 doc_info_pt->set_directory("RESLT");

 // Initialise counter for solutions
 doc_info_pt->number()=0;

 // Create and initialise trace file
 for(unsigned i=0;i<n_perturbed_problems;i++)
  {
   perturbed_problem_pt[i]->create_trace_file(doc_info_pt);
  }

 // Initialise timestep (also sets weights for all timesteppers)
 base_problem.initialise_dt(dt);
 for(unsigned i=0;i<n_perturbed_problems;i++)
  {
   perturbed_problem_pt[i]->initialise_dt(dt);
  }
 
 // Set initial conditions
 base_problem.set_initial_condition();
 for(unsigned i=0;i<n_perturbed_problems;i++)
  {
   perturbed_problem_pt[i]->set_initial_condition();
  }

 // Create info file
 ofstream info_file;
 char info_filename[100];
 sprintf(info_filename,"%s/info.dat",doc_info_pt->directory().c_str());
 info_file.open(info_filename);
 info_file.close();

 // Output base state solution
 base_problem.doc_solution(doc_info_pt);

 // Doc perturbed state initial conditions
 for(unsigned i=0;i<n_perturbed_problems;i++)
  {
   perturbed_problem_pt[i]->doc_solution(doc_info_pt,true);
  }

 // Increment counter for solutions
 doc_info_pt->number()++;

 // Determine number of timesteps
 const unsigned n_timestep = unsigned(t_max/dt);

 // Timestepping loop
 for(unsigned t=0;t<n_timestep;t++)
  {
   // Output timestep and global time to screen
   cout << "  Timestep " << t << " of " << n_timestep
        << "; Time is " << perturbed_problem_pt[0]->time_pt()->time()
        << endl;
   
   // Output timestep and global time to oomph output
   oomph_info << "================================"
              << "==============================\n"
              << "  Timestep " << t << " of " << n_timestep
              << "; Time is " << perturbed_problem_pt[0]->time_pt()->time()
              << "\n=============================="
              << "================================"
              << std::endl;

   // Solve perturbed problems
   for(unsigned i=0;i<n_perturbed_problems;i++)
    {
     oomph_info << "\nPerturbed problem (" << i+1 << " of " 
                << n_perturbed_problems 
                << "):\n-----------------------------\n" << endl;
     perturbed_problem_pt[i]->unsteady_newton_solve(dt);
    }
   
   // Doc all problems
   base_problem.doc_solution(doc_info_pt);
   for(unsigned i=0;i<n_perturbed_problems;i++)
    {
     perturbed_problem_pt[i]->doc_solution(doc_info_pt,true);
    }
   
   // Increment counter for solutions
   doc_info_pt->number()++;
   
  } // End of loop over timesteps
 
 // Clear and close trace files
 for(unsigned i=0;i<n_perturbed_problems;i++)
  {
   perturbed_problem_pt[i]->close_trace_file();
  }
 
} // End of main
