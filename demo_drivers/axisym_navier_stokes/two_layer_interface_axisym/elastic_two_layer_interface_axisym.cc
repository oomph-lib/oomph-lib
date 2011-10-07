//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
// Driver for an adaptive axisymmetric two fluid interface problem,
// where the mesh is deformed using a pseudo-solid node-update strategy
 
// Generic oomph-lib header
#include "generic.h"

// Axisymmetric Navier-Stokes headers
#include "axisym_navier_stokes.h"
#include "navier_stokes.h"

// Interface headers
#include "fluid_interface.h"

// Constitutive law headers
#include "constitutive.h"

// Solid headers
#include "solid.h"

// Bessel function headers
#include "oomph_crbond_bessel.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;


//==start_of_namespace====================================================
/// Namespace for physical parameters
//========================================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re = 5.0;

 /// Strouhal number
 double St = 1.0;

 /// Womersley number (Reynolds x Strouhal, computed automatically)
 double ReSt;
 
 /// Product of Reynolds number and inverse of Froude number
 double ReInvFr = 5.0; // (Fr = 1)

 /// \short Ratio of viscosity in upper fluid to viscosity in lower
 /// fluid. Reynolds number etc. is based on viscosity in lower fluid.
 double Viscosity_Ratio = 0.1;

 /// \short Ratio of density in upper fluid to density in lower
 /// fluid. Reynolds number etc. is based on density in lower fluid.
 double Density_Ratio = 0.5;

 /// Capillary number
 double Ca = 0.01;

 /// Direction of gravity
 Vector<double> G(3);

 /// Pseudo-solid Poisson ratio
 double Nu = 0.1;

} // End of namespace


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//==start_of_specific_mesh_class==========================================
/// Two layer mesh which employs a pseudo-solid node-update strategy.
/// This class is essentially a wrapper to an ElasticRectangularQuadMesh,
/// with additional members to keep track of which elements are in the
/// lower or upper fluid. The class also contains Vectors of (pointers to)
/// "root" elements in both fluids (for use with spatial adaptivity).
//========================================================================
template <class ELEMENT>
class ElasticRefineableTwoLayerMesh :
 public virtual ElasticRefineableRectangularQuadMesh<ELEMENT>
{

public:

 /// \short Constructor: Pass number of elements in x-direction, number of
 /// elements in y-direction in bottom and top layer, respectively,
 /// axial length and height of top and bottom layers, and pointer 
 /// to timestepper (defaults to Steady timestepper)
 ElasticRefineableTwoLayerMesh(const unsigned &nx, 
                               const unsigned &ny1,
                               const unsigned &ny2, 
                               const double &lx,
                               const double &h1,
                               const double &h2,
                               TimeStepper* time_stepper_pt=
                               &Mesh::Default_TimeStepper)
  : RectangularQuadMesh<ELEMENT>(nx,ny1+ny2,lx,h1+h2,
                                 false,time_stepper_pt),
    ElasticRectangularQuadMesh<ELEMENT>(nx,ny1+ny2,lx,h1+h2,
                                        false,time_stepper_pt),
    ElasticRefineableRectangularQuadMesh<ELEMENT>(nx,ny1+ny2,lx,h1+h2,
                                                  false,time_stepper_pt)
  {
   // ----------------------------------------------------
   // Convert all nodes on the interface to boundary nodes
   // ----------------------------------------------------

   // Set the number of boundaries to 5
   this->set_nboundary(5);

   // Loop over horizontal elements
   for(unsigned e=0;e<nx;e++)
    {
     // Get pointer to element in lower fluid adjacent to interface
     FiniteElement* el_pt = this->finite_element_pt(nx*(ny1-1)+e);

     // Determine number of nodes in this element
     const unsigned n_node = el_pt->nnode();

     // The last three nodes in this element are those on the interface.
     // Loop over these nodes and convert them to boundary nodes.
     for(unsigned n=0;n<3;n++)
      {
       Node* nod_pt = el_pt->node_pt(n_node-3+n);
       this->convert_to_boundary_node(nod_pt);
       this->add_boundary_node(4,nod_pt);
      }
    }

   // -------------------------------------------------
   // Populate Vectors of (pointers to) "root" elements
   // -------------------------------------------------

   // Determine number of elements in lower and upper fluids
   const unsigned long n_lower = nx*ny1;
   const unsigned long n_upper = nx*ny2;
   
   // Add all the bulk elements in the lower fluid to the vector of
   // "root" elements of the lower layer
   Lower_layer_root_element_pt.reserve(n_lower);
   for(unsigned e=0;e<n_lower;e++)
    {
     Lower_layer_root_element_pt.push_back(this->finite_element_pt(e));
    }
   
   // Add all the bulk elements in the upper fluid to the vector of
   // "root" elements of the upper layer
   Upper_layer_root_element_pt.reserve(n_upper);
   for(unsigned e=n_lower;e<(n_lower+n_upper);e++)
    {
     Upper_layer_root_element_pt.push_back(this->finite_element_pt(e));
    }
   
   // Update the "current" lower/upper element pt vectors
   this->update_lower_and_upper_element_pt_vectors();

  } // End of constructor

 /// \short Update the lower and upper layer element pointer vectors.
 /// This needs to be called when the mesh is adapted so that they
 /// contain pointers to all newly created elements (and pointers to
 /// now deleted elements are removed)
 void update_lower_and_upper_element_pt_vectors()
 {
  // Clear the element pt vectors
  Lower_layer_element_pt.clear();
  Upper_layer_element_pt.clear();

  // Determine number of "root" elements in lower/upper layers
  const unsigned n_root_lower = Lower_layer_root_element_pt.size();
  const unsigned n_root_upper = Upper_layer_root_element_pt.size();

  // Loop over the lower layer "root" elements and call function which
  // will add the "leaf" elements to Lower_layer_element_pt
  for(unsigned e=0;e<n_root_lower;e++)
   {
    get_sons_recursively(Lower_layer_root_element_pt[e],
                         Lower_layer_element_pt);
   }

  // Loop over the upper layer "root" elements and call function which
  // will add the "leaf" elements to Upper_layer_element_pt
  for(unsigned e=0;e<n_root_upper;e++)
   {
    get_sons_recursively(Upper_layer_root_element_pt[e],
                         Upper_layer_element_pt);
   }

 } // End of update_lower_and_upper_element_pt_vectors

 /// Access function for number of elements in lower layer
 unsigned long nlower() const { return Lower_layer_element_pt.size(); }

 /// Access function for number of elements in upper layer
 unsigned long nupper() const { return Upper_layer_element_pt.size(); }

 /// Access function for pointers to elements in bottom layer
 FiniteElement* &lower_layer_element_pt(const unsigned long &i) 
  {
   return Lower_layer_element_pt[i];
  }

 /// Access function for pointers to elements in upper layer
 FiniteElement* &upper_layer_element_pt(const unsigned long &i) 
  {
   return Upper_layer_element_pt[i];
  }


private:

 /// Vector of pointers to elements in the lower layer
 Vector<FiniteElement*> Lower_layer_element_pt;

 /// Vector of pointers to elements in the upper layer
 Vector<FiniteElement*> Upper_layer_element_pt;
 
 /// \short Vector of pointers to those elements in the lower
 /// layer which were created with the initial mesh. Since the mesh
 /// can never unrefine past these they can be used to update
 /// Lower_layer_element_pt after mesh adaptation.
 Vector<FiniteElement*> Lower_layer_root_element_pt;

 /// \short Vector of pointers to those elements in the upper
 /// layer which were created with the initial mesh. Since the mesh
 /// can never unrefine past these they can be used to update
 /// Upper_layer_element_pt after mesh adaptation.
 Vector<FiniteElement*> Upper_layer_root_element_pt;

 /// \short Helper function to recursively get the sons of the element
 /// pointed to by father_pt and, if they are a "leaf" of the quadtree,
 /// add them to the vector of pointers sons_which_are_leaves_pt
 void get_sons_recursively(
  FiniteElement* const &father_pt, 
  Vector<FiniteElement*> &sons_which_are_leaves_pt)
  {
   // Cast to local element
   ELEMENT* local_element_pt = dynamic_cast<ELEMENT*>(father_pt);
   
   // If the father has no sons, we're done
   if(local_element_pt->quadtree_pt()->nsons()==0)
    {
     sons_which_are_leaves_pt.push_back(local_element_pt);
    }
   // Otherwise, call the recursion on all four sons
   else
    {
     using namespace QuadTreeNames;
     
     get_sons_recursively(local_element_pt->quadtree_pt()
                          ->son_pt(NW)->object_pt(),
                          sons_which_are_leaves_pt);
     get_sons_recursively(local_element_pt->quadtree_pt()
                          ->son_pt(NE)->object_pt(),
                          sons_which_are_leaves_pt);
     get_sons_recursively(local_element_pt->quadtree_pt()
                          ->son_pt(SE)->object_pt(),
                          sons_which_are_leaves_pt);
     get_sons_recursively(local_element_pt->quadtree_pt()
                          ->son_pt(SW)->object_pt(),
                           sons_which_are_leaves_pt);
    }
  } // End of get_sons_recursively

}; // End of specific mesh class


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//==start_of_problem_class================================================
/// Axisymmetric two fluid interface problem in a rectangular domain
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
class InterfaceProblem : public Problem
{

public:
 
 /// Constructor: Pass the number of elements and the width of the
 /// domain in the r direction. Also pass the number of elements in
 /// the z direction of the bottom (fluid 1) and top (fluid 2) layers,
 /// along with the heights of both layers.
 InterfaceProblem(const unsigned &n_r, const unsigned &n_z1, 
                  const unsigned &n_z2, const double &l_r, 
                  const double &h1, const double &h2);

 /// Destructor (empty)
 ~InterfaceProblem() {}

 /// Set initial conditions
 void set_initial_condition();

 /// Set boundary conditions
 void set_boundary_conditions();

 /// Doc the solution
 void doc_solution(DocInfo &doc_info);

 /// Do unsteady run up to maximum time t_max with given timestep dt
 void unsteady_run(const double &t_max, const double &dt); 

private:
 
 /// No actions required before solve step
 void actions_before_newton_solve() {}
 
 /// No actions required after solve step
 void actions_after_newton_solve() {}

 /// \short Actions before the timestep: For maximum stability, reset
 /// the current nodal positions to be the "stress-free" ones.
 void actions_before_implicit_timestep()
  {
   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
  }

 /// Strip off the interface elements before adapting the bulk mesh
 void actions_before_adapt();

 /// Rebuild the mesh of interface elements after adapting the bulk mesh
 void actions_after_adapt();

 /// Create the 1d interface elements
 void create_interface_elements();

 /// Delete the 1d interface elements
 void delete_interface_elements();

 /// Deform the mesh/free surface to a prescribed function
 void deform_free_surface(const double &epsilon, const double &k);
 
 /// \short Helper function to recursively get the NW and NE sons of the
 /// element pointer to by father_pt and, if they are a "leaf" of the
 /// quadtree, add them to the vector of pointers sons_which_are_leaves_pt
 void get_neighbour_sons_recursively(
  FiniteElement* const &father_pt,
  Vector<FiniteElement*> &sons_which_are_leaves_pt);

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e,
                   const unsigned &pdof, 
                   const double &pvalue)
  {
   // Fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  }

 /// Pointer to the (specific) "bulk" mesh
 ElasticRefineableTwoLayerMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

 // Pointer to the constitutive law used to determine the mesh deformation
 ConstitutiveLaw* Constitutive_law_pt;

 /// \short Vector of pointers to "root" elements (in the lower layer)
 /// which are adjacent to the interface
 Vector<FiniteElement*> Root_neighbour_element_pt;

 /// Width of domain
 double Lr;

 /// Trace file
 ofstream Trace_file;

}; // End of problem class



//==start_of_constructor==================================================
/// Constructor for axisymmetric two fluid interface problem
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
InterfaceProblem<ELEMENT,TIMESTEPPER>::
InterfaceProblem(const unsigned &n_r, const unsigned &n_z1,
                 const unsigned &n_z2, const double &l_r,
                 const double& h1, const double &h2) : Lr(l_r)
{

 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new TIMESTEPPER); 

 // Build and assign "bulk" mesh
 Bulk_mesh_pt = new ElasticRefineableTwoLayerMesh<ELEMENT>
  (n_r,n_z1,n_z2,l_r,h1,h2,time_stepper_pt());

 // Create and set the error estimator for spatial adaptivity
 Bulk_mesh_pt->spatial_error_estimator_pt() = new Z2ErrorEstimator;

 // Set the maximum refinement level for the mesh to 4
 Bulk_mesh_pt->max_refinement_level() = 4;

 // Create the "surface" mesh that will contain only the interface
 // elements. The constructor just creates the mesh without giving
 // it any elements, nodes, etc.
 Surface_mesh_pt = new Mesh;
 
 // Create interface elements at the boundary between the two fluids,
 // and add them to the surface mesh
 create_interface_elements();

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all sub-meshes into a single mesh
 build_global_mesh();

 // --------------------------------------------
 // Set the boundary conditions for this problem
 // --------------------------------------------

 // All nodes are free by default -- just pin the ones that have
 // Dirichlet conditions here

 // Determine number of mesh boundaries
 const unsigned n_boundary = Bulk_mesh_pt->nboundary();

 // Loop over mesh boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine number of nodes on boundary b
   const unsigned n_node = Bulk_mesh_pt->nboundary_node(b);

   // Loop over nodes on boundary b
   for(unsigned n=0;n<n_node;n++)
    {
     // Fluid boundary conditions:
     // --------------------------

     // Pin radial and azimuthal velocities (no slip/penetration)
     // on all boundaries other than the interface (b=4)
     if(b!=4)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(2);
      }

     // Pin axial velocity on top (b=2) and bottom (b=0) boundaries
     // (no penetration). Because we have a slippery outer wall we do
     // NOT pin the axial velocity on this boundary (b=1); similarly,
     // we do not pin the axial velocity on the symmetry boundary (b=3).
     if(b==0 || b==2) { Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1); }

     // Solid boundary conditions:
     // --------------------------

     // Pin vertical displacement on solid boundaries
     if(b==0 || b==2) { Bulk_mesh_pt->boundary_node_pt(b,n)->pin_position(1); }

    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries

 // Loop over all nodes in mesh
 const unsigned n_node = Bulk_mesh_pt->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   // Pin horizontal displacement of all nodes
   Bulk_mesh_pt->node_pt(n)->pin_position(0);

   // Pin all azimuthal velocities throughout the bulk of the domain
   Bulk_mesh_pt->node_pt(n)->pin(2);
  }

 // Define a constitutive law for the solid equations: generalised Hookean
 Constitutive_law_pt = new GeneralisedHookean(&Global_Physical_Variables::Nu);

 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------

 // Determine number of bulk elements in lower/upper fluids
 const unsigned n_lower = Bulk_mesh_pt->nlower();
 const unsigned n_upper = Bulk_mesh_pt->nupper();

 // Loop over bulk elements in lower fluid
 for(unsigned e=0;e<n_lower;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                           lower_layer_element_pt(e));

   // Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;

   // Set the product of the Reynolds number and the inverse of the
   // Froude number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

   // Set the direction of gravity
   el_pt->g_pt() = &Global_Physical_Variables::G;

   // Set the constitutive law
   el_pt->constitutive_law_pt() = Constitutive_law_pt;

   // Assign the time pointer
   el_pt->time_pt() = time_pt();

  } // End of loop over bulk elements in lower fluid

 // Loop over bulk elements in upper fluid 
 for(unsigned e=0;e<n_upper;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                           upper_layer_element_pt(e));

   // Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;

   // Set the product of the Reynolds number and the inverse of the
   // Froude number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

   // Set the direction of gravity
   el_pt->g_pt() = &Global_Physical_Variables::G;

   // Set the viscosity ratio
   el_pt->viscosity_ratio_pt() = &Global_Physical_Variables::Viscosity_Ratio;

   // Set the density ratio
   el_pt->density_ratio_pt() = &Global_Physical_Variables::Density_Ratio;

   // Set the constitutive law
   el_pt->constitutive_law_pt() = Constitutive_law_pt;

   // Assign the time pointer
   el_pt->time_pt() = time_pt();

  } // End of loop over bulk elements in upper fluid

 // Set the pressure in the first element at 'node' 0 to 0.0
 fix_pressure(0,0,0.0);

 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  Bulk_mesh_pt->element_pt());

 // Apply the boundary conditions
 set_boundary_conditions();

 // Now we set up the information about which of the "root" elements (i.e.
 // those that were created with the original mesh) in the lower layer
 // are adjacent to the interface. Pointers to these elements are stored
 // in a Vector which is then used by create_interface_elements() as a
 // starting point for determining, at any particular state of mesh
 // adaptation, which of the current elements are adjacent to the interface.
 // We use this as a starting point since ONLY "root" elements are
 // guarenteed to exist at any particular state of mesh adaptation.
 // PATRICKFLAG MOVE ABOVE INTO DESCRIPTION IN DOCUMENTATION

 // ----------------------------------------------------------------
 // Populate vector which stores (pointers to) "root" elements in
 // the lower layer which are adjacent to the interface
 // ----------------------------------------------------------------

 // Loop over horizontal elements
 for(unsigned e=0;e<n_r;e++)
  {
   // Add those elements adjacent to the interface to the vector of
   // pointers to neighbouring "root" elements
   Root_neighbour_element_pt.push_back(Bulk_mesh_pt->
                                       finite_element_pt(n_r*(n_z1-1)+e));
  }
 
 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // End of constructor



//==start_of_set_initial_condition========================================
/// \short Set initial conditions: Set all nodal velocities to zero and
/// initialise the previous velocities and nodal positions to correspond
/// to an impulsive start
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::set_initial_condition()
{
 // Determine number of nodes in mesh
 const unsigned n_node = Bulk_mesh_pt->nnode();
 
 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Loop over the three velocity components
   for(unsigned i=0;i<3;i++)
    {
     // Set velocity component i of node n to zero
     Bulk_mesh_pt->node_pt(n)->set_value(i,0.0);
    }
  }
 
 // Initialise the previous velocity values and nodal positions
 // for timestepping corresponding to an impulsive start
 assign_initial_values_impulsive();
 
} // End of set_initial_condition



//==start_of_set_boundary_conditions======================================
/// \short Set boundary conditions: Set all velocity components to zero
/// on the top and bottom (solid) walls and the radial and azimuthal
/// components only to zero on the side boundaries
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::set_boundary_conditions()
{
 // Determine number of mesh boundaries
 const unsigned n_boundary = Bulk_mesh_pt->nboundary();
 
 // Loop over mesh boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine number of nodes on boundary b
   const unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
   
   // Loop over nodes on boundary b
   for(unsigned n=0;n<n_node;n++)
    {
     // Set radial component of the velocity to zero on all boundaries
     // other than the interface (b=4)
     if(b!=4) { Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(0,0.0); }

     // Set azimuthal component of the velocity to zero on all boundaries
     // other than the interface (b=4)
     if(b!=4) { Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(2,0.0); }

     // Set axial component of the velocity to zero on solid boundaries
     if(b==0 || b==2)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(1,0.0);
      }
    }
  }
} // End of set_boundary_conditions



//==start_of_actions_before_adapt=========================================
/// Strip off the interface elements before adapting the bulk mesh
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::actions_before_adapt()
{
 // Delete the interface elements and wipe the surface mesh
 delete_interface_elements();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();

} // End of actions_before_adapt



//==start_of_actions_after_adapt==========================================
/// Rebuild the mesh of interface elements after adapting the bulk mesh
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::actions_after_adapt()
{
 // Update the vectors of pointers to elements in the lower/upper
 // fluid layers so that they remain consistent after adaptation
 Bulk_mesh_pt->update_lower_and_upper_element_pt_vectors();
 
 // Create the interface elements
 this->create_interface_elements();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
 
 // Pin horizontal displacement of all nodes
 const unsigned n_node = Bulk_mesh_pt->nnode();
 for(unsigned n=0;n<n_node;n++) { Bulk_mesh_pt->node_pt(n)->pin_position(0); }

 // Unpin all fluid pressure dofs
 RefineableAxisymmetricNavierStokesEquations::
  unpin_all_pressure_dofs(Bulk_mesh_pt->element_pt());
 
 // Pin redudant fluid pressure dofs
 RefineableAxisymmetricNavierStokesEquations::
  pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());
 
 // Now set the pressure in the first element at 'node' 0 to 0.0
 fix_pressure(0,0,0.0);
 
 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  Bulk_mesh_pt->element_pt());

 // Reset the boundary conditions
 set_boundary_conditions();

} // End of actions_after_adapt



//==start_of_create_interface_elements====================================
/// \short Create interface elements between the two fluids in the mesh
/// pointed to by Bulk_mesh_pt and add the elements to the Mesh object
/// pointed to by Surface_mesh_pt.
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::create_interface_elements()
{
 // ---------------------------------------------------------------
 // Determine which bulk elements (in the lower fluid) are adjacent
 // to the interface (we will call these "neighbours")
 // ---------------------------------------------------------------

 // Create storage for the neighbours
 Vector<FiniteElement*> neighbour_pt;

 // Determine the number of "root" neighbours (i.e. those elements which
 // were neighbours in the original mesh, before mesh adaptation)
 const unsigned n_root_neighbour = Root_neighbour_element_pt.size();

 // Loop over the "root" neighbours to determine current neighbours
 for(unsigned e=0;e<n_root_neighbour;e++)
  {
   get_neighbour_sons_recursively(Root_neighbour_element_pt[e],neighbour_pt);
  }
  
 // --------------------------------------------------------------
 // Now loop over the neighbours and create the interface elements
 // --------------------------------------------------------------

 // Determine the number of neighbours
 const unsigned n_neighbour = neighbour_pt.size();

 // Loop over the neighbours
 for(unsigned e=0;e<n_neighbour;e++)
  {
   // Create the interface element (on face 2 of the bulk element)
   FiniteElement* interface_element_element_pt =
    new ElasticAxisymmetricFluidInterfaceElement<ELEMENT>(neighbour_pt[e],2);

   // Add the interface element to the surface mesh
   this->Surface_mesh_pt->add_element_pt(interface_element_element_pt);
  }

 // --------------------------------------------------------
 // Complete the setup to make the elements fully functional
 // --------------------------------------------------------

 // Determine number of 1D interface elements in mesh
 const unsigned n_interface_element = this->Surface_mesh_pt->nelement();

 // Loop over the interface elements
 for(unsigned e=0;e<n_interface_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ElasticAxisymmetricFluidInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<ElasticAxisymmetricFluidInterfaceElement<ELEMENT>*>
    (Surface_mesh_pt->element_pt(e));

   // Set the Strouhal number
   el_pt->st_pt() = &Global_Physical_Variables::St;

   // Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;

  } // End of loop over interface elements

} // End of create_interface_elements()



//==start_of_delete_interface_elements====================================
/// Delete the interface elements and wipe the surface mesh
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::delete_interface_elements()
{
 // Determine number of interface elements
 const unsigned n_interface_element = Surface_mesh_pt->nelement();
 
 // Loop over interface elements and delete
 for(unsigned e=0;e<n_interface_element;e++)
  {
   delete Surface_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 Surface_mesh_pt->flush_element_and_node_storage();
 
} // End of delete_interface_elements



//==start_of_deform_free_surface==========================================
/// Deform the mesh/free surface to a prescribed function
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::
deform_free_surface(const double &epsilon,const double &k)
{
 // Initialise Bessel functions (only need the first!)
 double j0, j1, y0, y1, j0p, j1p, y0p, y1p;

 // Determine number of nodes in the "bulk" mesh
 const unsigned n_node = Bulk_mesh_pt->nnode();
 
 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Determine eulerian position of node
   const double current_r_pos = Bulk_mesh_pt->node_pt(n)->x(0);
   const double current_z_pos = Bulk_mesh_pt->node_pt(n)->x(1);
   
   // Compute Bessel functions
   CRBond_Bessel::bessjy01a(k*current_r_pos,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
   
   // Determine new vertical position of node
   const double new_z_pos = current_z_pos
    + (1.0-fabs(1.0-current_z_pos))*epsilon*j0;
   
   // Set new position
   Bulk_mesh_pt->node_pt(n)->x(1) = new_z_pos;
  }
} // End of deform_free_surface



//==start_of_get_neighbour_sons_recursively===============================
/// \short Helper function to recursively get the NW and NE sons of the
/// element pointer to by father_pt and, if they are a "leaf" of the
/// quadtree, add them to the vector of pointers sons_which_are_leaves_pt
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::get_neighbour_sons_recursively(
 FiniteElement* const &father_pt,
 Vector<FiniteElement*> &sons_which_are_leaves_pt)
{
 // Upcast from GeneralisedElement to the present element
 ELEMENT* local_element_pt = dynamic_cast<ELEMENT*>(father_pt);

 // If the father has no sons, we're done
 if(local_element_pt->quadtree_pt()->nsons()==0)
  {
   sons_which_are_leaves_pt.push_back(local_element_pt);
  }
 // Otherwise call the recursion on the NW and NE sons
 else
  {
   using namespace QuadTreeNames;
   get_neighbour_sons_recursively(local_element_pt->quadtree_pt()
                                  ->son_pt(NW)->object_pt(),
                                  sons_which_are_leaves_pt);
   get_neighbour_sons_recursively(local_element_pt->quadtree_pt()
                                  ->son_pt(NE)->object_pt(),
                                  sons_which_are_leaves_pt);
  }
} // End of get_neighbour_sons_recursively



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::doc_solution(DocInfo &doc_info)
{ 

 // Output the time
 cout << "Time is now " << time_pt()->time() << std::endl;

 // Upcast from GeneralisedElement to the present element
 ElasticAxisymmetricFluidInterfaceElement<ELEMENT>* el_pt = 
  dynamic_cast<ElasticAxisymmetricFluidInterfaceElement<ELEMENT>*>
  (Surface_mesh_pt->element_pt(0));
 
 // Document time and vertical position of left hand side of interface
 // in trace file
 Trace_file << time_pt()->time() << " "
            << el_pt->node_pt(0)->x(1) << std::endl;
 
 ofstream some_file;
 char filename[100];
 
 // Set number of plot points (in each coordinate direction)
 const unsigned npts = 5;
 
 // Open solution output file
 sprintf(filename,"%s/soln%i.dat",
         doc_info.directory().c_str(),doc_info.number());
 some_file.open(filename);

 // Output solution to file
 Bulk_mesh_pt->output(some_file,npts);

 // Close solution output file
 some_file.close();

 // Open interface solution output file
 sprintf(filename,"%s/interface_soln%i.dat",
         doc_info.directory().c_str(),doc_info.number());
 some_file.open(filename);
 
 // Output solution to file
 Surface_mesh_pt->output(some_file,npts);
 
 // Close solution output file
 some_file.close();
 
} // End of doc_solution



//==start_of_unsteady_run=================================================
/// Perform run up to specified time t_max with given timestep dt
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::
unsteady_run(const double &t_max, const double &dt)
{

 // Set value of epsilon
 const double epsilon = 0.1;
 
 // Set value of k in Bessel function J_0(kr)
 const double k_bessel = 3.8317;

 // Deform the mesh/free surface
 deform_free_surface(epsilon,k_bessel);

 // Initialise DocInfo object
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Initialise counter for solutions
 doc_info.number()=0;

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 // Initialise trace file
 Trace_file << "time, interface height" << std::endl;

 // Initialise timestep
 initialise_dt(dt);

 // Set initial condition
 set_initial_condition();
 
 // Maximum number of spatial adaptations per timestep
 unsigned max_adapt = 2;

 // Call refine_uniformly twice
 for(unsigned i=0;i<2;i++) { refine_uniformly(); }

 // Doc initial solution
 doc_solution(doc_info);

 // Increment counter for solutions
 doc_info.number()++;

 // Determine number of timesteps
 const unsigned n_timestep = unsigned(t_max/dt);
 
 // Are we on the first timestep? At this point, yes!
 bool first_timestep = true;

 // Timestepping loop
 for(unsigned t=1;t<=n_timestep;t++)
  {
   // Output current timestep to screen
   cout << "\nTimestep " << t << " of " << n_timestep << std::endl;
   
   // Take one fixed timestep with spatial adaptivity
   unsteady_newton_solve(dt,max_adapt,first_timestep);

   // No longer on first timestep, so set first_timestep flag to false
   first_timestep = false; 

   // Reset maximum number of adaptations for all future timesteps
   max_adapt = 1;

   // Doc solution
   doc_solution(doc_info);

   // Increment counter for solutions 
   doc_info.number()++;

  } // End of timestepping loop

} // End of unsteady_run


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//==start_of_main=========================================================
/// Driver code for axisymmetric two fluid interface problem
//========================================================================
int main(int argc, char* argv[]) 
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Compute the Womersley number
 Global_Physical_Variables::ReSt =
  Global_Physical_Variables::Re*Global_Physical_Variables::St;

 /// Maximum time
 double t_max = 1.2;

 /// Duration of timestep
 const double dt = 0.005;

 // If we are doing validation run, use smaller number of timesteps
 if(CommandLineArgs::Argc>1) { t_max = 0.01; }

 // Number of elements in radial (r) direction
 const unsigned n_r = 3;
   
 // Number of elements in axial (z) direction in bottom fluid (fluid 1)
 const unsigned n_z1 = 3;
   
 // Number of elements in axial (z) direction in top fluid (fluid 2)
 const unsigned n_z2 = 3;

 // Width of domain
 const double l_r = 1.0;

 // Height of lower fluid layer
 const double h1 = 1.0;

 // Height of upper fluid layer
 const double h2 = 1.0;

 // Set direction of gravity (vertically downwards)
 Global_Physical_Variables::G[0] = 0.0;
 Global_Physical_Variables::G[1] = -1.0;
 Global_Physical_Variables::G[2] = 0.0;

 // Set up the spine test problem with AxisymmetricQCrouzeixRaviartElements,
 // using the BDF<2> timestepper
 InterfaceProblem<RefineablePseudoSolidNodeUpdateElement<
 RefineableAxisymmetricQCrouzeixRaviartElement,
  RefineableQPVDElement<2,3> >,BDF<2> >
  problem(n_r,n_z1,n_z2,l_r,h1,h2);
   
 // Run the unsteady simulation
 problem.unsteady_run(t_max,dt);

} // End of main
