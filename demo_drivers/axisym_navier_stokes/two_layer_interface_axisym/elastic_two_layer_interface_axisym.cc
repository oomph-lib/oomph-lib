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

 /// Womersley number (Reynolds x Strouhal)
 double ReSt = 5.0;

 /// Product of Reynolds number and inverse of Froude number
 double ReInvFr = 5.0;

 /// Ratio of viscosity in upper fluid to viscosity in lower
 /// fluid. Reynolds number etc. is based on viscosity in lower fluid.
 double Viscosity_Ratio = 0.1;

 /// Ratio of density in upper fluid to density in lower
 /// fluid. Reynolds number etc. is based on density in lower fluid.
 double Density_Ratio = 0.5;

 /// Capillary number
 double Ca = 0.01;

 /// Direction of gravity
 Vector<double> G(3);

 /// Pseudo-solid Poisson ratio
 double Nu = 0.1;

} // End of namespace


/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////


//==start_of_specific_mesh_class==========================================
/// Two layer mesh which employs a pseudo-solid node-update strategy.
/// This class is essentially a wrapper to an
/// ElasticRefineableRectangularQuadMesh, with an additional boundary
/// to represent the interface between the two fluid layers.
//========================================================================
template <class ELEMENT>
class ElasticRefineableTwoLayerMesh :
 public virtual ElasticRefineableRectangularQuadMesh<ELEMENT>
{

public:

 /// Constructor: Pass number of elements in x-direction, number of
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
    } // End of loop over horizontal elements

   // Set up the boundary element information
   this->setup_boundary_element_info();
  }
 
}; // End of specific mesh class


/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////


//==start_of_problem_class================================================
/// Axisymmetric two fluid interface problem in a rectangular domain
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
class InterfaceProblem : public Problem
{

public:
 
 /// Constructor
 InterfaceProblem();

 /// Destructor (empty)
 ~InterfaceProblem() {}

 /// Set initial conditions
 void set_initial_condition();

 /// Set boundary conditions
 void set_boundary_conditions();

 /// Document the solution
 void doc_solution(DocInfo &doc_info);

 /// Do unsteady run up to maximum time t_max with given timestep dt
 void unsteady_run(const double &t_max, const double &dt); 

private:
 
 /// No actions required before solve step
 void actions_before_newton_solve() {}
 
 /// No actions required after solve step
 void actions_after_newton_solve() {}

 /// Actions before the timestep: For maximum stability, reset
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
InterfaceProblem()
{
 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new TIMESTEPPER); 

 // Define number of elements in r direction
 const unsigned n_r = 3;
   
 // Define number of elements in z direction in lower fluid (fluid 1)
 const unsigned n_z1 = 3;

 // Define number of elements in z direction in upper fluid (fluid 2)
 const unsigned n_z2 = 3;

 // Define width of domain and store as class member data
 const double l_r = 1.0;
 this->Lr = l_r;

 // Define height of lower fluid layer
 const double h1 = 1.0;

 // Define height of upper fluid layer
 const double h2 = 1.0;

 // Build and assign the "bulk" mesh
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

 // Compute number of bulk elements in lower/upper fluids
 const unsigned n_lower = n_r*n_z1;
 const unsigned n_upper = n_r*n_z2;

 // Loop over bulk elements in lower fluid
 for(unsigned e=0;e<n_lower;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

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

  } // End of loop over bulk elements in lower fluid

 // Loop over bulk elements in upper fluid 
 for(unsigned e=n_lower;e<(n_lower+n_upper);e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

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

  } // End of loop over bulk elements in upper fluid

 // Set the pressure in the first element at 'node' 0 to 0.0
 fix_pressure(0,0,0.0);

 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  Bulk_mesh_pt->element_pt());

 // Apply the boundary conditions
 set_boundary_conditions();

 // Set up equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // End of constructor



//==start_of_set_initial_condition========================================
/// Set initial conditions: Set all nodal velocities to zero and
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
/// Set boundary conditions: Set all velocity components to zero
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
/// Create interface elements between the two fluids in the mesh
/// pointed to by Bulk_mesh_pt and add the elements to the Mesh object
/// pointed to by Surface_mesh_pt.
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::create_interface_elements()
{
 // Determine number of bulk elements adjacent to interface (boundary 4)
 const unsigned n_element = this->Bulk_mesh_pt->nboundary_element(4);
 
 // Loop over those elements adjacent to the interface
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to the interface
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    this->Bulk_mesh_pt->boundary_element_pt(4,e));

   // We only want to attach interface elements to the bulk elements
   // which are BELOW the interface, and so we filter out those above by
   // referring to the viscosity_ratio_pt
   if(bulk_elem_pt->viscosity_ratio_pt()
      !=&Global_Physical_Variables::Viscosity_Ratio)
    {
     // Find index of the face of element e that corresponds to the interface
     const int face_index = this->Bulk_mesh_pt->face_index_at_boundary(4,e);
     
     // Create the interface element
     FiniteElement* interface_element_element_pt =
      new ElasticAxisymmetricFluidInterfaceElement<ELEMENT>(bulk_elem_pt,
                                                            face_index);

     // Add the interface element to the surface mesh
     this->Surface_mesh_pt->add_element_pt(interface_element_element_pt);
    }
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


/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////


//==start_of_main=========================================================
/// Driver code for axisymmetric two fluid interface problem
//========================================================================
int main(int argc, char* argv[]) 
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // -----------------------------------------------------------------
 // Define possible command line arguments and parse the ones that
 // were actually specified
 // -----------------------------------------------------------------

 // Are we performing a validation run?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Check that definition of Womersley number is consistent with those
 // of the Reynolds and Strouhal numbers
 if(Global_Physical_Variables::ReSt !=
    Global_Physical_Variables::Re*Global_Physical_Variables::St)
  {
   std::ostringstream error_stream;
   error_stream << "Definition of Global_Physical_Variables::ReSt is "
                << "inconsistant with those\n"
                << "of Global_Physical_Variables::Re and "
                << "Global_Physical_Variables::St." << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
  }

 /// Maximum time
 double t_max = 1.2;

 /// Duration of timestep
 const double dt = 0.005;

 // If we are doing validation run, use smaller number of timesteps
 if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   t_max = 0.01;
  }

 // Set direction of gravity (vertically downwards)
 Global_Physical_Variables::G[0] = 0.0;
 Global_Physical_Variables::G[1] = -1.0;
 Global_Physical_Variables::G[2] = 0.0;

 // Set up the elastic test problem with AxisymmetricQCrouzeixRaviartElements,
 // using the BDF<2> timestepper
 InterfaceProblem<RefineablePseudoSolidNodeUpdateElement<
 RefineableAxisymmetricQCrouzeixRaviartElement,
  RefineableQPVDElement<2,3> >,BDF<2> > problem;
   
 // Run the unsteady simulation
 problem.unsteady_run(t_max,dt);

} // End of main
