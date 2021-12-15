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
//Driver for a multi-physics problem that couples the Axisymmetric Navier--Stokes
//equations to the axisymmetric advection diffusion equations, 
//in a cilynder jet with Robin boundary (temperature) at bound 1

//Oomph-lib headers, we require the generic, advection-diffusion
//and navier-stokes elements.

// Generic oomph-lib header
#include "generic.h"

// Axisymmetric Navier-Stokes headers
#include "axisym_navier_stokes.h"

// Interface headers
#include "fluid_interface.h"

// The Axisymetric Advection Diffusion equations
#include "steady_axisym_advection_diffusion.h"

// The mesh
#include "meshes/horizontal_single_layer_spine_mesh.h"

// Use the oomph and std namespaces 
using namespace oomph;

using namespace std;

// Megaclass header

#include "non_isothermal_axisym_nst_elements.h"

//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{
 /// ----- KINEMATIC AND DYNAMIC PROBLEM -----

 /// Reynolds number
 double Re = 0.01; // (Re = epsilon*R)

 /// Product of Reynolds number and inverse of Froude number
 double ReInvFr = 0.01; // (Re/Fr = epsilon^2*R/F -> 1/F = 0)

 /// Capillary number
 double Ca = 10.0; // (Ca = C/epsilon)

 /// External pressure
 double P_ext = 0.0;

 /// Gravity vector
 Vector<double> Direction_of_gravity(3);

 /// ----- ENERGY PROBLEM -----

 /// Peclet number 
 double Pe = 0.0; // (Pe = epsilon*P)

 /// Biot number
 double Bi = 0.01; // (Bi = epsilon^2*B)

 /// Exterior temperature
 double T_ext = 0.0;

 /// Inlet temperature
 double T_inlet = 1.0;

 /// ----- VISCOSITY PARAMETERS -----

 /// Pre-exponetial factor 
 double G = 1.0;

 /// Eta factor (exponent of the exponecial function)
 double eta = 0.0;

 /// Viscosity ratio function modeled following a Arrhenius fashion
 void viscosity_ratio_function(double& temperature,
                               double& result)
 {
  result = G*exp(eta*(T_inlet-temperature));
 }

 /// Beta on a boundary on which r is fixed
 void prescribed_beta_on_fixed_r_boundary(const Vector<double>& x_vector, 
                                          double& beta)
 {
  beta = -Bi*T_ext; 
 }

 /// Alfa on a boundary on which r is fixed
 void prescribed_alpha_on_fixed_r_boundary(const Vector<double>& x_vect, 
                                           double& alpha)
 {
  // Use alpha<0 for put T_inf<u(r,z)<T_inlet
  // alpha = -Bi and beta = -Bi*T_inf with Bi>0
  alpha = -Bi; 
 }

} // end_of_namespace

/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////

//====== start_of_problem_class=======================================
/// 2D Convection  problem on rectangular domain, discretised 
/// with refineable elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class AxisymFreeSurfaceNozzleAdvDiffRobinProblem : public Problem
{

public:

 /// Constructor: Pass the number of elements and the lengths of the
 /// domain in the r and z directions (h is the height of the fluid layer
 /// i.e. the length of the domain in the z direction)
 AxisymFreeSurfaceNozzleAdvDiffRobinProblem(const unsigned &n_r, 
                                            const unsigned &n_z, 
                                            const double &l_r, 
                                            const double &h);

 /// Destructor. Empty
 ~AxisymFreeSurfaceNozzleAdvDiffRobinProblem() {}

 /// Spine heights/lengths are unknowns in the problem so their values get
 /// corrected during each Newton step. However, changing their value does
 /// not automatically change the nodal positions, so we need to update all
 /// of them here.
 void actions_before_newton_convergence_check()
  {
   Bulk_mesh_pt->node_update();
  }

 /// Doc the solution.
 void doc_solution(DocInfo& doc_info);

 /// Do steady run up to maximum time t_max with given timestep dt
 void steady_run();

private:

 /// Update the problem specs before solve
 /// Re-set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve(); 

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}
  
 /// Create Axisymmetric Advection Diffusion flux elements on boundary b of 
 /// the Mesh pointed to by bulk_mesh_pt and add them to the Mesh 
 /// object pointed to by surface_mesh_pt
 void create_flux_elements(const unsigned &b, 
                           Mesh* const &bulk_mesh_pt,
                           Mesh* const &surface_mesh_pt);

 /// Delete Axisymmetric Advection Diffusion flux elements and wipe the surface mesh
 void delete_flux_elements(Mesh* const &surface_mesh_pt);

 /// Pointer to the "bulk" mesh
 HorizontalSingleLayerSpineMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh of flux elements
 Mesh* Surface_mesh_pt;

 /// Pointer to the mesh of interface elements
 Mesh* Interface_mesh_pt;

 /// Deform the mesh/free surface to a prescribed function
 void deform_free_surface(const double &Dr)
  {
   // Determine number of spines in mesh
   const unsigned n_spine = Bulk_mesh_pt->nspine();

   // Loop over spines in mesh
   for(unsigned i=0;i<n_spine;i++)
    {
     // Determine z coordinate of spine
     double z_value = Bulk_mesh_pt->boundary_node_pt(1,i)->x(1);

     if (z_value<=(Height/1.1))
      {
       // Set spine height
       Bulk_mesh_pt->spine_pt(i)->height() = sqrt( exp(-log(Dr)*(1.0 - (1.1*z_value/Height))) );
      }
    } // End of loop over spines
   
   // Update nodes in bulk mesh
   Bulk_mesh_pt->node_update();

  } // End of deform_free_surface

 /// Pointer to viscosity ratio function 
 //ViscosityRatioFctPt Viscosity_ratio_fct_pt;

 /// Width of domain
 double Lr;

 /// Height of the domain
 double Height;

}; // end of problem class

//===========start_of_constructor=========================================
/// Constructor for convection problem
//========================================================================
template<class ELEMENT>
AxisymFreeSurfaceNozzleAdvDiffRobinProblem<ELEMENT>::
AxisymFreeSurfaceNozzleAdvDiffRobinProblem(const unsigned &n_r, 
                                           const unsigned &n_z,
                                           const double &l_r, 
                                           const double& h) : Lr(l_r),
                                                              Height(h)
{
 // Build and assign mesh
 Bulk_mesh_pt = 
  new HorizontalSingleLayerSpineMesh<ELEMENT>(n_r,n_z,l_r,h);

 //Create "surface mesh" that will only contain the interface elements
 Interface_mesh_pt = new Mesh;
 {
  // How many bulk elements are adjacent to boundary b?
  unsigned n_element = Bulk_mesh_pt->nboundary_element(1);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(1,e));
   
   // Find the index of the face of element e along boundary b
   int face_index = Bulk_mesh_pt->face_index_at_boundary(1,e);

   // Build the corresponding free surface element
   SpineAxisymmetricFluidInterfaceElement<ELEMENT>* interface_element_pt = new 
   SpineAxisymmetricFluidInterfaceElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the prescribed-flux element to the surface mesh
   Interface_mesh_pt->add_element_pt(interface_element_pt);

  } //end of loop over bulk elements adjacent to boundary b
 }


 // Create "surface mesh" that will contain only the prescribed-flux 
 // elements. The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Surface_mesh_pt = new Mesh;

 // Create prescribed-flux elements from all elements that are 
 // adjacent to boundary 1, but add them to a separate mesh.
 // Note that this is exactly the same function as used in the 
 // single mesh version of the problem, we merely pass different Mesh pointers.
 create_flux_elements(1,Bulk_mesh_pt,Surface_mesh_pt);

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Interface_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();

 // --------------------------------------------
 // Set the boundary conditions for this problem
 // --------------------------------------------

 // All nodes are free by default -- only need to pin the ones 
 // that have Dirichlet conditions here

 //Loop over the boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Loop over the number of nodes on the boundary
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   switch (ibound) 
    {
    case 0:
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //Pin all velocities (Outlet zone)
       //Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(2);
      }
     break;
    case 1:
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       double z_value = Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(1);
       //Pin all velocities (Wall)
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(2);
       if (z_value>=(Height/1.1))
        {
         Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
         Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
         //Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(3);
        }
      }
     break;
    case 2:
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //Pin all velocities (Inlet zone) and temperature
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(2);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(3);
      }   
     break;
    default:
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //Pin all velocities (axis of simmetry)
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(2);
      }    
     break;
    }
  }

 // Pin spine height at top of domain
 unsigned first_spine = Bulk_mesh_pt->nspine()-5;
 unsigned last_spine = Bulk_mesh_pt->nspine()-1;
 for (unsigned num_spine=first_spine;num_spine<=last_spine;num_spine++)
  {
   Bulk_mesh_pt->spine_pt(num_spine)->spine_height_pt()->pin(0);
  }
 
 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------
 
 //Find number of elements in mesh
 unsigned n_element = Bulk_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 // Loop over the bulk elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   
   // Set the product of the Reynolds number and the inverse of the
   // Froude number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

   // Set the direction of gravity
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;

   // Set the Reynolds number
   el_pt->pe_pt() = &Global_Physical_Variables::Pe;

   // Set the viscosity ratio functino
   el_pt->viscosity_ratio_fct_pt() = 
    &Global_Physical_Variables::viscosity_ratio_function;

  } // End of loop over elements

  // Create a Data object whose single value stores the external pressure
  Data* external_pressure_data_pt = new Data(1);
 
  // Pin and set the external pressure to some arbitrary value
  double p_ext = Global_Physical_Variables::P_ext;

  external_pressure_data_pt->pin(0);
  external_pressure_data_pt->set_value(0,p_ext);

 // Determine number of 1D interface elements in mesh
 const unsigned n_interface_element = Interface_mesh_pt->nelement();

 // Loop over the interface elements
 for(unsigned e=0;e<n_interface_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   SpineAxisymmetricFluidInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<SpineAxisymmetricFluidInterfaceElement<ELEMENT>*>
    (Interface_mesh_pt->element_pt(e));

   // Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;

   // Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(external_pressure_data_pt);

  } // End of loop over interface elements

  // Loop over the flux elements to pass pointer to prescribed flux function
   n_element = Surface_mesh_pt->nelement()-2; //2 elements at the nozzle

  for(unsigned e=0;e<n_element;e++)
   {
    // Upcast from GeneralisedElement to AdvectionDiffusion flux element
    SteadyAxisymAdvectionDiffusionFluxElement<ELEMENT> *el_pt = 
     dynamic_cast< SteadyAxisymAdvectionDiffusionFluxElement<ELEMENT>*>(
     Surface_mesh_pt->element_pt(e));
    if (e<n_element) //Free surface
     {
      // Set the pointer to the prescribed beta function
      el_pt->beta_fct_pt() = 
       &Global_Physical_Variables::prescribed_beta_on_fixed_r_boundary;
      // Set the pointer to the prescribed alpha function
      el_pt->alpha_fct_pt() = 
       &Global_Physical_Variables::prescribed_alpha_on_fixed_r_boundary;
     }
   }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor

//===================start_of_actions_before_newton_solve================
/// Update the problem specs before solve: (Re-)set boundary conditions
/// just to be on the safe side...
//========================================================================
template<class ELEMENT>
void AxisymFreeSurfaceNozzleAdvDiffRobinProblem<ELEMENT>::
actions_before_newton_solve() 
{
 const double Dr = 10.0;
 const double t_inlet = Global_Physical_Variables::T_inlet;

 // Determine number of nodes in mesh
 const unsigned n_node = Bulk_mesh_pt->nnode();

 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Determine r coordinate of node
   double r_value = Bulk_mesh_pt->node_pt(n)->x(0);
   // Determine z coordinate of node
   double z_value = Bulk_mesh_pt->node_pt(n)->x(1);
      
   // Initial guess for ur (Multiply by epsilon)
   double ur_value = -0.5*r_value*(1.0/Height)*log(Dr)*exp(log(Dr)*(1.0-(1.1*z_value/Height)));
   //Initial guess for uz
   double uz_value = -exp(log(Dr)*(1.0-(1.1*z_value/Height)));

   // Set velocity component i of node n to guess
   Bulk_mesh_pt->node_pt(n)->set_value(0,ur_value);
   Bulk_mesh_pt->node_pt(n)->set_value(1,uz_value);
   // Set theta velocity component of node n to zero
   Bulk_mesh_pt->node_pt(n)->set_value(2,0.0);
   // Set temperature of  node n to one
   Bulk_mesh_pt->node_pt(n)->set_value(3,t_inlet);
  }

 // Correct the value in order to imposse boundary conditions
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Determine number of nodes in the bound mesh
   const unsigned n_node = Bulk_mesh_pt->nboundary_node(ibound);
   switch (ibound) 
    {
    case 0:
     for (unsigned inod=0;inod<n_node;inod++)
      {
       // Determine r coordinate of node
       double r_value = Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(0);
       double w_bound0 = -Dr*(1.0 - 0.0*pow(r_value,2.0));
       //Set the velocities at the bottom zone
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,0.0);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,w_bound0);
       //Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(2,0.0);
       //Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(3,t_bound0);
      }
     break;
    case 1:
     for (unsigned inod=0;inod<n_node;inod++)
      {
       // Determine r coordinate of node
       double z_value = Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(1);
       //Set the azimuthal velocity at the outer wall
       //Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(2,0.0);
       if (z_value>=(Height/1.1))
        {
         Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,0.0);
         Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);
        }
      }
     break;
    case 2:
     for (unsigned inod=0;inod<n_node;inod++)
      {
       // Determine r coordinate of node
       double r_value = Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(0);
       double w_bound2 = -2.0*(1.0 - 1.0*pow(r_value,2.0));
       double t_bound2 = t_inlet - 0.0*pow(r_value,2.0);
       //Set all of the magnitudes at the top zone
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,0.0);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,w_bound2);
       //Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(2,0.0);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(3,t_bound2);
      }   
     break;
    default: // simmetry axis
     for (unsigned inod=0;inod<n_node;inod++)
      {
       //Set the radial and azimuthal velocity at the axis
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,0.0);
       //Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(2,0.0);
      }    
     break;
    }
   } 
}


//===============start_doc_solution=======================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void AxisymFreeSurfaceNozzleAdvDiffRobinProblem<ELEMENT>::
doc_solution(DocInfo& doc_info)
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts = 5;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);

 //Find number of elements in mesh
 unsigned n_element = Bulk_mesh_pt->nelement();

 // Loop over the bulk elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Set the Reynolds number
   el_pt->output(some_file,npts);
  }


// Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

} // end of doc_solution

//============start_of_create_flux_elements==============================
/// Create AdvectionDiffusion Flux Elements on the b-th boundary of 
/// the Mesh object pointed to by bulk_mesh_pt and add the elements 
/// to the Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void AxisymFreeSurfaceNozzleAdvDiffRobinProblem<ELEMENT>::
create_flux_elements(const unsigned &b, 
                     Mesh* const &bulk_mesh_pt,
                     Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   // Find the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-flux element
   SteadyAxisymAdvectionDiffusionFluxElement<ELEMENT>* flux_element_pt = new 
   SteadyAxisymAdvectionDiffusionFluxElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_flux_elements

//============start_of_delete_flux_elements==============================
/// Delete Advection Diffusion Flux Elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void AxisymFreeSurfaceNozzleAdvDiffRobinProblem<ELEMENT>::
delete_flux_elements(Mesh* const &surface_mesh_pt)
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

} // end of delete_flux_elements

//==start_of_steady_run=================================================
/// Perform run
//========================================================================
template<class ELEMENT>
void AxisymFreeSurfaceNozzleAdvDiffRobinProblem<ELEMENT>::
steady_run()
{
 // Increase maximum residual and iteration number
 Problem::Max_residuals=3000.0;
 Problem::Max_newton_iterations=200;

 // Set value of Dr
 const double Dr = 10.0;

 // Deform the mesh/free surface
 deform_free_surface(Dr);

 // Initialise DocInfo object
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Initialise counter for solutions
 doc_info.number()=0;
 
 // Doc solution
 doc_solution(doc_info);

 // Step number
 doc_info.number()++;

 for (unsigned i=0;i<2;i++)
  {

   cout << "Solving for Re = " 
        << Global_Physical_Variables::Re
        << " Ca = "
        << Global_Physical_Variables::Ca
        << " Pe = "
        << Global_Physical_Variables::Pe
        << " Bi = "
        << Global_Physical_Variables::Bi
        << std::endl;

   // Solve the problem
   newton_solve();
   
   // Doc solution
   doc_solution(doc_info);
   
   // Step number
   doc_info.number()++;

   // Bump up parameter
   Global_Physical_Variables::Pe+=1.0;
   Global_Physical_Variables::eta+=5.0;

  }

} // End of steady_run

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

//=======start_of_main================================================
/// Driver code for 2D Boussinesq convection problem
//====================================================================
int main(int argc, char **argv)
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // # of elements in r-direction
 unsigned n_r = 8; //20;

 // # of elements in z-direction
 unsigned n_z = 22;

 // Domain length in r-direction
 double l_r = 1.0;

 // Domain length in z-direction
 double h = 11.0;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 //FiniteElement::Accept_negative_jacobian = true;

 // Set direction of gravity (vertically downwards)
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] =-1.0;
 Global_Physical_Variables::Direction_of_gravity[2] = 0.0;

 //Construct our problem
 AxisymFreeSurfaceNozzleAdvDiffRobinProblem<SpineElement<NonIsothermalAxisymmetricQCrouzeixRaviartElement> > 
 problem(n_r,n_z,l_r,h);
 
 // Run the steady simulation
 problem.steady_run();

} // end of main
