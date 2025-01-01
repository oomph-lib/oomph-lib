//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
// Driver for fibre spinning axisymmetric fluid problem exiting from a nozzle
 
// Generic oomph-lib header
#include "generic.h"

// Navier-Stokes headers
#include "navier_stokes.h"

// Axisymmetric Navier-Stokes headers
#include "axisym_navier_stokes.h"

// Interface headers
#include "fluid_interface.h"

// The mesh
#include "meshes/horizontal_single_layer_spine_mesh.h"

using namespace std;

using namespace oomph;


//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re = 0.1; // (Re = epsilon*R)

 /// Product of Reynolds number and inverse of Froude number
 double ReInvFr = 0.01; // (Re/Fr = epsilon^2*R/F)

 /// Capillary number
 double Ca = 10.0; // (Ca = C/epsilon)

 /// External pressure
 double P_ext = 0.0;

 /// Direction of gravity
 Vector<double> G(3);

} // End of namespace


//==start_of_problem_class===============================================
/// Melt spinning axisymmetric fluid interface problem in rectangular domain
//=======================================================================
template<class ELEMENT>
class MeltSpinningProblem : public Problem
{
 
public:
 
 /// Constructor: Pass the number of elements and the lengths of the
 /// domain in the r and z directions (h is the height of the fluid layer
 /// i.e. the length of the domain in the z direction)
 MeltSpinningProblem(const unsigned &n_r, 
                     const unsigned &n_z, 
                     const double &l_r, 
                     const double &h);
 
 /// Destructor (empty)
 ~MeltSpinningProblem() {}

 /// Spine heights/lengths are unknowns in the problem so their values get
 /// corrected during each Newton step. However, changing their value does
 /// not automatically change the nodal positions, so we need to update all
 /// of them here.
 void actions_before_newton_convergence_check()
  {
   Bulk_mesh_pt->node_update();
  }

 /// Doc the solution
 void doc_solution(DocInfo &doc_info);

 /// Do steady run 
 void steady_run();

private:

 /// Update the problem specs before solve
 /// Re-set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve(); 

 /// Update after solve can remain empty, because the update 
 /// is performed automatically after every Newton step.
 void actions_after_newton_solve() {}

 /// Access function for the specific mesh
 HorizontalSingleLayerSpineMesh<ELEMENT>* Bulk_mesh_pt; 

 /// Mesh for the interface elements
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
       Bulk_mesh_pt->spine_pt(i)->height() = sqrt(exp(-log(Dr)*(1.0-(1.1*z_value/Height))));
      }

    } // End of loop over spines
   
   // Update nodes in bulk mesh
   Bulk_mesh_pt->node_update();

  } // End of deform_free_surface

 /// Width of domain
 double Lr;

 /// Height of the domain
 double Height;

}; // End of problem class



//==start_of_constructor==================================================
/// Constructor for melt spinning axisymmetric fluid problem 
//========================================================================
template<class ELEMENT>
MeltSpinningProblem<ELEMENT>::
MeltSpinningProblem(const unsigned &n_r, 
                    const unsigned &n_z,
                    const double &l_r, 
                    const double& h) : Lr(l_r),
                                       Height(h)

{

 // Build and assign mesh (the "false" boolean flag tells the mesh
 // constructor that the domain is not periodic in r)
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

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Interface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();


 // --------------------------------------------
 // Set the boundary conditions for this problem
 // --------------------------------------------

 // All nodes are free by default -- just pin the ones that have
 // Dirichlet conditions here

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
       //Pin all velocities (Wall))
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
        }
      }
     break;
    case 2:
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //Pin all velocities (Wall)
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(2);
      }   
     break;
    default: // simmetry axis
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //Pin all velocities (Wall)
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
 
 cout << "First spine is " 
      << first_spine
      << " and last spine is  "
      << last_spine
      << std::endl;

 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------
 
 // Determine number of bulk elements in mesh
 const unsigned n_bulk = Bulk_mesh_pt->nelement();

 // Loop over the bulk elements
 for(unsigned e=0;e<n_bulk;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the product of the Reynolds number and the inverse of the
   // Froude number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

   // Set the direction of gravity
   el_pt->g_pt() = &Global_Physical_Variables::G;

  } // End of loop over bulk elements

 // Create a Data object whose single value stores the external pressure
 Data* external_pressure_data_pt = new Data(1);

// Pin and set the external pressure to some arbitrary value
 double p_ext = Global_Physical_Variables::P_ext;
 
 // Pin and set the external pressure to some arbitrary value
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

 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // End of constructor



//===================start_of_actions_before_newton_solve================
/// Update the problem specs before solve: (Re-)set boundary conditions
/// just to be on the safe side...
//========================================================================
template<class ELEMENT>
void MeltSpinningProblem<ELEMENT>::
actions_before_newton_solve() 
{
 const double Dr = 20.0;

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
   double ur_value = -0.5*r_value*(1.0/Height)*log(Dr)*exp(log(Dr)*(1.0-(z_value/Height)));
   //Initial guess for uz
   double uz_value = -exp(log(Dr)*(1.0-(z_value/Height)));

   // Set velocity component i of node n to guess
   Bulk_mesh_pt->node_pt(n)->set_value(0,ur_value);
   Bulk_mesh_pt->node_pt(n)->set_value(1,uz_value);
   // Set theta velocity component of node n to zero
   Bulk_mesh_pt->node_pt(n)->set_value(2,0.0);
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
      }
     break;
    case 1:
     for (unsigned inod=0;inod<n_node;inod++)
      {
       // Determine r coordinate of node
       double z_value = Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(1);
       //Set the zero velocity at the nozzle
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
       //Set all of the magnitudes at the top zone
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,0.0);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,w_bound2);
      }   
     break;
    default: // simmetry axis
     for (unsigned inod=0;inod<n_node;inod++)
      {
       //Set the radial velocity at the axis
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,0.0);
      }    
     break;
    }
   } 
} // End of actions_before_newton_solve

   
//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void MeltSpinningProblem<ELEMENT>::
doc_solution(DocInfo &doc_info)
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Set number of plot points (in each coordinate direction)
 const unsigned npts = 5;

 // Open solution output file
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);

 // Output solution to file
 Bulk_mesh_pt->output(some_file,npts);
 Interface_mesh_pt->output(some_file,npts);

 // Close solution output file
 some_file.close();
 
} // End of doc_solution

//==start_of_steady_run=================================================
/// Perform run
//========================================================================
 template<class ELEMENT>
void MeltSpinningProblem<ELEMENT>::
steady_run()
{
 // Increase maximum residual and iteration number
 Problem::Max_residuals=1000.0;
 Problem::Max_newton_iterations=200;

 // Set value of Dr
 const double Dr = 20.0;

 // Deform the mesh/free surface
 deform_free_surface(Dr);

 // Initialise DocInfo object
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Initialise counter for solutions
 doc_info.number()=0;
 
 // Doc initial solution
 doc_solution(doc_info);

 // Increment counter for solutions 
 doc_info.number()++;

 for(unsigned i=0;i<1;i++)
  {

   cout << "Solving for Re = " 
        << Global_Physical_Variables::Re
        << " Ca = "
        << Global_Physical_Variables::Ca
        << std::endl;

   // Solve the problem
   newton_solve();
   
   // Doc solution
   doc_solution(doc_info);
   
   // Step number
   doc_info.number()++;

   // Bump up parameter
   Global_Physical_Variables::Re=1.0*Global_Physical_Variables::Re;
   Global_Physical_Variables::Ca=1.0*Global_Physical_Variables::Ca;

  } 

} // End of steady_run



/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for isothermal fibre spinning axisymmetric fluid problem 
//=====================================================================
int main(int argc, char* argv[]) 
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Number of elements in radial (r) direction
 const unsigned n_r = 20;
   
 // Number of elements in axial (z) direction
 const unsigned n_z = 22;

 // Width of domain
 const double l_r = 1.0;

 // Height of fluid layer
 const double h = 11.0;
 
 // Accept negative jacobian 

 //FiniteElement::Accept_negative_jacobian = true;

 // Set direction of gravity (vertically downwards)
 Global_Physical_Variables::G[0] = 0.0;
 Global_Physical_Variables::G[1] =-1.0;
 Global_Physical_Variables::G[2] = 0.0;
 
 // Set up the spine test problem with AxisymmetricQCrouzeixRaviartElements,
 MeltSpinningProblem<SpineElement<AxisymmetricQCrouzeixRaviartElement > >
  problem(n_r,n_z,l_r,h);
 
 // Run the steady simulation
  problem.steady_run();
 
} // End of main
