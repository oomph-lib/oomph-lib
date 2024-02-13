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
// Driver for a 2D steady axisymetric advection diffusion problem with Robin 
// boundary condition using two separate meshes for the bulk and surface meshes.

//Generic routines
#include "generic.h"

// The Axisymetric Advection Diffusion equations
#include "steady_axisym_advection_diffusion.h"

// The mesh 
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;


//======start_of_namespace============================================
/// Namespace for exact solution for SteadyAxisymAdvectionDiffusion equation 
/// with "sharp" step 
//====================================================================
namespace TanhSolnForSteadyAxisymAdvectionDiffusion
{

 /// Peclet number
 double Peclet = 0.0;

 /// Parameter for angle of step
 double TanPhi = 0.5;

 /// Parameter for switch z-dependence
 double Omega = 0.5;

 /// Amplitude of boundary deflection (Epsilon<R!!!)
 double Epsilon = 0.1;

 /// Relation between alfa and beta (\f$\alpha = K_{\alpha}\cdot\beta $\f)
 double K_Alpha = 1.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, 
                  Vector<double>& u)
 {
  double f_rz = 1.0-(TanPhi*pow(x[0],2.0)-Omega*x[1]); 
  u[0] = tanh(f_rz);
 }

 /// Exact solution as a scalar
 void get_exact_u(const Vector<double>& x, 
                  double& u)
 {
  double f_rz = 1.0-(TanPhi*pow(x[0],2.0)-Omega*x[1]); 
  u = tanh(f_rz);

 }

 /// Source function required to make the solution above an exact solution 
 void source_function(const Vector<double>& x_vect, 
                      double& source)
 {
  double r = x_vect[0];
  double z = x_vect[1];

  double f_rz = 1.0-(TanPhi*pow(r,2.0)-Omega*z); 

   source =  
    (1.0-pow(tanh(f_rz),2.0))*
    ( Peclet*(-2*r*TanPhi*sin(6.0*z) + Omega*cos(6.0*r)) + 
      2.0*tanh(f_rz)*( pow(Omega,2.0) + 4.0*pow(r*TanPhi,2.0) ) +
      4.0*TanPhi); 
 }

/// Beta required by the exact solution on a boundary on which r is fixed
 void prescribed_beta_on_fixed_r_boundary(const Vector<double>& x_vect, 
                                          double& beta)
 {
  double r = x_vect[0];
  double z = x_vect[1];

  double f_rz = 1.0-(TanPhi*pow(r,2.0)-Omega*z);

  //The beta is (for the problem with bound1 = F(z) = R + epsilon*sin(z) with 0<z<H=2*Pi)

  beta = -1.0*( 1.0/( sqrt(1.0 + pow(Epsilon*cos(z),2.0) ) ) )*
              ( 1.0/( 1.0 -  K_Alpha*tanh(f_rz) ) )*
              ( 1.0 - pow(tanh(f_rz),2.0))*
              (-2.0*TanPhi*r + Epsilon*cos(z)*Omega);

 }

/// Alfa required by the exact solution on a boundary on which r is fixed
 void prescribed_alpha_on_fixed_r_boundary(const Vector<double>& x_vect, 
                                           double& alpha)
 {
  double r = x_vect[0];
  double z = x_vect[1];

  double f_rz = 1.0-(TanPhi*pow(r,2.0)-Omega*z);

  alpha = -K_Alpha*( 1.0/( sqrt(1.0 + pow(Epsilon*cos(z),2.0) ) ) )*
                   ( 1.0/( 1.0 -  K_Alpha*tanh(f_rz) ) )*
                   ( 1.0 - pow(tanh(f_rz),2.0))*
                   (-2.0*TanPhi*r + Epsilon*cos(z)*Omega);

 }

 /// Wind
 void wind_function(const Vector<double>& x, 
                    Vector<double>& wind)
 {
  double r = x[0];
  double z = x[1];

  wind[0] = sin(6.0*z);
  wind[1] = cos(6.0*r);

 }
 
} // end of namespace

/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////

//====== start_of_problem_class=======================================
/// 2D SteadyAxisymAdvectionDiffusion problem on rectangular domain, 
/// discretised with 2D QSteadyAxisymAdvectionDiffusion elements. The 
/// specific type of element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class TwoMeshFluxSteadyAxisymAdvectionDiffusionProblem : public Problem
{

public:

 TwoMeshFluxSteadyAxisymAdvectionDiffusionProblem(
  SteadyAxisymAdvectionDiffusionEquations::SteadyAxisymAdvectionDiffusionSourceFctPt source_fct_pt,
  SteadyAxisymAdvectionDiffusionEquations::SteadyAxisymAdvectionDiffusionWindFctPt wind_fct_pt);

 /// Destructor. Empty
 ~TwoMeshFluxSteadyAxisymAdvectionDiffusionProblem(){}

 /// Doc the solution.
 void doc_solution(DocInfo& doc_info);

private:

 /// Update the problem specs before solve: Reset boundary conditions
 /// to the values from the tanh solution.
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
 RectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

 /// Pointer to source function
 SteadyAxisymAdvectionDiffusionEquations::SteadyAxisymAdvectionDiffusionSourceFctPt Source_fct_pt;

 /// Pointer to wind function
 SteadyAxisymAdvectionDiffusionEquations::SteadyAxisymAdvectionDiffusionWindFctPt Wind_fct_pt;

}; //End of class TwoMeshFluxSteadyAxisymAdvectionDiffusionProblem



//=====start_of_constructor===============================================
/// Constructor for SteadyAxisymAdvectionDiffusion problem: Pass 
/// pointer to  source and wind functions.
//========================================================================
template<class ELEMENT>
TwoMeshFluxSteadyAxisymAdvectionDiffusionProblem<ELEMENT>::
TwoMeshFluxSteadyAxisymAdvectionDiffusionProblem(
 SteadyAxisymAdvectionDiffusionEquations::SteadyAxisymAdvectionDiffusionSourceFctPt source_fct_pt,
 SteadyAxisymAdvectionDiffusionEquations::SteadyAxisymAdvectionDiffusionWindFctPt wind_fct_pt)
       :  Source_fct_pt(source_fct_pt), Wind_fct_pt(wind_fct_pt)
{ 

 // Setup "bulk" mesh

 // # of elements in r-direction
 unsigned n_r = 20;

 // # of elements in z-direction
 unsigned n_z = 20;

 // Domain length in r-direction
 double l_r = 1.0;

 // Domain length in z-direction
 double l_z = 2.0*MathematicalConstants::Pi;

 // Build and assign "bulk" mesh

 Bulk_mesh_pt = new RectangularQuadMesh<ELEMENT>(n_r,n_z,l_r,l_z);

 // MODIFY THE MESH 
 // ---------------

 // Determine number of nodes in mesh
 const unsigned n_node = Bulk_mesh_pt->nnode();

 // Loop over nodes
 for(unsigned i=0;i<n_node;i++)
  {
     
   // Determine z and r coordinate of node
   double z_value = Bulk_mesh_pt->node_pt(i)->x(1);
   double r_value = Bulk_mesh_pt->node_pt(i)->x(0);

   // Set mesh using the external surface
   double epsilon = TanhSolnForSteadyAxisymAdvectionDiffusion::Epsilon;

   Bulk_mesh_pt->node_pt(i)->x(0) = r_value - (r_value/l_r)*epsilon*sin(z_value);

  } // End of loop over nodes
   
 // Update nodes in bulk mesh
 Bulk_mesh_pt->node_update(); 

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
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();
 
 // BOUNDARY CONDITIONS
 // ------------------- 
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here
 unsigned num_bound = Bulk_mesh_pt->nboundary();

 //Leave nodes on boundary 3 (symmetry axis) free
 for(unsigned ibound=0;ibound<num_bound-1;ibound++)
  {
   //Leave nodes on boundary 1 free
   if (ibound==0 || ibound==2)
    {
     unsigned num_nod= mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0); 
      }
    }
  } // end loop over boundaries
 
 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor: Pass pointer to source function
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;

   //Set the wind function pointer
   el_pt->wind_fct_pt() = Wind_fct_pt;

   // Set the Peclet number pointer
   el_pt->pe_pt() = &TanhSolnForSteadyAxisymAdvectionDiffusion::Peclet;
  }

 // Loop over the flux elements to pass pointer to prescribed flux function
 n_element = Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to AdvectionDiffusion flux element
   SteadyAxisymAdvectionDiffusionFluxElement<ELEMENT> *el_pt = 
    dynamic_cast< SteadyAxisymAdvectionDiffusionFluxElement<ELEMENT>*>(
     Surface_mesh_pt->element_pt(e));

   // Set the pointer to the prescribed beta function
   el_pt->beta_fct_pt() = 
    &TanhSolnForSteadyAxisymAdvectionDiffusion::prescribed_beta_on_fixed_r_boundary;
   // Set the pointer to the prescribed alpha function
   el_pt->alpha_fct_pt() = 
    &TanhSolnForSteadyAxisymAdvectionDiffusion::prescribed_alpha_on_fixed_r_boundary;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor


//===================start_of_actions_before_newton_solve================
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the tanh solution.
//========================================================================
template<class ELEMENT>
void TwoMeshFluxSteadyAxisymAdvectionDiffusionProblem<ELEMENT>::
actions_before_newton_solve()
{
 // How many boundaries are there?
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 
 //Loop over the boundaries (except the symmetry axis bound)
 for(unsigned ibound=0;ibound<num_bound-1;ibound++)
  {
   // Only update the Dirichlet conditions!
   if (ibound==0 || ibound==2)
    {
     // How many nodes are there on this boundary?
     unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);

     // Loop over the nodes on boundary
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get pointer to node
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound,inod);

       // Extract nodal coordinates from node:
       Vector<double> x(2);
       x[0] = nod_pt->x(0);
       x[1] = nod_pt->x(1);

       // Compute the value of the exact solution at the nodal point
       Vector<double> u(1);
       TanhSolnForSteadyAxisymAdvectionDiffusion::get_exact_u(x,u);

       // Assign the value to the one (and only) nodal value at this node
       nod_pt->set_value(0,u[0]);

      }
    }
  } 
}  // end of actions before solve



//===============start_of_doc=============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TwoMeshFluxSteadyAxisymAdvectionDiffusionProblem<ELEMENT>::
doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts = 5;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,
                          npts,
                          TanhSolnForSteadyAxisymAdvectionDiffusion::get_exact_u); 
 some_file.close();

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,
                             TanhSolnForSteadyAxisymAdvectionDiffusion::get_exact_u,
                             error,
                             norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

} // end of doc

//============start_of_create_flux_elements==============================
/// Create AdvectionDiffusion Flux Elements on the b-th boundary of 
/// the Mesh object pointed to by bulk_mesh_pt and add the elements 
/// to the Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TwoMeshFluxSteadyAxisymAdvectionDiffusionProblem<ELEMENT>::
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
void TwoMeshFluxSteadyAxisymAdvectionDiffusionProblem<ELEMENT>::
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


//===== start_of_main=====================================================
/// Driver code for 2D SteadyAxisymAdvectionDiffusion problem
//========================================================================
int main()
{

 //Set up the problem
 //------------------

 // Create the problem with 2D nine-node elements from the
 // QSteadyAxisymAdvectionDiffusionElement family. Pass pointer to 
 // source and wind function. 

 TwoMeshFluxSteadyAxisymAdvectionDiffusionProblem<QSteadyAxisymAdvectionDiffusionElement<3> > 
 problem(&TanhSolnForSteadyAxisymAdvectionDiffusion::source_function,
         &TanhSolnForSteadyAxisymAdvectionDiffusion::wind_function);

 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number() = 0;

 // Check if we're ready to go:
 //----------------------------
 cout << "\n\n\nProblem self-test ";
 if (problem.self_test()==0) 
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("Self test failed",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Solve the problem
 problem.newton_solve();

 //Output the solution
 problem.doc_solution(doc_info);
 
} // end of main
