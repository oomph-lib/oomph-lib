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

// Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;
using namespace oomph;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//======start_of_TanhSolnForPoisson=======================================
/// Namespace for exact solution for Poisson equation
//========================================================================
namespace TanhSolnForPoisson
{
 /// Parameter for steepness of "step"
 double Alpha=4.0;

 /// Parameter for angle Phi of "step"
 double TanPhi=1.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x,Vector<double>& u)
 {
  // Assign the solution value
  u[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 } // End of get_exact_u

 /// Source function required to make the solution above an exact solution
 void source_function(const Vector<double>& x,double& source)
 {
  // Return the value of the source function at x, i.e. f(x)
  source=2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*Alpha*TanPhi*TanPhi+2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*Alpha;
 } // End of source_function

 /// Flux required by the exact solution on a boundary on which x is fixed
 void prescribed_flux_on_fixed_x_boundary(const Vector<double>& x,
                                          double& flux)
 {
  // The outer unit normal to the boundary is (1,0)
  double N[2]={1.0,0.0};

  // Calculate the value of the flux in terms of the normal
  flux=-(1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*TanPhi*N[0]+(1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*N[1];
 } // End of prescribed_flux_on_fixed_x_boundary
} // End of TanhSolnForPoisson

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//======start_of_Smoother_Factory_Function_Helper=========================
/// Returns a pointer to a smoother of the appropriate type
//========================================================================
namespace Smoother_Factory_Function_Helper
{
 /// Returns a pointer to a Gauss-Seidel Smoother object which can
 /// be used as a pre- or post-smoother
 Smoother* set_smoother()
 {
  // Create a new GS object
  return new GS<CRDoubleMatrix>;
 } // End of set_pre_smoother
} // End of Smoother_Factory_Function_Helper

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//======start_of_problem_class============================================
/// 2D Poisson problem on rectangular domain, discretised with
/// 2D QPoisson elements. Flux boundary conditions are applied
/// along boundary 1 (the boundary where x=1). The specific type of
/// element and mesh is specified via the template parameter.
//========================================================================
template<class ELEMENT,class MESH>
class FluxPoissonMGProblem : public MGProblem
{
public:

 /// Constructor: Pass a pointer to the source function
 FluxPoissonMGProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt);

 /// Destructor
 ~FluxPoissonMGProblem();

 /// Actions before adapt: Wipe the mesh of prescribed flux elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of prescribed flux elements
 void actions_after_adapt();

 /// Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// Doc the solution. DocInfo object stores flags/labels for where
 /// the output gets written to
 void doc_solution(DocInfo& doc_info);

private:

 /// Create Poisson flux elements on the b-th boundary of the
 /// problem's mesh
 void create_flux_elements(const unsigned& b);

 /// Delete Poisson flux elements and wipe the surface mesh
 void delete_flux_elements();

 /// Build and set multgrid solver
 void set_multigrid_solver();

 /// Pointer to the bulk mesh. Overloads the pure virtual function in
 /// the abstract base class, MGProblem. Must be refineable to allow the
 /// use of refine_base_mesh_as_in_reference_mesh_minus_one() in make_copy()
 TreeBasedRefineableMeshBase* mg_bulk_mesh_pt()
  {
   // Return the pointer to the bulk mesh
   return Bulk_mesh_pt;
  } // End of mg_bulk_mesh_pt

 /// Return a pointer to a new instance of the same problem.
 MGProblem* make_new_problem()
  {
   // Make new problem of the FluxPoissonMGProblem class whose template
   // parameters are specified by the template parameters of the current
   // problem
   return new FluxPoissonMGProblem<ELEMENT,MESH>
    (&TanhSolnForPoisson::source_function);
  } // End of make_new_problem

 /// Pointer to the "bulk" mesh
 MESH* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;
}; // End of FluxPoissonMGProblem class

//=======start_of_constructor=============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT,class MESH>
FluxPoissonMGProblem<ELEMENT,MESH>::
FluxPoissonMGProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt)
 : Source_fct_pt(source_fct_pt)
{
 // Set up the multigrid solver
 set_multigrid_solver();

 // Build the mesh:
 //----------------
 // Number of elements in the x-direction
 unsigned n_x=4;

 // Number of elements in the y-direction
 unsigned n_y=4;

 // Domain length in the x-direction
 double l_x=1.0;

 // Domain length in the y-direction
 double l_y=2.0;

 // Build the "bulk" mesh
 Bulk_mesh_pt=new MESH(n_x,n_y,l_x,l_y);

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
 // Create "surface mesh" that will contain only the prescribed-flux
 // elements. The constructor just creates the mesh without giving it
 // any elements, nodes, etc.
 Surface_mesh_pt=new Mesh;

 // Set boundary conditions:
 //-------------------------
 // Create prescribed-flux elements from all elements that are adjacent to
 // boundary 1, but add them to a separate mesh
 create_flux_elements(1);

 // Find out how many boundaries the "bulk" mesh possesses
 unsigned n_bound=Bulk_mesh_pt->nboundary();

 // Loop over the boundaries of the "bulk" mesh
 for(unsigned b=0;b<n_bound;b++)
 {
  // Leave nodes on boundary 1 free
  if(b!=1)
  {
   // Find the number of nodes on the b-th boundary
   unsigned n_node=Bulk_mesh_pt->nboundary_node(b);

   // Loop over the boundary nodes
   for (unsigned n=0;n<n_node;n++)
   {
    // All nodes are free by default -- we just need to pin the ones that
    // have Dirichlet conditions here
    Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
   }
  } // if(b!=1)
 } // for(unsigned b=0;b<n_bound;b++)

 // Build the global mesh from its submeshes:
 //------------------------------------------
 // Add the "bulk" mesh to the Problem's submeshes
 add_sub_mesh(Bulk_mesh_pt);

 // Add the mesh containing flux elements to the Problem's submeshes
 add_sub_mesh(Surface_mesh_pt);

 // Rebuild the Problem's global mesh from its various sub-meshes
 build_global_mesh();

 // Complete the build of all elements so they are fully functional:
 //-----------------------------------------------------------------
 // Loop over the "bulk" mesh elements and pass in the source function pointer
 for(unsigned e=0;e<Bulk_mesh_pt->nelement();e++)
 {
  // Upcast from GeneralisedElement to the element type specified by the
  // template argument
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // Set the source function pointer
  el_pt->source_fct_pt()=Source_fct_pt;
 }

 // Set up the equation numbering scheme
 assign_eqn_numbers();
} // End of FluxPoissonMGProblem constructor

//======start_of_destructor===============================================
/// Destructor for Poisson problem.
//========================================================================
template<class ELEMENT,class MESH>
FluxPoissonMGProblem<ELEMENT,MESH>::~FluxPoissonMGProblem()
{
 // Delete the MG solver pointer
 delete linear_solver_pt();

 // Set the MG solver pointer to null
 linear_solver_pt()=0;
 
 // Delete the error estimator
 delete Bulk_mesh_pt->spatial_error_estimator_pt();

 // Set the error estimator pointer to null
 Bulk_mesh_pt->spatial_error_estimator_pt()=0;
 
 // Delete the "bulk" mesh
 delete Bulk_mesh_pt;

 // Set the "bulk" mesh pointer to null
 Bulk_mesh_pt=0;

 // Delete the mesh containing flux elements
 delete Surface_mesh_pt;

 // Set the surface mesh pointer to null
 Surface_mesh_pt=0;
} // End of destructor

//======start_of_actions_before_adapt=====================================
/// Actions before adapt: Wipe the mesh of prescribed flux elements
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::actions_before_adapt()
{
 // Kill the flux elements and wipe surface mesh
 delete_flux_elements();

 // Flush the submeshes
 flush_sub_meshes();

 // Add the bulk mesh back in
 add_sub_mesh(Bulk_mesh_pt);

 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
} // End of actions_before_adapt

//======start_of_actions_after_adapt======================================
/// Actions after adapt: Rebuild the mesh of prescribed flux elements
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::actions_after_adapt()
{
 // Create prescribed-flux elements from all elements that are
 // adjacent to boundary 1 and add them to surfac mesh
 create_flux_elements(1);

 // Add the flux mesh back into the global mesh
 add_sub_mesh(Surface_mesh_pt);

 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
} // End of actions_after_adapt

//======start_of_actions_before_newton_solve==============================
/// Update the problem specs before solve: Reset boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::actions_before_newton_solve()
{
 // How many boundaries are there?
 unsigned n_bound=Bulk_mesh_pt->nboundary();

 // Loop over the boundaries
 for(unsigned i=0;i<n_bound;i++)
 {
  // Only update Dirichlet nodes
  if (i!=1)
  {
   // How many nodes are there on this boundary?
   unsigned n_node=Bulk_mesh_pt->nboundary_node(i);

   // Loop over the nodes on boundary
   for (unsigned n=0;n<n_node;n++)
   {
    // Get pointer to node
    Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(i,n);

    // Extract nodal coordinates from node:
    Vector<double> x(2);
    x[0]=nod_pt->x(0);
    x[1]=nod_pt->x(1);

    // Create storage for the solution
    Vector<double> u(1);

    // Compute the value of the exact solution at the nodal point
    TanhSolnForPoisson::get_exact_u(x,u);

    // Assign the value to the one (and only) nodal value at this node
    nod_pt->set_value(0,u[0]);
   }
  } // if (i!=1)
 } // for(unsigned i=0;i<n_bound;i++)
} // End of actions_before_newton_solve

//======start_of_create_flux_elements=====================================
/// Create Poisson Flux Elements on the boundary 1
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::create_flux_elements(const unsigned& b)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_belement=Bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_belement;e++)
 {
  // Get pointer to the bulk element that is adjacent to boundary b
  ELEMENT* bulk_elem_pt=dynamic_cast<ELEMENT*>
   (Bulk_mesh_pt->boundary_element_pt(b,e));

  // Find the index of the face of the bulk element at the boundary
  int face_index=Bulk_mesh_pt->face_index_at_boundary(b,e);

  // Build the corresponding prescribed-flux element
  PoissonFluxElement<ELEMENT>* flux_element_pt=
   new PoissonFluxElement<ELEMENT>(bulk_elem_pt,face_index);

  // Add the prescribed flux element to the mesh
  Surface_mesh_pt->add_element_pt(flux_element_pt);
 }

 // Find out how many flux elements there are on boundary 1
 unsigned n_selement=Surface_mesh_pt->nelement();

 // Loop over the flux elements to pass pointer to prescribed flux function
 for(unsigned e=0;e<n_selement;e++)
 {
  // Upcast from GeneralisedElement to Poisson flux element
  PoissonFluxElement<ELEMENT>* el_pt=
   dynamic_cast< PoissonFluxElement<ELEMENT>*>(Surface_mesh_pt->element_pt(e));

  // Set the pointer to the prescribed flux function
  el_pt->flux_fct_pt()=&TanhSolnForPoisson::prescribed_flux_on_fixed_x_boundary;
 }
} // End of create_flux_elements

//======start_of_delete_flux_elements=====================================
/// Delete Poisson flux elements and wipe the surface mesh
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::delete_flux_elements()
{
 // Find out how many surface elements there are in the surface mesh
 unsigned n_element=Surface_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
 {
  // Kill the e-th surface element
  delete Surface_mesh_pt->element_pt(e);
 }

 // Wipe the mesh
 Surface_mesh_pt->flush_element_and_node_storage();
} // End of delete_flux_elements

//======start_of_set_multigrid_solver=====================================
/// Build and set the multigrid solver
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::set_multigrid_solver()
{
 // Make an object of the MGSolver class and get the pointer to it
 MGSolver<2>* mg_solver_pt=new MGSolver<2>(this);

 // Switch solver to MG
 linear_solver_pt()=mg_solver_pt;

 // Set the pre-smoother factory function
 mg_solver_pt->set_pre_smoother_factory_function
  (Smoother_Factory_Function_Helper::set_smoother);

 // Set the post-smoother factory function
 mg_solver_pt->set_post_smoother_factory_function
  (Smoother_Factory_Function_Helper::set_smoother);

 // Create a file to record the convergence history
 mg_solver_pt->open_convergence_history_file_stream("RESLT/conv.dat");
} // End of set_multigrid_solver

//======start_of_doc_solution=============================================
/// Document the solution: doc_info contains labels/output directory, etc.
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::doc_solution(DocInfo& doc_info)
{
 // Output file stream
 ofstream some_file;

 // Create an array of characters to store the filename
 char filename[100];

 // Number of plot points
 unsigned npts=3;

 // Output solution:
 //-----------------
 // Create the filename
 sprintf(filename,"%s/soln%i.dat",
	 doc_info.directory().c_str(),
	 doc_info.number());

 // Open a file with the filename given above
 some_file.open(filename);

 // Document the computed solution
 Bulk_mesh_pt->output(some_file,npts);

 // Close the file
 some_file.close();

 // Output exact solution:
 //-----------------------
 // Create the filename
 sprintf(filename,"%s/exact_soln%i.dat",
	 doc_info.directory().c_str(),
	 doc_info.number());

 // Open a file with the filename given above
 some_file.open(filename);

 // Document the exact solution
 Bulk_mesh_pt->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u);

 // Close the file
 some_file.close();

 // Document the error:
 //--------------------
 // Note: this is only calculated over the bulk mesh as we cannot properly
 // deal with two different types of elements

 // Create storage for the L2 norm of the error
 double error=0.0;

 // Create storage for the L2 norm of the solution
 double soln=0.0;

 // Create the filename
 sprintf(filename,"%s/error%i.dat",
	 doc_info.directory().c_str(),
	 doc_info.number());

 // Open a file with the filename given above
 some_file.open(filename);

 // Compute the error and document it
 Bulk_mesh_pt->compute_error(some_file,
			     TanhSolnForPoisson::get_exact_u,
			     error,soln);

 // Close the file
 some_file.close();

 // Output the L2 norm of the error and solution and the relative error:
 //---------------------------------------------------------------------
 oomph_info << "Norm of error    : " << sqrt(error) << std::endl;
 oomph_info << "Norm of solution : " << sqrt(soln) << std::endl;
 oomph_info << "Relative error   : " << sqrt(error)/sqrt(soln) << std::endl;

 // Increment the documentation counter
 doc_info.number()++;
} // End of doc_solution

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//======start_of_main=====================================================
/// Demonstrate how to solve 2D Poisson problem with flux boundary
/// conditions.
//========================================================================
int main(int argc, char **argv)
{
 //------------------------------------
 // Sort out documentation information:
 //------------------------------------
 // Create a DocInfo object to document the solution
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;

 //--------------------
 // Set up the problem:
 //--------------------
 // Use 2D quad elements with quadratic interpolation
 typedef RefineableQPoissonElement<2,3> ELEMENT;

 // Use a rectangular quad mesh templated by the chosen element type
 typedef RefineableRectangularQuadMesh<ELEMENT> MESH;

 // Set up the problem with the chosen element and mesh type
 FluxPoissonMGProblem<ELEMENT,MESH> problem
  (&TanhSolnForPoisson::source_function);

 //-------------------
 // Solve the problem:
 //-------------------
 // Set the number of times to initially refine the problem. Must be at
 // least 1 otherwise the multigrid solver uses SuperLU to solve the problem
 unsigned n_refine=1;

 // Keep refining until the minimum refinement level is reached
 for (unsigned i=0;i<n_refine;i++)
 {
  // Indicate the problem is about to be refined
  oomph_info << "\n===================="
	     << "Initial Refinement"
	     << "====================\n"
	     << std::endl;

  // Refine the problem
  problem.refine_uniformly();
 }
 
 // Set the orientation of the "step" to 45 degrees
 TanhSolnForPoisson::TanPhi=1.0;
 
 // Initial value for the steepness of the "step"
 TanhSolnForPoisson::Alpha=1.0;
 
 // Do a couple of solutions for different forcing functions
 //---------------------------------------------------------
 // Choose the number of choices of Alpha to use
 unsigned nstep=4;

 // Loop over the different choices of Alpha
 for (unsigned istep=0;istep<nstep;istep++)
 {
  // Increase the steepness of the step:
  TanhSolnForPoisson::Alpha+=2.0;

  // Tell the user what problem we're solving
  std::cout << "\n\nSolving for TanhSolnForPoisson::Alpha="
	    << TanhSolnForPoisson::Alpha << "\n" << std::endl;
   
  // Solve the problem
  problem.newton_solve();
   
  // Document the solution
  problem.doc_solution(doc_info);
 }
} // End of main
