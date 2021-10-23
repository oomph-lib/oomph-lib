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
//Driver for a 3D poisson problem

//Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/simple_cubic_mesh.h"

using namespace std;

using namespace oomph;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//======start_of_namespace================================================
/// Returns a pointer to a smoother of the appropriate type
//========================================================================
namespace Smoother_Factory_Function_Helper
{
 /// Returns a pointer to a Gauss-Seidel Smoother object which is
 /// to be used as a pre- or post-smoother
 Smoother* set_smoother()
 {
  // Create a new GS object
  return new GS<CRDoubleMatrix>;
 } // End of set_pre_smoother
} // End of Smoother_Factory_Function_Helper

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//=============start_of_namespace=====================================
/// Namespace for exact solution for Poisson equation with sharp step 
//====================================================================
namespace TanhSolnForPoisson
{
 /// Parameter for steepness of step
 double Alpha=1.0;

 /// Orientation (non-normalised x-component of unit vector in
 /// direction of step plane)
 double N_x=-1.0;

 /// Orientation (non-normalised y-component of unit vector in
 /// direction of step plane)
 double N_y=-1.0;

 /// Orientation (non-normalised z-component of unit vector in
 /// direction of step plane)
 double N_z=1.0;

 /// Orientation (x-coordinate of step plane) 
 double X_0=0.0;

 /// Orientation (y-coordinate of step plane) 
 double Y_0=0.0;

 /// Orientation (z-coordinate of step plane) 
 double Z_0=0.0;

 // Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  // Assign the solution value  
  u[0]=tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-Y_0)*
		   N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*
		   N_z/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)));
 }

 /// Exact solution as a scalar
 void get_exact_u(const Vector<double>& x, double& u)
 {
  // Assign the solution value
  u=tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-Y_0)*
		  N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*
		  N_z/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)));
 }


 /// Source function to make it an exact solution 
 void get_source(const Vector<double>& x, double& source)
 {
  // Create auxiliary variables
  double s1;
  double s2;
  double s3;

  // Assign the value of s1
  s1 = -2.0*tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z)))*(1.0-pow(tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z))),2.0))*Alpha*Alpha*N_x*N_x/(N_x*N_x+N_y*N_y+N_z*N_z);
  
  // Assign the value of s2
      s2 = -2.0*tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z)))*(1.0-pow(tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z))),2.0))*Alpha*Alpha*N_y*N_y/(N_x*N_x+N_y*N_y+N_z*N_z);
      
  // Assign the value of s3
      s3 = -2.0*tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z)))*(1.0-pow(tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z))),2.0))*Alpha*Alpha*N_z*N_z/(N_x*N_x+N_y*N_y+N_z*N_z);
      
      // Calculate the source value
      source=s1+s2+s3;
 }
} // End of TanhSolnForPoisson

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//======start_of_problem_class========================================
/// 3D Poisson problem in a unit cube with Dirichlet conditions
/// applied on all six faces. The specific type of element and mesh
/// used here is specified via the template parameter.
//====================================================================
template<class ELEMENT,class MESH>
class UnitCubePoissonMGProblem : public MGProblem
{
public:

 /// Constructor: Pass pointer to source function
 UnitCubePoissonMGProblem(PoissonEquations<3>::PoissonSourceFctPt source_fct_pt);

 /// Destructor: Empty
 ~UnitCubePoissonMGProblem()
  {
   // Delete the MG solver pointers
   delete linear_solver_pt();

   // Set the pointer to null
   linear_solver_pt()=0;
   
   // Delete the mesh that was created in the constructor
   delete mesh_pt();

   // Make the associated pointer a null pointer
   Problem::mesh_pt()=0;
  } // End of ~UnitCubePoissonMGProblem

 /// Update the problem specs before solve: 
 /// Set Dirichlet boundary conditions from exact solution.
 void actions_before_newton_solve();
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}
 
 /// Document the solution
 void doc_solution(DocInfo& doc_info);
 
private:
    
 /// Overload generic access function to the mesh by one that returns
 /// a pointer to the mesh with return type specified by the template parameter 
 MESH* mesh_pt() 
  {
   // Upcast to the appropriate mesh type
   return dynamic_cast<MESH*>(Problem::mesh_pt());
  } // End of mesh_pt
 
 /// Build and set multgrid solver
 void set_multigrid_solver();
 
 /// Pointer to the bulk mesh. Overloads the pure virtual function in
 /// the abstract base class, MGProblem. Must be refineable to allow the
 /// use of refine_base_mesh_as_in_reference_mesh_minus_one() in make_copy()
 TreeBasedRefineableMeshBase* mg_bulk_mesh_pt() 
  {
   // Return a handle to the global mesh
   return mesh_pt();
  } // End of mg_bulk_mesh_pt

 /// Return a pointer to a new instance of the same problem.
 MGProblem* make_new_problem()
  {
   // Return a pointer to a new object of the UnitCubePoissonMGProblem class
   return new UnitCubePoissonMGProblem<ELEMENT,MESH>
    (&TanhSolnForPoisson::get_source);
  } // End of make_new_problem
  
 /// Pointer to source function
 PoissonEquations<3>::PoissonSourceFctPt Source_fct_pt;
}; // End of UnitCubePoissonMGProblem class

//======start_of_constructor==============================================
/// Constructor for Poisson problem in a unit cube
//========================================================================
template<class ELEMENT,class MESH>
UnitCubePoissonMGProblem<ELEMENT,MESH>::UnitCubePoissonMGProblem
(PoissonEquations<3>::PoissonSourceFctPt source_fct_pt)
 : Source_fct_pt(source_fct_pt)
{ 
 // Set up the multigrid solver
 set_multigrid_solver();
 
 // Set the steepness of the step in the exact tanh solution
 TanhSolnForPoisson::Alpha=50.0;

 // Build the mesh:
 //---------------- 
 // The number of elements in the x-direction
 unsigned n_x=2;
 
 // The number of elements in the y-direction
 unsigned n_y=2;
 
 // The number of elements in the z-direction
 unsigned n_z=2;
 
 // The length of the mesh in the x-direction
 unsigned l_x=1.0;
 
 // The length of the mesh in the y-direction
 unsigned l_y=1.0;
 
 // The length of the mesh in the z-direction
 unsigned l_z=1.0;
  
 // Create the mesh
 Problem::mesh_pt()=new MESH(n_x,n_y,n_z,l_x,l_y,l_z);

 // Set boundary conditions:
 //-------------------------
 // Find the number of boundaries the mesh possesses
 unsigned n_bound=mesh_pt()->nboundary();

 // Loop over the boundaries of the mesh
 for(unsigned b=0;b<n_bound;b++)
 {
  // Find the number of nodes on the b-th boundary
  unsigned n_node=mesh_pt()->nboundary_node(b);

  // Loop over the boundary nodes
  for (unsigned n=0;n<n_node;n++)
  {
   // All nodes are free by default -- we just need to pin the ones that
   // have Dirichlet conditions here
   mesh_pt()->boundary_node_pt(b,n)->pin(0); 
  }
 } // for(unsigned b=0;b<n_bound;b++)

 // Find the number of elements in the mesh
 unsigned n_element=mesh_pt()->nelement();

 // Complete the build of all elements so they are fully functional:
 //-----------------------------------------------------------------
 // Loop over the elements in the mesh and pass in the source function pointer
 for(unsigned i=0;i<n_element;i++)
 {
  // Upcast the element to the type specified by the template parameter
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

  // Set the source function pointer
  el_pt->source_fct_pt()=Source_fct_pt;
 }

 // Set up the equation numbering scheme
 assign_eqn_numbers();
} // End of UnitCubePoissonMGProblem constructor

//======start_of_actions_before_newton_solve==============================
/// Update the problem specs before solve: Reset boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT,class MESH>
void UnitCubePoissonMGProblem<ELEMENT,MESH>::actions_before_newton_solve()
{
 // Get the number of boundaries in the mesh
 unsigned num_bound=mesh_pt()->nboundary();
   
 // Loop over the boundaries
 for(unsigned ibound=0;ibound<num_bound;ibound++)
 {
  // Get the number of nodes on the ibound-th boundary
  unsigned num_nod=mesh_pt()->nboundary_node(ibound);

  // Loop over the nodes on ibound-th boundary
  for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Grab the inod-th node on the ibound-th boundary
   Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

   // Create a vector to hold the (Eulerian) position of the node
   Vector<double> x(3,0.0);

   // Loop over each coordinate
   for (unsigned k=0;k<3;k++)
   {
    // Store the k-th coordinate of the node in x
    x[k]=nod_pt->x(k);
   }
     
   // Create a variable to store the solution at this position
   double u=0.0;

   // Calculate the exact solution at this point
   TanhSolnForPoisson::get_exact_u(x,u);

   // Assign the value of the solution at this node
   nod_pt->set_value(0,u);
  }
 } // for(unsigned ibound=0;ibound<num_bound;ibound++)
} // End of actions_before_newton_solve 

//======set_multigrid_solver==============================================
/// Build and set the multigrid solver
//========================================================================
template<class ELEMENT,class MESH>
void UnitCubePoissonMGProblem<ELEMENT,MESH>::set_multigrid_solver()
{
 // Make an object of the MGSolver class and get the pointer to it
 MGSolver<3>* mg_solver_pt=new MGSolver<3>(this);
 
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

//========================start_of_doc====================================
/// Doc the solution
//========================================================================
template<class ELEMENT,class MESH>
void UnitCubePoissonMGProblem<ELEMENT,MESH>::doc_solution(DocInfo& doc_info)
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
 mesh_pt()->output(some_file,npts);
 
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
 mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
 
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
 mesh_pt()->compute_error(some_file,TanhSolnForPoisson::get_exact_u,error,soln);

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

//=========start_of_main=============================================
/// Driver for 3D Poisson problem in a unit cube. Solution has
/// a sharp step.
//===================================================================
int main(int argc, char *argv[])
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
 // Set up the problem
 //--------------------
 // Using linear interpolation
 typedef RefineableQPoissonElement<3,3> ELEMENT;
  
 // Typedef the mesh and template it by the chosen element type
 typedef RefineableSimpleCubicMesh<ELEMENT> MESH;

 // Set the problem
 UnitCubePoissonMGProblem<ELEMENT,MESH> problem(&TanhSolnForPoisson::get_source);
  
 //-------------------
 // Solve the problem:
 //-------------------
 // Set the number of times to initially refine the problem. Must be at
 // least 1 otherwise the multigrid solver uses SuperLU to solve the problem
 unsigned n_refine=2;
 
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
 
 // Solve the problem
 problem.newton_solve();

 // Document the solution
 problem.doc_solution(doc_info);
} // End of main
