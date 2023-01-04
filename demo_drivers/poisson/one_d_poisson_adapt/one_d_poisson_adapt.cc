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
// Driver for a simple 1D poisson problem with adaptive mesh refinement

// Generic oomph-lib routines
#include "generic.h"

// Poisson elements/equations
#include "poisson.h"

// The mesh
#include "meshes/one_d_mesh.h"

using namespace std;

using namespace oomph;




//======start_of_namespace================================================
/// Namespace for exact solution for Poisson equation with "sharp step" 
//========================================================================
namespace ArcTanSolnForPoisson
{
 
 /// Parameter for steepness of "step"
 double Alpha=100.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = atan(Alpha*(x[0]-0.5));
 }

 /// Source function required to make the solution above an exact solution 
 void get_source(const Vector<double>& x, double& source)
 {
  // Cache tan( atan(Alpha*(x[0]-0.5)) ) term
  double tan_term = tan( atan(Alpha*(x[0]-0.5)) );
  
  // Compute source function
  source = -(2.0 * Alpha*Alpha * tan_term) /
   ( (1.0 + tan_term*tan_term)*(1.0 + tan_term*tan_term) );
 }
 
} // End of namespace




//======start_of_problem_class============================================
/// 1D Poisson problem discretised with refineable 1D QPoisson elements.
/// The specific type of element is specified via the template parameter.
//========================================================================
template<class ELEMENT> 
class RefineableOneDPoissonProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source function
 RefineableOneDPoissonProblem(PoissonEquations<1>::PoissonSourceFctPt 
                              source_fct_pt);
 
 /// Destructor (empty)
 ~RefineableOneDPoissonProblem() {}
 
 /// Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve() {}

 /// Doc the solution. DocInfo object stores flags/labels for where
 /// the output gets written to.
 void doc_solution(DocInfo& doc_info);

 /// Overloaded version of the Problem's access function to the mesh.
 /// Recasts the pointer to the base Mesh object to the actual mesh type.
 RefineableOneDMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableOneDMesh<ELEMENT>*>(Problem::mesh_pt());
  }

private:

 /// Pointer to source function
 PoissonEquations<1>::PoissonSourceFctPt Source_fct_pt;
 
}; // End of problem class




//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
RefineableOneDPoissonProblem<ELEMENT>::
RefineableOneDPoissonProblem(PoissonEquations<1>::PoissonSourceFctPt 
                             source_fct_pt)
 : Source_fct_pt(source_fct_pt)
{ 
 
 // Set up mesh
 // -----------

 // Number of elements
 const unsigned n = 2;

 // Domain length
 const double length = 1.0;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableOneDMesh<ELEMENT>(n,length);

 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt() = new Z2ErrorEstimator;
  
 // Set the boundary conditions for this problem. All nodes are free by
 // default so only need to pin the ones that have Dirichlet conditions here.
 const unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
    }
  }

 // Complete build of all elements so they are fully functional

 // Loop over the elements to set up element-specific things that cannot
 // be handled by the (argument-free!) ELEMENT constructor: Pass pointer
 // to source function
 const unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;
  }

 // Set up equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

} // End of constructor




//=====start_of_actions_before_newton_solve===============================
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void RefineableOneDPoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 // Determine the number of mesh boundaries
 const unsigned n_boundary = mesh_pt()->nboundary();
 
 // Loop over these boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine the number of nodes on this boundary b
   const unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
   
   // Loop over these nodes
   for (unsigned n=0;n<n_boundary_node;n++)
    {
     // Get pointer to node
     Node* nod_pt = mesh_pt()->boundary_node_pt(b,n);
     
     // Extract nodal coordinates from node:
     Vector<double> x(1);
     x[0] = nod_pt->x(0);

     // Compute the value of the exact solution at the nodal point
     Vector<double> u(1);
     ArcTanSolnForPoisson::get_exact_u(x,u);

     // Assign the value to the one (and only) nodal value at this node
     nod_pt->set_value(0,u[0]);
    }
  } 
}  // End of actions_before_newton_solve




//=====start_of_doc_solution==============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void RefineableOneDPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 
 // Declare output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 const unsigned npts = 5;

 // Output solution 
 // ---------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
 // Output exact solution 
 // ---------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,ArcTanSolnForPoisson::get_exact_u); 
 some_file.close();

 // Doc error and return of the square of the L2 error
 // --------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,ArcTanSolnForPoisson::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

} // End of doc_solution




//======start_of_main=====================================================
/// Driver code for 1D Poisson problem
//========================================================================
int main()
{
 
 //Set up the problem
 //------------------
 
 // Create the problem with 1D three-node refineable elements from the
 // RefineableLinePoissonElement family. Pass pointer to source function. 
 RefineableOneDPoissonProblem<RefineableQPoissonElement<1,3> > 
  problem(&ArcTanSolnForPoisson::get_source);
 
 // Create label for output
 // -----------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;
 
 // Check if we're ready to go:
 // ---------------------------
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

 // Choose a large value for the steepness of the "step"
 ArcTanSolnForPoisson::Alpha=100.0; 

 // Refine problem uniformly 5 times (to check automatic unrefinement)
 for(unsigned i=0;i<4;i++) { problem.refine_uniformly(); }

 // Solve the problem, performing up to 10 adaptive refinements
 problem.newton_solve(10);
 
 // Output the solution
 problem.doc_solution(doc_info);
 
} // End of main









