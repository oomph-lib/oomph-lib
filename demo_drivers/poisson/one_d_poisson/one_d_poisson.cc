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
//kruemelmonster
//Driver for a simple 1D poisson problem

// Generic oomph-lib routines
#include "generic.h"

// Include Poisson elements/equations
#include "poisson.h"

// Include the mesh
#include "meshes/one_d_mesh.h"

using namespace std;

using namespace oomph;

//==start_of_namespace================================================
/// Namespace for fish-shaped solution of 1D Poisson equation
//====================================================================
namespace FishSolnOneDPoisson
{

 /// Sign of the source function 
 /// (- gives the upper half of the fish, + the lower half)
 int Sign=-1;


 /// Exact, fish-shaped solution as a 1D vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = double(Sign)*((sin(sqrt(30.0))-1.0)*x[0]-sin(sqrt(30.0)*x[0]));
 }


 /// Source function required to make the fish shape an exact solution 
 void source_function(const Vector<double>& x, double& source)
 {
  source = double(Sign)*30.0*sin(sqrt(30.0)*x[0]);
 }

} // end of namespace







//==start_of_problem_class============================================
/// 1D Poisson problem in unit interval.
//====================================================================
template<class ELEMENT> 
class OneDPoissonProblem : public Problem
{

public:

 /// Constructor: Pass number of elements and pointer to source function
 OneDPoissonProblem(const unsigned& n_element, 
                    PoissonEquations<1>::PoissonSourceFctPt source_fct_pt);

 /// Destructor (empty)
 ~OneDPoissonProblem()
  {
   delete mesh_pt();
  }

 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve();

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution, pass the number of the case considered,
 /// so that output files can be distinguished.
 void doc_solution(const unsigned& label);

private:

 /// Pointer to source function
 PoissonEquations<1>::PoissonSourceFctPt Source_fct_pt;

}; // end of problem class





//=====start_of_constructor===============================================
/// Constructor for 1D Poisson problem in unit interval.
/// Discretise the 1D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function. 
//========================================================================
template<class ELEMENT>
OneDPoissonProblem<ELEMENT>::OneDPoissonProblem(const unsigned& n_element,
 PoissonEquations<1>::PoissonSourceFctPt source_fct_pt) : 
 Source_fct_pt(source_fct_pt)
{ 
 Problem::Sparse_assembly_method = Perform_assembly_using_two_arrays;

// Problem::Problem_is_nonlinear = false;
 // Set domain length 
 double L=1.0;

 // Build mesh and store pointer in Problem
 Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element,L);

 // Set the boundary conditions for this problem: By default, all nodal
 // values are free -- we only need to pin the ones that have 
 // Dirichlet conditions. 

 // Pin the single nodal value at the single node on mesh 
 // boundary 0 (= the left domain boundary at x=0)
 mesh_pt()->boundary_node_pt(0,0)->pin(0);
 
 // Pin the single nodal value at the single node on mesh 
 // boundary 1 (= the right domain boundary at x=1)
 mesh_pt()->boundary_node_pt(1,0)->pin(0);

 // Complete the setup of the 1D Poisson problem:

 // Loop over elements and set pointers to source function
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   //Set the source function pointer
   elem_pt->source_fct_pt() = Source_fct_pt;
  }

 // Setup equation numbering scheme
 assign_eqn_numbers();

} // end of constructor




//===start_of_actions_before_newton_solve========================================
/// Update the problem specs before solve: (Re)set boundary values
/// from the exact solution. 
//========================================================================
template<class ELEMENT>
void OneDPoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 
 // Assign boundary values for this problem by reading them out
 // from the exact solution.

 // Left boundary is node 0 in the mesh:
 Node* left_node_pt=mesh_pt()->node_pt(0);

 // Determine the position of the boundary node (the exact solution
 // requires the coordinate in a 1D vector!)
 Vector<double> x(1);
 x[0]=left_node_pt->x(0);
 
 // Boundary value (read in from exact solution which returns
 // the solution in a 1D vector)
 Vector<double> u(1);
 FishSolnOneDPoisson::get_exact_u(x,u);
 
 // Assign the boundary condition to one (and only) nodal value
 left_node_pt->set_value(0,u[0]);


 // Right boundary is last node in the mesh:
 unsigned last_node=mesh_pt()->nnode()-1;
 Node* right_node_pt=mesh_pt()->node_pt(last_node);

 // Determine the position of the boundary node
 x[0]=right_node_pt->x(0);
 
 // Boundary value (read in from exact solution which returns
 // the solution in a 1D vector)
 FishSolnOneDPoisson::get_exact_u(x,u);
 
 // Assign the boundary condition to one (and only) nodal value
 right_node_pt->set_value(0,u[0]);

 
} // end of actions before solve



//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT>
void OneDPoissonProblem<ELEMENT>::doc_solution(const unsigned& label)
{ 
 using namespace StringConversion;

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution with specified number of plot points per element
 ofstream solution_file(("soln" + to_string(label) + ".dat").c_str());
 mesh_pt()->output(solution_file,npts);
 solution_file.close();

 // Output exact solution at much higher resolution (so we can
 // see how well the solutions agree between nodal points)
 ofstream exact_file(("exact_soln" + to_string(label) + ".dat").c_str());
 mesh_pt()->output_fct(exact_file,20*npts,FishSolnOneDPoisson::get_exact_u); 
 exact_file.close();

 // Doc pointwise error and compute norm of error and of the solution
 double error,norm;
 ofstream error_file(("error" + to_string(label) + ".dat").c_str());
 mesh_pt()->compute_error(error_file,FishSolnOneDPoisson::get_exact_u,
                          error,norm); 
 error_file.close();

 // Doc error norm:
 cout << "\nNorm of error    : " << sqrt(error) << std::endl; 
 cout << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
 cout << std::endl;

} // end of doc

 

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//======start_of_main==================================================
/// Driver for 1D Poisson problem
//=====================================================================
int main()
{

 // Set up the problem: 
 // Solve a 1D Poisson problem using a source function that generates
 // a fish shaped exact solution
 unsigned n_element=40; //Number of elements
 OneDPoissonProblem<QPoissonElement<1,4> > //Element type as template parameter
  problem(n_element,FishSolnOneDPoisson::source_function);

 // Check whether the problem can be solved
 cout << "\n\n\nProblem self-test ";
 if (problem.self_test()==0)  
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("failed!",
                       OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
  }

 // Set the sign of the source function:
 cout << "\n\n\nSolving with negative sign:\n" << std::endl;
 FishSolnOneDPoisson::Sign=-1;

 // Solve the problem with this Sign
 problem.newton_solve();

 //Output solution for this case (label output files with "0")
 problem.doc_solution(0);


 // Change the sign of the source function:
 cout << "\n\n\nSolving with positive sign:\n" << std::endl;
 FishSolnOneDPoisson::Sign=1;

 // Re-solve the problem with this Sign (boundary conditions get
 // updated automatically when Problem::actions_before_newton_solve() is
 // called.
 problem.newton_solve();

 //Output solution for this case (label output files with "1")
 problem.doc_solution(1);

} // end of main









