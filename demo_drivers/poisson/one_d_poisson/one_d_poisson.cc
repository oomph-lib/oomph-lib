//kruemelmonster
//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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

 /// \short Sign of the source function 
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

 /// \short Doc the solution, pass the number of the case considered,
 /// so that output files can be distinguished.
 void doc_solution(const unsigned& label);

private:

 /// Pointer to source function
 PoissonEquations<1>::PoissonSourceFctPt Source_fct_pt;

}; // end of problem class





//=====start_of_constructor===============================================
/// \short Constructor for 1D Poisson problem in unit interval.
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
/// \short Update the problem specs before solve: (Re)set boundary values
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

 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// To generate (normally distributed) random numbers
#include <iostream>
#include <chrono>
#include <random>

//========================================================================
/// Helper namespace 
//========================================================================
namespace PowerIterationHelperNamespace
{
 /// \short Construct a trivial random generator engine from a time-based
 /// seed. Used to construct an initial guess vector (which isn't orthogonal
 /// to the dominant eigenvector)
 unsigned Seed=std::chrono::system_clock::now().time_since_epoch().count();
 
 /// \short The maximum number of iterations to use inside the power (and
 /// inverse power) iteration
 unsigned Max_iteration=100000;
 
 // The tolerance for convergence
 double Tolerance=1.0e-12;

 /// \short Boolean to indicate whether we want to compute the singular
 /// values (for the condition number) or the eigenvalues
 bool Compute_singular_values=true;
 
 /// Function that returns the sign of a double (returns zero if x=0)
 int sgn(const double& x)
 {
  // If x is less than zero then we get false-true=-1 and if x is
  // greater than zero we get true-false=1. The final case is if
  // x is in fact zero. In this case we get false-false=0.
  return (x>0.0)-(x<0.0);
 } // End of sgn
 
 /// \short Function that populates the entries of the DoubleVector 
 /// random_vector with (normally-distributed) random numbers
 void randn(DoubleVector& random_vector)
 {
  // Create a trivial random generator engine from a time-based seed
  std::default_random_engine generator(Seed);

  // Set the distribution to have zero mean and unit variance 
  std::normal_distribution<double> distribution(0.0,1.0);

  // How many rows?
  unsigned n_row=random_vector.nrow();
 
  // Loop over 
  for (unsigned i=0;i<n_row;i++)
  {
   // Generate a random number and store it
   random_vector[i]=distribution(generator);
  }
 } // End of randn
} // End of PowerIterationHelperNamespace

//========================================================================
/// \short Power method: used to compute the *largest* eigenvalue (and
/// corresponding eigenvector) of the matrix pointed to by matrix_pt
/// Based on implementation described here:
///             http://www.cs.huji.ac.il/~csip/tirgul2.pdf
//========================================================================
double power_method(CRDoubleMatrix* const matrix_pt, DoubleVector& eigenvector)
{
 // Grab the LinearAlgebraDistribution associated with the input matrix
 LinearAlgebraDistribution* dist_pt=matrix_pt->distribution_pt();
 
 // If the eigenvector hasn't been built by the user
 if (eigenvector.built()==0)
 {
  // Build it!
  eigenvector.build(dist_pt,false);
 } 

 // Boolean to indicate convergence
 bool has_converged=false;
 
 // Allocate space for the initial guess
 DoubleVector q_old(dist_pt,0.0);
 
 // Allocate space for the updated guess
 DoubleVector q_new(dist_pt,0.0);
 
 // Allocate space for the matrix-vector product A*q_old
 DoubleVector z_new(dist_pt,0.0);
 
 // Temporary vector for intermediate calculations (only needed if we're
 // computing singular values instead of the eigenvalues)
 DoubleVector sing_calc_vec(dist_pt,0.0);
 
 // Populate the entries of q_old
 PowerIterationHelperNamespace::randn(q_old);
 
 // Normalise it
 q_old/=q_old.norm();

 // Storage for the eigenvalue
 double lambda=0.0;
 
 // The difference between the previous and current guess (for the eigenvector)
 double q_diff_norm=0.0;

 // How many rows in the matrix?
 unsigned n_row=matrix_pt->nrow();
 
 // How many iterations did it take to converge?
 unsigned n_iteration=0;
 
 // Iterate until convergence
 for (unsigned i=0;i<PowerIterationHelperNamespace::Max_iteration;i++)
 {  
  // If we're computing eigenvalues compute q_new=A*q_old
  if (!PowerIterationHelperNamespace::Compute_singular_values)
  {
   // Compute q_new=A*q_old 
   matrix_pt->multiply(q_old,q_new);
  }
  // If we're computing singular values we actually need to compute
  // q_new=A^{T}*A*q_old so also multiply by A^{T} 
  else
  {
   // Compute sing_calc_vec=A*q_old 
   matrix_pt->multiply(q_old,sing_calc_vec);
   
   // Compute q_new=A^{T}*A*q_old=A^{T}*sing_calc_vec
   matrix_pt->multiply_transpose(sing_calc_vec,q_new);
  }

  // Normalise
  q_new/=q_new.norm();
  
  // If we're computing eigenvalues compute q_new=A*q_new
  if (!PowerIterationHelperNamespace::Compute_singular_values)
  {
   // Compute z_new=A*q_new 
   matrix_pt->multiply(q_new,z_new);
   
   // ...and now compute the eigenvalue
   lambda=q_new.dot(z_new);
  }
  // If we're computing singular values we actually need to compute
  // z_new=A^{T}*A*q_new so also multiply by A^{T} 
  else
  {
   // Compute sing_calc_vec=A*q_new 
   matrix_pt->multiply(q_new,sing_calc_vec);
   
   // Compute z_new=A^{T}*A*q_new=A^{T}*sing_calc_vec
   matrix_pt->multiply_transpose(sing_calc_vec,z_new);
   
   // ...and now compute the singular value
   lambda=std::sqrt(q_new.dot(z_new));
  }

  // Increment the iteration counter
  n_iteration++;

  // Scale q_old and q_new so they have the same sign in the first entry as
  // the eigenvectors can be invariant up to a sign change. Scale q_old first
  q_old*=PowerIterationHelperNamespace::sgn(q_old[0]);
  
  // Now scale q_new in the same way
  q_new*=PowerIterationHelperNamespace::sgn(q_new[0]);
  
  // The norm of the difference between the previous and current guess
  q_diff_norm=0.0;

  // Loop over the entries of the vector
  for (unsigned j=0;j<n_row;j++)
  {
   // Update the value of q_diff_norm
   q_diff_norm+=std::pow(q_new[j]-q_old[j],2);
  }

  // Now square root it
  q_diff_norm=std::sqrt(q_diff_norm);
  
  // Check if the convergence tolerance has been met
  if (q_diff_norm<PowerIterationHelperNamespace::Tolerance)
  {
   // We've converged
   has_converged=true;
   
   // We're done so jump out now
   break;   
  }

  // If we're not done yet, copy the updated vector into the old one
  q_old=q_new;
 } // for (unsigned i=0;i<PowerIterationHelperNamespace::Max_iteration;i++)

 // If we've actually converged
 if (has_converged)
 {  
  // Document the convergence
  oomph_info << "\nPower method converged within " << n_iteration
	     << " iterations to accuracy: " << q_diff_norm
	     << std::endl;
 }
 // Otherwise we've run out of iterations
 else
 {
  // Document the (non-)convergence
  oomph_info << "\nPower method has not converged within " << n_iteration
	     << " iterations (with\ntolerance: "
	     << PowerIterationHelperNamespace::Tolerance
	     << "). Difference between previous and current\nguess: "
	     << q_diff_norm << std::endl;
 }

 // Copy the approximation to the dominant eigenvector over
 eigenvector=q_new;
 
 // Return the largest eigenvalue
 return lambda; 
} // End of power_method


//========================================================================
/// \short Inverse power method: used to compute the *smallest* eigenvalue
/// (and corresponding eigenvector) of the matrix pointed to by matrix_pt
/// Based on implementation described here:
///             http://www.cs.huji.ac.il/~csip/tirgul2.pdf
//========================================================================
double inverse_power_method(CRDoubleMatrix* const matrix_pt,
			    DoubleVector& eigenvector)
{
 // Grab the LinearAlgebraDistribution associated with the input matrix
 LinearAlgebraDistribution* dist_pt=matrix_pt->distribution_pt();
 
 // If the eigenvector hasn't been built by the user
 if (eigenvector.built()==0)
 {
  // Build it!
  eigenvector.build(dist_pt,false);
 } 

 // Handle to the LinearSolver associated with the input matrix
 LinearSolver* solver_pt=matrix_pt->linear_solver_pt();
   
 // Enable resolves for later
 solver_pt->enable_resolve();

 // ...and be quiet!
 solver_pt->disable_doc_time();
 
 // Allocate space for the transpose of the matrix (only used if we're
 // computing singular values instead of the eigenvalues)
 CRDoubleMatrix* matrix_transpose_pt=new CRDoubleMatrix;

 // Handle to the LinearSolver associated with the transpose of the input matrix
 LinearSolver* transpose_solver_pt=0;
 
 // If we're computing singular values we need the matrix transpose
 if (PowerIterationHelperNamespace::Compute_singular_values)
 {
  // Compute the transpose of the input matrix and store it
  matrix_pt->get_matrix_transpose(matrix_transpose_pt);
  
  // Get the LinearSolver associated with the transpose of the input matrix
  transpose_solver_pt=matrix_transpose_pt->linear_solver_pt();
 
  // Enable resolves for later (using the transposed matrix solver)
  transpose_solver_pt->enable_resolve();

  // ...and, again, be quiet!
  transpose_solver_pt->disable_doc_time();
 }

 // Make sure everything has been set up to do resolves
 {
  // Dummy RHS vector 
  DoubleVector rhs_temp(dist_pt,0.0);
 
  // Dummy LHS vector
  DoubleVector lhs_temp(dist_pt,0.0);
  
  // Solve the system once to set everything necessary for a resolve up
  matrix_pt->solve(rhs_temp,lhs_temp);
  
  // If we're computing singular values we need the matrix transpose
  if (PowerIterationHelperNamespace::Compute_singular_values)
  {
   // Solve the system once to set everything necessary for a resolve up
   matrix_transpose_pt->solve(rhs_temp,lhs_temp);
  }
 } // Set up everything necessary for a resolve
 
 // Boolean to indicate convergence
 bool has_converged=false;
 
 // Allocate space for the initial guess
 DoubleVector q_old(dist_pt,0.0);
 
 // Allocate space for the updated guess
 DoubleVector q_new(dist_pt,0.0);
 
 // Allocate space for the matrix-vector product A^{-1}*q_old
 DoubleVector z_new(dist_pt,0.0);
 
 // Temporary vector for intermediate calculations (only needed if we're
 // computing singular values instead of the eigenvalues)
 DoubleVector sing_calc_vec(dist_pt,0.0);
 
 // Populate the entries of q_old
 PowerIterationHelperNamespace::randn(q_old);
 
 // Normalise it
 q_old/=q_old.norm();

 // Storage for the eigenvalue
 double lambda=0.0;
 
 // The difference between the previous and current guess (for the eigenvector)
 double q_diff_norm=0.0;

 // How many rows in the matrix?
 unsigned n_row=matrix_pt->nrow();

 // How many iterations did it take to converge?
 unsigned n_iteration=0;
 
 // Iterate until convergence
 for (unsigned i=0;i<PowerIterationHelperNamespace::Max_iteration;i++)
 {  
  // Reinitialise
  q_new.initialise(0.0);
  
  // If we're computing eigenvalues compute q_new=A*q_old
  if (!PowerIterationHelperNamespace::Compute_singular_values)
  {
   // Compute q_new=A^{-1}*q_old 
   solver_pt->resolve(q_old,q_new);
  }
  // If we're computing singular values we actually need to compute
  // q_new=A^{T}*A*q_old so also multiply by A^{T} 
  else
  {
   // Compute q_new=A^{-1}*q_old 
   solver_pt->resolve(q_old,sing_calc_vec);
   
   // Compute q_new=A^{-T}*A^{-1}*q_old=A^{-T}*sing_calc_vec
   transpose_solver_pt->resolve(sing_calc_vec,q_new);
  }
  
  // Normalise
  q_new/=q_new.norm();
  
  // Reinitialise
  z_new.initialise(0.0);
  
  // If we're computing eigenvalues compute q_new=A*q_old
  if (!PowerIterationHelperNamespace::Compute_singular_values)
  {
   // Compute z_new=A^{-1}*q_new
   solver_pt->resolve(q_new,z_new);
   
   // ...and now compute the singular value (recalling that we're actually
   // calculating 1/lambda with the inverse power method)
   lambda=1.0/q_new.dot(z_new);
  }
  // If we're computing singular values we actually need to compute
  // q_new=A^{-T}*A^{-1}*q_old so also multiply by A^{-T} 
  else
  {
   // Compute z_new=A^{-1}*q_new
   solver_pt->resolve(q_new,sing_calc_vec);
   
   // Compute z_new=A^{-T}*A^{-1}*q_old=A^{-T}*sing_calc_vec
   transpose_solver_pt->resolve(sing_calc_vec,z_new);
  
   // ...and now compute the singular value (recalling that we're actually
   // calculating 1/lambda with the inverse power method)
   lambda=std::sqrt(1.0/q_new.dot(z_new));
  }
  
  // Increment the iteration counter
  n_iteration++;

  // Scale q_old and q_new so they have the same sign in the first entry as
  // the eigenvectors can be invariant up to a sign change. Scale q_old first
  q_old*=PowerIterationHelperNamespace::sgn(q_old[0]);
  
  // Now scale q_new in the same way
  q_new*=PowerIterationHelperNamespace::sgn(q_new[0]);
  
  // The norm of the difference between the previous and current guess
  q_diff_norm=0.0;

  // Loop over the entries of the vector
  for (unsigned j=0;j<n_row;j++)
  {
   // Update the value of q_diff_norm
   q_diff_norm+=std::pow(q_new[j]-q_old[j],2);
  }

  // Now square root it
  q_diff_norm=std::sqrt(q_diff_norm);

  // Check if the convergence tolerance has been met
  if (q_diff_norm<PowerIterationHelperNamespace::Tolerance)
  {
   // We've converged
   has_converged=true;
   
   // We're done so jump out now
   break;   
  }

  // If we're not done yet, copy the updated vector into the old one
  q_old=q_new;
 } // for (unsigned i=0;i<PowerIterationHelperNamespace::Max_iteration;i++)

 // If we've actually converged
 if (has_converged)
 {
  // Document the convergence
  oomph_info << "\nInverse power method converged within " << n_iteration
	     << " iterations to accuracy: " << q_diff_norm
	     << std::endl;
 }
 // Otherwise we've run out of iterations
 else
 {
  // Document the (non-)convergence
  oomph_info << "\nInverse power method has not converged within "
	     << n_iteration << " iterations (with\ntolerance: "
	     << PowerIterationHelperNamespace::Tolerance
	     << "). Difference between previous and current\nguess: "
	     << q_diff_norm << std::endl;
 }

 // Copy the approximation to the dominant eigenvector over
 eigenvector=q_new;
 
 // Return the smallest eigenvalue
 return lambda; 
} // End of inverse_power_method

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//======start_of_main==================================================
/// Driver for 1D Poisson problem
//=====================================================================
int main()
{
 // The number of elements in the 1D mesh
 unsigned n_element=40; 
 
 // Solve a 1D Poisson Problem (using a source function that generates
 // a fish shaped exact solution)
 OneDPoissonProblem<QPoissonElement<1,2> > problem(
  n_element,FishSolnOneDPoisson::source_function);

 // Set the sign of the source function
 FishSolnOneDPoisson::Sign=-1;

 // Storage for the Jacobian
 CRDoubleMatrix matrix;

 // Storage for the (to be unused) residual vector
 DoubleVector residual;

 // Storage for the eigenvector associated with the largest eigenvalue
 DoubleVector eigenvector_max;

 // Storage for the eigenvector associated with the smallest eigenvalue
 DoubleVector eigenvector_min;
 
 // Get the Jacobian matrix and residual vector
 problem.get_jacobian(residual,matrix);
 
 // Tell the user
 oomph_info << "\n=================================================="
	    << "\nBeginning eigenvalue/singular value calculation..."
	    << "\n=================================================="
	    << std::endl;
 
 // Storage for the largest eigenvalue/singular value
 double lambda_max=power_method(&matrix,eigenvector_max);

 // Storage for the smallest eigenvalue/singular value
 double lambda_min=inverse_power_method(&matrix,eigenvector_min);

 // Set the output precision
 oomph_info.stream_pt()->precision(10);
  
 // If we're just computing eigenvalues
 if (!PowerIterationHelperNamespace::Compute_singular_values)
 {
  // Output the eigenvalues
  oomph_info << "\nLargest eigenvalue: " << lambda_max
	     << "\nSmallest eigenvalue: " << lambda_min
	     << std::endl;
 }
 // If we're computing singular values instead
 else
 {
  // Output the singular values
  oomph_info << "\nLargest singular value: " << lambda_max
	     << "\nSmallest singular value: " << lambda_min
	     << "\nCondition number: " << lambda_max/lambda_min
	     << std::endl;
 }

 // Create a method of outputting the Jacobian
 std::ofstream matrix_file("jacobian.dat");

 // Dump the matrix so we can check that we're doing the right thing...
 matrix.sparse_indexed_output(matrix_file);

 // Close the file
 matrix_file.close();

 // Output the eigenvector associated with the largest eigenvalue
 eigenvector_max.output("eigenvector_max.dat");
 
 // Also output the eigenvector associated with the smallest eigenvalue
 eigenvector_min.output("eigenvector_min.dat"); 
} // End of main









