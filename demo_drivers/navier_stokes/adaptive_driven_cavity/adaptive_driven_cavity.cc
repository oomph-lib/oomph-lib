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
// Driver for adaptive 2D rectangular driven cavity. Solved with black
// box adaptation, using Taylor Hood and Crouzeix Raviart elements.

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=100;
} // end_of_namespace



//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain, templated
/// by element type. 
//====================================================================
template<class ELEMENT>
class RefineableDrivenCavityProblem : public Problem
{

public:

 /// Constructor
 RefineableDrivenCavityProblem();

 /// Destructor: Empty
 ~RefineableDrivenCavityProblem() {}

 /// Update the after solve (empty)
 void actions_after_newton_solve() {}

 /// \short Update the problem specs before solve. 
 /// (Re-)set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
  { 
  // Setup tangential flow along boundary 0:
  unsigned ibound=0; 
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Tangential flow
    unsigned i=0;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,1.0);
    // No penetration
    i=1;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
   }
  
  // Overwrite with no flow along all other boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=1;ibound<num_bound;ibound++)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      for (unsigned i=0;i<2;i++)
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
       }
     }
   }
  } // end_of_actions_before_newton_solve


 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
    
    // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now set the first pressure dof in the first element to 0.0
   fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
 
private:

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to proper element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

}; // end_of_problem_class



//==start_of_constructor==================================================
/// Constructor for RefineableDrivenCavity problem 
///
//========================================================================
template<class ELEMENT>
RefineableDrivenCavityProblem<ELEMENT>::RefineableDrivenCavityProblem()
{ 

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=10;

 // # of elements in y-direction
 unsigned n_y=10;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(mesh_pt())->
  spatial_error_estimator_pt()=error_estimator_pt;
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here: All boundaries are Dirichlet boundaries.
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to Reynolds
 // number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  } // end loop over elements
 
 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
  // Now set the first pressure dof in the first element to 0.0
 fix_pressure(0,0,0.0);
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end_of_constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 


 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
} // end_of_doc_solution



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


//==start_of_main======================================================
/// Driver for RefineableDrivenCavity test problem 
//=====================================================================
int main()
{
 // Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");

 // Set max. number of black-box adaptation
 unsigned max_adapt=3;

 // Build problem
 RefineableDrivenCavityProblem<RefineableQCrouzeixRaviartElement<2> > problem;
  
 // Solve the problem with automatic adaptation
 problem.newton_solve(max_adapt);
  
 // Step number
 doc_info.number()=1;

 //Output solution
 problem.doc_solution(doc_info);

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













