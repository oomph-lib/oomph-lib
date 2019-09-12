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
//Driver for a simple 2D poisson problem with adaptive mesh refinement

//Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/simple_rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//===== start_of_namespace=============================================
/// Namespace for exact solution for Poisson equation with "sharp step" 
//=====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of "step"
 double Alpha=5.0;

 /// Parameter for angle Phi of "step"
 double TanPhi=0.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Exact gradient as a Vector
 void get_exact_gradient(const Vector<double>& x, Vector<double>& dudx)
 {
  //dudx[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
  double CoshTerm = cosh(1.0-Alpha*(TanPhi*x[0]-x[1]))
                   *cosh(1.0-Alpha*(TanPhi*x[0]-x[1]));
  dudx[0] = -Alpha*TanPhi/CoshTerm;
  dudx[1] = Alpha/CoshTerm;
 }

 /// Source function required to make the solution above an exact solution 
 void get_source(const Vector<double>& x, double& source)
 {
  source = 2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*Alpha*TanPhi*TanPhi+2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*Alpha;
 }
 
} // end of namespace



//==============================start_of_mesh======================
/// p-refineable equivalent of the SimpleRectangularQuadMesh.
/// Refinement is performed by the QuadTree-based procedures
/// implemented in the RefineableQuadMesh base class.
//=================================================================
template<class ELEMENT>
class SimpleRefineableRectangularQuadMesh : 
 public virtual SimpleRectangularQuadMesh<ELEMENT>,  
 public RefineableQuadMesh<ELEMENT>
{ 

public: 

 /// \short  Pass number of elements in the horizontal 
 /// and vertical directions, and the corresponding dimensions.
 /// Timestepper defaults to Static.
 SimpleRefineableRectangularQuadMesh(const unsigned &Nx,
                                     const unsigned &Ny, 
                                     const double &Lx, const double &Ly,
                                     TimeStepper* time_stepper_pt=
                                     &Mesh::Default_TimeStepper) :
  SimpleRectangularQuadMesh<ELEMENT>(Nx,Ny,Lx,Ly,time_stepper_pt)
  {
   // Nodal positions etc. were created in constructor for
   // SimpleRectangularQuadMesh<...> --> We only need to set up 
   // adaptivity information: Associate finite elements with their 
   // QuadTrees and plant them in a QuadTreeForest:
   this->setup_quadtree_forest();
  } // end of constructor
 

 /// Destructor: Empty
 virtual ~SimpleRefineableRectangularQuadMesh() {}

}; // end of mesh



//====== start_of_problem_class=======================================
/// hp-refineable 2D Poisson problem on rectangular domain, discretised
/// with hp-refineable 2D QPoisson elements. The specific type of
/// element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class RefineableTwoDPoissonProblem : public Problem
{
public:

 /// Constructor: Pass pointer to source function
 RefineableTwoDPoissonProblem(PoissonEquations<2>::PoissonSourceFctPt 
                              source_fct_pt);

 /// Destructor (empty)
 ~RefineableTwoDPoissonProblem()
  {
   delete mesh_pt()->spatial_error_estimator_pt();
   delete Problem::mesh_pt();
  }

 /// \short Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// \short Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

 /// \short Overloaded version of the Problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 SimpleRefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<SimpleRefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

private:

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

}; // end of problem class




//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
RefineableTwoDPoissonProblem<ELEMENT>::
      RefineableTwoDPoissonProblem(PoissonEquations<2>::PoissonSourceFctPt 
                               source_fct_pt)
       :  Source_fct_pt(source_fct_pt)
{ 

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=2;

 // # of elements in y-direction
 unsigned n_y=2;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=2.0;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new SimpleRefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
    }
  }

 // Complete the build of all elements so they are fully functional

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor: Pass pointer to source function
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor




//=================================start_of_actions_before_newton_solve===
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void RefineableTwoDPoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 // How many boundaries are there?
 unsigned num_bound = mesh_pt()->nboundary();
 
 //Loop over the boundaries
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // How many nodes are there on this boundary?
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);

   // Loop over the nodes on boundary
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
     
     // Pin the node
     nod_pt->pin(0);

     // Extract nodal coordinates from node:
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     // Compute the value of the exact solution at the nodal point
     Vector<double> u(1);
     TanhSolnForPoisson::get_exact_u(x,u);

     // Assign the value to the one (and only) nodal value at this node
     nod_pt->set_value(0,u[0]);
    }
  } 
}  // end of actions before solve



//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void RefineableTwoDPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();


 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
 some_file.close();

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

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
 bool Compute_singular_values=false;
 
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

//===== start_of_main=====================================================
/// Driver code for hp-adaptive solution to the 2D Poisson problem
//========================================================================
int main()
{

 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");
 
 //Set up the problem
 //------------------
 // Create the problem with 2D hp-refineable elements from the
 // PRefineableQPoissonElement family. Pass pointer to source function. 
 RefineableTwoDPoissonProblem<PRefineableQPoissonElement<2> > 
  problem(&TanhSolnForPoisson::get_source);

 // Set the orientation of the "step" to 45 degrees
 TanhSolnForPoisson::TanPhi=1.0;
 
 // Choose a large value for the steepness of the "step"
 TanhSolnForPoisson::Alpha=20.0;
 
 // Initially uniformly refine the mesh
 // for (unsigned p=0; p<0; p++)
 //  {
 //   cout << "p-refining:" << endl;
 //   problem.p_refine_uniformly();
 //  }
 for (unsigned h=0; h<1; h++)
 {
  cout << "h-refining:" << endl;
  problem.refine_uniformly();
 }
 
 problem.newton_solve();
 
 //Output the solution
 doc_info.number()=1;
 problem.doc_solution(doc_info);
 
 problem.p_adapt();
 problem.newton_solve();
 doc_info.number()=2;
 problem.doc_solution(doc_info);
 
 problem.adapt();
 problem.newton_solve();
 doc_info.number()=3;
 problem.doc_solution(doc_info);
 
 problem.p_adapt();
 problem.newton_solve();
 doc_info.number()=4;
 problem.doc_solution(doc_info);
 
 problem.adapt();
 problem.newton_solve();
 doc_info.number()=5;
 problem.doc_solution(doc_info);
 
 problem.p_adapt();
 problem.newton_solve();
 doc_info.number()=6;
 problem.doc_solution(doc_info);
 
 problem.adapt();
 problem.newton_solve();
 doc_info.number()=7;
 problem.doc_solution(doc_info);
 
 problem.p_adapt();
 problem.newton_solve();
 doc_info.number()=8;
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
 oomph_info.stream_pt()->precision(15);
  
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
	     << "\n" << std::endl;
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
} //end of main

