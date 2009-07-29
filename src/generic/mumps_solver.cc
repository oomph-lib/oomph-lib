//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
//Interface to MUMPS solver (fortran)

#include<iostream>
#include<vector>


//oomph-lib headers
#include "cfortran.h"
#include "mumps.h"
#include "mumps_solver.h"
#include "Vector.h"
#include "oomph_utilities.h"
#include "problem.h"

using namespace std;

namespace oomph
{

 

//====================================================================
/// \short Namespace for pool of fortran mumps solvers
//====================================================================
namespace MumpsSolverPool 
{

 ///Stack containing the IDs of available mumps solvers
 std::stack<int> Available_solver_ids;

 /// Bool indicating that the pool has been set up
 bool Pool_has_been_setup=false;
 
 /// Default max. number of mumps solvers
 int Max_n_solvers=10;
 
 /// Get new solver from pool and return its id
 void get_new_solver_id(unsigned& solver_id)
 {
  if (Available_solver_ids.size()==0)
   {
    std::ostringstream error_message_stream;
    error_message_stream 
     << "Sorry, all available mumps solvers are in use.\n" 
     << "Please increase MumpsSolverPool::Max_n_solvers from its \n" 
     << "current value of " << Max_n_solvers << " and try again." 
     << std::endl;
    
    throw OomphLibError(error_message_stream.str(),
                        "MumpsSolverPool::get_new_solver_id()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  else
   {
    // Get solver id from top of stack and pops it
    solver_id=Available_solver_ids.top();
    Available_solver_ids.pop();
   }
 }
 
 /// Return solver to pool
 void return_solver(const unsigned& solver_id)
 {
  // Put it back on stack
  Available_solver_ids.push(solver_id);
 }

 /// \short Setup namespace -- specify the max. number of solver
 /// instantiations required.
 void setup(const unsigned& max_n_solvers=Max_n_solvers)
 {
  if (!Pool_has_been_setup)
   {
    Pool_has_been_setup=true;

    // Setup pool in fortran
    mumps_setup_solver_pool(int(max_n_solvers));
    Max_n_solvers=max_n_solvers;

    // Declare available IDs -- Fortran 1-based indexing!
    for (unsigned i=1;i<=max_n_solvers;i++)
     {
      Available_solver_ids.push(i);
     }
   }
  else
   {
    std::ostringstream error_message_stream;
    error_message_stream 
     << "Mumps solver pool has already been set up. Can't do it again.\n";
    throw OomphLibError(error_message_stream.str(),
                        "MumpsSolverPool::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }

 }

}




//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////





 //=========================================================================
 /// \short Default factor for workspace -- static so it can be overwritten
 /// globally.
 //=========================================================================
 int MumpsSolver::Default_workspace_scaling_factor=1;



//=============================================================================
/// Constructor: Call setup
//=============================================================================
 MumpsSolver::MumpsSolver()
 {
  // Setup pool
  if (!MumpsSolverPool::Pool_has_been_setup)
   {
    // Setup with default number of solver instantiations
    MumpsSolverPool::setup();
   }

  // Get new solver ID from pool
  MumpsSolverPool::get_new_solver_id(Solver_ID_in_pool);

  oomph_info << "New solver id: " << Solver_ID_in_pool << std::endl;

  Doc_stats=false;
  Suppress_solve=false;
  Delete_matrix_data=false;
  mumps_setup(Solver_ID_in_pool,Default_workspace_scaling_factor);
 }
 
//=============================================================================
/// Destructor: Shutdown mumps
//=============================================================================
MumpsSolver::~MumpsSolver()
 {

  oomph_info << "Returning solver id: " << Solver_ID_in_pool << std::endl;

  /// Return solver to pool
  MumpsSolverPool::return_solver(Solver_ID_in_pool);
  
  /// Shut it down
  mumps_shutdown(Solver_ID_in_pool);
 }
 

//=============================================================================
/// Set scaling factor for workspace (defaults is 2)
//=============================================================================
void MumpsSolver::set_workspace_scaling_factor(const unsigned& s)
{
 mumps_set_workspace_scaling_factor(Solver_ID_in_pool,int(s));
}
 

//=============================================================================
/// LU decompose the matrix addressed by matrix_pt using
/// mumps. The resulting matrix factors are stored internally.
/// Note: if Delete_matrix_data is true the function 
/// matrix_pt->clean_up_memory() will be used to wipe the matrix data.
//=============================================================================
void MumpsSolver::factorise(DoubleMatrixBase* const &matrix_pt)
{

 //Check that we have a square matrix
#ifdef PARANOID
 int n = matrix_pt->nrow();
 int m = matrix_pt->ncol();
 if(n != m)
  {
   std::ostringstream error_message_stream;
   error_message_stream << "Can only solve for square matrices\n" 
                        << "N, M " << n << " " << m << std::endl;
   
   throw OomphLibError(error_message_stream.str(),
                       "MumpsSolver::factorise()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // Doc stats?
 if (Doc_stats) mumps_switch_on_doc(Solver_ID_in_pool);

 // number of processors
 unsigned nproc = MPI_Helpers::Nproc;
 if(dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt) != 0)
  {
   nproc = dynamic_cast<DistributableLinearAlgebraObject*>
    (matrix_pt)->distribution_pt()->communicator_pt()->nproc();
  }

 // Make sure any existing factors are deleted
 clean_up_memory();
   
 // Is it a CRDoubleMatrix?
 if(dynamic_cast<CRDoubleMatrix*>(matrix_pt) != 0)
  {
   // Get a cast pointer to the matrix
   CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);

#ifdef PARANOID
   // paranoid check that the matrix has been set up
   if (!cr_matrix_pt->built())
    {
     throw OomphLibError
      ("To apply MumpsSolver to a CRDoubleMatrix - it must be built",
       "MumpsSolver::factorise()",OOMPH_EXCEPTION_LOCATION);
    }
#endif
   
   // if the matrix is distributed then set up solver
   if ((cr_matrix_pt->distribution_pt()->communicator_pt()->nproc()==1)||
       (cr_matrix_pt->distributed()))
    {
     // Find the number of rows and non-zero entries in the matrix
     const int nnz_loc = int(cr_matrix_pt->nnz());
     const int n = matrix_pt->nrow();
     
     // Create mumps storage
     
     // Create vector for row numbers
     Vector<int> irn_loc(nnz_loc);
     
     // Create vector for column numbers
     Vector<int> jcn_loc(nnz_loc);
     
     // Vector for entries
     Vector<double> a_loc(nnz_loc);
     
     // First row held on this processor
     int first_row = cr_matrix_pt->first_row();

     // Copy into coordinate storage scheme using pointer arithmetic
     double* matrix_value_pt = cr_matrix_pt->value();
     int* matrix_index_pt = cr_matrix_pt->column_index();
     int* matrix_start_pt = cr_matrix_pt->row_start();
     int i_row=0;
     for(int count=0;count<nnz_loc;count++) 
      {
       a_loc[count] = matrix_value_pt[count];
       jcn_loc[count] = matrix_index_pt[count]+1;
       if (count<matrix_start_pt[i_row+1])
        {
         irn_loc[count] = first_row+i_row+1;
        }
       else
        {
         i_row++;
         irn_loc[count] = first_row+i_row+1;
        }
      }
      
     // Now delete the matrix if we are allowed
     if (Delete_matrix_data==true)
      {
       oomph_info << "MUMPS is clearing the matrix" << std::endl;
       cr_matrix_pt->clear();
      }

     // Call mumps factorisation
     mumps_factorise(Solver_ID_in_pool,n,nnz_loc,
                     &irn_loc[0],&jcn_loc[0],&a_loc[0]);
     
    }
   // else the CRDoubleMatrix is not distributed
   else
    {

     std::ostringstream error_message_stream;
     error_message_stream << "MumpsSolver only works for a "
                          << " distributed CRDoubleMatrix\n";
     throw OomphLibError(error_message_stream.str(),
                         "MumpsSolver::factorise()",
                         OOMPH_EXCEPTION_LOCATION);

    }
  }
 // Otherwise throw an error
 else
  {
   std::ostringstream error_message_stream;
   error_message_stream << "MumpsSolver implemented only for "
                        << "distributed CRDoubleMatrix. \n";
   throw OomphLibError(error_message_stream.str(),
                       "MumpsSolver::factorise()",
                       OOMPH_EXCEPTION_LOCATION);
  }


 // Switch off docing again
 mumps_switch_off_doc(Solver_ID_in_pool);
 

}

//=============================================================================
/// Do the backsubstitution for mumps solver. 
/// This does not make any assumption about the distribution of the
/// vectors
//=============================================================================
 void MumpsSolver::backsub(const DoubleVector &rhs,
                           DoubleVector &result)
 {
  
  // Doc stats?
  if (Doc_stats) mumps_switch_on_doc(Solver_ID_in_pool);

  // number of DOFs
  int ndof = Distribution_pt->nrow();
  
  // Make backup to avoid over-writing
  DoubleVector tmp_rhs;
  tmp_rhs=rhs;
  
  // Now turn this into a global (non-distributed) vector
  // because that's what mumps needs
  
  // Make a global distribution (i.e. one that isn't distributed)
  LinearAlgebraDistribution global_distribution(
   Distribution_pt->communicator_pt(),ndof,false);
   
  // Redistribute the tmp_rhs vector with this distribution -- it's
  // now "global", as required for mumps
  tmp_rhs.redistribute(&global_distribution);
  
  // Do the backsubsitution phase -- overwrites the tmp_rhs vector with the
  // solution
  mumps_backsub(Solver_ID_in_pool,ndof,&tmp_rhs[0]);
  
  // Broadcast the result which is only held on root
  MPI_Bcast(&tmp_rhs[0],ndof,MPI_DOUBLE,0,
            Distribution_pt->communicator_pt()->mpi_comm());

  // If the result vector is distributed, re-distribute the
  // non-distributed tmp_rhs vector to match
  if (result.distribution_setup()) 
   {
    tmp_rhs.redistribute(result.distribution_pt());
   }
  else
   {
    tmp_rhs.redistribute(Distribution_pt);    
   }
  
  // Now copy the tmp_rhs vector into the (matching) result
  result = tmp_rhs;
 

  // Switch off docing again
  mumps_switch_off_doc(Solver_ID_in_pool);

}


//===============================================================
/// Clean up the memory allocated for the MumpsSolver solver
//===============================================================
void MumpsSolver::clean_up_memory()
{
 //Cleanup
 mumps_cleanup_memory(Solver_ID_in_pool);
}


//=========================================================================
/// Linear-algebra-type solver: Takes pointer to a matrix and rhs 
/// vector and returns the solution of the linear system. Problem pointer 
/// defaults to NULL and can be omitted. The function returns the global 
/// result vector. Matrix must be CRDoubleMatrix.
/// Note: if Delete_matrix_data is true the function 
/// matrix_pt->clean_up_memory() will be used to wipe the matrix data.
//=========================================================================
void MumpsSolver::solve(DoubleMatrixBase* const &matrix_pt,
                        const DoubleVector &rhs,
                        DoubleVector &result)
{

 // Initialise timer
 double t_start = TimingHelpers::timer(); 
 
 // Doc stats?
 if (Doc_stats) mumps_switch_on_doc(Solver_ID_in_pool);

#ifdef PARANOID
 // check that the rhs vector is setup
 if (!rhs.distribution_pt()->setup())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The vectors rhs must be setup";
   throw OomphLibError(error_message_stream.str(),
                       "MumpsSolver::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // check that the matrix is square
 if (matrix_pt->nrow() != matrix_pt->ncol())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The matrix at matrix_pt must be square.";
   throw OomphLibError(error_message_stream.str(),
                       "MumpsSolver::solve()",
                       OOMPH_EXCEPTION_LOCATION);    
  }

 // check that the matrix and the rhs vector have the same nrow()
 if (matrix_pt->nrow() != rhs.nrow())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The matrix and the rhs vector must have the same number of rows.";
   throw OomphLibError(error_message_stream.str(),
                       "MumpsSolver::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 
 // Setup the distribution of the solver to match that of the matrix
 if (dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt))
  {
   Distribution_pt->rebuild(dynamic_cast<DistributableLinearAlgebraObject*>
                            (matrix_pt)->distribution_pt());
  }

 //Factorise the matrix
 factorise(matrix_pt);
 
 //Now do the back solve
 backsub(rhs,result);

 // Doc time for solve
 double t_end = TimingHelpers::timer(); 
 Solution_time = t_end-t_start;
 
 if ((Doc_time) && (MPI_Helpers::My_rank==0))
  {
   
   oomph_info << std::endl 
              << "Time for MumpsSolver solve [sec]       : "
              << t_end-t_start << std::endl << std::endl;
  }
 
 // If we are not storing the solver data for resolves, delete it
 if (!Enable_resolve) 
  {
   clean_up_memory();
  }


 // Switch off docing again
 mumps_switch_off_doc(Solver_ID_in_pool);
 
 }

//==================================================================
/// Solver: Takes pointer to problem and returns the results Vector
/// which contains the solution of the linear system defined by
/// the problem's fully assembled Jacobian and residual Vector.
//==================================================================
void MumpsSolver::solve(Problem* const &problem_pt, DoubleVector &result)
{

 // Initialise timer
 double t_start = TimingHelpers::timer();

 // Doc stats?
 if (Doc_stats) mumps_switch_on_doc(Solver_ID_in_pool);
  
 // number of dofs
 unsigned n_dof = problem_pt->ndof();

 // Set the distribution for the solver.
 Distribution_pt->rebuild(problem_pt->communicator_pt(),n_dof);
  
 // Take a copy of Delete_matrix_data
 bool copy_of_Delete_matrix_data = Delete_matrix_data;
 
 // Set Delete_matrix to true
 Delete_matrix_data = true;
 
 // Initialise timer
 t_start = TimingHelpers::timer();
 
 // Storage for the distributed residuals vector
 DoubleVector residuals(Distribution_pt,0.0);
 
 // Get the sparse jacobian and residuals of the problem
 CRDoubleMatrix jacobian(Distribution_pt);
 problem_pt->get_jacobian(residuals, jacobian);
 
 // Doc time for setup
 double t_end = TimingHelpers::timer();
 Jacobian_setup_time = t_end-t_start;
 if ((Doc_time) && (MPI_Helpers::My_rank==0))
  {
   oomph_info << "Time to set up CRDoubleMatrix Jacobian [sec]        : "
              << Jacobian_setup_time << std::endl;
  }

 //Now call the linear algebra solve, if desired
 if(!Suppress_solve) 
  {
   solve(&jacobian,residuals,result);
  }
 
 // Set Delete_matrix back to original value
 Delete_matrix_data = copy_of_Delete_matrix_data;
 
 // Switch off docing again
 mumps_switch_off_doc(Solver_ID_in_pool);

 // Finalise/doc timings
 if ((Doc_time) && (MPI_Helpers::My_rank==0))
  {
   double t_end = TimingHelpers::timer();
   oomph_info << std::endl << "Total time for MumpsSolver " << "(np=" 
              << MPI_Helpers::Nproc << ",N=" << problem_pt->ndof()
              <<") [sec] : " << t_end-t_start << std::endl << std::endl;
  }
 
}


//===============================================================
/// Resolve the system defined by the last assembled jacobian
/// and the specified rhs vector if resolve has been enabled.
/// Note: returns the global result Vector.
//===============================================================
void MumpsSolver::resolve(const DoubleVector &rhs, DoubleVector &result)
{
#ifdef PARANOID
 // check that the rhs vector is setup
 if (!rhs.distribution_pt()->setup())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The vectors rhs must be setup";
   throw OomphLibError(error_message_stream.str(),
                       "MumpsSolver::resolve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check that the rhs distribution is the same as the distribution as this 
 // solver
 if (!(*rhs.distribution_pt() == *Distribution_pt))
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The distribution of rhs vector must match the solver";
   throw OomphLibError(error_message_stream.str(),
                       "MumpsSolver::resolve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // if the result vector is setup then check it has the same distribution
 // as the rhs
 if (result.distribution_setup())
  {
   if (!(*result.distribution_pt() == *rhs.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The result vector distribution has been setup; it must have the "
      << "same distribution as the rhs vector.";
     throw OomphLibError(error_message_stream.str(),
                         "MumpsSolver::resolve()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }   
#endif
 
 // Store starting time for solve
 double t_start = TimingHelpers::timer();
 

 // Doc stats?
 if (Doc_stats) mumps_switch_on_doc(Solver_ID_in_pool);
  
 //Now do the back substitution phase
 backsub(rhs,result);

 // Doc time for solve
 double t_end = TimingHelpers::timer();
 Solution_time = t_end-t_start;

 // Switch off docing again
 mumps_switch_off_doc(Solver_ID_in_pool);
 
 if ((Doc_time) && (MPI_Helpers::My_rank==0))
  {
   oomph_info << "Time for MumpsSolver solve [sec]: "
              << t_end-t_start << std::endl;
  }
 
}




}
