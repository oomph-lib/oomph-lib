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
#include "hypre_solver.h"


namespace oomph
{
 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
 
//==================================================================
/// Default settings for various uses of the HYPRE solver
//==================================================================
 namespace Hypre_default_settings
 {
  /// \short Set default parameters for use as preconditioner in
  /// 2D Poisson-type problem.
  void set_defaults_for_2D_poisson_problem(
   HyprePreconditioner* hypre_preconditioner_pt)
  {
   // Use simple smoother
   hypre_preconditioner_pt->amg_using_simple_smoothing();
   
   // Smoother types:
   //           0=Jacobi
   //           1=Gauss-Seidel
   hypre_preconditioner_pt->amg_simple_smoother() = 1;
   
   // AMG preconditioner
   hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
   
   // Choose strength parameter for amg
   hypre_preconditioner_pt->amg_strength() = 0.25;
  }
  

  /// \short Set default parameters for use as preconditioner in
  /// 3D Poisson-type problem.
  void set_defaults_for_3D_poisson_problem(
   HyprePreconditioner* hypre_preconditioner_pt)
  {
   // Set default settings as for 2D Poisson
   set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
   
   // Change strength parameter for amg
   hypre_preconditioner_pt->amg_strength() = 0.7;
  }
  
  
  /// \short Set default parameters for use as preconditioner in
  /// for momentum block in Navier-Stokes problem
  void set_defaults_for_navier_stokes_momentum_block(
   HyprePreconditioner* hypre_preconditioner_pt)
  {
   // Set default settings as for 2D Poisson
   set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
   
   // Change smoother type:
   //           0=Jacobi
   //           1=Gauss-Seidel
   hypre_preconditioner_pt->amg_simple_smoother() = 0;
    
   // Set smoother damping
   hypre_preconditioner_pt->amg_damping() = 0.5;
   
   // Change strength parameter for amg
   hypre_preconditioner_pt->amg_strength() = 0.75;
  }
  
 }


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// functions for HypreHelpers namespace
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

 namespace HypreHelpers
 {
  
//========================================================================
/// Helper function to check the Hypre error flag, return the message
/// associated with any error, and reset the global error flag to zero.
/// This function also returns the error value.
//========================================================================
  int check_HYPRE_error_flag(std::ostringstream& message)
  {
   // get the Hypre error flag
   int err = HYPRE_GetError();
   
   // Tell us all about it...
   if (err)
    {
     oomph_info << "Hypre error flag=" << err << std::endl;
     char* error_message = new char[128];
     HYPRE_DescribeError(err, error_message);
     message << "WARNING: " << std::endl
             << "HYPRE error message: " << error_message 
             << std::endl;
     delete[] error_message;
    }
   
   // reset Hypre's global error flag
   hypre__global_error=0;
   
   return err;
 }
 
 
//========================================================================
/// Helper function to create a serial HYPRE_IJVector and
/// HYPRE_ParVector. The length of the vector created is
/// n_values. An array of values (of length n_values) can be
/// inserted, or if (values==0) an empty Hypre vector is created.
/// indices can define the index of the entries in values, or if
/// (indices==0) values are inserted in the order they occur.
//========================================================================
  void create_HYPRE_Vector(const int& n_values,
                           HYPRE_IJVector& hypre_ij_vector,
                           HYPRE_ParVector& hypre_par_vector,
                           const double* values,
                           int* indices)
  {
   // initialize Hypre vector
   HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, n_values-1, &hypre_ij_vector);
   HYPRE_IJVectorSetObjectType(hypre_ij_vector, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(hypre_ij_vector);

   // keep track of whether we need to delete any arrays
   bool delete_indices = false;

   // insert values if required
   if (values != 0)
    {
     // generate an array of indices if non exists
     if (indices == 0)
      {
       delete_indices = true;
       // set up array containing indices
       indices = new int[n_values];
       for (int i=0; i<n_values; i++)
        {
         indices[i] = i;
        }
      }

     // insert values
     HYPRE_IJVectorSetValues(hypre_ij_vector, n_values, indices, values);

     // tidy up memory if required
     if (delete_indices)
      {
       delete[] indices;
       indices = 0;
      }
    }

   // assemble vectors
   HYPRE_IJVectorAssemble(hypre_ij_vector);
   HYPRE_IJVectorGetObject(hypre_ij_vector, (void **) &hypre_par_vector);
  }
  

//========================================================================
/// Helper function to create a serial HYPRE_IJMatrix and
/// HYPRE_ParCSRMatrix from a CRDoubleMatrix
//========================================================================
  void create_HYPRE_Matrix(const CRDoubleMatrix& oomph_matrix,
                           HYPRE_IJMatrix& hypre_ij_matrix,
                           HYPRE_ParCSRMatrix& hypre_par_matrix)
  {
   // find number of rows/columns
   const int n_rows = int(oomph_matrix.nrow());
   const int n_cols = int(oomph_matrix.ncol());

   // get pointers to the matrix

   // column indicies of matrix
   const int* matrix_cols = oomph_matrix.column_index();

   // entries of matrix
   const double* matrix_vals = oomph_matrix.value();

   // row starts
   const int* matrix_row_start = oomph_matrix.row_start();

   // initialize hypre matrix
   HYPRE_IJMatrixCreate(MPI_COMM_WORLD,
                        0, n_rows-1, 0, n_cols-1,
                        &hypre_ij_matrix);
   HYPRE_IJMatrixSetObjectType(hypre_ij_matrix, HYPRE_PARCSR);
   HYPRE_IJMatrixInitialize(hypre_ij_matrix);

   // set up simple row_map for serial case
   int* ncols_per_row = new int[n_rows];  // Number of columns in each row
   int* row_map = new int[n_rows];        // Map of local to global rows
   for (int i=0; i<n_rows; i++)
    {
     ncols_per_row[i] = matrix_row_start[i+1] - matrix_row_start[i];
     // simple row map for serial matrix
     row_map[i] = i;
    }

   // put values in HYPRE matrix
   HYPRE_IJMatrixSetValues(hypre_ij_matrix,
                           n_rows,
                           ncols_per_row,
                           row_map,
                           matrix_cols,
                           matrix_vals);

   // tidy up memory
   delete[] ncols_per_row;
   delete[] row_map;

   // assemble matrix
   HYPRE_IJMatrixAssemble(hypre_ij_matrix);
   HYPRE_IJMatrixGetObject(hypre_ij_matrix, (void **) &hypre_par_matrix);
  }


#ifdef OOMPH_HAS_MPI
//========================================================================
/// Helper function to create a distributed HYPRE_IJVector and
/// HYPRE_ParVector. lower and upper define the ranges of the contiguous
/// partitioning of the vector. An array of local values
/// (of length 1+upper-lower) can be inserted, or if (values==0) an empty
/// vector is created. indices can be used to define the global index
/// of the entries in values, or if (indices==0) values are inserted in
/// the order they occur.
//========================================================================
  void create_HYPRE_Vector(const int& lower,
                           const int& upper,
                           HYPRE_IJVector& hypre_ij_vector,
                           HYPRE_ParVector& hypre_par_vector,
                           const double* values,
                           int* indices)
  {
   // initialize Hypre vector
   HYPRE_IJVectorCreate(MPI_COMM_WORLD, lower, upper, &hypre_ij_vector);
   HYPRE_IJVectorSetObjectType(hypre_ij_vector, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(hypre_ij_vector);
   
   // insert values if required
   if (values != 0)
    {
     // keep track of whether we need to point indices to zero at end
     bool nullify_indices = false;
     
     // find out how many values there are
     const int n_values = 1+upper-lower;
     
     // generate Vector of indices if non exists
     if (indices == 0)
      {
       Vector<int> vec_indices(n_values);
              
       // set up array containing indices
       indices = &vec_indices[0];
       nullify_indices = true;
       for (int i=0; i<n_values; i++)
        {
         indices[i] = lower+i;
        }
      }
     
     // insert values
     HYPRE_IJVectorSetValues(hypre_ij_vector, n_values, indices, values);
     
     // point indecies to null if required
     if (nullify_indices)
      {
        indices = 0;
      }
    }
   
   // assemble vectors
   HYPRE_IJVectorAssemble(hypre_ij_vector);
   HYPRE_IJVectorGetObject(hypre_ij_vector, (void **) &hypre_par_vector);
  }
    
  
//========================================================================
/// Helper function to create a distributed HYPRE_IJMatrix and
/// HYPRE_ParCSRMatrix from a square CRDoubleMatrix. lower and upper define
/// the range of the contiguous row partitioning of the Hypre matrix created.
//========================================================================
  void create_HYPRE_Matrix(const CRDoubleMatrix& oomph_matrix,
                           const int& lower,
                           const int& upper,
                           HYPRE_IJMatrix& hypre_ij_matrix,
                           HYPRE_ParCSRMatrix& hypre_par_matrix)
  {
#ifdef PARANOID
   // check the matrix is square
   if ( oomph_matrix.nrow() != oomph_matrix.ncol() )
    {
     std::ostringstream error_message;
     error_message << "create_HYPRE_Matrix require a square matrix. "
                   << "Matrix is " << oomph_matrix.nrow()
                   << " by " << oomph_matrix.ncol() << std::endl;
     throw OomphLibError(error_message.str(),
                         "HypreHelpers::create_HYPRE_Matrix()",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   
   // get pointers to the matrix
   
   // column indicies of matrix
   const int* matrix_cols = oomph_matrix.column_index();
   
   // entries of matrix
   const double* matrix_vals = oomph_matrix.value();
   
   // row starts
   const int* matrix_row_starts = oomph_matrix.row_start();
   
   // initialize hypre matrix
   HYPRE_IJMatrixCreate(MPI_COMM_WORLD,
                        lower,
                        upper,
                        lower,
                        upper,
                        &hypre_ij_matrix);
   HYPRE_IJMatrixSetObjectType(hypre_ij_matrix, HYPRE_PARCSR);
   HYPRE_IJMatrixInitialize(hypre_ij_matrix);
   
   // find number of local rows
   const int n_local_rows = 1+upper-lower;
   
   // Vector to store number of columns in each row
   Vector<int> ncols_per_row(n_local_rows);
   
   // Vector to store map of local to global rows
   Vector<int> row_map(n_local_rows);
   
   // find the number of entries per row and set up the row map
   for (int i=0; i<n_local_rows; i++)
    {
     int global_row = lower+i;
     ncols_per_row[i] = matrix_row_starts[global_row+1] - 
      matrix_row_starts[global_row];
     row_map[i] = global_row;
    }
   
   // find where this row block starts in matrix_cols and matrix_vals
   int local_start = matrix_row_starts[lower];
   
   // put values in HYPRE matrix
   HYPRE_IJMatrixSetValues(hypre_ij_matrix,
                           n_local_rows,
                           &ncols_per_row[0],
                           &row_map[0],
                           matrix_cols+local_start,
                           matrix_vals+local_start);
   
   // assemble matrix
   HYPRE_IJMatrixAssemble(hypre_ij_matrix);
   HYPRE_IJMatrixGetObject(hypre_ij_matrix, (void **) &hypre_par_matrix);
  }
  
  
//========================================================================
/// Helper function to create a distributed HYPRE_IJMatrix and
/// HYPRE_ParCSRMatrix from a DistributedCRDoubleMatrix.
//========================================================================
  void create_HYPRE_Matrix(DistributedCRDoubleMatrix& oomph_matrix,
                           HYPRE_IJMatrix& hypre_ij_matrix,
                           HYPRE_ParCSRMatrix& hypre_par_matrix)
  {
#ifdef PARANOID
   // check the matrix is square
   if ( oomph_matrix.nrow() != oomph_matrix.ncol() )
    {
     std::ostringstream error_message;
     error_message << "create_HYPRE_Matrix require a square matrix. "
                   << "Matrix is " << oomph_matrix.nrow()
                   << " by " << oomph_matrix.ncol() << std::endl;
     throw OomphLibError(error_message.str(),
                         "HypreHelpers::create_HYPRE_Matrix()",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   
   // get pointers to the matrix
   
   // column indicies of matrix
   const int* matrix_cols = oomph_matrix.column_index();
   
   // entries of matrix
   const double* matrix_vals = oomph_matrix.value();
   
   // row starts
   const int* matrix_row_starts = oomph_matrix.row_start();
   
   // find number of local rows
   const int n_local_rows = oomph_matrix.nrow_local();
   
   // find the first row on this processor
   int first_row = oomph_matrix.first_row();
   
   // initialize hypre matrix
   HYPRE_IJMatrixCreate(MPI_COMM_WORLD,
                        first_row,
                        first_row + n_local_rows - 1,
                        first_row,
                        first_row + n_local_rows - 1,
                        &hypre_ij_matrix);
   HYPRE_IJMatrixSetObjectType(hypre_ij_matrix, HYPRE_PARCSR);
   HYPRE_IJMatrixInitialize(hypre_ij_matrix);
   
   // Vector to store number of columns in each row
   Vector<int> ncols_per_row(n_local_rows);
   
   // Vector to store map of local to global rows
   Vector<int> row_map(n_local_rows);
   
   // find the number of entries per row and set up the row map
   for (int i=0; i<n_local_rows; i++)
    {
     ncols_per_row[i] = matrix_row_starts[i+1] - matrix_row_starts[i];
     row_map[i] = first_row+i;
    }
   
   // put values in HYPRE matrix
   HYPRE_IJMatrixSetValues(hypre_ij_matrix,
                           n_local_rows,
                           &ncols_per_row[0],
                           &row_map[0],
                           matrix_cols,
                           matrix_vals);
   
   // assemble matrix
   HYPRE_IJMatrixAssemble(hypre_ij_matrix);
   HYPRE_IJMatrixGetObject(hypre_ij_matrix, (void **) &hypre_par_matrix);
  }

#endif  
 }


 
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// functions for HypreInterface class
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//=============================================================================
/// Helper function which creates a Hypre matrix from a CRDoubleMatrix
/// If OOMPH-LIB has been set up for MPI use, the Hypre matrix is 
/// distributed over the available processors.
//=============================================================================
 void HypreInterface::hypre_matrix_setup(CRDoubleMatrix* matrix_pt)
 {
  // reset Hypre's global error flag
  hypre__global_error=0;
  
#ifdef OOMPH_HAS_MPI
  // MPI version
  double t_start = MPI_Wtime();
  
  // Nrows_local stores number of rows on each processor -
  // needed later to gather the distributed solution using MPI
  Nrows_local.resize(MPI_Helpers::Nproc);
  
  // Find the number of rows
  Nrow = matrix_pt->nrow();
  
  // issue warning if the matrix is small compared to the number of processors
  if ( 2*MPI_Helpers::Nproc>int(Nrow) )
   {
    oomph_info
     << "Warning: HYPRE based solvers may fail if 2*number of processors "
     << "is greater than the number of unknowns!" << std::endl;
   }
  
  // Set the row partitioning for the processors and store
  double size = double(Nrow)/MPI_Helpers::Nproc;
  int first_row = 0;
  int next_first_row = 0;
  for (int p=0;p<MPI_Helpers::Nproc; p++)
   {
    // find the first_row on the next processor unless there
    // are no more processors, in which case set it to be
    // the number of global rows
    if (p==MPI_Helpers::Nproc)
     {
      next_first_row = int(Nrow);
     }
    else
     {
      next_first_row = int((p+1)*size);
     }
    
    if(p==MPI_Helpers::My_rank)
     {
      First_row = first_row;
      Last_row = next_first_row-1;
     }
    // set Nrow_local
    Nrows_local[p] = next_first_row-first_row;
    
    // set first row for next processor
    first_row = next_first_row;
   }
  
  // generate the Hypre matrix
  HypreHelpers::create_HYPRE_Matrix(*matrix_pt,
                                    First_row,
                                    Last_row,
                                    Matrix_ij,
                                    Matrix_par);
#else
  // Non-MPI version
  clock_t t_start = clock();
  
  // find number of rows in matrix
  Nrow = matrix_pt->nrow();
  
  // set up hypre matrix
  // -------------------
  HypreHelpers::create_HYPRE_Matrix(*matrix_pt,
                                    Matrix_ij,
                                    Matrix_par);
#endif
  
  // Output error messages if required  
  if (Hypre_error_messages)
   {
    std::ostringstream message;
    int err = HypreHelpers::check_HYPRE_error_flag(message);
    if (err)
     {
      OomphLibWarning(message.str(),
                      "HypreSolver::hypre_matrix_setup()",
                      OOMPH_EXCEPTION_LOCATION);
     }
   }
  
  // delete CRDoubleMatrix if required
  if (Delete_input_data)
   {
    matrix_pt->clean_up_memory();
   }

  // Output times
  if (Output_time)
   {
    oomph_info << "Time to generate HYPRE matrix [s]          : ";
#ifdef OOMPH_HAS_MPI   
    double t_end = MPI_Wtime();
    oomph_info << t_end-t_start << std::endl; 
#else
    clock_t t_end = clock();
    oomph_info << double(t_end-t_start)/CLOCKS_PER_SEC << std::endl;
#endif
   }

 }



#ifdef OOMPH_HAS_MPI   
//=============================================================================
/// Helper function which creates a distributed Hypre matrix from an 
/// OOMPH-LIB DistributedCRDoubleMatrix
//=============================================================================
 void HypreInterface::hypre_matrix_setup(DistributedCRDoubleMatrix* matrix_pt)
 {
  double t_start = MPI_Wtime();

  // reset Hypre's global error flag
  hypre__global_error=0;
  
  // find the matrix partitioning and store
  First_row = matrix_pt->first_row();
  int nrow_loc = matrix_pt->nrow_local();
  Last_row = First_row + nrow_loc - 1;
  
  Nrow = matrix_pt->nrow();
  
  // set up Nrow_local - use MPI_Allgather to get the values
  Nrows_local.resize(MPI_Helpers::Nproc);
  MPI_Allgather(&nrow_loc,
                1,
                MPI_INT,
                &Nrows_local[0],
                1,
                MPI_INT,
                MPI_COMM_WORLD);
  
  // create the Hypre Matrix
  HypreHelpers::create_HYPRE_Matrix(*matrix_pt,
                                    Matrix_ij,
                                    Matrix_par);

  // check error flag
  if (Hypre_error_messages)
   {
    std::ostringstream message;
    int err = HypreHelpers::check_HYPRE_error_flag(message);
    if (err)
     {
      OomphLibWarning(message.str(),
                      "HypreSolver::hypre_matrix_setup()",
                      OOMPH_EXCEPTION_LOCATION);
     }
   }
 
  // delete DistributedCRDoubleMatrix matrix if required
  if (Delete_input_data)
   {
    matrix_pt->clean_up_memory();
   }

  // output times
  if (Output_time)
   {
    double t_end = MPI_Wtime();
    oomph_info << "Time to generate HYPRE matrix [s]          : "
               << t_end - t_start << std::endl;
   }

 }
#endif


//=============================================================================
/// Sets up the solver data required for use in an oomph-lib
/// LinearSolver or Preconditioner, once the Hypre matrix has been
/// generated using hypre_matrix_setup(...).
//=============================================================================
 void HypreInterface::hypre_solver_setup()
 {
  // Store time
#ifdef OOMPH_HAS_MPI   
    double t_start = MPI_Wtime();
    double t_end = 0;
#else
    clock_t t_start = clock();
    clock_t t_end = 0; 
#endif

  // reset Hypre's global error flag
  hypre__global_error=0;

  // create dummy Hypre vectors which are required for setup
  HYPRE_IJVector dummy_sol_ij;
  HYPRE_ParVector dummy_sol_par;
  HYPRE_IJVector dummy_rhs_ij;
  HYPRE_ParVector dummy_rhs_par;
#ifdef OOMPH_HAS_MPI 
  HypreHelpers::create_HYPRE_Vector(First_row,
                                    Last_row,
                                    dummy_sol_ij,
                                    dummy_sol_par);
  HypreHelpers::create_HYPRE_Vector(First_row,
                                    Last_row,
                                    dummy_rhs_ij,
                                    dummy_rhs_par);
#else
  HypreHelpers::create_HYPRE_Vector(int(Nrow),
                                    dummy_sol_ij,
                                    dummy_sol_par);
  HypreHelpers::create_HYPRE_Vector(int(Nrow),
                                    dummy_rhs_ij,
                                    dummy_rhs_par);
#endif

  // Set up internal preconditioner for CG, GMRES or BiCGSTAB
  // --------------------------------------------------------
  if ( (Hypre_method>=CG) && (Hypre_method<=BiCGStab) )
   {
    // AMG preconditioner
    if (Internal_preconditioner==BoomerAMG)
     {
      // set up BoomerAMG
      HYPRE_BoomerAMGCreate(&Preconditioner);
      HYPRE_BoomerAMGSetPrintLevel(Preconditioner, AMG_print_level);
      HYPRE_BoomerAMGSetMaxLevels(Preconditioner, AMG_max_levels);
      HYPRE_BoomerAMGSetMaxIter(Preconditioner, 1);
      HYPRE_BoomerAMGSetTol(Preconditioner, 0.0);
      HYPRE_BoomerAMGSetCoarsenType(Preconditioner, AMG_coarsening);
      HYPRE_BoomerAMGSetStrongThreshold(Preconditioner, AMG_strength);
      HYPRE_BoomerAMGSetMaxRowSum(Preconditioner, AMG_max_row_sum);
      HYPRE_BoomerAMGSetTruncFactor(Preconditioner, AMG_truncation);

      if (AMG_using_simple_smoothing)
       {
        HYPRE_BoomerAMGSetRelaxType(Preconditioner, AMG_simple_smoother);
        HYPRE_BoomerAMGSetNumSweeps(Preconditioner, AMG_smoother_iterations);

        // This one gives a memory leak
        //double * relaxweight = new double[AMG_max_levels];

        // This is how they do it in a hypre demo code
        double* relaxweight = hypre_CTAlloc(double, AMG_max_levels);

        for (unsigned i=0; i<AMG_max_levels; i++)
         {
          relaxweight[i] = AMG_damping;
         }
        HYPRE_BoomerAMGSetRelaxWeight(Preconditioner, relaxweight);
       }
      else
       {
        HYPRE_BoomerAMGSetSmoothType(Preconditioner, AMG_complex_smoother);
        HYPRE_BoomerAMGSetSmoothNumLevels(Preconditioner, AMG_max_levels);
        HYPRE_BoomerAMGSetSmoothNumSweeps(Preconditioner,
                                          AMG_smoother_iterations);
       }       

      Existing_preconditioner = BoomerAMG;
     }

    // Euclid preconditioner
    else if (Internal_preconditioner==Euclid)
     {
      HYPRE_EuclidCreate(MPI_COMM_WORLD, &Preconditioner);

      // set Euclid parameters using command line like array
      int n_args = 0;
      char* args[22];

      // first argument is empty string
      args[n_args++] = "";
  
      // switch on/off Block Jacobi ILU
      args[n_args++] = "-bj";
      if (Euclid_using_BJ)
      {
       args[n_args++] = "1";
      }
      else
      {
       args[n_args++] = "0";
      }

      // switch on/off row scaling
      args[n_args++] = "-rowScale";
      if (Euclid_rowScale)
      {
       args[n_args++] = "1";
      }
      else
      {
       args[n_args++] = "0";
      }

      // set level for ILU(k)
      args[n_args++] = "-level";
      char level_value[10];
      sprintf(level_value,"%d",Euclid_level);
      args[n_args++] = level_value;

      // set drop tol for ILU(k) factorization
      args[n_args++] = "-sparseA";
      char droptol[20];
      sprintf(droptol,"%f",Euclid_droptol);
      args[n_args++] = droptol;

      // set ILUT factorization if required
      if (Euclid_using_ILUT)
      {
       args[n_args++] = "-ilut";
       args[n_args++] = droptol;
      }

      // set printing of Euclid data
      if (Euclid_print_level == 0)
      {
       args[n_args++] = "-eu_stats";
       args[n_args++] = "0";
       args[n_args++] = "-eu_mem";
       args[n_args++] = "0";
      }
      if (Euclid_print_level == 1)
      {
       args[n_args++] = "-eu_stats";
       args[n_args++] = "1";
       args[n_args++] = "-eu_mem";
       args[n_args++] = "0";
      }
      if (Euclid_print_level == 2)
      {
       args[n_args++] = "-eu_stats";
       args[n_args++] = "1";
       args[n_args++] = "-eu_mem";
       args[n_args++] = "1";
      }

      // set next entry in array to null
      args[n_args] = 0;

      HYPRE_EuclidSetParams(Preconditioner, n_args, args);

      Existing_preconditioner = Euclid;
     }

    // ParaSails preconditioner
    else if (Internal_preconditioner==ParaSails)
     {
      HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &Preconditioner);
      HYPRE_ParaSailsSetSym(Preconditioner, ParaSails_symmetry);
      HYPRE_ParaSailsSetParams(Preconditioner,
                               ParaSails_thresh,
                               ParaSails_nlevel);
      HYPRE_ParaSailsSetFilter(Preconditioner, ParaSails_filter);
      Existing_preconditioner = ParaSails;
     }

    // check error flag
    if (Hypre_error_messages)
     {
      std::ostringstream message;
      int err = HypreHelpers::check_HYPRE_error_flag(message);
      if (err)
       {
        OomphLibWarning(message.str(),
                        "HypreSolver::hypre_setup()",
                        OOMPH_EXCEPTION_LOCATION);
       }
     }
   } // end of setting up internal preconditioner

  // Record preconditioner set up time
  double preconditioner_setup_time=0;
  if (Output_time)
   {
#ifdef OOMPH_HAS_MPI   
    t_end = MPI_Wtime();
    preconditioner_setup_time = t_end-t_start; 
#else
    t_end = clock();
    preconditioner_setup_time = double(t_end-t_start)/CLOCKS_PER_SEC;
#endif
   }

  // set up solver
  // -------------
  t_start = t_end;

  // AMG solver
  if (Hypre_method==BoomerAMG)
   {
    if (Output_time)
     {
      oomph_info << "Setting up BoomerAMG" << std::endl;
     }

    // set up BoomerAMG
    HYPRE_BoomerAMGCreate(&Solver);
    HYPRE_BoomerAMGSetPrintLevel(Solver, AMG_print_level);
    HYPRE_BoomerAMGSetMaxLevels(Solver, AMG_max_levels);
    HYPRE_BoomerAMGSetMaxIter(Solver, Max_iter);
    HYPRE_BoomerAMGSetTol(Solver, Tolerance);
    HYPRE_BoomerAMGSetCoarsenType(Solver, AMG_coarsening);
    HYPRE_BoomerAMGSetStrongThreshold(Solver, AMG_strength);
    HYPRE_BoomerAMGSetMaxRowSum(Solver, AMG_max_row_sum);
    HYPRE_BoomerAMGSetTruncFactor(Solver, AMG_truncation);
    
    if (AMG_using_simple_smoothing)
     {
      HYPRE_BoomerAMGSetRelaxType(Solver, AMG_simple_smoother);
      HYPRE_BoomerAMGSetNumSweeps(Solver, AMG_smoother_iterations);

      // This one gives a memory leak
      //double * relaxweight = new double[AMG_max_levels];

      // This is how they do it in a hypre demo code
      double* relaxweight = hypre_CTAlloc(double, AMG_max_levels);

      for (unsigned i=0; i<AMG_max_levels; i++)
       {
        relaxweight[i] = AMG_damping;
       }
      HYPRE_BoomerAMGSetRelaxWeight(Solver, relaxweight);
     }
    else
     {
      HYPRE_BoomerAMGSetSmoothType(Solver, AMG_complex_smoother);
      HYPRE_BoomerAMGSetSmoothNumLevels(Solver, AMG_max_levels);
      HYPRE_BoomerAMGSetSmoothNumSweeps(Solver, AMG_smoother_iterations);
     }  

    HYPRE_BoomerAMGSetup(Solver,
                         Matrix_par,
                         dummy_rhs_par,
                         dummy_sol_par);

    Existing_solver = BoomerAMG;
   }

  // Euclid solver
  else if (Hypre_method==Euclid)
   {
    if (Output_time)
     {
      oomph_info << "Setting up Euclid" << std::endl;
     }
    HYPRE_EuclidCreate(MPI_COMM_WORLD, &Solver);

    // set Euclid parameters using command line like array
    int n_args = 0;
    char* args[20];

    // first argument is empty string
    args[n_args++] = "";

    // switch on/off Block Jacobi ILU
    args[n_args++] = "-bj";
    if (Euclid_using_BJ)
     {
      args[n_args++] = "1";
     }
    else
     {
      args[n_args++] = "0";
     }

    // switch on/off row scaling
    args[n_args++] = "-rowScale";
    if (Euclid_rowScale)
    {
     args[n_args++] = "1";
    }
    else
    {
     args[n_args++] = "0";
    }

    // set level for ILU(k)
    args[n_args++] = "-level";
    char level_value[10];
    sprintf(level_value,"%d",Euclid_level);
    args[n_args++] = level_value;

    // set drop tol for ILU(k) factorization
    args[n_args++] = "-sparseA";
    char droptol[20];
    sprintf(droptol,"%f",Euclid_droptol);
    args[n_args++] = droptol;

    // set ILUT factorization if required
    if (Euclid_using_ILUT)
    {
     args[n_args++] = "-ilut";
     args[n_args++] = droptol;
    }

    // set printing of Euclid data
    if (Euclid_print_level == 0)
    {
     args[n_args++] = "-eu_stats";
     args[n_args++] = "0";
     args[n_args++] = "-eu_mem";
     args[n_args++] = "0";
    }

    if (Euclid_print_level == 1)
    {
     args[n_args++] = "-eu_stats";
     args[n_args++] = "1";
     args[n_args++] = "-eu_mem";
     args[n_args++] = "0";
    }
    if (Euclid_print_level == 2)
    {
     args[n_args++] = "-eu_stats";
     args[n_args++] = "1";
     args[n_args++] = "-eu_mem";
     args[n_args++] = "1";
    }

    // set next entry in array to null
    args[n_args] = 0;

    HYPRE_EuclidSetParams(Solver, n_args, args);

    HYPRE_EuclidSetup(Solver,
                      Matrix_par,
                      dummy_rhs_par,
                      dummy_sol_par);
    Existing_solver = Euclid;
   }

  // ParaSails preconditioner
  else if (Hypre_method==ParaSails)
   {
    if (Output_time) 
     {
      oomph_info << "Setting up ParaSails" << std::endl;
     }

    HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &Solver);
    HYPRE_ParaSailsSetSym(Solver, ParaSails_symmetry);
    HYPRE_ParaSailsSetParams(Solver,
                             ParaSails_thresh,
                             ParaSails_nlevel);
    HYPRE_ParaSailsSetFilter(Solver, ParaSails_filter);

    HYPRE_ParaSailsSetup(Solver,
                         Matrix_par,
                         dummy_rhs_par,
                         dummy_sol_par);
    Existing_solver = ParaSails;
   }

  // CG solver
  else if (Hypre_method==CG)
   {
    if (Output_time)
     {
      oomph_info << "Setting up CG";
     }

    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &Solver);
    HYPRE_PCGSetTol(Solver, Tolerance);
    HYPRE_PCGSetLogging(Solver, 0);
    HYPRE_PCGSetPrintLevel(Solver, Krylov_print_level);
    HYPRE_PCGSetMaxIter(Solver, Max_iter);

    // set preconditioner
    if (Internal_preconditioner==BoomerAMG)  // AMG
     {
      if (Output_time)
       {
        oomph_info << " with BoomerAMG preconditioner";
       }

      HYPRE_PCGSetPrecond(Solver,
                         (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                         (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                         Preconditioner);
     }
    else if (Internal_preconditioner==Euclid)  // Euclid
     {
      if (Output_time) 
       {
        oomph_info << " with Euclid ILU preconditioner";
       }
      
      HYPRE_PCGSetPrecond(Solver,
                          (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
                          (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup,
                          Preconditioner);
     }
    else if (Internal_preconditioner==ParaSails)  // ParaSails
     {
      if (Output_time)
       {
        oomph_info << " with ParaSails approximate inverse preconditioner";
       }

      HYPRE_PCGSetPrecond(Solver,
                          (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
                          (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup,
                          Preconditioner);
     }
    else
     {
      if (Output_time) 
       {
        oomph_info << " with no preconditioner";
       }
     }
    if (Output_time) 
     {
      oomph_info << std::endl;
     }
    
    HYPRE_PCGSetup(Solver,
                   (HYPRE_Matrix) Matrix_par,
                   (HYPRE_Vector) dummy_rhs_par,
                   (HYPRE_Vector) dummy_sol_par);
    

    Existing_solver = CG;
   }

  // GMRES solver
  else if (Hypre_method==GMRES)
   {
    if (Output_time) 
     {
      oomph_info << "Setting up GMRES";
     }
    
    HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &Solver);
    HYPRE_GMRESSetTol(Solver, Tolerance);
    HYPRE_GMRESSetKDim(Solver, Max_iter);
    HYPRE_GMRESSetLogging(Solver, 0);
    HYPRE_GMRESSetPrintLevel(Solver, Krylov_print_level);
    HYPRE_GMRESSetMaxIter(Solver, Max_iter);

    // set preconditioner
    if (Internal_preconditioner==BoomerAMG)  // AMG
     {
      if (Output_time) 
       {
        oomph_info << " with BoomerAMG preconditioner";
       }

      HYPRE_GMRESSetPrecond(Solver,
                            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                            Preconditioner);
     }
    else if (Internal_preconditioner==Euclid) // Euclid
     {
      if (Output_time)
       {
        oomph_info << " with Euclid ILU preconditioner";
       }

      HYPRE_GMRESSetPrecond(Solver,
                            (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
                            (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup,
                            Preconditioner);
     }
    else if (Internal_preconditioner==ParaSails)  // ParaSails
     {
      if (Output_time) 
       {
        oomph_info << " with ParaSails approximate inverse preconditioner";
       }
      
      HYPRE_GMRESSetPrecond(Solver,
                            (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
                            (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup,
                            Preconditioner);
     }
    else
     {
      if (Output_time) 
       {
        oomph_info << " with no preconditioner";
       }
     }
    if (Output_time) 
     {
      oomph_info << std::endl;
     }
    
    HYPRE_GMRESSetup(Solver,
                     (HYPRE_Matrix) Matrix_par,
                     (HYPRE_Vector) dummy_rhs_par,
                     (HYPRE_Vector) dummy_sol_par);
        
    Existing_solver = GMRES;
   }
  
  // BiCGStab solver
  else if (Hypre_method==BiCGStab)
   {
    if (Output_time) 
     {
      oomph_info << "Setting up BiCGStab";
     }
    
    HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &Solver);
    HYPRE_BiCGSTABSetTol(Solver, Tolerance);
    HYPRE_BiCGSTABSetLogging(Solver, 0);
    HYPRE_BiCGSTABSetPrintLevel(Solver, Krylov_print_level);
    HYPRE_BiCGSTABSetMaxIter(Solver, Max_iter);
 
    // set preconditioner
    if (Internal_preconditioner==BoomerAMG)  // AMG
     {
      if (Output_time)
       {
        oomph_info << " with BoomerAMG preconditioner";
       }
      
      HYPRE_BiCGSTABSetPrecond(Solver,
                               (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                               (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                               Preconditioner);
     }
    else if (Internal_preconditioner==Euclid)  // Euclid
     {
      if (Output_time)
       {
        oomph_info << " with Euclid ILU preconditioner";
       }
      
      HYPRE_BiCGSTABSetPrecond(Solver,
                               (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
                               (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup,
                               Preconditioner);
     }
    else if (Internal_preconditioner==ParaSails)  // ParaSails
     {
      if (Output_time) 
       {
        oomph_info << " with ParaSails approximate inverse preconditioner";
       }
      
      HYPRE_BiCGSTABSetPrecond(Solver,
                               (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
                               (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup,
                               Preconditioner);
     }
    else
     {
      if (Output_time)
       {
        oomph_info << " with no preconditioner";
       }
     }
    if (Output_time) 
     {
      oomph_info << std::endl;
     }
    
    HYPRE_BiCGSTABSetup(Solver,
                        (HYPRE_Matrix) Matrix_par,
                        (HYPRE_Vector) dummy_rhs_par,
                        (HYPRE_Vector) dummy_sol_par);
    
    Existing_solver = BiCGStab;
   }

  // no solver exists for this value of Solver flag
  else
   {
    std::ostringstream error_message;
    error_message << "Solver has been set to an invalid value. "
                  << "current value=" << Solver;
    throw OomphLibError(error_message.str(),
                     	"Hypre_Solver::hypre_solver_setup()",
                        OOMPH_EXCEPTION_LOCATION);

   }

#ifdef OOMPH_HAS_MPI   
  t_end = MPI_Wtime();
  double solver_setup_time = t_end-t_start; 
#else
  t_end = clock();
  double solver_setup_time = double(t_end-t_start)/CLOCKS_PER_SEC;
#endif
  
  // destroy dummy hypre vectors
  HYPRE_IJVectorDestroy(dummy_sol_ij);
  HYPRE_IJVectorDestroy(dummy_rhs_ij);
  
  // check error flag
  if (Hypre_error_messages)
   {
    std::ostringstream message;
    int err = HypreHelpers::check_HYPRE_error_flag(message);
    if (err)
     {
      OomphLibWarning(message.str(),
                      "HypreSolver::hypre_solver_setup()",
                      OOMPH_EXCEPTION_LOCATION);
     }
   }

  // output times
  if (Output_time)
   {
    if ( (Hypre_method>=CG) && (Hypre_method<=BiCGStab) )
     {
      if ((Internal_preconditioner >= BoomerAMG)
          && (Internal_preconditioner <= ParaSails))
       {
        oomph_info << "Time for internal preconditioner setup [s] : "
                   << preconditioner_setup_time << std::endl;
       }
     }
    oomph_info << "Time for HYPRE solver setup [s]            : "
               << solver_setup_time
               << "s" << std::endl;
   }

 }



//===================================================================
/// Helper function performs a solve if solver data has been set
/// up using hypre_solver_setup(...).
//====================================================================
 void HypreInterface::hypre_solve(const Vector<double> &rhs,
                               Vector<double> &solution)
 {
  // Record time
#ifdef OOMPH_HAS_MPI   
  double t_start = MPI_Wtime();
  double t_end = 0;
#else
  clock_t t_start = clock();
  clock_t t_end = 0;;
#endif

  // Set up hypre vectors
  // --------------------

  // Hypre vector for rhs values
  HYPRE_IJVector rhs_ij;
  HYPRE_ParVector rhs_par;

  // Hypre vector for solution
  HYPRE_IJVector solution_ij;
  HYPRE_ParVector solution_par;
  
  // find length of rhs vectors
#ifdef OOMPH_HAS_MPI   
  const int vec_length = Last_row-First_row+1;
#else
  const int vec_length = int(Nrow);
#endif

  // storage for vector indices required to make Hypre vector
  Vector<int> vec_indices(vec_length);

  // pointer to the array of rhs values
  const double* rhs_values=0; 

  // set up rhs_values and vec_indices
#ifdef OOMPH_HAS_MPI 
  for (int i=0; i<vec_length;i++)
   {
    vec_indices[i] = First_row+i;
   }
  if (Using_distributed_rhs)
   {
    rhs_values = &rhs[0];
   }
  else
   {
    rhs_values = &rhs[First_row];
   }

  HypreHelpers::create_HYPRE_Vector(First_row,
                                    Last_row,
                                    rhs_ij,
                                    rhs_par,
                                    rhs_values,
                                    &vec_indices[0]);
  
  HypreHelpers::create_HYPRE_Vector(First_row,
                                    Last_row,
                                    solution_ij,
                                    solution_par);
#else
  rhs_values = &rhs[0];
  for (int i=0; i<vec_length;i++)
   {
    vec_indices[i] = i;
   }
  HypreHelpers::create_HYPRE_Vector(vec_length,
                                    rhs_ij,
                                    rhs_par,
                                    rhs_values,
                                    &vec_indices[0]);
  
  HypreHelpers::create_HYPRE_Vector(vec_length,
                                    solution_ij,
                                    solution_par);
#endif
  
  // Record time
  double vector_setup_time=0;
  if (Output_time)
   {
#ifdef OOMPH_HAS_MPI   
    t_end = MPI_Wtime();
    vector_setup_time = t_end-t_start; 
#else
    t_end = clock();
    vector_setup_time = double(t_end-t_start)/CLOCKS_PER_SEC;
#endif
   }
  
  // check error flag
  if (Hypre_error_messages)
   {
    std::ostringstream message;
    int err = HypreHelpers::check_HYPRE_error_flag(message);
    if (err)
     {
      OomphLibWarning(message.str(),
                      "HypreSolver::hypre_solve()",
                      OOMPH_EXCEPTION_LOCATION);
     }
   }
  
  // solve
  // -----
  t_start = t_end;
  
  // for solver stats
  int iterations=0;
  double norm=0;
  
  if (Existing_solver==BoomerAMG)
   {
    HYPRE_BoomerAMGSolve(Solver, Matrix_par, rhs_par, solution_par);
    HYPRE_BoomerAMGGetNumIterations(Solver, &iterations);
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(Solver, &norm);
   }
  else if (Existing_solver==CG)
   {
    HYPRE_PCGSolve(Solver,
                   (HYPRE_Matrix) Matrix_par, 
                   (HYPRE_Vector) rhs_par, 
                   (HYPRE_Vector) solution_par);
    HYPRE_PCGGetNumIterations(Solver, &iterations);
    HYPRE_PCGGetFinalRelativeResidualNorm(Solver, &norm);
   }
  else if (Existing_solver==GMRES)
   {
    HYPRE_GMRESSolve(Solver, 
                     (HYPRE_Matrix) Matrix_par, 
                     (HYPRE_Vector) rhs_par, 
                     (HYPRE_Vector) solution_par);
    HYPRE_GMRESGetNumIterations(Solver, &iterations);
    HYPRE_GMRESGetFinalRelativeResidualNorm(Solver, &norm);
   }
  else if (Existing_solver==BiCGStab)
   {
    HYPRE_BiCGSTABSolve(Solver,
                        (HYPRE_Matrix) Matrix_par, 
                        (HYPRE_Vector) rhs_par, 
                        (HYPRE_Vector) solution_par);
    HYPRE_BiCGSTABGetNumIterations(Solver, &iterations);
    HYPRE_BiCGSTABGetFinalRelativeResidualNorm(Solver, &norm);
   }
  else if (Existing_solver==Euclid)
   {
    HYPRE_EuclidSolve(Solver, Matrix_par, rhs_par, solution_par);
   }
  else if (Existing_solver==ParaSails)
   {
    HYPRE_ParaSailsSolve(Solver, Matrix_par, rhs_par, solution_par);
   }

  // output any error message
  if (Hypre_error_messages)
   {
    std::ostringstream message;
    int err = HypreHelpers::check_HYPRE_error_flag(message);
    if (err)
     {
      OomphLibWarning(message.str(),
                      "HypreSolver::hypre_solve()",
                      OOMPH_EXCEPTION_LOCATION);
     }
   }

  // Copy result to solution
#ifdef OOMPH_HAS_MPI   
  // If we're returning the distributed vector or there is only one
  // processor simply get the solution vector
  if ((Returning_distributed_solution) || (MPI_Helpers::Nproc==1))
   {
    solution.resize(vec_length);
    HYPRE_IJVectorGetValues(solution_ij,
                            vec_length,
                            &vec_indices[0],
                            &solution[0]);
   }
  else
   {
    // Vector to hold local solution values
    Vector<double> sol_vals_local(vec_length);
    
    // Get the solution
    HYPRE_IJVectorGetValues(solution_ij,
                            vec_length,
                            &vec_indices[0],
                            &sol_vals_local[0]);
  
 
    // Gather the results - first set up displacements for received 
    //local solution values
    Vector<int> displacements(MPI_Helpers::Nproc);
    displacements[0]=0;
    for (int p = 0; p<MPI_Helpers::Nproc-1; p++)
     {
      displacements[p+1]=displacements[p]+Nrows_local[p];
     }
    
    // gather the solution values
    solution.resize(Nrow);
    MPI_Allgatherv(&sol_vals_local[0],
                   vec_length,
                   MPI_DOUBLE,
                   &solution[0],
                   &Nrows_local[0],
                   &displacements[0],
                   MPI_DOUBLE,
                   MPI_COMM_WORLD);
   }
#else
  solution.resize(Nrow);
  HYPRE_IJVectorGetValues(solution_ij, 
                          vec_length, 
                          &vec_indices[0], 
                          &solution[0]);
#endif
  
  // output any error message
  if (Hypre_error_messages)
   {
    std::ostringstream message;
    int err = HypreHelpers::check_HYPRE_error_flag(message);
    if (err)
     {
      OomphLibWarning(message.str(),
                      "HypreSolver::hypre_solve()",
                      OOMPH_EXCEPTION_LOCATION);
     }
   }
  
  // deallocation
  HYPRE_IJVectorDestroy(solution_ij);
  HYPRE_IJVectorDestroy(rhs_ij);
  
  // Record time
  double solve_time=0;
  if (Output_time)
   {
#ifdef OOMPH_HAS_MPI   
    double t_end = MPI_Wtime();
    solve_time = t_end-t_start; 
#else
    clock_t t_end = clock();
    solve_time = double(t_end-t_start)/CLOCKS_PER_SEC;
#endif
   }
  
  // output timings and info
  if (Output_time)
   {
    oomph_info << "Time to generate HYPRE vectors [s]         : "
               << vector_setup_time << std::endl;
    oomph_info << "Time for HYPRE solve [s]                   : "
               << solve_time << std::endl;
   }
  
  // for iterative solvers output iterations and final norm
  if ((Hypre_method>=CG) && (Hypre_method<=BoomerAMG))
   {
    if (iterations>1)
     {
      if (Output_time) oomph_info << "Number of iterations         : "
                               << iterations << std::endl;
      if (Output_time) oomph_info << "Final Relative Residual Norm : " 
                               << norm << std::endl;
     }
   }
 }
 

//===================================================================
/// hypre_clean_up_memory() deletes any existing Hypre solver and
/// Hypre matrix
//====================================================================
 void HypreInterface::hypre_clean_up_memory()
 {
  // is there an existing solver
  if (Existing_solver != None)
   {
#ifdef OOMPH_HAS_MPI   
    double t_start = MPI_Wtime();
    Nrows_local.resize(0);
#else
    clock_t t_start = clock();
#endif

    // delete matrix
    HYPRE_IJMatrixDestroy(Matrix_ij);

    // delete solver
    if (Existing_solver==BoomerAMG)
     {
      HYPRE_BoomerAMGDestroy(Solver);
     }
    else if (Existing_solver==CG)
     {
      HYPRE_ParCSRPCGDestroy(Solver);
     }
    else if (Existing_solver==GMRES)
     {
      HYPRE_ParCSRGMRESDestroy(Solver);
     }
    else if (Existing_solver==BiCGStab)
     {
      HYPRE_ParCSRBiCGSTABDestroy(Solver);
     }
    else if (Existing_solver==Euclid)
     {
      HYPRE_EuclidDestroy(Solver);
     }
    else if (Existing_solver==ParaSails)
     {
      HYPRE_ParaSailsDestroy(Solver);
     }
    Existing_solver = None;

    // delete preconditioner
    if (Existing_preconditioner==BoomerAMG)
     {
      HYPRE_BoomerAMGDestroy(Preconditioner);
     }
    else if (Existing_preconditioner==Euclid)
     {
      HYPRE_EuclidDestroy(Preconditioner);
     }
    else if (Existing_preconditioner==ParaSails)
     {
      HYPRE_ParaSailsDestroy(Preconditioner);
     }
    Existing_preconditioner = None;

#ifdef OOMPH_HAS_MPI   
    double t_end = MPI_Wtime();
#else
    clock_t t_end = clock();
#endif

    // check error flag
    if (Hypre_error_messages)
     {
      std::ostringstream message;
      int err = HypreHelpers::check_HYPRE_error_flag(message);
      if (err)
       {
        OomphLibWarning(message.str(),
                        "HypreSolver::clean_up_memory()",
                        OOMPH_EXCEPTION_LOCATION);
       }
     }

    if (Output_time)
     {
      oomph_info << "Time to deallocate HYPRE solver data [s]   : "
#ifdef OOMPH_HAS_MPI   
                 << t_end-t_start
#else
                 << double(t_end-t_start)/CLOCKS_PER_SEC
#endif
                 << "\n";
     }
   }
 }


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// functions for HypreSolver class
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//===================================================================
/// Problem-based solve function to generate the Jacobian matrix and
/// residual vector and use HypreInterface::hypre_solver_setup(...) 
/// and HypreInterface::hypre_solve(...) to solve the linear system. 
/// This function will delete any existing data.
/// Note: the returned solution vector is NOT distributed, i.e. all
/// processors hold all values because this is what the Newton solver
/// requires.
//====================================================================
 void HypreSolver::solve(Problem* const &problem_pt,
                         Vector<double> &solution)
 {
  // Set Output_time flag for HypreInterface
  Output_time = Doc_time;

  // Delete any existing solver data
  clean_up_memory();

  // Set flag to allow deletion of the oomphlib Jacobian matrix
  // (we're in control)
  Delete_input_data = true;

  //  Get Jacobian and rhs
#ifdef OOMPH_HAS_MPI    
  // Record start time
  double t_start = MPI_Wtime();
  
  // Storage for Jacobian
  DistributedCRDoubleMatrix matrix;  
  
  // set flag so solver expects distributed rhs
  Using_distributed_rhs = true;
  
  // Set flat to return global solution vector 
  Returning_distributed_solution = false;

  // Storage for the Residual
  DistributedVector<double> residual;
#else
  // Record start time
  clock_t t_start = clock();
  
  // Storage for Jacobian
  CRDoubleMatrix matrix;

  // Storage for the Residual
  Vector<double> residual;
#endif
  
  //Get oomph-lib Jacobian matrix and residual vector
  problem_pt->get_jacobian(residual,matrix);
  
  // Output times
  if(Doc_time)
   {
    oomph_info << "Time to generate Jacobian and residual [s] : ";
#ifdef OOMPH_HAS_MPI   
    double t_end = MPI_Wtime();
    oomph_info << t_end-t_start << std::endl; 
#else
    clock_t t_end = clock();
    oomph_info << double(t_end-t_start)/CLOCKS_PER_SEC << std::endl;
#endif
   }
  
  // generate hypre matrix
  hypre_matrix_setup(&matrix);
  
  // call hypre_solver_setup to generate linear solver data
  hypre_solver_setup();
  
  // perform hypre_solve
#ifdef OOMPH_HAS_MPI
  hypre_solve(residual.vector(), solution);
#else
  hypre_solve(residual, solution);
#endif  

  // delete solver data if required
  if (!Enable_resolve)
   {
    clean_up_memory();
   }
 }
 
 
 
//===================================================================
/// Uses HypreInterface::hypre_solve(...) to solve the linear system 
/// for a CRDoubleMatrix or a DistributedCRDoubleMatrix. 
/// In the latter case, the rhs needs to be distributed too, 
/// i.e. the length of the rhs vector must be equal to the number of 
/// rows stored locally. 
/// Note: the returned solution vector is never distributed, i.e. all
/// processors hold all values
//====================================================================
 void HypreSolver::solve(DoubleMatrixBase* const &matrix_pt,
                         const Vector<double> &rhs,
                         Vector<double> &solution)
 {
#ifdef PARANOID
  // check the matrix is square
  if ( matrix_pt->nrow() != matrix_pt->ncol() )
   {
    std::ostringstream error_message;
    error_message << "HypreSolver require a square matrix. "
                  << "Matrix is " << matrix_pt->nrow()
                  << " by " << matrix_pt->ncol() << std::endl;
    throw OomphLibError(error_message.str(),
                        "HypreSolver::solve()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Set Output_time flag for HypreInterface
  Output_time = Doc_time;

  // Clean up existing solver data
  clean_up_memory();

  // Set flag to decide if oomphlib matrix can be deleted
  // (Recall that Delete_matrix defaults to false).
  Delete_input_data = Delete_matrix;

  // Flag to record if matrix cast was successful
  bool successful_matrix_cast = false;

  // Try cast to a CRDoubleMatrix
  CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);
  
  // If cast successful set things up for a serial solve
  if (cr_matrix_pt)
   {
#ifdef PARANOID
    // check the matrix and rhs are of consistent sizes
    if ( matrix_pt->nrow() != rhs.size() )
     {
      std::ostringstream error_message;
      error_message << "HypreSolver require rhs vector length equals number of"
                    << " rows/columns in a CRDoubleMatrix. "
                    << "Vector length is " << rhs.size()
                    << " and number of rows/columns in matrix is "
                    << matrix_pt->ncol() << std::endl;
      throw OomphLibError(error_message.str(),
                          "HypreSolver::solve()",
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif
    successful_matrix_cast = true;
    hypre_matrix_setup(cr_matrix_pt);
    
#ifdef OOMPH_HAS_MPI
    // Set flag for non-distributed rhs in the solver
    Using_distributed_rhs = false;
#endif
   }
  
#ifdef OOMPH_HAS_MPI
  // Set flat to return global solution vector 
  Returning_distributed_solution = false;

  // If required try casting to DistributedCRDoubleMatrix
  DistributedCRDoubleMatrix* dist_cr_matrix_pt=0; 
  if (!successful_matrix_cast)
   {
    dist_cr_matrix_pt = dynamic_cast<DistributedCRDoubleMatrix*>(matrix_pt);
    
    // If cast successful set things up for a parallel solve
    if (dist_cr_matrix_pt)
     {
#ifdef PARANOID
      // check the matrix and rhs are of consistent sizes
      if ( unsigned(dist_cr_matrix_pt->nrow_local()) != rhs.size() )
       {
        std::ostringstream error_message;
        error_message
         << "HypreSolver require rhs vector length equals number of"
         << " local rows in a DistributedCRDoubleMatrix. "
         << "Vector length is " << rhs.size()
         << " and number of local rows in matrix is "
         << dist_cr_matrix_pt->nrow_local() << std::endl;
        throw OomphLibError(error_message.str(),
                            "HypreSolver::solve()",
                            OOMPH_EXCEPTION_LOCATION);
       }
#endif
      successful_matrix_cast = true;
      hypre_matrix_setup(dist_cr_matrix_pt);
      
      // Set flag for non-distributed rhs in the solver
      Using_distributed_rhs = true;
     }
   }
#endif
  
#ifdef PARANOID
  // Check casts were successful
  if (!successful_matrix_cast)
   {
    std::ostringstream error_message;
    error_message << "HypreSolver only work with "
#ifdef OOMPH_HAS_MPI
                  << "CRDoubleMatrix and DistributedCRDoubleMatrix matrices" 
#else
                  << "CRDoubleMatrix matrices"
#endif
                  << std::endl;
    throw OomphLibError(error_message.str(),
                     	"HypreSolver::solve()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // call hypre_setup to generate Hypre matrix and linear solver data
  hypre_solver_setup();
  
  // perform hypre_solve
  hypre_solve(rhs, solution);
  
  // delete solver data if required
  if (!Enable_resolve)
   {
    clean_up_memory();
   }
 }
 
 
//===================================================================
/// Resolve performs a linear solve using current solver data (if
/// such data exists).
//====================================================================
 void HypreSolver::resolve(const Vector<double> &rhs,
                           Vector<double> &solution)
 {
#ifdef PARANOID
  // check solver data exists
  if (existing_solver()==None)
   {
    std::ostringstream error_message;
    error_message << "resolve(...) requires that solver data has been "
                  << "set up by a previous call to solve(...) after "
                  << "a call to enable_resolve()" << std::endl;
    throw OomphLibError(error_message.str(),
                     	  "HypreSolver::resolve()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Set Output_time flag for HypreInterface
  Output_time = Doc_time;

  // solve
  hypre_solve(rhs, solution);

  // Note: never delete solver data as the preconditioner is typically
  // called repeatedly.
 }


//===================================================================
/// clean_up_memory() deletes any existing solver data
//====================================================================
 void HypreSolver::clean_up_memory()
 {
  hypre_clean_up_memory();
 }



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// functions for HyprePreconditioner class
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//=============================================================================
/// An interface to allow HypreSolver to be used as a Preconditioner
/// for the oomph-lib IterativeLinearSolver class.
/// Matrix has to be of type CRDoubleMatrix or DistributedCRDoubleMatrix.
//=============================================================================
 void HyprePreconditioner::setup(Problem* problem_pt, 
                                 DoubleMatrixBase* matrix_pt)
 {
  // Set Output_time flag for HypreInterface
  Output_time = Doc_time;
  
#ifdef PARANOID
  // check the matrix is square
  if ( matrix_pt->nrow() != matrix_pt->ncol() )
   {
    std::ostringstream error_message;
    error_message << "HyprePreconditioner require a square matrix. "
                  << "Matrix is " << matrix_pt->nrow()
                  << " by " << matrix_pt->ncol() << std::endl;
    throw OomphLibError(error_message.str(),
                        "HyprePreconditioner::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // clean up any previous solver data
  clean_up_memory();

  // set flag to decide if oomphlib matrix can be deleted
  // (Recall that Delete_matrix defaults to false).
  Delete_input_data = Delete_matrix;
  
  // Try casting to a CRDoubleMatrix
  bool successful_matrix_cast = false;
  CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);
  
  // If cast successful set things up for a serial solve
  if (cr_matrix_pt)
   {
    successful_matrix_cast = true;
    hypre_matrix_setup(cr_matrix_pt);
    
#ifdef OOMPH_HAS_MPI
    // Set flag for non-distributed rhs and solution in the solver
    Using_distributed_rhs = false;
    Returning_distributed_solution = false;
#endif
   }
  
#ifdef OOMPH_HAS_MPI
  // If required try casting to DistributedCRDoubleMatrix
  DistributedCRDoubleMatrix* dist_cr_matrix_pt=0; 
  if (!successful_matrix_cast)
   {
    dist_cr_matrix_pt = dynamic_cast<DistributedCRDoubleMatrix*>(matrix_pt);
    
    // If cast successful set things up for a parallel solve
    if (dist_cr_matrix_pt)
     {
      successful_matrix_cast = true;
      hypre_matrix_setup(dist_cr_matrix_pt);
      
      // Set flag for non-distributed rhs and soltion in the solver
      Using_distributed_rhs = true;
      Returning_distributed_solution = true;
     }
   }
#endif
  
#ifdef PARANOID
  // Check casts were successful
  if (!successful_matrix_cast)
   {
    std::ostringstream error_message;
    error_message << "HyprePreconditioner only work with "
#ifdef OOMPH_HAS_MPI
                  << "CRDoubleMatrix and DistributedCRDoubleMatrix matrices" 
#else
                  << "CRDoubleMatrix matrices"
#endif
                  << std::endl;
    throw OomphLibError(error_message.str(),
                     	"HyprePreconditioner::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
 
  // call hypre_solver_setup
  hypre_solver_setup();

 }

 
//===================================================================
/// Preconditioner_solve uses a hypre solver to precondition vector r
//====================================================================
 void HyprePreconditioner::preconditioner_solve(const Vector<double> &r,
                                        Vector<double> &z)
 {
#ifdef PARANOID
  // check solver data exists
  if (existing_solver()==None)
   {
    std::ostringstream error_message;
    error_message << "preconditioner_solve(...) requires that data has "
                  << "been set up using the function setup(...)" << std::endl;
    throw OomphLibError(error_message.str(),
                     	"HyprePreconditioner::preconditioner_solve()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Switch off any timings for the solve
  Output_time = false;
  
  // perform hypre_solve
  hypre_solve(r,z);

 }

 

//===================================================================
/// clean_up_memory() deletes any existing Hypre solver and
/// Hypre matrix
//====================================================================
 void HyprePreconditioner::clean_up_memory()
 {
  hypre_clean_up_memory();
 }

}
