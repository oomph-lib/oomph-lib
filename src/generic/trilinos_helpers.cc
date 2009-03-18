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
#include "trilinos_helpers.h"

namespace oomph
{
namespace TrilinosHelpers
{

 //============================================================================
 /// \short creates a distributed Epetra_Vector that has the same contents as
 /// oomph_v with distribution row_map_pt. \n
 /// NOTE 1. the Epetra_Map and the LinearAlgebraDistribution object in the
 /// vector must descibe the same distribution
 /// NOTE 2. if the bool view is false (default) then the values in the 
 /// oomph-lib vector (oomph_v) will be copied into the epetra vector, or if 
 /// it is true then the epetra vector will 'view' (i.e. point to) the 
 /// contents of the oomph-lib vector
 //============================================================================
 void create_epetra_vector(const DoubleVector& oomph_v,
                           const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt,
                           bool view)
 {
#ifdef PARANOID
  // check the the oomph lib vector is setup
  if (!oomph_v.distribution_pt()->setup())
   {
    std::ostringstream error_message;
    error_message << "The oomph-lib vector (oomph_v) must be setup.";
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpers::create_epetra_vector()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#ifdef OOMPH_HAS_MPI
  // if the oomph-lib vector is distributed, then check it has the same
  // distribution as the Epetra_Map
  if (oomph_v.distributed() && 
      !compare_maps(oomph_v.distribution_pt(),*row_map_pt))
   {
    std::ostringstream error_message;
    error_message << "The oomph-lib distributed vector (oomph_v) and the "
                  << "epetra map must describe the same distribution.";
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpers::create_epetra_vector()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  // if the oomph-lib vector is no distribited, the check is has the same 
  // number of global rows
  if (static_cast<int>(oomph_v.nrow()) != row_map_pt->NumGlobalElements())
   {
    std::ostringstream error_message;
    error_message
     << "The epetra map and the oomph-lib matrix have different numbers "
     << "of global rows.\n"
     << "NumGlobalElements() for the Epetra_Map: " 
     << row_map_pt->NumGlobalElements() << "\n" 
     << "nrow() for DoubleVector: " << oomph_v.nrow();
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpers::create_epetra_vector()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // first local row of oomph vector to be inserted into the Epetra_Vector
  unsigned offset = 0;
  // NOTE assumption that the first global row in the vector is the 
  // first_row. This is true of oomph-lib vectors, but not necessarily 
  // for epetra vectors
  if (!oomph_v.distributed() && 
      oomph_v.distribution_pt()->communicator_pt()->nproc() > 1)
   {
    int* global_rows = row_map_pt->MyGlobalElements();
    offset = static_cast<unsigned>(global_rows[0]);
   }

  //
  double* v_pt = oomph_v.values_pt();
  if (view)
   {
    epetra_v_pt = new Epetra_Vector(View,*row_map_pt,v_pt+offset);
   }
  else
   {
    epetra_v_pt = new Epetra_Vector(Copy,*row_map_pt,v_pt+offset);
   }
 }

 //============================================================================
 /// creates an empty trilinos vector with Epetra_Map row_map_pt
 //============================================================================
 void create_epetra_vector(const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt)
 {

  // create the Epetra_Vector
  epetra_v_pt = new Epetra_Vector(*row_map_pt,true);
 }

 //============================================================================
 /// \short Helper function to copy the contents of a Trilinos vector to an
 /// oomph-lib distributed vector. The distribution of the two vectors must
 /// be identical
 //============================================================================
 void copy_to_oomphlib_vector(const Epetra_Vector* epetra_v_pt,
                              DoubleVector& oomph_v)
 {
#ifdef PARANOID
  // check the the oomph lib vector is setup
  if (!oomph_v.distribution_pt()->setup())
   {
    std::ostringstream error_message;
    error_message << "The oomph-lib vector (oomph_v) must be setup.";
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpers::copy_to_oomphlib_vector()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // if the oomph-lib vector is distributed
  if (oomph_v.distributed())
   {
    // extract values from epetra_v_pt
    double* v_values;
    epetra_v_pt->ExtractView(&v_values);
    
    // copy the values
    unsigned nrow_local = oomph_v.nrow_local();
    for (unsigned i = 0; i < nrow_local; i++)
     {
      oomph_v[i] = v_values[i];
     }
   }

  // else teh oomph-lib vector is not distributed
  else
   {

    // number of global rows
    unsigned nrow = epetra_v_pt->GlobalLength();

    // values to be copied into the oomph-lib vector
    double* values = oomph_v.values_pt();
    
    // get the values vector
#ifdef OOMPH_HAS_MPI
    int nproc = epetra_v_pt->Map().Comm().NumProc();
    if (nproc == 1)
     {
      epetra_v_pt->ExtractView(&values);
     }
    else
     {
      // get the local values
      double* local_values;
      epetra_v_pt->ExtractView(&local_values);

      // my rank
      int my_rank =epetra_v_pt->Map().Comm().MyPID();
      
      // number of local rows
      Vector<int> nrow_local(nproc);
      nrow_local[my_rank] = epetra_v_pt->MyLength();
      

      // gather the First_row vector
      MPI_Allgather(&nrow_local[my_rank],1,MPI_INT,&nrow_local[0],1,MPI_INT,
                    oomph_v.distribution_pt()->communicator_pt()->mpi_comm());
      
      // number of local rows
      Vector<int> first_row(nproc);
      int* global_rows = epetra_v_pt->Map().MyGlobalElements();
      first_row[my_rank] = global_rows[0];
      
      // gather the First_row vector
      MPI_Allgather(&first_row[my_rank],1,MPI_INT,&first_row[0],1,MPI_INT,
                    oomph_v.distribution_pt()->communicator_pt()->mpi_comm());

      // gather the local solution values
      values = new double[nrow];
      MPI_Allgatherv(local_values,nrow_local[my_rank],MPI_DOUBLE,values,
                     &nrow_local[0],&first_row[0],MPI_DOUBLE,
                     oomph_v.distribution_pt()->communicator_pt()->mpi_comm());
    }
#else
    epetra_v_pt->ExtractView(&values);
#endif
    double* oomph_v_pt = oomph_v.values_pt();
    for (unsigned i=0; i<nrow; i++)
     {
      oomph_v_pt[i] = values[i];
     }
   }
 }

//=============================================================================
/// \short Helper function to create a distributed Epetra_CrsMatrix from a 
/// from a CRDoubleMatrix
/// \b NOTE 1. This function constructs an Epetra_CrsMatrix using new,
/// "delete trilinos_matrix_pt;" is NOT called.
//=============================================================================
void create_epetra_matrix(DoubleMatrixBase* matrix_pt,
                          const Epetra_Map* epetra_range_map_pt,
                          const Epetra_Map* epetra_domain_map_pt,
                          const Epetra_Map* epetra_col_map_pt,
                          Epetra_CrsMatrix* &epetra_matrix_pt,
                          bool view)
{

 // first try DistributedCRDoubleMatrix
 CRDoubleMatrix* cast_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);

#ifdef PARANOID
 // if the cast failed then throw error
 if (cast_matrix_pt == 0)
  {
   std::ostringstream error_message;
   error_message << "An create_epetra_matrix(...) is only compatiible with "
                 << "CRDoubleMatrix" << std::endl;
   throw OomphLibError(error_message.str(),
                       "TrilinosHelpers::create_epetra_matrix()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check the matrix is setup
 if (!cast_matrix_pt->built())
  {
    std::ostringstream error_message;
    error_message << "The oomph-lib matrix must be built.";
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpers::create_epetra_matrix()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

 // if this is a MPI build then attempt to build Trilinos Matrix from
 // a CRDoubleMatrix
#ifdef OOMPH_HAS_MPI
 if (cast_matrix_pt->distributed())
  {

   // paranoid check that the distribution of the epetra row map and the 
   // oomph-lib matrix are the same
#ifdef PARANOID
   compare_maps(cast_matrix_pt->distribution_pt(),
                *epetra_range_map_pt);
#endif
   // paranoid check that the distribution of the map and the vector are the
  // same
#ifdef PARANOID
  if (!compare_maps(cast_matrix_pt->distribution_pt(),
                    *epetra_range_map_pt))
  {
    std::ostringstream error_message;
    error_message << "The oomph-lib distributed vector (oomph_v) and the "
                  << "epetra map must describe the same distribution.";
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpers::create_epetra_matrix()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

   // get pointers to the matrix values, column indices etc
   int* column = cast_matrix_pt->column_index();
   double* value = cast_matrix_pt->value();
   int* row_start = cast_matrix_pt->row_start();
   
   // get my nrow_local and first_row
   unsigned nrow_local = cast_matrix_pt->nrow_local();
   unsigned first_row = cast_matrix_pt->first_row();
   
   // store the number of non zero entries per row
   int* nnz_per_row = new int[nrow_local];
   for (unsigned row=0;row<nrow_local;row++)
    {
     nnz_per_row[row] = row_start[row+1] - row_start[row];
    }
   
   // construct the new Epetra matrix
   if (view)
    {
     epetra_matrix_pt = new Epetra_CrsMatrix(View,*epetra_range_map_pt,
                                             *epetra_col_map_pt,
                                             nnz_per_row,true);
    }
   else
    {
     epetra_matrix_pt = new Epetra_CrsMatrix(Copy,*epetra_range_map_pt,
                                             nnz_per_row,true);
    }

   // insert the values
   for (unsigned row=0; row<nrow_local; row++)
    {
     // get pointer to this row in values/columns
     int ptr = row_start[row];
#ifdef PARANOID
     int err = epetra_matrix_pt->InsertGlobalValues(first_row+row,
                                                    nnz_per_row[row],
                                                    value+ptr,
                                                    column+ptr);
     if (err != 0)
      {
       std::ostringstream error_message;
       error_message 
        << "Epetra Matrix Insert Global Values : epetra_error_flag = " 
        << err;
       throw OomphLibError(error_message.str(),
                           "TrilinosHelpers::create_epetra_matrix(...)",
                           OOMPH_EXCEPTION_LOCATION);       
      }
#else
     epetra_matrix_pt->InsertGlobalValues(first_row+row,
                                          nnz_per_row[row],
                                          value+ptr,
                                          column+ptr);
#endif
    }

   // complete the build of the trilinos matrix
#ifdef PARANOID
   int err=0;
   err = epetra_matrix_pt->FillComplete(*epetra_domain_map_pt,
                                        *epetra_range_map_pt);
   if (err != 0)
    {
     std::ostringstream error_message;
     error_message 
      << "Epetra Matrix Fill Complete Error : epetra_error_flag = " 
      << err;
     throw OomphLibError(error_message.str(),
                         "TrilinosHelpers::create_epetra_matrix(...)",
                         OOMPH_EXCEPTION_LOCATION);
    }
#else
   epetra_matrix_pt->FillComplete(*epetra_domain_map_pt,*epetra_range_map_pt);
#endif

   // tidy up memory
   delete[] nnz_per_row;
  }

 // if the matrix is not distributed 
 else if (!cast_matrix_pt->distributed())
  {

   // get pointers to the matrix values, column indices etc
   int* column = cast_matrix_pt->column_index();
   double* value = cast_matrix_pt->value();
   int* row_start = cast_matrix_pt->row_start();
   
   // get my nrow_local and first_row
   unsigned nrow_local = epetra_range_map_pt->NumMyElements();
   int* global_rows = epetra_range_map_pt->MyGlobalElements();
   unsigned first_row = static_cast<unsigned>(global_rows[0]);

   // store the number of non zero entries per row
   int* nnz_per_row = new int[nrow_local];
   for (unsigned row=0;row<nrow_local;row++)
    {
     nnz_per_row[row] = row_start[first_row + row+1] 
      - row_start[first_row + row];
    }
   
   // construct the new Epetra matrix
   if (view)
    {
     epetra_matrix_pt = new Epetra_CrsMatrix(View,*epetra_range_map_pt,
                                             *epetra_col_map_pt,
                                             nnz_per_row,true);
    }
   else
    {
     epetra_matrix_pt = new Epetra_CrsMatrix(Copy,*epetra_range_map_pt,
                                             nnz_per_row,true);
    }

   // insert the values
   for (unsigned row=0; row< nrow_local; row++)
    {

     // get pointer to this row in values/columns
     int ptr = row_start[row+first_row];
     epetra_matrix_pt->InsertGlobalValues(first_row + row,
                                          nnz_per_row[row],
                                          value+ptr,
                                          column+ptr);
    }
   
   // tidy up memory
   delete[] nnz_per_row;
   
   // complete the build of the trilinos matrix
#ifdef PARANOID
   int err=0;
   err = epetra_matrix_pt->FillComplete(*epetra_domain_map_pt,
                                        *epetra_range_map_pt);
   if (err != 0)
    {
     std::ostringstream error_message;
     error_message 
      << "Epetra Matrix Fill Complete Error : epetra_error_flag = " 
      << err;
     throw OomphLibError(error_message.str(),
                         "TrilinosHelpers::create_epetra_matrix(...)",
                         OOMPH_EXCEPTION_LOCATION);
    }
#else
   epetra_matrix_pt->FillComplete(*epetra_domain_map_pt,*epetra_range_map_pt);
#endif

 }
#else

 // find number of rows and columns
 unsigned n_rows = cast_matrix_pt->nrow();

 // get pointers to the matrix values, column indices etc
 int* columns = cast_matrix_pt->column_index();
 double* values = cast_matrix_pt->value();
 int* row_starts = cast_matrix_pt->row_start();

 // store the number of non zero entries per row
 int* nnz_per_row = new int[n_rows];
 for (unsigned row=0;row<n_rows;row++)
  {
   nnz_per_row[row] = row_starts[row+1] - row_starts[row];
  }

 // construct the new Epetra matrix
 if (view)
  {
   epetra_matrix_pt = new Epetra_CrsMatrix(View,*epetra_range_map_pt,
                                           *epetra_col_map_pt,
                                           nnz_per_row);
  }
 else
  {
   epetra_matrix_pt = new Epetra_CrsMatrix(Copy,*epetra_range_map_pt,
                                           *epetra_col_map_pt,
                                           nnz_per_row);
  }

 // insert the values
 for (unsigned row=0;row<n_rows;row++)
  {
   // get pointer to this row in values/columns
   int ptr = row_starts[row];
   epetra_matrix_pt->InsertGlobalValues(row,
                                        nnz_per_row[row],
                                        values+ptr,
                                        columns+ptr);
  }
 delete[] nnz_per_row;

 // complete the build of the trilinos matrix
#ifdef PARANOID
 int err=0;
 err = epetra_matrix_pt->FillComplete(*epetra_domain_map_pt,
                                      *epetra_range_map_pt);
 if (err != 0)
  {
   std::ostringstream error_message;
     error_message
      << "Epetra Matrix Fill Complete Error : epetra_error_flag = "
      << err;
     throw OomphLibError(error_message.str(),
                         "TrilinosHelpers::create_epetra_matrix(...)",
                         OOMPH_EXCEPTION_LOCATION);
  }
#else
 epetra_matrix_pt->FillComplete(*epetra_domain_map_pt,*epetra_range_map_pt);
#endif
#endif
  }


//=============================================================================
/// Helper function to create a distributed Epetra_CrsMatrix from a 
/// from a DistributedCRDoubleMatrix or a CRDoubleMatrix
/// \b NOTE 1. This function constructs an Epetra_CrsMatrix using new,
/// "delete trilinos_matrix_pt;" is NOT called.
/// \b NOTE 2. Specialisation for SQUARE matrices.
//=============================================================================
void create_epetra_matrix(DoubleMatrixBase* oomph_matrix,
                          const Epetra_Map* epetra_range_map_pt,
                          const Epetra_Map* epetra_col_map_pt,
                          Epetra_CrsMatrix* &epetra_matrix_pt)
{
 create_epetra_matrix(oomph_matrix,epetra_range_map_pt,epetra_range_map_pt,
                      epetra_col_map_pt,epetra_matrix_pt);
}


//=============================================================================
/// Function to perform a matrix-vector multiplication on a 
/// distributed matrix and a distributed vector using Trilinos functionality.
/// \n
/// \b NOTE 1. the matrix (matrix) and the vectors (x and soln) must have the 
/// same communicator.
/// \b NOTE 2. The vector (soln) will be returned with the same distribution
/// as the matrix, unless a distribution is predefined in the solution
/// vector in which case the vector will be returned with that distribution
//=============================================================================
void multiply(CRDoubleMatrix &oomph_matrix,
              const DoubleVector &oomph_x,
              DoubleVector &oomph_soln)
{
#ifdef PARANOID
 // check that this matrix is built
 if (!oomph_matrix.built())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "This matrix has not been built";
   throw OomphLibError(error_message_stream.str(),
                       "TrilinosHelpers::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // check that the distribution of the matrix and the soln are the same
 if (oomph_soln.distribution_setup())
  {
   if (!(*oomph_matrix.distribution_pt() == *oomph_soln.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The soln vector and this matrix must have the same distribution.";
     throw OomphLibError(error_message_stream.str(),
                       "TrilinosHelpers::multiply()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
 
 // check that the distribution of the oomph-lib vector x is setup
 if (!oomph_x.distribution_pt()->setup())
  {
   std::ostringstream error_message_stream;
     error_message_stream 
      << "The x vector must be setup";
     throw OomphLibError(error_message_stream.str(),
                       "TrilinosHelpers::multiply()",
                         OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // setup the distribution
 if (!oomph_soln.distribution_pt()->setup())
  {
   oomph_soln.rebuild(oomph_matrix.distribution_pt());
  }

 // create the communicator
#ifdef OOMPH_HAS_MPI
 Epetra_MpiComm* comm_pt = new Epetra_MpiComm(
  oomph_matrix.distribution_pt()->communicator_pt()->mpi_comm());
#else
 Epetra_SerialComm* comm_pt = new Epetra_SerialComm();
#endif

 // create the maps corresponding the distribution of matrix1 and matrix2

 // matrix range map
 Epetra_Map* matrix_range_map_pt;

#ifdef OOMPH_HAS_MPI
 int* global_rows;
 create_epetra_map(oomph_matrix.distribution_pt(),comm_pt,matrix_range_map_pt,
                   global_rows);
#else
 create_epetra_map(oomph_matrix.distribution_pt(),comm_pt,matrix_range_map_pt);
#endif

 // matrix domain map
 Epetra_Map* matrix_domain_map_pt;

#ifdef OOMPH_HAS_MPI
 int* global_cols;
 create_epetra_map(oomph_x.distribution_pt(),comm_pt,matrix_domain_map_pt,
                   global_cols); 
#else
 create_epetra_map(oomph_x.distribution_pt(),comm_pt,matrix_domain_map_pt); 
#endif

 // matrix col map
 Epetra_Map* matrix_col_map_pt = new Epetra_Map(oomph_matrix.ncol(),
                                                oomph_matrix.ncol(),
                                                0,*comm_pt);
 
 // convert matrix1 to epetra matrix
 Epetra_CrsMatrix* epetra_matrix_pt;
 create_epetra_matrix(&oomph_matrix,matrix_range_map_pt,matrix_domain_map_pt,
                      matrix_col_map_pt,epetra_matrix_pt,false);

 // convert x to Trilinos vector
 Epetra_Vector* epetra_x_pt;
 create_epetra_vector(oomph_x,matrix_domain_map_pt,epetra_x_pt);

 // create Trilinos vector for soln ( 'viewing' the contents of the oomph-lib
 // matrix)
 Epetra_Vector* epetra_soln_pt;
 create_epetra_vector(oomph_soln,matrix_range_map_pt,epetra_soln_pt);

 // do the multiply
 int epetra_error_flag = epetra_matrix_pt->Multiply(false,*epetra_x_pt,
                                                    *epetra_soln_pt);

 // return the solution
 copy_to_oomphlib_vector(epetra_soln_pt,oomph_soln);

 // throw error if there is an epetra error
#if PARANOID
 if (epetra_error_flag != 0)
  {
   std::ostringstream error_message;
   error_message 
    << "Epetra Matrix Vector Multiply Error : epetra_error_flag = " 
    << epetra_error_flag;
   throw OomphLibError(error_message.str(),
                       "TrilinosHelpersMPI::multiply(...)",
                       OOMPH_EXCEPTION_LOCATION);
   }
#endif

 // do something with the error flag to keep the compiler from
 // complaining in non-paranoid mode
 epetra_error_flag=0;

 // clean up
 delete matrix_range_map_pt;
 delete matrix_domain_map_pt;
 delete matrix_col_map_pt;
 delete epetra_matrix_pt;
 delete epetra_x_pt;
 delete epetra_soln_pt;
 delete comm_pt;
#ifdef OOMPH_HAS_MPI
 delete[] global_rows;
 delete[] global_cols;
#endif
}

//=============================================================================
/// \short Function to perform a matrix-matrix multiplication on distributed 
/// matrices by using Trilinos functionality.\n
/// \b NOTE 1. There are two Trilinos matrix-matrix multiplication methods 
/// available, using either the EpetraExt::MatrixMatrix class (if use_ml == 
/// false) or using ML (Epetra_MatrixMult method)
/// \b NOTE 2. the solution matrix (matrix_soln) will be returned with the 
/// same distribution as matrix1
/// \b NOTE 3. All matrices must share the same communicator. 
//=============================================================================
void multiply(CRDoubleMatrix &matrix1,
              CRDoubleMatrix &matrix2,
              CRDoubleMatrix &matrix_soln,
              const bool& use_ml)
{


#ifdef PARANOID
 // check that matrix 1 is built
 if (!matrix1.built())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "This matrix matrix1 has not been built";
   throw OomphLibError(error_message_stream.str(),
                       "TrilinosHelpers::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check that matrix 2 is built
 if (!matrix2.built())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "This matrix matrix2 has not been built";
   throw OomphLibError(error_message_stream.str(),
                       "TrilinosHelpers::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check matrix dimensions are compatable
 if ( matrix1.ncol() != matrix2.nrow()  )
 {
  std::ostringstream error_message;
  error_message
   << "Matrix dimensions incompatible for matrix-matrix multiplication"
   << "ncol() for first matrix: " << matrix1.ncol()
   << "nrow() for second matrix: " << matrix2.nrow();
  throw OomphLibError(error_message.str(),
                      "TrilinosHelpers::multiply(...)",
                      OOMPH_EXCEPTION_LOCATION);
 }
 // check that the have the same communicator
 OomphCommunicator temp_comm(matrix1.distribution_pt()->communicator_pt());
 if (temp_comm != *matrix2.distribution_pt()->communicator_pt())
  {
   std::ostringstream error_message;
   error_message
    << "Matrix dimensions incompatible for matrix-matrix multiplication"
    << "ncol() for first matrix: " << matrix1.ncol()
    << "nrow() for second matrix: " << matrix2.nrow();
   throw OomphLibError(error_message.str(),
                       "TrilinosHelpers::multiply(...)",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check that the distribution of the matrix and the soln are the same
 if (matrix_soln.distribution_pt()->setup())
  {
   if (!(*matrix_soln.distribution_pt() == *matrix1.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The solution matrix and matrix1 must have the same distribution.";
     throw OomphLibError(error_message_stream.str(),
                       "TrilinosHelpers::multiply()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif

 // setup the distribution
 if (!matrix_soln.distribution_pt()->setup())
  {
   matrix_soln.rebuild(matrix1.distribution_pt());
  }

 // temporary fix
 // ML MM method only appears to work for square matrices
 // Should be investigated further.
 bool temp_use_ml = false;
 if ((*matrix1.distribution_pt() == *matrix2.distribution_pt()) &&
     (matrix1.ncol() == matrix2.ncol()))
  {
   temp_use_ml = use_ml;
  }

 // create the communicator
#ifdef OOMPH_HAS_MPI
 Epetra_MpiComm* comm_pt = new Epetra_MpiComm(
  matrix1.distribution_pt()->communicator_pt()->mpi_comm());
#else
 Epetra_SerialComm* comm_pt = new Epetra_SerialComm();
#endif
 
 // create the maps corresponding the distribution of matrix1 and matrix2

 // matrix 1
 Epetra_Map* matrix1_range_map_pt;

#ifdef OOMPH_HAS_MPI
 int* matrix1_global_rows;
 create_epetra_map(matrix1.distribution_pt(),comm_pt,matrix1_range_map_pt,
                   matrix1_global_rows);
#else
 create_epetra_map(matrix1.distribution_pt(),comm_pt,matrix1_range_map_pt);
#endif

 // matrix 2
 Epetra_Map* matrix2_range_map_pt;

#ifdef OOMPH_HAS_MPI
 int* matrix2_global_rows;
 create_epetra_map(matrix2.distribution_pt(),comm_pt,matrix2_range_map_pt,
                   matrix2_global_rows);
#else
 create_epetra_map(matrix2.distribution_pt(),comm_pt,matrix2_range_map_pt);
#endif

 // matrix 2 column distribution
 Epetra_Map* matrix2_domain_map_pt;
 LinearAlgebraDistribution 
  matrix2_col_distribution(matrix2.distribution_pt()->communicator_pt(),
                           matrix2.ncol(),true); 

#ifdef OOMPH_HAS_MPI
 int* matrix2_global_cols;
 create_epetra_map(&matrix2_col_distribution,comm_pt,matrix2_domain_map_pt,
                   matrix2_global_cols);
#else
 create_epetra_map(&matrix2_col_distribution,comm_pt,matrix2_domain_map_pt);
#endif
 // create matrix 1
 Epetra_Map* matrix1_col_map_pt = new Epetra_Map(matrix1.ncol(),
                                                 matrix1.ncol(),
                                                 0,*comm_pt);
 Epetra_CrsMatrix* epetra_matrix1_pt;
 create_epetra_matrix(&matrix1,matrix1_range_map_pt,matrix2_range_map_pt,
                      matrix1_col_map_pt,epetra_matrix1_pt,false);
 
 // create matrix 2
 Epetra_Map* matrix2_col_map_pt = new Epetra_Map(matrix2.ncol(),
                                                 matrix2.ncol(),
                                                 0,*comm_pt);
 Epetra_CrsMatrix* epetra_matrix2_pt;
 create_epetra_matrix(&matrix2,matrix2_range_map_pt,matrix2_domain_map_pt,
                      matrix2_col_map_pt,epetra_matrix2_pt,false);

 // create Trilinos matrix to hold solution - will have same map 
 // (and number of rows) as matrix1
 Epetra_CrsMatrix* solution_pt;

 // do the multiplication
 // ---------------------
 if (temp_use_ml)
 {
  // there is a problem using this function, many pages of
  // warning messages are issued....
  // "tmpresult->InsertGlobalValues returned 3"
  // and
  // "Result_epet->InsertGlobalValues returned 3"
  // from function ML_back_to_epetraCrs(...) in
  // file ml_epetra_utils.cpp unless the
  // relevant lines are commented out. However this function
  // is much faster (at least on Biowulf) than the alternative
  // below.
  solution_pt = Epetra_MatrixMult(epetra_matrix1_pt, epetra_matrix2_pt);
 }
 else
 {

  // this method requires us to pass in the solution matrix
  solution_pt = new Epetra_CrsMatrix(Copy,*matrix1_range_map_pt,1);
  int epetra_error_flag = 
   EpetraExt::MatrixMatrix::Multiply(*epetra_matrix1_pt,
                                     false,
                                     *epetra_matrix2_pt,
                                     false,
                                     *solution_pt);
#ifdef PARANOID
  if (epetra_error_flag != 0)
   {
    std::ostringstream error_message;
    error_message << "error flag from Multiply(): "
                  << epetra_error_flag
                  << " from TrilinosHelpers::multiply"
                  << std::endl;
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpers::multiply(...)",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif


  // do something with the error flag to keep the compiler from
  // complaining in non-paranoid mode
  epetra_error_flag=0;
  
 }

 // extract values and put into solution
 // ------------------------------------

 // find
 int nnz_local = solution_pt->NumMyNonzeros();
 int nrow_local = matrix1.nrow_local();

 // do some checks
#ifdef PARANOID
 // check number of global rows in soluton matches that in matrix1
 if ( (int)matrix1.nrow() != solution_pt->NumGlobalRows()  )
 {
  std::ostringstream error_message;
  error_message
     << "Incorrect number of global rows in solution matrix. "
     << "nrow() for first input matrix: " << matrix1.nrow()
     << " nrow() for solution: " << solution_pt->NumGlobalRows();
  throw OomphLibError(error_message.str(),
                      "TrilinosHelpers::multiply()",
                      OOMPH_EXCEPTION_LOCATION);
 }

 // check number of local rows in soluton matches that in matrix1
 if ( static_cast<int>(matrix1.nrow_local()) != solution_pt->NumMyRows()  )
 {
  std::ostringstream error_message;
  error_message
     << "Incorrect number of local rows in solution matrix. "
     << "nrow_local() for first input matrix: " << matrix1.nrow_local()
     << " nrow_local() for solution: " << solution_pt->NumMyRows();
  throw OomphLibError(error_message.str(),
                      "TrilinosHelpers::multiply()",
                      OOMPH_EXCEPTION_LOCATION);
 }

 // check number of global columns in soluton matches that in matrix2
 if ( (int)matrix2.ncol() != solution_pt->NumGlobalCols()  )
 {
  std::ostringstream error_message;
  error_message
     << "Incorrect number of global columns in solution matrix. "
     << "ncol() for second input matrix: " << matrix2.ncol()
     << " ncol() for solution: " << solution_pt->NumGlobalCols();
  throw OomphLibError(error_message.str(),
                      "TrilinosHelpers::multiply()",
                      OOMPH_EXCEPTION_LOCATION);
 }

 // check global index of the first row matches
 if ( static_cast<int>(matrix1.first_row()) != solution_pt->GRID(0)  )
 {
  std::ostringstream error_message;
  error_message
     << "Incorrect global index for first row of solution matrix. "
     << "first_row() for first input matrix : " << matrix1.first_row()
     << " first_row() for solution: " << solution_pt->GRID(0);
  throw OomphLibError(error_message.str(),
                      "TrilinosHelpers::multiply()",
                      OOMPH_EXCEPTION_LOCATION);
 }                     
#endif

 // extract values from Epetra matrix row by row
 Vector<double> value(nnz_local);
 Vector<int> column_index(nnz_local);
 Vector<int> row_start(nrow_local+1);
 
 int ptr = 0;
 int num_entries=0;
 double* values=0;
 int* local_column=0;

 for (int row=0; row<nrow_local; row++)
 {
  row_start[row] = ptr;
  solution_pt->ExtractMyRowView(row, num_entries, values, local_column);
  for (int i=0; i<num_entries; i++)
  {
   value[ptr] = *(values+i);
   // obtain global column id and store in column_index
   column_index[ptr] = solution_pt->GCID( *(local_column+i) );
   ptr++;
  }
 }
 row_start[nrow_local] = ptr;

 // delete Trilinos objects
 delete matrix1_range_map_pt;
 delete matrix2_range_map_pt;
 delete matrix2_domain_map_pt;
 delete matrix1_col_map_pt;
 delete matrix2_col_map_pt;
 delete epetra_matrix1_pt;
 delete epetra_matrix2_pt;
 delete solution_pt;
 delete comm_pt;
// delete column_indices1;
// delete column_indices2;
#ifdef OOMPH_HAS_MPI
 delete[] matrix1_global_rows;
 delete[] matrix2_global_rows; 
 delete[] matrix2_global_cols;
#endif

 // Build the Oomph-lib solution matrix using build function
 matrix_soln.rebuild(matrix1.distribution_pt(),
                     matrix2.ncol(),
                     value,
                     column_index,
                     row_start);
}

#ifdef OOMPH_HAS_MPI
//=============================================================================
/// \short Function to compare a distribution described by an oomph-lib 
/// DistributionInfo object to a distribution described by a Trilinos
/// Epetra_Map object. Returns true if they describe the same distribution.
/// Used as a PARANOID check in other TrilinosHelpers methods
//=============================================================================
bool compare_maps(const LinearAlgebraDistribution* oomph_distribution, 
                   const Epetra_Map& epetra_map)
{
 // create the epetra_map for oomph distribution
 Epetra_Map* temp_map_pt = 0;
 int* global_rows;
 Epetra_MpiComm* comm_pt = new Epetra_MpiComm(MPI_COMM_WORLD);
 create_epetra_map(oomph_distribution,comm_pt,temp_map_pt,global_rows);

 // compare
 bool same = epetra_map.SameAs(*temp_map_pt);

 // clean up
 delete temp_map_pt;
 delete[] global_rows;
 delete comm_pt;
 
 // return
 return same;
}
#endif

#ifdef OOMPH_HAS_MPI

 //============================================================================
 /// \short takes a oomph-lib distribution and returns a trilinos mpi map 
 /// (Epetra_Map)
 //============================================================================
 void create_epetra_map(const LinearAlgebraDistribution* distribution_pt,
                        const Epetra_MpiComm* epetra_comm_pt,
                        Epetra_Map* &epetra_map_pt, int* &my_global_rows)
{ 
  // get my first_row and local nrow                                
  unsigned long first_row;
  unsigned long nrow_local;
 
  // this method always creates a distributed epetra map if using
  // multiple procs
  if (distribution_pt->communicator_pt()->nproc() > 1 && 
      distribution_pt->distributed() == false)
    {
      // number of processors
      unsigned nproc = distribution_pt->communicator_pt()->nproc();
      
      // number of global rows
      unsigned global_nrow = distribution_pt->nrow();
      
      // temporay storage for the first rows and nrow locals
      Vector<int> first_row_v(nproc);
      Vector<int> nrow_local_v(nproc);
      
      // compute first row
      for (unsigned p=0;p<nproc;p++)
	{                                          
	  first_row_v[p] = unsigned(double(p*global_nrow)/ 
				    double(nproc));                  
	}             
      
      // compute local nrow
      for (unsigned p=0; p<nproc-1; p++) 
	{                                                  
	  nrow_local_v[p] = first_row_v[p+1] - first_row_v[p];
	}  
      nrow_local_v[nproc-1] = global_nrow - first_row_v[nproc-1];
      
      // store my first row
      unsigned my_rank = distribution_pt->communicator_pt()->my_rank();
      first_row = first_row_v[my_rank];
      nrow_local = nrow_local_v[my_rank];
    }
  else
    {
      first_row = distribution_pt->first_row();                     
      nrow_local = distribution_pt->nrow_local();
    }
  
 // create a distributed trilinos map 
 my_global_rows = new int[nrow_local];
 for (unsigned i=0; i < nrow_local; i++)
  {
   my_global_rows[i] = first_row + i;
  }  
 epetra_map_pt = new Epetra_Map(distribution_pt->nrow(),nrow_local,
                                my_global_rows,0,*epetra_comm_pt); 
}

#else
//=============================================================================
/// \short creates a serial epetra_map of size nrow
//=============================================================================
void create_epetra_map(const LinearAlgebraDistribution* dist_pt,
                       const Epetra_SerialComm* epetra_comm_pt,
                       Epetra_Map* &epetra_map_pt)
{
 // create a simple Trilinos map (contiguous)
 epetra_map_pt = new Epetra_Map(dist_pt->nrow(),0,*epetra_comm_pt);
}
#endif

}
}
