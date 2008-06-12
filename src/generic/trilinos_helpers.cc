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

#ifdef OOMPH_HAS_MPI



 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
 // Distributed Helpers
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
 


 //============================================================================
 /// \short creates a distributed Epetra_Vector that has the same contents as
 /// oomph_v with distribution row_map_pt. \n
 /// \b NOTE 1. the Epetra_Map and the DistributionInfo objects must descibe
 /// the same distribution
 /// \b NOTE 2. if the bool copy is true (default) then the values in the 
 /// oomph-lib vector (oomph_v) will be copied into the epetra vector, or if 
 /// it is false then the epetra vector will 'view' (i.e. point to) the 
 /// contents of the oomph-lib vector. (defult is copy = true)
 //============================================================================
 void create_epetra_vector(DistributedVector<double>& oomph_v,
                           const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt,
                           const bool& copy)
 {
  // paranoid check that the distribution of the map and the vector are the
  // same
#ifdef PARANOID
  if (!compare_maps(oomph_v.distribution(),*row_map_pt))
   {
    std::ostringstream error_message;
    error_message << "The oomph-lib distributed vector (oomph_v) and the "
                  << "epetra map must describe the same distribution.";
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpers::create_epetra_vector()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // create the Epetra_Vector
  if (copy)
   {
    epetra_v_pt = new Epetra_Vector(Copy,*row_map_pt,&oomph_v[0]);
   }
  else
   {
    epetra_v_pt = new Epetra_Vector(View,*row_map_pt,&oomph_v[0]);
   }
 }


 //============================================================================
 /// creates a distributed Epetra_Vector that has the same contents as
 /// oomph_v with distribution row_map_pt, the contents of the trilinos vector 
 /// are copied from the oomph-lib vector. \n
 /// \b NOTE 1. duplication of helper due to the const oomph-lib vector
 /// \b NOTE 2. the Epetra_Map and the DistributionInfo objects must descibe
 /// the same distribution
 //============================================================================
 void create_epetra_vector(const DistributedVector<double>& oomph_v,
                           const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt)
 {
  // paranoid check that the distribution of the map and the vector are the
  // same
#ifdef PARANOID
  if (!compare_maps(oomph_v.distribution(),*row_map_pt))
   {
    std::ostringstream error_message;
    error_message << "The oomph-lib distributed vector (oomph_v) and the "
                  << "epetra map must describe the same distribution.";
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpers::create_epetra_vector()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // create the Epetra_Vector
  unsigned nrow_local = oomph_v.distribution().nrow_local();
  double* v = new double[nrow_local];
  for (unsigned i = 0; i < nrow_local; i++)
   {
    v[i] = oomph_v[i];
   }
  epetra_v_pt = new Epetra_Vector(Copy,*row_map_pt,v);
  delete[] v;
 }

 //============================================================================
 /// Helper function to copy the contents of a Trilinos vector to an
 /// oomph-lib distributed vector. The my_distribution DistributionInfo should 
 /// be the same as the distribution of epetra_v_pt
 //============================================================================
 void copy_to_oomphlib_vector(const Epetra_Vector* epetra_v_pt,
                              const DistributionInfo my_distribution,
                              DistributedVector<double>& oomph_v)
 {
  // extract values from epetra_v_pt
  double* v_values;
  epetra_v_pt->ExtractView(&v_values);

  // distribute the oomph-lib vector
  oomph_v.distribute(my_distribution);

  // copy the values
  unsigned nrow_local = my_distribution.nrow_local();
  for (unsigned i = 0; i < nrow_local; i++)
   {
    oomph_v[i] = v_values[i];
   }
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
void multiply(DistributedCRDoubleMatrix &oomph_matrix,
              const DistributedVector< double > &oomph_x,
              DistributedVector< double > &oomph_soln)
{
#ifdef PARANOID
 // check matrix and vector are compatible
 if ( oomph_matrix.ncol() != oomph_x.nrow_global()  )
 {
  std::ostringstream error_message;
  error_message
   << "Matrix dimensions incompatible for matrix-matrix multiplication"
   << "ncol() for the matrix: " << oomph_matrix.ncol()
   << "nrow() for the vector: " << oomph_x.nrow_global();
  throw OomphLibError(error_message.str(),
                      "TrilinosHelpers::multiply(...)",
                      OOMPH_EXCEPTION_LOCATION);
 }
 
 // check that matrix and x vector have the same communicator
 int comm_compare_result;
 MPI_Comm_compare(oomph_matrix.distribution().communicator(), 
                 oomph_x.distribution().communicator(),&comm_compare_result);
 if (comm_compare_result != MPI_IDENT)
  {
   std::ostringstream error_message;
   error_message
    << "Matrix (oomph_matrix) and rhs vector (oomph_x) must have the same " 
    << "communicator.\n";
   throw OomphLibError(error_message.str(),
                      "TrilinosHelpers::multiply(...)",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 // create the communicator
 // note use of MPI_COMM_WORLD
 Epetra_MpiComm* comm_pt = new Epetra_MpiComm(MPI_COMM_WORLD);

 // create the maps corresponding the distribution of matrix1 and matrix2
 Epetra_Map* matrix_row_map_pt;
 int* global_rows;
 create_epetra_map(oomph_matrix.distribution(),comm_pt,matrix_row_map_pt,
                   global_rows);
 Epetra_Map* matrix_col_map_pt;
 int* global_cols;
 create_epetra_map(oomph_x.distribution(),comm_pt,matrix_col_map_pt,
                   global_cols); 

 // convert matrix1 to epetra matrix
 Epetra_CrsMatrix* epetra_matrix_pt;
 create_epetra_matrix(&oomph_matrix,matrix_row_map_pt,matrix_col_map_pt,
                      epetra_matrix_pt);

 // convert x to Trilinos vector
 Epetra_Vector* epetra_x_pt;
 create_epetra_vector(oomph_x,matrix_col_map_pt,epetra_x_pt);

 // distribute the vector
 oomph_soln.distribute(oomph_matrix.distribution());

 // create Trilinos vector for soln ( 'viewing' the contents of the oomph-lib
 // matrix )
 Epetra_Vector* epetra_soln_pt;
 create_epetra_vector(oomph_soln,matrix_row_map_pt,epetra_soln_pt,false);

 // do the multiply
 int epetra_error_flag = epetra_matrix_pt->Multiply(false,*epetra_x_pt,
                                                    *epetra_soln_pt);

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
 delete matrix_row_map_pt;
 delete[] global_rows;
 delete matrix_col_map_pt;
 delete[] global_cols;
 delete epetra_matrix_pt;
 delete epetra_x_pt;
 delete epetra_soln_pt;
 delete comm_pt;
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
void multiply(DistributedCRDoubleMatrix &matrix1,
              DistributedCRDoubleMatrix &matrix2,
              DistributedCRDoubleMatrix &matrix_soln,
              const bool& use_ml)
{

#ifdef PARANOID
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
 
 // check that matrices have the same communicator
  int comm_compare_result;
 MPI_Comm_compare(matrix1.distribution().communicator(), 
                 matrix2.distribution().communicator(),&comm_compare_result);
 if (comm_compare_result != MPI_IDENT)
  {
   std::ostringstream error_message;
   error_message
    << "Input matrices must have the same communicator.\n";
   throw OomphLibError(error_message.str(),
                      "TrilinosHelpers::multiply(...)",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 // create the communicator pt
 Epetra_MpiComm* comm_pt = new Epetra_MpiComm(MPI_COMM_WORLD);

 // create the maps corresponding the distribution of matrix1 and matrix2
 Epetra_Map* matrix1_map_pt;
 int* matrix1_global_rows;
 create_epetra_map(matrix1.distribution(),comm_pt,
                   matrix1_map_pt,matrix1_global_rows);
 Epetra_Map* matrix2_map_pt;
 int* matrix2_global_rows;
 create_epetra_map(matrix2.distribution(),comm_pt,
                   matrix2_map_pt,matrix2_global_rows); 
 Epetra_Map* matrix2_col_map_pt;
 int* matrix2_global_cols;
 DistributionInfo matrix2_col_distribution
  (matrix2.distribution().communicator(),matrix2.ncol()); 
 create_epetra_map(matrix2_col_distribution,comm_pt,matrix2_col_map_pt,
                   matrix2_global_cols);  

 // convert matrix1 to epetra matrix
 Epetra_CrsMatrix* epetra_matrix1_pt;
 create_epetra_matrix(&matrix1,matrix1_map_pt,matrix2_map_pt,epetra_matrix1_pt);

 // convert matrix2 to epetra matrix
 Epetra_CrsMatrix* epetra_matrix2_pt;
 create_epetra_matrix(&matrix2,matrix2_map_pt,matrix2_col_map_pt,
                      epetra_matrix2_pt);

 // create Trilinos matrix to hold solution - will have same map 
 // (and number of rows) as matrix1
 Epetra_CrsMatrix* solution_pt;

 // do the multiplication
 // ---------------------
 if (use_ml)
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
  solution_pt = new Epetra_CrsMatrix(Copy,*matrix1_map_pt,
                                     *matrix2_col_map_pt,1);
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
 delete matrix1_map_pt;
 delete[] matrix1_global_rows;
 delete matrix2_map_pt;
 delete[] matrix2_global_rows; 
 delete matrix2_col_map_pt;
 delete[] matrix2_global_cols;
 delete epetra_matrix1_pt;
 delete epetra_matrix2_pt;
 delete solution_pt;
 delete comm_pt;

 // Build the Oomph-lib solution matrix using build function
 matrix_soln.build(value,
                   column_index,
                   row_start,
                   matrix1.distribution(),
                   matrix2.ncol());
}


//=============================================================================
/// \short takes a oomph-lib distribution (DistributionInfo) and returns
/// a trilinos map (Epetra_Map)
//=============================================================================
void create_epetra_map(const DistributionInfo& my_distribution,
                       const Epetra_MpiComm* epetra_comm_pt,
                       Epetra_Map* &epetra_map_pt, int* &my_global_rows)

{ 
 // get my first_row and local nrow                                
 unsigned long first_row = my_distribution.first_row();                     
 unsigned long nrow_local = my_distribution.nrow_local(); 
 
 // create a distributed trilinos map 
 my_global_rows = new int[nrow_local];
 for (unsigned i=0; i < nrow_local; i++)
  {
   my_global_rows[i] = first_row + i;
  }  
 epetra_map_pt = new Epetra_Map(my_distribution.nrow_global(),nrow_local,
                                my_global_rows,0,*epetra_comm_pt); 
}

//=============================================================================
/// \short Function to compare a distribution described by an oomph-lib 
/// DistributionInfo object to a distribution described by a Trilinos
/// Epetra_Map object. Returns true if they describe the same distribution.
/// Used as a PARANOID check in other TrilinosHelpers methods
//=============================================================================
bool compare_maps(const DistributionInfo& oomph_distribution, 
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
#else




 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
 // SERIAL Helpers
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////



/// \short creates a serial epetra_map of size nrow
void create_epetra_map(const unsigned& nrow,
                       const Epetra_SerialComm* epetra_comm_pt,
                       Epetra_Map* &epetra_map_pt)
{
 // create a simple Trilinos map (contiguous)
 epetra_map_pt = new Epetra_Map(nrow,0,*epetra_comm_pt);
}

#endif



 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
 // SERIAL and DISTRIBUTED Helpers
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////



 //============================================================================
 /// \short creates an Epetra_Vector from a serial oomph-lib vector oomph_v
 /// with map row_map_pt.
 //============================================================================
 void create_epetra_vector(const Vector<double>& oomph_v,
                           const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt)
 {
  // paranoid check that the oomphlib vector and epetra map have the same 
  // number of global rows
#ifdef PARANOID
  if (static_cast<int>(oomph_v.size()) != row_map_pt->NumGlobalElements())
   {
    std::ostringstream error_message;
    error_message
     << "The epetra map and the oomph-lib matrix have different numbers "
     << "of global rows.\n"
     << "NumGlobalElements() for the Epetra_Map: " 
     << row_map_pt->NumGlobalElements() << "\n" 
     << "size() for Vector<double>: " << oomph_v.size();
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpers::create_epetra_vector()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // storage for data to be inserted
  double* values;
  int* indices;
  unsigned nrow = 0;

#ifdef OOMPH_HAS_MPI
  int nproc = row_map_pt->Comm().NumProc();
  if (nproc == 1)
   {
    // store the data for insertion
    nrow = oomph_v.size();
    values = new double[nrow];
    indices = new int[nrow];
    for (unsigned i = 0; i < nrow; i++)
     {
      values[i] = oomph_v[i];
      indices[i] = i;
     }
   }
  else
   {
    // get my nrow_local and first_row
    unsigned nrow_local = row_map_pt->NumMyElements();
    int* global_rows = row_map_pt->MyGlobalElements();
    unsigned first_row = static_cast<unsigned>(global_rows[0]);

    // store the data for insertion
    nrow = nrow_local;
    values = new double[nrow];
    indices = new int[nrow];
    for (unsigned i = 0; i < nrow; i++)
     {
      values[i] = oomph_v[i + first_row];
      indices[i] = i + first_row;
     }
   }
#else
  // store the data for insertion
  nrow = oomph_v.size();
  values = new double[nrow];
  indices = new int[nrow];
  for (unsigned i = 0; i < nrow; i++)
   {
    values[i] = oomph_v[i];
    indices[i] = i;
   }
#endif

  // insert
  epetra_v_pt = new Epetra_Vector(*row_map_pt);
  epetra_v_pt->ReplaceGlobalValues(nrow,values,indices);

  // clean up memory
  delete[] values;
  delete[] indices;
 }

 //============================================================================
 /// creates an empty trilinos vector with Epetra_Map row_map_pt
 //============================================================================
 void create_epetra_vector(const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt)
 {

  // create the Epetra_Vector
  epetra_v_pt = new Epetra_Vector(*row_map_pt);
 }


//=============================================================================
/// \short Helper function to copy the contents of a epetra vector to an
/// oomph-lib distributed vector.
//=============================================================================
void copy_to_oomphlib_vector(const Epetra_Vector* epetra_v_pt,
                             Vector<double>& oomph_lib_v)
{
 // number of global rows
 unsigned nrow = epetra_v_pt->GlobalLength();

 // values to be copied into the oomph-lib vector
 double* values;

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
  
   // get the communicator
   MPI_Comm my_mpi_comm = 
    dynamic_cast<const Epetra_MpiComm&>(epetra_v_pt->Map().Comm()).Comm();

   // gather the First_row vector
   MPI_Allgather(&nrow_local[my_rank],1,
                 MPI_INT,
                 &nrow_local[0],1,MPI_INT,my_mpi_comm);

   // number of local rows
   Vector<int> first_row(nproc);
   int* global_rows = epetra_v_pt->Map().MyGlobalElements();
   first_row[my_rank] = global_rows[0];

   // gather the First_row vector
   MPI_Allgather(&first_row[my_rank],1,
                 MPI_INT,
                 &first_row[0],1,MPI_INT,
                 my_mpi_comm);

   // gather the local solution values
   values = new double[nrow];
   MPI_Allgatherv(local_values,
                  nrow_local[my_rank],
                  MPI_DOUBLE,
                  values,
                  &nrow_local[0],
                  &first_row[0],
                  MPI_DOUBLE,
                  my_mpi_comm);
   
   // clean up memory
   MPI_Comm_free(&my_mpi_comm);
  }
#else
  epetra_v_pt->ExtractView(&values);
#endif

  // copy the vector to the oomphlib vector
  oomph_lib_v.resize(nrow);
  for (unsigned i=0; i<nrow; i++)
   {
    oomph_lib_v[i] = values[i];
   }

  // clean up memory
#ifdef OOMPH_HAS_MPI
  if (nproc > 1)
   {
    delete[] values;
   }
#endif
}

//=============================================================================
/// \short Helper function to create a distributed Epetra_CrsMatrix from a 
/// from a DistributedCRDoubleMatrix or a CRDoubleMatrix
/// \b NOTE 1. This function constructs an Epetra_CrsMatrix using new,
/// "delete trilinos_matrix_pt;" is NOT called.
//=============================================================================
void create_epetra_matrix(DoubleMatrixBase* matrix_pt,
                          const Epetra_Map* epetra_row_map_pt,
                          const Epetra_Map* epetra_col_map_pt,
                          Epetra_CrsMatrix* &epetra_matrix_pt)
{
 // if this is a MPI build then attempt to build Trilinos Matrix from
 // a CRDoubleMatrix or a DistributedCRDoubleMatrix
#ifdef OOMPH_HAS_MPI
 bool cast_failed = true;

 // first try DistributedCRDoubleMatrix
 DistributedCRDoubleMatrix* cast_dist_matrix_pt = 
  dynamic_cast<DistributedCRDoubleMatrix*>(matrix_pt);
 if (cast_dist_matrix_pt!=0)
  {

   // paranoid check that the distribution of the epetra row map and the 
   // oomph-lib matrix are the same
#ifdef PARANOID
   compare_maps(cast_dist_matrix_pt->distribution(),
                *epetra_row_map_pt);
#endif
     // paranoid check that the distribution of the map and the vector are the
  // same
#ifdef PARANOID
  if (!compare_maps(cast_dist_matrix_pt->distribution(),
                    *epetra_row_map_pt))
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
   int* column = cast_dist_matrix_pt->column_index();
   double* value = cast_dist_matrix_pt->value();
   int* row_start = cast_dist_matrix_pt->row_start();
   
   // get my nrow_local and first_row
   unsigned nrow_local = cast_dist_matrix_pt->nrow_local();
   unsigned first_row = cast_dist_matrix_pt->first_row();
   
   // store the number of non zero entries per row
   int* nnz_per_row = new int[nrow_local];
   for (unsigned row=0;row<nrow_local;row++)
    {
     nnz_per_row[row] = row_start[row+1] - row_start[row];
    }
   
   // construct the new Epetra matrix
   epetra_matrix_pt = new Epetra_CrsMatrix(Copy,*epetra_row_map_pt,
                                           nnz_per_row);
   
   // insert the values
   for (unsigned row=0; row<nrow_local; row++)
    {
     // get pointer to this row in values/columns
     int ptr = row_start[row];
     epetra_matrix_pt->InsertGlobalValues(first_row+row,
                                          nnz_per_row[row],
                                          value+ptr,
                                          column+ptr);
    }
   
   // tidy up memory
   delete[] nnz_per_row;
   
   // complete the build of the trilinos matrix
#ifdef PARANOID
   int err=0;
   err = epetra_matrix_pt->FillComplete(*epetra_col_map_pt,*epetra_row_map_pt);
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
   epetra_matrix_pt->FillComplete(*epetra_col_map_pt,*epetra_row_map_pt);
#endif
   
   // cast succeeded
   cast_failed = false;
  }

 // if the cast failed try cast to CRDoubleMatrix
 if (cast_failed)
 {
  CRDoubleMatrix* cast_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);
  if (cast_matrix_pt != 0)
  {

   // get pointers to the matrix values, column indices etc
   int* column = cast_matrix_pt->column_index();
   double* value = cast_matrix_pt->value();
   int* row_start = cast_matrix_pt->row_start();
   
   // get my nrow_local and first_row
   unsigned nrow_local = epetra_row_map_pt->NumMyElements();
   int* global_rows = epetra_row_map_pt->MyGlobalElements();
   unsigned first_row = static_cast<unsigned>(global_rows[0]);

   // store the number of non zero entries per row
   int* nnz_per_row = new int[nrow_local];
   for (unsigned row=0;row<nrow_local;row++)
    {
     nnz_per_row[row] = row_start[first_row + row+1] 
      - row_start[first_row + row];
    }
   
   // construct the new Epetra matrix
   epetra_matrix_pt = new Epetra_CrsMatrix(Copy,*epetra_row_map_pt,
                                           nnz_per_row);

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
   err = epetra_matrix_pt->FillComplete(*epetra_col_map_pt,*epetra_row_map_pt);
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
   epetra_matrix_pt->FillComplete(*epetra_col_map_pt,*epetra_row_map_pt);
#endif
   
   // cast succeeded
   cast_failed = false;
  }
 }
 
#ifdef PARANOID
 // check the matrix is a DistributedCRDoubleMatrix
 if (cast_failed)
  {
   std::ostringstream error_message;
   error_message << "An create_epetra_matrix(...) is only compatiible with "
                 << "DistributedCRDoubleMatrix or CRDoubleMatrix" << std::endl;
   throw OomphLibError(error_message.str(),
                       "TrilinosHelpers::create_epetra_matrix()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
#else

 CRDoubleMatrix* cast_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);
#ifdef PARANOID
 if (cast_matrix_pt == 0)
  {
   std::ostringstream error_message;
   error_message << "An create_epetra_matrix(...) is only compatiible with "
                 << "CRDoubleMatrix" << std::endl;
   throw OomphLibError(error_message.str(),
                       "TrilinosHelpers::create_epetra_matrix()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

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
 epetra_matrix_pt = new Epetra_CrsMatrix(Copy,*epetra_row_map_pt,
                                         nnz_per_row);

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
 err = epetra_matrix_pt->FillComplete(*epetra_col_map_pt,*epetra_row_map_pt);
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
 epetra_matrix_pt->FillComplete(*epetra_col_map_pt,*epetra_row_map_pt);
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
                          const Epetra_Map* epetra_row_map_pt,
                          Epetra_CrsMatrix* &epetra_matrix_pt)
{
 create_epetra_matrix(oomph_matrix,epetra_row_map_pt,epetra_row_map_pt,
                      epetra_matrix_pt);
}
}
}
