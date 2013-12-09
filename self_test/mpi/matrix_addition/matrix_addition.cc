//LIC//====================================================================
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

// Oomph-lib includes
#include "generic.h"

using namespace oomph;

// Helper function to construct a Vector given an array
template<typename myType>
void construct_vector(myType given_array[],unsigned given_arraysize,
                      Vector<myType> &result_vector)
{
  // Clear and reserve the required memory.
  result_vector.clear();
  result_vector.reserve(given_arraysize);

  for (unsigned i = 0; i < given_arraysize; i++) 
  {
    result_vector.push_back(given_array[i]);
  }
}

// Helper function to output a Vector.
template<typename myType>
void output_vector(Vector<myType> &given_vector)
{
  typename Vector<myType>::iterator it;

  for(it = given_vector.begin(); it != given_vector.end(); ++it)
  {
    oomph_info << *it << std::endl; 
  }
}

//===start_of_main=============================================================
/// Driver code: Testing CRDoubleMatrixHelpers::add(...)
///
/// The script validate.sh should run this self test on 1, 2, 3 and 4 cores.
//=============================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Get the global oomph-lib communicator 
  const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

  // my rank and number of processors. This is used later for putting the data.
  unsigned my_rank = comm_pt->my_rank();
  unsigned nproc = comm_pt->nproc();

  // Matrix 0
  // [1 2 0 0 0
  //  3 0 4 0 0
  //  0 0 0 5 6
  //  7 8 0 0 9
  //  0 0 0 0 0]
  //
  //  values = 1 2 3 4 5 6 7 8 9
  //  col i  = 0 1 0 2 3 4 0 1 4
  //  row s  = 0 2 4 6 9 9
  //
  //  Matrix 1
  //  [1  2  3  4  5
  //   6  7  8  9  10
  //   0  0  11 0  12
  //   13 14 15 0  0
  //   0  0  16 17 18]
  //  
  //  values = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
  //  col i  = 0 1 2 3 4 0 1 2 3 4  2  4  0  1  2  2  3  4
  //  row s  = 0 5 10 12 15 18
  //
  //  Result matrix
  //  [2  4  3  4  5 
  //   9  7  12 9  10
  //   0  0  11 5  18
  //   20 22 15 0  9 
  //   0  0  16 17 18]
  //
  //   values = 2 4 3 4 5 9 7 12 9 10 11 5 18 20 22 15 9 16 17 18
  //   col i  = 0 1 2 3 4 0 1 2 3 4 2 3 4 0 1 2 4 2 3 4
  //   row s  = 0 5 10 13 17 19
  //

  // First matrix (mat_zero)
  unsigned nrow_global_t1m1 = 5;
  unsigned ncol_t1m1 = 5;
  unsigned nnz_t1m1 = 9;
  Vector<double> val_t1m1;
  Vector<int> col_i_t1m1;
  Vector<int> row_s_t1m1;

  double val_array_t1m1[] = {1,2,3,4,5,6,7,8,9};
  int col_i_array_t1m1[] = {0,1,0,2,3,4,0,1,4};
  int row_s_array_t1m1[] = {0,2,4,6,9,9};

  construct_vector(val_array_t1m1,nnz_t1m1,val_t1m1);
  construct_vector(col_i_array_t1m1,nnz_t1m1,col_i_t1m1);
  construct_vector(row_s_array_t1m1,nrow_global_t1m1+1,row_s_t1m1);

  LinearAlgebraDistribution distri_zero(comm_pt,nrow_global_t1m1,true);

  CRDoubleMatrix mat_t1m1;
  CRDoubleMatrixHelpers::create_uniformly_distributed_matrix(
      nrow_global_t1m1,ncol_t1m1,comm_pt,
      val_t1m1,col_i_t1m1,row_s_t1m1,mat_t1m1);

  std::ostringstream mat_t1m1_stream;
  mat_t1m1_stream << "t1m1_NP"<<nproc<<"R"<<my_rank;
  mat_t1m1.sparse_indexed_output(mat_t1m1_stream.str());

  // Next matrix (mat_t1m2)
  unsigned nrow_global_t1m2 = 5;
  unsigned ncol_t1m2 = 5;
  unsigned nnz_t1m2 = 18;
  Vector<double> val_t1m2;
  Vector<int> col_i_t1m2;
  Vector<int> row_s_t1m2;

  double val_array_t1m2[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
  int col_i_array_t1m2[] = {0,1,2,3,4,0,1,2,3,4,2,4,0,1,2,2,3,4};
  int row_s_array_t1m2[] = {0,5,10,12,15,18};

  construct_vector(val_array_t1m2,nnz_t1m2,val_t1m2);
  construct_vector(col_i_array_t1m2,nnz_t1m2,col_i_t1m2);
  construct_vector(row_s_array_t1m2,nrow_global_t1m2+1,row_s_t1m2);

  LinearAlgebraDistribution distri_t1m2(comm_pt,nrow_global_t1m2,true);

  CRDoubleMatrix mat_t1m2;
  CRDoubleMatrixHelpers::create_uniformly_distributed_matrix(
      nrow_global_t1m2,ncol_t1m2,comm_pt,
      val_t1m2,col_i_t1m2,row_s_t1m2,mat_t1m2);

  std::ostringstream mat_t1m2_stream;

  mat_t1m2_stream << "t1m2_NP"<<nproc<<"R"<<my_rank;
  mat_t1m2.sparse_indexed_output(mat_t1m2_stream.str());

  CRDoubleMatrix mat_t1_result;

  mat_t1m1.add(mat_t1m2,mat_t1_result);

  std::ostringstream mat_t1_result_stream;
  mat_t1_result_stream << "t1_res_NP"<<nproc<<"P"<<my_rank;
  mat_t1_result.sparse_indexed_output(mat_t1_result_stream.str());


  // Now test with the following matrices: A + B = C where
  // A = [1  2  3  4  5  6
  //      7  8  9  10 11 12
  //      13 14 15 16 17 18]
  //
  // B = [1  2  3  4  5  6
  //      7  8  9  10 11 12
  //      13 14 15 16 17 18]
  //
  // B = [2  4  6  8  10 12
  //      14 16 18 20 22 24
  //      26 28 30 32 34 36]
  //
  



  //  // The number of block rows and columns, these will be set to the correct
  //  // value per test below.
  //  unsigned nblock_row = 0;
  //  unsigned nblock_col = 0;
  //
  //  // Supply the dimensions of the matrices to concatenate.
  //  // The matrices must be in the order: top row, then along the columns,
  //  // then second row, et cetera...
  //  // The number of elements in this must be 2*nblock_row*nblock_col.
  //  unsigned dimarray[] = {7,7,7,5,7,3,5,7,5,5,5,3,3,7,3,5,3,3};
  //
  //  /////////////////////////////////////////////////////////////////////////////
  //  // Test 0: Concatenate the sub matrices with the sizes
  //  // (7,7)(7,5)(7,3)
  //  // (5,7)(5,5)(5,3)
  //  // (3,7)(3,5)(3,3)
  //  
  //  // This is a 3 by 3 block matrix.
  //  nblock_row = 3;
  //  nblock_col = 3;
  //
  //  // The data structure to store the pointers to matrices.
  //  DenseMatrix<CRDoubleMatrix*> mat0_pt(nblock_row,nblock_col,0);
  //
  //  // Create the matrice to concatenate.
  //  create_matrices_to_cat(nblock_row,nblock_col,dimarray,
  //                         comm_pt,mat0_pt);
  //  
  //  // The result matrix.
  //  CRDoubleMatrix result_matrix0;
  //
  //  // Call the concatenate function.
  //  CRDoubleMatrixHelpers::concatenate(mat0_pt,result_matrix0);
  //
  //  // output the result matrix
  //  std::ostringstream result_mat0_stream;
  //  result_mat0_stream << "out0_NP" << nproc << "R" << my_rank;
  //  result_matrix0.sparse_indexed_output(result_mat0_stream.str().c_str());
  //
  //  // No longer need the result matrix.
  //  result_matrix0.clear();
  //
  //  /////////////////////////////////////////////////////////////////////////////
  //  // Test 1: Concatenate the sub matrices with the sizes
  //  // (7,7)(7,5)
  //  // (5,7)(5,5)
  //  // (3,7)(3,5)
  //  
  //  // This is a 3 by 3 block matrix.
  //  nblock_row = 3;
  //  nblock_col = 2;
  //
  //  // The data structure to store the pointers to matrices.
  //  DenseMatrix<CRDoubleMatrix*> mat1_pt(nblock_row,nblock_col,0);
  //  
  //  // Get the matrices from mat0_pt
  //  for(unsigned block_row_i = 0; block_row_i < nblock_row; block_row_i++) 
  //  {
  //    for (unsigned block_col_i = 0; block_col_i < nblock_col; block_col_i++) 
  //    {
  //      mat1_pt(block_row_i,block_col_i) = mat0_pt(block_row_i,block_col_i);
  //    }
  //  }
  //
  //  // The result matrix.
  //  CRDoubleMatrix result_matrix1;
  //
  //  // Call the concatenate function.
  //  CRDoubleMatrixHelpers::concatenate(mat1_pt,result_matrix1);
  //
  //  // output the result matrix
  //  std::ostringstream result_mat1_stream;
  //  result_mat1_stream << "out1_NP" << nproc << "R" << my_rank;
  //  result_matrix1.sparse_indexed_output(result_mat1_stream.str().c_str());
  //
  //  // No longer need the result matrix.
  //  result_matrix1.clear();
  //
  //  /////////////////////////////////////////////////////////////////////////////
  //  // Test 2: Concatenate the sub matrices with the sizes
  //  // (7,7)(7,5)(7,3)
  //  // (5,7)(5,5)(5,3)
  //  
  //  // This is a 3 by 3 block matrix.
  //  nblock_row = 2;
  //  nblock_col = 3;
  //
  //  // The data structure to store the pointers to matrices.
  //  DenseMatrix<CRDoubleMatrix*> mat2_pt(nblock_row,nblock_col,0);
  //  
  //  // Get the matrices from mat0_pt
  //  for(unsigned block_row_i = 0; block_row_i < nblock_row; block_row_i++) 
  //  {
  //    for (unsigned block_col_i = 0; block_col_i < nblock_col; block_col_i++) 
  //    {
  //      mat2_pt(block_row_i,block_col_i) = mat0_pt(block_row_i,block_col_i);
  //    }
  //  }
  //
  //  // The result matrix.
  //  CRDoubleMatrix result_matrix2;
  //
  //  // Call the concatenate function.
  //  CRDoubleMatrixHelpers::concatenate(mat2_pt,result_matrix2);
  //
  //  // output the result matrix
  //  std::ostringstream result_mat2_stream;
  //  result_mat2_stream << "out2_NP" << nproc << "R" << my_rank;
  //  result_matrix2.sparse_indexed_output(result_mat2_stream.str().c_str());
  //
  //  // No longer need the result matrix.
  //  result_matrix2.clear();
  //
  //  /////////////////////////////////////////////////////////////////////////////
  //  // Delete the sub block matrices, we only need to delete them from mat0_pt.
  //  nblock_row = mat0_pt.nrow();
  //  nblock_col = mat0_pt.ncol();
  //  for (unsigned block_row_i = 0; block_row_i < nblock_row; block_row_i++) 
  //  {
  //    for (unsigned block_col_i = 0; block_col_i < nblock_col; block_col_i++) 
  //    {
  //      delete mat0_pt(block_row_i,block_col_i);
  //    } // for col_i
  //  } // for row_i

#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main
