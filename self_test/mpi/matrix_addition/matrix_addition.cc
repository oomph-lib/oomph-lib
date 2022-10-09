//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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

  // First matrix (t1m1)
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

  CRDoubleMatrix mat_t1m1;
  CRDoubleMatrixHelpers::create_uniformly_distributed_matrix(
      nrow_global_t1m1,ncol_t1m1,comm_pt,
      val_t1m1,col_i_t1m1,row_s_t1m1,mat_t1m1);

  std::ostringstream mat_t1m1_stream;
  mat_t1m1_stream << "t1m1_NP"<<nproc<<"R"<<my_rank;
  mat_t1m1.sparse_indexed_output(mat_t1m1_stream.str());

  // Next matrix (t1m2)
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

  CRDoubleMatrix mat_t1m2;
  CRDoubleMatrixHelpers::create_uniformly_distributed_matrix(
      nrow_global_t1m2,ncol_t1m2,comm_pt,
      val_t1m2,col_i_t1m2,row_s_t1m2,mat_t1m2);

  std::ostringstream mat_t1m2_stream;

  mat_t1m2_stream << "t1m2_NP"<<nproc<<"R"<<my_rank;
  mat_t1m2.sparse_indexed_output(mat_t1m2_stream.str());

  // Now perform element-wise addition of t1m1 and t1m2, and put the result
  // in mat_t1_result.
  CRDoubleMatrix mat_t1_result;

  mat_t1m1.add(mat_t1m2,mat_t1_result);

  std::ostringstream mat_t1_result_stream;
  mat_t1_result_stream << "t1_res_NP"<<nproc<<"R"<<my_rank;
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
  
  // First matrix (t2m1)
  unsigned nrow_global_t2m1 = 3;
  unsigned ncol_t2m1 = 6;
  unsigned nnz_t2m1 = 18;
  Vector<double> val_t2m1;
  Vector<int> col_i_t2m1;
  Vector<int> row_s_t2m1;

  double val_array_t2m1[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
  int col_i_array_t2m1[] = {0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5};
  int row_s_array_t2m1[] = {0,6,12,18};

  construct_vector(val_array_t2m1,nnz_t2m1,val_t2m1);
  construct_vector(col_i_array_t2m1,nnz_t2m1,col_i_t2m1);
  construct_vector(row_s_array_t2m1,nrow_global_t2m1+1,row_s_t2m1);

  CRDoubleMatrix mat_t2m1;
  CRDoubleMatrixHelpers::create_uniformly_distributed_matrix(
      nrow_global_t2m1,ncol_t2m1,comm_pt,
      val_t2m1,col_i_t2m1,row_s_t2m1,mat_t2m1);

  std::ostringstream mat_t2m1_stream;
  mat_t2m1_stream << "t2m1_NP"<<nproc<<"R"<<my_rank;
  mat_t2m1.sparse_indexed_output(mat_t2m1_stream.str());
 
  // Second matrix (t2m2)
  unsigned nrow_global_t2m2 = 3;
  unsigned ncol_t2m2 = 6;
  unsigned nnz_t2m2 = 18;
  Vector<double> val_t2m2;
  Vector<int> col_i_t2m2;
  Vector<int> row_s_t2m2;

  double val_array_t2m2[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
  int col_i_array_t2m2[] = {0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5};
  int row_s_array_t2m2[] = {0,6,12,18};

  construct_vector(val_array_t2m2,nnz_t2m2,val_t2m2);
  construct_vector(col_i_array_t2m2,nnz_t2m2,col_i_t2m2);
  construct_vector(row_s_array_t2m2,nrow_global_t2m2+1,row_s_t2m2);

  CRDoubleMatrix mat_t2m2;
  CRDoubleMatrixHelpers::create_uniformly_distributed_matrix(
      nrow_global_t2m2,ncol_t2m2,comm_pt,
      val_t2m2,col_i_t2m2,row_s_t2m2,mat_t2m2);

  std::ostringstream mat_t2m2_stream;
  mat_t2m2_stream << "t2m2_NP"<<nproc<<"R"<<my_rank;
  mat_t2m2.sparse_indexed_output(mat_t2m2_stream.str());

  // Now perform element-wise addition of t2m1 and t2m2, and put the result
  // in mat_t2_result.
  CRDoubleMatrix mat_t2_result;

  mat_t2m1.add(mat_t2m2,mat_t2_result);

  std::ostringstream mat_t2_result_stream;
  mat_t2_result_stream << "t2_res_NP"<<nproc<<"R"<<my_rank;
  mat_t2_result.sparse_indexed_output(mat_t2_result_stream.str());


  // Now test with the following matrices: A + B = C where
  // A = [1  2  3
  //      4  5  6
  //      7  8  9
  //      10 11 12
  //      13 14 15
  //      16 17 18]
  //
  // B = [1  2  3
  //      4  5  6
  //      7  8  9
  //      10 11 12
  //      13 14 15
  //      16 17 18]
  //
  // C = [2  4  6
  //      8  10 12
  //      14 16 18
  //      20 22 24
  //      26 28 30
  //      32 34 36]
  
  // First matrix (t3m1)
  unsigned nrow_global_t3m1 = 6;
  unsigned ncol_t3m1 = 3;
  unsigned nnz_t3m1 = 18;
  Vector<double> val_t3m1;
  Vector<int> col_i_t3m1;
  Vector<int> row_s_t3m1;

  double val_array_t3m1[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
  int col_i_array_t3m1[] = {0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2};
  int row_s_array_t3m1[] = {0,3,6,9,12,15,18};

  construct_vector(val_array_t3m1,nnz_t3m1,val_t3m1);
  construct_vector(col_i_array_t3m1,nnz_t3m1,col_i_t3m1);
  construct_vector(row_s_array_t3m1,nrow_global_t3m1+1,row_s_t3m1);

  CRDoubleMatrix mat_t3m1;
  CRDoubleMatrixHelpers::create_uniformly_distributed_matrix(
      nrow_global_t3m1,ncol_t3m1,comm_pt,
      val_t3m1,col_i_t3m1,row_s_t3m1,mat_t3m1);

  std::ostringstream mat_t3m1_stream;
  mat_t3m1_stream << "t3m1_NP"<<nproc<<"R"<<my_rank;
  mat_t3m1.sparse_indexed_output(mat_t3m1_stream.str());
 
  // Second matrix (t3m2)
  unsigned nrow_global_t3m2 = 6;
  unsigned ncol_t3m2 = 3;
  unsigned nnz_t3m2 = 18;
  Vector<double> val_t3m2;
  Vector<int> col_i_t3m2;
  Vector<int> row_s_t3m2;

  double val_array_t3m2[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
  int col_i_array_t3m2[] = {0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2};
  int row_s_array_t3m2[] = {0,3,6,9,12,15,18};


  construct_vector(val_array_t3m2,nnz_t3m2,val_t3m2);
  construct_vector(col_i_array_t3m2,nnz_t3m2,col_i_t3m2);
  construct_vector(row_s_array_t3m2,nrow_global_t3m2+1,row_s_t3m2);

  CRDoubleMatrix mat_t3m2;
  CRDoubleMatrixHelpers::create_uniformly_distributed_matrix(
      nrow_global_t3m2,ncol_t3m2,comm_pt,
      val_t3m2,col_i_t3m2,row_s_t3m2,mat_t3m2);

  std::ostringstream mat_t3m2_stream;
  mat_t3m2_stream << "t3m2_NP"<<nproc<<"R"<<my_rank;
  mat_t3m2.sparse_indexed_output(mat_t3m2_stream.str());

  // Now perform element-wise addition of t3m1 and t3m2, and put the result
  // in mat_t3_result.
  CRDoubleMatrix mat_t3_result;

  mat_t3m1.add(mat_t3m2,mat_t3_result);

  std::ostringstream mat_t3_result_stream;
  mat_t3_result_stream << "t3_res_NP"<<nproc<<"R"<<my_rank;
  mat_t3_result.sparse_indexed_output(mat_t3_result_stream.str());


  // COLUMN INDICES NOT ORDERED. 
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
  // In this test, the column indices (per row) is not ordered.
  
  // First matrix (t4m1)
  unsigned nrow_global_t4m1 = 3;
  unsigned ncol_t4m1 = 6;
  unsigned nnz_t4m1 = 18;
  Vector<double> val_t4m1;
  Vector<int> col_i_t4m1;
  Vector<int> row_s_t4m1;

  double val_array_t4m1[] = {1,5,4,3,2,6,7,8,9,10,11,12,13,17,16,15,14,18};
  int col_i_array_t4m1[] = {0,4,3,2,1,5,0,1,2,3,4,5,0,4,3,2,1,5};
  int row_s_array_t4m1[] = {0,6,12,18};

  construct_vector(val_array_t4m1,nnz_t4m1,val_t4m1);
  construct_vector(col_i_array_t4m1,nnz_t4m1,col_i_t4m1);
  construct_vector(row_s_array_t4m1,nrow_global_t4m1+1,row_s_t4m1);

  CRDoubleMatrix mat_t4m1;
  CRDoubleMatrixHelpers::create_uniformly_distributed_matrix(
      nrow_global_t4m1,ncol_t4m1,comm_pt,
      val_t4m1,col_i_t4m1,row_s_t4m1,mat_t4m1);

  std::ostringstream mat_t4m1_stream;
  mat_t4m1_stream << "t4m1_NP"<<nproc<<"R"<<my_rank;
  mat_t4m1.sparse_indexed_output(mat_t4m1_stream.str());
 
  // Second matrix (t4m2)
  unsigned nrow_global_t4m2 = 3;
  unsigned ncol_t4m2 = 6;
  unsigned nnz_t4m2 = 18;
  Vector<double> val_t4m2;
  Vector<int> col_i_t4m2;
  Vector<int> row_s_t4m2;

  double val_array_t4m2[] = {1,2,3,4,5,6,7,11,10,9,8,12,13,14,15,16,17,18};
  int col_i_array_t4m2[] = {0,1,2,3,4,5,0,4,3,2,1,5,0,1,2,3,4,5};
  int row_s_array_t4m2[] = {0,6,12,18};

  construct_vector(val_array_t4m2,nnz_t4m2,val_t4m2);
  construct_vector(col_i_array_t4m2,nnz_t4m2,col_i_t4m2);
  construct_vector(row_s_array_t4m2,nrow_global_t4m2+1,row_s_t4m2);

  CRDoubleMatrix mat_t4m2;
  CRDoubleMatrixHelpers::create_uniformly_distributed_matrix(
      nrow_global_t4m2,ncol_t4m2,comm_pt,
      val_t4m2,col_i_t4m2,row_s_t4m2,mat_t4m2);

  std::ostringstream mat_t4m2_stream;
  mat_t4m2_stream << "t4m2_NP"<<nproc<<"R"<<my_rank;
  mat_t4m2.sparse_indexed_output(mat_t4m2_stream.str());

  // Now perform element-wise addition of t4m1 and t4m2, and put the result
  // in mat_t4_result.
  CRDoubleMatrix mat_t4_result;

  mat_t4m1.add(mat_t4m2,mat_t4_result);

  std::ostringstream mat_t4_result_stream;
  mat_t4_result_stream << "t4_res_NP"<<nproc<<"R"<<my_rank;
  mat_t4_result.sparse_indexed_output(mat_t4_result_stream.str());

#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main
