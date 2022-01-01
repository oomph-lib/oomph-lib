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

// Helper function to construct a Vector given an array. We use this because
// there is no better way to create an arbitrary Vector.
template<typename myType>
void construct_vector(const myType given_array[],
                      const unsigned given_arraysize,                  
                      Vector<myType> &result_vector)
{
  // Clear and reserve the required memory.
  result_vector.clear();

  std::copy(given_array,
            given_array + given_arraysize,
            std::back_inserter(result_vector));
}

// Helper function to output an oomph-lib Vector to oomph_info
template<typename myType>
void output_vector(const Vector<myType> &given_vector)
{
  typename Vector<myType>::const_iterator it;

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

  // my rank and number of processors. 
  // This is used later for putting the data.
  const unsigned my_rank = comm_pt->my_rank();
  const unsigned nproc = comm_pt->nproc();

  // Matrix 1
  // [1 2 0 0 0
  //  3 0 4 0 0
  //  0 0 0 5 6
  //  7 8 0 0 9
  //  0 0 0 0 0]
  //
  //  values = 1 2 3 4 5 6 7 8 9
  //  col i  = 0 1 0 2 3 4 0 1 4
  //  row s  = 0 2 4 6 9 9
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

  // Test the copy constructor
  CRDoubleMatrix copy_of_mat_t1m1(mat_t1m1);

  std::ostringstream copy_of_mat_t1m1_stream;
  copy_of_mat_t1m1_stream << "copy_of_t1m1_NP" << nproc << "R" << my_rank;
  copy_of_mat_t1m1.sparse_indexed_output(copy_of_mat_t1m1_stream.str());


  // Next matrix, testing non-square matrix where nrow < ncol.
  // A = [1  2  3  4  5  6
  //      7  8  9  10 11 12
  //      13 14 15 16 17 18]
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

  // Test the copy constructor
  CRDoubleMatrix copy_of_mat_t2m1(mat_t2m1);

  std::ostringstream copy_of_mat_t2m1_stream;
  copy_of_mat_t2m1_stream << "copy_of_t2m1_NP" << nproc << "R" << my_rank;
  copy_of_mat_t2m1.sparse_indexed_output(copy_of_mat_t2m1_stream.str());


  // Testing non square matrix where nrow > ncol
  // A = [1  2  3
  //      4  5  6
  //      7  8  9
  //      10 11 12
  //      13 14 15
  //      16 17 18]
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

  // Test the copy constructor
  CRDoubleMatrix copy_of_mat_t3m1(mat_t3m1);

  std::ostringstream copy_of_mat_t3m1_stream;
  copy_of_mat_t3m1_stream << "copy_of_t3m1_NP" << nproc << "R" << my_rank;
  copy_of_mat_t3m1.sparse_indexed_output(copy_of_mat_t3m1_stream.str());

#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main
