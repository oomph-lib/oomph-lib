//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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

// Given the number of blocks (in both dimensions), and an array containing
// the dimensions.
// Return: mat_info_vec, which contains the dimensions of the sub blocks 
// in a coherent manner. The number of global rows and columns of the block i,j
// is represented as: mat_info_vec[i][j][0] and mat_info_vec[i][j][1]
// respectively.
void fill_in_mat_info(unsigned const nblock_row, unsigned const nblock_col, 
                      unsigned dimarray[], 
                      Vector<Vector<Vector<unsigned> > > &mat_info_vec)
{
 // The total number of dimensions is two times the number of blocks
 // Since matrices are squares. This should be the same number of elements in
 // dimarray.
 unsigned ndims = 2*nblock_row*nblock_col;
  
 // Put dimarray into dimvec_all, since Vectors are nicer to work with.
 Vector<unsigned> dimvec_all;
 construct_vector(dimarray,ndims,dimvec_all);
  
 // index for the dimvec_all array.
 unsigned dimvec_all_i = 0;
 for (unsigned block_row_i = 0; block_row_i < nblock_row; block_row_i++) 
  {
   // create the Vector for the columns for the current row.
   Vector<Vector<unsigned> > current_block_col_vec;    
     
   // loop through the block columns
   for (unsigned block_col_i = 0; block_col_i < nblock_col; block_col_i++) 
    {
     // loop through the dimensions for this block.
     Vector<unsigned> current_dim_vec;
       
     for (unsigned dim_i = 0; dim_i < 2; dim_i++) 
      {
       current_dim_vec.push_back(dimvec_all[dimvec_all_i]);
       dimvec_all_i++;
      } // for the individual block dimensions.

     // push current_dim_vec into the current entry...
     current_block_col_vec.push_back(current_dim_vec);
    } // for the columns in this block row.

   // push the current block column vec onto the current row vec.
   mat_info_vec.push_back(current_block_col_vec);
  }// for the number of block rows.
}

// Given the dimensions of a matrix, returns a matrix where the elements 
// (starting from the number 1) are ascending per column, per row. 
// I.e., for a 3 by 3 matrix we have:
// [1 2 3
//  4 5 6
//  7 8 9]
void create_matrix_ascend_col_row(unsigned const nrow, unsigned const ncol,
                                  const OomphCommunicator* const comm_pt,
                                  bool const distributed, 
                                  CRDoubleMatrix &block)
{  
 // Clear the block
 block.clear();
  
 // Create the distribution.
 LinearAlgebraDistribution distri(comm_pt,nrow,distributed);

 // The number of rows this processor is responsible for.
 unsigned nrow_local = distri.nrow_local();
  
 // The number of values this processor will have to insert.
 unsigned nval = nrow_local*ncol;
  
 // The first_row will be used as an offset for the values to insert.
 unsigned first_row = distri.first_row();

 // Fill in values...
 Vector<double> values(nval,0);
 for (unsigned val_i = 0; val_i < nval; val_i++)
  {
   values[val_i] = (val_i+1) + (first_row)*ncol;
  }
  
 // Vectors for the column index and row_start.
 Vector<int> column_indicies(nval,0);
 Vector<int> row_start(nrow_local+1,0);

 unsigned column_indicies_i = 0;
 for (unsigned row_i = 0; row_i < nrow_local; row_i++)
  {
   for (unsigned col_i = 0; col_i < ncol; col_i++) 
    {
     column_indicies[column_indicies_i] = col_i;
     column_indicies_i++;
    }
   row_start[row_i] = ncol*row_i;
  }

 // Fill in the last row_start.
 row_start[nrow_local] = nval;

 block.build(&distri,ncol,values,column_indicies,row_start);
}


// Given mat_info_vec, which describes the sub blocks as so:
// The number of global rows and columns of the block i,j
// is represented as: mat_info_vec[i][j][0] and mat_info_vec[i][j][1]
// respectively. Pointers to the uniformly distribted matrices is returned.
void fill_in_sub_matrices(const OomphCommunicator* const comm_pt, 
                          const bool distributed,
                          Vector<Vector<Vector<unsigned> > > &mat_info,
                          DenseMatrix<CRDoubleMatrix* > &mat_pt)
{
 unsigned nblock_row = mat_pt.nrow();
 unsigned nblock_col = mat_pt.ncol();

 for (unsigned block_row_i = 0; block_row_i < nblock_row; block_row_i++) 
  {
   for (unsigned block_col_i = 0; block_col_i < nblock_col; block_col_i++) 
    {
     unsigned nrow = mat_info[block_row_i][block_col_i][0];
     unsigned ncol = mat_info[block_row_i][block_col_i][1];
      
     mat_pt(block_row_i,block_col_i)=new CRDoubleMatrix;

     create_matrix_ascend_col_row(nrow,ncol,comm_pt,distributed,
                                  *mat_pt(block_row_i,block_col_i));
    }
  }
}

// Helper funcion to call other helper functions
// 1) fill_in_mat_info(...) to fill in a data structure describing the matrices
// 2) fill_in_sub_matrices(...) uses the mat_info to fill in the sub matrices 
void create_matrices_to_cat
(unsigned dimarray[], const OomphCommunicator* const comm_pt, 
 DenseMatrix<CRDoubleMatrix*>&mat_pt)
{

 const unsigned nblock_row = mat_pt.nrow();
 const unsigned nblock_col = mat_pt.ncol();

 // We store the above array in a data structure such that
 // mat_info[0][0][0] gives us the number of rows in the block (0,0)
 // mat_info[0][0][1] gives us the number of columns in the block (0,0)
 Vector<Vector<Vector<unsigned> > > mat_info_vec;

 // Helper function to fill in mat_info given
 // nblock_row, nblock_col and dim_array.
 fill_in_mat_info(nblock_row,nblock_col,dimarray,mat_info_vec);
  
 // Fill in each sub matrix using create_matrix_ascend_col_row(..)
 bool distributed = true;
 fill_in_sub_matrices(comm_pt,distributed,mat_info_vec,mat_pt);
}

//===start_of_main======================================================
/// Driver code: Testing 
/// CRDoubleMatrixHelpers::concatenate_without_communication(...)
/// We concatenate uniformly distributed matrices.
/// 
/// Let (x,y) be a matrix with x rows and y columns, 
/// with entries increasing along the columns, then along the rows.
/// For example, (3,3) is
/// [1 2 3
///  4 5 6
///  7 8 9].
/// 
/// We concatenate the following matrices:
/// (7,7)(7,5)(7,3)
/// (5,7)(5,5)(5,3)
/// (3,7)(3,5)(3,3)
///
/// (7,7)(7,5)
/// (5,7)(5,5)
/// (3,7)(3,5)
/// 
/// (7,7)(7,5)(7,3)
/// (5,7)(5,5)(5,3)
///
/// Communication between processors is NOT required, but the block structure
/// of the sub blocks are NOT preserved in the result matrix. Please see
/// self_test/mpi/vector_concatenation_without_communication/ for more detail.
///
/// The script validate.sh should run this self test on 1, 2, 3 and 4 cores. 
//======================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
 // Initialise MPI
 MPI_Helpers::init(argc,argv);
#endif

 // Get the global oomph-lib communicator 
 const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

 unsigned nblock_row = 3;
 unsigned nblock_col = 3;

 // Supply the dimensions of the matrices to concatenate.
 // The matrices must be in the order: first row, then along the columns,
 // then second row, et cetera...
 // The number of elements in this must be 2*nblock_row*nblock_col.
 // In this test, we concatenate blocks with the following dimension and 
 // arrangements:
 //
 // (7,7) (7,5) (7,3)
 // (5,7) (5,5) (5,3)
 // (3,7) (3,5) (3,3)
 unsigned dimarray[] = {7,7,7,5,7,3,5,7,5,5,5,3,3,7,3,5,3,3};
  
 // The data structure to store the pointers to matrices.
 DenseMatrix<CRDoubleMatrix*> mat0_pt(nblock_row,nblock_col,0);

 // Create the matrice to concatenate.
 create_matrices_to_cat(dimarray,comm_pt,mat0_pt);

 /// //////////////////////////////////////////////////////////////////////////
  
 // stuff for the output
 unsigned my_rank = comm_pt->my_rank();
 unsigned nproc = comm_pt->nproc();
  
 // Test #0, we concaenate
 // (7,7) (7,5) (7,3)
 // (5,7) (5,5) (5,3)
 // (3,7) (3,5) (3,3)

 // We require the LinearAlgebraDistribution of the block rows.
 // We use the first column.
 Vector<LinearAlgebraDistribution*> row_distribution0_pt;
 for (unsigned block_row_i = 0; block_row_i < nblock_row; block_row_i++) 
  {
   row_distribution0_pt.push_back(mat0_pt(block_row_i,0)
                                  ->distribution_pt());
  }

 // Perform the concatenation.
 // Note: Because this is a square block matrix, we only need to pass one
 // Vector of distributions. This will be used for both row_distribution_pt
 // and col_distribution_pt in the concatenate_without_communication(...)
 // function.
 CRDoubleMatrix result_matrix0;
 CRDoubleMatrixHelpers::concatenate_without_communication(
  row_distribution0_pt,mat0_pt,result_matrix0);
 
 std::ostringstream result_mat0_stream;
 result_mat0_stream << "out0_NP" << nproc << "R" << my_rank;
 result_matrix0.sparse_indexed_output(result_mat0_stream.str().c_str());

 // Clear the result matrix.
 result_matrix0.clear();
  
 /// //////////////////////////////////////////////////////////////////////////
 // Test #1, we concaenate
 // (7,7) (7,5)
 // (5,7) (5,5)
 // (3,7) (3,5)

 nblock_row = 3;
 nblock_col = 2;
 DenseMatrix<CRDoubleMatrix*> mat1_pt(nblock_row,nblock_col,0);

 // Get the blocks.
 for (unsigned block_row_i = 0; block_row_i < nblock_row; block_row_i++) 
  {
   for (unsigned block_col_i = 0; block_col_i < nblock_col; block_col_i++) 
    {
     mat1_pt(block_row_i,block_col_i) = mat0_pt(block_row_i,block_col_i);
    }
  }

 // We require the LinearAlgebraDistribution of the block rows.
 // We use the first column.
 Vector<LinearAlgebraDistribution*> row_distribution1_pt;
 for (unsigned block_row_i = 0; block_row_i < nblock_row; block_row_i++) 
  {
   row_distribution1_pt.push_back(row_distribution0_pt[block_row_i]);
  }
 
 // We require the LinearAlgebraDistribution of the block "columns".
 // We use the first column.
 Vector<LinearAlgebraDistribution*> col_distribution1_pt;
 for (unsigned block_col_i = 0; block_col_i < nblock_col; block_col_i++) 
  {
   col_distribution1_pt.push_back(row_distribution0_pt[block_col_i]);
  }
 
 // The result matrix.
 CRDoubleMatrix result_matrix1;
  
 // Call the concatenate function.
 CRDoubleMatrixHelpers::concatenate_without_communication(
  row_distribution1_pt,col_distribution1_pt,mat1_pt,result_matrix1);
  
 std::ostringstream result_mat1_stream;
 result_mat1_stream << "out1_NP" << nproc << "R" << my_rank;
 result_matrix1.sparse_indexed_output(result_mat1_stream.str().c_str());

 // Clear the result matrix.
 result_matrix1.clear();
 
 /// //////////////////////////////////////////////////////////////////////////
 // Test #2, we concaenate
 // (7,7) (7,5) (7,3)
 // (5,7) (5,5) (5,3)
 nblock_row = 2;
 nblock_col = 3;
 DenseMatrix<CRDoubleMatrix*> mat2_pt(nblock_row,nblock_col,0);

 // Get the blocks.
 for (unsigned block_row_i = 0; block_row_i < nblock_row; block_row_i++) 
  {
   for (unsigned block_col_i = 0; block_col_i < nblock_col; block_col_i++) 
    {
     mat2_pt(block_row_i,block_col_i) = mat0_pt(block_row_i,block_col_i);
    }
  }

 // We require the LinearAlgebraDistribution of the block rows.
 // We use the first column.
 Vector<LinearAlgebraDistribution*> row_distribution2_pt;
 for (unsigned block_row_i = 0; block_row_i < nblock_row; block_row_i++) 
  {
   row_distribution2_pt.push_back(row_distribution0_pt[block_row_i]);
  }
 
 // We require the LinearAlgebraDistribution of the block "columns".
 // We use the first column.
 Vector<LinearAlgebraDistribution*> col_distribution2_pt;
 for (unsigned block_col_i = 0; block_col_i < nblock_col; block_col_i++) 
  {
   col_distribution2_pt.push_back(row_distribution0_pt[block_col_i]);
  }
 
 // The result matrix.
 CRDoubleMatrix result_matrix2;
  
 // Call the concatenate function.
 CRDoubleMatrixHelpers::concatenate_without_communication(
  row_distribution2_pt,col_distribution2_pt,mat2_pt,result_matrix2);
  
 std::ostringstream result_mat2_stream;
 result_mat2_stream << "out2_NP" << nproc << "R" << my_rank;
 result_matrix2.sparse_indexed_output(result_mat2_stream.str().c_str());

 // Clear the result matrix.
 result_matrix2.clear();

 // There are 3 by 3 block matrices to delete.
 nblock_row = 3;
 nblock_col = 3;
 // Delete the matrices.
 for (unsigned row_i = 0; row_i < nblock_row; row_i++) 
  {
   for (unsigned col_i = 0; col_i < nblock_col; col_i++) 
    {
     delete mat0_pt(row_i,col_i);
    } // for col_i
  } // for row_i

#ifdef OOMPH_HAS_MPI
 // finalize MPI
 MPI_Helpers::finalize();
#endif
 return(EXIT_SUCCESS);
} // end_of_main
