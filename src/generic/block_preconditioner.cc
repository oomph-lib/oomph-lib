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
#include "block_preconditioner.h"

namespace oomph
{

//=============================================================================
/// \short Gets block (i,j) from the original matrix and returns it in
/// block_matrix_pt (Specialisation for CCDoubleMatrix)
//=============================================================================
 template<> 
 void BlockPreconditioner<CCDoubleMatrix>:: 
 get_block(const unsigned& i, const unsigned& j, 
           CCDoubleMatrix* matrix_pt,
           CCDoubleMatrix*& block_pt)
 {

  // pointers for the jacobian matrix is compressed column sparse format 
  int* j_column_start;
  int* j_row_index;
  double* j_value;
  
  // sets pointers to jacobian matrix
  j_column_start = matrix_pt->column_start();
  j_row_index = matrix_pt->row_index();
  j_value = matrix_pt->value();

  // get the block dimensions
  unsigned block_nrow = this->block_dimension(i);
  unsigned block_ncol = this->block_dimension(j);

  // allocate temporary storage for the component vectors of block (i,j)
  // temp_ptr is used to point to an element in each column - required as
  // cannot assume that order of block's rows in jacobian and the block
  // matrix will be the same
  Vector<int> temp_row_index;
  Vector<int> temp_column_start(block_ncol+1);
  Vector<int> temp_ptr(block_ncol+1);
  Vector<double> temp_value;
  int block_nnz = 0;
  
  // get number of rows in source matrix
  unsigned master_nrow = this->master_nrow();
  
  // determine how many non zeros there are in the block (i,j)
  // also determines how many non zeros are stored in each row or column - 
  // stored in temp_ptr temporarily
  for (unsigned k = 0; k < master_nrow; k++)
   {
    if (block_number(k) == static_cast<int>(j))
     {
      for (int l = j_column_start[k]; 
           l < j_column_start[k+1]; l++)
       {
        if (block_number(j_row_index[l]) == 
            static_cast<int>(i))
         {
          block_nnz++;
          temp_ptr[index_in_block(k)+1]++;
         }
       }
     }
   }

  // if the block matrix is not empty
 if (block_nnz > 0)
   {

    // uses number of elements in each column of block to determine values
    // for the block column start (temp_column_start)
    temp_column_start[0] = 0;
    for (unsigned k = 1; k <= block_ncol; k++)
     {
      
      temp_column_start[k] = temp_column_start[k-1]+temp_ptr[k];
      temp_ptr[k] = temp_column_start[k];
     }
    
    // resizes the block row index and value to store all non zeros
    temp_row_index.resize(block_nnz);
    temp_value.resize(block_nnz);
    
    // copies the relevant elements of the jacobian to the correct entries 
    // of the block matrix
    for (unsigned k = 0; k < master_nrow; k++)
     {
      if (block_number(k) == static_cast<int>(j))
       {
        for (int l = j_column_start[k]; 
             l < j_column_start[k+1]; l++)
         {
          if (block_number(j_row_index[l]) == 
              static_cast<int>(i))
           {
            int kk = temp_ptr[index_in_block(k)]++;
            temp_value[kk] = j_value[l];
            temp_row_index[kk] = index_in_block(j_row_index[l]); 
           }
         }
       }
     }								
    
    // creates a new compressed column sparse matrix for the pointer for
    // the current block 
    block_pt = new CCDoubleMatrix(temp_value,temp_row_index,
                                  temp_column_start,block_nrow,block_ncol);

#ifdef PARANOID
    // checks to see if block matrix has been set up correctly 
    block_matrix_test(matrix_pt,i,j,block_pt);
#endif
   }
 
 // else the matrix is empty
 else
  {
   block_pt = 0;
  }
 }


//=============================================================================
/// \short Gets block (i,j) from the original matrix and returns it in
/// block_matrix_pt (Specialisation for CRDoubleMatrix)
//=============================================================================
 template<> 
 void BlockPreconditioner<CRDoubleMatrix>:: 
 get_block(const unsigned& i, const unsigned& j, CRDoubleMatrix* matrix_pt,
           CRDoubleMatrix*& block_pt)
 {

  // pointers for the jacobian matrix is compressed row sparse format 
  int* j_row_start;
  int* j_column_index;
  double* j_value;
  
  // sets pointers to jacobian matrix
  j_row_start = matrix_pt->row_start();
  j_column_index = matrix_pt->column_index();
  j_value = matrix_pt->value();

  // get the block dimensions
  unsigned block_nrow = this->block_dimension(i);
  unsigned block_ncol = this->block_dimension(j);

  // allocate temporary storage for the component vectors of block (i,j)
  // temp_ptr is used to point to an element in each column - required as
  // cannot assume that order of block's rows in jacobian and the block
  // matrix will be the same
  Vector<int> temp_column_index;
  Vector<int> temp_row_start(block_nrow+1);
  Vector<int> temp_ptr(block_nrow+1);
  Vector<double> temp_value;
  int block_nnz = 0;
  
  // get number of rows in source matrix
  unsigned master_nrow = this->master_nrow();
  
  // determine how many non zeros there are in the block (i,j)
  // also determines how many non zeros are stored in each row or column - 
  // stored in temp_ptr temporarily
  for (unsigned k = 0; k < master_nrow; k++)
   {
    if (block_number(k) == static_cast<int>(i))
     {
      for (int l = j_row_start[k]; 
           l < j_row_start[k+1]; l++)
       {
        if (block_number(j_column_index[l]) == 
            static_cast<int>(j))
         {
          block_nnz++;
          temp_ptr[index_in_block(k)+1]++;
         }
       }
     }
   }
  
  // if the matrix is not empty
  if (block_nnz > 0)
   {
    
    // uses number of elements in each column of block to determine values
    // for the block column start (temp_row_start)
    temp_row_start[0] = 0;
    for (unsigned k = 1; k <= block_nrow; k++)
     {
      temp_row_start[k] = temp_row_start[k-1]+temp_ptr[k];
      temp_ptr[k] = temp_row_start[k];
     }
  
    // resizes the block row index and value to store all non zeros
    temp_column_index.resize(block_nnz);
    temp_value.resize(block_nnz);
    
    // copies the relevant elements of the jacobian to the correct entries 
    // of the block matrix
    for (unsigned k = 0; k < master_nrow; k++)
     {
      if (block_number(k) == static_cast<int>(i))
       {
        for (int l = j_row_start[k]; 
             l < j_row_start[k+1]; l++)
         {
          if (block_number(j_column_index[l]) == 
              static_cast<int>(j))
           {
            int kk = temp_ptr[index_in_block(k)]++;
            temp_value[kk] = j_value[l];
            temp_column_index[kk] = 
             index_in_block(j_column_index[l]); 
           }
         }
       }
     }								
    
   
    // creates a new compressed column sparse matrix for the pointer for
    // the current block 
    block_pt = new CRDoubleMatrix(temp_value,temp_column_index,
                                  temp_row_start,block_nrow,block_ncol);
    
#ifdef PARANOID
    // checks to see if block matrix has been set up correctly 
    block_matrix_test(matrix_pt,i,j,block_pt);
#endif
    
   }
  // else the matrix is empty
  else
   {
    block_pt = 0;
   }
 }


//=============================================================================
/// \short test function to check that every element in the block matrix
/// (block_i,block_j) matches the corresponding element in the original matrix
//=============================================================================
 template<typename MATRIX>
 void BlockPreconditioner<MATRIX>::block_matrix_test(const MATRIX* matrix_pt,
                                                     const unsigned& block_i,
                                                     const unsigned& block_j,
                                                     const MATRIX*
                                                     block_matrix_pt)
 {

  // boolean flag to indicate whether test is passed
  bool check = true;
  
  // number of rows in matrix
  unsigned n_row = matrix_pt->nrow();
  
  // number of columns in matrix
  unsigned n_col = matrix_pt->ncol();
  
  // loop over rows of original matrix
  for (unsigned i = 0; i < n_row; i++)
   {
    
    // if this coefficient is associated with a block in this block 
    // preconditioner
    if (static_cast<int>(block_i) == this->block_number(i))
     {
      
      // loop over columns of original matrix
      for (unsigned j = 0; j < n_col; j++)
       {
        
        // if the coeeficient is associated with a block in this block
        // preconditioner
        if (static_cast<int>(block_j) == this->block_number(j))
         {
          
          // check whether elements in original matrix and matrix of block 
          // pointers match
          if ( matrix_pt->operator()(i,j) !=
               block_matrix_pt
               ->operator()(index_in_block(i),index_in_block(j)) )
           {
            check = false;
           }
         }
       }
     }
   }
  
  // throw error
  if (!check)
   {
    std::ostringstream error_message;
    error_message << "The require elements have not been successfully copied"
                  << " from the original matrix to the block matrices";
    throw OomphLibError(error_message.str(),
                        "BlockPreconditioner<MATRIX>::block_matrix_test()",
                        OOMPH_EXCEPTION_LOCATION);
   }
 }

//=============================================================================
/// Assemble the block preconditioner as a single matrix. 
/// In some cases the block preconditioner cannot be applied in 
/// individual components; this function takes the matrix
/// of block pointers and returns a single matrix containing all the 
/// blocks of the matrix of blocks in a single matrix that can
/// be solved directly. Specialised version for CCDoubleMatrix.
//=============================================================================
 template<> 
 void BlockPreconditioner<CCDoubleMatrix>::build_preconditioner_matrix( 
  DenseMatrix<CCDoubleMatrix*>& block_matrix_pt, CCDoubleMatrix*& 
  preconditioner_matrix_pt)
 {

  // number of non zeros in the final preconditioner matrix
  int p_nnz = 0;
  for (unsigned j = 0; j < Nblock_types; j++)
   {
    for (unsigned i = 0; i < Nblock_types; i++)
     {
      if (block_matrix_pt(i,j) != 0)
       {
        p_nnz += block_matrix_pt(i,j)->nnz();	
       }
     }
   }

  // vector indicating the offset required for a block
  Vector<unsigned> block_start(Nblock_types,0);
  for (unsigned a = 1; a < Nblock_types; a++)
   {
    block_start[a] = block_start[a-1]+block_dimension(a-1);
   }

  // determine p matrix n_row
  unsigned p_nrow = Nrow;

  // temporary storage for complete block preconditioner matrix
  Vector<int> p_column_start(p_nrow+1);
  Vector<int> p_row_index(p_nnz);
  Vector<double> p_value(p_nnz);
  int p_ptr = 0;
  

  p_column_start[0] = 0;
  
  // loops of the block columns 
  for (unsigned j = 0; j < Nblock_types; j++)
   {
    
    // determines the block column offset
    int j_block_start = block_start[j];
    
    // loop over the columns of the current block
    for (unsigned k = 0; k < block_dimension(j); k++)
     {
      
      //
      p_column_start[j_block_start + k+1] = p_column_start[j_block_start+k];
      
      // loop over the block rows
      for (unsigned i = 0; i < Nblock_types; i++)
       {
        
        // if block(i,j) pointer no null then
        if (block_matrix_pt(i,j) != 0)
         {
        
          // sets the block row offset
          unsigned i_block_start = block_start[i];			
          
          // creates pointers for the elements in the current block
          int* temp_column_start= block_matrix_pt(i,j)->column_start();
          int* temp_row_index = block_matrix_pt(i,j)->row_index();
          double* temp_value = block_matrix_pt(i,j)->value();
          
          // sets the next column start
          p_column_start[j_block_start + k + 1] += temp_column_start[k+1] -
           temp_column_start[k];
          
          // adds of the current row of the current block to the preconditioner
          // matrix
          for (int l = temp_column_start[k]; l < temp_column_start[k+1]; l++)
           {
            p_row_index[p_ptr] = temp_row_index[l] + i_block_start;
            p_value[p_ptr] = temp_value[l];
            p_ptr++;
           }
         }	
       }			
     }
   }

  // builds the preconditioner matrix 
  // ALH:: Assumed to be square
  preconditioner_matrix_pt = new CCDoubleMatrix(p_value, p_row_index, 
                                              p_column_start,p_nrow,p_nrow);
 }
   


//=============================================================================
/// Assemble the block preconditioner as a single matrix. 
/// In some cases the block preconditioner cannot be applied in 
/// individual components; this function takes the matrix
/// of block pointers and returns a single matrix containing all the 
/// blocks of the matrix of blocks in a single matrix that can
/// be solved directly. Specialised version for CRDoubleMatrix.
//=============================================================================
 template<>
 void BlockPreconditioner<CRDoubleMatrix>::build_preconditioner_matrix(
  DenseMatrix<CRDoubleMatrix*>& block_matrix_pt, CRDoubleMatrix*&
  preconditioner_matrix_pt)
 {

  // determine the minimum and maximum block type number for the blocks to be 
  // built into a single matrix
  unsigned p_nnz = 0;
  for (unsigned i = 0; i < Nblock_types; i++)
   {
    for (unsigned j = 0; j < Nblock_types; j++)
     {
      if (block_matrix_pt(i,j) !=0)
       {
        p_nnz += block_matrix_pt(i,j)->nnz();
       }
     }
   }

  // vector indicating the offset required for a block
  Vector<unsigned> block_start(Nblock_types,0);
  for (unsigned a = 1; a < Nblock_types; a++)
   {
    block_start[a] = block_start[a-1]+this->block_dimension(a-1);
   }
  
  // determine p matrix n_row
  unsigned p_nrow = Nrow;

  // temporary storage for complete block preconditioner matrix
  Vector<int> p_row_start(p_nrow+1); 
  Vector<int> p_column_index;
  Vector<double> p_value;
  
  p_row_start[0] = 0;
  
  // loops of the block columns 
  for (unsigned i = 0; i < Nblock_types; i++)
   {
    
    // determines the block row offset
    int i_block_start = block_start[i];
    
    // loop over the rows of the current block
    for (unsigned k = 0; k < this->block_dimension(i); k++)
     {
      
      //
      p_row_start[i_block_start + k+1] = p_row_start[i_block_start+k];
      
      // loop over the block rows
      for (unsigned j = 0; j < Nblock_types; j++)
       {
        
        // if block(i,j) pointer not null then
        if (block_matrix_pt(i,j) != 0)
         {
          // sets the block row offset
          unsigned j_block_start = block_start[j];
          
          // creates pointers for the elements in the current block
          int* temp_row_start = block_matrix_pt(i,j)->row_start();
          int* temp_column_index =block_matrix_pt(i,j)->
           column_index();
          double* temp_value = block_matrix_pt(i,j)->value();
          
          // adds of the current row of the current block to the preconditioner
          // matrix
          for (int l = temp_row_start[k]; l < temp_row_start[k+1]; l++)
           {
            p_row_start[i_block_start+k+1]++;
            p_column_index.push_back(temp_column_index[l] + j_block_start);
            p_value.push_back(temp_value[l]);
           }
         }	
       }			
     }
   }
  
  // builds the preconditioner matrix
  // ALH: Assumed to be square
  preconditioner_matrix_pt = new CRDoubleMatrix(p_value, p_column_index, 
                                                p_row_start,
                                                p_nrow,p_nrow);
 }
}

