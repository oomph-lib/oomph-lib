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

 /// \short Static boolean to allow block_matrix_test(...) to be run.
 /// Defaults to false.
 template<typename MATRIX> 
 bool BlockPreconditioner<MATRIX>::Run_block_matrix_test=false;
 
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
    if (Run_block_matrix_test)
     {
      // checks to see if block matrix has been set up correctly 
      block_matrix_test(matrix_pt,i,j,block_pt);
     }
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
 get_block(const unsigned& block_i, const unsigned& block_j, 
           CRDoubleMatrix* matrix_pt,
           CRDoubleMatrix*& block_pt)
 {

#ifdef PARANOID
  // the number of blocks
  unsigned n_blocks = this->nblock_types();

  // paranoid check that block i is in this block preconditioner
  if (block_i >= n_blocks || block_j >= n_blocks)
   {
    std::ostringstream error_message;
    error_message << "Requested block (" << block_i << "," << block_j   
                  << "), however this preconditioner has nblock_types() "
                  << "= " << nblock_types() << std::endl;
    throw OomphLibError(error_message.str(),
                        "BlockPreconditioner::get_block(...)",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // if + only one processor
  //    + more than one processor but matrix_pt is not distributed
  // then use the serial get_block method
  if (matrix_pt->distribution_pt()->communicator_pt()->nproc() == 1 ||
      !matrix_pt->distribution_pt()->distributed())
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
    unsigned block_nrow = this->block_dimension(block_i);
    unsigned block_ncol = this->block_dimension(block_j);
    
    // allocate temporary storage for the component vectors of block (i,j)
    // temp_ptr is used to point to an element in each column - required as
    // cannot assume that order of block's rows in jacobian and the block
    // matrix will be the same
    int* temp_row_start = new int[block_nrow+1];
    for (unsigned i = 0; i <= block_nrow; i++)
     {
      temp_row_start[i] = 0;
     }
    Vector<int> temp_ptr(block_nrow+1);
    int block_nnz = 0;
    
    // get number of rows in source matrix
    unsigned master_nrow = this->master_nrow();
    
    // determine how many non zeros there are in the block (i,j)
    // also determines how many non zeros are stored in each row or column - 
    // stored in temp_ptr temporarily
    for (unsigned k = 0; k < master_nrow; k++)
     {
      if (block_number(k) == static_cast<int>(block_i))
       {
        for (int l = j_row_start[k]; 
             l < j_row_start[k+1]; l++)
         {
          if (block_number(j_column_index[l]) == 
              static_cast<int>(block_j))
           {
            block_nnz++;
            temp_ptr[index_in_block(k)+1]++;
           }
         }
       }
     }
    
    // if the matrix is not empty
    int* temp_column_index = new int[block_nnz];
    double* temp_value = new double[block_nnz];
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
      
      // copies the relevant elements of the jacobian to the correct entries 
      // of the block matrix
      for (unsigned k = 0; k < master_nrow; k++)
       {
        if (block_number(k) == static_cast<int>(block_i))
         {
          for (int l = j_row_start[k]; 
               l < j_row_start[k+1]; l++)
           {
            if (block_number(j_column_index[l]) == 
                static_cast<int>(block_j))
             {
              int kk = temp_ptr[index_in_block(k)]++;
              temp_value[kk] = j_value[l];
              temp_column_index[kk] = 
               index_in_block(j_column_index[l]); 
             }
           }
         }
       }
     }
      
      
    // creates a new compressed column sparse matrix for the pointer for
    // the current block 
    if (block_pt == 0)
     {
      block_pt = new CRDoubleMatrix(Block_distribution_pt[block_i]);
     }
    block_pt->rebuild_matrix_without_copy(block_ncol,block_nnz,
                                          temp_value,temp_column_index,
                                          temp_row_start);
 
#ifdef PARANOID
    if (Run_block_matrix_test)
     {
      // checks to see if block matrix has been set up correctly 
      block_matrix_test(matrix_pt,block_i,block_j,block_pt);
     }
#endif 
   }


  // otherwise we are dealing with a distributed matrix
  else
   {
#ifdef OOMPH_HAS_MPI
    // number of processors
    unsigned nproc = Distribution_pt->communicator_pt()->nproc();

    // my rank
    unsigned my_rank = Distribution_pt->communicator_pt()->my_rank();

    // sets pointers to jacobian matrix
    int* j_row_start = matrix_pt->row_start();
    int* j_column_index = matrix_pt->column_index();
    double* j_value = matrix_pt->value();

    // number of non zeros in each row to be sent
    Vector<int*> nnz_send(nproc,0);

    // number of non zeros in each row to be received
    Vector<int*> nnz_recv(nproc,0);

    // storage for data to be sent
    Vector<int*> column_index_for_proc(nproc,0);
    Vector<double*> values_for_proc(nproc,0);

    // number of non zeros to be sent to each processor
    Vector<unsigned> total_nnz_send(nproc,0);

    // number of rows of the block matrix on this processor
    unsigned nrow_local = Block_distribution_pt[block_i]->nrow_local();

    // resize the nnz storage and compute nnz_send
    // and send and recv the nnz
    Vector<MPI_Request> send_req;
    Vector<MPI_Request> recv1_req;
    for (unsigned p = 0; p < nproc; p++)
     {
      int nrow_send = Nrows_to_send_for_get_block(block_i,p);
      int nrow_recv = Nrows_to_recv_for_get_block(block_i,p);

      // assemble nnz recv
      nnz_recv[p] = new int[nrow_recv];

      // assemble the storage to send
      if (nrow_send > 0 && p != my_rank)
       {
        nnz_send[p] = new int[nrow_send];
       }

      // compute the number of nnzs in each row and the total number
      // of nnzs
      for (int i = 0; i < nrow_send; i++)
       {
        unsigned row = Rows_to_send_for_get_block(block_i,p)[i];
        int c = 0;
        for (int r = j_row_start[row]; r < j_row_start[row+1]; r++)
         {
          if (block_number(j_column_index[r]) == int(block_j))
           {
            c++;
           }
         }
        if (p != my_rank)
         {
          nnz_send[p][i] = c;
         }
        else
         {
          nnz_recv[p][i] = c;
         }
        total_nnz_send[p] += c;
       }

      // send
      if (p != my_rank)
       {
        if (nrow_send)
         {
          MPI_Request req;
          MPI_Isend(nnz_send[p],nrow_send,MPI_INT,p,0,
                    Distribution_pt->communicator_pt()->mpi_comm(),&req);
          send_req.push_back(req);
         }

        // recv
        if (nrow_recv)
         {
          MPI_Request req;
          MPI_Irecv(nnz_recv[p],nrow_recv,MPI_INT,p,0,
                    Distribution_pt->communicator_pt()->mpi_comm(),&req);
          recv1_req.push_back(req);
         }
       }
     }

    // next assemble the values and row_start data to be sent for each
    // processor
    for (unsigned p = 0; p < nproc; p++)
     {
      int nrow_send = Nrows_to_send_for_get_block(block_i,p);

      // assemble the storage for the values and column indices to be sent
      if (p != my_rank)
       {
        if (total_nnz_send[p] > 0)
         {
          values_for_proc[p] = new double[total_nnz_send[p]];
          column_index_for_proc[p] = new int[total_nnz_send[p]];
          
          // copy the values and column indices to the storage
          unsigned ptr = 0;
          for (int i = 0; i < nrow_send; i++)
           {
            unsigned row = Rows_to_send_for_get_block(block_i,p)[i];
            for (int r = j_row_start[row]; r < j_row_start[row+1]; r++)
             {
              if (block_number(j_column_index[r]) == int(block_j))
               {
                values_for_proc[p][ptr] = j_value[r];
                column_index_for_proc[p][ptr] = 
                 index_in_block(j_column_index[r]);
                ptr++;
               }
             }
           }
       
          // create the datatypes
          MPI_Datatype types[2];
          MPI_Type_contiguous(total_nnz_send[p],MPI_DOUBLE,&types[0]);
          MPI_Type_commit(&types[0]);
          MPI_Type_contiguous(total_nnz_send[p],MPI_INT,&types[1]);
          MPI_Type_commit(&types[1]);
          
          // get the start address of the vectors
          MPI_Aint displacement[2];
          MPI_Address(values_for_proc[p],&displacement[0]);
          MPI_Address(column_index_for_proc[p],&displacement[1]);
          
          // compute the displacements
          displacement[1] -= displacement[0];
          displacement[0] -= displacement[0];

          // compute the block lengths
          int length[2];
          length[0] = length[1] = 1;

          // build the struct data type
          MPI_Datatype final_type;
          MPI_Type_struct(2,length,displacement,types,&final_type);
          MPI_Type_commit(&final_type);
          MPI_Type_free(&types[0]);
          MPI_Type_free(&types[1]);

          // and send
          MPI_Request req;
          MPI_Isend(values_for_proc[p],1,final_type,p,1,
                    Distribution_pt->communicator_pt()->mpi_comm(),&req);
          send_req.push_back(req);
          MPI_Type_free(&final_type);
         }
       }
     }

    // wait for the recv to complete (the row_start recv which actually
    // contains the number of nnzs in each row)
    int c_recv = recv1_req.size();
    if (c_recv != 0)
     {
      MPI_Waitall(c_recv,&recv1_req[0],MPI_STATUS_IGNORE);
     }

    // compute the total number of nnzs to be received
    Vector<int> total_nnz_recv_from_proc(nproc);
    int local_block_nnz = 0;
    for (unsigned p = 0; p < nproc; p++)
     {
      // compute the total nnzs
      for (unsigned i = 0; i < Nrows_to_recv_for_get_block(block_i,p); i++)
       {
        total_nnz_recv_from_proc[p] += nnz_recv[p][i];

       }
      local_block_nnz += total_nnz_recv_from_proc[p];
     }

    // compute the offset for each block of nnzs (a matrix row) in the 
    // values_recv and column_index_recv vectors

    // fisrt determine how many blocks of rows are to be recv
    Vector<int> n_recv_block(nproc,0);
    for (unsigned p = 0; p < nproc; p++)
     {
      if (Nrows_to_recv_for_get_block(block_i,p) > 0)
       {
        n_recv_block[p] = 1;
       }
      for (unsigned i = 1; i < Nrows_to_recv_for_get_block(block_i,p); i++)
       {
        if (Rows_to_recv_for_get_block(block_i,p)[i] !=
            Rows_to_recv_for_get_block(block_i,p)[i-1] + 1)
         {
          n_recv_block[p]++;
         }
       }
     }

    // next assemble row start recv
    int* row_start_recv = new int[nrow_local+1];
    for (unsigned i = 0; i <= nrow_local; i++)
     {
      row_start_recv[i] = 0;
     }
    for (unsigned p = 0; p < nproc; p++)
     {
      for (unsigned i = 0; i < Nrows_to_recv_for_get_block(block_i,p); i++)   
       {
        row_start_recv[Rows_to_recv_for_get_block(block_i,p)[i]] 
         = nnz_recv[p][i];
       }
     }
    int g = row_start_recv[0];
    row_start_recv[0] = 0;
    for (unsigned i = 1; i < nrow_local; i++)
     {
      int temp_g = g;
      g = row_start_recv[i];
      row_start_recv[i] = row_start_recv[i-1] + temp_g;
     }    
    row_start_recv[nrow_local] = row_start_recv[nrow_local-1] + g;

    // next assemble the offset and the number of nzs in each recv block
    Vector<int*> offset_recv_block(nproc,0);
    Vector<int*> nnz_recv_block(nproc,0);
    for (unsigned p = 0; p < nproc; p++)
     {
      if (Nrows_to_recv_for_get_block(block_i,p) > 0)
       {
        offset_recv_block[p] = new int[n_recv_block[p]];
        offset_recv_block[p][0] = 0;
        nnz_recv_block[p] = new int[n_recv_block[p]];
        for (int i = 0; i < n_recv_block[p]; i++)
         {
          nnz_recv_block[p][i] = 0;
         }
        unsigned ptr = 0;
        nnz_recv_block[p][ptr] += nnz_recv[p][0];
        offset_recv_block[p][0] 
         = row_start_recv[Rows_to_recv_for_get_block(block_i,p)[0]];
        for (unsigned i = 1; i < Nrows_to_recv_for_get_block(block_i,p); i++)
         {
          if (Rows_to_recv_for_get_block(block_i,p)[i] !=
              Rows_to_recv_for_get_block(block_i,p)[i-1] + 1)
           {
            ptr++;
            offset_recv_block[p][ptr] 
             = row_start_recv[Rows_to_recv_for_get_block(block_i,p)[i]];
           }
          nnz_recv_block[p][ptr] += nnz_recv[p][i];
         }
       }
      delete[] nnz_recv[p];
     }

    // post the receives
    int* column_index_recv = new int[local_block_nnz];
    double* values_recv = new double[local_block_nnz];
    Vector<MPI_Request> recv2_req;
    for (unsigned p = 0; p < nproc; p++)
     {
      if (p != my_rank)
       {
        if (total_nnz_recv_from_proc[p] != 0)
         {
          // create the datatypes
          MPI_Datatype types[2];
          MPI_Type_indexed(n_recv_block[p],nnz_recv_block[p],         
                           offset_recv_block[p],MPI_DOUBLE,&types[0]);
          MPI_Type_commit(&types[0]);
          MPI_Type_indexed(n_recv_block[p],nnz_recv_block[p],         
                           offset_recv_block[p],MPI_INT,&types[1]);
          MPI_Type_commit(&types[1]);
          
          // compute the displacements
          MPI_Aint displacements[2];
          MPI_Address(values_recv, &displacements[0]);
          MPI_Address(column_index_recv,&displacements[1]);
          displacements[1] -= displacements[0];
          displacements[0] -= displacements[0];
          
          // compute the block lengths
          int length[2];
          length[0] = length[1] = 1;
          
          // create the final datatype
          MPI_Datatype final_type;
          MPI_Type_struct(2,length,displacements,types,&final_type);
          MPI_Type_commit(&final_type);
	  MPI_Type_free(&types[0]);
	  MPI_Type_free(&types[1]);
          
          // and the recv
          MPI_Request req;
          MPI_Irecv(values_recv,1,final_type,p,1,
                    Distribution_pt->communicator_pt()->mpi_comm(),&req);
          recv2_req.push_back(req);
	  MPI_Type_free(&final_type);
         }
       }
      else
       {
        // next send the values and column indices to self
        unsigned block_ptr = 0;
        unsigned counter = 0;
        int nrow_send = Nrows_to_send_for_get_block(block_i,my_rank);
        if (nrow_send > 0)
         {
          unsigned offset = offset_recv_block[my_rank][0];
          for (int i = 0; i < nrow_send; i++)
           {
            if (i > 0)
             {
              if (Rows_to_recv_for_get_block(block_i,p)[i] !=
                  Rows_to_recv_for_get_block(block_i,p)[i-1] + 1)
               {
                counter = 0;
                block_ptr++;
                offset = offset_recv_block[my_rank][block_ptr];
               }
             }
            unsigned row = Rows_to_send_for_get_block(block_i,my_rank)[i];
            for (int r = j_row_start[row]; r < j_row_start[row+1]; r++)
             {
              if (block_number(j_column_index[r]) == int(block_j))
               {
                values_recv[offset+counter] = j_value[r];
                column_index_recv[offset + counter] = 
                 index_in_block(j_column_index[r]);
                counter++;
               }
             }
           }
         }
       }
     }
       
    // wait for the recv to complete (for the column_index and the values_
    c_recv = recv2_req.size();
    if (c_recv != 0)
     {
      MPI_Waitall(c_recv,&recv2_req[0],MPI_STATUS_IGNORE);   
     }

    // create the matrix
    // creates a new compressed column sparse matrix for the pointer for
    // the current block 
    if (block_pt == 0)
     {
      block_pt = new CRDoubleMatrix(Block_distribution_pt[block_i]);
     }
    block_pt->rebuild_matrix_without_copy(this->block_dimension(block_j),
                                          local_block_nnz,
                                          values_recv,column_index_recv,
                                          row_start_recv);

    // wait for the send to complete (nnz / row_start)
    int c_send = send_req.size();
    if (c_send)
     {
      MPI_Waitall(c_send,&send_req[0],MPI_STATUS_IGNORE);
     }

    // delete temp storage used for assembling data for communication
    for (unsigned p = 0; p < nproc; p++)
     {
      delete[] nnz_send[p];
      delete[] column_index_for_proc[p];
      delete[] values_for_proc[p];
      delete[] offset_recv_block[p];
      delete[] nnz_recv_block[p];
     }
#else
    // throw error
    std::ostringstream error_message;
    error_message << "The matrix is distributed and on more than one "
                  << "processor. MPI is required.";
    throw OomphLibError(error_message.str(),
                        "BlockPreconditioner<MATRIX>::get_block(...)",
                        OOMPH_EXCEPTION_LOCATION);
#endif 
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
    
  // creates a new compressed column sparse matrix for the pointer for
  // the current block 
  if (preconditioner_matrix_pt != 0)
   {
    delete preconditioner_matrix_pt;
   }
  preconditioner_matrix_pt = new CCDoubleMatrix(p_value, 
						p_row_index, 
						p_column_start,
						p_nrow,p_nrow);
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
#ifdef PARANOID
  // the number of blocks
  unsigned n_blocks = this->nblock_types();

  // paranoid check that block i is in this block preconditioner
  if (block_matrix_pt.nrow() != n_blocks || block_matrix_pt.ncol() != n_blocks)
   {
    std::ostringstream error_message;
    error_message << "Block_matrix_pt must be of dimension nblock x nblock "
		  << "(" << n_blocks << " x " << n_blocks << "), but it is "
		  << block_matrix_pt.nrow() << " x " << block_matrix_pt.ncol()
		  << std::endl;
    throw OomphLibError(error_message.str(),
                        "BlockPreconditioner::get_block(...)",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // if + only one processor
  //    + more than one processor but matrix_pt is not distributed
  // then use the serial get_block method
  if (Problem_pt->communicator_pt()->nproc() == 1 ||
      !Distribution_pt->distributed())
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
    int* p_row_start = new int[p_nrow+1]; 
    int* p_column_index = new int[p_nnz];
    double* p_value = new double[p_nnz];
  
    p_row_start[0] = 0;
  
    // loops of the block columns 
    unsigned p_coef_pt = 0;
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
              p_column_index[p_coef_pt] = temp_column_index[l] + j_block_start;
              p_value[p_coef_pt] = temp_value[l];
              p_coef_pt++;
             }
           }	
         }			
       }
     }

    // creates a new compressed column sparse matrix for the pointer for
    // the current block 
    if (preconditioner_matrix_pt == 0)
     {
      preconditioner_matrix_pt = new CRDoubleMatrix(Distribution_pt);
     }
    preconditioner_matrix_pt->rebuild_matrix_without_copy(p_nrow,p_nnz,
                                                          p_value,
                                                          p_column_index, 
                                                          p_row_start);
   }
  else
   {
#ifdef OOMPH_HAS_MPI

    // my rank
    unsigned my_rank = Distribution_pt->communicator_pt()->my_rank();
    
    // the number of processors
    unsigned nproc = Distribution_pt->communicator_pt()->nproc();

    // number of blocks
    unsigned nblock = this->Nblock_types;

    // begin by determining which rows of which matrix should be sent
    // and recv
    Vector< Vector<int> > first_row_to_send(nblock);
    Vector< Vector<int> > nrow_to_send(nblock);
    Vector< Vector<int> > first_row_to_recv(nblock);
    Vector< Vector<int> > nrow_to_recv(nblock);
    unsigned my_block_offset = 0;
    for (unsigned b = 0; b < nblock; b++)
     {
      first_row_to_send[b].resize(nproc,0);
      nrow_to_send[b].resize(nproc,0);
      first_row_to_recv[b].resize(nproc,0);
      nrow_to_recv[b].resize(nproc,0);
      for (unsigned p = 0; p < nproc; p++)
       {
        if ((Distribution_pt->first_row(p) < 
             (my_block_offset + Block_distribution_pt[b]->first_row(my_rank) +
              Block_distribution_pt[b]->nrow_local(my_rank))) &&
            (my_block_offset + Block_distribution_pt[b]->first_row(my_rank) <
             (Distribution_pt->first_row(p) + Distribution_pt->nrow_local(p))))
         {
          first_row_to_send[b][p] = 
           std::max(Block_distribution_pt[b]->first_row(my_rank)
                    + my_block_offset,
                    Distribution_pt->first_row(p))-my_block_offset;
          nrow_to_send[b][p] = 
           std::min(Block_distribution_pt[b]->first_row(my_rank) + 
                    Block_distribution_pt[b]->nrow_local(my_rank),
                    Distribution_pt->first_row(p) + Distribution_pt->nrow_local(p)
                    - my_block_offset) - first_row_to_send[b][p];
         }

        if ((Distribution_pt->first_row(my_rank) < 
             (Block_distribution_pt[b]->first_row(p) +
              Block_distribution_pt[b]->nrow_local(p) + 
              my_block_offset)) &&
            (Block_distribution_pt[b]->first_row(p) + my_block_offset <
             (Distribution_pt->first_row(my_rank) +
              Distribution_pt->nrow_local(my_rank))))
         {
          first_row_to_recv[b][p] = 
           std::max(Block_distribution_pt[b]->first_row(p) + my_block_offset,
                    Distribution_pt->first_row(my_rank)) - my_block_offset;
          nrow_to_recv[b][p] = 
           std::min(Block_distribution_pt[b]->first_row(p) + 
                    Block_distribution_pt[b]->nrow_local(p),
                    Distribution_pt->first_row(my_rank) + 
                    Distribution_pt->nrow_local(my_rank) - my_block_offset) -
           first_row_to_recv[b][p];
         }
       }
      my_block_offset += this->block_dimension(b);
     }

    // assemble the number of coefs in each row to be sent to each proc
    // and post the sends and recvs
    Vector<int*> nnz_to_be_sent(nproc);
    Vector<int*> nnz_to_be_recv(nproc);
    Vector<MPI_Request> send_reqs_nnz;
    Vector<MPI_Request> recv_reqs_nnz;
    Vector<unsigned> nrow_send_p(nproc);
    Vector<unsigned> nrow_recv_p(nproc);
    Vector<unsigned> nmatrix_send_p(nproc,0);
    for (unsigned p = 0; p < nproc; p++)
     {

      // compute the number of rows to send
      unsigned nrow_send = 0;
      unsigned nrow_recv = 0;
      unsigned nmatrix_send = 0;
      for (unsigned ib = 0; ib < nblock; ib++)
       {
        for (unsigned jb = 0; jb < nblock; jb++)
         {
          if (block_matrix_pt(ib,jb) != 0)
           {
            nrow_recv += nrow_to_recv[ib][p];
            nrow_send += nrow_to_send[ib][p];
	    if (nrow_to_send[ib][p] > 0)
             {
              nmatrix_send++;
             }
           }
         }
       }
      if (p != my_rank)
       {
        nnz_to_be_recv[p] = new int[nrow_recv];
       }
      nnz_to_be_sent[p] = new int[nrow_send];
      nrow_send_p[p] = nrow_send;
      nrow_recv_p[p] = nrow_recv;
      nmatrix_send_p[p] = nmatrix_send;

      // if data to send then compute and send
      if (nrow_send > 0)
       {
        unsigned ptr = 0;
        for (unsigned ib = 0; ib < nblock; ib++)
         {
          if (nrow_to_send[ib][p] > 0)
           {
            unsigned l = first_row_to_send[ib][p]-
             Block_distribution_pt[ib]->first_row();
            unsigned u = l + nrow_to_send[ib][p];
            for (unsigned jb = 0; jb < nblock; jb++)
             {
              if (block_matrix_pt(ib,jb) != 0)
               {
                int* row_start = block_matrix_pt(ib,jb)->row_start();
                for (unsigned k = l; k < u; k++)
                 {
                  nnz_to_be_sent[p][ptr] = row_start[k+1] - row_start[k];
                  ptr++;
                 }
               }
             }
           }
         }

        // and send
	if (p != my_rank)
         {
          MPI_Request req;
          MPI_Isend(nnz_to_be_sent[p],nrow_send,MPI_INT,p,0,
                    Problem_pt->communicator_pt()->mpi_comm(),
                    &req);
          send_reqs_nnz.push_back(req);
         }
	else
         {
          nnz_to_be_recv[p] = nnz_to_be_sent[p];
          nnz_to_be_sent[p] = 0;
         }
       }

      // if data to be recv then post recvs
      if (nrow_recv > 0 && p != my_rank)
       {
        MPI_Request req;
        MPI_Irecv(nnz_to_be_recv[p],nrow_recv,MPI_INT,p,0,
                  Problem_pt->communicator_pt()->mpi_comm(),
                  &req);
        recv_reqs_nnz.push_back(req);
       }
     }

    // next send the values and column indices

    // the base displacement
    MPI_Aint base_displacement;
    int* comm_pt = new int;
    MPI_Address(comm_pt,&base_displacement);
    Vector<MPI_Request> send_reqs_coefs;
    for (unsigned p = 0; p < nproc; p++)
     {
      if (p != my_rank)
       {
        if (nrow_send_p[p] > 0)
	 {
          MPI_Datatype types[nmatrix_send_p[p]*2];
          int sz[nmatrix_send_p[p]*2];
          MPI_Aint displacements[nmatrix_send_p[p]*2];
          unsigned ptr = 0;
          unsigned nmsp = nmatrix_send_p[p];
          for (unsigned ib = 0; ib < nblock; ib++)
           {
            if (nrow_to_send[ib][p] > 0)
             {
              unsigned l = first_row_to_send[ib][p]-
               Block_distribution_pt[ib]->first_row();
              unsigned u = l + nrow_to_send[ib][p];
              for (unsigned jb = 0; jb < nblock; jb++)
               {
                if (block_matrix_pt(ib,jb) != 0)
                 {
                  // the underlying storage
                  int* row_start = block_matrix_pt(ib,jb)->row_start();
                  int* column_index = block_matrix_pt(ib,jb)->column_index();
                  double* value = block_matrix_pt(ib,jb)->value();
		       
                  // number of coefs to be sent
                  unsigned ns = row_start[u]-row_start[l];
		       
                  // assemble the data for send
                  // column index
                  MPI_Type_contiguous(ns,MPI_INT,&types[ptr]);
                  MPI_Type_commit(&types[ptr]);
                  MPI_Address(&column_index[row_start[l]],&displacements[ptr]);
                  displacements[ptr] -= base_displacement;
                  sz[ptr] = 1;
		       
                  // values
                  MPI_Type_contiguous(ns,MPI_DOUBLE,&types[nmsp+ptr]); 
                  MPI_Type_commit(&types[nmsp+ptr]);
                  MPI_Address(&value[row_start[l]],&displacements[nmsp+ptr]);
                  displacements[nmsp+ptr] -= base_displacement;
                  sz[nmsp+ptr] = 1;
                  ptr++;
                 }
               }
             }
           }

          // build the final send datatype
          MPI_Datatype final_send_type;
          MPI_Type_struct(nmatrix_send_p[p]*2,sz,displacements,types,
                          &final_send_type);
          MPI_Type_commit(&final_send_type);
          for (unsigned i = 0; i < nmatrix_send_p[p]*2; i++)
           {
            MPI_Type_free(&types[i]);
           }

          // send
          MPI_Request req;
          MPI_Isend(comm_pt,1,final_send_type,p,1,
                    Problem_pt->communicator_pt()->mpi_comm(),&req);
          send_reqs_coefs.push_back(req);
          MPI_Type_free(&final_send_type);
	 }
       }
     }

    // before we can receive the column indices and values we
    // must create the storage and build the mpi data types 
    // hence the nnz receive must be complete

    // wait for recv to complete
    unsigned n_recv_req_nnz = recv_reqs_nnz.size();
    if (n_recv_req_nnz)
     {
    MPI_Waitall(n_recv_req_nnz,&recv_reqs_nnz[0],
                MPI_STATUS_IGNORE);
     }
    
    // first compute the total number of nnzs to recv
    unsigned total_nnz_recv = 0;
    for (unsigned p = 0; p < nproc; p++)
     {
      for (unsigned i = 0; i < nrow_recv_p[p]; i++)
       {
        total_nnz_recv += nnz_to_be_recv[p][i];
       }
     }

    // compute row start
    unsigned nrow_local = Distribution_pt->nrow_local();
    int* row_start = new int[nrow_local+1];
    for (unsigned i = 0; i < nrow_local; i++)
     {
      row_start[i] = 0;
     }

    // initially row start is populated with the number of non zeros in each 
    // row
    unsigned first_row = Distribution_pt->first_row();
    for (unsigned p = 0; p < nproc; p++)
     {
      unsigned block_offset = 0;
      unsigned ptr = 0;
      for (unsigned b = 0; b < nblock; b++)
       {
        if (nrow_to_recv[b][p] > 0)
         {
          unsigned l = first_row_to_recv[b][p]+block_offset-first_row;
          unsigned u = l + nrow_to_recv[b][p];
          for (unsigned c = 0; c < nblock; c++)
           {
            if (block_matrix_pt(b,c)!=0)
             {
              for (unsigned i = l; i < u; i++)
               {
                row_start[i] += nnz_to_be_recv[p][ptr];
                ptr++;
               }
             }
           }
         }
        block_offset += this->block_dimension(b);  
       }
     }




    // next we compute row start
    int g = row_start[0];
    row_start[0] = 0;
    for (unsigned i = 1; i < nrow_local; i++)
     {
      int h = row_start[i];
      row_start[i] = row_start[i-1] + g;
      g = h;
     }
    row_start[nrow_local] = row_start[nrow_local-1]+g;

    // storage for the index vector for the mpi
    // indexed datatype
    Vector<int*> offset_for_recv(nproc);
    for (unsigned p = 0; p < nproc; p++)
     {
      offset_for_recv[p] = new int[nrow_recv_p[p]];
     }

    // next we compute the offsets for the mpi indexed datatype
    // row_start will be altered, but can be easily fixed
    for (unsigned p = 0; p < nproc; p++)
     {
      unsigned block_offset = 0;
      unsigned ptr = 0;
      for (unsigned b = 0; b < nblock; b++)
       {
        if (nrow_to_recv[b][p] > 0)
         {
          unsigned l = first_row_to_recv[b][p]+block_offset-first_row;
          unsigned u = l + nrow_to_recv[b][p];
          for (unsigned c = 0; c < nblock; c++)
           {
            if (block_matrix_pt(b,c) != 0)
             {
              for (unsigned i = l; i < u; i++)
               {
                offset_for_recv[p][ptr] = row_start[i];
                row_start[i] += nnz_to_be_recv[p][ptr];
                ptr++;
               }
             }
           }
         }
        block_offset += this->block_dimension(b);  
       }
     }

    // resize the storage for the column indices and the values
    int* column_index = new int[total_nnz_recv];
    double* values = new double[total_nnz_recv];

    // recv the column indices and values
    Vector<MPI_Request> recv_reqs_coefs;
    Vector<unsigned> recv_reqs_coefs_proc;
    MPI_Aint displacements[2];
    MPI_Address(column_index,&displacements[0]);
    MPI_Address(values,&displacements[1]);
    displacements[1] -= displacements[0];
    displacements[0] = 0;
    int sz[2];
    sz[0] = 1;
    sz[1] = 1;
    for (unsigned p = 0; p < nproc; p++)
     {
      if (p != my_rank)
       {
        if (nrow_recv_p[p] > 0)
         {
          // the datatypes
          // index 1 - column index
          //       2 - values
          MPI_Datatype types[2];
		
          // build the type for the column indices
          MPI_Type_indexed(nrow_recv_p[p],nnz_to_be_recv[p],
                           offset_for_recv[p],MPI_INT,&types[0]);
          MPI_Type_commit(&types[0]);
		
          // build the type for the values
          MPI_Type_indexed(nrow_recv_p[p],nnz_to_be_recv[p],
                           offset_for_recv[p],MPI_DOUBLE,&types[1]);
          MPI_Type_commit(&types[1]);
		
          // build the combined datatype
          MPI_Datatype final_recv_type;
          MPI_Type_struct(2,sz,displacements,types,&final_recv_type);
          MPI_Type_commit(&final_recv_type);
          MPI_Type_free(&types[0]);
          MPI_Type_free(&types[1]);
		
          // and recv
          MPI_Request req;
          MPI_Irecv(column_index,1,final_recv_type,p,1,
                    Problem_pt->communicator_pt()->mpi_comm(),&req);
          MPI_Type_free(&final_recv_type);
          recv_reqs_coefs.push_back(req);
          recv_reqs_coefs_proc.push_back(p);
         }
       }
     }

    // revert row start
    for (unsigned i = nrow_local; i > 0; i--)
     {
      row_start[i] = row_start[i-1];
     }
    row_start[0] = 0;

    // receive from self
    if (nrow_recv_p[my_rank] > 0)
     {
      unsigned row_ptr = 0;
      for (unsigned b = 0; b < nblock; b++)
       {
	if (nrow_to_recv[b][my_rank] > 0)
         {
          unsigned nr = nrow_to_recv[b][my_rank];
          unsigned column_index_offset = 0;
          for (unsigned c = 0; c < nblock; c++)
           {
            if (block_matrix_pt(b,c) != 0)
             {
              double* block_values = block_matrix_pt(b,c)->value();
              int* block_column_indices 
               = block_matrix_pt(b,c)->column_index();
              int* block_row_start 
               = block_matrix_pt(b,c)->row_start(); 
              for (unsigned i = 0; i < nr; i++)
               {
                unsigned n = block_row_start[i+1]-block_row_start[i];
                for (unsigned j = 0; j < n; j++)
                 {
                  values[offset_for_recv[my_rank][row_ptr]+j] =
                   block_values[block_row_start[i]+j];
                  column_index[offset_for_recv[my_rank][row_ptr]+j] =
                   block_column_indices[block_row_start[i]+j]
                   + column_index_offset;
                 }
                row_ptr++;
               }
             }
            column_index_offset += this->block_dimension(c);
           }
         }
       }

      // clean
      delete[] offset_for_recv[my_rank];
      delete[] nnz_to_be_recv[my_rank];
     }

    // update the column indices as the recv is completed
    unsigned n_coef_recv_req = recv_reqs_coefs.size();
    while (n_coef_recv_req > 0)
     {
      // wait for a communication to complete
      int i;
      MPI_Waitany(n_coef_recv_req,&recv_reqs_coefs[0],&i,MPI_STATUS_IGNORE);
      int p = recv_reqs_coefs_proc[i];
      recv_reqs_coefs.erase(recv_reqs_coefs.begin()+i);
      recv_reqs_coefs_proc.erase(recv_reqs_coefs_proc.begin()+i);
      n_coef_recv_req--;

      // update the column indices
      unsigned ptr = 0;
      for (unsigned b = 0; b < nblock; b++)
       {
        unsigned nr = nrow_to_recv[b][p];
        if (nr > 0)
         {
          unsigned column_index_offset = 0;
          for (unsigned c = 0; c < nblock; c++)
           {
            if (block_matrix_pt(b,c) != 0)
             {
              for (unsigned i = 0; i < nr; i++)
               {
                unsigned n = nnz_to_be_recv[p][ptr];
                for (unsigned j = 0; j < n; j++)
                 {
                  column_index[offset_for_recv[p][ptr]+j]
                   += column_index_offset;
                 }
                ptr++;
               }
              column_index_offset += this->block_dimension(c);
             }
           }
         }
       }
	
      // clean
      delete[] offset_for_recv[p];
      delete[] nnz_to_be_recv[p];   
     }

    // build the matrix
    if (preconditioner_matrix_pt == 0)
     {
      preconditioner_matrix_pt = new CRDoubleMatrix(this->Distribution_pt);
     }
    preconditioner_matrix_pt->
     rebuild_matrix_without_copy(Distribution_pt->nrow(),
                                 total_nnz_recv,
                                 values,column_index,
                                 row_start);

    // wait for the nnz requests to complete
    unsigned n_send_req_nnz = send_reqs_nnz.size();
    if (n_send_req_nnz)
     {
      MPI_Waitall(n_send_req_nnz,&send_reqs_nnz[0],MPI_STATUS_IGNORE);
     }
    for (unsigned p = 0; p < nproc; p++)
     {
      delete[] nnz_to_be_sent[p];
     } 

    // wait for the indices/value send requests to complete
    unsigned n_send_reqs_coefs = send_reqs_coefs.size();
    if (n_send_reqs_coefs)
     {
      MPI_Waitall(n_send_reqs_coefs,&send_reqs_coefs[0],MPI_STATUS_IGNORE);
     }
#endif    
   }
 }
}

