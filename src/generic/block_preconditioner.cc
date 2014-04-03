//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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



 //============================================================================
 /// Determine the size of the matrix blocks and setup the
 /// lookup schemes relating the global degrees of freedom with
 /// their "blocks" and their indices (row/column numbers) in those
 /// blocks.
 /// The distributions of the preconditioner and the blocks are
 /// automatically specified (and assumed to be uniform) at this
 /// stage.
 /// This method should be used if any block contains more than one
 /// type of DOF. The argument vector dof_to_block_map should be of length
 /// ndof. Each element should contain an integer indicating the block number
 /// corresponding to that type of DOF.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 block_setup(Vector<unsigned>& dof_to_block_map)
 {

  if(is_subsidiary_block_preconditioner())
   {
    // Compute number of rows in this (sub) preconditioner using data from
    // the master.
    Nrow = 0;
    for (unsigned b = 0; b < Internal_ndof_types; b++)
     {
      Nrow += this->dof_block_dimension(b);
     }

#ifdef PARANOID
    if (Nrow==0)
     {
      std::ostringstream error_message;
      error_message
       << "Nrow=0 in subsidiary preconditioner. This seems fishy and\n"
       << "suggests that block_setup() was not called for the \n"
       << "master block preconditioner yet.";
      throw OomphLibWarning(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
     }
#endif
   }

  if(Preconditioner_blocks_have_been_precomputed)
   {
    // Check that the dof_to_block map "makes sense" for the 
    // Doftype_to_doftype_map.
    // The Doftype_to_doftype_map describes which doftypes should be considered
    // as a single doftype in THIS preconditioner.
    // For example, if this preconditioner is LSC in 3D, it expects 4 doftypes;
    // 3 velocity, [u, v, w] and 1 pressure [p], giving us the doftype ordering
    // [u v w p].
    // The LSC preconditioner groups the velocity and pressure doftypes 
    // separately, thus the dof_to_block_map will be:
    // [0 0 0 1]
    // 
    // The Doftype_to_doftype_map MUST have length 4, to identify which of the
    // OTHER doftypes should be grouped together to be one of THIS doftypes.
    //
    // For example, if the preconditioner above this one has the doftype 
    // ordering:
    // 0  1  2  3  4  5  6  7  8  9
    // ub vb wb up vp wp ut vt wt p
    // Then we want to tell THIS preconditioner which doftypes belongs to 
    // u, v, w and p, by giving the Doftype_to_doftype_map:
    // [0 3 6]
    // [1 4 7]
    // [2 5 8]
    // [9]
    //
    // So, it is important that the length of dof_to_block_map is the same as
    // the length of Doftype_to_doftype_map.
    unsigned dof_to_block_map_size = dof_to_block_map.size();
    if(dof_to_block_map_size != Doftype_to_doftype_map.size())
     {
      std::ostringstream err_msg;
      err_msg
       << "The size of dof_to_block_map and Doftype_to_doftype_map is not "
       << "the same.\n"
       << "One of the two list is incorrect, please look at the comments\n"
       << "in the source code for more details.";
      throw OomphLibWarning(err_msg.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
     }
    
    // Create the Block_to_block_map from 
    // the dof_to_block_map and Doftype_to_doftype_map.

    // find the maximum block number
    unsigned max_block_number = 0;
    for (unsigned i = 0; i < dof_to_block_map_size; i++)
     {
      if (dof_to_block_map[i] > max_block_number)
       {
        max_block_number = dof_to_block_map[i];
       }
     }

    // Now we do the following:
    // Lets say the Doftype_to_doftype_map is:
    // [0 3 6]
    // [1 4 7]
    // [2 5 8]
    // [9]
    //
    // and the dof_to_block_map is [0 0 0 1].
    //
    // Then we need to form the Block_to_block_map:
    // [0 3 6 1 4 7 2 5 8]
    // [9]
    Block_to_block_map.clear();
    for (unsigned blocktype_i = 0; blocktype_i <= max_block_number; 
         blocktype_i++)
     {
      // Temp vector to store the doftypes.
      Vector<unsigned> temp_doftype_vec;

      // Loop through the entires in dof_to_block_map
      for (unsigned i = 0; i < dof_to_block_map_size; i++) 
       {
        // If the entry in dof_to_block_map matches the current blocktype_i
        // push all entries of Doftype_to_doftype_map[i] into temp_doftype_vec.
        if(dof_to_block_map[i] == blocktype_i)
         {
          unsigned doftype_to_doftype_map_i_size 
            = Doftype_to_doftype_map[i].size();
          for (unsigned j = 0; j < doftype_to_doftype_map_i_size; j++) 
           {
            temp_doftype_vec.push_back(Doftype_to_doftype_map[i][j]);
           }
         }
       }
      Block_to_block_map.push_back(temp_doftype_vec);
     }

    // Now set the dof_to_block_map to the identify.
    dof_to_block_map.resize(Internal_ndof_types,0);
    for (unsigned i = 0; i < Internal_ndof_types; i++) 
     {
      dof_to_block_map[i] = i;
     }
   }

#ifdef PARANOID

  // Check that the meshes are ok. This only needs to be done in the master
  // because subsidiary preconditioners don't do anything with the meshes
  // here.
  if(is_master_block_preconditioner())
   {
    // This is declared as local_nmesh because there are other variables
    // called nmesh floating about... but this will not exist if PARANOID is
    // switched on.
    unsigned local_nmesh = nmesh();
  
    // Check that some mesh pointers have been assigned.
    if(local_nmesh == 0)
     {
      std::ostringstream error_msg;
      error_msg << "Cannot setup blocks because no meshes have been set.";
      throw OomphLibError(error_msg.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
 
    // Each mesh must contain elements with the same number of dof.
    // A stricter check is to ensure that the mesh contains only one type of
    // elements. This is relaxed in same cases.
    for (unsigned mesh_i = 0; mesh_i < local_nmesh; mesh_i++) 
     {
      // The number of elements in the current mesh.
      unsigned n_element = mesh_pt(mesh_i)->nelement();
      
      // When the bulk mesh is distributed, there may not be any elements
      // in the surface mesh(es).
      if(n_element > 0)
      {
        // The string of the first element in the current mesh.
        std::string first_element_string
         =typeid(*(mesh_pt(mesh_i)->element_pt(0))).name();
  
      // If there are multiple element types in the current mesh, 
      // we can at least make sure that they contain the same types of DOFs.
      if(bool(Allow_multiple_element_type_in_mesh[mesh_i]))
       {
        // The ndof types of the first element.
        unsigned first_element_ndof_type = 
         mesh_pt(mesh_i)->element_pt(0)->ndof_types();

        // Loop through the meshes and compare the number of types of DOFs.
        for (unsigned el_i = 1; el_i < n_element; el_i++) 
         {
          // The ndof type of the current element.
          unsigned current_element_ndof_type =
           mesh_pt(mesh_i)->element_pt(el_i)->ndof_types();

          // The string of the current element.
          std::string current_element_string
           =typeid(*(mesh_pt(mesh_i)->element_pt(el_i))).name();   

          // Compare against the first element.
          if(current_element_ndof_type != first_element_ndof_type)
           {
            std::ostringstream error_message;
            error_message 
             << "Elements in the same mesh MUST have the same number of types "
             << "of DOFs.\n"
             << "The element in mesh " << mesh_i << ", at position "
             << el_i << " is: \n"
             << current_element_string <<", \n"
             << "with ndof types: " << current_element_ndof_type << ".\n"
             << "The first element in the same mesh is: \n"
             << first_element_string << ", \n"
             << "with ndof types: " << first_element_ndof_type << ".\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
           }
         }
       }
      else
       // There should be only one type of elements in the current mesh. Check
       // that this is the case!
       {
        // Loop through the elements in the current mesh.
        for (unsigned el_i = 1; el_i < n_element; el_i++) 
         {
          // The string of the current element.
          std::string current_element_string
           =typeid(*(mesh_pt(mesh_i)->element_pt(el_i))).name();      
  
          // Compare against the first element.
          if(current_element_string.compare(first_element_string) != 0)
           {
            std::ostringstream error_message;
            error_message 
             << "By default, a mesh containing block preconditionable "
             << "elements must contain only one type of element.\n"
             << "The element in mesh " << mesh_i << ", at position "
             << el_i << " is: \n" << current_element_string << "\n"
             << "The first element in the same mesh is: \n"
             << first_element_string << "\n"
             << "If this is correct, consider calling the set_mesh(...) with\n"
             << "the optional argument set true to allow multiple element\n"
             << "types in the same mesh.\n"
             << "Note: A minimal requirement is that the elements in the same\n"
             << "mesh MUST have the same number types DOFs.";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
           }
         }
       }
       }// Here
     }
   }

#endif


  // clear the memory
  this->clear_block_preconditioner_base();

  // get my_rank and nproc
#ifdef OOMPH_HAS_MPI
  unsigned my_rank = comm_pt()->my_rank();
  unsigned nproc = comm_pt()->nproc();
#endif


  /////////////////////////////////////////////////////////////////////////////
  // start of master block preconditioner only operations
  /////////////////////////////////////////////////////////////////////////////
#ifdef OOMPH_HAS_MPI
  unsigned* nreq_sparse = new unsigned[nproc];
  unsigned* nreq_sparse_for_proc = new unsigned[nproc];
  unsigned** index_in_dof_block_sparse_send = new unsigned*[nproc];
  unsigned** dof_number_sparse_send = new unsigned*[nproc];
  for (unsigned p = 0; p < nproc; p++)
   {
    index_in_dof_block_sparse_send[p] = 0;
    dof_number_sparse_send[p] = 0;
   }
  Vector<MPI_Request> send_requests_sparse;
  Vector<MPI_Request> recv_requests_sparse;
#endif

  // if this preconditioner is the master preconditioner then we need
  // to assemble the vectors : Dof_number
  //                           Index_in_dof_block
  if (is_master_block_preconditioner())
   {

    // Get the number of dof types in each mesh.
    Ndof_types_in_mesh.assign(nmesh(),0);
    for(unsigned i=0; i<nmesh(); i++)
     {
      Ndof_types_in_mesh[i] = mesh_pt(i)->ndof_types();
     }

    // Setup the distribution of this preconditioner, assumed to be the same
    // as the matrix if the matrix is distributable.
    if (dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt()))
     {
      this->build_distribution
       (dynamic_cast<DistributableLinearAlgebraObject*>
        (matrix_pt())->distribution_pt());
     }
    else
     {
      LinearAlgebraDistribution dist(comm_pt(),
                                     matrix_pt()->nrow(),false);
      this->build_distribution(dist);
     }
    Nrow = matrix_pt()->nrow();

    // boolean to indicate whether the matrix is actually distributed
    // ie distributed and on more than one processor
    bool matrix_distributed =
     (this->distribution_pt()->distributed() &&
      this->distribution_pt()->communicator_pt()->nproc() > 1);


    // matrix must be CR
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());
    if (cr_matrix_pt == 0) {
     std::ostringstream error_message;
     error_message << "Block setup for distributed matrices only works "
                   << "for CRDoubleMatrices";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }




    // my distribution
    unsigned first_row = this->distribution_pt()->first_row();
    unsigned nrow_local = this->distribution_pt()->nrow_local();
    unsigned last_row = first_row+nrow_local-1;

#ifdef OOMPH_HAS_MPI

    // storage for the rows required by each processor in the dense
    // block lookup storage scheme
    // dense_required_rows(p,0) is the minimum global index required by proc p
    //                 ...(p,1) is the maximum global index required by proc p
    DenseMatrix<unsigned> dense_required_rows(nproc,2);
    for (unsigned p = 0; p < nproc; p++)
     {
      dense_required_rows(p,0) = this->distribution_pt()->first_row(p);
      dense_required_rows(p,1) = this->distribution_pt()->first_row(p)
       +this->distribution_pt()->nrow_local(p) - 1;
     }

    // determine the global rows That are not in the range first_row to
    // first_row+nrow_local for which we should store the
    // Dof_index and Index_in_dof_block for
    // then send the lists to other processors
    Vector<unsigned> sparse_global_rows_for_block_lookup;
    if (matrix_distributed)
     {
      unsigned nnz = cr_matrix_pt->nnz();
      int* column_index = cr_matrix_pt->column_index();
      for (unsigned i = 0; i < nnz; i++)
       {
        unsigned ci = column_index[i];
        if (ci<first_row || ci>last_row)
         {
          if (find(sparse_global_rows_for_block_lookup.begin(),
                   sparse_global_rows_for_block_lookup.end(),
                   ci) ==
              sparse_global_rows_for_block_lookup.end())
           {
            sparse_global_rows_for_block_lookup.push_back(ci);
           }
         }
       }
     }
    sort(sparse_global_rows_for_block_lookup.begin(),
         sparse_global_rows_for_block_lookup.end());
    int nsparse = sparse_global_rows_for_block_lookup.size();
    Global_index_sparse.resize(nsparse);
    Index_in_dof_block_sparse.resize(nsparse);
    Dof_number_sparse.resize(nsparse);
    // RAY - Why is this looping through a vector?...
    for (int i = 0; i < nsparse; i++)
     {
      Global_index_sparse[i]=sparse_global_rows_for_block_lookup[i];
     }
    sparse_global_rows_for_block_lookup.clear();


    Vector<MPI_Request> recv_requests_sparse_nreq;
    if (matrix_distributed)
     {
      MPI_Aint base_displacement_sparse;
      MPI_Address(nreq_sparse,&base_displacement_sparse);

      int zero = 0;
      for (unsigned p = 0; p < nproc; p++)
       {
        nreq_sparse[p] = 0;
        nreq_sparse_for_proc[p] = 0;

        // determine the global eqn numbers required by this processor
        // that can be classified by processor p
        int begin = 0;
        for (int i = 0; i < nsparse; ++i)
         {
          if (Global_index_sparse[i]<dense_required_rows(p,0))
           {
            ++begin;
           }
          else
           {
            if (Global_index_sparse[i]<=dense_required_rows(p,1))
             {
              ++nreq_sparse[p];
             }
            else
             {
              break;
             }
           }
         }

        // if this processor has rows to be classified by proc p
        if (nreq_sparse[p]>0)
         {

          // send the number of global eqn numbers
          MPI_Request req1;
          MPI_Isend(&nreq_sparse[p],1,MPI_UNSIGNED,p,31,
                    comm_pt()->mpi_comm(),&req1);
          send_requests_sparse.push_back(req1);

          // send the global eqn numbers
          MPI_Request req2;
          MPI_Isend(&Global_index_sparse[begin],
                    nreq_sparse[p],MPI_UNSIGNED,p,32,
                    comm_pt()->mpi_comm(),&req2);
          send_requests_sparse.push_back(req2);

          // post the recvs for the data that will be returned

          // the datatypes, displacements, lengths for the two datatypes
          MPI_Datatype types[2];
          MPI_Aint displacements[2];
          int lengths[2];

          // index in dof block
          MPI_Type_contiguous(nreq_sparse[p],MPI_UNSIGNED,&types[0]);
          MPI_Type_commit(&types[0]);
          MPI_Address(&Index_in_dof_block_sparse[begin],&displacements[0]);
          displacements[0] -= base_displacement_sparse;
          lengths[0] = 1;

          // dof number
          MPI_Type_contiguous(nreq_sparse[p],MPI_UNSIGNED,&types[1]);
          MPI_Type_commit(&types[1]);
          MPI_Address(&Dof_number_sparse[begin],&displacements[1]);
          displacements[1] -= base_displacement_sparse;
          lengths[1] = 1;

          // build the final type
          MPI_Datatype recv_type;
          MPI_Type_struct(2,lengths,displacements,types,&recv_type);
          MPI_Type_commit(&recv_type);
          MPI_Type_free(&types[0]);
          MPI_Type_free(&types[1]);

          // and recv
          MPI_Request req;
          MPI_Irecv(nreq_sparse,1,recv_type,p,33,
                    comm_pt()->mpi_comm(),&req);
          recv_requests_sparse.push_back(req);
          MPI_Type_free(&recv_type);
         }

        // if no communication required, confirm this
        if (nreq_sparse[p]==0)
         {
          MPI_Request req1;
          MPI_Isend(&zero,1,MPI_UNSIGNED,p,31,
                    comm_pt()->mpi_comm(),&req1);
          send_requests_sparse.push_back(req1);
         }

        //
        MPI_Request req;
        MPI_Irecv(&nreq_sparse_for_proc[p],1,MPI_UNSIGNED,p,31,
                  comm_pt()->mpi_comm(),&req);
        recv_requests_sparse_nreq.push_back(req);
       }
     }
#endif

    // resize the storage
    Dof_number_dense.resize(nrow_local);
    Index_in_dof_block_dense.resize(nrow_local);

    // zero the number of dof types
    Internal_ndof_types = 0;

#ifdef PARANOID
    // Vector to keep track of previously assigned block numbers
    // to check consistency between multiple assignments.
    Vector<int> previously_assigned_block_number(nrow_local,
                                                 Data::Is_unclassified);
#endif

    // determine whether the problem is distribution
    bool problem_distributed = false;

    // the problem method distributed() is only accessible with MPI
#ifdef OOMPH_HAS_MPI
    problem_distributed = any_mesh_distributed();
#endif

    // if the problem is not distributed
    if (!problem_distributed)
     {

      // Offset for the block type in the overall system.
      // Different meshes contain different block-preconditionable
      // elements -- their blocks are added one after the other...
      unsigned dof_offset=0;

      // Loop over all meshes
      for (unsigned m=0;m<nmesh();m++)
       {
        // Number of elements in this mesh
        unsigned n_element = mesh_pt(m)->nelement();

        // Find the number of block types that the elements in this mesh
        // are in charge of
        unsigned ndof_in_element = ndof_types_in_mesh(m);
        Internal_ndof_types += ndof_in_element;

        for (unsigned e=0;e<n_element;e++)
         {
          // List containing pairs of global equation number and
          // dof number for each global dof in an element
          std::list<std::pair<unsigned long,unsigned> > dof_lookup_list;

          // Get list of blocks associated with the element's global unknowns
          mesh_pt(m)->element_pt(e)->
           get_dof_numbers_for_unknowns(dof_lookup_list);

          // Loop over all entries in the list
          // and store the block number
          typedef std::list<std::pair<unsigned long,unsigned> >::iterator IT;
          for (IT it=dof_lookup_list.begin();
               it!=dof_lookup_list.end();it++)
           {

            unsigned long global_dof = it->first;
            if (global_dof >= unsigned(first_row) &&
                global_dof <= unsigned(last_row))
             {
              unsigned dof_number = (it->second)+dof_offset;
              Dof_number_dense[global_dof-first_row]
               = dof_number;

#ifdef PARANOID
              // Check consistency of block numbers if assigned multiple times
              if (previously_assigned_block_number[global_dof-
                                                   first_row]<0)
               {
                previously_assigned_block_number[global_dof-first_row]
                 =dof_number;
               }
#endif
             }
           }
         }

        // About to do the next mesh which contains block preconditionable
        // elements of a different type; all the dofs that these elements are
        // "in charge of" differ from the ones considered so far.
        // Bump up the block counter to make sure we're not overwriting
        // anything here
        dof_offset+=ndof_in_element;
       }

#ifdef PARANOID
      // check that every global equation number has been allocated a dof type
      for (unsigned i = 0; i < nrow_local; i++)
       {
        if (previously_assigned_block_number[i] < 0)
         {
          std::ostringstream error_message;
          error_message << "Not all degrees of freedom have had DOF type "
                        << "numbers allocated. Dof number " << i
                        << " is unallocated.";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
         }
       }
#endif
     }
    // else the problem is distributed
    else
     {
#ifdef OOMPH_HAS_MPI


      // Offset for the block type in the overall system.
      // Different meshes contain different block-preconditionable
      // elements -- their blocks are added one after the other...
      unsigned dof_offset=0;

      // the set of global degrees of freedom and their corresponding dof
      // number on this processor
      std::map<unsigned long,unsigned > my_dof_map;

      // Loop over all meshes
      for (unsigned m=0;m<nmesh();m++)
       {
        // Number of elements in this mesh
        unsigned n_element = this->mesh_pt(m)->nelement();

        // Find the number of block types that the elements in this mesh
        // are in charge of
        unsigned ndof_in_element = ndof_types_in_mesh(m);
        Internal_ndof_types += ndof_in_element;

        // Loop over all elements
        for (unsigned e=0;e<n_element;e++)
         {

          // if the element is not a halo element
          if (!this->mesh_pt(m)->element_pt(e)->is_halo())
           {
            // List containing pairs of global equation number and
            // dof number for each global dof in an element
            std::list<std::pair<unsigned long,unsigned> > dof_lookup_list;

            // Get list of blocks associated with the element's global
            // unknowns
            this->mesh_pt(m)->element_pt(e)->
             get_dof_numbers_for_unknowns(dof_lookup_list);

            // update the block numbers and put it in the map.
            typedef
             std::list<std::pair<unsigned long,unsigned> >::iterator IT;
            for (IT it=dof_lookup_list.begin();
                 it!=dof_lookup_list.end();it++)
             {
              it->second = (it->second)+dof_offset;
              my_dof_map[it->first] = it->second;
             }
           }
         }

        // About to do the next mesh which contains block preconditionable
        // elements of a different type; all the dofs that these elements are
        // "in charge of" differ from the ones considered so far.
        // Bump up the block counter to make sure we're not overwriting
        // anything here
        dof_offset+=ndof_in_element;
       }

      // next copy the map of my dofs to two vectors to send
      unsigned my_ndof = my_dof_map.size();
      unsigned long* my_global_dofs = new unsigned long[my_ndof];
      unsigned* my_dof_numbers = new unsigned[my_ndof];
      typedef
       std::map<unsigned long,unsigned >::iterator IT;
      unsigned pt = 0;
      for (IT it = my_dof_map.begin(); it != my_dof_map.end(); it++)
       {
        my_global_dofs[pt] = it->first;
        my_dof_numbers[pt] = it->second;
        pt++;
       }

      // and then clear the map
      my_dof_map.clear();

      // count up how many DOFs need to be sent to each processor
      int* first_dof_to_send = new int[nproc];
      int* ndof_to_send = new int[nproc];
      unsigned ptr = 0;
      for (unsigned p = 0; p < nproc; p++)
       {
        first_dof_to_send[p] = 0;
        ndof_to_send[p] = 0;
        while (ptr < my_ndof && my_global_dofs[ptr] < dense_required_rows(p,0))
         {
          ptr++;
         }
        first_dof_to_send[p] = ptr;
        while (ptr < my_ndof && my_global_dofs[ptr] <= dense_required_rows(p,1))
         {
          ndof_to_send[p]++;
          ptr++;
         }
       }

      // next communicate to each processor how many dofs it expects to recv
      int* ndof_to_recv = new int[nproc];
      MPI_Alltoall(ndof_to_send,1,MPI_INT,ndof_to_recv,1,MPI_INT,
                   comm_pt()->mpi_comm());

      // the base displacements for the sends
      MPI_Aint base_displacement;
      MPI_Address(my_global_dofs,&base_displacement);

#ifdef PARANOID
      // storage for paranoid check to ensure that every row as been
      // imported
      std::vector<bool> dof_recv(nrow_local,false);
#endif

      // next send and recv
      Vector<MPI_Request> send_requests;
      Vector<MPI_Request> recv_requests;
      Vector<unsigned long*> global_dofs_recv(nproc,0);
      Vector<unsigned*> dof_numbers_recv(nproc,0);
      Vector<unsigned> proc;
      for (unsigned p = 0; p < nproc; p++)
       {
        if (p != my_rank)
         {

          // send
          if (ndof_to_send[p] > 0)
           {
            // the datatypes, displacements, lengths for the two datatypes
            MPI_Datatype types[2];
            MPI_Aint displacements[2];
            int lengths[2];

            // my global dofs
            MPI_Type_contiguous(ndof_to_send[p],MPI_UNSIGNED_LONG,&types[0]);
            MPI_Type_commit(&types[0]);
            MPI_Address(my_global_dofs + first_dof_to_send[p],
                        &displacements[0]);
            displacements[0] -= base_displacement;
            lengths[0] = 1;

            // my dof numbers
            MPI_Type_contiguous(ndof_to_send[p],MPI_UNSIGNED,&types[1]);
            MPI_Type_commit(&types[1]);
            MPI_Address(my_dof_numbers + first_dof_to_send[p],
                        &displacements[1]);
            displacements[1] -= base_displacement;
            lengths[1] = 1;

            // build the final type
            MPI_Datatype send_type;
            MPI_Type_struct(2,lengths,displacements,types,&send_type);
            MPI_Type_commit(&send_type);
            MPI_Type_free(&types[0]);
            MPI_Type_free(&types[1]);

            // and send
            MPI_Request req;
            MPI_Isend(my_global_dofs,1,send_type,p,2,
                      comm_pt()->mpi_comm(),&req);
            send_requests.push_back(req);
            MPI_Type_free(&send_type);
           }

          // and recv
          if (ndof_to_recv[p] > 0)
           {
            // resize the storage
            global_dofs_recv[p] = new unsigned long[ndof_to_recv[p]];
            dof_numbers_recv[p] = new unsigned[ndof_to_recv[p]];
            proc.push_back(p);

            // the datatypes, displacements, lengths for the two datatypes
            MPI_Datatype types[2];
            MPI_Aint displacements[2];
            int lengths[2];

            // my global dofs
            MPI_Type_contiguous(ndof_to_recv[p],MPI_UNSIGNED_LONG,&types[0]);
            MPI_Type_commit(&types[0]);
            MPI_Address(global_dofs_recv[p],&displacements[0]);
            displacements[0] -= base_displacement;
            lengths[0] = 1;

            // my dof numbers
            MPI_Type_contiguous(ndof_to_recv[p],MPI_UNSIGNED,&types[1]);
            MPI_Type_commit(&types[1]);
            MPI_Address(dof_numbers_recv[p],&displacements[1]);
            displacements[1] -= base_displacement;
            lengths[1] = 1;

            // build the final type
            MPI_Datatype recv_type;
            MPI_Type_struct(2,lengths,displacements,types,&recv_type);
            MPI_Type_commit(&recv_type);
            MPI_Type_free(&types[0]);
            MPI_Type_free(&types[1]);

            // and recv
            MPI_Request req;
            MPI_Irecv(my_global_dofs,1,recv_type,p,2,
                      comm_pt()->mpi_comm(),&req);
            recv_requests.push_back(req);
            MPI_Type_free(&recv_type);
           }

         }
        // send to self
        else
         {
          unsigned u = first_dof_to_send[p] + ndof_to_recv[p];
          for (unsigned i = first_dof_to_send[p]; i < u; i++)
           {
#ifdef PARANOID
            // RAYRAY Legacy code, this checks that the DOFs have not been
            // re-classified. This has been changed so that we can re-classify
            // the DOFs. Example: Bulk element classifies it's own DOF types,
            // then a FaceElement will re-classify the bulk element's DOFs.
            // See ImposeParallelOutflowElement.
            // I have not deleted the below so we can easily "revert" back if
            // required.

            // paranoid check
            // check that if row has been imported the block number is the
            // same
            if (dof_recv[my_global_dofs[i]-first_row])
             {
              if (Dof_number_dense[my_global_dofs[i]-first_row]
                  != my_dof_numbers[i])
               {
                std::ostringstream error_message;
                error_message
                 << "Inconsistency in assigment of block numbers\n"
                 << "Global dof " <<  my_global_dofs[i]
                 << "was previously assigned to block "
                 <<  Dof_number_dense[my_global_dofs[i]-first_row]
                 << "\nNow it's been reassigned to block "
                 << my_dof_numbers[i] << ".\n"
                 << "This is most likely because one of your\n"
                 << "elements has classified its degrees of freedom\n"
                 << "wrongly. You should remember that \n"
                 << "GeneralisedElement::get_block_numbers_for_unknowns()\n"
                 << "should only classify the element's OWN dofs and not \n"
                 << "deal with dofs that were created by resizing nodes, \n"
                 << "say. Check that loops over nodal values in that \n"
                 << "function do not use Node::nvalue() to determine the\n"
                 << "dofs to be classified!\n";
                throw OomphLibWarning(error_message.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
               }
             }
            // indicate that this dof has ben recv
            dof_recv[my_global_dofs[i]-first_row] = true;
#endif
            Dof_number_dense[my_global_dofs[i]-first_row] =
             my_dof_numbers[i];
           }
         }
       }

      // recv and store the data
      unsigned c_recv = recv_requests.size();
      while (c_recv > 0)
       {

        // wait for any communication to finish
        int req_number;
        MPI_Waitany(c_recv,&recv_requests[0],&req_number,MPI_STATUS_IGNORE);
        recv_requests.erase(recv_requests.begin()+req_number);
        c_recv--;

        // determine the source processor
        unsigned p = proc[req_number];
        proc.erase(proc.begin()+req_number);

        // import the data
        for (int i  = 0; i < ndof_to_recv[p]; i++)
         {
#ifdef PARANOID
          // RAYRAY Legacy code, this checks that the DOFs have not been
          // re-classified. This has been changed so that we can re-classify
          // the DOFs. Example: Bulk element classifies it's own DOF types,
          // then a FaceElement will re-classify the bulk element's DOFs.
          // See ImposeParallelOutflowElement.
          // I have not deleted the below so we can easily "revert" back if
          // required.

          // paranoid check
          // check that if row has been imported the block number is the same
          if (dof_recv[global_dofs_recv[p][i]-first_row])
           {
            if (Dof_number_dense[global_dofs_recv[p][i]-first_row]
                != dof_numbers_recv[p][i])
             {
              std::ostringstream error_message;
              error_message
               << "Inconsistency in assignment of block numbers\n"
               << "Global dof "
               <<  global_dofs_recv[p][i]
               << " was previously assigned to block "
               <<  Dof_number_dense[global_dofs_recv[p][i]
                                    -first_row]
               << "\nNow it's been reassigned to block "
               << dof_numbers_recv[p][i] << ".\n"
               << "This is most likely because one of your\n"
               << "elements has classified its degrees of freedom\n"
               << "wrongly. You should remember that \n"
               << "GeneralisedElement::get_block_numbers_for_unknowns()\n"
               << "should only classify the element's OWN dofs and not \n"
               << "deal with dofs that were created by resizing nodes, \n"
               << "say. Check that loops over nodal values in that \n"
               << "function do not use Node::nvalue() to determine the\n"
               << "dofs to be classified!\n";
              throw OomphLibWarning(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
             }
           }
          // indicate that this dof has ben recv
          dof_recv[global_dofs_recv[p][i]-first_row] = true;
#endif
          Dof_number_dense[global_dofs_recv[p][i]-first_row]
           = dof_numbers_recv[p][i];
         }

        // delete the data
        delete[] global_dofs_recv[p];
        delete[] dof_numbers_recv[p];
       }

      // finally wait for the send requests to complete as we are leaving
      // an MPI block of code
      unsigned csr = send_requests.size();
      if (csr)
       {
        MPI_Waitall(csr,&send_requests[0],MPI_STATUS_IGNORE);
       }

      // clean up
      delete[] ndof_to_send;
      delete[] first_dof_to_send;
      delete[] ndof_to_recv;
      delete[] my_global_dofs;
      delete[] my_dof_numbers;
#ifdef PARANOID
      unsigned all_recv = true;
      for (unsigned i = 0; i < nrow_local; i++)
       {
        if (!dof_recv[i])
         {
          all_recv = false;
         }
       }
      if (!all_recv)
       {
        std::ostringstream error_message;
        error_message << "Not all the DOF numbers required were received";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
       }
#endif
#else
      std::ostringstream error_message;
      error_message
       << "The problem appears to be distributed, MPI is required";
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
#endif
     }

#ifdef OOMPH_HAS_MPI
    Vector<unsigned*> sparse_rows_for_proc(nproc,0);
    Vector<MPI_Request> sparse_rows_for_proc_requests;
    if (matrix_distributed)
     {
      // wait for number of sparse rows each processor requires
      // post recvs for that data
      if (recv_requests_sparse_nreq.size()>0)
       {
        MPI_Waitall(recv_requests_sparse_nreq.size(),
                    &recv_requests_sparse_nreq[0],
                    MPI_STATUS_IGNORE);
       }
      for (unsigned p = 0; p < nproc; ++p)
       {
        if (nreq_sparse_for_proc[p] > 0)
         {
          MPI_Request req;
          sparse_rows_for_proc[p] = new unsigned[nreq_sparse_for_proc[p]];
          MPI_Irecv(sparse_rows_for_proc[p],nreq_sparse_for_proc[p],
                    MPI_UNSIGNED,p,32,
                    comm_pt()->mpi_comm(),&req);
          sparse_rows_for_proc_requests.push_back(req);
         }
       }
     }
#endif


    // for every global degree of freedom required by this processor we now
    // have the corresponding dof number

    // clear the Ndof_in_dof_block storage
    Dof_dimension.assign(Internal_ndof_types,0);

    // first consider a non distributed matrix
    if (!matrix_distributed)
     {

      // set the Index_in_dof_block
      unsigned nrow  = this->distribution_pt()->nrow();
      Index_in_dof_block_dense.resize(nrow);
      Index_in_dof_block_dense.initialise(0);
      for (unsigned i = 0; i < nrow; i++)
       {
        Index_in_dof_block_dense[i] = Dof_dimension[Dof_number_dense[i]];
        Dof_dimension[Dof_number_dense[i]]++;
       }
     }

    // next a distributed matrix
    else
     {
#ifdef OOMPH_HAS_MPI

      // first compute how many instances of each dof are on this
      // processor
      unsigned* my_nrows_in_dof_block = new unsigned[Internal_ndof_types];
      for (unsigned i = 0; i < Internal_ndof_types; i++)
       {
        my_nrows_in_dof_block[i] = 0;
       }
      for (unsigned i = 0; i < nrow_local; i++)
       {
        my_nrows_in_dof_block[Dof_number_dense[i]]++;
       }

      // next share the data
      unsigned* nrow_in_dof_block_recv = new unsigned[Internal_ndof_types*nproc];
      MPI_Allgather(my_nrows_in_dof_block,Internal_ndof_types,MPI_UNSIGNED,
                    nrow_in_dof_block_recv,Internal_ndof_types,MPI_UNSIGNED,
                    comm_pt()->mpi_comm());
      delete[] my_nrows_in_dof_block;

      // compute my first dof index and Nrows_in_dof_block
      Vector<unsigned> my_first_dof_index(Internal_ndof_types,0);
      for (unsigned i = 0; i < Internal_ndof_types; i++)
       {
        for (unsigned p = 0; p < my_rank; p++)
         {
          my_first_dof_index[i] += nrow_in_dof_block_recv[p*Internal_ndof_types + i];
         }
        Dof_dimension[i] = my_first_dof_index[i];
        for (unsigned p = my_rank; p < nproc; p++)
         {
          Dof_dimension[i] += nrow_in_dof_block_recv[p*Internal_ndof_types + i];
         }
       }
      delete[] nrow_in_dof_block_recv;

      // next compute Index in dof block
      Index_in_dof_block_dense.resize(nrow_local);
      Index_in_dof_block_dense.initialise(0);
      Vector<unsigned> dof_counter(Internal_ndof_types,0);
      for (unsigned i = 0; i < nrow_local; i++)
       {
        Index_in_dof_block_dense[i] =
         my_first_dof_index[Dof_number_dense[i]] +
         dof_counter[Dof_number_dense[i]];
        dof_counter[Dof_number_dense[i]]++;
       }

      // the base displacements for the sends
      if (sparse_rows_for_proc_requests.size()>0)
       {
        MPI_Waitall(sparse_rows_for_proc_requests.size(),
                    &sparse_rows_for_proc_requests[0],
                    MPI_STATUS_IGNORE);
       }
      MPI_Aint base_displacement;
      MPI_Address(dof_number_sparse_send,&base_displacement);
      unsigned first_row = this->distribution_pt()->first_row();
      for (unsigned p = 0; p < nproc; ++p)
       {
        if (nreq_sparse_for_proc[p]>0)
         {
          // construct the data
          index_in_dof_block_sparse_send[p] =
           new unsigned[nreq_sparse_for_proc[p]];
          dof_number_sparse_send[p] =
           new unsigned[nreq_sparse_for_proc[p]];
          for (unsigned i = 0; i < nreq_sparse_for_proc[p]; ++i)
           {
            unsigned r = sparse_rows_for_proc[p][i];
            r -= first_row;
            index_in_dof_block_sparse_send[p][i]
             = Index_in_dof_block_dense[r];
            dof_number_sparse_send[p][i]
             = Dof_number_dense[r];
           }
          delete[] sparse_rows_for_proc[p];

          // send the data
          // the datatypes, displacements, lengths for the two datatypes
          MPI_Datatype types[2];
          MPI_Aint displacements[2];
          int lengths[2];

          // index in dof block
          MPI_Type_contiguous(nreq_sparse_for_proc[p],MPI_UNSIGNED,&types[0]);
          MPI_Type_commit(&types[0]);
          MPI_Address(index_in_dof_block_sparse_send[p],&displacements[0]);
          displacements[0] -= base_displacement;
          lengths[0] = 1;

          // dof number
          MPI_Type_contiguous(nreq_sparse_for_proc[p],MPI_UNSIGNED,&types[1]);
          MPI_Type_commit(&types[1]);
          MPI_Address(dof_number_sparse_send[p],&displacements[1]);
          displacements[1] -= base_displacement;
          lengths[1] = 1;

          // build the final type
          MPI_Datatype send_type;
          MPI_Type_struct(2,lengths,displacements,types,&send_type);
          MPI_Type_commit(&send_type);
          MPI_Type_free(&types[0]);
          MPI_Type_free(&types[1]);

          // and recv
          MPI_Request req;
          MPI_Isend(dof_number_sparse_send,1,send_type,p,33,
                    comm_pt()->mpi_comm(),&req);
          send_requests_sparse.push_back(req);
          MPI_Type_free(&send_type);
         }
        else
         {
          index_in_dof_block_sparse_send[p] = 0;
          dof_number_sparse_send[p] = 0;
         }
       }
#endif
     }
   }

  /////////////////////////////////////////////////////////////////////////////
  // end of master block preconditioner only operations
  /////////////////////////////////////////////////////////////////////////////

  // compute the number of rows in each block

#ifdef PARANOID
  //check the vector is the correct length
  if (dof_to_block_map.size() != Internal_ndof_types)
   {
    std::ostringstream error_message;
    error_message
     << "The dof_to_block_map vector (size="
     << dof_to_block_map.size() << ") must be of size Internal_ndof_types="
     << Internal_ndof_types;
    throw OomphLibError(
                        error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // find the maximum block number
  unsigned max_block_number = 0;
  for (unsigned i = 0; i < Internal_ndof_types; i++)
   {
    if (dof_to_block_map[i] > max_block_number)
     {
      max_block_number = dof_to_block_map[i];
     }
   }

  // resize the storage the the block to dof map
  Block_number_to_dof_number_lookup.clear();
  Block_number_to_dof_number_lookup.resize(max_block_number+1);
  Ndof_in_block.clear();
  Ndof_in_block.resize(max_block_number+1);

  // resize storage
  Dof_number_to_block_number_lookup.resize(Internal_ndof_types);

  // build the storage for the two maps (block to dof) and (dof to block)
  for (unsigned i = 0; i < Internal_ndof_types; i++)
   {
    Dof_number_to_block_number_lookup[i] = dof_to_block_map[i];
    Block_number_to_dof_number_lookup[dof_to_block_map[i]].push_back(i);
    Ndof_in_block[dof_to_block_map[i]]++;
   }

#ifdef PARANOID
  // paranoid check that every block number has at least one DOF associated
  // with it
  for (unsigned i = 0; i < max_block_number+1; i++)
   {
    if (Block_number_to_dof_number_lookup[i].size() == 0)
     {
      std::ostringstream error_message;
      error_message  << "block number " << i
                     << " does not have any DOFs associated with it";
      throw OomphLibWarning(
                            error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif

  // update the number of blocks types
  Internal_nblock_types = max_block_number+1;

  // distributed or not depends on if we have more than one processor
  bool distributed = this->master_distribution_pt()->distributed();

  // create the new block distributions
  Block_distribution_pt.resize(Internal_nblock_types);
  for (unsigned i = 0; i < Internal_nblock_types; i++)
   {
    unsigned block_dim = 0;
    for (unsigned j = 0; j < Ndof_in_block[i]; j++)
     {
      block_dim +=
       dof_block_dimension(Block_number_to_dof_number_lookup[i][j]);
     }
    Block_distribution_pt[i] = new
     LinearAlgebraDistribution(comm_pt(),
                               block_dim,distributed);
   }

  // If blocks have been precomputed, compute the precomputed 
  // block distributions.
  if(Preconditioner_blocks_have_been_precomputed)
   {
    // Delete any existing distributions in Precomputed_block_distribution_pt.
    unsigned n_existing_precom_block_dist 
     = Precomputed_block_distribution_pt.size();
    for (unsigned dist_i = 0; dist_i < n_existing_precom_block_dist; dist_i++) 
     {
      delete Precomputed_block_distribution_pt[dist_i];
     }

    // Work out the distributions of the concatenated blocks.
    unsigned super_block_size = Block_to_block_map.size();
    Precomputed_block_distribution_pt.resize(super_block_size,0);
    for (unsigned super_block_i = 0; 
         super_block_i < super_block_size; super_block_i++)
     {
      unsigned sub_block_size = Block_to_block_map[super_block_i].size();
      Vector<LinearAlgebraDistribution*> tmp_dist_pt(sub_block_size,0);
      
      for (unsigned sub_block_i = 0; 
           sub_block_i < sub_block_size; sub_block_i++) 
       {
        tmp_dist_pt[sub_block_i] 
         = Replacement_block_pt(
             Block_to_block_map[super_block_i][sub_block_i],0)
               ->distribution_pt();
       }
  
      Precomputed_block_distribution_pt[super_block_i] 
       = new LinearAlgebraDistribution;
  
      LinearAlgebraDistributionHelpers::concatenate(
        tmp_dist_pt,*Precomputed_block_distribution_pt[super_block_i]);
     }
   }

  // create the distribution of the preconditioner matrix
  // if this preconditioner is a subsidiary preconditioner then it stored
  // at Distribution_pt.
  // if this preconditioner is a master preconditioner then it is stored
  // at Preconditioner_matrix_distribution_pt
  LinearAlgebraDistribution dist;
  LinearAlgebraDistributionHelpers::concatenate(Block_distribution_pt,dist);
  
  // build the distribution
  if (is_subsidiary_block_preconditioner())
   {
    this->build_distribution(dist);
   }
  else
   {
    Preconditioner_matrix_distribution_pt = new
     LinearAlgebraDistribution(dist);
   }
  
  // clearing up after comm to assemble sparse lookup schemes
#ifdef OOMPH_HAS_MPI
  if (send_requests_sparse.size()>0)
   {
    MPI_Waitall(send_requests_sparse.size(),
                &send_requests_sparse[0],MPI_STATUS_IGNORE);
   }
  if (recv_requests_sparse.size()>0)
   {
    MPI_Waitall(recv_requests_sparse.size(),
                &recv_requests_sparse[0],MPI_STATUS_IGNORE);
   }
  for (unsigned p = 0; p < nproc; p++)
   {
    delete[] index_in_dof_block_sparse_send[p];
    delete[] dof_number_sparse_send[p];
   }
  delete[] index_in_dof_block_sparse_send;
  delete[] dof_number_sparse_send;
  delete[] nreq_sparse;
  delete[] nreq_sparse_for_proc;
#endif

  // next we assemble the lookup schemes for the rows
  // if the matrix is not distributed then we assemble Global_index
  // if the matrix is distributed then Rows_to_send_..., Rows_to_recv_... etc
  if (!distributed)
   {
    // resize the storage
    Global_index.resize(Internal_nblock_types);
    for (unsigned b = 0; b < Internal_nblock_types; b++)
     {
      Global_index[b].resize(Block_distribution_pt[b]->nrow());
     }

    // compute
    unsigned nrow = this->master_nrow();
    for (unsigned i = 0; i < nrow; i++)
     {
      // the dof type number
      int dof_number = this->dof_number(i);
      if (dof_number >= 0)
       {

        // the block number
        unsigned block_number = Dof_number_to_block_number_lookup[dof_number];

        // compute the index in the block
        unsigned index_in_block = 0;
        unsigned ptr = 0;
        while (int(Block_number_to_dof_number_lookup[block_number][ptr])
               != dof_number)
         {
          index_in_block +=
           dof_block_dimension(Block_number_to_dof_number_lookup[block_number]
                               [ptr]);
          ptr++;
         }
        index_in_block += index_in_dof(i);
        Global_index[block_number][index_in_block] = i;
       }
     }
   }
  // otherwise the matrix is distributed
  else
   {
#ifdef OOMPH_HAS_MPI

    // the pointer to the master distribution
    const LinearAlgebraDistribution* master_distribution_pt =
     this->master_distribution_pt();

    // resize the nrows... storage
    Nrows_to_send_for_get_block.resize(Internal_nblock_types,nproc);
    Nrows_to_send_for_get_block.initialise(0);
    Nrows_to_send_for_get_ordered.resize(nproc);
    Nrows_to_send_for_get_ordered.initialise(0);

    // loop over my rows
    unsigned nrow_local = master_distribution_pt->nrow_local();
    unsigned first_row = master_distribution_pt->first_row();
    for (unsigned i = 0; i < nrow_local; i++)
     {

      // the block number
      int b = this->block_number(first_row + i);

      // check that the DOF i is associated with this preconditioner
      if (b >= 0)
       {
        // the block index
        unsigned j = this->index_in_block(first_row + i);

        // the processor this row will be sent to
        unsigned block_p = 0;
        while(!(Block_distribution_pt[b]->first_row(block_p) <= j &&
                (Block_distribution_pt[b]->first_row(block_p) +
                 Block_distribution_pt[b]->nrow_local(block_p) > j)))
         {
          block_p++;
         }

        // and increment the counter
        Nrows_to_send_for_get_block(b,block_p)++;
        Nrows_to_send_for_get_ordered[block_p]++;
       }
     }

    // resize the storage for Nrows_to_recv
    Nrows_to_recv_for_get_block.resize(Internal_nblock_types,nproc);
    Nrows_to_recv_for_get_block.initialise(0);
    Nrows_to_recv_for_get_ordered.resize(nproc);
    Nrows_to_recv_for_get_ordered.initialise(0);

    // next we send the number of rows that will be sent by this processor
    Vector<unsigned*> nrows_to_send(nproc,0);
    Vector<unsigned*> nrows_to_recv(nproc,0);
    Vector<MPI_Request> send_requests_nrow;
    Vector<MPI_Request> recv_requests_nrow;
    Vector<unsigned> proc;
    for (unsigned p = 0; p < nproc; p++)
     {
      if (p != my_rank)
       {
        // send
        proc.push_back(p);
        nrows_to_send[p] = new unsigned[Internal_nblock_types];
        for (unsigned b = 0; b < Internal_nblock_types; b++)
         {
          nrows_to_send[p][b] =
           Nrows_to_send_for_get_block(b,p);
         }
        MPI_Request s_req;
        MPI_Isend(nrows_to_send[p],Internal_nblock_types,MPI_UNSIGNED,p,3,
                  comm_pt()->mpi_comm(),&s_req);
        send_requests_nrow.push_back(s_req);

        // recv
        nrows_to_recv[p] = new unsigned[Internal_nblock_types];
        MPI_Request r_req;
        MPI_Irecv(nrows_to_recv[p],Internal_nblock_types,MPI_UNSIGNED,p,3,
                  comm_pt()->mpi_comm(),&r_req);
        recv_requests_nrow.push_back(r_req);
       }
      // send to self
      else
       {
        for (unsigned b = 0; b < Internal_nblock_types; b++)
         {
          Nrows_to_recv_for_get_block(b,p) =
           Nrows_to_send_for_get_block(b,p);
         }
        Nrows_to_recv_for_get_ordered[p] = Nrows_to_send_for_get_ordered[p];
       }
     }

    // create some temporary storage for the global row indices that will
    // be received from another processor.
    DenseMatrix<int*> block_rows_to_send(Internal_nblock_types,nproc,0);
    Vector<int*> ordered_rows_to_send(nproc,0);

    // resize the rows... storage
    Rows_to_send_for_get_block.resize(Internal_nblock_types,nproc);
    Rows_to_send_for_get_block.initialise(0);
    Rows_to_send_for_get_ordered.resize(nproc);
    Rows_to_send_for_get_ordered.initialise(0);
    Rows_to_recv_for_get_block.resize(Internal_nblock_types,nproc);
    Rows_to_recv_for_get_block.initialise(0);

    // resize the storage
    for (unsigned p = 0; p < nproc; p++)
     {
      for (unsigned b = 0; b < Internal_nblock_types; b++)
       {
        Rows_to_send_for_get_block(b,p)
         = new int[Nrows_to_send_for_get_block(b,p)];
        if (p != my_rank)
         {
          block_rows_to_send(b,p)
           = new int[Nrows_to_send_for_get_block(b,p)];
         }
        else
         {
          Rows_to_recv_for_get_block(b,p)
           = new int[Nrows_to_send_for_get_block(b,p)];
         }
       }
      Rows_to_send_for_get_ordered[p]
       = new int [Nrows_to_send_for_get_ordered[p]];
     }



    // loop over my rows to allocate the nrows
    DenseMatrix<unsigned> ptr_block(Internal_nblock_types,nproc,0);
    for (unsigned i = 0; i < nrow_local; i++)
     {
      // the block number
      int b = this->block_number(first_row + i);

      // check that the DOF i is associated with this preconditioner
      if (b >= 0)
       {

        // the block index
        unsigned j = this->index_in_block(first_row + i);

        // the processor this row will be sent to
        unsigned block_p = 0;
        while(!(Block_distribution_pt[b]->first_row(block_p) <= j &&
                (Block_distribution_pt[b]->first_row(block_p) +
                 Block_distribution_pt[b]->nrow_local(block_p) > j)))
         {
          block_p++;
         }

        // and store the row
        Rows_to_send_for_get_block(b,block_p)[ptr_block(b,block_p)] = i;
        if (block_p != my_rank)
         {
          block_rows_to_send(b,block_p)[ptr_block(b,block_p)]
           = j - Block_distribution_pt[b]->first_row(block_p);
         }
        else
         {
          Rows_to_recv_for_get_block(b,block_p)[ptr_block(b,block_p)]
           = j - Block_distribution_pt[b]->first_row(block_p);
         }
        ptr_block(b,block_p)++;
       }
     }

    // next block ordered
    for (unsigned p = 0; p < nproc; ++p)
     {
      int pt = 0;
      for (unsigned b = 0; b < Internal_nblock_types; ++b)
       {

        for (unsigned i = 0; i < Nrows_to_send_for_get_block(b,p); ++i)
         {
          Rows_to_send_for_get_ordered[p][pt] =
           Rows_to_send_for_get_block(b,p)[i];
          pt++;
         }
       }
     }

    // next process the nrow recvs as they complete

    // recv and store the data
    unsigned c = recv_requests_nrow.size();
    while (c > 0)
     {

      // wait for any communication to finish
      int req_number;
      MPI_Waitany(c,&recv_requests_nrow[0],&req_number,MPI_STATUS_IGNORE);
      recv_requests_nrow.erase(recv_requests_nrow.begin()+req_number);
      c--;

      // determine the source processor
      unsigned p = proc[req_number];
      proc.erase(proc.begin()+req_number);

      // copy the data to its final storage
      Nrows_to_recv_for_get_ordered[p]=0;
      for (unsigned b = 0; b < Internal_nblock_types; b++)
       {
        Nrows_to_recv_for_get_block(b,p) = nrows_to_recv[p][b];
        Nrows_to_recv_for_get_ordered[p] += nrows_to_recv[p][b];
       }

      // and clear
      delete[] nrows_to_recv[p];
     }

    // resize the storage for the incoming rows data
    Rows_to_recv_for_get_ordered.resize(nproc,0);
    for (unsigned p = 0; p < nproc; p++)
     {
      if (p != my_rank)
       {
        for (unsigned b = 0; b < Internal_nblock_types; b++)
         {
          Rows_to_recv_for_get_block(b,p)
           = new int[Nrows_to_recv_for_get_block(b,p)];
         }
       }
     }

    // compute the number of sends and recv from this processor
    // to each other processor
    Vector<unsigned> nsend_for_rows(nproc,0);
    Vector<unsigned> nrecv_for_rows(nproc,0);
    for (unsigned p = 0; p < nproc; p++)
     {
      if (p != my_rank)
       {
        for (unsigned b = 0; b < Internal_nblock_types; b++)
         {
          if (Nrows_to_send_for_get_block(b,p) > 0)
           {
            nsend_for_rows[p]++;
           }
          if (Nrows_to_recv_for_get_block(b,p) > 0)
           {
            nrecv_for_rows[p]++;
           }
         }
       }
     }

    // finally post the sends and recvs
    MPI_Aint base_displacement;
    MPI_Address(matrix_pt(),&base_displacement);
    Vector<MPI_Request> req_rows;
    for (unsigned p = 0; p < nproc; p++)
     {
      if (p != my_rank)
       {
        // send
        if (nsend_for_rows[p] > 0)
         {
          MPI_Datatype send_types[nsend_for_rows[p]];
          MPI_Aint send_displacements[nsend_for_rows[p]];
          int send_sz[nsend_for_rows[p]];
          unsigned send_ptr = 0;
          for (unsigned b = 0; b < Internal_nblock_types; b++)
           {
            if (Nrows_to_send_for_get_block(b,p) > 0)
             {
              MPI_Type_contiguous(Nrows_to_send_for_get_block(b,p),
                                  MPI_INT,&send_types[send_ptr]);
              MPI_Type_commit(&send_types[send_ptr]);
              MPI_Address(block_rows_to_send(b,p),
                          &send_displacements[send_ptr]);
              send_displacements[send_ptr] -= base_displacement;
              send_sz[send_ptr] = 1;
              send_ptr++;
             }
           }
          MPI_Datatype final_send_type;
          MPI_Type_struct(nsend_for_rows[p],send_sz,send_displacements,
                          send_types,&final_send_type);
          MPI_Type_commit(&final_send_type);
          for (unsigned i = 0; i < nsend_for_rows[p]; i++)
           {
            MPI_Type_free(&send_types[i]);
           }
          MPI_Request send_req;
          MPI_Isend(matrix_pt(),1,final_send_type,p,4,
                    comm_pt()->mpi_comm(),&send_req);
          req_rows.push_back(send_req);
          MPI_Type_free(&final_send_type);
         }

        // recv
        if (nrecv_for_rows[p] > 0)
         {
          MPI_Datatype recv_types[nrecv_for_rows[p]];
          MPI_Aint recv_displacements[nrecv_for_rows[p]];
          int recv_sz[nrecv_for_rows[p]];
          unsigned recv_ptr = 0;
          for (unsigned b = 0; b < Internal_nblock_types; b++)
           {
            if (Nrows_to_recv_for_get_block(b,p) > 0)
             {
              MPI_Type_contiguous(Nrows_to_recv_for_get_block(b,p),
                                  MPI_INT,&recv_types[recv_ptr]);
              MPI_Type_commit(&recv_types[recv_ptr]);
              MPI_Address(Rows_to_recv_for_get_block(b,p),
                          &recv_displacements[recv_ptr]);
              recv_displacements[recv_ptr] -= base_displacement;
              recv_sz[recv_ptr] = 1;
              recv_ptr++;
             }
           }
          MPI_Datatype final_recv_type;
          MPI_Type_struct(nrecv_for_rows[p],recv_sz,recv_displacements,
                          recv_types,&final_recv_type);
          MPI_Type_commit(&final_recv_type);
          for (unsigned i = 0; i < nrecv_for_rows[p]; i++)
           {
            MPI_Type_free(&recv_types[i]);
           }
          MPI_Request recv_req;
          MPI_Irecv(matrix_pt(),1,final_recv_type,p,4,
                    comm_pt()->mpi_comm(),&recv_req);
          req_rows.push_back(recv_req);
          MPI_Type_free(&final_recv_type);
         }
       }
     }

    // cleaning up Waitalls


    // wait for the recv requests so we can compute
    // Nrows_to_recv_for_get_ordered
    unsigned n_req_rows = req_rows.size();
    if (n_req_rows)
     {
      MPI_Waitall(n_req_rows,&req_rows[0],MPI_STATUS_IGNORE);
     }

    // resize the storage
    Rows_to_recv_for_get_ordered.resize(nproc);
    Rows_to_recv_for_get_ordered.initialise(0);

    // construct block offset
    Vector<int> vec_offset(Internal_nblock_types,0);
    for (unsigned b = 1; b < Internal_nblock_types; ++b)
     {
      vec_offset[b]=vec_offset[b-1]+Block_distribution_pt[b-1]->nrow_local();
     }

    //
    for (unsigned p = 0; p < nproc; p++)
     {
      int pt = 0;
      Rows_to_recv_for_get_ordered[p]
       = new int[Nrows_to_recv_for_get_ordered[p]];
      for (unsigned b = 0; b < Internal_nblock_types; b++)
       {
        for (unsigned i = 0; i < Nrows_to_recv_for_get_block(b,p); i++)
         {
          Rows_to_recv_for_get_ordered[p][pt] =
           Rows_to_recv_for_get_block(b,p)[i]+vec_offset[b];
          pt++;
         }
       }
     }

    // clean up
    for (unsigned p = 0; p < nproc; p++)
     {
      if (p!= my_rank)
       {
        for (unsigned b = 0; b < Internal_nblock_types; b++)
         {
          delete[] block_rows_to_send(b,p);
         }
        if (Nrows_to_send_for_get_ordered[p] > 0)
         {
          delete[] ordered_rows_to_send[p];
         }
       }
     }

    // and the send reqs
    unsigned n_req_send_nrow = send_requests_nrow.size();
    if (n_req_send_nrow)
     {
      MPI_Waitall(n_req_send_nrow,&send_requests_nrow[0],MPI_STATUS_IGNORE);
     }
    for (unsigned p = 0; p < nproc; p++)
     {
      delete[] nrows_to_send[p];
     }
#endif
   }

  // If we asked for output of blocks to a file then do it.
  if(block_output_on())
   output_blocks_to_files(Output_base_filename);
 }

 //============================================================================
 //??ds
 /// \short Function to turn this preconditioner into a
 /// subsidiary preconditioner that operates within a bigger
 /// "master block preconditioner (e.g. a Navier-Stokes 2x2 block
 /// preconditioner dealing with the fluid sub-blocks within a
 /// 3x3 FSI preconditioner. Once this is done the master block
 /// preconditioner deals with the block setup etc.
 /// The vector block_map must specify the dof number in the
 /// master preconditioner that corresponds to a block number in this
 /// preconditioner. ??ds horribly misleading comment!
 /// The length of the vector is used to determine the number of
 /// blocks in this preconditioner therefore it must be correctly sized.
 /// This calls the other turn_into_subsidiary_block_preconditioner(...)
 /// function providing an empty doftype_to_doftype_map vector.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 turn_into_subsidiary_block_preconditioner
 (BlockPreconditioner<MATRIX>* master_block_prec_pt,
  Vector<unsigned>& block_map)
 {
   Vector<Vector<unsigned> > doftype_to_doftype_map;
   Preconditioner_doftypes_have_been_coarsened = false;
   turn_into_subsidiary_block_preconditioner(master_block_prec_pt,
       block_map, doftype_to_doftype_map);
 }


 //============================================================================
 //??ds
 /// \short Function to turn this preconditioner into a
 /// subsidiary preconditioner that operates within a bigger
 /// "master block preconditioner (e.g. a Navier-Stokes 2x2 block
 /// preconditioner dealing with the fluid sub-blocks within a
 /// 3x3 FSI preconditioner. Once this is done the master block
 /// preconditioner deals with the block setup etc.
 /// The vector block_map must specify the dof number in the
 /// master preconditioner that corresponds to a block number in this
 /// preconditioner. ??ds horribly misleading comment!
 /// The length of the vector is used to determine the number of
 /// blocks in this preconditioner therefore it must be correctly sized.
 /// The Vector doftype_to_doftype_map is provided if preconditioner 
 /// DOF types needs to be coarsened. 
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 turn_into_subsidiary_block_preconditioner
 (BlockPreconditioner<MATRIX>* master_block_prec_pt,
  Vector<unsigned>& block_map,
  Vector<Vector<unsigned> > &doftype_to_doftype_map)
 {
  
  if(doftype_to_doftype_map.size() != 0)
  {
#ifdef PARANOID
   // Checks for doftype_to_doftype_map
   // No more than ndof types described, 
   // and check that all entries are unique.
   std::set<unsigned> doftype_map_set;

   unsigned doftype_to_doftype_map_size = doftype_to_doftype_map.size();
   for (unsigned i = 0; i < doftype_to_doftype_map_size; i++)
    {
     unsigned doftype_to_doftype_map_i_size = doftype_to_doftype_map[i].size();
     for (unsigned j = 0; j < doftype_to_doftype_map_i_size; j++)
      {
       std::set<unsigned>::iterator doftype_map_it;
       std::pair<std::set<unsigned>::iterator,bool> doftype_map_ret;

       doftype_map_ret = doftype_map_set.insert(doftype_to_doftype_map[i][j]);
         
       if(!doftype_map_ret.second)
        {
         std::ostringstream error_message;
         error_message << "Error: the doftype number "
                       << doftype_to_doftype_map[i][j]
                       << " is already inserted."
                       << std::endl;
         throw OomphLibError(error_message.str(),
                             OOMPH_CURRENT_FUNCTION,
                             OOMPH_EXCEPTION_LOCATION);
        }
      }
    }
     
//   // All doftype described in doftype_to_doftype_map must be unique.
//   if(precomputed_block_nrow != doftype_map_set.size())
//    {
//     std::ostringstream error_message;
//     error_message << "Error: all doftypes must be assigned. \n"
//                   << "Only " << doftype_map_set.size()
//                   << " doftypes have been assigned."
//                   << std::endl;
//     throw OomphLibError(error_message.str(),
//                         OOMPH_CURRENT_FUNCTION,
//                         OOMPH_EXCEPTION_LOCATION);
//    }

#endif
    // We want to coarsen the preconditioner doftypes.
    Preconditioner_doftypes_have_been_coarsened = true;

    // Set the Doftype_to_doftype_map.
    Doftype_to_doftype_map = doftype_to_doftype_map;
  }
  else
  {
    Preconditioner_doftypes_have_been_coarsened = false;
  }

  // Set the master block preconditioner pointer
  Master_block_preconditioner_pt = master_block_prec_pt;

  // If the master block preconditioner has precomputed blocks
  // this preconditioner use these precomputed blocks.
  if(master_block_prec_pt->preconditioner_blocks_have_been_precomputed())
   {
    // We have precomputed preconditioner blocks.
    Preconditioner_blocks_have_been_precomputed = true;

    // Store the precomputed blocks.
    Replacement_block_pt
      = master_block_prec_pt->replacement_block_pt();
    
    // Only store the master's Doftype_to_doftype_map which is relevant
    // to this preconditioner.
    unsigned tmp_ndof_types_size = block_map.size();
    Vector<unsigned> tmp_dof_type_vec;
    Doftype_to_doftype_map.resize(block_map.size());
    for (unsigned doftype_i = 0; doftype_i < tmp_ndof_types_size; doftype_i++) 
     {
      Vector<unsigned> tmp_dof_type_subvec 
        = master_block_prec_pt->coarse_dof_type_subvec(doftype_i);
      Doftype_to_doftype_map[doftype_i] = tmp_dof_type_subvec;

      unsigned tmp_dof_type_subvec_size = tmp_dof_type_subvec.size();

      for (unsigned fine_dof_type_i = 0; 
           fine_dof_type_i < tmp_dof_type_subvec_size; fine_dof_type_i++) 
       {
        tmp_dof_type_vec.push_back(tmp_dof_type_subvec[fine_dof_type_i]);
       }
     }

    // Fill the block_map with the most fine grain dof types.
    block_map = tmp_dof_type_vec;
    std::sort(block_map.begin(),block_map.end());
      
    // Set the mapping from the master preconditioner blocks to the
    // subsidiary preconditioner blocks.
    Dof_number_in_master_preconditioner = block_map;
    
    // Get the number of block types (and dof types) in this preconditioner
    // from the length of the block_map vector.
    Internal_ndof_types = block_map.size();
    Internal_nblock_types = Internal_ndof_types;
   }
  // If the master block preconditioner does not have precomputed
  // preconditioner blocks.
  else
   {
    // Set the mapping from the master preconditioner blocks to the
    // subsidiary preconditioner blocks.
    Dof_number_in_master_preconditioner = block_map;
    
    // Get the number of block types (and dof types) in this preconditioner
    // from the length of the block_map vector.
    Internal_ndof_types = block_map.size();
    Internal_nblock_types = Internal_ndof_types;
   }
 } // end of turn_into_subsidiary_block_preconditioner(...)


 //============================================================================
 /// Determine the size of the matrix blocks and setup the
 /// lookup schemes relating the global degrees of freedom with
 /// their "blocks" and their indices (row/column numbers) in those
 /// blocks.
 /// The distributions of the preconditioner and the blocks are
 /// automatically specified (and assumed to be uniform) at this
 /// stage.
 /// This method should be used if each DOF type corresponds to a
 /// unique block type.
 //============================================================================
 template<typename MATRIX>
 void BlockPreconditioner<MATRIX>::block_setup()
 {
  // Get the number of dof types.
  unsigned internal_n_dof_types = internal_ndof_types();

  // Build the dof to block map - assume that each type of dof corresponds
  // to a different type of block.
  Vector<unsigned> dof_to_block_lookup(internal_n_dof_types);
  for (unsigned i = 0; i < internal_n_dof_types; i++)
   {
    dof_to_block_lookup[i] = i;
   }

  // call the block setup method
  this->block_setup(dof_to_block_lookup);
 }


 //============================================================================
 /// Get the block matrices required for the block preconditioner. Takes a
 /// pointer to a matrix of bools that indicate if a specified sub-block is
 /// required for the preconditioning operation. Computes the required block
 /// matrices, and stores pointers to them in the matrix block_matrix_pt. If an
 /// entry in block_matrix_pt is equal to NULL that sub-block has not been
 /// requested and is therefore not available.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 get_blocks(DenseMatrix<bool>& required_blocks,
            DenseMatrix<MATRIX*>& block_matrix_pt) const
 {

  // Cache number of block types
  const unsigned n_block_types=this->Internal_nblock_types;

#ifdef PARANOID
  // If required blocks matrix pointer is not the correct size then abort.
  if ((required_blocks.nrow() != n_block_types) ||
      (required_blocks.ncol() != n_block_types))
   {

    std::ostringstream error_message;
    error_message << "The size of the matrix of bools required_blocks"
                  << "(which indicates which blocks are required) is not the "
                  << "right size, required_blocks is "
                  << required_blocks.ncol()
                  << " x " << required_blocks.nrow() << ", whereas it should "
                  << "be " << n_block_types << " x " << n_block_types;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  // If block matrix pointer is not the correct size then abort.
  if ((block_matrix_pt.nrow() != n_block_types) ||
      (block_matrix_pt.ncol() != n_block_types))
   {
    std::ostringstream error_message;
    error_message << "The size of the block matrix pt is not the "
                  << "right size, block_matrix_pt is "
                  << block_matrix_pt.ncol()
                  << " x " << block_matrix_pt.nrow() << ", whereas it should "
                  << "be " << n_block_types << " x " << n_block_types;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

#endif

  // Loop over the blocks
  for (unsigned i = 0; i < n_block_types; i++)
   {
    for (unsigned j = 0; j < n_block_types; j++)
     {
      // If block(i,j) is required then create a matrix and fill it in.
      if (required_blocks(i,j))
       {
        //??ds might want to remove this use of new as well?
        block_matrix_pt(i,j) = new MATRIX;
        get_block(i, j, *block_matrix_pt(i,j));
       }

      // Otherwise set pointer to null.
      else
       {
        block_matrix_pt(i,j) = 0;
       }
     }
   }
 }

 //============================================================================
 /// \short Takes the naturally ordered vector and rearranges it into a
 /// vector of sub vectors corresponding to the blocks, so s[b][i] contains
 /// the i-th entry in the vector associated with block b.
 /// Note: If the preconditioner is a subsidiary preconditioner then only the
 /// sub-vectors associated with the blocks of the subsidiary preconditioner
 /// will be included. Hence the length of v is master_nrow() whereas the
 /// total length of the s s vectors is Nrow.
 /// NOTE: If the preconditioner blocks are precomputed, then this function
 /// calls get_block_vectors_with_precomputed_block_ordering(...),
 /// otherwise get_block_vector_with_original_matrix_ordering(...) is called.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 get_block_vectors(const DoubleVector& v, Vector<DoubleVector >& s) const
 {
  // If preconditioner_blocks have been precomputed then we call
  // get_block_vectors_with_precomputed_block_ordering(...) to ensure that
  // the distribution of the block vectors s are the same as the 
  // precomputed block distributions.
  if(Preconditioner_blocks_have_been_precomputed)
   {
    get_block_vectors_with_precomputed_block_ordering(v,s);
   }
  else
   // Otherwise the n-th block matrix came from the original matrix, so we
   // call get_block_vectors_with_original_matrix_ordering(...)
   {
    get_block_vectors_with_original_matrix_ordering(v,s);
   }
 }

 //============================================================================
 /// \short A helper function, takes the naturally ordered vector and 
 /// rearranges it into a vector of sub vectors corresponding to the blocks, 
 /// so s[b][i] contains the i-th entry in the vector associated with block b. 
 /// These blocks and vectors are those corresponding to the original block 
 /// matrix ordering, i.e. there are no precomputed blocks.
 /// Note: If the preconditioner is a subsidiary preconditioner then only the
 /// sub-vectors associated with the blocks of the subsidiary preconditioner
 /// will be included. Hence the length of v is master_nrow() whereas the
 /// total length of the s s vectors is Nrow.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 get_block_vectors_with_original_matrix_ordering(
                                                 const DoubleVector& v, Vector<DoubleVector >& s) const
 {
#ifdef PARANOID
  if (!v.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must match the "
                  << " specified master_distribution_pt(). \n"
                  << "i.e. Distribution_pt in the master preconditioner";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Number of block types
  const unsigned nblock = this->internal_nblock_types(true);

  // if + only one processor
  //    + more than one processor but matrix_pt is not distributed
  // then use the serial get_block method
  if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
      !this->distribution_pt()->distributed())
   {

    // Vector of vectors for each section of residual vector
    s.resize(nblock);

    // pointer to the data in v
    const double* v_pt = v.values_pt();

    // setup the block vector and then insert the data
    for (unsigned b = 0; b < nblock; b++)
     {
      s[b].build(Block_distribution_pt[b],0.0);
      double* s_pt = s[b].values_pt();
      unsigned nrow = s[b].nrow();
      for (unsigned i = 0; i < nrow; i++)
       {
        s_pt[i] = v_pt[this->Global_index[b][i]];
       }
     }
   }
  // otherwise use mpi
  else
   {
#ifdef OOMPH_HAS_MPI
    // my rank
    unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

    // the number of processors
    unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

    // build the vectors
    s.resize(nblock);
    for (unsigned b = 0; b < nblock; b++)
     {
      s[b].build(Block_distribution_pt[b],0.0);
     }

    // determine the maximum number of rows to be sent or recv
    // and determine the number of blocks each processor will send and recv
    // communication for
    Vector<int> nblock_send(nproc,0);
    Vector<int> nblock_recv(nproc,0);
    unsigned max_n_send_or_recv = 0;
    for (unsigned p = 0; p < nproc; p++)
     {
      for (unsigned b = 0; b < nblock; b++)
       {
        max_n_send_or_recv =
         std::max(max_n_send_or_recv,Nrows_to_send_for_get_block(b,p));
        max_n_send_or_recv =
         std::max(max_n_send_or_recv,Nrows_to_recv_for_get_block(b,p));
        if (Nrows_to_send_for_get_block(b,p) > 0)
         {
          nblock_send[p]++;
         }
        if (Nrows_to_recv_for_get_block(b,p) > 0)
         {
          nblock_recv[p]++;
         }
       }
     }

    // create a vectors of 1s the size of the nblock for the mpi indexed
    // data types
    int* block_lengths = new int[max_n_send_or_recv];
    for (unsigned i = 0; i < max_n_send_or_recv; i++)
     {
      block_lengths[i] = 1;
     }

    // perform the sends and receives
    Vector<MPI_Request> requests;
    for (unsigned p = 0; p < nproc; p++)
     {
      // send and recv with other processors
      if (p != my_rank)
       {
        // send
        if (nblock_send[p] > 0)
         {
          // create the datatypes vector
          MPI_Datatype block_send_types[nblock_send[p]];

          // create the datatypes
          unsigned ptr = 0;
          for (unsigned b = 0; b < nblock; b++)
           {
            if (Nrows_to_send_for_get_block(b,p) > 0)
             {
              MPI_Type_indexed(Nrows_to_send_for_get_block(b,p),block_lengths,
                               Rows_to_send_for_get_block(b,p),MPI_DOUBLE,
                               &block_send_types[ptr]);
              MPI_Type_commit(&block_send_types[ptr]);
              ptr++;
             }
           }

          // compute the displacements and lengths
          MPI_Aint displacements[nblock_send[p]];
          int lengths[nblock_send[p]];
          for (int i = 0; i < nblock_send[p]; i++)
           {
            lengths[i] = 1;
            displacements[i] = 0;
           }

          // build the final datatype
          MPI_Datatype type_send;
          MPI_Type_struct(nblock_send[p],lengths,displacements,
                          block_send_types,&type_send);
          MPI_Type_commit(&type_send);

          // send
          MPI_Request send_req;
          MPI_Isend(const_cast<double*>(v.values_pt()),1,type_send,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &send_req);
          MPI_Type_free(&type_send);
          for (int i = 0; i < nblock_send[p]; i++)
           {
            MPI_Type_free(&block_send_types[i]);
           }
          requests.push_back(send_req);
         }

        // recv
        if (nblock_recv[p] > 0)
         {
          // create the datatypes vector
          MPI_Datatype block_recv_types[nblock_recv[p]];

          // and the displacements
          MPI_Aint displacements[nblock_recv[p]];

          // and the lengths
          int lengths[nblock_recv[p]];

          // all displacements are computed relative to s[0] values
          MPI_Aint displacements_base;
          MPI_Address(s[0].values_pt(),&displacements_base);

          // now build
          unsigned ptr = 0;
          for (unsigned b = 0; b < nblock; b++)
           {
            if (Nrows_to_recv_for_get_block(b,p) > 0)
             {
              MPI_Type_indexed(Nrows_to_recv_for_get_block(b,p),block_lengths,
                               Rows_to_recv_for_get_block(b,p),MPI_DOUBLE,
                               &block_recv_types[ptr]);
              MPI_Type_commit(&block_recv_types[ptr]);
              MPI_Address(s[b].values_pt(),&displacements[ptr]);
              displacements[ptr] -= displacements_base;
              lengths[ptr] = 1;
              ptr++;
             }
           }

          // build the final data type
          MPI_Datatype type_recv;
          MPI_Type_struct(nblock_recv[p],lengths,displacements,
                          block_recv_types,&type_recv);
          MPI_Type_commit(&type_recv);

          // recv
          MPI_Request recv_req;
          MPI_Irecv(s[0].values_pt(),1,type_recv,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &recv_req);
          MPI_Type_free(&type_recv);
          for (int i = 0; i < nblock_recv[p]; i++)
           {
            MPI_Type_free(&block_recv_types[i]);
           }
          requests.push_back(recv_req);
         }
       }

      // communicate wih self
      else
       {
        const double* v_values_pt = v.values_pt();
        for (unsigned b = 0; b < nblock; b++)
         {
          double* w_values_pt = s[b].values_pt();
          for (unsigned i = 0; i < Nrows_to_send_for_get_block(b,p); i++)
           {
            w_values_pt[Rows_to_recv_for_get_block(b,p)[i]] =
             v_values_pt[Rows_to_send_for_get_block(b,p)[i]];
           }
         }
       }
     }

    // and then just wait
    unsigned c = requests.size();
    Vector<MPI_Status> stat(c);
    if (c)
     {
      MPI_Waitall(c,&requests[0],&stat[0]);
     }
    delete[] block_lengths;

#else
    // throw error
    std::ostringstream error_message;
    error_message << "The preconditioner is distributed and on more than one "
                  << "processor. MPI is required.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
#endif
   }
 }

 //============================================================================
 /// \short A helper function, takes the naturally ordered vector and 
 /// rearranges it into a vector of sub vectors corresponding to the 
 /// PRECOMPUTED blocks, so s[b][i] contains the i-th entry in the vector 
 /// associated with the precomputed block b.
 /// Note: If the preconditioner is a subsidiary preconditioner then only the
 /// sub-vectors associated with the blocks of the subsidiary preconditioner
 /// will be included. Hence the length of v is master_nrow() whereas the
 /// total length of the s vectors is Nrow. 
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 get_block_vectors_with_precomputed_block_ordering(
                                                   const DoubleVector& v, Vector<DoubleVector >& s) const
 {
#ifdef PARANOID
  if (!v.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must match the "
                  << " specified master_distribution_pt(). \n"
                  << "i.e. Distribution_pt in the master preconditioner";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if(!Preconditioner_blocks_have_been_precomputed)
   {
    std::ostringstream error_message;
    error_message << "You have not set precomputed blocks. It does not make "
                  << "sense to get vectors with precomputed block ordering.\n";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  const unsigned nprecomputed_blocks = Block_to_block_map.size();
  s.resize(nprecomputed_blocks);
  
  for (unsigned b_i = 0; b_i < nprecomputed_blocks; b_i++) 
   {
    get_block_vector_with_precomputed_block_ordering(b_i,v,s[b_i]);
   }

 }

 //============================================================================
 /// \short Takes the vector of block vectors, s, and copies its entries into
 /// the naturally ordered vector, v. If this is a subsidiary block
 /// preconditioner only those entries in v that are associated with its
 /// blocks are affected.
 /// NOTE: If the preconditioner blocks are precomputed, then this function
 /// calls return_block_vectors_with_precomputed_block_ordering(...),
 /// otherwise return_block_vector_with_original_matrix_ordering(...) 
 /// is called.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 return_block_vectors(const Vector<DoubleVector >& s, DoubleVector& v) const
 {
  // If the preconditioner blocks have been precomputed then it is assumed
  // that the block vectors s have the precomputed block ordering i.e. they 
  // have the same distributions as the block vectors from
  // get_block_vectors_with_precomputed_block_ordering(...). We call
  // return_block_vectors_with_precomputed_block_ordering(...)  ensure that
  // the entries are returned in the correct place.
  if(Preconditioner_blocks_have_been_precomputed)
   {
    return_block_vectors_with_precomputed_block_ordering(s,v);
   }
  else
   // Otherwise we assume that the block vectors s have the same distribution
   // as Block_distribution_pt.
   {
    return_block_vectors_with_original_matrix_ordering(s,v);
   }
 }


 //============================================================================
 /// \short A helper function, takes the vector of block vectors, s, and 
 /// copies its entries into the naturally ordered vector, v. 
 /// The block vectors are assumed to have the ordering of the original 
 /// block matrices. I.e. there are no precomputed blocks. 
 /// If this is a subsidiary block preconditioner only those entries in v 
 /// that are associated with its blocks are affected.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 return_block_vectors_with_original_matrix_ordering(
                                                    const Vector<DoubleVector >& s, DoubleVector& v) const
 {
  // the number of blocks
  unsigned nblock = this->internal_nblock_types(true);

#ifdef PARANOID
  if (!v.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must match the "
                  << " specified master_distribution_pt(). \n"
                  << "i.e. Distribution_pt in the master preconditioner";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  for (unsigned b = 0; b < nblock; b++)
   {
    if (!s[b].built())
     {
      std::ostringstream error_message;
      error_message << "The distribution of the block vector " << b
                    << " must be setup.";
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
    if (*(s[b].distribution_pt()) != *(Block_distribution_pt[b]))
     {
      std::ostringstream error_message;
      error_message << "The distribution of the block vector " << b
                    << "must match the"
                    << " specified distribution at Block_distribution_pt["
                    << b << "]";
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif

  // if + only one processor
  //    + more than one processor but matrix_pt is not distributed
  // then use the serial get_block method
  if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
      !this->distribution_pt()->distributed())
   {
    double* v_pt = v.values_pt();
    for (unsigned b = 0; b < nblock; b++)
     {
      const double* s_pt = s[b].values_pt();
      unsigned nrow = this->block_dimension(b);
      for (unsigned i = 0; i < nrow; i++)
       {
        v_pt[this->Global_index[b][i]] = s_pt[i];
       }
     }
   }
  // otherwise use mpi
  else
   {
#ifdef OOMPH_HAS_MPI

    // my rank
    unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

    // the number of processors
    unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

    // determine the maximum number of rows to be sent or recv
    // and determine the number of blocks each processor will send and recv
    // communication for
    Vector<int> nblock_send(nproc,0);
    Vector<int> nblock_recv(nproc,0);
    unsigned max_n_send_or_recv = 0;
    for (unsigned p = 0; p < nproc; p++)
     {
      for (unsigned b = 0; b < nblock; b++)
       {
        max_n_send_or_recv =
         std::max(max_n_send_or_recv,Nrows_to_send_for_get_block(b,p));
        max_n_send_or_recv =
         std::max(max_n_send_or_recv,Nrows_to_recv_for_get_block(b,p));
        if (Nrows_to_send_for_get_block(b,p) > 0)
         {
          nblock_recv[p]++;
         }
        if (Nrows_to_recv_for_get_block(b,p) > 0)
         {
          nblock_send[p]++;
         }
       }
     }

    // create a vectors of 1s the size of the nblock for the mpi indexed
    // data types
    int* block_lengths = new int[max_n_send_or_recv];
    for (unsigned i = 0; i < max_n_send_or_recv; i++)
     {
      block_lengths[i] = 1;
     }

    // perform the sends and receives
    Vector<MPI_Request> requests;
    for (unsigned p = 0; p < nproc; p++)
     {
      // send and recv with other processors
      if (p != my_rank)
       {
        // recv
        if (nblock_recv[p] > 0)
         {
          // create the datatypes vector
          MPI_Datatype block_recv_types[nblock_recv[p]];

          // create the datatypes
          unsigned ptr = 0;
          for (unsigned b = 0; b < nblock; b++)
           {
            if (Nrows_to_send_for_get_block(b,p) > 0)
             {
              MPI_Type_indexed(Nrows_to_send_for_get_block(b,p),block_lengths,
                               Rows_to_send_for_get_block(b,p),MPI_DOUBLE,
                               &block_recv_types[ptr]);
              MPI_Type_commit(&block_recv_types[ptr]);
              ptr++;
             }
           }

          // compute the displacements and lengths
          MPI_Aint displacements[nblock_recv[p]];
          int lengths[nblock_recv[p]];
          for (int i = 0; i < nblock_recv[p]; i++)
           {
            lengths[i] = 1;
            displacements[i] = 0;
           }

          // build the final datatype
          MPI_Datatype type_recv;
          MPI_Type_struct(nblock_recv[p],lengths,displacements,
                          block_recv_types,&type_recv);
          MPI_Type_commit(&type_recv);

          // recv
          MPI_Request recv_req;
          MPI_Irecv(v.values_pt(),1,type_recv,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &recv_req);
          MPI_Type_free(&type_recv);
          for (int i = 0; i < nblock_recv[p]; i++)
           {
            MPI_Type_free(&block_recv_types[i]);
           }
          requests.push_back(recv_req);
         }

        // send
        if (nblock_send[p] > 0)
         {
          // create the datatypes vector
          MPI_Datatype block_send_types[nblock_send[p]];

          // and the displacements
          MPI_Aint displacements[nblock_send[p]];

          // and the lengths
          int lengths[nblock_send[p]];

          // all displacements are computed relative to s[0] values
          MPI_Aint displacements_base;
          MPI_Address(const_cast<double*>(s[0].values_pt()),
                      &displacements_base);

          // now build
          unsigned ptr = 0;
          for (unsigned b = 0; b < nblock; b++)
           {
            if (Nrows_to_recv_for_get_block(b,p) > 0)
             {
              MPI_Type_indexed(Nrows_to_recv_for_get_block(b,p),block_lengths,
                               Rows_to_recv_for_get_block(b,p),MPI_DOUBLE,
                               &block_send_types[ptr]);
              MPI_Type_commit(&block_send_types[ptr]);
              MPI_Address(const_cast<double*>(s[b].values_pt()),
                          &displacements[ptr]);
              displacements[ptr] -= displacements_base;
              lengths[ptr] = 1;
              ptr++;
             }
           }

          // build the final data type
          MPI_Datatype type_send;
          MPI_Type_struct(nblock_send[p],lengths,displacements,
                          block_send_types,&type_send);
          MPI_Type_commit(&type_send);

          // send
          MPI_Request send_req;
          MPI_Isend(const_cast<double*>(s[0].values_pt()),1,type_send,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &send_req);
          MPI_Type_free(&type_send);
          for (int i = 0; i < nblock_send[p]; i++)
           {
            MPI_Type_free(&block_send_types[i]);
           }
          requests.push_back(send_req);
         }
       }

      // communicate wih self
      else
       {
        double* v_values_pt = v.values_pt();
        for (unsigned b = 0; b < nblock; b++)
         {
          const double* w_values_pt = s[b].values_pt();
          for (unsigned i = 0; i < Nrows_to_send_for_get_block(b,p); i++)
           {
            v_values_pt[Rows_to_send_for_get_block(b,p)[i]] =
             w_values_pt[Rows_to_recv_for_get_block(b,p)[i]];

           }
         }
       }
     }

    // and then just wait
    unsigned c = requests.size();
    Vector<MPI_Status> stat(c);
    if (c)
     {
      MPI_Waitall(c,&requests[0],&stat[0]);
     }
    delete[] block_lengths;

#else
    // throw error
    std::ostringstream error_message;
    error_message << "The preconditioner is distributed and on more than one "
                  << "processor. MPI is required.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
#endif
   }
 }

 //============================================================================
 /// \short A helper function, takes the vector of block vectors, s, 
 /// and copies its entries into the naturally ordered vector, v. 
 /// This function assume that there are nblocks_precomputed block vectors 
 /// and they have the precomputed block ordering. If this is a subsidiary 
 /// block preconditioner only those entries in v that are associated with 
 /// its blocks are affected.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 return_block_vectors_with_precomputed_block_ordering(
   const Vector<DoubleVector >& s, DoubleVector& v) const
 {
  // the number of blocks
  unsigned nprecomputedblock = Block_to_block_map.size();

#ifdef PARANOID
  if(!Preconditioner_blocks_have_been_precomputed)
   {
    std::ostringstream error_message;
    error_message << "Precomputed blocks are not set. I cannot return block "
                  << "vectors with precomputed ordering.\n";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (!v.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  
  if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must match the "
                  << " specified master_distribution_pt(). \n"
                  << "i.e. Distribution_pt in the master preconditioner";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  
  if(s.size() > nprecomputedblock)
   {
    std::ostringstream error_message;
    error_message << "You have supplied " << s.size()
                  << " block vectors. I require " << nprecomputedblock
                  << " block vectors.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  for (unsigned b = 0; b < nprecomputedblock; b++)
   {
    if (!s[b].built())
     {
      std::ostringstream error_message;
      error_message << "The distribution of the block vector " << b
                    << " must be setup.";
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
    
    if (*(s[b].distribution_pt()) != *(Precomputed_block_distribution_pt[b]))
     {
      std::ostringstream error_message;
      error_message << "The distribution of the block vector " << b
                    << "must match the"
                    << " specified distribution at "
                    << "Precomputed_block_distribution_pt["
                    << b << "]";
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif

  for (unsigned b_i = 0; b_i < nprecomputedblock; b_i++) 
   {
    return_block_vector_with_precomputed_block_ordering(b_i,s[b_i],v);
   }

 }

 //============================================================================
 /// \short A helper function, takes the naturally ordered vector, v, 
 /// and extracts the n-th block vector, b. 
 /// Here n is the block number in the current preconditioner. 
 /// NOTE: The ordering of the vector b is the same as the 
 /// ordering of the block matrix from get_precomputed_block(...).
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 get_block_vector_with_precomputed_block_ordering(const unsigned& b, 
                                                  const DoubleVector& v, 
                                                  DoubleVector& w)
  const
 {
#ifdef PARANOID
  // the number of blocks
  unsigned n_blocks = Block_to_block_map.size();

  // paranoid check that block i is in this block preconditioner
  if (b >= n_blocks)
   {
    std::ostringstream error_message;
    error_message << "Requested block  vector " << b
                  << ", however this preconditioner has merged blocks "
                  << "= " << n_blocks << std::endl;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (!v.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must match the "
                  << " specified master_distribution_pt(). \n"
                  << "i.e. Distribution_pt in the master preconditioner";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  // Are precomputed blocks set?
  if(!Preconditioner_blocks_have_been_precomputed)
   {
    std::ostringstream warning_message;
    warning_message << "There are no precomputed blocks set."
                    << "This function may not give the ordering you desire.";
    OomphLibWarning(warning_message.str(),
                    OOMPH_CURRENT_FUNCTION,
                    OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // How many block vectors do we need to concatenate?
  unsigned nblocks_to_cat = Block_to_block_map[b].size();
  if(nblocks_to_cat == 1)
   // If there is only one block vector, we simply extract it.
   {
    this->get_block_vector_with_original_matrix_ordering(
                                                         Block_to_block_map[b][0],v,w);
   }
  else
   // We need to concatenate multiple block vectors.
   {
    // Get the DoubleVectors to be concatenated.
    Vector<DoubleVector> tmp_block_vector(nblocks_to_cat);
    for (unsigned b_i = 0; b_i < nblocks_to_cat; b_i++) 
     {
      this->get_block_vector_with_original_matrix_ordering(
                                                           Block_to_block_map[b][b_i], v,tmp_block_vector[b_i]);
     }
     
    Vector<DoubleVector*> tmp_block_vector_pt(nblocks_to_cat,0);
    for (unsigned b_i = 0; b_i < nblocks_to_cat; b_i++) 
     {
      tmp_block_vector_pt[b_i] = &tmp_block_vector[b_i];
     }
  
    // Build w with the correct distribution.
    w.build(Precomputed_block_distribution_pt[b],0);
  
    // Concatenate the vectors
    DoubleVectorHelpers::concatenate_without_communication
     (tmp_block_vector_pt,w);
  
    // No longer need the sub vectors. Calling clear on the Vector will invoke
    // the destructors in the DoubleVectors.
    tmp_block_vector.clear();
   }
 }

 //============================================================================
 /// \short A helper function, takes the naturally ordered vector, v, 
 /// and extracts the n-th block vector, b. 
 /// Here n is the block number in the current preconditioner. 
 /// NOTE: The ordering of the vector b is the same as the 
 /// ordering of the block matrix from get_block_from_original_matrix(...).
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 get_block_vector_with_original_matrix_ordering(const unsigned& b, 
                                                const DoubleVector& v, 
                                                DoubleVector& w)
  const
 {
#ifdef PARANOID
  // the number of blocks
  unsigned n_blocks = this->internal_nblock_types(true);

  // paranoid check that block i is in this block preconditioner
  if (b >= n_blocks)
   {
    std::ostringstream error_message;
    error_message << "Requested block  vector " << b
                  << ", however this preconditioner has internal_nblock_types(true) "
                  << "= " << internal_nblock_types(true) << std::endl;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (!v.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must match the "
                  << " specified master_distribution_pt(). \n"
                  << "i.e. Distribution_pt in the master preconditioner";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // rebuild the block vector
  w.build(Block_distribution_pt[b],0.0);

  // if + only one processor
  //    + more than one processor but matrix_pt is not distributed
  // then use the serial get_block method
  if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
      !this->distribution_pt()->distributed())
   {
    double* w_pt = w.values_pt();
    const double* v_pt = v.values_pt();
    unsigned n_row = w.nrow();
    for (unsigned i = 0; i < n_row; i++)
     {
      w_pt[i] = v_pt[this->Global_index[b][i]];
     }
   }
  // otherwise use mpi
  else
   {
#ifdef OOMPH_HAS_MPI

    // my rank
    unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

    // the number of processors
    unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

    // determine the maximum number of rows to be sent or recv
    unsigned max_n_send_or_recv = 0;
    for (unsigned p = 0; p < nproc; p++)
     {
      max_n_send_or_recv =
       std::max(max_n_send_or_recv,Nrows_to_send_for_get_block(b,p));
      max_n_send_or_recv =
       std::max(max_n_send_or_recv,Nrows_to_recv_for_get_block(b,p));
     }

    // create a vectors of 1s (the size of the nblock for the mpi indexed
    // data types
    int* block_lengths = new int[max_n_send_or_recv];
    for (unsigned i = 0; i < max_n_send_or_recv; i++)
     {
      block_lengths[i] = 1;
     }

    // perform the sends and receives
    Vector<MPI_Request> requests;
    for (unsigned p = 0; p < nproc; p++)
     {
      // send and recv with other processors
      if (p != my_rank)
       {
        if (Nrows_to_send_for_get_block(b,p) > 0)
         {
          // create the send datatype
          MPI_Datatype type_send;
          MPI_Type_indexed(Nrows_to_send_for_get_block(b,p),block_lengths,
                           Rows_to_send_for_get_block(b,p),MPI_DOUBLE,
                           &type_send);
          MPI_Type_commit(&type_send);

          // send
          MPI_Request send_req;
          MPI_Isend(const_cast<double*>(v.values_pt()),1,type_send,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &send_req);
          MPI_Type_free(&type_send);
          requests.push_back(send_req);
         }

        if (Nrows_to_recv_for_get_block(b,p) > 0)
         {
          // create the recv datatype
          MPI_Datatype type_recv;
          MPI_Type_indexed(Nrows_to_recv_for_get_block(b,p),block_lengths,
                           Rows_to_recv_for_get_block(b,p),MPI_DOUBLE,
                           &type_recv);
          MPI_Type_commit(&type_recv);

          // recv
          MPI_Request recv_req;
          MPI_Irecv(w.values_pt(),1,type_recv,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &recv_req);
          MPI_Type_free(&type_recv);
          requests.push_back(recv_req);
         }
       }

      // communicate wih self
      else
       {
        double* w_values_pt = w.values_pt();
        const double* v_values_pt = v.values_pt();
        for (unsigned i = 0; i < Nrows_to_send_for_get_block(b,p); i++)
         {
          w_values_pt[Rows_to_recv_for_get_block(b,p)[i]] =
           v_values_pt[Rows_to_send_for_get_block(b,p)[i]];
         }
       }
     }

    // and then just wait
    unsigned c = requests.size();
    Vector<MPI_Status> stat(c);
    if (c)
     {
      MPI_Waitall(c,&requests[0],&stat[0]);
     }
    delete[] block_lengths;

#else
    // throw error
    std::ostringstream error_message;
    error_message << "The preconditioner is distributed and on more than one "
                  << "processor. MPI is required.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
#endif
   }
 }

 //============================================================================
 /// \short Takes the naturally ordered vector, v and returns the n-th
 /// block vector, b. Here n is the block number in the current
 /// preconditioner. If blocks for this preconditioner has been precomputed
 /// then this function calls the function get_precomputed_block_vector(...).
 /// Otherwise it calls get_block_vector_from_original_matrix(...).
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 get_block_vector(const unsigned& b, const DoubleVector& v, DoubleVector& w)
  const
 {
  // If preconditioner blocks have been precomputed then we call 
  // get_block_vector_with_precomputed_block_ordering(...) to ensure that the
  // distribution of w is the same as the distribution of the precomputed
  // block.
  if(Preconditioner_blocks_have_been_precomputed)
   {
    get_block_vector_with_precomputed_block_ordering(b,v,w);
   }
  // Otherwise the n-th block matrix came from the original matrix, so we call
  // get_block_vector_with_original_matrix_ordering(...)
  else
   {
    get_block_vector_with_original_matrix_ordering(b,v,w);
   }
 }

 //============================================================================
 /// \short A helper function to return a block vector if the preconditioner
 /// blocks have been precomputed.
 /// Takes the n-th block ordered vector, b, and copies its entries
 /// to the appropriate entries in the naturally ordered vector, v.
 /// Here n is the block number in the current block preconditioner.
 /// If the preconditioner is a subsidiary block preconditioner
 /// the other entries in v  that are not associated with it
 /// are left alone.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 return_block_vector_with_precomputed_block_ordering(
                                                     const unsigned& b, const DoubleVector& w, DoubleVector& v) const
 {
#ifdef PARANOID
  // the number of blocks
  unsigned n_blocks = nblock_types_precomputed();

  // paranoid check that block i is in this block preconditioner
  if (b >= n_blocks)
   {
    std::ostringstream error_message;
    error_message << "Requested block  vector " << b
                  << ", however this preconditioner has nblock_types_precomputed() "
                  << "= " << nblock_types_precomputed() << std::endl;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (!v.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (*v.distribution_pt() != *this->master_distribution_pt())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must match the "
                  << " specified master_distribution_pt(). \n"
                  << "i.e. Distribution_pt in the master preconditioner";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (!w.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the block vector w must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  // Are precomputed blocks set?
  if(!Preconditioner_blocks_have_been_precomputed)
   {
    std::ostringstream warning_message;
    warning_message << "There are no precomputed blocks set."
                    << "This function may not give the ordering you desire.";
    OomphLibWarning(warning_message.str(),
                    OOMPH_CURRENT_FUNCTION,
                    OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // How many dof types does this block type represent?
  unsigned nblocks_to_split_into = Block_to_block_map[b].size();

  if(nblocks_to_split_into == 1)
   // If there is one block vector to return, we simply return the one block
   // vector
   {
    this->return_block_vector_with_original_matrix_ordering(
                                                            Block_to_block_map[b][0],w,v);
   }
  else
   // Otherwise this block vector correspondings to more than one dof types.
   // We need to split this.
   {
    // Temp vector to be returned to v.
    Vector<DoubleVector> tmp_block_vector(nblocks_to_split_into);
  
    // Fill it with the corresponding ditributions.
    for (unsigned b_i = 0; b_i < nblocks_to_split_into; b_i++) 
     {
      tmp_block_vector[b_i].build(this->block_distribution_pt
                                  (Block_to_block_map[b][b_i]));
     }
  
    Vector<DoubleVector*> tmp_block_vector_pt(nblocks_to_split_into,0);
    for (unsigned b_i = 0; b_i < nblocks_to_split_into; b_i++) 
     {
      tmp_block_vector_pt[b_i] = &tmp_block_vector[b_i];
     }
  
    // Perform the splitting of w into tmp_block_vector.
    DoubleVectorHelpers::split_without_communication(w,tmp_block_vector_pt);
  
    // return to v
    for (unsigned b_i = 0; b_i < nblocks_to_split_into; b_i++) 
     {
      this->return_block_vector_with_original_matrix_ordering(
                                                              Block_to_block_map[b][b_i], tmp_block_vector[b_i],v);
     }
   }
 }


 //============================================================================
 /// \short A helper function to return a block if no preconditioner blocks
 /// were precomputed.
 /// Takes the n-th block ordered vector, b,  and copies its entries
 /// to the appropriate entries in the naturally ordered vector, v.
 /// Here n is the block number in the current block preconditioner.
 /// If the preconditioner is a subsidiary block preconditioner
 /// the other entries in v  that are not associated with it
 /// are left alone.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 return_block_vector_with_original_matrix_ordering(const unsigned& b, 
                                                   const DoubleVector& w, 
                                                   DoubleVector& v)
  const
 {
#ifdef PARANOID
  // the number of blocks
  unsigned n_blocks = this->internal_nblock_types(true);

  // paranoid check that block i is in this block preconditioner
  if (b >= n_blocks)
   {
    std::ostringstream error_message;
    error_message << "Requested block  vector " << b
                  << ", however this preconditioner has internal_nblock_types(true) "
                  << "= " << internal_nblock_types(true) << std::endl;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (!v.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (*v.distribution_pt() != *this->master_distribution_pt())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must match the "
                  << " specified master_distribution_pt(). \n"
                  << "i.e. Distribution_pt in the master preconditioner";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (!w.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the block vector w must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (*w.distribution_pt() != *Block_distribution_pt[b])
   {
    std::ostringstream error_message;
    error_message << "The distribution of the block vector w must match the "
                  << " specified distribution at Block_distribution_pt[b]";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // if + only one processor
  //    + more than one processor but matrix_pt is not distributed
  // then use the serial get_block method
  if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
      !this->distribution_pt()->distributed())
   {

    // length of vector
    unsigned n_row = this->block_dimension(b);

    // copy back from the block vector to the naturally ordered vector
    double* v_pt = v.values_pt();
    const double* w_pt = w.values_pt();
    for (unsigned i = 0; i < n_row; i++)
     {
      v_pt[this->Global_index[b][i]] = w_pt[i];
     }
   }
  // otherwise use mpi
  else
   {
#ifdef OOMPH_HAS_MPI

    // my rank
    unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

    // the number of processors
    unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

    // determine the maximum number of rows to be sent or recv
    unsigned max_n_send_or_recv = 0;
    for (unsigned p = 0; p < nproc; p++)
     {
      max_n_send_or_recv =
       std::max(max_n_send_or_recv,Nrows_to_send_for_get_block(b,p));
      max_n_send_or_recv =
       std::max(max_n_send_or_recv,Nrows_to_recv_for_get_block(b,p));
     }

    // create a vectors of 1s (the size of the nblock for the mpi indexed
    // data types
    int* block_lengths = new int[max_n_send_or_recv];
    for (unsigned i = 0; i < max_n_send_or_recv; i++)
     {
      block_lengths[i] = 1;
     }

    // perform the sends and receives
    Vector<MPI_Request> requests;
    for (unsigned p = 0; p < nproc; p++)
     {
      // send and recv with other processors
      if (p != my_rank)
       {
        if (Nrows_to_recv_for_get_block(b,p) > 0)
         {
          // create the send datatype
          MPI_Datatype type_send;
          MPI_Type_indexed(Nrows_to_recv_for_get_block(b,p),block_lengths,
                           Rows_to_recv_for_get_block(b,p),MPI_DOUBLE,
                           &type_send);
          MPI_Type_commit(&type_send);

          // send
          MPI_Request send_req;
          MPI_Isend(const_cast<double*>(w.values_pt()),1,type_send,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &send_req);
          MPI_Type_free(&type_send);
          requests.push_back(send_req);
         }

        if (Nrows_to_send_for_get_block(b,p) > 0)
         {
          // create the recv datatype
          MPI_Datatype type_recv;
          MPI_Type_indexed(Nrows_to_send_for_get_block(b,p),block_lengths,
                           Rows_to_send_for_get_block(b,p),MPI_DOUBLE,
                           &type_recv);
          MPI_Type_commit(&type_recv);

          // recv
          MPI_Request recv_req;
          MPI_Irecv(v.values_pt(),1,type_recv,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &recv_req);
          MPI_Type_free(&type_recv);
          requests.push_back(recv_req);
         }
       }

      // communicate wih self
      else
       {
        const double* w_values_pt = w.values_pt();
        double* v_values_pt = v.values_pt();
        for (unsigned i = 0; i < Nrows_to_send_for_get_block(b,p); i++)
         {
          v_values_pt[Rows_to_send_for_get_block(b,p)[i]] =
           w_values_pt[Rows_to_recv_for_get_block(b,p)[i]];
         }
       }
     }

    // and then just wait
    unsigned c = requests.size();
    Vector<MPI_Status> stat(c);
    if (c)
     {
      MPI_Waitall(c,&requests[0],&stat[0]);
     }
    delete[] block_lengths;

#else
    // throw error
    std::ostringstream error_message;
    error_message << "The preconditioner is distributed and on more than one "
                  << "processor. MPI is required.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
#endif
   }
 }

 //============================================================================
 /// \short Takes the n-th block vector, b, and copies its entries
 /// to the appropriate entries in the naturally ordered vector, v, either
 /// by calling return_block_vector_with_precomputed_block_ordering(...)
 /// if the preconditioner blocks have been precomputed or
 /// return_block_vector_with_original_matrix_ordering(...) otherwise.
 /// Here n is the block number in the current block preconditioner.
 /// If the preconditioner is a subsidiary block preconditioner
 /// the other entries in v  that are not associated with it
 /// are left alone.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 return_block_vector(const unsigned& b, const DoubleVector& w, DoubleVector& v)
  const
 {
  // If the preconditioner blocks have been precomputed, then it is assumed
  // that the ordering of b is the same as the ordering of the
  // precomputed blocks. We invoke the correct function for returning the 
  // block vector.
  if(Preconditioner_blocks_have_been_precomputed)
   {
    return_block_vector_with_precomputed_block_ordering(b,w,v);
   }
  else
   {
    return_block_vector_with_original_matrix_ordering(b,w,v);
   }
 }

 //============================================================================
 /// \short Given the naturally ordered vector, v, return the vector rearranged
 /// in block order in w.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 get_block_ordered_preconditioner_vector(const DoubleVector& v, 
                                         DoubleVector& w)
  const
 {
#ifdef PARANOID
  if (!v.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (*v.distribution_pt() != *this->master_distribution_pt())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must match the "
                  << " specified master_distribution_pt(). \n"
                  << "i.e. Distribution_pt in the master preconditioner";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  //  Cleared and resized w for reordered vector
  w.build(this->preconditioner_matrix_distribution_pt(),0.0);

  // if + only one processor
  //    + more than one processor but matrix_pt is not distributed
  // then use the serial get_block method
  if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
      !this->distribution_pt()->distributed())
   {

    // number of blocks
    unsigned nblock = this->Internal_nblock_types;

    // copy to w
    unsigned block_offset = 0;
    double* w_pt = w.values_pt();
    const double* v_pt = v.values_pt();
    for (unsigned b = 0; b < nblock;b++)
     {
      unsigned block_nrow = this->block_dimension(b);
      for (unsigned i = 0; i < block_nrow; i++)
       {
        w_pt[block_offset+i] = v_pt[this->Global_index[b][i]];
       }
      block_offset += block_nrow;
     }
   }
  // otherwise use mpi
  else
   {
#ifdef OOMPH_HAS_MPI

    // my rank
    unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

    // the number of processors
    unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

    // determine the maximum number of rows to be sent or recv
    unsigned max_n_send_or_recv = 0;
    for (unsigned p = 0; p < nproc; p++)
     {
      max_n_send_or_recv =
       std::max(max_n_send_or_recv,Nrows_to_send_for_get_ordered[p]);
      max_n_send_or_recv =
       std::max(max_n_send_or_recv,Nrows_to_recv_for_get_ordered[p]);
     }

    // create a vectors of 1s (the size of the nblock for the mpi indexed
    // data types
    int* block_lengths = new int[max_n_send_or_recv];
    for (unsigned i = 0; i < max_n_send_or_recv; i++)
     {
      block_lengths[i] = 1;
     }

    // perform the sends and receives
    Vector<MPI_Request> requests;
    for (unsigned p = 0; p < nproc; p++)
     {
      // send and recv with other processors
      if (p != my_rank)
       {
        if (Nrows_to_send_for_get_ordered[p] > 0)
         {
          // create the send datatype
          MPI_Datatype type_send;
          MPI_Type_indexed(Nrows_to_send_for_get_ordered[p],block_lengths,
                           Rows_to_send_for_get_ordered[p],MPI_DOUBLE,
                           &type_send);
          MPI_Type_commit(&type_send);

          // send
          MPI_Request send_req;
          MPI_Isend(const_cast<double*>(v.values_pt()),1,type_send,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &send_req);
          MPI_Type_free(&type_send);
          requests.push_back(send_req);
         }

        if (Nrows_to_recv_for_get_ordered[p] > 0)
         {
          // create the recv datatype
          MPI_Datatype type_recv;
          MPI_Type_indexed(Nrows_to_recv_for_get_ordered[p],block_lengths,
                           Rows_to_recv_for_get_ordered[p],MPI_DOUBLE,
                           &type_recv);
          MPI_Type_commit(&type_recv);

          // recv
          MPI_Request recv_req;
          MPI_Irecv(w.values_pt(),1,type_recv,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &recv_req);
          MPI_Type_free(&type_recv);
          requests.push_back(recv_req);
         }
       }

      // communicate wih self
      else
       {
        double* w_values_pt = w.values_pt();
        const double* v_values_pt = v.values_pt();
        for (unsigned i = 0; i < Nrows_to_send_for_get_ordered[p]; i++)
         {
          w_values_pt[Rows_to_recv_for_get_ordered[p][i]] =
           v_values_pt[Rows_to_send_for_get_ordered[p][i]];
         }
       }
     }

    // and then just wait
    unsigned c = requests.size();
    Vector<MPI_Status> stat(c);
    if (c)
     {
      MPI_Waitall(c,&requests[0],&stat[0]);
     }
    delete[] block_lengths;

#else
    // throw error
    std::ostringstream error_message;
    error_message << "The preconditioner is distributed and on more than one "
                  << "processor. MPI is required.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
#endif
   }
 }

 //============================================================================
 /// \short Takes the naturally ordered vector, w, and reorders it in block
 /// order. Reordered vector is returned in v. Note: If the preconditioner is a
 /// subsidiary preconditioner then only the components of the vector
 /// associated with the blocks of the subsidiary preconditioner will be
 /// included.
 //============================================================================
 template<typename MATRIX> void BlockPreconditioner<MATRIX>::
 return_block_ordered_preconditioner_vector(const DoubleVector& w,
                                            DoubleVector& v) const
 {
#ifdef PARANOID
  if (!v.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (*v.distribution_pt() != *this->master_distribution_pt())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the global vector v must match the "
                  << " specified master_distribution_pt(). \n"
                  << "i.e. Distribution_pt in the master preconditioner";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (!w.built())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the block vector w must be setup.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
  if (*w.distribution_pt() != *this->preconditioner_matrix_distribution_pt())
   {
    std::ostringstream error_message;
    error_message << "The distribution of the block vector w must match the "
                  << " specified distribution at Distribution_pt[b]";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif


  // if + only one processor
  //    + more than one processor but matrix_pt is not distributed
  // then use the serial get_block method
  if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
      !this->distribution_pt()->distributed())
   {
    // number of blocks
    unsigned nblock = this->Internal_nblock_types;

    // copy to w
    unsigned block_offset = 0;
    const double* w_pt = w.values_pt();
    double* v_pt = v.values_pt();
    for (unsigned b = 0; b < nblock;b++)
     {
      unsigned block_nrow = this->block_dimension(b);
      for (unsigned i = 0; i < block_nrow; i++)
       {
        v_pt[this->Global_index[b][i]] = w_pt[block_offset+i];
       }
      block_offset += block_nrow;
     }
   }
  // otherwise use mpi
  else
   {
#ifdef OOMPH_HAS_MPI

    // my rank
    unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

    // the number of processors
    unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

    // determine the maximum number of rows to be sent or recv
    unsigned max_n_send_or_recv = 0;
    for (unsigned p = 0; p < nproc; p++)
     {
      max_n_send_or_recv =
       std::max(max_n_send_or_recv,Nrows_to_send_for_get_ordered[p]);
      max_n_send_or_recv =
       std::max(max_n_send_or_recv,Nrows_to_recv_for_get_ordered[p]);
     }

    // create a vectors of 1s (the size of the nblock for the mpi indexed
    // data types
    int* block_lengths = new int[max_n_send_or_recv];
    for (unsigned i = 0; i < max_n_send_or_recv; i++)
     {
      block_lengths[i] = 1;
     }

    // perform the sends and receives
    Vector<MPI_Request> requests;
    for (unsigned p = 0; p < nproc; p++)
     {
      // send and recv with other processors
      if (p != my_rank)
       {
        if (Nrows_to_recv_for_get_ordered[p] > 0)
         {
          // create the send datatype
          MPI_Datatype type_send;
          MPI_Type_indexed(Nrows_to_recv_for_get_ordered[p],block_lengths,
                           Rows_to_recv_for_get_ordered[p],MPI_DOUBLE,
                           &type_send);
          MPI_Type_commit(&type_send);

          // send
          MPI_Request send_req;
          MPI_Isend(const_cast<double*>(w.values_pt()),1,type_send,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &send_req);
          MPI_Type_free(&type_send);
          requests.push_back(send_req);
         }

        if (Nrows_to_send_for_get_ordered[p] > 0)
         {
          // create the recv datatype
          MPI_Datatype type_recv;
          MPI_Type_indexed(Nrows_to_send_for_get_ordered[p],block_lengths,
                           Rows_to_send_for_get_ordered[p],MPI_DOUBLE,
                           &type_recv);
          MPI_Type_commit(&type_recv);

          // recv
          MPI_Request recv_req;
          MPI_Irecv(v.values_pt(),1,type_recv,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),&recv_req);
          MPI_Type_free(&type_recv);
          requests.push_back(recv_req);
         }
       }

      // communicate wih self
      else
       {
        const double* w_values_pt = w.values_pt();
        double* v_values_pt = v.values_pt();
        for (unsigned i = 0; i < Nrows_to_send_for_get_ordered[p]; i++)
         {
          v_values_pt[Rows_to_send_for_get_ordered[p][i]] =
           w_values_pt[Rows_to_recv_for_get_ordered[p][i]];
         }
       }
     }

    // and then just wait
    unsigned c = requests.size();
    Vector<MPI_Status> stat(c);
    if (c)
     {
      MPI_Waitall(c,&requests[0],&stat[0]);
     }
    delete[] block_lengths;

#else
    // throw error
    std::ostringstream error_message;
    error_message << "The preconditioner is distributed and on more than one "
                  << "processor. MPI is required.";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
#endif
   } // else use mpi
 } // function return_block_ordered_preconditioner_vector

//=============================================================================
/// \short Gets block (i,j) from the original matrix and returns it in
/// block_matrix_pt (Specialisation for CRDoubleMatrix)
//=============================================================================
 template<> 
 void BlockPreconditioner<CRDoubleMatrix>:: 
 get_block_from_original_matrix(const unsigned& block_i, const unsigned& block_j, 
                                CRDoubleMatrix& output_block) const
 {

#ifdef PARANOID
  // the number of blocks
  unsigned n_blocks = this->internal_nblock_types(true);

  // paranoid check that block i is in this block preconditioner
  if (block_i >= n_blocks || block_j >= n_blocks)
   {
    std::ostringstream error_message;
    error_message << "Requested block (" << block_i << "," << block_j   
                  << "), however this preconditioner has internal_nblock_types(true) "
                  << "= " << internal_nblock_types(true) << std::endl;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  // Check that the matrix is the same as that of the master
  if(is_subsidiary_block_preconditioner())
   {
    if(master_block_preconditioner_pt()->matrix_pt() != matrix_pt())
     {
        std::string err = "Master and subs should have same matrix.";
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif

    // Cast the pointer
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());

  // if + only one processor
  //    + more than one processor but matrix_pt is not distributed
  // then use the serial get_block method
    if (cr_matrix_pt->distribution_pt()->communicator_pt()->nproc() == 1 ||
	!cr_matrix_pt->distribution_pt()->distributed())
   {
    // pointers for the jacobian matrix is compressed row sparse format 
    int* j_row_start;
    int* j_column_index;
    double* j_value;
    
    // sets pointers to jacobian matrix
	j_row_start = cr_matrix_pt->row_start();
	j_column_index = cr_matrix_pt->column_index();
	j_value = cr_matrix_pt->value();
    
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
      
      
    // Fill in the compressed row matrix ??ds Note: I kept the calls to
    // build as close as I could to before (had to replace new(dist) with
    // .build(dist) ).
    output_block.build(Block_distribution_pt[block_i]);
    output_block.build_without_copy(block_ncol,block_nnz,
                                    temp_value,temp_column_index,
                                    temp_row_start);
 
#ifdef PARANOID
    // checks to see if block matrix has been set up correctly 
    //   block_matrix_test(matrix_pt,block_i,block_j,block_pt);
    if (Run_block_matrix_test)
     {
      // checks to see if block matrix has been set up correctly 
      block_matrix_test(block_i, block_j, &output_block);
     }
#endif
   }


  // otherwise we are dealing with a distributed matrix
  else
   {
#ifdef OOMPH_HAS_MPI
    // number of processors
    unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

    // my rank
    unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

    // sets pointers to jacobian matrix
	int* j_row_start = cr_matrix_pt->row_start();
	int* j_column_index = cr_matrix_pt->column_index();
	double* j_value = cr_matrix_pt->value();

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
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &req);
          send_req.push_back(req);
         }

        // recv
        if (nrow_recv)
         {
          MPI_Request req;
          MPI_Irecv(nnz_recv[p],nrow_recv,MPI_INT,p,0,
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &req);
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
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &req);
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
          MPI_Address(values_recv,&displacements[0]);
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
                    this->distribution_pt()->communicator_pt()->mpi_comm(),
                    &req);
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

    // Fill in the compressed row matrix
    output_block.build(Block_distribution_pt[block_i]);
    output_block.build_without_copy(this->block_dimension(block_j),
                                    local_block_nnz,
                                    values_recv,
                                    column_index_recv,
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
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
#endif 
   }
 }


//=============================================================================
/// \short Gets block (i,j) from Replacement_block_pt and returns it in
/// block_matrix_pt.
//=============================================================================
 template<> 
 void BlockPreconditioner<CRDoubleMatrix>:: 
 get_coarsened_block(const unsigned& block_i, const unsigned& block_j, 
                       CRDoubleMatrix& output_block) const
 {
#ifdef PARANOID
  // the number of blocks
  unsigned nblocks = nblock_types_precomputed();
  
  // paranoid check that block i is in this block preconditioner
  if (block_i >= nblocks || block_j >= nblocks)
   {
    std::ostringstream error_message;
    error_message << "Requested block (" << block_i << "," << block_j   
                  << "), however this preconditioner has nblock_types_precomputed() "
                  << "= " << nblocks << std::endl;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  if(!Preconditioner_blocks_have_been_precomputed)
  {
    std::ostringstream error_message;
    error_message << "There are no precomputed blocks. Please call "
                  << "set_replacement_block(...)  ";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Create the dense matrix required for the merge.
  // How many block rows and columns?
  const unsigned nblock_in_row = Block_to_block_map[block_i].size();
  const unsigned nblock_in_col = Block_to_block_map[block_j].size();
  
  if((nblock_in_row == 1) && (nblock_in_col == 1))
   {
     
    // Do not need to invoke concatenate function.
    unsigned prec_block_i = Block_to_block_map[block_i][0];
    unsigned prec_block_j = Block_to_block_map[block_j][0];

    // Cache the pointer to the precomputed block.
    CRDoubleMatrix* precom_block_pt = 0;
    if(Preconditioner_blocks_have_been_precomputed)
    {
      precom_block_pt = Replacement_block_pt(prec_block_i,prec_block_j);
    }
    else
    {
      precom_block_pt = new CRDoubleMatrix;
      get_block_from_original_matrix(prec_block_i,prec_block_j,*precom_block_pt);
    }

    // Create a new matrix with the precomputed block
    // (prec_block_i,prec_block_j)

    // Temp storage.
    Vector<double> tmp_values;
    Vector<int> tmp_column_indices;
    Vector<int> tmp_row_start;
    
    // The precomputed block nrow and nnz
    unsigned precom_nrow = precom_block_pt->nrow();
    unsigned long precom_nnz = precom_block_pt->nnz();
    
    // Reserve space.
    tmp_values.reserve(precom_nnz);
    tmp_column_indices.reserve(precom_nnz);
    tmp_row_start.reserve(precom_nrow + 1);

    // The data to copy over.
    double* precom_values = precom_block_pt->value();
    int* precom_column_indices = precom_block_pt->column_index();
    int* precom_row_start = precom_block_pt->row_start();
    
    // Copy the values and column indices.
    for (unsigned i = 0; i < precom_nnz; i++) 
     {
      tmp_values.push_back(precom_values[i]);
      tmp_column_indices.push_back(precom_column_indices[i]);
     }

    // Copy the row start
    for (unsigned i = 0; i < precom_nrow + 1; i++) 
     {
      tmp_row_start.push_back(precom_row_start[i]);
     }

    unsigned precom_ncol
      = precom_block_pt->ncol();

    output_block.build(precom_block_pt->distribution_pt(),
                       precom_ncol,
                       tmp_values,
                       tmp_column_indices,
                       tmp_row_start);
    if(!Preconditioner_blocks_have_been_precomputed)
    {
      delete precom_block_pt;
    }
   }
  else
   {
  DenseMatrix<CRDoubleMatrix*> tmp_block_pt(nblock_in_row,nblock_in_col,0);
  Vector<LinearAlgebraDistribution*> tmp_row_distribution_pt(nblock_in_row,0);
  Vector<LinearAlgebraDistribution*> tmp_col_distribution_pt(nblock_in_col,0);

  // Fill in the corresponding matrices.
  for (unsigned block_row_i = 0; block_row_i < nblock_in_row; block_row_i++) 
   {
    unsigned prec_block_i = Block_to_block_map[block_i][block_row_i];

    for (unsigned block_col_i = 0; block_col_i < nblock_in_col; block_col_i++) 
     {
      unsigned prec_block_j = Block_to_block_map[block_j][block_col_i];
    
      tmp_block_pt(block_row_i,block_col_i) = 0;
      if(Preconditioner_blocks_have_been_precomputed)
      { 
        tmp_block_pt(block_row_i,block_col_i)
          = Replacement_block_pt(prec_block_i, prec_block_j);
      }
      else
      {
        tmp_block_pt(block_row_i,block_col_i) = new CRDoubleMatrix;
        get_block_from_original_matrix(prec_block_i,prec_block_j,
                                       *tmp_block_pt(block_row_i,block_col_i));
      }
     } // for
   } // for

  // Fill in the row distributions, use the first block column.
  for (unsigned block_row_i = 0; block_row_i < nblock_in_row; block_row_i++) 
   {
    tmp_row_distribution_pt[block_row_i] 
      = tmp_block_pt(block_row_i,0)->distribution_pt();
   }

  // Fill in the col distributions, use the first block row.
  // This is a bit more tricky, we need the distributions of the block
  // rows that these block columns correspond to.
  for (unsigned block_col_i = 0; block_col_i < nblock_in_col; block_col_i++) 
   {
    unsigned prec_row_block_i = Block_to_block_map[block_j][block_col_i];

    tmp_col_distribution_pt[block_col_i] 
      = Replacement_block_pt(prec_row_block_i,0)->distribution_pt();
   }

  output_block.build(Precomputed_block_distribution_pt[block_i]);

  // Concatenate the matrix.
  // For now, we use concatenate_without_communication(...) since none of the
  // current preconditioners require the block matrix to be in a particular
  // arrangement. We could use concatenate(...) which requires communication.
  CRDoubleMatrixHelpers::concatenate_without_communication(
    tmp_row_distribution_pt,tmp_col_distribution_pt, tmp_block_pt, output_block);
   }
 }




//=============================================================================
/// \short test function to check that every element in the block matrix
/// (block_i,block_j) matches the corresponding element in the original matrix
//=============================================================================
  template<typename MATRIX> void BlockPreconditioner<MATRIX>::
  block_matrix_test(const unsigned& block_i, const unsigned& block_j,
		    const MATRIX* block_matrix_pt) const
 {

  // boolean flag to indicate whether test is passed
  bool check = true;
  
  // number of rows in matrix
    unsigned n_row = matrix_pt()->nrow();
  
  // number of columns in matrix
    unsigned n_col = matrix_pt()->ncol();
  
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
		    if ( matrix_pt()->operator()(i,j) !=
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
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
 }


 template class BlockPreconditioner<CRDoubleMatrix>;

} // Namespace: oomph

