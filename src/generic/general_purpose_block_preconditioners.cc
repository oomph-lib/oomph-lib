//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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

#include "oomph_definitions.h"
#include "preconditioner.h"
#include "block_preconditioner.h"
#include "matrices.h"
#include "general_purpose_block_preconditioners.h"

namespace oomph
{

 //============================================================================
 /// setup for the block diagonal preconditioner
 //============================================================================
 template<typename MATRIX> 
 void BlockDiagonalPreconditioner<MATRIX>::setup()
 {
  // clean the memory
  this->clean_up_memory();

  // Note: Subsidiary block preconditioners don't really use their
  // mesh pointers (since the lookup schemes inferred from them are
  // set up by their masters). Generally we insist that (for uniformity)
  // a mesh should be set, but sometimes subsidiary preconditioners are 
  // set up on the fly and we can't sensibly set meshes, so we're forgiving...)
  if (this->is_master_block_preconditioner())
   {
#ifdef PARANOID
    if( this->gp_nmesh() == 0)
     {
      std::ostringstream err_msg;
      err_msg << "There are no meshes set.\n"
              << "Did you remember to call add_mesh(...)?";
      throw OomphLibError(err_msg.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    // Set all meshes if this is master block preconditioner
    this->gp_preconditioner_set_all_meshes();
   }

  // Set up the block look up schemes
  this->gp_preconditioner_block_setup();

  // number of types of blocks.
  unsigned nblock_types = this->nblock_types();

  // Create any subsidiary preconditioners needed
  this->fill_in_subsidiary_preconditioners(nblock_types);

  // If using two level parallelisation then we need to use a
  // PrecondtionerArray which requires very different setup. ??ds possibly
  // it should have it's own class?
  if(Use_two_level_parallelisation)
   {
    // Get the blocks. We have to use new because you can't have containers
    // of matrices because the copy constructor is buggy (so we create a
    // container of pointers instead). ??ds
    Vector<CRDoubleMatrix*> block_diagonal_matrix_pt(nblock_types, 0);
    for(unsigned i=0; i<nblock_types; i++)
     {
      block_diagonal_matrix_pt[i] = new CRDoubleMatrix;
      this->get_block(i, get_other_diag_ds(i, nblock_types),
                      *block_diagonal_matrix_pt[i]);
     }

    // Construct the preconditioner array
    Preconditioner_array_pt = new PreconditionerArray;

    Preconditioner_array_pt->
     setup_preconditioners(block_diagonal_matrix_pt,
                           this->Subsidiary_preconditioner_pt,
                           this->comm_pt());

    // and delete the blocks
    for(unsigned i=0; i<nblock_types; i++)
     {
      delete block_diagonal_matrix_pt[i];
      block_diagonal_matrix_pt[i] = 0;
     }

    // Preconditioner array is weird, it calls delete on all the
    // preconditioners you give it and requires new ones each time!
    this->Subsidiary_preconditioner_pt.clear();
   }


  // Otherwise just set up each block's preconditioner in order
  else
   {
    for(unsigned i=0; i<nblock_types; i++)
     {
      // Get the block
      unsigned j = get_other_diag_ds(i, nblock_types);
      CRDoubleMatrix block = this->get_block(i, j);

      // Set up preconditioner (i.e. approximately solve the block + store solution)
      double superlusetup_start = TimingHelpers::timer();
      this->Subsidiary_preconditioner_pt[i]->setup(&block);
      double superlusetup_end = TimingHelpers::timer();
      oomph_info << "Took " << superlusetup_end - superlusetup_start
                 << "s to setup."<< std::endl;
     }
   }

 }
 
 //=============================================================================
 /// Preconditioner solve for the block diagonal preconditioner
 //=============================================================================
 template<typename MATRIX> 
 void BlockDiagonalPreconditioner<MATRIX>::
 preconditioner_solve(const DoubleVector& r, DoubleVector& z)
 {
  // Cache umber of block types
  unsigned n_block = this->nblock_types();

  // Get the right hand side vector (b in Ax = b) in block form.
  Vector<DoubleVector> block_r;
  this->get_block_vectors(r, block_r);

  // Make sure z vector is built
  if (!z.built()) { z.build(this->distribution_pt(),0.0); }

  // Vector of vectors for storage of solution in block form.
  Vector<DoubleVector> block_z(n_block);

  if (Use_two_level_parallelisation)
   {
    Preconditioner_array_pt->solve_preconditioners(block_r, block_z);
   }
  else
   {
    // solve each diagonal block
    for (unsigned i = 0; i < n_block; i++)
     {
      double t_start=0.0;
      if (Doc_time_during_preconditioner_solve)
       {
        t_start=TimingHelpers::timer();
       }
      this->Subsidiary_preconditioner_pt[i]->preconditioner_solve(block_r[i],
                                                                   block_z[i]);
      if (Doc_time_during_preconditioner_solve)
       {
        oomph_info << "Time for application of " << i 
                   << "-th block preconditioner: " 
                   << TimingHelpers::timer()-t_start 
                   << std::endl;
       }
     }
   }

  // copy solution in block vectors block_z back to z
  this->return_block_vectors(block_z,z);
 }


 //============================================================================
 /// setup for the block triangular preconditioner
 //============================================================================
 template<typename MATRIX> 
 void BlockTriangularPreconditioner<MATRIX>::
 setup()
 {
  // clean the memory
  this->clean_up_memory();

  // Subsidiary preconditioners don't really need the meshes
 if (this->is_master_block_preconditioner())
  {
#ifdef PARANOID
   if( this->gp_nmesh() == 0)
     {
      std::ostringstream err_msg;
      err_msg << "There are no meshes set.\n"
              << "Did you remember to call add_mesh(...)?";
      throw OomphLibError(err_msg.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    // Set all meshes if this is master block preconditioner
    this->gp_preconditioner_set_all_meshes();
  }

  // Set up the block look up schemes
  this->gp_preconditioner_block_setup();
   
  // number of block types
  unsigned nblock_types = this->nblock_types();

  // storage for the off diagonal matrix vector products
  Off_diagonal_matrix_vector_products.resize(nblock_types,nblock_types,0);

  // Fill in any null subsidiary preconditioners
  this->fill_in_subsidiary_preconditioners(nblock_types);

  // build the preconditioners and matrix vector products
  for(unsigned i = 0; i < nblock_types; i++)
   {
    // Get the block and set up the preconditioner.
    {
     CRDoubleMatrix block_matrix = this->get_block(i,i);
     this->Subsidiary_preconditioner_pt[i]
      ->setup(&block_matrix);
    }
     
    // next setup the off diagonal mat vec operators
    unsigned l, u;
    if(Upper_triangular)
     {
      l = i+1;
      u = nblock_types;
     }
    else
     {
      l = 0;
      u = i;
     }

    for(unsigned j = l; j < u; j++)
     {
      // Get the block
      CRDoubleMatrix block_matrix = this->get_block(i,j);

      // Copy the block into a "multiplier" class. If trilinos is being
      // used this should also be faster than oomph-lib's multiphys.
      Off_diagonal_matrix_vector_products(i,j) = new MatrixVectorProduct();

      this->setup_matrix_vector_product(
              Off_diagonal_matrix_vector_products(i,j),
              &block_matrix,j);
     }
   }
 }
 
 //=============================================================================
 /// Preconditioner solve for the block triangular preconditioner
 //=============================================================================
 template<typename MATRIX> void BlockTriangularPreconditioner<MATRIX>::
 preconditioner_solve(const DoubleVector& r, DoubleVector& z)
 {
  // Cache number of block types
  unsigned n_block = this->nblock_types();
  
  //
  int start, end, step;
  if(Upper_triangular)
   {
    start = n_block-1;
    end = -1;
    step = -1;
   }
  else
   {
    start = 0;
    end = n_block;
    step = 1;
   }

  // vector of vectors for each section of residual vector
  Vector<DoubleVector> block_r;
  
  // rearrange the vector r into the vector of block vectors block_r
  this->get_block_vectors(r,block_r);

  // vector of vectors for the solution block vectors
  Vector<DoubleVector> block_z(n_block);

  //
  for (int i = start; i != end; i+=step)
   {
    // ??ds ugly, fix this?
    if(dynamic_cast<BlockPreconditioner<CRDoubleMatrix>*>
       (this->Subsidiary_preconditioner_pt[i]) == 0)
     {
      // solve on the block
      this->Subsidiary_preconditioner_pt[i]->
       preconditioner_solve(block_r[i], block_z[i]);
     }
    else
     {
      // This is pretty inefficient but there is no other way to do this with
      // block sub preconditioners as far as I can tell because they demand
      // to have the full r and z vectors...

      DoubleVector block_z_with_size_of_full_z(r.distribution_pt()),
       r_updated(r.distribution_pt());

      // Construct the new r vector (with the update given by back subs
      // below).
      this->return_block_vectors(block_r, r_updated);

      // Solve (blocking is handled inside the block preconditioner).
      this->Subsidiary_preconditioner_pt[i]->
       preconditioner_solve(r_updated, block_z_with_size_of_full_z);

      // Extract this block's z vector because we need it for the next step
      // (and throw away block_z_with_size_of_full_z).
      this->get_block_vector(i, block_z_with_size_of_full_z,  block_z[i]);
     }

    // substitute
    for (int j = i + step; j !=end; j+=step)
     {
      DoubleVector temp;
      Off_diagonal_matrix_vector_products(j,i)->multiply(block_z[i],temp);
      block_r[j] -= temp;
     }
   }

  // copy solution in block vectors block_r back to z
  this->return_block_vectors(block_z,z);
 }




 //=============================================================================
 /// Setup for the exact block preconditioner
 //=============================================================================
 template<typename MATRIX> 
 void ExactBlockPreconditioner<MATRIX>::setup()
 {

  // If this is a master block preconditioner, 
  // we give the mesh pointers to the BlockPreconditioner base class.
  // Subsidiary preconditioners don't need meshes...
  if(this->is_master_block_preconditioner())
  {
#ifdef PARANOID
   if( this->gp_nmesh() == 0)
    {
      std::ostringstream err_msg;
      err_msg << "There are no meshes set.\n"
              << "Did you remember to call add_mesh(...)?";
      throw OomphLibError(err_msg.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Set all meshes if this is master block preconditioner
    this->gp_preconditioner_set_all_meshes();
  }

  // Set up the block look up schemes
  this->gp_preconditioner_block_setup();

  // get the number of DOF types
  unsigned nblock_types = this->nblock_types();

  // Build the preconditioner matrix
  VectorMatrix<BlockSelector> required_blocks(nblock_types,nblock_types);

  // boolean indicating if we want the block, stored for readability.
  const bool want_block = true;
  for (unsigned b_i = 0; b_i < nblock_types; b_i++) 
  {
    for (unsigned b_j = 0; b_j < nblock_types; b_j++) 
    {
      required_blocks[b_i][b_j].select_block(b_i,b_j,want_block);
    }
  }

  // Get the concatenated blocks
  CRDoubleMatrix exact_block_matrix 
    = this->get_concatenated_block(required_blocks);
  
  // Setup the preconditioner.
  this->fill_in_subsidiary_preconditioners(1); // Only need one preconditioner
  preconditioner_pt()->setup(&exact_block_matrix);
 }
 
 //=============================================================================
 /// Preconditioner solve for the exact block preconditioner
 //=============================================================================
 template<typename MATRIX> 
 void ExactBlockPreconditioner<MATRIX>::
 preconditioner_solve(const DoubleVector& r, DoubleVector& z)
 {
  // get  the block ordered components of the r vector for this preconditioner
  DoubleVector block_order_r;
  this->get_block_ordered_preconditioner_vector(r,block_order_r);

  // vector for solution
  DoubleVector block_order_z;

  // apply the preconditioner
  preconditioner_pt()->preconditioner_solve(block_order_r,block_order_z);

  // copy solution back to z vector
  this->return_block_ordered_preconditioner_vector(block_order_z,z);
 }

 template class BlockDiagonalPreconditioner<CRDoubleMatrix>;
 template class BlockTriangularPreconditioner<CRDoubleMatrix>;
 template class ExactBlockPreconditioner<CRDoubleMatrix>;
 template class BlockAntiDiagonalPreconditioner<CRDoubleMatrix>;
 template class DummyBlockPreconditioner<CRDoubleMatrix>;

}
