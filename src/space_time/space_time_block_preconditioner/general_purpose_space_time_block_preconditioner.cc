// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Oomph-lib headers
#include "generic/SuperLU_preconditioner.h"

// Header file
#include "general_purpose_space_time_block_preconditioner.h"

/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////

namespace oomph
{
  //============================================================================
  /// Set up the exact preconditioner
  //============================================================================
  template<typename MATRIX>
  void ExactDGPBlockPreconditioner<MATRIX>::setup()
  {
    // Clean the memory
    this->clean_up_memory();

    // Subsidiary preconditioners don't really need the meshes
    if (this->is_master_block_preconditioner())
    {
#ifdef PARANOID
      if (this->gp_nmesh() == 0)
      {
        std::ostringstream err_msg;
        err_msg << "There are no meshes set.\n"
                << "Did you remember to call add_mesh(...)?";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Set all meshes if this is master block preconditioner
      this->gp_preconditioner_set_all_meshes();
    }

    // Initialise the memory usage variable
    Memory_usage_in_bytes = 0.0;

    // If we're meant to build silently
    if (this->Silent_preconditioner_setup == true)
    {
      // Store the output stream pointer
      this->Stream_pt = oomph_info.stream_pt();

      // Now set the oomph_info stream pointer to the null stream to
      // disable all possible output
      oomph_info.stream_pt() = &oomph_nullstream;
    }

    // Set up the block look up schemes
    this->gp_preconditioner_block_setup();

    // Number of block types
    unsigned n_block_types = this->nblock_types();

    // Fill in any null subsidiary preconditioners
    this->fill_in_subsidiary_preconditioners(n_block_types);

    // The total time for setting up the subsidiary preconditioners
    double t_subsidiary_setup_total = 0.0;

    // Stores the blocks we want to extract
    VectorMatrix<BlockSelector> required_blocks(n_block_types, n_block_types);

    // Boolean indicating if we want the block, stored for readability.
    const bool want_block = true;

    // Loop over the block rows
    for (unsigned i = 0; i < n_block_types; i++)
    {
      // Loop over the block columns
      for (unsigned j = 0; j < n_block_types; j++)
      {
        // Indicate that we want this block
        required_blocks[i][j].select_block(i, j, want_block);
      }
    } // for (unsigned i=0;i<n_block_types;i++)

    // Get the start time
    double t_extract_start = TimingHelpers::timer();

    // Get the concatenated blocks
    CRDoubleMatrix full_matrix = this->get_concatenated_block(required_blocks);

    // The total time for extracting all the blocks from the "global" matrix
    double t_extraction_total = (TimingHelpers::timer() - t_extract_start);

    // Get the start time
    double t_subsidiary_setup_start = TimingHelpers::timer();

    // It's a block preconditioner so pass it a pointer to the global matrix!
    this->Subsidiary_preconditioner_pt[0]->setup(&full_matrix);

    // Update the timing total
    t_subsidiary_setup_total +=
      (TimingHelpers::timer() - t_subsidiary_setup_start);

    // Remember that the preconditioner has been set up
    Preconditioner_has_been_setup = true;

    // If we're meant to build silently, reassign the oomph stream pointer
    if (this->Silent_preconditioner_setup == true)
    {
      // Store the output stream pointer
      oomph_info.stream_pt() = this->Stream_pt;

      // Reset our own stream pointer
      this->Stream_pt = 0;
    }

    // Tell the user
    oomph_info << "Total block extraction time [sec]: " << t_extraction_total
               << "\nTotal subsidiary preconditioner setup time [sec]: "
               << t_subsidiary_setup_total << std::endl;

    // Do we need to doc. the memory statistics?
    // NOTE: We're going to assume that either:
    //    (1) GMRES is used as a subsidiary block preconditioner (preconditioned
    //        by the Navier-Stokes subsidiary block preconditioner), or;
    //    (2) SuperLU is used as the subsidiary preconditioner.
    if (Compute_memory_statistics)
    {
      // How many rows are there in the global Jacobian?
      unsigned n_row = this->matrix_pt()->nrow();

      // How many nonzeros are there in the global Jacobian?
      unsigned n_nnz = this->matrix_pt()->nnz();

      // Add in the subsidiary preconditioners contribution
      double total_memory_usage_for_setup_phase =
        dynamic_cast<SuperLUPreconditioner*>(
          this->Subsidiary_preconditioner_pt[0])
          ->get_total_memory_needed_for_superlu();

      // Add in the global Jacobian contribution
      total_memory_usage_for_setup_phase +=
        ((2 * ((n_row + 1) * sizeof(int))) +
         (n_nnz * (sizeof(int) + sizeof(double))));

      // How much memory have we used in the subsidiary preconditioners?
      oomph_info << "\nTotal amount of memory being used after setup (MB): "
                 << total_memory_usage_for_setup_phase / 1.0e+06 << "\n"
                 << std::endl;
    } // if (Compute_memory_statistics)
  } // End of setup


  //=============================================================================
  /// Preconditioner solve for the exact preconditioner
  //=============================================================================
  template<typename MATRIX>
  void ExactDGPBlockPreconditioner<MATRIX>::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
    // Vector of vectors for each section of residual vector
    Vector<DoubleVector> block_r;

    // Rearrange the vector r into the vector of block vectors block_r
    this->get_block_vectors(r, block_r);

    // Allocate space for the properly rearranged RHS vector
    DoubleVector rhs_reordered;

    // Concatenate the DoubleVectors together
    DoubleVectorHelpers::concatenate(block_r, rhs_reordered);

    // Allocate space for the rearranged solution vector
    DoubleVector z_reordered;

    // Solve the whole system exactly
    this->Subsidiary_preconditioner_pt[0]->preconditioner_solve(rhs_reordered,
                                                                z_reordered);

    // Vector of vectors for each section of solution vector (copy the RHS
    // block vector to get the sizing and distributions right)
    Vector<DoubleVector> block_z(block_r);

    // Split the solution vector into blocks
    DoubleVectorHelpers::split(z_reordered, block_z);

    // Copy the solution from the block vector block_z back into z
    this->return_block_vectors(block_z, z);
  } // End of preconditioner_solve

  //============================================================================
  /// Set up the block triangular preconditioner
  //============================================================================
  template<typename MATRIX>
  void BandedBlockTriangularPreconditioner<MATRIX>::setup()
  {
    // Clean the memory
    this->clean_up_memory();

    // Subsidiary preconditioners don't really need the meshes
    if (this->is_master_block_preconditioner())
    {
#ifdef PARANOID
      if (this->gp_nmesh() == 0)
      {
        std::ostringstream err_msg;
        err_msg << "There are no meshes set.\n"
                << "Did you remember to call add_mesh(...)?";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Set all meshes if this is master block preconditioner
      this->gp_preconditioner_set_all_meshes();
    }

    // Initialise the memory usage variable
    Memory_usage_in_bytes = 0.0;

    // If we're meant to build silently
    if (this->Silent_preconditioner_setup == true)
    {
      // Store the output stream pointer
      this->Stream_pt = oomph_info.stream_pt();

      // Now set the oomph_info stream pointer to the null stream to
      // disable all possible output
      oomph_info.stream_pt() = &oomph_nullstream;
    }

    // Set up the block look up schemes
    this->gp_preconditioner_block_setup();

    // Number of block types
    unsigned n_block_types = this->nblock_types();

    // Storage for the off diagonal matrix vector products
    Off_diagonal_matrix_vector_products.resize(n_block_types, n_block_types, 0);

    // Fill in any null subsidiary preconditioners
    this->fill_in_subsidiary_preconditioners(n_block_types);

    // The total time for extracting all the blocks from the "global" matrix
    double t_extraction_total = 0.0;

    // The total time for setting up the subsidiary preconditioners
    double t_subsidiary_setup_total = 0.0;

    // The total time for setting up the matrix-vector products
    double t_mvp_setup_total = 0.0;

    // Build the preconditioners and matrix vector products
    for (unsigned i = 0; i < n_block_types; i++)
    {
      // See if the i-th subsidiary preconditioner is a block preconditioner
      if (dynamic_cast<BlockPreconditioner<CRDoubleMatrix>*>(
            this->Subsidiary_preconditioner_pt[i]) != 0)
      {
        // If we need to compute the memory statistics
        if (Compute_memory_statistics)
        {
          // If we're dealing with a GMRES block preconditioner
          if (dynamic_cast<GMRESBlockPreconditioner*>(
                this->Subsidiary_preconditioner_pt[i]))
          {
            // Enable the doc-ing of memory usage in the GMRES block
            // preconditioner
            dynamic_cast<GMRESBlockPreconditioner*>(
              this->Subsidiary_preconditioner_pt[i])
              ->enable_doc_memory_usage();
          }
        } // if (Compute_memory_statistics)

        // Get the start time
        double t_subsidiary_setup_start = TimingHelpers::timer();

        // It's a block preconditioner so pass it a pointer to the global
        // matrix!
        this->Subsidiary_preconditioner_pt[i]->setup(this->matrix_pt());

        // Get the end time
        double t_subsidiary_setup_end = TimingHelpers::timer();

        // Update the timing total
        t_subsidiary_setup_total +=
          (t_subsidiary_setup_end - t_subsidiary_setup_start);
      }
      // It's just a preconditioner so pass it a deep copy of the block!
      else
      {
        // Get the start time
        double t_extract_start = TimingHelpers::timer();

        // Grab the i-th diagonal block
        CRDoubleMatrix block_matrix = this->get_block(i, i);

        // Get the end time
        double t_extract_end = TimingHelpers::timer();

        // Update the timing total
        t_extraction_total += (t_extract_end - t_extract_start);

        // Get the start time
        double t_subsidiary_setup_start = TimingHelpers::timer();

        // Set up the i-th subsidiary preconditioner with this block
        this->Subsidiary_preconditioner_pt[i]->setup(&block_matrix);

        // Get the end time
        double t_subsidiary_setup_end = TimingHelpers::timer();

        // Update the timing total
        t_subsidiary_setup_total +=
          (t_subsidiary_setup_end - t_subsidiary_setup_start);
      }

      // The index of the lower column bound
      unsigned l;

      // The index of the upper column bound
      unsigned u;

      // If we're using the upper triangular form of the preconditioner
      if (Upper_triangular)
      {
        // The lower bound column bound
        l = i + 1;

        // If the block bandwidth has been provided
        if (Block_bandwidth >= 0)
        {
          // The upper bound is the last nonzero block column
          u = std::min(n_block_types, l + Block_bandwidth);
        }
        // Default case; loop over all the block columns
        else
        {
          // The upper bound is the last column
          u = n_block_types;
        }
      }
      // If we're using the lower triangular form of the preconditioner
      else
      {
        // The upper bound column bound
        u = i;

        // If the block bandwidth has been provided
        if (Block_bandwidth >= 0)
        {
          // The lower bound is the last nonzero block column. The explicit
          // casting is a bit over the top but it's better to be safe than
          // sorry...
          l = std::max(0, int(int(u) - Block_bandwidth));
        }
        // Default case; loop over all the block columns
        else
        {
          // The lower bound is the first column
          l = 0;
        }
      } // for (unsigned i=0;i<nblock_types;i++)

      // Loop over the chosen block columns
      for (unsigned j = l; j < u; j++)
      {
        // Get the start time
        double t_extract_start = TimingHelpers::timer();

        // Get the (i,j)-th block matrix
        CRDoubleMatrix block_matrix = this->get_block(i, j);

        // Do we need to doc. the memory statistics?
        if (Compute_memory_statistics)
        {
          // How many rows are there in this block?
          unsigned n_row = block_matrix.nrow();

          // How many nonzeros are there in this block?
          unsigned n_nnz = block_matrix.nnz();

          // Compute the memory usage. The only memory overhead here (for the
          // setup phase) is the storage of the sub-diagonal blocks for the MVPs
          Memory_usage_in_bytes += ((2 * ((n_row + 1) * sizeof(int))) +
                                    (n_nnz * (sizeof(int) + sizeof(double))));
        }

        // Get the end time
        double t_extract_end = TimingHelpers::timer();

        // Update the timing total
        t_extraction_total += (t_extract_end - t_extract_start);

        // Get the start time
        double t_mvp_start = TimingHelpers::timer();

        // Copy the block into a "multiplier" class. If trilinos is being
        // used this should also be faster than oomph-lib's multiphys.
        Off_diagonal_matrix_vector_products(i, j) = new MatrixVectorProduct();

        // Set the damn thing up
        this->setup_matrix_vector_product(
          Off_diagonal_matrix_vector_products(i, j), &block_matrix, j);

        // Get the end time
        double t_mvp_end = TimingHelpers::timer();

        // Update the timing total
        t_mvp_setup_total += (t_mvp_end - t_mvp_start);
      }
    } // for (unsigned i=0;i<n_block_types;i++)

    // Remember that the preconditioner has been set up
    Preconditioner_has_been_setup = true;

    // If we're meant to build silently, reassign the oomph stream pointer
    if (this->Silent_preconditioner_setup == true)
    {
      // Store the output stream pointer
      oomph_info.stream_pt() = this->Stream_pt;

      // Reset our own stream pointer
      this->Stream_pt = 0;
    }

    // Tell the user
    oomph_info << "Total block extraction time [sec]: " << t_extraction_total
               << "\nTotal subsidiary preconditioner setup time [sec]: "
               << t_subsidiary_setup_total
               << "\nTotal matrix-vector product setup time [sec]: "
               << t_mvp_setup_total << std::endl;

    // Do we need to doc. the memory statistics?
    // NOTE: We're going to assume that either:
    //    (1) GMRES is used as a subsidiary block preconditioner (preconditioned
    //        by the Navier-Stokes subsidiary block preconditioner), or;
    //    (2) SuperLU is used as the subsidiary preconditioner.
    if (Compute_memory_statistics)
    {
      // Allocate space for the total memory usage
      double total_memory_usage_for_setup_phase = 0.0;

      // How many rows are there in the global Jacobian?
      unsigned n_row = this->matrix_pt()->nrow();

      // How many nonzeros are there in the global Jacobian?
      unsigned n_nnz = this->matrix_pt()->nnz();

      // Storage for the memory usage of the global system matrix.
      // NOTE: This calculation is done by taking into account the storage
      // scheme for a compressed row matrix
      double memory_usage_for_storage_of_global_jacobian_in_bytes =
        ((2 * ((n_row + 1) * sizeof(int))) +
         (n_nnz * (sizeof(int) + sizeof(double))));

      // Allocate storage for the memory usage of the subsidiary preconditioner
      double memory_usage_for_subsidiary_preconditioner_in_bytes = 0.0;

      // Allocate storage for the memory usage of the Navier-Stokes subsidiary
      // block preconditioner (if it's used, that is)
      double memory_usage_for_nssbp_preconditioner_in_bytes = 0.0;

      // Boolean to indicate whether we've issue a warning about not being
      // able to compute memory statistics. If it doesn't work with one
      // block, chances are it won't work with any of them, so don't be a
      // pain in the gluteus maximus and just spit out the warning once
      bool have_issued_warning_about_memory_calculations = false;

      // Build the preconditioners and matrix vector products
      for (unsigned i = 0; i < n_block_types; i++)
      {
        // See if the i-th subsidiary preconditioner the GMRES block
        // preconditioner
        if (dynamic_cast<GMRESBlockPreconditioner*>(
              this->Subsidiary_preconditioner_pt[i]) != 0)
        {
          // Upcast the subsidiary preconditioner pointer to a GMRES block
          // preconditioner pointer
          GMRESBlockPreconditioner* gmres_block_prec_pt =
            dynamic_cast<GMRESBlockPreconditioner*>(
              this->Subsidiary_preconditioner_pt[i]);

          // Update the subsidiary preconditioner memory usage variable
          memory_usage_for_subsidiary_preconditioner_in_bytes +=
            gmres_block_prec_pt->get_memory_usage_in_bytes();

          // Add on the memory usage for the NS subsidiary block preconditioner
          memory_usage_for_nssbp_preconditioner_in_bytes +=
            gmres_block_prec_pt->navier_stokes_subsidiary_preconditioner_pt()
              ->get_memory_usage_in_bytes();
        }
        // This block is preconditioned by the SuperLU preconditioner
        else if (dynamic_cast<SuperLUPreconditioner*>(
                   this->Subsidiary_preconditioner_pt[i]) != 0)
        {
          // Update the subsidiary preconditioner memory usage variable
          memory_usage_for_subsidiary_preconditioner_in_bytes +=
            dynamic_cast<SuperLUPreconditioner*>(
              this->Subsidiary_preconditioner_pt[i])
              ->get_total_memory_needed_for_superlu();
        }
        // Don't know what to do otherwise!
        else
        {
          // Have we already complained once?
          if (!have_issued_warning_about_memory_calculations)
          {
            // Remember that we've issued a warning now
            have_issued_warning_about_memory_calculations = true;

            // Allocate storage for an output stream
            std::ostringstream warning_message_stream;

            // Create a warning message
            warning_message_stream
              << "Can't compute the memory statistics for the " << i
              << "-th diagonal block in the banded\nblock "
              << "triangular preconditioner so I'm just going "
              << "to skip it." << std::endl;

            // Give the user a warning
            OomphLibWarning(warning_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
          }
        }
      } // for (unsigned i=0;i<nblock_types;i++)

      // Create some space
      oomph_info << std::endl;

      // If we've got a nonzero contribution from the NSSBP
      if (memory_usage_for_nssbp_preconditioner_in_bytes > 0.0)
      {
        // Doc. the memory usage for the NSSBP
        oomph_info << "Amount of memory used by Navier-Stokes subsidiary block "
                   << "preconditioner (MB): "
                   << memory_usage_for_nssbp_preconditioner_in_bytes / 1.0e+06
                   << std::endl;

        // Add in the NSSBP contribution
        total_memory_usage_for_setup_phase +=
          memory_usage_for_nssbp_preconditioner_in_bytes;
      }

      // Add in the subsidiary preconditioners contribution
      total_memory_usage_for_setup_phase +=
        memory_usage_for_subsidiary_preconditioner_in_bytes;

      // Add in the banded block triangular preconditioner contribution
      total_memory_usage_for_setup_phase += get_memory_usage_in_bytes();

      // Add in the global Jacobian contribution
      total_memory_usage_for_setup_phase +=
        memory_usage_for_storage_of_global_jacobian_in_bytes;

      // How much memory have we used in the subsidiary preconditioners?
      oomph_info
        << "Amount of memory used by the subsidiary preconditioners "
        << "(MB): "
        << memory_usage_for_subsidiary_preconditioner_in_bytes / 1.0e+06
        << "\nAmount of memory used by banded block triangular "
        << "preconditioner (MB): " << get_memory_usage_in_bytes() / 1.0e+06
        << "\nAmount of memory used to store (global) Jacobian (MB): "
        << memory_usage_for_storage_of_global_jacobian_in_bytes / 1.0e+06
        << "\nTotal amount of memory being used after setup (MB): "
        << total_memory_usage_for_setup_phase / 1.0e+06 << "\n"
        << std::endl;
    } // if (Compute_memory_statistics)
  } // End of setup

  //=============================================================================
  /// Preconditioner solve for the block triangular preconditioner
  //=============================================================================
  template<typename MATRIX>
  void BandedBlockTriangularPreconditioner<MATRIX>::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
    // Cache number of block types
    unsigned n_block = this->nblock_types();

    // The start index
    int start;

    // The end index
    int end;

    // The max. bandwidth of the block matrix
    int bandwidth_limit;

    // The steps to take (=+1/-1) which is dependent on whether we're using
    // the lower or upper triangular form of the preconditioner
    int step;

    // If we're using this as an block upper triangular preconditioner
    if (Upper_triangular)
    {
      start = n_block - 1;
      end = -1;
      step = -1;
    }
    else
    {
      start = 0;
      end = n_block;
      step = 1;
    }

    // Storage for the iteration count. The entries will be updated as we loop
    // over the blocks. Here -1 is used to indicate that a block was not solved
    // by GMRES so we ignore it
    Vector<int> iteration_count(n_block, -1);

    // Storage for the number of blocks solved with GMRES
    unsigned n_block_solved_with_gmres = 0;

    // Vector of vectors for each section of residual vector
    Vector<DoubleVector> block_r;

    // Rearrange the vector r into the vector of block vectors block_r
    this->get_block_vectors(r, block_r);

    // Vector of vectors for the solution block vectors
    Vector<DoubleVector> block_z(n_block);

    // Loop over the block rows
    for (int i = start; i != end; i += step)
    {
      // If the i-th subsidiary preconditioner is a block preconditioner
      if (dynamic_cast<BlockPreconditioner<CRDoubleMatrix>*>(
            this->Subsidiary_preconditioner_pt[i]) != 0)
      {
        // This is pretty inefficient but there is no other way to do this with
        // block sub preconditioners as far as I can tell because they demand
        // to have the full r and z vectors...
        DoubleVector block_z_with_size_of_full_z(r.distribution_pt());
        DoubleVector r_updated(r.distribution_pt());

        // Construct the new r vector (with the update given by back subs
        // below).
        this->return_block_vectors(block_r, r_updated);

        // Solve (blocking is handled inside the block preconditioner).
        this->Subsidiary_preconditioner_pt[i]->preconditioner_solve(
          r_updated, block_z_with_size_of_full_z);

        // Extract this block's z vector because we need it for the next step
        // (and throw away block_z_with_size_of_full_z).
        this->get_block_vector(i, block_z_with_size_of_full_z, block_z[i]);
      }
      // If the subsidiary preconditioner is just a regular preconditioner
      else
      {
        // Solve on the block
        this->Subsidiary_preconditioner_pt[i]->preconditioner_solve(block_r[i],
                                                                    block_z[i]);
      }

      // See if the i-th subsidiary preconditioner the GMRES block
      // preconditioner
      if (dynamic_cast<GMRESBlockPreconditioner*>(
            this->Subsidiary_preconditioner_pt[i]) != 0)
      {
        // Update the iteration count
        iteration_count[i] = dynamic_cast<GMRESBlockPreconditioner*>(
                               this->Subsidiary_preconditioner_pt[i])
                               ->iterations();

        // Increment the counter
        n_block_solved_with_gmres++;
      }

      // If the block bandwidth has been provided
      if (Block_bandwidth >= 0)
      {
        // If we're using this as an block upper triangular preconditioner
        if (Upper_triangular)
        {
          bandwidth_limit = std::max(end - i, step + step * Block_bandwidth);
        }
        else
        {
          bandwidth_limit = std::min(end - i, step + step * Block_bandwidth);
        }
      }
      // Default case; loop over all the block rows
      else
      {
        bandwidth_limit = end - i;
      }

      // Substitute (over all the rows)
      for (int j = i + step; j != i + bandwidth_limit; j += step)
      {
        // Allocate space for the matrix-vector product (MVP)
        DoubleVector temp;

        // Calculate the MVP
        Off_diagonal_matrix_vector_products(j, i)->multiply(block_z[i], temp);

        // Now update the RHS vector
        block_r[j] -= temp;
      }
    } // for (int i=start;i!=end;i+=step)

    // Copy the solution from the block vector block_z back into z
    this->return_block_vectors(block_z, z);

    // Only compute the iteration count statistics if we have some data
    if (n_block_solved_with_gmres > 0)
    {
      // Storage for the average iteration count
      double n_iter_mean = 0.0;

      // Storage for the variance of the iteration count
      double n_iter_var = 0.0;

      // Loop over the entries of the iteration count vector
      for (unsigned i = 0; i < n_block; i++)
      {
        // If we used GMRES on this block
        if (dynamic_cast<GMRESBlockPreconditioner*>(
              this->Subsidiary_preconditioner_pt[i]) != 0)
        {
          // Sanity check: make sure there's a non-negative iteration count
          if (iteration_count[i] < 0)
          {
            // Create an output stream
            std::ostringstream error_message_stream;

            // Create an error message
            error_message_stream
              << "Iteration count should be non-negative but "
              << "you have an iteration\ncount of " << iteration_count[i] << "!"
              << std::endl;

            // Throw an error
            throw OomphLibError(error_message_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        } // if (dynamic_cast<GMRESBlockPreconditioner*>...
      } // for (unsigned i=0;i<n_block;i++)

      // Loop over the entries of the iteration count vector
      for (unsigned i = 0; i < n_block; i++)
      {
        // If we used GMRES on this block
        if (dynamic_cast<GMRESBlockPreconditioner*>(
              this->Subsidiary_preconditioner_pt[i]) != 0)
        {
          // Update the mean
          n_iter_mean +=
            (double(iteration_count[i]) / double(n_block_solved_with_gmres));
        } // if (dynamic_cast<GMRESBlockPreconditioner*>...
      } // for (unsigned i=0;i<n_block;i++)

      // Loop over the entries of the iteration count vector (again)
      // NOTE: We don't need to do the error checks again, we did them on the
      // previous run through
      for (unsigned i = 0; i < n_block; i++)
      {
        // If we used GMRES on this block
        if (dynamic_cast<GMRESBlockPreconditioner*>(
              this->Subsidiary_preconditioner_pt[i]) != 0)
        {
          // Update the variance
          n_iter_var +=
            (std::pow(double(iteration_count[i]) - n_iter_mean, 2.0) /
             double(n_block_solved_with_gmres));
        } // if (dynamic_cast<GMRESBlockPreconditioner*>...
      } // for (unsigned i=0;i<n_block;i++)

      // Doc. the statistics
      oomph_info << "\nNumber of subsidiary blocks solved with GMRES: "
                 << n_block_solved_with_gmres
                 << "\nAverage subsidiary preconditioner iteration count: "
                 << n_iter_mean
                 << "\nStandard deviation of subsidiary preconditioner "
                 << "iteration count: " << std::sqrt(n_iter_var) << std::endl;
    } // if (n_block_solved_with_gmres>0)
  } // End of preconditioner_solve

  // Ensure build of required objects (BUT only for CRDoubleMatrix objects)
  template class ExactDGPBlockPreconditioner<CRDoubleMatrix>;
  template class BandedBlockTriangularPreconditioner<CRDoubleMatrix>;
} // End of namespace oomph
