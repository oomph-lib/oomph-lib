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
// Header file
#include "spacetime_navier_stokes_block_preconditioner.h"

/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////

namespace oomph
{
  /*
 #ifdef OOMPH_HAS_HYPRE
  //======start_of_HypreSubsidiaryPreconditionerHelper_namespace==============
  /// Helper method for the block diagonal F block preconditioner to allow
  /// hypre to be used as a subsidiary block preconditioner
  //==========================================================================
  namespace HypreSubsidiaryPreconditionerHelper
  {
   /// Assign the Hypre preconditioner pointer
   Preconditioner* set_hypre_preconditioner()
   {
    return new HyprePreconditioner;
   }
  } // End of HypreSubsidiaryPreconditionerHelper
 #endif
  */

  //============================================================================
  /// Setup for the block triangular preconditioner
  //============================================================================
  void SpaceTimeNavierStokesSubsidiaryPreconditioner::setup()
  {
    // For debugging...
    bool doc_block_matrices = false;

    // Clean up any previously created data
    clean_up_memory();

    // If we're meant to build silently
    if (this->Silent_preconditioner_setup == true)
    {
      // Store the output stream pointer
      this->Stream_pt = oomph_info.stream_pt();

      // Now set the oomph_info stream pointer to the null stream to
      // disable all possible output
      oomph_info.stream_pt() = &oomph_nullstream;
    }

#ifdef PARANOID
    // If the upcast was unsuccessful
    if (dynamic_cast<CRDoubleMatrix*>(this->matrix_pt()) == 0)
    {
      // Allocate space for an error message
      std::ostringstream error_message;

      // Create the error message
      error_message << "NavierStokesSchurComplementPreconditioner only works "
                    << "with CRDoubleMatrix matrices" << std::endl;

      // Throw the error message
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Start the timer for the block setup
    double t_block_setup_start = TimingHelpers::timer();

    // Set up block look up schemes (done automatically in the
    // BlockPreconditioner base class, based on the information
    // provided in the block-preconditionable elements in the problem)

    // This preconditioner has two types of block:
    //        Type 0: Velocity - Corresponding to DOF type 0 & 1
    //        Type 1: Pressure - Corresponding to DOF type 2
    unsigned n_dof_types = 0;

    // If this has been set up as a subsidiary preconditioner
    if (this->is_subsidiary_block_preconditioner())
    {
      // Get the number of dof types (will be told by the master preconditioner)
      n_dof_types = this->ndof_types();

      // DRAIG: Delete
      oomph_info << "Number of dof types: " << n_dof_types << std::endl;

      // Make sure there's only 3 dof types for now!
      if (n_dof_types != 3)
      {
        // Allocate space for an error message
        std::ostringstream error_message;

        // Create the error message
        error_message << "Should only be 3 dof types! You have " << n_dof_types
                      << " dof types!" << std::endl;

        // Throw an error
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
    else
    {
      // Throw an error
      throw OomphLibError("Can only be used as a subsidiary preconditioner!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    } // if (this->is_subsidiary_block_preconditioner())

    // Allocate storage for the dof-to-block mapping
    Vector<unsigned> dof_to_block_map(n_dof_types);

    // The last dof type (the pressure) should have it's own block type
    dof_to_block_map[n_dof_types - 1] = 1;

    // Call the generic block setup function
    this->block_setup(dof_to_block_map);

    // Stop the timer
    double t_block_setup_finish = TimingHelpers::timer();

    // Calculate the time difference
    double block_setup_time = t_block_setup_finish - t_block_setup_start;

    // Document the time taken
    oomph_info << "Time for block_setup(...) [sec]: " << block_setup_time
               << std::endl;

    // Get the number of block types
    unsigned n_block_type = this->nblock_types();

    // How many block types are there? Should only be 2...
    oomph_info << "Number of block types: " << n_block_type << std::endl;

    // Allocate storage for the momentum block
    CRDoubleMatrix* f_pt = new CRDoubleMatrix;

    // Allocate storage for the divergence block
    CRDoubleMatrix* d_pt = new CRDoubleMatrix;

    // Allocate storage for the gradient block
    CRDoubleMatrix* g_pt = new CRDoubleMatrix;

    // Allocate storage for the pressure Poisson matrix
    CRDoubleMatrix* p_matrix_pt = new CRDoubleMatrix;

    // Get the momentum block F
    this->get_block(0, 0, *f_pt);

    // Get the divergence block, D
    this->get_block(1, 0, *d_pt);

    // Get gradient matrix, G=D^T
    this->get_block(0, 1, *g_pt);

    // Calculate the pressure Poisson-like matrix
    d_pt->multiply(*g_pt, *p_matrix_pt);

    // If the user wants the block matrices themselves
    if (doc_block_matrices)
    {
      // Filename suffix (defaults should be an empty string)
      std::string suffix = "";

      // Create space for all of the blocks
      DenseMatrix<CRDoubleMatrix*> block_matrix_pt(
        n_block_type, n_block_type, 0);

      // Loop over the block rows
      for (unsigned i = 0; i < n_block_type; i++)
      {
        // Loop over the block columns
        for (unsigned j = 0; j < n_block_type; j++)
        {
          // Allocate space for the matrix
          block_matrix_pt(i, j) = new CRDoubleMatrix;

          // Get the block
          this->get_block(i, j, *block_matrix_pt(i, j));
        }
      } // for (unsigned i=0;i<n_block_type;i++)

      // Space for the full matrix
      CRDoubleMatrix* j_pt = new CRDoubleMatrix;

      // Concatenate the sub-blocks to get the full matrix
      CRDoubleMatrixHelpers::concatenate(block_matrix_pt, *j_pt);

      // Output the matrix in MATLAB form
      j_pt->sparse_indexed_output_with_offset("j_" + suffix + ".csv");

      // Tell the user
      oomph_info << "Done output of j_" + suffix + ".csv" << std::endl;

      // Loop over the block rows
      for (unsigned i = 0; i < n_block_type; i++)
      {
        // Loop over the block columns
        for (unsigned j = 0; j < n_block_type; j++)
        {
          // Clear up this block
          delete block_matrix_pt(i, j);

          // Make it a null pointer
          block_matrix_pt(i, j) = 0;
        }
      } // for (unsigned i=0;i<n_block_type;i++)

      // Now clear up j_pt
      delete j_pt;

      // Make it a null pointer
      j_pt = 0;

      // Output the matrix in MATLAB form
      f_pt->sparse_indexed_output_with_offset("f_matrix" + suffix + ".dat");

      // Tell the user
      oomph_info << "Done output of f_matrix" + suffix + ".dat" << std::endl;

      // Output the matrix in MATLAB form
      d_pt->sparse_indexed_output_with_offset("d_matrix" + suffix + ".dat");

      // Tell the user
      oomph_info << "Done output of d_matrix" + suffix + ".dat" << std::endl;

      // Output the matrix in MATLAB form
      g_pt->sparse_indexed_output_with_offset("g_matrix" + suffix + ".dat");

      // Tell the user
      oomph_info << "Done output of g_matrix" + suffix + ".dat" << std::endl;

      // Output the matrix in MATLAB form
      p_matrix_pt->sparse_indexed_output_with_offset("p_matrix" + suffix +
                                                     ".dat");

      // Tell the user
      oomph_info << "Done output of p_matrix" + suffix + ".dat" << std::endl;

      // Exit
      exit(0);
    }

    // Kill divergence matrix because we don't need it any more
    delete d_pt;

    // Make it a null pointer
    d_pt = 0;

    // Make a new matrix-vector product
    F_mat_vec_pt = new MatrixVectorProduct;

    // Make a new matrix-vector product
    G_mat_vec_pt = new MatrixVectorProduct;

    // Set the matrix-vector product up
    this->setup_matrix_vector_product(F_mat_vec_pt, f_pt, 0);

    // Set the matrix-vector product up
    this->setup_matrix_vector_product(G_mat_vec_pt, g_pt, 1);

    // Do we need to doc. the memory statistics?
    if (Compute_memory_statistics)
    {
      // Initialise the memory usage variable
      Memory_usage_in_bytes = 0.0;

      // How many rows are there in the momentum block?
      unsigned n_row_f = f_pt->nrow();

      // How many nonzeros are there in the momentum block?
      unsigned n_nnz_f = f_pt->nnz();

      // Update the memory usage variable. NOTE: This calculation is done by
      // taking into account the storage scheme for a compressed row matrix
      Memory_usage_in_bytes += ((2 * ((n_row_f + 1) * sizeof(int))) +
                                (n_nnz_f * (sizeof(int) + sizeof(double))));

      // How many rows are there in the gradient block?
      unsigned n_row_g = g_pt->nrow();

      // How many nonzeros are there in the gradient block?
      unsigned n_nnz_g = g_pt->nnz();

      // Update the memory usage variable
      Memory_usage_in_bytes += ((2 * ((n_row_g + 1) * sizeof(int))) +
                                (n_nnz_g * (sizeof(int) + sizeof(double))));
    }

    // Kill the gradient matrix because we don't need it any more
    delete g_pt;

    // Make it a null pointer
    g_pt = 0;

    // If the P preconditioner has not been set
    if (P_preconditioner_pt == 0)
    {
      // Just use SuperLU
      P_preconditioner_pt = new SuperLUPreconditioner;

      // Indicate that we're using the default preconditioner
      Using_default_p_preconditioner = true;

      // Do we need to doc. the memory statistics?
      if (Compute_memory_statistics)
      {
        // Add on the memory needed to store and calculate the LU factors of P
        Memory_usage_in_bytes +=
          dynamic_cast<SuperLUPreconditioner*>(P_preconditioner_pt)
            ->get_total_memory_needed_for_superlu();
      }
    } // if (P_preconditioner_pt==0)

    // Compute the LU decomposition of P
    P_preconditioner_pt->setup(p_matrix_pt);

    // Now delete the actual matrix P, we don't need it anymore
    delete p_matrix_pt;

    // Make it a null pointer
    p_matrix_pt = 0;

    // If the F preconditioner has not been setup
    if (F_preconditioner_pt == 0)
    {
      // Just use SuperLU
      F_preconditioner_pt = new SuperLUPreconditioner;

      // Indicate that we're using the default preconditioner
      Using_default_f_preconditioner = true;

      // Do we need to doc. the memory statistics?
      if (Compute_memory_statistics)
      {
        // Add on the memory needed to store and calculate the LU factors of P
        Memory_usage_in_bytes +=
          dynamic_cast<SuperLUPreconditioner*>(F_preconditioner_pt)
            ->get_total_memory_needed_for_superlu();
      }
    } // if (F_preconditioner_pt==0)

    // Compute the LU decomposition of F
    F_preconditioner_pt->setup(f_pt);

    // Now delete the matrix, we don't need it anymore
    delete f_pt;

    // Make it a null pointer
    f_pt = 0;

    // Remember that the preconditioner has been setup so the stored
    // information can be wiped when we come here next...
    Preconditioner_has_been_setup = true;

    // If we're meant to build silently, reassign the oomph stream pointer
    if (this->Silent_preconditioner_setup == true)
    {
      // Store the output stream pointer
      oomph_info.stream_pt() = this->Stream_pt;

      // Reset our own stream pointer
      this->Stream_pt = 0;
    }
  } // End of setup

  //=============================================================================
  /// Preconditioner solve for the block triangular preconditioner
  //=============================================================================
  void SpaceTimeNavierStokesSubsidiaryPreconditioner::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
    // If z is not set up then give it the same distribution
    if (!z.distribution_pt()->built())
    {
      // Build z with the same distribution as r
      z.build(r.distribution_pt(), 0.0);
    }

    // -----------------------------------------------------------------------
    // Step 1 - apply approximate Schur inverse to pressure unknowns (block 1)
    // -----------------------------------------------------------------------
    // First working vectors
    DoubleVector temp_vec;

    // Second working vector
    DoubleVector another_temp_vec;

    // Copy pressure values from residual vector to temp_vec:
    // Loop over all entries in the global vector (this one
    // includes velocity and pressure dofs in some random fashion)
    this->get_block_vector(1, r, temp_vec);

    // Compute the vector P^{-1}r_p
    P_preconditioner_pt->preconditioner_solve(temp_vec, another_temp_vec);

    // Now clear the vector storing r_p
    temp_vec.clear();

    // Now compute G(P^{-1}r_p)
    G_mat_vec_pt->multiply(another_temp_vec, temp_vec);

    // Clear the vector storing P^{-1}r_p
    another_temp_vec.clear();

    // Now update to get F(GP^{-1}r_p)
    F_mat_vec_pt->multiply(temp_vec, another_temp_vec);

    // Clear the vector storing GP^{-1}r_p
    temp_vec.clear();

    // Now compute D(FGP^{-1}r_p)=G^{T}(FGP^{-1}r_p)
    G_mat_vec_pt->multiply_transpose(another_temp_vec, temp_vec);

    // Clear the vector storing FGP^{-1}r_p
    another_temp_vec.clear();

    // Finally, compute the full approximation to -z_p;
    // -z_p=P^{-1}(DFG)P^{-1}r_p NOTE: We'll move the sign over in the next step
    P_preconditioner_pt->preconditioner_solve(temp_vec, another_temp_vec);

    // Rebuild temp_vec with another_temp_vec's distribution with all
    // entries initialised to zero
    temp_vec.build(another_temp_vec.distribution_pt(), 0.0);

    // Now fix the sign and store the result in temp_vec
    temp_vec -= another_temp_vec;

    // Now copy temp_vec (i.e. z_p) back into the global vector z
    return_block_vector(1, temp_vec, z);

    // ------------------------------------------------------------
    // Step 2 - apply preconditioner to velocity unknowns (block 0)
    // ------------------------------------------------------------
    // Clear the vector storing z_p=-P^{-1}(DFG)P^{-1}r_p
    temp_vec.clear();

    // Recall that another_temp_vec (computed above) contains the negative
    // of the solution of the Schur complement system, -z_p. Multiply by G
    // and store the result in temp_vec (vector resizes itself).
    G_mat_vec_pt->multiply(another_temp_vec, temp_vec);

    // Now clear the vector storing P^{-1}(DFG)P^{-1}r_p
    another_temp_vec.clear();

    // Get the residual vector entries associated with the velocities, r_u
    get_block_vector(0, r, another_temp_vec);

    // Update the result to get r_u-Gz_p
    another_temp_vec += temp_vec;

    // Compute the solution z_u which satisfies z_u=F^{-1}(r_u-Gz_p)
    F_preconditioner_pt->preconditioner_solve(another_temp_vec, temp_vec);

    // Now store the result in the global solution vector
    return_block_vector(0, temp_vec, z);
  } // End of preconditioner_solve


  //============================================================================
  /// Setup for the GMRES block preconditioner
  //============================================================================
  void GMRESBlockPreconditioner::setup()
  {
    // Clean up any previously created data
    clean_up_memory();

    // If we're meant to build silently
    if (this->Silent_preconditioner_setup == true)
    {
      // Store the output stream pointer
      this->Stream_pt = oomph_info.stream_pt();

      // Now set the oomph_info stream pointer to the null stream to
      // disable all possible output
      oomph_info.stream_pt() = &oomph_nullstream;
    }

    // Start the timer for the block setup
    double t_block_setup_start = TimingHelpers::timer();

    // Set up block look up schemes (done automatically in the
    // BlockPreconditioner base class, based on the information
    // provided in the block-preconditionable elements in the problem)

    // This preconditioner has two types of block:
    //        Type 0: Velocity - Corresponding to DOF type 0 & 1
    //        Type 1: Pressure - Corresponding to DOF type 2
    unsigned n_dof_types = 0;

    // If this has been set up as a subsidiary preconditioner
    if (this->is_subsidiary_block_preconditioner())
    {
      // Get the number of dof types (will be told by the master preconditioner)
      n_dof_types = this->ndof_types();

      // Output the number of dof types
      oomph_info << "Number of dof types: " << n_dof_types << std::endl;

      // Make sure there's only 3 dof types for now!
      if (n_dof_types != 3)
      {
        // Allocate space for an error message
        std::ostringstream error_message;

        // Create the error message
        error_message << "Should only be 3 dof types! You have " << n_dof_types
                      << " dof types!" << std::endl;

        // Throw an error
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
    // If it's being used as a master preconditioner
    else
    {
      // Throw an error
      throw OomphLibError("Currently only used as a subsidiary preconditioner!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    } // if (this->is_subsidiary_block_preconditioner())

    // Aggregrate everything together since GMRES itself isn't going to take
    // advantage of the block structure...
    Vector<unsigned> dof_to_block_map(n_dof_types, 0);

    // Call the generic block setup function
    this->block_setup(dof_to_block_map);

    // Stop the timer
    double t_block_setup_finish = TimingHelpers::timer();

    // Calculate the time difference
    double block_setup_time = t_block_setup_finish - t_block_setup_start;

    // Document the time taken
    oomph_info << "Time for block_setup(...) [sec]: " << block_setup_time
               << std::endl;

    // Get the number of block types
    unsigned n_block_type = this->nblock_types();

    // How many block types are there? Should only be 2...
    oomph_info << "Number of block types: " << n_block_type << std::endl;

    // Space for the concatenated block matrix
    Matrix_pt = new CRDoubleMatrix;

    // Get the (one and only) block to be used by this block preconditioner
    this->get_block(0, 0, *Matrix_pt);

    // Create a new instance of the NS subsidiary preconditioner
    Navier_stokes_subsidiary_preconditioner_pt =
      new SpaceTimeNavierStokesSubsidiaryPreconditioner;

    // Do we need to doc. the memory statistics?
    if (Compute_memory_statistics)
    {
      // How many rows are there in this block?
      unsigned n_row = Matrix_pt->nrow();

      // How many nonzeros are there in this block?
      unsigned n_nnz = Matrix_pt->nnz();

      // Compute the memory usage. The only memory overhead here (for the
      // setup phase) is storing the system matrix
      Memory_usage_in_bytes = ((2 * ((n_row + 1) * sizeof(int))) +
                               (n_nnz * (sizeof(int) + sizeof(double))));

      // We might as well doc. the memory statistics of the Navier-Stokes
      // subsidiary block preconditioner while we're at it...
      Navier_stokes_subsidiary_preconditioner_pt->enable_doc_memory_usage();
    }

    // Create a map to declare which of the 3 dof types in the GMRES
    // subsidiary preconditioner correspond to the dof types in the NS
    // subsidiary block preconditioner
    Vector<unsigned> dof_map(n_dof_types);

    // Loop over the dof types associated with the GMRES subsidiary block
    // preconditioner
    for (unsigned i = 0; i < n_dof_types; i++)
    {
      // There is, essentially, an identity mapping between the GMRES subsidiary
      // block preconditioner and the NSSP w.r.t. the dof types
      dof_map[i] = i;
    }

    // Now turn it into an "proper" subsidiary preconditioner
    Navier_stokes_subsidiary_preconditioner_pt
      ->turn_into_subsidiary_block_preconditioner(this, dof_map);

    // Setup: The NS subsidiary preconditioner is itself a block preconditioner
    // so we need to pass it a pointer to full-size matrix!
    Navier_stokes_subsidiary_preconditioner_pt->setup(this->matrix_pt());

    // Remember that the preconditioner has been setup so the stored
    // information can be wiped when we come here next...
    Preconditioner_has_been_setup = true;

    // If we're meant to build silently, reassign the oomph stream pointer
    if (this->Silent_preconditioner_setup == true)
    {
      // Store the output stream pointer
      oomph_info.stream_pt() = this->Stream_pt;

      // Reset our own stream pointer
      this->Stream_pt = 0;
    }
  } // End of setup

  //=============================================================================
  /// Preconditioner solve for the GMRES block preconditioner
  //=============================================================================
  void GMRESBlockPreconditioner::preconditioner_solve(const DoubleVector& rhs,
                                                      DoubleVector& solution)
  {
    // Get the number of block types
    unsigned n_block_types = this->nblock_types();

    // If there is more than one block type (there shouldn't be!)
    if (n_block_types != 1)
    {
      // Throw an error
      throw OomphLibError("Can only deal with one dof type!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Allocate space for the solution blocks
    Vector<DoubleVector> block_z(n_block_types);

    // Allocate space for the residual blocks
    Vector<DoubleVector> block_r(n_block_types);

    // Get the reordered RHS vector
    this->get_block_vectors(rhs, block_r);

    // Get the reordered solution vector
    this->get_block_vectors(solution, block_z);

    // Initialise the entries in the solution vector to zero
    block_z[0].initialise(0.0);

    // Get number of dofs
    unsigned n_dof = block_r[0].nrow();

    // Time solver
    double t_start = TimingHelpers::timer();

    // Relative residual
    double resid;

    // Iteration counter
    unsigned iter = 1;

    // Initialise vectors
    Vector<double> s(n_dof + 1, 0);
    Vector<double> cs(n_dof + 1);
    Vector<double> sn(n_dof + 1);
    DoubleVector w(block_r[0].distribution_pt(), 0.0);

    // Storage for the time taken for all preconditioner applications
    double t_prec_application_total = 0.0;

    // Storage for the RHS vector
    DoubleVector r(block_r[0].distribution_pt(), 0.0);

    // If we're using LHS preconditioning then solve Mr=b-Jx for r (assumes x=0)
    if (Preconditioner_LHS)
    {
      // Initialise timer
      double t_prec_application_start = TimingHelpers::timer();

      // This is pretty inefficient but there is no other way to do this with
      // block subsidiary preconditioners as far as I can tell because they
      // demand to have the full r and z vectors...
      DoubleVector block_z_with_size_of_full_z(rhs.distribution_pt(), 0.0);
      DoubleVector r_updated(rhs.distribution_pt(), 0.0);

      // Reorder the entries in block_r[0] and put them back into the
      // global-sized vector. The subsidiary preconditioner should only
      // operate on the dofs that this preconditioner operates on so
      // don't worry about all the other entries
      this->return_block_vector(0, block_r[0], r_updated);

      // Solve on the global-sized vectors
      Navier_stokes_subsidiary_preconditioner_pt->preconditioner_solve(
        r_updated, block_z_with_size_of_full_z);

      // Now grab the part of the solution associated with this preconditioner
      this->get_block_vector(0, block_z_with_size_of_full_z, block_z[0]);

      // Doc time for setup of preconditioner
      double t_prec_application_end = TimingHelpers::timer();

      // Update the preconditioner application time total
      t_prec_application_total +=
        (t_prec_application_end - t_prec_application_start);
    }
    // If we're using RHS preconditioning, do nothing to the actual RHS vector
    else
    {
      // Copy the entries into r
      r = block_r[0];
    }

    // Calculate the initial residual norm (assumes x=0 initially)
    double normb = r.norm();

    // Set beta (the initial residual)
    double beta = normb;

    // Compute initial relative residual
    if (normb == 0.0)
    {
      // To avoid dividing by zero, set normb to 1
      normb = 1;
    }

    // Now calculate the relative residual norm
    resid = beta / normb;

    // If required, document convergence history to screen or file (if
    // stream is open)
    if (Doc_convergence_history)
    {
      // If a file has not been provided to output to
      if (!Output_file_stream.is_open())
      {
        // Output the initial relative residual norm to screen
        oomph_info << 0 << " " << resid << std::endl;
      }
      // If there is a file provided to output to
      else
      {
        // Output the initial relative residual norm to file
        Output_file_stream << 0 << " " << resid << std::endl;
      }
    } // if (Doc_convergence_history)

    // If GMRES converges immediately
    if (resid <= Tolerance)
    {
      // If the user wishes GMRES to output information about the solve to
      // screen
      if (Doc_time)
      {
        // Tell the user
        oomph_info << "GMRES block preconditioner converged immediately. "
                   << "Normalised residual norm: " << resid << std::endl;
      }

      // Now reorder the values in block_z[0] and put them back into the
      // global solution vector
      this->return_block_vector(0, block_z[0], solution);

      // Doc time for solver
      double t_end = TimingHelpers::timer();
      Solution_time = t_end - t_start;

      // If the user wishes GMRES to output information about the solve to
      // screen
      if (Doc_time)
      {
        // Tell the user
        oomph_info << "Time for all preconditioner applications [sec]: "
                   << t_prec_application_total
                   << "\nTotal time for solve with GMRES block preconditioner "
                   << "[sec]: " << Solution_time << std::endl;
      }

      // We're finished; end here!
      return;
    } // if (resid<=Tolerance)

    // Initialise vector of orthogonal basis vectors, v
    Vector<DoubleVector> v(n_dof + 1);

    // Initialise the upper Hessenberg matrix H.
    // NOTE: For implementation purposes the upper hessenberg matrix indices
    // are swapped so the matrix is effectively transposed
    Vector<Vector<double>> H(n_dof + 1);

    // Loop as many times as we're allowed
    while (iter <= Max_iter)
    {
      // Copy the entries of the residual vector, r, into v[0]
      v[0] = r;

      // Now update the zeroth basis vector, v[0], to have values r/beta
      v[0] /= beta;

      // Storage for ...?
      s[0] = beta;

      // Inner iteration counter for restarted version (we don't actually use
      // a restart in this implementation of GMRES...)
      unsigned iter_restart;

      // Perform iterations
      for (iter_restart = 0; iter_restart < n_dof && iter <= Max_iter;
           iter_restart++, iter++)
      {
        // Resize next column of the upper Hessenberg matrix
        H[iter_restart].resize(iter_restart + 2);

        // Solve for w inside the scope of these braces
        {
          // Temporary vector
          DoubleVector temp(block_r[0].distribution_pt(), 0.0);

          // If we're using LHS preconditioning solve Mw=Jv[i] for w
          if (Preconditioner_LHS)
          {
            // Calculate temp=Jv[i]
            Matrix_pt->multiply(v[iter_restart], temp);

            // Get the current time
            double t_prec_application_start = TimingHelpers::timer();

            // This is pretty inefficient but there is no other way to do this
            // with block sub preconditioners as far as I can tell because they
            // demand to have the full r and z vectors...
            DoubleVector block_z_with_size_of_full_z(rhs.distribution_pt(),
                                                     0.0);
            DoubleVector r_updated(rhs.distribution_pt(), 0.0);

            // Put the values in temp into the global-sized vector
            this->return_block_vector(0, temp, r_updated);

            // Apply the NSSP to the global-sized vectors (this will extract
            // the appropriate sub-vectors and solve on the reduced system)
            Navier_stokes_subsidiary_preconditioner_pt->preconditioner_solve(
              r_updated, block_z_with_size_of_full_z);

            // Now grab the part of the solution associated with this
            // preconditioner
            this->get_block_vector(0, block_z_with_size_of_full_z, w);

            // End timer
            double t_prec_application_end = TimingHelpers::timer();

            // Update the preconditioner application time total
            t_prec_application_total +=
              (t_prec_application_end - t_prec_application_start);
          }
          // If we're using RHS preconditioning solve
          else
          {
            // Initialise timer
            double t_prec_application_start = TimingHelpers::timer();

            // This is pretty inefficient but there is no other way to do this
            // with block sub preconditioners as far as I can tell because they
            // demand to have the full r and z vectors...
            DoubleVector block_z_with_size_of_full_z(rhs.distribution_pt());
            DoubleVector r_updated(rhs.distribution_pt());

            // Put the values in v[iter_restart] back into the global-sized
            // vector
            this->return_block_vector(0, v[iter_restart], r_updated);

            // Solve on the global-sized vectors
            Navier_stokes_subsidiary_preconditioner_pt->preconditioner_solve(
              r_updated, block_z_with_size_of_full_z);

            // Now grab the part of the solution associated with this
            // preconditioner
            this->get_block_vector(0, block_z_with_size_of_full_z, temp);

            // Doc time for setup of preconditioner
            double t_prec_application_end = TimingHelpers::timer();

            // Update the preconditioner application time total
            t_prec_application_total +=
              (t_prec_application_end - t_prec_application_start);

            // Calculate w=Jv_m where v_m=M^{-1}v
            Matrix_pt->multiply(temp, w);
          }
        } // End of solve for w

        // Get a pointer to the values in w
        double* w_pt = w.values_pt();

        // Loop over the columns of the (transposed) Hessenberg matrix
        for (unsigned k = 0; k <= iter_restart; k++)
        {
          // Initialise the entry in the upper Hessenberg matrix
          H[iter_restart][k] = 0.0;

          // Get the pointer to the values of the k-th basis vector
          double* vk_pt = v[k].values_pt();

          // Update H with the dot product of w and v[k]
          H[iter_restart][k] += w.dot(v[k]);

          // Loop over the entries in w and v[k] again
          for (unsigned i = 0; i < n_dof; i++)
          {
            // Now update the values in w
            w_pt[i] -= H[iter_restart][k] * vk_pt[i];
          }
        } // for (unsigned k=0;k<=iter_restart;k++)

        // Calculate the (iter_restart+1,iter_restart)-th entry in the upper
        // Hessenberg matrix (remember, H is actually the transpose of the
        // proper upper Hessenberg matrix)
        H[iter_restart][iter_restart + 1] = w.norm();

        // Build the (iter_restart+1)-th basis vector by copying w
        v[iter_restart + 1] = w;

        // Now rescale the basis vector
        v[iter_restart + 1] /= H[iter_restart][iter_restart + 1];

        // Loop over the columns of the transposed upper Hessenberg matrix
        for (unsigned k = 0; k < iter_restart; k++)
        {
          // Apply a Givens rotation to entries of H
          apply_plane_rotation(
            H[iter_restart][k], H[iter_restart][k + 1], cs[k], sn[k]);
        }

        // Now generate the cs and sn entries for a plane rotation
        generate_plane_rotation(H[iter_restart][iter_restart],
                                H[iter_restart][iter_restart + 1],
                                cs[iter_restart],
                                sn[iter_restart]);

        // Use the newly generated entries in cs and sn to apply a Givens
        // rotation to those same entries
        apply_plane_rotation(H[iter_restart][iter_restart],
                             H[iter_restart][iter_restart + 1],
                             cs[iter_restart],
                             sn[iter_restart]);

        // Now apply the Givens rotation to the entries in s
        apply_plane_rotation(s[iter_restart],
                             s[iter_restart + 1],
                             cs[iter_restart],
                             sn[iter_restart]);

        // Compute the current residual
        beta = std::fabs(s[iter_restart + 1]);

        // Compute the relative residual norm
        resid = beta / normb;

        // If required, document the convergence history to screen or file
        if (Doc_convergence_history)
        {
          // If an output file has not been provided
          if (!Output_file_stream.is_open())
          {
            // Output the convergence history to screen
            oomph_info << iter << " " << resid << std::endl;
          }
          // An output file has been provided so output to that
          else
          {
            // Output the convergence history to file
            Output_file_stream << iter << " " << resid << std::endl;
          }
        } // if (Doc_convergence_history)

        // If the required tolerance was met
        if (resid < Tolerance)
        {
          // Storage for the global-sized solution vector. Strictly speaking
          // we could actually just use the vector, solution but this can be
          // done after we know we've implemented everything correctly...
          DoubleVector block_z_with_size_of_full_z(rhs.distribution_pt(), 0.0);

          // Take the current solution, reorder the entries appropriately
          // and stick them into the global sized vector
          this->return_block_vector(0, block_z[0], block_z_with_size_of_full_z);

          // This will update block_z[0] and won't touch the global-size vector
          // block_x_with_size_of_full_x. We're in charge of reordering the
          // entries and putting it in solution
          update(
            iter_restart, H, s, v, block_z_with_size_of_full_z, block_z[0]);

          // Now go into block_z[0], reorder the entries in the appropriate
          // manner and put them into the vector, solution
          this->return_block_vector(0, block_z[0], solution);

          // If the user wishes GMRES to document the convergence
          if (Doc_time)
          {
            // Document the convergence to screen
            oomph_info
              << "\nGMRES block preconditioner converged (1). Normalised "
              << "residual norm: " << resid
              << "\nNumber of iterations to convergence: " << iter << "\n"
              << std::endl;
          }

          // Get the current time
          double t_end = TimingHelpers::timer();

          // Calculate the time difference, i.e. the time taken for the whole
          // solve
          Solution_time = t_end - t_start;

          // Storage the iteration count
          Iterations = iter;

          // If the user wishes GMRES to document the convergence
          if (Doc_time)
          {
            // Document the convergence to screen
            oomph_info
              << "Time for all preconditioner applications [sec]: "
              << t_prec_application_total
              << "\nTotal time for solve with GMRES block preconditioner "
              << "[sec]: " << Solution_time << std::endl;
          }

          // We're done now; finish here
          return;
        } // if (resid<Tolerance)
      } // for (iter_restart=0;iter_restart<n_dof&&iter<=Max_iter;...

      // Update
      if (iter_restart > 0)
      {
        // Storage for the global-sized solution vector
        DoubleVector block_z_with_size_of_full_z(rhs.distribution_pt(), 0.0);

        // Take the current solution, reorder the entries appropriately
        // and stick them into the global sized vector
        this->return_block_vector(0, block_z[0], block_z_with_size_of_full_z);

        // Update the (reordered block) solution vector
        update(
          (iter_restart - 1), H, s, v, block_z_with_size_of_full_z, block_z[0]);

        // Now go into block_z[0], reorder the entries in the appropriate manner
        // and put them into the vector, solution
        this->return_block_vector(0, block_z[0], solution);
      }

      // Solve Mr=(b-Jx) for r
      {
        // Temporary vector to store b-Jx
        DoubleVector temp(block_r[0].distribution_pt(), 0.0);

        // Calculate temp=b-Jx
        Matrix_pt->residual(block_z[0], block_r[0], temp);

        // If we're using LHS preconditioning
        if (Preconditioner_LHS)
        {
          // Initialise timer
          double t_prec_application_start = TimingHelpers::timer();

          // This is pretty inefficient but there is no other way to do this
          // with block sub preconditioners as far as I can tell because they
          // demand to have the full r and z vectors...
          DoubleVector block_z_with_size_of_full_z(rhs.distribution_pt());
          DoubleVector r_updated(rhs.distribution_pt());

          // Reorder the entries of temp and put them into the global-sized
          // vector
          this->return_block_vector(0, temp, r_updated);

          // Solve on the global-sized vectors
          Navier_stokes_subsidiary_preconditioner_pt->preconditioner_solve(
            r_updated, block_z_with_size_of_full_z);

          // Now grab the part of the solution associated with this
          // preconditioner
          this->get_block_vector(0, block_z_with_size_of_full_z, r);

          // Doc time for setup of preconditioner
          double t_prec_application_end = TimingHelpers::timer();

          // Update the preconditioner application time total
          t_prec_application_total +=
            (t_prec_application_end - t_prec_application_start);
        }
      } // End of solve Mr=(b-Jx) for r

      // Compute the current residual norm
      beta = r.norm();

      // if relative residual within tolerance
      resid = beta / normb;
      if (resid < Tolerance)
      {
        // Now store the result in the global solution vector
        this->return_block_vector(0, block_z[0], solution);

        if (Doc_time)
        {
          oomph_info
            << "\nGMRES block preconditioner converged (2). Normalised "
            << "residual norm: " << resid
            << "\nNumber of iterations to convergence: " << iter << "\n"
            << std::endl;
        }

        // Doc time for solver
        double t_end = TimingHelpers::timer();
        Solution_time = t_end - t_start;

        Iterations = iter;

        if (Doc_time)
        {
          oomph_info
            << "Time for all preconditioner applications [sec]: "
            << t_prec_application_total
            << "\nTotal time for solve with GMRES block preconditioner "
            << "[sec]: " << Solution_time << std::endl;
        }
        return;
      }
    }


    // otherwise GMRES failed convergence
    oomph_info << "\nGMRES block preconditioner did not converge to required "
               << "tolerance! \nReturning with normalised residual norm: "
               << resid << "\nafter " << Max_iter << " iterations.\n"
               << std::endl;

    if (Throw_error_after_max_iter)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "Solver failed to converge and you requested an "
                           << "error on convergence failures.";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
    }
  } // End of solve_helper
} // End of namespace oomph
