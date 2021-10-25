// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
#include "hypre_solver.h"

// For problem->get_jacobian(...)
#include "problem.h"


namespace oomph
{
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //==================================================================
  /// Default settings for various uses of the HYPRE solver
  //==================================================================
  namespace Hypre_default_settings
  {
    /// Set default parameters for use as preconditioner in
    /// 2D Poisson-type problem.
    void set_defaults_for_2D_poisson_problem(
      HyprePreconditioner* hypre_preconditioner_pt)
    {
      // Set iterations to 1
      hypre_preconditioner_pt->set_amg_iterations(1);

      // Use simple smoother
      hypre_preconditioner_pt->amg_using_simple_smoothing();

      // Smoother types:
      //           0=Jacobi
      //           1=Gauss-Seidel
      hypre_preconditioner_pt->amg_simple_smoother() = 1;

      // AMG preconditioner
      hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

      // Choose strength parameter for amg
      hypre_preconditioner_pt->amg_strength() = 0.25;

      // Coarsening type
      hypre_preconditioner_pt->amg_coarsening() = 0;
    }


    /// Set default parameters for use as preconditioner in
    /// 3D Poisson-type problem.
    void set_defaults_for_3D_poisson_problem(
      HyprePreconditioner* hypre_preconditioner_pt)
    {
      // Set default settings as for 2D Poisson
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

      // Change strength parameter for amg
      hypre_preconditioner_pt->amg_strength() = 0.7;
    }


    /// Set default parameters for use as preconditioner in
    /// for momentum block in Navier-Stokes problem
    void set_defaults_for_navier_stokes_momentum_block(
      HyprePreconditioner* hypre_preconditioner_pt)
    {
      // Set default settings as for 2D Poisson
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

      // Change smoother type:
      //           0=Jacobi
      //           1=Gauss-Seidel
      hypre_preconditioner_pt->amg_simple_smoother() = 0;

      // Set smoother damping
      hypre_preconditioner_pt->amg_damping() = 0.5;

      // Change strength parameter for amg
      hypre_preconditioner_pt->amg_strength() = 0.75;
    }

  } // namespace Hypre_default_settings


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  // functions for HypreHelpers namespace
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  namespace HypreHelpers
  {
    //========================================================================
    /// Default for AMG strength (0.25 recommended for 2D problems;
    /// larger (0.5-0.75, say) for 3D
    //========================================================================
    double AMG_strength = 0.25;

    //========================================================================
    /// Default AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    ///  3 = modified RS with 3rd pass to add C points on the boundaries
    ///  6 = Falgout (uses 1 then CLJP using interior coarse points as
    ///      first independent set)
    ///  8 = PMIS (parallel coarsening using independent sets - lower
    ///      complexities than 0, maybe also slower convergence)
    ///  10= HMIS (one pass RS on each processor then PMIS on interior
    ///      coarse points as first independent set)
    ///  11= One pass RS on each processor (not recommended)
    //========================================================================
    unsigned AMG_coarsening = 6;


    //========================================================================
    /// AMG interpolation truncation factor
    //========================================================================
    double AMG_truncation = 0.0;

    //========================================================================
    /// Helper function to check the Hypre error flag, return the message
    /// associated with any error, and reset the global error flag to zero.
    /// This function also returns the error value.
    //========================================================================
    int check_HYPRE_error_flag(std::ostringstream& message)
    {
      // get the Hypre error flag
      int err = HYPRE_GetError();

      // Tell us all about it...
      if (err)
      {
        oomph_info << "Hypre error flag=" << err << std::endl;
        char* error_message = new char[128];
        HYPRE_DescribeError(err, error_message);
        message << "WARNING: " << std::endl
                << "HYPRE error message: " << error_message << std::endl;
        delete[] error_message;
      }

      // reset Hypre's global error flag
      hypre__global_error = 0;

      return err;
    }


    //========================================================================
    /// Helper function to create a HYPRE_IJVector and HYPRE_ParVector.
    /// + If no MPI then serial vectors are created
    /// + If MPI and serial input vector then distributed hypre vectors are
    /// created
    /// + If MPI and distributed input vector the distributed output vectors
    ///   are created.
    //========================================================================
    void create_HYPRE_Vector(const DoubleVector& oomph_vec,
                             const LinearAlgebraDistribution* dist_pt,
                             HYPRE_IJVector& hypre_ij_vector,
                             HYPRE_ParVector& hypre_par_vector)
    {
      // the lower and upper row of the vector on this processor
      unsigned lower = dist_pt->first_row();
      unsigned upper = dist_pt->first_row() + dist_pt->nrow_local() - 1;

      // number of local rows
      unsigned nrow_local = dist_pt->nrow_local();

      // initialize Hypre vector
#ifdef OOMPH_HAS_MPI
      HYPRE_IJVectorCreate(
        oomph_vec.distribution_pt()->communicator_pt()->mpi_comm(),
        lower,
        upper,
        &hypre_ij_vector);
#else
      HYPRE_IJVectorCreate(MPI_COMM_WORLD, lower, upper, &hypre_ij_vector);
#endif
      HYPRE_IJVectorSetObjectType(hypre_ij_vector, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(hypre_ij_vector);

      // set up array containing indices
      int* indices = new int[nrow_local];
      double* values = new double[nrow_local];
      const unsigned hypre_first_row = dist_pt->first_row();
      unsigned j = 0;
      if (!oomph_vec.distributed() && dist_pt->distributed())
      {
        j = hypre_first_row;
      }
      const double* o_pt = oomph_vec.values_pt();
      for (unsigned i = 0; i < nrow_local; i++)
      {
        indices[i] = hypre_first_row + i;
        values[i] = o_pt[j + i];
      }

      // insert values
      HYPRE_IJVectorSetValues(hypre_ij_vector, nrow_local, indices, values);

      // assemble vectors
      HYPRE_IJVectorAssemble(hypre_ij_vector);
      HYPRE_IJVectorGetObject(hypre_ij_vector, (void**)&hypre_par_vector);

      // clean up
      delete[] indices;
      delete[] values;
    }


    //========================================================================
    /// Helper function to create a HYPRE_IJVector and HYPRE_ParVector.
    /// + If no MPI then serial vectors are created
    /// + If MPI and serial input vector then distributed hypre vectors are
    /// created
    /// + If MPI and distributed input vector the distributed output vectors
    ///   are created.
    //========================================================================
    void create_HYPRE_Vector(const LinearAlgebraDistribution* dist_pt,
                             HYPRE_IJVector& hypre_ij_vector,
                             HYPRE_ParVector& hypre_par_vector)
    {
      // the lower and upper row of the vector on this processor
      unsigned lower = dist_pt->first_row();
      unsigned upper = dist_pt->first_row() + dist_pt->nrow_local() - 1;

      // initialize Hypre vector
#ifdef OOMPH_HAS_MPI
      HYPRE_IJVectorCreate(
        dist_pt->communicator_pt()->mpi_comm(), lower, upper, &hypre_ij_vector);
#else
      HYPRE_IJVectorCreate(MPI_COMM_WORLD, lower, upper, &hypre_ij_vector);
#endif
      HYPRE_IJVectorSetObjectType(hypre_ij_vector, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(hypre_ij_vector);

      // assemble vectors
      HYPRE_IJVectorAssemble(hypre_ij_vector);
      HYPRE_IJVectorGetObject(hypre_ij_vector, (void**)&hypre_par_vector);
    }


    //========================================================================
    /// Helper function to create a serial HYPRE_IJMatrix and
    /// HYPRE_ParCSRMatrix from a CRDoubleMatrix
    /// NOTE: dist_pt is rebuilt to match the distribution of the hypre solver
    /// which is not necassarily the same as the oomph lib matrix
    //========================================================================
    void create_HYPRE_Matrix(CRDoubleMatrix* oomph_matrix,
                             HYPRE_IJMatrix& hypre_ij_matrix,
                             HYPRE_ParCSRMatrix& hypre_par_matrix,
                             LinearAlgebraDistribution* dist_pt)
    {
#ifdef PARANOID
      // check that the matrix is built
      if (!oomph_matrix->built())
      {
        std::ostringstream error_message;
        error_message << "The matrix has not been built";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      // check the matrix is square
      if (oomph_matrix->nrow() != oomph_matrix->ncol())
      {
        std::ostringstream error_message;
        error_message << "create_HYPRE_Matrix require a square matrix. "
                      << "Matrix is " << oomph_matrix->nrow() << " by "
                      << oomph_matrix->ncol() << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

#ifdef OOMPH_HAS_MPI
      // Trap the case when we have compiled with MPI,
      // but are running in serial
      if (!MPI_Helpers::mpi_has_been_initialised())
      {
        std::ostringstream error_stream;
        error_stream
          << "Oomph-lib has been compiled with MPI support and "
          << "you are using HYPRE.\n"
          << "For this combination of flags, MPI must be initialised.\n"
          << "Call MPI_Helpers::init() in the "
          << "main() function of your driver code\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // find number of rows/columns
      const unsigned nrow = int(oomph_matrix->nrow());

      // get pointers to the matrix

      // column indices of matrix
      const int* matrix_cols = oomph_matrix->column_index();

      // entries of matrix
      const double* matrix_vals = oomph_matrix->value();

      // row starts
      const int* matrix_row_start = oomph_matrix->row_start();

      // build the distribution
      if (oomph_matrix->distribution_pt()->distributed())
      {
        dist_pt->build(oomph_matrix->distribution_pt());
      }
      else
      {
        bool distributed = true;
        if (oomph_matrix->distribution_pt()->communicator_pt()->nproc() == 1)
        {
          distributed = false;
        }
        dist_pt->build(oomph_matrix->distribution_pt()->communicator_pt(),
                       nrow,
                       distributed);
      }

      // initialize hypre matrix
      unsigned lower = dist_pt->first_row();
      unsigned upper = lower + dist_pt->nrow_local() - 1;

#ifdef OOMPH_HAS_MPI
      HYPRE_IJMatrixCreate(dist_pt->communicator_pt()->mpi_comm(),
                           lower,
                           upper,
                           lower,
                           upper,
                           &hypre_ij_matrix);
#else
      HYPRE_IJMatrixCreate(
        MPI_COMM_WORLD, lower, upper, lower, upper, &hypre_ij_matrix);
#endif
      HYPRE_IJMatrixSetObjectType(hypre_ij_matrix, HYPRE_PARCSR);
      HYPRE_IJMatrixInitialize(hypre_ij_matrix);

      // set up a row map
      // and first row / nrow_local
      const unsigned hypre_nrow_local = dist_pt->nrow_local();
      const unsigned hypre_first_row = dist_pt->first_row();
      int* ncols_per_row = new int[hypre_nrow_local];
      int* row_map = new int[hypre_nrow_local];
      for (unsigned i = 0; i < hypre_nrow_local; i++)
      {
        unsigned j = i;
        if (!oomph_matrix->distributed() && dist_pt->distributed())
        {
          j += hypre_first_row;
        }
        ncols_per_row[i] = matrix_row_start[j + 1] - matrix_row_start[j];
        row_map[i] = hypre_first_row + i;
      }

      // put values in HYPRE matrix
      int local_start = 0;
      if (!oomph_matrix->distributed() && dist_pt->distributed())
      {
        local_start += matrix_row_start[hypre_first_row];
      }


      HYPRE_IJMatrixSetValues(hypre_ij_matrix,
                              hypre_nrow_local,
                              ncols_per_row,
                              row_map,
                              matrix_cols + local_start,
                              matrix_vals + local_start);

      // assemble matrix
      HYPRE_IJMatrixAssemble(hypre_ij_matrix); // hierher leak?
      HYPRE_IJMatrixGetObject(hypre_ij_matrix, (void**)&hypre_par_matrix);

      // tidy up memory
      delete[] ncols_per_row;
      delete[] row_map;
    }

    //====================================================================
    /// Helper function to set Euclid options using a command line
    /// like array.
    //=====================================================================
    void euclid_settings_helper(const bool& use_block_jacobi,
                                const bool& use_row_scaling,
                                const bool& use_ilut,
                                const int& level,
                                const double& drop_tol,
                                const int& print_level,
                                HYPRE_Solver& euclid_object)
    {
      // Easier to use C-arrays rather than std::strings because Euclid takes
      // char** as argument and because C++ doesn't provide decent number to
      // std::string conversion functions.

      int n_args = 0;
      const char* args[22];

      // first argument is empty string
      args[n_args++] = "";

      // switch on/off Block Jacobi ILU
      args[n_args++] = "-bj";
      if (use_block_jacobi)
      {
        args[n_args++] = "1";
      }
      else
      {
        args[n_args++] = "0";
      }

      // switch on/off row scaling
      args[n_args++] = "-rowScale";
      if (use_row_scaling)
      {
        args[n_args++] = "1";
      }
      else
      {
        args[n_args++] = "0";
      }

      // set level for ILU(k)
      args[n_args++] = "-level";
      char level_value[10];
      sprintf(level_value, "%d", level);
      args[n_args++] = level_value;

      // // set drop tol for ILU(k) factorization
      // args[n_args++] = "-sparseA";
      // char droptol[20];
      // sprintf(droptol,"%f",drop_tol);
      // args[n_args++] = droptol;

      // // set ILUT factorization if required
      // if (use_ilut)
      //  {
      //   args[n_args++] = "-ilut";
      //   args[n_args++] = droptol;
      //  }

      // set printing of Euclid data
      if (print_level == 0)
      {
        args[n_args++] = "-eu_stats";
        args[n_args++] = "0";
        args[n_args++] = "-eu_mem";
        args[n_args++] = "0";
      }
      if (print_level == 1)
      {
        args[n_args++] = "-eu_stats";
        args[n_args++] = "1";
        args[n_args++] = "-eu_mem";
        args[n_args++] = "0";
      }
      if (print_level == 2)
      {
        args[n_args++] = "-eu_stats";
        args[n_args++] = "1";
        args[n_args++] = "-eu_mem";
        args[n_args++] = "1";
      }

      // set next entry in array to null
      args[n_args] = 0;

      // Send the parameters
      HYPRE_EuclidSetParams(euclid_object, n_args, const_cast<char**>(args));
    }

  } // namespace HypreHelpers


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // functions for HypreInterface class
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// Helper function which creates a Hypre matrix from a CRDoubleMatrix
  /// If OOMPH-LIB has been set up for MPI use, the Hypre matrix is
  /// distributed over the available processors.
  //=============================================================================
  void HypreInterface::hypre_matrix_setup(CRDoubleMatrix* matrix_pt)
  {
    // reset Hypre's global error flag
    hypre__global_error = 0;

    // issue warning if the matrix is small compared to the number of processors
    if (unsigned(2 * matrix_pt->distribution_pt()->communicator_pt()->nproc()) >
        matrix_pt->nrow())
    {
      oomph_info
        << "Warning: HYPRE based solvers may fail if 2*number of processors "
        << "is greater than the number of unknowns!" << std::endl;
    }

    // store the distribution
    // generate the Hypre matrix
    HypreHelpers::create_HYPRE_Matrix(
      matrix_pt, Matrix_ij, Matrix_par, Hypre_distribution_pt);

    // Output error messages if required
    if (Hypre_error_messages)
    {
      std::ostringstream message;
      int err = HypreHelpers::check_HYPRE_error_flag(message);
      if (err)
      {
        OomphLibWarning(message.str(),
                        "HypreSolver::hypre_matrix_setup()",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }

    // delete CRDoubleMatrix if required
    if (Delete_input_data)
    {
      matrix_pt->clear();
    }
  }


  //=============================================================================
  /// Sets up the solver data required for use in an oomph-lib
  /// LinearSolver or Preconditioner, once the Hypre matrix has been
  /// generated using hypre_matrix_setup(...).
  //=============================================================================
  void HypreInterface::hypre_solver_setup()
  {
    // Store time
    double t_start = TimingHelpers::timer();
    double t_end = 0;


    // reset Hypre's global error flag
    hypre__global_error = 0;

    // create dummy Hypre vectors which are required for setup
    HYPRE_IJVector dummy_sol_ij;
    HYPRE_ParVector dummy_sol_par;
    HYPRE_IJVector dummy_rhs_ij;
    HYPRE_ParVector dummy_rhs_par;
    HypreHelpers::create_HYPRE_Vector(
      Hypre_distribution_pt, dummy_sol_ij, dummy_sol_par);
    HypreHelpers::create_HYPRE_Vector(
      Hypre_distribution_pt, dummy_rhs_ij, dummy_rhs_par);

    // Set up internal preconditioner for CG, GMRES or BiCGSTAB
    // --------------------------------------------------------
    if ((Hypre_method >= CG) && (Hypre_method <= BiCGStab))
    {
      // AMG preconditioner
      if (Internal_preconditioner == BoomerAMG)
      {
        // set up BoomerAMG
        HYPRE_BoomerAMGCreate(&Preconditioner);
        HYPRE_BoomerAMGSetPrintLevel(Preconditioner, AMG_print_level);
        HYPRE_BoomerAMGSetMaxLevels(Preconditioner, AMG_max_levels);
        HYPRE_BoomerAMGSetMaxIter(Preconditioner, 1);
        HYPRE_BoomerAMGSetTol(Preconditioner, 0.0);
        HYPRE_BoomerAMGSetCoarsenType(Preconditioner, AMG_coarsening);
        HYPRE_BoomerAMGSetStrongThreshold(Preconditioner, AMG_strength);
        HYPRE_BoomerAMGSetMaxRowSum(Preconditioner, AMG_max_row_sum);
        HYPRE_BoomerAMGSetTruncFactor(Preconditioner, AMG_truncation);

        if (AMG_using_simple_smoothing)
        {
          HYPRE_BoomerAMGSetRelaxType(Preconditioner, AMG_simple_smoother);
          HYPRE_BoomerAMGSetNumSweeps(Preconditioner, AMG_smoother_iterations);

          // This one gives a memory leak
          // double * relaxweight = new double[AMG_max_levels];

          // This is how they do it in a hypre demo code
          double* relaxweight = hypre_CTAlloc(double, AMG_max_levels);

          for (unsigned i = 0; i < AMG_max_levels; i++)
          {
            relaxweight[i] = AMG_damping;
          }
          HYPRE_BoomerAMGSetRelaxWeight(Preconditioner, relaxweight);
        }
        else
        {
          HYPRE_BoomerAMGSetSmoothType(Preconditioner, AMG_complex_smoother);
          HYPRE_BoomerAMGSetSmoothNumLevels(Preconditioner, AMG_max_levels);
          HYPRE_BoomerAMGSetSmoothNumSweeps(Preconditioner,
                                            AMG_smoother_iterations);

          // If we are using Euclid then set up additional Euclid only options
          if (AMG_complex_smoother == 9)
          {
            HypreHelpers::euclid_settings_helper(
              AMGEuclidSmoother_use_block_jacobi,
              AMGEuclidSmoother_use_row_scaling,
              AMGEuclidSmoother_use_ilut,
              AMGEuclidSmoother_level,
              AMGEuclidSmoother_drop_tol,
              AMGEuclidSmoother_print_level,
              Preconditioner);
          }
        }

        Existing_preconditioner = BoomerAMG;
      }

      // Euclid preconditioner
      else if (Internal_preconditioner == Euclid)
      {
#ifdef OOMPH_HAS_MPI
        HYPRE_EuclidCreate(Hypre_distribution_pt->communicator_pt()->mpi_comm(),
                           &Preconditioner);
#else
        HYPRE_EuclidCreate(MPI_COMM_WORLD, &Preconditioner);
#endif

        // Set parameters
        HypreHelpers::euclid_settings_helper(Euclid_using_BJ,
                                             Euclid_rowScale,
                                             Euclid_using_ILUT,
                                             Euclid_level,
                                             Euclid_droptol,
                                             Euclid_print_level,
                                             Preconditioner);

        Existing_preconditioner = Euclid;
      }

      // ParaSails preconditioner
      else if (Internal_preconditioner == ParaSails)
      {
#ifdef OOMPH_HAS_MPI
        HYPRE_ParaSailsCreate(
          Hypre_distribution_pt->communicator_pt()->mpi_comm(),
          &Preconditioner);
#else
        HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &Preconditioner);
#endif
        HYPRE_ParaSailsSetSym(Preconditioner, ParaSails_symmetry);
        HYPRE_ParaSailsSetParams(
          Preconditioner, ParaSails_thresh, ParaSails_nlevel);
        HYPRE_ParaSailsSetFilter(Preconditioner, ParaSails_filter);
        Existing_preconditioner = ParaSails;
      }

      // check error flag
      if (Hypre_error_messages)
      {
        std::ostringstream message;
        int err = HypreHelpers::check_HYPRE_error_flag(message);
        if (err)
        {
          OomphLibWarning(message.str(),
                          "HypreSolver::hypre_setup()",
                          OOMPH_EXCEPTION_LOCATION);
        }
      }
    } // end of setting up internal preconditioner


    // set up solver
    // -------------
    t_start = TimingHelpers::timer();

    // AMG solver
    if (Hypre_method == BoomerAMG)
    {
      if (Output_info)
      {
        oomph_info << "Setting up BoomerAMG, ";
      }

      // set up BoomerAMG
      HYPRE_BoomerAMGCreate(&Solver);
      HYPRE_BoomerAMGSetPrintLevel(Solver, AMG_print_level);
      HYPRE_BoomerAMGSetMaxLevels(Solver, AMG_max_levels);
      HYPRE_BoomerAMGSetMaxIter(Solver, Max_iter);
      HYPRE_BoomerAMGSetTol(Solver, Tolerance);
      HYPRE_BoomerAMGSetCoarsenType(Solver, AMG_coarsening);
      HYPRE_BoomerAMGSetStrongThreshold(Solver, AMG_strength);
      HYPRE_BoomerAMGSetMaxRowSum(Solver, AMG_max_row_sum);
      HYPRE_BoomerAMGSetTruncFactor(Solver, AMG_truncation);

      if (AMG_using_simple_smoothing)
      {
        HYPRE_BoomerAMGSetRelaxType(Solver, AMG_simple_smoother);
        HYPRE_BoomerAMGSetNumSweeps(Solver, AMG_smoother_iterations);

        // This one gives a memory leak
        // double * relaxweight = new double[AMG_max_levels];

        // This is how they do it in a hypre demo code
        double* relaxweight = hypre_CTAlloc(double, AMG_max_levels);

        for (unsigned i = 0; i < AMG_max_levels; i++)
        {
          relaxweight[i] = AMG_damping;
        }
        HYPRE_BoomerAMGSetRelaxWeight(Solver, relaxweight);
      }
      else
      {
        HYPRE_BoomerAMGSetSmoothType(Solver, AMG_complex_smoother);
        HYPRE_BoomerAMGSetSmoothNumLevels(Solver, AMG_max_levels);
        HYPRE_BoomerAMGSetSmoothNumSweeps(Solver, AMG_smoother_iterations);

        /* Other settings
         * 6 &	Schwarz smoothers & HYPRE_BoomerAMGSetDomainType,
         * HYPRE_BoomerAMGSetOverlap, \\
         *  &  & HYPRE_BoomerAMGSetVariant, HYPRE_BoomerAMGSetSchwarzRlxWeight
         * \\
         * 7 &	Pilut & HYPRE_BoomerAMGSetDropTol, HYPRE_BoomerAMGSetMaxNzPerRow
         * \\
         * 8 &	ParaSails & HYPRE_BoomerAMGSetSym, HYPRE_BoomerAMGSetLevel, \\
         * &  &  HYPRE_BoomerAMGSetFilter, HYPRE_BoomerAMGSetThreshold \\
         * 9 &	Euclid & HYPRE_BoomerAMGSetEuclidFile \\
         */

        // If we are using Euclid then set up additional Euclid only options
        if (AMG_complex_smoother == 9)
        {
          HypreHelpers::euclid_settings_helper(
            AMGEuclidSmoother_use_block_jacobi,
            AMGEuclidSmoother_use_row_scaling,
            AMGEuclidSmoother_use_ilut,
            AMGEuclidSmoother_level,
            AMGEuclidSmoother_drop_tol,
            AMGEuclidSmoother_print_level,
            Preconditioner);
        }

        // Add any others here as required...
      }

      //     MemoryUsage::doc_memory_usage("before amg setup [solver]");
      //     MemoryUsage::insert_comment_to_continous_top("BEFORE AMG SETUP
      //     [SOLVER]");

      HYPRE_BoomerAMGSetup(Solver, Matrix_par, dummy_rhs_par, dummy_sol_par);

      //     MemoryUsage::doc_memory_usage("after amg setup [solver]");
      //     MemoryUsage::insert_comment_to_continous_top("AFTER AMG SETUP
      //     [SOLVER]");

      Existing_solver = BoomerAMG;
    }

    // Euclid solver
    else if (Hypre_method == Euclid)
    {
      if (Output_info)
      {
        oomph_info << "Setting up Euclid, ";
      }
#ifdef OOMPH_HAS_MPI
      HYPRE_EuclidCreate(Hypre_distribution_pt->communicator_pt()->mpi_comm(),
                         &Solver);
#else
      HYPRE_EuclidCreate(MPI_COMM_WORLD, &Solver);
#endif

      // Set parameters
      HypreHelpers::euclid_settings_helper(Euclid_using_BJ,
                                           Euclid_rowScale,
                                           Euclid_using_ILUT,
                                           Euclid_level,
                                           Euclid_droptol,
                                           Euclid_print_level,
                                           Preconditioner);

      HYPRE_EuclidSetup(Solver, Matrix_par, dummy_rhs_par, dummy_sol_par);
      Existing_solver = Euclid;
    }

    // ParaSails preconditioner
    else if (Hypre_method == ParaSails)
    {
      if (Output_info)
      {
        oomph_info << "Setting up ParaSails, ";
      }
#ifdef OOMPH_HAS_MPI
      HYPRE_ParaSailsCreate(
        Hypre_distribution_pt->communicator_pt()->mpi_comm(), &Solver);
#else
      HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &Solver);
#endif
      HYPRE_ParaSailsSetSym(Solver, ParaSails_symmetry);
      HYPRE_ParaSailsSetParams(Solver, ParaSails_thresh, ParaSails_nlevel);
      HYPRE_ParaSailsSetFilter(Solver, ParaSails_filter);

      HYPRE_ParaSailsSetup(Solver, Matrix_par, dummy_rhs_par, dummy_sol_par);
      Existing_solver = ParaSails;
    }

    // CG solver
    else if (Hypre_method == CG)
    {
      if (Output_info)
      {
        oomph_info << "Setting up CG, ";
      }

#ifdef OOMPH_HAS_MPI
      HYPRE_ParCSRPCGCreate(
        Hypre_distribution_pt->communicator_pt()->mpi_comm(), &Solver);
#else
      HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &Solver);
#endif
      HYPRE_PCGSetTol(Solver, Tolerance);
      HYPRE_PCGSetLogging(Solver, 0);
      HYPRE_PCGSetPrintLevel(Solver, Krylov_print_level);
      HYPRE_PCGSetMaxIter(Solver, Max_iter);

      // set preconditioner
      if (Internal_preconditioner == BoomerAMG) // AMG
      {
        if (Output_info)
        {
          oomph_info << " with BoomerAMG preconditioner, ";
        }

        HYPRE_PCGSetPrecond(Solver,
                            (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                            (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup,
                            Preconditioner);
      }
      else if (Internal_preconditioner == Euclid) // Euclid
      {
        if (Output_info)
        {
          oomph_info << " with Euclid ILU preconditioner, ";
        }

        HYPRE_PCGSetPrecond(Solver,
                            (HYPRE_PtrToSolverFcn)HYPRE_EuclidSolve,
                            (HYPRE_PtrToSolverFcn)HYPRE_EuclidSetup,
                            Preconditioner);
      }
      else if (Internal_preconditioner == ParaSails) // ParaSails
      {
        if (Output_info)
        {
          oomph_info << " with ParaSails approximate inverse preconditioner, ";
        }

        HYPRE_PCGSetPrecond(Solver,
                            (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSolve,
                            (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSetup,
                            Preconditioner);
      }
      else
      {
        if (Output_info)
        {
          oomph_info << " with no preconditioner";
        }
      }


      HYPRE_PCGSetup(Solver,
                     (HYPRE_Matrix)Matrix_par,
                     (HYPRE_Vector)dummy_rhs_par,
                     (HYPRE_Vector)dummy_sol_par);


      Existing_solver = CG;
    }

    // GMRES solver
    else if (Hypre_method == GMRES)
    {
      if (Output_info)
      {
        oomph_info << "Setting up GMRES";
      }

#ifdef OOMPH_HAS_MPI
      HYPRE_ParCSRGMRESCreate(
        Hypre_distribution_pt->communicator_pt()->mpi_comm(), &Solver);
#else
      HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &Solver);
#endif
      HYPRE_GMRESSetTol(Solver, Tolerance);
      HYPRE_GMRESSetKDim(Solver, Max_iter);
      HYPRE_GMRESSetLogging(Solver, 0);
      HYPRE_GMRESSetPrintLevel(Solver, Krylov_print_level);
      HYPRE_GMRESSetMaxIter(Solver, Max_iter);

      // set preconditioner
      if (Internal_preconditioner == BoomerAMG) // AMG
      {
        if (Output_info)
        {
          oomph_info << " with BoomerAMG preconditioner, ";
        }

        HYPRE_GMRESSetPrecond(Solver,
                              (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                              (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup,
                              Preconditioner);
      }
      else if (Internal_preconditioner == Euclid) // Euclid
      {
        if (Output_info)
        {
          oomph_info << " with Euclid ILU preconditioner, ";
        }

        HYPRE_GMRESSetPrecond(Solver,
                              (HYPRE_PtrToSolverFcn)HYPRE_EuclidSolve,
                              (HYPRE_PtrToSolverFcn)HYPRE_EuclidSetup,
                              Preconditioner);
      }
      else if (Internal_preconditioner == ParaSails) // ParaSails
      {
        if (Output_info)
        {
          oomph_info << " with ParaSails approximate inverse preconditioner, ";
        }

        HYPRE_GMRESSetPrecond(Solver,
                              (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSolve,
                              (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSetup,
                              Preconditioner);
      }
      else
      {
        if (Output_info)
        {
          oomph_info << " with no preconditioner";
        }
      }

      HYPRE_GMRESSetup(Solver,
                       (HYPRE_Matrix)Matrix_par,
                       (HYPRE_Vector)dummy_rhs_par,
                       (HYPRE_Vector)dummy_sol_par);

      Existing_solver = GMRES;
    }

    // BiCGStab solver
    else if (Hypre_method == BiCGStab)
    {
      if (Output_info)
      {
        oomph_info << "Setting up BiCGStab";
      }
#ifdef OOMPH_HAS_MPI
      HYPRE_ParCSRBiCGSTABCreate(
        Hypre_distribution_pt->communicator_pt()->mpi_comm(), &Solver);
#else
      HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &Solver);
#endif
      HYPRE_BiCGSTABSetTol(Solver, Tolerance);
      HYPRE_BiCGSTABSetLogging(Solver, 0);
      HYPRE_BiCGSTABSetPrintLevel(Solver, Krylov_print_level);
      HYPRE_BiCGSTABSetMaxIter(Solver, Max_iter);

      // set preconditioner
      if (Internal_preconditioner == BoomerAMG) // AMG
      {
        if (Output_info)
        {
          oomph_info << " with BoomerAMG preconditioner, ";
        }

        HYPRE_BiCGSTABSetPrecond(Solver,
                                 (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                                 (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup,
                                 Preconditioner);
      }
      else if (Internal_preconditioner == Euclid) // Euclid
      {
        if (Output_info)
        {
          oomph_info << " with Euclid ILU preconditioner, ";
        }

        HYPRE_BiCGSTABSetPrecond(Solver,
                                 (HYPRE_PtrToSolverFcn)HYPRE_EuclidSolve,
                                 (HYPRE_PtrToSolverFcn)HYPRE_EuclidSetup,
                                 Preconditioner);
      }
      else if (Internal_preconditioner == ParaSails) // ParaSails
      {
        if (Output_info)
        {
          oomph_info << " with ParaSails approximate inverse preconditioner, ";
        }

        HYPRE_BiCGSTABSetPrecond(Solver,
                                 (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSolve,
                                 (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSetup,
                                 Preconditioner);
      }
      else
      {
        if (Output_info)
        {
          oomph_info << " with no preconditioner, ";
        }
      }

      HYPRE_BiCGSTABSetup(Solver,
                          (HYPRE_Matrix)Matrix_par,
                          (HYPRE_Vector)dummy_rhs_par,
                          (HYPRE_Vector)dummy_sol_par);

      Existing_solver = BiCGStab;
    }

    // no solver exists for this value of Solver flag
    else
    {
      std::ostringstream error_message;
      error_message << "Solver has been set to an invalid value. "
                    << "current value=" << Solver;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    t_end = TimingHelpers::timer();
    double solver_setup_time = t_end - t_start;

    // destroy dummy hypre vectors
    HYPRE_IJVectorDestroy(dummy_sol_ij);
    HYPRE_IJVectorDestroy(dummy_rhs_ij);

    // check error flag
    if (Hypre_error_messages)
    {
      std::ostringstream message;
      int err = HypreHelpers::check_HYPRE_error_flag(message);
      if (err)
      {
        OomphLibWarning(message.str(),
                        "HypreSolver::hypre_solver_setup()",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }

    // output times
    if (Output_info)
    {
      oomph_info << "time for setup [s] : " << solver_setup_time << std::endl;
    }
  }


  //===================================================================
  /// Helper function performs a solve if solver data has been set
  /// up using hypre_solver_setup(...).
  //====================================================================
  void HypreInterface::hypre_solve(const DoubleVector& rhs,
                                   DoubleVector& solution)
  {
    // Record time
    double t_start = TimingHelpers::timer();

    // Set up hypre vectors
    // --------------------

    // Hypre vector for rhs values
    HYPRE_IJVector rhs_ij;
    HYPRE_ParVector rhs_par;

    // Hypre vector for solution
    HYPRE_IJVector solution_ij;
    HYPRE_ParVector solution_par;

    // set up rhs_values and vec_indices
    HypreHelpers::create_HYPRE_Vector(
      rhs, Hypre_distribution_pt, rhs_ij, rhs_par);

    HypreHelpers::create_HYPRE_Vector(
      Hypre_distribution_pt, solution_ij, solution_par);

    // check error flag
    if (Hypre_error_messages)
    {
      std::ostringstream message;
      int err = HypreHelpers::check_HYPRE_error_flag(message);
      if (err)
      {
        OomphLibWarning(message.str(),
                        "HypreSolver::hypre_solve()",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }

    // solve
    // -----

    // for solver stats
    int iterations = 0;
    double norm = 0;

    // Get the norm of rhs
    const double rhs_norm = rhs.norm();
    bool do_solving = false;
    if (rhs_norm > 0.0)
    {
      do_solving = true;
    }

#ifdef OOMPH_HAS_MPI
    // We need to check whether any processor requires to solve, if that
    // is the case then do the solving
    if (MPI_Helpers::mpi_has_been_initialised())
    {
      if (MPI_Helpers::communicator_pt()->nproc() > 1)
      {
        unsigned this_processor_do_solving = 0;
        unsigned all_processors_do_solving = 0;
        if (do_solving)
        {
          this_processor_do_solving = 1;
        }
        // Get the communicator
        OomphCommunicator* comm_pt = MPI_Helpers::communicator_pt();
        // Communicate with all procesoors
        MPI_Allreduce(&this_processor_do_solving,
                      &all_processors_do_solving,
                      1,
                      MPI_UNSIGNED,
                      MPI_SUM,
                      comm_pt->mpi_comm());
        if (all_processors_do_solving > 0)
        {
          do_solving = true;
        }
      }
    }
#endif

    if (do_solving)
    {
      if (Existing_solver == BoomerAMG)
      {
        HYPRE_BoomerAMGSolve(Solver, Matrix_par, rhs_par, solution_par);
        HYPRE_BoomerAMGGetNumIterations(Solver, &iterations);
        HYPRE_BoomerAMGGetFinalRelativeResidualNorm(Solver, &norm);
      }
      else if (Existing_solver == CG)
      {
        HYPRE_PCGSolve(Solver,
                       (HYPRE_Matrix)Matrix_par,
                       (HYPRE_Vector)rhs_par,
                       (HYPRE_Vector)solution_par);
        HYPRE_PCGGetNumIterations(Solver, &iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(Solver, &norm);
      }
      else if (Existing_solver == GMRES)
      {
        HYPRE_GMRESSolve(Solver,
                         (HYPRE_Matrix)Matrix_par,
                         (HYPRE_Vector)rhs_par,
                         (HYPRE_Vector)solution_par);
        HYPRE_GMRESGetNumIterations(Solver, &iterations);
        HYPRE_GMRESGetFinalRelativeResidualNorm(Solver, &norm);
      }
      else if (Existing_solver == BiCGStab)
      {
        HYPRE_BiCGSTABSolve(Solver,
                            (HYPRE_Matrix)Matrix_par,
                            (HYPRE_Vector)rhs_par,
                            (HYPRE_Vector)solution_par);
        HYPRE_BiCGSTABGetNumIterations(Solver, &iterations);
        HYPRE_BiCGSTABGetFinalRelativeResidualNorm(Solver, &norm);
      }
      else if (Existing_solver == Euclid)
      {
        HYPRE_EuclidSolve(Solver, Matrix_par, rhs_par, solution_par);
      }
      else if (Existing_solver == ParaSails)
      {
        HYPRE_ParaSailsSolve(Solver, Matrix_par, rhs_par, solution_par);
      }

      // output any error message
      if (Hypre_error_messages)
      {
        std::ostringstream message;
        int err = HypreHelpers::check_HYPRE_error_flag(message);
        if (err)
        {
          OomphLibWarning(message.str(),
                          "HypreSolver::hypre_solve()",
                          OOMPH_EXCEPTION_LOCATION);
        }
      }

    } // if (do_solving)

    // Copy result to solution
    unsigned nrow_local = Hypre_distribution_pt->nrow_local();
    unsigned first_row = Hypre_distribution_pt->first_row();
    int* indices = new int[nrow_local];
    for (unsigned i = 0; i < nrow_local; i++)
    {
      indices[i] = first_row + i;
    }
    LinearAlgebraDistribution* soln_dist_pt;
    if (solution.built())
    {
      soln_dist_pt = new LinearAlgebraDistribution(solution.distribution_pt());
    }
    else
    {
      soln_dist_pt = new LinearAlgebraDistribution(rhs.distribution_pt());
    }
    solution.build(Hypre_distribution_pt, 0.0);
    HYPRE_IJVectorGetValues(
      solution_ij, nrow_local, indices, solution.values_pt());
    solution.redistribute(soln_dist_pt);
    delete[] indices;
    delete soln_dist_pt;

    // output any error message
    if (Hypre_error_messages)
    {
      std::ostringstream message;
      int err = HypreHelpers::check_HYPRE_error_flag(message);
      if (err)
      {
        OomphLibWarning(message.str(),
                        "HypreSolver::hypre_solve()",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }

    // deallocation
    HYPRE_IJVectorDestroy(solution_ij);
    HYPRE_IJVectorDestroy(rhs_ij);

    // Record time
    double solve_time = 0;
    if (Output_info)
    {
      double t_end = TimingHelpers::timer();
      solve_time = t_end - t_start;
    }

    // output timings and info
    if (Output_info)
    {
      oomph_info << "Time for HYPRE solve [s] : " << solve_time << std::endl;
    }

    // for iterative solvers output iterations and final norm
    if ((Hypre_method >= CG) && (Hypre_method <= BoomerAMG))
    {
      if (iterations > 1)
      {
        if (Output_info)
          oomph_info << "Number of iterations         : " << iterations
                     << std::endl;
        if (Output_info)
          oomph_info << "Final Relative Residual Norm : " << norm << std::endl;
      }
    }
  }


  //===================================================================
  /// hypre_clean_up_memory() deletes any existing Hypre solver and
  /// Hypre matrix
  //====================================================================
  void HypreInterface::hypre_clean_up_memory()
  {
    // is there an existing solver
    if (Existing_solver != None)
    {
      // delete matrix
      HYPRE_IJMatrixDestroy(Matrix_ij);

      // delete solver
      if (Existing_solver == BoomerAMG)
      {
        HYPRE_BoomerAMGDestroy(Solver);
      }
      else if (Existing_solver == CG)
      {
        HYPRE_ParCSRPCGDestroy(Solver);
      }
      else if (Existing_solver == GMRES)
      {
        HYPRE_ParCSRGMRESDestroy(Solver);
      }
      else if (Existing_solver == BiCGStab)
      {
        HYPRE_ParCSRBiCGSTABDestroy(Solver);
      }
      else if (Existing_solver == Euclid)
      {
        HYPRE_EuclidDestroy(Solver);
      }
      else if (Existing_solver == ParaSails)
      {
        HYPRE_ParaSailsDestroy(Solver);
      }
      Existing_solver = None;

      // delete preconditioner
      if (Existing_preconditioner == BoomerAMG)
      {
        HYPRE_BoomerAMGDestroy(Preconditioner);
      }
      else if (Existing_preconditioner == Euclid)
      {
        HYPRE_EuclidDestroy(Preconditioner);
      }
      else if (Existing_preconditioner == ParaSails)
      {
        HYPRE_ParaSailsDestroy(Preconditioner);
      }
      Existing_preconditioner = None;

      // check error flag
      if (Hypre_error_messages)
      {
        std::ostringstream message;
        int err = HypreHelpers::check_HYPRE_error_flag(message);
        if (err)
        {
          OomphLibWarning(message.str(),
                          "HypreSolver::clean_up_memory()",
                          OOMPH_EXCEPTION_LOCATION);
        }
      }
    }
  }


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // functions for HypreSolver class
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //===================================================================
  /// Problem-based solve function to generate the Jacobian matrix and
  /// residual vector and use HypreInterface::hypre_solver_setup(...)
  /// and HypreInterface::hypre_solve(...) to solve the linear system.
  /// This function will delete any existing data.
  /// Note: the returned solution vector is NOT distributed, i.e. all
  /// processors hold all values because this is what the Newton solver
  /// requires.
  //====================================================================
  void HypreSolver::solve(Problem* const& problem_pt, DoubleVector& solution)
  {
    double t_start = TimingHelpers::timer();

    // Set Output_time flag for HypreInterface
    Output_info = Doc_time;

    // Delete any existing solver data
    clean_up_memory();

    // Set flag to allow deletion of the oomphlib Jacobian matrix
    // (we're in control)
    Delete_input_data = true;

    // Get oomph-lib Jacobian matrix and residual vector
    DoubleVector residual;
    CRDoubleMatrix* matrix = new CRDoubleMatrix;
    problem_pt->get_jacobian(residual, *matrix);

    // Output times
    if (Doc_time)
    {
      oomph_info << "Time to generate Jacobian and residual [s] : ";
      double t_end = TimingHelpers::timer();
      oomph_info << t_end - t_start << std::endl;
    }

    // generate hypre matrix
    hypre_matrix_setup(matrix);

    // call hypre_solver_setup to generate linear solver data
    hypre_solver_setup();

    // perform hypre_solve
    hypre_solve(residual, solution);

    // delete solver data if required
    if (!Enable_resolve)
    {
      clean_up_memory();
    }
  }


  //===================================================================
  /// Uses HypreInterface::hypre_solve(...) to solve the linear system
  /// for a CRDoubleMatrix or a DistributedCRDoubleMatrix.
  /// In the latter case, the rhs needs to be distributed too,
  /// i.e. the length of the rhs vector must be equal to the number of
  /// rows stored locally.
  /// Note: the returned solution vector is never distributed, i.e. all
  /// processors hold all values
  //====================================================================
  void HypreSolver::solve(DoubleMatrixBase* const& matrix_pt,
                          const DoubleVector& rhs,
                          DoubleVector& solution)
  {
#ifdef PARANOID
    // check the matrix is square
    if (matrix_pt->nrow() != matrix_pt->ncol())
    {
      std::ostringstream error_message;
      error_message << "HypreSolver require a square matrix. "
                    << "Matrix is " << matrix_pt->nrow() << " by "
                    << matrix_pt->ncol() << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Set Output_time flag for HypreInterface
    Output_info = Doc_time;

    // Clean up existing solver data
    clean_up_memory();

    // Set flag to decide if oomphlib matrix can be deleted
    // (Recall that Delete_matrix defaults to false).
    Delete_input_data = Delete_matrix;

    // Try cast to a CRDoubleMatrix
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);

    // If cast successful set things up for a serial solve
    if (cr_matrix_pt)
    {
      // rebuild the distribution
      this->build_distribution(cr_matrix_pt->distribution_pt());

#ifdef PARANOID
      // check that rhs has the same distribution as the matrix (now stored
      // in Distribution_pt)
      if (*this->distribution_pt() != *rhs.distribution_pt())
      {
        std::ostringstream error_message;
        error_message << "The distribution of the rhs vector and the matrix "
                      << " must be the same" << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      // if the solution is setup make sure it has the same distribution as
      // the matrix as well
      if (solution.built())
      {
        if (*this->distribution_pt() != *solution.distribution_pt())
        {
          std::ostringstream error_message;
          error_message << "The distribution of the solution vector is setup "
                        << "there it must be the same as the matrix."
                        << std::endl;
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif

      hypre_matrix_setup(cr_matrix_pt);
    }
    else
    {
#ifdef PARANOID
      std::ostringstream error_message;
      error_message << "HypreSolver only work with "
                    << "CRDoubleMatrix matrices" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    }

    // call hypre_setup to generate Hypre matrix and linear solver data
    hypre_solver_setup();

    // perform hypre_solve
    hypre_solve(rhs, solution);

    // delete solver data if required
    if (!Enable_resolve)
    {
      clean_up_memory();
    }
  }


  //===================================================================
  /// Resolve performs a linear solve using current solver data (if
  /// such data exists).
  //====================================================================
  void HypreSolver::resolve(const DoubleVector& rhs, DoubleVector& solution)
  {
#ifdef PARANOID
    // check solver data exists
    if (existing_solver() == None)
    {
      std::ostringstream error_message;
      error_message << "resolve(...) requires that solver data has been "
                    << "set up by a previous call to solve(...) after "
                    << "a call to enable_resolve()" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // check that rhs has the same distribution as the matrix (now stored
    // in Distribution_pt)
    if (*this->distribution_pt() != *rhs.distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the rhs vector and the matrix "
                    << " must be the same" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // if the solution is setup make sure it has the same distribution as
    // the matrix as well
    if (solution.built())
    {
      if (*this->distribution_pt() != *solution.distribution_pt())
      {
        std::ostringstream error_message;
        error_message << "The distribution of the solution vector is setup "
                      << "there it must be the same as the matrix."
                      << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Set Output_info flag for HypreInterface
    Output_info = Doc_time;

    // solve
    hypre_solve(rhs, solution);

    // Note: never delete solver data as the preconditioner is typically
    // called repeatedly.
  }


  //===================================================================
  /// clean_up_memory() deletes any existing solver data
  //====================================================================
  void HypreSolver::clean_up_memory()
  {
    hypre_clean_up_memory();
  }


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // functions for HyprePreconditioner class
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// Static double that accumulates the preconditioner
  /// solve time of all instantiations of this class. Reset
  /// it manually, e.g. after every Newton solve, using
  /// reset_cumulative_solve_times().
  //=============================================================================
  double HyprePreconditioner::Cumulative_preconditioner_solve_time = 0.0;

  //=============================================================================
  /// map of static doubles that accumulates the preconditioner
  /// solve time of all instantiations of this class, labeled by
  /// context string. Reset
  /// it manually, e.g. after every Newton solve, using
  /// reset_cumulative_solve_times().
  //=============================================================================
  std::map<std::string, double>
    HyprePreconditioner::Context_based_cumulative_solve_time;

  //=============================================================================
  /// Static unsigned that accumulates the number of preconditioner
  /// solves of all instantiations of this class. Reset
  /// it manually, e.g. after every Newton solve, using
  /// reset_cumulative_solve_times().
  //=============================================================================
  unsigned HyprePreconditioner::Cumulative_npreconditioner_solve = 0;

  //=============================================================================
  /// Static unsigned that accumulates the number of preconditioner
  /// solves of all instantiations of this class, labeled by
  /// context string. Reset
  /// it manually, e.g. after every Newton solve, using
  /// reset_cumulative_solve_times().
  //=============================================================================
  std::map<std::string, unsigned>
    HyprePreconditioner::Context_based_cumulative_npreconditioner_solve;

  //=============================================================================
  /// Static unsigned that stores nrow for the most recent
  /// instantiations of this class, labeled by
  /// context string. Reset
  /// it manually, e.g. after every Newton solve, using
  /// reset_cumulative_solve_times().
  //=============================================================================
  std::map<std::string, unsigned> HyprePreconditioner::Context_based_nrow;

  //=============================================================================
  /// Report cumulative solve times of all instantiations of this
  /// class
  //=============================================================================
  void HyprePreconditioner::report_cumulative_solve_times()
  {
    oomph_info << "\n\n=====================================================\n";
    oomph_info << "Cumulative HyprePreconditioner solve time "
               << HyprePreconditioner::Cumulative_preconditioner_solve_time
               << " for " << Cumulative_npreconditioner_solve << " solves";
    if (Cumulative_npreconditioner_solve != 0)
    {
      oomph_info << " ( "
                 << HyprePreconditioner::Cumulative_preconditioner_solve_time /
                      double(Cumulative_npreconditioner_solve)
                 << " per solve )";
    }
    oomph_info << std::endl << std::endl;
    if (Context_based_cumulative_solve_time.size() > 0)
    {
      oomph_info << "Breakdown by context: " << std::endl;
      for (std::map<std::string, double>::iterator it =
             Context_based_cumulative_solve_time.begin();
           it != Context_based_cumulative_solve_time.end();
           it++)
      {
        oomph_info
          << (*it).first << " " << (*it).second << " for "
          << Context_based_cumulative_npreconditioner_solve[(*it).first]
          << " solves";
        if (Context_based_cumulative_npreconditioner_solve[(*it).first] != 0)
        {
          oomph_info
            << " ( "
            << (*it).second /
                 double(
                   Context_based_cumulative_npreconditioner_solve[(*it).first])
            << " per solve; "
            << (*it).second /
                 double(
                   Context_based_cumulative_npreconditioner_solve[(*it)
                                                                    .first]) /
                 double(Context_based_nrow[(*it).first])
            << " per solve per dof )";
        }
        oomph_info << std::endl;
      }
    }
    oomph_info << "\n=====================================================\n";
    oomph_info << std::endl;
  }

  //=============================================================================
  /// Reset cumulative solve times
  //=============================================================================
  void HyprePreconditioner::reset_cumulative_solve_times()
  {
    Cumulative_preconditioner_solve_time = 0.0;
    Context_based_cumulative_solve_time.clear();
    Cumulative_npreconditioner_solve = 0;
    Context_based_cumulative_npreconditioner_solve.clear();
    Context_based_nrow.clear();
  }


  //=============================================================================
  /// An interface to allow HypreSolver to be used as a Preconditioner
  /// for the oomph-lib IterativeLinearSolver class.
  /// Matrix has to be of type CRDoubleMatrix or DistributedCRDoubleMatrix.
  //=============================================================================
  void HyprePreconditioner::setup()
  {
    // Set Output_info flag for HypreInterface
    Output_info = Doc_time;

#ifdef PARANOID
    // check the matrix is square
    if (matrix_pt()->nrow() != matrix_pt()->ncol())
    {
      std::ostringstream error_message;
      error_message << "HyprePreconditioner require a square matrix. "
                    << "Matrix is " << matrix_pt()->nrow() << " by "
                    << matrix_pt()->ncol() << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // clean up any previous solver data
    clean_up_memory();

    // set flag to decide if oomphlib matrix can be deleted
    // (Recall that Delete_matrix defaults to false).
    Delete_input_data = Delete_matrix;

    // Try casting to a CRDoubleMatrix
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());

    // If cast successful set things up for a serial solve
    if (cr_matrix_pt)
    {
      this->build_distribution(cr_matrix_pt->distribution_pt());
      hypre_matrix_setup(cr_matrix_pt);
    }
    else
    {
#ifdef PARANOID
      std::ostringstream error_message;
      error_message << "HyprePreconditioner only work with "
                    << "CRDoubleMatrix matrices" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    }

    if (Context_string != "")
    {
      oomph_info << "Setup of HyprePreconditioner in context \" "
                 << Context_string << "\": nrow, nrow_local, nnz "
                 << cr_matrix_pt->nrow() << " " << cr_matrix_pt->nrow_local()
                 << " " << cr_matrix_pt->nnz() << std::endl;
    }
    Context_based_nrow[Context_string] = cr_matrix_pt->nrow();

    // call hypre_solver_setup
    hypre_solver_setup();
  }

  //===================================================================
  /// Preconditioner_solve uses a hypre solver to precondition vector r
  //====================================================================
  void HyprePreconditioner::preconditioner_solve(const DoubleVector& r,
                                                 DoubleVector& z)
  {
    // Store time
    double t_start = TimingHelpers::timer();

#ifdef PARANOID
    // check solver data exists
    if (existing_solver() == None)
    {
      std::ostringstream error_message;
      error_message << "preconditioner_solve(...) requires that data has "
                    << "been set up using the function setup(...)" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // check that rhs has the same distribution as the matrix (now stored
    // in Distribution_pt)
    if (*this->distribution_pt() != *r.distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the rhs vector and the matrix "
                    << " must be the same" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // if the solution is setup make sure it has the same distribution as
    // the matrix as well
    if (z.built())
    {
      if (*this->distribution_pt() != *z.distribution_pt())
      {
        std::ostringstream error_message;
        error_message << "The distribution of the solution vector is setup "
                      << "there it must be the same as the matrix."
                      << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Switch off any timings for the solve
    Output_info = false;

    // perform hypre_solve
    hypre_solve(r, z);

    // Add to cumulative solve time
    double t_end = TimingHelpers::timer();
    Cumulative_preconditioner_solve_time += (t_end - t_start);
    Cumulative_npreconditioner_solve++;
    My_cumulative_preconditioner_solve_time += (t_end - t_start);
    if (Context_string != "")
    {
      Context_based_cumulative_solve_time[Context_string] += (t_end - t_start);
      Context_based_cumulative_npreconditioner_solve[Context_string]++;
    }
  }


  //===================================================================
  /// clean_up_memory() deletes any existing Hypre solver and
  /// Hypre matrix
  //====================================================================
  void HyprePreconditioner::clean_up_memory()
  {
    hypre_clean_up_memory();
  }

} // namespace oomph
