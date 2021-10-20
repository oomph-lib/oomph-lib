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
// Interface to MUMPS solver


#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif


#include <iostream>
#include <vector>


// oomph-lib headers
#include "mumps_solver.h"
#include "Vector.h"
#include "oomph_utilities.h"
#include "problem.h"


namespace oomph
{
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////


  //=========================================================================
  ///  Default factor for workspace -- static so it can be overwritten
  /// globally.
  //=========================================================================
  int MumpsSolver::Default_workspace_scaling_factor = 5;


  //=========================================================================
  /// Static warning to suppress warnings about incorrect distribution of
  /// RHS vector. Default is false
  //=========================================================================
  bool MumpsSolver::Suppress_incorrect_rhs_distribution_warning_in_resolve =
    false;

  //=============================================================================
  /// Constructor: Call setup
  //=============================================================================
  MumpsSolver::MumpsSolver()
    : Suppress_solve(false),
      Doc_stats(false),
      Suppress_warning_about_MPI_COMM_WORLD(false),
      Suppress_mumps_info_during_solve(false),
      Mumps_is_initialised(false),
      Workspace_scaling_factor(Default_workspace_scaling_factor),
      Delete_matrix_data(false),
      Mumps_struc_pt(nullptr),
      Jacobian_symmetry_flag(MumpsJacobianSymmetryFlags::Unsymmetric),
      Jacobian_ordering_flag(MumpsJacobianOrderingFlags::Metis_ordering)
  {
  }

  //=============================================================================
  /// Initialise instance of mumps data structure
  //=============================================================================
  void MumpsSolver::initialise_mumps()
  {
    // Make new instance of Mumps data structure
    Mumps_struc_pt = new DMUMPS_STRUC_C;

    // Mumps' hack to indicate that mpi_comm_world is used. Conversion of any
    // other communicators appears to be non-portable, so we simply
    // issue a warning if we later discover that oomph-lib's communicator
    // is not MPI_COMM_WORLD
    Mumps_struc_pt->comm_fortran = -987654;

    // Root processor participates in solution
    Mumps_struc_pt->par = 1;

    // Set matrix symmetry properties
    // (unsymmetric / symmetric positive definite / general symmetric)
    Mumps_struc_pt->sym = Jacobian_symmetry_flag;

    // First call does the initialise phase
    Mumps_struc_pt->job = -1;

    // Call c interface to double precision mumps to initialise
    dmumps_c(Mumps_struc_pt);
    Mumps_is_initialised = true;

    // Output stream for global info on host. Negative value suppresses printing
    Mumps_struc_pt->ICNTL(3) = -1;

    // Only show error messages and stats
    // (4 for full doc; creates huge amount of output)
    Mumps_struc_pt->ICNTL(4) = 2;

    // Write matrix
    // sprintf(Mumps_struc_pt->write_problem,"/work/e173/e173/mheil/matrix");

    // Assembled matrix (rather than element-by_element)
    Mumps_struc_pt->ICNTL(5) = 0;

    // set the package to use for ordering during (sequential) analysis
    Mumps_struc_pt->ICNTL(7) = Jacobian_ordering_flag;

    // Distributed problem with user-specified distribution
    Mumps_struc_pt->ICNTL(18) = 3;

    // Dense RHS
    Mumps_struc_pt->ICNTL(20) = 0;

    // Non-distributed solution
    Mumps_struc_pt->ICNTL(21) = 0;
  }


  //=============================================================================
  /// Shutdown mumps
  //=============================================================================
  void MumpsSolver::shutdown_mumps()
  {
    if (Mumps_is_initialised)
    {
      // Shut down flag
      Mumps_struc_pt->job = -2;

      // Call c interface to double precision mumps to shut down
      dmumps_c(Mumps_struc_pt);

      // Cleanup
      delete Mumps_struc_pt;
      Mumps_struc_pt = 0;

      Mumps_is_initialised = false;

      // Cleanup storage
      Irn_loc.clear();
      Jcn_loc.clear();
      A_loc.clear();
    }
  }


  //=============================================================================
  /// Destructor: Shutdown mumps
  //=============================================================================
  MumpsSolver::~MumpsSolver()
  {
    shutdown_mumps();
  }


  //=============================================================================
  /// LU decompose the matrix addressed by matrix_pt using
  /// mumps. The resulting matrix factors are stored internally.
  /// Note: if Delete_matrix_data is true the function
  /// matrix_pt->clean_up_memory() will be used to wipe the matrix data.
  //=============================================================================
  void MumpsSolver::factorise(DoubleMatrixBase* const& matrix_pt)
  {
    // Initialise timer
    double t_start = TimingHelpers::timer();

    // set the distribution
    DistributableLinearAlgebraObject* dist_matrix_pt =
      dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt);
    if (dist_matrix_pt)
    {
      // the solver has the same distribution as the matrix if possible
      this->build_distribution(dist_matrix_pt->distribution_pt());
    }
    else
    {
      throw OomphLibError("Matrix must be distributable",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }


    // Check that we have a square matrix
#ifdef PARANOID
    int n = matrix_pt->nrow();
    int m = matrix_pt->ncol();
    if (n != m)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "Can only solve for square matrices\n"
                           << "N, M " << n << " " << m << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    if (!Suppress_warning_about_MPI_COMM_WORLD)
    {
      if (this->distribution_pt()->communicator_pt()->mpi_comm() !=
          MPI_COMM_WORLD)
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "Warning: Mumps wrapper assumes that communicator is "
             "MPI_COMM_WORLD\n"
          << "         which is not the case, so mumps may die...\n"
          << "         If it does initialise oomph-lib's mpi without "
             "requesting\n"
          << "         the communicator to be a duplicate of MPI_COMM_WORLD\n"
          << "         (done via an optional boolean to "
             "MPI_Helpers::init(...)\n\n"
          << "         (You can suppress this warning by recompiling without\n"
          << "         paranoia or calling \n\n"
          << "         "
             "MumpsSolver::enable_suppress_warning_about_MPI_COMM_WORLD()\n\n"
          << "         \n";
        OomphLibWarning(error_message_stream.str(),
                        "MumpsSolver::factorise()",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }

#endif

    // Initialise
    if (Mumps_is_initialised)
    {
      shutdown_mumps();
    }
    initialise_mumps();


    // Doc stats?
    if (Doc_stats)
    {
      // Output stream for global info on host. Negative value suppresses
      // printing
      Mumps_struc_pt->ICNTL(3) = 6;
    }

    // number of processors
    unsigned nproc = 1;
    if (dist_matrix_pt != 0)
    {
      nproc = dist_matrix_pt->distribution_pt()->communicator_pt()->nproc();
    }

    // Is it a CRDoubleMatrix?
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);
    if (cr_matrix_pt != 0)
    {
#ifdef PARANOID
      // paranoid check that the matrix has been set up
      if (!cr_matrix_pt->built())
      {
        throw OomphLibError(
          "To apply MumpsSolver to a CRDoubleMatrix - it must be built",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // if the matrix is distributed then set up solver
      if ((nproc == 1) || (cr_matrix_pt->distributed()))
      {
        double t_start_copy = TimingHelpers::timer();

        // Find the number of rows and non-zero entries in the matrix
        const int nnz_loc = int(cr_matrix_pt->nnz());
        const int n = matrix_pt->nrow();

        // Create mumps storage

        // Create vector for row numbers
        Irn_loc.resize(nnz_loc);

        // Create vector for column numbers
        Jcn_loc.resize(nnz_loc);

        // Vector for entries
        A_loc.resize(nnz_loc);

        // First row held on this processor
        int first_row = cr_matrix_pt->first_row();

        // Copy into coordinate storage scheme using pointer arithmetic
        double* matrix_value_pt = cr_matrix_pt->value();
        int* matrix_index_pt = cr_matrix_pt->column_index();
        int* matrix_start_pt = cr_matrix_pt->row_start();
        int i_row = 0;

        // is the matrix symmetric? If so, we must only provide
        // MUMPS with the upper (or lower) triangular part of the Jacobian
        if (Mumps_struc_pt->sym != MumpsJacobianSymmetryFlags::Unsymmetric)
        {
          oomph_info << "Assuming Jacobian matrix is symmetric "
                     << "(passing MUMPS the upper triangular part)"
                     << std::endl;
        }

        for (int count = 0; count < nnz_loc; count++)
        {
          A_loc[count] = matrix_value_pt[count];
          Jcn_loc[count] = matrix_index_pt[count] + 1;
          if (count < matrix_start_pt[i_row + 1])
          {
            Irn_loc[count] = first_row + i_row + 1;
          }
          else
          {
            i_row++;
            Irn_loc[count] = first_row + i_row + 1;
          }

          // only pass the upper triangular part (and diagonal)
          // if we have a symmetric matrix (MUMPS sums values,
          // so need to set the lwoer triangle to zero)
          if ((Mumps_struc_pt->sym !=
               MumpsJacobianSymmetryFlags::Unsymmetric) &&
              (Irn_loc[count] > Jcn_loc[count]))
          {
            A_loc[count] = 0.0;
          }
        }

        // Now delete the matrix if we are allowed
        if (Delete_matrix_data == true)
        {
          cr_matrix_pt->clear();
        }

        if ((Doc_time) &&
            (this->distribution_pt()->communicator_pt()->my_rank() == 0))
        {
          double t_end_copy = TimingHelpers::timer();
          oomph_info << "Time for copying matrix into MumpsSolver data "
                        "structure [sec]       : "
                     << t_end_copy - t_start_copy << std::endl;
        }

        // Call mumps factorisation
        //-------------------------

        // Specify size of system
        Mumps_struc_pt->n = n;

        // Number of nonzeroes
        Mumps_struc_pt->nz_loc = nnz_loc;

        // The entries
        Mumps_struc_pt->irn_loc = &Irn_loc[0];
        Mumps_struc_pt->jcn_loc = &Jcn_loc[0];
        Mumps_struc_pt->a_loc = &A_loc[0];

        double t_start_analyse = TimingHelpers::timer();

        // Do analysis
        Mumps_struc_pt->job = 1;
        dmumps_c(Mumps_struc_pt);


        if ((Doc_time) &&
            (this->distribution_pt()->communicator_pt()->my_rank() == 0))
        {
          double t_end_analyse = TimingHelpers::timer();
          oomph_info
            << "Time for mumps analysis stage in MumpsSolver [sec]       : "
            << t_end_analyse - t_start_analyse
            << "\n(Ordering generated using: ";

          switch (Mumps_struc_pt->INFOG(7))
          {
            case MumpsJacobianOrderingFlags::Scotch_ordering:
              oomph_info << "SCOTCH";
              break;
            case MumpsJacobianOrderingFlags::Pord_ordering:
              oomph_info << "PORD";
              break;
            case MumpsJacobianOrderingFlags::Metis_ordering:
              oomph_info << "METIS";
              break;
            default:
              oomph_info << Mumps_struc_pt->INFOG(7);
          }

          oomph_info << ")" << std::endl;
        }


        int my_rank = this->distribution_pt()->communicator_pt()->my_rank();

        // Document estimate for size of workspace
        if (my_rank == 0)
        {
          if (!Suppress_mumps_info_during_solve)
          {
            oomph_info << "Estimated max. workspace in MB: "
                       << Mumps_struc_pt->INFOG(26) << std::endl;
          }
        }

        double t_start_factor = TimingHelpers::timer();

        // Loop until successfully factorised
        bool factorised = false;
        while (!factorised)
        {
          // Set workspace to multiple of that -- ought to be "significantly
          // larger than infog(26)", according to the manual :(
          Mumps_struc_pt->ICNTL(23) =
            Workspace_scaling_factor * Mumps_struc_pt->INFOG(26);

          if (!Suppress_mumps_info_during_solve)
          {
            oomph_info << "Attempting factorisation with workspace estimate: "
                       << Mumps_struc_pt->ICNTL(23) << " MB\n";
            oomph_info << "corresponding to Workspace_scaling_factor = "
                       << Workspace_scaling_factor << "\n";
          }

          // Do factorisation
          Mumps_struc_pt->job = 2;
          dmumps_c(Mumps_struc_pt);

          // Check for error
          if (Mumps_struc_pt->INFOG(1) != 0)
          {
            if (!Suppress_mumps_info_during_solve)
            {
              oomph_info << "Error during mumps factorisation!\n";
              oomph_info << "Error codes: " << Mumps_struc_pt->INFO(1) << " "
                         << Mumps_struc_pt->INFO(2) << std::endl;
            }

            // Increase scaling factor for workspace and run again
            Workspace_scaling_factor *= 2;

            if (!Suppress_mumps_info_during_solve)
            {
              oomph_info << "Increasing workspace_scaling_factor to "
                         << Workspace_scaling_factor << std::endl;
            }
          }
          else
          {
            factorised = true;
          }
        }


        if ((Doc_time) &&
            (this->distribution_pt()->communicator_pt()->my_rank() == 0))
        {
          double t_end_factor = TimingHelpers::timer();
          oomph_info << "Time for actual mumps factorisation in MumpsSolver "
                        "[sec]       : "
                     << t_end_factor - t_start_factor << std::endl;
        }
      }
      // else the CRDoubleMatrix is not distributed
      else
      {
        std::ostringstream error_message_stream;
        error_message_stream << "MumpsSolver only works for a "
                             << " distributed CRDoubleMatrix\n";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
    // Otherwise throw an error
    else
    {
      std::ostringstream error_message_stream;
      error_message_stream << "MumpsSolver implemented only for "
                           << "distributed CRDoubleMatrix. \n";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }


    if ((Doc_time) &&
        (this->distribution_pt()->communicator_pt()->my_rank() == 0))
    {
      double t_end = TimingHelpers::timer();
      oomph_info << "Time for MumpsSolver factorisation [sec]       : "
                 << t_end - t_start << std::endl;
    }

    // Switch off docing again by setting output stream for global info on
    // to negative number
    Mumps_struc_pt->ICNTL(3) = -1;
  }

  //=============================================================================
  /// Do the backsubstitution for mumps solver.
  /// This does not make any assumption about the distribution of the
  /// vectors
  //=============================================================================
  void MumpsSolver::backsub(const DoubleVector& rhs, DoubleVector& result)
  {
    double t_start = TimingHelpers::timer();

#ifdef PARANOID
    if (!Suppress_warning_about_MPI_COMM_WORLD)
    {
      if (this->distribution_pt()->communicator_pt()->mpi_comm() !=
          MPI_COMM_WORLD)
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "Warning: Mumps wrapper assumes that communicator is "
             "MPI_COMM_WORLD\n"
          << "         which is not the case, so mumps may die...\n"
          << "         If it does initialise oomph-lib's mpi without "
             "requesting\n"
          << "         the communicator to be a duplicate of MPI_COMM_WORLD\n"
          << "         (done via an optional boolean to "
             "MPI_Helpers::init(...)\n\n"
          << "         (You can suppress this warning by recompiling without\n"
          << "         paranoia or calling \n\n"
          << "         "
             "MumpsSolver::enable_suppress_warning_about_MPI_COMM_WORLD()\n\n"
          << "         \n";
        OomphLibWarning(error_message_stream.str(),
                        "MumpsSolver::backsub()",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Initialised?
    if (!Mumps_is_initialised)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "Mumps has not been initialised.";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // check that the rhs vector is setup
    if (!rhs.distribution_pt()->built())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The vectors rhs must be setup";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

#endif

    // Check that the rhs distribution is the same as the distribution as this
    // solver. If not redistribute and issue a warning
    LinearAlgebraDistribution rhs_distribution(rhs.distribution_pt());
    if (!(*rhs.distribution_pt() == *this->distribution_pt()))
    {
      if (!Suppress_incorrect_rhs_distribution_warning_in_resolve)
      {
        std::ostringstream warning_stream;
        warning_stream << "The distribution of rhs vector does not match that "
                          "of the solver.\n";
        warning_stream
          << "The rhs may have to be redistributed but we're not doing this "
             "because\n"
          << "I'm no longer convinced it's necessary. Keep an eye on this...\n";
        warning_stream
          << "To remove this warning you can either:\n"
          << "    i) Ensure that the rhs vector has the correct distribution\n"
          << "       before calling the resolve() function\n"
          << "or ii) Set the flag \n"
          << " MumpsSolver::Suppress_incorrect_rhs_distribution_warning_in_"
             "resolve\n"
          << "       to be true\n\n";

        OomphLibWarning(warning_stream.str(),
                        "MumpsSolver::resolve()",
                        OOMPH_EXCEPTION_LOCATION);
      }

      // //Have to cast away const-ness (which tells us that we shouldn't really
      // //be doing this!)
      // const_cast<DoubleVector&>(rhs).redistribute(this->distribution_pt());
    }


#ifdef PARANOID
    // if the result vector is setup then check it has the same distribution
    // as the rhs
    if (result.distribution_built())
    {
      if (!(*result.distribution_pt() == *rhs.distribution_pt()))
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The result vector distribution has been setup; it must have the "
          << "same distribution as the rhs vector.";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif


    // Doc stats?
    if (Doc_stats)
    {
      // Output stream for global info on host. Negative value suppresses
      // printing
      Mumps_struc_pt->ICNTL(3) = 6;
    }

    // number of DOFs
    int ndof = this->distribution_pt()->nrow();

    // Make backup to avoid over-writing
    DoubleVector tmp_rhs;
    tmp_rhs = rhs;

    // Now turn this into a global (non-distributed) vector
    // because that's what mumps needs

    // Make a global distribution (i.e. one that isn't distributed)
    LinearAlgebraDistribution global_distribution(
      this->distribution_pt()->communicator_pt(), ndof, false);

    // Redistribute the tmp_rhs vector with this distribution -- it's
    // now "global", as required for mumps
    tmp_rhs.redistribute(&global_distribution);

    // Do the backsubsitution phase -- overwrites the tmp_rhs vector with the
    // solution
    Mumps_struc_pt->rhs = &tmp_rhs[0];
    Mumps_struc_pt->job = 3;
    dmumps_c(Mumps_struc_pt);

    // Copy back
    unsigned n = Mumps_struc_pt->n;
    for (unsigned j = 0; j < n; j++)
    {
      tmp_rhs[j] = Mumps_struc_pt->rhs[j];
    }

    // Broadcast the result which is only held on root
    MPI_Bcast(&tmp_rhs[0],
              ndof,
              MPI_DOUBLE,
              0,
              this->distribution_pt()->communicator_pt()->mpi_comm());

    // If the result vector is distributed, re-distribute the
    // non-distributed tmp_rhs vector to match
    if (result.built())
    {
      tmp_rhs.redistribute(result.distribution_pt());
    }
    else
    {
      tmp_rhs.redistribute(this->distribution_pt());
    }

    // Now copy the tmp_rhs vector into the (matching) result
    result = tmp_rhs;

    if ((Doc_time) &&
        (this->distribution_pt()->communicator_pt()->my_rank() == 0))
    {
      double t_end = TimingHelpers::timer();
      oomph_info << "Time for MumpsSolver backsub [sec]       : "
                 << t_end - t_start << std::endl;
    }

    // Switch off docing again by setting output stream for global info on
    // to negative number
    Mumps_struc_pt->ICNTL(3) = -1;
  }

  //=========================================================================
  /// Clean up the memory allocated for the MumpsSolver solver
  //=========================================================================
  void MumpsSolver::clean_up_memory()
  {
    shutdown_mumps();
  }


  //=========================================================================
  /// Linear-algebra-type solver: Takes pointer to a matrix and rhs
  /// vector and returns the solution of the linear system. Problem pointer
  /// defaults to NULL and can be omitted. The function returns the global
  /// result vector. Matrix must be CRDoubleMatrix.
  /// Note: if Delete_matrix_data is true the function
  /// matrix_pt->clean_up_memory() will be used to wipe the matrix data.
  //=========================================================================
  void MumpsSolver::solve(DoubleMatrixBase* const& matrix_pt,
                          const DoubleVector& rhs,
                          DoubleVector& result)
  {
#ifdef PARANOID
    if (!Suppress_warning_about_MPI_COMM_WORLD)
    {
      if (dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt)
            ->distribution_pt()
            ->communicator_pt()
            ->mpi_comm() != MPI_COMM_WORLD)
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "Warning: Mumps wrapper assumes that communicator is "
             "MPI_COMM_WORLD\n"
          << "         which is not the case, so mumps may die...\n"
          << "         If it does initialise oomph-lib's mpi without "
             "requesting\n"
          << "         the communicator to be a duplicate of MPI_COMM_WORLD\n"
          << "         (done via an optional boolean to "
             "MPI_Helpers::init(...)\n\n"
          << "         (You can suppress this warning by recompiling without\n"
          << "         paranoia or calling \n\n"
          << "         "
             "MumpsSolver::enable_suppress_warning_about_MPI_COMM_WORLD()\n\n"
          << "         \n";
        OomphLibWarning(error_message_stream.str(),
                        "MumpsSolver::solve()",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Initialise timer
    double t_start = TimingHelpers::timer();

#ifdef PARANOID
    // check that the rhs vector is setup
    if (!rhs.distribution_pt()->built())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The vectors rhs must be setup";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // check that the matrix is square
    if (matrix_pt->nrow() != matrix_pt->ncol())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The matrix at matrix_pt must be square.";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // check that the matrix and the rhs vector have the same nrow()
    if (matrix_pt->nrow() != rhs.nrow())
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "The matrix and the rhs vector must have the same number of rows.";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }


    // if the matrix is distributable then should have the same distribution
    // as the rhs vector
    DistributableLinearAlgebraObject* ddist_matrix_pt =
      dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt);
    if (ddist_matrix_pt != 0)
    {
      if (!(*ddist_matrix_pt->distribution_pt() == *rhs.distribution_pt()))
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The matrix matrix_pt must have the same distribution as the "
          << "rhs vector.";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
    // if the matrix is not distributable then it the rhs vector should not be
    // distributed
    else
    {
      if (rhs.distribution_pt()->distributed())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The matrix (matrix_pt) is not distributable and therefore the rhs"
          << " vector must not be distributed";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
    // if the result vector is setup then check it has the same distribution
    // as the rhs
    if (result.built())
    {
      if (!(*result.distribution_pt() == *rhs.distribution_pt()))
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The result vector distribution has been setup; it must have the "
          << "same distribution as the rhs vector.";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }

#endif


    // set the distribution
    DistributableLinearAlgebraObject* dist_matrix_pt =
      dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt);
    if (dist_matrix_pt)
    {
      // the solver has the same distribution as the matrix if possible
      this->build_distribution(dist_matrix_pt->distribution_pt());
    }
    else
    {
      // the solver has the same distribution as the RHS
      this->build_distribution(rhs.distribution_pt());
    }

    // Factorise the matrix
    factorise(matrix_pt);

    // Now do the back solve
    backsub(rhs, result);

    // Doc time for solve
    double t_end = TimingHelpers::timer();
    Solution_time = t_end - t_start;

    if ((Doc_time) &&
        (this->distribution_pt()->communicator_pt()->my_rank() == 0))
    {
      oomph_info << "Time for MumpsSolver solve [sec]       : "
                 << t_end - t_start << std::endl;
    }

    // Switch off docing again by setting output stream for global info on
    // to negative number
    Mumps_struc_pt->ICNTL(3) = -1;

    // If we are not storing the solver data for resolves, delete it
    if (!Enable_resolve)
    {
      clean_up_memory();
    }
  }

  //==================================================================
  /// Solver: Takes pointer to problem and returns the results Vector
  /// which contains the solution of the linear system defined by
  /// the problem's fully assembled Jacobian and residual Vector.
  //==================================================================
  void MumpsSolver::solve(Problem* const& problem_pt, DoubleVector& result)
  {
#ifdef PARANOID
    if (!Suppress_warning_about_MPI_COMM_WORLD)
    {
      if (problem_pt->communicator_pt()->mpi_comm() != MPI_COMM_WORLD)
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "Warning: Mumps wrapper assumes that communicator is "
             "MPI_COMM_WORLD\n"
          << "         which is not the case, so mumps may die...\n"
          << "         If it does initialise oomph-lib's mpi without "
             "requesting\n"
          << "         the communicator to be a duplicate of MPI_COMM_WORLD\n"
          << "         (done via an optional boolean to "
             "MPI_Helpers::init(...)\n\n"
          << "         (You can suppress this warning by recompiling without\n"
          << "         paranoia or calling \n\n"
          << "         "
             "MumpsSolver::enable_suppress_warning_about_MPI_COMM_WORLD()\n\n"
          << "         \n";
        OomphLibWarning(error_message_stream.str(),
                        "MumpsSolver::solve()",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Initialise timer
    double t_start = TimingHelpers::timer();

    // number of dofs
    unsigned n_dof = problem_pt->ndof();

    // Set the distribution for the solver.
    this->distribution_pt()->build(problem_pt->communicator_pt(), n_dof);

    // Take a copy of Delete_matrix_data
    bool copy_of_Delete_matrix_data = Delete_matrix_data;

    // Set Delete_matrix to true
    Delete_matrix_data = true;

    // Initialise timer
    t_start = TimingHelpers::timer();

    // Storage for the distributed residuals vector
    DoubleVector residuals(this->distribution_pt(), 0.0);

    // Get the sparse jacobian and residuals of the problem
    CRDoubleMatrix jacobian(this->distribution_pt());
    problem_pt->get_jacobian(residuals, jacobian);

    // Doc time for setup
    double t_end = TimingHelpers::timer();
    Jacobian_setup_time = t_end - t_start;
    int my_rank = this->distribution_pt()->communicator_pt()->my_rank();
    if ((Doc_time) && (my_rank == 0))
    {
      oomph_info << "Time to set up CRDoubleMatrix Jacobian [sec]        : "
                 << Jacobian_setup_time << std::endl;
    }


    // Now call the linear algebra solve, if desired
    if (!Suppress_solve)
    {
      // If the distribution of the result has been build and
      // does not match that of
      // the solver then redistribute before the solve and return
      // to the incoming distribution afterwards.
      if ((result.built()) &&
          (!(*result.distribution_pt() == *this->distribution_pt())))
      {
        LinearAlgebraDistribution temp_global_dist(result.distribution_pt());
        result.build(this->distribution_pt(), 0.0);
        solve(&jacobian, residuals, result);
        result.redistribute(&temp_global_dist);
      }
      else
      {
        solve(&jacobian, residuals, result);
      }
    }

    // Set Delete_matrix back to original value
    Delete_matrix_data = copy_of_Delete_matrix_data;

    // Finalise/doc timings
    if ((Doc_time) && (my_rank == 0))
    {
      double t_end = TimingHelpers::timer();
      oomph_info << "Total time for MumpsSolver "
                 << "(np="
                 << this->distribution_pt()->communicator_pt()->nproc()
                 << ",N=" << problem_pt->ndof()
                 << ") [sec] : " << t_end - t_start << std::endl;
    }
  }


  //===============================================================
  /// Resolve the system defined by the last assembled jacobian
  /// and the specified rhs vector if resolve has been enabled.
  /// Note: returns the global result Vector.
  //===============================================================
  void MumpsSolver::resolve(const DoubleVector& rhs, DoubleVector& result)
  {
#ifdef PARANOID
    if (!Suppress_warning_about_MPI_COMM_WORLD)
    {
      if (this->distribution_pt()->communicator_pt()->mpi_comm() !=
          MPI_COMM_WORLD)
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "Warning: Mumps wrapper assumes communicator is MPI_COMM_WORLD\n"
          << "         which is not the case, so mumps may die...\n"
          << "         If it does initialise oomph-lib's mpi without "
             "requesting\n"
          << "         the communicator to be a duplicate of MPI_COMM_WORLD\n"
          << "         (done via an optional boolean to "
             "MPI_Helpers::init(...)\n\n"
          << "         (You can suppress this warning by recompiling without\n"
          << "         paranoia or calling\n\n"
          << "         "
             "MumpsSolver::enable_suppress_warning_about_MPI_COMM_WORLD()\n\n"
          << "         \n";
        OomphLibWarning(error_message_stream.str(),
                        "MumpsSolver::resolve()",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Store starting time for solve
    double t_start = TimingHelpers::timer();

    // Doc stats?
    if (Doc_stats)
    {
      // Output stream for global info on host. Negative value suppresses
      // printing
      Mumps_struc_pt->ICNTL(3) = 6;
    }

    // Now do the back substitution phase
    backsub(rhs, result);

    // Doc time for solve
    double t_end = TimingHelpers::timer();
    Solution_time = t_end - t_start;

    // Switch off docing again by setting output stream for global info on
    // to negative number
    Mumps_struc_pt->ICNTL(3) = -1;

    if ((Doc_time) &&
        (this->distribution_pt()->communicator_pt()->my_rank() == 0))
    {
      oomph_info << "Time for MumpsSolver solve [sec]: " << t_end - t_start
                 << std::endl;
    }
  }


} // namespace oomph
