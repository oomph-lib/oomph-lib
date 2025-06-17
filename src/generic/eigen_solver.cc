// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
// Non-inline functions for the eigensolvers

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif


// Include cfortran.h and the header for the LAPACK QZ routines
#include "cfortran.h"
#include "lapack_qz.h"

// Oomph-lib headers
#include "eigen_solver.h"
#include "linear_solver.h"
#include "problem.h"


namespace oomph
{
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////


  //==========================================================================
  /// Use LAPACK QZ to solve the real eigenproblem that is assembled by elements
  /// in a mesh in a Problem object. Note that the assembled matrices include
  /// the shift and are real. The eigenvalues and eigenvectors are, in general,
  /// complex.
  /// Eigenvalues may be infinite and are therefore returned as
  /// \f$ \lambda_i = \alpha_i / \beta_i \f$ where \f$ \alpha_i \f$ is complex
  /// while \f$ \beta_i \f$ is real. The actual eigenvalues may then be
  /// computed by doing the division, checking for zero betas to avoid NaNs.
  /// This is actually a helper function that stores re & imag parts of
  /// eigenvectors in a collection of real vectors; they
  /// are disentangled in the alternative version of this function that returns
  /// Vectors of complex vectors.
  /// At least n_eval eigenvalues are computed.
  //==========================================================================
  void LAPACK_QZ::solve_eigenproblem_helper(
    Problem* const& problem_pt,
    const int& n_eval,
    Vector<std::complex<double>>& alpha_eval,
    Vector<double>& beta_eval,
    Vector<DoubleVector>& eigenvector_aux)
  {
    // Some character identifiers for use in the LAPACK routine
    // Do not calculate the left eigenvectors
    char no_eigvecs[2] = "N";

    // Do caculate the eigenvectors
    char eigvecs[2] = "V";

    // Get the dimension of the matrix
    int n = problem_pt->ndof(); // Total size of matrix

    // Use padding?
    bool use_padding = false;

    // If the dimension of the matrix is even, then pad the arrays to
    // make the size odd. This somehow sorts out a strange run-time behaviour
    // identified by Rich Hewitt.
    // Actual size of matrix that will be allocated
    int padded_n = n;
    if (n % 2 == 0)
    {
      use_padding = true;
      padded_n += 1;
    }

    // Storage for the matrices in the packed form required by the LAPACK
    // routine
    double* M = new double[padded_n * padded_n];
    double* A = new double[padded_n * padded_n];

    // TEMPORARY
    // only use non-distributed matrices and vectors
    LinearAlgebraDistribution dist(problem_pt->communicator_pt(), n, false);
    this->build_distribution(dist);

    // Enclose in a separate scope so that memory is cleaned after assembly
    {
      // Allocated Row compressed matrices for the mass matrix and shifted main
      // matrix
      CRDoubleMatrix temp_M(this->distribution_pt()),
        temp_AsigmaM(this->distribution_pt());

      // Assemble the matrices; pass the shift into the assembly
      problem_pt->get_eigenproblem_matrices(temp_M, temp_AsigmaM, Sigma_real);

      // Now convert these matrices into the appropriate packed form
      unsigned index = 0;
      for (int i = 0; i < n; ++i)
      {
        for (int j = 0; j < n; ++j)
        {
          M[index] = temp_M(j, i);
          A[index] = temp_AsigmaM(j, i);
          ++index;
        }
        // If necessary pad the columns with zeros
        if (use_padding)
        {
          M[index] = 0.0;
          A[index] = 0.0;
          ++index;
        }
      }
      // No need to pad the final row because it is never accessed by the
      // routine.
    }

    // Temporary eigenvalue storage
    double* alpha_r = new double[n];
    double* alpha_i = new double[n];
    double* beta = new double[n];
    // Temporary eigenvector storage
    double* vec_left = new double[1];
    const int leading_dimension_vec_left = 1;
    double* vec_right = new double[n * n];
    const int leading_dimension_vec_right = n;

    // Workspace for the LAPACK routine
    std::vector<double> work(1, 0.0);
    // The dimension of the workspace array. Set to -1 so that a workspace
    // query is assumed and LAPACK calculates the optimum size which is
    // returned as the first value of work.
    const int query_workspace = -1;
    // Info flag for the LAPACK routine
    int info = 0;

    // Call FORTRAN LAPACK to get the required workspace
    // Note the use of the padding leading dimension for the matrices
    LAPACK_DGGEV(no_eigvecs,
                 eigvecs,
                 n,
                 &A[0],
                 padded_n,
                 &M[0],
                 padded_n,
                 alpha_r,
                 alpha_i,
                 beta,
                 vec_left,
                 leading_dimension_vec_left,
                 vec_right,
                 leading_dimension_vec_right,
                 &work[0],
                 query_workspace,
                 info);

    // Succesful completion?
    if (info != 0)
    {
      DGGEV_error(info, n);
    }

    // Get the amount of requires workspace
    int required_workspace = (int)work[0];
    // If we need it
    work.resize(required_workspace);

    // call FORTRAN LAPACK again with the optimum workspace
    LAPACK_DGGEV(no_eigvecs,
                 eigvecs,
                 n,
                 &A[0],
                 padded_n,
                 &M[0],
                 padded_n,
                 alpha_r,
                 alpha_i,
                 beta,
                 vec_left,
                 leading_dimension_vec_left,
                 vec_right,
                 leading_dimension_vec_right,
                 &work[0],
                 required_workspace,
                 info);

    // Succesful completion?
    if (info != 0)
    {
      DGGEV_error(info, n);
    }

    // Now resize storage for the eigenvalues and eigenvectors
    // We get them all!
    alpha_eval.resize(n);
    beta_eval.resize(n);
    eigenvector_aux.resize(n);

    // Now convert the output into our format
    for (int i = 0; i < n; ++i)
    {
      // Encode eigenvalues in alpha and beta vectors. Shift is undone here
      alpha_eval[i] = std::complex<double>(Sigma_real + alpha_r[i], alpha_i[i]);
      beta_eval[i] = beta[i];

      // Resize the eigenvector  storage
      eigenvector_aux[i].build(this->distribution_pt(), 0.0);

      // Load up the mangeled storage for the eigenvector. What's called
      // eigenvector here isn't necessarily an eigenvector but provides storage
      // for the real and imaginary parts of the real or complex conjugate
      // eigenvectors. They are translated into actual complex eigenvectors
      // (at the cost of some additional storage and the gain of considerable
      // clarity) in the public version of this function.
      for (int k = 0; k < n; ++k)
      {
        eigenvector_aux[i][k] = vec_right[i * n + k];
      }
    }

    // Delete all allocated storage
    delete[] vec_right;
    delete[] vec_left;
    delete[] beta;
    delete[] alpha_r;
    delete[] alpha_i;
    delete[] A;
    delete[] M;
  }


  //==========================================================================
  /// Solve the real eigenproblem that is assembled by elements in
  /// a mesh in a Problem object. Note that the assembled matrices include the
  /// shift and are real. The eigenvalues and eigenvectors are,
  /// in general, complex. Eigenvalues may be infinite and are therefore
  /// returned as
  /// \f$ \lambda_i = \alpha_i / \beta_i \f$ where \f$ \alpha_i \f$ is complex
  /// while \f$ \beta_i \f$ is real. The actual eigenvalues may then be
  /// computed by doing the division, checking for zero betas to avoid NaNs.
  /// There's a convenience wrapper to this function that simply computes
  /// these eigenvalues regardless. That version may die in NaN checking is
  /// enabled (via the fenv.h header and the associated feenable function).
  /// At least n_eval eigenvalues are computed.
  //==========================================================================
  void LAPACK_QZ::solve_eigenproblem(Problem* const& problem_pt,
                                     const int& n_eval,
                                     Vector<std::complex<double>>& alpha_eval,
                                     Vector<double>& beta_eval,
                                     Vector<DoubleVector>& eigenvector_real,
                                     Vector<DoubleVector>& eigenvector_imag,
                                     const bool& do_adjoint_problem)
  {
    if (do_adjoint_problem)
    {
      throw OomphLibError("Solving an adjoint eigenproblem is not currently "
                          "implemented for LAPACK_QZ.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    Vector<DoubleVector> eigenvector_aux;

    // Call raw interface to lapack qz
    solve_eigenproblem_helper(
      problem_pt, n_eval, alpha_eval, beta_eval, eigenvector_aux);

    // Now move auxiliary data into actual complex eigenvalues
    // Instrutions from lapack qz (where VR translates into eigenvector_aux):
    //        VR is DOUBLE PRECISION array, dimension (LDVR,N)
    //        If JOBVR = 'V', the right eigenvectors v(j) are stored one
    //        after another in the columns of VR, in the same order as
    //        their eigenvalues. If the j-th eigenvalue is real, then
    //        v(j) = VR(:,j), the j-th column of VR. If the j-th and
    //        (j+1)-th eigenvalues form a complex conjugate pair, then
    //        v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
    //        Each eigenvector is scaled so the largest component has
    //        abs(real part)+abs(imag. part)=1.
    unsigned n = problem_pt->ndof();
    eigenvector_real.resize(n);
    eigenvector_imag.resize(n);
    unsigned eval_count = 0;
    while (eval_count < n)
    {
      // i-th eigenvalue is real:
      if (alpha_eval[eval_count].imag() == 0.0)
      {
        // Resize the single eigenvector associated with this
        // single real eigenvalue
        eigenvector_real[eval_count].build(this->distribution_pt(), 0.0);
        eigenvector_imag[eval_count].build(this->distribution_pt(), 0.0);
        for (unsigned j = 0; j < n; ++j)
        {
          eigenvector_real[eval_count][j] = eigenvector_aux[eval_count][j];
          eigenvector_imag[eval_count][j] = 0.0;
        }
        eval_count++;
      }
      // Assume (and check!) that complex conjugate pairs follow each other
      // as implied by
      // http://www.netlib.org/lapack/explore-html/d9/d8e/group__double_g_eeigen_ga4f59e87e670a755b41cbdd7e97f36bea.html
      else
      {
#ifdef PARANOID

        // Are the eigenvalues finite? If not skip test.
        if ((beta_eval[eval_count] != 0.0) &&
            (beta_eval[eval_count + 1] != 0.0))
        {
          std::complex<double> lambda_this =
            alpha_eval[eval_count] / beta_eval[eval_count];

          std::complex<double> lambda_next =
            alpha_eval[eval_count + 1] / beta_eval[eval_count + 1];

          // Check failure of cc-ness to within tolerance
          if (fabs(lambda_this.imag() + lambda_next.imag()) >
              Tolerance_for_ccness_check)
          {
            std::ostringstream error_stream;
            error_stream << "Non-zero imaginary part of eigenvalue "
                         << eval_count << " : " << lambda_this.imag()
                         << std::endl;
            error_stream
              << "isn't the negative of its subsequent value       : "
              << lambda_next.imag() << std::endl
              << "Their sum " << (lambda_this.imag() + lambda_next.imag())
              << " is greater than Tolerance_for_ccness_check = "
              << Tolerance_for_ccness_check << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
          if (fabs(lambda_this.real() - lambda_next.real()) >
              Tolerance_for_ccness_check)
          {
            std::ostringstream error_stream;
            error_stream << "Real parts of complex eigenvalue  " << eval_count
                         << " : " << lambda_this.real() << std::endl;
            error_stream
              << " doesn't agree with its supposed-to-be cc counterpart     : "
              << lambda_next.real() << std::endl
              << "Their difference "
              << (lambda_this.real() - lambda_next.real())
              << " is greater than Tolerance_for_ccness_check = "
              << Tolerance_for_ccness_check << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }

#endif


        // Resize the two cc eigenvectors associated with the
        // two cc eigenvalues
        eigenvector_real[eval_count].build(this->distribution_pt(), 0.0);
        eigenvector_imag[eval_count].build(this->distribution_pt(), 0.0);
        eigenvector_real[eval_count + 1].build(this->distribution_pt(), 0.0);
        eigenvector_imag[eval_count + 1].build(this->distribution_pt(), 0.0);
        for (unsigned j = 0; j < n; ++j)
        {
          eigenvector_real[eval_count][j] = eigenvector_aux[eval_count][j];
          eigenvector_imag[eval_count][j] = eigenvector_aux[eval_count + 1][j];

          eigenvector_real[eval_count + 1][j] = eigenvector_aux[eval_count][j];
          eigenvector_imag[eval_count + 1][j] =
            -eigenvector_aux[eval_count + 1][j];
        }
        eval_count += 2;
      }
    }
  }

  //==========================================================================
  /// Use LAPACK to solve a complex eigen problem specified by the given
  /// matrices. Note: the (real) shift
  /// that's specifiable via the EigenSolver base class is ignored here.
  /// A warning gets issued if it's set to a nonzero value.
  //==========================================================================
  void LAPACK_QZ::find_eigenvalues(
    const ComplexMatrixBase& A,
    const ComplexMatrixBase& M,
    Vector<std::complex<double>>& eigenvalue,
    Vector<Vector<std::complex<double>>>& eigenvector)
  {
#ifdef PARANOID
    if (Sigma_real != 0.0)
    {
      std::stringstream error_stream;
      error_stream << "Non-zero shift Sigma_real = " << Sigma_real
                   << " ignored in LAPACK_QZ::find_eigenvalues\n";
      OomphLibWarning(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Some character identifiers for use in the LAPACK routine
    // Do not calculate the left eigenvectors
    char no_eigvecs[2] = "N";

    // Do caculate the eigenvectors
    char eigvecs[2] = "V";

    // Get the dimension of the matrix
    int n = A.nrow(); // Total size of matrix

    // Storage for the matrices in the packed form required by the LAPACK
    // routine
    double* M_linear = new double[2 * n * n];
    double* A_linear = new double[2 * n * n];

    // Now convert the matrices into the appropriate packed form
    unsigned index = 0;
    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < n; ++j)
      {
        M_linear[index] = M(j, i).real();
        A_linear[index] = A(j, i).real();
        ++index;
        M_linear[index] = M(j, i).imag();
        A_linear[index] = A(j, i).imag();
        ++index;
      }
    }

    // Temporary eigenvalue storage
    double* alpha = new double[2 * n];
    double* beta = new double[2 * n];
    // Temporary eigenvector storage
    double* vec_left = new double[2];
    const int leading_dimension_vec_left = 1;
    double* vec_right = new double[2 * n * n];
    const int leading_dimension_vec_right = n;

    // Workspace for the LAPACK routine
    std::vector<double> work(2, 0.0);
    // The dimension of the workspace array. Set to -1 so that a workspace
    // query is assumed and LAPACK calculates the optimum size which is
    // returned as the first value of work.
    const int query_workspace = -1;
    std::vector<double> rwork(8 * n, 0.0);

    // Info flag for the LAPACK routine
    int info = 0;

    // Call FORTRAN LAPACK to get the required workspace
    LAPACK_ZGGEV(no_eigvecs,
                 eigvecs,
                 n,
                 &A_linear[0],
                 n,
                 &M_linear[0],
                 n,
                 alpha,
                 beta,
                 vec_left,
                 leading_dimension_vec_left,
                 vec_right,
                 leading_dimension_vec_right,
                 &work[0],
                 query_workspace,
                 &rwork[0],
                 info);

    // Succesful completion?
    if (info != 0)
    {
      ZGGEV_error(info, n);
    }


    // Get the amount of requires workspace
    int required_workspace = (int)work[0];
    // If we need it
    work.resize(2 * required_workspace);

    // call FORTRAN LAPACK again with the optimum workspace
    LAPACK_ZGGEV(no_eigvecs,
                 eigvecs,
                 n,
                 &A_linear[0],
                 n,
                 &M_linear[0],
                 n,
                 alpha,
                 beta,
                 vec_left,
                 1,
                 vec_right,
                 n,
                 &work[0],
                 required_workspace,
                 &rwork[0],
                 info);

    // Succesful completion?
    if (info != 0)
    {
      ZGGEV_error(info, n);
    }

    // Now resize storage for the eigenvalues and eigenvectors
    // We get them all!
    eigenvalue.resize(n);
    eigenvector.resize(n);

    // Now convert the output into our format
    for (int i = 0; i < n; ++i)
    {
      // We have temporary complex numbers
      std::complex<double> num(alpha[2 * i], alpha[2 * i + 1]);
      std::complex<double> den(beta[2 * i], beta[2 * i + 1]);

      // N.B. This is naive and dangerous according to the documentation
      // beta could be very small giving over or under flow
      eigenvalue[i] = num / den;

      // Resize the eigenvector storage
      eigenvector[i].resize(n);

      // Copy across into the complex valued eigenvector
      for (int k = 0; k < n; ++k)
      {
        eigenvector[i][k] = std::complex<double>(
          vec_right[2 * i * n + 2 * k], vec_right[2 * i * n + 2 * k + 1]);
      }
    }

    // Delete all allocated storage
    delete[] vec_right;
    delete[] vec_left;
    delete[] beta;
    delete[] alpha;
    delete[] A_linear;
    delete[] M_linear;
  }
} // namespace oomph
