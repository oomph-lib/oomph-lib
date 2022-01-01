// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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

// Include cfortran.h and the header for the FORTRAN ARPACK routines
#include "cfortran.h"
#include "arpack.h"
#include "lapack_qz.h"

// Oomph-lib headers
#include "eigen_solver.h"
#include "linear_solver.h"
#include "problem.h"


namespace oomph
{
  //===============================================================
  /// Constructor, set default values and set the initial
  /// linear solver to be superlu
  //===============================================================
  ARPACK::ARPACK()
    : EigenSolver(),
      Spectrum(1),
      NArnoldi(30),
      Small(true),
      Compute_eigenvectors(true)
  {
    Default_linear_solver_pt = Linear_solver_pt = new SuperLUSolver;
  }

  //===============================================================
  /// Destructor, delete the default linear solver
  //===============================================================
  ARPACK::~ARPACK()
  {
    delete Default_linear_solver_pt;
  }


  //==========================================================================
  /// Use ARPACK to solve an eigen problem that is assembled by elements in
  /// a mesh in a Problem object.
  //==========================================================================
  void ARPACK::solve_eigenproblem(Problem* const& problem_pt,
                                  const int& n_eval,
                                  Vector<std::complex<double>>& eigenvalue,
                                  Vector<DoubleVector>& eigenvector)
  {
    bool Verbose = false;

    // Control parameters
    int ido = 0; // Reverse communication flag
    char bmat[2]; // Specifies the type of matrix on the RHS
    // Must be 'I' for standard eigenvalue problem
    //        'G' for generalised eigenvalue problem
    int n; // Dimension of the problem
    char which[3]; // Set which eigenvalues are required.
    int nev; // Number of eigenvalues desired
    double tol = 0.0; // Stopping criteria
    int ncv; // Number of columns of the matrix V (Dimension of Arnoldi
             // subspace)
    int iparam[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // Control parameters
    // Pointers to starting locations in the workspace vectors
    int ipntr[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int info =
      0; // Setting to zero gives random initial guess for vector to start
    // Arnoldi iteration

    // Set up the sizes of the matrix
    n = problem_pt->ndof(); // Total size of matrix
    nev = n_eval; // Number of desired eigenvalues
    ncv = NArnoldi < n ? NArnoldi : n; // Number of Arnoldi vectors allowed
                                       // Maximum possible value is max
                                       // dimension of matrix

    // If we don't have enough Arnoldi vectors to compute the desired number
    // of eigenvalues then complain
    if (nev > ncv)
    {
      std::ostringstream warning_stream;
      warning_stream << "Number of requested eigenvalues " << nev << "\n"
                     << "is greater than the number of Arnoldi vectors:" << ncv
                     << "\n";
      // Increase the number of Arnoldi vectors
      ncv = nev + 10;
      if (ncv > n)
      {
        ncv = n;
      }

      warning_stream << "Increasing number of Arnoldi vectors to " << ncv
                     << "\n but you may want to increase further using\n"
                     << "ARPACK::narnoldi()\n"
                     << "which will also get rid of this warning.\n";

      OomphLibWarning(
        warning_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Allocate some workspace
    int lworkl = 3 * ncv * ncv + 6 * ncv;

    // Array that holds the final set of Arnoldi basis vectors
    double** v = new double*;
    *v = new double[ncv * n];
    // Leading dimension of v (n)
    int ldv = n;

    // Residual vector
    double* resid = new double[n];
    // Workspace vector
    double* workd = new double[3 * n];
    // Workspace vector
    double* workl = new double[lworkl];

    bmat[0] = 'G';
    // If we require eigenvalues to the left of the shift
    if (Small)
    {
      which[0] = 'S';
    }
    // Otherwise we require eigenvalues to the right of the shift
    else
    {
      which[0] = 'L';
    }

    // Which part of the eigenvalue interests us
    switch (Spectrum)
    {
        // Imaginary part
      case -1:
        which[1] = 'I';
        break;

        // Magnitude
      case 0:
        which[1] = 'M';
        break;

        // Real part
      case 1:
        which[1] = 'R';
        break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "Spectrum is set to an invalid value. It must be 0, 1 or -1\n"
          << "Not " << Spectrum << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    //     %---------------------------------------------------%
    //     | This program uses exact shifts with respect to    |
    //     | the current Hessenberg matrix (IPARAM(1) = 1).    |
    //     | IPARAM(3) specifies the maximum number of Arnoldi |
    //     | iterations allowed.  Mode 3 of DNAUPD is used     |
    //     | (IPARAM(7) = 3). All these options can be changed |
    //     | by the user. For details see the documentation in |
    //     | DNAUPD.                                           |
    //     %---------------------------------------------------%

    int ishifts = 1;
    int maxitr = 300;
    int mode = 3; // M is symetric and semi-definite

    iparam[0] = ishifts; // Exact shifts wrt Hessenberg matrix H
    iparam[2] = maxitr; // Maximum number of allowed iterations
    iparam[3] = 1; // Bloacksize to be used in the recurrence
    iparam[6] = mode; // Mode is shift and invert in real arithmetic

    // Real and imaginary parts of the shift
    double sigmar = Sigma_real, sigmai = 0.0;

    // TEMPORARY
    // only use non distributed matrice and vectors
    LinearAlgebraDistribution dist(problem_pt->communicator_pt(), n, false);
    this->build_distribution(dist);

    // Before we get into the Arnoldi loop solve the shifted matrix problem
    // Allocated Row compressed matrices for the mass matrix and shifted main
    // matrix
    CRDoubleMatrix M(this->distribution_pt()), AsigmaM(this->distribution_pt());

    // Assemble the matrices
    problem_pt->get_eigenproblem_matrices(M, AsigmaM, sigmar);

    // Allocate storage for the vectors to be used in matrix vector products
    DoubleVector rhs(this->distribution_pt(), 0.0);
    DoubleVector x(this->distribution_pt(), 0.0);


    bool LOOP_FLAG = true;
    bool First = true;

    // Enable resolves for the linear solver
    Linear_solver_pt->enable_resolve();

    // Do not report the time taken
    Linear_solver_pt->disable_doc_time();

    do
    {
      //        %---------------------------------------------%
      //        | Repeatedly call the routine DNAUPD and take |
      //        | actions indicated by parameter IDO until    |
      //        | either convergence is indicated or maxitr   |
      //        | has been exceeded.                          |
      //        %---------------------------------------------%
      //

      DNAUPD(ido,
             bmat,
             n,
             which,
             nev,
             tol,
             resid,
             ncv,
             v,
             ldv,
             iparam,
             ipntr,
             workd,
             workl,
             lworkl,
             info);

      if (ido == -1)
      {
        //           %-------------------------------------------%
        //           | Perform matrix vector multiplication      |
        //           |                y <--- OP*x                |
        //           | The user should supply his/her own        |
        //           | matrix vector multiplication routine here |
        //           | that takes workd(ipntr(1)) as the input   |
        //           | vector, and return the matrix vector      |
        //           | product to workd(ipntr(2)).               |
        //           %-------------------------------------------%
        //

        // Firstly multiply by mx
        for (int i = 0; i < n; i++)
        {
          x[i] = workd[ipntr[0] - 1 + i];
        }
        M.multiply(x, rhs);

        // Now solve the system
        DoubleVector temp(rhs);
        if (First)
        {
          Linear_solver_pt->solve(&AsigmaM, temp, rhs);
          First = false;
        }
        else
        {
          Linear_solver_pt->resolve(temp, rhs);
        }
        temp.clear();
        // AsigmaM.lubksub(rhs);
        // Put the solution into the workspace
        for (int i = 0; i < n; i++)
        {
          workd[ipntr[1] - 1 + i] = rhs[i];
        }

        continue;
      }
      else if (ido == 1)
      {
        // Already done multiplication by Mx
        // Need to load the rhs vector
        for (int i = 0; i < n; i++)
        {
          rhs[i] = workd[ipntr[2] - 1 + i];
        }
        // Now solve the system
        // AsigmaM.lubksub(rhs);
        DoubleVector temp(rhs);
        if (First)
        {
          Linear_solver_pt->solve(&AsigmaM, temp, rhs);
          First = false;
        }
        else
        {
          Linear_solver_pt->resolve(temp, rhs);
        }
        // Put the solution into the workspace
        for (int i = 0; i < n; i++)
        {
          workd[ipntr[1] - 1 + i] = rhs[i];
        }
        continue;
      }
      else if (ido == 2)
      {
        // Need to multiply by mass matrix
        // Vector<double> x(n);
        for (int i = 0; i < n; i++)
        {
          x[i] = workd[ipntr[0] - 1 + i];
        }
        M.multiply(x, rhs);

        for (int i = 0; i < n; i++)
        {
          workd[ipntr[1] - 1 + i] = rhs[i];
        }
        continue;
      }

      if (info < 0)
      {
        oomph_info << "ERROR" << std::endl;
        oomph_info << "Error with _naupd, info = '" << info << std::endl;
        if (info == -7)
        {
          oomph_info << "Not enough memory " << std::endl;
        }
      }

      //        %-------------------------------------------%
      //        | No fatal errors occurred.                 |
      //        | Post-Process using DNEUPD.                |
      //        |                                           |
      //        | Computed eigenvalues may be extracted.    |
      //        |                                           |
      //        | Eigenvectors may also be computed now if  |
      //        | desired.  (indicated by rvec = .true.)    |
      //        %-------------------------------------------%
      //
      LOOP_FLAG = false;
    } while (LOOP_FLAG);

    // Report any problems
    if (info == 1)
    {
      oomph_info << "Maximum number of iterations reached." << std::endl;
    }
    else if (info == 3)
    {
      oomph_info << "No shifts could be applied during implicit Arnoldi "
                 << "update, try increasing NCV. " << std::endl;
    }

    // Disable resolves for the linear solver
    Linear_solver_pt->disable_resolve();

    // Compute Ritz or Schur vectors, if desired
    int rvec = Compute_eigenvectors;
    // Specify the form of the basis computed
    char howmany[2];
    howmany[0] = 'A';
    // Find out the number of converged eigenvalues
    int nconv = iparam[4];

    // Note that there is an error (feature) in ARPACK that
    // means in certain cases, if we request eigenvectors,
    // consequent Schur factorization of the matrix spanned
    // by the Arnoldi vectors leads to more converged eigenvalues
    // than reported by DNAPUD. This is a pain because it
    // causes a segmentation fault and in every case that I've found
    // the eigenvalues/vectors are not those desired.
    // At the moment, I'm just going to leave it, but I note
    // the problem here to remind myself about it.

    // Array to select which Ritz vectors are computed
    int select[ncv];
    // Array that holds the real part of the eigenvalues
    double* eval_real = new double[nconv + 1];
    // Array that holds the imaginary part of the eigenvalues
    double* eval_imag = new double[nconv + 1];

    // Workspace
    double* workev = new double[3 * ncv];
    // Array to hold the eigenvectors
    double** z = new double*;
    *z = new double[(nconv + 1) * n];

    // Error flag
    int ierr;

    DNEUPD(rvec,
           howmany,
           select,
           eval_real,
           eval_imag,
           z,
           ldv,
           sigmar,
           sigmai,
           workev,
           bmat,
           n,
           which,
           nev,
           tol,
           resid,
           ncv,
           v,
           ldv,
           iparam,
           ipntr,
           workd,
           workl,
           lworkl,
           ierr);
    //
    //        %-----------------------------------------------%
    //        | The real part of the eigenvalue is returned   |
    //        | in the first column of the two dimensional    |
    //        | array D, and the imaginary part is returned   |
    //        | in the second column of D.  The corresponding |
    //        | eigenvectors are returned in the first NEV    |
    //        | columns of the two dimensional array V if     |
    //        | requested.  Otherwise, an orthogonal basis    |
    //        | for the invariant subspace corresponding to   |
    //        | the eigenvalues in D is returned in V.        |
    //        %-----------------------------------------------%

    if (ierr != 0)
    {
      oomph_info << "Error with _neupd, info = " << ierr << std::endl;
    }
    else
    // Print the output anyway
    {
      // Resize the eigenvalue output vector
      eigenvalue.resize(nconv);
      for (int j = 0; j < nconv; j++)
      {
        // Turn the output from ARPACK into a complex number
        std::complex<double> eigen(eval_real[j], eval_imag[j]);
        // Add the eigenvalues to the output vector
        eigenvalue[j] = eigen;
      }

      if (Compute_eigenvectors)
      {
        // Load the eigenvectors
        eigenvector.resize(nconv);
        for (int j = 0; j < nconv; j++)
        {
          eigenvector[j].build(this->distribution_pt(), 0.0);
          for (int i = 0; i < n; i++)
          {
            eigenvector[j][i] = z[0][j * n + i];
          }
        }
      }
      else
      {
        eigenvector.resize(0);
      }

      // Report the information
      if (Verbose)
      {
        oomph_info << " ARPACK " << std::endl;
        oomph_info << " ====== " << std::endl;
        oomph_info << " Size of the matrix is " << n << std::endl;
        oomph_info << " The number of Ritz values requested is " << nev
                   << std::endl;
        oomph_info << " The number of Arnoldi vectors generated (NCV) is "
                   << ncv << std::endl;
        oomph_info << " What portion of the spectrum: " << which << std::endl;
        oomph_info << " The number of converged Ritz values is " << nconv
                   << std::endl;
        oomph_info
          << " The number of Implicit Arnoldi update iterations taken is "
          << iparam[2] << std::endl;
        oomph_info << " The number of OP*x is " << iparam[8] << std::endl;
        oomph_info << " The convergence criterion is " << tol << std::endl;
      }
    }

    // Clean up the allocated memory
    delete[] * z;
    *z = 0;
    delete z;
    z = 0;
    delete[] workev;
    workev = 0;
    delete[] eval_imag;
    eval_imag = 0;
    delete[] eval_real;
    eval_real = 0;

    // Clean up the memory
    delete[] workl;
    workl = 0;
    delete[] workd;
    workd = 0;
    delete[] resid;
    resid = 0;
    delete[] * v;
    *v = 0;
    delete v;
    v = 0;
  }


  //==========================================================================
  /// Use LAPACK to solve an eigen problem that is assembled by elements in
  /// a mesh in a Problem object.
  //==========================================================================
  void LAPACK_QZ::solve_eigenproblem(Problem* const& problem_pt,
                                     const int& n_eval,
                                     Vector<std::complex<double>>& eigenvalue,
                                     Vector<DoubleVector>& eigenvector)
  {
    // Some character identifiers for use in the LAPACK routine
    // Do not calculate the left eigenvectors
    char no_eigvecs[2] = "N";
    // Do caculate the eigenvectors
    char eigvecs[2] = "V";

    // Get the dimension of the matrix
    int n = problem_pt->ndof(); // Total size of matrix

    // Initialise a padding integer
    int padding = 0;
    // If the dimension of the matrix is even, then pad the arrays to
    // make the size odd. This somehow sorts out a strange run-time behaviour
    // identified by Rich Hewitt.
    if (n % 2 == 0)
    {
      padding = 1;
    }

    // Get the real and imaginary parts of the shift
    double sigmar = Sigma_real; // sigmai = 0.0;

    // Actual size of matrix that will be allocated
    int padded_n = n + padding;

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

      // Assemble the matrices
      problem_pt->get_eigenproblem_matrices(temp_M, temp_AsigmaM, sigmar);

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
        if (padding)
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
    double* vec_right = new double[n * n];

    // Workspace for the LAPACK routine
    std::vector<double> work(1, 0.0);
    // Info flag for the LAPACK routine
    int info = 0;

    // Call FORTRAN LAPACK to get the required workspace
    // Note the use of the padding flag
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
                 1,
                 vec_right,
                 n,
                 &work[0],
                 -1,
                 info);

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
                 1,
                 vec_right,
                 n,
                 &work[0],
                 required_workspace,
                 info);

    // Now resize storage for the eigenvalues and eigenvectors
    // We get them all!
    eigenvalue.resize(n);
    eigenvector.resize(n);

    // Now convert the output into our format
    for (int i = 0; i < n; ++i)
    {
      // N.B. This is naive and dangerous according to the documentation
      // beta could be very small giving over or under flow
      // Remember to fix the shift back again
      eigenvalue[i] = std::complex<double>(sigmar + alpha_r[i] / beta[i],
                                           alpha_i[i] / beta[i]);

      // Resize the eigenvector  storage
      eigenvector[i].build(this->distribution_pt(), 0.0);
      // Load up the eigenvector (assume that it's real)
      for (int k = 0; k < n; ++k)
      {
        eigenvector[i][k] = vec_right[i * n + k];
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
  /// Use LAPACK to solve a complex eigen problem specified by the given
  /// matrices.
  //==========================================================================
  void LAPACK_QZ::find_eigenvalues(
    const ComplexMatrixBase& A,
    const ComplexMatrixBase& M,
    Vector<std::complex<double>>& eigenvalue,
    Vector<Vector<std::complex<double>>>& eigenvector)
  {
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
    double* vec_right = new double[2 * n * n];

    // Workspace for the LAPACK routine
    std::vector<double> work(2, 0.0);
    std::vector<double> rwork(8 * n, 0.0);

    // Info flag for the LAPACK routine
    int info = 0;

    // Call FORTRAN LAPACK to get the required workspace
    // Note the use of the padding flag
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
                 -1,
                 &rwork[0],
                 info);

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
      // Load up the eigenvector (assume that it's real)
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
