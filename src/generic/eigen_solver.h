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
// Header file for a class that defines interfaces to Eigensolvers

// Include guard to prevent multiple inclusions of the header
#ifndef OOMPH_EIGEN_SOLVER_HEADER
#define OOMPH_EIGEN_SOLVER_HEADER

// Include the header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

#include <complex>
#include "Vector.h"
#include "complex_matrices.h"

namespace oomph
{
  // Forward definition of problem class
  class Problem;

  // Forward definition of matrix class
  class DoubleMatrixBase;

  // Forward definition of linear solver class
  class LinearSolver;

  //=======================================================================
  /// Base class for all EigenProblem solves. This simply defines standard
  /// interfaces so that different solvers can be used easily.
  //=======================================================================
  class EigenSolver : public DistributableLinearAlgebraObject
  {
  protected:
    /// Double value that represents the real part of the shift in
    /// shifted eigensolvers
    double Sigma_real;

  public:
    /// Empty constructor
    EigenSolver() : Sigma_real(0.0) {}

    /// Empty copy constructor
    EigenSolver(const EigenSolver&) {}

    /// Empty destructor
    virtual ~EigenSolver() {}

    /// Eigensolver. This takes a pointer to a problem and returns
    /// a vector of complex numbers representing the eigenvalues
    /// and a corresponding vector of eigenvectors. n_eval specifies the min.
    /// number of eigenvalues/vectors required. This is primarily used in
    /// Arnoldi type implementations; direct solvers such as QZ compute all the
    /// eigenvalues/vectors.
    /// Note: this is a legacy version of this function that stores re & imag
    /// parts of eigenvectors in some solver-specific collection of real
    /// vectors.
    virtual void solve_eigenproblem_legacy(
      Problem* const& problem_pt,
      const int& n_eval,
      Vector<std::complex<double>>& eigenvalue,
      Vector<DoubleVector>& eigenvector,
      const bool& do_adjoint_problem = false) = 0;


    /// Solve the real eigenproblem that is assembled by elements in
    /// a mesh in a Problem object. Note that the assembled matrices include the
    /// shift and are real. The eigenvalues and eigenvectors are,
    /// in general, complex, and the eigenvalues can be NaNs or Infs.
    /// This function is therefore merely provided as a convenience, to be
    /// used if the user is confident that NaNs don't arise (e.g. in Arnoldi
    /// based solvers where typically only a small number of (finite)
    /// eigenvalues are computed), or if the users is happy to deal with NaNs in
    /// the subsequent post-processing.
    /// Function is virtual so it can be overloaded for Arnoldi type solvers
    /// that compute the (finite) eigenvalues directly
    /// At least n_eval eigenvalues are computed.
    virtual void solve_eigenproblem(Problem* const& problem_pt,
                                    const int& n_eval,
                                    Vector<std::complex<double>>& eigenvalue,
                                    Vector<DoubleVector>& eigenvector_real,
                                    Vector<DoubleVector>& eigenvector_imag,
                                    const bool& do_adjoint_problem = false)
    {
      Vector<std::complex<double>> alpha;
      Vector<double> beta;

      // Call the "safe" version
      solve_eigenproblem(problem_pt,
                         n_eval,
                         alpha,
                         beta,
                         eigenvector_real,
                         eigenvector_imag,
                         do_adjoint_problem);

      // Now do the brute force conversion, possibly creating NaNs and Infs...
      unsigned n = alpha.size();
      eigenvalue.resize(n);
      for (unsigned i = 0; i < n; i++)
      {
        eigenvalue[i] = alpha[i] / beta[i];
      }
    }

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
    virtual void solve_eigenproblem(Problem* const& problem_pt,
                                    const int& n_eval,
                                    Vector<std::complex<double>>& alpha,
                                    Vector<double>& beta,
                                    Vector<DoubleVector>& eigenvector_real,
                                    Vector<DoubleVector>& eigenvector_imag,
                                    const bool& do_adjoint_problem = false) = 0;

    /// Set the value of the (real) shift
    void set_shift(const double& shift_value)
    {
      Sigma_real = shift_value;
    }

    /// Return the value of the (real) shift (const version)
    const double& get_shift() const
    {
      return Sigma_real;
    }
  };

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////


  //=====================================================================
  /// Class for the ARPACK eigensolver
  //=====================================================================
  class ARPACK : public EigenSolver
  {
  private:
    /// Pointer to a linear solver
    LinearSolver* Linear_solver_pt;

    /// Pointer to a default linear solver
    LinearSolver* Default_linear_solver_pt;

    /// Integer to set whether the real, imaginary or magnitude is
    /// required
    /// to be small or large.
    int Spectrum;

    /// Number of Arnoldi vectors to compute
    int NArnoldi;

    /// Boolean to set which part of the spectrum left (default) or right
    /// of the shifted value.
    bool Small;

    /// Boolean to indicate whether or not to compute the eigenvectors
    bool Compute_eigenvectors;


  public:
    /// Constructor
    ARPACK();

    /// Empty copy constructor
    ARPACK(const ARPACK&) {}

    /// Destructor, delete the linear solver
    virtual ~ARPACK();

    /// Access function for the number of Arnoldi vectors
    int& narnoldi()
    {
      return NArnoldi;
    }

    /// Access function for the number of Arnoldi vectors (const version)
    const int& narnoldi() const
    {
      return NArnoldi;
    }

    /// Set to enable the computation of the eigenvectors (default)
    void enable_compute_eigenvectors()
    {
      Compute_eigenvectors = true;
    }

    /// Set to disable the computation of the eigenvectors
    void disable_compute_eigenvectors()
    {
      Compute_eigenvectors = false;
    }

    /// Solve the eigen problem
    void solve_eigenproblem_legacy(Problem* const& problem_pt,
                                   const int& n_eval,
                                   Vector<std::complex<double>>& eigenvalue,
                                   Vector<DoubleVector>& eigenvector,
                                   const bool& do_adjoint_problem = false);


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
    void solve_eigenproblem(Problem* const& problem_pt,
                            const int& n_eval,
                            Vector<std::complex<double>>& alpha,
                            Vector<double>& beta,
                            Vector<DoubleVector>& eigenvector_real,
                            Vector<DoubleVector>& eigenvector_imag,
                            const bool& do_adjoint_problem = false)
    {
      oomph_info << "Broken, but then don't we want arpack to go anyway?\n";
      abort();
    }

    /// Set the desired eigenvalues to be left of the shift
    void get_eigenvalues_left_of_shift()
    {
      Small = true;
    }

    /// Set the desired eigenvalues to be right of the shift
    void get_eigenvalues_right_of_shift()
    {
      Small = false;
    }

    /// Set the real part to be the quantity of interest (default)
    void track_eigenvalue_real_part()
    {
      Spectrum = 1;
    }

    /// Set the imaginary part fo the quantity of interest
    void track_eigenvalue_imaginary_part()
    {
      Spectrum = -1;
    }

    /// Set the magnitude to be the quantity of interest
    void track_eigenvalue_magnitude()
    {
      Spectrum = 0;
    }

    /// Return a pointer to the linear solver object
    LinearSolver*& linear_solver_pt()
    {
      return Linear_solver_pt;
    }

    /// Return a pointer to the linear solver object (const version)
    LinearSolver* const& linear_solver_pt() const
    {
      return Linear_solver_pt;
    }
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  //=====================================================================
  /// Class for the LAPACK QZ eigensolver
  //=====================================================================
  class LAPACK_QZ : public EigenSolver
  {
  public:
    /// Empty constructor
    LAPACK_QZ() : EigenSolver(), Tolerance_for_ccness_check(1.0e-13) {}

    /// Broken copy constructor
    LAPACK_QZ(const LAPACK_QZ&) = delete;

    /// Broken assignment operator
    void operator=(const LAPACK_QZ&) = delete;

    /// Empty desctructor
    virtual ~LAPACK_QZ() {}

    /// Use LAPACK QZ to solve the real eigenproblem that is assembled
    /// by elements in a mesh in a Problem object. Note that the assembled
    /// matrices include the shift and are real. The eigenvalues and
    /// eigenvectors are, in general, complex.
    /// This is a legacy version of this function that stores re & imag parts of
    /// eigenvectors in some solver-specific collection of real vectors;
    /// they are disentangled in the alternative version of this function
    /// that returns Vectors of complex Vectors.
    /// At least n_eval eigenvalues are computed.
    void solve_eigenproblem_legacy(Problem* const& problem_pt,
                                   const int& n_eval,
                                   Vector<std::complex<double>>& eigenvalue,
                                   Vector<DoubleVector>& eigenvector,
                                   const bool& do_adjoint_problem = false);

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
    void solve_eigenproblem(Problem* const& problem_pt,
                            const int& n_eval,
                            Vector<std::complex<double>>& alpha,
                            Vector<double>& beta,
                            Vector<DoubleVector>& eigenvector_real,
                            Vector<DoubleVector>& eigenvector_imag,
                            const bool& do_adjoint_problem = false);

    /// Find the eigenvalues of a complex generalised eigenvalue problem
    /// specified by \f$ Ax = \lambda  Mx \f$. Note: the (real) shift
    /// that's specifiable via the EigenSolver base class is ignored here.
    /// A warning gets issued if it's set to a nonzero value.
    void find_eigenvalues(const ComplexMatrixBase& A,
                          const ComplexMatrixBase& M,
                          Vector<std::complex<double>>& eigenvalue,
                          Vector<Vector<std::complex<double>>>& eigenvector);


    /// Access to tolerance for checking complex conjugateness of eigenvalues
    /// (const version)
    double tolerance_for_ccness_check() const
    {
      return Tolerance_for_ccness_check;
    }


    /// Access to tolerance for checking complex conjugateness of eigenvalues
    double& tolerance_for_ccness_check()
    {
      return Tolerance_for_ccness_check;
    }

  private:
    /// Helper function called from legacy and updated version from "raw" lapack
    /// code
    void solve_eigenproblem_helper(Problem* const& problem_pt,
                                   const int& n_eval,
                                   Vector<std::complex<double>>& alpha,
                                   Vector<double>& beta,
                                   Vector<DoubleVector>& eigenvector);


    /// Provide diagonstic for DGGEV error return
    void DGGEV_error(const int& info, const int& n)
    {
      // Throw an error
      std::ostringstream error_stream;
      error_stream << "Failure in LAPACK_DGGEV(...).\n"
                   << "info = " << info << std::endl;
      error_stream
        << "Diagnostics below are from \n\n"
        << "http://www.netlib.org/lapack/explore-html/d9/d8e/"
           "group__double_g_eeigen_ga4f59e87e670a755b41cbdd7e97f36bea.html"
        << std::endl;
      if (info < 0)
      {
        error_stream << -info << "-th input arg had an illegal value\n";
      }
      else if (info <= n)
      {
        error_stream << "The QZ iteration failed.  No eigenvectors have been\n"
                     << "calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)\n"
                     << "should be correct for j=INFO+1,...,N, where \n"
                     << "info = " << info << " and N = " << n << std::endl;
      }
      else if (info == (n + 1))
      {
        error_stream << "QZ iteration failed in DHGEQZ.\n";
      }
      else if (info == (n + 2))
      {
        error_stream << "error return from DTGEVC.\n";
      }
      error_stream
        << "Aborting here; if you know how to proceed then\n"
        << "then implement ability to catch this error and continue\n";

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// Provide diagonstic for ZGGEV error return
    void ZGGEV_error(const int& info, const int& n)
    {
      // Throw an error
      std::ostringstream error_stream;
      error_stream << "Failure in LAPACK_ZGGEV(...).\n"
                   << "info = " << info << std::endl;
      error_stream << "Diagnostics below are from \n\n"
                   << "http://www.netlib.org/lapack/explore-html/"
                   << "db/d55/group__complex16_g_eeigen_ga79fcce20c"
                   << "617429ccf985e6f123a6171.html" << std::endl;
      if (info < 0)
      {
        error_stream << -info << "-th input arg had an illegal value\n";
      }
      else if (info <= n)
      {
        error_stream << "The QZ iteration failed.  No eigenvectors have been\n"
                     << "calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)\n"
                     << "should be correct for j=INFO+1,...,N, where \n"
                     << "info = " << info << " and N = " << n << std::endl;
      }
      else if (info == (n + 1))
      {
        error_stream << "QZ iteration failed in ZHGEQZ.\n";
      }
      else if (info == (n + 2))
      {
        error_stream << "error return from ZTGEVC.\n";
      }
      error_stream
        << "Aborting here; if you know how to proceed then\n"
        << "then implement ability to catch this error and continue\n";

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// Tolerance for checking complex conjugateness of eigenvalues
    double Tolerance_for_ccness_check;
  };

} // namespace oomph

#endif
