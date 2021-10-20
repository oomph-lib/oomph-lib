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
// Non-inline member functions for the matrix classes

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

#include <cstring>

#include <set>
#include <map>

//#include <valgrind/callgrind.h>

// oomph-lib headers
#include "matrices.h"
#include "linear_solver.h"


namespace oomph
{
  //============================================================================
  /// Complete LU solve (overwrites RHS with solution). This is the
  /// generic version which should not need to be over-written.
  //============================================================================
  void DoubleMatrixBase::solve(DoubleVector& rhs)
  {
#ifdef PARANOID
    if (Linear_solver_pt == 0)
    {
      throw OomphLibError("Linear_solver_pt not set in matrix",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Copy rhs vector into local storage so it doesn't get overwritten
    // if the linear solver decides to initialise the solution vector, say,
    // which it's quite entitled to do!
    DoubleVector actual_rhs(rhs);

    // Use the linear algebra interface to the linear solver
    Linear_solver_pt->solve(this, actual_rhs, rhs);
  }

  //============================================================================
  /// Complete LU solve (Nothing gets overwritten!). This generic
  /// version should never need to be overwritten
  //============================================================================
  void DoubleMatrixBase::solve(const DoubleVector& rhs, DoubleVector& soln)
  {
#ifdef PARANOID
    if (Linear_solver_pt == 0)
    {
      throw OomphLibError("Linear_solver_pt not set in matrix",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif
    // Use the linear algebra interface to the linear solver
    Linear_solver_pt->solve(this, rhs, soln);
  }

  //============================================================================
  /// Complete LU solve (overwrites RHS with solution). This is the
  /// generic version which should not need to be over-written.
  //============================================================================
  void DoubleMatrixBase::solve(Vector<double>& rhs)
  {
#ifdef PARANOID
    if (Linear_solver_pt == 0)
    {
      throw OomphLibError("Linear_solver_pt not set in matrix",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Copy rhs vector into local storage so it doesn't get overwritten
    // if the linear solver decides to initialise the solution vector, say,
    // which it's quite entitled to do!
    Vector<double> actual_rhs(rhs);

    // Use the linear algebra interface to the linear solver
    Linear_solver_pt->solve(this, actual_rhs, rhs);
  }

  //============================================================================
  /// Complete LU solve (Nothing gets overwritten!). This generic
  /// version should never need to be overwritten
  //============================================================================
  void DoubleMatrixBase::solve(const Vector<double>& rhs, Vector<double>& soln)
  {
#ifdef PARANOID
    if (Linear_solver_pt == 0)
    {
      throw OomphLibError("Linear_solver_pt not set in matrix",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif
    // Use the linear algebra interface to the linear solver
    Linear_solver_pt->solve(this, rhs, soln);
  }


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //===============================================================
  /// Constructor, set the default linear solver to be the DenseLU
  /// solver
  //===============================================================
  DenseDoubleMatrix::DenseDoubleMatrix() : DenseMatrix<double>()
  {
    Linear_solver_pt = Default_linear_solver_pt = new DenseLU;
  }

  //==============================================================
  /// Constructor to build a square n by n matrix.
  /// Set the default linear solver to be DenseLU
  //==============================================================
  DenseDoubleMatrix::DenseDoubleMatrix(const unsigned long& n)
    : DenseMatrix<double>(n)
  {
    Linear_solver_pt = Default_linear_solver_pt = new DenseLU;
  }


  //=================================================================
  /// Constructor to build a matrix with n rows and m columns.
  /// Set the default linear solver to be DenseLU
  //=================================================================
  DenseDoubleMatrix::DenseDoubleMatrix(const unsigned long& n,
                                       const unsigned long& m)
    : DenseMatrix<double>(n, m)
  {
    Linear_solver_pt = Default_linear_solver_pt = new DenseLU;
  }

  //=====================================================================
  /// Constructor to build a matrix with n rows and m columns,
  /// with initial value initial_val
  /// Set the default linear solver to be DenseLU
  //=====================================================================
  DenseDoubleMatrix::DenseDoubleMatrix(const unsigned long& n,
                                       const unsigned long& m,
                                       const double& initial_val)
    : DenseMatrix<double>(n, m, initial_val)
  {
    Linear_solver_pt = Default_linear_solver_pt = new DenseLU;
  }

  //=======================================================================
  /// Destructor delete the default linear solver
  //======================================================================
  DenseDoubleMatrix::~DenseDoubleMatrix()
  {
    // Delete the default linear solver
    delete Default_linear_solver_pt;
  }

  //============================================================================
  /// LU decompose a matrix, by using the default linear solver
  /// (DenseLU)
  //============================================================================
  void DenseDoubleMatrix::ludecompose()
  {
    // Use the default (DenseLU) solver to ludecompose the matrix
    static_cast<DenseLU*>(Default_linear_solver_pt)->factorise(this);
  }


  //============================================================================
  ///  Back substitute an LU decomposed matrix.
  //============================================================================
  void DenseDoubleMatrix::lubksub(DoubleVector& rhs)
  {
    // Use the default (DenseLU) solver to perform the backsubstitution
    static_cast<DenseLU*>(Default_linear_solver_pt)->backsub(rhs, rhs);
  }

  //============================================================================
  ///  Back substitute an LU decomposed matrix.
  //============================================================================
  void DenseDoubleMatrix::lubksub(Vector<double>& rhs)
  {
    // Use the default (DenseLU) solver to perform the backsubstitution
    static_cast<DenseLU*>(Default_linear_solver_pt)->backsub(rhs, rhs);
  }


  //============================================================================
  ///  Determine eigenvalues and eigenvectors, using
  /// Jacobi rotations. Only for symmetric matrices. Nothing gets overwritten!
  /// - \c eigen_vect(i,j) = j-th component of i-th eigenvector.
  /// - \c eigen_val[i] is the i-th eigenvalue; same ordering as in eigenvectors
  //============================================================================
  void DenseDoubleMatrix::eigenvalues_by_jacobi(
    Vector<double>& eigen_vals, DenseMatrix<double>& eigen_vect) const
  {
#ifdef PARANOID
    // Check Matrix is square
    if (N != M)
    {
      throw OomphLibError(
        "This matrix is not square, the matrix MUST be square!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif
    // Make a copy of the matrix & check that it's symmetric

    // Check that the sizes of eigen_vals and eigen_vect are correct. If not
    // correct them.
    if (eigen_vals.size() != N)
    {
      eigen_vals.resize(N);
    }
    if (eigen_vect.ncol() != N || eigen_vect.nrow() != N)
    {
      eigen_vect.resize(N);
    }

    DenseDoubleMatrix working_matrix(N);
    for (unsigned long i = 0; i < N; i++)
    {
      for (unsigned long j = 0; j < M; j++)
      {
#ifdef PARANOID
        if (Matrixdata[M * i + j] != Matrixdata[M * j + i])
        {
          throw OomphLibError(
            "Matrix needs to be symmetric for eigenvalues_by_jacobi()",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
#endif
        working_matrix(i, j) = (*this)(i, j);
      }
    }

    DenseDoubleMatrix aux_eigen_vect(N);

    throw OomphLibError("Sorry JacobiEigenSolver::jacobi() removed because of "
                        "licencing problems.",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);

    // // Call eigensolver
    // unsigned long nrot;
    // JacobiEigenSolver::jacobi(working_matrix, eigen_vals, aux_eigen_vect,
    //                           nrot);

    // Copy across (and transpose)
    for (unsigned long i = 0; i < N; i++)
    {
      for (unsigned long j = 0; j < M; j++)
      {
        eigen_vect(i, j) = aux_eigen_vect(j, i);
      }
    }
  }


  //============================================================================
  ///  Multiply the matrix by the vector x: soln=Ax
  //============================================================================
  void DenseDoubleMatrix::multiply(const DoubleVector& x,
                                   DoubleVector& soln) const
  {
#ifdef PARANOID
    // Check to see if x.size() = ncol().
    if (x.nrow() != this->ncol())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.nrow() << ", it should be " << this->ncol()
                           << std::endl;
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // check that x is not distributed
    if (x.distributed())
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "The x vector cannot be distributed for DenseDoubleMatrix "
        << "matrix-vector multiply" << std::endl;
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // if soln is setup...
    if (soln.built())
    {
      // check that soln is not distributed
      if (soln.distributed())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The x vector cannot be distributed for DenseDoubleMatrix "
          << "matrix-vector multiply" << std::endl;
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      if (soln.nrow() != this->nrow())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector is setup and therefore must have the same "
          << "number of rows as the matrix";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      if (*x.distribution_pt()->communicator_pt() !=
          *soln.distribution_pt()->communicator_pt())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector and the x vector must have the same communicator"
          << std::endl;
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if soln is not setup then setup the distribution
    if (!soln.built())
    {
      LinearAlgebraDistribution dist(
        x.distribution_pt()->communicator_pt(), this->nrow(), false);
      soln.build(&dist, 0.0);
    }

    // Initialise the solution
    soln.initialise(0.0);

    // Multiply the matrix A, by the vector x
    const double* x_pt = x.values_pt();
    double* soln_pt = soln.values_pt();
    for (unsigned long i = 0; i < N; i++)
    {
      for (unsigned long j = 0; j < M; j++)
      {
        soln_pt[i] += Matrixdata[M * i + j] * x_pt[j];
      }
    }
  }


  //=================================================================
  /// Multiply the transposed matrix by the vector x: soln=A^T x
  //=================================================================
  void DenseDoubleMatrix::multiply_transpose(const DoubleVector& x,
                                             DoubleVector& soln) const
  {
#ifdef PARANOID
    // Check to see if x.size() = ncol().
    if (x.nrow() != this->nrow())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.nrow() << ", it should be " << this->nrow()
                           << std::endl;
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // check that x is not distributed
    if (x.distributed())
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "The x vector cannot be distributed for DenseDoubleMatrix "
        << "matrix-vector multiply" << std::endl;
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // if soln is setup...
    if (soln.built())
    {
      // check that soln is not distributed
      if (soln.distributed())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The x vector cannot be distributed for DenseDoubleMatrix "
          << "matrix-vector multiply" << std::endl;
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      if (soln.nrow() != this->ncol())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector is setup and therefore must have the same "
          << "number of columns as the matrix";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      if (*soln.distribution_pt()->communicator_pt() !=
          *x.distribution_pt()->communicator_pt())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector and the x vector must have the same communicator"
          << std::endl;
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if soln is not setup then setup the distribution
    if (!soln.built())
    {
      LinearAlgebraDistribution* dist_pt = new LinearAlgebraDistribution(
        x.distribution_pt()->communicator_pt(), this->ncol(), false);
      soln.build(dist_pt, 0.0);
      delete dist_pt;
    }

    // Initialise the solution
    soln.initialise(0.0);

    // Matrix vector product
    double* soln_pt = soln.values_pt();
    const double* x_pt = x.values_pt();
    for (unsigned long i = 0; i < N; i++)
    {
      for (unsigned long j = 0; j < M; j++)
      {
        soln_pt[j] += Matrixdata[N * i + j] * x_pt[i];
      }
    }
  }


  //=================================================================
  /// For every row, find the maximum absolute value of the
  /// entries in this row. Set all values that are less than alpha times
  /// this maximum to zero and return the resulting matrix in
  /// reduced_matrix. Note: Diagonal entries are retained regardless
  /// of their size.
  //=================================================================
  void DenseDoubleMatrix::matrix_reduction(const double& alpha,
                                           DenseDoubleMatrix& reduced_matrix)
  {
    reduced_matrix.resize(N, M, 0.0);
    // maximum value in a row
    double max_row;

    // Loop over rows
    for (unsigned i = 0; i < N; i++)
    {
      // Initialise max value in row
      max_row = 0.0;

      // Loop over entries in columns
      for (unsigned long j = 0; j < M; j++)
      {
        // Find max. value in row
        if (std::fabs(Matrixdata[M * i + j]) > max_row)
        {
          max_row = std::fabs(Matrixdata[M * i + j]);
        }
      }

      // Decide if we need to retain the entries in the row
      for (unsigned long j = 0; j < M; j++)
      {
        // If we're on the diagonal or the value is sufficiently large: retain
        // i.e. copy across.
        if (i == j || std::fabs(Matrixdata[M * i + j]) > alpha * max_row)
        {
          reduced_matrix(i, j) = Matrixdata[M * i + j];
        }
      }
    }
  }


  //=============================================================================
  /// Function to multiply this matrix by the DenseDoubleMatrix  matrix_in.
  //=============================================================================
  void DenseDoubleMatrix::multiply(const DenseDoubleMatrix& matrix_in,
                                   DenseDoubleMatrix& result)
  {
#ifdef PARANOID
    // check matrix dimensions are compatable
    if (this->ncol() != matrix_in.nrow())
    {
      std::ostringstream error_message;
      error_message
        << "Matrix dimensions incompatable for matrix-matrix multiplication"
        << "ncol() for first matrix:" << this->ncol()
        << "nrow() for second matrix: " << matrix_in.nrow();

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // NB N is number of rows!
    unsigned long n_row = this->nrow();
    unsigned long m_col = matrix_in.ncol();

    // resize and intialize result
    result.resize(n_row, m_col, 0.0);

    // clock_t clock1 = clock();

    // do calculation
    unsigned long n_col = this->ncol();
    for (unsigned long k = 0; k < n_col; k++)
    {
      for (unsigned long i = 0; i < n_row; i++)
      {
        for (unsigned long j = 0; j < m_col; j++)
        {
          result(i, j) += Matrixdata[m_col * i + k] * matrix_in(k, j);
        }
      }
    }
  }


  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////


  //=======================================================================
  ///  Default constructor, set the default linear solver and
  /// matrix-matrix multiplication method.
  //========================================================================
  CCDoubleMatrix::CCDoubleMatrix() : CCMatrix<double>()
  {
    Linear_solver_pt = Default_linear_solver_pt = new SuperLUSolver;
    Matrix_matrix_multiply_method = 2;
  }

  //========================================================================
  ///  Constructor: Pass vector of values, vector of row indices,
  /// vector of column starts and number of rows (can be suppressed
  /// for square matrices). Number of nonzero entries is read
  /// off from value, so make sure the vector has been shrunk
  /// to its correct length.
  //=======================================================================
  CCDoubleMatrix::CCDoubleMatrix(const Vector<double>& value,
                                 const Vector<int>& row_index,
                                 const Vector<int>& column_start,
                                 const unsigned long& n,
                                 const unsigned long& m)
    : CCMatrix<double>(value, row_index, column_start, n, m)
  {
    Linear_solver_pt = Default_linear_solver_pt = new SuperLUSolver;
    Matrix_matrix_multiply_method = 2;
  }

  /// Destructor: delete the default linear solver
  CCDoubleMatrix::~CCDoubleMatrix()
  {
    delete Default_linear_solver_pt;
  }


  //===================================================================
  /// Perform LU decomposition. Return the sign of the determinant
  //===================================================================
  void CCDoubleMatrix::ludecompose()
  {
    static_cast<SuperLUSolver*>(Default_linear_solver_pt)->factorise(this);
  }

  //===================================================================
  /// Do the backsubstitution
  //===================================================================
  void CCDoubleMatrix::lubksub(DoubleVector& rhs)
  {
    static_cast<SuperLUSolver*>(Default_linear_solver_pt)->backsub(rhs, rhs);
  }

  //===================================================================
  ///  Multiply the matrix by the vector x
  //===================================================================
  void CCDoubleMatrix::multiply(const DoubleVector& x, DoubleVector& soln) const
  {
#ifdef PARANOID
    // Check to see if x.size() = ncol().
    if (x.nrow() != this->ncol())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.nrow() << ", it should be " << this->ncol()
                           << std::endl;
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // check that x is not distributed
    if (x.distributed())
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "The x vector cannot be distributed for CCDoubleMatrix "
        << "matrix-vector multiply" << std::endl;
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // if soln is setup...
    if (soln.built())
    {
      // check that soln is not distributed
      if (soln.distributed())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The x vector cannot be distributed for CCDoubleMatrix "
          << "matrix-vector multiply" << std::endl;
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      if (soln.nrow() != this->nrow())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector is setup and therefore must have the same "
          << "number of rows as the matrix";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      if (*soln.distribution_pt()->communicator_pt() !=
          *x.distribution_pt()->communicator_pt())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector and the x vector must have the same communicator"
          << std::endl;
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if soln is not setup then setup the distribution
    if (!soln.built())
    {
      LinearAlgebraDistribution* dist_pt = new LinearAlgebraDistribution(
        x.distribution_pt()->communicator_pt(), this->nrow(), false);
      soln.build(dist_pt, 0.0);
      delete dist_pt;
    }

    // zero
    soln.initialise(0.0);

    // multiply
    double* soln_pt = soln.values_pt();
    const double* x_pt = x.values_pt();
    for (unsigned long j = 0; j < N; j++)
    {
      for (long k = Column_start[j]; k < Column_start[j + 1]; k++)
      {
        unsigned long i = Row_index[k];
        double a_ij = Value[k];
        soln_pt[i] += a_ij * x_pt[j];
      }
    }
  }


  //=================================================================
  /// Multiply the  transposed matrix by the vector x: soln=A^T x
  //=================================================================
  void CCDoubleMatrix::multiply_transpose(const DoubleVector& x,
                                          DoubleVector& soln) const
  {
#ifdef PARANOID
    // Check to see if x.size() = ncol().
    if (x.nrow() != this->nrow())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.nrow() << ", it should be " << this->nrow()
                           << std::endl;
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // check that x is not distributed
    if (x.distributed())
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "The x vector cannot be distributed for CCDoubleMatrix "
        << "matrix-vector multiply" << std::endl;
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // if soln is setup...
    if (soln.built())
    {
      // check that soln is not distributed
      if (soln.distributed())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The x vector cannot be distributed for CCDoubleMatrix "
          << "matrix-vector multiply" << std::endl;
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      if (soln.nrow() != this->ncol())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector is setup and therefore must have the same "
          << "number of columns as the matrix";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      if (*soln.distribution_pt()->communicator_pt() !=
          *x.distribution_pt()->communicator_pt())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector and the x vector must have the same communicator"
          << std::endl;
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if soln is not setup then setup the distribution
    if (!soln.built())
    {
      LinearAlgebraDistribution* dist_pt = new LinearAlgebraDistribution(
        x.distribution_pt()->communicator_pt(), this->ncol(), false);
      soln.build(dist_pt, 0.0);
      delete dist_pt;
    }

    // zero
    soln.initialise(0.0);

    // Matrix vector product
    double* soln_pt = soln.values_pt();
    const double* x_pt = x.values_pt();
    for (unsigned long i = 0; i < N; i++)
    {
      for (long k = Column_start[i]; k < Column_start[i + 1]; k++)
      {
        unsigned long j = Row_index[k];
        double a_ij = Value[k];
        soln_pt[j] += a_ij * x_pt[i];
      }
    }
  }


  //===========================================================================
  /// Function to multiply this matrix by the CCDoubleMatrix matrix_in
  /// The multiplication method used can be selected using the flag
  /// Matrix_matrix_multiply_method. By default Method 2 is used.
  /// Method 1: First runs through this matrix and matrix_in to find the storage
  ///           requirements for result - arrays of the correct size are
  ///           then allocated before performing the calculation.
  ///           Minimises memory requirements but more costly.
  /// Method 2: Grows storage for values and column indices of result 'on the
  ///           fly' using an array of maps. Faster but more memory
  ///           intensive.
  /// Method 3: Grows storage for values and column indices of result 'on the
  ///           fly' using a vector of vectors. Not particularly impressive
  ///           on the platforms we tried...
  //=============================================================================
  void CCDoubleMatrix::multiply(const CCDoubleMatrix& matrix_in,
                                CCDoubleMatrix& result)
  {
#ifdef PARANOID
    // check matrix dimensions are compatible
    if (this->ncol() != matrix_in.nrow())
    {
      std::ostringstream error_message;
      error_message
        << "Matrix dimensions incompatable for matrix-matrix multiplication"
        << "ncol() for first matrix:" << this->ncol()
        << "nrow() for second matrix: " << matrix_in.nrow();

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // NB N is number of rows!
    unsigned long N = this->nrow();
    unsigned long M = matrix_in.ncol();
    unsigned long Nnz = 0;

    // pointers to arrays which store result
    int* Column_start;
    double* Value;
    int* Row_index;

    // get pointers to matrix_in
    const int* matrix_in_col_start = matrix_in.column_start();
    const int* matrix_in_row_index = matrix_in.row_index();
    const double* matrix_in_value = matrix_in.value();

    // get pointers to this matrix
    const double* this_value = this->value();
    const int* this_col_start = this->column_start();
    const int* this_row_index = this->row_index();

    // set method
    unsigned method = Matrix_matrix_multiply_method;

    // clock_t clock1 = clock();

    // METHOD 1
    // --------
    if (method == 1)
    {
      // allocate storage for column starts
      Column_start = new int[M + 1];
      Column_start[0] = 0;

      // a set to store number of non-zero rows in each column of result
      std::set<unsigned> rows;

      // run through columns of this matrix and matrix_in to find number of
      // non-zero entries in each column of result
      for (unsigned long this_col = 0; this_col < M; this_col++)
      {
        // run through non-zeros in this_col of this matrix
        for (int this_ptr = this_col_start[this_col];
             this_ptr < this_col_start[this_col + 1];
             this_ptr++)
        {
          // find row index for non-zero
          unsigned matrix_in_col = this_row_index[this_ptr];

          // run through corresponding column in matrix_in
          for (int matrix_in_ptr = matrix_in_col_start[matrix_in_col];
               matrix_in_ptr < matrix_in_col_start[matrix_in_col + 1];
               matrix_in_ptr++)
          {
            // find row index for non-zero in matrix_in and store in rows
            rows.insert(matrix_in_row_index[matrix_in_ptr]);
          }
        }
        // update Column_start
        Column_start[this_col + 1] = Column_start[this_col] + rows.size();

        // wipe values in rows
        rows.clear();
      }

      // set Nnz
      Nnz = Column_start[M];

      // allocate arrays for result
      Value = new double[Nnz];
      Row_index = new int[Nnz];

      // set all values of Row_index to -1
      for (unsigned long i = 0; i < Nnz; i++) Row_index[i] = -1;

      // Calculate values for result - first run through columns of this matrix
      for (unsigned long this_col = 0; this_col < M; this_col++)
      {
        // run through non-zeros in this_column
        for (int this_ptr = this_col_start[this_col];
             this_ptr < this_col_start[this_col + 1];
             this_ptr++)
        {
          // find value of non-zero
          double this_val = this_value[this_ptr];

          // find row associated with non-zero
          unsigned matrix_in_col = this_row_index[this_ptr];

          // run through corresponding column in matrix_in
          for (int matrix_in_ptr = matrix_in_col_start[matrix_in_col];
               matrix_in_ptr < matrix_in_col_start[matrix_in_col + 1];
               matrix_in_ptr++)
          {
            // find row index for non-zero in matrix_in
            int row = matrix_in_row_index[matrix_in_ptr];

            // find position in result to insert value
            for (int ptr = Column_start[this_col];
                 ptr <= Column_start[this_col + 1];
                 ptr++)
            {
              if (ptr == Column_start[this_col + 1])
              {
                // error - have passed end of column without finding
                // correct row index
                std::ostringstream error_message;
                error_message << "Error inserting value in result";

                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
              else if (Row_index[ptr] == -1)
              {
                // first entry for this row index
                Row_index[ptr] = row;
                Value[ptr] = this_val * matrix_in_value[matrix_in_ptr];
                break;
              }
              else if (Row_index[ptr] == row)
              {
                // row index already exists - add value
                Value[ptr] += this_val * matrix_in_value[matrix_in_ptr];
                break;
              }
            }
          }
        }
      }
    }

    // METHOD 2
    // --------
    else if (method == 2)
    {
      // generate array of maps to store values for result
      std::map<int, double>* result_maps = new std::map<int, double>[M];

      // run through columns of this matrix
      for (unsigned long this_col = 0; this_col < M; this_col++)
      {
        // run through non-zeros in this_col
        for (int this_ptr = this_col_start[this_col];
             this_ptr < this_col_start[this_col + 1];
             this_ptr++)
        {
          // find value of non-zero
          double this_val = this_value[this_ptr];

          // find row index associated with non-zero
          unsigned matrix_in_col = this_row_index[this_ptr];

          // run through corresponding column in matrix_in
          for (int matrix_in_ptr = matrix_in_col_start[matrix_in_col];
               matrix_in_ptr < matrix_in_col_start[matrix_in_col + 1];
               matrix_in_ptr++)
          {
            // find row index for non-zero in matrix_in
            int row = matrix_in_row_index[matrix_in_ptr];

            // insert value
            result_maps[this_col][row] +=
              this_val * matrix_in_value[matrix_in_ptr];
          }
        }
      }

      // allocate Column_start
      Column_start = new int[M + 1];

      // copy across column starts
      Column_start[0] = 0;
      for (unsigned long col = 0; col < M; col++)
      {
        int size = result_maps[col].size();
        Column_start[col + 1] = Column_start[col] + size;
      }

      // set Nnz
      Nnz = Column_start[M];

      // allocate other arrays
      Value = new double[Nnz];
      Row_index = new int[Nnz];

      // copy values and row indices
      for (unsigned long col = 0; col < M; col++)
      {
        unsigned ptr = Column_start[col];
        for (std::map<int, double>::iterator i = result_maps[col].begin();
             i != result_maps[col].end();
             i++)
        {
          Row_index[ptr] = i->first;
          Value[ptr] = i->second;
          ptr++;
        }
      }

      // tidy up memory
      delete[] result_maps;
    }

    // METHOD 3
    // --------
    else if (method == 3)
    {
      // vectors of vectors to store results
      std::vector<std::vector<int>> result_rows(N);
      std::vector<std::vector<double>> result_vals(N);

      // run through the columns of this matrix
      for (unsigned long this_col = 0; this_col < M; this_col++)
      {
        // run through non-zeros in this_col
        for (int this_ptr = this_col_start[this_col];
             this_ptr < this_col_start[this_col + 1];
             this_ptr++)
        {
          // find value of non-zero
          double this_val = this_value[this_ptr];

          // find row index associated with non-zero
          unsigned matrix_in_col = this_row_index[this_ptr];

          // run through corresponding column in matrix_in
          for (int matrix_in_ptr = matrix_in_col_start[matrix_in_col];
               matrix_in_ptr < matrix_in_col_start[matrix_in_col + 1];
               matrix_in_ptr++)
          {
            // find row index for non-zero in matrix_in
            int row = matrix_in_row_index[matrix_in_ptr];

            // insert value
            int size = result_rows[this_col].size();
            for (int i = 0; i <= size; i++)
            {
              if (i == size)
              {
                // first entry for this row index
                result_rows[this_col].push_back(row);
                result_vals[this_col].push_back(this_val *
                                                matrix_in_value[matrix_in_ptr]);
              }
              else if (row == result_rows[this_col][i])
              {
                // row index already exists
                result_vals[this_col][i] +=
                  this_val * matrix_in_value[matrix_in_ptr];
                break;
              }
            }
          }
        }
      }

      // allocate Column_start
      Column_start = new int[M + 1];

      // copy across column starts
      Column_start[0] = 0;
      for (unsigned long col = 0; col < M; col++)
      {
        int size = result_rows[col].size();
        Column_start[col + 1] = Column_start[col] + size;
      }

      // set Nnz
      Nnz = Column_start[M];

      // allocate other arrays
      Value = new double[Nnz];
      Row_index = new int[Nnz];

      // copy across values and row indices
      for (unsigned long col = 0; col < N; col++)
      {
        unsigned ptr = Column_start[col];
        unsigned n_rows = result_rows[col].size();
        for (unsigned i = 0; i < n_rows; i++)
        {
          Row_index[ptr] = result_rows[col][i];
          Value[ptr] = result_vals[col][i];
          ptr++;
        }
      }
    }

    // INCORRECT VALUE FOR METHOD
    else
    {
      std::ostringstream error_message;
      error_message << "Incorrect method set in matrix-matrix multiply"
                    << "method=" << method << " not allowed";

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    result.build_without_copy(Value, Row_index, Column_start, Nnz, N, M);
  }


  //=================================================================
  /// For every row, find the maximum absolute value of the
  /// entries in this row. Set all values that are less than alpha times
  /// this maximum to zero and return the resulting matrix in
  /// reduced_matrix. Note: Diagonal entries are retained regardless
  /// of their size.
  //=================================================================
  void CCDoubleMatrix::matrix_reduction(const double& alpha,
                                        CCDoubleMatrix& reduced_matrix)
  {
    // number of columns in matrix
    long n_coln = ncol();

    Vector<double> max_row(nrow(), 0.0);

    // Here's the packed format for the new matrix
    Vector<int> B_row_start(1);
    Vector<int> B_column_index;
    Vector<double> B_value;


    // k is counter for the number of entries in the reduced matrix
    unsigned k = 0;

    // Initialise row start
    B_row_start[0] = 0;

    // Loop over columns
    for (long i = 0; i < n_coln; i++)
    {
      // Loop over entries in columns
      for (long j = Column_start[i]; j < Column_start[i + 1]; j++)
      {
        // Find max. value in row
        if (std::fabs(Value[j]) > max_row[Row_index[j]])
        {
          max_row[Row_index[j]] = std::fabs(Value[j]);
        }
      }

      // Decide if we need to retain the entries in the row
      for (long j = Column_start[i]; j < Column_start[i + 1]; j++)
      {
        // If we're on the diagonal or the value is sufficiently large: retain
        // i.e. copy across.
        if (i == Row_index[j] ||
            std::fabs(Value[j]) > alpha * max_row[Row_index[j]])
        {
          B_value.push_back(Value[j]);
          B_column_index.push_back(Row_index[j]);
          k++;
        }
      }
      // This writes the row start for the next row -- equal to
      // to the number of entries written so far (C++ zero-based indexing!)
      B_row_start.push_back(k);
    }


    // Build the matrix from the compressed format
    dynamic_cast<CCDoubleMatrix&>(reduced_matrix)
      .build(B_value, B_column_index, B_row_start, nrow(), ncol());
  }


  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  //=============================================================================
  /// Default constructor
  //=============================================================================
  CRDoubleMatrix::CRDoubleMatrix()
  {
    // set the default solver
    Linear_solver_pt = Default_linear_solver_pt = new SuperLUSolver;

    // matrix not built
    Built = false;

    // set the serial matrix-matrix multiply method
#ifdef OOMPH_HAS_TRILINOS
    //    Serial_matrix_matrix_multiply_method = 4;
    Serial_matrix_matrix_multiply_method = 2;
#else
    Serial_matrix_matrix_multiply_method = 2;
#endif
  }

  //=============================================================================
  /// Copy constructor
  //=============================================================================
  CRDoubleMatrix::CRDoubleMatrix(const CRDoubleMatrix& other_matrix)
  {
    // copy the distribution
    this->build_distribution(other_matrix.distribution_pt());

    // copy coefficients
    const double* values_pt = other_matrix.value();
    const int* column_indices = other_matrix.column_index();
    const int* row_start = other_matrix.row_start();

    // This is the local nnz.
    const unsigned nnz = other_matrix.nnz();

    // Using number of local rows since the underlying CRMatrix is local to
    // each processor.
    const unsigned nrow_local = other_matrix.nrow_local();

    // Storage for the (yet to be copied) data.
    double* my_values_pt = new double[nnz];
    int* my_column_indices = new int[nnz];
    int* my_row_start = new int[nrow_local + 1];

    // Copying over the data.
    std::copy(values_pt, values_pt + nnz, my_values_pt);
    std::copy(column_indices, column_indices + nnz, my_column_indices);
    std::copy(row_start, row_start + nrow_local + 1, my_row_start);


    // Build without copy since we have made a deep copy of the data structure.
    this->build_without_copy(
      other_matrix.ncol(), nnz, my_values_pt, my_column_indices, my_row_start);

    // set the default solver
    Linear_solver_pt = Default_linear_solver_pt = new SuperLUSolver;

    // matrix is built
    Built = true;

    // set the serial matrix-matrix multiply method
#ifdef OOMPH_HAS_TRILINOS
    // Serial_matrix_matrix_multiply_method = 4;
    Serial_matrix_matrix_multiply_method = 2;
#else
    Serial_matrix_matrix_multiply_method = 2;
#endif
  }


  //=============================================================================
  /// Constructor: just stores the distribution but does not build the
  /// matrix
  //=============================================================================
  CRDoubleMatrix::CRDoubleMatrix(
    const LinearAlgebraDistribution* distribution_pt)
  {
    this->build_distribution(distribution_pt);

    // set the default solver
    Linear_solver_pt = Default_linear_solver_pt = new SuperLUSolver;

    // matrix not built
    Built = false;

// set the serial matrix-matrix multiply method
#ifdef OOMPH_HAS_TRILINOS
    //    Serial_matrix_matrix_multiply_method = 4;
    Serial_matrix_matrix_multiply_method = 2;
#else
    Serial_matrix_matrix_multiply_method = 2;
#endif
  }

  //=============================================================================
  ///  Constructor: Takes the distribution and the number of columns, as
  /// well as the vector of values, vector of column indices,vector of row
  /// starts.
  //=============================================================================
  CRDoubleMatrix::CRDoubleMatrix(const LinearAlgebraDistribution* dist_pt,
                                 const unsigned& ncol,
                                 const Vector<double>& value,
                                 const Vector<int>& column_index,
                                 const Vector<int>& row_start)
  {
    // build the compressed row matrix
    CR_matrix.build(
      value, column_index, row_start, dist_pt->nrow_local(), ncol);

    // store the Distribution
    this->build_distribution(dist_pt);

    // set the linear solver
    Linear_solver_pt = Default_linear_solver_pt = new SuperLUSolver;

    // set the serial matrix-matrix multiply method
#ifdef OOMPH_HAS_TRILINOS
    // Serial_matrix_matrix_multiply_method = 4;
    Serial_matrix_matrix_multiply_method = 2;
#else
    Serial_matrix_matrix_multiply_method = 2;
#endif

    // matrix has been built
    Built = true;
  }


  //=============================================================================
  /// Destructor
  //=============================================================================
  CRDoubleMatrix::~CRDoubleMatrix()
  {
    this->clear();
    delete Default_linear_solver_pt;
    Default_linear_solver_pt = 0;
  }

  //=============================================================================
  /// Rebuild the matrix - assembles an empty matrix with a defined distribution
  //=============================================================================
  void CRDoubleMatrix::build(const LinearAlgebraDistribution* distribution_pt)
  {
    this->clear();
    this->build_distribution(distribution_pt);
  }

  //=============================================================================
  ///  Runs through the column index vector and checks if the entries
  /// are arranged arbitrarily or if they follow the regular lexicographical
  /// of matrices. If a boolean argument is provided with the assignment
  /// TRUE then information on the first entry which is not in the correct
  /// position will also be given
  //=============================================================================
  bool CRDoubleMatrix::entries_are_sorted(
    const bool& doc_unordered_entries) const
  {
#ifdef OOMPH_HAS_MPI
    // We only need to produce a warning if the matrix is distributed
    if (this->distributed())
    {
      // Create an ostringstream object to store the warning message
      std::ostringstream warning_message;

      // Create the warning messsage
      warning_message << "This method currently works for serial but "
                      << "has not been implemented for use with MPI!\n";

      // Issue the warning
      OomphLibWarning(warning_message.str(),
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Get the number of rows in this matrix
    unsigned n_rows = this->nrow();

    // Acquire access to the value, row_start and column_index arrays from
    // the CR matrix. Since we do not change anything in row_start_pt we
    // give it the const prefix
    const int* column_index_pt = this->column_index();
    const int* row_start_pt = this->row_start();

    // Loop over the rows of matrix
    for (unsigned i = 0; i < n_rows; i++)
    {
      // Calculate the number of nonzeros in the i-th row
      unsigned nnz_row_i = *(row_start_pt + i + 1) - *(row_start_pt + i);

      // Get the index of the first entry in row i
      unsigned i_row_start = *(row_start_pt + i);

      // Loop over the entries of the i-th row
      for (unsigned j = 0; j < nnz_row_i - 1; j++)
      {
        // Check if the column index of the following entry is greater than the
        // current entry
        if ((*(column_index_pt + i_row_start + j + 1)) <
            (*(column_index_pt + i_row_start + j)))
        {
          // If the input argument was set to TRUE we document we output
          // information of the first entry which is not in the correct position
          if (doc_unordered_entries)
          {
            // Tell the user
            oomph_info << "Matrix has not been correctly sorted!"
                       << "\nOn row: " << i << "\nEntry: " << j
                       << "\nEntry 1 = " << *(column_index_pt + i_row_start + j)
                       << "\nEntry 2 = "
                       << *(column_index_pt + i_row_start + j + 1) << std::endl;
          }

          // It hasn't worked
          return false;
        } // if ((*(column_index_pt+i_row_start+j+1)) ...
      } // for (unsigned j=0;j<nnz_row_i-1;j++)
    } // for (unsigned i=0;i<n_rows;i++)

    // If it gets here without a warning then the entries in each row of
    // the matrix are ordered by increasing column index
    return true;
  } // End of entries_are_sorted()

  //=============================================================================
  ///  This helper function sorts the entries in the column index vector
  /// and the value vector. During the construction of the matrix the entries
  /// were most likely assigned in an arbitrary order. As a result, it cannot
  /// be assumed that the entries in the column index vector corresponding to
  /// each row of the matrix have been arranged in increasing order. During
  /// the setup an additional vector will be set up; Index_of_diagonal_entries.
  /// The i-th entry of this vector contains the index of the last entry
  /// below or on the diagonal. If there are no entries below or on the
  /// diagonal then the corresponding entry is -1. If, however, there are
  /// no entries in the row then the entry is irrelevant and is kept
  /// as the initialised value; 0.
  //=============================================================================
  void CRDoubleMatrix::sort_entries()
  {
#ifdef OOMPH_HAS_MPI
    // We only need to produce a warning if the matrix is distributed
    if (this->distributed())
    {
      // Create an ostringstream object to store the warning message
      std::ostringstream warning_message;

      // Create the warning messsage
      warning_message << "This method currently works for serial but "
                      << "has not been tested with MPI!\n";

      // Issue the warning
      OomphLibWarning(warning_message.str(),
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Get the number of rows in the matrix
    unsigned n_rows = this->nrow();

    // Acquire access to the value, row_start and column_index arrays from
    // the CR matrix. Since we do not change anything in row_start_pt we
    // give it the const prefix
    double* value_pt = this->value();
    int* column_index_pt = this->column_index();
    const int* row_start_pt = this->row_start();

    // Resize the Index_of_diagonal_entries vector
    Index_of_diagonal_entries.resize(n_rows, 0);

    // Vector of pairs to store the column_index of each value in the i-th row
    // and its corresponding matrix entry
    Vector<std::pair<int, double>> column_index_and_value_row_i;

    // Loop over the rows of the matrix
    for (unsigned i = 0; i < n_rows; i++)
    {
      // Find the number of nonzeros in the i-th row
      unsigned nnz_row_i = *(row_start_pt + i + 1) - *(row_start_pt + i);

      // Variable to store the start of the i-th row
      unsigned i_row_start = *(row_start_pt + i);

      // If there are no nonzeros in this row then the i-th entry of the vector
      // Index_of_diagonal_entries is irrelevant so we can simply let it be 0
      if (nnz_row_i == 0)
      {
        // Set the i-th entry
        Index_of_diagonal_entries[i] = 0;
      }
      // If there are nonzeros in the i-th row
      else
      {
        // If there is more than one entry in the row resize the vector
        // column_index_and_value_row_i
        column_index_and_value_row_i.resize(nnz_row_i);

        // Loop over the entries in the row
        for (unsigned j = 0; j < nnz_row_i; j++)
        {
          // Assign the appropriate entries to column_index_and_value_row_i
          column_index_and_value_row_i[j] =
            std::make_pair(*(column_index_pt + i_row_start + j),
                           *(value_pt + i_row_start + j));
        }

        // Sort the vector of pairs using the struct
        // CRDoubleMatrixComparisonHelper
        std::sort(column_index_and_value_row_i.begin(),
                  column_index_and_value_row_i.end(),
                  Comparison_struct);

        //-----------------------------------------------------------------------
        // Now that the entries of the i-th row have been sorted we can read
        // them back into value_pt and column_index_pt:
        //-----------------------------------------------------------------------

        // Create a boolean variable to indicate whether or not the i-th entry
        // of Index_of_diagonal_entries has been set
        bool is_ith_entry_set = false;

        // Loop over the entries in the vector column_index_and_value_row_i and
        // assign its entries to value_pt and column_index_pt
        for (unsigned j = 0; j < nnz_row_i; j++)
        {
          // Set the column index of the j-th nonzero value in the i-th row of
          // the matrix
          *(column_index_pt + i_row_start + j) =
            column_index_and_value_row_i[j].first;

          // Set the value of the j-th nonzero value in the i-th row of
          // the matrix
          *(value_pt + i_row_start + j) =
            column_index_and_value_row_i[j].second;

          // This if statement is used to set the i-th entry of the vector
          // Index_of_diagonal_entries if it has not yet been set
          if (!is_ith_entry_set)
          {
            // If the column index of the first entry in row i is greater than
            // the row number then the first entry must lie above the diagonal
            if (unsigned(*(column_index_pt + i_row_start)) > i)
            {
              // If the column index of the first entry in the row is greater
              // than the row number, i, then the i-th entry of
              // Index_of_diagonal_entries needs to be set to -1 to indicate
              // there are no entries below or on the diagonal
              Index_of_diagonal_entries[i] = -1;

              // Indicate that the i-th entry of Index_of_diagonal_entries has
              // been set
              is_ith_entry_set = true;
            }
            // If there are entries below or on the diagonal
            else
            {
              // If there is only one entry in the row then we know that this
              // will be the last entry below or on the diagonal because we have
              // eliminated the possibility that if there is only one entry,
              // that it lies above the diagonal
              if (nnz_row_i == 1)
              {
                // Set the index of the current entry to be the value of i-th
                // entry of Index_of_diagonal_entries
                Index_of_diagonal_entries[i] = i_row_start + j;

                // Indicate that the i-th entry of Index_of_diagonal_entries has
                // been set
                is_ith_entry_set = true;
              }
              // It remains to now check the case that there is more than one
              // entry in the row. If there is more than one entry in the row
              // and there are entries below or on the diagonal then we have
              // three cases:
              //          (1) The current entry lies on the diagonal;
              //          (2) The current entry lies above the diagonal;
              //          (3) The current entry lies below the diagonal;
              // The first case can easily be checked as done below. If the
              // second case occurs then we have just passed the last entry. We
              // know this because at least one entry lies on or below the
              // diagonal. If the second case it true then we need to assign the
              // previous entry to the vector Index_of_diagonal_entries.
              // Finally, we are left with case (3), which can be split into two
              // cases:
              //          (3.1) The current entry lies below the diagonal but it
              //                is not the last entry below or on the diagonal;
              //          (3.2) The current entry lies below the diagonal and is
              //                the last entry below or on the diagonal.
              // If case (3.1) holds then we can simply wait until we get to the
              // next entry in the row and examine that. If the next entry lies
              // on the diagonal then it will be swept up by case (1). If the
              // next entry lies above the diagonal then case (2) will sweep it
              // up and if neither is the case then we wait until the next entry
              // and so on. If, instead, case (3.2) holds then our last check
              // simply needs to check if the current entry is the last entry in
              // the row because if the last entry lies on the diagonal, case
              // (1) will sweep it up. If it lies above the diagonal, case (2)
              // will take care of it. Therefore, the only remaining case is
              // that it lies strictly below the diagonal and since it is the
              // last entry in the row it means the index of this entry needs to
              // be assigned to Index_of_diagonal_entries

              // Case (1) : The current entry lies on the diagonal
              else if (unsigned(*(column_index_pt + i_row_start + j)) == i)
              {
                // Set the index of the current entry to be the value of i-th
                // entry of Index_of_diagonal_entries
                Index_of_diagonal_entries[i] = i_row_start + j;

                // Indicate that the i-th entry of Index_of_diagonal_entries has
                // been set
                is_ith_entry_set = true;
              }
              // Case (2) : The current entry lies above the diagonal
              else if (unsigned(*(column_index_pt + i_row_start + j)) > i)
              {
                // Set the index of the current entry to be the value of i-th
                // entry of Index_of_diagonal_entries
                Index_of_diagonal_entries[i] = i_row_start + j - 1;

                // Indicate that the i-th entry of Index_of_diagonal_entries has
                // been set
                is_ith_entry_set = true;
              }
              // Case (3.2) : The current entry is the last entry in the row
              else if (j == nnz_row_i - 1)
              {
                // Set the index of the current entry to be the value of i-th
                // entry of Index_of_diagonal_entries
                Index_of_diagonal_entries[i] = i_row_start + j;

                // Indicate that the i-th entry of Index_of_diagonal_entries has
                // been set
                is_ith_entry_set = true;
              } // if (nnz_row_i==1) else if
            } // if (*(column_index_pt+i_row_start)>i)
          } // if (!is_ith_entry_set)
        } // for (unsigned j=0;j<nnz_row_i;j++)
      } // if (nnz_row_i==0) else
    } // for (unsigned i=0;i<n_rows;i++)
  } // End of sort_entries()

  //=============================================================================
  /// Clean method
  //=============================================================================
  void CRDoubleMatrix::clear()
  {
    this->clear_distribution();
    CR_matrix.clean_up_memory();
    Built = false;

    if (Linear_solver_pt != 0) // Only clean up if it exists
      Linear_solver_pt->clean_up_memory();
  }

  //=============================================================================
  ///  build method: Takes the distribution and the number of columns, as
  /// well as the vector of values, vector of column indices,vector of row
  /// starts.
  //=============================================================================
  void CRDoubleMatrix::build(const LinearAlgebraDistribution* distribution_pt,
                             const unsigned& ncol,
                             const Vector<double>& value,
                             const Vector<int>& column_index,
                             const Vector<int>& row_start)
  {
    // clear
    this->clear();

    // store the Distribution
    this->build_distribution(distribution_pt);

    // set the linear solver
    Default_linear_solver_pt = new SuperLUSolver;

    // now build the matrix
    this->build(ncol, value, column_index, row_start);
  }

  //=============================================================================
  ///  method to rebuild the matrix, but not the distribution
  //=============================================================================
  void CRDoubleMatrix::build(const unsigned& ncol,
                             const Vector<double>& value,
                             const Vector<int>& column_index,
                             const Vector<int>& row_start)
  {
    // call the underlying build method
    CR_matrix.clean_up_memory();
    CR_matrix.build(value, column_index, row_start, this->nrow_local(), ncol);

    // matrix has been build
    Built = true;
  }

  //=============================================================================
  ///  method to rebuild the matrix, but not the distribution
  //=============================================================================
  void CRDoubleMatrix::build_without_copy(const unsigned& ncol,
                                          const unsigned& nnz,
                                          double* value,
                                          int* column_index,
                                          int* row_start)
  {
    // call the underlying build method
    CR_matrix.clean_up_memory();
    CR_matrix.build_without_copy(
      value, column_index, row_start, nnz, this->nrow_local(), ncol);

    // matrix has been build
    Built = true;
  }

  //=============================================================================
  /// Do LU decomposition
  //=============================================================================
  void CRDoubleMatrix::ludecompose()
  {
#ifdef PARANOID
    // check that the this matrix is built
    if (!Built)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "This matrix has not been built.";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // factorise using superlu or superlu dist if we oomph has mpi
    static_cast<SuperLUSolver*>(Default_linear_solver_pt)->factorise(this);
  }

  //=============================================================================
  /// Do back-substitution
  //=============================================================================
  void CRDoubleMatrix::lubksub(DoubleVector& rhs)
  {
#ifdef PARANOID
    // check that the rhs vector is setup
    if (!rhs.built())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The vector rhs has not been setup";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // check that the rhs vector has the same distribution as this matrix
    if (!(*this->distribution_pt() == *rhs.distribution_pt()))
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "The vector rhs must have the same distribution as the matrix";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // backsub
    DoubleVector rhs_copy(rhs);
    static_cast<SuperLUSolver*>(Default_linear_solver_pt)
      ->backsub(rhs_copy, rhs);
  }

  //=============================================================================
  ///  Multiply the matrix by the vector x
  //=============================================================================
  void CRDoubleMatrix::multiply(const DoubleVector& x, DoubleVector& soln) const
  {
#ifdef PARANOID
    // check that this matrix is built
    if (!Built)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "This matrix has not been built";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // check that the distribution of x is setup
    if (!x.built())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The distribution of the vector x must be setup";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // Check to see if x.size() = ncol().
    if (this->ncol() != x.distribution_pt()->nrow())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The number of rows in the x vector and the "
                              "number of columns in the "
                           << "matrix must be the same";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // if the soln is distributed
    if (soln.built())
    {
      if (!(*soln.distribution_pt() == *this->distribution_pt()))
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector is setup and therefore must have the same "
          << "distribution as the matrix";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if soln is not setup then setup the distribution
    if (!soln.built())
    {
      // Resize and initialize the solution vector
      soln.build(this->distribution_pt(), 0.0);
    }

    // Initialise
    soln.initialise(0.0);

    // if distributed and on more than one processor use trilinos
    // otherwise use the oomph-lib methods
    if (this->distributed() &&
        this->distribution_pt()->communicator_pt()->nproc() > 1)
    {
#ifdef OOMPH_HAS_TRILINOS
      // This will only work if we have trilinos on board
      TrilinosEpetraHelpers::multiply(this, x, soln);
#else
      std::ostringstream error_message_stream;
      error_message_stream
        << "Matrix-vector product on multiple processors with distributed "
        << "CRDoubleMatrix requires Trilinos.";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
#endif
    }
    else
    {
      unsigned n = this->nrow();
      const int* row_start = CR_matrix.row_start();
      const int* column_index = CR_matrix.column_index();
      const double* value = CR_matrix.value();
      double* soln_pt = soln.values_pt();
      const double* x_pt = x.values_pt();
      for (unsigned long i = 0; i < n; i++)
      {
        soln_pt[i] = 0.0;
        for (long k = row_start[i]; k < row_start[i + 1]; k++)
        {
          unsigned long j = column_index[k];
          double a_ij = value[k];
          soln_pt[i] += a_ij * x_pt[j];
        }
      }
    }
  }

  //=================================================================
  /// Multiply the transposed matrix by the vector x: soln=A^T x
  //=================================================================
  void CRDoubleMatrix::multiply_transpose(const DoubleVector& x,
                                          DoubleVector& soln) const
  {
#ifdef PARANOID
    // check that this matrix is built
    if (!Built)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "This matrix has not been built";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // Check to see if x.size() = ncol().
    if (!(*this->distribution_pt() == *x.distribution_pt()))
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "The x vector and this matrix must have the same distribution.";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // if soln is setup then it should have the same distribution as x
    if (soln.built())
    {
      if (soln.distribution_pt()->nrow() != this->ncol())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector is setup and therefore must have the same "
          << "number of rows as the vector x";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if soln is not setup then setup the distribution
    if (!soln.built())
    {
      LinearAlgebraDistribution* dist_pt =
        new LinearAlgebraDistribution(x.distribution_pt()->communicator_pt(),
                                      this->ncol(),
                                      this->distributed());
      soln.build(dist_pt, 0.0);
      delete dist_pt;
    }

    // Initialise
    soln.initialise(0.0);

    if (this->distributed() &&
        this->distribution_pt()->communicator_pt()->nproc() > 1)
    {
#ifdef OOMPH_HAS_TRILINOS
      // This will only work if we have trilinos on board
      TrilinosEpetraHelpers::multiply(this, x, soln);
#else
      std::ostringstream error_message_stream;
      error_message_stream
        << "Matrix-vector product on multiple processors with distributed "
        << "CRDoubleMatrix requires Trilinos.";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
#endif
    }
    else
    {
      unsigned n = this->nrow();
      const int* row_start = CR_matrix.row_start();
      const int* column_index = CR_matrix.column_index();
      const double* value = CR_matrix.value();
      double* soln_pt = soln.values_pt();
      const double* x_pt = x.values_pt();
      // Matrix vector product
      for (unsigned long i = 0; i < n; i++)
      {
        for (long k = row_start[i]; k < row_start[i + 1]; k++)
        {
          unsigned long j = column_index[k];
          double a_ij = value[k];
          soln_pt[j] += a_ij * x_pt[i];
        }
      }
    }
  }

  //===========================================================================
  /// Function to multiply this matrix by the CRDoubleMatrix matrix_in.
  /// In a serial matrix, there are 4 methods available:
  /// Method 1: First runs through this matrix and matrix_in to find the storage
  ///           requirements for result - arrays of the correct size are
  ///           then allocated before performing the calculation.
  ///           Minimises memory requirements but more costly.
  /// Method 2: Grows storage for values and column indices of result 'on the
  ///           fly' using an array of maps. Faster but more memory
  ///           intensive.
  /// Method 3: Grows storage for values and column indices of result 'on the
  ///           fly' using a vector of vectors. Not particularly impressive
  ///           on the platforms we tried...
  /// Method 4: Trilinos Epetra Matrix Matrix multiply.
  /// Method 5: Trilinox Epetra Matrix Matrix Mulitply (ml based)
  /// If Trilinos is installed then Method 4 is employed by default, otherwise
  /// Method 2 is employed by default.
  /// In a distributed matrix, only Trilinos Epetra Matrix Matrix multiply
  /// is available.
  //=============================================================================
  void CRDoubleMatrix::multiply(const CRDoubleMatrix& matrix_in,
                                CRDoubleMatrix& result) const
  {
#ifdef PARANOID
    // check that this matrix is built
    if (!Built)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "This matrix has not been built";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // check that this matrix is built
    if (!matrix_in.built())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "This matrix matrix_in has not been built";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // if soln is setup then it should have the same distribution as x
    if (result.built())
    {
      if (!(*result.distribution_pt() == *this->distribution_pt()))
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The matrix result is setup and therefore must have the same "
          << "distribution as the vector x";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if the result has not been setup, then store the distribution
    if (!result.distribution_built())
    {
      result.build(this->distribution_pt());
    }

    // short name for Serial_matrix_matrix_multiply_method
    unsigned method = Serial_matrix_matrix_multiply_method;

    // if this matrix is not distributed and matrix in is not distributed
    if (!this->distributed() && !matrix_in.distributed() &&
        ((method == 1) || (method == 2) || (method == 3)))
    {
      // NB N is number of rows!
      unsigned long N = this->nrow();
      unsigned long M = matrix_in.ncol();
      unsigned long Nnz = 0;

      // pointers to arrays which store result
      int* Row_start = 0;
      double* Value = 0;
      int* Column_index = 0;

      // get pointers to matrix_in
      const int* matrix_in_row_start = matrix_in.row_start();
      const int* matrix_in_column_index = matrix_in.column_index();
      const double* matrix_in_value = matrix_in.value();

      // get pointers to this matrix
      const double* this_value = this->value();
      const int* this_row_start = this->row_start();
      const int* this_column_index = this->column_index();

      // clock_t clock1 = clock();

      // METHOD 1
      // --------
      if (method == 1)
      {
        // allocate storage for row starts
        Row_start = new int[N + 1];
        Row_start[0] = 0;

        // a set to store number of non-zero columns in each row of result
        std::set<unsigned> columns;

        // run through rows of this matrix and matrix_in to find number of
        // non-zero entries in each row of result
        for (unsigned long this_row = 0; this_row < N; this_row++)
        {
          // run through non-zeros in this_row of this matrix
          for (int this_ptr = this_row_start[this_row];
               this_ptr < this_row_start[this_row + 1];
               this_ptr++)
          {
            // find column index for non-zero
            int matrix_in_row = this_column_index[this_ptr];

            // run through corresponding row in matrix_in
            for (int matrix_in_ptr = matrix_in_row_start[matrix_in_row];
                 matrix_in_ptr < matrix_in_row_start[matrix_in_row + 1];
                 matrix_in_ptr++)
            {
              // find column index for non-zero in matrix_in and store in
              // columns
              columns.insert(matrix_in_column_index[matrix_in_ptr]);
            }
          }
          // update Row_start
          Row_start[this_row + 1] = Row_start[this_row] + columns.size();

          // wipe values in columns
          columns.clear();
        }

        // set Nnz
        Nnz = Row_start[N];

        // allocate arrays for result
        Value = new double[Nnz];
        Column_index = new int[Nnz];

        // set all values of Column_index to -1
        for (unsigned long i = 0; i < Nnz; i++)
        {
          Column_index[i] = -1;
        }

        // Calculate values for result - first run through rows of this matrix
        for (unsigned long this_row = 0; this_row < N; this_row++)
        {
          // run through non-zeros in this_row
          for (int this_ptr = this_row_start[this_row];
               this_ptr < this_row_start[this_row + 1];
               this_ptr++)
          {
            // find value of non-zero
            double this_val = this_value[this_ptr];

            // find column associated with non-zero
            int matrix_in_row = this_column_index[this_ptr];

            // run through corresponding row in matrix_in
            for (int matrix_in_ptr = matrix_in_row_start[matrix_in_row];
                 matrix_in_ptr < matrix_in_row_start[matrix_in_row + 1];
                 matrix_in_ptr++)
            {
              // find column index for non-zero in matrix_in
              int col = matrix_in_column_index[matrix_in_ptr];

              // find position in result to insert value
              for (int ptr = Row_start[this_row];
                   ptr <= Row_start[this_row + 1];
                   ptr++)
              {
                if (ptr == Row_start[this_row + 1])
                {
                  // error - have passed end of row without finding
                  // correct column
                  std::ostringstream error_message;
                  error_message << "Error inserting value in result";

                  throw OomphLibError(error_message.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
                else if (Column_index[ptr] == -1)
                {
                  // first entry for this column index
                  Column_index[ptr] = col;
                  Value[ptr] = this_val * matrix_in_value[matrix_in_ptr];
                  break;
                }
                else if (Column_index[ptr] == col)
                {
                  // column index already exists - add value
                  Value[ptr] += this_val * matrix_in_value[matrix_in_ptr];
                  break;
                }
              }
            }
          }
        }
      }

      // METHOD 2
      // --------
      else if (method == 2)
      {
        // generate array of maps to store values for result
        std::map<int, double>* result_maps = new std::map<int, double>[N];

        // run through rows of this matrix
        for (unsigned long this_row = 0; this_row < N; this_row++)
        {
          // run through non-zeros in this_row
          for (int this_ptr = this_row_start[this_row];
               this_ptr < this_row_start[this_row + 1];
               this_ptr++)
          {
            // find value of non-zero
            double this_val = this_value[this_ptr];

            // find column index associated with non-zero
            int matrix_in_row = this_column_index[this_ptr];

            // run through corresponding row in matrix_in
            for (int matrix_in_ptr = matrix_in_row_start[matrix_in_row];
                 matrix_in_ptr < matrix_in_row_start[matrix_in_row + 1];
                 matrix_in_ptr++)
            {
              // find column index for non-zero in matrix_in
              int col = matrix_in_column_index[matrix_in_ptr];

              // insert value
              result_maps[this_row][col] +=
                this_val * matrix_in_value[matrix_in_ptr];
            }
          }
        }

        // allocate Row_start
        Row_start = new int[N + 1];

        // copy across row starts
        Row_start[0] = 0;
        for (unsigned long row = 0; row < N; row++)
        {
          int size = result_maps[row].size();
          Row_start[row + 1] = Row_start[row] + size;
        }

        // set Nnz
        Nnz = Row_start[N];

        // allocate other arrays
        Value = new double[Nnz];
        Column_index = new int[Nnz];

        // copy values and column indices
        for (unsigned long row = 0; row < N; row++)
        {
          unsigned ptr = Row_start[row];
          for (std::map<int, double>::iterator i = result_maps[row].begin();
               i != result_maps[row].end();
               i++)
          {
            Column_index[ptr] = i->first;
            Value[ptr] = i->second;
            ptr++;
          }
        }

        // tidy up memory
        delete[] result_maps;
      }

      // METHOD 3
      // --------
      else if (method == 3)
      {
        // vectors of vectors to store results
        std::vector<std::vector<int>> result_cols(N);
        std::vector<std::vector<double>> result_vals(N);

        // run through the rows of this matrix
        for (unsigned long this_row = 0; this_row < N; this_row++)
        {
          // run through non-zeros in this_row
          for (int this_ptr = this_row_start[this_row];
               this_ptr < this_row_start[this_row + 1];
               this_ptr++)
          {
            // find value of non-zero
            double this_val = this_value[this_ptr];

            // find column index associated with non-zero
            int matrix_in_row = this_column_index[this_ptr];

            // run through corresponding row in matrix_in
            for (int matrix_in_ptr = matrix_in_row_start[matrix_in_row];
                 matrix_in_ptr < matrix_in_row_start[matrix_in_row + 1];
                 matrix_in_ptr++)
            {
              // find column index for non-zero in matrix_in
              int col = matrix_in_column_index[matrix_in_ptr];

              // insert value
              int size = result_cols[this_row].size();
              for (int i = 0; i <= size; i++)
              {
                if (i == size)
                {
                  // first entry for this column
                  result_cols[this_row].push_back(col);
                  result_vals[this_row].push_back(
                    this_val * matrix_in_value[matrix_in_ptr]);
                }
                else if (col == result_cols[this_row][i])
                {
                  // column already exists
                  result_vals[this_row][i] +=
                    this_val * matrix_in_value[matrix_in_ptr];
                  break;
                }
              }
            }
          }
        }

        // allocate Row_start
        Row_start = new int[N + 1];

        // copy across row starts
        Row_start[0] = 0;
        for (unsigned long row = 0; row < N; row++)
        {
          int size = result_cols[row].size();
          Row_start[row + 1] = Row_start[row] + size;
        }

        // set Nnz
        Nnz = Row_start[N];

        // allocate other arrays
        Value = new double[Nnz];
        Column_index = new int[Nnz];

        // copy across values and column indices
        for (unsigned long row = 0; row < N; row++)
        {
          unsigned ptr = Row_start[row];
          unsigned nnn = result_cols[row].size();
          for (unsigned i = 0; i < nnn; i++)
          {
            Column_index[ptr] = result_cols[row][i];
            Value[ptr] = result_vals[row][i];
            ptr++;
          }
        }
      }

      // build
      result.build_without_copy(M, Nnz, Value, Column_index, Row_start);
    }

    // else we have to use trilinos
    else
    {
#ifdef OOMPH_HAS_TRILINOS
      bool use_ml = false;
      if (method == 5)
      {
        use_ml = true;
      }
      TrilinosEpetraHelpers::multiply(*this, matrix_in, result, use_ml);
#else
      std::ostringstream error_message;
      error_message << "Serial_matrix_matrix_multiply_method = "
                    << Serial_matrix_matrix_multiply_method
                    << " requires trilinos.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    }
  }


  //=================================================================
  /// For every row, find the maximum absolute value of the
  /// entries in this row. Set all values that are less than alpha times
  /// this maximum to zero and return the resulting matrix in
  /// reduced_matrix. Note: Diagonal entries are retained regardless
  /// of their size.
  //=================================================================
  void CRDoubleMatrix::matrix_reduction(const double& alpha,
                                        CRDoubleMatrix& reduced_matrix)
  {
    // number of rows in matrix
    long n_row = nrow_local();
    double max_row;

    // Here's the packed format for the new matrix
    Vector<int> B_row_start(1);
    Vector<int> B_column_index;
    Vector<double> B_value;

    // get pointers to the underlying data
    const int* row_start = CR_matrix.row_start();
    const int* column_index = CR_matrix.column_index();
    const double* value = CR_matrix.value();

    // k is counter for the number of entries in the reduced matrix
    unsigned k = 0;

    // Initialise row start
    B_row_start[0] = 0;

    // Loop over rows
    for (long i = 0; i < n_row; i++)
    {
      // Initialise max value in row
      max_row = 0.0;

      // Loop over entries in columns
      for (long j = row_start[i]; j < row_start[i + 1]; j++)
      {
        // Find max. value in row
        if (std::fabs(value[j]) > max_row)
        {
          max_row = std::fabs(value[j]);
        }
      }

      // Decide if we need to retain the entries in the row
      for (long j = row_start[i]; j < row_start[i + 1]; j++)
      {
        // If we're on the diagonal or the value is sufficiently large: retain
        // i.e. copy across.
        if (i == column_index[j] || std::fabs(value[j]) > alpha * max_row)
        {
          B_value.push_back(value[j]);
          B_column_index.push_back(column_index[j]);
          k++;
        }
      }
      // This writes the row start for the next row -- equal to
      // to the number of entries written so far (C++ zero-based indexing!)
      B_row_start.push_back(k);
    }

    // Build the matrix from the compressed format
    dynamic_cast<CRDoubleMatrix&>(reduced_matrix)
      .build(this->ncol(), B_value, B_column_index, B_row_start);
  }

  //=============================================================================
  /// if this matrix is distributed then the equivalent global matrix is built
  /// using new and returned. The calling method is responsible for the
  /// destruction of the new matrix.
  //=============================================================================
  CRDoubleMatrix* CRDoubleMatrix::global_matrix() const
  {
#ifdef OOMPH_HAS_MPI
    // if this matrix is not distributed then this method is redundant
    if (!this->distributed() ||
        this->distribution_pt()->communicator_pt()->nproc() == 1)
    {
      return new CRDoubleMatrix(*this);
    }

    // nnz
    int nnz = this->nnz();

    // my nrow local
    unsigned nrow_local = this->nrow_local();

    // nrow global
    unsigned nrow = this->nrow();

    // cache nproc
    int nproc = this->distribution_pt()->communicator_pt()->nproc();

    // get the nnzs on the other processors
    int* dist_nnz_pt = new int[nproc];
    MPI_Allgather(&nnz,
                  1,
                  MPI_INT,
                  dist_nnz_pt,
                  1,
                  MPI_INT,
                  this->distribution_pt()->communicator_pt()->mpi_comm());

    // create a int vector of first rows and nrow local and compute nnz global
    int* dist_first_row = new int[nproc];
    int* dist_nrow_local = new int[nproc];
    int nnz_global = 0;
    for (int p = 0; p < nproc; p++)
    {
      nnz_global += dist_nnz_pt[p];
      dist_first_row[p] = this->first_row(p);
      dist_nrow_local[p] = this->nrow_local(p);
    }

    // compute the offset for the values and column index data
    int* nnz_offset = new int[nproc];
    nnz_offset[0] = 0;
    for (int p = 1; p < nproc; p++)
    {
      nnz_offset[p] = nnz_offset[p - 1] + dist_nnz_pt[p - 1];
    }

    // get pointers to the (current) distributed data
    // const_cast required because MPI requires non-const data when sending
    // data
    int* dist_row_start = const_cast<int*>(this->row_start());
    int* dist_column_index = const_cast<int*>(this->column_index());
    double* dist_value = const_cast<double*>(this->value());

    // space for the global matrix
    int* global_row_start = new int[nrow + 1];
    int* global_column_index = new int[nnz_global];
    double* global_value = new double[nnz_global];

    // get the row starts
    MPI_Allgatherv(dist_row_start,
                   nrow_local,
                   MPI_INT,
                   global_row_start,
                   dist_nrow_local,
                   dist_first_row,
                   MPI_INT,
                   this->distribution_pt()->communicator_pt()->mpi_comm());

    // get the column indexes
    MPI_Allgatherv(dist_column_index,
                   nnz,
                   MPI_INT,
                   global_column_index,
                   dist_nnz_pt,
                   nnz_offset,
                   MPI_INT,
                   this->distribution_pt()->communicator_pt()->mpi_comm());

    // get the values
    MPI_Allgatherv(dist_value,
                   nnz,
                   MPI_DOUBLE,
                   global_value,
                   dist_nnz_pt,
                   nnz_offset,
                   MPI_DOUBLE,
                   this->distribution_pt()->communicator_pt()->mpi_comm());

    // finally the last row start
    global_row_start[nrow] = nnz_global;

    // update the other row start
    for (int p = 0; p < nproc; p++)
    {
      for (int i = 0; i < dist_nrow_local[p]; i++)
      {
        unsigned j = dist_first_row[p] + i;
        global_row_start[j] += nnz_offset[p];
      }
    }

    // create the global distribution
    LinearAlgebraDistribution* dist_pt = new LinearAlgebraDistribution(
      this->distribution_pt()->communicator_pt(), nrow, false);

    // create the matrix
    CRDoubleMatrix* matrix_pt = new CRDoubleMatrix(dist_pt);

    // copy of distribution taken so delete
    delete dist_pt;

    // pass data into matrix
    matrix_pt->build_without_copy(this->ncol(),
                                  nnz_global,
                                  global_value,
                                  global_column_index,
                                  global_row_start);

    // clean up
    delete[] dist_first_row;
    delete[] dist_nrow_local;
    delete[] nnz_offset;
    delete[] dist_nnz_pt;

    // and return
    return matrix_pt;
#else
    return new CRDoubleMatrix(*this);
#endif
  }

  //============================================================================
  /// The contents of the matrix are redistributed to match the new
  /// distribution. In a non-MPI build this method does nothing.
  /// \b NOTE 1: The current distribution and the new distribution must have
  /// the same number of global rows.
  /// \b NOTE 2: The current distribution and the new distribution must have
  /// the same Communicator.
  //============================================================================
  void CRDoubleMatrix::redistribute(
    const LinearAlgebraDistribution* const& dist_pt)
  {
#ifdef OOMPH_HAS_MPI
#ifdef PARANOID
    // paranoid check that the nrows for both distributions is the
    // same
    if (dist_pt->nrow() != this->distribution_pt()->nrow())
    {
      std::ostringstream error_message;
      error_message << "The number of global rows in the new distribution ("
                    << dist_pt->nrow() << ") is not equal to the number"
                    << " of global rows in the current distribution ("
                    << this->distribution_pt()->nrow() << ").\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // paranoid check that the current distribution and the new distribution
    // have the same Communicator
    OomphCommunicator temp_comm(*dist_pt->communicator_pt());
    if (!(temp_comm == *this->distribution_pt()->communicator_pt()))
    {
      std::ostringstream error_message;
      error_message << "The new distribution and the current distribution must "
                    << "have the same communicator.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // paranoid check that the matrix is build
    if (!this->built())
    {
      std::ostringstream error_message;
      error_message << "The matrix must be build to be redistributed";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // if the two distributions are not the same
    // =========================================
    if (!((*this->distribution_pt()) == *dist_pt))
    {
      // Get the number of columns to build the matrix.
      unsigned long ncol = this->ncol();

      // current data
      int* current_row_start = this->row_start();
      int* current_column_index = this->column_index();
      double* current_value = this->value();

      // get the rank and the number of processors
      int my_rank = this->distribution_pt()->communicator_pt()->my_rank();
      int nproc = this->distribution_pt()->communicator_pt()->nproc();

      // if both distributions are distributed
      // =====================================
      if (this->distributed() && dist_pt->distributed())
      {
        // new nrow_local and first_row data
        Vector<unsigned> new_first_row(nproc);
        Vector<unsigned> new_nrow_local(nproc);
        Vector<unsigned> current_first_row(nproc);
        Vector<unsigned> current_nrow_local(nproc);
        for (int i = 0; i < nproc; i++)
        {
          new_first_row[i] = dist_pt->first_row(i);
          new_nrow_local[i] = dist_pt->nrow_local(i);
          current_first_row[i] = this->first_row(i);
          current_nrow_local[i] = this->nrow_local(i);
        }

        // compute which local rows are expected to be received from each
        // processor / sent to each processor
        Vector<unsigned> first_row_for_proc(nproc, 0);
        Vector<unsigned> nrow_local_for_proc(nproc, 0);
        Vector<unsigned> first_row_from_proc(nproc, 0);
        Vector<unsigned> nrow_local_from_proc(nproc, 0);

        // for every processor compute first_row and nrow_local that will
        // will sent and received by this processor
        for (int p = 0; p < nproc; p++)
        {
          // start with data to be sent
          if ((new_first_row[p] <
               (current_first_row[my_rank] + current_nrow_local[my_rank])) &&
              (current_first_row[my_rank] <
               (new_first_row[p] + new_nrow_local[p])))
          {
            first_row_for_proc[p] =
              std::max(current_first_row[my_rank], new_first_row[p]);
            nrow_local_for_proc[p] =
              std::min(
                (current_first_row[my_rank] + current_nrow_local[my_rank]),
                (new_first_row[p] + new_nrow_local[p])) -
              first_row_for_proc[p];
          }

          // and data to be received
          if ((new_first_row[my_rank] <
               (current_first_row[p] + current_nrow_local[p])) &&
              (current_first_row[p] <
               (new_first_row[my_rank] + new_nrow_local[my_rank])))
          {
            first_row_from_proc[p] =
              std::max(current_first_row[p], new_first_row[my_rank]);
            nrow_local_from_proc[p] =
              std::min((current_first_row[p] + current_nrow_local[p]),
                       (new_first_row[my_rank] + new_nrow_local[my_rank])) -
              first_row_from_proc[p];
          }
        }

        // determine how many nnzs to send to each processor
        Vector<unsigned> nnz_for_proc(nproc, 0);
        for (int p = 0; p < nproc; p++)
        {
          if (nrow_local_for_proc[p] > 0)
          {
            nnz_for_proc[p] = (current_row_start[first_row_for_proc[p] -
                                                 current_first_row[my_rank] +
                                                 nrow_local_for_proc[p]] -
                               current_row_start[first_row_for_proc[p] -
                                                 current_first_row[my_rank]]);
          }
        }

        // next post non-blocking sends and recv for the nnzs
        Vector<unsigned> nnz_from_proc(nproc, 0);
        Vector<MPI_Request> send_req;
        Vector<MPI_Request> nnz_recv_req;
        for (int p = 0; p < nproc; p++)
        {
          if (p != my_rank)
          {
            // send
            if (nrow_local_for_proc[p] > 0)
            {
              MPI_Request req;
              MPI_Isend(&nnz_for_proc[p],
                        1,
                        MPI_UNSIGNED,
                        p,
                        0,
                        this->distribution_pt()->communicator_pt()->mpi_comm(),
                        &req);
              send_req.push_back(req);
            }

            // recv
            if (nrow_local_from_proc[p] > 0)
            {
              MPI_Request req;
              MPI_Irecv(&nnz_from_proc[p],
                        1,
                        MPI_UNSIGNED,
                        p,
                        0,
                        this->distribution_pt()->communicator_pt()->mpi_comm(),
                        &req);
              nnz_recv_req.push_back(req);
            }
          }
          // "send to self"
          else
          {
            nnz_from_proc[p] = nnz_for_proc[p];
          }
        }

        // allocate new storage for the new row_start
        int* new_row_start = new int[new_nrow_local[my_rank] + 1];

        // wait for recvs to complete
        unsigned n_recv_req = nnz_recv_req.size();
        if (n_recv_req > 0)
        {
          Vector<MPI_Status> recv_status(n_recv_req);
          MPI_Waitall(n_recv_req, &nnz_recv_req[0], &recv_status[0]);
        }

        // compute the nnz offset for each processor
        unsigned next_row = 0;
        unsigned nnz_count = 0;
        Vector<unsigned> nnz_offset(nproc, 0);
        for (int p = 0; p < nproc; p++)
        {
          unsigned pp = 0;
          while (new_first_row[pp] != next_row)
          {
            pp++;
          }
          nnz_offset[pp] = nnz_count;
          nnz_count += nnz_from_proc[pp];
          next_row += new_nrow_local[pp];
        }

        // allocate storage for the values and column indices
        int* new_column_index = new int[nnz_count];
        double* new_value = new double[nnz_count];

        // post the sends and recvs for the matrix data
        Vector<MPI_Request> recv_req;
        MPI_Aint base_address;
        MPI_Get_address(new_value, &base_address);
        for (int p = 0; p < nproc; p++)
        {
          // communicated with other processors
          if (p != my_rank)
          {
            // SEND
            if (nrow_local_for_proc[p] > 0)
            {
              // array of datatypes
              MPI_Datatype types[3];

              // array of offsets
              MPI_Aint offsets[3];

              // array of lengths
              int len[3];

              // row start
              unsigned first_row_to_send =
                first_row_for_proc[p] - current_first_row[my_rank];
              MPI_Type_contiguous(nrow_local_for_proc[p], MPI_INT, &types[0]);
              MPI_Type_commit(&types[0]);
              len[0] = 1;
              MPI_Get_address(current_row_start + first_row_to_send,
                              &offsets[0]);
              offsets[0] -= base_address;

              // values
              unsigned first_coef_to_send =
                current_row_start[first_row_to_send];
              MPI_Type_contiguous(nnz_for_proc[p], MPI_DOUBLE, &types[1]);
              MPI_Type_commit(&types[1]);
              len[1] = 1;
              MPI_Get_address(current_value + first_coef_to_send, &offsets[1]);
              offsets[1] -= base_address;

              // column index
              MPI_Type_contiguous(nnz_for_proc[p], MPI_DOUBLE, &types[2]);
              MPI_Type_commit(&types[2]);
              len[2] = 1;
              MPI_Get_address(current_column_index + first_coef_to_send,
                              &offsets[2]);
              offsets[2] -= base_address;

              // build the combined datatype
              MPI_Datatype send_type;
              MPI_Type_create_struct(3, len, offsets, types, &send_type);
              MPI_Type_commit(&send_type);
              MPI_Type_free(&types[0]);
              MPI_Type_free(&types[1]);
              MPI_Type_free(&types[2]);

              // and send
              MPI_Request req;
              MPI_Isend(new_value,
                        1,
                        send_type,
                        p,
                        1,
                        this->distribution_pt()->communicator_pt()->mpi_comm(),
                        &req);
              send_req.push_back(req);
              MPI_Type_free(&send_type);
            }

            // RECV
            if (nrow_local_from_proc[p] > 0)
            {
              // array of datatypes
              MPI_Datatype types[3];

              // array of offsets
              MPI_Aint offsets[3];

              // array of lengths
              int len[3];

              // row start
              unsigned first_row_to_recv =
                first_row_from_proc[p] - new_first_row[my_rank];
              MPI_Type_contiguous(nrow_local_from_proc[p], MPI_INT, &types[0]);
              MPI_Type_commit(&types[0]);
              len[0] = 1;
              MPI_Get_address(new_row_start + first_row_to_recv, &offsets[0]);
              offsets[0] -= base_address;

              // values
              unsigned first_coef_to_recv = nnz_offset[p];
              MPI_Type_contiguous(nnz_from_proc[p], MPI_DOUBLE, &types[1]);
              MPI_Type_commit(&types[1]);
              len[1] = 1;
              MPI_Get_address(new_value + first_coef_to_recv, &offsets[1]);
              offsets[1] -= base_address;

              // column index
              MPI_Type_contiguous(nnz_from_proc[p], MPI_INT, &types[2]);
              MPI_Type_commit(&types[2]);
              len[2] = 1;
              MPI_Get_address(new_column_index + first_coef_to_recv,
                              &offsets[2]);
              offsets[2] -= base_address;

              // build the combined datatype
              MPI_Datatype recv_type;
              MPI_Type_create_struct(3, len, offsets, types, &recv_type);
              MPI_Type_commit(&recv_type);
              MPI_Type_free(&types[0]);
              MPI_Type_free(&types[1]);
              MPI_Type_free(&types[2]);

              // and send
              MPI_Request req;
              MPI_Irecv(new_value,
                        1,
                        recv_type,
                        p,
                        1,
                        this->distribution_pt()->communicator_pt()->mpi_comm(),
                        &req);
              recv_req.push_back(req);
              MPI_Type_free(&recv_type);
            }
          }
          // other wise transfer data internally
          else
          {
            unsigned j =
              first_row_for_proc[my_rank] - current_first_row[my_rank];
            unsigned k = first_row_from_proc[my_rank] - new_first_row[my_rank];
            for (unsigned i = 0; i < nrow_local_for_proc[my_rank]; i++)
            {
              new_row_start[k + i] = current_row_start[j + i];
            }
            unsigned first_coef_to_send = current_row_start[j];
            for (unsigned i = 0; i < nnz_for_proc[my_rank]; i++)
            {
              new_value[nnz_offset[p] + i] =
                current_value[first_coef_to_send + i];
              new_column_index[nnz_offset[p] + i] =
                current_column_index[first_coef_to_send + i];
            }
          }
        }

        // wait for all recvs to complete
        n_recv_req = recv_req.size();
        if (n_recv_req > 0)
        {
          Vector<MPI_Status> recv_status(n_recv_req);
          MPI_Waitall(n_recv_req, &recv_req[0], &recv_status[0]);
        }

        // next we need to update the row starts
        for (int p = 0; p < nproc; p++)
        {
          if (nrow_local_from_proc[p] > 0)
          {
            unsigned first_row =
              first_row_from_proc[p] - new_first_row[my_rank];
            unsigned last_row = first_row + nrow_local_from_proc[p] - 1;
            int update = nnz_offset[p] - new_row_start[first_row];
            for (unsigned i = first_row; i <= last_row; i++)
            {
              new_row_start[i] += update;
            }
          }
        }
        new_row_start[dist_pt->nrow_local()] = nnz_count;

        // wait for sends to complete
        unsigned n_send_req = send_req.size();
        if (n_recv_req > 0)
        {
          Vector<MPI_Status> send_status(n_recv_req);
          MPI_Waitall(n_send_req, &send_req[0], &send_status[0]);
        }
        // if (my_rank == 0)
        //  {
        //   CRDoubleMatrix* m_pt = this->global_matrix();
        //   m_pt->sparse_indexed_output("m1.dat");
        //  }

        //
        this->build(dist_pt);
        this->build_without_copy(
          ncol, nnz_count, new_value, new_column_index, new_row_start);
        // if (my_rank == 0)
        //  {
        //   CRDoubleMatrix* m_pt = this->global_matrix();
        //   m_pt->sparse_indexed_output("m2.dat");
        //  }
        //       this->sparse_indexed_output(oomph_info);
        abort();
      }

      // if this matrix is distributed but the new distributed matrix is global
      // ======================================================================
      else if (this->distributed() && !dist_pt->distributed())
      {
        // nnz
        int nnz = this->nnz();

        // nrow global
        unsigned nrow = this->nrow();

        // cache nproc
        int nproc = this->distribution_pt()->communicator_pt()->nproc();

        // get the nnzs on the other processors
        int* dist_nnz_pt = new int[nproc];
        MPI_Allgather(&nnz,
                      1,
                      MPI_INT,
                      dist_nnz_pt,
                      1,
                      MPI_INT,
                      this->distribution_pt()->communicator_pt()->mpi_comm());

        // create an int array of first rows and nrow local and
        // compute nnz global
        int* dist_first_row = new int[nproc];
        int* dist_nrow_local = new int[nproc];
        for (int p = 0; p < nproc; p++)
        {
          dist_first_row[p] = this->first_row(p);
          dist_nrow_local[p] = this->nrow_local(p);
        }

        // conpute the offset for the values and column index data
        // compute the nnz offset for each processor
        int next_row = 0;
        unsigned nnz_count = 0;
        Vector<unsigned> nnz_offset(nproc, 0);
        for (int p = 0; p < nproc; p++)
        {
          unsigned pp = 0;
          while (dist_first_row[pp] != next_row)
          {
            pp++;
          }
          nnz_offset[pp] = nnz_count;
          nnz_count += dist_nnz_pt[pp];
          next_row += dist_nrow_local[pp];
        }

        // get pointers to the (current) distributed data
        int* dist_row_start = this->row_start();
        int* dist_column_index = this->column_index();
        double* dist_value = this->value();

        // space for the global matrix
        int* global_row_start = new int[nrow + 1];
        int* global_column_index = new int[nnz_count];
        double* global_value = new double[nnz_count];

        // post the sends and recvs for the matrix data
        Vector<MPI_Request> recv_req;
        Vector<MPI_Request> send_req;
        MPI_Aint base_address;
        MPI_Get_address(global_value, &base_address);

        // SEND
        if (dist_nrow_local[my_rank] > 0)
        {
          // types
          MPI_Datatype types[3];

          // offsets
          MPI_Aint offsets[3];

          // lengths
          int len[3];

          // row start
          MPI_Type_contiguous(dist_nrow_local[my_rank], MPI_INT, &types[0]);
          MPI_Type_commit(&types[0]);
          MPI_Get_address(dist_row_start, &offsets[0]);
          offsets[0] -= base_address;
          len[0] = 1;

          // value
          MPI_Type_contiguous(nnz, MPI_DOUBLE, &types[1]);
          MPI_Type_commit(&types[1]);
          MPI_Get_address(dist_value, &offsets[1]);
          offsets[1] -= base_address;
          len[1] = 1;

          // column indices
          MPI_Type_contiguous(nnz, MPI_INT, &types[2]);
          MPI_Type_commit(&types[2]);
          MPI_Get_address(dist_column_index, &offsets[2]);
          offsets[2] -= base_address;
          len[2] = 1;

          // build the send type
          MPI_Datatype send_type;
          MPI_Type_create_struct(3, len, offsets, types, &send_type);
          MPI_Type_commit(&send_type);
          MPI_Type_free(&types[0]);
          MPI_Type_free(&types[1]);
          MPI_Type_free(&types[2]);

          // and send
          for (int p = 0; p < nproc; p++)
          {
            if (p != my_rank)
            {
              MPI_Request req;
              MPI_Isend(global_value,
                        1,
                        send_type,
                        p,
                        1,
                        this->distribution_pt()->communicator_pt()->mpi_comm(),
                        &req);
              send_req.push_back(req);
            }
          }
          MPI_Type_free(&send_type);
        }

        // RECV
        for (int p = 0; p < nproc; p++)
        {
          // communicated with other processors
          if (p != my_rank)
          {
            // RECV
            if (dist_nrow_local[p] > 0)
            {
              // types
              MPI_Datatype types[3];

              // offsets
              MPI_Aint offsets[3];

              // lengths
              int len[3];

              // row start
              MPI_Type_contiguous(dist_nrow_local[p], MPI_INT, &types[0]);
              MPI_Type_commit(&types[0]);
              MPI_Get_address(global_row_start + dist_first_row[p],
                              &offsets[0]);
              offsets[0] -= base_address;
              len[0] = 1;

              // value
              MPI_Type_contiguous(dist_nnz_pt[p], MPI_DOUBLE, &types[1]);
              MPI_Type_commit(&types[1]);
              MPI_Get_address(global_value + nnz_offset[p], &offsets[1]);
              offsets[1] -= base_address;
              len[1] = 1;

              // column indices
              MPI_Type_contiguous(dist_nnz_pt[p], MPI_INT, &types[2]);
              MPI_Type_commit(&types[2]);
              MPI_Get_address(global_column_index + nnz_offset[p], &offsets[2]);
              offsets[2] -= base_address;
              len[2] = 1;

              // build the send type
              MPI_Datatype recv_type;
              MPI_Type_create_struct(3, len, offsets, types, &recv_type);
              MPI_Type_commit(&recv_type);
              MPI_Type_free(&types[0]);
              MPI_Type_free(&types[1]);
              MPI_Type_free(&types[2]);

              // and send
              MPI_Request req;
              MPI_Irecv(global_value,
                        1,
                        recv_type,
                        p,
                        1,
                        this->distribution_pt()->communicator_pt()->mpi_comm(),
                        &req);
              recv_req.push_back(req);
              MPI_Type_free(&recv_type);
            }
          }
          // otherwise send to self
          else
          {
            unsigned nrow_local = dist_nrow_local[my_rank];
            unsigned first_row = dist_first_row[my_rank];
            for (unsigned i = 0; i < nrow_local; i++)
            {
              global_row_start[first_row + i] = dist_row_start[i];
            }
            unsigned offset = nnz_offset[my_rank];
            for (int i = 0; i < nnz; i++)
            {
              global_value[offset + i] = dist_value[i];
              global_column_index[offset + i] = dist_column_index[i];
            }
          }
        }

        // wait for all recvs to complete
        unsigned n_recv_req = recv_req.size();
        if (n_recv_req > 0)
        {
          Vector<MPI_Status> recv_status(n_recv_req);
          MPI_Waitall(n_recv_req, &recv_req[0], &recv_status[0]);
        }

        // finally the last row start
        global_row_start[nrow] = nnz_count;

        // update the other row start
        for (int p = 0; p < nproc; p++)
        {
          for (int i = 0; i < dist_nrow_local[p]; i++)
          {
            unsigned j = dist_first_row[p] + i;
            global_row_start[j] += nnz_offset[p];
          }
        }

        // wait for sends to complete
        unsigned n_send_req = send_req.size();
        if (n_recv_req > 0)
        {
          Vector<MPI_Status> send_status(n_recv_req);
          MPI_Waitall(n_send_req, &send_req[0], &send_status[0]);
        }

        // rebuild the matrix
        LinearAlgebraDistribution* dist_pt = new LinearAlgebraDistribution(
          this->distribution_pt()->communicator_pt(), nrow, false);
        this->build(dist_pt);
        this->build_without_copy(
          ncol, nnz_count, global_value, global_column_index, global_row_start);

        // clean up
        delete dist_pt;
        delete[] dist_first_row;
        delete[] dist_nrow_local;
        delete[] dist_nnz_pt;
      }

      // other the matrix is not distributed but it needs to be turned
      // into a distributed matrix
      // =============================================================
      else if (!this->distributed() && dist_pt->distributed())
      {
        // cache the new nrow_local
        unsigned nrow_local = dist_pt->nrow_local();

        // and first_row
        unsigned first_row = dist_pt->first_row();

        // get pointers to the (current) distributed data
        int* global_row_start = this->row_start();
        int* global_column_index = this->column_index();
        double* global_value = this->value();

        // determine the number of non zeros required by this processor
        unsigned nnz = global_row_start[first_row + nrow_local] -
                       global_row_start[first_row];

        // allocate
        int* dist_row_start = new int[nrow_local + 1];
        int* dist_column_index = new int[nnz];
        double* dist_value = new double[nnz];

        // copy
        int offset = global_row_start[first_row];
        for (unsigned i = 0; i <= nrow_local; i++)
        {
          dist_row_start[i] = global_row_start[first_row + i] - offset;
        }
        for (unsigned i = 0; i < nnz; i++)
        {
          dist_column_index[i] = global_column_index[offset + i];
          dist_value[i] = global_value[offset + i];
        }

        // rebuild
        this->build(dist_pt);
        this->build_without_copy(
          ncol, nnz, dist_value, dist_column_index, dist_row_start);
      }
    }
#endif
  }

  //=============================================================================
  /// Compute transpose of matrix
  //=============================================================================
  void CRDoubleMatrix::get_matrix_transpose(CRDoubleMatrix* result) const
  {
    // Get the number of non_zeros
    unsigned long nnon_zeros = this->nnz();

    // Find the number of rows and columns in the transposed
    // matrix. We differentiate these from those associated
    // with the original matrix by appending the characters
    // '_t' onto the end
    unsigned long n_rows = this->nrow();
    unsigned long n_rows_t = this->ncol();

#ifdef OOMPH_HAS_MPI
    // We only need to produce a warning if the matrix is distributed
    if (this->distributed())
    {
      // Create an ostringstream object to store the warning message
      std::ostringstream warning_message;

      // Create the warning messsage
      warning_message << "This method currently works for serial but "
                      << "has not been tested with MPI!\n";

      // Issue the warning
      OomphLibWarning(warning_message.str(),
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Set up the distribution for the transposed matrix
    result->distribution_pt()->build(
      this->distribution_pt()->communicator_pt(), n_rows_t, false);

    // Acquire access to the value, row_start and column_index
    // arrays from the CR matrix
    const double* value_pt = this->value();
    const int* row_start_pt = this->row_start();
    const int* column_index_pt = this->column_index();

    // Allocate space for the row_start and column_index vectors
    // associated with the transpose of the current matrix.
    Vector<double> value_t(nnon_zeros, 0.0);
    Vector<int> column_index_t(nnon_zeros, 0);
    Vector<int> row_start_t(n_rows_t + 1, 0);

    // Loop over the column index vector and count how many times
    // each column number occurs and increment the appropriate
    // entry in row_start_t whose i+1'th entry will contain the
    // number of non-zeros in the i'th column of the matrix (this
    // is done so that the cumulative sum done next returns the
    // correct result)
    for (unsigned i = 0; i < nnon_zeros; i++)
    {
      // Assign entries to row_start_t (noting the first
      // entry is left as 0 for the cumulative sum done later)
      row_start_t[*(column_index_pt + i) + 1]++;
    }

    // Calculate the sum of the first i entries in the row_start_t
    // vector and store the values in the i'th entry of row_start_t
    for (unsigned i = 1; i < n_rows_t + 1; i++)
    {
      // Calculate the cumulative sum
      row_start_t[i] += row_start_t[i - 1];
    }

    // Allocate space for variables to store the indices of the
    // start and end of the non-zeros in a given row of the matrix
    unsigned i_row_s = 0;
    unsigned i_row_e = 0;

    // Initialise 3 extra variables for readability of the
    // code in the subsequent piece of code
    unsigned a = 0;
    unsigned b = 0;
    unsigned c = 0;

    // Vector needed to count the number of entries added
    // to each segment in column_index_t where each segment
    // is associated with each row in the transposed matrix
    Vector<int> counter(n_rows_t, 0);

    // Set the entries in column_index_t. To do this we loop
    // over each row of the original matrix (equivalently
    // the number of columns in the transpose)
    for (unsigned i_row = 0; i_row < n_rows; i_row++)
    {
      // Here we find the indices of the start and end
      // of the non-zeros in i_row'th row of the matrix.
      // [Note, there should be a -1 on i_row_e but this
      // is ignored so that we use a strict inequality
      // in the subsequent for-loop!]
      i_row_s = *(row_start_pt + i_row);
      i_row_e = *(row_start_pt + i_row + 1);

      // Loop over the entries in the i_row'th row
      // of the matrix
      for (unsigned j = i_row_s; j < i_row_e; j++)
      {
        // The value of a is the column index of the j'th
        // element in the i_row'th row of the original matrix
        // (which is also the row index of the associated
        // non-zero in the transposed matrix)
        a = *(column_index_pt + j);

        // The value of b will be used to find the start
        // of the appropriate segment of column_index_t
        // that the non-zero belongs in
        b = row_start_t[a];

        // Find the number of elements already added to
        // this segment (to find which entry of the segment
        // the value i_row, the column index of the non-zero
        // in the transposed matrix, should be assigned to)
        c = counter[*(column_index_pt + j)];

        // Assign the value i_row to the appropriate entry
        // in column_index_t
        column_index_t[b + c] = i_row;
        value_t[b + c] = *(value_pt + j);

        // Increment the j'th value of the vector counter
        // to indicate another non-zero index has been
        // added into the
        counter[*(column_index_pt + j)]++;

      } // End of for-loop over i_row'th row of the matrix

    } // End of for-loop over rows of the matrix

    // Build the matrix (note: the value of n_cols for the
    // transposed matrix is n_rows for the original matrix)
    result->build(n_rows, value_t, column_index_t, row_start_t);

  } // End of the function


  //=============================================================================
  /// Compute infinity (maximum) norm of matrix
  //=============================================================================
  double CRDoubleMatrix::inf_norm() const
  {
#ifdef PARANOID
    // paranoid check that the vector is setup
    if (!this->distribution_built())
    {
      std::ostringstream error_message;
      error_message << "This matrix must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // compute the local norm
    unsigned nrow_local = this->nrow_local();
    double n = 0;
    const int* row_start = CR_matrix.row_start();
    const double* value = CR_matrix.value();
    for (unsigned i = 0; i < nrow_local; i++)
    {
      double a = 0;
      for (int j = row_start[i]; j < row_start[i + 1]; j++)
      {
        a += fabs(value[j]);
      }
      n = std::max(n, a);
    }

    // if this vector is distributed and on multiple processors then gather
#ifdef OOMPH_HAS_MPI
    double n2 = n;
    if (this->distributed() &&
        this->distribution_pt()->communicator_pt()->nproc() > 1)
    {
      MPI_Allreduce(&n,
                    &n2,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    this->distribution_pt()->communicator_pt()->mpi_comm());
    }
    n = n2;
#endif

    // and return
    return n;
  }

  //=============================================================================
  /// Return the diagonal entries of the matrix.
  /// This only works with square matrices. This condition may be relaxed
  /// in the future if need be.
  //=============================================================================
  Vector<double> CRDoubleMatrix::diagonal_entries() const
  {
#ifdef PARANOID
    // Check if the matrix has been built.
    if (!this->built())
    {
      std::ostringstream error_message;
      error_message << "The matrix has not been built.\n"
                    << "Please build it...\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check if this is a square matrix.
    if (this->nrow() != this->ncol())
    {
      std::ostringstream error_message;
      error_message << "The matrix is not square. Can only get the diagonal\n"
                    << "entries of a square matrix.\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // We get the diagonal entries on this processor.
    unsigned nrow_local = this->nrow_local();

    // Create storage for the diagonal entries.
    Vector<double> result_vec;
    result_vec.reserve(nrow_local);

    // Get the first row for the column offset.
    unsigned first_row = this->first_row();

    // Loop through the local rows.
    for (unsigned i = 0; i < nrow_local; i++)
    {
      // The column entries are globally indexed.
      unsigned diag_entry_col = first_row + i;

      // Push back the diagonal entry.
      result_vec.push_back(CR_matrix.get_entry(i, diag_entry_col));
    }

    return result_vec;
  }

  //=============================================================================
  /// Element-wise addition of this matrix with matrix_in.
  //=============================================================================
  void CRDoubleMatrix::add(const CRDoubleMatrix& matrix_in,
                           CRDoubleMatrix& result_matrix) const
  {
#ifdef PARANOID
    // Check if this matrix is built.
    if (!this->built())
    {
      std::ostringstream error_message;
      error_message << "The matrix is not built.\n"
                    << "Please build the matrix!\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check if this matrix_in is built.
    if (!matrix_in.built())
    {
      std::ostringstream error_message;
      error_message << "The matrix matrix_in is not built.\n"
                    << "Please build the matrix!\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check if the dimensions of this matrix and matrix_in are the same.
    unsigned long this_nrow = this->nrow();
    unsigned long matrix_in_nrow = matrix_in.nrow();
    if (this_nrow != matrix_in_nrow)
    {
      std::ostringstream error_message;
      error_message << "matrix_in has a different number of rows than\n"
                    << "this matrix.\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    unsigned long this_ncol = this->ncol();
    unsigned long matrix_in_ncol = matrix_in.ncol();
    if (this_ncol != matrix_in_ncol)
    {
      std::ostringstream error_message;
      error_message << "matrix_in has a different number of columns than\n"
                    << "this matrix.\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check if the distribution is the same (Otherwise we may have to send and
    // receive data from other processors - which is not implemented!)
    if (*(this->distribution_pt()) != *(matrix_in.distribution_pt()))
    {
      std::ostringstream error_message;
      error_message << "matrix_in must have the same distribution as\n"
                    << "this matrix.\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // If the matrix is built, check that it's existing distribution is the
    // same as the in matrix. Since we'll use the same distribution instead
    // of completely rebuilding it.
    if (result_matrix.built() &&
        (*result_matrix.distribution_pt() != *matrix_in.distribution_pt()))
    {
      std::ostringstream error_message;
      error_message << "The result_matrix is built. "
                    << "But has a different distribution from matrix_in \n"
                    << "They need to be the same.\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // To add the elements of two CRDoubleMatrices, we need to know the union of
    // the sparsity patterns. This is determined by the column indices.
    // We add the column indices and values (entries) as a key-value pair in
    // a map (per row). We then read these out into two column indices and
    // values vector for the result matrix.

    unsigned nrow_local = this->nrow_local();
    Vector<int> res_column_indices;
    Vector<double> res_values;
    Vector<int> res_row_start;
    res_row_start.reserve(nrow_local + 1);

    // The row_start and column_indices
    const int* this_column_indices = this->column_index();
    const int* this_row_start = this->row_start();
    const int* in_column_indices = matrix_in.column_index();
    const int* in_row_start = matrix_in.row_start();

    // Values from this matrix and matrix_in.
    const double* this_values = this->value();
    const double* in_values = matrix_in.value();


    // The first entry in row_start is always zero.
    res_row_start.push_back(0);

    // Loop through the rows of both matrices and insert the column indices and
    // values as a key-value pair.
    for (unsigned row_i = 0; row_i < nrow_local; row_i++)
    {
      // Create the map for this row.
      std::map<int, double> res_row_map;

      // Insert the column and value pair for this matrix.
      for (int i = this_row_start[row_i]; i < this_row_start[row_i + 1]; i++)
      {
        res_row_map[this_column_indices[i]] = this_values[i];
      }

      // Insert the column and value pair for in matrix.
      for (int i = in_row_start[row_i]; i < in_row_start[row_i + 1]; i++)
      {
        res_row_map[in_column_indices[i]] += in_values[i];
      }

      // Fill in the row start
      res_row_start.push_back(res_row_start.back() + res_row_map.size());

      // Now insert the key into res_column_indices and value into res_values
      for (std::map<int, double>::iterator it = res_row_map.begin();
           it != res_row_map.end();
           ++it)
      {
        res_column_indices.push_back(it->first);
        res_values.push_back(it->second);
      }
    }

    // Finally build the result_matrix.
    if (result_matrix.distribution_pt()->built())
    {
      // Use the existing distribution.
      result_matrix.build(
        this->ncol(), res_values, res_column_indices, res_row_start);
    }
    else
    {
      // Build with THIS distribution
      result_matrix.build(this->distribution_pt(),
                          this->ncol(),
                          res_values,
                          res_column_indices,
                          res_row_start);
    }
  }

  //=================================================================
  /// Namespace for helper functions for CRDoubleMatrices
  //=================================================================
  namespace CRDoubleMatrixHelpers
  {
    //============================================================================
    ///  Builds a uniformly distributed matrix.
    /// A locally replicated matrix is constructed then redistributed using
    /// OOMPH-LIB's default uniform row distribution.
    /// This is memory intensive thus should be used for
    /// testing or small problems only.
    //============================================================================
    void create_uniformly_distributed_matrix(
      const unsigned& nrow,
      const unsigned& ncol,
      const OomphCommunicator* const comm_pt,
      const Vector<double>& values,
      const Vector<int>& column_indices,
      const Vector<int>& row_start,
      CRDoubleMatrix& matrix_out)
    {
#ifdef PARANOID
      // Check if the communicator exists.
      if (comm_pt == 0)
      {
        std::ostringstream error_message;
        error_message << "Please supply the communicator.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      // Is the out matrix built? We need an empty matrix!
      if (matrix_out.built())
      {
        std::ostringstream error_message;
        error_message << "The result matrix has been built.\n"
                      << "Please clear the matrix.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Create the locally replicated distribution.
      bool distributed = false;
      LinearAlgebraDistribution locally_replicated_distribution(
        comm_pt, nrow, distributed);

      // Create the matrix.
      matrix_out.build(&locally_replicated_distribution,
                       ncol,
                       values,
                       column_indices,
                       row_start);

      // Create the distributed distribution.
      distributed = true;
      LinearAlgebraDistribution distributed_distribution(
        comm_pt, nrow, distributed);

      // Redistribute the matrix.
      matrix_out.redistribute(&distributed_distribution);
    }

    //============================================================================
    /// Compute infinity (maximum) norm of sub blocks as if it was one matrix
    //============================================================================
    double inf_norm(const DenseMatrix<CRDoubleMatrix*>& matrix_pt)
    {
      // The number of block rows and columns
      const unsigned nblockrow = matrix_pt.nrow();
      const unsigned nblockcol = matrix_pt.ncol();

#ifdef PARANOID
      // Check that tehre is at least one matrix.
      if (matrix_pt.nrow() == 0)
      {
        std::ostringstream error_message;
        error_message << "There are no matrices... \n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }


      // Check that all matrix_pt pointers are not null
      // and the matrices are built.
      for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
      {
        for (unsigned block_col_i = 0; block_col_i < nblockcol; block_col_i++)
        {
          if (matrix_pt(block_row_i, block_col_i) == 0)
          {
            std::ostringstream error_message;
            error_message << "The pointer martrix_pt(" << block_row_i << ","
                          << block_col_i << ") is null.\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          if (!matrix_pt(block_row_i, block_col_i)->built())
          {
            std::ostringstream error_message;
            error_message << "The matrix at martrix_pt(" << block_row_i << ","
                          << block_col_i << ") is not built.\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
#endif

#ifdef OOMPH_HAS_MPI

      // The communicator pointer from block (0,0)
      const OomphCommunicator* const comm_pt =
        matrix_pt(0, 0)->distribution_pt()->communicator_pt();

#ifdef PARANOID


      // Check that all communicators are the same
      for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
      {
        for (unsigned block_col_i = 0; block_col_i < nblockcol; block_col_i++)
        {
          // Communicator for this block matrix.
          const OomphCommunicator current_block_comm =
            *(matrix_pt(block_row_i, block_col_i)
                ->distribution_pt()
                ->communicator_pt());
          if (*comm_pt != current_block_comm)
          {
            std::ostringstream error_message;
            error_message << "The communicator of block martrix_pt("
                          << block_row_i << "," << block_col_i
                          << ") is not the same as block "
                          << "matrix_pt(0,0).\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // Check that all distributed boolean are the same (if on more than 1
      // core)
      if (comm_pt->nproc() > 1)
      {
        // Get the distributed boolean from matrix_pt(0,0)
        bool first_distributed = matrix_pt(0, 0)->distributed();

        for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
        {
          for (unsigned block_col_i = 0; block_col_i < nblockcol; block_col_i++)
          {
            // Is the current block distributed?
            bool current_distributed =
              matrix_pt(block_row_i, block_col_i)->distributed();

            if (first_distributed != current_distributed)
            {
              std::ostringstream error_message;
              error_message << "Block matrix_pt(" << block_row_i << ","
                            << block_col_i << ") and block matrix_pt(0,0) "
                            << "have a different distributed boolean.\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }

      // Check that all sub matrix dimensions "make sense"
      // We need to check that all the matrices in the same row has the same
      // nrow. Then repeat for the columns.

      // Check the nrow of each block row.
      for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
      {
        // Get the nrow to compare against from the first column.
        const unsigned first_block_nrow = matrix_pt(block_row_i, 0)->nrow();

        // Loop through the block columns.
        for (unsigned block_col_i = 1; block_col_i < nblockcol; block_col_i++)
        {
          // If the nrow of this block is not the same as the nrow from the
          // first block in this block row, throw an error.
          const unsigned current_block_nrow =
            matrix_pt(block_row_i, block_col_i)->nrow();

          if (first_block_nrow != current_block_nrow)
          {
            std::ostringstream error_message;
            error_message << "First block has nrow = " << current_block_nrow
                          << ". But martrix_pt(" << block_row_i << ","
                          << block_col_i
                          << ") has nrow = " << current_block_nrow << ".\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // Check the ncol of each block column.
      for (unsigned block_col_i = 0; block_col_i < nblockcol; block_col_i++)
      {
        // Get the ncol from the first block row to compare against.
        const unsigned first_block_ncol = matrix_pt(0, block_col_i)->ncol();

        for (unsigned block_row_i = 1; block_row_i < nblockrow; block_row_i++)
        {
          // Get the ncol for the current block.
          const unsigned current_block_ncol =
            matrix_pt(block_row_i, block_col_i)->ncol();

          if (first_block_ncol != current_block_ncol)
          {
            std::ostringstream error_message;
            error_message << "First block has ncol = " << current_block_ncol
                          << ". But martrix_pt(" << block_row_i << ","
                          << block_col_i
                          << ") has ncol = " << current_block_ncol << ".\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // Check that the distribution for each block row is the same.
      for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
      {
        // The first distribution of this block row.
        const LinearAlgebraDistribution first_dist =
          *(matrix_pt(block_row_i, 0)->distribution_pt());

        // Loop through the rest of the block columns.
        for (unsigned block_col_i = 1; block_col_i < nblockcol; block_col_i++)
        {
          // Get the distribution from the current block.
          const LinearAlgebraDistribution current_dist =
            matrix_pt(block_row_i, block_col_i)->distribution_pt();

          // Compare the first distribution against the current.
          if (first_dist != current_dist)
          {
            std::ostringstream error_message;
            error_message << "First distribution of block row " << block_row_i
                          << " is different from the distribution from "
                          << "martrix_pt(" << block_row_i << "," << block_col_i
                          << ").\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
#endif

#endif

      // Loop thrpugh the block rows, then block columns to
      // compute the local inf norm
      double inf_norm = 0;
      for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
      {
        // Get the number of local rows from the first block.
        unsigned block_nrow_local = matrix_pt(block_row_i, 0)->nrow_local();

        // Loop through the block_nrow_local in this block row
        for (unsigned local_row_i = 0; local_row_i < block_nrow_local;
             local_row_i++)
        {
          double abs_sum_of_row = 0;
          // Loop through the block columns
          for (unsigned block_col_i = 0; block_col_i < nblockcol; block_col_i++)
          {
            // Locally cache the pointer to the current block.
            CRDoubleMatrix* block_pt = matrix_pt(block_row_i, block_col_i);

            const int* row_start = block_pt->row_start();
            const double* value = block_pt->value();

            // Loop through the values
            for (int val_i = row_start[local_row_i];
                 val_i < row_start[local_row_i + 1];
                 val_i++)
            {
              abs_sum_of_row += fabs(value[val_i]);
            }
          }
          // Store the max row
          inf_norm = std::max(inf_norm, abs_sum_of_row);
        }
      }

      // if this vector is distributed and on multiple processors then gather
#ifdef OOMPH_HAS_MPI
      double inf_norm_local = inf_norm;
      if (matrix_pt(0, 0)->distributed() && comm_pt->nproc() > 1)
      {
        MPI_Allreduce(&inf_norm,
                      &inf_norm_local,
                      1,
                      MPI_DOUBLE,
                      MPI_MAX,
                      comm_pt->mpi_comm());
      }
      inf_norm = inf_norm_local;
#endif

      // and return
      return inf_norm;
    }

    //============================================================================
    ///  Calculates the largest Gershgorin disc whilst preserving the
    /// sign. Let A be an n by n matrix, with entries aij. For \f$ i \in \{
    /// 1,...,n \} \f$ let \f$ R_i = \sum_{i\neq j} |a_{ij}| \f$ be the sum of
    /// the absolute values of the non-diagonal entries in the i-th row. Let \f$
    /// D(a_{ii},R_i) \f$ be the closed disc centered at \f$ a_{ii} \f$ with
    /// radius \f$ R_i \f$, such a disc is called a Gershgorin disc.
    ///
    /// \n
    ///
    /// We calculate \f$ |D(a_{ii},R_i)|_{max} \f$ and multiply by the sign of
    /// the diagonal entry.
    ///
    /// \n
    ///
    /// The DenseMatrix of CRDoubleMatrices are treated as if they are one
    /// large matrix. Therefore the dimensions of the sub matrices has to
    /// "make sense", there is a paranoid check for this.
    //============================================================================
    double gershgorin_eigenvalue_estimate(
      const DenseMatrix<CRDoubleMatrix*>& matrix_pt)
    {
      // The number of block rows and columns
      const unsigned nblockrow = matrix_pt.nrow();
      const unsigned nblockcol = matrix_pt.ncol();

#ifdef PARANOID
      // Check that tehre is at least one matrix.
      if (matrix_pt.nrow() == 0)
      {
        std::ostringstream error_message;
        error_message << "There are no matrices... \n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }


      // Check that all matrix_pt pointers are not null
      // and the matrices are built.
      for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
      {
        for (unsigned block_col_i = 0; block_col_i < nblockcol; block_col_i++)
        {
          if (matrix_pt(block_row_i, block_col_i) == 0)
          {
            std::ostringstream error_message;
            error_message << "The pointer martrix_pt(" << block_row_i << ","
                          << block_col_i << ") is null.\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          if (!matrix_pt(block_row_i, block_col_i)->built())
          {
            std::ostringstream error_message;
            error_message << "The matrix at martrix_pt(" << block_row_i << ","
                          << block_col_i << ") is not built.\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
#endif


#ifdef OOMPH_HAS_MPI

      // The communicator pointer from block (0,0)
      // All communicators should be the same, we check this next.
      const OomphCommunicator* const comm_pt =
        matrix_pt(0, 0)->distribution_pt()->communicator_pt();

#ifdef PARANOID

      // Check that all communicators are the same
      for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
      {
        for (unsigned block_col_i = 0; block_col_i < nblockcol; block_col_i++)
        {
          // Communicator for this block matrix.
          const OomphCommunicator current_block_comm =
            *(matrix_pt(block_row_i, block_col_i)
                ->distribution_pt()
                ->communicator_pt());
          if (*comm_pt != current_block_comm)
          {
            std::ostringstream error_message;
            error_message << "The communicator of block martrix_pt("
                          << block_row_i << "," << block_col_i
                          << ") is not the same as block "
                          << "matrix_pt(0,0).\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // Check that all distributed boolean are the same (if on more than 1
      // core)
      if (comm_pt->nproc() > 1)
      {
        // Get the distributed boolean from matrix_pt(0,0)
        bool first_distributed = matrix_pt(0, 0)->distributed();

        for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
        {
          for (unsigned block_col_i = 0; block_col_i < nblockcol; block_col_i++)
          {
            // Is the current block distributed?
            bool current_distributed =
              matrix_pt(block_row_i, block_col_i)->distributed();

            if (first_distributed != current_distributed)
            {
              std::ostringstream error_message;
              error_message << "Block matrix_pt(" << block_row_i << ","
                            << block_col_i << ") and block matrix_pt(0,0) "
                            << "have a different distributed boolean.\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }

      // Check that all sub matrix dimensions "make sense"
      // We need to check that all the matrices in the same row has the same
      // nrow. Then repeat for the columns.

      // Check the nrow of each block row.
      for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
      {
        // Get the nrow to compare against from the first column.
        const unsigned first_block_nrow = matrix_pt(block_row_i, 0)->nrow();

        // Loop through the block columns.
        for (unsigned block_col_i = 1; block_col_i < nblockcol; block_col_i++)
        {
          // If the nrow of this block is not the same as the nrow from the
          // first block in this block row, throw an error.
          const unsigned current_block_nrow =
            matrix_pt(block_row_i, block_col_i)->nrow();

          if (first_block_nrow != current_block_nrow)
          {
            std::ostringstream error_message;
            error_message << "First block has nrow = " << current_block_nrow
                          << ". But martrix_pt(" << block_row_i << ","
                          << block_col_i
                          << ") has nrow = " << current_block_nrow << ".\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // Check the ncol of each block column.
      for (unsigned block_col_i = 0; block_col_i < nblockcol; block_col_i++)
      {
        // Get the ncol from the first block row to compare against.
        const unsigned first_block_ncol = matrix_pt(0, block_col_i)->ncol();

        for (unsigned block_row_i = 1; block_row_i < nblockrow; block_row_i++)
        {
          // Get the ncol for the current block.
          const unsigned current_block_ncol =
            matrix_pt(block_row_i, block_col_i)->ncol();

          if (first_block_ncol != current_block_ncol)
          {
            std::ostringstream error_message;
            error_message << "First block has ncol = " << current_block_ncol
                          << ". But martrix_pt(" << block_row_i << ","
                          << block_col_i
                          << ") has ncol = " << current_block_ncol << ".\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // Check that the distribution for each block row is the same.
      for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
      {
        // The first distribution of this block row.
        const LinearAlgebraDistribution first_dist =
          *(matrix_pt(block_row_i, 0)->distribution_pt());

        // Loop through the rest of the block columns.
        for (unsigned block_col_i = 1; block_col_i < nblockcol; block_col_i++)
        {
          // Get the distribution from the current block.
          const LinearAlgebraDistribution current_dist =
            matrix_pt(block_row_i, block_col_i)->distribution_pt();

          // Compare the first distribution against the current.
          if (first_dist != current_dist)
          {
            std::ostringstream error_message;
            error_message << "First distribution of block row " << block_row_i
                          << " is different from the distribution from "
                          << "martrix_pt(" << block_row_i << "," << block_col_i
                          << ").\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

#endif
#endif

      // Loop thrpugh the block rows, then block columns to
      // compute the local inf norm
      double extreme_disc = 0;
      for (unsigned block_row_i = 0; block_row_i < nblockrow; block_row_i++)
      {
        // Get the number of local rows from the first block.
        unsigned block_nrow_local = matrix_pt(block_row_i, 0)->nrow_local();

        // Loop through the block_nrow_local in this block row
        for (unsigned local_row_i = 0; local_row_i < block_nrow_local;
             local_row_i++)
        {
          double abs_sum_of_row = 0;
          // Loop through the block columns
          for (unsigned block_col_i = 0; block_col_i < nblockcol; block_col_i++)
          {
            // Locally cache the pointer to the current block.
            CRDoubleMatrix* block_pt = matrix_pt(block_row_i, block_col_i);

            const int* row_start = block_pt->row_start();
            const double* value = block_pt->value();

            // Loop through the values
            for (int val_i = row_start[local_row_i];
                 val_i < row_start[local_row_i + 1];
                 val_i++)
            {
              abs_sum_of_row += fabs(value[val_i]);
            }
          }

          // Now minus the diagonal entry...
          // Locate the diagonal block matrix.
          double* s_values = matrix_pt(block_row_i, block_row_i)->value();
          int* s_column_index =
            matrix_pt(block_row_i, block_row_i)->column_index();
          int* s_row_start = matrix_pt(block_row_i, block_row_i)->row_start();
          // int s_nrow_local =
          // matrix_pt(block_row_i,block_row_i)->nrow_local();
          int s_first_row = matrix_pt(block_row_i, block_row_i)->first_row();

          // Get the diagonal value...
          double diagonal_value = 0;
          bool found = false;
          for (int j = s_row_start[local_row_i];
               j < s_row_start[local_row_i + 1] && !found;
               j++)
          {
            if (s_column_index[j] == int(local_row_i + s_first_row))
            {
              diagonal_value = s_values[j];
              found = true;
            }
          }

          // Check if the diagonal entry is found.
          if (!found)
          {
            std::ostringstream error_message;
            error_message << "The diagonal entry for the block(" << block_row_i
                          << "," << block_row_i << ")\n"
                          << "on local row " << local_row_i
                          << " does not exist." << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          // This is the disc.
          abs_sum_of_row -= fabs(diagonal_value);

          // Now we have to check if the diagonal entry is
          // on the left or right side of zero.
          if (diagonal_value > 0)
          {
            double extreme_disc_local = diagonal_value + abs_sum_of_row;
            extreme_disc = std::max(extreme_disc_local, extreme_disc);
          }
          else
          {
            double extreme_disc_local = diagonal_value - abs_sum_of_row;
            extreme_disc = std::min(extreme_disc_local, extreme_disc);
          }
        } // Loop through local row (of all block column)
      } // Loop through block row

      // if this vector is distributed and on multiple processors then gather
#ifdef OOMPH_HAS_MPI
      double extreme_disc_local = extreme_disc;
      if (matrix_pt(0, 0)->distributed() && comm_pt->nproc() > 1)
      {
        if (extreme_disc > 0)
        {
          MPI_Allreduce(&extreme_disc,
                        &extreme_disc_local,
                        1,
                        MPI_DOUBLE,
                        MPI_MAX,
                        comm_pt->mpi_comm());
        }
        else
        {
          MPI_Allreduce(&extreme_disc,
                        &extreme_disc_local,
                        1,
                        MPI_DOUBLE,
                        MPI_MIN,
                        comm_pt->mpi_comm());
        }
      }
      extreme_disc = extreme_disc_local;
#endif

      // and return
      return extreme_disc;
    }

    //============================================================================
    ///  Concatenate CRDoubleMatrix matrices.
    /// The in matrices are concatenated such that the block structure of the
    /// in matrices are preserved in the result matrix. Communication between
    /// processors is required. If the block structure of the sub matrices does
    /// not need to be preserved, consider using
    /// CRDoubleMatrixHelpers::concatenate_without_communication(...).
    ///
    /// The matrix manipulation functions
    /// CRDoubleMatrixHelpers::concatenate(...) and
    /// CRDoubleMatrixHelpers::concatenate_without_communication(...)
    /// are analogous to the Vector manipulation functions
    /// DoubleVectorHelpers::concatenate(...) and
    /// DoubleVectorHelpers::concatenate_without_communication(...).
    /// Please look at the DoubleVector functions for an illustration of the
    /// differences between concatenate(...) and
    /// concatenate_without_communication(...).
    ///
    /// Distribution of the result matrix:
    /// If the result matrix does not have a distribution built, then it will be
    /// given a uniform row distribution. Otherwise we use the existing
    /// distribution. This gives the user the ability to define their own
    /// distribution, or save computing power if a distribution has
    /// been pre-built.
    ///
    /// NOTE: ALL the matrices pointed to by matrix_pt has to be built. This is
    /// not the case with concatenate_without_communication(...)
    //============================================================================
    void concatenate(const DenseMatrix<CRDoubleMatrix*>& matrix_pt,
                     CRDoubleMatrix& result_matrix)
    {
      // The number of block rows and block columns.
      unsigned matrix_nrow = matrix_pt.nrow();
      unsigned matrix_ncol = matrix_pt.ncol();

      // PARANOID checks involving only the in matrices.
#ifdef PARANOID
      // Are there matrices to concatenate?
      if (matrix_nrow == 0)
      {
        std::ostringstream error_message;
        error_message << "There are no matrices to concatenate.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Does this matrix need concatenating?
      if ((matrix_nrow == 1) && (matrix_ncol == 1))
      {
        std::ostringstream warning_message;
        warning_message << "There is only one matrix to concatenate...\n"
                        << "This does not require concatenating...\n";
        OomphLibWarning(warning_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
      }

      // Are all sub matrices built?
      for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
      {
        for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
        {
          if (!(matrix_pt(block_row_i, block_col_i)->built()))
          {
            std::ostringstream error_message;
            error_message << "The sub matrix (" << block_row_i << ","
                          << block_col_i << ")\n"
                          << "is not built. \n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // Do all dimensions of sub matrices "make sense"?
      // Compare the number of rows of each block matrix in a block row.
      for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
      {
        // Use the first column to compare against the rest.
        unsigned long current_block_nrow = matrix_pt(block_row_i, 0)->nrow();

        // Compare against columns 1 to matrix_ncol - 1
        for (unsigned block_col_i = 1; block_col_i < matrix_ncol; block_col_i++)
        {
          // Get the nrow for this sub block.
          unsigned long subblock_nrow =
            matrix_pt(block_row_i, block_col_i)->nrow();

          if (current_block_nrow != subblock_nrow)
          {
            std::ostringstream error_message;
            error_message << "The sub matrix (" << block_row_i << ","
                          << block_col_i << ")\n"
                          << "requires nrow = " << current_block_nrow
                          << ", but has nrow = " << subblock_nrow << ".\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // Compare the number of columns of each block matrix in a block column.
      for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
      {
        // Use the first row to compare against the rest.
        unsigned long current_block_ncol = matrix_pt(0, block_col_i)->ncol();

        // Compare against rows 1 to matrix_nrow - 1
        for (unsigned block_row_i = 1; block_row_i < matrix_nrow; block_row_i++)
        {
          // Get the ncol for this sub block.
          unsigned long subblock_ncol =
            matrix_pt(block_row_i, block_col_i)->ncol();

          if (current_block_ncol != subblock_ncol)
          {
            std::ostringstream error_message;
            error_message << "The sub matrix (" << block_row_i << ","
                          << block_col_i << ")\n"
                          << "requires ncol = " << current_block_ncol
                          << ", but has ncol = " << subblock_ncol << ".\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
#endif

      // The communicator pointer from block (0,0)
      const OomphCommunicator* const comm_pt =
        matrix_pt(0, 0)->distribution_pt()->communicator_pt();

      // Check if the block (0,0) is distributed or not.
      bool distributed = matrix_pt(0, 0)->distributed();

      // If the result matrix does not have a distribution, we create a uniform
      // distribution.
      if (!result_matrix.distribution_pt()->built())
      {
        // Sum of sub matrix nrow. We use the first column.
        unsigned tmp_nrow = 0;
        for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
        {
          tmp_nrow += matrix_pt(block_row_i, 0)->nrow();
        }

        LinearAlgebraDistribution tmp_distribution(
          comm_pt, tmp_nrow, distributed);

        result_matrix.build(&tmp_distribution);
      }
      else
      // A distribution is supplied for the result matrix.
      {
#ifdef PARANOID
        // Check if the sum of the nrow from the sub matrices is the same as the
        // the nrow from the result matrix.

        // Sum of sub matrix nrow. We use the first column.
        unsigned tmp_nrow = 0;
        for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
        {
          tmp_nrow += matrix_pt(block_row_i, 0)->nrow();
        }

        if (tmp_nrow != result_matrix.nrow())
        {
          std::ostringstream error_message;
          error_message << "The total number of rows from the matrices to\n"
                        << "concatenate does not match the nrow from the\n"
                        << "result matrix\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }

#ifdef PARANOID

      // Are all the communicators the same?
      // Compare the communicator for sub matrices (against the result matrix).
      {
        const OomphCommunicator communicator =
          *(result_matrix.distribution_pt()->communicator_pt());

        // Are all communicator pointers the same?
        for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
        {
          for (unsigned block_col_i = 0; block_col_i < matrix_ncol;
               block_col_i++)
          {
            const OomphCommunicator another_communicator =
              *(matrix_pt(block_row_i, block_col_i)
                  ->distribution_pt()
                  ->communicator_pt());

            if (!(communicator == another_communicator))
            {
              std::ostringstream error_message;
              error_message << "The OomphCommunicator of the sub matrix ("
                            << block_row_i << "," << block_col_i << ")\n"
                            << "does not have the same communicator as the "
                               "result matrix. \n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }

      // Are all the distributed boolean the same? This only applies if we have
      // more than one processor. If there is only one processor, then it does
      // not matter if it is distributed or not - they are conceptually the
      // same.
      if (comm_pt->nproc() != 1)
      {
        // Compare distributed for sub matrices (against the result matrix).
        const bool res_distributed = result_matrix.distributed();

        // Loop over all sub blocks.
        for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
        {
          for (unsigned block_col_i = 0; block_col_i < matrix_ncol;
               block_col_i++)
          {
            const bool another_distributed =
              matrix_pt(block_row_i, block_col_i)->distributed();

            if (res_distributed != another_distributed)
            {
              std::ostringstream error_message;
              error_message << "The distributed boolean of the sub matrix ("
                            << block_row_i << "," << block_col_i << ")\n"
                            << "is not the same as the result matrix. \n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }
#endif


      // Get the number of columns up to each block for offset
      // in calculating the result column indices.
      // Since the number of columns in each block column is the same,
      // we only loop through the first block row (row zero).
      Vector<unsigned long> sum_of_ncol_up_to_block(matrix_ncol);

      // Also compute the total number of columns to build the resulting matrix.
      unsigned long res_ncol = 0;

      for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
      {
        sum_of_ncol_up_to_block[block_col_i] = res_ncol;
        res_ncol += matrix_pt(0, block_col_i)->ncol();
      }

      // We begin the process of extracting and ordering the values,
      // column_indices and row_start of all the sub blocks.
      if ((comm_pt->nproc() == 1) || !distributed)
      // Serial version of the code.
      {
        // Get the total number of non zero entries so we can reserve storage
        // for the values and column_indices vectors.
        unsigned long res_nnz = 0;
        for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
        {
          for (unsigned block_col_i = 0; block_col_i < matrix_ncol;
               block_col_i++)
          {
            res_nnz += matrix_pt(block_row_i, block_col_i)->nnz();
          }
        }

        // Declare the vectors required to build a CRDoubleMatrix
        Vector<double> res_values;
        Vector<int> res_column_indices;
        Vector<int> res_row_start;

        // Reserve space for the vectors.
        res_values.reserve(res_nnz);
        res_column_indices.reserve(res_nnz);
        res_row_start.reserve(result_matrix.nrow() + 1);

        // Now we fill in the data.

        // Running sum of nnz per row.
        int nnz_running_sum = 0;

        // Loop through the block rows.
        for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
        {
          // Get the number of rows in this block row, from the first block.
          unsigned long block_row_nrow = matrix_pt(block_row_i, 0)->nrow();

          // Loop through the number of rows in this block row
          for (unsigned row_i = 0; row_i < block_row_nrow; row_i++)
          {
            // The row start is the nnz at the start of the row.
            res_row_start.push_back(nnz_running_sum);

            // Loop through the block columns
            for (unsigned block_col_i = 0; block_col_i < matrix_ncol;
                 block_col_i++)
            {
              // Get the current block.
              CRDoubleMatrix* current_block_pt =
                matrix_pt(block_row_i, block_col_i);

              // Get the values, column_indices and row_start for this block.
              double* current_block_values = current_block_pt->value();
              int* current_block_column_indices =
                current_block_pt->column_index();
              int* current_block_row_start = current_block_pt->row_start();

              for (int val_i = current_block_row_start[row_i];
                   val_i < current_block_row_start[row_i + 1];
                   val_i++)
              {
                res_values.push_back(current_block_values[val_i]);
                res_column_indices.push_back(
                  current_block_column_indices[val_i] +
                  sum_of_ncol_up_to_block[block_col_i]);
              }

              // Update the running sum of nnz per row
              nnz_running_sum += current_block_row_start[row_i + 1] -
                                 current_block_row_start[row_i];
            } // for block cols
          } // for rows
        } // for block rows

        // Fill in the last row start
        res_row_start.push_back(res_nnz);

        // Build the matrix
        result_matrix.build(
          res_ncol, res_values, res_column_indices, res_row_start);
      }
      // Otherwise we are dealing with a distributed matrix.
      else
      {
#ifdef OOMPH_HAS_MPI

        // Flag to enable timing. This is for debugging
        // and/or testing purposes only.
        bool enable_timing = false;

        // Get the number of processors
        unsigned nproc = comm_pt->nproc();

        // My rank
        unsigned my_rank = comm_pt->my_rank();

        // Storage for the data (per processor) to send.
        Vector<Vector<unsigned>> column_indices_to_send(nproc);
        Vector<Vector<double>> values_to_send(nproc);

        // The sum of the nrow for the sub blocks (so far). This is used as an
        // offset to calculate the global equation number in the result matrix.
        unsigned long sum_of_block_nrow = 0;

        double t_prep_data_start;
        if (enable_timing)
        {
          t_prep_data_start = TimingHelpers::timer();
        }

        // Get the pointer to the result distribution, for convenience...
        LinearAlgebraDistribution* res_distribution_pt =
          result_matrix.distribution_pt();

        // loop over the sub blocks to calculate the global_eqn, get the values
        // and column indices.
        for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
        {
          // Get the number of local rows in this block_row from the first
          // block.
          unsigned current_block_nrow_local =
            matrix_pt(block_row_i, 0)->nrow_local();

          // Get the first_row for this block_row
          unsigned current_block_row_first_row =
            matrix_pt(block_row_i, 0)->first_row();

          // Loop through the number of local rows
          for (unsigned sub_local_eqn = 0;
               sub_local_eqn < current_block_nrow_local;
               sub_local_eqn++)
          {
            // Calculate the corresponding (res_global_eqn) equation number
            // for this local row number in this block.
            unsigned long res_global_eqn =
              sub_local_eqn + current_block_row_first_row + sum_of_block_nrow;

            // Get the processor that this global row belongs to.
            // The rank_of_global_row(...) function loops through all the
            // processors and does two unsigned comparisons. Since we have to do
            // this for every row, it may be better to store a list mapping for
            // very large number of processors.
            unsigned res_p =
              res_distribution_pt->rank_of_global_row(res_global_eqn);

            // With the res_p, we get the res_first_row to
            // work out the res_local_eqn
            unsigned res_first_row = res_distribution_pt->first_row(res_p);
            unsigned res_local_eqn = res_global_eqn - res_first_row;

            // Loop through the block columns, calculate the nnz. This is used
            // to reserve space for the value and column_indices Vectors.
            unsigned long current_row_nnz = 0;
            for (unsigned block_col_i = 0; block_col_i < matrix_ncol;
                 block_col_i++)
            {
              // Get the row_start
              int* current_block_row_start =
                matrix_pt(block_row_i, block_col_i)->row_start();

              // Update the nnz for this row.
              current_row_nnz += current_block_row_start[sub_local_eqn + 1] -
                                 current_block_row_start[sub_local_eqn];
            } // for block column, get nnz.

            // Reserve space for efficiency.
            // unsigned capacity_in_res_p_vec
            //  = column_indices_to_send[res_p].capacity();

            // Reserve memory for nnz+2, since we need to store the
            // res_local_eqn and nnz as well as the data (values/column
            // indices). Note: The two reserve functions are called per row. If
            // the matrix is very sparse (just a few elements per row), it will
            // be more efficient to not reserve and let the STL vector handle
            // this. On average, this implementation is more efficient.
            // column_indices_to_send[res_p].reserve(capacity_in_res_p_vec
            //    + current_row_nnz+2);
            // values_to_send[res_p].reserve(capacity_in_res_p_vec
            //    + current_row_nnz+2);

            // Push back the res_local_eqn and nnz
            column_indices_to_send[res_p].push_back(res_local_eqn);
            column_indices_to_send[res_p].push_back(current_row_nnz);
            values_to_send[res_p].push_back(res_local_eqn);
            values_to_send[res_p].push_back(current_row_nnz);

            // Loop through the block columns again and get the values
            // and column_indices
            for (unsigned block_col_i = 0; block_col_i < matrix_ncol;
                 block_col_i++)
            {
              // Cache the pointer to the current block for convenience.
              CRDoubleMatrix* current_block_pt =
                matrix_pt(block_row_i, block_col_i);

              // Values, column indices and row_start for the current block.
              double* current_block_values = current_block_pt->value();
              int* current_block_column_indices =
                current_block_pt->column_index();
              int* current_block_row_start = current_block_pt->row_start();

              // Loop though the values and column_indices
              for (int val_i = current_block_row_start[sub_local_eqn];
                   val_i < current_block_row_start[sub_local_eqn + 1];
                   val_i++)
              {
                // Push back the value.
                values_to_send[res_p].push_back(current_block_values[val_i]);

                // Push back the (offset) column index.
                column_indices_to_send[res_p].push_back(
                  current_block_column_indices[val_i] +
                  sum_of_ncol_up_to_block[block_col_i]);
              } // for block columns
            } // for block column, get values and column_indices.
          } // for sub_local_eqn

          // update the sum_of_block_nrow
          sum_of_block_nrow += matrix_pt(block_row_i, 0)->nrow();

        } // for block_row

        if (enable_timing)
        {
          double t_prep_data_finish = TimingHelpers::timer();
          double t_prep_data_time = t_prep_data_finish - t_prep_data_start;
          oomph_info << "Time for prep data: " << t_prep_data_time << std::endl;
        }

        // Prepare to send data!

        // Storage for the number of data to be sent to each processor.
        Vector<int> send_n(nproc, 0);

        // Storage for all the values/column indices to be sent
        // to each processor.
        Vector<double> send_values_data;
        Vector<unsigned> send_column_indices_data;

        // Storage location within send_values_data
        // (and send_column_indices_data) for data to be sent to each processor.
        Vector<int> send_displacement(nproc, 0);

        double t_total_ndata_start;
        if (enable_timing) t_total_ndata_start = TimingHelpers::timer();

        // Get the total amount of data which needs to be sent, so we can
        // reserve space for it.
        unsigned total_ndata = 0;
        for (unsigned rank = 0; rank < nproc; rank++)
        {
          if (rank != my_rank)
          {
            total_ndata += values_to_send[rank].size();
          }
        }

        if (enable_timing)
        {
          double t_total_ndata_finish = TimingHelpers::timer();
          double t_total_ndata_time =
            t_total_ndata_finish - t_total_ndata_start;
          oomph_info << "Time for total_ndata: " << t_total_ndata_time
                     << std::endl;
        }

        double t_flat_pack_start;
        if (enable_timing) t_flat_pack_start = TimingHelpers::timer();

        // Now we don't have to re-allocate data/memory when push_back is
        // called. Nb. Using push_back without reserving memory may cause
        // multiple re-allocation behind the scenes, this is expensive.
        send_values_data.reserve(total_ndata);
        send_column_indices_data.reserve(total_ndata);

        // Loop over all the processors to "flat pack" the data for sending.
        for (unsigned rank = 0; rank < nproc; rank++)
        {
          // Set the offset for the current processor
          // This only has to be done once for both values and column indices.
          send_displacement[rank] = send_values_data.size();

          // Don't bother to do anything if
          // the processor in the loop is the current processor.
          if (rank != my_rank)
          {
            // Put the values into the send data vector.
            // n_data is the same for both values and column indices.
            unsigned n_data = values_to_send[rank].size();
            for (unsigned j = 0; j < n_data; j++)
            {
              send_values_data.push_back(values_to_send[rank][j]);
              send_column_indices_data.push_back(
                column_indices_to_send[rank][j]);
            } // for
          } // if rank != my_rank

          // Find the number of data to be added to the vector.
          // send_n is the same for both values and column indices.
          send_n[rank] = send_values_data.size() - send_displacement[rank];
        } // loop over processors

        if (enable_timing)
        {
          double t_flat_pack_finish = TimingHelpers::timer();
          double t_flat_pack_time = t_flat_pack_finish - t_flat_pack_start;
          oomph_info << "t_flat_pack_time: " << t_flat_pack_time << std::endl;
        }

        double t_sendn_start;
        if (enable_timing) t_sendn_start = TimingHelpers::timer();

        // Strorage for the number of data to be received from each processor
        Vector<int> receive_n(nproc, 0);

        MPI_Alltoall(&send_n[0],
                     1,
                     MPI_INT,
                     &receive_n[0],
                     1,
                     MPI_INT,
                     comm_pt->mpi_comm());

        if (enable_timing)
        {
          double t_sendn_finish = TimingHelpers::timer();
          double t_sendn_time = t_sendn_finish - t_sendn_start;
          oomph_info << "t_sendn_time: " << t_sendn_time << std::endl;
        }


        // Prepare the data to be received
        // by working out the displacement from the received data
        // receive_displacement is the same for both values and column indices.
        Vector<int> receive_displacement(nproc, 0);
        int receive_data_count = 0;
        for (unsigned rank = 0; rank < nproc; rank++)
        {
          receive_displacement[rank] = receive_data_count;
          receive_data_count += receive_n[rank];
        }

        // Now resize the receive buffer for all data from all processors.
        // Make sure that it has a size of at least one.
        if (receive_data_count == 0)
        {
          receive_data_count++;
        }
        Vector<double> receive_values_data(receive_data_count);
        Vector<unsigned> receive_column_indices_data(receive_data_count);

        // Make sure that the send buffer has size at least one
        // so that we don't get a segmentation fault.
        if (send_values_data.size() == 0)
        {
          send_values_data.resize(1);
        }

        double t_send_data_start;
        if (enable_timing) t_send_data_start = TimingHelpers::timer();

        // Now send the data between all processors
        MPI_Alltoallv(&send_values_data[0],
                      &send_n[0],
                      &send_displacement[0],
                      MPI_DOUBLE,
                      &receive_values_data[0],
                      &receive_n[0],
                      &receive_displacement[0],
                      MPI_DOUBLE,
                      comm_pt->mpi_comm());

        // Now send the data between all processors
        MPI_Alltoallv(&send_column_indices_data[0],
                      &send_n[0],
                      &send_displacement[0],
                      MPI_UNSIGNED,
                      &receive_column_indices_data[0],
                      &receive_n[0],
                      &receive_displacement[0],
                      MPI_UNSIGNED,
                      comm_pt->mpi_comm());

        if (enable_timing)
        {
          double t_send_data_finish = TimingHelpers::timer();
          double t_send_data_time = t_send_data_finish - t_send_data_start;
          oomph_info << "t_send_data_time: " << t_send_data_time << std::endl;
        }

        // All the rows for this processor are stored in:
        // from other processors:
        // receive_column_indices_data and receive_values_data
        // from this processor:
        // column_indices_to_send[my_rank] and values_to_send[my_rank]
        //
        // They are in some order determined by the distribution.
        // We need to re-arrange them. To do this, we do some pre-processing.

        // nrow_local for this processor.
        unsigned res_nrow_local = res_distribution_pt->nrow_local();

        // Per row, store:
        // 1) where this row came from, 0 - this proc, 1 - other procs.
        // 2) the nnz,
        // 3) the offset - where the values/columns in the receive data vectors
        //                 begins. This is different from the offset of where
        //                 the data from a certain processor starts.
        Vector<Vector<unsigned>> value_column_locations(res_nrow_local,
                                                        Vector<unsigned>(3, 0));

        // Store the local nnz so we can reserve space for
        // the values and column indices.
        unsigned long res_nnz_local = 0;

        double t_locations_start;
        if (enable_timing) t_locations_start = TimingHelpers::timer();

        // Loop through the data currently on this processor.
        unsigned location_i = 0;
        unsigned my_column_indices_to_send_size =
          column_indices_to_send[my_rank].size();
        while (location_i < my_column_indices_to_send_size)
        {
          unsigned current_local_eqn =
            column_indices_to_send[my_rank][location_i++];

          unsigned current_nnz = column_indices_to_send[my_rank][location_i++];

          // No need to fill [*][0] with 0 since it is already initialised to 0.

          // Store the nnz.
          value_column_locations[current_local_eqn][1] = current_nnz;

          // Also increment the res_local_nnz
          res_nnz_local += current_nnz;

          // Store the offset.
          value_column_locations[current_local_eqn][2] = location_i;

          // Update the location_i so it starts at the next row.
          location_i += current_nnz;
        }

        // Loop through the data from different processors.

        // Check to see if data has been received.
        bool data_has_been_received = false;
        unsigned send_rank = 0;
        while (send_rank < nproc)
        {
          if (receive_n[send_rank] > 0)
          {
            data_has_been_received = true;
            break;
          }
          send_rank++;
        }

        location_i = 0; // start at 0.
        if (data_has_been_received)
        {
          unsigned receive_column_indices_data_size =
            receive_column_indices_data.size();
          while (location_i < receive_column_indices_data_size)
          {
            unsigned current_local_eqn =
              receive_column_indices_data[location_i++];
            unsigned current_nnz = receive_column_indices_data[location_i++];

            // These comes from other processors.
            value_column_locations[current_local_eqn][0] = 1;

            // Store the nnz.
            value_column_locations[current_local_eqn][1] = current_nnz;

            // Also increment the res_local_nnz
            res_nnz_local += current_nnz;

            // Store the offset.
            value_column_locations[current_local_eqn][2] = location_i;

            // Update the location_i so it starts at the next row.
            location_i += current_nnz;
          }
        }

        if (enable_timing)
        {
          double t_locations_finish = TimingHelpers::timer();
          double t_locations_time = t_locations_finish - t_locations_start;
          oomph_info << "t_locations_time: " << t_locations_time << std::endl;
        }

        double t_fillvecs_start;
        if (enable_timing) t_fillvecs_start = TimingHelpers::timer();

        // Now loop through the locations and store the values
        // the column indices in the correct order.
        Vector<int> res_column_indices;
        Vector<double> res_values;
        Vector<int> res_row_start;

        res_column_indices.reserve(res_nnz_local);
        res_values.reserve(res_nnz_local);
        res_row_start.reserve(res_nrow_local + 1);

        // Running sum of nnz for the row_start. Must be int because
        // res_row_start is templated with int.
        int nnz_running_sum = 0;

        // Now insert the rows.
        for (unsigned local_row_i = 0; local_row_i < res_nrow_local;
             local_row_i++)
        {
          // Fill the res_row_start with the nnz so far.
          res_row_start.push_back(nnz_running_sum);

          bool data_is_from_other_proc =
            bool(value_column_locations[local_row_i][0]);

          unsigned row_i_nnz = value_column_locations[local_row_i][1];

          unsigned row_i_offset = value_column_locations[local_row_i][2];

          if (data_is_from_other_proc)
          {
            // Insert range [offset, offset+nnz) from
            // receive_column_indices_data and receive_values_data into
            // res_column_indices and res_values respectively.
            res_column_indices.insert(
              res_column_indices.end(),
              receive_column_indices_data.begin() + row_i_offset,
              receive_column_indices_data.begin() + row_i_offset + row_i_nnz);

            res_values.insert(res_values.end(),
                              receive_values_data.begin() + row_i_offset,
                              receive_values_data.begin() + row_i_offset +
                                row_i_nnz);
          }
          else
          {
            res_column_indices.insert(res_column_indices.end(),
                                      column_indices_to_send[my_rank].begin() +
                                        row_i_offset,
                                      column_indices_to_send[my_rank].begin() +
                                        row_i_offset + row_i_nnz);

            res_values.insert(res_values.end(),
                              values_to_send[my_rank].begin() + row_i_offset,
                              values_to_send[my_rank].begin() + row_i_offset +
                                row_i_nnz);
          }

          // Update the running sum of nnz
          nnz_running_sum += row_i_nnz;
        }

        // Insert the last row_start value
        res_row_start.push_back(res_nnz_local);

        if (enable_timing)
        {
          double t_fillvecs_finish = TimingHelpers::timer();
          double t_fillvecs_time = t_fillvecs_finish - t_fillvecs_start;
          oomph_info << "t_fillvecs_time: " << t_fillvecs_time << std::endl;
        }

        double t_buildres_start;
        if (enable_timing) t_buildres_start = TimingHelpers::timer();

        // build the matrix.
        result_matrix.build(
          res_ncol, res_values, res_column_indices, res_row_start);

        if (enable_timing)
        {
          double t_buildres_finish = TimingHelpers::timer();
          double t_buildres_time = t_buildres_finish - t_buildres_start;
          oomph_info << "t_buildres_time: " << t_buildres_time << std::endl;
        }
        //  */
#endif
      }
    }

    //============================================================================
    ///  Concatenate CRDoubleMatrix matrices.
    ///
    /// The Vector row_distribution_pt contains the LinearAlgebraDistribution
    /// of each block row.
    /// The Vector col_distribution_pt contains the LinearAlgebraDistribution
    /// of each block column.
    /// The DenseMatrix matrix_pt contains pointers to the CRDoubleMatrices
    /// to concatenate.
    /// The CRDoubleMatrix result_matrix is the result matrix.
    ///
    /// The result matrix is a permutation of the sub matrices such that the
    /// data stays on the same processor when the result matrix is built, there
    /// is no communication between processors. Thus the block structure of the
    /// sub matrices are NOT preserved in the result matrix. The rows are
    /// block-permuted, defined by the concatenation of the distributions in
    /// row_distribution_pt. Similarly, the columns are block-permuted, defined
    /// by the concatenation of the distributions in col_distribution_pt. For
    /// more details on the block-permutation, see
    /// LinearAlgebraDistributionHelpers::concatenate(...).
    ///
    /// If one wishes to preserve the block structure of the sub matrices in the
    /// result matrix, consider using CRDoubleMatrixHelpers::concatenate(...),
    /// which uses communication between processors to ensure that the block
    /// structure of the sub matrices are preserved.
    ///
    /// The matrix manipulation functions
    /// CRDoubleMatrixHelpers::concatenate(...) and
    /// CRDoubleMatrixHelpers::concatenate_without_communication(...)
    /// are analogous to the Vector manipulation functions
    /// DoubleVectorHelpers::concatenate(...) and
    /// DoubleVectorHelpers::concatenate_without_communication(...).
    /// Please look at the DoubleVector functions for an illustration of the
    /// differences between concatenate(...) and
    /// concatenate_without_communication(...).
    ///
    /// Distribution of the result matrix:
    /// If the result matrix does not have a distribution built, then it will be
    /// given a distribution built from the concatenation of the distributions
    /// from row_distribution_pt, see
    /// LinearAlgebraDistributionHelpers::concatenate(...) for more detail.
    /// Otherwise we use the existing distribution.
    /// If there is an existing distribution then it must be the same as the
    /// distribution from the concatenation of row distributions as described
    /// above.
    /// Why don't we always compute the distribution "on the fly"?
    /// Because a non-uniform distribution requires communication.
    /// All block preconditioner distributions are concatenations of the
    /// distributions of the individual blocks.
    //============================================================================
    void concatenate_without_communication(
      const Vector<LinearAlgebraDistribution*>& row_distribution_pt,
      const Vector<LinearAlgebraDistribution*>& col_distribution_pt,
      const DenseMatrix<CRDoubleMatrix*>& matrix_pt,
      CRDoubleMatrix& result_matrix)
    {
      // The number of block rows and block columns.
      unsigned matrix_nrow = matrix_pt.nrow();
      unsigned matrix_ncol = matrix_pt.ncol();

      // PARANOID checks involving in matrices and block_distribution only.
      // PARANOID checks involving the result matrix will come later since
      // we have to create the result matrix distribution from the in
      // distribution if it does not already exist.
#ifdef PARANOID

      // Are there matrices to concatenate?
      if (matrix_nrow == 0 || matrix_ncol == 0)
      {
        std::ostringstream error_message;
        error_message << "There are no matrices to concatenate.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Does this matrix need concatenating?
      if ((matrix_nrow == 1) && (matrix_ncol == 1))
      {
        std::ostringstream warning_message;
        warning_message << "There is only one matrix to concatenate...\n"
                        << "This does not require concatenating...\n";
        OomphLibWarning(warning_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
      }


      // The distribution for each block row is stored in row_distribution_pt.
      // So the number of distributions in row_distribution_pt must be the
      // same as matrix_nrow.
      if (matrix_nrow != row_distribution_pt.size())
      {
        std::ostringstream error_message;
        error_message << "The number of row distributions must be the same as\n"
                      << "the number of block rows.";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // The number of distributions for the columns must match the number of
      // block columns.
      if (matrix_ncol != col_distribution_pt.size())
      {
        std::ostringstream error_message;
        error_message
          << "The number of column distributions must be the same as\n"
          << "the number of block columns.";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Check that all pointers in row_distribution_pt is not null.
      for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
      {
        if (row_distribution_pt[block_row_i] == 0)
        {
          std::ostringstream error_message;
          error_message << "The row distribution pointer in position "
                        << block_row_i << " is null.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Check that all pointers in row_distribution_pt is not null.
      for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
      {
        if (col_distribution_pt[block_col_i] == 0)
        {
          std::ostringstream error_message;
          error_message << "The column distribution pointer in position "
                        << block_col_i << " is null.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Check that all distributions are built.
      // First the row distributions
      for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
      {
        if (!row_distribution_pt[block_row_i]->built())
        {
          std::ostringstream error_message;
          error_message << "The distribution pointer in position "
                        << block_row_i << " is not built.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
      // Now the column distributions
      for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
      {
        if (!col_distribution_pt[block_col_i]->built())
        {
          std::ostringstream error_message;
          error_message << "The distribution pointer in position "
                        << block_col_i << " is not built.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Check that all communicators in row_distribution_pt are the same.
      const OomphCommunicator first_row_comm =
        *(row_distribution_pt[0]->communicator_pt());

      for (unsigned block_row_i = 1; block_row_i < matrix_nrow; block_row_i++)
      {
        const OomphCommunicator current_comm =
          *(row_distribution_pt[block_row_i]->communicator_pt());

        if (first_row_comm != current_comm)
        {
          std::ostringstream error_message;
          error_message
            << "The communicator from the row distribution in position "
            << block_row_i << " is not the same as the first "
            << "communicator from row_distribution_pt";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Check that all communicators in col_distribution_pt are the same as the
      // first row communicator from above.
      for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
      {
        const OomphCommunicator current_comm =
          *(col_distribution_pt[block_col_i]->communicator_pt());

        if (first_row_comm != current_comm)
        {
          std::ostringstream error_message;
          error_message
            << "The communicator from the col distribution in position "
            << block_col_i << " is not the same as the first "
            << "communicator from row_distribution_pt";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Are all sub matrices built? If the matrix_pt is not null, make sure
      // that it is built.
      for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
      {
        for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
        {
          if (matrix_pt(block_row_i, block_col_i) != 0 &&
              !(matrix_pt(block_row_i, block_col_i)->built()))
          {
            std::ostringstream error_message;
            error_message << "The sub matrix_pt(" << block_row_i << ","
                          << block_col_i << ")\n"
                          << "is not built.\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // For the matrices which are built, do they have the same communicator as
      // the first communicator from row_distribution_pt?
      for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
      {
        for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
        {
          if (matrix_pt(block_row_i, block_col_i) != 0)
          {
            const OomphCommunicator current_comm =
              *(matrix_pt(block_row_i, block_col_i)
                  ->distribution_pt()
                  ->communicator_pt());
            if (first_row_comm != current_comm)
            {
              std::ostringstream error_message;
              error_message
                << "The sub matrix_pt(" << block_row_i << "," << block_col_i
                << ")\n"
                << "does not have the same communicator pointer as those in\n"
                << "(row|col)_distribution_pt.\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }

      // Do all dimensions of sub matrices "make sense"?
      // Compare the number of rows of each block matrix in a block row.
      for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
      {
        // Use the first column to compare against the rest.
        unsigned long current_block_nrow =
          row_distribution_pt[block_row_i]->nrow();

        // Compare against columns 0 to matrix_ncol - 1
        for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
        {
          // Perform the check if the matrix_pt is not null.
          if (matrix_pt(block_row_i, block_col_i) != 0)
          {
            // Get the nrow for this sub block.
            unsigned long subblock_nrow =
              matrix_pt(block_row_i, block_col_i)->nrow();

            if (current_block_nrow != subblock_nrow)
            {
              std::ostringstream error_message;
              error_message << "The sub matrix (" << block_row_i << ","
                            << block_col_i << ")\n"
                            << "requires nrow = " << current_block_nrow
                            << ", but has nrow = " << subblock_nrow << ".\n"
                            << "Either the row_distribution_pt is incorrect or "
                            << "the sub matrices are incorrect.\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }

      // Compare the number of columns of each block matrix in a block column.
      for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
      {
        // Get the current block ncol from the linear algebra distribution.
        // Note that we assume that the dimensions are symmetrical.
        unsigned current_block_ncol = col_distribution_pt[block_col_i]->nrow();

        for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
        {
          if (matrix_pt(block_row_i, block_col_i) != 0)
          {
            // Get the ncol for this sub block.
            unsigned subblock_ncol =
              matrix_pt(block_row_i, block_col_i)->ncol();

            if (current_block_ncol != subblock_ncol)
            {
              std::ostringstream error_message;
              error_message << "The sub matrix (" << block_row_i << ","
                            << block_col_i << ")\n"
                            << "requires ncol = " << current_block_ncol
                            << ", but has ncol = " << subblock_ncol << ".\n"
                            << "Either the col_distribution_pt is incorrect or "
                            << "the sub matrices are incorrect.\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }

      // Ensure that the distributions for all sub matrices in the same block
      // row are the same. This is because we permute the row across several
      // matrices.

      // Loop through each block row.
      for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
      {
        // Get the distribution from the first block in this row.
        LinearAlgebraDistribution* block_row_distribution_pt =
          row_distribution_pt[block_row_i];

        // Loop through the block columns
        for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
        {
          if (matrix_pt(block_row_i, block_col_i) != 0)
          {
            // Get the distribution for this block.
            LinearAlgebraDistribution* current_block_distribution_pt =
              matrix_pt(block_row_i, block_col_i)->distribution_pt();

            // Ensure that the in matrices is a square block matrix.
            if ((*block_row_distribution_pt) !=
                (*current_block_distribution_pt))
            {
              std::ostringstream error_message;
              error_message
                << "Sub block(" << block_row_i << "," << block_col_i << ")"
                << "does not have the same distributoin as the first"
                << "block in this block row.\n"
                << "All distributions on a block row must be the same"
                << "for this function to concatenate matrices.\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }
#endif

      // The communicator pointer from the first row_distribution_pt
      const OomphCommunicator* const comm_pt =
        row_distribution_pt[0]->communicator_pt();

      // Renamed for so it makes more sense.
      unsigned nblock_row = matrix_nrow;

      // If the result matrix does not have a distribution, then we concatenate
      // the sub matrix distributions.
      if (!result_matrix.distribution_pt()->built())
      {
        // The result distribution
        LinearAlgebraDistribution tmp_distribution;
        LinearAlgebraDistributionHelpers::concatenate(row_distribution_pt,
                                                      tmp_distribution);

        result_matrix.build(&tmp_distribution);
      }
      else
      // A distribution is supplied for the result matrix.
      {
#ifdef PARANOID
        // Check that the result distribution is a concatenation of the
        // distributions of the sub matrices.

        LinearAlgebraDistribution wanted_distribution;

        LinearAlgebraDistributionHelpers::concatenate(row_distribution_pt,
                                                      wanted_distribution);

        if (*(result_matrix.distribution_pt()) != wanted_distribution)
        {
          std::ostringstream error_message;
          error_message
            << "The result distribution is not correct.\n"
            << "Please call the function without a result\n"
            << "distribution (clear the result matrix) or check the\n"
            << "distribution of the result matrix.\n"
            << "The result distribution must be the same as the one \n"
            << "created by\n"
            << "LinearAlgebraDistributionHelpers::concatenate(...)";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }

      // The rest of the paranoid checks.
#ifdef PARANOID

      // Make sure that the communicator from the result matrix is the same as
      // all the others. This test is redundant if this function created the
      // result matrix distribution, since then it is guaranteed that the
      // communicators are the same.
      {
        // Communicator from the result matrix.
        const OomphCommunicator res_comm =
          *(result_matrix.distribution_pt()->communicator_pt());

        // Is the result communicator pointer the same as the others?
        // Since we have already tested the others, we only need to compare
        // against one of them. Say the first communicator from
        // row_distribution_pt.
        const OomphCommunicator first_comm =
          *(row_distribution_pt[0]->communicator_pt());

        if (res_comm != first_comm)
        {
          std::ostringstream error_message;
          error_message << "The OomphCommunicator of the result matrix is not "
                           "the same as the "
                        << "others!";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Are all the distributed boolean the same? This only applies if we have
      // more than one processor. If there is only one processor, then it does
      // not matter if it is distributed or not - they are conceptually the
      // same.
      if (comm_pt->nproc() != 1)
      {
        // Compare distributed for sub matrices (against the result matrix).
        const bool res_distributed = result_matrix.distributed();

        // Loop over all sub blocks.
        for (unsigned block_row_i = 0; block_row_i < matrix_nrow; block_row_i++)
        {
          for (unsigned block_col_i = 0; block_col_i < matrix_ncol;
               block_col_i++)
          {
            if (matrix_pt(block_row_i, block_col_i) != 0)
            {
              const bool another_distributed =
                matrix_pt(block_row_i, block_col_i)->distributed();

              if (res_distributed != another_distributed)
              {
                std::ostringstream error_message;
                error_message << "The distributed boolean of the sub matrix ("
                              << block_row_i << "," << block_col_i << ")\n"
                              << "is not the same as the result matrix. \n";
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            }
          }
        }

        // Do this test for row_distribution_pt
        const bool first_row_distribution_distributed =
          row_distribution_pt[0]->distributed();

        for (unsigned block_row_i = 1; block_row_i < matrix_nrow; block_row_i++)
        {
          const bool another_distributed =
            row_distribution_pt[block_row_i]->distributed();

          if (first_row_distribution_distributed != another_distributed)
          {
            std::ostringstream error_message;
            error_message
              << "The distributed boolean of row_distribution_pt["
              << block_row_i << "]\n"
              << "is not the same as the one from row_distribution_pt[0]. \n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }

        // Repeat for col_distribution_pt
        for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
        {
          const bool another_distributed =
            col_distribution_pt[block_col_i]->distributed();

          if (first_row_distribution_distributed != another_distributed)
          {
            std::ostringstream error_message;
            error_message
              << "The distributed boolean of col_distribution_pt["
              << block_col_i << "]\n"
              << "is not the same as the one from row_distribution_pt[0]. \n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
#endif

      ////////////////// END OF PARANOID TESTS
      ///////////////////////////////////////

      // The number of processors.
      unsigned nproc = comm_pt->nproc();

      // Cache the result distribution pointer for convenience.
      LinearAlgebraDistribution* res_distribution_pt =
        result_matrix.distribution_pt();

      // nrow_local for the result matrix
      unsigned res_nrow_local = res_distribution_pt->nrow_local();

      // renamed for readability.
      unsigned nblock_col = matrix_ncol;

      // construct the block offset
      //  DenseMatrix<unsigned> col_offset(nproc,nblock_col,0);
      std::vector<std::vector<unsigned>> col_offset(
        nproc, std::vector<unsigned>(nblock_col));
      unsigned off = 0;
      for (unsigned proc_i = 0; proc_i < nproc; proc_i++)
      {
        for (unsigned block_i = 0; block_i < nblock_col; block_i++)
        {
          col_offset[proc_i][block_i] = off;
          off += col_distribution_pt[block_i]->nrow_local(proc_i);
        }
      }

      // Do some pre-processing for the processor number a global row number is
      // on. This is required when permuting the column entries.
      // We need to do this for each distribution, so we have a vector of
      // vectors. First index corresponds to the distribution, the second is
      // the processor number.
      std::vector<std::vector<unsigned>> p_for_rows(nblock_col,
                                                    std::vector<unsigned>());
      // initialise 2D vector
      for (unsigned blocki = 0; blocki < nblock_col; blocki++)
      {
        int blockinrow = col_distribution_pt[blocki]->nrow();
        p_for_rows[blocki].resize(blockinrow);
        // FOR each global index in the block, work out the corresponding proc.
        for (int rowi = 0; rowi < blockinrow; rowi++)
        {
          unsigned p = 0;
          int b_first_row = col_distribution_pt[blocki]->first_row(p);
          int b_nrow_local = col_distribution_pt[blocki]->nrow_local(p);

          while (rowi < b_first_row || rowi >= b_nrow_local + b_first_row)
          {
            p++;
            b_first_row = col_distribution_pt[blocki]->first_row(p);
            b_nrow_local = col_distribution_pt[blocki]->nrow_local(p);
          }
          p_for_rows[blocki][rowi] = p;
        }
      }

      // determine nnz of all blocks on this processor only.
      // This is used to create storage space.
      unsigned long res_nnz = 0;
      for (unsigned row_i = 0; row_i < nblock_row; row_i++)
      {
        for (unsigned col_i = 0; col_i < nblock_col; col_i++)
        {
          if (matrix_pt(row_i, col_i) != 0)
          {
            res_nnz += matrix_pt(row_i, col_i)->nnz();
          }
        }
      }

      // My rank
      //    unsigned my_rank = comm_pt->my_rank();
      //    my_rank = my_rank;

      // Turn the above into a string.
      //      std::ostringstream myrankstream;
      //      myrankstream << "THISDOESNOTHINGnp" << my_rank << std::endl;
      //      std::string myrankstring = myrankstream.str();


      // CALLGRIND_ZERO_STATS;
      // CALLGRIND_START_INSTRUMENTATION;

      // storage for the result matrix.
      int* res_row_start = new int[res_nrow_local + 1];
      int* res_column_index = new int[res_nnz];
      double* res_value = new double[res_nnz];

      // initialise the zero-th entry
      res_row_start[0] = 0;

      // loop over the block rows
      unsigned long res_i = 0; // index for the result matrix.
      unsigned long res_row_i = 0; // index for the row
      for (unsigned i = 0; i < nblock_row; i++)
      {
        // loop over the rows of the current block local rows.
        unsigned block_nrow = row_distribution_pt[i]->nrow_local();
        for (unsigned k = 0; k < block_nrow; k++)
        {
          // initialise res_row_start
          res_row_start[res_row_i + 1] = res_row_start[res_row_i];

          // Loop over the block columns
          for (unsigned j = 0; j < nblock_col; j++)
          {
            // if block(i,j) pointer is not null then
            if (matrix_pt(i, j) != 0)
            {
              // get pointers for the elements in the current block
              int* b_row_start = matrix_pt(i, j)->row_start();
              int* b_column_index = matrix_pt(i, j)->column_index();
              double* b_value = matrix_pt(i, j)->value();

              // memcpy( &dst[dstIdx], &src[srcIdx], numElementsToCopy * sizeof(
              // Element ) );
              // no ele to copy
              int numEleToCopy = b_row_start[k + 1] - b_row_start[k];
              memcpy(res_value + res_i,
                     b_value + b_row_start[k],
                     numEleToCopy * sizeof(double));
              // Loop through the current local row.
              for (int l = b_row_start[k]; l < b_row_start[k + 1]; l++)
              {
                // if b_column_index[l] was a row index, what processor
                // would it be on
                //            unsigned p = col_distribution_pt[j]
                //              ->rank_of_global_row_map(b_column_index[l]);
                unsigned p = p_for_rows[j][b_column_index[l]];

                int b_first_row = col_distribution_pt[j]->first_row(p);
                //            int b_nrow_local =
                //            col_distribution_pt[j]->nrow_local(p);

                //            while (b_column_index[l] < b_first_row ||
                //                   b_column_index[l] >=
                //                   b_nrow_local+b_first_row)
                //             {
                //              p++;
                //              b_first_row =
                //              col_distribution_pt[j]->first_row(p);
                //              b_nrow_local =
                //              col_distribution_pt[j]->nrow_local(p);
                //             }

                // determine the local equation number in the block j/processor
                // p "column block"
                int eqn = b_column_index[l] - b_first_row;

                // add to the result matrix
                //            res_value[res_i] = b_value[l];
                res_column_index[res_i] = col_offset[p][j] + eqn;
                res_row_start[res_row_i + 1]++;
                res_i++;
              }
            }
          }

          // increment the row pt
          res_row_i++;
        }
      }
      // CALLGRIND_STOP_INSTRUMENTATION;
      // CALLGRIND_DUMP_STATS_AT(myrankstring.c_str());


      // Get the number of columns of the result matrix.
      unsigned res_ncol = 0;
      for (unsigned block_col_i = 0; block_col_i < matrix_ncol; block_col_i++)
      {
        res_ncol += col_distribution_pt[block_col_i]->nrow();
      }

      // Build the result matrix.
      result_matrix.build_without_copy(
        res_ncol, res_nnz, res_value, res_column_index, res_row_start);
    }


    //============================================================================
    ///  Concatenate CRDoubleMatrix matrices.
    /// This calls the other concatenate_without_communication(...) function,
    /// passing block_distribution_pt as both the row_distribution_pt and
    /// col_distribution_pt. This should only be called for block square
    /// matrices.
    //============================================================================
    void concatenate_without_communication(
      const Vector<LinearAlgebraDistribution*>& block_distribution_pt,
      const DenseMatrix<CRDoubleMatrix*>& matrix_pt,
      CRDoubleMatrix& result_matrix)
    {
#ifdef PARANOID
      // The number of block rows and block columns.
      unsigned matrix_nrow = matrix_pt.nrow();
      unsigned matrix_ncol = matrix_pt.ncol();

      // Are there matrices to concatenate?
      if (matrix_nrow == 0)
      {
        std::ostringstream error_message;
        error_message << "There are no matrices to concatenate.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Ensure that the sub matrices is a square block matrix.
      if (matrix_nrow != matrix_ncol)
      {
        std::ostringstream error_message;
        error_message
          << "The number of block rows and block columns\n"
          << "must be the same. Otherwise, call the other\n"
          << "concatenate_without_communication function, passing in\n"
          << "a Vector of distributions describing how to permute the\n"
          << "columns.";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      concatenate_without_communication(
        block_distribution_pt, block_distribution_pt, matrix_pt, result_matrix);
    }

  } // namespace CRDoubleMatrixHelpers

} // namespace oomph
