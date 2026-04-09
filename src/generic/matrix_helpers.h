#pragma once
#include "generic/matrices.h"

namespace oomph
{

  // ===========================================================================
  /// \short This class contains static functions, which can help with
  /// manipulating matrices.
  // ===========================================================================
  class MatrixHelpers
  {
  public:
    // =========================================================================
    /// \short Creates a formatted string representation of the matrix.
    /// Useful for debugging and logging.
    /// \param M The matrix to convert to a string.
    /// \param width The character width per number (default 12).
    /// \param prec The precision for floating point numbers (default 5).
    // =========================================================================
    static std::string format(const DenseMatrix<double>& M,
                              int width = 12,
                              int prec = 5)
    {
      std::ostringstream output;
      unsigned long nrow = M.nrow();
      unsigned long ncol = M.ncol();

      // Set fixed formatting for alignment
      output << std::scientific << std::setprecision(prec);

      for (unsigned long i = 0; i < nrow; i++)
      {
        // Optional: Add a visual bracket at start of row
        output << "| ";

        for (unsigned long j = 0; j < ncol; j++)
        {
          output << std::setw(width) << M(i, j);
          // Add a small spacer between columns if not using the last one
          if (j < ncol - 1)
          {
            output << " ";
          }
        }

        // Optional: Add a visual bracket at end of row
        output << " |" << std::endl;
      }

      return output.str();
    }

    // =========================================================================
    /// \short Multiplies two matrices
    // =========================================================================
    static void multiply(const DenseMatrix<double>& A,
                         const DenseMatrix<double>& B,
                         DenseMatrix<double>& result)
    {
      unsigned long rows_A = A.nrow();
      unsigned long cols_A = A.ncol();
      unsigned long rows_B = B.nrow();
      unsigned long cols_B = B.ncol();

#ifdef PARANOID
      // check matrix dimensions are compatible
      if (cols_A != rows_B)
      {
        std::ostringstream error_message;
        error_message
          << "Matrix dimensions incompatible for matrix-matrix multiplication. "
          << "Cols of A: " << cols_A << ", Rows of B: " << rows_B;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      result.resize(rows_A, cols_B, 0.0);

      for (unsigned long i = 0; i < rows_A; i++)
      {
        for (unsigned long k = 0; k < cols_A; k++)
        {
          double a_val = A(i, k);
          for (unsigned long j = 0; j < cols_B; j++)
          {
            result(i, j) += a_val * B(k, j);
          }
        }
      }
    }

    // =========================================================================
    /// \short Computes result = A^T B
    // =========================================================================
    static void transpose_multiply(const DenseMatrix<double>& A,
                                   const DenseMatrix<double>& B,
                                   DenseMatrix<double>& result)
    {
      unsigned long rows_A = A.nrow();
      unsigned long cols_A = A.ncol();
      unsigned long rows_B = B.nrow();
      unsigned long cols_B = B.ncol();

#ifdef PARANOID
      // check matrix dimensions are compatible
      if (rows_A != rows_B)
      {
        std::ostringstream error_message;
        error_message
          << "Matrix dimensions incompatible for matrix-matrix multiplication. "
          << "Rows of A: " << rows_A << ", Rows of B: " << rows_B;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      result.resize(cols_A, cols_B, 0.0);

      for (unsigned long k = 0; k < rows_A; k++)
      {
        for (unsigned long i = 0; i < cols_A; i++)
        {
          double a_val = A(k, i);
          for (unsigned long j = 0; j < cols_B; j++)
          {
            result(i, j) += a_val * B(k, j);
          }
        }
      }
    }

    // =========================================================================
    /// \short Computes result = A B^T
    // =========================================================================
    static void multiply_transpose(const DenseMatrix<double>& A,
                                   const DenseMatrix<double>& B,
                                   DenseMatrix<double>& result)
    {
      unsigned long rows_A = A.nrow();
      unsigned long cols_A = A.ncol();
      unsigned long rows_B = B.nrow();
      unsigned long cols_B = B.ncol();

#ifdef PARANOID
      // check matrix dimensions are compatible
      if (cols_A != cols_B)
      {
        std::ostringstream error_message;
        error_message
          << "Matrix dimensions incompatible for matrix-matrix multiplication. "
          << "Cols of A: " << cols_A << ", Cols of B: " << cols_B;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      result.resize(rows_A, rows_B, 0.0);

      for (unsigned long i = 0; i < rows_A; i++)
      {
        for (unsigned long j = 0; j < rows_B; j++)
        {
          double sum = 0.0;
          for (unsigned long k = 0; k < cols_A; k++)
          {
            sum += A(i, k) * B(j, k);
          }
          result(i, j) = sum;
        }
      }
    }

    // =========================================================================
    /// \short Computes the squared magnitude of a matrix:
    ///                                                |M|^2 = sum_i,j M(i, j)^2
    // =========================================================================
    static double magnitude_squared(const DenseMatrix<double>& M)
    {
      const unsigned long nrow = M.nrow();
      const unsigned long ncol = M.ncol();

      double squared_magnitude = 0.0;
      for (unsigned long i = 0; i < nrow; i++)
      {
        for (unsigned long j = 0; j < ncol; j++)
        {
          squared_magnitude += M(i, j) * M(i, j);
        }
      }

      return squared_magnitude;
    }

    // =========================================================================
    /// \short Computes the magnitude of a matrix: |M| = sqrt{sum_i,j M(i, j)^2}
    // =========================================================================
    static double magnitude(const DenseMatrix<double>& M)
    {
      return std::sqrt(magnitude_squared(M));
    }

    // =========================================================================
    /// \short Computes the dot product between two matrices:
    //                                         A dot B = sum_i,j A(i, j) B(i, j)
    // =========================================================================
    static double dot(const DenseMatrix<double>& A,
                      const DenseMatrix<double>& B)
    {
#ifdef PARANOID
      if (A.nrow() != B.nrow() || A.ncol() != B.ncol())
      {
        std::ostringstream error_message;
        error_message << "Matrix dimensions incompatible for dot product. "
                      << "nrow, ncol for first matrix: " << A.nrow() << ", "
                      << A.ncol()
                      << " nrow, ncol for second matrix: " << B.nrow() << ", "
                      << B.ncol();

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      double sum = 0.0;
      for (unsigned long i = 0; i < A.nrow(); i++)
      {
        for (unsigned long j = 0; j < A.ncol(); j++)
        {
          sum += A(i, j) * B(i, j);
        }
      }

      return sum;
    }

    // =========================================================================
    /// \short Normalises a matrix inplace: hat{M} = M / |M|
    // =========================================================================
    static void normalise(DenseMatrix<double>& M)
    {
      double mag = magnitude(M);
      if (mag == 0)
      {
        return;
      }

      double inv_mag = 1.0 / mag;

      unsigned long nrow = M.nrow();
      unsigned long ncol = M.ncol();
      for (unsigned long i = 0; i < nrow; i++)
      {
        for (unsigned long j = 0; j < ncol; j++)
        {
          M(i, j) *= inv_mag;
        }
      }
    }

    // =========================================================================
    /// \short Symmetrizes a matrix: Msym = 0.5 (M + M^T)
    // =========================================================================
    static void symmetrize(const DenseMatrix<double>& M,
                           DenseMatrix<double>& Msym)
    {
      unsigned long nrow = M.nrow();
      unsigned long ncol = M.ncol();

#ifdef PARANOID
      if (nrow != ncol)
      {
        throw OomphLibError("Matrix must be square to be symmetrized.",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      Msym.resize(nrow, ncol, 0.0);

      for (unsigned long i = 0; i < nrow; i++)
      {
        for (unsigned long j = 0; j < ncol; j++)
        {
          Msym(i, j) = 0.5 * (M(i, j) + M(j, i));
        }
      }
    }

    // =========================================================================
    /// \short Anitsymmetrizes a matrix: Mskew = 0.5 (M - M^T)
    // =========================================================================
    static void antisymmetrize(const DenseMatrix<double>& M,
                               DenseMatrix<double>& Mskew)
    {
      unsigned long nrow = M.nrow();
      unsigned long ncol = M.ncol();

#ifdef PARANOID
      if (nrow != ncol)
      {
        throw OomphLibError("Matrix must be square to be antisymmetrized.",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      Mskew.resize(nrow, ncol, 0.0);

      for (unsigned long i = 0; i < nrow; i++)
      {
        for (unsigned long j = 0; j < ncol; j++)
        {
          Mskew(i, j) = 0.5 * (M(i, j) - M(j, i));
        }
      }
    }

    // =========================================================================
    /// \short A convenience alias for dot
    // =========================================================================
    static double reduce(const DenseMatrix<double>& A,
                         const DenseMatrix<double>& B)
    {
      return dot(A, B);
    }

    // =========================================================================
    /// \short Computes a tensor matrix reduction / multiplication:
    ///   M_ret(i, j) = T(i, j, k, l) M(k, l)
    // =========================================================================
    static void reduce(const RankFourTensor<double>& T,
                       const DenseMatrix<double>& M,
                       DenseMatrix<double>& M_ret)
    {
      unsigned long n1 = T.nindex1();
      unsigned long n2 = T.nindex2();
      unsigned long n3 = T.nindex3();
      unsigned long n4 = T.nindex4();

#ifdef PARANOID
      if (M.nrow() != n3 || M.ncol() != n4)
      {
        std::ostringstream error_message;
        error_message << "Matrix dimensions incompatible for tensor reduction. "
                      << "Expected matrix of size " << n3 << "x" << n4
                      << ", but got " << M.nrow() << "x" << M.ncol() << ".";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      M_ret.resize(n1, n2, 0.0);

      for (unsigned long i = 0; i < n1; i++)
      {
        for (unsigned long j = 0; j < n2; j++)
        {
          double sum = 0.0;

          // Contract over the 3rd and 4th indices of the Tensor
          // against the rows and columns of the Matrix
          for (unsigned long k = 0; k < n3; k++)
          {
            for (unsigned long l = 0; l < n4; l++)
            {
              // T_ijkl * M_kl
              sum += T(i, j, k, l) * M(k, l);
            }
          }

          M_ret(i, j) = sum;
        }
      }
    }

    // =========================================================================
    /// \short Inverts a matrix. Note: This is only implemented for DIM <= 3
    // =========================================================================
    template<unsigned DIM>
    static double invert_matrix(const DenseMatrix<double>& M,
                                DenseMatrix<double>& inverse_M)
    {
      throw OomphLibError("invert_matrix is not implemented for DIM > 3",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
      return 0;
    }

    // =========================================================================
    /// \short Computes the determinant of a matrix. Note: This is only
    /// implemented for DIM <= 3
    // =========================================================================
    template<unsigned DIM>
    static double determinant(const DenseMatrix<double>& M)
    {
      throw OomphLibError("determinant is not implemented for DIM > 3",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
      return 0;
    }

    // Used to check if indices make sense when accessing a matrix which is of
    // size DIM x DIM. Should be optimised out when not compiled with PARANOID
    template<unsigned DIM>
    static void check_matrix_indices(const unsigned& i, const unsigned& j)
    {
      // If the indices don't make sense then throw an error
#ifdef PARANOID
      if (i >= DIM || j >= DIM)
      {
        throw OomphLibError("Indices must be less than dimension size",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }
  };

  template<>
  inline double MatrixHelpers::invert_matrix<1>(const DenseMatrix<double>& M,
                                                DenseMatrix<double>& inverse_M)
  {
    // Calculate the determinant of the matrix
    const double det = M(0, 0);

    //     // Report if Matrix is singular or negative
    // #ifdef PARANOID
    //     check_jacobian(det);
    // #endif

    // Calculate the inverse --- trivial in 1D
    inverse_M(0, 0) = 1.0 / M(0, 0);

    // Return the determinant
    return (det);
  }

  //===========================================================================
  /// Two-d specialisation of function to calculate inverse of M mapping
  //===========================================================================
  template<>
  inline double MatrixHelpers::invert_matrix<2>(const DenseMatrix<double>& M,
                                                DenseMatrix<double>& inverse_M)
  {
    // Calculate the determinant of the matrix
    const double det = M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);

    // Report if Matrix is singular or negative
    // #ifdef PARANOID
    //     check_jacobian(det);
    // #endif

    // Calculate the inverset of the 2x2 matrix
    inverse_M(0, 0) = M(1, 1) / det;
    inverse_M(0, 1) = -M(0, 1) / det;
    inverse_M(1, 0) = -M(1, 0) / det;
    inverse_M(1, 1) = M(0, 0) / det;

    // Return the M
    return (det);
  }

  //=============================================================================
  /// Three-d specialisation of function to calculate inverse of M
  /// mapping
  //=============================================================================
  template<>
  inline double MatrixHelpers::invert_matrix<3>(const DenseMatrix<double>& M,
                                                DenseMatrix<double>& inverse_M)
  {
    // Calculate the determinant of the matrix
    const double det =
      M(0, 0) * M(1, 1) * M(2, 2) + M(0, 1) * M(1, 2) * M(2, 0) +
      M(0, 2) * M(1, 0) * M(2, 1) - M(0, 0) * M(1, 2) * M(2, 1) -
      M(0, 1) * M(1, 0) * M(2, 2) - M(0, 2) * M(1, 1) * M(2, 0);

    //     // Report if Matrix is singular or negative
    // #ifdef PARANOID
    //     check_jacobian(det);
    // #endif

    // Calculate the inverse of the 3x3 matrix
    inverse_M(0, 0) = (M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1)) / det;
    inverse_M(0, 1) = -(M(0, 1) * M(2, 2) - M(0, 2) * M(2, 1)) / det;
    inverse_M(0, 2) = (M(0, 1) * M(1, 2) - M(0, 2) * M(1, 1)) / det;
    inverse_M(1, 0) = -(M(1, 0) * M(2, 2) - M(1, 2) * M(2, 0)) / det;
    inverse_M(1, 1) = (M(0, 0) * M(2, 2) - M(0, 2) * M(2, 0)) / det;
    inverse_M(1, 2) = -(M(0, 0) * M(1, 2) - M(0, 2) * M(1, 0)) / det;
    inverse_M(2, 0) = (M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0)) / det;
    inverse_M(2, 1) = -(M(0, 0) * M(2, 1) - M(0, 1) * M(2, 0)) / det;
    inverse_M(2, 2) = (M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0)) / det;

    // Return the determinant
    return (det);
  }

  //===========================================================================
  /// One-d specialisation of function to calculate determinant
  //===========================================================================
  template<>
  inline double MatrixHelpers::determinant<1>(const DenseMatrix<double>& M)
  {
    return M(0, 0);
  }

  //===========================================================================
  /// Two-d specialisation of function to calculate determinant
  //===========================================================================
  template<>
  inline double MatrixHelpers::determinant<2>(const DenseMatrix<double>& M)
  {
    return M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
  }

  //===========================================================================
  /// Three-d specialisation of function to calculate determinant
  //===========================================================================
  template<>
  inline double MatrixHelpers::determinant<3>(const DenseMatrix<double>& M)
  {
    return M(0, 0) * M(1, 1) * M(2, 2) + M(0, 1) * M(1, 2) * M(2, 0) +
           M(0, 2) * M(1, 0) * M(2, 1) - M(0, 0) * M(1, 2) * M(2, 1) -
           M(0, 1) * M(1, 0) * M(2, 2) - M(0, 2) * M(1, 1) * M(2, 0);
  }

} // End of namespace oomph