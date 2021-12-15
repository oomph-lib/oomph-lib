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
#include "complex_matrices.h"

namespace oomph
{
  //============================================================================
  /// Complete LU solve (overwrites RHS with solution). This is the
  /// generic version which should not need to be over-written.
  //============================================================================
  void ComplexMatrixBase::solve(Vector<std::complex<double>>& rhs)
  {
#ifdef PARANOID
    // Check Matrix is square
    if (nrow() != ncol())
    {
      throw OomphLibError(
        "This matrix is not square, the matrix MUST be square!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
    // Check that size of rhs = nrow()
    if (rhs.size() != nrow())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The rhs vector is not the right size. It is "
                           << rhs.size() << ", it should be " << nrow()
                           << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Call the LU decomposition
    ludecompose();

    // Call the back substitution
    lubksub(rhs);
  }


  //============================================================================
  /// Complete LU solve (Nothing gets overwritten!). This generic
  /// version should never need to be overwritten
  //============================================================================
  void ComplexMatrixBase::solve(const Vector<std::complex<double>>& rhs,
                                Vector<std::complex<double>>& soln)
  {
    // Set the solution vector equal to the rhs
    // N.B. This won't work if we change to another vector format
    soln = rhs;

    // Overwrite the solution vector (rhs is unchanged)
    solve(soln);
  }


  //=======================================================================
  /// Delete the storage that has been allocated for the LU factors, if
  /// the matrix data is not itself being overwritten.
  //======================================================================
  void DenseComplexMatrix::delete_lu_factors()
  {
    // Clean up the LU factor storage, if it has been allocated
    // and it's not the same as the matrix storage
    if ((!Overwrite_matrix_storage) && (LU_factors != 0))
    {
      // Delete the pointer to the LU factors
      delete[] LU_factors;
      // Null out the vector
      LU_factors = 0;
    }
  }


  //=======================================================================
  /// Destructor clean up the LU factors that have been allocated
  //======================================================================
  DenseComplexMatrix::~DenseComplexMatrix()
  {
    // Delete the storage for the index vector
    delete Index;

    // Delete the pointer to the LU factors
    delete[] LU_factors;

    // Null out the vector
    LU_factors = 0;
  }

  //============================================================================
  /// LU decompose a matrix, over-writing it and recording the pivots
  /// numbers in the Index vector.
  /// Returns the sign of the determinant.
  //============================================================================
  int DenseComplexMatrix::ludecompose()
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

    // Small constant
    const double small_number = 1.0e-20;

    // If the Index vector has not already been allocated,
    // allocated storage for the index of permutations
    if (Index == 0)
    {
      Index = new Vector<long>;
    }

    // Resize the index vector to the correct length
    Index->resize(N);

    // Vector scaling stores the implicit scaling of each row
    Vector<double> scaling(N);

    // Integer to store the sign that must multiply the determinant as
    // a consequence of the row/column interchanges
    int signature = 1;

    // Loop over rows to get implicit scaling information
    for (unsigned long i = 0; i < N; i++)
    {
      double big = 0.0;
      for (unsigned long j = 0; j < M; j++)
      {
        double tmp = std::abs((*this)(i, j));
        if (tmp > big) big = tmp;
      }
      if (big == 0.0)
      {
        throw OomphLibError(
          "Singular Matrix", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      // Save the scaling
      scaling[i] = 1.0 / big;
    }

    // Firsly, we shall delete any previous LU storage.
    // If the user calls this function twice without changing the matrix
    // then it is their own inefficiency, not ours (this time).
    delete_lu_factors();

    // Now if we're saving the LU factors separately (the default).
    if (!Overwrite_matrix_storage)
    {
      // Allocate storage for the LU factors
      // Assign space for the n rows
      LU_factors = new std::complex<double>[N * M];

      // Now we know that memory has been allocated, copy over
      // the matrix values
      for (unsigned long i = 0; i < (N * M); i++)
      {
        LU_factors[i] = Matrixdata[i];
      }
    }
    // Otherwise the pointer to the LU_factors is the same as the
    // matrix data
    else
    {
      LU_factors = Matrixdata;
    }

    // Loop over columns
    for (unsigned long j = 0; j < M; j++)
    {
      // Initialise imax
      unsigned long imax = 0;

      for (unsigned long i = 0; i < j; i++)
      {
        std::complex<double> sum = LU_factors[M * i + j];
        for (unsigned long k = 0; k < i; k++)
        {
          sum -= LU_factors[M * i + k] * LU_factors[M * k + j];
        }
        LU_factors[M * i + j] = sum;
      }

      // Initialise search for largest pivot element
      double largest_entry = 0.0;
      for (unsigned long i = j; i < N; i++)
      {
        std::complex<double> sum = LU_factors[M * i + j];
        for (unsigned long k = 0; k < j; k++)
        {
          sum -= LU_factors[M * i + k] * LU_factors[M * k + j];
        }
        LU_factors[M * i + j] = sum;
        // Set temporary
        double tmp = scaling[i] * std::abs(sum);
        if (tmp >= largest_entry)
        {
          largest_entry = tmp;
          imax = i;
        }
      }

      // Test to see if we need to interchange rows
      if (j != imax)
      {
        for (unsigned long k = 0; k < N; k++)
        {
          std::complex<double> tmp = LU_factors[M * imax + k];
          LU_factors[M * imax + k] = LU_factors[M * j + k];
          LU_factors[M * j + k] = tmp;
        }
        // Change the parity of signature
        signature = -signature;
        // Interchange scale factor
        scaling[imax] = scaling[j];
      }

      // Set the index
      (*Index)[j] = imax;
      if (LU_factors[M * j + j] == 0.0)
      {
        LU_factors[M * j + j] = small_number;
      }
      // Divide by pivot element
      if (j != N - 1)
      {
        std::complex<double> tmp = 1.0 / LU_factors[M * j + j];
        for (unsigned long i = j + 1; i < N; i++)
        {
          LU_factors[M * i + j] *= tmp;
        }
      }

    } // End of loop over columns


    // Now multiply all the diagonal terms together to get the determinant
    // Note that we need to use the mantissa, exponent formulation to
    // avoid underflow errors. This is way more tedious in complex arithmetic.
    std::complex<double> determinant_mantissa(1.0, 0.0);
    std::complex<int> determinant_exponent(0, 0);
    for (unsigned i = 0; i < N; i++)
    {
      // Get the next entry in mantissa exponent form
      std::complex<double> entry;
      int iexp_real = 0, iexp_imag = 0;
      entry =
        std::complex<double>(frexp(LU_factors[M * i + i].real(), &iexp_real),
                             frexp(LU_factors[M * i + i].imag(), &iexp_imag));

      // Now calculate the four terms that contribute to the real
      // and imaginary parts
      // i.e. A + Bi  * C + Di = AC - BD + i(AD + BC)
      double temp_mantissa[4];
      int temp_exp[4];

      // Get the first contribution to the real part, AC
      temp_mantissa[0] = determinant_mantissa.real() * entry.real();
      temp_exp[0] = determinant_exponent.real() + iexp_real;
      // Get the second contribution to the real part, DB
      temp_mantissa[1] = determinant_mantissa.imag() * entry.imag();
      temp_exp[1] = determinant_exponent.imag() + iexp_imag;

      // Get the first contribution to the imaginary part, AD
      temp_mantissa[2] = determinant_mantissa.real() * entry.imag();
      temp_exp[2] = determinant_exponent.real() + iexp_imag;
      // Get the second contribution to the imaginary part, BC
      temp_mantissa[3] = determinant_mantissa.imag() * entry.real();
      temp_exp[3] = determinant_exponent.imag() + iexp_real;

      // Align the exponents in the two terms for each sum (real or imaginary)
      // We always align up to the larger exponent
      for (unsigned i = 0; i < 3; i += 2)
      {
        // If the first exponent is smaller, move it up
        if (temp_exp[i] < temp_exp[i + 1])
        {
          // The smaller term
          int lower = temp_exp[i];
          // Loop over the difference in the exponents,
          // scaling down the mantissa each time
          for (int index = lower; index < temp_exp[i + 1]; ++index)
          {
            temp_mantissa[i] /= 2.0;
            ++temp_exp[i];
          }
        }
        // Otherwise the second exponent is smaller
        else
        {
          // The smaller term is the other
          int lower = temp_exp[i + 1];
          // Loop over the difference in the exponents,
          // scaling down the mantissa each time
          for (int index = lower; index < temp_exp[i]; ++index)
          {
            temp_mantissa[i + 1] /= 2.0;
            ++temp_exp[i + 1];
          }
        }
      }

      // Now combine the terms AC - BD
      // and Combine the terms AD + BC
      determinant_mantissa =
        std::complex<double>(temp_mantissa[0] - temp_mantissa[1],
                             temp_mantissa[2] + temp_mantissa[3]);
      // The exponents will be the same, so pick one
      determinant_exponent = std::complex<int>(temp_exp[0], temp_exp[2]);

      // Finally normalise the result
      determinant_mantissa =
        std::complex<double>(frexp(determinant_mantissa.real(), &iexp_real),
                             frexp(determinant_mantissa.imag(), &iexp_imag));

      int temp_real = determinant_exponent.real() + iexp_real;
      int temp_imag = determinant_exponent.imag() + iexp_imag;

      determinant_exponent = std::complex<int>(temp_real, temp_imag);
    }

    // Integer to store the sign of the determinant
    int sign = 0;
    // Find the sign of the determinant (left or right half-plane)
    if (std::abs(determinant_mantissa) > 0.0)
    {
      sign = 1;
    }
    if (std::abs(determinant_mantissa) < 0.0)
    {
      sign = -1;
    }

    // Multiply the sign by the signature
    sign *= signature;

    // Return the sign of the determinant
    return sign;
  }

  //============================================================================
  ///  Back substitute an LU decomposed matrix.
  //============================================================================
  void DenseComplexMatrix::lubksub(Vector<std::complex<double>>& rhs)
  {
#ifdef PARANOID
    // Check Matrix is square
    if (N != M)
    {
      throw OomphLibError(
        "This matrix is not square,  the matrix MUST be square!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
    // Check that size of rhs=nrow() and index=nrow() are correct!
    if (rhs.size() != N)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The rhs vector is not the right size. It is "
                           << rhs.size() << ", it should be " << N << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    if (Index == 0)
    {
      throw OomphLibError(
        "Index vector has not been allocated. Have you called ludecompse()\n",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
    if (Index->size() != N)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The Index vector is not the right size. It is "
                           << Index->size() << ", it should be " << N
                           << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    unsigned long ii = 0;
    for (unsigned i = 0; i < N; i++)
    {
      unsigned long ip = (*Index)[i];
      std::complex<double> sum = rhs[ip];
      rhs[ip] = rhs[i];

      if (ii != 0)
      {
        for (unsigned long j = ii - 1; j < i; j++)
        {
          sum -= LU_factors[M * i + j] * rhs[j];
        }
      }
      else if (sum != 0.0)
      {
        ii = i + 1;
      }
      rhs[i] = sum;
    }

    // Now do the back substitution
    for (long i = long(N) - 1; i >= 0; i--)
    {
      std::complex<double> sum = rhs[i];
      for (long j = i + 1; j < long(N); j++)
        sum -= LU_factors[M * i + j] * rhs[j];
      rhs[i] = sum / LU_factors[M * i + i];
    }
  }

  //============================================================================
  ///  Find the residual of Ax=b, i.e. r=b-Ax
  //============================================================================
  void DenseComplexMatrix::residual(const Vector<std::complex<double>>& x,
                                    const Vector<std::complex<double>>& rhs,
                                    Vector<std::complex<double>>& residual)
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
    // Check that size of rhs = nrow()
    if (rhs.size() != N)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The rhs vector is not the right size. It is "
                           << rhs.size() << ", it should be " << N << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // Check that the size of x is correct
    if (x.size() != N)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.size() << ", it should be " << N << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif
    // If size of residual is wrong, correct it!
    if (residual.size() != N)
    {
      residual.resize(N);
    }

    // Multiply the matrix by the vector x in residual vector
    for (unsigned long i = 0; i < N; i++)
    {
      residual[i] = rhs[i];
      for (unsigned long j = 0; j < M; j++)
      {
        residual[i] -= Matrixdata[M * i + j] * x[j];
      }
    }
  }


  //============================================================================
  ///  Multiply the matrix by the vector x: soln=Ax
  //============================================================================
  void DenseComplexMatrix::multiply(const Vector<std::complex<double>>& x,
                                    Vector<std::complex<double>>& soln)
  {
#ifdef PARANOID
    // Check to see if x.size() = ncol().
    if (x.size() != M)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.size() << ", it should be " << M << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    if (soln.size() != N)
    {
      // Resize and initialize the solution vector
      soln.resize(N);
    }

    // Multiply the matrix A, by the vector x
    for (unsigned long i = 0; i < N; i++)
    {
      soln[i] = 0.0;
      for (unsigned long j = 0; j < M; j++)
      {
        soln[i] += Matrixdata[M * i + j] * x[j];
      }
    }
  }


  //=================================================================
  /// Multiply the  transposed matrix by the vector x: soln=A^T x
  //=================================================================
  void DenseComplexMatrix::multiply_transpose(
    const Vector<std::complex<double>>& x, Vector<std::complex<double>>& soln)
  {
#ifdef PARANOID
    // Check to see x.size() = nrow()
    if (x.size() != N)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.size() << ", it should be " << N << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    if (soln.size() != M)
    {
      // Resize and initialize the solution vector
      soln.resize(M);
    }

    // Initialise the solution
    for (unsigned long i = 0; i < M; i++)
    {
      soln[i] = 0.0;
    }


    // Matrix vector product
    for (unsigned long i = 0; i < N; i++)
    {
      for (unsigned long j = 0; j < M; j++)
      {
        soln[j] += Matrixdata[N * i + j] * x[i];
      }
    }
  }

  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //===================================================================
  // Interface to SuperLU wrapper
  //===================================================================
  extern "C"
  {
    int superlu_complex(int*,
                        int*,
                        int*,
                        int*,
                        std::complex<double>*,
                        int*,
                        int*,
                        std::complex<double>*,
                        int*,
                        int*,
                        int*,
                        void*,
                        int*);
  }


  //===================================================================
  /// Perform LU decomposition. Return the sign of the determinant
  //===================================================================
  int CCComplexMatrix::ludecompose()
  {
#ifdef PARANOID
    if (N != M)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "Can only solve for quadratic matrices\n"
                           << "N, M " << N << " " << M << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // SuperLU expects compressed column format: Set flag for
    // transpose (0/1) = (false/true)
    int transpose = 0;

    // Doc (0/1) = (false/true)
    int doc = 0;
    if (Doc_stats_during_solve) doc = 1;

    // Cast to integers for stupid SuperLU
    int n_aux = (int)N;
    int nnz_aux = (int)Nnz;

    // Integer to hold the sign of the determinant
    int sign = 0;

    // Just do the lu decompose phase (i=1)
    int i = 1;
    sign = superlu_complex(&i,
                           &n_aux,
                           &nnz_aux,
                           0,
                           Value,
                           Row_index,
                           Column_start,
                           0,
                           &n_aux,
                           &transpose,
                           &doc,
                           &F_factors,
                           &Info);

    // Return the sign of the determinant
    return sign;
  }


  //===================================================================
  /// Do the backsubstitution
  //===================================================================
  void CCComplexMatrix::lubksub(Vector<std::complex<double>>& rhs)
  {
#ifdef PARANOID
    if (N != rhs.size())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "Size of RHS doesn't match matrix size\n"
                           << "N, rhs.size() " << N << " " << rhs.size()
                           << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    if (N != M)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "Can only solve for quadratic matrices\n"
                           << "N, M " << N << " " << M << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    /// RHS vector
    std::complex<double>* b = new std::complex<double>[N];

    // Copy across
    for (unsigned long i = 0; i < N; i++)
    {
      b[i] = rhs[i];
    }

    // SuperLU expects compressed column format: Set flag for
    // transpose (0/1) = (false/true)
    int transpose = 0;

    // Doc (0/1) = (false/true)
    int doc = 0;
    if (Doc_stats_during_solve) doc = 1;

    // Number of RHSs
    int nrhs = 1;


    // Cast to integers for stupid SuperLU
    int n_aux = (int)N;
    int nnz_aux = (int)Nnz;

    // SuperLU in three steps (lu decompose, back subst, cleanup)
    // Do the solve phase
    int i = 2;
    superlu_complex(&i,
                    &n_aux,
                    &nnz_aux,
                    &nrhs,
                    Value,
                    Row_index,
                    Column_start,
                    b,
                    &n_aux,
                    &transpose,
                    &doc,
                    &F_factors,
                    &Info);

    // Copy b into rhs vector
    for (unsigned long i = 0; i < N; i++)
    {
      rhs[i] = b[i];
    }

    // Cleanup
    delete[] b;
  }


  //===================================================================
  /// Cleanup storage
  //===================================================================
  void CCComplexMatrix::clean_up_memory()
  {
    // Only bother if we've done an LU solve
    if (F_factors != 0)
    {
      // SuperLU expects compressed column format: Set flag for
      // transpose (0/1) = (false/true)
      int transpose = 0;

      // Doc (0/1) = (false/true)
      int doc = 0;
      if (Doc_stats_during_solve) doc = 1;


      // Cast to integers for stupid SuperLU
      int n_aux = (int)N;
      int nnz_aux = (int)Nnz;

      // SuperLU in three steps (lu decompose, back subst, cleanup)
      // Flag to indicate which solve step to do (1, 2 or 3)
      int i = 3;
      superlu_complex(&i,
                      &n_aux,
                      &nnz_aux,
                      0,
                      Value,
                      Row_index,
                      Column_start,
                      0,
                      &n_aux,
                      &transpose,
                      &doc,
                      &F_factors,
                      &Info);
    }
  }


  //===================================================================
  /// Work out residual vector r = b-Ax for candidate solution x
  //===================================================================
  void CCComplexMatrix::residual(const Vector<std::complex<double>>& x,
                                 const Vector<std::complex<double>>& rhs,
                                 Vector<std::complex<double>>& residual)
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
    // Check that size of rhs = nrow()
    if (rhs.size() != N)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The rhs vector is not the right size. It is "
                           << rhs.size() << ", it should be " << N << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // Check that the size of x is correct
    if (x.size() != N)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.size() << ", it should be " << N << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    unsigned long r_n = residual.size();
    if (r_n != N)
    {
      residual.resize(N);
    }

    // Need to do this in loop over rows
    for (unsigned i = 0; i < N; i++)
    {
      residual[i] = rhs[i];
    }
    // Now loop over columns
    for (unsigned long j = 0; j < N; j++)
    {
      for (long k = Column_start[j]; k < Column_start[j + 1]; k++)
      {
        unsigned long i = Row_index[k];
        std::complex<double> a_ij = Value[k];
        residual[i] -= a_ij * x[j];
      }
    }
  }

  //===================================================================
  ///  Multiply the matrix by the vector x
  //===================================================================
  void CCComplexMatrix::multiply(const Vector<std::complex<double>>& x,
                                 Vector<std::complex<double>>& soln)
  {
#ifdef PARANOID
    // Check to see if x.size() = ncol()
    if (x.size() != M)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.size() << ", it should be " << M << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    if (soln.size() != N)
    {
      // Resize and initialize the solution vector
      soln.resize(N);
    }
    for (unsigned i = 0; i < N; i++)
    {
      soln[i] = 0.0;
    }

    for (unsigned long j = 0; j < N; j++)
    {
      for (long k = Column_start[j]; k < Column_start[j + 1]; k++)
      {
        unsigned long i = Row_index[k];
        std::complex<double> a_ij = Value[k];
        soln[i] += a_ij * x[j];
      }
    }
  }


  //=================================================================
  /// Multiply the  transposed matrix by the vector x: soln=A^T x
  //=================================================================
  void CCComplexMatrix::multiply_transpose(
    const Vector<std::complex<double>>& x, Vector<std::complex<double>>& soln)
  {
#ifdef PARANOID
    // Check to see x.size() = nrow()
    if (x.size() != N)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.size() << ", it should be " << N << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    if (soln.size() != M)
    {
      // Resize and initialize the solution vector
      soln.resize(M);
    }

    // Initialise the solution
    for (unsigned long i = 0; i < M; i++)
    {
      soln[i] = 0.0;
    }

    // Matrix vector product
    for (unsigned long i = 0; i < N; i++)
    {
      for (long k = Column_start[i]; k < Column_start[i + 1]; k++)
      {
        unsigned long j = Row_index[k];
        std::complex<double> a_ij = Value[k];
        soln[j] += a_ij * x[i];
      }
    }
  }


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //===================================================================
  /// Do LU decomposition and return sign of determinant
  //===================================================================
  int CRComplexMatrix::ludecompose()
  {
#ifdef PARANOID
    if (N != M)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "Can only solve for quadratic matrices\n"
                           << "N, M " << N << " " << M << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // SuperLU expects compressed column format: Set flag for
    // transpose (0/1) = (false/true)
    int transpose = 1;

    // Doc (0/1) = (false/true)
    int doc = 0;
    if (Doc_stats_during_solve) doc = 1;

    // Copies so that const-ness is maintained
    int n_aux = int(N);
    int nnz_aux = int(Nnz);

    // Integer to hold the sign of the determinant
    int sign = 0;

    // Just do the lu decompose phase (i=1)
    int i = 1;
    sign = superlu_complex(&i,
                           &n_aux,
                           &nnz_aux,
                           0,
                           Value,
                           Column_index,
                           Row_start,
                           0,
                           &n_aux,
                           &transpose,
                           &doc,
                           &F_factors,
                           &Info);
    // Return sign of determinant
    return sign;
  }


  //===================================================================
  /// Do back-substitution
  //===================================================================
  void CRComplexMatrix::lubksub(Vector<std::complex<double>>& rhs)
  {
#ifdef PARANOID
    if (N != rhs.size())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "Size of RHS doesn't match matrix size\n"
                           << "N, rhs.size() " << N << " " << rhs.size()
                           << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    if (N != M)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "Can only solve for quadratic matrices\n"
                           << "N, M " << N << " " << M << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    /// RHS vector
    std::complex<double>* b = new std::complex<double>[N];

    // Copy across
    for (unsigned long i = 0; i < N; i++)
    {
      b[i] = rhs[i];
    }

    // SuperLU expects compressed column format: Set flag for
    // transpose (0/1) = (false/true)
    int transpose = 1;

    // Doc (0/1) = (false/true)
    int doc = 0;
    if (Doc_stats_during_solve) doc = 1;

    // Number of RHSs
    int nrhs = 1;

    // Copies so that const-ness is maintained
    int n_aux = int(N);
    int nnz_aux = int(Nnz);

    // SuperLU in three steps (lu decompose, back subst, cleanup)
    // Do the solve phase
    int i = 2;
    superlu_complex(&i,
                    &n_aux,
                    &nnz_aux,
                    &nrhs,
                    Value,
                    Column_index,
                    Row_start,
                    b,
                    &n_aux,
                    &transpose,
                    &doc,
                    &F_factors,
                    &Info);

    // Copy b into rhs vector
    for (unsigned long i = 0; i < N; i++)
    {
      rhs[i] = b[i];
    }

    // Cleanup
    delete[] b;
  }


  //===================================================================
  /// Cleanup memory
  //===================================================================
  void CRComplexMatrix::clean_up_memory()
  {
    // Only bother if we've done an LU solve
    if (F_factors != 0)
    {
      // SuperLU expects compressed column format: Set flag for
      // transpose (0/1) = (false/true)
      int transpose = 1;

      // Doc (0/1) = (false/true)
      int doc = 0;
      if (Doc_stats_during_solve) doc = 1;

      // Copies so that const-ness is maintained
      int n_aux = int(N);
      int nnz_aux = int(Nnz);

      // SuperLU in three steps (lu decompose, back subst, cleanup)
      // Flag to indicate which solve step to do (1, 2 or 3)
      int i = 3;
      superlu_complex(&i,
                      &n_aux,
                      &nnz_aux,
                      0,
                      Value,
                      Column_index,
                      Row_start,
                      0,
                      &n_aux,
                      &transpose,
                      &doc,
                      &F_factors,
                      &Info);
    }
  }


  //=================================================================
  ///  Find the residulal to x of Ax=b, ie r=b-Ax
  //=================================================================
  void CRComplexMatrix::residual(const Vector<std::complex<double>>& x,
                                 const Vector<std::complex<double>>& rhs,
                                 Vector<std::complex<double>>& residual)
  {
#ifdef PARANOID
    // Check that size of rhs = nrow()
    if (rhs.size() != N)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The rhs vector is not the right size. It is "
                           << rhs.size() << ", it should be " << N << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // Check that the size of x is correct
    if (x.size() != M)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.size() << ", it should be " << M << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    if (residual.size() != N)
    {
      residual.resize(N);
    }

    for (unsigned long i = 0; i < N; i++)
    {
      residual[i] = rhs[i];
      for (long k = Row_start[i]; k < Row_start[i + 1]; k++)
      {
        unsigned long j = Column_index[k];
        std::complex<double> a_ij = Value[k];
        residual[i] -= a_ij * x[j];
      }
    }
  }

  //=================================================================
  ///  Multiply the matrix by the vector x
  //=================================================================
  void CRComplexMatrix::multiply(const Vector<std::complex<double>>& x,
                                 Vector<std::complex<double>>& soln)
  {
#ifdef PARANOID
    // Check to see x.size() = ncol()
    if (x.size() != M)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.size() << ", it should be " << M << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    if (soln.size() != N)
    {
      // Resize and initialize the solution vector
      soln.resize(N);
    }
    for (unsigned long i = 0; i < N; i++)
    {
      soln[i] = 0.0;
      for (long k = Row_start[i]; k < Row_start[i + 1]; k++)
      {
        unsigned long j = Column_index[k];
        std::complex<double> a_ij = Value[k];
        soln[i] += a_ij * x[j];
      }
    }
  }


  //=================================================================
  /// Multiply the  transposed matrix by the vector x: soln=A^T x
  //=================================================================
  void CRComplexMatrix::multiply_transpose(
    const Vector<std::complex<double>>& x, Vector<std::complex<double>>& soln)
  {
#ifdef PARANOID
    // Check to see x.size() = nrow()
    if (x.size() != N)
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector is not the right size. It is "
                           << x.size() << ", it should be " << N << std::endl;

      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    if (soln.size() != M)
    {
      // Resize and initialize the solution vector
      soln.resize(M);
    }

    // Initialise the solution
    for (unsigned long i = 0; i < M; i++)
    {
      soln[i] = 0.0;
    }

    // Matrix vector product
    for (unsigned long i = 0; i < N; i++)
    {
      for (long k = Row_start[i]; k < Row_start[i + 1]; k++)
      {
        unsigned long j = Column_index[k];
        std::complex<double> a_ij = Value[k];
        soln[j] += a_ij * x[i];
      }
    }
  }
} // namespace oomph
