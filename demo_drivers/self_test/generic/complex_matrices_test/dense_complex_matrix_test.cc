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

// Oomph-lib includes
#include "generic.h"
#include "complex_matrices_test_utility_functions.h"

using namespace std;
using namespace oomph;

/// Main test function for the dense complex matrix class
int main()
{
  // test default constructor
  DenseComplexMatrix matrix_default;
  oomph_info << matrix_default.nrow() << endl;
  oomph_info << matrix_default.ncol() << endl;

  constexpr unsigned long n_row = 2;
  constexpr unsigned long n_col = 3;
  // Number of non-zero entries
  constexpr unsigned n_nz = 4;

  // test square matrix constructor
  DenseComplexMatrix matrix_square(n_row);
  oomph_info << matrix_square.nrow() << endl;
  oomph_info << matrix_square.ncol() << endl;


  // test rectangular matrix constructor
  DenseComplexMatrix matrix_rect(n_row, n_col);
  oomph_info << matrix_rect.nrow() << endl;
  oomph_info << matrix_rect.ncol() << endl;

  complex<double> inital_value(0.0, 1.0);

  // test constructor with every entry equal to a single initial value
  DenseComplexMatrix matrix(n_row, n_col, inital_value);
  oomph_info << matrix.nrow() << endl;
  oomph_info << matrix.ncol() << endl;
  print_complex_matrix(matrix);

  Vector<complex<double>> values(n_nz);
  values[0] = complex<double>(1.0, 1.0);
  values[1] = complex<double>(0.0, 1.0);
  values[2] = complex<double>(2.0, 1.0);
  values[3] = complex<double>(1.0, 1.0);

  // test round bracket access overloaded operator
  matrix_square(0, 0) = values[0];
  matrix_square(0, 1) = values[1];
  matrix_square(1, 0) = values[2];
  matrix_square(1, 1) = values[3];
  print_complex_matrix(matrix_square);

  // test LU decomposition
  oomph_info << matrix_square.ludecompose() << endl;

  constexpr unsigned long vector_length = n_row;

  Vector<complex<double>> rhs(vector_length);
  rhs[0] = complex<double>(1.0, 0.0);
  rhs[1] = complex<double>(1.4, 1.0);

  // test LU back substitution
  // note: overwrites rhs
  matrix_square.lubksub(rhs);
  print_complex_vector(rhs);

  Vector<complex<double>> x(vector_length);
  x[0] = complex<double>(1.0, 1.5);
  x[1] = complex<double>(-2.0, 3.0);

  Vector<complex<double>> residual(vector_length);

  // test residual
  matrix_square.residual(x, rhs, residual);
  print_complex_vector(residual);

  Vector<complex<double>> soln(vector_length);

  // test multiply
  matrix_square.multiply(x, soln);
  print_complex_vector(soln);

  // test multiply transposed
  matrix_square.multiply_transpose(x, soln);
  print_complex_vector(soln);

  return (EXIT_SUCCESS);
}
