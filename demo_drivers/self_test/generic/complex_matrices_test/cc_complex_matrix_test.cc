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

/// Main test function for the compressed-column complex matrix class
int main()
{
  // test default constructor
  CCComplexMatrix matrix_default;
  oomph_info << matrix_default.nrow() << endl;
  oomph_info << matrix_default.ncol() << endl;

  constexpr unsigned long n_col_rect = 4;
  constexpr unsigned long n_row_rect = 2;
  // Number of non-zero entries
  constexpr unsigned n_nz_rect = 4;

  Vector<complex<double>> value_rect(n_nz_rect);
  value_rect[0] = complex<double>(1.0, 1.0);
  value_rect[1] = complex<double>(2.0, 1.0);
  value_rect[2] = complex<double>(3.0, 1.0);
  value_rect[3] = complex<double>(4.0, 1.0);

  Vector<int> row_index_rect(n_nz_rect);
  row_index_rect[0] = 0;
  row_index_rect[1] = 2;
  row_index_rect[2] = 1;
  row_index_rect[3] = 3;

  // Final entry stores a fictitous index equal to the number of non-zeros
  Vector<int> column_start_rect(n_col_rect + 1);
  column_start_rect[0] = 0;
  column_start_rect[1] = 2;
  column_start_rect[2] = n_nz_rect;

  // test full matrix constructor with a rectangular matrix
  CCComplexMatrix matrix_rect(
    value_rect, row_index_rect, column_start_rect, n_col_rect, n_row_rect);
  oomph_info << matrix_rect.nrow() << endl;
  oomph_info << matrix_rect.ncol() << endl;
  print_complex_matrix(matrix_rect);

  constexpr unsigned long n_col_square = 4;
  constexpr unsigned long n_row_square = 4;
  constexpr unsigned n_nz_square = 7;

  Vector<complex<double>> values_square(n_nz_square);
  values_square[0] = complex<double>(1.0, 1.0);
  values_square[1] = complex<double>(2.0, -1.0);
  values_square[2] = complex<double>(3.0, 0.1);
  values_square[3] = complex<double>(4.0, 3.0);
  values_square[4] = complex<double>(5.0, 1.0);
  values_square[5] = complex<double>(6.0, -2.0);
  values_square[6] = complex<double>(0.5, -1.0);

  Vector<int> row_index_square(n_nz_square);
  row_index_square[0] = 0;
  row_index_square[1] = 2;
  row_index_square[2] = 0;
  row_index_square[3] = 1;
  row_index_square[4] = 1;
  row_index_square[5] = 2;
  row_index_square[6] = 3;

  // Final entry stores a fictitous index equal to the number of non-zeros
  Vector<int> column_start_square(n_col_square + 1);
  column_start_square[0] = 0;
  column_start_square[1] = 2;
  column_start_square[2] = 4;
  column_start_square[3] = 6;
  column_start_square[4] = n_nz_square;

  // test full matrix constructor with a square matrix
  CCComplexMatrix matrix_square(values_square,
                                row_index_square,
                                column_start_square,
                                n_col_square,
                                n_row_square);
  oomph_info << matrix_square.nrow() << endl;
  oomph_info << matrix_square.ncol() << endl;
  print_complex_matrix(matrix_square);

  // test LU decomposition
  oomph_info << matrix_square.ludecompose() << endl;

  constexpr unsigned long vector_length = n_row_square;

  Vector<complex<double>> rhs(vector_length);
  rhs[0] = complex<double>(1.0, -0.5);
  rhs[1] = complex<double>(0.5, 1.2);
  rhs[2] = complex<double>(2.0, 0.3);
  rhs[3] = complex<double>(-3.0, -0.4);

  // test residual
  matrix_square.lubksub(rhs);
  print_complex_vector(rhs);

  Vector<complex<double>> x(vector_length);
  x[0] = complex<double>(2.0, 2.0);
  x[1] = complex<double>(-2.0, 3.0);
  x[2] = complex<double>(1.0, -2.0);
  x[3] = complex<double>(-2.0, 1.0);

  Vector<complex<double>> soln(vector_length);

  // test multiply
  matrix_square.multiply(x, soln);
  print_complex_vector(soln);

  // test multiply transposed
  matrix_square.multiply_transpose(x, soln);
  print_complex_vector(soln);

  return (EXIT_SUCCESS);
}
