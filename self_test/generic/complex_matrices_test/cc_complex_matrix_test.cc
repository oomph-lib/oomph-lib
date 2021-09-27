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

// Oomph-lib includes
#include "generic.h"

using namespace std;
using namespace oomph;

void print_complex_vector(Vector<complex<double>>& x)
{
  unsigned vector_length = x.size();
  for (unsigned i = 0; i < vector_length; i++)
  {
    cout << x[i] << endl;
  }
}

void print_cc_complex_matrix(CCComplexMatrix& matrix)
{
  unsigned n_row = matrix.nrow();
  unsigned n_col = matrix.ncol();
  for (unsigned i = 0; i < n_row; i++)
  {
    for (unsigned j = 0; j < n_col; j++)
    {
      cout << matrix(i, j) << ", ";
    }
    cout << endl;
  }
}

int main()
{
  Vector<complex<double>> values(4);
  values[0] = complex<double>(1, 1);
  values[1] = complex<double>(2, 1);
  values[2] = complex<double>(3, 1);
  values[3] = complex<double>(4, 1);

  Vector<int> row_index(4);
  row_index[0] = 0;
  row_index[1] = 2;
  row_index[2] = 1;
  row_index[3] = 3;

  // Final entry stores a fictitous index equal to the number of non-zeros, NNz
  constexpr unsigned n_column = 2;
  Vector<int> column_start(n_column + 1);
  column_start[0] = 0;
  column_start[1] = 2;
  column_start[2] = 4;

  Vector<complex<double>> values_square(7);
  values_square[0] = complex<double>(1, 1);
  values_square[1] = complex<double>(2, -1);
  values_square[2] = complex<double>(3, 0.1);
  values_square[3] = complex<double>(4, 3);
  values_square[4] = complex<double>(5, 1);
  values_square[5] = complex<double>(6, -2);
  values_square[6] = complex<double>(0.5, -1);

  Vector<int> row_index_square(7);
  row_index_square[0] = 0;
  row_index_square[1] = 2;
  row_index_square[2] = 0;
  row_index_square[3] = 1;
  row_index_square[4] = 1;
  row_index_square[5] = 2;
  row_index_square[6] = 3;

  // Final entry stores a fictitous index equal to the number of non-zeros, NNz
  constexpr unsigned n_column_square = 4;
  Vector<int> column_start_square(n_column_square + 1);
  column_start_square[0] = 0;
  column_start_square[1] = 2;
  column_start_square[2] = 4;
  column_start_square[3] = 6;
  column_start_square[4] = 7;

  Vector<std::complex<double>> x(4);
  x[0] = complex<double>(2, 2);
  x[1] = complex<double>(-2, 3);
  x[2] = complex<double>(1, -2);
  x[3] = complex<double>(-2, 1);

  Vector<std::complex<double>> rhs(4);
  rhs[0] = complex<double>(1, -0.5);
  rhs[1] = complex<double>(0.5, 1.2);
  rhs[2] = complex<double>(2, 0.3);
  rhs[3] = complex<double>(-3, -0.4);

  Vector<std::complex<double>> soln(4);

  // test default constructor
  CCComplexMatrix matrix_default;
  cout << matrix_default.nrow() << endl;
  cout << matrix_default.ncol() << endl;

  // test full matrix constructor
  constexpr unsigned long n = 2;
  constexpr unsigned long m = 4;
  CCComplexMatrix matrix(values, row_index, column_start, n, m);
  cout << matrix.nrow() << endl;
  cout << matrix.ncol() << endl;
  print_cc_complex_matrix(matrix);

  CCComplexMatrix matrix_square(
    values_square, row_index_square, column_start_square, 4, 4);
  cout << matrix_square.nrow() << endl;
  cout << matrix_square.ncol() << endl;
  print_cc_complex_matrix(matrix_square);

  // test LU decomposition
  cout << matrix_square.ludecompose() << endl;

  // test residual
  matrix_square.lubksub(rhs);
  print_complex_vector(rhs);

  // test multiply
  matrix_square.multiply(x, soln);
  print_complex_vector(soln);

  // test multiply transposed
  matrix_square.multiply_transpose(x, soln);
  print_complex_vector(soln);

  return (EXIT_SUCCESS);
}
