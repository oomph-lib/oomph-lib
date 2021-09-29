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

void print_complex_vector(const Vector<complex<double>>& x)
{
  unsigned vector_length = x.size();
  for (unsigned i = 0; i < vector_length; i++)
  {
    cout << x[i] << endl;
  }
}

void print_cr_complex_matrix(const CRComplexMatrix& matrix)
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
  // test default constructor
  CRComplexMatrix matrix_default;
  cout << matrix_default.nrow() << endl;
  cout << matrix_default.ncol() << endl;

  constexpr unsigned long n_col_rect = 2;
  constexpr unsigned long n_row_rect = 4;
  // Number of non-zero entries
  constexpr unsigned n_nz_rect = 4;

  Vector<complex<double>> value_rect(n_nz_rect);
  value_rect[0] = complex<double>(1.0, 1.0);
  value_rect[1] = complex<double>(2.0, 1.0);
  value_rect[2] = complex<double>(3.0, 1.0);
  value_rect[3] = complex<double>(4.0, 1.0);

  Vector<int> column_index_rect(n_nz_rect);
  column_index_rect[0] = 0;
  column_index_rect[1] = 2;
  column_index_rect[2] = 1;
  column_index_rect[3] = 3;

  // Final entry stores a fictitous index equal to the number of non-zeros
  Vector<int> row_start_rect(n_col_rect + 1);
  row_start_rect[0] = 0;
  row_start_rect[1] = 2;
  row_start_rect[2] = n_nz_rect;

  // test full matrix constructor with a rectangular matrix
  CRComplexMatrix matrix_rect(
    value_rect, column_index_rect, row_start_rect, n_col_rect, n_row_rect);
  cout << matrix_rect.nrow() << endl;
  cout << matrix_rect.ncol() << endl;
  print_cr_complex_matrix(matrix_rect);

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

  Vector<int> column_index_square(n_nz_square);
  column_index_square[0] = 0;
  column_index_square[1] = 2;
  column_index_square[2] = 0;
  column_index_square[3] = 1;
  column_index_square[4] = 1;
  column_index_square[5] = 2;
  column_index_square[6] = 3;

  // Final entry stores a fictitous index equal to the number of non-zeros
  Vector<int> row_start_square(n_col_square + 1);
  row_start_square[0] = 0;
  row_start_square[1] = 2;
  row_start_square[2] = 4;
  row_start_square[3] = 6;
  row_start_square[4] = n_nz_square;

  // test full matrix constructor with a square matrix
  CRComplexMatrix matrix_square(values_square,
                                column_index_square,
                                row_start_square,
                                n_col_square,
                                n_row_square);
  cout << matrix_square.nrow() << endl;
  cout << matrix_square.ncol() << endl;
  print_cr_complex_matrix(matrix_square);

  // test LU decomposition
  cout << matrix_square.ludecompose() << endl;

  constexpr unsigned long vector_length = n_row_square;

  Vector<std::complex<double>> rhs(vector_length);
  rhs[0] = complex<double>(1.0, 0.0);
  rhs[1] = complex<double>(0.0, 1.0);
  rhs[2] = complex<double>(2.0, 1.0);
  rhs[3] = complex<double>(-3.0, 1.0);

  // test residual
  matrix_square.lubksub(rhs);
  print_complex_vector(rhs);

  Vector<std::complex<double>> x(vector_length);
  x[0] = complex<double>(2.0, 2.0);
  x[1] = complex<double>(-2.0, 3.0);
  x[2] = complex<double>(0.5, 2.2);
  x[3] = complex<double>(-1.4, 3.0);

  Vector<std::complex<double>> soln(vector_length);

  // test multiply
  matrix_square.multiply(x, soln);
  print_complex_vector(soln);

  // test multiply transposed
  matrix_square.multiply_transpose(x, soln);
  print_complex_vector(soln);

  LinearAlgebraDistribution* dist_pt = new LinearAlgebraDistribution();

  constexpr unsigned n_nz_double = 5;

  Vector<double> value_double(5);
  value_double[0] = double(1.0);
  value_double[1] = double(2.0);
  value_double[2] = double(3.0);
  value_double[3] = double(4.0);
  value_double[4] = double(5.0);

  Vector<int> column_index_double(n_nz_double);
  column_index_double[0] = 0;
  column_index_double[1] = 2;
  column_index_double[2] = 0;
  column_index_double[3] = 1;
  column_index_double[4] = 1;

  Vector<int> row_start_double(n_col_square + 1);
  row_start_double[0] = 0;
  row_start_double[1] = 2;
  row_start_double[2] = 3;
  row_start_double[3] = 4;
  row_start_double[4] = n_nz_double;

  CRDoubleMatrix matrix_in_double(
    dist_pt, n_col_square, value_double, column_index_double, row_start_double);

  // test add (CRDoubleMatrix)
  matrix_square.add(matrix_in_double, matrix_square);
  print_cr_complex_matrix(matrix_square);

  // test add (CRComplexMatrix)
  matrix_square.add(matrix_square, matrix_square);
  print_cr_complex_matrix(matrix_square);

  // Print input matrices
  cout << "A" << endl;
  print_cr_complex_matrix(matrix_square);

  CRComplexMatrix matrix_square_2(values_square,
                                column_index_square,
                                row_start_square,
                                n_col_square,
                                n_row_square);

  cout << "B" << endl;
  print_cr_complex_matrix(matrix_square_2);

  CRComplexMatrix matrix_result(values_square,
                                column_index_square,
                                row_start_square,
                                n_col_square,
                                n_row_square);

  // test multiply method 2 (default)
  matrix_square.multiply(matrix_square_2, matrix_result);
  cout << "A B" << endl;
  print_cr_complex_matrix(matrix_result);

  // Create a reference to Serial_matrix_matrix_multiply_method
  unsigned& method = matrix_square.serial_matrix_matrix_multiply_method();
  cout << method << endl;

  // test multiply method 1
  method = 1;
  matrix_square.multiply(matrix_square_2, matrix_result);
  cout << method << endl;
  cout << "A B" << endl;
  print_cr_complex_matrix(matrix_result);

  // test multiply method 3
  method = 3;
  matrix_square.multiply(matrix_square_2, matrix_result);
  cout << method << endl;
  cout << "A B" << endl;
  print_cr_complex_matrix(matrix_result);

  // Print cleanup messages for easier debugging
  cout << "Begin cleanup" << endl;
  delete dist_pt;
  dist_pt = nullptr;
  cout << "Finished cleanup" << endl;

  return (EXIT_SUCCESS);
}
