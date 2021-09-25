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
  for (unsigned i = 0; i < x.size(); i++)
  {
    cout << x[i] << endl;
  }
}

void print_cr_complex_matrix(CRComplexMatrix& matrix)
{
  for (unsigned i = 0; i < matrix.nrow(); i++)
  {
    for (unsigned j = 0; j < matrix.ncol(); j++)
    {
      cout << matrix(i, j) << ", ";
    }
    cout << endl;
  }
}

int main()
{
  unsigned long n = 2;
  unsigned long m = 4;

  Vector<complex<double>> values(4);
  values[0] = complex<double>(1, 1);
  values[1] = complex<double>(2, 1);
  values[2] = complex<double>(3, 1);
  values[3] = complex<double>(4, 1);

  Vector<int> column_index(4);
  column_index[0] = 0;
  column_index[1] = 2;
  column_index[2] = 1;
  column_index[3] = 3;

  Vector<int> row_start(3); // length = n+1
  row_start[0] = 0;
  row_start[1] = 2;
  row_start[2] = 4; // Fictitous index = NNz

  Vector<complex<double>> values_square(5);
  values_square[0] = complex<double>(1, 1);
  values_square[1] = complex<double>(2, 1);
  values_square[2] = complex<double>(3, 1);
  values_square[3] = complex<double>(4, 1);
  values_square[4] = complex<double>(5, 1);

  Vector<int> column_index_square(5);
  column_index_square[0] = 0;
  column_index_square[1] = 2;
  column_index_square[2] = 1;
  column_index_square[3] = 0;
  column_index_square[4] = 3;

  Vector<int> row_start_square(5); // length = n+1
  row_start_square[0] = 0;
  row_start_square[1] = 2;
  row_start_square[2] = 3;
  row_start_square[3] = 4;
  row_start_square[4] = 5; // Fictitous index = NNz

  Vector<std::complex<double>> x(4);
  x[0] = complex<double>(2, 2);
  x[1] = complex<double>(-2, 3);
  x[2] = complex<double>(0.5, 2.2);
  x[3] = complex<double>(-1.4, 3);

  Vector<std::complex<double>> rhs(4);
  rhs[0] = complex<double>(1, 0);
  rhs[1] = complex<double>(0, 1);
  rhs[2] = complex<double>(2, 1);
  rhs[3] = complex<double>(-3, 1);

  Vector<std::complex<double>> soln(4);

  Vector<double> value_double(5);
  value_double[0] = double(1);
  value_double[1] = double(2);
  value_double[2] = double(3);
  value_double[3] = double(4);
  value_double[4] = double(5);

  unsigned* method_ptr = new unsigned;

  // test default constructor
  CRComplexMatrix matrix_default;
  cout << matrix_default.nrow() << endl;
  cout << matrix_default.ncol() << endl;

  // test full matrix constructor
  CRComplexMatrix matrix(values, column_index, row_start, n, m);
  cout << matrix.nrow() << endl;
  cout << matrix.ncol() << endl;
  print_cr_complex_matrix(matrix);

  CRComplexMatrix matrix_square(
    values_square, column_index_square, row_start_square, 4, 4);
  cout << matrix_square.nrow() << endl;
  cout << matrix_square.ncol() << endl;
  print_cr_complex_matrix(matrix_square);

  CRComplexMatrix matrix_result(
    values_square, column_index_square, row_start_square, 4, 4);
  print_cr_complex_matrix(matrix_result);

  unsigned ncol_double = matrix_square.ncol();
  LinearAlgebraDistribution* dist_ptr = new LinearAlgebraDistribution();
  CRDoubleMatrix matrix_in_double(
    dist_ptr, ncol_double, value_double, column_index, row_start);

  CRComplexMatrix matrix_square_2(
    values_square, column_index_square, row_start_square, 4, 4);

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


  // test add (CRDoubleMatrix)
  matrix_square.add(matrix_in_double, matrix_square);
  print_cr_complex_matrix(matrix_square);

  // test add (CRComplexMatrix)
  matrix_square.add(matrix_square, matrix_square);
  print_cr_complex_matrix(matrix_square);

  // Print input matrices
  cout << "A" << endl;
  print_cr_complex_matrix(matrix_square);
  cout << "B" << endl;
  print_cr_complex_matrix(matrix_square_2);

  // test multiply method 2 (default)
  matrix_square.multiply(matrix_square_2, matrix_result);
  cout << "A B" << endl;
  print_cr_complex_matrix(matrix_result);

  // test multiply set method and default value
  method_ptr = &matrix_square.serial_matrix_matrix_multiply_method();
  cout << *method_ptr << endl;

  // test multiply method 1
  *method_ptr = 1;
  matrix_square.multiply(matrix_square_2, matrix_result);
  cout << *method_ptr << endl;
  cout << "A B" << endl;
  print_cr_complex_matrix(matrix_result);

  // test multiply method 3
  *method_ptr = 3;
  matrix_square.multiply(matrix_square_2, matrix_result);
  cout << *method_ptr << endl;
  cout << "A B" << endl;
  print_cr_complex_matrix(matrix_result);

  // Print cleanup messages for easier debugging
  cout << "Begin cleanup" << endl;
  delete dist_ptr;
  dist_ptr = 0;
  method_ptr = 0;
  delete method_ptr;
  method_ptr = 0;
  cout << "Finished cleanup" << endl;

  return (EXIT_SUCCESS);
}
