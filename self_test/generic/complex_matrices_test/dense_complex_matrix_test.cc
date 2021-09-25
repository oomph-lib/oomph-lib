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

void print_dense_complex_matrix(DenseComplexMatrix& matrix)
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
  Vector<complex<double>> values(4);
  values[0] = complex<double>(1, 1);
  values[1] = complex<double>(0, 1);
  values[2] = complex<double>(2, 1);
  values[3] = complex<double>(1, 1);

  complex<double> single_value(0, 1);

  Vector<std::complex<double>> x(2);
  x[0] = complex<double>(1, 1.5);
  x[1] = complex<double>(-2, 3);

  Vector<std::complex<double>> rhs(2);
  rhs[0] = complex<double>(1, 0);
  rhs[1] = complex<double>(1.4, 1);

  Vector<std::complex<double>> residual(2);
  Vector<std::complex<double>> soln(2);

  // test default constructor
  DenseComplexMatrix matrix_default;
  cout << matrix_default.nrow() << endl;
  cout << matrix_default.ncol() << endl;

  // test square matrix constructor
  DenseComplexMatrix matrix_square(2);
  cout << matrix_square.nrow() << endl;
  cout << matrix_square.ncol() << endl;

  // test rectangular matrix constructor
  DenseComplexMatrix matrix_rect(3, 2);
  cout << matrix_rect.nrow() << endl;
  cout << matrix_rect.ncol() << endl;

  // test full constructor
  DenseComplexMatrix matrix(2, 2, single_value);
  cout << matrix.nrow() << endl;
  cout << matrix.ncol() << endl;
  print_dense_complex_matrix(matrix);

  // test round bracket access overloaded operator
  matrix_square(0, 0) = values[0];
  matrix_square(0, 1) = values[1];
  matrix_square(1, 0) = values[2];
  matrix_square(1, 1) = values[3];
  print_dense_complex_matrix(matrix_square);

  // test LU decomposition
  cout << matrix_square.ludecompose() << endl;

  // test LU back substitution
  // note: overwrites rhs
  matrix_square.lubksub(rhs);
  print_complex_vector(rhs);

  // test residual
  matrix_square.residual(x, rhs, residual);
  print_complex_vector(residual);

  // test multiply
  matrix_square.multiply(x, soln);
  print_complex_vector(soln);

  // test multiply transposed
  matrix_square.multiply_transpose(x, soln);
  print_complex_vector(soln);

  return (EXIT_SUCCESS);
}
