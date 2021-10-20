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
#include "complex_matrices_test_utility_functions.h"

using namespace std;
using namespace oomph;

/// Print complex vector to console out
void print_complex_vector(const Vector<complex<double>>& x)
{
  unsigned vector_length = x.size();
  for (unsigned i = 0; i < vector_length; i++)
  {
    cout << x[i] << endl;
  }
}

/// Print a complex matrix in its dense form to console out
void print_complex_matrix(const ComplexMatrixBase& matrix)
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
