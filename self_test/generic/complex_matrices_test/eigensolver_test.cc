//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================

// Oomph-lib includes
#include "generic.h"

using namespace std;
using namespace oomph;

int main()
{
  Vector<int> row_start(3);
  row_start[0] = 0;
  row_start[1] = 2;
  row_start[2] = 4;

  Vector<int> column_index(4);
  column_index[0] = 0;
  column_index[1] = 1;
  column_index[2] = 0;
  column_index[3] = 1;

  Vector<complex<double>> values(4);
  values[0] = complex<double>(1, 1);
  values[1] = complex<double>(0, 1);
  values[2] = complex<double>(0, 1);
  values[3] = complex<double>(1, 1);

  Vector<complex<double>> rhs(2);
  rhs[0] = complex<double>(1, 0);
  rhs[1] = complex<double>(2, 0);

  CRComplexMatrix A(values, column_index, row_start, 2, 2);

  values[0] = complex<double>(1, 0);
  values[1] = complex<double>(0, 1);
  values[2] = complex<double>(0.2, 0);
  values[3] = complex<double>(1, 0);

  CCComplexMatrix M(values, column_index, row_start, 2, 2);

  LAPACK_QZ eigen_solver;

  Vector<complex<double>> eval;
  Vector<Vector<complex<double>>> evec;

  eigen_solver.find_eigenvalues(A, M, eval, evec);

  // test eigen_solver complex output
  for (int i = 0; i < 2; i++)
  {
    cout << eval[i] << ": ";
    for (int j = 0; j < 2; j++)
    {
      cout << evec[i][j] << ", ";
    }
    cout << endl;
  }

  return(EXIT_SUCCESS);
}
