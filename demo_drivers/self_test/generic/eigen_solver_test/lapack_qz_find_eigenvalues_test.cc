// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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

/// Main test function for the eigensolver of the compressed form of a complex
/// matrix eigenvalue problem
int main()
{
  // Create DocInfo and add an output directory
  DocInfo doc_info;
  doc_info.set_directory("RESLT/");

  // Create eigensolver
  LAPACK_QZ eigen_solver;

  constexpr unsigned long n_row = 2;
  constexpr unsigned long n_col = 2;
  // Number of non zero entries for the different matrices
  constexpr unsigned long n_nz_cr = 4;
  constexpr unsigned long n_nz_cc = 4;

  Vector<int> row_start_cr(n_row + 1);
  row_start_cr[0] = 0;
  row_start_cr[1] = 2;
  row_start_cr[2] = n_nz_cr;

  Vector<int> column_index_cr(n_nz_cr);
  column_index_cr[0] = 0;
  column_index_cr[1] = 1;
  column_index_cr[2] = 0;
  column_index_cr[3] = 1;

  Vector<complex<double>> values_cr(n_nz_cr);
  values_cr[0] = complex<double>(1.0, 1.0);
  values_cr[1] = complex<double>(0.0, 1.0);
  values_cr[2] = complex<double>(0.0, 1.0);
  values_cr[3] = complex<double>(1.0, 1.0);

  CRComplexMatrix matrix_cr(
    values_cr, column_index_cr, row_start_cr, n_row, n_col);

  Vector<int> row_start_cc(n_row + 1);
  row_start_cc[0] = 0;
  row_start_cc[1] = 2;
  row_start_cc[2] = n_nz_cc;

  Vector<int> column_index_cc(n_nz_cc);
  column_index_cc[0] = 0;
  column_index_cc[1] = 1;
  column_index_cc[2] = 0;
  column_index_cc[3] = 1;

  Vector<complex<double>> values_cc(n_nz_cc);
  values_cc[0] = complex<double>(1.0, 0.0);
  values_cc[1] = complex<double>(0.0, 1.0);
  values_cc[2] = complex<double>(0.2, 0.0);
  values_cc[3] = complex<double>(1.0, 0.0);

  CCComplexMatrix matrix_cc(
    values_cc, column_index_cc, row_start_cc, n_row, n_col);

  Vector<complex<double>> eval;
  Vector<Vector<complex<double>>> evec;

  // Test eigen_solver with complex matrices
  eigen_solver.find_eigenvalues(matrix_cr, matrix_cc, eval, evec);

  ofstream output_stream;
  output_stream.open(doc_info.directory() +
                     "lapack_qz_find_eigenvalues_test.dat");

  // Output real and imaginary parts of the eigenvalues and eigenvectors
  for (int i = 0; i < 2; i++)
  {
    output_stream << eval[i].real() << " , " << eval[i].imag() << endl << endl;
    for (int j = 0; j < 2; j++)
    {
      output_stream << evec[i][j].real() << " , " << evec[i][j].imag() << endl;
    }
    output_stream << endl;
  }
  output_stream.close();

  return (EXIT_SUCCESS);
}
