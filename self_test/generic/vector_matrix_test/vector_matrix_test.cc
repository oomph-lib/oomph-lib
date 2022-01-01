//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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

using namespace oomph;

//===start_of_main======================================================
/// Driver code: Testing VectorMatrix class
//======================================================================
int main(int argc, char* argv[])
{
  std::ostringstream out_stream;
  out_stream << "OUTPUT";

  std::ofstream out_file;
  out_file.open(out_stream.str().c_str());

  VectorMatrix<int> vector_matrix_1(3,4,-1);

  // test 1: nrow(), should be 3
  out_file << vector_matrix_1.nrow() << "\n";

  // test 2: ncol(), should be 4
  out_file << vector_matrix_1.ncol() << "\n";

  // test 3: [] operator, should be -1
  out_file << vector_matrix_1[2][2] << "\n";

  // test 4: [] operator
  vector_matrix_1[2][2] = 42;
  // should be 42
  out_file << vector_matrix_1[2][2] << "\n";

  // test 5: resize
  vector_matrix_1.resize(6,8,-13);

  // should be 6
  out_file << vector_matrix_1.nrow() << "\n";
  // should be 8
  out_file << vector_matrix_1.ncol() << "\n";

  // elements [0][0] and [2][2] is not changed since it is within the 
  // resize bound.
  
  // should be -1
  out_file << vector_matrix_1[0][0] << "\n";
  // should be 42
  out_file << vector_matrix_1[2][2] << "\n";
  // should be -13 (the value given to resize)
  out_file << vector_matrix_1[5][7] << "\n";

  // test 6: assign
  vector_matrix_1.assign(24,32,-42);


  // All elements should be destroyed and replaced, so all elements should be
  // -42. Check a random one.
  out_file << vector_matrix_1[24 * (rand() / (RAND_MAX + 1.0))]
                             [32 * (rand() / (RAND_MAX + 1.0))]<<"\n";

  // test clear
  vector_matrix_1.clear();
  // should be 0
  out_file << vector_matrix_1.nrow() << "\n";
  // should be 0
  out_file << vector_matrix_1.ncol() << "\n";

  /// /////////////////////////////////////////////////////////
  
  // Lastly, check other constructors.
  const unsigned nrow = 101;
  const unsigned ncol = 97;
  const double val = 77.7;

  // Default constructor
  VectorMatrix<double> vector_matrix_2;
  // should be 0
  out_file << vector_matrix_2.nrow() << "\n";
  // should be 0
  out_file << vector_matrix_2.ncol() << "\n";

  // constructor: n by m, given val
  VectorMatrix<double> vector_matrix_3(nrow,ncol,val);
  // should be 101
  out_file << vector_matrix_3.nrow() << "\n";
  // should be 97
  out_file << vector_matrix_3.ncol() << "\n";
  // should be 77.7
  out_file << vector_matrix_3[nrow * (rand() / (RAND_MAX + 1.0))]
                             [ncol * (rand() / (RAND_MAX + 1.0))]<<"\n";

  // constructor: n by n, given val
  VectorMatrix<double>vector_matrix_4(nrow,ncol);
  // should be 101
  out_file << vector_matrix_4.nrow() << "\n";
  // should be 97
  out_file << vector_matrix_4.ncol() << "\n";
  // should be 0.0
  out_file << vector_matrix_4[nrow * (rand() / (RAND_MAX + 1.0))]
                             [nrow * (rand() / (RAND_MAX + 1.0))]<<"\n";

  out_file.close();

  return(EXIT_SUCCESS);
} // end_of_main
