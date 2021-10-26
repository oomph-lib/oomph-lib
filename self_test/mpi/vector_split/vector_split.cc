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

using namespace oomph;

// Helper function to construct a Vector given an array
// in C++11, we can inistialise Vectors properly...
template<typename myType>
void construct_vector(myType given_array[],unsigned given_arraysize,
                      Vector<myType> &result_vector)
{
  // Clear and reserve the required memory.
  result_vector.clear();
  result_vector.reserve(given_arraysize);

  for (unsigned i = 0; i < given_arraysize; i++) 
  {
    result_vector.push_back(given_array[i]);
  }
}

// Helper function to output a Vector.
template<typename myType>
void output_vector(Vector<myType> &given_vector)
{
  typename Vector<myType>::iterator it;

  for(it = given_vector.begin(); it != given_vector.end(); ++it)
  {
    oomph_info << *it << std::endl; 
  }
}

/// ////////////////////////////////////////////////////////////////////////////

// Given the dimensions of a vector, returns a vector where the elements
// are increasing along (or rather, down), the rows.
void create_vector_ascend_row(unsigned const nrow, 
                              const OomphCommunicator* const comm_pt,
                              bool const distributed, 
                              DoubleVector &my_vec)
{  
  // Clear the block
  my_vec.clear();
  
  // Create the distribution.
  LinearAlgebraDistribution distri(comm_pt,nrow,distributed);

  // Build the vector
  my_vec.build(distri,0.0);

  // The number of rows this processor is responsible for.
  unsigned nrow_local = distri.nrow_local();
  
  // The first_row will be used as an offset for the values to insert.
  unsigned first_row = distri.first_row();

  // Fill in values...
  for (unsigned row_i = 0; row_i < nrow_local; row_i++)
  {
    my_vec[row_i] = first_row + row_i + 1; // Use natural numbers
  }
}

// Helper function to create DoubleVectors.
void create_vectors_to_split
(unsigned nrowarray[], const OomphCommunicator* const comm_pt, 
 const bool distributed, Vector<DoubleVector>&out_vector)
{
  // The number of out vectors.
  unsigned nvectors = out_vector.size();

  // Build the each sub DoubleVector with just a uniform distribution.
  // By default, each entry is initialised to 0.
  for (unsigned vec_i = 0; vec_i < nvectors; vec_i++) 
  {
    // Create the distribution
    LinearAlgebraDistribution dist(comm_pt,nrowarray[vec_i],distributed);
    out_vector[vec_i].build(dist);
  }
}

//===start_of_main======================================================
/// Driver code: Testing DoubleVectorHelpers::split(...)
/// This is the reverse of DoubleVectorsHelpers::concatenate(...),
/// which is demonstrated in self_test/mpi/vector_concatenation/
/// Given a vector with length = 7+5+3 = 15
/// [1
///  2
///  3
///  4
///  5
///  6
///  7
///  8
///  9
///  10
///  11
///  12
///  13
///  14
///  15],
/// 
/// we split it into vectors v1, v2 and v3 of lengths 7, 5 and 3 respectively.
/// The script validate.sh should run this on 1, 2, 3 and 4 cores.
/// 
/// Communication is required and the order of the entries is preserved
/// across the vectors. We demonstrate this on two cores, p0 and p1:
/// v1, p0: [1
///          2
///          3]
/// v1, p1: [4
///          5
///          6
///          7]
///
/// v2, p0: [8
///          9]
/// v2, p1: [10
///          11
///          12]
///
/// v2, p0: [13]
/// v2, p1: [14
///          15]
//======================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Get the global oomph-lib communicator 
  const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();
  
  // The number of sub vectors
  unsigned nvectors = 3;

  // Supply the nrow for the sub distributions. In C++11, we can initialise
  // Vectors like this, but for now, we use arrays.
  unsigned nrowarray[] = {7,5,3};
  
  // The vectors to split the in vector in to.
  Vector<DoubleVector> out_vector(nvectors);
  bool distributed = true;

  // Create the vectors and get the pointers to them.
  // This may be a bit long winded but we do not want to create pointers
  // to objects using "new" since then we will have to remember to call
  // delete and null the pointers. This is considered bad practice, memory
  // management should be done automatically with smart pointers.
  // We could use smart pointers but it is only available in C++11,
  // or use boost smart pointers, but we do not have boost...
  // /rant.
  create_vectors_to_split(nrowarray,comm_pt,distributed,out_vector);

  // Global row for the in vector (must match the sum of the global rows for 
  // the sub vectors).
  unsigned in_nrow = 0;
  for (unsigned vec_i = 0; vec_i < nvectors; vec_i++) 
  {
    in_nrow += out_vector[vec_i].nrow();
  }
  
  // The in vector
  DoubleVector in_vector;
  create_vector_ascend_row(in_nrow, comm_pt, distributed, in_vector);
  
  // Call the split function.
  DoubleVectorHelpers::split(in_vector,out_vector);
  
  // The split is done, now we output the results.
  // We do not use the output function from DoubleVector
  // because it outputs the whole DoubleVector per processor 
  // (requiring communication). We want to check that the values per processor
  // is correct.

  // My rank
  unsigned my_rank = comm_pt->my_rank();
  unsigned nproc = comm_pt->nproc();

  // The output file name.
  std::ostringstream outfile_stream;
  outfile_stream << "out_NP"<<nproc<<"R"<< my_rank;

  // The output file.
  std::ofstream out_file;
  out_file.open(outfile_stream.str().c_str());

  // Output data from all out vectors:
  // nrow()
  // first_row()
  // nrow_local()
  // distributed()
  // values
  for (unsigned vec_i = 0; vec_i < nvectors; vec_i++) 
  {
    // The values and nrow local.
    double* out_values = out_vector[vec_i].values_pt();
    unsigned out_nrow_local = out_vector[vec_i].nrow_local();

    out_file << out_vector[vec_i].nrow() << "\n";
    out_file << out_vector[vec_i].first_row() << "\n";
    out_file << out_nrow_local << "\n";
    out_file << out_vector[vec_i].distributed() << "\n";
    
    for (unsigned val_i = 0; val_i < out_nrow_local; val_i++) 
    {
      out_file << out_values[val_i] << "\n";
    }
  }

  out_file.close();

#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main
