//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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

// Given dimvector, which describes the sizes of each DoubleVector,
// returns a Vector of DoubleVector (as described by dimvector).
void fill_in_sub_vectors(const OomphCommunicator* const comm_pt, 
                         const bool distributed,
                         Vector<unsigned > &dimvector,
                         Vector<DoubleVector> &my_vecs)
{
  unsigned nvectors = my_vecs.size();

  for (unsigned vec_i = 0; vec_i < nvectors; vec_i++) 
  {
    unsigned vec_size = dimvector[vec_i];
    create_vector_ascend_row(vec_size,comm_pt,distributed,
                             my_vecs[vec_i]);
  }
}

// Helper funcion to call other helper functions
// 1) converts dimarray to dimvector for robustness and ease of management.
// 2) fill_in_sub_vectors(...) creates a distributed vector for each dim
//    in dimvector.
void create_vectors_to_cat
(const unsigned nvectors, unsigned dimarray[], 
 const OomphCommunicator* const comm_pt, Vector<DoubleVector>&my_vecs)
{
  // Convert dimarray into a Vector so it is easier to manage.
  Vector<unsigned> dimvector;
  construct_vector(dimarray,nvectors,dimvector);
  
  // Fill in each sub DoubleVector using create_double_vector_ascend_row(...)
  bool distributed = true;
  fill_in_sub_vectors(comm_pt,distributed,dimvector,my_vecs);
}

//===start_of_main======================================================
/// Driver code: Testing 
/// DoubleVectorHelpers::concatenate_without_communication(...)
/// 
/// We concatenate uniformly distributed DoubleVectors.
/// The vectors v1, v2, and v3 have increasing entries 
/// with lengths 7, 5 and 3 as depicted below:
///
/// v1   v2   v3
/// [1   [1   [1
///  2    2    2
///  3    3    3]
///  4    4
///  5    5]
///  6
///  7]
///
/// The script validate.sh should run this program on 1, 2, 3 and 4 cores. 
/// 
/// Communication is NOT required but the order of the entries is NOT 
/// preserved. We demonstrate this on two cores, p0 and p1:
/// v1 p0   p1
///    [1   [4
///     2    5
///     3]   6
///          7]
///
/// v2 p0   p1
///    [1   [3
///     2]   4
///          5]
///
/// v3 p0   p1
///    [1]  [2
///          3]
/// 
/// Result vector:
///    p0   p1
///    [1   [4
///     2    5
///     3    6
///     1    7
///     2    3
///     1]   4
///          5
///          2
///          3]
//======================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);
#endif
  
  // Get the global oomph-lib communicator 
  const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

  //  My rank
  unsigned my_rank = comm_pt->my_rank();
  
  // number of sub vectors
  unsigned nvectors = 3;

  // Supply the length of the vectors
  unsigned nrowarray[] = {7,5,3};
  
  // The data structure to store the DoubleVectors to concatenate
  Vector<DoubleVector> in_vector(nvectors);
  
  // Create the vectors to concatenate.
  create_vectors_to_cat(nvectors,nrowarray,comm_pt,in_vector);
  
  // The out vector.
  DoubleVector out_vector;
  
  // Call the concatenate function.
  DoubleVectorHelpers::concatenate_without_communication(in_vector,
                                                         out_vector);
  
  // Output data from out_vector:
  // nrow()
  // first_row
  // nrow_local()
  // distributed()
  // values()
  unsigned nproc = comm_pt->nproc();
  std::ostringstream out_stream;
  out_stream << "out_NP"<<nproc<<"R"<< my_rank;

  std::ofstream out_file;
  out_file.open(out_stream.str().c_str());

  double* out_values = out_vector.values_pt();
  unsigned out_nrow_local = out_vector.nrow_local();
  
  out_file << out_vector.nrow() << "\n";
  out_file << out_vector.first_row() << "\n";
  out_file << out_nrow_local << "\n";
  out_file << out_vector.distributed() << "\n";

  for (unsigned val_i = 0; val_i < out_nrow_local; val_i++) 
  {
    out_file << out_values[val_i] << "\n";
  }

  out_file.close();

#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main
