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


// Helper function to create a Vector of DoubleVectors 
// with uniform distribution and fills them with natural numbers.
void create_vectors_to_split
(unsigned nrowarray[], const OomphCommunicator* const comm_pt, 
 const bool distributed, Vector<DoubleVector>&out_vector)
{
  // The number of out vectors.
  unsigned nvectors = out_vector.size();

  // Fill in each sub DoubleVector using create_double_vector_ascend_row(...)
  for (unsigned vec_i = 0; vec_i < nvectors; vec_i++) 
  {
    // Create the distribution
    LinearAlgebraDistribution dist(comm_pt,nrowarray[vec_i],distributed);
    out_vector[vec_i].build(dist);
  }
}

//===start_of_main======================================================
/// Driver code: Testing 
/// DoubleVectorHelpers::split_without_communication(...)
/// This is the reverse of 
/// DoubleVectorHelpers::concatenate_without_communication(...), demonstrated
/// in self_test/mpi/vector_concatenation_without_communication/
///
/// There is a "strong" relationship between the distributions of the in_vector 
/// and the out vectors: The distribution of the in_vector must be the same
/// as the concatenation of the out vectors as defined by
/// LinearAlgebraDistributionHelpers::concatenate(...).
///
/// To see why, we demonstrate this on two cores, p0 and p1, and 
/// out vectors v1, v2 and v3 with lengths 7, 5 and 3 respectively.
/// As seen in the driver code distribution_concatenation, we have:
///
/// v1: nrow = 7
///     nrow_local p0 = 3
///     nrow_local p1 = 4
///
/// v2: nrow = 5
///     nrow_local p0 = 2
///     nrow_local p1 = 3
///
/// v3: nrow = 3
///     nrow_local p0 = 1
///     nrow_local p1 = 2
///
/// If no communication takes place, i.e. all data on p0 stays on p0 
/// (and similarly for p1), the in vector MUST have the following distribution:
/// in_vec: nrow = 15
///         nrow_local p0 = 3+2+1 = 6
///         nrow_local p1 = 4+3+2 = 9
///
/// In this driver test, we split a vector of length 15 into vectors of
/// length 7, 5 and 3, ensuring that the above relationship is maintained. 
/// The script validate.sh should run this program on 1, 2, 3 and 4 cores.
//======================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);
#endif
  
  // Get the global oomph-lib communicator 
  const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

  //  my rank
  unsigned my_rank = comm_pt->my_rank();
  
  // number of out vectors
  unsigned nvectors = 3;

  // Supply the number of global rows of the vectors
  unsigned nrowarray[] = {7,5,3};
  
  // To create an in_vector with the correct distribution, we first create
  // the out vectors, then store their distribution pointers in a Vector.

  // The out vectors
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

  // Create the vectors to split the vector into.
  create_vectors_to_split(nrowarray,comm_pt,distributed,out_vector);
  
  // Pointers to the out distributions
  Vector<LinearAlgebraDistribution*> out_distribution_pt(nvectors,0);
  for (unsigned vec_i = 0; vec_i < nvectors; vec_i++) 
  {
    out_distribution_pt[vec_i] = out_vector[vec_i].distribution_pt();
  }

  // Create the in vector distribution by concatenating the distributions
  // from out vectors.
  LinearAlgebraDistribution in_distribution;
  LinearAlgebraDistributionHelpers::concatenate(out_distribution_pt,
                                                in_distribution);

  // Create the in vector with natural numbers.
  DoubleVector in_vector(in_distribution,0.0);
  
  double* in_values_pt = in_vector.values_pt();
  unsigned in_nrow_local = in_vector.nrow_local();
  unsigned in_first_row = in_vector.first_row();

  // Fill in the values:
  for (unsigned row_i = 0; row_i < in_nrow_local; row_i++) 
  {
    in_values_pt[row_i] = in_first_row + row_i + 1;
  }

  // Now split the vector.
  DoubleVectorHelpers::split_without_communication(in_vector,
                                                   out_vector);
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
