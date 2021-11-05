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

//===start_of_main======================================================
/// Driver code: Testing LinearAlgebraDistributionHelpers::concatenate(...)
///
/// We concatenate three uniformly distributed distributions 
/// with nrow = 7, 5, and 3 respectively.
///
/// On one core:
/// dist_0: nrow = 7
///         nrow_local rank0 = 7
///
/// dist_1: nrow = 5
///         nrow_local rank0 = 5
///
/// dist_2: nrow = 3
///         nrow_local rank0 = 3
///
/// The result distribution should have:
/// dist_r: nrow = 15
///         nrow_local rank0 = 15
///
/// ////////////////////////////////////
/// 
/// On two cores:
/// dist_0: nrow = 7
///         nrow_local rank0 = 3
///         nrow_local rank1 = 4
///
/// dist_1: nrow = 5
///         nrow_local rank0 = 2
///         nrow_local rank1 = 3
///
/// dist_2: nrow = 3
///         nrow_local rank0 = 1
///         nrow_local rank1 = 2
///
/// The result distribution should have:
/// dist_r: nrow = 15
///         nrow_local rank0 = 6
///         nrow_local rank1 = 9
///
/// ////////////////////////////////////
///
/// On three cores:
/// dist_0: nrow = 7
///         nrow_local rank0 = 2
///         nrow_local rank1 = 2
///         nrow_local rank2 = 3
///
/// dist_1: nrow = 5
///         nrow_local rank0 = 1
///         nrow_local rank1 = 2
///         nrow_local rank2 = 2
///
/// dist_2: nrow = 3
///         nrow_local rank0 = 1
///         nrow_local rank1 = 1
///         nrow_local rank2 = 1
///
/// The result distribution should have:
/// dist_r: nrow = 15
///         nrow_local rank0 = 4
///         nrow_local rank1 = 5
///         nrow_local rank2 = 6
///
/// ////////////////////////////////////
///
/// On four cores:
/// dist_0: nrow = 7
///         nrow_local rank0 = 1
///         nrow_local rank1 = 2
///         nrow_local rank2 = 2
///         nrow_local rank2 = 2
///
/// dist_1: nrow = 5
///         nrow_local rank0 = 1
///         nrow_local rank1 = 1
///         nrow_local rank1 = 1
///         nrow_local rank1 = 2
///
/// dist_2: nrow = 3
///         nrow_local rank0 = 0
///         nrow_local rank1 = 1
///         nrow_local rank1 = 1
///         nrow_local rank1 = 1
///
/// The result distribution should have:
/// dist_r: nrow = 15
///         nrow_local rank0 = 2
///         nrow_local rank1 = 4
///         nrow_local rank1 = 4 
///         nrow_local rank1 = 5 
///
/// The script validate.sh should run this test on 1, 2, 3 and4 cores.
//======================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);
#endif
  
  // Get the global oomph-lib communicator 
  const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

  // How many distributions do we want to generate?
  unsigned ndistributions = 3;

  // What are the lengths of the above distributions?
  unsigned distlengtharray[] = {7,5,3};

  // The distributions to concatenate.
  Vector<LinearAlgebraDistribution> dist_to_cat;
  bool distributed = true;

  // Create the distribution and get the pointers to them.
  // This may be a bit long winded but we do not want to create pointers
  // to distributions using "new", then we will have to remember to call
  // delete and null the pointers. This is considered bad practice, memory
  // management should be done automatically with smart pointers.
  // We could use smart pointers but it is only available in C++11,
  // or use boost smart pointers, but we do not have boost...
  // /rant.
  for (unsigned dist_i = 0; dist_i < ndistributions; dist_i++) 
  {
    // Create the distribution.
    dist_to_cat.push_back(LinearAlgebraDistribution(
          comm_pt,
          distlengtharray[dist_i],
          distributed));
  }
  
  // The pointers to distributions to concatenate.
  Vector<LinearAlgebraDistribution*> dist_to_cat_pt;
  for (unsigned dist_i = 0; dist_i < ndistributions; dist_i++) 
  {
    dist_to_cat_pt.push_back(&dist_to_cat[dist_i]);
  }

  // The result distribution.
  LinearAlgebraDistribution result_distribution;
  
  // Call the concatenate function.
  LinearAlgebraDistributionHelpers::concatenate(dist_to_cat_pt,
                                                result_distribution);
  
  // Output data about the result distribution:
  // nrow()
  // first_row()
  // nrow_local()
  // distributed()
  unsigned my_rank = comm_pt->my_rank();
  unsigned nproc = comm_pt->nproc();
  std::ostringstream result_stream;
  result_stream << "out_NP" << nproc << "R" << my_rank;
  
  std::ofstream result_file;
  result_file.open(result_stream.str().c_str());
  result_file << result_distribution.nrow() << "\n";
  result_file << result_distribution.first_row() << "\n";
  result_file << result_distribution.nrow_local() << "\n";
  result_file << result_distribution.distributed() << "\n";
  result_file.close();
  
#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main
