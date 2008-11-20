//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
#include "distribution_info.h"

namespace oomph
{

 //============================================================================
 /// Sets the distribution. Takes first_row, local_nrow and 
 /// global_nrow as arguments. If global_nrow is not provided or equal to
 /// 0 then it is computed automatically
 //============================================================================
 void DistributionInfo::distribute(MPI_Comm comm,
                                   const long unsigned& first_row, 
                                   const long unsigned& local_nrow,
                                   const long unsigned& global_nrow)
 {
  // copy the communicator
  Communicator = comm;

  // get the rank and the number of processors
  int my_rank;
  MPI_Comm_rank(Communicator,&my_rank);
  int n_proc;
  MPI_Comm_size(Communicator,&n_proc);

  // resize the storage
  First_row.resize(n_proc);
  Nrow_local.resize(n_proc);
  
  // set first row and local nrow for this processor
  First_row[my_rank] = first_row;
  Nrow_local[my_rank] = local_nrow;
  
  // gather the First_row vector
  MPI_Allgather(&Nrow_local[my_rank],1,
                MPI_UNSIGNED_LONG,
                &Nrow_local[0],1,MPI_UNSIGNED_LONG,
                Communicator);

  // gather the Nrow_local vector
  MPI_Allgather(&First_row[my_rank],1,
                MPI_UNSIGNED_LONG,
                &First_row[0],1,MPI::UNSIGNED_LONG,
                Communicator);    

  // if global nrow is not provided then compute by summing local_nrow over
  // all processors
  if (global_nrow == 0)
   {
    Nrow_global = 0;
    for (int p = 0; p < n_proc; p++)
     {
      Nrow_global += Nrow_local[p];
     }
   }
  else
   {
    Nrow_global = global_nrow;
   }
  
  // distribution has been setup
  Setup = true;
  
#ifdef PARANOID
  // paranoid check that the distribution works

  // check that none of the processors partition overlap
  for (int p = 0; p < n_proc; p++)
   {
    for (int pp = p+1; pp < n_proc; pp++)
     {
      if ((First_row[p] >= First_row[pp] &&
           First_row[p] < First_row[pp] + Nrow_local[pp]) ||
          (First_row[p] + Nrow_local[p] -1 >= First_row[pp] &&
           First_row[p] + Nrow_local[p] -1 <
           First_row[pp] + Nrow_local[pp]))
       {
        std::ostringstream error_message;
        error_message << "The distributed rows on processor " << p
                      << " and processor " << pp << " overlap.\n"
                      << "Processor " << p << " : first_row = "
                      << First_row[p] << ", nrow = " 
                      << Nrow_local[p] << ".\n"
                      << "Processor " << pp << " : first_row = "
                      << First_row[pp] << ", nrow = " 
                      << Nrow_local[pp] << ".\n";
        throw OomphLibWarning(error_message.str(),
                              "DistributionInfo::distribute(...)",
                              OOMPH_EXCEPTION_LOCATION);
       }
     }
   }
  
  // check that no processor has a row with a global row index greater than
  // the number of global rows
  for (int p = 0; p < n_proc; p++)
   {
    if (First_row[p] + Nrow_local[p] > Nrow_global)
     {
      std::ostringstream error_message;
      error_message << "Processor " << p << " contains rows "
                    << First_row[p] << " to " 
                    << First_row[p]+Nrow_local[p]-1
                    << " but there are only " << Nrow_global 
                    << " to be distributed." << std::endl;
      throw OomphLibWarning(error_message.str(),
                            "DistributionInfo::distribute(...)",
                            OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif
 }
 
 //============================================================================
 /// \short Uniformly distribute global_nrow over all processors where 
 /// processors 0 holds approximately the first 
 /// global_nrow/n_proc, processor 1 holds the next 
 /// global_nrow/n_proc and so on...
 //============================================================================
 void DistributionInfo::distribute(MPI_Comm comm,
                                   const long unsigned& global_nrow)
 {
  // copy the communicator
  Communicator = comm;

  // the number of processors in this distribution
  int n_proc;
  MPI_Comm_size(Communicator,&n_proc);

  // set global nrow
  Nrow_global = global_nrow;
  
  // resize the vectors
  Nrow_local.clear();
  First_row.clear();
  Nrow_local.resize(n_proc);                               
  First_row.resize(n_proc);
  
  // compute first row
  for (int p=0;p<n_proc;p++)
   {                                          
    First_row[p] = unsigned(double(p*global_nrow)/ 
                            double(n_proc));                  
   }             
  
  // compute local nrow
  for (int p=0; p<n_proc-1; p++) 
   {                                                  
    Nrow_local[p] = First_row[p+1] - First_row[p];
   }  
  Nrow_local[n_proc-1] = global_nrow -  
   First_row[n_proc-1];          
  
  // distribution has been setup
  Setup = true;
 }
}
