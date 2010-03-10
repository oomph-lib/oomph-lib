//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
#include "linear_algebra_distribution.h"

namespace oomph
{


 //============================================================================
 /// Sets the distribution. Takes first_row, local_nrow and 
 /// global_nrow as arguments. If global_nrow is not provided or equal to
 /// 0 then it is computed automatically
 //============================================================================
 void LinearAlgebraDistribution::build(const OomphCommunicator* 
				       const comm_pt,
				       const unsigned& first_row, 
				       const unsigned& local_nrow,
				       const unsigned& global_nrow)
 {
  // copy the communicator
  delete Comm_pt;
  Comm_pt = new OomphCommunicator(*comm_pt);

  // get the rank and the number of processors
  int my_rank = Comm_pt->my_rank();
  int nproc = Comm_pt->nproc();

  // resize the storage
  First_row.clear();
  First_row.resize(nproc);
  Nrow_local.clear();
  Nrow_local.resize(nproc);
  
  // set first row and local nrow for this processor
  First_row[my_rank] = first_row;
  Nrow_local[my_rank] = local_nrow;
  
#ifdef OOMPH_HAS_MPI
  // gather the First_row vector
  unsigned my_nrow_local = Nrow_local[my_rank];
  MPI_Allgather(&my_nrow_local,1,MPI_UNSIGNED,
                &Nrow_local[0],1,MPI_UNSIGNED,
                Comm_pt->mpi_comm());

  // gather the Nrow_local vector
  unsigned my_first_row = First_row[my_rank];
  MPI_Allgather(&my_first_row,1,MPI_UNSIGNED,
                &First_row[0],1,MPI_UNSIGNED,
                Comm_pt->mpi_comm());    
#endif

  // if global nrow is not provided then compute by summing local_nrow over
  // all processors
  if (global_nrow == 0)
   {
    if (nproc == 1)
     {
      Nrow = local_nrow;
     }
    else
     {
      Nrow = 0;
      for (int p = 0; p < nproc; p++)
       {
        Nrow += Nrow_local[p];
       }
     }
   }
  else
   {
    Nrow = global_nrow;
   }

  // the distribution is true
  Distributed = true;

#ifdef OOMPH_HAS_MPI  
#ifdef PARANOID
  // paranoid check that the distribution works


  // check that none of the processors partition overlap
  for (int p = 0; p < nproc; p++)
   {
    if (Nrow_local[p] > 0)
     {
      for (int pp = p+1; pp < nproc; pp++)
       {
        if (Nrow_local[pp] > 0)
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
                                  "LinearAlgebraDistribution::distribute(...)",
                                  OOMPH_EXCEPTION_LOCATION);
           }
         }
       }
     }
   }
  
  // check that no processor has a row with a global row index greater than
  // the number of global rows
  for (int p = 0; p < nproc; p++)
   {
    if (First_row[p] + Nrow_local[p] > Nrow)
     {
      std::ostringstream error_message;
      error_message << "Processor " << p << " contains rows "
                    << First_row[p] << " to " 
                    << First_row[p]+Nrow_local[p]-1
                    << " but there are only " << Nrow 
                    << " to be distributed." << std::endl;
      throw OomphLibWarning(error_message.str(),
                            "LinearAlgebraDistribution::distribute(...)",
                            OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif
#endif
 }
 
 //============================================================================
 /// \short Uniformly distribute global_nrow over all processors where 
 /// processors 0 holds approximately the first 
 /// global_nrow/n_proc, processor 1 holds the next 
 /// global_nrow/n_proc and so on...
 //============================================================================
 void LinearAlgebraDistribution::build(const OomphCommunicator* 
				       const comm_pt,
				       const unsigned& global_nrow,
				       const bool& distribute)
 {
  // copy the communicator
  delete Comm_pt;
  Comm_pt = new OomphCommunicator(*comm_pt);

  // delete existing storage
  First_row.clear();
  Nrow_local.clear();

  // set global nrow
  Nrow = global_nrow;
  
  // store the distributed flag
  Distributed = distribute;

#ifdef OOMPH_HAS_MPI

  // if distributed object then compute uniform distribution
  if (distribute == true)
   {

    // get the number of processors
    int nproc = Comm_pt->nproc();

    // resize the storage
    First_row.resize(nproc);
    Nrow_local.resize(nproc);

    // compute first row
    for (int p=0;p<nproc;p++)
     {                                          
      First_row[p] = unsigned(double(p*global_nrow)/ 
                              double(nproc));                  
     }             
    
    // compute local nrow
    for (int p=0; p<nproc-1; p++) 
     {                                                  
      Nrow_local[p] = First_row[p+1] - First_row[p];
     }  
    Nrow_local[nproc-1] = global_nrow -  
     First_row[nproc-1];         
   }
#endif
 }

 //============================================================================
 /// helper method for the =assignment operator and copy constructor
 //============================================================================
 void LinearAlgebraDistribution::build(const LinearAlgebraDistribution& 
				       new_dist)
  {

   // delete the existing storage
   First_row.clear();
   Nrow_local.clear();

   // if new_dist is not setup
   if (new_dist.communicator_pt() == 0)
    {
     delete Comm_pt;
     Comm_pt = 0;
     Distributed = true;
     Nrow = 0;
     }
   else
    {
     // copy the communicator
     delete Comm_pt;
     Comm_pt = new OomphCommunicator(*new_dist.communicator_pt());
         
     // the new distribution is distributed
     if (new_dist.distributed())
      {
       
       // # of processors
       int nproc = Comm_pt->nproc();
 
       // resize storage
       First_row.clear();
       First_row.resize(nproc);
       Nrow_local.clear();
       Nrow_local.resize(nproc);
       
       // copy contents of first_row and nrow_local
       for (int i = 0; i < nproc; i++)
        {
         First_row[i] = new_dist.first_row(i);
         Nrow_local[i] = new_dist.nrow_local(i); 
        }

       Distributed = true;
      }
     // else if the new ditribution is not distributed
     else
      {
       Distributed = false;
      }
     Nrow = new_dist.nrow();
    }
  }


 //============================================================================
 /// operator==
 //============================================================================
 bool LinearAlgebraDistribution::operator==
 (const LinearAlgebraDistribution& other_dist) const
  {   
#ifdef OOMPH_HAS_MPI
   // compare the communcators
   if (!((*Comm_pt) == (*other_dist.communicator_pt())))
    {
     return false;
    }

   // compare Distributed
   if (Distributed != other_dist.distributed())
    {
     return false;
    }

   // if not distributed compare nrow
   if (!Distributed)
    {
     if (other_dist.nrow() == Nrow)
      {
       return true;
      }
     return false;
    }
   
   // compare
   bool flag = true;
   int nproc = Comm_pt->nproc();
   for (int i = 0; i < nproc && flag == true; i++)
    {
     if (other_dist.first_row(i) != First_row[i] ||
         other_dist.nrow_local(i) != Nrow_local[i])
      {
       flag = false;
      }
    }
   return flag;
#else
   if (other_dist.nrow() == Nrow)
    {
     return true;
    }
   return false;
#endif
  }

 //=============================================================================
 /// output operator
 //=============================================================================
 std::ostream& operator<<(std::ostream& stream, 
                          LinearAlgebraDistribution dist)
 {
  stream << "nrow()=" << dist.nrow() 
         << ", first_row()=" << dist.first_row()
         << ", nrow_local()=" << dist.nrow_local()
         << ", distributed()=" << dist.distributed()
         << std::endl;
  return stream;
 }



}//end of oomph namespace
