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
#include "linear_algebra_distribution.h"

namespace oomph
{
  //============================================================================
  /// Sets the distribution. Takes first_row, local_nrow and
  /// global_nrow as arguments. If global_nrow is not provided or equal to
  /// 0 then it is computed automatically
  //============================================================================
  void LinearAlgebraDistribution::build(const OomphCommunicator* const comm_pt,
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
    MPI_Allgather(&my_nrow_local,
                  1,
                  MPI_UNSIGNED,
                  &Nrow_local[0],
                  1,
                  MPI_UNSIGNED,
                  Comm_pt->mpi_comm());

    // gather the Nrow_local vector
    unsigned my_first_row = First_row[my_rank];
    MPI_Allgather(&my_first_row,
                  1,
                  MPI_UNSIGNED,
                  &First_row[0],
                  1,
                  MPI_UNSIGNED,
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

    // the distribution is true, unless we only have 1 processor
    if (nproc != 1)
    {
      Distributed = true;
    }
    else
    {
      Distributed = false;
    }

#ifdef OOMPH_HAS_MPI
#ifdef PARANOID
    // paranoid check that the distribution works


    // check that none of the processors partition overlap
    for (int p = 0; p < nproc; p++)
    {
      if (Nrow_local[p] > 0)
      {
        for (int pp = p + 1; pp < nproc; pp++)
        {
          if (Nrow_local[pp] > 0)
          {
            if ((First_row[p] >= First_row[pp] &&
                 First_row[p] < First_row[pp] + Nrow_local[pp]) ||
                (First_row[p] + Nrow_local[p] - 1 >= First_row[pp] &&
                 First_row[p] + Nrow_local[p] - 1 <
                   First_row[pp] + Nrow_local[pp]))
            {
              std::ostringstream error_message;
              error_message
                << "The distributed rows on processor " << p
                << " and processor " << pp << " overlap.\n"
                << "Processor " << p << " : first_row = " << First_row[p]
                << ", nrow = " << Nrow_local[p] << ".\n"
                << "Processor " << pp << " : first_row = " << First_row[pp]
                << ", nrow = " << Nrow_local[pp] << ".\n";
              throw OomphLibWarning(
                error_message.str(),
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
        error_message << "Processor " << p << " contains rows " << First_row[p]
                      << " to " << First_row[p] + Nrow_local[p] - 1
                      << " but there are only " << Nrow << " to be distributed."
                      << std::endl;
        throw OomphLibWarning(error_message.str(),
                              "LinearAlgebraDistribution::distribute(...)",
                              OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif
#endif
  }

  //============================================================================
  /// Uniformly distribute global_nrow over all processors where
  /// processors 0 holds approximately the first
  /// global_nrow/n_proc, processor 1 holds the next
  /// global_nrow/n_proc and so on...
  //============================================================================
  void LinearAlgebraDistribution::build(const OomphCommunicator* const comm_pt,
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
      for (int p = 0; p < nproc; p++)
      {
        First_row[p] =
          unsigned((double(p) / double(nproc)) * double(global_nrow));
      }

      // compute local nrow
      for (int p = 0; p < nproc - 1; p++)
      {
        Nrow_local[p] = First_row[p + 1] - First_row[p];
      }
      Nrow_local[nproc - 1] = global_nrow - First_row[nproc - 1];
    }
#endif
  }

  //============================================================================
  /// helper method for the =assignment operator and copy constructor
  //============================================================================
  void LinearAlgebraDistribution::build(
    const LinearAlgebraDistribution& new_dist)
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
        First_row = new_dist.first_row_vector();
        Nrow_local = new_dist.nrow_local_vector();

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
  bool LinearAlgebraDistribution::operator==(
    const LinearAlgebraDistribution& other_dist) const
  {
#ifdef OOMPH_HAS_MPI
    // compare the communicators
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
                           LinearAlgebraDistribution& dist)
  {
    stream << "nrow()=" << dist.nrow() << ", first_row()=" << dist.first_row()
           << ", nrow_local()=" << dist.nrow_local()
           << ", distributed()=" << dist.distributed() << std::endl;
    return stream;
  }

  //=============================================================================
  /// Namespace for helper functions for LinearAlgebraDistributions
  //=============================================================================
  namespace LinearAlgebraDistributionHelpers
  {
    //===========================================================================
    /// Takes a vector of LinearAlgebraDistribution objects and
    /// concatenates them such that the nrow_local of the out_distribution
    /// is the sum of the nrow_local of all the in_distributions and the number
    /// of global rows of the out_distribution is the sum of the number of
    /// global rows of all the in_distributions. This results in a permutation
    /// of the rows in the out_distribution. Think of this in terms of
    /// DoubleVectors, if we have DoubleVectors with distributions A and B,
    /// distributed across two processors (p0 and p1), A: [a0] (on p0)    B:
    /// [b0] (on p0)
    ///    [a1] (on p1)       [b1] (on P1),
    ///
    /// then the out_distribution is
    /// [a0  (on p0)
    ///  b0] (on p0)
    /// [a1  (on p1)
    ///  b1] (on p1),
    ///
    /// as opposed to
    /// [a0  (on p0)
    ///  a1] (on p0)
    /// [b0  (on p1)
    ///  b1] (on p1).
    ///
    /// Note (1): The out_distribution may not be uniformly distributed even
    /// if the in_distributions are uniform distributions.
    /// Try this out with two distributions of global rows 3 and 5, uniformly
    /// distributed across two processors. Compare this against a distribution
    /// of global row 8 distributed across two processors.
    ///
    /// Note (2): There is no equivalent function which takes a Vector of
    /// LinearAlgebraDistribution objects (as opposed to pointers), there should
    /// not be one since we do not want to invoke the assignment operator when
    /// creating the Vector of LinearAlgebraDistribution objects.
    //===========================================================================
    void concatenate(
      const Vector<LinearAlgebraDistribution*>& in_distribution_pt,
      LinearAlgebraDistribution& out_distribution)
    {
      // How many distributions are in in_distribution?
      unsigned ndistributions = in_distribution_pt.size();

#ifdef PARANOID

      // Are any of the in_distribution pointers null?
      // We do not want null pointers.
      for (unsigned dist_i = 0; dist_i < ndistributions; dist_i++)
      {
        if (in_distribution_pt[dist_i] == 0)
        {
          std::ostringstream error_message;
          error_message << "The pointer in_distribution_pt[" << dist_i
                        << "] is null.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      /// /////////////////////////////

      // Check that all in distributions are built.
      for (unsigned dist_i = 0; dist_i < ndistributions; dist_i++)
      {
        if (!in_distribution_pt[dist_i]->built())
        {
          std::ostringstream error_message;
          error_message << "The in_distribution_pt[" << dist_i
                        << "] is not built.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      /// /////////////////////////////

      // Check that all communicators to concatenate are the same
      // by comparing all communicators against the first one.
      const OomphCommunicator first_comm =
        *(in_distribution_pt[0]->communicator_pt());

      for (unsigned dist_i = 0; dist_i < ndistributions; dist_i++)
      {
        // Get the Communicator for the current vector.
        const OomphCommunicator another_comm =
          *(in_distribution_pt[dist_i]->communicator_pt());

        if (!(first_comm == another_comm))
        {
          std::ostringstream error_message;
          error_message << "The communicator in position " << dist_i << " is \n"
                        << "not the same as the first one.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      /// /////////////////////////////

      // Ensure that all distributions are either distributed or not.
      // This is because we use the distributed() function from the first
      // distribution to indicate if the result distribution should be
      // distributed or not.
      bool first_distributed = in_distribution_pt[0]->distributed();
      for (unsigned dist_i = 0; dist_i < ndistributions; dist_i++)
      {
        // Is the current distribution distributed?
        bool another_distributed = in_distribution_pt[dist_i]->distributed();
        if (first_distributed != another_distributed)
        {
          std::ostringstream error_message;
          error_message
            << "The distribution in position " << dist_i << " has a different\n"
            << "different distributed boolean than the distribution.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      /// /////////////////////////////

      // Check that the out distribution is not built.
      if (out_distribution.built())
      {
        std::ostringstream error_message;
        error_message << "The out distribution is built.\n"
                      << "Please clear it.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif

      // Get the communicator pointer
      const OomphCommunicator* const comm_pt =
        in_distribution_pt[0]->communicator_pt();

      // Number of processors
      unsigned nproc = comm_pt->nproc();

      // First determine the out_nrow_local and the out_nrow (the global nrow)
      unsigned out_nrow_local = 0;
      unsigned out_nrow = 0;
      for (unsigned in_dist_i = 0; in_dist_i < ndistributions; in_dist_i++)
      {
        out_nrow_local += in_distribution_pt[in_dist_i]->nrow_local();
        out_nrow += in_distribution_pt[in_dist_i]->nrow();
      }

      // Now calculate the first row for this processor. We need to know the
      // out_nrow_local for all the other processors, MPI_Allgather must be
      // used.
      unsigned out_first_row = 0;

      // Distributed case: We need to gather all the out_nrow_local from all
      // processors, the out_first_row for this processors is the partial sum up
      // of the out_nrow_local to my_rank.
      bool distributed = in_distribution_pt[0]->distributed();
      if (distributed)
      {
#ifdef OOMPH_HAS_MPI
        // My rank
        unsigned my_rank = comm_pt->my_rank();

        unsigned* out_nrow_local_all = new unsigned[nproc];
        MPI_Allgather(&out_nrow_local,
                      1,
                      MPI_UNSIGNED,
                      out_nrow_local_all,
                      1,
                      MPI_UNSIGNED,
                      comm_pt->mpi_comm());
        for (unsigned proc_i = 0; proc_i < my_rank; proc_i++)
        {
          out_first_row += out_nrow_local_all[proc_i];
        }
        delete[] out_nrow_local_all;
#endif
      }
      // Non-distributed case
      else
      {
        out_first_row = 0;
      }

      // Build the distribution.
      if (nproc == 1 || !distributed)
      {
        // Some ambiguity here -- on a single processor it doesn't really
        // matter if we think of ourselves as distributed or not but this
        // follows the pattern used elsewhere.
        out_distribution.build(comm_pt, out_nrow, false);
      }
      else
      {
        out_distribution.build(
          comm_pt, out_first_row, out_nrow_local, out_nrow);
      }
    } // function concatenate

  } // end of namespace LinearAlgebraDistributionHelpers

} // namespace oomph
