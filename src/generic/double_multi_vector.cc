// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
#include "double_multi_vector.h"
#include "matrices.h"

namespace oomph
{
  /// The contents of the vector are redistributed to match the new
  /// distribution. In a non-MPI rebuild this method works, but does nothing.
  /// \b NOTE 1: The current distribution and the new distribution must have
  /// the same number of global rows.
  /// \b NOTE 2: The current distribution and the new distribution must have
  /// the same Communicator.
  void DoubleMultiVector::redistribute(
    const LinearAlgebraDistribution* const& dist_pt)
  {
#ifdef OOMPH_HAS_MPI
#ifdef PARANOID
    if (!Internal_values)
    {
      // if this vector does not own the double* values then it cannot be
      // distributed.
      // note: this is not stictly necessary - would just need to be careful
      // with delete[] below.
      std::ostringstream error_message;
      error_message
        << "This multi vector does not own its data (i.e. data has been "
        << "passed in via set_external_values() and therefore "
        << "cannot be redistributed";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // paranoid check that the nrows for both distributions is the
    // same
    if (dist_pt->nrow() != this->nrow())
    {
      std::ostringstream error_message;
      error_message << "The number of global rows in the new distribution ("
                    << dist_pt->nrow() << ") is not equal to the number"
                    << " of global rows in the current distribution ("
                    << this->nrow() << ").\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // paranoid check that the current distribution and the new distribution
    // have the same Communicator
    OomphCommunicator temp_comm(*dist_pt->communicator_pt());
    if (!(temp_comm == *this->distribution_pt()->communicator_pt()))
    {
      std::ostringstream error_message;
      error_message << "The new distribution and the current distribution must "
                    << "have the same communicator.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // check the distributions are not the same
    if (!((*this->distribution_pt()) == *dist_pt))
    {
      // Cache the number of vectors
      const unsigned n_vector = this->Nvector;

      // get the rank and the number of processors
      int my_rank = this->distribution_pt()->communicator_pt()->my_rank();
      int nproc = this->distribution_pt()->communicator_pt()->nproc();

      // if both vectors are distributed
      if (this->distributed() && dist_pt->distributed())
      {
        // new nrow_local and first_row data
        Vector<unsigned> new_first_row_data(nproc);
        Vector<unsigned> new_nrow_local_data(nproc);
        Vector<unsigned> current_first_row_data(nproc);
        Vector<unsigned> current_nrow_local_data(nproc);
        for (int i = 0; i < nproc; i++)
        {
          new_first_row_data[i] = dist_pt->first_row(i);
          new_nrow_local_data[i] = dist_pt->nrow_local(i);
          current_first_row_data[i] = this->first_row(i);
          current_nrow_local_data[i] = this->nrow_local(i);
        }

        // compute which local rows are expected to be received from each
        // processor / sent to each processor
        Vector<unsigned> new_first_row_for_proc(nproc);
        Vector<unsigned> new_nrow_local_for_proc(nproc);
        Vector<unsigned> new_first_row_from_proc(nproc);
        Vector<unsigned> new_nrow_local_from_proc(nproc);

        // for every processor compute first_row and nrow_local that will
        // will sent and received by this processor
        for (int p = 0; p < nproc; p++)
        {
          // start with data to be sent
          if ((new_first_row_data[p] < (current_first_row_data[my_rank] +
                                        current_nrow_local_data[my_rank])) &&
              (current_first_row_data[my_rank] <
               (new_first_row_data[p] + new_nrow_local_data[p])))
          {
            new_first_row_for_proc[p] =
              std::max(current_first_row_data[my_rank], new_first_row_data[p]);
            new_nrow_local_for_proc[p] =
              std::min((current_first_row_data[my_rank] +
                        current_nrow_local_data[my_rank]),
                       (new_first_row_data[p] + new_nrow_local_data[p])) -
              new_first_row_for_proc[p];
          }

          // and data to be received
          if ((new_first_row_data[my_rank] <
               (current_first_row_data[p] + current_nrow_local_data[p])) &&
              (current_first_row_data[p] <
               (new_first_row_data[my_rank] + new_nrow_local_data[my_rank])))
          {
            new_first_row_from_proc[p] =
              std::max(current_first_row_data[p], new_first_row_data[my_rank]);
            new_nrow_local_from_proc[p] =
              std::min(
                (current_first_row_data[p] + current_nrow_local_data[p]),
                (new_first_row_data[my_rank] + new_nrow_local_data[my_rank])) -
              new_first_row_from_proc[p];
          }
        }

        // Storage for the new data
        double** temp_data = new double*[n_vector];
        double* contiguous_temp_data =
          new double[n_vector * new_nrow_local_data[my_rank]];
        for (unsigned v = 0; v < n_vector; ++v)
        {
          temp_data[v] =
            &contiguous_temp_data[v * new_nrow_local_data[my_rank]];
        }

        // "send to self" or copy Data that does not need to be sent else where
        // to temp_data
        if (new_nrow_local_for_proc[my_rank] != 0)
        {
          unsigned j =
            new_first_row_for_proc[my_rank] - current_first_row_data[my_rank];
          unsigned k =
            new_first_row_for_proc[my_rank] - new_first_row_data[my_rank];
          for (unsigned i = 0; i < new_nrow_local_for_proc[my_rank]; i++)
          {
            for (unsigned v = 0; v < n_vector; ++v)
            {
              temp_data[v][k + i] = Values[v][j + i];
            }
          }
        }

        // send and receive circularly
        for (int p = 1; p < nproc; p++)
        {
          // next processor to send to
          unsigned dest_p = (my_rank + p) % nproc;

          // next processor to receive from
          unsigned source_p = (nproc + my_rank - p) % nproc;

          // send and receive the value
          MPI_Status status;
          for (unsigned v = 0; v < n_vector; v++)
          {
            MPI_Sendrecv(Values[v] + new_first_row_for_proc[dest_p] -
                           current_first_row_data[my_rank],
                         new_nrow_local_for_proc[dest_p],
                         MPI_DOUBLE,
                         dest_p,
                         1,
                         temp_data[v] + new_first_row_from_proc[source_p] -
                           new_first_row_data[my_rank],
                         new_nrow_local_from_proc[source_p],
                         MPI_DOUBLE,
                         source_p,
                         1,
                         this->distribution_pt()->communicator_pt()->mpi_comm(),
                         &status);
          }
        }

        // copy from temp data to Values_pt
        delete[] Values[0];
        delete[] Values;
        Values = temp_data;
      }
      // if this vector is distributed but the new distributed is global
      else if (this->distributed() && !dist_pt->distributed())
      {
        // copy existing Values_pt to temp_data
        unsigned n_local_data = this->nrow_local();
        double** temp_data = new double*[n_vector];
        // New continguous data
        double* contiguous_temp_data = new double[n_vector * n_local_data];
        for (unsigned v = 0; v < n_vector; ++v)
        {
          temp_data[v] = &contiguous_temp_data[v * n_local_data];
          for (unsigned i = 0; i < n_local_data; i++)
          {
            temp_data[v][i] = Values[v][i];
          }
        }

        // clear and resize Values_pt
        delete[] Values[0];
        double* values = new double[this->nrow() * n_vector];
        for (unsigned v = 0; v < n_vector; v++)
        {
          Values[v] = &values[v * this->nrow()];
        }

        // create a int vector of first rows
        int* dist_first_row = new int[nproc];
        int* dist_nrow_local = new int[nproc];
        for (int p = 0; p < nproc; p++)
        {
          dist_first_row[p] = this->first_row(p);
          dist_nrow_local[p] = this->nrow_local(p);
        }

        // gather the local vectors from all processors on all processors
        int my_local_data(this->nrow_local());

        // Loop over all vectors
        for (unsigned v = 0; v < n_vector; v++)
        {
          MPI_Allgatherv(
            temp_data[v],
            my_local_data,
            MPI_DOUBLE,
            Values[v],
            dist_nrow_local,
            dist_first_row,
            MPI_DOUBLE,
            this->distribution_pt()->communicator_pt()->mpi_comm());
        }

        // update the distribution
        this->build_distribution(dist_pt);

        // delete the temp_data
        delete[] temp_data[0];
        delete[] temp_data;

        // clean up
        delete[] dist_first_row;
        delete[] dist_nrow_local;
      }

      // if this vector is not distrubted but the target vector is
      else if (!this->distributed() && dist_pt->distributed())
      {
        // cache the new nrow_local
        unsigned nrow_local = dist_pt->nrow_local();

        // and first_row
        unsigned first_row = dist_pt->first_row();

        const unsigned n_local_data = nrow_local;
        double** temp_data = new double*[n_vector];
        double* contiguous_temp_data = new double[n_vector * n_local_data];

        // copy the data
        for (unsigned v = 0; v < n_vector; v++)
        {
          temp_data[v] = &contiguous_temp_data[v * n_local_data];
          for (unsigned i = 0; i < n_local_data; i++)
          {
            temp_data[v][i] = Values[v][first_row + i];
          }
        }

        // copy to Values_pt
        delete[] Values[0];
        delete[] Values;
        Values = temp_data;

        // update the distribution
        this->build_distribution(dist_pt);
      }

      // copy the Distribution
      this->build_distribution(dist_pt);
    }
#endif

    // Update the doublevector representation
    this->setup_doublevector_representation();
  }


} // namespace oomph
