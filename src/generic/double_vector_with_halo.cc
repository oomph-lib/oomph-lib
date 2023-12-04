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
#include "double_vector_with_halo.h"

namespace oomph
{
  //============================================================================
  /// Constructor that sets up the required information communicating
  /// between all processors. Requires two "all to all" communications.
  /// Arguments are the distribution of the DoubleVector and a
  /// Vector of global unknowns required on this processor.
  //===========================================================================
  DoubleVectorHaloScheme::DoubleVectorHaloScheme(
    LinearAlgebraDistribution* const& dist_pt,
    const Vector<unsigned>& required_global_eqn)
    : Distribution_pt(dist_pt)
  {
#ifdef OOMPH_HAS_MPI
    // Only bother to do anything if the vector is distributed
    if (dist_pt->distributed())
    {
      // First create temporary maps for halo requests.
      // Using the map structure ensures that the data will be sorted
      // into processor order: the order of the unsigned integer key

      // These are the halo requests indexed by master processor and the
      // local equation number that is to be haloed on that processor
      std::map<unsigned, Vector<unsigned>> to_be_haloed;
      // These are the halo requests indexed by master processor and the
      // index in the additional storage
      // that corresponds to the halo data on this processor
      std::map<unsigned, Vector<unsigned>> halo_entries;

      // Find rank of current processor
      const unsigned my_rank =
        static_cast<int>(dist_pt->communicator_pt()->my_rank());
      // Loop over the desired equations (on this processor)
      // and work out whether we need to halo them according to the
      // given distribution
      const unsigned n_global_eqn = required_global_eqn.size();
      // Index for the locally stored halo values
      unsigned index = 0;
      for (unsigned n = 0; n < n_global_eqn; n++)
      {
        // Cache the required GLOBAL equation number
        const unsigned i_global = required_global_eqn[n];
        // Where is it stored?
        unsigned rank_of_global = dist_pt->rank_of_global_row(i_global);
        // If the equation is not stored locally then
        // populate the two maps
        if (my_rank != rank_of_global)
        {
          // Work out the local entry on the appropriate processor
          unsigned i_local = i_global - dist_pt->first_row(rank_of_global);
          // Mark the local storage index as halo with rank_of_global as master
          halo_entries[rank_of_global].push_back(index);
          // Mark the local equation of the rank_of_global as to be
          // haloed
          to_be_haloed[rank_of_global].push_back(i_local);
          // Store the local index corresponding to the global equation
          Local_index[i_global] = index;
          // Increment the index
          ++index;
        }
      }

      // We now need to tell the other processors which of their data are
      // haloed on this processor

      // First find out how many processors there are!
      const int n_proc = dist_pt->communicator_pt()->nproc();

      // Setup storage for number of data to be sent to each processor
      Vector<int> send_n(n_proc, 0);
      Vector<int> send_displacement(n_proc, 0);
      int send_data_count = 0;

      // Iterate over the entries in the map
      // This will be in rank order because of the ordering of the map
      for (std::map<unsigned, Vector<unsigned>>::iterator it =
             to_be_haloed.begin();
           it != to_be_haloed.end();
           ++it)
      {
        const unsigned rank = it->first;
        const unsigned size_ = it->second.size();
        // The displacement is the current number of data
        send_displacement[rank] = send_data_count;
        // The number to send is just the size of the array
        send_n[rank] = static_cast<int>(size_);
        send_data_count += size_;
      }

      // Now send the number of haloed entries from every processor
      // to every processor

      // Receive the data directly into the storage for haloed daa
      Haloed_n.resize(n_proc, 0);
      MPI_Alltoall(&send_n[0],
                   1,
                   MPI_INT,
                   &Haloed_n[0],
                   1,
                   MPI_INT,
                   dist_pt->communicator_pt()->mpi_comm());

      // Prepare the data to be sent
      // Always resize to at least one
      if (send_data_count == 0)
      {
        ++send_data_count;
      }
      Vector<unsigned> send_data(send_data_count);
      // Iterate over the entries in the map
      unsigned count = 0;
      for (std::map<unsigned, Vector<unsigned>>::iterator it =
             to_be_haloed.begin();
           it != to_be_haloed.end();
           ++it)
      {
        // Iterate over the vector
        for (Vector<unsigned>::iterator it2 = it->second.begin();
             it2 != it->second.end();
             ++it2)
        {
          send_data[count] = (*it2);
          ++count;
        }
      }

      // Prepare the data to be received,
      // Again this can go directly into Haloed storage
      int receive_data_count = 0;
      Haloed_displacement.resize(n_proc);
      for (int d = 0; d < n_proc; d++)
      {
        // The displacement is the amount of data received so far
        Haloed_displacement[d] = receive_data_count;
        receive_data_count += Haloed_n[d];
      }

      // Now resize the receive buffer
      // Always make sure that it has size of at least one
      if (receive_data_count == 0)
      {
        ++receive_data_count;
      }
      Haloed_eqns.resize(receive_data_count);
      // Send the data between all the processors
      MPI_Alltoallv(&send_data[0],
                    &send_n[0],
                    &send_displacement[0],
                    MPI_UNSIGNED,
                    &Haloed_eqns[0],
                    &Haloed_n[0],
                    &Haloed_displacement[0],
                    MPI_UNSIGNED,
                    dist_pt->communicator_pt()->mpi_comm());

      // Finally, we translate the map of halo entries into the permanent
      // storage
      Halo_n.resize(n_proc, 0);
      Halo_displacement.resize(n_proc, 0);

      // Loop over all the entries in the map
      unsigned receive_haloed_count = 0;
      for (int d = 0; d < n_proc; d++)
      {
        // Pointer to the map entry
        std::map<unsigned, Vector<unsigned>>::iterator it =
          halo_entries.find(d);
        // If we don't have it in the map, skip
        if (it == halo_entries.end())
        {
          Halo_displacement[d] = receive_haloed_count;
          Halo_n[d] = 0;
        }
        else
        {
          Halo_displacement[d] = receive_haloed_count;
          const int size_ = it->second.size();
          Halo_n[d] = size_;
          // Resize the equations to be sent
          Halo_eqns.resize(receive_haloed_count + size_);
          for (int i = 0; i < size_; i++)
          {
            Halo_eqns[receive_haloed_count + i] = it->second[i];
          }
          receive_haloed_count += size_;
        }
      }
    }
#endif
  }

  //=====================================================================
  /// Function that sets up a vector of pointers to halo
  /// data, index using the scheme in Local_index. The first arguement
  /// is a map of pointers to all halo data index by the global equation
  /// number
  //====================================================================
  void DoubleVectorHaloScheme::setup_halo_dofs(
    const std::map<unsigned, double*>& halo_data_pt,
    Vector<double*>& halo_dof_pt)
  {
    // How many entries are there in the map
    unsigned n_halo = Local_index.size();
    // Resize the vector
    halo_dof_pt.resize(n_halo);

    // Loop over all the entries in the map
    for (std::map<unsigned, unsigned>::iterator it = Local_index.begin();
         it != Local_index.end();
         ++it)
    {
      // Find the pointer in the halo_data_pt map
      std::map<unsigned, double*>::const_iterator it2 =
        halo_data_pt.find(it->first);
      // Did we find it
      if (it2 != halo_data_pt.end())
      {
        // Now set the entry
        halo_dof_pt[it->second] = it2->second;
      }
      else
      {
        std::ostringstream error_stream;
        error_stream << "Global equation " << it->first
                     << " reqired as halo is not stored in halo_data_pt\n";

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
  }

  //------------------------------------------------------------------
  // Member functions for the DoubleVectorWithHaloEntries
  //-------------------------------------------------------------------


  //=========================================================================
  /// Synchronise the halo data within the vector. This requires one
  /// "all to all" communnication.
  //====================================================================
  void DoubleVectorWithHaloEntries::synchronise()
  {
#ifdef OOMPH_HAS_MPI
    // Only need to do anything if the DoubleVector is distributed
    if (this->distributed())
    {
      // Read out the number of entries to send
      const unsigned n_send = Halo_scheme_pt->Haloed_eqns.size();
      Vector<double> send_data(n_send);
      // Read out the data values
      for (unsigned i = 0; i < n_send; i++)
      {
        send_data[i] = (*this)[Halo_scheme_pt->Haloed_eqns[i]];
      }

      // Read out the number of entries to receive
      const unsigned n_receive = Halo_scheme_pt->Halo_eqns.size();
      Vector<double> receive_data(n_receive);

      // Make sure that the send and receive data have size at least one
      if (n_send == 0)
      {
        send_data.resize(1);
      }
      if (n_receive == 0)
      {
        receive_data.resize(1);
      }
      // Communicate
      MPI_Alltoallv(&send_data[0],
                    &Halo_scheme_pt->Haloed_n[0],
                    &Halo_scheme_pt->Haloed_displacement[0],
                    MPI_DOUBLE,
                    &receive_data[0],
                    &Halo_scheme_pt->Halo_n[0],
                    &Halo_scheme_pt->Halo_displacement[0],
                    MPI_DOUBLE,
                    this->distribution_pt()->communicator_pt()->mpi_comm());


      // Now I need simply to update my local values
      for (unsigned i = 0; i < n_receive; i++)
      {
        Halo_value[Halo_scheme_pt->Halo_eqns[i]] = receive_data[i];
      }
    }
#endif
  }

  //=========================================================================
  /// Gather all ther data from multiple processors and sum the result
  /// which will be stored in the master copy and then synchronised to
  /// all copies. This requires two "all to all" communications
  //====================================================================
  void DoubleVectorWithHaloEntries::sum_all_halo_and_haloed_values()
  {
#ifdef OOMPH_HAS_MPI
    // Only need to do anything if the DoubleVector is distributed
    if (this->distributed())
    {
      // Send the Halo entries to the master processor
      const unsigned n_send = Halo_scheme_pt->Halo_eqns.size();
      Vector<double> send_data(n_send);
      // Read out the data values
      for (unsigned i = 0; i < n_send; i++)
      {
        send_data[i] = Halo_value[Halo_scheme_pt->Halo_eqns[i]];
      }

      // Read out the number of entries to receive
      const unsigned n_receive = Halo_scheme_pt->Haloed_eqns.size();
      Vector<double> receive_data(n_receive);

      // Make sure that the send and receive data have size at least one
      if (n_send == 0)
      {
        send_data.resize(1);
      }
      if (n_receive == 0)
      {
        receive_data.resize(1);
      }
      // Communicate
      MPI_Alltoallv(&send_data[0],
                    &Halo_scheme_pt->Halo_n[0],
                    &Halo_scheme_pt->Halo_displacement[0],
                    MPI_DOUBLE,
                    &receive_data[0],
                    &Halo_scheme_pt->Haloed_n[0],
                    &Halo_scheme_pt->Haloed_displacement[0],
                    MPI_DOUBLE,
                    this->distribution_pt()->communicator_pt()->mpi_comm());


      // Now I need simply to update and sum my  local values
      for (unsigned i = 0; i < n_receive; i++)
      {
        (*this)[Halo_scheme_pt->Haloed_eqns[i]] += receive_data[i];
      }

      // Then synchronise
      this->synchronise();
    }
#endif
  }


  //===================================================================
  /// Construct the halo scheme and storage for the halo data
  //=====================================================================
  void DoubleVectorWithHaloEntries::build_halo_scheme(
    DoubleVectorHaloScheme* const& halo_scheme_pt)
  {
    Halo_scheme_pt = halo_scheme_pt;

    if (Halo_scheme_pt != 0)
    {
      // Need to set up the halo data
      unsigned n_halo_data = halo_scheme_pt->Local_index.size();

      // Resize the halo storage
      Halo_value.resize(n_halo_data);

      // Now let's get the initial values from the other processors
      this->synchronise();
    }
  }


} // namespace oomph
