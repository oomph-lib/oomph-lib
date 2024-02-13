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
// Non-inline functions for BlockPrecSpaceTimeNavierStokes elements
#include "spacetime_navier_stokes_block_preconditionable_elements.h"

/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////

namespace oomph
{
  //=====start_of_setup=========================================================
  /// The setup function...
  //============================================================================
  void BlockPrecQTaylorHoodSpaceTimeElement::get_dof_numbers_for_unknowns(
    std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // Number of nodes
    unsigned n_node = this->nnode();

    // Temporary pair (used to store dof lookup prior to being added to list)
    std::pair<unsigned, unsigned> dof_lookup;

    // The number of nodes in the each direction
    unsigned n_node_1d = this->nnode_1d();

    // Make sure that we're using quadratic interpolation otherwise this might
    // go wrong...
    if (n_node_1d != 3)
    {
      // Create an output stream
      std::ostringstream error_message_stream;

      // Create an error message
      error_message_stream
        << "Can only deal with Navier-Stokes elements which "
        << "use quadratic interpolation at the moment. Using " << n_node_1d
        << " nodes in each direction" << std::endl;

      // Throw an error
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Loop over the nodes
    for (unsigned j = 0; j < n_node; j++)
    {
      // Only bother if the node isn't a copy...
      if (!(this->node_pt(j)->is_a_copy()))
      {
        // Storage for the local time slice ID (0<i_local<n_node_1d)
        unsigned i_local = 0;

        // If we're on the first time slice
        if (j < n_node_1d * n_node_1d)
        {
          // In the first local time slice
          i_local = 0;
        }
        // If we're in the second time slice (already excluded the first time
        // slice)
        else if (j < 2 * n_node_1d * n_node_1d)
        {
          // In the second local time slice
          i_local = 1;
        }
        // If we're in the final time slice (already excluded the first/second
        // slice)
        else if (j < 3 * n_node_1d * n_node_1d)
        {
          // In the final local time slice
          i_local = 2;
        }
        // We should never get here
        else
        {
          // Create an output stream
          std::ostringstream error_message_stream;

          // Create an error message
          error_message_stream << "Looping over too many nodes!" << std::endl;

          // Throw an error
          throw OomphLibError(error_message_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        // The number of dofs at this node
        unsigned n_dof = this->node_pt(j)->nvalue();

        // Loop over components
        for (unsigned i = 0; i < n_dof; i++)
        {
          // Determine local eqn number
          int local_eqn_number = this->nodal_local_eqn(j, i);

          // Ignore pinned values - far away degrees of freedom resulting
          // from hanging nodes can be ignored since these are be dealt
          // with by the element containing their master nodes
          if (local_eqn_number >= 0)
          {
            // The id shift (depending on which elemental time slice we're in)
            unsigned elemental_id_shift = 5 * Time_slab_id;

            // Store dof lookup in temporary pair: Global equation number
            // is the first entry in pair
            dof_lookup.first = this->eqn_number(local_eqn_number);

            // If we're on the first time slice
            if (i_local == 0)
            {
              // Here i can range only between 0 and 2 (but it's also an
              // unsigned so we only need to check if it's strictly greater than
              // 2)
              if (i > 2)
              {
                // Create an output stream
                std::ostringstream error_message_stream;

                // Create an error message
                error_message_stream
                  << "Don't know what to do when i=" << i
                  << ". Can only handle when i is between 0 and 2!"
                  << std::endl;

                // Throw an error
                throw OomphLibError(error_message_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }

              // Set dof numbers: Dof number is the second entry in pair
              dof_lookup.second = elemental_id_shift + i;
            }
            // If we're on the second time slice
            else if (i_local == 1)
            {
              // Here i can range only between 0 and 1
              if (i > 1)
              {
                // Create an output stream
                std::ostringstream error_message_stream;

                // Create an error message
                error_message_stream
                  << "Don't know what to do when i=" << i
                  << ". Can only handle when i is either 0 or 1!" << std::endl;

                // Throw an error
                throw OomphLibError(error_message_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }

              // The local dof shift given which local time slice we're on (i.e.
              // which node we're on in the time-direction). We've already
              // covered 2 velocities and 1 pressure so the dof shift is 3.
              unsigned local_id_shift = 3;

              // Set dof numbers: Dof number is the second entry in pair
              dof_lookup.second = elemental_id_shift + local_id_shift + i;
            }
            // If we're on the final time slice
            else if (i_local == 2)
            {
              // Here i can range only between 0 and 2
              if (i > 2)
              {
                // Create an output stream
                std::ostringstream error_message_stream;

                // Create an error message
                error_message_stream
                  << "Don't know what to do when i=" << i
                  << ". Can only handle when i is between 0 and 2!"
                  << std::endl;

                // Throw an error
                throw OomphLibError(error_message_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }

              // The local dof shift given which local time slice we're on (i.e.
              // which node we're on in the time-direction). We've already
              // covered 4 velocities and 1 pressure so the dof shift is 5.
              unsigned local_id_shift = 5;

              // Set dof numbers: Dof number is the second entry in pair
              dof_lookup.second = elemental_id_shift + local_id_shift + i;
            } // if (i_local==0)

            // Add to list
            dof_lookup_list.push_front(dof_lookup);
          }
        } // for (unsigned i=0;i<DIM;i++)
      } // if (!(this->node_pt(j)->is_a_copy()))
    } // for (unsigned j=0;j<n_node;j++)
  } // End of get_dof_numbers_for_unknowns
} // End of namespace oomph
