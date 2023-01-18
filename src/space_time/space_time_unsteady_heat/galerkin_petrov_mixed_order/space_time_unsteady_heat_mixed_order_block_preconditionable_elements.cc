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
// Non-inline functions for BlockPrecSpaceTimeUnsteadyHeat elements
#include "space_time_unsteady_heat_mixed_order_block_preconditionable_elements.h"

/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////

namespace oomph
{
  //=====start_of_setup=========================================================
  /// The setup function...
  //============================================================================
  template<unsigned SPATIAL_DIM, unsigned NNODE_1D>
  void BlockPrecQUnsteadyHeatMixedOrderSpaceTimeElement<SPATIAL_DIM, NNODE_1D>::
    get_dof_numbers_for_unknowns(
      std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // Number of nodes
    unsigned n_node = this->nnode();

    // Temporary pair (used to store dof lookup prior to being added to list)
    std::pair<unsigned, unsigned> dof_lookup;

    // Loop over the nodes but ignore the nodes on the first local time-slice
    // because their dof number will either be set by elements from the
    // previous space-time slab (including on time-periodic meshes) or they'll
    // be pinned through an initial condition.
    for (unsigned i = NNODE_1D * NNODE_1D; i < n_node; i++)
    {
      // Storage for the local time slice ID (0<=i_temporal<=NNODE_1D-1)
      unsigned i_temporal = 0;

      // The spatial node number
      unsigned i_spatial = i % (NNODE_1D * NNODE_1D);

      // Which local time slice are we in?
      i_temporal = (i - i_spatial) / (NNODE_1D * NNODE_1D);

      // Find the index at which the variable is stored
      unsigned u_nodal_index = this->u_index_ust_heat();

      // Determine local eqn number
      int local_eqn_number = this->nodal_local_eqn(i, u_nodal_index);

      // Ignore pinned values - far away degrees of freedom resulting
      // from hanging nodes can be ignored since these are be dealt
      // with by the element containing their master nodes
      if (local_eqn_number >= 0)
      {
        // The only unpinned dofs we should encounter at this point should be
        // those on the final temporal boundary of the space-time element
        if (i_temporal != NNODE_1D - 1)
        {
          // Create an output stream
          std::ostringstream error_message_stream;

          // Create an error message
          error_message_stream << "All nodes strictly on the interior of the "
                               << "elements temporal boundaries\nshould be "
                               << "pinned in the mixed order element!"
                               << std::endl;

          // Throw an error
          throw OomphLibError(error_message_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        // Store dof lookup in temporary pair: Global equation number
        // is the first entry in pair
        dof_lookup.first = this->eqn_number(local_eqn_number);

        // Set dof numbers: dof number is the second entry in pair
        dof_lookup.second = (this->Time_slab_id);

        // Add to list
        dof_lookup_list.push_front(dof_lookup);
      } // if (local_eqn_number>=0)
    } // for (unsigned i=NNODE_1D*NNODE_1D; i<n_node; i++)
  } // End of get_dof_numbers_for_unknowns


  //=====start_of_setup=========================================================
  /// The setup function...
  //============================================================================
  template<unsigned SPATIAL_DIM, unsigned NNODE_1D>
  void BlockPrecRefineableQUnsteadyHeatMixedOrderSpaceTimeElement<SPATIAL_DIM,
                                                                  NNODE_1D>::
    get_dof_numbers_for_unknowns(
      std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // Number of nodes
    unsigned n_node = this->nnode();

    // Temporary pair (used to store dof lookup prior to being added to list)
    std::pair<unsigned, unsigned> dof_lookup;

    // Loop over the nodes but ignore the nodes on the first local time-slice
    // because their dof number will either be set by elements from the
    // previous space-time slab (including on time-periodic meshes) or they'll
    // be pinned through an initial condition.
    for (unsigned i = NNODE_1D * NNODE_1D; i < n_node; i++)
    {
      // Storage for the local time slice ID (0<=i_temporal<=NNODE_1D-1)
      unsigned i_temporal = 0;

      // The spatial node number
      unsigned i_spatial = i % (NNODE_1D * NNODE_1D);

      // Which local time slice are we in?
      i_temporal = (i - i_spatial) / (NNODE_1D * NNODE_1D);

      // Find the index at which the variable is stored
      unsigned u_nodal_index = this->u_index_ust_heat();

      // Determine local eqn number
      int local_eqn_number = this->nodal_local_eqn(i, u_nodal_index);

      // Ignore pinned values - far away degrees of freedom resulting
      // from hanging nodes can be ignored since these are be dealt
      // with by the element containing their master nodes
      if (local_eqn_number >= 0)
      {
        // The only unpinned dofs we should encounter at this point should be
        // those on the final temporal boundary of the space-time element
        if (i_temporal != NNODE_1D - 1)
        {
          // Create an output stream
          std::ostringstream error_message_stream;

          // Create an error message
          error_message_stream << "All nodes strictly on the interior of the "
                               << "elements temporal boundaries\nshould be "
                               << "pinned in the mixed order element!"
                               << std::endl;

          // Throw an error
          throw OomphLibError(error_message_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        // Store dof lookup in temporary pair: Global equation number
        // is the first entry in pair
        dof_lookup.first = this->eqn_number(local_eqn_number);

        // Set dof numbers: dof number is the second entry in pair
        dof_lookup.second = (this->Time_slab_id);

        // Add to list
        dof_lookup_list.push_front(dof_lookup);
      } // if (local_eqn_number>=0)
    } // for (unsigned i=NNODE_1D*NNODE_1D; i<n_node; i++)
  } // End of get_dof_numbers_for_unknowns


  //====================================================================
  // Force build of templates
  //====================================================================
  template class BlockPrecQUnsteadyHeatMixedOrderSpaceTimeElement<2, 2>;
  template class BlockPrecQUnsteadyHeatMixedOrderSpaceTimeElement<2, 3>;
  template class BlockPrecQUnsteadyHeatMixedOrderSpaceTimeElement<2, 4>;
  template class BlockPrecRefineableQUnsteadyHeatMixedOrderSpaceTimeElement<2,
                                                                            2>;
  template class BlockPrecRefineableQUnsteadyHeatMixedOrderSpaceTimeElement<2,
                                                                            3>;
  template class BlockPrecRefineableQUnsteadyHeatMixedOrderSpaceTimeElement<2,
                                                                            4>;
} // End of namespace oomph
