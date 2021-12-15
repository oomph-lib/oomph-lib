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
// Non-inline functions for triangle/tet NS elements

#include "generalised_newtonian_Taxisym_navier_stokes_elements.h"


namespace oomph
{
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //========================================================================
  /// Unpin all internal pressure dofs.
  //========================================================================
  void GeneralisedNewtonianAxisymmetricTCrouzeixRaviartElement::
    unpin_all_internal_pressure_dofs()
  {
    unsigned n_pres = this->npres_axi_nst();
    // loop over pressure dofs
    for (unsigned l = 0; l < n_pres; l++)
    {
      // unpin internal pressure
      this->internal_data_pt(P_axi_nst_internal_index)->unpin(l);
    }
  }


  //=========================================================================
  ///  Add to the set \c paired_load_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for all values (pressures, velocities) that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  void GeneralisedNewtonianAxisymmetricTCrouzeixRaviartElement::
    identify_load_data(std::set<std::pair<Data*, unsigned>>& paired_load_data)
  {
    // Find the index at which the velocity is stored
    unsigned u_index[3];
    for (unsigned i = 0; i < 3; i++)
    {
      u_index[i] = this->u_index_axi_nst(i);
    }

    // Loop over the nodes
    unsigned n_node = this->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      // Loop over the velocity components and add pointer to their data
      // and indices to the vectors
      for (unsigned i = 0; i < 3; i++)
      {
        paired_load_data.insert(std::make_pair(this->node_pt(n), u_index[i]));
      }
    }

    // Identify the pressure data
    identify_pressure_data(paired_load_data);
  }


  //=========================================================================
  ///  Add to the set \c paired_pressue_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for all pressures values that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  void GeneralisedNewtonianAxisymmetricTCrouzeixRaviartElement::
    identify_pressure_data(
      std::set<std::pair<Data*, unsigned>>& paired_pressure_data)
  {
    // Loop over the internal data
    unsigned n_internal = this->ninternal_data();
    for (unsigned l = 0; l < n_internal; l++)
    {
      unsigned nval = this->internal_data_pt(l)->nvalue();
      // Add internal data
      for (unsigned j = 0; j < nval; j++)
      {
        paired_pressure_data.insert(
          std::make_pair(this->internal_data_pt(l), j));
      }
    }
  }


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////

  // Set the data for the number of Variables at each node
  const unsigned
    GeneralisedNewtonianAxisymmetricTTaylorHoodElement::Initial_Nvalue[6] = {
      4, 4, 4, 3, 3, 3};

  // Set the data for the pressure conversion array
  const unsigned GeneralisedNewtonianAxisymmetricTTaylorHoodElement::Pconv[3] =
    {0, 1, 2};


  //========================================================================
  /// Unpin all pressure dofs, incl the mid-face/side ones where
  /// they have been allocated (e.g. in the refineable version of this
  /// element).
  //========================================================================
  void GeneralisedNewtonianAxisymmetricTTaylorHoodElement::
    unpin_all_nodal_pressure_dofs()
  {
    unsigned n_node = this->nnode();
    // loop over nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      if (this->node_pt(l)->nvalue() == 3 + 1)
      {
        // unpin pressure dof
        this->node_pt(l)->unpin(3);
      }
    }
  }

  //========================================================================
  /// Pin all nodal pressure dofs, incl the mid-face/side ones where
  /// they have been allocated (e.g. in the refineable version of this
  /// element).
  //========================================================================
  void GeneralisedNewtonianAxisymmetricTTaylorHoodElement::
    pin_all_nodal_pressure_dofs()
  {
    // Loop over all nodes and pin pressure
    unsigned n_node = this->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      if (this->node_pt(n)->nvalue() == 3 + 1)
      {
        this->node_pt(n)->pin(3);
      }
    }
  }

  //========================================================================
  /// Unpin the proper nodal pressure dofs which are not hanging.
  //========================================================================
  void GeneralisedNewtonianAxisymmetricTTaylorHoodElement::
    unpin_proper_nodal_pressure_dofs()
  {
    // Loop over all pressure nodes and unpin if they're not hanging
    unsigned n_pres = npres_axi_nst();
    for (unsigned l = 0; l < n_pres; l++)
    {
      Node* nod_pt = this->node_pt(Pconv[l]);
      if (!nod_pt->is_hanging(3))
      {
        nod_pt->unpin(3);
      }
    }
  }


  //=========================================================================
  ///  Add to the set \c paired_load_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for all values (pressures, velocities) that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  void GeneralisedNewtonianAxisymmetricTTaylorHoodElement::identify_load_data(
    std::set<std::pair<Data*, unsigned>>& paired_load_data)
  {
    // Loop over the nodes
    unsigned n_node = this->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      // Loop over the velocity components and add pointer to their data
      // and indices to the vectors
      for (unsigned i = 0; i < 3; i++)
      {
        paired_load_data.insert(std::make_pair(this->node_pt(n), i));
      }
    }

    // Add the pressure data
    identify_pressure_data(paired_load_data);
  }

  //=========================================================================
  ///  Add to the set \c paired_load_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for all values (pressures, velocities) that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  void GeneralisedNewtonianAxisymmetricTTaylorHoodElement::
    identify_pressure_data(
      std::set<std::pair<Data*, unsigned>>& paired_load_data)
  {
    // Loop over the pressure data
    unsigned n_pres = npres_axi_nst();
    for (unsigned l = 0; l < n_pres; l++)
    {
      // The DIMth entry in each nodal data is the pressure, which
      // affects the traction
      paired_load_data.insert(std::make_pair(this->node_pt(Pconv[l]), 3));
    }
  }

} // namespace oomph
