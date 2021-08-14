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
#include "macro_element_node_update_element.h"

namespace oomph
{
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // MacroElementNodeUpdateNodes
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Excute the update function: Update the current (and if
  /// update_all_time_levels_for_new_node==true also the previous)
  /// nodal position. Also update the current nodal values if
  /// an auxiliary update function is defined.
  /// Note: Updating of previous positions is only required (and should
  /// only be performed for) newly created MacroElementNodeUpdateNodes
  /// i.e. when this function is called from...
  /// If a node is hanging, its position is determined via its hanging
  /// node constraints after updating the position of its master nodes.
  //========================================================================
  void MacroElementNodeUpdateNode::node_update(
    const bool& update_all_time_levels_for_new_node)
  {
    // Number of time levels that need to be updated
    unsigned ntime;
    if (update_all_time_levels_for_new_node)
    {
      // Present value plus the previous values
      ntime = 1 + Position_time_stepper_pt->nprev_values();
    }
    else
    {
      // Just the present value
      ntime = 1;
    }

    // Is it a hanging node?
    if (is_hanging())
    {
      // Loop over all master nodes and update their position
      // That's all we need to update the position of hanging nodes!
      // (Recall that for hanging nodes Node::x(...) is not
      // guaranteed to be kept up-to-date; the (constrained!) nodal
      // position of hanging nodes must be determined via
      // Node::position() which determines the position
      // via the hanging node constraints from the position of
      // the master nodes)
      unsigned nmaster = hanging_pt()->nmaster();
      for (unsigned imaster = 0; imaster < nmaster; imaster++)
      {
        dynamic_cast<MacroElementNodeUpdateNode*>(
          hanging_pt()->master_node_pt(imaster))
          ->node_update(update_all_time_levels_for_new_node);
      }
    }
    // Node isn't hanging --> update it directly
    else
    {
      // If no update element is defined, keep the nodal positions where
      // they were (i.e. don't do anything), else update
      if (Node_update_element_pt != 0)
      {
        // Vector of local coordinates
        unsigned n_dim = ndim();
        Vector<double> x_new(n_dim);

        // Loop over time levels
        for (unsigned t = 0; t < ntime; t++)
        {
          // Update via macro element representation
          Node_update_element_pt->get_x(t, S_in_node_update_element, x_new);
          for (unsigned i = 0; i < n_dim; i++)
          {
            x(t, i) = x_new[i];
          }
        }
      }
    }

    // Perform auxiliary update of function values? Node passes itself
    // to this function so its position etc. is available to the auxiliary
    // node update function.
    if (Aux_node_update_fct_pt != 0)
    {
      Aux_node_update_fct_pt(this);
    }
  }

} // namespace oomph
