// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
// Non-inline functions for pml meshes and associated (generic) elements

#include "pml_meshes.h"

namespace oomph
{
  //===================================================================
  /// Namespace with functions that allow the construction of
  /// PML layers on axis aligned boundaries
  //===================================================================
  namespace TwoDimensionalPMLHelper
  {
    /// helper function for sorting the right boundary nodes
    bool sorter_right_boundary(Node* nod_i_pt, Node* nod_j_pt)
    {
      return (nod_i_pt->x(1) < nod_j_pt->x(1));
    }

    /// helper function for sorting the top boundary nodes
    bool sorter_top_boundary(Node* nod_i_pt, Node* nod_j_pt)
    {
      return (nod_i_pt->x(0) < nod_j_pt->x(0));
    }

    /// helper function for sorting the left boundary nodes
    bool sorter_left_boundary(Node* nod_i_pt, Node* nod_j_pt)
    {
      return (nod_i_pt->x(1) > nod_j_pt->x(1));
    }

    /// helper function for sorting the bottom boundary nodes
    bool sorter_bottom_boundary(Node* nod_i_pt, Node* nod_j_pt)
    {
      return (nod_i_pt->x(0) > nod_j_pt->x(0));
    }
  } // namespace TwoDimensionalPMLHelper
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////

} // namespace oomph
