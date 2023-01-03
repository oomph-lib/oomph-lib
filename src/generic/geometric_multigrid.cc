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
// Include the corresponding header file
#include "geometric_multigrid.h"

namespace oomph
{
  //=========================================================================
  /// Given the son type of the element and the local coordinate s of
  /// a given node in the son element, return the local coordinate s in its
  /// father element. 2D case.
  //=========================================================================
  template<>
  void MGSolver<2>::level_up_local_coord_of_node(const int& son_type,
                                                 Vector<double>& s)
  {
    // If the element is unrefined between the levels the local coordinate
    // of the node in one element is the same as that in the other element
    // therefore we only need to perform calculations if the levels are
    // different (i.e. son_type is not OMEGA)
    if (son_type != Tree::OMEGA)
    {
      // Scale the local coordinate from the range [-1,1]x[-1,1] to the range
      // [0,1]x[0,1] to match the width of the local coordinate range of the
      // fine element from the perspective of the father element. This
      // then simply requires a shift of the coordinates to match which type
      // of son element we're dealing with
      s[0] = (s[0] + 1.0) / 2.0;
      s[1] = (s[1] + 1.0) / 2.0;

      // Cases: The son_type determines how the local coordinates should be
      // shifted to give the local coordinates in the coarse mesh element
      switch (son_type)
      {
          // If we're dealing with the bottom-left element we need to shift
          // the range [0,1]x[0,1] to [-1,0]x[-1,0]
        case QuadTreeNames::SW:
          s[0] -= 1;
          s[1] -= 1;
          break;

          // If we're dealing with the bottom-right element we need to shift
          // the range [0,1]x[0,1] to [0,1]x[-1,0]
        case QuadTreeNames::SE:
          s[1] -= 1;
          break;

          // If we're dealing with the top-right element we need to shift the
          // range [0,1]x[0,1] to [0,1]x[0,1], i.e. nothing needs to be done
        case QuadTreeNames::NE:
          break;

          // If we're dealing with the top-left element we need to shift
          // the range [0,1]x[0,1] to [-1,0]x[0,1]
        case QuadTreeNames::NW:
          s[0] -= 1;
          break;
      }
    } // if son_type!=Tree::OMEGA
  } // End of level_up_local_coord_of_node

  //=========================================================================
  /// Given the son type of the element and the local coordinate s of
  /// a given node in the son element, return the local coordinate s in its
  /// father element. 3D case.
  //=========================================================================
  template<>
  void MGSolver<3>::level_up_local_coord_of_node(const int& son_type,
                                                 Vector<double>& s)
  {
    // If the element is unrefined between the levels the local coordinate
    // of the node in one element is the same as that in the other element
    // therefore we only need to perform calculations if the levels are
    // different (i.e. son_type is not OMEGA)
    if (son_type != Tree::OMEGA)
    {
      // Scale the local coordinate from the range [-1,1]x[-1,1]x[-1,1]
      // to the range [0,1]x[0,1]x[0,1] to match the width of the local
      // coordinate range of the fine element from the perspective of
      // the father element. This then simply requires a shift of the
      // coordinates to match which type of son element we're dealing with
      s[0] = (s[0] + 1.0) / 2.0;
      s[1] = (s[1] + 1.0) / 2.0;
      s[2] = (s[2] + 1.0) / 2.0;

      // Cases: The son_type determines how the local coordinates should be
      // shifted to give the local coordinates in the coarse mesh element
      switch (son_type)
      {
        case OcTreeNames::LDF:
          s[0] -= 1;
          s[1] -= 1;
          break;

        case OcTreeNames::LDB:
          s[0] -= 1;
          s[1] -= 1;
          s[2] -= 1;
          break;

        case OcTreeNames::LUF:
          s[0] -= 1;
          break;

        case OcTreeNames::LUB:
          s[0] -= 1;
          s[2] -= 1;
          break;

        case OcTreeNames::RDF:
          s[1] -= 1;
          break;

        case OcTreeNames::RDB:
          s[1] -= 1;
          s[2] -= 1;
          break;

        case OcTreeNames::RUF:
          break;

        case OcTreeNames::RUB:
          s[2] -= 1;
          break;
      }
    } // if son_type!=Tree::OMEGA
  } // End of level_up_local_coord_of_node
} // End of namespace oomph
