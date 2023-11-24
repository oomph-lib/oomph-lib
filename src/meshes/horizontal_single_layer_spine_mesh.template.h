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
#ifndef OOMPH_HORIZONTAL_SINGLE_LAYER_SPINE_MESH_HEADER
#define OOMPH_HORIZONTAL_SINGLE_LAYER_SPINE_MESH_HEADER

// oomph-lib includes
#include "../generic/spines.h"
#include "rectangular_quadmesh.template.h"

// Created by Francisco

namespace oomph
{
  //======================================================================
  /// Horizontal Single-layer spine mesh class derived from standard 2D mesh.
  /// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
  /// e.g  SpineElement<QCrouzeixRaviartElement<2>)
  /// and the information required to update their position. Additional
  /// equations must be specified in order to determine how the spines move.
  //======================================================================
  template<class ELEMENT>
  class HorizontalSingleLayerSpineMesh : public RectangularQuadMesh<ELEMENT>,
                                         public SpineMesh
  {
  public:
    /// Constructor: Pass number of elements in x-direction, number of
    /// elements in y-direction, axial length, height of layer, and pointer
    /// to timestepper (defaults to Steady timestepper)
    HorizontalSingleLayerSpineMesh(
      const unsigned& nx,
      const unsigned& ny,
      const double& lx,
      const double& h,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);


    /// General node update function implements pure virtual function
    /// defined in SpineMesh base class and performs specific node update
    /// actions:  along vertical spines
    virtual void spine_node_update(SpineNode* spine_node_pt)
    {
      // Get fraction along the spine
      double W = spine_node_pt->fraction();
      // Get spine height
      double H = spine_node_pt->h();
      // Set the value of y
      // spine_node_pt->x(0) = this->Xmin + W*H;
      spine_node_pt->x(0) = W * H;
    }


  protected:
    /// Helper function to actually build the single-layer spine mesh
    /// (called from various constructors)
    virtual void build_horizontal_single_layer_mesh(
      TimeStepper* time_stepper_pt);
  };

} // namespace oomph

#endif
