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
// Header file for a relatively simple Quad Meshe
#ifndef OOMPH_BACKWARD_STEP_MESH_HEADER
#define OOMPH_BACKWARD_STEP_MESH_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// OOMPH-LIB headers
#include "rectangular_quadmesh.h"
#include "generic/refineable_quad_mesh.h"

namespace oomph
{
  //=================================================================
  /// Backward step mesh
  //=================================================================
  template<class ELEMENT>
  class BackwardStepQuadMesh : public virtual RectangularQuadMesh<ELEMENT>
  {
  public:
    /// Pass overall number of elements in the horizontal
    /// and vertical directions, nx and ny, and the corresponding
    /// dimensions, lx and ly. nx_cut_out and ny_cut_out elements
    /// are cut out from the lower right corner to create the
    /// (reversed) backward step geometry. Timestepper defaults
    /// to Steady.
    BackwardStepQuadMesh(
      const unsigned& nx,
      const unsigned& ny,
      const unsigned& nx_cut_out,
      const unsigned& ny_cut_out,
      const double& lx,
      const double& ly,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : RectangularQuadMesh<ELEMENT>(nx, ny, lx, ly, time_stepper_pt)
    {
      // Do the actual build
      build_mesh(nx, ny, nx_cut_out, ny_cut_out, lx, ly);
    }

    /// Destructor: Empty
    virtual ~BackwardStepQuadMesh() {}

  private:
    /// Actual build function
    void build_mesh(const unsigned& nx,
                    const unsigned& ny,
                    const unsigned& nx_cut_out,
                    const unsigned& ny_cut_out,
                    const double& lx,
                    const double& ly);

  }; // end of mesh


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  //=================================================================
  /// Refineable backward step mesh
  //=================================================================
  template<class ELEMENT>
  class RefineableBackwardStepQuadMesh
    : public virtual BackwardStepQuadMesh<ELEMENT>,
      public RefineableQuadMesh<ELEMENT>
  {
  public:
    /// Pass overall number of elements in the horizontal
    /// and vertical directions, nx and ny, and the corresponding
    /// dimensions, lx and ly. nx_cut_out and ny_cut_out elements
    /// are cut out from the lower right corner to create the
    /// (reversed) backward step geometry. Timestepper defaults
    /// to Steady.
    RefineableBackwardStepQuadMesh(
      const unsigned& nx,
      const unsigned& ny,
      const unsigned& nx_cut_out,
      const unsigned& ny_cut_out,
      const double& lx,
      const double& ly,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : RectangularQuadMesh<ELEMENT>(nx, ny, lx, ly, time_stepper_pt),
        BackwardStepQuadMesh<ELEMENT>(
          nx, ny, nx_cut_out, ny_cut_out, lx, ly, time_stepper_pt)
    {
      // Nodal positions etc. were created in constructor for
      // SimpleRectangularQuadMesh<...> --> We only need to set up
      // adaptivity information: Associate finite elements with their
      // QuadTrees and plant them in a QuadTreeForest:
      this->setup_quadtree_forest();

    } // end of constructor

    /// Destructor: Empty
    virtual ~RefineableBackwardStepQuadMesh() {}

  }; // end of mesh

} // namespace oomph
#include "backward_step_mesh.tpp"
#endif
