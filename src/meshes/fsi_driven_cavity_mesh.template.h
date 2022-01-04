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

// Include guards
#ifndef OOMPH_FSI_DRIVEN_CAVITY_MESH_HEADER
#define OOMPH_FSI_DRIVEN_CAVITY_MESH_HEADER

// Generic includes
#include "../generic/refineable_quad_mesh.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"
#include "../generic/quad_mesh.h"

// Mesh is based on simple_rectangular_quadmesh
#include "simple_rectangular_quadmesh.template.h"
#include "simple_rectangular_quadmesh.template.cc"

// Include algebraic elements
#include "../generic/algebraic_elements.h"


namespace oomph
{
  //========================================================================
  /// Mesh for W. Wall's FSI driven cavity problem.
  /// The mesh is derived from the \c SimpleRectangularQuadMesh
  /// so it's node and element numbering scheme is the same
  /// as in that mesh. Only the boundaries are numbered differently
  /// to allow the easy identification of the "collapsible" segment.
  /// Boundary coordinates are set up for all nodes
  /// located on boundary 3 (the collapsible segment).
  /// The curvilinear ("collapsible") segment is defined by
  /// a \c GeomObject.
  /// - Boundary 0 is the moving lid.
  /// - Boundary 1 is the gap above the moving lid on the right wall
  /// - Boundary 2 is the rigid part of the right wall
  /// - Boundary 3 is the moving (elastic) wall
  /// - Boundary 4 is the rigid part of the left wall
  /// - Boundary 5 is the gap above the moving lid on the left wall
  //========================================================================
  template<class ELEMENT>
  class FSIDrivenCavityMesh : public SimpleRectangularQuadMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass number of elements, number of elements,
    /// fractional height of the gap above the moving wall,
    /// pointer to GeomObject that defines the collapsible segment and pointer
    /// to TimeStepper (defaults to the default timestepper, Steady).
    FSIDrivenCavityMesh(
      const unsigned& nx,
      const unsigned& ny,
      const double& lx,
      const double& ly,
      const double& gap_fraction,
      GeomObject* wall_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

  protected:
    /// Number of elements in x direction
    unsigned Nx;

    /// Number of elements in y direction
    unsigned Ny;

    /// Fraction of the gap next to moving lid, relative to the
    /// height of the domain
    double Gap_fraction;

    /// Pointer to geometric object that represents the moving wall
    GeomObject* Wall_pt;
  };


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //===================================================================
  /// Refineable version of FSIDrivenCavityMesh.
  /// The mesh is derived from the \c SimpleRectangularQuadMesh
  /// so it's node and element numbering scheme is the same
  /// as in that mesh. Only the boundaries are numbered differently
  /// to allow the easy identification of the "collapsible" segment.
  ///  Boundary coordinates are set up for all nodes
  /// located on boundary 3 (the collapsible segment).
  /// The curvilinear ("collapsible") segment is defined by
  /// a \c GeomObject.
  /// - Boundary 0 is the moving lid.
  /// - Boundary 1 is the gap above the moving lid on the right wall
  /// - Boundary 2 is the rigid part of the right wall
  /// - Boundary 3 is the moving (elastic) wall
  /// - Boundary 4 is the rigid part of the left wall
  /// - Boundary 5 is the gap above the moving lid on the left wall
  //====================================================================
  template<class ELEMENT>
  class RefineableFSIDrivenCavityMesh
    : public virtual FSIDrivenCavityMesh<ELEMENT>,
      public RefineableQuadMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass number of elements, lengths, pointer to
    /// geometric object that describes the wall and timestepper
    RefineableFSIDrivenCavityMesh(
      const unsigned& nx,
      const unsigned& ny,
      const double& lx,
      const double& ly,
      const double& gap_fraction,
      GeomObject* wall_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : FSIDrivenCavityMesh<ELEMENT>(
          nx, ny, lx, ly, gap_fraction, wall_pt, time_stepper_pt)
    {
      // Build quadtree forest
      this->setup_quadtree_forest();
    }


    /// Destructor(empty)
    ~RefineableFSIDrivenCavityMesh() {}
  };


  //=================================================================
  /// / Alebraic node update version of FSIDrivenCavityMesh
  /// - Boundary 0 is the moving lid.
  /// - Boundary 1 is the gap above the moving lid on the right wall
  /// - Boundary 2 is the rigid part of the right wall
  /// - Boundary 3 is the moving (elastic) wall
  /// - Boundary 4 is the rigid part of the left wall
  /// - Boundary 5 is the gap above the moving lid on the left wall
  //=================================================================
  template<class ELEMENT>
  class AlgebraicFSIDrivenCavityMesh
    : public virtual FSIDrivenCavityMesh<ELEMENT>,
      public AlgebraicMesh
  {
  public:
    /// Constructor: Pass number of elements, lengths, pointer to
    /// GeomObject that defines the collapsible segment and pointer to
    /// TimeStepper (defaults to the default timestepper, Steady).
    AlgebraicFSIDrivenCavityMesh(
      const unsigned& nx,
      const unsigned& ny,
      const double& lx,
      const double& ly,
      const double& gap_fraction,
      GeomObject* wall_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : FSIDrivenCavityMesh<ELEMENT>(
          nx, ny, lx, ly, gap_fraction, wall_pt, time_stepper_pt)
    {
      // Add the geometric object to the list associated with this AlgebraicMesh
      AlgebraicMesh::add_geom_object_list_pt(wall_pt);

      // Setup algebraic node update operations
      setup_algebraic_node_update();
    }

    /// Destructor: empty
    virtual ~AlgebraicFSIDrivenCavityMesh() {}

    /// Update nodal position at time level t (t=0: present;
    /// t>0: previous)
    void algebraic_node_update(const unsigned& t, AlgebraicNode*& node_pt);

    /// Update the node-udate data after mesh adaptation.
    /// Empty -- no update of node update required as this is
    /// non-refineable mesh.
    void update_node_update(AlgebraicNode*& node_pt) {}

  protected:
    /// Function to setup the algebraic node update
    void setup_algebraic_node_update();
  };


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Refineable version algebraic FSIDrivenCavityMesh.
  /// - Boundary 0 is the moving lid.
  /// - Boundary 1 is the gap above the moving lid on the right wall
  /// - Boundary 2 is the rigid part of the right wall
  /// - Boundary 3 is the moving (elastic) wall
  /// - Boundary 4 is the rigid part of the left wall
  /// - Boundary 5 is the gap above the moving lid on the left wall
  //=================================================================
  template<class ELEMENT>
  class RefineableAlgebraicFSIDrivenCavityMesh
    : public RefineableQuadMesh<ELEMENT>,
      public virtual AlgebraicFSIDrivenCavityMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass number of elements, lengths, pointer to
    /// GeomObject that defines the collapsible segment and pointer to
    /// TimeStepper (defaults to the default timestepper, Steady).
    RefineableAlgebraicFSIDrivenCavityMesh(
      const unsigned& nx,
      const unsigned& ny,
      const double& lx,
      const double& ly,
      const double& gap_fraction,
      GeomObject* wall_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : FSIDrivenCavityMesh<ELEMENT>(
          nx, ny, lx, ly, gap_fraction, wall_pt, time_stepper_pt),
        AlgebraicFSIDrivenCavityMesh<ELEMENT>(
          nx, ny, lx, ly, gap_fraction, wall_pt, time_stepper_pt)
    {
      // Build quadtree forest
      this->setup_quadtree_forest();
    }

    /// Update the node update data for specified node following
    /// any mesh adapation
    void update_node_update(AlgebraicNode*& node_pt);
  };


} // namespace oomph

#endif
