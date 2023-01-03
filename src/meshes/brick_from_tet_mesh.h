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
#ifndef OOMPH_BRICK_FROM_TET_MESH_HEADER
#define OOMPH_BRICK_FROM_TET_MESH_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
// mpi headers
#include "mpi.h"
#endif

#include <iterator>
#include <algorithm>

#include "generic/mesh.h"
#include "generic/tet_mesh.h"
#include "generic/brick_mesh.h"
#include "generic/refineable_brick_mesh.h"
#include "generic/Telements.h"
#include "xda_tet_mesh.h"
#include "tetgen_mesh.h"

namespace oomph
{
  //=====================================================================
  /// Brick mesh built by brickifying an existing tet mesh -- each
  /// tet gets split into four bricks. Can only be built with
  /// quadratic (27 node) elements.
  //=====================================================================
  template<class ELEMENT>
  class BrickFromTetMesh : public virtual BrickMeshBase
  {
  public:
    /// Constructor: Pass xda file name.
    BrickFromTetMesh(const std::string xda_file_name,
                     TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    {
      // Mesh can only be built with 3D Qelements.
      MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3,
                                                                            3);

      // Build temporary tet mesh
      XdaTetMesh<TElement<3, 3>>* tmp_mesh_pt =
        new XdaTetMesh<TElement<3, 3>>(xda_file_name, time_stepper_pt);

      // Actually build the mesh
      build_mesh(tmp_mesh_pt, time_stepper_pt);

      // Now kill the temporary mesh
      delete tmp_mesh_pt;
    }

    /// Constructor: Pass the files required for the tetgen mesh.
    BrickFromTetMesh(const std::string& node_file_name,
                     const std::string& element_file_name,
                     const std::string& face_file_name,
                     const bool& split_corner_elements,
                     TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
                     const bool& use_attributes = false)
    {
      // Mesh can only be built with 3D Qelements.
      MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3,
                                                                            3);

      // Build temporary tet mesh
      TetgenMesh<TElement<3, 3>>* tmp_mesh_pt =
        new TetgenMesh<TElement<3, 3>>(node_file_name,
                                       element_file_name,
                                       face_file_name,
                                       split_corner_elements,
                                       time_stepper_pt,
                                       use_attributes);

      // Actually build the mesh
      build_mesh(tmp_mesh_pt, time_stepper_pt);

      // Now kill the temporary mesh
      delete tmp_mesh_pt;
    }


    /// Constructor: Pass xda file name. This returns a pointer to the
    /// internally built XdaTetMesh for external use. Note that YOU
    /// are responsible for deleting this mesh.
    BrickFromTetMesh(const std::string xda_file_name,
                     XdaTetMesh<TElement<3, 3>>*& xda_tet_mesh_pt,
                     TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    {
      // Mesh can only be built with 3D Qelements.
      MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3,
                                                                            3);

      // Build temporary tet mesh
      xda_tet_mesh_pt =
        new XdaTetMesh<TElement<3, 3>>(xda_file_name, time_stepper_pt);

      // Actually build the brick mesh
      build_mesh(xda_tet_mesh_pt, time_stepper_pt);

      // Note that we're keeping the tet mesh alive for external use...
    }

    /// Access functions to the Vector of oomph-lib boundary ids
    /// that make up boundary b in the original xda enumeration
    Vector<unsigned> oomph_lib_boundary_ids(const unsigned& xda_boundary_id)
    {
      return Boundary_id[xda_boundary_id];
    }

  private:
    /// Build fct: Pass pointer to existing tet mesh.
    void build_mesh(XdaTetMesh<TElement<3, 3>>* tet_mesh_pt,
                    TimeStepper* time_stepper_pt);

    /// Build fct: Pass pointer to existing tet mesh.
    void build_mesh(TetgenMesh<TElement<3, 3>>* tet_mesh_pt,
                    TimeStepper* time_stepper_pt);

    /// Vector of vectors containing the boundary IDs of
    /// the overall boundary specified in the xda file.
    Vector<Vector<unsigned>> Boundary_id;
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  //=====================================================================
  /// Solid brick mesh built by brickifying an existing tet mesh -- each
  /// tet gets split into four bricks. Can only be built with
  /// quadratic (27 node) elements.
  //=====================================================================
  template<class ELEMENT>
  class SolidBrickFromTetMesh : public virtual BrickFromTetMesh<ELEMENT>,
                                public SolidMesh
  {
  public:
    /// Constructor: Pass xda file name.
    SolidBrickFromTetMesh(
      const std::string xda_file_name,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : BrickFromTetMesh<ELEMENT>(xda_file_name, time_stepper_pt)
    {
      // Make the current configuration the undeformed one by
      // setting the nodal Lagrangian coordinates to their current
      // Eulerian ones
      set_lagrangian_nodal_coordinates();
    }


    /// Constructor: Pass xda file name. This returns a pointer to the
    /// internally built XdaTetMesh for external use. Note that YOU
    /// are responsible for deleting this mesh.
    SolidBrickFromTetMesh(
      const std::string xda_file_name,
      XdaTetMesh<TElement<3, 3>>*& xda_tet_mesh_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : BrickFromTetMesh<ELEMENT>(
          xda_file_name, xda_tet_mesh_pt, time_stepper_pt)
    {
      // Make the current configuration the undeformed one by
      // setting the nodal Lagrangian coordinates to their current
      // Eulerian ones
      set_lagrangian_nodal_coordinates();
    }
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  //=====================================================================
  /// Refineable brick mesh built by brickifying an existing tet mesh -- each
  /// tet gets split into four bricks. Can only be built with
  /// quadratic (27 node) elements.
  //=====================================================================
  template<class ELEMENT>
  class RefineableBrickFromTetMesh : public virtual BrickFromTetMesh<ELEMENT>,
                                     public virtual RefineableBrickMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass xda file name.
    RefineableBrickFromTetMesh(
      const std::string xda_file_name,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : BrickFromTetMesh<ELEMENT>(xda_file_name, time_stepper_pt)
    {
      // Nodal positions etc. were created in constructor for
      // nonrefineable mesh. Only need to setup quadtree forest
      this->setup_octree_forest();
    }


    /// Constructor: Pass xda file name. This returns a pointer to the
    /// internally built XdaTetMesh for external use. Note that YOU
    /// are responsible for deleting this mesh.
    RefineableBrickFromTetMesh(
      const std::string xda_file_name,
      XdaTetMesh<TElement<3, 3>>*& xda_tet_mesh_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : BrickFromTetMesh<ELEMENT>(
          xda_file_name, xda_tet_mesh_pt, time_stepper_pt)
    {
      // Nodal positions etc. were created in constructor for
      // RectangularMesh<...>. Only need to setup quadtree forest
      this->setup_octree_forest();
    }
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  //=====================================================================
  /// Refineable solid brick mesh built by brickifying an existing tet
  /// mesh -- each tet gets split into four bricks. Can only be built with
  /// quadratic (27 node) elements.
  //=====================================================================
  template<class ELEMENT>
  class RefineableSolidBrickFromTetMesh
    : public virtual BrickFromTetMesh<ELEMENT>,
      public virtual RefineableBrickMesh<ELEMENT>,
      public SolidMesh

  {
  public:
    /// Constructor: Pass xda file name.
    RefineableSolidBrickFromTetMesh(
      const std::string xda_file_name,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : BrickFromTetMesh<ELEMENT>(xda_file_name, time_stepper_pt)
    {
      // Make the current configuration the undeformed one by
      // setting the nodal Lagrangian coordinates to their current
      // Eulerian ones
      set_lagrangian_nodal_coordinates();

      // Nodal positions etc. were created in constructor for
      // nonrefineable mesh. Only need to setup quadtree forest
      this->setup_octree_forest();
    }


    /// Constructor: Pass xda file name. This returns a pointer to the
    /// internally built XdaTetMesh for external use. Note that YOU
    /// are responsible for deleting this mesh.
    RefineableSolidBrickFromTetMesh(
      const std::string xda_file_name,
      XdaTetMesh<TElement<3, 3>>*& xda_tet_mesh_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : BrickFromTetMesh<ELEMENT>(
          xda_file_name, xda_tet_mesh_pt, time_stepper_pt)
    {
      // Make the current configuration the undeformed one by
      // setting the nodal Lagrangian coordinates to their current
      // Eulerian ones
      set_lagrangian_nodal_coordinates();

      // Nodal positions etc. were created in constructor for
      // RectangularMesh<...>. Only need to setup quadtree forest
      this->setup_octree_forest();
    }
  };

} // namespace oomph

#include "brick_from_tet_mesh.tpp"
#endif
