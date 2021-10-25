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
#ifndef OOMPH_THIN_LAYER_BRICK_ON_TET_MESH_HEADER
#define OOMPH_THIN_LAYER_BRICK_ON_TET_MESH_HEADER


#include "../generic/brick_mesh.h"
#include "../generic/refineable_brick_mesh.h"

namespace oomph
{
  //========================================================================
  /// Brick mesh layer built on top of a given tet mesh. Typically
  /// used in FSI problems where the tet mesh is the fluid mesh and this
  /// mesh acts as the solid mesh that surrounds the FSI interface.
  //========================================================================
  template<class ELEMENT>
  class ThinLayerBrickOnTetMesh : public virtual BrickMeshBase
  {
  public:
    /// Function pointer to function that specifies the wall thickness
    /// as a fct of the coordinates of the inner surface
    typedef void (*ThicknessFctPt)(const Vector<double>& x, double& h_thick);

    /// Constructor: Specify (quadratic) tet mesh, boundary IDs of
    /// boundary on which the current mesh is to be erected (in an FSI context
    /// this boundary tends to be the FSI boundary of the fluid mesh. Also
    /// specify the uniform thickness of layer, and the number of element
    /// layers. The vectors stored in in_out_boundary_ids contain the boundary
    /// IDs of the other boundaries in the tet mesh. In an FSI context
    /// these typically identify the in/outflow boundaries in the fluid
    /// mesh. The boundary enumeration of the current mesh follows the
    /// one of the underlying fluid mesh: The enumeration of the FSI boundary
    /// matches (to enable the setup of the FSI matching); the "in/outflow"
    /// faces in this mesh inherit the same enumeration as the in/outflow
    /// faces in the underlying fluid mesh. Finally, the "outer" boundary
    /// gets its own boundary ID.
    /// Timestepper defaults to steady pseudo-timestepper.
    ThinLayerBrickOnTetMesh(
      Mesh* tet_mesh_pt,
      const Vector<unsigned>& boundary_ids,
      ThicknessFctPt thickness_fct_pt,
      const unsigned& nlayer,
      const Vector<Vector<unsigned>>& in_out_boundary_id,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);


    /// Access functions to the Vector of oomph-lib boundary ids
    /// that make up boundary on which the mesh was erected (typically
    /// the FSI interface in an FSI problem)
    Vector<unsigned> fsi_boundary_id()
    {
      return FSI_boundary_id;
    }

    /// Boundary ID of the "outer" surface -- in an FSI context
    /// this is the non-wetted tube surface at a distance h_thick from
    /// the FSI surface
    unsigned outer_boundary_id()
    {
      return Outer_boundary_id;
    }

    /// Access function to the vector containing the ids of the oomph-lib
    /// mesh boundaries that make up the specified in/outflow boundaries
    /// as specified in constructor.
    Vector<unsigned> in_out_boundary_id(const unsigned& boundary_id)
    {
      return In_out_boundary_id[boundary_id];
    }


  private:
    /// Vector of oomph-lib boundary ids
    /// that make up boundary on which the mesh was erected (typically
    /// the FSI interface in an FSI problem)
    Vector<unsigned> FSI_boundary_id;

    /// Boundary ID of the "outer" surface -- the non-wetted
    /// tube surface at a distance h_thick from the FSI surface
    unsigned Outer_boundary_id;

    /// Vector of vectors containing the ids of the oomph-lib
    /// mesh boundaries that make up the specified in/outflow boundaries
    Vector<Vector<unsigned>> In_out_boundary_id;

    /// Function pointer to function that specifies the wall thickness
    /// as a fct of the coordinates of the inner surface
    ThicknessFctPt Thickness_fct_pt;
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Refineable brick mesh layer built on top of a given tet mesh. Typically
  /// used in FSI problems where the tet mesh is the fluid mesh and this
  /// mesh acts as the solid mesh that surrounds the FSI interface.
  //========================================================================
  template<class ELEMENT>
  class RefineableThinLayerBrickOnTetMesh
    : public virtual ThinLayerBrickOnTetMesh<ELEMENT>,
      public virtual RefineableBrickMesh<ELEMENT>
  {
  public:
    /// Function pointer to function that specifies the wall thickness
    /// as a fct of the coordinates of the inner surface
    typedef void (*ThicknessFctPt)(const Vector<double>& x, double& h_thick);

    /// Constructor: Specify (quadratic) tet mesh, boundary IDs of
    /// boundary on which the current mesh is to be erected (in an FSI context
    /// this boundary tends to be the FSI boundary of the fluid mesh. Also
    /// specify the uniform thickness of layer, and the number of element
    /// layers. The vectors stored in in_out_boundary_ids contain the boundary
    /// IDs of the other boundaries in the tet mesh. In an FSI context
    /// these typically identify the in/outflow boundaries in the fluid
    /// mesh. The boundary enumeration of the current mesh follows the
    /// one of the underlying fluid mesh: The enumeration of the FSI boundary
    /// matches (to enable the setup of the FSI matching); the "in/outflow"
    /// faces in this mesh inherit the same enumeration as the in/outflow
    /// faces in the underlying fluid mesh. Finally, the "outer" boundary
    /// gets its own boundary ID.
    /// Timestepper defaults to steady pseudo-timestepper.
    RefineableThinLayerBrickOnTetMesh(
      Mesh* tet_mesh_pt,
      const Vector<unsigned>& boundary_ids,
      ThicknessFctPt thickness_fct_pt,
      const unsigned& nlayer,
      const Vector<Vector<unsigned>>& in_out_boundary_id,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : ThinLayerBrickOnTetMesh<ELEMENT>(tet_mesh_pt,
                                         boundary_ids,
                                         thickness_fct_pt,
                                         nlayer,
                                         in_out_boundary_id,
                                         time_stepper_pt)
    {
      // Nodal positions etc. were created in constructor for
      // nonrefineable mesh. Only need to setup quadtree forest
      this->setup_octree_forest();
    }
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Solid brick mesh layer built on top of a given tet mesh. Typically
  /// used in FSI problems where the tet mesh is the fluid mesh and this
  /// mesh acts as the solid mesh that surrounds the FSI interface.
  //========================================================================
  template<class ELEMENT>
  class SolidThinLayerBrickOnTetMesh
    : public virtual ThinLayerBrickOnTetMesh<ELEMENT>,
      public virtual SolidMesh
  {
  public:
    /// Function pointer to function that specifies the wall thickness
    /// as a fct of the coordinates of the inner surface
    typedef void (*ThicknessFctPt)(const Vector<double>& x, double& h_thick);

    /// Constructor: Specify (quadratic) tet mesh, boundary IDs of
    /// boundary on which the current mesh is to be erected (in an FSI context
    /// this boundary tends to be the FSI boundary of the fluid mesh. Also
    /// specify the uniform thickness of layer, and the number of element
    /// layers. The vectors stored in in_out_boundary_ids contain the boundary
    /// IDs of the other boundaries in the tet mesh. In an FSI context
    /// these typically identify the in/outflow boundaries in the fluid
    /// mesh. The boundary enumeration of the current mesh follows the
    /// one of the underlying fluid mesh: The enumeration of the FSI boundary
    /// matches (to enable the setup of the FSI matching); the "in/outflow"
    /// faces in this mesh inherit the same enumeration as the in/outflow
    /// faces in the underlying fluid mesh. Finally, the "outer" boundary
    /// gets its own boundary ID.
    /// Timestepper defaults to steady pseudo-timestepper.
    SolidThinLayerBrickOnTetMesh(
      Mesh* tet_mesh_pt,
      const Vector<unsigned>& boundary_ids,
      ThicknessFctPt thickness_fct_pt,
      const unsigned& nlayer,
      const Vector<Vector<unsigned>>& in_out_boundary_id,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : ThinLayerBrickOnTetMesh<ELEMENT>(tet_mesh_pt,
                                         boundary_ids,
                                         thickness_fct_pt,
                                         nlayer,
                                         in_out_boundary_id,
                                         time_stepper_pt)
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


  //========================================================================
  /// Refineable solid brick mesh layer built on top of a given tet mesh.
  /// Typically used in FSI problems where the tet mesh is the fluid mesh and
  /// this mesh acts as the solid mesh that surrounds the FSI interface.
  //========================================================================
  template<class ELEMENT>
  class RefineableSolidThinLayerBrickOnTetMesh
    : public virtual ThinLayerBrickOnTetMesh<ELEMENT>,
      public virtual RefineableBrickMesh<ELEMENT>,
      public virtual SolidMesh
  {
  public:
    /// Function pointer to function that specifies the wall thickness
    /// as a fct of the coordinates of the inner surface
    typedef void (*ThicknessFctPt)(const Vector<double>& x, double& h_thick);

    /// Constructor: Specify (quadratic) tet mesh, boundary IDs of
    /// boundary on which the current mesh is to be erected (in an FSI context
    /// this boundary tends to be the FSI boundary of the fluid mesh. Also
    /// specify the uniform thickness of layer, and the number of element
    /// layers. The vectors stored in in_out_boundary_ids contain the boundary
    /// IDs of the other boundaries in the tet mesh. In an FSI context
    /// these typically identify the in/outflow boundaries in the fluid
    /// mesh. The boundary enumeration of the current mesh follows the
    /// one of the underlying fluid mesh: The enumeration of the FSI boundary
    /// matches (to enable the setup of the FSI matching); the "in/outflow"
    /// faces in this mesh inherit the same enumeration as the in/outflow
    /// faces in the underlying fluid mesh. Finally, the "outer" boundary
    /// gets its own boundary ID.
    /// Timestepper defaults to steady pseudo-timestepper.
    RefineableSolidThinLayerBrickOnTetMesh(
      Mesh* tet_mesh_pt,
      const Vector<unsigned>& boundary_ids,
      ThicknessFctPt thickness_fct_pt,
      const unsigned& nlayer,
      const Vector<Vector<unsigned>>& in_out_boundary_id,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : ThinLayerBrickOnTetMesh<ELEMENT>(tet_mesh_pt,
                                         boundary_ids,
                                         thickness_fct_pt,
                                         nlayer,
                                         in_out_boundary_id,
                                         time_stepper_pt)
    {
      // Make the current configuration the undeformed one by
      // setting the nodal Lagrangian coordinates to their current
      // Eulerian ones
      set_lagrangian_nodal_coordinates();

      // Nodal positions etc. were created in constructor for
      // nonrefineable mesh. Only need to setup quadtree forest
      this->setup_octree_forest();
    }
  };


} // namespace oomph

#endif
