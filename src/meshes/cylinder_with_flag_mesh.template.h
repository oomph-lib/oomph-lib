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
#ifndef OOMPH_CYLINDER_WITH_FLAG_MESH_HEADER
#define OOMPH_CYLINDER_WITH_FLAG_MESH_HEADER


// Include the headers file for domain
#include "cylinder_with_flag_domain.h"

// Generic includes
#include "../generic/refineable_quad_mesh.h"
#include "../generic/quad_mesh.h"

// Include algebraic elements
#include "../generic/algebraic_elements.h"


namespace oomph
{
  //=============================================================
  /// Domain-based mesh for cylinder with flag as in Turek
  /// benchmark.
  //=============================================================
  template<class ELEMENT>
  class CylinderWithFlagMesh : public virtual Mesh, public virtual QuadMeshBase
  {
  public:
    /// Constructor. Pass the pointers to the GeomObjects that
    /// parametrise the cylinder, the three edges of the flag, the length and
    /// height of the domain, the length and height of the flag, the coordinates
    /// of the centre of the cylinder and its radius. Timestepper defaults to
    /// Steady default timestepper.
    CylinderWithFlagMesh(
      Circle* cylinder_pt,
      GeomObject* top_flag_pt,
      GeomObject* bottom_flag_pt,
      GeomObject* tip_flag_pt,
      const double& length,
      const double& height,
      const double& flag_length,
      const double& flag_height,
      const double& centre_x,
      const double& centre_y,
      const double& a,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Destructor: Kill the domain
    virtual ~CylinderWithFlagMesh()
    {
      delete Domain_pt;
    }

    /// Access function to the domain
    CylinderWithFlagDomain* domain_pt()
    {
      return Domain_pt;
    }

  protected:
    /// Pointer to the domain
    CylinderWithFlagDomain* Domain_pt;

  }; // end of mesh class


  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


  //===================================================================
  /// Refineable version of CylinderWithFlagMesh.
  //===================================================================
  template<class ELEMENT>
  class RefineableCylinderWithFlagMesh : public CylinderWithFlagMesh<ELEMENT>,
                                         public RefineableQuadMesh<ELEMENT>
  {
  public:
    /// Constructor. Pass the pointers to the GeomObjects that
    /// parametrise the cylinder, the three edges of the flag, the length and
    /// height of the domain, the length and height of the flag, the coordinates
    /// of the centre of the cylinder and its radius. Timestepper defaults to
    /// Steady default timestepper.
    RefineableCylinderWithFlagMesh(
      Circle* cylinder_pt,
      GeomObject* top_flag_pt,
      GeomObject* bottom_flag_pt,
      GeomObject* tip_flag_pt,
      const double& length,
      const double& height,
      const double& flag_length,
      const double& flag_height,
      const double& centre_x,
      const double& centre_y,
      const double& a,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : CylinderWithFlagMesh<ELEMENT>(cylinder_pt,
                                      top_flag_pt,
                                      bottom_flag_pt,
                                      tip_flag_pt,
                                      length,
                                      height,
                                      flag_length,
                                      flag_height,
                                      centre_x,
                                      centre_y,
                                      a,
                                      time_stepper_pt)
    {
      // Nodal positions etc. were created in constructor for
      // Cylinder...<...>. Need to setup adaptive information.

      // Setup quadtree forest for mesh refinement
      this->setup_quadtree_forest();
    }


    /// Destructor: Empty
    virtual ~RefineableCylinderWithFlagMesh() {}
  };


  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


  //===================================================================
  /// Algebraic version of CylinderWithFlagMesh.
  //===================================================================
  template<class ELEMENT>
  class AlgebraicCylinderWithFlagMesh
    : public AlgebraicMesh,
      public virtual CylinderWithFlagMesh<ELEMENT>
  {
  public:
    /// Constructor. Pass the pointers to the GeomObjects that
    /// parametrise the cylinder, the three edges of the flag, the length and
    /// height of the domain, the length and height of the flag, the coordinates
    /// of the centre of the cylinder and its radius. Timestepper defaults to
    /// Steady default timestepper.
    AlgebraicCylinderWithFlagMesh(
      Circle* cylinder_pt,
      GeomObject* top_flag_pt,
      GeomObject* bottom_flag_pt,
      GeomObject* tip_flag_pt,
      const double& length,
      const double& height,
      const double& flag_length,
      const double& flag_height,
      const double& centre_x,
      const double& centre_y,
      const double& a,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : CylinderWithFlagMesh<ELEMENT>(cylinder_pt,
                                      top_flag_pt,
                                      bottom_flag_pt,
                                      tip_flag_pt,
                                      length,
                                      height,
                                      flag_length,
                                      flag_height,
                                      centre_x,
                                      centre_y,
                                      a,
                                      time_stepper_pt),
        Cylinder_pt(cylinder_pt),
        Top_flag_pt(top_flag_pt),
        Bottom_flag_pt(bottom_flag_pt),
        Tip_flag_pt(tip_flag_pt),
        Length(length),
        Height(height),
        Flag_length(flag_length),
        Flag_height(flag_height),
        Centre_x(centre_x),
        Centre_y(centre_y),
        A(a)
    {
      // Add the geometric objects to the list associated with this
      // AlgebraicMesh
      AlgebraicMesh::add_geom_object_list_pt(cylinder_pt);
      AlgebraicMesh::add_geom_object_list_pt(top_flag_pt);
      AlgebraicMesh::add_geom_object_list_pt(bottom_flag_pt);
      AlgebraicMesh::add_geom_object_list_pt(tip_flag_pt);

      // Setup algebraic node update operations
      setup_algebraic_node_update();
    }

    /// Destructor: empty
    virtual ~AlgebraicCylinderWithFlagMesh() {}


    /// Set geometric object that defines the
    /// bottom face of the flag
    void set_bottom_flag_pt(GeomObject* bottom_flag_pt)
    {
      // Need to alter the domain's bottom_flag_pt too
      this->domain_pt()->bottom_flag_pt() = bottom_flag_pt;
      Bottom_flag_pt = bottom_flag_pt;
    }

    /// Set the geometric object that defines the
    /// top face of the flag
    void set_top_flag_pt(GeomObject* top_flag_pt)
    {
      this->domain_pt()->top_flag_pt() = top_flag_pt;
      Top_flag_pt = top_flag_pt;
    }


    /// Set the geometric object that defines the
    /// tip of the flag
    void set_tip_flag_pt(GeomObject* tip_flag_pt)
    {
      this->domain_pt()->tip_flag_pt() = tip_flag_pt;
      Tip_flag_pt = tip_flag_pt;
    }


    /// Read-only access to geometric object that defines the
    /// bottom face of the flag
    GeomObject* bottom_flag_pt() const
    {
      return Bottom_flag_pt;
    }


    /// Read-only access to geometric object that defines the
    /// top face of the flag
    GeomObject* top_flag_pt() const
    {
      return Top_flag_pt;
    }

    /// Read-only access to geometric object that defines the
    /// tip of the flag
    GeomObject* tip_flag_pt() const
    {
      return Tip_flag_pt;
    }


    /// Update the geometric references that are used
    /// to update node after mesh adaptation.
    /// Empty -- no update of node update required without adaptativity
    void update_node_update(AlgebraicNode*& node_pt) {}


    /// Update nodal position at time level t (t=0: present;
    /// t>0: previous)
    void algebraic_node_update(const unsigned& t, AlgebraicNode*& node_pt);

  protected:
    /// Function to setup the algebraic node update
    void setup_algebraic_node_update();

    /// Helper function
    void node_update_I(const unsigned& t, AlgebraicNode*& node_pt);

    /// Helper function
    void node_update_II(const unsigned& t, AlgebraicNode*& node_pt);

    /// Helper function
    void node_update_III(const unsigned& t, AlgebraicNode*& node_pt);

    /// Helper function
    void node_update_IV(const unsigned& t, AlgebraicNode*& node_pt);

    /// Helper function
    void node_update_V(const unsigned& t, AlgebraicNode*& node_pt);

    /// Helper function
    void node_update_VI(const unsigned& t, AlgebraicNode*& node_pt);

    /// Helper function
    void node_update_VII(const unsigned& t, AlgebraicNode*& node_pt);

    /// Helper function
    void node_update_VIII(const unsigned& t, AlgebraicNode*& node_pt);

    /// Helper function
    void node_update_IX(const unsigned& t, AlgebraicNode*& node_pt);

    /// Cylinder
    GeomObject* Cylinder_pt;

    /// Top flag
    GeomObject* Top_flag_pt;

    /// Bottom flag
    GeomObject* Bottom_flag_pt;

    /// Tip flag
    GeomObject* Tip_flag_pt;

    /// Length of the domain
    double Length;

    /// Height of the domain
    double Height;

    /// Flag length
    double Flag_length;

    /// Flag thickness
    double Flag_height;

    /// x position of the centre of the cylinder
    double Centre_x;

    /// x position of the centre of the cylinder
    double Centre_y;

    /// radius of the cylinder
    double A;
  };


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////

  //===================================================================
  /// Refineable version of  AlgebraicCylinderWithFlagMesh
  //===================================================================
  template<class ELEMENT>
  class RefineableAlgebraicCylinderWithFlagMesh
    : public RefineableQuadMesh<ELEMENT>,
      public virtual AlgebraicCylinderWithFlagMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass the pointers to the GeomObjects that
    /// parametrise the cylinder, the three edges of the flag, the length and
    /// height of the domain, the length and height of the flag, the coordinates
    /// of the centre of the cylinder and its radius. Timestepper defaults to
    /// Steady default timestepper.
    RefineableAlgebraicCylinderWithFlagMesh(
      Circle* cylinder_pt,
      GeomObject* top_flag_pt,
      GeomObject* bottom_flag_pt,
      GeomObject* tip_flag_pt,
      const double& length,
      const double& height,
      const double& flag_length,
      const double& flag_height,
      const double& centre_x,
      const double& centre_y,
      const double& a,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : CylinderWithFlagMesh<ELEMENT>(cylinder_pt,
                                      top_flag_pt,
                                      bottom_flag_pt,
                                      tip_flag_pt,
                                      length,
                                      height,
                                      flag_length,
                                      flag_height,
                                      centre_x,
                                      centre_y,
                                      a,
                                      time_stepper_pt),
        AlgebraicCylinderWithFlagMesh<ELEMENT>(cylinder_pt,
                                               top_flag_pt,
                                               bottom_flag_pt,
                                               tip_flag_pt,
                                               length,
                                               height,
                                               flag_length,
                                               flag_height,
                                               centre_x,
                                               centre_y,
                                               a,
                                               time_stepper_pt)
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
