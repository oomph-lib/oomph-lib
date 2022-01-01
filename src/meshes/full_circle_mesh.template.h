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

#ifndef OOMPH_FULL_CIRCLE_MESH_HEADER
#define OOMPH_FULL_CIRCLE_MESH_HEADER

// Headers
#include "../generic/refineable_quad_mesh.h"

// Include the headers file for domain
#include "full_circle_domain.h"

namespace oomph
{
  //====================================================================
  /// Full circle mesh class.
  /// The domain is specified by the GeomObject that identifies
  /// the entire area. Non-refineable base version!
  ///
  /// The mesh boundaries are numbered as follows:
  /// - Boundary 0: The outer wall, represented by \f$\xi_1 = 1\f$.
  /// .
  ///
  //====================================================================
  template<class ELEMENT>
  class FullCircleMesh : public virtual QuadMeshBase
  {
  public:
    /// Constructor: Pass pointer to geometric object that
    /// specifies the area; values of theta at which dividing lines
    /// are to be placed, fractions of the radius for the central box
    /// at the dividing lines and the timestepper.
    /// Timestepper defaults to Steady dummy timestepper.
    FullCircleMesh(GeomObject* wall_pt,
                   const Vector<double>& theta_positions,
                   const Vector<double>& radius_box,
                   TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Destructor: empty
    virtual ~FullCircleMesh()
    {
      delete Domain_pt;
    }

    /// Access function to GeomObject representing wall
    GeomObject*& area_pt()
    {
      return Area_pt;
    }

    /// Access function to domain
    FullCircleDomain* domain_pt()
    {
      return Domain_pt;
    }

    /// Access function to underlying domain
    FullCircleDomain* domain_pt() const
    {
      return Domain_pt;
    }

  protected:
    /// Pointer to domain
    FullCircleDomain* Domain_pt;

    /// Pointer to the geometric object that represents the entire domain
    GeomObject* Area_pt;
  };


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////


  //=============================================================
  /// Adaptative version of the FullCircleMesh base mesh.
  /// The domain is specified by the GeomObject that identifies
  /// the entire area
  ///
  /// The mesh boundaries are numbered as follows:
  /// - Boundary 1: The outer wall, represetned by \f$\xi_1 = 1\f$.
  /// .
  ///
  //=============================================================
  template<class ELEMENT>
  class RefineableFullCircleMesh : public FullCircleMesh<ELEMENT>,
                                   public RefineableQuadMesh<ELEMENT>

  {
  public:
    /// Constructor for adaptive deformable quarter tube mesh class.
    /// Pass pointer to geometric object that
    /// specifies the volume, start and end coordinates for the centreline
    /// on the geometric object. Values of theta at which dividing lines
    /// are to be placed, fractions of the radius for the central box
    /// at the dividing lines, the number of layers
    /// and the timestepper.
    /// Timestepper defaults to Steady dummy timestepper.
    RefineableFullCircleMesh(
      GeomObject* wall_pt,
      const Vector<double>& theta_positions,
      const Vector<double>& radius_box,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : FullCircleMesh<ELEMENT>(
          wall_pt, theta_positions, radius_box, time_stepper_pt)
    {
      // Loop over all elements and set macro element pointer
      for (unsigned ielem = 0; ielem < FullCircleMesh<ELEMENT>::nelement();
           ielem++)
      {
        dynamic_cast<RefineableQElement<2>*>(
          FullCircleMesh<ELEMENT>::element_pt(ielem))
          ->set_macro_elem_pt(this->Domain_pt->macro_element_pt(ielem));
      }


      // Setup Quadtree forest: Turn elements into individual octrees
      // and plant in forest
      Vector<TreeRoot*> trees_pt;
      for (unsigned iel = 0; iel < FullCircleMesh<ELEMENT>::nelement(); iel++)
      {
        FiniteElement* el_pt = FullCircleMesh<ELEMENT>::finite_element_pt(iel);
        ELEMENT* ref_el_pt = dynamic_cast<ELEMENT*>(el_pt);
        QuadTreeRoot* quadtree_root_pt = new QuadTreeRoot(ref_el_pt);
        trees_pt.push_back(quadtree_root_pt);
      }

      this->Forest_pt = new QuadTreeForest(trees_pt);

#ifdef PARANOID
      // Run self test
      unsigned success_flag =
        dynamic_cast<QuadTreeForest*>(this->Forest_pt)->self_test();
      if (success_flag == 0)
      {
        oomph_info << "Successfully built quadtree forest " << std::endl;
      }
      else
      {
        throw OomphLibError("Trouble in building quadtree forest ",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    /// Destructor: empty
    virtual ~RefineableFullCircleMesh() {}
  };

} // namespace oomph
#endif
