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
#ifndef OOMPH_QUARTER_PIPE_MESH_TEMPLATE_HEADER
#define OOMPH_QUARTER_PIPE_MESH_TEMPLATE_HEADER

// Generic oomph-lib includes
#include "../generic/mesh.h"
#include "../generic/brick_mesh.h"
#include "../generic/refineable_brick_mesh.h"
#include "simple_cubic_mesh.template.h"
#include "simple_cubic_mesh.template.cc"
#include "quarter_pipe_domain.h"

#include "../generic/macro_element.h"
#include "../generic/domain.h"


namespace oomph
{
  //================================================================
  /// Non refineable quarter pipe mesh class
  /// Deform a simple cubic mesh into a quarter pipe
  /// r: radial direction
  /// theta: azimuthal direction
  /// z: axis direction
  //================================================================
  template<class ELEMENT>
  class QuarterPipeMesh : public virtual SimpleCubicMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass number of elements in various directions,
    /// the inner and outer radius and the length of the tube
    QuarterPipeMesh(const unsigned& ntheta,
                    const unsigned& nr,
                    const unsigned& nz,
                    const double& rmin,
                    const double& rmax,
                    const double& length,
                    TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);


    /// Empty Destructor
    virtual ~QuarterPipeMesh()
    {
      delete Domain_pt;
    }

    /// Access function to domain
    QuarterPipeDomain* domain_pt()
    {
      return Domain_pt;
    }

    /// Access function to underlying domain
    QuarterPipeDomain* domain_pt() const
    {
      return Domain_pt;
    }

  protected:
    /// Number of elements azimuthal direction
    unsigned Ntheta;

    /// Number of elements radial direction
    unsigned Nr;

    /// Number of elements axial direction
    unsigned Nz;

    /// Inner radius
    double Rmin;

    /// Outer radius
    double Rmax;

    /// Length
    double Length;

    /// Pointer to domain
    QuarterPipeDomain* Domain_pt;

  }; // endofclass


  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


  //================================================================
  /// Refineable quarter pipe mesh class
  //================================================================
  template<class ELEMENT>
  class RefineableQuarterPipeMesh : public virtual QuarterPipeMesh<ELEMENT>,
                                    public RefineableBrickMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass number of elements in various directions,
    /// the inner and outer radius and the length of the tube
    RefineableQuarterPipeMesh(
      const unsigned& ntheta,
      const unsigned& nr,
      const unsigned& nz,
      const double& rmin,
      const double& rmax,
      const double& length,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : SimpleCubicMesh<ELEMENT>(
          ntheta, nr, nz, 1.0, 1.0, 1.0, time_stepper_pt),
        QuarterPipeMesh<ELEMENT>(
          ntheta, nr, nz, rmin, rmax, length, time_stepper_pt)
    {
      // Setup Octree forest: Turn elements into individual octrees
      // and plant in forest
      Vector<TreeRoot*> trees_pt;
      for (unsigned iel = 0; iel < (nr * ntheta * nz); iel++)
      {
        FiniteElement* el_pt = QuarterPipeMesh<ELEMENT>::finite_element_pt(iel);
        ELEMENT* ref_el_pt = dynamic_cast<ELEMENT*>(el_pt);
        OcTreeRoot* octree_root_pt = new OcTreeRoot(ref_el_pt);
        trees_pt.push_back(octree_root_pt);
      }

      this->Forest_pt = new OcTreeForest(trees_pt);
    }


    /// Destructor -- delete forest
    virtual ~RefineableQuarterPipeMesh()
    {
      delete this->Forest_pt;
    }

  }; // endofclass


  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


  //================================================================
  /// Non refineable elastic quarter pipe mesh class
  /// setup lagrangian coordinates for solid mechanics problems
  //================================================================
  template<class ELEMENT>
  class ElasticQuarterPipeMesh : public virtual QuarterPipeMesh<ELEMENT>,
                                 public virtual SolidMesh
  {
  public:
    /// Constructor:  Pass number of elements in various directions,
    /// the inner and outer radius and the length of the tube.
    /// Builds mesh and copies Eulerian coords to Lagrangian
    /// ones so that the initial configuration is the stress-free one.
    ElasticQuarterPipeMesh(
      const unsigned& ntheta,
      const unsigned& nr,
      const unsigned& nz,
      const double& rmin,
      const double& rmax,
      const double& length,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : SimpleCubicMesh<ELEMENT>(
          ntheta, nr, nz, 1.0, 1.0, 1.0, time_stepper_pt),
        QuarterPipeMesh<ELEMENT>(
          ntheta, nr, nz, rmin, rmax, length, time_stepper_pt)
    {
      /// Make the current configuration the undeformed one by
      /// setting the nodal Lagrangian coordinates to their current
      /// Eulerian ones
      set_lagrangian_nodal_coordinates();
    }
  };


  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


  //================================================================
  /// Refineable elastic quarter pipe mesh class
  //================================================================
  template<class ELEMENT>
  class ElasticRefineableQuarterPipeMesh
    : public virtual ElasticQuarterPipeMesh<ELEMENT>,
      public RefineableBrickMesh<ELEMENT>
  {
  public:
    /// Constructor:  Pass number of elements in various directions,
    /// the inner and outer radius and the length of the tube.
    /// Builds mesh and copies Eulerian coords to Lagrangian
    /// ones so that the initial configuration is the stress-free one.
    ElasticRefineableQuarterPipeMesh(
      const unsigned& ntheta,
      const unsigned& nr,
      const unsigned& nz,
      const double& rmin,
      const double& rmax,
      const double& length,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : SimpleCubicMesh<ELEMENT>(
          ntheta, nr, nz, 1.0, 1.0, 1.0, time_stepper_pt),
        QuarterPipeMesh<ELEMENT>(
          ntheta, nr, nz, rmin, rmax, length, time_stepper_pt),
        ElasticQuarterPipeMesh<ELEMENT>(
          ntheta, nr, nz, rmin, rmax, length, time_stepper_pt)
    {
      // Setup Octree forest: Turn elements into individual octrees
      // and plant in forest
      Vector<TreeRoot*> trees_pt;
      for (unsigned iel = 0; iel < (nr * ntheta * nz); iel++)
      {
        FiniteElement* el_pt = QuarterPipeMesh<ELEMENT>::finite_element_pt(iel);
        ELEMENT* ref_el_pt = dynamic_cast<ELEMENT*>(el_pt);
        OcTreeRoot* octree_root_pt = new OcTreeRoot(ref_el_pt);
        trees_pt.push_back(octree_root_pt);
      }
      this->Forest_pt = new OcTreeForest(trees_pt);

      // Loop over all elements and set the undeformed macro element pointer
      unsigned n_element = this->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        // Get pointer to full element type
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(this->element_pt(e));

        // Set pointer to macro element so the curvlinear boundaries
        // of the undeformed mesh/domain get picked up during adaptive
        // mesh refinement
        el_pt->set_undeformed_macro_elem_pt(
          this->Domain_pt->macro_element_pt(e));

        // Use MacroElement representation for
        // Lagrangian coordinates of newly created
        // nodes
        el_pt
          ->enable_use_of_undeformed_macro_element_for_new_lagrangian_coords();
      }
    }
  };


} // namespace oomph


#endif
