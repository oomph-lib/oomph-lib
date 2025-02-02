// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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

#ifndef OOMPH_ANNULAR_MESH_HEADER
#define OOMPH_ANNULAR_MESH_HEADER

// Headers
#include "rectangular_quadmesh.template.h"
#include "rectangular_quadmesh.template.cc"


// Include the headers file for domain
#include "annular_domain.h"

namespace oomph
{
  //===================================================================
  /// 2D annular mesh with a unit circle in the middle and a layer
  /// of thickness h surrounding it.
  //==================================================================
  template<class ELEMENT>
  class TwoDAnnularMesh : public virtual RectangularQuadMesh<ELEMENT>
  {
  public:
    /// Constructor
    TwoDAnnularMesh(const bool& periodic,
                    const double& azimuthal_fraction,
                    const unsigned& ntheta,
                    const unsigned& nr,
                    const double& a,
                    const double& h,
                    TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : RectangularQuadMesh<ELEMENT>(
          ntheta, nr, 1.0, 1.0, periodic, time_stepper_pt)
    {
      // Mesh can only be built with 2D Qelements.
      MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

      // Wrap mesh into annular shape
      double phi = 0.0;
      wrap_into_annular_shape(a, h, azimuthal_fraction, phi);
    }

    /// Constructor; rotate mesh by angle phi.
    TwoDAnnularMesh(const bool& periodic,
                    const double& azimuthal_fraction,
                    const unsigned& ntheta,
                    const unsigned& nr,
                    const double& a,
                    const double& h,
                    const double& phi,
                    TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : RectangularQuadMesh<ELEMENT>(
          ntheta, nr, 1.0, 1.0, periodic, time_stepper_pt)
    {
      // Mesh can only be built with 2D Qelements.
      MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

      // Wrap mesh into annular shape
      wrap_into_annular_shape(a, h, azimuthal_fraction, phi);
    }


  private:
    /// Wrap mesh into annular shape
    void wrap_into_annular_shape(const double& a,
                                 const double& h,
                                 const double& azimuthal_fraction,
                                 const double& phi);
  };


  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //===================================================================
  /// Refineable 2D annular mesh with a unit circle in the middle and a layer
  /// of thickness h surrounding it.
  //==================================================================
  template<class ELEMENT>
  class RefineableTwoDAnnularMesh : public virtual TwoDAnnularMesh<ELEMENT>,
                                    public virtual RefineableQuadMesh<ELEMENT>
  {
  public:
    /// Constructor
    RefineableTwoDAnnularMesh(
      const bool& periodic,
      const double& azimuthal_fraction,
      const unsigned& ntheta,
      const unsigned& nr,
      const double& a,
      const double& h,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : RectangularQuadMesh<ELEMENT>(
          ntheta, nr, 1.0, 1.0, periodic, time_stepper_pt),
        TwoDAnnularMesh<ELEMENT>(
          periodic, azimuthal_fraction, ntheta, nr, a, h, time_stepper_pt)
    {
      // Mesh can only be built with 2D Qelements.
      MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

      // Build macro element-based domain
      double phi = 0.0;
      Domain_pt = new AnnularDomain(azimuthal_fraction, ntheta, nr, a, h, phi);

      // Loop over all elements and set macro element pointer
      unsigned nel = this->nelement();
      for (unsigned ielem = 0; ielem < nel; ielem++)
      {
        dynamic_cast<RefineableQElement<2>*>(this->element_pt(ielem))
          ->set_macro_elem_pt(this->Domain_pt->macro_element_pt(ielem));
      }

      // Update nodal positions based on macro-element representation
      this->node_update();

      // Nodal positions etc. were created in constructor for
      // RectangularMesh<...>. Only need to setup quadtree forest
      this->setup_quadtree_forest();

      // Setup the periodic neighbour information of the TreeRoots
      // Cast to specific elements
      if (periodic)
      {
        Vector<TreeRoot*> left_root_pt(nr);
        Vector<TreeRoot*> right_root_pt(nr);
        for (unsigned i = 0; i < nr; i++)
        {
          left_root_pt[i] = dynamic_cast<ELEMENT*>(this->element_pt(i * ntheta))
                              ->tree_pt()
                              ->root_pt();

          right_root_pt[i] =
            dynamic_cast<ELEMENT*>(this->element_pt((i + 1) * ntheta - 1))
              ->tree_pt()
              ->root_pt();
        }

        // Set the neighbour and periodicity
        using namespace QuadTreeNames;
        for (unsigned i = 0; i < nr; i++)
        {
          left_root_pt[i]->neighbour_pt(W) = right_root_pt[i];
          left_root_pt[i]->set_neighbour_periodic(W);

          right_root_pt[i]->neighbour_pt(E) = left_root_pt[i];
          right_root_pt[i]->set_neighbour_periodic(E);
        }
      }
    }


    /// Constructor; rotate mesh by angle phi
    RefineableTwoDAnnularMesh(
      const bool& periodic,
      const double& azimuthal_fraction,
      const unsigned& ntheta,
      const unsigned& nr,
      const double& a,
      const double& h,
      const double& phi,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : RectangularQuadMesh<ELEMENT>(
          ntheta, nr, 1.0, 1.0, periodic, time_stepper_pt),
        TwoDAnnularMesh<ELEMENT>(
          periodic, azimuthal_fraction, ntheta, nr, a, h, phi, time_stepper_pt)
    {
      // Mesh can only be built with 2D Qelements.
      MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

      // Build macro element-based domain
      Domain_pt = new AnnularDomain(azimuthal_fraction, ntheta, nr, a, h, phi);

      // Loop over all elements and set macro element pointer
      unsigned nel = this->nelement();
      for (unsigned ielem = 0; ielem < nel; ielem++)
      {
        dynamic_cast<RefineableQElement<2>*>(this->element_pt(ielem))
          ->set_macro_elem_pt(this->Domain_pt->macro_element_pt(ielem));
      }

      // Update nodal positions based on macro-element representation
      this->node_update();

      // Nodal positions etc. were created in constructor for
      // RectangularMesh<...>. Only need to setup quadtree forest
      this->setup_quadtree_forest();

      // Setup the periodic neighbour information of the TreeRoots
      // Cast to specific elements
      if (periodic)
      {
        Vector<TreeRoot*> left_root_pt(nr);
        Vector<TreeRoot*> right_root_pt(nr);
        for (unsigned i = 0; i < nr; i++)
        {
          left_root_pt[i] = dynamic_cast<ELEMENT*>(this->element_pt(i * ntheta))
                              ->tree_pt()
                              ->root_pt();

          right_root_pt[i] =
            dynamic_cast<ELEMENT*>(this->element_pt((i + 1) * ntheta - 1))
              ->tree_pt()
              ->root_pt();
        }

        // Set the neighbour and periodicity
        using namespace QuadTreeNames;
        for (unsigned i = 0; i < nr; i++)
        {
          left_root_pt[i]->neighbour_pt(W) = right_root_pt[i];
          left_root_pt[i]->set_neighbour_periodic(W);

          right_root_pt[i]->neighbour_pt(E) = left_root_pt[i];
          right_root_pt[i]->set_neighbour_periodic(E);
        }
      }
    }

  private:
    /// Pointer to domain
    AnnularDomain* Domain_pt;
  };

} // namespace oomph

#endif
