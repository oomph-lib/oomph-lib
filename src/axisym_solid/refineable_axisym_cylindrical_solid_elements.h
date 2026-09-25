// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2026 Matthias Heil and Andrew Hazel
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

// Include guard to prevent multiple inclusions of this header
#ifndef OOMPH_REFINEABLE_AXISYM_CYLINDRICAL_ELASTICITY_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_AXISYM_CYLINDRICAL_ELASTICITY_ELEMENTS_HEADER

// oomph-lib headers
#include "axisym_cylindrical_solid_elements.h"
#include "../generic/refineable_quad_element.h"
#include "../generic/refineable_brick_element.h"
#include "../generic/error_estimator.h"

namespace oomph
{
  //========================================================================
  /// Class for Refineable axisymmetric PVD equations in cylindrical coords
  //========================================================================
  class RefineableAxisymmetricCylindricalPVDEquations
    : public virtual AxisymmetricCylindricalPVDEquations,
      public virtual RefineableSolidElement,
      public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// Constructor
    RefineableAxisymmetricCylindricalPVDEquations()
      : AxisymmetricCylindricalPVDEquations(),
        RefineableElement(),
        RefineableSolidElement(),
        ElementWithZ2ErrorEstimator()
    {
    }

    /// Call the residuals including hanging node cases
    void fill_in_contribution_to_residuals_axisym_pvd(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag);

    /// No values are interpolated in this element (pure solid)
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      values.clear();
    }

    /// No values are interpolated in this element (pure solid)
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      values.clear();
    }

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      return 4;
    }

    /// Get 'flux' for Z2 error recovery:   Upper triangular entries
    /// in strain tensor.
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
#ifdef PARANOID
      unsigned num_entries = 4;
      if (flux.size() != num_entries)
      {
        std::ostringstream error_message;
        error_message << "The flux vector has the wrong number of entries, "
                      << flux.size() << ", whereas it should be " << num_entries
                      << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get strain matrix
      DenseMatrix<double> strain(3);
      this->get_strain(s, strain);

      // Pack into flux Vector
      unsigned icount = 0;

      // Start with diagonal terms
      for (unsigned i = 0; i < 3; i++)
      {
        flux[icount] = strain(i, i);
        icount++;
      }

      // Only non-zero off-diagonal
      flux[icount] = strain(0, 1);
    }

    /// Number of continuously interpolated values: 0 (pure solid problem)
    unsigned ncont_interpolated_values() const
    {
      return 0;
    }

    /// Further build function, pass the pointers down to the sons
    void further_build()
    {
      RefineableAxisymmetricCylindricalPVDEquations* cast_father_element_pt =
        dynamic_cast<RefineableAxisymmetricCylindricalPVDEquations*>(
          this->father_element_pt());

      // Do whatever needs to be done in the base class
      RefineableSolidElement::further_build();

      // Set pointer to body force function
      this->Body_force_fct_pt = cast_father_element_pt->body_force_fct_pt();

      // Set pointer to the contitutive law
      this->Constitutive_law_pt = cast_father_element_pt->constitutive_law_pt();

      // Set the timescale ratio (non-dim. density)
      this->Lambda_sq_pt = cast_father_element_pt->lambda_sq_pt();

      // Set the mass damping parameter
      this->Eta_mass_pt = cast_father_element_pt->eta_mass_pt();
    }
  };

  //========================================================================
  /// Class for refineable QPVDElement elements
  //========================================================================
  template<unsigned NNODE_1D>
  class RefineableQAxisymmetricCylindricalPVDElement
    : public virtual QAxisymmetricCylindricalPVDElement<NNODE_1D>,
      public virtual RefineableAxisymmetricCylindricalPVDEquations,
      public virtual RefineableSolidQElement<2>
  {
  public:
    /// Constructor:
    RefineableQAxisymmetricCylindricalPVDElement()
      : QAxisymmetricCylindricalPVDElement<NNODE_1D>(),
        RefineableElement(),
        RefineableSolidElement(),
        RefineableAxisymmetricCylindricalPVDEquations(),
        RefineableSolidQElement<2>()
    {
    }

    /// Empty rebuild from sons, no need to reconstruct anything here
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QAxisymmetricCylindricalPVDElement<NNODE_1D>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QAxisymmetricCylindricalPVDElement<NNODE_1D>::vertex_node_pt(j);
    }

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return NNODE_1D - 1;
    }

    ///  No additional hanging node procedures are required
    /// for the solid elements.
    void further_setup_hanging_nodes() {}
  };

  //==============================================================
  /// FaceGeometry of the
  /// RefineableQAxisymmetricPVDElement
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<RefineableQAxisymmetricCylindricalPVDElement<NNODE_1D>>
    : public virtual SolidQElement<1, NNODE_1D>
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<1, NNODE_1D>() {}
  };

  //==============================================================
  /// FaceGeometry of the FaceGeometry of the
  /// RefineableQAxisymmetricPVDElement
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<
    FaceGeometry<RefineableQAxisymmetricCylindricalPVDElement<NNODE_1D>>>
    : public virtual PointElement
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : PointElement() {}
  };

} // namespace oomph

#endif
