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
// Header file for refineable linear elasticity elements

// Include guard to prevent multiple inclusions of this header
#ifndef OOMPH_REFINEABLE_TIME_HARMONIC_LINEAR_ELASTICITY_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_TIME_HARMONIC_LINEAR_ELASTICITY_ELEMENTS_HEADER

// oomph-lib headers
#include "time_harmonic_linear_elasticity_elements.h"
#include "../generic/refineable_quad_element.h"
#include "../generic/refineable_brick_element.h"
#include "../generic/error_estimator.h"

namespace oomph
{
  //========================================================================
  /// Class for Refineable TimeHarmonicLinearElasticity equations
  //========================================================================
  template<unsigned DIM>
  class RefineableTimeHarmonicLinearElasticityEquations
    : public virtual TimeHarmonicLinearElasticityEquations<DIM>,
      public virtual RefineableElement,
      public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// Constructor
    RefineableTimeHarmonicLinearElasticityEquations()
      : TimeHarmonicLinearElasticityEquations<DIM>(),
        RefineableElement(),
        ElementWithZ2ErrorEstimator()
    {
    }


    /// Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Create enough initialised storage
      values.resize(2 * DIM, 0.0);

      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Shape functions
      Shape psi(n_node);
      this->shape(s, psi);

      // Calculate displacements
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get the index at which the i-th velocity is stored
        std::complex<unsigned> u_nodal_index =
          this->u_index_time_harmonic_linear_elasticity(i);
        for (unsigned l = 0; l < n_node; l++)
        {
          values[i] += this->nodal_value(t, l, u_nodal_index.real()) * psi(l);
          values[i + DIM] +=
            this->nodal_value(t, l, u_nodal_index.imag()) * psi(l);
        }
      }
    }

    /// Get the current interpolated values (displacements).
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines) ,the values Vector sets its own size in here.
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      this->get_interpolated_values(0, s, values);
    }

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      // DIM Diagonal strain rates and DIM*(DIM-1)/2 off diagonal terms
      return 2 * (DIM + DIM * (DIM - 1) / 2);
    }

    /// Get 'flux' for Z2 error recovery:   Upper triangular entries
    /// in strain tensor.
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
#ifdef PARANOID
      unsigned num_entries = 2 * (DIM + ((DIM * DIM) - DIM) / 2);
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
      DenseMatrix<std::complex<double>> strain(DIM);
      this->get_strain(s, strain);

      // Pack into flux Vector
      unsigned icount = 0;

      // Start with diagonal terms
      for (unsigned i = 0; i < DIM; i++)
      {
        flux[icount] = strain(i, i).real();
        icount++;
        flux[icount] = strain(i, i).imag();
        icount++;
      }

      // Off diagonals row by row
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = i + 1; j < DIM; j++)
        {
          flux[icount] = strain(i, j).real();
          icount++;
          flux[icount] = strain(i, j).imag();
          icount++;
        }
      }
    }

    /// Number of continuously interpolated values: 2*DIM
    unsigned ncont_interpolated_values() const
    {
      return 2 * DIM;
    }

    /// Further build function, pass the pointers down to the sons
    void further_build()
    {
      RefineableTimeHarmonicLinearElasticityEquations<DIM>*
        cast_father_element_pt =
          dynamic_cast<RefineableTimeHarmonicLinearElasticityEquations<DIM>*>(
            this->father_element_pt());

      // Set pointer to body force function
      this->Body_force_fct_pt = cast_father_element_pt->body_force_fct_pt();

      // Set pointer to the contitutive law
      this->Elasticity_tensor_pt =
        cast_father_element_pt->elasticity_tensor_pt();

      // Set the frequency
      this->Omega_sq_pt = cast_father_element_pt->omega_sq_pt();
    }


  private:
    /// Overloaded helper function to take hanging nodes into account
    void fill_in_generic_contribution_to_residuals_time_harmonic_linear_elasticity(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag);
  };


  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////

  //========================================================================
  /// Class for refineable QTimeHarmonicLinearElasticityElement elements
  //========================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class RefineableQTimeHarmonicLinearElasticityElement
    : public virtual QTimeHarmonicLinearElasticityElement<DIM, NNODE_1D>,
      public virtual RefineableTimeHarmonicLinearElasticityEquations<DIM>,
      public virtual RefineableQElement<DIM>
  {
  public:
    /// Constructor:
    RefineableQTimeHarmonicLinearElasticityElement()
      : QTimeHarmonicLinearElasticityElement<DIM, NNODE_1D>(),
        RefineableElement(),
        RefineableTimeHarmonicLinearElasticityEquations<DIM>(),
        RefineableQElement<DIM>()
    {
    }

    /// Empty rebuild from sons, no need to reconstruct anything here
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QTimeHarmonicLinearElasticityElement<DIM,
                                                  NNODE_1D>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QTimeHarmonicLinearElasticityElement<DIM,
                                                  NNODE_1D>::vertex_node_pt(j);
    }

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return NNODE_1D - 1;
    }

    ///  No additional hanging node procedures are required
    void further_setup_hanging_nodes() {}
  };


  //==============================================================
  /// FaceGeometry of the 2D
  /// RefineableQTimeHarmonicLinearElasticityElement elements
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<
    RefineableQTimeHarmonicLinearElasticityElement<2, NNODE_1D>>
    : public virtual QElement<1, NNODE_1D>
  {
  public:
    // Make sure that we call the constructor of the QElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : QElement<1, NNODE_1D>() {}
  };

  //==============================================================
  /// FaceGeometry of the FaceGeometry of the 2D
  /// RefineableQTimeHarmonicLinearElasticityElement
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<
    FaceGeometry<RefineableQTimeHarmonicLinearElasticityElement<2, NNODE_1D>>>
    : public virtual PointElement
  {
  public:
    // Make sure that we call the constructor of the QElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : PointElement() {}
  };


  //==============================================================
  /// FaceGeometry of the 3D RefineableQTimeHarmonicLinearElasticityElement
  /// elements
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<
    RefineableQTimeHarmonicLinearElasticityElement<3, NNODE_1D>>
    : public virtual QElement<2, NNODE_1D>
  {
  public:
    // Make sure that we call the constructor of the QElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : QElement<2, NNODE_1D>() {}
  };

  //==============================================================
  /// FaceGeometry of the FaceGeometry of the 3D
  /// RefineableQTimeHarmonicLinearElasticityElement
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<
    FaceGeometry<RefineableQTimeHarmonicLinearElasticityElement<3, NNODE_1D>>>
    : public virtual QElement<1, NNODE_1D>
  {
  public:
    // Make sure that we call the constructor of the QElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : QElement<1, NNODE_1D>() {}
  };


} // namespace oomph

#endif
