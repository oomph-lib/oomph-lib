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
// Header file for Tri/Tet linear elasticity elements
#ifndef OOMPH_TPML_TIME_HARMONIC_LINEAR_ELASTICITY_ELEMENTS_HEADER
#define OOMPH_TPML_TIME_HARMONIC_LINEAR_ELASTICITY_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// OOMPH-LIB headers
#include "generic/nodes.h"
#include "generic/oomph_utilities.h"
#include "generic/Telements.h"
#include "generic/error_estimator.h"

#include "./pml_time_harmonic_linear_elasticity_elements.h"

namespace oomph
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // TPMLTimeHarmonicLinearElasticityElement
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// TPMLTimeHarmonicLinearElasticityElement<DIM,NNODE_1D> elements are
  /// isoparametric triangular
  /// DIM-dimensional PMLTimeHarmonicLinearElasticity elements with
  /// NNODE_1D nodal points along each
  /// element edge. Inherits from TElement and
  /// PMLTimeHarmonicLinearElasticityEquations
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class TPMLTimeHarmonicLinearElasticityElement
    : public TElement<DIM, NNODE_1D>,
      public PMLTimeHarmonicLinearElasticityEquations<DIM>,
      public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// Constructor: Call constructors for TElement and
    /// PMLTimeHarmonicLinearElasticity equations
    TPMLTimeHarmonicLinearElasticityElement()
      : TElement<DIM, NNODE_1D>(),
        PMLTimeHarmonicLinearElasticityEquations<DIM>()
    {
    }


    /// Broken copy constructor
    TPMLTimeHarmonicLinearElasticityElement(
      const TPMLTimeHarmonicLinearElasticityElement<DIM, NNODE_1D>& dummy) =
      delete;

    /// Broken assignment operator
    void operator=(
      const TPMLTimeHarmonicLinearElasticityElement<DIM, NNODE_1D>&) = delete;

    /// Output function:
    void output(std::ostream& outfile)
    {
      PMLTimeHarmonicLinearElasticityEquations<DIM>::output(outfile);
    }

    ///  Output function:
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      PMLTimeHarmonicLinearElasticityEquations<DIM>::output(outfile, nplot);
    }


    /// C-style output function:
    void output(FILE* file_pt)
    {
      PMLTimeHarmonicLinearElasticityEquations<DIM>::output(file_pt);
    }

    ///  C-style output function:
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      PMLTimeHarmonicLinearElasticityEquations<DIM>::output(file_pt, n_plot);
    }


    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return TElement<DIM, NNODE_1D>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return TElement<DIM, NNODE_1D>::vertex_node_pt(j);
    }

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return NNODE_1D - 1;
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
  };

  //=======================================================================
  /// Face geometry for the TPMLTimeHarmonicLinearElasticityElement elements:
  /// The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<TPMLTimeHarmonicLinearElasticityElement<DIM, NNODE_1D>>
    : public virtual TElement<DIM - 1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : TElement<DIM - 1, NNODE_1D>() {}
  };

  //=======================================================================
  /// Face geometry for the 1D TPMLTimeHarmonicLinearElasticityElement
  /// elements:  Point elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<TPMLTimeHarmonicLinearElasticityElement<1, NNODE_1D>>
    : public virtual PointElement
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : PointElement() {}
  };


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Policy class defining the elements to be used in the actual
  /// PML layers. Same spatial dimension and nnode_1d but quads
  /// rather than triangles.
  //=======================================================================
  template<unsigned NNODE_1D>
  class PMLLayerElement<TPMLTimeHarmonicLinearElasticityElement<2, NNODE_1D>>
    : public virtual QPMLTimeHarmonicLinearElasticityElement<2, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate QElement
    PMLLayerElement() : QPMLTimeHarmonicLinearElasticityElement<2, NNODE_1D>()
    {
    }
  };

  //=======================================================================
  /// Face geometry for the TPMLTimeHarmonicLinearElasticityElement elements:
  /// The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<
    PMLLayerElement<TPMLTimeHarmonicLinearElasticityElement<DIM, NNODE_1D>>>
    : public virtual QElement<DIM - 1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QElement<DIM - 1, NNODE_1D>() {}
  };


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Policy class defining the elements to be used in the actual
  /// PML layers. Same spatial dimension and nnode_1d but quads
  /// rather than triangles.
  //=======================================================================
  template<unsigned NNODE_1D>
  class PMLLayerElement<ProjectablePMLTimeHarmonicLinearElasticityElement<
    TPMLTimeHarmonicLinearElasticityElement<2, NNODE_1D>>>
    : public virtual QPMLTimeHarmonicLinearElasticityElement<2, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate QElement
    PMLLayerElement() : QPMLTimeHarmonicLinearElasticityElement<2, NNODE_1D>()
    {
    }
  };


} // namespace oomph

#endif
