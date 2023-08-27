//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#ifndef TINTEGRAL_ELEMENTS_HEADER
#define TINTEGRAL_ELEMENTS_HEADER

#include "../generic/Telements.h"
#include "../generic/error_estimator.h"
#include "integral_equations.h"

namespace oomph
{
  template<unsigned DIM, unsigned NNODE_1D>
  class TIntegralElement : public virtual TElement<2, NNODE_1D>,
                           public virtual IntegralEquations<2>,
                           public virtual ElementWithZ2ErrorEstimator
  {
  public:
    // Constructor: Call constructors for TElement and
    /// Integral equations
    TIntegralElement() : TElement<2, NNODE_1D>(), IntegralEquations<2>() {}

    ///  Access function for Nvalue: # of `values' (pinned or dofs)
    /// at node n (always returns the same value at every node, 1)
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return (NNODE_1D - 1);
    }

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      return 2;
    }

    /// Get 'flux' for Z2 error recovery:  Standard flux.from Poisson equations
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      this->get_integral_flux(s, flux);
    }

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return TElement<2, NNODE_1D>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return TElement<2, NNODE_1D>::vertex_node_pt(j);
    }

    virtual void output(std::ostream& outfile)
    {
      output(outfile, NNODE_1D);
    }

    void output(std::ostream& outfile, const unsigned& nplot)
    {
      IntegralEquations<2>::output(outfile, nplot);
    }

  private:
    /// Static unsigned that holds the (same) number of variables at every node
    static const unsigned Initial_Nvalue;
  };

  //=======================================================================
  /// Face geometry for the TIntegralElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<TIntegralElement<DIM, NNODE_1D>>
    : public virtual TElement<DIM - 1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : TElement<DIM - 1, NNODE_1D>() {}
  };

  //======================================================================
  // Set the data for the number of Variables at each node, always 1
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned TIntegralElement<DIM, NNODE_1D>::Initial_Nvalue = 0;

} // namespace oomph
#endif
