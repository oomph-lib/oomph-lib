// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
// Header file for TFourierDecomposedHelmholtz elements
#ifndef OOMPH_TFOURIER_DECOMPOSED_HELMHOLTZ_ELEMENTS_HEADER
#define OOMPH_TFOURIER_DECOMPOSED_HELMHOLTZ_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// OOMPH-LIB headers
#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"
#include "../generic/Telements.h"
#include "../generic/error_estimator.h"
#include "fourier_decomposed_helmholtz_elements.h"

namespace oomph
{
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  // TFourierDecomposedHelmholtzElement
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //======================================================================
  /// TFourierDecomposedHelmholtzElement<NNODE_1D> elements are
  /// isoparametric triangular
  /// FourierDecomposedHelmholtz elements with  NNODE_1D nodal points along each
  /// element edge. Inherits from TElement and
  /// FourierDecomposedHelmholtzEquations
  //======================================================================
  template<unsigned NNODE_1D>
  class TFourierDecomposedHelmholtzElement
    : public TElement<2, NNODE_1D>,
      public FourierDecomposedHelmholtzEquations,
      public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// Constructor: Call constructors for TElement and
    /// FourierDecomposedHelmholtz equations
    TFourierDecomposedHelmholtzElement()
      : TElement<2, NNODE_1D>(), FourierDecomposedHelmholtzEquations()
    {
    }


    /// Broken copy constructor
    TFourierDecomposedHelmholtzElement(
      const TFourierDecomposedHelmholtzElement<NNODE_1D>& dummy) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const TFourierDecomposedHelmholtzElement<NNODE_1D>&) =
     * delete;*/

    ///  Access function for Nvalue: # of `values' (pinned or dofs)
    /// at node n (always returns the same value at every node, 2)
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }

    /// Output function:
    ///  r,z,u
    void output(std::ostream& outfile)
    {
      FourierDecomposedHelmholtzEquations::output(outfile);
    }

    ///  Output function:
    ///   r,z,u  n_plot^2 plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      FourierDecomposedHelmholtzEquations::output(outfile, n_plot);
    }


    /// C-style output function:
    ///  r,z,u   or    x,y,z,u
    void output(FILE* file_pt)
    {
      FourierDecomposedHelmholtzEquations::output(file_pt);
    }


    ///  C-style output function:
    ///   r,z,u  at n_plot^2 plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FourierDecomposedHelmholtzEquations::output(file_pt, n_plot);
    }


    /// Output function for an exact solution:
    ///  r,z,u_exact
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      FourierDecomposedHelmholtzEquations::output_fct(
        outfile, n_plot, exact_soln_pt);
    }


    /// Output function for a time-dependent exact solution.
    ///  x,y,u_exact (calls the steady version)
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      FourierDecomposedHelmholtzEquations::output_fct(
        outfile, n_plot, time, exact_soln_pt);
    }

  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_fourier_decomposed_helmholtz(
      const Vector<double>& s,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;


    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_fourier_decomposed_helmholtz(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;


    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return (NNODE_1D - 1);
    }

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      return 2 * 2;
    }

    /// Get 'flux' for Z2 error recovery:  Standard flux from
    /// UnsteadyHeat equations
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      Vector<std::complex<double>> complex_flux(2);
      this->get_flux(s, complex_flux);
      unsigned count = 0;
      for (unsigned i = 0; i < 2; i++)
      {
        flux[count++] = complex_flux[i].real();
        flux[count++] = complex_flux[i].imag();
      }
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

  private:
    /// Static unsigned that holds the (same) number of variables at every node
    static const unsigned Initial_Nvalue;
  };


  // Inline functions:


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned NNODE_1D>
  double TFourierDecomposedHelmholtzElement<NNODE_1D>::
    dshape_and_dtest_eulerian_fourier_decomposed_helmholtz(
      const Vector<double>& s,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const
  {
    unsigned n_node = this->nnode();

    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian(s, psi, dpsidx);

    // Loop over the test functions and derivatives and set them equal to the
    // shape functions
    for (unsigned i = 0; i < n_node; i++)
    {
      test[i] = psi[i];
      dtestdx(i, 0) = dpsidx(i, 0);
      dtestdx(i, 1) = dpsidx(i, 1);
    }

    // Return the jacobian
    return J;
  }


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned NNODE_1D>
  double TFourierDecomposedHelmholtzElement<NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_fourier_decomposed_helmholtz(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

    // Set the pointers of the test functions
    test = psi;
    dtestdx = dpsidx;

    // Return the jacobian
    return J;
  }


  //=======================================================================
  /// Face geometry for the TFourierDecomposedHelmholtzElement elements:
  /// The spatial dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<TFourierDecomposedHelmholtzElement<NNODE_1D>>
    : public virtual TElement<1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : TElement<1, NNODE_1D>() {}
  };


} // namespace oomph

#endif
