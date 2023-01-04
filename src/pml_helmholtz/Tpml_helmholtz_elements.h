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
// Header file for THelmholtz elements
#ifndef OOMPH_TPML_HELMHOLTZ_ELEMENTS_HEADER
#define OOMPH_TPML_HELMHOLTZ_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// OOMPH-LIB headers
#include "generic/nodes.h"
#include "generic/oomph_utilities.h"
#include "generic/Telements.h"
#include "generic/error_estimator.h"
#include "pml_helmholtz_elements.h"

namespace oomph
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // TPMLHelmholtzElement
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// TPMLHelmholtzElement<DIM,NNODE_1D> elements are
  /// isoparametric triangular DIM-dimensional PMLHelmholtz elements
  /// with  NNODE_1D nodal points along each element edge. Inherits from
  /// TElement and PMLHelmholtzEquations
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class TPMLHelmholtzElement : public TElement<DIM, NNODE_1D>,
                               public PMLHelmholtzEquations<DIM>,
                               public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// Constructor: Call constructors for TElement and
    /// PMLHelmholtz equations
    TPMLHelmholtzElement()
      : TElement<DIM, NNODE_1D>(), PMLHelmholtzEquations<DIM>()
    {
    }


    /// Broken copy constructor
    TPMLHelmholtzElement(const TPMLHelmholtzElement<DIM, NNODE_1D>& dummy) =
      delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const TPMLHelmholtzElement<DIM,NNODE_1D>&) = delete;*/

    ///  Access function for Nvalue: # of `values' (pinned or dofs)
    /// at node n (always returns the same value at every node, 1)
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }

    /// Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream& outfile)
    {
      PMLHelmholtzEquations<DIM>::output(outfile);
    }

    ///  Output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      PMLHelmholtzEquations<DIM>::output(outfile, n_plot);
    }


    /// C-style output function:
    ///  x,y,u   or    x,y,z,u
    void output(FILE* file_pt)
    {
      PMLHelmholtzEquations<DIM>::output(file_pt);
    }


    ///  C-style output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      PMLHelmholtzEquations<DIM>::output(file_pt, n_plot);
    }


    /// Output function for an exact solution:
    ///  x,y,u_exact
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      PMLHelmholtzEquations<DIM>::output_fct(outfile, n_plot, exact_soln_pt);
    }


    /// Output function for a time-dependent exact solution.
    ///  x,y,u_exact (calls the steady version)
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      PMLHelmholtzEquations<DIM>::output_fct(
        outfile, n_plot, time, exact_soln_pt);
    }

  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_helmholtz(const Vector<double>& s,
                                                      Shape& psi,
                                                      DShape& dpsidx,
                                                      Shape& test,
                                                      DShape& dtestdx) const;


    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_helmholtz(
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
      return 2 * DIM;
    }

    /// Get 'flux' for Z2 error recovery:  Standard flux from
    /// UnsteadyHeat equations
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      Vector<std::complex<double>> complex_flux(DIM);
      this->get_flux(s, complex_flux);
      unsigned count = 0;
      for (unsigned i = 0; i < DIM; i++)
      {
        flux[count++] = complex_flux[i].real();
        flux[count++] = complex_flux[i].imag();
      }
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

  private:
    /// Static unsigned that holds the (same) number of variables at every node
    static const unsigned Initial_Nvalue;
  };


  //!! Cleanup - this was not here before!
  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned TPMLHelmholtzElement<DIM, NNODE_1D>::Initial_Nvalue = 2;

  // Inline functions:


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double TPMLHelmholtzElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_helmholtz(const Vector<double>& s,
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
  template<unsigned DIM, unsigned NNODE_1D>
  double TPMLHelmholtzElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_helmholtz(const unsigned& ipt,
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
  /// Face geometry for the TPMLHelmholtzElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<TPMLHelmholtzElement<DIM, NNODE_1D>>
    : public virtual TElement<DIM - 1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : TElement<DIM - 1, NNODE_1D>() {}
  };

  //=======================================================================
  /// Face geometry for the 1D TPMLHelmholtzElement elements:
  /// Point elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<TPMLHelmholtzElement<1, NNODE_1D>>
    : public virtual PointElement
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : PointElement() {}
  };


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Policy class defining the elements to be used in the actual
  /// PML layers. It's the corresponding quads.
  //=======================================================================
  template<unsigned NNODE_1D>
  class PMLLayerElement<TPMLHelmholtzElement<2, NNODE_1D>>
    : public virtual QPMLHelmholtzElement<2, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate QElement
    PMLLayerElement() : QPMLHelmholtzElement<2, NNODE_1D>() {}
  };

  //=======================================================================
  /// Face geometry for the TPMLHelmholtzElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<PMLLayerElement<TPMLHelmholtzElement<DIM, NNODE_1D>>>
    : public virtual QElement<DIM - 1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : QElement<DIM - 1, NNODE_1D>() {}
  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Policy class defining the elements to be used in the actual
  /// PML layers. It's the corresponding quads.
  //=======================================================================
  template<unsigned NNODE_1D>
  class PMLLayerElement<
    ProjectablePMLHelmholtzElement<TPMLHelmholtzElement<2, NNODE_1D>>>
    : public virtual QPMLHelmholtzElement<2, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate QElement
    PMLLayerElement() : QPMLHelmholtzElement<2, NNODE_1D>() {}
  };

} // namespace oomph


#endif
