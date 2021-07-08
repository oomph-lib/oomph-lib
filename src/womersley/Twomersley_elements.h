// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
// LIC//
// LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
// Header file for TWomersley elements
#ifndef OOMPH_TWOMERSLEY_ELEMENTS_HEADER
#define OOMPH_TWOMERSLEY_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// OOMPH-LIB headers
#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"
#include "../generic/Telements.h"
#include "../generic/error_estimator.h"
#include "womersley_elements.h"

namespace oomph
{
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  //======================================================================
  /// TWomersleyElement elements are linear/triangular/tetrahedral
  /// Womersley elements with isoparametric interpolation for the function.
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class TWomersleyElement :
    public virtual TElement<DIM, NNODE_1D>,
    public virtual WomersleyEquations<DIM>
  {
  private:
    /// \short Static array of ints to hold number of variables at
    /// nodes: Initial_Nvalue[n]
    static const unsigned Initial_Nvalue;

  public:
    ///\short  Constructor: Call constructors for TElement and
    /// Womersley equations
    TWomersleyElement() : TElement<DIM, NNODE_1D>(), WomersleyEquations<DIM>()
    {
    }

    /// Broken copy constructor
    TWomersleyElement(const TWomersleyElement<DIM, NNODE_1D> &dummy)
    {
      BrokenCopy::broken_copy("TWomersleyElement");
    }

    /// Broken assignment operator
    void operator=(const TWomersleyElement<DIM, NNODE_1D> &)
    {
      BrokenCopy::broken_assign("TWomersleyElement");
    }

    /// \short  Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned &n) const
    {
      return Initial_Nvalue;
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream &outfile)
    {
      WomersleyEquations<DIM>::output(outfile);
    }

    ///  \short Output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot)
    {
      WomersleyEquations<DIM>::output(outfile, n_plot);
    }

    /// \short C-style output function:
    ///  x,y,u   or    x,y,z,u
    void output(FILE *file_pt)
    {
      WomersleyEquations<DIM>::output(file_pt);
    }

    ///  \short C-style output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(FILE *file_pt, const unsigned &n_plot)
    {
      WomersleyEquations<DIM>::output(file_pt, n_plot);
    }

    /// \short Output function for an exact solution:
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    void output_fct(std::ostream &outfile,
                    const unsigned &n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      WomersleyEquations<DIM>::output_fct(outfile, n_plot, exact_soln_pt);
    }

    /// \short Output function for a time-dependent exact solution.
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    /// (Calls the steady version)
    void output_fct(std::ostream &outfile,
                    const unsigned &n_plot,
                    const double &time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      WomersleyEquations<DIM>::output_fct(outfile, n_plot, time, exact_soln_pt);
    }

  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_womersley(const Vector<double> &s,
                                                      Shape &psi,
                                                      DShape &dpsidx,
                                                      Shape &test,
                                                      DShape &dtestdx) const;

    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return  Jacobian of mapping
    inline double dshape_and_dtest_eulerian_at_knot_womersley(
      const unsigned &ipt,
      Shape &psi,
      DShape &dpsidx,
      Shape &test,
      DShape &dtestdx) const;
  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double TWomersleyElement<DIM, NNODE_1D>::dshape_and_dtest_eulerian_womersley(
    const Vector<double> &s,
    Shape &psi,
    DShape &dpsidx,
    Shape &test,
    DShape &dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian(s, psi, dpsidx);

    // Loop over the test functions and derivatives and set them equal to the
    // shape functions
    for (unsigned i = 0; i < NNODE_1D; i++)
    {
      test[i] = psi[i];
      for (unsigned j = 0; j < DIM; j++)
      {
        dtestdx(i, j) = dpsidx(i, j);
      }
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
  double TWomersleyElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_womersley(const unsigned &ipt,
                                                Shape &psi,
                                                DShape &dpsidx,
                                                Shape &test,
                                                DShape &dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

    // Set the test functions equal to the shape functions
    //(sets internal pointers)
    test = psi;
    dtestdx = dpsidx;

    // Return the jacobian
    return J;
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Face geometry for the TWomersleyElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<TWomersleyElement<DIM, NNODE_1D>> :
    public virtual TElement<DIM - 1, NNODE_1D>
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : TElement<DIM - 1, NNODE_1D>() {}
  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Face geometry for the 1D TWomersleyElement elements: Point elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<TWomersleyElement<1, NNODE_1D>> :
    public virtual PointElement
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : PointElement() {}
  };

} // namespace oomph

#endif
