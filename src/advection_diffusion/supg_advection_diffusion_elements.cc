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
#include "supg_advection_diffusion_elements.h"


namespace oomph
{
  //======================================================================
  /// \short QSUPGAdvectionDiffusionElement<DIM,NNODE_1D> elements are
  /// SUPG-stabilised Advection Diffusion elements with
  /// NNODE_1D nodal points in each coordinate direction. Inherits
  /// from QAdvectionDiffusionElement and overwrites their
  /// test functions
  ///
  //======================================================================

  //======================================================================
  /// \short Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// SUPG stabilisation: Petrov-Galerkin, i.e. test functions \f$ \ne \f$
  /// shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QSUPGAdvectionDiffusionElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_adv_diff(const Vector<double>& s,
                                       Shape& psi,
                                       DShape& dpsidx,
                                       Shape& test,
                                       DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = QElement<DIM, NNODE_1D>::dshape_eulerian(s, psi, dpsidx);

    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Calculate Eulerian coordinates
    Vector<double> interpolated_x(DIM, 0.0);

    // Loop over nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over directions
      for (unsigned j = 0; j < DIM; j++)
      {
        interpolated_x[j] += this->nodal_position(l, j) * psi[l];
      }
    }

    // Get wind
    Vector<double> wind(DIM);
    // Dummy ipt argument
    unsigned ipt = 0;
    this->get_wind_adv_diff(ipt, s, interpolated_x, wind);

    // Loop over the test functions and derivatives and set them equal to the
    // shape functions + add stabilisation
    for (unsigned j = 0; j < n_node; j++)
    {
      test[j] = psi[j];

      for (unsigned i = 0; i < DIM; i++)
      {
        dtestdx(j, i) = dpsidx(j, i);
        test[j] += Tau_SUPG * wind[i] * dpsidx(j, i);
      }
    }

    // Return the jacobian
    return J;
  }


  //======================================================================
  /// \short Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// SUPG stabilisation: Petrov-Galerkin, i.e. test functions \f$ \ne \f$
  /// shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QSUPGAdvectionDiffusionElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_adv_diff(const unsigned& ipt,
                                               Shape& psi,
                                               DShape& dpsidx,
                                               Shape& test,
                                               DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Calculate Eulerian coordinates
    Vector<double> interpolated_x(DIM, 0.0);

    // Loop over nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over directions
      for (unsigned j = 0; j < DIM; j++)
      {
        interpolated_x[j] += this->nodal_position(l, j) * psi(l);
      }
    }

    // Find the dimension of the element
    unsigned Dim = this->dim();
    // Storage for the local coordinates of the integration point
    Vector<double> s(Dim);
    // Set the local coordinate
    for (unsigned i = 0; i < Dim; i++)
    {
      s[i] = this->integral_pt()->knot(ipt, i);
    }

    // Get wind
    Vector<double> wind(DIM);
    this->get_wind_adv_diff(ipt, s, interpolated_x, wind);

    // Loop over the test functions and derivatives and set them equal to the
    // shape functions + add stabilisation
    for (unsigned j = 0; j < n_node; j++)
    {
      test(j) = psi(j);
      for (unsigned i = 0; i < DIM; i++)
      {
        dtestdx(j, i) = dpsidx(j, i);
        test(j) += Tau_SUPG * wind[i] * dpsidx(j, i);
      }
    }


    // Return the jacobian
    return J;
  }


  // Force template instantiation.
  template class QSUPGAdvectionDiffusionElement<2, 2>;
  template class QSUPGAdvectionDiffusionElement<2, 3>;
  template class QSUPGAdvectionDiffusionElement<2, 4>;


  template class QSUPGAdvectionDiffusionElement<3, 2>;
  template class QSUPGAdvectionDiffusionElement<3, 3>;
  template class QSUPGAdvectionDiffusionElement<3, 4>;


  // Force template instantiation.
  template class RefineableQSUPGAdvectionDiffusionElement<2, 2>;
  template class RefineableQSUPGAdvectionDiffusionElement<2, 3>;
  template class RefineableQSUPGAdvectionDiffusionElement<2, 4>;


  template class RefineableQSUPGAdvectionDiffusionElement<3, 2>;
  template class RefineableQSUPGAdvectionDiffusionElement<3, 3>;
  template class RefineableQSUPGAdvectionDiffusionElement<3, 4>;


} // namespace oomph
