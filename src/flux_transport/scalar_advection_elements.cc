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
// Non-inline member function of the flux transport elements class

#include "scalar_advection_elements.h"

namespace oomph
{
  //=========================================================
  /// Return the flux as a function of the unknowns
  //=========================================================
  template<unsigned DIM>
  void ScalarAdvectionEquations<DIM>::flux(const Vector<double>& u,
                                           DenseMatrix<double>& f)
  {
    // Get the wind
    Vector<double> W(DIM);
    Vector<double> s(DIM), x(DIM);
    // Dummy integration point
    unsigned ipt = 0;
    this->get_wind_scalar_adv(ipt, s, x, W);

    // Flux is the wind multiplied by the flux
    for (unsigned j = 0; j < DIM; j++)
    {
      f(0, j) = W[j] * u[0];
    }
  }

  //======================================================================
  /// Return the flux derivatives as a function of the unknowns
  //=====================================================================
  template<unsigned DIM>
  void ScalarAdvectionEquations<DIM>::dflux_du(const Vector<double>& u,
                                               RankThreeTensor<double>& df_du)
  {
    const unsigned n_flux = this->nflux();

    // Get the wind
    Vector<double> W(DIM);
    Vector<double> s(DIM), x(DIM);
    // Dummy integration point
    unsigned ipt = 0;
    this->get_wind_scalar_adv(ipt, s, x, W);

    df_du.initialise(0.0);

    for (unsigned i = 0; i < n_flux; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        df_du(i, j, i) = W[j];
      }
    }
  }

  template class ScalarAdvectionEquations<1>;
  template class ScalarAdvectionEquations<2>;
  template class ScalarAdvectionEquations<3>;

} // namespace oomph
