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

#include "euler_elements.h"

namespace oomph
{
  //===========================================================
  /// Set the default value of Gamma to be 1.4 in all three
  /// dimension specialisations. This form seems to be required for
  /// gcc 4.4.4, rather than a more general templated version.
  //===========================================================
  template<>
  double EulerEquations<1>::Default_Gamma_Value = 1.4;
  template<>
  double EulerEquations<2>::Default_Gamma_Value = 1.4;
  template<>
  double EulerEquations<3>::Default_Gamma_Value = 1.4;


  //======================================================
  /// Calculate the pressure value from the unknowns
  //=====================================================
  template<unsigned DIM>
  double EulerEquations<DIM>::pressure(const Vector<double>& u) const
  {
    // Initialise the pressure to zero
    double p = 0.0;
    // Subtract off the momentum components
    for (unsigned j = 0; j < DIM; j++)
    {
      p -= u[2 + j] * u[2 + j];
    }
    // Multiply by half and divide by the extra density component
    p *= 0.5 / u[0];
    // Now add on the energy
    p += u[1];
    // Finaly multiply by gamma minus 1
    p *= (this->gamma() - 1);

    // return the pressure
    return p;
  }


  //=========================================================
  /// \short Return the flux as a function of the unknowns
  /// The unknowns are stored as density, energy and then
  /// the velocity components
  //=========================================================
  template<unsigned DIM>
  void EulerEquations<DIM>::flux(const Vector<double>& u,
                                 DenseMatrix<double>& f)
  {
    // The density flux is the momentum
    for (unsigned j = 0; j < DIM; j++)
    {
      f(0, j) = u[2 + j];
    }
    // The energy flux is given by the velocity component multiplied by
    // E + p
    // Find the pressure
    double p = pressure(u);

    // The we can do the energy fluxes
    for (unsigned j = 0; j < DIM; j++)
    {
      f(1, j) = u[2 + j] * (u[1] + p) / u[0];
    }

    // Now the momentum fluxes
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        f(2 + i, j) = u[2 + j] * u[2 + i] / u[0];
      }
      // Add the additional diagonal terms
      f(2 + i, i) += p;
    }
  }


  //====================================================================
  /// Output function, print the values of all unknowns
  //==================================================================
  template<unsigned DIM>
  void EulerEquations<DIM>::output(std::ostream& outfile, const unsigned& nplot)
  {
    // Find the number of fluxes
    const unsigned n_flux = this->nflux();

    // Vector of local coordinates
    Vector<double> s(DIM);
    // Vector of values
    Vector<double> u(n_flux, 0.0);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << this->interpolated_x(s, i) << " ";
      }

      // Values
      for (unsigned i = 0; i < n_flux; i++)
      {
        u[i] = this->interpolated_u_flux_transport(s, i);
        outfile << u[i] << " ";
      }

      // Now output the velocity
      for (unsigned j = 0; j < DIM; j++)
      {
        outfile << u[2 + j] / u[0] << " ";
      }

      // Also the pressure
      outfile << pressure(u);

      outfile << std::endl;
    }
    outfile << std::endl;

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// \short Return the flux derivatives as a function of the unknowns
  //=====================================================================
  /*template<unsigned DIM>
  void EulerEquations<DIM>::
  dflux_du(const Vector<double> &u, RankThreeTensor<double> &df_du)
  {
  }*/

  template class EulerEquations<1>;
  template class EulerEquations<2>;
  template class EulerEquations<3>;

  template<unsigned NNODE_1D>
  Gauss<1, NNODE_1D>
    DGSpectralEulerElement<2, NNODE_1D>::Default_face_integration_scheme;

  template class DGSpectralEulerElement<2, 2>;
  template class DGSpectralEulerElement<2, 3>;
  template class DGSpectralEulerElement<2, 4>;

} // namespace oomph
