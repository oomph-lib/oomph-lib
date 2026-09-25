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
// Non-inline functions for Poisson elements
#include "benchmark_elements.h"


namespace oomph
{
  //======================================================================
  /// Compute element residual Vector and/or element Jacobian matrix
  //======================================================================
  template<unsigned NNODE_1D>
  void ShellBenchmarkElement<NNODE_1D>::
    fill_in_generic_residual_contribution_shell_benchmark(
      Vector<double>& residuals)
  {
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node, 1);

    // Index at which the poisson unknown is stored
    const unsigned u_nodal_index = 0; // u_index_shell_benchmark();

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Integer to store local equation number
    int local_eqn = 0;

    // Get the Poisson ratio
    double nu_local = nu();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_eulerian_at_knot(ipt, psi, dpsidx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of unknown
      // Allocate and initialise to zero
      double interpolated_u = 0.0;
      double interpolated_x = 0.0;
      double interpolated_dudx = 0.0;

      // Calculate function value and derivatives:
      //-----------------------------------------
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the nodal value of the poisson unknown
        double u_value = raw_nodal_value(l, u_nodal_index);
        interpolated_u += u_value * psi(l);
        interpolated_dudx += u_value * dpsidx(l, 0);
        interpolated_x += raw_nodal_position(l, 0) * psi(l);
      }

      // Assign the principal stretches, as well as the J term
      double lambda_r = interpolated_dudx;
      double lambda_t = interpolated_u / interpolated_x;
      double J_vol = lambda_r * lambda_t * lambda_t;
      double tmp = 0.5 / (1.0 + nu_local);
      // Assemble the A term
      double A = tmp *
                 (lambda_r * lambda_r - 1.0 +
                  ((2.0 * nu_local) / (1.0 - 2.0 * nu_local)) * log(J_vol)) /
                 J_vol;
      // Assemble B term
      double B =
        tmp * (2.0 / interpolated_u) * (pow(lambda_r / lambda_t, 2.0) - 1.0);

      // Assemble residuals (Jacobian is done by finite
      // differences)

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the local equation
        local_eqn = nodal_local_eqn(l, u_nodal_index);
        // IF it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Add the contribution
          residuals[local_eqn] += (-1.0 * A * dpsidx(l, 0) + B * psi(l)) * W;
        }
      } // End of loop over test functions
    } // End of loop over integration points
  }


  //======================================================================
  /// Output function
  //======================================================================
  template<unsigned NNODE_1D>
  void ShellBenchmarkElement<NNODE_1D>::output(std::ostream& outfile,
                                               const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(1);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);
      outfile << interpolated_x(s, 0) << " ";
      outfile << interpolated_u_shell_benchmark(s) << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }

  //===========================================================================
  /// Compute the element's residual vector and the (zero) Jacobian matrix.
  //===========================================================================
  void ShellBenchmarkPressureElement::
    fill_in_generic_residual_contribution_shell_benchmark_pressure(
      Vector<double>& residuals)
  {
    // Get the pressure
    double pressure_local = pressure();

    // Get the local equation
    int local_eqn = nodal_local_eqn(0, 0);

    // If a Dirichlet condition has not been set
    if (local_eqn >= 0)
    {
      // Add the contribution
      residuals[local_eqn] -= pressure_local;
    }
  }

  // Instantiate required elements
  template class ShellBenchmarkElement<2>;
  template class ShellBenchmarkElement<3>;
  template class ShellBenchmarkElement<4>;

} // namespace oomph