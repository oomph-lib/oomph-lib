// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
// Non-inline functions for FoepplvonKarman displacement elements

#include "displacement_based_foeppl_von_karman_elements.h"

#include <iostream>

namespace oomph
{
  // Foeppl von Karman displacement equations static data

  /// Default value for Poisson's ratio
  double DisplacementBasedFoepplvonKarmanEquations::Default_Nu_Value = 0.5;

  /// Default value physical constants
  double
    DisplacementBasedFoepplvonKarmanEquations::Default_Physical_Constant_Value =
      0.0;


  //======================================================================
  /// Self-test:  Return 0 for OK
  //======================================================================
  unsigned DisplacementBasedFoepplvonKarmanEquations::self_test()
  {
    bool passed = true;

    // Check lower-level stuff
    if (FiniteElement::self_test() != 0)
    {
      passed = false;
    }

    // Return verdict
    if (passed)
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }

  //======================================================================
  /// Output function:
  ///
  ///   x,y,w
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  void DisplacementBasedFoepplvonKarmanEquations::output(std::ostream& outfile,
                                                         const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);
    // Vector of Eulerian coordinates
    Vector<double> x(2);

    // Sigma_{\alpha \beta}
    DenseMatrix<double> sigma(2, 2);

    // E_{\alpha \beta}
    DenseMatrix<double> strain(2, 2);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get the Eulerian coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        x[i] = interpolated_x(s, i);
        outfile << x[i] << " ";
      }

      outfile << interpolated_w_fvk(s, 0) << " "; // w
      outfile << interpolated_w_fvk(s, 1) << " "; // Laplacian W
      outfile << interpolated_w_fvk(s, 2) << " "; // Ux
      outfile << interpolated_w_fvk(s, 3) << " "; // Uy
      //      outfile << std::endl; // New line

      // Get the stresses at local coordinate s
      this->get_stress_and_strain_for_output(s, sigma, strain);

      // Output stress
      outfile << sigma(0, 0) << " "; // sigma_xx
      outfile << sigma(0, 1) << " "; // sigma_xy
      outfile << sigma(1, 1) << " "; // sigma_yy

      // Output strain
      outfile << strain(0, 0) << " "; // E_xx
      outfile << strain(0, 1) << " "; // E_xy
      outfile << strain(1, 1) << " "; // E_yy
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// C-style output function:
  ///
  ///   x,y,w
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  void DisplacementBasedFoepplvonKarmanEquations::output(FILE* file_pt,
                                                         const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Tecplot header info
    fprintf(file_pt, "%s", tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < 2; i++)
      {
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }
      fprintf(file_pt, "%g \n", interpolated_w_fvk(s));
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, nplot);
  }


  //======================================================================
  /// Output exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  ///
  ///   x,y,w_exact
  //======================================================================
  void DisplacementBasedFoepplvonKarmanEquations::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordintes
    Vector<double> x(2);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Exact solution Vector (here a scalar)
    // Vector<double> exact_soln(1);
    Vector<double> exact_soln(3);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,w_exact
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << " ";
      outfile << exact_soln[1] << " ";
      outfile << exact_soln[2] << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// Validate against exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  ///
  //======================================================================
  void DisplacementBasedFoepplvonKarmanEquations::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    // Initialise
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordintes
    Vector<double> x(2);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Tecplot
    outfile << "ZONE" << std::endl;

    // Exact solution vector
    Vector<double> exact_soln(3);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      interpolated_x(s, x);

      // Get FE function value
      double w_fe = interpolated_w_fvk(s);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,error
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << " " << exact_soln[0] - w_fe << std::endl;

      // Add to error and norm
      norm += exact_soln[0] * exact_soln[0] * W;
      error += (exact_soln[0] - w_fe) * (exact_soln[0] - w_fe) * W;
    }
  }

} // namespace oomph
