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
// Non-inline functions for axisym FoepplvonKarman elements

#include "axisym_fvk_elements.h"

namespace oomph
{
  //======================================================================
  /// Set the data for the number of Variables at each node: 6
  //======================================================================
  template<unsigned NNODE_1D>
  const unsigned AxisymFoepplvonKarmanElement<NNODE_1D>::Initial_Nvalue = 6;


  //======================================================================
  /// Default value physical constants
  //======================================================================
  double AxisymFoepplvonKarmanEquations::Default_Physical_Constant_Value = 0.0;


  //======================================================================
  /// Compute contribution to element residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  void AxisymFoepplvonKarmanEquations::fill_in_contribution_to_residuals(
    Vector<double>& residuals)
  {
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidr(n_node, 1), dtestdr(n_node, 1);

    // Indices at which the unknowns are stored
    const unsigned w_nodal_index = nodal_index_fvk(0);
    const unsigned laplacian_w_nodal_index = nodal_index_fvk(1);
    const unsigned phi_nodal_index = nodal_index_fvk(2);
    const unsigned laplacian_phi_nodal_index = nodal_index_fvk(3);
    const unsigned smooth_dwdr_nodal_index = nodal_index_fvk(4);
    const unsigned smooth_dphidr_nodal_index = nodal_index_fvk(5);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Integers to store the local equation numbers
    int local_eqn = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_axisym_fvk(
        ipt, psi, dpsidr, test, dtestdr);

      // Allocate and initialise to zero storage for the interpolated values
      double interpolated_r = 0;

      double interpolated_w = 0;
      double interpolated_laplacian_w = 0;
      double interpolated_phi = 0;
      double interpolated_laplacian_phi = 0;

      double interpolated_dwdr = 0;
      double interpolated_dlaplacian_wdr = 0;
      double interpolated_dphidr = 0;
      double interpolated_dlaplacian_phidr = 0;

      double interpolated_smooth_dwdr = 0;
      double interpolated_smooth_dphidr = 0;
      double interpolated_continuous_d2wdr2 = 0;
      double interpolated_continuous_d2phidr2 = 0;

      // Calculate function values and derivatives:
      //-----------------------------------------
      Vector<double> nodal_value(6, 0.0);
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the nodal values
        nodal_value[0] = raw_nodal_value(l, w_nodal_index);
        nodal_value[1] = raw_nodal_value(l, laplacian_w_nodal_index);

        if (!Linear_bending_model)
        {
          nodal_value[2] = raw_nodal_value(l, phi_nodal_index);
          nodal_value[3] = raw_nodal_value(l, laplacian_phi_nodal_index);
          nodal_value[4] = raw_nodal_value(l, smooth_dwdr_nodal_index);
          nodal_value[5] = raw_nodal_value(l, smooth_dphidr_nodal_index);
        }

        // Add contributions from current node/shape function
        interpolated_w += nodal_value[0] * psi(l);
        interpolated_laplacian_w += nodal_value[1] * psi(l);

        if (!Linear_bending_model)
        {
          interpolated_phi += nodal_value[2] * psi(l);
          interpolated_laplacian_phi += nodal_value[3] * psi(l);
          interpolated_smooth_dwdr += nodal_value[4] * psi(l);
          interpolated_smooth_dphidr += nodal_value[5] * psi(l);

          interpolated_continuous_d2wdr2 += nodal_value[4] * dpsidr(l, 0);
          interpolated_continuous_d2phidr2 += nodal_value[5] * dpsidr(l, 0);
        }

        interpolated_r += raw_nodal_position(l, 0) * psi(l);
        interpolated_dwdr += nodal_value[0] * dpsidr(l, 0);
        interpolated_dlaplacian_wdr += nodal_value[1] * dpsidr(l, 0);

        if (!Linear_bending_model)
        {
          interpolated_dphidr += nodal_value[2] * dpsidr(l, 0);
          interpolated_dlaplacian_phidr += nodal_value[3] * dpsidr(l, 0);
        }
      } // End of loop over the nodes


      // Premultiply the weights and the Jacobian
      double W = w * interpolated_r * J;

      // Get pressure function
      //-------------------
      double pressure;
      get_pressure_fvk(ipt, interpolated_r, pressure);

      double airy_forcing;
      get_airy_forcing_fvk(ipt, interpolated_r, airy_forcing);


      // Assemble residuals and Jacobian
      //--------------------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the local equation
        local_eqn = nodal_local_eqn(l, w_nodal_index);

        // IF it's not a boundary condition
        if (local_eqn >= 0)
        {
          residuals[local_eqn] += pressure * test(l) * W;

          // Reduced order biharmonic operator
          residuals[local_eqn] +=
            interpolated_dlaplacian_wdr * dtestdr(l, 0) * W;

          if (!Linear_bending_model)
          {
            // Monge-Ampere part
            residuals[local_eqn] +=
              eta() *
              (interpolated_continuous_d2wdr2 * interpolated_dphidr +
               interpolated_dwdr * interpolated_continuous_d2phidr2) /
              interpolated_r * test(l) * W;
          }
        }

        // Get the local equation
        local_eqn = nodal_local_eqn(l, laplacian_w_nodal_index);

        // IF it's not a boundary condition
        if (local_eqn >= 0)
        {
          // The coupled Poisson equations for the biharmonic operator
          residuals[local_eqn] += interpolated_laplacian_w * test(l) * W;
          residuals[local_eqn] += interpolated_dwdr * dtestdr(l, 0) * W;
        }

        // Get the local equation
        local_eqn = nodal_local_eqn(l, phi_nodal_index);

        // IF it's not a boundary condition
        if (local_eqn >= 0)
        {
          residuals[local_eqn] += airy_forcing * test(l) * W;

          // Reduced order biharmonic operator
          residuals[local_eqn] +=
            interpolated_dlaplacian_phidr * dtestdr(l, 0) * W;

          // Monge-Ampere part
          residuals[local_eqn] -= interpolated_continuous_d2wdr2 *
                                  interpolated_dwdr / interpolated_r * test(l) *
                                  W;
        }

        // Get the local equation
        local_eqn = nodal_local_eqn(l, laplacian_phi_nodal_index);

        // IF it's not a boundary condition
        if (local_eqn >= 0)
        {
          // The coupled Poisson equations for the biharmonic operator
          residuals[local_eqn] += interpolated_laplacian_phi * test(l) * W;
          residuals[local_eqn] += interpolated_dphidr * dtestdr(l, 0) * W;
        }

        // Residuals for the smooth derivatives
        local_eqn = nodal_local_eqn(l, smooth_dwdr_nodal_index);
        if (local_eqn >= 0)
        {
          residuals[local_eqn] +=
            (interpolated_dwdr - interpolated_smooth_dwdr) * test(l) * W;
        }

        local_eqn = nodal_local_eqn(l, smooth_dphidr_nodal_index);
        if (local_eqn >= 0)
        {
          residuals[local_eqn] +=
            (interpolated_dphidr - interpolated_smooth_dphidr) * test(l) * W;
        }

      } // End of loop over test functions
    } // End of loop over integration points
  }


  //======================================================================
  /// Self-test:  Return 0 for OK
  //======================================================================
  unsigned AxisymFoepplvonKarmanEquations::self_test()
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
  /// Compute in-plane stresses. Return boolean to indicate success
  /// (false if attempt to evaluate stresses at zero radius)
  //======================================================================
  bool AxisymFoepplvonKarmanEquations::interpolated_stress(
    const Vector<double>& s, double& sigma_r_r, double& sigma_phi_phi)

  {
    // No in plane stresses if linear bending
    if (Linear_bending_model)
    {
      sigma_r_r = 0.0;
      sigma_phi_phi = 0.0;

      // Success!
      return true;
    }
    else
    {
      // Get shape fcts and derivs
      unsigned n_dim = this->dim();
      unsigned n_node = this->nnode();
      Shape psi(n_node);
      DShape dpsi_dr(n_node, n_dim);

      // Check if we're dividing by zero
      Vector<double> r(1);
      this->interpolated_x(s, r);
      if (r[0] == 0.0)
      {
        sigma_r_r = 0.0;
        sigma_phi_phi = 0.0;
        return false;
      }

      // Get shape fcts and derivs
      dshape_eulerian(s, psi, dpsi_dr);

      // Indices at which the unknowns are stored
      const unsigned smooth_dphi_dr_nodal_index = nodal_index_fvk(5);

      // Allocate and initialise to zero storage for the interpolated values
      double interpolated_r = 0;
      double interpolated_dphi_dr = 0;
      double interpolated_continuous_d2phi_dr2 = 0;


      // Calculate function values and derivatives:
      //-----------------------------------------
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Add contributions from current node/shape function
        interpolated_r += raw_nodal_position(l, 0) * psi(l);
        interpolated_dphi_dr +=
          this->nodal_value(l, smooth_dphi_dr_nodal_index) * psi(l);
        interpolated_continuous_d2phi_dr2 +=
          this->nodal_value(l, smooth_dphi_dr_nodal_index) * dpsi_dr(l, 0);
      } // End of loop over nodes

      // Compute stresses
      sigma_r_r = interpolated_dphi_dr / interpolated_r;
      sigma_phi_phi = interpolated_continuous_d2phi_dr2;

      // Success!
      return true;

    } // End if

  } // End of interpolated_stress function


  //======================================================================
  /// Output function:
  ///   r, w, sigma_r_r, sigma_phi_phi
  /// nplot points
  //======================================================================
  void AxisymFoepplvonKarmanEquations::output(std::ostream& outfile,
                                              const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(1);

    // Tecplot header info
    outfile << "ZONE\n";

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get stress
      double sigma_r_r = 0.0;
      double sigma_phi_phi = 0.0;
      bool success = interpolated_stress(s, sigma_r_r, sigma_phi_phi);
      if (success)
      {
        outfile << interpolated_x(s, 0) << " " << interpolated_w_fvk(s) << " "
                << sigma_r_r << " " << sigma_phi_phi << std::endl;
      }
    }
  }


  //======================================================================
  /// C-style output function:
  ///   r,w
  /// nplot points
  //======================================================================
  void AxisymFoepplvonKarmanEquations::output(FILE* file_pt,
                                              const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(1);

    // Tecplot header info
    fprintf(file_pt, "%s", tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      fprintf(file_pt, "%g ", interpolated_x(s, 0));
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
  ///   r,w_exact
  //======================================================================
  void AxisymFoepplvonKarmanEquations::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(1);

    // Vector for coordinates
    Vector<double> r(1);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Exact solution Vector (here a scalar)
    // Vector<double> exact_soln(1);
    Vector<double> exact_soln(1);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get r position as Vector
      interpolated_x(s, r);

      // Get exact solution at this point
      (*exact_soln_pt)(r, exact_soln);

      // Output r,w_exact
      outfile << r[0] << " ";
      outfile << exact_soln[0] << std::endl;
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
  void AxisymFoepplvonKarmanEquations::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    // Initialise
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(1);

    // Vector for coordintes
    Vector<double> r(1);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Tecplot
    outfile << "ZONE" << std::endl;

    // Exact solution Vector (here a scalar)
    // Vector<double> exact_soln(1);
    Vector<double> exact_soln(1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      s[0] = integral_pt()->knot(ipt, 0);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get r position as Vector
      interpolated_x(s, r);

      // Get FE function value
      double w_fe = interpolated_w_fvk(s);

      // Get exact solution at this point
      (*exact_soln_pt)(r, exact_soln);

      // Output r,error
      outfile << r[0] << " ";
      outfile << exact_soln[0] << " " << exact_soln[0] - w_fe << std::endl;

      // Add to error and norm
      norm += exact_soln[0] * exact_soln[0] * W;
      error += (exact_soln[0] - w_fe) * (exact_soln[0] - w_fe) * W;
    }

    {
      // Initialise
      error = 0.0;
      norm = 0.0;

      // Vector of local coordinates
      Vector<double> s(1);

      // Vector for coordintes
      Vector<double> r(1);

      // Find out how many nodes there are in the element
      unsigned n_node = nnode();

      Shape psi(n_node);

      // Set the value of n_intpt
      unsigned n_intpt = integral_pt()->nweight();

      // Tecplot
      outfile << "ZONE" << std::endl;

      // Exact solution Vector (here a scalar)
      // Vector<double> exact_soln(1);
      Vector<double> exact_soln(1);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s
        s[0] = integral_pt()->knot(ipt, 0);

        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Get jacobian of mapping
        double J = J_eulerian(s);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get r position as Vector
        interpolated_x(s, r);

        // Get FE function value
        double w_fe = interpolated_w_fvk(s);

        // Get exact solution at this point
        (*exact_soln_pt)(r, exact_soln);

        // Output r error
        outfile << r[0] << " ";
        outfile << exact_soln[0] << " " << exact_soln[0] - w_fe << std::endl;

        // Add to error and norm
        norm += exact_soln[0] * exact_soln[0] * W;
        error += (exact_soln[0] - w_fe) * (exact_soln[0] - w_fe) * W;
      }
    }
  }


  //====================================================================
  // Force build of templates
  //====================================================================
  template class AxisymFoepplvonKarmanElement<2>;
  template class AxisymFoepplvonKarmanElement<3>;
  template class AxisymFoepplvonKarmanElement<4>;

} // namespace oomph
