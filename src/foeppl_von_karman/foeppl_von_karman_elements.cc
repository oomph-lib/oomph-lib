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
// Non-inline functions for FoepplvonKarman elements

#include "foeppl_von_karman_elements.h"

#include <iostream>

namespace oomph
{
  //======================================================================
  /// Set the data for the number of Variables at each node - 8
  //======================================================================
  template<unsigned NNODE_1D>
  const unsigned QFoepplvonKarmanElement<NNODE_1D>::Initial_Nvalue = 8;

  // Foeppl von Karman equations static data

  /// Default value physical constants
  double FoepplvonKarmanEquations::Default_Physical_Constant_Value = 0.0;

  //======================================================================
  /// Compute contribution to element residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  void FoepplvonKarmanEquations::fill_in_contribution_to_residuals(
    Vector<double>& residuals)
  {
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, 2), dtestdx(n_node, 2);

    // Indices at which the unknowns are stored
    const unsigned w_nodal_index = nodal_index_fvk(0);
    const unsigned laplacian_w_nodal_index = nodal_index_fvk(1);
    const unsigned phi_nodal_index = nodal_index_fvk(2);
    const unsigned laplacian_phi_nodal_index = nodal_index_fvk(3);
    const unsigned smooth_dwdx_nodal_index = nodal_index_fvk(4);
    const unsigned smooth_dwdy_nodal_index = nodal_index_fvk(5);
    const unsigned smooth_dphidx_nodal_index = nodal_index_fvk(6);
    const unsigned smooth_dphidy_nodal_index = nodal_index_fvk(7);

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
      double J =
        dshape_and_dtest_eulerian_at_knot_fvk(ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Allocate and initialise to zero storage for the interpolated values
      Vector<double> interpolated_x(2, 0.0);

      double interpolated_w = 0;
      double interpolated_laplacian_w = 0;
      double interpolated_phi = 0;
      double interpolated_laplacian_phi = 0;

      Vector<double> interpolated_dwdx(2, 0.0);
      Vector<double> interpolated_dlaplacian_wdx(2, 0.0);
      Vector<double> interpolated_dphidx(2, 0.0);
      Vector<double> interpolated_dlaplacian_phidx(2, 0.0);

      Vector<double> interpolated_smooth_dwdx(2, 0.0);
      Vector<double> interpolated_smooth_dphidx(2, 0.0);
      double interpolated_continuous_d2wdx2 = 0;
      double interpolated_continuous_d2wdy2 = 0;
      double interpolated_continuous_d2phidx2 = 0;
      double interpolated_continuous_d2phidy2 = 0;
      double interpolated_continuous_d2wdxdy = 0;
      double interpolated_continuous_d2phidxdy = 0;

      // Calculate function values and derivatives:
      //-----------------------------------------
      Vector<double> nodal_value(8, 0.0);
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
          nodal_value[4] = raw_nodal_value(l, smooth_dwdx_nodal_index);
          nodal_value[5] = raw_nodal_value(l, smooth_dwdy_nodal_index);
          nodal_value[6] = raw_nodal_value(l, smooth_dphidx_nodal_index);
          nodal_value[7] = raw_nodal_value(l, smooth_dphidy_nodal_index);
        }

        // Add contributions from current node/shape function
        interpolated_w += nodal_value[0] * psi(l);
        interpolated_laplacian_w += nodal_value[1] * psi(l);

        if (!Linear_bending_model)
        {
          interpolated_phi += nodal_value[2] * psi(l);
          interpolated_laplacian_phi += nodal_value[3] * psi(l);

          interpolated_smooth_dwdx[0] += nodal_value[4] * psi(l);
          interpolated_smooth_dwdx[1] += nodal_value[5] * psi(l);
          interpolated_smooth_dphidx[0] += nodal_value[6] * psi(l);
          interpolated_smooth_dphidx[1] += nodal_value[7] * psi(l);

          interpolated_continuous_d2wdx2 += nodal_value[4] * dpsidx(l, 0);
          interpolated_continuous_d2wdy2 += nodal_value[5] * dpsidx(l, 1);
          interpolated_continuous_d2phidx2 += nodal_value[6] * dpsidx(l, 0);
          interpolated_continuous_d2phidy2 += nodal_value[7] * dpsidx(l, 1);
          // mjr CHECK THESE
          interpolated_continuous_d2wdxdy +=
            0.5 *
            (nodal_value[4] * dpsidx(l, 1) + nodal_value[5] * dpsidx(l, 0));
          interpolated_continuous_d2phidxdy +=
            0.5 *
            (nodal_value[6] * dpsidx(l, 1) + nodal_value[7] * dpsidx(l, 0));
        }

        // Loop over directions
        for (unsigned j = 0; j < 2; j++)
        {
          interpolated_x[j] += raw_nodal_position(l, j) * psi(l);
          interpolated_dwdx[j] += nodal_value[0] * dpsidx(l, j);
          interpolated_dlaplacian_wdx[j] += nodal_value[1] * dpsidx(l, j);

          if (!Linear_bending_model)
          {
            interpolated_dphidx[j] += nodal_value[2] * dpsidx(l, j);
            interpolated_dlaplacian_phidx[j] += nodal_value[3] * dpsidx(l, j);
          }
        }
      }

      // Get pressure function
      //-------------------
      double pressure;
      get_pressure_fvk(ipt, interpolated_x, pressure);

      double airy_forcing;
      get_airy_forcing_fvk(ipt, interpolated_x, airy_forcing);

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
          for (unsigned k = 0; k < 2; k++)
          {
            residuals[local_eqn] +=
              interpolated_dlaplacian_wdx[k] * dtestdx(l, k) * W;
          }

          if (!Linear_bending_model)
          {
            // Monge-Ampere part
            residuals[local_eqn] += eta() *
                                    (interpolated_continuous_d2wdx2 *
                                       interpolated_continuous_d2phidy2 +
                                     interpolated_continuous_d2wdy2 *
                                       interpolated_continuous_d2phidx2 -
                                     2.0 * interpolated_continuous_d2wdxdy *
                                       interpolated_continuous_d2phidxdy) *
                                    test(l) * W;
          }
        }

        // Get the local equation
        local_eqn = nodal_local_eqn(l, laplacian_w_nodal_index);

        // IF it's not a boundary condition
        if (local_eqn >= 0)
        {
          // The coupled Poisson equations for the biharmonic operator
          residuals[local_eqn] += interpolated_laplacian_w * test(l) * W;

          for (unsigned k = 0; k < 2; k++)
          {
            residuals[local_eqn] += interpolated_dwdx[k] * dtestdx(l, k) * W;
          }
        }

        // Get the local equation
        local_eqn = nodal_local_eqn(l, phi_nodal_index);

        // IF it's not a boundary condition
        if (local_eqn >= 0)
        {
          residuals[local_eqn] += airy_forcing * test(l) * W;

          // Reduced order biharmonic operator
          for (unsigned k = 0; k < 2; k++)
          {
            residuals[local_eqn] +=
              interpolated_dlaplacian_phidx[k] * dtestdx(l, k) * W;
          }

          // Monge-Ampere part
          residuals[local_eqn] -=
            (interpolated_continuous_d2wdx2 * interpolated_continuous_d2wdy2 -
             interpolated_continuous_d2wdxdy *
               interpolated_continuous_d2wdxdy) *
            test(l) * W;
        }

        // Get the local equation
        local_eqn = nodal_local_eqn(l, laplacian_phi_nodal_index);

        // IF it's not a boundary condition
        if (local_eqn >= 0)
        {
          // The coupled Poisson equations for the biharmonic operator
          residuals[local_eqn] += interpolated_laplacian_phi * test(l) * W;

          for (unsigned k = 0; k < 2; k++)
          {
            residuals[local_eqn] += interpolated_dphidx[k] * dtestdx(l, k) * W;
          }
        }

        // Residuals for the smooth derivatives
        local_eqn = nodal_local_eqn(l, smooth_dwdx_nodal_index);

        if (local_eqn >= 0)
        {
          residuals[local_eqn] +=
            (interpolated_dwdx[0] - interpolated_smooth_dwdx[0]) * test(l) * W;
        }

        local_eqn = nodal_local_eqn(l, smooth_dwdy_nodal_index);

        if (local_eqn >= 0)
        {
          residuals[local_eqn] +=
            (interpolated_dwdx[1] - interpolated_smooth_dwdx[1]) * test(l) * W;
        }

        local_eqn = nodal_local_eqn(l, smooth_dphidx_nodal_index);

        if (local_eqn >= 0)
        {
          residuals[local_eqn] +=
            (interpolated_dphidx[0] - interpolated_smooth_dphidx[0]) * test(l) *
            W;
        }

        local_eqn = nodal_local_eqn(l, smooth_dphidy_nodal_index);

        if (local_eqn >= 0)
        {
          residuals[local_eqn] +=
            (interpolated_dphidx[1] - interpolated_smooth_dphidx[1]) * test(l) *
            W;
        }
      }

    } // End of loop over integration points


    // Finally: My contribution to the volume constraint equation
    // (if any). Note this must call get_bounded_volume since the
    // definition of the bounded volume can be overloaded in derived
    // elements.
    if (Volume_constraint_pressure_external_data_index >= 0)
    {
      local_eqn =
        external_local_eqn(Volume_constraint_pressure_external_data_index, 0);
      if (local_eqn >= 0)
      {
        residuals[local_eqn] += get_bounded_volume();
      }
    }
  }

  /*
  void FoepplvonKarmanEquations::fill_in_contribution_to_jacobian(Vector<double>
  &residuals, DenseMatrix<double> &jacobian)
  {
   //Add the contribution to the residuals
   FoepplvonKarmanEquations::fill_in_contribution_to_residuals(residuals);
   //Allocate storage for the full residuals (residuals of entire element)
   unsigned n_dof = ndof();
   Vector<double> full_residuals(n_dof);
   //Get the residuals for the entire element
   FoepplvonKarmanEquations::get_residuals(full_residuals);
   //Calculate the contributions from the internal dofs
   //(finite-difference the lot by default)
   fill_in_jacobian_from_internal_by_fd(full_residuals,jacobian,true);
   //Calculate the contributions from the external dofs
   //(finite-difference the lot by default)
   fill_in_jacobian_from_external_by_fd(full_residuals,jacobian,true);
   //Calculate the contributions from the nodal dofs
   fill_in_jacobian_from_nodal_by_fd(full_residuals,jacobian);
  }
  */


  //======================================================================
  /// Self-test:  Return 0 for OK
  //======================================================================
  unsigned FoepplvonKarmanEquations::self_test()
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
  /// Compute in-plane stresses
  //======================================================================
  void FoepplvonKarmanEquations::interpolated_stress(const Vector<double>& s,
                                                     double& sigma_xx,
                                                     double& sigma_yy,
                                                     double& sigma_xy)
  {
    // No in plane stresses if linear bending
    if (Linear_bending_model)
    {
      sigma_xx = 0.0;
      sigma_yy = 0.0;
      sigma_xy = 0.0;
      return;
    }

    // Get shape fcts and derivs
    unsigned n_dim = this->dim();
    unsigned n_node = this->nnode();
    Shape psi(n_node);
    DShape dpsidx(n_node, n_dim);
    dshape_eulerian(s, psi, dpsidx);
    double interpolated_continuous_d2phidx2 = 0;
    double interpolated_continuous_d2phidy2 = 0;
    double interpolated_continuous_d2phidxdy = 0;

    const unsigned smooth_dphidx_nodal_index = nodal_index_fvk(6);
    const unsigned smooth_dphidy_nodal_index = nodal_index_fvk(7);

    // Loop over nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      interpolated_continuous_d2phidx2 +=
        raw_nodal_value(l, smooth_dphidx_nodal_index) * dpsidx(l, 0);
      interpolated_continuous_d2phidy2 +=
        raw_nodal_value(l, smooth_dphidy_nodal_index) * dpsidx(l, 1);
      interpolated_continuous_d2phidxdy +=
        0.5 * (raw_nodal_value(l, smooth_dphidx_nodal_index) * dpsidx(l, 1) +
               raw_nodal_value(l, smooth_dphidy_nodal_index) * dpsidx(l, 0));
    }

    sigma_xx = interpolated_continuous_d2phidy2;
    sigma_yy = interpolated_continuous_d2phidx2;
    sigma_xy = -interpolated_continuous_d2phidxdy;
  }


  //======================================================================
  /// Output function:
  ///
  ///   x,y,w
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  void FoepplvonKarmanEquations::output(std::ostream& outfile,
                                        const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }
      outfile << interpolated_w_fvk(s) << std::endl;
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
  void FoepplvonKarmanEquations::output(FILE* file_pt, const unsigned& nplot)
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
  void FoepplvonKarmanEquations::output_fct(
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
    Vector<double> exact_soln(2);

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
  void FoepplvonKarmanEquations::compute_error(
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

    // Exact solution Vector (here a scalar)
    // Vector<double> exact_soln(1);
    Vector<double> exact_soln(2);

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


  //====================================================================
  // Force build of templates
  //====================================================================
  template class QFoepplvonKarmanElement<2>;
  template class QFoepplvonKarmanElement<3>;
  template class QFoepplvonKarmanElement<4>;

} // namespace oomph
