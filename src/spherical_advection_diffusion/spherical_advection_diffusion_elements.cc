// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
// Non-inline functions for Advection Diffusion elements in a spherical
// coordinate system
#include "spherical_advection_diffusion_elements.h"

namespace oomph
{
  /// 2D Steady Axisymmetric Advection Diffusion elements

  /// Default value for Peclet number
  double SphericalAdvectionDiffusionEquations::Default_peclet_number = 0.0;


  //======================================================================
  /// Compute element residual Vector and/or element Jacobian matrix
  ///
  /// flag=1: compute both
  /// flag=0: compute only residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  void SphericalAdvectionDiffusionEquations::
    fill_in_generic_residual_contribution_spherical_adv_diff(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Get the nodal index at which the unknown is stored
    const unsigned u_nodal_index = u_index_spherical_adv_diff();

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, 2), dtestdx(n_node, 2);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Get Physical Variables from Element
    const double scaled_peclet = pe();

    // Get peclet number multiplied by Strouhal number
    const double scaled_peclet_st = pe_st();

    // Integers used to store the local equation number and local unknown
    // indices for the residuals and jacobians
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++) s[i] = integral_pt()->knot(ipt, i);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_spherical_adv_diff(
        ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the solution and its derivatives
      // Allocate
      double interpolated_u = 0.0;
      double dudt = 0.0;
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> interpolated_dudx(2, 0.0);

      // Calculate function value and derivatives:
      //-----------------------------------------
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the value at the node
        double u_value = raw_nodal_value(l, u_nodal_index);
        interpolated_u += u_value * psi(l);
        dudt += du_dt_spherical_adv_diff(l) * psi(l);
        // Loop over the two coordinates directions
        for (unsigned j = 0; j < 2; j++)
        {
          interpolated_x[j] += raw_nodal_position(l, j) * psi(l);
          interpolated_dudx[j] += u_value * dpsidx(l, j);
        }
      }

      // Get source function
      //-------------------
      double source;
      get_source_spherical_adv_diff(ipt, interpolated_x, source);


      // Get wind
      //--------
      Vector<double> wind(3);
      get_wind_spherical_adv_diff(ipt, s, interpolated_x, wind);

      // r is the first position component
      double r = interpolated_x[0];
      // theta is the second position component
      double sin_th = sin(interpolated_x[1]);
      // dS is the area weighting
      double dS = r * r * sin_th;

      // TEMPERATURE EQUATION (Neglected viscous dissipation)
      // Assemble residuals and Jacobian
      //--------------------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Set the local equation number
        local_eqn = nodal_local_eqn(l, u_nodal_index);

        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          // Add body force/source term
          residuals[local_eqn] -=
            (scaled_peclet_st * dudt + source) * dS * test(l) * W;

          // The Advection Diffusion bit itself
          residuals[local_eqn] -=
            // radial terms
            (dS * interpolated_dudx[0] *
               (scaled_peclet * wind[0] * test(l) + dtestdx(l, 0)) +
             // azimuthal terms
             (sin_th * interpolated_dudx[1] *
              (r * scaled_peclet * wind[1] * test(l) + dtestdx(l, 1)))) *
            W;


          // Calculate the jacobian
          //-----------------------
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Set the number of the unknown
              local_unknown = nodal_local_eqn(l2, u_nodal_index);

              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                // Mass matrix term
                jacobian(local_eqn, local_unknown) -=
                  scaled_peclet_st * test(l) * psi(l2) *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * dS * W;

                // add the mass matrix term
                if (flag == 2)
                {
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_peclet_st * test(l) * psi(l2) * dS * W;
                }

                // Assemble Jacobian term
                jacobian(local_eqn, local_unknown) -=
                  // radial terms
                  (dS * dpsidx(l2, 0) *
                     (scaled_peclet * wind[0] * test(l) + dtestdx(l, 0)) +
                   // azimuthal terms
                   (sin_th * dpsidx(l2, 1) *
                    (r * scaled_peclet * wind[1] * test(l) + dtestdx(l, 1)))) *
                  W;
              }
            }
          }
        }
      }


    } // End of loop over integration points
  }


  //======================================================================
  /// Self-test:  Return 0 for OK
  //======================================================================
  // template <unsigned DIM>
  unsigned SphericalAdvectionDiffusionEquations::self_test()
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
  ///   r,z,u,w_r,w_z
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  void SphericalAdvectionDiffusionEquations::output(std::ostream& outfile,
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

      // Get Eulerian coordinate of plot point
      Vector<double> x(2);
      interpolated_x(s, x);

      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Convert to Cartesians
      outfile << x[0] * sin(x[1]) << " " << x[0] * cos(x[1]) << " ";
      // Output concentration
      outfile << this->interpolated_u_spherical_adv_diff(s) << " ";

      // Get the wind
      Vector<double> wind(2);
      // Dummy ipt argument needed... ?
      unsigned ipt = 0;
      get_wind_spherical_adv_diff(ipt, s, x, wind);
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << wind[i] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// C-style output function:
  ///
  ///   r,z,u
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  // template <unsigned DIM>
  void SphericalAdvectionDiffusionEquations::output(FILE* file_pt,
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
      // Write Cartesian coordinates
      double r = interpolated_x(s, 0);
      double theta = interpolated_x(s, 1);
      fprintf(file_pt, "%g %g ", r * sin(theta), r * cos(theta));

      fprintf(file_pt, "%g \n", interpolated_u_spherical_adv_diff(s));
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, nplot);
  }

  //======================================================================
  ///  Output exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  ///
  ///   r,z,u_exact
  //======================================================================
  // template <unsigned DIM>
  void SphericalAdvectionDiffusionEquations::output_fct(
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
    Vector<double> exact_soln(1);

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

      // Output x,y,...,u_exact
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
  // template <unsigned DIM>
  void SphericalAdvectionDiffusionEquations::compute_error(
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

    // Tecplot header info
    outfile << "ZONE" << std::endl;

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(1);

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
      double u_fe = interpolated_u_spherical_adv_diff(s);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,...,error
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << " " << exact_soln[0] - u_fe << std::endl;

      // Add to error and norm
      norm += exact_soln[0] * exact_soln[0] * W;
      error += (exact_soln[0] - u_fe) * (exact_soln[0] - u_fe) * W;
    }
  }


  //======================================================================
  // Set the data for the number of Variables at each node
  //======================================================================
  template<unsigned NNODE_1D>
  const unsigned QSphericalAdvectionDiffusionElement<NNODE_1D>::Initial_Nvalue =
    1;

  //====================================================================
  // Force build of templates
  //====================================================================

  template class QSphericalAdvectionDiffusionElement<2>;
  template class QSphericalAdvectionDiffusionElement<3>;
  template class QSphericalAdvectionDiffusionElement<4>;


} // namespace oomph
