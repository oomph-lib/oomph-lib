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
#include "integral_equations.h"

namespace oomph
{
  template<unsigned DIM>
  void IntegralEquations<DIM>::set_output_data_pt(Data* const& data_pt)
  {
    bool use_fd_for_jacobian = true;
    this->add_external_data(data_pt, use_fd_for_jacobian);
  }

  template<unsigned DIM>
  void IntegralEquations<DIM>::setup_integrand(
    Data* const& output_data_pt,
    const IntegralEquations<DIM>::IntegrandFctPt& integrand_fct_pt)
  {
    bool use_fd_for_jacobian = true;
    this->add_external_data(output_data_pt, use_fd_for_jacobian);

    Vector_integrand_fct_pt.push_back(integrand_fct_pt);
  }

  template<unsigned DIM>
  void IntegralEquations<DIM>::fill_in_generic_residual_contribution(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    const unsigned& flag)
  {
    //  if (this->nexternal_data() != 1)
    //  {
    //    throw OomphLibError("Should only be one external data",
    //                        OOMPH_CURRENT_FUNCTION,
    //                        OOMPH_EXCEPTION_LOCATION);
    //  }

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Set up memory for the shape functions
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Get continuous time from timestepper of first node
    const double t = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Integers to store the local equation and unknown numbers
    int local_eqn = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and functions
      double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      Vector<double> interpolated_x(DIM, 0.0);
      // Calculate function value
      //-----------------------------------------
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_x[j] += raw_nodal_position(l, j) * psi(l);
        }
      }

      const unsigned n_integrand = Vector_integrand_fct_pt.size();
      for (unsigned n = 0; n < n_integrand; n++)
      {
        // Get source function
        //-------------------
        double f;
        get_integrand(n, t, interpolated_x, f);

        // Assemble residuals and Jacobian
        //--------------------------------

        unsigned i_value = 0;

        // Get the local equation
        local_eqn = external_local_eqn(n, i_value);

        // Loop over the test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Add body force/source term here
          residuals[local_eqn] += f * psi(l) * W;
        }
      }
    }
  }

  template<unsigned DIM>
  void IntegralEquations<DIM>::output(std::ostream& outfile,
                                      const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Get continuous time from timestepper of first node
    const double t = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      Vector<double> x(DIM);
      double f;
      // Coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        x[i] = interpolated_x(s, i);
        outfile << x[i] << " ";
      }

      const unsigned n_integrand = Vector_integrand_fct_pt.size();
      for (unsigned n = 0; n < n_integrand; n++)
      {
        get_integrand(n, t, x, f);

        // Integrand
        outfile << f << " ";
      }

      outfile << std::endl;
    }
    outfile << std::endl;

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }

  template<unsigned DIM>
  void IntegralEquations<DIM>::output(std::ostream& outfile)
  {
    const unsigned default_n_plot = 5;
    output(outfile, default_n_plot);
  }

  /// Get flux: flux[i] = du/dx_i
  template<unsigned DIM>
  void IntegralEquations<DIM>::get_integral_flux(const Vector<double>& s,
                                                 Vector<double>& flux) const
  {
    // Find out how many nodes there are in the element
    const unsigned n_node = nnode();

    // Vector for coordinates
    Vector<double> x(DIM);

    // Get x position as Vector
    interpolated_x(s, x);

    // Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

    // Call the derivatives of the shape and test functions
    dshape_eulerian(s, psi, dpsidx);

    // Initialise to zero
    for (unsigned j = 0; j < DIM; j++)
    {
      flux[j] = 0.0;
    }

    double f;
    // Loop over nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over derivative directions
      for (unsigned j = 0; j < DIM; j++)
      {
        (*Vector_integrand_fct_pt[0])(0, x, f);
        flux[j] += f * dpsidx(l, j);
      }
    }
  }

  template<unsigned DIM>
  inline void IntegralEquations<DIM>::get_integrand(const unsigned& n,
                                                    const double& t,
                                                    const Vector<double>& x,
                                                    double& f) const
  {
    if (Vector_integrand_fct_pt[n] == 0)
    {
      f = 1.0;
    }
    else
    {
      (*Vector_integrand_fct_pt[n])(t, x, f);
    }
  }

  //====================================================================
  /// / Force build of templates
  //====================================================================
  template class IntegralEquations<1>;
  template class IntegralEquations<2>;
  template class IntegralEquations<3>;
} // namespace oomph
