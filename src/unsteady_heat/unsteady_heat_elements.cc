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
// Non-inline functions for UnsteadyHeat elements
#include "unsteady_heat_elements.h"


namespace oomph
{
  /// 2D UnsteadyHeat elements


  /// Default value for Alpha parameter (thermal inertia)
  template<unsigned DIM>
  double UnsteadyHeatEquations<DIM>::Default_alpha_parameter = 1.0;

  /// Default value for Beta parameter (thermal conductivity)
  template<unsigned DIM>
  double UnsteadyHeatEquations<DIM>::Default_beta_parameter = 1.0;


  //======================================================================
  // Set the data for the number of Variables at each node
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned QUnsteadyHeatElement<DIM, NNODE_1D>::Initial_Nvalue = 1;

  //======================================================================
  /// Compute element residual Vector and/or element Jacobian matrix
  ///
  /// flag=1: compute both
  /// flag=0: compute only residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  template<unsigned DIM>
  void UnsteadyHeatEquations<DIM>::
    fill_in_generic_residual_contribution_ust_heat(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Find the index at which the variable is stored
    unsigned u_nodal_index = u_index_ust_heat();

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get Alpha and beta parameters number
    double alpha_local = alpha();
    double beta_local = beta();

    // Integers to hold the local equation and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++) s[i] = integral_pt()->knot(ipt, i);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_ust_heat(
        ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Allocate memory for local quantities and initialise to zero
      double interpolated_u = 0.0;
      double dudt = 0.0;
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> interpolated_dudx(DIM, 0.0);
      Vector<double> mesh_velocity(DIM, 0.0);

      // Calculate function value and derivatives:
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Calculate the value at the nodes
        double u_value = raw_nodal_value(l, u_nodal_index);
        interpolated_u += u_value * psi(l);
        dudt += du_dt_ust_heat(l) * psi(l);
        // Loop over directions
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_x[j] += raw_nodal_position(l, j) * psi(l);
          interpolated_dudx[j] += u_value * dpsidx(l, j);
        }
      }

      // Mesh velocity?
      if (!ALE_is_disabled)
      {
        for (unsigned l = 0; l < n_node; l++)
        {
          for (unsigned j = 0; j < DIM; j++)
          {
            mesh_velocity[j] += raw_dnodal_position_dt(l, j) * psi(l);
          }
        }
      }

      // Get source function
      //-------------------
      double source = 0.0;
      get_source_ust_heat(time, ipt, interpolated_x, source);

      // Assemble residuals and Jacobian
      //--------------------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        local_eqn = nodal_local_eqn(l, u_nodal_index);
        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          // Add body force/source term and time derivative
          residuals[local_eqn] += (alpha_local * dudt + source) * test(l) * W;

          // The mesh velocity bit
          if (!ALE_is_disabled)
          {
            for (unsigned k = 0; k < DIM; k++)
            {
              residuals[local_eqn] -= alpha_local * mesh_velocity[k] *
                                      interpolated_dudx[k] * test(l) * W;
            }
          }

          // Laplace operator
          for (unsigned k = 0; k < DIM; k++)
          {
            residuals[local_eqn] +=
              beta_local * interpolated_dudx[k] * dtestdx(l, k) * W;
          }


          // Calculate the jacobian
          //-----------------------
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              local_unknown = nodal_local_eqn(l2, u_nodal_index);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                // Mass matrix
                jacobian(local_eqn, local_unknown) +=
                  alpha_local * test(l) * psi(l2) *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * W;

                // Laplace operator & mesh velocity bit
                for (unsigned i = 0; i < DIM; i++)
                {
                  double tmp = beta_local * dtestdx(l, i);
                  if (!ALE_is_disabled)
                    tmp -= alpha_local * mesh_velocity[i] * test(l);
                  jacobian(local_eqn, local_unknown) += dpsidx(l2, i) * tmp * W;
                }
              }
            }
          }
        }
      }


    } // End of loop over integration points
  }


  //======================================================================
  /// Compute norm of fe solution
  //======================================================================
  template<unsigned DIM>
  void UnsteadyHeatEquations<DIM>::compute_norm(double& norm)
  {
    // Initialise
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get FE function value
      double u = interpolated_u_ust_heat(s);

      // Add to  norm
      norm += u * u * W;
    }
  }

  //======================================================================
  /// Self-test:  Return 0 for OK
  //======================================================================
  template<unsigned DIM>
  unsigned UnsteadyHeatEquations<DIM>::self_test()
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
  ///   x,y,u   or    x,y,z,u
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  template<unsigned DIM>
  void UnsteadyHeatEquations<DIM>::output(std::ostream& outfile,
                                          const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }
      outfile << interpolated_u_ust_heat(s) << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// C-style output function:
  ///
  ///   x,y,u   or    x,y,z,u
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  template<unsigned DIM>
  void UnsteadyHeatEquations<DIM>::output(FILE* file_pt, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    fprintf(file_pt, "%s", tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < DIM; i++)
      {
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }
      fprintf(file_pt, "%g \n", interpolated_u_ust_heat(s));
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
  ///   x,y,u_exact    or    x,y,z,u_exact
  //======================================================================
  template<unsigned DIM>
  void UnsteadyHeatEquations<DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

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
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// Output exact solution at time t
  ///
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  ///
  ///   x,y,u_exact    or    x,y,z,u_exact
  //======================================================================
  template<unsigned DIM>
  void UnsteadyHeatEquations<DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)

  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);


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
      (*exact_soln_pt)(time, x, exact_soln);

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < DIM; i++)
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
  template<unsigned DIM>
  void UnsteadyHeatEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    // Initialise
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

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
      for (unsigned i = 0; i < DIM; i++)
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
      double u_fe = interpolated_u_ust_heat(s);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,...,error
      for (unsigned i = 0; i < DIM; i++)
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
  /// Validate against exact solution at time t.
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  ///
  //======================================================================
  template<unsigned DIM>
  void UnsteadyHeatEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
    double& error,
    double& norm)

  {
    // Initialise
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);


    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Tecplot
    outfile << "ZONE" << std::endl;

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
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
      double u_fe = interpolated_u_ust_heat(s);

      // Get exact solution at this point
      (*exact_soln_pt)(time, x, exact_soln);

      // Output x,y,...,error
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << " " << exact_soln[0] - u_fe << std::endl;

      // Add to error and norm
      norm += exact_soln[0] * exact_soln[0] * W;
      error += (exact_soln[0] - u_fe) * (exact_soln[0] - u_fe) * W;
    }
  }


  //====================================================================
  // Force build of templates
  //====================================================================
  template class QUnsteadyHeatElement<1, 2>;
  template class QUnsteadyHeatElement<1, 3>;
  template class QUnsteadyHeatElement<1, 4>;

  template class QUnsteadyHeatElement<2, 2>;
  template class QUnsteadyHeatElement<2, 3>;
  template class QUnsteadyHeatElement<2, 4>;

  template class QUnsteadyHeatElement<3, 2>;
  template class QUnsteadyHeatElement<3, 3>;
  template class QUnsteadyHeatElement<3, 4>;


} // namespace oomph
