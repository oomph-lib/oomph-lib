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
// Non-inline functions for Advection Diffusion elements
#include "advection_diffusion_reaction_elements.h"

namespace oomph
{
  // Specify the number of reagents
  template<unsigned NREAGENT, unsigned DIM>
  const unsigned AdvectionDiffusionReactionEquations<NREAGENT, DIM>::N_reagent =
    NREAGENT;

  /// 2D Advection Diffusion elements


  /// Default value for Peclet number
  template<unsigned NREAGENT, unsigned DIM>
  Vector<double> AdvectionDiffusionReactionEquations<NREAGENT, DIM>::
    Default_dimensionless_number(NREAGENT, 1.0);

  //======================================================================
  ///  Compute element residual Vector and/or element Jacobian matrix
  /// and/or the mass matrix
  ///
  /// flag=2: compute Jacobian, residuals and mass matrix
  /// flag=1: compute Jacobian and residuals
  /// flag=0: compute only residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  template<unsigned NREAGENT, unsigned DIM>
  void AdvectionDiffusionReactionEquations<NREAGENT, DIM>::
    fill_in_generic_residual_contribution_adv_diff_react(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Get the nodal index at which the unknown is stored
    unsigned c_nodal_index[NREAGENT];
    for (unsigned r = 0; r < NREAGENT; r++)
    {
      c_nodal_index[r] = c_index_adv_diff_react(r);
    }

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get diffusion coefficients
    Vector<double> D = diff();

    // Get the timescales
    Vector<double> T = tau();

    // Integers used to store the local equation number and local unknown
    // indices for the residuals and jacobians
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++) s[i] = integral_pt()->knot(ipt, i);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_adv_diff_react(
        ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the solution and its derivatives
      // Allocate
      Vector<double> interpolated_c(NREAGENT, 0.0);
      Vector<double> dcdt(NREAGENT, 0.0);
      Vector<double> interpolated_x(DIM, 0.0);
      DenseMatrix<double> interpolated_dcdx(NREAGENT, DIM, 0.0);
      Vector<double> mesh_velocity(DIM, 0.0);


      // Calculate function value and derivatives:
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions to calculate the position
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_x[j] += raw_nodal_position(l, j) * psi(l);
        }

        // Loop over the unknown reagents
        for (unsigned r = 0; r < NREAGENT; r++)
        {
          // Get the value at the node
          const double c_value = raw_nodal_value(l, c_nodal_index[r]);

          // Calculate the interpolated value
          interpolated_c[r] += c_value * psi(l);
          dcdt[r] += dc_dt_adv_diff_react(l, r) * psi(l);

          // Loop over directions to calculate the derivatie
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_dcdx(r, j) += c_value * dpsidx(l, j);
          }
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
      Vector<double> source(NREAGENT);
      get_source_adv_diff_react(ipt, interpolated_x, source);


      // Get wind
      Vector<double> wind(DIM);
      get_wind_adv_diff_react(ipt, s, interpolated_x, wind);

      // Get reaction terms
      Vector<double> R(NREAGENT);
      get_reaction_adv_diff_react(ipt, interpolated_c, R);

      // If we are getting the jacobian the get the derivative terms
      DenseMatrix<double> dRdC(NREAGENT);
      if (flag)
      {
        get_reaction_deriv_adv_diff_react(ipt, interpolated_c, dRdC);
      }


      // Assemble residuals and Jacobian
      //--------------------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the reagents
        for (unsigned r = 0; r < NREAGENT; r++)
        {
          // Set the local equation number
          local_eqn = nodal_local_eqn(l, c_nodal_index[r]);

          /*IF it's not a boundary condition*/
          if (local_eqn >= 0)
          {
            // Add body force/source/reaction term and time derivative
            residuals[local_eqn] -=
              (T[r] * dcdt[r] + source[r] + R[r]) * test(l) * W;

            // The Advection Diffusion bit itself
            for (unsigned k = 0; k < DIM; k++)
            {
              // Terms that multiply the test function
              double tmp = wind[k];
              // If the mesh is moving need to subtract the mesh velocity
              if (!ALE_is_disabled)
              {
                tmp -= T[r] * mesh_velocity[k];
              }
              // Now construct the contribution to the residuals
              residuals[local_eqn] -= interpolated_dcdx(r, k) *
                                      (tmp * test(l) + D[r] * dtestdx(l, k)) *
                                      W;
            }

            // Calculate the jacobian
            //-----------------------
            if (flag)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Loop over the reagents again
                for (unsigned r2 = 0; r2 < NREAGENT; r2++)
                {
                  // Set the number of the unknown
                  local_unknown = nodal_local_eqn(l2, c_nodal_index[r2]);

                  // If at a non-zero degree of freedom add in the entry
                  if (local_unknown >= 0)
                  {
                    // Diagonal terms (i.e. the basic equations)
                    if (r2 == r)
                    {
                      // Mass matrix term
                      jacobian(local_eqn, local_unknown) -=
                        T[r] * test(l) * psi(l2) *
                        node_pt(l2)->time_stepper_pt()->weight(1, 0) * W;

                      // Add the mass matrix term
                      if (flag == 2)
                      {
                        mass_matrix(local_eqn, local_unknown) +=
                          T[r] * test(l) * psi(l2) * W;
                      }

                      // Add contribution to Elemental Matrix
                      for (unsigned i = 0; i < DIM; i++)
                      {
                        // Temporary term used in assembly
                        double tmp = wind[i];
                        if (!ALE_is_disabled) tmp -= T[r] * mesh_velocity[i];
                        // Now assemble Jacobian term
                        jacobian(local_eqn, local_unknown) -=
                          dpsidx(l2, i) *
                          (tmp * test(l) + D[r] * dtestdx(l, i)) * W;
                      }

                    } // End of diagonal terms

                    // Now add the cross-reaction terms
                    jacobian(local_eqn, local_unknown) -=
                      dRdC(r, r2) * psi(l2) * test(l) * W;
                  }
                }
              }
            } // End of jacobian
          }
        } // End of loop over reagents
      } // End of loop over nodes
    } // End of loop over integration points
  }


  //=======================================================================
  /// Compute norm of the solution: sum of squares of L2 norms for reagents
  //=======================================================================
  template<unsigned NREAGENT, unsigned DIM>
  void AdvectionDiffusionReactionEquations<NREAGENT, DIM>::compute_norm(
    double& norm)
  {
    // Initialise
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Solution
    double c = 0.0;

    // Find out how many nodes there are in the element
    unsigned n_node = this->nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = this->J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Loop over all reagents
      for (unsigned r = 0; r < NREAGENT; ++r)
      {
        // Get FE function value
        c = this->interpolated_c_adv_diff_react(s, r);
        // Add to  norm
        norm += c * c * W;
      }
    }
  }


  //======================================================================
  /// Self-test:  Return 0 for OK
  //======================================================================
  template<unsigned NREAGENT, unsigned DIM>
  unsigned AdvectionDiffusionReactionEquations<NREAGENT, DIM>::self_test()
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


  //=========================================================================
  /// Integrate the reagent concentrations over the element
  //========================================================================
  template<unsigned NREAGENT, unsigned DIM>
  void AdvectionDiffusionReactionEquations<NREAGENT, DIM>::integrate_reagents(
    Vector<double>& C) const
  {
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Get the nodal index at which the unknown is stored
    unsigned c_nodal_index[NREAGENT];
    for (unsigned r = 0; r < NREAGENT; r++)
    {
      c_nodal_index[r] = c_index_adv_diff_react(r);
    }

    // Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_eulerian_at_knot(ipt, psi, dpsidx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the solution and its derivatives
      // Allocate
      Vector<double> interpolated_c(NREAGENT, 0.0);

      // Calculate function value and derivatives:
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the unknown reagents
        for (unsigned r = 0; r < NREAGENT; r++)
        {
          // Get the value at the node
          const double c_value = raw_nodal_value(l, c_nodal_index[r]);

          // Calculate the interpolated value
          interpolated_c[r] += c_value * psi(l);
        }
      }

      for (unsigned r = 0; r < NREAGENT; r++)
      {
        C[r] += interpolated_c[r] * W;
      }
    } // End of loop over integration points
  }


  //======================================================================
  ///  Output function:
  ///
  ///   x,y,u,w_x,w_y   or    x,y,z,u,w_x,w_y,w_z
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  template<unsigned NREAGENT, unsigned DIM>
  void AdvectionDiffusionReactionEquations<NREAGENT, DIM>::output(
    std::ostream& outfile, const unsigned& nplot)
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

      // Get Eulerian coordinate of plot point
      Vector<double> x(DIM);
      interpolated_x(s, x);

      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      for (unsigned i = 0; i < NREAGENT; i++)
      {
        outfile << interpolated_c_adv_diff_react(s, i) << " ";
      }

      // Get the wind
      Vector<double> wind(DIM);
      // Dummy integration point needed
      unsigned ipt = 0;
      get_wind_adv_diff_react(ipt, s, x, wind);
      for (unsigned i = 0; i < DIM; i++)
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
  ///   x,y,u   or    x,y,z,u
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  template<unsigned NREAGENT, unsigned DIM>
  void AdvectionDiffusionReactionEquations<NREAGENT, DIM>::output(
    FILE* file_pt, const unsigned& nplot)
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
      for (unsigned i = 0; i < NREAGENT; i++)
      {
        fprintf(file_pt, "%g \n", interpolated_c_adv_diff_react(s, i));
      }
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, nplot);
  }


  //======================================================================
  ///   Output exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  ///
  ///   x,y,u_exact    or    x,y,z,u_exact
  //======================================================================
  template<unsigned NREAGENT, unsigned DIM>
  void AdvectionDiffusionReactionEquations<NREAGENT, DIM>::output_fct(
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
  ///  Validate against exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  ///
  //======================================================================
  template<unsigned NREAGENT, unsigned DIM>
  void AdvectionDiffusionReactionEquations<NREAGENT, DIM>::compute_error(
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
      double u_fe = interpolated_c_adv_diff_react(s, 0);

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

  //====================================================================
  // Force build of templates, only building binary reactions at present
  //====================================================================
  /// One reagent only (!)
  template class AdvectionDiffusionReactionEquations<1, 1>;
  template class AdvectionDiffusionReactionEquations<1, 2>;
  template class AdvectionDiffusionReactionEquations<1, 3>;

  template class QAdvectionDiffusionReactionElement<1, 1, 2>;
  template class QAdvectionDiffusionReactionElement<1, 1, 3>;
  template class QAdvectionDiffusionReactionElement<1, 1, 4>;

  template class QAdvectionDiffusionReactionElement<1, 2, 2>;
  template class QAdvectionDiffusionReactionElement<1, 2, 3>;
  template class QAdvectionDiffusionReactionElement<1, 2, 4>;

  template class QAdvectionDiffusionReactionElement<1, 3, 2>;
  template class QAdvectionDiffusionReactionElement<1, 3, 3>;
  template class QAdvectionDiffusionReactionElement<1, 3, 4>;

  // Two reagents
  template class AdvectionDiffusionReactionEquations<2, 1>;
  template class AdvectionDiffusionReactionEquations<2, 2>;
  template class AdvectionDiffusionReactionEquations<2, 3>;

  template class QAdvectionDiffusionReactionElement<2, 1, 2>;
  template class QAdvectionDiffusionReactionElement<2, 1, 3>;
  template class QAdvectionDiffusionReactionElement<2, 1, 4>;

  template class QAdvectionDiffusionReactionElement<2, 2, 2>;
  template class QAdvectionDiffusionReactionElement<2, 2, 3>;
  template class QAdvectionDiffusionReactionElement<2, 2, 4>;

  template class QAdvectionDiffusionReactionElement<2, 3, 2>;
  template class QAdvectionDiffusionReactionElement<2, 3, 3>;
  template class QAdvectionDiffusionReactionElement<2, 3, 4>;

} // namespace oomph
