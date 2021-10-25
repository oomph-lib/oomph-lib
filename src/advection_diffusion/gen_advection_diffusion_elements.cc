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
// Non-inline functions for GeneralisedAdvection Diffusion elements
#include "gen_advection_diffusion_elements.h"

namespace oomph
{
  /// 2D GeneralisedAdvection Diffusion elements


  /// Default value for Peclet number
  template<unsigned DIM>
  double GeneralisedAdvectionDiffusionEquations<DIM>::Default_peclet_number =
    0.0;

  //======================================================================
  /// Compute element residual Vector and/or element Jacobian matrix
  ///
  /// flag=1: compute both
  /// flag=0: compute only residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  template<unsigned DIM>
  void GeneralisedAdvectionDiffusionEquations<DIM>::
    fill_in_generic_residual_contribution_cons_adv_diff(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Get the nodal index at which the unknown is stored
    const unsigned u_nodal_index = u_index_cons_adv_diff();

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get Peclet number
    const double peclet = pe();

    // Get the Peclet*Strouhal number
    const double peclet_st = pe_st();

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
      double J = dshape_and_dtest_eulerian_at_knot_cons_adv_diff(
        ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the solution and its derivatives
      // Allocate
      double interpolated_u = 0.0;
      double dudt = 0.0;
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> interpolated_dudx(DIM, 0.0);
      Vector<double> mesh_velocity(DIM, 0.0);


      // Calculate function value and derivatives:
      //-----------------------------------------
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the value at the node
        double u_value = raw_nodal_value(l, u_nodal_index);
        interpolated_u += u_value * psi(l);
        dudt += du_dt_cons_adv_diff(l) * psi(l);
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
      double source;
      get_source_cons_adv_diff(ipt, interpolated_x, source);


      // Get wind
      //--------
      Vector<double> wind(DIM);
      get_wind_cons_adv_diff(ipt, s, interpolated_x, wind);

      // Get the conserved wind (non-divergence free)
      Vector<double> conserved_wind(DIM);
      get_conserved_wind_cons_adv_diff(ipt, s, interpolated_x, conserved_wind);


      // Get diffusivity tensor
      DenseMatrix<double> D(DIM, DIM);
      get_diff_cons_adv_diff(ipt, s, interpolated_x, D);

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
          // Add body force/source term and time derivative
          residuals[local_eqn] -= (peclet_st * dudt + source) * test(l) * W;

          // The Generalised Advection Diffusion bit itself
          for (unsigned k = 0; k < DIM; k++)
          {
            // Terms that multiply the test function
            // divergence-free wind
            double tmp = peclet * wind[k];
            // If the mesh is moving need to subtract the mesh velocity
            if (!ALE_is_disabled)
            {
              tmp -= peclet_st * mesh_velocity[k];
            }
            tmp *= interpolated_dudx[k];

            // Terms that multiply the derivative of the test function
            // Advective term
            double tmp2 = -conserved_wind[k] * interpolated_u;
            // Now the diuffusive term
            for (unsigned j = 0; j < DIM; j++)
            {
              tmp2 += D(k, j) * interpolated_dudx[j];
            }
            // Now construct the contribution to the residuals
            residuals[local_eqn] -= (tmp * test(l) + tmp2 * dtestdx(l, k)) * W;
          }

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
                  peclet_st * test(l) * psi(l2) *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * W;

                // Add the mass matrix term
                if (flag == 2)
                {
                  mass_matrix(local_eqn, local_unknown) +=
                    peclet_st * test(l) * psi(l2) * W;
                }

                // Add contribution to Elemental Matrix
                for (unsigned k = 0; k < DIM; k++)
                {
                  // Temporary term used in assembly
                  double tmp = -peclet * wind[k];
                  if (!ALE_is_disabled)
                  {
                    tmp -= peclet_st * mesh_velocity[k];
                  }
                  tmp *= dpsidx(l2, k);

                  double tmp2 = -conserved_wind[k] * psi(l2);
                  // Now the diffusive term
                  for (unsigned j = 0; j < DIM; j++)
                  {
                    tmp2 += D(k, j) * dpsidx(l2, j);
                  }

                  // Now assemble Jacobian term
                  jacobian(local_eqn, local_unknown) -=
                    (tmp * test(l) + tmp2 * dtestdx(l, k)) * W;
                }
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
  template<unsigned DIM>
  unsigned GeneralisedAdvectionDiffusionEquations<DIM>::self_test()
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
  ///   x,y,u,w_x,w_y   or    x,y,z,u,w_x,w_y,w_z
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  template<unsigned DIM>
  void GeneralisedAdvectionDiffusionEquations<DIM>::output(
    std::ostream& outfile, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);


    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    const unsigned n_node = this->nnode();
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

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
      outfile << interpolated_u_cons_adv_diff(s) << " ";

      // Get the gradients
      (void)this->dshape_eulerian(s, psi, dpsidx);
      Vector<double> interpolated_dudx(DIM, 0.0);
      for (unsigned n = 0; n < n_node; n++)
      {
        const double u_ = this->nodal_value(n, 0);
        for (unsigned i = 0; i < DIM; i++)
        {
          interpolated_dudx[i] += u_ * dpsidx(n, i);
        }
      }

      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_dudx[i] << " ";
      }

      // Get the wind
      Vector<double> wind(DIM);
      // Dummy integration point variable needed
      unsigned ipt = 0;
      get_wind_cons_adv_diff(ipt, s, x, wind);
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
  template<unsigned DIM>
  void GeneralisedAdvectionDiffusionEquations<DIM>::output(
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
      fprintf(file_pt, "%g \n", interpolated_u_cons_adv_diff(s));
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
  ///   x,y,u_exact    or    x,y,z,u_exact
  //======================================================================
  template<unsigned DIM>
  void GeneralisedAdvectionDiffusionEquations<DIM>::output_fct(
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
  /// Validate against exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  ///
  //======================================================================
  template<unsigned DIM>
  void GeneralisedAdvectionDiffusionEquations<DIM>::compute_error(
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
      double u_fe = interpolated_u_cons_adv_diff(s);

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
  /// Calculate the integrated value of the unknown over the element
  ///
  //======================================================================
  template<unsigned DIM>
  double GeneralisedAdvectionDiffusionEquations<DIM>::integrate_u()
  {
    // Initialise
    double sum = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Find out how many nodes there are in the element
    const unsigned n_node = nnode();

    // Find the index at which the concentration is stored
    const unsigned u_nodal_index = this->u_index_cons_adv_diff();

    // Allocate memory for the shape functions
    Shape psi(n_node);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      const double w = integral_pt()->weight(ipt);

      // Get the shape functions
      this->shape_at_knot(ipt, psi);

      // Calculate the concentration
      double interpolated_u = 0.0;
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += this->nodal_value(l, u_nodal_index) * psi(l);
      }

      // Get jacobian of mapping
      const double J = J_eulerian_at_knot(ipt);

      // Add the values to the sum
      sum += interpolated_u * w * J;
    }

    // return the sum
    return sum;
  }


  //======================================================================
  // Set the data for the number of Variables at each node
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned
    QGeneralisedAdvectionDiffusionElement<DIM, NNODE_1D>::Initial_Nvalue = 1;

  //====================================================================
  // Force build of templates
  //====================================================================
  template class GeneralisedAdvectionDiffusionEquations<1>;
  template class GeneralisedAdvectionDiffusionEquations<2>;
  template class GeneralisedAdvectionDiffusionEquations<3>;

  template class QGeneralisedAdvectionDiffusionElement<1, 2>;
  template class QGeneralisedAdvectionDiffusionElement<1, 3>;
  template class QGeneralisedAdvectionDiffusionElement<1, 4>;

  template class QGeneralisedAdvectionDiffusionElement<2, 2>;
  template class QGeneralisedAdvectionDiffusionElement<2, 3>;
  template class QGeneralisedAdvectionDiffusionElement<2, 4>;

  template class QGeneralisedAdvectionDiffusionElement<3, 2>;
  template class QGeneralisedAdvectionDiffusionElement<3, 3>;
  template class QGeneralisedAdvectionDiffusionElement<3, 4>;


} // namespace oomph
