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
#include "axisym_poroelasticity_elements.h"

namespace oomph
{
  //===================================================================
  /// Static default value for Young's modulus (1.0 -- for natural
  /// scaling, i.e. all stresses have been non-dimensionalised by
  /// the same reference Young's modulus. Setting the "non-dimensional"
  /// Young's modulus (obtained by de-referencing Youngs_modulus_pt)
  /// to a number larger than one means that the material is stiffer
  /// than assumed in the non-dimensionalisation.
  //===================================================================
  double AxisymmetricPoroelasticityEquations::Default_youngs_modulus_value =
    1.0;

  //========================================================================
  /// Static default value for timescale ratio (1.0)
  //===================================================================
  double AxisymmetricPoroelasticityEquations::Default_lambda_sq_value = 1.0;

  //===================================================================
  /// Static default value for the density ratio (fluid to solid)
  //===================================================================
  double AxisymmetricPoroelasticityEquations::Default_density_ratio_value = 1.0;

  //===================================================================
  /// Static default value for permeability (1.0 for natural scaling
  /// i.e. timescale is given by L^2/(k^* E)
  //===================================================================
  double AxisymmetricPoroelasticityEquations::Default_permeability_value = 1.0;


  //===================================================================
  /// Static default value for ratio of the material's actual
  /// permeability to the permeability used in the non-dimensionalisastion
  /// of the equations
  //===================================================================
  double AxisymmetricPoroelasticityEquations::Default_permeability_ratio_value =
    1.0;

  //===================================================================
  /// Static default value for alpha, the Biot parameter
  //===================================================================
  double AxisymmetricPoroelasticityEquations::Default_alpha_value = 1.0;

  //===================================================================
  /// Static default value for the porosity
  //===================================================================
  double AxisymmetricPoroelasticityEquations::Default_porosity_value = 1.0;

  //======================================================================
  /// Performs a div-conserving transformation of the vector basis
  /// functions from the reference element to the actual element
  //======================================================================
  double AxisymmetricPoroelasticityEquations::transform_basis(
    const Vector<double>& s,
    const Shape& q_basis_local,
    Shape& psi,
    DShape& dpsi,
    Shape& q_basis) const
  {
    // Call the (geometric) shape functions and their derivatives
    this->dshape_local(s, psi, dpsi);

    // Storage for the (geometric) jacobian and its inverse
    DenseMatrix<double> jacobian(2), inverse_jacobian(2);

    // Get the jacobian of the geometric mapping and its determinant
    double det = local_to_eulerian_mapping(dpsi, jacobian, inverse_jacobian);

    // Transform the derivative of the geometric basis so that it's w.r.t.
    // global coordinates
    this->transform_derivatives(inverse_jacobian, dpsi);

    // Get the number of computational basis vectors
    const unsigned n_q_basis = this->nq_basis();

    // Loop over the basis vectors
    for (unsigned l = 0; l < n_q_basis; l++)
    {
      // Loop over the spatial components
      for (unsigned i = 0; i < 2; i++)
      {
        // Zero the basis
        q_basis(l, i) = 0.0;
      }
    }

    // Loop over the spatial components
    for (unsigned i = 0; i < 2; i++)
    {
      // And again
      for (unsigned j = 0; j < 2; j++)
      {
        // Get the element of the jacobian (must transpose it due to different
        // conventions) and divide by the determinant
        double jac_trans = jacobian(j, i) / det;

        // Loop over the computational basis vectors
        for (unsigned l = 0; l < n_q_basis; l++)
        {
          // Apply the appropriate div-conserving mapping
          q_basis(l, i) += jac_trans * q_basis_local(l, j);
        }
      }
    }

    // Scale the basis by the ratio of the length of the edge to the length of
    // the corresponding edge on the reference element
    scale_basis(q_basis);

    return det;
  }


  //========================================================================
  /// Output FE representation of soln:
  /// x,y,u1,u2,q1,q2,div_q,p at Nplot^2 plot points
  //========================================================================
  void AxisymmetricPoroelasticityEquations::output(std::ostream& outfile,
                                                   const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Skeleton velocity
    Vector<double> du_dt(2);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);

    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Output the components of the position
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      // Output the components of the FE representation of skeleton displ. at s
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_u(s, i) << " "; // soln 0 and 1
      }

      // Output the components of the FE representation of q at s
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_q(s, i) << " "; // soln 2 and 3
      }

      // Output FE representation of div u at s
      outfile << interpolated_div_q(s) << " "; // soln 4

      // Output FE representation of p at s
      outfile << interpolated_p(s) << " "; // soln 5

      // Skeleton velocity
      interpolated_du_dt(s, du_dt);
      outfile << du_dt[0] << " "; // soln 6
      outfile << du_dt[1] << " "; // soln 7

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  //============================================================================
  /// Output FE representation of exact soln at
  /// Nplot^2 plot points
  //============================================================================
  void AxisymmetricPoroelasticityEquations::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordintes
    Vector<double> x(2);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Exact solution Vector
    Vector<double> exact_soln(13);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);

    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }
      for (unsigned i = 0; i < 13; i++)
      {
        outfile << exact_soln[i] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  //========================================================================
  /// Output FE representation of exact soln at
  /// Nplot^2 plot points. Unsteady version
  //========================================================================
  void AxisymmetricPoroelasticityEquations::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordintes
    Vector<double> x(2);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Exact solution Vector
    Vector<double> exact_soln(13);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);

    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(time, x, exact_soln);

      // Output
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }
      for (unsigned i = 0; i < 13; i++)
      {
        outfile << exact_soln[i] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  //========================================================================
  /// Compute the error between the FE solution and the exact solution
  /// using the H(div) norm for q and L^2 norm for u and p
  //========================================================================
  void AxisymmetricPoroelasticityEquations::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    Vector<double>& error,
    Vector<double>& norm)
  {
    for (unsigned i = 0; i < 3; i++)
    {
      error[i] = 0.0;
      norm[i] = 0.0;
    }

    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordinates
    Vector<double> x(2);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    outfile << "ZONE" << std::endl;

    // Exact solution Vector
    Vector<double> exact_soln(13);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = this->J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Displacement error
      for (unsigned i = 0; i < 2; i++)
      {
        norm[0] += exact_soln[i] * exact_soln[i] * W;
        // Error due to q_i
        error[0] += (exact_soln[i] - this->interpolated_u(s, i)) *
                    (exact_soln[i] - this->interpolated_u(s, i)) * W;
      }

      // Flux error
      for (unsigned i = 0; i < 2; i++)
      {
        norm[1] += exact_soln[2 + i] * exact_soln[2 + i] * W;
        // Error due to q_i
        error[1] += (exact_soln[2 + i] - this->interpolated_q(s, i)) *
                    (exact_soln[2 + i] - this->interpolated_q(s, i)) * W;
      }

      // Flux divergence error
      norm[1] += exact_soln[2 * 2] * exact_soln[2 * 2] * W;
      error[1] += (exact_soln[2 * 2] - interpolated_div_q(s)) *
                  (exact_soln[2 * 2] - interpolated_div_q(s)) * W;

      // Pressure error
      norm[2] += exact_soln[2 * 2 + 1] * exact_soln[2 * 2 + 1] * W;
      error[2] += (exact_soln[2 * 2 + 1] - this->interpolated_p(s)) *
                  (exact_soln[2 * 2 + 1] - this->interpolated_p(s)) * W;

      // Output x,y,[z]
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output u_1_error,u_2_error,...
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << exact_soln[i] - this->interpolated_u(s, i) << " ";
      }

      // Output q_1_error,q_2_error,...
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << exact_soln[2 + i] - this->interpolated_q(s, i) << " ";
      }

      // Output p_error
      outfile << exact_soln[2 * 2 + 1] - this->interpolated_p(s) << " ";

      outfile << std::endl;
    }
  }

  //========================================================================
  /// Compute the error between the FE solution and the exact solution
  /// using the H(div) norm for u and L^2 norm for p. Unsteady version
  //========================================================================
  void AxisymmetricPoroelasticityEquations::compute_error(
    std::ostream& outfile,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
    Vector<double>& error,
    Vector<double>& norm)
  {
    for (unsigned i = 0; i < 3; i++)
    {
      error[i] = 0.0;
      norm[i] = 0.0;
    }

    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordinates
    Vector<double> x(2);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    outfile << "ZONE" << std::endl;

    // Exact solution Vector
    Vector<double> exact_soln(13);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = this->J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(time, x, exact_soln);

      // Displacement error
      for (unsigned i = 0; i < 2; i++)
      {
        norm[0] += exact_soln[i] * exact_soln[i] * W;
        // Error due to q_i
        error[0] += (exact_soln[i] - this->interpolated_u(s, i)) *
                    (exact_soln[i] - this->interpolated_u(s, i)) * W;
      }

      // Flux error
      for (unsigned i = 0; i < 2; i++)
      {
        norm[1] += exact_soln[2 + i] * exact_soln[2 + i] * W;
        // Error due to q_i
        error[1] += (exact_soln[2 + i] - this->interpolated_q(s, i)) *
                    (exact_soln[2 + i] - this->interpolated_q(s, i)) * W;
      }

      // Flux divergence error
      norm[1] += exact_soln[2 * 2] * exact_soln[2 * 2] * W;
      error[1] += (exact_soln[2 * 2] - interpolated_div_q(s)) *
                  (exact_soln[2 * 2] - interpolated_div_q(s)) * W;

      // Pressure error
      norm[2] += exact_soln[2 * 2 + 1] * exact_soln[2 * 2 + 1] * W;
      error[2] += (exact_soln[2 * 2 + 1] - this->interpolated_p(s)) *
                  (exact_soln[2 * 2 + 1] - this->interpolated_p(s)) * W;

      // Output x,y,[z]
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output u_1_error,u_2_error,...
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << exact_soln[i] - this->interpolated_u(s, i) << " ";
      }

      // Output q_1_error,q_2_error,...
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << exact_soln[2 + i] - this->interpolated_q(s, i) << " ";
      }

      // Output p_error
      outfile << exact_soln[2 * 2 + 1] - this->interpolated_p(s) << " ";

      outfile << std::endl;
    }
  }

  //========================================================================
  /// Fill in residuals and, if flag==true, jacobian
  //========================================================================
  void AxisymmetricPoroelasticityEquations::
    fill_in_generic_residual_contribution(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian,
                                          bool flag)
  {
    // Get the number of geometric nodes, total number of basis functions,
    // and number of edges basis functions
    const unsigned n_node = nnode();
    const unsigned n_q_basis = nq_basis();
    const unsigned n_q_basis_edge = nq_basis_edge();
    const unsigned n_p_basis = np_basis();

    // Storage for the geometric and computational bases and the test functions
    Shape psi(n_node), u_basis(n_node), u_test(n_node), q_basis(n_q_basis, 2),
      q_test(n_q_basis, 2), p_basis(n_p_basis), p_test(n_p_basis),
      div_q_basis_ds(n_q_basis), div_q_test_ds(n_q_basis);

    DShape dpsidx(n_node, 2), du_basis_dx(n_node, 2), du_test_dx(n_node, 2);

    // Get the number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Storage for the local coordinates
    Vector<double> s(2);

    // Storage for the elasticity source function
    Vector<double> f_solid(2);

    // Storage for the source function
    Vector<double> f_fluid(2);

    // Storage for the mass source function
    double mass_source_local = 0.0;

    // Get elastic parameters
    double nu_local = this->nu();
    double youngs_modulus_local = this->youngs_modulus();

    // Obtain Lame parameters from Young's modulus and Poisson's ratio
    double lambda = youngs_modulus_local * nu_local / (1.0 + nu_local) /
                    (1.0 - 2.0 * nu_local);

    double mu = youngs_modulus_local / 2.0 / (1.0 + nu_local);

    // Storage for Lambda_sq
    double lambda_sq = this->lambda_sq();

    // Get the value of permeability
    double local_permeability = this->permeability();

    // Ratio of the material's permeability to the permeability used
    // to non-dimensionalise the equations
    double local_permeability_ratio = this->permeability_ratio();

    // Get the value of alpha
    double alpha = this->alpha();

    // Get the value of the porosity
    double porosity = this->porosity();

    // Get the density ratio
    double density_ratio = this->density_ratio();

    // Precompute the ratio of fluid density to combined density
    double rho_f_over_rho =
      density_ratio / (1.0 + porosity * (density_ratio - 1.0));

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Find the local coordinates at the integration point
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the weight of the intetgration point
      double w = integral_pt()->weight(ipt);

      // Call the basis functions and test functions and get the
      // (geometric) jacobian of the current element
      double J = shape_basis_test_local_at_knot(ipt,
                                                psi,
                                                dpsidx,
                                                u_basis,
                                                u_test,
                                                du_basis_dx,
                                                du_test_dx,
                                                q_basis,
                                                q_test,
                                                p_basis,
                                                p_test,
                                                div_q_basis_ds,
                                                div_q_test_ds);

      // Storage for interpolated values
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> interpolated_u(2, 0.0);
      DenseMatrix<double> interpolated_du_dx(2, 2, 0.0);
      double interpolated_div_du_dt_dx = 0.0;
      double interpolated_du_r_dt = 0.0;
      Vector<double> interpolated_d2u_dt2(2, 0.0);
      Vector<double> interpolated_q(2, 0.0);
      double interpolated_div_q_ds = 0.0;
      Vector<double> interpolated_dq_dt(2, 0.0);
      double interpolated_p = 0.0;

      // loop over geometric basis functions to find interpolated x
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the geometric basis functions
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += nodal_position(l, i) * psi(l);
          interpolated_d2u_dt2[i] += this->d2u_dt2(l, i) * u_basis(l);

          // Get the nodal displacements
          const double u_value =
            this->raw_nodal_value(l, u_index_axisym_poroelasticity(i));
          interpolated_u[i] += u_value * u_basis(l);

          // Loop over derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_du_dx(i, j) += u_value * du_basis_dx(l, j);
          }

          // divergence of the time derivative of the solid displacement
          interpolated_div_du_dt_dx += this->du_dt(l, i) * du_basis_dx(l, i);
        }

        // r-component of the solid velocity
        interpolated_du_r_dt += du_dt(l, 0) * u_basis(l);
      }

      // loop over the nodes and use the vector basis functions to find the
      // interpolated flux and its time derivative
      for (unsigned i = 0; i < 2; i++)
      {
        // Loop over the edge basis vectors
        for (unsigned l = 0; l < n_q_basis_edge; l++)
        {
          interpolated_q[i] += q_edge(l) * q_basis(l, i);
          interpolated_dq_dt[i] += dq_edge_dt(l) * q_basis(l, i);
        }
        // Loop over the internal basis vectors
        for (unsigned l = n_q_basis_edge; l < n_q_basis; l++)
        {
          interpolated_q[i] += q_internal(l - n_q_basis_edge) * q_basis(l, i);
          interpolated_dq_dt[i] +=
            dq_internal_dt(l - n_q_basis_edge) * q_basis(l, i);
        }
      }

      // loop over the pressure basis and find the interpolated pressure
      for (unsigned l = 0; l < n_p_basis; l++)
      {
        interpolated_p += p_value(l) * p_basis(l);
      }

      // loop over the q edge divergence basis and the q internal divergence
      // basis to find interpolated div q
      for (unsigned l = 0; l < n_q_basis_edge; l++)
      {
        interpolated_div_q_ds += q_edge(l) * div_q_basis_ds(l);
      }
      for (unsigned l = n_q_basis_edge; l < n_q_basis; l++)
      {
        interpolated_div_q_ds +=
          q_internal(l - n_q_basis_edge) * div_q_basis_ds(l);
      }

      // Get the solid body force
      this->solid_body_force(time, interpolated_x, f_solid);

      // Get the fluid nody force
      this->fluid_body_force(time, interpolated_x, f_fluid);

      // Get the mass source function
      this->mass_source(time, interpolated_x, mass_source_local);

      double r = interpolated_x[0];

      // Linear elasticity:
      //-------------------

      double u_r = interpolated_u[0];
      double du_r_dr = interpolated_du_dx(0, 0);
      double du_r_dz = interpolated_du_dx(0, 1);
      double du_z_dr = interpolated_du_dx(1, 0);
      double du_z_dz = interpolated_du_dx(1, 1);

      // Storage for terms of Jacobian
      double G_r = 0, G_z = 0;

      for (unsigned l = 0; l < n_node; l++)
      {
        for (unsigned a = 0; a < 2; a++)
        {
          local_eqn =
            this->nodal_local_eqn(l, u_index_axisym_poroelasticity(a));

          if (local_eqn >= 0)
          {
            residuals[local_eqn] +=
              (lambda_sq *
                 (interpolated_d2u_dt2[a] +
                  rho_f_over_rho * local_permeability * interpolated_dq_dt[a]) -
               f_solid[a]) *
              u_test(l) * r * w * J;

            // r-equation
            if (a == 0)
            {
              residuals[local_eqn] +=
                (mu * (2.0 * du_r_dr * du_test_dx(l, 0) +
                       du_test_dx(l, 1) * (du_r_dz + du_z_dr) +
                       2.0 * u_test(l) / pow(r, 2) * (u_r)) +
                 (lambda * (du_r_dr + u_r / r + du_z_dz) -
                  alpha * interpolated_p) *
                   (du_test_dx(l, 0) + u_test(l) / r)) *
                r * w * J;
            }
            else if (a == 1)
            {
              residuals[local_eqn] +=
                (mu * (du_test_dx(l, 0) * (du_r_dz + du_z_dr) +
                       2.0 * du_z_dz * du_test_dx(l, 1)) +
                 (lambda * (du_r_dr + u_r / r + du_z_dz) -
                  alpha * interpolated_p) *
                   du_test_dx(l, 1)) *
                r * w * J;
            }
            // error: a should be 0 or 1
            else
            {
              throw OomphLibError("a should equal 0 or 1",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }

            // Jacobian entries
            if (flag)
            {
              // d(u_eqn_l,a)/d(U_l2,c)
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                if (a == 0)
                {
                  G_r =
                    (mu * (2.0 * du_basis_dx(l2, 0) * du_test_dx(l, 0) +
                           du_test_dx(l, 1) * du_basis_dx(l2, 1) +
                           2 * u_test(l) / pow(r, 2) * u_basis(l2)) +
                     (lambda * (du_basis_dx(l2, 0) + u_basis(l2) / r)) *
                       (du_test_dx(l, 0) + u_test(l) / r) +
                     lambda_sq * node_pt(l2)->time_stepper_pt()->weight(2, 0) *
                       u_basis(l2) * u_test(l)) *
                    r * w * J;

                  G_z = (mu * du_test_dx(l, 1) * du_basis_dx(l2, 0) +
                         lambda * du_basis_dx(l2, 1) *
                           (du_test_dx(l, 0) + u_test(l) / r)) *
                        r * w * J;
                }
                else if (a == 1)
                {
                  G_r = (mu * du_test_dx(l, 0) * du_basis_dx(l2, 1) +
                         lambda * (du_basis_dx(l2, 0) + u_basis(l2) / r) *
                           du_test_dx(l, 1)) *
                        r * w * J;
                  G_z =
                    (mu * (du_test_dx(l, 0) * du_basis_dx(l2, 0) +
                           2.0 * du_basis_dx(l2, 1) * du_test_dx(l, 1)) +
                     lambda * du_basis_dx(l2, 1) * du_test_dx(l, 1) +
                     lambda_sq * node_pt(l2)->time_stepper_pt()->weight(2, 0) *
                       u_basis(l2) * u_test(l)) *
                    r * w * J;
                }

                for (unsigned c = 0; c < 2; c++)
                {
                  local_unknown =
                    this->nodal_local_eqn(l2, u_index_axisym_poroelasticity(c));
                  if (local_unknown >= 0)
                  {
                    if (c == 0)
                    {
                      jacobian(local_eqn, local_unknown) += G_r;
                    }
                    else if (c == 1)
                    {
                      jacobian(local_eqn, local_unknown) += G_z;
                    }
                  }
                }
              }

              // d(u_eqn_l,a)/d(Q_l2)
              for (unsigned l2 = 0; l2 < n_q_basis; l2++)
              {
                TimeStepper* timestepper_pt = 0;

                if (l2 < n_q_basis_edge)
                {
                  local_unknown = q_edge_local_eqn(l2);
                  timestepper_pt =
                    this->node_pt(q_edge_node_number(l2))->time_stepper_pt();
                }
                else // n_q_basis_edge <= l < n_basis
                {
                  local_unknown = q_internal_local_eqn(l2 - n_q_basis_edge);
                  timestepper_pt = this->internal_data_pt(q_internal_index())
                                     ->time_stepper_pt();
                }

                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    lambda_sq * rho_f_over_rho * local_permeability *
                    timestepper_pt->weight(1, 0) * q_basis(l2, a) * u_test(l) *
                    r * w * J;
                }
              }

              // d(u_eqn_l,a)/d(P_l2)
              for (unsigned l2 = 0; l2 < n_p_basis; l2++)
              {
                local_unknown = p_local_eqn(l2);
                if (local_unknown >= 0)
                {
                  if (a == 0)
                  {
                    jacobian(local_eqn, local_unknown) -=
                      alpha * p_basis(l2) * (du_test_dx(l, 0) + u_test(l) / r) *
                      r * w * J;
                  }
                  else if (a == 1)
                  {
                    jacobian(local_eqn, local_unknown) -=
                      alpha * p_basis(l2) * du_test_dx(l, 1) * r * w * J;
                  }
                }
              }
            } // End of Jacobian entries
          } // End of if not boundary condition
        } // End of loop over dimensions
      } // End of loop over u test functions


      // Darcy:
      //-------

      // Loop over the test functions
      for (unsigned l = 0; l < n_q_basis; l++)
      {
        if (l < n_q_basis_edge)
        {
          local_eqn = q_edge_local_eqn(l);
        }
        else // n_q_basis_edge <= l < n_basis
        {
          local_eqn = q_internal_local_eqn(l - n_q_basis_edge);
        }

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          for (unsigned i = 0; i < 2; i++)
          {
            residuals[local_eqn] +=
              (rho_f_over_rho * lambda_sq *
                 (interpolated_d2u_dt2[i] +
                  (local_permeability / porosity) * interpolated_dq_dt[i]) +
               interpolated_q[i] / local_permeability_ratio -
               rho_f_over_rho * f_fluid[i]) *
              q_test(l, i) * r * w * J;
          }

          // deliberately no jacobian factor in this integral
          residuals[local_eqn] -= interpolated_p * div_q_test_ds(l) * r * w;

          // deliberately no r factor in this integral
          residuals[local_eqn] -= interpolated_p * q_test(l, 0) * w * J;

          // Jacobian entries
          if (flag)
          {
            // d(q_eqn_l)/d(U_l2,c)
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              for (unsigned c = 0; c < 2; c++)
              {
                local_unknown =
                  this->nodal_local_eqn(l2, u_index_axisym_poroelasticity(c));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    rho_f_over_rho * lambda_sq *
                    this->node_pt(l2)->time_stepper_pt()->weight(2, 0) *
                    u_basis(l2) * q_test(l, c) * r * w * J;
                }
              }
            }

            // d(q_eqn_l)/d(Q_l2)
            for (unsigned l2 = 0; l2 < n_q_basis; l2++)
            {
              TimeStepper* timestepper_pt = 0;

              if (l2 < n_q_basis_edge)
              {
                local_unknown = q_edge_local_eqn(l2);
                timestepper_pt =
                  this->node_pt(q_edge_node_number(l2))->time_stepper_pt();
              }
              else // n_q_basis_edge <= l < n_basis
              {
                local_unknown = q_internal_local_eqn(l2 - n_q_basis_edge);
                timestepper_pt =
                  this->internal_data_pt(q_internal_index())->time_stepper_pt();
              }

              if (local_unknown >= 0)
              {
                for (unsigned c = 0; c < 2; c++)
                {
                  jacobian(local_eqn, local_unknown) +=
                    q_basis(l2, c) * q_test(l, c) *
                    (1.0 / local_permeability_ratio +
                     rho_f_over_rho * lambda_sq * local_permeability *
                       timestepper_pt->weight(1, 0) / porosity) *
                    r * w * J;
                }
              }
            }

            // d(q_eqn_l)/d(P_l2)
            for (unsigned l2 = 0; l2 < n_p_basis; l2++)
            {
              local_unknown = p_local_eqn(l2);

              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) -=
                  p_basis(l2) * (div_q_test_ds(l) * r + q_test(l, 0) * J) * w;
              }
            }
          } // End of Jacobian entries
        } // End of if not boundary condition
      } // End of loop over q test functions

      // loop over pressure test functions
      for (unsigned l = 0; l < n_p_basis; l++)
      {
        // get the local equation number
        local_eqn = p_local_eqn(l);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // solid divergence term
          residuals[local_eqn] +=
            alpha * interpolated_div_du_dt_dx * p_test(l) * r * w * J;
          residuals[local_eqn] +=
            alpha * interpolated_du_r_dt * p_test(l) * w * J;
          // deliberately no jacobian factor in this integral
          residuals[local_eqn] +=
            local_permeability * interpolated_div_q_ds * p_test(l) * r * w;
          // deliberately no r factor in this integral
          residuals[local_eqn] +=
            local_permeability * interpolated_q[0] * p_test(l) * w * J;
          residuals[local_eqn] -= mass_source_local * p_test(l) * r * w * J;

          // Jacobian entries
          if (flag)
          {
            // d(p_eqn_l)/d(U_l2,c)
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              for (unsigned c = 0; c < 2; c++)
              {
                local_unknown =
                  this->nodal_local_eqn(l2, u_index_axisym_poroelasticity(c));

                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    alpha * this->node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                    du_basis_dx(l2, c) * p_test(l) * r * w * J;

                  // Extra term due to cylindrical coordinate system
                  if (c == 0)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      alpha *
                      this->node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                      u_basis(l2) * p_test(l) * w * J;
                  }
                }
              }
            }

            // d(p_eqn_l)/d(Q_l2)
            for (unsigned l2 = 0; l2 < n_q_basis; l2++)
            {
              if (l2 < n_q_basis_edge)
              {
                local_unknown = q_edge_local_eqn(l2);
              }
              else // n_q_basis_edge <= l < n_basis
              {
                local_unknown = q_internal_local_eqn(l2 - n_q_basis_edge);
              }

              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  (div_q_basis_ds(l2) * r + q_basis(l2, 0) * J) *
                  local_permeability * p_test(l) * w;
              }
            }
          } // End of Jacobian entries
        } // End of if not boundary condition
      } // End of loop over p test functions
    } // End of loop over integration points
  }


} // namespace oomph
