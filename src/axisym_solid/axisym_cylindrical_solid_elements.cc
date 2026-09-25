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

// Non-inline functions for axisymmetric solid mechanics elements

#include "axisym_cylindrical_solid_elements.h"
#include "axisym_cylindrical_solid_with_pressure_elements.h"

namespace oomph
{
  /// Static default value for the timescale ratio
  double AxisymmetricCylindricalPVDEquations::Default_lambda_sq_value = 1.0;

  /// Static default value for the damping parameters
  double AxisymmetricCylindricalPVDEquations::Default_eta_value = 0.0;

  /// Static default value for the timescale ratio
  double
    AxisymmetricCylindricalPVDWithPressureEquations::Default_lambda_sq_value =
      1.0;

  /// Static default value for the damping parameters
  double AxisymmetricCylindricalPVDWithPressureEquations::Default_eta_value =
    0.0;

  //====================================================================
  /// Residuals for AxisymCylindricalPVDWithPressureElements
  //====================================================================
  void AxisymmetricCylindricalPVDEquations::
    fill_in_contribution_to_residuals_axisym_pvd(Vector<double>& residuals,
                                                 DenseMatrix<double>& jacobian,
                                                 const unsigned& flag)
  {
    // Set the number of Lagrangian coordinates
    unsigned n_lagrangian = 2;
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Integers to store local equation number and local unknown
    int local_eqn = 0, local_unknown = 0;

    // Set up memory for the shape functions
    Shape psi(n_node);
    DShape dpsidxi(n_node, n_lagrangian);

    // Timescale ratio (non-dim density)
    double lambda_sq = this->lambda_sq();

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Get the mass damping parameter
    const double eta_M = eta_mass();

    // Time factor
    double time_factor = 0.0;
    if (lambda_sq > 0)
    {
      time_factor = this->node_pt(0)->position_time_stepper_pt()->weight(2, 0);
    }

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);
      // Call the derivatives of the shape functions
      double J = dshape_lagrangian_at_knot(ipt, psi, dpsidxi);
      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate the local Lagrangian coordinates, position components
      // and the derivatives of global position components
      // wrt lagrangian coordinates, as well as acceleration
      Vector<double> interpolated_xi(2, 0.0);
      Vector<double> interpolated_X(2, 0.0);
      Vector<double> accel(2, 0.0);
      Vector<double> interpolated_dXdt(2, 0.0);
      DenseMatrix<double> interpolated_dXdxi(2, 2, 0.0);

      // Calculate displacements and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over displacement components (deformed position)
        for (unsigned i = 0; i < 2; i++)
        {
          // Set the value of the lagrangian coordinate
          interpolated_xi[i] += lagrangian_position(l, i) * psi(l);
          // Set the value of the position component
          interpolated_X[i] += nodal_position(l, i) * psi(l);
          // Set accel.
          accel[i] += dnodal_position_dt(l, 2, i) * psi(l);
          // Set velocity
          interpolated_dXdt[i] += dnodal_position_dt(l, 1, i) * psi(l);
          // Loop over Lagrangian derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            // Calculate dX[i]/dxi_{j}
            interpolated_dXdxi(i, j) += nodal_position(l, i) * dpsidxi(l, j);
          }
        }
      }

      // We are now in a position to calculate the undeformed metric tensor
      DenseMatrix<double> g(3);
      // r row
      g(0, 0) = 1.0;
      g(0, 1) = 0.0;
      g(0, 2) = 0.0;
      // z row
      g(1, 0) = 0.0;
      g(1, 1) = 1.0;
      g(1, 2) = 0.0;
      // phi row
      g(2, 0) = 0.0;
      g(2, 1) = 0.0;
      g(2, 2) = interpolated_xi[0] * interpolated_xi[0];

      // Now multiply the weight by the square-root of the undeformed metric
      // tensor r
      double detg = g(0, 0) * g(1, 1) * g(2, 2);
      W *= sqrt(detg);

      // Now calculate the deformed metric tensor
      DenseMatrix<double> G(3);
      // r row
      G(0, 0) = interpolated_dXdxi(0, 0) * interpolated_dXdxi(0, 0) +
                interpolated_dXdxi(1, 0) * interpolated_dXdxi(1, 0);
      G(0, 1) = interpolated_dXdxi(0, 0) * interpolated_dXdxi(0, 1) +
                interpolated_dXdxi(1, 0) * interpolated_dXdxi(1, 1);
      G(0, 2) = 0.0;
      // z row
      G(1, 0) = G(0, 1);
      G(1, 1) = interpolated_dXdxi(0, 1) * interpolated_dXdxi(0, 1) +
                interpolated_dXdxi(1, 1) * interpolated_dXdxi(1, 1);
      G(1, 2) = 0.0;
      // phi row
      G(2, 0) = 0.0;
      G(2, 1) = 0.0;
      G(2, 2) = interpolated_X[0] * interpolated_X[0];

      // Now calculate the stress tensor from the constitutive law
      DenseMatrix<double> sigma(3, 3, 0.0);
      get_stress(g, G, sigma);

      // If we're calculating the Jacobian, will need derivative of stress
      // tensor w.r.t. the deformed metric tensor
      RankFourTensor<double> d_stress_dG(3, 3, 3, 3, 0.0);
      RankFourTensor<double> dG_dX(n_node, 2, 3, 3, 0.0);

      if (flag == 1)
      {
        // Get the "upper triangular"
        // entries of the derivatives of the stress tensor with
        // respect to G
        this->get_d_stress_dG_upper(g, G, sigma, d_stress_dG);

        // Construct upper triangle of dGdX
        // Loop over nodes
        for (unsigned m = 0; m < n_node; m++)
        {
          // Loops over directions
          for (unsigned i = 0; i < 2; i++)
          {
            for (unsigned a = 0; a < 2; a++)
            {
              for (unsigned b = a; b < 2; b++)
              {
                dG_dX(m, i, a, b) = interpolated_dXdxi(i, a) * dpsidxi(m, b) +
                                    interpolated_dXdxi(i, b) * dpsidxi(m, a);
              }
            }
          }

          // Accounting for axisymmetric coord system
          dG_dX(m, 0, 2, 2) = 2.0 * interpolated_X[0] * psi(m);
        }
      }

      // Get body force at current time
      Vector<double> b(2, 0.0);
      this->body_force(interpolated_xi, b);

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over displacement components
        for (unsigned i = 0; i < 2; i++)
        {
          // Get the local eqn
          local_eqn = position_local_eqn(l, 0, i);

          // If not a boundary condition
          if (local_eqn >= 0)
          {
            // Forcing/inertial contributions
            residuals[local_eqn] +=
              (lambda_sq * (accel[i] + eta_M * interpolated_dXdt[i]) - b[i]) *
              psi(l) * W;

            // Stress term
            for (unsigned a = 0; a < 2; a++)
            {
              for (unsigned b = 0; b < 2; b++)
              {
                residuals[local_eqn] +=
                  sigma(a, b) * dpsidxi(l, b) * interpolated_dXdxi(i, a) * W;
              }
            }

            // Additional stress term if it's the r component
            if (i == 0)
            {
              residuals[local_eqn] +=
                sigma(2, 2) * interpolated_X[0] * psi(l) * W;
            }

            // Get Jacobian too?
            if (flag == 1)
            {
              // Loop over the nodes of the element again
              for (unsigned ll = 0; ll < n_node; ll++)
              {
                // Loop over the displacement components again
                for (unsigned ii = 0; ii < 2; ii++)
                {
                  // Get the number of the unknown
                  local_unknown = position_local_eqn(ll, 0, ii);

                  /*IF it's not a boundary condition*/
                  if (local_unknown >= 0)
                  {
                    // General stress term
                    //--------------------
                    double sum = 0.0;
                    for (unsigned a = 0; a < 2; a++)
                    {
                      for (unsigned b = a; b < 2; b++)
                      {
                        double factor = dG_dX(l, i, a, b);
                        if (a == b) factor *= 0.5;

                        for (unsigned aa = 0; aa < 3; aa++)
                        {
                          // Only upper half of derivatives w.r.t.
                          // symm tensor
                          for (unsigned bb = aa; bb < 3; bb++)
                          {
                            sum += factor * d_stress_dG(a, b, aa, bb) *
                                   dG_dX(ll, ii, aa, bb);
                          }
                        }
                      }
                    }

                    // Contribution accounting for axisymmetry
                    if (i == 0)
                    {
                      double factor = 0.5 * dG_dX(l, i, 2, 2);
                      for (unsigned aa = 0; aa < 3; aa++)
                      {
                        // Only upper half of derivatives w.r.t.
                        // symm tensor
                        for (unsigned bb = aa; bb < 3; bb++)
                        {
                          sum += factor * d_stress_dG(2, 2, aa, bb) *
                                 dG_dX(ll, ii, aa, bb);
                        }
                      }
                    }

                    // Multiply by weight and add contribution
                    // (Add directly because this bit is nonsymmetric)
                    jacobian(local_eqn, local_unknown) += sum * W;

                    // Only upper triangle (no separate test for bc as
                    // local_eqn is already nonnegative). Can be done
                    // for remaining terms as they are symmetric
                    if ((i == ii) && (local_unknown >= local_eqn))
                    {
                      // Initialise contribution
                      double sum = 0.0;

                      // Inertia term
                      sum += lambda_sq * time_factor * psi(ll) * psi(l);

                      // Stress term
                      for (unsigned a = 0; a < 2; a++)
                      {
                        for (unsigned b = 0; b < 2; b++)
                        {
                          sum += sigma(a, b) * dpsidxi(ll, a) * dpsidxi(l, b);
                        }
                      }

                      // Accounting for axisymmetry
                      if (i == 0)
                      {
                        sum += sigma(2, 2) * psi(l) * psi(ll);
                      }

                      // Multiply by weights to form contribution
                      double sym_entry = sum * W;
                      // Add contribution to jacobian
                      jacobian(local_eqn, local_unknown) += sym_entry;
                      // Add to lower triangular entries
                      if (local_eqn != local_unknown)
                      {
                        jacobian(local_unknown, local_eqn) += sym_entry;
                      }
                    }
                  } // End of if not boundary condition
                } // End of second displacement component loop
              } // End of second node loop
            } // End of Jacobian
          } // End of if not a boundary condition
        } // End of loop over displacement components
      } // End of loop over test functions
    } // End of loop over integration points
  }

  //====================================================================
  /// Residuals for Refineable AxisymCylindricalPVDWithPressureElements
  //====================================================================
  void AxisymmetricCylindricalPVDWithPressureEquations::
    fill_in_contribution_to_residuals_axisym_pvd_with_pressure(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
  {
    // Set the number of Lagrangian coordinates
    unsigned n_lagrangian = 2;
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Integers to store local equation number and local unknown
    int local_eqn = 0, local_unknown = 0;

    // Find out how many pressure dofs there are
    const unsigned n_solid_pres = this->nsolid_pres();

    // Set up memory for the shape functions
    Shape psi(n_node);
    DShape dpsidxi(n_node, n_lagrangian);

    // Set up memory for the pressure shape functions
    Shape psisp(n_solid_pres);

    // Timescale ratio (non-dim density)
    double lambda_sq = this->lambda_sq();

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Get the mass damping parameter
    const double eta_M = eta_mass();

    // Time factor
    double time_factor = 0.0;
    if (lambda_sq > 0)
    {
      time_factor = this->node_pt(0)->position_time_stepper_pt()->weight(2, 0);
    }

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);
      // Call the derivatives of the shape functions
      double J = dshape_lagrangian_at_knot(ipt, psi, dpsidxi);
      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Call the pressure shape functions
      this->solid_pshape_at_knot(ipt, psisp);

      // Calculate the local Lagrangian coordinates, position components
      // and the derivatives of global position components
      // wrt lagrangian coordinates, as well as acceleration
      Vector<double> interpolated_xi(2, 0.0);
      Vector<double> interpolated_X(2, 0.0);
      Vector<double> accel(2, 0.0);
      Vector<double> interpolated_dXdt(2, 0.0);
      DenseMatrix<double> interpolated_dXdxi(2, 2, 0.0);
      double interpolated_solid_p = 0.0;

      // Calculate displacements and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over displacement components (deformed position)
        for (unsigned i = 0; i < 2; i++)
        {
          // Set the value of the lagrangian coordinate
          interpolated_xi[i] += lagrangian_position(l, i) * psi(l);
          // Set the value of the position component
          interpolated_X[i] += nodal_position(l, i) * psi(l);
          // Set accel.
          accel[i] += dnodal_position_dt(l, 2, i) * psi(l);
          // Set velocity
          interpolated_dXdt[i] += dnodal_position_dt(l, 1, i) * psi(l);
          // Loop over Lagrangian derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            // Calculate dX[i]/dxi_{j}
            interpolated_dXdxi(i, j) += nodal_position(l, i) * dpsidxi(l, j);
          }
        }
      }

      // Calculate the local internal pressure
      for (unsigned l = 0; l < n_solid_pres; l++)
      {
        interpolated_solid_p += solid_p(l) * psisp(l);
      }

      // We are now in a position to calculate the undeformed metric tensor
      DenseMatrix<double> g(3);
      // r row
      g(0, 0) = 1.0;
      g(0, 1) = 0.0;
      g(0, 2) = 0.0;
      // z row
      g(1, 0) = 0.0;
      g(1, 1) = 1.0;
      g(1, 2) = 0.0;
      // phi row
      g(2, 0) = 0.0;
      g(2, 1) = 0.0;
      g(2, 2) = interpolated_xi[0] * interpolated_xi[0];

      // Now multiply the weight by the square-root of the undeformed metric
      // tensor r
      double detg = g(0, 0) * g(1, 1) * g(2, 2);
      W *= sqrt(detg);

      // Now calculate the deformed metric tensor
      DenseMatrix<double> G(3);
      // r row
      G(0, 0) = interpolated_dXdxi(0, 0) * interpolated_dXdxi(0, 0) +
                interpolated_dXdxi(1, 0) * interpolated_dXdxi(1, 0);
      G(0, 1) = interpolated_dXdxi(0, 0) * interpolated_dXdxi(0, 1) +
                interpolated_dXdxi(1, 0) * interpolated_dXdxi(1, 1);
      G(0, 2) = 0.0;
      // z row
      G(1, 0) = G(0, 1);
      G(1, 1) = interpolated_dXdxi(0, 1) * interpolated_dXdxi(0, 1) +
                interpolated_dXdxi(1, 1) * interpolated_dXdxi(1, 1);
      G(1, 2) = 0.0;
      // phi row
      G(2, 0) = 0.0;
      G(2, 1) = 0.0;
      G(2, 2) = interpolated_X[0] * interpolated_X[0];

      // Now calculate the deviatoric stress tensor from the constitutive law
      DenseMatrix<double> sigma_dev(3), Gup(3);
      double detG = 0.0, gen_dil = 0.0, inv_kappa = 0.0;
      // If it's incompressible call one form of the constitutive law
      if (Incompressible)
      {
        get_stress(g, G, sigma_dev, Gup, detG);
      }
      // Otherwise call another form
      else
      {
        get_stress(g, G, sigma_dev, Gup, gen_dil, inv_kappa);
      }

      // Build the stress tensor up from its pressure and deviatoric
      // components
      DenseMatrix<double> sigma(3, 3, 0.0);
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = 0; j < 3; j++)
        {
          sigma(i, j) =
            -1.0 * interpolated_solid_p * Gup(i, j) + sigma_dev(i, j);
        }
      }

      // If we're calculating the Jacobian, will need derivative of stress
      // tensor w.r.t. the deformed metric tensor
      RankFourTensor<double> d_stress_dG(3, 3, 3, 3, 0.0);
      RankFourTensor<double> dG_dX(n_node, 2, 3, 3, 0.0);
      DenseMatrix<double> d_detG_dG(3, 3, 0.0);
      DenseMatrix<double> d_gen_dil_dG(3, 3, 0.0);

      if (flag == 1)
      {
        // If incompressible, call the incompressible form
        if (Incompressible)
        {
          this->get_d_stress_dG_upper(
            g, G, sigma, detG, interpolated_solid_p, d_stress_dG, d_detG_dG);
        }
        else
        // Otherwise call the near-incompressible form
        {
          this->get_d_stress_dG_upper(g,
                                      G,
                                      sigma,
                                      gen_dil,
                                      inv_kappa,
                                      interpolated_solid_p,
                                      d_stress_dG,
                                      d_gen_dil_dG);
        }

        // Construct upper triangle of dGdX
        // Loop over nodes
        for (unsigned m = 0; m < n_node; m++)
        {
          // Loops over directions
          for (unsigned i = 0; i < 2; i++)
          {
            for (unsigned a = 0; a < 2; a++)
            {
              for (unsigned b = a; b < 2; b++)
              {
                dG_dX(m, i, a, b) = interpolated_dXdxi(i, a) * dpsidxi(m, b) +
                                    interpolated_dXdxi(i, b) * dpsidxi(m, a);
              }
            }
          }

          // Accounting for axisymmetric coord system
          dG_dX(m, 0, 2, 2) = 2.0 * interpolated_X[0] * psi(m);
        }
      }

      // Get body force at current time
      Vector<double> b(2, 0.0);
      this->body_force(interpolated_xi, b);

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over displacement components
        for (unsigned i = 0; i < 2; i++)
        {
          // Get the local eqn
          local_eqn = position_local_eqn(l, 0, i);

          // If not a boundary condition
          if (local_eqn >= 0)
          {
            // Forcing/inertial contributions
            residuals[local_eqn] +=
              (lambda_sq * (accel[i] + eta_M * interpolated_dXdt[i]) - b[i]) *
              psi(l) * W;

            // Stress term
            for (unsigned a = 0; a < 2; a++)
            {
              for (unsigned b = 0; b < 2; b++)
              {
                residuals[local_eqn] +=
                  sigma(a, b) * dpsidxi(l, b) * interpolated_dXdxi(i, a) * W;
              }
            }

            // Additional stress term if it's the r component
            if (i == 0)
            {
              residuals[local_eqn] +=
                sigma(2, 2) * interpolated_X[0] * psi(l) * W;
            }

            // Get Jacobian too?
            if (flag == 1)
            {
              // Loop over the nodes of the element again
              for (unsigned ll = 0; ll < n_node; ll++)
              {
                // Loop over the displacement components again
                for (unsigned ii = 0; ii < 2; ii++)
                {
                  // Get the number of the unknown
                  local_unknown = position_local_eqn(ll, 0, ii);

                  /*IF it's not a boundary condition*/
                  if (local_unknown >= 0)
                  {
                    // General stress term
                    //--------------------
                    double sum = 0.0;
                    for (unsigned a = 0; a < 2; a++)
                    {
                      for (unsigned b = a; b < 2; b++)
                      {
                        double factor = dG_dX(l, i, a, b);
                        if (a == b) factor *= 0.5;

                        for (unsigned aa = 0; aa < 3; aa++)
                        {
                          // Only upper half of derivatives w.r.t.
                          // symm tensor
                          for (unsigned bb = aa; bb < 3; bb++)
                          {
                            sum += factor * d_stress_dG(a, b, aa, bb) *
                                   dG_dX(ll, ii, aa, bb);
                          }
                        }
                      }
                    }

                    // Contribution accounting for axisymmetry
                    if (i == 0)
                    {
                      double factor = 0.5 * dG_dX(l, i, 2, 2);
                      for (unsigned aa = 0; aa < 3; aa++)
                      {
                        // Only upper half of derivatives w.r.t.
                        // symm tensor
                        for (unsigned bb = aa; bb < 3; bb++)
                        {
                          sum += factor * d_stress_dG(2, 2, aa, bb) *
                                 dG_dX(ll, ii, aa, bb);
                        }
                      }
                    }

                    // Multiply by weight and add contribution
                    // (Add directly because this bit is nonsymmetric)
                    jacobian(local_eqn, local_unknown) += sum * W;

                    // Only upper triangle (no separate test for bc as
                    // local_eqn is already nonnegative). Can be done
                    // for remaining terms as they are symmetric
                    if ((i == ii) && (local_unknown >= local_eqn))
                    {
                      // Initialise contribution
                      double sum = 0.0;

                      // Inertia term
                      sum += lambda_sq * time_factor * psi(ll) * psi(l);

                      // Stress term
                      for (unsigned a = 0; a < 2; a++)
                      {
                        for (unsigned b = 0; b < 2; b++)
                        {
                          sum += sigma(a, b) * dpsidxi(ll, a) * dpsidxi(l, b);
                        }
                      }

                      // Accounting for axisymmetry
                      if (i == 0)
                      {
                        sum += sigma(2, 2) * psi(l) * psi(ll);
                      }

                      // Multiply by weights to form contribution
                      double sym_entry = sum * W;
                      // Add contribution to jacobian
                      jacobian(local_eqn, local_unknown) += sym_entry;
                      // Add to lower triangular entries
                      if (local_eqn != local_unknown)
                      {
                        jacobian(local_unknown, local_eqn) += sym_entry;
                      }
                    }
                  } // End of if not boundary condition
                }
              }

              // Can add in the pressure jacobian terms
              // Loop over the pressure nodes
              for (unsigned l2 = 0; l2 < n_solid_pres; l2++)
              {
                local_unknown = this->solid_p_local_eqn(l2);

                // If it's not a boundary condition
                if (local_unknown >= 0)
                {
                  // Add the pressure terms to the jacobian
                  for (unsigned a = 0; a < 2; a++)
                  {
                    for (unsigned b = 0; b < 2; b++)
                    {
                      jacobian(local_eqn, local_unknown) -=
                        psisp(l2) * Gup(a, b) * interpolated_dXdxi(i, a) *
                        dpsidxi(l, b) * W;
                    }
                  }

                  // Additional contribution for axisymmetry
                  if (i == 0)
                  {
                    jacobian(local_eqn, local_unknown) -=
                      psisp(l2) * Gup(2, 2) * interpolated_X[0] * psi(l) * W;
                  }
                }
              } // End of loop over pressure dofs
            } // End of Jacobian
          }
        } // End of loop over displacement components
      } // End of loop over test functions

      // Now loop over the pressure dofs
      for (unsigned l = 0; l < n_solid_pres; l++)
      {
        local_eqn = this->solid_p_local_eqn(l);

        // Pinned (unlikely, actually) or real dof?
        if (local_eqn >= 0)
        {
          // For true incompressibility we need to conserve volume
          // so the determinant of the deformed metric tensor
          // needs to be equal to that of the undeformed one, which
          // is equal to the volumetric growth factor
          if (this->Incompressible)
          {
            residuals[local_eqn] += (detG / detg - 1.0) * psisp[l] * W;

            // Add in Jacobian terms
            if (flag == 1)
            {
              // Loop over the nodes of the element again
              for (unsigned ll = 0; ll < n_node; ll++)
              {
                // Loop over the displacement components again
                for (unsigned ii = 0; ii < 2; ii++)
                {
                  // Get the number of the unknown
                  int local_unknown = this->position_local_eqn(ll, 0, ii);

                  /*IF it's not a boundary condition*/
                  if (local_unknown >= 0)
                  {
                    // General stress term
                    //--------------------
                    double sum = 0.0;

                    for (unsigned aa = 0; aa < 3; aa++)
                    {
                      // Only upper half of derivatives w.r.t.
                      // symm tensor
                      for (unsigned bb = aa; bb < 3; bb++)
                      {
                        sum += d_detG_dG(aa, bb) * dG_dX(ll, ii, aa, bb) *
                               psisp(l) / detg;
                      }
                    }
                    jacobian(local_eqn, local_unknown) += sum * W;
                  }
                }
              }
              // No Jacobian terms due to pressure since it does not feature
              // in the incompressibility constraint
            }
          }
          else
          {
            residuals[local_eqn] +=
              (inv_kappa * interpolated_solid_p + gen_dil) * psisp[l] * W;

            // Add in the jacobian terms
            if (flag == 1)
            {
              // Loop over the nodes of the element again
              for (unsigned ll = 0; ll < n_node; ll++)
              {
                // Loop over the displacement components again
                for (unsigned ii = 0; ii < 2; ii++)
                {
                  // Get the number of the unknown
                  local_unknown = position_local_eqn(ll, 0, ii);

                  /*IF it's not a boundary condition*/
                  if (local_unknown >= 0)
                  {
                    // General stress term
                    double sum = 0.0;
                    for (unsigned aa = 0; aa < 3; aa++)
                    {
                      // Only upper half
                      for (unsigned bb = aa; bb < 3; bb++)
                      {
                        sum += d_gen_dil_dG(aa, bb) * dG_dX(ll, ii, aa, bb) *
                               psisp(l);
                      }
                    }
                    jacobian(local_eqn, local_unknown) += sum * W;
                  }
                }
              }

              // Loop over the pressure nodes again
              for (unsigned l2 = 0; l2 < n_solid_pres; l2++)
              {
                local_unknown = this->solid_p_local_eqn(l2);

                // If it's not a boundary condition
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    inv_kappa * psisp(l2) * psisp(l) * W;
                }
              } // End of loop over pressure dofs
            } // End of Jacobian
          }
        } // End of if not boundary condition
      } // End of loop over pressure dofs
    } // End of loop over integration points
  }

  //====================================================================
  /// Data for the number of Variables at each node
  //====================================================================
  const unsigned
    QAxisymmetricCylindricalPVDWithPressureElement::Initial_Nvalue[9] = {
      1, 0, 1, 0, 0, 0, 1, 0, 1};

  //==========================================================================
  /// Conversion from pressure dof to Node number at which pressure is stored
  //==========================================================================
  const unsigned QAxisymmetricCylindricalPVDWithPressureElement::Pconv[4] = {
    0, 2, 6, 8};

  //=======================================================================
  /// Data for the number of variables at each node
  //=======================================================================
  const unsigned TAxisymCylindricalPVDWithPressureElement::Initial_Nvalue[6] = {
    1, 1, 1, 0, 0, 0};

  //=======================================================================
  /// Data for the pressure conversion array
  //=======================================================================
  const unsigned TAxisymCylindricalPVDWithPressureElement::Pconv[3] = {0, 1, 2};

  //==================================================================
  /// Solid pressure shape function evaluated at integration point
  //==================================================================
  void AxisymmetricCylindricalPVDWithPressureEquations::solid_pshape_at_knot(
    const unsigned& ipt, Shape& psi) const
  {
    // Storage for local coordinates of the integration point
    Vector<double> s(2);
    // Set the local coordinates
    for (unsigned i = 0; i < 2; i++)
    {
      s[i] = this->integral_pt()->knot(ipt, i);
    }
    // Get the shape function
    solid_pshape(s, psi);
  }

  // Build required templates
  template class TAxisymCylindricalPVDElement<3>;
  template class QAxisymmetricCylindricalPVDElement<3>;

} // namespace oomph
