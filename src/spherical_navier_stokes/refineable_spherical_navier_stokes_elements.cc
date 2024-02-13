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
#include "refineable_spherical_navier_stokes_elements.h"


namespace oomph
{
  //========================================================================
  /// Add element's contribution to the elemental
  /// residual vector and/or Jacobian matrix.
  /// flag=1: compute both
  /// flag=0: compute only residual vector
  //=======================================================================
  void RefineableSphericalNavierStokesEquations::
    fill_in_generic_residual_contribution_spherical_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Find out how many pressure dofs there are
    const unsigned n_pres = npres_spherical_nst();

    // Get the local indices of the nodal coordinates
    unsigned u_nodal_index[3];
    for (unsigned i = 0; i < 3; ++i)
    {
      u_nodal_index[i] = this->u_index_spherical_nst(i);
    }

    // Which nodal value represents the pressure? (Negative if pressure
    // is not based on nodal interpolation).
    const int p_index = this->p_nodal_index_spherical_nst();

    // Local array of booleans that are true if the l-th pressure value is
    // hanging (avoid repeated virtual function calls)
    bool pressure_dof_is_hanging[n_pres];
    // If the pressure is stored at a node
    if (p_index >= 0)
    {
      // Read out whether the pressure is hanging
      for (unsigned l = 0; l < n_pres; ++l)
      {
        pressure_dof_is_hanging[l] = pressure_node_pt(l)->is_hanging(p_index);
      }
    }
    // Otherwise the pressure is not stored at a node and so cannot hang
    else
    {
      for (unsigned l = 0; l < n_pres; ++l)
      {
        pressure_dof_is_hanging[l] = false;
      }
    }


    // Set up memory for the shape and test functions
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, 2), dtestfdx(n_node, 2);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Get Physical Variables from Element
    // Reynolds number must be multiplied by the density ratio
    const double dens_ratio = density_ratio();
    const double scaled_re = re() * dens_ratio;
    const double scaled_re_st = re_st() * dens_ratio;
    const double scaled_re_inv_ro = re_invro() * dens_ratio;
    // const double scaled_re_inv_fr = re_invfr()*dens_ratio;
    // const double visc_ratio = viscosity_ratio();
    Vector<double> G = g();

    // Integers to store the local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0;

    // Local storage for pointers to hang info objects
    HangInfo *hang_info_pt = 0, *hang_info2_pt = 0;

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

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_spherical_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Call the pressure shape and test functions
      pshape_spherical_nst(s, psip, testp);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the pressure and velocity components
      //--------------------------------------------------------------
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(3, 0.0);
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> mesh_velocity(2, 0.0);
      Vector<double> dudt(3, 0.0);
      DenseMatrix<double> interpolated_dudx(3, 2, 0.0);

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
      {
        interpolated_p += this->p_spherical_nst(l) * psip(l);
      }

      // Calculate velocities and derivatives

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the shape function
        double psi_ = psif(l);
        // Loop over positions separately
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += nodal_position(l, i) * psi_;
        }

        // Loop over velocity directions (three of these)
        for (unsigned i = 0; i < 3; i++)
        {
          const double u_value = this->nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psi_;
          dudt[i] += du_dt_spherical_nst(l, i) * psi_;

          // Loop over derivative directions for gradients
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }

        // Only bother to calculate the mesh velocity if we are using an ALE
        // method
        if (!ALE_is_disabled)
        {
          throw OomphLibError(
            "ALE is not properly implemented for Refineable Spherical NS yet",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);

          // Loop over directions (only DIM (2) of these)
          for (unsigned i = 0; i < 2; i++)
          {
            mesh_velocity[i] += this->dnodal_position_dt(l, i) * psi_;
          }
        }
      } // End of loop over the nodes

      // Get the user-defined body force terms
      Vector<double> body_force(3);
      this->get_body_force_spherical_nst(
        time, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      // double
      // source=this->get_source_spherical_nst(time(),ipt,interpolated_x);

      // r is the first postition component
      const double r = interpolated_x[0];
      const double r2 = r * r;
      // const double theta = interpolated_x[1];
      const double sin_theta = sin(interpolated_x[1]);
      const double cos_theta = cos(interpolated_x[1]);
      const double cot_theta = cos_theta / sin_theta;

      const double u_r = interpolated_u[0];
      const double u_theta = interpolated_u[1];
      const double u_phi = interpolated_u[2];

      // Pre-calculate the scaling factor of the jacobian
      // double W_pure = W;

      // W *= r*r*sin(theta);

      // MOMENTUM EQUATIONS
      //==================
      // Number of master nodes and storage for the weight of the shape function
      unsigned n_master = 1;
      double hang_weight = 1.0;

      // Loop over the nodes for the test functions/equations
      //----------------------------------------------------
      for (unsigned l = 0; l < n_node; l++)
      {
        // Local boolean that indicates the hanging status of the node
        const bool is_node_hanging = node_pt(l)->is_hanging();

        // If the node is hanging
        if (is_node_hanging)
        {
          // Get the hanging pointer
          hang_info_pt = node_pt(l)->hanging_pt();
          // Read out number of master nodes from hanging data
          n_master = hang_info_pt->nmaster();
        }
        // Otherwise the node is its own master
        else
        {
          n_master = 1;
        }

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Loop over velocity components for equations
          for (unsigned i = 0; i < 2 + 1; i++)
          {
            // Get the equation number
            // If the node is hanging
            if (is_node_hanging)
            {
              // Get the equation number from the master node
              local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                               u_nodal_index[i]);
              // Get the hang weight from the master node
              hang_weight = hang_info_pt->master_weight(m);
            }
            // If the node is not hanging
            else
            {
              // Local equation number
              local_eqn = this->nodal_local_eqn(l, u_nodal_index[i]);
              // Node contributes with full weight
              hang_weight = 1.0;
            }

            // If it's not a boundary condition...
            if (local_eqn >= 0)
            {
              // initialise for residual calculation
              double sum = 0.0;

              switch (i)
              {
                  // RADIAL MOMENTUM EQUATION
                case 0:
                {
                  // Convective r-terms
                  double conv = r * u_r * interpolated_dudx(0, 0);

                  // Convective theta-terms
                  conv += u_theta * interpolated_dudx(0, 1);

                  // Remaining convective terms
                  conv -= (u_theta * u_theta + u_phi * u_phi);

                  // Add the time-derivative and convective terms
                  sum += (scaled_re_st * dudt[0] * r2 + scaled_re * r * conv) *
                         sin_theta * testf(l);

                  // Add the user-defined body force term
                  sum -= r2 * sin_theta * body_force[0] * testf(l);

                  // Add the Coriolis term
                  sum -= 2.0 * scaled_re_inv_ro * sin_theta * u_phi * r2 *
                         sin_theta * testf(l);

                  // r-derivative test function component of stress tensor
                  sum += (-interpolated_p + 2 * interpolated_dudx(0, 0)) * r2 *
                         sin_theta * dtestfdx(l, 0);

                  // theta-derivative test function component of stress tensor
                  sum += (r * interpolated_dudx(1, 0) - u_theta +
                          interpolated_dudx(0, 1)) *
                         sin_theta * dtestfdx(l, 1);

                  // Undifferentiated test function component of stress tensor
                  sum += 2.0 *
                         ((-r * interpolated_p + +interpolated_dudx(1, 1) +
                           2.0 * u_r) *
                            sin_theta +
                          u_theta * cos_theta) *
                         testf(l);
                }
                break;

                // THETA-COMPONENT MOMENTUM EQUATION
                case 1:
                {
                  // All convective terms
                  double conv =
                    (u_r * interpolated_dudx(1, 0) * r +
                     u_theta * interpolated_dudx(1, 1) + u_r * u_theta) *
                      sin_theta -
                    u_phi * u_phi * cos_theta;

                  // Add the time-derivative and convective terms to the
                  // residual
                  sum += (scaled_re_st * r2 * sin_theta * dudt[1] +
                          scaled_re * r * conv) *
                         testf(l);

                  // Add the user-defined body force term
                  sum -= r2 * sin_theta * body_force[1] * testf(l);

                  // Add the Coriolis term
                  sum -= 2.0 * scaled_re_inv_ro * cos_theta * u_phi * r2 *
                         sin_theta * testf(l);

                  // r-derivative test function component of stress tensor
                  sum += (r * interpolated_dudx(1, 0) - u_theta +
                          interpolated_dudx(0, 1)) *
                         r * sin_theta * dtestfdx(l, 0);

                  // theta-derivative test function component of stress tensor
                  sum += (-r * interpolated_p + 2.0 * interpolated_dudx(1, 1) +
                          2 * u_r) *
                         sin_theta * dtestfdx(l, 1);

                  // Undifferentiated test function component of stress tensor
                  sum -= ((r * interpolated_dudx(1, 0) - u_theta +
                           interpolated_dudx(0, 1)) *
                            sin_theta -
                          (-r * interpolated_p + 2.0 * u_r +
                           2.0 * u_theta * cot_theta) *
                            cos_theta) *
                         testf(l);
                }
                break;

                // PHI-COMPONENT MOMENTUM EQUATION
                case 2:

                {
                  // Convective r-terms
                  double conv = u_r * interpolated_dudx(2, 0) * r * sin_theta;

                  // Convective theta-terms
                  conv += u_theta * interpolated_dudx(2, 1) * sin_theta;

                  // Remaining convective terms
                  conv += u_phi * (u_r * sin_theta + u_theta * cos_theta);

                  // Add the time-derivative and convective terms
                  sum += (scaled_re_st * r2 * dudt[2] * sin_theta +
                          scaled_re * conv * r) *
                         testf(l);

                  sum -= r2 * sin_theta * body_force[2] * testf(l);

                  // Add the Coriolis term
                  sum += 2.0 * scaled_re_inv_ro *
                         (cos_theta * u_theta + sin_theta * u_r) * r2 *
                         sin_theta * testf(l);


                  // r-derivative test function component of stress tensor
                  sum += (r2 * interpolated_dudx(2, 0) - r * u_phi) *
                         dtestfdx(l, 0) * sin_theta;

                  // theta-derivative test function component of stress tensor
                  sum +=
                    (interpolated_dudx(2, 1) * sin_theta - u_phi * cos_theta) *
                    dtestfdx(l, 1);

                  // Undifferentiated test function component of stress tensor
                  sum -= ((r * interpolated_dudx(2, 0) - u_phi) * sin_theta +
                          (interpolated_dudx(2, 1) - u_phi * cot_theta) *
                            cos_theta) *
                         testf(l);
                }
                break;
              }

              // Add contribution
              // Sign changed to be consistent with other NS implementations
              residuals[local_eqn] -= sum * W * hang_weight;

              // CALCULATE THE JACOBIAN
              if (flag)
              {
                // Number of master nodes and weights
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                // Loop over the velocity nodes for columns
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Local boolean for hanging status
                  const bool is_node2_hanging = node_pt(l2)->is_hanging();

                  // If the node is hanging
                  if (is_node2_hanging)
                  {
                    hang_info2_pt = node_pt(l2)->hanging_pt();
                    // Read out number of master nodes from hanging data
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise the node is its own master
                  else
                  {
                    n_master2 = 1;
                  }

                  // Loop over the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    // Loop over the velocity components
                    for (unsigned i2 = 0; i2 < 2 + 1; i2++)
                    {
                      // Get the number of the unknown
                      // If the node is hanging
                      if (is_node2_hanging)
                      {
                        // Get the equation number from the master node
                        local_unknown = this->local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), u_nodal_index[i2]);
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      else
                      {
                        local_unknown =
                          this->nodal_local_eqn(l2, u_nodal_index[i2]);
                        hang_weight2 = 1.0;
                      }

                      // If the unknown is non-zero
                      if (local_unknown >= 0)
                      {
                        // Different results for i and i2
                        switch (i)
                        {
                            // RADIAL MOMENTUM EQUATION
                          case 0:
                            switch (i2)
                            {
                                // radial component
                              case 0:

                                // Add the mass matrix entries
                                if (flag == 2)
                                {
                                  mass_matrix(local_eqn, local_unknown) +=
                                    scaled_re_st * psif(l2) * testf(l) * r2 *
                                    sin_theta * W * hang_weight * hang_weight2;
                                }

                                {
                                  // Convective r-term contribution
                                  double jac_conv =
                                    r * (psif(l2) * interpolated_dudx(0, 0) +
                                         u_r * dpsifdx(l2, 0));

                                  // Convective theta-term contribution
                                  jac_conv += u_theta * dpsifdx(l2, 1);

                                  // Add the time-derivative contribution and
                                  // the convective
                                  // contribution to the sum
                                  double jac_sum =
                                    (scaled_re_st *
                                       node_pt(l2)->time_stepper_pt()->weight(
                                         1, 0) *
                                       psif(l2) * r2 +
                                     scaled_re * jac_conv * r) *
                                    testf(l);


                                  // Contribution from the r-derivative test
                                  // function
                                  // part of stress tensor
                                  jac_sum +=
                                    2.0 * dpsifdx(l2, 0) * dtestfdx(l, 0) * r2;

                                  // Contribution from the theta-derivative
                                  // test function part  of stress tensor
                                  jac_sum += dpsifdx(l2, 1) * dtestfdx(l, 1);


                                  // Contribution from the undifferentiated
                                  // test function part
                                  // of stress tensor
                                  jac_sum += 4.0 * psif[l2] * testf(l);

                                  // Add the total contribution to the
                                  // jacobian multiplied
                                  // by the jacobian weight
                                  jacobian(local_eqn, local_unknown) -=
                                    jac_sum * sin_theta * W * hang_weight *
                                    hang_weight2;
                                }

                                break;

                                // axial component
                              case 1:
                              {
                                // No time derivative contribution

                                // Convective contribution
                                double jac_sum =
                                  scaled_re *
                                  (interpolated_dudx(0, 1) - 2.0 * u_theta) *
                                  psif(l2) * r * sin_theta * testf(l);

                                // Contribution from the theta-derivative
                                // test function
                                // part of stress tensor
                                jac_sum += (r * dpsifdx(l2, 0) - psif(l2)) *
                                           dtestfdx(l, 1) * sin_theta;

                                // Contribution from the undifferentiated
                                // test function
                                // part of stress tensor
                                jac_sum += 2.0 *
                                           (dpsifdx(l2, 1) * sin_theta +
                                            psif(l2) * cos_theta) *
                                           testf(l);

                                // Add the full contribution to the jacobian
                                // matrix
                                jacobian(local_eqn, local_unknown) -=
                                  jac_sum * W * hang_weight * hang_weight2;

                              } // End of i2 = 1 section

                              break;

                                // azimuthal component
                              case 2:
                              {
                                // Single convective-term contribution
                                jacobian(local_eqn, local_unknown) +=
                                  2.0 * scaled_re * u_phi * psif[l2] * r *
                                  sin_theta * testf[l] * W * hang_weight *
                                  hang_weight2;

                                // Add the Coriolis term
                                jacobian(local_eqn, local_unknown) +=
                                  2.0 * scaled_re_inv_ro * sin_theta *
                                  psif(l2) * r2 * sin_theta * testf[l] * W *
                                  hang_weight * hang_weight2;
                              }

                              break;
                            } // End of contribution radial momentum eqn
                            break;

                            // AXIAL MOMENTUM EQUATION
                          case 1:
                            switch (i2)
                            {
                                // radial component
                              case 0:
                              {
                                // Convective terms
                                double jac_sum =
                                  scaled_re *
                                  (r2 * interpolated_dudx(1, 0) + r * u_theta) *
                                  psif(l2) * sin_theta * testf(l);

                                // Contribution from the r-derivative
                                // test function part of stress tensor
                                jac_sum += dpsifdx(l2, 1) * dtestfdx(l, 0) *
                                           sin_theta * r;

                                // Contribution from the theta-derivative
                                // test function
                                // part of stress tensor
                                jac_sum +=
                                  2.0 * psif(l2) * dtestfdx(l, 1) * sin_theta;

                                // Contribution from the undifferentiated
                                // test function
                                // part of stress tensor
                                jac_sum -= (dpsifdx(l2, 1) * sin_theta -
                                            2.0 * psif(l2) * cos_theta) *
                                           testf(l);

                                // Add the sum to the jacobian
                                jacobian(local_eqn, local_unknown) -=
                                  jac_sum * W * hang_weight * hang_weight2;
                              }

                              break;

                                // axial component
                              case 1:

                                // Add the mass matrix terms
                                if (flag == 2)
                                {
                                  mass_matrix(local_eqn, local_unknown) +=
                                    scaled_re_st * psif[l2] * testf[l] * W *
                                    r2 * sin_theta * hang_weight * hang_weight2;
                                }


                                {
                                  // Contribution from the convective terms
                                  double jac_conv =
                                    r * u_r * dpsifdx(l2, 0) +
                                    u_theta * dpsifdx(l2, 1) +
                                    (interpolated_dudx(1, 1) + u_r) * psif(l2);

                                  // Add the time-derivative term and the
                                  // convective terms
                                  double jac_sum =
                                    (scaled_re_st *
                                       node_pt(l2)->time_stepper_pt()->weight(
                                         1, 0) *
                                       psif(l2) * r2 +
                                     scaled_re * r * jac_conv) *
                                    testf(l) * sin_theta;


                                  // Contribution from the r-derivative test
                                  // function
                                  // part of stress tensor
                                  jac_sum += (r * dpsifdx(l2, 0) - psif(l2)) *
                                             r * dtestfdx(l, 0) * sin_theta;

                                  // Contribution from the theta-derivative
                                  // test function
                                  // part of stress tensor
                                  jac_sum += 2.0 * dpsifdx(l2, 1) *
                                             dtestfdx(l, 1) * sin_theta;

                                  // Contribution from the undifferentiated
                                  // test function
                                  // part of stress tensor
                                  jac_sum -=
                                    ((r * dpsifdx(l2, 0) - psif(l2)) *
                                       sin_theta -
                                     2.0 * psif(l2) * cot_theta * cos_theta) *
                                    testf(l);

                                  // Add the contribution to the jacobian
                                  jacobian(local_eqn, local_unknown) -=
                                    jac_sum * W * hang_weight * hang_weight2;
                                }

                                break;

                                // azimuthal component
                              case 2:
                              {
                                // Only a convective term contribution
                                jacobian(local_eqn, local_unknown) +=
                                  scaled_re * 2.0 * r * psif(l2) * u_phi *
                                  cos_theta * testf(l) * W * hang_weight *
                                  hang_weight2;

                                // Add the Coriolis term
                                jacobian(local_eqn, local_unknown) +=
                                  2.0 * scaled_re_inv_ro * cos_theta *
                                  psif(l2) * r2 * sin_theta * testf[l] * W *
                                  hang_weight * hang_weight2;
                              }

                              break;
                            }
                            break;

                            // AZIMUTHAL MOMENTUM EQUATION
                          case 2:
                            switch (i2)
                            {
                                // radial component
                              case 0:
                              {
                                // Contribution from convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re *
                                  (r * interpolated_dudx(2, 0) + u_phi) *
                                  psif(l2) * testf(l) * r * sin_theta * W *
                                  hang_weight * hang_weight2;

                                // Coriolis term
                                jacobian(local_eqn, local_unknown) -=
                                  2.0 * scaled_re_inv_ro * sin_theta *
                                  psif(l2) * r2 * sin_theta * testf[l] * W *
                                  hang_weight * hang_weight2;
                              }
                              break;

                              // axial component
                              case 1:
                              {
                                // Contribution from convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re *
                                  (interpolated_dudx(2, 1) * sin_theta +
                                   u_phi * cos_theta) *
                                  r * psif(l2) * testf(l) * W * hang_weight *
                                  hang_weight2;

                                // Coriolis term
                                jacobian(local_eqn, local_unknown) -=
                                  2.0 * scaled_re_inv_ro * cos_theta *
                                  psif(l2) * r2 * sin_theta * testf[l] * W *
                                  hang_weight * hang_weight2;
                              }

                              break;

                                // azimuthal component
                              case 2:

                                // Add the mass matrix terms
                                if (flag == 2)
                                {
                                  mass_matrix(local_eqn, local_unknown) +=
                                    scaled_re_st * psif[l2] * testf[l] * W *
                                    r2 * sin_theta * hang_weight * hang_weight2;
                                }

                                {
                                  // Convective terms
                                  double jac_conv =
                                    r * u_r * dpsifdx(l2, 0) * sin_theta;

                                  // Convective theta-term contribution
                                  jac_conv +=
                                    u_theta * dpsifdx(l2, 1) * sin_theta;

                                  // Contribution from the remaining convective
                                  // terms
                                  jac_conv +=
                                    (u_r * sin_theta + u_theta * cos_theta) *
                                    psif(l2);

                                  // Add the convective and time derivatives
                                  double jac_sum =
                                    (scaled_re_st *
                                       node_pt(l2)->time_stepper_pt()->weight(
                                         1, 0) *
                                       psif(l2) * r2 * sin_theta +
                                     scaled_re * r * jac_conv) *
                                    testf(l);


                                  // Contribution from the r-derivative test
                                  // function
                                  // part of stress tensor
                                  jac_sum += (r * dpsifdx(l2, 0) - psif(l2)) *
                                             dtestfdx(l, 0) * r * sin_theta;

                                  // Contribution from the theta-derivative
                                  // test function
                                  // part of stress tensor
                                  jac_sum += (dpsifdx(l2, 1) * sin_theta -
                                              psif(l2) * cos_theta) *
                                             dtestfdx(l, 1);

                                  // Contribution from the undifferentiated
                                  // test function
                                  // part of stress tensor
                                  jac_sum -=
                                    ((r * dpsifdx(l2, 0) - psif(l2)) *
                                       sin_theta +
                                     (dpsifdx(l2, 1) - psif(l2) * cot_theta) *
                                       cos_theta) *
                                    testf(l);

                                  // Add to the jacobian
                                  jacobian(local_eqn, local_unknown) -=
                                    jac_sum * W * hang_weight * hang_weight2;
                                }

                                break;
                            }
                            break;
                        }
                      }
                    }
                  }
                } // End of loop over the nodes


                // Loop over the pressure shape functions
                for (unsigned l2 = 0; l2 < n_pres; l2++)
                {
                  // If the pressure dof is hanging
                  if (pressure_dof_is_hanging[l2])
                  {
                    // Pressure dof is hanging so it must be nodal-based
                    hang_info2_pt = pressure_node_pt(l2)->hanging_pt(p_index);

                    // Get the number of master nodes from the pressure node
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise the node is its own master
                  else
                  {
                    n_master2 = 1;
                  }

                  // Loop over the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    // Get the number of the unknown
                    // If the pressure dof is hanging
                    if (pressure_dof_is_hanging[l2])
                    {
                      // Get the unknown from the node
                      local_unknown = local_hang_eqn(
                        hang_info2_pt->master_node_pt(m2), p_index);
                      // Get the weight from the hanging object
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    else
                    {
                      local_unknown = p_local_eqn(l2);
                      hang_weight2 = 1.0;
                    }

                    // If the unknown is not pinned
                    if (local_unknown >= 0)
                    {
                      // Add in contributions to different equations
                      switch (i)
                      {
                          // RADIAL MOMENTUM EQUATION
                        case 0:
                          jacobian(local_eqn, local_unknown) +=
                            psip(l2) *
                            (r2 * dtestfdx(l, 0) + 2.0 * r * testf[l]) * W *
                            sin_theta * hang_weight * hang_weight2;


                          break;

                          // AXIAL MOMENTUM EQUATION
                        case 1:
                          jacobian(local_eqn, local_unknown) +=
                            psip(l2) * r *
                            (dtestfdx(l, 1) * sin_theta +
                             cos_theta * testf(l)) *
                            W * hang_weight * hang_weight2;

                          break;

                          // AZIMUTHAL MOMENTUM EQUATION
                        case 2:
                          break;
                      }
                    }
                  }
                } // End of loop over pressure dofs
              } // End of Jacobian calculation
            } // End of if not boundary condition statement
          } // End of loop over velocity components
        } // End of loop over master nodes
      } // End of loop over nodes


      // CONTINUITY EQUATION
      //===================

      // Loop over the pressure shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        // If the pressure dof is hanging
        if (pressure_dof_is_hanging[l])
        {
          // Pressure dof is hanging so it must be nodal-based
          hang_info_pt = pressure_node_pt(l)->hanging_pt(p_index);
          // Get the number of master nodes from the pressure node
          n_master = hang_info_pt->nmaster();
        }
        // Otherwise the node is its own master
        else
        {
          n_master = 1;
        }

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Get the number of the unknown
          // If the pressure dof is hanging
          if (pressure_dof_is_hanging[l])
          {
            local_eqn =
              local_hang_eqn(hang_info_pt->master_node_pt(m), p_index);
            hang_weight = hang_info_pt->master_weight(m);
          }
          else
          {
            local_eqn = p_local_eqn(l);
            hang_weight = 1.0;
          }

          // If the equation is not pinned
          if (local_eqn >= 0)
          {
            // The entire continuity equation
            residuals[local_eqn] += ((2.0 * u_r + r * interpolated_dudx(0, 0) +
                                      interpolated_dudx(1, 1)) *
                                       sin_theta +
                                     u_theta * cos_theta) *
                                    r * testp(l) * W * hang_weight;

            // CALCULATE THE JACOBIAN
            //======================
            if (flag)
            {
              unsigned n_master2 = 1;
              double hang_weight2 = 1.0;
              // Loop over the velocity nodes for columns
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Local boolean to indicate whether the node is hanging
                bool is_node2_hanging = node_pt(l2)->is_hanging();

                // If the node is hanging
                if (is_node2_hanging)
                {
                  hang_info2_pt = node_pt(l2)->hanging_pt();
                  // Read out number of master nodes from hanging data
                  n_master2 = hang_info2_pt->nmaster();
                }
                // Otherwise the node is its own master
                else
                {
                  n_master2 = 1;
                }

                // Loop over the master nodes
                for (unsigned m2 = 0; m2 < n_master2; m2++)
                {
                  // Loop over the velocity components
                  for (unsigned i2 = 0; i2 < 2 + 1; i2++)
                  {
                    // Get the number of the unknown
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Get the equation number from the master node
                      local_unknown = local_hang_eqn(
                        hang_info2_pt->master_node_pt(m2), u_nodal_index[i2]);
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    else
                    {
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_index[i2]);
                      hang_weight2 = 1.0;
                    }

                    // If the unknown is not pinned
                    if (local_unknown >= 0)
                    {
                      switch (i2)
                      {
                          // radial component
                        case 0:
                          jacobian(local_eqn, local_unknown) +=
                            (2.0 * psif(l2) + r * dpsifdx(l2, 0)) * r *
                            sin_theta * testp(l) * W * hang_weight *
                            hang_weight2;


                          break;

                          // axial component
                        case 1:
                          jacobian(local_eqn, local_unknown) +=
                            r *
                            (dpsifdx(l2, 1) * sin_theta +
                             psif(l2) * cos_theta) *
                            testp(l) * W * hang_weight * hang_weight2;


                          break;

                          // azimuthal component
                        case 2:
                          break;
                      }
                    }
                  }
                }
              } // End of loop over nodes

              // NO PRESSURE CONTRIBUTIONS TO CONTINUITY EQUATION
            } // End of Jacobian calculation
          }
        }
      } // End of loop over pressure nodes

    } // end of loop over integration points
  }

} // namespace oomph
