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
#include "refineable_axisym_navier_stokes_elements.h"


namespace oomph
{
  //========================================================================
  /// Add element's contribution to the elemental
  /// residual vector and/or Jacobian matrix.
  /// flag=1: compute both
  /// flag=0: compute only residual vector
  //=======================================================================
  void RefineableAxisymmetricNavierStokesEquations::
    fill_in_generic_residual_contribution_axi_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // The dimension is actually two
    unsigned DIM = 2;

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Find out how many pressure dofs there are
    unsigned n_pres = npres_axi_nst();

    // Get the local indices of the nodal coordinates
    unsigned u_nodal_index[3];
    for (unsigned i = 0; i < 3; ++i)
    {
      u_nodal_index[i] = u_index_axi_nst(i);
    }

    // Which nodal value represents the pressure? (Negative if pressure
    // is not based on nodal interpolation).
    int p_index = this->p_nodal_index_axi_nst();

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
    DShape dpsifdx(n_node, DIM), dtestfdx(n_node, DIM);


    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Set the value of Nintpt
    unsigned Nintpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get Physical Variables from Element
    // Reynolds number must be multiplied by the density ratio
    double scaled_re = re() * density_ratio();
    double scaled_re_st = re_st() * density_ratio();
    double scaled_re_inv_fr = re_invfr() * density_ratio();
    double scaled_re_inv_ro = re_invro() * density_ratio();
    double visc_ratio = viscosity_ratio(); // hierher -- rewrite and
                                           // make consistent with
                                           // non-refineable version
    Vector<double> G = g();

    // Integers to store the local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0;

    // Local storage for pointers to hang info objects
    HangInfo *hang_info_pt = 0, *hang_info2_pt = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < Nintpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_axi_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Call the pressure shape and test functions
      pshape_axi_nst(s, psip, testp);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the pressure and velocity components
      //--------------------------------------------------------------
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(DIM + 1, 0.0);
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> mesh_veloc(DIM, 0.0);
      Vector<double> dudt(DIM + 1, 0.0);
      DenseMatrix<double> interpolated_dudx(DIM + 1, DIM, 0.0);

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
      {
        interpolated_p += p_axi_nst(l) * psip[l];
      }


      // Calculate velocities and derivatives

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the shape function
        const double psif_ = psif(l);
        // Loop over directions
        for (unsigned i = 0; i < DIM; i++)
        {
          interpolated_x[i] += nodal_position(l, i) * psif_;
          // mesh_veloc[i] +=dnodal_position_dt(l,i)*psif(l);
        }

        for (unsigned i = 0; i < DIM + 1; i++)
        {
          const double u_value = nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif_;
          dudt[i] += du_dt_axi_nst(l, i) * psif_;
          // Loop over derivative directions for gradients
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      // Get the mesh velocity if ALE is enabled
      if (!ALE_is_disabled)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the two coordinate directions
          for (unsigned i = 0; i < 2; i++)
          {
            mesh_veloc[i] += this->dnodal_position_dt(l, i) * psif(l);
          }
        }
      }


      // Get the user-defined body force terms
      Vector<double> body_force(DIM + 1);
      get_body_force_axi_nst(time, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      double source = get_source_fct(time, ipt, interpolated_x);

      // r is the first postition component
      double r = interpolated_x[0];

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
        bool is_node_hanging = node_pt(l)->is_hanging();

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
          for (unsigned i = 0; i < DIM + 1; i++)
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
                  // Add the user-defined body force terms
                  sum += r * body_force[0] * testf[l] * W * hang_weight;

                  // Add the gravitational body force term
                  sum +=
                    r * scaled_re_inv_fr * testf[l] * G[0] * W * hang_weight;

                  // Add the pressure gradient term
                  sum += interpolated_p * (testf[l] + r * dtestfdx(l, 0)) * W *
                         hang_weight;

                  // Add in the stress tensor terms
                  // The viscosity ratio needs to go in here to ensure
                  // continuity of normal stress is satisfied even in flows
                  // with zero pressure gradient!
                  sum -= visc_ratio * r * (1.0 + Gamma[0]) *
                         interpolated_dudx(0, 0) * dtestfdx(l, 0) * W *
                         hang_weight;

                  sum -= visc_ratio * r *
                         (interpolated_dudx(0, 1) +
                          Gamma[0] * interpolated_dudx(1, 0)) *
                         dtestfdx(l, 1) * W * hang_weight;

                  sum -= visc_ratio * (1.0 + Gamma[0]) * interpolated_u[0] *
                         testf[l] * W * hang_weight / r;

                  // Add in the inertial terms
                  // du/dt term
                  sum -=
                    scaled_re_st * r * dudt[0] * testf[l] * W * hang_weight;

                  // Convective terms
                  sum -= scaled_re *
                         (r * interpolated_u[0] * interpolated_dudx(0, 0) -
                          interpolated_u[2] * interpolated_u[2] +
                          r * interpolated_u[1] * interpolated_dudx(0, 1)) *
                         testf[l] * W * hang_weight;

                  // Mesh velocity terms
                  if (!ALE_is_disabled)
                  {
                    for (unsigned k = 0; k < 2; k++)
                    {
                      sum += scaled_re_st * r * mesh_veloc[k] *
                             interpolated_dudx(0, k) * testf[l] * W *
                             hang_weight;
                    }
                  }

                  // Coriolis term
                  sum += 2.0 * r * scaled_re_inv_ro * interpolated_u[2] *
                         testf[l] * W * hang_weight;

                  break;

                  // AXIAL MOMENTUM EQUATION
                case 1:
                  // If it's not a boundary condition
                  // Add the user-defined body force terms
                  // Remember to multiply by the density ratio!
                  sum += r * body_force[1] * testf[l] * W * hang_weight;

                  // Add the gravitational body force term
                  sum +=
                    r * scaled_re_inv_fr * testf[l] * G[1] * W * hang_weight;

                  // Add the pressure gradient term
                  sum += r * interpolated_p * dtestfdx(l, 1) * W * hang_weight;

                  // Add in the stress tensor terms
                  // The viscosity ratio needs to go in here to ensure
                  // continuity of normal stress is satisfied even in flows
                  // with zero pressure gradient!
                  sum -= visc_ratio * r *
                         (interpolated_dudx(1, 0) +
                          Gamma[1] * interpolated_dudx(0, 1)) *
                         dtestfdx(l, 0) * W * hang_weight;

                  sum -= visc_ratio * r * (1.0 + Gamma[1]) *
                         interpolated_dudx(1, 1) * dtestfdx(l, 1) * W *
                         hang_weight;

                  // Add in the inertial terms
                  // du/dt term
                  sum -=
                    scaled_re_st * r * dudt[1] * testf[l] * W * hang_weight;

                  // Convective terms
                  sum -= scaled_re *
                         (r * interpolated_u[0] * interpolated_dudx(1, 0) +
                          r * interpolated_u[1] * interpolated_dudx(1, 1)) *
                         testf[l] * W * hang_weight;

                  // Mesh velocity terms
                  if (!ALE_is_disabled)
                  {
                    for (unsigned k = 0; k < 2; k++)
                    {
                      sum += scaled_re_st * r * mesh_veloc[k] *
                             interpolated_dudx(1, k) * testf[l] * W *
                             hang_weight;
                    }
                  }
                  break;

                  // AZIMUTHAL MOMENTUM EQUATION
                case 2:
                  // Add the user-defined body force terms
                  // Remember to multiply by the density ratio!
                  sum += r * body_force[2] * testf[l] * W * hang_weight;

                  // Add the gravitational body force term
                  sum +=
                    r * scaled_re_inv_fr * testf[l] * G[2] * W * hang_weight;

                  // There is NO pressure gradient term

                  // Add in the stress tensor terms
                  // The viscosity ratio needs to go in here to ensure
                  // continuity of normal stress is satisfied even in flows
                  // with zero pressure gradient!
                  sum -= visc_ratio *
                         (r * interpolated_dudx(2, 0) -
                          Gamma[0] * interpolated_u[2]) *
                         dtestfdx(l, 0) * W * hang_weight;

                  sum -= visc_ratio * r * interpolated_dudx(2, 1) *
                         dtestfdx(l, 1) * W * hang_weight;

                  sum -= visc_ratio *
                         ((interpolated_u[2] / r) -
                          Gamma[0] * interpolated_dudx(2, 0)) *
                         testf[l] * W * hang_weight;


                  // Add in the inertial terms
                  // du/dt term
                  sum -=
                    scaled_re_st * r * dudt[2] * testf[l] * W * hang_weight;

                  // Convective terms
                  sum -= scaled_re *
                         (r * interpolated_u[0] * interpolated_dudx(2, 0) +
                          interpolated_u[0] * interpolated_u[2] +
                          r * interpolated_u[1] * interpolated_dudx(2, 1)) *
                         testf[l] * W * hang_weight;

                  // Mesh velocity terms
                  if (!ALE_is_disabled)
                  {
                    for (unsigned k = 0; k < 2; k++)
                    {
                      sum += scaled_re_st * r * mesh_veloc[k] *
                             interpolated_dudx(2, k) * testf[l] * W *
                             hang_weight;
                    }
                  }

                  // Coriolis term
                  sum -= 2.0 * r * scaled_re_inv_ro * interpolated_u[0] *
                         testf[l] * W * hang_weight;

                  break;
              }

              // Add contribution
              residuals[local_eqn] += sum;

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
                    for (unsigned i2 = 0; i2 < DIM + 1; i2++)
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
                                    scaled_re_st * r * psif[l2] * testf[l] * W *
                                    hang_weight * hang_weight2;
                                }

                                // Add contribution to the Jacobian matrix
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * (1.0 + Gamma[0]) *
                                  dpsifdx(l2, 0) * dtestfdx(l, 0) * W *
                                  hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * dpsifdx(l2, 1) *
                                  dtestfdx(l, 1) * W * hang_weight *
                                  hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * (1.0 + Gamma[0]) * psif[l2] *
                                  testf[l] * W * hang_weight * hang_weight2 / r;

                                // Add in the inertial terms
                                // du/dt term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re_st * r *
                                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                                  psif[l2] * testf[l] * W * hang_weight *
                                  hang_weight2;

                                // Convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re *
                                  (r * psif[l2] * interpolated_dudx(0, 0) +
                                   r * interpolated_u[0] * dpsifdx(l2, 0) +
                                   r * interpolated_u[1] * dpsifdx(l2, 1)) *
                                  testf[l] * W * hang_weight * hang_weight2;

                                // Mesh velocity terms
                                if (!ALE_is_disabled)
                                {
                                  for (unsigned k = 0; k < 2; k++)
                                  {
                                    jacobian(local_eqn, local_unknown) +=
                                      scaled_re_st * r * mesh_veloc[k] *
                                      dpsifdx(l2, k) * testf[l] * W *
                                      hang_weight * hang_weight2;
                                  }
                                }
                                break;

                                // axial component
                              case 1:
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * Gamma[0] * dpsifdx(l2, 0) *
                                  dtestfdx(l, 1) * W * hang_weight *
                                  hang_weight2;

                                // Convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r * psif[l2] *
                                  interpolated_dudx(0, 1) * testf[l] * W *
                                  hang_weight * hang_weight2;
                                break;

                                // azimuthal component
                              case 2:
                                // Convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  -scaled_re * 2.0 * interpolated_u[2] *
                                  psif[l2] * testf[l] * W * hang_weight *
                                  hang_weight2;

                                // Coriolis terms
                                jacobian(local_eqn, local_unknown) +=
                                  2.0 * r * scaled_re_inv_ro * psif[l2] *
                                  testf[l] * W * hang_weight * hang_weight2;

                                break;
                            } /*End of contribution radial momentum eqn*/
                            break;

                            // AXIAL MOMENTUM EQUATION
                          case 1:
                            switch (i2)
                            {
                                // radial component
                              case 0:
                                // Add in the stress tensor terms
                                // The viscosity ratio needs to go in here to
                                // ensure continuity of normal stress is
                                // satisfied even in flows with zero pressure
                                // gradient!
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * Gamma[1] * dpsifdx(l2, 1) *
                                  dtestfdx(l, 0) * W * hang_weight *
                                  hang_weight2;

                                // Convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r * psif[l2] *
                                  interpolated_dudx(1, 0) * testf[l] * W *
                                  hang_weight * hang_weight2;
                                break;

                                // axial component
                              case 1:

                                // Add the mass matrix terms
                                if (flag == 2)
                                {
                                  mass_matrix(local_eqn, local_unknown) +=
                                    scaled_re_st * r * psif[l2] * testf[l] * W *
                                    hang_weight * hang_weight2;
                                }

                                // Add in the stress tensor terms
                                // The viscosity ratio needs to go in here to
                                // ensure continuity of normal stress is
                                // satisfied even in flows with zero pressure
                                // gradient!
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * dpsifdx(l2, 0) *
                                  dtestfdx(l, 0) * W * hang_weight *
                                  hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * (1.0 + Gamma[1]) *
                                  dpsifdx(l2, 1) * dtestfdx(l, 1) * W *
                                  hang_weight * hang_weight2;

                                // Add in the inertial terms
                                // du/dt term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re_st * r *
                                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                                  psif[l2] * testf[l] * W * hang_weight *
                                  hang_weight2;

                                // Convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re *
                                  (r * interpolated_u[0] * dpsifdx(l2, 0) +
                                   r * psif[l2] * interpolated_dudx(1, 1) +
                                   r * interpolated_u[1] * dpsifdx(l2, 1)) *
                                  testf[l] * W * hang_weight * hang_weight2;

                                // Mesh velocity terms
                                if (!ALE_is_disabled)
                                {
                                  for (unsigned k = 0; k < 2; k++)
                                  {
                                    jacobian(local_eqn, local_unknown) +=
                                      scaled_re_st * r * mesh_veloc[k] *
                                      dpsifdx(l2, k) * testf[l] * W *
                                      hang_weight * hang_weight2;
                                  }
                                }
                                break;

                                // azimuthal component
                              case 2:
                                // There are no azimithal terms in the axial
                                // momentum equation
                                break;
                            }
                            break;

                            // AZIMUTHAL MOMENTUM EQUATION
                          case 2:
                            switch (i2)
                            {
                                // radial component
                              case 0:
                                // Convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re *
                                  (r * psif[l2] * interpolated_dudx(2, 0) +
                                   psif[l2] * interpolated_u[2]) *
                                  testf[l] * W * hang_weight * hang_weight2;

                                // Coriolis term
                                jacobian(local_eqn, local_unknown) -=
                                  2.0 * r * scaled_re_inv_ro * psif[l2] *
                                  testf[l] * W * hang_weight * hang_weight2;

                                break;

                                // axial component
                              case 1:
                                // Convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r * psif[l2] *
                                  interpolated_dudx(2, 1) * testf[l] * W *
                                  hang_weight * hang_weight2;
                                break;

                                // azimuthal component
                              case 2:

                                // Add the mass matrix terms
                                if (flag == 2)
                                {
                                  mass_matrix(local_eqn, local_unknown) +=
                                    scaled_re_st * r * psif[l2] * testf[l] * W *
                                    hang_weight * hang_weight2;
                                }

                                // Add in the stress tensor terms
                                // The viscosity ratio needs to go in here to
                                // ensure continuity of normal stress is
                                // satisfied even in flows with zero pressure
                                // gradient!
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio *
                                  (r * dpsifdx(l2, 0) - Gamma[0] * psif[l2]) *
                                  dtestfdx(l, 0) * W * hang_weight *
                                  hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * dpsifdx(l2, 1) *
                                  dtestfdx(l, 1) * W * hang_weight *
                                  hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio *
                                  ((psif[l2] / r) - Gamma[0] * dpsifdx(l2, 0)) *
                                  testf[l] * W * hang_weight * hang_weight2;

                                // Add in the inertial terms
                                // du/dt term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re_st * r *
                                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                                  psif[l2] * testf[l] * W * hang_weight *
                                  hang_weight2;

                                // Convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re *
                                  (r * interpolated_u[0] * dpsifdx(l2, 0) +
                                   interpolated_u[0] * psif[l2] +
                                   r * interpolated_u[1] * dpsifdx(l2, 1)) *
                                  testf[l] * W * hang_weight * hang_weight2;

                                // Mesh velocity terms
                                if (!ALE_is_disabled)
                                {
                                  for (unsigned k = 0; k < 2; k++)
                                  {
                                    jacobian(local_eqn, local_unknown) +=
                                      scaled_re_st * r * mesh_veloc[k] *
                                      dpsifdx(l2, k) * testf[l] * W *
                                      hang_weight * hang_weight2;
                                  }
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
                            psip[l2] * (testf[l] + r * dtestfdx(l, 0)) * W *
                            hang_weight * hang_weight2;
                          break;

                          // AXIAL MOMENTUM EQUATION
                        case 1:
                          jacobian(local_eqn, local_unknown) +=
                            r * psip[l2] * dtestfdx(l, 1) * W * hang_weight *
                            hang_weight2;
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
            // Source term
            residuals[local_eqn] -= r * source * testp[l] * W * hang_weight;

            // Gradient terms
            residuals[local_eqn] +=
              (interpolated_u[0] + r * interpolated_dudx(0, 0) +
               r * interpolated_dudx(1, 1)) *
              testp[l] * W * hang_weight;

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
                  for (unsigned i2 = 0; i2 < DIM + 1; i2++)
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
                            (psif[l2] + r * dpsifdx(l2, 0)) * testp[l] * W *
                            hang_weight * hang_weight2;
                          break;

                          // axial component
                        case 1:
                          jacobian(local_eqn, local_unknown) +=
                            r * dpsifdx(l2, 1) * testp[l] * W * hang_weight *
                            hang_weight2;
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

    } // End of loop over integration points


  } // End of fill_in_generic_residual_contribution_axi_nst(...)


  //======================================================================
  /// Compute derivatives of elemental residual vector with respect
  /// to nodal coordinates.
  /// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
  /// Overloads the FD-based version in the FiniteElement base class.
  //======================================================================
  void RefineableAxisymmetricNavierStokesEquations::
    get_dresidual_dnodal_coordinates(
      RankThreeTensor<double>& dresidual_dnodal_coordinates)
  {
    // Create an Oomph Lib warning
    std::string warning_message = "This function has not been tested.\n";
    std::string function = "RefineableAxisymmetricNavierStokesEquations::\n";
    function += "get_dresidual_dnodal_coordinates(...)";
    OomphLibWarning(warning_message, function, OOMPH_EXCEPTION_LOCATION);

    // Return immediately if there are no dofs
    if (ndof() == 0)
    {
      return;
    }

    // Determine number of nodes in element
    const unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Determine number of pressure dofs in element
    const unsigned n_pres = this->npres_axi_nst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[3];
    for (unsigned i = 0; i < 3; i++)
    {
      u_nodal_index[i] = this->u_index_axi_nst(i);
    }

    // Which nodal value represents the pressure? (Negative if pressure
    // is not based on nodal interpolation).
    const int p_index = this->p_nodal_index_axi_nst();

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
    // Note that there are only two dimensions, r and z, in this problem
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, 2), dtestfdx(n_node, 2);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Determine number of shape controlling nodes
    const unsigned n_shape_controlling_node = nshape_controlling_nodes();

    // Deriatives of shape fct derivatives w.r.t. nodal coords
    RankFourTensor<double> d_dpsifdx_dX(2, n_shape_controlling_node, n_node, 2);
    RankFourTensor<double> d_dtestfdx_dX(
      2, n_shape_controlling_node, n_node, 2);

    // Derivative of Jacobian of mapping w.r.t. to nodal coords
    DenseMatrix<double> dJ_dX(2, n_shape_controlling_node);

    // Derivatives of derivative of u w.r.t. nodal coords
    // Note that the entry d_dudx_dX(p,q,i,k) contains the derivative w.r.t.
    // nodal coordinate X_pq of du_i/dx_k. Since there are three velocity
    // components, i can take on the values 0, 1 and 2, while k and p only
    // take on the values 0 and 1 since there are only two spatial dimensions.
    RankFourTensor<double> d_dudx_dX(2, n_shape_controlling_node, 3, 2);

    // Derivatives of nodal velocities w.r.t. nodal coords:
    // Assumption: Interaction only local via no-slip so
    // X_pq only affects U_iq.
    // Note that the entry d_U_dX(p,q,i) contains the derivative w.r.t. nodal
    // coordinate X_pq of nodal value U_iq. In other words this entry is the
    // derivative of the i-th velocity component w.r.t. the p-th spatial
    // coordinate, taken at the q-th local node.
    RankThreeTensor<double> d_U_dX(2, n_shape_controlling_node, 3, 0.0);

    // Determine the number of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Vector to hold local coordinates (two dimensions)
    Vector<double> s(2);

    // Get physical variables from element
    // (Reynolds number must be multiplied by the density ratio)
    const double scaled_re = this->re() * this->density_ratio();
    const double scaled_re_st = this->re_st() * this->density_ratio();
    const double scaled_re_inv_fr = this->re_invfr() * this->density_ratio();
    const double scaled_re_inv_ro = this->re_invro() * this->density_ratio();
    const double visc_ratio = this->viscosity_ratio();
    const Vector<double> G = this->g();

    // FD step
    double eps_fd = GeneralisedElement::Default_fd_jacobian_step;

    // Pre-compute derivatives of nodal velocities w.r.t. nodal coords:
    // Assumption: Interaction only local via no-slip so
    // X_ij only affects U_ij.
    bool element_has_node_with_aux_node_update_fct = false;

    std::map<Node*, unsigned> local_shape_controlling_node_lookup =
      shape_controlling_node_lookup();

    // FD loop over shape-controlling nodes
    for (std::map<Node*, unsigned>::iterator it =
           local_shape_controlling_node_lookup.begin();
         it != local_shape_controlling_node_lookup.end();
         it++)
    {
      // Get pointer to q-th local shape-controlling node
      Node* nod_pt = it->first;

      // Get its number
      unsigned q = it->second;

      // Only compute if there's a node-update fct involved
      if (nod_pt->has_auxiliary_node_update_fct_pt())
      {
        element_has_node_with_aux_node_update_fct = true;

        // This functionality has not been tested so issue a warning
        // to this effect
        std::ostringstream warning_stream;
        warning_stream << "\nThe functionality to evaluate the additional"
                       << "\ncontribution to the deriv of the residual eqn"
                       << "\nw.r.t. the nodal coordinates which comes about"
                       << "\nif a node's values are updated using an auxiliary"
                       << "\nnode update function has NOT been tested for"
                       << "\nrefineable axisymmetric Navier-Stokes elements."
                       << "\nUse at your own risk" << std::endl;
        OomphLibWarning(warning_stream.str(),
                        "RefineableAxisymmetricNavierStokesEquations::get_"
                        "dresidual_dnodal_coordinates",
                        OOMPH_EXCEPTION_LOCATION);

        // Current nodal velocity
        Vector<double> u_ref(3);
        for (unsigned i = 0; i < 3; i++)
        {
          u_ref[i] = *(nod_pt->value_pt(u_nodal_index[i]));
        }

        // FD
        for (unsigned p = 0; p < 2; p++)
        {
          // Make backup
          double backup = nod_pt->x(p);

          // Do FD step. No node update required as we're
          // attacking the coordinate directly...
          nod_pt->x(p) += eps_fd;

          // Do auxiliary node update (to apply no slip)
          nod_pt->perform_auxiliary_node_update_fct();

          // Loop over velocity components
          for (unsigned i = 0; i < 3; i++)
          {
            // Evaluate
            d_U_dX(p, q, i) =
              (*(nod_pt->value_pt(u_nodal_index[p])) - u_ref[p]) / eps_fd;
          }

          // Reset
          nod_pt->x(p) = backup;

          // Do auxiliary node update (to apply no slip)
          nod_pt->perform_auxiliary_node_update_fct();
        }
      }
    }

    // Integer to store the local equation number
    int local_eqn = 0;

    // Pointers to hang info object
    HangInfo* hang_info_pt = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      const double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      const double J =
        this->dshape_and_dtest_eulerian_at_knot_axi_nst(ipt,
                                                        psif,
                                                        dpsifdx,
                                                        d_dpsifdx_dX,
                                                        testf,
                                                        dtestfdx,
                                                        d_dtestfdx_dX,
                                                        dJ_dX);

      // Call the pressure shape and test functions
      this->pshape_axi_nst(s, psip, testp);

      // Allocate storage for the position and the derivative of the
      // mesh positions w.r.t. time
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> mesh_velocity(2, 0.0);

      // Allocate storage for the pressure, velocity components and their
      // time and spatial derivatives
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(3, 0.0);
      Vector<double> dudt(3, 0.0);
      DenseMatrix<double> interpolated_dudx(3, 2, 0.0);

      // Calculate pressure at integration point
      for (unsigned l = 0; l < n_pres; l++)
      {
        interpolated_p += this->p_axi_nst(l) * psip[l];
      }

      // Calculate velocities and derivatives at integration point
      // ---------------------------------------------------------

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the shape function
        const double psif_ = psif(l);

        // Loop over the two coordinate directions
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += nodal_position(l, i) * psif_;
        }

        // Loop over the three velocity directions
        for (unsigned i = 0; i < 3; i++)
        {
          // Get the nodal value
          const double u_value = nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif_;
          dudt[i] += this->du_dt_axi_nst(l, i) * psif_;

          // Loop over derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      // Get the mesh velocity if ALE is enabled
      if (!this->ALE_is_disabled)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the two coordinate directions
          for (unsigned i = 0; i < 2; i++)
          {
            mesh_velocity[i] += this->dnodal_position_dt(l, i) * psif[l];
          }
        }
      }

      // Calculate derivative of du_i/dx_k w.r.t. nodal positions X_{pq}

      // Loop over shape-controlling nodes
      for (unsigned q = 0; q < n_shape_controlling_node; q++)
      {
        // Loop over the two coordinate directions
        for (unsigned p = 0; p < 2; p++)
        {
          // Loop over the three velocity components
          for (unsigned i = 0; i < 3; i++)
          {
            // Loop over the two coordinate directions
            for (unsigned k = 0; k < 2; k++)
            {
              double aux = 0.0;

              // Loop over nodes and add contribution
              for (unsigned j = 0; j < n_node; j++)
              {
                aux +=
                  nodal_value(j, u_nodal_index[i]) * d_dpsifdx_dX(p, q, j, k);
              }
              d_dudx_dX(p, q, i, k) = aux;
            }
          }
        }
      }

      // Get weight of actual nodal position/value in computation of mesh
      // velocity from positional/value time stepper
      const double pos_time_weight =
        node_pt(0)->position_time_stepper_pt()->weight(1, 0);
      const double val_time_weight =
        node_pt(0)->time_stepper_pt()->weight(1, 0);

      // Get the user-defined body force terms
      Vector<double> body_force(3);
      this->get_body_force_axi_nst(time, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      const double source = this->get_source_fct(time, ipt, interpolated_x);

      // Get gradient of body force function
      DenseMatrix<double> d_body_force_dx(3, 2, 0.0);
      this->get_body_force_gradient_axi_nst(
        time, ipt, s, interpolated_x, d_body_force_dx);

      // Get gradient of source function
      Vector<double> source_gradient(2, 0.0);
      this->get_source_fct_gradient(time, ipt, interpolated_x, source_gradient);

      // Cache r (the first position component)
      const double r = interpolated_x[0];

      // Assemble shape derivatives
      // --------------------------

      // ==================
      // MOMENTUM EQUATIONS
      // ==================

      // Number of master nodes and storage for the weight of the shape function
      unsigned n_master = 1;
      double hang_weight = 1.0;

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the test function
        const double testf_ = testf[l];

        // Local boolean to indicate whether the node is hanging
        bool is_node_hanging = node_pt(l)->is_hanging();

        // If the node is hanging
        if (is_node_hanging)
        {
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
          // --------------------------------
          // FIRST (RADIAL) MOMENTUM EQUATION
          // --------------------------------

          // Get the equation number
          // If the node is hanging
          if (is_node_hanging)
          {
            // Get the equation number from the master node
            local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                             u_nodal_index[0]);
            // Get the hang weight from the master node
            hang_weight = hang_info_pt->master_weight(m);
          }
          // If the node is not hanging
          else
          {
            // Local equation number
            local_eqn = this->nodal_local_eqn(l, u_nodal_index[0]);

            // Node contributes with full weight
            hang_weight = 1.0;
          }

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            // Loop over the two coordinate directions
            for (unsigned p = 0; p < 2; p++)
            {
              // Loop over shape controlling nodes
              for (unsigned q = 0; q < n_shape_controlling_node; q++)
              {
                // Residual x deriv of Jacobian
                // ----------------------------

                // Add the user-defined body force terms
                double sum = r * body_force[0] * testf_;

                // Add the gravitational body force term
                sum += r * scaled_re_inv_fr * testf_ * G[0];

                // Add the pressure gradient term
                sum += interpolated_p * (testf_ + r * dtestfdx(l, 0));

                // Add the stress tensor terms
                // The viscosity ratio needs to go in here to ensure
                // continuity of normal stress is satisfied even in flows
                // with zero pressure gradient!
                sum -= visc_ratio * r * (1.0 + Gamma[0]) *
                       interpolated_dudx(0, 0) * dtestfdx(l, 0);

                sum -= visc_ratio * r *
                       (interpolated_dudx(0, 1) +
                        Gamma[0] * interpolated_dudx(1, 0)) *
                       dtestfdx(l, 1);

                sum -= visc_ratio * (1.0 + Gamma[0]) * interpolated_u[0] *
                       testf_ / r;

                // Add the du/dt term
                sum -= scaled_re_st * r * dudt[0] * testf_;

                // Add the convective terms
                sum -= scaled_re *
                       (r * interpolated_u[0] * interpolated_dudx(0, 0) -
                        interpolated_u[2] * interpolated_u[2] +
                        r * interpolated_u[1] * interpolated_dudx(0, 1)) *
                       testf_;

                // Add the mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned k = 0; k < 2; k++)
                  {
                    sum += scaled_re_st * r * mesh_velocity[k] *
                           interpolated_dudx(0, k) * testf_;
                  }
                }

                // Add the Coriolis term
                sum += 2.0 * r * scaled_re_inv_ro * interpolated_u[2] * testf_;

                // Multiply through by deriv of Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  sum * dJ_dX(p, q) * w * hang_weight;

                // Derivative of residual x Jacobian
                // ---------------------------------

                // Body force
                sum = r * d_body_force_dx(0, p) * psif[q] * testf_;
                if (p == 0)
                {
                  sum += body_force[0] * psif[q] * testf_;
                }

                // Gravitational body force
                if (p == 0)
                {
                  sum += scaled_re_inv_fr * G[0] * psif[q] * testf_;
                }

                // Pressure gradient term
                sum += r * interpolated_p * d_dtestfdx_dX(p, q, l, 0);
                if (p == 0)
                {
                  sum += interpolated_p * psif[q] * dtestfdx(l, 0);
                }

                // Viscous terms
                sum -=
                  r * visc_ratio *
                  ((1.0 + Gamma[0]) *
                     (d_dudx_dX(p, q, 0, 0) * dtestfdx(l, 0) +
                      interpolated_dudx(0, 0) * d_dtestfdx_dX(p, q, l, 0)) +
                   (d_dudx_dX(p, q, 0, 1) + Gamma[0] * d_dudx_dX(p, q, 1, 0)) *
                     dtestfdx(l, 1) +
                   (interpolated_dudx(0, 1) +
                    Gamma[0] * interpolated_dudx(1, 0)) *
                     d_dtestfdx_dX(p, q, l, 1));
                if (p == 0)
                {
                  sum -=
                    visc_ratio *
                    ((1.0 + Gamma[0]) *
                       (interpolated_dudx(0, 0) * psif[q] * dtestfdx(l, 0) -
                        interpolated_u[0] * psif[q] * testf_ / (r * r)) +
                     (interpolated_dudx(0, 1) +
                      Gamma[0] * interpolated_dudx(1, 0)) *
                       psif[q] * dtestfdx(l, 1));
                }

                // Convective terms, including mesh velocity
                for (unsigned k = 0; k < 2; k++)
                {
                  double tmp = scaled_re * interpolated_u[k];
                  if (!ALE_is_disabled)
                  {
                    tmp -= scaled_re_st * mesh_velocity[k];
                  }
                  sum -= r * tmp * d_dudx_dX(p, q, 0, k) * testf_;
                  if (p == 0)
                  {
                    sum -= tmp * interpolated_dudx(0, k) * psif[q] * testf_;
                  }
                }
                if (!ALE_is_disabled)
                {
                  sum += scaled_re_st * r * pos_time_weight *
                         interpolated_dudx(0, p) * psif[q] * testf_;
                }

                // du/dt term
                if (p == 0)
                {
                  sum -= scaled_re_st * dudt[0] * psif[q] * testf_;
                }

                // Coriolis term
                if (p == 0)
                {
                  sum += 2.0 * scaled_re_inv_ro * interpolated_u[2] * psif[q] *
                         testf_;
                }

                // Multiply through by Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  sum * J * w * hang_weight;

              } // End of loop over shape controlling nodes q
            } // End of loop over coordinate directions p

            // Derivs w.r.t. to nodal velocities
            // ---------------------------------
            if (element_has_node_with_aux_node_update_fct)
            {
              // Loop over local nodes
              for (unsigned q_local = 0; q_local < n_node; q_local++)
              {
                // Number of master nodes and storage for the weight of
                // the shape function
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                HangInfo* hang_info2_pt = 0;

                // Local boolean to indicate whether the node is hanging
                bool is_node_hanging2 = node_pt(q_local)->is_hanging();

                Node* actual_shape_controlling_node_pt = node_pt(q_local);

                // If the node is hanging
                if (is_node_hanging2)
                {
                  hang_info2_pt = node_pt(q_local)->hanging_pt();

                  // Read out number of master nodes from hanging data
                  n_master2 = hang_info2_pt->nmaster();
                }
                // Otherwise the node is its own master
                else
                {
                  n_master2 = 1;
                }

                // Loop over the master nodes
                for (unsigned mm = 0; mm < n_master2; mm++)
                {
                  if (is_node_hanging2)
                  {
                    actual_shape_controlling_node_pt =
                      hang_info2_pt->master_node_pt(mm);
                    hang_weight2 = hang_info2_pt->master_weight(mm);
                  }

                  // Find the corresponding number
                  unsigned q = local_shape_controlling_node_lookup
                    [actual_shape_controlling_node_pt];

                  // Loop over the two coordinate directions
                  for (unsigned p = 0; p < 2; p++)
                  {
                    // Contribution from deriv of first velocity component
                    double tmp = scaled_re_st * r * val_time_weight *
                                 psif[q_local] * testf_;
                    tmp += r * scaled_re * interpolated_dudx(0, 0) *
                           psif[q_local] * testf_;

                    for (unsigned k = 0; k < 2; k++)
                    {
                      double tmp2 = scaled_re * interpolated_u[k];
                      if (!ALE_is_disabled)
                      {
                        tmp2 -= scaled_re_st * mesh_velocity[k];
                      }
                      tmp += r * tmp2 * dpsifdx(q_local, k) * testf_;
                    }

                    tmp += r * visc_ratio * (1.0 + Gamma[0]) *
                           dpsifdx(q_local, 0) * dtestfdx(l, 0);
                    tmp +=
                      r * visc_ratio * dpsifdx(q_local, 1) * dtestfdx(l, 1);
                    tmp += visc_ratio * (1.0 + Gamma[0]) * psif[q_local] *
                           testf_ / r;

                    // Multiply through by dU_0q/dX_pq
                    double sum = -d_U_dX(p, q_local, 0) * tmp;

                    // Contribution from deriv of second velocity component
                    tmp = scaled_re * r * interpolated_dudx(0, 1) *
                          psif[q_local] * testf_;
                    tmp += r * visc_ratio * Gamma[0] * dpsifdx(q_local, 0) *
                           dtestfdx(l, 1);

                    // Multiply through by dU_1q/dX_pq
                    sum -= d_U_dX(p, q, 1) * tmp;

                    // Contribution from deriv of third velocity component
                    tmp =
                      2.0 *
                      (r * scaled_re_inv_ro + scaled_re * interpolated_u[2]) *
                      psif[q_local] * testf_;

                    // Multiply through by dU_2q/dX_pq
                    sum += d_U_dX(p, q, 2) * tmp;

                    // Multiply through by Jacobian and integration weight
                    dresidual_dnodal_coordinates(local_eqn, p, q) +=
                      sum * J * w * hang_weight * hang_weight2;
                  }
                } // End of loop over master nodes
              } // End of loop over local nodes
            } // End of if(element_has_node_with_aux_node_update_fct)
          } // End of if it's not a boundary condition

          // --------------------------------
          // SECOND (AXIAL) MOMENTUM EQUATION
          // --------------------------------

          // Get the equation number
          // If the node is hanging
          if (is_node_hanging)
          {
            // Get the equation number from the master node
            local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                             u_nodal_index[1]);
            // Get the hang weight from the master node
            hang_weight = hang_info_pt->master_weight(m);
          }
          // If the node is not hanging
          else
          {
            // Local equation number
            local_eqn = this->nodal_local_eqn(l, u_nodal_index[1]);

            // Node contributes with full weight
            hang_weight = 1.0;
          }

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            // Loop over the two coordinate directions
            for (unsigned p = 0; p < 2; p++)
            {
              // Loop over shape controlling nodes
              for (unsigned q = 0; q < n_shape_controlling_node; q++)
              {
                // Residual x deriv of Jacobian
                // ----------------------------

                // Add the user-defined body force terms
                double sum = r * body_force[1] * testf_;

                // Add the gravitational body force term
                sum += r * scaled_re_inv_fr * testf_ * G[1];

                // Add the pressure gradient term
                sum += r * interpolated_p * dtestfdx(l, 1);

                // Add the stress tensor terms
                // The viscosity ratio needs to go in here to ensure
                // continuity of normal stress is satisfied even in flows
                // with zero pressure gradient!
                sum -= visc_ratio * r *
                       (interpolated_dudx(1, 0) +
                        Gamma[1] * interpolated_dudx(0, 1)) *
                       dtestfdx(l, 0);

                sum -= visc_ratio * r * (1.0 + Gamma[1]) *
                       interpolated_dudx(1, 1) * dtestfdx(l, 1);

                // Add the du/dt term
                sum -= scaled_re_st * r * dudt[1] * testf_;

                // Add the convective terms
                sum -= scaled_re *
                       (r * interpolated_u[0] * interpolated_dudx(1, 0) +
                        r * interpolated_u[1] * interpolated_dudx(1, 1)) *
                       testf_;

                // Add the mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned k = 0; k < 2; k++)
                  {
                    sum += scaled_re_st * r * mesh_velocity[k] *
                           interpolated_dudx(1, k) * testf_;
                  }
                }

                // Multiply through by deriv of Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  sum * dJ_dX(p, q) * w * hang_weight;

                // Derivative of residual x Jacobian
                // ---------------------------------

                // Body force
                sum = r * d_body_force_dx(1, p) * psif[q] * testf_;
                if (p == 0)
                {
                  sum += body_force[1] * psif[q] * testf_;
                }

                // Gravitational body force
                if (p == 0)
                {
                  sum += scaled_re_inv_fr * G[1] * psif[q] * testf_;
                }

                // Pressure gradient term
                sum += r * interpolated_p * d_dtestfdx_dX(p, q, l, 1);
                if (p == 0)
                {
                  sum += interpolated_p * psif[q] * dtestfdx(l, 1);
                }

                // Viscous terms
                sum -=
                  r * visc_ratio *
                  ((d_dudx_dX(p, q, 1, 0) + Gamma[1] * d_dudx_dX(p, q, 0, 1)) *
                     dtestfdx(l, 0) +
                   (interpolated_dudx(1, 0) +
                    Gamma[1] * interpolated_dudx(0, 1)) *
                     d_dtestfdx_dX(p, q, l, 0) +
                   (1.0 + Gamma[1]) * d_dudx_dX(p, q, 1, 1) * dtestfdx(l, 1) +
                   (1.0 + Gamma[1]) * interpolated_dudx(1, 1) *
                     d_dtestfdx_dX(p, q, l, 1));
                if (p == 0)
                {
                  sum -=
                    visc_ratio * ((interpolated_dudx(1, 0) +
                                   Gamma[1] * interpolated_dudx(0, 1)) *
                                    psif[q] * dtestfdx(l, 0) +
                                  (1.0 + Gamma[1]) * interpolated_dudx(1, 1) *
                                    psif[q] * dtestfdx(l, 1));
                }

                // Convective terms, including mesh velocity
                for (unsigned k = 0; k < 2; k++)
                {
                  double tmp = scaled_re * interpolated_u[k];
                  if (!ALE_is_disabled)
                  {
                    tmp -= scaled_re_st * mesh_velocity[k];
                  }
                  sum -= r * tmp * d_dudx_dX(p, q, 1, k) * testf_;
                  if (p == 0)
                  {
                    sum -= tmp * interpolated_dudx(1, k) * psif[q] * testf_;
                  }
                }
                if (!ALE_is_disabled)
                {
                  sum += scaled_re_st * r * pos_time_weight *
                         interpolated_dudx(1, p) * psif[q] * testf_;
                }

                // du/dt term
                if (p == 0)
                {
                  sum -= scaled_re_st * dudt[1] * psif[q] * testf_;
                }

                // Multiply through by Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  sum * J * w * hang_weight;

              } // End of loop over shape controlling nodes q
            } // End of loop over coordinate directions p

            // Derivs w.r.t. to nodal velocities
            // ---------------------------------
            if (element_has_node_with_aux_node_update_fct)
            {
              // Loop over local nodes
              for (unsigned q_local = 0; q_local < n_node; q_local++)
              {
                // Number of master nodes and storage for the weight of
                // the shape function
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                HangInfo* hang_info2_pt = 0;

                // Local boolean to indicate whether the node is hanging
                bool is_node_hanging2 = node_pt(q_local)->is_hanging();

                Node* actual_shape_controlling_node_pt = node_pt(q_local);

                // If the node is hanging
                if (is_node_hanging2)
                {
                  hang_info2_pt = node_pt(q_local)->hanging_pt();

                  // Read out number of master nodes from hanging data
                  n_master2 = hang_info2_pt->nmaster();
                }
                // Otherwise the node is its own master
                else
                {
                  n_master2 = 1;
                }

                // Loop over the master nodes
                for (unsigned mm = 0; mm < n_master2; mm++)
                {
                  if (is_node_hanging2)
                  {
                    actual_shape_controlling_node_pt =
                      hang_info2_pt->master_node_pt(mm);
                    hang_weight2 = hang_info2_pt->master_weight(mm);
                  }

                  // Find the corresponding number
                  unsigned q = local_shape_controlling_node_lookup
                    [actual_shape_controlling_node_pt];

                  // Loop over the two coordinate directions
                  for (unsigned p = 0; p < 2; p++)
                  {
                    // Contribution from deriv of first velocity component
                    double tmp = scaled_re * r * interpolated_dudx(1, 0) *
                                 psif[q_local] * testf_;
                    tmp += r * visc_ratio * Gamma[1] * dpsifdx(q_local, 1) *
                           dtestfdx(l, 0);

                    // Multiply through by dU_0q/dX_pq
                    double sum = -d_U_dX(p, q, 0) * tmp;

                    // Contribution from deriv of second velocity component
                    tmp = scaled_re_st * r * val_time_weight * psif[q_local] *
                          testf_;
                    tmp += r * scaled_re * interpolated_dudx(1, 1) *
                           psif[q_local] * testf_;

                    for (unsigned k = 0; k < 2; k++)
                    {
                      double tmp2 = scaled_re * interpolated_u[k];
                      if (!ALE_is_disabled)
                      {
                        tmp2 -= scaled_re_st * mesh_velocity[k];
                      }
                      tmp += r * tmp2 * dpsifdx(q_local, k) * testf_;
                    }
                    tmp +=
                      r * visc_ratio *
                      (dpsifdx(q_local, 0) * dtestfdx(l, 0) +
                       (1.0 + Gamma[1]) * dpsifdx(q_local, 1) * dtestfdx(l, 1));

                    // Multiply through by dU_1q/dX_pq
                    sum -= d_U_dX(p, q, 1) * tmp;

                    // Multiply through by Jacobian and integration weight
                    dresidual_dnodal_coordinates(local_eqn, p, q) +=
                      sum * J * w * hang_weight * hang_weight2;
                  }
                } // End of loop over master nodes
              } // End of loop over local nodes
            } // End of if(element_has_node_with_aux_node_update_fct)
          } // End of if it's not a boundary condition

          // -----------------------------------
          // THIRD (AZIMUTHAL) MOMENTUM EQUATION
          // -----------------------------------

          // Get the equation number
          // If the node is hanging
          if (is_node_hanging)
          {
            // Get the equation number from the master node
            local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                             u_nodal_index[2]);
            // Get the hang weight from the master node
            hang_weight = hang_info_pt->master_weight(m);
          }
          // If the node is not hanging
          else
          {
            // Local equation number
            local_eqn = this->nodal_local_eqn(l, u_nodal_index[2]);

            // Node contributes with full weight
            hang_weight = 1.0;
          }

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            // Loop over the two coordinate directions
            for (unsigned p = 0; p < 2; p++)
            {
              // Loop over shape controlling nodes
              for (unsigned q = 0; q < n_shape_controlling_node; q++)
              {
                // Residual x deriv of Jacobian
                // ----------------------------

                // Add the user-defined body force terms
                double sum = r * body_force[2] * testf_;

                // Add the gravitational body force term
                sum += r * scaled_re_inv_fr * testf_ * G[2];

                // There is NO pressure gradient term

                // Add in the stress tensor terms
                // The viscosity ratio needs to go in here to ensure
                // continuity of normal stress is satisfied even in flows
                // with zero pressure gradient!
                sum -=
                  visc_ratio *
                  (r * interpolated_dudx(2, 0) - Gamma[0] * interpolated_u[2]) *
                  dtestfdx(l, 0);

                sum -=
                  visc_ratio * r * interpolated_dudx(2, 1) * dtestfdx(l, 1);

                sum -= visc_ratio *
                       ((interpolated_u[2] / r) -
                        Gamma[0] * interpolated_dudx(2, 0)) *
                       testf_;

                // Add the du/dt term
                sum -= scaled_re_st * r * dudt[2] * testf_;

                // Add the convective terms
                sum -= scaled_re *
                       (r * interpolated_u[0] * interpolated_dudx(2, 0) +
                        interpolated_u[0] * interpolated_u[2] +
                        r * interpolated_u[1] * interpolated_dudx(2, 1)) *
                       testf_;

                // Add the mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned k = 0; k < 2; k++)
                  {
                    sum += scaled_re_st * r * mesh_velocity[k] *
                           interpolated_dudx(2, k) * testf_;
                  }
                }

                // Add the Coriolis term
                sum -= 2.0 * r * scaled_re_inv_ro * interpolated_u[0] * testf_;

                // Multiply through by deriv of Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  sum * dJ_dX(p, q) * w * hang_weight;

                // Derivative of residual x Jacobian
                // ---------------------------------

                // Body force
                sum = r * d_body_force_dx(2, p) * psif[q] * testf_;
                if (p == 0)
                {
                  sum += body_force[2] * psif[q] * testf_;
                }

                // Gravitational body force
                if (p == 0)
                {
                  sum += scaled_re_inv_fr * G[2] * psif[q] * testf_;
                }

                // There is NO pressure gradient term

                // Viscous terms
                sum -= visc_ratio * r *
                       (d_dudx_dX(p, q, 2, 0) * dtestfdx(l, 0) +
                        interpolated_dudx(2, 0) * d_dtestfdx_dX(p, q, l, 0) +
                        d_dudx_dX(p, q, 2, 1) * dtestfdx(l, 1) +
                        interpolated_dudx(2, 1) * d_dtestfdx_dX(p, q, l, 1));

                sum += visc_ratio * Gamma[0] * d_dudx_dX(p, q, 2, 0) * testf_;

                if (p == 0)
                {
                  sum -= visc_ratio *
                         (interpolated_dudx(2, 0) * psif[q] * dtestfdx(l, 0) +
                          interpolated_dudx(2, 1) * psif[q] * dtestfdx(l, 1) +
                          interpolated_u[2] * psif[q] * testf_ / (r * r));
                }

                // Convective terms, including mesh velocity
                for (unsigned k = 0; k < 2; k++)
                {
                  double tmp = scaled_re * interpolated_u[k];
                  if (!ALE_is_disabled)
                  {
                    tmp -= scaled_re_st * mesh_velocity[k];
                  }
                  sum -= r * tmp * d_dudx_dX(p, q, 2, k) * testf_;
                  if (p == 0)
                  {
                    sum -= tmp * interpolated_dudx(2, k) * psif[q] * testf_;
                  }
                }
                if (!ALE_is_disabled)
                {
                  sum += scaled_re_st * r * pos_time_weight *
                         interpolated_dudx(2, p) * psif[q] * testf_;
                }

                // du/dt term
                if (p == 0)
                {
                  sum -= scaled_re_st * dudt[2] * psif[q] * testf_;
                }

                // Coriolis term
                if (p == 0)
                {
                  sum -= 2.0 * scaled_re_inv_ro * interpolated_u[0] * psif[q] *
                         testf_;
                }

                // Multiply through by Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  sum * J * w * hang_weight;

              } // End of loop over shape controlling nodes q
            } // End of loop over coordinate directions p

            // Derivs w.r.t. to nodal velocities
            // ---------------------------------
            if (element_has_node_with_aux_node_update_fct)
            {
              // Loop over local nodes
              for (unsigned q_local = 0; q_local < n_node; q_local++)
              {
                // Number of master nodes and storage for the weight of
                // the shape function
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                HangInfo* hang_info2_pt = 0;

                // Local boolean to indicate whether the node is hanging
                bool is_node_hanging2 = node_pt(q_local)->is_hanging();

                Node* actual_shape_controlling_node_pt = node_pt(q_local);

                // If the node is hanging
                if (is_node_hanging2)
                {
                  hang_info2_pt = node_pt(q_local)->hanging_pt();

                  // Read out number of master nodes from hanging data
                  n_master2 = hang_info2_pt->nmaster();
                }
                // Otherwise the node is its own master
                else
                {
                  n_master2 = 1;
                }

                // Loop over the master nodes
                for (unsigned mm = 0; mm < n_master2; mm++)
                {
                  if (is_node_hanging2)
                  {
                    actual_shape_controlling_node_pt =
                      hang_info2_pt->master_node_pt(mm);
                    hang_weight2 = hang_info2_pt->master_weight(mm);
                  }

                  // Find the corresponding number
                  unsigned q = local_shape_controlling_node_lookup
                    [actual_shape_controlling_node_pt];

                  // Loop over the two coordinate directions
                  for (unsigned p = 0; p < 2; p++)
                  {
                    // Contribution from deriv of first velocity component
                    double tmp = (2.0 * r * scaled_re_inv_ro +
                                  r * scaled_re * interpolated_dudx(2, 0) -
                                  scaled_re * interpolated_u[2]) *
                                 psif[q_local] * testf_;

                    // Multiply through by dU_0q/dX_pq
                    double sum = -d_U_dX(p, q, 0) * tmp;

                    // Contribution from deriv of second velocity component
                    // Multiply through by dU_1q/dX_pq
                    sum -= d_U_dX(p, q, 1) * scaled_re * r *
                           interpolated_dudx(2, 1) * psif[q_local] * testf_;

                    // Contribution from deriv of third velocity component
                    tmp = scaled_re_st * r * val_time_weight * psif[q_local] *
                          testf_;
                    tmp -=
                      scaled_re * interpolated_u[0] * psif[q_local] * testf_;

                    for (unsigned k = 0; k < 2; k++)
                    {
                      double tmp2 = scaled_re * interpolated_u[k];
                      if (!ALE_is_disabled)
                      {
                        tmp2 -= scaled_re_st * mesh_velocity[k];
                      }
                      tmp += r * tmp2 * dpsifdx(q_local, k) * testf_;
                    }
                    tmp +=
                      visc_ratio *
                      (r * dpsifdx(q_local, 0) - Gamma[0] * psif[q_local]) *
                      dtestfdx(l, 0);
                    tmp +=
                      r * visc_ratio * dpsifdx(q_local, 1) * dtestfdx(l, 1);
                    tmp +=
                      visc_ratio *
                      ((psif[q_local] / r) - Gamma[0] * dpsifdx(q_local, 0)) *
                      testf_;

                    // Multiply through by dU_2q/dX_pq
                    sum -= d_U_dX(p, q, 2) * tmp;

                    // Multiply through by Jacobian and integration weight
                    dresidual_dnodal_coordinates(local_eqn, p, q) +=
                      sum * J * w * hang_weight * hang_weight2;
                  }
                } // End of loop over master nodes
              } // End of loop over local nodes
            } // End of if(element_has_node_with_aux_node_update_fct)
          } // End of if it's not a boundary condition
        } // End of loop over master nodes
      } // End of loop over test functions


      // ===================
      // CONTINUITY EQUATION
      // ===================

      // Loop over the shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        // If the pressure dof is hanging
        if (pressure_dof_is_hanging[l])
        {
          // Pressure dof is hanging so it must be nodal-based
          // Get the hang info object
          hang_info_pt = this->pressure_node_pt(l)->hanging_pt(p_index);

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
            // Get the local equation from the master node
            local_eqn =
              this->local_hang_eqn(hang_info_pt->master_node_pt(m), p_index);
            // Get the weight for the node
            hang_weight = hang_info_pt->master_weight(m);
          }
          else
          {
            local_eqn = this->p_local_eqn(l);
            hang_weight = 1.0;
          }

          // Cache the test function
          const double testp_ = testp[l];

          // If not a boundary conditions
          if (local_eqn >= 0)
          {
            // Loop over the two coordinate directions
            for (unsigned p = 0; p < 2; p++)
            {
              // Loop over shape controlling nodes
              for (unsigned q = 0; q < n_shape_controlling_node; q++)
              {
                // Residual x deriv of Jacobian
                //-----------------------------

                // Source term
                double aux = -r * source;

                // Gradient terms
                aux += (interpolated_u[0] + r * interpolated_dudx(0, 0) +
                        r * interpolated_dudx(1, 1));

                // Multiply through by deriv of Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  aux * dJ_dX(p, q) * testp_ * w * hang_weight;

                // Derivative of residual x Jacobian
                // ---------------------------------

                // Gradient of source function
                aux = -r * source_gradient[p] * psif[q];
                if (p == 0)
                {
                  aux -= source * psif[q];
                }

                // Gradient terms
                aux += r * (d_dudx_dX(p, q, 0, 0) + d_dudx_dX(p, q, 1, 1));
                if (p == 0)
                {
                  aux += (interpolated_dudx(0, 0) + interpolated_dudx(1, 1)) *
                         psif[q];
                }

                // Multiply through by Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  aux * testp_ * J * w * hang_weight;
              }
            }

            // Derivs w.r.t. to nodal velocities
            // ---------------------------------
            if (element_has_node_with_aux_node_update_fct)
            {
              // Loop over local nodes
              for (unsigned q_local = 0; q_local < n_node; q_local++)
              {
                // Number of master nodes and storage for the weight of
                // the shape function
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                HangInfo* hang_info2_pt = 0;

                // Local boolean to indicate whether the node is hanging
                bool is_node_hanging2 = node_pt(q_local)->is_hanging();

                Node* actual_shape_controlling_node_pt = node_pt(q_local);

                // If the node is hanging
                if (is_node_hanging2)
                {
                  hang_info2_pt = node_pt(q_local)->hanging_pt();

                  // Read out number of master nodes from hanging data
                  n_master2 = hang_info2_pt->nmaster();
                }
                // Otherwise the node is its own master
                else
                {
                  n_master2 = 1;
                }

                // Loop over the master nodes
                for (unsigned mm = 0; mm < n_master2; mm++)
                {
                  if (is_node_hanging2)
                  {
                    actual_shape_controlling_node_pt =
                      hang_info2_pt->master_node_pt(mm);
                    hang_weight2 = hang_info2_pt->master_weight(mm);
                  }

                  // Find the corresponding number
                  unsigned q = local_shape_controlling_node_lookup
                    [actual_shape_controlling_node_pt];

                  // Loop over the two coordinate directions
                  for (unsigned p = 0; p < 2; p++)
                  {
                    double aux = d_U_dX(p, q, 0) *
                                   (psif[q_local] + r * dpsifdx(q_local, 0)) +
                                 d_U_dX(p, q, 1) * r * dpsifdx(q_local, 1);

                    // Multiply through by Jacobian, test function and
                    // integration weight
                    dresidual_dnodal_coordinates(local_eqn, p, q) +=
                      aux * testp_ * J * w * hang_weight * hang_weight2;
                  }
                } // End of loop over (mm) master nodes
              } // End of loop over local nodes q_local
            } // End of if(element_has_node_with_aux_node_update_fct)
          } // End of if it's not a boundary condition
        } // End of loop over (m) master nodes
      } // End of loop over shape functions for continuity eqn

    } // End of loop over integration points

  } // End of get_dresidual_dnodal_coordinates(...)


} // namespace oomph
