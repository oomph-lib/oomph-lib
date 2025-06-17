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
// Non-inline functions for the refineable linearised axisymmetric
// Navier-Stokes elements

// oomph-lib includes
#include "refineable_linearised_axisym_navier_stokes_elements.h"

namespace oomph
{
  //=======================================================================
  /// Compute the residuals for the refineable linearised axisymmetric
  /// Navier--Stokes equations; flag=1(or 0): do (or don't) compute the
  /// Jacobian as well.
  //=======================================================================
  void RefineableLinearisedAxisymmetricNavierStokesEquations::
    fill_in_generic_residual_contribution_linearised_axi_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Get the time from the first node in the element
    const double time = this->node_pt(0)->time_stepper_pt()->time();

    // Determine number of nodes in the element
    const unsigned n_node = nnode();

    // Determine how many pressure values there are associated with
    // a single pressure component
    const unsigned n_pres = npres_linearised_axi_nst();

    // Get the nodal indices at which the velocity is stored
    unsigned u_nodal_index[6];
    for (unsigned i = 0; i < 6; ++i)
    {
      u_nodal_index[i] = u_index_linearised_axi_nst(i);
    }

    // Which nodal values represent the two pressure components?
    // (Negative if pressure is not based on nodal interpolation).
    Vector<int> p_index(2);
    for (unsigned i = 0; i < 2; i++)
    {
      p_index[i] = this->p_index_linearised_axi_nst(i);
    }

    // Local array of booleans that are true if the l-th pressure value is
    // hanging (avoid repeated virtual function calls)
    bool pressure_dof_is_hanging[n_pres];

    // If the pressure is stored at a node
    if (p_index[0] >= 0)
    {
      // Read out whether the pressure is hanging
      for (unsigned l = 0; l < n_pres; ++l)
      {
        pressure_dof_is_hanging[l] =
          pressure_node_pt(l)->is_hanging(p_index[0]);
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

    // Set up memory for the fluid shape and test functions
    // Note that there are two dimensions, r and z, in this problem
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, 2), dtestfdx(n_node, 2);

    // Set up memory for the pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Determine number of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Set up memory for the vector to hold local coordinates (two dimensions)
    Vector<double> s(2);

    // Get physical variables from the element
    // (Reynolds number must be multiplied by the density ratio)
    const double scaled_re = re() * density_ratio();
    const double scaled_re_st = re_st() * density_ratio();
    const double visc_ratio = viscosity_ratio();
    const int k = azimuthal_mode_number();

    // Integers used to store the local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0;

    // Local storage for pointers to hang info objects
    HangInfo *hang_info_pt = 0, *hang_info2_pt = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of the local coordinates s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      const double w = integral_pt()->weight(ipt);

      // Calculate the fluid shape and test functions, and their derivatives
      // w.r.t. the global coordinates
      const double J = dshape_and_dtest_eulerian_at_knot_linearised_axi_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Calculate the pressure shape and test functions
      pshape_linearised_axi_nst(s, psip, testp);

      // Premultiply the weights and the Jacobian of the mapping between
      // local and global coordinates
      const double W = w * J;

      // Allocate storage for the position and the derivative of the
      // mesh positions w.r.t. time
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> mesh_velocity(2, 0.0);

      // Allocate storage for the velocity components (six of these)
      // and their derivatives w.r.t. time
      Vector<double> interpolated_u(6, 0.0);
      Vector<double> dudt(6, 0.0);

      // Allocate storage for the pressure components (two of these)
      Vector<double> interpolated_p(2, 0.0);

      // Allocate storage for the derivatives of the velocity components
      // w.r.t. global coordinates (r and z)
      DenseMatrix<double> interpolated_dudx(6, 2, 0.0);

      // Calculate pressure at the integration point
      // -------------------------------------------

      // Loop over pressure degrees of freedom (associated with a single
      // pressure component) in the element
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Cache the shape function
        const double psip_ = psip(l);

        // Loop over the two pressure components
        for (unsigned i = 0; i < 2; i++)
        {
          // Get the value
          const double p_value = this->p_linearised_axi_nst(l, i);

          // Add contribution
          interpolated_p[i] += p_value * psip_;
        }
      } // End of loop over the pressure degrees of freedom in the element

      // Calculate velocities and their derivatives at the integration point
      // -------------------------------------------------------------------

      // Loop over the element's nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the shape function
        const double psif_ = psif(l);

        // Loop over the two coordinate directions
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += this->nodal_position(l, i) * psif_;
        }

        // Loop over the six velocity components
        for (unsigned i = 0; i < 6; i++)
        {
          // Get the value
          const double u_value = this->nodal_value(l, u_nodal_index[i]);

          // Add contribution
          interpolated_u[i] += u_value * psif_;

          // Add contribution to dudt
          dudt[i] += du_dt_linearised_axi_nst(l, i) * psif_;

          // Loop over two coordinate directions (for derivatives)
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      } // End of loop over the element's nodes

      // Get the mesh velocity if ALE is enabled
      if (!ALE_is_disabled)
      {
        // Loop over the element's nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the two coordinate directions
          for (unsigned i = 0; i < 2; i++)
          {
            mesh_velocity[i] += this->dnodal_position_dt(l, i) * psif(l);
          }
        }
      }

      // Get velocities and their derivatives from base flow problem
      // -----------------------------------------------------------

      // Allocate storage for the velocity components of the base state
      // solution (initialise to zero)
      Vector<double> base_flow_u(3, 0.0);

      // Get the user-defined base state solution velocity components
      get_base_flow_u(time, ipt, interpolated_x, base_flow_u);

      // Allocate storage for the derivatives of the base state solution's
      // velocity components w.r.t. global coordinate (r and z)
      // N.B. the derivatives of the base flow components w.r.t. the
      // azimuthal coordinate direction (theta) are always zero since the
      // base flow is axisymmetric
      DenseMatrix<double> base_flow_dudx(3, 2, 0.0);

      // Get the derivatives of the user-defined base state solution
      // velocity components w.r.t. global coordinates
      get_base_flow_dudx(time, ipt, interpolated_x, base_flow_dudx);

      // Cache base flow velocities and their derivatives
      const double interpolated_ur = base_flow_u[0];
      const double interpolated_uz = base_flow_u[1];
      const double interpolated_utheta = base_flow_u[2];
      const double interpolated_durdr = base_flow_dudx(0, 0);
      const double interpolated_durdz = base_flow_dudx(0, 1);
      const double interpolated_duzdr = base_flow_dudx(1, 0);
      const double interpolated_duzdz = base_flow_dudx(1, 1);
      const double interpolated_duthetadr = base_flow_dudx(2, 0);
      const double interpolated_duthetadz = base_flow_dudx(2, 1);

      // Cache r-component of position
      const double r = interpolated_x[0];

      // Cache unknowns
      const double interpolated_UC = interpolated_u[0];
      const double interpolated_US = interpolated_u[1];
      const double interpolated_WC = interpolated_u[2];
      const double interpolated_WS = interpolated_u[3];
      const double interpolated_VC = interpolated_u[4];
      const double interpolated_VS = interpolated_u[5];
      const double interpolated_PC = interpolated_p[0];
      const double interpolated_PS = interpolated_p[1];

      // Cache derivatives of the unknowns
      const double interpolated_dUCdr = interpolated_dudx(0, 0);
      const double interpolated_dUCdz = interpolated_dudx(0, 1);
      const double interpolated_dUSdr = interpolated_dudx(1, 0);
      const double interpolated_dUSdz = interpolated_dudx(1, 1);
      const double interpolated_dWCdr = interpolated_dudx(2, 0);
      const double interpolated_dWCdz = interpolated_dudx(2, 1);
      const double interpolated_dWSdr = interpolated_dudx(3, 0);
      const double interpolated_dWSdz = interpolated_dudx(3, 1);
      const double interpolated_dVCdr = interpolated_dudx(4, 0);
      const double interpolated_dVCdz = interpolated_dudx(4, 1);
      const double interpolated_dVSdr = interpolated_dudx(5, 0);
      const double interpolated_dVSdz = interpolated_dudx(5, 1);

      // ==================
      // MOMENTUM EQUATIONS
      // ==================

      // Number of master nodes
      unsigned n_master = 1;

      // Storage for the weight of the shape function
      double hang_weight = 1.0;

      // Loop over the nodes for the test functions/equations
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache test functions and their derivatives
        const double testf_ = testf(l);
        const double dtestfdr = dtestfdx(l, 0);
        const double dtestfdz = dtestfdx(l, 1);

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
          for (unsigned i = 0; i < 6; i++)
          {
            // Get the equation number
            // -----------------------

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
                  // ---------------------------------------------
                  // FIRST (RADIAL) MOMENTUM EQUATION: COSINE PART
                  // ---------------------------------------------

                case 0:

                  // Pressure gradient term
                  sum +=
                    interpolated_PC * (testf_ + r * dtestfdr) * W * hang_weight;

                  // Stress tensor terms
                  sum -= visc_ratio * r * (1.0 + Gamma[0]) *
                         interpolated_dUCdr * dtestfdr * W * hang_weight;

                  sum -= visc_ratio * r *
                         (interpolated_dUCdz + Gamma[0] * interpolated_dWCdr) *
                         dtestfdz * W * hang_weight;

                  sum += visc_ratio *
                         ((k * Gamma[0] * interpolated_dVSdr) -
                          (k * (2.0 + Gamma[0]) * interpolated_VS / r) -
                          ((1.0 + Gamma[0] + (k * k)) * interpolated_UC / r)) *
                         testf_ * W * hang_weight;

                  // Inertial terms (du/dt)
                  sum -= scaled_re_st * r * dudt[0] * testf_ * W * hang_weight;

                  // Inertial terms (convective)
                  sum -= scaled_re *
                         (r * interpolated_ur * interpolated_dUCdr +
                          r * interpolated_durdr * interpolated_UC +
                          k * interpolated_utheta * interpolated_US -
                          2 * interpolated_utheta * interpolated_VC +
                          r * interpolated_uz * interpolated_dUCdz +
                          r * interpolated_durdz * interpolated_WC) *
                         testf_ * W * hang_weight;

                  // Mesh velocity terms
                  if (!ALE_is_disabled)
                  {
                    for (unsigned j = 0; j < 2; j++)
                    {
                      sum += scaled_re_st * r * mesh_velocity[j] *
                             interpolated_dudx(0, j) * testf_ * W * hang_weight;
                    }
                  }

                  break;

                  // --------------------------------------------
                  // SECOND (RADIAL) MOMENTUM EQUATION: SINE PART
                  // --------------------------------------------

                case 1:

                  // Pressure gradient term
                  sum +=
                    interpolated_PS * (testf_ + r * dtestfdr) * W * hang_weight;

                  // Stress tensor terms
                  sum -= visc_ratio * r * (1.0 + Gamma[0]) *
                         interpolated_dUSdr * dtestfdr * W * hang_weight;

                  sum -= visc_ratio * r *
                         (interpolated_dUSdz + Gamma[0] * interpolated_dWSdr) *
                         dtestfdz * W * hang_weight;

                  sum -= visc_ratio *
                         ((k * Gamma[0] * interpolated_dVCdr) -
                          (k * (2.0 + Gamma[0]) * interpolated_VC / r) +
                          ((1.0 + Gamma[0] + (k * k)) * interpolated_US / r)) *
                         testf_ * W * hang_weight;

                  // Inertial terms (du/dt)
                  sum -= scaled_re_st * r * dudt[1] * testf_ * W * hang_weight;

                  // Inertial terms (convective)
                  sum -= scaled_re *
                         (r * interpolated_ur * interpolated_dUSdr +
                          r * interpolated_durdr * interpolated_US -
                          k * interpolated_utheta * interpolated_UC -
                          2 * interpolated_utheta * interpolated_VS +
                          r * interpolated_uz * interpolated_dUSdz +
                          r * interpolated_durdz * interpolated_WS) *
                         testf_ * W * hang_weight;

                  // Mesh velocity terms
                  if (!ALE_is_disabled)
                  {
                    for (unsigned j = 0; j < 2; j++)
                    {
                      sum += scaled_re_st * r * mesh_velocity[j] *
                             interpolated_dudx(1, j) * testf_ * W * hang_weight;
                    }
                  }

                  break;

                  // --------------------------------------------
                  // THIRD (AXIAL) MOMENTUM EQUATION: COSINE PART
                  // --------------------------------------------

                case 2:

                  // Pressure gradient term
                  sum += r * interpolated_PC * dtestfdz * W * hang_weight;

                  // Stress tensor terms
                  sum -= visc_ratio * r *
                         (interpolated_dWCdr + Gamma[1] * interpolated_dUCdz) *
                         dtestfdr * W * hang_weight;

                  sum -= visc_ratio * r * (1.0 + Gamma[1]) *
                         interpolated_dWCdz * dtestfdz * W * hang_weight;

                  sum +=
                    visc_ratio * k *
                    (Gamma[1] * interpolated_dVSdz - k * interpolated_WC / r) *
                    testf_ * W * hang_weight;

                  // Inertial terms (du/dt)
                  sum -= scaled_re_st * r * dudt[2] * testf_ * W * hang_weight;

                  // Inertial terms (convective)
                  sum -= scaled_re *
                         (r * interpolated_ur * interpolated_dWCdr +
                          r * interpolated_duzdr * interpolated_UC +
                          k * interpolated_utheta * interpolated_WS +
                          r * interpolated_uz * interpolated_dWCdz +
                          r * interpolated_duzdz * interpolated_WC) *
                         testf_ * W * hang_weight;

                  // Mesh velocity terms
                  if (!ALE_is_disabled)
                  {
                    for (unsigned j = 0; j < 2; j++)
                    {
                      sum += scaled_re_st * r * mesh_velocity[j] *
                             interpolated_dudx(2, j) * testf_ * W * hang_weight;
                    }
                  }

                  break;

                  // -------------------------------------------
                  // FOURTH (AXIAL) MOMENTUM EQUATION: SINE PART
                  // -------------------------------------------

                case 3:

                  // Pressure gradient term
                  sum += r * interpolated_PS * dtestfdz * W * hang_weight;

                  // Stress tensor terms
                  sum -= visc_ratio * r *
                         (interpolated_dWSdr + Gamma[1] * interpolated_dUSdz) *
                         dtestfdr * W * hang_weight;

                  sum -= visc_ratio * r * (1.0 + Gamma[1]) *
                         interpolated_dWSdz * dtestfdz * W * hang_weight;

                  sum -=
                    visc_ratio * k *
                    (Gamma[1] * interpolated_dVCdz + k * interpolated_WS / r) *
                    testf_ * W * hang_weight;

                  // Inertial terms (du/dt)
                  sum -= scaled_re_st * r * dudt[3] * testf_ * W * hang_weight;

                  // Inertial terms (convective)
                  sum -= scaled_re *
                         (r * interpolated_ur * interpolated_dWSdr +
                          r * interpolated_duzdr * interpolated_US -
                          k * interpolated_utheta * interpolated_WC +
                          r * interpolated_uz * interpolated_dWSdz +
                          r * interpolated_duzdz * interpolated_WS) *
                         testf_ * W * hang_weight;

                  // Mesh velocity terms
                  if (!ALE_is_disabled)
                  {
                    for (unsigned j = 0; j < 2; j++)
                    {
                      sum += scaled_re_st * r * mesh_velocity[j] *
                             interpolated_dudx(3, j) * testf_ * W * hang_weight;
                    }
                  }

                  break;

                  // ------------------------------------------------
                  // FIFTH (AZIMUTHAL) MOMENTUM EQUATION: COSINE PART
                  // ------------------------------------------------

                case 4:

                  // Pressure gradient term
                  sum -= k * interpolated_PS * testf_ * W * hang_weight;

                  // Stress tensor terms
                  sum +=
                    visc_ratio *
                    (-r * interpolated_dVCdr - k * Gamma[0] * interpolated_US +
                     Gamma[0] * interpolated_VC) *
                    dtestfdr * W * hang_weight;

                  sum -=
                    visc_ratio *
                    (k * Gamma[0] * interpolated_WS + r * interpolated_dVCdz) *
                    dtestfdz * W * hang_weight;

                  sum += visc_ratio *
                         (Gamma[0] * interpolated_dVCdr +
                          k * (2.0 + Gamma[0]) * interpolated_US / r -
                          (1.0 + (k * k) + (k * k * Gamma[0])) *
                            interpolated_VC / r) *
                         testf_ * W * hang_weight;

                  // Inertial terms (du/dt)
                  sum -= scaled_re_st * r * dudt[4] * testf_ * W * hang_weight;

                  // Inertial terms (convective)
                  sum -= scaled_re *
                         (r * interpolated_ur * interpolated_dVCdr +
                          r * interpolated_duthetadr * interpolated_UC +
                          k * interpolated_utheta * interpolated_VS +
                          interpolated_utheta * interpolated_UC +
                          interpolated_ur * interpolated_VC +
                          r * interpolated_uz * interpolated_dVCdz +
                          r * interpolated_duthetadz * interpolated_WC) *
                         testf_ * W * hang_weight;

                  // Mesh velocity terms
                  if (!ALE_is_disabled)
                  {
                    for (unsigned j = 0; j < 2; j++)
                    {
                      sum += scaled_re_st * r * mesh_velocity[j] *
                             interpolated_dudx(4, j) * testf_ * W * hang_weight;
                    }
                  }

                  break;

                  // ----------------------------------------------
                  // SIXTH (AZIMUTHAL) MOMENTUM EQUATION: SINE PART
                  // ----------------------------------------------

                case 5:

                  // Pressure gradient term
                  sum += k * interpolated_PC * testf_ * W * hang_weight;

                  // Stress tensor terms
                  sum +=
                    visc_ratio *
                    (-r * interpolated_dVSdr + k * Gamma[0] * interpolated_UC +
                     Gamma[0] * interpolated_VS) *
                    dtestfdr * W * hang_weight;

                  sum +=
                    visc_ratio *
                    (k * Gamma[0] * interpolated_WC - r * interpolated_dVSdz) *
                    dtestfdz * W * hang_weight;

                  sum += visc_ratio *
                         (Gamma[0] * interpolated_dVSdr -
                          k * (2.0 + Gamma[0]) * interpolated_UC / r -
                          (1.0 + (k * k) + (k * k * Gamma[0])) *
                            interpolated_VS / r) *
                         testf_ * W * hang_weight;

                  // Inertial terms (du/dt)
                  sum -= scaled_re_st * r * dudt[5] * testf_ * W * hang_weight;

                  // Inertial terms (convective)
                  sum -= scaled_re *
                         (r * interpolated_ur * interpolated_dVSdr +
                          r * interpolated_duthetadr * interpolated_US -
                          k * interpolated_utheta * interpolated_VC +
                          interpolated_utheta * interpolated_US +
                          interpolated_ur * interpolated_VS +
                          r * interpolated_uz * interpolated_dVSdz +
                          r * interpolated_duthetadz * interpolated_WS) *
                         testf_ * W * hang_weight;

                  // Mesh velocity terms
                  if (!ALE_is_disabled)
                  {
                    for (unsigned j = 0; j < 2; j++)
                    {
                      sum += scaled_re_st * r * mesh_velocity[j] *
                             interpolated_dudx(5, j) * testf_ * W * hang_weight;
                    }
                  }

                  break;

              } // End of switch statement for momentum equations

              // Add contribution to elemental residual vector
              residuals[local_eqn] += sum;

              // ======================
              // Calculate the Jacobian
              // ======================
              if (flag)
              {
                // Number of master nodes
                unsigned n_master2 = 1;

                // Storage for the weight of the shape function
                double hang_weight2 = 1.0;

                // Loop over the velocity nodes for columns
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Cache velocity shape functions and their derivatives
                  const double psif_ = psif[l2];
                  const double dpsifdr = dpsifdx(l2, 0);
                  const double dpsifdz = dpsifdx(l2, 1);

                  // Local boolean that indicates the hanging status of the node
                  bool is_node2_hanging = node_pt(l2)->is_hanging();

                  // If the node is hanging
                  if (is_node2_hanging)
                  {
                    // Get the hanging pointer
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
                    for (unsigned i2 = 0; i2 < 6; i2++)
                    {
                      // Get the number of the unknown
                      // -----------------------------

                      // If the node is hanging
                      if (is_node2_hanging)
                      {
                        // Get the equation number from the master node
                        local_unknown = this->local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), u_nodal_index[i2]);

                        // Get the hang weight from the master node
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      // If the node is not hanging
                      else
                      {
                        // Local equation number
                        local_unknown =
                          this->nodal_local_eqn(l2, u_nodal_index[i2]);

                        // Node contributes with full weight
                        hang_weight2 = 1.0;
                      }

                      // If the unknown is not pinned
                      if (local_unknown >= 0)
                      {
                        // Different results for i and i2
                        switch (i)
                        {
                            // ---------------------------------------------
                            // FIRST (RADIAL) MOMENTUM EQUATION: COSINE PART
                            // ---------------------------------------------

                          case 0:

                            switch (i2)
                            {
                                // Radial velocity component (cosine part) U_k^C

                              case 0:

                                // Add the mass matrix entries
                                if (flag == 2)
                                {
                                  mass_matrix(local_eqn, local_unknown) +=
                                    scaled_re_st * r * psif_ * testf_ * W *
                                    hang_weight * hang_weight2 * hang_weight *
                                    hang_weight2;
                                }

                                // Add contributions to the Jacobian matrix

                                // Stress tensor terms
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * (1.0 + Gamma[0]) * dpsifdr *
                                  dtestfdr * W * hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * dpsifdz * dtestfdz * W *
                                  hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * (1.0 + Gamma[0] + (k * k)) *
                                  psif_ * testf_ * W * hang_weight *
                                  hang_weight2 / r;

                                // Inertial terms (du/dt)
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re_st * r *
                                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                                  psif_ * testf_ * W * hang_weight *
                                  hang_weight2;

                                // Inertial terms (convective)
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r *
                                  (psif_ * interpolated_durdr +
                                   interpolated_ur * dpsifdr +
                                   interpolated_uz * dpsifdz) *
                                  testf_ * W * hang_weight * hang_weight2;

                                // Mesh velocity terms
                                if (!ALE_is_disabled)
                                {
                                  for (unsigned j = 0; j < 2; j++)
                                  {
                                    jacobian(local_eqn, local_unknown) +=
                                      scaled_re_st * r * mesh_velocity[j] *
                                      dpsifdx(l2, j) * testf_ * W *
                                      hang_weight * hang_weight2;
                                  }
                                }

                                break;

                                // Radial velocity component (sine part) U_k^S

                              case 1:

                                // Convective term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * k * interpolated_utheta * psif_ *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                                // Axial velocity component (cosine part) W_k^C

                              case 2:

                                // Stress tensor term
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * Gamma[0] * r * dpsifdr *
                                  dtestfdz * W * hang_weight * hang_weight2;

                                // Convective term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r * interpolated_durdz * psif_ *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                                // Axial velocity component (sine part) W_k^S
                                // has no contribution

                              case 3:
                                break;

                                // Azimuthal velocity component (cosine part)
                                // V_k^C

                              case 4:

                                // Convective term
                                jacobian(local_eqn, local_unknown) +=
                                  scaled_re * 2 * interpolated_utheta * psif_ *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                                // Azimuthal velocity component (sine part)
                                // V_k^S

                              case 5:

                                // Stress tensor term
                                jacobian(local_eqn, local_unknown) +=
                                  visc_ratio *
                                  ((Gamma[0] * k * dpsifdr) -
                                   (k * (2.0 + Gamma[0]) * psif_ / r)) *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                            } // End of first (radial) momentum equation
                            break;

                            // --------------------------------------------
                            // SECOND (RADIAL) MOMENTUM EQUATION: SINE PART
                            // --------------------------------------------

                          case 1:

                            switch (i2)
                            {
                                // Radial velocity component (cosine part) U_k^C

                              case 0:

                                // Convective term
                                jacobian(local_eqn, local_unknown) +=
                                  scaled_re * k * interpolated_utheta * psif_ *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                                // Radial velocity component (sine part) U_k^S

                              case 1:

                                // Add the mass matrix entries
                                if (flag == 2)
                                {
                                  mass_matrix(local_eqn, local_unknown) +=
                                    scaled_re_st * r * psif_ * testf_ * W *
                                    hang_weight * hang_weight2;
                                }

                                // Add contributions to the Jacobian matrix

                                // Stress tensor terms
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * (1.0 + Gamma[0]) * dpsifdr *
                                  dtestfdr * W * hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * dpsifdz * dtestfdz * W *
                                  hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * (1.0 + Gamma[0] + (k * k)) *
                                  psif_ * testf_ * W * hang_weight *
                                  hang_weight2 / r;

                                // Inertial terms (du/dt)
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re_st * r *
                                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                                  psif_ * testf_ * W * hang_weight *
                                  hang_weight2;

                                // Inertial terms (convective)
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r *
                                  (psif_ * interpolated_durdr +
                                   interpolated_ur * dpsifdr +
                                   interpolated_uz * dpsifdz) *
                                  testf_ * W * hang_weight * hang_weight2;

                                // Mesh velocity terms
                                if (!ALE_is_disabled)
                                {
                                  for (unsigned j = 0; j < 2; j++)
                                  {
                                    jacobian(local_eqn, local_unknown) +=
                                      scaled_re_st * r * mesh_velocity[j] *
                                      dpsifdx(l2, j) * testf_ * W *
                                      hang_weight * hang_weight2;
                                  }
                                }

                                break;

                                // Axial velocity component (cosine part) W_k^C
                                // has no contribution

                              case 2:
                                break;

                                // Axial velocity component (sine part) W_k^S

                              case 3:

                                // Stress tensor term
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * Gamma[0] * r * dpsifdr *
                                  dtestfdz * W * hang_weight * hang_weight2;

                                // Convective term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r * interpolated_durdz * psif_ *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                                // Azimuthal velocity component (cosine part)
                                // V_k^C

                              case 4:

                                // Stress tensor terms
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio *
                                  ((Gamma[0] * k * dpsifdr) -
                                   (k * (2.0 + Gamma[0]) * psif_ / r)) *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                                // Azimuthal velocity component (sine part)
                                // V_k^S

                              case 5:

                                // Convective term
                                jacobian(local_eqn, local_unknown) +=
                                  scaled_re * 2 * interpolated_utheta * psif_ *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                            } // End of second (radial) momentum equation
                            break;

                            // --------------------------------------------
                            // THIRD (AXIAL) MOMENTUM EQUATION: COSINE PART
                            // --------------------------------------------

                          case 2:

                            switch (i2)
                            {
                                // Radial velocity component (cosine part) U_k^C

                              case 0:

                                // Stress tensor term
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * Gamma[1] * dpsifdz *
                                  dtestfdr * W * hang_weight * hang_weight2;

                                // Convective term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r * psif_ * interpolated_duzdr *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                                // Radial velocity component (sine part) U_k^S
                                // has no contribution

                              case 1:
                                break;

                                // Axial velocity component (cosine part) W_k^C

                              case 2:

                                // Add the mass matrix entries
                                if (flag == 2)
                                {
                                  mass_matrix(local_eqn, local_unknown) +=
                                    scaled_re_st * r * psif_ * testf_ * W *
                                    hang_weight * hang_weight2;
                                }

                                // Add contributions to the Jacobian matrix

                                // Stress tensor terms
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * dpsifdr * dtestfdr * W *
                                  hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * (1.0 + Gamma[1]) * dpsifdz *
                                  dtestfdz * W * hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * k * k * psif_ * testf_ * W *
                                  hang_weight * hang_weight2 / r;

                                // Inertial terms (du/dt)
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re_st * r *
                                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                                  psif_ * testf_ * W * hang_weight *
                                  hang_weight2;

                                // Inertial terms (convective)
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r *
                                  (interpolated_ur * dpsifdr +
                                   psif_ * interpolated_duzdz +
                                   interpolated_uz * dpsifdz) *
                                  testf_ * W * hang_weight * hang_weight2;


                                // Mesh velocity terms
                                if (!ALE_is_disabled)
                                {
                                  for (unsigned j = 0; j < 2; j++)
                                  {
                                    jacobian(local_eqn, local_unknown) +=
                                      scaled_re_st * r * mesh_velocity[j] *
                                      dpsifdx(l2, j) * testf_ * W *
                                      hang_weight * hang_weight2;
                                  }
                                }

                                break;

                                // Axial velocity component (sine part) W_k^S

                              case 3:

                                // Convective term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * k * interpolated_utheta * psif_ *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                                // Azimuthal velocity component (cosine part)
                                // V_k^C has no contribution

                              case 4:
                                break;

                                // Azimuthal velocity component (sine part)
                                // V_k^S

                              case 5:

                                // Stress tensor term
                                jacobian(local_eqn, local_unknown) +=
                                  visc_ratio * Gamma[1] * k * dpsifdz * testf_ *
                                  W * hang_weight * hang_weight2;

                                break;

                            } // End of third (axial) momentum equation
                            break;

                            // -------------------------------------------
                            // FOURTH (AXIAL) MOMENTUM EQUATION: SINE PART
                            // -------------------------------------------

                          case 3:

                            switch (i2)
                            {
                                // Radial velocity component (cosine part) U_k^C
                                // has no contribution

                              case 0:
                                break;

                                // Radial velocity component (sine part) U_k^S

                              case 1:

                                // Stress tensor term
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * Gamma[1] * dpsifdz *
                                  dtestfdr * W * hang_weight * hang_weight2;

                                // Convective term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r * psif_ * interpolated_duzdr *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                                // Axial velocity component (cosine part) W_k^S

                              case 2:

                                // Convective term
                                jacobian(local_eqn, local_unknown) +=
                                  scaled_re * k * interpolated_utheta * psif_ *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                                // Axial velocity component (sine part) W_k^S

                              case 3:

                                // Add the mass matrix entries
                                if (flag == 2)
                                {
                                  mass_matrix(local_eqn, local_unknown) +=
                                    scaled_re_st * r * psif_ * testf_ * W *
                                    hang_weight * hang_weight2;
                                }

                                // Add contributions to the Jacobian matrix

                                // Stress tensor terms
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * dpsifdr * dtestfdr * W *
                                  hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * (1.0 + Gamma[1]) * dpsifdz *
                                  dtestfdz * W * hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * k * k * psif_ * testf_ * W *
                                  hang_weight * hang_weight2 / r;

                                // Inertial terms (du/dt)
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re_st * r *
                                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                                  psif_ * testf_ * W * hang_weight *
                                  hang_weight2;

                                // Inertial terms (convective)
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r *
                                  (interpolated_ur * dpsifdr +
                                   psif_ * interpolated_duzdz +
                                   interpolated_uz * dpsifdz) *
                                  testf_ * W * hang_weight * hang_weight2;


                                // Mesh velocity terms
                                if (!ALE_is_disabled)
                                {
                                  for (unsigned j = 0; j < 2; j++)
                                  {
                                    jacobian(local_eqn, local_unknown) +=
                                      scaled_re_st * r * mesh_velocity[j] *
                                      dpsifdx(l2, j) * testf_ * W *
                                      hang_weight * hang_weight2;
                                  }
                                }

                                break;

                                // Azimuthal velocity component (cosine part)
                                // V_k^C

                              case 4:

                                // Stress tensor term
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * Gamma[1] * k * dpsifdz * testf_ *
                                  W * hang_weight * hang_weight2;

                                break;

                                // Azimuthal velocity component (sine part)
                                // V_k^S has no contribution

                              case 5:
                                break;

                            } // End of fourth (axial) momentum equation
                            break;

                            // ------------------------------------------------
                            // FIFTH (AZIMUTHAL) MOMENTUM EQUATION: COSINE PART
                            // ------------------------------------------------

                          case 4:

                            switch (i2)
                            {
                                // Radial velocity component (cosine part) U_k^C

                              case 0:

                                // Convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re *
                                  (r * interpolated_duthetadr +
                                   interpolated_utheta) *
                                  psif_ * testf_ * W * hang_weight *
                                  hang_weight2;

                                break;

                                // Radial velocity component (sine part) U_k^S

                              case 1:

                                // Stress tensor terms
                                jacobian(local_eqn, local_unknown) +=
                                  visc_ratio * k * psif_ *
                                  (((2.0 + Gamma[0]) * testf_ / r) -
                                   (Gamma[0] * dtestfdr)) *
                                  W * hang_weight * hang_weight2;

                                break;

                                // Axial velocity component (cosine part) W_k^C

                              case 2:

                                // Convective term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r * psif_ *
                                  interpolated_duthetadz * testf_ * W *
                                  hang_weight * hang_weight2;

                                break;

                                // Axial velocity component (sine part) W_k^S

                              case 3:

                                // Stress tensor term
                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * k * Gamma[0] * psif_ * dtestfdz *
                                  W * hang_weight * hang_weight2;

                                break;

                                // Azimuthal velocity component (cosine part)
                                // V_k^C

                              case 4:

                                // Add the mass matrix entries
                                if (flag == 2)
                                {
                                  mass_matrix(local_eqn, local_unknown) +=
                                    scaled_re_st * r * psif_ * testf_ * W *
                                    hang_weight * hang_weight2;
                                }

                                // Add contributions to the Jacobian matrix

                                // Stress tensor terms
                                jacobian(local_eqn, local_unknown) +=
                                  visc_ratio *
                                  (Gamma[0] * psif_ - r * dpsifdr) * dtestfdr *
                                  W * hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * dpsifdz * dtestfdz * W *
                                  hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) +=
                                  visc_ratio *
                                  (Gamma[0] * dpsifdr -
                                   (1.0 + (k * k) + (k * k * Gamma[0])) *
                                     psif_ / r) *
                                  testf_ * W * hang_weight * hang_weight2;

                                // Inertial terms (du/dt)
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re_st * r *
                                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                                  psif_ * testf_ * W * hang_weight *
                                  hang_weight2;

                                // Inertial terms (convective)
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re *
                                  (r * interpolated_ur * dpsifdr +
                                   interpolated_ur * psif_ +
                                   r * interpolated_uz * dpsifdz) *
                                  testf_ * W * hang_weight * hang_weight2;

                                // Mesh velocity terms
                                if (!ALE_is_disabled)
                                {
                                  for (unsigned j = 0; j < 2; j++)
                                  {
                                    jacobian(local_eqn, local_unknown) +=
                                      scaled_re_st * r * mesh_velocity[j] *
                                      dpsifdx(l2, j) * testf_ * W *
                                      hang_weight * hang_weight2;
                                  }
                                }

                                break;

                                // Azimuthal velocity component (sine part)
                                // V_k^S

                              case 5:

                                // Convective term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * k * interpolated_utheta * psif_ *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                            } // End of fifth (azimuthal) momentum equation
                            break;

                            // ----------------------------------------------
                            // SIXTH (AZIMUTHAL) MOMENTUM EQUATION: SINE PART
                            // ----------------------------------------------

                          case 5:

                            switch (i2)
                            {
                                // Radial velocity component (cosine part) U_k^C

                              case 0:

                                // Stress tensor terms
                                jacobian(local_eqn, local_unknown) +=
                                  visc_ratio * k * psif_ *
                                  ((Gamma[0] * dtestfdr) -
                                   ((2.0 + Gamma[0]) * testf_ / r)) *
                                  W * hang_weight * hang_weight2;

                                break;

                                // Radial velocity component (sine part) U_k^S

                              case 1:

                                // Convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re *
                                  (r * interpolated_duthetadr +
                                   interpolated_utheta) *
                                  psif_ * testf_ * W * hang_weight *
                                  hang_weight2;

                                break;

                                // Axial velocity component (cosine part) W_k^C

                              case 2:

                                // Stress tensor term
                                jacobian(local_eqn, local_unknown) +=
                                  visc_ratio * k * Gamma[0] * psif_ * dtestfdz *
                                  W * hang_weight * hang_weight2;

                                break;

                                // Axial velocity component (sine part) W_k^S

                              case 3:

                                // Convective term
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re * r * psif_ *
                                  interpolated_duthetadz * testf_ * W *
                                  hang_weight * hang_weight2;

                                break;

                                // Azimuthal velocity component (cosine part)
                                // V_k^C

                              case 4:

                                // Convective term
                                jacobian(local_eqn, local_unknown) +=
                                  scaled_re * k * interpolated_utheta * psif_ *
                                  testf_ * W * hang_weight * hang_weight2;

                                break;

                                // Azimuthal velocity component (sine part)
                                // V_k^S

                              case 5:

                                // Add the mass matrix entries
                                if (flag == 2)
                                {
                                  mass_matrix(local_eqn, local_unknown) +=
                                    scaled_re_st * r * psif_ * testf_ * W *
                                    hang_weight * hang_weight2;
                                }

                                // Add contributions to the Jacobian matrix

                                // Stress tensor terms
                                jacobian(local_eqn, local_unknown) +=
                                  visc_ratio *
                                  (Gamma[0] * psif_ - r * dpsifdr) * dtestfdr *
                                  W * hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) -=
                                  visc_ratio * r * dpsifdz * dtestfdz * W *
                                  hang_weight * hang_weight2;

                                jacobian(local_eqn, local_unknown) +=
                                  visc_ratio *
                                  (Gamma[0] * dpsifdr -
                                   (1.0 + (k * k) + (k * k * Gamma[0])) *
                                     psif_ / r) *
                                  testf_ * W * hang_weight * hang_weight2;

                                // Inertial terms (du/dt)
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re_st * r *
                                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                                  psif_ * testf_ * W * hang_weight *
                                  hang_weight2;

                                // Convective terms
                                jacobian(local_eqn, local_unknown) -=
                                  scaled_re *
                                  (r * interpolated_ur * dpsifdr +
                                   interpolated_ur * psif_ +
                                   r * interpolated_uz * dpsifdz) *
                                  testf_ * W * hang_weight * hang_weight2;

                                // Mesh velocity terms
                                if (!ALE_is_disabled)
                                {
                                  for (unsigned j = 0; j < 2; j++)
                                  {
                                    jacobian(local_eqn, local_unknown) +=
                                      scaled_re_st * r * mesh_velocity[j] *
                                      dpsifdx(l2, j) * testf_ * W *
                                      hang_weight * hang_weight2;
                                  }
                                }

                                break;

                            } // End of sixth (azimuthal) momentum equation
                            break;
                        }
                      } // End of if not boundary condition statement
                    } // End of loop over velocity components
                  } // End of loop over master (m2) nodes
                } // End of loop over the velocity nodes (l2)


                // Loop over the pressure shape functions
                for (unsigned l2 = 0; l2 < n_pres; l2++)
                {
                  // If the pressure dof is hanging
                  if (pressure_dof_is_hanging[l2])
                  {
                    // Pressure dof is hanging so it must be nodal-based
                    hang_info2_pt =
                      pressure_node_pt(l2)->hanging_pt(p_index[0]);

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
                    // Loop over the two pressure components
                    for (unsigned i2 = 0; i2 < 2; i2++)
                    {
                      // Get the number of the unknown
                      // -----------------------------

                      // If the pressure dof is hanging
                      if (pressure_dof_is_hanging[l2])
                      {
                        // Get the unknown from the node
                        local_unknown = local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), p_index[i2]);

                        // Get the weight from the hanging object
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      else
                      {
                        local_unknown = p_local_eqn(l2, i2);
                        hang_weight2 = 1.0;
                      }

                      // If the unknown is not pinned
                      if (local_unknown >= 0)
                      {
                        // Add in contributions to momentum equations
                        switch (i)
                        {
                            // ---------------------------------------------
                            // FIRST (RADIAL) MOMENTUM EQUATION: COSINE PART
                            // ---------------------------------------------

                          case 0:

                            switch (i2)
                            {
                                // Cosine part P_k^C
                              case 0:

                                jacobian(local_eqn, local_unknown) +=
                                  psip[l2] * (testf_ + r * dtestfdr) * W *
                                  hang_weight * hang_weight2;

                                break;

                                // Sine part P_k^S has no contribution
                              case 1:
                                break;

                            } // End of first (radial) momentum equation
                            break;

                            // --------------------------------------------
                            // SECOND (RADIAL) MOMENTUM EQUATION: SINE PART
                            // --------------------------------------------

                          case 1:

                            switch (i2)
                            {
                                // Cosine part P_k^C has no contribution
                              case 0:
                                break;

                                // Sine part P_k^S
                              case 1:

                                jacobian(local_eqn, local_unknown) +=
                                  psip[l2] * (testf_ + r * dtestfdr) * W *
                                  hang_weight * hang_weight2;

                                break;

                            } // End of second (radial) momentum equation
                            break;

                            // --------------------------------------------
                            // THIRD (AXIAL) MOMENTUM EQUATION: COSINE PART
                            // --------------------------------------------

                          case 2:

                            switch (i2)
                            {
                                // Cosine part P_k^C
                              case 0:

                                jacobian(local_eqn, local_unknown) +=
                                  r * psip[l2] * dtestfdz * W * hang_weight *
                                  hang_weight2;

                                break;

                                // Sine part P_k^S has no contribution
                              case 1:
                                break;

                            } // End of third (axial) momentum equation
                            break;

                            // -------------------------------------------
                            // FOURTH (AXIAL) MOMENTUM EQUATION: SINE PART
                            // -------------------------------------------

                          case 3:

                            switch (i2)
                            {
                                // Cosine part P_k^C has no contribution
                              case 0:
                                break;

                                // Sine part P_k^S
                              case 1:

                                jacobian(local_eqn, local_unknown) +=
                                  r * psip[l2] * dtestfdz * W * hang_weight *
                                  hang_weight2;

                                break;

                            } // End of fourth (axial) momentum equation
                            break;

                            // ------------------------------------------------
                            // FIFTH (AZIMUTHAL) MOMENTUM EQUATION: COSINE PART
                            // ------------------------------------------------

                          case 4:

                            switch (i2)
                            {
                                // Cosine part P_k^C has no contribution
                              case 0:
                                break;

                                // Sine part P_k^S
                              case 1:

                                jacobian(local_eqn, local_unknown) -=
                                  k * psip[l2] * testf_ * W * hang_weight *
                                  hang_weight2;

                                break;

                            } // End of fifth (azimuthal) momentum equation
                            break;

                            // ----------------------------------------------
                            // SIXTH (AZIMUTHAL) MOMENTUM EQUATION: SINE PART
                            // ----------------------------------------------

                          case 5:

                            switch (i2)
                            {
                                // Cosine part P_k^C
                              case 0:

                                jacobian(local_eqn, local_unknown) +=
                                  k * psip[l2] * testf_ * W * hang_weight *
                                  hang_weight2;

                                break;

                                // Sine part P_k^S has no contribution
                              case 1:
                                break;
                            } // End of sixth (azimuthal) momentum equation
                        } // End of add in contributions to different equations
                      } // End of if unknown is not pinned statement
                    } // End of loop over pressure components
                  } // End of loop over master (m2) nodes
                } // End of loop over pressure "nodes"
              } // End of Jacobian calculation
            } // End of if not boundary condition statement
          } // End of loop over velocity components
        } // End of loop over master (m) nodes
      } // End of loop over nodes

      // ====================
      // CONTINUITY EQUATIONS
      // ====================

      // Loop over the pressure shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Cache test function
        const double testp_ = testp[l];

        // If the pressure dof is hanging
        if (pressure_dof_is_hanging[l])
        {
          // Pressure dof is hanging so it must be nodal-based
          hang_info_pt = pressure_node_pt(l)->hanging_pt(p_index[0]);

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
          // Loop over the two pressure components
          for (unsigned i = 0; i < 2; i++)
          {
            // Get the equation number
            // -----------------------

            // If the pressure dof is hanging
            if (pressure_dof_is_hanging[l])
            {
              // Get the equation number from the master node
              local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                               p_index[i]);

              // Get the hang weight from the master node
              hang_weight = hang_info_pt->master_weight(m);
            }
            else
            {
              // Local equation number
              local_eqn = this->p_local_eqn(l, i);

              // Node contributes with full weight
              hang_weight = 1.0;
            }

            // If it's not a boundary condition...
            if (local_eqn >= 0)
            {
              switch (i)
              {
                  // --------------------------------------
                  // FIRST CONTINUITY EQUATION: COSINE PART
                  // --------------------------------------

                case 0:

                  // Gradient terms
                  residuals[local_eqn] +=
                    (interpolated_UC + r * interpolated_dUCdr +
                     k * interpolated_VS + r * interpolated_dWCdz) *
                    testp_ * W * hang_weight;

                  break;

                  // -------------------------------------
                  // SECOND CONTINUITY EQUATION: SINE PART
                  // -------------------------------------

                case 1:

                  // Gradient terms
                  residuals[local_eqn] +=
                    (interpolated_US + r * interpolated_dUSdr -
                     k * interpolated_VC + r * interpolated_dWSdz) *
                    testp_ * W * hang_weight;

                  break;
              }

              // ======================
              // Calculate the Jacobian
              // ======================
              if (flag)
              {
                // Number of master nodes
                unsigned n_master2 = 1;

                // Storage for the weight of the shape function
                double hang_weight2 = 1.0;

                // Loop over the velocity nodes for columns
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Cache velocity shape functions and their derivatives
                  const double psif_ = psif[l2];
                  const double dpsifdr = dpsifdx(l2, 0);
                  const double dpsifdz = dpsifdx(l2, 1);

                  // Local boolean to indicate whether the node is hanging
                  bool is_node2_hanging = node_pt(l2)->is_hanging();

                  // If the node is hanging
                  if (is_node2_hanging)
                  {
                    // Get the hanging pointer
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
                    for (unsigned i2 = 0; i2 < 6; i2++)
                    {
                      // Get the number of the unknown
                      // -----------------------------

                      // If the node is hanging
                      if (is_node2_hanging)
                      {
                        // Get the equation number from the master node
                        local_unknown = local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), u_nodal_index[i2]);

                        // Get the hang weight from the master node
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      // If the node is not hanging
                      else
                      {
                        // Local equation number
                        local_unknown =
                          this->nodal_local_eqn(l2, u_nodal_index[i2]);

                        // Node contributes with full weight
                        hang_weight2 = 1.0;
                      }

                      // If the unknown is not pinned
                      if (local_unknown >= 0)
                      {
                        // Different results for i and i2
                        switch (i)
                        {
                            // --------------------------------------
                            // FIRST CONTINUITY EQUATION: COSINE PART
                            // --------------------------------------

                          case 0:

                            switch (i2)
                            {
                                // Radial velocity component (cosine part) U_k^C

                              case 0:

                                jacobian(local_eqn, local_unknown) +=
                                  (psif_ + r * dpsifdr) * testp_ * W *
                                  hang_weight * hang_weight2;

                                break;

                                // Radial velocity component (sine part) U_k^S
                                // has no contribution

                              case 1:
                                break;

                                // Axial velocity component (cosine part) W_k^C

                              case 2:

                                jacobian(local_eqn, local_unknown) +=
                                  r * dpsifdz * testp_ * W * hang_weight *
                                  hang_weight2;

                                break;

                                // Axial velocity component (sine part) W_k^S
                                // has no contribution

                              case 3:
                                break;

                                // Azimuthal velocity component (cosine part)
                                // V_k^C has no contribution

                              case 4:
                                break;

                                // Azimuthal velocity component (sine part)
                                // V_k^S

                              case 5:

                                jacobian(local_eqn, local_unknown) +=
                                  k * psif_ * testp_ * W * hang_weight *
                                  hang_weight2;

                                break;

                            } // End of first continuity equation
                            break;

                            // -------------------------------------
                            // SECOND CONTINUITY EQUATION: SINE PART
                            // -------------------------------------

                          case 1:

                            switch (i2)
                            {
                                // Radial velocity component (cosine part) U_k^C
                                // has no contribution

                              case 0:
                                break;

                                // Radial velocity component (sine part) U_k^S

                              case 1:

                                jacobian(local_eqn, local_unknown) +=
                                  (psif_ + r * dpsifdr) * testp_ * W *
                                  hang_weight * hang_weight2;

                                break;

                                // Axial velocity component (cosine part) W_k^C
                                // has no contribution

                              case 2:
                                break;

                                // Axial velocity component (sine part) W_k^S

                              case 3:

                                jacobian(local_eqn, local_unknown) +=
                                  r * dpsifdz * testp_ * W * hang_weight *
                                  hang_weight2;

                                break;

                                // Azimuthal velocity component (cosine part)
                                // V_k^C

                              case 4:

                                jacobian(local_eqn, local_unknown) -=
                                  k * psif_ * testp_ * W * hang_weight *
                                  hang_weight2;

                                break;

                                // Azimuthal velocity component (sine part)
                                // V_k^S has no contribution

                              case 5:
                                break;

                            } // End of second continuity equation
                            break;
                        }
                      } // End of if unknown is not pinned statement
                    } // End of loop over velocity components
                  } // End of loop over master (m2) nodes
                } // End of loop over velocity nodes

                // Real and imaginary pressure components, P_k^C and P_k^S,
                // have no contribution to Jacobian

              } // End of Jacobian calculation
            } // End of if not boundary condition statement
          } // End of loop over the two pressure components
        } // End of loop over master nodes
      } // End of loop over pressure nodes

    } // End of loop over integration points

  } // End of fill_in_generic_residual_contribution_linearised_axi_nst

} // End of namespace oomph
