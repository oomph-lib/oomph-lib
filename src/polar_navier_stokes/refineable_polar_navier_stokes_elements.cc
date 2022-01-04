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
#include "refineable_polar_navier_stokes_elements.h"

namespace oomph
{
  /// ///////////////////////////////////////////////////////////////////////
  //======================================================================//
  /// Start of what would've been refineable_navier_stokes_elements.cc    //
  //======================================================================//
  /// ///////////////////////////////////////////////////////////////////////

  //==============================================================
  ///  Compute the residuals for the Navier--Stokes
  ///  equations; flag=1(or 0): do (or don't) compute the
  ///  Jacobian as well.
  ///  flag=2 for Residuals, Jacobian and mass_matrix
  ///
  ///  This is now my new version with Jacobian and
  ///  dimensionless phi
  ///
  ///  This version supports hanging nodes
  //==============================================================
  void RefineablePolarNavierStokesEquations::
    fill_in_generic_residual_contribution(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian,
                                          DenseMatrix<double>& mass_matrix,
                                          unsigned flag)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many pressure dofs there are
    unsigned n_pres = npres_pnst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[2];
    for (unsigned i = 0; i < 2; i++)
    {
      u_nodal_index[i] = u_index_pnst(i);
    }

    // Which nodal value represents the pressure? (Negative if pressure
    // is not based on nodal interpolation).
    int p_index = this->p_nodal_index_pnst();

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

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Get the reynolds number and Alpha
    const double Re = re();
    const double Alpha = alpha();
    const double Re_St = re_st();

    // Integers to store the local equations and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Pointers to hang info objects
    HangInfo *hang_info_pt = 0, *hang_info2_pt = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++) s[i] = integral_pt()->knot(ipt, i);
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_pnst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Call the pressure shape and test functions
      this->pshape_pnst(s, psip, testp);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the pressure and velocity components
      // Allocate storage initialised to zero
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(2, 0.0);
      Vector<double> interpolated_x(2, 0.0);
      // Vector<double> dudt(2);
      DenseMatrix<double> interpolated_dudx(2, 2, 0.0);

      // Initialise to zero
      for (unsigned i = 0; i < 2; i++)
      {
        // dudt[i] = 0.0;
        interpolated_u[i] = 0.0;
        interpolated_x[i] = 0.0;
        for (unsigned j = 0; j < 2; j++)
        {
          interpolated_dudx(i, j) = 0.0;
        }
      }

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
        interpolated_p += this->p_pnst(l) * psip[l];

      // Calculate velocities and derivatives:

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned i = 0; i < 2; i++)
        {
          // Get the nodal value
          double u_value = this->nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif[l];
          interpolated_x[i] += this->nodal_position(l, i) * psif[l];
          // dudt[i]+=du_dt_pnst(l,i)*psif[l];

          // Loop over derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      // MOMENTUM EQUATIONS
      //------------------
      // Number of master nodes and storage for the weight of the shape function
      unsigned n_master = 1;
      double hang_weight = 1.0;
      // Loop over the nodes for the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
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

        // Now add in a new loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Can't loop over velocity components as don't have identical
          // contributions Do seperately for i = {0,1} instead
          unsigned i = 0;
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

            /*IF it's not a boundary condition*/
            if (local_eqn >= 0)
            {
              // Add the testf[l] term of the stress tensor
              residuals[local_eqn] +=
                ((interpolated_p / interpolated_x[0]) -
                 ((1. + Gamma[i]) / pow(interpolated_x[0], 2.)) *
                   ((1. / Alpha) * interpolated_dudx(1, 1) +
                    interpolated_u[0])) *
                testf[l] * interpolated_x[0] * Alpha * W * hang_weight;

              // Add the dtestfdx(l,0) term of the stress tensor
              residuals[local_eqn] +=
                (interpolated_p - (1. + Gamma[i]) * interpolated_dudx(0, 0)) *
                dtestfdx(l, 0) * interpolated_x[0] * Alpha * W * hang_weight;

              // Add the dtestfdx(l,1) term of the stress tensor
              residuals[local_eqn] -=
                ((1. / (interpolated_x[0] * Alpha)) * interpolated_dudx(0, 1) -
                 (interpolated_u[1] / interpolated_x[0]) +
                 Gamma[i] * interpolated_dudx(1, 0)) *
                (1. / (interpolated_x[0] * Alpha)) * dtestfdx(l, 1) *
                interpolated_x[0] * Alpha * W * hang_weight;

              // Convective terms
              residuals[local_eqn] -=
                Re *
                (interpolated_u[0] * interpolated_dudx(0, 0) +
                 (interpolated_u[1] / (interpolated_x[0] * Alpha)) *
                   interpolated_dudx(0, 1) -
                 (pow(interpolated_u[1], 2.) / interpolated_x[0])) *
                testf[l] * interpolated_x[0] * Alpha * W * hang_weight;


              // CALCULATE THE JACOBIAN
              if (flag)
              {
                // Number of master nodes and weights
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                // Loop over the velocity shape functions again
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
                    // Again can't loop over velocity components due to loss of
                    // symmetry
                    unsigned i2 = 0;
                    {
                      // Get the number of the unknown
                      // If the node is hanging
                      if (is_node2_hanging)
                      {
                        // Get the equation number from the master node
                        local_unknown = this->local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), u_nodal_index[i2]);
                        // Get the hang weights
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      else
                      {
                        local_unknown =
                          this->nodal_local_eqn(l2, u_nodal_index[i2]);
                        hang_weight2 = 1.0;
                      }

                      // If at a non-zero degree of freedom add in the entry
                      if (local_unknown >= 0)
                      {
                        // Add contribution to Elemental Matrix
                        jacobian(local_eqn, local_unknown) -=
                          (1. + Gamma[i]) *
                          (psif[l2] / pow(interpolated_x[0], 2.)) * testf[l] *
                          interpolated_x[0] * Alpha * W * hang_weight *
                          hang_weight2;

                        jacobian(local_eqn, local_unknown) -=
                          (1. + Gamma[i]) * dpsifdx(l2, 0) * dtestfdx(l, 0) *
                          interpolated_x[0] * Alpha * W * hang_weight *
                          hang_weight2;

                        jacobian(local_eqn, local_unknown) -=
                          (1. / (interpolated_x[0] * Alpha)) * dpsifdx(l2, 1) *
                          (1. / (interpolated_x[0] * Alpha)) * dtestfdx(l, 1) *
                          interpolated_x[0] * Alpha * W * hang_weight *
                          hang_weight2;

                        // Now add in the inertial terms
                        jacobian(local_eqn, local_unknown) -=
                          Re *
                          (psif[l2] * interpolated_dudx(0, 0) +
                           interpolated_u[0] * dpsifdx(l2, 0) +
                           (interpolated_u[1] / (interpolated_x[0] * Alpha)) *
                             dpsifdx(l2, 1)) *
                          testf[l] * interpolated_x[0] * Alpha * W *
                          hang_weight * hang_weight2;

                        // extra bit for mass matrix
                        if (flag == 2)
                        {
                          mass_matrix(local_eqn, local_unknown) +=
                            Re_St * psif[l2] * testf[l] * interpolated_x[0] *
                            Alpha * W * hang_weight * hang_weight2;
                        }

                      } // End of (Jacobian's) if not boundary condition
                        // statement
                    } // End of i2=0 section

                    i2 = 1;
                    {
                      // Get the number of the unknown
                      // If the node is hanging
                      if (is_node2_hanging)
                      {
                        // Get the equation number from the master node
                        local_unknown = this->local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), u_nodal_index[i2]);
                        // Get the hang weights
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      else
                      {
                        local_unknown =
                          this->nodal_local_eqn(l2, u_nodal_index[i2]);
                        hang_weight2 = 1.0;
                      }

                      // If at a non-zero degree of freedom add in the entry
                      if (local_unknown >= 0)
                      {
                        // Add contribution to Elemental Matrix
                        jacobian(local_eqn, local_unknown) -=
                          ((1. + Gamma[i]) /
                           (pow(interpolated_x[0], 2.) * Alpha)) *
                          dpsifdx(l2, 1) * testf[l] * interpolated_x[0] *
                          Alpha * W * hang_weight * hang_weight2;

                        jacobian(local_eqn, local_unknown) -=
                          (-(psif[l2] / interpolated_x[0]) +
                           Gamma[i] * dpsifdx(l2, 0)) *
                          (1. / (interpolated_x[0] * Alpha)) * dtestfdx(l, 1) *
                          interpolated_x[0] * Alpha * W * hang_weight *
                          hang_weight2;

                        // Now add in the inertial terms
                        jacobian(local_eqn, local_unknown) -=
                          Re *
                          ((psif[l2] / (interpolated_x[0] * Alpha)) *
                             interpolated_dudx(0, 1) -
                           2 * interpolated_u[1] *
                             (psif[l2] / interpolated_x[0])) *
                          testf[l] * interpolated_x[0] * Alpha * W *
                          hang_weight * hang_weight2;

                      } // End of (Jacobian's) if not boundary condition
                        // statement
                    } // End of i2=1 section

                  } // End of loop over master nodes m2
                } // End of l2 loop

                /*Now loop over pressure shape functions*/
                /*This is the contribution from pressure gradient*/
                for (unsigned l2 = 0; l2 < n_pres; l2++)
                {
                  // If the pressure dof is hanging
                  if (pressure_dof_is_hanging[l2])
                  {
                    hang_info2_pt =
                      this->pressure_node_pt(l2)->hanging_pt(p_index);
                    // Pressure dof is hanging so it must be nodal-based
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
                      // Get the unknown from the master node
                      local_unknown = this->local_hang_eqn(
                        hang_info2_pt->master_node_pt(m2), p_index);
                      // Get the weight from the hanging object
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    else
                    {
                      local_unknown = this->p_local_eqn(l2);
                      hang_weight2 = 1.0;
                    }

                    /*If we are at a non-zero degree of freedom in the entry*/
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) +=
                        (psip[l2] / interpolated_x[0]) * testf[l] *
                        interpolated_x[0] * Alpha * W * hang_weight *
                        hang_weight2;

                      jacobian(local_eqn, local_unknown) +=
                        psip[l2] * dtestfdx(l, 0) * interpolated_x[0] * Alpha *
                        W * hang_weight * hang_weight2;

                    } // End of Jacobian pressure if not a boundary condition
                      // statement

                  } // End of loop over master nodes m2
                } // End of loop over pressure shape functions l2

              } /*End of Jacobian calculation*/

            } // End of if not boundary condition statement
          } // End of i=0 section

          i = 1;
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

            /*IF it's not a boundary condition*/
            if (local_eqn >= 0)
            {
              // Add the testf[l] term of the stress tensor
              residuals[local_eqn] +=
                ((1. / (pow(interpolated_x[0], 2.) * Alpha)) *
                   interpolated_dudx(0, 1) -
                 (interpolated_u[1] / pow(interpolated_x[0], 2.)) +
                 Gamma[i] * (1. / interpolated_x[0]) *
                   interpolated_dudx(1, 0)) *
                testf[l] * interpolated_x[0] * Alpha * W * hang_weight;

              // Add the dtestfdx(l,0) term of the stress tensor
              residuals[local_eqn] -=
                (interpolated_dudx(1, 0) +
                 Gamma[i] * ((1. / (interpolated_x[0] * Alpha)) *
                               interpolated_dudx(0, 1) -
                             (interpolated_u[1] / interpolated_x[0]))) *
                dtestfdx(l, 0) * interpolated_x[0] * Alpha * W * hang_weight;

              // Add the dtestfdx(l,1) term of the stress tensor
              residuals[local_eqn] +=
                (interpolated_p -
                 (1. + Gamma[i]) * ((1. / (interpolated_x[0] * Alpha)) *
                                      interpolated_dudx(1, 1) +
                                    (interpolated_u[0] / interpolated_x[0]))) *
                (1. / (interpolated_x[0] * Alpha)) * dtestfdx(l, 1) *
                interpolated_x[0] * Alpha * W * hang_weight;

              // Convective terms
              residuals[local_eqn] -=
                Re *
                (interpolated_u[0] * interpolated_dudx(1, 0) +
                 (interpolated_u[1] / (interpolated_x[0] * Alpha)) *
                   interpolated_dudx(1, 1) +
                 ((interpolated_u[0] * interpolated_u[1]) /
                  interpolated_x[0])) *
                testf[l] * interpolated_x[0] * Alpha * W * hang_weight;

              // CALCULATE THE JACOBIAN
              if (flag)
              {
                // Number of master nodes and weights
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;

                // Loop over the velocity shape functions again
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
                    // Again can't loop over velocity components due to loss of
                    // symmetry
                    unsigned i2 = 0;
                    {
                      // Get the number of the unknown
                      // If the node is hanging
                      if (is_node2_hanging)
                      {
                        // Get the equation number from the master node
                        local_unknown = this->local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), u_nodal_index[i2]);
                        // Get the hang weights
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      else
                      {
                        local_unknown =
                          this->nodal_local_eqn(l2, u_nodal_index[i2]);
                        hang_weight2 = 1.0;
                      }

                      // If at a non-zero degree of freedom add in the entry
                      if (local_unknown >= 0)
                      {
                        // Add contribution to Elemental Matrix
                        jacobian(local_eqn, local_unknown) +=
                          (1. / (pow(interpolated_x[0], 2.) * Alpha)) *
                          dpsifdx(l2, 1) * testf[l] * interpolated_x[0] *
                          Alpha * W * hang_weight * hang_weight2;

                        jacobian(local_eqn, local_unknown) -=
                          Gamma[i] * (1. / (interpolated_x[0] * Alpha)) *
                          dpsifdx(l2, 1) * dtestfdx(l, 0) * interpolated_x[0] *
                          Alpha * W * hang_weight * hang_weight2;

                        jacobian(local_eqn, local_unknown) -=
                          (1 + Gamma[i]) * (psif[l2] / interpolated_x[0]) *
                          (1. / (interpolated_x[0] * Alpha)) * dtestfdx(l, 1) *
                          interpolated_x[0] * Alpha * W * hang_weight *
                          hang_weight2;

                        // Now add in the inertial terms
                        jacobian(local_eqn, local_unknown) -=
                          Re *
                          (psif[l2] * interpolated_dudx(1, 0) +
                           (psif[l2] * interpolated_u[1] / interpolated_x[0])) *
                          testf[l] * interpolated_x[0] * Alpha * W *
                          hang_weight * hang_weight2;

                      } // End of (Jacobian's) if not boundary condition
                        // statement
                    } // End of i2=0 section

                    i2 = 1;
                    {
                      // Get the number of the unknown
                      // If the node is hanging
                      if (is_node2_hanging)
                      {
                        // Get the equation number from the master node
                        local_unknown = this->local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), u_nodal_index[i2]);
                        // Get the hang weights
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      else
                      {
                        local_unknown =
                          this->nodal_local_eqn(l2, u_nodal_index[i2]);
                        hang_weight2 = 1.0;
                      }

                      // If at a non-zero degree of freedom add in the entry
                      if (local_unknown >= 0)
                      {
                        // Add contribution to Elemental Matrix
                        jacobian(local_eqn, local_unknown) +=
                          (-(psif[l2] / pow(interpolated_x[0], 2.)) +
                           Gamma[i] * (1. / interpolated_x[0]) *
                             dpsifdx(l2, 0)) *
                          testf[l] * interpolated_x[0] * Alpha * W *
                          hang_weight * hang_weight2;

                        jacobian(local_eqn, local_unknown) -=
                          (dpsifdx(l2, 0) -
                           Gamma[i] * (psif[l2] / interpolated_x[0])) *
                          dtestfdx(l, 0) * interpolated_x[0] * Alpha * W *
                          hang_weight * hang_weight2;

                        jacobian(local_eqn, local_unknown) -=
                          (1. + Gamma[i]) * (1. / (interpolated_x[0] * Alpha)) *
                          dpsifdx(l2, 1) * (1. / (interpolated_x[0] * Alpha)) *
                          dtestfdx(l, 1) * interpolated_x[0] * Alpha * W *
                          hang_weight * hang_weight2;

                        // Now add in the inertial terms
                        jacobian(local_eqn, local_unknown) -=
                          Re *
                          (interpolated_u[0] * dpsifdx(l2, 0) +
                           (psif[l2] / (interpolated_x[0] * Alpha)) *
                             interpolated_dudx(1, 1) +
                           (interpolated_u[1] / (interpolated_x[0] * Alpha)) *
                             dpsifdx(l2, 1) +
                           (interpolated_u[0] * psif[l2] / interpolated_x[0])) *
                          testf[l] * interpolated_x[0] * Alpha * W *
                          hang_weight * hang_weight2;

                        // extra bit for mass matrix
                        if (flag == 2)
                        {
                          mass_matrix(local_eqn, local_unknown) +=
                            Re_St * psif[l2] * testf[l] * interpolated_x[0] *
                            Alpha * W * hang_weight * hang_weight2;
                        }

                      } // End of (Jacobian's) if not boundary condition
                        // statement
                    } // End of i2=1 section

                  } // End of loop over master nodes m2
                } // End of l2 loop

                /*Now loop over pressure shape functions*/
                /*This is the contribution from pressure gradient*/
                for (unsigned l2 = 0; l2 < n_pres; l2++)
                {
                  // If the pressure dof is hanging
                  if (pressure_dof_is_hanging[l2])
                  {
                    hang_info2_pt =
                      this->pressure_node_pt(l2)->hanging_pt(p_index);
                    // Pressure dof is hanging so it must be nodal-based
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
                      // Get the unknown from the master node
                      local_unknown = this->local_hang_eqn(
                        hang_info2_pt->master_node_pt(m2), p_index);
                      // Get the weight from the hanging object
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    else
                    {
                      local_unknown = this->p_local_eqn(l2);
                      hang_weight2 = 1.0;
                    }


                    /*If we are at a non-zero degree of freedom in the entry*/
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) +=
                        (psip[l2] / interpolated_x[0]) * (1. / Alpha) *
                        dtestfdx(l, 1) * interpolated_x[0] * Alpha * W *
                        hang_weight * hang_weight2;

                    } // End of if not boundary condition for pressure in
                      // jacobian

                  } // End of loop over master nodes m2
                } // End of loop over pressure test functions l2

              } /*End of Jacobian calculation*/

            } // End of if not boundary condition statement
          } // End of i=1 section

        } // End of loop over master nodes m
      } // End of loop over shape functions


      // CONTINUITY EQUATION
      //-------------------

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

          // If not a boundary conditions
          if (local_eqn >= 0)
          {
            residuals[local_eqn] +=
              (interpolated_dudx(0, 0) +
               (interpolated_u[0] / interpolated_x[0]) +
               (1. / (interpolated_x[0] * Alpha)) * interpolated_dudx(1, 1)) *
              testp[l] * interpolated_x[0] * Alpha * W * hang_weight;

            /*CALCULATE THE JACOBIAN*/
            if (flag)
            {
              unsigned n_master2 = 1;
              double hang_weight2 = 1.0;
              /*Loop over the velocity shape functions*/
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
                  unsigned i2 = 0;
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
                    /*If we're at a non-zero degree of freedom add it in*/
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) +=
                        (dpsifdx(l2, 0) + (psif[l2] / interpolated_x[0])) *
                        testp[l] * interpolated_x[0] * Alpha * W * hang_weight *
                        hang_weight2;
                    }
                  } // End of i2=0 section

                  i2 = 1;
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

                    /*If we're at a non-zero degree of freedom add it in*/
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) +=
                        (1. / (interpolated_x[0] * Alpha)) * dpsifdx(l2, 1) *
                        testp[l] * interpolated_x[0] * Alpha * W * hang_weight *
                        hang_weight2;
                    }
                  } // End of i2=1 section

                } // End of loop over master nodes m2
              } /*End of loop over l2*/
            } /*End of Jacobian calculation*/

          } // End of if not boundary condition
        } // End of loop over master nodes m
      } // End of loop over pressure test functions l

    } // End of loop over integration points

  } // End of add_generic_residual_contribution


} // End of namespace oomph
