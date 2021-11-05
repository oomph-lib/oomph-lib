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
// Source file for refineable unsteady heat elements
#include "refineable_discontinuous_galerkin_space_time_unsteady_heat_elements.h"

/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////

namespace oomph
{
  //========================================================================
  /// Add element's contribution to the elemental residual vector and/or
  /// Jacobian matrix.
  /// flag=0: compute residual vector only
  /// flag=1: compute both
  //========================================================================
  template<unsigned SPATIAL_DIM>
  void RefineableSpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::
    fill_in_generic_residual_contribution_ust_heat(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
  {
    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Find the index at which the unknown is stored
    unsigned u_nodal_index = this->u_index_ust_heat();

    // Set up memory for the shape functions
    Shape psi(n_node);

    // Set up memory for the test functions
    Shape test(n_node);

    // Allocate space for the derivatives of the shape functions
    DShape dpsidx(n_node, SPATIAL_DIM + 1);

    // Allocate space for the derivatives of the test functions
    DShape dtestdx(n_node, SPATIAL_DIM + 1);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Storage for the local coordinates
    Vector<double> s(SPATIAL_DIM + 1, 0.0);

    // Get the Alpha parameter
    double alpha_local = this->alpha();

    // Get the Beta parameter
    double beta_local = this->beta();

    // Integer to hold the local equation
    int local_eqn = 0;

    // Integer to hold the local unknowns
    int local_unknown = 0;

    // Local storage for pointer to hang info object
    HangInfo* hang_info_pt = 0;

    // Local storage for pointer to (another) hang info object
    HangInfo* hang_info2_pt = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < SPATIAL_DIM + 1; i++)
      {
        // Calculate the i-th local coordinate
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = this->dshape_and_dtest_eulerian_at_knot_ust_heat(
        ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Storage for the interpolated time value
      double interpolated_t = 0.0;

      // Storage for the interpolated solution value
      double interpolated_u = 0.0;

      // Storage for the interpolated time-derivative of the solution
      double interpolated_dudt = 0.0;

      // Storage for the spatial coordinates
      Vector<double> interpolated_x(SPATIAL_DIM, 0.0);

      // Storage for the spatial derivatives of the solution
      Vector<double> interpolated_dudx(SPATIAL_DIM + 1, 0.0);

      // Storage for the mesh velocity
      Vector<double> mesh_velocity(SPATIAL_DIM, 0.0);

      //-------------------------------------------------
      // Calculate derivatives and source function value:
      //-------------------------------------------------
      // Loop over the nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the nodal value at the l-th node
        double u_value = nodal_value(l, u_nodal_index);

        // Update the interpolated u value
        interpolated_u += u_value * psi(l);

        // Loop over the coordinate directions (both spatial AND time)
        for (unsigned j = 0; j < SPATIAL_DIM; j++)
        {
          // Update the interpolated x value
          interpolated_x[j] += nodal_position(l, j) * psi(l);

          // Update the interpolated du/dx_j value
          interpolated_dudx[j] += u_value * dpsidx(l, j);
        }

        // Update the interpolated time value
        interpolated_t += nodal_position(l, SPATIAL_DIM) * psi(l);

        // Update the interpolated du/dt value
        interpolated_dudt += u_value * dpsidx(l, SPATIAL_DIM);
      } // for (unsigned l=0;l<n_node;l++)

      // If ALE is enabled
      if (!(this->ALE_is_disabled))
      {
        // Loop over the coordinate directions
        for (unsigned j = 0; j < SPATIAL_DIM; j++)
        {
          // Loop over the nodes
          for (unsigned l = 0; l < n_node; l++)
          {
            // Update the mesh velocity
            mesh_velocity[j] += nodal_position(l, j) * dpsidx(l, SPATIAL_DIM);
          }
        } // for (unsigned j=0;j<SPATIAL_DIM;j++)
      } // if (!ALE_is_disabled)

      // Initialise the source term value
      double source = 0.0;

      // Get the interpolated source term value
      this->get_source_ust_heat(interpolated_t, ipt, interpolated_x, source);

      //--------------------------------
      // Assemble residuals and Jacobian
      //--------------------------------
      // Loop over the nodes (or equivalently the test functions)
      for (unsigned l = 0; l < n_node; l++)
      {
        // Storage for the number of master nodes
        unsigned n_master = 1;

        // Storage for the hang weight associated with the shape function
        double hang_weight = 1.0;

        // Check if the node is hanging
        bool is_node_hanging = this->node_pt(l)->is_hanging();

        // If the node is hanging, get the number of master nodes
        if (is_node_hanging)
        {
          // Get the hang info pointer associated with the l-th node in the
          // element
          hang_info_pt = this->node_pt(l)->hanging_pt();

          // Get the number of master nodes associated with the node
          n_master = hang_info_pt->nmaster();
        }
        // If it's not hanging there is just one master node; the node itself
        else
        {
          // Set the number of master nodes to one
          n_master = 1;
        }

        // Loop over the number of master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Check if the node is hanging
          if (is_node_hanging)
          {
            // Read out the local equation from the master node
            local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                             u_nodal_index);

            // Read out the weight from the master node
            hang_weight = hang_info_pt->master_weight(m);
          }
          // If the node is not hanging
          else
          {
            // The local equation number comes from the node itself
            local_eqn = this->nodal_local_eqn(l, u_nodal_index);

            // The hang weight is one
            hang_weight = 1.0;
          }

          // If the nodal equation is not a boundary condition
          if (local_eqn >= 0)
          {
            // Add source term and time derivative
            residuals[local_eqn] +=
              ((source + alpha_local * interpolated_dudt) * test(l) * W *
               hang_weight);

            // If ALE is enabled
            if (!(this->ALE_is_disabled))
            {
              // Loop over the coordinate directions
              for (unsigned k = 0; k < SPATIAL_DIM; k++)
              {
                // Add in the mesh velocity contribution
                residuals[local_eqn] -=
                  (alpha_local * mesh_velocity[k] * interpolated_dudx[k] *
                   test(l) * W * hang_weight);
              }
            } // if (!ALE_is_disabled)

            // Loop over the coordinate directions
            for (unsigned k = 0; k < SPATIAL_DIM; k++)
            {
              // Add in the contribution from the Laplace operator
              residuals[local_eqn] += (beta_local * interpolated_dudx[k] *
                                       dtestdx(l, k) * W * hang_weight);
            }

            //-----------------------
            // Calculate the Jacobian
            //-----------------------
            // If we also need to construct the Jacobian
            if (flag)
            {
              // Storage for the number of master nodes
              unsigned n_master2 = 1;

              // Storage for the hang weight associated with the shape function
              double hang_weight2 = 1.0;

              // Loop over the nodes for the variables
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Check if the node is hanging
                bool is_node2_hanging = this->node_pt(l2)->is_hanging();

                // If the node is hanging, get the number of master nodes
                if (is_node2_hanging)
                {
                  // Get the hang info pointer associated with the l2-th node
                  hang_info2_pt = this->node_pt(l2)->hanging_pt();

                  // Get the number of master nodes associated with the node
                  n_master2 = hang_info2_pt->nmaster();
                }
                // If it's not hanging there is just one master node; the node
                // itself
                else
                {
                  // Set the number of master nodes to one
                  n_master2 = 1;
                }

                // Loop over the master nodes
                for (unsigned m2 = 0; m2 < n_master2; m2++)
                {
                  // Check if the node is hanging
                  if (is_node2_hanging)
                  {
                    // Read out the local equation from the master node
                    local_unknown = this->local_hang_eqn(
                      hang_info2_pt->master_node_pt(m2), u_nodal_index);

                    // Read out the weight from the master node
                    hang_weight2 = hang_info2_pt->master_weight(m2);
                  }
                  // If the node is not hanging
                  else
                  {
                    // The local equation number comes from the node itself
                    local_unknown = this->nodal_local_eqn(l2, u_nodal_index);

                    // The hang weight is one
                    hang_weight2 = 1.0;
                  }

                  // If the unknown is not pinned
                  if (local_unknown >= 0)
                  {
                    // Time-derivative contribution
                    jacobian(local_eqn, local_unknown) +=
                      (alpha_local * test(l) * dpsidx(l2, SPATIAL_DIM) * W *
                       hang_weight * hang_weight2);

                    // Laplace operator contribution
                    for (unsigned k = 0; k < SPATIAL_DIM; k++)
                    {
                      // Add the Laplacian contribution to the elemental
                      // Jacobian
                      jacobian(local_eqn, local_unknown) +=
                        (beta_local * dpsidx(l2, k) * dtestdx(l, k) * W *
                         hang_weight * hang_weight2);
                    }

                    // If ALE is enabled
                    if (!(this->ALE_is_disabled))
                    {
                      // Loop over the spatial coordinates
                      for (unsigned k = 0; k < SPATIAL_DIM; k++)
                      {
                        // Add the ALE contribution to the Jacobian
                        jacobian(local_eqn, local_unknown) -=
                          (alpha_local * mesh_velocity[k] * dpsidx(l2, k) *
                           test(l) * W * hang_weight * hang_weight2);
                      }
                    } // if (!(this->ALE_is_disabled))
                  } // if (local_unknown>=0)
                } // for (unsigned m2=0;m2<n_master;m2++)
              } // for (unsigned l2=0;l2<n_node;l2++)
            } // if (flag)
          } // if (local_eqn>=0)
        } // for (unsigned m=0;m<n_master;m++)
      } // for (unsigned l=0;l<n_node;l++)
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of fill_in_generic_residual_contribution_ust_heat

  //=========================================================================
  // Force build of templates:
  // Build 2D (in space) elements --> 3D (space-time) elements
  //=========================================================================
  template class RefineableQUnsteadyHeatSpaceTimeElement<2, 2>;
  template class RefineableQUnsteadyHeatSpaceTimeElement<2, 3>;
  template class RefineableQUnsteadyHeatSpaceTimeElement<2, 4>;
} // namespace oomph
