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
#include "refineable_gen_advection_diffusion_elements.h"

namespace oomph
{
  //==========================================================================
  /// Add  the element's contribution to the elemental residual vector
  /// and/or elemental jacobian matrix.
  /// This function overloads the standard version so that the possible
  /// presence of hanging nodes is taken into account.
  //=========================================================================
  template<unsigned DIM>
  void RefineableGeneralisedAdvectionDiffusionEquations<DIM>::
    fill_in_generic_residual_contribution_cons_adv_diff(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Get the nodal index at which the unknown is stored
    unsigned u_nodal_index = this->u_index_cons_adv_diff();

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get Peclet number
    double peclet = this->pe();

    // Get the Peclet multiplied by the Strouhal number
    double peclet_st = this->pe_st();

    // Integers used to store the local equation number and local unknown
    // indices for the residuals and jacobians
    int local_eqn = 0, local_unknown = 0;

    // Local storage for pointers to hang_info objects
    HangInfo *hang_info_pt = 0, *hang_info2_pt = 0;

    // Local variable to determine the ALE stuff
    bool ALE_is_disabled_flag = this->ALE_is_disabled;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++) s[i] = integral_pt()->knot(ipt, i);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = this->dshape_and_dtest_eulerian_at_knot_cons_adv_diff(
        ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the function, initialise to zero
      double dudt = 0.0;
      double interpolated_u = 0.0;

      // These need to be a Vector to be ANSI C++, initialise to zero
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> interpolated_dudx(DIM, 0.0);
      Vector<double> mesh_velocity(DIM, 0.0);

      // Calculate function value and derivatives:
      //-----------------------------------------

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the value at the node
        double u_value = this->nodal_value(l, u_nodal_index);
        interpolated_u += u_value * psi(l);
        dudt += this->du_dt_cons_adv_diff(l) * psi(l);
        // Loop over directions
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_x[j] += nodal_position(l, j) * psi(l);
          interpolated_dudx[j] += u_value * dpsidx(l, j);
        }
      }

      // Get the mesh velocity, if required
      if (!ALE_is_disabled_flag)
      {
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned j = 0; j < DIM; j++)
          {
            mesh_velocity[j] += dnodal_position_dt(l, j) * psi(l);
          }
        }
      }


      // Get body force
      double source;
      this->get_source_cons_adv_diff(ipt, interpolated_x, source);


      // Get wind
      //--------
      Vector<double> wind(DIM);
      this->get_wind_cons_adv_diff(ipt, s, interpolated_x, wind);

      // Get the conserved wind (non-divergence free)
      Vector<double> conserved_wind(DIM);
      this->get_conserved_wind_cons_adv_diff(
        ipt, s, interpolated_x, conserved_wind);

      // Get diffusivity tensor
      DenseMatrix<double> D(DIM, DIM);
      this->get_diff_cons_adv_diff(ipt, s, interpolated_x, D);

      // Assemble residuals and Jacobian
      //================================

      // Loop over the nodes for the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Local variables to store the number of master nodes and
        // the weight associated with the shape function if the node is hanging
        unsigned n_master = 1;
        double hang_weight = 1.0;
        // Local bool (is the node hanging)
        bool is_node_hanging = this->node_pt(l)->is_hanging();

        // If the node is hanging, get the number of master nodes
        if (is_node_hanging)
        {
          hang_info_pt = this->node_pt(l)->hanging_pt();
          n_master = hang_info_pt->nmaster();
        }
        // Otherwise there is just one master node, the node itself
        else
        {
          n_master = 1;
        }

        // Loop over the number of master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Get the local equation number and hang_weight
          // If the node is hanging
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

          // If the nodal equation is not a boundary conditino
          if (local_eqn >= 0)
          {
            // Add du/dt and body force/source term here
            residuals[local_eqn] -=
              (peclet_st * dudt + source) * test(l) * W * hang_weight;

            // The Mesh velocity and GeneralisedAdvection--Diffusion bit
            for (unsigned k = 0; k < DIM; k++)
            {
              // Terms that multiply the test function
              double tmp = peclet * wind[k];
              // If the mesh is moving need to subtract the mesh velocity
              if (!ALE_is_disabled_flag)
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
              residuals[local_eqn] -=
                (tmp * test(l) + tmp2 * dtestdx(l, k)) * W * hang_weight;
            }

            // Calculate the Jacobian
            if (flag)
            {
              // Local variables to store the number of master nodes
              // and the weights associated with each hanging node
              unsigned n_master2 = 1;
              double hang_weight2 = 1.0;
              // Loop over the nodes for the variables
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Local bool (is the node hanging)
                bool is_node2_hanging = this->node_pt(l2)->is_hanging();
                // If the node is hanging, get the number of master nodes
                if (is_node2_hanging)
                {
                  hang_info2_pt = this->node_pt(l2)->hanging_pt();
                  n_master2 = hang_info2_pt->nmaster();
                }
                // Otherwise there is one master node, the node itself
                else
                {
                  n_master2 = 1;
                }

                // Loop over the master nodes
                for (unsigned m2 = 0; m2 < n_master2; m2++)
                {
                  // Get the local unknown and weight
                  // If the node is hanging
                  if (is_node2_hanging)
                  {
                    // Read out the local unknown from the master node
                    local_unknown = this->local_hang_eqn(
                      hang_info2_pt->master_node_pt(m2), u_nodal_index);
                    // Read out the hanging weight from the master node
                    hang_weight2 = hang_info2_pt->master_weight(m2);
                  }
                  // If the node is not hanging
                  else
                  {
                    // The local unknown number comes from the node
                    local_unknown = this->nodal_local_eqn(l2, u_nodal_index);
                    // The hang weight is one
                    hang_weight2 = 1.0;
                  }

                  // If the unknown is not pinned
                  if (local_unknown >= 0)
                  {
                    // Add contribution to Elemental Matrix
                    // Mass matrix du/dt term
                    jacobian(local_eqn, local_unknown) -=
                      peclet_st * test(l) * psi(l2) *
                      this->node_pt(l2)->time_stepper_pt()->weight(1, 0) * W *
                      hang_weight * hang_weight2;

                    // Add contribution to mass matrix
                    if (flag == 2)
                    {
                      mass_matrix(local_eqn, local_unknown) +=
                        peclet_st * test(l) * psi(l2) * W * hang_weight *
                        hang_weight2;
                    }

                    // Add contribution to Elemental Matrix
                    for (unsigned k = 0; k < DIM; k++)
                    {
                      // Temporary term used in assembly
                      double tmp = peclet * wind[k];
                      if (!ALE_is_disabled_flag)
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
                        (tmp * test(l) + tmp2 * dtestdx(l, k)) * W *
                        hang_weight * hang_weight2;
                    }
                  }
                } // End of loop over master nodes
              } // End of loop over nodes
            } // End of Jacobian calculation

          } // End of non-zero equation

        } // End of loop over the master nodes for residual
      } // End of loop over nodes

    } // End of loop over integration points
  }


  //====================================================================
  // Force build of templates
  //====================================================================
  template class RefineableQGeneralisedAdvectionDiffusionElement<2, 2>;
  template class RefineableQGeneralisedAdvectionDiffusionElement<2, 3>;
  template class RefineableQGeneralisedAdvectionDiffusionElement<2, 4>;

  template class RefineableQGeneralisedAdvectionDiffusionElement<3, 2>;
  template class RefineableQGeneralisedAdvectionDiffusionElement<3, 3>;
  template class RefineableQGeneralisedAdvectionDiffusionElement<3, 4>;

} // namespace oomph
