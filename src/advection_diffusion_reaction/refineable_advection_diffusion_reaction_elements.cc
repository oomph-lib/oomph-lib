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
#include "refineable_advection_diffusion_reaction_elements.h"

namespace oomph
{
  //==========================================================================
  /// Add  the element's contribution to the elemental residual vector
  /// and/or elemental jacobian matrix.
  /// This function overloads the standard version so that the possible
  /// presence of hanging nodes is taken into account.
  //=========================================================================
  template<unsigned NREAGENT, unsigned DIM>
  void RefineableAdvectionDiffusionReactionEquations<NREAGENT, DIM>::
    fill_in_generic_residual_contribution_adv_diff_react(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Find out how many nodes there are in the element
    const unsigned n_node = nnode();

    // Get the nodal index at which the unknown is stored
    unsigned c_nodal_index[NREAGENT];
    for (unsigned r = 0; r < NREAGENT; r++)
    {
      c_nodal_index[r] = this->c_index_adv_diff_react(r);
    }

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get diffusion coefficients
    Vector<double> D = this->diff();

    // Get the timescales
    Vector<double> T = this->tau();

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
      double J = this->dshape_and_dtest_eulerian_at_knot_adv_diff_react(
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
          interpolated_x[j] += nodal_position(l, j) * psi(l);
        }

        // Loop over the unknown reagents
        for (unsigned r = 0; r < NREAGENT; r++)
        {
          // Get the value at the node
          const double c_value = nodal_value(l, c_nodal_index[r]);

          // Calculate the interpolated value
          interpolated_c[r] += c_value * psi(l);
          dcdt[r] += this->dc_dt_adv_diff_react(l, r) * psi(l);

          // Loop over directions to calculate the derivatie
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_dcdx(r, j) += c_value * dpsidx(l, j);
          }
        }
      }

      // Mesh velocity?
      if (!ALE_is_disabled_flag)
      {
        for (unsigned l = 0; l < n_node; l++)
        {
          for (unsigned j = 0; j < DIM; j++)
          {
            mesh_velocity[j] += dnodal_position_dt(l, j) * psi(l);
          }
        }
      }

      // Get source function
      Vector<double> source(NREAGENT);
      this->get_source_adv_diff_react(ipt, interpolated_x, source);


      // Get wind
      Vector<double> wind(DIM);
      this->get_wind_adv_diff_react(ipt, s, interpolated_x, wind);

      // Get reaction terms
      Vector<double> R(NREAGENT);
      this->get_reaction_adv_diff_react(ipt, interpolated_c, R);

      // If we are getting the jacobian the get the derivative terms
      DenseMatrix<double> dRdC(NREAGENT);
      if (flag)
      {
        this->get_reaction_deriv_adv_diff_react(ipt, interpolated_c, dRdC);
      }

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
          // Loop over the number of reagents
          for (unsigned r = 0; r < NREAGENT; r++)
          {
            // Get the local equation number and hang_weight
            // If the node is hanging
            if (is_node_hanging)
            {
              // Read out the local equation from the master node
              local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                               c_nodal_index[r]);
              // Read out the weight from the master node
              hang_weight = hang_info_pt->master_weight(m);
            }
            // If the node is not hanging
            else
            {
              // The local equation number comes from the node itself
              local_eqn = this->nodal_local_eqn(l, c_nodal_index[r]);
              // The hang weight is one
              hang_weight = 1.0;
            }

            // If the nodal equation is not a boundary conditino
            if (local_eqn >= 0)
            {
              // Add body force/source/reaction term and time derivative
              residuals[local_eqn] -=
                (T[r] * dcdt[r] + source[r] + R[r]) * test(l) * W * hang_weight;

              // The Advection Diffusion bit itself
              for (unsigned k = 0; k < DIM; k++)
              {
                // Terms that multiply the test function
                double tmp = wind[k];
                // If the mesh is moving need to subtract the mesh velocity
                if (!ALE_is_disabled_flag)
                {
                  tmp -= T[r] * mesh_velocity[k];
                }
                // Now construct the contribution to the residuals
                residuals[local_eqn] -= interpolated_dcdx(r, k) *
                                        (tmp * test(l) + D[r] * dtestdx(l, k)) *
                                        W * hang_weight;
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
                    // Loop over the reagents again
                    for (unsigned r2 = 0; r2 < NREAGENT; r2++)
                    {
                      // Get the local unknown and weight
                      // If the node is hanging
                      if (is_node2_hanging)
                      {
                        // Read out the local unknown from the master node
                        local_unknown = this->local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), c_nodal_index[r2]);
                        // Read out the hanging weight from the master node
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      // If the node is not hanging
                      else
                      {
                        // The local unknown number comes from the node
                        local_unknown =
                          this->nodal_local_eqn(l2, c_nodal_index[r2]);
                        // The hang weight is one
                        hang_weight2 = 1.0;
                      }

                      // If the unknown is not pinned
                      if (local_unknown >= 0)
                      {
                        // Diagonal terms (i.e. the basic equations)
                        if (r2 == r)
                        {
                          // Mass matrix term
                          jacobian(local_eqn, local_unknown) -=
                            T[r] * test(l) * psi(l2) *
                            node_pt(l2)->time_stepper_pt()->weight(1, 0) * W *
                            hang_weight * hang_weight2;

                          // Add the mass matrix term
                          if (flag == 2)
                          {
                            mass_matrix(local_eqn, local_unknown) +=
                              T[r] * test(l) * psi(l2) * W * hang_weight *
                              hang_weight2;
                          }

                          // Add contribution to Elemental Matrix
                          for (unsigned i = 0; i < DIM; i++)
                          {
                            // Temporary term used in assembly
                            double tmp = wind[i];
                            if (!ALE_is_disabled_flag)
                            {
                              tmp -= T[r] * mesh_velocity[i];
                            }
                            // Now assemble Jacobian term
                            jacobian(local_eqn, local_unknown) -=
                              dpsidx(l2, i) *
                              (tmp * test(l) + D[r] * dtestdx(l, i)) * W *
                              hang_weight * hang_weight2;
                          }

                        } // End of diagonal terms

                        // Now add the cross-reaction terms
                        jacobian(local_eqn, local_unknown) -=
                          dRdC(r, r2) * psi(l2) * test(l) * W * hang_weight *
                          hang_weight2;
                      }
                    } // End of loop over reagents
                  } // End of loop over master nodes
                } // End of loop over nodes
              } // End of Jacobian calculation

            } // End of non-zero equation

          } // End of loop over reagents
        } // End of loop over the master nodes for residual
      } // End of loop over nodes

    } // End of loop over integration points
  }


  //====================================================================
  // Force build of templates
  //====================================================================
  /// One reagent
  template class RefineableQAdvectionDiffusionReactionElement<1, 2, 2>;
  template class RefineableQAdvectionDiffusionReactionElement<1, 2, 3>;
  template class RefineableQAdvectionDiffusionReactionElement<1, 2, 4>;

  template class RefineableQAdvectionDiffusionReactionElement<1, 3, 2>;
  template class RefineableQAdvectionDiffusionReactionElement<1, 3, 3>;
  template class RefineableQAdvectionDiffusionReactionElement<1, 3, 4>;

  // Two reagents
  template class RefineableQAdvectionDiffusionReactionElement<2, 1, 2>;
  template class RefineableQAdvectionDiffusionReactionElement<2, 1, 3>;
  template class RefineableQAdvectionDiffusionReactionElement<2, 1, 4>;

  template class RefineableQAdvectionDiffusionReactionElement<2, 2, 2>;
  template class RefineableQAdvectionDiffusionReactionElement<2, 2, 3>;
  template class RefineableQAdvectionDiffusionReactionElement<2, 2, 4>;

  template class RefineableQAdvectionDiffusionReactionElement<2, 3, 2>;
  template class RefineableQAdvectionDiffusionReactionElement<2, 3, 3>;
  template class RefineableQAdvectionDiffusionReactionElement<2, 3, 4>;

} // namespace oomph
