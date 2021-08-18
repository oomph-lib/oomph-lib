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
#include "refineable_helmholtz_elements.h"


namespace oomph
{
  //========================================================================
  /// Add element's contribution to the elemental
  /// residual vector and/or Jacobian matrix.
  /// flag=1: compute both
  /// flag=0: compute only residual vector
  //========================================================================
  template<unsigned DIM>
  void RefineableHelmholtzEquations<DIM>::
    fill_in_generic_residual_contribution_helmholtz(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
  {
    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Integers to store the local equation and unknown numbers
    int local_eqn_real = 0, local_unknown_real = 0;
    int local_eqn_imag = 0, local_unknown_imag = 0;

    // Local storage for pointers to hang_info objects
    HangInfo *hang_info_pt = 0, *hang_info2_pt = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = this->dshape_and_dtest_eulerian_at_knot_helmholtz(
        ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Position and gradient
      std::complex<double> interpolated_u(0.0, 0.0);
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<std::complex<double>> interpolated_dudx(DIM);

      // Calculate function value and derivatives:
      //-----------------------------------------

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_x[j] += nodal_position(l, j) * psi(l);
        }
        // Get the nodal value of the helmholtz unknown
        const std::complex<double> u_value(
          nodal_value(l, this->u_index_helmholtz().real()),
          nodal_value(l, this->u_index_helmholtz().imag()));
        // Add to the interpolated value
        interpolated_u += u_value * psi(l);

        // Loop over directions
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_dudx[j] += u_value * dpsidx(l, j);
        }
      }

      // Get body force
      std::complex<double> source(0.0, 0.0);
      this->get_source_helmholtz(ipt, interpolated_x, source);


      // Assemble residuals and Jacobian

      // Loop over the nodes for the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Local variables used to store the number of master nodes and the
        // weight associated with the shape function if the node is hanging
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

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Get the local equation number and hang_weight
          // If the node is hanging
          if (is_node_hanging)
          {
            // Read out the local equation number from the m-th master node
            local_eqn_real =
              this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                   this->u_index_helmholtz().real());
            local_eqn_imag =
              this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                   this->u_index_helmholtz().imag());

            // Read out the weight from the master node
            hang_weight = hang_info_pt->master_weight(m);
          }
          // If the node is not hanging
          else
          {
            // The local equation number comes from the node itself
            local_eqn_real =
              this->nodal_local_eqn(l, this->u_index_helmholtz().real());
            local_eqn_imag =
              this->nodal_local_eqn(l, this->u_index_helmholtz().imag());

            // The hang weight is one
            hang_weight = 1.0;
          }

          // If the nodal equation is not a boundary condition
          if (local_eqn_real >= 0)
          {
            // Add body force/source term here and Helmholtz bit
            residuals[local_eqn_real] +=
              (source.real() - this->k_squared() * interpolated_u.real()) *
              test(l) * W * hang_weight;

            // The Helmholtz bit itself
            for (unsigned k = 0; k < DIM; k++)
            {
              residuals[local_eqn_real] +=
                interpolated_dudx[k].real() * dtestdx(l, k) * W * hang_weight;
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
                    local_unknown_real =
                      this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                           this->u_index_helmholtz().real());

                    // Read out the hanging weight from the master node
                    hang_weight2 = hang_info2_pt->master_weight(m2);
                  }
                  // If the node is not hanging
                  else
                  {
                    // The local unknown number comes from the node
                    local_unknown_real = this->nodal_local_eqn(
                      l2, this->u_index_helmholtz().real());

                    // The hang weight is one
                    hang_weight2 = 1.0;
                  }

                  // If the unknown is not pinned
                  if (local_unknown_real >= 0)
                  {
                    // Add contribution to Elemental Matrix
                    for (unsigned i = 0; i < DIM; i++)
                    {
                      jacobian(local_eqn_real, local_unknown_real) +=
                        dpsidx(l2, i) * dtestdx(l, i) * W * hang_weight *
                        hang_weight2;
                    }

                    // Add the helmholtz contribution
                    jacobian(local_eqn_real, local_unknown_real) -=
                      this->k_squared() * psi(l2) * test(l) * W * hang_weight *
                      hang_weight2;
                  }
                } // End of loop over master nodes
              } // End of loop over nodes
            } // End of Jacobian calculation
          } // End of case when residual equation is not pinned


          if (local_eqn_imag >= 0)
          {
            // Add body force/source term here and Helmholtz bit
            residuals[local_eqn_imag] +=
              (source.imag() - this->k_squared() * interpolated_u.imag()) *
              test(l) * W * hang_weight;

            // The Helmholtz bit itself
            for (unsigned k = 0; k < DIM; k++)
            {
              residuals[local_eqn_imag] +=
                interpolated_dudx[k].imag() * dtestdx(l, k) * W * hang_weight;
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
                    local_unknown_imag =
                      this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                           this->u_index_helmholtz().imag());

                    // Read out the hanging weight from the master node
                    hang_weight2 = hang_info2_pt->master_weight(m2);
                  }
                  // If the node is not hanging
                  else
                  {
                    // The local unknown number comes from the node

                    local_unknown_imag = this->nodal_local_eqn(
                      l2, this->u_index_helmholtz().imag());

                    // The hang weight is one
                    hang_weight2 = 1.0;
                  }

                  if (local_unknown_imag >= 0)
                  {
                    // Add contribution to Elemental Matrix
                    for (unsigned i = 0; i < DIM; i++)
                    {
                      jacobian(local_eqn_imag, local_unknown_imag) +=
                        dpsidx(l2, i) * dtestdx(l, i) * W * hang_weight *
                        hang_weight2;
                    }

                    // Add the helmholtz contribution
                    jacobian(local_eqn_imag, local_unknown_imag) -=
                      this->k_squared() * psi(l2) * test(l) * W * hang_weight *
                      hang_weight2;
                  }
                } // End of loop over master nodes
              } // End of loop over nodes
            } // End of Jacobian calculation
          }
        } // End of loop over master nodes for residual
      } // End of loop over nodes

    } // End of loop over integration points
  }


  //====================================================================
  // Force build of templates
  //====================================================================
  template class RefineableQHelmholtzElement<1, 2>;
  template class RefineableQHelmholtzElement<1, 3>;
  template class RefineableQHelmholtzElement<1, 4>;

  template class RefineableQHelmholtzElement<2, 2>;
  template class RefineableQHelmholtzElement<2, 3>;
  template class RefineableQHelmholtzElement<2, 4>;

  template class RefineableQHelmholtzElement<3, 2>;
  template class RefineableQHelmholtzElement<3, 3>;
  template class RefineableQHelmholtzElement<3, 4>;

} // namespace oomph
