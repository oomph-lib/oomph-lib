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
#include "refineable_young_laplace_elements.h"


namespace oomph
{
  //======================================================================
  /// Compute element residual vector taking hanging nodes into account
  //======================================================================
  void RefineableYoungLaplaceEquations::fill_in_contribution_to_residuals(
    Vector<double>& residuals)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Set up memory for the shape (test) functions
    Shape psi(n_node);
    DShape dpsidzeta(n_node, 2);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Integers to store the local equation numbers
    int local_eqn = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape (test) functions
      double J = dshape_eulerian_at_knot(ipt, psi, dpsidzeta);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the pressure and velocity components
      // Allocate and initialise to zero
      double interpolated_u = 0.0;
      Vector<double> interpolated_zeta(2, 0.0);
      Vector<double> interpolated_du_dzeta(2, 0.0);

      // Calculate function value and derivatives:
      //-----------------------------------------
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += u(l) * psi(l);
        // Loop over directions
        for (unsigned j = 0; j < 2; j++)
        {
          interpolated_zeta[j] += nodal_position(l, j) * psi(l);
          interpolated_du_dzeta[j] += u(l) * dpsidzeta(l, j);
        }
      }

      // Allocation and definition of variables necessary for
      // further calculations

      /// "Simple" case
      /// --------------
      double nonlinearterm = 1.0;
      double sqnorm = 0.0;

      /// Spine case
      /// -----------

      // Derivs of position vector w.r.t. global intrinsic coords
      Vector<Vector<double>> dRdzeta;
      allocate_vector_of_vectors(2, 3, dRdzeta);

      // Unnormalised normal
      Vector<double> N_unnormalised(3, 0.0);

      // Spine and spine basis vectors, entries initialised to zero
      Vector<double> spine_base(3, 0.0), spine(3, 0.0);

      // Derivative of spine basis vector w.r.t to the intrinsic
      // coordinates: dspine_base[i,j] = j-th component of the deriv.
      // of the spine basis vector w.r.t. to the i-th global intrinsic
      // coordinate
      Vector<Vector<double>> dspine_base;
      allocate_vector_of_vectors(2, 3, dspine_base);

      // Derivative of spine vector w.r.t to the intrinsic
      // coordinates: dspine[i,j] = j-th component of the deriv.
      // of the spine vector w.r.t. to the i-th global intrinsic
      // coordinate
      Vector<Vector<double>> dspine;
      allocate_vector_of_vectors(2, 3, dspine);

      // Vector v_\alpha contains the numerator of the variations of the
      // area element {\cal A}^{1/2} w.r.t. the components of dR/d\zeta_\alpha
      Vector<double> area_variation_numerator_0(3, 0.0);
      Vector<double> area_variation_numerator_1(3, 0.0);

      // Vector position
      Vector<double> r(3, 0.0);

      // No spines
      //---------
      if (!use_spines())
      {
        for (unsigned j = 0; j < 2; j++)
          sqnorm += interpolated_du_dzeta[j] * interpolated_du_dzeta[j];

        nonlinearterm = 1.0 / sqrt(1.0 + sqnorm);
      }

      // Spines
      //------
      else
      {
        // Get the spines
        get_spine_base(interpolated_zeta, spine_base, dspine_base);
        get_spine(interpolated_zeta, spine, dspine);

        // calculation of dR/d\zeta_\alpha
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          // Product rule for d(u {\bf S} ) / d \zeta_\alpha
          Vector<double> dudzeta_times_spine(3, 0.0);
          scalar_times_vector(
            interpolated_du_dzeta[alpha], spine, dudzeta_times_spine);

          Vector<double> u_times_dspinedzeta(3, 0.0);
          scalar_times_vector(
            interpolated_u, dspine[alpha], u_times_dspinedzeta);

          Vector<double> d_u_times_spine_dzeta(3, 0.0);
          vector_sum(
            dudzeta_times_spine, u_times_dspinedzeta, d_u_times_spine_dzeta);

          // Add derivative of spine base
          vector_sum(d_u_times_spine_dzeta, dspine_base[alpha], dRdzeta[alpha]);
        }

        /// Get the unnormalized normal
        cross_product(dRdzeta[0], dRdzeta[1], N_unnormalised);

        /// Tmp storage
        Vector<double> v_tmp_1(3, 0.0);
        Vector<double> v_tmp_2(3, 0.0);

        // Calculation of
        // |dR/d\zeta_1|^2 dR/d\zeta_0 - <dR/d\zeta_0,dR/d\zeta_1>dR/d\zeta_1
        scalar_times_vector(pow(two_norm(dRdzeta[1]), 2), dRdzeta[0], v_tmp_1);
        scalar_times_vector(
          -1 * scalar_product(dRdzeta[0], dRdzeta[1]), dRdzeta[1], v_tmp_2);
        vector_sum(v_tmp_1, v_tmp_2, area_variation_numerator_0);

        // Calculation of
        // |dR/d\zeta_0|^2 dR/d\zeta_1 - <dR/d\zeta_0,dR/d\zeta_1>dR/d\zeta_0
        scalar_times_vector(pow(two_norm(dRdzeta[0]), 2), dRdzeta[1], v_tmp_1);
        scalar_times_vector(
          -1 * scalar_product(dRdzeta[0], dRdzeta[1]), dRdzeta[0], v_tmp_2);
        vector_sum(v_tmp_1, v_tmp_2, area_variation_numerator_1);

        // Global Eulerian cooordinates
        for (unsigned j = 0; j < 3; j++)
        {
          r[j] = spine_base[j] + interpolated_u * spine[j];
        }
      }


      // Assemble residuals
      //-------------------

      // Loop over the nodes for the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Local variables used to store the number of master nodes and the
        // weight associated with the shape function if the node is hanging
        unsigned n_master = 1;
        double hang_weight = 1.0;

        // If the node is hanging, get the number of master nodes
        if (node_pt(l)->is_hanging())
        {
          n_master = node_pt(l)->hanging_pt()->nmaster();
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
          if (node_pt(l)->is_hanging())
          {
            // Read out the local equation from the master node
            local_eqn =
              local_hang_eqn(node_pt(l)->hanging_pt()->master_node_pt(m), 0);
            // Read out the weight from the master node
            hang_weight = node_pt(l)->hanging_pt()->master_weight(m);
          }
          // If the node is not hanging
          else
          {
            // The local equation number comes from the node itself
            local_eqn = this->u_local_eqn(l);
            // The hang weight is one
            hang_weight = 1.0;
          }

          /*IF it's not a boundary condition*/
          if (local_eqn >= 0)
          {
            // "simple" calculation case
            if (!use_spines())
            {
              // Add source term: The curvature
              residuals[local_eqn] += get_kappa() * psi(l) * W * hang_weight;

              // The YoungLaplace bit itself
              for (unsigned k = 0; k < 2; k++)
              {
                residuals[local_eqn] += nonlinearterm *
                                        interpolated_du_dzeta[k] *
                                        dpsidzeta(l, k) * W * hang_weight;
              }
            }

            // Spine calculation case
            else
            {
              // Calculation of d(u S)/d\zeta_0
              //-------------------------------
              Vector<double> v_tmp_1(3, 0.0);
              scalar_times_vector(dpsidzeta(l, 0), spine, v_tmp_1);

              Vector<double> v_tmp_2(3, 0.0);
              scalar_times_vector(psi(l), dspine[0], v_tmp_2);

              Vector<double> d_uS_dzeta0(3, 0.0);
              vector_sum(v_tmp_1, v_tmp_2, d_uS_dzeta0);

              // Add contribution to residual
              residuals[local_eqn] +=
                W * hang_weight *
                scalar_product(area_variation_numerator_0, d_uS_dzeta0) /
                two_norm(N_unnormalised);

              // Calculation of d(u S)/d\zeta_1
              scalar_times_vector(dpsidzeta(l, 1), spine, v_tmp_1);
              scalar_times_vector(psi(l), dspine[1], v_tmp_2);
              Vector<double> d_uS_dzeta1(3, 0.0);
              vector_sum(v_tmp_1, v_tmp_2, d_uS_dzeta1);

              // Add contribution to residual
              residuals[local_eqn] +=
                W * hang_weight *
                scalar_product(area_variation_numerator_1, d_uS_dzeta1) /
                two_norm(N_unnormalised);

              // Curvature contribution to the residual  : kappa N S test
              residuals[local_eqn] += W * hang_weight * (get_kappa()) *
                                      scalar_product(N_unnormalised, spine) *
                                      psi(l);
            }
          }

        } // End of loop over master nodes for residual

      } // End of loop over nodes

    } // End of loop over integration points
  }

  //====================================================================
  // Force build of templates
  //====================================================================
  template class RefineableQYoungLaplaceElement<2>;
  template class RefineableQYoungLaplaceElement<3>;
  template class RefineableQYoungLaplaceElement<4>;

} // namespace oomph
