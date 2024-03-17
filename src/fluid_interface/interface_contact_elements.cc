// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
// Non-inline functions for fluid free surface elements

// OOMPH-LIB headers
#include "interface_contact_elements.h"


namespace oomph
{
  //=========================================================================
  /// Add contribution to element's residual vector and Jacobian
  //=========================================================================
  void ContactLineElement::
    fill_in_generic_residual_contribution_interface_boundary(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
    // Let's get the info from the parent
    FaceElement* parent_pt = dynamic_cast<FaceElement*>(bulk_element_pt());

    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Set up memory for the shape functions
    Shape psi(n_node);
    DShape dpsids(n_node, dim());
    Vector<double> s_local(dim());

    // Find the dimension of the problem
    unsigned spatial_dim = this->nodal_dimension();

    // Set up the identity matrix
    DenseMatrix<double> identity(spatial_dim, spatial_dim, 0.0);
    for (unsigned i = 0; i < spatial_dim; i++)
    {
      identity(i, i) = 1.0;
    }

    // Storage for the coordinate
    Vector<double> x(spatial_dim);
    Vector<double> dx_dt(spatial_dim);
    Vector<double> u(spatial_dim);

    // Outer unit normal to the wall and surface
    Vector<double> wall_normal(spatial_dim);
    Vector<double> surface_normal(spatial_dim);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Jacobian of mapping
      double J = 1.0;
      if (dim() >= 1)
      {
        J = J_eulerian_at_knot(ipt);
      }


      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate the shape functions
      shape_at_knot(ipt, psi);

      // Find the dimension of the element
      const unsigned el_dim = dim();
      // Storage for the local coordinates of the integration point
      Vector<double> s_local(el_dim);
      // Set the local coordinate
      for (unsigned i = 0; i < el_dim; i++)
      {
        s_local[i] = integral_pt()->knot(ipt, i);
      }

      // Get the x coordinate
      this->interpolated_x(s_local, x);

      // Get the dx/dt of the coordinate
      const unsigned t_deriv = 1;
      this->interpolated_dxdt(s_local, t_deriv, dx_dt);

      // Assemble interpolated values
      // Loop over the nodes
      double lambda = 0.0;
      for (unsigned n = 0; n < n_node; n++)
      {
        for (unsigned i = 0; i < spatial_dim; i++)
        {
          u[i] += this->nst_u(n, i) * psi(n);
        }
        lambda += fsi_lagrange_multiplier(n) * psi(n);
      }
      // Get the unit normal to the wall
      wall_unit_normal(x, wall_normal);

      // Find the local coordinates in the parent
      // Find the dimension of the parent
      unsigned n_dim = parent_pt->dim();
      Vector<double> s_parent(n_dim);
      this->get_local_coordinate_in_bulk(s_local, s_parent);

      // Just get the outer unit normal
      parent_pt->outer_unit_normal(s_parent, surface_normal);

      double contact_angle_local = get_contact_angle(x, dx_dt);

      // Get the value of sigma from the parent
      double sigma_local = sigma();

      // Get the capillary number
      double ca_local = ca();

      // Get the wall tangent vector
      Vector<double> wall_tangent(spatial_dim);
      // Rotate the normal vector clockwise 90 degrees
      wall_tangent[0] = wall_normal[1];
      wall_tangent[1] = -wall_normal[0];

      Vector<double> project_to_interface_normal_on_wall(spatial_dim, 0.0);
      for (unsigned d = 0; d < spatial_dim; d++)
      {
        for (unsigned i = 0; i < spatial_dim; i++)
        {
          project_to_interface_normal_on_wall[d] +=
            (identity(d, i) - wall_normal[d] * wall_normal[i]) *
            surface_normal[d];
        }
      }

      // Loop over the nodes
      for (unsigned n = 0; n < n_node; n++)
      {
        // Add the tension component to the momentum equations
        // Used in both the strong and weak forms
        for (unsigned i = 0; i < spatial_dim; i++)
        {
          int local_eqn = nst_momentum_local_eqn(n, i);
          if (local_eqn >= 0)
          {
            residuals[local_eqn] += (sigma_local / ca_local) *
                                    (sin(contact_angle_local) * wall_normal[i] +
                                     cos(contact_angle_local) *
                                       project_to_interface_normal_on_wall[i]);
          }
        }

        // STRONG contact angle imposition
        if (Contact_angle_flag[n] == STRONG)
        {
          /// Strong form of the contact angle constraint by hijacking of the
          /// free surface kinematic equation

          // Read out the kinematic equation number
          int local_eqn = this->fsi_kinematic_local_eqn(n);

          // Note that because we have outer unit normals for the free surface
          // and the wall, the cosine of the contact angle is equal to
          // MINUS the dot product computed above
          if (local_eqn >= 0)
          {
            // Find the dot product of the two vectors
            double dot = 0.0;
            for (unsigned i = 0; i < spatial_dim; i++)
            {
              dot += surface_normal[i] * wall_normal[i];
            }

            // Reset the residuals, as we are hijacking.
            // TODO Check that we need to do this.
            if (n == 0)
            {
              residuals[local_eqn] = 0.0;
            }

            residuals[local_eqn] +=
              (cos(contact_angle_local) + dot) * psi(n) * W;
          }
        }
        // WEAK contact angle imposition
        else if (Contact_angle_flag[n] == WEAK)
        {
          //  int local_eqn = wall_bounded_kinematic_local_eqn(n);
          //  std::cout << "local_eqn: " << local_eqn << std::endl;
          //  std::cout << "mom 0 eqn: " << nst_momentum_local_eqn(0, 0)
          //            << std::endl;
          //  std::cout << "mom 1 eqn: " << nst_momentum_local_eqn(0, 1)
          //            << std::endl;
          //  std::cout << "mom 2 eqn: " << nst_momentum_local_eqn(0, 2)
          //            << std::endl;
          //  std::cout << "cont eqn: "
          //            << nodal_local_eqn(0, nst_continuity_index(0)) <<
          //            std::endl;
          //  std::cout << "kin eqn: " << fsi_kinematic_local_eqn(0) <<
          //  std::endl; std::cout << "nvalue: " << node_pt(n)->nvalue() <<
          //  std::endl;
          //  // std::cout << "nodal index: " <<
          //  // cl_lagrange_multiplier_nodal_index(n)
          //  //          << std::endl;
          //  if (local_eqn >= 0)
          //  {
          //    // std::cout << "Add to wall bounded kinematic equations, "
          //    //          << cl_lagrange_multiplier_nodal_index(n) <<
          //    std::endl;
          //    // Point kinematic equation
          //    double st_local = st();

          //    // Wall-bounded kinematic equation
          //    for (unsigned k = 0; k < spatial_dim; k++)
          //    {
          //      residuals[local_eqn] += (u[k] - st_local * dx_dt[k]) *
          //                              project_to_interface_normal_on_wall[k]
          //                              * psi(n) * W;

          //      // Lagrange multiplier contribution to momentum equations
          //      local_eqn = this->nst_momentum_local_eqn(n, k);
          //      if (local_eqn >= 0)
          //      {
          //        residuals[local_eqn] +=
          //          lambda * project_to_interface_normal_on_wall[k] * psi(n) *
          //          W;
          //      }

          //      // Lagrange multiplier contribution to solid displacement
          //      // equations
          //      local_eqn = this->fsi_kinematic_local_eqn(n);
          //      if (local_eqn >= 0)
          //      {
          //        double residual_contribution =
          //          -lambda * st_local *
          //          project_to_interface_normal_on_wall[k] * psi(n) * W;
          //        // TODO Check if we need further contributions from
          //        // project_to_interface_normal_on_wall
          //        // If we are timestepping,...
          //        const unsigned time_derivative = 1;
          //        const unsigned weight_number = 0;
          //        const double time_stepper_weight =
          //          this->node_pt(n)->time_stepper_pt()->weight(time_derivative,
          //                                                      weight_number);
          //        if (time_stepper_weight > 0)
          //        {
          //          // ... scale the contribution by the weight.
          //          residual_contribution *= time_stepper_weight;
          //        }

          //        residuals[local_eqn] += residual_contribution;
          //      }
          //    }
          //  }
        }
        // NOTE: The jacobian entries will be computed automatically
        // by finite differences.

      } // Node loop

      // Now add the additional contributions
      add_additional_residual_contributions_interface_boundary(
        residuals, jacobian, flag, psi, dpsids, surface_normal, W);
    } // Integration points loop
    // std::cout << "Out" << std::endl;
  }
} // namespace oomph
