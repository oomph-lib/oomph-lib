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
// Non-inline functions for fluid free surface elements

// OOMPH-LIB headers
#include "surfactant_transport_elements.h"

namespace oomph
{
  // Define the default physical value to be one
  double SurfactantTransportInterfaceElement::Default_Physical_Constant_Value =
    1.0;

  //=====================================================================
  /// Get the surfactant concentration
  //=====================================================================
  double SurfactantTransportInterfaceElement::interpolated_C(
    const Vector<double>& s)
  {
    // Find number of nodes
    unsigned n_node = this->nnode();

    // Local shape function
    Shape psi(n_node);

    // Find values of shape function
    this->shape(s, psi);

    // Initialise value of C
    double C = 0.0;

    // Loop over the local nodes and sum
    for (unsigned l = 0; l < n_node; l++)
    {
      C += this->nodal_value(l, this->C_index[l]) * psi(l);
    }

    return (C);
  }

  //=====================================================================
  /// The time derivative of the surface concentration
  //=====================================================================
  double SurfactantTransportInterfaceElement::dcdt_surface(
    const unsigned& l) const
  {
    // Get the data's timestepper
    TimeStepper* time_stepper_pt = this->node_pt(l)->time_stepper_pt();

    // Initialise dudt
    double dcdt = 0.0;
    // Loop over the timesteps, if there is a non Steady timestepper
    if (time_stepper_pt->type() != "Steady")
    {
      // Number of timsteps (past & present)
      const unsigned n_time = time_stepper_pt->ntstorage();

      for (unsigned t = 0; t < n_time; t++)
      {
        dcdt += time_stepper_pt->weight(1, t) *
                this->nodal_value(t, l, this->C_index[l]);
      }
    }
    return dcdt;
  }

  //=====================================================================
  /// The surface tension function is linear in the
  /// concentration with constant of proportionality equal
  /// to the elasticity  number.
  //=====================================================================
  double SurfactantTransportInterfaceElement::sigma(const Vector<double>& s)
  {
    // Find the number of shape functions
    const unsigned n_node = this->nnode();
    // Now get the shape fuctions at the local coordinate
    Shape psi(n_node);
    this->shape(s, psi);

    // Now interpolate the temperature and surfactant concentration
    double C = 0.0;
    for (unsigned l = 0; l < n_node; l++)
    {
      C += this->nodal_value(l, this->C_index[l]) * psi(l);
    }

    // Get the Elasticity numbers
    double Beta = this->beta();
    // Return the variable surface tension
    return (1.0 - Beta * (C - 1.0));
  }

  //=======================================================================
  /// \short Overload the Helper function to calculate the residuals and
  /// jacobian entries. This particular function ensures that the
  /// additional entries are calculated inside the integration loop
  void SurfactantTransportInterfaceElement::
    add_additional_residual_contributions_interface(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag,
      const Shape& psif,
      const DShape& dpsifds,
      const DShape& dpsifdS,
      const DShape& dpsifdS_div,
      const Vector<double>& s,
      const Vector<double>& interpolated_x,
      const Vector<double>& interpolated_n,
      const double& W,
      const double& J)
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Storage for the local equation numbers and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Surface advection-diffusion equation

    // Find the index at which the concentration is stored
    Vector<unsigned> u_index = this->U_index_interface;

    // Read out the surface peclect number
    const double Pe_s = this->peclet_s();
    // const double PeSt_s = this->peclet_strouhal_s();

    // Now calculate the concentration and derivatives at this point
    // Assuming the same shape functions are used (which they are)
    double interpolated_C = 0.0;
    double dCdt = 0.0;
    // The tangent vectors and velocity
    const unsigned n_dim = this->node_pt(0)->ndim();
    Vector<double> interpolated_u(n_dim, 0.0);
    Vector<double> mesh_velocity(n_dim, 0.0);
    Vector<double> interpolated_grad_C(n_dim, 0.0);
    double interpolated_div_u = 0.0;

    // Loop over the shape functions
    for (unsigned l = 0; l < n_node; l++)
    {
      const double psi = psif(l);
      const double C_ = this->nodal_value(l, this->C_index[l]);

      interpolated_C += C_ * psi;
      dCdt += dcdt_surface(l) * psi;
      // Velocity and Mesh Velocity
      for (unsigned j = 0; j < n_dim; j++)
      {
        const double u_ = this->nodal_value(l, u_index[j]);
        interpolated_u[j] += u_ * psi;
        mesh_velocity[j] += this->dnodal_position_dt(l, j) * psi;
        interpolated_grad_C[j] += C_ * dpsifdS(l, j);
        interpolated_div_u += u_ * dpsifdS_div(l, j);
      }
    }

    // Pre-compute advection term
    double interpolated_advection_term = interpolated_C * interpolated_div_u;
    for (unsigned i = 0; i < n_dim; i++)
    {
      interpolated_advection_term +=
        (interpolated_u[i] - mesh_velocity[i]) * interpolated_grad_C[i];
    }


    // Now we add the flux term to the appropriate residuals
    for (unsigned l = 0; l < n_node; l++)
    {
      // Read out the apprporiate local equation
      local_eqn = this->nodal_local_eqn(l, this->C_index[l]);

      // If not a boundary condition
      if (local_eqn >= 0)
      {
        // Time derivative term
        residuals[local_eqn] += dCdt * psif(l) * W * J;

        // First Advection term
        residuals[local_eqn] += interpolated_advection_term * psif(l) * W * J;

        // Diffusion term
        double diffusion_term = 0.0;
        for (unsigned i = 0; i < n_dim; i++)
        {
          diffusion_term += interpolated_grad_C[i] * dpsifdS(l, i);
        }
        residuals[local_eqn] += (1.0 / Pe_s) * diffusion_term * W * J;

        // We also need to worry about the jacobian terms
        if (flag)
        {
          // Loop over the nodes again
          for (unsigned l2 = 0; l2 < n_node; l2++)
          {
            // Get the time stepper
            TimeStepper* time_stepper_pt = this->node_pt(l2)->time_stepper_pt();

            // Get the unknown c_index
            local_unknown = this->nodal_local_eqn(l2, this->C_index[l2]);

            if (local_unknown >= 0)
            {
              jacobian(local_eqn, local_unknown) +=
                time_stepper_pt->weight(1, 0) * psif(l2) * psif(l) * W * J;

              jacobian(local_eqn, local_unknown) +=
                psif(l2) * interpolated_div_u * psif(l) * W * J;

              for (unsigned i = 0; i < n_dim; i++)
              {
                jacobian(local_eqn, local_unknown) +=
                  (interpolated_u[i] - mesh_velocity[i]) * dpsifdS(l2, i) *
                  psif(l) * W * J;
              }

              for (unsigned i = 0; i < n_dim; i++)
              {
                jacobian(local_eqn, local_unknown) +=
                  (1.0 / Pe_s) * dpsifdS(l2, i) * dpsifdS(l, i) * W * J;
              }
            }


            // Loop over the velocity components
            for (unsigned i2 = 0; i2 < n_dim; i2++)
            {
              // Get the unknown
              local_unknown = this->nodal_local_eqn(l2, u_index[i2]);


              // If not a boundary condition
              if (local_unknown >= 0)
              {
                // Bits from the advection term
                jacobian(local_eqn, local_unknown) +=
                  (interpolated_C * dpsifdS_div(l2, i2) +
                   psif(l2) * interpolated_grad_C[i2]) *
                  psif(l) * W * J;
              }
            }
          }
        }
      }

      if (flag)
      {
        const double dsigma = this->dsigma_dC(s);
        const double Ca = this->ca();
        for (unsigned l2 = 0; l2 < n_node; l2++)
        {
          local_unknown = this->nodal_local_eqn(l2, this->C_index[l2]);
          if (local_unknown >= 0)
          {
            const double psi_ = psif(l2);
            for (unsigned i = 0; i < n_dim; i++)
            {
              // Add the Jacobian contribution from the surface tension
              local_eqn = this->nodal_local_eqn(l, u_index[i]);
              if (local_eqn >= 0)
              {
                jacobian(local_eqn, local_unknown) -=
                  (dsigma / Ca) * psi_ * dpsifdS_div(l, i) * J * W;
              }
            }
          }
        }
      }

    } // End of loop over the nodes
  }

  //=======================================================
  /// Overload the output function
  //=======================================================
  void SurfactantTransportInterfaceElement::output(std::ostream& outfile,
                                                   const unsigned& n_plot)
  {
    outfile.precision(16);

    const unsigned el_dim = this->dim();
    const unsigned n_dim = this->nodal_dimension();
    const unsigned n_velocity = this->U_index_interface.size();

    // Set output Vector
    Vector<double> s(el_dim);
    Vector<double> n(n_dim);
    Vector<double> u(n_velocity);


    outfile << this->tecplot_zone_string(n_plot);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(n_plot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, n_plot, s);
      // Get the outer unit normal
      this->outer_unit_normal(s, n);

      double u_n = 0.0;
      for (unsigned i = 0; i < n_velocity; i++)
      {
        u[i] = this->interpolated_u(s, i);
      }

      // Not the same as above for axisymmetric case
      for (unsigned i = 0; i < n_dim; i++)
      {
        u_n += u[i] * n[i];
      }

      // Output the x,y,u,v
      for (unsigned i = 0; i < n_dim; i++)
        outfile << this->interpolated_x(s, i) << " ";
      for (unsigned i = 0; i < n_dim; i++) outfile << u[i] << " ";
      // Output a dummy pressure
      outfile << 0.0 << " ";
      // Output the concentration
      outfile << interpolated_C(s) << " ";
      // Output the interfacial tension
      outfile << sigma(s) << " ";
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << u[i] - u_n * n[i] << " ";
      }
      outfile << std::endl;
    }
    outfile << std::endl;
  }

  //=======================================================================
  /// \short Compute the concentration intergated over the surface area
  //=======================================================================
  double SurfactantTransportInterfaceElement::integrate_c()
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    const unsigned el_dim = this->dim();
    const unsigned n_dim = this->nodal_dimension();

    // Set up memeory for the shape functions
    Shape psif(n_node);
    DShape dpsifds(n_node, el_dim);
    DShape dpsifdS(n_node, n_dim);
    DShape dpsifdS_div(n_node, n_dim);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    // Storage for the local coordinate
    Vector<double> s(el_dim);

    // Storage for the total mass
    double mass = 0.0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the local coordinate at the integration point
      for (unsigned i = 0; i < el_dim; i++)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double W = this->integral_pt()->weight(ipt);

      // Call the derivatives of the shape function
      this->dshape_local_at_knot(ipt, psif, dpsifds);

      Vector<double> interpolated_x(n_dim, 0.0);
      DenseMatrix<double> interpolated_t(el_dim, n_dim, 0.0);
      double interpolated_c = 0.0;

      // Loop over the shape functions
      for (unsigned l = 0; l < n_node; l++)
      {
        const double psi_ = psif(l);
        interpolated_c += this->nodal_value(l, this->C_index[l]) * psi_;
        // Loop over directional components
        for (unsigned i = 0; i < n_dim; i++)
        {
          // Coordinate
          interpolated_x[i] += this->nodal_position(l, i) * psi_;

          // Calculate the tangent vectors
          for (unsigned j = 0; j < el_dim; j++)
          {
            interpolated_t(j, i) += this->nodal_position(l, i) * dpsifds(l, j);
          }
        }
      }

      // Calculate the surface gradient and divergence
      double J = this->compute_surface_derivatives(
        psif, dpsifds, interpolated_t, interpolated_x, dpsifdS, dpsifdS_div);

      mass += interpolated_c * W * J;
    }
    return mass;
  }


} // namespace oomph
