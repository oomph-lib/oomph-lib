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
// Non-inline member function of the flux transport elements class

#include "flux_transport_elements.h"
#include "../generic/shape.h"
#include "../generic/timesteppers.h"

namespace oomph
{
  //======================================================================
  /// Return the derivatives of the flux components as functions of the
  /// unknowns. The default implementation is to use second-order
  /// finite-differences, but these can be overloaded in specific equations
  //========================================================================
  template<unsigned DIM>
  void FluxTransportEquations<DIM>::dflux_du(const Vector<double>& u,
                                             RankThreeTensor<double>& df_du)
  {
    // Find the number of fluxes
    const unsigned n_flux = nflux();
    // Local copy of the unknowns
    Vector<double> u_local = u;
    // Finite differences
    DenseMatrix<double> F(n_flux, DIM), F_plus(n_flux, DIM),
      F_minus(n_flux, DIM);
    const double fd_step = GeneralisedElement::Default_fd_jacobian_step;
    // Now loop over all the fluxes
    for (unsigned p = 0; p < n_flux; p++)
    {
      // Store the old value
      const double old_var = u_local[p];
      // Increment the value
      u_local[p] += fd_step;
      // Get the new values
      flux(u_local, F_plus);

      // Reset the value
      u_local[p] = old_var;
      // Decrement the value
      u_local[p] -= fd_step;
      // Get the new values
      flux(u_local, F_minus);

      // Assemble the entries of the jacobian
      for (unsigned r = 0; r < n_flux; r++)
      {
        for (unsigned j = 0; j < DIM; j++)
        {
          df_du(r, j, p) = (F_plus(r, j) - F_minus(r, j)) / (2.0 * fd_step);
        }
      }

      // Reset the value
      u_local[p] = old_var;
    }
  }

  //=====================================================================
  /// Compute the residuals for the flux transport equations;
  /// flag=1 compute jacobian as well
  /// flag=2 compute mass matrix and jacobian as well
  /// flag=3 compute mass matrix as well
  //=======================================================================
  template<unsigned DIM>
  void FluxTransportEquations<DIM>::
    fill_in_generic_residual_contribution_flux_transport(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Find the number of fluxes
    const unsigned n_flux = this->nflux();
    // Find the number of nodes
    const unsigned n_node = this->nnode();
    // Storage for the shape function and derivatives of shape function
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

    // Cache the nodal indices at which the unknowns are stored
    unsigned u_nodal_index[n_flux];
    for (unsigned i = 0; i < n_flux; i++)
    {
      u_nodal_index[i] = this->u_index_flux_transport(i);
    }

    // Find the number of integration points
    unsigned n_intpt = this->integral_pt()->nweight();

    // Integers to store the local equations and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the shape functions at the knot
      double J = this->dshape_and_dtest_eulerian_at_knot_flux_transport(
        ipt, psi, dpsidx, test, dtestdx);

      // Get the integral weight
      double W = this->integral_pt()->weight(ipt) * J;

      // Storage for the local unknowns
      Vector<double> interpolated_u(n_flux, 0.0);
      Vector<double> interpolated_dudt(n_flux, 0.0);

      // Loop over the shape functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Locally cache the shape function
        const double psi_ = psi(l);
        for (unsigned i = 0; i < n_flux; i++)
        {
          // Calculate the velocity and tangent vector
          interpolated_u[i] += this->nodal_value(l, u_nodal_index[i]) * psi_;
          interpolated_dudt[i] += this->du_dt_flux_transport(l, i) * psi_;
        }
      }

      // Calculate the flux at the integration point
      DenseMatrix<double> F(n_flux, DIM);
      this->flux(interpolated_u, F);

      RankThreeTensor<double> dF_du(n_flux, DIM, n_flux);
      if ((flag) && (flag != 3))
      {
        this->dflux_du(interpolated_u, dF_du);
      }

      // We need to assemble the contributions to the volume integral
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the unknowns
        for (unsigned i = 0; i < n_flux; i++)
        {
          // Get the local equation number
          local_eqn = this->nodal_local_eqn(l, i);

          // if it's not a boundary condition
          if (local_eqn >= 0)
          {
            // Add the time derivatives
            residuals[local_eqn] -= interpolated_dudt[i] * test(l) * W;


            // Calculate the flux dot product
            double flux_dot = 0.0;
            for (unsigned k = 0; k < DIM; k++)
            {
              flux_dot += F(i, k) * dtestdx(l, k);
            }

            // Add the contribution to the residuals
            residuals[local_eqn] += flux_dot * W;

            // Now worry about the jacobian and mass matrix terms
            if (flag)
            {
              // If we are assembling the jacobian
              if (flag < 3)
              {
                // Loop over the shape functions again
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Loop over the unknowns again
                  for (unsigned i2 = 0; i2 < n_flux; i2++)
                  {
                    // Get the local unknowns
                    local_unknown = this->nodal_local_eqn(l2, i2);
                    // If not pinned
                    if (local_unknown >= 0)
                    {
                      // Add the time derivative terms
                      if (i2 == i)
                      {
                        jacobian(local_eqn, local_unknown) -=
                          node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                          psi(l2) * test(l) * W;

                        // Add the mass matrix if we are computing
                        // it and the jacobian
                        if (flag == 2)
                        {
                          mass_matrix(local_eqn, local_unknown) +=
                            psi(l2) * test(l) * W;
                        }
                      }

                      // Add the flux derivative terms
                      double dflux_dot = 0.0;
                      for (unsigned k = 0; k < DIM; k++)
                      {
                        dflux_dot += dF_du(i, k, i2) * dtestdx(l, k);
                      }
                      jacobian(local_eqn, local_unknown) +=
                        dflux_dot * psi(l2) * W;
                    }
                  }
                }
              }
              // End of jacobian assembly, here we are just assembling the
              // mass matrix
              else
              {
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  local_unknown = this->nodal_local_eqn(l2, i);
                  if (local_unknown >= 0)
                  {
                    mass_matrix(local_eqn, local_unknown) +=
                      psi(l2) * test(l) * W;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  //==================================================================
  /// Return the i-th unknown at the local coordinate s
  //==================================================================
  template<unsigned DIM>
  double FluxTransportEquations<DIM>::interpolated_u_flux_transport(
    const Vector<double>& s, const unsigned& i)
  {
    // Find the number of nodes
    const unsigned n_node = this->nnode();
    // Get the shape functions at the local coordinate
    Shape psi(n_node);
    this->shape(s, psi);

    // Now interpolate each unknown
    double u = 0.0;

    // Loop over the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      const double psi_ = psi[n];
      u += this->nodal_value(n, u_index_flux_transport(i)) * psi_;
    }
    return u;
  }

  //======================================================================
  /// i-th component of du/dt at local node n.
  /// Uses suitably interpolated value for hanging nodes.
  //======================================================================
  template<unsigned DIM>
  double FluxTransportEquations<DIM>::du_dt_flux_transport(
    const unsigned& n, const unsigned& i) const
  {
    // Get the data's timestepper
    TimeStepper* time_stepper_pt = this->node_pt(n)->time_stepper_pt();

    // Initialise dudt
    double dudt = 0.0;

    // Loop over the timesteps, if there is a non Steady timestepper
    if (!time_stepper_pt->is_steady())
    {
      // Find the index at which the dof is stored
      const unsigned u_nodal_index = this->u_index_flux_transport(i);

      // Number of timsteps (past & present)
      const unsigned n_time = time_stepper_pt->ntstorage();
      // Loop over the timesteps
      for (unsigned t = 0; t < n_time; t++)
      {
        dudt +=
          time_stepper_pt->weight(1, t) * nodal_value(t, n, u_nodal_index);
      }
    }

    return dudt;
  }

  //==================================================================
  /// Calculate the averages of the flux variables
  //=================================================================
  template<unsigned DIM>
  void FluxTransportEquations<DIM>::calculate_element_averages(
    double*& average_value)
  {
    // Find the number of fluxes
    const unsigned n_flux = this->nflux();
    // Resize the memory if necessary
    if (average_value == 0)
    {
      average_value = new double[n_flux + 1];
    }

    // Initialise the averages to zero
    for (unsigned i = 0; i < n_flux + 1; i++)
    {
      average_value[i] = 0.0;
    }

    // Find the number of nodes
    const unsigned n_node = this->nnode();
    // Storage for the shape functions
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

    // Cache the nodal indices at which the unknowns are stored
    unsigned u_nodal_index[n_flux];
    for (unsigned i = 0; i < n_flux; i++)
    {
      u_nodal_index[i] = this->u_index_flux_transport(i);
    }

    // Find the number of integration points
    unsigned n_intpt = this->integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the shape functions and derivatives at the knot
      double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

      // Get the integral weight
      double W = this->integral_pt()->weight(ipt) * J;

      // Storage for the local unknowns
      Vector<double> interpolated_u(n_flux, 0.0);

      // Loop over the shape functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Locally cache the shape function
        const double psi_ = psi(l);
        for (unsigned i = 0; i < n_flux; i++)
        {
          // Calculate the velocity and tangent vector
          interpolated_u[i] += this->nodal_value(l, u_nodal_index[i]) * psi_;
        }
      }

      average_value[n_flux] += W;
      // Loop over the values
      for (unsigned i = 0; i < n_flux; i++)
      {
        average_value[i] += interpolated_u[i] * W;
      }
    }

    // Divide the results by the size of the element
    for (unsigned i = 0; i < n_flux; i++)
    {
      average_value[i] /= average_value[n_flux];
    }
  }


  //====================================================================
  /// Output function, print the values of all unknowns
  //==================================================================
  template<unsigned DIM>
  void FluxTransportEquations<DIM>::output(std::ostream& outfile,
                                           const unsigned& nplot)
  {
    // Find the number of fluxes
    const unsigned n_flux = this->nflux();

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      // Values
      for (unsigned i = 0; i < n_flux; i++)
      {
        outfile << interpolated_u_flux_transport(s, i) << " ";
      }

      outfile << std::endl;
    }
    outfile << std::endl;

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }

  template class FluxTransportEquations<1>;
  template class FluxTransportEquations<2>;
  template class FluxTransportEquations<3>;

} // namespace oomph
