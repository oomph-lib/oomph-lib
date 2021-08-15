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
#include "periodic_orbit_handler.h"

namespace oomph
{
  void PeriodicOrbitEquations::orbit_output(GeneralisedElement* const& elem_pt,
                                            std::ostream& outfile,
                                            const unsigned& n_plot)
  {
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidt(n_node, 1);

    Vector<double> s(1);

    // Cast the element
    PeriodicOrbitBaseElement* const base_el_pt =
      dynamic_cast<PeriodicOrbitBaseElement*>(elem_pt);

    const double inverse_timescale = this->omega();

    // Loop over the plot point
    for (unsigned i = 0; i < n_plot; i++)
    {
      // Get the coordinate
      s[0] = -1.0 + 2.0 * i / ((double)(n_plot - 1));

      (void)dshape_eulerian(s, psi, dpsidt);

      // Sort out the timestepper weights and the time in here
      // Calculate the time
      double interpolated_time = 0.0;
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_time += raw_nodal_position(l, 0) * psi(l);
      }


      // Need to multiply by period and the set global time
      this->time_pt()->time() = interpolated_time;

      // Set the weights of the timestepper
      this->set_timestepper_weights(psi, dpsidt);

      base_el_pt->spacetime_output(
        outfile, n_plot, interpolated_time / inverse_timescale);
    }
  }

  void PeriodicOrbitEquations::fill_in_generic_residual_contribution_orbit(
    PeriodicOrbitAssemblyHandlerBase* const& assembly_handler_pt,
    GeneralisedElement* const& elem_pt,
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    const unsigned& flag)
  {
    // A simple integration loop
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidt(n_node, 1), dtestdt(n_node, 1);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Integers to store the local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0;

    // Storage for the underlying element's residuals
    const unsigned n_elem_dof = elem_pt->ndof();
    Vector<double> el_residuals(n_elem_dof);
    DenseMatrix<double> el_mass;
    DenseMatrix<double> el_jacobian;

    if (flag)
    {
      el_mass.resize(n_elem_dof, n_elem_dof);
      el_jacobian.resize(n_elem_dof, n_elem_dof);
    }

    // Storage for the current value of the unkowns
    Vector<double> u(n_elem_dof);
    // and the derivatives
    Vector<double> du_dt(n_elem_dof);
    // And the previous derivative
    Vector<double> du_dt_old(n_elem_dof);
    // Storage for the inner product matrix
    DenseMatrix<double> inner_product(n_elem_dof);
    // Cast the element
    PeriodicOrbitBaseElement* const base_el_pt =
      dynamic_cast<PeriodicOrbitBaseElement*>(elem_pt);
    // Get the inner product matrix
    base_el_pt->get_inner_product_matrix(inner_product);

    // Storage for "all" unknowns
    Vector<double> all_current_unknowns;
    Vector<double> all_previous_unknowns;

    // Get the value of omega
    const double inverse_timescale = this->omega();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_orbit(
        ipt, psi, dpsidt, test, dtestdt);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Sort out the timestepper weights and the time in here
      // Calculate the time
      double interpolated_time = 0.0;
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_time += raw_nodal_position(l, 0) * psi(l);
      }


      // Need to multiply by period and the set global time
      this->time_pt()->time() = interpolated_time / inverse_timescale;
      // Set the weights of the timestepper
      this->set_timestepper_weights(psi, dpsidt);

      // Get the residuals of the element or residuals and jacobian
      if (flag)
      {
        elem_pt->get_jacobian_and_mass_matrix(
          el_residuals, el_jacobian, el_mass);
      }
      else
      {
        elem_pt->get_residuals(el_residuals);
      }

      // Get the current value of the unknown time derivatives
      base_el_pt->get_non_external_ddofs_dt(du_dt);

      // Multiply the elemental residuals by the appropriate test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the local equation number
        local_eqn = nodal_local_eqn(l, 0);
        // If it's not a boundary condition (it will never be)
        if (local_eqn >= 0)
        {
          // Work out the offset which is the GLOBAL time unknown
          // multiplied by the number of elemental unknowns
          unsigned offset = this->eqn_number(local_eqn) * n_elem_dof;
          // Add to the appropriate residuals
          for (unsigned i = 0; i < n_elem_dof; i++)
          {
            residuals[i + offset] += el_residuals[i] * psi(l) * W;
          }


          // Add the jacobian contributions
          if (flag)
          {
            // The form of the equations is -M du/dt + R = 0

            // Now add the contribution to the jacobian
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              local_unknown = nodal_local_eqn(l2, 0);
              if (local_unknown >= 0)
              {
                // Work out the second offset
                unsigned offset2 = this->eqn_number(local_unknown) * n_elem_dof;
                // Add to the appropriate jacobian terms
                for (unsigned i = 0; i < n_elem_dof; i++)
                {
                  for (unsigned j = 0; j < n_elem_dof; j++)
                  {
                    // Add in the Jacobian terms
                    jacobian(i + offset, j + offset2) +=
                      el_jacobian(i, j) * psi(l2) * psi(l) * W;

                    // Add the time derivative terms,
                    jacobian(i + offset, j + offset2) -=
                      el_mass(i, j) * dpsidt(l2, 0) * inverse_timescale *
                      psi(l) * W;
                  }
                }
              }
            }


            // Add the variation of the period
            for (unsigned i = 0; i < n_elem_dof; i++)
            {
              for (unsigned j = 0; j < n_elem_dof; j++)
              {
                jacobian(i + offset, Ntstorage * n_elem_dof) -=
                  el_mass(i, j) * (du_dt[j] / inverse_timescale) * psi(l) * W;
              }
            }
          } // End of Jacobian flag
        }
      }

      // Sort out the phase condition

      // Get the current value of the unknowns
      base_el_pt->get_non_external_dofs(u);

      // Now get the unknowns required by the assembly handler
      // i.e. including all values throughout the period for backup
      assembly_handler_pt->get_dofs_for_element(elem_pt, all_current_unknowns);
      // Get the previous values as stored by the assembly handler
      assembly_handler_pt->get_previous_dofs_for_element(elem_pt,
                                                         all_previous_unknowns);

      // Now set the elemental values to the previous values
      assembly_handler_pt->set_dofs_for_element(elem_pt, all_previous_unknowns);
      // Get the previous time derivatives
      base_el_pt->get_non_external_ddofs_dt(du_dt_old);
      // Reset the element's values to the current
      assembly_handler_pt->set_dofs_for_element(elem_pt, all_current_unknowns);


      // Assemble the inner product
      double sum = 0.0;
      for (unsigned i = 0; i < n_elem_dof; i++)
      {
        for (unsigned j = 0; j < n_elem_dof; j++)
        {
          sum += u[i] * inner_product(i, j) * du_dt_old[j];
        }
      }

      // Add to the residuals
      residuals[Ntstorage * n_elem_dof] += sum * W;

      // Sort out the jacobian
      if (flag)
      {
        // Loop over the unknown time points
        for (unsigned l2 = 0; l2 < n_node; l2++)
        {
          // Get the local unknown
          local_unknown = nodal_local_eqn(l2, 0);
          if (local_unknown >= 0)
          {
            // Work out the offset
            unsigned offset2 = this->eqn_number(local_unknown) * n_elem_dof;
            // Now add in the appropriate jacobian terms
            for (unsigned i2 = 0; i2 < n_elem_dof; i2++)
            {
              double sum2 = 0.0;
              for (unsigned j = 0; j < n_elem_dof; j++)
              {
                sum2 += inner_product(i2, j) * du_dt_old[j];
              }
              jacobian(Ntstorage * n_elem_dof, i2 + offset2) +=
                psi(l2) * sum2 * W;
            }
          }
        }
      } // End of jacobian calculation

    } // End of loop over time period integration points
  }


} // namespace oomph
