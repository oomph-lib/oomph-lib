// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
// Non-inline functions for space-time Navier-Stokes elements
#include "refineable_discontinuous_galerkin_spacetime_navier_stokes_elements.h"

/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////

namespace oomph
{
  //===================================================================
  /// Compute the diagonal of the velocity/pressure mass matrices.
  /// If which one=0, both are computed, otherwise only the pressure
  /// (which_one=1) or the velocity mass matrix (which_one=2 -- the
  /// LSC version of the preconditioner only needs that one)
  //===================================================================
  template<unsigned DIM>
  void RefineableSpaceTimeNavierStokesEquations<DIM>::
    get_pressure_and_velocity_mass_matrix_diagonal(
      Vector<double>& press_mass_diag,
      Vector<double>& veloc_mass_diag,
      const unsigned& which_one)
  {
    // Resize and initialise
    unsigned n_dof = ndof();

    // If the pressure mass matrix needs to be constructed
    if ((which_one == 0) || (which_one == 1))
    {
      // Assign the appropriate amount of space
      press_mass_diag.assign(n_dof, 0.0);
    }

    // If the velocity mass matrix needs to be constructed
    if ((which_one == 0) || (which_one == 2))
    {
      // Assign the appropriate amount of space
      veloc_mass_diag.assign(n_dof, 0.0);
    }

    // Pointer to hang info object
    HangInfo* hang_info_pt = 0;

    // Number of master nodes
    unsigned n_master = 1;

    // Hang weight for shape functions
    double hang_weight = 1.0;

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Storage for the local velocity indices
    Vector<unsigned> u_nodal_index(DIM, 0.0);

    // Find the indices at which the local velocities are stored
    for (unsigned i = 0; i < DIM; i++)
    {
      // Calculate the i-th local velocity component
      u_nodal_index[i] = this->u_index_nst(i);
    }

    // Set up memory for velocity shape functions
    Shape psi(n_node);

    // Find number of pressure dofs
    unsigned n_pres = this->npres_nst();

    // Pressure shape function
    Shape psi_p(n_pres);

    // Which nodal value represents the pressure? (Negative if pressure
    // is not based on nodal interpolation).
    int p_index = this->p_nodal_index_nst();

    // Local array of booleans that are true if the l-th pressure value is
    // hanging (avoid repeated virtual function calls)
    bool pressure_dof_is_hanging[n_pres];

    // If the pressure is stored at a node
    if (p_index >= 0)
    {
      // Loop over the pressure nodes
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Read out the hang status of the pressure node
        pressure_dof_is_hanging[l] = pressure_node_pt(l)->is_hanging(p_index);
      }
    }
    // Otherwise the pressure is not stored at a node and so cannot hang
    else
    {
      // Loop over the pressures nodes
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Indicate that the pressure node is not hanging
        pressure_dof_is_hanging[l] = false;
      }
    } // if (p_index>=0)

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Integer to store the local equations no
    int local_eqn = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get determinant of Jacobian of the mapping
      double J = J_eulerian_at_knot(ipt);

      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        // Calculate the i-th local coordinate at the ipt-th integration point
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Premultiply weights and Jacobian
      double W = w * J;

      // Do we want the velocity one?
      if ((which_one == 0) || (which_one == 2))
      {
        // Get the velocity shape functions
        shape_at_knot(ipt, psi);

        // Loop over the velocity shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Local boolean to indicate whether or not the node is hanging
          bool is_node_hanging = node_pt(l)->is_hanging();

          // If the node is hanging
          if (is_node_hanging)
          {
            // Get the HangInfo pointer from the node
            hang_info_pt = node_pt(l)->hanging_pt();

            // Read out number of master nodes from hanging data
            n_master = hang_info_pt->nmaster();
          }
          // Otherwise the node is its own master
          else
          {
            // There is only one node to consider; the node itself
            n_master = 1;
          }

          // Loop over the master nodes
          for (unsigned m = 0; m < n_master; m++)
          {
            // Loop over velocity components for equations
            for (unsigned i = 0; i < DIM; i++)
            {
              // If the node is hanging
              if (is_node_hanging)
              {
                // Get the equation number from the master node
                local_eqn = this->local_hang_eqn(
                  hang_info_pt->master_node_pt(m), u_nodal_index[i]);

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

              // If it's not a boundary condition
              if (local_eqn >= 0)
              {
                // Add the contribution
                veloc_mass_diag[local_eqn] += pow(psi[l] * hang_weight, 2) * W;
              }
            } // for (unsigned i=0;i<DIM;i++)
          } // for (unsigned m=0;m<n_master;m++)
        } // if ((which_one==0)||(which_one==2))
      } // for (unsigned ipt=0;ipt<n_intpt;ipt++)

      // Do we want the pressure one?
      if ((which_one == 0) || (which_one == 1))
      {
        // Get the pressure shape functions
        this->pshape_nst(s, psi_p);

        // Loop over the pressure shape functions
        for (unsigned l = 0; l < n_pres; l++)
        {
          // If the pressure dof is hanging
          if (pressure_dof_is_hanging[l])
          {
            // Get the HangInfo object
            hang_info_pt = this->pressure_node_pt(l)->hanging_pt(p_index);

            // Get the number of master nodes from the pressure node
            n_master = hang_info_pt->nmaster();
          }
          // Otherwise the node is its own master
          else
          {
            // There is only one node to consider; the node itself
            n_master = 1;
          }

          // Loop over the master nodes
          for (unsigned m = 0; m < n_master; m++)
          {
            // If the pressure dof is hanging
            if (pressure_dof_is_hanging[l])
            {
              // Get the local equation from the master node
              local_eqn =
                this->local_hang_eqn(hang_info_pt->master_node_pt(m), p_index);

              // Get the weight for the node
              hang_weight = hang_info_pt->master_weight(m);
            }
            // If the pressure dof is not hanging
            else
            {
              // Get the local equation number of this pressure dof
              local_eqn = this->p_local_eqn(l);

              // Get the hang weight associated with this pressure node
              hang_weight = 1.0;
            }

            // If the equation is not pinned
            if (local_eqn >= 0)
            {
              // Add the contribution
              press_mass_diag[local_eqn] += pow(psi_p[l] * hang_weight, 2) * W;
            }
          } // for (unsigned m=0;m<n_master;m++)
        } // for (unsigned l=0;l<n_pres;l++)
      } // if ((which_one==0)||(which_one==1))
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of get_pressure_and_velocity_mass_matrix_diagonal


  //==============================================================
  /// Compute the residuals for the associated pressure advection
  /// diffusion problem. Used by the Fp preconditioner.
  /// flag=1(or 0): do (or don't) compute the Jacobian as well.
  //==============================================================
  template<unsigned DIM>
  void RefineableSpaceTimeNavierStokesEquations<DIM>::
    fill_in_generic_pressure_advection_diffusion_contribution_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
  {
    OomphLibWarning("I'm not sure this is correct yet...",
                    OOMPH_CURRENT_FUNCTION,
                    OOMPH_EXCEPTION_LOCATION);

    // Return immediately if there are no dofs
    if (ndof() == 0) return;

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many pressure dofs there are
    unsigned n_pres = this->npres_nst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_nodal_index[i] = this->u_index_nst(i);
    }


    // Which nodal value represents the pressure? (Negative if pressure
    // is not based on nodal interpolation).
    int p_index = this->p_nodal_index_nst();

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
      // pressure advection diffusion doesn't work for this one!
      throw OomphLibError(
        "Pressure advection diffusion does not work in this case\n",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);

      for (unsigned l = 0; l < n_pres; ++l)
      {
        pressure_dof_is_hanging[l] = false;
      }
    }

    // Set up memory for the velocity shape fcts
    Shape psif(n_node);
    DShape dpsidx(n_node, DIM);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);
    DShape dpsip(n_pres, DIM);
    DShape dtestp(n_pres, DIM);

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get Physical Variables from Element
    // Reynolds number must be multiplied by the density ratio
    double scaled_re = this->re() * this->density_ratio();

    // Integers to store the local equations and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Pointers to hang info objects
    HangInfo *hang_info_pt = 0, *hang_info2_pt = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++) s[i] = integral_pt()->knot(ipt, i);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the veloc shape functions
      // (Derivs not needed but they are free)
      double J = this->dshape_eulerian_at_knot(ipt, psif, dpsidx);

      // Call the pressure shape and test functions
      this->dpshape_and_dptest_eulerian_nst(s, psip, dpsip, testp, dtestp);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the pressure and velocity components
      // Allocate
      Vector<double> interpolated_u(DIM, 0.0);
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> interpolated_dpdx(DIM, 0.0);

      // Calculate pressure gradient
      for (unsigned l = 0; l < n_pres; l++)
      {
        for (unsigned i = 0; i < DIM; i++)
        {
          interpolated_dpdx[i] += this->p_nst(l) * dpsip(l, i);
        }
      }

      // Calculate velocities

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the nodal value
          double u_value = nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif[l];
          interpolated_x[i] += nodal_position(l, i) * psif[l];
        }
      }

      // Source function (for validaton only)
      double source = 0.0;
      if (this->Press_adv_diff_source_fct_pt != 0)
      {
        source = this->Press_adv_diff_source_fct_pt(interpolated_x);
      }


      // Number of master nodes and storage for the weight of the shape function
      unsigned n_master = 1;
      double hang_weight = 1.0;


      // Loop over the pressure shape functions
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

          // If the equation is not pinned
          if (local_eqn >= 0)
          {
            residuals[local_eqn] -= source * testp[l] * W * hang_weight;
            for (unsigned k = 0; k < DIM; k++)
            {
              residuals[local_eqn] +=
                interpolated_dpdx[k] *
                (scaled_re * interpolated_u[k] * testp[l] + dtestp(l, k)) * W *
                hang_weight;
            }

            // Jacobian too?
            if (flag)
            {
              // Number of master nodes and weights
              unsigned n_master2 = 1;
              double hang_weight2 = 1.0;

              // Loop over the pressure shape functions
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

                  // If the unknown is not pinned
                  if (local_unknown >= 0)
                  {
                    if ((int(eqn_number(local_eqn)) !=
                         this->Pinned_fp_pressure_eqn) &&
                        (int(eqn_number(local_unknown)) !=
                         this->Pinned_fp_pressure_eqn))
                    {
                      for (unsigned k = 0; k < DIM; k++)
                      {
                        jacobian(local_eqn, local_unknown) +=
                          dtestp(l2, k) *
                          (scaled_re * interpolated_u[k] * testp[l] +
                           dtestp(l, k)) *
                          W * hang_weight * hang_weight2;
                      }
                    }
                    else
                    {
                      if ((int(eqn_number(local_eqn)) ==
                           this->Pinned_fp_pressure_eqn) &&
                          (int(eqn_number(local_unknown)) ==
                           this->Pinned_fp_pressure_eqn))
                      {
                        jacobian(local_eqn, local_unknown) = 1.0;
                      }
                    }
                  }
                }
              }
            } /*End of Jacobian calculation*/
          } // End of if not boundary condition
        } // End of loop over master nodes
      } // End of loop over l
    } // end of integration loop

    // Now add boundary contributions from Robin BCs
    unsigned nrobin =
      this->Pressure_advection_diffusion_robin_element_pt.size();
    for (unsigned e = 0; e < nrobin; e++)
    {
      this->Pressure_advection_diffusion_robin_element_pt[e]
        ->fill_in_generic_residual_contribution_fp_press_adv_diff_robin_bc(
          residuals, jacobian, flag);
    }
  }


  //========================================================================
  /// Add element's contribution to the elemental
  /// residual vector and/or Jacobian matrix.
  /// flag=1: compute both
  /// flag=0: compute only residual vector
  //========================================================================
  template<unsigned DIM>
  void RefineableSpaceTimeNavierStokesEquations<DIM>::
    fill_in_generic_residual_contribution_nst(Vector<double>& residuals,
                                              DenseMatrix<double>& jacobian,
                                              DenseMatrix<double>& mass_matrix,
                                              const unsigned& flag)
  {
    // Return immediately if there are no dofs
    if (ndof() == 0) return;

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many pressure dofs there are
    unsigned n_pres = this->npres_nst();

    // Allocate storage for the indices of the velocity components
    unsigned u_nodal_index[DIM];

    // Loop over the velocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Find the index at which the i-th local velocity is stored
      u_nodal_index[i] = this->u_index_nst(i);
    }

    // Which nodal value represents the pressure? (Negative if pressure
    // is not based on nodal interpolation).
    int p_index = this->p_nodal_index_nst();

    // Local array of booleans that are true if the l-th pressure value is
    // hanging (avoid repeated virtual function calls)
    bool pressure_dof_is_hanging[n_pres];

    // If the pressure is stored at a node
    if (p_index >= 0)
    {
      // Loop over the pressure nodes
      for (unsigned l = 0; l < n_pres; ++l)
      {
        // Get the hang status of the l-th pressure node
        pressure_dof_is_hanging[l] = pressure_node_pt(l)->is_hanging(p_index);
      }
    }
    // Otherwise the pressure is not stored at a node and so cannot hang
    else
    {
      // Loop over the pressure nodes
      for (unsigned l = 0; l < n_pres; ++l)
      {
        // Indicate that the l-th pressure node is not hangings
        pressure_dof_is_hanging[l] = false;
      }
    } // if (p_index>=0)

    // Set up memory for the shape functions
    Shape psif(n_node);

    // Set up memory for the test functions
    Shape testf(n_node);

    // Allocate space for the derivatives of the shape functions
    DShape dpsifdx(n_node, DIM + 1);

    // Allocate space for the derivatives of the test functions
    DShape dtestfdx(n_node, DIM + 1);

    // Set up memory for pressure shape functions
    Shape psip(n_pres);

    // Set up memory for pressure test functions
    Shape testp(n_pres);

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM + 1, 0.0);

    //-------------------------------------
    // Get physical variables from element:
    //-------------------------------------
    // Reynolds number must be multiplied by the density ratio
    double scaled_re = this->re() * this->density_ratio();

    // Get the scaled Womersley value
    double scaled_re_st = this->re() * this->st() * this->density_ratio();

    // Get the scaled Womersley value differentiated w.r.t. the Strouhal number
    double scaled_dre_st_dst = this->re() * this->density_ratio();

    // Get the scaled Reynolds / Froude number
    double scaled_re_inv_fr = this->re_invfr() * this->density_ratio();

    // Get the viscosity ratio
    double visc_ratio = this->viscosity_ratio();

    // Get the gravity vector
    Vector<double> G = this->g();

    // Integer to store the local equation number
    int local_eqn = 0;

    // Integer to store the local unknown number
    int local_unknown = 0;

    // Pointers to first hang info object
    HangInfo* hang_info_pt = 0;

    // Pointers to second hang info object
    HangInfo* hang_info2_pt = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        // Calculate the i-th local coordinate
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = this->dshape_and_dtest_eulerian_at_knot_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Call the pressure shape and test functions
      this->pshape_nst(s, psip, testp);

      // Pre-multiply the weights and the Jacobian
      double W = w * J;

      // Storage for the interpolated time value
      double interpolated_t = 0.0;

      // Storage for the interpolated pressure value
      double interpolated_p = 0.0;

      // Storage for the spatial coordinates
      Vector<double> interpolated_x(DIM, 0.0);

      // Storage for the interpolated velocity components
      Vector<double> interpolated_u(DIM, 0.0);

      // Storage for the interpolated time-derivative of the velocities
      Vector<double> interpolated_dudt(DIM, 0.0);

      // Storage for the mesh velocity
      Vector<double> mesh_velocity(DIM, 0.0);

      // Storage for the spatial derivatives of the velocity components
      DenseMatrix<double> interpolated_dudx(DIM, DIM, 0.0);

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Update the current approximation to the interpolated pressure
        interpolated_p += this->p_nst(l) * psip[l];
      }

      //-------------------------------------------------------------------
      // Calculate velocities, derivatives, source function and body force:
      //-------------------------------------------------------------------
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Update the interpolated time value
        interpolated_t += this->nodal_position(l, DIM) * psif(l);

        // Loop over coordinate directions
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the nodal value
          double u_value = this->nodal_value(l, u_nodal_index[i]);

          // Update the i-th interpolated velocity component
          interpolated_u[i] += u_value * psif[l];

          // Update the i-th interpolated coordinate value
          interpolated_x[i] += this->nodal_position(l, i) * psif[l];

          // Update the interpolated du_i/dt value
          interpolated_dudt[i] += u_value * dpsifdx(l, DIM);

          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            // Update the interpolated du_i/dx_j value
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        } // for (unsigned i=0;i<DIM;i++)
      } // for (unsigned l=0;l<n_node;l++)

      // If ALE is enabled
      if (!(this->ALE_is_disabled))
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < DIM; i++)
          {
            // Update the i-th mesh velocity component
            mesh_velocity[i] += this->nodal_position(l, i) * dpsifdx(l, DIM);
          }
        } // for (unsigned l=0;l<n_node;l++)
      } // if (!(this->ALE_is_disabled))

      // Allocate space for the body force
      Vector<double> body_force(DIM, 0.0);

      // Get the user-defined body force term
      this->get_body_force_nst(
        interpolated_t, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      double source = this->get_source_nst(interpolated_t, ipt, interpolated_x);

      //---------------------------------
      // Assemble residuals and Jacobian:
      //---------------------------------
      //--------------------
      // Momentum equations:
      //--------------------
      // Number of master nodes
      unsigned n_master = 1;

      // Storage for the weight of the shape function
      double hang_weight = 1.0;

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Local boolean to indicate whether the node is hanging
        bool is_node_hanging = node_pt(l)->is_hanging();

        // If the node is hanging
        if (is_node_hanging)
        {
          // Get the HangInfo object associated with the l-th node
          hang_info_pt = node_pt(l)->hanging_pt();

          // Read out number of master nodes from hanging data
          n_master = hang_info_pt->nmaster();
        }
        // Otherwise the node is its own master
        else
        {
          // There is only one node to consider; the node itself
          n_master = 1;
        }

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Loop over velocity components
          for (unsigned i = 0; i < DIM; i++)
          {
            // Check if the node is hanging
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

            // If it's not a boundary condition
            if (local_eqn >= 0)
            {
              // Add the user-defined body force terms
              residuals[local_eqn] +=
                body_force[i] * testf[l] * W * hang_weight;

              // Add the gravitational body force term
              residuals[local_eqn] +=
                scaled_re_inv_fr * testf[l] * G[i] * W * hang_weight;

              // Add the pressure gradient term
              residuals[local_eqn] +=
                interpolated_p * dtestfdx(l, i) * W * hang_weight;

              // Add in the contribution from the time derivative
              residuals[local_eqn] -= (scaled_re_st * interpolated_dudt[i] *
                                       testf[l] * W * hang_weight);

              // If ALE is enabled
              if (!(this->ALE_is_disabled))
              {
                // Loop over the coordinate directions
                for (unsigned k = 0; k < DIM; k++)
                {
                  // Add in the mesh velocity contribution
                  residuals[local_eqn] +=
                    (scaled_re_st * mesh_velocity[k] * interpolated_dudx(i, k) *
                     testf[l] * W * hang_weight);
                }
              } // if (!(this->ALE_is_disabled))

              // Loop over the coordinate directions
              for (unsigned k = 0; k < DIM; k++)
              {
                // Add in the convective term contribution
                residuals[local_eqn] -=
                  (scaled_re * interpolated_u[k] * interpolated_dudx(i, k) *
                   testf[l] * W * hang_weight);
              }

              // Loop over the coordinate directions
              for (unsigned k = 0; k < DIM; k++)
              {
                // Add in the stress tensor terms:
                // NOTE: The viscosity ratio needs to go in here to ensure
                // continuity of normal stress is satisfied even in flows
                // with zero pressure gradient!
                residuals[local_eqn] -=
                  ((interpolated_dudx(i, k) +
                    this->Gamma[i] * interpolated_dudx(k, i)) *
                   visc_ratio * dtestfdx(l, k) * W * hang_weight);
              }

              //------------------------
              // Calculate the Jacobian:
              //------------------------
              // If we also need to construct the Jacobian
              if (flag)
              {
                // Number of master nodes
                unsigned n_master2 = 1;

                // Storage for the weight of the shape function
                double hang_weight2 = 1.0;

                // Loop over the velocity nodes for columns
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Local boolean to indicate whether the node is hanging
                  bool is_node2_hanging = node_pt(l2)->is_hanging();

                  // If the node is hanging
                  if (is_node2_hanging)
                  {
                    // Get the HangInfo object associated with the l-th node
                    hang_info2_pt = node_pt(l2)->hanging_pt();

                    // Read out number of master nodes from hanging data
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise the node is its own master
                  else
                  {
                    // There is only one node to consider; the node itself
                    n_master2 = 1;
                  }

                  // Loop over the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    // Loop over the velocity components
                    for (unsigned i2 = 0; i2 < DIM; i2++)
                    {
                      // Check if the node is hanging
                      if (is_node2_hanging)
                      {
                        // Get the equation number from the master node
                        local_unknown = this->local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), u_nodal_index[i2]);

                        // Get the hang weight from the master node
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      // If the node is not hanging
                      else
                      {
                        // Get the local equation number
                        local_unknown =
                          this->nodal_local_eqn(l2, u_nodal_index[i2]);

                        // Node contributes with full weight
                        hang_weight2 = 1.0;
                      }

                      // If the unknown is non-zero
                      if (local_unknown >= 0)
                      {
                        // Add contribution to elemental matrix
                        jacobian(local_eqn, local_unknown) -=
                          (visc_ratio * this->Gamma[i] * dpsifdx(l2, i) *
                           dtestfdx(l, i2) * W * hang_weight * hang_weight2);

                        // Now add in the inertial terms
                        jacobian(local_eqn, local_unknown) -=
                          (scaled_re * psif[l2] * interpolated_dudx(i, i2) *
                           testf[l] * W * hang_weight * hang_weight2);

                        // Extra component if i2=i
                        if (i2 == i)
                        {
                          // If we also need to construct the mass matrix (only
                          // diagonal entries)
                          if (flag == 2)
                          {
                            // NOTE: This is positive because the mass matrix is
                            // taken to the other side of the equation when
                            // formulating the generalised eigenproblem.
                            mass_matrix(local_eqn, local_unknown) +=
                              (scaled_re_st * psif[l2] * testf[l] * W *
                               hang_weight * hang_weight2);
                          }

                          // Add in the time-derivative contribution
                          jacobian(local_eqn, local_unknown) -=
                            (scaled_re_st * dpsifdx(l2, DIM) * testf[l] * W *
                             hang_weight * hang_weight2);

                          // Loop over the velocity components
                          for (unsigned k = 0; k < DIM; k++)
                          {
                            // Add in the convective term contribution
                            jacobian(local_eqn, local_unknown) -=
                              (scaled_re * interpolated_u[k] * dpsifdx(l2, k) *
                               testf[l] * W * hang_weight * hang_weight2);
                          }

                          // If ALE is enabled
                          if (!(this->ALE_is_disabled))
                          {
                            // Loop over the velocity components
                            for (unsigned k = 0; k < DIM; k++)
                            {
                              // Add in the mesh velocity contribution
                              jacobian(local_eqn, local_unknown) +=
                                (scaled_re_st * mesh_velocity[k] *
                                 dpsifdx(l2, k) * testf[l] * W * hang_weight *
                                 hang_weight2);
                            }
                          } // if (!(this->ALE_is_disabled))

                          // Loop over velocity components
                          for (unsigned k = 0; k < DIM; k++)
                          {
                            // Add in the velocity gradient terms
                            jacobian(local_eqn, local_unknown) -=
                              (visc_ratio * dpsifdx(l2, k) * dtestfdx(l, k) *
                               W * hang_weight * hang_weight2);
                          }
                        } // if (i2==i)
                      } // if (local_unknown>=0)
                    } // for (unsigned i2=0;i2<DIM;i2++)
                  } // for (unsigned m2=0;m2<n_master2;m2++)
                } // for (unsigned l2=0;l2<n_node;l2++)

                // Loop over the pressure shape functions
                for (unsigned l2 = 0; l2 < n_pres; l2++)
                {
                  // If the pressure dof is hanging
                  if (pressure_dof_is_hanging[l2])
                  {
                    // Get the HangInfo object associated with the l2-th
                    // pressure node
                    hang_info2_pt =
                      this->pressure_node_pt(l2)->hanging_pt(p_index);

                    // Get the number of master nodes from the pressure node
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise the node is its own master
                  else
                  {
                    // There is only one node to consider; the node itself
                    n_master2 = 1;
                  }

                  // Loop over the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    // If the pressure dof is hanging
                    if (pressure_dof_is_hanging[l2])
                    {
                      // Get the unknown from the master node
                      local_unknown = this->local_hang_eqn(
                        hang_info2_pt->master_node_pt(m2), p_index);

                      // Get the weight from the HangInfo object
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // Get the local equation number
                      local_unknown = this->p_local_eqn(l2);

                      // Node contributes with full weight
                      hang_weight2 = 1.0;
                    }

                    // If the unknown is not pinned
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) +=
                        (psip[l2] * dtestfdx(l, i) * W * hang_weight *
                         hang_weight2);
                    }
                  } // for (unsigned m2=0;m2<n_master2;m2++)
                } // for (unsigned l2=0;l2<n_pres;l2++)

                //------------------------------------
                // Calculate external data information
                //------------------------------------
                // If we're storing the Strouhal number as external data then we
                // need to calculate dr/d(St) in the elemental Jacobian
                if (this->Strouhal_is_stored_as_external_data)
                {
                  // The index of the external data (there's only one!)
                  unsigned data_index = 0;

                  // The index of the unknown value stored in the external data
                  unsigned value_index = 0;

                  // Get the local equation number associated with the extra
                  // unknown
                  local_unknown =
                    this->external_local_eqn(data_index, value_index);

                  // If we're at a non-zero degree of freedom add in the entry
                  if (local_unknown >= 0)
                  {
                    // Add in the contribution from the time derivative
                    jacobian(local_eqn, local_unknown) -=
                      (scaled_dre_st_dst * interpolated_dudt[i] * testf[l] * W *
                       hang_weight);

                    // If ALE is enabled
                    if (!this->ALE_is_disabled)
                    {
                      // Loop over the coordinate directions
                      for (unsigned k = 0; k < DIM; k++)
                      {
                        // Add in the mesh velocity contribution
                        jacobian(local_eqn, local_unknown) +=
                          (scaled_dre_st_dst * mesh_velocity[k] *
                           interpolated_dudx(i, k) * testf[l] * W *
                           hang_weight);
                      }
                    } // if (!ALE_is_disabled)
                  } // if (local_unknown>=0)
                } // if (Strouhal_is_stored_as_external_data)
              } // if (flag)
            } // if (local_eqn>=0)
          } // for (unsigned i=0;i<DIM;i++)
        } // for (unsigned m=0;m<n_master;m++)
      } // for (unsigned l=0;l<n_node;l++)

      //---------------------
      // Continuity equation:
      //---------------------
      // Loop over the pressure shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        // If the pressure dof is hanging
        if (pressure_dof_is_hanging[l])
        {
          // Get the HangInfo object associated with the l2-th pressure node
          hang_info_pt = this->pressure_node_pt(l)->hanging_pt(p_index);

          // Get the number of master nodes from the pressure node
          n_master = hang_info_pt->nmaster();
        }
        // Otherwise the node is its own master
        else
        {
          // There is only one node to consider; the node itself
          n_master = 1;
        }

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // If the pressure dof is hanging
          if (pressure_dof_is_hanging[l])
          {
            // Get the unknown from the master node
            local_eqn =
              this->local_hang_eqn(hang_info_pt->master_node_pt(m), p_index);

            // Get the weight from the HangInfo object
            hang_weight = hang_info_pt->master_weight(m);
          }
          // If the node is not hanging
          else
          {
            // Get the local equation number
            local_eqn = this->p_local_eqn(l);

            // Node contributes with full weight
            hang_weight = 1.0;
          }

          // If the equation is not pinned
          if (local_eqn >= 0)
          {
            // Add in the source term contribution
            residuals[local_eqn] -= source * testp[l] * W * hang_weight;

            // Loop over velocity components
            for (unsigned k = 0; k < DIM; k++)
            {
              // Add in the velocity gradient terms
              residuals[local_eqn] +=
                interpolated_dudx(k, k) * testp[l] * W * hang_weight;
            }

            //------------------------
            // Calculate the Jacobian:
            //------------------------
            // If we also need to construct the Jacobian
            if (flag)
            {
              // Number of master nodes
              unsigned n_master2 = 1;

              // Storage for the weight of the shape function
              double hang_weight2 = 1.0;

              // Loop over the velocity nodes for columns
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Local boolean to indicate whether the node is hanging
                bool is_node2_hanging = node_pt(l2)->is_hanging();

                // If the node is hanging
                if (is_node2_hanging)
                {
                  // Get the HangInfo object associated with the l2-th node
                  hang_info2_pt = node_pt(l2)->hanging_pt();

                  // Read out number of master nodes from hanging data
                  n_master2 = hang_info2_pt->nmaster();
                }
                // Otherwise the node is its own master
                else
                {
                  // There is only one node to consider; the node itself
                  n_master2 = 1;
                }

                // Loop over the master nodes
                for (unsigned m2 = 0; m2 < n_master2; m2++)
                {
                  // Loop over the velocity components
                  for (unsigned i2 = 0; i2 < DIM; i2++)
                  {
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Get the equation number from the master node
                      local_unknown = this->local_hang_eqn(
                        hang_info2_pt->master_node_pt(m2), u_nodal_index[i2]);

                      // Get the weight from the HangInfo object
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // Get the local equation number
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_index[i2]);

                      // Node contributes with full weight
                      hang_weight2 = 1.0;
                    }

                    // If the unknown is not pinned
                    if (local_unknown >= 0)
                    {
                      // Add in the velocity gradient contribution
                      jacobian(local_eqn, local_unknown) +=
                        (dpsifdx(l2, i2) * testp[l] * W * hang_weight *
                         hang_weight2);
                    }
                  } // for (unsigned i2=0;i2<DIM;i2++)s
                } // for (unsigned m2=0;m2<n_master2;m2++)
              } // for (unsigned l2=0;l2<n_node;l2++)
            } // if (flag)
          } // if (local_eqn>=0)
        } // for (unsigned m=0;m<n_master;m++)
      } // for (unsigned l=0;l<n_pres;l++)
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of fill_in_generic_residual_contribution_nst


  //======================================================================
  /// Compute derivatives of elemental residual vector with respect
  /// to nodal coordinates.
  /// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
  /// Overloads the FD-based version in the FE base class.
  //======================================================================
  template<unsigned DIM>
  void RefineableSpaceTimeNavierStokesEquations<
    DIM>::get_dresidual_dnodal_coordinates(RankThreeTensor<double>&
                                             dresidual_dnodal_coordinates)
  {
    // Throw a warning
    throw OomphLibError("Space-time update needs to be checked!",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);

    // Return immediately if there are no dofs
    if (ndof() == 0)
    {
      return;
    }

    // Determine number of nodes in element
    const unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Determine number of pressure dofs in element
    const unsigned n_pres = this->npres_nst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_nodal_index[i] = this->u_index_nst(i);
    }

    // Which nodal value represents the pressure? (Negative if pressure
    // is not based on nodal interpolation).
    const int p_index = this->p_nodal_index_nst();

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
    DShape dpsifdx(n_node, DIM), dtestfdx(n_node, DIM);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Determine number of shape controlling nodes
    const unsigned n_shape_controlling_node = nshape_controlling_nodes();

    // Deriatives of shape fct derivatives w.r.t. nodal coords
    RankFourTensor<double> d_dpsifdx_dX(
      DIM, n_shape_controlling_node, n_node, DIM);
    RankFourTensor<double> d_dtestfdx_dX(
      DIM, n_shape_controlling_node, n_node, DIM);

    // Derivative of Jacobian of mapping w.r.t. to nodal coords
    DenseMatrix<double> dJ_dX(DIM, n_shape_controlling_node);

    // Derivatives of derivative of u w.r.t. nodal coords
    RankFourTensor<double> d_dudx_dX(DIM, n_shape_controlling_node, DIM, DIM);

    // Derivatives of nodal velocities w.r.t. nodal coords:
    // Assumption: Interaction only local via no-slip so
    // X_ij only affects U_ij.
    DenseMatrix<double> d_U_dX(DIM, n_shape_controlling_node, 0.0);

    // Determine the number of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get physical variables from element
    // (Reynolds number must be multiplied by the density ratio)
    double scaled_re = this->re() * this->density_ratio();
    double scaled_re_st = this->re_st() * this->density_ratio();
    double scaled_re_inv_fr = this->re_invfr() * this->density_ratio();
    double visc_ratio = this->viscosity_ratio();
    Vector<double> G = this->g();

    // FD step
    double eps_fd = GeneralisedElement::Default_fd_jacobian_step;

    // Pre-compute derivatives of nodal velocities w.r.t. nodal coords:
    // Assumption: Interaction only local via no-slip so
    // X_ij only affects U_ij.
    bool element_has_node_with_aux_node_update_fct = false;

    std::map<Node*, unsigned> local_shape_controlling_node_lookup =
      shape_controlling_node_lookup();

    // FD loop over shape-controlling nodes
    for (std::map<Node*, unsigned>::iterator it =
           local_shape_controlling_node_lookup.begin();
         it != local_shape_controlling_node_lookup.end();
         it++)
    {
      // Get node
      Node* nod_pt = it->first;

      // Get its number
      unsigned q = it->second;

      // Only compute if there's a node-update fct involved
      if (nod_pt->has_auxiliary_node_update_fct_pt())
      {
        element_has_node_with_aux_node_update_fct = true;

        // Current nodal velocity
        Vector<double> u_ref(DIM);
        for (unsigned i = 0; i < DIM; i++)
        {
          u_ref[i] = *(nod_pt->value_pt(u_nodal_index[i]));
        }

        // FD
        for (unsigned p = 0; p < DIM; p++)
        {
          // Make backup
          double backup = nod_pt->x(p);

          // Do FD step. No node update required as we're
          // attacking the coordinate directly...
          nod_pt->x(p) += eps_fd;

          // Do auxiliary node update (to apply no slip)
          nod_pt->perform_auxiliary_node_update_fct();

          // Evaluate
          d_U_dX(p, q) =
            (*(nod_pt->value_pt(u_nodal_index[p])) - u_ref[p]) / eps_fd;

          // Reset
          nod_pt->x(p) = backup;

          // Do auxiliary node update (to apply no slip)
          nod_pt->perform_auxiliary_node_update_fct();
        }
      }
    }

    // Integer to store the local equation number
    int local_eqn = 0;

    // Pointers to hang info object
    HangInfo* hang_info_pt = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      const double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      const double J =
        this->dshape_and_dtest_eulerian_at_knot_nst(ipt,
                                                    psif,
                                                    dpsifdx,
                                                    d_dpsifdx_dX,
                                                    testf,
                                                    dtestfdx,
                                                    d_dtestfdx_dX,
                                                    dJ_dX);

      // Call the pressure shape and test functions
      this->pshape_nst(s, psip, testp);

      // Calculate local values of the pressure and velocity components
      // Allocate
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(DIM, 0.0);
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> mesh_velocity(DIM, 0.0);
      Vector<double> dudt(DIM, 0.0);
      DenseMatrix<double> interpolated_dudx(DIM, DIM, 0.0);

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
      {
        interpolated_p += this->p_nst(l) * psip[l];
      }

      // Calculate velocities and derivatives:

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the nodal value
          const double u_value = nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif[l];
          interpolated_x[i] += nodal_position(l, i) * psif[l];
          dudt[i] += this->du_dt_nst(l, i) * psif[l];

          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      if (!this->ALE_is_disabled)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < DIM; i++)
          {
            mesh_velocity[i] += this->dnodal_position_dt(l, i) * psif[l];
          }
        }
      }

      // Calculate derivative of du_i/dx_k w.r.t. nodal positions X_{pq}

      // Loop over shape-controlling nodes
      for (unsigned q = 0; q < n_shape_controlling_node; q++)
      {
        // Loop over coordinate directions
        for (unsigned p = 0; p < DIM; p++)
        {
          for (unsigned i = 0; i < DIM; i++)
          {
            for (unsigned k = 0; k < DIM; k++)
            {
              double aux = 0.0;
              for (unsigned j = 0; j < n_node; j++)
              {
                aux +=
                  nodal_value(j, u_nodal_index[i]) * d_dpsifdx_dX(p, q, j, k);
              }
              d_dudx_dX(p, q, i, k) = aux;
            }
          }
        }
      }

      // Get weight of actual nodal position/value in computation of mesh
      // velocity from positional/value time stepper
      const double pos_time_weight =
        node_pt(0)->position_time_stepper_pt()->weight(1, 0);
      const double val_time_weight =
        node_pt(0)->time_stepper_pt()->weight(1, 0);

      // Get the user-defined body force terms
      Vector<double> body_force(DIM);
      this->get_body_force_nst(time, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      const double source = this->get_source_nst(time, ipt, interpolated_x);

      // Get gradient of body force function
      DenseMatrix<double> d_body_force_dx(DIM, DIM, 0.0);
      this->get_body_force_gradient_nst(
        time, ipt, s, interpolated_x, d_body_force_dx);

      // Get gradient of source function
      Vector<double> source_gradient(DIM, 0.0);
      this->get_source_gradient_nst(time, ipt, interpolated_x, source_gradient);


      // Assemble shape derivatives
      //---------------------------

      // MOMENTUM EQUATIONS
      // ------------------

      // Number of master nodes and storage for the weight of the shape function
      unsigned n_master = 1;
      double hang_weight = 1.0;

      // Loop over the test functions
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

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Loop over coordinate directions
          for (unsigned i = 0; i < DIM; i++)
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

            // IF it's not a boundary condition
            if (local_eqn >= 0)
            {
              // Loop over coordinate directions
              for (unsigned p = 0; p < DIM; p++)
              {
                // Loop over shape controlling nodes
                for (unsigned q = 0; q < n_shape_controlling_node; q++)
                {
                  // Residual x deriv of Jacobian
                  // ----------------------------

                  // Add the user-defined body force terms
                  double sum = body_force[i] * testf[l];

                  // Add the gravitational body force term
                  sum += scaled_re_inv_fr * testf[l] * G[i];

                  // Add the pressure gradient term
                  sum += interpolated_p * dtestfdx(l, i);

                  // Add in the stress tensor terms
                  // The viscosity ratio needs to go in here to ensure
                  // continuity of normal stress is satisfied even in flows
                  // with zero pressure gradient!
                  for (unsigned k = 0; k < DIM; k++)
                  {
                    sum -= visc_ratio *
                           (interpolated_dudx(i, k) +
                            this->Gamma[i] * interpolated_dudx(k, i)) *
                           dtestfdx(l, k);
                  }

                  // Add in the inertial terms

                  // du/dt term
                  sum -= scaled_re_st * dudt[i] * testf[l];

                  // Convective terms, including mesh velocity
                  for (unsigned k = 0; k < DIM; k++)
                  {
                    double tmp = scaled_re * interpolated_u[k];
                    if (!this->ALE_is_disabled)
                    {
                      tmp -= scaled_re_st * mesh_velocity[k];
                    }
                    sum -= tmp * interpolated_dudx(i, k) * testf[l];
                  }

                  // Multiply throsugh by deriv of Jacobian and integration
                  // weight
                  dresidual_dnodal_coordinates(local_eqn, p, q) +=
                    sum * dJ_dX(p, q) * w * hang_weight;

                  // Derivative of residual x Jacobian
                  // ---------------------------------

                  // Body force
                  sum = d_body_force_dx(i, p) * psif(q) * testf(l);

                  // Pressure gradient term
                  sum += interpolated_p * d_dtestfdx_dX(p, q, l, i);

                  // Viscous term
                  for (unsigned k = 0; k < DIM; k++)
                  {
                    sum -=
                      visc_ratio * ((interpolated_dudx(i, k) +
                                     this->Gamma[i] * interpolated_dudx(k, i)) *
                                      d_dtestfdx_dX(p, q, l, k) +
                                    (d_dudx_dX(p, q, i, k) +
                                     this->Gamma[i] * d_dudx_dX(p, q, k, i)) *
                                      dtestfdx(l, k));
                  }

                  // Convective terms, including mesh velocity
                  for (unsigned k = 0; k < DIM; k++)
                  {
                    double tmp = scaled_re * interpolated_u[k];
                    if (!this->ALE_is_disabled)
                    {
                      tmp -= scaled_re_st * mesh_velocity[k];
                    }
                    sum -= tmp * d_dudx_dX(p, q, i, k) * testf(l);
                  }
                  if (!this->ALE_is_disabled)
                  {
                    sum += scaled_re_st * pos_time_weight * psif(q) *
                           interpolated_dudx(i, p) * testf(l);
                  }

                  // Multiply through by Jacobian and integration weight
                  dresidual_dnodal_coordinates(local_eqn, p, q) +=
                    sum * J * w * hang_weight;

                } // End of loop over shape controlling nodes q
              } // End of loop over coordinate directions p


              // Derivs w.r.t. to nodal velocities
              // ---------------------------------
              if (element_has_node_with_aux_node_update_fct)
              {
                // Loop over local nodes
                for (unsigned q_local = 0; q_local < n_node; q_local++)
                {
                  // Number of master nodes and storage for the weight of
                  // the shape function
                  unsigned n_master2 = 1;
                  double hang_weight2 = 1.0;
                  HangInfo* hang_info2_pt = 0;

                  // Local boolean to indicate whether the node is hanging
                  bool is_node_hanging2 = node_pt(q_local)->is_hanging();

                  Node* actual_shape_controlling_node_pt = node_pt(q_local);

                  // If the node is hanging
                  if (is_node_hanging2)
                  {
                    hang_info2_pt = node_pt(q_local)->hanging_pt();

                    // Read out number of master nodes from hanging data
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise the node is its own master
                  else
                  {
                    n_master2 = 1;
                  }

                  // Loop over the master nodes
                  for (unsigned mm = 0; mm < n_master2; mm++)
                  {
                    if (is_node_hanging2)
                    {
                      actual_shape_controlling_node_pt =
                        hang_info2_pt->master_node_pt(mm);
                      hang_weight2 = hang_info2_pt->master_weight(mm);
                    }

                    // Find the corresponding number
                    unsigned q = local_shape_controlling_node_lookup
                      [actual_shape_controlling_node_pt];

                    // Loop over coordinate directions
                    for (unsigned p = 0; p < DIM; p++)
                    {
                      double sum = -visc_ratio * this->Gamma[i] *
                                     dpsifdx(q_local, i) * dtestfdx(l, p) -
                                   scaled_re * psif(q_local) *
                                     interpolated_dudx(i, p) * testf(l);
                      if (i == p)
                      {
                        sum -= scaled_re_st * val_time_weight * psif(q_local) *
                               testf(l);
                        for (unsigned k = 0; k < DIM; k++)
                        {
                          sum -=
                            visc_ratio * dpsifdx(q_local, k) * dtestfdx(l, k);
                          double tmp = scaled_re * interpolated_u[k];
                          if (!this->ALE_is_disabled)
                          {
                            tmp -= scaled_re_st * mesh_velocity[k];
                          }
                          sum -= tmp * dpsifdx(q_local, k) * testf(l);
                        }
                      }

                      dresidual_dnodal_coordinates(local_eqn, p, q) +=
                        sum * d_U_dX(p, q) * J * w * hang_weight * hang_weight2;
                    }
                  } // End of loop over master nodes
                } // End of loop over local nodes
              } // End of if(element_has_node_with_aux_node_update_fct)


            } // local_eqn>=0
          }
        }
      } // End of loop over test functions


      // CONTINUITY EQUATION
      // -------------------

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
            // Loop over coordinate directions
            for (unsigned p = 0; p < DIM; p++)
            {
              // Loop over nodes
              for (unsigned q = 0; q < n_shape_controlling_node; q++)
              {
                // Residual x deriv of Jacobian
                //-----------------------------

                // Source term
                double aux = -source;

                // Loop over velocity components
                for (unsigned k = 0; k < DIM; k++)
                {
                  aux += interpolated_dudx(k, k);
                }

                // Multiply through by deriv of Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  aux * dJ_dX(p, q) * testp[l] * w * hang_weight;


                // Derivative of residual x Jacobian
                // ---------------------------------

                // Loop over velocity components
                aux = -source_gradient[p] * psif(q);
                for (unsigned k = 0; k < DIM; k++)
                {
                  aux += d_dudx_dX(p, q, k, k);
                }
                // Multiply through by Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  aux * testp[l] * J * w * hang_weight;
              }
            }


            // Derivs w.r.t. to nodal velocities
            // ---------------------------------
            if (element_has_node_with_aux_node_update_fct)
            {
              // Loop over local nodes
              for (unsigned q_local = 0; q_local < n_node; q_local++)
              {
                // Number of master nodes and storage for the weight of
                // the shape function
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                HangInfo* hang_info2_pt = 0;

                // Local boolean to indicate whether the node is hanging
                bool is_node_hanging2 = node_pt(q_local)->is_hanging();

                Node* actual_shape_controlling_node_pt = node_pt(q_local);

                // If the node is hanging
                if (is_node_hanging2)
                {
                  hang_info2_pt = node_pt(q_local)->hanging_pt();

                  // Read out number of master nodes from hanging data
                  n_master2 = hang_info2_pt->nmaster();
                }
                // Otherwise the node is its own master
                else
                {
                  n_master2 = 1;
                }

                // Loop over the master nodes
                for (unsigned mm = 0; mm < n_master2; mm++)
                {
                  if (is_node_hanging2)
                  {
                    actual_shape_controlling_node_pt =
                      hang_info2_pt->master_node_pt(mm);
                    hang_weight2 = hang_info2_pt->master_weight(mm);
                  }

                  // Find the corresponding number
                  unsigned q = local_shape_controlling_node_lookup
                    [actual_shape_controlling_node_pt];

                  // Loop over coordinate directions
                  for (unsigned p = 0; p < DIM; p++)
                  {
                    double aux = dpsifdx(q_local, p) * testp(l);
                    dresidual_dnodal_coordinates(local_eqn, p, q) +=
                      aux * d_U_dX(p, q) * J * w * hang_weight * hang_weight2;
                  }
                } // End of loop over (mm) master nodes
              } // End of loop over local nodes q_local
            } // End of if(element_has_node_with_aux_node_update_fct)
          } // End of if it's not a boundary condition
        } // End of loop over (m) master nodes
      } // End of loop over shape functions for continuity eqn

    } // End of loop over integration points
  }


  //====================================================================
  /// Force build of templates
  //====================================================================
  template class RefineableSpaceTimeNavierStokesEquations<2>;
  template class RefineableQTaylorHoodSpaceTimeElement<2>;
} // namespace oomph
