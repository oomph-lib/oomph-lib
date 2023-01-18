// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Non-inline functions for the refineable linearised axisymmetric
// Navier-Stokes elements

// oomph-lib includes
#include "refineable_linearised_navier_stokes_elements.h"

namespace oomph
{
  //=======================================================================
  /// Compute the residuals for the refineable linearised axisymmetric
  /// Navier--Stokes equations; flag=1(or 0): do (or don't) compute the
  /// Jacobian as well.
  //=======================================================================
  void RefineableLinearisedNavierStokesEquations::
    fill_in_generic_residual_contribution_linearised_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Get the time from the first node in the element
    const double time = this->node_pt(0)->time_stepper_pt()->time();

    // Determine number of nodes in the element
    const unsigned n_node = nnode();

    // Determine how many pressure values there are associated with
    // a single pressure component
    const unsigned n_pres = npres_linearised_nst();

    const unsigned n_veloc = 4 * DIM;

    // Get the nodal indices at which the velocity is stored
    unsigned u_nodal_index[n_veloc];
    for (unsigned i = 0; i < n_veloc; ++i)
    {
      u_nodal_index[i] = u_index_linearised_nst(i);
    }

    // Which nodal values represent the two pressure components?
    // (Negative if pressure is not based on nodal interpolation).
    Vector<int> p_index(2);
    for (unsigned i = 0; i < 2; i++)
    {
      p_index[i] = this->p_index_linearised_nst(i);
    }

    // Local array of booleans that are true if the l-th pressure value is
    // hanging (avoid repeated virtual function calls)
    bool pressure_dof_is_hanging[n_pres];

    // If the pressure is stored at a node
    if (p_index[0] >= 0)
    {
      // Read out whether the pressure is hanging
      for (unsigned l = 0; l < n_pres; ++l)
      {
        pressure_dof_is_hanging[l] =
          pressure_node_pt(l)->is_hanging(p_index[0]);
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


    // Set up memory for the fluid shape and test functions
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, DIM), dtestfdx(n_node, DIM);

    // Set up memory for the pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Determine number of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Set up memory for the vector to hold local coordinates (two dimensions)
    Vector<double> s(DIM);

    // Get physical variables from the element
    // (Reynolds number must be multiplied by the density ratio)
    const double scaled_re = re() * density_ratio();
    const double scaled_re_st = re_st() * density_ratio();
    const double visc_ratio = viscosity_ratio();

    const double eval_real = lambda();
    const double eval_imag = omega();

    const std::complex<double> eigenvalue(eval_real, eval_imag);

    // Integers used to store the local equation numbers
    int local_eqn = 0;

    // Local storage for pointers to hang info objects
    HangInfo* hang_info_pt = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of the local coordinates s
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      const double w = integral_pt()->weight(ipt);

      // Calculate the fluid shape and test functions, and their derivatives
      // w.r.t. the global coordinates
      const double J = dshape_and_dtest_eulerian_at_knot_linearised_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Calculate the pressure shape and test functions
      pshape_linearised_nst(s, psip, testp);

      // Premultiply the weights and the Jacobian of the mapping between
      // local and global coordinates
      const double W = w * J;

      // Allocate storage for the position and the derivative of the
      // mesh positions w.r.t. time
      Vector<double> interpolated_x(DIM, 0.0);
      // Vector<double> mesh_velocity(2,0.0);

      // Allocate storage for the velocity components (six of these)
      // and their derivatives w.r.t. time
      Vector<std::complex<double>> interpolated_u(DIM);
      // Vector<double> dudt(6,0.0);
      // Allocate storage for the eigen function normalisation
      Vector<std::complex<double>> interpolated_u_normalisation(DIM);
      for (unsigned i = 0; i < DIM; ++i)
      {
        interpolated_u[i].real(0.0);
        interpolated_u[i].imag(0.0);
        interpolated_u_normalisation[i].real(0.0);
        interpolated_u_normalisation[i].imag(0.0);
      }

      // Allocate storage for the pressure components (two of these
      std::complex<double> interpolated_p(0.0, 0.0);
      std::complex<double> interpolated_p_normalisation(0.0, 0.0);

      // Allocate storage for the derivatives of the velocity components
      // w.r.t. global coordinates (r and z)
      Vector<Vector<std::complex<double>>> interpolated_dudx(DIM);
      for (unsigned i = 0; i < DIM; ++i)
      {
        interpolated_dudx[i].resize(DIM);
        for (unsigned j = 0; j < DIM; ++j)
        {
          interpolated_dudx[i][j].real(0.0);
          interpolated_dudx[i][j].imag(0.0);
        }
      }

      // Calculate pressure at the integration point
      // -------------------------------------------

      // Loop over pressure degrees of freedom (associated with a single
      // pressure component) in the element
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Cache the shape function
        const double psip_ = psip(l);

        // Get the complex nodal pressure values
        const std::complex<double> p_value(this->p_linearised_nst(l, 0),
                                           this->p_linearised_nst(l, 1));

        // Add contribution
        interpolated_p += p_value * psip_;

        // Repeat for normalisation
        const std::complex<double> p_norm_value(this->p_linearised_nst(l, 2),
                                                this->p_linearised_nst(l, 3));
        interpolated_p_normalisation += p_norm_value * psip_;
      }
      // End of loop over the pressure degrees of freedom in the element

      // Calculate velocities and their derivatives at the integration point
      // -------------------------------------------------------------------

      // Loop over the element's nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the shape function
        const double psif_ = psif(l);

        // Loop over the DIM coordinate directions
        for (unsigned i = 0; i < DIM; i++)
        {
          interpolated_x[i] += this->nodal_position(l, i) * psif_;
        }

        // Loop over the DIM complex velocity components
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the value
          const std::complex<double> u_value(
            this->nodal_value(l, u_nodal_index[2 * i + 0]),
            this->nodal_value(l, u_nodal_index[2 * i + 1]));

          // Add contribution
          interpolated_u[i] += u_value * psif_;

          // Add contribution to dudt
          // dudt[i] += du_dt_linearised_nst(l,i)*psif_;

          // Loop over two coordinate directions (for derivatives)
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_dudx[i][j] += u_value * dpsifdx(l, j);
          }

          // Interpolate the normalisation function
          const std::complex<double> normalisation_value(
            this->nodal_value(l, u_nodal_index[2 * (DIM + i)]),
            this->nodal_value(l, u_nodal_index[2 * (DIM + i) + 1]));
          interpolated_u_normalisation[i] += normalisation_value * psif_;
        }
      } // End of loop over the element's nodes

      // Get the mesh velocity if ALE is enabled
      /*if(!ALE_is_disabled)
       {
        // Loop over the element's nodes
        for(unsigned l=0;l<n_node;l++)
         {
          // Loop over the two coordinate directions
          for(unsigned i=0;i<2;i++)
           {
            mesh_velocity[i] += this->raw_dnodal_position_dt(l,i)*psif(l);
           }
         }
         }*/

      // Get velocities and their derivatives from base flow problem
      // -----------------------------------------------------------

      // Allocate storage for the velocity components of the base state
      // solution (initialise to zero)
      Vector<double> base_flow_u(DIM, 0.0);

      // Get the user-defined base state solution velocity components
      get_base_flow_u(time, ipt, interpolated_x, base_flow_u);

      // Allocate storage for the derivatives of the base state solution's
      // velocity components w.r.t. global coordinate (r and z)
      // N.B. the derivatives of the base flow components w.r.t. the
      // azimuthal coordinate direction (theta) are always zero since the
      // base flow is axisymmetric
      DenseMatrix<double> base_flow_dudx(DIM, DIM, 0.0);

      // Get the derivatives of the user-defined base state solution
      // velocity components w.r.t. global coordinates
      get_base_flow_dudx(time, ipt, interpolated_x, base_flow_dudx);


      // MOMENTUM EQUATIONS
      //------------------

      // Number of master nodes
      unsigned n_master = 1;

      // Storage for the weight of the shape function
      double hang_weight = 1.0;

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Local boolean to indicate whether the node is hanging
        bool is_node_hanging = node_pt(l)->is_hanging();

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
          // Loop over the velocity components
          for (unsigned i = 0; i < DIM; i++)
          {
            // Assemble the residuals
            // Time dependent term
            std::complex<double> residual_contribution =
              -scaled_re_st * eigenvalue * interpolated_u[i] * testf[l] * W;
            // Pressure term
            residual_contribution += interpolated_p * dtestfdx(l, i) * W;
            // Viscous terms
            for (unsigned k = 0; k < DIM; ++k)
            {
              residual_contribution -=
                visc_ratio *
                (interpolated_dudx[i][k] + Gamma[i] * interpolated_dudx[k][i]) *
                dtestfdx(l, k) * W;
            }

            // Advective terms
            for (unsigned k = 0; k < DIM; ++k)
            {
              residual_contribution -=
                scaled_re *
                (base_flow_u[k] * interpolated_dudx[i][k] +
                 interpolated_u[k] * base_flow_dudx(i, k)) *
                testf[l] * W;
            }

            // Now separate real and imaginary parts

            if (is_node_hanging)
            {
              local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                               u_nodal_index[2 * i]);
              hang_weight = hang_info_pt->master_weight(m);
            }
            // If node is not hanging number or not
            else
            {
              local_eqn = nodal_local_eqn(l, u_nodal_index[2 * i]);
              hang_weight = 1.0;
            }

            if (local_eqn >= 0)
            {
              residuals[local_eqn] +=
                residual_contribution.real() * hang_weight;
            }


            if (is_node_hanging)
            {
              local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                               u_nodal_index[2 * i + 1]);
              hang_weight = hang_info_pt->master_weight(m);
            }
            // If node is not hanging number or not
            else
            {
              local_eqn = nodal_local_eqn(l, u_nodal_index[2 * i + 1]);
              hang_weight = 1.0;
            }
            if (local_eqn >= 0)
            {
              residuals[local_eqn] +=
                residual_contribution.imag() * hang_weight;
            }


            // CALCULATE THE JACOBIAN
            /*if(flag)
              {
              //Loop over the velocity shape functions again
              for(unsigned l2=0;l2<n_node;l2++)
              {
              //Loop over the velocity components again
              for(unsigned i2=0;i2<DIM;i2++)
              {
              //If at a non-zero degree of freedom add in the entry
              local_unknown = nodal_local_eqn(l2,u_nodal_index[i2]);
              if(local_unknown >= 0)
              {
              //Add contribution to Elemental Matrix
              jacobian(local_eqn,local_unknown)
              -= visc_ratio*Gamma[i]*dpsifdx(l2,i)*dtestfdx(l,i2)*W;

              //Extra component if i2 = i
              if(i2 == i)
              {
              //Loop over velocity components
              for(unsigned k=0;k<DIM;k++)
              {
              jacobian(local_eqn,local_unknown)
              -= visc_ratio*dpsifdx(l2,k)*dtestfdx(l,k)*W;
              }
              }

              //Now add in the inertial terms
              jacobian(local_eqn,local_unknown)
              -= scaled_re*psif[l2]*interpolated_dudx(i,i2)*testf[l]*W;

              //Extra component if i2=i
              if(i2 == i)
              {
              //Add the mass matrix term (only diagonal entries)
              //Note that this is positive because the mass matrix
              //is taken to the other side of the equation when
              //formulating the generalised eigenproblem.
              if(flag==2)
              {
              mass_matrix(local_eqn,local_unknown) +=
              scaled_re_st*psif[l2]*testf[l]*W;
              }

              //du/dt term
              jacobian(local_eqn,local_unknown)
              -= scaled_re_st*
              node_pt(l2)->time_stepper_pt()->weight(1,0)*
              psif[l2]*testf[l]*W;

              //Loop over the velocity components
              for(unsigned k=0;k<DIM;k++)
              {
              double tmp=scaled_re*interpolated_u[k];
              if (!ALE_is_disabled) tmp-=scaled_re_st*mesh_velocity[k];
              jacobian(local_eqn,local_unknown) -=
              tmp*dpsifdx(l2,k)*testf[l]*W;
              }
              }

              }
              }
              }

              //Now loop over pressure shape functions
              //This is the contribution from pressure gradient
              for(unsigned l2=0;l2<n_pres;l2++)
              {
              //If we are at a non-zero degree of freedom in the entry
              local_unknown = p_local_eqn(l2);
              if(local_unknown >= 0)
              {
              jacobian(local_eqn,local_unknown)
              += psip[l2]*dtestfdx(l,i)*W;
              }
              }
              } //End of Jacobian calculation

              }*/ //End of if not boundary condition statement

          } // End of loop over dimension
        } // End of loop over master nodes
      } // End of loop over shape functions


      // CONTINUITY EQUATION
      //-------------------

      // Loop over the shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        if (pressure_dof_is_hanging[l])
        {
          hang_info_pt = this->pressure_node_pt(l)->hanging_pt(p_index[0]);
          n_master = hang_info_pt->nmaster();
        }
        else
        {
          n_master = 1;
        }

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; ++m)
        {
          // Assemble the residuals
          std::complex<double> residual_contribution = interpolated_dudx[0][0];
          for (unsigned k = 1; k < DIM; ++k)
          {
            residual_contribution += interpolated_dudx[k][k];
          }

          if (pressure_dof_is_hanging[l])
          {
            local_eqn =
              this->local_hang_eqn(hang_info_pt->master_node_pt(m), p_index[0]);
            hang_weight = hang_info_pt->master_weight(m);
          }
          else
          {
            local_eqn = this->p_local_eqn(l, 0);
            hang_weight = 1.0;
          }

          // If not a boundary conditions
          if (local_eqn >= 0)
          {
            residuals[local_eqn] +=
              residual_contribution.real() * testp[l] * W * hang_weight;
          }

          if (pressure_dof_is_hanging[l])
          {
            local_eqn =
              this->local_hang_eqn(hang_info_pt->master_node_pt(m), p_index[1]);
            hang_weight = hang_info_pt->master_weight(m);
          }
          else
          {
            local_eqn = this->p_local_eqn(l, 1);
            hang_weight = 1.0;
          }

          // If not a boundary conditions
          if (local_eqn >= 0)
          {
            residuals[local_eqn] +=
              residual_contribution.imag() * testp[l] * W * hang_weight;
          }

        } // End of loop over master nodes
      } // End of loop over l

      // Normalisation condition. Leave this alone because there is
      // no test function involved.
      std::complex<double> residual_contribution =
        interpolated_p_normalisation * interpolated_p;
      for (unsigned k = 0; k < DIM; ++k)
      {
        residual_contribution +=
          interpolated_u_normalisation[k] * interpolated_u[k];
      }

      local_eqn = this->eigenvalue_local_eqn(0);
      if (local_eqn >= 0)
      {
        residuals[local_eqn] += residual_contribution.real() * W;
      }

      local_eqn = this->eigenvalue_local_eqn(1);
      if (local_eqn >= 0)
      {
        residuals[local_eqn] += residual_contribution.imag() * W;
      }

    } // End of loop over the integration points

  } // End of fill_in_generic_residual_contribution_linearised_nst

} // End of namespace oomph
