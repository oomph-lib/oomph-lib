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
// Non-inline functions for linearised axisymmetric Navier-Stokes elements

// oomph-lib includes
#include "linearised_navier_stokes_elements.h"

namespace oomph
{
  //=======================================================================
  /// Linearised axisymmetric Navier--Stokes equations static data
  //=======================================================================

  // Use the stress-divergence form by default (Gamma=1)
  Vector<double> LinearisedNavierStokesEquations::Gamma(2, 1.0);

  // "Magic" number to indicate pressure is not stored at node
  int LinearisedNavierStokesEquations::Pressure_not_stored_at_node = -100;

  // Physical constants default to zero
  double LinearisedNavierStokesEquations::Default_Physical_Constant_Value = 0.0;

  // Density/viscosity ratios default to one
  double LinearisedNavierStokesEquations::Default_Physical_Ratio_Value = 1.0;


  //=======================================================================
  /// Output function in tecplot format: Velocities only
  /// r, z, U^C, U^S, V^C, V^S, W^C, W^S
  /// at specified previous timestep (t=0 present; t>0 previous timestep).
  /// Specified number of plot points in each coordinate direction.
  //=======================================================================
  void LinearisedNavierStokesEquations::output_veloc(std::ostream& outfile,
                                                     const unsigned& nplot,
                                                     const unsigned& t)
  {
    // Determine number of nodes in element
    const unsigned n_node = nnode();

    // Provide storage for local shape functions
    Shape psi(n_node);

    // Provide storage for vectors of local coordinates and
    // global coordinates and velocities
    Vector<double> s(DIM);
    Vector<double> interpolated_x(DIM);
    Vector<double> interpolated_u(2 * DIM);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Determine number of plot points
    const unsigned n_plot_points = nplot_points(nplot);

    // Loop over plot points
    for (unsigned iplot = 0; iplot < n_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get shape functions
      shape(s, psi);

      // Loop over coordinate directions
      for (unsigned i = 0; i < DIM; i++)
      {
        // Initialise global coordinate
        interpolated_x[i] = 0.0;

        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          interpolated_x[i] += nodal_position(t, l, i) * psi[l];
        }
      }

      // Loop over the velocity components
      for (unsigned i = 0; i < (2 * DIM); i++)
      {
        // Get the index at which the velocity is stored
        const unsigned u_nodal_index = u_index_linearised_nst(i);

        // Initialise velocity
        interpolated_u[i] = 0.0;

        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          interpolated_u[i] += nodal_value(t, l, u_nodal_index) * psi[l];
        }
      }

      // Output global coordinates to file
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x[i] << " ";
      }

      // Output velocities to file
      for (unsigned i = 0; i < (2 * DIM); i++)
      {
        outfile << interpolated_u[i] << " ";
      }

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);

  } // End of output_veloc


  //=======================================================================
  /// Output function in tecplot format:
  /// r, z, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //=======================================================================
  void LinearisedNavierStokesEquations::output(std::ostream& outfile,
                                               const unsigned& nplot)
  {
    // Provide storage for vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Determine number of plot points
    const unsigned n_plot_points = nplot_points(nplot);

    // Loop over plot points
    for (unsigned iplot = 0; iplot < n_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Output global coordinates to file
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      //  Output velocities (and normalisation) to file
      for (unsigned i = 0; i < (4 * DIM); i++)
      {
        outfile << interpolated_u_linearised_nst(s, i) << " ";
      }

      // Output pressure to file
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_p_linearised_nst(s, i) << " ";
      }

      outfile << std::endl;
    }
    outfile << std::endl;

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);

  } // End of output


  //=======================================================================
  /// Output function in tecplot format:
  /// r, z, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S
  /// Specified number of plot points in each coordinate direction.
  //=======================================================================
  void LinearisedNavierStokesEquations::output(FILE* file_pt,
                                               const unsigned& nplot)
  {
    // Provide storage for vector of local coordinates
    Vector<double> s(2);

    // Tecplot header info
    fprintf(file_pt, "%s ", tecplot_zone_string(nplot).c_str());

    // Determine number of plot points
    const unsigned n_plot_points = nplot_points(nplot);

    // Loop over plot points
    for (unsigned iplot = 0; iplot < n_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Output global coordinates to file
      for (unsigned i = 0; i < 2; i++)
      {
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }

      //  Output velocities to file
      for (unsigned i = 0; i < (2 * DIM); i++)
      {
        fprintf(file_pt, "%g ", interpolated_u_linearised_nst(s, i));
      }

      // Output pressure to file
      for (unsigned i = 0; i < 2; i++)
      {
        fprintf(file_pt, "%g ", interpolated_p_linearised_nst(s, i));
      }

      fprintf(file_pt, "\n");
    }

    fprintf(file_pt, "\n");

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, nplot);

  } // End of output


  //=======================================================================
  /// Get strain-rate tensor: \f$ e_{ij} \f$  where
  /// \f$ i,j = r,z,\theta \f$ (in that order). We evaluate this tensor
  /// at a value of theta such that the product of theta and the azimuthal
  /// mode number (k) gives \f$ \pi/4 \f$. Therefore
  /// \f$ \cos(k \theta) = \sin(k \theta) = 1/\sqrt{2} \f$.
  //=======================================================================
  void LinearisedNavierStokesEquations::strain_rate(
    const Vector<double>& s,
    DenseMatrix<double>& strainrate,
    const unsigned& real) const
  {
#ifdef PARANOID
    if ((strainrate.ncol() != DIM) || (strainrate.nrow() != DIM))
    {
      std::ostringstream error_message;
      error_message << "The strain rate has incorrect dimensions "
                    << strainrate.ncol() << " x " << strainrate.nrow()
                    << " Not " << DIM << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // The question is what to do about the real and imaginary parts
    // The answer is that we know that the "real" velocity is
    // U_r cos(omega t) - U_i sin(omega t) and we can choose omega t
    // to be 7pi/4 so that cos = 1/sqrt(2) and sin = -1/sqrt(2)
    // to get equal contributions
    double cosomegat = 1.0 / sqrt(2.0);
    double sinomegat = cosomegat;

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psif(n_node);
    DShape dpsifdx(n_node, DIM);

    // Call the derivatives of the shape functions
    dshape_eulerian(s, psif, dpsifdx);

    // Velocity gradient matrix
    Vector<Vector<std::complex<double>>> dudx(DIM);
    // Initialise to zero
    for (unsigned i = 0; i < DIM; ++i)
    {
      dudx[i].resize(DIM);
      for (unsigned j = 0; j < DIM; ++j)
      {
        dudx[i][j].real(0.0);
        dudx[i][j].imag(0.0);
      }
    }

    // Get the nodal indices at which the velocity is stored
    unsigned n_veloc = 2 * DIM;
    unsigned u_nodal_index[n_veloc];
    for (unsigned i = 0; i < n_veloc; ++i)
    {
      u_nodal_index[i] = u_index_linearised_nst(i);
    }

    // Loop over the element's nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over the DIM complex velocity components
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get the value
        const std::complex<double> u_value(
          this->raw_nodal_value(l, u_nodal_index[2 * i + 0]),
          this->raw_nodal_value(l, u_nodal_index[2 * i + 1]));

        // Loop over two coordinate directions (for derivatives)
        for (unsigned j = 0; j < DIM; j++)
        {
          dudx[i][j] += u_value * dpsifdx(l, j);
        }
      }
    }

    // Take the real part of the velocity gradient matrix
    DenseMatrix<double> real_dudx(DIM);
    if (real == 0)
    {
      cosomegat = 1.0;
      sinomegat = 0.0;
    }
    else
    {
      cosomegat = 0.0;
      sinomegat = 1.0;
    }

    for (unsigned i = 0; i < DIM; ++i)
    {
      for (unsigned j = 0; j < DIM; ++j)
      {
        real_dudx(i, j) =
          dudx[i][j].real() * cosomegat + dudx[i][j].imag() * sinomegat;
      }
    }


    // Now set the strain rate
    // Loop over veclocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Loop over derivative directions
      for (unsigned j = 0; j < DIM; j++)
      {
        strainrate(i, j) = 0.5 * (real_dudx(i, j) + real_dudx(j, i));
      }
    }

  } // End of strain_rate


  //=======================================================================
  /// Compute the residuals for the linearised axisymmetric Navier--Stokes
  /// equations; flag=1(or 0): do (or don't) compute the Jacobian as well.
  //=======================================================================
  void LinearisedNavierStokesEquations::
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

    // Integers used to store the local equation and unknown numbers
    int local_eqn = 0;

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
          interpolated_x[i] += this->raw_nodal_position(l, i) * psif_;
        }

        // Loop over the DIM complex velocity components
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the value
          const std::complex<double> u_value(
            this->raw_nodal_value(l, u_nodal_index[2 * i + 0]),
            this->raw_nodal_value(l, u_nodal_index[2 * i + 1]));

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
            this->raw_nodal_value(l, u_nodal_index[2 * (DIM + i)]),
            this->raw_nodal_value(l, u_nodal_index[2 * (DIM + i) + 1]));
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

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
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

          /*IF it's not a boundary condition*/
          // Here assume that we're only going to pin entire complex
          // number or not
          local_eqn = nodal_local_eqn(l, u_nodal_index[2 * i]);
          if (local_eqn >= 0)
          {
            residuals[local_eqn] += residual_contribution.real();
          }

          local_eqn = nodal_local_eqn(l, u_nodal_index[2 * i + 1]);
          if (local_eqn >= 0)
          {
            residuals[local_eqn] += residual_contribution.imag();
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
      } // End of loop over shape functions


      // CONTINUITY EQUATION
      //-------------------

      // Loop over the shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Assemble the residuals
        std::complex<double> residual_contribution = interpolated_dudx[0][0];
        for (unsigned k = 1; k < DIM; ++k)
        {
          residual_contribution += interpolated_dudx[k][k];
        }

        local_eqn = p_local_eqn(l, 0);
        // If not a boundary conditions
        if (local_eqn >= 0)
        {
          residuals[local_eqn] += residual_contribution.real() * testp[l] * W;
        }

        local_eqn = p_local_eqn(l, 1);
        // If not a boundary conditions
        if (local_eqn >= 0)
        {
          residuals[local_eqn] += residual_contribution.imag() * testp[l] * W;
        }

      } // End of loop over l

      // Normalisation condition
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


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


  /// Linearised axisymmetric Crouzeix-Raviart elements
  /// -------------------------------------------------


  //=======================================================================
  /// Set the data for the number of variables at each node
  //=======================================================================
  const unsigned LinearisedQCrouzeixRaviartElement::Initial_Nvalue[9] = {
    8, 8, 8, 8, 8, 8, 8, 8, 8};


  //========================================================================
  /// Number of values (pinned or dofs) required at node n
  //========================================================================
  unsigned LinearisedQCrouzeixRaviartElement::required_nvalue(
    const unsigned& n) const
  {
    return Initial_Nvalue[n];
  }


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


  /// Linearised axisymmetric Taylor-Hood elements
  /// --------------------------------------------


  //=======================================================================
  /// Set the data for the number of variables at each node
  //=======================================================================
  const unsigned LinearisedQTaylorHoodElement::Initial_Nvalue[9] = {
    12, 8, 12, 8, 8, 8, 12, 8, 12};


  //=======================================================================
  /// Set the data for the pressure conversion array
  //=======================================================================
  const unsigned LinearisedQTaylorHoodElement::Pconv[4] = {0, 2, 6, 8};


} // namespace oomph
