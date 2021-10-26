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
// Non-inline functions for linearised axisymmetric Navier-Stokes elements

// oomph-lib includes
#include "linearised_axisym_navier_stokes_elements.h"

namespace oomph
{
  //=======================================================================
  /// Linearised axisymmetric Navier--Stokes equations static data
  //=======================================================================

  // Use the stress-divergence form by default (Gamma=1)
  Vector<double> LinearisedAxisymmetricNavierStokesEquations::Gamma(2, 1.0);

  // "Magic" number to indicate pressure is not stored at node
  int LinearisedAxisymmetricNavierStokesEquations::Pressure_not_stored_at_node =
    -100;

  // Physical constants default to zero
  double LinearisedAxisymmetricNavierStokesEquations::
    Default_Physical_Constant_Value = 0.0;

  // Azimuthal mode number defaults to zero
  int LinearisedAxisymmetricNavierStokesEquations::
    Default_Azimuthal_Mode_Number_Value = 0;

  // Density/viscosity ratios default to one
  double
    LinearisedAxisymmetricNavierStokesEquations::Default_Physical_Ratio_Value =
      1.0;


  //=======================================================================
  /// Output function in tecplot format: Velocities only
  /// r, z, U^C, U^S, V^C, V^S, W^C, W^S
  /// at specified previous timestep (t=0 present; t>0 previous timestep).
  /// Specified number of plot points in each coordinate direction.
  //=======================================================================
  void LinearisedAxisymmetricNavierStokesEquations::output_veloc(
    std::ostream& outfile, const unsigned& nplot, const unsigned& t)
  {
    // Determine number of nodes in element
    const unsigned n_node = nnode();

    // Provide storage for local shape functions
    Shape psi(n_node);

    // Provide storage for vectors of local coordinates and
    // global coordinates and velocities
    Vector<double> s(2);
    Vector<double> interpolated_x(2);
    Vector<double> interpolated_u(6);

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
      for (unsigned i = 0; i < 2; i++)
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
      for (unsigned i = 0; i < 6; i++)
      {
        // Get the index at which the velocity is stored
        const unsigned u_nodal_index = u_index_linearised_axi_nst(i);

        // Initialise velocity
        interpolated_u[i] = 0.0;

        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          interpolated_u[i] += nodal_value(t, l, u_nodal_index) * psi[l];
        }
      }

      // Output global coordinates to file
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_x[i] << " ";
      }

      // Output velocities to file
      for (unsigned i = 0; i < 6; i++)
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
  void LinearisedAxisymmetricNavierStokesEquations::output(
    std::ostream& outfile, const unsigned& nplot)
  {
    // Provide storage for vector of local coordinates
    Vector<double> s(2);

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
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      //  Output velocities to file
      for (unsigned i = 0; i < 6; i++)
      {
        outfile << interpolated_u_linearised_axi_nst(s, i) << " ";
      }

      // Output pressure to file
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_p_linearised_axi_nst(s, i) << " ";
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
  void LinearisedAxisymmetricNavierStokesEquations::output(
    FILE* file_pt, const unsigned& nplot)
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
      for (unsigned i = 0; i < 6; i++)
      {
        fprintf(file_pt, "%g ", interpolated_u_linearised_axi_nst(s, i));
      }

      // Output pressure to file
      for (unsigned i = 0; i < 2; i++)
      {
        fprintf(file_pt, "%g ", interpolated_p_linearised_axi_nst(s, i));
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
  void LinearisedAxisymmetricNavierStokesEquations::strain_rate(
    const Vector<double>& s, DenseMatrix<double>& strainrate) const
  {
#ifdef PARANOID
    if ((strainrate.ncol() != 3) || (strainrate.nrow() != 3))
    {
      std::ostringstream error_message;
      error_message << "The strain rate has incorrect dimensions "
                    << strainrate.ncol() << " x " << strainrate.nrow()
                    << " Not 3" << std::endl;

      throw OomphLibError(
        error_message.str(),
        "LinearisedAxisymmetricNavierStokeEquations::strain_rate()",
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Determine number of nodes in element
    const unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node, 2);

    // Call the derivatives of the shape functions
    dshape_eulerian(s, psi, dpsidx);

    // Radius
    double interpolated_r = 0.0;

    // Velocity components and their derivatives
    double UC = 0.0, US = 0.0;
    double dUCdr = 0.0, dUSdr = 0.0;
    double dUCdz = 0.0, dUSdz = 0.0;
    double WC = 0.0, WS = 0.0;
    double dWCdr = 0.0, dWSdr = 0.0;
    double dWCdz = 0.0, dWSdz = 0.0;
    double VC = 0.0, VS = 0.0;
    double dVCdr = 0.0, dVSdr = 0.0;
    double dVCdz = 0.0, dVSdz = 0.0;

    // Get the local storage for the indices
    unsigned u_nodal_index[6];
    for (unsigned i = 0; i < 6; ++i)
    {
      u_nodal_index[i] = u_index_linearised_axi_nst(i);
    }

    // Loop over nodes to assemble velocities and their derivatives
    // w.r.t. r and z (x_0 and x_1)
    for (unsigned l = 0; l < n_node; l++)
    {
      interpolated_r += nodal_position(l, 0) * psi[l];

      UC += nodal_value(l, u_nodal_index[0]) * psi[l];
      US += nodal_value(l, u_nodal_index[1]) * psi[l];
      WC += nodal_value(l, u_nodal_index[2]) * psi[l];
      WS += nodal_value(l, u_nodal_index[3]) * psi[l];
      VC += nodal_value(l, u_nodal_index[4]) * psi[l];
      VS += nodal_value(l, u_nodal_index[4]) * psi[l];

      dUCdr += nodal_value(l, u_nodal_index[0]) * dpsidx(l, 0);
      dUSdr += nodal_value(l, u_nodal_index[1]) * dpsidx(l, 0);
      dWCdr += nodal_value(l, u_nodal_index[2]) * dpsidx(l, 0);
      dWSdr += nodal_value(l, u_nodal_index[3]) * dpsidx(l, 0);
      dVCdr += nodal_value(l, u_nodal_index[4]) * dpsidx(l, 0);
      dVSdr += nodal_value(l, u_nodal_index[5]) * dpsidx(l, 0);

      dUCdz += nodal_value(l, u_nodal_index[0]) * dpsidx(l, 1);
      dUSdz += nodal_value(l, u_nodal_index[1]) * dpsidx(l, 1);
      dWCdz += nodal_value(l, u_nodal_index[2]) * dpsidx(l, 1);
      dWSdz += nodal_value(l, u_nodal_index[3]) * dpsidx(l, 1);
      dVCdz += nodal_value(l, u_nodal_index[4]) * dpsidx(l, 1);
      dVSdz += nodal_value(l, u_nodal_index[5]) * dpsidx(l, 1);
    }

    // Cache azimuthal mode number
    const int k = this->azimuthal_mode_number();

    // We wish to evaluate the strain-rate tensor at a value of theta
    // such that k*theta = pi/4 radians. That way we pick up equal
    // contributions from the real and imaginary parts of the velocities.
    // sin(pi/4) = cos(pi/4) = 1/sqrt(2)
    const double cosktheta = 1.0 / sqrt(2);
    const double sinktheta = cosktheta;

    // Assemble velocities and their derivatives w.r.t. r, z and theta
    // from real and imaginary parts
    const double ur = UC * cosktheta + US * sinktheta;
    const double utheta = VC * cosktheta + VS * sinktheta;

    const double durdr = dUCdr * cosktheta + dUSdr * sinktheta;
    const double durdz = dUCdz * cosktheta + dUSdz * sinktheta;
    const double durdtheta = k * US * cosktheta - k * UC * sinktheta;

    const double duzdr = dWCdr * cosktheta + dWSdr * sinktheta;
    const double duzdz = dWCdz * cosktheta + dWSdz * sinktheta;
    const double duzdtheta = k * WS * cosktheta - k * WC * sinktheta;

    const double duthetadr = dVCdr * cosktheta + dVSdr * sinktheta;
    const double duthetadz = dVCdz * cosktheta + dVSdz * sinktheta;
    const double duthetadtheta = k * VS * cosktheta - k * VC * sinktheta;

    // Assign strain rates without negative powers of the radius
    // and zero those with:
    strainrate(0, 0) = durdr;
    strainrate(0, 1) = 0.5 * (durdz + duzdr);
    strainrate(1, 0) = strainrate(0, 1);
    strainrate(0, 2) = 0.5 * duthetadr;
    strainrate(2, 0) = strainrate(0, 2);
    strainrate(1, 1) = duzdz;
    strainrate(1, 2) = 0.5 * duthetadz;
    strainrate(2, 1) = strainrate(1, 2);
    strainrate(2, 2) = 0.0;

    // Overwrite the strain rates with negative powers of the radius
    // unless we're at the origin
    if (std::abs(interpolated_r) > 1.0e-16)
    {
      double inverse_radius = 1.0 / interpolated_r;
      strainrate(0, 2) =
        0.5 * (duthetadr + inverse_radius * (durdtheta - utheta));
      strainrate(2, 0) = strainrate(0, 2);
      strainrate(2, 2) = inverse_radius * (ur + duthetadtheta);
      strainrate(1, 2) = 0.5 * (duthetadz + inverse_radius * duzdtheta);
      strainrate(2, 1) = strainrate(1, 2);
    }

  } // End of strain_rate


  //=======================================================================
  /// Compute the residuals for the linearised axisymmetric Navier--Stokes
  /// equations; flag=1(or 0): do (or don't) compute the Jacobian as well.
  //=======================================================================
  void LinearisedAxisymmetricNavierStokesEquations::
    fill_in_generic_residual_contribution_linearised_axi_nst(
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
    const unsigned n_pres = npres_linearised_axi_nst();

    // Get the nodal indices at which the velocity is stored
    unsigned u_nodal_index[6];
    for (unsigned i = 0; i < 6; ++i)
    {
      u_nodal_index[i] = u_index_linearised_axi_nst(i);
    }

    // Set up memory for the fluid shape and test functions
    // Note that there are two spatial dimensions, r and z, in this problem
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, 2), dtestfdx(n_node, 2);

    // Set up memory for the pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Determine number of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Set up memory for the vector to hold local coordinates (two dimensions)
    Vector<double> s(2);

    // Get physical variables from the element
    // (Reynolds number must be multiplied by the density ratio)
    const double scaled_re = re() * density_ratio();
    const double scaled_re_st = re_st() * density_ratio();
    const double visc_ratio = viscosity_ratio();
    const int k = azimuthal_mode_number();

    // Integers used to store the local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of the local coordinates s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      const double w = integral_pt()->weight(ipt);

      // Calculate the fluid shape and test functions, and their derivatives
      // w.r.t. the global coordinates
      const double J = dshape_and_dtest_eulerian_at_knot_linearised_axi_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Calculate the pressure shape and test functions
      pshape_linearised_axi_nst(s, psip, testp);

      // Premultiply the weights and the Jacobian of the mapping between
      // local and global coordinates
      const double W = w * J;

      // Allocate storage for the position and the derivative of the
      // mesh positions w.r.t. time
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> mesh_velocity(2, 0.0);

      // Allocate storage for the velocity components (six of these)
      // and their derivatives w.r.t. time
      Vector<double> interpolated_u(6, 0.0);
      Vector<double> dudt(6, 0.0);

      // Allocate storage for the pressure components (two of these)
      Vector<double> interpolated_p(2, 0.0);

      // Allocate storage for the derivatives of the velocity components
      // w.r.t. global coordinates (r and z)
      DenseMatrix<double> interpolated_dudx(6, 2, 0.0);

      // Calculate pressure at the integration point
      // -------------------------------------------

      // Loop over pressure degrees of freedom (associated with a single
      // pressure component) in the element
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Cache the shape function
        const double psip_ = psip(l);

        // Loop over the two pressure components
        for (unsigned i = 0; i < 2; i++)
        {
          // Get the value
          const double p_value = this->p_linearised_axi_nst(l, i);

          // Add contribution
          interpolated_p[i] += p_value * psip_;
        }
      } // End of loop over the pressure degrees of freedom in the element

      // Calculate velocities and their derivatives at the integration point
      // -------------------------------------------------------------------

      // Loop over the element's nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the shape function
        const double psif_ = psif(l);

        // Loop over the two coordinate directions
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += this->raw_nodal_position(l, i) * psif_;
        }

        // Loop over the six velocity components
        for (unsigned i = 0; i < 6; i++)
        {
          // Get the value
          const double u_value = this->raw_nodal_value(l, u_nodal_index[i]);

          // Add contribution
          interpolated_u[i] += u_value * psif_;

          // Add contribution to dudt
          dudt[i] += du_dt_linearised_axi_nst(l, i) * psif_;

          // Loop over two coordinate directions (for derivatives)
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      } // End of loop over the element's nodes

      // Get the mesh velocity if ALE is enabled
      if (!ALE_is_disabled)
      {
        // Loop over the element's nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the two coordinate directions
          for (unsigned i = 0; i < 2; i++)
          {
            mesh_velocity[i] += this->raw_dnodal_position_dt(l, i) * psif(l);
          }
        }
      }

      // Get velocities and their derivatives from base flow problem
      // -----------------------------------------------------------

      // Allocate storage for the velocity components of the base state
      // solution (initialise to zero)
      Vector<double> base_flow_u(3, 0.0);

      // Get the user-defined base state solution velocity components
      get_base_flow_u(time, ipt, interpolated_x, base_flow_u);

      // Allocate storage for the derivatives of the base state solution's
      // velocity components w.r.t. global coordinate (r and z)
      // N.B. the derivatives of the base flow components w.r.t. the
      // azimuthal coordinate direction (theta) are always zero since the
      // base flow is axisymmetric
      DenseMatrix<double> base_flow_dudx(3, 2, 0.0);

      // Get the derivatives of the user-defined base state solution
      // velocity components w.r.t. global coordinates
      get_base_flow_dudx(time, ipt, interpolated_x, base_flow_dudx);

      // Cache base flow velocities and their derivatives
      const double interpolated_ur = base_flow_u[0];
      const double interpolated_uz = base_flow_u[1];
      const double interpolated_utheta = base_flow_u[2];
      const double interpolated_durdr = base_flow_dudx(0, 0);
      const double interpolated_durdz = base_flow_dudx(0, 1);
      const double interpolated_duzdr = base_flow_dudx(1, 0);
      const double interpolated_duzdz = base_flow_dudx(1, 1);
      const double interpolated_duthetadr = base_flow_dudx(2, 0);
      const double interpolated_duthetadz = base_flow_dudx(2, 1);

      // Cache r-component of position
      const double r = interpolated_x[0];

      // Cache unknowns
      const double interpolated_UC = interpolated_u[0];
      const double interpolated_US = interpolated_u[1];
      const double interpolated_WC = interpolated_u[2];
      const double interpolated_WS = interpolated_u[3];
      const double interpolated_VC = interpolated_u[4];
      const double interpolated_VS = interpolated_u[5];
      const double interpolated_PC = interpolated_p[0];
      const double interpolated_PS = interpolated_p[1];

      // Cache derivatives of the unknowns
      const double interpolated_dUCdr = interpolated_dudx(0, 0);
      const double interpolated_dUCdz = interpolated_dudx(0, 1);
      const double interpolated_dUSdr = interpolated_dudx(1, 0);
      const double interpolated_dUSdz = interpolated_dudx(1, 1);
      const double interpolated_dWCdr = interpolated_dudx(2, 0);
      const double interpolated_dWCdz = interpolated_dudx(2, 1);
      const double interpolated_dWSdr = interpolated_dudx(3, 0);
      const double interpolated_dWSdz = interpolated_dudx(3, 1);
      const double interpolated_dVCdr = interpolated_dudx(4, 0);
      const double interpolated_dVCdz = interpolated_dudx(4, 1);
      const double interpolated_dVSdr = interpolated_dudx(5, 0);
      const double interpolated_dVSdz = interpolated_dudx(5, 1);

      // ==================
      // MOMENTUM EQUATIONS
      // ==================

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache test functions and their derivatives
        const double testf_ = testf(l);
        const double dtestfdr = dtestfdx(l, 0);
        const double dtestfdz = dtestfdx(l, 1);

        // ---------------------------------------------
        // FIRST (RADIAL) MOMENTUM EQUATION: COSINE PART
        // ---------------------------------------------

        // Get local equation number of first velocity value at this node
        local_eqn = nodal_local_eqn(l, u_nodal_index[0]);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Pressure gradient term
          residuals[local_eqn] += interpolated_PC * (testf_ + r * dtestfdr) * W;

          // Stress tensor terms
          residuals[local_eqn] -= visc_ratio * r * (1.0 + Gamma[0]) *
                                  interpolated_dUCdr * dtestfdr * W;

          residuals[local_eqn] -=
            visc_ratio * r *
            (interpolated_dUCdz + Gamma[0] * interpolated_dWCdr) * dtestfdz * W;

          residuals[local_eqn] +=
            visc_ratio *
            ((k * Gamma[0] * interpolated_dVSdr) -
             (k * (2.0 + Gamma[0]) * interpolated_VS / r) -
             ((1.0 + Gamma[0] + (k * k)) * interpolated_UC / r)) *
            testf_ * W;

          // Inertial terms (du/dt)
          residuals[local_eqn] -= scaled_re_st * r * dudt[0] * testf_ * W;

          // Inertial terms (convective)
          residuals[local_eqn] -= scaled_re *
                                  (r * interpolated_ur * interpolated_dUCdr +
                                   r * interpolated_durdr * interpolated_UC +
                                   k * interpolated_utheta * interpolated_US -
                                   2 * interpolated_utheta * interpolated_VC +
                                   r * interpolated_uz * interpolated_dUCdz +
                                   r * interpolated_durdz * interpolated_WC) *
                                  testf_ * W;

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            for (unsigned j = 0; j < 2; j++)
            {
              residuals[local_eqn] += scaled_re_st * r * mesh_velocity[j] *
                                      interpolated_dudx(0, j) * testf_ * W;
            }
          }

          // Calculate the Jacobian
          // ----------------------

          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Cache velocity shape functions and their derivatives
              const double psif_ = psif[l2];
              const double dpsifdr = dpsifdx(l2, 0);
              const double dpsifdz = dpsifdx(l2, 1);

              // Radial velocity component (cosine part) U_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  // Add the mass matrix
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif_ * testf_ * W;
                }

                // Add contributions to the Jacobian matrix

                // Stress tensor terms
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * (1.0 + Gamma[0]) * dpsifdr * dtestfdr * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * dpsifdz * dtestfdz * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * (1.0 + Gamma[0] + (k * k)) * psif_ * testf_ * W /
                  r;

                // Inertial terms (du/dt)
                jacobian(local_eqn, local_unknown) -=
                  scaled_re_st * r *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif_ *
                  testf_ * W;

                // Inertial terms (convective)
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * r *
                  (psif_ * interpolated_durdr + interpolated_ur * dpsifdr +
                   interpolated_uz * dpsifdz) *
                  testf_ * W;

                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned j = 0; j < 2; j++)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      scaled_re_st * r * mesh_velocity[j] * dpsifdx(l2, j) *
                      testf_ * W;
                  }
                }
              }

              // Radial velocity component (sine part) U_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                // Convective term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * k * interpolated_utheta * psif_ * testf_ * W;
              }

              // Axial velocity component (cosine part) W_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              if (local_unknown >= 0)
              {
                // Stress tensor term
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * Gamma[0] * r * dpsifdr * dtestfdz * W;

                // Convective term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * r * interpolated_durdz * psif_ * testf_ * W;
              }

              // Axial velocity component (sine part) W_k^S
              // has no contribution

              // Azimuthal velocity component (cosine part) V_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[4]);
              if (local_unknown >= 0)
              {
                // Convective term
                jacobian(local_eqn, local_unknown) +=
                  scaled_re * 2 * interpolated_utheta * psif_ * testf_ * W;
              }

              // Azimuthal velocity component (sine part) V_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[5]);
              if (local_unknown >= 0)
              {
                // Stress tensor term
                jacobian(local_eqn, local_unknown) +=
                  visc_ratio *
                  ((Gamma[0] * k * dpsifdr) -
                   (k * (2.0 + Gamma[0]) * psif_ / r)) *
                  testf_ * W;
              }
            } // End of loop over velocity shape functions

            // Now loop over pressure shape functions
            // (This is the contribution from pressure gradient)
            for (unsigned l2 = 0; l2 < n_pres; l2++)
            {
              // Cosine part P_k^C
              local_unknown = p_local_eqn(l2, 0);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  psip[l2] * (testf_ + r * dtestfdr) * W;
              }

              // Sine part P_k^S has no contribution

            } // End of loop over pressure shape functions
          } // End of Jacobian calculation

        } // End of if not boundary condition statement

        // --------------------------------------------
        // SECOND (RADIAL) MOMENTUM EQUATION: SINE PART
        // --------------------------------------------

        // Get local equation number of second velocity value at this node
        local_eqn = nodal_local_eqn(l, u_nodal_index[1]);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Pressure gradient term
          residuals[local_eqn] += interpolated_PS * (testf_ + r * dtestfdr) * W;

          // Stress tensor terms
          residuals[local_eqn] -= visc_ratio * r * (1.0 + Gamma[0]) *
                                  interpolated_dUSdr * dtestfdr * W;

          residuals[local_eqn] -=
            visc_ratio * r *
            (interpolated_dUSdz + Gamma[0] * interpolated_dWSdr) * dtestfdz * W;

          residuals[local_eqn] -=
            visc_ratio *
            ((k * Gamma[0] * interpolated_dVCdr) -
             (k * (2.0 + Gamma[0]) * interpolated_VC / r) +
             ((1.0 + Gamma[0] + (k * k)) * interpolated_US / r)) *
            testf_ * W;

          // Inertial terms (du/dt)
          residuals[local_eqn] -= scaled_re_st * r * dudt[1] * testf_ * W;

          // Inertial terms (convective)
          residuals[local_eqn] -= scaled_re *
                                  (r * interpolated_ur * interpolated_dUSdr +
                                   r * interpolated_durdr * interpolated_US -
                                   k * interpolated_utheta * interpolated_UC -
                                   2 * interpolated_utheta * interpolated_VS +
                                   r * interpolated_uz * interpolated_dUSdz +
                                   r * interpolated_durdz * interpolated_WS) *
                                  testf_ * W;

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            for (unsigned j = 0; j < 2; j++)
            {
              residuals[local_eqn] += scaled_re_st * r * mesh_velocity[j] *
                                      interpolated_dudx(1, j) * testf_ * W;
            }
          }

          // Calculate the Jacobian
          // ----------------------

          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Cache velocity shape functions and their derivatives
              const double psif_ = psif[l2];
              const double dpsifdr = dpsifdx(l2, 0);
              const double dpsifdz = dpsifdx(l2, 1);

              // Radial velocity component (cosine part) U_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              if (local_unknown >= 0)
              {
                // Convective term
                jacobian(local_eqn, local_unknown) +=
                  scaled_re * k * interpolated_utheta * psif_ * testf_ * W;
              }

              // Radial velocity component (sine part) U_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  // Add the mass matrix
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif_ * testf_ * W;
                }

                // Add contributions to the Jacobian matrix

                // Stress tensor terms
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * (1.0 + Gamma[0]) * dpsifdr * dtestfdr * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * dpsifdz * dtestfdz * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * (1.0 + Gamma[0] + (k * k)) * psif_ * testf_ * W /
                  r;

                // Inertial terms (du/dt)
                jacobian(local_eqn, local_unknown) -=
                  scaled_re_st * r *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif_ *
                  testf_ * W;

                // Inertial terms (convective)
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * r *
                  (psif_ * interpolated_durdr + interpolated_ur * dpsifdr +
                   interpolated_uz * dpsifdz) *
                  testf_ * W;

                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned j = 0; j < 2; j++)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      scaled_re_st * r * mesh_velocity[j] * dpsifdx(l2, j) *
                      testf_ * W;
                  }
                }
              }

              // Axial velocity component (cosine part) W_k^C
              // has no contribution

              // Axial velocity component (sine part) W_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[3]);
              if (local_unknown >= 0)
              {
                // Stress tensor term
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * Gamma[0] * r * dpsifdr * dtestfdz * W;

                // Convective term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * r * interpolated_durdz * psif_ * testf_ * W;
              }

              // Azimuthal velocity component (cosine part) V_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[4]);
              if (local_unknown >= 0)
              {
                // Stress tensor terms
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio *
                  ((Gamma[0] * k * dpsifdr) -
                   (k * (2.0 + Gamma[0]) * psif_ / r)) *
                  testf_ * W;
              }

              // Azimuthal velocity component (sine part) V_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[5]);
              if (local_unknown >= 0)
              {
                // Convective term
                jacobian(local_eqn, local_unknown) +=
                  scaled_re * 2 * interpolated_utheta * psif_ * testf_ * W;
              }
            } // End of loop over velocity shape functions

            // Now loop over pressure shape functions
            // (This is the contribution from pressure gradient)
            for (unsigned l2 = 0; l2 < n_pres; l2++)
            {
              // Cosine part P_k^C has no contribution

              // Sine part P_k^S
              local_unknown = p_local_eqn(l2, 1);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  psip[l2] * (testf_ + r * dtestfdr) * W;
              }
            } // End of loop over pressure shape functions
          } // End of Jacobian calculation

        } // End of if not boundary condition statement

        // --------------------------------------------
        // THIRD (AXIAL) MOMENTUM EQUATION: COSINE PART
        // --------------------------------------------

        // Get local equation number of third velocity value at this node
        local_eqn = nodal_local_eqn(l, u_nodal_index[2]);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Pressure gradient term
          residuals[local_eqn] += r * interpolated_PC * dtestfdz * W;

          // Stress tensor terms
          residuals[local_eqn] -=
            visc_ratio * r *
            (interpolated_dWCdr + Gamma[1] * interpolated_dUCdz) * dtestfdr * W;

          residuals[local_eqn] -= visc_ratio * r * (1.0 + Gamma[1]) *
                                  interpolated_dWCdz * dtestfdz * W;

          residuals[local_eqn] +=
            visc_ratio * k *
            (Gamma[1] * interpolated_dVSdz - k * interpolated_WC / r) * testf_ *
            W;

          // Inertial terms (du/dt)
          residuals[local_eqn] -= scaled_re_st * r * dudt[2] * testf_ * W;

          // Inertial terms (convective)
          residuals[local_eqn] -= scaled_re *
                                  (r * interpolated_ur * interpolated_dWCdr +
                                   r * interpolated_duzdr * interpolated_UC +
                                   k * interpolated_utheta * interpolated_WS +
                                   r * interpolated_uz * interpolated_dWCdz +
                                   r * interpolated_duzdz * interpolated_WC) *
                                  testf_ * W;

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            for (unsigned j = 0; j < 2; j++)
            {
              residuals[local_eqn] += scaled_re_st * r * mesh_velocity[j] *
                                      interpolated_dudx(2, j) * testf_ * W;
            }
          }

          // Calculate the Jacobian
          // ----------------------

          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Cache velocity shape functions and their derivatives
              const double psif_ = psif[l2];
              const double dpsifdr = dpsifdx(l2, 0);
              const double dpsifdz = dpsifdx(l2, 1);

              // Radial velocity component (cosine part) U_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              if (local_unknown >= 0)
              {
                // Stress tensor term
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * Gamma[1] * dpsifdz * dtestfdr * W;

                // Convective term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * r * psif_ * interpolated_duzdr * testf_ * W;
              }

              // Radial velocity component (sine part) U_k^S
              // has no contribution

              // Axial velocity component (cosine part) W_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  // Add the mass matrix
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif_ * testf_ * W;
                }

                // Add contributions to the Jacobian matrix

                // Stress tensor terms
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * dpsifdr * dtestfdr * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * (1.0 + Gamma[1]) * dpsifdz * dtestfdz * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * k * k * psif_ * testf_ * W / r;

                // Inertial terms (du/dt)
                jacobian(local_eqn, local_unknown) -=
                  scaled_re_st * r *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif_ *
                  testf_ * W;

                // Inertial terms (convective)
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * r *
                  (interpolated_ur * dpsifdr + psif_ * interpolated_duzdz +
                   interpolated_uz * dpsifdz) *
                  testf_ * W;


                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned j = 0; j < 2; j++)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      scaled_re_st * r * mesh_velocity[j] * dpsifdx(l2, j) *
                      testf_ * W;
                  }
                }
              }

              // Axial velocity component (sine part) W_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[3]);
              if (local_unknown >= 0)
              {
                // Convective term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * k * interpolated_utheta * psif_ * testf_ * W;
              }

              // Azimuthal velocity component (cosine part) V_k^C
              // has no contribution

              // Azimuthal velocity component (sine part) V_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[5]);
              if (local_unknown >= 0)
              {
                // Stress tensor term
                jacobian(local_eqn, local_unknown) +=
                  visc_ratio * Gamma[1] * k * dpsifdz * testf_ * W;
              }
            } // End of loop over velocity shape functions

            // Now loop over pressure shape functions
            // (This is the contribution from pressure gradient)
            for (unsigned l2 = 0; l2 < n_pres; l2++)
            {
              // Cosine part P_k^C
              local_unknown = p_local_eqn(l2, 0);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  r * psip[l2] * dtestfdz * W;
              }

              // Sine part P_k^S has no contribution

            } // End of loop over pressure shape functions
          } // End of Jacobian calculation

        } // End of if not boundary condition statement

        // -------------------------------------------
        // FOURTH (AXIAL) MOMENTUM EQUATION: SINE PART
        // -------------------------------------------

        // Get local equation number of fourth velocity value at this node
        local_eqn = nodal_local_eqn(l, u_nodal_index[3]);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Pressure gradient term
          residuals[local_eqn] += r * interpolated_PS * dtestfdz * W;

          // Stress tensor terms
          residuals[local_eqn] -=
            visc_ratio * r *
            (interpolated_dWSdr + Gamma[1] * interpolated_dUSdz) * dtestfdr * W;

          residuals[local_eqn] -= visc_ratio * r * (1.0 + Gamma[1]) *
                                  interpolated_dWSdz * dtestfdz * W;

          residuals[local_eqn] -=
            visc_ratio * k *
            (Gamma[1] * interpolated_dVCdz + k * interpolated_WS / r) * testf_ *
            W;

          // Inertial terms (du/dt)
          residuals[local_eqn] -= scaled_re_st * r * dudt[3] * testf_ * W;

          // Inertial terms (convective)
          residuals[local_eqn] -= scaled_re *
                                  (r * interpolated_ur * interpolated_dWSdr +
                                   r * interpolated_duzdr * interpolated_US -
                                   k * interpolated_utheta * interpolated_WC +
                                   r * interpolated_uz * interpolated_dWSdz +
                                   r * interpolated_duzdz * interpolated_WS) *
                                  testf_ * W;

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            for (unsigned j = 0; j < 2; j++)
            {
              residuals[local_eqn] += scaled_re_st * r * mesh_velocity[j] *
                                      interpolated_dudx(3, j) * testf_ * W;
            }
          }

          // Calculate the Jacobian
          // ----------------------

          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Cache velocity shape functions and their derivatives
              const double psif_ = psif[l2];
              const double dpsifdr = dpsifdx(l2, 0);
              const double dpsifdz = dpsifdx(l2, 1);

              // Radial velocity component (cosine part) U_k^C
              // has no contribution

              // Radial velocity component (sine part) U_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                // Stress tensor term
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * Gamma[1] * dpsifdz * dtestfdr * W;

                // Convective term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * r * psif_ * interpolated_duzdr * testf_ * W;
              }

              // Axial velocity component (cosine part) W_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              if (local_unknown >= 0)
              {
                // Convective term
                jacobian(local_eqn, local_unknown) +=
                  scaled_re * k * interpolated_utheta * psif_ * testf_ * W;
              }

              // Axial velocity component (sine part) W_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[3]);
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  // Add the mass matrix
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif_ * testf_ * W;
                }

                // Add contributions to the Jacobian matrix

                // Stress tensor terms
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * dpsifdr * dtestfdr * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * (1.0 + Gamma[1]) * dpsifdz * dtestfdz * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * k * k * psif_ * testf_ * W / r;

                // Inertial terms (du/dt)
                jacobian(local_eqn, local_unknown) -=
                  scaled_re_st * r *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif_ *
                  testf_ * W;

                // Inertial terms (convective)
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * r *
                  (interpolated_ur * dpsifdr + psif_ * interpolated_duzdz +
                   interpolated_uz * dpsifdz) *
                  testf_ * W;


                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned j = 0; j < 2; j++)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      scaled_re_st * r * mesh_velocity[j] * dpsifdx(l2, j) *
                      testf_ * W;
                  }
                }
              }

              // Azimuthal velocity component (cosine part) V_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[4]);
              if (local_unknown >= 0)
              {
                // Stress tensor term
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * Gamma[1] * k * dpsifdz * testf_ * W;
              }

              // Azimuthal velocity component (sine part) V_k^S
              // has no contribution

            } // End of loop over velocity shape functions

            // Now loop over pressure shape functions
            // (This is the contribution from pressure gradient)
            for (unsigned l2 = 0; l2 < n_pres; l2++)
            {
              // Cosine part P_k^C has no contribution

              // Sine part P_k^S
              local_unknown = p_local_eqn(l2, 1);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  r * psip[l2] * dtestfdz * W;
              }
            } // End of loop over pressure shape functions
          } // End of Jacobian calculation

        } // End of if not boundary condition statement

        // ------------------------------------------------
        // FIFTH (AZIMUTHAL) MOMENTUM EQUATION: COSINE PART
        // ------------------------------------------------

        // Get local equation number of fifth velocity value at this node
        local_eqn = nodal_local_eqn(l, u_nodal_index[4]);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Pressure gradient term
          residuals[local_eqn] -= k * interpolated_PS * testf_ * W;

          // Stress tensor terms
          residuals[local_eqn] +=
            visc_ratio *
            (-r * interpolated_dVCdr - k * Gamma[0] * interpolated_US +
             Gamma[0] * interpolated_VC) *
            dtestfdr * W;

          residuals[local_eqn] -=
            visc_ratio *
            (k * Gamma[0] * interpolated_WS + r * interpolated_dVCdz) *
            dtestfdz * W;

          residuals[local_eqn] +=
            visc_ratio *
            (Gamma[0] * interpolated_dVCdr +
             k * (2.0 + Gamma[0]) * interpolated_US / r -
             (1.0 + (k * k) + (k * k * Gamma[0])) * interpolated_VC / r) *
            testf_ * W;

          // Inertial terms (du/dt)
          residuals[local_eqn] -= scaled_re_st * r * dudt[4] * testf_ * W;

          // Inertial terms (convective)
          residuals[local_eqn] -=
            scaled_re *
            (r * interpolated_ur * interpolated_dVCdr +
             r * interpolated_duthetadr * interpolated_UC +
             k * interpolated_utheta * interpolated_VS +
             interpolated_utheta * interpolated_UC +
             interpolated_ur * interpolated_VC +
             r * interpolated_uz * interpolated_dVCdz +
             r * interpolated_duthetadz * interpolated_WC) *
            testf_ * W;

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            for (unsigned j = 0; j < 2; j++)
            {
              residuals[local_eqn] += scaled_re_st * r * mesh_velocity[j] *
                                      interpolated_dudx(4, j) * testf_ * W;
            }
          }

          // Calculate the Jacobian
          // ----------------------

          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Cache velocity shape functions and their derivatives
              const double psif_ = psif[l2];
              const double dpsifdr = dpsifdx(l2, 0);
              const double dpsifdz = dpsifdx(l2, 1);

              // Radial velocity component (cosine part) U_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              if (local_unknown >= 0)
              {
                // Convective terms
                jacobian(local_eqn, local_unknown) -=
                  scaled_re *
                  (r * interpolated_duthetadr + interpolated_utheta) * psif_ *
                  testf_ * W;
              }

              // Radial velocity component (sine part) U_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                // Stress tensor terms
                jacobian(local_eqn, local_unknown) +=
                  visc_ratio * k * psif_ *
                  (((2.0 + Gamma[0]) * testf_ / r) - (Gamma[0] * dtestfdr)) * W;
              }

              // Axial velocity component (cosine part) W_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              if (local_unknown >= 0)
              {
                // Convective term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * r * psif_ * interpolated_duthetadz * testf_ * W;
              }

              // Axial velocity component (sine part) W_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[3]);
              if (local_unknown >= 0)
              {
                // Stress tensor term
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * k * Gamma[0] * psif_ * dtestfdz * W;
              }

              // Azimuthal velocity component (cosine part) V_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[4]);
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  // Add the mass matrix
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif_ * testf_ * W;
                }

                // Add contributions to the Jacobian matrix

                // Stress tensor terms
                jacobian(local_eqn, local_unknown) +=
                  visc_ratio * (Gamma[0] * psif_ - r * dpsifdr) * dtestfdr * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * dpsifdz * dtestfdz * W;

                jacobian(local_eqn, local_unknown) +=
                  visc_ratio *
                  (Gamma[0] * dpsifdr -
                   (1.0 + (k * k) + (k * k * Gamma[0])) * psif_ / r) *
                  testf_ * W;

                // Inertial terms (du/dt)
                jacobian(local_eqn, local_unknown) -=
                  scaled_re_st * r *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif_ *
                  testf_ * W;

                // Inertial terms (convective)
                jacobian(local_eqn, local_unknown) -=
                  scaled_re *
                  (r * interpolated_ur * dpsifdr + interpolated_ur * psif_ +
                   r * interpolated_uz * dpsifdz) *
                  testf_ * W;

                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned j = 0; j < 2; j++)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      scaled_re_st * r * mesh_velocity[j] * dpsifdx(l2, j) *
                      testf_ * W;
                  }
                }
              }

              // Azimuthal velocity component (sine part) V_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[5]);
              if (local_unknown >= 0)
              {
                // Convective term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * k * interpolated_utheta * psif_ * testf_ * W;
              }
            } // End of loop over velocity shape functions

            // Now loop over pressure shape functions
            // (This is the contribution from pressure gradient)
            for (unsigned l2 = 0; l2 < n_pres; l2++)
            {
              // Cosine part P_k^C has no contribution

              // Sine part P_k^S
              local_unknown = p_local_eqn(l2, 1);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) -= k * psip[l2] * testf_ * W;
              }
            } // End of loop over pressure shape functions
          } // End of Jacobian calculation

        } // End of if not boundary condition statement

        // ----------------------------------------------
        // SIXTH (AZIMUTHAL) MOMENTUM EQUATION: SINE PART
        // ----------------------------------------------

        // Get local equation number of sixth velocity value at this node
        local_eqn = nodal_local_eqn(l, u_nodal_index[5]);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Pressure gradient term
          residuals[local_eqn] += k * interpolated_PC * testf_ * W;

          // Stress tensor terms
          residuals[local_eqn] +=
            visc_ratio *
            (-r * interpolated_dVSdr + k * Gamma[0] * interpolated_UC +
             Gamma[0] * interpolated_VS) *
            dtestfdr * W;

          residuals[local_eqn] +=
            visc_ratio *
            (k * Gamma[0] * interpolated_WC - r * interpolated_dVSdz) *
            dtestfdz * W;

          residuals[local_eqn] +=
            visc_ratio *
            (Gamma[0] * interpolated_dVSdr -
             k * (2.0 + Gamma[0]) * interpolated_UC / r -
             (1.0 + (k * k) + (k * k * Gamma[0])) * interpolated_VS / r) *
            testf_ * W;

          // Inertial terms (du/dt)
          residuals[local_eqn] -= scaled_re_st * r * dudt[5] * testf_ * W;

          // Inertial terms (convective)
          residuals[local_eqn] -=
            scaled_re *
            (r * interpolated_ur * interpolated_dVSdr +
             r * interpolated_duthetadr * interpolated_US -
             k * interpolated_utheta * interpolated_VC +
             interpolated_utheta * interpolated_US +
             interpolated_ur * interpolated_VS +
             r * interpolated_uz * interpolated_dVSdz +
             r * interpolated_duthetadz * interpolated_WS) *
            testf_ * W;

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            for (unsigned j = 0; j < 2; j++)
            {
              residuals[local_eqn] += scaled_re_st * r * mesh_velocity[j] *
                                      interpolated_dudx(5, j) * testf_ * W;
            }
          }

          // Calculate the Jacobian
          // ----------------------

          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Cache velocity shape functions and their derivatives
              const double psif_ = psif[l2];
              const double dpsifdr = dpsifdx(l2, 0);
              const double dpsifdz = dpsifdx(l2, 1);

              // Radial velocity component (cosine part) U_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              if (local_unknown >= 0)
              {
                // Stress tensor terms
                jacobian(local_eqn, local_unknown) +=
                  visc_ratio * k * psif_ *
                  ((Gamma[0] * dtestfdr) - ((2.0 + Gamma[0]) * testf_ / r)) * W;
              }

              // Radial velocity component (sine part) U_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                // Convective terms
                jacobian(local_eqn, local_unknown) -=
                  scaled_re *
                  (r * interpolated_duthetadr + interpolated_utheta) * psif_ *
                  testf_ * W;
              }

              // Axial velocity component (cosine part) W_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              if (local_unknown >= 0)
              {
                // Stress tensor term
                jacobian(local_eqn, local_unknown) +=
                  visc_ratio * k * Gamma[0] * psif_ * dtestfdz * W;
              }

              // Axial velocity component (sine part) W_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[3]);
              if (local_unknown >= 0)
              {
                // Convective term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * r * psif_ * interpolated_duthetadz * testf_ * W;
              }

              // Azimuthal velocity component (cosine part) V_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[4]);
              if (local_unknown >= 0)
              {
                // Convective term
                jacobian(local_eqn, local_unknown) +=
                  scaled_re * k * interpolated_utheta * psif_ * testf_ * W;
              }

              // Azimuthal velocity component (sine part) V_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[5]);
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  // Add the mass matrix
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif_ * testf_ * W;
                }

                // Add contributions to the Jacobian matrix

                // Stress tensor terms
                jacobian(local_eqn, local_unknown) +=
                  visc_ratio * (Gamma[0] * psif_ - r * dpsifdr) * dtestfdr * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * dpsifdz * dtestfdz * W;

                jacobian(local_eqn, local_unknown) +=
                  visc_ratio *
                  (Gamma[0] * dpsifdr -
                   (1.0 + (k * k) + (k * k * Gamma[0])) * psif_ / r) *
                  testf_ * W;

                // Inertial terms (du/dt)
                jacobian(local_eqn, local_unknown) -=
                  scaled_re_st * r *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif_ *
                  testf_ * W;

                // Convective terms
                jacobian(local_eqn, local_unknown) -=
                  scaled_re *
                  (r * interpolated_ur * dpsifdr + interpolated_ur * psif_ +
                   r * interpolated_uz * dpsifdz) *
                  testf_ * W;

                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned j = 0; j < 2; j++)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      scaled_re_st * r * mesh_velocity[j] * dpsifdx(l2, j) *
                      testf_ * W;
                  }
                }
              }
            } // End of loop over velocity shape functions

            // Now loop over pressure shape functions
            // (This is the contribution from pressure gradient)
            for (unsigned l2 = 0; l2 < n_pres; l2++)
            {
              // Cosine part P_k^C
              local_unknown = p_local_eqn(l2, 0);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) += k * psip[l2] * testf_ * W;
              }

              // Sine part P_k^S has no contribution

            } // End of loop over pressure shape functions
          } // End of Jacobian calculation

        } // End of if not boundary condition statement

      } // End of loop over shape functions


      // ====================
      // CONTINUITY EQUATIONS
      // ====================

      // Loop over the pressure shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Cache test function
        const double testp_ = testp[l];

        // --------------------------------------
        // FIRST CONTINUITY EQUATION: COSINE PART
        // --------------------------------------

        // Get local equation number of first pressure value at this node
        local_eqn = p_local_eqn(l, 0);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Gradient terms
          residuals[local_eqn] +=
            (interpolated_UC + r * interpolated_dUCdr + k * interpolated_VS +
             r * interpolated_dWCdz) *
            testp_ * W;

          // Calculate the Jacobian
          // ----------------------

          if (flag)
          {
            // Loop over the velocity shape functions
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Cache velocity shape functions and their derivatives
              const double psif_ = psif[l2];
              const double dpsifdr = dpsifdx(l2, 0);
              const double dpsifdz = dpsifdx(l2, 1);

              // Radial velocity component (cosine part) U_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  (psif_ + r * dpsifdr) * testp_ * W;
              }

              // Radial velocity component (sine part) U_k^S
              // has no contribution

              // Axial velocity component (cosine part) W_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) += r * dpsifdz * testp_ * W;
              }

              // Axial velocity component (sine part) W_k^S
              // has no contribution

              // Azimuthal velocity component (cosine part) V_k^C
              // has no contribution

              // Azimuthal velocity component (sine part) V_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[5]);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) += k * psif_ * testp_ * W;
              }
            } // End of loop over velocity shape functions

            // Real and imaginary pressure components, P_k^C and P_k^S,
            // have no contribution

          } // End of Jacobian calculation

        } // End of if not boundary condition statement

        // -------------------------------------
        // SECOND CONTINUITY EQUATION: SINE PART
        // -------------------------------------

        // Get local equation number of second pressure value at this node
        local_eqn = p_local_eqn(l, 1);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Gradient terms
          residuals[local_eqn] +=
            (interpolated_US + r * interpolated_dUSdr - k * interpolated_VC +
             r * interpolated_dWSdz) *
            testp_ * W;

          // Calculate the Jacobian
          // ----------------------

          if (flag)
          {
            // Loop over the velocity shape functions
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Cache velocity shape functions and their derivatives
              const double psif_ = psif[l2];
              const double dpsifdr = dpsifdx(l2, 0);
              const double dpsifdz = dpsifdx(l2, 1);

              // Radial velocity component (cosine part) U_k^C
              // has no contribution

              // Radial velocity component (sine part) U_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  (psif_ + r * dpsifdr) * testp_ * W;
              }

              // Axial velocity component (cosine part) W_k^C
              // has no contribution

              // Axial velocity component (sine part) W_k^S
              local_unknown = nodal_local_eqn(l2, u_nodal_index[3]);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) += r * dpsifdz * testp_ * W;
              }

              // Azimuthal velocity component (cosine part) V_k^C
              local_unknown = nodal_local_eqn(l2, u_nodal_index[4]);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) -= k * psif_ * testp_ * W;
              }

              // Azimuthal velocity component (sine part) V_k^S
              // has no contribution

            } // End of loop over velocity shape functions

            // Real and imaginary pressure components, P_k^C and P_k^S,
            // have no contribution

          } // End of Jacobian calculation

        } // End of if not boundary condition statement

      } // End of loop over pressure shape functions

    } // End of loop over the integration points

  } // End of fill_in_generic_residual_contribution_linearised_axi_nst


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


  /// Linearised axisymmetric Crouzeix-Raviart elements
  /// -------------------------------------------------


  //=======================================================================
  /// Set the data for the number of variables at each node
  //=======================================================================
  const unsigned
    LinearisedAxisymmetricQCrouzeixRaviartElement::Initial_Nvalue[9] = {
      6, 6, 6, 6, 6, 6, 6, 6, 6};


  //========================================================================
  /// Number of values (pinned or dofs) required at node n
  //========================================================================
  unsigned LinearisedAxisymmetricQCrouzeixRaviartElement::required_nvalue(
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
  const unsigned LinearisedAxisymmetricQTaylorHoodElement::Initial_Nvalue[9] = {
    8, 6, 8, 6, 6, 6, 8, 6, 8};


  //=======================================================================
  /// Set the data for the pressure conversion array
  //=======================================================================
  const unsigned LinearisedAxisymmetricQTaylorHoodElement::Pconv[4] = {
    0, 2, 6, 8};


} // namespace oomph
