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
// Non-inline functions for elements that solve the equations of linear
// elasticity in cartesian coordinates

#include "time_harmonic_fourier_decomposed_linear_elasticity_elements.h"


namespace oomph
{
  /// Static default value for square of angular frequency: Zero
  std::complex<double> TimeHarmonicFourierDecomposedLinearElasticityEquationsBase::
    Default_omega_sq_value(0.0, 0.0);

  /// Static default value for Young's modulus (1.0 -- for natural
  /// scaling, i.e. all stresses have been non-dimensionalised by
  /// the same reference Young's modulus. Setting the "non-dimensional"
  /// Young's modulus (obtained by de-referencing Youngs_modulus_pt)
  /// to a number larger than one means that the material is stiffer
  /// than assumed in the non-dimensionalisation.
  std::complex<double> TimeHarmonicFourierDecomposedLinearElasticityEquationsBase::
    Default_youngs_modulus_value(1.0, 0.0);


  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////


  //=======================================================================
  /// Compute norm of the solution
  //=======================================================================
  void TimeHarmonicFourierDecomposedLinearElasticityEquations::compute_norm(
    double& norm)
  {
    // Initialise
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordintes
    Vector<double> x(2);

    // Displacement vector
    Vector<std::complex<double>> disp(3);

    // Find out how many nodes there are in the element
    unsigned n_node = this->nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = this->J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get FE function value
      this->interpolated_u_time_harmonic_fourier_decomposed_linear_elasticity(
        s, disp);

      // Add to  norm
      for (unsigned ii = 0; ii < 3; ii++)
      {
        norm += (disp[ii].real() * disp[ii].real() +
                 disp[ii].imag() * disp[ii].imag()) *
                W;
      }
    }
  }

  //=======================================================================
  /// Get strain tensor
  //=======================================================================
  void TimeHarmonicFourierDecomposedLinearElasticityEquations::get_strain(
    const Vector<double>& s, DenseMatrix<std::complex<double>>& strain)
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

#ifdef PARANOID
    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();

    if (n_position_type != 1)
    {
      throw OomphLibError("TimeHarmonicFourierDecomposedLinearElasticity is "
                          "not yet implemented for more than one position type",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the indices at which the local displacements are stored
    std::complex<unsigned> u_nodal_index[3];
    for (unsigned i = 0; i < 3; i++)
    {
      u_nodal_index[i] =
        this->u_index_time_harmonic_fourier_decomposed_linear_elasticity(i);
    }

    // Fourier wavenumber as double
    double n = double(this->fourier_wavenumber());

    // Set up memory for the shape functions
    Shape psi(n_node);

    // Derivs only w.r.t. r [0] and z [1]
    DShape dpsidx(n_node, 2);

    // Call the derivatives of the shape functions (ignore Jacobian)
    this->dshape_eulerian(s, psi, dpsidx);

    // Storage for Eulerian coordinates (r,z; initialised to zero)
    Vector<double> interpolated_x(2, 0.0);

    // Displacements u_r,u_z,u_theta
    Vector<std::complex<double>> interpolated_u(3,
                                                std::complex<double>(0.0, 0.0));

    // Calculate interpolated values of the derivatives w.r.t.
    // Eulerian coordinates
    DenseMatrix<std::complex<double>> interpolated_dudx(
      3, 2, std::complex<double>(0.0, 0.0));

    // Calculate displacements and derivatives
    for (unsigned l = 0; l < n_node; l++)
    {
      // Calculate the coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        interpolated_x[i] += this->nodal_position(l, i) * psi(l);
      }
      // Get the nodal displacements
      for (unsigned i = 0; i < 3; i++)
      {
        const std::complex<double> u_value =
          std::complex<double>(this->nodal_value(l, u_nodal_index[i].real()),
                               this->nodal_value(l, u_nodal_index[i].imag()));

        interpolated_u[i] += u_value * psi(l);

        // Loop over derivative directions
        for (unsigned j = 0; j < 2; j++)
        {
          interpolated_dudx(i, j) += u_value * dpsidx(l, j);
        }
      }
    }

    // define shorthand notation for regularly-occurring terms
    double r = interpolated_x[0];
    const std::complex<double> I(0.0, 1.0);

    // r component of displacement
    std::complex<double> ur = interpolated_u[0];

    // z component of displacement
    std::complex<double> uz = interpolated_u[1];

    // theta component of displacement
    std::complex<double> uth = interpolated_u[2];

    // derivatives w.r.t. r and z:
    std::complex<double> durdr = interpolated_dudx(0, 0);
    std::complex<double> durdz = interpolated_dudx(0, 1);
    std::complex<double> duzdr = interpolated_dudx(1, 0);
    std::complex<double> duzdz = interpolated_dudx(1, 1);
    std::complex<double> duthdr = interpolated_dudx(2, 0);
    std::complex<double> duthdz = interpolated_dudx(2, 1);

    strain(0, 0) = durdr;
    strain(2, 2) = I * double(n) * uth / r + ur / r;
    strain(1, 1) = duzdz;
    strain(0, 2) = 0.5 * (I * double(n) * ur / r - uth / r + duthdr);
    strain(2, 0) = 0.5 * (I * double(n) * ur / r - uth / r + duthdr);
    strain(0, 1) = 0.5 * (durdz + duzdr);
    strain(1, 0) = 0.5 * (durdz + duzdr);
    strain(2, 1) = 0.5 * (duthdz + I * double(n) * uz / r);
    strain(1, 2) = 0.5 * (duthdz + I * double(n) * uz / r);
  }

  /// ////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////

  //=======================================================================
  /// Compute the residuals for the Fourier decomposed (in cyl. polars)
  /// time harmonic linear elasticity equations in. Flag indicates if
  /// we want Jacobian too.
  //=======================================================================
  void TimeHarmonicFourierDecomposedLinearElasticityEquations::
    fill_in_generic_contribution_to_residuals_fourier_decomp_time_harmonic_linear_elasticity(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

#ifdef PARANOID
    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();

    if (n_position_type != 1)
    {
      throw OomphLibError("TimeHarmonicFourierDecomposedLinearElasticity is "
                          "not yet implemented for more than one position type",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the indices at which the local displacements are stored
    std::complex<unsigned> u_nodal_index[3];
    for (unsigned i = 0; i < 3; i++)
    {
      u_nodal_index[i] =
        this->u_index_time_harmonic_fourier_decomposed_linear_elasticity(i);
    }

    // Get (complex) elastic parameters
    std::complex<double> nu_local = this->nu();
    std::complex<double> youngs_modulus_local = this->youngs_modulus();

    // Obtain Lame parameters from Young's modulus and Poisson's ratio
    std::complex<double> lambda = youngs_modulus_local * nu_local /
                                  (1.0 + nu_local) / (1.0 - 2.0 * nu_local);
    std::complex<double> mu = youngs_modulus_local / 2.0 / (1.0 + nu_local);


    // Fourier wavenumber as double
    double n = double(this->fourier_wavenumber());

    // Square of non-dimensional frequency
    const std::complex<double> omega_sq = this->omega_sq();

    // Set up memory for the shape functions
    Shape psi(n_node);

    // Derivs only w.r.t. r [0] and z [1]
    DShape dpsidx(n_node, 2);

    // Set the value of Nintpt -- the number of integration points
    unsigned n_intpt = this->integral_pt()->nweight();

    // Set the vector to hold the local coordinates in the element
    Vector<double> s(2);

    // Integers to store the local equation numbers
    int local_eqn_real = 0, local_eqn_imag = 0, local_unknown_real = 0,
        local_unknown_imag = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign the values of s
      for (unsigned i = 0; i < 2; ++i)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Call the derivatives of the shape functions (and get Jacobian)
      double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

      // Storage for Eulerian coordinates (r,z; initialised to zero)
      Vector<double> interpolated_x(2, 0.0);

      // Displacements u_r,u_z,u_theta
      Vector<std::complex<double>> interpolated_u(
        3, std::complex<double>(0.0, 0.0));

      // Calculate interpolated values of the derivatives w.r.t.
      // Eulerian coordinates
      DenseMatrix<std::complex<double>> interpolated_dudx(
        3, 2, std::complex<double>(0.0, 0.0));

      // Calculate displacements and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Calculate the coordinates
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += this->raw_nodal_position(l, i) * psi(l);
        }
        // Get the nodal displacements
        for (unsigned i = 0; i < 3; i++)
        {
          const std::complex<double> u_value = std::complex<double>(
            this->raw_nodal_value(l, u_nodal_index[i].real()),
            this->raw_nodal_value(l, u_nodal_index[i].imag()));

          interpolated_u[i] += u_value * psi(l);

          // Loop over derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsidx(l, j);
          }
        }
      }

      // Get body force
      Vector<std::complex<double>> b(3);
      this->body_force(interpolated_x, b);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      //=====EQUATIONS OF FOURIER DECOMPOSED TIME HARMONIC LINEAR ELASTICITY
      //========

      // define shorthand notation for regularly-occurring terms
      double r = interpolated_x[0];
      double rsq = pow(r, 2);
      double nsq = pow(n, 2);
      const std::complex<double> I(0.0, 1.0);

      // r component of displacement
      std::complex<double> ur = interpolated_u[0];

      // z component of displacement
      std::complex<double> uz = interpolated_u[1];

      // theta component of displacement
      std::complex<double> uth = interpolated_u[2];

      // derivatives w.r.t. r and z:
      std::complex<double> durdr = interpolated_dudx(0, 0);
      std::complex<double> durdz = interpolated_dudx(0, 1);
      std::complex<double> duzdr = interpolated_dudx(1, 0);
      std::complex<double> duzdz = interpolated_dudx(1, 1);
      std::complex<double> duthdr = interpolated_dudx(2, 0);
      std::complex<double> duthdz = interpolated_dudx(2, 1);

      // storage for (complex) terms required for analytic Jacobian
      std::complex<double> G_r, G_z, G_theta;

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the displacement components
        for (unsigned a = 0; a < 3; a++)
        {
          // Get the REAL equation number
          local_eqn_real = this->nodal_local_eqn(l, u_nodal_index[a].real());

          /*IF it's not a boundary condition*/
          if (local_eqn_real >= 0)
          {
            // Acceleration and body force
            residuals[local_eqn_real] +=
              (-omega_sq.real() * interpolated_u[a].real() +
               omega_sq.imag() * interpolated_u[a].imag() - b[a].real()) *
              psi(l) * r * W;

            // Three components of the stress divergence term:
            // a=0: r; a=1: z; a=2: theta

            // Real part of r-equation
            if (a == 0)
            {
              residuals[local_eqn_real] +=
                (mu.real() *
                   (2.0 * durdr.real() * dpsidx(l, 0) +
                    n * psi(l) / r *
                      (n * ur.real() / r + duthdr.imag() - uth.imag() / r) +
                    dpsidx(l, 1) * (durdz.real() + duzdr.real()) +
                    2.0 * psi(l) / pow(r, 2) * (ur.real() - n * uth.imag())) +
                 mu.imag() *
                   (-2.0 * durdr.imag() * dpsidx(l, 0) +
                    n * psi(l) / r *
                      (-n * ur.imag() / r + duthdr.real() - uth.real() / r) -
                    dpsidx(l, 1) * (durdz.imag() + duzdr.imag()) -
                    2.0 * psi(l) / pow(r, 2) * (ur.imag() + n * uth.real())) +
                 lambda.real() *
                   (durdr.real() + ur.real() / r - n / r * uth.imag() +
                    duzdz.real()) *
                   (dpsidx(l, 0) + psi(l) / r) -
                 lambda.imag() *
                   (durdr.imag() + ur.imag() / r + n / r * uth.real() +
                    duzdz.imag()) *
                   (dpsidx(l, 0) + psi(l) / r)) *
                r * W;
            }
            // Real part of z-equation
            else if (a == 1)
            {
              residuals[local_eqn_real] +=
                (mu.real() *
                   (dpsidx(l, 0) * (durdz.real() + duzdr.real()) +
                    n * psi(l) / r * (n * uz.real() / r + duthdz.imag()) +
                    2.0 * duzdz.real() * dpsidx(l, 1)) +
                 mu.imag() *
                   (-dpsidx(l, 0) * (durdz.imag() + duzdr.imag()) +
                    n * psi(l) / r * (-n * uz.imag() / r + duthdz.real()) -
                    2.0 * duzdz.imag() * dpsidx(l, 1)) +
                 lambda.real() *
                   (durdr.real() + ur.real() / r - n / r * uth.imag() +
                    duzdz.real()) *
                   dpsidx(l, 1) -
                 lambda.imag() *
                   (durdr.imag() + ur.imag() / r + n / r * uth.real() +
                    duzdz.imag()) *
                   dpsidx(l, 1)) *
                r * W;
            }
            // Real part of theta-equation
            else if (a == 2)
            {
              residuals[local_eqn_real] +=
                (mu.real() *
                   ((duthdr.real() - uth.real() / r - n * ur.imag() / r) *
                      (dpsidx(l, 0) - psi(l) / r) +
                    2.0 * n * psi(l) / pow(r, 2) *
                      (n * uth.real() + ur.imag()) +
                    dpsidx(l, 1) * (duthdz.real() - n / r * uz.imag())) +
                 mu.imag() *
                   ((-duthdr.imag() + uth.imag() / r - n * ur.real() / r) *
                      (dpsidx(l, 0) - psi(l) / r) +
                    2.0 * n * psi(l) / pow(r, 2) *
                      (-n * uth.imag() + ur.real()) +
                    dpsidx(l, 1) * (-duthdz.imag() - n / r * uz.real())) +
                 lambda.real() *
                   (durdr.imag() + ur.imag() / r + n / r * uth.real() +
                    duzdz.imag()) *
                   n * psi(l) / r +
                 lambda.imag() *
                   (durdr.real() + ur.real() / r - n / r * uth.imag() +
                    duzdz.real()) *
                   n * psi(l) / r) *
                r * W;
            }
            // error: a should be 0, 1 or 2
            else
            {
              throw OomphLibError("a should equal 0, 1 or 2",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }

            // Jacobian entries
            if (flag)
            {
              // Loop over the displacement basis functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // define complex terms used to obtain entries of current row in
                // the Jacobian:

                // terms for rows of Jacobian matrix corresponding to r-equation
                if (a == 0)
                {
                  G_r = (mu * (2.0 * dpsidx(l2, 0) * dpsidx(l, 0) +
                               (nsq + 2.0) / rsq * psi(l2) * psi(l) +
                               dpsidx(l2, 1) * dpsidx(l, 1)) +
                         lambda * (dpsidx(l2, 0) + psi(l2) / r) *
                           (dpsidx(l, 0) + psi(l) / r) -
                         omega_sq * psi(l2) * psi(l)) *
                        r * W;

                  G_z = (mu * dpsidx(l2, 0) * dpsidx(l, 1) +
                         lambda * dpsidx(l2, 1) * (dpsidx(l, 0) + psi(l) / r)) *
                        r * W;

                  G_theta = (-I * (mu * (n / r * dpsidx(l2, 0) * psi(l) -
                                         3.0 * n / rsq * psi(l2) * psi(l)) -
                                   lambda * n / r * psi(l2) *
                                     (dpsidx(l, 0) + psi(l) / r))) *
                            r * W;
                }
                // terms for rows of Jacobian matrix corresponding to z-equation
                else if (a == 1)
                {
                  G_r =
                    (mu * dpsidx(l2, 1) * dpsidx(l, 0) +
                     lambda * (dpsidx(l2, 0) + psi(l2) / r) * dpsidx(l, 1)) *
                    r * W;

                  G_z = (mu * (dpsidx(l2, 0) * dpsidx(l, 0) +
                               nsq / rsq * psi(l2) * psi(l) +
                               2.0 * dpsidx(l2, 1) * dpsidx(l, 1)) +
                         lambda * dpsidx(l2, 1) * dpsidx(l, 1) -
                         omega_sq * psi(l2) * psi(l)) *
                        r * W;

                  G_theta = (-I * (mu * n / r * dpsidx(l2, 1) * psi(l) -
                                   lambda * n / r * dpsidx(l, 1) * psi(l2))) *
                            r * W;
                }
                // terms for rows of Jacobian matrix corresponding to
                // theta-equation
                else if (a == 2)
                {
                  G_r = (-I * (mu * (-n / r * psi(l2) * dpsidx(l, 0) +
                                     3.0 * n / rsq * psi(l2) * psi(l)) +
                               lambda * n / r * (dpsidx(l2, 0) + psi(l2) / r) *
                                 psi(l))) *
                        r * W;

                  G_z = (-I * (-mu * n / r * psi(l2) * dpsidx(l, 1) +
                               lambda * n / r * dpsidx(l2, 1) * psi(l))) *
                        r * W;

                  G_theta = (mu * ((dpsidx(l2, 0) - psi(l2) / r) *
                                     (dpsidx(l, 0) - psi(l) / r) +
                                   2.0 * nsq / rsq * psi(l2) * psi(l) +
                                   dpsidx(l2, 1) * dpsidx(l, 1)) +
                             lambda * nsq / rsq * psi(l2) * psi(l) -
                             omega_sq * psi(l2) * psi(l)) *
                            r * W;
                }

                // Loop over the displacement components
                for (unsigned c = 0; c < 3; c++)
                {
                  // Get real and imaginary parts of local unknown
                  local_unknown_real =
                    this->nodal_local_eqn(l2, u_nodal_index[c].real());
                  local_unknown_imag =
                    this->nodal_local_eqn(l2, u_nodal_index[c].imag());

                  // If the real part of the local unknown is not pinned
                  if (local_unknown_real >= 0)
                  {
                    if (c == 0)
                    {
                      jacobian(local_eqn_real, local_unknown_real) +=
                        G_r.real();
                    }
                    else if (c == 1)
                    {
                      jacobian(local_eqn_real, local_unknown_real) +=
                        G_z.real();
                    }
                    else if (c == 2)
                    {
                      jacobian(local_eqn_real, local_unknown_real) +=
                        G_theta.real();
                    }
                  }

                  // If the imaginary part of the local unknown is not pinned
                  if (local_unknown_imag >= 0)
                  {
                    if (c == 0)
                    {
                      jacobian(local_eqn_real, local_unknown_imag) +=
                        -G_r.imag();
                    }
                    else if (c == 1)
                    {
                      jacobian(local_eqn_real, local_unknown_imag) +=
                        -G_z.imag();
                    }
                    else if (c == 2)
                    {
                      jacobian(local_eqn_real, local_unknown_imag) +=
                        -G_theta.imag();
                    }
                  }
                }
              }
            } // End of jacobian calculation

          } // End of if not boundary condition for real eqn


          // Get the IMAG equation number
          local_eqn_imag = this->nodal_local_eqn(l, u_nodal_index[a].imag());

          /*IF it's not a boundary condition*/
          if (local_eqn_imag >= 0)
          {
            // Acceleration and body force
            residuals[local_eqn_imag] +=
              (-omega_sq.real() * interpolated_u[a].imag() -
               omega_sq.imag() * interpolated_u[a].real() - b[a].imag()) *
              psi(l) * r * W;

            // Three components of the stress divergence term:
            // a=0: r; a=1: z; a=2: theta

            // Imag part of r-equation
            if (a == 0)
            {
              residuals[local_eqn_imag] +=
                (mu.real() *
                   (2 * durdr.imag() * dpsidx(l, 0) +
                    n * psi(l) / r *
                      (n * ur.imag() / r - duthdr.real() + uth.real() / r) +
                    dpsidx(l, 1) * (durdz.imag() + duzdr.imag()) +
                    2 * psi(l) / pow(r, 2) * (ur.imag() + n * uth.real())) +
                 mu.imag() *
                   (2.0 * durdr.real() * dpsidx(l, 0) +
                    n * psi(l) / r *
                      (n * ur.real() / r + duthdr.imag() - uth.imag() / r) +
                    dpsidx(l, 1) * (durdz.real() + duzdr.real()) +
                    2.0 * psi(l) / pow(r, 2) * (ur.real() - n * uth.imag())) +
                 lambda.real() *
                   (durdr.imag() + ur.imag() / r + n / r * uth.real() +
                    duzdz.imag()) *
                   (dpsidx(l, 0) + psi(l) / r) +
                 lambda.imag() *
                   (durdr.real() + ur.real() / r - n / r * uth.imag() +
                    duzdz.real()) *
                   (dpsidx(l, 0) + psi(l) / r)) *
                r * W;
            }
            // Imag part of z-equation
            else if (a == 1)
            {
              residuals[local_eqn_imag] +=
                (mu.real() *
                   (dpsidx(l, 0) * (durdz.imag() + duzdr.imag()) +
                    n * psi(l) / r * (n * uz.imag() / r - duthdz.real()) +
                    2 * duzdz.imag() * dpsidx(l, 1)) +
                 mu.imag() *
                   (dpsidx(l, 0) * (durdz.real() + duzdr.real()) +
                    n * psi(l) / r * (n * uz.real() / r + duthdz.imag()) +
                    2.0 * duzdz.real() * dpsidx(l, 1)) +
                 lambda.real() *
                   (durdr.imag() + ur.imag() / r + n / r * uth.real() +
                    duzdz.imag()) *
                   dpsidx(l, 1) +
                 lambda.imag() *
                   (durdr.real() + ur.real() / r - n / r * uth.imag() +
                    duzdz.real()) *
                   dpsidx(l, 1)) *
                r * W;
            }
            // Imag part of theta-equation
            else if (a == 2)
            {
              residuals[local_eqn_imag] +=
                (mu.real() *
                   ((duthdr.imag() - uth.imag() / r + n * ur.real() / r) *
                      (dpsidx(l, 0) - psi(l) / r) +
                    2.0 * n * psi(l) / pow(r, 2.) *
                      (n * uth.imag() - ur.real()) +
                    dpsidx(l, 1) * (duthdz.imag() + n / r * uz.real())) +
                 mu.imag() *
                   ((duthdr.real() - uth.real() / r - n * ur.imag() / r) *
                      (dpsidx(l, 0) - psi(l) / r) +
                    2.0 * n * psi(l) / pow(r, 2) *
                      (n * uth.real() + ur.imag()) +
                    dpsidx(l, 1) * (duthdz.real() - n / r * uz.imag())) +
                 lambda.real() *
                   (-durdr.real() - ur.real() / r + n / r * uth.imag() -
                    duzdz.real()) *
                   n * psi(l) / r +
                 lambda.imag() *
                   (durdr.imag() + ur.imag() / r + n / r * uth.real() +
                    duzdz.imag()) *
                   n * psi(l) / r) *
                r * W;
            }
            // error: a should be 0, 1 or 2
            else
            {
              throw OomphLibError("a should equal 0, 1 or 2",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }

            // Jacobian entries
            if (flag)
            {
              // Loop over the displacement basis functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // define complex terms used to obtain entries of current row in
                // the Jacobian:

                // terms for rows of Jacobian matrix corresponding to r-equation
                if (a == 0)
                {
                  G_r = (mu * (2.0 * dpsidx(l2, 0) * dpsidx(l, 0) +
                               (nsq + 2.0) / rsq * psi(l2) * psi(l) +
                               dpsidx(l2, 1) * dpsidx(l, 1)) +
                         lambda * (dpsidx(l2, 0) + psi(l2) / r) *
                           (dpsidx(l, 0) + psi(l) / r) -
                         omega_sq * psi(l2) * psi(l)) *
                        r * W;

                  G_z = (mu * dpsidx(l2, 0) * dpsidx(l, 1) +
                         lambda * dpsidx(l2, 1) * (dpsidx(l, 0) + psi(l) / r)) *
                        r * W;

                  G_theta = (-I * (mu * (n / r * dpsidx(l2, 0) * psi(l) -
                                         3.0 * n / rsq * psi(l2) * psi(l)) -
                                   lambda * n / r * psi(l2) *
                                     (dpsidx(l, 0) + psi(l) / r))) *
                            r * W;
                }
                // terms for rows of Jacobian matrix corresponding to z-equation
                else if (a == 1)
                {
                  G_r =
                    (mu * dpsidx(l2, 1) * dpsidx(l, 0) +
                     lambda * (dpsidx(l2, 0) + psi(l2) / r) * dpsidx(l, 1)) *
                    r * W;

                  G_z = (mu * (dpsidx(l2, 0) * dpsidx(l, 0) +
                               nsq / rsq * psi(l2) * psi(l) +
                               2.0 * dpsidx(l2, 1) * dpsidx(l, 1)) +
                         lambda * dpsidx(l2, 1) * dpsidx(l, 1) -
                         omega_sq * psi(l2) * psi(l)) *
                        r * W;

                  G_theta = (-I * (mu * n / r * dpsidx(l2, 1) * psi(l) -
                                   lambda * n / r * dpsidx(l, 1) * psi(l2))) *
                            r * W;
                }
                // terms for rows of Jacobian matrix corresponding to
                // theta-equation
                else if (a == 2)
                {
                  G_r = (-I * (mu * (-n / r * psi(l2) * dpsidx(l, 0) +
                                     3.0 * n / rsq * psi(l2) * psi(l)) +
                               lambda * n / r * (dpsidx(l2, 0) + psi(l2) / r) *
                                 psi(l))) *
                        r * W;

                  G_z = (-I * (-mu * n / r * psi(l2) * dpsidx(l, 1) +
                               lambda * n / r * dpsidx(l2, 1) * psi(l))) *
                        r * W;

                  G_theta = (mu * ((dpsidx(l2, 0) - psi(l2) / r) *
                                     (dpsidx(l, 0) - psi(l) / r) +
                                   2.0 * nsq / rsq * psi(l2) * psi(l) +
                                   dpsidx(l2, 1) * dpsidx(l, 1)) +
                             lambda * nsq / rsq * psi(l2) * psi(l) -
                             omega_sq * psi(l2) * psi(l)) *
                            r * W;
                }

                // Loop over the displacement components
                for (unsigned c = 0; c < 3; c++)
                {
                  // Get real and imaginary parts of local unknown
                  local_unknown_real =
                    this->nodal_local_eqn(l2, u_nodal_index[c].real());
                  local_unknown_imag =
                    this->nodal_local_eqn(l2, u_nodal_index[c].imag());

                  // If the real part of the local unknown is not pinned
                  if (local_unknown_real >= 0)
                  {
                    if (c == 0)
                    {
                      jacobian(local_eqn_imag, local_unknown_real) +=
                        G_r.imag();
                    }
                    else if (c == 1)
                    {
                      jacobian(local_eqn_imag, local_unknown_real) +=
                        G_z.imag();
                    }
                    else if (c == 2)
                    {
                      jacobian(local_eqn_imag, local_unknown_real) +=
                        G_theta.imag();
                    }
                  }

                  // If the imaginary part of the local unknown is not pinned
                  if (local_unknown_imag >= 0)
                  {
                    if (c == 0)
                    {
                      jacobian(local_eqn_imag, local_unknown_imag) +=
                        G_r.real();
                    }
                    else if (c == 1)
                    {
                      jacobian(local_eqn_imag, local_unknown_imag) +=
                        G_z.real();
                    }
                    else if (c == 2)
                    {
                      jacobian(local_eqn_imag, local_unknown_imag) +=
                        G_theta.real();
                    }
                  }
                }
              }
            } // End of jacobian calculation

          } // End of if not boundary condition for imag eqn

        } // End of loop over coordinate directions
      } // End of loop over shape functions
    } // End of loop over integration points
  }

  //=======================================================================
  /// Output exact solution  r,z, u_r_real, u_z_real, ..., u_theta_imag
  //=======================================================================
  void TimeHarmonicFourierDecomposedLinearElasticityEquations::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordintes
    Vector<double> x(2);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Exact solution Vector
    Vector<double> exact_soln(6);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output r,z,...,u_exact,...
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }
      for (unsigned i = 0; i < 6; i++)
      {
        outfile << exact_soln[i] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  //=======================================================================
  /// Output: r,z, u_r_real, u_z_real, ..., u_theta_imag
  //=======================================================================
  void TimeHarmonicFourierDecomposedLinearElasticityEquations::output(
    std::ostream& outfile, const unsigned& nplot)
  {
    // Set output Vector
    Vector<double> s(2);
    Vector<double> x(2);
    Vector<std::complex<double>> u(3);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Get Eulerian coordinates and displacements
      this->interpolated_x(s, x);
      this->interpolated_u_time_harmonic_fourier_decomposed_linear_elasticity(
        s, u);

      // Output the r,z,..
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output real part of displacement
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << u[i].real() << " ";
      }

      // Output imag part of displacement
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << u[i].imag() << " ";
      }

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }


  //=======================================================================
  /// C-style output:r,z, u_r_real, u_z_real, ..., u_theta_imag
  //=======================================================================
  void TimeHarmonicFourierDecomposedLinearElasticityEquations::output(
    FILE* file_pt, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Tecplot header info
    fprintf(file_pt, "%s", this->tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        fprintf(file_pt, "%g ", this->interpolated_x(s, i));
      }

      // Displacement
      for (unsigned i = 0; i < 3; i++)
      {
        fprintf(
          file_pt,
          "%g ",
          this
            ->interpolated_u_time_harmonic_fourier_decomposed_linear_elasticity(
              s, i)
            .real());
      }
      for (unsigned i = 0; i < 3; i++)
      {
        fprintf(
          file_pt,
          "%g ",
          this
            ->interpolated_u_time_harmonic_fourier_decomposed_linear_elasticity(
              s, i)
            .imag());
      }
    }
    fprintf(file_pt, "\n");

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(file_pt, nplot);
  }

  //======================================================================
  /// Validate against exact solution
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points and compute L2 error
  /// and L2 norm of velocity solution over element.
  //=======================================================================
  void TimeHarmonicFourierDecomposedLinearElasticityEquations::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordinates
    Vector<double> x(2);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    outfile << "ZONE" << std::endl;

    // Exact solution Vector (u_r_real, u_z_real, ..., u_theta_imag)
    Vector<double> exact_soln(6);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = this->J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Displacement error
      for (unsigned i = 0; i < 3; i++)
      {
        norm += (exact_soln[i] * exact_soln[i] +
                 exact_soln[i + 3] * exact_soln[i + 3]) *
                W;
        error +=
          ((exact_soln[i] -
            this
              ->interpolated_u_time_harmonic_fourier_decomposed_linear_elasticity(
                s, i)
              .real()) *
             (exact_soln[i] -
              this
                ->interpolated_u_time_harmonic_fourier_decomposed_linear_elasticity(
                  s, i)
                .real()) +
           (exact_soln[i + 3] -
            this
              ->interpolated_u_time_harmonic_fourier_decomposed_linear_elasticity(
                s, i)
              .imag()) *
             (exact_soln[i + 3] -
              this
                ->interpolated_u_time_harmonic_fourier_decomposed_linear_elasticity(
                  s, i)
                .imag())) *
          W;
      }


      // Output r,z coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output ur_error, uz_error, uth_error
      for (unsigned i = 0; i < 3; i++)
      {
        outfile
          << exact_soln[i] -
               this
                 ->interpolated_u_time_harmonic_fourier_decomposed_linear_elasticity(
                   s, i)
                 .real()
          << " ";
      }
      for (unsigned i = 0; i < 3; i++)
      {
        outfile
          << exact_soln[i + 3] -
               this
                 ->interpolated_u_time_harmonic_fourier_decomposed_linear_elasticity(
                   s, i)
                 .imag()
          << " ";
      }
      outfile << std::endl;
    }
  }

  // Instantiate required elements
  template class QTimeHarmonicFourierDecomposedLinearElasticityElement<2>;
  template class QTimeHarmonicFourierDecomposedLinearElasticityElement<3>;
  template class QTimeHarmonicFourierDecomposedLinearElasticityElement<4>;


} // namespace oomph
