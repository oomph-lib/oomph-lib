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
// Non-inline functions for Kirchhoff Love shell elements
#include "shell_elements.h"


namespace oomph
{
  //====================================================================
  /// Default value for the Poisson ratio: 0.49 for no particular reason
  //====================================================================
  double KirchhoffLoveShellEquations::Default_nu_value = 0.49;

  //=======================================================================
  /// Static default value for timescale ratio (1.0 for natural scaling)
  //=======================================================================
  double KirchhoffLoveShellEquations::Default_lambda_sq_value = 1.0;

  //====================================================================
  /// Static default value for non-dim wall thickness (1/20)
  //===================================================================
  double KirchhoffLoveShellEquations::Default_h_value = 0.05;

  //=======================================================================
  /// Default load function (zero traction)
  //=======================================================================
  void KirchhoffLoveShellEquations::Zero_traction_fct(const Vector<double>& xi,
                                                      const Vector<double>& x,
                                                      const Vector<double>& N,
                                                      Vector<double>& load)
  {
    unsigned n_dim = load.size();
    for (unsigned i = 0; i < n_dim; i++)
    {
      load[i] = 0.0;
    }
  }


  //=======================================================================
  /// Static default for prestress (set to zero)
  //=======================================================================
  double KirchhoffLoveShellEquations::Zero_prestress = 0.0;


  //======================================================================
  /// Get normal to the shell
  //=====================================================================
  void KirchhoffLoveShellEquations::get_normal(const Vector<double>& s,
                                               Vector<double>& N)
  {
    // Set the dimension of the coordinates
    const unsigned n_dim = 3;

    // Set the number of lagrangian coordinates
    const unsigned n_lagrangian = 2;

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Find out how many positional dofs there are
    const unsigned n_position_type = nnodal_position_type();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, n_lagrangian);

    // Could/Should we make this more general?
    DShape d2psidxi(n_node, n_position_type, 3);


    // Call the derivatives of the shape functions:
    // d2psidxi(i,0) = \f$ \partial^2 \psi_j / \partial^2 \xi_0^2 \f$
    // d2psidxi(i,1) = \f$ \partial^2 \psi_j / \partial^2 \xi_1^2 \f$
    // d2psidxi(i,2) = \f$ \partial^2 \psi_j/\partial \xi_0 \partial \xi_1 \f$
    (void)dshape_lagrangian(s, psi, dpsidxi);

    // Calculate local values of lagrangian position and
    // the derivative of global position wrt lagrangian coordinates
    DenseMatrix<double> interpolated_A(2, 3);

    // Initialise to zero
    for (unsigned i = 0; i < n_dim; i++)
    {
      for (unsigned j = 0; j < n_lagrangian; j++)
      {
        interpolated_A(j, i) = 0.0;
      }
    }

    // Calculate displacements and derivatives
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over positional dofs
      for (unsigned k = 0; k < n_position_type; k++)
      {
        // Loop over displacement components (deformed position)
        for (unsigned i = 0; i < n_dim; i++)
        {
          double x_value = raw_nodal_position_gen(l, k, i);

          // Loop over derivative directions
          for (unsigned j = 0; j < n_lagrangian; j++)
          {
            interpolated_A(j, i) += x_value * dpsidxi(l, k, j);
          }
        }
      }
    }

    // Unscaled normal
    N[0] = (interpolated_A(0, 1) * interpolated_A(1, 2) -
            interpolated_A(0, 2) * interpolated_A(1, 1));
    N[1] = (interpolated_A(0, 2) * interpolated_A(1, 0) -
            interpolated_A(0, 0) * interpolated_A(1, 2));
    N[2] = (interpolated_A(0, 0) * interpolated_A(1, 1) -
            interpolated_A(0, 1) * interpolated_A(1, 0));

    // Normalise
    double norm = 1.0 / sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
    for (unsigned i = 0; i < 3; i++)
    {
      N[i] *= norm;
    }
  }


  //======================================================================
  /// Get strain and bending tensors
  //=====================================================================
  std::pair<double, double> KirchhoffLoveShellEquations::get_strain_and_bend(
    const Vector<double>& s, DenseDoubleMatrix& strain, DenseDoubleMatrix& bend)
  {
    // Set the dimension of the coordinates
    const unsigned n_dim = 3;

    // Set the number of lagrangian coordinates
    const unsigned n_lagrangian = 2;

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Find out how many positional dofs there are
    const unsigned n_position_type = nnodal_position_type();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, n_lagrangian);

    // Could/Should we make this more general?
    DShape d2psidxi(n_node, n_position_type, 3);


    // Call the derivatives of the shape functions:
    // d2psidxi(i,0) = \f$ \partial^2 \psi_j / \partial^2 \xi_0^2 \f$
    // d2psidxi(i,1) = \f$ \partial^2 \psi_j / \partial^2 \xi_1^2 \f$
    // d2psidxi(i,2) = \f$ \partial^2 \psi_j/\partial \xi_0 \partial \xi_1 \f$
    (void)d2shape_lagrangian(s, psi, dpsidxi, d2psidxi);

    // Calculate local values of lagrangian position and
    // the derivative of global position wrt lagrangian coordinates
    Vector<double> interpolated_xi(2, 0.0), interpolated_x(3, 0.0);
    Vector<double> accel(3, 0.0);
    DenseMatrix<double> interpolated_A(2, 3);
    DenseMatrix<double> interpolated_dAdxi(3);

    // Initialise to zero
    for (unsigned i = 0; i < n_dim; i++)
    {
      for (unsigned j = 0; j < n_lagrangian; j++)
      {
        interpolated_A(j, i) = 0.0;
      }
      for (unsigned j = 0; j < 3; j++)
      {
        interpolated_dAdxi(i, j) = 0.0;
      }
    }

    // Calculate displacements and derivatives
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over positional dofs
      for (unsigned k = 0; k < n_position_type; k++)
      {
        // Loop over lagrangian coordinate directions
        for (unsigned i = 0; i < n_lagrangian; i++)
        {
          interpolated_xi[i] +=
            raw_lagrangian_position_gen(l, k, i) * psi(l, k);
        }

        // Loop over displacement components (deformed position)
        for (unsigned i = 0; i < n_dim; i++)
        {
          double x_value = raw_nodal_position_gen(l, k, i);

          // Calculate interpolated (deformed) position
          interpolated_x[i] += x_value * psi(l, k);

          // Calculate the acceleration
          accel[i] += raw_dnodal_position_gen_dt(2, l, k, i) * psi(l, k);

          // Loop over derivative directions
          for (unsigned j = 0; j < n_lagrangian; j++)
          {
            interpolated_A(j, i) += x_value * dpsidxi(l, k, j);
          }
          // Loop over the second derivative directions
          for (unsigned j = 0; j < 3; j++)
          {
            interpolated_dAdxi(j, i) += x_value * d2psidxi(l, k, j);
          }
        }
      }
    }

    // Get the values of the undeformed parameters
    Vector<double> interpolated_R(3);
    DenseMatrix<double> interpolated_a(2, 3);
    RankThreeTensor<double> dadxi(2, 2, 3);

    // Read out the undeformed position from the geometric object
    Undeformed_midplane_pt->d2position(
      interpolated_xi, interpolated_R, interpolated_a, dadxi);

    // Copy the values of the second derivatives into a slightly
    // different data structure, where the mixed derivatives [0][1] and [1][0]
    // are given by the single index [2]
    DenseMatrix<double> interpolated_dadxi(3);
    for (unsigned i = 0; i < 3; i++)
    {
      interpolated_dadxi(0, i) = dadxi(0, 0, i);
      interpolated_dadxi(1, i) = dadxi(1, 1, i);
      interpolated_dadxi(2, i) = dadxi(0, 1, i);
    }

    // Declare and calculate the undeformed and deformed metric tensor
    double a[2][2], A[2][2], aup[2][2], Aup[2][2];

    // Assign values of A and gamma
    for (unsigned al = 0; al < 2; al++)
    {
      for (unsigned be = 0; be < 2; be++)
      {
        // Initialise a(al,be) and A(al,be) to zero
        a[al][be] = 0.0;
        A[al][be] = 0.0;
        // Now calculate the dot product
        for (unsigned k = 0; k < 3; k++)
        {
          a[al][be] += interpolated_a(al, k) * interpolated_a(be, k);
          A[al][be] += interpolated_A(al, k) * interpolated_A(be, k);
        }
        // Calculate strain tensor
        strain(al, be) = 0.5 * (A[al][be] - a[al][be]);
      }
    }

    // Calculate the contravariant metric tensor
    double adet = calculate_contravariant(a, aup);
    double Adet = calculate_contravariant(A, Aup);

    // Square roots are expensive, so let's do them once here
    double sqrt_adet = sqrt(adet);
    double sqrt_Adet = sqrt(Adet);

    // Calculate the normal Vectors
    Vector<double> n(3), N(3);
    n[0] = (interpolated_a(0, 1) * interpolated_a(1, 2) -
            interpolated_a(0, 2) * interpolated_a(1, 1)) /
           sqrt_adet;
    n[1] = (interpolated_a(0, 2) * interpolated_a(1, 0) -
            interpolated_a(0, 0) * interpolated_a(1, 2)) /
           sqrt_adet;
    n[2] = (interpolated_a(0, 0) * interpolated_a(1, 1) -
            interpolated_a(0, 1) * interpolated_a(1, 0)) /
           sqrt_adet;
    N[0] = (interpolated_A(0, 1) * interpolated_A(1, 2) -
            interpolated_A(0, 2) * interpolated_A(1, 1)) /
           sqrt_Adet;
    N[1] = (interpolated_A(0, 2) * interpolated_A(1, 0) -
            interpolated_A(0, 0) * interpolated_A(1, 2)) /
           sqrt_Adet;
    N[2] = (interpolated_A(0, 0) * interpolated_A(1, 1) -
            interpolated_A(0, 1) * interpolated_A(1, 0)) /
           sqrt_Adet;


    // Calculate the curvature tensors
    double b[2][2], B[2][2];

    b[0][0] = n[0] * interpolated_dadxi(0, 0) +
              n[1] * interpolated_dadxi(0, 1) + n[2] * interpolated_dadxi(0, 2);

    // Off-diagonal terms are the same
    b[0][1] = b[1][0] = n[0] * interpolated_dadxi(2, 0) +
                        n[1] * interpolated_dadxi(2, 1) +
                        n[2] * interpolated_dadxi(2, 2);

    b[1][1] = n[0] * interpolated_dadxi(1, 0) +
              n[1] * interpolated_dadxi(1, 1) + n[2] * interpolated_dadxi(1, 2);

    // Deformed curvature tensor
    B[0][0] = N[0] * interpolated_dAdxi(0, 0) +
              N[1] * interpolated_dAdxi(0, 1) + N[2] * interpolated_dAdxi(0, 2);

    // Off-diagonal terms are the same
    B[0][1] = B[1][0] = N[0] * interpolated_dAdxi(2, 0) +
                        N[1] * interpolated_dAdxi(2, 1) +
                        N[2] * interpolated_dAdxi(2, 2);

    B[1][1] = N[0] * interpolated_dAdxi(1, 0) +
              N[1] * interpolated_dAdxi(1, 1) + N[2] * interpolated_dAdxi(1, 2);

    // Set up the change of curvature tensor
    for (unsigned i = 0; i < 2; i++)
    {
      for (unsigned j = 0; j < 2; j++)
      {
        bend(i, j) = b[i][j] - B[i][j];
      }
    }

    // Return undeformed and deformed determinant of in-plane metric
    // tensor
    return std::make_pair(adet, Adet);
  }


  //====================================================================
  /// Return the residuals for the equations of KL shell theory
  //====================================================================
  void KirchhoffLoveShellEquations::fill_in_contribution_to_residuals_shell(
    Vector<double>& residuals)
  {
    // Set the dimension of the coordinates
    const unsigned n_dim = 3;

    // Set the number of lagrangian coordinates
    const unsigned n_lagrangian = 2;

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Find out how many positional dofs there are
    const unsigned n_position_type = nnodal_position_type();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, n_lagrangian);

    // Could/Should we make this more general?
    DShape d2psidxi(n_node, n_position_type, 3);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Get Physical Variables from Element
    // External pressure, Poisson ratio , and non-dim thickness h
    const double nu_cached = nu();
    const double h_cached = h();
    const double lambda_sq_cached = lambda_sq();

    // Integers used to store the local equation numbers
    int local_eqn = 0;

    const double bending_scale = (1.0 / 12.0) * h_cached * h_cached;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape functions:
      // d2psidxi(i,0) = \f$ \partial^2 \psi_j / \partial^2 \xi_0^2 \f$
      // d2psidxi(i,1) = \f$ \partial^2 \psi_j / \partial^2 \xi_1^2 \f$
      // d2psidxi(i,2) = \f$ \partial^2 \psi_j/\partial \xi_0 \partial \xi_1 \f$
      double J = d2shape_lagrangian_at_knot(ipt, psi, dpsidxi, d2psidxi);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of lagrangian position and
      // the derivative of global position wrt lagrangian coordinates
      Vector<double> interpolated_xi(2, 0.0), interpolated_x(3, 0.0);
      Vector<double> accel(3, 0.0);
      DenseMatrix<double> interpolated_A(2, 3);
      DenseMatrix<double> interpolated_dAdxi(3);

      // Initialise to zero
      for (unsigned i = 0; i < n_dim; i++)
      {
        for (unsigned j = 0; j < n_lagrangian; j++)
        {
          interpolated_A(j, i) = 0.0;
        }
        for (unsigned j = 0; j < 3; j++)
        {
          interpolated_dAdxi(i, j) = 0.0;
        }
      }

      // Calculate displacements and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Loop over lagrangian coordinate directions
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            interpolated_xi[i] +=
              raw_lagrangian_position_gen(l, k, i) * psi(l, k);
          }

          // Loop over displacement components (deformed position)
          for (unsigned i = 0; i < n_dim; i++)
          {
            double x_value = raw_nodal_position_gen(l, k, i);

            // Calculate interpolated (deformed) position
            interpolated_x[i] += x_value * psi(l, k);

            // Calculate the acceleration
            accel[i] += raw_dnodal_position_gen_dt(2, l, k, i) * psi(l, k);

            // Loop over derivative directions
            for (unsigned j = 0; j < n_lagrangian; j++)
            {
              interpolated_A(j, i) += x_value * dpsidxi(l, k, j);
            }
            // Loop over the second derivative directions
            for (unsigned j = 0; j < 3; j++)
            {
              interpolated_dAdxi(j, i) += x_value * d2psidxi(l, k, j);
            }
          }
        }
      }

      // Get the values of the undeformed parameters
      Vector<double> interpolated_R(3);
      DenseMatrix<double> interpolated_a(2, 3);
      RankThreeTensor<double> dadxi(2, 2, 3);

      // Read out the undeformed position from the geometric object
      Undeformed_midplane_pt->d2position(
        interpolated_xi, interpolated_R, interpolated_a, dadxi);

      // Copy the values of the second derivatives into a slightly
      // different data structure, where the mixed derivatives [0][1] and [1][0]
      // are given by the single index [2]
      DenseMatrix<double> interpolated_dadxi(3);
      for (unsigned i = 0; i < 3; i++)
      {
        interpolated_dadxi(0, i) = dadxi(0, 0, i);
        interpolated_dadxi(1, i) = dadxi(1, 1, i);
        interpolated_dadxi(2, i) = dadxi(0, 1, i);
      }

      // Declare and calculate the undeformed and deformed metric tensor,
      // the strain tensor
      double a[2][2], A[2][2], aup[2][2], Aup[2][2], gamma[2][2];

      // Assign values of A and gamma
      for (unsigned al = 0; al < 2; al++)
      {
        for (unsigned be = 0; be < 2; be++)
        {
          // Initialise a(al,be) and A(al,be) to zero
          a[al][be] = 0.0;
          A[al][be] = 0.0;
          // Now calculate the dot product
          for (unsigned k = 0; k < 3; k++)
          {
            a[al][be] += interpolated_a(al, k) * interpolated_a(be, k);
            A[al][be] += interpolated_A(al, k) * interpolated_A(be, k);
          }
          // Calculate strain tensor
          gamma[al][be] = 0.5 * (A[al][be] - a[al][be]);
        }
      }

      // Calculate the contravariant metric tensor
      double adet = calculate_contravariant(a, aup);
      double Adet = calculate_contravariant(A, Aup);

      // Square roots are expensive, so let's do them once here
      double sqrt_adet = sqrt(adet);
      double sqrt_Adet = sqrt(Adet);

      // Calculate the normal Vectors
      Vector<double> n(3), N(3);
      n[0] = (interpolated_a(0, 1) * interpolated_a(1, 2) -
              interpolated_a(0, 2) * interpolated_a(1, 1)) /
             sqrt_adet;
      n[1] = (interpolated_a(0, 2) * interpolated_a(1, 0) -
              interpolated_a(0, 0) * interpolated_a(1, 2)) /
             sqrt_adet;
      n[2] = (interpolated_a(0, 0) * interpolated_a(1, 1) -
              interpolated_a(0, 1) * interpolated_a(1, 0)) /
             sqrt_adet;
      N[0] = (interpolated_A(0, 1) * interpolated_A(1, 2) -
              interpolated_A(0, 2) * interpolated_A(1, 1)) /
             sqrt_Adet;
      N[1] = (interpolated_A(0, 2) * interpolated_A(1, 0) -
              interpolated_A(0, 0) * interpolated_A(1, 2)) /
             sqrt_Adet;
      N[2] = (interpolated_A(0, 0) * interpolated_A(1, 1) -
              interpolated_A(0, 1) * interpolated_A(1, 0)) /
             sqrt_Adet;


      // Calculate the curvature tensors
      double b[2][2], B[2][2];

      b[0][0] = n[0] * interpolated_dadxi(0, 0) +
                n[1] * interpolated_dadxi(0, 1) +
                n[2] * interpolated_dadxi(0, 2);

      // Off-diagonal terms are the same
      b[0][1] = b[1][0] = n[0] * interpolated_dadxi(2, 0) +
                          n[1] * interpolated_dadxi(2, 1) +
                          n[2] * interpolated_dadxi(2, 2);

      b[1][1] = n[0] * interpolated_dadxi(1, 0) +
                n[1] * interpolated_dadxi(1, 1) +
                n[2] * interpolated_dadxi(1, 2);

      // Deformed curvature tensor
      B[0][0] = N[0] * interpolated_dAdxi(0, 0) +
                N[1] * interpolated_dAdxi(0, 1) +
                N[2] * interpolated_dAdxi(0, 2);

      // Off-diagonal terms are the same
      B[0][1] = B[1][0] = N[0] * interpolated_dAdxi(2, 0) +
                          N[1] * interpolated_dAdxi(2, 1) +
                          N[2] * interpolated_dAdxi(2, 2);

      B[1][1] = N[0] * interpolated_dAdxi(1, 0) +
                N[1] * interpolated_dAdxi(1, 1) +
                N[2] * interpolated_dAdxi(1, 2);

      // Set up the change of curvature tensor
      double kappa[2][2];

      for (unsigned i = 0; i < 2; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          kappa[i][j] = b[i][j] - B[i][j];
        }
      }

      // Calculate the plane stress stiffness tensor
      double Et[2][2][2][2];

      // Premultiply some constants
      // NOTE C2 has 1.0-Nu in the denominator
      // double C1 = 1.0/(2.0*(1.0+nu_cached)), C2
      // = 2.0*nu_cached/(1.0-nu_cached);

      // Some constants
      double C1 = 0.5 * (1.0 - nu_cached);
      double C2 = nu_cached;

      // Loop over first index
      for (unsigned i = 0; i < 2; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          for (unsigned k = 0; k < 2; k++)
          {
            for (unsigned l = 0; l < 2; l++)
            {
              // Et[i][j][k][l] = C1*(aup[i][k]*aup[j][l] + aup[i][l]*aup[j][k]
              //                     + C2*aup[i][j]*aup[k][l]);
              Et[i][j][k][l] =
                C1 * (aup[i][k] * aup[j][l] + aup[i][l] * aup[j][k]) +
                C2 * aup[i][j] * aup[k][l];
            }
          }
        }
      }

      // Define the load Vector
      Vector<double> f(3);
      // Determine the load
      load_vector(ipt, interpolated_xi, interpolated_x, N, f);

      // Scale the load by the bending stiffness scale
      // for(unsigned i=0;i<3;i++){f[i] *=
      // bending_scale/(1.0-nu_cached*nu_cached);}

      // Allocate storage for variations of the the normal vector
      double normal_var[3][3];

      //=====EQUATIONS OF ELASTICITY FROM PRINCIPLE OF VIRTUAL
      // DISPLACEMENTS========

      // Little tensor to handle the mixed derivative terms
      unsigned mix[2][2] = {{0, 2}, {2, 1}};

      // Loop over the number of nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the type of degree of freedom
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Setup components of variations in the normal vector
          normal_var[0][0] = 0.0;
          normal_var[0][1] = (-interpolated_A(1, 2) * dpsidxi(l, k, 0) +
                              interpolated_A(0, 2) * dpsidxi(l, k, 1)) /
                             sqrt_Adet;
          normal_var[0][2] = (interpolated_A(1, 1) * dpsidxi(l, k, 0) -
                              interpolated_A(0, 1) * dpsidxi(l, k, 1)) /
                             sqrt_Adet;

          normal_var[1][0] = (interpolated_A(1, 2) * dpsidxi(l, k, 0) -
                              interpolated_A(0, 2) * dpsidxi(l, k, 1)) /
                             sqrt_Adet;
          normal_var[1][1] = 0.0;
          normal_var[1][2] = (-interpolated_A(1, 0) * dpsidxi(l, k, 0) +
                              interpolated_A(0, 0) * dpsidxi(l, k, 1)) /
                             sqrt_Adet;

          normal_var[2][0] = (-interpolated_A(1, 1) * dpsidxi(l, k, 0) +
                              interpolated_A(0, 1) * dpsidxi(l, k, 1)) /
                             sqrt_Adet;
          normal_var[2][1] = (interpolated_A(1, 0) * dpsidxi(l, k, 0) -
                              interpolated_A(0, 0) * dpsidxi(l, k, 1)) /
                             sqrt_Adet;
          normal_var[2][2] = 0.0;

          // Loop over the coordinate direction
          for (unsigned i = 0; i < 3; i++)
          {
            // Get the local equation number
            local_eqn = position_local_eqn(l, k, i);

            // If it's not a boundary condition
            if (local_eqn >= 0)
            {
              // Temporary to hold a common part of the bending terms
              double bending_dot_product_multiplier =
                ((A[1][1] * interpolated_A(0, i) -
                  A[1][0] * interpolated_A(1, i)) *
                   dpsidxi(l, k, 0) +
                 (A[0][0] * interpolated_A(1, i) -
                  A[0][1] * interpolated_A(0, i)) *
                   dpsidxi(l, k, 1)) /
                Adet;

              // Combine the cross and dot product bending terms
              for (unsigned i2 = 0; i2 < 3; i2++)
              {
                normal_var[i][i2] -= N[i2] * bending_dot_product_multiplier;
              }

              // Add in external forcing
              residuals[local_eqn] -=
                f[i] / h_cached * psi(l, k) * W * sqrt_Adet;

              // Storage for all residual terms that are multipled by the
              // Jacobian and area of the undeformed shell
              // Start with the acceleration term
              double other_residual_terms =
                lambda_sq_cached * accel[i] * psi(l, k);

              // Now need the loops over the Greek indicies
              for (unsigned al = 0; al < 2; al++)
              {
                for (unsigned be = 0; be < 2; be++)
                {
                  // Membrane prestress
                  other_residual_terms += *Prestress_pt(al, be) *
                                          interpolated_A(al, i) *
                                          dpsidxi(l, k, be);

                  // Remaining terms
                  for (unsigned ga = 0; ga < 2; ga++)
                  {
                    for (unsigned de = 0; de < 2; de++)
                    {
                      if (!Ignore_membrane_terms)
                      {
                        // Pure membrane term
                        other_residual_terms +=
                          Et[al][be][ga][de] * gamma[al][be] *
                          interpolated_A(ga, i) * dpsidxi(l, k, de);
                      }

                      // Bending terms
                      other_residual_terms -=
                        bending_scale * Et[al][be][ga][de] * kappa[al][be] *
                        (N[i] * d2psidxi(l, k, mix[ga][de])
                         // Cross and dot product terms
                         +
                         normal_var[i][0] * interpolated_dAdxi(mix[ga][de], 0) +
                         normal_var[i][1] * interpolated_dAdxi(mix[ga][de], 1) +
                         normal_var[i][2] * interpolated_dAdxi(mix[ga][de], 2));
                    }
                  }
                }
              }

              residuals[local_eqn] += other_residual_terms * W * sqrt_adet;
            }
          }
        }
      }

    } // End of loop over the integration points
  }

  //=========================================================================
  /// Return the jacobian is calculated by finite differences by default,
  //=========================================================================
  void KirchhoffLoveShellEquations::fill_in_contribution_to_jacobian(
    Vector<double>& residuals, DenseMatrix<double>& jacobian)
  {
    // Call the element's residuals vector
    fill_in_contribution_to_residuals_shell(residuals);

    // Solve for the consistent acceleration in Newmark scheme?
    if (Solve_for_consistent_newmark_accel_flag)
    {
      fill_in_jacobian_for_newmark_accel(jacobian);
      return;
    }

    // Allocate storage for the full residuals of the element
    unsigned n_dof = ndof();
    Vector<double> full_residuals(n_dof);

    // Call the full residuals
    get_residuals(full_residuals);

    // Get the solid entries in the jacobian using finite differences
    SolidFiniteElement::fill_in_jacobian_from_solid_position_by_fd(
      full_residuals, jacobian);

    // Get the entries from the external data, usually load terms
    SolidFiniteElement::fill_in_jacobian_from_external_by_fd(full_residuals,
                                                             jacobian);
  }

  //=======================================================================
  /// Compute the potential (strain) and kinetic energy of the
  /// element.
  //=======================================================================
  void KirchhoffLoveShellEquations::get_energy(double& pot_en, double& kin_en)
  {
    // Initialise
    pot_en = 0.0;
    kin_en = 0.0;

    // Set the dimension of the coordinates
    const unsigned n_dim = Undeformed_midplane_pt->ndim();

    // Set the number of lagrangian coordinates
    const unsigned n_lagrangian = Undeformed_midplane_pt->nlagrangian();

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Find out how many positional dofs there are
    const unsigned n_position_type = nnodal_position_type();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, n_lagrangian);
    DShape d2psidxi(n_node, n_position_type, 3);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Get Physical Variables from Element
    // External pressure, Poisson ratio , and non-dim thickness h
    const double nu_cached = nu();
    const double h_cached = h();
    const double lambda_sq_cached = lambda_sq();

#ifdef PARANOID
    // Check for non-zero prestress
    if ((prestress(0, 0) != 0) || (prestress(1, 0) != 0) ||
        (prestress(1, 1) != 0))
    {
      std::string error_message =
        "Warning: Not sure if the energy is computed correctly\n";
      error_message += " for nonzero prestress\n";

      oomph_info << error_message << std::endl;

      //    throw OomphLibWarning(error_message,
      //                          "KirchhoffLoveShellEquations::get_energy()",
      //                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Setup storage for various quantities
    Vector<double> veloc(n_dim);
    Vector<double> interpolated_xi(n_lagrangian);
    DenseMatrix<double> interpolated_A(n_lagrangian, n_dim);
    DenseMatrix<double> interpolated_dAdxi(3);

    const double bending_scale = (1.0 / 12.0) * h_cached * h_cached;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape functions:
      double J = d2shape_lagrangian_at_knot(ipt, psi, dpsidxi, d2psidxi);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Initialise values to zero
      for (unsigned i = 0; i < n_dim; i++)
      {
        veloc[i] = 0.0;
        for (unsigned j = 0; j < n_lagrangian; j++)
        {
          interpolated_A(j, i) = 0.0;
        }
        for (unsigned j = 0; j < 3; j++)
        {
          interpolated_dAdxi(i, j) = 0.0;
        }
      }
      for (unsigned j = 0; j < n_lagrangian; j++)
      {
        interpolated_xi[j] = 0.0;
      }

      // Calculate velocity and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Loop over lagrangian coordinate directions
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            interpolated_xi[i] += lagrangian_position_gen(l, k, i) * psi(l, k);
          }

          // Loop over displacement components
          for (unsigned i = 0; i < n_dim; i++)
          {
            veloc[i] += dnodal_position_gen_dt(1, l, k, i) * psi(l, k);

            double x_value = nodal_position_gen(l, k, i);

            // Loop over derivative directions
            for (unsigned j = 0; j < n_lagrangian; j++)
            {
              interpolated_A(j, i) += x_value * dpsidxi(l, k, j);
            }
            // Loop over the second derivative directions
            for (unsigned j = 0; j < 3; j++)
            {
              interpolated_dAdxi(j, i) += x_value * d2psidxi(l, k, j);
            }
          }
        }
      }

      // Get square of veloc
      double veloc_sq = 0;
      for (unsigned i = 0; i < n_dim; i++)
      {
        veloc_sq += veloc[i] * veloc[i];
      }

      // Get the values of the undeformed parameters
      Vector<double> interpolated_r(n_dim);
      DenseMatrix<double> interpolated_a(n_lagrangian, n_dim);
      RankThreeTensor<double> dadxi(n_lagrangian, n_lagrangian, n_dim);

      // Read out the undeformed position from the geometric object
      Undeformed_midplane_pt->d2position(
        interpolated_xi, interpolated_r, interpolated_a, dadxi);

      // Copy the values of the second derivatives into a slightly
      // different data structure, where the mixed derivatives [0][1] and [1][0]
      // are given by the single index [2]
      DenseMatrix<double> interpolated_dadxi(3);
      for (unsigned i = 0; i < 3; i++)
      {
        interpolated_dadxi(0, i) = dadxi(0, 0, i);
        interpolated_dadxi(1, i) = dadxi(1, 1, i);
        interpolated_dadxi(2, i) = dadxi(0, 1, i);
      }

      // Declare and calculate the undeformed and deformed metric tensor,
      // the strain tensor
      double a[2][2], A[2][2], aup[2][2], Aup[2][2], gamma[2][2];

      // Assign values of A and gamma
      for (unsigned al = 0; al < 2; al++)
      {
        for (unsigned be = 0; be < 2; be++)
        {
          // Initialise a(al,be) and A(al,be) to zero
          a[al][be] = 0.0;
          A[al][be] = 0.0;
          // Now calculate the dot product
          for (unsigned k = 0; k < n_dim; k++)
          {
            a[al][be] += interpolated_a(al, k) * interpolated_a(be, k);
            A[al][be] += interpolated_A(al, k) * interpolated_A(be, k);
          }
          // Calculate strain tensor
          gamma[al][be] = 0.5 * (A[al][be] - a[al][be]);
        }
      }

      // Calculate the contravariant metric tensor
      double adet = calculate_contravariant(a, aup);
      double Adet = calculate_contravariant(A, Aup);

      // Square roots are expensive, so let's do them once here
      double sqrt_adet = sqrt(adet);
      double sqrt_Adet = sqrt(Adet);

      // Calculate the normal Vectors
      Vector<double> n(3), N(3);
      n[0] = (interpolated_a(0, 1) * interpolated_a(1, 2) -
              interpolated_a(0, 2) * interpolated_a(1, 1)) /
             sqrt_adet;
      n[1] = (interpolated_a(0, 2) * interpolated_a(1, 0) -
              interpolated_a(0, 0) * interpolated_a(1, 2)) /
             sqrt_adet;
      n[2] = (interpolated_a(0, 0) * interpolated_a(1, 1) -
              interpolated_a(0, 1) * interpolated_a(1, 0)) /
             sqrt_adet;
      N[0] = (interpolated_A(0, 1) * interpolated_A(1, 2) -
              interpolated_A(0, 2) * interpolated_A(1, 1)) /
             sqrt_Adet;
      N[1] = (interpolated_A(0, 2) * interpolated_A(1, 0) -
              interpolated_A(0, 0) * interpolated_A(1, 2)) /
             sqrt_Adet;
      N[2] = (interpolated_A(0, 0) * interpolated_A(1, 1) -
              interpolated_A(0, 1) * interpolated_A(1, 0)) /
             sqrt_Adet;


      // Calculate the curvature tensors
      double b[2][2], B[2][2];

      b[0][0] = n[0] * interpolated_dadxi(0, 0) +
                n[1] * interpolated_dadxi(0, 1) +
                n[2] * interpolated_dadxi(0, 2);

      // Off-diagonal terms are the same
      b[0][1] = b[1][0] = n[0] * interpolated_dadxi(2, 0) +
                          n[1] * interpolated_dadxi(2, 1) +
                          n[2] * interpolated_dadxi(2, 2);

      b[1][1] = n[0] * interpolated_dadxi(1, 0) +
                n[1] * interpolated_dadxi(1, 1) +
                n[2] * interpolated_dadxi(1, 2);

      // Deformed curvature tensor
      B[0][0] = N[0] * interpolated_dAdxi(0, 0) +
                N[1] * interpolated_dAdxi(0, 1) +
                N[2] * interpolated_dAdxi(0, 2);

      // Off-diagonal terms are the same
      B[0][1] = B[1][0] = N[0] * interpolated_dAdxi(2, 0) +
                          N[1] * interpolated_dAdxi(2, 1) +
                          N[2] * interpolated_dAdxi(2, 2);

      B[1][1] = N[0] * interpolated_dAdxi(1, 0) +
                N[1] * interpolated_dAdxi(1, 1) +
                N[2] * interpolated_dAdxi(1, 2);

      // Set up the change of curvature tensor
      double kappa[2][2];

      for (unsigned i = 0; i < 2; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          kappa[i][j] = b[i][j] - B[i][j];
        }
      }

      // Calculate the plane stress stiffness tensor
      double Et[2][2][2][2];

      // Some constants
      double C1 = 0.5 * (1.0 - nu_cached);
      double C2 = nu_cached;

      // Loop over first index
      for (unsigned i = 0; i < 2; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          for (unsigned k = 0; k < 2; k++)
          {
            for (unsigned l = 0; l < 2; l++)
            {
              Et[i][j][k][l] =
                C1 * (aup[i][k] * aup[j][l] + aup[i][l] * aup[j][k]) +
                C2 * aup[i][j] * aup[k][l];
            }
          }
        }
      }

      // Add contributions to potential energy
      double pot_en_cont = 0.0;

      for (unsigned i = 0; i < 2; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          for (unsigned k = 0; k < 2; k++)
          {
            for (unsigned l = 0; l < 2; l++)
            {
              pot_en_cont +=
                Et[i][j][k][l] * (gamma[i][j] * gamma[k][l] +
                                  bending_scale * kappa[i][j] * kappa[k][l]);
            }
          }
        }
      }
      pot_en += pot_en_cont * 0.5 * W * sqrt_adet;

      // Add contribution to kinetic energy
      kin_en += 0.5 * lambda_sq_cached * veloc_sq * W * sqrt_adet;

    } // End of loop over the integration points
  }


  //===================================================================
  /// Get integral of instantaneous rate of work done on
  /// the wall due to the load returned by the virtual
  /// function load_vector_for_rate_of_work_computation(...).
  /// In the current class
  /// the latter function simply de-references the external
  /// load but this can be overloaded in derived classes
  /// (e.g. in FSI elements) to determine the rate of work done
  /// by individual constituents of this load (e.g. the fluid load
  /// in an FSI problem).
  //===================================================================
  double KirchhoffLoveShellEquations::load_rate_of_work()
  {
    // Initialise
    double rate_of_work_integral = 0.0;

    // The number of dimensions
    const unsigned n_dim = 3;

    // The number of lagrangian coordinates
    const unsigned n_lagrangian = 2;

    // The number of nodes
    const unsigned n_node = nnode();

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Find out how many positional dofs there are
    unsigned n_position_type = nnodal_position_type();

    // Vector to hold local coordinate in element
    Vector<double> s(n_lagrangian);

    // Set up storage for shape functions and derivatives of shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, n_lagrangian);

    // Storage for jacobian of mapping from local to Lagrangian coords
    DenseMatrix<double> jacobian(n_lagrangian);

    // Storage for velocity vector
    Vector<double> velocity(n_dim);

    // Storage for load vector
    Vector<double> load(n_dim);

    // Storage for normal vector
    Vector<double> normal(n_dim);

    // Storage for Lagrangian and Eulerian coordinates
    Vector<double> xi(n_lagrangian);
    Vector<double> x(n_dim);

    // Storage for covariant vectors
    DenseMatrix<double> interpolated_A(n_lagrangian, n_dim);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get local coords at knot
      for (unsigned i = 0; i < n_lagrangian; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get shape functions, derivatives, determinant of Jacobian of mapping
      double J = dshape_lagrangian(s, psi, dpsidxi);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Initialise velocity, x and interpolated_A to zero
      for (unsigned i = 0; i < n_dim; i++)
      {
        velocity[i] = 0.0;
        x[i] = 0.0;
        for (unsigned j = 0; j < n_lagrangian; j++)
        {
          interpolated_A(j, i) = 0.0;
        }
      }

      // Initialise xi to zero
      for (unsigned i = 0; i < n_lagrangian; i++)
      {
        xi[i] = 0.0;
      }

      // Calculate velocity, coordinates and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Loop over lagrangian coordinate directions
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            xi[i] += lagrangian_position_gen(l, k, i) * psi(l, k);
          }

          // Loop over displacement components
          for (unsigned i = 0; i < n_dim; i++)
          {
            velocity[i] += dnodal_position_gen_dt(1, l, k, i) * psi(l, k);

            double x_value = nodal_position_gen(l, k, i);
            x[i] += x_value * psi(l, k);

            // Loop over derivative directions
            for (unsigned j = 0; j < n_lagrangian; j++)
            {
              interpolated_A(j, i) += x_value * dpsidxi(l, k, j);
            }
          }
        }
      }

      // Get deformed metric tensor
      double A[2][2];
      for (unsigned al = 0; al < 2; al++)
      {
        for (unsigned be = 0; be < 2; be++)
        {
          // Initialise A(al,be) to zero
          A[al][be] = 0.0;
          // Calculate the dot product
          for (unsigned k = 0; k < n_dim; k++)
          {
            A[al][be] += interpolated_A(al, k) * interpolated_A(be, k);
          }
        }
      }

      // Get determinant of metric tensor
      double A_det = A[0][0] * A[1][1] - A[0][1] * A[1][0];

      // Get outer unit normal
      get_normal(s, normal);

      // Get load vector
      load_vector_for_rate_of_work_computation(ipt, xi, x, normal, load);

      // Local rate of work:
      double rate_of_work = 0.0;
      for (unsigned i = 0; i < n_dim; i++)
      {
        rate_of_work += load[i] * velocity[i];
      }

      // Add rate of work
      rate_of_work_integral += rate_of_work / h() * W * A_det;
    }

    return rate_of_work_integral;
  }


  //===================================================================
  /// The output function: position, veloc and accel.
  //===================================================================
  void HermiteShellElement::output_with_time_dep_quantities(
    std::ostream& outfile, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(2);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;

    // Loop over plot points
    for (unsigned l2 = 0; l2 < n_plot; l2++)
    {
      s[1] = -1.0 + l2 * 2.0 / (n_plot - 1);
      for (unsigned l1 = 0; l1 < n_plot; l1++)
      {
        s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

        // Output
        outfile << interpolated_x(s, 0) << " " << interpolated_x(s, 1) << " "
                << interpolated_x(s, 2) << " " << interpolated_dxdt(s, 0, 1)
                << " " << interpolated_dxdt(s, 1, 1) << " "
                << interpolated_dxdt(s, 2, 1) << " "
                << interpolated_dxdt(s, 0, 2) << " "
                << interpolated_dxdt(s, 1, 2) << " "
                << interpolated_dxdt(s, 2, 2) << " " << std::endl;
      }
    }
    outfile << std::endl;
  }


  //===================================================================
  /// The output function
  //===================================================================
  void HermiteShellElement::output(std::ostream& outfile,
                                   const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(2);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;

    // Loop over plot points
    for (unsigned l2 = 0; l2 < n_plot; l2++)
    {
      s[1] = -1.0 + l2 * 2.0 / (n_plot - 1);
      for (unsigned l1 = 0; l1 < n_plot; l1++)
      {
        s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

        // Output the x and y positions
        outfile << interpolated_x(s, 0) << " " << interpolated_x(s, 1) << " "
                << interpolated_x(s, 2) << " ";
        outfile << interpolated_xi(s, 0) << " " << interpolated_xi(s, 1)
                << std::endl;
      }
    }
    outfile << std::endl;
  }


  //===================================================================
  /// The output function
  //===================================================================
  void HermiteShellElement::output(FILE* file_pt, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(2);

    // Tecplot header info
    fprintf(file_pt, "ZONE I=%i, J=%i \n", n_plot, n_plot);

    // Loop over plot points
    for (unsigned l2 = 0; l2 < n_plot; l2++)
    {
      s[1] = -1.0 + l2 * 2.0 / (n_plot - 1);
      for (unsigned l1 = 0; l1 < n_plot; l1++)
      {
        s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

        // Output the x and y positions
        fprintf(file_pt,
                "%g %g %g ",
                interpolated_x(s, 0),
                interpolated_x(s, 1),
                interpolated_x(s, 2));
        fprintf(
          file_pt, "%g %g \n ", interpolated_xi(s, 0), interpolated_xi(s, 1));
      }
    }
    fprintf(file_pt, "\n");
  }


  //=========================================================================
  /// Define the dposition function. This is used to set no-slip boundary
  /// conditions in some FSI problems.
  //=========================================================================
  void FSIDiagHermiteShellElement::dposition_dlagrangian_at_local_coordinate(
    const Vector<double>& s, DenseMatrix<double>& drdxi) const
  {
#ifdef PARANOID
    if (Undeformed_midplane_pt == 0)
    {
      throw OomphLibError("Undeformed_midplane_pt has not been set",
                          "FSIDiagHermiteShellElement::dposition_dlagrangian_"
                          "at_local_coordinate()",
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the dimension of the global coordinate
    unsigned n_dim = Undeformed_midplane_pt->ndim();

    // Find the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_midplane_pt->nlagrangian();

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, n_lagrangian);

    // Get the derivatives of the shape functions w.r.t local coordinates
    // at this point
    dshape_lagrangian(s, psi, dpsidxi);

    // Initialise the derivatives to zero
    drdxi.initialise(0.0);

    // Loop over the nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over the positional types
      for (unsigned k = 0; k < n_position_type; k++)
      {
        // Loop over the dimensions
        for (unsigned i = 0; i < n_dim; i++)
        {
          // Loop over the lagrangian coordinates
          for (unsigned j = 0; j < n_lagrangian; j++)
          {
            // Add the contribution to the derivative
            drdxi(j, i) += nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
          }
        }
      }
    }
  }


  //=============================================================================
  /// Create a list of pairs for all unknowns in this element,
  /// so that the first entry in each pair contains the global equation
  /// number of the unknown, while the second one contains the number
  /// of the "DOF types" that this unknown is associated with.
  /// (Function can obviously only be called if the equation numbering
  /// scheme has been set up.)
  /// This element is only in charge of the solid dofs.
  //=============================================================================
  void FSIDiagHermiteShellElement::get_dof_numbers_for_unknowns(
    std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // temporary pair (used to store dof lookup prior to being added to list)
    std::pair<unsigned long, unsigned> dof_lookup;

    // number of nodes
    const unsigned n_node = this->nnode();

    // Get the number of position dofs and dimensions at the node
    const unsigned n_position_type = nnodal_position_type();
    const unsigned nodal_dim = nodal_dimension();

    // Integer storage for local unknown
    int local_unknown = 0;

    // Loop over the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // Loop over position dofs
      for (unsigned k = 0; k < n_position_type; k++)
      {
        // Loop over dimension
        for (unsigned i = 0; i < nodal_dim; i++)
        {
          // If the variable is free
          local_unknown = position_local_eqn(n, k, i);

          // ignore pinned values
          if (local_unknown >= 0)
          {
            // store dof lookup in temporary pair: First entry in pair
            // is global equation number; second entry is dof type
            dof_lookup.first = this->eqn_number(local_unknown);
            dof_lookup.second = 0;

            // add to list
            dof_lookup_list.push_front(dof_lookup);
          }
        }
      }
    }
  }


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  ClampedHermiteShellBoundaryConditionElement::
    ClampedHermiteShellBoundaryConditionElement(
      FiniteElement* const& bulk_el_pt, const int& face_index)
    : FaceGeometry<HermiteShellElement>(), FaceElement()
  {
    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index, this);

    // Overwrite default assignments:

    // Number of position types (for automatic local equation numbering)
    // must be re-set to the shell case because the constraint
    // equation makes contributions to all 3 (coordinate directions) x
    // 4 (types of positional dofs) = 12 dofs associated with each node in the
    // shell element. But keep reading...
    set_nnodal_position_type(4);


    // The trouble with the above is that the parametrisation of the
    // element geometry (used in J_eulerian and elsewhere)
    // is based on the 4 (2 nodes x 2 types -- value and slope) shape functions
    // in the 1D element. In general this is done by looping over the
    // nodal position types (2 per node) and referring to the
    // relevant position Data at the nodes indirectly via the
    // bulk_position_type(...) function. This is now screwed up
    // because we've bumped up the number of nodal position types
    // to four! The slightly hacky way around this is resize
    // the lookup scheme underlying the bulk_position_type(...)
    // to four types of degree of freedom per node so it's consistent
    // with the above assignement. However, we only have shape
    // functions associated with the two types of degree of freedom,
    // so (below) we overload shape(...) and dshape_local(...)
    // to that it returns zero for all "superfluous" shape functions.
    bulk_position_type_resize(4);

    // Make some dummy assigments to the position types that
    // correspond to the two superfluous ones -- the values associated
    // with these will never be used since we're setting the
    // associated shape functions to zero!
    bulk_position_type(2) = 0;
    bulk_position_type(3) = 0;

    // Resize normal to clamping plane and initialise for clamping
    // being applied in a plane where z=const.
    Normal_to_clamping_plane.resize(3);
    Normal_to_clamping_plane[0] = 0.0;
    Normal_to_clamping_plane[1] = 0.0;
    Normal_to_clamping_plane[2] = 1.0;

    // We need two additional values (for the two types of the Lagrange
    // multiplier interpolated by the 1D Hermite basis functions)
    // at the element's two nodes
    Vector<unsigned> nadditional_data_values(2);
    nadditional_data_values[0] = 2;
    nadditional_data_values[1] = 2;
    resize_nodes(nadditional_data_values);
  }


  //===========================================================================
  /// Add the element's contribution to its residual vector
  //===========================================================================
  void ClampedHermiteShellBoundaryConditionElement::
    fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
    // Get position vector to and normal vector on wall (from bulk element)
    HermiteShellElement* bulk_el_pt =
      dynamic_cast<HermiteShellElement*>(bulk_element_pt());

    // Local coordinates in bulk
    Vector<double> s_bulk(2);

    // Local coordinate in FaceElement:
    Vector<double> s(1);

    // Normal to shell
    Vector<double> shell_normal(3);
    Vector<double> shell_normal_plus(3);

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Types of (nonzero) shape/test fcts (for Lagrange multiplier) in
    // this element
    unsigned n_type = 2;

    // Integer to store the local equation number
    int local_eqn = 0;

    // Set up memory for the shape functions:

    // Shape fcts
    Shape psi(n_node, n_type);

    // Set # of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get integration point in local coordinates
      s[0] = integral_pt()->knot(ipt, 0);

      // Get shape functions
      shape(s, psi);

      // Get jacobian
      double J = J_eulerian(s);

      // Assemble the Lagrange multiplier
      double lambda = 0.0;
      for (unsigned j = 0; j < n_node; j++)
      {
        for (unsigned k = 0; k < n_type; k++)
        {
          // The (additional) Lagrange multiplier values are stored
          // after those that were created by the bulk elements:
          lambda += node_pt(j)->value(Nbulk_value[j] + k) * psi(j, k);
        }
      }

      // Get vector of coordinates in bulk element
      s_bulk = local_coordinate_in_bulk(s);

      // Get unit normal on bulk shell element
      bulk_el_pt->get_normal(s_bulk, shell_normal);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Here's the actual constraint
      double constraint_residual =
        shell_normal[0] * Normal_to_clamping_plane[0] +
        shell_normal[1] * Normal_to_clamping_plane[1] +
        shell_normal[2] * Normal_to_clamping_plane[2];


      // Assemble residual for Lagrange multiplier:
      //-------------------------------------------

      // Loop over the number of nodes
      for (unsigned j = 0; j < n_node; j++)
      {
        // Loop over the type of degree of freedom
        for (unsigned k = 0; k < n_type; k++)
        {
          // Local eqn number. Recall that the
          // (additional) Lagrange multiplier values are stored
          // after those that were created by the bulk elements:
          local_eqn = nodal_local_eqn(j, Nbulk_value[j] + k);
          if (local_eqn >= 0)
          {
            residuals[local_eqn] += constraint_residual * psi(j, k) * W;
          }
        }
      }


      // Add Lagrange multiplier contribution to bulk equations
      //-------------------------------------------------------

      // Derivatives of constraint equations w.r.t. positional values
      // is done by finite differencing
      double fd_step = 1.0e-8;

      // Loop over the number of nodes
      for (unsigned j = 0; j < n_node; j++)
      {
        Node* nod_pt = node_pt(j);

        // Loop over directions
        for (unsigned i = 0; i < 3; i++)
        {
          // Loop over the type of degree of freedom
          for (unsigned k = 0; k < 4; k++)
          {
            // Local eqn number: Node, type, direction
            local_eqn = position_local_eqn(j, k, i);
            if (local_eqn >= 0)
            {
              // Backup
              double backup = nod_pt->x_gen(k, i);

              // Increment
              nod_pt->x_gen(k, i) += fd_step;

              // Re-evaluate unit normal on bulk shell element
              bulk_el_pt->get_normal(s_bulk, shell_normal_plus);

              // Re-evaluate the constraint
              double constraint_residual_plus =
                shell_normal_plus[0] * Normal_to_clamping_plane[0] +
                shell_normal_plus[1] * Normal_to_clamping_plane[1] +
                shell_normal_plus[2] * Normal_to_clamping_plane[2];

              // Derivarive of constraint w.r.t. to current discrete
              // displacement
              double dres_dx =
                (constraint_residual_plus - constraint_residual) / fd_step;

              // Add to residual
              residuals[local_eqn] += lambda * dres_dx * W;

              // Reset
              nod_pt->x_gen(k, i) = backup;
            }
          }
        }
      }
    } // End of loop over the integration points

    // Loop over integration points
  }


  //===========================================================================
  /// Output function
  //===========================================================================
  void ClampedHermiteShellBoundaryConditionElement::output(
    std::ostream& outfile, const unsigned& n_plot)
  {
    // Get bulk element
    HermiteShellElement* bulk_el_pt =
      dynamic_cast<HermiteShellElement*>(bulk_element_pt());

    // Local coord
    Vector<double> s(1);

    // Coordinate in bulk
    Vector<double> s_bulk(2);

    // Normal vector
    Vector<double> shell_normal(3);

    // # of nodes, # of dof types
    Shape psi(2, 2);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << std::endl;

    // Loop over plot points
    for (unsigned l = 0; l < n_plot; l++)
    {
      s[0] = -1.0 + l * 2.0 / (n_plot - 1);

      // Get vector of coordinates in bulk element
      s_bulk = local_coordinate_in_bulk(s);

      // Get unit normal on bulk shell element
      bulk_el_pt->get_normal(s_bulk, shell_normal);

      // Get shape function
      shape(s, psi);

      // Assemble the Lagrange multiplier
      double lambda = 0.0;
      for (unsigned j = 0; j < 2; j++)
      {
        for (unsigned k = 0; k < 2; k++)
        {
          // The (additional) Lagrange multiplier values are stored
          // after those that were created by the bulk elements:
          lambda += node_pt(j)->value(Nbulk_value[j] + k) * psi(j, k);
        }
      }

      double dot_product = shell_normal[0] * Normal_to_clamping_plane[0] +
                           shell_normal[1] * Normal_to_clamping_plane[1] +
                           shell_normal[2] * Normal_to_clamping_plane[2];

      // Output stuff
      outfile << s[0] << " " << interpolated_xi(s, 0) << " "
              << interpolated_xi(s, 1) << " " << interpolated_x(s, 0) << " "
              << interpolated_x(s, 1) << " " << interpolated_x(s, 2) << " "
              << J_eulerian(s) << " " << bulk_el_pt->interpolated_xi(s_bulk, 0)
              << " " << bulk_el_pt->interpolated_xi(s_bulk, 1) << " "
              << bulk_el_pt->interpolated_x(s_bulk, 0) << " "
              << bulk_el_pt->interpolated_x(s_bulk, 1) << " "
              << bulk_el_pt->interpolated_x(s_bulk, 2) << " " << lambda << " "
              << shell_normal[0] << " " << shell_normal[1] << " "
              << shell_normal[2] << " " << Normal_to_clamping_plane[0] << " "
              << Normal_to_clamping_plane[1] << " "
              << Normal_to_clamping_plane[2] << " " << dot_product << " "
              << std::endl;
    }


    // Initialise error
    const bool check_error = false;
    if (check_error)
    {
      double max_error = 0.0;

      // Set # of integration points
      const unsigned n_intpt = integral_pt()->nweight();

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get integration point in local coordinates
        s[0] = integral_pt()->knot(ipt, 0);

        // Get eulerian jacobian via s
        double J_via_s = J_eulerian(s);

        // Get eulerian jacobian via ipt
        double J_via_ipt = J_eulerian_at_knot(ipt);

        double error = std::fabs(J_via_s - J_via_ipt);
        if (error > max_error) max_error = error;
      }
      if (max_error > 1.0e-14)
      {
        oomph_info << "Warning: Max. error between J_via_s J_via_ipt "
                   << max_error << std::endl;
      }
    }
  }

  //=============================================================================
  /// Create a list of pairs for all unknowns in this element,
  /// so that the first entry in each pair contains the global equation
  /// number of the unknown, while the second one contains the number
  /// of the "DOF types" that this unknown is associated with.
  /// (Function can obviously only be called if the equation numbering
  /// scheme has been set up.)
  /// This element is only in charge of the solid dofs.
  //=============================================================================
  void ClampedHermiteShellBoundaryConditionElement::
    get_dof_numbers_for_unknowns(
      std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // temporary pair (used to store dof lookup prior to being added to list)
    std::pair<unsigned long, unsigned> dof_lookup;

    // number of nodes
    const unsigned n_node = this->nnode();

    // Get the number of position dofs and dimensions at the node
    const unsigned n_position_type = nnodal_position_type();
    const unsigned nodal_dim = nodal_dimension();

    // Integer storage for local unknown
    int local_unknown = 0;

    // Loop over the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      Node* nod_pt = this->node_pt(n);

      // Loop over position dofs
      for (unsigned k = 0; k < n_position_type; k++)
      {
        // Loop over dimension
        for (unsigned i = 0; i < nodal_dim; i++)
        {
          // If the variable is free
          local_unknown = position_local_eqn(n, k, i);

          // ignore pinned values
          if (local_unknown >= 0)
          {
            // store dof lookup in temporary pair: First entry in pair
            // is global equation number; second entry is dof type
            dof_lookup.first = this->eqn_number(local_unknown);
            dof_lookup.second = 0;

            // add to list
            dof_lookup_list.push_front(dof_lookup);
          }
        }
      }

      // Loop over values
      unsigned n_val = nod_pt->nvalue();
      for (unsigned j = 0; j < n_val; j++)
      {
        // If the variable is free
        local_unknown = nodal_local_eqn(n, j);

        // ignore pinned values
        if (local_unknown >= 0)
        {
          // store dof lookup in temporary pair: First entry in pair
          // is global equation number; second entry is dof type
          dof_lookup.first = this->eqn_number(local_unknown);
          dof_lookup.second = 0;

          // add to list
          dof_lookup_list.push_front(dof_lookup);
        }
      }
    }
  }


} // namespace oomph
