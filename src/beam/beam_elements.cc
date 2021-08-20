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
// Non-inline functions for Kirchhoff Love beam elements

#include "beam_elements.h"

namespace oomph
{
  //================================================================
  /// Static default value for 2nd Piola Kirchhoff prestress (zero)
  //================================================================
  double KirchhoffLoveBeamEquations::Default_sigma0_value = 0.0;

  //=====================================================================
  /// Static default value for timescale ratio (1.0 for natural scaling)
  //=====================================================================
  double KirchhoffLoveBeamEquations::Default_lambda_sq_value = 1.0;

  //=========================================================
  /// Static default value for non-dim wall thickness (1/20)
  //=========================================================
  // i.e. the reference thickness 'h_0'
  double KirchhoffLoveBeamEquations::Default_h_value = 0.05;

  //=======================================================================
  /// Default load function (zero traction)
  //=======================================================================
  void KirchhoffLoveBeamEquations::Zero_traction_fct(const Vector<double>& xi,
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
  /// Default wall profile function (constant thickness h_0)
  //=======================================================================
  void KirchhoffLoveBeamEquations::Unit_profile_fct(const Vector<double>& xi,
                                                    const Vector<double>& x,
                                                    double& h_ratio)
  {
    h_ratio = 1.0;
  }


  //=======================================================================
  /// Get position vector to and normal vector on wall
  //=======================================================================
  void KirchhoffLoveBeamEquations::get_normal(const Vector<double>& s,
                                              Vector<double>& r,
                                              Vector<double>& N)
  {
#ifdef PARANOID
    if (N.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "Normal vector should have dimension 2, not" << N.size()
                    << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (r.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "Position vector should have dimension 2, not"
                    << r.size() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    if (s.size() != 1)
    {
      std::ostringstream error_message;
      error_message << "Local coordinate should have dimension 1, not"
                    << s.size() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Set the dimension of the global coordinates
    unsigned n_dim = Undeformed_beam_pt->ndim();

    // Set the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_type = nnodal_position_type();

    // Set up memory for the shape functions:

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_type);

    // # of nodes, # of positional dofs, # of lagrangian coords (for deriv)
    DShape dpsidxi(n_node, n_position_type, n_lagrangian);

    // Call the derivatives of the shape functions w.r.t. Lagrangian coords
    dshape_lagrangian(s, psi, dpsidxi);

    // Base Vector
    DenseMatrix<double> interpolated_A(n_lagrangian, n_dim);

    // Initialise to zero

    // Loop over coordinate directions/components of Vector
    for (unsigned i = 0; i < n_dim; i++)
    {
      r[i] = 0.0;
      // Loop over derivatives/base Vectors (just one here)
      for (unsigned j = 0; j < n_lagrangian; j++)
      {
        interpolated_A(j, i) = 0.0;
      }
    }

    // Loop over directions
    for (unsigned i = 0; i < n_dim; i++)
    {
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_type; k++)
        {
          r[i] += raw_nodal_position_gen(l, k, i) * psi(l, k);

          // Loop over derivative directions (just one here)
          for (unsigned j = 0; j < n_lagrangian; j++)
          {
            interpolated_A(j, i) +=
              raw_nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
          }
        }
      }
    }

    // Calculate the length of the normal vector
    double length = pow(interpolated_A(0, 0), 2) + pow(interpolated_A(0, 1), 2);

    // Calculate the normal
    N[0] = -interpolated_A(0, 1) / sqrt(length);
    N[1] = interpolated_A(0, 0) / sqrt(length);
  }


  //=======================================================================
  /// Get position vector to and non-unit tangent vector on wall: dr/ds
  //=======================================================================
  void KirchhoffLoveBeamEquations::get_non_unit_tangent(const Vector<double>& s,
                                                        Vector<double>& r,
                                                        Vector<double>& drds)
  {
#ifdef PARANOID
    if (drds.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "Tangent vector should have dimension 2, not"
                    << drds.size() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (r.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "Position vector should have dimension 2, not"
                    << r.size() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    if (s.size() != 1)
    {
      std::ostringstream error_message;
      error_message << "Local coordinate should have dimension 1, not"
                    << s.size() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Set the dimension of the global coordinates
    unsigned n_dim = Undeformed_beam_pt->ndim();

    // Set the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_type = nnodal_position_type();

    // Set up memory for the shape functions:

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_type);

    // # of nodes, # of positional dofs, # of lagrangian coords (for deriv)
    DShape dpsids(n_node, n_position_type, n_lagrangian);

    // Call the derivatives of the shape functions w.r.t. local coords
    dshape_local(s, psi, dpsids);

    // Initialise to zero

    // Loop over coordinate directions/components of Vector
    for (unsigned i = 0; i < n_dim; i++)
    {
      r[i] = 0.0;
      drds[i] = 0.0;
    }

    // Loop over directions
    for (unsigned i = 0; i < n_dim; i++)
    {
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_type; k++)
        {
          r[i] += raw_nodal_position_gen(l, k, i) * psi(l, k);
          // deriv w.r.t. to zero-th (and only) local coordiate
          drds[i] += raw_nodal_position_gen(l, k, i) * dpsids(l, k, 0);
        }
      }
    }
  }


  //=======================================================================
  /// Return the residuals for the equations of Kirchhoff-Love beam
  /// theory with linear constitutive equations; if  Solid_ic_pt!=0, we
  /// assign residuals which force the assignement of an initial shape/
  /// veloc/accel to the dofs.
  //=======================================================================
  void KirchhoffLoveBeamEquations::fill_in_contribution_to_residuals_beam(
    Vector<double>& residuals)
  {
    // Set up the initial conditions, if an IC pointer has been set
    if (Solid_ic_pt != 0)
    {
      fill_in_residuals_for_solid_ic(residuals);
      return;
    }

    // Set the dimension of the global coordinates
    const unsigned n_dim = Undeformed_beam_pt->ndim();

    // Set the number of lagrangian coordinates
    const unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Find out how many positional dofs there are
    const unsigned n_position_type = nnodal_position_type();

    // Integer to store the local equation number
    int local_eqn = 0;

    // Setup memory for accelerations
    Vector<double> accel(2);

    // Set up memory for the shape functions:

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_type);

    // # of nodes, # of positional dofs, # of lagrangian coords (for deriv)
    DShape dpsidxi(n_node, n_position_type, n_lagrangian);

    // # of nodes, # of positional dofs, # of derivs)
    DShape d2psidxi(n_node, n_position_type, n_lagrangian);

    // Set # of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Get Physical Variables from Element

    // Thickness h_0/R, axial prestress, timescale ratio (density)
    const double HoR_0 = h(); // i.e. refers to reference thickness 'h_0'
    const double sigma_0 = sigma0();
    const double Lambda_sq = lambda_sq();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape functions w.r.t. Lagrangian coords
      double J = d2shape_lagrangian_at_knot(ipt, psi, dpsidxi, d2psidxi);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of lagrangian position and
      // the derivative of global position wrt lagrangian coordinates
      Vector<double> interpolated_xi(n_lagrangian, 0.0),
        interpolated_x(n_dim, 0.0);
      DenseMatrix<double> interpolated_A(n_lagrangian, n_dim);
      DenseMatrix<double> interpolated_dAdxi(n_lagrangian, n_dim);

      // Initialise to zero
      // Loop over coordinate directions/components of Vector
      for (unsigned i = 0; i < n_dim; i++)
      {
        // Initialise acclerations
        accel[i] = 0.0;
        // Loop over derivatives/base Vectors (just one here)
        for (unsigned j = 0; j < n_lagrangian; j++)
        {
          interpolated_A(j, i) = 0.0;
        }
        // Loop over derivatives of base Vector (just one here)
        for (unsigned j = 0; j < n_lagrangian; j++)
        {
          interpolated_dAdxi(j, i) = 0.0;
        }
      }

      // Calculate displacements, accelerations and spatial derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Loop over Lagrangian coordinate directions [xi_gen[] are the
          // the *gen*eralised Lagrangian coordinates: node, type, direction]
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            interpolated_xi[i] +=
              raw_lagrangian_position_gen(l, k, i) * psi(l, k);
          }

          // Loop over components of the deformed position Vector
          for (unsigned i = 0; i < n_dim; i++)
          {
            interpolated_x[i] += raw_nodal_position_gen(l, k, i) * psi(l, k);

            accel[i] += raw_dnodal_position_gen_dt(2, l, k, i) * psi(l, k);

            // Loop over derivative directions (just one here)
            for (unsigned j = 0; j < n_lagrangian; j++)
            {
              interpolated_A(j, i) +=
                raw_nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
            }

            // Loop over the second derivative directions (just one here)
            for (unsigned j = 0; j < n_lagrangian; j++)
            {
              interpolated_dAdxi(j, i) +=
                raw_nodal_position_gen(l, k, i) * d2psidxi(l, k, j);
            }
          }
        }
      }

      // Setup position Vector and derivatives of undeformed config
      Vector<double> R(n_dim);
      DenseMatrix<double> a(n_lagrangian, n_dim);
      RankThreeTensor<double> dadxi(n_lagrangian, n_lagrangian, n_dim);

      // Get the undeformed geometry
      Undeformed_beam_pt->d2position(interpolated_xi, R, a, dadxi);

      // Declare and calculate the undeformed and deformed metric tensor
      // and the strain tensor (these are just scalars)
      double amet = 0.0, Amet = 0.0;

      // Now calculate the dot product
      for (unsigned k = 0; k < n_dim; k++)
      {
        amet += a(0, k) * a(0, k);
        Amet += interpolated_A(0, k) * interpolated_A(0, k);
      }

      double gamma = 0.5 * (Amet - amet);

      // Calculate the contravariant metric tensors
      double adet = amet; // double aup = 1.0/amet;
      double Adet = Amet; // double Aup = 1.0/Amet;

      // Calculate the unit normal Vectors
      Vector<double> n(2);
      Vector<double> N(2);
      n[0] = -a(0, 1) / sqrt(adet);
      n[1] = a(0, 0) / sqrt(adet);

      N[0] = -interpolated_A(0, 1) / sqrt(Adet);
      N[1] = interpolated_A(0, 0) / sqrt(Adet);

      // Square root of deformed determinant
      double sqrt_Adet = sqrt(Adet);

      // Calculate the curvature tensors
      double b = n[0] * dadxi(0, 0, 0) + n[1] * dadxi(0, 0, 1);
      double B =
        N[0] * interpolated_dAdxi(0, 0) + N[1] * interpolated_dAdxi(0, 1);

      // Set up the change of curvature tensor
      double kappa = b - B;

      // Define the load Vector
      Vector<double> f(n_dim);

      // Get load Vector
      load_vector(ipt, interpolated_xi, interpolated_x, N, f);

      // Define the wall thickness ratio profile
      double h_ratio = 0;

      // Get wall thickness ratio profile
      wall_profile(interpolated_xi, interpolated_x, h_ratio);

      // Thickness h/R
      double HoR = HoR_0 * h_ratio;

      // Allocate storage for variations of normal Vector
      double normal_var[n_dim][n_dim];

      // Geometrically nonlinear beam equations from
      //--------------------------------------------
      // principle of virtual displacements
      //-----------------------------------

      // Loop over the number of nodes
      for (unsigned n = 0; n < n_node; n++)
      {
        // Loop over the type of degree of freedom
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Setup components of variation of normal Vector:
          // ...(variation w.r.t. direction, component)

          normal_var[0][0] = interpolated_A(0, 0) * interpolated_A(0, 1) *
                             dpsidxi(n, k, 0) /
                             (sqrt_Adet * sqrt_Adet * sqrt_Adet);

          normal_var[1][0] =
            (interpolated_A(0, 1) * interpolated_A(0, 1) / Adet - 1.0) *
            dpsidxi(n, k, 0) / (sqrt_Adet);

          normal_var[0][1] =
            (1.0 - interpolated_A(0, 0) * interpolated_A(0, 0) / Adet) *
            dpsidxi(n, k, 0) / (sqrt_Adet);

          normal_var[1][1] = -interpolated_A(0, 0) * interpolated_A(0, 1) *
                             dpsidxi(n, k, 0) /
                             (sqrt_Adet * sqrt_Adet * sqrt_Adet);

          // Loop over the coordinate directions
          for (unsigned i = 0; i < n_dim; i++)
          {
            // Find the equation number
            local_eqn = position_local_eqn(n, k, i);

            // If it's not a boundary condition
            if (local_eqn >= 0)
            {
              // Premultiply by thickness profile

              // External forcing
              residuals[local_eqn] -=
                h_ratio * (1.0 / HoR) * f[i] * psi(n, k) * W * sqrt(Adet);

              residuals[local_eqn] +=
                h_ratio * Lambda_sq * accel[i] * psi(n, k) * W * sqrt(adet);

              // Membrane term with axial prestress
              residuals[local_eqn] += h_ratio * (sigma_0 + gamma) *
                                      interpolated_A(0, i) * dpsidxi(n, k, 0) *
                                      W * sqrt(adet);

              // Bending term: Minus sign because \delta \kappa = - \delta B
              residuals[local_eqn] -=
                h_ratio * (1.0 / 12.0) * HoR * HoR * kappa *
                (N[i] * d2psidxi(n, k, 0) +
                 normal_var[i][0] * interpolated_dAdxi(0, 0) +
                 normal_var[i][1] * interpolated_dAdxi(0, 1)) *
                W * sqrt(adet);
            }
          }
        }
      }

    } // End of loop over the integration points
  }


  //=======================================================================
  /// Get FE jacobian and residuals (Jacobian done by finite differences)
  //=======================================================================
  void KirchhoffLoveBeamEquations::fill_in_contribution_to_jacobian(
    Vector<double>& residuals, DenseMatrix<double>& jacobian)
  {
    // Call the element's residuals vector
    fill_in_contribution_to_residuals_beam(residuals);

    // Solve for the consistent acceleration in Newmark scheme?
    if (Solve_for_consistent_newmark_accel_flag)
    {
      SolidFiniteElement::fill_in_jacobian_for_newmark_accel(jacobian);
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
  /// element (wrapper function).
  //=======================================================================
  void KirchhoffLoveBeamEquations::get_energy(double& pot_en, double& kin_en)
  {
    // Initialise
    double stretch = 0.0;
    double bend = 0.0;
    kin_en = 0.0;

    // Compute energy
    get_energy(stretch, bend, kin_en);

    // Sum components of potential energy
    pot_en = stretch + bend;
  }


  //=======================================================================
  /// Compute the potential (strain) and kinetic energy of the
  /// element, breaking down the potential energy into stretching and
  /// bending components.
  //=======================================================================
  void KirchhoffLoveBeamEquations::get_energy(double& stretch,
                                              double& bend,
                                              double& kin_en)
  {
    // Initialise
    stretch = 0.0;
    bend = 0.0;
    kin_en = 0.0;

    // Set the dimension of the global coordinates
    unsigned n_dim = Undeformed_beam_pt->ndim();

    // Set the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_dofs = nnodal_position_type();

    // Setup memory for veloc
    Vector<double> veloc(2);

    // Set up memory for the shape functions:

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_dofs);

    // # of nodes, # of positional dofs, # of lagrangian coords (for deriv)
    DShape dpsidxi(n_node, n_position_dofs, n_lagrangian);

    // # of nodes, # of positional dofs, # of derivs)
    DShape d2psidxi(n_node, n_position_dofs, n_lagrangian);

    // Set # of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Get Physical Variables from Element

    // Thickness h_0/R, axial prestress, timescale ratio (density)
    double HoR_0 = h(); // i.e. refers to reference thickness 'h_0'
    double sigma_0 = sigma0();

    double Lambda_sq = lambda_sq();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape functions w.r.t. Lagrangian coords
      double J = d2shape_lagrangian_at_knot(ipt, psi, dpsidxi, d2psidxi);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of lagrangian position and
      // the derivative of global position wrt lagrangian coordinates
      Vector<double> interpolated_xi(n_lagrangian, 0.0),
        interpolated_x(n_dim, 0.0);

      DenseMatrix<double> interpolated_A(n_lagrangian, n_dim);
      DenseMatrix<double> interpolated_dAdxi(n_lagrangian, n_dim);

      // Initialise to zero

      // Loop over coordinate directions/components of Vector
      for (unsigned i = 0; i < n_dim; i++)
      {
        // Initialise veloc
        veloc[i] = 0.0;

        // Loop over derivatives/base Vectors (just one here)
        for (unsigned j = 0; j < n_lagrangian; j++)
        {
          interpolated_A(j, i) = 0.0;
        }
        // Loop over derivatives of base Vector (just one here)
        for (unsigned j = 0; j < n_lagrangian; j++)
        {
          interpolated_dAdxi(j, i) = 0.0;
        }
      }

      // Calculate displacements, accelerations and spatial derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_dofs; k++)
        {
          // Loop over Lagrangian coordinate directions [xi_gen[] are the
          // the *gen*eralised Lagrangian coordinates: node, type, direction]
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            interpolated_xi[i] +=
              raw_lagrangian_position_gen(l, k, i) * psi(l, k);
          }

          // Loop over components of the deformed position Vector
          for (unsigned i = 0; i < n_dim; i++)
          {
            // Need this for wall thickness ratio profile
            interpolated_x[i] += raw_nodal_position_gen(l, k, i) * psi(l, k);

            veloc[i] += raw_dnodal_position_gen_dt(1, l, k, i) * psi(l, k);

            // Loop over derivative directions (just one here)
            for (unsigned j = 0; j < n_lagrangian; j++)
            {
              interpolated_A(j, i) +=
                raw_nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
            }

            // Loop over the second derivative directions (just one here)
            for (unsigned j = 0; j < n_lagrangian; j++)
            {
              interpolated_dAdxi(j, i) +=
                raw_nodal_position_gen(l, k, i) * d2psidxi(l, k, j);
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

      // Setup position Vector and derivatives of undeformed config
      Vector<double> R(n_dim);
      DenseMatrix<double> a(n_lagrangian, n_dim);
      RankThreeTensor<double> dadxi(n_lagrangian, n_lagrangian, n_dim);

      // Get the undeformed geometry
      Undeformed_beam_pt->d2position(interpolated_xi, R, a, dadxi);

      // Declare and calculate the undeformed and deformed metric tensor
      // and the strain tensor (these are 1d tensors, i.e. scalars)
      double amet = 0.0, Amet = 0.0;

      // Work out metric and strain tensors
      // Now calculate the dot product
      for (unsigned k = 0; k < n_dim; k++)
      {
        amet += a(0, k) * a(0, k);
        Amet += interpolated_A(0, k) * interpolated_A(0, k);
      }

      // Calculate strain tensor
      double gamma = 0.5 * (Amet - amet);

      // Calculate the contravariant metric tensors
      double adet = amet; // aup = 1.0/amet;
      double Adet = Amet; // Aup = 1.0/Amet;

      // Calculate the unit normal Vectors
      Vector<double> n(2);
      Vector<double> N(2);
      n[0] = -a(0, 1) / sqrt(adet);
      n[1] = a(0, 0) / sqrt(adet);

      N[0] = -interpolated_A(0, 1) / sqrt(Adet);
      N[1] = interpolated_A(0, 0) / sqrt(Adet);

      // Calculate the curvature tensors
      double b = n[0] * dadxi(0, 0, 0) + n[1] * dadxi(0, 0, 1);
      double B =
        N[0] * interpolated_dAdxi(0, 0) + N[1] * interpolated_dAdxi(0, 1);

      // Set up the change of curvature tensor
      double kappa = b - B;

      // Define the wall thickness ratio profile
      double h_ratio = 0;

      // Get wall thickness ratio profile
      wall_profile(interpolated_xi, interpolated_x, h_ratio);

      // Thickness h/R
      double HoR = HoR_0 * h_ratio;

      // Add contributions
      stretch +=
        h_ratio * 0.5 * (gamma + sigma_0) * (gamma + sigma_0) * W * sqrt(adet);
      bend += h_ratio * 0.5 * (1.0 / 12.0) * HoR * HoR * kappa * kappa * W *
              sqrt(adet);
      kin_en += h_ratio * 0.5 * Lambda_sq * veloc_sq * W * sqrt(adet);
    } // End of loop over the integration points
  }

  //=======================================================================
  /// Output position at previous time (t=0: present; t>0: previous)
  /// with specified number of plot points.
  //=======================================================================
  void HermiteBeamElement::output(const unsigned& t,
                                  std::ostream& outfile,
                                  const unsigned& n_plot) const
  {
#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "HermiteBeamElement::output(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << std::endl;

    // Find the dimension of the first node
    unsigned n_dim = this->node_pt(0)->ndim();

    // Loop over plot points
    for (unsigned l = 0; l < n_plot; l++)
    {
      s[0] = -1.0 + l * 2.0 / (n_plot - 1);

      // Output the Eulerian coordinates
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << interpolated_x(t, s, i) << " ";
      }
      outfile << std::endl;
    }
  }


  //=======================================================================
  /// Output function
  //=======================================================================
  void HermiteBeamElement::output(std::ostream& outfile, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << std::endl;

    // Set the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // Set the dimension of the global coordinates
    unsigned n_dim = Undeformed_beam_pt->ndim();

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_dofs = nnodal_position_type();

    Vector<double> veloc(n_dim);
    Vector<double> accel(n_dim);
    Vector<double> posn(n_dim);

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_dofs);

    // Loop over element plot points
    for (unsigned l1 = 0; l1 < n_plot; l1++)
    {
      s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

      // Get shape functions
      shape(s, psi);

      Vector<double> interpolated_xi(n_lagrangian);
      interpolated_xi[0] = 0.0;

      // Loop over coordinate directions/components of Vector
      for (unsigned i = 0; i < n_dim; i++)
      {
        // Initialise acclerations and veloc
        accel[i] = 0.0;
        veloc[i] = 0.0;
        posn[i] = 0.0;
      }


      // Calculate displacements, accelerations and spatial derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_dofs; k++)
        {
          // Loop over Lagrangian coordinate directions [xi_gen[] are the
          // the *gen*eralised Lagrangian coordinates: node, type, direction]
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            interpolated_xi[i] +=
              raw_lagrangian_position_gen(l, k, i) * psi(l, k);
          }

          // Loop over components of the deformed position Vector
          for (unsigned i = 0; i < n_dim; i++)
          {
            accel[i] += raw_dnodal_position_gen_dt(2, l, k, i) * psi(l, k);
            veloc[i] += raw_dnodal_position_gen_dt(1, l, k, i) * psi(l, k);
            posn[i] += raw_dnodal_position_gen_dt(0, l, k, i) * psi(l, k);
          }
        }
      }

      double scalar_accel = 0.0;
      double scalar_veloc = 0.0;
      double scalar_posn = 0.0;

      // Output position etc.
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << posn[i] << " ";
        scalar_posn += pow(posn[i], 2);
      }
      outfile << interpolated_xi[0] << " ";
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << veloc[i] << " ";
        scalar_veloc += pow(veloc[i], 2);
      }
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << accel[i] << " ";
        scalar_accel += pow(accel[i], 2);
      }
      outfile << sqrt(scalar_posn) << " ";
      outfile << sqrt(scalar_veloc) << " ";
      outfile << sqrt(scalar_accel) << " ";
      outfile << std::endl;
    }
  }


  //=======================================================================
  /// Output function
  //=======================================================================
  void HermiteBeamElement::output(std::ostream& outfile)
  {
    unsigned n_plot = 5;
    output(outfile, n_plot);
  }


  //=======================================================================
  /// C-style output position at previous time (t=0: present; t>0: previous)
  /// with specified number of plot points.
  //=======================================================================
  void HermiteBeamElement::output(const unsigned& t,
                                  FILE* file_pt,
                                  const unsigned& n_plot) const
  {
#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "HermiteBeamElement::output(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    // outfile << "ZONE I=" << n_plot << std::endl;
    fprintf(file_pt, "ZONE I=%i \n", n_plot);

    // Find the dimension of the first node
    unsigned n_dim = this->node_pt(0)->ndim();

    // Loop over plot points
    for (unsigned l = 0; l < n_plot; l++)
    {
      s[0] = -1.0 + l * 2.0 / (n_plot - 1);

      // Output the Eulerian coordinates
      for (unsigned i = 0; i < n_dim; i++)
      {
        // outfile << interpolated_x(s,t,i) << " " ;
        fprintf(file_pt, "%g ", interpolated_x(t, s, i));
      }
      // outfile << std::endl;
      fprintf(file_pt, "\n");
    }
  }


  //=======================================================================
  /// Output function
  //=======================================================================
  void HermiteBeamElement::output(FILE* file_pt, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    // outfile << "ZONE I=" << n_plot << std::endl;
    fprintf(file_pt, "ZONE I=%i \n", n_plot);

    // Set the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // Set the dimension of the global coordinates
    unsigned n_dim = Undeformed_beam_pt->ndim();

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_dofs = nnodal_position_type();

    Vector<double> veloc(n_dim);
    Vector<double> accel(n_dim);
    Vector<double> posn(n_dim);

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_dofs);

    // Loop over element plot points
    for (unsigned l1 = 0; l1 < n_plot; l1++)
    {
      s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

      // Get shape functions
      shape(s, psi);

      Vector<double> interpolated_xi(n_lagrangian);
      interpolated_xi[0] = 0.0;

      // Loop over coordinate directions/components of Vector
      for (unsigned i = 0; i < n_dim; i++)
      {
        // Initialise acclerations and veloc
        accel[i] = 0.0;
        veloc[i] = 0.0;
        posn[i] = 0.0;
      }


      // Calculate displacements, accelerations and spatial derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_dofs; k++)
        {
          // Loop over Lagrangian coordinate directions [xi_gen[] are the
          // the *gen*eralised Lagrangian coordinates: node, type, direction]
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            interpolated_xi[i] +=
              raw_lagrangian_position_gen(l, k, i) * psi(l, k);
          }

          // Loop over components of the deformed position Vector
          for (unsigned i = 0; i < n_dim; i++)
          {
            accel[i] += raw_dnodal_position_gen_dt(2, l, k, i) * psi(l, k);
            veloc[i] += raw_dnodal_position_gen_dt(1, l, k, i) * psi(l, k);
            posn[i] += raw_dnodal_position_gen_dt(0, l, k, i) * psi(l, k);
          }
        }
      }

      double scalar_accel = 0.0;
      double scalar_veloc = 0.0;
      double scalar_posn = 0.0;

      // Output position etc.
      for (unsigned i = 0; i < n_dim; i++)
      {
        // outfile << posn[i] << " " ;
        fprintf(file_pt, "%g ", posn[i]);
        scalar_posn += pow(posn[i], 2);
      }
      // outfile << interpolated_xi[0] << " ";
      fprintf(file_pt, "%g ", interpolated_xi[0]);
      for (unsigned i = 0; i < n_dim; i++)
      {
        // outfile << veloc[i] << " " ;
        fprintf(file_pt, "%g ", veloc[i]);
        scalar_veloc += pow(veloc[i], 2);
      }
      for (unsigned i = 0; i < n_dim; i++)
      {
        // outfile << accel[i] << " " ;
        fprintf(file_pt, "%g ", accel[i]);
        scalar_accel += pow(accel[i], 2);
      }
      // outfile << sqrt(scalar_posn)  << " ";
      fprintf(file_pt, "%g ", sqrt(scalar_posn));
      // outfile << sqrt(scalar_veloc) << " ";
      fprintf(file_pt, "%g ", sqrt(scalar_veloc));
      // outfile << sqrt(scalar_accel) << " ";
      fprintf(file_pt, "%g ", sqrt(scalar_accel));
      // outfile << std::endl;
      fprintf(file_pt, "\n");
    }
  }


  //=======================================================================
  /// C-style output function
  //=======================================================================
  void HermiteBeamElement::output(FILE* file_pt)
  {
    unsigned n_plot = 5;
    output(file_pt, n_plot);
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Function used to find the local coordinate s that corresponds to the
  /// "global" intrinsic coordinate zeta (in the element's incarnation
  /// as a GeomObject). For this element, zeta is equal to the Lagrangian
  /// coordinate xi. If the required zeta is located within this
  /// element, geom_object_pt points to "this" element. If zeta
  /// is not located within the element it is set to NULL.
  //======================================================================
  void FSIHermiteBeamElement::locate_zeta(
    const Vector<double>& zeta,
    GeomObject*& geom_object_pt,
    Vector<double>& s,
    const bool& use_coordinate_as_initial_guess)
  {
    // Assumed that the first node has a lower xi coordinate than the second
    unsigned lo = 0, hi = 1;

    // If the first node has a higher xi then swap
    if (raw_lagrangian_position(0, 0) > raw_lagrangian_position(1, 0))
    {
      lo = 1;
      hi = 0;
    }

    // Tolerance for finding zeta
    double epsilon = 1.0e-13;

    // If zeta is not in the element, then return a null pointer
    // Correct for rounding errors here
    if ((zeta[0] - raw_lagrangian_position(lo, 0) < -epsilon) ||
        (zeta[0] - raw_lagrangian_position(hi, 0) > epsilon))
    {
      geom_object_pt = 0;
      return;
    }

    // Otherwise, zeta is located in this element. For now assume
    // that the relationship between zeta and s is linear as it will
    // be for uniform node spacing... We'll trap this further down.
    // In general we'll need a Newton method here and we'll implement
    // this as soon as we have an example with non-uniformly spaced
    // FSIHermiteBeamElements...

    // The GeomObject that contains the zeta coordinate is "this":
    geom_object_pt = this;


    // Find the fraction along the element
    double zeta_fraction =
      (zeta[0] - raw_lagrangian_position(lo, 0)) /
      (raw_lagrangian_position(hi, 0) - raw_lagrangian_position(lo, 0));

    s[0] = -1.0 + zeta_fraction * 2.0;

#ifdef PARANOID
    // Check if we've really ended where we wanted to be...
    Vector<double> zeta_test(1);
    interpolated_zeta(s, zeta_test);
    if (std::fabs(zeta[0] - zeta_test[0]) > epsilon)
    {
      std::ostringstream error_stream;
      error_stream
        << "The zeta coordinate " << zeta_test[0] << " \n"
        << "computed by interpolated_zeta() for s[0]=" << s[0] << " \n"
        << "differs by more than the tolerance (" << epsilon << ") from \n "
        << "the required value " << zeta_test[0] << " \n\n"
        << "You're probably using a mesh with non-uniformly \n "
        << "spaced  FSIHermiteBeamElements. For such cases the root finding"
        << "in FSIHermiteBeamElement::locate_zeta() must be replaced "
        << "by a proper Newton method or some such thing...\n";
      OomphLibError(error_stream.str(),
                    "FSIHermiteBeamElement::locate_zeta()",
                    OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Check for rounding error
    if (s[0] > 1.0)
    {
      s[0] = 1.0;
    }
    if (s[0] < -1.0)
    {
      s[0] = -1.0;
    }
  }


  //=========================================================================
  /// Define the dposition function. This is used to set no-slip boundary
  /// conditions in FSI problems.
  //=========================================================================
  void FSIHermiteBeamElement::dposition_dlagrangian_at_local_coordinate(
    const Vector<double>& s, DenseMatrix<double>& drdxi) const
  {
#ifdef PARANOID
    if (Undeformed_beam_pt == 0)
    {
      throw OomphLibError(
        "Undeformed_beam_pt has not been set",
        "FSIHermiteBeamElement::dposition_dlagrangian_at_local_coordinate()",
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the dimension of the global coordinate
    unsigned n_dim = Undeformed_beam_pt->ndim();
    // Find the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();
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
  /// of the "DOF type" that this unknown is associated with.
  /// (Function can obviously only be called if the equation numbering
  /// scheme has been set up.)
  /// This element is only in charge of the solid dofs.
  //=============================================================================
  void FSIHermiteBeamElement::get_dof_numbers_for_unknowns(
    std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // temporary pair (used to store dof lookup prior to being added to list)
    std::pair<unsigned, unsigned> dof_lookup;

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


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////


  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  ClampedSlidingHermiteBeamBoundaryConditionElement::
    ClampedSlidingHermiteBeamBoundaryConditionElement(
      FiniteElement* const& bulk_el_pt, const int& face_index)
    : FaceGeometry<HermiteBeamElement>(), FaceElement()
  {
    // Number of nodal position types: 2
    set_nnodal_position_type(2);

    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index, this);

    // Resize vector to some point on the symmetry line along which the
    // end of the beam is sliding. Initialise for symmetry line = y-axis
    Vector_to_symmetry_line.resize(2);
    Vector_to_symmetry_line[0] = 0.0;
    Vector_to_symmetry_line[1] = 0.0;

    // Resize normal vector to the symmetry line along which the
    // end of the beam is sliding. Initialise for symmetry line = y-axis
    Normal_to_symmetry_line.resize(2);
    Normal_to_symmetry_line[0] = -1.0;
    Normal_to_symmetry_line[1] = 0.0;

    // Resize number of dofs at the element's one and only node by two
    // to accomodate the Lagrange multipliers

    // We need two additional values at the one-and-only node in this element
    Vector<unsigned> nadditional_data_values(1);
    nadditional_data_values[0] = 2;
    resize_nodes(nadditional_data_values);
  }


  //===========================================================================
  /// Add the element's contribution to its residual vector
  //===========================================================================
  void ClampedSlidingHermiteBeamBoundaryConditionElement::
    fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
    // Get position vector to and normal & tangent vectors on wall (from
    // bulk element)
    HermiteBeamElement* bulk_el_pt =
      dynamic_cast<HermiteBeamElement*>(bulk_element_pt());

    // Get the value of the constant local coordinate in the bulk element
    // Dummy local coordinate of size zero
    Vector<double> s(0);
    Vector<double> s_bulk(1);
    this->get_local_coordinate_in_bulk(s, s_bulk);

    // Get position vector to and normal on wall
    Vector<double> r(2);
    Vector<double> n(2);
    bulk_el_pt->get_normal(s_bulk, r, n);

    // Get (non-unit) tangent vector
    Vector<double> drds(2);
    bulk_el_pt->get_non_unit_tangent(s_bulk, r, drds);

    // Residual corresponding to: Point must be on symmetry line
    double res0 =
      (r[0] - Vector_to_symmetry_line[0]) * Normal_to_symmetry_line[0] +
      (r[1] - Vector_to_symmetry_line[1]) * Normal_to_symmetry_line[1];

    // Work out tangent to symmetry line
    Vector<double> tangent_to_symmetry_line(2);
    tangent_to_symmetry_line[0] = -Normal_to_symmetry_line[1];
    tangent_to_symmetry_line[1] = Normal_to_symmetry_line[0];

    // Residual corresponding to: Beam must meet symmetry line at right angle
    double res1 = drds[0] * tangent_to_symmetry_line[0] +
                  drds[1] * tangent_to_symmetry_line[1];

    // The first Lagrange multiplier enforces the along-the-beam displacement:
    int j_local = nodal_local_eqn(0, Nbulk_value[0]);
    residuals[j_local] += res0;

    // Second Lagrange multiplier enforces the derivative of the
    // transverse displacement:
    j_local = nodal_local_eqn(0, Nbulk_value[0] + 1);
    residuals[j_local] += res1;

    // Add Lagrange multiplier contribution to the bulk equations

    // Lagrange multipliers
    double lambda0 = node_pt(0)->value(Nbulk_value[0]);
    double lambda1 = node_pt(0)->value(Nbulk_value[0] + 1);

    // Loop over the solid dofs

    // How many positional values are there?
    unsigned n_dim = nodal_dimension();
    unsigned n_type = nnodal_position_type();

    /// Loop over directions
    for (unsigned i = 0; i < n_dim; i++)
    {
      // Loop over types
      for (unsigned k = 0; k < n_type; k++)
      {
        // Get local equation number
        int j_local = position_local_eqn(0, k, i);

        // Real dof?
        if (j_local >= 0)
        {
          // Derivative of first residual w.r.t. to discrete displacements:
          // Only the type-zero dofs make a contribution:
          double dres0dxik = 0.0;
          if (k == 0)
          {
            dres0dxik = Normal_to_symmetry_line[i];
          }

          // Derivative of second residual w.r.t. to discrete displacements:
          // Only the type-one dofs make a contribution:
          double dres1dxik = 0.0;
          if (k == 1)
          {
            dres1dxik = tangent_to_symmetry_line[i];
          }

          // Add to the residual
          residuals[j_local] += lambda0 * dres0dxik + lambda1 * dres1dxik;
        }
      }
    }
  }


} // namespace oomph
