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
// Non-inline functions for elements that solve the equations of linear
// elasticity in cartesian coordinates

#include "axisym_linear_elasticity_elements.h"


namespace oomph
{
  /// \short Static default value for Young's modulus (1.0 -- for natural
  /// scaling, i.e. all stresses have been non-dimensionalised by
  /// the same reference Young's modulus. Setting the "non-dimensional"
  /// Young's modulus (obtained by de-referencing Youngs_modulus_pt)
  /// to a number larger than one means that the material is stiffer
  /// than assumed in the non-dimensionalisation.
  double
    AxisymmetricLinearElasticityEquationsBase::Default_youngs_modulus_value =
      1.0;

  /// Static default value for timescale ratio (1.0 -- for natural scaling)
  double AxisymmetricLinearElasticityEquationsBase::Default_lambda_sq_value =
    1.0;

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Get strain (3x3 entries; r, z, phi)
  //=======================================================================
  void AxisymmetricLinearElasticityEquations::get_strain(
    const Vector<double>& s, DenseMatrix<double>& strain)
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();


#ifdef PARANOID
    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();

    if (n_position_type != 1)
    {
      throw OomphLibError("AxisymmetricLinearElasticity is not yet implemented "
                          "for more than one position type",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the indices at which the local displacements are stored
    unsigned u_nodal_index[3];
    for (unsigned i = 0; i < 3; i++)
    {
      u_nodal_index[i] = this->u_index_axisymmetric_linear_elasticity(i);
    }

    // Set up memory for the shape functions
    Shape psi(n_node);

    // Derivs only w.r.t. r [0] and z [1]
    DShape dpsidx(n_node, 2);

    // Storage for Eulerian coordinates (r,z; initialised to zero)
    Vector<double> interpolated_x(2, 0.0);

    // Displacements u_r,u_z,u_theta
    Vector<double> interpolated_u(3, 0.0);

    // Calculate interpolated values of the derivatives w.r.t.
    // Eulerian coordinates
    DenseMatrix<double> interpolated_dudx(3, 2, 0.0);


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
        const double u_value = this->raw_nodal_value(l, u_nodal_index[i]);
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

    // r component of displacement
    double ur = interpolated_u[0];

    // theta component of displacement
    double uth = interpolated_u[2];

    // derivatives w.r.t. r and z:
    double durdr = interpolated_dudx(0, 0);
    double durdz = interpolated_dudx(0, 1);
    double duzdr = interpolated_dudx(1, 0);
    double duzdz = interpolated_dudx(1, 1);
    double duthdr = interpolated_dudx(2, 0);
    double duthdz = interpolated_dudx(2, 1);


    // e_rr
    strain(0, 0) = durdr;
    // e_rz
    strain(0, 1) = 0.5 * (durdz + duzdr);
    strain(1, 0) = 0.5 * (durdz + duzdr);
    // e_rphi
    strain(0, 2) = 0.5 * (duthdr - uth / r);
    strain(2, 0) = 0.5 * (duthdr - uth / r);
    // e_zz
    strain(1, 1) = duzdz;
    // e_zphi
    strain(1, 2) = 0.5 * duthdz;
    strain(2, 1) = 0.5 * duthdz;
    // e_phiphi
    strain(2, 2) = ur / r;
  }

  //=======================================================================
  /// Compute the residuals for the axisymmetric (in cyl. polars)
  /// linear elasticity equations in. Flag indicates if we want Jacobian too.
  //=======================================================================
  void AxisymmetricLinearElasticityEquations::
    fill_in_generic_contribution_to_residuals_axisymmetric_linear_elasticity(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

#ifdef PARANOID
    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();

    if (n_position_type != 1)
    {
      throw OomphLibError("AxisymmetricLinearElasticity is not yet implemented "
                          "for more than one position type",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the indices at which the local displacements are stored
    unsigned u_nodal_index[3];
    for (unsigned i = 0; i < 3; i++)
    {
      u_nodal_index[i] = this->u_index_axisymmetric_linear_elasticity(i);
    }

    // Get elastic parameters
    double nu_local = this->nu();
    double youngs_modulus_local = this->youngs_modulus();

    // Obtain Lame parameters from Young's modulus and Poisson's ratio
    double lambda = youngs_modulus_local * nu_local / (1.0 + nu_local) /
                    (1.0 - 2.0 * nu_local);
    double mu = youngs_modulus_local / 2.0 / (1.0 + nu_local);


    // Lambda squared --- time scaling, NOT sqaure of Lame parameter lambda
    const double lambda_sq = this->lambda_sq();

    // Set up memory for the shape functions
    Shape psi(n_node);

    // Derivs only w.r.t. r [0] and z [1]
    DShape dpsidx(n_node, 2);

    // Set the value of Nintpt -- the number of integration points
    unsigned n_intpt = this->integral_pt()->nweight();

    // Set the vector to hold the local coordinates in the element
    Vector<double> s(2);

    // Integers to store the local equation numbers
    int local_eqn = 0, local_unknown = 0;

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
      Vector<double> interpolated_u(3, 0.0);

      // Calculate interpolated values of the derivatives w.r.t.
      // Eulerian coordinates
      DenseMatrix<double> interpolated_dudx(3, 2, 0.0);

      Vector<double> d2u_dt2(3, 0.0);

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
          const double u_value = this->raw_nodal_value(l, u_nodal_index[i]);

          interpolated_u[i] += u_value * psi(l);

          d2u_dt2[i] += d2u_dt2_axisymmetric_linear_elasticity(l, i) * psi(l);

          // Loop over derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsidx(l, j);
          }
        }
      }

      // Get body force
      Vector<double> b(3);
      this->body_force(time, interpolated_x, b);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      //=====EQUATIONS OF AXISYMMETRIC LINEAR ELASTICITY ========

      // define shorthand notation for regularly-occurring terms
      double r = interpolated_x[0];
      double rsq = pow(r, 2);

      // r component of displacement
      double ur = interpolated_u[0];

      // theta component of displacement
      double uth = interpolated_u[2];

      // derivatives w.r.t. r and z:
      double durdr = interpolated_dudx(0, 0);
      double durdz = interpolated_dudx(0, 1);
      double duzdr = interpolated_dudx(1, 0);
      double duzdz = interpolated_dudx(1, 1);
      double duthdr = interpolated_dudx(2, 0);
      double duthdz = interpolated_dudx(2, 1);

      // storage for terms required for analytic Jacobian
      double G_r, G_z, G_theta;

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the displacement components
        for (unsigned a = 0; a < 3; a++)
        {
          // Get the equation number
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[a]);

          /*IF it's not a boundary condition*/
          if (local_eqn >= 0)
          {
            // Acceleration and body force
            residuals[local_eqn] +=
              (lambda_sq * d2u_dt2[a] - b[a]) * psi(l) * r * W;

            // Three components of the stress divergence term:
            // a=0: r; a=1: z; a=2: theta

            // r-equation
            if (a == 0)
            {
              residuals[local_eqn] += (mu * (2.0 * durdr * dpsidx(l, 0) +
                                             +dpsidx(l, 1) * (durdz + duzdr) +
                                             2.0 * psi(l) / pow(r, 2) * (ur)) +
                                       lambda * (durdr + ur / r + duzdz) *
                                         (dpsidx(l, 0) + psi(l) / r)) *
                                      r * W;
            }
            // z-equation
            else if (a == 1)
            {
              residuals[local_eqn] +=
                (mu * (dpsidx(l, 0) * (durdz + duzdr) +
                       2.0 * duzdz * dpsidx(l, 1)) +
                 lambda * (durdr + ur / r + duzdz) * dpsidx(l, 1)) *
                r * W;
            }
            // theta-equation
            else if (a == 2)
            {
              residuals[local_eqn] +=
                (mu * ((duthdr - uth / r) * (dpsidx(l, 0) - psi(l) / r) +
                       dpsidx(l, 1) * (duthdz))) *
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
                // define terms used to obtain entries of current row in the
                // Jacobian:

                // terms for rows of Jacobian matrix corresponding to r-equation
                if (a == 0)
                {
                  G_r =
                    (mu * (2.0 * dpsidx(l2, 0) * dpsidx(l, 0) +
                           2.0 / rsq * psi(l2) * psi(l) +
                           dpsidx(l2, 1) * dpsidx(l, 1)) +
                     lambda * (dpsidx(l2, 0) + psi(l2) / r) *
                       (dpsidx(l, 0) + psi(l) / r) +
                     lambda_sq * node_pt(l2)->time_stepper_pt()->weight(2, 0) *
                       psi(l2) * psi(l)) *
                    r * W;

                  G_z = (mu * dpsidx(l2, 0) * dpsidx(l, 1) +
                         lambda * dpsidx(l2, 1) * (dpsidx(l, 0) + psi(l) / r)) *
                        r * W;

                  G_theta = 0;
                }

                // terms for rows of Jacobian matrix corresponding to z-equation
                else if (a == 1)
                {
                  G_r =
                    (mu * dpsidx(l2, 1) * dpsidx(l, 0) +
                     lambda * (dpsidx(l2, 0) + psi(l2) / r) * dpsidx(l, 1)) *
                    r * W;

                  G_z =
                    (mu * (dpsidx(l2, 0) * dpsidx(l, 0) +
                           2.0 * dpsidx(l2, 1) * dpsidx(l, 1)) +
                     lambda * dpsidx(l2, 1) * dpsidx(l, 1) +
                     lambda_sq * node_pt(l2)->time_stepper_pt()->weight(2, 0) *
                       psi(l2) * psi(l)) *
                    r * W;

                  G_theta = 0;
                }

                // terms for rows of Jacobian matrix corresponding to
                // theta-equation
                else if (a == 2)
                {
                  G_r = 0;

                  G_z = 0;

                  G_theta =
                    (mu * ((dpsidx(l2, 0) - psi(l2) / r) *
                             (dpsidx(l, 0) - psi(l) / r) +
                           dpsidx(l2, 1) * dpsidx(l, 1)) +
                     lambda_sq * node_pt(l2)->time_stepper_pt()->weight(2, 0) *
                       psi(l2) * psi(l)) *
                    r * W;
                }

                // Loop over the displacement components
                for (unsigned c = 0; c < 3; c++)
                {
                  // Get local unknown
                  local_unknown = this->nodal_local_eqn(l2, u_nodal_index[c]);

                  // If the local unknown is not pinned
                  if (local_unknown >= 0)
                  {
                    if (c == 0)
                    {
                      jacobian(local_eqn, local_unknown) += G_r;
                    }
                    else if (c == 1)
                    {
                      jacobian(local_eqn, local_unknown) += G_z;
                    }
                    else if (c == 2)
                    {
                      jacobian(local_eqn, local_unknown) += G_theta;
                    }
                  }
                }
              }
            } // End of jacobian calculation

          } // End of if not boundary condition
        } // End of loop over coordinate directions
      } // End of loop over shape functions
    } // End of loop over integration points
  }

  //=======================================================================
  /// Output exact solution  r,z, u_r, u_z, u_theta
  //=======================================================================
  void AxisymmetricLinearElasticityEquations::output_fct(
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
    Vector<double> exact_soln(9);

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

      for (unsigned i = 0; i < 3; i++)
      {
        outfile << exact_soln[i] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  //=======================================================================
  /// Output exact solution  r,z, u_r, u_z, u_theta
  /// Time dependent version
  //=======================================================================
  void AxisymmetricLinearElasticityEquations::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordintes
    Vector<double> x(2);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Exact solution Vector
    Vector<double> exact_soln(9);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(time, x, exact_soln);

      // Output r,z,...,u_exact,...
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      for (unsigned i = 0; i < 9; i++)
      {
        outfile << exact_soln[i] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  //=======================================================================
  /// Output: r,z, u_r, u_z, u_theta
  //=======================================================================
  void AxisymmetricLinearElasticityEquations::output(std::ostream& outfile,
                                                     const unsigned& nplot)
  {
    // Set output Vector
    Vector<double> s(2);
    Vector<double> x(2);
    Vector<double> u(3);
    Vector<double> du_dt(3);
    Vector<double> d2u_dt2(3);


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
      this->interpolated_u_axisymmetric_linear_elasticity(s, u);
      this->interpolated_du_dt_axisymmetric_linear_elasticity(s, du_dt);
      this->interpolated_d2u_dt2_axisymmetric_linear_elasticity(s, d2u_dt2);

      // Output the r,z,..
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output displacement
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << u[i] << " ";
      }

      // Output veloc
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << du_dt[i] << " ";
      }

      // Output accel
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << d2u_dt2[i] << " ";
      }

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }


  //=======================================================================
  /// C-style output:r,z, u_r, u_z, u_theta
  //=======================================================================
  void AxisymmetricLinearElasticityEquations::output(FILE* file_pt,
                                                     const unsigned& nplot)
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
        fprintf(file_pt,
                "%g ",
                this->interpolated_u_axisymmetric_linear_elasticity(s, i));
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
  void AxisymmetricLinearElasticityEquations::compute_error(
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

    // Exact solution Vector (u_r, u_z, u_theta)
    Vector<double> exact_soln(9);

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
        norm += (exact_soln[i] * exact_soln[i]) * W;
        error += ((exact_soln[i] -
                   this->interpolated_u_axisymmetric_linear_elasticity(s, i)) *
                  (exact_soln[i] -
                   this->interpolated_u_axisymmetric_linear_elasticity(s, i))) *
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
        outfile << exact_soln[i] -
                     this->interpolated_u_axisymmetric_linear_elasticity(s, i)
                << " ";
      }
      outfile << std::endl;
    }
  }

  //======================================================================
  /// Validate against exact solution
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points and compute L2 error
  /// and L2 norm of velocity solution over element.
  //=======================================================================
  void AxisymmetricLinearElasticityEquations::compute_error(
    std::ostream& outfile,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
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

    // Exact solution Vector (u_r, u_z, u_theta)
    Vector<double> exact_soln(9);

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
      (*exact_soln_pt)(time, x, exact_soln);

      // Displacement error
      for (unsigned i = 0; i < 3; i++)
      {
        norm += (exact_soln[i] * exact_soln[i]) * W;
        error += ((exact_soln[i] -
                   this->interpolated_u_axisymmetric_linear_elasticity(s, i)) *
                  (exact_soln[i] -
                   this->interpolated_u_axisymmetric_linear_elasticity(s, i))) *
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
        outfile << exact_soln[i] -
                     this->interpolated_u_axisymmetric_linear_elasticity(s, i)
                << " ";
      }
      outfile << std::endl;
    }
  }

  // Instantiate required elements
  template class QAxisymmetricLinearElasticityElement<2>;
  template class QAxisymmetricLinearElasticityElement<3>;
  template class QAxisymmetricLinearElasticityElement<4>;


} // namespace oomph
