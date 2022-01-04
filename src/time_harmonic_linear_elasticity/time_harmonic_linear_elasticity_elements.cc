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

#include "time_harmonic_linear_elasticity_elements.h"


namespace oomph
{
  /// Static default value for square of frequency
  template<unsigned DIM>
  double
    TimeHarmonicLinearElasticityEquationsBase<DIM>::Default_omega_sq_value =
      1.0;


  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////

  //======================================================================
  /// Compute the strain tensor at local coordinate s
  //======================================================================
  template<unsigned DIM>
  void TimeHarmonicLinearElasticityEquationsBase<DIM>::get_strain(
    const Vector<double>& s, DenseMatrix<std::complex<double>>& strain) const
  {
#ifdef PARANOID
    if ((strain.ncol() != DIM) || (strain.nrow() != DIM))
    {
      std::ostringstream error_message;
      error_message << "Strain matrix is " << strain.ncol() << " x "
                    << strain.nrow() << ", but dimension of the equations is "
                    << DIM << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Find out how many position types there are
    unsigned n_position_type = this->nnodal_position_type();

    if (n_position_type != 1)
    {
      throw OomphLibError("TimeHarmonicLinearElasticity is not yet implemented "
                          "for more than one position type",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Find the indices at which the local velocities are stored
    std::complex<unsigned> u_nodal_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_nodal_index[i] = u_index_time_harmonic_linear_elasticity(i);
    }

    // Set up memory for the shape and derivative functions
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

    // Call the derivatives of the shape functions
    (void)dshape_eulerian(s, psi, dpsidx);

    // Calculate interpolated values of the derivative of global position
    DenseMatrix<std::complex<double>> interpolated_dudx(
      DIM, DIM, std::complex<double>(0.0, 0.0));

    // Loop over nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over velocity components
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get the nodal value
        const std::complex<double> u_value =
          std::complex<double>(this->nodal_value(l, u_nodal_index[i].real()),
                               this->nodal_value(l, u_nodal_index[i].imag()));

        // Loop over derivative directions
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_dudx(i, j) += u_value * dpsidx(l, j);
        }
      }
    }

    /// Now fill in the entries of the strain tensor
    for (unsigned i = 0; i < DIM; i++)
    {
      // Do upper half of matrix
      // Note that j must be signed here for the comparison test to work
      // Also i must be cast to an int
      for (int j = (DIM - 1); j >= static_cast<int>(i); j--)
      {
        // Off diagonal terms
        if (static_cast<int>(i) != j)
        {
          strain(i, j) =
            0.5 * (interpolated_dudx(i, j) + interpolated_dudx(j, i));
        }
        // Diagonal terms will including growth factor when it comes back in
        else
        {
          strain(i, i) = interpolated_dudx(i, i);
        }
      }
      // Matrix is symmetric so just copy lower half
      for (int j = (i - 1); j >= 0; j--)
      {
        strain(i, j) = strain(j, i);
      }
    }
  }


  //======================================================================
  /// Compute the Cauchy stress tensor at local coordinate s for
  /// displacement formulation.
  //======================================================================
  template<unsigned DIM>
  void TimeHarmonicLinearElasticityEquations<DIM>::get_stress(
    const Vector<double>& s, DenseMatrix<std::complex<double>>& stress) const
  {
#ifdef PARANOID
    if ((stress.ncol() != DIM) || (stress.nrow() != DIM))
    {
      std::ostringstream error_message;
      error_message << "Stress matrix is " << stress.ncol() << " x "
                    << stress.nrow() << ", but dimension of the equations is "
                    << DIM << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Get strain
    DenseMatrix<std::complex<double>> strain(DIM, DIM);
    this->get_strain(s, strain);

    // Now fill in the entries of the stress tensor without exploiting
    // symmetry -- sorry too lazy. This fct is only used for
    // postprocessing anyway...
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        stress(i, j) = 0.0;
        for (unsigned k = 0; k < DIM; k++)
        {
          for (unsigned l = 0; l < DIM; l++)
          {
            stress(i, j) += this->E(i, j, k, l) * strain(k, l);
          }
        }
      }
    }
  }


  //=======================================================================
  /// Compute the residuals for the linear elasticity equations in
  /// cartesian coordinates. Flag indicates if we want Jacobian too.
  //=======================================================================
  template<unsigned DIM>
  void TimeHarmonicLinearElasticityEquations<DIM>::
    fill_in_generic_contribution_to_residuals_time_harmonic_linear_elasticity(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

#ifdef PARANOID
    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();

    if (n_position_type != 1)
    {
      throw OomphLibError("TimeHarmonicLinearElasticity is not yet implemented "
                          "for more than one position type",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Throw and error if an elasticity tensor has not been set
    if (this->Elasticity_tensor_pt == 0)
    {
      throw OomphLibError("No elasticity tensor set.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the indices at which the local velocities are stored
    std::complex<unsigned> u_nodal_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_nodal_index[i] = this->u_index_time_harmonic_linear_elasticity(i);
    }


    // Square of non-dimensional frequency
    const double omega_sq_local = this->omega_sq();

    // Set up memory for the shape functions
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

    // Set the value of Nintpt -- the number of integration points
    unsigned n_intpt = this->integral_pt()->nweight();

    // Set the vector to hold the local coordinates in the element
    Vector<double> s(DIM);

    // Integer to store the local equation number
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign the values of s
      for (unsigned i = 0; i < DIM; ++i)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Call the derivatives of the shape functions (and get Jacobian)
      double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

      // Storage for Eulerian coordinates (initialised to zero)
      Vector<double> interpolated_x(DIM, 0.0);

      // Displacement
      Vector<std::complex<double>> interpolated_u(
        DIM, std::complex<double>(0.0, 0.0));

      // Calculate interpolated values of the derivative of global position
      // wrt lagrangian coordinates
      DenseMatrix<std::complex<double>> interpolated_dudx(
        DIM, DIM, std::complex<double>(0.0, 0.0));

      // Calculate displacements and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over displacement components (deformed position)
        for (unsigned i = 0; i < DIM; i++)
        {
          // Calculate the Lagrangian coordinates and the accelerations
          interpolated_x[i] += this->raw_nodal_position(l, i) * psi(l);


          // Get the nodal displacements
          const std::complex<double> u_value = std::complex<double>(
            this->raw_nodal_value(l, u_nodal_index[i].real()),
            this->raw_nodal_value(l, u_nodal_index[i].imag()));

          interpolated_u[i] += u_value * psi(l);

          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsidx(l, j);
          }
        }
      }

      // Get body force at current time
      Vector<std::complex<double>> b(DIM);
      this->body_force(interpolated_x, b);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      //=====EQUATIONS OF LINEAR ELASTICITY ========

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the displacement components
        for (unsigned a = 0; a < DIM; a++)
        {
          // Get the REAL equation number
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[a].real());

          /*IF it's not a boundary condition*/
          if (local_eqn >= 0)
          {
            // Acceleration and body force
            residuals[local_eqn] +=
              (-omega_sq_local * interpolated_u[a].real() - b[a].real()) *
              psi(l) * W;

            // Stress term
            for (unsigned b = 0; b < DIM; b++)
            {
              for (unsigned c = 0; c < DIM; c++)
              {
                for (unsigned d = 0; d < DIM; d++)
                {
                  // Add the stress terms to the residuals
                  residuals[local_eqn] += this->E(a, b, c, d) *
                                          interpolated_dudx(c, d).real() *
                                          dpsidx(l, b) * W;
                }
              }
            }

            // Jacobian entries
            if (flag)
            {
              // Loop over the displacement basis functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Loop over the displacement components again
                for (unsigned c = 0; c < DIM; c++)
                {
                  local_unknown =
                    this->nodal_local_eqn(l2, u_nodal_index[c].real());
                  // If it's not pinned
                  if (local_unknown >= 0)
                  {
                    // Inertial term
                    if (a == c)
                    {
                      jacobian(local_eqn, local_unknown) -=
                        omega_sq_local * psi(l) * psi(l2) * W;
                    }

                    // Stress term
                    for (unsigned b = 0; b < DIM; b++)
                    {
                      for (unsigned d = 0; d < DIM; d++)
                      {
                        // Add the contribution to the Jacobian matrix
                        jacobian(local_eqn, local_unknown) +=
                          this->E(a, b, c, d) * dpsidx(l2, d) * dpsidx(l, b) *
                          W;
                      }
                    }
                  } // End of if not boundary condition
                }
              }
            } // End of jacobian calculation

          } // End of if not boundary condition for real eqn


          // Get the IMAG equation number
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[a].imag());

          /*IF it's not a boundary condition*/
          if (local_eqn >= 0)
          {
            // Acceleration and body force
            residuals[local_eqn] +=
              (-omega_sq_local * interpolated_u[a].imag() - b[a].imag()) *
              psi(l) * W;

            // Stress term
            for (unsigned b = 0; b < DIM; b++)
            {
              for (unsigned c = 0; c < DIM; c++)
              {
                for (unsigned d = 0; d < DIM; d++)
                {
                  // Add the stress terms to the residuals
                  residuals[local_eqn] += this->E(a, b, c, d) *
                                          interpolated_dudx(c, d).imag() *
                                          dpsidx(l, b) * W;
                }
              }
            }

            // Jacobian entries
            if (flag)
            {
              // Loop over the displacement basis functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Loop over the displacement components again
                for (unsigned c = 0; c < DIM; c++)
                {
                  local_unknown =
                    this->nodal_local_eqn(l2, u_nodal_index[c].imag());
                  // If it's not pinned
                  if (local_unknown >= 0)
                  {
                    // Inertial term
                    if (a == c)
                    {
                      jacobian(local_eqn, local_unknown) -=
                        omega_sq_local * psi(l) * psi(l2) * W;
                    }

                    // Stress term
                    for (unsigned b = 0; b < DIM; b++)
                    {
                      for (unsigned d = 0; d < DIM; d++)
                      {
                        // Add the contribution to the Jacobian matrix
                        jacobian(local_eqn, local_unknown) +=
                          this->E(a, b, c, d) * dpsidx(l2, d) * dpsidx(l, b) *
                          W;
                      }
                    }
                  } // End of if not boundary condition
                }
              }
            } // End of jacobian calculation

          } // End of if not boundary condition for imag eqn

        } // End of loop over coordinate directions
      } // End of loop over shape functions
    } // End of loop over integration points
  }

  //=======================================================================
  /// Output exact solution x,y,[z],u_r,v_r,[w_r],u_i,v_i,[w_i]
  //=======================================================================
  template<unsigned DIM>
  void TimeHarmonicLinearElasticityEquations<DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Exact solution Vector
    Vector<double> exact_soln(2 * DIM);

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

      // Output x,y,...,u_exact,...
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      for (unsigned i = 0; i < 2 * DIM; i++)
      {
        outfile << exact_soln[i] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }


  //=======================================================================
  /// Output: x,y,[z],u,v,[w]
  //=======================================================================
  template<unsigned DIM>
  void TimeHarmonicLinearElasticityEquations<DIM>::output(std::ostream& outfile,
                                                          const unsigned& nplot)
  {
    // Set output Vector
    Vector<double> s(DIM);
    Vector<double> x(DIM);
    Vector<std::complex<double>> u(DIM);

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
      this->interpolated_u_time_harmonic_linear_elasticity(s, u);

      // Output the x,y,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }

      // Output u,v,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << u[i].real() << " ";
      }

      // Output u,v,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << u[i].imag() << " ";
      }

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }


  //=======================================================================
  /// C-style output: x,y,[z],u,v,[w]
  //=======================================================================
  template<unsigned DIM>
  void TimeHarmonicLinearElasticityEquations<DIM>::output(FILE* file_pt,
                                                          const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    fprintf(file_pt, "%s", this->tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        fprintf(file_pt, "%g ", this->interpolated_x(s, i));
      }

      // Displacement
      for (unsigned i = 0; i < DIM; i++)
      {
        fprintf(
          file_pt,
          "%g ",
          this->interpolated_u_time_harmonic_linear_elasticity(s, i).real());
      }
      for (unsigned i = 0; i < DIM; i++)
      {
        fprintf(
          file_pt,
          "%g ",
          this->interpolated_u_time_harmonic_linear_elasticity(s, i).imag());
      }
    }
    fprintf(file_pt, "\n");

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(file_pt, nplot);
  }


  //=======================================================================
  /// Compute norm of the solution
  //=======================================================================
  template<unsigned DIM>
  void TimeHarmonicLinearElasticityEquations<DIM>::compute_norm(double& norm)
  {
    // Initialise
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Displacement vector
    Vector<std::complex<double>> disp(DIM);

    // Find out how many nodes there are in the element
    unsigned n_node = this->nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
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
      this->interpolated_u_time_harmonic_linear_elasticity(s, disp);

      // Add to  norm
      for (unsigned ii = 0; ii < DIM; ii++)
      {
        norm += (disp[ii].real() * disp[ii].real() +
                 disp[ii].imag() * disp[ii].imag()) *
                W;
      }
    }
  }


  // Instantiate the required elements
  template class TimeHarmonicLinearElasticityEquationsBase<2>;
  template class TimeHarmonicLinearElasticityEquations<2>;

  template class QTimeHarmonicLinearElasticityElement<3, 3>;
  template class TimeHarmonicLinearElasticityEquationsBase<3>;
  template class TimeHarmonicLinearElasticityEquations<3>;


} // namespace oomph
