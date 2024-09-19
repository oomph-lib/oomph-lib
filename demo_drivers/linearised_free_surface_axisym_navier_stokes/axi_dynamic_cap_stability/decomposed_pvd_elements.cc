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
// Non-inline functions for elements that solve the equations of linear
// elasticity in cartesian coordinates

#include "decomposed_pvd_elements.h"


namespace oomph
{
  /// Static default value for timescale ratio (1.0 -- for natural scaling)
  template<unsigned DIM>
  double DecomposedPVDEquationsBase<DIM>::Default_lambda_sq_value = 1.0;


  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////

  //=======================================================================
  /// Compute the residuals for the discretised principle of
  /// virtual displacements.
  //=======================================================================
  template<unsigned DIM>
  void DecomposedPVDEquations<DIM>::
    fill_in_generic_contribution_to_residuals_pvd(Vector<double>& residuals,
                                                  DenseMatrix<double>& jacobian,
                                                  const unsigned& flag)
  {
#ifdef PARANOID
    // Check if the constitutive equation requires the explicit imposition of an
    // incompressibility constraint
    if (this->Constitutive_law_pt->requires_incompressibility_constraint())
    {
      throw OomphLibError(
        "PVDEquations cannot be used with incompressible constitutive laws.",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find out how many nodes there are
    const unsigned n_node = this->nnode();

    // Find out how many positional dofs there are
    const unsigned n_position_type = this->nnodal_position_type();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, DIM);

    // Set the value of Nintpt -- the number of integration points
    const unsigned n_intpt = this->integral_pt()->nweight();

    // Set the vector to hold the local coordinates in the element
    Vector<double> s(DIM);

    // Timescale ratio (non-dim density)
    double lambda_sq = this->lambda_sq();

    // Time factor
    double time_factor = 0.0;
    if (lambda_sq > 0)
    {
      //  time_factor = this->node_pt(0)->position_time_stepper_pt()->weight(2,
      //  0);
    }

    // Integer to store the local equation number
    int local_eqn = 0;

    for (unsigned sin_test_function = 0; sin_test_function < 2;
         sin_test_function++)
    {
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
        double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidxi);

        // Calculate interpolated values of the derivative of global position
        // wrt lagrangian coordinates
        DenseMatrix<double> interpolated_G(DIM);

        // Setup memory for accelerations
        Vector<double> accel(DIM);

        // Initialise to zero
        for (unsigned i = 0; i < DIM; i++)
        {
          // Initialise acclerations
          accel[i] = 0.0;
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_G(i, j) = 0.0;
          }
        }

        // Storage for Lagrangian coordinates (initialised to zero)
        Vector<double> interpolated_xi(DIM, 0.0);

        // Calculate displacements and derivatives and lagrangian coordinates
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over positional dofs
          for (unsigned k = 0; k < n_position_type; k++)
          {
            double psi_ = psi(l, k);
            // Loop over displacement components (deformed position)
            for (unsigned i = 0; i < DIM; i++)
            {
              // Calculate the Lagrangian coordinates and the accelerations
              // The lagrangian position is the same as the nodal after we
              // reset
              interpolated_xi[i] += this->nodal_position(l, i) * psi_;

              // Only compute accelerations if inertia is switched on
              if ((lambda_sq > 0.0) && (this->Unsteady))
              {
                accel[i] += this->du_pvd_dt(2, l, i, sin_test_function) * psi_;
              }

              // Loop over derivative directions
              for (unsigned j = 0; j < DIM; j++)
              {
                interpolated_G(j, i) +=
                  this->u_pvd(l, i, sin_test_function) * dpsidxi(l, k, j);
              }
            }
          }
        }

        // Get isotropic growth factor
        double gamma = 1.0;
        this->get_isotropic_growth(ipt, s, interpolated_xi, gamma);


        // Get body force at current time
        Vector<double> b(DIM);
        this->body_force(interpolated_xi, b);

        // We use Cartesian coordinates as the reference coordinate
        // system. In this case the undeformed metric tensor is always
        // the identity matrix -- stretched by the isotropic growth
        double diag_entry = pow(gamma, 2.0 / double(DIM));
        DenseMatrix<double> g(DIM);
        for (unsigned i = 0; i < DIM; i++)
        {
          for (unsigned j = 0; j < DIM; j++)
          {
            if (i == j)
            {
              g(i, j) = diag_entry;
            }
            else
            {
              g(i, j) = 0.0;
            }

            // The G currently just has the displacement component.
            // Add the undeformed part. This works because X = x + d and G = g
            // + dd/dxi
            interpolated_G(i, j) += g(i, j);
          }
        }

        // Premultiply the undeformed volume ratio (from the isotropic
        // growth), the weights and the Jacobian
        double W = gamma * w * J;

        // Declare and calculate the deformed metric tensor
        DenseMatrix<double> G(DIM);

        // Assign values of G
        for (unsigned i = 0; i < DIM; i++)
        {
          // Do upper half of matrix
          for (unsigned j = i; j < DIM; j++)
          {
            // Initialise G(i,j) to zero
            G(i, j) = 0.0;
            // Now calculate the dot product
            for (unsigned k = 0; k < DIM; k++)
            {
              G(i, j) += interpolated_G(i, k) * interpolated_G(j, k);
            }
          }
          // Matrix is symmetric so just copy lower half
          for (unsigned j = 0; j < i; j++)
          {
            G(i, j) = G(j, i);
          }
        }

        // Now calculate the stress tensor from the constitutive law
        DenseMatrix<double> sigma(DIM);
        get_stress(g, G, sigma);

        // Add pre-stress
        for (unsigned i = 0; i < DIM; i++)
        {
          for (unsigned j = 0; j < DIM; j++)
          {
            sigma(i, j) += this->prestress(i, j, interpolated_xi);
          }
        }

        // Get stress derivative by FD only needed for Jacobian
        //-----------------------------------------------------

        // Stress derivative
        RankFourTensor<double> d_stress_dG(DIM, DIM, DIM, DIM, 0.0);
        // Derivative of metric tensor w.r.t. to nodal coords
        RankFiveTensor<double> d_G_dX(
          n_node, n_position_type, DIM, DIM, DIM, 0.0);

        // Get Jacobian too?
        if (flag == 1)
        {
          // Derivative of metric tensor w.r.t. to discrete positional dofs
          // NOTE: Since G is symmetric we only compute the upper triangle
          //       and DO NOT copy the entries across. Subsequent computations
          //       must (and, in fact, do) therefore only operate with upper
          //       triangular entries
          for (unsigned ll = 0; ll < n_node; ll++)
          {
            for (unsigned kk = 0; kk < n_position_type; kk++)
            {
              for (unsigned ii = 0; ii < DIM; ii++)
              {
                for (unsigned aa = 0; aa < DIM; aa++)
                {
                  for (unsigned bb = aa; bb < DIM; bb++)
                  {
                    d_G_dX(ll, kk, ii, aa, bb) =
                      interpolated_G(aa, ii) * dpsidxi(ll, kk, bb) +
                      interpolated_G(bb, ii) * dpsidxi(ll, kk, aa);
                  }
                }
              }
            }
          }

          // Get the "upper triangular" entries of the derivatives of the stress
          // tensor with respect to G
          this->get_d_stress_dG_upper(g, G, sigma, d_stress_dG);
        }

        //=====EQUATIONS OF ELASTICITY FROM PRINCIPLE OF VIRTUAL
        // DISPLACEMENTS========

        // Loop over the test functions, nodes of the element
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop of types of dofs
          for (unsigned k = 0; k < n_position_type; k++)
          {
            // Offset for faster access
            const unsigned offset5 = dpsidxi.offset(l, k);

            // Loop over the displacement components
            for (unsigned i = 0; i < DIM; i++)
            {
              // Get the equation number
              local_eqn = this->nodal_local_eqn(
                l, this->u_index_pvd(l, i, sin_test_function));

              /*IF it's not a boundary condition*/
              if (local_eqn >= 0)
              {
                // Initialise contribution to sum
                double sum = 0.0;

                // Acceleration and body force
                sum += (lambda_sq * accel[i] - b[i]) * psi(l, k);

                // Stress term
                for (unsigned a = 0; a < DIM; a++)
                {
                  unsigned count = offset5;
                  for (unsigned b = 0; b < DIM; b++)
                  {
                    // Add the stress terms to the residuals
                    sum += sigma(a, b) * interpolated_G(a, i) *
                           dpsidxi.raw_direct_access(count);
                    ++count;
                  }
                }
                residuals[local_eqn] += W * sum;

                // Get Jacobian too?
                if (flag == 1)
                {
                  // Offset for faster access in general stress loop
                  const unsigned offset1 = d_G_dX.offset(l, k, i);

                  // Loop over the nodes of the element again
                  for (unsigned ll = 0; ll < n_node; ll++)
                  {
                    // Loop of types of dofs again
                    for (unsigned kk = 0; kk < n_position_type; kk++)
                    {
                      // Loop over the displacement components again
                      for (unsigned ii = 0; ii < DIM; ii++)
                      {
                        // Get the number of the unknown
                        int local_unknown = this->nodal_local_eqn(
                          ll, this->u_index_pvd(ll, ii, sin_test_function));

                        /*IF it's not a boundary condition*/
                        if (local_unknown >= 0)
                        {
                          // Offset for faster access in general stress loop
                          const unsigned offset2 = d_G_dX.offset(ll, kk, ii);
                          const unsigned offset4 = dpsidxi.offset(ll, kk);

                          // General stress term
                          //--------------------
                          double sum = 0.0;
                          unsigned count1 = offset1;
                          for (unsigned a = 0; a < DIM; a++)
                          {
                            // Bump up direct access because we're only
                            // accessing upper triangle
                            count1 += a;
                            for (unsigned b = a; b < DIM; b++)
                            {
                              double factor = d_G_dX.raw_direct_access(count1);
                              if (a == b) factor *= 0.5;

                              // Offset for faster access
                              unsigned offset3 = d_stress_dG.offset(a, b);
                              unsigned count2 = offset2;
                              unsigned count3 = offset3;

                              for (unsigned aa = 0; aa < DIM; aa++)
                              {
                                // Bump up direct access because we're only
                                // accessing upper triangle
                                count2 += aa;
                                count3 += aa;

                                // Only upper half of derivatives w.r.t. symm
                                // tensor
                                for (unsigned bb = aa; bb < DIM; bb++)
                                {
                                  sum += factor *
                                         d_stress_dG.raw_direct_access(count3) *
                                         d_G_dX.raw_direct_access(count2);
                                  ++count2;
                                  ++count3;
                                }
                              }
                              ++count1;
                            }
                          }

                          // Multiply by weight and add contribution
                          // (Add directly because this bit is nonsymmetric)
                          jacobian(local_eqn, local_unknown) += sum * W;

                          // Only upper triangle (no separate test for bc as
                          // local_eqn is already nonnegative)
                          if ((i == ii) && (local_unknown >= local_eqn))
                          {
                            // Initialise contribution
                            double sum = 0.0;

                            // Inertia term
                            sum +=
                              lambda_sq * time_factor * psi(ll, kk) * psi(l, k);

                            // Stress term
                            unsigned count4 = offset4;
                            for (unsigned a = 0; a < DIM; a++)
                            {
                              // Cache term
                              const double factor =
                                dpsidxi.raw_direct_access(count4); // ll ,kk
                              ++count4;

                              unsigned count5 = offset5;
                              for (unsigned b = 0; b < DIM; b++)
                              {
                                sum +=
                                  sigma(a, b) * factor *
                                  dpsidxi.raw_direct_access(count5); // l  ,k
                                ++count5;
                              }
                            }
                            // Add contribution to jacobian
                            jacobian(local_eqn, local_unknown) += sum * W;
                            // Add to lower triangular section
                            if (local_eqn != local_unknown)
                            {
                              jacobian(local_unknown, local_eqn) += sum * W;
                            }
                          }

                        } // End of if not boundary condition
                      }
                    }
                  }
                }

              } // End of if not boundary condition
            } // End of loop over coordinate directions
          } // End of loop over type of dof
        } // End of loop over shape functions
      } // End of loop over integration points
    }
  }


  //=======================================================================
  /// Output: x,y,[z],u,v,[w]
  //=======================================================================
  template<unsigned DIM>
  void DecomposedPVDEquations<DIM>::output(std::ostream& outfile,
                                           const unsigned& nplot)
  {
    // Set output Vector
    Vector<double> s(DIM);
    Vector<double> x(DIM);
    Vector<double> u1(DIM);
    Vector<double> u2(DIM);

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
      this->interpolated_u_pvd(s, 0, u1);
      this->interpolated_u_pvd(s, 1, u2);

      // Output the x,y,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }

      // Output u,v,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << u1[i] << " ";
      }
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << u2[i] << " ";
      }

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////

  //=====================================================================
  /// "Magic" number that indicates that the solid pressure is not stored
  /// at a node. It is a negative number that cannot be -1 because that is
  /// used to represent the positional hanging scheme in Hanging_pt objects
  //======================================================================
  template<unsigned DIM>
  int DecomposedPVDEquationsBase<DIM>::Solid_pressure_not_stored_at_node = -100;


  // Instantiate the required elements
  template class DecomposedPVDEquationsBase<1>;
  template class DecomposedPVDEquations<1>;

  template class DecomposedPVDEquationsBase<2>;
  template class DecomposedPVDEquations<2>;

  // template class QPVDElement<3, 3>;
  template class DecomposedPVDEquationsBase<3>;
  template class DecomposedPVDEquations<3>;


} // namespace oomph
