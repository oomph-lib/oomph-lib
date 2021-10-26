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
// Non-inline functions for elements that solve the principle of virtual
// equations of solid mechanics

#include "solid_elements.h"


namespace oomph
{
  /// Static default value for timescale ratio (1.0 -- for natural scaling)
  template<unsigned DIM>
  double PVDEquationsBase<DIM>::Default_lambda_sq_value = 1.0;


  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////

  //======================================================================
  /// Compute the strain tensor at local coordinate s
  //======================================================================
  template<unsigned DIM>
  void PVDEquationsBase<DIM>::get_strain(const Vector<double>& s,
                                         DenseMatrix<double>& strain) const
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
#endif

    // Find out how many nodes there are in the element
    const unsigned n_node = nnode();

    // Find out how many position types there are
    const unsigned n_position_type = this->nnodal_position_type();

    // Set up memory for the shape and test functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, DIM);

    // Call the derivatives of the shape functions
    (void)dshape_lagrangian(s, psi, dpsidxi);

    // Calculate interpolated values of the derivative of global position
    DenseMatrix<double> interpolated_G(DIM);

    // Initialise to zero
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        interpolated_G(i, j) = 0.0;
      }
    }

    // Storage for Lagrangian coordinates (initialised to zero)
    Vector<double> interpolated_xi(DIM, 0.0);

    // Loop over nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over the positional dofs
      for (unsigned k = 0; k < n_position_type; k++)
      {
        // Loop over velocity components
        for (unsigned i = 0; i < DIM; i++)
        {
          // Calculate the Lagrangian coordinates
          interpolated_xi[i] +=
            this->lagrangian_position_gen(l, k, i) * psi(l, k);

          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_G(j, i) +=
              this->nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
          }
        }
      }
    }

    // Get isotropic growth factor
    double gamma = 1.0;
    // Dummy integration point
    unsigned ipt = 0;
    get_isotropic_growth(ipt, s, interpolated_xi, gamma);

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
      }
    }


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

    // Fill in the strain tensor
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        strain(i, j) = 0.5 * (G(i, j) - g(i, j));
      }
    }
  }

  //=======================================================================
  /// Compute the residuals for the discretised principle of
  /// virtual displacements.
  //=======================================================================
  template<unsigned DIM>
  void PVDEquations<DIM>::fill_in_generic_contribution_to_residuals_pvd(
    Vector<double>& residuals,
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

    // Simply set up initial condition?
    if (this->Solid_ic_pt != 0)
    {
      this->fill_in_residuals_for_solid_ic(residuals);
      return;
    }

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
      time_factor = this->node_pt(0)->position_time_stepper_pt()->weight(2, 0);
    }

    // Integer to store the local equation number
    int local_eqn = 0;

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
      double J = this->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);

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
            interpolated_xi[i] += this->lagrangian_position_gen(l, k, i) * psi_;

            // Only compute accelerations if inertia is switched on
            if ((lambda_sq > 0.0) && (this->Unsteady))
            {
              accel[i] += this->dnodal_position_gen_dt(2, l, k, i) * psi_;
            }

            // Loop over derivative directions
            for (unsigned j = 0; j < DIM; j++)
            {
              interpolated_G(j, i) +=
                this->nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
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
            local_eqn = this->position_local_eqn(l, k, i);

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
                      int local_unknown = this->position_local_eqn(ll, kk, ii);

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
                              sum += sigma(a, b) * factor *
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


  //=======================================================================
  /// Output: x,y,[z],xi0,xi1,[xi2],gamma
  //=======================================================================
  template<unsigned DIM>
  void PVDEquations<DIM>::output(std::ostream& outfile, const unsigned& n_plot)
  {
    Vector<double> x(DIM);
    Vector<double> xi(DIM);
    Vector<double> s(DIM);

    // Tecplot header info
    outfile << this->tecplot_zone_string(n_plot);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(n_plot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, n_plot, s);

      // Get Eulerian and Lagrangian coordinates
      this->interpolated_x(s, x);
      this->interpolated_xi(s, xi);

      // Get isotropic growth
      double gamma;
      // Dummy integration point
      unsigned ipt = 0;
      this->get_isotropic_growth(ipt, s, xi, gamma);

      // Output the x,y,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }

      // Output xi0,xi1,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << xi[i] << " ";
      }

      // Output growth
      outfile << gamma;
      outfile << std::endl;
    }


    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, n_plot);
    outfile << std::endl;
  }


  //=======================================================================
  /// C-style output: x,y,[z],xi0,xi1,[xi2],gamma
  //=======================================================================
  template<unsigned DIM>
  void PVDEquations<DIM>::output(FILE* file_pt, const unsigned& n_plot)
  {
    // Set output Vector
    Vector<double> s(DIM);
    Vector<double> x(DIM);
    Vector<double> xi(DIM);

    switch (DIM)
    {
      case 2:

        // Tecplot header info
        // outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
        fprintf(file_pt, "ZONE I=%i, J=%i\n", n_plot, n_plot);

        // Loop over element nodes
        for (unsigned l2 = 0; l2 < n_plot; l2++)
        {
          s[1] = -1.0 + l2 * 2.0 / (n_plot - 1);
          for (unsigned l1 = 0; l1 < n_plot; l1++)
          {
            s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

            // Get Eulerian and Lagrangian coordinates
            this->interpolated_x(s, x);
            this->interpolated_xi(s, xi);

            // Get isotropic growth
            double gamma;
            // Dummy integration point
            unsigned ipt = 0;
            this->get_isotropic_growth(ipt, s, xi, gamma);

            // Output the x,y,..
            for (unsigned i = 0; i < DIM; i++)
            {
              // outfile << x[i] << " ";
              fprintf(file_pt, "%g ", x[i]);
            }
            // Output xi0,xi1,..
            for (unsigned i = 0; i < DIM; i++)
            {
              // outfile << xi[i] << " ";
              fprintf(file_pt, "%g ", xi[i]);
            }
            // Output growth
            // outfile << gamma << " ";
            // outfile << std::endl;
            fprintf(file_pt, "%g \n", gamma);
          }
        }
        // outfile << std::endl;
        fprintf(file_pt, "\n");

        break;

      case 3:

        // Tecplot header info
        // outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
        fprintf(file_pt, "ZONE I=%i, J=%i, K=%i \n", n_plot, n_plot, n_plot);

        // Loop over element nodes
        for (unsigned l3 = 0; l3 < n_plot; l3++)
        {
          s[2] = -1.0 + l3 * 2.0 / (n_plot - 1);
          for (unsigned l2 = 0; l2 < n_plot; l2++)
          {
            s[1] = -1.0 + l2 * 2.0 / (n_plot - 1);
            for (unsigned l1 = 0; l1 < n_plot; l1++)
            {
              s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

              // Get Eulerian and Lagrangian coordinates
              this->interpolated_x(s, x);
              this->interpolated_xi(s, xi);

              // Get isotropic growth
              double gamma;
              // Dummy integration point
              unsigned ipt = 0;
              this->get_isotropic_growth(ipt, s, xi, gamma);

              // Output the x,y,z
              for (unsigned i = 0; i < DIM; i++)
              {
                // outfile << x[i] << " ";
                fprintf(file_pt, "%g ", x[i]);
              }
              // Output xi0,xi1,xi2
              for (unsigned i = 0; i < DIM; i++)
              {
                // outfile << xi[i] << " ";
                fprintf(file_pt, "%g ", xi[i]);
              }
              // Output growth
              // outfile << gamma << " ";
              // outfile << std::endl;
              fprintf(file_pt, "%g \n", gamma);
            }
          }
        }
        // outfile << std::endl;
        fprintf(file_pt, "\n");

        break;

      default:
        std::ostringstream error_message;
        error_message << "No output routine for PVDEquations<" << DIM
                      << "> elements --  write it yourself!" << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }


  //=======================================================================
  /// Output: x,y,[z],xi0,xi1,[xi2],gamma strain and stress components
  //=======================================================================
  template<unsigned DIM>
  void PVDEquations<DIM>::extended_output(std::ostream& outfile,
                                          const unsigned& n_plot)
  {
    Vector<double> x(DIM);
    Vector<double> xi(DIM);
    Vector<double> s(DIM);
    DenseMatrix<double> stress_or_strain(DIM, DIM);

    // Tecplot header info
    outfile << this->tecplot_zone_string(n_plot);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(n_plot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, n_plot, s);

      // Get Eulerian and Lagrangian coordinates
      this->interpolated_x(s, x);
      this->interpolated_xi(s, xi);

      // Get isotropic growth
      double gamma;
      // Dummy integration point
      unsigned ipt = 0;
      this->get_isotropic_growth(ipt, s, xi, gamma);

      // Output the x,y,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }

      // Output xi0,xi1,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << xi[i] << " ";
      }

      // Output growth
      outfile << gamma << " ";

      // get the strain
      this->get_strain(s, stress_or_strain);
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = 0; j <= i; j++)
        {
          outfile << stress_or_strain(j, i) << " ";
        }
      }

      // get the stress
      this->get_stress(s, stress_or_strain);
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = 0; j <= i; j++)
        {
          outfile << stress_or_strain(j, i) << " ";
        }
      }


      outfile << std::endl;
    }


    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, n_plot);
    outfile << std::endl;
  }


  //=======================================================================
  /// Get potential (strain) and kinetic energy
  //=======================================================================
  template<unsigned DIM>
  void PVDEquationsBase<DIM>::get_energy(double& pot_en, double& kin_en)
  {
    // Initialise
    pot_en = 0;
    kin_en = 0;

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Find out how many nodes there are
    const unsigned n_node = this->nnode();

    // Find out how many positional dofs there are
    const unsigned n_position_type = this->nnodal_position_type();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, DIM);

    // Timescale ratio (non-dim density)
    double lambda_sq = this->lambda_sq();

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

      // Call the derivatives of the shape functions and get Jacobian
      double J = this->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);

      // Storage for Lagrangian coordinates and velocity (initialised to zero)
      Vector<double> interpolated_xi(DIM, 0.0);
      Vector<double> veloc(DIM, 0.0);

      // Calculate lagrangian coordinates
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Loop over displacement components (deformed position)
          for (unsigned i = 0; i < DIM; i++)
          {
            // Calculate the Lagrangian coordinates
            interpolated_xi[i] +=
              this->lagrangian_position_gen(l, k, i) * psi(l, k);

            // Calculate the velocity components (if unsteady solve)
            if (this->Unsteady)
            {
              veloc[i] += this->dnodal_position_gen_dt(l, k, i) * psi(l, k);
            }
          }
        }
      }

      // Get isotropic growth factor
      double gamma = 1.0;
      this->get_isotropic_growth(ipt, s, interpolated_xi, gamma);

      // Premultiply the undeformed volume ratio (from the isotropic
      // growth), the weights and the Jacobian
      double W = gamma * w * J;

      DenseMatrix<double> sigma(DIM, DIM);
      DenseMatrix<double> strain(DIM, DIM);

      // Now calculate the stress tensor from the constitutive law
      this->get_stress(s, sigma);

      // Add pre-stress
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = 0; j < DIM; j++)
        {
          sigma(i, j) += prestress(i, j, interpolated_xi);
        }
      }

      // get the strain
      this->get_strain(s, strain);

      // Initialise
      double local_pot_en = 0;
      double veloc_sq = 0;

      // Compute integrals
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = 0; j < DIM; j++)
        {
          local_pot_en += sigma(i, j) * strain(i, j);
        }
        veloc_sq += veloc[i] * veloc[i];
      }

      pot_en += 0.5 * local_pot_en * W;
      kin_en += lambda_sq * 0.5 * veloc_sq * W;
    }
  }


  //=======================================================================
  /// Compute the contravariant second Piola Kirchoff stress at a given local
  /// coordinate. Note: this replicates a lot of code that is already
  /// coontained in get_residuals() but without sacrificing efficiency
  /// (re-computing the shape functions several times) or creating
  /// helper functions with horrendous interfaces (to pass all the
  /// functions which shouldn't be recomputed) about this is
  /// unavoidable.
  //=======================================================================
  template<unsigned DIM>
  void PVDEquations<DIM>::get_stress(const Vector<double>& s,
                                     DenseMatrix<double>& sigma)
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, DIM);

    // Call the derivatives of the shape functions (ignore Jacobian)
    (void)this->dshape_lagrangian(s, psi, dpsidxi);

    // Lagrangian coordinates
    Vector<double> xi(DIM);
    this->interpolated_xi(s, xi);

    // Get isotropic growth factor
    double gamma;
    // Dummy integration point
    unsigned ipt = 0;
    this->get_isotropic_growth(ipt, s, xi, gamma);

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
      }
    }


    // Calculate interpolated values of the derivative of global position
    // wrt lagrangian coordinates
    DenseMatrix<double> interpolated_G(DIM);

    // Initialise to zero
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        interpolated_G(i, j) = 0.0;
      }
    }

    // Calculate displacements and derivatives
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over positional dofs
      for (unsigned k = 0; k < n_position_type; k++)
      {
        // Loop over displacement components (deformed position)
        for (unsigned i = 0; i < DIM; i++)
        {
          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_G(j, i) +=
              this->nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
          }
        }
      }
    }

    // Declare and calculate the deformed metric tensor
    DenseMatrix<double> G(DIM);
    // Assign values of G
    for (unsigned i = 0; i < DIM; i++)
    {
      // Do upper half of matrix
      // Note that j must be signed here for the comparison test to work
      // Also i must be cast to an int
      for (int j = (DIM - 1); j >= static_cast<int>(i); j--)
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
      for (int j = (i - 1); j >= 0; j--)
      {
        G(i, j) = G(j, i);
      }
    }

    // Now calculate the stress tensor from the constitutive law
    get_stress(g, G, sigma);
  }


  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////


  //=======================================================================
  /// Compute principal stress vectors and (scalar) principal stresses
  /// at specified local coordinate: \c  principal_stress_vector(i,j)
  /// is the j-th component of the i-th principal stress vector.
  //=======================================================================
  template<unsigned DIM>
  void PVDEquationsBase<DIM>::get_principal_stress(
    const Vector<double>& s,
    DenseMatrix<double>& principal_stress_vector,
    Vector<double>& principal_stress)
  {
    // Compute contravariant ("upper") 2nd Piola Kirchhoff stress
    DenseDoubleMatrix sigma(DIM, DIM);
    get_stress(s, sigma);

    // Lagrangian coordinates
    Vector<double> xi(DIM);
    this->interpolated_xi(s, xi);

    // Add pre-stress
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        sigma(i, j) += this->prestress(i, j, xi);
      }
    }

    // Get covariant base vectors in deformed configuration
    DenseMatrix<double> lower_deformed_basis(DIM);
    get_deformed_covariant_basis_vectors(s, lower_deformed_basis);

    // Work out covariant ("lower") metric tensor
    DenseDoubleMatrix lower_metric(DIM);
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        lower_metric(i, j) = 0.0;
        for (unsigned k = 0; k < DIM; k++)
        {
          lower_metric(i, j) +=
            lower_deformed_basis(i, k) * lower_deformed_basis(j, k);
        }
      }
    }

    // Work out cartesian components of contravariant ("upper") basis vectors
    DenseMatrix<double> upper_deformed_basis(DIM);

    // Loop over RHSs
    Vector<double> rhs(DIM);
    Vector<double> aux(DIM);
    for (unsigned k = 0; k < DIM; k++)
    {
      for (unsigned l = 0; l < DIM; l++)
      {
        rhs[l] = lower_deformed_basis(l, k);
      }

      lower_metric.solve(rhs, aux);

      for (unsigned l = 0; l < DIM; l++)
      {
        upper_deformed_basis(l, k) = aux[l];
      }
    }

    // Eigenvalues (=principal stresses) and eigenvectors
    DenseMatrix<double> ev(DIM);

    // Get eigenvectors of contravariant 2nd Piola Kirchoff stress
    sigma.eigenvalues_by_jacobi(principal_stress, ev);

    // ev(j,i) is the i-th component of the j-th eigenvector
    // relative to the deformed "lower variance" basis!
    // Work out cartesian components of eigenvectors by multiplying
    // the "lower variance components" by these "upper variance" basis
    // vectors

    // Loop over cartesian compnents
    for (unsigned i = 0; i < DIM; i++)
    {
      // Initialise the row
      for (unsigned j = 0; j < DIM; j++)
      {
        principal_stress_vector(j, i) = 0.0;
      }

      // Loop over basis vectors
      for (unsigned j = 0; j < DIM; j++)
      {
        for (unsigned k = 0; k < DIM; k++)
        {
          principal_stress_vector(j, i) +=
            upper_deformed_basis(k, i) * ev(j, k);
        }
      }
    }

    // Scaling factor to turn these vectors into unit vectors
    Vector<double> norm(DIM);
    for (unsigned i = 0; i < DIM; i++)
    {
      norm[i] = 0.0;
      for (unsigned j = 0; j < DIM; j++)
      {
        norm[i] += pow(principal_stress_vector(i, j), 2);
      }
      norm[i] = sqrt(norm[i]);
    }


    // Scaling and then multiplying by eigenvalue gives the principal stress
    // vectors
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        principal_stress_vector(j, i) =
          ev(j, i) / norm[j] * principal_stress[j];
      }
    }
  }


  //=======================================================================
  /// Return the deformed covariant basis vectors
  /// at specified local coordinate:  \c def_covariant_basis(i,j)
  /// is the j-th component of the i-th basis vector.
  //=======================================================================
  template<unsigned DIM>
  void PVDEquationsBase<DIM>::get_deformed_covariant_basis_vectors(
    const Vector<double>& s, DenseMatrix<double>& def_covariant_basis)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, DIM);


    // Call the derivatives of the shape functions (ignore Jacobian)
    (void)dshape_lagrangian(s, psi, dpsidxi);


    // Initialise to zero
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        def_covariant_basis(i, j) = 0.0;
      }
    }

    // Calculate displacements and derivatives
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over positional dofs
      for (unsigned k = 0; k < n_position_type; k++)
      {
        // Loop over displacement components (deformed position)
        for (unsigned i = 0; i < DIM; i++)
        {
          // Loop over derivative directions (i.e. base vectors)
          for (unsigned j = 0; j < DIM; j++)
          {
            def_covariant_basis(j, i) +=
              nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
          }
        }
      }
    }
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
  int PVDEquationsBase<DIM>::Solid_pressure_not_stored_at_node = -100;


  //=======================================================================
  /// Fill in element's contribution to the elemental
  /// residual vector and/or Jacobian matrix.
  /// flag=0: compute only residual vector
  /// flag=1: compute both, fully analytically
  /// flag=2: compute both, using FD for the derivatives w.r.t. to the
  ///         discrete displacment dofs.
  /// flag=3: compute residuals, jacobian (full analytic) and mass matrix
  /// flag=4: compute residuals, jacobian (FD for derivatives w.r.t.
  ///          displacements) and mass matrix
  //=======================================================================
  template<unsigned DIM>
  void PVDEquationsWithPressure<DIM>::
    fill_in_generic_residual_contribution_pvd_with_pressure(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      const unsigned& flag)
  {
#ifdef PARANOID
    // Check if the constitutive equation requires the explicit imposition of an
    // incompressibility constraint
    if (this->Constitutive_law_pt->requires_incompressibility_constraint() &&
        (!Incompressible))
    {
      throw OomphLibError("The constitutive law requires the use of the "
                          "incompressible formulation by setting the element's "
                          "member function set_incompressible()",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Simply set up initial condition?
    if (this->Solid_ic_pt != 0)
    {
      this->get_residuals_for_solid_ic(residuals);
      return;
    }

    // Find out how many nodes there are
    const unsigned n_node = this->nnode();

    // Find out how many position types of dof there are
    const unsigned n_position_type = this->nnodal_position_type();

    // Find out how many pressure dofs there are
    const unsigned n_solid_pres = this->npres_solid();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, DIM);

    // Set up memory for the pressure shape functions
    Shape psisp(n_solid_pres);

    // Set the value of n_intpt
    const unsigned n_intpt = this->integral_pt()->nweight();

    // Set the vector to hold the local coordinates in the element
    Vector<double> s(DIM);

    // Timescale ratio (non-dim density)
    double lambda_sq = this->lambda_sq();

    // Time factor
    double time_factor = 0.0;
    if (lambda_sq > 0)
    {
      time_factor = this->node_pt(0)->position_time_stepper_pt()->weight(2, 0);
    }

    // Integers to hold the local equation and unknown numbers
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

      // Call the derivatives of the shape functions
      double J = this->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);

      // Call the pressure shape functions
      solid_pshape_at_knot(ipt, psisp);

      // Storage for Lagrangian coordinates (initialised to zero)
      Vector<double> interpolated_xi(DIM, 0.0);

      // Deformed tangent vectors
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

      // Calculate displacements and derivatives and lagrangian coordinates
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Loop over displacement components (deformed position)
          for (unsigned i = 0; i < DIM; i++)
          {
            // Calculate the lagrangian coordinates and the accelerations
            interpolated_xi[i] +=
              this->lagrangian_position_gen(l, k, i) * psi(l, k);

            // Only compute accelerations if inertia is switched on
            // otherwise the timestepper might not be able to
            // work out dx_gen_dt(2,...)
            if ((lambda_sq > 0.0) && (this->Unsteady))
            {
              accel[i] += this->dnodal_position_gen_dt(2, l, k, i) * psi(l, k);
            }

            // Loop over derivative directions
            for (unsigned j = 0; j < DIM; j++)
            {
              interpolated_G(j, i) +=
                this->nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
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
        }
      }

      // Premultiply the undeformed volume ratio (from the isotropic
      // growth), the weights and the Jacobian
      double W = gamma * w * J;

      // Calculate the interpolated solid pressure
      double interpolated_solid_p = 0.0;
      for (unsigned l = 0; l < n_solid_pres; l++)
      {
        interpolated_solid_p += solid_p(l) * psisp[l];
      }


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

      // Now calculate the deviatoric stress and all pressure-related
      // quantitites
      DenseMatrix<double> sigma(DIM, DIM), sigma_dev(DIM, DIM), Gup(DIM, DIM);
      double detG = 0.0;
      double gen_dil = 0.0;
      double inv_kappa = 0.0;

      // Get stress derivative by FD only needed for Jacobian

      // Stress etc derivatives
      RankFourTensor<double> d_stress_dG(DIM, DIM, DIM, DIM, 0.0);
      DenseMatrix<double> d_detG_dG(DIM, DIM, 0.0);
      DenseMatrix<double> d_gen_dil_dG(DIM, DIM, 0.0);

      // Derivative of metric tensor w.r.t. to nodal coords
      RankFiveTensor<double> d_G_dX(
        n_node, n_position_type, DIM, DIM, DIM, 0.0);

      // Get Jacobian too?
      if ((flag == 1) || (flag == 3))
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
      }


      // Incompressible: Compute the deviatoric part of the stress tensor, the
      // contravariant deformed metric tensor and the determinant
      // of the deformed covariant metric tensor.
      if (Incompressible)
      {
        get_stress(g, G, sigma_dev, Gup, detG);

        // Get full stress
        for (unsigned a = 0; a < DIM; a++)
        {
          for (unsigned b = 0; b < DIM; b++)
          {
            sigma(a, b) = sigma_dev(a, b) - interpolated_solid_p * Gup(a, b);
          }
        }

        // Get Jacobian too?
        if ((flag == 1) || (flag == 3))
        {
          // Get the "upper triangular" entries of the derivatives of the stress
          // tensor with respect to G
          this->get_d_stress_dG_upper(
            g, G, sigma, detG, interpolated_solid_p, d_stress_dG, d_detG_dG);
        }
      }
      // Nearly incompressible: Compute the deviatoric part of the
      // stress tensor, the contravariant deformed metric tensor,
      // the generalised dilatation and the inverse bulk modulus.
      else
      {
        get_stress(g, G, sigma_dev, Gup, gen_dil, inv_kappa);

        // Get full stress
        for (unsigned a = 0; a < DIM; a++)
        {
          for (unsigned b = 0; b < DIM; b++)
          {
            sigma(a, b) = sigma_dev(a, b) - interpolated_solid_p * Gup(a, b);
          }
        }

        // Get Jacobian too?
        if ((flag == 1) || (flag == 3))
        {
          // Get the "upper triangular" entries of the derivatives of the stress
          // tensor with respect to G
          this->get_d_stress_dG_upper(g,
                                      G,
                                      sigma,
                                      gen_dil,
                                      inv_kappa,
                                      interpolated_solid_p,
                                      d_stress_dG,
                                      d_gen_dil_dG);
        }
      }

      // Add pre-stress
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = 0; j < DIM; j++)
        {
          sigma(i, j) += this->prestress(i, j, interpolated_xi);
        }
      }

      //=====EQUATIONS OF ELASTICITY FROM PRINCIPLE OF VIRTUAL
      // DISPLACEMENTS========

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the types of dof
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Offset for faster access
          const unsigned offset5 = dpsidxi.offset(l, k);

          // Loop over the displacement components
          for (unsigned i = 0; i < DIM; i++)
          {
            // Get the equation number
            local_eqn = this->position_local_eqn(l, k, i);

            /*IF it's not a boundary condition*/
            if (local_eqn >= 0)
            {
              // Initialise the contribution
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

              // Add in the mass matrix terms
              if (flag > 2)
              {
                // Loop over the nodes of the element again
                for (unsigned ll = 0; ll < n_node; ll++)
                {
                  // Loop of types of dofs again
                  for (unsigned kk = 0; kk < n_position_type; kk++)
                  {
                    // Get the number of the unknown
                    int local_unknown = this->position_local_eqn(ll, kk, i);

                    /*IF it's not a boundary condition*/
                    if (local_unknown >= 0)
                    {
                      mass_matrix(local_eqn, local_unknown) +=
                        lambda_sq * psi(l, k) * psi(ll, kk) * W;
                    }
                  }
                }
              }

              // Add in the jacobian terms
              if ((flag == 1) || (flag == 3))
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
                      int local_unknown = this->position_local_eqn(ll, kk, ii);

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
                          // Initialise the contribution
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
                              sum += sigma(a, b) * factor *
                                     dpsidxi.raw_direct_access(count5); // l  ,k
                              ++count5;
                            }
                          }

                          // Add to jacobian
                          jacobian(local_eqn, local_unknown) += sum * W;
                          // Add to lower triangular parts
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

              // Derivatives w.r.t. pressure dofs
              if (flag > 0)
              {
                // Loop over the pressure dofs for unknowns
                for (unsigned l2 = 0; l2 < n_solid_pres; l2++)
                {
                  local_unknown = this->solid_p_local_eqn(l2);

                  // If it's not a boundary condition
                  if (local_unknown >= 0)
                  {
                    // Add the pressure terms to the jacobian
                    for (unsigned a = 0; a < DIM; a++)
                    {
                      for (unsigned b = 0; b < DIM; b++)
                      {
                        jacobian(local_eqn, local_unknown) -=
                          psisp[l2] * Gup(a, b) * interpolated_G(a, i) *
                          dpsidxi(l, k, b) * W;
                      }
                    }
                  }
                }
              } // End of Jacobian terms

            } // End of if not boundary condition
          } // End of loop over coordinate directions
        } // End of loop over types of dof
      } // End of loop over shape functions

      //==============CONSTRAINT EQUATIONS FOR PRESSURE=====================

      // Now loop over the pressure degrees of freedom
      for (unsigned l = 0; l < n_solid_pres; l++)
      {
        local_eqn = this->solid_p_local_eqn(l);

        // Pinned (unlikely, actually) or real dof?
        if (local_eqn >= 0)
        {
          // For true incompressibility we need to conserve volume
          // so the determinant of the deformed metric tensor
          // needs to be equal to that of the undeformed one, which
          // is equal to the volumetric growth factor
          if (Incompressible)
          {
            residuals[local_eqn] += (detG - gamma) * psisp[l] * W;


            // Get Jacobian too?
            if ((flag == 1) || (flag == 3))
            {
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
                    int local_unknown = this->position_local_eqn(ll, kk, ii);

                    /*IF it's not a boundary condition*/
                    if (local_unknown >= 0)
                    {
                      // Offset for faster access
                      const unsigned offset = d_G_dX.offset(ll, kk, ii);

                      // General stress term
                      double sum = 0.0;
                      unsigned count = offset;
                      for (unsigned aa = 0; aa < DIM; aa++)
                      {
                        // Bump up direct access because we're only
                        // accessing upper triangle
                        count += aa;

                        // Only upper half
                        for (unsigned bb = aa; bb < DIM; bb++)
                        {
                          sum += d_detG_dG(aa, bb) *
                                 d_G_dX.raw_direct_access(count) * psisp(l);
                          ++count;
                        }
                      }
                      jacobian(local_eqn, local_unknown) += sum * W;
                    }
                  }
                }
              }

              // No Jacobian terms due to pressure since it does not feature
              // in the incompressibility constraint
            }
          }
          // Nearly incompressible: (Neg.) pressure given by product of
          // bulk modulus and generalised dilatation
          else
          {
            residuals[local_eqn] +=
              (inv_kappa * interpolated_solid_p + gen_dil) * psisp[l] * W;

            // Add in the jacobian terms
            if ((flag == 1) || (flag == 3))
            {
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
                    int local_unknown = this->position_local_eqn(ll, kk, ii);

                    /*IF it's not a boundary condition*/
                    if (local_unknown >= 0)
                    {
                      // Offset for faster access
                      const unsigned offset = d_G_dX.offset(ll, kk, ii);

                      // General stress term
                      double sum = 0.0;
                      unsigned count = offset;
                      for (unsigned aa = 0; aa < DIM; aa++)
                      {
                        // Bump up direct access because we're only
                        // accessing upper triangle
                        count += aa;

                        // Only upper half
                        for (unsigned bb = aa; bb < DIM; bb++)
                        {
                          sum += d_gen_dil_dG(aa, bb) *
                                 d_G_dX.raw_direct_access(count) * psisp(l);
                          ++count;
                        }
                      }
                      jacobian(local_eqn, local_unknown) += sum * W;
                    }
                  }
                }
              }
            }

            // Derivatives w.r.t. pressure dofs
            if (flag > 0)
            {
              // Loop over the pressure nodes again
              for (unsigned l2 = 0; l2 < n_solid_pres; l2++)
              {
                local_unknown = this->solid_p_local_eqn(l2);
                // If not pinnned
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    inv_kappa * psisp[l2] * psisp[l] * W;
                }
              }
            } // End of jacobian terms

          } // End of else

        } // End of if not boundary condition
      }

    } // End of loop over integration points
  }


  //=======================================================================
  /// Output: x,y,[z],xi0,xi1,[xi2],p,gamma
  //=======================================================================
  template<unsigned DIM>
  void PVDEquationsWithPressure<DIM>::output(std::ostream& outfile,
                                             const unsigned& n_plot)
  {
    // Set output Vector
    Vector<double> s(DIM);
    Vector<double> x(DIM);
    Vector<double> xi(DIM);

    // Tecplot header info
    outfile << this->tecplot_zone_string(n_plot);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(n_plot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, n_plot, s);

      // Get Eulerian and Lagrangian coordinates
      this->interpolated_x(s, x);
      this->interpolated_xi(s, xi);

      // Get isotropic growth
      double gamma;
      // Dummy integration point
      unsigned ipt = 0;
      this->get_isotropic_growth(ipt, s, xi, gamma);

      // Output the x,y,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }

      // Output xi0,xi1,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << xi[i] << " ";
      }

      // Output growth
      outfile << gamma << " ";
      // Output pressure
      outfile << interpolated_solid_p(s) << " ";
      outfile << "\n";
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, n_plot);
  }


  //=======================================================================
  /// C-stsyle output: x,y,[z],xi0,xi1,[xi2],p,gamma
  //=======================================================================
  template<unsigned DIM>
  void PVDEquationsWithPressure<DIM>::output(FILE* file_pt,
                                             const unsigned& n_plot)
  {
    // Set output Vector
    Vector<double> s(DIM);
    Vector<double> x(DIM);
    Vector<double> xi(DIM);

    switch (DIM)
    {
      case 2:
        // Tecplot header info
        // outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
        fprintf(file_pt, "ZONE I=%i, J=%i\n", n_plot, n_plot);

        // Loop over element nodes
        for (unsigned l2 = 0; l2 < n_plot; l2++)
        {
          s[1] = -1.0 + l2 * 2.0 / (n_plot - 1);
          for (unsigned l1 = 0; l1 < n_plot; l1++)
          {
            s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

            // Get Eulerian and Lagrangian coordinates
            this->interpolated_x(s, x);
            this->interpolated_xi(s, xi);

            // Get isotropic growth
            double gamma;
            // Dummy integration point
            unsigned ipt = 0;
            this->get_isotropic_growth(ipt, s, xi, gamma);

            // Output the x,y,..
            for (unsigned i = 0; i < DIM; i++)
            {
              // outfile << x[i] << " ";
              fprintf(file_pt, "%g ", x[i]);
            }
            // Output xi0,xi1,..
            for (unsigned i = 0; i < DIM; i++)
            {
              // outfile << xi[i] << " ";
              fprintf(file_pt, "%g ", xi[i]);
            }
            // Output growth
            // outfile << gamma << " ";
            fprintf(file_pt, "%g ", gamma);

            // Output pressure
            // outfile << interpolated_solid_p(s) << " ";
            // outfile << std::endl;
            fprintf(file_pt, "%g \n", interpolated_solid_p(s));
          }
        }

        break;

      case 3:
        // Tecplot header info
        // outfile << "ZONE I=" << n_plot
        //        << ", J=" << n_plot
        //        << ", K=" << n_plot << std::endl;
        fprintf(file_pt, "ZONE I=%i, J=%i, K=%i \n", n_plot, n_plot, n_plot);

        // Loop over element nodes
        for (unsigned l3 = 0; l3 < n_plot; l3++)
        {
          s[2] = -1.0 + l3 * 2.0 / (n_plot - 1);
          for (unsigned l2 = 0; l2 < n_plot; l2++)
          {
            s[1] = -1.0 + l2 * 2.0 / (n_plot - 1);
            for (unsigned l1 = 0; l1 < n_plot; l1++)
            {
              s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

              // Get Eulerian and Lagrangian coordinates
              this->interpolated_x(s, x);
              this->interpolated_xi(s, xi);

              // Get isotropic growth
              double gamma;
              // Dummy integration point
              unsigned ipt = 0;
              this->get_isotropic_growth(ipt, s, xi, gamma);

              // Output the x,y,..
              for (unsigned i = 0; i < DIM; i++)
              {
                // outfile << x[i] << " ";
                fprintf(file_pt, "%g ", x[i]);
              }
              // Output xi0,xi1,..
              for (unsigned i = 0; i < DIM; i++)
              {
                // outfile << xi[i] << " ";
                fprintf(file_pt, "%g ", xi[i]);
              }
              // Output growth
              // outfile << gamma << " ";
              fprintf(file_pt, "%g ", gamma);

              // Output pressure
              // outfile << interpolated_solid_p(s) << " ";
              // outfile << std::endl;
              fprintf(file_pt, "%g \n", interpolated_solid_p(s));
            }
          }
        }
        break;

      default:
        std::ostringstream error_message;
        error_message << "No output routine for PVDEquationsWithPressure<"
                      << DIM << "> elements. Write it yourself!" << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }


  //=======================================================================
  /// Output: x,y,[z],xi0,xi1,[xi2],gamma strain and stress components
  //=======================================================================
  template<unsigned DIM>
  void PVDEquationsWithPressure<DIM>::extended_output(std::ostream& outfile,
                                                      const unsigned& n_plot)
  {
    Vector<double> x(DIM);
    Vector<double> xi(DIM);
    Vector<double> s(DIM);
    DenseMatrix<double> stress_or_strain(DIM, DIM);

    // Tecplot header info
    outfile << this->tecplot_zone_string(n_plot);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(n_plot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, n_plot, s);

      // Get Eulerian and Lagrangian coordinates
      this->interpolated_x(s, x);
      this->interpolated_xi(s, xi);

      // Get isotropic growth
      double gamma;
      // Dummy integration point
      unsigned ipt = 0;
      this->get_isotropic_growth(ipt, s, xi, gamma);

      // Output the x,y,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }

      // Output xi0,xi1,..
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << xi[i] << " ";
      }

      // Output growth
      outfile << gamma << " ";

      // Output pressure
      outfile << interpolated_solid_p(s) << " ";

      // get the strain
      this->get_strain(s, stress_or_strain);
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = 0; j <= i; j++)
        {
          outfile << stress_or_strain(j, i) << " ";
        }
      }

      // get the stress
      this->get_stress(s, stress_or_strain);
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = 0; j <= i; j++)
        {
          outfile << stress_or_strain(j, i) << " ";
        }
      }


      outfile << std::endl;
    }


    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, n_plot);
    outfile << std::endl;
  }


  //=======================================================================
  /// Compute the diagonal of the velocity mass matrix for LSC
  /// preconditioner.
  //=======================================================================
  template<unsigned DIM>
  void PVDEquationsWithPressure<DIM>::get_mass_matrix_diagonal(
    Vector<double>& mass_diag)
  {
    // Resize and initialise
    mass_diag.assign(this->ndof(), 0.0);

    // find out how many nodes there are
    unsigned n_node = this->nnode();

    // Find out how many position types of dof there are
    const unsigned n_position_type = this->nnodal_position_type();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, DIM);

    // Number of integration points
    unsigned n_intpt = this->integral_pt()->nweight();

    // Integer to store the local equations no
    int local_eqn = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Call the derivatives of the shape functions
      double J = this->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);

      // Premultiply weights and Jacobian
      double W = w * J;

      // Loop over the nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the types of dof
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Loop over the directions
          for (unsigned i = 0; i < DIM; i++)
          {
            // Get the equation number
            local_eqn = this->position_local_eqn(l, k, i);

            // If not a boundary condition
            if (local_eqn >= 0)
            {
              // Add the contribution
              mass_diag[local_eqn] += pow(psi(l, k), 2) * W;
            } // End of if not boundary condition statement
          } // End of loop over dimension
        } // End of dof type
      } // End of loop over basis functions
    }
  }


  //=======================================================================
  /// Compute the contravariant second Piola Kirchoff stress at a given local
  /// coordinate. Note: this replicates a lot of code that is already
  /// coontained in get_residuals() but without sacrificing efficiency
  /// (re-computing the shape functions several times) or creating
  /// helper functions with horrendous interfaces (to pass all the
  /// functions which shouldn't be recomputed) about this is
  /// unavoidable.
  //=======================================================================
  template<unsigned DIM>
  void PVDEquationsWithPressure<DIM>::get_stress(const Vector<double>& s,
                                                 DenseMatrix<double>& sigma)
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();

    // Find out how many pressure dofs there are
    unsigned n_solid_pres = this->npres_solid();

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, DIM);

    // Set up memory for the pressure shape functions
    Shape psisp(n_solid_pres);

    // Find values of shape function
    solid_pshape(s, psisp);

    // Call the derivatives of the shape functions (ignore Jacobian)
    (void)this->dshape_lagrangian(s, psi, dpsidxi);

    // Lagrangian coordinates
    Vector<double> xi(DIM);
    this->interpolated_xi(s, xi);

    // Get isotropic growth factor
    double gamma;
    // Dummy integration point
    unsigned ipt = 0;
    this->get_isotropic_growth(ipt, s, xi, gamma);

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
      }
    }


    // Calculate interpolated values of the derivative of global position
    // wrt lagrangian coordinates
    DenseMatrix<double> interpolated_G(DIM);

    // Initialise to zero
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        interpolated_G(i, j) = 0.0;
      }
    }

    // Calculate displacements and derivatives
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over positional dofs
      for (unsigned k = 0; k < n_position_type; k++)
      {
        // Loop over displacement components (deformed position)
        for (unsigned i = 0; i < DIM; i++)
        {
          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_G(j, i) +=
              this->nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
          }
        }
      }
    }

    // Declare and calculate the deformed metric tensor
    DenseMatrix<double> G(DIM);

    // Assign values of G
    for (unsigned i = 0; i < DIM; i++)
    {
      // Do upper half of matrix
      // Note that j must be signed here for the comparison test to work
      // Also i must be cast to an int
      for (int j = (DIM - 1); j >= static_cast<int>(i); j--)
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
      for (int j = (i - 1); j >= 0; j--)
      {
        G(i, j) = G(j, i);
      }
    }


    // Calculate the interpolated solid pressure
    double interpolated_solid_p = 0.0;
    for (unsigned l = 0; l < n_solid_pres; l++)
    {
      interpolated_solid_p += solid_p(l) * psisp[l];
    }

    // Now calculate the deviatoric stress and all pressure-related
    // quantitites
    DenseMatrix<double> sigma_dev(DIM), Gup(DIM);
    double detG = 0.0;
    double gen_dil = 0.0;
    double inv_kappa = 0.0;

    // Incompressible: Compute the deviatoric part of the stress tensor, the
    // contravariant deformed metric tensor and the determinant
    // of the deformed covariant metric tensor.

    if (Incompressible)
    {
      get_stress(g, G, sigma_dev, Gup, detG);
    }
    // Nearly incompressible: Compute the deviatoric part of the
    // stress tensor, the contravariant deformed metric tensor,
    // the generalised dilatation and the inverse bulk modulus.
    else
    {
      get_stress(g, G, sigma_dev, Gup, gen_dil, inv_kappa);
    }

    // Get complete stress
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        sigma(i, j) = -interpolated_solid_p * Gup(i, j) + sigma_dev(i, j);
      }
    }
  }


  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //====================================================================
  /// Data for the number of Variables at each node
  //====================================================================
  template<>
  const unsigned QPVDElementWithContinuousPressure<2>::Initial_Nvalue[9] = {
    1, 0, 1, 0, 0, 0, 1, 0, 1};

  //==========================================================================
  /// Conversion from pressure dof to Node number at which pressure is stored
  //==========================================================================
  template<>
  const unsigned QPVDElementWithContinuousPressure<2>::Pconv[4] = {0, 2, 6, 8};

  //====================================================================
  /// Data for the number of Variables at each node
  //====================================================================
  template<>
  const unsigned QPVDElementWithContinuousPressure<3>::Initial_Nvalue[27] = {
    1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1};

  //==========================================================================
  /// Conversion from pressure dof to Node number at which pressure is stored
  //==========================================================================
  template<>
  const unsigned QPVDElementWithContinuousPressure<3>::Pconv[8] = {
    0, 2, 6, 8, 18, 20, 24, 26};


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Data for the number of variables at each node
  //=======================================================================
  template<>
  const unsigned TPVDElementWithContinuousPressure<2>::Initial_Nvalue[6] = {
    1, 1, 1, 0, 0, 0};

  //=======================================================================
  /// Data for the pressure conversion array
  //=======================================================================
  template<>
  const unsigned TPVDElementWithContinuousPressure<2>::Pconv[3] = {0, 1, 2};

  //=======================================================================
  /// Data for the number of variables at each node
  //=======================================================================
  template<>
  const unsigned TPVDElementWithContinuousPressure<3>::Initial_Nvalue[10] = {
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0};

  //=======================================================================
  /// Data for the pressure conversion array
  //=======================================================================
  template<>
  const unsigned TPVDElementWithContinuousPressure<3>::Pconv[4] = {0, 1, 2, 3};


  // Instantiate the required elements
  template class QPVDElementWithPressure<2>;
  template class QPVDElementWithContinuousPressure<2>;
  template class PVDEquationsBase<2>;
  template class PVDEquations<2>;
  template class PVDEquationsWithPressure<2>;

  template class QPVDElementWithPressure<3>;
  template class QPVDElementWithContinuousPressure<3>;
  template class PVDEquationsBase<3>;
  template class PVDEquations<3>;
  template class PVDEquationsWithPressure<3>;

  template class TPVDElementWithContinuousPressure<2>;
  template class TPVDElementWithContinuousPressure<3>;

} // namespace oomph
