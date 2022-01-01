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
// Non-inline member functions and static member data for refineable solid
// mechanics elements

#include "refineable_solid_elements.h"

namespace oomph
{
  //====================================================================
  /// Residuals for Refineable QPVDElements
  //====================================================================
  template<unsigned DIM>
  void RefineablePVDEquations<DIM>::
    fill_in_generic_contribution_to_residuals_pvd(Vector<double>& residuals,
                                                  DenseMatrix<double>& jacobian,
                                                  const unsigned& flag)
  {
#ifdef PARANOID
    // Check if the constitutive equation requires the explicit imposition of an
    // incompressibility constraint
    if (this->Constitutive_law_pt->requires_incompressibility_constraint())
    {
      throw OomphLibError("RefineablePVDEquations cannot be used with "
                          "incompressible constitutive laws.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Simply set up initial condition?
    if (Solid_ic_pt != 0)
    {
      get_residuals_for_solid_ic(residuals);
      return;
    }

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Find out how many positional dofs there are
    const unsigned n_position_type = this->nnodal_position_type();

    // Integers to store local equation numbers
    int local_eqn = 0;

    // Timescale ratio (non-dim density)
    double lambda_sq = this->lambda_sq();

    // Time factor
    double time_factor = 0.0;
    if (lambda_sq > 0)
    {
      time_factor = this->node_pt(0)->position_time_stepper_pt()->weight(2, 0);
    }

    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, DIM);

    // Set the value of n_intpt -- the number of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Set the vector to hold the local coordinates in the element
    Vector<double> s(DIM);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign the values of s
      for (unsigned i = 0; i < DIM; ++i)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape functions (and get Jacobian)
      double J = dshape_lagrangian_at_knot(ipt, psi, dpsidxi);

      // Calculate interpolated values of the derivative of global position
      // wrt lagrangian coordinates
      DenseMatrix<double> interpolated_G(DIM, DIM, 0.0);

      // Setup memory for accelerations
      Vector<double> accel(DIM, 0.0);

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
            interpolated_xi[i] += lagrangian_position_gen(l, k, i) * psi_;

            // Only compute accelerations if inertia is switched on
            // otherwise the timestepper might not be able to
            // work out dx_gen_dt(2,...)
            if ((lambda_sq > 0.0) && (this->Unsteady))
            {
              accel[i] += dnodal_position_gen_dt(2, l, k, i) * psi_;
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
      this->get_stress(g, G, sigma);

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

        // Get the "upper triangular"
        // entries of the derivatives of the stress tensor with
        // respect to G
        this->get_d_stress_dG_upper(g, G, sigma, d_stress_dG);
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


      // Default setting for non-hanging node
      unsigned n_master = 1;
      double hang_weight = 1.0;

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get pointer to local node l
        Node* local_node_pt = node_pt(l);

        // Cache hang status
        bool is_hanging = local_node_pt->is_hanging();

        // If the node is a hanging node
        if (is_hanging)
        {
          n_master = local_node_pt->hanging_pt()->nmaster();
        }
        // Otherwise the node is its own master
        else
        {
          n_master = 1;
        }


        // Storage for local equation numbers at node indexed by
        // type and direction
        DenseMatrix<int> position_local_eqn_at_node(n_position_type, DIM);

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          if (is_hanging)
          {
            // Find the equation numbers
            position_local_eqn_at_node = local_position_hang_eqn(
              local_node_pt->hanging_pt()->master_node_pt(m));

            // Find the hanging node weight
            hang_weight = local_node_pt->hanging_pt()->master_weight(m);
          }
          else
          {
            // Loop of types of dofs
            for (unsigned k = 0; k < n_position_type; k++)
            {
              // Loop over the displacement components
              for (unsigned i = 0; i < DIM; i++)
              {
                position_local_eqn_at_node(k, i) = position_local_eqn(l, k, i);
              }
            }

            // Hang weight is one
            hang_weight = 1.0;
          }

          // Loop of types of dofs
          for (unsigned k = 0; k < n_position_type; k++)
          {
            // Offset for faster access
            const unsigned offset5 = dpsidxi.offset(l, k);

            // Loop over the displacement components
            for (unsigned i = 0; i < DIM; i++)
            {
              local_eqn = position_local_eqn_at_node(k, i);

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
                residuals[local_eqn] += W * sum * hang_weight;


                // Get Jacobian too?
                if (flag == 1)
                {
                  // Offset for faster access in general stress loop
                  const unsigned offset1 = d_G_dX.offset(l, k, i);

                  // Default setting for non-hanging node
                  unsigned nn_master = 1;
                  double hhang_weight = 1.0;

                  // Loop over the nodes of the element again
                  for (unsigned ll = 0; ll < n_node; ll++)
                  {
                    // Get pointer to local node ll
                    Node* llocal_node_pt = node_pt(ll);

                    // Cache hang status
                    bool iis_hanging = llocal_node_pt->is_hanging();

                    // If the node is a hanging node
                    if (iis_hanging)
                    {
                      nn_master = llocal_node_pt->hanging_pt()->nmaster();
                    }
                    // Otherwise the node is its own master
                    else
                    {
                      nn_master = 1;
                    }


                    // Storage for local unknown numbers at node indexed by
                    // type and direction
                    DenseMatrix<int> position_local_unk_at_node(n_position_type,
                                                                DIM);

                    // Loop over the master nodes
                    for (unsigned mm = 0; mm < nn_master; mm++)
                    {
                      if (iis_hanging)
                      {
                        // Find the unknown numbers
                        position_local_unk_at_node = local_position_hang_eqn(
                          llocal_node_pt->hanging_pt()->master_node_pt(mm));

                        // Find the hanging node weight
                        hhang_weight =
                          llocal_node_pt->hanging_pt()->master_weight(mm);
                      }
                      else
                      {
                        // Loop of types of dofs
                        for (unsigned kk = 0; kk < n_position_type; kk++)
                        {
                          // Loop over the displacement components
                          for (unsigned ii = 0; ii < DIM; ii++)
                          {
                            position_local_unk_at_node(kk, ii) =
                              position_local_eqn(ll, kk, ii);
                          }
                        }

                        // Hang weight is one
                        hhang_weight = 1.0;
                      }


                      // Loop of types of dofs again
                      for (unsigned kk = 0; kk < n_position_type; kk++)
                      {
                        // Loop over the displacement components again
                        for (unsigned ii = 0; ii < DIM; ii++)
                        {
                          // Get the number of the unknown
                          int local_unknown =
                            position_local_unk_at_node(kk, ii);


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
                                double factor =
                                  d_G_dX.raw_direct_access(count1);
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

                                  // Only upper half of derivatives w.r.t.
                                  // symm tensor
                                  for (unsigned bb = aa; bb < DIM; bb++)
                                  {
                                    sum +=
                                      factor *
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
                            jacobian(local_eqn, local_unknown) +=
                              sum * W * hang_weight * hhang_weight;

                            // Only upper triangle (no separate test for bc as
                            // local_eqn is already nonnegative)
                            if ((i == ii) && (local_unknown >= local_eqn))
                            {
                              // Initialise the contribution
                              double sum = 0.0;

                              // Inertia term
                              sum += lambda_sq * time_factor * psi(ll, kk) *
                                     psi(l, k);

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

                              // Multiply by weights to form contribution
                              double sym_entry =
                                sum * W * hang_weight * hhang_weight;
                              // Add contribution to jacobian
                              jacobian(local_eqn, local_unknown) += sym_entry;
                              // Add to lower triangular entries
                              if (local_eqn != local_unknown)
                              {
                                jacobian(local_unknown, local_eqn) += sym_entry;
                              }
                            }
                          } // End of if not boundary condition
                        }
                      }
                    }
                  }
                }

              } // End of if not boundary condition

            } // End of loop over coordinate directions
          } // End of loop over type of dof
        } // End of loop over master nodes
      } // End of loop over nodes
    } // End of loop over integration points
  }


  //=======================================================================
  /// Compute the diagonal of the velocity mass matrix for LSC
  /// preconditioner.
  //=======================================================================
  template<unsigned DIM>
  void RefineablePVDEquationsWithPressure<DIM>::get_mass_matrix_diagonal(
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

      unsigned n_master = 1;
      double hang_weight = 1.0;

      // Loop over the nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get pointer to local node l
        Node* local_node_pt = node_pt(l);

        // Cache hang status
        bool is_hanging = local_node_pt->is_hanging();

        // If the node is a hanging node
        if (is_hanging)
        {
          n_master = local_node_pt->hanging_pt()->nmaster();
        }
        // Otherwise the node is its own master
        else
        {
          n_master = 1;
        }

        // Storage for local equation numbers at node indexed by
        // type and direction
        DenseMatrix<int> position_local_eqn_at_node(n_position_type, DIM);

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          if (is_hanging)
          {
            // Find the equation numbers
            position_local_eqn_at_node = local_position_hang_eqn(
              local_node_pt->hanging_pt()->master_node_pt(m));

            // Find the hanging node weight
            hang_weight = local_node_pt->hanging_pt()->master_weight(m);
          }
          else
          {
            // Loop of types of dofs
            for (unsigned k = 0; k < n_position_type; k++)
            {
              // Loop over the displacement components
              for (unsigned i = 0; i < DIM; i++)
              {
                position_local_eqn_at_node(k, i) = position_local_eqn(l, k, i);
              }
            }

            // Hang weight is one
            hang_weight = 1.0;
          }

          // Loop over the types of dof
          for (unsigned k = 0; k < n_position_type; k++)
          {
            // Loop over the directions
            for (unsigned i = 0; i < DIM; i++)
            {
              // Get the equation number
              local_eqn = position_local_eqn_at_node(k, i);

              // If not a boundary condition
              if (local_eqn >= 0)
              {
                // Add the contribution
                mass_diag[local_eqn] += pow(psi(l, k) * hang_weight, 2) * W;
              } // End of if not boundary condition statement
            } // End of loop over dimension
          } // End of dof type
        } // End of loop over master nodes
      } // End of loop over basis functions (nodes)
    } // End integration loop
  }


  //===========================================================================
  /// Fill in element's contribution to the elemental
  /// residual vector and/or Jacobian matrix.
  /// flag=0: compute only residual vector
  /// flag=1: compute both, fully analytically
  /// flag=2: compute both, using FD for the derivatives w.r.t. to the
  ///         discrete displacment dofs.
  /// flag=3: compute residuals, jacobian (full analytic) and mass matrix
  /// flag=4: compute residuals, jacobian (FD for derivatives w.r.t.
  ///          displacements) and mass matrix
  //==========================================================================
  template<unsigned DIM>
  void RefineablePVDEquationsWithPressure<DIM>::
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
        (!this->Incompressible))
    {
      throw OomphLibError("The constitutive law requires the use of the "
                          "incompressible formulation by setting the element's "
                          "member function set_incompressible()",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Simply set up initial condition?
    if (Solid_ic_pt != 0)
    {
      get_residuals_for_solid_ic(residuals);
      return;
    }

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Find out how many position types of dof there are
    const unsigned n_position_type = this->nnodal_position_type();

    // Find out how many pressure dofs there are
    const unsigned n_solid_pres = this->npres_solid();

    // Find out the index of the solid dof
    const int solid_p_index = this->solid_p_nodal_index();

    // Local array of booleans that is true if the l-th pressure value is
    // hanging This is an optimization because it avoids repeated virtual
    // function calls
    bool solid_pressure_dof_is_hanging[n_solid_pres];

    // If the solid pressure is stored at a node
    if (solid_p_index >= 0)
    {
      // Read out whether the solid pressure is hanging
      for (unsigned l = 0; l < n_solid_pres; ++l)
      {
        solid_pressure_dof_is_hanging[l] =
          solid_pressure_node_pt(l)->is_hanging(solid_p_index);
      }
    }
    // Otherwise the pressure is not stored at a node and so
    // it cannot hang
    else
    {
      for (unsigned l = 0; l < n_solid_pres; ++l)
      {
        solid_pressure_dof_is_hanging[l] = false;
      }
    }

    // Integer for storage of local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0;

    // Timescale ratio (non-dim density)
    double lambda_sq = this->lambda_sq();


    // Time factor
    double time_factor = 0.0;
    if (lambda_sq > 0)
    {
      time_factor = this->node_pt(0)->position_time_stepper_pt()->weight(2, 0);
    }


    // Set up memory for the shape functions
    Shape psi(n_node, n_position_type);
    DShape dpsidxi(n_node, n_position_type, DIM);

    // Set up memory for the pressure shape functions
    Shape psisp(n_solid_pres);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Set the vector to hold the local coordinates in the element
    Vector<double> s(DIM);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign the values of s
      for (unsigned i = 0; i < DIM; ++i)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape functions
      double J = dshape_lagrangian_at_knot(ipt, psi, dpsidxi);

      // Call the pressure shape functions
      this->solid_pshape_at_knot(ipt, psisp);

      // Storage for Lagrangian coordinates (initialised to zero)
      Vector<double> interpolated_xi(DIM, 0.0);

      // Deformed tangent vectors
      DenseMatrix<double> interpolated_G(DIM, DIM, 0.0);

      // Setup memory for accelerations
      Vector<double> accel(DIM, 0.0);

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
            // Calculate the lagrangian coordinates and the accelerations
            interpolated_xi[i] += lagrangian_position_gen(l, k, i) * psi_;

            // Only compute accelerations if inertia is switched on
            // otherwise the timestepper might not be able to
            // work out dx_gen_dt(2,...)
            if ((lambda_sq > 0.0) && (this->Unsteady))
            {
              accel[i] += dnodal_position_gen_dt(2, l, k, i) * psi_;
            }

            // Loop over derivative directions
            for (unsigned j = 0; j < DIM; j++)
            {
              interpolated_G(j, i) +=
                nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
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
        interpolated_solid_p += this->solid_p(l) * psisp[l];
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
      // RankFourTensor<double> d_Gup_dG(DIM,DIM,DIM,DIM,0.0);
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
      if (this->Incompressible)
      {
        this->get_stress(g, G, sigma_dev, Gup, detG);

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
          // Get the "upper triangular" entries of the
          // derivatives of the stress tensor with
          // respect to G
          this->get_d_stress_dG_upper(
            g, G, sigma, detG, interpolated_solid_p, d_stress_dG, d_detG_dG);
        }
      }
      // Nearly incompressible: Compute the deviatoric part of the
      // stress tensor, the contravariant deformed metric tensor,
      // the generalised dilatation and the inverse bulk modulus.
      else
      {
        this->get_stress(g, G, sigma_dev, Gup, gen_dil, inv_kappa);

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
          // Get the "upper triangular" entries of the derivatives
          // of the stress tensor with
          // respect to G
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

      unsigned n_master = 1;
      double hang_weight = 1.0;

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get pointer to local node l
        Node* local_node_pt = node_pt(l);

        // Cache hang status
        bool is_hanging = local_node_pt->is_hanging();

        // If the node is a hanging node
        if (is_hanging)
        {
          n_master = local_node_pt->hanging_pt()->nmaster();
        }
        // Otherwise the node is its own master
        else
        {
          n_master = 1;
        }


        // Storage for local equation numbers at node indexed by
        // type and direction
        DenseMatrix<int> position_local_eqn_at_node(n_position_type, DIM);

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          if (is_hanging)
          {
            // Find the equation numbers
            position_local_eqn_at_node = local_position_hang_eqn(
              local_node_pt->hanging_pt()->master_node_pt(m));

            // Find the hanging node weight
            hang_weight = local_node_pt->hanging_pt()->master_weight(m);
          }
          else
          {
            // Loop of types of dofs
            for (unsigned k = 0; k < n_position_type; k++)
            {
              // Loop over the displacement components
              for (unsigned i = 0; i < DIM; i++)
              {
                position_local_eqn_at_node(k, i) = position_local_eqn(l, k, i);
              }
            }

            // Hang weight is one
            hang_weight = 1.0;
          }


          // Loop of types of dofs
          for (unsigned k = 0; k < n_position_type; k++)
          {
            // Offset for faster access
            const unsigned offset5 = dpsidxi.offset(l, k);

            // Loop over the displacement components
            for (unsigned i = 0; i < DIM; i++)
            {
              local_eqn = position_local_eqn_at_node(k, i);

              /*IF it's not a boundary condition*/
              if (local_eqn >= 0)
              {
                // Initialise contribution to zero
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
                residuals[local_eqn] += W * sum * hang_weight;


                // Get the mass matrix
                // This involves another loop over the points
                // because the jacobian may NOT be being calculated analytically
                // It could be made more efficient in th event that
                // we eventually decide not (never) to
                // use finite differences.
                if (flag > 2)
                {
                  // Default setting for non-hanging node
                  unsigned nn_master = 1;
                  double hhang_weight = 1.0;

                  // Loop over the nodes of the element again
                  for (unsigned ll = 0; ll < n_node; ll++)
                  {
                    // Get pointer to local node ll
                    Node* llocal_node_pt = node_pt(ll);

                    // Cache hang status
                    bool iis_hanging = llocal_node_pt->is_hanging();

                    // If the node is a hanging node
                    if (iis_hanging)
                    {
                      nn_master = llocal_node_pt->hanging_pt()->nmaster();
                    }
                    // Otherwise the node is its own master
                    else
                    {
                      nn_master = 1;
                    }


                    // Storage for local unknown numbers at node indexed by
                    // type and direction
                    DenseMatrix<int> position_local_unk_at_node(n_position_type,
                                                                DIM);

                    // Loop over the master nodes
                    for (unsigned mm = 0; mm < nn_master; mm++)
                    {
                      if (iis_hanging)
                      {
                        // Find the unknown numbers
                        position_local_unk_at_node = local_position_hang_eqn(
                          llocal_node_pt->hanging_pt()->master_node_pt(mm));

                        // Find the hanging node weight
                        hhang_weight =
                          llocal_node_pt->hanging_pt()->master_weight(mm);
                      }
                      else
                      {
                        // Loop of types of dofs
                        for (unsigned kk = 0; kk < n_position_type; kk++)
                        {
                          // Loop over the displacement components
                          for (unsigned ii = 0; ii < DIM; ii++)
                          {
                            position_local_unk_at_node(kk, ii) =
                              position_local_eqn(ll, kk, ii);
                          }
                        }

                        // Hang weight is one
                        hhang_weight = 1.0;
                      }


                      // Loop of types of dofs again
                      for (unsigned kk = 0; kk < n_position_type; kk++)
                      {
                        // Get the number of the unknown
                        int local_unknown = position_local_unk_at_node(kk, i);

                        /*IF it's not a boundary condition*/
                        if (local_unknown >= 0)
                        {
                          mass_matrix(local_eqn, local_unknown) +=
                            lambda_sq * psi(l, k) * psi(ll, kk) * hang_weight *
                            hhang_weight * W;
                        }
                      }
                    }
                  }
                }


                // Get Jacobian too?
                if ((flag == 1) || (flag == 3))
                {
                  // Offset for faster access in general stress loop
                  const unsigned offset1 = d_G_dX.offset(l, k, i);

                  // Default setting for non-hanging node
                  unsigned nn_master = 1;
                  double hhang_weight = 1.0;

                  // Loop over the nodes of the element again
                  for (unsigned ll = 0; ll < n_node; ll++)
                  {
                    // Get pointer to local node ll
                    Node* llocal_node_pt = node_pt(ll);

                    // Cache hang status
                    bool iis_hanging = llocal_node_pt->is_hanging();

                    // If the node is a hanging node
                    if (iis_hanging)
                    {
                      nn_master = llocal_node_pt->hanging_pt()->nmaster();
                    }
                    // Otherwise the node is its own master
                    else
                    {
                      nn_master = 1;
                    }


                    // Storage for local unknown numbers at node indexed by
                    // type and direction
                    DenseMatrix<int> position_local_unk_at_node(n_position_type,
                                                                DIM);

                    // Loop over the master nodes
                    for (unsigned mm = 0; mm < nn_master; mm++)
                    {
                      if (iis_hanging)
                      {
                        // Find the unknown numbers
                        position_local_unk_at_node = local_position_hang_eqn(
                          llocal_node_pt->hanging_pt()->master_node_pt(mm));

                        // Find the hanging node weight
                        hhang_weight =
                          llocal_node_pt->hanging_pt()->master_weight(mm);
                      }
                      else
                      {
                        // Loop of types of dofs
                        for (unsigned kk = 0; kk < n_position_type; kk++)
                        {
                          // Loop over the displacement components
                          for (unsigned ii = 0; ii < DIM; ii++)
                          {
                            position_local_unk_at_node(kk, ii) =
                              position_local_eqn(ll, kk, ii);
                          }
                        }

                        // Hang weight is one
                        hhang_weight = 1.0;
                      }


                      // Loop of types of dofs again
                      for (unsigned kk = 0; kk < n_position_type; kk++)
                      {
                        // Loop over the displacement components again
                        for (unsigned ii = 0; ii < DIM; ii++)
                        {
                          // Get the number of the unknown
                          int local_unknown =
                            position_local_unk_at_node(kk, ii);

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
                                double factor =
                                  d_G_dX.raw_direct_access(count1);
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

                                  // Only upper half of derivatives w.r.t.
                                  // symm tensor
                                  for (unsigned bb = aa; bb < DIM; bb++)
                                  {
                                    sum +=
                                      factor *
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
                            jacobian(local_eqn, local_unknown) +=
                              sum * W * hang_weight * hhang_weight;

                            // Only upper triangle (no separate test for bc as
                            // local_eqn is already nonnegative)
                            if ((i == ii) && (local_unknown >= local_eqn))
                            {
                              // Initialise contribution
                              double sum = 0.0;

                              // Inertia term
                              sum += lambda_sq * time_factor * psi(ll, kk) *
                                     psi(l, k);

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

                              // Multiply by weights to form contribution
                              double sym_entry =
                                sum * W * hang_weight * hhang_weight;
                              // Add contribution to jacobian
                              jacobian(local_eqn, local_unknown) += sym_entry;
                              // Add to lower triangular entries
                              if (local_eqn != local_unknown)
                              {
                                jacobian(local_unknown, local_eqn) += sym_entry;
                              }
                            }
                          } // End of if not boundary condition
                        }
                      }
                    }
                  }
                }

                // Can add in the pressure jacobian terms
                if (flag > 0)
                {
                  // Loop over the pressure nodes
                  for (unsigned l2 = 0; l2 < n_solid_pres; l2++)
                  {
                    unsigned n_master2 = 1;
                    double hang_weight2 = 1.0;
                    HangInfo* hang_info2_pt = 0;

                    bool is_hanging2 = solid_pressure_dof_is_hanging[l2];
                    if (is_hanging2)
                    {
                      // Get the HangInfo object associated with the
                      // hanging solid pressure
                      hang_info2_pt =
                        solid_pressure_node_pt(l2)->hanging_pt(solid_p_index);

                      n_master2 = hang_info2_pt->nmaster();
                    }
                    else
                    {
                      n_master2 = 1;
                    }

                    // Loop over all the master nodes
                    for (unsigned m2 = 0; m2 < n_master2; m2++)
                    {
                      if (is_hanging2)
                      {
                        // Get the equation numbers at the master node
                        local_unknown = local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), solid_p_index);

                        // Find the hanging node weight at the node
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      else
                      {
                        local_unknown = this->solid_p_local_eqn(l2);
                        hang_weight2 = 1.0;
                      }

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
                              dpsidxi(l, k, b) * W * hang_weight * hang_weight2;
                          }
                        }
                      }
                    } // End of loop over master nodes
                  } // End of loop over pressure dofs
                } // End of Jacobian terms

              } // End of if not boundary condition
            }
          }
        } // End of loop of over master nodes

      } // End of loop over nodes

      //==============CONSTRAINT EQUATIONS FOR PRESSURE=====================

      // Now loop over the pressure degrees of freedom
      for (unsigned l = 0; l < n_solid_pres; l++)
      {
        bool is_hanging = solid_pressure_dof_is_hanging[l];

        unsigned n_master = 1;
        double hang_weight = 1.0;
        HangInfo* hang_info_pt = 0;

        // If the node is a hanging node
        if (is_hanging)
        {
          // Get a pointer to the HangInfo object associated with the
          // solid pressure (stored at solid_p_index)
          hang_info_pt = solid_pressure_node_pt(l)->hanging_pt(solid_p_index);

          // Number of master nodes
          n_master = hang_info_pt->nmaster();
        }
        // Otherwise the node is its own master
        else
        {
          n_master = 1;
        }

        // Loop over all the master nodes
        // Note that the pressure is stored at the inded solid_p_index
        for (unsigned m = 0; m < n_master; m++)
        {
          if (is_hanging)
          {
            // Get the equation numbers at the master node
            local_eqn =
              local_hang_eqn(hang_info_pt->master_node_pt(m), solid_p_index);

            // Find the hanging node weight at the node
            hang_weight = hang_info_pt->master_weight(m);
          }
          else
          {
            local_eqn = this->solid_p_local_eqn(l);
          }

          // Pinned (unlikely, actually) or real dof?
          if (local_eqn >= 0)
          {
            // For true incompressibility we need to conserve volume
            // so the determinant of the deformed metric tensor
            // needs to be equal to that of the undeformed one, which
            // is equal to the volumetric growth factor
            if (this->Incompressible)
            {
              residuals[local_eqn] +=
                (detG - gamma) * psisp[l] * W * hang_weight;

              // Get Jacobian too?
              if ((flag == 1) || (flag == 3))
              {
                // Default setting for non-hanging node
                unsigned nn_master = 1;
                double hhang_weight = 1.0;

                // Loop over the nodes of the element again
                for (unsigned ll = 0; ll < n_node; ll++)
                {
                  // Get pointer to local node ll
                  Node* llocal_node_pt = node_pt(ll);

                  // Cache hang status
                  bool iis_hanging = llocal_node_pt->is_hanging();

                  // If the node is a hanging node
                  if (iis_hanging)
                  {
                    nn_master = llocal_node_pt->hanging_pt()->nmaster();
                  }
                  // Otherwise the node is its own master
                  else
                  {
                    nn_master = 1;
                  }

                  // Storage for local unknown numbers at node indexed by
                  // type and direction
                  DenseMatrix<int> position_local_unk_at_node(n_position_type,
                                                              DIM);

                  // Loop over the master nodes
                  for (unsigned mm = 0; mm < nn_master; mm++)
                  {
                    if (iis_hanging)
                    {
                      // Find the unknown numbers
                      position_local_unk_at_node = local_position_hang_eqn(
                        llocal_node_pt->hanging_pt()->master_node_pt(mm));

                      // Find the hanging node weight
                      hhang_weight =
                        llocal_node_pt->hanging_pt()->master_weight(mm);
                    }
                    else
                    {
                      // Loop of types of dofs
                      for (unsigned kk = 0; kk < n_position_type; kk++)
                      {
                        // Loop over the displacement components
                        for (unsigned ii = 0; ii < DIM; ii++)
                        {
                          position_local_unk_at_node(kk, ii) =
                            position_local_eqn(ll, kk, ii);
                        }
                      }

                      // Hang weight is one
                      hhang_weight = 1.0;
                    }


                    // Loop of types of dofs again
                    for (unsigned kk = 0; kk < n_position_type; kk++)
                    {
                      // Loop over the displacement components again
                      for (unsigned ii = 0; ii < DIM; ii++)
                      {
                        // Get the number of the unknown
                        int local_unknown = position_local_unk_at_node(kk, ii);

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
                          jacobian(local_eqn, local_unknown) +=
                            sum * W * hang_weight * hhang_weight;
                        }
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
                (inv_kappa * interpolated_solid_p + gen_dil) * psisp[l] * W *
                hang_weight;

              // Add in the jacobian terms
              if ((flag == 1) || (flag == 3))
              {
                // Default setting for non-hanging node
                unsigned nn_master = 1;
                double hhang_weight = 1.0;

                // Loop over the nodes of the element again
                for (unsigned ll = 0; ll < n_node; ll++)
                {
                  // Get pointer to local node ll
                  Node* llocal_node_pt = node_pt(ll);

                  // Cache hang status
                  bool iis_hanging = llocal_node_pt->is_hanging();

                  // If the node is a hanging node
                  if (iis_hanging)
                  {
                    nn_master = llocal_node_pt->hanging_pt()->nmaster();
                  }
                  // Otherwise the node is its own master
                  else
                  {
                    nn_master = 1;
                  }


                  // Storage for local unknown numbers at node indexed by
                  // type and direction
                  DenseMatrix<int> position_local_unk_at_node(n_position_type,
                                                              DIM);

                  // Loop over the master nodes
                  for (unsigned mm = 0; mm < nn_master; mm++)
                  {
                    if (iis_hanging)
                    {
                      // Find the unknown numbers
                      position_local_unk_at_node = local_position_hang_eqn(
                        llocal_node_pt->hanging_pt()->master_node_pt(mm));

                      // Find the hanging node weight
                      hhang_weight =
                        llocal_node_pt->hanging_pt()->master_weight(mm);
                    }
                    else
                    {
                      // Loop of types of dofs
                      for (unsigned kk = 0; kk < n_position_type; kk++)
                      {
                        // Loop over the displacement components
                        for (unsigned ii = 0; ii < DIM; ii++)
                        {
                          position_local_unk_at_node(kk, ii) =
                            position_local_eqn(ll, kk, ii);
                        }
                      }

                      // Hang weight is one
                      hhang_weight = 1.0;
                    }


                    // Loop of types of dofs again
                    for (unsigned kk = 0; kk < n_position_type; kk++)
                    {
                      // Loop over the displacement components again
                      for (unsigned ii = 0; ii < DIM; ii++)
                      {
                        // Get the number of the unknown
                        int local_unknown = position_local_unk_at_node(kk, ii);

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
                          jacobian(local_eqn, local_unknown) +=
                            sum * W * hang_weight * hhang_weight;
                        }
                      }
                    }
                  }
                }
              }


              // Add in the pressure jacobian terms
              if (flag > 0)
              {
                // Loop over the pressure nodes again
                for (unsigned l2 = 0; l2 < n_solid_pres; l2++)
                {
                  bool is_hanging2 = solid_pressure_dof_is_hanging[l2];

                  unsigned n_master2 = 1;
                  double hang_weight2 = 1.0;
                  HangInfo* hang_info2_pt = 0;

                  if (is_hanging2)
                  {
                    // Get pointer to hang info object
                    // Note that the pressure is stored at
                    // the index solid_p_index
                    hang_info2_pt =
                      solid_pressure_node_pt(l2)->hanging_pt(solid_p_index);

                    n_master2 = hang_info2_pt->nmaster();
                  }
                  else
                  {
                    n_master2 = 1;
                  }

                  // Loop over all the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    if (is_hanging2)
                    {
                      // Get the equation numbers at the master node
                      local_unknown = local_hang_eqn(
                        hang_info2_pt->master_node_pt(m2), solid_p_index);

                      // Find the hanging node weight at the node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    else
                    {
                      local_unknown = this->solid_p_local_eqn(l2);
                      hang_weight2 = 1.0;
                    }

                    // If it's not a boundary condition
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) +=
                        inv_kappa * psisp[l2] * psisp[l] * W * hang_weight *
                        hang_weight2;
                    }

                  } // End of loop over master nodes
                } // End of loop over pressure dofs
              } // End of pressure Jacobian


            } // End of nearly incompressible case
          } // End of if not boundary condition
        } // End of loop over master nodes
      } // End of loop over pressure dofs
    } // End of loop over integration points
  }


  //====================================================================
  /// Forcing building of required templates
  //====================================================================
  template class RefineablePVDEquations<2>;
  template class RefineablePVDEquations<3>;

  template class RefineablePVDEquationsWithPressure<2>;
  template class RefineablePVDEquationsWithPressure<3>;

} // namespace oomph
