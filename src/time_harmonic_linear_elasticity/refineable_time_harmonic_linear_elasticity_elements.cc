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
// Non-inline member functions and static member data for refineable
// linear elasticity elements


#include "refineable_time_harmonic_linear_elasticity_elements.h"

namespace oomph
{
  //====================================================================
  /// Residuals for Refineable QTimeHarmonicLinearElasticityElements
  //====================================================================
  template<unsigned DIM>
  void RefineableTimeHarmonicLinearElasticityEquations<DIM>::
    fill_in_generic_contribution_to_residuals_time_harmonic_linear_elasticity(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
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
      throw OomphLibError("No elasticity tensor set",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

#endif

    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Find the indices at which the local displacements are stored
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

    // Local storage for pointers to hang_info objects
    HangInfo *hang_info_pt = 0, *hang_info2_pt = 0;

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

      // Calculate interpolated values of the derivative of displacements
      DenseMatrix<std::complex<double>> interpolated_dudx(
        DIM, DIM, std::complex<double>(0.0, 0.0));

      // Calculate displacements and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over displacement components
        for (unsigned i = 0; i < DIM; i++)
        {
          // Calculate the coordinates and the accelerations
          interpolated_x[i] += this->nodal_position(l, i) * psi(l);

          // Get the nodal displacements
          const std::complex<double> u_value =
            std::complex<double>(nodal_value(l, u_nodal_index[i].real()),
                                 nodal_value(l, u_nodal_index[i].imag()));

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

      // Number of master nodes and storage for the weight of the shape function
      unsigned n_master = 1;
      double hang_weight = 1.0;

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Local boolean to indicate whether the node is hanging
        bool is_node_hanging = node_pt(l)->is_hanging();

        // If the node is hanging
        if (is_node_hanging)
        {
          hang_info_pt = node_pt(l)->hanging_pt();
          // Read out number of master nodes from hanging data
          n_master = hang_info_pt->nmaster();
        }
        // Otherwise the node is its own master
        else
        {
          n_master = 1;
        }

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Loop over the displacement components
          for (unsigned a = 0; a < DIM; a++)
          {
            // Real
            //-----

            // Get the equation number
            if (is_node_hanging)
            {
              // Get the equation number from the master node
              local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                               u_nodal_index[a].real());
              // Get the hang weight from the master node
              hang_weight = hang_info_pt->master_weight(m);
            }
            // Otherwise the node is not hanging
            else
            {
              local_eqn = this->nodal_local_eqn(l, u_nodal_index[a].real());
              hang_weight = 1.0;
            }

            /*IF it's not a boundary condition*/
            if (local_eqn >= 0)
            {
              // Acceleration and body force
              residuals[local_eqn] +=
                (-omega_sq_local * interpolated_u[a].real() - b[a].real()) *
                psi(l) * W * hang_weight;

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
                                            dpsidx(l, b) * W * hang_weight;
                  }
                }
              }

              // Jacobian entries
              if (flag)
              {
                // Number of master nodes and weights
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                // Loop over the displacement basis functions again
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Local boolean to indicate whether the node is hanging
                  bool is_node2_hanging = node_pt(l2)->is_hanging();

                  // If the node is hanging
                  if (is_node2_hanging)
                  {
                    hang_info2_pt = node_pt(l2)->hanging_pt();
                    // Read out number of master nodes from hanging data
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise the node is its own master
                  else
                  {
                    n_master2 = 1;
                  }

                  // Loop over the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    // Loop over the displacement components again
                    for (unsigned c = 0; c < DIM; c++)
                    {
                      // Get the number of the unknown
                      // If the node is hanging
                      if (is_node2_hanging)
                      {
                        // Get the equation number from the master node
                        local_unknown = this->local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2),
                          u_nodal_index[c].real());
                        // Get the hang weights
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      else
                      {
                        local_unknown =
                          this->nodal_local_eqn(l2, u_nodal_index[c].real());
                        hang_weight2 = 1.0;
                      }

                      // If it's not pinned
                      if (local_unknown >= 0)
                      {
                        // Inertial term
                        if (a == c)
                        {
                          jacobian(local_eqn, local_unknown) -=
                            omega_sq_local * psi(l) * psi(l2) * W *
                            hang_weight * hang_weight2;
                        }

                        // Stress term
                        for (unsigned b = 0; b < DIM; b++)
                        {
                          for (unsigned d = 0; d < DIM; d++)
                          {
                            // Add the contribution to the Jacobian matrix
                            jacobian(local_eqn, local_unknown) +=
                              this->E(a, b, c, d) * dpsidx(l2, d) *
                              dpsidx(l, b) * W * hang_weight * hang_weight2;
                          }
                        }
                      } // End of if not boundary condition
                    }
                  }
                }
              } // End of jacobian calculation

            } // End of if not boundary condition


            // Imag
            //-----

            // Get the equation number
            if (is_node_hanging)
            {
              // Get the equation number from the master node
              local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                               u_nodal_index[a].imag());
              // Get the hang weight from the master node
              hang_weight = hang_info_pt->master_weight(m);
            }
            // Otherwise the node is not hanging
            else
            {
              local_eqn = this->nodal_local_eqn(l, u_nodal_index[a].imag());
              hang_weight = 1.0;
            }

            /*IF it's not a boundary condition*/
            if (local_eqn >= 0)
            {
              // Acceleration and body force
              residuals[local_eqn] +=
                (-omega_sq_local * interpolated_u[a].imag() - b[a].imag()) *
                psi(l) * W * hang_weight;

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
                                            dpsidx(l, b) * W * hang_weight;
                  }
                }
              }

              // Jacobian entries
              if (flag)
              {
                // Number of master nodes and weights
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                // Loop over the displacement basis functions again
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Local boolean to indicate whether the node is hanging
                  bool is_node2_hanging = node_pt(l2)->is_hanging();

                  // If the node is hanging
                  if (is_node2_hanging)
                  {
                    hang_info2_pt = node_pt(l2)->hanging_pt();
                    // Read out number of master nodes from hanging data
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise the node is its own master
                  else
                  {
                    n_master2 = 1;
                  }

                  // Loop over the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    // Loop over the displacement components again
                    for (unsigned c = 0; c < DIM; c++)
                    {
                      // Get the number of the unknown
                      // If the node is hanging
                      if (is_node2_hanging)
                      {
                        // Get the equation number from the master node
                        local_unknown = this->local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2),
                          u_nodal_index[c].imag());
                        // Get the hang weights
                        hang_weight2 = hang_info2_pt->master_weight(m2);
                      }
                      else
                      {
                        local_unknown =
                          this->nodal_local_eqn(l2, u_nodal_index[c].imag());
                        hang_weight2 = 1.0;
                      }

                      // If it's not pinned
                      if (local_unknown >= 0)
                      {
                        // Inertial term
                        if (a == c)
                        {
                          jacobian(local_eqn, local_unknown) -=
                            omega_sq_local * psi(l) * psi(l2) * W *
                            hang_weight * hang_weight2;
                        }

                        // Stress term
                        for (unsigned b = 0; b < DIM; b++)
                        {
                          for (unsigned d = 0; d < DIM; d++)
                          {
                            // Add the contribution to the Jacobian matrix
                            jacobian(local_eqn, local_unknown) +=
                              this->E(a, b, c, d) * dpsidx(l2, d) *
                              dpsidx(l, b) * W * hang_weight * hang_weight2;
                          }
                        }
                      } // End of if not boundary condition
                    }
                  }
                }
              } // End of jacobian calculation

            } // End of if not boundary condition


          } // End of loop over coordinate directions
        }
      } // End of loop over shape functions
    } // End of loop over integration points
  }


  //====================================================================
  /// Force building of required templates
  //====================================================================
  template class RefineableTimeHarmonicLinearElasticityEquations<2>;
  template class RefineableTimeHarmonicLinearElasticityEquations<3>;

} // namespace oomph
