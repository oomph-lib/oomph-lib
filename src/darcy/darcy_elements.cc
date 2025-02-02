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
#include "darcy_elements.h"

namespace oomph
{
  //=====================================================================
  /// Performs a div-conserving transformation of the vector basis
  /// functions from the reference element to the actual element
  //=====================================================================
  template<unsigned DIM>
  double DarcyEquations<DIM>::transform_basis(const Vector<double>& s,
                                              const Shape& q_basis_local,
                                              Shape& psi,
                                              Shape& q_basis) const
  {
    // Get the number of nodes in the element
    const unsigned n_node = this->nnode();

    // Storage for derivatives of the (geometric) shape functions, and call
    // the shape functions
    DShape dpsi(n_node, DIM);
    this->dshape_local(s, psi, dpsi);

    // Storage for the (geometric) jacobian and its inverse
    DenseMatrix<double> jacobian(DIM), inverse_jacobian(DIM);

    // Get the jacobian of the geometric mapping and its determinant
    double det = local_to_eulerian_mapping(dpsi, jacobian, inverse_jacobian);

    // Get the number of computational basis vectors
    const unsigned n_q_basis = this->nq_basis();

    // Loop over the basis vectors
    for (unsigned l = 0; l < n_q_basis; l++)
    {
      // Loop over the spatial components
      for (unsigned i = 0; i < DIM; i++)
      {
        // Zero the basis
        q_basis(l, i) = 0.0;
      }
    }

    // Loop over the spatial components
    for (unsigned i = 0; i < DIM; i++)
    {
      // And again
      for (unsigned j = 0; j < DIM; j++)
      {
        // Get the element of the jacobian (must transpose it due to different
        // conventions) and divide by the determinant
        double jac_trans = jacobian(j, i) / det;

        // Loop over the computational basis vectors
        for (unsigned l = 0; l < n_q_basis; l++)
        {
          // Apply the appropriate div-conserving mapping
          q_basis(l, i) += jac_trans * q_basis_local(l, j);
        }
      }
    }

    // Scale the basis by the ratio of the length of the edge to the length of
    // the corresponding edge on the reference element
    scale_basis(q_basis);

    return det;
  }

  //=====================================================================
  /// Output FE representation of soln: x,y,q1,q2,div_q,p,q \cdot n
  //=====================================================================
  template<unsigned DIM>
  void DarcyEquations<DIM>::output_with_projected_flux(
    std::ostream& outfile,
    const unsigned& nplot,
    const Vector<double> unit_normal)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);

    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Output the components of the position
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << this->interpolated_x(s, i) << " ";
      }

      // Output the components of the FE representation of q at s
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << this->interpolated_q(s, i) << " ";
      }

      // Output FE representation of div q at s
      outfile << this->interpolated_div_q(s) << " ";

      // Output FE representation of p at s
      outfile << this->interpolated_p(s) << " ";

      // Fluxes projected into the direction of the face normal
      double flux = 0.0;
      for (unsigned i = 0; i < 2; i++)
      {
        flux += this->interpolated_q(s, i) * unit_normal[i];
      }
      outfile << flux << " ";

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }


  //=====================================================================
  /// Output FE representation of soln: x,y,q1,q2,div_q,p at
  /// Nplot^DIM plot points
  //=====================================================================
  template<unsigned DIM>
  void DarcyEquations<DIM>::output(std::ostream& outfile, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);

    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Output the components of the position
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      // Output the components of the FE representation of q at s
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_q(s, i) << " ";
      }

      // Output FE representation of div q at s
      outfile << interpolated_div_q(s) << " ";

      // Output FE representation of p at s
      outfile << interpolated_p(s) << " ";

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  //=====================================================================
  /// Output FE representation of exact soln: x,y,q1,q2,div_q,p at
  /// Nplot^DIM plot points
  //=====================================================================
  template<unsigned DIM>
  void DarcyEquations<DIM>::output_fct(
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
    Vector<double> exact_soln(DIM + 2);

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

      // Output x,y,q_exact,p_exact,div_q_exact
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      for (unsigned i = 0; i < DIM + 2; i++)
      {
        outfile << exact_soln[i] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  //=====================================================================
  /// Compute the error between the FE solution and the exact solution
  /// using the H(div) norm for q and L^2 norm for p
  //=====================================================================
  template<unsigned DIM>
  void DarcyEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    Vector<double>& error,
    Vector<double>& norm)
  {
    for (unsigned i = 0; i < 2; i++)
    {
      error[i] = 0.0;
      norm[i] = 0.0;
    }

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordinates
    Vector<double> x(DIM);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    outfile << "ZONE" << std::endl;

    // Exact solution Vector (u,v,[w])
    Vector<double> exact_soln(DIM + 2);

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

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Flux error
      for (unsigned i = 0; i < DIM; i++)
      {
        norm[0] += exact_soln[i] * exact_soln[i] * W;
        // Error due to q_i
        error[0] += (exact_soln[i] - this->interpolated_q(s, i)) *
                    (exact_soln[i] - this->interpolated_q(s, i)) * W;
      }

      // // Flux divergence error
      // norm[0]+=exact_soln[DIM]*exact_soln[DIM]*W;
      // error[0]+=(exact_soln[DIM]-interpolated_div_q(s))*
      //           (exact_soln[DIM]-interpolated_div_q(s))*W;

      // Pressure error
      norm[1] += exact_soln[DIM + 1] * exact_soln[DIM + 1] * W;
      error[1] += (exact_soln[DIM + 1] - this->interpolated_p(s)) *
                  (exact_soln[DIM + 1] - this->interpolated_p(s)) * W;

      // Output x,y,[z]
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }

      // Output q_1_error,q_2_error,...
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << exact_soln[i] - this->interpolated_q(s, i) << " ";
      }

      // Output p_error
      outfile << exact_soln[DIM + 1] - this->interpolated_p(s) << " ";

      outfile << std::endl;
    }
  }

  //=====================================================================
  /// Fill in residuals and, if flag==true, jacobian
  //=====================================================================
  template<unsigned DIM>
  void DarcyEquations<DIM>::fill_in_generic_residual_contribution(
    Vector<double>& residuals, DenseMatrix<double>& jacobian, bool flag)
  {
    // Get the number of geometric nodes, total number of basis functions,
    // and number of edges basis functions
    const unsigned n_node = nnode();
    const unsigned n_q_basis = nq_basis();
    const unsigned n_q_basis_edge = nq_basis_edge();
    const unsigned n_p_basis = np_basis();

    // Storage for the geometric and computational bases and the test functions
    Shape psi(n_node), q_basis(n_q_basis, DIM), q_test(n_q_basis, DIM),
      p_basis(n_p_basis), p_test(n_p_basis), div_q_basis_ds(n_q_basis),
      div_q_test_ds(n_q_basis);

    // Get the number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Storage for the local coordinates
    Vector<double> s(DIM);

    // Storage for the source function
    Vector<double> f(DIM);

    // Storage for the mass source function
    double mass_source_local;

    // Local equation and unknown numbers
    int local_eqn = 0; //, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Find the local coordinates at the integration point
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the weight of the intetgration point
      double w = integral_pt()->weight(ipt);

      // Call the basis functions and test functions and get the
      // (geometric) jacobian of the current element
      double J = shape_basis_test_local_at_knot(ipt,
                                                psi,
                                                q_basis,
                                                q_test,
                                                p_basis,
                                                p_test,
                                                div_q_basis_ds,
                                                div_q_test_ds);

      // Storage for interpolated values
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> interpolated_q(DIM, 0.0);
      double interpolated_div_q_ds = 0.0;
      double interpolated_p = 0.0;

      // loop over the geometric basis functions to find interpolated x
      for (unsigned l = 0; l < n_node; l++)
      {
        for (unsigned i = 0; i < DIM; i++)
        {
          interpolated_x[i] += nodal_position(l, i) * psi[l];
        }
      }

      // loop over the nodes and use the vector basis functions to find the
      // interpolated flux
      for (unsigned i = 0; i < DIM; i++)
      {
        // Loop over the edge basis vectors
        for (unsigned l = 0; l < n_q_basis_edge; l++)
        {
          interpolated_q[i] += q_edge(l) * q_basis(l, i);
        }
        // Loop over the internal basis vectors
        for (unsigned l = n_q_basis_edge; l < n_q_basis; l++)
        {
          interpolated_q[i] += q_internal(l - n_q_basis_edge) * q_basis(l, i);
        }
      }

      // loop over the pressure basis and find the interpolated pressure
      for (unsigned l = 0; l < n_p_basis; l++)
      {
        interpolated_p += p_value(l) * p_basis(l);
      }

      // loop over the q edge divergence basis and the q internal divergence
      // basis to find interpolated div q
      for (unsigned l = 0; l < n_q_basis_edge; l++)
      {
        interpolated_div_q_ds += q_edge(l) * div_q_basis_ds(l);
      }
      for (unsigned l = n_q_basis_edge; l < n_q_basis; l++)
      {
        interpolated_div_q_ds +=
          q_internal(l - n_q_basis_edge) * div_q_basis_ds(l);
      }

      // Get the source function
      this->source(interpolated_x, f);

      // Get the mass source function
      this->mass_source(interpolated_x, mass_source_local);

      // Loop over the test functions
      for (unsigned l = 0; l < n_q_basis; l++)
      {
        if (l < n_q_basis_edge)
        {
          local_eqn = q_edge_local_eqn(l);
        }
        else // n_q_basis_edge <= l < n_basis
        {
          local_eqn = q_internal_local_eqn(l - n_q_basis_edge);
        }

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          for (unsigned i = 0; i < DIM; i++)
          {
            residuals[local_eqn] +=
              (interpolated_q[i] - f[i]) * q_test(l, i) * w * J;
          }

          // deliberately no jacobian factor in this integral
          residuals[local_eqn] -= (interpolated_p * div_q_test_ds(l)) * w;
        }
      } // End of loop over test functions

      // loop over pressure test functions
      for (unsigned l = 0; l < n_p_basis; l++)
      {
        // get the local equation number
        local_eqn = p_local_eqn(l);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // deliberately no jacobian factor in this integral
          residuals[local_eqn] += interpolated_div_q_ds * p_test(l) * w;
          residuals[local_eqn] -= mass_source_local * p_test(l) * w * J;
        }
      }
    } // End of loop over integration points
  }


  //=====================================================================
  // Force building of templates
  //=====================================================================
  template class DarcyEquations<2>;

} // namespace oomph
