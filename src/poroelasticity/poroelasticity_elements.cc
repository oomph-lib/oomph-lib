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
#include "poroelasticity_elements.h"

namespace oomph
{
  /// Static default value for timescale ratio (1.0 -- for natural scaling)
  template<unsigned DIM>
  double PoroelasticityEquations<DIM>::Default_lambda_sq_value = 1.0;

  /// Static default value for the density ratio
  template<unsigned DIM>
  double PoroelasticityEquations<DIM>::Default_density_ratio_value = 1.0;

  /// Static default value for 1/k
  template<unsigned DIM>
  double PoroelasticityEquations<DIM>::Default_k_inv_value = 1.0;

  /// Static default value for alpha
  template<unsigned DIM>
  double PoroelasticityEquations<DIM>::Default_alpha_value = 1.0;

  /// Static default value for the porosity
  template<unsigned DIM>
  double PoroelasticityEquations<DIM>::Default_porosity_value = 1.0;

  //======================================================================
  /// Compute the strain tensor at local coordinate s
  //======================================================================
  template<unsigned DIM>
  void PoroelasticityEquations<DIM>::get_strain(
    const Vector<double>& s, DenseMatrix<double>& strain) const
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

    /*
    //Find out how many position types there are
    unsigned n_position_type = this->nnodal_position_type();

    if(n_position_type != 1)
     {
      throw OomphLibError(
       "PoroElasticity is not yet implemented for more than one position type",
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
    */
#endif


    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_nodal_index[i] = u_index(i);
    }

    // Set up memory for the shape and derivative functions
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

    // Call the derivatives of the shape functions
    (void)dshape_eulerian(s, psi, dpsidx);

    // Calculate interpolated values of the derivative of global position
    DenseMatrix<double> interpolated_dudx(DIM, DIM, 0.0);

    // Loop over nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over velocity components
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get the nodal value
        const double u_value = this->nodal_value(l, u_nodal_index[i]);

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
  void PoroelasticityEquations<DIM>::get_stress(
    const Vector<double>& s, DenseMatrix<double>& stress) const
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
    DenseMatrix<double> strain(DIM, DIM);
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
          for (unsigned l = 0; k < DIM; k++)
          {
            stress(i, j) += this->E(i, j, k, l) * strain(k, l);
          }
        }
      }
    }
  }


  /// \short Performs a div-conserving transformation of the vector basis
  /// functions from the reference element to the actual element
  template<unsigned DIM>
  double PoroelasticityEquations<DIM>::transform_basis(
    const Vector<double>& s,
    const Shape& q_basis_local,
    Shape& psi,
    DShape& dpsi,
    Shape& q_basis) const
  {
    // Call the (geometric) shape functions and their derivatives
    //(void)this->dshape_eulerian(s,psi,dpsi);
    this->dshape_local(s, psi, dpsi);

    // Storage for the (geometric) jacobian and its inverse
    DenseMatrix<double> jacobian(DIM), inverse_jacobian(DIM);

    // Get the jacobian of the geometric mapping and its determinant
    double det = local_to_eulerian_mapping(dpsi, jacobian, inverse_jacobian);

    // Transform the derivative of the geometric basis so that it's w.r.t.
    // global coordinates
    this->transform_derivatives(inverse_jacobian, dpsi);

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

  /// \short Output FE representation of soln: x,y,u1,u2,q1,q2,div_q,p at
  /// Nplot^DIM plot points
  template<unsigned DIM>
  void PoroelasticityEquations<DIM>::output(std::ostream& outfile,
                                            const unsigned& nplot)
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

      // Output the components of the FE representation of u at s
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_u(s, i) << " ";
      }

      // Output the components of the FE representation of q at s
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_q(s, i) << " ";
      }

      // Output FE representation of div u at s
      outfile << interpolated_div_q(s) << " ";

      // Output FE representation of p at s
      outfile << interpolated_p(s) << " ";

      const unsigned n_node = this->nnode();

      Shape psi(n_node);
      shape(s, psi);

      Vector<double> interpolated_du_dt(DIM, 0.0);

      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned l = 0; l < n_node; l++)
        {
          interpolated_du_dt[i] += du_dt(l, i) * psi(l);
        }
        outfile << interpolated_du_dt[i] << " ";
      }

      Vector<double> interpolated_d2u_dt2(DIM, 0.0);

      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned l = 0; l < n_node; l++)
        {
          interpolated_d2u_dt2[i] += d2u_dt2(l, i) * psi(l);
        }
        outfile << interpolated_d2u_dt2[i] << " ";
      }

      const unsigned n_q_basis = this->nq_basis();
      const unsigned n_q_basis_edge = this->nq_basis_edge();

      Shape q_basis(n_q_basis, DIM), q_basis_local(n_q_basis, DIM);
      this->get_q_basis_local(s, q_basis_local);
      (void)this->transform_basis(s, q_basis_local, psi, q_basis);

      Vector<double> interpolated_dq_dt(DIM, 0.0);

      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned l = 0; l < n_q_basis_edge; l++)
        {
          interpolated_dq_dt[i] += dq_edge_dt(l) * q_basis(l, i);
        }
        for (unsigned l = n_q_basis_edge; l < n_q_basis; l++)
        {
          interpolated_dq_dt[i] +=
            dq_internal_dt(l - n_q_basis_edge) * q_basis(l, i);
        }
        outfile << interpolated_dq_dt[i] << " ";
      }

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  /// \short Output FE representation of exact soln: x,y,u1,u2,q1,q2,div_q,p at
  /// Nplot^DIM plot points
  template<unsigned DIM>
  void PoroelasticityEquations<DIM>::output_fct(
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
    Vector<double> exact_soln(2 * DIM + 2);

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
      for (unsigned i = 0; i < 2 * DIM + 2; i++)
      {
        outfile << exact_soln[i] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  /// \short Output FE representation of exact soln: x,y,u1,u2,div_q,p at
  /// Nplot^DIM plot points. Unsteady version
  template<unsigned DIM>
  void PoroelasticityEquations<DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Exact solution Vector
    Vector<double> exact_soln(2 * DIM + 2);

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

      // Output x,y,q_exact,p_exact,div_q_exact
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      for (unsigned i = 0; i < 2 * DIM + 2; i++)
      {
        outfile << exact_soln[i] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  /// \short Compute the error between the FE solution and the exact solution
  /// using the H(div) norm for q and L^2 norm for u and p
  template<unsigned DIM>
  void PoroelasticityEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    Vector<double>& error,
    Vector<double>& norm)
  {
    for (unsigned i = 0; i < 3; i++)
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
    Vector<double> exact_soln(2 * DIM + 2);

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

      // Displacement error
      for (unsigned i = 0; i < DIM; i++)
      {
        norm[0] += exact_soln[i] * exact_soln[i] * W;
        // Error due to q_i
        error[0] += (exact_soln[i] - this->interpolated_u(s, i)) *
                    (exact_soln[i] - this->interpolated_u(s, i)) * W;
      }

      // Flux error
      for (unsigned i = 0; i < DIM; i++)
      {
        norm[1] += exact_soln[DIM + i] * exact_soln[DIM + i] * W;
        // Error due to q_i
        error[1] += (exact_soln[DIM + i] - this->interpolated_q(s, i)) *
                    (exact_soln[DIM + i] - this->interpolated_q(s, i)) * W;
      }

      // Flux divergence error
      norm[1] += exact_soln[2 * DIM] * exact_soln[2 * DIM] * W;
      error[1] += (exact_soln[2 * DIM] - interpolated_div_q(s)) *
                  (exact_soln[2 * DIM] - interpolated_div_q(s)) * W;

      // Pressure error
      norm[2] += exact_soln[2 * DIM + 1] * exact_soln[2 * DIM + 1] * W;
      error[2] += (exact_soln[2 * DIM + 1] - this->interpolated_p(s)) *
                  (exact_soln[2 * DIM + 1] - this->interpolated_p(s)) * W;

      // Output x,y,[z]
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }

      // Output u_1_error,u_2_error,...
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << exact_soln[i] - this->interpolated_u(s, i) << " ";
      }

      // Output q_1_error,q_2_error,...
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << exact_soln[DIM + i] - this->interpolated_q(s, i) << " ";
      }

      // Output p_error
      outfile << exact_soln[2 * DIM + 1] - this->interpolated_p(s) << " ";

      outfile << std::endl;
    }
  }

  /// \short Compute the error between the FE solution and the exact solution
  /// using the H(div) norm for u and L^2 norm for p. Unsteady version
  template<unsigned DIM>
  void PoroelasticityEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
    Vector<double>& error,
    Vector<double>& norm)
  {
    for (unsigned i = 0; i < 3; i++)
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
    Vector<double> exact_soln(2 * DIM + 2);

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
      (*exact_soln_pt)(time, x, exact_soln);

      // Displacement error
      for (unsigned i = 0; i < DIM; i++)
      {
        norm[0] += exact_soln[i] * exact_soln[i] * W;
        // Error due to q_i
        error[0] += (exact_soln[i] - this->interpolated_u(s, i)) *
                    (exact_soln[i] - this->interpolated_u(s, i)) * W;
      }

      // Flux error
      for (unsigned i = 0; i < DIM; i++)
      {
        norm[1] += exact_soln[DIM + i] * exact_soln[DIM + i] * W;
        // Error due to q_i
        error[1] += (exact_soln[DIM + i] - this->interpolated_q(s, i)) *
                    (exact_soln[DIM + i] - this->interpolated_q(s, i)) * W;
      }

      // Flux divergence error
      norm[1] += exact_soln[2 * DIM] * exact_soln[2 * DIM] * W;
      error[1] += (exact_soln[2 * DIM] - interpolated_div_q(s)) *
                  (exact_soln[2 * DIM] - interpolated_div_q(s)) * W;

      // Pressure error
      norm[2] += exact_soln[2 * DIM + 1] * exact_soln[2 * DIM + 1] * W;
      error[2] += (exact_soln[2 * DIM + 1] - this->interpolated_p(s)) *
                  (exact_soln[2 * DIM + 1] - this->interpolated_p(s)) * W;

      // Output x,y,[z]
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }

      // Output u_1_error,u_2_error,...
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << exact_soln[i] - this->interpolated_u(s, i) << " ";
      }

      // Output q_1_error,q_2_error,...
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << exact_soln[DIM + i] - this->interpolated_q(s, i) << " ";
      }

      // Output p_error
      outfile << exact_soln[2 * DIM + 1] - this->interpolated_p(s) << " ";

      outfile << std::endl;
    }
  }

  /// Fill in residuals and, if flag==true, jacobian
  template<unsigned DIM>
  void PoroelasticityEquations<DIM>::fill_in_generic_residual_contribution(
    Vector<double>& residuals, DenseMatrix<double>& jacobian, bool flag)
  {
    // Get the number of geometric nodes, total number of basis functions,
    // and number of edges basis functions
    const unsigned n_node = nnode();
    const unsigned n_q_basis = nq_basis();
    const unsigned n_q_basis_edge = nq_basis_edge();
    const unsigned n_p_basis = np_basis();

    // Storage for the geometric and computational bases and the test functions
    Shape psi(n_node), u_basis(n_node), u_test(n_node), q_basis(n_q_basis, DIM),
      q_test(n_q_basis, DIM), p_basis(n_p_basis), p_test(n_p_basis),
      div_q_basis_ds(n_q_basis), div_q_test_ds(n_q_basis);

    DShape dpsidx(n_node, DIM), du_basis_dx(n_node, DIM),
      du_test_dx(n_node, DIM);

    // Get the number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Storage for the local coordinates
    Vector<double> s(DIM);

    // Storage for the elasticity source function
    Vector<double> f_solid(DIM);

    // Storage for the source function
    Vector<double> f_fluid(DIM);

    // Storage for the mass source function
    double mass_source_local;

    // Storage for Lambda_sq
    double lambda_sq = this->lambda_sq();

    // Get the value of 1/k
    double k_inv = this->k_inv();

    // Get the value of alpha
    double alpha = this->alpha();

    // Get the value of the porosity
    double porosity = this->porosity();

    // Get the density ratio
    double density_ratio = this->density_ratio();

    // Precompute the ratio of fluid density to combined density
    double rho_f_over_rho =
      density_ratio / (1 + porosity * (density_ratio - 1));

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0;

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
                                                dpsidx,
                                                u_basis,
                                                u_test,
                                                du_basis_dx,
                                                du_test_dx,
                                                q_basis,
                                                q_test,
                                                p_basis,
                                                p_test,
                                                div_q_basis_ds,
                                                div_q_test_ds);

      // Storage for interpolated values
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> interpolated_u(DIM, 0.0);
      DenseMatrix<double> interpolated_du_dx(DIM, DIM, 0.0);
      double interpolated_div_du_dt_dx = 0.0;
      Vector<double> interpolated_d2u_dt2(DIM, 0.0);
      Vector<double> interpolated_q(DIM, 0.0);
      double interpolated_div_q_ds = 0.0;
      Vector<double> interpolated_dq_dt(DIM, 0.0);
      double interpolated_p = 0.0;

      // loop over the vector basis functions to find interpolated x
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the geometric basis functions
        for (unsigned i = 0; i < DIM; i++)
        {
          interpolated_x[i] += nodal_position(l, i) * psi(l);

          interpolated_d2u_dt2[i] += this->d2u_dt2(l, i) * u_basis(l);

          // Get the nodal displacements
          const double u_value = this->raw_nodal_value(l, u_index(i));
          interpolated_u[i] += u_value * u_basis(l);

          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_du_dx(i, j) += u_value * du_basis_dx(l, j);
          }

          // divergence of the time derivative of the solid displacement
          interpolated_div_du_dt_dx += this->du_dt(l, i) * du_basis_dx(l, i);
        }
      }

      // loop over the nodes and use the geometric basis functions to find the
      // interpolated flux and its time derivative
      for (unsigned i = 0; i < DIM; i++)
      {
        // Loop over the edge basis vectors
        for (unsigned l = 0; l < n_q_basis_edge; l++)
        {
          interpolated_q[i] += q_edge(l) * q_basis(l, i);
          interpolated_dq_dt[i] += dq_edge_dt(l) * q_basis(l, i);
        }
        // Loop over the internal basis vectors
        for (unsigned l = n_q_basis_edge; l < n_q_basis; l++)
        {
          interpolated_q[i] += q_internal(l - n_q_basis_edge) * q_basis(l, i);
          interpolated_dq_dt[i] +=
            dq_internal_dt(l - n_q_basis_edge) * q_basis(l, i);
        }
      }

      // loop over the pressure basis and find the interpolated pressure
      for (unsigned l = 0; l < n_p_basis; l++)
      {
        interpolated_p += p_value(l) * p_basis(l);
      }

      // loop over the u edge divergence basis and the u internal divergence
      // basis to find interpolated div u
      for (unsigned l = 0; l < n_q_basis_edge; l++)
      {
        interpolated_div_q_ds += q_edge(l) * div_q_basis_ds(l);
      }
      for (unsigned l = n_q_basis_edge; l < n_q_basis; l++)
      {
        interpolated_div_q_ds +=
          q_internal(l - n_q_basis_edge) * div_q_basis_ds(l);
      }

      // Get the solid forcing
      this->force_solid(time, interpolated_x, f_solid);

      // Get the fluid forcing
      this->force_fluid(time, interpolated_x, f_fluid);

      // Get the mass source function
      this->mass_source(time, interpolated_x, mass_source_local);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Linear elasticity:

      for (unsigned l = 0; l < n_node; l++)
      {
        for (unsigned a = 0; a < DIM; a++)
        {
          local_eqn = this->nodal_local_eqn(l, u_index(a));

          if (local_eqn >= 0)
          {
            residuals[local_eqn] +=
              (lambda_sq * (interpolated_d2u_dt2[a] +
                            rho_f_over_rho * interpolated_dq_dt[a]) -
               f_solid[a]) *
              u_test(l) * W;

            // Stress term
            for (unsigned b = 0; b < DIM; b++)
            {
              if (a == b)
              {
                residuals[local_eqn] -=
                  alpha * interpolated_p * du_test_dx(l, b) * W;
              }

              for (unsigned c = 0; c < DIM; c++)
              {
                for (unsigned d = 0; d < DIM; d++)
                {
                  // Add the stress terms to the residuals
                  residuals[local_eqn] += this->E(a, b, c, d) *
                                          interpolated_du_dx(c, d) *
                                          du_test_dx(l, b) * W;
                }
              }
            }
            // Jacobian entries
            if (flag)
            {
              // d(u_eqn_l,a)/d(U_l2,c)
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                for (unsigned c = 0; c < DIM; c++)
                {
                  local_unknown = this->nodal_local_eqn(l2, u_index(c));
                  if (local_unknown >= 0)
                  {
                    if (a == c)
                    {
                      jacobian(local_eqn, local_unknown) +=
                        lambda_sq *
                        this->node_pt(l2)->time_stepper_pt()->weight(2, 0) *
                        u_basis(l2) * u_test(l) * W;
                    }

                    for (unsigned b = 0; b < DIM; b++)
                    {
                      for (unsigned d = 0; d < DIM; d++)
                      {
                        jacobian(local_eqn, local_unknown) +=
                          this->E(a, b, c, d) * du_basis_dx(l2, d) *
                          du_test_dx(l, b) * W;
                      }
                    }
                  }
                }
              }

              // d(u_eqn_l,a)/d(Q_l2)
              for (unsigned l2 = 0; l2 < n_q_basis; l2++)
              {
                TimeStepper* timestepper_pt = 0;

                if (l2 < n_q_basis_edge)
                {
                  local_unknown = q_edge_local_eqn(l2);
                  timestepper_pt =
                    this->node_pt(q_edge_node_number(l2))->time_stepper_pt();
                }
                else // n_q_basis_edge <= l < n_basis
                {
                  local_unknown = q_internal_local_eqn(l2 - n_q_basis_edge);
                  timestepper_pt = this->internal_data_pt(q_internal_index())
                                     ->time_stepper_pt();
                }

                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    lambda_sq * rho_f_over_rho * timestepper_pt->weight(1, 0) *
                    q_basis(l2, a) * u_test(l) * W;
                }
              }

              // d(u_eqn_l,a)/d(P_l2)
              for (unsigned l2 = 0; l2 < n_p_basis; l2++)
              {
                local_unknown = p_local_eqn(l2);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    alpha * p_basis(l2) * du_test_dx(l, a) * W;
                }
              }
            } // End of Jacobian entries
          } // End of if not boundary condition
        } // End of loop over dimensions
      } // End of loop over u test functions


      // Darcy:
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
              (k_inv * interpolated_q[i] - rho_f_over_rho * f_fluid[i] +
               rho_f_over_rho * lambda_sq *
                 (interpolated_d2u_dt2[i] +
                  (1.0 / porosity) * interpolated_dq_dt[i])) *
              q_test(l, i) * w * J;
          }

          // deliberately no jacobian factor in this integral
          residuals[local_eqn] -= (interpolated_p * div_q_test_ds(l)) * w;

          // Jacobian entries
          if (flag)
          {
            // d(q_eqn_l)/d(U_l2,c)
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              for (unsigned c = 0; c < DIM; c++)
              {
                local_unknown = this->nodal_local_eqn(l2, u_index(c));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    rho_f_over_rho * lambda_sq *
                    this->node_pt(l2)->time_stepper_pt()->weight(2, 0) *
                    u_basis(l2) * q_test(l, c) * W;
                }
              }
            }

            // d(q_eqn_l)/d(Q_l2)
            for (unsigned l2 = 0; l2 < n_q_basis; l2++)
            {
              TimeStepper* timestepper_pt = 0;

              if (l2 < n_q_basis_edge)
              {
                local_unknown = q_edge_local_eqn(l2);
                timestepper_pt =
                  this->node_pt(q_edge_node_number(l2))->time_stepper_pt();
              }
              else // n_q_basis_edge <= l < n_basis
              {
                local_unknown = q_internal_local_eqn(l2 - n_q_basis_edge);
                timestepper_pt =
                  this->internal_data_pt(q_internal_index())->time_stepper_pt();
              }

              if (local_unknown >= 0)
              {
                for (unsigned c = 0; c < DIM; c++)
                {
                  jacobian(local_eqn, local_unknown) +=
                    q_basis(l2, c) * q_test(l, c) *
                    (k_inv + rho_f_over_rho * lambda_sq *
                               timestepper_pt->weight(1, 0) / porosity) *
                    W;
                }
              }
            }

            // d(q_eqn_l)/d(P_l2)
            for (unsigned l2 = 0; l2 < n_p_basis; l2++)
            {
              local_unknown = p_local_eqn(l2);

              // deliberately no jacobian factor in this term
              jacobian(local_eqn, local_unknown) -=
                p_basis(l2) * div_q_test_ds(l) * w;
            }
          } // End of Jacobian entries
        } // End of if not boundary condition
      } // End of loop over q test functions

      // loop over pressure test functions
      for (unsigned l = 0; l < n_p_basis; l++)
      {
        // get the local equation number
        local_eqn = p_local_eqn(l);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // solid divergence term
          residuals[local_eqn] +=
            alpha * interpolated_div_du_dt_dx * p_test(l) * w * J;
          // fluid divergence term - deliberately no jacobian factor in this
          // term
          residuals[local_eqn] += interpolated_div_q_ds * p_test(l) * w;
          // mass source term
          residuals[local_eqn] -= mass_source_local * p_test(l) * w * J;

          // Jacobian entries
          if (flag)
          {
            // d(p_eqn_l)/d(U_l2,c)
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              for (unsigned c = 0; c < DIM; c++)
              {
                local_unknown = this->nodal_local_eqn(l2, u_index(c));

                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    alpha * this->node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                    du_basis_dx(l2, c) * p_test(l) * W;
                }
              }
            }

            // d(p_eqn_l)/d(Q_l2)
            for (unsigned l2 = 0; l2 < n_q_basis; l2++)
            {
              if (l2 < n_q_basis_edge)
              {
                local_unknown = q_edge_local_eqn(l2);
              }
              else // n_q_basis_edge <= l < n_basis
              {
                local_unknown = q_internal_local_eqn(l2 - n_q_basis_edge);
              }

              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  p_test(l) * div_q_basis_ds(l2) * w;
              }
            }
          } // End of Jacobian entries
        } // End of if not boundary condition
      } // End of loop over p test functions
    } // End of loop over integration points
  }

  // Force building of templates
  template class PoroelasticityEquations<2>;

} // namespace oomph
