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
// Non-inline functions for NS elements

#include "spherical_navier_stokes_elements.h"


namespace oomph
{
  /// Navier--Stokes equations static data
  Vector<double> SphericalNavierStokesEquations::Gamma(3, 1.0);

  //=================================================================
  /// "Magic" negative number that indicates that the pressure is
  /// not stored at a node. This cannot be -1 because that represents
  /// the positional hanging scheme in the hanging_pt object of nodes
  //=================================================================
  int SphericalNavierStokesEquations::Pressure_not_stored_at_node = -100;

  /// Navier--Stokes equations static data
  double SphericalNavierStokesEquations::Default_Physical_Constant_Value = 0.0;

  /// Navier--Stokes equations static data
  double SphericalNavierStokesEquations::Default_Physical_Ratio_Value = 1.0;

  /// Navier-Stokes equations default gravity vector
  Vector<double> SphericalNavierStokesEquations::Default_Gravity_vector(3, 0.0);


  //================================================================
  /// Compute the diagonal of the velocity/pressure mass matrices.
  /// If which one=0, both are computed, otherwise only the pressure
  /// (which_one=1) or the velocity mass matrix (which_one=2 -- the
  /// LSC version of the preconditioner only needs that one)
  /// NOTE: pressure versions isn't implemented yet because this
  ///       element has never been tried with Fp preconditoner.
  //================================================================
  void SphericalNavierStokesEquations::
    get_pressure_and_velocity_mass_matrix_diagonal(
      Vector<double>& press_mass_diag,
      Vector<double>& veloc_mass_diag,
      const unsigned& which_one)
  {
#ifdef PARANOID
    if ((which_one == 0) || (which_one == 1))
    {
      throw OomphLibError("Computation of diagonal of pressure mass matrix is "
                          "not impmented yet.\n",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Resize and initialise
    veloc_mass_diag.assign(ndof(), 0.0);

    // find out how many nodes there are
    const unsigned n_node = nnode();

    // find number of coordinates
    const unsigned n_value = 3;

    // find the indices at which the local velocities are stored
    Vector<unsigned> u_nodal_index(n_value);
    for (unsigned i = 0; i < n_value; i++)
    {
      u_nodal_index[i] = this->u_index_spherical_nst(i);
    }

    // Set up memory for test functions
    Shape test(n_node);

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Integer to store the local equations no
    int local_eqn = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get determinant of Jacobian of the mapping
      double J = J_eulerian_at_knot(ipt);

      // Premultiply weights and Jacobian
      double W = w * J;

      // Get the velocity test functions - there is no explicit
      // function to give the test function so use shape
      shape_at_knot(ipt, test);

      // Need to get the position to sort out the jacobian properly
      double r = 0.0;
      double theta = 0.0;
      for (unsigned l = 0; l < n_node; l++)
      {
        r += this->nodal_position(l, 0) * test(l);
        theta += this->nodal_position(l, 1) * test(l);
      }

      // Multiply by the geometric factor
      W *= r * r * sin(theta);

      // Loop over the veclocity test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the velocity components
        for (unsigned i = 0; i < n_value; i++)
        {
          local_eqn = nodal_local_eqn(l, u_nodal_index[i]);

          // If not a boundary condition
          if (local_eqn >= 0)
          {
            // Add the contribution
            veloc_mass_diag[local_eqn] += test[l] * test[l] * W;
          } // End of if not boundary condition statement
        } // End of loop over dimension
      } // End of loop over test functions
    }
  }


  //======================================================================
  /// Validate against exact velocity solution at given time.
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points and compute L2 error
  /// and L2 norm of velocity solution over element.
  //=======================================================================
  void SphericalNavierStokesEquations::compute_error(
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

    // Vector for coordintes
    Vector<double> x(2);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    outfile << "ZONE" << std::endl;

    // Exact solution Vector (u,v,[w],p)
    Vector<double> exact_soln(4);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(time, x, exact_soln);

      // Velocity error
      for (unsigned i = 0; i < 3; i++)
      {
        norm += exact_soln[i] * exact_soln[i] * W;
        error += (exact_soln[i] - interpolated_u_spherical_nst(s, i)) *
                 (exact_soln[i] - interpolated_u_spherical_nst(s, i)) * W;
      }

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output x,y,[z],u_error,v_error,[w_error]
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << exact_soln[i] - interpolated_u_spherical_nst(s, i) << " ";
      }

      outfile << std::endl;
    }
  }


  //======================================================================
  /// Validate against exact velocity solution
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points and compute L2 error
  /// and L2 norm of velocity solution over element.
  //=======================================================================
  void SphericalNavierStokesEquations::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordintes
    Vector<double> x(2);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();


    outfile << "ZONE" << std::endl;

    // Exact solution Vector (u,v,w,p)
    Vector<double> exact_soln(4);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      interpolated_x(s, x);

      // Get the exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Velocity error
      for (unsigned i = 0; i < 3; i++)
      {
        norm += exact_soln[i] * exact_soln[i] * W;
        error += (exact_soln[i] - interpolated_u_spherical_nst(s, i)) *
                 (exact_soln[i] - interpolated_u_spherical_nst(s, i)) * W;
      }

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output x,y,[z],u_error,v_error,[w_error]
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << exact_soln[i] - interpolated_u_spherical_nst(s, i) << " ";
      }
      outfile << std::endl;
    }
  }

  //======================================================================
  /// Validate against exact velocity solution
  /// Solution is provided via direct coding.
  /// Plot at a given number of plot points and compute energy error
  /// and energy norm of velocity solution over element.
  //=======================================================================
  void SphericalNavierStokesEquations::compute_error_e(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_dr_pt,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_dtheta_pt,
    double& u_error,
    double& u_norm,
    double& p_error,
    double& p_norm)
  {
    u_error = 0.0;
    p_error = 0.0;
    u_norm = 0.0;
    p_norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordinates
    Vector<double> x(2);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();


    /// Exact solution Vectors (u,v,[w],p) -
    /// - need two extras to deal with the derivatives in the energy norm
    Vector<double> exact_soln(4);
    Vector<double> exact_soln_dr(4);
    Vector<double> exact_soln_dth(4);


    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      interpolated_x(s, x);

      // Get exact solution at the point x using the 'actual' functions
      (*exact_soln_pt)(x, exact_soln);
      (*exact_soln_dr_pt)(x, exact_soln_dr);
      (*exact_soln_dtheta_pt)(x, exact_soln_dth);

      // exact_soln = actual(x);
      // exact_soln_dr = actual_dr(x);
      // exact_soln_dth = actual_dth(x);


      double r = x[0];
      double theta = x[1];

      W *= r * r * sin(theta);


      // Velocity error in energy norm
      // Owing to the additional terms arising from the phi-components, we can't
      // loop over the velocities. Calculate each section separately instead.

      // Norm calculation involves the double contraction of two tensors, so we
      // take the summation of nine separate components for clarity. Fortunately
      // most of the components zero out in the axisymmetric case.
      //---------------------------------------------
      // Entry (1,1) in dyad matrix
      // No contribution
      // Entry (1,2) in dyad matrix
      // No contribution
      // Entry (1,3) in dyad matrix
      u_norm += (1 / (r * r)) * exact_soln[2] * exact_soln[2] * W;
      //----------------------------------------------
      // Entry (2,1) in dyad matrix
      // No contribution
      // Entry (2,2) in dyad matrix
      // No contribution
      // Entry (2,3) in dyad matrix
      u_norm += cot(theta) * cot(theta) * (1 / (r * r)) * exact_soln[2] *
                exact_soln[2] * W;
      //-----------------------------------------------
      // Entry (3,1) in dyad matrix
      u_norm += exact_soln_dr[2] * exact_soln_dr[2] * W;
      // Entry (3,2) in dyad matrix
      u_norm += (1 / (r * r)) * exact_soln_dth[2] * exact_soln_dth[2] * W;
      // Entry (3,3) in dyad matrix
      // No contribution
      //-----------------------------------------------
      // End of norm calculation


      // Error calculation involves the summation of 9 squared components,
      // stated separately here for clarity
      //---------------------------------------------
      // Entry (1,1) in dyad matrix
      u_error += interpolated_dudx_spherical_nst(s, 0, 0) *
                 interpolated_dudx_spherical_nst(s, 0, 0) * W;
      // Entry (1,2) in dyad matrix
      u_error += (1 / (r * r)) *
                 (interpolated_u_spherical_nst(s, 1) -
                  interpolated_dudx_spherical_nst(s, 0, 1)) *
                 (interpolated_u_spherical_nst(s, 1) -
                  interpolated_dudx_spherical_nst(s, 0, 1)) *
                 W;
      // Entry (1,3) in dyad matrix
      u_error += (1 / (r * r)) *
                 (interpolated_u_spherical_nst(s, 2) - exact_soln[2]) *
                 (interpolated_u_spherical_nst(s, 2) - exact_soln[2]) * W;
      //---------------------------------------------
      // Entry (2,1) in dyad matrix
      u_error += interpolated_dudx_spherical_nst(s, 1, 0) *
                 interpolated_dudx_spherical_nst(s, 1, 0) * W;
      // Entry (2,2) in dyad matrix
      u_error += (1 / (r * r)) *
                 (interpolated_dudx_spherical_nst(s, 1, 1) +
                  interpolated_u_spherical_nst(s, 0)) *
                 (interpolated_dudx_spherical_nst(s, 1, 1) +
                  interpolated_u_spherical_nst(s, 0)) *
                 W;
      // Entry (2,3) in dyad matrix
      u_error += (1 / (r * r)) * cot(theta) * cot(theta) *
                 (interpolated_u_spherical_nst(s, 2) - exact_soln[2]) *
                 (interpolated_u_spherical_nst(s, 2) - exact_soln[2]) * W;
      //---------------------------------------------
      // Entry (3,1) in dyad matrix
      u_error += (exact_soln_dr[2] - interpolated_dudx_spherical_nst(s, 2, 0)) *
                 (exact_soln_dr[2] - interpolated_dudx_spherical_nst(s, 2, 0)) *
                 W;
      // Entry (3,2) in dyad matrix
      u_error +=
        (1 / (r * r)) *
        (exact_soln_dth[2] - interpolated_dudx_spherical_nst(s, 2, 1)) *
        (exact_soln_dth[2] - interpolated_dudx_spherical_nst(s, 2, 1)) * W;
      // Entry (3,3) in dyad matrix
      u_error += (1 / (r * r)) *
                 (interpolated_u_spherical_nst(s, 0) -
                  cot(theta) * interpolated_u_spherical_nst(s, 1)) *
                 (interpolated_u_spherical_nst(s, 0) -
                  cot(theta) * interpolated_u_spherical_nst(s, 1)) *
                 W;
      //---------------------------------------------
      // End of velocity error calculation

      // Pressure error in L2 norm - energy norm not required here as the
      // pressure only appears undifferentiated in the NS SPC weak form.

      p_norm += exact_soln[3] * exact_soln[3] * W;

      p_error += (exact_soln[3] - interpolated_p_spherical_nst(s)) *
                 (exact_soln[3] - interpolated_p_spherical_nst(s)) * W;


      // Output r and theta coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output the velocity error details in the energy norm
      outfile << u_error << "  " << u_norm << "  ";
      // Output the pressure error details in the L2 norm
      outfile << p_error << "  " << p_norm << "  ";
      outfile << std::endl;

    } // End loop over integration points

  } // End of function statement


  //======================================================================
  // Output the shear stress at the outer wall of a rotating chamber.
  //======================================================================
  void SphericalNavierStokesEquations::compute_shear_stress(
    std::ostream& outfile)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordinates
    Vector<double> x(2);

    // Number of plot points
    unsigned Npts = 3;

    // Declaration of the shear stress variable
    double shear = 0.0;


    // Assign values of s
    for (unsigned i = 0; i < Npts; i++)
    {
      s[0] = 1.0;
      s[1] = -1 + 2.0 * i / (Npts - 1);


      // Get position as global coordinate
      interpolated_x(s, x);

      // Output r and theta coordinates
      outfile << x[0] << " ";
      outfile << x[1] << " ";

      // Output the shear stress at the current coordinate
      double r = x[0];
      shear = interpolated_dudx_spherical_nst(s, 2, 0) -
              (1 / r) * interpolated_u_spherical_nst(s, 2);
      outfile << shear << "  ";

      outfile << std::endl;
    }
  }


  //======================================================================
  // Output given velocity values and shear stress along a specific
  // section of a rotating chamber.
  //======================================================================
  void SphericalNavierStokesEquations::extract_velocity(std::ostream& outfile)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordinates
    Vector<double> x(2);

    // Number of plot points
    unsigned Npts = 3;


    // Assign values of s
    for (unsigned i = 0; i < Npts; i++)
    {
      s[1] = 1.0;
      s[0] = -1 + 2.0 * i / (Npts - 1);


      // Get position as global coordinate
      interpolated_x(s, x);
      double r = x[0];
      double theta = x[1];

      // Output r and theta coordinates
      outfile << r << " ";
      outfile << theta << " ";

      // Output the velocity at the current coordinate
      outfile << interpolated_u_spherical_nst(s, 0) << "  ";
      outfile << interpolated_u_spherical_nst(s, 1) << "  ";
      outfile << interpolated_u_spherical_nst(s, 2) << "  ";
      outfile << interpolated_p_spherical_nst(s) << "  ";

      // Output shear stress along the given radius
      double shear = interpolated_dudx_spherical_nst(s, 2, 0) -
                     (1 / r) * interpolated_u_spherical_nst(s, 2);
      outfile << shear << "  ";

      // Coordinate setup for norm
      double J = J_eulerian(s);
      if (i != 1)
      {
        double w = integral_pt()->weight(i);
        double W = w * J;
        W *= r * r * sin(theta);

        // Output the velocity integral at current point
        double u_int = 0.0;
        for (unsigned j = 0; j < 3; j++)
        {
          u_int += interpolated_u_spherical_nst(s, j) * W;
        }
        outfile << u_int << " ";

        // Output the pressure integral at current point
        double p_int = interpolated_p_spherical_nst(s) * W;
        outfile << p_int;
        outfile << std::endl;
      }

      else
      {
        double w = integral_pt()->weight(5);
        double W = w * J;
        W *= r * r * sin(theta);

        // Output the velocity integral at current point
        double u_int = 0.0;
        for (unsigned j = 0; j < 3; j++)
        {
          u_int += interpolated_u_spherical_nst(s, j) * W;
        }
        outfile << u_int << " ";

        // Output the pressure integral at current point
        double p_int = interpolated_p_spherical_nst(s) * W;
        outfile << p_int;
        outfile << std::endl;
      }
    }
  }


  //======================================================================
  /// Output "exact" solution
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  /// Function prints as many components as are returned in solution Vector.
  //=======================================================================
  void SphericalNavierStokesEquations::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordintes
    Vector<double> x(2);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Exact solution Vector
    Vector<double> exact_soln;


    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,...
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << x[i] << " ";
      }

      // Output "exact solution"
      for (unsigned i = 0; i < exact_soln.size(); i++)
      {
        outfile << exact_soln[i] << " ";
      }

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }

  //======================================================================
  /// Output "exact" solution at a given time
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  /// Function prints as many components as are returned in solution Vector.
  //=======================================================================
  void SphericalNavierStokesEquations::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Vector for coordinates
    Vector<double> x(2);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Exact solution Vector
    Vector<double> exact_soln;

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(time, x, exact_soln);

      // Output x,y,...
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << x[i] << " ";
      }

      // Output "exact solution"
      for (unsigned i = 0; i < exact_soln.size(); i++)
      {
        outfile << exact_soln[i] << " ";
      }

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }

  //==============================================================
  /// Output function: Velocities only
  /// x,y,[z],u,v,[w]
  /// in tecplot format at specified previous timestep (t=0: present;
  /// t>0: previous timestep). Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  void SphericalNavierStokesEquations::output_veloc(std::ostream& outfile,
                                                    const unsigned& nplot,
                                                    const unsigned& t)
  {
    // Find number of nodes
    unsigned n_node = nnode();

    // Local shape function
    Shape psi(n_node);

    // Vectors of local coordinates and coords and velocities
    Vector<double> s(2);
    Vector<double> interpolated_x(2);
    Vector<double> interpolated_u(3);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get shape functions
      shape(s, psi);

      // Loop over directions
      for (unsigned i = 0; i < 3; i++)
      {
        interpolated_x[i] = 0.0;
        interpolated_u[i] = 0.0;
        // Get the index at which velocity i is stored
        unsigned u_nodal_index = u_index_spherical_nst(i);
        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          interpolated_u[i] += nodal_value(t, l, u_nodal_index) * psi[l];
          interpolated_x[i] += nodal_position(t, l, i) * psi[l];
        }
      }

      // Coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_x[i] << " ";
      }

      // Velocities
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << interpolated_u[i] << " ";
      }

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }

  //==============================================================
  /// Output function:
  /// x,y,[z],u,v,[w],p
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  void SphericalNavierStokesEquations::output(std::ostream& outfile,
                                              const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Coordinates in Cartesian Coordinates
      // for(unsigned i=0;i<DIM;i++)
      // {
      //  outfile << interpolated_x(s,i) << " ";
      // }

      // Output the coordinates in Cartesian coordintes
      double r = interpolated_x(s, 0);
      double theta = interpolated_x(s, 1);
      // The actual Cartesian values
      outfile << r * sin(theta) << " " << r * cos(theta) << " ";

      // Now the x and y velocities
      double u_r = interpolated_u_spherical_nst(s, 0);
      double u_theta = interpolated_u_spherical_nst(s, 1);

      outfile << u_r * sin(theta) + u_theta * cos(theta) << " "
              << u_r * cos(theta) - u_theta * sin(theta) << " ";

      // Now the polar coordinates
      // outfile << r << " " << theta << " ";

      // Velocities r, theta, phi
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << interpolated_u_spherical_nst(s, i) << " ";
      }

      // Pressure
      outfile << interpolated_p_spherical_nst(s) << " ";

      // Dissipation
      outfile << dissipation(s) << " ";

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //==============================================================
  /// C-style output function:
  /// x,y,[z],u,v,[w],p
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  void SphericalNavierStokesEquations::output(FILE* file_pt,
                                              const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Tecplot header info
    fprintf(file_pt, "%s", tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }

      // Velocities
      for (unsigned i = 0; i < 3; i++)
      {
        fprintf(file_pt, "%g ", interpolated_u_spherical_nst(s, i));
      }

      // Pressure
      fprintf(file_pt, "%g \n", interpolated_p_spherical_nst(s));
    }
    fprintf(file_pt, "\n");

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, nplot);
  }


  //==============================================================
  /// Full output function:
  /// x,y,[z],u,v,[w],p,du/dt,dv/dt,[dw/dt],dissipation
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction
  //==============================================================
  void SphericalNavierStokesEquations::full_output(std::ostream& outfile,
                                                   const unsigned& nplot)
  {
    throw OomphLibError(
      "Probably Broken", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);

    // Vector of local coordinates
    Vector<double> s(2);

    // Acceleration
    Vector<double> dudt(3);

    // Mesh elocity
    Vector<double> mesh_veloc(3, 0.0);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Set up memory for the shape functions
    Shape psif(n_node);
    DShape dpsifdx(n_node, 2);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psif, dpsifdx);

      // Allocate storage
      Vector<double> mesh_veloc(3);
      Vector<double> dudt(3);
      Vector<double> dudt_ALE(3);
      DenseMatrix<double> interpolated_dudx(3, 2);

      // Initialise everything to zero
      for (unsigned i = 0; i < 3; i++)
      {
        mesh_veloc[i] = 0.0;
        dudt[i] = 0.0;
        dudt_ALE[i] = 0.0;
        for (unsigned j = 0; j < 2; j++)
        {
          interpolated_dudx(i, j) = 0.0;
        }
      }

      // Calculate velocities and derivatives


      // Loop over directions
      for (unsigned i = 0; i < 3; i++)
      {
        // Get the index at which velocity i is stored
        unsigned u_nodal_index = u_index_spherical_nst(i);
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          dudt[i] += du_dt_spherical_nst(l, u_nodal_index) * psif(l);
          mesh_veloc[i] += dnodal_position_dt(l, i) * psif(l);

          // Loop over derivative directions for velocity gradients
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) +=
              nodal_value(l, u_nodal_index) * dpsifdx(l, j);
          }
        }
      }


      // Get dudt in ALE form (incl mesh veloc)
      for (unsigned i = 0; i < 3; i++)
      {
        dudt_ALE[i] = dudt[i];
        for (unsigned k = 0; k < 2; k++)
        {
          dudt_ALE[i] -= mesh_veloc[k] * interpolated_dudx(i, k);
        }
      }


      // Coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      // Velocities
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << interpolated_u_spherical_nst(s, i) << " ";
      }

      // Pressure
      outfile << interpolated_p_spherical_nst(s) << " ";

      // Accelerations
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << dudt_ALE[i] << " ";
      }

      // Dissipation
      outfile << dissipation(s) << " ";


      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //==============================================================
  /// Output function for vorticity.
  /// x,y,[z],[omega_x,omega_y,[and/or omega_z]]
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  void SphericalNavierStokesEquations::output_vorticity(std::ostream& outfile,
                                                        const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Create vorticity vector of the required size
    Vector<double> vorticity(1);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      // Get vorticity
      get_vorticity(s, vorticity);

      outfile << vorticity[0] << std::endl;
      ;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //==============================================================
  /// Return integral of dissipation over element
  //==============================================================
  double SphericalNavierStokesEquations::dissipation() const
  {
    // Initialise
    double diss = 0.0;

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Get strain rate matrix
      DenseMatrix<double> strainrate(3, 3);
      strain_rate(s, strainrate);

      // Initialise
      double local_diss = 0.0;
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = 0; j < 3; j++)
        {
          local_diss += 2.0 * strainrate(i, j) * strainrate(i, j);
        }
      }

      diss += local_diss * w * J;
    }

    return diss;
  }

  //==============================================================
  /// Compute traction (on the viscous scale) exerted onto
  /// the fluid at local coordinate s. N has to be outer unit normal
  /// to the fluid.
  //==============================================================
  void SphericalNavierStokesEquations::get_traction(const Vector<double>& s,
                                                    const Vector<double>& N,
                                                    Vector<double>& traction)
  {
    // Get velocity gradients
    DenseMatrix<double> strainrate(3, 3);
    strain_rate(s, strainrate);

    // Get pressure
    double press = interpolated_p_spherical_nst(s);

    // Loop over traction components
    for (unsigned i = 0; i < 3; i++)
    {
      // If we are in the r and theta direction
      // initialise the traction
      if (i < 2)
      {
        traction[i] = -press * N[i];
      }
      // Otherwise it's zero because the normal cannot have a component
      // in the phi direction
      else
      {
        traction[i] = 0.0;
      }
      // Loop over the possible traction directions
      for (unsigned j = 0; j < 2; j++)
      {
        traction[i] += 2.0 * strainrate(i, j) * N[j];
      }
    }
  }

  //==============================================================
  /// Return dissipation at local coordinate s
  //==============================================================
  double SphericalNavierStokesEquations::dissipation(
    const Vector<double>& s) const
  {
    // Get strain rate matrix
    DenseMatrix<double> strainrate(3, 3);
    strain_rate(s, strainrate);

    // Initialise
    double local_diss = 0.0;
    for (unsigned i = 0; i < 3; i++)
    {
      for (unsigned j = 0; j < 3; j++)
      {
        local_diss += 2.0 * strainrate(i, j) * strainrate(i, j);
      }
    }

    return local_diss;
  }

  //==============================================================
  /// Get strain-rate tensor: 1/2 (du_i/dx_j + du_j/dx_i)
  //==============================================================
  void SphericalNavierStokesEquations::strain_rate(
    const Vector<double>& s, DenseMatrix<double>& strainrate) const
  {
    // Velocity components
    Vector<double> interpolated_u(3, 0.0);
    // Coordinates
    double interpolated_r = 0.0;
    double interpolated_theta = 0.0;
    // Velocity gradient matrix
    DenseMatrix<double> dudx(3, 2, 0.0);


    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node, 2);

    // Call the derivatives of the shape functions
    dshape_eulerian(s, psi, dpsidx);

    // Get the values of the positions
    for (unsigned l = 0; l < n_node; l++)
    {
      interpolated_r += this->nodal_position(l, 0) * psi(l);
      interpolated_theta += this->nodal_position(l, 1) * psi(l);
    }


    // Loop over velocity components
    for (unsigned i = 0; i < 3; i++)
    {
      // Get the index at which the i-th velocity is stored
      unsigned u_nodal_index = u_index_spherical_nst(i);
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the velocity value
        const double nodal_u = this->nodal_value(l, u_nodal_index);
        // Add in the velocities
        interpolated_u[i] += nodal_u * psi(l);
        // Loop over the derivative directions
        for (unsigned j = 0; j < 2; j++)
        {
          dudx(i, j) += nodal_u * dpsidx(l, j);
        }
      }
    }

    // Assign those strain rates without any negative powers of the radius
    strainrate(0, 0) = dudx(0, 0);
    strainrate(0, 1) = 0.0;
    strainrate(0, 2) = 0.0;
    ;
    strainrate(1, 1) = 0.0;
    strainrate(1, 2) = 0.0;
    strainrate(2, 2) = 0.0;


    // Add the negative powers of the radius, unless we are
    // at the origin
    if (std::fabs(interpolated_r) > 1.0e-15)
    {
      const double pi = MathematicalConstants::Pi;
      double inverse_r = 1.0 / interpolated_r;

      // Should we include the cot terms (default no)
      bool include_cot_terms = false;
      double cot_theta = 0.0;
      // If we in the legal range then include cot
      if ((std::fabs(interpolated_theta) > 1.0e-15) &&
          (std::fabs(pi - interpolated_theta) > 1.0e-15))
      {
        include_cot_terms = true;
        cot_theta = this->cot(interpolated_theta);
      }

      strainrate(0, 1) =
        0.5 * (dudx(1, 0) + (dudx(0, 1) - interpolated_u[1]) * inverse_r);
      strainrate(0, 2) = 0.5 * (dudx(2, 0) - interpolated_u[2] * inverse_r);
      strainrate(1, 1) = (dudx(1, 1) + interpolated_u[0]) * inverse_r;

      if (include_cot_terms)
      {
        strainrate(1, 2) =
          0.5 * (dudx(2, 1) - interpolated_u[2] * cot_theta) * inverse_r;
        strainrate(2, 2) =
          (interpolated_u[0] + interpolated_u[1] * cot_theta) * inverse_r;
      }
    }

    // Sort out the symmetries
    strainrate(1, 0) = strainrate(0, 1);
    strainrate(2, 0) = strainrate(0, 2);
    strainrate(2, 1) = strainrate(1, 2);
  }


  //==============================================================
  /// Compute 2D vorticity vector at local coordinate s (return in
  /// one and only component of vorticity vector
  //==============================================================
  void SphericalNavierStokesEquations::get_vorticity(
    const Vector<double>& s, Vector<double>& vorticity) const
  {
    throw OomphLibError("Not tested for spherical Navier-Stokes elements",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);

#ifdef PARANOID
    if (vorticity.size() != 1)
    {
      std::ostringstream error_message;
      error_message << "Computation of vorticity in 2D requires a 1D vector\n"
                    << "which contains the only non-zero component of the\n"
                    << "vorticity vector. You've passed a vector of size "
                    << vorticity.size() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Specify spatial dimension
    unsigned DIM = 2;

    // Velocity gradient matrix
    DenseMatrix<double> dudx(DIM);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

    // Call the derivatives of the shape functions
    dshape_eulerian(s, psi, dpsidx);

    // Initialise to zero
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        dudx(i, j) = 0.0;
      }
    }

    // Loop over veclocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Get the index at which the i-th velocity is stored
      unsigned u_nodal_index = u_index_spherical_nst(i);
      // Loop over derivative directions
      for (unsigned j = 0; j < DIM; j++)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          dudx(i, j) += nodal_value(l, u_nodal_index) * dpsidx(l, j);
        }
      }
    }

    // Z-component of vorticity
    vorticity[0] = dudx(1, 0) - dudx(0, 1);
  }


  //==============================================================
  ///  Get integral of kinetic energy over element:
  //==============================================================
  double SphericalNavierStokesEquations::kin_energy() const
  {
    // Initialise
    double kin_en = 0.0;

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Loop over directions
      double veloc_squared = 0.0;
      for (unsigned i = 0; i < 3; i++)
      {
        veloc_squared += interpolated_u_spherical_nst(s, i) *
                         interpolated_u_spherical_nst(s, i);
      }

      kin_en += 0.5 * veloc_squared * w * J;
    }

    return kin_en;
  }


  //==========================================================================
  ///  Get integral of time derivative of kinetic energy over element:
  //==========================================================================
  double SphericalNavierStokesEquations::d_kin_energy_dt() const
  {
    throw OomphLibError("Not tested for spherical Navier-Stokes elements",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);


    // Initialise
    double d_kin_en_dt = 0.0;

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Get the number of nodes
    const unsigned n_node = this->nnode();

    // Storage for the shape function
    Shape psi(n_node);
    DShape dpsidx(n_node, 2);

    // Get the value at which the velocities are stored
    unsigned u_index[3];
    for (unsigned i = 0; i < 3; i++)
    {
      u_index[i] = this->u_index_spherical_nst(i);
    }

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the jacobian of the mapping
      double J = dshape_eulerian_at_knot(ipt, psi, dpsidx);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Now work out the velocity and the time derivative
      Vector<double> interpolated_u(3, 0.0), interpolated_dudt(3, 0.0);

      // Loop over the shape functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the velocity components
        for (unsigned i = 0; i < 3; i++)
        {
          interpolated_u[i] += nodal_value(l, u_index[i]) * psi(l);
          interpolated_dudt[i] += du_dt_spherical_nst(l, u_index[i]) * psi(l);
        }
      }

      // Get mesh velocity bit
      if (!ALE_is_disabled)
      {
        Vector<double> mesh_velocity(2, 0.0);
        DenseMatrix<double> interpolated_dudx(3, 2, 0.0);

        // Loop over nodes again
        for (unsigned l = 0; l < n_node; l++)
        {
          // Mesh can only move in the first two directions
          for (unsigned i = 0; i < 2; i++)
          {
            mesh_velocity[i] += this->dnodal_position_dt(l, i) * psi(l);
          }

          // Now du/dx bit
          for (unsigned i = 0; i < 3; i++)
          {
            for (unsigned j = 0; j < 2; j++)
            {
              interpolated_dudx(i, j) +=
                this->nodal_value(l, u_index[i]) * dpsidx(l, j);
            }
          }
        }

        // Subtract mesh velocity from du_dt
        for (unsigned i = 0; i < 3; i++)
        {
          for (unsigned k = 0; k < 2; k++)
          {
            interpolated_dudt[i] -= mesh_velocity[k] * interpolated_dudx(i, k);
          }
        }
      }


      // Loop over directions and add up u du/dt  terms
      double sum = 0.0;
      for (unsigned i = 0; i < 3; i++)
      {
        sum += interpolated_u[i] * interpolated_dudt[i];
      }

      d_kin_en_dt += sum * w * J;
    }

    return d_kin_en_dt;
  }


  //==============================================================
  /// Return pressure integrated over the element
  //==============================================================
  double SphericalNavierStokesEquations::pressure_integral() const
  {
    // Initialise
    double press_int = 0;

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get pressure
      double press = interpolated_p_spherical_nst(s);

      // Add
      press_int += press * W;
    }

    return press_int;
  }

  //==============================================================
  ///  Compute the residuals for the Navier--Stokes
  ///  equations; flag=1(or 0): do (or don't) compute the
  ///  Jacobian as well.
  //==============================================================
  void SphericalNavierStokesEquations::
    fill_in_generic_residual_contribution_spherical_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Return immediately if there are no dofs
    if (ndof() == 0) return;

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Find out how many pressure dofs there are
    const unsigned n_pres = npres_spherical_nst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[3];
    for (unsigned i = 0; i < 3; i++)
    {
      u_nodal_index[i] = u_index_spherical_nst(i);
    }

    // Set up memory for the shape and test functions
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, 2), dtestfdx(n_node, 2);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Number of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Get Physical Variables from Element
    // Reynolds number must be multiplied by the density ratio
    const double scaled_re = re() * density_ratio();
    const double scaled_re_st = re_st() * density_ratio();
    const double scaled_re_inv_ro = re_invro() * density_ratio();
    // double scaled_re_inv_fr = re_invfr()*density_ratio();
    // double visc_ratio = viscosity_ratio();
    Vector<double> G = g();

    // Integers to store the local equations and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++) s[i] = integral_pt()->knot(ipt, i);
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_spherical_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Call the pressure shape and test functions
      pshape_spherical_nst(s, psip, testp);

      // Premultiply the weights and the Jacobian
      const double W = w * J;

      // Calculate local values of the pressure and velocity components
      // Allocate
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(3, 0.0);
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> mesh_velocity(2, 0.0);
      Vector<double> dudt(3, 0.0);
      DenseMatrix<double> interpolated_dudx(3, 2, 0.0);

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
      {
        interpolated_p += p_spherical_nst(l) * psip(l);
      }


      // Calculate velocities and derivatives:

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // cache the shape functions
        double psi_ = psif(l);
        // Loop over the positions separately
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += raw_nodal_position(l, i) * psi_;
        }

        // Loop over velocity directions (three of these)
        for (unsigned i = 0; i < 3; i++)
        {
          // Get the nodal value
          const double u_value = raw_nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psi_;
          dudt[i] += du_dt_spherical_nst(l, i) * psi_;

          // Loop over derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      if (!ALE_is_disabled)
      {
        OomphLibWarning("ALE is not tested for spherical Navier Stokes",
                        "SphericalNS::fill_in...",
                        OOMPH_EXCEPTION_LOCATION);
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions (only DIM (2) of these)
          for (unsigned i = 0; i < 2; i++)
          {
            mesh_velocity[i] += this->raw_dnodal_position_dt(l, i) * psif(l);
          }
        }
      }

      // Get the user-defined body force terms
      Vector<double> body_force(3);
      get_body_force_spherical_nst(time, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      // double source = get_source_spherical_nst(time(),ipt,interpolated_x);


      // MOMENTUM EQUATIONS
      //------------------

      const double r = interpolated_x[0];
      const double r2 = r * r;
      const double sin_theta = sin(interpolated_x[1]);
      const double cos_theta = cos(interpolated_x[1]);
      const double cot_theta = cos_theta / sin_theta;

      const double u_r = interpolated_u[0];
      const double u_theta = interpolated_u[1];
      const double u_phi = interpolated_u[2];

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cannot loop over the velocity components,
        // so address each of them separately

        // FIRST: r-component momentum equation

        unsigned i = 0;

        // If it's not a boundary condition
        local_eqn = nodal_local_eqn(l, u_nodal_index[i]);
        if (local_eqn >= 0)
        {
          // Convective r-terms
          double conv = u_r * interpolated_dudx(0, 0) * r;

          // Convective theta-terms
          conv += u_theta * interpolated_dudx(0, 1);

          // Remaining convective terms
          conv -= (u_theta * u_theta + u_phi * u_phi);

          // Add the time-derivative term and the convective terms
          // to the sum
          double sum = (scaled_re_st * r2 * dudt[0] + r * scaled_re * conv) *
                       sin_theta * testf[l];

          // Subtract the mesh velocity terms
          if (!ALE_is_disabled)
          {
            sum -= scaled_re_st *
                   (mesh_velocity[0] * interpolated_dudx(0, 0) * r +
                    mesh_velocity[1] * (interpolated_dudx(0, 1) - u_theta)) *
                   r * sin_theta * testf[l];
          }

          // Add the Coriolis term
          sum -= 2.0 * scaled_re_inv_ro * sin_theta * u_phi * r2 * sin_theta *
                 testf[l];

          // r-derivative test function component of stress tensor
          sum += (-interpolated_p + 2.0 * interpolated_dudx(0, 0)) * r2 *
                 sin_theta * dtestfdx(l, 0);

          // theta-derivative test function component of stress tensor
          sum +=
            (r * interpolated_dudx(1, 0) - u_theta + interpolated_dudx(0, 1)) *
            sin_theta * dtestfdx(l, 1);

          // Undifferentiated test function component of stress tensor
          sum += 2.0 *
                 ((-r * interpolated_p + interpolated_dudx(1, 1) + 2.0 * u_r) *
                    sin_theta +
                  u_theta * cos_theta) *
                 testf(l);

          // Make the residuals negative for consistency with the
          // other Navier-Stokes equations
          residuals[local_eqn] -= sum * W;

          // CALCULATE THE JACOBIAN
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Cannot loop over the velocity components
              // --- need to compute these separately instead

              unsigned i2 = 0;

              // If at a non-zero degree of freedom add in the entry
              local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
              if (local_unknown >= 0)
              {
                // Convective r-term contribution
                double jac_conv = r * (psif(l2) * interpolated_dudx(0, 0) +
                                       u_r * dpsifdx(l2, 0));

                // Convective theta-term contribution
                jac_conv += u_theta * dpsifdx(l2, 1);

                // Add the time-derivative contribution and the convective
                // contribution to the sum
                double jac_sum =
                  (scaled_re_st * node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                     psif(l2) * r2 +
                   scaled_re * jac_conv * r) *
                  testf(l);

                // Subtract the mesh-velocity terms
                if (!ALE_is_disabled)
                {
                  jac_sum -= scaled_re_st *
                             (mesh_velocity[0] * dpsifdx(l2, 0) * r +
                              mesh_velocity[1] * dpsifdx(l2, 1)) *
                             r * sin_theta * testf(l);
                }


                // Contribution from the r-derivative test function part
                // of stress tensor
                jac_sum += 2.0 * dpsifdx(l2, 0) * dtestfdx(l, 0) * r2;

                // Contribution from the theta-derivative test function part
                // of stress tensor
                jac_sum += dpsifdx(l2, 1) * dtestfdx(l, 1);


                // Contribution from the undifferentiated test function part
                // of stress tensor
                jac_sum += 4.0 * psif[l2] * testf(l);

                // Add the total contribution to the jacobian multiplied
                // by the jacobian weight
                jacobian(local_eqn, local_unknown) -= jac_sum * sin_theta * W;

                // Mass matrix term
                if (flag == 2)
                {
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * testf[l] * r2 * sin_theta * W;
                }
              }
              // End of i2=0 section


              i2 = 1;

              // If at a non-zero degree of freedom add in the entry
              local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
              if (local_unknown >= 0)
              {
                // No time derivative contribution

                // Convective contribution
                double jac_sum = scaled_re *
                                 (interpolated_dudx(0, 1) - 2.0 * u_theta) *
                                 psif(l2) * r * sin_theta * testf(l);

                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  jac_sum += scaled_re_st * mesh_velocity[1] * psif(l2) * r *
                             sin_theta * testf(l);
                }

                // Contribution from the theta-derivative test function
                // part of stress tensor
                jac_sum +=
                  (r * dpsifdx(l2, 0) - psif(l2)) * dtestfdx(l, 1) * sin_theta;

                // Contribution from the undifferentiated test function
                // part of stress tensor
                jac_sum += 2.0 *
                           (dpsifdx(l2, 1) * sin_theta + psif(l2) * cos_theta) *
                           testf(l);

                // Add the full contribution to the jacobian matrix
                jacobian(local_eqn, local_unknown) -= jac_sum * W;

              } // End of i2=1 section


              i2 = 2;

              // If at a non-zero degree of freedom add in the entry
              local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
              if (local_unknown >= 0)
              {
                // Single convective-term contribution
                jacobian(local_eqn, local_unknown) += 2.0 * scaled_re * u_phi *
                                                      psif(l2) * r * sin_theta *
                                                      testf[l] * W;

                // Coriolis term
                jacobian(local_eqn, local_unknown) +=
                  2.0 * scaled_re_inv_ro * sin_theta * psif(l2) * r2 *
                  sin_theta * testf[l] * W;
              }
              // End of i2=2 section
            }
            // End of the l2 loop over the test functions

            // Now loop over pressure shape functions
            // This is the contribution from pressure gradient
            for (unsigned l2 = 0; l2 < n_pres; l2++)
            {
              /*If we are at a non-zero degree of freedom in the entry*/
              local_unknown = p_local_eqn(l2);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  psip(l2) * (r2 * dtestfdx(l, 0) + 2.0 * r * testf(l)) *
                  sin_theta * W;
              }
            }
          } // End of Jacobian calculation

        } // End of if not boundary condition statement
          // End of r-component momentum equation


        // SECOND: theta-component momentum equation

        i = 1;

        // IF it's not a boundary condition
        local_eqn = nodal_local_eqn(l, u_nodal_index[i]);
        if (local_eqn >= 0)
        {
          // All convective terms
          double conv = (u_r * interpolated_dudx(1, 0) * r +
                         u_theta * interpolated_dudx(1, 1) + u_r * u_theta) *
                          sin_theta -
                        u_phi * u_phi * cos_theta;

          // Add the time-derivative and convective terms to the
          // residuals sum
          double sum =
            (scaled_re_st * dudt[1] * r2 * sin_theta + scaled_re * r * conv) *
            testf(l);

          // Subtract the mesh velocity terms
          if (!ALE_is_disabled)
          {
            sum -= scaled_re_st *
                   (mesh_velocity[0] * interpolated_dudx(1, 0) * r +
                    mesh_velocity[1] * (interpolated_dudx(1, 1) + u_r)) *
                   r * sin_theta * testf(l);
          }

          // Add the Coriolis term
          sum -= 2.0 * scaled_re_inv_ro * cos_theta * u_phi * r2 * sin_theta *
                 testf[l];

          // r-derivative test function component of stress tensor
          sum +=
            (r * interpolated_dudx(1, 0) - u_theta + interpolated_dudx(0, 1)) *
            r * sin_theta * dtestfdx(l, 0);

          // theta-derivative test function component of stress tensor
          sum +=
            (-r * interpolated_p + 2.0 * interpolated_dudx(1, 1) + 2.0 * u_r) *
            dtestfdx(l, 1) * sin_theta;

          // Undifferentiated test function component of stress tensor
          sum -=
            ((r * interpolated_dudx(1, 0) - u_theta + interpolated_dudx(0, 1)) *
               sin_theta -
             (-r * interpolated_p + 2.0 * u_r + 2.0 * u_theta * cot_theta) *
               cos_theta) *
            testf(l);

          // Add the sum to the residuals,
          //(Negative for consistency)
          residuals[local_eqn] -= sum * W;

          // CALCULATE THE JACOBIAN
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Cannot loop over the velocity components
              // --- need to compute these separately instead

              unsigned i2 = 0;

              // If at a non-zero degree of freedom add in the entry
              local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
              if (local_unknown >= 0)
              {
                // No time derivative contribution

                // Convective terms
                double jac_sum = scaled_re *
                                 (r2 * interpolated_dudx(1, 0) + r * u_theta) *
                                 psif(l2) * sin_theta * testf(l);

                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  jac_sum -= scaled_re_st * mesh_velocity[1] * psif(l2) * r *
                             sin_theta * testf(l);
                }

                // Contribution from the r-derivative
                // test function part of stress tensor
                jac_sum += dpsifdx(l2, 1) * dtestfdx(l, 0) * sin_theta * r;

                // Contribution from the theta-derivative test function
                // part of stress tensor
                jac_sum += 2.0 * psif(l2) * dtestfdx(l, 1) * sin_theta;

                // Contribution from the undifferentiated test function
                // part of stress tensor
                jac_sum -=
                  (dpsifdx(l2, 1) * sin_theta - 2.0 * psif(l2) * cos_theta) *
                  testf(l);

                // Add the sum to the jacobian
                jacobian(local_eqn, local_unknown) -= jac_sum * W;
              }
              // End of i2=0 section


              i2 = 1;

              // If at a non-zero degree of freedom add in the entry
              local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
              if (local_unknown >= 0)
              {
                // Contribution from the convective terms
                double jac_conv = r * u_r * dpsifdx(l2, 0) +
                                  u_theta * dpsifdx(l2, 1) +
                                  (interpolated_dudx(1, 1) + u_r) * psif(l2);

                // Add the time-derivative term and the
                // convective terms
                double jac_sum =
                  (scaled_re_st * node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                     psif(l2) * r2 +
                   scaled_re * r * jac_conv) *
                  testf(l) * sin_theta;


                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  jac_sum -= scaled_re_st *
                             (mesh_velocity[0] * dpsifdx(l2, 0) * r +
                              mesh_velocity[1] * dpsifdx(l2, 1)) *
                             r * sin_theta * testf(l);
                }

                // Contribution from the r-derivative test function
                // part of stress tensor
                jac_sum += (r * dpsifdx(l2, 0) - psif(l2)) * r *
                           dtestfdx(l, 0) * sin_theta;

                // Contribution from the theta-derivative test function
                // part of stress tensor
                jac_sum += 2.0 * dpsifdx(l2, 1) * dtestfdx(l, 1) * sin_theta;

                // Contribution from the undifferentiated test function
                // part of stress tensor
                jac_sum -= ((r * dpsifdx(l2, 0) - psif(l2)) * sin_theta -
                            2.0 * psif(l2) * cot_theta * cos_theta) *
                           testf(l);

                // Add the contribution to the jacobian
                jacobian(local_eqn, local_unknown) -= jac_sum * W;

                // Mass matrix term
                if (flag == 2)
                {
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * testf[l] * r2 * sin_theta * W;
                }
              }
              // End of i2=1 section


              i2 = 2;

              // If at a non-zero degree of freedom add in the entry
              local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
              if (local_unknown >= 0)
              {
                // Only a convective term contribution
                jacobian(local_eqn, local_unknown) += scaled_re * 2.0 * r *
                                                      psif(l2) * u_phi *
                                                      cos_theta * testf(l) * W;

                // Coriolis term
                jacobian(local_eqn, local_unknown) +=
                  2.0 * scaled_re_inv_ro * cos_theta * psif(l2) * r2 *
                  sin_theta * testf[l] * W;
              }
              // End of i2=2 section
            }
            // End of the l2 loop over the test functions

            /*Now loop over pressure shape functions*/
            /*This is the contribution from pressure gradient*/
            for (unsigned l2 = 0; l2 < n_pres; l2++)
            {
              /*If we are at a non-zero degree of freedom in the entry*/
              local_unknown = p_local_eqn(l2);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  psip(l2) * r *
                  (dtestfdx(l, 1) * sin_theta + cos_theta * testf[l]) * W;
              }
            }
          } /*End of Jacobian calculation*/

        } // End of if not boundary condition statement
          // End of theta-component of momentum


        // THIRD: phi-component momentum equation

        i = 2;

        // IF it's not a boundary condition
        local_eqn = nodal_local_eqn(l, u_nodal_index[i]);
        if (local_eqn >= 0)
        {
          // Convective r-terms
          double conv = u_r * interpolated_dudx(2, 0) * r * sin_theta;

          // Convective theta-terms
          conv += u_theta * interpolated_dudx(2, 1) * sin_theta;

          // Remaining convective terms
          conv += u_phi * (u_r * sin_theta + u_theta * cos_theta);

          // Add the time-derivative term and the convective terms to the sum
          double sum =
            (scaled_re_st * r2 * dudt[2] * sin_theta + scaled_re * conv * r) *
            testf(l);

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            sum -= scaled_re_st *
                   (mesh_velocity[0] * interpolated_dudx(2, 0) * r +
                    mesh_velocity[1] * interpolated_dudx(2, 1)) *
                   r * sin_theta * testf(l);
          }

          // Add the Coriolis term
          sum += 2.0 * scaled_re_inv_ro *
                 (cos_theta * u_theta + sin_theta * u_r) * r2 * sin_theta *
                 testf[l];

          // r-derivative test function component of stress tensor
          sum += (r2 * interpolated_dudx(2, 0) - r * u_phi) * dtestfdx(l, 0) *
                 sin_theta;

          // theta-derivative test function component of stress tensor
          sum += (interpolated_dudx(2, 1) * sin_theta - u_phi * cos_theta) *
                 dtestfdx(l, 1);

          // Undifferentiated test function component of stress tensor
          sum -= ((r * interpolated_dudx(2, 0) - u_phi) * sin_theta +
                  (interpolated_dudx(2, 1) - u_phi * cot_theta) * cos_theta) *
                 testf(l);

          // Add the sum to the residuals
          residuals[local_eqn] -= sum * W;


          // CALCULATE THE JACOBIAN
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Cannot loop over the velocity components
              // -- need to compute these separately instead

              unsigned i2 = 0;

              // If at a non-zero degree of freedom add in the entry
              local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
              if (local_unknown >= 0)
              {
                // Contribution from convective terms
                jacobian(local_eqn, local_unknown) -=
                  scaled_re * (r * interpolated_dudx(2, 0) + u_phi) * psif(l2) *
                  testf(l) * r * sin_theta * W;

                // Coriolis term
                jacobian(local_eqn, local_unknown) -=
                  2.0 * scaled_re_inv_ro * sin_theta * psif(l2) * r2 *
                  sin_theta * testf[l] * W;
              }
              // End of i2=0 section

              i2 = 1;

              // If at a non-zero degree of freedom add in the entry
              local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
              if (local_unknown >= 0)
              {
                // Contribution from convective terms
                jacobian(local_eqn, local_unknown) -=
                  scaled_re *
                  (interpolated_dudx(2, 1) * sin_theta + u_phi * cos_theta) *
                  r * psif(l2) * testf(l) * W;


                // Coriolis term
                jacobian(local_eqn, local_unknown) -=
                  2.0 * scaled_re_inv_ro * cos_theta * psif(l2) * r2 *
                  sin_theta * testf[l] * W;
              }
              // End of i2=1 section


              i2 = 2;

              // If at a non-zero degree of freedom add in the entry
              local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
              if (local_unknown >= 0)
              {
                // Convective terms
                double jac_conv = r * u_r * dpsifdx(l2, 0) * sin_theta;

                // Convective theta-term contribution
                jac_conv += u_theta * dpsifdx(l2, 1) * sin_theta;

                // Contribution from the remaining convective terms
                jac_conv += (u_r * sin_theta + u_theta * cos_theta) * psif(l2);

                // Add the convective and time derivatives
                double jac_sum =
                  (scaled_re_st * node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                     psif(l2) * r2 * sin_theta +
                   scaled_re * r * jac_conv) *
                  testf(l);


                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  jac_sum -= scaled_re_st *
                             (mesh_velocity[0] * dpsifdx(l2, 0) * r +
                              mesh_velocity[1] * dpsifdx(l2, 1)) *
                             r * sin_theta * testf(l);
                }

                // Contribution from the r-derivative test function
                // part of stress tensor
                jac_sum += (r * dpsifdx(l2, 0) - psif(l2)) * dtestfdx(l, 0) *
                           r * sin_theta;

                // Contribution from the theta-derivative test function
                // part of stress tensor
                jac_sum += (dpsifdx(l2, 1) * sin_theta - psif(l2) * cos_theta) *
                           dtestfdx(l, 1);

                // Contribution from the undifferentiated test function
                // part of stress tensor
                jac_sum -=
                  ((r * dpsifdx(l2, 0) - psif(l2)) * sin_theta +
                   (dpsifdx(l2, 1) - psif(l2) * cot_theta) * cos_theta) *
                  testf(l);

                // Add to the jacobian
                jacobian(local_eqn, local_unknown) -= jac_sum * W;

                // Mass matrix term
                if (flag == 2)
                {
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * psif(l2) * testf[l] * r2 * sin_theta * W;
                }
              }
              // End of i2=2 section
            }
            // End of the l2 loop over the test functions

            // We assume phi-derivatives to be zero
            // No phi-contribution to the pressure gradient to include here

          } // End of Jacobian calculation

        } // End of if not boundary condition statement
          // End of phi-component of momentum

      } // End of loop over shape functions


      // CONTINUITY EQUATION
      //-------------------

      // Loop over the shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        local_eqn = p_local_eqn(l);
        // If not a boundary conditions
        // Can't loop over the velocity components, so put these in separately
        if (local_eqn >= 0)
        {
          residuals[local_eqn] += ((2.0 * u_r + r * interpolated_dudx(0, 0) +
                                    interpolated_dudx(1, 1)) *
                                     sin_theta +
                                   u_theta * cos_theta) *
                                  r * testp(l) * W;

          // CALCULATE THE JACOBIAN
          if (flag)
          {
            // Loop over the velocity shape functions
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Can't loop over velocity components, so put these in separately

              // Contribution from the r-component
              unsigned i2 = 0;
              /*If we're at a non-zero degree of freedom add it in*/
              local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);

              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  (2.0 * psif(l2) + r * dpsifdx(l2, 0)) * r * sin_theta *
                  testp(l) * W;
              }
              // End of contribution from r-component


              // Contribution from the theta-component
              i2 = 1;
              /*If we're at a non-zero degree of freedom add it in*/
              local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);

              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  r * (dpsifdx(l2, 1) * sin_theta + psif(l2) * cos_theta) *
                  testp(l) * W;
              }
              // End of contribution from theta-component

            } // End of loop over l2
          } // End of Jacobian calculation

        } // End of if not boundary condition
      } // End of loop over l
    }
  }

  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //=========================================================================
  ///  Add to the set \c paired_load_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for all values (pressures, velocities) that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  void QSphericalCrouzeixRaviartElement::identify_load_data(
    std::set<std::pair<Data*, unsigned>>& paired_load_data)
  {
    // Find the index at which the velocity is stored
    unsigned u_index[3];
    for (unsigned i = 0; i < 3; i++)
    {
      u_index[i] = this->u_index_spherical_nst(i);
    }

    // Loop over the nodes
    unsigned n_node = this->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      // Loop over the velocity components and add pointer to their data
      // and indices to the vectors
      for (unsigned i = 0; i < 3; i++)
      {
        paired_load_data.insert(std::make_pair(this->node_pt(n), u_index[i]));
      }
    }

    // Identify the pressure data
    this->identify_pressure_data(paired_load_data);
  }

  //=========================================================================
  ///  Add to the set \c paired_pressure_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for pressures values that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  void QSphericalCrouzeixRaviartElement::identify_pressure_data(
    std::set<std::pair<Data*, unsigned>>& paired_pressure_data)
  {
    // Loop over the internal data
    unsigned n_internal = this->ninternal_data();
    for (unsigned l = 0; l < n_internal; l++)
    {
      unsigned nval = this->internal_data_pt(l)->nvalue();
      // Add internal data
      for (unsigned j = 0; j < nval; j++)
      {
        paired_pressure_data.insert(
          std::make_pair(this->internal_data_pt(l), j));
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
  //=============================================================================
  void QSphericalCrouzeixRaviartElement::get_dof_numbers_for_unknowns(
    std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // number of nodes
    unsigned n_node = this->nnode();

    // number of pressure values
    unsigned n_press = this->npres_spherical_nst();

    // temporary pair (used to store dof lookup prior to being added to list)
    std::pair<unsigned, unsigned> dof_lookup;

    // pressure dof number
    unsigned pressure_dof_number = 3;

    // loop over the pressure values
    for (unsigned n = 0; n < n_press; n++)
    {
      // determine local eqn number
      int local_eqn_number = this->p_local_eqn(n);

      // ignore pinned values - far away degrees of freedom resulting
      // from hanging nodes can be ignored since these are be dealt
      // with by the element containing their master nodes
      if (local_eqn_number >= 0)
      {
        // store dof lookup in temporary pair: First entry in pair
        // is global equation number; second entry is dof type
        dof_lookup.first = this->eqn_number(local_eqn_number);
        dof_lookup.second = pressure_dof_number;

        // add to list
        dof_lookup_list.push_front(dof_lookup);
      }
    }

    // loop over the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // find the number of values at this node
      unsigned nv = this->node_pt(n)->nvalue();

      // loop over these values
      for (unsigned v = 0; v < nv; v++)
      {
        // determine local eqn number
        int local_eqn_number = this->nodal_local_eqn(n, v);

        // ignore pinned values
        if (local_eqn_number >= 0)
        {
          // store dof lookup in temporary pair: First entry in pair
          // is global equation number; second entry is dof type
          dof_lookup.first = this->eqn_number(local_eqn_number);
          dof_lookup.second = v;

          // add to list
          dof_lookup_list.push_front(dof_lookup);
        }
      }
    }
  }


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////

  // 2D Taylor--Hood
  // Set the data for the number of Variables at each node
  const unsigned QSphericalTaylorHoodElement::Initial_Nvalue[9] = {
    4, 3, 4, 3, 3, 3, 4, 3, 4};

  // Set the data for the pressure conversion array
  const unsigned QSphericalTaylorHoodElement::Pconv[4] = {0, 2, 6, 8};

  //=========================================================================
  ///  Add to the set \c paired_load_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for all values (pressures, velocities) that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  void QSphericalTaylorHoodElement::identify_load_data(
    std::set<std::pair<Data*, unsigned>>& paired_load_data)
  {
    // Find the index at which the velocity is stored
    unsigned u_index[3];
    for (unsigned i = 0; i < 3; i++)
    {
      u_index[i] = this->u_index_spherical_nst(i);
    }

    // Loop over the nodes
    unsigned n_node = this->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      // Loop over the velocity components and add pointer to their data
      // and indices to the vectors
      for (unsigned i = 0; i < 3; i++)
      {
        paired_load_data.insert(std::make_pair(this->node_pt(n), u_index[i]));
      }
    }

    // Identify the pressure data
    this->identify_pressure_data(paired_load_data);
  }

  //=========================================================================
  ///  Add to the set \c paired_load_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for pressure values that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  void QSphericalTaylorHoodElement::identify_pressure_data(
    std::set<std::pair<Data*, unsigned>>& paired_pressure_data)
  {
    // Find the index at which the pressure is stored
    unsigned p_index =
      static_cast<unsigned>(this->p_nodal_index_spherical_nst());

    // Loop over the pressure data
    unsigned n_pres = npres_spherical_nst();
    for (unsigned l = 0; l < n_pres; l++)
    {
      // The DIMth entry in each nodal data is the pressure, which
      // affects the traction
      paired_pressure_data.insert(
        std::make_pair(this->node_pt(Pconv[l]), p_index));
    }
  }


  //============================================================================
  /// Create a list of pairs for all unknowns in this element,
  /// so the first entry in each pair contains the global equation
  /// number of the unknown, while the second one contains the number
  /// of the "DOF type" that this unknown is associated with.
  /// (Function can obviously only be called if the equation numbering
  /// scheme has been set up.)
  //============================================================================
  void QSphericalTaylorHoodElement::get_dof_numbers_for_unknowns(
    std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // number of nodes
    unsigned n_node = this->nnode();

    // local eqn no for pressure unknown
    // unsigned p_index = this->p_nodal_index_spherical_nst();

    // temporary pair (used to store dof lookup prior to being added to list)
    std::pair<unsigned, unsigned> dof_lookup;

    // loop over the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // find the number of values at this node
      unsigned nv = this->node_pt(n)->nvalue();

      // loop over these values
      for (unsigned v = 0; v < nv; v++)
      {
        // determine local eqn number
        int local_eqn_number = this->nodal_local_eqn(n, v);

        // ignore pinned values - far away degrees of freedom resulting
        // from hanging nodes can be ignored since these are be dealt
        // with by the element containing their master nodes
        if (local_eqn_number >= 0)
        {
          // store dof lookup in temporary pair: Global equation number
          // is the first entry in pair
          dof_lookup.first = this->eqn_number(local_eqn_number);

          // set dof numbers: Dof number is the second entry in pair
          dof_lookup.second = v;

          // add to list
          dof_lookup_list.push_front(dof_lookup);
        }
      }
    }
  }


} // namespace oomph
