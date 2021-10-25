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
// Non-inline functions for NS elements
#include "spacetime_navier_stokes_elements.h"

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

namespace oomph
{
  /// Navier-Stokes equations static data
  template<unsigned DIM>
  Vector<double> SpaceTimeNavierStokesEquations<DIM>::Gamma(DIM, 1.0);

  //===================================================================
  /// "Magic" negative number that indicates that the pressure is
  /// not stored at a node. This cannot be -1 because that represents
  /// the positional hanging scheme in the hanging_pt object of nodes.
  //===================================================================
  template<unsigned DIM>
  int SpaceTimeNavierStokesEquations<DIM>::Pressure_not_stored_at_node = -100;

  /// Navier-Stokes equations static data
  template<unsigned DIM>
  double SpaceTimeNavierStokesEquations<DIM>::Default_Physical_Constant_Value =
    0.0;

  /// Navier-Stokes equations static data
  template<unsigned DIM>
  double SpaceTimeNavierStokesEquations<DIM>::Default_Physical_Ratio_Value =
    1.0;

  /// Navier-Stokes equations default gravity vector
  template<unsigned DIM>
  Vector<double> SpaceTimeNavierStokesEquations<DIM>::Default_Gravity_vector(
    DIM, 0.0);

  //===================================================================
  /// Compute the diagonal of the velocity/pressure mass matrices.
  /// If which one=0, both are computed, otherwise only the pressure
  /// (which_one=1) or the velocity mass matrix (which_one=2 -- the
  /// LSC version of the preconditioner only needs that one).
  //===================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::
    get_pressure_and_velocity_mass_matrix_diagonal(
      Vector<double>& press_mass_diag,
      Vector<double>& veloc_mass_diag,
      const unsigned& which_one)
  {
    // Resize and initialise
    unsigned n_dof = ndof();

    // If space needs to be assigned for the pressure mass matrix
    if ((which_one == 0) || (which_one == 1))
    {
      // Assign the appropriate amount of space
      press_mass_diag.assign(n_dof, 0.0);
    }

    // If space needs to be assigned for the velocity mass matrix
    if ((which_one == 0) || (which_one == 2))
    {
      // Assign the appropriate amount of space
      veloc_mass_diag.assign(n_dof, 0.0);
    }

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Storage for the local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Storage for the local velocity indices
    Vector<unsigned> u_nodal_index(DIM, 0.0);

    // Find the indices at which the local velocities are stored
    for (unsigned i = 0; i < DIM; i++)
    {
      // Calculate the i-th local velocity component
      u_nodal_index[i] = this->u_index_nst(i);
    }

    // Set up memory for velocity shape functions
    Shape psi(n_node);

    // Find number of pressure dofs
    unsigned n_pres = npres_nst();

    // Pressure shape function
    Shape psi_p(n_pres);

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

      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        // Calculate the i-th local coordinate at the ipt-th integration point
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Premultiply weights and Jacobian
      double W = w * J;

      // Do we want the velocity one?
      if ((which_one == 0) || (which_one == 2))
      {
        // Get the velocity shape functions
        shape_at_knot(ipt, psi);

        // Loop over the velocity shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the velocity components
          for (unsigned i = 0; i < DIM; i++)
          {
            // Get the local equation number
            local_eqn = nodal_local_eqn(l, u_nodal_index[i]);

            // If not a boundary condition
            if (local_eqn >= 0)
            {
              // Add the contribution
              veloc_mass_diag[local_eqn] += pow(psi[l], 2) * W;
            }
          } // for (unsigned i=0;i<n_dim;i++)
        } // for (unsigned l=0;l<n_node;l++)
      } // if ((which_one==0)||(which_one==2))

      // Do we want the pressure one?
      if ((which_one == 0) || (which_one == 1))
      {
        // Get the pressure shape functions
        pshape_nst(s, psi_p);

        // Loop over the veclocity shape functions
        for (unsigned l = 0; l < n_pres; l++)
        {
          // Get equation number
          local_eqn = p_local_eqn(l);

          // If not a boundary condition
          if (local_eqn >= 0)
          {
            // Add the contribution
            press_mass_diag[local_eqn] += pow(psi_p[l], 2) * W;
          }
        } // for (unsigned l=0;l<n_pres;l++)
      } // if ((which_one==0)||(which_one==1))
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of get_pressure_and_velocity_mass_matrix_diagonal


  //======================================================================
  /// Compute the vector norm of FEM solution
  //======================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::compute_norm(Vector<double>& norm)
  {
    // Resize the solution norm vector
    norm.resize(DIM + 1, 0.0);

    // Vector of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Compute the velocity norm
      for (unsigned i = 0; i < DIM; i++)
      {
        // Update the solution norm value
        norm[i] += pow(interpolated_u_nst(s, i), 2) * W;
      }

      // Update the pressure norm value
      norm[DIM] += pow(interpolated_p_nst(s), 2) * W;
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of compute_norm


  //======================================================================
  /// Validate against exact velocity solution at given time. Solution is
  /// provided via function pointer. Plot at a given number of plot points
  /// and compute L2 error and L2 norm of velocity solution over element.
  //=======================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
    double& error,
    double& norm)
  {
    // Initialise the error norm value
    error = 0.0;

    // Initialise the solution norm value
    norm = 0.0;

    // Storage for the time value
    double interpolated_t = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Vector for the spatial coordinates
    Vector<double> x(DIM, 0.0);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Exact solution Vector (u,v,[w],p)
    Vector<double> exact_soln(DIM + 1, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Assign values of x
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get x position as Vector
        x[i] = interpolated_x(s, i);
      }

      // Get the time value
      interpolated_t = interpolated_x(s, DIM);

      // Get exact solution at this point
      (*exact_soln_pt)(interpolated_t, x, exact_soln);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Velocity error
      for (unsigned i = 0; i < DIM; i++)
      {
        // Update the solution norm value
        norm += exact_soln[i] * exact_soln[i] * W;

        // Update the error norm value
        error += pow(exact_soln[i] - interpolated_u_nst(s, i), 2) * W;
      }

      // ------ DRAIG: REMOVE ----------------------------------------
      // Update the solution norm value
      norm += pow(exact_soln[DIM], 2) * W;

      // Update the error norm value
      error += pow(exact_soln[DIM] - interpolated_p_nst(s), 2) * W;
      // ------ DRAIG: REMOVE ----------------------------------------
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)

    //------------------------------------------------
    // Output the error at the appropriate plot points
    //------------------------------------------------
    // Just output at the default number of plot points
    unsigned n_plot = 5;

    // Tecplot header info
    outfile << tecplot_zone_string(n_plot);

    // How many plot points are there in total?
    unsigned num_plot_points = nplot_points(n_plot);

    // Loop over plot points
    for (unsigned i_plot = 0; i_plot < num_plot_points; i_plot++)
    {
      // Get local coordinates of plot point
      get_s_plot(i_plot, n_plot, s);

      // Loop over the spatial coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        // Assign the i-th spatial coordinate
        x[i] = interpolated_x(s, i);
      }

      // Get the time value
      interpolated_t = interpolated_x(s, DIM);

      // Get exact solution at this point
      (*exact_soln_pt)(interpolated_t, x, exact_soln);

      // Output x,y,...
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the i-th spatial coordinate
        outfile << x[i] << " ";
      }

      // Output the time value at this point
      outfile << interpolated_t << " ";

      // Output u_error,v_error
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the error in the i-th velocity component at this point
        outfile << exact_soln[i] - interpolated_u_nst(s, i) << " ";
      }

      // Output the error in the pressure field at this point
      outfile << exact_soln[DIM] - interpolated_p_nst(s) << " ";

      // End the line
      outfile << std::endl;
    } // for (unsigned i_plot=0;i_plot<num_plot_points;i_plot++)

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, n_plot);
  } // End of compute_error


  //======================================================================
  /// Validate against exact velocity solution at given time. Solution is
  /// provided via function pointer. Plot at a given number of plot points
  /// and compute L2 error and L2 norm of velocity solution over element.
  //=======================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
    Vector<double>& error,
    Vector<double>& norm)
  {
    // Resize the error norm vector
    error.resize(DIM + 1, 0.0);

    // Resize the solution norm vector
    norm.resize(DIM + 1, 0.0);

    // Storage for the time value
    double interpolated_t = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Vector for the spatial coordinates
    Vector<double> x(DIM, 0.0);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Exact solution Vector (u,v,[w],p)
    Vector<double> exact_soln(DIM + 1, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Assign values of x
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get x position as Vector
        x[i] = interpolated_x(s, i);
      }

      // Get the time value
      interpolated_t = interpolated_x(s, DIM);

      // Get exact solution at this point
      (*exact_soln_pt)(interpolated_t, x, exact_soln);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Velocity error
      for (unsigned i = 0; i < DIM; i++)
      {
        // Update the solution norm value
        norm[i] += exact_soln[i] * exact_soln[i] * W;

        // Update the error norm value
        error[i] += pow(exact_soln[i] - interpolated_u_nst(s, i), 2) * W;
      }

      // ------ DRAIG: REMOVE ----------------------------------------
      // Update the solution norm value
      norm[DIM] += pow(exact_soln[DIM], 2) * W;

      // Update the error norm value
      error[DIM] += pow(exact_soln[DIM] - interpolated_p_nst(s), 2) * W;
      // ------ DRAIG: REMOVE ----------------------------------------
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)

    //------------------------------------------------
    // Output the error at the appropriate plot points
    //------------------------------------------------
    // Just output at the default number of plot points
    unsigned n_plot = 5;

    // Tecplot header info
    outfile << tecplot_zone_string(n_plot);

    // How many plot points are there in total?
    unsigned num_plot_points = nplot_points(n_plot);

    // Loop over plot points
    for (unsigned i_plot = 0; i_plot < num_plot_points; i_plot++)
    {
      // Get local coordinates of plot point
      get_s_plot(i_plot, n_plot, s);

      // Loop over the spatial coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        // Assign the i-th spatial coordinate
        x[i] = interpolated_x(s, i);
      }

      // Get the time value
      interpolated_t = interpolated_x(s, DIM);

      // Get exact solution at this point
      (*exact_soln_pt)(interpolated_t, x, exact_soln);

      // Output x,y,...
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the i-th spatial coordinate
        outfile << x[i] << " ";
      }

      // Output the time value at this point
      outfile << interpolated_t << " ";

      // Output u_error,v_error
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the error in the i-th velocity component at this point
        outfile << exact_soln[i] - interpolated_u_nst(s, i) << " ";
      }

      // Output the error in the pressure field at this point
      outfile << exact_soln[DIM] - interpolated_p_nst(s) << " ";

      // End the line
      outfile << std::endl;
    } // for (unsigned i_plot=0;i_plot<num_plot_points;i_plot++)

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, n_plot);
  } // End of compute_error


  //======================================================================
  /// Validate against exact velocity solution at given time.
  /// Solution is provided via function pointer.
  /// Compute L2 error and L2 norm of velocity solution over element.
  //=======================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::compute_error(
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
    double& error,
    double& norm)
  {
    // Initialise the error norm value
    error = 0.0;

    // Initialise the solution norm value
    norm = 0.0;

    // Storage for the time value
    double interpolated_t = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Vector for the spatial coordinates
    Vector<double> x(DIM, 0.0);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Exact solution Vector (u,v,[w],p)
    Vector<double> exact_soln(DIM + 1, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Assign values of x
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get x position as Vector
        x[i] = interpolated_x(s, i);
      }

      // Get the time value
      interpolated_t = interpolated_x(s, DIM);

      // Get exact solution at this point
      (*exact_soln_pt)(interpolated_t, x, exact_soln);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Velocity error
      for (unsigned i = 0; i < DIM; i++)
      {
        // Update the solution norm value
        norm += exact_soln[i] * exact_soln[i] * W;

        // Update the error norm value
        error += pow(exact_soln[i] - interpolated_u_nst(s, i), 2) * W;
      }
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of compute_error

  //======================================================================
  /// Validate against exact velocity solution Solution is provided via a
  /// function pointer. Compute L2 error and L2 norm of velocity solution
  /// over element.
  //=======================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::compute_error(
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    // Initialise the error norm value
    error = 0.0;

    // Initialise the solution norm value
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Vector for the spatial coordinates
    Vector<double> x(DIM, 0.0);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Exact solution Vector (u,v,[w],p)
    Vector<double> exact_soln(DIM + 1, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Assign values of x
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get x position as Vector
        x[i] = interpolated_x(s, i);
      }

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Velocity error
      for (unsigned i = 0; i < DIM; i++)
      {
        // Update the solution norm value
        norm += exact_soln[i] * exact_soln[i] * W;

        // Update the error norm value
        error += pow(exact_soln[i] - interpolated_u_nst(s, i), 2) * W;
      }
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of compute_error

  //======================================================================
  /// Validate against exact velocity solution
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points and compute L2 error
  /// and L2 norm of velocity solution over element.
  //=======================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    // Initialise the error norm value
    error = 0.0;

    // Initialise the solution norm value
    norm = 0.0;

    // Storage for the time value
    double interpolated_t = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Vector for the spatial coordinates
    Vector<double> x(DIM, 0.0);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Output the tecplot header
    outfile << "ZONE" << std::endl;

    // Exact solution Vector (u,v,[w],p)
    Vector<double> exact_soln(DIM + 1, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Assign values of x
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get x position as Vector
        x[i] = interpolated_x(s, i);
      }

      // Get the time value
      interpolated_t = interpolated_x(s, DIM);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Velocity error
      for (unsigned i = 0; i < DIM; i++)
      {
        // Update the solution norm value
        norm += exact_soln[i] * exact_soln[i] * W;

        // Update the error norm value
        error += pow(exact_soln[i] - interpolated_u_nst(s, i), 2) * W;
      }

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the i-th coordinate value
        outfile << x[i] << " ";
      }

      // Output the time value at this point
      outfile << interpolated_t << " ";

      // Output x,y,u_error,v_error
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the error in the i-th velocity component at this point
        outfile << exact_soln[i] - interpolated_u_nst(s, i) << " ";
      }

      // Finish the line off
      outfile << std::endl;
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of compute_error

  //======================================================================
  /// Output "exact" solution
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  /// Function prints as many components as are returned in solution Vector.
  //=======================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& n_plot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Storage for the time value
    double interpolated_t = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Vector for coordinates
    Vector<double> x(DIM, 0.0);

    // Tecplot header info
    outfile << tecplot_zone_string(n_plot);

    // Exact solution Vector
    Vector<double> exact_soln;

    // How many plot points are there in total?
    unsigned num_plot_points = nplot_points(n_plot);

    // Loop over plot points
    for (unsigned i_plot = 0; i_plot < num_plot_points; i_plot++)
    {
      // Get local coordinates of plot point
      get_s_plot(i_plot, n_plot, s);

      // Loop over the spatial coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        // Assign the i-th spatial coordinate
        x[i] = interpolated_x(s, i);
      }

      // Get the time value
      interpolated_t = interpolated_x(s, DIM);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,...
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the i-th spatial coordinate
        outfile << x[i] << " ";
      }

      // Output the time value at this point
      outfile << interpolated_t << " ";

      // Output "exact solution"
      for (unsigned i = 0; i < exact_soln.size(); i++)
      {
        // Output the i-th (exact) velocity component value
        outfile << exact_soln[i] << " ";
      }

      // End the line
      outfile << std::endl;
    } // for (unsigned i_plot=0;i_plot<num_plot_points;i_plot++)

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, n_plot);
  } // End of output_fct

  //======================================================================
  /// Output "exact" solution at a given time
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  /// Function prints as many components as are returned in solution Vector.
  //=======================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& n_plot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
    // Storage for the time value
    double interpolated_t = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Vector for coordinates
    Vector<double> x(DIM, 0.0);

    // Tecplot header info
    outfile << tecplot_zone_string(n_plot);

    // Exact solution Vector
    Vector<double> exact_soln;

    // How many plot points are there in total?
    unsigned num_plot_points = nplot_points(n_plot);

    // Loop over plot points
    for (unsigned i_plot = 0; i_plot < num_plot_points; i_plot++)
    {
      // Get local coordinates of plot point
      get_s_plot(i_plot, n_plot, s);

      // Loop over the spatial coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        // Assign the i-th spatial coordinate
        x[i] = interpolated_x(s, i);
      }

      // Get the time value
      interpolated_t = interpolated_x(s, DIM);

      // Get exact solution at this point
      (*exact_soln_pt)(interpolated_t, x, exact_soln);

      // Output x,y,...
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the i-th spatial coordinate value
        outfile << x[i] << " ";
      }

      // Output the time value at this point
      outfile << interpolated_t << " ";

      // Output "exact solution"
      for (unsigned i = 0; i < exact_soln.size(); i++)
      {
        // Output the i-th (exact) velocity component value
        outfile << exact_soln[i] << " ";
      }

      // End the line
      outfile << std::endl;
    } // for (unsigned i_plot=0;i_plot<num_plot_points;i_plot++)

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, n_plot);
  } // End of output_fct

  //==============================================================
  /// Output function: Velocities only
  /// x,y,[z],u,v,[w]
  /// in tecplot format at specified previous timestep (t=0: present;
  /// t>0: previous timestep). Specified number of plot points in each
  /// coordinate direction.
  /// DRAIG: Should be broken!
  //==============================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::output_veloc(std::ostream& outfile,
                                                         const unsigned& n_plot,
                                                         const unsigned& t)
  {
    // Find number of nodes
    unsigned n_node = nnode();

    // Local shape function
    Shape psi(n_node);

    // Vectors of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Tecplot header info
    outfile << tecplot_zone_string(n_plot);

    // How many plot points are there?
    unsigned num_plot_points = nplot_points(n_plot);

    // Loop over plot points
    for (unsigned i_plot = 0; i_plot < num_plot_points; i_plot++)
    {
      // Get the local coordinates of plot point
      get_s_plot(i_plot, n_plot, s);

      // Coordinates
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        // Output the i-th spatial coordinate value
        outfile << interpolated_x(s, i) << " ";
      }

      // Velocities
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the i-th velocity component
        outfile << interpolated_u_nst(s, i) << " ";
      }

      // End the line
      outfile << std::endl;
    } // for (unsigned i_plot=0;i_plot<num_plot_points;i_plot++)

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, n_plot);
  } // End of output_veloc

  //==============================================================
  /// Output function:
  /// x,y,[z],u,v,[w],p
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::output(std::ostream& outfile,
                                                   const unsigned& n_plot)
  {
    // Find number of nodes
    unsigned n_node = nnode();

    // Local shape function
    Shape psi(n_node);

    // Vectors of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Tecplot header info
    outfile << tecplot_zone_string(n_plot);

    // How many plot points are there?
    unsigned num_plot_points = nplot_points(n_plot);

    // Loop over plot points
    for (unsigned i_plot = 0; i_plot < num_plot_points; i_plot++)
    {
      // Get the local coordinates of plot point
      get_s_plot(i_plot, n_plot, s);

      // Coordinates
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        // Output the i-th spatial coordinate value
        outfile << interpolated_x(s, i) << " ";
      }

      // Velocities
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the i-th velocity component
        outfile << interpolated_u_nst(s, i) << " ";
      }

      // Pressure
      outfile << interpolated_p_nst(s) << " ";

      // End the line
      outfile << std::endl;
    } // for (unsigned i_plot=0;i_plot<num_plot_points;i_plot++)

    // Add an extra line
    outfile << std::endl;

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, n_plot);
  } // End of output

  //==============================================================
  /// C-style output function:
  /// x,y,[z],u,v,[w],p
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::output(FILE* file_pt,
                                                   const unsigned& n_plot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Tecplot header info
    fprintf(file_pt, "%s", tecplot_zone_string(n_plot).c_str());

    // How many plot points in total?
    unsigned num_plot_points = nplot_points(n_plot);

    // Loop over plot points
    for (unsigned i_plot = 0; i_plot < num_plot_points; i_plot++)
    {
      // Get the local coordinates associated with this plot point
      get_s_plot(i_plot, n_plot, s);

      // Loop over the coordinate directions
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        // Output the i-th coordinate value to file
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }

      // Loop over the velocity components
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the i-th velocity component to file
        fprintf(file_pt, "%g ", interpolated_u_nst(s, i));
      }

      // Pressure
      fprintf(file_pt, "%g \n", interpolated_p_nst(s));
    }

    // Finish the line off
    fprintf(file_pt, "\n");

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, n_plot);
  } // End of output

  //==============================================================
  /// Full output function:
  /// x,y,t,u,v,p,du/dt,dv/dt,dissipation
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction
  //==============================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::full_output(std::ostream& outfile,
                                                        const unsigned& n_plot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Tecplot header info
    outfile << tecplot_zone_string(n_plot);

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Set up memory for the shape functions
    Shape psif(n_node);

    // Set up memory for the shape function derivatives
    DShape dpsifdx(n_node, DIM);

    // How many plot points in total?
    unsigned num_plot_points = nplot_points(n_plot);

    // Loop over plot points
    for (unsigned i_plot = 0; i_plot < num_plot_points; i_plot++)
    {
      // Get the local coordinates of plot point
      get_s_plot(i_plot, n_plot, s);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psif, dpsifdx);

      // Allocate storage for the mesh velocity
      Vector<double> mesh_velocity(DIM, 0.0);

      // Allocate storage for the acceleration
      Vector<double> du_dt(DIM, 0.0);

      // Allocate storage for the ALE acceleration
      Vector<double> du_dt_ALE(DIM, 0.0);

      // Allocate storage for the velocity derivatives
      DenseMatrix<double> interpolated_dudx(DIM, DIM, 0.0);

      //--------------------------------------
      // Calculate velocities and derivatives:
      //--------------------------------------
      // Loop over directions
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get the index at which velocity i is stored
        unsigned u_nodal_index = u_index_nst(i);

        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Get the nodal value
          double u_value = nodal_value(l, u_nodal_index);

          // Update the i-th acceleration component
          du_dt[i] += u_value * dpsifdx(l, DIM);

          // Update the i-th mesh velocity component
          mesh_velocity[i] += nodal_position(l, i) * dpsifdx(l, DIM);

          // Loop over derivative directions for velocity gradients
          for (unsigned j = 0; j < DIM; j++)
          {
            // Update the value of du_i/dx_j
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        } // for (unsigned l=0;l<n_node;l++)
      } // for (unsigned i=0;i<DIM;i++)

      //---------------------------------------------
      // Get du/dt in ALE form (incl. mesh velocity):
      //---------------------------------------------
      // Loop over the coordinate directions
      for (unsigned i = 0; i < DIM; i++)
      {
        // Store the i-th acceleration component
        du_dt_ALE[i] = du_dt[i];

        // Loop over the coordinate directions
        for (unsigned k = 0; k < DIM; k++)
        {
          // Take the mesh velocity into account
          du_dt_ALE[i] -= mesh_velocity[k] * interpolated_dudx(i, k);
        }
      } // for (unsigned i=0;i<DIM;i++)

      // Output the coordinates
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        // Output the i-th coordinate value
        outfile << interpolated_x(s, i) << " ";
      }

      // Output the velocity components
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the i-th velocity component
        outfile << interpolated_u_nst(s, i) << " ";
      }

      // Output the pressure
      outfile << interpolated_p_nst(s) << " ";

      // Output the acceleration
      for (unsigned i = 0; i < DIM; i++)
      {
        // Output the ALE acceleration term
        outfile << du_dt_ALE[i] << " ";
      }

      // Dissipation
      outfile << dissipation(s) << " ";

      // End the line
      outfile << std::endl;
    } // for (unsigned i_plot=0;i_plot<num_plot_points;i_plot++)

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, n_plot);
  } // End of full_output


  //==============================================================
  /// Output function for vorticity.
  /// x,y,[z],[omega_x,omega_y,[and/or omega_z]]
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::output_vorticity(
    std::ostream& outfile, const unsigned& n_plot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Create vorticity vector of the required size
    Vector<double> vorticity;

    // If we're in 2D the vorticity field is a scalar field
    if (DIM == 2)
    {
      // Resize the vorticity vector
      vorticity.resize(1);
    }
    // If we're in 3D the vorticity field is a vector field
    else if (DIM == 3)
    {
      // Resize the vorticity vector
      vorticity.resize(3);
    }
    // If we're in 1D
    else
    {
      std::string error_message =
        "Can't output vorticity in 1D - in fact, what ";
      error_message += "do you mean\nby the 1D Navier-Stokes equations?\n";
      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Tecplot header info
    outfile << tecplot_zone_string(n_plot);

    // How many plot points are there in total?
    unsigned num_plot_points = nplot_points(n_plot);

    // Loop over plot points
    for (unsigned i_plot = 0; i_plot < num_plot_points; i_plot++)
    {
      // Get local coordinates of plot point
      get_s_plot(i_plot, n_plot, s);

      // Coordinates
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        // Output the i-th coordinate value
        outfile << interpolated_x(s, i) << " ";
      }

      // Get vorticity
      get_vorticity(s, vorticity);

      // If we're in 2D
      if (DIM == 2)
      {
        // Output the vorticity field value
        outfile << vorticity[0];
      }
      // If we're in 3D
      else
      {
        // Output the vorticity field
        outfile << vorticity[0] << " " << vorticity[1] << " " << vorticity[2]
                << " ";
      }

      // Finish the line off
      outfile << std::endl;
    }

    // Add in an extra line
    outfile << std::endl;

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, n_plot);
  } // End of output_vorticity

  //==============================================================
  /// Return integral of dissipation over element
  //==============================================================
  template<unsigned DIM>
  double SpaceTimeNavierStokesEquations<DIM>::dissipation() const
  {
    // Initialise the elemental dissipation value
    double diss = 0.0;

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        // Calculate the i-th local coordinate value
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Storage for the strain rate matrix
      DenseMatrix<double> strainrate(DIM, DIM);

      // Get strain rate matrix
      strain_rate(s, strainrate);

      // Initialise the local dissipation
      double local_diss = 0.0;

      // Loop over the coordinate directions
      for (unsigned i = 0; i < DIM; i++)
      {
        // Loop over the coordinate directions
        for (unsigned j = 0; j < DIM; j++)
        {
          // Update the local dissipation value
          local_diss += 2.0 * strainrate(i, j) * strainrate(i, j);
        }
      } // for (unsigned i=0;i<DIM;i++)

      // Update the elemental dissipation value
      diss += local_diss * w * J;
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)

    // Return the elemental dissipation value
    return diss;
  } // End of dissipation

  //==============================================================
  /// Compute traction (on the viscous scale) exerted onto
  /// the fluid at local coordinate s. N has to be outer unit normal
  /// to the fluid.
  //==============================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::get_traction(
    const Vector<double>& s, const Vector<double>& N, Vector<double>& traction)
  {
    // Allocate space for the strain rate matrix
    DenseMatrix<double> strainrate(DIM, DIM);

    // Get velocity gradients
    strain_rate(s, strainrate);

    // Get pressure
    double press = interpolated_p_nst(s);

    // Loop over traction components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Add in the pressure contribution
      traction[i] = -press * N[i];

      // Loop over the strain rate entries
      for (unsigned j = 0; j < DIM; j++)
      {
        // Add in the strain rate contribution
        traction[i] += 2.0 * strainrate(i, j) * N[j];
      }
    } // for (unsigned i=0;i<DIM;i++)
  } // End of get_traction

  //==============================================================
  /// Compute traction (on the viscous scale) exerted onto
  /// the fluid at local coordinate s, decomposed into pressure and
  /// normal and tangential viscous components.
  /// N has to be outer unit normal to the fluid.
  //==============================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::get_traction(
    const Vector<double>& s,
    const Vector<double>& N,
    Vector<double>& traction_p,
    Vector<double>& traction_visc_n,
    Vector<double>& traction_visc_t)
  {
    // Allocate space for the traction components
    Vector<double> traction_visc(DIM);

    // Allocate space for the velocity gradients
    DenseMatrix<double> strainrate(DIM, DIM);

    // Calculate the strain rate at the local coordinate s
    strain_rate(s, strainrate);

    // Get pressure
    double press = interpolated_p_nst(s);

    // Loop over traction components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Get the pressure contribution
      traction_p[i] = -press * N[i];

      // Loop over the coordinate directions
      for (unsigned j = 0; j < DIM; j++)
      {
        // Add in the viscous stress contribution
        traction_visc[i] += 2.0 * strainrate(i, j) * N[j];
      }

      // Get the normal component of the viscous stress
      traction_visc_n[i] = traction_visc[i] * N[i];

      // Get the tangential component of the viscous stress
      traction_visc_t[i] = traction_visc[i] - traction_visc_n[i];
    } // for (unsigned i=0;i<DIM;i++)
  } // End of get_traction

  //==============================================================
  /// Return dissipation at local coordinate s
  //==============================================================
  template<unsigned DIM>
  double SpaceTimeNavierStokesEquations<DIM>::dissipation(
    const Vector<double>& s) const
  {
    // Get strain rate matrix
    DenseMatrix<double> strainrate(DIM, DIM);
    strain_rate(s, strainrate);

    // Initialise
    double local_diss = 0.0;
    for (unsigned i = 0; i < DIM; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        local_diss += 2.0 * strainrate(i, j) * strainrate(i, j);
      }
    }

    return local_diss;
  }

  //==============================================================
  /// Get strain-rate tensor: (1/2)*(du_i/dx_j+du_j/dx_i)
  //==============================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::strain_rate(
    const Vector<double>& s, DenseMatrix<double>& strainrate) const
  {
#ifdef PARANOID
    if ((strainrate.ncol() != DIM) || (strainrate.nrow() != DIM))
    {
      std::ostringstream error_message;
      error_message << "The strain rate has incorrect dimensions "
                    << strainrate.ncol() << " x " << strainrate.nrow()
                    << " not " << DIM << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Velocity gradient matrix
    DenseMatrix<double> dudx(DIM);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node);

    // Set up memory for the shape and test function derivatives
    DShape dpsidx(n_node, DIM + 1);

    // Evaluate the shape functions and shape function derivatives at this point
    dshape_eulerian(s, psi, dpsidx);

    // Initialise to zero
    dudx.initialise(0.0);

    // Loop over velocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Get the index at which the i-th velocity is stored
      unsigned u_nodal_index = u_index_nst(i);

      // Loop over derivative directions
      for (unsigned j = 0; j < DIM; j++)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Update the value of du_i/dx_j
          dudx(i, j) += nodal_value(l, u_nodal_index) * dpsidx(l, j);
        }
      } // for (unsigned j=0;j<DIM;j++)
    } // for (unsigned i=0;i<DIM;i++)

    // Loop over velocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Loop over derivative directions
      for (unsigned j = 0; j < DIM; j++)
      {
        // Calculate the (i,j)-th strain rate entry
        strainrate(i, j) = 0.5 * (dudx(i, j) + dudx(j, i));
      }
    } // for (unsigned i=0;i<DIM;i++)
  } // End of strain_rate

  //==============================================================
  /// Compute 2D vorticity vector at local coordinate s (return in
  /// one and only component of vorticity vector
  //==============================================================
  template<>
  void SpaceTimeNavierStokesEquations<2>::get_vorticity(
    const Vector<double>& s, Vector<double>& vorticity) const
  {
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
    DenseMatrix<double> dudx(DIM, DIM, 0.0);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node);

    // Set up memory for the shape and test function derivatives
    DShape dpsidx(n_node, DIM + 1);

    // Evaluate the shape functions and shape function derivatives at this point
    dshape_eulerian(s, psi, dpsidx);

    // Initialise to zero
    dudx.initialise(0.0);

    // Loop over velocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Get the index at which the i-th velocity is stored
      unsigned u_nodal_index = u_index_nst(i);

      // Loop over derivative directions
      for (unsigned j = 0; j < DIM; j++)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Update the value of du_i/dx_j
          dudx(i, j) += nodal_value(l, u_nodal_index) * dpsidx(l, j);
        }
      } // for (unsigned j=0;j<DIM;j++)
    } // for (unsigned i=0;i<DIM;i++)

    // Z-component of vorticity
    vorticity[0] = dudx(1, 0) - dudx(0, 1);
  } // End of get_vorticity

  //==============================================================
  /// Compute 2D vorticity vector at local coordinate s (return in
  /// one and only component of vorticity vector as a double
  //==============================================================
  template<>
  void SpaceTimeNavierStokesEquations<2>::get_vorticity(const Vector<double>& s,
                                                        double& vorticity) const
  {
    // Create a vector to store the vorticity
    Vector<double> vort(1, 0.0);

    // Calculate the vorticity
    get_vorticity(s, vort);

    // Assign the vorticity
    vorticity = vort[0];
  } // End of get_vorticity

  //==============================================================
  ///  Get integral of kinetic energy over element:
  /// Note that this is the "raw" kinetic energy in the sense
  /// that the density ratio has not been included. In problems
  /// with two or more fluids the user will have to remember
  /// to premultiply certain elements by the appropriate density
  /// ratio.
  //==============================================================
  template<unsigned DIM>
  double SpaceTimeNavierStokesEquations<DIM>::kin_energy() const
  {
    // Initialise the elemental kinetic energy value
    double kin_en = 0.0;

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Loop over directions
      double veloc_squared = 0.0;
      for (unsigned i = 0; i < DIM; i++)
      {
        veloc_squared += interpolated_u_nst(s, i) * interpolated_u_nst(s, i);
      }

      kin_en += 0.5 * veloc_squared * w * J;
    }

    return kin_en;

  } // End of kin_energy


  //==========================================================================
  /// Get integral of time derivative of kinetic energy over element:
  //==========================================================================
  template<unsigned DIM>
  double SpaceTimeNavierStokesEquations<DIM>::d_kin_energy_dt() const
  {
    // Initialise
    double d_kin_en_dt = 0.0;

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Get the number of nodes
    const unsigned n_node = this->nnode();

    // Storage for the shape function
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM + 1);

    // Get the value at which the velocities are stored
    unsigned u_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_index[i] = this->u_index_nst(i);
    }

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the jacobian of the mapping
      double J = dshape_eulerian_at_knot(ipt, psi, dpsidx);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Now work out the velocity and the time derivative
      Vector<double> interpolated_u(DIM, 0.0);
      Vector<double> interpolated_dudt(DIM, 0.0);

      // Loop over the shape functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the dimensions
        for (unsigned i = 0; i < DIM; i++)
        {
          interpolated_u[i] += nodal_value(l, u_index[i]) * psi(l);
          interpolated_dudt[i] += du_dt_nst(l, u_index[i]) * psi(l);
        }
      }

      // Get mesh velocity bit
      if (!ALE_is_disabled)
      {
        Vector<double> mesh_velocity(DIM, 0.0);
        DenseMatrix<double> interpolated_dudx(DIM, DIM, 0.0);

        // Loop over nodes again
        for (unsigned l = 0; l < n_node; l++)
        {
          for (unsigned i = 0; i < DIM; i++)
          {
            mesh_velocity[i] += this->dnodal_position_dt(l, i) * psi(l);

            for (unsigned j = 0; j < DIM; j++)
            {
              interpolated_dudx(i, j) +=
                this->nodal_value(l, u_index[i]) * dpsidx(l, j);
            }
          }
        }

        // Subtract mesh velocity from du_dt
        for (unsigned i = 0; i < DIM; i++)
        {
          for (unsigned k = 0; k < DIM; k++)
          {
            interpolated_dudt[i] -= mesh_velocity[k] * interpolated_dudx(i, k);
          }
        }
      }


      // Loop over directions and add up u du/dt  terms
      double sum = 0.0;
      for (unsigned i = 0; i < DIM; i++)
      {
        sum += interpolated_u[i] * interpolated_dudt[i];
      }

      d_kin_en_dt += sum * w * J;
    }

    return d_kin_en_dt;

  } // End of d_kin_energy_dt


  //==============================================================
  /// Return pressure integrated over the element
  //==============================================================
  template<unsigned DIM>
  double SpaceTimeNavierStokesEquations<DIM>::pressure_integral() const
  {
    // Initialise
    double press_int = 0;

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
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
      double press = interpolated_p_nst(s);

      // Add
      press_int += press * W;
    }

    return press_int;
  }


  //==============================================================
  /// Compute the residuals for the associated pressure advection
  /// diffusion problem. Used by the Fp preconditioner.
  /// flag=1(or 0): do (or don't) compute the Jacobian as well.
  //==============================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::
    fill_in_generic_pressure_advection_diffusion_contribution_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
  {
    OomphLibWarning("I'm not sure this is correct yet...",
                    OOMPH_CURRENT_FUNCTION,
                    OOMPH_EXCEPTION_LOCATION);

    // Return immediately if there are no dofs
    if (ndof() == 0) return;

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many pressure dofs there are
    unsigned n_pres = npres_nst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_nodal_index[i] = u_index_nst(i);
    }

    // Set up memory for the velocity shape fcts
    Shape psif(n_node);
    DShape dpsidx(n_node, DIM + 1);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres);
    Shape testp(n_pres);
    DShape dpsip(n_pres, DIM);
    DShape dtestp(n_pres, DIM);

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Get Physical Variables from Element
    // Reynolds number must be multiplied by the density ratio
    double scaled_re = re() * density_ratio();

    // Integers to store the local equations and unknowns
    int local_eqn = 0;
    int local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the velocity shape functions
      // (Derivs not needed but they are free)
      double J = this->dshape_eulerian_at_knot(ipt, psif, dpsidx);

      // Call the pressure shape and test functions
      this->dpshape_and_dptest_eulerian_nst(s, psip, dpsip, testp, dtestp);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Storage for the (Eulerian) coordinates
      Vector<double> x(DIM + 1, 0.0);

      // Calculate local values of the pressure and velocity components
      Vector<double> interpolated_u(DIM, 0.0);
      Vector<double> interpolated_dpdx(DIM, 0.0);

      // Calculate pressure gradient
      for (unsigned l = 0; l < n_pres; l++)
      {
        for (unsigned i = 0; i < DIM; i++)
        {
          interpolated_dpdx[i] += p_nst(l) * dpsip(l, i);
        }
      }

      // Loop over the velocity components
      for (unsigned i = 0; i < DIM; i++)
      {
        interpolated_u[i] = interpolated_u_nst(s, i);
      }

      // Loop over coordinate directions
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        x[i] = interpolated_x(s, i);
      }

      // Source function (for validaton only)
      double source = 0.0;
      if (Press_adv_diff_source_fct_pt != 0)
      {
        source = Press_adv_diff_source_fct_pt(x);
      }

      // Loop over the shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        local_eqn = p_local_eqn(l);

        // If not a boundary conditions
        if (local_eqn >= 0)
        {
          residuals[local_eqn] -= source * testp[l] * W;

          for (unsigned k = 0; k < DIM; k++)
          {
            residuals[local_eqn] +=
              interpolated_dpdx[k] *
              (scaled_re * interpolated_u[k] * testp[l] + dtestp(l, k)) * W;
          }

          // Jacobian too?
          if (flag)
          {
            // Loop over the shape functions
            for (unsigned l2 = 0; l2 < n_pres; l2++)
            {
              local_unknown = p_local_eqn(l2);

              // If not a boundary conditions
              if (local_unknown >= 0)
              {
                if ((int(eqn_number(local_eqn)) != Pinned_fp_pressure_eqn) &&
                    (int(eqn_number(local_unknown)) != Pinned_fp_pressure_eqn))
                {
                  for (unsigned k = 0; k < DIM; k++)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      dtestp(l2, k) *
                      (scaled_re * interpolated_u[k] * testp[l] +
                       dtestp(l, k)) *
                      W;
                  }
                }
                else
                {
                  if ((int(eqn_number(local_eqn)) == Pinned_fp_pressure_eqn) &&
                      (int(eqn_number(local_unknown)) ==
                       Pinned_fp_pressure_eqn))
                  {
                    //
                    jacobian(local_eqn, local_unknown) = 1.0;
                  }
                } // if
                  // ((int(eqn_number(local_eqn))!=Pinned_fp_pressure_eqn)&&...
              } // if (local_unknown>=0)
            } // for (unsigned l2=0;l2<n_pres;l2++)
          } // if (flag)
        } // if (local_eqn>=0)
      } // for (unsigned l=0;l<n_pres;l++)
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)

    // Now add boundary contributions from Robin BCs
    unsigned nrobin = Pressure_advection_diffusion_robin_element_pt.size();
    for (unsigned e = 0; e < nrobin; e++)
    {
      Pressure_advection_diffusion_robin_element_pt[e]
        ->fill_in_generic_residual_contribution_fp_press_adv_diff_robin_bc(
          residuals, jacobian, flag);
    }
  } // End of fill_in_generic_pressure_advection_diffusion_contribution_nst


  //==============================================================
  /// Compute the residuals for the Navier-Stokes
  /// equations; flag=1 (or 0): do (or don't) compute the
  /// Jacobian as well.
  //==============================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::
    fill_in_generic_residual_contribution_nst(Vector<double>& residuals,
                                              DenseMatrix<double>& jacobian,
                                              DenseMatrix<double>& mass_matrix,
                                              const unsigned& flag)
  {
    // Return immediately if there are no dofs
    if (ndof() == 0) return;

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many pressure dofs there are
    unsigned n_pres = npres_nst();

    // Allocate storage for the indices of the velocity components
    unsigned u_nodal_index[DIM];

    // Loop over the velocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Find the index at which the i-th local velocity is stored
      u_nodal_index[i] = u_index_nst(i);
    }

    // Set up memory for the shape functions
    Shape psif(n_node);

    // Set up memory for the test functions
    Shape testf(n_node);

    // Allocate space for the derivatives of the shape functions
    DShape dpsifdx(n_node, DIM + 1);

    // Allocate space for the derivatives of the test functions
    DShape dtestfdx(n_node, DIM + 1);

    // Set up memory for pressure shape functions
    Shape psip(n_pres);

    // Set up memory for pressure test functions
    Shape testp(n_pres);

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM + 1, 0.0);

    //-------------------------------------
    // Get physical variables from element:
    //-------------------------------------
    // Reynolds number must be multiplied by the density ratio
    double scaled_re = re() * density_ratio();

    // Get the scaled Strouhal value
    double scaled_re_st = re() * st() * density_ratio();

    // Get the scaled Strouhal value differentiated w.r.t. the Strouhal number
    double scaled_dre_st_dst = re() * density_ratio();

    // Get the scaled Reynolds / Froude number
    double scaled_re_inv_fr = re_invfr() * density_ratio();

    // Get the viscosity ratio
    double visc_ratio = viscosity_ratio();

    // Get the gravity vector
    Vector<double> G = g();

    // Integer to store the local equation number
    int local_eqn = 0;

    // Integer to store the local unknown number
    int local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        // Calculate the i-th local coordinate
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Call the pressure shape and test functions
      pshape_nst(s, psip, testp);

      // Pre-multiply the weights and the Jacobian
      double W = w * J;

      // Storage for the interpolated time value
      double interpolated_t = 0.0;

      // Storage for the interpolated pressure value
      double interpolated_p = 0.0;

      // Storage for the spatial coordinates
      Vector<double> interpolated_x(DIM, 0.0);

      // Storage for the interpolated velocity components
      Vector<double> interpolated_u(DIM, 0.0);

      // Storage for the interpolated time-derivative of the velocities
      Vector<double> interpolated_dudt(DIM, 0.0);

      // Storage for the mesh velocity
      Vector<double> mesh_velocity(DIM, 0.0);

      // Storage for the spatial derivatives of the velocity components
      DenseMatrix<double> interpolated_dudx(DIM, DIM, 0.0);

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Update the current approximation to the interpolated pressure
        interpolated_p += p_nst(l) * psip[l];
      }

      //-------------------------------------------------------------------
      // Calculate velocities, derivatives, source function and body force:
      //-------------------------------------------------------------------
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Update the interpolated time value
        interpolated_t += raw_nodal_position(l, DIM) * psif(l);

        // Loop over coordinate directions
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the nodal value
          double u_value = raw_nodal_value(l, u_nodal_index[i]);

          // Update the i-th interpolated velocity component
          interpolated_u[i] += u_value * psif[l];

          // Update the i-th interpolated coordinate value
          interpolated_x[i] += raw_nodal_position(l, i) * psif[l];

          // Update the interpolated du_i/dt value
          interpolated_dudt[i] += u_value * dpsifdx(l, DIM);

          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            // Update the interpolated du_i/dx_j value
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        } // for (unsigned i=0;i<DIM;i++)
      } // for (unsigned l=0;l<n_node;l++)

      // If ALE is enabled
      if (!ALE_is_disabled)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < DIM; i++)
          {
            // Update the i-th mesh velocity component
            mesh_velocity[i] += raw_nodal_position(l, i) * dpsifdx(l, DIM);
          }
        } // for (unsigned l=0;l<n_node;l++)
      } // if (!ALE_is_disabled)

      // Allocate space for the body force
      Vector<double> body_force(DIM, 0.0);

      // Get the user-defined body force term
      get_body_force_nst(interpolated_t, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      double source = get_source_nst(interpolated_t, ipt, interpolated_x);

      //---------------------------------
      // Assemble residuals and Jacobian:
      //---------------------------------
      //--------------------
      // Momentum equations:
      //--------------------
      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the velocity components
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the local equation number associated with this unknown and node
          local_eqn = nodal_local_eqn(l, u_nodal_index[i]);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            // Add the user-defined body force terms
            residuals[local_eqn] += body_force[i] * testf[l] * W;

            // Add the gravitational body force term
            residuals[local_eqn] += scaled_re_inv_fr * testf[l] * G[i] * W;

            // Add the pressure gradient term
            residuals[local_eqn] += interpolated_p * dtestfdx(l, i) * W;

            // Add in the contribution from the time derivative
            residuals[local_eqn] -=
              scaled_re_st * interpolated_dudt[i] * testf[l] * W;

            // If ALE is enabled
            if (!ALE_is_disabled)
            {
              // Loop over the coordinate directions
              for (unsigned k = 0; k < DIM; k++)
              {
                // Add in the mesh velocity contribution
                residuals[local_eqn] +=
                  (scaled_re_st * mesh_velocity[k] * interpolated_dudx(i, k) *
                   testf[l] * W);
              }
            } // if (!ALE_is_disabled)

            // Loop over the coordinate directions
            for (unsigned k = 0; k < DIM; k++)
            {
              // Add in the convective term contribution
              residuals[local_eqn] -= (scaled_re * interpolated_u[k] *
                                       interpolated_dudx(i, k) * testf[l] * W);
            }

            // Loop over the coordinate directions
            for (unsigned k = 0; k < DIM; k++)
            {
              // Add in the stress tensor terms:
              // NOTE: The viscosity ratio needs to go in here to ensure
              // continuity of normal stress is satisfied even in flows
              // with zero pressure gradient!
              residuals[local_eqn] -= ((interpolated_dudx(i, k) +
                                        Gamma[i] * interpolated_dudx(k, i)) *
                                       visc_ratio * dtestfdx(l, k) * W);
            }

            //------------------------
            // Calculate the Jacobian:
            //------------------------
            // If we also need to construct the Jacobian
            if (flag)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Loop over the velocity components again
                for (unsigned i2 = 0; i2 < DIM; i2++)
                {
                  // Get the local equation number associated with this unknown
                  local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);

                  // If we're at a proper degree of freedom
                  if (local_unknown >= 0)
                  {
                    // Add the contribution to the elemental matrix
                    jacobian(local_eqn, local_unknown) -=
                      (visc_ratio * Gamma[i] * dpsifdx(l2, i) *
                       dtestfdx(l, i2) * W);

                    // Now add in the inertial terms
                    jacobian(local_eqn, local_unknown) -=
                      (scaled_re * psif[l2] * interpolated_dudx(i, i2) *
                       testf[l] * W);

                    // Extra component if i2=i
                    if (i2 == i)
                    {
                      // If we also need to construct the mass matrix (only
                      // diagonal entries)
                      if (flag == 2)
                      {
                        // NOTE: This is positive because the mass matrix is
                        // taken to the other side of the equation when
                        // formulating the generalised eigenproblem.
                        mass_matrix(local_eqn, local_unknown) +=
                          (scaled_re_st * psif[l2] * testf[l] * W);
                      }

                      // Add in the time-derivative contribution
                      jacobian(local_eqn, local_unknown) -=
                        (scaled_re_st * dpsifdx(l2, DIM) * testf[l] * W);

                      // If ALE is enabled
                      if (!ALE_is_disabled)
                      {
                        // Loop over the velocity components
                        for (unsigned k = 0; k < DIM; k++)
                        {
                          // Add in the mesh velocity contribution
                          jacobian(local_eqn, local_unknown) +=
                            (scaled_re_st * mesh_velocity[k] * dpsifdx(l2, k) *
                             testf[l] * W);
                        }
                      } // if (!ALE_is_disabled)

                      // Loop over the velocity components
                      for (unsigned k = 0; k < DIM; k++)
                      {
                        // Add in the convective term contribution
                        jacobian(local_eqn, local_unknown) -=
                          (scaled_re * interpolated_u[k] * dpsifdx(l2, k) *
                           testf[l] * W);
                      }

                      // Loop over the velocity components
                      for (unsigned k = 0; k < DIM; k++)
                      {
                        // Add in the velocity gradient terms
                        jacobian(local_eqn, local_unknown) -=
                          (visc_ratio * dpsifdx(l2, k) * dtestfdx(l, k) * W);
                      }
                    } // if (i2==i)
                  } // if (local_unknown>=0)
                } // for (unsigned i2=0;i2<DIM;i2++)
              } // for (unsigned l2=0;l2<n_node;l2++)

              // Loop over pressure shape functions
              for (unsigned l2 = 0; l2 < n_pres; l2++)
              {
                // Get the local equation number associated with this degree of
                // freedom
                local_unknown = p_local_eqn(l2);

                // If we are at a proper degree of freedom
                if (local_unknown >= 0)
                {
                  // Add in the pressure gradient contribution
                  jacobian(local_eqn, local_unknown) +=
                    psip[l2] * dtestfdx(l, i) * W;
                }
              } // for (unsigned l2=0;l2<n_pres;l2++)

              //------------------------------------
              // Calculate external data information
              //------------------------------------
              // If we're storing the Strouhal number as external data then we
              // need to calculate dr/d(St) in the elemental Jacobian
              if (Strouhal_is_stored_as_external_data)
              {
                // The index of the external data (there's only one!)
                unsigned data_index = 0;

                // The index of the unknown value stored in the external data
                unsigned value_index = 0;

                // Get the local equation number associated with the extra
                // unknown
                local_unknown =
                  this->external_local_eqn(data_index, value_index);

                // If we're at a non-zero degree of freedom add in the entry
                if (local_unknown >= 0)
                {
                  // Add in the contribution from the time derivative
                  jacobian(local_eqn, local_unknown) -=
                    (scaled_dre_st_dst * interpolated_dudt[i] * testf[l] * W);

                  // If ALE is enabled
                  if (!ALE_is_disabled)
                  {
                    // Loop over the coordinate directions
                    for (unsigned k = 0; k < DIM; k++)
                    {
                      // Add in the mesh velocity contribution
                      jacobian(local_eqn, local_unknown) +=
                        (scaled_dre_st_dst * mesh_velocity[k] *
                         interpolated_dudx(i, k) * testf[l] * W);
                    }
                  } // if (!ALE_is_disabled)
                } // if (local_unknown>=0)
              } // if (Strouhal_is_stored_as_external_data)
            } // if (flag)
          } // if (local_eqn>=0)
        } // for (unsigned i=0;i<DIM;i++)
      } // for (unsigned l=0;l<n_node;l++)

      //---------------------
      // Continuity equation:
      //---------------------
      // Loop over the shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Get the local equation number associated with the pressure dof
        local_eqn = p_local_eqn(l);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Add in the source term contribution
          residuals[local_eqn] -= source * testp[l] * W;

          // Loop over velocity components
          for (unsigned k = 0; k < DIM; k++)
          {
            // Add in the velocity gradient terms
            residuals[local_eqn] += interpolated_dudx(k, k) * testp[l] * W;
          }

          //------------------------
          // Calculate the Jacobian:
          //------------------------
          // If we also need to construct the Jacobian
          if (flag)
          {
            // Loop over the velocity shape functions
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Loop over velocity components
              for (unsigned i2 = 0; i2 < DIM; i2++)
              {
                // Get the local equation number associated with this node
                local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);

                // If we're at a proper degree of freedom
                if (local_unknown >= 0)
                {
                  // Add in the velocity gradient contribution
                  jacobian(local_eqn, local_unknown) +=
                    dpsifdx(l2, i2) * testp[l] * W;
                }
              } // for (unsigned i2=0;i2<DIM;i2++)
            } // for (unsigned l2=0;l2<n_node;l2++)
          } // if (flag)
        } // if (local_eqn>=0)
      } // for (unsigned l=0;l<n_pres;l++)
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of fill_in_generic_residual_contribution_nst

  //================================================================
  /// Compute the derivatives of the residuals for the Navier-Stokes
  /// equations with respect to a parameter;
  /// flag=2 or 1 or 0: do (or don't) compute the Jacobian and
  /// mass matrix as well
  //================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::
    fill_in_generic_dresidual_contribution_nst(
      double* const& parameter_pt,
      Vector<double>& dres_dparam,
      DenseMatrix<double>& djac_dparam,
      DenseMatrix<double>& dmass_matrix_dparam,
      const unsigned& flag)
  {
    // Throw an error
    throw OomphLibError("Not yet implemented\n",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  } // End of fill_in_generic_dresidual_contribution_nst

  //==================================================================
  /// Compute the hessian tensor vector products required to
  /// perform continuation of bifurcations analytically
  //==================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::
    fill_in_contribution_to_hessian_vector_products(
      Vector<double> const& Y,
      DenseMatrix<double> const& C,
      DenseMatrix<double>& product)
  {
    // Throw an error
    throw OomphLibError("Not yet implemented\n",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  } // End of fill_in_contribution_to_hessian_vector_products

  //======================================================================
  /// Compute derivatives of elemental residual vector with respect
  /// to nodal coordinates.
  /// dresidual_dnodal_coordinates(l,i,j)=d res(l) / dX_{ij}
  /// Overloads the FD-based version in the FE base class.
  /// DRAIG: This needs doing carefully if the ALE nodes aren't fixed!!!
  //======================================================================
  template<unsigned DIM>
  void SpaceTimeNavierStokesEquations<DIM>::get_dresidual_dnodal_coordinates(
    RankThreeTensor<double>& dresidual_dnodal_coordinates)
  {
    // Throw a warning
    throw OomphLibError("Space-time update needs to be checked!",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);

    // Return immediately if there are no dofs
    if (ndof() == 0)
    {
      return;
    }

    // Determine number of nodes in element
    const unsigned n_node = nnode();

    // Determine number of pressure dofs in element
    const unsigned n_pres = npres_nst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_nodal_index[i] = u_index_nst(i);
    }

    // Set up memory for the shape and test functions
    Shape psif(n_node);
    Shape testf(n_node);
    DShape dpsifdx(n_node, DIM + 1);
    DShape dtestfdx(n_node, DIM + 1);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Deriatives of shape fct derivatives w.r.t. nodal coords
    RankFourTensor<double> d_dpsifdx_dX(DIM, n_node, n_node, DIM);
    RankFourTensor<double> d_dtestfdx_dX(DIM, n_node, n_node, DIM);

    // Derivative of Jacobian of mapping w.r.t. to nodal coords
    DenseMatrix<double> dJ_dX(DIM, n_node);

    // Derivatives of derivative of u w.r.t. nodal coords
    RankFourTensor<double> d_dudx_dX(DIM, n_node, DIM, DIM);

    // Derivatives of nodal velocities w.r.t. nodal coords:
    // Assumption: Interaction only local via no-slip so
    // X_ij only affects U_ij.
    DenseMatrix<double> d_U_dX(DIM, n_node, 0.0);

    // Determine the number of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Vector to hold local coordinates
    Vector<double> s(DIM + 1, 0.0);

    // Get physical variables from element
    // (Reynolds number must be multiplied by the density ratio)
    double scaled_re = re() * density_ratio();
    double scaled_re_st = re_st() * density_ratio();
    double scaled_re_inv_fr = re_invfr() * density_ratio();
    double visc_ratio = viscosity_ratio();
    Vector<double> G = g();

    // FD step
    double eps_fd = GeneralisedElement::Default_fd_jacobian_step;

    // Pre-compute derivatives of nodal velocities w.r.t. nodal coords:
    // Assumption: Interaction only local via no-slip so
    // X_ij only affects U_ij.
    bool element_has_node_with_aux_node_update_fct = false;
    for (unsigned q = 0; q < n_node; q++)
    {
      Node* nod_pt = node_pt(q);

      // Only compute if there's a node-update fct involved
      if (nod_pt->has_auxiliary_node_update_fct_pt())
      {
        element_has_node_with_aux_node_update_fct = true;

        // Current nodal velocity
        Vector<double> u_ref(DIM);
        for (unsigned i = 0; i < DIM; i++)
        {
          u_ref[i] = *(nod_pt->value_pt(u_nodal_index[i]));
        }

        // FD
        for (unsigned p = 0; p < DIM; p++)
        {
          // Make backup
          double backup = nod_pt->x(p);

          // Do FD step. No node update required as we're
          // attacking the coordinate directly...
          nod_pt->x(p) += eps_fd;

          // Do auxiliary node update (to apply no-slip)
          nod_pt->perform_auxiliary_node_update_fct();

          // Evaluate
          d_U_dX(p, q) =
            (*(nod_pt->value_pt(u_nodal_index[p])) - u_ref[p]) / eps_fd;

          // Reset
          nod_pt->x(p) = backup;

          // Do auxiliary node update (to apply no slip)
          nod_pt->perform_auxiliary_node_update_fct();
        }
      }
    }

    // Integer to store the local equation number
    int local_eqn = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      const double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      const double J = dshape_and_dtest_eulerian_at_knot_nst(ipt,
                                                             psif,
                                                             dpsifdx,
                                                             d_dpsifdx_dX,
                                                             testf,
                                                             dtestfdx,
                                                             d_dtestfdx_dX,
                                                             dJ_dX);

      // Call the pressure shape and test functions
      pshape_nst(s, psip, testp);

      // Storage for the interpolated time value
      double interpolated_t = 0.0;

      // Storage for the interpolated pressure value
      double interpolated_p = 0.0;

      // Storage for the spatial coordinates
      Vector<double> interpolated_x(DIM, 0.0);

      // Storage for the interpolated velocity components
      Vector<double> interpolated_u(DIM, 0.0);

      // Storage for the interpolated time-derivative of the velocities
      Vector<double> interpolated_dudt(DIM, 0.0);

      // Storage for the mesh velocity
      Vector<double> mesh_velocity(DIM, 0.0);

      // Storage for the spatial derivatives of the velocity components
      DenseMatrix<double> interpolated_dudx(DIM, DIM, 0.0);

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Update the current approximation to the interpolated pressure
        interpolated_p += p_nst(l) * psip[l];
      }

      //-------------------------------------------------------------------
      // Calculate velocities, derivatives, source function and body force:
      //-------------------------------------------------------------------
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Update the interpolated time value
        interpolated_t += raw_nodal_position(l, DIM) * psif(l);

        // Loop over coordinate directions
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the nodal value
          double u_value = raw_nodal_value(l, u_nodal_index[i]);

          // Update the i-th interpolated velocity component
          interpolated_u[i] += u_value * psif[l];

          // Update the i-th interpolated coordinate value
          interpolated_x[i] += raw_nodal_position(l, i) * psif[l];

          // Update the interpolated du_i/dt value
          interpolated_dudt[i] += u_value * dpsifdx(l, DIM);

          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            // Update the interpolated du_i/dx_j value
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        } // for (unsigned i=0;i<DIM;i++)
      } // for (unsigned l=0;l<n_node;l++)

      // If ALE is enabled
      if (!ALE_is_disabled)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < DIM; i++)
          {
            // Update the i-th mesh velocity component
            mesh_velocity[i] += raw_nodal_position(l, i) * dpsifdx(l, DIM);
          }
        } // for (unsigned l=0;l<n_node;l++)
      } // if (!ALE_is_disabled)

      // Calculate derivative of du_i/dx_k w.r.t. nodal positions X_{pq}
      // DRAIG: CHECK
      for (unsigned q = 0; q < n_node; q++)
      {
        // Loop over coordinate directions
        for (unsigned p = 0; p < DIM; p++)
        {
          for (unsigned i = 0; i < DIM; i++)
          {
            for (unsigned k = 0; k < DIM; k++)
            {
              double aux = 0.0;

              for (unsigned j = 0; j < n_node; j++)
              {
                aux += raw_nodal_value(j, u_nodal_index[i]) *
                       d_dpsifdx_dX(p, q, j, k);
              }

              //
              d_dudx_dX(p, q, i, k) = aux;
            }
          }
        }
      }

      // Allocate space for the body force
      Vector<double> body_force(DIM);

      // Get the user-defined body force term
      get_body_force_nst(interpolated_t, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      const double source = get_source_nst(interpolated_t, ipt, interpolated_x);

      // Note: Can use raw values and avoid bypassing hanging information
      // because this is the non-refineable version!

      // Allocate space for the gradient of the body force function
      DenseMatrix<double> d_body_force_dx(DIM, DIM, 0.0);

      // Get gradient of body force function
      get_body_force_gradient_nst(
        interpolated_t, ipt, s, interpolated_x, d_body_force_dx);

      // Allocate space for the gradient of the source function
      Vector<double> source_gradient(DIM, 0.0);

      // Get gradient of source function
      get_source_gradient_nst(
        interpolated_t, ipt, interpolated_x, source_gradient);

      // Get weight of actual nodal position in computation of mesh velocity
      // from positional time stepper
      const double pos_time_weight =
        (node_pt(0)->position_time_stepper_pt()->weight(1, 0));

      // Get weight of actual nodal value in computation of mesh velocity
      // from value time stepper
      const double val_time_weight =
        node_pt(0)->time_stepper_pt()->weight(1, 0);

      //---------------------------------
      // Assemble residuals and Jacobian:
      //---------------------------------
      //--------------------
      // Momentum equations:
      //--------------------
      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over coordinate directions
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the local equation number associated with this degree of
          // freedom
          local_eqn = nodal_local_eqn(l, u_nodal_index[i]);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            // Loop over coordinate directions
            for (unsigned p = 0; p < DIM; p++)
            {
              // Loop over nodes
              for (unsigned q = 0; q < n_node; q++)
              {
                //-----------------------------------
                // Residual x derivative of Jacobian:
                //-----------------------------------
                // Add the user-defined body force terms
                double sum = body_force[i] * testf[l];

                // Add the gravitational body force term
                sum += scaled_re_inv_fr * testf[l] * G[i];

                // Add the pressure gradient term
                sum += interpolated_p * dtestfdx(l, i);

                // Loop over the coordinate directions
                for (unsigned k = 0; k < DIM; k++)
                {
                  // Add in the stress tensor contribution
                  // NOTE: The viscosity ratio needs to go in here to ensure
                  // continuity of normal stress is satisfied even in flows
                  // with zero pressure gradient!
                  sum -= (visc_ratio *
                          (interpolated_dudx(i, k) +
                           Gamma[i] * interpolated_dudx(k, i)) *
                          dtestfdx(l, k));
                }

                // Add in the contribution from the time derivative
                sum -= scaled_re_st * interpolated_dudt[i] * testf[l];

                // Loop over the coordinate directions
                for (unsigned k = 0; k < DIM; k++)
                {
                  // Add in convective term contribution
                  sum -= scaled_re * interpolated_u[k] *
                         interpolated_dudx(i, k) * testf[l];
                }

                // If ALE is enabled
                if (!ALE_is_disabled)
                {
                  // Loop over the coordinate directions
                  for (unsigned k = 0; k < DIM; k++)
                  {
                    // Add in the mesh velocity contribution
                    sum += scaled_re_st * mesh_velocity[k] *
                           interpolated_dudx(i, k) * testf[l];
                  }
                }

                // Multiply through by derivative of Jacobian and integration
                // weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  sum * dJ_dX(p, q) * w;

                //----------------------------------
                // Derivative of residual x Jacobian
                //----------------------------------
                // Body force
                sum = d_body_force_dx(i, p) * psif(q) * testf(l);

                // Pressure gradient term
                sum += interpolated_p * d_dtestfdx_dX(p, q, l, i);

                // Viscous term
                for (unsigned k = 0; k < DIM; k++)
                {
                  sum -= (visc_ratio * ((interpolated_dudx(i, k) +
                                         Gamma[i] * interpolated_dudx(k, i)) *
                                          d_dtestfdx_dX(p, q, l, k) +
                                        (d_dudx_dX(p, q, i, k) +
                                         Gamma[i] * d_dudx_dX(p, q, k, i)) *
                                          dtestfdx(l, k)));
                }

                // Convective terms, including mesh velocity
                for (unsigned k = 0; k < DIM; k++)
                {
                  double tmp = scaled_re * interpolated_u[k];
                  if (!ALE_is_disabled)
                  {
                    tmp -= scaled_re_st * mesh_velocity[k];
                  }
                  sum -= tmp * d_dudx_dX(p, q, i, k) * testf(l);
                }

                if (!ALE_is_disabled)
                {
                  sum += scaled_re_st * pos_time_weight * psif(q) *
                         interpolated_dudx(i, p) * testf(l);
                }

                // Multiply through by Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) += sum * J * w;

                // Derivs w.r.t. to nodal velocities
                // ---------------------------------
                if (element_has_node_with_aux_node_update_fct)
                {
                  sum =
                    -visc_ratio * Gamma[i] * dpsifdx(q, i) * dtestfdx(l, p) -
                    scaled_re * psif(q) * interpolated_dudx(i, p) * testf(l);
                  if (i == p)
                  {
                    sum -= scaled_re_st * val_time_weight * psif(q) * testf(l);
                    for (unsigned k = 0; k < DIM; k++)
                    {
                      sum -= visc_ratio * dpsifdx(q, k) * dtestfdx(l, k);
                      double tmp = scaled_re * interpolated_u[k];
                      if (!ALE_is_disabled)
                        tmp -= scaled_re_st * mesh_velocity[k];
                      sum -= tmp * dpsifdx(q, k) * testf(l);
                    }
                  }
                  dresidual_dnodal_coordinates(local_eqn, p, q) +=
                    sum * d_U_dX(p, q) * J * w;
                }
              }
            }
          }
        }
      } // End of loop over test functions


      // CONTINUITY EQUATION
      // -------------------

      // Loop over the shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        local_eqn = p_local_eqn(l);

        // If not a boundary conditions
        if (local_eqn >= 0)
        {
          // Loop over coordinate directions
          for (unsigned p = 0; p < DIM; p++)
          {
            // Loop over nodes
            for (unsigned q = 0; q < n_node; q++)
            {
              // Residual x deriv of Jacobian
              //-----------------------------

              // Source term
              double aux = -source;

              // Loop over velocity components
              for (unsigned k = 0; k < DIM; k++)
              {
                aux += interpolated_dudx(k, k);
              }

              // Multiply through by deriv of Jacobian and integration weight
              dresidual_dnodal_coordinates(local_eqn, p, q) +=
                aux * dJ_dX(p, q) * testp[l] * w;

              // Derivative of residual x Jacobian
              //----------------------------------

              // Loop over velocity components
              aux = -source_gradient[p] * psif(q);
              for (unsigned k = 0; k < DIM; k++)
              {
                aux += d_dudx_dX(p, q, k, k);
              }
              // Multiply through by Jacobian and integration weight
              dresidual_dnodal_coordinates(local_eqn, p, q) +=
                aux * testp[l] * J * w;

              // Derivs w.r.t. to nodal velocities
              //---------------------------------
              if (element_has_node_with_aux_node_update_fct)
              {
                aux = dpsifdx(q, p) * testp(l);
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  aux * d_U_dX(p, q) * J * w;
              }
            }
          }
        }
      }
    } // End of loop over integration points
  } // End of get_dresidual_dnodal_coordinates

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  //---------------------------------
  // Space-time Taylor-Hood elements:
  //---------------------------------
  // Set the data for the number of Variables at each node
  template<>
  const unsigned QTaylorHoodSpaceTimeElement<2>::Initial_Nvalue[27] = {
    3, 2, 3, 2, 2, 2, 3, 2, 3, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 3, 2, 3};

  // Set the data for the pressure conversion array
  template<>
  const unsigned QTaylorHoodSpaceTimeElement<2>::Pconv[8] = {
    0, 2, 6, 8, 18, 20, 24, 26};

  //=========================================================================
  /// Add to the set \c paired_load_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for all values (pressures, velocities) that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  template<unsigned DIM>
  void QTaylorHoodSpaceTimeElement<DIM>::identify_load_data(
    std::set<std::pair<Data*, unsigned>>& paired_load_data)
  {
    // Allocate storage for the indices of the velocity components
    unsigned u_index[DIM];

    // Loop over the velocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Get the index at which the i-th velocity component is stored
      u_index[i] = this->u_index_nst(i);
    }

    // Get the number of nodes in this element
    unsigned n_node = this->nnode();

    // Loop over the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // Loop over the velocity components and add pointer to their data
      // and indices to the vectors
      for (unsigned i = 0; i < DIM; i++)
      {
        // Add in the node and equation number pair
        paired_load_data.insert(std::make_pair(this->node_pt(n), u_index[i]));
      }
    } // for (unsigned n=0;n<n_node;n++)

    // Identify the pressure data
    this->identify_pressure_data(paired_load_data);
  } // End of identify_load_data

  //=========================================================================
  ///  Add to the set \c paired_pressure_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for pressure values that affect the
  /// load computed in the \c get_load(...) function.,
  //=========================================================================
  template<unsigned DIM>
  void QTaylorHoodSpaceTimeElement<DIM>::identify_pressure_data(
    std::set<std::pair<Data*, unsigned>>& paired_pressure_data)
  {
    // Find the index at which the pressure is stored
    unsigned p_index = static_cast<unsigned>(this->p_nodal_index_nst());

    // Get the number of pressure degrees of freedom
    unsigned n_pres = npres_nst();

    // Loop over the pressure data
    for (unsigned l = 0; l < n_pres; l++)
    {
      // The DIM-th entry in each nodal data is the pressure, which affects
      // the traction
      paired_pressure_data.insert(
        std::make_pair(this->node_pt(Pconv[l]), p_index));
    }
  } // End of identify_pressure_data

  /*
  //============================================================================
  /// Create a list of pairs for all unknowns in this element,
  /// so the first entry in each pair contains the global equation
  /// number of the unknown, while the second one contains the number
  /// of the "DOF type" that this unknown is associated with.
  /// (Function can obviously only be called if the equation numbering
  /// scheme has been set up.)
  //============================================================================
  template<unsigned DIM>
  void QTaylorHoodSpaceTimeElement<DIM>::get_dof_numbers_for_unknowns(
   std::list<std::pair<unsigned long, unsigned> >& dof_lookup_list) const
  {
   // Get the number of nodes in this element
   unsigned n_node=this->nnode();

   // Temporary pair (used to store dof lookup prior to being added to list)
   std::pair<unsigned,unsigned> dof_lookup;

   // Loop over the nodes
   for (unsigned n=0;n<n_node;n++)
   {
    // Find the number of values at this node
    unsigned n_value=this->required_nvalue(n);

    // Loop over these values
    for (unsigned i_value=0;i_value<n_value;i_value++)
    {
     // Determine the local eqn number associated with this value
     int local_eqn_number=this->nodal_local_eqn(n,i_value);

     // Ignore pinned values - far away degrees of freedom resulting
     // from hanging nodes can be ignored since these are be dealt
     // with by the element containing their master nodes
     if (local_eqn_number>=0)
     {
      // Store dof lookup in temporary pair; the global equation number
      // is the first entry in pair
      dof_lookup.first=this->eqn_number(local_eqn_number);

      // Set dof numbers; dof number is the second entry in pair
      // DRAIG: Uncomment whichever one you want to use. Setting
      // all dof numbers to 0 means you don't distinguish between
      // velocity and pressure component; just aggregate everything.
      // In contrast, setting the dof number to i_value means the
      // dofs associated with the first velocity component, second
      // velocity component and the pressure are all separated
      dof_lookup.second=i_value;

      // Add to list
      dof_lookup_list.push_front(dof_lookup);
     }
    } // for (unsigned v=0;v<nv;v++)
   } // for (unsigned n=0;n<n_node;n++)
  } // End of get_dof_numbers_for_unknowns
  */

  //====================================================================
  /// Force build of templates
  //====================================================================
  template class SpaceTimeNavierStokesEquations<2>;
  template class QTaylorHoodSpaceTimeElement<2>;
} // namespace oomph
