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

#include "navier_stokes_elements.h"


namespace oomph
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


  /// Navier--Stokes equations static data
  template<unsigned DIM>
  Vector<double> NavierStokesEquations<DIM>::Gamma(DIM, 1.0);

  //=================================================================
  /// "Magic" negative number that indicates that the pressure is
  /// not stored at a node. This cannot be -1 because that represents
  /// the positional hanging scheme in the hanging_pt object of nodes
  //=================================================================
  template<unsigned DIM>
  int NavierStokesEquations<DIM>::Pressure_not_stored_at_node = -100;

  /// Navier--Stokes equations static data
  template<unsigned DIM>
  double NavierStokesEquations<DIM>::Default_Physical_Constant_Value = 0.0;

  /// Navier--Stokes equations static data
  template<unsigned DIM>
  double NavierStokesEquations<DIM>::Default_Physical_Ratio_Value = 1.0;

  /// Navier-Stokes equations default gravity vector
  template<unsigned DIM>
  Vector<double> NavierStokesEquations<DIM>::Default_Gravity_vector(DIM, 0.0);


  //===================================================================
  /// Compute the diagonal of the velocity/pressure mass matrices.
  /// If which one=0, both are computed, otherwise only the pressure
  /// (which_one=1) or the velocity mass matrix (which_one=2 -- the
  /// LSC version of the preconditioner only needs that one)
  //===================================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::
    get_pressure_and_velocity_mass_matrix_diagonal(
      Vector<double>& press_mass_diag,
      Vector<double>& veloc_mass_diag,
      const unsigned& which_one)
  {
    // Resize and initialise
    unsigned n_dof = ndof();

    if ((which_one == 0) || (which_one == 1))
      press_mass_diag.assign(n_dof, 0.0);
    if ((which_one == 0) || (which_one == 2))
      veloc_mass_diag.assign(n_dof, 0.0);

    // find out how many nodes there are
    unsigned n_node = nnode();

    // find number of spatial dimensions
    unsigned n_dim = this->dim();

    // Local coordinates
    Vector<double> s(n_dim);

    // find the indices at which the local velocities are stored
    Vector<unsigned> u_nodal_index(n_dim);
    for (unsigned i = 0; i < n_dim; i++)
    {
      u_nodal_index[i] = this->u_index_nst(i);
    }

    // Set up memory for veloc shape functions
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
      for (unsigned i = 0; i < n_dim; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Premultiply weights and Jacobian
      double W = w * J;


      // Do we want the velocity one?
      if ((which_one == 0) || (which_one == 2))
      {
        // Get the velocity shape functions
        shape_at_knot(ipt, psi);

        // Loop over the veclocity shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the velocity components
          for (unsigned i = 0; i < n_dim; i++)
          {
            local_eqn = nodal_local_eqn(l, u_nodal_index[i]);

            // If not a boundary condition
            if (local_eqn >= 0)
            {
              // Add the contribution
              veloc_mass_diag[local_eqn] += pow(psi[l], 2) * W;
            }
          }
        }
      }

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
        }
      }
    }
  }


  //=======================================================================
  /// Compute norm of the solution
  //=======================================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::compute_norm(double& norm)
  {
    // Initialise
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Solution
    double u = 0.0;

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

      for (unsigned i = 0; i < DIM; i++)
      {
        // Get FE function value
        u = this->interpolated_u_nst(s, i);

        // Add to  norm
        norm += u * u * W;
      }
    }
  }


  //======================================================================
  /// Compute the vector norm of FEM solution
  //======================================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::compute_norm(Vector<double>& norm)
  {
    // Resize the solution norm vector
    norm.resize(DIM + 1, 0.0);

    // Vector of local coordinates
    Vector<double> s(DIM, 0.0);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
      {
        // Compute the i-th local coordinate of this knot point
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
        norm[i] += std::pow(interpolated_u_nst(s, i), 2) * W;
      }

      // Update the pressure norm value
      norm[DIM] += std::pow(interpolated_p_nst(s), 2) * W;
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of compute_norm


  //======================================================================
  /// Validate against exact velocity solution at given time.
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points and compute L2 error
  /// and L2 norm of velocity solution over element.
  //=======================================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
    double& error,
    double& norm)
  {
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    outfile << "ZONE" << std::endl;

    // Exact solution Vector (u,v,[w],p)
    Vector<double> exact_soln(DIM + 1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
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
      for (unsigned i = 0; i < DIM; i++)
      {
        norm += exact_soln[i] * exact_soln[i] * W;
        error += (exact_soln[i] - interpolated_u_nst(s, i)) *
                 (exact_soln[i] - interpolated_u_nst(s, i)) * W;
      }

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }

      // Output x,y,[z],u_error,v_error,[w_error]
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << exact_soln[i] - interpolated_u_nst(s, i) << " ";
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
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();


    outfile << "ZONE" << std::endl;

    // Exact solution Vector (u,v,[w],p)
    Vector<double> exact_soln(DIM + 1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
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
      (*exact_soln_pt)(x, exact_soln);

      // Velocity error
      for (unsigned i = 0; i < DIM; i++)
      {
        norm += exact_soln[i] * exact_soln[i] * W;
        error += (exact_soln[i] - interpolated_u_nst(s, i)) *
                 (exact_soln[i] - interpolated_u_nst(s, i)) * W;
      }

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }

      // Output x,y,[z],u_error,v_error,[w_error]
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << exact_soln[i] - interpolated_u_nst(s, i) << " ";
      }
      outfile << std::endl;
    }
  }


  //======================================================================
  /// Validate against exact velocity solution at given time.
  /// Solution is provided via function pointer.
  /// Compute L2 error and L2 norm of velocity solution over element.
  //=======================================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::compute_error(
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
    double& error,
    double& norm)
  {
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Exact solution Vector (u,v,[w],p)
    Vector<double> exact_soln(DIM + 1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
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
      for (unsigned i = 0; i < DIM; i++)
      {
        norm += exact_soln[i] * exact_soln[i] * W;
        error += (exact_soln[i] - interpolated_u_nst(s, i)) *
                 (exact_soln[i] - interpolated_u_nst(s, i)) * W;
      }
    }
  } // End of compute_error


  //======================================================================
  /// Validate against exact velocity solution Solution is provided via a
  /// function pointer. Compute L2 error and L2 norm of velocity solution
  /// over element.
  //=======================================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::compute_error(
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    // Initialise error norm value
    error = 0.0;

    // Initialise exact solution norm value
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Exact solution Vector (u,v,[w],p)
    Vector<double> exact_soln(DIM + 1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
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
      (*exact_soln_pt)(x, exact_soln);

      // Velocity error
      for (unsigned i = 0; i < DIM; i++)
      {
        norm += exact_soln[i] * exact_soln[i] * W;
        error += (exact_soln[i] - interpolated_u_nst(s, i)) *
                 (exact_soln[i] - interpolated_u_nst(s, i)) * W;
      }
    } // for(unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of compute_error


  //======================================================================
  /// Output "exact" solution
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  /// Function prints as many components as are returned in solution Vector.
  //=======================================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

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
      for (unsigned i = 0; i < DIM; i++)
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
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::output_fct(
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
      for (unsigned i = 0; i < DIM; i++)
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
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::output_veloc(std::ostream& outfile,
                                                const unsigned& nplot,
                                                const unsigned& t)
  {
    // Find number of nodes
    unsigned n_node = nnode();

    // Local shape function
    Shape psi(n_node);

    // Vectors of local coordinates and coords and velocities
    Vector<double> s(DIM);
    Vector<double> interpolated_x(DIM);
    Vector<double> interpolated_u(DIM);

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
      for (unsigned i = 0; i < DIM; i++)
      {
        interpolated_x[i] = 0.0;
        interpolated_u[i] = 0.0;
        // Get the index at which velocity i is stored
        unsigned u_nodal_index = u_index_nst(i);
        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          interpolated_u[i] += nodal_value(t, l, u_nodal_index) * psi[l];
          interpolated_x[i] += nodal_position(t, l, i) * psi[l];
        }
      }

      // Coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x[i] << " ";
      }

      // Velocities
      for (unsigned i = 0; i < DIM; i++)
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
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::output(std::ostream& outfile,
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

      // Coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      // Velocities
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_u_nst(s, i) << " ";
      }

      // Pressure
      outfile << interpolated_p_nst(s) << " ";

      outfile << std::endl;
    }
    outfile << std::endl;

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //==============================================================
  /// C-style output function:
  /// x,y,[z],u,v,[w],p
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::output(FILE* file_pt, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    fprintf(file_pt, "%s", tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }

      // Velocities
      for (unsigned i = 0; i < DIM; i++)
      {
        fprintf(file_pt, "%g ", interpolated_u_nst(s, i));
      }

      // Pressure
      fprintf(file_pt, "%g \n", interpolated_p_nst(s));
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
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::full_output(std::ostream& outfile,
                                               const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Acceleration
    Vector<double> dudt(DIM);

    // Mesh elocity
    Vector<double> mesh_veloc(DIM, 0.0);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Set up memory for the shape functions
    Shape psif(n_node);
    DShape dpsifdx(n_node, DIM);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psif, dpsifdx);

      // Allocate storage
      Vector<double> mesh_veloc(DIM);
      Vector<double> dudt(DIM);
      Vector<double> dudt_ALE(DIM);
      DenseMatrix<double> interpolated_dudx(DIM, DIM);

      // Initialise everything to zero
      for (unsigned i = 0; i < DIM; i++)
      {
        mesh_veloc[i] = 0.0;
        dudt[i] = 0.0;
        dudt_ALE[i] = 0.0;
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_dudx(i, j) = 0.0;
        }
      }

      // Calculate velocities and derivatives


      // Loop over directions
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get the index at which velocity i is stored
        unsigned u_nodal_index = u_index_nst(i);
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          dudt[i] += du_dt_nst(l, u_nodal_index) * psif[l];
          mesh_veloc[i] += dnodal_position_dt(l, i) * psif[l];

          // Loop over derivative directions for velocity gradients
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_dudx(i, j) +=
              nodal_value(l, u_nodal_index) * dpsifdx(l, j);
          }
        }
      }


      // Get dudt in ALE form (incl mesh veloc)
      for (unsigned i = 0; i < DIM; i++)
      {
        dudt_ALE[i] = dudt[i];
        for (unsigned k = 0; k < DIM; k++)
        {
          dudt_ALE[i] -= mesh_veloc[k] * interpolated_dudx(i, k);
        }
      }


      // Coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      // Velocities
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_u_nst(s, i) << " ";
      }

      // Pressure
      outfile << interpolated_p_nst(s) << " ";

      // Accelerations
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << dudt_ALE[i] << " ";
      }

      //    // Mesh velocity
      //    for(unsigned i=0;i<DIM;i++)
      //     {
      //      outfile << mesh_veloc[i] << " ";
      //     }

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
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::output_vorticity(std::ostream& outfile,
                                                    const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Create vorticity vector of the required size
    Vector<double> vorticity;
    if (DIM == 2)
    {
      vorticity.resize(1);
    }
    else if (DIM == 3)
    {
      vorticity.resize(3);
    }
    else
    {
      std::string error_message =
        "Can't output vorticity in 1D - in fact, what do you mean\n";
      error_message += "by the 1D Navier-Stokes equations?\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      // Get vorticity
      get_vorticity(s, vorticity);

      if (DIM == 2)
      {
        outfile << vorticity[0];
      }
      else
      {
        outfile << vorticity[0] << " " << vorticity[1] << " " << vorticity[2]
                << " ";
      }

      outfile << std::endl;
    }
    outfile << std::endl;

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //==============================================================
  /// Return integral of dissipation over element
  //==============================================================
  template<unsigned DIM>
  double NavierStokesEquations<DIM>::dissipation() const
  {
    // Initialise
    double diss = 0.0;

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

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

      diss += local_diss * w * J;
    }

    return diss;
  }

  //==============================================================
  /// Compute traction (on the viscous scale) exerted onto
  /// the fluid at local coordinate s. N has to be outer unit normal
  /// to the fluid.
  //==============================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::get_traction(const Vector<double>& s,
                                                const Vector<double>& N,
                                                Vector<double>& traction)
  {
    // Get velocity gradients
    DenseMatrix<double> strainrate(DIM, DIM);
    strain_rate(s, strainrate);

    // Get pressure
    double press = interpolated_p_nst(s);

    // Loop over traction components
    for (unsigned i = 0; i < DIM; i++)
    {
      traction[i] = -press * N[i];
      for (unsigned j = 0; j < DIM; j++)
      {
        traction[i] += 2.0 * strainrate(i, j) * N[j];
      }
    }
  }

  //==============================================================
  /// Compute traction (on the viscous scale) exerted onto
  /// the fluid at local coordinate s, decomposed into pressure and
  /// normal and tangential viscous components.
  /// N has to be outer unit normal to the fluid.
  //==============================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::get_traction(const Vector<double>& s,
                                                const Vector<double>& N,
                                                Vector<double>& traction_p,
                                                Vector<double>& traction_visc_n,
                                                Vector<double>& traction_visc_t)
  {
    Vector<double> traction_visc(DIM);

    // Get velocity gradients
    DenseMatrix<double> strainrate(DIM, DIM);
    strain_rate(s, strainrate);

    // Get pressure
    double press = interpolated_p_nst(s);

    // Loop over traction components
    for (unsigned i = 0; i < DIM; i++)
    {
      // pressure
      traction_p[i] = -press * N[i];
      for (unsigned j = 0; j < DIM; j++)
      {
        // viscous stress
        traction_visc[i] += 2.0 * strainrate(i, j) * N[j];
      }
      // Decompose viscous stress into normal and tangential components
      traction_visc_n[i] = traction_visc[i] * N[i];
      traction_visc_t[i] = traction_visc[i] - traction_visc_n[i];
    }
  }

  //==============================================================
  /// Return dissipation at local coordinate s
  //==============================================================
  template<unsigned DIM>
  double NavierStokesEquations<DIM>::dissipation(const Vector<double>& s) const
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
  /// Get strain-rate tensor: 1/2 (du_i/dx_j + du_j/dx_i)
  //==============================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::strain_rate(
    const Vector<double>& s, DenseMatrix<double>& strainrate) const
  {
#ifdef PARANOID
    if ((strainrate.ncol() != DIM) || (strainrate.nrow() != DIM))
    {
      std::ostringstream error_message;
      error_message << "The strain rate has incorrect dimensions "
                    << strainrate.ncol() << " x " << strainrate.nrow()
                    << " Not " << DIM << std::endl;

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
      unsigned u_nodal_index = u_index_nst(i);
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

    // Loop over veclocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Loop over derivative directions
      for (unsigned j = 0; j < DIM; j++)
      {
        strainrate(i, j) = 0.5 * (dudx(i, j) + dudx(j, i));
      }
    }
  }


  //==============================================================
  /// Compute 2D vorticity vector at local coordinate s (return in
  /// one and only component of vorticity vector
  //==============================================================
  template<>
  void NavierStokesEquations<2>::get_vorticity(const Vector<double>& s,
                                               Vector<double>& vorticity) const
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
      unsigned u_nodal_index = u_index_nst(i);
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


  //==================================================================
  /// Compute vorticity vector at local coordinate s and return the
  /// one and only component of vorticity vector (only makes sense
  /// when solving the 2D N.St. equations)
  //==================================================================
  template<>
  void NavierStokesEquations<2>::get_vorticity(const Vector<double>& s,
                                               double& vorticity) const
  {
    // Create a vector to store the vorticity
    Vector<double> vort(1, 0.0);

    // Calculate the vorticity
    get_vorticity(s, vort);

    // Assign the vorticity
    vorticity = vort[0];
  }


  //==============================================================
  /// Compute 3D vorticity vector at local coordinate s
  //==============================================================
  template<>
  void NavierStokesEquations<3>::get_vorticity(const Vector<double>& s,
                                               Vector<double>& vorticity) const
  {
#ifdef PARANOID
    if (vorticity.size() != 3)
    {
      std::ostringstream error_message;
      error_message << "Computation of vorticity in 3D requires a 3D vector\n"
                    << "which contains the only non-zero component of the\n"
                    << "vorticity vector. You've passed a vector of size "
                    << vorticity.size() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Specify spatial dimension
    unsigned DIM = 3;

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
      // Get the index at which the i-th velocity component is stored
      unsigned u_nodal_index = u_index_nst(i);
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

    vorticity[0] = dudx(2, 1) - dudx(1, 2);
    vorticity[1] = dudx(0, 2) - dudx(2, 0);
    vorticity[2] = dudx(1, 0) - dudx(0, 1);
  }


  //==============================================================
  ///  \short Get integral of kinetic energy over element:
  /// Note that this is the "raw" kinetic energy in the sense
  /// that the density ratio has not been included. In problems
  /// with two or more fluids the user will have to remember
  /// to premultiply certain elements by the appropriate density
  /// ratio.
  //==============================================================
  template<unsigned DIM>
  double NavierStokesEquations<DIM>::kin_energy() const
  {
    // Initialise
    double kin_en = 0.0;

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
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
  }


  //==========================================================================
  ///  \short Get integral of time derivative of kinetic energy over element:
  //==========================================================================
  template<unsigned DIM>
  double NavierStokesEquations<DIM>::d_kin_energy_dt() const
  {
    // Initialise
    double d_kin_en_dt = 0.0;

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get the number of nodes
    const unsigned n_node = this->nnode();

    // Storage for the shape function
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

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
      Vector<double> interpolated_u(DIM, 0.0), interpolated_dudt(DIM, 0.0);

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
  }


  //==============================================================
  /// Return pressure integrated over the element
  //==============================================================
  template<unsigned DIM>
  double NavierStokesEquations<DIM>::pressure_integral() const
  {
    // Initialise
    double press_int = 0;

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
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
  void NavierStokesEquations<DIM>::
    fill_in_generic_pressure_advection_diffusion_contribution_nst(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
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
    DShape dpsidx(n_node, DIM);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);
    DShape dpsip(n_pres, DIM);
    DShape dtestp(n_pres, DIM);

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get Physical Variables from Element
    // Reynolds number must be multiplied by the density ratio
    double scaled_re = re() * density_ratio();

    // Integers to store the local equations and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++) s[i] = integral_pt()->knot(ipt, i);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the veloc shape functions
      // (Derivs not needed but they are free)
      double J = this->dshape_eulerian_at_knot(ipt, psif, dpsidx);

      // Call the pressure shape and test functions
      this->dpshape_and_dptest_eulerian_nst(s, psip, dpsip, testp, dtestp);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the pressure and velocity components
      // Allocate
      Vector<double> interpolated_u(DIM, 0.0);
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> interpolated_dpdx(DIM, 0.0);

      // Calculate pressure gradient
      for (unsigned l = 0; l < n_pres; l++)
      {
        for (unsigned i = 0; i < DIM; i++)
        {
          interpolated_dpdx[i] += p_nst(l) * dpsip(l, i);
        }
      }

      // Calculate velocities

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the nodal value
          double u_value = raw_nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif[l];
          interpolated_x[i] += raw_nodal_position(l, i) * psif[l];
        }
      }

      // Source function (for validaton only)
      double source = 0.0;
      if (Press_adv_diff_source_fct_pt != 0)
      {
        source = Press_adv_diff_source_fct_pt(interpolated_x);
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
                    jacobian(local_eqn, local_unknown) = 1.0;
                  }
                }
              }
            }
          } /*End of Jacobian calculation*/
        } // End of if not boundary condition
      } // End of loop over l
    }

    // Now add boundary contributions from Robin BCs
    unsigned nrobin = Pressure_advection_diffusion_robin_element_pt.size();
    for (unsigned e = 0; e < nrobin; e++)
    {
      Pressure_advection_diffusion_robin_element_pt[e]
        ->fill_in_generic_residual_contribution_fp_press_adv_diff_robin_bc(
          residuals, jacobian, flag);
    }
  }


  //==============================================================
  ///  Compute the residuals for the Navier--Stokes
  ///  equations; flag=1(or 0): do (or don't) compute the
  ///  Jacobian as well.
  //==============================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::fill_in_generic_residual_contribution_nst(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix,
    unsigned flag)
  {
    // Return immediately if there are no dofs
    if (ndof() == 0) return;

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Find out how many pressure dofs there are
    unsigned n_pres = npres_nst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_nodal_index[i] = u_index_nst(i);
    }

    // Set up memory for the shape and test functions
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, DIM), dtestfdx(n_node, DIM);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get Physical Variables from Element
    // Reynolds number must be multiplied by the density ratio
    double scaled_re = re() * density_ratio();
    double scaled_re_st = re_st() * density_ratio();
    double scaled_re_inv_fr = re_invfr() * density_ratio();
    double visc_ratio = viscosity_ratio();
    Vector<double> G = g();

    // Integers to store the local equations and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++) s[i] = integral_pt()->knot(ipt, i);
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Call the pressure shape and test functions
      pshape_nst(s, psip, testp);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the pressure and velocity components
      // Allocate
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(DIM, 0.0);
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> mesh_velocity(DIM, 0.0);
      Vector<double> dudt(DIM, 0.0);
      DenseMatrix<double> interpolated_dudx(DIM, DIM, 0.0);

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
        interpolated_p += p_nst(l) * psip[l];

      // Calculate velocities and derivatives:

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the nodal value
          double u_value = raw_nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif[l];
          interpolated_x[i] += raw_nodal_position(l, i) * psif[l];
          dudt[i] += du_dt_nst(l, i) * psif[l];

          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      if (!ALE_is_disabled)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < DIM; i++)
          {
            mesh_velocity[i] += this->raw_dnodal_position_dt(l, i) * psif[l];
          }
        }
      }

      // Get the user-defined body force terms
      Vector<double> body_force(DIM);
      get_body_force_nst(time, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      double source = get_source_nst(time, ipt, interpolated_x);


      // MOMENTUM EQUATIONS
      //------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the velocity components
        for (unsigned i = 0; i < DIM; i++)
        {
          /*IF it's not a boundary condition*/
          local_eqn = nodal_local_eqn(l, u_nodal_index[i]);
          if (local_eqn >= 0)
          {
            // Add the user-defined body force terms
            residuals[local_eqn] += body_force[i] * testf[l] * W;

            // Add the gravitational body force term
            residuals[local_eqn] += scaled_re_inv_fr * testf[l] * G[i] * W;

            // Add the pressure gradient term
            residuals[local_eqn] += interpolated_p * dtestfdx(l, i) * W;

            // Add in the stress tensor terms
            // The viscosity ratio needs to go in here to ensure
            // continuity of normal stress is satisfied even in flows
            // with zero pressure gradient!
            for (unsigned k = 0; k < DIM; k++)
            {
              residuals[local_eqn] -=
                visc_ratio *
                (interpolated_dudx(i, k) + Gamma[i] * interpolated_dudx(k, i)) *
                dtestfdx(l, k) * W;
            }

            // Add in the inertial terms
            // du/dt term
            residuals[local_eqn] -= scaled_re_st * dudt[i] * testf[l] * W;


            // Convective terms, including mesh velocity
            for (unsigned k = 0; k < DIM; k++)
            {
              double tmp = scaled_re * interpolated_u[k];
              if (!ALE_is_disabled) tmp -= scaled_re_st * mesh_velocity[k];
              residuals[local_eqn] -=
                tmp * interpolated_dudx(i, k) * testf[l] * W;
            }

            // CALCULATE THE JACOBIAN
            if (flag)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Loop over the velocity components again
                for (unsigned i2 = 0; i2 < DIM; i2++)
                {
                  // If at a non-zero degree of freedom add in the entry
                  local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
                  if (local_unknown >= 0)
                  {
                    // Add contribution to Elemental Matrix
                    jacobian(local_eqn, local_unknown) -=
                      visc_ratio * Gamma[i] * dpsifdx(l2, i) * dtestfdx(l, i2) *
                      W;

                    // Extra component if i2 = i
                    if (i2 == i)
                    {
                      /*Loop over velocity components*/
                      for (unsigned k = 0; k < DIM; k++)
                      {
                        jacobian(local_eqn, local_unknown) -=
                          visc_ratio * dpsifdx(l2, k) * dtestfdx(l, k) * W;
                      }
                    }

                    // Now add in the inertial terms
                    jacobian(local_eqn, local_unknown) -=
                      scaled_re * psif[l2] * interpolated_dudx(i, i2) *
                      testf[l] * W;

                    // Extra component if i2=i
                    if (i2 == i)
                    {
                      // Add the mass matrix term (only diagonal entries)
                      // Note that this is positive because the mass matrix
                      // is taken to the other side of the equation when
                      // formulating the generalised eigenproblem.
                      if (flag == 2)
                      {
                        mass_matrix(local_eqn, local_unknown) +=
                          scaled_re_st * psif[l2] * testf[l] * W;
                      }

                      // du/dt term
                      jacobian(local_eqn, local_unknown) -=
                        scaled_re_st *
                        node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                        psif[l2] * testf[l] * W;

                      // Loop over the velocity components
                      for (unsigned k = 0; k < DIM; k++)
                      {
                        double tmp = scaled_re * interpolated_u[k];
                        if (!ALE_is_disabled)
                          tmp -= scaled_re_st * mesh_velocity[k];
                        jacobian(local_eqn, local_unknown) -=
                          tmp * dpsifdx(l2, k) * testf[l] * W;
                      }
                    }
                  }
                }
              }

              /*Now loop over pressure shape functions*/
              /*This is the contribution from pressure gradient*/
              for (unsigned l2 = 0; l2 < n_pres; l2++)
              {
                /*If we are at a non-zero degree of freedom in the entry*/
                local_unknown = p_local_eqn(l2);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    psip[l2] * dtestfdx(l, i) * W;
                }
              }
            } /*End of Jacobian calculation*/

          } // End of if not boundary condition statement

        } // End of loop over dimension
      } // End of loop over shape functions


      // CONTINUITY EQUATION
      //-------------------

      // Loop over the shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        local_eqn = p_local_eqn(l);
        // If not a boundary conditions
        if (local_eqn >= 0)
        {
          // Source term
          // residuals[local_eqn] -=source*testp[l]*W;
          double aux = -source;

          // Loop over velocity components
          for (unsigned k = 0; k < DIM; k++)
          {
            // residuals[local_eqn] += interpolated_dudx(k,k)*testp[l]*W;
            aux += interpolated_dudx(k, k);
          }

          residuals[local_eqn] += aux * testp[l] * W;

          /*CALCULATE THE JACOBIAN*/
          if (flag)
          {
            /*Loop over the velocity shape functions*/
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              /*Loop over velocity components*/
              for (unsigned i2 = 0; i2 < DIM; i2++)
              {
                /*If we're at a non-zero degree of freedom add it in*/
                local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    dpsifdx(l2, i2) * testp[l] * W;
                }
              } /*End of loop over i2*/
            } /*End of loop over l2*/
          } /*End of Jacobian calculation*/

        } // End of if not boundary condition
      } // End of loop over l
    }
  }

  //==============================================================
  ///  Compute the derivatives of the residuals for the Navier--Stokes
  ///  equations with respect to a parameter;
  /// flag=2 or 1 or 0: do (or don't) compute the
  ///  Jacobian and mass matrix as well
  //==============================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::fill_in_generic_dresidual_contribution_nst(
    double* const& parameter_pt,
    Vector<double>& dres_dparam,
    DenseMatrix<double>& djac_dparam,
    DenseMatrix<double>& dmass_matrix_dparam,
    unsigned flag)
  {
    throw OomphLibError("Not yet implemented\n",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

  //==================================================================
  /// \short Compute the hessian tensor vector products required to
  /// perform continuation of bifurcations analytically
  //================================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::
    fill_in_contribution_to_hessian_vector_products(
      Vector<double> const& Y,
      DenseMatrix<double> const& C,
      DenseMatrix<double>& product)
  {
    throw OomphLibError("Not yet implemented\n",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }


  //======================================================================
  /// Compute derivatives of elemental residual vector with respect
  /// to nodal coordinates.
  /// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
  /// Overloads the FD-based version in the FE base class.
  //======================================================================
  template<unsigned DIM>
  void NavierStokesEquations<DIM>::get_dresidual_dnodal_coordinates(
    RankThreeTensor<double>& dresidual_dnodal_coordinates)
  {
    // Return immediately if there are no dofs
    if (ndof() == 0)
    {
      return;
    }

    // Determine number of nodes in element
    const unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Determine number of pressure dofs in element
    const unsigned n_pres = npres_nst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_nodal_index[i] = u_index_nst(i);
    }

    // Set up memory for the shape and test functions
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, DIM), dtestfdx(n_node, DIM);

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
    Vector<double> s(DIM);

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

          // Do auxiliary node update (to apply no slip)
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
      for (unsigned i = 0; i < DIM; i++)
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

      // Calculate local values of the pressure and velocity components
      // Allocate
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(DIM, 0.0);
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> mesh_velocity(DIM, 0.0);
      Vector<double> dudt(DIM, 0.0);
      DenseMatrix<double> interpolated_dudx(DIM, DIM, 0.0);

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
      {
        interpolated_p += p_nst(l) * psip[l];
      }

      // Calculate velocities and derivatives:

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the nodal value
          double u_value = raw_nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif[l];
          interpolated_x[i] += raw_nodal_position(l, i) * psif[l];
          dudt[i] += du_dt_nst(l, i) * psif[l];

          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      if (!ALE_is_disabled)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < DIM; i++)
          {
            mesh_velocity[i] += this->raw_dnodal_position_dt(l, i) * psif[l];
          }
        }
      }

      // Calculate derivative of du_i/dx_k w.r.t. nodal positions X_{pq}
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
              d_dudx_dX(p, q, i, k) = aux;
            }
          }
        }
      }

      // Get weight of actual nodal position/value in computation of mesh
      // velocity from positional/value time stepper
      const double pos_time_weight =
        node_pt(0)->position_time_stepper_pt()->weight(1, 0);
      const double val_time_weight =
        node_pt(0)->time_stepper_pt()->weight(1, 0);

      // Get the user-defined body force terms
      Vector<double> body_force(DIM);
      get_body_force_nst(time, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      const double source = get_source_nst(time, ipt, interpolated_x);

      // Note: Can use raw values and avoid bypassing hanging information
      // because this is the non-refineable version!

      // Get gradient of body force function
      DenseMatrix<double> d_body_force_dx(DIM, DIM, 0.0);
      get_body_force_gradient_nst(
        time, ipt, s, interpolated_x, d_body_force_dx);

      // Get gradient of source function
      Vector<double> source_gradient(DIM, 0.0);
      get_source_gradient_nst(time, ipt, interpolated_x, source_gradient);


      // Assemble shape derivatives
      // --------------------------

      // MOMENTUM EQUATIONS
      // ------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over coordinate directions
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the local equation
          local_eqn = nodal_local_eqn(l, u_nodal_index[i]);
          ;

          // IF it's not a boundary condition
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

                // Add the user-defined body force terms
                double sum = body_force[i] * testf[l];

                // Add the gravitational body force term
                sum += scaled_re_inv_fr * testf[l] * G[i];

                // Add the pressure gradient term
                sum += interpolated_p * dtestfdx(l, i);

                // Add in the stress tensor terms
                // The viscosity ratio needs to go in here to ensure
                // continuity of normal stress is satisfied even in flows
                // with zero pressure gradient!
                for (unsigned k = 0; k < DIM; k++)
                {
                  sum -= visc_ratio *
                         (interpolated_dudx(i, k) +
                          Gamma[i] * interpolated_dudx(k, i)) *
                         dtestfdx(l, k);
                }

                // Add in the inertial terms

                // du/dt term
                sum -= scaled_re_st * dudt[i] * testf[l];

                // Convective terms, including mesh velocity
                for (unsigned k = 0; k < DIM; k++)
                {
                  double tmp = scaled_re * interpolated_u[k];
                  if (!ALE_is_disabled)
                  {
                    tmp -= scaled_re_st * mesh_velocity[k];
                  }
                  sum -= tmp * interpolated_dudx(i, k) * testf[l];
                }

                // Multiply through by deriv of Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) +=
                  sum * dJ_dX(p, q) * w;

                // Derivative of residual x Jacobian
                //----------------------------------

                // Body force
                sum = d_body_force_dx(i, p) * psif(q) * testf(l);

                // Pressure gradient term
                sum += interpolated_p * d_dtestfdx_dX(p, q, l, i);

                // Viscous term
                for (unsigned k = 0; k < DIM; k++)
                {
                  sum -= visc_ratio * ((interpolated_dudx(i, k) +
                                        Gamma[i] * interpolated_dudx(k, i)) *
                                         d_dtestfdx_dX(p, q, l, k) +
                                       (d_dudx_dX(p, q, i, k) +
                                        Gamma[i] * d_dudx_dX(p, q, k, i)) *
                                         dtestfdx(l, k));
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
  }


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  /// 2D Crouzeix-Raviart elements
  // Set the data for the number of Variables at each node
  template<>
  const unsigned QCrouzeixRaviartElement<2>::Initial_Nvalue[9] = {
    2, 2, 2, 2, 2, 2, 2, 2, 2};

  /// 3D Crouzeix-Raviart elements
  // Set the data for the number of Variables at each node
  template<>
  const unsigned QCrouzeixRaviartElement<3>::Initial_Nvalue[27] = {
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};


  //========================================================================
  /// Number of values (pinned or dofs) required at node n.
  //========================================================================
  template<unsigned DIM>
  unsigned QCrouzeixRaviartElement<DIM>::required_nvalue(
    const unsigned& n) const
  {
    return Initial_Nvalue[n];
  }


  //=========================================================================
  ///  Add to the set \c paired_load_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for all values (pressures, velocities) that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  template<unsigned DIM>
  void QCrouzeixRaviartElement<DIM>::identify_load_data(
    std::set<std::pair<Data*, unsigned>>& paired_load_data)
  {
    // Find the index at which the velocity is stored
    unsigned u_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_index[i] = this->u_index_nst(i);
    }

    // Loop over the nodes
    unsigned n_node = this->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      // Loop over the velocity components and add pointer to their data
      // and indices to the vectors
      for (unsigned i = 0; i < DIM; i++)
      {
        paired_load_data.insert(std::make_pair(this->node_pt(n), u_index[i]));
      }
    }

    // Identify the pressure data
    identify_pressure_data(paired_load_data);
  }


  //=========================================================================
  ///  Add to the set \c paired_pressue_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for all pressures values that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  template<unsigned DIM>
  void QCrouzeixRaviartElement<DIM>::identify_pressure_data(
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
  /// of the DOF that this unknown is associated with.
  /// (Function can obviously only be called if the equation numbering
  /// scheme has been set up.)
  //=============================================================================
  template<unsigned DIM>
  void QCrouzeixRaviartElement<DIM>::get_dof_numbers_for_unknowns(
    std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // number of nodes
    unsigned n_node = this->nnode();

    // number of pressure values
    unsigned n_press = this->npres_nst();

    // temporary pair (used to store dof lookup prior to being added to list)
    std::pair<unsigned, unsigned> dof_lookup;

    // pressure dof number
    unsigned pressure_dof_number = DIM;

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


  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  // 2D Taylor--Hood
  // Set the data for the number of Variables at each node
  template<>
  const unsigned QTaylorHoodElement<2>::Initial_Nvalue[9] = {
    3, 2, 3, 2, 2, 2, 3, 2, 3};

  // Set the data for the pressure conversion array
  template<>
  const unsigned QTaylorHoodElement<2>::Pconv[4] = {0, 2, 6, 8};

  // 3D Taylor--Hood
  // Set the data for the number of Variables at each node
  template<>
  const unsigned QTaylorHoodElement<3>::Initial_Nvalue[27] = {
    4, 3, 4, 3, 3, 3, 4, 3, 4, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 4, 3, 4, 3, 3, 3, 4, 3, 4};

  // Set the data for the pressure conversion array
  template<>
  const unsigned QTaylorHoodElement<3>::Pconv[8] = {0, 2, 6, 8, 18, 20, 24, 26};

  //=========================================================================
  ///  Add to the set \c paired_load_data pairs containing
  /// - the pointer to a Data object
  /// and
  /// - the index of the value in that Data object
  /// .
  /// for all values (pressures, velocities) that affect the
  /// load computed in the \c get_load(...) function.
  //=========================================================================
  template<unsigned DIM>
  void QTaylorHoodElement<DIM>::identify_load_data(
    std::set<std::pair<Data*, unsigned>>& paired_load_data)
  {
    // Find the index at which the velocity is stored
    unsigned u_index[DIM];
    for (unsigned i = 0; i < DIM; i++)
    {
      u_index[i] = this->u_index_nst(i);
    }

    // Loop over the nodes
    unsigned n_node = this->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      // Loop over the velocity components and add pointer to their data
      // and indices to the vectors
      for (unsigned i = 0; i < DIM; i++)
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
  /// for pressure values that affect the
  /// load computed in the \c get_load(...) function.,
  //=========================================================================
  template<unsigned DIM>
  void QTaylorHoodElement<DIM>::identify_pressure_data(
    std::set<std::pair<Data*, unsigned>>& paired_pressure_data)
  {
    // Find the index at which the pressure is stored
    unsigned p_index = static_cast<unsigned>(this->p_nodal_index_nst());

    // Loop over the pressure data
    unsigned n_pres = npres_nst();
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
  template<unsigned DIM>
  void QTaylorHoodElement<DIM>::get_dof_numbers_for_unknowns(
    std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // number of nodes
    unsigned n_node = this->nnode();

    // local eqn no for pressure unknown
    // unsigned p_index = this->p_nodal_index_nst();

    // temporary pair (used to store dof lookup prior to being added to list)
    std::pair<unsigned, unsigned> dof_lookup;

    // loop over the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // find the number of values at this node
      unsigned nv = this->required_nvalue(n);

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


  //====================================================================
  //// Force build of templates
  //====================================================================
  template class NavierStokesEquations<2>;
  template class QCrouzeixRaviartElement<2>;
  template class QTaylorHoodElement<2>;

  template class NavierStokesEquations<3>;
  template class QCrouzeixRaviartElement<3>;
  template class QTaylorHoodElement<3>;

} // namespace oomph
