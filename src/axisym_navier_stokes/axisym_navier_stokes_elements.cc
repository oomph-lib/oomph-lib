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
#include "axisym_navier_stokes_elements.h"


namespace oomph
{
  /// Navier--Stokes equations static data
  Vector<double> AxisymmetricNavierStokesEquations::Gamma(2, 1.0);

  //=================================================================
  /// "Magic" negative number that indicates that the pressure is
  /// not stored at a node. This cannot be -1 because that represents
  /// the positional hanging scheme in the hanging_pt object of nodes
  //=================================================================
  int AxisymmetricNavierStokesEquations::Pressure_not_stored_at_node = -100;

  /// Navier--Stokes equations static data
  double AxisymmetricNavierStokesEquations::Default_Physical_Constant_Value =
    0.0;

  // Navier--Stokes equations static data
  double AxisymmetricNavierStokesEquations::Default_Physical_Ratio_Value = 1.0;

  /// Navier-Stokes equations default gravity vector
  Vector<double> AxisymmetricNavierStokesEquations::Default_Gravity_vector(3,
                                                                           0.0);


  //================================================================
  /// Compute the diagonal of the velocity/pressure mass matrices.
  /// If which one=0, both are computed, otherwise only the pressure
  /// (which_one=1) or the velocity mass matrix (which_one=2 -- the
  /// LSC version of the preconditioner only needs that one)
  /// NOTE: pressure versions isn't implemented yet because this
  ///       element has never been tried with Fp preconditoner.
  //================================================================
  void AxisymmetricNavierStokesEquations::
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
      u_nodal_index[i] = this->u_index_axi_nst(i);
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
      for (unsigned l = 0; l < n_node; l++)
      {
        r += this->nodal_position(l, 0) * test(l);
      }

      // Multiply by the geometric factor
      W *= r;

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
  void AxisymmetricNavierStokesEquations::compute_error(
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

    // Set the value of Nintpt
    unsigned Nintpt = integral_pt()->nweight();

    outfile << "ZONE" << std::endl;

    // Exact solution Vector (u,v,w,p)
    Vector<double> exact_soln(4);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < Nintpt; ipt++)
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

      // Get x position as Vector
      interpolated_x(s, x);

      // Premultiply the weights and the Jacobian and r
      double W = w * J * x[0];

      // Get exact solution at this point
      (*exact_soln_pt)(time, x, exact_soln);

      // Velocity error
      for (unsigned i = 0; i < 3; i++)
      {
        norm += exact_soln[i] * exact_soln[i] * W;
        error += (exact_soln[i] - interpolated_u_axi_nst(s, i)) *
                 (exact_soln[i] - interpolated_u_axi_nst(s, i)) * W;
      }

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output x,y,[z],u_error,v_error,[w_error]
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << exact_soln[i] - interpolated_u_axi_nst(s, i) << " ";
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
  void AxisymmetricNavierStokesEquations::compute_error(
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

    // Set the value of Nintpt
    unsigned Nintpt = integral_pt()->nweight();

    outfile << "ZONE" << std::endl;

    // Exact solution Vector (u,v,w,p)
    Vector<double> exact_soln(4);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < Nintpt; ipt++)
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

      // Get x position as Vector
      interpolated_x(s, x);

      // Premultiply the weights and the Jacobian and r
      double W = w * J * x[0];

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Velocity error
      for (unsigned i = 0; i < 3; i++)
      {
        norm += exact_soln[i] * exact_soln[i] * W;
        error += (exact_soln[i] - interpolated_u_axi_nst(s, i)) *
                 (exact_soln[i] - interpolated_u_axi_nst(s, i)) * W;
      }

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output x,y,u_error,v_error,w_error
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << exact_soln[i] - interpolated_u_axi_nst(s, i) << " ";
      }

      outfile << std::endl;
    }
  }

  //======================================================================
  /// Output "exact" solution
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  /// Function prints as many components as are returned in solution Vector.
  //=======================================================================
  void AxisymmetricNavierStokesEquations::output_fct(
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
      for (unsigned i = 0; i < 2; i++)
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
  void AxisymmetricNavierStokesEquations::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
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
      (*exact_soln_pt)(time, x, exact_soln);

      // Output x,y,...
      for (unsigned i = 0; i < 2; i++)
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
  void AxisymmetricNavierStokesEquations::output_veloc(std::ostream& outfile,
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

      // Loop over coordinate directions
      for (unsigned i = 0; i < 2; i++)
      {
        interpolated_x[i] = 0.0;
        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          interpolated_x[i] += nodal_position(t, l, i) * psi[l];
        }
      }

      // Loop over the velocity components
      for (unsigned i = 0; i < 3; i++)
      {
        // Get the index at which the velocity is stored
        unsigned u_nodal_index = u_index_axi_nst(i);
        interpolated_u[i] = 0.0;
        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          interpolated_u[i] += nodal_value(t, l, u_nodal_index) * psi[l];
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
  /// r,z,u,v,w,p
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  void AxisymmetricNavierStokesEquations::output(std::ostream& outfile,
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

      // Coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      // Velocities
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << interpolated_u_axi_nst(s, i) << " ";
      }

      // Pressure
      outfile << interpolated_p_axi_nst(s) << " ";

      outfile << std::endl;
    }
    outfile << std::endl;

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //==============================================================
  /// Output function:
  /// r,z,u,v,w,p
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  void AxisymmetricNavierStokesEquations::output(FILE* file_pt,
                                                 const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Tecplot header info
    fprintf(file_pt, "%s ", tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        // outfile << interpolated_x(s,i) << " ";
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }

      // Velocities
      for (unsigned i = 0; i < 3; i++)
      {
        // outfile << interpolated_u(s,i) << " ";
        fprintf(file_pt, "%g ", interpolated_u_axi_nst(s, i));
      }

      // Pressure
      // outfile << interpolated_p(s)  << " ";
      fprintf(file_pt, "%g ", interpolated_p_axi_nst(s));

      // outfile << std::endl;
      fprintf(file_pt, "\n");
    }
    // outfile << std::endl;
    fprintf(file_pt, "\n");

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, nplot);
  }


  //==============================================================
  /// Return integral of dissipation over element
  //==============================================================
  double AxisymmetricNavierStokesEquations::dissipation() const
  {
    throw OomphLibError(
      "Check the dissipation calculation for axisymmetric NSt",
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);

    // Initialise
    double diss = 0.0;

    // Set the value of Nintpt
    unsigned Nintpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < Nintpt; ipt++)
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
  /// \short Compute traction (on the viscous scale) at local
  /// coordinate s for outer unit normal N
  //==============================================================
  void AxisymmetricNavierStokesEquations::traction(
    const Vector<double>& s,
    const Vector<double>& N,
    Vector<double>& traction) const
  {
    // throw OomphLibError(
    // "Check the traction calculation for axisymmetric NSt",
    // OOMPH_CURRENT_FUNCTION,
    // OOMPH_EXCEPTION_LOCATION);

    // Pad out normal vector if required
    Vector<double> n_local(3, 0.0);
    n_local[0] = N[0];
    n_local[1] = N[1];

#ifdef PARANOID
    if ((N.size() == 3) && (N[2] != 0.0))
    {
      throw OomphLibError(
        "Unit normal passed into this fct should either be 2D (r,z) or have a "
        "zero component in the theta-direction",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Get velocity gradients
    DenseMatrix<double> strainrate(3, 3);
    strain_rate(s, strainrate);

    // Get pressure
    double press = interpolated_p_axi_nst(s);

    // Loop over traction components
    for (unsigned i = 0; i < 3; i++)
    {
      traction[i] = -press * n_local[i];
      for (unsigned j = 0; j < 3; j++)
      {
        traction[i] += 2.0 * strainrate(i, j) * n_local[j];
      }
    }
  }

  //==============================================================
  /// Return dissipation at local coordinate s
  //==============================================================
  double AxisymmetricNavierStokesEquations::dissipation(
    const Vector<double>& s) const
  {
    throw OomphLibError(
      "Check the dissipation calculation for axisymmetric NSt",
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);

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
  /// Get strain-rate tensor: \f$ e_{ij} \f$  where
  /// \f$ i,j = r,z,\theta \f$ (in that order)
  //==============================================================
  void AxisymmetricNavierStokesEquations::strain_rate(
    const Vector<double>& s, DenseMatrix<double>& strainrate) const
  {
#ifdef PARANOID
    if ((strainrate.ncol() != 3) || (strainrate.nrow() != 3))
    {
      std::ostringstream error_message;
      error_message << "The strain rate has incorrect dimensions "
                    << strainrate.ncol() << " x " << strainrate.nrow()
                    << " Not 3" << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node, 2);

    // Call the derivatives of the shape functions
    dshape_eulerian(s, psi, dpsidx);

    // Radius
    double interpolated_r = 0.0;

    // Velocity components and their derivatives
    double ur = 0.0;
    double durdr = 0.0;
    double durdz = 0.0;
    double uz = 0.0;
    double duzdr = 0.0;
    double duzdz = 0.0;
    double uphi = 0.0;
    double duphidr = 0.0;
    double duphidz = 0.0;

    // Get the local storage for the indices
    unsigned u_nodal_index[3];
    for (unsigned i = 0; i < 3; ++i)
    {
      u_nodal_index[i] = u_index_axi_nst(i);
    }

    // Loop over nodes to assemble velocities and their derivatives
    // w.r.t. to r and z (x_0 and x_1)
    for (unsigned l = 0; l < n_node; l++)
    {
      interpolated_r += nodal_position(l, 0) * psi[l];

      ur += nodal_value(l, u_nodal_index[0]) * psi[l];
      uz += nodal_value(l, u_nodal_index[1]) * psi[l];
      uphi += nodal_value(l, u_nodal_index[2]) * psi[l];

      durdr += nodal_value(l, u_nodal_index[0]) * dpsidx(l, 0);
      durdz += nodal_value(l, u_nodal_index[0]) * dpsidx(l, 1);

      duzdr += nodal_value(l, u_nodal_index[1]) * dpsidx(l, 0);
      duzdz += nodal_value(l, u_nodal_index[1]) * dpsidx(l, 1);

      duphidr += nodal_value(l, u_nodal_index[2]) * dpsidx(l, 0);
      duphidz += nodal_value(l, u_nodal_index[2]) * dpsidx(l, 1);
    }


    // Assign strain rates without negative powers of the radius
    // and zero those with:
    strainrate(0, 0) = durdr;
    strainrate(0, 1) = 0.5 * (durdz + duzdr);
    strainrate(1, 0) = strainrate(0, 1);
    strainrate(0, 2) = 0.0;
    strainrate(2, 0) = strainrate(0, 2);
    strainrate(1, 1) = duzdz;
    strainrate(1, 2) = 0.5 * duphidz;
    strainrate(2, 1) = strainrate(1, 2);
    strainrate(2, 2) = 0.0;


    // Overwrite the strain rates with negative powers of the radius
    // unless we're at the origin
    if (std::fabs(interpolated_r) > 1.0e-16)
    {
      double inverse_radius = 1.0 / interpolated_r;
      strainrate(0, 2) = 0.5 * (duphidr - inverse_radius * uphi);
      strainrate(2, 0) = strainrate(0, 2);
      strainrate(2, 2) = inverse_radius * ur;
    }
  }


  //==============================================================
  ///  \short Get integral of kinetic energy over element:
  //==============================================================
  double AxisymmetricNavierStokesEquations::kin_energy() const
  {
    throw OomphLibError(
      "Check the kinetic energy calculation for axisymmetric NSt",
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);

    // Initialise
    double kin_en = 0.0;

    // Set the value of Nintpt
    unsigned Nintpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < Nintpt; ipt++)
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
        veloc_squared +=
          interpolated_u_axi_nst(s, i) * interpolated_u_axi_nst(s, i);
      }

      kin_en += 0.5 * veloc_squared * w * J * interpolated_x(s, 0);
    }

    return kin_en;
  }

  //==============================================================
  /// Return pressure integrated over the element
  //==============================================================
  double AxisymmetricNavierStokesEquations::pressure_integral() const
  {
    // Initialise
    double press_int = 0;

    // Set the value of Nintpt
    unsigned Nintpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < Nintpt; ipt++)
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
      double W = w * J * interpolated_x(s, 0);

      // Get pressure
      double press = interpolated_p_axi_nst(s);

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
  void AxisymmetricNavierStokesEquations::
    fill_in_generic_residual_contribution_axi_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Find out how many pressure dofs there are
    unsigned n_pres = npres_axi_nst();

    // Get the nodal indices at which the velocity is stored
    unsigned u_nodal_index[3];
    for (unsigned i = 0; i < 3; ++i)
    {
      u_nodal_index[i] = u_index_axi_nst(i);
    }

    // Set up memory for the shape and test functions
    // Note that there are only two dimensions, r and z in this problem
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, 2), dtestfdx(n_node, 2);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Number of integration points
    unsigned Nintpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates (two dimensions)
    Vector<double> s(2);

    // Get Physical Variables from Element
    // Reynolds number must be multiplied by the density ratio
    double scaled_re = re() * density_ratio();
    double scaled_re_st = re_st() * density_ratio();
    double scaled_re_inv_fr = re_invfr() * density_ratio();
    double scaled_re_inv_ro = re_invro() * density_ratio();
    // double visc_ratio = viscosity_ratio();
    Vector<double> G = g();

    // Integers used to store the local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < Nintpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++) s[i] = integral_pt()->knot(ipt, i);
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_axi_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Call the pressure shape and test functions
      pshape_axi_nst(s, psip, testp);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Allocate storage for the position and the derivative of the
      // mesh positions wrt time
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> mesh_velocity(2, 0.0);
      // Allocate storage for the pressure, velocity components and their
      // time and spatial derivatives
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(3, 0.0);
      Vector<double> dudt(3, 0.0);
      DenseMatrix<double> interpolated_dudx(3, 2, 0.0);

      // Calculate pressure at integration point
      for (unsigned l = 0; l < n_pres; l++)
      {
        interpolated_p += p_axi_nst(l) * psip[l];
      }

      // Calculate velocities and derivatives at integration point

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the shape function
        const double psif_ = psif(l);
        // Loop over the two coordinate directions
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += this->raw_nodal_position(l, i) * psif_;
        }
        // mesh_velocity[i]  += dnodal_position_dt(l,i)*psif[l];

        // Loop over the three velocity directions
        for (unsigned i = 0; i < 3; i++)
        {
          // Get the u_value
          const double u_value = this->raw_nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif_;
          dudt[i] += du_dt_axi_nst(l, i) * psif_;
          // Loop over derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      // Get the mesh velocity if ALE is enabled
      if (!ALE_is_disabled)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the two coordinate directions
          for (unsigned i = 0; i < 2; i++)
          {
            mesh_velocity[i] += this->raw_dnodal_position_dt(l, i) * psif(l);
          }
        }
      }


      // Get the user-defined body force terms
      Vector<double> body_force(3);
      get_body_force_axi_nst(time, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      double source = get_source_fct(time, ipt, interpolated_x);

      // Get the user-defined viscosity function
      double visc_ratio;
      get_viscosity_ratio_axisym_nst(ipt, s, interpolated_x, visc_ratio);

      // r is the first position component
      double r = interpolated_x[0];

      // MOMENTUM EQUATIONS
      //------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // FIRST (RADIAL) MOMENTUM EQUATION
        local_eqn = nodal_local_eqn(l, u_nodal_index[0]);
        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Add the user-defined body force terms
          residuals[local_eqn] += r * body_force[0] * testf[l] * W;

          // Add the gravitational body force term
          residuals[local_eqn] += r * scaled_re_inv_fr * testf[l] * G[0] * W;

          // Add the pressure gradient term
          residuals[local_eqn] +=
            interpolated_p * (testf[l] + r * dtestfdx(l, 0)) * W;

          // Add in the stress tensor terms
          // The viscosity ratio needs to go in here to ensure
          // continuity of normal stress is satisfied even in flows
          // with zero pressure gradient!
          residuals[local_eqn] -= visc_ratio * r * (1.0 + Gamma[0]) *
                                  interpolated_dudx(0, 0) * dtestfdx(l, 0) * W;

          residuals[local_eqn] -=
            visc_ratio * r *
            (interpolated_dudx(0, 1) + Gamma[0] * interpolated_dudx(1, 0)) *
            dtestfdx(l, 1) * W;

          residuals[local_eqn] -= visc_ratio * (1.0 + Gamma[0]) *
                                  interpolated_u[0] * testf[l] * W / r;

          // Add in the inertial terms
          // du/dt term
          residuals[local_eqn] -= scaled_re_st * r * dudt[0] * testf[l] * W;

          // Convective terms
          residuals[local_eqn] -=
            scaled_re *
            (r * interpolated_u[0] * interpolated_dudx(0, 0) -
             interpolated_u[2] * interpolated_u[2] +
             r * interpolated_u[1] * interpolated_dudx(0, 1)) *
            testf[l] * W;

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            for (unsigned k = 0; k < 2; k++)
            {
              residuals[local_eqn] += scaled_re_st * r * mesh_velocity[k] *
                                      interpolated_dudx(0, k) * testf[l] * W;
            }
          }

          // Add the Coriolis term
          residuals[local_eqn] +=
            2.0 * r * scaled_re_inv_ro * interpolated_u[2] * testf[l] * W;

          // CALCULATE THE JACOBIAN
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              // Radial velocity component
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  // Add the mass matrix
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] * testf[l] * W;
                }

                // Add contribution to the Jacobian matrix
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * (1.0 + Gamma[0]) * dpsifdx(l2, 0) *
                  dtestfdx(l, 0) * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * dpsifdx(l2, 1) * dtestfdx(l, 1) * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * (1.0 + Gamma[0]) * psif[l2] * testf[l] * W / r;

                // Add in the inertial terms
                // du/dt term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re_st * r *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif[l2] *
                  testf[l] * W;

                // Convective terms
                jacobian(local_eqn, local_unknown) -=
                  scaled_re *
                  (r * psif[l2] * interpolated_dudx(0, 0) +
                   r * interpolated_u[0] * dpsifdx(l2, 0) +
                   r * interpolated_u[1] * dpsifdx(l2, 1)) *
                  testf[l] * W;

                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned k = 0; k < 2; k++)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      scaled_re_st * r * mesh_velocity[k] * dpsifdx(l2, k) *
                      testf[l] * W;
                  }
                }
              }


              // Axial velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * Gamma[0] * dpsifdx(l2, 0) * dtestfdx(l, 1) *
                  W;

                // Convective terms
                jacobian(local_eqn, local_unknown) -= scaled_re * r * psif[l2] *
                                                      interpolated_dudx(0, 1) *
                                                      testf[l] * W;
              }

              // Azimuthal velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              if (local_unknown >= 0)
              {
                // Convective terms
                jacobian(local_eqn, local_unknown) -= -scaled_re * 2.0 *
                                                      interpolated_u[2] *
                                                      psif[l2] * testf[l] * W;

                // Coriolis terms
                jacobian(local_eqn, local_unknown) +=
                  2.0 * r * scaled_re_inv_ro * psif[l2] * testf[l] * W;
              }
            }

            /*Now loop over pressure shape functions*/
            /*This is the contribution from pressure gradient*/
            for (unsigned l2 = 0; l2 < n_pres; l2++)
            {
              local_unknown = p_local_eqn(l2);
              /*If we are at a non-zero degree of freedom in the entry*/
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  psip[l2] * (testf[l] + r * dtestfdx(l, 0)) * W;
              }
            }
          } /*End of Jacobian calculation*/

        } // End of if not boundary condition statement

        // SECOND (AXIAL) MOMENTUM EQUATION
        local_eqn = nodal_local_eqn(l, u_nodal_index[1]);
        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Add the user-defined body force terms
          // Remember to multiply by the density ratio!
          residuals[local_eqn] += r * body_force[1] * testf[l] * W;

          // Add the gravitational body force term
          residuals[local_eqn] += r * scaled_re_inv_fr * testf[l] * G[1] * W;

          // Add the pressure gradient term
          residuals[local_eqn] += r * interpolated_p * dtestfdx(l, 1) * W;

          // Add in the stress tensor terms
          // The viscosity ratio needs to go in here to ensure
          // continuity of normal stress is satisfied even in flows
          // with zero pressure gradient!
          residuals[local_eqn] -=
            visc_ratio * r *
            (interpolated_dudx(1, 0) + Gamma[1] * interpolated_dudx(0, 1)) *
            dtestfdx(l, 0) * W;

          residuals[local_eqn] -= visc_ratio * r * (1.0 + Gamma[1]) *
                                  interpolated_dudx(1, 1) * dtestfdx(l, 1) * W;

          // Add in the inertial terms
          // du/dt term
          residuals[local_eqn] -= scaled_re_st * r * dudt[1] * testf[l] * W;

          // Convective terms
          residuals[local_eqn] -=
            scaled_re *
            (r * interpolated_u[0] * interpolated_dudx(1, 0) +
             r * interpolated_u[1] * interpolated_dudx(1, 1)) *
            testf[l] * W;

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            for (unsigned k = 0; k < 2; k++)
            {
              residuals[local_eqn] += scaled_re_st * r * mesh_velocity[k] *
                                      interpolated_dudx(1, k) * testf[l] * W;
            }
          }

          // CALCULATE THE JACOBIAN
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              // Radial velocity component
              if (local_unknown >= 0)
              {
                // Add in the stress tensor terms
                // The viscosity ratio needs to go in here to ensure
                // continuity of normal stress is satisfied even in flows
                // with zero pressure gradient!
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * Gamma[1] * dpsifdx(l2, 1) * dtestfdx(l, 0) *
                  W;

                // Convective terms
                jacobian(local_eqn, local_unknown) -= scaled_re * r * psif[l2] *
                                                      interpolated_dudx(1, 0) *
                                                      testf[l] * W;
              }

              // Axial velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  // Add the mass matrix
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] * testf[l] * W;
                }


                // Add in the stress tensor terms
                // The viscosity ratio needs to go in here to ensure
                // continuity of normal stress is satisfied even in flows
                // with zero pressure gradient!
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * dpsifdx(l2, 0) * dtestfdx(l, 0) * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * (1.0 + Gamma[1]) * dpsifdx(l2, 1) *
                  dtestfdx(l, 1) * W;

                // Add in the inertial terms
                // du/dt term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re_st * r *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif[l2] *
                  testf[l] * W;

                // Convective terms
                jacobian(local_eqn, local_unknown) -=
                  scaled_re *
                  (r * interpolated_u[0] * dpsifdx(l2, 0) +
                   r * psif[l2] * interpolated_dudx(1, 1) +
                   r * interpolated_u[1] * dpsifdx(l2, 1)) *
                  testf[l] * W;


                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned k = 0; k < 2; k++)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      scaled_re_st * r * mesh_velocity[k] * dpsifdx(l2, k) *
                      testf[l] * W;
                  }
                }
              }

              // There are no azimithal terms in the axial momentum equation
            } // End of loop over velocity shape functions

            /*Now loop over pressure shape functions*/
            /*This is the contribution from pressure gradient*/
            for (unsigned l2 = 0; l2 < n_pres; l2++)
            {
              local_unknown = p_local_eqn(l2);
              /*If we are at a non-zero degree of freedom in the entry*/
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  r * psip[l2] * dtestfdx(l, 1) * W;
              }
            }
          } /*End of Jacobian calculation*/

        } // End of AXIAL MOMENTUM EQUATION

        // THIRD (AZIMUTHAL) MOMENTUM EQUATION
        local_eqn = nodal_local_eqn(l, u_nodal_index[2]);
        if (local_eqn >= 0)
        {
          // Add the user-defined body force terms
          // Remember to multiply by the density ratio!
          residuals[local_eqn] += r * body_force[2] * testf[l] * W;

          // Add the gravitational body force term
          residuals[local_eqn] += r * scaled_re_inv_fr * testf[l] * G[2] * W;

          // There is NO pressure gradient term

          // Add in the stress tensor terms
          // The viscosity ratio needs to go in here to ensure
          // continuity of normal stress is satisfied even in flows
          // with zero pressure gradient!
          residuals[local_eqn] -=
            visc_ratio *
            (r * interpolated_dudx(2, 0) - Gamma[0] * interpolated_u[2]) *
            dtestfdx(l, 0) * W;

          residuals[local_eqn] -=
            visc_ratio * r * interpolated_dudx(2, 1) * dtestfdx(l, 1) * W;

          residuals[local_eqn] -=
            visc_ratio *
            ((interpolated_u[2] / r) - Gamma[0] * interpolated_dudx(2, 0)) *
            testf[l] * W;


          // Add in the inertial terms
          // du/dt term
          residuals[local_eqn] -= scaled_re_st * r * dudt[2] * testf[l] * W;

          // Convective terms
          residuals[local_eqn] -=
            scaled_re *
            (r * interpolated_u[0] * interpolated_dudx(2, 0) +
             interpolated_u[0] * interpolated_u[2] +
             r * interpolated_u[1] * interpolated_dudx(2, 1)) *
            testf[l] * W;

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            for (unsigned k = 0; k < 2; k++)
            {
              residuals[local_eqn] += scaled_re_st * r * mesh_velocity[k] *
                                      interpolated_dudx(2, k) * testf[l] * W;
            }
          }

          // Add the Coriolis term
          residuals[local_eqn] -=
            2.0 * r * scaled_re_inv_ro * interpolated_u[0] * testf[l] * W;

          // CALCULATE THE JACOBIAN
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Radial velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              if (local_unknown >= 0)
              {
                // Convective terms
                jacobian(local_eqn, local_unknown) -=
                  scaled_re *
                  (r * psif[l2] * interpolated_dudx(2, 0) +
                   psif[l2] * interpolated_u[2]) *
                  testf[l] * W;

                // Coriolis term
                jacobian(local_eqn, local_unknown) -=
                  2.0 * r * scaled_re_inv_ro * psif[l2] * testf[l] * W;
              }

              // Axial velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                // Convective terms
                jacobian(local_eqn, local_unknown) -= scaled_re * r * psif[l2] *
                                                      interpolated_dudx(2, 1) *
                                                      testf[l] * W;
              }

              // Azimuthal velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  // Add the mass matrix
                  mass_matrix(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] * testf[l] * W;
                }

                // Add in the stress tensor terms
                // The viscosity ratio needs to go in here to ensure
                // continuity of normal stress is satisfied even in flows
                // with zero pressure gradient!
                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * (r * dpsifdx(l2, 0) - Gamma[0] * psif[l2]) *
                  dtestfdx(l, 0) * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * r * dpsifdx(l2, 1) * dtestfdx(l, 1) * W;

                jacobian(local_eqn, local_unknown) -=
                  visc_ratio * ((psif[l2] / r) - Gamma[0] * dpsifdx(l2, 0)) *
                  testf[l] * W;

                // Add in the inertial terms
                // du/dt term
                jacobian(local_eqn, local_unknown) -=
                  scaled_re_st * r *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif[l2] *
                  testf[l] * W;

                // Convective terms
                jacobian(local_eqn, local_unknown) -=
                  scaled_re *
                  (r * interpolated_u[0] * dpsifdx(l2, 0) +
                   interpolated_u[0] * psif[l2] +
                   r * interpolated_u[1] * dpsifdx(l2, 1)) *
                  testf[l] * W;

                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned k = 0; k < 2; k++)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      scaled_re_st * r * mesh_velocity[k] * dpsifdx(l2, k) *
                      testf[l] * W;
                  }
                }
              }
            }

            // There are no pressure terms
          } // End of Jacobian

        } // End of AZIMUTHAL EQUATION

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
          residuals[local_eqn] -= r * source * testp[l] * W;

          // Gradient terms
          residuals[local_eqn] +=
            (interpolated_u[0] + r * interpolated_dudx(0, 0) +
             r * interpolated_dudx(1, 1)) *
            testp[l] * W;

          /*CALCULATE THE JACOBIAN*/
          if (flag)
          {
            /*Loop over the velocity shape functions*/
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Radial velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  (psif[l2] + r * dpsifdx(l2, 0)) * testp[l] * W;
              }

              // Axial velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  r * dpsifdx(l2, 1) * testp[l] * W;
              }
            } /*End of loop over l2*/
          } /*End of Jacobian calculation*/

        } // End of if not boundary condition

      } // End of loop over l
    }
  }


  //======================================================================
  /// Compute derivatives of elemental residual vector with respect
  /// to nodal coordinates.
  /// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
  /// Overloads the FD-based version in the FiniteElement base class.
  //======================================================================
  void AxisymmetricNavierStokesEquations::get_dresidual_dnodal_coordinates(
    RankThreeTensor<double>& dresidual_dnodal_coordinates)
  {
    // Return immediately if there are no dofs
    if (ndof() == 0)
    {
      return;
    }

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Determine number of nodes in element
    const unsigned n_node = nnode();

    // Determine number of pressure dofs in element
    const unsigned n_pres = npres_axi_nst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[3];
    for (unsigned i = 0; i < 3; i++)
    {
      u_nodal_index[i] = u_index_axi_nst(i);
    }

    // Set up memory for the shape and test functions
    // Note that there are only two dimensions, r and z, in this problem
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, 2), dtestfdx(n_node, 2);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Deriatives of shape fct derivatives w.r.t. nodal coords
    RankFourTensor<double> d_dpsifdx_dX(2, n_node, n_node, 2);
    RankFourTensor<double> d_dtestfdx_dX(2, n_node, n_node, 2);

    // Derivative of Jacobian of mapping w.r.t. to nodal coords
    DenseMatrix<double> dJ_dX(2, n_node);

    // Derivatives of derivative of u w.r.t. nodal coords
    // Note that the entry d_dudx_dX(p,q,i,k) contains the derivative w.r.t.
    // nodal coordinate X_pq of du_i/dx_k. Since there are three velocity
    // components, i can take on the values 0, 1 and 2, while k and p only
    // take on the values 0 and 1 since there are only two spatial dimensions.
    RankFourTensor<double> d_dudx_dX(2, n_node, 3, 2);

    // Derivatives of nodal velocities w.r.t. nodal coords:
    // Assumption: Interaction only local via no-slip so
    // X_pq only affects U_iq.
    // Note that the entry d_U_dX(p,q,i) contains the derivative w.r.t. nodal
    // coordinate X_pq of nodal value U_iq. In other words this entry is the
    // derivative of the i-th velocity component w.r.t. the p-th spatial
    // coordinate, taken at the q-th local node.
    RankThreeTensor<double> d_U_dX(2, n_node, 3, 0.0);

    // Determine the number of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Vector to hold local coordinates (two dimensions)
    Vector<double> s(2);

    // Get physical variables from element
    // (Reynolds number must be multiplied by the density ratio)
    const double scaled_re = re() * density_ratio();
    const double scaled_re_st = re_st() * density_ratio();
    const double scaled_re_inv_fr = re_invfr() * density_ratio();
    const double scaled_re_inv_ro = re_invro() * density_ratio();
    const double visc_ratio = viscosity_ratio();
    const Vector<double> G = g();

    // FD step
    double eps_fd = GeneralisedElement::Default_fd_jacobian_step;

    // Pre-compute derivatives of nodal velocities w.r.t. nodal coords:
    // Assumption: Interaction only local via no-slip so
    // X_ij only affects U_ij.
    bool element_has_node_with_aux_node_update_fct = false;
    for (unsigned q = 0; q < n_node; q++)
    {
      // Get pointer to q-th local node
      Node* nod_pt = node_pt(q);

      // Only compute if there's a node-update fct involved
      if (nod_pt->has_auxiliary_node_update_fct_pt())
      {
        element_has_node_with_aux_node_update_fct = true;

        // This functionality has not been tested so issue a warning
        // to this effect
        std::ostringstream warning_stream;
        warning_stream << "\nThe functionality to evaluate the additional"
                       << "\ncontribution to the deriv of the residual eqn"
                       << "\nw.r.t. the nodal coordinates which comes about"
                       << "\nif a node's values are updated using an auxiliary"
                       << "\nnode update function has NOT been tested for"
                       << "\naxisymmetric Navier-Stokes elements. Use at your"
                       << "\nown risk" << std::endl;
        OomphLibWarning(
          warning_stream.str(),
          "AxisymmetricNavierStokesEquations::get_dresidual_dnodal_coordinates",
          OOMPH_EXCEPTION_LOCATION);

        // Current nodal velocity
        Vector<double> u_ref(3);
        for (unsigned i = 0; i < 3; i++)
        {
          u_ref[i] = *(nod_pt->value_pt(u_nodal_index[i]));
        }

        // FD
        for (unsigned p = 0; p < 2; p++)
        {
          // Make backup
          double backup = nod_pt->x(p);

          // Do FD step. No node update required as we're
          // attacking the coordinate directly...
          nod_pt->x(p) += eps_fd;

          // Do auxiliary node update (to apply no slip)
          nod_pt->perform_auxiliary_node_update_fct();

          // Loop over velocity components
          for (unsigned i = 0; i < 3; i++)
          {
            // Evaluate
            d_U_dX(p, q, i) =
              (*(nod_pt->value_pt(u_nodal_index[i])) - u_ref[i]) / eps_fd;
          }

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
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      const double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      const double J = dshape_and_dtest_eulerian_at_knot_axi_nst(ipt,
                                                                 psif,
                                                                 dpsifdx,
                                                                 d_dpsifdx_dX,
                                                                 testf,
                                                                 dtestfdx,
                                                                 d_dtestfdx_dX,
                                                                 dJ_dX);

      // Call the pressure shape and test functions
      pshape_axi_nst(s, psip, testp);

      // Allocate storage for the position and the derivative of the
      // mesh positions w.r.t. time
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> mesh_velocity(2, 0.0);

      // Allocate storage for the pressure, velocity components and their
      // time and spatial derivatives
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(3, 0.0);
      Vector<double> dudt(3, 0.0);
      DenseMatrix<double> interpolated_dudx(3, 2, 0.0);

      // Calculate pressure at integration point
      for (unsigned l = 0; l < n_pres; l++)
      {
        interpolated_p += p_axi_nst(l) * psip[l];
      }

      // Calculate velocities and derivatives at integration point
      // ---------------------------------------------------------

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the shape function
        const double psif_ = psif(l);

        // Loop over the two coordinate directions
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += this->raw_nodal_position(l, i) * psif_;
        }

        // Loop over the three velocity directions
        for (unsigned i = 0; i < 3; i++)
        {
          // Get the nodal value
          const double u_value = this->raw_nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif_;
          dudt[i] += du_dt_axi_nst(l, i) * psif_;

          // Loop over derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      // Get the mesh velocity if ALE is enabled
      if (!ALE_is_disabled)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the two coordinate directions
          for (unsigned i = 0; i < 2; i++)
          {
            mesh_velocity[i] += this->raw_dnodal_position_dt(l, i) * psif[l];
          }
        }
      }

      // Calculate derivative of du_i/dx_k w.r.t. nodal positions X_{pq}
      for (unsigned q = 0; q < n_node; q++)
      {
        // Loop over the two coordinate directions
        for (unsigned p = 0; p < 2; p++)
        {
          // Loop over the three velocity components
          for (unsigned i = 0; i < 3; i++)
          {
            // Loop over the two coordinate directions
            for (unsigned k = 0; k < 2; k++)
            {
              double aux = 0.0;

              // Loop over nodes and add contribution
              for (unsigned j = 0; j < n_node; j++)
              {
                aux += this->raw_nodal_value(j, u_nodal_index[i]) *
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
      Vector<double> body_force(3);
      get_body_force_axi_nst(time, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      const double source = get_source_fct(time, ipt, interpolated_x);

      // Note: Can use raw values and avoid bypassing hanging information
      // because this is the non-refineable version!

      // Get gradient of body force function
      DenseMatrix<double> d_body_force_dx(3, 2, 0.0);
      get_body_force_gradient_axi_nst(
        time, ipt, s, interpolated_x, d_body_force_dx);

      // Get gradient of source function
      Vector<double> source_gradient(2, 0.0);
      get_source_fct_gradient(time, ipt, interpolated_x, source_gradient);

      // Cache r (the first position component)
      const double r = interpolated_x[0];

      // Assemble shape derivatives
      // --------------------------

      // ==================
      // MOMENTUM EQUATIONS
      // ==================

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the test function
        const double testf_ = testf[l];

        // --------------------------------
        // FIRST (RADIAL) MOMENTUM EQUATION
        // --------------------------------
        local_eqn = nodal_local_eqn(l, u_nodal_index[0]);
        ;

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Loop over the two coordinate directions
          for (unsigned p = 0; p < 2; p++)
          {
            // Loop over nodes
            for (unsigned q = 0; q < n_node; q++)
            {
              // Residual x deriv of Jacobian
              // ----------------------------

              // Add the user-defined body force terms
              double sum = r * body_force[0] * testf_;

              // Add the gravitational body force term
              sum += r * scaled_re_inv_fr * testf_ * G[0];

              // Add the pressure gradient term
              sum += interpolated_p * (testf_ + r * dtestfdx(l, 0));

              // Add the stress tensor terms
              // The viscosity ratio needs to go in here to ensure
              // continuity of normal stress is satisfied even in flows
              // with zero pressure gradient!
              sum -= visc_ratio * r * (1.0 + Gamma[0]) *
                     interpolated_dudx(0, 0) * dtestfdx(l, 0);

              sum -=
                visc_ratio * r *
                (interpolated_dudx(0, 1) + Gamma[0] * interpolated_dudx(1, 0)) *
                dtestfdx(l, 1);

              sum -=
                visc_ratio * (1.0 + Gamma[0]) * interpolated_u[0] * testf_ / r;

              // Add the du/dt term
              sum -= scaled_re_st * r * dudt[0] * testf_;

              // Add the convective terms
              sum -= scaled_re *
                     (r * interpolated_u[0] * interpolated_dudx(0, 0) -
                      interpolated_u[2] * interpolated_u[2] +
                      r * interpolated_u[1] * interpolated_dudx(0, 1)) *
                     testf_;

              // Add the mesh velocity terms
              if (!ALE_is_disabled)
              {
                for (unsigned k = 0; k < 2; k++)
                {
                  sum += scaled_re_st * r * mesh_velocity[k] *
                         interpolated_dudx(0, k) * testf_;
                }
              }

              // Add the Coriolis term
              sum += 2.0 * r * scaled_re_inv_ro * interpolated_u[2] * testf_;

              // Multiply through by deriv of Jacobian and integration weight
              dresidual_dnodal_coordinates(local_eqn, p, q) +=
                sum * dJ_dX(p, q) * w;

              // Derivative of residual x Jacobian
              // ---------------------------------

              // Body force
              sum = r * d_body_force_dx(0, p) * psif[q] * testf_;
              if (p == 0)
              {
                sum += body_force[0] * psif[q] * testf_;
              }

              // Gravitational body force
              if (p == 0)
              {
                sum += scaled_re_inv_fr * G[0] * psif[q] * testf_;
              }

              // Pressure gradient term
              sum += r * interpolated_p * d_dtestfdx_dX(p, q, l, 0);
              if (p == 0)
              {
                sum += interpolated_p * psif[q] * dtestfdx(l, 0);
              }

              // Viscous terms
              sum -=
                r * visc_ratio *
                ((1.0 + Gamma[0]) *
                   (d_dudx_dX(p, q, 0, 0) * dtestfdx(l, 0) +
                    interpolated_dudx(0, 0) * d_dtestfdx_dX(p, q, l, 0)) +
                 (d_dudx_dX(p, q, 0, 1) + Gamma[0] * d_dudx_dX(p, q, 1, 0)) *
                   dtestfdx(l, 1) +
                 (interpolated_dudx(0, 1) +
                  Gamma[0] * interpolated_dudx(1, 0)) *
                   d_dtestfdx_dX(p, q, l, 1));
              if (p == 0)
              {
                sum -= visc_ratio *
                       ((1.0 + Gamma[0]) *
                          (interpolated_dudx(0, 0) * psif[q] * dtestfdx(l, 0) -
                           interpolated_u[0] * psif[q] * testf_ / (r * r)) +
                        (interpolated_dudx(0, 1) +
                         Gamma[0] * interpolated_dudx(1, 0)) *
                          psif[q] * dtestfdx(l, 1));
              }

              // Convective terms, including mesh velocity
              for (unsigned k = 0; k < 2; k++)
              {
                double tmp = scaled_re * interpolated_u[k];
                if (!ALE_is_disabled)
                {
                  tmp -= scaled_re_st * mesh_velocity[k];
                }
                sum -= r * tmp * d_dudx_dX(p, q, 0, k) * testf_;
                if (p == 0)
                {
                  sum -= tmp * interpolated_dudx(0, k) * psif[q] * testf_;
                }
              }
              if (!ALE_is_disabled)
              {
                sum += scaled_re_st * r * pos_time_weight *
                       interpolated_dudx(0, p) * psif[q] * testf_;
              }

              // du/dt term
              if (p == 0)
              {
                sum -= scaled_re_st * dudt[0] * psif[q] * testf_;
              }

              // Coriolis term
              if (p == 0)
              {
                sum +=
                  2.0 * scaled_re_inv_ro * interpolated_u[2] * psif[q] * testf_;
              }

              // Multiply through by Jacobian and integration weight
              dresidual_dnodal_coordinates(local_eqn, p, q) += sum * J * w;

              // Derivs w.r.t. nodal velocities
              // ------------------------------
              if (element_has_node_with_aux_node_update_fct)
              {
                // Contribution from deriv of first velocity component
                double tmp =
                  scaled_re_st * r * val_time_weight * psif[q] * testf_;
                tmp +=
                  r * scaled_re * interpolated_dudx(0, 0) * psif[q] * testf_;
                for (unsigned k = 0; k < 2; k++)
                {
                  double tmp2 = scaled_re * interpolated_u[k];
                  if (!ALE_is_disabled)
                  {
                    tmp2 -= scaled_re_st * mesh_velocity[k];
                  }
                  tmp += r * tmp2 * dpsifdx(q, k) * testf_;
                }

                tmp += r * visc_ratio * (1.0 + Gamma[0]) * dpsifdx(q, 0) *
                       dtestfdx(l, 0);
                tmp += r * visc_ratio * dpsifdx(q, 1) * dtestfdx(l, 1);
                tmp += visc_ratio * (1.0 + Gamma[0]) * psif[q] * testf_ / r;

                // Multiply through by dU_0q/dX_pq
                sum = -d_U_dX(p, q, 0) * tmp;

                // Contribution from deriv of second velocity component
                tmp =
                  scaled_re * r * interpolated_dudx(0, 1) * psif[q] * testf_;
                tmp +=
                  r * visc_ratio * Gamma[0] * dpsifdx(q, 0) * dtestfdx(l, 1);

                // Multiply through by dU_1q/dX_pq
                sum -= d_U_dX(p, q, 1) * tmp;

                // Contribution from deriv of third velocity component
                tmp = 2.0 *
                      (r * scaled_re_inv_ro + scaled_re * interpolated_u[2]) *
                      psif[q] * testf_;

                // Multiply through by dU_2q/dX_pq
                sum += d_U_dX(p, q, 2) * tmp;

                // Multiply through by Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) += sum * J * w;
              }
            } // End of loop over q
          } // End of loop over p
        } // End of if it's not a boundary condition

        // --------------------------------
        // SECOND (AXIAL) MOMENTUM EQUATION
        // --------------------------------
        local_eqn = nodal_local_eqn(l, u_nodal_index[1]);
        ;

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Loop over the two coordinate directions
          for (unsigned p = 0; p < 2; p++)
          {
            // Loop over nodes
            for (unsigned q = 0; q < n_node; q++)
            {
              // Residual x deriv of Jacobian
              // ----------------------------

              // Add the user-defined body force terms
              double sum = r * body_force[1] * testf_;

              // Add the gravitational body force term
              sum += r * scaled_re_inv_fr * testf_ * G[1];

              // Add the pressure gradient term
              sum += r * interpolated_p * dtestfdx(l, 1);

              // Add the stress tensor terms
              // The viscosity ratio needs to go in here to ensure
              // continuity of normal stress is satisfied even in flows
              // with zero pressure gradient!
              sum -=
                visc_ratio * r *
                (interpolated_dudx(1, 0) + Gamma[1] * interpolated_dudx(0, 1)) *
                dtestfdx(l, 0);

              sum -= visc_ratio * r * (1.0 + Gamma[1]) *
                     interpolated_dudx(1, 1) * dtestfdx(l, 1);

              // Add the du/dt term
              sum -= scaled_re_st * r * dudt[1] * testf_;

              // Add the convective terms
              sum -= scaled_re *
                     (r * interpolated_u[0] * interpolated_dudx(1, 0) +
                      r * interpolated_u[1] * interpolated_dudx(1, 1)) *
                     testf_;

              // Add the mesh velocity terms
              if (!ALE_is_disabled)
              {
                for (unsigned k = 0; k < 2; k++)
                {
                  sum += scaled_re_st * r * mesh_velocity[k] *
                         interpolated_dudx(1, k) * testf_;
                }
              }

              // Multiply through by deriv of Jacobian and integration weight
              dresidual_dnodal_coordinates(local_eqn, p, q) +=
                sum * dJ_dX(p, q) * w;

              // Derivative of residual x Jacobian
              // ---------------------------------

              // Body force
              sum = r * d_body_force_dx(1, p) * psif[q] * testf_;
              if (p == 0)
              {
                sum += body_force[1] * psif[q] * testf_;
              }

              // Gravitational body force
              if (p == 0)
              {
                sum += scaled_re_inv_fr * G[1] * psif[q] * testf_;
              }

              // Pressure gradient term
              sum += r * interpolated_p * d_dtestfdx_dX(p, q, l, 1);
              if (p == 0)
              {
                sum += interpolated_p * psif[q] * dtestfdx(l, 1);
              }

              // Viscous terms
              sum -=
                r * visc_ratio *
                ((d_dudx_dX(p, q, 1, 0) + Gamma[1] * d_dudx_dX(p, q, 0, 1)) *
                   dtestfdx(l, 0) +
                 (interpolated_dudx(1, 0) +
                  Gamma[1] * interpolated_dudx(0, 1)) *
                   d_dtestfdx_dX(p, q, l, 0) +
                 (1.0 + Gamma[1]) * d_dudx_dX(p, q, 1, 1) * dtestfdx(l, 1) +
                 (1.0 + Gamma[1]) * interpolated_dudx(1, 1) *
                   d_dtestfdx_dX(p, q, l, 1));
              if (p == 0)
              {
                sum -=
                  visc_ratio * ((interpolated_dudx(1, 0) +
                                 Gamma[1] * interpolated_dudx(0, 1)) *
                                  psif[q] * dtestfdx(l, 0) +
                                (1.0 + Gamma[1]) * interpolated_dudx(1, 1) *
                                  psif[q] * dtestfdx(l, 1));
              }

              // Convective terms, including mesh velocity
              for (unsigned k = 0; k < 2; k++)
              {
                double tmp = scaled_re * interpolated_u[k];
                if (!ALE_is_disabled)
                {
                  tmp -= scaled_re_st * mesh_velocity[k];
                }
                sum -= r * tmp * d_dudx_dX(p, q, 1, k) * testf_;
                if (p == 0)
                {
                  sum -= tmp * interpolated_dudx(1, k) * psif[q] * testf_;
                }
              }
              if (!ALE_is_disabled)
              {
                sum += scaled_re_st * r * pos_time_weight *
                       interpolated_dudx(1, p) * psif[q] * testf_;
              }

              // du/dt term
              if (p == 0)
              {
                sum -= scaled_re_st * dudt[1] * psif[q] * testf_;
              }

              // Multiply through by Jacobian and integration weight
              dresidual_dnodal_coordinates(local_eqn, p, q) += sum * J * w;

              // Derivs w.r.t. to nodal velocities
              // ---------------------------------
              if (element_has_node_with_aux_node_update_fct)
              {
                // Contribution from deriv of first velocity component
                double tmp =
                  scaled_re * r * interpolated_dudx(1, 0) * psif[q] * testf_;
                tmp +=
                  r * visc_ratio * Gamma[1] * dpsifdx(q, 1) * dtestfdx(l, 0);

                // Multiply through by dU_0q/dX_pq
                sum = -d_U_dX(p, q, 0) * tmp;

                // Contribution from deriv of second velocity component
                tmp = scaled_re_st * r * val_time_weight * psif[q] * testf_;
                tmp +=
                  r * scaled_re * interpolated_dudx(1, 1) * psif[q] * testf_;
                for (unsigned k = 0; k < 2; k++)
                {
                  double tmp2 = scaled_re * interpolated_u[k];
                  if (!ALE_is_disabled)
                  {
                    tmp2 -= scaled_re_st * mesh_velocity[k];
                  }
                  tmp += r * tmp2 * dpsifdx(q, k) * testf_;
                }
                tmp += r * visc_ratio *
                       (dpsifdx(q, 0) * dtestfdx(l, 0) +
                        (1.0 + Gamma[1]) * dpsifdx(q, 1) * dtestfdx(l, 1));

                // Multiply through by dU_1q/dX_pq
                sum -= d_U_dX(p, q, 1) * tmp;

                // Multiply through by Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) += sum * J * w;
              }
            } // End of loop over q
          } // End of loop over p
        } // End of if it's not a boundary condition

        // -----------------------------------
        // THIRD (AZIMUTHAL) MOMENTUM EQUATION
        // -----------------------------------
        local_eqn = nodal_local_eqn(l, u_nodal_index[2]);
        ;

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Loop over the two coordinate directions
          for (unsigned p = 0; p < 2; p++)
          {
            // Loop over nodes
            for (unsigned q = 0; q < n_node; q++)
            {
              // Residual x deriv of Jacobian
              // ----------------------------

              // Add the user-defined body force terms
              double sum = r * body_force[2] * testf_;

              // Add the gravitational body force term
              sum += r * scaled_re_inv_fr * testf_ * G[2];

              // There is NO pressure gradient term

              // Add in the stress tensor terms
              // The viscosity ratio needs to go in here to ensure
              // continuity of normal stress is satisfied even in flows
              // with zero pressure gradient!
              sum -=
                visc_ratio *
                (r * interpolated_dudx(2, 0) - Gamma[0] * interpolated_u[2]) *
                dtestfdx(l, 0);

              sum -= visc_ratio * r * interpolated_dudx(2, 1) * dtestfdx(l, 1);

              sum -=
                visc_ratio *
                ((interpolated_u[2] / r) - Gamma[0] * interpolated_dudx(2, 0)) *
                testf_;

              // Add the du/dt term
              sum -= scaled_re_st * r * dudt[2] * testf_;

              // Add the convective terms
              sum -= scaled_re *
                     (r * interpolated_u[0] * interpolated_dudx(2, 0) +
                      interpolated_u[0] * interpolated_u[2] +
                      r * interpolated_u[1] * interpolated_dudx(2, 1)) *
                     testf_;

              // Add the mesh velocity terms
              if (!ALE_is_disabled)
              {
                for (unsigned k = 0; k < 2; k++)
                {
                  sum += scaled_re_st * r * mesh_velocity[k] *
                         interpolated_dudx(2, k) * testf_;
                }
              }

              // Add the Coriolis term
              sum -= 2.0 * r * scaled_re_inv_ro * interpolated_u[0] * testf_;

              // Multiply through by deriv of Jacobian and integration weight
              dresidual_dnodal_coordinates(local_eqn, p, q) +=
                sum * dJ_dX(p, q) * w;

              // Derivative of residual x Jacobian
              // ---------------------------------

              // Body force
              sum = r * d_body_force_dx(2, p) * psif[q] * testf_;
              if (p == 0)
              {
                sum += body_force[2] * psif[q] * testf_;
              }

              // Gravitational body force
              if (p == 0)
              {
                sum += scaled_re_inv_fr * G[2] * psif[q] * testf_;
              }

              // There is NO pressure gradient term

              // Viscous terms
              sum -= visc_ratio * r *
                     (d_dudx_dX(p, q, 2, 0) * dtestfdx(l, 0) +
                      interpolated_dudx(2, 0) * d_dtestfdx_dX(p, q, l, 0) +
                      d_dudx_dX(p, q, 2, 1) * dtestfdx(l, 1) +
                      interpolated_dudx(2, 1) * d_dtestfdx_dX(p, q, l, 1));

              sum += visc_ratio * Gamma[0] * d_dudx_dX(p, q, 2, 0) * testf_;

              if (p == 0)
              {
                sum -= visc_ratio *
                       (interpolated_dudx(2, 0) * psif[q] * dtestfdx(l, 0) +
                        interpolated_dudx(2, 1) * psif[q] * dtestfdx(l, 1) +
                        interpolated_u[2] * psif[q] * testf_ / (r * r));
              }

              // Convective terms, including mesh velocity
              for (unsigned k = 0; k < 2; k++)
              {
                double tmp = scaled_re * interpolated_u[k];
                if (!ALE_is_disabled)
                {
                  tmp -= scaled_re_st * mesh_velocity[k];
                }
                sum -= r * tmp * d_dudx_dX(p, q, 2, k) * testf_;
                if (p == 0)
                {
                  sum -= tmp * interpolated_dudx(2, k) * psif[q] * testf_;
                }
              }
              if (!ALE_is_disabled)
              {
                sum += scaled_re_st * r * pos_time_weight *
                       interpolated_dudx(2, p) * psif[q] * testf_;
              }

              // du/dt term
              if (p == 0)
              {
                sum -= scaled_re_st * dudt[2] * psif[q] * testf_;
              }

              // Coriolis term
              if (p == 0)
              {
                sum -=
                  2.0 * scaled_re_inv_ro * interpolated_u[0] * psif[q] * testf_;
              }

              // Multiply through by Jacobian and integration weight
              dresidual_dnodal_coordinates(local_eqn, p, q) += sum * J * w;

              // Derivs w.r.t. to nodal velocities
              // ---------------------------------
              if (element_has_node_with_aux_node_update_fct)
              {
                // Contribution from deriv of first velocity component
                double tmp = (2.0 * r * scaled_re_inv_ro +
                              r * scaled_re * interpolated_dudx(2, 0) -
                              scaled_re * interpolated_u[2]) *
                             psif[q] * testf_;

                // Multiply through by dU_0q/dX_pq
                sum = -d_U_dX(p, q, 0) * tmp;

                // Contribution from deriv of second velocity component
                // Multiply through by dU_1q/dX_pq
                sum -= d_U_dX(p, q, 1) * scaled_re * r *
                       interpolated_dudx(2, 1) * psif[q] * testf_;

                // Contribution from deriv of third velocity component
                tmp = scaled_re_st * r * val_time_weight * psif[q] * testf_;
                tmp -= scaled_re * interpolated_u[0] * psif[q] * testf_;
                for (unsigned k = 0; k < 2; k++)
                {
                  double tmp2 = scaled_re * interpolated_u[k];
                  if (!ALE_is_disabled)
                  {
                    tmp2 -= scaled_re_st * mesh_velocity[k];
                  }
                  tmp += r * tmp2 * dpsifdx(q, k) * testf_;
                }
                tmp += visc_ratio * (r * dpsifdx(q, 0) - Gamma[0] * psif[q]) *
                       dtestfdx(l, 0);
                tmp += r * visc_ratio * dpsifdx(q, 1) * dtestfdx(l, 1);
                tmp += visc_ratio * ((psif[q] / r) - Gamma[0] * dpsifdx(q, 0)) *
                       testf_;

                // Multiply through by dU_2q/dX_pq
                sum -= d_U_dX(p, q, 2) * tmp;

                // Multiply through by Jacobian and integration weight
                dresidual_dnodal_coordinates(local_eqn, p, q) += sum * J * w;
              }
            } // End of loop over q
          } // End of loop over p
        } // End of if it's not a boundary condition

      } // End of loop over test functions


      // ===================
      // CONTINUITY EQUATION
      // ===================

      // Loop over the shape functions
      for (unsigned l = 0; l < n_pres; l++)
      {
        local_eqn = p_local_eqn(l);

        // Cache the test function
        const double testp_ = testp[l];

        // If not a boundary conditions
        if (local_eqn >= 0)
        {
          // Loop over the two coordinate directions
          for (unsigned p = 0; p < 2; p++)
          {
            // Loop over nodes
            for (unsigned q = 0; q < n_node; q++)
            {
              // Residual x deriv of Jacobian
              //-----------------------------

              // Source term
              double aux = -r * source;

              // Gradient terms
              aux += (interpolated_u[0] + r * interpolated_dudx(0, 0) +
                      r * interpolated_dudx(1, 1));

              // Multiply through by deriv of Jacobian and integration weight
              dresidual_dnodal_coordinates(local_eqn, p, q) +=
                aux * dJ_dX(p, q) * testp_ * w;

              // Derivative of residual x Jacobian
              // ---------------------------------

              // Gradient of source function
              aux = -r * source_gradient[p] * psif[q];
              if (p == 0)
              {
                aux -= source * psif[q];
              }

              // Gradient terms
              aux += r * (d_dudx_dX(p, q, 0, 0) + d_dudx_dX(p, q, 1, 1));
              if (p == 0)
              {
                aux +=
                  (interpolated_dudx(0, 0) + interpolated_dudx(1, 1)) * psif[q];
              }

              // Derivs w.r.t. nodal velocities
              // ------------------------------
              if (element_has_node_with_aux_node_update_fct)
              {
                aux += d_U_dX(p, q, 0) * (psif[q] + r * dpsifdx(q, 0));
                aux += d_U_dX(p, q, 1) * r * dpsifdx(q, 1);
              }

              // Multiply through by Jacobian, test fct and integration weight
              dresidual_dnodal_coordinates(local_eqn, p, q) +=
                aux * testp_ * J * w;
            }
          }
        }
      } // End of loop over shape functions for continuity eqn

    } // End of loop over integration points
  }


  //==============================================================
  ///  Compute the residuals for the Navier--Stokes
  ///  equations; flag=1(or 0): do (or don't) compute the
  ///  Jacobian as well.
  //==============================================================
  void AxisymmetricNavierStokesEquations::
    fill_in_generic_dresidual_contribution_axi_nst(
      double* const& parameter_pt,
      Vector<double>& dres_dparam,
      DenseMatrix<double>& djac_dparam,
      DenseMatrix<double>& dmass_matrix_dparam,
      unsigned flag)
  {
    // Die if the parameter is not the Reynolds number
    if (parameter_pt != this->re_pt())
    {
      std::ostringstream error_stream;
      error_stream
        << "Cannot compute analytic jacobian for parameter addressed by "
        << parameter_pt << "\n";
      error_stream << "Can only compute derivatives wrt Re (" << Re_pt << ")\n";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Which parameters are we differentiating with respect to
    bool diff_re = false;
    bool diff_re_st = false;
    bool diff_re_inv_fr = false;
    bool diff_re_inv_ro = false;

    // Set the boolean flags according to the parameter pointer
    if (parameter_pt == this->re_pt())
    {
      diff_re = true;
    }
    if (parameter_pt == this->re_st_pt())
    {
      diff_re_st = true;
    }
    if (parameter_pt == this->re_invfr_pt())
    {
      diff_re_inv_fr = true;
    }
    if (parameter_pt == this->re_invro_pt())
    {
      diff_re_inv_ro = true;
    }


    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many pressure dofs there are
    unsigned n_pres = npres_axi_nst();

    // Get the nodal indices at which the velocity is stored
    unsigned u_nodal_index[3];
    for (unsigned i = 0; i < 3; ++i)
    {
      u_nodal_index[i] = u_index_axi_nst(i);
    }

    // Set up memory for the shape and test functions
    // Note that there are only two dimensions, r and z in this problem
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, 2), dtestfdx(n_node, 2);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Number of integration points
    unsigned Nintpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates (two dimensions)
    Vector<double> s(2);

    // Get Physical Variables from Element
    // Reynolds number must be multiplied by the density ratio
    // double scaled_re = re()*density_ratio();
    // double scaled_re_st = re_st()*density_ratio();
    // double scaled_re_inv_fr = re_invfr()*density_ratio();
    // double scaled_re_inv_ro = re_invro()*density_ratio();
    double dens_ratio = this->density_ratio();
    // double visc_ratio = viscosity_ratio();
    Vector<double> G = g();

    // Integers used to store the local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < Nintpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++) s[i] = integral_pt()->knot(ipt, i);
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_axi_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Call the pressure shape and test functions
      pshape_axi_nst(s, psip, testp);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Allocate storage for the position and the derivative of the
      // mesh positions wrt time
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> mesh_velocity(2, 0.0);
      // Allocate storage for the pressure, velocity components and their
      // time and spatial derivatives
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(3, 0.0);
      Vector<double> dudt(3, 0.0);
      DenseMatrix<double> interpolated_dudx(3, 2, 0.0);

      // Calculate pressure at integration point
      for (unsigned l = 0; l < n_pres; l++)
      {
        interpolated_p += p_axi_nst(l) * psip[l];
      }

      // Calculate velocities and derivatives at integration point

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the shape function
        const double psif_ = psif(l);
        // Loop over the two coordinate directions
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += this->raw_nodal_position(l, i) * psif_;
        }
        // mesh_velocity[i]  += dnodal_position_dt(l,i)*psif[l];

        // Loop over the three velocity directions
        for (unsigned i = 0; i < 3; i++)
        {
          // Get the u_value
          const double u_value = this->raw_nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif_;
          dudt[i] += du_dt_axi_nst(l, i) * psif_;
          // Loop over derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      // Get the mesh velocity if ALE is enabled
      if (!ALE_is_disabled)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the two coordinate directions
          for (unsigned i = 0; i < 2; i++)
          {
            mesh_velocity[i] += this->raw_dnodal_position_dt(l, i) * psif(l);
          }
        }
      }


      // Get the user-defined body force terms
      // Vector<double> body_force(3);
      // get_body_force(time(),ipt,interpolated_x,body_force);

      // Get the user-defined source function
      // double source = get_source_fct(time(),ipt,interpolated_x);

      // Get the user-defined viscosity function
      double visc_ratio;
      get_viscosity_ratio_axisym_nst(ipt, s, interpolated_x, visc_ratio);

      // r is the first position component
      double r = interpolated_x[0];


      // MOMENTUM EQUATIONS
      //------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // FIRST (RADIAL) MOMENTUM EQUATION
        local_eqn = nodal_local_eqn(l, u_nodal_index[0]);
        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // No user-defined body force terms
          // dres_dparam[local_eqn] +=
          // r*body_force[0]*testf[l]*W;

          // Add the gravitational body force term if the reynolds number
          // is equal to re_inv_fr
          if (diff_re_inv_fr)
          {
            dres_dparam[local_eqn] += r * dens_ratio * testf[l] * G[0] * W;
          }

          // No pressure gradient term
          // residuals[local_eqn]  +=
          // interpolated_p*(testf[l] + r*dtestfdx(l,0))*W;

          // No in the stress tensor terms
          // The viscosity ratio needs to go in here to ensure
          // continuity of normal stress is satisfied even in flows
          // with zero pressure gradient!
          // residuals[local_eqn] -= visc_ratio*
          // r*(1.0+Gamma[0])*interpolated_dudx(0,0)*dtestfdx(l,0)*W;

          // residuals[local_eqn] -= visc_ratio*r*
          // (interpolated_dudx(0,1) + Gamma[0]*interpolated_dudx(1,0))*
          // dtestfdx(l,1)*W;

          // residuals[local_eqn] -=
          // visc_ratio*(1.0 + Gamma[0])*interpolated_u[0]*testf[l]*W/r;

          // Add in the inertial terms
          // du/dt term
          if (diff_re_st)
          {
            dres_dparam[local_eqn] -= dens_ratio * r * dudt[0] * testf[l] * W;
          }

          // Convective terms
          if (diff_re)
          {
            dres_dparam[local_eqn] -=
              dens_ratio *
              (r * interpolated_u[0] * interpolated_dudx(0, 0) -
               interpolated_u[2] * interpolated_u[2] +
               r * interpolated_u[1] * interpolated_dudx(0, 1)) *
              testf[l] * W;
          }

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            if (diff_re_st)
            {
              for (unsigned k = 0; k < 2; k++)
              {
                dres_dparam[local_eqn] += dens_ratio * r * mesh_velocity[k] *
                                          interpolated_dudx(0, k) * testf[l] *
                                          W;
              }
            }
          }

          // Add the Coriolis term
          if (diff_re_inv_ro)
          {
            dres_dparam[local_eqn] +=
              2.0 * r * dens_ratio * interpolated_u[2] * testf[l] * W;
          }

          // CALCULATE THE JACOBIAN
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              // Radial velocity component
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  if (diff_re_st)
                  {
                    // Add the mass matrix
                    dmass_matrix_dparam(local_eqn, local_unknown) +=
                      dens_ratio * r * psif[l2] * testf[l] * W;
                  }
                }

                // Add contribution to the Jacobian matrix
                // jacobian(local_eqn,local_unknown)
                // -= visc_ratio*r*(1.0+Gamma[0])
                //*dpsifdx(l2,0)*dtestfdx(l,0)*W;

                // jacobian(local_eqn,local_unknown)
                // -= visc_ratio*r*dpsifdx(l2,1)*dtestfdx(l,1)*W;

                // jacobian(local_eqn,local_unknown)
                // -= visc_ratio*(1.0 + Gamma[0])*psif[l2]*testf[l]*W/r;

                // Add in the inertial terms
                // du/dt term
                if (diff_re_st)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    dens_ratio * r *
                    node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif[l2] *
                    testf[l] * W;
                }

                // Convective terms
                if (diff_re)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    dens_ratio *
                    (r * psif[l2] * interpolated_dudx(0, 0) +
                     r * interpolated_u[0] * dpsifdx(l2, 0) +
                     r * interpolated_u[1] * dpsifdx(l2, 1)) *
                    testf[l] * W;
                }

                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned k = 0; k < 2; k++)
                  {
                    if (diff_re_st)
                    {
                      djac_dparam(local_eqn, local_unknown) +=
                        dens_ratio * r * mesh_velocity[k] * dpsifdx(l2, k) *
                        testf[l] * W;
                    }
                  }
                }
              }

              // Axial velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                // jacobian(local_eqn,local_unknown) -=
                // visc_ratio*r*Gamma[0]*dpsifdx(l2,0)*dtestfdx(l,1)*W;

                // Convective terms
                if (diff_re)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    dens_ratio * r * psif[l2] * interpolated_dudx(0, 1) *
                    testf[l] * W;
                }
              }

              // Azimuthal velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              if (local_unknown >= 0)
              {
                // Convective terms
                if (diff_re)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    -dens_ratio * 2.0 * interpolated_u[2] * psif[l2] *
                    testf[l] * W;
                }
                // Coriolis terms
                if (diff_re_inv_ro)
                {
                  djac_dparam(local_eqn, local_unknown) +=
                    2.0 * r * dens_ratio * psif[l2] * testf[l] * W;
                }
              }
            }

            /*Now loop over pressure shape functions*/
            /*This is the contribution from pressure gradient*/
            // for(unsigned l2=0;l2<n_pres;l2++)
            // {
            //  local_unknown = p_local_eqn(l2);
            //  /*If we are at a non-zero degree of freedom in the entry*/
            //  if(local_unknown >= 0)
            //   {
            //    jacobian(local_eqn,local_unknown)
            //     += psip[l2]*(testf[l] + r*dtestfdx(l,0))*W;
            //   }
            // }
          } /*End of Jacobian calculation*/

        } // End of if not boundary condition statement

        // SECOND (AXIAL) MOMENTUM EQUATION
        local_eqn = nodal_local_eqn(l, u_nodal_index[1]);
        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Add the user-defined body force terms
          // Remember to multiply by the density ratio!
          // residuals[local_eqn] +=
          // r*body_force[1]*testf[l]*W;

          // Add the gravitational body force term
          if (diff_re_inv_fr)
          {
            dres_dparam[local_eqn] += r * dens_ratio * testf[l] * G[1] * W;
          }

          // Add the pressure gradient term
          // residuals[local_eqn]  += r*interpolated_p*dtestfdx(l,1)*W;

          // Add in the stress tensor terms
          // The viscosity ratio needs to go in here to ensure
          // continuity of normal stress is satisfied even in flows
          // with zero pressure gradient!
          // residuals[local_eqn] -= visc_ratio*
          //  r*(interpolated_dudx(1,0) + Gamma[1]*interpolated_dudx(0,1))
          // *dtestfdx(l,0)*W;

          // residuals[local_eqn] -= visc_ratio*r*
          // (1.0 + Gamma[1])*interpolated_dudx(1,1)*dtestfdx(l,1)*W;

          // Add in the inertial terms
          // du/dt term
          if (diff_re_st)
          {
            dres_dparam[local_eqn] -= dens_ratio * r * dudt[1] * testf[l] * W;
          }

          // Convective terms
          if (diff_re)
          {
            dres_dparam[local_eqn] -=
              dens_ratio *
              (r * interpolated_u[0] * interpolated_dudx(1, 0) +
               r * interpolated_u[1] * interpolated_dudx(1, 1)) *
              testf[l] * W;
          }

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            if (diff_re_st)
            {
              for (unsigned k = 0; k < 2; k++)
              {
                dres_dparam[local_eqn] += dens_ratio * r * mesh_velocity[k] *
                                          interpolated_dudx(1, k) * testf[l] *
                                          W;
              }
            }
          }

          // CALCULATE THE JACOBIAN
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              // Radial velocity component
              if (local_unknown >= 0)
              {
                // Add in the stress tensor terms
                // The viscosity ratio needs to go in here to ensure
                // continuity of normal stress is satisfied even in flows
                // with zero pressure gradient!
                // jacobian(local_eqn,local_unknown) -=
                // visc_ratio*r*Gamma[1]*dpsifdx(l2,1)*dtestfdx(l,0)*W;

                // Convective terms
                if (diff_re)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    dens_ratio * r * psif[l2] * interpolated_dudx(1, 0) *
                    testf[l] * W;
                }
              }

              // Axial velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  if (diff_re_st)
                  {
                    // Add the mass matrix
                    dmass_matrix_dparam(local_eqn, local_unknown) +=
                      dens_ratio * r * psif[l2] * testf[l] * W;
                  }
                }


                // Add in the stress tensor terms
                // The viscosity ratio needs to go in here to ensure
                // continuity of normal stress is satisfied even in flows
                // with zero pressure gradient!
                // jacobian(local_eqn,local_unknown) -=
                // visc_ratio*r*dpsifdx(l2,0)*dtestfdx(l,0)*W;

                // jacobian(local_eqn,local_unknown) -=
                // visc_ratio*r*(1.0 + Gamma[1])*dpsifdx(l2,1)*
                // dtestfdx(l,1)*W;

                // Add in the inertial terms
                // du/dt term
                if (diff_re_st)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    dens_ratio * r *
                    node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif[l2] *
                    testf[l] * W;
                }

                // Convective terms
                if (diff_re)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    dens_ratio *
                    (r * interpolated_u[0] * dpsifdx(l2, 0) +
                     r * psif[l2] * interpolated_dudx(1, 1) +
                     r * interpolated_u[1] * dpsifdx(l2, 1)) *
                    testf[l] * W;
                }

                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  for (unsigned k = 0; k < 2; k++)
                  {
                    if (diff_re_st)
                    {
                      djac_dparam(local_eqn, local_unknown) +=
                        dens_ratio * r * mesh_velocity[k] * dpsifdx(l2, k) *
                        testf[l] * W;
                    }
                  }
                }
              }

              // There are no azimithal terms in the axial momentum equation
            } // End of loop over velocity shape functions

          } /*End of Jacobian calculation*/

        } // End of AXIAL MOMENTUM EQUATION

        // THIRD (AZIMUTHAL) MOMENTUM EQUATION
        local_eqn = nodal_local_eqn(l, u_nodal_index[2]);
        if (local_eqn >= 0)
        {
          // Add the user-defined body force terms
          // Remember to multiply by the density ratio!
          // residuals[local_eqn] +=
          // r*body_force[2]*testf[l]*W;

          // Add the gravitational body force term
          if (diff_re_inv_fr)
          {
            dres_dparam[local_eqn] += r * dens_ratio * testf[l] * G[2] * W;
          }

          // There is NO pressure gradient term

          // Add in the stress tensor terms
          // The viscosity ratio needs to go in here to ensure
          // continuity of normal stress is satisfied even in flows
          // with zero pressure gradient!
          // residuals[local_eqn] -= visc_ratio*
          // (r*interpolated_dudx(2,0) -
          //  Gamma[0]*interpolated_u[2])*dtestfdx(l,0)*W;

          // residuals[local_eqn] -= visc_ratio*r*
          // interpolated_dudx(2,1)*dtestfdx(l,1)*W;

          // residuals[local_eqn] -= visc_ratio*
          // ((interpolated_u[2]/r) -
          // Gamma[0]*interpolated_dudx(2,0))*testf[l]*W;


          // Add in the inertial terms
          // du/dt term
          if (diff_re_st)
          {
            dres_dparam[local_eqn] -= dens_ratio * r * dudt[2] * testf[l] * W;
          }

          // Convective terms
          if (diff_re)
          {
            dres_dparam[local_eqn] -=
              dens_ratio *
              (r * interpolated_u[0] * interpolated_dudx(2, 0) +
               interpolated_u[0] * interpolated_u[2] +
               r * interpolated_u[1] * interpolated_dudx(2, 1)) *
              testf[l] * W;
          }

          // Mesh velocity terms
          if (!ALE_is_disabled)
          {
            if (diff_re_st)
            {
              for (unsigned k = 0; k < 2; k++)
              {
                dres_dparam[local_eqn] += dens_ratio * r * mesh_velocity[k] *
                                          interpolated_dudx(2, k) * testf[l] *
                                          W;
              }
            }
          }

          // Add the Coriolis term
          if (diff_re_inv_ro)
          {
            dres_dparam[local_eqn] -=
              2.0 * r * dens_ratio * interpolated_u[0] * testf[l] * W;
          }

          // CALCULATE THE JACOBIAN
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Radial velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              if (local_unknown >= 0)
              {
                // Convective terms
                if (diff_re)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    dens_ratio *
                    (r * psif[l2] * interpolated_dudx(2, 0) +
                     psif[l2] * interpolated_u[2]) *
                    testf[l] * W;
                }

                // Coriolis term
                if (diff_re_inv_ro)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    2.0 * r * dens_ratio * psif[l2] * testf[l] * W;
                }
              }

              // Axial velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              if (local_unknown >= 0)
              {
                // Convective terms
                if (diff_re)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    dens_ratio * r * psif[l2] * interpolated_dudx(2, 1) *
                    testf[l] * W;
                }
              }

              // Azimuthal velocity component
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              if (local_unknown >= 0)
              {
                if (flag == 2)
                {
                  // Add the mass matrix
                  if (diff_re_st)
                  {
                    dmass_matrix_dparam(local_eqn, local_unknown) +=
                      dens_ratio * r * psif[l2] * testf[l] * W;
                  }
                }

                // Add in the stress tensor terms
                // The viscosity ratio needs to go in here to ensure
                // continuity of normal stress is satisfied even in flows
                // with zero pressure gradient!
                // jacobian(local_eqn,local_unknown) -=
                // visc_ratio*(r*dpsifdx(l2,0) -
                //                  Gamma[0]*psif[l2])*dtestfdx(l,0)*W;

                // jacobian(local_eqn,local_unknown) -=
                // visc_ratio*r*dpsifdx(l2,1)*dtestfdx(l,1)*W;

                // jacobian(local_eqn,local_unknown) -=
                // visc_ratio*((psif[l2]/r) - Gamma[0]*dpsifdx(l2,0))
                // *testf[l]*W;

                // Add in the inertial terms
                // du/dt term
                if (diff_re_st)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    dens_ratio * r *
                    node_pt(l2)->time_stepper_pt()->weight(1, 0) * psif[l2] *
                    testf[l] * W;
                }

                // Convective terms
                if (diff_re)
                {
                  djac_dparam(local_eqn, local_unknown) -=
                    dens_ratio *
                    (r * interpolated_u[0] * dpsifdx(l2, 0) +
                     interpolated_u[0] * psif[l2] +
                     r * interpolated_u[1] * dpsifdx(l2, 1)) *
                    testf[l] * W;
                }

                // Mesh velocity terms
                if (!ALE_is_disabled)
                {
                  if (diff_re_st)
                  {
                    for (unsigned k = 0; k < 2; k++)
                    {
                      djac_dparam(local_eqn, local_unknown) +=
                        dens_ratio * r * mesh_velocity[k] * dpsifdx(l2, k) *
                        testf[l] * W;
                    }
                  }
                }
              }
            }

            // There are no pressure terms
          } // End of Jacobian

        } // End of AZIMUTHAL EQUATION

      } // End of loop over shape functions


      // CONTINUITY EQUATION NO PARAMETERS
      //-------------------
    }
  }

  //=========================================================================
  /// \short Compute the hessian tensor vector products required to
  /// perform continuation of bifurcations analytically
  //=========================================================================
  void AxisymmetricNavierStokesEquations::
    fill_in_contribution_to_hessian_vector_products(
      Vector<double> const& Y,
      DenseMatrix<double> const& C,
      DenseMatrix<double>& product)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Get the nodal indices at which the velocity is stored
    unsigned u_nodal_index[3];
    for (unsigned i = 0; i < 3; ++i)
    {
      u_nodal_index[i] = u_index_axi_nst(i);
    }

    // Set up memory for the shape and test functions
    // Note that there are only two dimensions, r and z in this problem
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, 2), dtestfdx(n_node, 2);

    // Number of integration points
    unsigned Nintpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates (two dimensions)
    Vector<double> s(2);

    // Get Physical Variables from Element
    // Reynolds number must be multiplied by the density ratio
    double scaled_re = re() * density_ratio();
    // double visc_ratio = viscosity_ratio();
    Vector<double> G = g();

    // Integers used to store the local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0, local_freedom = 0;

    // How may dofs are there
    const unsigned n_dof = this->ndof();

    // Create a local matrix eigenvector product contribution
    DenseMatrix<double> jac_y(n_dof, n_dof, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < Nintpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 2; i++) s[i] = integral_pt()->knot(ipt, i);
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_axi_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Allocate storage for the position and the derivative of the
      // mesh positions wrt time
      Vector<double> interpolated_x(2, 0.0);
      // Vector<double> mesh_velocity(2,0.0);
      // Allocate storage for the pressure, velocity components and their
      // time and spatial derivatives
      Vector<double> interpolated_u(3, 0.0);
      // Vector<double> dudt(3,0.0);
      DenseMatrix<double> interpolated_dudx(3, 2, 0.0);

      // Calculate velocities and derivatives at integration point

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Cache the shape function
        const double psif_ = psif(l);
        // Loop over the two coordinate directions
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += this->raw_nodal_position(l, i) * psif_;
        }

        // Loop over the three velocity directions
        for (unsigned i = 0; i < 3; i++)
        {
          // Get the u_value
          const double u_value = this->raw_nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif_;
          // dudt[i]+= du_dt_axi_nst(l,i)*psif_;
          // Loop over derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      // Get the mesh velocity if ALE is enabled
      if (!ALE_is_disabled)
      {
        throw OomphLibError("Moving nodes not implemented\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // r is the first position component
      double r = interpolated_x[0];


      // MOMENTUM EQUATIONS
      //------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // FIRST (RADIAL) MOMENTUM EQUATION
        local_eqn = nodal_local_eqn(l, u_nodal_index[0]);
        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Loop over the velocity shape functions yet again
          for (unsigned l3 = 0; l3 < n_node; l3++)
          {
            // Derivative of jacobian terms with respect to radial velocity
            local_freedom = nodal_local_eqn(l3, u_nodal_index[0]);
            if (local_freedom >= 0)
            {
              // Storage for the sums
              double temp = 0.0;

              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
                // Radial velocity component
                if (local_unknown >= 0)
                {
                  // Remains of convective terms
                  temp -= scaled_re *
                          (r * psif[l2] * dpsifdx(l3, 0) +
                           r * psif[l3] * dpsifdx(l2, 0)) *
                          Y[local_unknown] * testf[l] * W;
                }

                // Axial velocity component
                local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
                if (local_unknown >= 0)
                {
                  // Convective terms
                  temp -= scaled_re * r * psif[l2] * dpsifdx(l3, 1) *
                          Y[local_unknown] * testf[l] * W;
                }
              }
              // Add the summed contribution to the product matrix
              jac_y(local_eqn, local_freedom) += temp;
            } // End of derivative wrt radial coordinate


            // Derivative of jacobian terms with respect to axial velocity
            local_freedom = nodal_local_eqn(l3, u_nodal_index[1]);
            if (local_freedom >= 0)
            {
              double temp = 0.0;

              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
                // Radial velocity component
                if (local_unknown >= 0)
                {
                  // Remains of convective terms
                  temp -= scaled_re * (r * psif[l3] * dpsifdx(l2, 1)) *
                          Y[local_unknown] * testf[l] * W;
                }
              }
              // Add the summed contribution to the product matrix
              jac_y(local_eqn, local_freedom) += temp;
            } // End of derivative wrt axial coordinate

            // Derivative of jacobian terms with respect to azimuthal velocity
            local_freedom = nodal_local_eqn(l3, u_nodal_index[2]);
            if (local_freedom >= 0)
            {
              double temp = 0.0;

              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Azimuthal velocity component
                local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
                if (local_unknown >= 0)
                {
                  // Convective terms
                  temp -= -scaled_re * 2.0 * psif[l3] * psif[l2] *
                          Y[local_unknown] * testf[l] * W;
                }
              }
              // Add the summed contibution
              jac_y(local_eqn, local_freedom) += temp;

            } // End of if not boundary condition statement
          } // End of loop over freedoms
        } // End of RADIAL MOMENTUM EQUATION


        // SECOND (AXIAL) MOMENTUM EQUATION
        local_eqn = nodal_local_eqn(l, u_nodal_index[1]);
        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Loop over the velocity shape functions yet again
          for (unsigned l3 = 0; l3 < n_node; l3++)
          {
            // Derivative of jacobian terms with respect to radial velocity
            local_freedom = nodal_local_eqn(l3, u_nodal_index[0]);
            if (local_freedom >= 0)
            {
              double temp = 0.0;

              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Axial velocity component
                local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
                if (local_unknown >= 0)
                {
                  // Convective terms
                  temp -= scaled_re * (r * psif[l3] * dpsifdx(l2, 0)) *
                          Y[local_unknown] * testf[l] * W;
                }
              }
              jac_y(local_eqn, local_freedom) += temp;

              // There are no azimithal terms in the axial momentum equation
            } // End of loop over velocity shape functions


            // Derivative of jacobian terms with respect to axial velocity
            local_freedom = nodal_local_eqn(l3, u_nodal_index[1]);
            if (local_freedom >= 0)
            {
              double temp = 0.0;

              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
                // Radial velocity component
                if (local_unknown >= 0)
                {
                  // Convective terms
                  temp -= scaled_re * r * psif[l2] * dpsifdx(l3, 0) *
                          Y[local_unknown] * testf[l] * W;
                }

                // Axial velocity component
                local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
                if (local_unknown >= 0)
                {
                  // Convective terms
                  temp -= scaled_re *
                          (r * psif[l2] * dpsifdx(l3, 1) +
                           r * psif[l3] * dpsifdx(l2, 1)) *
                          Y[local_unknown] * testf[l] * W;
                }

                // There are no azimithal terms in the axial momentum equation
              } // End of loop over velocity shape functions

              // Add summed contributiont to jacobian product matrix
              jac_y(local_eqn, local_freedom) += temp;
            }
          } // End of loop over local freedoms

        } // End of AXIAL MOMENTUM EQUATION

        // THIRD (AZIMUTHAL) MOMENTUM EQUATION
        local_eqn = nodal_local_eqn(l, u_nodal_index[2]);
        if (local_eqn >= 0)
        {
          // Loop over the velocity shape functions yet again
          for (unsigned l3 = 0; l3 < n_node; l3++)
          {
            // Derivative of jacobian terms with respect to radial velocity
            local_freedom = nodal_local_eqn(l3, u_nodal_index[0]);
            if (local_freedom >= 0)
            {
              double temp = 0.0;

              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Azimuthal velocity component
                local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
                if (local_unknown >= 0)
                {
                  // Convective terms
                  temp -=
                    scaled_re *
                    (r * psif[l3] * dpsifdx(l2, 0) + psif[l3] * psif[l2]) *
                    Y[local_unknown] * testf[l] * W;
                }
              }
              // Add the summed contribution to the jacobian eigenvector sum
              jac_y(local_eqn, local_freedom) += temp;
            }

            // Derivative of jacobian terms with respect to axial velocity
            local_freedom = nodal_local_eqn(l3, u_nodal_index[1]);
            if (local_freedom >= 0)
            {
              double temp = 0.0;

              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Azimuthal velocity component
                local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
                if (local_unknown >= 0)
                {
                  // Convective terms
                  temp -= scaled_re * (r * psif[l3] * dpsifdx(l2, 1)) *
                          Y[local_unknown] * testf[l] * W;
                }
              }
              // Add the summed contribution to the jacobian eigenvector sum
              jac_y(local_eqn, local_freedom) += temp;
            }


            // Derivative of jacobian terms with respect to azimuthal velocity
            local_freedom = nodal_local_eqn(l3, u_nodal_index[2]);
            if (local_freedom >= 0)
            {
              double temp = 0.0;


              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Radial velocity component
                local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
                if (local_unknown >= 0)
                {
                  // Convective terms
                  temp -=
                    scaled_re *
                    (r * psif[l2] * dpsifdx(l3, 0) + psif[l2] * psif[l3]) *
                    Y[local_unknown] * testf[l] * W;
                }

                // Axial velocity component
                local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
                if (local_unknown >= 0)
                {
                  // Convective terms
                  temp -= scaled_re * r * psif[l2] * dpsifdx(l3, 1) *
                          Y[local_unknown] * testf[l] * W;
                }
              }
              // Add the fully summed contribution
              jac_y(local_eqn, local_freedom) += temp;
            }
          } // End of loop over freedoms

          // There are no pressure terms
        } // End of AZIMUTHAL EQUATION

      } // End of loop over shape functions
    }

    // We have now assembled the matrix (d J_{ij} Y_j)/d u_{k}
    // and simply need to sum over the vectors
    const unsigned n_vec = C.nrow();
    for (unsigned i = 0; i < n_dof; i++)
    {
      for (unsigned k = 0; k < n_dof; k++)
      {
        // Cache the value of the hessian y product
        const double j_y = jac_y(i, k);
        // Loop over the possible vectors
        for (unsigned v = 0; v < n_vec; v++)
        {
          product(v, i) += j_y * C(v, k);
        }
      }
    }
  }


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  //=============================================================================
  /// Create a list of pairs for all unknowns in this element,
  /// so that the first entry in each pair contains the global equation
  /// number of the unknown, while the second one contains the number
  /// of the "DOF type" that this unknown is associated with.
  /// (Function can obviously only be called if the equation numbering
  /// scheme has been set up.)
  //=============================================================================
  void AxisymmetricQCrouzeixRaviartElement::get_dof_numbers_for_unknowns(
    std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // number of nodes
    unsigned n_node = this->nnode();

    // number of pressure values
    unsigned n_press = this->npres_axi_nst();

    // temporary pair (used to store dof lookup prior to being added to list)
    std::pair<unsigned, unsigned> dof_lookup;

    // pressure dof number (is this really OK?)
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


  /// Axisymmetric Crouzeix-Raviart elements
  // Set the data for the number of Variables at each node
  const unsigned AxisymmetricQCrouzeixRaviartElement::Initial_Nvalue[9] = {
    3, 3, 3, 3, 3, 3, 3, 3, 3};

  //========================================================================
  /// Number of values (pinned or dofs) required at node n.
  //========================================================================
  unsigned AxisymmetricQCrouzeixRaviartElement::required_nvalue(
    const unsigned& n) const
  {
    return Initial_Nvalue[n];
  }

  //========================================================================
  /// Compute traction at local coordinate s for outer unit normal N
  //========================================================================
  void AxisymmetricQCrouzeixRaviartElement::get_traction(
    const Vector<double>& s,
    const Vector<double>& N,
    Vector<double>& traction) const
  {
    AxisymmetricNavierStokesEquations::traction(s, N, traction);
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  //============================================================================
  /// Create a list of pairs for all unknowns in this element,
  /// so the first entry in each pair contains the global equation
  /// number of the unknown, while the second one contains the number
  /// of the "DOF type" that this unknown is associated with.
  /// (Function can obviously only be called if the equation numbering
  /// scheme has been set up.)
  //============================================================================
  void AxisymmetricQTaylorHoodElement::get_dof_numbers_for_unknowns(
    std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
  {
    // number of nodes
    unsigned n_node = this->nnode();

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


  // Axisymmetric Taylor--Hood
  // Set the data for the number of Variables at each node
  const unsigned AxisymmetricQTaylorHoodElement::Initial_Nvalue[9] = {
    4, 3, 4, 3, 3, 3, 4, 3, 4};

  // Set the data for the pressure conversion array
  const unsigned AxisymmetricQTaylorHoodElement::Pconv[4] = {0, 2, 6, 8};


  //========================================================================
  /// Compute traction at local coordinate s for outer unit normal N
  //========================================================================
  void AxisymmetricQTaylorHoodElement::get_traction(
    const Vector<double>& s,
    const Vector<double>& N,
    Vector<double>& traction) const
  {
    AxisymmetricNavierStokesEquations::traction(s, N, traction);
  }

} // namespace oomph
