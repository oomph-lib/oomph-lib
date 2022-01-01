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
#include "polar_navier_stokes_elements.h"

namespace oomph
{
  /// ///////////////////////////////////////////////////////////////////////
  //======================================================================//
  /// Start of what would've been navier_stokes_elements.cc               //
  //======================================================================//
  /// ///////////////////////////////////////////////////////////////////////

  /// Navier--Stokes equations static data
  Vector<double> PolarNavierStokesEquations::Gamma(2, 0.0);

  //=================================================================
  /// "Magic" negative number that indicates that the pressure is
  /// not stored at a node. This cannot be -1 because that represents
  /// the positional hanging scheme in the hanging_pt object of nodes
  //=================================================================
  int PolarNavierStokesEquations::Pressure_not_stored_at_node = -100;

  /// Navier--Stokes equations static data
  double PolarNavierStokesEquations::Default_Physical_Constant_Value = 0.0;

  /// Navier--Stokes equations static data
  double PolarNavierStokesEquations::Default_Physical_Ratio_Value = 1.0;

  /// Navier-Stokes equations default gravity vector
  Vector<double> PolarNavierStokesEquations::Default_Gravity_vector(2, 0.0);

  //======================================================================
  /// Validate against exact velocity solution at given time.
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points and compute L2 error
  /// and L2 norm of velocity solution over element.
  //=======================================================================
  void PolarNavierStokesEquations::compute_error(
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
    Vector<double> exact_soln(2 + 1);

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
      for (unsigned i = 0; i < 2; i++)
      {
        norm += exact_soln[i] * exact_soln[i] * W;
        error += (exact_soln[i] - interpolated_u_pnst(s, i)) *
                 (exact_soln[i] - interpolated_u_pnst(s, i)) * W;
      }

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output x,y,[z],u_error,v_error,[w_error]
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << exact_soln[i] - interpolated_u_pnst(s, i) << " ";
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
  void PolarNavierStokesEquations::compute_error(
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

    // Exact solution Vector (u,v,[w],p)
    Vector<double> exact_soln(2 + 1);

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
      (*exact_soln_pt)(x, exact_soln);

      // Velocity error
      for (unsigned i = 0; i < 2; i++)
      {
        norm += exact_soln[i] * exact_soln[i] * W;
        error += (exact_soln[i] - interpolated_u_pnst(s, i)) *
                 (exact_soln[i] - interpolated_u_pnst(s, i)) * W;
      }

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << x[i] << " ";
      }

      // Output x,y,[z],u_error,v_error,[w_error]
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << exact_soln[i] - interpolated_u_pnst(s, i) << " ";
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
  void PolarNavierStokesEquations::output_fct(
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
  void PolarNavierStokesEquations::output_fct(
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
  void PolarNavierStokesEquations::output_veloc(std::ostream& outfile,
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
    Vector<double> interpolated_u(2);

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
      for (unsigned i = 0; i < 2; i++)
      {
        interpolated_x[i] = 0.0;
        interpolated_u[i] = 0.0;
        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          interpolated_u[i] += u_pnst(t, l, i) * psi[l];
          interpolated_x[i] += nodal_position(t, l, i) * psi[l];
        }
      }

      // Coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_x[i] << " ";
      }

      // Velocities
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_u[i] << " ";
      }

      outfile << std::endl;
    }
    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }

  //  keyword: primary
  //==============================================================
  /// Output function:
  /// x,y,[z],u,v,[w],p
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  void PolarNavierStokesEquations::output(std::ostream& outfile,
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

      // Work out global physical coordinate
      const double Alpha = alpha();
      double r = interpolated_x(s, 0);
      double phi = interpolated_x(s, 1);
      double theta = Alpha * phi;

      // Coordinates
      outfile << r * cos(theta) << " " << r * sin(theta) << " ";

      // Velocities
      outfile << interpolated_u_pnst(s, 0) * cos(theta) -
                   interpolated_u_pnst(s, 1) * sin(theta)
              << " ";
      outfile << interpolated_u_pnst(s, 0) * sin(theta) +
                   interpolated_u_pnst(s, 1) * cos(theta)
              << " ";

      // Pressure
      outfile << interpolated_p_pnst(s) << " ";

      // Radial and Azimuthal velocities
      outfile << interpolated_u_pnst(s, 0) << " " << interpolated_u_pnst(s, 1)
              << " ";
      // comment start here
      /*
      double similarity_solution,dsimilarity_solution;
      get_similarity_solution(theta,Alpha,similarity_solution,dsimilarity_solution);
      // Similarity solution:
      outfile << similarity_solution/(Alpha*r) << " ";
      // Error from similarity solution:
      outfile << interpolated_u_pnst(s,0) - similarity_solution/(Alpha*r) << "
      ";
      */
      /*
      //Work out the Stokes flow similarity solution
      double mult = (Alpha/(sin(2.*Alpha)-2.*Alpha*cos(2.*Alpha)));
      double similarity_solution = mult*(cos(2.*Alpha*phi)-cos(2.*Alpha));
      // Similarity solution:
      outfile << similarity_solution/(Alpha*r) << " ";
      // Error from similarity solution:
      outfile << interpolated_u_pnst(s,0) - similarity_solution/(Alpha*r) << "
      ";
      //comment end here
      */
      outfile << 0 << " " << 1 << " ";

      // r and theta for better plotting
      outfile << r << " " << phi << " ";

      outfile << std::endl;
    }
    outfile << std::endl;
  }

  /*
  //Full_convergence_checks output:
  //==============================================================
  /// Output the exact similarity solution and error
  /// Output function in tecplot format.
  /// Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  void PolarNavierStokesEquations::output(std::ostream &outfile,
                                 const unsigned &nplot)
  {

   //Vector of local coordinates
   Vector<double> s(2);

   // Tecplot header info
   outfile << tecplot_zone_string(nplot);

   // Loop over plot points
   unsigned num_plot_points=nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {

     // Get local coordinates of plot point
     get_s_plot(iplot,nplot,s);

     //Work out global physical coordinate
     const double Alpha = alpha();
     double r = interpolated_x(s,0);
     double phi = interpolated_x(s,1);
     double theta = Alpha*phi;

     // Coordinates
     outfile << r*cos(theta) << " " << r*sin(theta) << " ";

     // Velocities
     outfile << interpolated_u_pnst(s,0)*cos(theta) -
  interpolated_u_pnst(s,1)*sin(theta)
       << " ";
     outfile << interpolated_u_pnst(s,0)*sin(theta) +
  interpolated_u_pnst(s,1)*cos(theta)
       << " ";

     // Pressure
     outfile << interpolated_p_pnst(s)  << " ";

     // Radial and Azimuthal velocities
     outfile << interpolated_u_pnst(s,0) << " " << interpolated_u_pnst(s,1) << "
  ";

     // I need to add exact pressure, radial velocity and radial velocity
  derivatives here
     // Plus the error from these
     double m=2.*Alpha;
     double mu=(1./(sin(m)-m*cos(m)));
     double exact_u=(mu/r)*(cos(m*phi)-cos(m));
     double exact_p=(2.*mu)*((cos(m*phi)/(r*r))-cos(m));
     double exact_dudr=-(exact_u/r);
     double exact_dudphi=-mu*m*sin(m*phi)/r;

     // If we don't have Stokes flow then we need to overwrite the Stokes
  solution by
     // reading in the correct similarity solution from file
     const double Re = re();
     if(Re>1.e-8)
      {
       double similarity_solution,dsimilarity_solution;
       get_similarity_solution(theta,Alpha,similarity_solution,dsimilarity_solution);
       double A = Global_Physical_Variables::P2/Alpha;

       exact_u = similarity_solution/(r*Alpha);
       exact_p = 2.*(exact_u/r)-(A/(2.*Alpha*Alpha))*((1./(r*r))-1.);
       exact_dudr = -(exact_u/r);
       exact_dudphi = dsimilarity_solution/(r*Alpha);
      }

     // exact pressure and error
     outfile << exact_p << " " << (interpolated_p_pnst(s)-exact_p) << " ";

     // r and theta for better plotting
     outfile << r << " " << phi << " ";

     // exact and oomph du/dr and du/dphi?
     outfile << exact_u << " " << (interpolated_u_pnst(s,0)-exact_u) << " ";
     outfile << exact_dudr << " " << (interpolated_dudx_pnst(s,0,0)-exact_dudr)
  << " "; outfile << exact_dudphi << " " <<
  (interpolated_dudx_pnst(s,0,1)-exact_dudphi) << " ";

     outfile << std::endl;
    }
   outfile << std::endl;
  }
  */
  //==============================================================
  /// C-style output function:
  /// x,y,[z],u,v,[w],p
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction.
  //==============================================================
  void PolarNavierStokesEquations::output(FILE* file_pt, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Tecplot header info
    // outfile << tecplot_zone_string(nplot);
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
        // outfile << interpolated_x(s,i) << " ";
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }

      // Velocities
      for (unsigned i = 0; i < 2; i++)
      {
        // outfile << interpolated_u_pnst(s,i) << " ";
        fprintf(file_pt, "%g ", interpolated_u_pnst(s, i));
      }

      // Pressure
      // outfile << interpolated_p_pnst(s)  << " ";
      fprintf(file_pt, "%g \n", interpolated_p_pnst(s));
    }
    fprintf(file_pt, "\n");
  }


  //==============================================================
  /// Full output function:
  /// x,y,[z],u,v,[w],p,du/dt,dv/dt,[dw/dt],dissipation
  /// in tecplot format. Specified number of plot points in each
  /// coordinate direction
  //==============================================================
  void PolarNavierStokesEquations::full_output(std::ostream& outfile,
                                               const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Acceleration
    Vector<double> dudt(2);

    // Mesh elocity
    Vector<double> mesh_veloc(2);

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
      Vector<double> mesh_veloc(2);
      Vector<double> dudt(2);
      Vector<double> dudt_ALE(2);
      DenseMatrix<double> interpolated_dudx(2, 2);
      // DenseDoubleMatrix interpolated_dudx(2,2);

      // Initialise everything to zero
      for (unsigned i = 0; i < 2; i++)
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

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned i = 0; i < 2; i++)
        {
          dudt[i] += du_dt_pnst(l, i) * psif[l];
          mesh_veloc[i] += dnodal_position_dt(l, i) * psif[l];

          // Loop over derivative directions for velocity gradients
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_pnst(l, i) * dpsifdx(l, j);
          }
        }
      }


      // Get dudt in ALE form (incl mesh veloc)
      for (unsigned i = 0; i < 2; i++)
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
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_u_pnst(s, i) << " ";
      }

      // Pressure
      outfile << interpolated_p_pnst(s) << " ";

      // Accelerations
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << dudt_ALE[i] << " ";
      }

      // Dissipation
      outfile << dissipation(s) << " ";


      outfile << std::endl;
    }
  }

  //==============================================================
  /// Return integral of dissipation over element
  //==============================================================
  double PolarNavierStokesEquations::dissipation() const
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
      DenseMatrix<double> strainrate(2, 2);
      // DenseDoubleMatrix strainrate(2,2);
      strain_rate(s, strainrate);

      // Initialise
      double local_diss = 0.0;
      for (unsigned i = 0; i < 2; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          local_diss += 2.0 * strainrate(i, j) * strainrate(i, j);
        }
      }

      diss += local_diss * w * J;
    }

    return diss;
  }

  //==============================================================
  /// Compute traction (on the viscous scale) at local
  /// coordinate s for outer unit normal N
  //==============================================================
  void PolarNavierStokesEquations::get_traction(const Vector<double>& s,
                                                const Vector<double>& N,
                                                Vector<double>& traction)
  {
    // Get velocity gradients
    DenseMatrix<double> strainrate(2, 2);
    // DenseDoubleMatrix strainrate(2,2);
    strain_rate(s, strainrate);

    // Get pressure
    double press = interpolated_p_pnst(s);

    // Loop over traction components
    for (unsigned i = 0; i < 2; i++)
    {
      traction[i] = -press * N[i];
      for (unsigned j = 0; j < 2; j++)
      {
        traction[i] += 2.0 * strainrate(i, j) * N[j];
      }
    }
  }

  //==============================================================
  /// Return dissipation at local coordinate s
  //==============================================================
  double PolarNavierStokesEquations::dissipation(const Vector<double>& s) const
  {
    // Get strain rate matrix
    DenseMatrix<double> strainrate(2, 2);
    // DenseDoubleMatrix strainrate(2,2);
    strain_rate(s, strainrate);

    // Initialise
    double local_diss = 0.0;
    for (unsigned i = 0; i < 2; i++)
    {
      for (unsigned j = 0; j < 2; j++)
      {
        local_diss += 2.0 * strainrate(i, j) * strainrate(i, j);
      }
    }

    return local_diss;
  }

  //==============================================================
  /// Get strain-rate tensor:
  /// Slightly more complicated in polar coordinates (see eg. Aris)
  //==============================================================
  void PolarNavierStokesEquations::strain_rate(
    const Vector<double>& s, DenseMatrix<double>& strainrate) const
  {
#ifdef PARANOID
    if ((strainrate.ncol() != 2) || (strainrate.nrow() != 2))
    {
      std::ostringstream error_stream;
      error_stream << "Wrong size " << strainrate.ncol() << " "
                   << strainrate.nrow() << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Velocity gradient matrix
    DenseMatrix<double> dudx(2);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psif(n_node);
    DShape dpsifdx(n_node, 2);

    // Call the derivatives of the shape functions
    dshape_eulerian(s, psif, dpsifdx);

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[2];
    for (unsigned i = 0; i < 2; i++)
    {
      u_nodal_index[i] = u_index_pnst(i);
    }

    // Calculate local values of the velocity components
    // Allocate
    Vector<double> interpolated_u(2);
    Vector<double> interpolated_x(2);
    DenseMatrix<double> interpolated_dudx(2, 2);

    // Initialise to zero
    for (unsigned i = 0; i < 2; i++)
    {
      interpolated_u[i] = 0.0;
      interpolated_x[i] = 0.0;
      for (unsigned j = 0; j < 2; j++)
      {
        interpolated_dudx(i, j) = 0.0;
      }
    }

    // Calculate velocities and derivatives:
    // Loop over nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over directions
      for (unsigned i = 0; i < 2; i++)
      {
        // Get the nodal value
        double u_value = this->nodal_value(l, u_nodal_index[i]);
        interpolated_u[i] += u_value * psif[l];
        interpolated_x[i] += this->nodal_position(l, i) * psif[l];

        // Loop over derivative directions
        for (unsigned j = 0; j < 2; j++)
        {
          interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
        }
      }
    }

    const double Alpha = alpha();
    // Can no longer loop over strain components
    strainrate(0, 0) = interpolated_dudx(0, 0);
    strainrate(1, 1) =
      (1. / (Alpha * interpolated_x[0])) * interpolated_dudx(1, 1) +
      (interpolated_u[0] / interpolated_x[0]);
    strainrate(1, 0) =
      0.5 * (interpolated_dudx(1, 0) - (interpolated_u[1] / interpolated_x[0]) +
             (1. / (Alpha * interpolated_x[0])) * interpolated_dudx(0, 1));
    strainrate(0, 1) = strainrate(1, 0);
  }

  //==============================================================
  /// Return polar strain-rate tensor multiplied by r
  /// Slightly more complicated in polar coordinates (see eg. Aris)
  //==============================================================
  void PolarNavierStokesEquations::strain_rate_by_r(
    const Vector<double>& s, DenseMatrix<double>& strainrate) const
  {
#ifdef PARANOID
    if ((strainrate.ncol() != 2) || (strainrate.nrow() != 2))
    {
      std::ostringstream error_stream;
      error_stream << "Wrong size " << strainrate.ncol() << " "
                   << strainrate.nrow() << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Velocity gradient matrix
    DenseMatrix<double> dudx(2);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psif(n_node);
    DShape dpsifdx(n_node, 2);

    // Call the derivatives of the shape functions
    dshape_eulerian(s, psif, dpsifdx);

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[2];
    for (unsigned i = 0; i < 2; i++)
    {
      u_nodal_index[i] = u_index_pnst(i);
    }

    // Calculate local values of the velocity components
    // Allocate
    Vector<double> interpolated_u(2);
    Vector<double> interpolated_x(2);
    DenseMatrix<double> interpolated_dudx(2, 2);

    // Initialise to zero
    for (unsigned i = 0; i < 2; i++)
    {
      interpolated_u[i] = 0.0;
      interpolated_x[i] = 0.0;
      for (unsigned j = 0; j < 2; j++)
      {
        interpolated_dudx(i, j) = 0.0;
      }
    }

    // Calculate velocities and derivatives:
    // Loop over nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over directions
      for (unsigned i = 0; i < 2; i++)
      {
        // Get the nodal value
        double u_value = this->nodal_value(l, u_nodal_index[i]);
        interpolated_u[i] += u_value * psif[l];
        interpolated_x[i] += this->nodal_position(l, i) * psif[l];

        // Loop over derivative directions
        for (unsigned j = 0; j < 2; j++)
        {
          interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
        }
      }
    }

    const double Alpha = alpha();
    // Can no longer loop over strain components
    strainrate(0, 0) = interpolated_dudx(0, 0);
    strainrate(1, 1) =
      (1. / (Alpha * interpolated_x[0])) * interpolated_dudx(1, 1) +
      (interpolated_u[0] / interpolated_x[0]);
    strainrate(1, 0) =
      0.5 * (interpolated_dudx(1, 0) - (interpolated_u[1] / interpolated_x[0]) +
             (1. / (Alpha * interpolated_x[0])) * interpolated_dudx(0, 1));
    strainrate(0, 1) = strainrate(1, 0);

    strainrate(0, 0) *= interpolated_x[0];
    strainrate(1, 1) *= interpolated_x[0];
    strainrate(1, 0) *= interpolated_x[0];
    strainrate(0, 1) *= interpolated_x[0];
    /*
    strainrate(0,0)*=interpolated_x[0];
    strainrate(1,1)*=interpolated_x[0];
    strainrate(1,0)*=interpolated_x[0];
    strainrate(0,1)*=interpolated_x[0];
    */
  }

  //==============================================================
  ///  Get integral of kinetic energy over element:
  //==============================================================
  double PolarNavierStokesEquations::kin_energy() const
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
      for (unsigned i = 0; i < 2; i++)
      {
        veloc_squared += interpolated_u_pnst(s, i) * interpolated_u_pnst(s, i);
      }

      kin_en += 0.5 * veloc_squared * w * J;
    }

    return kin_en;
  }

  //==============================================================
  /// Return pressure integrated over the element
  //==============================================================
  double PolarNavierStokesEquations::pressure_integral() const
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
      double press = interpolated_p_pnst(s);

      // Add
      press_int += press * W;
    }

    return press_int;
  }

  /// ///////////////////////////////////////////////////////////////////
  /// / The finished version of the new equations                   /////
  /// ///////////////////////////////////////////////////////////////////

  //==============================================================
  ///  Compute the residuals for the Navier--Stokes
  ///  equations; flag=1(or 0): do (or don't) compute the
  ///  Jacobian as well.
  ///  flag=2 for Residuals, Jacobian and mass_matrix
  ///
  ///  This is now my new version with Jacobian and
  ///  dimensionless phi
  //==============================================================
  void PolarNavierStokesEquations::fill_in_generic_residual_contribution(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix,
    unsigned flag)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many pressure dofs there are
    unsigned n_pres = npres_pnst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[2];
    for (unsigned i = 0; i < 2; i++)
    {
      u_nodal_index[i] = u_index_pnst(i);
    }

    // Set up memory for the shape and test functions
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, 2), dtestfdx(n_node, 2);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(2);

    // Get the reynolds number and Alpha
    const double Re = re();
    const double Alpha = alpha();
    const double Re_St = re_st();

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
      double J = dshape_and_dtest_eulerian_at_knot_pnst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Call the pressure shape and test functions
      pshape_pnst(s, psip, testp);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the pressure and velocity components
      // Allocate
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(2);
      Vector<double> interpolated_x(2);
      // Vector<double> mesh_velocity(2);
      // Vector<double> dudt(2);
      DenseMatrix<double> interpolated_dudx(2, 2);

      // Initialise to zero
      for (unsigned i = 0; i < 2; i++)
      {
        // dudt[i] = 0.0;
        // mesh_velocity[i] = 0.0;
        interpolated_u[i] = 0.0;
        interpolated_x[i] = 0.0;
        for (unsigned j = 0; j < 2; j++)
        {
          interpolated_dudx(i, j) = 0.0;
        }
      }

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
        interpolated_p += p_pnst(l) * psip[l];

      // Calculate velocities and derivatives:

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned i = 0; i < 2; i++)
        {
          // Get the nodal value
          double u_value = this->nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif[l];
          interpolated_x[i] += this->nodal_position(l, i) * psif[l];
          // dudt[i]+=du_dt_pnst(l,i)*psif[l];
          // mesh_velocity[i] +=dx_dt(l,i)*psif[l];

          // Loop over derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      // MOMENTUM EQUATIONS
      //------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Can't loop over velocity components as don't have identical
        // contributions Do seperately for i = {0,1} instead
        unsigned i = 0;
        {
          /*IF it's not a boundary condition*/
          local_eqn = nodal_local_eqn(l, u_nodal_index[i]);
          if (local_eqn >= 0)
          {
            // Add the testf[l] term of the stress tensor
            residuals[local_eqn] +=
              ((interpolated_p / interpolated_x[0]) -
               ((1. + Gamma[i]) / pow(interpolated_x[0], 2.)) *
                 ((1. / Alpha) * interpolated_dudx(1, 1) + interpolated_u[0])) *
              testf[l] * interpolated_x[0] * Alpha * W;

            // Add the dtestfdx(l,0) term of the stress tensor
            residuals[local_eqn] +=
              (interpolated_p - (1. + Gamma[i]) * interpolated_dudx(0, 0)) *
              dtestfdx(l, 0) * interpolated_x[0] * Alpha * W;

            // Add the dtestfdx(l,1) term of the stress tensor
            residuals[local_eqn] -=
              ((1. / (interpolated_x[0] * Alpha)) * interpolated_dudx(0, 1) -
               (interpolated_u[1] / interpolated_x[0]) +
               Gamma[i] * interpolated_dudx(1, 0)) *
              (1. / (interpolated_x[0] * Alpha)) * dtestfdx(l, 1) *
              interpolated_x[0] * Alpha * W;

            // Convective terms
            residuals[local_eqn] -=
              Re *
              (interpolated_u[0] * interpolated_dudx(0, 0) +
               (interpolated_u[1] / (interpolated_x[0] * Alpha)) *
                 interpolated_dudx(0, 1) -
               (pow(interpolated_u[1], 2.) / interpolated_x[0])) *
              testf[l] * interpolated_x[0] * Alpha * W;


            // CALCULATE THE JACOBIAN
            if (flag)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Again can't loop over velocity components due to loss of
                // symmetry
                unsigned i2 = 0;
                {
                  // If at a non-zero degree of freedom add in the entry
                  local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
                  if (local_unknown >= 0)
                  {
                    // Add contribution to Elemental Matrix
                    jacobian(local_eqn, local_unknown) -=
                      (1. + Gamma[i]) *
                      (psif[l2] / pow(interpolated_x[0], 2.)) * testf[l] *
                      interpolated_x[0] * Alpha * W;

                    jacobian(local_eqn, local_unknown) -=
                      (1. + Gamma[i]) * dpsifdx(l2, 0) * dtestfdx(l, 0) *
                      interpolated_x[0] * Alpha * W;

                    jacobian(local_eqn, local_unknown) -=
                      (1. / (interpolated_x[0] * Alpha)) * dpsifdx(l2, 1) *
                      (1. / (interpolated_x[0] * Alpha)) * dtestfdx(l, 1) *
                      interpolated_x[0] * Alpha * W;

                    // Now add in the inertial terms
                    jacobian(local_eqn, local_unknown) -=
                      Re *
                      (psif[l2] * interpolated_dudx(0, 0) +
                       interpolated_u[0] * dpsifdx(l2, 0) +
                       (interpolated_u[1] / (interpolated_x[0] * Alpha)) *
                         dpsifdx(l2, 1)) *
                      testf[l] * interpolated_x[0] * Alpha * W;

                    // extra bit for mass matrix
                    if (flag == 2)
                    {
                      mass_matrix(local_eqn, local_unknown) +=
                        Re_St * psif[l2] * testf[l] * interpolated_x[0] *
                        Alpha * W;
                    }

                  } // End of (Jacobian's) if not boundary condition statement
                } // End of i2=0 section

                i2 = 1;
                {
                  // If at a non-zero degree of freedom add in the entry
                  local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
                  if (local_unknown >= 0)
                  {
                    // Add contribution to Elemental Matrix
                    jacobian(local_eqn, local_unknown) -=
                      ((1. + Gamma[i]) / (pow(interpolated_x[0], 2.) * Alpha)) *
                      dpsifdx(l2, 1) * testf[l] * interpolated_x[0] * Alpha * W;

                    jacobian(local_eqn, local_unknown) -=
                      (-(psif[l2] / interpolated_x[0]) +
                       Gamma[i] * dpsifdx(l2, 0)) *
                      (1. / (interpolated_x[0] * Alpha)) * dtestfdx(l, 1) *
                      interpolated_x[0] * Alpha * W;

                    // Now add in the inertial terms
                    jacobian(local_eqn, local_unknown) -=
                      Re *
                      ((psif[l2] / (interpolated_x[0] * Alpha)) *
                         interpolated_dudx(0, 1) -
                       2 * interpolated_u[1] * (psif[l2] / interpolated_x[0])) *
                      testf[l] * interpolated_x[0] * Alpha * W;

                  } // End of (Jacobian's) if not boundary condition statement
                } // End of i2=1 section

              } // End of l2 loop

              /*Now loop over pressure shape functions*/
              /*This is the contribution from pressure gradient*/
              for (unsigned l2 = 0; l2 < n_pres; l2++)
              {
                /*If we are at a non-zero degree of freedom in the entry*/
                local_unknown = p_local_eqn(l2);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    (psip[l2] / interpolated_x[0]) * testf[l] *
                    interpolated_x[0] * Alpha * W;

                  jacobian(local_eqn, local_unknown) +=
                    psip[l2] * dtestfdx(l, 0) * interpolated_x[0] * Alpha * W;
                }
              }
            } /*End of Jacobian calculation*/

          } // End of if not boundary condition statement
        } // End of i=0 section

        i = 1;
        {
          /*IF it's not a boundary condition*/
          local_eqn = nodal_local_eqn(l, u_nodal_index[i]);
          if (local_eqn >= 0)
          {
            // Add the testf[l] term of the stress tensor
            residuals[local_eqn] +=
              ((1. / (pow(interpolated_x[0], 2.) * Alpha)) *
                 interpolated_dudx(0, 1) -
               (interpolated_u[1] / pow(interpolated_x[0], 2.)) +
               Gamma[i] * (1. / interpolated_x[0]) * interpolated_dudx(1, 0)) *
              testf[l] * interpolated_x[0] * Alpha * W;

            // Add the dtestfdx(l,0) term of the stress tensor
            residuals[local_eqn] -=
              (interpolated_dudx(1, 0) +
               Gamma[i] *
                 ((1. / (interpolated_x[0] * Alpha)) * interpolated_dudx(0, 1) -
                  (interpolated_u[1] / interpolated_x[0]))) *
              dtestfdx(l, 0) * interpolated_x[0] * Alpha * W;

            // Add the dtestfdx(l,1) term of the stress tensor
            residuals[local_eqn] +=
              (interpolated_p -
               (1. + Gamma[i]) *
                 ((1. / (interpolated_x[0] * Alpha)) * interpolated_dudx(1, 1) +
                  (interpolated_u[0] / interpolated_x[0]))) *
              (1. / (interpolated_x[0] * Alpha)) * dtestfdx(l, 1) *
              interpolated_x[0] * Alpha * W;

            // Convective terms
            residuals[local_eqn] -=
              Re *
              (interpolated_u[0] * interpolated_dudx(1, 0) +
               (interpolated_u[1] / (interpolated_x[0] * Alpha)) *
                 interpolated_dudx(1, 1) +
               ((interpolated_u[0] * interpolated_u[1]) / interpolated_x[0])) *
              testf[l] * interpolated_x[0] * Alpha * W;

            // CALCULATE THE JACOBIAN
            if (flag)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Again can't loop over velocity components due to loss of
                // symmetry
                unsigned i2 = 0;
                {
                  // If at a non-zero degree of freedom add in the entry
                  local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
                  if (local_unknown >= 0)
                  {
                    // Add contribution to Elemental Matrix
                    jacobian(local_eqn, local_unknown) +=
                      (1. / (pow(interpolated_x[0], 2.) * Alpha)) *
                      dpsifdx(l2, 1) * testf[l] * interpolated_x[0] * Alpha * W;

                    jacobian(local_eqn, local_unknown) -=
                      Gamma[i] * (1. / (interpolated_x[0] * Alpha)) *
                      dpsifdx(l2, 1) * dtestfdx(l, 0) * interpolated_x[0] *
                      Alpha * W;

                    jacobian(local_eqn, local_unknown) -=
                      (1 + Gamma[i]) * (psif[l2] / interpolated_x[0]) *
                      (1. / (interpolated_x[0] * Alpha)) * dtestfdx(l, 1) *
                      interpolated_x[0] * Alpha * W;

                    // Now add in the inertial terms
                    jacobian(local_eqn, local_unknown) -=
                      Re *
                      (psif[l2] * interpolated_dudx(1, 0) +
                       (psif[l2] * interpolated_u[1] / interpolated_x[0])) *
                      testf[l] * interpolated_x[0] * Alpha * W;

                  } // End of (Jacobian's) if not boundary condition statement
                } // End of i2=0 section

                i2 = 1;
                {
                  // If at a non-zero degree of freedom add in the entry
                  local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
                  if (local_unknown >= 0)
                  {
                    // Add contribution to Elemental Matrix
                    jacobian(local_eqn, local_unknown) +=
                      (-(psif[l2] / pow(interpolated_x[0], 2.)) +
                       Gamma[i] * (1. / interpolated_x[0]) * dpsifdx(l2, 0)) *
                      testf[l] * interpolated_x[0] * Alpha * W;

                    jacobian(local_eqn, local_unknown) -=
                      (dpsifdx(l2, 0) -
                       Gamma[i] * (psif[l2] / interpolated_x[0])) *
                      dtestfdx(l, 0) * interpolated_x[0] * Alpha * W;

                    jacobian(local_eqn, local_unknown) -=
                      (1. + Gamma[i]) * (1. / (interpolated_x[0] * Alpha)) *
                      dpsifdx(l2, 1) * (1. / (interpolated_x[0] * Alpha)) *
                      dtestfdx(l, 1) * interpolated_x[0] * Alpha * W;

                    // Now add in the inertial terms
                    jacobian(local_eqn, local_unknown) -=
                      Re *
                      (interpolated_u[0] * dpsifdx(l2, 0) +
                       (psif[l2] / (interpolated_x[0] * Alpha)) *
                         interpolated_dudx(1, 1) +
                       (interpolated_u[1] / (interpolated_x[0] * Alpha)) *
                         dpsifdx(l2, 1) +
                       (interpolated_u[0] * psif[l2] / interpolated_x[0])) *
                      testf[l] * interpolated_x[0] * Alpha * W;

                    // extra bit for mass matrix
                    if (flag == 2)
                    {
                      mass_matrix(local_eqn, local_unknown) +=
                        Re_St * psif[l2] * testf[l] * interpolated_x[0] *
                        Alpha * W;
                    }

                  } // End of (Jacobian's) if not boundary condition statement
                } // End of i2=1 section

              } // End of l2 loop

              /*Now loop over pressure shape functions*/
              /*This is the contribution from pressure gradient*/
              for (unsigned l2 = 0; l2 < n_pres; l2++)
              {
                /*If we are at a non-zero degree of freedom in the entry*/
                local_unknown = p_local_eqn(l2);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    (psip[l2] / interpolated_x[0]) * (1. / Alpha) *
                    dtestfdx(l, 1) * interpolated_x[0] * Alpha * W;
                }
              }
            } /*End of Jacobian calculation*/

          } // End of if not boundary condition statement

        } // End of i=1 section

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
          residuals[local_eqn] +=
            (interpolated_dudx(0, 0) + (interpolated_u[0] / interpolated_x[0]) +
             (1. / (interpolated_x[0] * Alpha)) * interpolated_dudx(1, 1)) *
            testp[l] * interpolated_x[0] * Alpha * W;


          /*CALCULATE THE JACOBIAN*/
          if (flag)
          {
            /*Loop over the velocity shape functions*/
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              unsigned i2 = 0;
              {
                /*If we're at a non-zero degree of freedom add it in*/
                local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    (dpsifdx(l2, 0) + (psif[l2] / interpolated_x[0])) *
                    testp[l] * interpolated_x[0] * Alpha * W;
                }
              } // End of i2=0 section

              i2 = 1;
              {
                /*If we're at a non-zero degree of freedom add it in*/
                local_unknown = nodal_local_eqn(l2, u_nodal_index[i2]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    (1. / (interpolated_x[0] * Alpha)) * dpsifdx(l2, 1) *
                    testp[l] * interpolated_x[0] * Alpha * W;
                }
              } // End of i2=1 section

            } /*End of loop over l2*/
          } /*End of Jacobian calculation*/

        } // End of if not boundary condition
      } // End of loop over l


    } // End of loop over integration points

  } // End of add_generic_residual_contribution

  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////

  /// 2D Crouzeix-Raviart elements
  // Set the data for the number of Variables at each node
  const unsigned PolarCrouzeixRaviartElement::Initial_Nvalue[9] = {
    2, 2, 2, 2, 2, 2, 2, 2, 2};

  //========================================================================
  /// Number of values (pinned or dofs) required at node n.
  //========================================================================
  unsigned PolarCrouzeixRaviartElement::required_nvalue(
    const unsigned int& n) const
  {
    return Initial_Nvalue[n];
  }

  //=========================================================================
  /// Add pointers to Data and indices of the values
  /// that affect the potential load (traction) applied
  /// to the SolidElements to the set paired_load_data
  //=========================================================================
  void PolarCrouzeixRaviartElement::get_load_data(
    std::set<std::pair<Data*, unsigned>>& paired_load_data)
  {
    // Loop over the nodes
    unsigned n_node = this->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      // Loop over the velocity components and add pointer to their data
      // and indices to the vectors
      for (unsigned i = 0; i < 2; i++)
      {
        paired_load_data.insert(std::make_pair(this->node_pt(n), i));
      }
    }

    // Loop over the internal data
    unsigned n_internal = this->ninternal_data();
    for (unsigned l = 0; l < n_internal; l++)
    {
      unsigned nval = this->internal_data_pt(l)->nvalue();
      // Add internal data
      for (unsigned j = 0; j < nval; j++)
      {
        paired_load_data.insert(std::make_pair(this->internal_data_pt(l), j));
      }
    }
  }

  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////

  // 2D Taylor--Hood
  // Set the data for the number of Variables at each node
  const unsigned PolarTaylorHoodElement::Initial_Nvalue[9] = {
    3, 2, 3, 2, 2, 2, 3, 2, 3};

  // Set the data for the pressure conversion array
  const unsigned PolarTaylorHoodElement::Pconv[4] = {0, 2, 6, 8};

  //=========================================================================
  /// Add pointers to Data and the indices of the values
  /// that affect the potential load (traction) applied
  /// to the SolidElements to the set  paired_load_data.
  //=========================================================================
  void PolarTaylorHoodElement::get_load_data(
    std::set<std::pair<Data*, unsigned>>& paired_load_data)
  {
    // Loop over the nodes
    unsigned n_node = this->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      // Loop over the velocity components and add pointer to their data
      // and indices to the vectors
      for (unsigned i = 0; i < 2; i++)
      {
        paired_load_data.insert(std::make_pair(this->node_pt(n), i));
      }
    }

    // Loop over the pressure data
    unsigned n_pres = npres_pnst();
    for (unsigned l = 0; l < n_pres; l++)
    {
      // The 2th entry in each nodal data is the pressure, which
      // affects the traction
      paired_load_data.insert(std::make_pair(this->node_pt(Pconv[l]), 2));
    }
  }

  //========================================================================
  /// Compute traction at local coordinate s for outer unit normal N
  //========================================================================
  // template<unsigned 2>
  // void PolarTaylorHoodElement::get_traction(const Vector<double>& s,
  //                                      const Vector<double>& N,
  //                                      Vector<double>& traction)
  //{
  // PolarNavierStokesEquations::traction(s,N,traction);
  //}

} // namespace oomph
