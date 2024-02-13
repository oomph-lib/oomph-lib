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
// Non-inline functions for SpaceTimeUnsteadyHeat elements
#include "discontinuous_galerkin_space_time_unsteady_heat_elements.h"

/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////

namespace oomph
{
  //======================================================================
  // Default parameters
  //======================================================================
  /// Default value for Alpha parameter (thermal inertia)
  template<unsigned SPATIAL_DIM>
  double SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::Default_alpha_parameter =
    1.0;

  /// Default value for Beta parameter (thermal conductivity)
  template<unsigned SPATIAL_DIM>
  double SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::Default_beta_parameter =
    1.0;

  //======================================================================
  // Set the data for the number of variables at each node
  //======================================================================
  template<unsigned SPATIAL_DIM, unsigned NNODE_1D>
  const unsigned
    QUnsteadyHeatSpaceTimeElement<SPATIAL_DIM, NNODE_1D>::Initial_Nvalue = 1;

  //======================================================================
  /// Compute element residual vector and/or element Jacobian matrix
  ///
  /// flag=0: compute only residual vector
  /// flag=1: compute both
  ///
  /// Pure version without hanging nodes
  //======================================================================
  template<unsigned SPATIAL_DIM>
  void SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::
    fill_in_generic_residual_contribution_ust_heat(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find the index at which the variable is stored
    unsigned u_nodal_index = u_index_ust_heat();

    // Set up memory for the shape functions
    Shape psi(n_node);

    // Set up memory for the test functions
    Shape test(n_node);

    // Allocate space for the derivatives of the shape functions
    DShape dpsidx(n_node, SPATIAL_DIM + 1);

    // Allocate space for the derivatives of the test functions
    DShape dtestdx(n_node, SPATIAL_DIM + 1);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Storage for the local coordinates
    Vector<double> s(SPATIAL_DIM + 1, 0.0);

    // Get the Alpha parameter
    double alpha_local = alpha();

    // Get the Beta parameter
    double beta_local = beta();

    // Integer to hold the local equation
    int local_eqn = 0;

    // Integer to hold the local unknowns
    int local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < SPATIAL_DIM + 1; i++)
      {
        // Calculate the i-th local coordinate
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_ust_heat(
        ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Storage for the interpolated time value
      double interpolated_t = 0.0;

      // Storage for the interpolated solution value
      double interpolated_u = 0.0;

      // Storage for the interpolated time-derivative of the solution
      double interpolated_dudt = 0.0;

      // Storage for the spatial coordinates
      Vector<double> interpolated_x(SPATIAL_DIM, 0.0);

      // Storage for the spatial derivatives of the solution
      Vector<double> interpolated_dudx(SPATIAL_DIM, 0.0);

      // Storage for the mesh velocity
      Vector<double> mesh_velocity(SPATIAL_DIM, 0.0);

      //-------------------------------------------------
      // Calculate derivatives and source function value:
      //-------------------------------------------------
      // Loop over the nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the nodal value at the l-th node
        double u_value = raw_nodal_value(l, u_nodal_index);

        // Update the interpolated time value
        interpolated_t += raw_nodal_position(l, SPATIAL_DIM) * psi(l);

        // Loop over the coordinate directions (both spatial AND time)
        for (unsigned j = 0; j < SPATIAL_DIM; j++)
        {
          // Update the interpolated x value
          interpolated_x[j] += raw_nodal_position(l, j) * psi(l);

          // Update the interpolated du/dx_j value
          interpolated_dudx[j] += u_value * dpsidx(l, j);
        }

        // Update the interpolated u value
        interpolated_u += u_value * psi(l);

        // Update the interpolated du/dt value
        interpolated_dudt += u_value * dpsidx(l, SPATIAL_DIM);
      } // for (unsigned l=0;l<n_node;l++)

      // Initialise the source term value
      double source = 0.0;

      // Get the interpolated source term value
      get_source_ust_heat(interpolated_t, ipt, interpolated_x, source);

      //---------------------------------
      // Assemble residuals and Jacobian:
      //---------------------------------
      // Loop over the nodes (or equivalently the test functions)
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the local equation number
        local_eqn = nodal_local_eqn(l, u_nodal_index);

        // If it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Add source term and time derivative
          residuals[local_eqn] +=
            (source + alpha_local * interpolated_dudt) * test(l) * W;

          // Loop over the coordinate directions
          for (unsigned k = 0; k < SPATIAL_DIM; k++)
          {
            // Add in the contribution from the Laplace operator
            residuals[local_eqn] +=
              beta_local * interpolated_dudx[k] * dtestdx(l, k) * W;
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
              // Get the local equation number
              local_unknown = nodal_local_eqn(l2, u_nodal_index);

              // If we're at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                // Add in the time derivative contribution
                jacobian(local_eqn, local_unknown) +=
                  (alpha_local * test(l) * dpsidx(l2, SPATIAL_DIM) * W);

                // Laplace operator
                for (unsigned i = 0; i < SPATIAL_DIM; i++)
                {
                  // Add the test function contribution to the Jacobian
                  jacobian(local_eqn, local_unknown) +=
                    (beta_local * dpsidx(l2, i) * dtestdx(l, i) * W);
                }
              } // if (local_unknown>=0)
            } // for (unsigned l2=0;l2<n_node;l2++)
          } // if (flag)
        } // if (local_eqn>=0)
      } // for (unsigned l=0;l<n_node;l++)
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of fill_in_generic_residual_contribution_ust_heat


  //======================================================================
  /// Compute norm of FE solution
  //======================================================================
  template<unsigned SPATIAL_DIM>
  void SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::compute_norm(double& norm)
  {
    // Initialise
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(SPATIAL_DIM + 1, 0.0);

    // Vector for coordinates
    Vector<double> x(SPATIAL_DIM + 1, 0.0);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Allocate memory for the shape and test functions
    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < SPATIAL_DIM + 1; i++)
      {
        // Get the i-th local coordinate at the ipt-th integration point
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get Jacobian of mapping
      double J = J_eulerian(s);

      // Pre-multiply the weights and the Jacobian
      double W = w * J;

      // Get FE function value
      double u = interpolated_u_ust_heat(s);

      // Update the norm value
      norm += u * u * W;
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of compute_norm


  //======================================================================
  /// Self-test: Return 0 for OK
  //======================================================================
  template<unsigned SPATIAL_DIM>
  unsigned SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::self_test()
  {
    // Initialise the boolean variable
    bool passed = true;

    // Check lower-level stuff
    if (FiniteElement::self_test() != 0)
    {
      // If we get here then the lower-level self-tests did not pass
      passed = false;
    }

    // If the self-tests passed
    if (passed)
    {
      // Return the value zero
      return 0;
    }
    // If the self-tests didn't pass
    else
    {
      // Return the value one
      return 1;
    }
  } // End of self_test


  //======================================================================
  /// Output function:
  ///   x,t,u   or    x,y,t,u
  /// at nplot points in each coordinate direction
  //======================================================================
  template<unsigned SPATIAL_DIM>
  void SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::output(
    std::ostream& outfile, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(SPATIAL_DIM + 1, 0.0);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Get the number of plot points
    unsigned num_plot_points = nplot_points(nplot);

    // Loop over plot points
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Loop over the coordinate directions
      for (unsigned i = 0; i < SPATIAL_DIM + 1; i++)
      {
        // Output the interpolated coordinate
        outfile << interpolated_x(s, i) << " ";
      }

      // Calculate the interpolated solution value
      outfile << interpolated_u_ust_heat(s) << std::endl;
    } // for (unsigned iplot=0;iplot<num_plot_points;iplot++)

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  } // End of output


  //======================================================================
  /// C-style output function:
  ///   x,t,u   or    x,y,t,u
  /// at nplot points in each coordinate direction
  //======================================================================
  template<unsigned SPATIAL_DIM>
  void SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::output(
    FILE* file_pt, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(SPATIAL_DIM + 1, 0.0);

    // Tecplot header info
    fprintf(file_pt, "%s", tecplot_zone_string(nplot).c_str());

    // Get the number of plot points
    unsigned num_plot_points = nplot_points(nplot);

    // Loop over plot points
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Loop over the coordinate directions
      for (unsigned i = 0; i < SPATIAL_DIM + 1; i++)
      {
        // Print the i-th coordinate value at local coordinate s
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }

      // Output the interpolated solution value at local coordinate s
      fprintf(file_pt, "%g \n", interpolated_u_ust_heat(s));
    } // for (unsigned iplot=0;iplot<num_plot_points;iplot++)

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, nplot);
  } // End of output


  //======================================================================
  /// Output exact solution at a given number of plot points:
  ///   x,t,u_exact    or    x,y,t,u_exact
  /// Solution is provided via function pointer.
  //======================================================================
  template<unsigned SPATIAL_DIM>
  void SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(SPATIAL_DIM + 1, 0.0);

    // Vector for spatial coordinates
    Vector<double> spatial_coordinates(SPATIAL_DIM, 0.0);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Exact solution vector (here it's simply a scalar)
    Vector<double> exact_soln(1, 0.0);

    // Get the number of plot points
    unsigned num_plot_points = nplot_points(nplot);

    // Loop over plot points
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Loop over the spatial coordinates
      for (unsigned i = 0; i < SPATIAL_DIM; i++)
      {
        // Assign the i-th spatial coordinate
        spatial_coordinates[i] = interpolated_x(s, i);

        // Output the i-th coordinate at the point
        outfile << spatial_coordinates[i] << " ";
      }

      // Output the time value at this point
      outfile << interpolated_x(s, SPATIAL_DIM) << " ";

      // Get the exact solution at this point
      (*exact_soln_pt)(spatial_coordinates, exact_soln);

      // Output the exact solution at this point
      outfile << exact_soln[0] << std::endl;
    } // for (unsigned iplot=0;iplot<num_plot_points;iplot++)

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  } // End of output_fct


  //======================================================================
  /// Output exact solution at a given number of plot points:
  ///   x,t,u_exact    or    x,y,t,u_exact
  /// Solution is provided via function pointer.
  //======================================================================
  template<unsigned SPATIAL_DIM>
  void SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
    // Storage for the time value
    double interpolated_t = 0.0;

    // Vector of local coordinates
    Vector<double> s(SPATIAL_DIM + 1, 0.0);

    // Vector for spatial coordinates
    Vector<double> spatial_coordinates(SPATIAL_DIM, 0.0);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Exact solution vector (here it's simply a scalar)
    Vector<double> exact_soln(1, 0.0);

    // Get the number of plot points
    unsigned num_plot_points = nplot_points(nplot);

    // Loop over plot points
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Loop over the spatial coordinates
      for (unsigned i = 0; i < SPATIAL_DIM; i++)
      {
        // Assign the i-th spatial coordinate
        spatial_coordinates[i] = interpolated_x(s, i);

        // Output the i-th coordinate at the point
        outfile << spatial_coordinates[i] << " ";
      }

      // Get the time value
      interpolated_t = interpolated_x(s, SPATIAL_DIM);

      // Output the time value at this point
      outfile << interpolated_t << " ";

      // Get the exact solution at this point
      (*exact_soln_pt)(interpolated_t, spatial_coordinates, exact_soln);

      // Output the exact solution at this point
      outfile << exact_soln[0] << std::endl;
    } // for (unsigned iplot=0;iplot<num_plot_points;iplot++)

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  } // End of output_fct


  //======================================================================
  /// Validate against exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  //======================================================================
  template<unsigned SPATIAL_DIM>
  void SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    // Initialise error value
    error = 0.0;

    // Initialise norm value
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(SPATIAL_DIM + 1, 0.0);

    // Vector for spatial coordinates
    Vector<double> spatial_coordinates(SPATIAL_DIM, 0.0);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Initialise shape functions
    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Tecplot header info
    outfile << "ZONE" << std::endl;

    // Exact solution vector (here it's simply a scalar)
    Vector<double> exact_soln(1, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < SPATIAL_DIM + 1; i++)
      {
        // Get the i-th local coordinate at the ipt-th integration point
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get FE function value
      double u_fe = interpolated_u_ust_heat(s);

      // Loop over the spatial coordinates
      for (unsigned i = 0; i < SPATIAL_DIM; i++)
      {
        // Assign the i-th spatial coordinate
        spatial_coordinates[i] = interpolated_x(s, i);

        // Output the i-th coordinate at the point
        outfile << spatial_coordinates[i] << " ";
      }

      // Output the i-th coordinate at this point
      outfile << interpolated_x(s, SPATIAL_DIM) << " ";

      // Get exact solution at this point
      (*exact_soln_pt)(spatial_coordinates, exact_soln);

      // Output the error
      outfile << exact_soln[0] << " " << exact_soln[0] - u_fe << std::endl;

      // Add to (exact) solution norm value
      norm += exact_soln[0] * exact_soln[0] * W;

      // Update the error norm value
      error += (exact_soln[0] - u_fe) * (exact_soln[0] - u_fe) * W;
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of compute_error


  //======================================================================
  /// Validate against exact solution at time t.
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  //======================================================================
  template<unsigned SPATIAL_DIM>
  void SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
    double& error,
    double& norm)
  {
    // Initialise error value
    error = 0.0;

    // Initialise norm value
    norm = 0.0;

    // Storage for the time value
    double interpolated_t = 0.0;

    // Vector of local coordinates
    Vector<double> s(SPATIAL_DIM + 1, 0.0);

    // Vector for spatial coordinates
    Vector<double> spatial_coordinates(SPATIAL_DIM, 0.0);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    // Initialise shape functions
    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Tecplot header info
    outfile << "ZONE" << std::endl;

    // Exact solution vector (here it's simply a scalar)
    Vector<double> exact_soln(1, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < SPATIAL_DIM + 1; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get FE function value
      double u_fe = interpolated_u_ust_heat(s);

      // Loop over the spatial coordinates
      for (unsigned i = 0; i < SPATIAL_DIM; i++)
      {
        // Assign the i-th spatial coordinate
        spatial_coordinates[i] = interpolated_x(s, i);

        // Output the i-th coordinate at the point
        outfile << spatial_coordinates[i] << " ";
      }

      // Get the time value
      interpolated_t = interpolated_x(s, SPATIAL_DIM);

      // Output the time value at this point
      outfile << interpolated_t << " ";

      // Get the exact solution at this point
      (*exact_soln_pt)(interpolated_t, spatial_coordinates, exact_soln);

      // Output the error
      outfile << exact_soln[0] << " " << exact_soln[0] - u_fe << std::endl;

      // Add to (exact) solution norm value
      norm += exact_soln[0] * exact_soln[0] * W;

      // Update the error norm value
      error += (exact_soln[0] - u_fe) * (exact_soln[0] - u_fe) * W;
    } // for (unsigned ipt=0;ipt<n_intpt;ipt++)
  } // End of compute_error


  //======================================================================
  /// Output function:
  ///   x,t,u   or    x,y,t,u
  /// at nplot points in each coordinate direction
  //======================================================================
  template<unsigned SPATIAL_DIM>
  void SpaceTimeUnsteadyHeatEquations<SPATIAL_DIM>::output_element_paraview(
    std::ofstream& file_out, const unsigned& nplot)
  {
    // Change the scientific format so that E is used rather than e
    file_out.setf(std::ios_base::uppercase);

    // Make variables to hold the number of nodes and elements
    unsigned number_of_nodes = this->nplot_points_paraview(nplot);

    // Make variables to hold the number of elements
    unsigned total_number_of_elements = this->nsub_elements_paraview(nplot);

    //------------------
    // File Declaration:
    //------------------
    // Insert the necessary lines plus header of file, and
    // number of nodes and elements
    file_out << "<?xml version=\"1.0\"?>\n"
             << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
             << "byte_order=\"LittleEndian\">\n"
             << "<UnstructuredGrid>\n"
             << "<Piece NumberOfPoints=\"" << number_of_nodes
             << "\" NumberOfCells=\"" << total_number_of_elements << "\">\n";

    //------------
    // Point Data:
    //------------
    // Check the number of degrees of freedom
    unsigned ndof = this->nscalar_paraview();

    // Point data is going in here
    file_out << "<PointData ";

    // Insert just the first scalar name, since paraview reads everything
    // else after that as being of the same type. Get information from
    // first element.
    file_out << "Scalars=\"" << this->scalar_name_paraview(0) << "\">\n";

    // Loop over i scalar fields and j number of elements
    for (unsigned i = 0; i < ndof; i++)
    {
      file_out << "<DataArray type=\"Float32\" "
               << "Name=\"" << this->scalar_name_paraview(i) << "\" "
               << "format=\"ascii\""
               << ">\n";

      // Output the i-th scalar field with nplot plot points
      this->scalar_value_paraview(file_out, i, nplot);

      // Close of the DataArray
      file_out << "</DataArray>\n";
    }

    // Close off the PointData set
    file_out << "</PointData>\n";

    //------------------
    // Geometric Points:
    //------------------
    // Always has to be 3 components for an unstructured grid
    file_out << "<Points>\n"
             << "<DataArray type=\"Float32\""
             << " NumberOfComponents=\"" << 3 << "\" "
             << "format=\"ascii\">\n";

    // Print the plot points
    this->output_paraview(file_out, nplot);

    // Close off the geometric points set
    file_out << "</DataArray>\n"
             << "</Points>\n";

    //-------
    // Cells:
    //-------
    file_out << "<Cells>\n"
             << "<DataArray type=\"Int32\" Name=\""
             << "connectivity\" format=\"ascii\">\n";

    // Make counter for keeping track of all the local elements,
    // because Paraview requires global coordinates
    unsigned counter = 0;

    // Write connectivity with the local elements
    this->write_paraview_output_offset_information(file_out, nplot, counter);

    // Output header stuff
    file_out << "</DataArray>\n"
             << "<DataArray type=\"Int32\" "
             << "Name=\"offsets\" format=\"ascii\">\n";

    // Make variable that holds the current offset number
    unsigned offset_sum = 0;

    // Write the offset for the specific elements
    this->write_paraview_offsets(file_out, nplot, offset_sum);

    // Add in header information
    file_out << "</DataArray>\n"
             << "<DataArray type=\"UInt8\" Name=\"types\">\n";

    // Get the type the element has
    this->write_paraview_type(file_out, nplot);

    // Finish off the data set
    file_out << "</DataArray>\n"
             << "</Cells>\n";

    //--------------
    // File Closure:
    //--------------
    file_out << "</Piece>\n"
             << "</UnstructuredGrid>\n"
             << "</VTKFile>";
  } // End of output_element_paraview


  //====================================================================
  // Force build of templates
  //====================================================================
  template class QUnsteadyHeatSpaceTimeElement<2, 2>;
  template class QUnsteadyHeatSpaceTimeElement<2, 3>;
  template class QUnsteadyHeatSpaceTimeElement<2, 4>;
} // End of namespace oomph
