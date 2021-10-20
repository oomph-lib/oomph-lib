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
// Non-inline member functions for discontinuous galerkin elements

// oomph-lib includes
#include "dg_elements.h"
#include "shape.h"
#include <iomanip>

namespace oomph
{
  //===================================================================
  /// Find pointers to neighbouring faces and the local coordinates
  /// in those faces that correspond to the integration points in the
  /// present face.
  /// This is achieved by moving up to the bulk element and thence
  /// the mesh which MUST have implemented a neighbour finding scheme
  //==================================================================
  void DGFaceElement::setup_neighbour_info(
    const bool& add_neighbour_data_to_bulk)
  {
    // Cache the pointer to the bulk element
    DGElement* const bulk_element_pt =
      dynamic_cast<DGElement*>(this->bulk_element_pt());

    // Find the number of points in the integration scheme
    const unsigned n_intpt = integral_pt()->nweight();
    // Resize the storage in the element
    Neighbour_face_pt.resize(n_intpt);
    Neighbour_local_coordinate.resize(n_intpt);

    // If we are adding the neighbour data to the bulk element
    // then resize this storage
    if (add_neighbour_data_to_bulk)
    {
      Neighbour_external_data.resize(n_intpt);
    }


    // Get the dimension of the present element
    const unsigned el_dim = this->dim();
    // Local coordinate in the face element
    Vector<double> s(el_dim);

    // Get the dimension of the bulk element
    const unsigned n_dim = bulk_element_pt->dim();
    // Local coordinate in the bulk element
    Vector<double> s_bulk(n_dim);

    // Storage for the interpolation data in the neighbour
    Vector<Data*> neighbour_data;

    // Now loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the local coordinate of the integration points
      for (unsigned i = 0; i < el_dim; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Now get the bulk coordinate
      this->get_local_coordinate_in_bulk(s, s_bulk);

      // Now we find the neighbouring face via the bulk element's member
      // function which calls the Mesh's member function
      bulk_element_pt->get_neighbouring_face_and_local_coordinate(
        this->face_index(),
        s_bulk,
        Neighbour_face_pt[ipt],
        Neighbour_local_coordinate[ipt]);

      // If we are adding the external data to the bulk
      if (add_neighbour_data_to_bulk)
      {
        // Throw the appropriate data from the neighbour face to the set
        dynamic_cast<DGFaceElement*>(Neighbour_face_pt[ipt])
          ->get_interpolation_data(neighbour_data);

        // Find the number of data
        unsigned n_neighbour_data = neighbour_data.size();
        // Resize the storage accordingly
        Neighbour_external_data.resize(n_neighbour_data);

        // Add the data to the external data of the bulk element
        for (unsigned n = 0; n < n_neighbour_data; n++)
        {
          Neighbour_external_data[ipt][n] =
            bulk_element_pt->add_external_data(neighbour_data[n]);
        }
      }
    }
  }


  //======================================================================
  // Report the global coordinates corresponding to the integration points
  // and the corresponding coordinates in the neighbouring faces
  //======================================================================
  void DGFaceElement::report_info()
  {
    // Find the number of nodes
    const unsigned n_node = this->nnode();
    // Allocate storage for the shape functions
    Shape psi(n_node);

    // Find the dimension of the problem
    const unsigned dim = this->nodal_dimension();
    // Storage for local coordinates in this face and its neighbour
    Vector<double> x(dim), face_x(dim);

    // Find the dimension of the element
    const unsigned el_dim = this->dim();
    // Storage for local coordinate
    Vector<double> s(el_dim);

    // Calculate the number of integration points from the array
    const unsigned n_intpt = this->integral_pt()->nweight();
    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Find the local coordinate at the knot point
      for (unsigned i = 0; i < el_dim; i++)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }
      // Get the value of the global coordinate in the present face
      this->interpolated_x(s, x);

      // Get the value of x in the neighbouring face
      Neighbour_face_pt[ipt]->interpolated_x(Neighbour_local_coordinate[ipt],
                                             face_x);

      // Let's compare
      oomph_info << "In Face                   In Neighbour\n";
      for (unsigned i = 0; i < dim; i++)
      {
        if (i == 0)
        {
          oomph_info << "(";
        }
        else
        {
          oomph_info << ", ";
        }
        oomph_info << std::setw(5) << std::left << x[i];
      }
      oomph_info << ")";

      oomph_info << "                   ";

      for (unsigned i = 0; i < dim; i++)
      {
        if (i == 0)
        {
          oomph_info << "(";
        }
        else
        {
          oomph_info << ", ";
        }
        oomph_info << std::setw(5) << std::left << face_x[i];
      }
      oomph_info << ")";
      oomph_info << std::endl;
    }
  }


  //=====================================================================
  /// Return the interpolated values of the unknown fluxes
  //=====================================================================
  void DGFaceElement::interpolated_u(const Vector<double>& s, Vector<double>& u)
  {
    // Find the number of nodes
    const unsigned n_node = nnode();
    // If there are no nodes then return immediately
    if (n_node == 0)
    {
      return;
    }

    // Get the shape functions at the local coordinate
    Shape psi(n_node);
    this->shape(s, psi);

    // Find the number of fluxes
    const unsigned n_flux = this->required_nflux();

    // Find the indices at which the local fluxes are stored
    Vector<unsigned> flux_nodal_index(n_flux);
    for (unsigned i = 0; i < n_flux; i++)
    {
      flux_nodal_index[i] = this->flux_index(i);
    }
    // Initialise the fluxes to zero
    for (unsigned i = 0; i < n_flux; i++)
    {
      u[i] = 0.0;
    }

    // Loop over the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      const double psi_ = psi[n];
      for (unsigned i = 0; i < n_flux; i++)
      {
        u[i] += this->node_pt(n)->value(flux_nodal_index[i]) * psi_;
      }
    }
  }

  //=====================================================================
  /// Add all the data that are used to interpolate the unknowns. This must
  /// be consistent with the interpolated_u function above and the default
  /// implementation assumes pure nodal interpolaton.
  //=====================================================================
  void DGFaceElement::get_interpolation_data(Vector<Data*>& interpolation_data)
  {
    // Find the number of nodes
    const unsigned n_node = nnode();
    // Now resize the vector
    interpolation_data.resize(n_node);

    // If there are no nodes then return immediately
    if (n_node == 0)
    {
      return;
    }

    // Otherwise loop over the nodes and add to the vector in order
    for (unsigned n = 0; n < n_node; n++)
    {
      interpolation_data[n] = this->node_pt(n);
    }
  }

  //====================================================================
  /// Calculate the numerical flux at the knot point ipt. This is the
  /// most general interface than can be overloaded if desired. The shape
  /// functions at the knot point will be passed into this function.
  //====================================================================
  void DGFaceElement::numerical_flux_at_knot(const unsigned& ipt,
                                             const Shape& psi,
                                             Vector<double>& flux,
                                             DenseMatrix<double>& dflux_du_int,
                                             DenseMatrix<double>& dflux_du_ext,
                                             unsigned flag)
  {
    // Find the number of nodes
    const unsigned n_node = this->nnode();
    // Find the nodal dimension
    const unsigned nodal_dim = this->nodal_dimension();
    // Number of fluxes
    const unsigned n_flux = this->required_nflux();
    // Find the indices at which the local fluxes are stored
    Vector<unsigned> flux_nodal_index(n_flux);
    for (unsigned i = 0; i < n_flux; i++)
    {
      flux_nodal_index[i] = this->flux_index(i);
    }

    // Calculate the local unknowns
    Vector<double> interpolated_u(n_flux, 0.0);

    // Loop over the shape functions
    for (unsigned l = 0; l < n_node; l++)
    {
      // Cache the shape functions
      const double psi_ = psi(l);
      // Loop over the fluxes
      for (unsigned i = 0; i < n_flux; i++)
      {
        // Calculate the velocity from the most recent nodal values
        interpolated_u[i] += this->nodal_value(l, flux_nodal_index[i]) * psi_;
      }
    }

    // Now calculate the outer unit normal Vector
    Vector<double> interpolated_n(nodal_dim);
    this->outer_unit_normal(ipt, interpolated_n);

    // Get the pointer to the neighbour
    DGFaceElement* neighbour_element_pt =
      dynamic_cast<DGFaceElement*>(Neighbour_face_pt[ipt]);

    // Get the neighbour's fluxes
    Vector<double> interpolated_u_neigh(n_flux);

    neighbour_element_pt->interpolated_u(Neighbour_local_coordinate[ipt],
                                         interpolated_u_neigh);

    // Call the "standard" numerical flux function
    this->numerical_flux(
      interpolated_n, interpolated_u, interpolated_u_neigh, flux);

    // If we are calculating the jacobian add this term
    if (flag && (flag < 3))
    {
      this->dnumerical_flux_du(interpolated_n,
                               interpolated_u,
                               interpolated_u_neigh,
                               dflux_du_int,
                               dflux_du_ext);
    }
  }

  //===============================================================
  /// Calculate the derivative of the
  /// normal flux, which is the dot product of our
  /// approximation to the flux with the outer unit normal,
  /// with respect to the interior and exterior variables
  /// Default is to use finite differences
  //=============================================================
  void DGFaceElement::dnumerical_flux_du(const Vector<double>& n_out,
                                         const Vector<double>& u_int,
                                         const Vector<double>& u_ext,
                                         DenseMatrix<double>& dflux_du_int,
                                         DenseMatrix<double>& dflux_du_ext)
  {
    // Find the number of fluxes
    const unsigned n_flux = this->required_nflux();

    // Get a local copy of the unknowns
    Vector<double> u_int_local = u_int;
    Vector<double> u_ext_local = u_ext;

    // Storage for incremented and decremented fluxes
    Vector<double> flux_plus(n_flux), flux_minus(n_flux);

    const double fd_step = GeneralisedElement::Default_fd_jacobian_step;

    // Now loop over all the fluxes
    for (unsigned n = 0; n < n_flux; n++)
    {
      // Increase internal value

      // Store the old value
      double old_var = u_int_local[n];
      // Increment the value
      u_int_local[n] += fd_step;
      // Get the new values
      this->numerical_flux(n_out, u_int_local, u_ext_local, flux_plus);

      // Reset the value
      u_int_local[n] = old_var;
      // Decrement the value
      u_int_local[n] -= fd_step;
      // Get the new values
      this->numerical_flux(n_out, u_int_local, u_ext_local, flux_minus);

      // Assemble the column of the jacobian
      for (unsigned m = 0; m < n_flux; m++)
      {
        dflux_du_int(m, n) = (flux_plus[m] - flux_minus[m]) / (2.0 * fd_step);
      }

      // Reset the value
      u_int_local[n] = old_var;

      // Increase external value

      // Store the old value
      old_var = u_ext_local[n];
      // Increment the value
      u_ext_local[n] += fd_step;
      // Get the new values
      this->numerical_flux(n_out, u_int_local, u_ext_local, flux_plus);

      // Reset the value
      u_ext_local[n] = old_var;
      // Decrement the value
      u_ext_local[n] -= fd_step;
      // Get the new values
      this->numerical_flux(n_out, u_int_local, u_ext_local, flux_minus);

      // Assemble the column of the jacobian
      for (unsigned m = 0; m < n_flux; m++)
      {
        dflux_du_ext(m, n) = (flux_plus[m] - flux_minus[m]) / (2.0 * fd_step);
      }

      // Reset the value
      u_ext_local[n] = old_var;
    }
  }


  //===================================================================
  /// Calculate the integrated (numerical) flux out of the face and add
  /// it to the residuals vector
  //===================================================================
  void DGFaceElement::add_flux_contributions(Vector<double>& residuals,
                                             DenseMatrix<double>& jacobian,
                                             unsigned flag)
  {
// Check that we have set up the coupling data if we are computing the jacobin
#ifdef PARANOID
    if (flag && (flag < 3))
    {
      if (Neighbour_external_data.size() == 0)
      {
        std::ostringstream error_stream;
        error_stream
          << "Coupling data between elements not included in jacobian\n"
          << "You should call DGMesh::setup_face_neighbour_info(true) to "
             "ensure\n"
          << "that this information is included in the jacobian\n";

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif


    // Find the number of nodes
    const unsigned n_node = nnode();
    // Dimension of the face
    const unsigned el_dim = dim();
    // Storage for the shape functions
    Shape psi(n_node);

    // Number of integration points
    const unsigned n_intpt = this->integral_pt()->nweight();
    // Number of fluxes
    const unsigned n_flux = this->required_nflux();
    // Find the indices at which the local fluxes are stored
    Vector<unsigned> flux_nodal_index(n_flux);
    for (unsigned i = 0; i < n_flux; i++)
    {
      flux_nodal_index[i] = this->flux_index(i);
    }

    // Cache the bulk element
    DGElement* bulk_elem_pt = dynamic_cast<DGElement*>(this->bulk_element_pt());

    // Storage for the flux and its derivatives
    Vector<double> F(n_flux);
    DenseMatrix<double> dF_du_int(n_flux, n_flux);
    DenseMatrix<double> dF_du_ext(n_flux, n_flux);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double W = this->integral_pt()->weight(ipt);
      // Get the shape functions at the knot
      this->shape_at_knot(ipt, psi);

      // Calculate the Jacobian
      // For a point element, it's one
      double J = W;
      // Otherwise calculate for the element
      if (el_dim != 0)
      {
        J *= this->J_eulerian_at_knot(ipt);
      }

      // Now calculate the numerical flux (and derivatives)
      this->numerical_flux_at_knot(ipt, psi, F, dF_du_int, dF_du_ext, flag);

      // Limit if desired here

      // Cache the pointer to the neighbour
      DGFaceElement* neighbour_element_pt =
        dynamic_cast<DGFaceElement*>(Neighbour_face_pt[ipt]);

      // Finally  we need to assemble the appropriate contributions
      // to the residuals
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the fluxes
        for (unsigned i = 0; i < n_flux; i++)
        {
          // Get the local equation number in the bulk
          int local_eqn = bulk_elem_pt->nodal_local_eqn(bulk_node_number(l),
                                                        flux_nodal_index[i]);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            // Add the flux multiplied by the shape function and the jacobian
            residuals[local_eqn] -= psi(l) * F[i] * J;

            // Add the jacobian contributions
            if (flag)
            {
              // If we are assembling the jacobian
              if (flag < 3)
              {
                // Loop over the internal nodes and fluxes again
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  for (unsigned i2 = 0; i2 < n_flux; i2++)
                  {
                    // Get the local unknown equation number in the bulk
                    int local_unknown = bulk_elem_pt->nodal_local_eqn(
                      bulk_node_number(l2), flux_nodal_index[i2]);

                    // If it's not a boundary condition
                    if (local_unknown >= 0)
                    {
                      // Add the flux multiplied by the shape function a
                      // nd the jacobian
                      jacobian(local_eqn, local_unknown) -=
                        psi(l) * dF_du_int(i, i2) * psi(l2) * J;
                    }
                  }
                }

                // How many nodes does it have
                unsigned neigh_n_node = neighbour_element_pt->nnode();

                // Get the flux indices from the neighbour
                Vector<unsigned> neigh_flux_index(n_flux);
                for (unsigned i2 = 0; i2 < n_flux; i2++)
                {
                  neigh_flux_index[i2] = neighbour_element_pt->flux_index(i2);
                }

                // Loop over the neighbours nodes
                for (unsigned l2 = 0; l2 < neigh_n_node; l2++)
                {
                  // Loop over the fluxed
                  for (unsigned i2 = 0; i2 < n_flux; i2++)
                  {
                    // Get the local unknown equation number in the bulk
                    int local_unknown =
                      dynamic_cast<DGElement*>(this->bulk_element_pt())
                        ->external_local_eqn(Neighbour_external_data[ipt][l2],
                                             flux_nodal_index[i2]);

                    // If it's not a boundary condition
                    if (local_unknown >= 0)
                    {
                      // Add the flux multiplied by the shape function
                      // and the jacobian
                      jacobian(local_eqn, local_unknown) -=
                        psi(l) * dF_du_ext(i, i2) * psi(l2) * J;
                    }
                  }
                }
              }
            } // End of jacobian calculation
          }
        }
      }
    }
  }

  //========================================================================
  /// Function that computes and stores the (inverse) mass matrix
  //========================================================================
  void DGElement::pre_compute_mass_matrix()
  {
    // Now let's assemble stuff
    const unsigned n_dof = this->ndof();
    // Allocate storage for the local mass matrix (if required)
    if (M_pt == 0)
    {
      M_pt = new DenseDoubleMatrix;
    }

    // Resize and initialise the vector that will holds the residuals
    Vector<double> dummy(n_dof, 0.0);

    // Resize the mass matrix
    M_pt->resize(n_dof, n_dof);
    // Initialise the entries to zero
    M_pt->initialise(0.0);
    // Get the local mass matrix and residuals
    this->fill_in_contribution_to_mass_matrix(dummy, *M_pt);

    // Now invert the mass matrix it will always be small
    // This can possibly be streamlined (for example in spectral
    // elements the mass matrix is diagonal)
    M_pt->ludecompose();

    // The mass matrix has been computed
    Mass_matrix_has_been_computed = true;
  }


  //============================================================================
  /// Function that returns the current value of the residuals
  /// multiplied by the inverse mass matrix (virtual so that it can be
  /// overloaded specific elements in which time and memory saving tricks can be
  /// applied)
  //============================================================================
  void DGElement::get_inverse_mass_matrix_times_residuals(
    Vector<double>& minv_res)
  {
    // If there are external data this is not going to work
    if (nexternal_data() > 0)
    {
      std::ostringstream error_stream;
      error_stream
        << "Cannot use a discontinuous formulation for the mass matrix when\n"
        << "there are external data.\n "
        << "Do not call Problem::enable_discontinuous_formulation()\n";

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Now let's assemble stuff
    const unsigned n_dof = this->ndof();
    // Allocate storage for the local mass matrix (if required)
    if (M_pt == 0)
    {
      M_pt = new DenseDoubleMatrix;
    }

    // Resize and initialise the vector that will holds the residuals
    minv_res.resize(n_dof);
    for (unsigned n = 0; n < n_dof; n++)
    {
      minv_res[n] = 0.0;
    }

    // If we are recycling the mass matrix
    if (Mass_matrix_reuse_is_enabled && Mass_matrix_has_been_computed)
    {
      // Get the residuals
      this->fill_in_contribution_to_residuals(minv_res);
    }
    // Otherwise
    else
    {
      // Resize the mass matrix
      M_pt->resize(n_dof, n_dof);
      // Initialise the entries to zero
      M_pt->initialise(0.0);
      // Get the local mass matrix and residuals
      this->fill_in_contribution_to_mass_matrix(minv_res, *M_pt);

      // Now invert the mass matrix it will always be small
      // This can possibly be streamlined (for example in spectral
      // elements the mass matrix is diagonal)
      M_pt->ludecompose();

      // The mass matrix has been computed
      Mass_matrix_has_been_computed = true;
    }

    // Always do the backsubstitution
    M_pt->lubksub(minv_res);
  }


  void DGElement::get_neighbouring_face_and_local_coordinate(
    const int& face_index,
    const Vector<double>& s,
    FaceElement*& face_element_pt,
    Vector<double>& s_face)
  {
    DG_mesh_pt->neighbour_finder(this, face_index, s, face_element_pt, s_face);
  }


  /// Limit the slope within an element
  void DGElement::slope_limit(SlopeLimiter* const& slope_limiter_pt)
  {
    // Firstly find the dimension
    const unsigned n_dim = this->dim();
    // Find the number of fluxes
    const unsigned n_flux = this->required_nflux();

    switch (n_dim)
    {
        // One dimensional (easy-ish) case
      case 1:
      {
        // Storage for the element and its neighbours
        Vector<DGElement*> required_element_pt(3);
        required_element_pt[0] = this;

        // Get the pointer to the element on the left
        required_element_pt[1] = dynamic_cast<DGElement*>(
          dynamic_cast<DGFaceElement*>(this->face_element_pt(0))
            ->neighbour_face_pt(0)
            ->bulk_element_pt());
        // Get the pointer to the element on the right
        required_element_pt[2] = dynamic_cast<DGElement*>(
          dynamic_cast<DGFaceElement*>(this->face_element_pt(1))
            ->neighbour_face_pt(0)
            ->bulk_element_pt());

        // Loop over the fluxed
        for (unsigned i = 0; i < n_flux; i++)
        {
          // Call our limiter, which will take as it's arguments, the current
          // element and the required neighbours
          slope_limiter_pt->limit(i, required_element_pt);
        }
      }
      break;

      default:
      {
        std::ostringstream error_stream;
        error_stream << "Slope limiting is not implemented for this dimension: "
                     << n_dim;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
  }

  double DGMesh::FaceTolerance = 1.0e-10;


  //====================================================
  /// Helper minmod function
  //====================================================
  double MinModLimiter::minmod(Vector<double>& args)
  {
    const unsigned n_arg = args.size();
    // If no arguments return zero
    if (n_arg == 0)
    {
      return 0.0;
    }

    // Initialise the sign from the sign of the first entry
    int sign = 0;
    if (args[0] < 0.0)
    {
      sign = -1;
    }
    else if (args[0] > 0.0)
    {
      sign = 1;
    }
    else
    {
      return 0.0;
    }

    // Initialise the minimum value
    double min = std::fabs(args[0]);

    // Now loop over the rest of the values
    for (unsigned i = 1; i < n_arg; i++)
    {
      if (sign == 1)
      {
        if (args[i] < 0.0)
        {
          return 0.0;
        }
        else if (args[i] < min)
        {
          min = args[i];
        }
      }
      else if (sign == -1)
      {
        if (args[i] > 0.0)
        {
          return 0.0;
        }
        else if (std::fabs(args[i]) < min)
        {
          min = std::fabs(args[i]);
        }
      }
    }

    // Make sure to return the sign multiplied by the minimum value
    return sign * min;
  }


  /// Modified minmod limiter to fix behaviour in smooth regions
  double MinModLimiter::minmodB(Vector<double>& args, const double& h)
  {
    const unsigned n_arg = args.size();
    // If no arguments return zero
    if (n_arg == 0)
    {
      return 0.0;
    }

    // Modification to fix extrema
    if (std::fabs(args[0]) < this->M * h * h)
    {
      return args[0];
    }
    // Otherwise just return the usual minmod
    else
    {
      return minmod(args);
    }
  }

  /// Implement the limiter function for the basic MinModlimiter
  void MinModLimiter::limit(const unsigned& i,
                            const Vector<DGElement*>& required_element_pt)
  {
    // Set the tolerance
    const double tol = 1.0e-16;

    // Find the geometric parameters of the element
    const unsigned n_node = required_element_pt[0]->nnode();
    const double x_l = required_element_pt[0]->node_pt(0)->x(0);
    const double x_r = required_element_pt[0]->node_pt(n_node - 1)->x(0);
    const double h = x_r - x_l;
    const double x0 = 0.5 * (x_l + x_r);
    // Find the average value
    const double u_av = required_element_pt[0]->average_value(i);

    // Storage for the gradients to the minmod function
    Vector<double> arg;
    arg.reserve(3);

    // Add the approximation for the element's left flux
    arg.push_back(u_av - required_element_pt[0]->node_pt(0)->value(i));

    // If there is a left element add the approximate gradient
    // to the argument vector
    if (required_element_pt[1] != required_element_pt[0])
    {
      arg.push_back(u_av - required_element_pt[1]->average_value(i));
    }
    // If there is a right element add the approximate gradient
    // to the argument vector
    if (required_element_pt[2] != required_element_pt[0])
    {
      arg.push_back(required_element_pt[2]->average_value(i) - u_av);
    }

    // Set the left value
    const double u_l = u_av - this->minmod(arg);

    // Now replace the first term in the argument list with the
    // approximation for the element's right flux
    arg.front() = required_element_pt[0]->node_pt(n_node - 1)->value(i) - u_av;

    // Set the right value
    const double u_r = u_av + this->minmod(arg);

    // If the basic limited values are different from
    // the unlimited values then limit
    if ((std::fabs(u_l - required_element_pt[0]->node_pt(0)->value(i)) > tol) &&
        (std::fabs(
           u_r - required_element_pt[0]->node_pt(n_node - 1)->value(i)) > tol))
    {
      // Find the centre of the element on the left
      const double x0_l =
        0.5 * (required_element_pt[1]
                 ->node_pt(required_element_pt[1]->nnode() - 1)
                 ->x(0) +
               required_element_pt[1]->node_pt(0)->x(0));
      // Find the centre of the element on the right
      const double x0_r =
        0.5 * (required_element_pt[2]
                 ->node_pt(required_element_pt[2]->nnode() - 1)
                 ->x(0) +
               required_element_pt[2]->node_pt(0)->x(0));

      // Clear the argument list and reserve its size to
      arg.clear();
      arg.reserve(3);

      // Approximate the gradient over the whole cell
      arg.push_back((required_element_pt[0]->node_pt(n_node - 1)->value(i) -
                     required_element_pt[0]->node_pt(0)->value(i)) /
                    h);

      // Adjust the estimates for the gradient calculated from neighbouring
      // averages for the straight min-mod and MUSCL limiters
      double gradient_factor = 0.5;
      if (MUSCL)
      {
        gradient_factor = 1.0;
      }

      // If there is a left element, form the gradient
      if (required_element_pt[0] != required_element_pt[1])
      {
        arg.push_back((u_av - required_element_pt[1]->average_value(i)) /
                      (gradient_factor * (x0 - x0_l)));
      }

      // If there is a right element, form the gradient
      if (required_element_pt[0] != required_element_pt[2])
      {
        arg.push_back((required_element_pt[2]->average_value(i) - u_av) /
                      (gradient_factor * (x0_r - x0)));
      }

      // Calculate the limited gradient of these three
      double limited_gradient = this->minmodB(arg, h);

      // Loop over the nodes and limit
      for (unsigned n = 0; n < n_node; n++)
      {
        double x = required_element_pt[0]->node_pt(n)->x(0) - x0;
        required_element_pt[0]->node_pt(n)->set_value(
          i, u_av + x * limited_gradient);
      }
    }
  }

} // namespace oomph
