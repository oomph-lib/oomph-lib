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
#ifndef OOMPH_GEN_HH_TIME_HARMONIC_LIN_ELAST_HEADER
#define OOMPH_GEN_HH_TIME_HARMONIC_LIN_ELAST_HEADER


// Oomph-lib includes
#include "generic.h"
#include "pml_helmholtz.h"
#include "time_harmonic_linear_elasticity.h"

namespace oomph
{
  //======================================================================
  /// A class for elements that allow the imposition of an applied traction
  /// in the equations of time-harmonic linear elasticity from a
  /// PMLHelmholtz potential (interpreted as a displacement potential
  /// for the fluid in a  quasi-steady, linearised FSI problem.)
  /// The geometrical information can be read from the FaceGeometry<ELEMENT>
  /// class and thus, we can be generic enough without the need to have
  /// a separate equations class.
  //======================================================================
  template<class ELASTICITY_BULK_ELEMENT, class HELMHOLTZ_BULK_ELEMENT>
  class TimeHarmonicLinElastLoadedByPMLHelmholtzPressureBCElement
    : public virtual FaceGeometry<ELASTICITY_BULK_ELEMENT>,
      public virtual FaceElement,
      public virtual ElementWithExternalElement
  {
  protected:
    /// \short Pointer to the ratio, \f$ Q \f$ , of the stress used to
    /// non-dimensionalise the fluid stresses to the stress used to
    /// non-dimensionalise the solid stresses.
    double* Q_pt;

    /// \short Static default value for the ratio of stress scales
    /// used in the fluid and solid equations (default is 1.0)
    static double Default_Q_Value;

    /// Index at which the i-th displacement component is stored
    Vector<std::complex<unsigned>>
      U_index_time_harmonic_linear_elasticity_helmholtz_traction;

    /// \short Helper function that actually calculates the residuals
    // This small level of indirection is required to avoid calling
    // fill_in_contribution_to_residuals in fill_in_contribution_to_jacobian
    // which causes all kinds of pain if overloading later on
    void fill_in_contribution_to_residuals_helmholtz_traction(
      Vector<double>& residuals);

  public:
    /// \short Constructor, which takes a "bulk" element and the
    /// value of the index and its limit
    TimeHarmonicLinElastLoadedByPMLHelmholtzPressureBCElement(
      FiniteElement* const& element_pt, const int& face_index)
      : FaceGeometry<ELASTICITY_BULK_ELEMENT>(),
        FaceElement(),
        Q_pt(&Default_Q_Value)
    {
      // Attach the geometrical information to the element. N.B. This function
      // also assigns nbulk_value from the required_nvalue of the bulk element
      element_pt->build_face_element(face_index, this);

#ifdef PARANOID
      {
        // Check that the element is not a refineable 3d element
        ELASTICITY_BULK_ELEMENT* elem_pt =
          dynamic_cast<ELASTICITY_BULK_ELEMENT*>(element_pt);
        // If it's three-d
        if (elem_pt->dim() == 3)
        {
          // Is it refineable
          RefineableElement* ref_el_pt =
            dynamic_cast<RefineableElement*>(elem_pt);
          if (ref_el_pt != 0)
          {
            if (this->has_hanging_nodes())
            {
              throw OomphLibError("This flux element will not work correctly "
                                  "if nodes are hanging\n",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }
#endif

      // Set source element storage: one interaction with an external
      // element -- the PMLHelmholtz bulk element that provides
      // the pressure
      this->set_ninteraction(1);

      // Find the dimension of the problem
      unsigned n_dim = element_pt->nodal_dimension();

      // Find the index at which the displacement unknowns are stored
      ELASTICITY_BULK_ELEMENT* cast_element_pt =
        dynamic_cast<ELASTICITY_BULK_ELEMENT*>(element_pt);
      this->U_index_time_harmonic_linear_elasticity_helmholtz_traction.resize(
        n_dim);
      for (unsigned i = 0; i < n_dim; i++)
      {
        this->U_index_time_harmonic_linear_elasticity_helmholtz_traction[i] =
          cast_element_pt->u_index_time_harmonic_linear_elasticity(i);
      }
    }


    /// Return the residuals
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      fill_in_contribution_to_residuals_helmholtz_traction(residuals);
    }


    /// Fill in contribution from Jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the residuals
      fill_in_contribution_to_residuals_helmholtz_traction(residuals);

      // Derivatives w.r.t. external data
      fill_in_jacobian_from_external_interaction_by_fd(residuals, jacobian);
    }

    /// \short Return the ratio of the stress scales used to non-dimensionalise
    /// the fluid and elasticity equations. E.g.
    /// \f$ Q = (\omega a)^2 \rho/E \f$, i.e. the
    /// ratio between the inertial fluid stress and the solid's elastic
    /// modulus E.
    const double& q() const
    {
      return *Q_pt;
    }

    /// \short Return a pointer the ratio of stress scales used to
    /// non-dimensionalise the fluid and solid equations.
    double*& q_pt()
    {
      return Q_pt;
    }


    /// \short Output function
    void output(std::ostream& outfile)
    {
      /// Dummy
      unsigned nplot = 0;
      output(outfile, nplot);
    }

    /// \short Output function: Plot traction etc at Gauss points
    /// nplot is ignored.
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      // Dimension
      unsigned n_dim = this->nodal_dimension();

      // Storage for traction
      Vector<std::complex<double>> traction(n_dim);

      // Get FSI parameter
      const double q_value = q();

      outfile << "ZONE\n";

      // Set the value of n_intpt
      unsigned n_intpt = integral_pt()->nweight();

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        Vector<double> s_int(n_dim - 1);
        for (unsigned i = 0; i < n_dim - 1; i++)
        {
          s_int[i] = integral_pt()->knot(ipt, i);
        }

        // Get the outer unit normal
        Vector<double> interpolated_normal(n_dim);
        outer_unit_normal(ipt, interpolated_normal);

        // Boundary coordinate
        Vector<double> zeta(1);
        interpolated_zeta(s_int, zeta);

        // Get external element for potential
        HELMHOLTZ_BULK_ELEMENT* ext_el_pt =
          dynamic_cast<HELMHOLTZ_BULK_ELEMENT*>(external_element_pt(0, ipt));
        Vector<double> s_ext(external_element_local_coord(0, ipt));
        std::complex<double> u_helmholtz =
          ext_el_pt->interpolated_u_pml_helmholtz(s_ext);

        // Traction: Pressure is proportional to POSITIVE potential
        ext_el_pt->interpolated_u_pml_helmholtz(s_ext);
        traction[0] = -q_value * interpolated_normal[0] * u_helmholtz;
        traction[1] = -q_value * interpolated_normal[1] * u_helmholtz;

        outfile << ext_el_pt->interpolated_x(s_ext, 0) << " "
                << ext_el_pt->interpolated_x(s_ext, 1) << " "
                << traction[0].real() << " " << traction[1].real() << " "
                << traction[0].imag() << " " << traction[1].imag() << " "
                << interpolated_normal[0] << " " << interpolated_normal[1]
                << " " << u_helmholtz.real() << " " << u_helmholtz.imag() << " "
                << interpolated_x(s_int, 0) << " " << interpolated_x(s_int, 1)
                << " "
                << sqrt(pow(ext_el_pt->interpolated_x(s_ext, 0) -
                              interpolated_x(s_int, 0),
                            2) +
                        pow(ext_el_pt->interpolated_x(s_ext, 1) -
                              interpolated_x(s_int, 1),
                            2))
                << " " << zeta[0] << std::endl;
      }
    }


    /// \short C_style output function
    void output(FILE* file_pt)
    {
      FaceGeometry<ELASTICITY_BULK_ELEMENT>::output(file_pt);
    }

    /// \short C-style output function
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FaceGeometry<ELASTICITY_BULK_ELEMENT>::output(file_pt, n_plot);
    }

    /// \short Compute the global_flux_contribution through the template
    /// elasticity elements : we compute u.n
    std::complex<double> global_flux_contribution_from_solid()
    {
      // Dummy output file
      std::ofstream outfile;
      return global_flux_contribution_from_solid(outfile);
    }

    /// \short  Compute the element's contribution to the integral of
    /// the flux over the artificial boundary. Also output the
    /// flux  in the specified output file if it's open.
    std::complex<double> global_flux_contribution_from_solid(
      std::ofstream& outfile)
    {
      // pointer to the corresponding elasticity bulk element
      ELASTICITY_BULK_ELEMENT* bulk_elem_pt =
        dynamic_cast<ELASTICITY_BULK_ELEMENT*>(this->bulk_element_pt());

      // get the dim of the bulk and local nodes
      const unsigned bulk_dim = bulk_elem_pt->dim();
      const unsigned local_dim = this->dim();

      // Set up memory for the outer unit normal
      Vector<double> unit_normal(bulk_dim);

      // Set the value of n_intpt
      const unsigned n_intpt = integral_pt()->nweight();

      // Set the Vector to hold local coordinates
      Vector<double> s(local_dim);
      std::complex<double> flux(0.0, 0.0);
      Vector<std::complex<double>> u(bulk_dim);

      // Output?
      if (outfile.is_open())
      {
        outfile << "ZONE\n";
      }

      // Loop over the integration points
      //--------------------------------
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s
        for (unsigned i = 0; i < local_dim; i++)
        {
          s[i] = integral_pt()->knot(ipt, i);
        }
        // get the outer_unit_ext vector
        this->outer_unit_normal(s, unit_normal);

        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Get jacobian of mapping
        double J = J_eulerian(s);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get local coordinates in bulk element by copy construction
        Vector<double> s_bulk(local_coordinate_in_bulk(s));

        // Get the displacement from the bulk element
        bulk_elem_pt->interpolated_u_time_harmonic_linear_elasticity(s_bulk, u);

        // Compute normal displacement u.n
        std::complex<double> u_dot_n(0.0, 0.0);
        for (unsigned i = 0; i < bulk_dim; i++)
        {
          u_dot_n += u[i] * unit_normal[i];
        }

        // Output?
        if (outfile.is_open())
        {
          Vector<double> x(bulk_dim);
          interpolated_x(s, x);
          outfile << x[0] << " " << x[1] << " " << u_dot_n.real() << " "
                  << u_dot_n.imag() << "\n";
        }

        // ...add to integral
        flux += u_dot_n * W;
      }

      return flux;
    }


    /// \short Compute the global_flux_contribution through the helmholtz
    /// elements : we compute dphidn.n
    std::complex<double> global_flux_contribution_from_helmholtz()
    {
      // Dummy output file
      std::ofstream outfile;
      return global_flux_contribution_from_helmholtz(outfile);
    }

    /// \short  Compute the element's contribution to the integral of
    /// the flux over the artificial boundary. Also output the
    /// flux  in the specified output file if it's open.
    std::complex<double> global_flux_contribution_from_helmholtz(
      std::ofstream& outfile)
    {
      // Get the dim of this element
      const unsigned local_dim = this->dim();

      // Set the value of n_intpt
      const unsigned n_intpt = integral_pt()->nweight();

      // Set the Vector to hold local coordinates
      Vector<double> s(local_dim);
      std::complex<double> flux(0.0, 0.0);

      // Output?
      if (outfile.is_open())
      {
        outfile << "ZONE\n";
      }

      // Loop over the integration points
      //--------------------------------
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s
        for (unsigned i = 0; i < local_dim; i++)
        {
          s[i] = integral_pt()->knot(ipt, i);
        }

        // Get local coordinates in bulk element by copy construction
        Vector<double> s_bulk(local_coordinate_in_bulk(s));

        // Get external element for potential
        HELMHOLTZ_BULK_ELEMENT* ext_el_pt =
          dynamic_cast<HELMHOLTZ_BULK_ELEMENT*>(external_element_pt(0, ipt));

        // number of nodes in ext element
        unsigned nnode_ext = ext_el_pt->nnode();

        // get the dim of the external nodes
        const unsigned ext_dim = ext_el_pt->dim();

        // Set up memory for the external shape function
        Shape psi_ext(nnode_ext);
        DShape dpsi_ext_dx(nnode_ext, ext_dim);

        // Call the derivatives of the shape  functions
        // in the external element -- must do this via s because this point
        // is not an integration point the bulk element!
        (void)ext_el_pt->dshape_eulerian(s_bulk, psi_ext, dpsi_ext_dx);

        // Set up memory for the outer unit normal
        Vector<double> unit_normal(ext_dim);

        // get the outer_unit_ext vector
        this->outer_unit_normal(s, unit_normal);

        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Get jacobian of mapping
        double J = J_eulerian(s);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Global gradient
        Vector<std::complex<double>> interpolated_dphidx(
          ext_dim, std::complex<double>(0.0, 0.0));
        for (unsigned l = 0; l < nnode_ext; l++)
        {
          // Get the nodal value of the helmholtz unknown
          const std::complex<double> phi_value(
            ext_el_pt->nodal_value(l, ext_el_pt->u_index_helmholtz().real()),
            ext_el_pt->nodal_value(l, ext_el_pt->u_index_helmholtz().imag()));

          // Loop over directions
          for (unsigned i = 0; i < ext_dim; i++)
          {
            interpolated_dphidx[i] += phi_value * dpsi_ext_dx(l, i);
          }
        }

        // define dphi_dn
        std::complex<double> dphi_dn(0.0, 0.0);
        for (unsigned i = 0; i < ext_dim; i++)
        {
          dphi_dn += interpolated_dphidx[i] * unit_normal[i];
        }


#ifdef PARANOID
        double max_permitted_discrepancy = 1.0e-10;
        double diff = sqrt(
          pow(ext_el_pt->interpolated_x(s_bulk, 0) - interpolated_x(s, 0), 2) +
          pow(ext_el_pt->interpolated_x(s_bulk, 1) - interpolated_x(s, 1), 2));
        if (diff > max_permitted_discrepancy)
        {
          std::ostringstream error_stream;
          error_stream
            << "Position computed by external el and face element differ by "
            << diff << " which is more than the acceptable tolerance "
            << max_permitted_discrepancy << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Output?
        if (outfile.is_open())
        {
          Vector<double> x(ext_dim);
          interpolated_x(s, x);
          outfile << x[0] << " " << x[1] << " " << dphi_dn.real() << " "
                  << dphi_dn.imag() << " "
                  << "\n";
        }

        // ...add to integral
        flux += dphi_dn * W;
      }

      return flux;
    }
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Static default value for the ragoogletio of stress scales
  /// used in the fluid and solid equations (default is 1.0)
  //=================================================================
  template<class ELASTICITY_BULK_ELEMENT, class HELMHOLTZ_BULK_ELEMENT>
  double TimeHarmonicLinElastLoadedByPMLHelmholtzPressureBCElement<
    ELASTICITY_BULK_ELEMENT,
    HELMHOLTZ_BULK_ELEMENT>::Default_Q_Value = 1.0;


  //=====================================================================
  /// Return the residuals
  //=====================================================================
  template<class ELASTICITY_BULK_ELEMENT, class HELMHOLTZ_BULK_ELEMENT>
  void TimeHarmonicLinElastLoadedByPMLHelmholtzPressureBCElement<
    ELASTICITY_BULK_ELEMENT,
    HELMHOLTZ_BULK_ELEMENT>::
    fill_in_contribution_to_residuals_helmholtz_traction(
      Vector<double>& residuals)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

#ifdef PARANOID
    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();
    if (n_position_type != 1)
    {
      throw OomphLibError("LinearElasticity is not yet implemented for more "
                          "than one position type.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find out the dimension of the node
    const unsigned n_dim = this->nodal_dimension();

    // Cache the nodal indices at which the displacement components are stored
    std::complex<unsigned> u_nodal_index[n_dim];
    for (unsigned i = 0; i < n_dim; i++)
    {
      u_nodal_index[i] =
        this->U_index_time_harmonic_linear_elasticity_helmholtz_traction[i];
    }

    // Integer to hold the local equation number
    int local_eqn = 0;

    // Set up memory for the shape functions
    Shape psi(n_node);
    DShape dpsids(n_node, n_dim - 1);

    // Get FSI parameter
    const double q_value = q();

    // Storage for traction
    Vector<std::complex<double>> traction(2);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Only need to call the local derivatives
      dshape_local_at_knot(ipt, psi, dpsids);

      // Calculate the coordinates
      Vector<double> interpolated_x(n_dim, 0.0);

      // Also calculate the surface tangent vectors
      DenseMatrix<double> interpolated_A(n_dim - 1, n_dim, 0.0);

      // Calculate displacements and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned i = 0; i < n_dim; i++)
        {
          // Calculate the Eulerian coords
          const double x_local = nodal_position(l, i);
          interpolated_x[i] += x_local * psi(l);

          // Loop over LOCAL derivative directions, to calculate the tangent(s)
          for (unsigned j = 0; j < n_dim - 1; j++)
          {
            interpolated_A(j, i) += x_local * dpsids(l, j);
          }
        }
      }

      // Now find the local metric tensor from the tangent Vectors
      DenseMatrix<double> A(n_dim - 1);
      for (unsigned i = 0; i < n_dim - 1; i++)
      {
        for (unsigned j = 0; j < n_dim - 1; j++)
        {
          // Initialise surface metric tensor to zero
          A(i, j) = 0.0;

          // Take the dot product
          for (unsigned k = 0; k < n_dim; k++)
          {
            A(i, j) += interpolated_A(i, k) * interpolated_A(j, k);
          }
        }
      }

      // Get the outer unit normal
      Vector<double> interpolated_normal(n_dim);
      outer_unit_normal(ipt, interpolated_normal);

      // Find the determinant of the metric tensor
      double Adet = 0.0;
      switch (n_dim)
      {
        case 2:
          Adet = A(0, 0);
          break;
        case 3:
          Adet = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
          break;
        default:
          throw OomphLibError(
            "Wrong dimension in TimeHarmonicLinElastLoadedByPressureElement",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      // Premultiply the weights and the square-root of the determinant of
      // the metric tensor
      double W = w * sqrt(Adet);

      // Get bulk element for potential
      HELMHOLTZ_BULK_ELEMENT* ext_el_pt =
        dynamic_cast<HELMHOLTZ_BULK_ELEMENT*>(external_element_pt(0, ipt));
      Vector<double> s_ext(external_element_local_coord(0, ipt));

      // Traction: Pressure is proportional to POSITIVE potential
      std::complex<double> u_helmholtz =
        ext_el_pt->interpolated_u_pml_helmholtz(s_ext);
      traction[0] = -q_value * interpolated_normal[0] * u_helmholtz;
      traction[1] = -q_value * interpolated_normal[1] * u_helmholtz;

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the displacement components
        for (unsigned i = 0; i < n_dim; i++)
        {
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[i].real());
          /*IF it's not a boundary condition*/
          if (local_eqn >= 0)
          {
            // Add the loading terms to the residuals
            residuals[local_eqn] -= traction[i].real() * psi(l) * W;
          }

          local_eqn = this->nodal_local_eqn(l, u_nodal_index[i].imag());
          /*IF it's not a boundary condition*/
          if (local_eqn >= 0)
          {
            // Add the loading terms to the residuals
            residuals[local_eqn] -= traction[i].imag() * psi(l) * W;
          }

        } // End of if not boundary condition
      } // End of loop over shape functions
    } // End of loop over integration points
  }


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //======================================================================
  /// \short A class for elements that allow the imposition of an
  /// prescribed flux (determined from the normal displacements of an
  /// adjacent linearly elastic solid. Normal derivative for displacement
  /// potential is given by normal displacement of adjacent solid multiplies
  /// by FSI parameter (q = k^2 B/E).
  /// The element geometry is obtained from the FaceGeometry<ELEMENT>
  /// policy class.
  //======================================================================
  template<class HELMHOLTZ_BULK_ELEMENT, class ELASTICITY_BULK_ELEMENT>
  class PMLHelmholtzFluxFromNormalDisplacementBCElement
    : public virtual FaceGeometry<HELMHOLTZ_BULK_ELEMENT>,
      public virtual FaceElement,
      public virtual ElementWithExternalElement
  {
  public:
    /// \short Constructor, takes the pointer to the "bulk" element and the
    /// face index identifying the face to which the element is attached.
    PMLHelmholtzFluxFromNormalDisplacementBCElement(
      FiniteElement* const& bulk_el_pt, const int& face_index);

    /// Broken copy constructor
    PMLHelmholtzFluxFromNormalDisplacementBCElement(
      const PMLHelmholtzFluxFromNormalDisplacementBCElement& dummy)
    {
      BrokenCopy::broken_copy(
        "PMLHelmholtzFluxFromNormalDisplacementBCElement");
    }

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const
                   PMLHelmholtzFluxFromNormalDisplacementBCElement&)
     {
      BrokenCopy::broken_assign
       ("PMLHelmholtzFluxFromNormalDisplacementBCElement");
       }*/


    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      fill_in_generic_residual_contribution_helmholtz_flux_from_displacement(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }


    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                                 DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_helmholtz_flux_from_displacement(
        residuals, jacobian, 1);

      // Derivatives w.r.t. external data
      fill_in_jacobian_from_external_interaction_by_fd(residuals, jacobian);
    }

    /// Output function
    void output(std::ostream& outfile)
    {
      // Dummy
      unsigned nplot = 0;
      output(outfile, nplot);
    }

    /// Output function: flux etc at Gauss points; n_plot is ignored.
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      outfile << "ZONE\n";

      // Get the value of Nintpt
      const unsigned n_intpt = integral_pt()->nweight();

      // Set the Vector to hold local coordinates
      Vector<double> s_int(Dim - 1);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s
        for (unsigned i = 0; i < (Dim - 1); i++)
        {
          s_int[i] = integral_pt()->knot(ipt, i);
        }

        // Get the unit normal
        Vector<double> interpolated_normal(Dim);
        outer_unit_normal(ipt, interpolated_normal);

        Vector<double> zeta(1);
        interpolated_zeta(s_int, zeta);

        // Get displacements
        ELASTICITY_BULK_ELEMENT* ext_el_pt =
          dynamic_cast<ELASTICITY_BULK_ELEMENT*>(external_element_pt(0, ipt));
        Vector<double> s_ext(external_element_local_coord(0, ipt));
        Vector<std::complex<double>> displ(2);
        ext_el_pt->interpolated_u_time_harmonic_linear_elasticity(s_ext, displ);

        // Convert into flux BC: This takes the dot product of the
        // actual displacement with the flux element's own outer
        // unit normal so the plus sign is OK.
        std::complex<double> flux = (displ[0] * interpolated_normal[0] +
                                     displ[1] * interpolated_normal[1]);

        // Output
        outfile << ext_el_pt->interpolated_x(s_ext, 0) << " "
                << ext_el_pt->interpolated_x(s_ext, 1) << " "
                << flux.real() * interpolated_normal[0] << " "
                << flux.real() * interpolated_normal[1] << " "
                << flux.imag() * interpolated_normal[0] << " "
                << flux.imag() * interpolated_normal[1] << " "
                << interpolated_normal[0] << " " << interpolated_normal[1]
                << " " << flux.real() << " " << flux.imag() << " "
                << interpolated_x(s_int, 0) << " " << interpolated_x(s_int, 1)
                << " "
                << sqrt(pow(ext_el_pt->interpolated_x(s_ext, 0) -
                              interpolated_x(s_int, 0),
                            2) +
                        pow(ext_el_pt->interpolated_x(s_ext, 1) -
                              interpolated_x(s_int, 1),
                            2))
                << " " << zeta[0] << std::endl;
      }
    }


    /// C-style output function
    void output(FILE* file_pt)
    {
      FaceGeometry<HELMHOLTZ_BULK_ELEMENT>::output(file_pt);
    }

    /// C-style output function
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FaceGeometry<HELMHOLTZ_BULK_ELEMENT>::output(file_pt, n_plot);
    }


  protected:
    /// \short Function to compute the shape and test functions and to return
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test(const Vector<double>& s,
                                 Shape& psi,
                                 Shape& test) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Get the shape functions
      shape(s, psi);

      // Set the test functions to be the same as the shape functions
      for (unsigned i = 0; i < n_node; i++)
      {
        test[i] = psi[i];
      }

      // Return the value of the jacobian
      return J_eulerian(s);
    }


    /// \short Function to compute the shape and test functions and to return
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test_at_knot(const unsigned& ipt,
                                         Shape& psi,
                                         Shape& test) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Get the shape functions
      shape_at_knot(ipt, psi);

      // Set the test functions to be the same as the shape functions
      for (unsigned i = 0; i < n_node; i++)
      {
        test[i] = psi[i];
      }

      // Return the value of the jacobian
      return J_eulerian_at_knot(ipt);
    }


  private:
    /// \short Add the element's contribution to its residual vector.
    /// flag=1(or 0): do (or don't) compute the contribution to the
    /// Jacobian as well.
    void fill_in_generic_residual_contribution_helmholtz_flux_from_displacement(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag);

    /// The spatial dimension of the problem
    unsigned Dim;

    /// The index at which the unknown is stored at the nodes
    std::complex<unsigned> U_index_helmholtz_from_displacement;
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////


  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, and the
  /// face index that identifies the face of the bulk element to which
  /// this face element is to be attached.
  //===========================================================================
  template<class HELMHOLTZ_BULK_ELEMENT, class ELASTICITY_BULK_ELEMENT>
  PMLHelmholtzFluxFromNormalDisplacementBCElement<HELMHOLTZ_BULK_ELEMENT,
                                                  ELASTICITY_BULK_ELEMENT>::
    PMLHelmholtzFluxFromNormalDisplacementBCElement(
      FiniteElement* const& bulk_el_pt, const int& face_index)
    : FaceGeometry<HELMHOLTZ_BULK_ELEMENT>(), FaceElement()
  {
    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index, this);

#ifdef PARANOID
    {
      // Check that the element is not a refineable 3d element
      HELMHOLTZ_BULK_ELEMENT* elem_pt =
        dynamic_cast<HELMHOLTZ_BULK_ELEMENT*>(bulk_el_pt);
      // If it's three-d
      if (elem_pt->dim() == 3)
      {
        // Is it refineable
        RefineableElement* ref_el_pt =
          dynamic_cast<RefineableElement*>(elem_pt);
        if (ref_el_pt != 0)
        {
          if (this->has_hanging_nodes())
          {
            throw OomphLibError("This flux element will not work correctly if "
                                "nodes are hanging\n",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
    }
#endif

    // Set source element storage: one interaction with an external element
    // that provides the displacement of the adjacent linear elasticity
    // element
    this->set_ninteraction(1);

    // Extract the dimension of the problem from the dimension of
    // the first node
    Dim = this->node_pt(0)->ndim();

    // Set up U_index_helmholtz_displacement. Initialise to zero, which
    // probably won't change in most cases, oh well, the price we
    // pay for generality
    U_index_helmholtz_from_displacement = std::complex<unsigned>(0, 0);

    // Cast to the appropriate PMLHelmholtzEquation so that we can
    // find the index at which the variable is stored
    // We assume that the dimension of the full problem is the same
    // as the dimension of the node, if this is not the case you will have
    // to write custom elements, sorry
    switch (Dim)
    {
        // One dimensional problem
      case 1:
      {
        PMLHelmholtzEquations<1>* eqn_pt =
          dynamic_cast<PMLHelmholtzEquations<1>*>(bulk_el_pt);
        // If the cast has failed die
        if (eqn_pt == 0)
        {
          std::string error_string =
            "Bulk element must inherit from PMLHelmholtzEquations.";
          error_string +=
            "Nodes are one dimensional, but cannot cast the bulk element to\n";
          error_string += "PMLHelmholtzEquations<1>\n.";
          error_string += "If you desire this functionality, you must "
                          "implement it yourself\n";

          throw OomphLibError(
            error_string, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
        // Otherwise read out the value
        else
        {
          // Read the index from the (cast) bulk element
          U_index_helmholtz_from_displacement = eqn_pt->u_index_helmholtz();
        }
      }
      break;

      // Two dimensional problem
      case 2:
      {
        PMLHelmholtzEquations<2>* eqn_pt =
          dynamic_cast<PMLHelmholtzEquations<2>*>(bulk_el_pt);
        // If the cast has failed die
        if (eqn_pt == 0)
        {
          std::string error_string =
            "Bulk element must inherit from PMLHelmholtzEquations.";
          error_string +=
            "Nodes are two dimensional, but cannot cast the bulk element to\n";
          error_string += "PMLHelmholtzEquations<2>\n.";
          error_string += "If you desire this functionality, you must "
                          "implement it yourself\n";

          throw OomphLibError(
            error_string, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
        else
        {
          // Read the index from the (cast) bulk element.
          U_index_helmholtz_from_displacement = eqn_pt->u_index_helmholtz();
        }
      }
      break;

      // Three dimensional problem
      case 3:
      {
        PMLHelmholtzEquations<3>* eqn_pt =
          dynamic_cast<PMLHelmholtzEquations<3>*>(bulk_el_pt);
        // If the cast has failed die
        if (eqn_pt == 0)
        {
          std::string error_string =
            "Bulk element must inherit from PMLHelmholtzEquations.";
          error_string += "Nodes are three dimensional, but cannot cast the "
                          "bulk element to\n";
          error_string += "PMLHelmholtzEquations<3>\n.";
          error_string += "If you desire this functionality, you must "
                          "implement it yourself\n";

          throw OomphLibError(
            error_string, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
        else
        {
          // Read the index from the (cast) bulk element.
          U_index_helmholtz_from_displacement = eqn_pt->u_index_helmholtz();
        }
      }
      break;

      // Any other case is an error
      default:
        std::ostringstream error_stream;
        error_stream << "Dimension of node is " << Dim
                     << ". It should be 1,2, or 3!" << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        break;
    }
  }


  //===========================================================================
  /// \short Helper function to compute the element's residual vector and
  /// the Jacobian matrix.
  //===========================================================================
  template<class HELMHOLTZ_BULK_ELEMENT, class ELASTICITY_BULK_ELEMENT>
  void PMLHelmholtzFluxFromNormalDisplacementBCElement<
    HELMHOLTZ_BULK_ELEMENT,
    ELASTICITY_BULK_ELEMENT>::
    fill_in_generic_residual_contribution_helmholtz_flux_from_displacement(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
  {
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psif(n_node), testf(n_node);

    // Set the value of Nintpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(Dim - 1);

    // Integers to hold the local equation and unknown numbers
    int local_eqn = 0;

    // Locally cache the index at which the variable is stored
    const std::complex<unsigned> u_index_helmholtz =
      U_index_helmholtz_from_displacement;

    // Loop over the integration points
    //--------------------------------
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < (Dim - 1); i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Find the shape and test functions and return the Jacobian
      // of the mapping
      double J = shape_and_test(s, psif, testf);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Need to find position to feed into flux function, initialise to zero
      Vector<double> interpolated_x(Dim, 0.0);

      // Calculate velocities and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over velocity components
        for (unsigned i = 0; i < Dim; i++)
        {
          interpolated_x[i] += nodal_position(l, i) * psif[l];
        }
      }

      // Get the outer unit normal
      Vector<double> interpolated_normal(2);
      outer_unit_normal(ipt, interpolated_normal);

      // Get displacements
      ELASTICITY_BULK_ELEMENT* ext_el_pt =
        dynamic_cast<ELASTICITY_BULK_ELEMENT*>(external_element_pt(0, ipt));
      Vector<double> s_ext(external_element_local_coord(0, ipt));
      Vector<std::complex<double>> displ(2);
      ext_el_pt->interpolated_u_time_harmonic_linear_elasticity(s_ext, displ);

      // Convert into flux BC: This takes the dot product of the
      // actual displacement with the flux element's own outer
      // unit normal so the plus sign is OK.
      std::complex<double> flux =
        (displ[0] * interpolated_normal[0] + displ[1] * interpolated_normal[1]);

      // Now add to the appropriate equations

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        local_eqn = nodal_local_eqn(l, u_index_helmholtz.real());
        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          // Add the prescribed flux terms
          residuals[local_eqn] -= flux.real() * testf[l] * W;

          // Imposed traction doesn't depend upon the solution,
          // --> the Jacobian is always zero, so no Jacobian
          // terms are required
        }

        local_eqn = nodal_local_eqn(l, u_index_helmholtz.imag());
        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          // Add the prescribed flux terms
          residuals[local_eqn] -= flux.imag() * testf[l] * W;

          // Imposed traction doesn't depend upon the solution,
          // --> the Jacobian is always zero, so no Jacobian
          // terms are required
        }
      }
    }
  }
} // namespace oomph

#endif
