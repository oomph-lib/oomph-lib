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


// Oomph-lib includes
#include "generic.h"
#include "helmholtz.h"
#include "time_harmonic_linear_elasticity.h"

namespace oomph
{
  //======================================================================
  /// A class for elements that allow the imposition of an applied traction
  /// in the equations of time-harmonic linear elasticity from a Helmholtz
  /// potential (interpreted as a displacement potential for the fluid in a
  /// quasi-steady, linearised FSI problem.)
  /// The geometrical information can be read from the FaceGeometry<ELEMENT>
  /// class and and thus, we can be generic enough without the need to have
  /// a separate equations class.
  //======================================================================
  template<class ELASTICITY_BULK_ELEMENT, class HELMHOLTZ_BULK_ELEMENT>
  class TimeHarmonicLinElastLoadedByHelmholtzPressureBCElement
    : public virtual FaceGeometry<ELASTICITY_BULK_ELEMENT>,
      public virtual FaceElement,
      public virtual ElementWithExternalElement
  {
  protected:
    /// Pointer to the ratio, \f$ Q \f$ , of the stress used to
    /// non-dimensionalise the fluid stresses to the stress used to
    /// non-dimensionalise the solid stresses.
    double* Q_pt;

    /// Static default value for the ratio of stress scales
    /// used in the fluid and solid equations (default is 1.0)
    static double Default_Q_Value;

    /// Index at which the i-th displacement component is stored
    Vector<std::complex<unsigned>>
      U_index_time_harmonic_linear_elasticity_helmholtz_traction;

    /// Helper function that actually calculates the residuals
    // This small level of indirection is required to avoid calling
    // fill_in_contribution_to_residuals in fill_in_contribution_to_jacobian
    // which causes all kinds of pain if overloading later on
    void fill_in_contribution_to_residuals_helmholtz_traction(
      Vector<double>& residuals);

  public:
    /// Constructor, which takes a "bulk" element and the
    /// value of the index and its limit
    TimeHarmonicLinElastLoadedByHelmholtzPressureBCElement(
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
      // element -- the Helmholtz bulk element that provides the pressure
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

    /// Return the ratio of the stress scales used to non-dimensionalise
    /// the fluid and elasticity equations. E.g.
    /// \f$ Q = (\omega a)^2 \rho/E \f$, i.e. the
    /// ratio between the inertial fluid stress and the solid's elastic
    /// modulus E.
    const double& q() const
    {
      return *Q_pt;
    }

    /// Return a pointer the ratio of stress scales used to
    /// non-dimensionalise the fluid and solid equations.
    double*& q_pt()
    {
      return Q_pt;
    }


    /// Output function
    void output(std::ostream& outfile)
    {
      /// Dummy
      unsigned nplot = 0;
      output(outfile, nplot);
    }

    /// Output function: Plot traction etc at Gauss points
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

        // Get bulk element for potential
        HELMHOLTZ_BULK_ELEMENT* ext_el_pt =
          dynamic_cast<HELMHOLTZ_BULK_ELEMENT*>(external_element_pt(0, ipt));
        Vector<double> s_ext(external_element_local_coord(0, ipt));
        std::complex<double> u_helmholtz =
          ext_el_pt->interpolated_u_helmholtz(s_ext);

        // Traction: Pressure is proportional to POSITIVE potential
        ext_el_pt->interpolated_u_helmholtz(s_ext);
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

    /// C_style output function
    void output(FILE* file_pt)
    {
      FaceGeometry<ELASTICITY_BULK_ELEMENT>::output(file_pt);
    }

    /// C-style output function
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FaceGeometry<ELASTICITY_BULK_ELEMENT>::output(file_pt, n_plot);
    }
  };

  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Static default value for the ratio of stress scales
  /// used in the fluid and solid equations (default is 1.0)
  //=================================================================
  template<class ELASTICITY_BULK_ELEMENT, class HELMHOLTZ_BULK_ELEMENT>
  double TimeHarmonicLinElastLoadedByHelmholtzPressureBCElement<
    ELASTICITY_BULK_ELEMENT,
    HELMHOLTZ_BULK_ELEMENT>::Default_Q_Value = 1.0;


  //=====================================================================
  /// Return the residuals
  //=====================================================================
  template<class ELASTICITY_BULK_ELEMENT, class HELMHOLTZ_BULK_ELEMENT>
  void TimeHarmonicLinElastLoadedByHelmholtzPressureBCElement<
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
    std::vector<std::complex<unsigned>> u_nodal_index(n_dim);
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
            "TimeHarmonicLinElastLoadedByPressureElement::fill_in_contribution_"
            "to_residuals()",
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
        ext_el_pt->interpolated_u_helmholtz(s_ext);
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


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //======================================================================
  /// A class for elements that allow the imposition of an
  /// prescribed flux (determined from the normal displacements of an
  /// adjacent linearly elastic solid. Normal derivative for displacement
  /// potential is given by normal displacement of adjacent solid multiplies
  /// by FSI parameter (q = k^2 B/E).
  /// The element geometry is obtained from the FaceGeometry<ELEMENT>
  /// policy class.
  //======================================================================
  template<class HELMHOLTZ_BULK_ELEMENT, class ELASTICITY_BULK_ELEMENT>
  class HelmholtzFluxFromNormalDisplacementBCElement
    : public virtual FaceGeometry<HELMHOLTZ_BULK_ELEMENT>,
      public virtual FaceElement,
      public virtual ElementWithExternalElement
  {
  public:
    /// Constructor, takes the pointer to the "bulk" element and the
    /// face index identifying the face to which the element is attached.
    HelmholtzFluxFromNormalDisplacementBCElement(
      FiniteElement* const& bulk_el_pt, const int& face_index);

    /// Broken copy constructor
    HelmholtzFluxFromNormalDisplacementBCElement(
      const HelmholtzFluxFromNormalDisplacementBCElement& dummy) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const HelmholtzFluxFromNormalDisplacementBCElement&) =
      delete;*/


    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      fill_in_generic_residual_contribution_helmholtz_flux_from_displacement(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }


    /// Add the element's contribution to its residual vector and its
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
    /// Function to compute the shape and test functions and to return
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


    /// Function to compute the shape and test functions and to return
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
    /// Add the element's contribution to its residual vector.
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

  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, and the
  /// face index that identifies the face of the bulk element to which
  /// this face element is to be attached.
  //===========================================================================
  template<class HELMHOLTZ_BULK_ELEMENT, class ELASTICITY_BULK_ELEMENT>
  HelmholtzFluxFromNormalDisplacementBCElement<HELMHOLTZ_BULK_ELEMENT,
                                               ELASTICITY_BULK_ELEMENT>::
    HelmholtzFluxFromNormalDisplacementBCElement(
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

    // Cast to the appropriate HelmholtzEquation so that we can
    // find the index at which the variable is stored
    // We assume that the dimension of the full problem is the same
    // as the dimension of the node, if this is not the case you will have
    // to write custom elements, sorry
    switch (Dim)
    {
        // One dimensional problem
      case 1:
      {
        HelmholtzEquations<1>* eqn_pt =
          dynamic_cast<HelmholtzEquations<1>*>(bulk_el_pt);
        // If the cast has failed die
        if (eqn_pt == 0)
        {
          std::string error_string =
            "Bulk element must inherit from HelmholtzEquations.";
          error_string +=
            "Nodes are one dimensional, but cannot cast the bulk element to\n";
          error_string += "HelmholtzEquations<1>\n.";
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
        HelmholtzEquations<2>* eqn_pt =
          dynamic_cast<HelmholtzEquations<2>*>(bulk_el_pt);
        // If the cast has failed die
        if (eqn_pt == 0)
        {
          std::string error_string =
            "Bulk element must inherit from HelmholtzEquations.";
          error_string +=
            "Nodes are two dimensional, but cannot cast the bulk element to\n";
          error_string += "HelmholtzEquations<2>\n.";
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
        HelmholtzEquations<3>* eqn_pt =
          dynamic_cast<HelmholtzEquations<3>*>(bulk_el_pt);
        // If the cast has failed die
        if (eqn_pt == 0)
        {
          std::string error_string =
            "Bulk element must inherit from HelmholtzEquations.";
          error_string += "Nodes are three dimensional, but cannot cast the "
                          "bulk element to\n";
          error_string += "HelmholtzEquations<3>\n.";
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
  /// Helper function to compute the element's residual vector and
  /// the Jacobian matrix.
  //===========================================================================
  template<class HELMHOLTZ_BULK_ELEMENT, class ELASTICITY_BULK_ELEMENT>
  void HelmholtzFluxFromNormalDisplacementBCElement<HELMHOLTZ_BULK_ELEMENT,
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
