// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision$
// LIC//
// LIC// $LastChangedDate$
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
// Header file for elements that are used to apply prescribed flux
// boundary conditions to the Fourier decomposed Helmholtz equations
#ifndef OOMPH_FOURIER_DECOMPOSED_HELMHOLTZ_FLUX_ELEMENTS_HEADER
#define OOMPH_FOURIER_DECOMPOSED_HELMHOLTZ_FLUX_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include "math.h"
#include <complex>

// oomph-lib includes
#include "generic/Qelements.h"

namespace oomph
{
  //======================================================================
  /// \short A class for elements that allow the imposition of an
  /// applied flux on the boundaries of Fourier decomposed Helmholtz elements.
  /// The element geometry is obtained from the  FaceGeometry<ELEMENT>
  /// policy class.
  //======================================================================
  template<class ELEMENT>
  class FourierDecomposedHelmholtzFluxElement :
    public virtual FaceGeometry<ELEMENT>,
    public virtual FaceElement
  {
  public:
    /// \short Function pointer to the prescribed-flux function fct(x,f(x)) --
    /// x is a Vector and  the flux is a complex

    typedef void (*FourierDecomposedHelmholtzPrescribedFluxFctPt)(
      const Vector<double>& x, std::complex<double>& flux);

    /// \short Constructor, takes the pointer to the "bulk" element and the
    /// index of the face to which the element is attached.
    FourierDecomposedHelmholtzFluxElement(FiniteElement* const& bulk_el_pt,
                                          const int& face_index);

    ///\short  Broken empty constructor
    FourierDecomposedHelmholtzFluxElement()
    {
      throw OomphLibError("Don't call empty constructor for "
                          "FourierDecomposedHelmholtzFluxElement",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    /// Broken copy constructor
    FourierDecomposedHelmholtzFluxElement(
      const FourierDecomposedHelmholtzFluxElement& dummy)
    {
      BrokenCopy::broken_copy("FourierDecomposedHelmholtzFluxElement");
    }

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const FourierDecomposedHelmholtzFluxElement&)
     {
      BrokenCopy::broken_assign("FourierDecomposedHelmholtzFluxElement");
      }*/

    /// Access function for the prescribed-flux function pointer
    FourierDecomposedHelmholtzPrescribedFluxFctPt& flux_fct_pt()
    {
      return Flux_fct_pt;
    }

    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      fill_in_generic_residual_contribution_fourier_decomposed_helmholtz_flux(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                                 DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_fourier_decomposed_helmholtz_flux(
        residuals, jacobian, 1);
    }

    /// Output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(std::ostream& outfile)
    {
      FiniteElement::output(outfile);
    }

    /// \short Output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      FiniteElement::output(outfile, n_plot);
    }

    /// C-style output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(FILE* file_pt)
    {
      FiniteElement::output(file_pt);
    }

    /// \short C-style output function -- forward to broken version in
    /// FiniteElement until somebody decides what exactly they want to plot
    /// here...
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FiniteElement::output(file_pt, n_plot);
    }

    /// \short Return the index at which the unknown value
    /// is stored. (Real/imag part gives real index of real/imag part).
    virtual inline std::complex<unsigned> u_index_fourier_decomposed_helmholtz()
      const
    {
      return std::complex<unsigned>(
        U_index_fourier_decomposed_helmholtz.real(),
        U_index_fourier_decomposed_helmholtz.imag());
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

    /// Function to calculate the prescribed flux at a given spatial
    /// position
    void get_flux(const Vector<double>& x, std::complex<double>& flux)
    {
      // If the function pointer is zero return zero
      if (Flux_fct_pt == 0)
      {
        flux = std::complex<double>(0.0, 0.0);
      }
      // Otherwise call the function
      else
      {
        (*Flux_fct_pt)(x, flux);
      }
    }

    /// \short The index at which the real and imag part of the
    /// unknown is stored at the nodes
    std::complex<unsigned> U_index_fourier_decomposed_helmholtz;

    /// \short Add the element's contribution to its residual vector.
    /// flag=1(or 0): do (or don't) compute the contribution to the
    /// Jacobian as well.
    virtual void fill_in_generic_residual_contribution_fourier_decomposed_helmholtz_flux(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag);

    /// Function pointer to the (global) prescribed-flux function
    FourierDecomposedHelmholtzPrescribedFluxFctPt Flux_fct_pt;
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  template<class ELEMENT>
  FourierDecomposedHelmholtzFluxElement<ELEMENT>::
    FourierDecomposedHelmholtzFluxElement(FiniteElement* const& bulk_el_pt,
                                          const int& face_index) :
    FaceGeometry<ELEMENT>(), FaceElement()
  {
    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index, this);

    // Initialise the prescribed-flux function pointer to zero
    Flux_fct_pt = 0;

    // Initialise index at which real and imag unknowns are stored
    U_index_fourier_decomposed_helmholtz = std::complex<unsigned>(0, 1);

    // Now read out indices from bulk element
    FourierDecomposedHelmholtzEquations* eqn_pt =
      dynamic_cast<FourierDecomposedHelmholtzEquations*>(bulk_el_pt);
    // If the cast has failed die
    if (eqn_pt == 0)
    {
      std::string error_string =
        "Bulk element must inherit from FourierDecomposedHelmholtzEquations.";
      throw OomphLibError(
        error_string, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      // Read the index from the (cast) bulk element
      U_index_fourier_decomposed_helmholtz =
        eqn_pt->u_index_fourier_decomposed_helmholtz();
    }
  }

  //===========================================================================
  /// Compute the element's residual vector and the (zero) Jacobian matrix.
  //===========================================================================
  template<class ELEMENT>
  void FourierDecomposedHelmholtzFluxElement<ELEMENT>::
    fill_in_generic_residual_contribution_fourier_decomposed_helmholtz_flux(
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
    Vector<double> s(1);

    // Integers to hold the local equation and unknown numbers
    int local_eqn_real = 0, local_eqn_imag = 0;

    // Loop over the integration points
    //--------------------------------
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < 1; i++)
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
      Vector<double> interpolated_x(2, 0.0);

      // Calculate coordinates
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over coordinates
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += nodal_position(l, i) * psif[l];
        }
      }

      // first component
      double r = interpolated_x[0];

      // Get the imposed flux
      std::complex<double> flux(0.0, 0.0);
      get_flux(interpolated_x, flux);

      // Now add to the appropriate equations
      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        local_eqn_real =
          nodal_local_eqn(l, U_index_fourier_decomposed_helmholtz.real());

        /*IF it's not a boundary condition*/
        if (local_eqn_real >= 0)
        {
          // Add the prescribed flux terms
          residuals[local_eqn_real] -= flux.real() * testf[l] * r * W;

          // Imposed traction doesn't depend upon the solution,
          // --> the Jacobian is always zero, so no Jacobian
          // terms are required
        }
        local_eqn_imag =
          nodal_local_eqn(l, U_index_fourier_decomposed_helmholtz.imag());

        /*IF it's not a boundary condition*/
        if (local_eqn_imag >= 0)
        {
          // Add the prescribed flux terms
          residuals[local_eqn_imag] -= flux.imag() * testf[l] * r * W;

          // Imposed traction doesn't depend upon the solution,
          // --> the Jacobian is always zero, so no Jacobian
          // terms are required
        }
      }
    }
  }

} // namespace oomph

#endif
