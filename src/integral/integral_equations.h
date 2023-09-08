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
#ifndef INTEGRAL_EQUATIONS_HEADER
#define INTEGRAL_EQUATIONS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// OOMPH-LIB headers
#include "../generic/elements.h"
#include "../generic/Qelements.h"

namespace oomph
{
  //=============================================================
  // A class for all elements that provide the contributions to
  // integrate over a finite element domain.
  // This contains the generic maths. Shape functions, geometric
  // mapping etc. must get implemented in the derived classes.
  //=============================================================
  template<unsigned DIM>
  class IntegralEquations : public virtual FiniteElement
  {
  public:
    // Typedef for the integrand function
    // This can be a function of time and space in general
    typedef void (*IntegrandFctPt)(const double t,
                                   const Vector<double>& x,
                                   double& f);

  private:
    // Provide storage for the (maybe multiple) integrand_fct_pt.
    Vector<IntegrandFctPt> Vector_integrand_fct_pt;

  public:
    // Empty Contructor
    IntegralEquations() {}

    /// Broken copy constructor
    IntegralEquations(const IntegralEquations& dummy) = delete;

    // Empty Destructor
    ~IntegralEquations() {}

    /// Broken assignment operator
    void operator=(const IntegralEquations&) = delete;

    // Set the pointer for the output Data and integrand function.
    // The integrand_fct_pt defaults to not set, which results in a volume
    // integral.
    void setup_integrand(Data* const& output_data_pt,
                         const IntegrandFctPt& integrand_fct_pt = 0);

    // Specific fill in contribution to residuals
    // Calls the generic version which includes the Jacobian.
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      const bool compute_jacobian = false;
      fill_in_generic_residual_contribution(
        residuals, GeneralisedElement::Dummy_matrix, compute_jacobian);
    }

    // Get the value of the n-th integrand at time, t, and location x.
    inline virtual void get_integrand(const unsigned& n,
                                      const double& t,
                                      const Vector<double>& x,
                                      double& f) const;

    // Get the flux of the integrand: flux[i] = du/dx_i at local location s.
    // Useful for the Z2 error estimator refinement
    void get_integral_flux(const Vector<double>& s, Vector<double>& flux) const;

    // Override the FiniteElement output
    virtual void output(std::ostream& outfile);

    // Override the FiniteElement output
    virtual void output(std::ostream& outfile, const unsigned& nplot);

  private:
    // Generic fill in residual contribution
    void fill_in_generic_residual_contribution(Vector<double>& residuals,
                                               DenseMatrix<double>& jacobian,
                                               const unsigned& flag);
  };
} // namespace oomph
#endif
