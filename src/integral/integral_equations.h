//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#ifndef INTEGRAL_EQUATIONS_HEADER
#define INTEGRAL_EQUATIONS_HEADER

#include "../generic/elements.h"
#include "../generic/Qelements.h"

namespace oomph
{
  template<unsigned DIM>
  class IntegralEquations : public virtual FiniteElement
  {
  public:
    // Typedef the integrand function
    typedef void (*IntegrandFctPt)(const double t,
                                   const Vector<double>& x,
                                   double& f);

  private:
    Vector<IntegrandFctPt> Vector_integrand_fct_pt;

  public:
    // Contructor
    IntegralEquations() {}

    // Destructor
    ~IntegralEquations() {}

    /// Broken copy constructor
    IntegralEquations(const IntegralEquations& dummy) = delete;

    /// Broken assignment operator
    void operator=(const IntegralEquations&) = delete;

    void set_output_data_pt(Data* const& data_pt);

    void setup_integrand(Data* const& output_data_pt,
                         const IntegrandFctPt& integrand_fct_pt = 0);

    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      const bool compute_jacobian = false;
      fill_in_generic_residual_contribution(
        residuals, GeneralisedElement::Dummy_matrix, compute_jacobian);
    }

    inline virtual void get_integrand(const unsigned& n,
                                      const double& t,
                                      const Vector<double>& x,
                                      double& f) const;

    // /// Access function: Pointer to integrand function
    // IntegrandFctPt& integrand_fct_pt()
    // {
    //   return Integrand_fct_pt;
    // }

    // /// Access function: Pointer to integrand function. Const version
    // IntegrandFctPt integrand_fct_pt() const
    // {
    //   return Integrand_fct_pt;
    // }

    /// Get flux: flux[i] = du/dx_i
    void get_integral_flux(const Vector<double>& s, Vector<double>& flux) const;

    virtual void output(std::ostream& outfile);

    virtual void output(std::ostream& outfile, const unsigned& nplot);

  private:
    void fill_in_generic_residual_contribution(Vector<double>& residuals,
                                               DenseMatrix<double>& jacobian,
                                               const unsigned& flag);
  };
} // namespace oomph
#endif
