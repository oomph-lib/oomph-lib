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
#ifndef DATA_ELEMENT_HEADER
#define DATA_ELEMENT_HEADER

#include "generic.h"
#include "integral_equations.h"

namespace oomph
{
  class InfoElement : public virtual GeneralisedElement
  {
  public:
    // Constructor
    InfoElement(const unsigned n_internal_data);

    ~InfoElement() {}

    /// Broken copy constructor
    InfoElement(const InfoElement& dummy) = delete;

    /// Broken assignment operator
    void operator=(const InfoElement&) = delete;

    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      bool compute_jacobian = false;
      fill_in_generic_residual_contribution(
        residuals, GeneralisedElement::Dummy_matrix, compute_jacobian);
    }

  private:
    void fill_in_generic_residual_contribution(Vector<double>& residuals,
                                               DenseMatrix<double>& jacobian,
                                               const unsigned& flag);
  };

  InfoElement::InfoElement(const unsigned n_internal_data)
    : GeneralisedElement()
  {
    bool use_fd_for_jacobian = true;
    for (unsigned n = 0; n < n_internal_data; n++)
    {
      this->add_internal_data(new Data(1), use_fd_for_jacobian);
    }
  }


  void InfoElement::fill_in_generic_residual_contribution(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    const unsigned& flag)
  {
    unsigned n_data = this->ninternal_data();
    for (unsigned i_data = 0; i_data < n_data; i_data++)
    {
      unsigned n_value = this->internal_data_pt(i_data)->nvalue();
      for (unsigned i_value = 0; i_value < n_value; i_value++)
      {
        int local_eqn = this->internal_local_eqn(i_data, i_value);

        /// The contribution to the data is minus the value as it is moved from
        /// the right hand side of the equation
        residuals[local_eqn] += -this->internal_data_pt(i_data)->value(i_value);
      }
    }
  }
} // namespace oomph
#endif
