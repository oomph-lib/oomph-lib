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

#ifndef OOMPH_INFO_ELEMENTS_HEADER
#define OOMPH_INFO_ELEMENTS_HEADER

#include "generic.h"

namespace oomph
{

  class InfoElement : public virtual GeneralisedElement
  {
  public:
    Data* new_internal_data_pt(const unsigned& n_value = 1)
    {
      Data* data_pt = new Data(n_value);
      add_internal_data(data_pt);
      return data_pt;
    }

    void output(std::ostream& out)
    {
      // Loop over internal data and add it's value to the output stream
      for (unsigned i = 0; i < this->ninternal_data(); i++)
      {
        Data* data_pt = this->internal_data_pt(i);
        for (unsigned j = 0; j < data_pt->nvalue(); j++)
        {
          out << data_pt->value(j) << " ";
        }
        out << std::endl;
      }
    }

    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      this->fill_in_generic_residual_contribution(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1
      this->fill_in_generic_residual_contribution(residuals, jacobian, 1);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Call the generic routine with the flag set to 1
      this->fill_in_generic_residual_contribution(residuals, jacobian, 1);
    }

    void fill_in_generic_residual_contribution(Vector<double>& residuals,
                                               DenseMatrix<double>& jacobian,
                                               const unsigned& flag)
    {
      // Loop over internal data and subtract it's value to all the residuals
      for (unsigned i = 0; i < this->ninternal_data(); i++)
      {
        Data* data_pt = this->internal_data_pt(i);
        for (unsigned j = 0; j < data_pt->nvalue(); j++)
        {
          const unsigned local_eqn = this->internal_local_eqn(i, j);
          residuals[local_eqn] -= data_pt->value(j);
        }
      }
      if (flag)
      {
        // Loop over internal data and add minus one to all the jacobians
        for (unsigned i = 0; i < this->ninternal_data(); i++)
        {
          Data* data_pt = this->internal_data_pt(i);
          for (unsigned j = 0; j < data_pt->nvalue(); j++)
          {
            const unsigned local_eqn = this->internal_local_eqn(i, j);
            const unsigned local_unknown = local_eqn;
            jacobian(local_eqn, local_unknown) -= 1.0;
          }
        }
      }
    }
  };

}; // namespace oomph

#endif
