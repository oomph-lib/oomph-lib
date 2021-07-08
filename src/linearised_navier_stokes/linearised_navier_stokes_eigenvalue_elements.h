// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
// LIC//
// LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
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
// Header file for elements that allow the imposition of a "constant volume"
// constraint in free surface problems.

// Include guards, to prevent multiple includes
#ifndef LINEARISED_NAVIER_STOKES_EIGENVALUE_ELEMENTS_HEADER
#define LINEARISED_NAVIER_STOKES_EIGENVALUE_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// OOMPH-LIB headers
#include "generic/elements.h"

namespace oomph
{
  //==========================================================================
  /// A class that is used to implement the constraint that the eigenfunction
  /// has a particular normalisation. This element stores the two components
  /// of the eigenvalue.
  //=========================================================================
  class LinearisedNavierStokesEigenfunctionNormalisationElement :
    public GeneralisedElement
  {
  private:
    /// Pointer to the desired normalisation
    std::complex<double>* Normalisation_pt;

    /// Storage for the initial index of the eigenvalue
    unsigned External_or_internal_data_index_of_eigenvalue;

    /// Storage for the offset index of the eigenvalue
    unsigned Index_of_eigenvalue;

    /// \short The local eqn number for the traded pressure
    inline int eigenvalue_local_eqn(const unsigned& i)
    {
      return this->internal_local_eqn(
        External_or_internal_data_index_of_eigenvalue, Index_of_eigenvalue + i);
    }

    /// \short Fill in the residuals for the volume constraint
    void fill_in_generic_contribution_to_residuals_normalisation(
      Vector<double>& residuals);

  public:
    /// \short Constructor: Pass pointer to target volume. "Pressure" value that
    /// "traded" for the volume contraint is created internally (as a Data
    /// item with a single pressure value)
    LinearisedNavierStokesEigenfunctionNormalisationElement(
      std::complex<double>* const& normalisation_pt);

    /// \short Constructor: Pass pointer to target volume, pointer to Data
    /// item whose value specified by index_of_traded_pressure represents
    /// the "Pressure" value that "traded" for the volume contraint.
    /// The Data is stored as external Data for this element.
    /*LinearisedNavierStokesEigenfunctionNormalisationElement(double*
       prescribed_volume_pt, Data* p_traded_data_pt, const unsigned&
       index_of_traded_pressure);*/

    /// \short Empty destructor
    ~LinearisedNavierStokesEigenfunctionNormalisationElement() {}

    /// Access to Data that contains the traded pressure
    inline Data* eigenvalue_data_pt()
    {
      return internal_data_pt(External_or_internal_data_index_of_eigenvalue);
    }

    /// Return the traded pressure value
    inline double eigenvalue(const unsigned& i)
    {
      return eigenvalue_data_pt()->value(Index_of_eigenvalue + i);
    }

    /// Return the index of Data object  at which the traded pressure is stored
    inline unsigned index_of_eigenvalue()
    {
      return Index_of_eigenvalue;
    }

    /// \short Fill in the residuals for the volume constraint
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      this->fill_in_generic_contribution_to_residuals_normalisation(residuals);
    }

    /// \short Fill in the residuals and jacobian for the volume constraint
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // One contribution to jacobian; see comment in that function
      this->fill_in_generic_contribution_to_residuals_normalisation(residuals);

      const int local_eqn = this->eigenvalue_local_eqn(2);
      if (local_eqn >= 0)
      {
        const int local_unknown = this->eigenvalue_local_eqn(0);
        jacobian(local_eqn, local_unknown) += 1.0;
      }
    }

    /// \short Fill in the residuals, jacobian and mass matrix for the volume
    /// constraint.
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // No contribution to jacobian or mass matrix; see comment in that
      // function
      this->fill_in_generic_contribution_to_residuals_normalisation(residuals);
    }
  };

} // namespace oomph
#endif
