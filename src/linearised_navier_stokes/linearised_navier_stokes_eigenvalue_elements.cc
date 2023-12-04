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
// The element-independent guts for imposition of "constant volume"
// constraints in free surface/interface problems.


#include "linearised_navier_stokes_eigenvalue_elements.h"

namespace oomph
{
  //=====================================================================
  /// Fill in the residuals for the volume constraint
  //====================================================================
  void LinearisedNavierStokesEigenfunctionNormalisationElement::
    fill_in_generic_contribution_to_residuals_normalisation(
      Vector<double>& residuals)
  {
    // Note that this element can only be used with our linearised navier
    // stokes elements
    // Read part
    {
      const int local_eqn = this->eigenvalue_local_eqn(0);
      if (local_eqn >= 0)
      {
        residuals[local_eqn] -= (*Normalisation_pt).real();
      }
    }
    // Imaginary part
    {
      const int local_eqn = this->eigenvalue_local_eqn(1);
      if (local_eqn >= 0)
      {
        residuals[local_eqn] -= (*Normalisation_pt).imag();
      }
    }
    // Bifurcation constraint
    const int local_eqn = this->eigenvalue_local_eqn(2);
    if (local_eqn >= 0)
    {
      residuals[local_eqn] += this->eigenvalue(0);
    }
  }

  //===========================================================================
  /// Constructor: Pass pointer to target volume. "Pressure" value that
  /// "traded" for the volume contraint is created internally (as a Data
  /// item with a single pressure value)
  //===========================================================================
  LinearisedNavierStokesEigenfunctionNormalisationElement::
    LinearisedNavierStokesEigenfunctionNormalisationElement(
      std::complex<double>* const& normalisation_pt)
  {
    // Store pointer to normalisation
    Normalisation_pt = normalisation_pt;

    // Create data, add as internal data and record the index
    // (gets deleted automatically in destructor of GeneralisedElement)
    External_or_internal_data_index_of_eigenvalue =
      add_internal_data(new Data(3));

    // ...and stored the "traded pressure" value as first value
    Index_of_eigenvalue = 0;
  }

  //======================================================================
  /// Constructor: Pass pointer to target volume, pointer to Data
  /// item whose value specified by index_of_traded_pressure represents
  /// the "Pressure" value that "traded" for the volume contraint.
  /// The Data is stored as external Data for this element.
  //======================================================================
  /* LinearisedNavierStokesEigenfunctionNormalisationElement::LinearisedNavierStokesEigenfunctionNormalisationElement(
    double* prescribed_volume_pt,
    Data* p_traded_data_pt,
    const unsigned& index_of_traded_pressure)
   {
    // Store pointer to prescribed volume
    Prescribed_volume_pt = prescribed_volume_pt;

    // Add as external data and record the index
    External_or_internal_data_index_of_traded_pressure=
     add_external_data(p_traded_data_pt);

    // Record that it is external data
    Traded_pressure_stored_as_internal_data=false;

    // Record index
    Index_of_traded_pressure_value=index_of_traded_pressure;
    } */


} // namespace oomph
