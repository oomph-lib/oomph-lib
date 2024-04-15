// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// OOMPH-LIB headers
#include "axisym_fluid_traction_elements.h"

namespace oomph
{
  //=======================================================================
  /// Namespace containing the zero traction function for axisymmetric
  /// Navier Stokes traction elements
  //=======================================================================
  namespace AxisymmetricNavierStokesTractionElementHelper
  {
    //=======================================================================
    /// Default load function (zero traction)
    //=======================================================================
    void Zero_traction_fct(const double& time,
                           const Vector<double>& x,
                           const Vector<double>& N,
                           Vector<double>& load)
    {
      unsigned n_dim = load.size();
      for (unsigned i = 0; i < n_dim; i++)
      {
        load[i] = 0.0;
      }
    }

  } // namespace AxisymmetricNavierStokesTractionElementHelper

  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Namespace containing the default Strouhal number of axisymmetric
  /// linearised FSI.
  //=======================================================================
  namespace LinearisedFSIAxisymmetricNStNoSlipBCHelper
  {
    /// Default for fluid Strouhal number
    double Default_strouhal_number = 1.0;

  } // namespace LinearisedFSIAxisymmetricNStNoSlipBCHelper

} // namespace oomph
