// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
// Non-inline functions and static data for spectral poisson elements
#include "spectral_poisson_elements.h"


namespace oomph
{
  //==========================================================================
  // Static value that returns the number of variables at each node, always 1
  //==========================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned QSpectralPoissonElement<DIM, NNODE_1D>::Initial_Nvalue = 1;

  template class QSpectralPoissonElement<1, 2>;
  template class QSpectralPoissonElement<1, 3>;
  template class QSpectralPoissonElement<1, 4>;
  template class QSpectralPoissonElement<1, 5>;
  template class QSpectralPoissonElement<1, 6>;
  template class QSpectralPoissonElement<1, 7>;
  template class QSpectralPoissonElement<1, 8>;
  template class QSpectralPoissonElement<1, 9>;

  template class QSpectralPoissonElement<2, 2>;
  template class QSpectralPoissonElement<2, 3>;
  template class QSpectralPoissonElement<2, 4>;
  template class QSpectralPoissonElement<2, 5>;
  template class QSpectralPoissonElement<2, 6>;
  template class QSpectralPoissonElement<2, 7>;
  template class QSpectralPoissonElement<2, 8>;


  template class QSpectralPoissonElement<3, 2>;
  template class QSpectralPoissonElement<3, 3>;
  template class QSpectralPoissonElement<3, 4>;
  template class QSpectralPoissonElement<3, 5>;
  template class QSpectralPoissonElement<3, 6>;
  template class QSpectralPoissonElement<3, 7>;
  template class QSpectralPoissonElement<3, 8>;

} // namespace oomph
