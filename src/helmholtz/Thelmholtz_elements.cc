//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
//Non-inline functions for Helmholtz elements
#include "Thelmholtz_elements.h"


namespace oomph
{



/////////////////////////////////////////////////////////////////////////
// THelmholtzElement
/////////////////////////////////////////////////////////////////////////



//======================================================================
// Set the data for the number of Variables at each node, always two
/// (real and imag part) in every case
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
const unsigned THelmholtzElement<DIM,NNODE_1D>::Initial_Nvalue = 2;

//====================================================================
// Force build of templates
//====================================================================
template class THelmholtzElement<1,2>;
template class THelmholtzElement<1,3>;
template class THelmholtzElement<1,4>;

template class THelmholtzElement<2,2>;
template class THelmholtzElement<2,3>;
template class THelmholtzElement<2,4>;

template class THelmholtzElement<3,2>;
template class THelmholtzElement<3,3>;

}
