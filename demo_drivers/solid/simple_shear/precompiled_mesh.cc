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
#include "solid/solid_elements.h"
#include "solid/refineable_solid_elements.h"
#include "meshes/simple_cubic_mesh.h"

template class oomph::SimpleCubicMesh<oomph::QPVDElement<3,3> >;
template class oomph::SimpleCubicMesh<oomph::QPVDElementWithPressure<3> >;
template class 
oomph::SimpleCubicMesh<oomph::QPVDElementWithContinuousPressure<3> >;
template class 
oomph::SimpleCubicMesh<oomph::RefineableQPVDElement<3,3> >;
template class 
oomph::SimpleCubicMesh<oomph::RefineableQPVDElementWithPressure<3> >;
template class 
oomph::SimpleCubicMesh<oomph::RefineableQPVDElementWithContinuousPressure<3> >;
