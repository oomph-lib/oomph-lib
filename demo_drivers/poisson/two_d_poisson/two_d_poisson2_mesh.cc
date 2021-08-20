// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
// Mesh builder: This builds the mesh with a specific element
// so that the mesh sources don't have to be re-compiled over
// and over again in the driver code... Useful strategy if
// the driver code is under development and uses a complicated
// mesh with long compile times

// Include Poisson elements/equations
#include "poisson.h"

// Include the full source for mesh (this is the automatically generated
// combined header file that contains the templated header and the
// templated sources)
#include "meshes/simple_rectangular_quadmesh.h"

namespace oomph
{
  // Force build of specific mesh
  template class SimpleRectangularQuadMesh<QPoissonElement<2, 3>>;
} // namespace oomph
