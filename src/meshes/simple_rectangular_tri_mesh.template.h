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
// Header file for simple 2D triangle mesh class

// Include guards to prevent multiple inclusion of the header
#ifndef OOMPH_SIMPLE_RECTANGULAR_TRIMESH_HEADER
#define OOMPH_SIMPLE_RECTANGULAR_TRIMESH_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// Oomph-lib includes
#include "../generic/mesh.h"
#include "../generic/triangle_mesh.h"

namespace oomph
{
  //===================================================================
  /// Simple 2D triangular mesh for TElements
  //===================================================================
  template<class ELEMENT>
  class SimpleRectangularTriMesh : public virtual TriangleMeshBase
  {
  public:
    /// Constructor
    /// n_x  : number of elements in the x direction;
    /// n_y  : number of elements in the y direction;
    /// l_x  : length in the x direction
    /// l_y  : length in the y direction
    /// Ordering of elements: 'lower left' to 'lower right' then 'upwards'
    SimpleRectangularTriMesh(
      const unsigned& n_x,
      const unsigned& n_y,
      const double& l_x,
      const double& l_y,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Access function for number of elements in x directions
    const unsigned& nx() const
    {
      return Nx;
    }

    /// Access function for number of elements in y directions
    const unsigned& ny() const
    {
      return Ny;
    }

  private:
    /// Number of elements in x direction
    unsigned Nx;

    /// Number of elements in y directions
    unsigned Ny;

    /// Length of mesh in x-direction
    double Lx;

    /// Length of mesh in y-direction
    double Ly;
  };

} // namespace oomph

#endif
