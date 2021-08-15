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
#ifndef OOMPH_ONE_D_LAGRANGIAN_MESH_HEADER
#define OOMPH_ONE_D_LAGRANGIAN_MESH_HEADER

// OOMPH-LIB headers
#include "../generic/mesh.h"
#include "../generic/geom_objects.h"
#include "../generic/fsi.h"

// Include the header for the one dimensional mesh
#include "one_d_mesh.template.h"


namespace oomph
{
  //=======================================================================
  /// 1D mesh parametrised in terms of a 1D Lagrangian coordinate.
  /// The Eulerian positions of the nodes are determined by the GeomObject.
  //=======================================================================
  template<class ELEMENT>
  class OneDLagrangianMesh : public OneDMesh<ELEMENT>, public SolidMesh
  {
  private:
    /// Undeformed Eulerian shape
    GeomObject* Undef_eulerian_posn_pt;

    /// Set the default gradients of the elements
    void assign_default_element_gradients();

  public:
    /// \short Constructor: Pass number of elements, length,
    /// pointer to GeomObject that defines the undeformed Eulerian position,
    /// and the timestepper -- defaults to (Steady) default timestepper defined
    /// in the Mesh base class
    OneDLagrangianMesh(
      const unsigned& n_element,
      const double& length,
      GeomObject* undef_eulerian_posn_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// \short Constructor: Pass number of elements, xmin, xmax
    /// pointer to GeomObject that defines the undeformed Eulerian position,
    /// and the timestepper -- defaults to (Steady) default timestepper defined
    /// in the Mesh base class
    OneDLagrangianMesh(
      const unsigned& n_element,
      const double& xmin,
      const double& xmax,
      GeomObject* undef_eulerian_posn_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);


    /// Assign the undeformed Eulerian positions to the nodes
    void assign_undeformed_positions();
  };

} // namespace oomph

#endif
