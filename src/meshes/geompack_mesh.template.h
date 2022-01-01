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
#ifndef OOMPH_GEOMPACK_MESH_HEADER
#define OOMPH_GEOMPACK_MESH_HEADER

#include "../generic/geompack_scaffold_mesh.h"

namespace oomph
{
  //=========start_of_geompackquadmesh_class================================
  /// Quadrilateral mesh generator; Uses input from Geompack++.
  /// See: http://members.shaw.ca/bjoe/
  /// Currently only for four-noded quads -- extension to higher-order
  /// quads should be trivial (see the corresponding classes for
  /// triangular meshes).
  //========================================================================
  template<class ELEMENT>
  class GeompackQuadMesh : public Mesh
  {
  public:
    /// Constructor with the input files
    GeompackQuadMesh(const std::string& mesh_file_name,
                     const std::string& curve_file_name,
                     TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    {
      // Mesh can only be built with four-noded 2D Qelements.
      MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2,
                                                                            2);

      // Build scaffold
      Tmp_mesh_pt =
        new GeompackQuadScaffoldMesh(mesh_file_name, curve_file_name);

      // Convert mesh from scaffold to actual mesh
      build_from_scaffold(time_stepper_pt);

      // Kill the scaffold
      delete Tmp_mesh_pt;
      Tmp_mesh_pt = 0;
    }

    /// Empty destructor
    ~GeompackQuadMesh() {}

  private:
    /// Temporary scaffold mesh
    GeompackQuadScaffoldMesh* Tmp_mesh_pt;

    /// Build mesh from scaffold
    void build_from_scaffold(TimeStepper* time_stepper_pt);
  };

} // namespace oomph

#endif
