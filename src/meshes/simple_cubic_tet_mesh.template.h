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
#ifndef OOMPH_SIMPLE_CUBIC_TET_MESH_HEADER
#define OOMPH_SIMPLE_CUBIC_TET_MESH_HEADER

#include "../generic/Telements.h"
#include "../generic/tet_mesh.h"
#include "../generic/simple_cubic_scaffold_tet_mesh.h"

namespace oomph
{
  //===================================================================
  /// MySimple 3D tet mesh for TElements
  //===================================================================
  template<class ELEMENT>
  class SimpleCubicTetMesh : public TetMeshBase
  {
  public:
    /// Constructor: Pass number of element blocks
    /// in the x, y and z directions and the corresponding dimensions.
    /// Timestepper defaults to Steady.
    SimpleCubicTetMesh(
      const unsigned& n_x,
      const unsigned& n_y,
      const unsigned& n_z,
      const double& l_x,
      const double& l_y,
      const double& l_z,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    {
      // Mesh can only be built with 3D Telements.
      MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(3);


      std::ostringstream warn_message;
      warn_message
        << "Note: The SimpleCubicTetMesh() is quite inefficient.\n"
        << "      If your code takes a long time in the constructor\n"
        << "      consider using another tet mesh\n";
      OomphLibWarning(warn_message.str(),
                      "SimpleCubicTetMesh::SimpleCubicTetMesh()",
                      OOMPH_EXCEPTION_LOCATION);
      oomph_info << "Starting mesh construction..." << std::endl;
      double start_t = TimingHelpers::timer();

      // Build scaffold mesh
      Tmp_mesh_pt =
        new SimpleCubicScaffoldTetMesh(n_x, n_y, n_z, l_x, l_y, l_z);

      // Build actual mesh from scaffold mesh
      build_from_scaffold(time_stepper_pt);

      delete Tmp_mesh_pt;

      double end_t = TimingHelpers::timer();
      oomph_info << "...finished mesh construction. Total time [sec] "
                 << end_t - start_t << std::endl;
    }


  private:
    /// Build mesh from scaffold mesh
    void build_from_scaffold(TimeStepper* time_stepper_pt);

    /// Temporary scaffold mesh
    Mesh* Tmp_mesh_pt;
  };

} // namespace oomph

#endif
