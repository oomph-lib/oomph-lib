// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision$
// LIC//
// LIC// $LastChangedDate$
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_GEOMPACK_SCAFFOLD_MESH_HEADER
#define OOMPH_GEOMPACK_SCAFFOLD_MESH_HEADER

#include "mesh.h"
#include "Qelements.h"

namespace oomph
{
  //=====================================================================
  /// \short Mesh that is based on input files generated by the
  /// quadrilateral mesh generator Geompack.
  //=====================================================================
  class GeompackQuadScaffoldMesh : public virtual Mesh
  {
  public:
    /// Empty constructor
    GeompackQuadScaffoldMesh() {}

    /// \short Constructor: Pass the filename of the mesh files
    GeompackQuadScaffoldMesh(const std::string &mesh_file_name,
                             const std::string &curve_file_name);

    /// Broken copy constructor
    GeompackQuadScaffoldMesh(const GeompackQuadScaffoldMesh &)
    {
      BrokenCopy::broken_copy("GeompackQuadScaffoldMesh");
    }

    /// Broken assignment operator
    void operator=(const GeompackQuadScaffoldMesh &)
    {
      BrokenCopy::broken_assign("GeompackQuadScaffoldMesh");
    }

    /// Empty destructor
    ~GeompackQuadScaffoldMesh() {}

  }; // end class

} // namespace oomph

#endif
