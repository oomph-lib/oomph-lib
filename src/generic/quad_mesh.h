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
// Common base class for all QuadMeshes
#ifndef OOMPH_QUAD_MESH_HEADER
#define OOMPH_QUAD_MESH_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
// mpi headers
#include "mpi.h"
#endif

// oomph-lib includes
#include "Vector.h"
#include "nodes.h"
#include "matrices.h"
#include "mesh.h"

namespace oomph
{
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //================================================================
  /// Base class for quad meshes (meshes made of 2D quad elements).
  //================================================================
  class QuadMeshBase : public virtual Mesh
  {
  public:
    /// Constructor (empty)
    QuadMeshBase() {}

    /// Broken copy constructor
    QuadMeshBase(const QuadMeshBase& node) = delete;

    /// Broken assignment operator
    void operator=(const QuadMeshBase&) = delete;

    /// Destructor (empty)
    virtual ~QuadMeshBase() {}

    /// Setup lookup schemes which establish whic elements are located
    /// next to mesh's boundaries (wrapper to suppress doc).
    void setup_boundary_element_info()
    {
      std::ofstream outfile;
      setup_boundary_element_info(outfile);
    }

    /// Setup lookup schemes which establish whic elements are located
    /// next to mesh's boundaries. Doc in outfile (if it's open).
    void setup_boundary_element_info(std::ostream& outfile);
  };

} // namespace oomph

#endif
