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

#include "../generic/mesh.h"
#include "../generic/matrices.h"
#include "../generic/gmsh_brick_scaffold_mesh.h"
#include "../generic/brick_mesh.h"


namespace oomph
{
    //=========start_of_gmshquadmesh_class================================
    /// Quadrilateral mesh generator; Uses input from Reader.
    /// See: http://members.shaw.ca/bjoe/
    /// Currently only for four-noded quads -- extension to higher-order
    /// quads should be trivial (see the corresponding classes for
    /// triangular meshes).
    //========================================================================
    template<class ELEMENT>
    class GmshBrickMesh : public virtual BrickMeshBase
    {
    public:
        /// Default Constructor
        GmshBrickMesh()
        {
            // Mesh can only be built with 3D Telements.
            MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);
        }

        /// Constructor with the input files
        explicit GmshBrickMesh(const std::string& filename, bool verbose,
                          TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper )
        {
            // Mesh can only be built with 3D Qelements.
            MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);


            // Build scaffold
            Tmp_mesh_pt = new GmshBrickScaffoldMesh(filename, verbose);
            Boundaries = Tmp_mesh_pt->Boundaries;
            Ordered_boundary_conditions = Tmp_mesh_pt->Ordered_boundary_conditions;

            // Convert mesh from scaffold to actual mesh
            build_from_scaffold(time_stepper_pt);

            // Kill the scaffold
            delete Tmp_mesh_pt;
            Tmp_mesh_pt = nullptr;
        }

        /// Empty destructor
        ~GmshBrickMesh() override = default;


        std::vector<GMSH::Boundary> Boundaries;

        std::unordered_map<int, int> Ordered_boundary_conditions;
    protected:
        /// Temporary scaffold mesh
        GmshBrickScaffoldMesh* Tmp_mesh_pt = nullptr;

        /// Build mesh from scaffold
        void build_from_scaffold(TimeStepper* time_stepper_pt);
    };


    //==============start_mesh=================================================
    /// Gmsh-based mesh upgraded to become a solid mesh. Automatically
    /// enumerates all boundaries.
    //=========================================================================
    template<class ELEMENT>
    class SolidGmshBrickMesh : public virtual GmshBrickMesh<ELEMENT>,
                          public virtual SolidMesh
    {
    public:
        /// Constructor. Boundary coordinates are setup
        /// automatically.
        explicit SolidGmshBrickMesh(const std::string& file_name, bool verbose,
                        TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
                : GmshBrickMesh<ELEMENT>(file_name,verbose, time_stepper_pt )
        {
            // Assign the Lagrangian coordinates
            set_lagrangian_nodal_coordinates();
        }


        /// Empty Destructor
        ~SolidGmshBrickMesh() override = default;
    };


} // namespace oomph

#endif
