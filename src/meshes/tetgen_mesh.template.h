//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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

#ifndef OOMPH_TETGEN_MESH_HEADER
#define OOMPH_TETGEN_MESH_HEADER

#include "../generic/tetgen_scaffold_mesh.h"
#include "../generic/tet_mesh.h"

namespace oomph
{

//=========start of TetgenMesh class======================================
/// \short  Unstructured tet mesh based on output from Tetgen:
/// http://tetgen.berlios.de/
//========================================================================
template <class ELEMENT>
class TetgenMesh : public virtual TetMeshBase
{

public:

 /// \short Constructor with the input files
 TetgenMesh(const std::string& node_file_name,
            const std::string& element_file_name,
            const std::string& face_file_name,
            TimeStepper* time_stepper_pt=
            &Mesh::Default_TimeStepper)
              
  {
   // Build scaffold
   Tmp_mesh_pt= new 
    TetgenScaffoldMesh(node_file_name,element_file_name,face_file_name);
 
   // Convert mesh from scaffold to actual mesh
   build_from_scaffold(time_stepper_pt);

   // Kill the scaffold
   delete Tmp_mesh_pt;
   Tmp_mesh_pt=0;
  }


 /// \short Constructor with the input files. Setting the boolean
 /// flag to true splits "corner" elements, i.e. elements that
 /// that have at least three faces on a domain boundary. The 
 /// relevant elements are split without introducing hanging
 /// nodes so the sons have a "worse" shape than their fathers.
 /// However, this step avoids otherwise-hard-to-diagnose
 /// problems in fluids problems where the application of
 /// boundary conditions at such "corner" elements can
 /// overconstrain the solution. 
 TetgenMesh(const std::string& node_file_name,
            const std::string& element_file_name,
            const std::string& face_file_name,
            const bool& split_corner_elements,
            TimeStepper* time_stepper_pt=
            &Mesh::Default_TimeStepper)
              
  {
   // Build scaffold
   Tmp_mesh_pt= new 
    TetgenScaffoldMesh(node_file_name,element_file_name,face_file_name);
 
   // Convert mesh from scaffold to actual mesh
   build_from_scaffold(time_stepper_pt);

   // Kill the scaffold
   delete Tmp_mesh_pt;
   Tmp_mesh_pt=0;

   // Split corner elements
   if (split_corner_elements)
    {
     split_elements_in_corners();
    }
  }

 /// Empty destructor 
 ~TetgenMesh() {}


 /// \short Setup boundary coordinate on boundary b which is
 /// assumed to be planar. Boundary coordinates are the
 /// x-y coordinates in the plane of that boundary with the
 /// x-axis along the line from the (lexicographically)
 /// "lower left" to the "upper right" node. The y axis
 /// is obtained by taking the cross-product of the positive
 /// x direction with the outer unit normal computed by
 /// the face elements.
 void setup_boundary_coordinates(const unsigned& b)
 {
  // Dummy file
  std::ofstream some_file;
  
  // Don't switch the normal
  bool switch_normal=false;

  setup_boundary_coordinates(b,switch_normal,some_file);
 }
 

 /// \short Setup boundary coordinate on boundary b which is
 /// assumed to be planar. Boundary coordinates are the
 /// x-y coordinates in the plane of that boundary with the
 /// x-axis along the line from the (lexicographically)
 /// "lower left" to the "upper right" node. The y axis
 /// is obtained by taking the cross-product of the positive
 /// x direction with the outer unit normal computed by
 /// the face elements. Doc faces in output file.
 void setup_boundary_coordinates(const unsigned& b,
                                 std::ofstream& outfile)
  {
   // Don't switch the normal
   bool switch_normal=false;
   
   setup_boundary_coordinates(b,switch_normal,outfile);
  }


 /// \short Setup boundary coordinate on boundary b which is
 /// assumed to be planar. Boundary coordinates are the
 /// x-y coordinates in the plane of that boundary with the
 /// x-axis along the line from the (lexicographically)
 /// "lower left" to the "upper right" node. The y axis
 /// is obtained by taking the cross-product of the positive
 /// x direction with the outer unit normal computed by
 /// the face elements (or its negative if switch_normal is set
 /// to true). Doc faces in output file.
 void setup_boundary_coordinates(const unsigned& b,
                                 const bool& switch_normal,
                                 std::ofstream& outfile);

 
 /// \short Snap boundaries specified by the IDs listed in boundary_id to
 /// a quadratric surface, specified in the file 
 /// quadratic_surface_file_name. This is usually used with vmtk-based
 /// meshes for which oomph-lib's xda to poly conversion code produces the files
 /// "quadratic_fsi_boundary.dat" and "quadratic_outer_solid_boundary.dat"
 /// which specify the quadratic FSI boundary (for the fluid and the solid)
 /// and the quadratic representation of the outer boundary of the solid. 
 /// When used with these files, the flag switch_normal should be
 /// set to true when calling the function for the outer boundary of the
 /// solid. The DocInfo object can be used to label optional output
 /// files. (Uses directory and label).
 void snap_to_quadratic_surface(const Vector<unsigned>& boundary_id,
                                const std::string& quadratic_surface_file_name,
                                const bool& switch_normal,
                                DocInfo& doc_info);
 

 /// \short Non-Delaunay split elements that have three faces on a boundary
 /// into sons. Timestepper species timestepper for new nodes; defaults
 /// to to steady timestepper.
 void split_elements_in_corners(TimeStepper* time_stepper_pt=
                                &Mesh::Default_TimeStepper);
 

  private:

 /// Temporary scaffold mesh
 TetgenScaffoldMesh* Tmp_mesh_pt;
 
 /// Build mesh from scaffold
 void build_from_scaffold(TimeStepper* time_stepper_pt);
 
}; //end class

}

#endif
