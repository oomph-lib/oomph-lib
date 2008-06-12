//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_TRIANGLE_MESH_HEADER
#define OOMPH_TRIANGLE_MESH_HEADER

#include "../generic/triangle_scaffold_mesh.h" 

namespace oomph
{

//============start_of_triangle_class===================================
/// Triangle mesh build with the help of the scaffold mesh coming  
/// from the triangle mesh generator Triangle.
/// http://www.cs.cmu.edu/~quake/triangle.html
//======================================================================
template <class ELEMENT>
class TriangleMesh : public Mesh
{

public:

 /// \short Constructor with the input files
 TriangleMesh(const std::string& node_file_name,
              const std::string& element_file_name,
              const std::string& poly_file_name,
              TimeStepper* time_stepper_pt=
              &Mesh::Default_TimeStepper,
              const bool &use_attributes=false)
              
  {
   // Build scaffold
   Tmp_mesh_pt= new 
    TriangleScaffoldMesh(node_file_name,element_file_name,poly_file_name);
 
   // Convert mesh from scaffold to actual mesh
   build_from_scaffold(time_stepper_pt,use_attributes);

   // Kill the scaffold
   delete Tmp_mesh_pt;
   Tmp_mesh_pt=0;
  }


 /// Broken copy constructor
 TriangleMesh(const TriangleMesh& dummy) 
  { 
   BrokenCopy::broken_copy("TriangleMesh");
  } 
 
 /// Broken assignment operator
 void operator=(const TriangleMesh&) 
  {
   BrokenCopy::broken_assign("TriangleMesh");
  }


 /// Empty destructor 
 ~TriangleMesh() {}

 /// Return the number of regions specified by attributes
 unsigned nregion() {return Region_element_pt.size();}

 /// Return the number of elements in region i
 unsigned nregion_element(const unsigned &i) 
  {return Region_element_pt[i].size();}

 /// Return the attribute associated with region i
 double region_attribute(const unsigned &i)
  {return Region_attribute[i];}

 /// Return the e-th element in the i-th region
 FiniteElement* region_element_pt(const unsigned &i,
                                  const unsigned &e)
  {return Region_element_pt[i][e];}

  private:

 /// Temporary scaffold mesh
 TriangleScaffoldMesh* Tmp_mesh_pt;
 
 /// Build mesh from scaffold
 void build_from_scaffold(TimeStepper* time_stepper_pt,
                          const bool &use_attributes);


 /// Vectors of elements in each region differentiated by attribute
 Vector<Vector<FiniteElement* > > Region_element_pt;

 /// Vector of attributes associated with the elements in each region
 Vector<double> Region_attribute;
};

}

#endif
