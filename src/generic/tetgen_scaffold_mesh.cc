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
#include "mesh.h"
#include "Telements.h"
#include "tetgen_scaffold_mesh.h"

namespace oomph
{

//======================================================================
/// Constructor: Pass the filename of the tetrahedra file
//======================================================================
TetgenScaffoldMesh::TetgenScaffoldMesh(const std::string& node_file_name,
                                       const std::string& element_file_name,
                                       const std::string& face_file_name)
{

 // Read the number of nodes in the input file
 std::ifstream infile(node_file_name.c_str(),std::ios_base::in);
 unsigned  n_node;
 infile>>n_node;
 infile.close();

 // Create a vector of boolean so as not to create the same node twice
 std::vector<bool> done (n_node);
 for(unsigned i=0;i<n_node;i++){ done[i]=false;}
 
 // Resize the Node vector
 Node_pt.resize(n_node);

 // Read the number of elements, the number of local nodes per element,
 // and the global number of the nodes of each element in the input file
 std::ifstream infile2(element_file_name.c_str(),std::ios_base::in);
 unsigned n_element;
 infile2>>n_element;
 unsigned n_local_node;
 infile2>>n_local_node;

 if (n_local_node!=4)
  {
   std::ostringstream error_stream;
   error_stream
    << "Tetgen should only be used to generate 4-noded tetrahedra!\n" 
    << "Your tetgen input file, contains data for " 
    << n_local_node << "-noded tetrahedra" << std::endl;

   throw OomphLibError(error_stream.str(),
                       "TetgenScaffoldMesh::TetgenScaffoldMesh()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // Attributes can be provided but aren't used..
 unsigned attribute;
   
 // Dummy storage for element numbers
 unsigned dummy;

 // Storage for the local global node numbers for the nodes in an element
 Vector<unsigned> node(n_local_node);
   
 // Temporary storage for the global node numbers listed element-by-element
 // hierher unsigned b=n_element*n_local_node;
 Vector<unsigned> global_node(n_element*n_local_node);
   
 unsigned k=0;
 infile2>>attribute;
 if(attribute==0)
  {
   for(unsigned i=0;i<n_element;i++)
    {
     infile2>>dummy;
     for(unsigned j=0;j<n_local_node;j++)
      {
       infile2>>node[j];
       global_node[k]=node[j];
       k++;
      }
    }
  }
 else
  {
   for(unsigned i=0;i<n_element;i++)
    {
     infile2>>dummy;
     for(unsigned j=0;j<n_local_node;j++)
      {
       infile2>>node[j];
       global_node[k]=node[j];
       k++;
      }
     infile2>>dummy;
    }
  }
 infile2.close();
   
 // Resize the Element vector
 Element_pt.resize(n_element);

 // Read the coordinates, the attributes if it exists and the boundary marker 
 // if it exists of each node in the input file
 std::ifstream infile3(node_file_name.c_str(),std::ios_base::in);
 infile3>>n_node;
 unsigned dimension;
 infile3>>dimension;
 unsigned attribute2;
 infile3>>attribute2;
 unsigned bound_markers;
 infile3>>bound_markers;

 // Create storage
 Vector<double> x_node(n_node);
 Vector<double> y_node(n_node);
 Vector<double> z_node(n_node);
 Vector<unsigned> bound(n_node);

#ifdef PARANOID
 if(dimension!=3)
  {
   throw OomphLibError("The dimension must be 3\n",
                       "TetgenScaffoldMesh::TetgenScaffoldMesh()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // Check if there are attributes
 if (attribute2==0)
  {
   if(bound_markers==1)
    {
     for(unsigned i=0;i<n_node;i++)
      {
       infile3>>dummy;
       infile3>>x_node[i];
       infile3>>y_node[i];
       infile3>>z_node[i];
       infile3>>bound[i];
      }
     infile3.close();
    }
   else
    {
     for(unsigned i=0;i<n_node;i++)
      {
       infile3>>dummy;
       infile3>>x_node[i];
       infile3>>y_node[i];
       infile3>>z_node[i];
       bound[i]=0;
      }
     infile3.close();
    }
  }
 else
  { 
   if(bound_markers==1)
    {
     for(unsigned i=0;i<n_node;i++)
      {
       infile3>>dummy;
       infile3>>x_node[i];
       infile3>>y_node[i];
       infile3>>z_node[i];
       infile3>>dummy;
       infile3>>bound[i];
      }
     infile3.close();
    }
   else
    {
     for(unsigned i=0;i<n_node;i++)
      { 
       infile3>>dummy;
       infile3>>x_node[i];
       infile3>>y_node[i];
       infile3>>z_node[i];
       infile3>>dummy;
       bound[i]=0;
      }
     infile3.close();
    }
  }

 // Process face file to extract boundary faces
 //--------------------------------------------
   
 // Open face file
 std::ifstream face_file(face_file_name.c_str(),std::ios_base::in);
  
 // Number of faces in face file
 unsigned n_face;
 face_file>>n_face;

 // Boundary markers flag
 unsigned boundary_markers_flag;
 face_file>>boundary_markers_flag;

 // Storage for the global node numbers (in the tetgen 1-based
 // numbering scheme!) of the first, second and third  node in
 //  each segment
 Vector<unsigned> first_node(n_face);
 Vector<unsigned> second_node(n_face);
 Vector<unsigned> third_node(n_face);

 // Storage for the boundary marker for each face
 Vector<unsigned> face_boundary_id(n_face);

 // Dummy for global face number
 unsigned dummy_face_number;


 // Extract information for each segment
 for(unsigned i=0;i<n_face;i++)
  {
   face_file >> dummy_face_number;
   face_file >> first_node[i];
   face_file >> second_node[i];
   face_file >> third_node[i];
   face_file >> face_boundary_id[i];
  } 
 face_file.close();



 // Determine the number of distinct boundaries
 unsigned d=0;
 if(bound_markers==1)
  { 
   d=bound[0];
   for(unsigned i=1;i<n_node;i++)
    {
     if (bound[i]>d)
      {
       d=bound[i];
      }
    }
  }


 // Set number of boundaries
 if(d>0)
  {
   set_nboundary(d);
  }

 // Create the elements
 unsigned l=0;
 for(unsigned e=0;e<n_element;e++)
  {
   Element_pt[e]=new TElement<3,2>;
   unsigned c;
   c=global_node[l];
   if(done[c-1]==false) //c-1 because node number begins at 1 in triangle
    { 
     //If the node is on a boundary, construct a boundary node
     if((bound_markers==1) && (bound[c-1] > 0))
      {
       //Construct the boundary ndoe
       Node_pt[c-1] = finite_element_pt(e)->construct_boundary_node(3);
       //Add the boundary node
       add_boundary_node(bound[c-1]-1,Node_pt[c-1]);
      }
     else
      {
       Node_pt[c-1]=finite_element_pt(e)->construct_node(3);
      }
     done[c-1]=true;
     Node_pt[c-1]->x(0)=x_node[c-1];
     Node_pt[c-1]->x(1)=y_node[c-1];
     Node_pt[c-1]->x(2)=z_node[c-1];
    }
   else
    {
     finite_element_pt(e)->node_pt(3) = Node_pt[c-1];
    }
   l++;
   
   for(unsigned j=0;j<(n_local_node-1);j++)
    {
     c=global_node[l];
     if(done[c-1]==false) //c-1 because node number begins at 1 in triangle
      { 
       //If we're on a boundary
       if((bound_markers==1) && (bound[c-1] > 0))
        {
         //Construct the boundary node
         Node_pt[c-1] = finite_element_pt(e)->construct_boundary_node(j);
         //Add the boundary node look-up scheme
         add_boundary_node(bound[c-1]-1,Node_pt[c-1]);
        }
       else
        {
         Node_pt[c-1]=finite_element_pt(e)->construct_node(j);
        }
       done[c-1]=true;
       Node_pt[c-1]->x(0)=x_node[c-1];
       Node_pt[c-1]->x(1)=y_node[c-1];
       Node_pt[c-1]->x(2)=z_node[c-1];
      }
     else
      {
       finite_element_pt(e)->node_pt(j) = Node_pt[c-1];
      }
     l++;
    }
  }
   

 // Resize the "matrix" that stores the boundary id for each
 // face in each element.
 Face_boundary.resize(n_element);
  
 // Storage for the global node numbers (in tetgen's 1-based 
 // numbering scheme) for the zero-th, 1st, 2nd and 3rd node 
 // in each tetrahedron.
 unsigned zeroth_glob_num=0;
 unsigned first_glob_num=0;
 unsigned second_glob_num=0;
 unsigned third_glob_num=0;
  
 // Loop over the elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Each element has four faces
   Face_boundary[e].resize(4);
   // By default each face is NOT on a boundary
   for(unsigned i=0;i<4;i++)
    {
     Face_boundary[e][i]=0;
    }
  
   // Loop over all the nodes in the mesh to find out the
   // global node numbers of the element's four nodes.
   // Only search until all four have been found
   unsigned found=0;
   for(unsigned i=0;i<n_node;i++)
    {
     Node* nod_pt=Node_pt[i];
      
     // Is the 0th node in the element the same as the
     // candidate node in the mesh?
     if (finite_element_pt(e)->node_pt(0)==nod_pt)
      {
       // The global node number of the zero-th node in that element is
       // (in tetgen's 1-based numbering!):
       zeroth_glob_num=i+1; 
       found++;
      }
      
     // Is the 1st node in the element the same as the
     // candidate node in the mesh?
     if (finite_element_pt(e)->node_pt(1)==nod_pt)
      {
       // The global node number of the first node in that element is
       // (in tetgen's 1-based numbering!):
       first_glob_num=i+1; 
       found++;
      }
      
      
     // Is the 2nd node in the element the same as the
     // candidate node in the mesh?
     if (finite_element_pt(e)->node_pt(2)==nod_pt)
      {
       // The global node number of the second node in that element is
       // (in tetgen's 1-based numbering!):
       second_glob_num=i+1; 
       found++;
      }

     // Is the 3rd node in the element the same as the
     // candidate node in the mesh?
     if (finite_element_pt(e)->node_pt(3)==nod_pt)
      {
       // The global node number of the third node in that element is
       // (in tetgen's 1-based numbering!):
       third_glob_num=i+1; 
       found++;
      }
      
      
     // We've found four -- bail out...
     if (found==4) break;
    }
    
   //If we haven't found all the nodes the complain
   if(found < 4)
    {
     std::ostringstream error_stream;
     error_stream << "Only found global node numbers of " 
                  << found << " nodes in element " << e << std::endl; 
     throw OomphLibError(error_stream.str(),
                         "TetgenScaffoldMesh::TetgenScaffoldMesh()",
                         OOMPH_EXCEPTION_LOCATION);
    }
    

    
   // Now we know the global node numbers of the elements' four nodes
   // in tetgen's 1-based numbering. 
    
   // Loop over all the boundary faces and check if 
   // the face in the element coincide with it. If so,
   // copy the boundary id across.
   for(unsigned i=0;i<n_face;i++)
    {
     // Zero-th face
     if ( ( (zeroth_glob_num== first_node[i]) || 
            (zeroth_glob_num==second_node[i]) ||
            (zeroth_glob_num== third_node[i])    ) &&
          ( ( first_glob_num== first_node[i]) || 
            ( first_glob_num==second_node[i]) ||
            ( first_glob_num== third_node[i])    ) &&
          ( (second_glob_num== first_node[i]) || 
            (second_glob_num==second_node[i]) ||
            (second_glob_num== third_node[i])    )       )
      {
       // Copy boundary id across
       Face_boundary[e][0]=face_boundary_id[i];
      }
     
     // First face
     if ( ( (zeroth_glob_num== first_node[i]) || 
            (zeroth_glob_num==second_node[i]) ||
            (zeroth_glob_num== third_node[i])    ) &&
          ( ( first_glob_num== first_node[i]) || 
            ( first_glob_num==second_node[i]) ||
            ( first_glob_num== third_node[i])    ) &&
          ( ( third_glob_num== first_node[i]) || 
            ( third_glob_num==second_node[i]) ||
            ( third_glob_num== third_node[i])    )       )
      {
       // Copy boundary id across
       Face_boundary[e][1]=face_boundary_id[i];
      }

     // Second face
     if ( ( (zeroth_glob_num== first_node[i]) || 
            (zeroth_glob_num==second_node[i]) ||
            (zeroth_glob_num== third_node[i])    ) &&
          ( (second_glob_num== first_node[i]) || 
            (second_glob_num==second_node[i]) ||
            (second_glob_num== third_node[i])    ) &&
          ( ( third_glob_num== first_node[i]) || 
            ( third_glob_num==second_node[i]) ||
            ( third_glob_num== third_node[i])    )       )
      {
       // Copy boundary id across
       Face_boundary[e][2]=face_boundary_id[i];
      }

     // Third face
     if ( ( ( first_glob_num== first_node[i]) || 
            ( first_glob_num==second_node[i]) ||
            ( first_glob_num== third_node[i])    ) &&
          ( (second_glob_num== first_node[i]) || 
            (second_glob_num==second_node[i]) ||
            (second_glob_num== third_node[i])    ) &&
          ( ( third_glob_num== first_node[i]) || 
            ( third_glob_num==second_node[i]) ||
            ( third_glob_num== third_node[i])    )       )
      {
       // Copy boundary id across
       Face_boundary[e][3]=face_boundary_id[i];
      }



    } //end for i
  } // end for e


} //end of constructor
   

}


