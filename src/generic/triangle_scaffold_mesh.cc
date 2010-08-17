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
#include "triangle_scaffold_mesh.h"


namespace oomph
{


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//==================================================================
/// Helper namespace for triangle meshes
//==================================================================
namespace TriangleHelper
{
 /// Clear TriangulateIO structure
 void clear_triangulateio(TriangulateIO& triangulate_io,
                          const bool& clear_hole_data)
 {    
  // Clear the point,attribute and marker list 
  free(triangulate_io.pointlist); 
  free(triangulate_io.pointattributelist);
  free(triangulate_io.pointmarkerlist); 
  triangulate_io.numberofpoints = 0;
  triangulate_io.numberofpointattributes = 0;

  // Clear the triangle, attribute,neighbor and area list 
  free(triangulate_io.trianglelist);    
  free(triangulate_io.triangleattributelist);
  free(triangulate_io.trianglearealist);
  free(triangulate_io.neighborlist);
  triangulate_io.numberoftriangles = 0; 
  triangulate_io.numberofcorners = 0;
  triangulate_io.numberoftriangleattributes = 0;

  // Clear the segment and marker list 
  free(triangulate_io.segmentlist);
  free(triangulate_io.segmentmarkerlist);
  triangulate_io.numberofsegments = 0;

  // Clear hole list
  if (clear_hole_data) free(triangulate_io.holelist);
  triangulate_io.numberofholes = 0;

  // Clear region list 
  free(triangulate_io.regionlist);
  triangulate_io.numberofregions = 0;

  // Clear edge, marker and norm list 
  free(triangulate_io.edgelist);
  free(triangulate_io.edgemarkerlist);
  free(triangulate_io.normlist);
  triangulate_io.numberofedges = 0;

  // Now null it all out again 
  initialise_triangulateio(triangulate_io);
 }


 /// Initialise TriangulateIO structure
 void initialise_triangulateio(TriangulateIO& triangle_io)
 {
  // Initialize the point list 
  triangle_io.pointlist = (double *) NULL; 
  triangle_io.pointattributelist = (double *) NULL;  
  triangle_io.pointmarkerlist = (int *) NULL; 
  triangle_io.numberofpoints = 0;
  triangle_io.numberofpointattributes = 0;

  // Initialize the triangle list 
  triangle_io.trianglelist = (int *) NULL;    
  triangle_io.triangleattributelist = (double *) NULL;
  triangle_io.trianglearealist = (double *) NULL;
  triangle_io.neighborlist = (int *) NULL;
  triangle_io.numberoftriangles = 0; 
  triangle_io.numberofcorners = 0;
  triangle_io.numberoftriangleattributes = 0;

  // Initialize the segment list 
  triangle_io.segmentlist = (int *) NULL;
  triangle_io.segmentmarkerlist = (int *) NULL;
  triangle_io.numberofsegments = 0;
  
  // Initialise hole list
  triangle_io.holelist = (double *) NULL;
  triangle_io.numberofholes = 0;

  // Initialize region list 
  triangle_io.regionlist = (double *) NULL;
  triangle_io.numberofregions = 0;

  // Initialize edge list 
  triangle_io.edgelist = (int *) NULL;
  triangle_io.edgemarkerlist = (int *) NULL;
  triangle_io.normlist = (double *) NULL;
  triangle_io.numberofedges = 0;
  

 }


 /// \short Make (partial) deep copy of TriangulateIO object. We only copy
 /// those items we need within oomph-lib's adaptation procedures.
 /// Warnings are issued if triangulate_io contains data that is not
 /// not copied, unless quiet=true;
 TriangulateIO deep_copy_of_triangulateio_representation(
  TriangulateIO& triangle_io, const bool& quiet)
 {
  // Create the struct
  TriangulateIO triangle_out;
  
  // Initialise
  initialise_triangulateio(triangle_out);

  // Point data
  triangle_out.numberofpoints = triangle_io.numberofpoints;
  triangle_out.pointlist = (double *) 
   malloc(triangle_out.numberofpoints * 2 * sizeof(double));
  for (int j=0;j<triangle_out.numberofpoints*2;j++)
   {
    triangle_out.pointlist[j]=triangle_io.pointlist[j];
   }
  
  triangle_out.pointmarkerlist = 
   (int *) malloc(triangle_out.numberofpoints * sizeof(int));
  for (int j=0;j<triangle_out.numberofpoints;j++)
   {
    triangle_out.pointmarkerlist[j]=triangle_io.pointmarkerlist[j];
   }
  
  // Warn about laziness...
  if (!quiet)
   {
    if ((triangle_io.pointattributelist!=0)||
        (triangle_io.numberofpointattributes!=0))
     {
      OomphLibWarning(
       "Point attributes are not currently copied across",
       "TriangleHelper::deep_copy_of_triangulateio_representation",
       OOMPH_EXCEPTION_LOCATION);  
     }
   }


  // Triangle data
  triangle_out.numberoftriangles=triangle_io.numberoftriangles;
  triangle_out.trianglelist = 
   (int *) malloc(triangle_out.numberoftriangles * 3 * sizeof(int)); 
  for (int j=0;j<triangle_out.numberoftriangles*3;j++)
   {
    triangle_out.trianglelist[j]=triangle_io.trianglelist[j];
   }


  
  // Warn about laziness...
  if (!quiet)
   {
    if ((triangle_io.triangleattributelist!=0)||
        (triangle_io.numberoftriangleattributes!=0))
     {
      OomphLibWarning(
       "Triangle attributes are not currently copied across",
       "TriangleHelper::deep_copy_of_triangulateio_representation",
       OOMPH_EXCEPTION_LOCATION);  
     }

    if ((triangle_io.trianglearealist!=0))
     {
      OomphLibWarning(
       "Triangle areas are not currently copied across",
       "TriangleHelper::deep_copy_of_triangulateio_representation",
       OOMPH_EXCEPTION_LOCATION);  
     }
    
    if ((triangle_io.neighborlist!=0))
     {
      OomphLibWarning(
       "Triangle neighbours are not currently copied across",
       "TriangleHelper::deep_copy_of_triangulateio_representation",
       OOMPH_EXCEPTION_LOCATION);  
     }
   }
  
  
  triangle_out.numberofcorners=triangle_io.numberofcorners;

  // Segment data
  triangle_out.numberofsegments=triangle_io.numberofsegments;
  triangle_out.segmentlist = 
   (int *) malloc(triangle_out.numberofsegments * 2 * sizeof(int));
  for (int j=0;j<triangle_out.numberofsegments*2;j++)
   {
    triangle_out.segmentlist[j]=triangle_io.segmentlist[j];
   }
  triangle_out.segmentmarkerlist = 
   (int *) malloc(triangle_out.numberofsegments * sizeof(int));
  for (int j=0;j<triangle_out.numberofsegments;j++)
   {
    triangle_out.segmentmarkerlist[j]=triangle_io.segmentmarkerlist[j];
   }
  
  
  // Hole data
  triangle_out.numberofholes=triangle_io.numberofholes;
  triangle_out.holelist =
   (double*) malloc(triangle_out.numberofholes * 2 * sizeof(double));
  for (int j=0;j<triangle_out.numberofholes*2;j++)
   {
    triangle_out.holelist[j]=triangle_io.holelist[j];
   }
  
  
  // Warn about laziness...
  if (!quiet)
   {
    if ((triangle_io.regionlist!=0)||
        (triangle_io.numberofregions!=0))
     {
      OomphLibWarning(
       "Regions are not currently copied across",
       "TriangleHelper::deep_copy_of_triangulateio_representation",
       OOMPH_EXCEPTION_LOCATION);  
     }
    
    if ((triangle_io.edgelist!=0)||
        (triangle_io.numberofedges!=0))
     {
      OomphLibWarning(
       "Edges are not currently copied across",
       "TriangleHelper::deep_copy_of_triangulateio_representation",
       OOMPH_EXCEPTION_LOCATION);  
     } 

    if ((triangle_io.edgemarkerlist!=0))
     {
      OomphLibWarning(
       "Edge markers are not currently copied across",
       "TriangleHelper::deep_copy_of_triangulateio_representation",
       OOMPH_EXCEPTION_LOCATION);  
     } 
    
    if ((triangle_io.normlist!=0))
     {
      OomphLibWarning(
       "Normals are not currently copied across",
       "TriangleHelper::deep_copy_of_triangulateio_representation",
       OOMPH_EXCEPTION_LOCATION);  
     } 
   }

  // Give it back!
  return triangle_out;
  
 }

}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



//=====================================================================
/// \short Constructor: Pass the filenames of the triangle files
//=====================================================================
 TriangleScaffoldMesh::TriangleScaffoldMesh(const std::string& node_file_name,
                                            const std::string& ele_file_name,
                                            const std::string& poly_file_name)
 {
  
  // Process element file
  //---------------------
  std::ifstream element_file(ele_file_name.c_str(),std::ios_base::in);
   
  // Number of elements
  unsigned n_element;
  element_file>>n_element;
   
  // Number of nodes per element
  unsigned n_local_node;
  element_file>>n_local_node;
  if (n_local_node!=3)
   {
    std::ostringstream error_stream;
    error_stream 
     << "Triangle should only be used to generate 3-noded triangles!\n"
     << "Your triangle input file, contains data for " 
     << n_local_node << "-noded triangles" << std::endl;

    throw OomphLibError(error_stream.str(),
                        "TriangleScaffoldMesh::TriangleScaffoldMesh()",
                        OOMPH_EXCEPTION_LOCATION);
   }
   
  //Dummy nodal attributes
  unsigned dummy_attribute;

  //Element attributes may be used if we have internal boundaries
  Element_attribute.resize(n_element,0.0);

  // Dummy storage for element numbers
  unsigned dummy_element_number;
   
  // Temporary stoorage for global node numbers listed element-by-element
  Vector<unsigned> global_node(n_element*n_local_node);
   
  // Initialise counter
  unsigned k=0;

  // Are attributes specified?
  unsigned attribute_flag;
  element_file>>attribute_flag;

  // Read global node numbers for all elements
  if(attribute_flag==0)
   {
    for(unsigned i=0;i<n_element;i++)
     {
      element_file>>dummy_element_number;
      for(unsigned j=0;j<n_local_node;j++)
       {
        element_file>>global_node[k];
        k++;
       }
     }
   }
  else
   {
    for(unsigned i=0;i<n_element;i++)
     {
      element_file>>dummy_element_number;
      for(unsigned j=0;j<n_local_node;j++)
       {
        element_file>>global_node[k];
        k++;
       }
      element_file>> Element_attribute[i];
     }
   }
  element_file.close();
   
  // Resize the Element vector
  Element_pt.resize(n_element);
   
   
   
   
  // Process node file
  // -----------------
  std::ifstream node_file(node_file_name.c_str(),std::ios_base::in);

  // Read number of nodes
  unsigned n_node;
  node_file>>n_node;
   
  // Create a vector of boolean so as not to create the same node twice
  std::vector<bool> done(n_node);
  for (unsigned i=0;i<n_node;i++){done[i]=false;}
   
  // Resize the Node vector
  Node_pt.resize(n_node);
   
  // Spatial dimension of nodes
  unsigned dimension;
  node_file>>dimension;
   
#ifdef PARANOID
  if(dimension!=2)
   {
    throw OomphLibError("The dimension must be 2\n",
                        "TriangleScaffoldMesh::TriangleScaffoldMesh()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
   
  // Flag for attributes
  node_file>> attribute_flag;
   
  // Flag for boundary markers
  unsigned boundary_markers_flag;
  node_file>>boundary_markers_flag;
   
  // Dummy for node number
  unsigned dummy_node_number;
   
  // Create storage for nodal posititions and boundary markers
  Vector<double> x_node(n_node);
  Vector<double> y_node(n_node); 
  Vector<unsigned> bound(n_node);
   
  // Check if there are attributes
  if (attribute_flag==0)
   {
    if(boundary_markers_flag==1)
     {
      for(unsigned i=0;i<n_node;i++)
       {
        node_file>>dummy_node_number;
        node_file>>x_node[i];
        node_file>>y_node[i];
        node_file>>bound[i];
       }
      node_file.close();
     }
    else
     {
      for(unsigned i=0;i<n_node;i++)
       {
        node_file>>dummy_node_number;
        node_file>>x_node[i];
        node_file>>y_node[i];
        bound[i]=0;
       }
      node_file.close();
     }
   }
  else
   {
    if(boundary_markers_flag==1)
     {
      for(unsigned i=0;i<n_node;i++)
       {
        node_file>>dummy_node_number;
        node_file>>x_node[i];
        node_file>>y_node[i];
        node_file>>dummy_attribute;
        node_file>>bound[i];
       }
      node_file.close();
     }
    else
     {
      for(unsigned i=0;i<n_node;i++)
       {
        node_file>>dummy_node_number;
        node_file>>x_node[i];
        node_file>>y_node[i];
        node_file>>dummy_attribute;
        bound[i]=0;
       }
      node_file.close();
     }
   } //end
   
   
  // Determine highest boundary index
  // --------------------------------
  unsigned d=0;
  if(boundary_markers_flag==1)
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
   
   
  // Process poly file to extract edges
  //-----------------------------------
   
  // Open poly file
  std::ifstream poly_file(poly_file_name.c_str(),std::ios_base::in);

  // Number of nodes in poly file
  unsigned n_node_poly;
  poly_file>>n_node_poly;

  // Dimension
  poly_file>>dimension;

  // Attribute flag
  poly_file>> attribute_flag;

  // Boundary markers flag
  poly_file>>boundary_markers_flag;


  // Ignore node information: Note: No, we can't extract the
  // actual nodes themselves from here!
  unsigned dummy;
  if(n_node_poly>0)
   {
    if (attribute_flag==0)
     {
      if(boundary_markers_flag==1)
       {
        for(unsigned i=0;i<n_node_poly;i++)
         {
          poly_file>>dummy;
          poly_file>>dummy;
          poly_file>>dummy;
          poly_file>>dummy;
         }
       }
      else
       {
        for(unsigned i=0;i<n_node_poly;i++)
         {
          poly_file>>dummy;
          poly_file>>dummy;
          poly_file>>dummy;
         }
       }
     }
    else
     {
      if(boundary_markers_flag==1)
       {
        for(unsigned i=0;i<n_node_poly;i++)
         {
          poly_file>>dummy;
          poly_file>>dummy;
          poly_file>>dummy;
          poly_file>>dummy;
          poly_file>>dummy;
         }
       }
      else
       {
        for(unsigned i=0;i<n_node_poly;i++)
         {
          poly_file>>dummy;
          poly_file>>dummy;
          poly_file>>dummy;
          poly_file>>dummy;
         }
       }
     }
   } 
 
  // Now extract the segment information
  //------------------------------------

  // Number of segments
  unsigned n_segment;
  poly_file>>n_segment;

  // Boundary marker flag
  poly_file>>boundary_markers_flag;

  // Storage for the global node numbers (in the triangle 1-based
  // numbering scheme!) of the first and second node in each segment
  Vector<unsigned> first_node(n_segment);
  Vector<unsigned> second_node(n_segment);

  // Storage for the boundary marker for each segment
  Vector<unsigned> segment_boundary(n_segment);

  // Dummy for global segment number
  unsigned dummy_segment_number;

  // Extract information for each segment
  for(unsigned i=0;i<n_segment;i++)
   {
    poly_file >> dummy_segment_number;
    poly_file >> first_node[i];
    poly_file >> second_node[i];
    poly_file >> segment_boundary[i];
   } 
  
  // Extract hole center information
  unsigned nhole=0;
  poly_file>>nhole;
  if(nhole!=0)
   {
    Hole_centre.resize(nhole);
    
    // Dummy for hole number
    unsigned dummy_hole;
    // Loop over the holes to get centre coords
    for(unsigned ihole=0;ihole<nhole;ihole++)
     {
      Hole_centre[ihole].resize(2);
      
      // Read the centre value
      poly_file >> dummy_hole;
      poly_file >> Hole_centre[ihole][0];
      poly_file >> Hole_centre[ihole][1];
     }
   }
  else
   {
    Hole_centre.resize(0);
   }
  poly_file.close();


  
  // Create the elements
  //--------------------
  
  // Counter for nodes in the vector that lists
  // the global node numbers of the elements' local nodes 
  unsigned counter=0;
  for(unsigned e=0;e<n_element;e++)
   { 
    Element_pt[e]=new TElement<2,2>;
    for(unsigned j=0;j<n_local_node;j++)
     {
      unsigned global_node_number=global_node[counter];
      if(done[global_node_number-1]==false) //... -1 because node number
       // begins at 1 in triangle
       {
        //If we are on the boundary
        if((boundary_markers_flag==1) &&
           (bound[global_node_number-1]>0))
         {
          //Construct a boundary node
          Node_pt[global_node_number-1] = 
           finite_element_pt(e)->construct_boundary_node(j);
          //Add to the boundary node look-up scheme
          add_boundary_node(bound[global_node_number-1]-1,
                            Node_pt[global_node_number-1]);
         }
        //Otherwise make an ordinary node
        else
         {
          Node_pt[global_node_number-1]=finite_element_pt(e)->construct_node(j);
         }
        done[global_node_number-1]=true;
        Node_pt[global_node_number-1]->x(0)=x_node[global_node_number-1];
        Node_pt[global_node_number-1]->x(1)=y_node[global_node_number-1];
       }
      else
       {
        finite_element_pt(e)->node_pt(j) = Node_pt[global_node_number-1];
       }
      counter++;
     }
   }
  
  // Resize the "matrix" that stores the boundary id for each
  // edge in each element.
  Edge_boundary.resize(n_element);
  
  // Storage for the global node numbers (in triangle's 1-based 
  // numbering scheme) for the zero-th, 1st, and 2nd node in each
  // triangle.
  unsigned zeroth_glob_num=0;
  unsigned first_glob_num=0;
  unsigned second_glob_num=0;
  
  // Loop over the elements
  for(unsigned e=0;e<n_element;e++)
   {
    // Each element has three edges
    Edge_boundary[e].resize(3);
    // By default each edge is NOT on a boundary
    for(unsigned i=0;i<3;i++)
     {
      Edge_boundary[e][i]=0;
     }
    
    // Loop over all the nodes in the mesh to find out the
    // global node numbers of the element's three nodes.
    // Only search until all three have been found
    unsigned found=0;
    for(unsigned i=0;i<n_node;i++)
     {
      Node* nod_pt=Node_pt[i];
      
      // Is the 0th node in the element the same as the
      // candidate node in the mesh?
      if (finite_element_pt(e)->node_pt(0)==nod_pt)
       {
        // The global node number of the zero-th node in that element is
        // (in triangle's 1-based numbering!):
        zeroth_glob_num=i+1; 
        found++;
       }
      
      // Is the 1st node in the element the same as the
      // candidate node in the mesh?
      if (finite_element_pt(e)->node_pt(1)==nod_pt)
       {
        // The global node number of the first node in that element is
        // (in triangle's 1-based numbering!):
        first_glob_num=i+1; 
        found++;
       }
      
      
      // Is the 2nd node in the element the same as the
      // candidate node in the mesh?
      if (finite_element_pt(e)->node_pt(2)==nod_pt)
       {
        // The global node number of the second node in that element is
        // (in triangle's 1-based numbering!):
        second_glob_num=i+1; 
        found++;
       }
      
      // We've found three -- bail out...
      if (found==3) break;
     }
    
    //If we haven't found all the nodes then complain
    if(found < 3)
     {
      std::ostringstream error_stream;
      error_stream << "Only found global node numbers of " 
                   << found << " nodes in element " << e << std::endl; 
      throw OomphLibError(error_stream.str(),
                          "TriangleScaffoldMesh::TriangleScaffoldMesh()",
                          OOMPH_EXCEPTION_LOCATION);
     }
    

    // Now we know the global node numbers of the elements' three nodes
    // in triangle's 1-based numbering. 
    
    // Loop over all the boundary segments and check if 
    // the edges in the element coincide with it. If so,
    // copy the boundary id across.
    for(unsigned i=0;i<n_segment;i++)
     {
      // Zero-th edge
      if ( ( (zeroth_glob_num== first_node[i]) || 
             (zeroth_glob_num==second_node[i])   ) &&
           ( ( first_glob_num== first_node[i]) || 
             ( first_glob_num==second_node[i]) )     )
       {
        // Copy boundary id across
        Edge_boundary[e][0]=segment_boundary[i];

        //Add to the boundary node look-up scheme
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[zeroth_glob_num-1]);
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[first_glob_num-1]); 
       }
      
      // First edge
      if ( ( ( first_glob_num== first_node[i]) ||
             ( first_glob_num==second_node[i])   ) &&
           ( (second_glob_num== first_node[i]) || 
             (second_glob_num==second_node[i]) )     )
       {
        // Copy boundary id across
        Edge_boundary[e][1]=segment_boundary[i];

        //Add to the boundary node look-up scheme
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[second_glob_num-1]);
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[first_glob_num-1]);
       }
      // Second edge
      if ( ( (zeroth_glob_num==first_node[i]) ||
             (zeroth_glob_num==second_node[i])  ) && 
           ( (second_glob_num==first_node[i]) ||
             (second_glob_num==second_node[i]) )    )
       {
        // Copy boundary id across
        Edge_boundary[e][2]=segment_boundary[i];

        //Add to the boundary node look-up scheme
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[second_glob_num-1]);
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[zeroth_glob_num-1]);
       }
     }
   }
 } 


//=====================================================================
/// \short Constructor: Pass a data structure obtained from the triangulate
/// function
//=====================================================================
 TriangleScaffoldMesh::TriangleScaffoldMesh(TriangulateIO& triangle_data)
 {

  
  // Number of elements
  unsigned n_element = static_cast<unsigned>(triangle_data.numberoftriangles);
  
  // Number of nodes per element
  unsigned n_local_node = static_cast<unsigned>(triangle_data.numberofcorners);
  if (n_local_node!=3)
   {
    std::ostringstream error_stream;
    error_stream 
     << "Triangle should only be used to generate 3-noded triangles!\n"
     << "Your triangle input file, contains data for " 
     << n_local_node << "-noded triangles" << std::endl;
    
    throw OomphLibError(error_stream.str(),
                        "TriangleScaffoldMesh::TriangleScaffoldMesh()",
                        OOMPH_EXCEPTION_LOCATION);
   }

  //Element attributes may be used if we have internal boundaries
  Element_attribute.resize(n_element,0.0);
  
  // Temporary stoorage for global node numbers listed element-by-element
  Vector<unsigned> global_node(n_element*n_local_node);
  
  // Initialise counter
  unsigned k=0;
  
  // Are attributes specified?
  unsigned attribute_flag = 
   static_cast<unsigned>(triangle_data.numberoftriangleattributes);
  
  // Read global node numbers for all elements
  if(attribute_flag==0)
   {
    for(unsigned i=0;i<n_element;i++)
     {
      for(unsigned j=0;j<n_local_node;j++)
       {
        global_node[k] = static_cast<unsigned>(triangle_data.trianglelist[k]);
        k++;
       }
     }
   }
  else
   {
    for(unsigned i=0;i<n_element;i++)
     {
      for(unsigned j=0;j<n_local_node;j++)
       {
        global_node[k] = static_cast<unsigned>(triangle_data.trianglelist[k]);
        k++;
       }
      Element_attribute[i] = triangle_data.triangleattributelist[i];
     }
   }
  
  // Resize the Element vector
  Element_pt.resize(n_element);
  
  // Read number of nodes
  unsigned n_node = triangle_data.numberofpoints;
  
  // Create a vector of boolean so as not to create the same node twice
  std::vector<bool> done(n_node);
  for (unsigned i=0;i<n_node;i++){done[i]=false;}
  
  // Resize the Node vector
  Node_pt.resize(n_node);
    
  // Flag for boundary markers
  unsigned boundary_markers_flag = 0;
  if(triangle_data.pointmarkerlist!=0) {boundary_markers_flag=1;}
  
  // Create storage for nodal posititions and boundary markers
  Vector<double> x_node(n_node);
  Vector<double> y_node(n_node); 
  Vector<unsigned> bound(n_node);
  
  // We shall ingnore all point attributes
  if(boundary_markers_flag==1)
   {
    for(unsigned i=0;i<n_node;i++)
     {
      x_node[i] = triangle_data.pointlist[2*i];
      y_node[i] = triangle_data.pointlist[2*i+1];
      bound[i] = static_cast<unsigned>(triangle_data.pointmarkerlist[i]);
     }
   }
  else
   {
    for(unsigned i=0;i<n_node;i++)
     {
      x_node[i] = triangle_data.pointlist[2*i];
      y_node[i] = triangle_data.pointlist[2*i+1];
      bound[i]=0;
     }
   }
  
  // Now extract the segment information
  //------------------------------------
  
  // Number of segments
  unsigned n_segment = triangle_data.numberofsegments;
  
  // Storage for the global node numbers (in the triangle 1-based
  // numbering scheme!) of the first and second node in each segment
  Vector<unsigned> first_node(n_segment);
  Vector<unsigned> second_node(n_segment);
 
  // Storage for the boundary marker for each segment
  Vector<unsigned> segment_boundary(n_segment);
  
  // Extract information for each segment
  for(unsigned i=0;i<n_segment;i++)
   {
    first_node[i] = static_cast<unsigned>(triangle_data.segmentlist[2*i]);
    second_node[i] = static_cast<unsigned>(triangle_data.segmentlist[2*i+1]);
    segment_boundary[i] =  
     static_cast<unsigned>(triangle_data.segmentmarkerlist[i]);
   } 

  // Extract hole center information
  unsigned nhole=triangle_data.numberofholes;
  if(nhole!=0)
   {
    Hole_centre.resize(nhole);
    
    // Coords counter
    unsigned count_coords=0;
    
    // Loop over the holes to get centre coords
    for(unsigned ihole=0;ihole<nhole;ihole++)
     {
      Hole_centre[ihole].resize(2);
      
      // Read the centre value
      Hole_centre[ihole][0]=triangle_data.holelist[count_coords];
      Hole_centre[ihole][1]=triangle_data.holelist[count_coords+1];
      
      // Increment coords counter
      count_coords+=2;
     }
   }
  else
   {
    Hole_centre.resize(0);
   }

  // Determine highest boundary index using the segment_boundary_markers
  // (and not the vertex_boundary_markers!) segment marker value may be 
  // greater than the vertex one.
  // --------------------------------
  unsigned d=0;
  if(boundary_markers_flag==1)
   {  
    d=segment_boundary[0];
    for(unsigned i=1;i<n_segment;i++)
     {
      if (segment_boundary[i]>d)
       {
        d=segment_boundary[i];
       }
     }
   }
   
  // Set number of boundaries
  if(d>0)
   {
    set_nboundary(d);
   }

  
  // Create the elements
  //--------------------
  
  // Counter for nodes in the vector that lists
  // the global node numbers of the elements' local nodes 
  unsigned counter=0;
  for(unsigned e=0;e<n_element;e++)
   { 
    Element_pt[e]=new TElement<2,2>;
    for(unsigned j=0;j<n_local_node;j++)
     {
      unsigned global_node_number=global_node[counter];
      if(done[global_node_number-1]==false) //... -1 because node number
       // begins at 1 in triangle
       {
        //If we are on the boundary
        if((boundary_markers_flag==1) &&
           (bound[global_node_number-1]>0))
         {
          //Construct a boundary node
          Node_pt[global_node_number-1] = 
           finite_element_pt(e)->construct_boundary_node(j);
          //Add to the boundary node look-up scheme
          add_boundary_node(bound[global_node_number-1]-1,
                            Node_pt[global_node_number-1]);
         }
        //Otherwise make an ordinary node
        else
         {
          Node_pt[global_node_number-1]=finite_element_pt(e)->construct_node(j);
         }
        done[global_node_number-1]=true;
        Node_pt[global_node_number-1]->x(0)=x_node[global_node_number-1];
        Node_pt[global_node_number-1]->x(1)=y_node[global_node_number-1];
       }
      else
       {
        finite_element_pt(e)->node_pt(j) = Node_pt[global_node_number-1];
       }
      counter++;
     }
   }
  
  // Resize the "matrix" that stores the boundary id for each
  // edge in each element.
  Edge_boundary.resize(n_element);
  
  // Storage for the global node numbers (in triangle's 1-based 
  // numbering scheme) for the zero-th, 1st, and 2nd node in each
  // triangle.
  unsigned zeroth_glob_num=0;
  unsigned first_glob_num=0;
  unsigned second_glob_num=0;
  
  // Loop over the elements
  for(unsigned e=0;e<n_element;e++)
   {
    // Each element has three edges
    Edge_boundary[e].resize(3);
    // By default each edge is NOT on a boundary
    for(unsigned i=0;i<3;i++)
     {
      Edge_boundary[e][i]=0;
     }
    
    // Loop over all the nodes in the mesh to find out the
    // global node numbers of the element's three nodes.
    // Only search until all three have been found
    unsigned found=0;
    for(unsigned i=0;i<n_node;i++)
     {
      Node* nod_pt=Node_pt[i];
      
      // Is the 0th node in the element the same as the
      // candidate node in the mesh?
      if (finite_element_pt(e)->node_pt(0)==nod_pt)
       {
        // The global node number of the zero-th node in that element is
        // (in triangle's 1-based numbering!):
        zeroth_glob_num=i+1; 
        found++;
       }
      
      // Is the 1st node in the element the same as the
      // candidate node in the mesh?
      if (finite_element_pt(e)->node_pt(1)==nod_pt)
       {
        // The global node number of the first node in that element is
        // (in triangle's 1-based numbering!):
        first_glob_num=i+1; 
        found++;
       }
      
      
      // Is the 2nd node in the element the same as the
      // candidate node in the mesh?
      if (finite_element_pt(e)->node_pt(2)==nod_pt)
       {
        // The global node number of the second node in that element is
        // (in triangle's 1-based numbering!):
        second_glob_num=i+1; 
        found++;
       }
      
      // We've found three -- bail out...
      if (found==3) break;
     }
    
    //If we haven't found all the nodes then complain
    if(found < 3)
     {
      std::ostringstream error_stream;
      error_stream << "Only found global node numbers of " 
                   << found << " nodes in element " << e << std::endl; 
      throw OomphLibError(error_stream.str(),
                          "TriangleScaffoldMesh::TriangleScaffoldMesh()",
                          OOMPH_EXCEPTION_LOCATION);
     }
    

    // Now we know the global node numbers of the elements' three nodes
    // in triangle's 1-based numbering. 
    
    // Loop over all the boundary segments and check if 
    // the edges in the element coincide with it. If so,
    // copy the boundary id across.
    for(unsigned i=0;i<n_segment;i++)
     {
      // Zero-th edge
      if ( ( (zeroth_glob_num== first_node[i]) || 
             (zeroth_glob_num==second_node[i])   ) &&
           ( ( first_glob_num== first_node[i]) || 
             ( first_glob_num==second_node[i]) )     )
       {
        // Copy boundary id across
        Edge_boundary[e][0]=segment_boundary[i];

        //Add to the boundary node look-up scheme
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[zeroth_glob_num-1]);
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[first_glob_num-1]);
       }
      // First edge
      if ( ( ( first_glob_num== first_node[i]) ||
             ( first_glob_num==second_node[i])   ) &&
           ( (second_glob_num== first_node[i]) || 
             (second_glob_num==second_node[i]) )     )
       {
        // Copy boundary id across
        Edge_boundary[e][1]=segment_boundary[i];

        //Add to the boundary node look-up scheme
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[second_glob_num-1]);
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[first_glob_num-1]);
       }
      // Second edge
      if ( ( (zeroth_glob_num==first_node[i]) ||
             (zeroth_glob_num==second_node[i])  ) && 
           ( (second_glob_num==first_node[i]) ||
             (second_glob_num==second_node[i]) )    )
       {
        // Copy boundary id across
        Edge_boundary[e][2]=segment_boundary[i];

        //Add to the boundary node look-up scheme
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[second_glob_num-1]);
        add_boundary_node(segment_boundary[i]-1,
                          Node_pt[zeroth_glob_num-1]);
       }
     }
   }
 } 

}
