//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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

#include <algorithm>
#include "map_matrix.h"
#include "triangle_mesh.h"


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

#ifdef OOMPH_HAS_TRIANGLE_LIB

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
  if(clear_hole_data) {free(triangulate_io.regionlist);}
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


  //Copy over the triangle attribute data
  triangle_out.numberoftriangleattributes =
   triangle_io.numberoftriangleattributes;
  triangle_out.triangleattributelist =
   (double *) malloc(triangle_out.numberoftriangles * 
                     triangle_out.numberoftriangleattributes * sizeof(double));
  for(int j=0;
      j<(triangle_out.numberoftriangles*
         triangle_out.numberoftriangleattributes);++j)
   {
    triangle_out.triangleattributelist[j] = 
     triangle_io.triangleattributelist[j];
   }

  
  // Warn about laziness...
  if (!quiet)
   {
    /* if ((triangle_io.triangleattributelist!=0)||
        (triangle_io.numberoftriangleattributes!=0))
     {
      OomphLibWarning(
       "Triangle attributes are not currently copied across",
       "TriangleHelper::deep_copy_of_triangulateio_representation",
       OOMPH_EXCEPTION_LOCATION);  
       }*/

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
  

  //Region data
  triangle_out.numberofregions=triangle_io.numberofregions;
  triangle_out.regionlist =
   (double*) malloc(triangle_out.numberofregions * 4 * sizeof(double));
  for(int j=0;j<triangle_out.numberofregions*4;++j)
   {
    triangle_out.regionlist[j] = triangle_io.regionlist[j];
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
    /* if ((triangle_io.regionlist!=0)||
        (triangle_io.numberofregions!=0))
     {
      OomphLibWarning(
       "Regions are not currently copied across",
       "TriangleHelper::deep_copy_of_triangulateio_representation",
       OOMPH_EXCEPTION_LOCATION);  
       }*/
    
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

 /// \short Write the triangulateio data to disk as a poly file,
 /// mainly used for debugging
 void write_triangulateio_to_polyfile(TriangulateIO &triangle_io,
                                      std::ostream &poly_file)
 {
  //Up the precision dramatiacally
  poly_file.precision(20);

  //Output the number of points and their attributes
  //Store the number of attributes
  const int n_attr = triangle_io.numberofpointattributes;
  poly_file << triangle_io.numberofpoints << "  " << 2 << " " 
            << n_attr << " " ;
  //Determine whether there are point markers
  bool point_markers=true;
  if(triangle_io.pointmarkerlist==NULL) {point_markers=false;}
  //Indicate this in the file
  poly_file << point_markers << "\n";
  
  //Now output the point data
  poly_file << "#Points\n";
  for(int n=0;n<triangle_io.numberofpoints;++n)
   {
    //Output the point number and x and y coordinates
    poly_file << n+1 << " " 
              << triangle_io.pointlist[2*n] << " "
              << triangle_io.pointlist[2*n+1] << " ";
    //Output any attributes
    for(int i=0;i<n_attr;++i)
     {
      poly_file << triangle_io.pointattributelist[n_attr*n+i] << " ";
     }
    //Output the boundary marker
    if(point_markers)
     {
      poly_file << triangle_io.pointmarkerlist[n] << " ";
     }
    poly_file << "\n";
   }

  //Now move onto the segments
  poly_file << "#Lines\n";
  poly_file << triangle_io.numberofsegments << " ";
  //Determine whether there are segment markers
  bool seg_markers=true;
  if(triangle_io.segmentmarkerlist==NULL) {seg_markers=false;}
  //Output this info in the file
  poly_file << seg_markers << "\n";

  //Now output the segment data
  for(int n=0;n<triangle_io.numberofsegments;++n)
   {
    poly_file << n+1 << " "  
              << triangle_io.segmentlist[2*n] << " "
              << triangle_io.segmentlist[2*n+1] << " ";
    //If there is a boundary marker output
    if(seg_markers)
     {
      poly_file << triangle_io.segmentmarkerlist[n] << " ";
     }
    poly_file << "\n";
   }

  //Now output the number of holes
  poly_file << "#No holes\n";
  poly_file << triangle_io.numberofholes << "\n";
  //Output the hole data
  for(int h=0;h<triangle_io.numberofholes;++h)
   {
    poly_file << h+1 << " " 
              << triangle_io.holelist[2*h] << " "
              << triangle_io.holelist[2*h+1] << "\n";
   }

  //Now output the number of regions
  poly_file << "#Assignment of attributes to regions\n";
  poly_file << triangle_io.numberofregions << "\n";
  
  //Loop over the regions
  for(int r=0;r<triangle_io.numberofregions;++r)
   {
    poly_file << r+1 << " ";
    for(unsigned i=0;i<4;i++)
     {
      poly_file << triangle_io.regionlist[4*r+i] << " ";
     }
    poly_file << "\n";
   }
 }



 /// Create a triangulateio data file from ele node and poly
 /// files. This is used if the mesh is generated by using Triangle externally.
 /// The triangulateio structure is required to dump the mesh topology for
 /// restarts.
 void create_triangulateio_from_polyfiles(
  const std::string& node_file_name,
  const std::string& element_file_name,
  const std::string& poly_file_name, TriangulateIO &triangle_io,
  bool &use_attributes)
 {
  //Initialise the TriangulateIO data structure
  initialise_triangulateio(triangle_io);
  
  // Process element file
  std::ifstream element_file(element_file_name.c_str(),std::ios_base::in);

  // Check that the file actually opened correctly
  if(!element_file.is_open())
   {
    std::string error_msg("Failed to open element file: ");
    error_msg += "\"" + element_file_name + "\".";
    throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
   
  // Read in the number of elements
  element_file >> triangle_io.numberoftriangles;
  const unsigned n_element = 
   static_cast<unsigned>(triangle_io.numberoftriangles);

  //Read in the number of nodes per element
  element_file >> triangle_io.numberofcorners;
  const unsigned n_local_node =
   static_cast<unsigned>(triangle_io.numberofcorners);

  //Read in the element attributes
  element_file >> triangle_io.numberoftriangleattributes;
  const unsigned n_attributes =
   static_cast<unsigned>(triangle_io.numberoftriangleattributes);

  //Allocate storage in the data structure
  triangle_io.trianglelist =
   (int *) malloc(triangle_io.numberoftriangles *
                     triangle_io.numberofcorners * sizeof(int));

  if(n_attributes > 0)
   {
    triangle_io.triangleattributelist =
     (double *) malloc(triangle_io.numberoftriangles * 
                       triangle_io.numberoftriangleattributes 
                       * sizeof(double));
   }

  //Dummy storage
  int dummy_element_number;

  //Initialise counter
  unsigned counter=0;
  unsigned counter2=0;

  // Read global node numbers for all elements
  for(unsigned e=0;e<n_element;e++)
   {
    element_file >> dummy_element_number;
    for(unsigned j=0;j<n_local_node;j++)
     {
      element_file >> triangle_io.trianglelist[counter];
      ++counter;
     }
    for(unsigned j=0;j<n_attributes;j++)
     {
      element_file >> triangle_io.triangleattributelist[counter2];
      ++counter2;
     }
   }
  //Close the element file
  element_file.close();
   
  // Process node file
  // -----------------
  std::ifstream node_file(node_file_name.c_str(),std::ios_base::in);

  // Check that the file actually opened correctly
  if(!node_file.is_open())
   {
    std::string error_msg("Failed to open node file: ");
    error_msg += "\"" + node_file_name + "\".";
    throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  // Read number of nodes
  node_file >> triangle_io.numberofpoints;
  unsigned n_node = 
   static_cast<unsigned>(triangle_io.numberofpoints);
   
  // Spatial dimension of nodes
  unsigned dimension;
  node_file>>dimension;
   
#ifdef PARANOID
  if(dimension!=2)
   {
    throw OomphLibError("The dimension must be 2\n",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
   
  // Flag for attributes
  node_file >> triangle_io.numberofpointattributes;
  unsigned n_point_attributes = 
   static_cast<unsigned>(triangle_io.numberofpointattributes);
  
  // Flag for boundary markers
  unsigned boundary_markers_flag;
  node_file >> boundary_markers_flag;

  //Allocate storage
  triangle_io.pointlist = 
   (double*) malloc(triangle_io.numberofpoints * 2 * sizeof(double));
  triangle_io.pointattributelist =
   (double*) malloc(triangle_io.numberofpoints * 
                    triangle_io.numberofpointattributes * sizeof(double));
  if(boundary_markers_flag)
   {
    triangle_io.pointmarkerlist =
     (int*) malloc(triangle_io.numberofpoints * sizeof(int));
   }
   
  // Dummy for node number
  unsigned dummy_node_number;
   
  //Reset counter
  counter=0;
  // Load in nodal posititions, point attributes
  // and boundary markers
  for(unsigned i=0;i<n_node;i++)
   {
    node_file>>dummy_node_number;
    node_file>>triangle_io.pointlist[2*i];
    node_file>>triangle_io.pointlist[2*i+1];
    for(unsigned j=0;j<n_point_attributes;++j)
     {
      node_file>>triangle_io.pointattributelist[counter];
      ++counter;
     }
    if(boundary_markers_flag)
     {
      node_file>>triangle_io.pointmarkerlist[i];
     }
   }
  node_file.close();
   
   
  // Process poly file to extract edges
  //-----------------------------------
   
  // Open poly file
  std::ifstream poly_file(poly_file_name.c_str(),std::ios_base::in);

  // Check that the file actually opened correctly
  if(!poly_file.is_open())
   {
    std::string error_msg("Failed to open poly file: ");
    error_msg += "\"" + poly_file_name + "\".";
    throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  // Number of nodes in poly file --- these will be ignore
  unsigned n_node_poly;
  poly_file >> n_node_poly;

  // Dimension
  poly_file >> dimension;

  // Attribute flag
  unsigned attribute_flag;
  poly_file >> attribute_flag;

  // Boundary markers flag
  poly_file >> boundary_markers_flag;


  // Ignore node information: Note: No, we can't extract the
  // actual nodes themselves from here!
  unsigned dummy;
  for(unsigned i=0;i<n_node_poly;i++)
   {
    //Read in (and discard) node number and x and y coordinates
    poly_file>>dummy;
    poly_file>>dummy;
    poly_file>>dummy;
    //read in the attributes
    for(unsigned j=0;j<attribute_flag;++j)
     {
      poly_file >> dummy;
     }
    //read in the boundary marker
    if(boundary_markers_flag==1)
     {
      poly_file>>dummy;
     }
   }
 
  // Now extract the segment information
  //------------------------------------

  // Number of segments
  poly_file>> triangle_io.numberofsegments;
  unsigned n_segment = 
   static_cast<unsigned>(triangle_io.numberofsegments);

  // Boundary marker flag
  poly_file >> boundary_markers_flag;

  //Allocate storage
  triangle_io.segmentlist = 
   (int *) malloc(triangle_io.numberofsegments * 2 * sizeof(int));
  if(boundary_markers_flag)
   {
    triangle_io.segmentmarkerlist =
     (int *) malloc(triangle_io.numberofsegments * sizeof(int));
   }

  // Dummy for global segment number
  unsigned dummy_segment_number;

  // Extract information for each segment
  for(unsigned i=0;i<n_segment;i++)
   {
    poly_file >> dummy_segment_number;
    poly_file >> triangle_io.segmentlist[2*i];
    poly_file >> triangle_io.segmentlist[2*i+1];
    if(boundary_markers_flag)
     {
      poly_file >> triangle_io.segmentmarkerlist[i];
     }
   }
  
  // Extract hole center information
  poly_file >> triangle_io.numberofholes;
  unsigned n_hole = static_cast<unsigned>(triangle_io.numberofholes);

  //Allocate memory
  triangle_io.holelist = 
   (double*) malloc(triangle_io.numberofholes * 2 * sizeof(double));


  // Dummy for hole number
  unsigned dummy_hole;
  // Loop over the holes to get centre coords
  for(unsigned ihole=0;ihole<n_hole;ihole++)
   {
    // Read the centre value
    poly_file >> dummy_hole;
    poly_file >> triangle_io.holelist[2*ihole];
    poly_file >> triangle_io.holelist[2*ihole+1];
   }

  // Extract regions information
  poly_file >> triangle_io.numberofregions;
  unsigned n_regions = static_cast<unsigned>(triangle_io.numberofregions);

  //Allocate memory
  triangle_io.regionlist = 
   (double*) malloc(triangle_io.numberofregions * 4 * sizeof(double));

  // Check for using regions
  if (n_regions > 0)
   {use_attributes = true;}

  // Dummy for regions number
  unsigned dummy_region;

  // Loop over the regions to get their coords
  for(unsigned iregion=0;iregion<n_regions;iregion++)
   {
    // Read the regions coordinates
    poly_file >> dummy_region;
    poly_file >> triangle_io.regionlist[4*iregion];
    poly_file >> triangle_io.regionlist[4*iregion+1];
    poly_file >> triangle_io.regionlist[4*iregion+2];
    triangle_io.regionlist[4*iregion+3] = 0.0;
    
    // Skip the rest of the line, there is no need to read the size of
    // the elements in the region since that value is no longer used
    poly_file.ignore(80,'\n');
    
    // Verify if not using the default region number (zero)
    if (triangle_io.regionlist[4*iregion+2] == 0) 
     {
      std::ostringstream error_message;
      error_message << "Please use another region id different from zero.\n"
                    << "It is internally used as the default region number.\n";
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
 
   }

  poly_file.close();

 }

  


 /// \short Write all the triangulateio data to disk in a dump file
 /// that can then be used to restart simulations
 void dump_triangulateio(TriangulateIO &triangle_io,
                         std::ostream &dump_file)
 {
  //Dump the triangles first
  dump_file << triangle_io.numberoftriangles 
            << " # number of elements in TriangulateIO" << std::endl;

  dump_file << triangle_io.numberofcorners
            << " # number of nodes in each triangle" << std::endl;

  dump_file << triangle_io.numberoftriangleattributes 
            << " # number of triangle attributes" << std::endl;

  //Loop over and dump the triangle information
  const int n_element = triangle_io.numberoftriangles;
  const int n_local_node = triangle_io.numberofcorners;
  const int n_attribute = triangle_io.numberoftriangleattributes;
  unsigned counter=0, counter2=0;
  for(int e=0;e<n_element;++e)
   {
    //Dump the corners
    dump_file << e 
              << " # element number " << std::endl;
    for(int n=0;n<n_local_node;++n)
     {
      dump_file << triangle_io.trianglelist[counter] << std::endl;
      ++counter;
     }
    //Dump the attributes
    dump_file << n_attribute
              << " # number of attributes" << std::endl;
    for(int n=0;n<n_attribute;++n)
     {
      dump_file << triangle_io.triangleattributelist[counter2] << std::endl;
      ++counter2;
     }
   }
  

  //Dump the points (nodes) next
  dump_file << triangle_io.numberofpoints 
            << " # number of points in TriangulateIO" << std::endl; 
  dump_file << triangle_io.numberofpointattributes 
            << " # number of point attributes" << std::endl;
  //Test whether there are point markers
  bool point_marker_flag = true;
  if(triangle_io.pointmarkerlist==NULL) {point_marker_flag=false;}
  dump_file << point_marker_flag
              << " # point marker flag" << std::endl;
  

  //Now output the point data
  const int n_nodes = triangle_io.numberofpoints;
  const int n_point_attributes = triangle_io.numberofpointattributes;
  counter=0; counter2=0;
  for(int n=0;n<n_nodes;++n)
   {
    dump_file << n << " # point number " << std::endl;
    for(int i=0;i<2;++i)
     {
      dump_file << triangle_io.pointlist[counter] << std::endl;
      ++counter;
     }
    dump_file << n_point_attributes << " # number of point attributes " 
              << std::endl;
    //Output any attributes
    for(int i=0;i<n_point_attributes;++i)
     {
      dump_file << 
       triangle_io.pointattributelist[counter2] 
                << std::endl;
      ++counter2;
     }
    dump_file << point_marker_flag << " # point marker flag "
              << std::endl;
    //Output the boundary marker
    if(point_marker_flag)
     {
      dump_file << triangle_io.pointmarkerlist[n] << std::endl;
     }
   }
  
  //Now move onto the segments
  dump_file << triangle_io.numberofsegments 
            << " # Number of segments in TriangulateIO " << std::endl;
  
  //Determine whether there are segment markers
  bool seg_marker_flag=true;
  if(triangle_io.segmentmarkerlist==NULL) {seg_marker_flag=false;}
  //Output this info in the file
  dump_file << seg_marker_flag << " # segment marker flag " << std::endl;

  const int n_segments = triangle_io.numberofsegments;
  counter=0;
  //Now output the segment data
  for(int n=0;n<n_segments;++n)
   {
    dump_file << n << " # segment number " << std::endl;
    for(int i=0;i<2;++i)
     {
      dump_file << triangle_io.segmentlist[counter] << std::endl;
      ++counter;
     }

    //If there is a boundary marker output
    dump_file << seg_marker_flag << " # segment marker flag " << std::endl;
    if(seg_marker_flag)
     {
      dump_file << triangle_io.segmentmarkerlist[n] << std::endl;
     }
   }

  //Now output the number of holes
  dump_file << triangle_io.numberofholes << " # number of holes " << std::endl;
  const int n_hole = triangle_io.numberofholes;
  //Output the hole data
  for(int h=0;h<n_hole;++h)
   {
    dump_file << h << " # hole number " << std::endl;
    dump_file << triangle_io.holelist[2*h] << std::endl;
    dump_file << triangle_io.holelist[2*h+1] << std::endl;
   }

  //Now output the number of regions
  dump_file << triangle_io.numberofregions 
            << " # number of regions " << std::endl;

  const int n_region = triangle_io.numberofregions;
  //Loop over the regions
  counter=0;
  for(int r=0;r<n_region;++r)
   {
    dump_file << r << " # region number " << std::endl;
    for(unsigned i=0;i<4;i++)
     {
      dump_file << triangle_io.regionlist[counter] << std::endl;
      ++counter;
     }
   }
 }

 /// \short Read the triangulateio data from a dump file on 
 /// disk, which can then be used to restart simulations
void read_triangulateio(std::istream &restart_file,TriangulateIO &triangle_io)
{
 //String for reading
 std::string input_string;

 //Initialise the triangulate data structure
 initialise_triangulateio(triangle_io);

 //Read the first line up to termination sign
 getline(restart_file,input_string,'#');
 //Ignore the rest of the line
 restart_file.ignore(80,'\n');
 //Convert the number
 triangle_io.numberoftriangles = atoi(input_string.c_str());

 //Read the next line up to termination sign
 getline(restart_file,input_string,'#');
 //Ignore the rest of the line
 restart_file.ignore(80,'\n');
 //Convert the number
 triangle_io.numberofcorners = atoi(input_string.c_str());
 
//Read the next line up to termination sign
 getline(restart_file,input_string,'#');
 //Ignore the rest of the line
 restart_file.ignore(80,'\n');
 //Convert the number
 triangle_io.numberoftriangleattributes = atoi(input_string.c_str());
 
 //Convert numbers into register variables
 const int n_element = triangle_io.numberoftriangles;
 const int n_local_node = triangle_io.numberofcorners;
 const int n_attribute = triangle_io.numberoftriangleattributes;

 //Allocate storage in the data structure
 triangle_io.trianglelist =
  (int *) malloc(triangle_io.numberoftriangles *
                 triangle_io.numberofcorners * sizeof(int));
 
 if(n_attribute > 0)
  {
   triangle_io.triangleattributelist =
    (double *) malloc(triangle_io.numberoftriangles * 
                      triangle_io.numberoftriangleattributes 
                      * sizeof(double));
  }
 
 //Loop over elements and load in data
 unsigned counter=0, counter2=0;
 for(int e=0;e<n_element;++e)
  {
   //Read the next line and ignore it
   getline(restart_file,input_string);
   for(int n=0;n<n_local_node;++n)
    {
     getline(restart_file,input_string);
     triangle_io.trianglelist[counter] = atoi(input_string.c_str());
     ++counter;
    }
   //Read the attributes
   getline(restart_file,input_string);
   for(int n=0;n<n_attribute;++n)
    {
     getline(restart_file,input_string);
     triangle_io.triangleattributelist[counter2] = atof(input_string.c_str());
     ++counter2;
    }
  }
 

  //Read the points (nodes) next up to termination string
  getline(restart_file,input_string,'#');
  //ignore the rest
  restart_file.ignore(80,'\n');
  triangle_io.numberofpoints = atoi(input_string.c_str());

  //Read the point attributes next up to termination string
  getline(restart_file,input_string,'#');
  //ignore the rest
  restart_file.ignore(80,'\n');
  triangle_io.numberofpointattributes = atoi(input_string.c_str());

  //Read whether there are point markers
  getline(restart_file,input_string,'#');
  //ignore the rest
  restart_file.ignore(80,'\n');
  int point_marker_flag = atoi(input_string.c_str());

  //Allocate storage
  triangle_io.pointlist = 
   (double*) malloc(triangle_io.numberofpoints * 2 * sizeof(double));
  triangle_io.pointattributelist =
   (double*) malloc(triangle_io.numberofpoints * 
                    triangle_io.numberofpointattributes * sizeof(double));
  if(point_marker_flag)
   {
    triangle_io.pointmarkerlist =
     (int*) malloc(triangle_io.numberofpoints * sizeof(int));
   }
  

  //Now read the point data
  const int n_nodes = triangle_io.numberofpoints;
  const int n_point_attributes = triangle_io.numberofpointattributes;
  counter=0; counter2=0;
  for(int n=0;n<n_nodes;++n)
   {
    //Ignore the first line
    getline(restart_file,input_string);
    //Get the positions
    for(int i=0;i<2;++i)
     {
      getline(restart_file,input_string);
      triangle_io.pointlist[counter] = atof(input_string.c_str());
      ++counter;
     }
    
    //Ignore the next line about point attributes
    getline(restart_file,input_string);

    //Read any attributes
    for(int i=0;i<n_point_attributes;++i)
     {
      getline(restart_file,input_string);
      triangle_io.pointattributelist[counter2] =
       atof(input_string.c_str());
      ++counter2;
     }

    //Ignore the next line
    getline(restart_file,input_string);
    //Output the boundary marker
    if(point_marker_flag)
     {
      getline(restart_file,input_string);
      triangle_io.pointmarkerlist[n] = atoi(input_string.c_str());
     }
   }
  
  //Next read the segments
  getline(restart_file,input_string,'#');
  restart_file.ignore(80,'\n');
  triangle_io.numberofsegments  = atoi(input_string.c_str());
  
  //Determine whether there are segment markers
  getline(restart_file,input_string,'#');
  //ignore the rest
  restart_file.ignore(80,'\n');
  int seg_marker_flag = atoi(input_string.c_str());
  
  //Allocate storage
  triangle_io.segmentlist = 
   (int *) malloc(triangle_io.numberofsegments * 2 * sizeof(int));
  if(seg_marker_flag)
   {
    triangle_io.segmentmarkerlist =
     (int *) malloc(triangle_io.numberofsegments * sizeof(int));
   }

  const int n_segments = triangle_io.numberofsegments;
  counter=0;
  //Now output the segment data
  for(int n=0;n<n_segments;++n)
   {
    getline(restart_file,input_string);
    //get input
    for(int i=0;i<2;++i)
     {
      getline(restart_file,input_string);
      triangle_io.segmentlist[counter] = atoi(input_string.c_str());
      ++counter;
     }

    //Ignore the next line
    getline(restart_file,input_string);
    if(seg_marker_flag)
     {
      getline(restart_file,input_string);
      triangle_io.segmentmarkerlist[n] = atoi(input_string.c_str());
     }
   }
  
  //Now read the holes
  getline(restart_file,input_string,'#');
  restart_file.ignore(80,'\n');
  triangle_io.numberofholes = atoi(input_string.c_str());
  
  //Allocate memory
  triangle_io.holelist = 
   (double*) malloc(triangle_io.numberofholes * 2 * sizeof(double));
  
  const int n_hole = triangle_io.numberofholes;
  //Output the hole data
  for(int h=0;h<n_hole;++h)
   {
    //Ignore the first line
    getline(restart_file,input_string);
    //get the centre
    getline(restart_file,input_string);
    triangle_io.holelist[2*h] = atof(input_string.c_str());
    getline(restart_file,input_string);
    triangle_io.holelist[2*h+1] = atof(input_string.c_str());
   }

  //Now read the number of regions
  getline(restart_file,input_string,'#');
  restart_file.ignore(80,'\n');
  triangle_io.numberofregions = atoi(input_string.c_str());

  const int n_region = triangle_io.numberofregions;

  //Allocate memory
  triangle_io.regionlist = 
   (double*) malloc(triangle_io.numberofregions * 4 * sizeof(double));

  //Loop over the regions
  counter=0;
  for(int r=0;r<n_region;++r)
   {
    getline(restart_file,input_string);
    for(unsigned i=0;i<4;i++)
     {
      getline(restart_file,input_string);
      triangle_io.regionlist[counter] = atof(input_string.c_str());
      ++counter;
     }
   }
}
                         
#endif

} //End of TriangleHelper namespace







////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////






/// Namespace that allows the specification of a tolerance 
/// between vertices at the ends of polylines that are supposed
/// to be at the same position.
namespace ToleranceForVertexMismatchInPolygons
{
 
 /// Acceptable discrepancy for mismatch in vertex coordinates.
 /// In paranoid mode, the code will die if the beginning/end of
 /// two adjacent polylines differ by more than that. If the
 /// discrepancy is smaller (but nonzero) one of the vertex coordinates
 /// get adjusted to match perfectly; without paranoia the vertex
 /// coordinates are taken as they come...
 double Tolerable_error=1.0e-14;

}





////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////





// =======================================================================
// \short Connects the initial vertex of the curve section to a desired
/// target polyline by specifying the vertex number. There is a checking
/// which verifies that the initial vertex is close enough to the
/// destination vertex on the target polyline by no more than the specified
/// tolerance
// =======================================================================
void TriangleMeshCurveSection::connect_initial_vertex_to_polyline(
  TriangleMeshPolyLine *polyline_pt,
  const unsigned &vertex_number,
  const double &tolerance_for_connection)
{

#ifdef PARANOID
 unsigned n_vertices = polyline_pt->nvertex();

 if (n_vertices <= vertex_number)
  {
   std::ostringstream error_stream;
   error_stream
   << "The vertex number you provided (" << vertex_number
   << ") is greater\n than the number of vertices ("
   << n_vertices << "in the specified TriangleMeshPolyLine.\n"
   << "Remember that the vertex index starts at 0" << std::endl
   << "Source boundary (" << boundary_id() << ") wants to connect "
   << "to destination boundary (" << polyline_pt->boundary_id()
   << ")" << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
  }

 // Verify if there is really a match point in the specified
 // connection values
 Vector<double> v_src(2);
 Vector<double> v_dst(2);

 this->initial_vertex_coordinate(v_src);
 v_dst = polyline_pt->vertex_coordinate(vertex_number);

 double error = sqrt((v_src[0] - v_dst[0])*(v_src[0] - v_dst[0]) +
   (v_src[1] - v_dst[1])*(v_src[1] - v_dst[1]));

 if (error > tolerance_for_connection)
  {
   std::ostringstream error_stream;
   error_stream
   << "The associated vertices for the connection"
   << "\nare not close enough. Their respective values are:\n"
   << "Source boundary id:(" << this->boundary_id() << ")\n"
   << "Source vertex x:(" << v_src[0] << ") y:(" << v_src[1] << ")\n"
   << "Destination boundary id:(" << polyline_pt->boundary_id() << ")"
   << "\nAssociated vertex x:(" << v_dst[0] << ") y:(" << v_dst[1] << ")"
   << "\nThe corresponding distance is: ("<<  error << ") but the "
   << "allowed\ntolerance is: (" << tolerance_for_connection <<")"
   << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
  }

#endif

 Initial_vertex_connected = true;
 Initial_vertex_connected_bnd_id = polyline_pt->boundary_id();
 Initial_vertex_connected_n_vertex = vertex_number;
 Initial_vertex_connected_n_chunk = polyline_pt->boundary_chunk();

}

// =======================================================================
// \short Connects the final vertex of the curve section to a desired
/// target polyline by specifying the vertex number. There is a checking
/// which verifies that the final vertex is close enough to the
/// destination vertex on the target polyline by no more than the specified
/// tolerance
// =======================================================================
void TriangleMeshCurveSection::connect_final_vertex_to_polyline(
  TriangleMeshPolyLine *polyline_pt,
  const unsigned &vertex_number,
  const double &tolerance_for_connection)
{

#ifdef PARANOID
 unsigned n_vertices = polyline_pt->nvertex();

 if (n_vertices <= vertex_number)
  {
   std::ostringstream error_stream;
   error_stream
   << "The vertex number you provided (" << vertex_number
   << ") is greater\n than the number of vertices ("
   << n_vertices << "in the specified TriangleMeshPolyLine.\n"
   << "Remember that the vertex index starts at 0" << std::endl
   << "Source boundary (" << boundary_id() << ") wants to connect "
   << "to destination boundary (" << polyline_pt->boundary_id()
   << ")" << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
  }

 // Verify if there is really a match point in the specified
 // connection values
 Vector<double> v_src(2);
 Vector<double> v_dst(2);

 this->final_vertex_coordinate(v_src);
 v_dst = polyline_pt->vertex_coordinate(vertex_number);

 double error = sqrt((v_src[0] - v_dst[0])*(v_src[0] - v_dst[0]) +
   (v_src[1] - v_dst[1])*(v_src[1] - v_dst[1]));

 if (error > tolerance_for_connection)
  {
   std::ostringstream error_stream;
   error_stream
   << "The associated vertices for the connection"
   << "\nare not close enough. Their respective values are:\n"
   << "Source boundary id:(" << this->boundary_id() << ")\n"
   << "Source vertex x:(" << v_src[0] << ") y:(" << v_src[1] << ")\n"
   << "Destination boundary id:(" << polyline_pt->boundary_id() << ")"
   << "\nAssociated vertex x:(" << v_dst[0] << ") y:(" << v_dst[1] << ")"
   << "\nThe corresponding distance is: ("<<  error << ") but the "
   << "allowed\ntolerance is: (" << tolerance_for_connection <<")"
   << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
  }

#endif

 Final_vertex_connected = true;
 Final_vertex_connected_bnd_id = polyline_pt->boundary_id();
 Final_vertex_connected_n_vertex = vertex_number;
 Final_vertex_connected_n_chunk = polyline_pt->boundary_chunk();

}

// =======================================================================
// \short Connects the initial vertex of the curve section to a desired
/// target curviline by specifying the s value (intrinsic value on the
/// geometric object of the curviline) where to connect on the target
/// curviline. There is a checking which verifies that the initial vertex
/// and the coordinates on the given s value are close enough by no more
/// than the given tolerance
// =======================================================================
void TriangleMeshCurveSection::connect_initial_vertex_to_curviline(
  TriangleMeshCurviLine *curviline_pt,
  const double &s_value,
  const double &tolerance_for_connection)
{

#ifdef PARANOID
 double z_initial = curviline_pt->zeta_start();
 double z_final = curviline_pt->zeta_end();
 double z_max = std::max(z_initial,z_final);
 double z_min = std::min(z_initial,z_final);
 if (s_value < z_min || z_max < s_value)
  {
   std::ostringstream error_stream;
   error_stream
   << "The s value you provided for connection (" << s_value
   << ") is out\nof the limits of the specified "
   << "TriangleMeshCurviLine.\nThe limits are [" << z_initial
   << ", " << z_final << "]" << std::endl
   << "Source boundary (" << boundary_id() << ") wants to connect "
   << "to destination boundary (" << curviline_pt->boundary_id()
   << ")" << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
  }

 // Verify if there is really a match point in the specified
 // connection values
 Vector<double> v_src(2);
 Vector<double> v_dst(2);
 Vector<double> z(1);

 this->initial_vertex_coordinate(v_src);
 z[0] = s_value;
 curviline_pt->geom_object_pt()->position(z, v_dst);
 double error = sqrt((v_src[0] - v_dst[0])*(v_src[0] - v_dst[0]) +
   (v_src[1] - v_dst[1])*(v_src[1] - v_dst[1]));
 if (error >= tolerance_for_connection)
  {
   std::ostringstream error_stream;
   error_stream
   << "The associated vertex for the provided connection s value\n"
   << "are not close enough. Their respective values are:\n"
   << "Source boundary id:(" << this->boundary_id() << ")\n"
   << "Source vertex x:(" << v_src[0] << ") y:(" << v_src[1] << ")\n"
   << "Destination boundary id:(" << curviline_pt->boundary_id() << ")"
   << "\nDestination s value (" << s_value << ")\n"
   << "Associated vertex x:(" << v_dst[0] << ") y:(" << v_dst[1] << ")"
   << "\nThe corresponding distance is: ("<<  error << ") but the "
   << "allowed\ntolerance is: (" << tolerance_for_connection <<")"
   << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
  }

#endif

 Initial_vertex_connected_to_curviline = true;
 Initial_vertex_connected = true;
 Initial_vertex_connected_bnd_id = curviline_pt->boundary_id();
 Initial_vertex_connected_n_chunk = curviline_pt->boundary_chunk();

 // We are still not able to compute the vertex number but we can
 // at least store the corresponding s value
 // The corresponding vertex will be computed when the curviline be
 // changed into a polyline
 Initial_s_connection_value = s_value;
 Tolerance_for_s_connection = tolerance_for_connection;

 curviline_pt->add_connection_point(s_value, tolerance_for_connection);

}

// =======================================================================
// \short Connects the final vertex of the curve section to a desired
/// target curviline by specifying the s value (intrinsic value on the
/// geometric object of the curviline) where to connect on the target
/// curviline. There is a checking which verifies that the final vertex
/// and the coordinates on the given s value are close enough by no more
/// than the given tolerance
// =======================================================================
void TriangleMeshCurveSection::connect_final_vertex_to_curviline(
  TriangleMeshCurviLine *curviline_pt,
  const double &s_value,
  const double &tolerance_for_connection)
{

#ifdef PARANOID
 double z_initial = curviline_pt->zeta_start();
 double z_final = curviline_pt->zeta_end();
 double z_max = std::max(z_initial,z_final);
 double z_min = std::min(z_initial,z_final);
 if (s_value < z_min || z_max < s_value)
  {
   std::ostringstream error_stream;
   error_stream
   << "The s value you provided for connection (" << s_value
   << ") is out\nof the limits of the specified "
   << "TriangleMeshCurviLine.\nThe limits are [" << z_initial
   << ", " << z_final << "]" << std::endl
   << "Source boundary (" << boundary_id() << ") wants to connect "
   << "to destination boundary (" << curviline_pt->boundary_id()
   << ")" << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
  }

 // Verify if there is really a match point in the specified
 // connection values
 Vector<double> v_src(2);
 Vector<double> v_dst(2);
 Vector<double> z(1);

 this->final_vertex_coordinate(v_src);
 z[0] = s_value;
 curviline_pt->geom_object_pt()->position(z, v_dst);

 double error = sqrt((v_src[0] - v_dst[0])*(v_src[0] - v_dst[0]) +
   (v_src[1] - v_dst[1])*(v_src[1] - v_dst[1]));

 if (error >= tolerance_for_connection)
  {
   std::ostringstream error_stream;
   error_stream
   << "The associated vertex for the provided connection s value\n"
   << "are not close enough. Their respective values are:\n"
   << "Source boundary id:(" << this->boundary_id() << ")\n"
   << "Source vertex x:(" << v_src[0] << ") y:(" << v_src[1] << ")\n"
   << "Destination boundary id:(" << curviline_pt->boundary_id() << ")"
   << "\nDestination s value (" << s_value << ")\n"
   << "Associated vertex x:(" << v_dst[0] << ") y:(" << v_dst[1] << ")"
   << "\nThe corresponding distance is: ("<<  error << ") but the "
   << "allowed\ntolerance is: (" << tolerance_for_connection <<")"
   << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
  }

#endif

 Final_vertex_connected_to_curviline = true;
 Final_vertex_connected = true;
 Final_vertex_connected_bnd_id = curviline_pt->boundary_id();
 Final_vertex_connected_n_chunk = curviline_pt->boundary_chunk();

 // We are still not able to compute the vertex number but we can
 // at least store the corresponding s value.
 // The corresponding vertex will be computed when the curviline be
 // transformed into a polyline
 Final_s_connection_value = s_value;
 Tolerance_for_s_connection = tolerance_for_connection;

 curviline_pt->add_connection_point(s_value, tolerance_for_connection);

}







///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////







//=====================================================================
/// Class defining a closed curve for the Triangle mesh generation
//=====================================================================
TriangleMeshClosedCurve::TriangleMeshClosedCurve(
  const Vector<TriangleMeshCurveSection*> &curve_section_pt,
  const Vector<double>& internal_point_pt) :
  TriangleMeshCurve(curve_section_pt),
  Internal_point_pt(internal_point_pt)
{
 
#ifdef PARANOID
 
 // Matching of curve sections i.e. the last vertex of the i curve
 // section should match with the first vertex of the i+1 curve
 // section
 
 // Total number of boundaries
 const unsigned n_boundaries = Curve_section_pt.size();
 
 // Need at least two
 if (n_boundaries<2)
  {
   std::ostringstream error_stream;
   error_stream
    << "Sorry -- I'm afraid we insist that a closed curve is\n"
    << "specified by at least two separate CurveSections.\n"
    << "You've only specified (" << n_boundaries << ")" << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Check last point of each boundary bit coincides with first point
 // on next one
 for (unsigned i=0;i<n_boundaries-1;i++)
  {
   
   // Auxiliary vertex for storing the vertex values of contiguous curves
   Vector<double> v1(2);
   
   // This is for getting the final coordinates of the i curve section
   curve_section_pt[i]->final_vertex_coordinate(v1);
   
   // Auxiliary vertex for storing the vertex values of contiguous curves
   Vector<double> v2(2);
   
   // This is for the start coordinates of the i+1 curve section
   curve_section_pt[i+1]->initial_vertex_coordinate(v2);
   
   // Work out error
   double error = sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2));
   
   if (error>ToleranceForVertexMismatchInPolygons::Tolerable_error)
    {
     std::ostringstream error_stream;
     error_stream
     << "The start and end points of curve section boundary parts\n" << i
     << " and " <<i+1<< " don't match when judged with the tolerance of "
     << ToleranceForVertexMismatchInPolygons::Tolerable_error
     << " which\nis specified in the namespace variable\n"
     << "ToleranceForVertexMismatchInPolygons::Tolerable_error.\n\n"
     << "These are the vertices coordinates:\n"
     << "Curve section ("<<i<<") final vertex: ("
     << v1[0] <<", "<< v1[1] <<")\n"
     << "Curve section ("<<i+1<<") initial vertex: ("
     << v2[0] <<", "<< v2[1] <<")\n"
     << "The distance between the vertices is ("<<error<< ")\n"
     << "Feel free to adjust this or to recompile the code without\n"
     << "paranoia if you think this is OK...\n"
     << std::endl;
     throw OomphLibError(
       error_stream.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
    }
   else
    {
     // Aligns (only implemented for polylines)
     TriangleMeshPolyLine *current_polyline =
       dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i]);
     TriangleMeshPolyLine *next_polyline =
       dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i+1]);

     if (current_polyline && next_polyline)
      {
       unsigned last_vertex = current_polyline->nvertex() - 1;
       next_polyline->vertex_coordinate(0) =
         current_polyline->vertex_coordinate(last_vertex);
      }
    }

  } // For n_boundaries - 1

 // Check wrap around
 // Auxiliary vertex for storing the vertex values of contiguous curves
 Vector<double> v1(2);

 // This is for getting the first coordinates of the first curve section
 Curve_section_pt[0]->initial_vertex_coordinate(v1);

 // Auxiliary vertex for storing the vertex values of contiguous curves
 Vector<double> v2(2);

 // This is for getting the last coordinates of the last curve section
 Curve_section_pt[n_boundaries-1]->final_vertex_coordinate(v2);

 double error = sqrt(pow(v2[0]-v1[0],2)+pow(v2[1]-v1[1],2));

 if (error>ToleranceForVertexMismatchInPolygons::Tolerable_error)
  {
   std::ostringstream error_stream;
   error_stream
   << "The start and end points of the first and last curve segment\n"
   << "boundary parts don't match when judged \nwith the tolerance of "
   << ToleranceForVertexMismatchInPolygons::Tolerable_error
   << " which is specified in the namespace \nvariable "
   << "ToleranceForVertexMismatchInPolygons::Tolerable_error.\n\n"
   << "Feel free to adjust this or to recompile the code without\n"
   << "paranoia if you think this is OK...\n"
   << std::endl;
   throw OomphLibError(
     error_stream.str(),
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
  }
 else
  {
   // Aligns (only implemented for polylines)
   TriangleMeshPolyLine *first_polyline =
     dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[0]);
   TriangleMeshPolyLine *last_polyline =
     dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[n_boundaries-1]);

   if (first_polyline && last_polyline)
    {
     unsigned last_vertex = last_polyline->nvertex() - 1;
     first_polyline->vertex_coordinate(0) =
       last_polyline->vertex_coordinate(last_vertex);
    }
  }

#endif

}



/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//=========================================================================
/// Constructor: Specify vector of pointers to TriangleMeshCurveSection
/// that define the boundary of the segments of the polygon.
/// Each TriangleMeshCurveSection has its own boundary ID and can contain
/// multiple (straight-line) segments. For consistency across the
/// various uses of this class, we insist that the closed boundary
/// is represented by at least two separate TriangleMeshCurveSection
/// whose joint vertices must be specified in both.
/// (This is to allow the setup of unique boundary coordinate(s)
/// around the polygon.) This may seem slightly annoying
/// in cases where a polygon really only represents a single
/// boundary, but...
/// Note: The specified vector of pointers must consist of only
/// TriangleMeshPolyLine elements. There is a checking on the PARANOID
/// mode for this constraint
//=========================================================================
TriangleMeshPolygon::TriangleMeshPolygon(
    const Vector<TriangleMeshCurveSection*>& boundary_polyline_pt,
    const Vector<double>& internal_point_pt) :
    TriangleMeshCurve(boundary_polyline_pt),
    TriangleMeshClosedCurve(boundary_polyline_pt, internal_point_pt),
    Enable_redistribution_of_segments_between_polylines(false),
    Can_update_configuration(false),
    Polygon_fixed(false)
{
 
#ifdef PARANOID
 
 // Get the number of polylines
 const unsigned n_bound = boundary_polyline_pt.size();
 
 // Check that all the constituent TriangleMeshCurveSection are
 // instance of TriangleMeshPolyLine
 for (unsigned p = 0; p < n_bound; p++)
  {
   TriangleMeshPolyLine *tmp_polyline_pt =
    dynamic_cast<TriangleMeshPolyLine*>(boundary_polyline_pt[p]);
   if (tmp_polyline_pt == 0)
    {
     std::ostringstream error_stream;
     error_stream
      << "The (" << p << ") TriangleMeshCurveSection is not a "
      << "TriangleMeshPolyLine.\nThe TriangleMeshPolygon object"
      << "is constituent of only TriangleMeshPolyLine objects.\n"
      << "Verify that all the constituent elements of the "
      << "TriangleMeshPolygon\nare instantiated as "
      << "TriangleMeshPolyLines." << std::endl;
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
 
 // Check that the polylines are contiguous
 bool contiguous=true;
 unsigned i_offensive=0;
 
 // Need at least two
 if (n_bound<2)
  {
   std::ostringstream error_stream;
   error_stream
    << "Sorry -- I'm afraid we insist that a closed curve is\n"
    << "specified by at least two separate TriangleMeshPolyLines.\n"
    << "You've only specified (" << n_bound << ")" << std::endl;
   throw OomphLibError(error_stream.str(),
                       "TriangleMeshPolygon::TriangleMeshPolygon()",
                       OOMPH_EXCEPTION_LOCATION);
  } // if (n_bound<2)
 
 // Does the last node of the polyline connect to the first one of the
 // next one (only up the last but one!)
 for(unsigned i=0;i<n_bound-1;i++)
  {
   // Get vector last vertex in current polyline
   unsigned last_vertex = (polyline_pt(i)->nvertex())-1;
   Vector<double> v1=polyline_pt(i)->
    vertex_coordinate(last_vertex);
   
   // Get vector to first vertex in next polyline
   Vector<double> v2=polyline_pt(i+1)->vertex_coordinate(0);

   // Work out error
   double error=sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2));

   // Is error accetable?
   if (error>ToleranceForVertexMismatchInPolygons::Tolerable_error)
    {
     contiguous=false;
     i_offensive=i;
     break;
    }
   // Align
   else
    {
     polyline_pt(i+1)->vertex_coordinate(0)=
      polyline_pt(i)->vertex_coordinate(last_vertex);
    }
  }
 
 // Does the last one connect to the first one?
 
 // Get vector last vertex last polyline
 unsigned last_vertex = (polyline_pt(n_bound-1)->nvertex())-1;
 Vector<double> v1=polyline_pt(n_bound-1)->
  vertex_coordinate(last_vertex);
 
 // Get vector first vertex first polyline
 Vector<double> v2=polyline_pt(0)->vertex_coordinate(0);
 double error=sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2));
 
 if (error>ToleranceForVertexMismatchInPolygons::Tolerable_error)
  {
   contiguous=false;
   i_offensive=n_bound-1;
  }
 else
  {
   polyline_pt(0)->vertex_coordinate(0)=
    polyline_pt(n_bound-1)->vertex_coordinate(last_vertex);
  }
 
 if (!contiguous)
  {
   std::ostringstream error_stream;
   error_stream
    << "The polylines specified \n"
    << "should define a closed geometry, i.e. the first/last vertex of\n"
    << "adjacent polylines should match.\n\n"
    << "Your polyline number "<< i_offensive
    <<" has no contiguous neighbour, when judged \nwith the tolerance of "
    << ToleranceForVertexMismatchInPolygons::Tolerable_error
    << " which is specified in the namespace \nvariable "
    << "ToleranceForVertexMismatchInPolygons::Tolerable_error.\n\n"
    << "Feel free to adjust this or to recompile the code without\n"
    << "paranoia if you think this is OK...\n"
    << std::endl;
   throw OomphLibError(error_stream.str(),
                       "TriangleMeshPolygon::TriangleMeshPolygon()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Check if internal point is actually located in bounding polygon
 // Reference: http://paulbourke.net/geometry/insidepoly/
 
 // Only checked if there is an internal hole
 if (!Internal_point_pt.empty())
  {
   
   // Vertex coordinates
   Vector<Vector<double> > polygon_vertex;
   
   // Total number of vertices
   unsigned nvertex=0;
   
   // Storage for first/last point on polyline for contiguousness check
   Vector<double> last_vertex(2);
   Vector<double> first_vertex(2);

   // Get vertices
   unsigned npolyline=boundary_polyline_pt.size();
   for (unsigned i=0;i<npolyline;i++)
    {
     // Number of vertices
     unsigned nvert=boundary_polyline_pt[i]->nvertex();
     for (unsigned j=0;j<nvert;j++)
      {
       // Check contiguousness
       if ((i>1)&&(j==0))
        {
         first_vertex=polyline_pt(i)->vertex_coordinate(j);
         double dist=sqrt(pow(first_vertex[0]-last_vertex[0],2)+
                          pow(first_vertex[1]-last_vertex[1],2));
         if (dist
             >ToleranceForVertexMismatchInPolygons::Tolerable_error)
          {
           std::ostringstream error_stream;
           error_stream
            << "The start and end points of polylines " << i
            << " and " <<i+1<< " don't match when judged\n"
            << "with the tolerance ("
            << ToleranceForVertexMismatchInPolygons::Tolerable_error
            << ") which is specified in the namespace \nvariable "
            << "ToleranceForVertexMismatchInPolygons::"
            << "Tolerable_error.\n\n"
            << "Feel free to adjust this or to recompile the "
            << "code without\n"
            << "paranoia if you think this is OK...\n"
            << std::endl;
           throw OomphLibError(
            error_stream.str(),
            "TriangleMeshPolygon::TriangleMeshPolygon()",
            OOMPH_EXCEPTION_LOCATION);
          }
        }
       // Get vertex (ignore end point)
       if (j<nvert-1)
        {
         polygon_vertex.push_back(
          polyline_pt(i)->vertex_coordinate(j));
        }
       // Prepare for check of contiguousness
       else
        {
         last_vertex=polyline_pt(i)->vertex_coordinate(j);
        }
      }
    }

   // Total number of vertices
   nvertex=polygon_vertex.size();

   // Counter for number of intersections
   unsigned intersect_counter=0;

   //Get first vertex
   Vector<double> p1=polygon_vertex[0];
   for (unsigned i=1;i<=nvertex;i++)
    {
     // Get second vertex by wrap-around
     Vector<double> p2 = polygon_vertex[i%nvertex];
     
     if (Internal_point_pt[1] > std::min(p1[1],p2[1]))
      {
       if (Internal_point_pt[1] <= std::max(p1[1],p2[1]))
        {
         if (Internal_point_pt[0] <= std::max(p1[0],p2[0]))
          {
           if (p1[1] != p2[1])
            {
             double xintersect =
              (Internal_point_pt[1]-p1[1])*(p2[0]-p1[0])/
              (p2[1]-p1[1])+p1[0];
             if ( (p1[0] == p2[0]) ||
                  (Internal_point_pt[0] <= xintersect) )
              {
               intersect_counter++;
              }
            }
          }
        }
      }
     p1 = p2;
    }
   
   // Even number of intersections: outside
   if (intersect_counter%2==0)
    {
     std::ostringstream error_stream;
     error_stream
      << "The internal point at "
      << Internal_point_pt[0]<< " "
      << Internal_point_pt[1]
      << " isn't in the polygon that describes the internal closed "
      << "curve!\nPolygon vertices are at: \n";
     for (unsigned i=0;i<nvertex;i++)
      {
       error_stream << polygon_vertex[i][0] << " "
                    << polygon_vertex[i][1] << "\n";
      }
     error_stream
      << "This may be because the internal point is defined by a\n"
      << "GeomObject that has deformed so much that it's \n"
      << "swept over the (initial) internal point.\n"
      << "If so, you should update the position of the internal point. \n"
      << "This could be done automatically by generating \n"
      << "an internal mesh inside the polygon and using one\n"
      << "of its internal nodes as the internal point. Actually not \n"
      << "why triangle doesn't do that automatically....\n";
     throw OomphLibError(
      error_stream.str(),
      "TriangleMeshPolygon::TriangleMeshPolygon()",
      OOMPH_EXCEPTION_LOCATION);

    }

  }
 
#endif
 
}





/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////





//=====================================================================
/// Class defining an open curve for the Triangle mesh generation
//=====================================================================
TriangleMeshOpenCurve::TriangleMeshOpenCurve(
   const Vector<TriangleMeshCurveSection*> &curve_section_pt)
 : TriangleMeshCurve(curve_section_pt)
 {

#ifdef PARANOID

 // Matching of curve sections i.e. the last vertex of
 // the i curve section should match with the first
 // vertex of the i+1 curve section

 // Total number of boundaries
 unsigned n_boundaries = Curve_section_pt.size();

 // Check last point of each boundary bit coincides with first point
 // on next one
 for (unsigned i=0;i<n_boundaries-1;i++)
  {

   // Auxiliary vertex for storing the vertex values of contiguous curves
   Vector<double> v1(2);
   Vector<double> v2(2);

   // This is for getting the final coordinates of the i curve section
   Curve_section_pt[i]->final_vertex_coordinate(v1);

   // This is for the start coordinates of the i+1 curve section
   Curve_section_pt[i+1]->initial_vertex_coordinate(v2);

   // Work out error
   double error = sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2));
   if (error>ToleranceForVertexMismatchInPolygons::Tolerable_error)
    {
     std::ostringstream error_stream;
     error_stream
     << "The start and end points of curve section boundary parts " << i
     << " and " <<i+1<< " don't match when judged \nwith the tolerance of "
     << ToleranceForVertexMismatchInPolygons::Tolerable_error
     << " which is specified in the namespace \nvariable "
     << "ToleranceForVertexMismatchInPolygons::Tolerable_error.\n\n"
     << "Feel free to adjust this or to recompile the code without\n"
     << "paranoia if you think this is OK...\n"
     << std::endl;
     throw OomphLibError(
       error_stream.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
    }
   else
    {
     // Aligns (only implemented for polylines)
     TriangleMeshPolyLine *current_polyline =
       dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i]);
     TriangleMeshPolyLine *next_polyline =
       dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i+1]);

     if (current_polyline && next_polyline)
      {
       unsigned last_vertex = current_polyline->nvertex() - 1;
       next_polyline->vertex_coordinate(0) =
         current_polyline->vertex_coordinate(last_vertex);
      }
    }

  } // For n_boundaries - 1

#endif

 }






//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////






 //=================================================
 ///  Helper namespace for BCInfo object used
 /// in the identification of boundary elements.
 //=================================================
namespace TriangleBoundaryHelper
{

 /// Structure for Boundary Informations
 struct BCInfo
 {

  /// Face ID
  unsigned Face_id;
 
  /// Boundary ID
  unsigned Boundary;
 
  /// Pointer to bulk finite element
  FiniteElement* FE_pt;
  
 };
 
}

#ifdef OOMPH_HAS_TRIANGLE_LIB



//==============================================================
/// Dump the triangulateio structure to a dump file and
/// record boundary coordinates of boundary nodes
//==============================================================
 void TriangleMeshBase::dump_triangulateio(std::ostream &dump_file)
 {
  TriangleHelper::dump_triangulateio(Triangulateio,dump_file);

#ifdef OOMPH_HAS_MPI
  // If the mesh is not distributed then process what follows
  if (!this->is_mesh_distributed())
   {
#endif // #ifdef OOMPH_HAS_MPI
    
    // Loop over all boundary nodes and dump out boundary coordinates
    // if they exist
    Vector<double> zeta(1);
    unsigned nb=nboundary();
    for (unsigned b=0;b<nb;b++)
     {
      if (Boundary_coordinate_exists[b])
       {
        dump_file << "1 # Boundary coordinate for boundary " << b 
                  << " does exist\n";
        unsigned nnod=nboundary_node(b);
        dump_file << nnod << " # Number of dumped boundary nodes\n";
        for (unsigned j=0;j<nnod;j++)
         {
          Node* nod_pt=boundary_node_pt(b,j);
          nod_pt->get_coordinates_on_boundary(b,zeta);
          dump_file << zeta[0] << std::endl;
         }
        dump_file << "-999 # Done boundary coords for boundary " << b << "\n";
       }
      else
       {
        dump_file << "0 # Boundary coordinate for boundary " << b 
                  << " does not exist\n";
      
       }
     }
  
#ifdef OOMPH_HAS_MPI
   }
#endif // #ifdef OOMPH_HAS_MPI
  
 }


//==============================================================
/// Regenerate the mesh from a dumped triangulateio file
/// and dumped boundary coordinates of boundary nodes
//==============================================================
 void TriangleMeshBase::remesh_from_triangulateio(std::istream &restart_file)
 {


#ifdef PARANOID
  // Record number of boundaries
  unsigned nbound_old=nboundary();
#endif

  //Clear the existing triangulate io
  TriangleHelper::clear_triangulateio(Triangulateio);

  //Read the data into the file
  TriangleHelper::read_triangulateio(restart_file,Triangulateio);

  //Now remesh from the new data structure
  this->remesh_from_internal_triangulateio();
   
#ifdef OOMPH_HAS_MPI
  // If the mesh is not distributed then process what follows
  if (!this->is_mesh_distributed())
   {
#endif // #ifdef OOMPH_HAS_MPI
    
#ifdef PARANOID
    // Record number of boundary nodes after remesh
    unsigned nbound_new=nboundary();
    if (nbound_new!=nbound_old)
     {
      std::ostringstream error_stream;
      error_stream << "Number of boundaries before remesh from triangulateio, " 
                   << nbound_new 
                   << ",\ndoesn't match number boundaries afterwards, " 
                   << nbound_old 
                   << ". Have you messed \naround with boundary nodes in the "
                   << "derived mesh constructor (or after calling \nit)? If so,"
                   << " the dump/restart won't work as written at the moment.";
      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    
    // Loop over all boundary nodes and read boundary coordinates 
    // if they exist
    Vector<double> zeta(1);
    std::string input_string;
    unsigned nb=nboundary();
    for (unsigned b=0;b<nb;b++)
     {
      // Read line up to termination sign
      getline(restart_file,input_string,'#');

      // Ignore rest of line
      restart_file.ignore(80,'\n');

      // Did boundary coordinate exist?
      const unsigned bound_coord_exists=atoi(input_string.c_str());    
      if (bound_coord_exists==1)
       {
        // Remember it!
        Boundary_coordinate_exists[b]=true;

        // Read line up to termination sign
        getline(restart_file,input_string,'#');
      
        // Ignore rest of line
        restart_file.ignore(80,'\n');
      
        // How many nodes did we dump?
        const unsigned nnod_dumped=atoi(input_string.c_str());

        // Does it match?
        unsigned nnod=nboundary_node(b);
        if (nnod!=nnod_dumped)
         {
          std::ostringstream error_stream;
          error_stream << "Number of dumped boundary nodes " 
                       << nnod_dumped 
                       << " doesn't match number of nodes on boundary " 
                       << b << ": " << nnod << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
         }

        // Loop over all nodes
        for (unsigned j=0;j<nnod;j++)
         {
          // Read line up to termination sign
          getline(restart_file,input_string);
        
          // Boundary coordinate
          zeta[0]=atof(input_string.c_str());

          // Set it
          Node* nod_pt=boundary_node_pt(b,j);
          nod_pt->set_coordinates_on_boundary(b,zeta);
         }

        // Read line up to termination sign
        getline(restart_file,input_string,'#');
      
        // Ignore rest of line
        restart_file.ignore(80,'\n');
      
        // Have we reached the end?
        const int check=atoi(input_string.c_str());
        if (check!=-999)
         {
          std::ostringstream error_stream;
          error_stream << "Haven't read all nodes on boundary "<< b 
                       << std::endl; 
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
         }
       }
      else
       {
        oomph_info << "Restart: Boundary coordinate for boundary " << b 
                   << " does not exist.\n";
        
       }
     }
#ifdef OOMPH_HAS_MPI  
   } // if (!this->is_mesh_distributed())
#endif // #ifdef OOMPH_HAS_MPI
  
 }





//==============================================================
/// Write a Triangulateio_object file of the TriangulateIO object
/// String s is add to assign a different value for
/// input and/or output structure.
/// The function give the same result of the "report" function
/// included in the tricall.c, esternal_src.
//==============================================================
 void TriangleMeshBase::write_triangulateio(TriangulateIO& triangle, 
                                            std::string &s)
 {

  std::ofstream outfile;
  char filename[100];

  sprintf(filename,"Triangulateio_object_%s.dat",s.c_str());
  outfile.open(filename);
  outfile <<"# Triangulateio object values:\n\n"<<std::endl;

  // Write points coordinates
  if(triangle.numberofpoints!=0)
   {
    outfile <<"# Triangulateio number of points is:"
            <<triangle.numberofpoints <<std::endl;
   }
  if(triangle.pointlist != NULL)
   {
    outfile <<"# Vertex coordinates are:"<<std::endl;
    for(int k=0;k<triangle.numberofpoints*2;k+=2)
     {
      outfile << (k*0.5)+1 << " "
              << triangle.pointlist[k] << " "
              << triangle.pointlist[k+1] << std::endl;
     }
   }

  // Write points attribute list
  if(triangle.numberofpointattributes!=0)
   {
    outfile <<"# Triangulateio number of points attributelist is:"
            <<triangle.numberofpointattributes <<std::endl;
   }
  if(triangle.pointattributelist != NULL)
   {
       
    outfile <<"# Vertex attribute are:"<<std::endl;
    for(int k=0;k<triangle.numberofpointattributes;k++)
     {
      outfile << triangle.pointattributelist[k] << std::endl;
     }
   }

  // Write point markers list
  if(triangle.pointmarkerlist != NULL)
   {
    outfile <<"# Vertex Markers are:"<<std::endl;
    for(int k=0;k<triangle.numberofpoints;k++)
     {
      outfile << triangle.pointmarkerlist[k] << std::endl;
     }
   }
  
  // Write the 1.node file used by the showme function  
  std::ofstream nodefile;
  char nodename[100];

  sprintf(nodename,"file_%s.1.node",s.c_str());
  nodefile.open(nodename);
  nodefile <<triangle.numberofpoints<<" 2 "
           <<triangle.numberofpointattributes
           <<" 0"<<std::endl;
  for(int j=0;j<triangle.numberofpoints*2;j+=2)
   {
    nodefile <<(j/2)+1<<" " << triangle.pointlist[j] << " "
            << triangle.pointlist[j+1] << std::endl;
   }
  nodefile.close();
  

  // Write segments edge elements
  if(triangle.numberofsegments!=0)
   {
    outfile <<"# The number of segments is:"
            <<triangle.numberofsegments<<std::endl;
   }
  if(triangle.segmentlist != NULL)
   {
    outfile <<"# Segments are:"<<std::endl;
    for(int k=0;k<triangle.numberofsegments*2;k+=2)
     {
      outfile << triangle.segmentlist[k] << "  "
              << triangle.segmentlist[k+1] <<std::endl;
     }
   }

  // Write segments markers list
  if(triangle.segmentmarkerlist != NULL)
   {
    outfile <<"# Segments Markers are:"<<std::endl;
    for(int k=0;k<triangle.numberofsegments;k++)
     {
      outfile << triangle.segmentmarkerlist[k] << std::endl;
     }
   }

  // Write regions
  if(triangle.numberofregions!=0)
   {
    outfile <<"# The number of region is:"
            <<triangle.numberofregions<<std::endl;
   }

  // Write holes
  if(triangle.numberofholes!=0)
   {
    outfile <<"# The number of holes is:"
            <<triangle.numberofholes<<std::endl;
   }
  if(triangle.holelist != NULL)
   {
    outfile <<"#  Holes are:"<<std::endl;
    for(int k=0;k<triangle.numberofholes*2;k+=2)
     {
      outfile << triangle.holelist[k] << "  "
              << triangle.holelist[k+1] <<std::endl;
     }
   }

  // Write triangles
  if(triangle.numberoftriangles!=0)
   {
    outfile <<"# Triangulateio number of triangles:"
            <<triangle.numberoftriangles <<std::endl;
   }
  if(triangle.numberofcorners!=0)
   {
    outfile <<"# Triangulateio number of corners:"
            <<triangle.numberofcorners<<std::endl;
   }
  if(triangle.numberoftriangleattributes!=0)
   {
    outfile <<"# Triangulateio number of triangles attributes:"
            <<triangle.numberoftriangleattributes <<std::endl;
   }
  if(triangle.trianglelist != NULL)
   {
       
    outfile <<"# Traingles are:"<<std::endl;
    for(int k=0;k<triangle.numberoftriangles*3;k+=3)
     {
      outfile << triangle.trianglelist[k] << " "
              << triangle.trianglelist[k+1] << " "
              << triangle.trianglelist[k+2] << std::endl;
     }
   }

  if(triangle.trianglearealist != NULL)
   {
    outfile <<"# Triangle's areas are:"<<std::endl;
     for(int k=0;k<triangle.numberoftriangles;k++)
     {
      outfile << triangle.trianglearealist[k] << std::endl;
     }
   }

  if(triangle.trianglelist != NULL)
   {
    
    // Write the 1.ele file used by the showme function  
    std::ofstream elefile;
    char elename[100];
  
    sprintf(elename,"file_%s.1.ele",s.c_str());
    elefile.open(elename);
    elefile <<triangle.numberoftriangles<<" 3 0"<<std::endl;
    for(int j=0;j<triangle.numberoftriangles*3;j+=3)
     {
      elefile <<(j/3)+1<<" "<< triangle.trianglelist[j] << " "
              << triangle.trianglelist[j+1] << " "
              << triangle.trianglelist[j+2] << std::endl;
     }
    elefile.close();
   }
  
  outfile.close();
 }

#endif

//================================================================
/// Setup lookup schemes which establish which elements are located
/// next to which boundaries (Doc to outfile if it's open).
//================================================================
void TriangleMeshBase::setup_boundary_element_info(std::ostream &outfile)
{

 //Should we document the output here
 bool doc = false;

 if(outfile) doc = true;

 // Number of boundaries
 unsigned nbound=nboundary();
 
 // Wipe/allocate storage for arrays
 Boundary_element_pt.clear();
 Face_index_at_boundary.clear();
 Boundary_element_pt.resize(nbound);
 Face_index_at_boundary.resize(nbound);
 
 // Temporary vector of vectors of pointers to elements on the boundaries: 
 // This is a vector to ensure that order is strictly preserved
 Vector<Vector<FiniteElement*> > vector_of_boundary_element_pt;
 vector_of_boundary_element_pt.resize(nbound);
 
 // Matrix map for working out the fixed face for elements on boundary
 MapMatrixMixed<unsigned,FiniteElement*, int > 
  face_identifier;
 
 // Loop over elements
 //-------------------
 unsigned nel=nelement();
 
 // Get pointer to vector of boundaries that the
 // node lives on
 Vector<std::set<unsigned>*> boundaries_pt(3,0);
 
 // Data needed to deal with edges through the
 // interior of the domain
 std::map<Edge,unsigned > edge_count;
 std::map<Edge,TriangleBoundaryHelper::BCInfo> 
  edge_bcinfo;
 std::map<Edge,TriangleBoundaryHelper::BCInfo>
  face_info;
 MapMatrixMixed<unsigned,FiniteElement*, int > face_count;
 Vector<unsigned> bonus(nbound);

 // When using internal boundaries, an edge can be related to more than
 // one element (because of both sides of the internal boundaries)
 std::map<Edge,Vector<TriangleBoundaryHelper::BCInfo> > edge_internal_bnd;
 
 for (unsigned e=0;e<nel;e++)
  {
   // Get pointer to element
   FiniteElement* fe_pt=finite_element_pt(e);
   
   if (doc) {outfile << "Element: " << e << " " << fe_pt << std::endl;}
   
   // Only include 2D elements! Some meshes contain interface elements too.
   if (fe_pt->dim()==2)
    {
     // Loop over the element's nodes and find out which boundaries they're on
     // ----------------------------------------------------------------------

     //We need only loop over the corner nodes
     for(unsigned i=0;i<3;i++)
      {
       fe_pt->node_pt(i)->get_boundaries_pt(boundaries_pt[i]);
      }

     //Find the common boundaries of each edge
     Vector<std::set<unsigned> > edge_boundary(3);
    
     //Edge 0 connects points 1 and 2
     //-----------------------------

     if(boundaries_pt[1] && boundaries_pt[2])
      {
       // Create the corresponding edge
       Edge edge0(fe_pt->node_pt(1),fe_pt->node_pt(2));

       // Update infos about this edge
       TriangleBoundaryHelper::BCInfo info;
       info.Face_id=0;
       info.FE_pt = fe_pt;
	
       std::set_intersection(boundaries_pt[1]->begin(),boundaries_pt[1]->end(),
                             boundaries_pt[2]->begin(),boundaries_pt[2]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge_boundary[0],edge_boundary[0].begin()));
       std::set<unsigned>::iterator it0=edge_boundary[0].begin();

       // Edge does exist:
       if ( edge_boundary[0].size() > 0 )
        {
         info.Boundary=*it0;
         
         // How many times this edge has been visited
         edge_count[edge0]++;
         
         // Update edge_bcinfo
         edge_bcinfo.insert(std::make_pair(edge0,info));
         
         // ... and also update the info associated with internal bnd
         edge_internal_bnd[edge0].push_back(info);
        }
      }
     
     //Edge 1 connects points 0 and 2
     //-----------------------------

     if(boundaries_pt[0] && boundaries_pt[2])
      {
       std::set_intersection(boundaries_pt[0]->begin(),boundaries_pt[0]->end(),
                             boundaries_pt[2]->begin(),boundaries_pt[2]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge_boundary[1],edge_boundary[1].begin()));

       // Create the corresponding edge
       Edge edge1(fe_pt->node_pt(0),fe_pt->node_pt(2));

       // Update infos about this edge
       TriangleBoundaryHelper::BCInfo info;
       info.Face_id=1;
       info.FE_pt = fe_pt;
       std::set<unsigned>::iterator it1=edge_boundary[1].begin();

       // Edge does exist:
       if ( edge_boundary[1].size() > 0)
        {
         info.Boundary=*it1;
         
         // How many times this edge has been visited
         edge_count[edge1]++;  
         
         // Update edge_bcinfo              
         edge_bcinfo.insert(std::make_pair(edge1,info));
         
         // ... and also update the info associated with internal bnd
         edge_internal_bnd[edge1].push_back(info);
        }
      }
     
     //Edge 2 connects points 0 and 1
     //-----------------------------

     if(boundaries_pt[0] && boundaries_pt[1])
      {
       std::set_intersection(boundaries_pt[0]->begin(),boundaries_pt[0]->end(),
                             boundaries_pt[1]->begin(),boundaries_pt[1]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge_boundary[2],edge_boundary[2].begin()));

       // Create the corresponding edge
       Edge edge2(fe_pt->node_pt(0),fe_pt->node_pt(1));

       // Update infos about this edge
       TriangleBoundaryHelper::BCInfo info;
       info.Face_id=2;
       info.FE_pt = fe_pt;
       std::set<unsigned>::iterator it2=edge_boundary[2].begin();

       // Edge does exist:
       if ( edge_boundary[2].size() > 0)
        {
         info.Boundary=*it2;
         
         // How many times this edge has been visited
         edge_count[edge2]++;
         
         // Update edge_bcinfo
         edge_bcinfo.insert(std::make_pair(edge2,info));
         
         // ... and also update the info associated with internal bnd
         edge_internal_bnd[edge2].push_back(info);
        }
      }
     
     
#ifdef PARANOID
     
     // Check if edge is associated with multiple boundaries
     
     //We now know whether any edges lay on the boundaries
     for(unsigned i=0;i<3;i++)
      {
       //How many boundaries are there
       unsigned count = 0;

       //Loop over all the members of the set and add to the count
       //and set the boundary
       for(std::set<unsigned>::iterator it=edge_boundary[i].begin();
           it!=edge_boundary[i].end();++it)
        {
         ++count;
        }

       //If we're on more than one boundary, this is weird, so die
       if(count > 1)
        {
         std::ostringstream error_stream;
         error_stream << "Edge " << i << " is located on " << 
          count << " boundaries.\n";
         error_stream << "This is rather strange, so I'm going to die\n";
         throw OomphLibError(
          error_stream.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
        }
      }
     
#endif
     
     // Now we set the pointers to the boundary sets to zero
     for(unsigned i=0;i<3;i++) {boundaries_pt[i] = 0;}
     
    }
  } //end of loop over all elements
 
 // Loop over all edges that are located on a boundary
 typedef std::map<Edge,TriangleBoundaryHelper::BCInfo>::iterator ITE;
 for (ITE it=edge_bcinfo.begin();
      it!=edge_bcinfo.end();
      it++)
  {
   Edge current_edge = it->first;
   unsigned  bound=it->second.Boundary;
   
   // If the edge has been visited only once
   if(edge_count[current_edge]==1)
    {
     // Count the edges that are on the same element and on the same boundary
     face_count(static_cast<unsigned>(bound),it->second.FE_pt)=  
      face_count(static_cast<unsigned>(bound),it->second.FE_pt) + 1;
     
     //If such edges exist, let store the corresponding element
     if( face_count(bound,it->second.FE_pt) > 1)
      {
       // Update edge's infos
       TriangleBoundaryHelper::BCInfo info;
       info.Face_id=it->second.Face_id;
       info.FE_pt = it->second.FE_pt;
       info.Boundary=it->second.Boundary;
       
       // Add it to FIinfo, that stores infos of problematic elements
       face_info.insert(std::make_pair(current_edge,info));
       
       //How many edges on which boundary have to be added
       bonus[bound]++;
      }
     else
      {
       //Add element and face to the appropriate vectors
       // Does the pointer already exits in the vector
       Vector<FiniteElement*>::iterator b_el_it =
        std::find(vector_of_boundary_element_pt[
                   static_cast<unsigned>(bound)].begin(),
                  vector_of_boundary_element_pt[
                   static_cast<unsigned>(bound)].end(),
                  it->second.FE_pt);
       
       //Only insert if we have not found it (i.e. got to the end)
       if(b_el_it == vector_of_boundary_element_pt[
           static_cast<unsigned>(bound)].end())
        {
         vector_of_boundary_element_pt[static_cast<unsigned>(bound)].
          push_back(it->second.FE_pt);
        }
       
       //set_of_boundary_element_pt[static_cast<unsigned>(bound)].insert(
       // it->second.FE_pt);
       face_identifier(static_cast<unsigned>(bound),it->second.FE_pt) = 
        it->second.Face_id;
      }
     
    }
   
  }  //End of "adding-boundaries"-loop
 
  
 // Now copy everything across into permanent arrays
 //-------------------------------------------------

 // Loop over boundaries
 for (unsigned i=0;i<nbound;i++)
  {
   // Number of elements on this boundary that have to be added 
   // in addition to other elements
   unsigned bonus1=bonus[i];
   
   // Number of elements on this boundary
   unsigned nel=vector_of_boundary_element_pt[i].size() + bonus1;

   // Allocate storage for the coordinate identifiers
   Face_index_at_boundary[i].resize(nel);

   unsigned e_count=0;
   typedef Vector<FiniteElement*>::iterator IT;
   for (IT it=vector_of_boundary_element_pt[i].begin();
        it!=vector_of_boundary_element_pt[i].end();
        it++)
    {    
     // Recover pointer to element
     FiniteElement* fe_pt=*it;

     // Add to permanent storage
     Boundary_element_pt[i].push_back(fe_pt);

     Face_index_at_boundary[i][e_count] = face_identifier(i,fe_pt);

     // Increment counter
     e_count++;

    }
   // We add the elements that have two or more edges on this boundary
   for (ITE itt= face_info.begin(); itt!= face_info.end(); itt++)
    {
     if (itt->second.Boundary==i)
      {
       // Add to permanent storage
       Boundary_element_pt[i].push_back(itt->second.FE_pt);

       Face_index_at_boundary[i][e_count] = itt->second.Face_id;

       e_count++;
      }

    }

  } //End of loop over boundaries
 


 // Doc?
 //-----
 if (doc)
  {
   // Loop over boundaries
   for (unsigned i=0;i<nbound;i++)
    {
     unsigned nel=Boundary_element_pt[i].size();
     outfile << "Boundary: " << i
             << " is adjacent to " << nel
             << " elements" << std::endl;
     
     // Loop over elements on given boundary
     for (unsigned e=0;e<nel;e++)
      {
       FiniteElement* fe_pt=Boundary_element_pt[i][e];
       outfile << "Boundary element:" <<  fe_pt
               << " Face index of boundary is " 
               <<  Face_index_at_boundary[i][e] << std::endl;
      }
    }
  }
 

 // Lookup scheme has now been setup yet
 Lookup_for_elements_next_boundary_is_setup=true;

}


}
