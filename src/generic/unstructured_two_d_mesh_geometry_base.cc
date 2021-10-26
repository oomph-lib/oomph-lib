// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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

#include "unstructured_two_d_mesh_geometry_base.h"


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

namespace oomph
{
#ifdef OOMPH_HAS_TRIANGLE_LIB

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
      if (clear_hole_data)
      {
        free(triangulate_io.regionlist);
      }
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
      triangle_io.pointlist = (double*)NULL;
      triangle_io.pointattributelist = (double*)NULL;
      triangle_io.pointmarkerlist = (int*)NULL;
      triangle_io.numberofpoints = 0;
      triangle_io.numberofpointattributes = 0;

      // Initialize the triangle list
      triangle_io.trianglelist = (int*)NULL;
      triangle_io.triangleattributelist = (double*)NULL;
      triangle_io.trianglearealist = (double*)NULL;
      triangle_io.neighborlist = (int*)NULL;
      triangle_io.numberoftriangles = 0;
      triangle_io.numberofcorners = 0;
      triangle_io.numberoftriangleattributes = 0;

      // Initialize the segment list
      triangle_io.segmentlist = (int*)NULL;
      triangle_io.segmentmarkerlist = (int*)NULL;
      triangle_io.numberofsegments = 0;

      // Initialise hole list
      triangle_io.holelist = (double*)NULL;
      triangle_io.numberofholes = 0;

      // Initialize region list
      triangle_io.regionlist = (double*)NULL;
      triangle_io.numberofregions = 0;

      // Initialize edge list
      triangle_io.edgelist = (int*)NULL;
      triangle_io.edgemarkerlist = (int*)NULL;
      triangle_io.normlist = (double*)NULL;
      triangle_io.numberofedges = 0;
    }


    /// Make (partial) deep copy of TriangulateIO object. We only copy
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
      triangle_out.pointlist =
        (double*)malloc(triangle_out.numberofpoints * 2 * sizeof(double));
      for (int j = 0; j < triangle_out.numberofpoints * 2; j++)
      {
        triangle_out.pointlist[j] = triangle_io.pointlist[j];
      }

      triangle_out.pointmarkerlist =
        (int*)malloc(triangle_out.numberofpoints * sizeof(int));
      for (int j = 0; j < triangle_out.numberofpoints; j++)
      {
        triangle_out.pointmarkerlist[j] = triangle_io.pointmarkerlist[j];
      }

      // Warn about laziness...
      if (!quiet)
      {
        if ((triangle_io.pointattributelist != 0) ||
            (triangle_io.numberofpointattributes != 0))
        {
          OomphLibWarning(
            "Point attributes are not currently copied across",
            "TriangleHelper::deep_copy_of_triangulateio_representation",
            OOMPH_EXCEPTION_LOCATION);
        }
      }


      // Triangle data
      triangle_out.numberoftriangles = triangle_io.numberoftriangles;
      triangle_out.trianglelist =
        (int*)malloc(triangle_out.numberoftriangles * 3 * sizeof(int));
      for (int j = 0; j < triangle_out.numberoftriangles * 3; j++)
      {
        triangle_out.trianglelist[j] = triangle_io.trianglelist[j];
      }


      // Copy over the triangle attribute data
      triangle_out.numberoftriangleattributes =
        triangle_io.numberoftriangleattributes;
      triangle_out.triangleattributelist = (double*)malloc(
        triangle_out.numberoftriangles *
        triangle_out.numberoftriangleattributes * sizeof(double));
      for (int j = 0; j < (triangle_out.numberoftriangles *
                           triangle_out.numberoftriangleattributes);
           ++j)
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

        if ((triangle_io.trianglearealist != 0))
        {
          OomphLibWarning(
            "Triangle areas are not currently copied across",
            "TriangleHelper::deep_copy_of_triangulateio_representation",
            OOMPH_EXCEPTION_LOCATION);
        }

        if ((triangle_io.neighborlist != 0))
        {
          OomphLibWarning(
            "Triangle neighbours are not currently copied across",
            "TriangleHelper::deep_copy_of_triangulateio_representation",
            OOMPH_EXCEPTION_LOCATION);
        }
      }


      triangle_out.numberofcorners = triangle_io.numberofcorners;

      // Segment data
      triangle_out.numberofsegments = triangle_io.numberofsegments;
      triangle_out.segmentlist =
        (int*)malloc(triangle_out.numberofsegments * 2 * sizeof(int));
      for (int j = 0; j < triangle_out.numberofsegments * 2; j++)
      {
        triangle_out.segmentlist[j] = triangle_io.segmentlist[j];
      }
      triangle_out.segmentmarkerlist =
        (int*)malloc(triangle_out.numberofsegments * sizeof(int));
      for (int j = 0; j < triangle_out.numberofsegments; j++)
      {
        triangle_out.segmentmarkerlist[j] = triangle_io.segmentmarkerlist[j];
      }


      // Region data
      triangle_out.numberofregions = triangle_io.numberofregions;
      triangle_out.regionlist =
        (double*)malloc(triangle_out.numberofregions * 4 * sizeof(double));
      for (int j = 0; j < triangle_out.numberofregions * 4; ++j)
      {
        triangle_out.regionlist[j] = triangle_io.regionlist[j];
      }

      // Hole data
      triangle_out.numberofholes = triangle_io.numberofholes;
      triangle_out.holelist =
        (double*)malloc(triangle_out.numberofholes * 2 * sizeof(double));
      for (int j = 0; j < triangle_out.numberofholes * 2; j++)
      {
        triangle_out.holelist[j] = triangle_io.holelist[j];
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

        if ((triangle_io.edgelist != 0) || (triangle_io.numberofedges != 0))
        {
          OomphLibWarning(
            "Edges are not currently copied across",
            "TriangleHelper::deep_copy_of_triangulateio_representation",
            OOMPH_EXCEPTION_LOCATION);
        }

        if ((triangle_io.edgemarkerlist != 0))
        {
          OomphLibWarning(
            "Edge markers are not currently copied across",
            "TriangleHelper::deep_copy_of_triangulateio_representation",
            OOMPH_EXCEPTION_LOCATION);
        }

        if ((triangle_io.normlist != 0))
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

    /// Write the triangulateio data to disk as a poly file,
    /// mainly used for debugging
    void write_triangulateio_to_polyfile(TriangulateIO& triangle_io,
                                         std::ostream& poly_file)
    {
      // Up the precision dramatiacally
      poly_file.precision(20);

      // Output the number of points and their attributes
      // Store the number of attributes
      const int n_attr = triangle_io.numberofpointattributes;
      poly_file << triangle_io.numberofpoints << "  " << 2 << " " << n_attr
                << " ";
      // Determine whether there are point markers
      bool point_markers = true;
      if (triangle_io.pointmarkerlist == NULL)
      {
        point_markers = false;
      }
      // Indicate this in the file
      poly_file << point_markers << "\n";

      // Now output the point data
      poly_file << "#Points\n";
      for (int n = 0; n < triangle_io.numberofpoints; ++n)
      {
        // Output the point number and x and y coordinates
        poly_file << n + 1 << " " << triangle_io.pointlist[2 * n] << " "
                  << triangle_io.pointlist[2 * n + 1] << " ";
        // Output any attributes
        for (int i = 0; i < n_attr; ++i)
        {
          poly_file << triangle_io.pointattributelist[n_attr * n + i] << " ";
        }
        // Output the boundary marker
        if (point_markers)
        {
          poly_file << triangle_io.pointmarkerlist[n] << " ";
        }
        poly_file << "\n";
      }

      // Now move onto the segments
      poly_file << "#Lines\n";
      poly_file << triangle_io.numberofsegments << " ";
      // Determine whether there are segment markers
      bool seg_markers = true;
      if (triangle_io.segmentmarkerlist == NULL)
      {
        seg_markers = false;
      }
      // Output this info in the file
      poly_file << seg_markers << "\n";

      // Now output the segment data
      for (int n = 0; n < triangle_io.numberofsegments; ++n)
      {
        poly_file << n + 1 << " " << triangle_io.segmentlist[2 * n] << " "
                  << triangle_io.segmentlist[2 * n + 1] << " ";
        // If there is a boundary marker output
        if (seg_markers)
        {
          poly_file << triangle_io.segmentmarkerlist[n] << " ";
        }
        poly_file << "\n";
      }

      // Now output the number of holes
      poly_file << "#No holes\n";
      poly_file << triangle_io.numberofholes << "\n";
      // Output the hole data
      for (int h = 0; h < triangle_io.numberofholes; ++h)
      {
        poly_file << h + 1 << " " << triangle_io.holelist[2 * h] << " "
                  << triangle_io.holelist[2 * h + 1] << "\n";
      }

      // Now output the number of regions
      poly_file << "#Assignment of attributes to regions\n";
      poly_file << triangle_io.numberofregions << "\n";

      // Loop over the regions
      for (int r = 0; r < triangle_io.numberofregions; ++r)
      {
        poly_file << r + 1 << " ";
        for (unsigned i = 0; i < 4; i++)
        {
          poly_file << triangle_io.regionlist[4 * r + i] << " ";
        }
        poly_file << "\n";
      }
    }


    /// Create a triangulateio data file from ele node and poly
    /// files. This is used if the mesh is generated by using Triangle
    /// externally. The triangulateio structure is required to dump the mesh
    /// topology for restarts.
    void create_triangulateio_from_polyfiles(
      const std::string& node_file_name,
      const std::string& element_file_name,
      const std::string& poly_file_name,
      TriangulateIO& triangle_io,
      bool& use_attributes)
    {
      // Initialise the TriangulateIO data structure
      initialise_triangulateio(triangle_io);

      // Process element file
      std::ifstream element_file(element_file_name.c_str(), std::ios_base::in);

      // Check that the file actually opened correctly
      if (!element_file.is_open())
      {
        std::string error_msg("Failed to open element file: ");
        error_msg += "\"" + element_file_name + "\".";
        throw OomphLibError(
          error_msg, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Read in the number of elements
      element_file >> triangle_io.numberoftriangles;
      const unsigned n_element =
        static_cast<unsigned>(triangle_io.numberoftriangles);

      // Read in the number of nodes per element
      element_file >> triangle_io.numberofcorners;
      const unsigned n_local_node =
        static_cast<unsigned>(triangle_io.numberofcorners);

      // Read in the element attributes
      element_file >> triangle_io.numberoftriangleattributes;
      const unsigned n_attributes =
        static_cast<unsigned>(triangle_io.numberoftriangleattributes);

      // Allocate storage in the data structure
      triangle_io.trianglelist =
        (int*)malloc(triangle_io.numberoftriangles *
                     triangle_io.numberofcorners * sizeof(int));

      if (n_attributes > 0)
      {
        triangle_io.triangleattributelist = (double*)malloc(
          triangle_io.numberoftriangles *
          triangle_io.numberoftriangleattributes * sizeof(double));
      }

      // Dummy storage
      int dummy_element_number;

      // Initialise counter
      unsigned counter = 0;
      unsigned counter2 = 0;

      // Read global node numbers for all elements
      for (unsigned e = 0; e < n_element; e++)
      {
        element_file >> dummy_element_number;
        for (unsigned j = 0; j < n_local_node; j++)
        {
          element_file >> triangle_io.trianglelist[counter];
          ++counter;
        }
        for (unsigned j = 0; j < n_attributes; j++)
        {
          element_file >> triangle_io.triangleattributelist[counter2];
          ++counter2;
        }
      }
      // Close the element file
      element_file.close();

      // Process node file
      // -----------------
      std::ifstream node_file(node_file_name.c_str(), std::ios_base::in);

      // Check that the file actually opened correctly
      if (!node_file.is_open())
      {
        std::string error_msg("Failed to open node file: ");
        error_msg += "\"" + node_file_name + "\".";
        throw OomphLibError(
          error_msg, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Read number of nodes
      node_file >> triangle_io.numberofpoints;
      unsigned n_node = static_cast<unsigned>(triangle_io.numberofpoints);

      // Spatial dimension of nodes
      unsigned dimension;
      node_file >> dimension;

#ifdef PARANOID
      if (dimension != 2)
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

      // Allocate storage
      triangle_io.pointlist =
        (double*)malloc(triangle_io.numberofpoints * 2 * sizeof(double));
      triangle_io.pointattributelist =
        (double*)malloc(triangle_io.numberofpoints *
                        triangle_io.numberofpointattributes * sizeof(double));
      if (boundary_markers_flag)
      {
        triangle_io.pointmarkerlist =
          (int*)malloc(triangle_io.numberofpoints * sizeof(int));
      }

      // Dummy for node number
      unsigned dummy_node_number;

      // Reset counter
      counter = 0;
      // Load in nodal posititions, point attributes
      // and boundary markers
      for (unsigned i = 0; i < n_node; i++)
      {
        node_file >> dummy_node_number;
        node_file >> triangle_io.pointlist[2 * i];
        node_file >> triangle_io.pointlist[2 * i + 1];
        for (unsigned j = 0; j < n_point_attributes; ++j)
        {
          node_file >> triangle_io.pointattributelist[counter];
          ++counter;
        }
        if (boundary_markers_flag)
        {
          node_file >> triangle_io.pointmarkerlist[i];
        }
      }
      node_file.close();


      // Process poly file to extract edges
      //-----------------------------------

      // Open poly file
      std::ifstream poly_file(poly_file_name.c_str(), std::ios_base::in);

      // Check that the file actually opened correctly
      if (!poly_file.is_open())
      {
        std::string error_msg("Failed to open poly file: ");
        error_msg += "\"" + poly_file_name + "\".";
        throw OomphLibError(
          error_msg, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
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
      for (unsigned i = 0; i < n_node_poly; i++)
      {
        // Read in (and discard) node number and x and y coordinates
        poly_file >> dummy;
        poly_file >> dummy;
        poly_file >> dummy;
        // read in the attributes
        for (unsigned j = 0; j < attribute_flag; ++j)
        {
          poly_file >> dummy;
        }
        // read in the boundary marker
        if (boundary_markers_flag == 1)
        {
          poly_file >> dummy;
        }
      }

      // Now extract the segment information
      //------------------------------------

      // Number of segments
      poly_file >> triangle_io.numberofsegments;
      unsigned n_segment = static_cast<unsigned>(triangle_io.numberofsegments);

      // Boundary marker flag
      poly_file >> boundary_markers_flag;

      // Allocate storage
      triangle_io.segmentlist =
        (int*)malloc(triangle_io.numberofsegments * 2 * sizeof(int));
      if (boundary_markers_flag)
      {
        triangle_io.segmentmarkerlist =
          (int*)malloc(triangle_io.numberofsegments * sizeof(int));
      }

      // Dummy for global segment number
      unsigned dummy_segment_number;

      // Extract information for each segment
      for (unsigned i = 0; i < n_segment; i++)
      {
        poly_file >> dummy_segment_number;
        poly_file >> triangle_io.segmentlist[2 * i];
        poly_file >> triangle_io.segmentlist[2 * i + 1];
        if (boundary_markers_flag)
        {
          poly_file >> triangle_io.segmentmarkerlist[i];
        }
      }

      // Extract hole center information
      poly_file >> triangle_io.numberofholes;
      unsigned n_hole = static_cast<unsigned>(triangle_io.numberofholes);

      // Allocate memory
      triangle_io.holelist =
        (double*)malloc(triangle_io.numberofholes * 2 * sizeof(double));


      // Dummy for hole number
      unsigned dummy_hole;
      // Loop over the holes to get centre coords
      for (unsigned ihole = 0; ihole < n_hole; ihole++)
      {
        // Read the centre value
        poly_file >> dummy_hole;
        poly_file >> triangle_io.holelist[2 * ihole];
        poly_file >> triangle_io.holelist[2 * ihole + 1];
      }

      // Extract regions information
      poly_file >> triangle_io.numberofregions;
      unsigned n_regions = static_cast<unsigned>(triangle_io.numberofregions);

      // Allocate memory
      triangle_io.regionlist =
        (double*)malloc(triangle_io.numberofregions * 4 * sizeof(double));

      // Check for using regions
      if (n_regions > 0)
      {
        use_attributes = true;
      }

      // Dummy for regions number
      unsigned dummy_region;

      // Loop over the regions to get their coords
      for (unsigned iregion = 0; iregion < n_regions; iregion++)
      {
        // Read the regions coordinates
        poly_file >> dummy_region;
        poly_file >> triangle_io.regionlist[4 * iregion];
        poly_file >> triangle_io.regionlist[4 * iregion + 1];
        poly_file >> triangle_io.regionlist[4 * iregion + 2];
        triangle_io.regionlist[4 * iregion + 3] = 0.0;

        // Skip the rest of the line, there is no need to read the size of
        // the elements in the region since that value is no longer used
        poly_file.ignore(80, '\n');

        // Verify if not using the default region number (zero)
        if (triangle_io.regionlist[4 * iregion + 2] == 0)
        {
          std::ostringstream error_message;
          error_message
            << "Please use another region id different from zero.\n"
            << "It is internally used as the default region number.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      poly_file.close();
    }


    /// Write all the triangulateio data to disk in a dump file
    /// that can then be used to restart simulations
    void dump_triangulateio(TriangulateIO& triangle_io, std::ostream& dump_file)
    {
      // Dump the triangles first
      dump_file << triangle_io.numberoftriangles
                << " # number of elements in TriangulateIO" << std::endl;

      dump_file << triangle_io.numberofcorners
                << " # number of nodes in each triangle" << std::endl;

      dump_file << triangle_io.numberoftriangleattributes
                << " # number of triangle attributes" << std::endl;

      // Loop over and dump the triangle information
      const int n_element = triangle_io.numberoftriangles;
      const int n_local_node = triangle_io.numberofcorners;
      const int n_attribute = triangle_io.numberoftriangleattributes;
      unsigned counter = 0, counter2 = 0;
      for (int e = 0; e < n_element; ++e)
      {
        // Dump the corners
        dump_file << e << " # element number " << std::endl;
        for (int n = 0; n < n_local_node; ++n)
        {
          dump_file << triangle_io.trianglelist[counter] << std::endl;
          ++counter;
        }
        // Dump the attributes
        dump_file << n_attribute << " # number of attributes" << std::endl;
        for (int n = 0; n < n_attribute; ++n)
        {
          dump_file << triangle_io.triangleattributelist[counter2] << std::endl;
          ++counter2;
        }
      }


      // Dump the points (nodes) next
      dump_file << triangle_io.numberofpoints
                << " # number of points in TriangulateIO" << std::endl;
      dump_file << triangle_io.numberofpointattributes
                << " # number of point attributes" << std::endl;
      // Test whether there are point markers
      bool point_marker_flag = true;
      if (triangle_io.pointmarkerlist == NULL)
      {
        point_marker_flag = false;
      }
      dump_file << point_marker_flag << " # point marker flag" << std::endl;


      // Now output the point data
      const int n_nodes = triangle_io.numberofpoints;
      const int n_point_attributes = triangle_io.numberofpointattributes;
      counter = 0;
      counter2 = 0;
      for (int n = 0; n < n_nodes; ++n)
      {
        dump_file << n << " # point number " << std::endl;
        for (int i = 0; i < 2; ++i)
        {
          dump_file << triangle_io.pointlist[counter] << std::endl;
          ++counter;
        }
        dump_file << n_point_attributes << " # number of point attributes "
                  << std::endl;
        // Output any attributes
        for (int i = 0; i < n_point_attributes; ++i)
        {
          dump_file << triangle_io.pointattributelist[counter2] << std::endl;
          ++counter2;
        }
        dump_file << point_marker_flag << " # point marker flag " << std::endl;
        // Output the boundary marker
        if (point_marker_flag)
        {
          dump_file << triangle_io.pointmarkerlist[n] << std::endl;
        }
      }

      // Now move onto the segments
      dump_file << triangle_io.numberofsegments
                << " # Number of segments in TriangulateIO " << std::endl;

      // Determine whether there are segment markers
      bool seg_marker_flag = true;
      if (triangle_io.segmentmarkerlist == NULL)
      {
        seg_marker_flag = false;
      }
      // Output this info in the file
      dump_file << seg_marker_flag << " # segment marker flag " << std::endl;

      const int n_segments = triangle_io.numberofsegments;
      counter = 0;
      // Now output the segment data
      for (int n = 0; n < n_segments; ++n)
      {
        dump_file << n << " # segment number " << std::endl;
        for (int i = 0; i < 2; ++i)
        {
          dump_file << triangle_io.segmentlist[counter] << std::endl;
          ++counter;
        }

        // If there is a boundary marker output
        dump_file << seg_marker_flag << " # segment marker flag " << std::endl;
        if (seg_marker_flag)
        {
          dump_file << triangle_io.segmentmarkerlist[n] << std::endl;
        }
      }

      // Now output the number of holes
      dump_file << triangle_io.numberofholes << " # number of holes "
                << std::endl;
      const int n_hole = triangle_io.numberofholes;
      // Output the hole data
      for (int h = 0; h < n_hole; ++h)
      {
        dump_file << h << " # hole number " << std::endl;
        dump_file << triangle_io.holelist[2 * h] << std::endl;
        dump_file << triangle_io.holelist[2 * h + 1] << std::endl;
      }

      // Now output the number of regions
      dump_file << triangle_io.numberofregions << " # number of regions "
                << std::endl;

      const int n_region = triangle_io.numberofregions;
      // Loop over the regions
      counter = 0;
      for (int r = 0; r < n_region; ++r)
      {
        dump_file << r << " # region number " << std::endl;
        for (unsigned i = 0; i < 4; i++)
        {
          dump_file << triangle_io.regionlist[counter] << std::endl;
          ++counter;
        }
      }
    }

    /// Read the triangulateio data from a dump file on
    /// disk, which can then be used to restart simulations
    void read_triangulateio(std::istream& restart_file,
                            TriangulateIO& triangle_io)
    {
      // String for reading
      std::string input_string;

      // Initialise the triangulate data structure
      initialise_triangulateio(triangle_io);

      // Read the first line up to termination sign
      getline(restart_file, input_string, '#');
      // Ignore the rest of the line
      restart_file.ignore(80, '\n');
      // Convert the number
      triangle_io.numberoftriangles = atoi(input_string.c_str());

      // Read the next line up to termination sign
      getline(restart_file, input_string, '#');
      // Ignore the rest of the line
      restart_file.ignore(80, '\n');
      // Convert the number
      triangle_io.numberofcorners = atoi(input_string.c_str());

      // Read the next line up to termination sign
      getline(restart_file, input_string, '#');
      // Ignore the rest of the line
      restart_file.ignore(80, '\n');
      // Convert the number
      triangle_io.numberoftriangleattributes = atoi(input_string.c_str());

      // Convert numbers into register variables
      const int n_element = triangle_io.numberoftriangles;
      const int n_local_node = triangle_io.numberofcorners;
      const int n_attribute = triangle_io.numberoftriangleattributes;

      // Allocate storage in the data structure
      triangle_io.trianglelist =
        (int*)malloc(triangle_io.numberoftriangles *
                     triangle_io.numberofcorners * sizeof(int));

      if (n_attribute > 0)
      {
        triangle_io.triangleattributelist = (double*)malloc(
          triangle_io.numberoftriangles *
          triangle_io.numberoftriangleattributes * sizeof(double));
      }

      // Loop over elements and load in data
      unsigned counter = 0, counter2 = 0;
      for (int e = 0; e < n_element; ++e)
      {
        // Read the next line and ignore it
        getline(restart_file, input_string);
        for (int n = 0; n < n_local_node; ++n)
        {
          getline(restart_file, input_string);
          triangle_io.trianglelist[counter] = atoi(input_string.c_str());
          ++counter;
        }
        // Read the attributes
        getline(restart_file, input_string);
        for (int n = 0; n < n_attribute; ++n)
        {
          getline(restart_file, input_string);
          triangle_io.triangleattributelist[counter2] =
            atof(input_string.c_str());
          ++counter2;
        }
      }


      // Read the points (nodes) next up to termination string
      getline(restart_file, input_string, '#');
      // ignore the rest
      restart_file.ignore(80, '\n');
      triangle_io.numberofpoints = atoi(input_string.c_str());

      // Read the point attributes next up to termination string
      getline(restart_file, input_string, '#');
      // ignore the rest
      restart_file.ignore(80, '\n');
      triangle_io.numberofpointattributes = atoi(input_string.c_str());

      // Read whether there are point markers
      getline(restart_file, input_string, '#');
      // ignore the rest
      restart_file.ignore(80, '\n');
      int point_marker_flag = atoi(input_string.c_str());

      // Allocate storage
      triangle_io.pointlist =
        (double*)malloc(triangle_io.numberofpoints * 2 * sizeof(double));
      triangle_io.pointattributelist =
        (double*)malloc(triangle_io.numberofpoints *
                        triangle_io.numberofpointattributes * sizeof(double));
      if (point_marker_flag)
      {
        triangle_io.pointmarkerlist =
          (int*)malloc(triangle_io.numberofpoints * sizeof(int));
      }


      // Now read the point data
      const int n_nodes = triangle_io.numberofpoints;
      const int n_point_attributes = triangle_io.numberofpointattributes;
      counter = 0;
      counter2 = 0;
      for (int n = 0; n < n_nodes; ++n)
      {
        // Ignore the first line
        getline(restart_file, input_string);
        // Get the positions
        for (int i = 0; i < 2; ++i)
        {
          getline(restart_file, input_string);
          triangle_io.pointlist[counter] = atof(input_string.c_str());
          ++counter;
        }

        // Ignore the next line about point attributes
        getline(restart_file, input_string);

        // Read any attributes
        for (int i = 0; i < n_point_attributes; ++i)
        {
          getline(restart_file, input_string);
          triangle_io.pointattributelist[counter2] = atof(input_string.c_str());
          ++counter2;
        }

        // Ignore the next line
        getline(restart_file, input_string);
        // Output the boundary marker
        if (point_marker_flag)
        {
          getline(restart_file, input_string);
          triangle_io.pointmarkerlist[n] = atoi(input_string.c_str());
        }
      }

      // Next read the segments
      getline(restart_file, input_string, '#');
      restart_file.ignore(80, '\n');
      triangle_io.numberofsegments = atoi(input_string.c_str());

      // Determine whether there are segment markers
      getline(restart_file, input_string, '#');
      // ignore the rest
      restart_file.ignore(80, '\n');
      int seg_marker_flag = atoi(input_string.c_str());

      // Allocate storage
      triangle_io.segmentlist =
        (int*)malloc(triangle_io.numberofsegments * 2 * sizeof(int));
      if (seg_marker_flag)
      {
        triangle_io.segmentmarkerlist =
          (int*)malloc(triangle_io.numberofsegments * sizeof(int));
      }

      const int n_segments = triangle_io.numberofsegments;
      counter = 0;
      // Now output the segment data
      for (int n = 0; n < n_segments; ++n)
      {
        getline(restart_file, input_string);
        // get input
        for (int i = 0; i < 2; ++i)
        {
          getline(restart_file, input_string);
          triangle_io.segmentlist[counter] = atoi(input_string.c_str());
          ++counter;
        }

        // Ignore the next line
        getline(restart_file, input_string);
        if (seg_marker_flag)
        {
          getline(restart_file, input_string);
          triangle_io.segmentmarkerlist[n] = atoi(input_string.c_str());
        }
      }

      // Now read the holes
      getline(restart_file, input_string, '#');
      restart_file.ignore(80, '\n');
      triangle_io.numberofholes = atoi(input_string.c_str());

      // Allocate memory
      triangle_io.holelist =
        (double*)malloc(triangle_io.numberofholes * 2 * sizeof(double));

      const int n_hole = triangle_io.numberofholes;
      // Output the hole data
      for (int h = 0; h < n_hole; ++h)
      {
        // Ignore the first line
        getline(restart_file, input_string);
        // get the centre
        getline(restart_file, input_string);
        triangle_io.holelist[2 * h] = atof(input_string.c_str());
        getline(restart_file, input_string);
        triangle_io.holelist[2 * h + 1] = atof(input_string.c_str());
      }

      // Now read the number of regions
      getline(restart_file, input_string, '#');
      restart_file.ignore(80, '\n');
      triangle_io.numberofregions = atoi(input_string.c_str());

      const int n_region = triangle_io.numberofregions;

      // Allocate memory
      triangle_io.regionlist =
        (double*)malloc(triangle_io.numberofregions * 4 * sizeof(double));

      // Loop over the regions
      counter = 0;
      for (int r = 0; r < n_region; ++r)
      {
        getline(restart_file, input_string);
        for (unsigned i = 0; i < 4; i++)
        {
          getline(restart_file, input_string);
          triangle_io.regionlist[counter] = atof(input_string.c_str());
          ++counter;
        }
      }
    }
  } // namespace TriangleHelper

#endif


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////


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
    double Tolerable_error = 1.0e-14;

  } // namespace ToleranceForVertexMismatchInPolygons


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////


  // =======================================================================
  // Connects the initial vertex of the curve section to a desired
  /// target polyline by specifying the vertex number. There is a checking
  /// which verifies that the initial vertex is close enough to the
  /// destination vertex on the target polyline by no more than the specified
  /// tolerance
  // =======================================================================
  void TriangleMeshCurveSection::connect_initial_vertex_to_polyline(
    TriangleMeshPolyLine* polyline_pt,
    const unsigned& vertex_number,
    const double& tolerance_for_connection)
  {
#ifdef PARANOID
    unsigned n_vertices = polyline_pt->nvertex();

    if (n_vertices <= vertex_number)
    {
      std::ostringstream error_stream;
      error_stream << "The vertex number you provided (" << vertex_number
                   << ") is greater\n than the number of vertices ("
                   << n_vertices << "in the specified TriangleMeshPolyLine.\n"
                   << "Remember that the vertex index starts at 0" << std::endl
                   << "Source boundary (" << boundary_id()
                   << ") wants to connect "
                   << "to destination boundary (" << polyline_pt->boundary_id()
                   << ")" << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Verify if there is really a match point in the specified
    // connection values
    Vector<double> v_src(2);
    Vector<double> v_dst(2);

    this->initial_vertex_coordinate(v_src);
    v_dst = polyline_pt->vertex_coordinate(vertex_number);

    double error = sqrt((v_src[0] - v_dst[0]) * (v_src[0] - v_dst[0]) +
                        (v_src[1] - v_dst[1]) * (v_src[1] - v_dst[1]));

    if (error > tolerance_for_connection)
    {
      std::ostringstream error_stream;
      error_stream << "The associated vertices for the connection"
                   << "\nare not close enough. Their respective values are:\n"
                   << "Source boundary id:(" << this->boundary_id() << ")\n"
                   << "Source vertex x:(" << v_src[0] << ") y:(" << v_src[1]
                   << ")\n"
                   << "Destination boundary id:(" << polyline_pt->boundary_id()
                   << ")"
                   << "\nAssociated vertex x:(" << v_dst[0] << ") y:("
                   << v_dst[1] << ")"
                   << "\nThe corresponding distance is: (" << error
                   << ") but the "
                   << "allowed\ntolerance is: (" << tolerance_for_connection
                   << ")" << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

#endif

    Initial_vertex_connected = true;
    Initial_vertex_connected_bnd_id = polyline_pt->boundary_id();
    Initial_vertex_connected_n_vertex = vertex_number;
    Initial_vertex_connected_n_chunk = polyline_pt->boundary_chunk();
  }

  // =======================================================================
  // Connects the final vertex of the curve section to a desired
  /// target polyline by specifying the vertex number. There is a checking
  /// which verifies that the final vertex is close enough to the
  /// destination vertex on the target polyline by no more than the specified
  /// tolerance
  // =======================================================================
  void TriangleMeshCurveSection::connect_final_vertex_to_polyline(
    TriangleMeshPolyLine* polyline_pt,
    const unsigned& vertex_number,
    const double& tolerance_for_connection)
  {
#ifdef PARANOID
    unsigned n_vertices = polyline_pt->nvertex();

    if (n_vertices <= vertex_number)
    {
      std::ostringstream error_stream;
      error_stream << "The vertex number you provided (" << vertex_number
                   << ") is greater\n than the number of vertices ("
                   << n_vertices << "in the specified TriangleMeshPolyLine.\n"
                   << "Remember that the vertex index starts at 0" << std::endl
                   << "Source boundary (" << boundary_id()
                   << ") wants to connect "
                   << "to destination boundary (" << polyline_pt->boundary_id()
                   << ")" << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Verify if there is really a match point in the specified
    // connection values
    Vector<double> v_src(2);
    Vector<double> v_dst(2);

    this->final_vertex_coordinate(v_src);
    v_dst = polyline_pt->vertex_coordinate(vertex_number);

    double error = sqrt((v_src[0] - v_dst[0]) * (v_src[0] - v_dst[0]) +
                        (v_src[1] - v_dst[1]) * (v_src[1] - v_dst[1]));

    if (error > tolerance_for_connection)
    {
      std::ostringstream error_stream;
      error_stream << "The associated vertices for the connection"
                   << "\nare not close enough. Their respective values are:\n"
                   << "Source boundary id:(" << this->boundary_id() << ")\n"
                   << "Source vertex x:(" << v_src[0] << ") y:(" << v_src[1]
                   << ")\n"
                   << "Destination boundary id:(" << polyline_pt->boundary_id()
                   << ")"
                   << "\nAssociated vertex x:(" << v_dst[0] << ") y:("
                   << v_dst[1] << ")"
                   << "\nThe corresponding distance is: (" << error
                   << ") but the "
                   << "allowed\ntolerance is: (" << tolerance_for_connection
                   << ")" << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

#endif

    Final_vertex_connected = true;
    Final_vertex_connected_bnd_id = polyline_pt->boundary_id();
    Final_vertex_connected_n_vertex = vertex_number;
    Final_vertex_connected_n_chunk = polyline_pt->boundary_chunk();
  }

  // =======================================================================
  // Connects the initial vertex of the curve section to a desired
  /// target curviline by specifying the s value (intrinsic value on the
  /// geometric object of the curviline) where to connect on the target
  /// curviline. There is a checking which verifies that the initial vertex
  /// and the coordinates on the given s value are close enough by no more
  /// than the given tolerance
  // =======================================================================
  void TriangleMeshCurveSection::connect_initial_vertex_to_curviline(
    TriangleMeshCurviLine* curviline_pt,
    const double& s_value,
    const double& tolerance_for_connection)
  {
#ifdef PARANOID
    double z_initial = curviline_pt->zeta_start();
    double z_final = curviline_pt->zeta_end();
    double z_max = std::max(z_initial, z_final);
    double z_min = std::min(z_initial, z_final);
    if (s_value < z_min || z_max < s_value)
    {
      std::ostringstream error_stream;
      error_stream << "The s value you provided for connection (" << s_value
                   << ") is out\nof the limits of the specified "
                   << "TriangleMeshCurviLine.\nThe limits are [" << z_initial
                   << ", " << z_final << "]" << std::endl
                   << "Source boundary (" << boundary_id()
                   << ") wants to connect "
                   << "to destination boundary (" << curviline_pt->boundary_id()
                   << ")" << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Verify if there is really a match point in the specified
    // connection values
    Vector<double> v_src(2);
    Vector<double> v_dst(2);
    Vector<double> z(1);

    this->initial_vertex_coordinate(v_src);
    z[0] = s_value;
    curviline_pt->geom_object_pt()->position(z, v_dst);
    double error = sqrt((v_src[0] - v_dst[0]) * (v_src[0] - v_dst[0]) +
                        (v_src[1] - v_dst[1]) * (v_src[1] - v_dst[1]));
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
        << "\nThe corresponding distance is: (" << error << ") but the "
        << "allowed\ntolerance is: (" << tolerance_for_connection << ")"
        << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
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
  // Connects the final vertex of the curve section to a desired
  /// target curviline by specifying the s value (intrinsic value on the
  /// geometric object of the curviline) where to connect on the target
  /// curviline. There is a checking which verifies that the final vertex
  /// and the coordinates on the given s value are close enough by no more
  /// than the given tolerance
  // =======================================================================
  void TriangleMeshCurveSection::connect_final_vertex_to_curviline(
    TriangleMeshCurviLine* curviline_pt,
    const double& s_value,
    const double& tolerance_for_connection)
  {
#ifdef PARANOID
    double z_initial = curviline_pt->zeta_start();
    double z_final = curviline_pt->zeta_end();
    double z_max = std::max(z_initial, z_final);
    double z_min = std::min(z_initial, z_final);
    if (s_value < z_min || z_max < s_value)
    {
      std::ostringstream error_stream;
      error_stream << "The s value you provided for connection (" << s_value
                   << ") is out\nof the limits of the specified "
                   << "TriangleMeshCurviLine.\nThe limits are [" << z_initial
                   << ", " << z_final << "]" << std::endl
                   << "Source boundary (" << boundary_id()
                   << ") wants to connect "
                   << "to destination boundary (" << curviline_pt->boundary_id()
                   << ")" << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Verify if there is really a match point in the specified
    // connection values
    Vector<double> v_src(2);
    Vector<double> v_dst(2);
    Vector<double> z(1);

    this->final_vertex_coordinate(v_src);
    z[0] = s_value;
    curviline_pt->geom_object_pt()->position(z, v_dst);

    double error = sqrt((v_src[0] - v_dst[0]) * (v_src[0] - v_dst[0]) +
                        (v_src[1] - v_dst[1]) * (v_src[1] - v_dst[1]));

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
        << "\nThe corresponding distance is: (" << error << ") but the "
        << "allowed\ntolerance is: (" << tolerance_for_connection << ")"
        << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
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


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //=====================================================================
  /// Class defining a closed curve for the Triangle mesh generation
  //=====================================================================
  TriangleMeshClosedCurve::TriangleMeshClosedCurve(
    const Vector<TriangleMeshCurveSection*>& curve_section_pt,
    const Vector<double>& internal_point_pt,
    const bool& is_internal_point_fixed)
    : TriangleMeshCurve(curve_section_pt),
      Internal_point_pt(internal_point_pt),
      Is_internal_point_fixed(is_internal_point_fixed)
  {
    // Matching of curve sections i.e. the last vertex of the i curve
    // section should match with the first vertex of the i+1 curve
    // section

    // Total number of boundaries
    const unsigned n_boundaries = Curve_section_pt.size();

    // Need at least two
    if (n_boundaries < 2)
    {
      std::ostringstream error_stream;
      error_stream << "Sorry -- I'm afraid we insist that a closed curve is\n"
                   << "specified by at least two separate CurveSections.\n"
                   << "You've only specified (" << n_boundaries << ")"
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check last point of each boundary bit coincides with first point
    // on next one
    for (unsigned i = 0; i < n_boundaries - 1; i++)
    {
      // Auxiliary vertex for storing the vertex values of contiguous curves
      Vector<double> v1(2);

      // This is for getting the final coordinates of the i curve section
      curve_section_pt[i]->final_vertex_coordinate(v1);

      // Auxiliary vertex for storing the vertex values of contiguous curves
      Vector<double> v2(2);

      // This is for the start coordinates of the i+1 curve section
      curve_section_pt[i + 1]->initial_vertex_coordinate(v2);

      // Work out error
      double error = sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2));

      if (error > ToleranceForVertexMismatchInPolygons::Tolerable_error)
      {
        std::ostringstream error_stream;
        error_stream
          << "The start and end points of curve section boundary parts\n"
          << i << " and " << i + 1
          << " don't match when judged with the tolerance of "
          << ToleranceForVertexMismatchInPolygons::Tolerable_error
          << " which\nis specified in the namespace variable\n"
          << "ToleranceForVertexMismatchInPolygons::Tolerable_error.\n\n"
          << "These are the vertices coordinates:\n"
          << "Curve section (" << i << ") final vertex: (" << v1[0] << ", "
          << v1[1] << ")\n"
          << "Curve section (" << i + 1 << ") initial vertex: (" << v2[0]
          << ", " << v2[1] << ")\n"
          << "The distance between the vertices is (" << error << ")\n"
          << "Feel free to adjust this or to recompile the code without\n"
          << "paranoia if you think this is OK...\n"
          << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        // Aligns (only implemented for polylines)
        TriangleMeshPolyLine* current_polyline =
          dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i]);
        TriangleMeshPolyLine* next_polyline =
          dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i + 1]);

        // Was it able to do the cast?
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
    Curve_section_pt[n_boundaries - 1]->final_vertex_coordinate(v2);

    double error = sqrt(pow(v2[0] - v1[0], 2) + pow(v2[1] - v1[1], 2));

    if (error > ToleranceForVertexMismatchInPolygons::Tolerable_error)
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
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      // Aligns (only implemented for polylines)
      TriangleMeshPolyLine* first_polyline =
        dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[0]);
      TriangleMeshPolyLine* last_polyline =
        dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[n_boundaries - 1]);

      // Was it able to do the cast?
      if (first_polyline && last_polyline)
      {
        unsigned last_vertex = last_polyline->nvertex() - 1;
        first_polyline->vertex_coordinate(0) =
          last_polyline->vertex_coordinate(last_vertex);
      }
    }
  }


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


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
    const Vector<double>& internal_point_pt,
    const bool& is_internal_point_fixed)
    : TriangleMeshCurve(boundary_polyline_pt),
      TriangleMeshClosedCurve(
        boundary_polyline_pt, internal_point_pt, is_internal_point_fixed),
      Enable_redistribution_of_segments_between_polylines(false),
      Can_update_configuration(false),
      Polygon_fixed(false)
  {
    // Get the number of polylines
    const unsigned n_bound = boundary_polyline_pt.size();

    // Check that all the constituent TriangleMeshCurveSection are
    // instance of TriangleMeshPolyLine
    for (unsigned p = 0; p < n_bound; p++)
    {
      TriangleMeshPolyLine* tmp_polyline_pt =
        dynamic_cast<TriangleMeshPolyLine*>(boundary_polyline_pt[p]);
      if (tmp_polyline_pt == 0)
      {
        std::ostringstream error_stream;
        error_stream << "The (" << p << ") TriangleMeshCurveSection is not a "
                     << "TriangleMeshPolyLine.\nThe TriangleMeshPolygon object"
                     << "is constituent of only TriangleMeshPolyLine objects.\n"
                     << "Verify that all the constituent elements of the "
                     << "TriangleMeshPolygon\nare instantiated as "
                     << "TriangleMeshPolyLines." << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Check that the polylines are contiguous
    bool contiguous = true;
    unsigned i_offensive = 0;

    // Need at least two
    if (n_bound < 2)
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
    for (unsigned i = 0; i < n_bound - 1; i++)
    {
      // Get vector last vertex in current polyline
      unsigned last_vertex = (polyline_pt(i)->nvertex()) - 1;
      Vector<double> v1 = polyline_pt(i)->vertex_coordinate(last_vertex);

      // Get vector to first vertex in next polyline
      Vector<double> v2 = polyline_pt(i + 1)->vertex_coordinate(0);

      // Work out error
      double error = sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2));

      // Is error accetable?
      if (error > ToleranceForVertexMismatchInPolygons::Tolerable_error)
      {
        contiguous = false;
        i_offensive = i;
        break;
      }
      // Align
      else
      {
        polyline_pt(i + 1)->vertex_coordinate(0) =
          polyline_pt(i)->vertex_coordinate(last_vertex);
      }
    }

    // Does the last one connect to the first one?

    // Get vector last vertex last polyline
    unsigned last_vertex = (polyline_pt(n_bound - 1)->nvertex()) - 1;
    Vector<double> v1 =
      polyline_pt(n_bound - 1)->vertex_coordinate(last_vertex);

    // Get vector first vertex first polyline
    Vector<double> v2 = polyline_pt(0)->vertex_coordinate(0);
    double error = sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2));

    if (error > ToleranceForVertexMismatchInPolygons::Tolerable_error)
    {
      contiguous = false;
      i_offensive = n_bound - 1;
    }
    else
    {
      polyline_pt(0)->vertex_coordinate(0) =
        polyline_pt(n_bound - 1)->vertex_coordinate(last_vertex);
    }

    if (!contiguous)
    {
      std::ostringstream error_stream;
      error_stream
        << "The polylines specified \n"
        << "should define a closed geometry, i.e. the first/last vertex of\n"
        << "adjacent polylines should match.\n\n"
        << "Your polyline number " << i_offensive
        << " has no contiguous neighbour, when judged \nwith the tolerance of "
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
      Vector<Vector<double>> polygon_vertex;

      // Total number of vertices
      unsigned nvertex = 0;

      // Storage for first/last point on polyline for contiguousness check
      Vector<double> last_vertex(2);
      Vector<double> first_vertex(2);

      // Get vertices
      unsigned npolyline = boundary_polyline_pt.size();
      for (unsigned i = 0; i < npolyline; i++)
      {
        // Number of vertices
        unsigned nvert = boundary_polyline_pt[i]->nvertex();
        for (unsigned j = 0; j < nvert; j++)
        {
          // Check contiguousness
          if ((i > 1) && (j == 0))
          {
            first_vertex = polyline_pt(i)->vertex_coordinate(j);
            double dist = sqrt(pow(first_vertex[0] - last_vertex[0], 2) +
                               pow(first_vertex[1] - last_vertex[1], 2));
            if (dist > ToleranceForVertexMismatchInPolygons::Tolerable_error)
            {
              std::ostringstream error_stream;
              error_stream
                << "The start and end points of polylines " << i << " and "
                << i + 1 << " don't match when judged\n"
                << "with the tolerance ("
                << ToleranceForVertexMismatchInPolygons::Tolerable_error
                << ") which is specified in the namespace \nvariable "
                << "ToleranceForVertexMismatchInPolygons::"
                << "Tolerable_error.\n\n"
                << "Feel free to adjust this or to recompile the "
                << "code without\n"
                << "paranoia if you think this is OK...\n"
                << std::endl;
              throw OomphLibError(error_stream.str(),
                                  "TriangleMeshPolygon::TriangleMeshPolygon()",
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
          // Get vertex (ignore end point)
          if (j < nvert - 1)
          {
            polygon_vertex.push_back(polyline_pt(i)->vertex_coordinate(j));
          }
          // Prepare for check of contiguousness
          else
          {
            last_vertex = polyline_pt(i)->vertex_coordinate(j);
          }
        }
      }

      // Total number of vertices
      nvertex = polygon_vertex.size();

      // Counter for number of intersections
      unsigned intersect_counter = 0;

      // Get first vertex
      Vector<double> p1 = polygon_vertex[0];
      for (unsigned i = 1; i <= nvertex; i++)
      {
        // Get second vertex by wrap-around
        Vector<double> p2 = polygon_vertex[i % nvertex];

        if (Internal_point_pt[1] > std::min(p1[1], p2[1]))
        {
          if (Internal_point_pt[1] <= std::max(p1[1], p2[1]))
          {
            if (Internal_point_pt[0] <= std::max(p1[0], p2[0]))
            {
              if (p1[1] != p2[1])
              {
                double xintersect = (Internal_point_pt[1] - p1[1]) *
                                      (p2[0] - p1[0]) / (p2[1] - p1[1]) +
                                    p1[0];
                if ((p1[0] == p2[0]) || (Internal_point_pt[0] <= xintersect))
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
      if (intersect_counter % 2 == 0)
      {
        std::ostringstream error_stream;
        error_stream
          << "The internal point at " << Internal_point_pt[0] << " "
          << Internal_point_pt[1]
          << " isn't in the polygon that describes the internal closed "
          << "curve!\nPolygon vertices are at: \n";
        for (unsigned i = 0; i < nvertex; i++)
        {
          error_stream << polygon_vertex[i][0] << " " << polygon_vertex[i][1]
                       << "\n";
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
        throw OomphLibError(error_stream.str(),
                            "TriangleMeshPolygon::TriangleMeshPolygon()",
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
  }


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //=====================================================================
  /// Class defining an open curve for the Triangle mesh generation
  //=====================================================================
  TriangleMeshOpenCurve::TriangleMeshOpenCurve(
    const Vector<TriangleMeshCurveSection*>& curve_section_pt)
    : TriangleMeshCurve(curve_section_pt)
  {
    // Matching of curve sections i.e. the last vertex of
    // the i curve section should match with the first
    // vertex of the i+1 curve section

    // Total number of boundaries
    unsigned n_boundaries = Curve_section_pt.size();

    // Check last point of each boundary bit coincides with first point
    // on next one
    for (unsigned i = 0; i < n_boundaries - 1; i++)
    {
      // Auxiliary vertex for storing the vertex values of contiguous curves
      Vector<double> v1(2);
      Vector<double> v2(2);

      // This is for getting the final coordinates of the i curve section
      Curve_section_pt[i]->final_vertex_coordinate(v1);

      // This is for the start coordinates of the i+1 curve section
      Curve_section_pt[i + 1]->initial_vertex_coordinate(v2);

      // Work out error
      double error = sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2));
      if (error > ToleranceForVertexMismatchInPolygons::Tolerable_error)
      {
        std::ostringstream error_stream;
        error_stream
          << "The start and end points of curve section boundary parts " << i
          << " and " << i + 1
          << " don't match when judged \nwith the tolerance of "
          << ToleranceForVertexMismatchInPolygons::Tolerable_error
          << " which is specified in the namespace \nvariable "
          << "ToleranceForVertexMismatchInPolygons::Tolerable_error.\n\n"
          << "Feel free to adjust this or to recompile the code without\n"
          << "paranoia if you think this is OK...\n"
          << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        // Aligns (only implemented for polylines)
        TriangleMeshPolyLine* current_polyline =
          dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i]);
        TriangleMeshPolyLine* next_polyline =
          dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i + 1]);

        if (current_polyline && next_polyline)
        {
          unsigned last_vertex = current_polyline->nvertex() - 1;
          next_polyline->vertex_coordinate(0) =
            current_polyline->vertex_coordinate(last_vertex);
        }
      }

    } // For n_boundaries - 1
  }

#ifdef OOMPH_HAS_TRIANGLE_LIB

  //======================================================================
  /// Create TriangulateIO object from outer boundaries, internal
  /// boundaries, and open curves. Add the holes and regions
  /// information as well
  //======================================================================
  void UnstructuredTwoDMeshGeometryBase::build_triangulateio(
    Vector<TriangleMeshPolygon*>& outer_polygons_pt,
    Vector<TriangleMeshPolygon*>& internal_polygons_pt,
    Vector<TriangleMeshOpenCurve*>& open_curves_pt,
    Vector<Vector<double>>& extra_holes_coordinates,
    std::map<unsigned, Vector<double>>& regions_coordinates,
    std::map<unsigned, double>& regions_areas,
    TriangulateIO& triangulate_io)
  {
    // These are the general stages of the algorithm
    // --------------------------------------------------------------------
    // 1) Create the boundary connections structure, every entry states
    // if the initial or final vertex is connected to other vertex
    // --------------------------------------------------------------------
    // 2) Create the structure for base vertices, the base vertices are
    // those that are not connected
    // --------------------------------------------------------------------
    // 3) Assign a unique id to each vertex in the polylines (creates a
    // look up scheme used to create the segments)
    // --------------------------------------------------------------------
    // 4) Create the segments using the unique id of the vertices and
    // assign a segment id to each segment (the one from the
    // polyline-boundary)
    // ------------------------------------------------------------------
    // 5) Fill the triangulateio containers with the collected
    // information
    // --------------------------------------------------------------

    // ------------------------------------------------------------------
    // 1st- stage

    // 1) Create the boundary connections structure for the outer
    // boundaries, internal boundaries and open curves

    // Store the number of vertices on each boundary (which is
    // represented by a polyline -- for quick access--). We may have
    // more than one polyline with the same boundary id, that is why we
    // need a vector to represent the set of polylines associated to a
    // boundary (the chunks). The mentioned case is quite common when
    // working in parallel, when a boundary may be split because of the
    // distribution strategy
    std::map<unsigned, std::map<unsigned, unsigned>> boundary_chunk_n_vertices;

    // Note: For each polyline, we only consider (v-1) vertices since
    // the first vertex of the p-th polyline (when p > 1) is the same as
    // the last vertex of the (p-1)-th polyline (when p > 1). KEEP THIS
    // ---ALWAYS--- IN MIND WHEN REVIEWING THE CODE

    // The connections matrix

    // Stores the vertex_connection_info:
    // - is_connected
    // - boundary_id_to_connect
    // - boundary_chunk_to_connect
    // - vertex_number_to_connect
    // of the initial and final vertex of each polyline
    // -----------------------
    // (boundary, chunk#, vertex (initial[0] or final[1])) ->
    // vertex_connection_info
    // -----------------------
    // map[x][][] = boundary_id
    // map[][x][] = chunk_id
    // Vector[][][x] = vertex#, only initial or final (that is why only
    // two indexes)
    std::map<unsigned, std::map<unsigned, Vector<vertex_connection_info>>>
      connection_matrix;

    // Initialize the base vertex structure (every vertex is a not base
    // vertex by default)

    // The data structure that stores the base vertex information
    // Stores the base_vertex_info:
    // - done
    // - base_vertex
    // - boundary_id
    // - boundary_chunk
    // - vertex_number
    // of the initial and final vertex of each polyline
    // -----------------------
    // (boundary, chunk#, vertex (initial[0] or final[1])) -> base_vertex_info
    // -----------------------
    // map[x][][] = boundary_id
    // map[][x][] = chunk_id
    // Vector[][][x] = vertex#, only initial or final (that is why only
    // two indexes)
    std::map<unsigned, std::map<unsigned, Vector<base_vertex_info>>>
      base_vertices;

    // Get the number of outer polygons
    const unsigned n_outer_polygons = outer_polygons_pt.size();

    for (unsigned i = 0; i < n_outer_polygons; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = outer_polygons_pt[i]->npolyline();

      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get a pointer to the current polyline
        TriangleMeshPolyLine* tmp_polyline_pt =
          outer_polygons_pt[i]->polyline_pt(p);

        // Get the boundary id of the current polyline
        const unsigned bound_id = tmp_polyline_pt->boundary_id();

        // The number of vertices in the current polyline
        const unsigned n_vertices = tmp_polyline_pt->nvertex();

        // Get the chunk number associated with this boundary
        const unsigned bound_chunk = tmp_polyline_pt->boundary_chunk();

        // Store the number of vertices of the polyline in the global
        // storage
        boundary_chunk_n_vertices[bound_id][bound_chunk] = n_vertices;

        // Get the next polyline (or the initial polyline)
        TriangleMeshPolyLine* next_polyline_pt = 0;

        // Is there next polyline
        if (p < n_polylines - 1)
        {
          // Set the next polyline
          next_polyline_pt = outer_polygons_pt[i]->polyline_pt(p + 1);
        }
        else
        {
          // The next polyline is the initial polyline
          next_polyline_pt = outer_polygons_pt[i]->polyline_pt(0);
        }

        // Add the information to the connections matrix
        add_connection_matrix_info_helper(
          tmp_polyline_pt, connection_matrix, next_polyline_pt);

        // Initialise the base vertices structure for the current
        // polyline
        initialise_base_vertex(tmp_polyline_pt, base_vertices);

      } // for (p < n_polylines)

    } // for (i < n_outer_polygons)

    // Get the number of internal polygons
    const unsigned n_internal_polygons = internal_polygons_pt.size();

    for (unsigned i = 0; i < n_internal_polygons; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = internal_polygons_pt[i]->npolyline();

      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get a pointer to the current polyline
        TriangleMeshPolyLine* tmp_polyline_pt =
          internal_polygons_pt[i]->polyline_pt(p);

        // Get the boundary id of the current polyline
        const unsigned bound_id = tmp_polyline_pt->boundary_id();

        // The number of vertices in the current polyline
        const unsigned n_vertices = tmp_polyline_pt->nvertex();

        // Get the chunk number associated with this boundary
        const unsigned bound_chunk = tmp_polyline_pt->boundary_chunk();

        // Store the number of vertices of the polyline in the global
        // storage
        boundary_chunk_n_vertices[bound_id][bound_chunk] = n_vertices;

        // Get the next polyline (or the initial polyline)
        TriangleMeshPolyLine* next_polyline_pt = 0;

        // Is there next polyline
        if (p < n_polylines - 1)
        {
          // Set the next polyline
          next_polyline_pt = internal_polygons_pt[i]->polyline_pt(p + 1);
        }
        else
        {
          // The next polyline is the initial polyline
          next_polyline_pt = internal_polygons_pt[i]->polyline_pt(0);
        }

        // Add the information to the connections matrix
        add_connection_matrix_info_helper(
          tmp_polyline_pt, connection_matrix, next_polyline_pt);

        // Initialise the base vertices structure for the current
        // polyline
        initialise_base_vertex(tmp_polyline_pt, base_vertices);

      } // for (p < n_polylines)

    } // for (i < n_internal_polygons)

    // Get the number of open curves
    const unsigned n_open_curves = open_curves_pt.size();

    for (unsigned i = 0; i < n_open_curves; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = open_curves_pt[i]->ncurve_section();

      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get a pointer to the current polyline
        TriangleMeshPolyLine* tmp_polyline_pt =
          open_curves_pt[i]->polyline_pt(p);

        // Get the boundary id of the current polyline
        const unsigned bound_id = tmp_polyline_pt->boundary_id();

        // The number of vertices in the current polyline
        const unsigned n_vertices = tmp_polyline_pt->nvertex();

        // Get the chunk number associated with this boundary
        const unsigned bound_chunk = tmp_polyline_pt->boundary_chunk();

        // Store the number of vertices of the polyline in the global
        // storage
        boundary_chunk_n_vertices[bound_id][bound_chunk] = n_vertices;

        // Get the next polyline (or the initial polyline)
        TriangleMeshPolyLine* next_polyline_pt = 0;

        // Is there next polyline
        if (p < n_polylines - 1)
        {
          // Set the next polyline
          next_polyline_pt = open_curves_pt[i]->polyline_pt(p + 1);
        }
        else
        {
          // If we are in the last polyline of the open curve there is
          // no actual next polyline
        }

        // Add the information to the connections matrix
        add_connection_matrix_info_helper(
          tmp_polyline_pt, connection_matrix, next_polyline_pt);

        // Initialise the base vertices structure for the current
        // polyline
        initialise_base_vertex(tmp_polyline_pt, base_vertices);

      } // for (p < n_polylines)

    } // for (i < n_open_curves)

    // ------------------------------------------------------------------

    // ------------------------------------------------------------------
    // 2) Create the structure for base vertices, the base vertices are
    // those that are not connected
    // ------------------------------------------------------------------

    // Loop over the polylines in the outer polygons and indentify the
    // base vertices
    for (unsigned i = 0; i < n_outer_polygons; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = outer_polygons_pt[i]->npolyline();

      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get a pointer to the current polyline
        TriangleMeshPolyLine* tmp_polyline_pt =
          outer_polygons_pt[i]->polyline_pt(p);

        // Identify the base vertices in the current polyline
        add_base_vertex_info_helper(tmp_polyline_pt,
                                    base_vertices,
                                    connection_matrix,
                                    boundary_chunk_n_vertices);

      } // for (p < n_polylines)

    } // for (i < n_outer_polygons)

    // Loop over the polylines in the internal polygons and indentify the
    // base vertices
    for (unsigned i = 0; i < n_internal_polygons; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = internal_polygons_pt[i]->npolyline();

      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get a pointer to the current polyline
        TriangleMeshPolyLine* tmp_polyline_pt =
          internal_polygons_pt[i]->polyline_pt(p);

        // Identify the base vertices in the current polyline
        add_base_vertex_info_helper(tmp_polyline_pt,
                                    base_vertices,
                                    connection_matrix,
                                    boundary_chunk_n_vertices);

      } // for (p < n_polylines)

    } // for (i < n_internal_polygons)

    // Loop over the polylines in the open curves and indentify the base
    // vertices
    for (unsigned i = 0; i < n_open_curves; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = open_curves_pt[i]->ncurve_section();

      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get a pointer to the current polyline
        TriangleMeshPolyLine* tmp_polyline_pt =
          open_curves_pt[i]->polyline_pt(p);

        // Identify the base vertices in the current polyline
        add_base_vertex_info_helper(tmp_polyline_pt,
                                    base_vertices,
                                    connection_matrix,
                                    boundary_chunk_n_vertices);

      } // for (p < n_polylines)

    } // for (i < n_open_curves)

    // ------------------------------------------------------------------

    // ------------------------------------------------------------------
    // 3) Assign a unique id to each vertex in the polylines (creates a
    // look up scheme used to create the segments)
    // ------------------------------------------------------------------

    // Create the storage for the look-up scheme
    // (boundary_local, chunk#, vertex#) -> global_vertex_id
    // map[x][][] = boundary
    // map[][x][] = chunk_id
    // Vector[][][x] = vertex#
    std::map<unsigned, std::map<unsigned, Vector<int>>>
      boundary_chunk_vertex_to_global_vertex_id;

    // Create an entry in the map for each boundary, then do the same
    // for the chunks and finally resize the container (Vector) to store
    // the vertices
    for (std::map<unsigned, std::map<unsigned, unsigned>>::iterator it =
           boundary_chunk_n_vertices.begin();
         it != boundary_chunk_n_vertices.end();
         it++)
    {
      // Get the boundary id
      const unsigned b = (*it).first;

      // Now loop over the chunks
      for (std::map<unsigned, unsigned>::iterator itt = (*it).second.begin();
           itt != (*it).second.end();
           itt++)
      {
        // Get the chunk id
        const unsigned c = (*itt).first;

        // Get the number of vertices associated with the boundary-chunk
        // and resize the container
        const unsigned n_vertices = boundary_chunk_n_vertices[b][c];

        // Now create storage in the container and resize the vector that
        // stores the vertices ids. Initialize the entries to -1
        boundary_chunk_vertex_to_global_vertex_id[b][c].resize(n_vertices, -1);

      } // Loop over the chunks

    } // Loop over the boundaries

    // Counter for the numbering of the global vertices
    unsigned global_vertex_number = 0;

    // Container for the vertices
    Vector<Vector<double>> global_vertices;

    // Go through all the vertices in the polylines and assign a unique
    // id only to the base vertices, any other vertex copy the unique id
    // from its base vertex

    // The total number of added vertices in the outer polygons
    unsigned n_vertices_outer_polygons = 0;

    // Start with the outer polygons
    for (unsigned i = 0; i < n_outer_polygons; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = outer_polygons_pt[i]->npolyline();

      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get a pointer to the current polyline
        TriangleMeshPolyLine* tmp_polyline_pt =
          outer_polygons_pt[i]->polyline_pt(p);

        // Get the boundary id of the current polyline
        const unsigned bound_id = tmp_polyline_pt->boundary_id();

        // The number of vertices in the current polyline
        const unsigned n_vertices = tmp_polyline_pt->nvertex();

        // Get the current chunk number of the polyline
        const unsigned bnd_chunk = tmp_polyline_pt->boundary_chunk();

        // Assign a global vertex id to the initial vertex
        // -----------------------------------------------

        // Get the base vertex information of the initial vertex
        base_vertex_info initial_base_vertex_info =
          base_vertices[bound_id][bnd_chunk][0];

#ifdef PARANOID
        if (!initial_base_vertex_info.has_base_vertex_assigned)
        {
          std::ostringstream error_message;
          std::string output_string =
            "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
          error_message
            << "The initial vertex of the current polyline has no base\n"
            << "vertex assigned\n"
            << "Outer polygon number: (" << i << ")\n\n"
            << "Polyline number: (" << p << ")\n"
            << "Boundary id: (" << bound_id << ")\n"
            << "Boundary chunk: (" << bnd_chunk << ")\n";
          throw OomphLibError(
            error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
        } // if (!initial_base_vertex_info.has_base_vertex_assigned)
#endif

        // The base vertex boundary id
        unsigned bvbi = initial_base_vertex_info.boundary_id;
        // The base vertex boundary chunkx
        unsigned bvbc = initial_base_vertex_info.boundary_chunk;
        // The vertex number of the base vertex
        unsigned bvvn = initial_base_vertex_info.vertex_number;

        // Get the global vertex id of the base vertex
        int global_vertex_id =
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn];

        // Check if the global vertex id has been already assigned
        if (global_vertex_id != -1)
        {
          // If that is the case then copy the global vertex id in the
          // current slot
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][0] =
            global_vertex_id;
        } // if (global_vertex_id != -1)
        else
        {
          // Assign a global vertex id to the base vertex
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn] =
            global_vertex_number;

          // ... and copy the value to the initial vertex
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][0] =
            global_vertex_number;

          // ... increase the counter for the "global_vertex_number"
          global_vertex_number++;

          // Get the vertex
          Vector<double> vertex = tmp_polyline_pt->vertex_coordinate(0);

          // ... and add it to the global vertex container
          global_vertices.push_back(vertex);
        }

        // Is the initial vertex a base vertex
        if (initial_base_vertex_info.is_base_vertex)
        {
          // Increase the independent vertices counter
          n_vertices_outer_polygons++;
        }

        // ------------------------------------------------------------
        // Now loop over the intermediate vertices and assign a unique
        // vertex id (all intermediate vertices are base vertices)
        for (unsigned v = 1; v < n_vertices - 1; v++)
        {
          // Get the global vertex id
          global_vertex_id =
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][v];

          // Check if the global vertex id has been already assigned.
          // We do nothing if it has been already assigned
          // (global_vertex_id!=-1).

          // If it has not been already assigned (global_vertex_id==-1)
          // then set a new global vertex number, and add the vertex to
          // the global vertices container
          if (global_vertex_id == -1)
          {
            // Set a value for the global vertex
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][v] =
              global_vertex_number;

            // ... increase the counter for the "global_vertex_number"
            global_vertex_number++;

            // Add the vertex to the global vertex container
            Vector<double> vertex = tmp_polyline_pt->vertex_coordinate(v);
            // Add the vertex to the global vertex container
            global_vertices.push_back(vertex);

          } // if (global_vertex_id == -1)

          // Increase the independent vertices counter
          n_vertices_outer_polygons++;

        } // for (n_vertices-1)

        // Assign a global vertex id to the final vertex
        // -----------------------------------------------

        // Get the base vertex information of the final vertex
        base_vertex_info final_base_vertex_info =
          base_vertices[bound_id][bnd_chunk][1];

#ifdef PARANOID
        if (!final_base_vertex_info.has_base_vertex_assigned)
        {
          std::ostringstream error_message;
          std::string output_string =
            "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
          error_message
            << "The final vertex of the current polyline has no base\n"
            << "vertex assigned\n"
            << "Outer polygon number: (" << i << ")\n\n"
            << "Polyline number: (" << p << ")\n"
            << "Boundary id: (" << bound_id << ")\n"
            << "Boundary chunk: (" << bnd_chunk << ")\n";
          throw OomphLibError(
            error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
        } // if (!final_base_vertex_info.has_base_vertex_assigned)
#endif

        // The base vertex boundary id
        bvbi = final_base_vertex_info.boundary_id;
        // The base vertex boundary chunkx
        bvbc = final_base_vertex_info.boundary_chunk;
        // The vertex number of the base vertex
        bvvn = final_base_vertex_info.vertex_number;

        // Get the global vertex id of the base vertex
        global_vertex_id =
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn];

        // Check if the global vertex id has been already assigned
        if (global_vertex_id != -1)
        {
          // If that is the case then copy the global vertex id in the
          // current slot
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk]
                                                   [n_vertices - 1] =
                                                     global_vertex_id;
        } // if (global_vertex_id != -1)
        else
        {
          // Assign a global vertex id to the base vertex
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn] =
            global_vertex_number;

          // ... and copy the value to the final vertex
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk]
                                                   [n_vertices - 1] =
                                                     global_vertex_number;

          // ... increase the counter for the "global_vertex_number"
          global_vertex_number++;

          // Get the vertex
          Vector<double> vertex =
            tmp_polyline_pt->vertex_coordinate(n_vertices - 1);

          // ... and add it to the global vertex container
          global_vertices.push_back(vertex);
        }

        // Is the final vertex a base vertex
        if (final_base_vertex_info.is_base_vertex)
        {
          // Increase the independent vertices counter
          n_vertices_outer_polygons++;
        }

      } // for (p < n_polylines)

    } // for (i < n_outer_polygons)

    // The total number of added vertices in the internal polygons
    unsigned n_vertices_internal_polygons = 0;

    // Do the internal polygons
    for (unsigned i = 0; i < n_internal_polygons; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = internal_polygons_pt[i]->npolyline();

      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get a pointer to the current polyline
        TriangleMeshPolyLine* tmp_polyline_pt =
          internal_polygons_pt[i]->polyline_pt(p);

        // Get the boundary id of the current polyline
        const unsigned bound_id = tmp_polyline_pt->boundary_id();

        // The number of vertices in the current polyline
        const unsigned n_vertices = tmp_polyline_pt->nvertex();

        // Get the current chunk number of the polyline
        const unsigned bnd_chunk = tmp_polyline_pt->boundary_chunk();

        // Assign a global vertex id to the initial vertex
        // -----------------------------------------------

        // Get the base vertex information of the initial vertex
        base_vertex_info initial_base_vertex_info =
          base_vertices[bound_id][bnd_chunk][0];

#ifdef PARANOID
        if (!initial_base_vertex_info.has_base_vertex_assigned)
        {
          std::ostringstream error_message;
          std::string output_string =
            "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
          error_message
            << "The initial vertex of the current polyline has no base\n"
            << "vertex assigned\n"
            << "Internal polygon number: (" << i << ")\n\n"
            << "Polyline number: (" << p << ")\n"
            << "Boundary id: (" << bound_id << ")\n"
            << "Boundary chunk: (" << bnd_chunk << ")\n";
          throw OomphLibError(
            error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
        } // if (!initial_base_vertex_info.has_base_vertex_assigned)
#endif

        // The base vertex boundary id
        unsigned bvbi = initial_base_vertex_info.boundary_id;
        // The base vertex boundary chunkx
        unsigned bvbc = initial_base_vertex_info.boundary_chunk;
        // The vertex number of the base vertex
        unsigned bvvn = initial_base_vertex_info.vertex_number;

        // Get the global vertex id of the base vertex
        int global_vertex_id =
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn];

        // Check if the global vertex id has been already assigned
        if (global_vertex_id != -1)
        {
          // If that is the case then copy the global vertex id in the
          // current slot
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][0] =
            global_vertex_id;
        } // if (global_vertex_id != -1)
        else
        {
          // Assign a global vertex id to the base vertex
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn] =
            global_vertex_number;

          // ... and copy the value to the initial vertex
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][0] =
            global_vertex_number;

          // ... increase the counter for the "global_vertex_number"
          global_vertex_number++;

          // Get the vertex
          Vector<double> vertex = tmp_polyline_pt->vertex_coordinate(0);

          // ... and add it to the global vertex container
          global_vertices.push_back(vertex);
        }

        // Is the initial vertex a base vertex
        if (initial_base_vertex_info.is_base_vertex)
        {
          // Increase the independent vertices counter
          n_vertices_internal_polygons++;
        }

        // ------------------------------------------------------------
        // Now loop over the intermediate vertices and assign a unique
        // vertex id (all intermediate vertices are base vertices)
        for (unsigned v = 1; v < n_vertices - 1; v++)
        {
          // Get the global vertex id
          global_vertex_id =
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][v];

          // Check if the global vertex id has been already assigned.
          // We do nothing if it has been already assigned
          // (global_vertex_id!=-1).

          // If it has not been already assigned (global_vertex_id==-1)
          // then set a new global vertex number, and add the vertex to
          // the global vertices container
          if (global_vertex_id == -1)
          {
            // Set a value for the global vertex
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][v] =
              global_vertex_number;

            // ... increase the counter for the "global_vertex_number"
            global_vertex_number++;

            // Add the vertex to the global vertex container
            Vector<double> vertex = tmp_polyline_pt->vertex_coordinate(v);
            // Add the vertex to the global vertex container
            global_vertices.push_back(vertex);

          } // if (global_vertex_id == -1)

          // Increase the independent vertices counter
          n_vertices_internal_polygons++;

        } // for (n_vertices-1)

        // Assign a global vertex id to the final vertex
        // -----------------------------------------------

        // Get the base vertex information of the final vertex
        base_vertex_info final_base_vertex_info =
          base_vertices[bound_id][bnd_chunk][1];

#ifdef PARANOID
        if (!final_base_vertex_info.has_base_vertex_assigned)
        {
          std::ostringstream error_message;
          std::string output_string =
            "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
          error_message
            << "The final vertex of the current polyline has no base\n"
            << "vertex assigned\n"
            << "Internal polygon number: (" << i << ")\n\n"
            << "Polyline number: (" << p << ")\n"
            << "Boundary id: (" << bound_id << ")\n"
            << "Boundary chunk: (" << bnd_chunk << ")\n";
          throw OomphLibError(
            error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
        } // if (!final_base_vertex_info.has_base_vertex_assigned)
#endif

        // The base vertex boundary id
        bvbi = final_base_vertex_info.boundary_id;
        // The base vertex boundary chunkx
        bvbc = final_base_vertex_info.boundary_chunk;
        // The vertex number of the base vertex
        bvvn = final_base_vertex_info.vertex_number;

        // Get the global vertex id of the base vertex
        global_vertex_id =
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn];

        // Check if the global vertex id has been already assigned
        if (global_vertex_id != -1)
        {
          // If that is the case then copy the global vertex id in the
          // current slot
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk]
                                                   [n_vertices - 1] =
                                                     global_vertex_id;
        } // if (global_vertex_id != -1)
        else
        {
          // Assign a global vertex id to the base vertex
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn] =
            global_vertex_number;

          // ... and copy the value to the final vertex
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk]
                                                   [n_vertices - 1] =
                                                     global_vertex_number;

          // ... increase the counter for the "global_vertex_number"
          global_vertex_number++;

          // Get the vertex
          Vector<double> vertex =
            tmp_polyline_pt->vertex_coordinate(n_vertices - 1);

          // ... and add it to the global vertex container
          global_vertices.push_back(vertex);
        }

        // Is the final vertex a base vertex
        if (final_base_vertex_info.is_base_vertex)
        {
          // Increase the independent vertices counter
          n_vertices_internal_polygons++;
        }

      } // for (p < n_polylines)

    } // for (i < n_internal_polygons)

    // The total number of added vertices in the open curves
    unsigned n_vertices_open_curves = 0;

    // Do the open curves
    for (unsigned i = 0; i < n_open_curves; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = open_curves_pt[i]->ncurve_section();

      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get a pointer to the current polyline
        TriangleMeshPolyLine* tmp_polyline_pt =
          open_curves_pt[i]->polyline_pt(p);

        // Get the boundary id of the current polyline
        const unsigned bound_id = tmp_polyline_pt->boundary_id();

        // The number of vertices in the current polyline
        const unsigned n_vertices = tmp_polyline_pt->nvertex();

        // Get the current chunk number of the polyline
        const unsigned bnd_chunk = tmp_polyline_pt->boundary_chunk();

        // Assign a global vertex id to the initial vertex
        // -----------------------------------------------

        // Get the base vertex information of the initial vertex
        base_vertex_info initial_base_vertex_info =
          base_vertices[bound_id][bnd_chunk][0];

#ifdef PARANOID
        if (!initial_base_vertex_info.has_base_vertex_assigned)
        {
          std::ostringstream error_message;
          std::string output_string =
            "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
          error_message
            << "The initial vertex of the current polyline has no base\n"
            << "vertex assigned\n"
            << "Open curve number: (" << i << ")\n\n"
            << "Polyline number: (" << p << ")\n"
            << "Boundary id: (" << bound_id << ")\n"
            << "Boundary chunk: (" << bnd_chunk << ")\n";
          throw OomphLibError(
            error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
        } // if (!initial_base_vertex_info.has_base_vertex_assigned)
#endif

        // The base vertex boundary id
        unsigned bvbi = initial_base_vertex_info.boundary_id;
        // The base vertex boundary chunkx
        unsigned bvbc = initial_base_vertex_info.boundary_chunk;
        // The vertex number of the base vertex
        unsigned bvvn = initial_base_vertex_info.vertex_number;

        // Get the global vertex id of the base vertex
        int global_vertex_id =
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn];

        // Check if the global vertex id has been already assigned
        if (global_vertex_id != -1)
        {
          // If that is the case then copy the global vertex id in the
          // current slot
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][0] =
            global_vertex_id;
        } // if (global_vertex_id != -1)
        else
        {
          // Assign a global vertex id to the base vertex
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn] =
            global_vertex_number;

          // ... and copy the value to the initial vertex
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][0] =
            global_vertex_number;

          // ... increase the counter for the "global_vertex_number"
          global_vertex_number++;

          // Get the vertex
          Vector<double> vertex = tmp_polyline_pt->vertex_coordinate(0);

          // ... and add it to the global vertex container
          global_vertices.push_back(vertex);
        }

        // Is the initial vertex a base vertex
        if (initial_base_vertex_info.is_base_vertex)
        {
          // Increase the independent vertices counter
          n_vertices_open_curves++;
        }

        // ------------------------------------------------------------
        // Now loop over the intermediate vertices and assign a unique
        // vertex id (all intermediate vertices are base vertices)
        for (unsigned v = 1; v < n_vertices - 1; v++)
        {
          // Get the global vertex id
          global_vertex_id =
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][v];

          // Check if the global vertex id has been already assigned.
          // We do nothing if it has been already assigned
          // (global_vertex_id!=-1).

          // If it has not been already assigned (global_vertex_id==-1)
          // then set a new global vertex number, and add the vertex to
          // the global vertices container
          if (global_vertex_id == -1)
          {
            // Set a value for the global vertex
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][v] =
              global_vertex_number;

            // ... increase the counter for the "global_vertex_number"
            global_vertex_number++;

            // Add the vertex to the global vertex container
            Vector<double> vertex = tmp_polyline_pt->vertex_coordinate(v);
            // Add the vertex to the global vertex container
            global_vertices.push_back(vertex);

          } // if (global_vertex_id == -1)

          // Increase the independent vertices counter
          n_vertices_open_curves++;

        } // for (n_vertices-1)

        // Assign a global vertex id to the final vertex
        // -----------------------------------------------

        // Get the base vertex information of the final vertex
        base_vertex_info final_base_vertex_info =
          base_vertices[bound_id][bnd_chunk][1];

#ifdef PARANOID
        if (!final_base_vertex_info.has_base_vertex_assigned)
        {
          std::ostringstream error_message;
          std::string output_string =
            "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
          error_message
            << "The final vertex of the current polyline has no base\n"
            << "vertex assigned\n"
            << "Open curve number: (" << i << ")\n\n"
            << "Polyline number: (" << p << ")\n"
            << "Boundary id: (" << bound_id << ")\n"
            << "Boundary chunk: (" << bnd_chunk << ")\n";
          throw OomphLibError(
            error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
        } // if (!final_base_vertex_info.has_base_vertex_assigned)
#endif

        // The base vertex boundary id
        bvbi = final_base_vertex_info.boundary_id;
        // The base vertex boundary chunkx
        bvbc = final_base_vertex_info.boundary_chunk;
        // The vertex number of the base vertex
        bvvn = final_base_vertex_info.vertex_number;

        // Get the global vertex id of the base vertex
        global_vertex_id =
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn];

        // Check if the global vertex id has been already assigned
        if (global_vertex_id != -1)
        {
          // If that is the case then copy the global vertex id in the
          // current slot
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk]
                                                   [n_vertices - 1] =
                                                     global_vertex_id;
        } // if (global_vertex_id != -1)
        else
        {
          // Assign a global vertex id to the base vertex
          boundary_chunk_vertex_to_global_vertex_id[bvbi][bvbc][bvvn] =
            global_vertex_number;

          // ... and copy the value to the final vertex
          boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk]
                                                   [n_vertices - 1] =
                                                     global_vertex_number;

          // ... increase the counter for the "global_vertex_number"
          global_vertex_number++;

          // Get the vertex
          Vector<double> vertex =
            tmp_polyline_pt->vertex_coordinate(n_vertices - 1);

          // ... and add it to the global vertex container
          global_vertices.push_back(vertex);
        }

        // Is the final vertex a base vertex
        if (final_base_vertex_info.is_base_vertex)
        {
          // Increase the independent vertices counter
          n_vertices_open_curves++;
        }

      } // for (p < n_polylines)

    } // for (i < n_open_curves)

    // We have already assigned a unique id for each vertex in the
    // boundaries, and we have collected all the vertices in a container

    // Store the global number of vertices. Add the number of vertices
    // of all the polygons (outer and internal), and the open curves to
    // get the total number of vertices. NB This is the number of NON
    // REPEATED vertices, the sum of all the vertices of each individual
    // polyline can be computed from the boundary_chunk_nvertices
    // container
    const unsigned n_global_vertices = n_vertices_outer_polygons +
                                       n_vertices_internal_polygons +
                                       n_vertices_open_curves;

#ifdef PARANOID
    // Check that the total number of vertices be equal to the
    // pre-computed number of vertices
    const unsigned added_global_vertices = global_vertices.size();
    if (added_global_vertices != n_global_vertices)
    {
      std::ostringstream error_stream;
      std::string output_string =
        "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
      error_stream
        << "The number of added vertices to the global vertices container\n"
        << "is different from the pre-computed number of vertices\n"
        << "Added number of vertices: (" << added_global_vertices << ")\n"
        << "Pre-computed number of global vertices: (" << n_global_vertices
        << ")\n";
      throw OomphLibError(
        error_stream.str(), output_string, OOMPH_EXCEPTION_LOCATION);
    } // if (added_global_vertices != n_global_vertices)

    // Check that the counter for the global number of vertices is the same as
    // the pre-computed number of vertices
    if (global_vertex_number != n_global_vertices)
    {
      std::ostringstream error_stream;
      std::string output_string =
        "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
      error_stream
        << "The counter for the global number of vertices is different from\n"
        << "the pre-computed number of vertices\n"
        << "Global vertices counter: (" << global_vertex_number << ")\n"
        << "Pre-computed number of global vertices: (" << n_global_vertices
        << ")\n";
      throw OomphLibError(
        error_stream.str(), output_string, OOMPH_EXCEPTION_LOCATION);
    } // if (global_vertex_number != n_global_vertices)
#endif


    // ------------------------------------------------------------------

    // ------------------------------------------------------------------
    // 4th- stage
    // Create the segments using the unique id of the vertices and
    // assign a segment id to each segment (the one from the boundary)

    // Create the segment storage, each segment is composed of two
    // vertices (the vertices id is obtained from the global vertex ids
    // container)
    Vector<Vector<int>> global_segments;

    // Assign a segment marked to each segment (the boundary id to which
    // the segment belongs)
    Vector<int> global_segment_markers;

    // Loop over the boundaries again, but know get the global vertex id
    // from the container and create the segments

    // Start with the outer polygons
    for (unsigned i = 0; i < n_outer_polygons; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = outer_polygons_pt[i]->npolyline();

      // Loop over the polylines
      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get the boundary id of the current polyline
        const unsigned bound_id =
          outer_polygons_pt[i]->polyline_pt(p)->boundary_id();

        // The number of vertices in the current polyline
        const unsigned n_vertices =
          outer_polygons_pt[i]->polyline_pt(p)->nvertex();

        // Get the current chunk number of the polyline
        const unsigned bnd_chunk =
          outer_polygons_pt[i]->polyline_pt(p)->boundary_chunk();

        // Loop over the vertices in the polyline
        for (unsigned v = 1; v < n_vertices; v++)
        {
          // Get the global vertex id
          const int global_vertex_id_v_1 =
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk]
                                                     [v - 1];

#ifdef PARANOID
          if (global_vertex_id_v_1 == -1)
          {
            // Get the current vertex
            Vector<double> current_vertex =
              outer_polygons_pt[i]->polyline_pt(p)->vertex_coordinate(v - 1);
            std::ostringstream error_message;
            std::string output_string =
              "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
            error_message
              << "The global vertex number for the current vertex has not\n"
              << "been assigned\n."
              << "Outer polygon number: (" << i << ")\n\n"
              << "Polyline number: (" << p << ")\n"
              << "Boundary id: (" << bound_id << ")\n"
              << "Boundary chunk: (" << bnd_chunk << ")\n"
              << "Vertex[" << v - 1 << "]: (" << current_vertex[0] << ", "
              << current_vertex[1] << ")\n";
            throw OomphLibError(
              error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
          } // if (global_vertex_id_v_1 == -1)
#endif

          // Get the global vertex id
          const int global_vertex_id_v =
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][v];

#ifdef PARANOID
          if (global_vertex_id_v == -1)
          {
            // Get the current vertex
            Vector<double> current_vertex =
              outer_polygons_pt[i]->polyline_pt(p)->vertex_coordinate(v);
            std::ostringstream error_message;
            std::string output_string =
              "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
            error_message
              << "The global vertex number for the current vertex has not\n"
              << "been assigned\n."
              << "Outer polygon number: (" << i << ")\n\n"
              << "Polyline number: (" << p << ")\n"
              << "Boundary id: (" << bound_id << ")\n"
              << "Boundary chunk: (" << bnd_chunk << ")\n"
              << "Vertex[" << v << "]: (" << current_vertex[0] << ", "
              << current_vertex[1] << ")\n";
            throw OomphLibError(
              error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
          } // if (global_vertex_id_v == -1)
#endif

          // Add the vertices to the segments container
          Vector<int> current_segment(2);
          current_segment[0] = global_vertex_id_v_1;
          current_segment[1] = global_vertex_id_v;
          global_segments.push_back(current_segment);

          // ... and associate the segments marker
          global_segment_markers.push_back(bound_id);

        } // for (v < n_vertices)

      } // for (p < n_polylines)

    } // for (i < n_outer_polygons)

    // Now work the internal polygons
    for (unsigned i = 0; i < n_internal_polygons; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = internal_polygons_pt[i]->npolyline();

      // Loop over the polylines
      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get the boundary id of the current polyline
        const unsigned bound_id =
          internal_polygons_pt[i]->polyline_pt(p)->boundary_id();

        // The number of vertices in the current polyline
        const unsigned n_vertices =
          internal_polygons_pt[i]->polyline_pt(p)->nvertex();

        // Get the current chunk number of the polyline
        const unsigned bnd_chunk =
          internal_polygons_pt[i]->polyline_pt(p)->boundary_chunk();

        // Loop over the vertices in the polyline
        for (unsigned v = 1; v < n_vertices; v++)
        {
          // Get the global vertex id
          const int global_vertex_id_v_1 =
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk]
                                                     [v - 1];

#ifdef PARANOID
          if (global_vertex_id_v_1 == -1)
          {
            // Get the current vertex
            Vector<double> current_vertex =
              internal_polygons_pt[i]->polyline_pt(p)->vertex_coordinate(v - 1);
            std::ostringstream error_message;
            std::string output_string =
              "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
            error_message
              << "The global vertex number for the current vertex has not\n"
              << "been assigned\n."
              << "Internal polygon number: (" << i << ")\n\n"
              << "Polyline number: (" << p << ")\n"
              << "Boundary id: (" << bound_id << ")\n"
              << "Boundary chunk: (" << bnd_chunk << ")\n"
              << "Vertex[" << v - 1 << "]: (" << current_vertex[0] << ", "
              << current_vertex[1] << ")\n";
            throw OomphLibError(
              error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
          } // if (global_vertex_id_v_1 == -1)
#endif

          // Get the global vertex id
          const int global_vertex_id_v =
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][v];

#ifdef PARANOID
          if (global_vertex_id_v == -1)
          {
            // Get the current vertex
            Vector<double> current_vertex =
              internal_polygons_pt[i]->polyline_pt(p)->vertex_coordinate(v);
            std::ostringstream error_message;
            std::string output_string =
              "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
            error_message
              << "The global vertex number for the current vertex has not\n"
              << "been assigned\n."
              << "Internal polygon number: (" << i << ")\n\n"
              << "Polyline number: (" << p << ")\n"
              << "Boundary id: (" << bound_id << ")\n"
              << "Boundary chunk: (" << bnd_chunk << ")\n"
              << "Vertex[" << v << "]: (" << current_vertex[0] << ", "
              << current_vertex[1] << ")\n";
            throw OomphLibError(
              error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
          } // if (global_vertex_id_v == -1)
#endif

          // Add the vertices to the segments container
          Vector<int> current_segment(2);
          current_segment[0] = global_vertex_id_v_1;
          current_segment[1] = global_vertex_id_v;
          global_segments.push_back(current_segment);

          // ... and associate the segments marker
          global_segment_markers.push_back(bound_id);

        } // for (v < n_vertices)

      } // for (p < n_polylines)

    } // for (i < n_internal_polygons)

    // Now do the open curves
    for (unsigned i = 0; i < n_open_curves; i++)
    {
      // The number of polylines in the current polygon
      const unsigned n_polylines = open_curves_pt[i]->ncurve_section();

      // Loop over the polylines
      for (unsigned p = 0; p < n_polylines; p++)
      {
        // Get the boundary id of the current polyline
        const unsigned bound_id =
          open_curves_pt[i]->polyline_pt(p)->boundary_id();

        // The number of vertices in the current polyline
        const unsigned n_vertices =
          open_curves_pt[i]->polyline_pt(p)->nvertex();

        // Get the current chunk number of the polyline
        const unsigned bnd_chunk =
          open_curves_pt[i]->polyline_pt(p)->boundary_chunk();

        // Loop over the vertices in the polyline
        for (unsigned v = 1; v < n_vertices; v++)
        {
          // Get the global vertex id
          const int global_vertex_id_v_1 =
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk]
                                                     [v - 1];

#ifdef PARANOID
          if (global_vertex_id_v_1 == -1)
          {
            // Get the current vertex
            Vector<double> current_vertex =
              open_curves_pt[i]->polyline_pt(p)->vertex_coordinate(v - 1);
            std::ostringstream error_message;
            std::string output_string =
              "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
            error_message
              << "The global vertex number for the current vertex has not\n"
              << "been assigned\n."
              << "Open curve number: (" << i << ")\n\n"
              << "Polyline number: (" << p << ")\n"
              << "Boundary id: (" << bound_id << ")\n"
              << "Boundary chunk: (" << bnd_chunk << ")\n"
              << "Vertex[" << v - 1 << "]: (" << current_vertex[0] << ", "
              << current_vertex[1] << ")\n";
            throw OomphLibError(
              error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
          } // if (global_vertex_id_v_1 == -1)
#endif

          // Get the global vertex id
          const int global_vertex_id_v =
            boundary_chunk_vertex_to_global_vertex_id[bound_id][bnd_chunk][v];

#ifdef PARANOID
          if (global_vertex_id_v == -1)
          {
            // Get the current vertex
            Vector<double> current_vertex =
              open_curves_pt[i]->polyline_pt(p)->vertex_coordinate(v);
            std::ostringstream error_message;
            std::string output_string =
              "UnstructuredTwoDMeshGeometryBase::build_triangulateio()";
            error_message
              << "The global vertex number for the current vertex has not\n"
              << "been assigned\n."
              << "Open curve number: (" << i << ")\n\n"
              << "Polyline number: (" << p << ")\n"
              << "Boundary id: (" << bound_id << ")\n"
              << "Boundary chunk: (" << bnd_chunk << ")\n"
              << "Vertex[" << v << "]: (" << current_vertex[0] << ", "
              << current_vertex[1] << ")\n";
            throw OomphLibError(
              error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
          } // if (global_vertex_id_v == -1)
#endif

          // Add the vertices to the segments container
          Vector<int> current_segment(2);
          current_segment[0] = global_vertex_id_v_1;
          current_segment[1] = global_vertex_id_v;
          global_segments.push_back(current_segment);

          // ... and associate the segments marker
          global_segment_markers.push_back(bound_id);

        } // for (v < n_vertices)

      } // for (p < n_polylines)

    } // for (i < n_open_curves)

    // ------------------------------------------------------------------
    // 5th- stage
    // Fill the triangulateio containers with the collected information

    // Initialize the triangulateIO structure
    TriangleHelper::initialise_triangulateio(triangulate_io);

    // Set the number of vertices
    triangulate_io.numberofpoints = n_global_vertices;

    // Get the number of segments
    const unsigned n_global_segments = global_segments.size();
    // Set the number of segments
    triangulate_io.numberofsegments = n_global_segments;

    // Allocate memory in the triangulateIO structure to store the values
    triangulate_io.pointlist =
      (double*)malloc(triangulate_io.numberofpoints * 2 * sizeof(double));
    triangulate_io.segmentlist =
      (int*)malloc(triangulate_io.numberofsegments * 2 * sizeof(int));
    triangulate_io.segmentmarkerlist =
      (int*)malloc(triangulate_io.numberofsegments * sizeof(int));

    // Fill triangulateIO data structure
    // ---------------------------------------

    // Fill the vertices first
    // An external counter
    unsigned ii = 0;
    for (unsigned i = 0; i < n_global_vertices; i++, ii += 2)
    {
      triangulate_io.pointlist[ii] = global_vertices[i][0];
      triangulate_io.pointlist[ii + 1] = global_vertices[i][1];
    } // for (i < n_global_vertices)

    // Then fill the segments, and segments markers
    // Reset the external counter
    ii = 0;
    for (unsigned i = 0; i < n_global_segments; i++, ii += 2)
    {
      // The segment marker should start in 1 (our enumeration started
      // in 0, therefore we add one to every entry)
      triangulate_io.segmentlist[ii] = global_segments[i][0] + 1;
      triangulate_io.segmentlist[ii + 1] = global_segments[i][1] + 1;
      // Set the segment marker as the boundary id + 1
      triangulate_io.segmentmarkerlist[i] = global_segment_markers[i] + 1;
    }

    // Store any regions information
    // ---------------------------------------

    // Get the number of regions
    unsigned n_regions = regions_coordinates.size();

    // Check if there are regions to include
    if (n_regions > 0)
    {
      triangulate_io.numberofregions = n_regions;
      triangulate_io.regionlist =
        (double*)malloc(triangulate_io.numberofregions * 4 * sizeof(double));

      // Loop over the regions map
      unsigned p = 1;
      for (std::map<unsigned, Vector<double>>::iterator it_regions =
             regions_coordinates.begin();
           it_regions != regions_coordinates.end();
           it_regions++)
      {
        // Get the region id
        unsigned region_id = (*it_regions).first;
        // Set the x-coordinate
        triangulate_io.regionlist[4 * p - 4] = ((*it_regions).second)[0];
        // Set the y-coordinate
        triangulate_io.regionlist[4 * p - 3] = ((*it_regions).second)[1];
        // Set the region id
        triangulate_io.regionlist[4 * p - 2] = static_cast<double>(region_id);
        // This is used to define a target area for the region, zero
        // means no target area defined
        triangulate_io.regionlist[4 * p - 1] = regions_areas[region_id];
        // Increase the auxiliary counter
        p++;
      } // Loop over the regions map

    } // if (n_regions > 0)

    // Holes information
    // ---------------------------------------

    // Get the number of any extra defined holes
    const unsigned n_extra_holes = extra_holes_coordinates.size();
    // The internal polygon marked as a hole
    Vector<unsigned> internal_polygon_marked_as_hole;
    // Count how many internal polygons are holes
    for (unsigned h = 0; h < n_internal_polygons; h++)
    {
      if (!internal_polygons_pt[h]->internal_point().empty())
      {
        internal_polygon_marked_as_hole.push_back(h);
      } // Is internal polygon a hole?
    } // for (h < n_internal_polygons)

    // Get the number of internal polygons that should be treated as
    // holes
    const unsigned n_internal_polygons_are_holes =
      internal_polygon_marked_as_hole.size();

    // Set the number of holes in the triangulateIO structure
    triangulate_io.numberofholes =
      n_extra_holes + n_internal_polygons_are_holes;

    // Allocate memory for the holes coordinates
    triangulate_io.holelist =
      (double*)malloc(triangulate_io.numberofholes * 2 * sizeof(double));

    // Store the holes coordinates
    unsigned count_hole = 0;
    // Add first the internal polygons marked as holes
    for (; count_hole < n_internal_polygons_are_holes * 2; count_hole += 2)
    {
      const unsigned index_intenal_polygon_is_hole =
        internal_polygon_marked_as_hole[count_hole / 2];
      triangulate_io.holelist[count_hole] =
        internal_polygons_pt[index_intenal_polygon_is_hole]
          ->internal_point()[0];
      triangulate_io.holelist[count_hole + 1] =
        internal_polygons_pt[index_intenal_polygon_is_hole]
          ->internal_point()[1];
    } // for (count_hole < n_internal_polygons_are_holes*2)

    // Add the extra holes coordinates
    for (unsigned i_eh = 0;
         count_hole < 2 * (n_extra_holes + n_internal_polygons_are_holes);
         count_hole += 2, i_eh++)
    {
      triangulate_io.holelist[count_hole] = extra_holes_coordinates[i_eh][0];
      triangulate_io.holelist[count_hole + 1] =
        extra_holes_coordinates[i_eh][1];
    } // for (count_hole < 2*(n_extra_holes + n_internal_polygons_are_holes))

    // ------------------------------------------------------------------
  }

  //========================================================================
  /// Helps to add information to the connection matrix of the
  /// given polyline
  //========================================================================
  void UnstructuredTwoDMeshGeometryBase::add_connection_matrix_info_helper(
    TriangleMeshPolyLine* polyline_pt,
    std::map<unsigned, std::map<unsigned, Vector<vertex_connection_info>>>&
      connection_matrix,
    TriangleMeshPolyLine* next_polyline_pt)
  {
    // Get the boundary id of the current polyline
    const unsigned bound_id = polyline_pt->boundary_id();

    // Get the chunk number associated with this boundary
    const unsigned bound_chunk = polyline_pt->boundary_chunk();

    // Check if the ends of the polyline are connected
    const bool connected_init_end = polyline_pt->is_initial_vertex_connected();

    // ... check for both connections
    const bool connected_final_end = polyline_pt->is_final_vertex_connected();

    // Prepare the connection matrix (resize the vector to store the
    // initial and final vertices info.)
    connection_matrix[bound_id][bound_chunk].resize(2);

    // The default information for the matrix
    vertex_connection_info default_vertex_connection_info;
    default_vertex_connection_info.is_connected = false;
    default_vertex_connection_info.boundary_id_to_connect = 0;
    default_vertex_connection_info.boundary_chunk_to_connect = 0;
    default_vertex_connection_info.vertex_number_to_connect = 0;

    // Set the default information for the initial vertex
    connection_matrix[bound_id][bound_chunk][0] =
      default_vertex_connection_info;
    // Set the default information for the final vertex
    connection_matrix[bound_id][bound_chunk][1] =
      default_vertex_connection_info;

    // Set the info. for the connection matrix
    if (connected_init_end)
    {
      // The boundary id to connect
      const unsigned bc = polyline_pt->initial_vertex_connected_bnd_id();

      // The vertex number to which is connected
      const unsigned vc = polyline_pt->initial_vertex_connected_n_vertex();

      // The boundary chunk to which is connected
      const unsigned cc = polyline_pt->initial_vertex_connected_n_chunk();

      // Set the connection info
      vertex_connection_info initial_vertex_info;
      // Set as connected
      initial_vertex_info.is_connected = true;
      // The boundary id to connect
      initial_vertex_info.boundary_id_to_connect = bc;
      // The chunk number to connect
      initial_vertex_info.boundary_chunk_to_connect = cc;
      // The vertex number to connect
      initial_vertex_info.vertex_number_to_connect = vc;

      // Add the connection information to the connection matrix
      connection_matrix[bound_id][bound_chunk][0] = initial_vertex_info;

    } // if (connected_init_end)

    if (connected_final_end)
    {
      // The boundary id to connect
      const unsigned bc = polyline_pt->final_vertex_connected_bnd_id();

      // The vertex number to which is connected
      const unsigned vc = polyline_pt->final_vertex_connected_n_vertex();

      // The boundary chunk to which is connected
      const unsigned cc = polyline_pt->final_vertex_connected_n_chunk();

      // Set the connection info
      vertex_connection_info final_vertex_info;
      // Set as connected
      final_vertex_info.is_connected = true;
      // The boundary id to connect
      final_vertex_info.boundary_id_to_connect = bc;
      // The chunk number to connect
      final_vertex_info.boundary_chunk_to_connect = cc;
      // The vertex number to connect
      final_vertex_info.vertex_number_to_connect = vc;

      // Add the connection information to the connection matrix
      connection_matrix[bound_id][bound_chunk][1] = final_vertex_info;

    } // if (connected_final_end)
    else
    {
      // If the vertex is not connected but there is a next polyline in
      // the polygon (this only applies for polygons, for open curves
      // the next polyline is set to null)

      if (next_polyline_pt != 0)
      {
        // Get the info of the next polyline
        const unsigned next_bound_id = next_polyline_pt->boundary_id();

        // Get the chunk number associated with this boundary
        const unsigned next_bound_chunk = next_polyline_pt->boundary_chunk();

        // Set the connection info
        vertex_connection_info final_vertex_info;
        // Set as connected
        final_vertex_info.is_connected = true;
        // The boundary id to connect
        final_vertex_info.boundary_id_to_connect = next_bound_id;
        // The chunk number to connect
        final_vertex_info.boundary_chunk_to_connect = next_bound_chunk;
        // The vertex number to connect
        final_vertex_info.vertex_number_to_connect = 0;

        // Add the connection information to the connection matrix
        connection_matrix[bound_id][bound_chunk][1] = final_vertex_info;

      } // if (next_polyline_pt != 0)

    } // else if (connected_final_end)
  }

  //========================================================================
  // Initialize the base vertex structure, set every vertex to
  // no visited and not being a base vertex
  //========================================================================
  void UnstructuredTwoDMeshGeometryBase::initialise_base_vertex(
    TriangleMeshPolyLine* polyline_pt,
    std::map<unsigned, std::map<unsigned, Vector<base_vertex_info>>>&
      base_vertices)
  {
    // Get the boundary id of the current polyline
    const unsigned bound_id = polyline_pt->boundary_id();

    // Get the chunk number associated with this boundary
    const unsigned bound_chunk = polyline_pt->boundary_chunk();

    // Create the default data structure
    base_vertex_info default_base_vertex_info;
    // Set as not done
    default_base_vertex_info.has_base_vertex_assigned = false;
    // Set as not base vertex
    default_base_vertex_info.is_base_vertex = false;
    // Set the default base boundary id
    default_base_vertex_info.boundary_id = 0;
    // Set the default base boundary chunk
    default_base_vertex_info.boundary_chunk = 0;
    // Set the default base vertex number
    default_base_vertex_info.vertex_number = 0;

    // Initialize the base vertex info. for the current polyline
    // Allocate memory for the initial and final vertex
    base_vertices[bound_id][bound_chunk].resize(2);
    // Set the initial vertex info.
    base_vertices[bound_id][bound_chunk][0] = default_base_vertex_info;
    // Set the final vertex info.
    base_vertices[bound_id][bound_chunk][1] = default_base_vertex_info;
  }

  //========================================================================
  // Helps to identify the base vertex of the given polyline
  //========================================================================
  void UnstructuredTwoDMeshGeometryBase::add_base_vertex_info_helper(
    TriangleMeshPolyLine* polyline_pt,
    std::map<unsigned, std::map<unsigned, Vector<base_vertex_info>>>&
      base_vertices,
    std::map<unsigned, std::map<unsigned, Vector<vertex_connection_info>>>&
      connection_matrix,
    std::map<unsigned, std::map<unsigned, unsigned>>& boundary_chunk_nvertices)
  {
    // Get the general data of the polyline

    // Get the boundary id of the current polyline
    const unsigned bound_id = polyline_pt->boundary_id();

    // The number of vertices in the current polyline
    const unsigned n_vertices = polyline_pt->nvertex();

    // Get the chunk number associated with this boundary
    const unsigned bound_chunk = polyline_pt->boundary_chunk();

    // Keep track of the done vertices
    std::map<Vector<unsigned>, bool> done;

    // Loop over the vertices in the polyline
    for (unsigned v = 0; v < n_vertices; v++)
    {
      // For each vertex find its base vertex

      // Flag to state if repeat the search with the new vertex
      // reference
      bool repeat = false;

      // Store the info. of the vertex currently serching for base
      // vertex
      unsigned ibnd_id = bound_id;
      unsigned ibnd_chunk = bound_chunk;
      unsigned ivertex_number = v;

      // Store the info. of the base vertex of the current vertex
      unsigned base_bnd_id = 0;
      unsigned base_bnd_chunk = 0;
      unsigned base_vertex_number = 0;

      // Store all the found references to the vertex
      Vector<unsigned> iter_bnd_id;
      Vector<unsigned> iter_bnd_chunk;
      Vector<unsigned> iter_vertex_number;

      do
      {
        // Reset the flag
        repeat = false;

        // Get the number of vertices of the polyline where the current
        // reference index lives on
        const unsigned i_n_vertices =
          boundary_chunk_nvertices[ibnd_id][ibnd_chunk];
        // Is the vertex an initial or final vertex on the (ibnd_id,
        // ibnd_chunk)
        bool is_initial = false;
        if (ivertex_number == 0)
        {
          is_initial = true;
        }
        bool is_final = false;
        if (ivertex_number == i_n_vertices - 1)
        {
          is_final = true;
        }

        // Is initial or final?
        if (is_initial || is_final)
        {
          // Stores the base vertex info. of the current vertex
          base_vertex_info current_vertex_base_info;

          if (is_initial)
          {
            // Get the base vertex info. of the current vertex
            current_vertex_base_info = base_vertices[ibnd_id][ibnd_chunk][0];
          }
          else if (is_final)
          {
            // Get the base vertex info. of the current vertex
            current_vertex_base_info = base_vertices[ibnd_id][ibnd_chunk][1];
          }

          // Has the vertex assigned a base vertex?
          if (!current_vertex_base_info.has_base_vertex_assigned)
          {
            // Is the current vertex a base vertex?
            if (!current_vertex_base_info.is_base_vertex)
            {
              // Get the connection info. of the vertex
              vertex_connection_info vertex_connect_info =
                connection_matrix[ibnd_id][ibnd_chunk][0];

              if (is_initial)
              {
                vertex_connect_info = connection_matrix[ibnd_id][ibnd_chunk][0];
              }
              else if (is_final)
              {
                vertex_connect_info = connection_matrix[ibnd_id][ibnd_chunk][1];
              }

              // Get the complete connection information of the vertex

              // The new connection info.
              const bool current_is_connected =
                vertex_connect_info.is_connected;

              // Is the current vertex connected to other vertex
              if (current_is_connected)
              {
                // The new boundary id to connect
                const unsigned iibnd_id =
                  vertex_connect_info.boundary_id_to_connect;

                // The new boundary chunk to connect
                const unsigned iibnd_chunk =
                  vertex_connect_info.boundary_chunk_to_connect;

                // The new vertex number to connect
                const unsigned iivertex_number =
                  vertex_connect_info.vertex_number_to_connect;

                // Get the number of vertices in the new boundary to connect
                const unsigned ii_n_vertices =
                  boundary_chunk_nvertices[iibnd_id][iibnd_chunk];

                // Is the new vertex an initial or final vertex on the
                // (iibnd_id, iibnd_chunk)
                bool new_is_initial = false;
                if (iivertex_number == 0)
                {
                  new_is_initial = true;
                }
                bool new_is_final = false;
                if (iivertex_number == ii_n_vertices - 1)
                {
                  new_is_final = true;
                }

                // Is new initial or final?
                if (new_is_initial || new_is_final)
                {
                  // Stores the base vertex info. of the new vertex
                  base_vertex_info new_vertex_base_info;

                  if (new_is_initial)
                  {
                    // Get the base vertex info. of the current vertex
                    new_vertex_base_info =
                      base_vertices[iibnd_id][iibnd_chunk][0];
                  }
                  else if (new_is_final)
                  {
                    // Get the base vertex info. of the current vertex
                    new_vertex_base_info =
                      base_vertices[iibnd_id][iibnd_chunk][1];
                  }

                  // Verify if the new vertex has been done
                  Vector<unsigned> check_done(3);
                  check_done[0] = iibnd_id;
                  check_done[1] = iibnd_chunk;
                  check_done[2] = iivertex_number;
                  // Has the new vertex been already visited?
                  if (!done[check_done])
                  {
                    // Add the i-reference to the iteration vectors
                    // since it needs to be assigned the base node when
                    // it be found
                    iter_bnd_id.push_back(ibnd_id);
                    iter_bnd_chunk.push_back(ibnd_chunk);
                    iter_vertex_number.push_back(ivertex_number);

                    Vector<unsigned> current_done(3);
                    current_done[0] = ibnd_id;
                    current_done[1] = ibnd_chunk;
                    current_done[2] = ivertex_number;

                    // Mark the i-th reference of the vertex as done
                    done[current_done] = true;

                    // Change the vertex reference and repeat
                    ibnd_id = iibnd_id;
                    ibnd_chunk = iibnd_chunk;
                    ivertex_number = iivertex_number;

                    // Set the flag to repeat the search with the new
                    // reference of the vertex
                    repeat = true;

                  } // if (!done[check_done])
                  else
                  {
                    // We have found a base vertex, we only need to
                    // decide if it is the current vertex or the new
                    // vertex

                    // Add the i-reference to the iteration vectors to
                    // be assigned the found base vertex
                    iter_bnd_id.push_back(ibnd_id);
                    iter_bnd_chunk.push_back(ibnd_chunk);
                    iter_vertex_number.push_back(ivertex_number);

                    // Has the new vertex a base vertes assigned
                    if (!new_vertex_base_info.has_base_vertex_assigned)
                    {
                      // The current vertex is a base vertex (we are
                      // looping over the current and new vertex)

                      Vector<unsigned> current_done(3);
                      current_done[0] = ibnd_id;
                      current_done[1] = ibnd_chunk;
                      current_done[2] = ivertex_number;

                      // Mark the i-th reference of the vertex as done
                      done[current_done] = true;

                      // Mark the i-th reference of the vertex as base
                      // vertex
                      if (is_initial)
                      {
                        // Mark the vertex as base vertex
                        base_vertices[ibnd_id][ibnd_chunk][0].is_base_vertex =
                          true;
                      }
                      else if (is_final)
                      {
                        // Mark the vertex as base vertex
                        base_vertices[ibnd_id][ibnd_chunk][1].is_base_vertex =
                          true;
                      }

                      // Set as base vertex the current vertex info.
                      // The base boundary id
                      base_bnd_id = ibnd_id;
                      // The base boundary chunk number
                      base_bnd_chunk = ibnd_chunk;
                      // The base vertex number
                      base_vertex_number = ivertex_number;

                    } // if (!new_vertex_base_info.is_base_vertex)
                    else
                    {
                      // The new vertex is a base vertex (we are looping
                      // over the current and new vertex)

                      Vector<unsigned> current_done(3);
                      current_done[0] = ibnd_id;
                      current_done[1] = ibnd_chunk;
                      current_done[2] = ivertex_number;

                      // Mark the i-th reference of the vertex as done
                      done[current_done] = true;

                      // Set as base vertex the new vertex info.
                      // The base boundary id
                      base_bnd_id = iibnd_id;
                      // The base boundary chunk number
                      base_bnd_chunk = iibnd_chunk;
                      // The base vertex number
                      base_vertex_number = iivertex_number;

                    } // else if (!new_vertex_base_info.is_base_vertex)

                  } // else if (!new_vertex_base_info.done)

                } // if (new_is_initial || new_is_final)
                else
                {
                  // Add the i-reference to the iteration vectors since
                  // it needs to be assigned the just found base vertex
                  iter_bnd_id.push_back(ibnd_id);
                  iter_bnd_chunk.push_back(ibnd_chunk);
                  iter_vertex_number.push_back(ivertex_number);

                  Vector<unsigned> current_done(3);
                  current_done[0] = ibnd_id;
                  current_done[1] = ibnd_chunk;
                  current_done[2] = ivertex_number;

                  // Mark the i-th reference of the vertex as done
                  done[current_done] = true;

                  // Set as base node the new node info.
                  // The base boundary id
                  base_bnd_id = iibnd_id;
                  // The base boundary chunk number
                  base_bnd_chunk = iibnd_chunk;
                  // The base vertex number
                  base_vertex_number = iivertex_number;
                } // else if (new_is_initial || new_is_final)

              } // if (current_is_connected)
              else
              {
                // The current vertex is a base vertex

                // Add the i-reference to the iteration vectors to be
                // assigned the found base vertex
                iter_bnd_id.push_back(ibnd_id);
                iter_bnd_chunk.push_back(ibnd_chunk);
                iter_vertex_number.push_back(ivertex_number);

                Vector<unsigned> current_done(3);
                current_done[0] = ibnd_id;
                current_done[1] = ibnd_chunk;
                current_done[2] = ivertex_number;

                // Mark the i-th reference of the vertex as done
                done[current_done] = true;

                // Mark the i-th reference of the vertex as base vertex
                if (is_initial)
                {
                  // Mark the vertex as base vertex
                  base_vertices[ibnd_id][ibnd_chunk][0].is_base_vertex = true;
                }
                else if (is_final)
                {
                  // Mark the vertex as base vertex
                  base_vertices[ibnd_id][ibnd_chunk][1].is_base_vertex = true;
                }

                // Set as base vertex the current vertex info.
                // The base boundary id
                base_bnd_id = ibnd_id;
                // The base boundary chunk number
                base_bnd_chunk = ibnd_chunk;
                // The base vertex number
                base_vertex_number = ivertex_number;

              } // else if (current_is_connected)

            } // if (!current_vertex_base_info.is_base_vertex)
            else
            {
              // Copy the base vertex info. and assign it to all the
              // found vertex references

              // The base boundary id
              base_bnd_id = ibnd_id;
              // The base boundary chunk number
              base_bnd_chunk = ibnd_chunk;
              // The base vertex number
              base_vertex_number = ivertex_number;

            } // else if (!current_vertex_base_info.is_base_vertex)

          } // if (!current_vertex_base_info.has_base_vertex_assigned)
          else
          {
            // Copy the base vertex info. and assign it to all the found
            // vertex references

            // The base boundary id
            base_bnd_id = current_vertex_base_info.boundary_id;
            // The base boundary chunk number
            base_bnd_chunk = current_vertex_base_info.boundary_chunk;
            // The base vertex number
            base_vertex_number = current_vertex_base_info.vertex_number;
          } // else if (!current_vertex_base_info.has_base_vertex_assigned)

        } // if (is_initial || is_final)

      } while (repeat);

      // Assign the base vertex to all the found references of the
      // vertex

      // Get the number of found references
      const unsigned n_vertex_references = iter_bnd_id.size();
      // Loop over the found references and assign the base vertex
      for (unsigned r = 0; r < n_vertex_references; r++)
      {
        // Get the r-th reference of the vertex
        const unsigned rbnd_id = iter_bnd_id[r];
        const unsigned rbnd_chunk = iter_bnd_chunk[r];
        const unsigned rvertex_number = iter_vertex_number[r];

        // Get the number of vertices in the r-th boundary r-th boundary
        // chunk
        const unsigned r_n_vertices =
          boundary_chunk_nvertices[rbnd_id][rbnd_chunk];

        // Is the r-th vertex an initial or final vertex
        bool is_initial = false;
        if (rvertex_number == 0)
        {
          is_initial = true;
        }
        bool is_final = false;
        if (rvertex_number == r_n_vertices - 1)
        {
          is_final = true;
        }

        if (is_initial)
        {
          // Set the base vertex info. of the r-th vertex reference

          // Mark the vertex as having assigned a base vertex
          base_vertices[rbnd_id][rbnd_chunk][0].has_base_vertex_assigned = true;
          // The base boundary id
          base_vertices[rbnd_id][rbnd_chunk][0].boundary_id = base_bnd_id;
          // The base boundary chunk number
          base_vertices[rbnd_id][rbnd_chunk][0].boundary_chunk = base_bnd_chunk;
          // The base vertex number
          base_vertices[rbnd_id][rbnd_chunk][0].vertex_number =
            base_vertex_number;
        }
        else if (is_final)
        {
          // Set the base vertex info. of the r-th vertex reference

          // Mark the vertex as having assigned a base vertex
          base_vertices[rbnd_id][rbnd_chunk][1].has_base_vertex_assigned = true;
          // The base boundary id
          base_vertices[rbnd_id][rbnd_chunk][1].boundary_id = base_bnd_id;
          // The base boundary chunk number
          base_vertices[rbnd_id][rbnd_chunk][1].boundary_chunk = base_bnd_chunk;
          // The base vertex number
          base_vertices[rbnd_id][rbnd_chunk][1].vertex_number =
            base_vertex_number;
        }

      } // for (i < n_vertex_references)

    } // for (v < n_vertices)
  }

#endif


  //======================================================================
  /// Move the nodes on boundaries with associated Geometric Objects so
  /// that they exactly coincide with the geometric object. This requires
  /// that the boundary coordinates are set up consistently
  //======================================================================
  void UnstructuredTwoDMeshGeometryBase::snap_nodes_onto_geometric_objects()
  {
    // Loop over all boundaries
    unsigned n_bound = this->nboundary();
    for (unsigned b = 0; b < n_bound; b++)
    {
      // Find the geometric object
      GeomObject* const geom_object_pt = this->boundary_geom_object_pt(b);

      // If there is one
      if (geom_object_pt != 0)
      {
        Vector<double> b_coord(1);
        Vector<double> new_x(2);
        const unsigned n_boundary_node = this->nboundary_node(b);
        for (unsigned n = 0; n < n_boundary_node; ++n)
        {
          // Get the boundary node and coordinates
          Node* const nod_pt = this->boundary_node_pt(b, n);
          nod_pt->get_coordinates_on_boundary(b, b_coord);

          // Get the position and time history according to the underlying
          // geometric object, assuming that it has the same timestepper
          // as the nodes....
          unsigned n_tvalues =
            1 + nod_pt->position_time_stepper_pt()->nprev_values();
          for (unsigned t = 0; t < n_tvalues; ++t)
          {
            // Get the position according to the underlying geometric object
            geom_object_pt->position(t, b_coord, new_x);

            // Move the node
            for (unsigned i = 0; i < 2; i++)
            {
              nod_pt->x(t, i) = new_x[i];
            }
          }
        }
      }
    }
  }

  //======================================================================
  /// Helper function that checks if a given point is inside a polygon
  //======================================================================
  bool UnstructuredTwoDMeshGeometryBase::is_point_inside_polygon_helper(
    Vector<Vector<double>> polygon_vertices, Vector<double> point)
  {
    // Total number of vertices (the first and last vertex should be the same)
    const unsigned n_vertex = polygon_vertices.size();

    // Check if internal point is actually located in bounding polygon
    // Reference: http://paulbourke.net/geometry/insidepoly/

    // Counter for number of intersections
    unsigned intersect_counter = 0;

    // Get first vertex
    Vector<double> p1 = polygon_vertices[0];
    for (unsigned i = 1; i <= n_vertex; i++)
    {
      // Get second vertex by wrap-around
      Vector<double> p2 = polygon_vertices[i % n_vertex];

      if (point[1] > std::min(p1[1], p2[1]))
      {
        if (point[1] <= std::max(p1[1], p2[1]))
        {
          if (point[0] <= std::max(p1[0], p2[0]))
          {
            if (p1[1] != p2[1])
            {
              double xintersect =
                (point[1] - p1[1]) * (p2[0] - p1[0]) / (p2[1] - p1[1]) + p1[0];
              if ((p1[0] == p2[0]) || (point[0] <= xintersect))
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
    if (intersect_counter % 2 == 0)
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  //======================================================================
  /// Gets the vertex number on the destination polyline (used /
  /// to create the connections among shared boundaries)
  //======================================================================
  const bool UnstructuredTwoDMeshGeometryBase::
    get_connected_vertex_number_on_destination_polyline(
      TriangleMeshPolyLine* dst_polyline_pt,
      Vector<double>& vertex_coordinates,
      unsigned& vertex_number)
  {
    // The returning flag indicating that the vertex was found in the
    // destination boundary
    bool found_vertex_number = false;

    // Get the number of vertices in the destination polyline
    const unsigned n_vertices = dst_polyline_pt->nvertex();

    // Loop over the vertices and return the closest vertex
    // to the given vertex coordinates
    for (unsigned i = 0; i < n_vertices; i++)
    {
      // Get the i-vertex on the polyline
      Vector<double> current_vertex = dst_polyline_pt->vertex_coordinate(i);

      // Compute the distance to the input vertex
      double dist = (vertex_coordinates[0] - current_vertex[0]) *
                      (vertex_coordinates[0] - current_vertex[0]) +
                    (vertex_coordinates[1] - current_vertex[1]) *
                      (vertex_coordinates[1] - current_vertex[1]);
      dist = sqrt(dist);

      // Have we found the vertex?
      if (dist < ToleranceForVertexMismatchInPolygons::Tolerable_error)
      {
        // Set the vertex index
        vertex_number = i;

        // Set the flag to found
        found_vertex_number = true;

        // Break the serching
        break;
      }
    } // for (i < n_vertices)

    // Return if the connection was found
    return found_vertex_number;
  }


#ifdef OOMPH_HAS_TRIANGLE_LIB

  //======================================================================
  /// Helper function to copy the connection information from
  /// the input curve(polyline or curviline) to the output polyline
  //======================================================================
  void UnstructuredTwoDMeshGeometryBase::copy_connection_information(
    TriangleMeshCurveSection* input_curve_pt,
    TriangleMeshCurveSection* output_curve_pt)
  {
    // Is there a connection to the initial vertex
    const bool initial_connection =
      input_curve_pt->is_initial_vertex_connected();

    // Is there a connection to the initial vertex
    const bool final_connection = input_curve_pt->is_final_vertex_connected();

    // If there are any connection at the ends then copy the information
    if (initial_connection || final_connection)
    {
      // Get the initial and final vertex from the input polyline
      Vector<double> input_initial_vertex(2);
      input_curve_pt->initial_vertex_coordinate(input_initial_vertex);
      Vector<double> input_final_vertex(2);
      input_curve_pt->final_vertex_coordinate(input_final_vertex);

      // Get the initial and final vertex from the output polyline
      Vector<double> output_initial_vertex(2);
      output_curve_pt->initial_vertex_coordinate(output_initial_vertex);
      Vector<double> output_final_vertex(2);
      output_curve_pt->final_vertex_coordinate(output_final_vertex);

      // Check if the polyline is in the same order or if it is
      // reversed. We need to know to copy the connection information
      // correctly

      // The flag to state if the copy is in the same order or in
      // reversed order
      bool copy_reversed = false;

      // Start with the initial vertices
      double dist_initial =
        ((output_initial_vertex[0] - input_initial_vertex[0]) *
         (output_initial_vertex[0] - input_initial_vertex[0])) +
        ((output_initial_vertex[1] - input_initial_vertex[1]) *
         (output_initial_vertex[1] - input_initial_vertex[1]));
      dist_initial = sqrt(dist_initial);

      // Compute the distance of the final vertices
      double dist_final = ((output_final_vertex[0] - input_final_vertex[0]) *
                           (output_final_vertex[0] - input_final_vertex[0])) +
                          ((output_final_vertex[1] - input_final_vertex[1]) *
                           (output_final_vertex[1] - input_final_vertex[1]));
      dist_final = sqrt(dist_final);

      // Get the distace from the input vertices to the output vertices
      // in the same order
      const double dist = dist_initial + dist_final;

      // Compute the distance between the input initial vertex and the
      // output final vertex
      double dist_initial_rev =
        ((input_initial_vertex[0] - output_final_vertex[0]) *
         (input_initial_vertex[0] - output_final_vertex[0])) +
        ((input_initial_vertex[1] - output_final_vertex[1]) *
         (input_initial_vertex[1] - output_final_vertex[1]));
      dist_initial_rev = sqrt(dist_initial_rev);

      // Compute the distance between the input final vertex and the
      // output initial vertex
      double dist_final_rev =
        ((input_final_vertex[0] - output_initial_vertex[0]) *
         (input_final_vertex[0] - output_initial_vertex[0])) +
        ((input_final_vertex[1] - output_initial_vertex[1]) *
         (input_final_vertex[1] - output_initial_vertex[1]));
      dist_final_rev = sqrt(dist_final_rev);

      // Get the distace from the input to the output vertices in
      // reversed order
      const double dist_rev = dist_initial_rev + dist_final_rev;

      // If the distance is smaller than the reversed distance then the
      // order of the vertices is the same
      if (dist <= dist_rev)
      {
        // Do nothing
        copy_reversed = false;
      }
      else if (dist_rev < dist)
      {
        // The connection information is copied in reversed order
        copy_reversed = true;
      } // else if (dist_rev < dist)

      // Copy the connection information
      // ----------------------------------------------------------------
      // Copy in reversed order? (copy normal)
      if (!copy_reversed)
      {
        // Initial vertex
        if (input_curve_pt->is_initial_vertex_connected())
        {
          output_curve_pt->set_initial_vertex_connected();

          output_curve_pt->initial_vertex_connected_bnd_id() =
            input_curve_pt->initial_vertex_connected_bnd_id();

          output_curve_pt->initial_vertex_connected_n_chunk() =
            input_curve_pt->initial_vertex_connected_n_chunk();

          // We need to know if we have to compute the vertex number or
          // if we need just to copy it
          if (input_curve_pt->is_initial_vertex_connected_to_curviline())
          {
            double initial_s_connection =
              input_curve_pt->initial_s_connection_value();

            unsigned bnd_id =
              output_curve_pt->initial_vertex_connected_bnd_id();

            double s_tolerance = input_curve_pt->tolerance_for_s_connection();

            output_curve_pt->initial_vertex_connected_n_vertex() =
              get_associated_vertex_to_svalue(
                initial_s_connection, bnd_id, s_tolerance);

            // The output polyline is not longer connected to a curviline
            // because we will be using the polyline representation of the
            // curviline
            output_curve_pt->unset_initial_vertex_connected_to_curviline();
          }
          else
          {
            output_curve_pt->initial_vertex_connected_n_vertex() =
              input_curve_pt->initial_vertex_connected_n_vertex();
          }

        } // if (input_curve_pt->is_initial_vertex_connected())

        // Final vertex
        if (input_curve_pt->is_final_vertex_connected())
        {
          output_curve_pt->set_final_vertex_connected();

          output_curve_pt->final_vertex_connected_bnd_id() =
            input_curve_pt->final_vertex_connected_bnd_id();

          output_curve_pt->final_vertex_connected_n_chunk() =
            input_curve_pt->final_vertex_connected_n_chunk();

          // We need to know if we have to compute the vertex number or
          // if we need just to copy it
          if (input_curve_pt->is_final_vertex_connected_to_curviline())
          {
            double final_s_connection =
              input_curve_pt->final_s_connection_value();

            unsigned bnd_id = input_curve_pt->final_vertex_connected_bnd_id();

            double s_tolerance = input_curve_pt->tolerance_for_s_connection();

            output_curve_pt->final_vertex_connected_n_vertex() =
              get_associated_vertex_to_svalue(
                final_s_connection, bnd_id, s_tolerance);

            // The output polyline is not longer connected to a curviline
            // because we will be using the polyline representation of the
            // curviline
            output_curve_pt->unset_final_vertex_connected_to_curviline();
          }
          else
          {
            output_curve_pt->final_vertex_connected_n_vertex() =
              input_curve_pt->final_vertex_connected_n_vertex();
          }

        } // if (input_curve_pt->is_final_vertex_connected())

      } // if (!copy_reversed)
      else
      {
        // Copy the connections in reversed order

        // Initial vertex
        if (input_curve_pt->is_initial_vertex_connected())
        {
          output_curve_pt->set_final_vertex_connected();

          output_curve_pt->final_vertex_connected_bnd_id() =
            input_curve_pt->initial_vertex_connected_bnd_id();

          output_curve_pt->final_vertex_connected_n_chunk() =
            input_curve_pt->initial_vertex_connected_n_chunk();

          // We need to know if we have to compute the vertex number or
          // if we need just to copy it
          if (input_curve_pt->is_initial_vertex_connected_to_curviline())
          {
            double initial_s_connection =
              input_curve_pt->initial_s_connection_value();

            unsigned bnd_id = input_curve_pt->initial_vertex_connected_bnd_id();

            double s_tolerance = input_curve_pt->tolerance_for_s_connection();

            output_curve_pt->final_vertex_connected_n_vertex() =
              get_associated_vertex_to_svalue(
                initial_s_connection, bnd_id, s_tolerance);

            // The output polyline is not longer connected to a curviline
            // because we will be using the polyline representation of the
            // curviline
            output_curve_pt->unset_final_vertex_connected_to_curviline();
          }
          else
          {
            output_curve_pt->final_vertex_connected_n_vertex() =
              input_curve_pt->initial_vertex_connected_n_vertex();
          }

        } // if (input_curve_pt->is_initial_vertex_connected())

        // Final vertex
        if (input_curve_pt->is_final_vertex_connected())
        {
          output_curve_pt->set_initial_vertex_connected();

          output_curve_pt->initial_vertex_connected_bnd_id() =
            input_curve_pt->final_vertex_connected_bnd_id();

          output_curve_pt->initial_vertex_connected_n_chunk() =
            input_curve_pt->final_vertex_connected_n_chunk();

          // We need to know if we have to compute the vertex number or
          // if we need just to copy it
          if (input_curve_pt->is_final_vertex_connected_to_curviline())
          {
            double final_s_connection =
              input_curve_pt->final_s_connection_value();

            unsigned bnd_id = input_curve_pt->final_vertex_connected_bnd_id();

            double s_tolerance = input_curve_pt->tolerance_for_s_connection();

            output_curve_pt->initial_vertex_connected_n_vertex() =
              get_associated_vertex_to_svalue(
                final_s_connection, bnd_id, s_tolerance);

            // The output polyline is not longer connected to a curviline
            // because we will be using the polyline representation of the
            // curviline
            output_curve_pt->unset_initial_vertex_connected_to_curviline();
          }
          else
          {
            output_curve_pt->initial_vertex_connected_n_vertex() =
              input_curve_pt->final_vertex_connected_n_vertex();
          }
        } //  if (input_curve_pt->is_final_vertex_connected())
      } // else if (!copy_reversed)
    } // if (initial_connection || final_connection)
  }

  //======================================================================
  /// Helper function to copy the connection information from
  /// the input curve(polyline or curviline) to the output sub-polyline
  //======================================================================
  void UnstructuredTwoDMeshGeometryBase::
    copy_connection_information_to_sub_polylines(
      TriangleMeshCurveSection* input_curve_pt,
      TriangleMeshCurveSection* output_curve_pt)
  {
    // Is there a connection to the initial vertex
    const bool initial_connection =
      input_curve_pt->is_initial_vertex_connected();

    // Is there a connection to the initial vertex
    const bool final_connection = input_curve_pt->is_final_vertex_connected();

    // If there are any connection at the ends then copy the information
    if (initial_connection || final_connection)
    {
      // Get the initial and final vertex from the input polyline
      Vector<double> input_initial_vertex(2);
      input_curve_pt->initial_vertex_coordinate(input_initial_vertex);
      Vector<double> input_final_vertex(2);
      input_curve_pt->final_vertex_coordinate(input_final_vertex);

      // Get the initial and final vertex from the output polyline
      Vector<double> output_initial_vertex(2);
      output_curve_pt->initial_vertex_coordinate(output_initial_vertex);
      Vector<double> output_final_vertex(2);
      output_curve_pt->final_vertex_coordinate(output_final_vertex);

      // Check if the polyline is in the same order or if it is
      // reversed. We need to know to copy the connection information
      // correctly

      // The flag to state if the copy is in the same order or in
      // reversed order
      bool copy_reversed = false;

      // Start with the initial vertices
      double dist_initial =
        ((output_initial_vertex[0] - input_initial_vertex[0]) *
         (output_initial_vertex[0] - input_initial_vertex[0])) +
        ((output_initial_vertex[1] - input_initial_vertex[1]) *
         (output_initial_vertex[1] - input_initial_vertex[1]));
      dist_initial = sqrt(dist_initial);

      // Compute the distance of the final vertices
      double dist_final = ((output_final_vertex[0] - input_final_vertex[0]) *
                           (output_final_vertex[0] - input_final_vertex[0])) +
                          ((output_final_vertex[1] - input_final_vertex[1]) *
                           (output_final_vertex[1] - input_final_vertex[1]));
      dist_final = sqrt(dist_final);

      // Get the distace from the input vertices to the output vertices
      // in the same order
      const double dist = dist_initial + dist_final;

      // Compute the distance between the input initial vertex and the
      // output final vertex
      double dist_initial_rev =
        ((input_initial_vertex[0] - output_final_vertex[0]) *
         (input_initial_vertex[0] - output_final_vertex[0])) +
        ((input_initial_vertex[1] - output_final_vertex[1]) *
         (input_initial_vertex[1] - output_final_vertex[1]));
      dist_initial_rev = sqrt(dist_initial_rev);

      // Compute the distance between the input final vertex and the
      // output initial vertex
      double dist_final_rev =
        ((input_final_vertex[0] - output_initial_vertex[0]) *
         (input_final_vertex[0] - output_initial_vertex[0])) +
        ((input_final_vertex[1] - output_initial_vertex[1]) *
         (input_final_vertex[1] - output_initial_vertex[1]));
      dist_final_rev = sqrt(dist_final_rev);

      // Get the distace from the input to the output vertices in
      // reversed order
      const double dist_rev = dist_initial_rev + dist_final_rev;

      // If the distance is smaller than the reversed distance then the
      // order of the vertices is the same
      if (dist <= dist_rev)
      {
        // Do nothing
        copy_reversed = false;
      }
      else if (dist_rev < dist)
      {
        // The connection information is copied in reversed order
        copy_reversed = true;
      } // else if (dist_rev < dist)

      // Copy the connection information
      // ----------------------------------------------------------------
      // Copy in reversed order? (copy normal)
      if (!copy_reversed)
      {
        // Initial vertex. We need to double check that the initial
        // vertex of the input curve section correspond with the
        // initial vertex of the output sub-polyline. It may be the
        // case that this sub-polyline has not the initial vertex
        if (input_curve_pt->is_initial_vertex_connected() &&
            dist_initial <
              ToleranceForVertexMismatchInPolygons::Tolerable_error)
        {
          output_curve_pt->set_initial_vertex_connected();

          output_curve_pt->initial_vertex_connected_bnd_id() =
            input_curve_pt->initial_vertex_connected_bnd_id();

          output_curve_pt->initial_vertex_connected_n_chunk() =
            input_curve_pt->initial_vertex_connected_n_chunk();

          // We need to know if we have to compute the vertex
          // number or if we need just to copy it
          if (input_curve_pt->is_initial_vertex_connected_to_curviline())
          {
            double initial_s_connection =
              input_curve_pt->initial_s_connection_value();

            unsigned bnd_id =
              output_curve_pt->initial_vertex_connected_bnd_id();

            double s_tolerance = input_curve_pt->tolerance_for_s_connection();

            output_curve_pt->initial_vertex_connected_n_vertex() =
              get_associated_vertex_to_svalue(
                initial_s_connection, bnd_id, s_tolerance);

            // The output polyline is not longer connected to a
            // curviline because we will be using the polyline
            // representation of the curviline
            output_curve_pt->unset_initial_vertex_connected_to_curviline();
          }
          else
          {
            output_curve_pt->initial_vertex_connected_n_vertex() =
              input_curve_pt->initial_vertex_connected_n_vertex();
          }
        } // input_curve_pt->is_initial_vertex_connected() &&
          // dist_initial<Tolerance

        // Final vertex. We need to double check that the final vertex
        // of the input curve section correspond with the final vertex
        // of the output sub-polyline. It may be the case that this
        // sub-polyline has not the final vertex
        if (input_curve_pt->is_final_vertex_connected() &&
            dist_final < ToleranceForVertexMismatchInPolygons::Tolerable_error)
        {
          output_curve_pt->set_final_vertex_connected();

          output_curve_pt->final_vertex_connected_bnd_id() =
            input_curve_pt->final_vertex_connected_bnd_id();

          output_curve_pt->final_vertex_connected_n_chunk() =
            input_curve_pt->final_vertex_connected_n_chunk();

          // We need to know if we have to compute the vertex number or
          // if we need just to copy it
          if (input_curve_pt->is_final_vertex_connected_to_curviline())
          {
            double final_s_connection =
              input_curve_pt->final_s_connection_value();

            unsigned bnd_id = input_curve_pt->final_vertex_connected_bnd_id();

            double s_tolerance = input_curve_pt->tolerance_for_s_connection();

            output_curve_pt->final_vertex_connected_n_vertex() =
              get_associated_vertex_to_svalue(
                final_s_connection, bnd_id, s_tolerance);

            // The output polyline is not longer connected to a
            // curviline because we will be using the polyline
            // representation of the curviline
            output_curve_pt->unset_final_vertex_connected_to_curviline();
          }
          else
          {
            output_curve_pt->final_vertex_connected_n_vertex() =
              input_curve_pt->final_vertex_connected_n_vertex();
          }
        } // input_curve_pt->is_final_vertex_connected() && dist_final<Tolerance
      } // if (!copy_reversed)
      else
      {
        // Copy the connections in reversed order

        // Initial vertex. We need to double check that the initial
        // vertex of the input curve section correspond with the
        // initial vertex of the output sub-polyline. It may be the
        // case that this sub-polyline has not the initial vertex
        if (input_curve_pt->is_initial_vertex_connected() &&
            dist_initial_rev <
              ToleranceForVertexMismatchInPolygons::Tolerable_error)
        {
          output_curve_pt->set_final_vertex_connected();

          output_curve_pt->final_vertex_connected_bnd_id() =
            input_curve_pt->initial_vertex_connected_bnd_id();

          output_curve_pt->final_vertex_connected_n_chunk() =
            input_curve_pt->initial_vertex_connected_n_chunk();

          // We need to know if we have to compute the vertex number or
          // if we need just to copy it
          if (input_curve_pt->is_initial_vertex_connected_to_curviline())
          {
            double initial_s_connection =
              input_curve_pt->initial_s_connection_value();

            unsigned bnd_id = input_curve_pt->initial_vertex_connected_bnd_id();

            double s_tolerance = input_curve_pt->tolerance_for_s_connection();

            output_curve_pt->final_vertex_connected_n_vertex() =
              get_associated_vertex_to_svalue(
                initial_s_connection, bnd_id, s_tolerance);

            // The output polyline is not longer connected to a
            // curviline because we will be using the polyline
            // representation of the curviline
            output_curve_pt->unset_final_vertex_connected_to_curviline();
          }
          else
          {
            output_curve_pt->final_vertex_connected_n_vertex() =
              input_curve_pt->initial_vertex_connected_n_vertex();
          }
        } // input_curve_pt->is_final_vertex_connected() && dist_final<Tolerance

        // Final vertex. We need to double check that the final vertex
        // of the input curve section correspond with the final vertex
        // of the output sub-polyline. It may be the case that this
        // sub-polyline has not the final vertex
        if (input_curve_pt->is_final_vertex_connected() &&
            dist_final_rev <
              ToleranceForVertexMismatchInPolygons::Tolerable_error)
        {
          output_curve_pt->set_initial_vertex_connected();

          output_curve_pt->initial_vertex_connected_bnd_id() =
            input_curve_pt->final_vertex_connected_bnd_id();

          output_curve_pt->initial_vertex_connected_n_chunk() =
            input_curve_pt->final_vertex_connected_n_chunk();

          // We need to know if we have to compute the vertex number or
          // if we need just to copy it
          if (input_curve_pt->is_final_vertex_connected_to_curviline())
          {
            double final_s_connection =
              input_curve_pt->final_s_connection_value();

            unsigned bnd_id = input_curve_pt->final_vertex_connected_bnd_id();

            double s_tolerance = input_curve_pt->tolerance_for_s_connection();

            output_curve_pt->initial_vertex_connected_n_vertex() =
              get_associated_vertex_to_svalue(
                final_s_connection, bnd_id, s_tolerance);

            // The output polyline is not longer connected to a
            // curviline because we will be using the polyline
            // representation of the curviline
            output_curve_pt->unset_initial_vertex_connected_to_curviline();
          }
          else
          {
            output_curve_pt->initial_vertex_connected_n_vertex() =
              input_curve_pt->final_vertex_connected_n_vertex();
          }
        } // input_curve_pt->is_final_vertex_connected() && dist_final<Tolerance
      } // else if (!copy_reversed)
    } // if (initial_connection || final_connection)
  }


  /// Public static flag to suppress warning; defaults to true because
  /// it really clutters up the output. It's unlikely to ever be a
  /// genuine error...
  bool UnstructuredTwoDMeshGeometryBase::
    Suppress_warning_about_regions_and_boundaries = true;


#endif // OOMPH_HAS_TRIANGLE_LIB


} // namespace oomph
