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
#include "triangle_scaffold_mesh.h"


namespace oomph
{
  //=====================================================================
  /// Check mesh integrity -- performs some internal consistency checks
  /// and throws error if violated.
  //=====================================================================
  void TriangleScaffoldMesh::check_mesh_integrity()
  {
    // Check that all nodes are attached to elements
    std::map<Node*, bool> done;

    // Number of nodes per element
    unsigned nnod_el = finite_element_pt(0)->nnode();

    // Loop over elements to visit their nodes
    unsigned nelem = nelement();
    for (unsigned e = 0; e < nelem; e++)
    {
      // Loop over all nodes in element
      for (unsigned j = 0; j < nnod_el; j++)
      {
        // Pointer to node in the scaffold mesh
        Node* node_pt = finite_element_pt(e)->node_pt(j);

        done[node_pt] = true;
      }
    }

    // Now check if all nodes have been visited
    std::ostringstream error_stream;
    bool broken = false;
    unsigned nnod = nnode();
    for (unsigned j = 0; j < nnod; j++)
    {
      if (!done[Node_pt[j]])
      {
        error_stream << "Scaffold node " << j
                     << " does not seem to be attached to any element";
        broken = true;
      }
    }
    if (broken)
    {
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  //=====================================================================
  /// Constructor: Pass the filenames of the triangle files
  /// The assumptions are that the nodes have been assigned boundary
  /// information which is used in the nodal construction to ensure that
  /// BoundaryNodes are indeed constructed when necessary. Additional
  /// boundary information is added from the segment boundaries.
  //=====================================================================
  TriangleScaffoldMesh::TriangleScaffoldMesh(const std::string& node_file_name,
                                             const std::string& ele_file_name,
                                             const std::string& poly_file_name)
  {
    // Process element file
    //---------------------
    std::ifstream element_file(ele_file_name.c_str(), std::ios_base::in);

    // Check that the file actually opened correctly
    if (!element_file.is_open())
    {
      std::string error_msg("Failed to open element file: ");
      error_msg += "\"" + ele_file_name + "\".";
      throw OomphLibError(
        error_msg, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Number of elements
    unsigned n_element;
    element_file >> n_element;

    // Number of nodes per element
    unsigned n_local_node;
    element_file >> n_local_node;
    if (n_local_node != 3)
    {
      std::ostringstream error_stream;
      error_stream
        << "Triangle should only be used to generate 3-noded triangles!\n"
        << "Your triangle input file, contains data for " << n_local_node
        << "-noded triangles" << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Dummy nodal attributes
    unsigned dummy_attribute;

    // Element attributes may be used if we have internal boundaries
    Element_attribute.resize(n_element, 0.0);

    // Dummy storage for element numbers
    unsigned dummy_element_number;

    // Resize stoorage for global node numbers listed element-by-element
    Global_node.resize(n_element * n_local_node);

#ifdef PARANOID
    std::map<unsigned, bool> global_node_done;
#endif

    // Initialise counter
    unsigned k = 0;

    // Are attributes specified?
    unsigned attribute_flag;
    element_file >> attribute_flag;

    // Read global node numbers for all elements
    if (attribute_flag == 0)
    {
      for (unsigned i = 0; i < n_element; i++)
      {
        element_file >> dummy_element_number;
        for (unsigned j = 0; j < n_local_node; j++)
        {
          element_file >> Global_node[k];
#ifdef PARANOID
          global_node_done[Global_node[k] - 1] = true;
#endif
          k++;
        }
      }
    }
    else
    {
      for (unsigned i = 0; i < n_element; i++)
      {
        element_file >> dummy_element_number;
        for (unsigned j = 0; j < n_local_node; j++)
        {
          element_file >> Global_node[k];
#ifdef PARANOID
          global_node_done[Global_node[k] - 1] = true;
#endif
          k++;
        }
        element_file >> Element_attribute[i];
      }
    }
    element_file.close();

#ifdef PARANOID
    // Determine if the node numbering starts at 0 or 1 (triangle can do
    // either depending on input arguments). We can't (currently) handle
    // 0-based indexing so throw an error if it is being used.
    if (*std::min_element(Global_node.begin(), Global_node.end()) != 1)
    {
      std::string err("Triangle is using 0-based indexing, ");
      err += "however the oomph-lib interface can only handle 1-based indexing "
             "in Triangle.\n";
      err +=
        "This is likely to be due to the input file using 0-based indexing.";
      err += "Alternatively you may have specified the -z flag when running "
             "Triangle.";
      throw OomphLibError(
        err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Resize the Element vector
    Element_pt.resize(n_element);

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
    unsigned n_node;
    node_file >> n_node;

    // Create a vector of boolean so as not to create the same node twice
    std::vector<bool> done(n_node, false);

    // Resize the Node vector
    Node_pt.resize(n_node);

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
    node_file >> attribute_flag;

    // Flag for boundary markers
    unsigned boundary_markers_flag;
    node_file >> boundary_markers_flag;

    // Dummy for node number
    unsigned dummy_node_number;

    // Create storage for nodal posititions and boundary markers
    Vector<double> x_node(n_node);
    Vector<double> y_node(n_node);
    Vector<unsigned> bound(n_node);

    // Check if there are attributes
    if (attribute_flag == 0)
    {
      if (boundary_markers_flag == 1)
      {
        for (unsigned i = 0; i < n_node; i++)
        {
          node_file >> dummy_node_number;
          node_file >> x_node[i];
          node_file >> y_node[i];
          node_file >> bound[i];
        }
        node_file.close();
      }
      else
      {
        for (unsigned i = 0; i < n_node; i++)
        {
          node_file >> dummy_node_number;
          node_file >> x_node[i];
          node_file >> y_node[i];
          bound[i] = 0;
        }
        node_file.close();
      }
    }
    else
    {
      if (boundary_markers_flag == 1)
      {
        for (unsigned i = 0; i < n_node; i++)
        {
          node_file >> dummy_node_number;
          node_file >> x_node[i];
          node_file >> y_node[i];
          node_file >> dummy_attribute;
          node_file >> bound[i];
        }
        node_file.close();
      }
      else
      {
        for (unsigned i = 0; i < n_node; i++)
        {
          node_file >> dummy_node_number;
          node_file >> x_node[i];
          node_file >> y_node[i];
          node_file >> dummy_attribute;
          bound[i] = 0;
        }
        node_file.close();
      }
    } // end


    // Determine highest boundary index
    // --------------------------------
    unsigned n_bound = 0;
    if (boundary_markers_flag == 1)
    {
      n_bound = bound[0];
      for (unsigned i = 1; i < n_node; i++)
      {
        if (bound[i] > n_bound)
        {
          n_bound = bound[i];
        }
      }
    }

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

    // Number of nodes in poly file
    unsigned n_node_poly;
    poly_file >> n_node_poly;

    // Dimension
    poly_file >> dimension;

    // Attribute flag
    poly_file >> attribute_flag;

    // Boundary markers flag
    poly_file >> boundary_markers_flag;


    // Ignore node information: Note: No, we can't extract the
    // actual nodes themselves from here!
    unsigned dummy;
    if (n_node_poly > 0)
    {
      if (attribute_flag == 0)
      {
        if (boundary_markers_flag == 1)
        {
          for (unsigned i = 0; i < n_node_poly; i++)
          {
            poly_file >> dummy;
            poly_file >> dummy;
            poly_file >> dummy;
            poly_file >> dummy;
          }
        }
        else
        {
          for (unsigned i = 0; i < n_node_poly; i++)
          {
            poly_file >> dummy;
            poly_file >> dummy;
            poly_file >> dummy;
          }
        }
      }
      else
      {
        if (boundary_markers_flag == 1)
        {
          for (unsigned i = 0; i < n_node_poly; i++)
          {
            poly_file >> dummy;
            poly_file >> dummy;
            poly_file >> dummy;
            poly_file >> dummy;
            poly_file >> dummy;
          }
        }
        else
        {
          for (unsigned i = 0; i < n_node_poly; i++)
          {
            poly_file >> dummy;
            poly_file >> dummy;
            poly_file >> dummy;
            poly_file >> dummy;
          }
        }
      }
    }

    // Now extract the segment information
    // Segements are lines that lie on boundaries of the domain
    //----------------------------------------------------------

    // Number of segments
    unsigned n_segment;
    poly_file >> n_segment;

    // Boundary marker flag
    poly_file >> boundary_markers_flag;

    // Storage for the global node numbers (in the triangle 1-based
    // numbering scheme!) of the first and second node in each segment
    Vector<unsigned> first_node(n_segment);
    Vector<unsigned> second_node(n_segment);

    // Storage for the boundary marker for each segment
    Vector<unsigned> segment_boundary(n_segment);

    // Dummy for global segment number
    unsigned dummy_segment_number;

    // Storage for the edges associated with each node. Nodes are indexed
    // using the Triangle 1-based index which is why there is a +1 here.
    Vector<std::set<unsigned>> node_on_edges(n_node + 1);

    // Extract information for each segment
    for (unsigned i = 0; i < n_segment; i++)
    {
      poly_file >> dummy_segment_number;
      poly_file >> first_node[i];
      poly_file >> second_node[i];
      poly_file >> segment_boundary[i];
      // Check that we don't have a higher segment boundary number
      if (segment_boundary[i] > n_bound)
      {
        n_bound = segment_boundary[i];
      }
      // Add the segment index to each node
      node_on_edges[first_node[i]].insert(i);
      node_on_edges[second_node[i]].insert(i);
    }

    // Extract hole center information
    unsigned nhole = 0;
    poly_file >> nhole;
    if (nhole != 0)
    {
      Hole_centre.resize(nhole);

      // Dummy for hole number
      unsigned dummy_hole;
      // Loop over the holes to get centre coords
      for (unsigned ihole = 0; ihole < nhole; ihole++)
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

    // Set the number of boundaries
    if (n_bound > 0)
    {
      this->set_nboundary(n_bound);
    }


    // Create the elements
    //--------------------

    // Counter for nodes in the vector that lists
    // the global node numbers of the elements' local nodes
    unsigned counter = 0;
    for (unsigned e = 0; e < n_element; e++)
    {
      Element_pt[e] = new TElement<2, 2>;
      for (unsigned j = 0; j < n_local_node; j++)
      {
        unsigned global_node_number = Global_node[counter];
        if (done[global_node_number - 1] == false) //... -1 because node number
        // begins at 1 in triangle
        {
          // If we are on the boundary
          if ((boundary_markers_flag == 1) &&
              (bound[global_node_number - 1] > 0))
          {
            // Construct a boundary node
            Node_pt[global_node_number - 1] =
              finite_element_pt(e)->construct_boundary_node(j);
            // Add to the boundary node look-up scheme
            add_boundary_node(bound[global_node_number - 1] - 1,
                              Node_pt[global_node_number - 1]);
          }
          // Otherwise make an ordinary node
          else
          {
            Node_pt[global_node_number - 1] =
              finite_element_pt(e)->construct_node(j);
          }
          done[global_node_number - 1] = true;
          Node_pt[global_node_number - 1]->x(0) =
            x_node[global_node_number - 1];
          Node_pt[global_node_number - 1]->x(1) =
            y_node[global_node_number - 1];
        }
        else
        {
          finite_element_pt(e)->node_pt(j) = Node_pt[global_node_number - 1];
        }
        counter++;
      }
    }

    // Resize the "matrix" that stores the boundary id for each
    // edge in each element.
    Edge_boundary.resize(n_element);
    Edge_index.resize(n_element);

    // Storage for the global node numbers (in triangle's 1-based
    // numbering scheme) for the zero-th, 1st, and 2nd node in each
    // triangle.
    unsigned glob_num[3] = {0, 0, 0};

    // 0-based index used to construct a global index-based lookup scheme
    // for each edge that will be used to uniquely construct mid-side
    // nodes.
    // The segments (edges that lie on boundaries) have already
    // been added to the scheme, so we start with the number of segments.
    Nglobal_edge = n_segment;

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Each element has three edges
      Edge_boundary[e].resize(3);
      Edge_index[e].resize(3);
      // By default each edge is NOT on a boundary
      for (unsigned i = 0; i < 3; i++)
      {
        Edge_boundary[e][i] = 0;
      }

      // Read out the global node numbers from the triangle data structure
      const unsigned element_offset = e * n_local_node;
      for (unsigned i = 0; i < 3; i++)
      {
        glob_num[i] = Global_node[element_offset + i];
      }

      // Now we know the global node numbers of the elements' three nodes
      // in triangle's 1-based numbering.

      // Determine whether any of the elements edges have already been
      // allocated an index. This may be because they are on boundaries
      // (segments) or because they have already occured.
      // The global index for the i-th edge will be stored in edge_index[i]
      for (unsigned i = 0; i < 3; i++)
      {
        std::vector<unsigned> local_edge_index;

        // Find the common global edge index for the nodes on
        // the i-th element edge (note the use of moular arithmetic here)
        std::set_intersection(node_on_edges[glob_num[i]].begin(),
                              node_on_edges[glob_num[i]].end(),
                              node_on_edges[glob_num[(i + 1) % 3]].begin(),
                              node_on_edges[glob_num[(i + 1) % 3]].end(),
                              std::insert_iterator<std::vector<unsigned>>(
                                local_edge_index, local_edge_index.begin()));

        // If the nodes share more than one global edge index, then
        // we have a problem
        if (local_edge_index.size() > 1)
        {
          throw OomphLibError(
            "Nodes in scaffold mesh share more than one global edge",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }

        // If the element's edge is not already allocated, the intersection
        // will be empty
        if (local_edge_index.size() == 0)
        {
          // Allocate the next global index
          Edge_index[e][i] = Nglobal_edge;
          // Associate the new edge index with the nodes
          node_on_edges[glob_num[i]].insert(Nglobal_edge);
          node_on_edges[glob_num[(i + 1) % 3]].insert(Nglobal_edge);
          // Increment the global edge index
          ++Nglobal_edge;
        }
        // Otherwise we already have an edge
        else if (local_edge_index.size() == 1)
        {
          // Set the edge index
          Edge_index[e][i] = local_edge_index[0];
          // Allocate the boundary index, if it is a segment
          if (local_edge_index[0] < n_segment)
          {
            Edge_boundary[e][i] = segment_boundary[local_edge_index[0]];
            // Add the nodes to the boundary look-up scheme in
            // oomph-lib (0-based) index
            add_boundary_node(segment_boundary[local_edge_index[0]] - 1,
                              Node_pt[glob_num[i] - 1]);
            add_boundary_node(segment_boundary[local_edge_index[0]] - 1,
                              Node_pt[glob_num[(i + 1) % 3] - 1]);
          }
        }
      }
    }

#ifdef PARANOID

    std::ostringstream error_stream;
    bool broken = false;
    unsigned nnod = nnode();
    error_stream << "Checking presence of " << nnod << " global nodes\n";
    for (unsigned j = 0; j < nnod; j++)
    {
      if (!global_node_done[j])
      {
        error_stream << "Global node " << j
                     << " was not listed in *.ele file\n";
        broken = true;
      }
    }
    if (broken)
    {
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check things and throw if mesh is broken...
    check_mesh_integrity();
#endif
  }

#ifdef OOMPH_HAS_TRIANGLE_LIB

  //=====================================================================
  ///  Constructor: Pass a data structure obtained from the triangulate
  /// function
  //=====================================================================
  TriangleScaffoldMesh::TriangleScaffoldMesh(TriangulateIO& triangle_data)
  {
    // Number of elements
    unsigned n_element = static_cast<unsigned>(triangle_data.numberoftriangles);

    // Number of nodes per element
    unsigned n_local_node =
      static_cast<unsigned>(triangle_data.numberofcorners);
    if (n_local_node != 3)
    {
      std::ostringstream error_stream;
      error_stream
        << "Triangle should only be used to generate 3-noded triangles!\n"
        << "Your triangle input file, contains data for " << n_local_node
        << "-noded triangles" << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Element attributes may be used if we have internal boundaries
    Element_attribute.resize(n_element, 0.0);

    // Resizestoorage for global node numbers listed element-by-element
    Global_node.resize(n_element * n_local_node);

#ifdef PARANOID
    std::map<unsigned, bool> global_node_done;
#endif

    // Initialise counter
    unsigned k = 0;

    // Are attributes specified?
    unsigned attribute_flag =
      static_cast<unsigned>(triangle_data.numberoftriangleattributes);

    // Read global node numbers for all elements
    if (attribute_flag == 0)
    {
      for (unsigned i = 0; i < n_element; i++)
      {
        for (unsigned j = 0; j < n_local_node; j++)
        {
          Global_node[k] = static_cast<unsigned>(triangle_data.trianglelist[k]);
#ifdef PARANOID
          global_node_done[Global_node[k] - 1] = true;
#endif
          k++;
        }
      }
    }
    else
    {
      for (unsigned i = 0; i < n_element; i++)
      {
        for (unsigned j = 0; j < n_local_node; j++)
        {
          Global_node[k] = static_cast<unsigned>(triangle_data.trianglelist[k]);
#ifdef PARANOID
          global_node_done[Global_node[k] - 1] = true;
#endif
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
    std::vector<bool> done(n_node, false);

    // Resize the Node vector
    Node_pt.resize(n_node);

    // Flag for boundary markers
    unsigned boundary_markers_flag = 0;
    if (triangle_data.pointmarkerlist != 0)
    {
      boundary_markers_flag = 1;
    }

    // Create storage for nodal posititions and boundary markers
    Vector<double> x_node(n_node);
    Vector<double> y_node(n_node);
    Vector<unsigned> bound(n_node);

    // We shall ingnore all point attributes
    if (boundary_markers_flag == 1)
    {
      for (unsigned i = 0; i < n_node; i++)
      {
        x_node[i] = triangle_data.pointlist[2 * i];
        y_node[i] = triangle_data.pointlist[2 * i + 1];
        bound[i] = static_cast<unsigned>(triangle_data.pointmarkerlist[i]);
      }
    }
    else
    {
      for (unsigned i = 0; i < n_node; i++)
      {
        x_node[i] = triangle_data.pointlist[2 * i];
        y_node[i] = triangle_data.pointlist[2 * i + 1];
        bound[i] = 0;
      }
    }

    // Determine highest boundary index
    // --------------------------------
    unsigned n_bound = 0;
    if (boundary_markers_flag == 1)
    {
      n_bound = bound[0];
      for (unsigned i = 1; i < n_node; i++)
      {
        if (bound[i] > n_bound)
        {
          n_bound = bound[i];
        }
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

    // Storage for the edges associated with each node. Nodes are indexed
    // using the Triangle 1-based index which is why there is a +1 here.
    Vector<std::set<unsigned>> node_on_edges(n_node + 1);

    // Extract information for each segment
    for (unsigned i = 0; i < n_segment; i++)
    {
      first_node[i] = static_cast<unsigned>(triangle_data.segmentlist[2 * i]);
      second_node[i] =
        static_cast<unsigned>(triangle_data.segmentlist[2 * i + 1]);
      segment_boundary[i] =
        static_cast<unsigned>(triangle_data.segmentmarkerlist[i]);
      // Check that we don't have a higher segment boundary number
      if (segment_boundary[i] > n_bound)
      {
        n_bound = segment_boundary[i];
      }
      // Add the segment index to each node
      node_on_edges[first_node[i]].insert(i);
      node_on_edges[second_node[i]].insert(i);
    }

    // Extract hole center information
    unsigned nhole = triangle_data.numberofholes;
    if (nhole != 0)
    {
      Hole_centre.resize(nhole);

      // Coords counter
      unsigned count_coords = 0;

      // Loop over the holes to get centre coords
      for (unsigned ihole = 0; ihole < nhole; ihole++)
      {
        Hole_centre[ihole].resize(2);

        // Read the centre value
        Hole_centre[ihole][0] = triangle_data.holelist[count_coords];
        Hole_centre[ihole][1] = triangle_data.holelist[count_coords + 1];

        // Increment coords counter
        count_coords += 2;
      }
    }
    else
    {
      Hole_centre.resize(0);
    }

    // Set number of boundaries
    if (n_bound > 0)
    {
      set_nboundary(n_bound);
    }


    // Create the elements
    //--------------------

    // Counter for nodes in the vector that lists
    // the global node numbers of the elements' local nodes
    unsigned counter = 0;
    for (unsigned e = 0; e < n_element; e++)
    {
      Element_pt[e] = new TElement<2, 2>;
      for (unsigned j = 0; j < n_local_node; j++)
      {
        unsigned global_node_number = Global_node[counter];
        if (done[global_node_number - 1] == false) //... -1 because node number
        // begins at 1 in triangle
        {
          // If we are on the boundary
          if ((boundary_markers_flag == 1) &&
              (bound[global_node_number - 1] > 0))
          {
            // Construct a boundary node
            Node_pt[global_node_number - 1] =
              finite_element_pt(e)->construct_boundary_node(j);

            // Add to the boundary node look-up scheme
            add_boundary_node(bound[global_node_number - 1] - 1,
                              Node_pt[global_node_number - 1]);
          }
          // Otherwise make an ordinary node
          else
          {
            Node_pt[global_node_number - 1] =
              finite_element_pt(e)->construct_node(j);
          }
          done[global_node_number - 1] = true;
          Node_pt[global_node_number - 1]->x(0) =
            x_node[global_node_number - 1];
          Node_pt[global_node_number - 1]->x(1) =
            y_node[global_node_number - 1];
        }
        else
        {
          finite_element_pt(e)->node_pt(j) = Node_pt[global_node_number - 1];
        }
        counter++;
      }
    }

    // Resize the "matrix" that stores the boundary id for each
    // edge in each element.
    Edge_boundary.resize(n_element);
    Edge_index.resize(n_element);

    // Storage for the global node numbers (in triangle's 1-based
    // numbering scheme) for the zero-th, 1st, and 2nd node in each
    // triangle.
    unsigned glob_num[3] = {0, 0, 0};

    // 0-based index used to construct a global index-based lookup scheme
    // for each edge that will be used to uniquely construct mid-side
    // nodes.
    // The segments (edges that lie on boundaries) have already
    // been added to the scheme, so we start with the number of segments.
    Nglobal_edge = n_segment;

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Each element has three edges
      Edge_boundary[e].resize(3);
      Edge_index[e].resize(3);
      // By default each edge is NOT on a boundary
      for (unsigned i = 0; i < 3; i++)
      {
        Edge_boundary[e][i] = 0;
      }

      // Read out the global node numbers from the triangle data structure
      const unsigned element_offset = e * n_local_node;
      for (unsigned i = 0; i < 3; i++)
      {
        glob_num[i] = Global_node[element_offset + i];
      }

      // Now we know the global node numbers of the elements' three nodes
      // in triangle's 1-based numbering.

      // Determine whether any of the elements edges have already been
      // allocated an index. This may be because they are on boundaries
      // (segments) or because they have already occured.
      // The global index for the i-th edge will be stored in edge_index[i]
      for (unsigned i = 0; i < 3; i++)
      {
        std::vector<unsigned> local_edge_index;

        // Find the common global edge index for the nodes on
        // the i-th element edge (note the use of moular arithmetic here)
        std::set_intersection(node_on_edges[glob_num[i]].begin(),
                              node_on_edges[glob_num[i]].end(),
                              node_on_edges[glob_num[(i + 1) % 3]].begin(),
                              node_on_edges[glob_num[(i + 1) % 3]].end(),
                              std::insert_iterator<std::vector<unsigned>>(
                                local_edge_index, local_edge_index.begin()));

        // If the nodes share more than one global edge index, then
        // we have a problem
        if (local_edge_index.size() > 1)
        {
          throw OomphLibError(
            "Nodes in scaffold mesh share more than one global edge",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }


        // If the element's edge is not already allocated, the intersection
        // will be empty
        if (local_edge_index.size() == 0)
        {
          // Allocate the next global index
          Edge_index[e][i] = Nglobal_edge;
          // Associate the new edge index with the nodes
          node_on_edges[glob_num[i]].insert(Nglobal_edge);
          node_on_edges[glob_num[(i + 1) % 3]].insert(Nglobal_edge);
          // Increment the global edge index
          ++Nglobal_edge;
        }
        // Otherwise we already have an edge
        else if (local_edge_index.size() == 1)
        {
          // Set the edge index
          Edge_index[e][i] = local_edge_index[0];
          // Allocate the boundary index, if it is a segment
          if (local_edge_index[0] < n_segment)
          {
            Edge_boundary[e][i] = segment_boundary[local_edge_index[0]];
            // Add the nodes to the boundary look-up scheme in
            // oomph-lib (0-based) index
            add_boundary_node(segment_boundary[local_edge_index[0]] - 1,
                              Node_pt[glob_num[i] - 1]);
            add_boundary_node(segment_boundary[local_edge_index[0]] - 1,
                              Node_pt[glob_num[(i + 1) % 3] - 1]);
          }
        }
      }
    }

#ifdef PARANOID

    std::ostringstream error_stream;
    bool broken = false;
    unsigned nnod = nnode();
    error_stream << "Checking presence of " << nnod << " global nodes\n";
    for (unsigned j = 0; j < nnod; j++)
    {
      if (!global_node_done[j])
      {
        error_stream << "Global node " << j
                     << " was not listed in *.ele file\n";
        broken = true;
      }
    }
    if (broken)
    {
      error_stream
        << "This error means that some of the nodes are not connected \n"
        << " to (bulk) elements. This can happen if there is an isolated\n"
        << " boundary line in the mesh. One possible cause for this is\n"
        << " specifying a hole coordinate in the wrong place so that there\n"
        << " a gap between the mesh and the outer boundary.\n";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check things and throw if mesh is broken...
    check_mesh_integrity();
#endif
  }

#endif

} // namespace oomph
