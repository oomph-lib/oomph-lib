// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
#include "mesh.h"
#include "Telements.h"
#include "tetgen_scaffold_mesh.h"

namespace oomph
{
  //======================================================================
  /// Constructor: Pass the filename of the tetrahedra file
  /// The assumptions are that the nodes have been assigned boundary
  /// information which is used in the nodal construction to make sure
  /// that BoundaryNodes are constructed. Any additional boundaries are
  /// determined from the face boundary information.
  //======================================================================
  TetgenScaffoldMesh::TetgenScaffoldMesh(const std::string& node_file_name,
                                         const std::string& element_file_name,
                                         const std::string& face_file_name)
  {
    // Process the element file
    // --------------------------
    std::ifstream element_file(element_file_name.c_str(), std::ios_base::in);

    // Check that the file actually opened correctly
#ifdef PARANOID
    if (!element_file.is_open())
    {
      std::string error_msg("Failed to open element file: ");
      error_msg += "\"" + element_file_name + "\".";
      throw OomphLibError(
        error_msg, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Read in number of elements
    unsigned n_element;
    element_file >> n_element;

    // Read in number of nodes per element
    unsigned n_local_node;
    element_file >> n_local_node;

    // Throw an error if we have anything but linear simplices
    if (n_local_node != 4)
    {
      std::ostringstream error_stream;
      error_stream
        << "Tetgen should only be used to generate 4-noded tetrahedra!\n"
        << "Your tetgen input file, contains data for " << n_local_node
        << "-noded tetrahedra" << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Dummy nodal attributes
    unsigned dummy_attribute;

    // Element attributes may be used to distinguish internal regions
    // NOTE: This stores doubles because tetgen forces us to!
    Element_attribute.resize(n_element, 0.0);

    // Dummy storage for element numbers
    unsigned dummy_element_number;

    // Resize storage for the global node numbers listed element-by-element
    Global_node.resize(n_element * n_local_node);

    // Initialise (global) node counter
    unsigned k = 0;
    // Are there attributes?
    unsigned attribute_flag;
    element_file >> attribute_flag;

    // If no attributes
    if (attribute_flag == 0)
    {
      // Read in global node numbers
      for (unsigned i = 0; i < n_element; i++)
      {
        element_file >> dummy_element_number;
        for (unsigned j = 0; j < n_local_node; j++)
        {
          element_file >> Global_node[k];
          k++;
        }
      }
    }
    // Otherwise read in the attributes as well
    else
    {
      for (unsigned i = 0; i < n_element; i++)
      {
        element_file >> dummy_element_number;
        for (unsigned j = 0; j < n_local_node; j++)
        {
          element_file >> Global_node[k];
          k++;
        }
        element_file >> Element_attribute[i];
      }
    }
    element_file.close();

    // Resize the Element vector
    Element_pt.resize(n_element);

    // Process node file
    //--------------------
    std::ifstream node_file(node_file_name.c_str(), std::ios_base::in);

    // Check that the file actually opened correctly
#ifdef PARANOID
    if (!node_file.is_open())
    {
      std::string error_msg("Failed to open node file: ");
      error_msg += "\"" + node_file_name + "\".";
      throw OomphLibError(
        error_msg, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Read in the number of nodes
    unsigned n_node;
    node_file >> n_node;

    // Create a vector of boolean so as not to create the same node twice
    std::vector<bool> done(n_node, false);

    // Resize the Node vector
    Node_pt.resize(n_node);

    // Set the spatial dimension of the nodes
    unsigned dimension;
    node_file >> dimension;

#ifdef PARANOID
    if (dimension != 3)
    {
      throw OomphLibError("The dimesion of the nodes must be 3\n",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Flag for attributes
    node_file >> attribute_flag;

    // Flag for boundary markers
    unsigned boundary_markers_flag;
    node_file >> boundary_markers_flag;

    // Dummy storage for the node number
    unsigned dummy_node_number;

    // Create storage for nodal positions and boundary markers
    Vector<double> x_node(n_node);
    Vector<double> y_node(n_node);
    Vector<double> z_node(n_node);
    Vector<unsigned> bound(n_node);

    // Check if there are attributes
    // If not
    if (attribute_flag == 0)
    {
      // Are there boundary markers
      if (boundary_markers_flag == 1)
      {
        for (unsigned i = 0; i < n_node; i++)
        {
          node_file >> dummy_node_number;
          node_file >> x_node[i];
          node_file >> y_node[i];
          node_file >> z_node[i];
          node_file >> bound[i];
        }
        node_file.close();
      }
      // Otherwise no boundary markers
      else
      {
        for (unsigned i = 0; i < n_node; i++)
        {
          node_file >> dummy_node_number;
          node_file >> x_node[i];
          node_file >> y_node[i];
          node_file >> z_node[i];
          bound[i] = 0;
        }
        node_file.close();
      }
    }
    // Otherwise there are attributes
    else
    {
      if (boundary_markers_flag == 1)
      {
        for (unsigned i = 0; i < n_node; i++)
        {
          node_file >> dummy_node_number;
          node_file >> x_node[i];
          node_file >> y_node[i];
          node_file >> z_node[i];
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
          node_file >> z_node[i];
          node_file >> dummy_attribute;
          bound[i] = 0;
        }
        node_file.close();
      }
    }

    // Determine highest boundary index
    //------------------------------------
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

    // Process face file to extract boundary faces
    //--------------------------------------------

    // Open face file
    std::ifstream face_file(face_file_name.c_str(), std::ios_base::in);

    // Check that the file actually opened correctly
#ifdef PARANOID
    if (!face_file.is_open())
    {
      std::string error_msg("Failed to open face file: ");
      error_msg += "\"" + face_file_name + "\".";
      throw OomphLibError(
        error_msg, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Number of faces in face file
    unsigned n_face;
    face_file >> n_face;

    // Boundary markers flag
    face_file >> boundary_markers_flag;

    // Storage for the global node numbers (in the tetgen 1-based
    // numbering scheme!) of the first, second and third  node in
    //  each segment
    Vector<unsigned> first_node(n_face);
    Vector<unsigned> second_node(n_face);
    Vector<unsigned> third_node(n_face);

    // Storage for the boundary marker for each face
    Vector<unsigned> face_boundary(n_face);

    // Dummy for global face number
    unsigned dummy_face_number;

    // Storage for the (boundary) faces associated with each node.
    // Nodes are indexed using Tetgen's 1-based scheme, which is why
    // there is a +1 here
    Vector<std::set<unsigned>> node_on_faces(n_node + 1);

    // Extract information for each segment
    for (unsigned i = 0; i < n_face; i++)
    {
      face_file >> dummy_face_number;
      face_file >> first_node[i];
      face_file >> second_node[i];
      face_file >> third_node[i];
      face_file >> face_boundary[i];
      if (face_boundary[i] > n_bound)
      {
        n_bound = face_boundary[i];
      }
      // Add the face index to each node
      node_on_faces[first_node[i]].insert(i);
      node_on_faces[second_node[i]].insert(i);
      node_on_faces[third_node[i]].insert(i);
    }
    face_file.close();

    // Set number of boundaries
    if (n_bound > 0)
    {
      this->set_nboundary(n_bound);
    }

    // Create the elements
    unsigned counter = 0;
    for (unsigned e = 0; e < n_element; e++)
    {
      Element_pt[e] = new TElement<3, 2>;
      unsigned global_node_number = Global_node[counter];
      if (done[global_node_number - 1] == false)
      // ... -1 because node number begins at 1 in tetgen
      {
        // If the node is on a boundary, construct a boundary node
        if ((boundary_markers_flag == 1) && (bound[global_node_number - 1] > 0))
        {
          // Construct the boundary ndoe
          Node_pt[global_node_number - 1] =
            finite_element_pt(e)->construct_boundary_node(3);

          // Add to the boundary lookup scheme
          add_boundary_node(bound[global_node_number - 1] - 1,
                            Node_pt[global_node_number - 1]);
        }
        // Otherwise just construct a normal node
        else
        {
          Node_pt[global_node_number - 1] =
            finite_element_pt(e)->construct_node(3);
        }

        done[global_node_number - 1] = true;
        Node_pt[global_node_number - 1]->x(0) = x_node[global_node_number - 1];
        Node_pt[global_node_number - 1]->x(1) = y_node[global_node_number - 1];
        Node_pt[global_node_number - 1]->x(2) = z_node[global_node_number - 1];
      }
      // Otherwise just copy the node numbr accross
      else
      {
        finite_element_pt(e)->node_pt(3) = Node_pt[global_node_number - 1];
      }
      counter++;

      // Loop over the other nodes
      for (unsigned j = 0; j < (n_local_node - 1); j++)
      {
        global_node_number = Global_node[counter];
        if (done[global_node_number - 1] == false)
        // ... -1 because node number begins at 1 in tetgen
        {
          // If we're on a boundary
          if ((boundary_markers_flag == 1) &&
              (bound[global_node_number - 1] > 0))
          {
            // Construct the boundary ndoe
            Node_pt[global_node_number - 1] =
              finite_element_pt(e)->construct_boundary_node(j);

            // Add to the boundary lookup scheme
            add_boundary_node(bound[global_node_number - 1] - 1,
                              Node_pt[global_node_number - 1]);
          }
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
          Node_pt[global_node_number - 1]->x(2) =
            z_node[global_node_number - 1];
        }
        // Otherwise copy the pointer over
        else
        {
          finite_element_pt(e)->node_pt(j) = Node_pt[global_node_number - 1];
        }
        counter++;
      }
    }


    // Resize the "matrix" that stores the boundary id for each
    // face in each element.
    Face_boundary.resize(n_element);
    Face_index.resize(n_element);
    Edge_index.resize(n_element);


    // 0-based index scheme used to construct a global lookup for
    // each face that will be used to uniquely construct interior
    // face nodes.
    // The faces that lie on boundaries will have already been allocated so
    // we can start from there
    Nglobal_face = n_face;

    // Storage for the edges associated with each node.
    // Nodes are indexed using Tetgen's 1-based scheme, which is why
    // there is a +1 here
    Vector<std::set<unsigned>> node_on_edges(n_node + 1);

    // 0-based index scheme used to construct a global lookup for each
    // edge that will be used to uniquely construct interior edge nodes
    Nglobal_edge = 0;

    // Conversion from the edge numbers to the nodes at the end
    // of each each edge
    const unsigned first_local_edge_node[6] = {0, 0, 0, 1, 2, 1};
    const unsigned second_local_edge_node[6] = {1, 2, 3, 2, 3, 3};

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Each element has four faces
      Face_boundary[e].resize(4);
      Face_index[e].resize(4);
      // By default each face is NOT on a boundary
      for (unsigned i = 0; i < 4; i++)
      {
        Face_boundary[e][i] = 0;
      }

      Edge_index[e].resize(6);

      // Read out the global node numbers of the nodes from
      // tetgen's 1-based numbering.
      // The offset is to match the offset used above
      const unsigned element_offset = e * n_local_node;
      unsigned glob_num[4] = {0, 0, 0, 0};
      for (unsigned i = 0; i < 4; ++i)
      {
        glob_num[i] = Global_node[element_offset + ((i + 1) % 4)];
      }

      // Now we know the global node numbers of the elements' four nodes
      // in tetgen's 1-based numbering.

      // Determine whether any of the element faces have already been allocated
      // an index. This may be because they are located on boundaries, or have
      // already been visited The global index for the i-th face will be stored
      // in Face_index[i]

      // Loop over the local faces in the element
      for (unsigned i = 0; i < 4; ++i)
      {
        // On the i-th face, our numbering convention is such that
        // it is the (3-i)th node of the element that is omitted
        const unsigned omitted_node = 3 - i;

        // Start with the set of global faces associated with the i-th node
        // which is always on the i-th face
        std::set<unsigned> input = node_on_faces[glob_num[i]];

        // If there is no data yet, not point doing intersections
        // the face has not been set up
        if (input.size() > 0)
        {
          // Loop over the other nodes on the face
          unsigned local_node = i;
          for (unsigned i2 = 0; i2 < 3; ++i2)
          {
            // Create the next node index (mod 4)
            local_node = (local_node + 1) % 4;
            // Only take the intersection of nodes on the face
            // i.e. not 3-i
            if (local_node != omitted_node)
            {
              // Local storage for the intersection
              std::set<unsigned> intersection_result;
              // Calculate the intersection of the current input
              // and next node
              std::set_intersection(
                input.begin(),
                input.end(),
                node_on_faces[glob_num[local_node]].begin(),
                node_on_faces[glob_num[local_node]].end(),
                std::insert_iterator<std::set<unsigned>>(
                  intersection_result, intersection_result.begin()));
              // Let the result be the next input
              input = intersection_result;
            }
          }
        }

        // If the nodes share more than one global face index, then
        // we have a problem
        if (input.size() > 1)
        {
          throw OomphLibError(
            "Nodes in scaffold mesh share more than one global face",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }

        // If the element's face is not already allocated, the intersection
        // will be empty
        if (input.size() == 0)
        {
          // Allocate the next global index
          Face_index[e][i] = Nglobal_face;
          // Associate the new face index with the nodes
          for (unsigned i2 = 0; i2 < 4; ++i2)
          {
            // Don't add the omitted node
            if (i2 != omitted_node)
            {
              node_on_faces[glob_num[i2]].insert(Nglobal_face);
            }
          }
          // Increment the global face index
          ++Nglobal_face;
        }
        // Otherwise we already have a face
        else if (input.size() == 1)
        {
          const unsigned global_face_index = *input.begin();
          // Set the face index
          Face_index[e][i] = global_face_index;
          // Allocate the boundary index, if it's a boundary
          if (global_face_index < n_face)
          {
            Face_boundary[e][i] = face_boundary[global_face_index];
            // Add the nodes to the boundary look-up scheme in
            // oomph-lib (0-based) index
            for (unsigned i2 = 0; i2 < 4; ++i2)
            {
              // Don't add the omitted node
              if (i2 != omitted_node)
              {
                add_boundary_node(face_boundary[global_face_index] - 1,
                                  Node_pt[glob_num[i2] - 1]);
              }
            }
          }
        }
      } // end of loop over the faces


      // Loop over the element edges and assign global edge numbers
      for (unsigned i = 0; i < 6; ++i)
      {
        // Storage for the result of the intersection
        std::vector<unsigned> local_edge_index;

        const unsigned first_global_num = glob_num[first_local_edge_node[i]];
        const unsigned second_global_num = glob_num[second_local_edge_node[i]];

        // Find the common global edge index for the nodes on
        // the i-th element edge
        std::set_intersection(node_on_edges[first_global_num].begin(),
                              node_on_edges[first_global_num].end(),
                              node_on_edges[second_global_num].begin(),
                              node_on_edges[second_global_num].end(),
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
          node_on_edges[first_global_num].insert(Nglobal_edge);
          node_on_edges[second_global_num].insert(Nglobal_edge);
          // Increment the global edge index
          ++Nglobal_edge;
        }
        // Otherwise we already have an edge
        else if (local_edge_index.size() == 1)
        {
          // Set the edge index from the result of the intersection
          Edge_index[e][i] = local_edge_index[0];
        }
      }

    } // end for e


    // Now determine whether any edges lie on boundaries by using the
    // face boundary scheme and

    // Resize the storage
    Edge_boundary.resize(Nglobal_edge, false);

    // Now loop over all the boundary faces and mark that all edges
    // must also lie on the bounadry
    for (unsigned i = 0; i < n_face; ++i)
    {
      {
        // Storage for the result of the intersection
        std::vector<unsigned> local_edge_index;

        // Find the common global edge index for first edge of the face
        std::set_intersection(node_on_edges[first_node[i]].begin(),
                              node_on_edges[first_node[i]].end(),
                              node_on_edges[second_node[i]].begin(),
                              node_on_edges[second_node[i]].end(),
                              std::insert_iterator<std::vector<unsigned>>(
                                local_edge_index, local_edge_index.begin()));

        // If the nodes do not share exactly one global edge index, then
        // we have a problem
        if (local_edge_index.size() != 1)
        {
          throw OomphLibError(
            "Nodes in scaffold mesh face do not share exactly one global edge",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
        else
        {
          Edge_boundary[local_edge_index[0]] = true;
        }
      }

      {
        // Storage for the result of the intersection
        std::vector<unsigned> local_edge_index;

        // Find the common global edge index for second edge of the face
        std::set_intersection(node_on_edges[second_node[i]].begin(),
                              node_on_edges[second_node[i]].end(),
                              node_on_edges[third_node[i]].begin(),
                              node_on_edges[third_node[i]].end(),
                              std::insert_iterator<std::vector<unsigned>>(
                                local_edge_index, local_edge_index.begin()));

        // If the nodes do not share exactly one global edge index, then
        // we have a problem
        if (local_edge_index.size() != 1)
        {
          throw OomphLibError(
            "Nodes in scaffold mesh face do not share exactly one global edge",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
        else
        {
          Edge_boundary[local_edge_index[0]] = true;
        }
      }

      {
        // Storage for the result of the intersection
        std::vector<unsigned> local_edge_index;

        // Find the common global edge index for third edge of the face
        std::set_intersection(node_on_edges[first_node[i]].begin(),
                              node_on_edges[first_node[i]].end(),
                              node_on_edges[third_node[i]].begin(),
                              node_on_edges[third_node[i]].end(),
                              std::insert_iterator<std::vector<unsigned>>(
                                local_edge_index, local_edge_index.begin()));

        // If the nodes do not share exactly one global edge index, then
        // we have a problem
        if (local_edge_index.size() != 1)
        {
          throw OomphLibError(
            "Nodes in scaffold mesh face do not share exactly one global edge",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
        else
        {
          Edge_boundary[local_edge_index[0]] = true;
        }
      }
    }


  } // end of constructor


  //======================================================================
  /// Constructor: Pass a tetgenio data structure that represents
  /// the input data of the mesh.
  /// The assumptions are that the nodes have been assigned boundary
  /// information which is used in the nodal construction to make sure
  /// that BoundaryNodes are constructed. Any additional boundaries are
  /// determined from the face boundary information.
  //======================================================================
  TetgenScaffoldMesh::TetgenScaffoldMesh(tetgenio& tetgen_data)
  {
    // Find the number of elements
    unsigned n_element = static_cast<unsigned>(tetgen_data.numberoftetrahedra);

    // Read in number of nodes per element
    unsigned n_local_node = static_cast<unsigned>(tetgen_data.numberofcorners);

    // Throw an error if we have anything but linear simplices
    if (n_local_node != 4)
    {
      std::ostringstream error_stream;
      error_stream
        << "Tetgen should only be used to generate 4-noded tetrahedra!\n"
        << "Your tetgen input data, contains data for " << n_local_node
        << "-noded tetrahedra" << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Element attributes may be used to distinguish internal regions
    // NOTE: This stores doubles because tetgen forces us to!
    Element_attribute.resize(n_element, 0.0);

    // Resize storage for the global node numbers listed element-by-element
    Global_node.resize(n_element * n_local_node);

    // Initialise (global) node counter
    unsigned k = 0;
    // Are there attributes?
    unsigned attribute_flag =
      static_cast<unsigned>(tetgen_data.numberoftetrahedronattributes);

    // Read global node numbers for all elements
    if (attribute_flag == 0)
    {
      for (unsigned i = 0; i < n_element; i++)
      {
        for (unsigned j = 0; j < n_local_node; j++)
        {
          Global_node[k] =
            static_cast<unsigned>(tetgen_data.tetrahedronlist[k]);
          k++;
        }
      }
    }
    // Otherwise read in the attributes as well
    else
    {
      for (unsigned i = 0; i < n_element; i++)
      {
        for (unsigned j = 0; j < n_local_node; j++)
        {
          Global_node[k] =
            static_cast<unsigned>(tetgen_data.tetrahedronlist[k]);
          k++;
        }
        Element_attribute[i] = tetgen_data.tetrahedronattributelist[i];
      }
    }

    // Resize the Element vector
    Element_pt.resize(n_element);

    // Read in the number of nodes
    unsigned n_node = tetgen_data.numberofpoints;

    // Create a vector of boolean so as not to create the same node twice
    std::vector<bool> done(n_node, false);

    // Resize the Node vector
    Node_pt.resize(n_node);

    // Flag for boundary markers
    unsigned boundary_markers_flag = 0;
    if (tetgen_data.pointmarkerlist != 0)
    {
      boundary_markers_flag = 1;
    }

    // Create storage for nodal positions and boundary markers
    Vector<double> x_node(n_node);
    Vector<double> y_node(n_node);
    Vector<double> z_node(n_node);
    Vector<unsigned> bound(n_node);

    // We shall ignore all point attributes
    if (boundary_markers_flag == 1)
    {
      for (unsigned i = 0; i < n_node; i++)
      {
        x_node[i] = tetgen_data.pointlist[3 * i];
        y_node[i] = tetgen_data.pointlist[3 * i + 1];
        z_node[i] = tetgen_data.pointlist[3 * i + 2];
        bound[i] = static_cast<unsigned>(tetgen_data.pointmarkerlist[i]);
      }
    }
    // Otherwise no boundary markers
    else
    {
      for (unsigned i = 0; i < n_node; i++)
      {
        x_node[i] = tetgen_data.pointlist[3 * i];
        y_node[i] = tetgen_data.pointlist[3 * i + 1];
        z_node[i] = tetgen_data.pointlist[3 * i + 2];
        bound[i] = 0;
      }
    }


    // Determine highest boundary index
    //------------------------------------
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

    // Now extract face information
    //---------------------------------

    // Number of faces in face file
    unsigned n_face = tetgen_data.numberoftrifaces;

    // Storage for the global node numbers (in the tetgen 1-based
    // numbering scheme!) of the first, second and third  node in
    //  each segment
    Vector<unsigned> first_node(n_face);
    Vector<unsigned> second_node(n_face);
    Vector<unsigned> third_node(n_face);

    // Storage for the boundary marker for each face
    Vector<unsigned> face_boundary(n_face);

    // Storage for the (boundary) faces associated with each node.
    // Nodes are indexed using Tetgen's 1-based scheme, which is why
    // there is a +1 here
    Vector<std::set<unsigned>> node_on_faces(n_node + 1);

    // Extract information for each segment
    for (unsigned i = 0; i < n_face; i++)
    {
      first_node[i] = static_cast<unsigned>(tetgen_data.trifacelist[3 * i]);
      second_node[i] =
        static_cast<unsigned>(tetgen_data.trifacelist[3 * i + 1]);
      third_node[i] = static_cast<unsigned>(tetgen_data.trifacelist[3 * i + 2]);
      face_boundary[i] =
        static_cast<unsigned>(tetgen_data.trifacemarkerlist[i]);
      if (face_boundary[i] > n_bound)
      {
        n_bound = face_boundary[i];
      }
      // Add the face index to each node
      node_on_faces[first_node[i]].insert(i);
      node_on_faces[second_node[i]].insert(i);
      node_on_faces[third_node[i]].insert(i);
    }

    // Extract hole center information
    // unsigned nhole=tetgen_data.numberofholes;
    /*if(nhole!=0)
     {
      Hole_centre.resize(nhole);

      // Coords counter
      unsigned count_coords=0;

      // Loop over the holes to get centre coords
      for(unsigned ihole=0;ihole<nhole;ihole++)
       {
        Hole_centre[ihole].resize(3);

        // Read the centre value
        Hole_centre[ihole][0]=triangle_data.holelist[count_coords];
        Hole_centre[ihole][1]=triangle_data.holelist[count_coords+1];
        Hole_centre[ihole][2]=triangle_data.holelist[count_coords+2];

        // Increment coords counter
        count_coords+=3;
       }
     }
    else
     {
      Hole_centre.resize(0);
      }*/


    // Set number of boundaries
    if (n_bound > 0)
    {
      this->set_nboundary(n_bound);
    }

    // Create the elements
    unsigned counter = 0;
    for (unsigned e = 0; e < n_element; e++)
    {
      Element_pt[e] = new TElement<3, 2>;
      unsigned global_node_number = Global_node[counter];
      if (done[global_node_number - 1] == false)
      // ... -1 because node number begins at 1 in tetgen
      {
        // If the node is on a boundary, construct a boundary node
        if ((boundary_markers_flag == 1) && (bound[global_node_number - 1] > 0))
        {
          // Construct the boundary ndoe
          Node_pt[global_node_number - 1] =
            finite_element_pt(e)->construct_boundary_node(3);

          // Add to the boundary lookup scheme
          add_boundary_node(bound[global_node_number - 1] - 1,
                            Node_pt[global_node_number - 1]);
        }
        // Otherwise just construct a normal node
        else
        {
          Node_pt[global_node_number - 1] =
            finite_element_pt(e)->construct_node(3);
        }

        done[global_node_number - 1] = true;
        Node_pt[global_node_number - 1]->x(0) = x_node[global_node_number - 1];
        Node_pt[global_node_number - 1]->x(1) = y_node[global_node_number - 1];
        Node_pt[global_node_number - 1]->x(2) = z_node[global_node_number - 1];
      }
      // Otherwise just copy the node numbr accross
      else
      {
        finite_element_pt(e)->node_pt(3) = Node_pt[global_node_number - 1];
      }
      counter++;

      // Loop over the other nodes
      for (unsigned j = 0; j < (n_local_node - 1); j++)
      {
        global_node_number = Global_node[counter];
        if (done[global_node_number - 1] == false)
        // ... -1 because node number begins at 1 in tetgen
        {
          // If we're on a boundary
          if ((boundary_markers_flag == 1) &&
              (bound[global_node_number - 1] > 0))
          {
            // Construct the boundary ndoe
            Node_pt[global_node_number - 1] =
              finite_element_pt(e)->construct_boundary_node(j);

            // Add to the boundary lookup scheme
            add_boundary_node(bound[global_node_number - 1] - 1,
                              Node_pt[global_node_number - 1]);
          }
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
          Node_pt[global_node_number - 1]->x(2) =
            z_node[global_node_number - 1];
        }
        // Otherwise copy the pointer over
        else
        {
          finite_element_pt(e)->node_pt(j) = Node_pt[global_node_number - 1];
        }
        counter++;
      }
    }


    // Resize the "matrix" that stores the boundary id for each
    // face in each element.
    Face_boundary.resize(n_element);
    Face_index.resize(n_element);
    Edge_index.resize(n_element);


    // 0-based index scheme used to construct a global lookup for
    // each face that will be used to uniquely construct interior
    // face nodes.
    // The faces that lie on boundaries will have already been allocated so
    // we can start from there
    Nglobal_face = n_face;

    // Storage for the edges associated with each node.
    // Nodes are indexed using Tetgen's 1-based scheme, which is why
    // there is a +1 here
    Vector<std::set<unsigned>> node_on_edges(n_node + 1);

    // 0-based index scheme used to construct a global lookup for each
    // edge that will be used to uniquely construct interior edge nodes
    Nglobal_edge = 0;

    // Conversion from the edge numbers to the nodes at the end
    // of each each edge
    const unsigned first_local_edge_node[6] = {0, 0, 0, 1, 2, 1};
    const unsigned second_local_edge_node[6] = {1, 2, 3, 2, 3, 3};

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Each element has four faces
      Face_boundary[e].resize(4);
      Face_index[e].resize(4);
      // By default each face is NOT on a boundary
      for (unsigned i = 0; i < 4; i++)
      {
        Face_boundary[e][i] = 0;
      }

      Edge_index[e].resize(6);

      // Read out the global node numbers of the nodes from
      // tetgen's 1-based numbering.
      // The offset is to match the offset used above
      const unsigned element_offset = e * n_local_node;
      unsigned glob_num[4] = {0, 0, 0, 0};
      for (unsigned i = 0; i < 4; ++i)
      {
        glob_num[i] = Global_node[element_offset + ((i + 1) % 4)];
      }

      // Now we know the global node numbers of the elements' four nodes
      // in tetgen's 1-based numbering.

      // Determine whether any of the element faces have already been allocated
      // an index. This may be because they are located on boundaries, or have
      // already been visited The global index for the i-th face will be stored
      // in Face_index[i]

      // Loop over the local faces in the element
      for (unsigned i = 0; i < 4; ++i)
      {
        // On the i-th face, our numbering convention is such that
        // it is the (3-i)th node of the element that is omitted
        const unsigned omitted_node = 3 - i;

        // Start with the set of global faces associated with the i-th node
        // which is always on the i-th face
        std::set<unsigned> input = node_on_faces[glob_num[i]];

        // If there is no data yet, not point doing intersections
        // the face has not been set up
        if (input.size() > 0)
        {
          // Loop over the other nodes on the face
          unsigned local_node = i;
          for (unsigned i2 = 0; i2 < 3; ++i2)
          {
            // Create the next node index (mod 4)
            local_node = (local_node + 1) % 4;
            // Only take the intersection of nodes on the face
            // i.e. not 3-i
            if (local_node != omitted_node)
            {
              // Local storage for the intersection
              std::set<unsigned> intersection_result;
              // Calculate the intersection of the current input
              // and next node
              std::set_intersection(
                input.begin(),
                input.end(),
                node_on_faces[glob_num[local_node]].begin(),
                node_on_faces[glob_num[local_node]].end(),
                std::insert_iterator<std::set<unsigned>>(
                  intersection_result, intersection_result.begin()));
              // Let the result be the next input
              input = intersection_result;
            }
          }
        }

        // If the nodes share more than one global face index, then
        // we have a problem
        if (input.size() > 1)
        {
          throw OomphLibError(
            "Nodes in scaffold mesh share more than one global face",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }

        // If the element's face is not already allocated, the intersection
        // will be empty
        if (input.size() == 0)
        {
          // Allocate the next global index
          Face_index[e][i] = Nglobal_face;
          // Associate the new face index with the nodes
          for (unsigned i2 = 0; i2 < 4; ++i2)
          {
            // Don't add the omitted node
            if (i2 != omitted_node)
            {
              node_on_faces[glob_num[i2]].insert(Nglobal_face);
            }
          }
          // Increment the global face index
          ++Nglobal_face;
        }
        // Otherwise we already have a face
        else if (input.size() == 1)
        {
          const unsigned global_face_index = *input.begin();
          // Set the face index
          Face_index[e][i] = global_face_index;
          // Allocate the boundary index, if it's a boundary
          if (global_face_index < n_face)
          {
            Face_boundary[e][i] = face_boundary[global_face_index];
            // Add the nodes to the boundary look-up scheme in
            // oomph-lib (0-based) index
            for (unsigned i2 = 0; i2 < 4; ++i2)
            {
              // Don't add the omitted node
              if (i2 != omitted_node)
              {
                add_boundary_node(face_boundary[global_face_index] - 1,
                                  Node_pt[glob_num[i2] - 1]);
              }
            }
          }
        }
      } // end of loop over the faces


      // Loop over the element edges and assign global edge numbers
      for (unsigned i = 0; i < 6; ++i)
      {
        // Storage for the result of the intersection
        std::vector<unsigned> local_edge_index;

        const unsigned first_global_num = glob_num[first_local_edge_node[i]];
        const unsigned second_global_num = glob_num[second_local_edge_node[i]];

        // Find the common global edge index for the nodes on
        // the i-th element edge
        std::set_intersection(node_on_edges[first_global_num].begin(),
                              node_on_edges[first_global_num].end(),
                              node_on_edges[second_global_num].begin(),
                              node_on_edges[second_global_num].end(),
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
          node_on_edges[first_global_num].insert(Nglobal_edge);
          node_on_edges[second_global_num].insert(Nglobal_edge);
          // Increment the global edge index
          ++Nglobal_edge;
        }
        // Otherwise we already have an edge
        else if (local_edge_index.size() == 1)
        {
          // Set the edge index from the result of the intersection
          Edge_index[e][i] = local_edge_index[0];
        }
      }

    } // end for e


    // Now determine whether any edges lie on boundaries by using the
    // face boundary scheme and

    // Resize the storage
    Edge_boundary.resize(Nglobal_edge, false);

    // Now loop over all the boundary faces and mark that all edges
    // must also lie on the bounadry
    for (unsigned i = 0; i < n_face; ++i)
    {
      {
        // Storage for the result of the intersection
        std::vector<unsigned> local_edge_index;

        // Find the common global edge index for first edge of the face
        std::set_intersection(node_on_edges[first_node[i]].begin(),
                              node_on_edges[first_node[i]].end(),
                              node_on_edges[second_node[i]].begin(),
                              node_on_edges[second_node[i]].end(),
                              std::insert_iterator<std::vector<unsigned>>(
                                local_edge_index, local_edge_index.begin()));

        // If the nodes do not share exactly one global edge index, then
        // we have a problem
        if (local_edge_index.size() != 1)
        {
          throw OomphLibError(
            "Nodes in scaffold mesh face do not share exactly one global edge",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
        else
        {
          Edge_boundary[local_edge_index[0]] = true;
        }
      }

      {
        // Storage for the result of the intersection
        std::vector<unsigned> local_edge_index;

        // Find the common global edge index for second edge of the face
        std::set_intersection(node_on_edges[second_node[i]].begin(),
                              node_on_edges[second_node[i]].end(),
                              node_on_edges[third_node[i]].begin(),
                              node_on_edges[third_node[i]].end(),
                              std::insert_iterator<std::vector<unsigned>>(
                                local_edge_index, local_edge_index.begin()));

        // If the nodes do not share exactly one global edge index, then
        // we have a problem
        if (local_edge_index.size() != 1)
        {
          throw OomphLibError(
            "Nodes in scaffold mesh face do not share exactly one global edge",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
        else
        {
          Edge_boundary[local_edge_index[0]] = true;
        }
      }

      {
        // Storage for the result of the intersection
        std::vector<unsigned> local_edge_index;

        // Find the common global edge index for third edge of the face
        std::set_intersection(node_on_edges[first_node[i]].begin(),
                              node_on_edges[first_node[i]].end(),
                              node_on_edges[third_node[i]].begin(),
                              node_on_edges[third_node[i]].end(),
                              std::insert_iterator<std::vector<unsigned>>(
                                local_edge_index, local_edge_index.begin()));

        // If the nodes do not share exactly one global edge index, then
        // we have a problem
        if (local_edge_index.size() != 1)
        {
          throw OomphLibError(
            "Nodes in scaffold mesh face do not share exactly one global edge",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
        else
        {
          Edge_boundary[local_edge_index[0]] = true;
        }
      }
    }


  } // end of constructor


} // namespace oomph
