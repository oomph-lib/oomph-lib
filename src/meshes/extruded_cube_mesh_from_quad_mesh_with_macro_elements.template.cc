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
// Oomph-lib headers
#include "extruded_cube_mesh_from_quad_mesh_with_macro_elements.template.h"

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

namespace oomph
{
  //========================================================================
  /// \short Get all the boundary information of an element using the
  /// input (quad_mesh_pt) mesh. If the element lies on a boundary then
  /// the user will be given the corresponding boundary index and the
  /// index of the face of quad_el_pt attached to the boundary. If the
  /// element does NOT lie on any boundaries, this function simply
  /// returns a vector of size zero.
  //========================================================================
  template<class ELEMENT>
  Vector<std::pair<unsigned, int>> ExtrudedCubeMeshFromQuadMesh<
    ELEMENT>::get_element_boundary_information(QuadMeshBase* quad_mesh_pt,
                                               FiniteElement* quad_el_pt)
  {
    // Allocate space for the boundary info (initialised to empty)
    Vector<std::pair<unsigned, int>> boundary_info;

    // Find the number of boundaries in the input mesh
    unsigned n_boundary = quad_mesh_pt->nboundary();

    // Loop over the boundaries of the mesh
    for (unsigned b = 0; b < n_boundary; b++)
    {
      // How many elements are there on the b-th boundary?
      unsigned n_boundary_element = quad_mesh_pt->nboundary_element(b);

      // Loop over the elements on the b-th boundary
      for (unsigned e = 0; e < n_boundary_element; e++)
      {
        // Is the e-th element on the b-th boundary the input element?
        if (quad_el_pt == quad_mesh_pt->boundary_element_pt(b, e))
        {
          // Create a pair to hold the boundary index and element face index
          std::pair<unsigned, int> boundary_and_face_index;

          // Set the first entry (boundary index)
          boundary_and_face_index.first = b;

          // Set the second entry (face index)
          boundary_and_face_index.second =
            quad_mesh_pt->face_index_at_boundary(b, e);

          // Add the boundary index to the boundary_info vector
          boundary_info.push_back(boundary_and_face_index);
        }
      } // for (unsigned e=0;e<n_boundary_element;e++)
    } // for (unsigned b=0;b<n_boundary;b++)

    // Return the result
    return boundary_info;
  } // End of get_element_boundary_information

  //========================================================================
  /// Generic mesh construction. This function contains all the details of
  /// the mesh generation process, including all the tedious loops, counting
  /// spacing and boundary functions.
  /// NOTE: The boundary number of the extruded mesh will follow the same
  /// numbering as the input quad mesh. The newly created "front" and "back"
  /// face of the 3D mesh will added after those boundaries. For example,
  /// if the input mesh has 4 boundaries; b_0, b_1, b_2 & b_3 then the
  /// extruded mesh will have 6 boundaries; eb_0, eb_1, eb_2, eb_3, eb_4 &
  /// eb_5 where eb_4 is the "front" face and eb5 is the "back" face. The
  /// boundaries eb_i here satisfy eb_i = (b_i) x [Zmin,Zmax].
  //========================================================================
  template<class ELEMENT>
  void ExtrudedCubeMeshFromQuadMesh<ELEMENT>::build_mesh(
    QuadMeshBase* quad_mesh_pt, TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with 3D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

    // Need at least 2 elements in the z-direction
    if (Nz == 0)
    {
      // Create an ostringstream object to create an error message
      std::ostringstream error_message;

      // Create the error message
      error_message << "Extrude2DQuadMeshTo3DCubeMesh needs at least two "
                    << "elements in each,\ncoordinate direction. You have "
                    << "specified Nz=" << Nz << std::endl;

      // Throw an error
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Get the current (i.e. start) time
    double t_mesh_setup_start = TimingHelpers::timer();

    // Get the current (i.e. start) time
    double t_mesh_construct_start = TimingHelpers::timer();

    // Get the number of boundaries in the input mesh
    unsigned n_boundary_in_quad_mesh = quad_mesh_pt->nboundary();

    // The number of boundaries in the extruded mesh
    unsigned n_boundary_in_extruded_mesh = n_boundary_in_quad_mesh + 2;

    // Set the number of boundaries in the extruded mesh (just need
    // to add the front and back face to the number of boundaries)
    set_nboundary(n_boundary_in_extruded_mesh);

    // Get the number of elements in the input mesh
    unsigned n_element_in_quad_mesh = quad_mesh_pt->nelement();

    // Allocate storage for the boundary information of each quad element.
    // Takes up more memory but avoids recomputing it many, many times
    Vector<Vector<std::pair<unsigned, int>>> boundary_information(
      n_element_in_quad_mesh);

    // Loop over all the elements in the quad mesh
    for (unsigned e = 0; e < n_element_in_quad_mesh; e++)
    {
      // Get a pointer to the j-th element in the mesh
      FiniteElement* quad_el_pt = quad_mesh_pt->finite_element_pt(e);

      // Get the boundary information
      boundary_information[e] =
        get_element_boundary_information(quad_mesh_pt, quad_el_pt);
    }

    // Allocate storage for the number of elements in the extruded mesh
    Element_pt.resize(n_element_in_quad_mesh * Nz);

    // Create the first element
    ELEMENT* temp_el_pt = new ELEMENT;

    // Read out the number of linear points in the element (in one direction)
    unsigned n_node_1d = temp_el_pt->nnode_1d();

    // Store the variable n_node_1d as the private variable, Np
    N_node_1d = n_node_1d;

    // Delete the element
    delete temp_el_pt;

    // Make it a null pointer
    temp_el_pt = 0;

    // Need the same number of nodes in one direction in the 2D and 3D element
    if ((n_node_1d != quad_mesh_pt->finite_element_pt(0)->nnode_1d()) ||
        (n_node_1d * n_node_1d != quad_mesh_pt->finite_element_pt(0)->nnode()))
    {
      // Create an ostringstream object to create an error message
      std::ostringstream error_message;

      // Create the error message
      error_message
        << "Extrude2DQuadMeshTo3DCubeMesh needs the number of "
        << "nodes (in the 2D element) in one direction\n to match "
        << "the number of nodes in one direction in the 3D element!";

      // Throw an error
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Get the total number of nodes in the input mesh
    unsigned long n_node_in_quad_mesh = quad_mesh_pt->nnode();

    // Loop over the nodes in the spatial mesh
    for (unsigned i = 0; i < n_node_in_quad_mesh; i++)
    {
      // Check to see if the i-th node is hanging
      if (quad_mesh_pt->node_pt(i)->is_hanging())
      {
        // Create an output stream
        std::ostringstream warning_message_stream;

        // Create an error message
        warning_message_stream << "Extrusion machinery is not currently able "
                               << "to deal with hanging nodes!" << std::endl;

        // Throw an error
        throw OomphLibError(warning_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    } // for (unsigned i=0;i<n_node_in_quad_mesh;i++)

    // Number of nodes in the z-direction
    unsigned n_node_in_z_direction = (1 + (n_node_1d - 1) * Nz);

    // Get the total number of nodes in the input mesh
    unsigned long n_node_in_cube_mesh =
      n_node_in_quad_mesh * n_node_in_z_direction;

    // Allocate storage for the number of nodes in the extruded mesh
    Node_pt.resize(n_node_in_cube_mesh);

    // Counter for the number of nodes in the extruded mesh
    unsigned long node_count = 0;

    // Loop over each element slice in the z-direction
    for (unsigned i = 0; i < Nz; i++)
    {
      // Need a map at each element slice to tell us whether or not we've
      // already created an edge (and thus the nodes on that edge). If the
      // nodes have already been created then all we have to do is copy them
      // over in reverse order
      std::map<std::pair<Edge, double>, Vector<Node*>> edge_to_nodes_map;

      // As several edges can attach to a single (corner) node we can
      // accidentally create duplicate nodes (i.e. nodes that lie at
      // the same spatial position) because the edges aren't covered in
      // the required order (something we shouldn't rely on anyway!). To
      // get around this we store the corner nodes of each quad element
      // with its corresponding z-value. This will in turn return the
      // corresponding node in the extruded mesh
      std::map<std::pair<Node*, double>, Node*>
        quad_corner_node_to_cube_node_map;

      // Loop over all the elements in one slice
      for (unsigned e = 0; e < n_element_in_quad_mesh; e++)
      {
        // Get a pointer to the e-th element in the input mesh
        FiniteElement* quad_el_pt = quad_mesh_pt->finite_element_pt(e);

        // Get the index of the cube element in the extruded mesh
        unsigned element_index = i * n_element_in_quad_mesh + e;

        // Create the e-th element in the i-th slice of the extruded mesh
        ELEMENT* cube_el_pt = new ELEMENT;

        // Put it in global storage!
        Element_pt[element_index] = cube_el_pt;

        // Index of the local node number (in the quad element)
        unsigned quad_local_node_number = 0;

        // Index of the local node number (in the cube element)
        unsigned cube_local_node_number = 0;

        // Cache the boundary information associated with this element
        Vector<std::pair<unsigned, int>> el_boundary_information =
          boundary_information[e];

        // Loop over the nodes in the third local coordinate direction
        for (unsigned n_2 = 0; n_2 < n_node_1d; n_2++)
        {
          // Calculate the z-value at the n_2-th (z-)node in the i-th element
          // slice
          double z_value = z_spacing_function(i, n_2);

          //---------------------------------------------------------------
          // At a given node slice, there are three options:
          //     (1) If i>0 and n_2=0 then we already created the nodes
          //         (in the current slice) when we created the nodes in
          //         the previous slice;
          //     (2) If (1) doesn't hold we may have still have created
          //         certain nodes along a face/edge of an element that
          //         we previously created (in this node slice);
          //     (3) No nodes in this node slice have been created yet so
          //         set them all up!
          // NOTE: if (1) occurs, we can skip to the next node slice but
          // we can't if (2) occurs (or (3), obviously!).
          //---------------------------------------------------------------
          //---------------------------------------------------------------
          // Case (1): If we're past the first element slice and we're on
          // the first node slice (in the element) then copy the nodes from
          // the appropriate element.
          //---------------------------------------------------------------
          if ((i > 0) && (n_2 == 0))
          {
            // Get the index of the cube element (to copy) in the extruded mesh
            unsigned copy_element_index = (i - 1) * n_element_in_quad_mesh + e;

            // Get a pointer to the element in the previous slice
            // containing the node we want
            FiniteElement* copy_cube_el_pt =
              finite_element_pt(copy_element_index);

            // Storage for the local node number in the element that we're
            // copying
            unsigned copy_cube_local_node_number = 0;

            // Loop over the nodes in the second local coordinate direction
            for (unsigned n_1 = 0; n_1 < n_node_1d; n_1++)
            {
              // Loop over the nodes in the first local coordinate direction
              for (unsigned n_0 = 0; n_0 < n_node_1d; n_0++)
              {
                // Calculate the local node number (in the extruded mesh)
                cube_local_node_number = n_0 + (n_node_1d * n_1);

                // Calculate the local node number (in the input mesh)
                quad_local_node_number = n_0 + (n_node_1d * n_1);

                // Calculate the local node number (in the element we're
                // copying it from)
                copy_cube_local_node_number =
                  cube_local_node_number +
                  n_node_1d * n_node_1d * (n_node_1d - 1);

                // Copy the node over
                cube_el_pt->node_pt(cube_local_node_number) =
                  copy_cube_el_pt->node_pt(copy_cube_local_node_number);
              } // for (unsigned n_0=0;n_0<n_node_1d;n_0++)
            } // for (unsigned n_1=0;n_1<n_node_1d;n_1++)
          }
          //---------------------------------------------------------------
          // Case (2) & (3): Loop over the edges in the current node slice
          // and check if any of them have already been set up (and thus
          // stored in edge_to_nodes_map). If they have then copy them over
          // and construct all the other nodes. If they haven't then simply
          // construct all of the nodes.
          // NOTE: if an edge has been created already, we can't make the
          // nodes first, assign the coordinates THEN check as it would
          // be a wasteful operation. Instead, we get the nodes describing
          // an edge from the input mesh (quad_mesh_pt) and pair it with
          // the current z-value.
          //---------------------------------------------------------------
          else
          {
            // Create a vector of bools to tell us if we've set the nodes up
            // in the current node slice
            std::vector<bool> has_node_been_setup(n_node_1d * n_node_1d, false);

            // Number of edges in a slice (N/E/S/W)
            unsigned n_edge = 4;

            // Storage for the edges
            Vector<std::pair<Edge, double>> edges_and_z_value;

            // Vector to tell us if we are copying any edges (to make sure
            // we don't store the edge again later!)
            std::vector<bool> has_edge_been_done(n_edge, false);

            // What is the last edge? Convert this to an int so the compiler
            // doesn't complain about comparisons between unsigned and signed
            // integers...
            int last_edge = int(QuadTreeNames::N + n_edge);

            // Loop over the edges
            for (int i_edge = QuadTreeNames::N; i_edge < last_edge; i_edge++)
            {
              // Pointer for the first node associated with this edge
              Node* node1_pt = 0;

              // Pointer for the second node associated with this edge
              Node* node2_pt = 0;

              // Find the corner nodes associated with this edge
              switch (i_edge)
              {
                case QuadTreeNames::N:
                  // The first node
                  node1_pt = quad_el_pt->node_pt((n_node_1d * n_node_1d) - 1);

                  // The second node
                  node2_pt = quad_el_pt->node_pt(n_node_1d * (n_node_1d - 1));

                  // Break
                  break;
                case QuadTreeNames::E:
                  // The first node
                  node1_pt = quad_el_pt->node_pt(n_node_1d - 1);

                  // The second node
                  node2_pt = quad_el_pt->node_pt((n_node_1d * n_node_1d) - 1);

                  // Break
                  break;
                case QuadTreeNames::S:
                  // The first node
                  node1_pt = quad_el_pt->node_pt(0);

                  // The second node
                  node2_pt = quad_el_pt->node_pt(n_node_1d - 1);

                  // Break
                  break;
                case QuadTreeNames::W:
                  // The first node
                  node1_pt = quad_el_pt->node_pt(n_node_1d * (n_node_1d - 1));

                  // The second node
                  node2_pt = quad_el_pt->node_pt(0);

                  // Break
                  break;
                default:
                  // Create an ostringstream object to create an error message
                  std::ostringstream error_message_stream;

                  // Create the error message
                  error_message_stream
                    << "Input is " << i_edge << " but it can only\n"
                    << "be either N/E/S/W, i.e. " << QuadTreeNames::N << ","
                    << QuadTreeNames::E << "," << QuadTreeNames::S << " or "
                    << QuadTreeNames::W << std::endl;

                  // Throw an error
                  throw OomphLibError(error_message_stream.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
              }

              // Create the edge associated with the i_edge-th edge of this
              // element
              Edge edge(node1_pt, node2_pt);

              // Pair it up with the z-value at the current node slice
              std::pair<Edge, double> edge_and_z_value_pair(edge, z_value);

              // Store the edge (with its corresponding z-value)
              edges_and_z_value.push_back(edge_and_z_value_pair);

              // Is there is a matching entry in the map?
              std::map<std::pair<Edge, double>, Vector<Node*>>::iterator it =
                edge_to_nodes_map.find(edge_and_z_value_pair);

              // If we're not at the end then it has already been created
              if (it != edge_to_nodes_map.end())
              {
                // Get the vector of node pointers
                Vector<Node*> edge_nodes_pt = it->second;

#ifdef PARANOID
                // Sanity check; make sure enough nodes have been stored
                if (edge_nodes_pt.size() != n_node_1d)
                {
                  // Throw an error
                  throw OomphLibError("Not enough nodes in edge_nodes_pt!",
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
#endif

                // Remember that the i_edge-th edge has already been created!
                switch (i_edge)
                {
                  case QuadTreeNames::N:
                    // Indicate that the Northern edge has been done
                    has_edge_been_done[0] = true;

                    // Break
                    break;
                  case QuadTreeNames::E:
                    // Indicate that the Eastern edge has been done
                    has_edge_been_done[1] = true;

                    // Break
                    break;
                  case QuadTreeNames::S:
                    // Indicate that the Southern edge has been done
                    has_edge_been_done[2] = true;

                    // Break
                    break;
                  case QuadTreeNames::W:
                    // Indicate that the Western edge has been done
                    has_edge_been_done[3] = true;

                    // Break
                    break;
                  default:
                    // Throw an error
                    throw OomphLibError("Incorrect edge type!",
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                }

                //---------------------------------------------------------------
                // Assign the nodes appropriately depending on which edge has
                // already been created.
                // NOTE: the nodes will have been created in anti-clockwise
                // order which means we will have to assign them in clockwise
                // order (from the perspective of the current element).
                //---------------------------------------------------------------
                // Loop over the nodes
                for (unsigned i_node = 0; i_node < n_node_1d; i_node++)
                {
                  // Calculate the value of quad_local_node_number depending on
                  // which edge the node lies on
                  switch (i_edge)
                  {
                    case QuadTreeNames::N:
                      // Calculate this Node's quad local node number
                      quad_local_node_number =
                        ((n_node_1d * n_node_1d) - 1) - i_node;

                      // Break
                      break;
                    case QuadTreeNames::E:
                      // Calculate this Node's quad local node number
                      quad_local_node_number =
                        (n_node_1d - 1) + (i_node * n_node_1d);

                      // Break
                      break;
                    case QuadTreeNames::S:
                      // Calculate this Node's quad local node number
                      quad_local_node_number = i_node;

                      // Break
                      break;
                    case QuadTreeNames::W:
                      // Calculate this Node's quad local node number
                      quad_local_node_number =
                        (n_node_1d * (n_node_1d - 1)) - (i_node * n_node_1d);

                      // Break
                      break;
                  }

                  // Calculate this Node's cube local node number
                  cube_local_node_number =
                    quad_local_node_number + (n_node_1d * n_node_1d * n_2);

                  // Assign the i_node-th node
                  cube_el_pt->node_pt(cube_local_node_number) =
                    edge_nodes_pt[(n_node_1d - 1) - i_node];

                  // Indicate that the appropriate node has been set up
                  has_node_been_setup[quad_local_node_number] = true;
                } // for (unsigned i_node=0;i_node<n_node_1d;i_node++)
              } // if (it!=edge_to_nodes_map.end())
            } // for (unsigned i_edge=0;i<n_edge;i_edge++)

            // The number of corners
            unsigned n_corner = 4;

            // Loop over the corner nodes
            for (unsigned i_corner = 0; i_corner < n_corner; i_corner++)
            {
              // Calculate the value of quad_local_node_number depending on
              // which corner the node lies on
              switch (i_corner)
              {
                case 0:
                  // Bottom-left corner
                  quad_local_node_number = 0;

                  // Break
                  break;
                case 1:
                  // Bottom-right corner
                  quad_local_node_number = n_node_1d - 1;

                  // Break
                  break;
                case 2:
                  // Top-left corner
                  quad_local_node_number = n_node_1d * (n_node_1d - 1);

                  // Break
                  break;
                case 3:
                  // Top-right corner
                  quad_local_node_number = (n_node_1d * n_node_1d) - 1;

                  // Break
                  break;
              }

              // If the node hasn't already been copied from some other edge
              if (!has_node_been_setup[quad_local_node_number])
              {
                // Get a pointer to the first corner node (in the quad mesh)
                Node* quad_corner_node_pt =
                  quad_el_pt->node_pt(quad_local_node_number);

                // Pair up the corner node and z-value
                std::pair<Node*, double> quad_corner_node_and_z_value(
                  quad_corner_node_pt, z_value);

                // Have we made this Node yet?
                std::map<std::pair<Node*, double>, Node*>::iterator it =
                  quad_corner_node_to_cube_node_map.find(
                    quad_corner_node_and_z_value);

                // If we found a corresponding entry
                if (it != quad_corner_node_to_cube_node_map.end())
                {
                  // Index of the local node number
                  cube_local_node_number =
                    quad_local_node_number + (n_node_1d * n_node_1d * n_2);

                  // Get the node and store it in the element
                  cube_el_pt->node_pt(cube_local_node_number) = it->second;

                  // Indicate that the node has now been created
                  has_node_been_setup[quad_local_node_number] = true;
                }
              } // if (!has_node_been_setup[quad_local_node_number])
            } // for (unsigned i_corner=0;i_corner<n_corner;i_corner++)

            //------------------------------------------------------------------
            // By this point, all nodes which have already been created have
            // been copied over. All the other nodes now need to be constructed.
            //------------------------------------------------------------------
            // Loop over the nodes in the second local coordinate direction
            for (unsigned n_1 = 0; n_1 < n_node_1d; n_1++)
            {
              // Loop over the nodes in the first local coordinate direction
              for (unsigned n_0 = 0; n_0 < n_node_1d; n_0++)
              {
                // Calculate the local node number (in the input mesh)
                quad_local_node_number = n_0 + (n_node_1d * n_1);

                // Calculate the local node number (in the extruded mesh)
                cube_local_node_number =
                  quad_local_node_number + (n_node_1d * n_node_1d * n_2);

                // Check if the node has already been created (i.e. by copying)
                if (has_node_been_setup[quad_local_node_number])
                {
                  // If this node has already been set up then skip it
                  continue;
                }

                // Pointer to point to the node we are about to create
                Node* cube_node_pt = 0;

                //--------------------------------------------------------------
                //// Check if the node that we're about to construct lies on the
                // Front or back face of the mesh. The node can only lie on the
                // front face if we're on the first element slice and first node
                // slice. In contrast, the node can only lie on the back face if
                // we're on the final element slice and final node slice.
                //--------------------------------------------------------------
                // The node lies on the front face or the back face
                if (((i == 0) && (n_2 == 0)) ||
                    ((i == Nz - 1) && (n_2 == n_node_1d - 1)))
                {
                  // The node lies on the front face
                  if ((i == 0) && (n_2 == 0))
                  {
                    // Create the boundary node
                    cube_node_pt = cube_el_pt->construct_boundary_node(
                      cube_local_node_number, time_stepper_pt);

                    // Indicate that the node has been created now
                    has_node_been_setup[quad_local_node_number] = true;

                    // Add the node to the appropriate boundary
                    add_boundary_node(n_boundary_in_quad_mesh, cube_node_pt);
                  }
                  // The node lies on the back face
                  else if ((i == Nz - 1) && (n_2 == n_node_1d - 1))
                  {
                    // Create the boundary node
                    cube_node_pt = cube_el_pt->construct_boundary_node(
                      cube_local_node_number, time_stepper_pt);

                    // Indicate that the node has been created now
                    has_node_been_setup[quad_local_node_number] = true;

                    // Add the node to the appropriate boundary
                    add_boundary_node(n_boundary_in_quad_mesh + 1,
                                      cube_node_pt);
                  } // if ((i==0)&&(n_2==0))

                  // Get a pointer to the matching quad mesh node
                  Node* quad_node_pt =
                    quad_el_pt->node_pt(quad_local_node_number);

                  // If this Node lies on a boundary
                  if (quad_node_pt->is_on_boundary())
                  {
                    // Storage for the set of boundaries that this Node lies on
                    std::set<unsigned>* boundaries_pt;

                    // Get the boundaries that this Node lies on
                    quad_node_pt->get_boundaries_pt(boundaries_pt);

                    // Create an iterator to loop over the entries of the set
                    std::set<unsigned>::const_iterator const_it;

                    // Iterate over the boundaries
                    for (const_it = boundaries_pt->begin();
                         const_it != boundaries_pt->end();
                         const_it++)
                    {
                      // Add the node to the appropriate boundary (if it hasn't
                      // already been added)
                      add_boundary_node(*const_it, cube_node_pt);
                    }
                  } // quad_node_pt->is_on_boundary()
                } // if (((i==0)&&(n_2==0))||((i==Nz-1)&&(n_2==n_node_1d-1)))


                // Check if there is any other boundary information
                if (el_boundary_information.size() > 0)
                {
                  // How many boundaries are there to cover?
                  unsigned n_boundaries = el_boundary_information.size();

                  // Loop over the boundaries
                  for (unsigned i_bound = 0; i_bound < n_boundaries; i_bound++)
                  {
                    // If the node still hasn't been created yet
                    if (!has_node_been_setup[quad_local_node_number])
                    {
                      // Create the boundary node
                      cube_node_pt = cube_el_pt->construct_boundary_node(
                        cube_local_node_number, time_stepper_pt);

                      // Indicate that the node has been created now
                      has_node_been_setup[quad_local_node_number] = true;
                    }

                    // Get the boundary information associated with the
                    // i_bound-th entry
                    std::pair<unsigned, int> boundary_and_face_index =
                      el_boundary_information[i_bound];

                    // Decide whether or not the node lies on the boundary. If
                    // it does then add it to the boundary
                    switch (boundary_and_face_index.second)
                    {
                      case -1:
                        // If the s[0]=-1 boundary of the element lies on the
                        // boundary the node can only lie on the boundary if
                        // n_0=0.
                        if (n_0 == 0)
                        {
                          // Now add the node to the appropriate boundary
                          add_boundary_node(boundary_and_face_index.first,
                                            cube_node_pt);
                        }

                        // Break
                        break;
                      case 1:
                        // If the s[0]=1 boundary of the element lies on the
                        // boundary the node can only lie on the boundary if
                        // n_0=n_node_1d-1
                        if (n_0 == n_node_1d - 1)
                        {
                          // Now add the node to the appropriate boundary
                          add_boundary_node(boundary_and_face_index.first,
                                            cube_node_pt);
                        }

                        // Break
                        break;
                      case -2:
                        // If the s[1]=-1 boundary of the element lies on the
                        // boundary the node can only lie on the boundary if
                        // n_1=0.
                        if (n_1 == 0)
                        {
                          // Now add the node to the appropriate boundary
                          add_boundary_node(boundary_and_face_index.first,
                                            cube_node_pt);
                        }

                        // Break
                        break;
                      case 2:
                        // If the s[1]=1 boundary of the element lies on the
                        // boundary the node can only lie on the boundary if
                        // n_1=n_node_1d-1.
                        if (n_1 == n_node_1d - 1)
                        {
                          // Now add the node to the appropriate boundary
                          add_boundary_node(boundary_and_face_index.first,
                                            cube_node_pt);
                        }

                        // Break
                        break;
                    } // switch (boundary_and_face_index.second)
                  } // for (unsigned i_bound=0;i_bound<n_boundaries;i_bound++)
                } // if (el_boundary_information.size()>0)

                // Okay, so this Node does not lie on a face of a quad element
                // that lies on a boundary. However(!), there is still another
                // case in which it might lie on a mesh boundary; when the Node
                // is the ONLY Node in the element that lies on the boundary.
                // It's a slightly odd scenario but can easily occur on meshes
                // created using the QuadFromTriangleMesh class.
                Node* quad_node_pt =
                  quad_el_pt->node_pt(quad_local_node_number);

                // If this Node lies on a boundary
                if (quad_node_pt->is_on_boundary())
                {
                  // If we haven't created the Node yet AND the corresponding
                  // Node in the quad mesh definitely lies on the boundary
                  if (cube_node_pt == 0)
                  {
                    // Create the boundary node
                    cube_node_pt = cube_el_pt->construct_boundary_node(
                      cube_local_node_number, time_stepper_pt);

                    // Indicate that the node has been created now
                    has_node_been_setup[quad_local_node_number] = true;
                  }

                  // Storage for the set of boundaries that this Node lies on
                  std::set<unsigned>* boundaries_pt;

                  // Get the boundaries that this Node lies on
                  quad_node_pt->get_boundaries_pt(boundaries_pt);

                  // Create an iterator to loop over the entries of the set
                  std::set<unsigned>::const_iterator const_it;

                  // Iterate over the boundaries
                  for (const_it = boundaries_pt->begin();
                       const_it != boundaries_pt->end();
                       const_it++)
                  {
                    // Add the node to the appropriate boundary (if it hasn't
                    // already been added)
                    add_boundary_node(*const_it, cube_node_pt);
                  }
                } // quad_node_pt->is_on_boundary()

                // If the node still hasn't been created yet
                if (!has_node_been_setup[quad_local_node_number])
                {
                  // The only other possibility is that it is a regular node
                  cube_node_pt = cube_el_pt->construct_node(
                    cube_local_node_number, time_stepper_pt);

                  // Indicate that the node has now been set up
                  has_node_been_setup[quad_local_node_number] = true;
                }

                // Set the x-coordinate of the node
                cube_node_pt->x(0) = quad_node_pt->x(0);

                // Set the y-coordinate of the node
                cube_node_pt->x(1) = quad_node_pt->x(1);

                // Set the z-coordinate of the node
                cube_node_pt->x(2) = z_value;

                // Add it to global storage
                Node_pt[node_count] = cube_node_pt;

                // Increment the counter
                node_count++;
              } // for (unsigned n_0=0;n_0<n_node_1d;n_0++)
            } // for (unsigned n_1=0;n_1<n_node_1d;n_1++)

            //------------------------------------------------------------------
            // We've set up the nodes in one slice, now set up the edge-to-node
            // information and store it in the map edge_to_nodes_map.
            //------------------------------------------------------------------
            for (unsigned i_edge = 0; i_edge < n_edge; i_edge++)
            {
              // If it hasn't already been done
              if (!has_edge_been_done[i_edge])
              {
                // Create a vector to store the nodes on this edge
                Vector<Node*> edge_nodes_pt(n_node_1d, 0);

                // Loop over the nodes and add them
                for (unsigned i_node = 0; i_node < n_node_1d; i_node++)
                {
                  // Calculate the value of quad_local_node_number depending on
                  // which corner the node lies on
                  switch (i_edge)
                  {
                    case 0:
                      // Northern edge
                      quad_local_node_number =
                        ((n_node_1d * n_node_1d) - 1) - i_node;

                      // Break
                      break;
                    case 1:
                      // Eastern edge
                      quad_local_node_number =
                        (n_node_1d - 1) + (i_node * n_node_1d);

                      // Break
                      break;
                    case 2:
                      // Southern edge
                      quad_local_node_number = i_node;

                      // Break
                      break;
                    case 3:
                      // Western edge
                      quad_local_node_number =
                        (n_node_1d * (n_node_1d - 1)) - (i_node * n_node_1d);

                      // Break
                      break;
                  }

                  // Index of the local node number
                  cube_local_node_number =
                    quad_local_node_number + (n_node_1d * n_node_1d * n_2);

                  // Store the i_node-th node along this edge
                  edge_nodes_pt[i_node] =
                    cube_el_pt->node_pt(cube_local_node_number);
                } // for (unsigned i_node=0;i_node<n_node_1d;i_node++)

                // Store the nodes along the north edge of the element
                edge_to_nodes_map[edges_and_z_value[i_edge]] = edge_nodes_pt;
              } // if (!has_edge_been_done[i_edge])
            } // for (unsigned i_edge=0;i_edge<n_edge;i_edge++)

            //--------------------------------------------------------------------
            // Now store the corner nodes of the element at the current node
            // slice in the map quad_corner_node_to_cube_node_map.
            //--------------------------------------------------------------------
            // Loop over the corner nodes
            for (unsigned i_corner = 0; i_corner < n_corner; i_corner++)
            {
              // Calculate the value of quad_local_node_number depending on
              // which corner the node lies on
              switch (i_corner)
              {
                case 0:
                  // Bottom-left corner
                  quad_local_node_number = 0;

                  // Break
                  break;
                case 1:
                  // Bottom-right corner
                  quad_local_node_number = n_node_1d - 1;

                  // Break
                  break;
                case 2:
                  // Top-left corner
                  quad_local_node_number = n_node_1d * (n_node_1d - 1);

                  // Break
                  break;
                case 3:
                  // Top-right corner
                  quad_local_node_number = (n_node_1d * n_node_1d) - 1;

                  // Break
                  break;
              }

              // Index of the corresponding corner node (in the cube element)
              cube_local_node_number =
                quad_local_node_number + (n_node_1d * n_node_1d * n_2);

              // Get a pointer to the first corner node (in the quad mesh)
              Node* quad_corner_node_pt =
                quad_el_pt->node_pt(quad_local_node_number);

              // Get a pointer to the first corner node (in the cube mesh)
              Node* cube_corner_node_pt =
                cube_el_pt->node_pt(cube_local_node_number);

              // Pair up the corner node and z-value
              std::pair<Node*, double> quad_corner_node_and_z_value(
                quad_corner_node_pt, z_value);

              // Store it in the map
              quad_corner_node_to_cube_node_map[quad_corner_node_and_z_value] =
                cube_corner_node_pt;
            } // for (unsigned i_corner=0;i_corner<n_corner;i_corner++)

#ifdef PARANOID
            // Loop over all the nodes in the current node slice
            for (unsigned i_node = 0; i_node < (n_node_1d * n_node_1d);
                 i_node++)
            {
              // Make sure every node in the current slice has been set up
              if (!has_node_been_setup[i_node])
              {
                // Create an ostringstream object to create an error message
                std::ostringstream error_message_stream;

                // Create the error message
                error_message_stream << "There are nodes in element " << e
                                     << " which have not been constructed!"
                                     << std::endl;

                // Throw an error
                throw OomphLibError(error_message_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            } // for (unsigned i_node=0;i_node<(n_node_1d*n_node_1d);i_node++)
#endif
          } // if ((i>0)&&(n_2==0))
        } // for (unsigned n_2=0;n_2<n_node_1d;n_2++)
      } // for (unsigned e=0;e<n_element_in_quad_mesh;e++)

      // Sanity check; make sure enough nodes have been stored
      if (node_count != (1 + (n_node_1d - 1) * (i + 1)) * n_node_in_quad_mesh)
      {
        // Create an ostringstream object to create an error message
        std::ostringstream error_message_stream;

        // The number of nodes we expect there to be in the mesh
        int expected_n_node =
          (1 + (n_node_1d - 1) * (i + 1)) * n_node_in_quad_mesh;

        // Difference in node count
        int node_diff = expected_n_node - node_count;

        // Create the error message
        error_message_stream
          << "There are meant to be " << expected_n_node
          << " nodes in the extruded mesh by the\nend of slice " << i
          << " but you have " << node_count << ". "
          << "That's a difference of " << node_diff << "!";

        // Throw an error
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    } // for (unsigned i=0;i<Nz;i++)

    //------------------------------------------------------------------------
    // Count the number of nodes on the boundaries of the quad mesh; there
    // should really be a helper function that does this (hint hint)
    //------------------------------------------------------------------------
    // Use a set to make sure we don't count edge/corner nodes several times
    std::set<const Node*> all_quad_boundary_nodes_pt;

    // The number of boundary nodes in our extruded mesh
    int n_quad_mesh_boundary_node = 0;

    // The number of boundaries in the quad mesh
    for (unsigned b = 0; b < n_boundary_in_quad_mesh; b++)
    {
      // Add on the number of nodes on the b-th boundary
      unsigned n_boundary_node = quad_mesh_pt->nboundary_node(b);

      // Loop over the nodes on this boundary
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        // We haven't come across this boundary node yet
        if (all_quad_boundary_nodes_pt.find(quad_mesh_pt->boundary_node_pt(
              b, n)) == all_quad_boundary_nodes_pt.end())
        {
          // Count this node
          n_quad_mesh_boundary_node++;

          // Now remember it so we don't count it twice
          all_quad_boundary_nodes_pt.insert(
            quad_mesh_pt->boundary_node_pt(b, n));
        }
      } // for (unsigned n=0; n<n_boundary_node; n++)
    } // for (unsigned b=0; b<n_boundary_in_quad_mesh; b++)

    //------------------------------------------------------------------------
    // Count the number of nodes on the boundaries of the extruded cube mesh;
    // (recall hint hint from above...)
    //------------------------------------------------------------------------
    // Number of boundary nodes expected in the extruded mesh; initialise to
    // the number of boundary nodes on the front and back faces
    int n_expected_boundary_node = 2 * quad_mesh_pt->nnode();

    // Add on the nodes expected on the extruded mesh boundaries (but ignoring
    // the boundary nodes on the front and back faces)
    n_expected_boundary_node +=
      (((n_node_1d - 1) * Nz) - 1) * n_quad_mesh_boundary_node;

    // The number of boundary nodes in our extruded mesh
    int n_extruded_mesh_boundary_node = 0;

    // Use a set to make sure we don't count edge/corner nodes several times
    std::set<const Node*> all_extruded_boundary_nodes_pt;

    // Loop over the boundaries of the mesh
    for (unsigned b = 0; b < n_boundary_in_extruded_mesh; b++)
    {
      // Add on the number of nodes on the b-th boundary
      unsigned n_boundary_node = nboundary_node(b);

      // Loop over the nodes on this boundary
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        // We haven't come across this boundary node yet
        if (all_extruded_boundary_nodes_pt.find(this->boundary_node_pt(b, n)) ==
            all_extruded_boundary_nodes_pt.end())
        {
          // Count this node
          n_extruded_mesh_boundary_node++;

          // Now remember it so we don't count it twice
          all_extruded_boundary_nodes_pt.insert(this->boundary_node_pt(b, n));
        }
      } // for (unsigned n=0; n<n_boundary_node; n++)
    } // for (unsigned b=0; b<n_boundary_in_extruded_mesh; b++)

    //------------------------------------------------------------------------
    // Check that there is the correct number of boundary nodes
    //------------------------------------------------------------------------
    // Error: we don't have the correct number of boundary nodes
    if (n_extruded_mesh_boundary_node != n_expected_boundary_node)
    {
      // Create an ostringstream object to create an error message
      std::ostringstream error_message_stream;

      // Difference in node count
      int node_diff = n_expected_boundary_node - n_extruded_mesh_boundary_node;

      // Create the error message
      error_message_stream << "There should be " << n_expected_boundary_node
                           << " boundary nodes in the extruded mesh but there"
                           << "\nare only " << n_extruded_mesh_boundary_node
                           << " boundary nodes. That's a difference of "
                           << node_diff << "!" << std::endl;

      // Throw an error
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // If the user wishes the mesh setup time to be doc-ed
    if (MeshExtrusionHelpers::Mesh_extrusion_helper.doc_mesh_setup_time())
    {
      // Tell the user
      oomph_info << "Time taken to extrude mesh [sec]: "
                 << TimingHelpers::timer() - t_mesh_construct_start
                 << std::endl;
    }

    // Get the current (i.e. start) time
    double t_mesh_macro_start = TimingHelpers::timer();

    //----------------------------------------------------------------
    // If an element has a macro-element representation in the quad
    // mesh then it need an extruded macro-element representation in
    // the extruded mesh so set that up here. Then, use the extruded
    // macro-element representation of the elements to move all the
    // nodes that need to be moved.
    //----------------------------------------------------------------
    // Create a map that takes a Domain and returns the corresponding
    // ExtrudedDomain object
    std::map<Domain*, ExtrudedDomain*> domain_to_extruded_domain_map;

    // Loop over the elements in the input mesh
    for (unsigned e = 0; e < n_element_in_quad_mesh; e++)
    {
      // Get a pointer to the e-th element in the input mesh.
      // NOTE: The element is upcast to the QElementBase class so we
      // can access the s_macro_ll() and s_macro_ur() functions.
      QElementBase* quad_el_pt =
        dynamic_cast<QElementBase*>(quad_mesh_pt->finite_element_pt(e));

      // Check if the element has a macro-element representation
      if (quad_el_pt->macro_elem_pt() != 0)
      {
        // Get a pointer to the corresponding Domain object
        Domain* domain_pt = quad_el_pt->macro_elem_pt()->domain_pt();

        // How many macro elements are there in this Domain?
        unsigned n_macro_element = domain_pt->nmacro_element();

        // Create a pointer to point to the extruded version of the above Domain
        ExtrudedDomain* extruded_domain_pt = 0;

        // Create an iterator for the map
        std::map<Domain*, ExtrudedDomain*>::iterator it;

        // Check if there is a corresponding entry in the map
        it = domain_to_extruded_domain_map.find(domain_pt);

        // Check if the extruded version of this domain has already been created
        if (it != domain_to_extruded_domain_map.end())
        {
          // Get the associated ExtrudedDomain pointer
          extruded_domain_pt = it->second;
        }
        // If there is no entry then we need to create the extruded domain
        else
        {
          // Create an ExtrudedDomain object associated with the mesh
          extruded_domain_pt = new ExtrudedDomain(domain_pt, Nz, Zmin, Zmax);

          // Now store it in the map
          domain_to_extruded_domain_map[domain_pt] = extruded_domain_pt;

          // Add it to the mesh's list of extruded domains
          Extruded_domain_pt.push_back(extruded_domain_pt);
        }

        // Loop over each element slice in the z-direction
        for (unsigned i = 0; i < Nz; i++)
        {
          // What is the ID number of the macro element associated with
          // this quad element?
          unsigned macro_id =
            (i * n_macro_element +
             quad_el_pt->macro_elem_pt()->macro_element_number());

          // Get a pointer to the extruded macro element
          ExtrudedMacroElement* extruded_macro_element_pt =
            extruded_domain_pt->macro_element_pt(macro_id);

          // Get the index of the cube element in the extruded mesh
          unsigned element_index = i * n_element_in_quad_mesh + e;

          // Get the corresponding extruded element from the global storage
          ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Element_pt[element_index]);

          // Set its macro element pointer
          el_pt->set_macro_elem_pt(extruded_macro_element_pt);

          // The number of dimensions in the extruded mesh
          unsigned n_dim = 3;

          // Loop over the coordinate directions
          for (unsigned i = 0; i < n_dim - 1; i++)
          {
            // The i-th local coordinate value of the "lower-left" corner of
            // the current element
            el_pt->s_macro_ll(i) = quad_el_pt->s_macro_ll(i);

            // The i-th local coordinate value of the "upper-right" corner of
            // the current element
            el_pt->s_macro_ur(i) = quad_el_pt->s_macro_ur(i);
          }

          // The local coordinate value of the "lower-left" corner which
          // corresponds to the time direction
          el_pt->s_macro_ll(n_dim - 1) = -1.0;

          // The local coordinate value of the "upper-right" corner which
          // corresponds to the time direction
          el_pt->s_macro_ur(n_dim - 1) = 1.0;
        } // for (unsigned i=0;i<Nz;i++)
      } // if (quad_el_pt->macro_element_pt()!=0)
    } // for (unsigned e=0;e<n_element_in_quad_mesh;e++)

    // If the user wishes the mesh setup time to be doc-ed
    if (MeshExtrusionHelpers::Mesh_extrusion_helper.doc_mesh_setup_time())
    {
      // Tell the user
      oomph_info << "Time taken to set up extruded macro-elements [sec]: "
                 << TimingHelpers::timer() - t_mesh_macro_start << std::endl;
    }

    // Call the node update function to do all the relevant node moving
    node_update();

#ifdef PARANOID
    // Make sure there are no duplicate nodes (i.e. two or more nodes that
    // occupy the same Eulerian position)
    if (check_for_repeated_nodes())
    {
      // Throw a warning
      OomphLibWarning("The mesh contains repeated nodes!",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
    }

    // Also make sure enough nodes have been been created
    if (node_count != n_node_in_cube_mesh)
    {
      // Create an ostringstream object to create an error message
      std::ostringstream error_message_stream;

      // Create the error message
      error_message_stream << "An incorrect number of nodes have been created! "
                           << "There are\nmeant to be " << n_node_in_cube_mesh
                           << " nodes in the extruded mesh but\nonly "
                           << node_count << " have actually been constructed!"
                           << std::endl;

      // Throw a warning
      OomphLibWarning(error_message_stream.str(),
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Setup lookup scheme that establishes which elements are located
    // on the mesh boundaries
    setup_boundary_element_info();

    // If the user wishes the mesh setup time to be doc-ed
    if (MeshExtrusionHelpers::Mesh_extrusion_helper.doc_mesh_setup_time())
    {
      // Tell the user
      oomph_info << "Total time for extruded mesh creation/setup [sec]: "
                 << TimingHelpers::timer() - t_mesh_setup_start << std::endl;
    }
  } // End of build_mesh
} // End of namespace oomph
