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
// Driver for 2D moving block
#ifndef OOMPH_QUAD_FROM_TRIANGLE_MESH_TEMPLATE_CC
#define OOMPH_QUAD_FROM_TRIANGLE_MESH_TEMPLATE_CC

// The mesh
#include "quad_from_triangle_mesh.template.h"

using namespace std;
using namespace oomph;


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


namespace oomph
{
  //======================================================================
  /// \short Build the full mesh with the help of the scaffold mesh coming
  /// from the triangle mesh generator, Triangle. To build this quad
  /// element based mesh we make use of the fact that a triangle element
  /// can be split as shown in the diagram below:
  ///
  ///                  N2
  ///                 | |                 N0 : 1st scaffold element node
  ///                |   |                N1 : 2nd scaffold element node
  ///               |     |               N2 : 3rd scaffold element node
  ///              |       |
  ///           C |   Q_2   | B           Edge 0 : N0 --> N1
  ///            | |       | |            Edge 1 : N1 --> N2
  ///           |   |     |   |           Edge 2 : N2 --> N0
  ///          |     |   |     |
  ///         |       | |       |         A : Midpoint of edge 0
  ///        |   Q_0   |   Q_1   |        B : Midpoint of edge 1
  ///       |          |          |       C : Midpoint of edge 2
  ///      |           |           |
  ///     N0 __________|__________ N1
  ///                  A
  ///
  /// The intersection of all three quad elements is the centroid. Using
  /// this splitting, the subsequent mesh will consist of quadrilaterals
  /// whose shape which depend on the structure of the underlying mesh.
  //======================================================================
  template<class ELEMENT>
  void QuadFromTriangleMesh<ELEMENT>::build_from_scaffold(
    TriangleScaffoldMesh* tmp_mesh_pt,
    TimeStepper* time_stepper_pt,
    const bool& use_attributes)
  {
    // Create space for elements
    unsigned nelem = tmp_mesh_pt->nelement();

    // We will have 3 quad elements per scaffold element
    Element_pt.resize(3 * nelem);

    // Set number of boundaries
    unsigned nbound = tmp_mesh_pt->nboundary();

    // Resize the boundary information (the number of boundaries doesn't
    // change)
    set_nboundary(nbound);

    // Stores each element attached to a boundary and the index of the
    // face of the given element attached to the boundary
    Boundary_element_pt.resize(nbound);
    Face_index_at_boundary.resize(nbound);

    // Create a quad element for nodal data
    ELEMENT* temp_el_pt = new ELEMENT;

    // Get the number of nodes in a quad element
    unsigned nnode_el = temp_el_pt->nnode();

    // Find the number of nodes along one edge of a quad element
    unsigned nnode_1d = temp_el_pt->nnode_1d();

    // Calculate the number of nodes that will lie along an edge of a
    // triangle element in the scaffold mesh
    unsigned nnode_edge = 2 * nnode_1d - 1;

    // Delete the element pointer
    delete temp_el_pt;

    // Make it a null pointer
    temp_el_pt = 0;

    // Create dummy linear quad for geometry
    QElement<2, 2> dummy_element;

    // The dimension of the element
    unsigned n_dim = 2;

    // The position type
    unsigned n_position_type = 1;

    // Don't assign memory for any values
    unsigned initial_n_value = 0;

    // Loop over the nodes of the element and make them
    for (unsigned j = 0; j < 4; j++)
    {
      dummy_element.node_pt(j) =
        new Node(n_dim, n_position_type, initial_n_value);
    }

    // Local node number of each quad element corner
    unsigned corner_0 = 0;
    unsigned corner_1 = nnode_1d - 1;
    unsigned corner_2 = nnode_el - nnode_1d;
    unsigned corner_3 = nnode_el - 1;

    // Create a map to return a vector of pointers to nnode_1d nodes where
    // the input is an edge. If the edge hasn't been set up then this will
    // return a null pointer. Note: all node pointers on an edge will be
    // stored in clockwise ordering. Therefore, to copy the data of an
    // edge into the adjoining element we must proceed through the vector
    // backwards (as progressing through an edge of an element in a clockwise
    // manner is equivalent to proceeding through the edge of the neighbouring
    // element in an anti-clockwise manner)
    std::map<Edge, Vector<Node*>> edge_nodes_map;

    // Set up a map to check if the scaffold mesh node has been set up in the
    // quad mesh. If the node has been set up this map will return a pointer
    // to it otherwise it will return a null pointer
    std::map<Node*, Node*> scaffold_to_quad_mesh_node;

    // Loop over elements in scaffold mesh
    unsigned new_el_count = 0;

    // Create storage for the coordinates of the centroid
    Vector<double> centroid(2);

    // Create storage for the coordinates of the vertices of the triangle
    Vector<Vector<double>> triangle_vertex(3);

    // Loop over all of the elements in the scaffold mesh
    for (unsigned e = 0; e < nelem; e++)
    {
      // Initialise centroid values for the e-th triangle element
      centroid[0] = 0.0;
      centroid[1] = 0.0;

      // Loop over the scaffold element nodes
      for (unsigned j = 0; j < 3; j++)
      {
        // Resize the j-th triangle_vertex entry to contain the x and
        // y-coordinate
        triangle_vertex[j].resize(2);

        // Get the coordinates
        double x = tmp_mesh_pt->finite_element_pt(e)->node_pt(j)->x(0);
        double y = tmp_mesh_pt->finite_element_pt(e)->node_pt(j)->x(1);

        // Increment the centroid coordinates
        centroid[0] += x;
        centroid[1] += y;

        // Assign the triangle_vertex coordinates
        triangle_vertex[j][0] = x;
        triangle_vertex[j][1] = y;
      }

      // Divide the centroid entries by 3 to get the centroid coordinates
      centroid[0] /= 3.0;
      centroid[1] /= 3.0;

      // Create element pointers and assign them to a vector
      //----------------------------------------------------
      // Make new quad elements of the type specified by the template parameter
      ELEMENT* el0_pt = new ELEMENT;
      ELEMENT* el1_pt = new ELEMENT;
      ELEMENT* el2_pt = new ELEMENT;

      // Create a vector of ELEMENTs to store el0_pt, el1_pt and el2_pt
      Vector<ELEMENT*> el_vector_pt(3, 0);

      // Assign the entries to el_vector_pt
      el_vector_pt[0] = el0_pt;
      el_vector_pt[1] = el1_pt;
      el_vector_pt[2] = el2_pt;


      // Create the first node in each quad element and store in Node_pt.
      // These correspond to the nodes of the simplex triangle stored in
      // Tmp_mesh_pt. If they have already been set up then we do nothing:
      //------------------------------------------------------------------
      // Loop over the scaffold element nodes and check to see if they have
      // been set up
      for (unsigned j = 0; j < 3; j++)
      {
        // Pointer to node in the scaffold mesh
        Node* scaffold_node_pt = tmp_mesh_pt->finite_element_pt(e)->node_pt(j);

        // Check if the node has been set up yet
        Node* qmesh_node_pt = scaffold_to_quad_mesh_node[scaffold_node_pt];

        // Haven't done this one yet
        if (qmesh_node_pt == 0)
        {
          // Get pointer to set of mesh boundaries that this
          // scaffold node occupies; NULL if the node is not on any boundary
          std::set<unsigned>* boundaries_pt;
          scaffold_node_pt->get_boundaries_pt(boundaries_pt);

          // Check to see if it's on any boundaries
          if (boundaries_pt != 0)
          {
            // Create new boundary node. The scaffold element nodes are the
            // corners of a simplex triangle and thus always correspond to the
            // first node in each quad element
            qmesh_node_pt = el_vector_pt[j]->construct_boundary_node(
              corner_0, time_stepper_pt);

            // Add to boundaries
            for (std::set<unsigned>::iterator it = boundaries_pt->begin();
                 it != boundaries_pt->end();
                 ++it)
            {
              add_boundary_node(*it, qmesh_node_pt);
            }
          }
          // Build normal node
          else
          {
            // Create new normal node
            qmesh_node_pt =
              el_vector_pt[j]->construct_node(corner_0, time_stepper_pt);
          }

          // Add the mapping from the scaffold mesh node to the quad mesh node
          scaffold_to_quad_mesh_node[scaffold_node_pt] = qmesh_node_pt;

          // Copy new node, created using the NEW element's construct_node
          // function into global storage, using the same global
          // node number that we've just associated with the
          // corresponding node in the scaffold mesh
          Node_pt.push_back(qmesh_node_pt);
        }
        // If this node has already been done we need to copy the data across
        else
        {
          el_vector_pt[j]->node_pt(corner_0) = qmesh_node_pt;
        }


        // Set global coordinate
        el_vector_pt[j]->node_pt(corner_0)->x(0) = triangle_vertex[j][0];
        el_vector_pt[j]->node_pt(corner_0)->x(1) = triangle_vertex[j][1];
      }


      // Create the edges of the scaffold element and check to see if
      // they've been set up yet or not. If they haven't:
      //--------------------------------------------------------------
      // Make the three edges of the triangle
      Edge edge0(tmp_mesh_pt->finite_element_pt(e)->node_pt(0),
                 tmp_mesh_pt->finite_element_pt(e)->node_pt(1));
      Edge edge1(tmp_mesh_pt->finite_element_pt(e)->node_pt(1),
                 tmp_mesh_pt->finite_element_pt(e)->node_pt(2));
      Edge edge2(tmp_mesh_pt->finite_element_pt(e)->node_pt(2),
                 tmp_mesh_pt->finite_element_pt(e)->node_pt(0));

      // Check if the edges have been set up (each will have size nnode_1d).
      // If they have not been set up yet, this will
      Vector<Node*> edge0_vector_pt = edge_nodes_map[edge0];
      Vector<Node*> edge1_vector_pt = edge_nodes_map[edge1];
      Vector<Node*> edge2_vector_pt = edge_nodes_map[edge2];

      // Bools to indicate whether or not the edges have been set up
      bool edge0_setup = (edge0_vector_pt.size() != 0);
      bool edge1_setup = (edge1_vector_pt.size() != 0);
      bool edge2_setup = (edge2_vector_pt.size() != 0);

      // If edge 0 hasn't been set up (node 0 to node 1)
      if (!edge0_setup)
      {
        // Resize the vector to have length nnode_1d
        edge0_vector_pt.resize(nnode_edge, 0);

        // First node along edge 0 is the first node of element 0
        edge0_vector_pt[0] = el_vector_pt[0]->node_pt(0);

        // Last node along edge 0 is the first node of element 1
        edge0_vector_pt[nnode_edge - 1] = el_vector_pt[1]->node_pt(0);
      }

      // If edge 1 hasn't been set up (node 1 to node 2)
      if (!edge1_setup)
      {
        // Resize the vector to have length nnode_1d
        edge1_vector_pt.resize(nnode_edge, 0);

        // First node along edge 1 is the first node of element 1
        edge1_vector_pt[0] = el_vector_pt[1]->node_pt(0);

        // Last node along edge 1 is the first node of element 2
        edge1_vector_pt[nnode_edge - 1] = el_vector_pt[2]->node_pt(0);
      }

      // If edge 2 hasn't been set up (node 2 to node 0)
      if (!edge2_setup)
      {
        // Resize the vector to have length nnode_1d
        edge2_vector_pt.resize(nnode_edge, 0);

        // First node along edge 2 is the first node of element 2
        edge2_vector_pt[0] = el_vector_pt[2]->node_pt(0);

        // Last node along edge 2 is the first node of element 0
        edge2_vector_pt[nnode_edge - 1] = el_vector_pt[0]->node_pt(0);
      }


#ifdef PARANOID
      // If any of the edges have been set up, make sure that that the endpoints
      // in the returned vectors have the same address as those on the vertices

      // Come back and finish this off.
      // To check:
      //    - If two edges which have been set up have the same node in the
      //    middle
      //    - If an edge has already been set up then the map will return the
      //      same node as in the vector
#endif

      // Boundary IDs for bottom and left edge of quad
      // from scaffold mesh (if these remain zero the edges
      // are not on a boundary)
      unsigned q0_lower_boundary_id = 0;
      unsigned q0_left_boundary_id = 0;
      unsigned q1_lower_boundary_id = 0;
      unsigned q1_left_boundary_id = 0;
      unsigned q2_lower_boundary_id = 0;
      unsigned q2_left_boundary_id = 0;

      // Lower/left boundary IDs for quad 0; the lower edge in quad 0 is on
      // edge 0 of the scaffold triangle and the left edge in quad is on edge
      // 2 in scaffold triangle
      q0_lower_boundary_id = tmp_mesh_pt->edge_boundary(e, 0);
      q0_left_boundary_id = tmp_mesh_pt->edge_boundary(e, 2);

      // Lower/left boundary IDs for quad 1; the lower edge in quad 1 is on
      // edge 1 of the scaffold triangle and the left edge in quad is on edge
      // 0 of the scaffold triangle
      q1_lower_boundary_id = tmp_mesh_pt->edge_boundary(e, 1);
      q1_left_boundary_id = tmp_mesh_pt->edge_boundary(e, 0);

      // Lower/left boundary IDs for quad 2; the lower edge in quad 2 is on
      // edge 2 of the scaffold triangle and the left edge in quad is on edge
      // 1 of the scaffold triangle
      q2_lower_boundary_id = tmp_mesh_pt->edge_boundary(e, 2);
      q2_left_boundary_id = tmp_mesh_pt->edge_boundary(e, 1);

      // Store the boundary IDs as a vector; allows us to loop over them easily
      Vector<unsigned> boundary_id_vector(6, 0);
      boundary_id_vector[0] = q0_lower_boundary_id;
      boundary_id_vector[1] = q0_left_boundary_id;
      boundary_id_vector[2] = q1_lower_boundary_id;
      boundary_id_vector[3] = q1_left_boundary_id;
      boundary_id_vector[4] = q2_lower_boundary_id;
      boundary_id_vector[5] = q2_left_boundary_id;

      // Loop over the quad elements and store the boundary elements in the
      // vector Boundary_element_pt
      for (unsigned j = 0; j < 3; j++)
      {
        // Loop over the lower and the left boundary ID in the j'th element
        for (unsigned k = 0; k < 2; k++)
        {
          // The quad element lies on a boundary of the mesh
          if (boundary_id_vector[2 * j + k] > 0)
          {
            // Since the j'th quad element lies on a boundary of the mesh we add
            // a pointer to the element to the appropriate entry of
            // Boundary_element_pt
            Boundary_element_pt[boundary_id_vector[2 * j + k] - 1].push_back(
              el_vector_pt[j]);

            // If k=0 then the lower boundary of the quad element lies on
            // the boundary of the mesh and if k=1 then the left boundary
            // of the quad element lies on the mesh boundary. For quad elements
            // the indices are as follows:
            //       North face: 2
            //       East face: 1
            //       South face: -2
            //       West face: -1
            if (k == 0)
            {
              Face_index_at_boundary[boundary_id_vector[2 * j + k] - 1]
                .push_back(-2);
            }
            else
            {
              Face_index_at_boundary[boundary_id_vector[2 * j + k] - 1]
                .push_back(-1);
            } // if (k==0)
          } // if (boundary_id_vector[2*j+k]>0)
        } // for (unsigned k=0;k<2;k++)
      } // for (unsigned j=0;j<3;j++)


      // The upper right node is always the centroid. Note: The 'corner_3' node
      // lies within each of the three quad elements so we simply share the
      // pointer to it with each element:
      //---------------------------------------------------------------------------
      // Construct the centroid node
      Node* nod_pt = el0_pt->construct_node(corner_3, time_stepper_pt);

      // Add the pointer to the vector of nodal pointers
      Node_pt.push_back(nod_pt);

      // Quad 0
      el0_pt->node_pt(corner_3)->x(0) = centroid[0];
      el0_pt->node_pt(corner_3)->x(1) = centroid[1];

      // Quad 1
      el1_pt->node_pt(corner_3) = el0_pt->node_pt(corner_3);

      // Quad 2
      el2_pt->node_pt(corner_3) = el0_pt->node_pt(corner_3);


      // Set the nodal positions of the dummy element to emulate the FIRST
      // quad element (this allows for simple linear interpolation later):
      //------------------------------------------------------------------
      // Bottom-left corner
      dummy_element.node_pt(0)->x(0) = triangle_vertex[0][0];
      dummy_element.node_pt(0)->x(1) = triangle_vertex[0][1];

      // Bottom-right corner
      dummy_element.node_pt(1)->x(0) =
        0.5 * (triangle_vertex[0][0] + triangle_vertex[1][0]);
      dummy_element.node_pt(1)->x(1) =
        0.5 * (triangle_vertex[0][1] + triangle_vertex[1][1]);

      // Top-left corner
      dummy_element.node_pt(2)->x(0) =
        0.5 * (triangle_vertex[0][0] + triangle_vertex[2][0]);
      dummy_element.node_pt(2)->x(1) =
        0.5 * (triangle_vertex[0][1] + triangle_vertex[2][1]);

      // Top-right corner
      dummy_element.node_pt(3)->x(0) = centroid[0];
      dummy_element.node_pt(3)->x(1) = centroid[1];


      // Set up all of the nodes in the first quad element (Q0):
      //--------------------------------------------------------
      // Local and global coordinate vectors for the nodes
      Vector<double> s(2);
      Vector<double> x(2);

      // Loop over all of nodes in Q0 noting that the lower left corner node
      // (node 0) and the upper right corner node (centroid) have already
      // been set up
      for (unsigned j = 1; j < corner_3; j++)
      {
        // Indicates whether or not the node has been set up yet
        bool done = false;

        // On the lower edge
        if (j < nnode_1d)
        {
          // If the lower edge has already been set up then we already know the
          // nodes along this edge
          if (edge0_setup)
          {
            // The j'th node along this edge is the (nnode_1d-j)'th node in the
            // vector (remembering that the ordering is backwards since it has
            // already been set up)
            el0_pt->node_pt(j) = edge0_vector_pt[(nnode_edge - 1) - j];

            // Since the node has already been set up we do not need to sort
            // out its global coordinate data so skip to the next node
            continue;
          }
          // If the edge hasn't been set up yet
          else
          {
            // If the node lies on a boundary too then we need to construct a
            // boundary node
            if (q0_lower_boundary_id > 0)
            {
              // Construct a boundary node
              Node* nod_pt =
                el0_pt->construct_boundary_node(j, time_stepper_pt);

              // Add it to the list of boundary nodes
              add_boundary_node(q0_lower_boundary_id - 1, nod_pt);

              // Add the node into the vector of nodes on edge 0
              edge0_vector_pt[j] = nod_pt;

              // Add it to the list of nodes in the mesh
              Node_pt.push_back(nod_pt);

              // Indicate the j'th node has been set up
              done = true;
            }
          }

          // Node is not on a mesh boundary but on the lower edge
          if (!done)
          {
            // Construct a normal node
            Node* nod_pt = el0_pt->construct_node(j, time_stepper_pt);

            // Add the node into the vector of nodes on edge 0
            edge0_vector_pt[j] = nod_pt;

            // Add it to the list of nodes in the mesh
            Node_pt.push_back(nod_pt);

            // Indicate the j'th node has been set up
            done = true;
          }
        }
        // On the left edge
        else if (j % nnode_1d == 0)
        {
          // If the left edge has already been set up then we already know the
          // nodes along this edge
          if (edge2_setup)
          {
            // The j'th node is the (j/nnode_1d)'th node along this edge and
            // thus the (j/nnode_1d)'th entry in the edge vector
            el0_pt->node_pt(j) = edge2_vector_pt[j / nnode_1d];

            // Since the node has already been set up we do not need to sort
            // out its global coordinate data
            continue;
          }
          // If the edge hasn't been set up yet
          else
          {
            if (q0_left_boundary_id > 0)
            {
              // Construct a boundary node
              Node* nod_pt =
                el0_pt->construct_boundary_node(j, time_stepper_pt);

              // Add it to the list of boundary nodes
              add_boundary_node(q0_left_boundary_id - 1, nod_pt);

              // Add the node into the vector of nodes on edge 2 in clockwise
              // order
              edge2_vector_pt[(nnode_edge - 1) - (j / nnode_1d)] = nod_pt;

              // Add it to the list of nodes in the mesh
              Node_pt.push_back(nod_pt);

              // Indicate that the j'th node has been set up
              done = true;
            }
          }

          // Node is not on a mesh boundary but on the left edge
          if (!done)
          {
            // Construct a normal node
            Node* nod_pt = el0_pt->construct_node(j, time_stepper_pt);

            // Add the node into the vector of nodes on edge 2 in clockwise
            // order
            edge2_vector_pt[(nnode_edge - 1) - (j / nnode_1d)] = nod_pt;

            // Add it to the list of nodes in the mesh
            Node_pt.push_back(nod_pt);

            // Indicate the j'th node has been set up
            done = true;
          }
        }

        // Node is not on a mesh boundary or on the edge of the scaffold element
        if (!done)
        {
          // Construct a normal node
          Node* nod_pt = el0_pt->construct_node(j, time_stepper_pt);

          // Add it to the list of nodes in the mesh
          Node_pt.push_back(nod_pt);
        }

        // Get local coordinate
        el0_pt->local_coordinate_of_node(j, s);

        // Interpolate position linearly between vertex nodes
        dummy_element.interpolated_x(s, x);
        el0_pt->node_pt(j)->x(0) = x[0];
        el0_pt->node_pt(j)->x(1) = x[1];
      }


      // Set the nodal positions of the dummy element to now emulate the
      // SECOND quad element:
      //------------------------------------------------------------------
      // Note: we do not need to change the top-right corner since it always
      // coincides with the centroid of the triangle element

      // Bottom-left corner
      dummy_element.node_pt(0)->x(0) = triangle_vertex[1][0];
      dummy_element.node_pt(0)->x(1) = triangle_vertex[1][1];

      // Bottom-right corner
      dummy_element.node_pt(1)->x(0) =
        0.5 * (triangle_vertex[1][0] + triangle_vertex[2][0]);
      dummy_element.node_pt(1)->x(1) =
        0.5 * (triangle_vertex[1][1] + triangle_vertex[2][1]);

      // Top-left corner
      dummy_element.node_pt(2)->x(0) =
        0.5 * (triangle_vertex[0][0] + triangle_vertex[1][0]);
      dummy_element.node_pt(2)->x(1) =
        0.5 * (triangle_vertex[0][1] + triangle_vertex[1][1]);


      // Set up all of the nodes in the second quad element (Q1):
      //--------------------------------------------------------
      // Here we need to notice that we have already set up the final nnode_1d
      // nodes (the upper edge of Q1 coincides with the right edge of Q0)

      // Loop over nodes 1 to (corner_2-1) in Q1 noting that the lower left
      // corner node (node 0) and the upper edge of Q1 contains nodes
      // corner_2 to corner_3
      for (unsigned j = 1; j < corner_2; j++)
      {
        // Indicates whether or not the node has been set up yet
        bool done = false;

        // On the lower edge
        if (j < nnode_1d)
        {
          // If the lower edge has already been set up then we already know the
          // nodes along this edge
          if (edge1_setup)
          {
            // The j'th node along this edge is the (nnode_1d-j)'th node in the
            // vector (remembering that the ordering is backwards if it has
            // already been set up)
            el1_pt->node_pt(j) = edge1_vector_pt[(nnode_edge - 1) - j];

            // Since the node has already been set up we do not need to sort
            // out its global coordinate data
            continue;
          }
          // If the edge hasn't been set up yet
          else
          {
            // If the node lies on a boundary too then we need to construct a
            // boundary node
            if (q1_lower_boundary_id > 0)
            {
              // Construct a boundary node
              Node* nod_pt =
                el1_pt->construct_boundary_node(j, time_stepper_pt);

              // Add it to the list of boundary nodes
              add_boundary_node(q1_lower_boundary_id - 1, nod_pt);

              // Add the node into the vector of nodes on edge 1
              edge1_vector_pt[j] = nod_pt;

              // Add it to the list of nodes in the mesh
              Node_pt.push_back(nod_pt);

              // Indicate the j'th node has been set up
              done = true;
            }
          }

          // Node is not on a mesh boundary but on the lower edge
          if (!done)
          {
            // Construct a normal node
            Node* nod_pt = el1_pt->construct_node(j, time_stepper_pt);

            // Add the node into the vector of nodes on edge 1
            edge1_vector_pt[j] = nod_pt;

            // Add it to the list of nodes in the mesh
            Node_pt.push_back(nod_pt);

            // Indicate the j'th node has been set up
            done = true;
          }
        }
        // On the left edge
        else if (j % nnode_1d == 0)
        {
          // If the left edge has already been set up then we already know the
          // nodes along this edge
          if (edge0_setup)
          {
            // The j'th node along this edge is the (j/nnode_1d)'th node in the
            // vector
            el1_pt->node_pt(j) = edge0_vector_pt[j / nnode_1d];

            // Since the node has already been set up we do not need to sort
            // out its global coordinate data
            continue;
          }
          // If the edge hasn't been set up yet
          else
          {
            if (q1_left_boundary_id > 0)
            {
              // Construct a boundary node
              Node* nod_pt =
                el1_pt->construct_boundary_node(j, time_stepper_pt);

              // Add it to the list of boundary nodes
              add_boundary_node(q1_left_boundary_id - 1, nod_pt);

              // Add the node into the vector of nodes on edge 0 in clockwise
              // order
              edge0_vector_pt[(nnode_edge - 1) - (j / nnode_1d)] = nod_pt;

              // Add it to the list of nodes in the mesh
              Node_pt.push_back(nod_pt);

              // Indicate that the j'th node has been set up
              done = true;
            }
          }

          // Node is not on a mesh boundary but on the left edge
          if (!done)
          {
            // Construct a normal node
            Node* nod_pt = el1_pt->construct_node(j, time_stepper_pt);

            // Add the node into the vector of nodes on edge 0 in clockwise
            // order
            edge0_vector_pt[(nnode_edge - 1) - (j / nnode_1d)] = nod_pt;

            // Add it to the list of nodes in the mesh
            Node_pt.push_back(nod_pt);

            // Indicate the j'th node has been set up
            done = true;
          }
        }

        // Node is not on a mesh boundary or the scaffold element boundary
        if (!done)
        {
          // Construct a normal node
          Node* nod_pt = el1_pt->construct_node(j, time_stepper_pt);

          // Add it to the list of nodes in the mesh
          Node_pt.push_back(nod_pt);
        }

        // Get local coordinate
        el1_pt->local_coordinate_of_node(j, s);

        // Interpolate position linearly between vertex nodes
        dummy_element.interpolated_x(s, x);
        el1_pt->node_pt(j)->x(0) = x[0];
        el1_pt->node_pt(j)->x(1) = x[1];
      }


      // We now need to loop over nodes corner_2 to (corner_3-1) and copy the
      // given information from Q0. We do not need to set up the (corner_3)'th
      // node since it coincides with the centroid which has already been set up
      for (unsigned j = corner_2; j < corner_3; j++)
      {
        // The nodes along the upper edge of Q1 go from corner_2 to corner_3-1
        // while the nodes along the right edge of Q0 go from corner_1 to
        // (corner_3-nnode_1d) in increments of nnode_1d
        el1_pt->node_pt(j) =
          el0_pt->node_pt(corner_1 + (j - corner_2) * nnode_1d);
      }


      // Set the nodal positions of the dummy element to now emulate the
      // THIRD quad element:
      //------------------------------------------------------------------
      // Note: we do not need to change the top-right corner since it always
      // coincides with the centroid of the triangle element

      // Bottom-left corner
      dummy_element.node_pt(0)->x(0) = triangle_vertex[2][0];
      dummy_element.node_pt(0)->x(1) = triangle_vertex[2][1];

      // Bottom-right corner
      dummy_element.node_pt(1)->x(0) =
        0.5 * (triangle_vertex[0][0] + triangle_vertex[2][0]);
      dummy_element.node_pt(1)->x(1) =
        0.5 * (triangle_vertex[0][1] + triangle_vertex[2][1]);

      // Top-left corner
      dummy_element.node_pt(2)->x(0) =
        0.5 * (triangle_vertex[1][0] + triangle_vertex[2][0]);
      dummy_element.node_pt(2)->x(1) =
        0.5 * (triangle_vertex[1][1] + triangle_vertex[2][1]);


      // Set up all of the nodes in the third quad element (Q2):
      //--------------------------------------------------------
      // Here we need to notice that we have already set up the final nnode_1d
      // nodes (the upper edge of Q2 coincides with the right edge of Q1).
      // We have also already set up the nodes on the right edge of Q2 (the
      // right edge of Q2 coincides with the upper edge of Q0)

      // Loop over nodes 1 to (corner_2-1)
      for (unsigned j = 1; j < corner_2; j++)
      {
        // Indicates whether or not the node has been set up yet
        bool done = false;

        // On the lower edge
        if (j < nnode_1d - 1)
        {
          // If the lower edge has already been set up then we already know the
          // nodes along this edge
          if (edge2_setup)
          {
            // The j'th node along this edge is the (nnode_1d-j)'th node in the
            // vector (remembering that the ordering is backwards if it has
            // already been set up)
            el2_pt->node_pt(j) = edge2_vector_pt[(nnode_edge - 1) - j];

            // Since the node has already been set up we do not need to sort
            // out its global coordinate data
            continue;
          }
          // If the edge hasn't been set up yet
          else
          {
            // If the node lies on a boundary too then we need to construct a
            // boundary node
            if (q2_lower_boundary_id > 0)
            {
              // Construct a boundary node
              Node* nod_pt =
                el2_pt->construct_boundary_node(j, time_stepper_pt);

              // Add it to the list of boundary nodes
              add_boundary_node(q2_lower_boundary_id - 1, nod_pt);

              // Add the node into the vector of nodes on edge 2
              edge2_vector_pt[j] = nod_pt;

              // Add it to the list of nodes in the mesh
              Node_pt.push_back(nod_pt);

              // Indicate the j'th node has been set up
              done = true;
            }
          }

          // Node is not on a mesh boundary but on the lower edge
          if (!done)
          {
            // Construct a normal node
            Node* nod_pt = el2_pt->construct_node(j, time_stepper_pt);

            // Add the node into the vector of nodes on edge 2
            edge2_vector_pt[j] = nod_pt;

            // Add it to the list of nodes in the mesh
            Node_pt.push_back(nod_pt);

            // Indicate the j'th node has been set up
            done = true;
          }
        }
        // On the right edge
        else if ((j + 1) % nnode_1d == 0)
        {
          // Copy the data from the top edge of element 0 to element 2
          el2_pt->node_pt(j) =
            el0_pt->node_pt((corner_2 - 1) + (j + 1) / nnode_1d);

          // We don't need to set up the global coordinate data so just
          // skip to the next node in the element
          continue;
        }
        // On the left edge
        else if (j % nnode_1d == 0)
        {
          // If the left edge has already been set up then we already know the
          // nodes along this edge
          if (edge1_setup)
          {
            // The j'th node along this edge is the (j/nnode_1d)'th node in the
            // vector
            el2_pt->node_pt(j) = edge1_vector_pt[j / nnode_1d];

            // Since the node has already been set up we do not need to sort
            // out its global coordinate data
            continue;
          }
          // If the edge hasn't been set up yet
          else
          {
            if (q2_left_boundary_id > 0)
            {
              // Construct a boundary node
              Node* nod_pt =
                el2_pt->construct_boundary_node(j, time_stepper_pt);

              // Add it to the list of boundary nodes
              add_boundary_node(q2_left_boundary_id - 1, nod_pt);

              // Add the node into the vector of nodes on edge 1 in clockwise
              // order
              edge1_vector_pt[(nnode_edge - 1) - (j / nnode_1d)] = nod_pt;

              // Add it to the list of nodes in the mesh
              Node_pt.push_back(nod_pt);

              // Indicate that the j'th node has been set up
              done = true;
            }
          }

          // Node is not on a mesh boundary but on the left edge
          if (!done)
          {
            // Construct a normal node
            Node* nod_pt = el2_pt->construct_node(j, time_stepper_pt);

            // Add the node into the vector of nodes on edge 1 in clockwise
            // order
            edge1_vector_pt[(nnode_edge - 1) - (j / nnode_1d)] = nod_pt;

            // Add it to the list of nodes in the mesh
            Node_pt.push_back(nod_pt);

            // Indicate the j'th node has been set up
            done = true;
          }
        }

        // Node is not on a mesh boundary
        if (!done)
        {
          // Construct a normal node
          Node* nod_pt = el2_pt->construct_node(j, time_stepper_pt);

          // Add it to the list of nodes in the mesh
          Node_pt.push_back(nod_pt);
        }

        // Get local coordinate
        el2_pt->local_coordinate_of_node(j, s);

        // Interpolate position linearly between vertex nodes
        dummy_element.interpolated_x(s, x);
        el2_pt->node_pt(j)->x(0) = x[0];
        el2_pt->node_pt(j)->x(1) = x[1];
      }

      // We now need to loop over nodes corner_2 to (corner_3-1) and copy the
      // given information from Q1. We do not need to set up the (corner_3)'th
      // node since it coincides with the centroid which has already been set up
      for (unsigned j = corner_2; j < corner_3; j++)
      {
        // The nodes along the upper edge of Q2 go from corner_2 to corner_3-1
        // while the nodes along the right edge of Q1 go from corner_1 to
        // (corner_3-nnode_1d) in increments of nnode_1d
        el2_pt->node_pt(j) =
          el1_pt->node_pt(corner_1 + (j - corner_2) * nnode_1d);
      }

      // Add the element pointers to Element_pt and then increment the counter
      Element_pt[new_el_count] = el0_pt;
      Element_pt[new_el_count + 1] = el1_pt;
      Element_pt[new_el_count + 2] = el2_pt;
      new_el_count += 3;

      // Assign the edges to the edge map
      edge_nodes_map[edge0] = edge0_vector_pt;
      edge_nodes_map[edge1] = edge1_vector_pt;
      edge_nodes_map[edge2] = edge2_vector_pt;
    }

    // Indicate that the look up scheme has been set up
    Lookup_for_elements_next_boundary_is_setup = true;

    // Clean the dummy element nodes
    for (unsigned j = 0; j < 4; j++)
    {
      // Delete the j-th dummy element node
      delete dummy_element.node_pt(j);

      // Make it a null pointer
      dummy_element.node_pt(j) = 0;
    }
  }


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Adapt problem based on specified elemental error estimates
  //======================================================================
  template<class ELEMENT>
  void RefineableQuadFromTriangleMesh<ELEMENT>::adapt(
    const Vector<double>& elem_error)
  {
    // Call the adapt function from the TreeBasedRefineableMeshBase base class
    TreeBasedRefineableMeshBase::adapt(elem_error);

#ifdef OOMPH_HAS_TRIANGLE_LIB
    // Move the nodes on the new boundary onto the old curvilinear
    // boundary. If the boundary is straight this will do precisely
    // nothing but will be somewhat inefficient
    this->snap_nodes_onto_geometric_objects();
#endif
  } // End of adapt
} // End of namespace oomph

#endif // OOMPH_QUAD_FROM_TRIANGLE_MESH_TEMPLATE_CC
