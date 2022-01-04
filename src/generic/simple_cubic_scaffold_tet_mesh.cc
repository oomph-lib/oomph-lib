// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
#include "Telements.h"
#include "simple_cubic_scaffold_tet_mesh.h"


namespace oomph
{
  //===========================================================
  /// Scaffold mesh for cubic tet mesh.
  //===========================================================
  SimpleCubicScaffoldTetMesh::SimpleCubicScaffoldTetMesh(
    const unsigned& n_x,
    const unsigned& n_y,
    const unsigned& n_z,
    const double& l_x,
    const double& l_y,
    const double& l_z,
    TimeStepper* time_stepper_pt)
  {
    // Storage for the pointers to the nodes
    RankThreeTensor<Node*> vertex_node_pt(n_x + 1, n_y + 1, n_z + 1);
    RankThreeTensor<Node*> left_face_node_pt(n_x + 1, n_y + 1, n_z + 1);
    RankThreeTensor<Node*> front_face_node_pt(n_x + 1, n_y + 1, n_z + 1);
    RankThreeTensor<Node*> down_face_node_pt(n_x + 1, n_y + 1, n_z + 1);
    RankThreeTensor<Node*> central_node_pt(n_x, n_y, n_z);

    // "Lower left" corner
    double x_0 = 0.0;
    double y_0 = 0.0;
    double z_0 = 0.0;

    // Increments
    double dx = l_x / double(n_x);
    double dy = l_y / double(n_y);
    double dz = l_z / double(n_z);

    // Total number of nodes
    unsigned nnod = (n_x + 1) * (n_y + 1) * (n_z + 1) + // vertex nodes
                    n_x * n_y * n_z + // central nodes
                    (n_x * n_y) * (n_z + 1) + // faces
                    (n_x * n_z) * (n_y + 1) + (n_y * n_z) * (n_x + 1);
    Node_pt.reserve(nnod);

    // Total number of elements
    unsigned nelem = 24 * n_x * n_y * n_z;
    Element_pt.reserve(nelem);

    // Six boundaries
    set_nboundary(6);


    // Generate the vertex nodes of all cells
    //=======================================
    Node* node_pt;
    for (unsigned k = 0; k < n_z + 1; k++)
    {
      for (unsigned j = 0; j < n_y + 1; j++)
      {
        for (unsigned i = 0; i < n_x + 1; i++)
        {
          // Need to work out whether to create a boundary node or not
          if (((i < n_x + 1) && (j < n_y + 1) && (k < n_z + 1)) &&
              ((i == 0) || (i == n_x) || (j == 0) || (j == n_y) || (k == 0) ||
               (k == n_z)))
          {
            node_pt = new BoundaryNode<Node>(3, 1, 0);
          }
          else
          {
            node_pt = new Node(3, 1, 0);
          }

          vertex_node_pt(i, j, k) = node_pt;

          // Coordinates
          node_pt->x(0) = x_0 + double(i) * dx;
          node_pt->x(1) = y_0 + double(j) * dy;
          node_pt->x(2) = z_0 + double(k) * dz;

          // Real node?
          if ((i < n_x + 1) && (j < n_y + 1) && (k < n_z + 1))
          {
            // Add node
            Node_pt.push_back(node_pt);

            // Left boundary is 0
            if (i == 0)
            {
              add_boundary_node(0, node_pt);
            }
            // Right boundary is 1
            else if (i == n_x)
            {
              add_boundary_node(1, node_pt);
            }

            // Front boundary is 2
            if (j == 0)
            {
              add_boundary_node(2, node_pt);
            }
            // Back boundary is 3
            else if (j == n_y)
            {
              add_boundary_node(3, node_pt);
            }

            // Down boundary is 4
            if (k == 0)
            {
              add_boundary_node(4, node_pt);
            }
            // Up boundary is 5
            else if (k == n_z)
            {
              add_boundary_node(5, node_pt);
            }
          }
        }
      }
    }


    // Generate the face nodes of all cells
    //=====================================
    for (unsigned k = 0; k < n_z + 1; k++)
    {
      for (unsigned j = 0; j < n_y + 1; j++)
      {
        for (unsigned i = 0; i < n_x + 1; i++)
        {
          // Node on front face of cell
          //---------------------------
          // Need to work out whether to create boundary node or not
          if (((i < n_x) && (k < n_z)) && ((j == 0) || (j == n_y)))
          {
            node_pt = new BoundaryNode<Node>(3, 1, 0);
          }
          else
          {
            node_pt = new Node(3, 1, 0);
          }
          front_face_node_pt(i, j, k) = node_pt;


          // Coordinates
          node_pt->x(0) = x_0 + 0.5 * dx + double(i) * dx;
          node_pt->x(1) = y_0 + double(j) * dy;
          node_pt->x(2) = z_0 + 0.5 * dz + double(k) * dz;

          // Real node?
          if ((i < n_x) && (k < n_z))
          {
            // Add node to mesh
            Node_pt.push_back(node_pt);

            // Front boundary is 2
            if (j == 0)
            {
              add_boundary_node(2, node_pt);
            }
            // Back boundary is 3
            else if (j == n_y)
            {
              add_boundary_node(3, node_pt);
            }
          }

          // Node on left face of cell
          //--------------------------
          // Work out whether to construct boundary node or not
          if (((j < n_y) && (k < n_z)) && ((i == 0) || (i == n_x)))
          {
            node_pt = new BoundaryNode<Node>(3, 1, 0);
          }
          else
          {
            node_pt = new Node(3, 1, 0);
          }
          left_face_node_pt(i, j, k) = node_pt;

          // Coordinates
          node_pt->x(0) = x_0 + double(i) * dx;
          node_pt->x(1) = y_0 + 0.5 * dy + double(j) * dy;
          node_pt->x(2) = z_0 + 0.5 * dz + double(k) * dz;


          // Real node?
          if ((j < n_y) && (k < n_z))
          {
            // Add node to mesh
            Node_pt.push_back(node_pt);

            // Left boundary is 0
            if (i == 0)
            {
              add_boundary_node(0, node_pt);
            }
            // Right boundary is 1
            else if (i == n_x)
            {
              add_boundary_node(1, node_pt);
            }
          }

          // Node on down face of cell
          //--------------------------
          if (((i < n_x) && (j < n_y)) && ((k == 0) || (k == n_z)))
          {
            node_pt = new BoundaryNode<Node>(3, 1, 0);
          }
          else
          {
            node_pt = new Node(3, 1, 0);
          }
          down_face_node_pt(i, j, k) = node_pt;

          // Coordinates
          node_pt->x(0) = x_0 + 0.5 * dx + double(i) * dx;
          node_pt->x(1) = y_0 + 0.5 * dy + double(j) * dy;
          node_pt->x(2) = z_0 + double(k) * dz;


          // Real node?
          if ((i < n_x) && (j < n_y))
          {
            // Add node to mesh
            Node_pt.push_back(node_pt);

            // Down boundary is 4
            if (k == 0)
            {
              add_boundary_node(4, node_pt);
            }
            // Up boundary is 5
            else if (k == n_z)
            {
              add_boundary_node(5, node_pt);
            }
          }
        }
      }
    }


    // Central nodes for all cells
    for (unsigned k = 0; k < n_z; k++)
    {
      for (unsigned j = 0; j < n_y; j++)
      {
        for (unsigned i = 0; i < n_x; i++)
        {
          node_pt = new Node(3, 1, 0);
          central_node_pt(i, j, k) = node_pt;
          Node_pt.push_back(node_pt);
          node_pt->x(0) = x_0 + 0.5 * dx + double(i) * dx;
          node_pt->x(1) = y_0 + 0.5 * dy + double(j) * dy;
          node_pt->x(2) = z_0 + 0.5 * dz + double(k) * dz;
        }
      }
    }


    // Loop over blocks and create elements
    TElement<3, 2>* el_pt = 0;
    for (unsigned k = 0; k < n_z; k++)
    {
      for (unsigned j = 0; j < n_y; j++)
      {
        for (unsigned i = 0; i < n_x; i++)
        {
          // FRONT FACE
          //===========

          // Left element on front face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // LFD
          el_pt->node_pt(0) = vertex_node_pt(i, j, k);

          // Central front face
          el_pt->node_pt(1) = front_face_node_pt(i, j, k);

          // LFU
          el_pt->node_pt(2) = vertex_node_pt(i, j, k + 1);

          // central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Down element on front face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // LFD
          el_pt->node_pt(0) = vertex_node_pt(i, j, k);

          // RFD
          el_pt->node_pt(1) = vertex_node_pt(i + 1, j, k);

          // Central front face
          el_pt->node_pt(2) = front_face_node_pt(i, j, k);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Right element on front face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RFD
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j, k);

          // RFU
          el_pt->node_pt(1) = vertex_node_pt(i + 1, j, k + 1);

          // Central front face
          el_pt->node_pt(2) = front_face_node_pt(i, j, k);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Up element on front face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RFU
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j, k + 1);

          // LFU
          el_pt->node_pt(1) = vertex_node_pt(i, j, k + 1);

          // Central front face
          el_pt->node_pt(2) = front_face_node_pt(i, j, k);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // RIGHT FACE
          //===========

          // Front element on right face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RFD
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j, k);

          // Central right face
          el_pt->node_pt(1) = left_face_node_pt(i + 1, j, k);

          // RFU
          el_pt->node_pt(2) = vertex_node_pt(i + 1, j, k + 1);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Down element on right face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RFD
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j, k);

          // RBD
          el_pt->node_pt(1) = vertex_node_pt(i + 1, j + 1, k);

          // Central right face
          el_pt->node_pt(2) = left_face_node_pt(i + 1, j, k);

          // central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Back element on right face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RBD
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j + 1, k);

          // RBU
          el_pt->node_pt(1) = vertex_node_pt(i + 1, j + 1, k + 1);

          // Central right face
          el_pt->node_pt(2) = left_face_node_pt(i + 1, j, k);

          // central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Up element on right face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RBU
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j + 1, k + 1);

          // RFU
          el_pt->node_pt(1) = vertex_node_pt(i + 1, j, k + 1);

          // Central right face
          el_pt->node_pt(2) = left_face_node_pt(i + 1, j, k);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // UP FACE
          //===========

          // Front element on up face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // LFU
          el_pt->node_pt(0) = vertex_node_pt(i, j, k + 1);

          // RFU
          el_pt->node_pt(1) = vertex_node_pt(i + 1, j, k + 1);

          // Central up face
          el_pt->node_pt(2) = down_face_node_pt(i, j, k + 1);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Front element on up face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RFU
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j, k + 1);

          // RBU
          el_pt->node_pt(1) = vertex_node_pt(i + 1, j + 1, k + 1);

          // Central up face
          el_pt->node_pt(2) = down_face_node_pt(i, j, k + 1);

          // central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Back element on up face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RBU
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j + 1, k + 1);

          // LBU
          el_pt->node_pt(1) = vertex_node_pt(i, j + 1, k + 1);

          // Central up face
          el_pt->node_pt(2) = down_face_node_pt(i, j, k + 1);

          // central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Left element on up face
          //---------------------------
          el_pt = new TElement<3, 2>;


          // LBU
          el_pt->node_pt(0) = vertex_node_pt(i, j + 1, k + 1);

          // RFU
          el_pt->node_pt(1) = vertex_node_pt(i, j, k + 1);

          // Central up face
          el_pt->node_pt(2) = down_face_node_pt(i, j, k + 1);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // DOWN FACE
          //===========

          // Front element on down face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // LFD
          el_pt->node_pt(0) = vertex_node_pt(i, j, k);

          // RFD
          el_pt->node_pt(2) = vertex_node_pt(i + 1, j, k);

          // Central down face
          el_pt->node_pt(1) = down_face_node_pt(i, j, k);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Front element on down face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RFD
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j, k);

          // RBD
          el_pt->node_pt(2) = vertex_node_pt(i + 1, j + 1, k);

          // Central down face
          el_pt->node_pt(1) = down_face_node_pt(i, j, k);

          // central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Back element on down face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RBD
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j + 1, k);

          // LBD
          el_pt->node_pt(2) = vertex_node_pt(i, j + 1, k);

          // Central down face
          el_pt->node_pt(1) = down_face_node_pt(i, j, k);

          // central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Left element on down face
          //---------------------------
          el_pt = new TElement<3, 2>;


          // LBD
          el_pt->node_pt(0) = vertex_node_pt(i, j + 1, k);

          // RFD
          el_pt->node_pt(2) = vertex_node_pt(i, j, k);

          // Central down face
          el_pt->node_pt(1) = down_face_node_pt(i, j, k);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // BACK FACE
          //===========

          // Left element on back face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // LBD
          el_pt->node_pt(0) = vertex_node_pt(i, j + 1, k);

          // Central back face
          el_pt->node_pt(2) = front_face_node_pt(i, j + 1, k);

          // LBU
          el_pt->node_pt(1) = vertex_node_pt(i, j + 1, k + 1);

          // central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Down element on back face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // LBD
          el_pt->node_pt(0) = vertex_node_pt(i, j + 1, k);

          // RBD
          el_pt->node_pt(2) = vertex_node_pt(i + 1, j + 1, k);

          // Central back face
          el_pt->node_pt(1) = front_face_node_pt(i, j + 1, k);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Right element on back face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RBD
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j + 1, k);

          // RBU
          el_pt->node_pt(2) = vertex_node_pt(i + 1, j + 1, k + 1);

          // Central back face
          el_pt->node_pt(1) = front_face_node_pt(i, j + 1, k);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Up element on back face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // RBU
          el_pt->node_pt(0) = vertex_node_pt(i + 1, j + 1, k + 1);

          // LBU
          el_pt->node_pt(2) = vertex_node_pt(i, j + 1, k + 1);

          // Central back face
          el_pt->node_pt(1) = front_face_node_pt(i, j + 1, k);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // LEFT FACE
          //===========

          // Front element on left face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // LFD
          el_pt->node_pt(0) = vertex_node_pt(i, j, k);

          // Central left face
          el_pt->node_pt(2) = left_face_node_pt(i, j, k);

          // LFU
          el_pt->node_pt(1) = vertex_node_pt(i, j, k + 1);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Down element on left face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // LFD
          el_pt->node_pt(0) = vertex_node_pt(i, j, k);

          // LBD
          el_pt->node_pt(2) = vertex_node_pt(i, j + 1, k);

          // Central left face
          el_pt->node_pt(1) = left_face_node_pt(i, j, k);

          // central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Back element on left face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // LBD
          el_pt->node_pt(0) = vertex_node_pt(i, j + 1, k);

          // LBU
          el_pt->node_pt(2) = vertex_node_pt(i, j + 1, k + 1);

          // Central left face
          el_pt->node_pt(1) = left_face_node_pt(i, j, k);

          // central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);


          // Up element on left face
          //---------------------------
          el_pt = new TElement<3, 2>;

          // LBU
          el_pt->node_pt(0) = vertex_node_pt(i, j + 1, k + 1);

          // LFU
          el_pt->node_pt(2) = vertex_node_pt(i, j, k + 1);

          // Central left face
          el_pt->node_pt(1) = left_face_node_pt(i, j, k);

          // Central node
          el_pt->node_pt(3) = central_node_pt(i, j, k);

          Element_pt.push_back(el_pt);
        }
      }
    }

    if (Node_pt.size() != nnod)
    {
      std::ostringstream error_stream;
      error_stream << "Some internal error in the constructor\n"
                   << "Actual number of nodes          : " << Node_pt.size()
                   << "\ndoesn't match the prediction    : " << nnod
                   << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    if (Element_pt.size() != nelem)
    {
      std::ostringstream error_stream;
      error_stream << "Some internal error in the constructor\n"
                   << "Actual number of elements       : " << Element_pt.size()
                   << "\ndoesn't match the prediction    : " << nelem
                   << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    //  // Deliberately switch nodes to invert elements -- test for self_test()
    //  nelem=nelement();
    //  for (unsigned e=0;e<nelem;e++)
    //   {
    //    Node* node_pt=finite_element_pt(e)->node_pt(0);
    //    finite_element_pt(e)->node_pt(0)=
    //     finite_element_pt(e)->node_pt(1);
    //    finite_element_pt(e)->node_pt(1)=node_pt;
    //   }
  }


} // namespace oomph
