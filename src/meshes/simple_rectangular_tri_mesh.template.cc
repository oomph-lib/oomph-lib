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
#ifndef OOMPH_SIMPLE_RECTANGULAR_TRIMESH_TEMPLATE_CC
#define OOMPH_SIMPLE_RECTANGULAR_TRIMESH_TEMPLATE_CC


#include <algorithm>

#include "../generic/Telements.h"

// Simple 2D triangle mesh class
#include "simple_rectangular_tri_mesh.template.h"


namespace oomph
{
  //====================================================================
  /// Constructor for simple 2D triangular mesh class:
  ///
  /// n_x  : number of elements in the x direction
  ///
  /// n_y  : number of elements in the y direction
  ///
  /// l_x  : length in the x direction
  ///
  /// l_y  : length in the y direction
  ///
  /// Ordering of elements: 'lower left' to 'lower right' then 'upwards'
  //====================================================================
  template<class ELEMENT>
  SimpleRectangularTriMesh<ELEMENT>::SimpleRectangularTriMesh(
    const unsigned& n_x,
    const unsigned& n_y,
    const double& l_x,
    const double& l_y,
    TimeStepper* time_stepper_pt)
    : Nx(n_x), Ny(n_y), Lx(l_x), Ly(l_y)
  {
    using namespace MathematicalConstants;

    // Mesh can only be built with 2D Telements.
    MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(2);

    // Set number of boundaries
    this->set_nboundary(4);

    // Allocate the store for the elements
    unsigned n_element = (Nx) * (Ny)*2;
    Element_pt.resize(n_element, 0);

    // Create first element
    Element_pt[0] = new ELEMENT;

    // Currently this mesh only works for 3 and 6 noded triangles
    if ((finite_element_pt(0)->nnode_1d() != 2) &&
        (finite_element_pt(0)->nnode_1d() != 3) &&
        (finite_element_pt(0)->nnode_1d() != 4))
    {
      throw OomphLibError(
        "Currently this mesh only works for 3, 6 & 10-noded triangles",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    unsigned n_node = 0;
    // Unless nnode_1d returned as !=2 default no extra nodes
    unsigned additional_nnode = 0;

    // Allocate storage if no extra nodes
    if (finite_element_pt(0)->nnode_1d() == 2)
    {
      n_node = (Nx + 1) * (Ny + 1);
    }

    if (finite_element_pt(0)->nnode_1d() == 3)
    {
      additional_nnode = 3;
      // Allocate storage for nodes (including extra)
      n_node = (2 * Nx + 1) * (2 * Ny + 1);
    }

    if (finite_element_pt(0)->nnode_1d() == 4)
    {
      additional_nnode = 7;
      // Allocate storage for notes (including extra)
      n_node = (3 * Nx + 1) * (3 * Ny + 1);
    }

    Node_pt.resize(n_node);

    // Set up geometrical data
    //------------------------
    unsigned long node_count = 0;
    unsigned long element_count = 0;
    double xinit = 0.0, yinit = 0.0;
    // Set the values of the increments
    double xstep = Lx / (Nx);
    double ystep = Ly / (Ny);


    // FIRST NODE (lower left corner)
    //-------------------------------

    // Set the corner node

    // Allocate memory for the corner node of the first element
    //(which is not created loop as all subsequent bottom corners already exist)
    Node_pt[node_count] =
      finite_element_pt(0)->construct_node(1, time_stepper_pt);

    // Set the pointer from the element
    finite_element_pt(0)->node_pt(1) = Node_pt[node_count];

    // Don't increment node_count yet as we still need its containing box
    //(position and boundaries for first node are allocated in the main loop)

    // CREATE THE ELEMENTS (each has 3 local nodes)
    //--------------------------------------------------
    // Elements are created two at a time, the first being lower triangular
    // and the second upper triangular, so that they form a containing box.
    // Local nodes are numbered anti-clockwise with node_pt(1) being the
    // right-angle corner of the triangle.
    // The global node number, node_count, is considered to define
    // a containing box, with node_count defined as the node number
    // of the bottom left corner of the box.
    for (unsigned j = 0; j < Ny + 1; ++j)
    {
      for (unsigned i = 0; i < Nx + 1; ++i)
      {
        // CURRENT BOX
        //(nodes on RHS or top edge of domain do not define a box)
        if (j < Ny && i < Nx)
        {
          // Create two elements
          //------------------------------
          if (element_count != 0) // 0th element already exists
          {
            Element_pt[element_count] = new ELEMENT;
          }

          Element_pt[element_count + 1] = new ELEMENT;


          // Allocate memory for nodes in the current box
          //--------------------------------------------
          // If node on LHS of domain, allocate the top left corner node
          if (i == 0)
          {
            Node_pt[node_count + Nx + 1] =
              finite_element_pt(element_count)
                ->construct_node(0, time_stepper_pt);
          }

          // If on bottom row, allocate the bottom right corner node
          if (j == 0)
          {
            Node_pt[node_count + 1] = finite_element_pt(element_count)
                                        ->construct_node(2, time_stepper_pt);
          }

          // Always need to allocate top right corner node
          Node_pt[node_count + Nx + 2] = finite_element_pt(element_count + 1)
                                           ->construct_node(1, time_stepper_pt);

          // Bottom left corner node is already allocated


          // Set the pointers from the elements
          //----------------------------------
          finite_element_pt(element_count)->node_pt(0) =
            Node_pt[node_count + Nx + 1];
          finite_element_pt(element_count)->node_pt(1) = Node_pt[node_count];
          finite_element_pt(element_count)->node_pt(2) =
            Node_pt[node_count + 1];
          finite_element_pt(element_count + 1)->node_pt(0) =
            Node_pt[node_count + 1];
          finite_element_pt(element_count + 1)->node_pt(1) =
            Node_pt[node_count + Nx + 2];
          finite_element_pt(element_count + 1)->node_pt(2) =
            Node_pt[node_count + Nx + 1];

          element_count += 2;
        }

        // CURRENT (GLOBAL) NODE
        // Set the position
        Node_pt[node_count]->x(0) = xinit + i * xstep;
        Node_pt[node_count]->x(1) = yinit + j * ystep;

        // Add node to any relevant boundaries
        if (j == 0)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          add_boundary_node(0, Node_pt[node_count]);
        }
        if (j == Ny)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          add_boundary_node(2, Node_pt[node_count]);
        }
        if (i == 0)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          add_boundary_node(3, Node_pt[node_count]);
        }
        if (i == Nx)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          add_boundary_node(1, Node_pt[node_count]);
        }

        // Increment counter
        node_count++;
      }
    }


    if (additional_nnode == 3)
    {
      // Reset element counter to allow looping over existing elements
      // to add extra nodes.
      // Note that the node_count is not reset so no entries are overwritten
      element_count = 0;
      for (unsigned j = 0; j < Ny + 1; ++j)
      {
        // Note: i counter reduced by 1 since i axis runs through middle of
        // elements on x-axis
        for (unsigned i = 0; i < Nx; ++i)
        {
          // The additional nodes will be added in stages for ease of
          // application. Note: local numbering follows the anti-clockwise
          // pattern, node 3 halfway between  nodes 0-1, 4 between 1-2 and the
          // 5th local node between 2-0.

          // Restricted to stop creation of additional nodes outside the mesh
          if (j < Ny)
          {
            // If on the bottom row middle-node at bottom needs allocating
            if (j == 0)
            {
              Node_pt[node_count] = finite_element_pt(element_count)
                                      ->construct_node(4, time_stepper_pt);
            }

            // Due to directions of iteration node at middle of top box edge
            // (in next element) always needs allocating
            Node_pt[node_count + Nx] = finite_element_pt(element_count + 1)
                                         ->construct_node(4, time_stepper_pt);

            // Set pointers to the top/bottom middle nodes
            // Bottom node
            finite_element_pt(element_count)->node_pt(4) = Node_pt[node_count];
            // Top node
            finite_element_pt(element_count + 1)->node_pt(4) =
              Node_pt[node_count + Nx];

            // Increase the element counter to add top/bottom nodes to
            // next pair of element on next pass
            element_count += 2;
          } // End extra top/bottom node construction

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = xinit + double(i + 0.5) * xstep;
          Node_pt[node_count]->x(1) = yinit + j * ystep;

          // Add node to any applicable boundaries (node 4's can only be top
          // or bottom boundaries)
          if (j == 0)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(0, Node_pt[node_count]);
          }
          if (j == Ny)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(2, Node_pt[node_count]);
          }

          // Update node_count
          node_count++;
        }
      }

      // Next stage of additional node implementation involes the middle left
      // and right nodes (local number 3 on each triangle)

      // Element counter reset for second loop over existing elements
      element_count = 0;
      // Note: j counter reduced by 1 since j axis runs through middle of
      // elements on y-axis
      for (unsigned j = 0; j < Ny; ++j)
      {
        for (unsigned i = 0; i < Nx + 1; ++i)
        {
          if (j < Ny && i < Nx)
          {
            // If on the left hand boundary the node at middle of the left
            // side needs allocating
            if (i == 0)
            {
              Node_pt[node_count] = finite_element_pt(element_count)
                                      ->construct_node(3, time_stepper_pt);
            }

            // The mid node on the right hand side always needs constructing
            // within the bounds of the mesh
            Node_pt[node_count + 1] = finite_element_pt(element_count + 1)
                                        ->construct_node(3, time_stepper_pt);

            // Set pointers from the elements to new nodes
            finite_element_pt(element_count)->node_pt(3) = Node_pt[node_count];
            finite_element_pt(element_count + 1)->node_pt(3) =
              Node_pt[node_count + 1];

            // Increase element_count by 2
            element_count += 2;
          } // End extra left/right node construction

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = xinit + i * xstep;
          Node_pt[node_count]->x(1) = yinit + double(j + 0.5) * ystep;

          // Add node to any applicable boundaries again - only be left/right
          if (i == 0)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(3, Node_pt[node_count]);
          }
          if (i == Nx)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(1, Node_pt[node_count]);
          }

          // Update node_count
          node_count++;
        }
      }

      // Final stage of inserting extra nodes - inclusion of the local
      // number 5 (middle of hypot. edge)

      element_count = 0;
      // Note: both i,j counters reduced by 1 since j axis runs through middle
      // of elements in both x,y
      for (unsigned j = 0; j < Ny; ++j)
      {
        for (unsigned i = 0; i < Nx; ++i)
        {
          // The middle node always needs constructing
          Node_pt[node_count] = finite_element_pt(element_count)
                                  ->construct_node(5, time_stepper_pt);

          // Set pointers from the elements to new nodes
          finite_element_pt(element_count)->node_pt(5) = Node_pt[node_count];
          finite_element_pt(element_count + 1)->node_pt(5) =
            Node_pt[node_count];

          // Increase element_count by 2
          element_count += 2;
          // End extra left/right node construction

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = xinit + double(i + 0.5) * xstep;
          Node_pt[node_count]->x(1) = yinit + double(j + 0.5) * ystep;

          // All nodes are internal in this structure so no boundaries
          // applicable

          // Update node_count
          node_count++;
        }
      }
    }
    // End of extra nodes for 6 noded trianglur elements


    if (additional_nnode == 7)
    {
      // Reset element counter to allow looping over existing elements
      // to add extra nodes.
      // Note that the node_count is not reset so no entries are overwritten
      element_count = 0;
      for (unsigned j = 0; j < Ny + 1; ++j)
      {
        // Note: i counter reduced by 1 for implementation of lower element
        // node 5 and upper element node 6's.
        for (unsigned i = 0; i < Nx; ++i)
        {
          // The additional nodes will be added in stages for ease of
          // application. Note: local numbering follows the anti-clockwise
          // pattern, nodes 3 and 4 equidistant between nodes 0-1, 5 and 6
          // between 1-2, 7 and 8  between 2-0 and last node 9 located
          // (internally) at the centroid of the triangle.

          // Restricted to stop creation of additional nodes outside the mesh
          if (j < Ny)
          {
            // If on the bottom row middle-left-node at bottom needs allocating
            if (j == 0)
            {
              Node_pt[node_count] = finite_element_pt(element_count)
                                      ->construct_node(5, time_stepper_pt);
            }

            // Due to directions of iteration node at middle left of top box
            // edge (in next element) always needs allocating
            Node_pt[node_count + Nx] = finite_element_pt(element_count + 1)
                                         ->construct_node(6, time_stepper_pt);

            // Set pointers to the top/bottom middle nodes
            // Bottom node
            finite_element_pt(element_count)->node_pt(5) = Node_pt[node_count];
            // Top node
            finite_element_pt(element_count + 1)->node_pt(6) =
              Node_pt[node_count + Nx];

            // Increase the element counter to add top/bottom nodes to
            // next pair of element on next pass
            element_count += 2;
          } // End extra top/bottom node construction

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = xinit + double(i + 1.0 / 3.0) * xstep;
          Node_pt[node_count]->x(1) = yinit + j * ystep;

          // Add node to any applicable boundaries (exactly as previous)
          if (j == 0)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(0, Node_pt[node_count]);
          }
          if (j == Ny)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(2, Node_pt[node_count]);
          }

          // Update node_count
          node_count++;
        }
      }


      // Next stage: bottom-middle-right node (node 6 in lower tri.el.)
      // and top-middle-right node (node 5 in upper tri.el.)

      element_count = 0;
      for (unsigned j = 0; j < Ny + 1; ++j)
      {
        // Note: i counter as for above 5/6's
        for (unsigned i = 0; i < Nx; ++i)
        {
          // Restricted to stop creation of additional nodes outside the mesh
          if (j < Ny)
          {
            // If on the bottom row middle-right-node at bottom needs allocating
            if (j == 0)
            {
              Node_pt[node_count] = finite_element_pt(element_count)
                                      ->construct_node(6, time_stepper_pt);
            }

            // Due to directions of iteration node at middle left of top box
            // edge (in next element) always needs allocating
            Node_pt[node_count + Nx] = finite_element_pt(element_count + 1)
                                         ->construct_node(5, time_stepper_pt);

            // Set pointers to the top/bottom middle nodes
            // Bottom node
            finite_element_pt(element_count)->node_pt(6) = Node_pt[node_count];
            // Top node
            finite_element_pt(element_count + 1)->node_pt(5) =
              Node_pt[node_count + Nx];

            // Increase the element counter to add top/bottom nodes to
            // next pair of element on next pass
            element_count += 2;
          } // End extra top/bottom node construction

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = xinit + double(i + 2.0 / 3.0) * xstep;
          Node_pt[node_count]->x(1) = yinit + j * ystep;

          // Add node to any applicable boundaries (exactly as previous)
          if (j == 0)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(0, Node_pt[node_count]);
          }
          if (j == Ny)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(2, Node_pt[node_count]);
          }

          // Update node_count
          node_count++;
        }
      }


      // Next stage of additional node implementation involes node 4 on lower
      // tri. el. and node 3 on upper tri. el.

      // Element counter reset for next loop over existing elements
      element_count = 0;
      // Note: j counter reduced by 1 similar to others above
      for (unsigned j = 0; j < Ny; ++j)
      {
        for (unsigned i = 0; i < Nx + 1; ++i)
        {
          if (j < Ny && i < Nx)
          {
            // If on the left hand boundary the lower middle node of the left
            // side needs allocating
            if (i == 0)
            {
              Node_pt[node_count] = finite_element_pt(element_count)
                                      ->construct_node(4, time_stepper_pt);
            }

            // The mid node on the right hand side always needs constructing
            // within the bounds of the mesh
            Node_pt[node_count + 1] = finite_element_pt(element_count + 1)
                                        ->construct_node(3, time_stepper_pt);

            // Set pointers from the elements to new nodes
            finite_element_pt(element_count)->node_pt(4) = Node_pt[node_count];
            finite_element_pt(element_count + 1)->node_pt(3) =
              Node_pt[node_count + 1];

            // Increase element_count by 2
            element_count += 2;
          } // End extra left/right node construction

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = xinit + i * xstep;
          Node_pt[node_count]->x(1) = yinit + double(j + 1.0 / 3.0) * ystep;

          // Add node to any applicable boundaries again - only be left/right
          if (i == 0)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(3, Node_pt[node_count]);
          }
          if (i == Nx)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(1, Node_pt[node_count]);
          }

          // Update node_count
          node_count++;
        }
      }


      // Next stage of additional node implementation involes node 3 on lower
      // tri. el. and node 4 on upper tri. el.

      // Element counter reset for next loop over existing elements
      element_count = 0;
      // Note: j counter reduced by 1 similar to others above
      for (unsigned j = 0; j < Ny; ++j)
      {
        for (unsigned i = 0; i < Nx + 1; ++i)
        {
          if (j < Ny && i < Nx)
          {
            // If on the left hand boundary the lower middle node of the left
            // side needs allocating
            if (i == 0)
            {
              Node_pt[node_count] = finite_element_pt(element_count)
                                      ->construct_node(3, time_stepper_pt);
            }

            // The mid node on the right hand side always needs constructing
            // within the bounds of the mesh
            Node_pt[node_count + 1] = finite_element_pt(element_count + 1)
                                        ->construct_node(4, time_stepper_pt);

            // Set pointers from the elements to new nodes
            finite_element_pt(element_count)->node_pt(3) = Node_pt[node_count];
            finite_element_pt(element_count + 1)->node_pt(4) =
              Node_pt[node_count + 1];

            // Increase element_count by 2
            element_count += 2;
          } // End extra left/right node construction

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = xinit + i * xstep;
          Node_pt[node_count]->x(1) = yinit + double(j + 2.0 / 3.0) * ystep;

          // Add node to any applicable boundaries again - only be left/right
          if (i == 0)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(3, Node_pt[node_count]);
          }
          if (i == Nx)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(1, Node_pt[node_count]);
          }

          // Update node_count
          node_count++;
        }
      }


      // Next is the lower tri. el. totally internal node with local number 9
      element_count = 0;
      // Note: both i,j counters reduced by 1
      for (unsigned j = 0; j < Ny; ++j)
      {
        for (unsigned i = 0; i < Nx; ++i)
        {
          // The middle node always needs constructing
          Node_pt[node_count] = finite_element_pt(element_count)
                                  ->construct_node(9, time_stepper_pt);

          // Set pointers from the elements to new nodes
          finite_element_pt(element_count)->node_pt(9) = Node_pt[node_count];

          // Increase element_count by 2
          element_count += 2;

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = xinit + double(i + 1.0 / 3.0) * xstep;
          Node_pt[node_count]->x(1) = yinit + double(j + 1.0 / 3.0) * ystep;

          // All nodes are internal in this structure so no boundaries
          // applicable

          // Update node_count
          node_count++;
        }
      }


      // Next is the bottom hyp. node -
      // lower tri. el. number 7, upper tri. el. number 8
      element_count = 0;
      // Note: both i,j counters reduced by 1
      for (unsigned j = 0; j < Ny; ++j)
      {
        for (unsigned i = 0; i < Nx; ++i)
        {
          // The node always needs constructing
          Node_pt[node_count] = finite_element_pt(element_count)
                                  ->construct_node(7, time_stepper_pt);

          // Set pointers from the elements to new nodes
          finite_element_pt(element_count)->node_pt(7) = Node_pt[node_count];
          finite_element_pt(element_count + 1)->node_pt(8) =
            Node_pt[node_count];

          // Increase element_count by 2
          element_count += 2;

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = xinit + double(i + 2.0 / 3.0) * xstep;
          Node_pt[node_count]->x(1) = yinit + double(j + 1.0 / 3.0) * ystep;

          // All nodes are internal in this structure so no boundaries
          // applicable

          // Update node_count
          node_count++;
        }
      }


      // Next is the upper hyp. node -
      // lower tri. el. number 8, upper tri. el. number 7
      element_count = 0;
      // Note: both i,j counters reduced by 1
      for (unsigned j = 0; j < Ny; ++j)
      {
        for (unsigned i = 0; i < Nx; ++i)
        {
          // The node always needs constructing
          Node_pt[node_count] = finite_element_pt(element_count)
                                  ->construct_node(8, time_stepper_pt);

          // Set pointers from the elements to new nodes
          finite_element_pt(element_count)->node_pt(8) = Node_pt[node_count];
          finite_element_pt(element_count + 1)->node_pt(7) =
            Node_pt[node_count];

          // Increase element_count by 2
          element_count += 2;

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = xinit + double(i + 1.0 / 3.0) * xstep;
          Node_pt[node_count]->x(1) = yinit + double(j + 2.0 / 3.0) * ystep;

          // All nodes are internal in this structure so no boundaries
          // applicable

          // Update node_count
          node_count++;
        }
      }


      // Next is the upper tri. el. totally internal node with local number 9
      element_count = 0;
      // Note: both i,j counters reduced by 1
      for (unsigned j = 0; j < Ny; ++j)
      {
        for (unsigned i = 0; i < Nx; ++i)
        {
          // The middle node always needs constructing
          Node_pt[node_count] = finite_element_pt(element_count + 1)
                                  ->construct_node(9, time_stepper_pt);

          // Set pointers from the elements to new nodes
          finite_element_pt(element_count + 1)->node_pt(9) =
            Node_pt[node_count];

          // Increase element_count by 2
          element_count += 2;

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = xinit + double(i + 2.0 / 3.0) * xstep;
          Node_pt[node_count]->x(1) = yinit + double(j + 2.0 / 3.0) * ystep;

          // All nodes are internal in this structure so no boundaries
          // applicable

          // Update node_count
          node_count++;
        }
      }
    }
  }

} // namespace oomph
#endif
