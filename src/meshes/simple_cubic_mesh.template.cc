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
#ifndef OOMPH_SIMPLE_CUBIC_MESH_TEMPLATE_CC
#define OOMPH_SIMPLE_CUBIC_MESH_TEMPLATE_CC

// OOMPH-LIB Header
#include "simple_cubic_mesh.template.h"

namespace oomph
{
  //=========================================================================
  /// Generic mesh construction. This function contains all the details of
  /// the mesh generation process, including all the tedious loops, counting
  /// spacing and boundary functions.
  //========================================================================
  template<class ELEMENT>
  void SimpleCubicMesh<ELEMENT>::build_mesh(TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with 3D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

    if ((Nx == 1) || (Ny == 1) || (Nz == 1))
    {
      std::ostringstream error_message;
      error_message << "SimpleCubicMesh needs at least two elements in each,\n"
                    << "coordinate direction. You have specified \n"
                    << "Nx=" << Nx << "; Ny=" << Ny << "; Nz=" << Nz
                    << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Set the number of boundaries
    set_nboundary(6);

    // Allocate the store for the elements
    Element_pt.resize(Nx * Ny * Nz);
    // Create first element
    unsigned element_num = 0;
    Element_pt[element_num] = new ELEMENT;

    // Read out the number of linear points in the element
    unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

    // Can now allocate the store for the nodes
    Node_pt.resize((1 + (n_p - 1) * Nx) * (1 + (n_p - 1) * Ny) *
                   (1 + (n_p - 1) * Nz));

    // Set up geometrical data
    //------------------------

    unsigned long node_count = 0;

    // Set the length of the elments
    double el_length[3] = {(Xmax - Xmin) / double(Nx),
                           (Ymax - Ymin) / double(Ny),
                           (Zmax - Zmin) / double(Nz)};

    // Storage for local coordinate in element
    Vector<double> s_fraction;

    // Now assign the topology
    // Boundaries are numbered:
    //  0 is at the bottom
    // 1 2 3 4 from the front  proceeding anticlockwise
    // 5 is at the top
    // Pinned value are denoted by an integer value 1
    // Thus if a node is on two boundaries, ORing the values of the
    // boundary conditions will give the most restrictive case (pinning)

    // FIRST ELEMENT (lower left corner at z = 0 plane )
    //----------------------------------

    // Set the corner node
    // Storage for the local node number
    unsigned local_node_num = 0;
    // Allocate memory for the node
    Node_pt[node_count] =
      finite_element_pt(element_num)
        ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
      Node_pt[node_count];

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmin;

    // Add the node to boundaries 0,1 and 4
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(4, Node_pt[node_count]);
    // Increment the node number
    node_count++;

    // Loop over the other nodes in the first row in the x direction in
    // the z=0 plane
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // Set the local node number
      local_node_num = l2;

      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the local fraction of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundary 0 and 1
      add_boundary_node(0, Node_pt[node_count]);
      add_boundary_node(1, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // Loop over the other node columns in the z = 0 plane
    for (unsigned l1 = 1; l1 < n_p; l1++)
    {
      // Set the local node number
      local_node_num = l1 * n_p;

      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the first node of the row
      //(with boundaries with 0 and 4)
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundaries 0 and 4
      add_boundary_node(4, Node_pt[node_count]);
      add_boundary_node(0, Node_pt[node_count]);
      // Increment the node number
      node_count++;

      // Loop over the other nodes in the row
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Set the local node number
        local_node_num = l1 * n_p + l2;

        // Allocate the memory for the node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the boundary 0
        add_boundary_node(0, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }
    }


    //---------------------------------------------------------------------

    // Loop over the other node columns in the z direction
    for (unsigned l3 = 1; l3 < n_p; l3++)
    {
      // Set the local node number
      local_node_num = n_p * n_p * l3;

      // Allocate memory for the node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      // Add the node to boundaries 1 and 4
      add_boundary_node(1, Node_pt[node_count]);
      add_boundary_node(4, Node_pt[node_count]);
      // Increment the node number
      node_count++;

      // Loop over the other nodes in the first row in the x direction
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Set the local node number
        local_node_num = l2 + n_p * n_p * l3;

        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // Add the node to the boundary 1
        add_boundary_node(1, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // Loop over the other node columns
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // Set the local node number
        local_node_num = l1 * n_p + n_p * n_p * l3;

        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the first node of the row (with boundary 4)
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // Add the node to the boundary 4
        add_boundary_node(4, Node_pt[node_count]);
        // Increment the node number
        node_count++;

        // Loop over the other nodes in the row
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Set the local node number
          local_node_num = l2 + l1 * n_p + n_p * n_p * l3;

          // Allocate the memory for the node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // No boundary

          // Increment the node number
          node_count++;
        }
      }
    }


    //----------------------------------------------------------------------

    // CENTRE OF FIRST ROW OF ELEMENTS (PLANE Z = 0)
    //--------------------------------

    // Now loop over the first row of elements, apart from final element
    for (unsigned j = 1; j < (Nx - 1); j++)
    {
      // Allocate memory for new element
      element_num = j;
      Element_pt[element_num] = new ELEMENT;

      // We will begin with all the nodes at the plane z = 0

      // Do first row of nodes

      // First column of nodes is same as neighbouring element
      finite_element_pt(element_num)->node_pt(0) =
        finite_element_pt(element_num - 1)->node_pt((n_p - 1));

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Set the local node number
        local_node_num = l2;

        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the boundaries 0 an 1
        add_boundary_node(0, Node_pt[node_count]);
        add_boundary_node(1, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // Do the rest of the nodes at the plane z = 0
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // First column of nodes is same as neighbouring element
        finite_element_pt(element_num)->node_pt(l1 * n_p) =
          finite_element_pt(element_num - 1)->node_pt(l1 * n_p + (n_p - 1));

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Set the local node number
          local_node_num = l2 + l1 * n_p;

          // Allocate memory for the nodes
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) = Zmin;

          // Add the node to the Boundary
          add_boundary_node(0, Node_pt[node_count]);
          // Increment the node number
          node_count++;
        }
      }

      // Loop over the other node columns in the z direction
      for (unsigned l3 = 1; l3 < n_p; l3++)
      {
        // First column of nodes is same as neighbouring element
        finite_element_pt(j)->node_pt(l3 * n_p * n_p) =
          finite_element_pt(j - 1)->node_pt(l3 * n_p * n_p + (n_p - 1));

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Set the local node number
          local_node_num = l2 + l3 * n_p * n_p;

          // Allocate memory for the nodes
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin;
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // Add the node to the boundary 1
          add_boundary_node(1, Node_pt[node_count]);
          // Increment the node number
          node_count++;
        }

        // Do the rest of the nodes
        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // First column of nodes is same as neighbouring element
          finite_element_pt(j)->node_pt(l1 * n_p + l3 * n_p * n_p) =
            finite_element_pt(j - 1)->node_pt(l1 * n_p + (n_p - 1) +
                                              l3 * n_p * n_p);

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Set the local node number
            local_node_num = l2 + l1 * n_p + l3 * n_p * n_p;

            // Allocate memory for the nodes
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer from the element
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
            Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

            // No boundaries

            // Increment the node number
            node_count++;
          }
        }
      }
    }

    // MY FINAL ELEMENT IN FIRST ROW (lower right corner)
    //-----------------------------------------------

    // Allocate memory for new element
    element_num = Nx - 1;
    Element_pt[element_num] = new ELEMENT;

    // We will begin with all the nodes at the plane z = 0

    // Do first row of nodes

    // First node is same as neighbouring element
    finite_element_pt(element_num)->node_pt(0) =
      finite_element_pt(element_num - 1)->node_pt((n_p - 1));

    // New nodes for other columns
    for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
    {
      // Set the local node number
      local_node_num = l2;

      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
        Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundaries 0 an 1
      add_boundary_node(0, Node_pt[node_count]);
      add_boundary_node(1, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // Last node (corner)
    // Set the local node number
    local_node_num = n_p - 1;

    // Allocate memory for the node
    Node_pt[node_count] =
      finite_element_pt(element_num)
        ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
      Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmin;

    // Add the node to the boundaries
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(2, Node_pt[node_count]);
    // Increment the node number
    node_count++;

    // Do the middle nodes at the plane z = 0
    for (unsigned l1 = 1; l1 < n_p; l1++)
    {
      // First column of nodes is same as neighbouring element
      finite_element_pt(element_num)->node_pt(l1 * n_p) =
        finite_element_pt(element_num - 1)->node_pt(l1 * n_p + (n_p - 1));

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
      {
        // Set the local node number
        local_node_num = l2 + l1 * n_p;

        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the boundary 0
        add_boundary_node(0, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // New node for final column
      // Set the local node number
      local_node_num = l1 * n_p + (n_p - 1);

      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundaries 0 and 2
      add_boundary_node(2, Node_pt[node_count]);
      add_boundary_node(0, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // Loop over the other node columns in the z direction
    for (unsigned l3 = 1; l3 < n_p; l3++)
    {
      // First column of nodes is same as neighbouring element
      finite_element_pt(element_num)->node_pt(l3 * n_p * n_p) =
        finite_element_pt(element_num - 1)->node_pt(l3 * n_p * n_p + (n_p - 1));

      // New nodes for other rows (y=0)
      for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
      {
        // Set the local node number
        local_node_num = l2 + l3 * n_p * n_p;

        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // Add the node to the boundary 1
        add_boundary_node(1, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // Last node of the row (y=0)
      // Set the local node number
      local_node_num = n_p - 1 + l3 * n_p * n_p;

      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      // Add the node to the boundary 1 and 2
      add_boundary_node(1, Node_pt[node_count]);
      add_boundary_node(2, Node_pt[node_count]);
      // Increment the node number
      node_count++;

      // Do the rest of the nodes
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // First column of nodes is same as neighbouring element
        finite_element_pt(element_num)->node_pt(l1 * n_p + l3 * n_p * n_p) =
          finite_element_pt(element_num - 1)
            ->node_pt(l1 * n_p + (n_p - 1) + l3 * n_p * n_p);

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
        {
          // Set the local node number
          local_node_num = l2 + l1 * n_p + l3 * n_p * n_p;

          // Allocate memory for the nodes
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // No boundaries

          // Increment the node number
          node_count++;
        }

        // Last nodes (at the surface x = Lx)
        // Set the local nod number
        local_node_num = l1 * n_p + (n_p - 1) + l3 * n_p * n_p;
        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // Add the node to the boundary 2
        add_boundary_node(2, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }
    }


    // ALL CENTRAL ELEMENT ROWS (WE ARE STILL IN THE LAYER z=0)
    //------------------------

    // Loop over remaining element rows
    for (unsigned i = 1; i < (Ny - 1); i++)
    {
      // Set the first element in the row

      // Allocate memory for element (z = 0) (x =0)
      element_num = Nx * i;
      Element_pt[element_num] = new ELEMENT;

      // The first row of nodes is copied from the element below (at z=0)
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(element_num)->node_pt(l2) =
          finite_element_pt(element_num - Nx)->node_pt((n_p - 1) * n_p + l2);
      }

      // Other rows are new nodes
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // First column of nodes
        // Set the local node number
        local_node_num = l1 * n_p;

        // Allocate memory for the fist node in the x direction (x=0)
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the boundaries 4 and 0
        add_boundary_node(0, Node_pt[node_count]);
        add_boundary_node(4, Node_pt[node_count]);
        // Increment the node number
        node_count++;

        // Now do the other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Set the local node number
          local_node_num = l2 + l1 * n_p;

          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);


          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin;

          // Add the node to the boundary  and 0
          add_boundary_node(0, Node_pt[node_count]);
          // Increment the node number
          node_count++;
        }
      }

      // As always we extend now this element to the third direction
      for (unsigned l3 = 1; l3 < n_p; l3++)
      {
        // The first row of nodes is copied from the element below (at z=0)
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(element_num)->node_pt(l2 + l3 * n_p * n_p) =
            finite_element_pt(element_num - Nx)
              ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
        }

        // Other rows are new nodes (first the nodes for which x=0)
        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // First column of nodes
          // Set the local node number
          local_node_num = l1 * n_p + l3 * n_p * n_p;

          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin;
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // Add the node to the boundary 4
          add_boundary_node(4, Node_pt[node_count]);

          // Increment the node number
          node_count++;

          // Now do the other columns (we extend to the rest of the nodes)
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Set the local node number
            local_node_num = l2 + l1 * n_p + n_p * n_p * l3;

            // Allocate memory for node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer from the element
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

            // No boundaries

            // Increment the node number
            node_count++;
          }
        }
      }

      // Now loop over the rest of the elements in the row, apart from the last
      for (unsigned j = 1; j < (Nx - 1); j++)
      {
        // Allocate memory for new element
        element_num = Nx * i + j;
        Element_pt[element_num] = new ELEMENT;

        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(element_num)->node_pt(l2) =
            finite_element_pt(element_num - Nx)->node_pt((n_p - 1) * n_p + l2);
        }

        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(element_num)->node_pt(l1 * n_p) =
            finite_element_pt(element_num - 1)->node_pt(l1 * n_p + (n_p - 1));

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Set local node number
            local_node_num = l1 * n_p + l2;

            // Allocate memory for the nodes
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) = Zmin;

            // Add the node to the boundary  and 0
            add_boundary_node(0, Node_pt[node_count]);
            // Increment the node number
            node_count++;
          }
        }

        // We extend to the third dimension
        for (unsigned l3 = 1; l3 < n_p; l3++)
        {
          // The first row is copied from the lower element

          for (unsigned l2 = 0; l2 < n_p; l2++)
          {
            finite_element_pt(element_num)->node_pt(l2 + l3 * n_p * n_p) =
              finite_element_pt(element_num - Nx)
                ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
          }

          for (unsigned l1 = 1; l1 < n_p; l1++)
          {
            // First column is same as neighbouring element
            finite_element_pt(element_num)->node_pt(l1 * n_p + l3 * n_p * n_p) =
              finite_element_pt(element_num - 1)
                ->node_pt(l1 * n_p + l3 * n_p * n_p + (n_p - 1));

            // New nodes for other columns
            for (unsigned l2 = 1; l2 < n_p; l2++)
            {
              // Set the local node number
              local_node_num = l1 * n_p + l2 + l3 * n_p * n_p;

              // Allocate memory for the nodes
              Node_pt[node_count] =
                finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
              // Set the pointer
              finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];

              // Get the fractional position of the node
              finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);

              // Set the position of the node
              Node_pt[node_count]->x(0) =
                Xmin + el_length[0] * (j + s_fraction[0]);
              Node_pt[node_count]->x(1) =
                Ymin + el_length[1] * (i + s_fraction[1]);
              Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

              // No boundaries

              // Increment the node number
              node_count++;
            }
          }
        }

      } // End of loop over elements in row


      // Do final element in row

      // Allocate memory for element
      element_num = Nx * i + Nx - 1;
      Element_pt[element_num] = new ELEMENT;

      // We begin with z =0
      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(element_num)->node_pt(l2) =
          finite_element_pt(element_num - Nx)->node_pt((n_p - 1) * n_p + l2);
      }

      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // First column is same as neighbouring element
        finite_element_pt(element_num)->node_pt(l1 * n_p) =
          finite_element_pt(element_num - 1)->node_pt(l1 * n_p + (n_p - 1));

        // Middle nodes
        for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
        {
          // Set local node number
          local_node_num = l1 * n_p + l2;

          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin;

          // Add the node to the boundary  and 0
          add_boundary_node(0, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }

        // Final node

        // Set the local node number
        local_node_num = l1 * n_p + (n_p - 1);

        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);


        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the boundaries - and 1
        add_boundary_node(0, Node_pt[node_count]);
        add_boundary_node(2, Node_pt[node_count]);

        // Increment the node number
        node_count++;

      } // End of loop over rows of nodes in the element

      // We go to the third dimension
      for (unsigned l3 = 1; l3 < n_p; l3++)
      {
        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(element_num)->node_pt(l2 + l3 * n_p * n_p) =
            finite_element_pt(element_num - Nx)
              ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
        }

        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(element_num)->node_pt(l1 * n_p + l3 * n_p * n_p) =
            finite_element_pt(element_num - 1)
              ->node_pt(l1 * n_p + (n_p - 1) + l3 * n_p * n_p);

          // Middle nodes
          for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
          {
            // Set the local node number
            local_node_num = l1 * n_p + l2 + l3 * n_p * n_p;

            // Allocate memory for node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

            // No boundaries

            // Increment the node number
            node_count++;
          }

          // Final node
          // Set the local node number
          local_node_num = l1 * n_p + (n_p - 1) + l3 * n_p * n_p;

          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmax;
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // Add the node to the boundary 2
          add_boundary_node(2, Node_pt[node_count]);
          // Increment the node number
          node_count++;
        } // End of loop over rows of nodes in the element


      } // End of the 3dimension loop


    } // End of loop over rows of elements


    // FINAL ELEMENT ROW (IN THE z=0 LAYER)
    //=================


    // FIRST ELEMENT IN UPPER ROW (upper left corner)
    //----------------------------------------------

    // Allocate memory for element
    element_num = Nx * (Ny - 1);
    Element_pt[element_num] = new ELEMENT;
    // We begin with all the nodes for which z=0
    // The first row of nodes is copied from the element below
    for (unsigned l2 = 0; l2 < n_p; l2++)
    {
      finite_element_pt(element_num)->node_pt(l2) =
        finite_element_pt(element_num - Nx)->node_pt((n_p - 1) * n_p + l2);
    }

    // Second row of  nodes
    // First column of nodes
    for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
    {
      // Set the local node number
      local_node_num = n_p * l1;

      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) =
        Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundaries 4 and 0
      add_boundary_node(0, Node_pt[node_count]);
      add_boundary_node(4, Node_pt[node_count]);
      // Increment the node number
      node_count++;

      // Now do the other columns
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Set the local node number
        local_node_num = n_p * l1 + l2;

        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the boundary 0
        add_boundary_node(0, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }
    }

    // Final row of nodes
    // First column of nodes
    // Top left node
    // Set local node number
    local_node_num = n_p * (n_p - 1);
    // Allocate memory for node
    Node_pt[node_count] =
      finite_element_pt(element_num)
        ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
      Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmin;

    // Add the node to the boundaries 0,3 and 4
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(3, Node_pt[node_count]);
    add_boundary_node(4, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Now do the other columns
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // Set the local node number
      local_node_num = n_p * (n_p - 1) + l2;
      // Allocate memory for the node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundaries
      add_boundary_node(0, Node_pt[node_count]);
      add_boundary_node(3, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // We jump to the third dimension
    for (unsigned l3 = 1; l3 < n_p; l3++)
    {
      // The first row of nodes is copied from the element below
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(element_num)->node_pt(l2 + l3 * n_p * n_p) =
          finite_element_pt(element_num - Nx)
            ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
      }

      // Second row of  nodes
      // First column of nodes
      for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
      {
        // Set the local node number
        local_node_num = n_p * l1 + l3 * n_p * n_p;

        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // Add the node to the boundary 4
        add_boundary_node(4, Node_pt[node_count]);
        // Increment the node number
        node_count++;

        // Now do the other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Set local node number
          local_node_num = n_p * l1 + l2 + l3 * n_p * n_p;

          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // No boundaries

          // Increment the node number
          node_count++;
        }
      }

      // Final row of nodes
      // First column of nodes
      // Top left node
      local_node_num = n_p * (n_p - 1) + l3 * n_p * n_p;
      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      // Add the node to the boundaries
      add_boundary_node(3, Node_pt[node_count]);
      add_boundary_node(4, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Now do the other columns
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        local_node_num = n_p * (n_p - 1) + l2 + l3 * n_p * n_p;
        // Allocate memory for the node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // Add the node to the boundary 3
        add_boundary_node(3, Node_pt[node_count]);


        // Increment the node number
        node_count++;
      }
    }


    // Now loop over the rest of the elements in the row, apart from the last
    // (first plane z = 0)
    for (unsigned j = 1; j < (Nx - 1); j++)
    {
      // Allocate memory for element
      element_num = Nx * (Ny - 1) + j;
      Element_pt[element_num] = new ELEMENT;
      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(element_num)->node_pt(l2) =
          finite_element_pt(element_num - Nx)->node_pt((n_p - 1) * n_p + l2);
      }

      // Second rows
      for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
      {
        // First column is same as neighbouring element
        finite_element_pt(element_num)->node_pt(n_p * l1) =
          finite_element_pt(element_num - 1)->node_pt(n_p * l1 + (n_p - 1));

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          local_node_num = n_p * l1 + l2;
          // Allocate memory for the node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);

          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin;

          // Add the node to the boundary
          add_boundary_node(0, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }
      }

      // Top row
      // First column is same as neighbouring element
      finite_element_pt(element_num)->node_pt(n_p * (n_p - 1)) =
        finite_element_pt(element_num - 1)
          ->node_pt(n_p * (n_p - 1) + (n_p - 1));
      // New nodes for other columns
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        local_node_num = n_p * (n_p - 1) + l2;
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the boundary
        add_boundary_node(3, Node_pt[node_count]);
        add_boundary_node(0, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }


      // Jump in the third dimension

      for (unsigned l3 = 1; l3 < n_p; l3++)
      {
        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(element_num)->node_pt(l2 + l3 * n_p * n_p) =
            finite_element_pt(element_num - Nx)
              ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
        }

        // Second rows
        for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(element_num)->node_pt(n_p * l1 + l3 * n_p * n_p) =
            finite_element_pt(element_num - 1)
              ->node_pt(n_p * l1 + (n_p - 1) + l3 * n_p * n_p);

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            local_node_num = n_p * l1 + l2 + l3 * n_p * n_p;
            // Allocate memory for the node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);

            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
            Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

            // No boundaries

            // Increment the node number
            // add_boundary_node(0,Node_pt[node_count]);
            node_count++;
          }
        }

        // Top row
        // First column is same as neighbouring element
        finite_element_pt(element_num)
          ->node_pt(n_p * (n_p - 1) + l3 * n_p * n_p) =
          finite_element_pt(element_num - 1)
            ->node_pt(n_p * (n_p - 1) + (n_p - 1) + l3 * n_p * n_p);
        // New nodes for other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          local_node_num = n_p * (n_p - 1) + l2 + l3 * n_p * n_p;
          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymax;
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // Add the node to the boundary
          add_boundary_node(3, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }
      }

    } // End of loop over central elements in row


    // FINAL ELEMENT IN ROW (upper right corner) IN LAYER z = 0
    //--------------------------------------------------------

    // Allocate memory for element
    element_num = Nx * (Ny - 1) + Nx - 1;
    Element_pt[element_num] = new ELEMENT;

    // We work first in the plane z =0
    // The first row is copied from the lower element
    for (unsigned l2 = 0; l2 < n_p; l2++)
    {
      finite_element_pt(element_num)->node_pt(l2) =
        finite_element_pt(element_num - Nx)->node_pt((n_p - 1) * n_p + l2);
    }

    // Second rows
    for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
    {
      // First column is same as neighbouring element
      finite_element_pt(element_num)->node_pt(n_p * l1) =
        finite_element_pt(element_num - 1)->node_pt(n_p * l1 + (n_p - 1));

      // Middle nodes
      for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
      {
        local_node_num = n_p * l1 + l2;
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the boundary
        add_boundary_node(0, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }

      // Final node
      local_node_num = n_p * l1 + (n_p - 1);
      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) =
        Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundaries 0 and 2
      add_boundary_node(0, Node_pt[node_count]);
      add_boundary_node(2, Node_pt[node_count]);


      // Increment the node number
      node_count++;

    } // End of loop over middle rows

    // Final row
    // First column is same as neighbouring element
    finite_element_pt(element_num)->node_pt(n_p * (n_p - 1)) =
      finite_element_pt(element_num - 1)->node_pt(n_p * (n_p - 1) + (n_p - 1));

    // Middle nodes
    for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
    {
      local_node_num = n_p * (n_p - 1) + l2;
      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
        Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundaries 0 nd 3
      add_boundary_node(0, Node_pt[node_count]);
      add_boundary_node(3, Node_pt[node_count]);


      // Increment the node number
      node_count++;
    }

    // Final node
    // Determine number of values
    local_node_num = n_p * (n_p - 1) + (n_p - 1);
    // Allocate memory for node
    Node_pt[node_count] =
      finite_element_pt(element_num)
        ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer
    finite_element_pt(element_num)->node_pt(local_node_num) =
      Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmin;

    // Add the node to the boundaries 0,2 and 3
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(2, Node_pt[node_count]);
    add_boundary_node(3, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // We jump to the third dimension

    for (unsigned l3 = 1; l3 < n_p; l3++)
    {
      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(element_num)->node_pt(l2 + l3 * n_p * n_p) =
          finite_element_pt(element_num - Nx)
            ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
      }

      // Second rows
      for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
      {
        // First column is same as neighbouring element
        finite_element_pt(element_num)->node_pt(n_p * l1 + l3 * n_p * n_p) =
          finite_element_pt(element_num - 1)
            ->node_pt(n_p * l1 + (n_p - 1) + l3 * n_p * n_p);

        // Middle nodes
        for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
        {
          // Determine number of values
          local_node_num = n_p * l1 + l2 + l3 * n_p * n_p;
          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // No boundaries

          // Increment the node number
          node_count++;
        }

        // Final node
        // Determine number of values
        local_node_num = n_p * l1 + (n_p - 1) + l3 * n_p * n_p;
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // Add the node to the boundary 2
        add_boundary_node(2, Node_pt[node_count]);
        // Increment the node number
        node_count++;

      } // End of loop over middle rows

      // Final row
      // First column is same as neighbouring element
      finite_element_pt(element_num)
        ->node_pt(n_p * (n_p - 1) + l3 * n_p * n_p) =
        finite_element_pt(element_num - 1)
          ->node_pt(n_p * (n_p - 1) + (n_p - 1) + l3 * n_p * n_p);

      // Middle nodes
      for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
      {
        // Determine number of values
        local_node_num = n_p * (n_p - 1) + l2 + l3 * n_p * n_p;
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // Add the node to the boundary 3
        add_boundary_node(3, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // Final node
      // Determine number of values
      local_node_num = n_p * (n_p - 1) + (n_p - 1) + l3 * n_p * n_p;
      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      // Add the node to the boundaries 2 and 3
      add_boundary_node(2, Node_pt[node_count]);
      add_boundary_node(3, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // END OF THE FIRST LAYER


    //----------------------------------------------------------------------------------------------------------------------------
    // ***************************************NOW WE MAKE ALL THE INTREMEDIATE
    // LAYERS**********************************************
    //----------------------------------------------------------------------------------------------------------------------------


    for (unsigned k = 1; k < (Nz - 1); k++) // good loop for the diferent layers
    // for(unsigned k=1;k<Nz;k++)  // bad loop for testing the hole mesh but the
    // last layer
    {
      // FIRST ELEMENT OF THE LAYER
      //----------------------------------

      element_num = k * Nx * Ny;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < n_p; l1++)
      {
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(element_num)->node_pt(l2 + n_p * l1) =
            finite_element_pt(element_num - Nx * Ny)
              ->node_pt(l2 + n_p * l1 + n_p * n_p * (n_p - 1));
        }
      }


      // Loop over the other node columns in the z direction

      for (unsigned l3 = 1; l3 < n_p; l3++)
      {
        // Set the corner node
        // Determine number of values at this node
        local_node_num = n_p * n_p * l3;

        // Allocate memory for the node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);

        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        // Add the node to boundaries 1 and 4
        add_boundary_node(1, Node_pt[node_count]);
        add_boundary_node(4, Node_pt[node_count]);

        // Increment the node number
        node_count++;

        // Loop over the other nodes in the first row in the x direction
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Determine the number of values at this node
          local_node_num = l2 + n_p * n_p * l3;

          // Allocate memory for the nodes
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) = Ymin;
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // Add the node to the boundary 1
          add_boundary_node(1, Node_pt[node_count]);
          // Increment the node number
          node_count++;
        }

        // Loop over the other node columns
        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // Determine the number of values
          local_node_num = l1 * n_p + n_p * n_p * l3;

          // Allocate memory for the nodes
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the first node of the row (with boundary 4)
          Node_pt[node_count]->x(0) = Xmin;
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // Add the node to the boundary 4
          add_boundary_node(4, Node_pt[node_count]);
          // Increment the node number
          node_count++;

          // Loop over the other nodes in the row
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Set the number of values
            local_node_num = l1 * n_p + l2 + n_p * n_p * l3;

            // Allocate the memory for the node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer from the element
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
            Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (k + s_fraction[2]);


            // No boundary

            // Increment the node number
            node_count++;
          }
        }
      }


      //----------------------------------------------------------------------

      // CENTRE OF FIRST ROW OF ELEMENTS
      //--------------------------------

      // Now loop over the first row of elements, apart from final element
      for (unsigned j = 1; j < (Nx - 1); j++)
      {
        // Allocate memory for new element
        element_num = j + k * Nx * Ny;
        Element_pt[element_num] = new ELEMENT;

        // The lowest nodes are copied from the lower element
        for (unsigned l1 = 0; l1 < n_p; l1++)
        {
          for (unsigned l2 = 0; l2 < n_p; l2++)
          {
            finite_element_pt(j + k * Nx * Ny)->node_pt(l2 + n_p * l1) =
              finite_element_pt(j + (k - 1) * Nx * Ny)
                ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
          }
        }

        // Loop over the other node columns in the z direction
        for (unsigned l3 = 1; l3 < n_p; l3++)
        {
          // First column of nodes is same as neighbouring element
          finite_element_pt(j + k * Nx * Ny)->node_pt(l3 * n_p * n_p) =
            finite_element_pt(j - 1 + k * Nx * Ny)
              ->node_pt(l3 * n_p * n_p + (n_p - 1));

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Determine number of values
            local_node_num = l2 + l3 * n_p * n_p;

            // Allocate memory for the nodes
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
            // Set the pointer from the element
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) = Ymin;
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (k + s_fraction[2]);

            // Add the node to the boundary 1
            add_boundary_node(1, Node_pt[node_count]);
            // Increment the node number
            node_count++;
          }

          // Do the rest of the nodes
          for (unsigned l1 = 1; l1 < n_p; l1++)
          {
            // First column of nodes is same as neighbouring element
            finite_element_pt(j + k * Nx * Ny)
              ->node_pt(l1 * n_p + l3 * n_p * n_p) =
              finite_element_pt(j - 1 + k * Nx * Ny)
                ->node_pt(l1 * n_p + (n_p - 1) + l3 * n_p * n_p);

            // New nodes for other columns
            for (unsigned l2 = 1; l2 < n_p; l2++)
            {
              // Determine number of values
              local_node_num = l1 * n_p + l2 + l3 * n_p * n_p;

              // Allocate memory for the nodes
              Node_pt[node_count] =
                finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
              // Set the pointer from the element
              finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];

              // Get the fractional position of the node
              finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);

              // Set the position of the node
              Node_pt[node_count]->x(0) =
                Xmin + el_length[0] * (j + s_fraction[0]);
              Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
              Node_pt[node_count]->x(2) =
                Zmin + el_length[2] * (k + s_fraction[2]);

              // No boundaries

              // Increment the node number
              node_count++;
            }
          }
        }
      }

      // MY FINAL ELEMENT IN FIRST ROW (right corner)
      //-----------------------------------------------


      // Allocate memory for new element
      element_num = Nx - 1 + k * Nx * Ny;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < n_p; l1++)
      {
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx - 1 + k * Nx * Ny)->node_pt(l2 + n_p * l1) =
            finite_element_pt(Nx - 1 + (k - 1) * Nx * Ny)
              ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
        }
      }


      // Loop over the other node columns in the z direction
      for (unsigned l3 = 1; l3 < n_p; l3++)
      {
        // First column of nodes is same as neighbouring element
        finite_element_pt(Nx - 1 + k * Nx * Ny)->node_pt(l3 * n_p * n_p) =
          finite_element_pt(Nx - 2 + k * Nx * Ny)
            ->node_pt(l3 * n_p * n_p + (n_p - 1));

        // New nodes for other rows (y=0)
        for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
        {
          // Determine number of values
          local_node_num = l2 + l3 * n_p * n_p;

          // Allocate memory for the nodes
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin;
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // Add the node to the boundary 1
          add_boundary_node(1, Node_pt[node_count]);
          // Increment the node number
          node_count++;
        }

        // Last node of the row (y=0)

        // Determine number of values
        local_node_num = (n_p - 1) + l3 * n_p * n_p;

        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        // Add the node to the boundary 1 and 2
        add_boundary_node(1, Node_pt[node_count]);
        add_boundary_node(2, Node_pt[node_count]);
        // Increment the node number
        node_count++;

        // Do the rest of the nodes
        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // First column of nodes is same as neighbouring element
          finite_element_pt(Nx - 1 + k * Nx * Ny)
            ->node_pt(l1 * n_p + l3 * n_p * n_p) =
            finite_element_pt(Nx - 2 + k * Nx * Ny)
              ->node_pt(l1 * n_p + (n_p - 1) + l3 * n_p * n_p);

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
          {
            // Determine number of values
            local_node_num = l1 * n_p + l2 + l3 * n_p * n_p;

            // Allocate memory for the nodes
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer from the element
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
            Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (k + s_fraction[2]);

            // No boundaries

            // Increment the node number
            node_count++;
          }

          // Last nodes (at the surface x = Lx)
          // Determine number of values
          local_node_num = l1 * n_p + (n_p - 1) + l3 * n_p * n_p;
          // Allocate memory for the nodes
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);


          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmax;
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // Add the node to the boundary 2
          add_boundary_node(2, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }
      }


      // ALL CENTRAL ELEMENT ROWS
      //------------------------

      // Loop over remaining element rows
      for (unsigned i = 1; i < (Ny - 1); i++)
      {
        // Set the first element in the row

        // Allocate memory for element (x =0)
        element_num = Nx * i + Nx * Ny * k;
        Element_pt[element_num] = new ELEMENT;

        // The lowest nodes are copied from the lower element
        for (unsigned l1 = 0; l1 < n_p; l1++)
        {
          for (unsigned l2 = 0; l2 < n_p; l2++)
          {
            finite_element_pt(Nx * i + k * Nx * Ny)->node_pt(l2 + n_p * l1) =
              finite_element_pt(Nx * i + (k - 1) * Nx * Ny)
                ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
          }
        }

        // We extend now this element to the third direction

        for (unsigned l3 = 1; l3 < n_p; l3++)
        {
          // The first row of nodes is copied from the element below (at z=0)
          for (unsigned l2 = 0; l2 < n_p; l2++)
          {
            finite_element_pt(Nx * i + k * Nx * Ny)
              ->node_pt(l2 + l3 * n_p * n_p) =
              finite_element_pt(Nx * (i - 1) + k * Nx * Ny)
                ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
          }

          // Other rows are new nodes (first the nodes for which x=0)
          for (unsigned l1 = 1; l1 < n_p; l1++)
          {
            // First column of nodes

            // Determine number of values
            local_node_num = l1 * n_p + l3 * n_p * n_p;

            // Allocate memory for node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
            // Set the pointer from the element
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) = Xmin;
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (k + s_fraction[2]);

            // Add the node to the boundary 4
            add_boundary_node(4, Node_pt[node_count]);

            // Increment the node number
            node_count++;

            // Now do the other columns (we extend to the rest of the nodes)
            for (unsigned l2 = 1; l2 < n_p; l2++)
            {
              // Determine number of values
              local_node_num = l1 * n_p + l2 + n_p * n_p * l3;

              // Allocate memory for node
              Node_pt[node_count] =
                finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
              // Set the pointer from the element
              finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];

              // Get the fractional position of the node
              finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);

              // Set the position of the node
              Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
              Node_pt[node_count]->x(1) =
                Ymin + el_length[1] * (i + s_fraction[1]);
              Node_pt[node_count]->x(2) =
                Zmin + el_length[2] * (k + s_fraction[2]);

              // No boundaries

              // Increment the node number
              node_count++;
            }
          }
        }

        // Now loop over the rest of the elements in the row, apart from the
        // last
        for (unsigned j = 1; j < (Nx - 1); j++)
        {
          // Allocate memory for new element
          element_num = Nx * i + j + k * Nx * Ny;
          Element_pt[element_num] = new ELEMENT;

          // The lowest nodes are copied from the lower element
          for (unsigned l1 = 0; l1 < n_p; l1++)
          {
            for (unsigned l2 = 0; l2 < n_p; l2++)
            {
              finite_element_pt(Nx * i + j + k * Nx * Ny)
                ->node_pt(l2 + n_p * l1) =
                finite_element_pt(Nx * i + j + (k - 1) * Nx * Ny)
                  ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
            }
          }

          // We extend to the third dimension

          for (unsigned l3 = 1; l3 < n_p; l3++)
          {
            // The first row is copied from the lower element

            for (unsigned l2 = 0; l2 < n_p; l2++)
            {
              finite_element_pt(Nx * i + j + k * Nx * Ny)
                ->node_pt(l2 + l3 * n_p * n_p) =
                finite_element_pt(Nx * (i - 1) + j + k * Nx * Ny)
                  ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
            }

            for (unsigned l1 = 1; l1 < n_p; l1++)
            {
              // First column is same as neighbouring element
              finite_element_pt(Nx * i + j + k * Nx * Ny)
                ->node_pt(l1 * n_p + l3 * n_p * n_p) =
                finite_element_pt(Nx * i + (j - 1) + k * Nx * Ny)
                  ->node_pt(l1 * n_p + l3 * n_p * n_p + (n_p - 1));

              // New nodes for other columns
              for (unsigned l2 = 1; l2 < n_p; l2++)
              {
                // Determine number of values
                local_node_num = l1 * n_p + l2 + l3 * n_p * n_p;

                // Allocate memory for the nodes
                Node_pt[node_count] =
                  finite_element_pt(element_num)
                    ->construct_node(local_node_num, time_stepper_pt);
                // Set the pointer
                finite_element_pt(element_num)->node_pt(local_node_num) =
                  Node_pt[node_count];
                // Get the fractional position of the node
                finite_element_pt(element_num)
                  ->local_fraction_of_node(local_node_num, s_fraction);


                // Set the position of the node
                Node_pt[node_count]->x(0) =
                  Xmin + el_length[0] * (j + s_fraction[0]);
                Node_pt[node_count]->x(1) =
                  Ymin + el_length[1] * (i + s_fraction[1]);
                Node_pt[node_count]->x(2) =
                  Zmin + el_length[2] * (k + s_fraction[2]);


                // No boundaries
                // Increment the node number
                node_count++;
              }
            }
          }

        } // End of loop over elements in row


        // Do final element in row

        // Allocate memory for element
        element_num = Nx * i + Nx - 1 + k * Nx * Ny;
        Element_pt[element_num] = new ELEMENT;

        // The lowest nodes are copied from the lower element
        for (unsigned l1 = 0; l1 < n_p; l1++)
        {
          for (unsigned l2 = 0; l2 < n_p; l2++)
          {
            finite_element_pt(Nx * i + Nx - 1 + k * Nx * Ny)
              ->node_pt(l2 + n_p * l1) =
              finite_element_pt(Nx * i + Nx - 1 + (k - 1) * Nx * Ny)
                ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
          }
        }


        // We go to the third dimension
        for (unsigned l3 = 1; l3 < n_p; l3++)
        {
          // The first row is copied from the lower element
          for (unsigned l2 = 0; l2 < n_p; l2++)
          {
            finite_element_pt(Nx * i + Nx - 1 + k * Nx * Ny)
              ->node_pt(l2 + l3 * n_p * n_p) =
              finite_element_pt(Nx * (i - 1) + Nx - 1 + k * Nx * Ny)
                ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
          }

          for (unsigned l1 = 1; l1 < n_p; l1++)
          {
            // First column is same as neighbouring element
            finite_element_pt(Nx * i + Nx - 1 + k * Nx * Ny)
              ->node_pt(l1 * n_p + l3 * n_p * n_p) =
              finite_element_pt(Nx * i + Nx - 2 + k * Nx * Ny)
                ->node_pt(l1 * n_p + (n_p - 1) + l3 * n_p * n_p);

            // Middle nodes
            for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
            {
              // Determine number of values
              local_node_num = l1 * n_p + l2 + l3 * n_p * n_p;

              // Allocate memory for node
              Node_pt[node_count] =
                finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
              // Set the pointer
              finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];

              // Get the fractional position of the node
              finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);

              // Set the position of the node
              Node_pt[node_count]->x(0) =
                Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
              Node_pt[node_count]->x(1) =
                Ymin + el_length[1] * (i + s_fraction[1]);
              Node_pt[node_count]->x(2) =
                Zmin + el_length[2] * (k + s_fraction[2]);

              // No boundaries

              // Increment the node number
              node_count++;
            }

            // Final node

            // Determine number of values
            local_node_num = l1 * n_p + (n_p - 1) + l3 * n_p * n_p;

            // Allocate memory for node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) = Xmax;
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (k + s_fraction[2]);

            // Add the node to the boundary 2
            add_boundary_node(2, Node_pt[node_count]);

            // Increment the node number
            node_count++;
          } // End of loop over rows of nodes in the element


        } // End of the 3dimension loop


      } // End of loop over rows of elements


      // FINAL ELEMENT ROW (IN INTERMIDIATE  LAYERS)
      //=================


      // FIRST ELEMENT IN UPPER ROW (upper left corner)
      //----------------------------------------------

      // Allocate memory for element
      element_num = Nx * (Ny - 1) + k * Nx * Ny;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < n_p; l1++)
      {
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * (Ny - 1) + k * Nx * Ny)
            ->node_pt(l2 + n_p * l1) =
            finite_element_pt(Nx * (Ny - 1) + (k - 1) * Nx * Ny)
              ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
        }
      }

      // We jump to the third dimension
      for (unsigned l3 = 1; l3 < n_p; l3++)
      {
        // The first row of nodes is copied from the element below
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * (Ny - 1) + k * Nx * Ny)
            ->node_pt(l2 + l3 * n_p * n_p) =
            finite_element_pt(Nx * (Ny - 2) + k * Nx * Ny)
              ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
        }

        // Second row of  nodes
        // First column of nodes
        for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
        {
          // Determine number of values
          local_node_num = n_p * l1 + l3 * n_p * n_p;

          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin;
          Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // Add the node to the boundary 4
          add_boundary_node(4, Node_pt[node_count]);

          // Increment the node number
          node_count++;

          // Now do the other columns
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Determine number of values
            local_node_num = n_p * l1 + l2 + l3 * n_p * n_p;

            // Allocate memory for node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer from the element
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (k + s_fraction[2]);

            // No boundaries

            // Increment the node number
            node_count++;
          }
        }

        // Final row of nodes
        // First column of nodes
        // Top left node
        // Determine number of values
        local_node_num = n_p * (n_p - 1) + l3 * n_p * n_p;
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        // Add the node to the boundaries
        add_boundary_node(3, Node_pt[node_count]);
        add_boundary_node(4, Node_pt[node_count]);

        // Increment the node number
        node_count++;

        // Now do the other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Determine number of values
          local_node_num = n_p * (n_p - 1) + l2 + l3 * n_p * n_p;
          // Allocate memory for the node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) = Ymax;
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // Add the node to the boundary 3
          add_boundary_node(3, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }
      }


      // Now loop over the rest of the elements in the row, apart from the last
      for (unsigned j = 1; j < (Nx - 1); j++)
      {
        // Allocate memory for element
        element_num = Nx * (Ny - 1) + j + k * Nx * Ny;
        Element_pt[element_num] = new ELEMENT;

        // The lowest nodes are copied from the lower element
        for (unsigned l1 = 0; l1 < n_p; l1++)
        {
          for (unsigned l2 = 0; l2 < n_p; l2++)
          {
            finite_element_pt(Nx * (Ny - 1) + j + k * Nx * Ny)
              ->node_pt(l2 + n_p * l1) =
              finite_element_pt(Nx * (Ny - 1) + j + (k - 1) * Nx * Ny)
                ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
          }
        }

        // Jump in the third dimension

        for (unsigned l3 = 1; l3 < n_p; l3++)
        {
          // The first row is copied from the lower element
          for (unsigned l2 = 0; l2 < n_p; l2++)
          {
            finite_element_pt(Nx * (Ny - 1) + j + k * Nx * Ny)
              ->node_pt(l2 + l3 * n_p * n_p) =
              finite_element_pt(Nx * (Ny - 2) + j + k * Nx * Ny)
                ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
          }

          // Second rows
          for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
          {
            // First column is same as neighbouring element
            finite_element_pt(Nx * (Ny - 1) + j + k * Nx * Ny)
              ->node_pt(n_p * l1 + l3 * n_p * n_p) =
              finite_element_pt(Nx * (Ny - 1) + (j - 1) + k * Nx * Ny)
                ->node_pt(n_p * l1 + (n_p - 1) + l3 * n_p * n_p);

            // New nodes for other columns
            for (unsigned l2 = 1; l2 < n_p; l2++)
            {
              // Determine number of values
              local_node_num = n_p * l1 + l2 + l3 * n_p * n_p;
              // Allocate memory for the node
              Node_pt[node_count] =
                finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);

              // Set the pointer
              finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];

              // Get the fractional position of the node
              finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);

              // Set the position of the node
              Node_pt[node_count]->x(0) =
                Xmin + el_length[0] * (j + s_fraction[0]);
              Node_pt[node_count]->x(1) =
                Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
              Node_pt[node_count]->x(2) =
                Zmin + el_length[2] * (k + s_fraction[2]);

              // No boundaries

              // Increment the node number
              // add_boundary_node(0,Node_pt[node_count]);
              node_count++;
            }
          }

          // Top row
          // First column is same as neighbouring element
          finite_element_pt(Nx * (Ny - 1) + j + k * Nx * Ny)
            ->node_pt(n_p * (n_p - 1) + l3 * n_p * n_p) =
            finite_element_pt(Nx * (Ny - 1) + (j - 1) + k * Nx * Ny)
              ->node_pt(n_p * (n_p - 1) + (n_p - 1) + l3 * n_p * n_p);
          // New nodes for other columns
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Determine number of values
            local_node_num = n_p * (n_p - 1) + l2 + l3 * n_p * n_p;
            // Allocate memory for node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) = Ymax;
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (k + s_fraction[2]);

            // Add the node to the boundary
            add_boundary_node(3, Node_pt[node_count]);

            // Increment the node number
            node_count++;
          }
        }

      } // End of loop over central elements in row


      // FINAL ELEMENT IN ROW (upper right corner)
      //-----------------------------------------

      // Allocate memory for element
      element_num = Nx * (Ny - 1) + Nx - 1 + k * Nx * Ny;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < n_p; l1++)
      {
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * (Ny - 1) + Nx - 1 + k * Nx * Ny)
            ->node_pt(l2 + n_p * l1) =
            finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (k - 1) * Nx * Ny)
              ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
        }
      }

      // We jump to the third dimension


      for (unsigned l3 = 1; l3 < n_p; l3++)
      {
        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * (Ny - 1) + Nx - 1 + k * Nx * Ny)
            ->node_pt(l2 + l3 * n_p * n_p) =
            finite_element_pt(Nx * (Ny - 2) + Nx - 1 + k * Nx * Ny)
              ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
        }

        // Second rows
        for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(Nx * (Ny - 1) + Nx - 1 + k * Nx * Ny)
            ->node_pt(n_p * l1 + l3 * n_p * n_p) =
            finite_element_pt(Nx * (Ny - 1) + Nx - 2 + k * Nx * Ny)
              ->node_pt(n_p * l1 + (n_p - 1) + l3 * n_p * n_p);

          // Middle nodes
          for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
          {
            // Determine number of values
            local_node_num = n_p * l1 + l2 + l3 * n_p * n_p;
            // Allocate memory for node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];
            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);


            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (k + s_fraction[2]);

            // No boundaries

            // Increment the node number
            node_count++;
          }

          // Final node
          // Determine number of values
          local_node_num = n_p * l1 + (n_p - 1) + l3 * n_p * n_p;
          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmax;
          Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // Add the node to the boundary 2
          add_boundary_node(2, Node_pt[node_count]);

          // Increment the node number
          node_count++;

        } // End of loop over middle rows

        // Final row
        // First column is same as neighbouring element
        finite_element_pt(Nx * (Ny - 1) + Nx - 1 + k * Nx * Ny)
          ->node_pt(n_p * (n_p - 1) + l3 * n_p * n_p) =
          finite_element_pt(Nx * (Ny - 1) + Nx - 2 + k * Nx * Ny)
            ->node_pt(n_p * (n_p - 1) + (n_p - 1) + l3 * n_p * n_p);

        // Middle nodes
        for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
        {
          // Determine number of values
          local_node_num = n_p * (n_p - 1) + l2 + l3 * n_p * n_p;
          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymax;
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // Add the node to the boundary 3
          add_boundary_node(3, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }

        // Final node
        // Determine number of values
        local_node_num = n_p * (n_p - 1) + (n_p - 1) + l3 * n_p * n_p;
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        // Add the node to the boundaries 2 and 3
        add_boundary_node(2, Node_pt[node_count]);
        add_boundary_node(3, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }

    } // End loop of the elements layer


    // END OF THE INTERMEDIATE LAYERS

    // ----------------------------------------------------------------------------------
    //  ****************BEGINNING OF THE LAST
    //  LAYER**************************************
    // ----------------------------------------------------------------------------------


    // FIRST ELEMENT OF THE UPPER LAYER
    //----------------------------------

    element_num = (Nz - 1) * Nx * Ny;
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < n_p; l1++)
    {
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt((Nz - 1) * Nx * Ny)->node_pt(l2 + n_p * l1) =
          finite_element_pt((Nz - 2) * Nx * Ny)
            ->node_pt(l2 + n_p * l1 + n_p * n_p * (n_p - 1));
      }
    }


    // Loop over the other node columns in the z direction but the last

    for (unsigned l3 = 1; l3 < (n_p - 1); l3++)
    {
      // Set the corner nodes
      // Determine number of values at this node
      local_node_num = n_p * n_p * l3;

      // Allocate memory for the node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);

      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];
      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);


      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) =
        Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      // Add the node to boundaries 1 and 4
      add_boundary_node(1, Node_pt[node_count]);
      add_boundary_node(4, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Loop over the other nodes in the first row in the x direction
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Determine the number of values at this node
        local_node_num = l2 + n_p * n_p * l3;

        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];
        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);


        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // Add the node to the boundary 1
        add_boundary_node(1, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // Loop over the other node columns
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // Determine the number of values
        local_node_num = l1 * n_p + n_p * n_p * l3;

        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the first node of the row (with boundary 4)
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // Add the node to the boundary 4
        add_boundary_node(4, Node_pt[node_count]);
        // Increment the node number
        node_count++;

        // Loop over the other nodes in the row
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Set the number of values
          local_node_num = l1 * n_p + l2 + n_p * n_p * l3;

          // Allocate the memory for the node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // No boundary

          // Increment the node number
          node_count++;
        }
      }
    }

    // Make the top nodes

    // Set the corner nodes
    // Determine number of values at this node
    local_node_num = n_p * n_p * (n_p - 1);

    // Allocate memory for the node
    Node_pt[node_count] =
      finite_element_pt(element_num)
        ->construct_boundary_node(local_node_num, time_stepper_pt);

    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
      Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmax;

    // Add the node to boundaries 1, 4 and 5
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(4, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Loop over the other nodes in the first row in the x direction
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // Determine the number of values at this node
      local_node_num = l2 + n_p * n_p * (n_p - 1);

      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundaries 1 and 5
      add_boundary_node(1, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // Loop over the other node columns
    for (unsigned l1 = 1; l1 < n_p; l1++)
    {
      // Determine the number of values
      local_node_num = l1 * n_p + n_p * n_p * (n_p - 1);

      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the first node of the row (with boundary 4)
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundary 4
      add_boundary_node(4, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);
      // Increment the node number
      node_count++;

      // Loop over the other nodes in the row
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Set the number of values
        local_node_num = l1 * n_p + l2 + n_p * n_p * (n_p - 1);

        // Allocate the memory for the node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmax;

        // Add the node to the boundary 5
        add_boundary_node(5, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    }

    //----------------------------------------------------------------------

    // CENTRE OF FIRST ROW OF ELEMENTS OF THE UPPER LAYER
    //--------------------------------------------------

    // Now loop over the first row of elements, apart from final element
    for (unsigned j = 1; j < (Nx - 1); j++)
    {
      // Allocate memory for new element
      element_num = j + (Nz - 1) * Nx * Ny;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < n_p; l1++)
      {
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(j + (Nz - 1) * Nx * Ny)->node_pt(l2 + n_p * l1) =
            finite_element_pt(j + (Nz - 2) * Nx * Ny)
              ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
        }
      }

      // Loop over the other node columns in the z direction but the last
      for (unsigned l3 = 1; l3 < (n_p - 1); l3++)
      {
        // First column of nodes is same as neighbouring element
        finite_element_pt(j + (Nz - 1) * Nx * Ny)->node_pt(l3 * n_p * n_p) =
          finite_element_pt(j - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt(l3 * n_p * n_p + (n_p - 1));

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Determine number of values
          local_node_num = l2 + l3 * n_p * n_p;

          // Allocate memory for the nodes
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin;
          Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // Add the node to the boundary 1
          add_boundary_node(1, Node_pt[node_count]);
          // Increment the node number
          node_count++;
        }

        // Do the rest of the nodes
        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // First column of nodes is same as neighbouring element
          finite_element_pt(j + (Nz - 1) * Nx * Ny)
            ->node_pt(l1 * n_p + l3 * n_p * n_p) =
            finite_element_pt(j - 1 + (Nz - 1) * Nx * Ny)
              ->node_pt(l1 * n_p + (n_p - 1) + l3 * n_p * n_p);

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Determine number of values
            local_node_num = l1 * n_p + l2 + l3 * n_p * n_p;

            // Allocate memory for the nodes
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer from the element
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

            // No boundaries

            // Increment the node number
            node_count++;
          }
        }
      }


      // Top nodes

      // First node is same as neighbouring element
      finite_element_pt(j + (Nz - 1) * Nx * Ny)
        ->node_pt((n_p - 1) * n_p * n_p) =
        finite_element_pt(j - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt((n_p - 1) * n_p * n_p + (n_p - 1));

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Determine number of values
        local_node_num = l2 + (n_p - 1) * n_p * n_p;

        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) = Zmax;

        // Add the node to the boundaries 1 and 5
        add_boundary_node(1, Node_pt[node_count]);
        add_boundary_node(5, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // Do the rest of the nodes
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // First column of nodes is same as neighbouring element
        finite_element_pt(j + (Nz - 1) * Nx * Ny)
          ->node_pt(l1 * n_p + (n_p - 1) * n_p * n_p) =
          finite_element_pt(j - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt(l1 * n_p + (n_p - 1) + (n_p - 1) * n_p * n_p);

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Determine number of values
          local_node_num = l1 * n_p + l2 + (n_p - 1) * n_p * n_p;

          // Allocate memory for the nodes
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) = Zmax;

          // Add to the boundary 5
          add_boundary_node(5, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }
      }
    } // End loop of first row of elements


    // MY FINAL ELEMENT IN THE FIRST ROW OF THE UPPER LAYER (right corner)
    //--------------------------------------------------------------------


    // Allocate memory for new element
    element_num = Nx - 1 + (Nz - 1) * Nx * Ny;
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < n_p; l1++)
    {
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(Nx - 1 + (Nz - 1) * Nx * Ny)->node_pt(l2 + n_p * l1) =
          finite_element_pt(Nx - 1 + (Nz - 2) * Nx * Ny)
            ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
      }
    }


    // Loop over the other node columns in the z direction but the last
    for (unsigned l3 = 1; l3 < (n_p - 1); l3++)
    {
      // First column of nodes is same as neighbouring element
      finite_element_pt(Nx - 1 + (Nz - 1) * Nx * Ny)->node_pt(l3 * n_p * n_p) =
        finite_element_pt(Nx - 2 + (Nz - 1) * Nx * Ny)
          ->node_pt(l3 * n_p * n_p + (n_p - 1));

      // New nodes for other rows (y=0)
      for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
      {
        // Determine number of values
        local_node_num = l2 + l3 * n_p * n_p;

        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // Add the node to the boundary 1
        add_boundary_node(1, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // Last node of the row (y=0)

      // Determine number of values
      local_node_num = (n_p - 1) + l3 * n_p * n_p;

      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) =
        Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      // Add the node to the boundary 1 and 2
      add_boundary_node(1, Node_pt[node_count]);
      add_boundary_node(2, Node_pt[node_count]);
      // Increment the node number
      node_count++;

      // Do the rest of the nodes
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // First column of nodes is same as neighbouring element
        finite_element_pt(Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(l1 * n_p + l3 * n_p * n_p) =
          finite_element_pt(Nx - 2 + (Nz - 1) * Nx * Ny)
            ->node_pt(l1 * n_p + (n_p - 1) + l3 * n_p * n_p);

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
        {
          // Determine number of values
          local_node_num = l1 * n_p + l2 + l3 * n_p * n_p;

          // Allocate memory for the nodes
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // No boundaries

          // Increment the node number
          node_count++;
        }

        // Last nodes (at the surface x = Lx)
        // Determine number of values
        local_node_num = l1 * n_p + (n_p - 1) + l3 * n_p * n_p;
        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // Add the node to the boundary 2
        add_boundary_node(2, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    }

    // We make the top nodes
    // First node is same as neighbouring element
    finite_element_pt(Nx - 1 + (Nz - 1) * Nx * Ny)
      ->node_pt((n_p - 1) * n_p * n_p) =
      finite_element_pt(Nx - 2 + (Nz - 1) * Nx * Ny)
        ->node_pt((n_p - 1) * n_p * n_p + (n_p - 1));

    // New nodes for other rows (y=0)
    for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
    {
      // Determine number of values
      local_node_num = l2 + (n_p - 1) * n_p * n_p;

      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
        Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundaries 1 and 5
      add_boundary_node(1, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // Last node of the row (y=0)

    // Determine number of values
    local_node_num = (n_p - 1) + (n_p - 1) * n_p * n_p;

    // Allocate memory for the nodes
    Node_pt[node_count] =
      finite_element_pt(element_num)
        ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
      Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmax;

    // Add the node to the boundary 1 and 2
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(2, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);
    // Increment the node number
    node_count++;

    // Do the rest of the nodes
    for (unsigned l1 = 1; l1 < n_p; l1++)
    {
      // First column of nodes is same as neighbouring element
      finite_element_pt(Nx - 1 + (Nz - 1) * Nx * Ny)
        ->node_pt(l1 * n_p + (n_p - 1) * n_p * n_p) =
        finite_element_pt(Nx - 2 + (Nz - 1) * Nx * Ny)
          ->node_pt(l1 * n_p + (n_p - 1) + (n_p - 1) * n_p * n_p);

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
      {
        // Determine number of values
        local_node_num = l1 * n_p + l2 + (n_p - 1) * n_p * n_p;

        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmax;

        // Add node to the boundary 5
        add_boundary_node(5, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }

      // Last nodes (at the surface x = Lx)
      // Determine number of values
      local_node_num = l1 * n_p + (n_p - 1) + (n_p - 1) * n_p * n_p;
      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundaries 2 and 5
      add_boundary_node(2, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);


      // Increment the node number
      node_count++;
    }


    // ALL CENTRAL ELEMENT ROWS OF THE TOP  LAYER
    //------------------------------------------

    // Loop over remaining element rows
    for (unsigned i = 1; i < (Ny - 1); i++)
    {
      // Set the first element in the row

      // Allocate memory for element (x =0)
      element_num = Nx * i + Nx * Ny * (Nz - 1);
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < n_p; l1++)
      {
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * i + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + n_p * l1) =
            finite_element_pt(Nx * i + (Nz - 2) * Nx * Ny)
              ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
        }
      }

      // We extend now this element to the third dimension

      for (unsigned l3 = 1; l3 < (n_p - 1); l3++)
      {
        // The first row of nodes is copied from the element below (at z=0)
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * i + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + l3 * n_p * n_p) =
            finite_element_pt(Nx * (i - 1) + (Nz - 1) * Nx * Ny)
              ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
        }

        // Other rows are new nodes (first the nodes for which x=0)
        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // First column of nodes

          // Determine number of values
          local_node_num = l1 * n_p + l3 * n_p * n_p;

          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin;
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // Add the node to the boundary 4
          add_boundary_node(4, Node_pt[node_count]);

          // Increment the node number
          node_count++;

          // Now do the other columns (we extend to the rest of the nodes)
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Determine number of values
            local_node_num = l1 * n_p + l2 + n_p * n_p * l3;

            // Allocate memory for node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer from the element
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

            // No boundaries

            // Increment the node number
            node_count++;
          }
        }
      }

      // Now the top nodes

      // The first row of nodes is copied from the element below (at z=0)
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(Nx * i + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + (n_p - 1) * n_p * n_p) =
          finite_element_pt(Nx * (i - 1) + (Nz - 1) * Nx * Ny)
            ->node_pt((n_p - 1) * n_p + l2 + (n_p - 1) * n_p * n_p);
      }

      // Other rows are new nodes (first the nodes for which x=0)
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // First column of nodes

        // Determine number of values
        local_node_num = l1 * n_p + (n_p - 1) * n_p * n_p;

        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmax;

        // Add the node to the boundaries 4 and 5
        add_boundary_node(4, Node_pt[node_count]);
        add_boundary_node(5, Node_pt[node_count]);

        // Increment the node number
        node_count++;

        // Now do the other columns (we extend to the rest of the nodes)
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Determine number of values
          local_node_num = l1 * n_p + l2 + n_p * n_p * (n_p - 1);

          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmax;

          // Add the node to boundary 5
          add_boundary_node(5, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }
      }


      // Now loop over the rest of the elements in the row, apart from the last
      for (unsigned j = 1; j < (Nx - 1); j++)
      {
        // Allocate memory for new element
        element_num = Nx * i + j + (Nz - 1) * Nx * Ny;
        Element_pt[element_num] = new ELEMENT;

        // The lowest nodes are copied from the lower element
        for (unsigned l1 = 0; l1 < n_p; l1++)
        {
          for (unsigned l2 = 0; l2 < n_p; l2++)
          {
            finite_element_pt(Nx * i + j + (Nz - 1) * Nx * Ny)
              ->node_pt(l2 + n_p * l1) =
              finite_element_pt(Nx * i + j + (Nz - 2) * Nx * Ny)
                ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
          }
        }

        // We extend to the third dimension but the last layer of nodes

        for (unsigned l3 = 1; l3 < (n_p - 1); l3++)
        {
          // The first row is copied from the lower element

          for (unsigned l2 = 0; l2 < n_p; l2++)
          {
            finite_element_pt(Nx * i + j + (Nz - 1) * Nx * Ny)
              ->node_pt(l2 + l3 * n_p * n_p) =
              finite_element_pt(Nx * (i - 1) + j + (Nz - 1) * Nx * Ny)
                ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
          }

          for (unsigned l1 = 1; l1 < n_p; l1++)
          {
            // First column is same as neighbouring element
            finite_element_pt(Nx * i + j + (Nz - 1) * Nx * Ny)
              ->node_pt(l1 * n_p + l3 * n_p * n_p) =
              finite_element_pt(Nx * i + (j - 1) + (Nz - 1) * Nx * Ny)
                ->node_pt(l1 * n_p + l3 * n_p * n_p + (n_p - 1));

            // New nodes for other columns
            for (unsigned l2 = 1; l2 < n_p; l2++)
            {
              // Determine number of values
              local_node_num = l1 * n_p + l2 + l3 * n_p * n_p;

              // Allocate memory for the nodes
              Node_pt[node_count] =
                finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
              // Set the pointer
              finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];
              // Get the fractional position of the node
              finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);


              // Set the position of the node
              Node_pt[node_count]->x(0) =
                Xmin + el_length[0] * (j + s_fraction[0]);
              Node_pt[node_count]->x(1) =
                Ymin + el_length[1] * (i + s_fraction[1]);
              Node_pt[node_count]->x(2) =
                Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);
              // No boundaries

              // Increment the node number
              node_count++;
            }
          }
        }

        // Now the top nodes

        // The first row is copied from the lower element

        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * i + j + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + (n_p - 1) * n_p * n_p) =
            finite_element_pt(Nx * (i - 1) + j + (Nz - 1) * Nx * Ny)
              ->node_pt((n_p - 1) * n_p + l2 + (n_p - 1) * n_p * n_p);
        }

        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(Nx * i + j + (Nz - 1) * Nx * Ny)
            ->node_pt(l1 * n_p + (n_p - 1) * n_p * n_p) =
            finite_element_pt(Nx * i + (j - 1) + (Nz - 1) * Nx * Ny)
              ->node_pt(l1 * n_p + (n_p - 1) * n_p * n_p + (n_p - 1));

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Determine number of values
            local_node_num = l1 * n_p + l2 + (n_p - 1) * n_p * n_p;

            // Allocate memory for the nodes
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) = Zmax;

            // Add to boundary 5
            add_boundary_node(5, Node_pt[node_count]);

            // Increment the node number
            node_count++;
          }
        }


      } // End of loop over elements in row


      // Do final element in row

      // Allocate memory for element
      element_num = Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < n_p; l1++)
      {
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + n_p * l1) =
            finite_element_pt(Nx * i + Nx - 1 + (Nz - 2) * Nx * Ny)
              ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
        }
      }


      // We go to the third dimension but we do not reach the top
      for (unsigned l3 = 1; l3 < (n_p - 1); l3++)
      {
        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + l3 * n_p * n_p) =
            finite_element_pt(Nx * (i - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
              ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
        }

        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt(l1 * n_p + l3 * n_p * n_p) =
            finite_element_pt(Nx * i + Nx - 2 + (Nz - 1) * Nx * Ny)
              ->node_pt(l1 * n_p + (n_p - 1) + l3 * n_p * n_p);

          // Middle nodes
          for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
          {
            // Determine number of values
            local_node_num = l1 * n_p + l2 + l3 * n_p * n_p;

            // Allocate memory for node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

            // No boundaries

            // Increment the node number
            node_count++;
          }

          // Final node

          // Determine number of values
          local_node_num = l1 * n_p + (n_p - 1) + l3 * n_p * n_p;

          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmax;
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // Add the node to the boundary 2
          add_boundary_node(2, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        } // End of loop over rows of nodes in the element


      } // End of the 3dimension loop

      // Now the top nodes

      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + (n_p - 1) * n_p * n_p) =
          finite_element_pt(Nx * (i - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt((n_p - 1) * n_p + l2 + (n_p - 1) * n_p * n_p);
      }

      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // First column is same as neighbouring element
        finite_element_pt(Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(l1 * n_p + (n_p - 1) * n_p * n_p) =
          finite_element_pt(Nx * i + Nx - 2 + (Nz - 1) * Nx * Ny)
            ->node_pt(l1 * n_p + (n_p - 1) + (n_p - 1) * n_p * n_p);

        // Middle nodes
        for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
        {
          // Determine number of values
          local_node_num = l1 * n_p + l2 + (n_p - 1) * n_p * n_p;

          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmax;

          // Add to boundary 5
          add_boundary_node(5, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }

        // Final node

        // Determine number of values
        local_node_num = l1 * n_p + (n_p - 1) + (n_p - 1) * n_p * n_p;

        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmax;

        // Add the node to the boundaries 2 and 5
        add_boundary_node(2, Node_pt[node_count]);
        add_boundary_node(5, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      } // End of loop over rows of nodes in the element


    } // End of loop over rows of elements


    // FINAL ELEMENT ROW (IN TOP  LAYER)
    //===========================================


    // FIRST ELEMENT IN UPPER ROW (upper left corner)
    //----------------------------------------------

    // Allocate memory for element
    element_num = Nx * (Ny - 1) + (Nz - 1) * Nx * Ny;
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < n_p; l1++)
    {
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(Nx * (Ny - 1) + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + n_p * l1) =
          finite_element_pt(Nx * (Ny - 1) + (Nz - 2) * Nx * Ny)
            ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
      }
    }

    // We jump to the third dimension but the last layer of nodes
    for (unsigned l3 = 1; l3 < (n_p - 1); l3++)
    {
      // The first row of nodes is copied from the element below
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(Nx * (Ny - 1) + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + l3 * n_p * n_p) =
          finite_element_pt(Nx * (Ny - 2) + (Nz - 1) * Nx * Ny)
            ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
      }

      // Second row of  nodes
      // First column of nodes
      for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
      {
        // Determine number of values
        local_node_num = n_p * l1 + l3 * n_p * n_p;

        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // Add the node to the boundary 4
        add_boundary_node(4, Node_pt[node_count]);

        // Increment the node number
        node_count++;

        // Now do the other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Determine number of values
          local_node_num = n_p * l1 + l2 + l3 * n_p * n_p;

          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // No boundaries

          // Increment the node number
          node_count++;
        }
      }

      // Final row of nodes
      // First column of nodes
      // Top left node
      // Determine number of values
      local_node_num = n_p * (n_p - 1) + l3 * n_p * n_p;
      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) =
        Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      // Add the node to the boundaries
      add_boundary_node(3, Node_pt[node_count]);
      add_boundary_node(4, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Now do the other columns
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Determine number of values
        local_node_num = n_p * (n_p - 1) + l2 + l3 * n_p * n_p;
        // Allocate memory for the node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // Add the node to the boundary 3
        add_boundary_node(3, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    }
    // Now the top nodes

    // The first row of nodes is copied from the element below
    for (unsigned l2 = 0; l2 < n_p; l2++)
    {
      finite_element_pt(Nx * (Ny - 1) + (Nz - 1) * Nx * Ny)
        ->node_pt(l2 + (n_p - 1) * n_p * n_p) =
        finite_element_pt(Nx * (Ny - 2) + (Nz - 1) * Nx * Ny)
          ->node_pt((n_p - 1) * n_p + l2 + (n_p - 1) * n_p * n_p);
    }

    // Second row of  nodes
    // First column of nodes
    for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
    {
      // Determine number of values
      local_node_num = n_p * l1 + (n_p - 1) * n_p * n_p;

      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) =
        Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundaries 4 and 5
      add_boundary_node(4, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Now do the other columns
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Determine number of values
        local_node_num = n_p * l1 + l2 + (n_p - 1) * n_p * n_p;

        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];
        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmax;

        // Add the node to the boundary 5
        add_boundary_node(5, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    }

    // Final row of nodes
    // First column of nodes
    // Top left node
    // Determine number of values
    local_node_num = n_p * (n_p - 1) + (n_p - 1) * n_p * n_p;
    // Allocate memory for node
    Node_pt[node_count] =
      finite_element_pt(element_num)
        ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
      Node_pt[node_count];
    // Get the fractional position of the node
    finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmax;

    // Add the node to the boundaries
    add_boundary_node(3, Node_pt[node_count]);
    add_boundary_node(4, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Now do the other columns
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // Determine number of values
      local_node_num = n_p * (n_p - 1) + l2 + (n_p - 1) * n_p * n_p;
      // Allocate memory for the node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundaries 3 and 5
      add_boundary_node(3, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }


    // Now loop over the rest of the elements in the row, apart from the last
    for (unsigned j = 1; j < (Nx - 1); j++)
    {
      // Allocate memory for element
      element_num = Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < n_p; l1++)
      {
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + n_p * l1) =
            finite_element_pt(Nx * (Ny - 1) + j + (Nz - 2) * Nx * Ny)
              ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
        }
      }

      // Jump in the third dimension but the top nodes

      for (unsigned l3 = 1; l3 < (n_p - 1); l3++)
      {
        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + l3 * n_p * n_p) =
            finite_element_pt(Nx * (Ny - 2) + j + (Nz - 1) * Nx * Ny)
              ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
        }

        // Second rows
        for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
            ->node_pt(n_p * l1 + l3 * n_p * n_p) =
            finite_element_pt(Nx * (Ny - 1) + (j - 1) + (Nz - 1) * Nx * Ny)
              ->node_pt(n_p * l1 + (n_p - 1) + l3 * n_p * n_p);

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Determine number of values
            local_node_num = n_p * l1 + l2 + l3 * n_p * n_p;
            // Allocate memory for the node
            Node_pt[node_count] =
              finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);

            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
            Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

            // No boundaries

            // Increment the node number
            // add_boundary_node(0,Node_pt[node_count]);
            node_count++;
          }
        }

        // Top row
        // First column is same as neighbouring element
        finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
          ->node_pt(n_p * (n_p - 1) + l3 * n_p * n_p) =
          finite_element_pt(Nx * (Ny - 1) + (j - 1) + (Nz - 1) * Nx * Ny)
            ->node_pt(n_p * (n_p - 1) + (n_p - 1) + l3 * n_p * n_p);
        // New nodes for other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Determine number of values
          local_node_num = n_p * (n_p - 1) + l2 + l3 * n_p * n_p;
          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymax;
          Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // Add the node to the boundary
          add_boundary_node(3, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }
      }

      // Now the top nodes

      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + (n_p - 1) * n_p * n_p) =
          finite_element_pt(Nx * (Ny - 2) + j + (Nz - 1) * Nx * Ny)
            ->node_pt((n_p - 1) * n_p + l2 + (n_p - 1) * n_p * n_p);
      }

      // Second rows
      for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
      {
        // First column is same as neighbouring element
        finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
          ->node_pt(n_p * l1 + (n_p - 1) * n_p * n_p) =
          finite_element_pt(Nx * (Ny - 1) + (j - 1) + (Nz - 1) * Nx * Ny)
            ->node_pt(n_p * l1 + (n_p - 1) + (n_p - 1) * n_p * n_p);

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Determine number of values
          local_node_num = n_p * l1 + l2 + (n_p - 1) * n_p * n_p;
          // Allocate memory for the node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);

          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmax;

          // Add the node to the boundary 5
          add_boundary_node(5, Node_pt[node_count]);

          // Increment the node number add_boundary_node(0,Node_pt[node_count]);
          node_count++;
        }
      }

      // Top row
      // First column is same as neighbouring element
      finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
        ->node_pt(n_p * (n_p - 1) + (n_p - 1) * n_p * n_p) =
        finite_element_pt(Nx * (Ny - 1) + (j - 1) + (Nz - 1) * Nx * Ny)
          ->node_pt(n_p * (n_p - 1) + (n_p - 1) + (n_p - 1) * n_p * n_p);
      // New nodes for other columns
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Determine number of values
        local_node_num = n_p * (n_p - 1) + l2 + (n_p - 1) * n_p * n_p;
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) = Zmax;

        // Add the node to the boundaries 3 and 5
        add_boundary_node(3, Node_pt[node_count]);
        add_boundary_node(5, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }


    } // End of loop over central elements in row


    // LAST ELEMENT (Really!!)
    //-----------------------------------------

    // Allocate memory for element
    element_num = Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny;
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < n_p; l1++)
    {
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + n_p * l1) =
          finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 2) * Nx * Ny)
            ->node_pt(l2 + n_p * l1 + (n_p - 1) * n_p * n_p);
      }
    }

    // We jump to the third dimension but the top nodes

    for (unsigned l3 = 1; l3 < (n_p - 1); l3++)
    {
      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + l3 * n_p * n_p) =
          finite_element_pt(Nx * (Ny - 2) + Nx - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt((n_p - 1) * n_p + l2 + l3 * n_p * n_p);
      }

      // Second rows
      for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
      {
        // First column is same as neighbouring element
        finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(n_p * l1 + l3 * n_p * n_p) =
          finite_element_pt(Nx * (Ny - 1) + Nx - 2 + (Nz - 1) * Nx * Ny)
            ->node_pt(n_p * l1 + (n_p - 1) + l3 * n_p * n_p);

        // Middle nodes
        for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
        {
          // Determine number of values
          local_node_num = n_p * l1 + l2 + l3 * n_p * n_p;
          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(element_num)
              ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // No boundaries

          // Increment the node number
          node_count++;
        }

        // Final node
        // Determine number of values
        local_node_num = n_p * l1 + (n_p - 1) + l3 * n_p * n_p;
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // Add the node to the boundary 2
        add_boundary_node(2, Node_pt[node_count]);


        // Increment the node number
        node_count++;

      } // End of loop over middle rows

      // Final row
      // First column is same as neighbouring element
      finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
        ->node_pt(n_p * (n_p - 1) + l3 * n_p * n_p) =
        finite_element_pt(Nx * (Ny - 1) + Nx - 2 + (Nz - 1) * Nx * Ny)
          ->node_pt(n_p * (n_p - 1) + (n_p - 1) + l3 * n_p * n_p);

      // Middle nodes
      for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
      {
        // Determine number of values
        local_node_num = n_p * (n_p - 1) + l2 + l3 * n_p * n_p;
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // Add the node to the boundary 3
        add_boundary_node(3, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }

      // Final node
      // Determine number of values
      local_node_num = n_p * (n_p - 1) + (n_p - 1) + l3 * n_p * n_p;
      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      // In fluid 2
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) =
        Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      // Add the node to the boundaries 2 and 3
      add_boundary_node(2, Node_pt[node_count]);
      add_boundary_node(3, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }


    // Now the top nodes


    // The first row is copied from the lower element
    for (unsigned l2 = 0; l2 < n_p; l2++)
    {
      finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
        ->node_pt(l2 + (n_p - 1) * n_p * n_p) =
        finite_element_pt(Nx * (Ny - 2) + Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt((n_p - 1) * n_p + l2 + (n_p - 1) * n_p * n_p);
    }

    // Second rows
    for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
    {
      // First column is same as neighbouring element
      finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
        ->node_pt(n_p * l1 + (n_p - 1) * n_p * n_p) =
        finite_element_pt(Nx * (Ny - 1) + Nx - 2 + (Nz - 1) * Nx * Ny)
          ->node_pt(n_p * l1 + (n_p - 1) + (n_p - 1) * n_p * n_p);

      // Middle nodes
      for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
      {
        // Determine number of values
        local_node_num = n_p * l1 + l2 + (n_p - 1) * n_p * n_p;
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmax;

        // Add to boundary 5
        add_boundary_node(5, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }

      // Final node
      // Determine number of values
      local_node_num = n_p * l1 + (n_p - 1) + (n_p - 1) * n_p * n_p;
      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) =
        Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundaries 2 and 5
      add_boundary_node(2, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);

      // Increment the node number
      node_count++;

    } // End of loop over middle rows

    // Final row
    // First node is same as neighbouring element
    finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
      ->node_pt(n_p * (n_p - 1) + (n_p - 1) * n_p * n_p) =
      finite_element_pt(Nx * (Ny - 1) + Nx - 2 + (Nz - 1) * Nx * Ny)
        ->node_pt(n_p * (n_p - 1) + (n_p - 1) + (n_p - 1) * n_p * n_p);

    // Middle nodes
    for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
    {
      // Determine number of values
      local_node_num = n_p * (n_p - 1) + l2 + (n_p - 1) * n_p * n_p;
      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
        Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundary 3
      add_boundary_node(3, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);


      // Increment the node number
      node_count++;
    }

    // Final node (really!!)
    // Determine number of values
    local_node_num = n_p * (n_p - 1) + (n_p - 1) + (n_p - 1) * n_p * n_p;
    // Allocate memory for node
    Node_pt[node_count] =
      finite_element_pt(element_num)
        ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer
    finite_element_pt(element_num)->node_pt(local_node_num) =
      Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmax;

    // Add the node to the boundaries 2, 3 and 5
    add_boundary_node(2, Node_pt[node_count]);
    add_boundary_node(3, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);

    // Increment the node number
    node_count++;


    // Setup lookup scheme that establishes which elements are located
    // on the mesh boundaries
    setup_boundary_element_info();
  }


} // namespace oomph

#endif
