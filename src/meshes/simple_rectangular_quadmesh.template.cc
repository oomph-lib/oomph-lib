// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_SIMPLE_RECTANGULAR_QUADMESH_TEMPLATE_CC
#define OOMPH_SIMPLE_RECTANGULAR_QUADMESH_TEMPLATE_CC

#include "../generic/Qelements.h"
#include "simple_rectangular_quadmesh.template.h"


namespace oomph
{
  //=======================================================================
  /// Constructor for 2D Quad mesh class:
  ///
  /// Nx  : number of elements in the x direction
  ///
  /// Ny  : number of elements in the y direction
  ///
  /// Lx  : length in the x direction
  ///
  /// Ly  : length in the y direction
  ///
  /// Ordering of elements: 'Lower left' to 'lower right' then 'upwards'
  ///
  /// Timestepper defaults to Steady.
  //=======================================================================
  template<class ELEMENT>
  SimpleRectangularQuadMesh<ELEMENT>::SimpleRectangularQuadMesh(
    const unsigned& Nx,
    const unsigned& Ny,
    const double& Lx,
    const double& Ly,
    TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Set the internal values
    NX = Nx;
    NY = Ny;

    // Set the number of boundaries
    set_nboundary(4);

    // Allocate the store for the elements
    Element_pt.resize(Nx * Ny);

    // Create first element
    Element_pt[0] = new ELEMENT;

    // Read out the number of linear points in the element
    unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

    // Can now allocate the store for the nodes
    Node_pt.resize((1 + (n_p - 1) * Nx) * (1 + (n_p - 1) * Ny));

    // Set up geometrical data
    //------------------------

    unsigned long node_count = 0;
    double xinit = 0.0, yinit = 0.0;

    // Set the values of the increments
    // double xstep = Lx/((n_p-1)*Nx);
    // double ystep = Ly/((n_p-1)*Ny);

    // Set the length of the element
    double el_length_x = Lx / double(Nx);
    double el_length_y = Ly / double(Ny);

    // Storage for local coordinate in element
    Vector<double> s_fraction;

    // Now assign the topology
    // Boundaries are numbered 0 1 2 3 from the bottom proceeding anticlockwise
    // Pinned value are denoted by an integer value 1
    // Thus if a node is on two boundaries, ORing the values of the
    // boundary conditions will give the most restrictive case (pinning)


    // FIRST ELEMENT (lower left corner)
    //----------------------------------

    // Set the corner node

    // Allocate memory for the node
    Node_pt[node_count] =
      finite_element_pt(0)->construct_boundary_node(0, time_stepper_pt);

    // Set the pointer from the element
    finite_element_pt(0)->node_pt(0) = Node_pt[node_count];

    // Set the position of the node
    Node_pt[node_count]->x(0) = xinit;
    Node_pt[node_count]->x(1) = yinit;

    // Add the node to boundaries 0 and 3
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(3, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Loop over the other nodes in the first row
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(0)->construct_boundary_node(l2, time_stepper_pt);

      // Set the pointer from the element
      finite_element_pt(0)->node_pt(l2) = Node_pt[node_count];

      // Get the local fraction of the node
      finite_element_pt(0)->local_fraction_of_node(l2, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = xinit + el_length_x * s_fraction[0];
      Node_pt[node_count]->x(1) = yinit;

      // Add the node to the boundary
      add_boundary_node(0, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // Loop over the other node columns
    for (unsigned l1 = 1; l1 < n_p; l1++)
    {
      // Allocate memory for the nodes
      Node_pt[node_count] = finite_element_pt(0)->construct_boundary_node(
        l1 * n_p, time_stepper_pt);

      // Set the pointer from the element
      finite_element_pt(0)->node_pt(l1 * n_p) = Node_pt[node_count];

      // Get the fractional position
      finite_element_pt(0)->local_fraction_of_node(l1 * n_p, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = xinit;
      Node_pt[node_count]->x(1) = yinit + el_length_y * s_fraction[1];

      // Add the node to the boundary
      add_boundary_node(3, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Loop over the other nodes in the row
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Allocate the memory for the node
        Node_pt[node_count] =
          finite_element_pt(0)->construct_node(l1 * n_p + l2, time_stepper_pt);

        // Set the pointer from the element
        finite_element_pt(0)->node_pt(l1 * n_p + l2) = Node_pt[node_count];

        // Get the fractional position
        finite_element_pt(0)->local_fraction_of_node(l1 * n_p + l2, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = xinit + el_length_x * s_fraction[0];
        Node_pt[node_count]->x(1) = yinit + el_length_y * s_fraction[1];

        // Increment the node number
        node_count++;
      }
    }


    // CENTRE OF FIRST ROW OF ELEMENTS
    //--------------------------------
    // Now loop over the first row of elements, apart from final element
    for (unsigned j = 1; j < (Nx - 1); j++)
    {
      // Allocate memory for new element
      Element_pt[j] = new ELEMENT;

      // Do first row of nodes

      // First column of nodes is same as neighbouring element
      finite_element_pt(j)->node_pt(0) =
        finite_element_pt(j - 1)->node_pt((n_p - 1));

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(j)->construct_boundary_node(l2, time_stepper_pt);

        // Set the pointer from the element
        finite_element_pt(j)->node_pt(l2) = Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(j)->local_fraction_of_node(l2, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = xinit + el_length_x * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = yinit;

        // Add the node to the boundary
        add_boundary_node(0, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }

      // Do the rest of the nodes
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // First column of nodes is same as neighbouring element
        finite_element_pt(j)->node_pt(l1 * n_p) =
          finite_element_pt(j - 1)->node_pt(l1 * n_p + (n_p - 1));

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Allocate memory for the nodes
          Node_pt[node_count] = finite_element_pt(j)->construct_node(
            l1 * n_p + l2, time_stepper_pt);

          // Set the pointer from the element
          finite_element_pt(j)->node_pt(l1 * n_p + l2) = Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(j)->local_fraction_of_node(l1 * n_p + l2,
                                                       s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = xinit + el_length_x * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = yinit + el_length_y * s_fraction[1];

          // Increment the node number
          node_count++;
        }
      }
    }


    // FINAL ELEMENT IN FIRST ROW (lower right corner)
    //-----------------------------------------------

    // Allocate memory for element
    Element_pt[Nx - 1] = new ELEMENT;
    // First column of nodes is same as neighbouring element
    finite_element_pt(Nx - 1)->node_pt(0) =
      finite_element_pt(Nx - 2)->node_pt(n_p - 1);

    // New middle nodes
    for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
    {
      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(Nx - 1)->construct_boundary_node(l2, time_stepper_pt);

      // Set the pointer from the element
      finite_element_pt(Nx - 1)->node_pt(l2) = Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(Nx - 1)->local_fraction_of_node(l2, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
        xinit + el_length_x * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = yinit;

      // Add the node to the boundary
      add_boundary_node(0, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // New final node

    // Allocate memory for the node
    Node_pt[node_count] = finite_element_pt(Nx - 1)->construct_boundary_node(
      n_p - 1, time_stepper_pt);

    // Set the pointer from the element
    finite_element_pt(Nx - 1)->node_pt(n_p - 1) = Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(Nx - 1)->local_fraction_of_node(n_p - 1, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = xinit + el_length_x * (Nx - 1 + s_fraction[0]);
    Node_pt[node_count]->x(1) = yinit;

    // Add the node to the boundaries
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(1, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Do the rest of the nodes
    for (unsigned l1 = 1; l1 < n_p; l1++)
    {
      // First column of nodes is same as neighbouring element
      finite_element_pt(Nx - 1)->node_pt(l1 * n_p) =
        finite_element_pt(Nx - 2)->node_pt(l1 * n_p + (n_p - 1));

      // New node for middle column
      for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
      {
        // Allocate memory for node
        Node_pt[node_count] = finite_element_pt(Nx - 1)->construct_node(
          l1 * n_p + l2, time_stepper_pt);

        // Set the pointer from the element
        finite_element_pt(Nx - 1)->node_pt(l1 * n_p + l2) = Node_pt[node_count];

        // Get the fractional position
        finite_element_pt(Nx - 1)->local_fraction_of_node(l1 * n_p + l2,
                                                          s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          xinit + el_length_x * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = yinit + el_length_y * s_fraction[1];

        // Increment the node number
        node_count++;
      }

      // New node for final column

      // Allocate memory for node
      Node_pt[node_count] = finite_element_pt(Nx - 1)->construct_boundary_node(
        l1 * n_p + (n_p - 1), time_stepper_pt);

      // Set the pointer from the element
      finite_element_pt(Nx - 1)->node_pt(l1 * n_p + (n_p - 1)) =
        Node_pt[node_count];

      // Get the fractional position
      finite_element_pt(Nx - 1)->local_fraction_of_node(l1 * n_p + (n_p - 1),
                                                        s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
        xinit + el_length_x * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = yinit + el_length_y * s_fraction[1];

      // Add the node to the boundary
      add_boundary_node(1, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }


    // ALL CENTRAL ELEMENT ROWS
    //------------------------

    // Loop over remaining element rows
    for (unsigned i = 1; i < (Ny - 1); i++)
    {
      // Set the first element in the row

      // Allocate memory for element
      Element_pt[Nx * i] = new ELEMENT;

      // The first row of nodes is copied from the element below
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(Nx * i)->node_pt(l2) =
          finite_element_pt(Nx * (i - 1))->node_pt((n_p - 1) * n_p + l2);
      }

      // Other rows are new nodes
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // First column of nodes

        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(Nx * i)->construct_boundary_node(l1 * n_p,
                                                             time_stepper_pt);

        // Set the pointer from the element
        finite_element_pt(Nx * i)->node_pt(l1 * n_p) = Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(Nx * i)->local_fraction_of_node(l1 * n_p, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = xinit;
        Node_pt[node_count]->x(1) = yinit + el_length_y * (i + s_fraction[1]);

        // Add the node to the boundary
        add_boundary_node(3, Node_pt[node_count]);

        // Increment the node number
        node_count++;

        // Now do the other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Allocate memory for node
          Node_pt[node_count] = finite_element_pt(Nx * i)->construct_node(
            l1 * n_p + l2, time_stepper_pt);

          // Set the pointer from the element
          finite_element_pt(Nx * i)->node_pt(l1 * n_p + l2) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(Nx * i)->local_fraction_of_node(l1 * n_p + l2,
                                                            s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = xinit + el_length_x * s_fraction[0];
          Node_pt[node_count]->x(1) = yinit + el_length_y * (i + s_fraction[1]);

          // Increment the node number
          node_count++;
        }
      }

      // Now loop over the rest of the elements in the row, apart from the last
      for (unsigned j = 1; j < (Nx - 1); j++)
      {
        // Allocate memory for new element
        Element_pt[Nx * i + j] = new ELEMENT;

        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < n_p; l2++)
        {
          finite_element_pt(Nx * i + j)->node_pt(l2) =
            finite_element_pt(Nx * (i - 1) + j)->node_pt((n_p - 1) * n_p + l2);
        }

        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(Nx * i + j)->node_pt(l1 * n_p) =
            finite_element_pt(Nx * i + (j - 1))->node_pt(l1 * n_p + (n_p - 1));

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < n_p; l2++)
          {
            // Allocate memory for the nodes
            Node_pt[node_count] =
              finite_element_pt(Nx * i + j)
                ->construct_node(l1 * n_p + l2, time_stepper_pt);

            // Set the pointer
            finite_element_pt(Nx * i + j)->node_pt(l1 * n_p + l2) =
              Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(Nx * i + j)
              ->local_fraction_of_node(l1 * n_p + l2, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              xinit + el_length_x * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) =
              yinit + el_length_y * (i + s_fraction[1]);

            // Increment the node number
            node_count++;
          }
        }
      } // End of loop over elements in row

      // Do final element in row

      // Allocate memory for element
      Element_pt[Nx * i + Nx - 1] = new ELEMENT;

      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(Nx * i + Nx - 1)->node_pt(l2) =
          finite_element_pt(Nx * (i - 1) + Nx - 1)
            ->node_pt((n_p - 1) * n_p + l2);
      }

      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // First column is same as neighbouring element
        finite_element_pt(Nx * i + Nx - 1)->node_pt(l1 * n_p) =
          finite_element_pt(Nx * i + Nx - 2)->node_pt(l1 * n_p + (n_p - 1));

        // Middle nodes
        for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
        {
          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(Nx * i + Nx - 1)
              ->construct_node(l1 * n_p + l2, time_stepper_pt);

          // Set the pointer
          finite_element_pt(Nx * i + Nx - 1)->node_pt(l1 * n_p + l2) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(Nx * i + Nx - 1)
            ->local_fraction_of_node(l1 * n_p + l2, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
            xinit + el_length_x * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) = yinit + el_length_y * (i + s_fraction[1]);

          // Increment the node number
          node_count++;
        }

        // Final node

        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(Nx * i + Nx - 1)
            ->construct_boundary_node(l1 * n_p + (n_p - 1), time_stepper_pt);

        // Set the pointer
        finite_element_pt(Nx * i + Nx - 1)->node_pt(l1 * n_p + (n_p - 1)) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(Nx * i + Nx - 1)
          ->local_fraction_of_node(l1 * n_p + (n_p - 1), s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          xinit + el_length_x * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = yinit + el_length_y * (i + s_fraction[1]);

        // Add the node to the boundary
        add_boundary_node(1, Node_pt[node_count]);

        // Increment the node number
        node_count++;

      } // End of loop over rows of nodes in the element
    } // End of loop over rows of elements


    // FINAL ELEMENT ROW
    //=================


    // FIRST ELEMENT IN UPPER ROW (upper left corner)
    //----------------------------------------------

    // Allocate memory for element
    Element_pt[Nx * (Ny - 1)] = new ELEMENT;

    // The first row of nodes is copied from the element below
    for (unsigned l2 = 0; l2 < n_p; l2++)
    {
      finite_element_pt(Nx * (Ny - 1))->node_pt(l2) =
        finite_element_pt(Nx * (Ny - 2))->node_pt((n_p - 1) * n_p + l2);
    }

    // Second row of  nodes
    // First column of nodes
    for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
    {
      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(Nx * (Ny - 1))
          ->construct_boundary_node(n_p * l1, time_stepper_pt);

      // Set the pointer from the element
      finite_element_pt(Nx * (Ny - 1))->node_pt(n_p * l1) = Node_pt[node_count];

      // Get the fractional position of the element
      finite_element_pt(Nx * (Ny - 1))
        ->local_fraction_of_node(n_p * l1, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = xinit;
      Node_pt[node_count]->x(1) =
        yinit + el_length_y * (Ny - 1 + s_fraction[1]);

      // Add the node to the boundary
      add_boundary_node(3, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Now do the other columns
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(Nx * (Ny - 1))
            ->construct_node(n_p * l1 + l2, time_stepper_pt);

        // Set the pointer from the element
        finite_element_pt(Nx * (Ny - 1))->node_pt(n_p * l1 + l2) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(Nx * (Ny - 1))
          ->local_fraction_of_node(n_p * l1 + l2, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = xinit + el_length_x * s_fraction[0];
        Node_pt[node_count]->x(1) =
          yinit + el_length_y * (Ny - 1 + s_fraction[1]);

        // Increment the node number
        node_count++;
      }
    }

    // Final row of nodes
    // First column of nodes
    // Top left node

    // Allocate memory for node
    Node_pt[node_count] =
      finite_element_pt(Nx * (Ny - 1))
        ->construct_boundary_node(n_p * (n_p - 1), time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(Nx * (Ny - 1))->node_pt(n_p * (n_p - 1)) =
      Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(Nx * (Ny - 1))
      ->local_fraction_of_node(n_p * (n_p - 1), s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = xinit;
    Node_pt[node_count]->x(1) = yinit + el_length_y * Ny;

    // Add the node to the boundaries
    add_boundary_node(2, Node_pt[node_count]);
    add_boundary_node(3, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Now do the other columns
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // Allocate memory for the node
      Node_pt[node_count] =
        finite_element_pt(Nx * (Ny - 1))
          ->construct_boundary_node(n_p * (n_p - 1) + l2, time_stepper_pt);

      // Set the pointer from the element
      finite_element_pt(Nx * (Ny - 1))->node_pt(n_p * (n_p - 1) + l2) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(Nx * (Ny - 1))
        ->local_fraction_of_node(n_p * (n_p - 1) + l2, s_fraction);


      // Set the position of the node
      Node_pt[node_count]->x(0) = xinit + el_length_x * s_fraction[0];
      Node_pt[node_count]->x(1) = yinit + el_length_y * Ny;

      // Add the node to the boundary
      add_boundary_node(2, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // Now loop over the rest of the elements in the row, apart from the last
    for (unsigned j = 1; j < (Nx - 1); j++)
    {
      // Allocate memory for element
      Element_pt[Nx * (Ny - 1) + j] = new ELEMENT;
      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < n_p; l2++)
      {
        finite_element_pt(Nx * (Ny - 1) + j)->node_pt(l2) =
          finite_element_pt(Nx * (Ny - 2) + j)->node_pt((n_p - 1) * n_p + l2);
      }

      // Second rows
      for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
      {
        // First column is same as neighbouring element
        finite_element_pt(Nx * (Ny - 1) + j)->node_pt(n_p * l1) =
          finite_element_pt(Nx * (Ny - 1) + (j - 1))
            ->node_pt(n_p * l1 + (n_p - 1));

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < n_p; l2++)
        {
          // Allocate memory for the node
          Node_pt[node_count] =
            finite_element_pt(Nx * (Ny - 1) + j)
              ->construct_node(n_p * l1 + l2, time_stepper_pt);
          // Set the pointer
          finite_element_pt(Nx * (Ny - 1) + j)->node_pt(n_p * l1 + l2) =
            Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(Nx * (Ny - 1) + j)
            ->local_fraction_of_node(n_p * l1 + l2, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = xinit + el_length_x * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) =
            yinit + el_length_y * (Ny - 1 + s_fraction[1]);

          // Increment the node number
          node_count++;
        }
      }

      // Top row
      // First column is same as neighbouring element
      finite_element_pt(Nx * (Ny - 1) + j)->node_pt(n_p * (n_p - 1)) =
        finite_element_pt(Nx * (Ny - 1) + (j - 1))
          ->node_pt(n_p * (n_p - 1) + (n_p - 1));
      // New nodes for other columns
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(Nx * (Ny - 1) + j)
            ->construct_boundary_node(n_p * (n_p - 1) + l2, time_stepper_pt);
        // Set the pointer
        finite_element_pt(Nx * (Ny - 1) + j)->node_pt(n_p * (n_p - 1) + l2) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(Nx * (Ny - 1) + j)
          ->local_fraction_of_node(n_p * (n_p - 1) + l2, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = xinit + el_length_x * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = yinit + el_length_y * Ny;

        // Add the node to the boundary
        add_boundary_node(2, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    } // End of loop over central elements in row


    // FINAL ELEMENT IN ROW (upper right corner)
    //-----------------------------------------

    // Allocate memory for element
    Element_pt[Nx * (Ny - 1) + Nx - 1] = new ELEMENT;
    // The first row is copied from the lower element
    for (unsigned l2 = 0; l2 < n_p; l2++)
    {
      finite_element_pt(Nx * (Ny - 1) + Nx - 1)->node_pt(l2) =
        finite_element_pt(Nx * (Ny - 2) + Nx - 1)
          ->node_pt((n_p - 1) * n_p + l2);
    }

    // Second rows
    for (unsigned l1 = 1; l1 < (n_p - 1); l1++)
    {
      // First column is same as neighbouring element
      finite_element_pt(Nx * (Ny - 1) + Nx - 1)->node_pt(n_p * l1) =
        finite_element_pt(Nx * (Ny - 1) + Nx - 2)
          ->node_pt(n_p * l1 + (n_p - 1));

      // Middle nodes
      for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
      {
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(Nx * (Ny - 1) + Nx - 1)
            ->construct_node(n_p * l1 + l2, time_stepper_pt);
        // Set the pointer
        finite_element_pt(Nx * (Ny - 1) + Nx - 1)->node_pt(n_p * l1 + l2) =
          Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(Nx * (Ny - 1) + Nx - 1)
          ->local_fraction_of_node(n_p * l1 + l2, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          xinit + el_length_x * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) =
          yinit + el_length_y * (Ny - 1 + s_fraction[1]);

        // Increment the node number
        node_count++;
      }

      // Final node

      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(Nx * (Ny - 1) + Nx - 1)
          ->construct_boundary_node(n_p * l1 + (n_p - 1), time_stepper_pt);
      // Set the pointer
      finite_element_pt(Nx * (Ny - 1) + Nx - 1)->node_pt(n_p * l1 + (n_p - 1)) =
        Node_pt[node_count];

      // Get the fractional position
      finite_element_pt(Nx * (Ny - 1) + Nx - 1)
        ->local_fraction_of_node(n_p * l1 + (n_p - 1), s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = xinit + el_length_x * Nx;
      Node_pt[node_count]->x(1) =
        yinit + el_length_y * (Ny - 1 + s_fraction[1]);

      // Add the node to the boundary
      add_boundary_node(1, Node_pt[node_count]);

      // Increment the node number
      node_count++;

    } // End of loop over middle rows

    // Final row
    // First column is same as neighbouring element
    finite_element_pt(Nx * (Ny - 1) + Nx - 1)->node_pt(n_p * (n_p - 1)) =
      finite_element_pt(Nx * (Ny - 1) + Nx - 2)
        ->node_pt(n_p * (n_p - 1) + (n_p - 1));

    // Middle nodes
    for (unsigned l2 = 1; l2 < (n_p - 1); l2++)
    {
      // Allocate memory for node
      Node_pt[node_count] =
        finite_element_pt(Nx * (Ny - 1) + Nx - 1)
          ->construct_boundary_node(n_p * (n_p - 1) + l2, time_stepper_pt);
      // Set the pointer
      finite_element_pt(Nx * (Ny - 1) + Nx - 1)->node_pt(n_p * (n_p - 1) + l2) =
        Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(Nx * (Ny - 1) + Nx - 1)
        ->local_fraction_of_node(n_p * (n_p - 1) + l2, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
        xinit + el_length_x * (Nx - 1 + s_fraction[0]);
      // In fluid 2
      Node_pt[node_count]->x(1) = yinit + el_length_y * Ny;

      // Add the node to the boundary
      add_boundary_node(2, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // Final node

    // Allocate memory for node
    Node_pt[node_count] =
      finite_element_pt(Nx * (Ny - 1) + Nx - 1)
        ->construct_boundary_node(n_p * (n_p - 1) + (n_p - 1), time_stepper_pt);
    // Set the pointer
    finite_element_pt(Nx * (Ny - 1) + Nx - 1)
      ->node_pt(n_p * (n_p - 1) + (n_p - 1)) = Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(Nx * (Ny - 1) + Nx - 1)
      ->local_fraction_of_node(n_p * (n_p - 1) + (n_p - 1), s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = xinit + el_length_x * Nx;
    Node_pt[node_count]->x(1) = yinit + el_length_y * Ny;

    // Add the node to the boundaries
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(2, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Setup lookup scheme that establishes which elements are located
    // on the mesh boundaries
    setup_boundary_element_info();
  }

} // namespace oomph
#endif
