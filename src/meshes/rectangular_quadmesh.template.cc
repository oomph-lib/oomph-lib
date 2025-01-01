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
#ifndef OOMPH_RECTANGULAR_QUADMESH_TEMPLATE_CC
#define OOMPH_RECTANGULAR_QUADMESH_TEMPLATE_CC


// OOMPH-LIB Headers
#include "rectangular_quadmesh.template.h"


namespace oomph
{
  //===========================================================================
  /// Generic mesh construction. This function contains the "guts" of the
  /// mesh generation process, including all the tedious loops, counting
  /// and spacing functions. The function should be called in all constuctors
  /// of any derived classes.
  //===========================================================================
  template<class ELEMENT>
  void RectangularQuadMesh<ELEMENT>::build_mesh(TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Set the number of boundaries
    set_nboundary(4);

    // Allocate the store for the elements
    Element_pt.resize(Nx * Ny);
    // Allocate the memory for the first element
    Element_pt[0] = new ELEMENT;

    // Read out the number of linear points in the element
    Np = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

    // Can now allocate the store for the nodes
    Node_pt.resize((1 + (Np - 1) * Nx) * (1 + (Np - 1) * Ny));

    // Set up geometrical data
    unsigned long node_count = 0;

    // Now assign the topology
    // Boundaries are numbered 0 1 2 3 from the bottom proceeding anticlockwise
    // Pinned value are denoted by an integer value 1
    // Thus if a node is on two boundaries, ORing the values of the
    // boundary conditions will give the most restrictive case (pinning)

    // FIRST ELEMENT (lower left corner)

    // Set the corner node
    // Allocate memory for the node
    Node_pt[node_count] =
      finite_element_pt(0)->construct_boundary_node(0, time_stepper_pt);

    // Set the position of the node
    Node_pt[node_count]->x(0) = x_spacing_function(0, 0, 0, 0);
    Node_pt[node_count]->x(1) = y_spacing_function(0, 0, 0, 0);

    // Push the node back onto boundaries
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(3, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Loop over the other nodes in the first row
    for (unsigned l2 = 1; l2 < Np; l2++)
    {
      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(0)->construct_boundary_node(l2, time_stepper_pt);

      // Set the position of the node
      Node_pt[node_count]->x(0) = x_spacing_function(0, l2, 0, 0);
      Node_pt[node_count]->x(1) = y_spacing_function(0, l2, 0, 0);

      // Push the node back onto boundaries
      add_boundary_node(0, Node_pt[node_count]);

      // If we only have one column then the RHS node is on the right boundary
      if ((Nx == 1) && (l2 == (Np - 1)))
      {
        add_boundary_node(1, Node_pt[node_count]);
      }

      // Increment the node number
      node_count++;
    }

    // Loop over the other node columns
    for (unsigned l1 = 1; l1 < Np; l1++)
    {
      // Allocate memory for the nodes
      Node_pt[node_count] =
        finite_element_pt(0)->construct_boundary_node(l1 * Np, time_stepper_pt);

      // Set the position of the node
      Node_pt[node_count]->x(0) = x_spacing_function(0, 0, 0, l1);
      Node_pt[node_count]->x(1) = y_spacing_function(0, 0, 0, l1);

      // Push the node back onto boundaries
      add_boundary_node(3, Node_pt[node_count]);

      // If we only have one row, then the top node is on the top boundary
      if ((Ny == 1) && (l1 == (Np - 1)))
      {
        add_boundary_node(2, Node_pt[node_count]);
      }

      // Increment the node number
      node_count++;

      // Loop over the other nodes in the row
      for (unsigned l2 = 1; l2 < Np; l2++)
      {
        // Allocate the memory for the node
        // If it lies on a boundary make a boundary node
        if (((Nx == 1) && (l2 == (Np - 1))) || ((Ny == 1) && (l1 == (Np - 1))))
        {
          Node_pt[node_count] = finite_element_pt(0)->construct_boundary_node(
            l1 * Np + l2, time_stepper_pt);
        }
        // Otherwise its a normal node
        else
        {
          Node_pt[node_count] =
            finite_element_pt(0)->construct_node(l1 * Np + l2, time_stepper_pt);
        }

        // Set the position of the node
        Node_pt[node_count]->x(0) = x_spacing_function(0, l2, 0, l1);
        Node_pt[node_count]->x(1) = y_spacing_function(0, l2, 0, l1);

        // If we only have one column then the RHS node is on the right boundary
        if ((Nx == 1) && l2 == (Np - 1))
        {
          add_boundary_node(1, Node_pt[node_count]);
        }
        // If we only have one row, then the top node is on the top boundary
        if ((Ny == 1) && (l1 == (Np - 1)))
        {
          add_boundary_node(2, Node_pt[node_count]);
        }

        // Increment the node number
        node_count++;
      }
    }
    // END OF FIRST ELEMENT

    // CENTRE OF FIRST ROW OF ELEMENTS
    // Now loop over the first row of elements, apart from final element
    for (unsigned j = 1; j < (Nx - 1); j++)
    {
      // Allocate memory for new element
      Element_pt[j] = new ELEMENT;
      // Do first row of nodes
      // First column of nodes is same as neighbouring element
      finite_element_pt(j)->node_pt(0) =
        finite_element_pt(j - 1)->node_pt((Np - 1));
      // New nodes for other columns
      for (unsigned l2 = 1; l2 < Np; l2++)
      {
        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(j)->construct_boundary_node(l2, time_stepper_pt);

        // Set the position of the node
        Node_pt[node_count]->x(0) = x_spacing_function(j, l2, 0, 0);
        Node_pt[node_count]->x(1) = y_spacing_function(j, l2, 0, 0);

        // Push the node back onto boundaries
        add_boundary_node(0, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }

      // Do the rest of the nodes
      for (unsigned l1 = 1; l1 < Np; l1++)
      {
        // First column of nodes is same as neighbouring element
        finite_element_pt(j)->node_pt(l1 * Np) =
          finite_element_pt(j - 1)->node_pt(l1 * Np + (Np - 1));
        // New nodes for other columns
        for (unsigned l2 = 1; l2 < Np; l2++)
        {
          // Allocate memory for the nodes
          // If we only have one row, this node could be on the boundary
          if ((Ny == 1) && (l1 == (Np - 1)))
          {
            Node_pt[node_count] = finite_element_pt(j)->construct_boundary_node(
              l1 * Np + l2, time_stepper_pt);
          }
          // Otherwise create a normal node
          else
          {
            Node_pt[node_count] = finite_element_pt(j)->construct_node(
              l1 * Np + l2, time_stepper_pt);
          }

          // Set the position of the node
          Node_pt[node_count]->x(0) = x_spacing_function(j, l2, 0, l1);
          Node_pt[node_count]->x(1) = y_spacing_function(j, l2, 0, l1);

          // If we only have one row, then the top node is on the top boundary
          if ((Ny == 1) && (l1 == (Np - 1)))
          {
            add_boundary_node(2, Node_pt[node_count]);
          }

          // Increment the node number
          node_count++;
        }
      }
    }
    // END OF CENTRE OF FIRST ROW OF ELEMENTS

    // FINAL ELEMENT IN FIRST ROW (lower right corner)
    // Only allocate if there is more than one element in the row
    if (Nx > 1)
    {
      // Allocate memory for element
      Element_pt[Nx - 1] = new ELEMENT;

      // First column of nodes is same as neighbouring element
      finite_element_pt(Nx - 1)->node_pt(0) =
        finite_element_pt(Nx - 2)->node_pt(Np - 1);

      // New middle nodes
      for (unsigned l2 = 1; l2 < (Np - 1); l2++)
      {
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(Nx - 1)->construct_boundary_node(l2,
                                                             time_stepper_pt);

        // Set the position of the node
        Node_pt[node_count]->x(0) = x_spacing_function(Nx - 1, l2, 0, 0);
        Node_pt[node_count]->x(1) = y_spacing_function(Nx - 1, l2, 0, 0);

        // Push the node back onto boundaries
        add_boundary_node(0, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }

      // Allocate memory for a new boundary node
      Node_pt[node_count] = finite_element_pt(Nx - 1)->construct_boundary_node(
        Np - 1, time_stepper_pt);

      // If required make it periodic from the node on the other side
      if (Xperiodic == true)
      {
        Node_pt[node_count]->make_periodic(finite_element_pt(0)->node_pt(0));
      }

      // Set the position of the node
      Node_pt[node_count]->x(0) = x_spacing_function(Nx - 1, Np - 1, 0, 0);
      Node_pt[node_count]->x(1) = y_spacing_function(Nx - 1, Np - 1, 0, 0);

      // Push the node back onto boundaries
      add_boundary_node(0, Node_pt[node_count]);
      add_boundary_node(1, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Do the rest of the nodes
      for (unsigned l1 = 1; l1 < Np; l1++)
      {
        // First column of nodes is same as neighbouring element
        finite_element_pt(Nx - 1)->node_pt(l1 * Np) =
          finite_element_pt(Nx - 2)->node_pt(l1 * Np + (Np - 1));

        // New node for middle column
        for (unsigned l2 = 1; l2 < (Np - 1); l2++)
        {
          // Allocate memory for node
          // If it lies on the boundary, create a boundary node
          if ((Ny == 1) && (l1 == (Np - 1)))
          {
            Node_pt[node_count] =
              finite_element_pt(Nx - 1)->construct_boundary_node(
                l1 * Np + l2, time_stepper_pt);
          }
          // Otherwise create a normal node
          else
          {
            Node_pt[node_count] = finite_element_pt(Nx - 1)->construct_node(
              l1 * Np + l2, time_stepper_pt);
          }

          // Set the position of the node
          Node_pt[node_count]->x(0) = x_spacing_function(Nx - 1, l2, 0, l1);
          Node_pt[node_count]->x(1) = y_spacing_function(Nx - 1, l2, 0, l1);

          // If we only have one row, then the top node is on the top boundary
          if ((Ny == 1) && (l1 == (Np - 1)))
          {
            add_boundary_node(2, Node_pt[node_count]);
          }

          // Increment the node number
          node_count++;
        }

        // Allocate memory for a new boundary node
        Node_pt[node_count] =
          finite_element_pt(Nx - 1)->construct_boundary_node(l1 * Np + (Np - 1),
                                                             time_stepper_pt);

        // If required make it periodic from the corresponding node on the other
        // side
        if (Xperiodic == true)
        {
          Node_pt[node_count]->make_periodic(
            finite_element_pt(0)->node_pt(l1 * Np));
        }

        // Set the position of the node
        Node_pt[node_count]->x(0) = x_spacing_function(Nx - 1, Np - 1, 0, l1);
        Node_pt[node_count]->x(1) = y_spacing_function(Nx - 1, Np - 1, 0, l1);

        // Push the node back onto boundaries
        add_boundary_node(1, Node_pt[node_count]);

        // If we only have one row, then the top node is on the top boundary
        if ((Ny == 1) && (l1 == (Np - 1)))
        {
          add_boundary_node(2, Node_pt[node_count]);
        }

        // Increment the node number
        node_count++;
      }
    }
    // END OF FIRST ROW OF ELEMENTS

    // ALL CENTRAL ELEMENT ROWS
    // Loop over remaining element rows
    for (unsigned i = 1; i < (Ny - 1); i++)
    {
      // Set the first element in the row
      // Allocate memory for element
      Element_pt[Nx * i] = new ELEMENT;

      // The first row of nodes is copied from the element below
      for (unsigned l2 = 0; l2 < Np; l2++)
      {
        finite_element_pt(Nx * i)->node_pt(l2) =
          finite_element_pt(Nx * (i - 1))->node_pt((Np - 1) * Np + l2);
      }

      // Other rows are new nodes
      for (unsigned l1 = 1; l1 < Np; l1++)
      {
        // First column of nodes
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(Nx * i)->construct_boundary_node(l1 * Np,
                                                             time_stepper_pt);

        // Set the position of the node
        Node_pt[node_count]->x(0) = x_spacing_function(0, 0, i, l1);
        Node_pt[node_count]->x(1) = y_spacing_function(0, 0, i, l1);

        // Push the node back onto boundaries
        add_boundary_node(3, Node_pt[node_count]);

        // Increment the node number
        node_count++;

        // Now do the other columns
        for (unsigned l2 = 1; l2 < Np; l2++)
        {
          // Allocate memory for node
          // If we only have one column, the node could be on the boundary
          if ((Nx == 1) && (l2 == (Np - 1)))
          {
            Node_pt[node_count] =
              finite_element_pt(Nx * i)->construct_boundary_node(
                l1 * Np + l2, time_stepper_pt);
          }
          else
          {
            Node_pt[node_count] = finite_element_pt(Nx * i)->construct_node(
              l1 * Np + l2, time_stepper_pt);
          }

          // Set the position of the node
          Node_pt[node_count]->x(0) = x_spacing_function(0, l2, i, l1);
          Node_pt[node_count]->x(1) = y_spacing_function(0, l2, i, l1);

          // If we only have one column then the RHS node is on the
          // right boundary
          if ((Nx == 1) && (l2 == (Np - 1)))
          {
            add_boundary_node(1, Node_pt[node_count]);
          }

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
        for (unsigned l2 = 0; l2 < Np; l2++)
        {
          finite_element_pt(Nx * i + j)->node_pt(l2) =
            finite_element_pt(Nx * (i - 1) + j)->node_pt((Np - 1) * Np + l2);
        }

        for (unsigned l1 = 1; l1 < Np; l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(Nx * i + j)->node_pt(l1 * Np) =
            finite_element_pt(Nx * i + (j - 1))->node_pt(l1 * Np + (Np - 1));
          // New nodes for other columns
          for (unsigned l2 = 1; l2 < Np; l2++)
          {
            // Allocate memory for the nodes
            Node_pt[node_count] =
              finite_element_pt(Nx * i + j)
                ->construct_node(l1 * Np + l2, time_stepper_pt);

            // Set the position of the node
            Node_pt[node_count]->x(0) = x_spacing_function(j, l2, i, l1);
            Node_pt[node_count]->x(1) = y_spacing_function(j, l2, i, l1);

            // Increment the node number
            node_count++;
          }
        }
      } // End of loop over elements in row

      // Do final element in row
      // Only if there is more than one column
      if (Nx > 1)
      {
        // Allocate memory for element
        Element_pt[Nx * i + Nx - 1] = new ELEMENT;
        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < Np; l2++)
        {
          finite_element_pt(Nx * i + Nx - 1)->node_pt(l2) =
            finite_element_pt(Nx * (i - 1) + Nx - 1)
              ->node_pt((Np - 1) * Np + l2);
        }

        for (unsigned l1 = 1; l1 < Np; l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(Nx * i + Nx - 1)->node_pt(l1 * Np) =
            finite_element_pt(Nx * i + Nx - 2)->node_pt(l1 * Np + (Np - 1));

          // Middle nodes
          for (unsigned l2 = 1; l2 < (Np - 1); l2++)
          {
            // Allocate memory for node
            Node_pt[node_count] =
              finite_element_pt(Nx * i + Nx - 1)
                ->construct_node(l1 * Np + l2, time_stepper_pt);

            // Set the position of the node
            Node_pt[node_count]->x(0) = x_spacing_function(Nx - 1, l2, i, l1);
            Node_pt[node_count]->x(1) = y_spacing_function(Nx - 1, l2, i, l1);

            // Increment the node number
            node_count++;
          }

          // Allocate memory for a new boundary node
          Node_pt[node_count] =
            finite_element_pt(Nx * i + Nx - 1)
              ->construct_boundary_node(l1 * Np + (Np - 1), time_stepper_pt);

          // If required make it periodic from the corresponding node on
          // the other side
          if (Xperiodic == true)
          {
            Node_pt[node_count]->make_periodic(
              finite_element_pt(Nx * i)->node_pt(l1 * Np));
          }

          // Set the position of the node
          Node_pt[node_count]->x(0) = x_spacing_function(Nx - 1, Np - 1, i, l1);
          Node_pt[node_count]->x(1) = y_spacing_function(Nx - 1, Np - 1, i, l1);

          // Push the node back onto boundaries
          add_boundary_node(1, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        } // End of loop over rows of nodes in the element
      } // End of if more than one column
    } // End of loop over rows of elements

    // END OF LOOP OVER CENTRAL ELEMENT ROWS

    // FINAL ELEMENT ROW
    // ONLY NECESSARY IF THERE IS MORE THAN ONE ROW
    if (Ny > 1)
    {
      // FIRST ELEMENT IN UPPER ROW (upper left corner)
      // Allocate memory for element
      Element_pt[Nx * (Ny - 1)] = new ELEMENT;
      // The first row of nodes is copied from the element below
      for (unsigned l2 = 0; l2 < Np; l2++)
      {
        finite_element_pt(Nx * (Ny - 1))->node_pt(l2) =
          finite_element_pt(Nx * (Ny - 2))->node_pt((Np - 1) * Np + l2);
      }

      // Second row of  nodes
      // First column of nodes
      for (unsigned l1 = 1; l1 < (Np - 1); l1++)
      {
        // Allocate memory for node
        Node_pt[node_count] =
          finite_element_pt(Nx * (Ny - 1))
            ->construct_boundary_node(Np * l1, time_stepper_pt);

        // Set the position of the node
        Node_pt[node_count]->x(0) = x_spacing_function(0, 0, Ny - 1, l1);
        Node_pt[node_count]->x(1) = y_spacing_function(0, 0, Ny - 1, l1);

        // Push the node back onto boundaries
        add_boundary_node(3, Node_pt[node_count]);

        // Increment the node number
        node_count++;

        // Now do the other columns
        for (unsigned l2 = 1; l2 < Np; l2++)
        {
          // Allocate memory for node
          if ((Nx == 1) && (l2 == Np - 1))
          {
            Node_pt[node_count] =
              finite_element_pt(Nx * (Ny - 1))
                ->construct_boundary_node(Np * l1 + l2, time_stepper_pt);
          }
          else
          {
            Node_pt[node_count] =
              finite_element_pt(Nx * (Ny - 1))
                ->construct_node(Np * l1 + l2, time_stepper_pt);
          }

          // Set the position of the node
          Node_pt[node_count]->x(0) = x_spacing_function(0, l2, Ny - 1, l1);
          Node_pt[node_count]->x(1) = y_spacing_function(0, l2, Ny - 1, l1);

          // Push the node back onto boundaries
          if ((Nx == 1) && (l2 == Np - 1))
          {
            add_boundary_node(1, Node_pt[node_count]);
          }

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
          ->construct_boundary_node(Np * (Np - 1), time_stepper_pt);

      // Set the position of the node
      Node_pt[node_count]->x(0) = x_spacing_function(0, 0, Ny - 1, Np - 1);
      Node_pt[node_count]->x(1) = y_spacing_function(0, 0, Ny - 1, Np - 1);

      // Push the node back onto boundaries
      add_boundary_node(2, Node_pt[node_count]);
      add_boundary_node(3, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Now do the other columns
      for (unsigned l2 = 1; l2 < Np; l2++)
      {
        // Allocate memory for the node
        Node_pt[node_count] =
          finite_element_pt(Nx * (Ny - 1))
            ->construct_boundary_node(Np * (Np - 1) + l2, time_stepper_pt);

        // Set the position of the node
        Node_pt[node_count]->x(0) = x_spacing_function(0, l2, Ny - 1, Np - 1);
        Node_pt[node_count]->x(1) = y_spacing_function(0, l2, Ny - 1, Np - 1);

        // Push the node back onto boundaries
        add_boundary_node(2, Node_pt[node_count]);

        // Push the node back onto boundaries
        if ((Nx == 1) && (l2 == Np - 1))
        {
          add_boundary_node(1, Node_pt[node_count]);
        }

        // Increment the node number
        node_count++;
      }

      // Now loop over the rest of the elements in the row, apart from the last
      for (unsigned j = 1; j < (Nx - 1); j++)
      {
        // Allocate memory for element
        Element_pt[Nx * (Ny - 1) + j] = new ELEMENT;
        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < Np; l2++)
        {
          finite_element_pt(Nx * (Ny - 1) + j)->node_pt(l2) =
            finite_element_pt(Nx * (Ny - 2) + j)->node_pt((Np - 1) * Np + l2);
        }

        // Second rows
        for (unsigned l1 = 1; l1 < (Np - 1); l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(Nx * (Ny - 1) + j)->node_pt(Np * l1) =
            finite_element_pt(Nx * (Ny - 1) + (j - 1))
              ->node_pt(Np * l1 + (Np - 1));

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < Np; l2++)
          {
            // Allocate memory for the node
            Node_pt[node_count] =
              finite_element_pt(Nx * (Ny - 1) + j)
                ->construct_node(Np * l1 + l2, time_stepper_pt);

            // Set the position of the node
            Node_pt[node_count]->x(0) = x_spacing_function(j, l2, Ny - 1, l1);
            Node_pt[node_count]->x(1) = y_spacing_function(j, l2, Ny - 1, l1);

            // Increment the node number
            node_count++;
          }
        }

        // Top row
        // First column is same as neighbouring element
        finite_element_pt(Nx * (Ny - 1) + j)->node_pt(Np * (Np - 1)) =
          finite_element_pt(Nx * (Ny - 1) + (j - 1))
            ->node_pt(Np * (Np - 1) + (Np - 1));
        // New nodes for other columns
        for (unsigned l2 = 1; l2 < Np; l2++)
        {
          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(Nx * (Ny - 1) + j)
              ->construct_boundary_node(Np * (Np - 1) + l2, time_stepper_pt);

          // Set the position of the node
          Node_pt[node_count]->x(0) = x_spacing_function(j, l2, Ny - 1, Np - 1);
          Node_pt[node_count]->x(1) = y_spacing_function(j, l2, Ny - 1, Np - 1);

          // Push the node back onto boundaries
          add_boundary_node(2, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }
      } // End of loop over central elements in row

      // FINAL ELEMENT IN ROW (upper right corner)
      // Only if there is more than one column
      if (Nx > 1)
      {
        // Allocate memory for element
        Element_pt[Nx * (Ny - 1) + Nx - 1] = new ELEMENT;
        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < Np; l2++)
        {
          finite_element_pt(Nx * (Ny - 1) + Nx - 1)->node_pt(l2) =
            finite_element_pt(Nx * (Ny - 2) + Nx - 1)
              ->node_pt((Np - 1) * Np + l2);
        }

        // Second rows
        for (unsigned l1 = 1; l1 < (Np - 1); l1++)
        {
          // First column is same as neighbouring element
          finite_element_pt(Nx * (Ny - 1) + Nx - 1)->node_pt(Np * l1) =
            finite_element_pt(Nx * (Ny - 1) + Nx - 2)
              ->node_pt(Np * l1 + (Np - 1));

          // Middle nodes
          for (unsigned l2 = 1; l2 < (Np - 1); l2++)
          {
            // Allocate memory for node
            Node_pt[node_count] =
              finite_element_pt(Nx * (Ny - 1) + Nx - 1)
                ->construct_node(Np * l1 + l2, time_stepper_pt);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
              x_spacing_function(Nx - 1, l2, Ny - 1, l1);
            Node_pt[node_count]->x(1) =
              y_spacing_function(Nx - 1, l2, Ny - 1, l1);

            // Increment the node number
            node_count++;
          }

          // Final node
          // Allocate new memory for a boundary node
          Node_pt[node_count] =
            finite_element_pt(Nx * (Ny - 1) + Nx - 1)
              ->construct_boundary_node(Np * l1 + (Np - 1), time_stepper_pt);

          // If required make it periodic
          if (Xperiodic == true)
          {
            Node_pt[node_count]->make_periodic(
              finite_element_pt(Nx * (Ny - 1))->node_pt(Np * l1));
          }

          // Set the position of the node
          Node_pt[node_count]->x(0) =
            x_spacing_function(Nx - 1, Np - 1, Ny - 1, l1);
          Node_pt[node_count]->x(1) =
            y_spacing_function(Nx - 1, Np - 1, Ny - 1, l1);

          // Push the node back onto boundaries
          add_boundary_node(1, Node_pt[node_count]);

          // Increment the node number
          node_count++;

        } // End of loop over middle rows

        // Final row
        // First column is same as neighbouring element
        finite_element_pt(Nx * (Ny - 1) + Nx - 1)->node_pt(Np * (Np - 1)) =
          finite_element_pt(Nx * (Ny - 1) + Nx - 2)
            ->node_pt(Np * (Np - 1) + (Np - 1));

        // Middle nodes
        for (unsigned l2 = 1; l2 < (Np - 1); l2++)
        {
          // Allocate memory for node
          Node_pt[node_count] =
            finite_element_pt(Nx * (Ny - 1) + Nx - 1)
              ->construct_boundary_node(Np * (Np - 1) + l2, time_stepper_pt);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
            x_spacing_function(Nx - 1, l2, Ny - 1, Np - 1);
          Node_pt[node_count]->x(1) =
            y_spacing_function(Nx - 1, l2, Ny - 1, Np - 1);

          // Push the node back onto boundaries
          add_boundary_node(2, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }

        // Final node
        // Allocate new memory for a periodic node
        Node_pt[node_count] = finite_element_pt(Nx * (Ny - 1) + Nx - 1)
                                ->construct_boundary_node(
                                  Np * (Np - 1) + (Np - 1), time_stepper_pt);

        // If required make it periodic
        if (Xperiodic == true)
        {
          Node_pt[node_count]->make_periodic(
            finite_element_pt(Nx * (Ny - 1))->node_pt(Np * (Np - 1)));
        }

        // Set the position of the node
        Node_pt[node_count]->x(0) =
          x_spacing_function(Nx - 1, Np - 1, Ny - 1, Np - 1);
        Node_pt[node_count]->x(1) =
          y_spacing_function(Nx - 1, Np - 1, Ny - 1, Np - 1);

        // Push the node back onto boundaries
        add_boundary_node(1, Node_pt[node_count]);
        add_boundary_node(2, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
      // END OF FINAL ELEMENT IN MESH
    }

    // Setup boundary element lookup schemes
    setup_boundary_element_info();
  }

  //============================================================================
  /// Reorder the elements so they are listed in vertical slices
  /// (more efficient during the frontal solution if the domain
  /// is long in the x-direction.
  //============================================================================
  template<class ELEMENT>
  void RectangularQuadMesh<ELEMENT>::element_reorder()
  {
    // Find out how many elements there are
    unsigned long Nelement = nelement();
    // Create a dummy array of elements
    Vector<FiniteElement*> dummy;

    // Loop over the elements in horizontal order
    for (unsigned long j = 0; j < Nx; j++)
    {
      // Loop over the elements vertically
      for (unsigned long i = 0; i < Ny; i++)
      {
        // Push back onto the new stack
        dummy.push_back(finite_element_pt(Nx * i + j));
      }
    }

    // Now copy the reordered elements into the element_pt
    for (unsigned long e = 0; e < Nelement; e++)
    {
      Element_pt[e] = dummy[e];
    }
  }

} // namespace oomph
#endif
