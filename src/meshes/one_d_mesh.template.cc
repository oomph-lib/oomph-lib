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
#ifndef OOMPH_ONE_D_MESH_TEMPLATE_CC
#define OOMPH_ONE_D_MESH_TEMPLATE_CC


// oomph-lib headers
#include "one_d_mesh.template.h"


namespace oomph
{
  /// The generic mesh construction routine --- this contains all the hard
  /// work and is called by all constructors
  template<class ELEMENT>
  void OneDMesh<ELEMENT>::build_mesh(TimeStepper* time_stepper_pt)
  {
    // Set the length of the domain
    Length = Xmax - Xmin;

    // Set the number of boundaries -- there are 2 boundaries in a 1D mesh
    set_nboundary(2);

    // Allocate storage for the pointers to the elements
    Element_pt.resize(N);

    // Allocate memory for the first element
    Element_pt[0] = new ELEMENT;

    // Read out the number of nodes in the element (the member function
    // nnode_1d() is implemented in QElement)
    const unsigned n_node =
      dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

    // We can now allocate storage for the pointers to the nodes in the mesh
    Node_pt.resize(1 + (n_node - 1) * N);

    // Initialise the node counter
    unsigned node_count = 0;

    // Initialise the minimum x coordinate in the mesh
    const double xinit = Xmin;

    // Calculate the length of the element
    const double el_length = Length / double(N);

    // Allocate storage for the local coordinate in the element
    Vector<double> s_fraction;

    // If the number of elements is 1, the first element is also the
    // last element
    if (N == 1)
    {
      // Set the first node
      // ------------------

      // Allocate memory for the node, using the element's own construct_node
      // function -- only the element knows what type of nodes it needs!
      Node_pt[node_count] =
        finite_element_pt(0)->construct_boundary_node(0, time_stepper_pt);

      // Set the position of the node
      node_pt(node_count)->x(0) = xinit;

      // Add the node to the boundary 0
      add_boundary_node(0, Node_pt[node_count]);

      // Increment the counter for the nodes
      node_count++;

      // Now build central nodes (ignore last one which needs special
      // ------------------------------------------------------------
      // treatment because it's on the boundary)
      // ---------------------------------------
      for (unsigned jnod = 1; jnod < (n_node - 1); jnod++)
      {
        // Allocate memory for nodes, as before
        Node_pt[node_count] =
          finite_element_pt(0)->construct_node(jnod, time_stepper_pt);

        // Get the local coordinate of the node
        finite_element_pt(0)->local_fraction_of_node(jnod, s_fraction);

        // Set the position of the node (linear mapping)
        node_pt(node_count)->x(0) = xinit + el_length * s_fraction[0];

        // Increment the node counter
        node_count++;
      }

      // New final node
      // --------------

      // Allocate memory for the node, as before
      Node_pt[node_count] = finite_element_pt(0)->construct_boundary_node(
        n_node - 1, time_stepper_pt);

      // Set the position of the node
      node_pt(node_count)->x(0) = xinit + Length;

      // Add the node to the boundary 1
      add_boundary_node(1, Node_pt[node_count]);

      // Increment the node counter
      node_count++;
    }

    // Otherwise, i.e. if there is more than one element, build all elements
    else
    {
      // -------------
      // FIRST ELEMENT
      // -------------

      // Set the first node
      // ------------------

      // Allocate memory for the node, using the element's own construct_node
      // function -- only the element knows what type of nodes it needs!
      Node_pt[node_count] =
        finite_element_pt(0)->construct_boundary_node(0, time_stepper_pt);

      // Set the position of the node
      node_pt(node_count)->x(0) = xinit;

      // Add the node to the boundary 0
      add_boundary_node(0, Node_pt[node_count]);

      // Increment the counter for the nodes
      node_count++;

      // Now build the other nodes in the first element
      // ----------------------------------------------

      // Loop over the other nodes in the first element
      for (unsigned jnod = 1; jnod < n_node; jnod++)
      {
        // Allocate memory for the nodes
        Node_pt[node_count] =
          finite_element_pt(0)->construct_node(jnod, time_stepper_pt);

        // Get the local coordinate of the node
        finite_element_pt(0)->local_fraction_of_node(jnod, s_fraction);

        // Set the position of the node (linear mapping)
        node_pt(node_count)->x(0) = xinit + el_length * s_fraction[0];

        // Increment the node counter
        node_count++;
      }

      // ----------------
      // CENTRAL ELEMENTS
      // ----------------

      // Loop over central elements in mesh
      for (unsigned e = 1; e < (N - 1); e++)
      {
        // Allocate memory for the new element
        Element_pt[e] = new ELEMENT;

        // The first node is the same as the last node in the neighbouring
        // element on the left
        finite_element_pt(e)->node_pt(0) =
          finite_element_pt(e - 1)->node_pt((n_node - 1));

        // Loop over the remaining nodes in the element
        for (unsigned jnod = 1; jnod < n_node; jnod++)
        {
          // Allocate memory for the nodes, as before
          Node_pt[node_count] =
            finite_element_pt(e)->construct_node(jnod, time_stepper_pt);

          // Get the local coordinate of the nodes
          finite_element_pt(e)->local_fraction_of_node(jnod, s_fraction);

          // Set the position of the node (linear mapping)
          node_pt(node_count)->x(0) = xinit + el_length * (e + s_fraction[0]);

          // Increment the node counter
          node_count++;
        }
      } // End of loop over central elements


      // FINAL ELEMENT
      //--------------

      // Allocate memory for element
      Element_pt[N - 1] = new ELEMENT;

      // The first node is the same as the last node in the neighbouring
      // element on the left
      finite_element_pt(N - 1)->node_pt(0) =
        finite_element_pt(N - 2)->node_pt(n_node - 1);

      // New central nodes (ignore last one which needs special treatment
      // because it's on the boundary)
      for (unsigned jnod = 1; jnod < (n_node - 1); jnod++)
      {
        // Allocate memory for nodes, as before
        Node_pt[node_count] =
          finite_element_pt(N - 1)->construct_node(jnod, time_stepper_pt);

        // Get the local coordinate of the node
        finite_element_pt(N - 1)->local_fraction_of_node(jnod, s_fraction);

        // Set the position of the node
        node_pt(node_count)->x(0) = xinit + el_length * (N - 1 + s_fraction[0]);

        // Increment the node counter
        node_count++;
      }

      // New final node
      // --------------

      // Allocate memory for the node, as before
      Node_pt[node_count] = finite_element_pt(N - 1)->construct_boundary_node(
        n_node - 1, time_stepper_pt);

      // Set the position of the node
      node_pt(node_count)->x(0) = xinit + Length;

      // Add the node to the boundary 1
      add_boundary_node(1, Node_pt[node_count]);

      // Increment the node counter
      node_count++;
    }
  }

} // namespace oomph

#endif
