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
#ifndef OOMPH_GEOMPACK_MESH_TEMPLATE_CC
#define OOMPH_GEOMPACK_MESH_TEMPLATE_CC

#include "geompack_mesh.template.h"


namespace oomph
{
  //========================================================================
  /// Quadrilateral mesh generator; Uses input from Geompack++.
  /// See: http://members.shaw.ca/bjoe/
  /// Currently only for four-noded quads -- extension to higher-order
  /// quads should be trivial (see the corresponding classes for
  /// triangular meshes).
  //========================================================================
  template<class ELEMENT>
  void GeompackQuadMesh<ELEMENT>::build_from_scaffold(
    TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with four-noded 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2, 2);

    // Create space for elements
    unsigned nelem = Tmp_mesh_pt->nelement();
    Element_pt.resize(nelem);

    // Create space for nodes
    unsigned nnod = Tmp_mesh_pt->nnode();
    Node_pt.resize(nnod);

    // Set number of boundaries
    unsigned nbound = Tmp_mesh_pt->nboundary();
    set_nboundary(nbound);

    // Loop over elements in scaffold mesh, visit their nodes
    for (unsigned e = 0; e < nelem; e++)
    {
      Element_pt[e] = new ELEMENT;
    }


    // At the moment we can only do four-noded quads
    if (finite_element_pt(0)->nnode() != 4)
    {
      std::string error_message =
        "GeompackQuadMesh can currently only deal with\n";
      error_message += "four noded quads! \n";
      error_message += "Generalisation to higher-order elements should be \n";
      error_message += "straightforward but hasn't been done yet.\n";
      error_message += "Do you want to volunteer? \n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    // In the first instance build all nodes from within all the elements
    unsigned nnod_el = Tmp_mesh_pt->finite_element_pt(0)->nnode();

    // Loop over elements in scaffold mesh, visit their nodes
    for (unsigned e = 0; e < nelem; e++)
    {
      // Loop over all nodes in element
      for (unsigned j = 0; j < nnod_el; j++)
      {
        // Create new node, using the NEW element's construct_node
        // member function
        finite_element_pt(e)->construct_node(j, time_stepper_pt);
      }
    }


    std::map<Node*, unsigned> global_number;
    unsigned global_count = 0;
    // Loop over elements in scaffold mesh, visit their nodes
    for (unsigned e = 0; e < nelem; e++)
    {
      // Loop over all nodes in element
      for (unsigned j = 0; j < nnod_el; j++)
      {
        // Pointer to node in the scaffold mesh
        Node* scaffold_node_pt = Tmp_mesh_pt->finite_element_pt(e)->node_pt(j);

        // Get the (pseudo-)global node number in scaffold mesh
        // (It's zero [=default] if not visited this one yet)
        unsigned j_global = global_number[scaffold_node_pt];

        // Haven't done this one yet
        if (j_global == 0)
        {
          // Give it a number (not necessarily the global node
          // number in the scaffold mesh -- we just need something
          // to keep track...)
          global_count++;
          global_number[scaffold_node_pt] = global_count;

          // Copy new node, created using the NEW element's construct_node
          // function into global storage, using the same global
          // node number that we've just associated with the
          // corresponding node in the scaffold mesh
          Node_pt[global_count - 1] = finite_element_pt(e)->node_pt(j);

          // Assign coordinates
          Node_pt[global_count - 1]->x(0) = scaffold_node_pt->x(0);
          Node_pt[global_count - 1]->x(1) = scaffold_node_pt->x(1);


          // Get pointer to set of mesh boundaries that this
          // scaffold node occupies; NULL if the node is not on any boundary
          std::set<unsigned>* boundaries_pt;
          scaffold_node_pt->get_boundaries_pt(boundaries_pt);

          // Loop over the mesh boundaries that the node in the scaffold mesh
          // occupies and assign new node to the same ones.
          if (boundaries_pt != 0)
          {
            // Convert it to a boundary node
            this->convert_to_boundary_node(Node_pt[global_count - 1]);
            // Add the node to the boundaries
            for (std::set<unsigned>::iterator it = boundaries_pt->begin();
                 it != boundaries_pt->end();
                 ++it)
            {
              add_boundary_node(*it, Node_pt[global_count - 1]);
            }
          }
        }
        // This one has already been done: Kill it
        else
        {
          // Kill it
          delete finite_element_pt(e)->node_pt(j);

          // Overwrite the element's pointer to local node by
          // pointer to the already existing node -- identified by
          // the number kept in the map
          finite_element_pt(e)->node_pt(j) = Node_pt[j_global - 1];
        }
      }
    }
  }

} // namespace oomph
#endif
