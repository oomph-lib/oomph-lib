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
#ifndef OOMPH_SIMPLE_CUBIC_TET_MESH_TEMPLATE_CC
#define OOMPH_SIMPLE_CUBIC_TET_MESH_TEMPLATE_CC

#include <algorithm>

// Simple 3D tetrahedral mesh class
#include "simple_cubic_tet_mesh.template.h"
#include "../generic/map_matrix.h"
#include <algorithm>


namespace oomph
{
  //====================================================================
  /// Simple tetrahedral mesh - with 24 tet elements constructed within a
  /// "brick" form for each element block.
  //====================================================================
  template<class ELEMENT>
  void SimpleCubicTetMesh<ELEMENT>::build_from_scaffold(
    TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with 3D Telements.
    MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(3);

    // Create space for elements
    unsigned nelem = Tmp_mesh_pt->nelement();
    Element_pt.resize(nelem);

    // Create space for nodes
    unsigned nnode_scaffold = Tmp_mesh_pt->nnode();
    Node_pt.resize(nnode_scaffold);

    // Set number of boundaries
    unsigned nbound = Tmp_mesh_pt->nboundary();
    set_nboundary(nbound);

    // Loop over elements in scaffold mesh, visit their nodes
    for (unsigned e = 0; e < nelem; e++)
    {
      Element_pt[e] = new ELEMENT;
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


    // Setup map to check the (pseudo-)global node number
    // Nodes whose number is zero haven't been copied across
    // into the mesh yet.
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
          Node_pt[global_count - 1]->x(2) = scaffold_node_pt->x(2);

          // Get pointer to set of mesh boundaries that this
          // scaffold node occupies; NULL if the node is not on any boundary
          std::set<unsigned>* boundaries_pt;
          scaffold_node_pt->get_boundaries_pt(boundaries_pt);

          // Loop over the mesh boundaries that the node in the scaffold mesh
          // occupies and assign new node to the same ones.
          if (boundaries_pt != 0)
          {
            this->convert_to_boundary_node(Node_pt[global_count - 1]);
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


    // At this point we've created all the elements and
    // created their vertex nodes. Now we need to create
    // the additional (midside and internal) nodes!


    // We'll first create all local nodes for all elements
    // and then delete the superfluous ones that have
    // a matching node in an adjacent element.

    // Get number of nodes along element edge and dimension of element (3)
    // from first element
    unsigned nnode_1d = finite_element_pt(0)->nnode_1d();

    // At the moment we're only able to deal with nnode_1d=2 or 3.
    if ((nnode_1d != 2) && (nnode_1d != 3))
    {
      std::ostringstream error_message;
      error_message << "Mesh generation currently only works\n";
      error_message << "for nnode_1d = 2 or 3. You're trying to use it\n";
      error_message << "for nnode_1d=" << nnode_1d << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Spatial dimension of element = number of local coordinates
    unsigned dim = finite_element_pt(0)->dim();

    // Storage for the local coordinate of the new node
    Vector<double> s(dim);

    // Get number of nodes in the element from first element
    unsigned nnode = finite_element_pt(0)->nnode();

    // Loop over all elements
    for (unsigned e = 0; e < nelem; e++)
    {
      // Loop over the new nodes in the element and create them.
      for (unsigned j = 4; j < nnode; j++)
      {
        // Create new node
        Node* new_node_pt =
          finite_element_pt(e)->construct_node(j, time_stepper_pt);

        // What are the node's local coordinates?
        finite_element_pt(e)->local_coordinate_of_node(j, s);

        // Find the coordinates of the new node from the existing
        // and fully-functional element in the scaffold mesh
        for (unsigned i = 0; i < dim; i++)
        {
          new_node_pt->x(i) =
            Tmp_mesh_pt->finite_element_pt(e)->interpolated_x(s, i);
        }
      } // end of loop over new nodes
    } // end of loop over elements


    // Bracket this away so the edge map goes out of scope
    // when we're done
    {
      // Storage for pointer to mid-edge node
      MapMatrix<Node*, Node*> central_edge_node_pt;
      Node* edge_node1_pt = 0;
      Node* edge_node2_pt = 0;

      // Loop over elements
      for (unsigned e = 0; e < nelem; e++)
      {
        // Loop over new local nodes
        for (unsigned j = 4; j < nnode; j++)
        {
          // Pointer to the element's local node
          Node* node_pt = finite_element_pt(e)->node_pt(j);

          // By default, we assume the node is not new
          bool is_new = false;

          // This will have to be changed for higher-order elements
          //=======================================================

          // Switch on local node number (all located on edges)
          switch (j)
          {
              // Node 4 is located between nodes 0 and 1
            case 4:

              edge_node1_pt = finite_element_pt(e)->node_pt(0);
              edge_node2_pt = finite_element_pt(e)->node_pt(1);
              if (central_edge_node_pt(edge_node1_pt, edge_node2_pt) == 0)
              {
                is_new = true;
                central_edge_node_pt(edge_node1_pt, edge_node2_pt) = node_pt;
                central_edge_node_pt(edge_node2_pt, edge_node1_pt) = node_pt;
              }
              break;


              // Node 5 is located between nodes 0 and 2
            case 5:

              edge_node1_pt = finite_element_pt(e)->node_pt(0);
              edge_node2_pt = finite_element_pt(e)->node_pt(2);
              if (central_edge_node_pt(edge_node1_pt, edge_node2_pt) == 0)
              {
                is_new = true;
                central_edge_node_pt(edge_node1_pt, edge_node2_pt) = node_pt;
                central_edge_node_pt(edge_node2_pt, edge_node1_pt) = node_pt;
              }
              break;


              // Node 6 is located between nodes 0 and 3
            case 6:

              edge_node1_pt = finite_element_pt(e)->node_pt(0);
              edge_node2_pt = finite_element_pt(e)->node_pt(3);
              if (central_edge_node_pt(edge_node1_pt, edge_node2_pt) == 0)
              {
                is_new = true;
                central_edge_node_pt(edge_node1_pt, edge_node2_pt) = node_pt;
                central_edge_node_pt(edge_node2_pt, edge_node1_pt) = node_pt;
              }
              break;


              // Node 7 is located between nodes 1 and 2
            case 7:

              edge_node1_pt = finite_element_pt(e)->node_pt(1);
              edge_node2_pt = finite_element_pt(e)->node_pt(2);
              if (central_edge_node_pt(edge_node1_pt, edge_node2_pt) == 0)
              {
                is_new = true;
                central_edge_node_pt(edge_node1_pt, edge_node2_pt) = node_pt;
                central_edge_node_pt(edge_node2_pt, edge_node1_pt) = node_pt;
              }
              break;


              // Node 8 is located between nodes 2 and 3
            case 8:

              edge_node1_pt = finite_element_pt(e)->node_pt(2);
              edge_node2_pt = finite_element_pt(e)->node_pt(3);
              if (central_edge_node_pt(edge_node1_pt, edge_node2_pt) == 0)
              {
                is_new = true;
                central_edge_node_pt(edge_node1_pt, edge_node2_pt) = node_pt;
                central_edge_node_pt(edge_node2_pt, edge_node1_pt) = node_pt;
              }
              break;


              // Node 9 is located between nodes 1 and 3
            case 9:

              edge_node1_pt = finite_element_pt(e)->node_pt(1);
              edge_node2_pt = finite_element_pt(e)->node_pt(3);
              if (central_edge_node_pt(edge_node1_pt, edge_node2_pt) == 0)
              {
                is_new = true;
                central_edge_node_pt(edge_node1_pt, edge_node2_pt) = node_pt;
                central_edge_node_pt(edge_node2_pt, edge_node1_pt) = node_pt;
              }
              break;

            default:
              throw OomphLibError("More than ten nodes in Tet Element",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
          }

          if (is_new)
          {
            // New node: Add it to mesh
            Node_pt.push_back(node_pt);
          }
          else
          {
            // Delete local node in element...
            delete finite_element_pt(e)->node_pt(j);

            // ... and reset pointer to the existing node
            finite_element_pt(e)->node_pt(j) =
              central_edge_node_pt(edge_node1_pt, edge_node2_pt);
          }
        }
      }
    }


    // Boundary conditions

    // Matrix map to check if a node has already been added to
    // the boundary number b (NOTE: Enumerated by pointer to ORIGINAL
    // node before transfer to boundary node)
    MapMatrixMixed<Node*, unsigned, bool> bound_node_done;

    // Create lookup scheme for pointers to original local nodes
    // in elements
    Vector<Vector<Node*>> orig_node_pt(nelem);

    // Loop over elements
    for (unsigned e = 0; e < nelem; e++)
    {
      orig_node_pt[e].resize(nnode, 0);

      // Loop over new local nodes
      for (unsigned j = 4; j < nnode; j++)
      {
        // Pointer to the element's local node
        orig_node_pt[e][j] = finite_element_pt(e)->node_pt(j);
      }
    }


    // Loop over elements
    for (unsigned e = 0; e < nelem; e++)
    {
      // Loop over new local nodes
      for (unsigned j = 4; j < nnode; j++)
      {
        // Loop over the boundaries
        for (unsigned bo = 0; bo < nbound; bo++)
        {
          // Pointer to the element's local node
          Node* loc_node_pt = finite_element_pt(e)->node_pt(j);

          // Pointer to original node
          Node* orig_loc_node_pt = orig_node_pt[e][j];

          // value of the map for the node and boundary specified
          bool bound_test = bound_node_done(orig_loc_node_pt, bo);

          if (bound_test == false)
          {
            bound_node_done(orig_loc_node_pt, bo) = true;

            // This will have to be changed for higher-order elements
            //=======================================================

            // Switch on local node number (all located on edges)
            switch (j)
            {
                // Node 4 is located between nodes 0 and 1
              case 4:

                if (finite_element_pt(e)->node_pt(0)->is_on_boundary(bo) &&
                    finite_element_pt(e)->node_pt(1)->is_on_boundary(bo))
                {
                  this->convert_to_boundary_node(loc_node_pt);
                  add_boundary_node(bo, loc_node_pt);
                }
                break;


                // Node 5 is located between nodes 0 and 2
              case 5:

                if (finite_element_pt(e)->node_pt(0)->is_on_boundary(bo) &&
                    finite_element_pt(e)->node_pt(2)->is_on_boundary(bo))
                {
                  this->convert_to_boundary_node(loc_node_pt);
                  add_boundary_node(bo, loc_node_pt);
                }
                break;


                // Node 6 is located between nodes 0 and 3
              case 6:

                if (finite_element_pt(e)->node_pt(0)->is_on_boundary(bo) &&
                    finite_element_pt(e)->node_pt(3)->is_on_boundary(bo))
                {
                  this->convert_to_boundary_node(loc_node_pt);
                  add_boundary_node(bo, loc_node_pt);
                }
                break;


                // Node 7 is located between nodes 1 and 2
              case 7:

                if (finite_element_pt(e)->node_pt(1)->is_on_boundary(bo) &&
                    finite_element_pt(e)->node_pt(2)->is_on_boundary(bo))
                {
                  this->convert_to_boundary_node(loc_node_pt);
                  add_boundary_node(bo, loc_node_pt);
                }
                break;


                // Node 8 is located between nodes 2 and 3
              case 8:

                if (finite_element_pt(e)->node_pt(2)->is_on_boundary(bo) &&
                    finite_element_pt(e)->node_pt(3)->is_on_boundary(bo))
                {
                  this->convert_to_boundary_node(loc_node_pt);
                  add_boundary_node(bo, loc_node_pt);
                }
                break;


                // Node 9 is located between nodes 1 and 3
              case 9:

                if (finite_element_pt(e)->node_pt(1)->is_on_boundary(bo) &&
                    finite_element_pt(e)->node_pt(3)->is_on_boundary(bo))
                {
                  this->convert_to_boundary_node(loc_node_pt);
                  add_boundary_node(bo, loc_node_pt);
                }
                break;


              default:
                throw OomphLibError("More than ten nodes in Tet Element",
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
            }
          }

        } // end for bo
      } // end for j
    } // end for e
  }

} // namespace oomph

#endif
