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
#ifndef OOMPH_TETGEN_MESH_TEMPLATE_CC
#define OOMPH_TETGEN_MESH_TEMPLATE_CC


#include <algorithm>

#include "tetgen_mesh.template.h"
#include "../generic/Telements.h"
#include "../generic/map_matrix.h"


namespace oomph
{
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Build unstructured tet mesh based on output from scaffold
  //========================================================================
  template<class ELEMENT>
  void TetgenMesh<ELEMENT>::build_from_scaffold(TimeStepper* time_stepper_pt,
                                                const bool& use_attributes)
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
    Boundary_element_pt.resize(nbound);
    Face_index_at_boundary.resize(nbound);

    // If we have different regions, then resize the region
    // information
    if (use_attributes)
    {
      Boundary_region_element_pt.resize(nbound);
      Face_index_region_at_boundary.resize(nbound);
    }

    // Build elements
    for (unsigned e = 0; e < nelem; e++)
    {
      Element_pt[e] = new ELEMENT;
    }

    // Number of nodes per element
    unsigned nnod_el = Tmp_mesh_pt->finite_element_pt(0)->nnode();

    // Setup map to check the (pseudo-)global node number
    // Nodes whose number is zero haven't been copied across
    // into the mesh yet.
    std::map<Node*, unsigned> global_number;
    unsigned global_count = 0;

    // Map of element attribute pairs
    std::map<double, Vector<FiniteElement*>> element_attribute_map;

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
          // Get pointer to set of mesh boundaries that this
          // scaffold node occupies; NULL if the node is not on any boundary
          std::set<unsigned>* boundaries_pt;
          scaffold_node_pt->get_boundaries_pt(boundaries_pt);

          // Is it on boundaries?
          if (boundaries_pt != 0)
          {
            // Create new boundary node
            Node* new_node_pt =
              finite_element_pt(e)->construct_boundary_node(j, time_stepper_pt);

            // Give it a number (not necessarily the global node
            // number in the scaffold mesh -- we just need something
            // to keep track...)
            global_count++;
            global_number[scaffold_node_pt] = global_count;

            // Add to boundaries
            for (std::set<unsigned>::iterator it = boundaries_pt->begin();
                 it != boundaries_pt->end();
                 ++it)
            {
              add_boundary_node(*it, new_node_pt);
            }
          }
          // Build normal node
          else
          {
            // Create new normal node
            finite_element_pt(e)->construct_node(j, time_stepper_pt);

            // Give it a number (not necessarily the global node
            // number in the scaffold mesh -- we just need something
            // to keep track...)
            global_count++;
            global_number[scaffold_node_pt] = global_count;
          }

          // Copy new node, created using the NEW element's construct_node
          // function into global storage, using the same global
          // node number that we've just associated with the
          // corresponding node in the scaffold mesh
          Node_pt[global_count - 1] = finite_element_pt(e)->node_pt(j);

          // Assign coordinates
          Node_pt[global_count - 1]->x(0) = scaffold_node_pt->x(0);
          Node_pt[global_count - 1]->x(1) = scaffold_node_pt->x(1);
          Node_pt[global_count - 1]->x(2) = scaffold_node_pt->x(2);
        }
        // This one has already been done: Copy across
        else
        {
          finite_element_pt(e)->node_pt(j) = Node_pt[j_global - 1];
        }
      }

      // Store the attributes in the map
      if (use_attributes)
      {
        element_attribute_map[Tmp_mesh_pt->element_attribute(e)].push_back(
          finite_element_pt(e));
      }
    }

    // Now let's construct lists
    // Find the number of attributes
    if (use_attributes)
    {
      unsigned n_attribute = element_attribute_map.size();
      // There are n_attribute different regions
      Region_element_pt.resize(n_attribute);
      Region_attribute.resize(n_attribute);
      // Copy the vectors in the map over to our internal storage
      unsigned count = 0;
      for (std::map<double, Vector<FiniteElement*>>::iterator it =
             element_attribute_map.begin();
           it != element_attribute_map.end();
           ++it)
      {
        Region_attribute[count] = it->first;
        Region_element_pt[count] = it->second;
        ++count;
      }
    }

    // At this point we've created all the elements and
    // created their vertex nodes. Now we need to create
    // the additional (midside and internal) nodes!
    unsigned boundary_id;

    // Get number of nodes along element edge and dimension of element (3)
    // from first element
    unsigned n_node_1d = finite_element_pt(0)->nnode_1d();

    // At the moment we're only able to deal with nnode_1d=2 or 3.
    if ((n_node_1d != 2) && (n_node_1d != 3))
    {
      std::ostringstream error_message;
      error_message << "Mesh generation from tetgen currently only works\n";
      error_message << "for nnode_1d = 2 or 3. You're trying to use it\n";
      error_message << "for nnode_1d=" << n_node_1d << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Spatial dimension of element = number of local coordinates
    unsigned dim = finite_element_pt(0)->dim();

    // Storage for the local coordinate of the new node
    Vector<double> s(dim);

    // Get number of nodes in the element from first element
    unsigned n_node = finite_element_pt(0)->nnode();

    // Storage for each global edge of the mesh
    unsigned n_global_edge = Tmp_mesh_pt->nglobal_edge();
    Vector<Vector<Node*>> nodes_on_global_edge(n_global_edge);

    // Storage for each global face of the mesh
    unsigned n_global_face = Tmp_mesh_pt->nglobal_face();
    Vector<Vector<Node*>> nodes_on_global_face(n_global_face);

    // Map storing the mid-side of an edge; edge identified by
    // pointers to vertex nodes in scaffold mesh
    // MapMatrix<Node*,Node*> central_edge_node_pt;
    // Node* edge_node1_pt=0;
    // Node* edge_node2_pt=0;

    // Map storing the mid point of a face; face identified by
    // set of pointers to vertex nodes in scaffold mesh
    // std::map<std::set<Node*>,Node*> central_face_node_pt;
    // std::set<Node*> face_nodes_pt;

    // Mapping of Tetgen faces to face nodes in the enriched element
    unsigned face_map[4] = {1, 0, 2, 3};

    // Storage for the faces shared by the edges
    const unsigned faces_on_edge[6][2] = {
      {0, 1}, {0, 2}, {1, 2}, {0, 3}, {2, 3}, {1, 3}};

    // Loop over all elements
    for (unsigned e = 0; e < nelem; e++)
    {
      // Cache pointers to the elements
      FiniteElement* const elem_pt = this->finite_element_pt(e);
      FiniteElement* const tmp_elem_pt = Tmp_mesh_pt->finite_element_pt(e);

      // The number of edge nodes is 4 + 6*(n_node1d-2)
      unsigned n_edge_node = 4 + 6 * (n_node_1d - 2);

      // Now loop over the edge nodes
      // assuming that the numbering scheme is the same as the triangles
      // which puts edge nodes in ascending order.
      // We don't have any higher than quadratic at the moment, so it's
      // a bit academic really.

      // Only bother if we actually have some edge nodes
      if (n_edge_node > 4)
      {
        // Start from node number 4
        unsigned n = 4;

        // Loop over the edges
        for (unsigned j = 0; j < 6; ++j)
        {
          // Find the global edge index
          unsigned edge_index = Tmp_mesh_pt->edge_index(e, j);

          // Use the intersection of the appropriate faces to determine
          // whether the boundaries on which an edge lies
          std::set<unsigned> edge_boundaries;
          for (unsigned i = 0; i < 2; ++i)
          {
            unsigned face_boundary_id =
              Tmp_mesh_pt->face_boundary(e, faces_on_edge[j][i]);
            if (face_boundary_id > 0)
            {
              edge_boundaries.insert(face_boundary_id);
            }
          }

          // If the nodes on the edge have not been allocated, construct them
          if (nodes_on_global_edge[edge_index].size() == 0)
          {
            // Now loop over the nodes on the edge
            for (unsigned j2 = 0; j2 < n_node_1d - 2; ++j2)
            {
              // Storage for the new node
              Node* new_node_pt = 0;

              // If the edge is on a boundary, determined from the
              // scaffold mesh, construct a boundary node
              // The use of the scaffold mesh's edge_boundary data structure
              // ensures that a boundary node is created, even if the faces of
              // the current element do not lie on boundaries
              if (Tmp_mesh_pt->edge_boundary(edge_index) == true)
              {
                new_node_pt =
                  elem_pt->construct_boundary_node(n, time_stepper_pt);
                // Add it to the boundaries in the set,
                // remembering to subtract one to get to the oomph-lib numbering
                // scheme
                for (std::set<unsigned>::iterator it = edge_boundaries.begin();
                     it != edge_boundaries.end();
                     ++it)
                {
                  this->add_boundary_node((*it) - 1, new_node_pt);
                }
              }
              // Otherwise construct a normal node
              else
              {
                new_node_pt = elem_pt->construct_node(n, time_stepper_pt);
              }

              // Find the local coordinates of the node
              elem_pt->local_coordinate_of_node(n, s);

              // Find the coordinates of the new node from the exiting and
              // fully-functional element in the scaffold mesh
              for (unsigned i = 0; i < dim; ++i)
              {
                new_node_pt->x(i) = tmp_elem_pt->interpolated_x(s, i);
              }

              // Add the newly created node to the global node list
              Node_pt.push_back(new_node_pt);
              // Add to the edge index
              nodes_on_global_edge[edge_index].push_back(new_node_pt);
              // Increment the local node number
              ++n;
            } // end of loop over edge nodes
          }
          // Otherwise just set the pointers (orientation assumed the same)
          else
          {
            for (unsigned j2 = 0; j2 < n_node_1d - 2; ++j2)
            {
              elem_pt->node_pt(n) = nodes_on_global_edge[edge_index][j2];
              // It is possible that the edge may be on additional boundaries
              // through another element
              // So add again (note that this function will not add to
              // boundaries twice)
              for (std::set<unsigned>::iterator it = edge_boundaries.begin();
                   it != edge_boundaries.end();
                   ++it)
              {
                this->add_boundary_node((*it) - 1, elem_pt->node_pt(n));
              }
              ++n;
            }
          }
        } // End of loop over edges

        // Deal with enriched elements
        if (n_node == 15)
        {
          // Now loop over the faces
          for (unsigned j = 0; j < 4; ++j)
          {
            // Find the boundary id of the face (need to map between node
            // numbers and the face)
            boundary_id = Tmp_mesh_pt->face_boundary(e, face_map[j]);

            // Find the global face index (need to map between node numbers and
            // the face)
            unsigned face_index = Tmp_mesh_pt->face_index(e, face_map[j]);

            // If the nodes on the face have not been allocated
            if (nodes_on_global_face[face_index].size() == 0)
            {
              // Storage for the new node
              Node* new_node_pt = 0;

              // If the face is on a boundary, construct a boundary node
              if (boundary_id > 0)
              {
                new_node_pt =
                  elem_pt->construct_boundary_node(n, time_stepper_pt);
                // Add it to the boundary
                this->add_boundary_node(boundary_id - 1, new_node_pt);
              }
              // Otherwise construct a normal node
              else
              {
                new_node_pt = elem_pt->construct_node(n, time_stepper_pt);
              }

              // Find the local coordinates of the node
              elem_pt->local_coordinate_of_node(n, s);

              // Find the coordinates of the new node from the exiting and
              // fully-functional element in the scaffold mesh
              for (unsigned i = 0; i < dim; ++i)
              {
                new_node_pt->x(i) = tmp_elem_pt->interpolated_x(s, i);
              }

              // Add the newly created node to the global node list
              Node_pt.push_back(new_node_pt);
              // Add to the face index
              nodes_on_global_face[face_index].push_back(new_node_pt);
              // Increment the local node number
              ++n;
            }
            // Otherwise just set the single node from the face element
            else
            {
              elem_pt->node_pt(n) = nodes_on_global_face[face_index][0];
              ++n;
            }
          } // end of loop over faces

          // Construct the element's central node, which is not on a boundary
          {
            Node* new_node_pt =
              finite_element_pt(e)->construct_node(n, time_stepper_pt);
            Node_pt.push_back(new_node_pt);

            // Find the local coordinates of the node
            elem_pt->local_coordinate_of_node(n, s);

            // Find the coordinates of the new node from the existing
            // and fully-functional element in the scaffold mesh
            for (unsigned i = 0; i < dim; i++)
            {
              new_node_pt->x(i) = tmp_elem_pt->interpolated_x(s, i);
            }
          }
        } // End of enriched case

      } // End of case for edge nodes


      // Now loop over the faces to setup the information about elements
      // adjacent to the boundary
      for (unsigned j = 0; j < 4; ++j)
      {
        // Find the boundary id of the face
        boundary_id = Tmp_mesh_pt->face_boundary(e, j);

        if (boundary_id > 0)
        {
          Boundary_element_pt[boundary_id - 1].push_back(elem_pt);
          // Need to put a shift in here because of an inconsistent naming
          // convention between tetgen and our faces
          // Tetgen Face 0 is our Face 3
          // Tetgen Face 1 is our Face 2
          // Tetgen Face 2 is our Face 1
          // Tetgen Face 3 is our Face 0
          Face_index_at_boundary[boundary_id - 1].push_back(3 - j);

          // If using regions set up the boundary information
          if (use_attributes)
          {
            // Element adjacent to boundary
            Boundary_region_element_pt[boundary_id - 1]
                                      [static_cast<unsigned>(
                                         Tmp_mesh_pt->element_attribute(e))]
                                        .push_back(elem_pt);
            // Need to put a shift in here because of an inconsistent naming
            // convention between triangle and face elements
            Face_index_region_at_boundary[boundary_id - 1]
                                         [static_cast<unsigned>(
                                            Tmp_mesh_pt->element_attribute(e))]
                                           .push_back(3 - j);
          }
        }
      } // End of loop over faces


      // Lookup scheme has now been setup
      Lookup_for_elements_next_boundary_is_setup = true;


      //   /*

      //   //For standard quadratic elements all nodes are edge nodes
      //   unsigned n_vertex_and_edge_node = n_node;
      //   //If we have enriched elements, there are only 10 vertex and edge
      //   nodes if(n_node==15)
      //    {
      //     //There are only 10 vertex and edge nodes
      //     n_vertex_and_edge_node = 10;
      //    }

      //   // Loop over the new (edge) nodes in the element and create them.
      //   for(unsigned j=4;j<n_vertex_and_edge_node;j++)
      //    {

      //     // Figure out edge nodes
      //     switch(j)
      //      {

      //       // Node 4 is between nodes 0 and 1
      //      case 4:

      //       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(0);
      //       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(1);
      //       break;


      //       // Node 5 is between nodes 0 and 2
      //      case 5:

      //       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(0);
      //       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(2);
      //       break;

      //       // Node 6 is between nodes 0 and 3
      //      case 6:

      //       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(0);
      //       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(3);
      //       break;

      //       // Node 7 is between nodes 1 and 2
      //      case 7:

      //       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(1);
      //       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(2);
      //       break;

      //       // Node 8 is between nodes 2 and 3
      //      case 8:

      //       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(2);
      //       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(3);
      //       break;

      //       // Node 9 is between nodes 1 and 3
      //      case 9:

      //       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(1);
      //       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(3);
      //       break;

      //       break;

      //       //Error
      //      default:

      //       throw OomphLibError("More than ten edge nodes in Tet Element",
      //       OOMPH_CURRENT_FUNCTION,
      //                           OOMPH_EXCEPTION_LOCATION);
      //      }


      //     // Do we need a boundary node?
      //     bool need_boundary_node=false;

      //     // hierher At some point fine tune via intersection of boundary
      //     sets if
      //     (edge_node1_pt->is_on_boundary()&&edge_node2_pt->is_on_boundary())
      //      {
      //       need_boundary_node=true;
      //      }

      //     // Do we need a new node?
      //     if (central_edge_node_pt(edge_node1_pt,edge_node2_pt)==0)
      //      {
      //       Node* new_node_pt=0;

      //       // Create new  boundary node
      //       if (need_boundary_node)
      //        {
      //         new_node_pt=finite_element_pt(e)->
      //          construct_boundary_node(j,time_stepper_pt);
      //        }
      //       // Create new normal node
      //       else
      //        {
      //         new_node_pt=finite_element_pt(e)->
      //          construct_node(j,time_stepper_pt);
      //        }
      //       Node_pt.push_back(new_node_pt);

      //       // Now indicate existence of newly created mideside node in map
      //       central_edge_node_pt(edge_node1_pt,edge_node2_pt)=new_node_pt;
      //       central_edge_node_pt(edge_node2_pt,edge_node1_pt)=new_node_pt;

      //       // What are the node's local coordinates?
      //       finite_element_pt(e)->local_coordinate_of_node(j,s);

      //       // Find the coordinates of the new node from the existing
      //       // and fully-functional element in the scaffold mesh
      //       for(unsigned i=0;i<dim;i++)
      //        {
      //         new_node_pt->x(i)=
      //          Tmp_mesh_pt->finite_element_pt(e)->interpolated_x(s,i);
      //        }
      //      }
      //     else
      //      {
      //       // Set pointer to the existing node
      //       finite_element_pt(e)->node_pt(j)=
      //        central_edge_node_pt(edge_node1_pt,edge_node2_pt);
      //      }

      //    } // end of loop over new nodes

      //   //Need to sort out the face nodes
      //   if(n_node==15)
      //    {
      //     // Loop over the new (face) nodes in the element and create them.
      //     for(unsigned j=10;j<14;j++)
      //      {
      //       //Empty the set of face nodes
      //       face_nodes_pt.clear();
      //       // Figure out face nodes
      //       switch(j)
      //        {

      //         // Node 10 is between nodes 0 and 1 and 3
      //        case 10:

      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(0));
      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(1));
      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(3));
      //         break;

      //         // Node 11 is between nodes 0 and 1 and 2
      //        case 11:

      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(0));
      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(1));
      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(2));
      //         break;

      //         // Node 12 is between nodes 0 and 2 and 3
      //        case 12:

      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(0));
      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(2));
      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(3));
      //         break;

      //         // Node 13 is between nodes 1 and 2 and 3
      //        case 13:

      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(1));
      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(2));
      //         face_nodes_pt.insert(Tmp_mesh_pt->finite_element_pt(e)->node_pt(3));
      //         break;


      //         //Error
      //        default:

      //         throw OomphLibError("More than four face nodes in Tet Element",
      //         OOMPH_CURRENT_FUNCTION,
      //                             OOMPH_EXCEPTION_LOCATION);
      //        }

      //       // Do we need a boundary node?
      //       bool need_boundary_node=false;

      //       //Work it out from the face boundary
      //       boundary_id = Tmp_mesh_pt->face_boundary(e,face_map[j-10]);
      //       //If it's non-zero then we do need to create a boundary node
      //       if(boundary_id!=0) {need_boundary_node=true;}

      //       // Do we need a new node?
      //       if (central_face_node_pt[face_nodes_pt]==0)
      //        {
      //         Node* new_node_pt=0;

      //         // Create a new  boundary node
      //         if (need_boundary_node)
      //          {
      //           new_node_pt=finite_element_pt(e)->
      //            construct_boundary_node(j,time_stepper_pt);
      //           //Add it to the boundary
      //           add_boundary_node(boundary_id-1,new_node_pt);
      //          }
      //         // Create new normal node
      //         else
      //          {
      //           new_node_pt=finite_element_pt(e)->
      //            construct_node(j,time_stepper_pt);
      //          }
      //         Node_pt.push_back(new_node_pt);

      //       // Now indicate existence of newly created mideside node in map
      //       central_face_node_pt[face_nodes_pt]=new_node_pt;

      //       // What are the node's local coordinates?
      //       finite_element_pt(e)->local_coordinate_of_node(j,s);

      //       // Find the coordinates of the new node from the existing
      //       // and fully-functional element in the scaffold mesh
      //       for(unsigned i=0;i<dim;i++)
      //        {
      //         new_node_pt->x(i)=
      //          Tmp_mesh_pt->finite_element_pt(e)->interpolated_x(s,i);
      //        }
      //        }
      //       else
      //        {
      //         // Set pointer to the existing node
      //         finite_element_pt(e)->node_pt(j)=
      //          central_face_node_pt[face_nodes_pt];
      //        }
      //      } //End of loop over face nodes

      //     //Construct the element's central node, which is not on a boundary
      //     {
      //      Node* new_node_pt=
      //       finite_element_pt(e)->construct_node(14,time_stepper_pt);
      //      Node_pt.push_back(new_node_pt);

      //      // What are the node's local coordinates?
      //      finite_element_pt(e)->local_coordinate_of_node(14,s);
      //      // Find the coordinates of the new node from the existing
      //      // and fully-functional element in the scaffold mesh
      //      for(unsigned i=0;i<dim;i++)
      //       {
      //        new_node_pt->x(i)=
      //         Tmp_mesh_pt->finite_element_pt(e)->interpolated_x(s,i);
      //       }
      //     }
      //    } //End of enriched case

      //  } //end of loop over elements


      // //Boundary conditions

      // // Matrix map to check if a node has already been added to
      // // the boundary number b
      // MapMatrixMixed<Node*,unsigned,bool> bound_node_done;

      // // Loop over elements
      // for (unsigned e=0;e<nelem;e++)
      //  {
      //   // Loop over new local nodes
      //   for(unsigned j=4;j<n_node;j++)
      //    {
      //     // Pointer to the element's local node
      //     Node* loc_node_pt=finite_element_pt(e)->node_pt(j);

      //     // This will have to be changed for higher-order elements
      //     //=======================================================

      //     //  These are the face nodes on the element's face 0:
      //     if ( (j==4) || (j==5) || (j==7) )
      //      {
      //       boundary_id=Tmp_mesh_pt->face_boundary(e,0);
      //       if (boundary_id!=0)
      //        {
      //         if (!bound_node_done(loc_node_pt,boundary_id-1))
      //          {
      //           add_boundary_node(boundary_id-1,loc_node_pt);
      //           bound_node_done(loc_node_pt,boundary_id-1)=true;
      //          }
      //        }
      //      }


      //     // These are the face nodes on the element's face 1:
      //     if ( (j==4) || (j==6) || (j==9) )
      //      {
      //       boundary_id=Tmp_mesh_pt->face_boundary(e,1);
      //       if (boundary_id!=0)
      //        {
      //         if (!bound_node_done(loc_node_pt,boundary_id-1))
      //          {
      //           add_boundary_node(boundary_id-1,loc_node_pt);
      //           bound_node_done(loc_node_pt,boundary_id-1)=true;
      //          }
      //        }
      //      }

      //     // These are the face nodes on the element's face 2:
      //     if ( (j==5) || (j==6) || (j==8) )
      //      {
      //       boundary_id=Tmp_mesh_pt->face_boundary(e,2);
      //       if (boundary_id!=0)
      //        {
      //         if (!bound_node_done(loc_node_pt,boundary_id-1))
      //          {
      //           add_boundary_node(boundary_id-1,loc_node_pt);
      //           bound_node_done(loc_node_pt,boundary_id-1)=true;
      //          }
      //        }
      //      }


      //     // These are the face nodes on the element's face 3:
      //     if  ( (j==7) || (j==8) || (j==9) )
      //      {
      //       boundary_id=Tmp_mesh_pt->face_boundary(e,3);
      //       if (boundary_id!=0)
      //        {
      //         if (!bound_node_done(loc_node_pt,boundary_id-1))
      //          {
      //           add_boundary_node(boundary_id-1,loc_node_pt);
      //           bound_node_done(loc_node_pt,boundary_id-1)=true;
      //          }
      //        }
      //      }

      //    } //end for j

      // */

    } // end for e


  } // end function

  //=====================================================================
  /// Helper function to set up the reverse look up scheme for facets.
  /// This is used to set up boundary coordinates.
  //=====================================================================
  template<class ELEMENT>
  void TetgenMesh<ELEMENT>::setup_reverse_lookup_schemes_for_faceted_surface(
    TetMeshFacetedSurface* const& faceted_surface_pt)
  {
    // Set up reverse lookup scheme for a given faceted surface
    // Assumption is that there there is one boundary ID per facet.
    unsigned n_facet = faceted_surface_pt->nfacet();
    for (unsigned f = 0; f < n_facet; f++)
    {
      unsigned b = faceted_surface_pt->one_based_facet_boundary_id(f);
      if (b != 0)
      {
        this->Tet_mesh_faceted_surface_pt[b - 1] = faceted_surface_pt;
        this->Tet_mesh_facet_pt[b - 1] = faceted_surface_pt->facet_pt(f);
      }
      else
      {
        std::ostringstream error_message;
        error_message << "Boundary IDs have to be one-based. Yours is " << b
                      << "\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
  }


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //=========================================================================
  /// Transfer tetgenio data from the input to the output
  /// The output is assumed to have been constructed and "empty"
  //=========================================================================
  template<class ELEMENT>
  void TetgenMesh<ELEMENT>::deep_copy_of_tetgenio(tetgenio* const& input_pt,
                                                  tetgenio*& output_pt)
  {
    // Copy data values
    output_pt->firstnumber = input_pt->firstnumber;
    output_pt->mesh_dim = input_pt->mesh_dim;
    output_pt->useindex = input_pt->useindex;

    // Copy the number of points
    output_pt->numberofpoints = input_pt->numberofpoints;
    output_pt->numberofpointattributes = input_pt->numberofpointattributes;
    output_pt->numberofpointmtrs = input_pt->numberofpointmtrs;

    const int n_point = output_pt->numberofpoints;
    if (n_point > 0)
    {
      output_pt->pointlist = new double[n_point * 3];
      // Copy the values
      for (int n = 0; n < (n_point * 3); ++n)
      {
        output_pt->pointlist[n] = input_pt->pointlist[n];
      }

      // If there are point markers copy as well
      if (input_pt->pointmarkerlist != (int*)NULL)
      {
        output_pt->pointmarkerlist = new int[n_point];
        for (int n = 0; n < n_point; ++n)
        {
          output_pt->pointmarkerlist[n] = input_pt->pointmarkerlist[n];
        }
      }

      const int n_attr = output_pt->numberofpointattributes;
      if (n_attr > 0)
      {
        output_pt->pointattributelist = new double[n_point * n_attr];
        for (int n = 0; n < (n_point * n_attr); ++n)
        {
          output_pt->pointattributelist[n] = input_pt->pointattributelist[n];
        }
      }
    } // End of point case

    // Now add in metric tensors (if there are any)
    const int n_point_mtr = output_pt->numberofpointmtrs;
    if (n_point_mtr > 0)
    {
      output_pt->pointmtrlist = new double[n_point * n_point_mtr];
      for (int n = 0; n < (n_point * n_point_mtr); ++n)
      {
        output_pt->pointmtrlist[n] = input_pt->pointmtrlist[n];
      }
    }

    // Next tetrahedrons
    output_pt->numberoftetrahedra = input_pt->numberoftetrahedra;
    output_pt->numberofcorners = input_pt->numberofcorners;
    output_pt->numberoftetrahedronattributes =
      input_pt->numberoftetrahedronattributes;

    const int n_tetra = output_pt->numberoftetrahedra;
    const int n_corner = output_pt->numberofcorners;
    if (n_tetra > 0)
    {
      output_pt->tetrahedronlist = new int[n_tetra * n_corner];
      for (int n = 0; n < (n_tetra * n_corner); ++n)
      {
        output_pt->tetrahedronlist[n] = input_pt->tetrahedronlist[n];
      }

      // Add in the volume constriants
      if (input_pt->tetrahedronvolumelist != (double*)NULL)
      {
        output_pt->tetrahedronvolumelist = new double[n_tetra];
        for (int n = 0; n < n_tetra; ++n)
        {
          output_pt->tetrahedronvolumelist[n] =
            input_pt->tetrahedronvolumelist[n];
        }
      }

      // Add in the attributes
      const int n_tetra_attr = output_pt->numberoftetrahedronattributes;
      if (n_tetra_attr > 0)
      {
        output_pt->tetrahedronattributelist =
          new double[n_tetra * n_tetra_attr];
        for (int n = 0; n < (n_tetra * n_tetra_attr); ++n)
        {
          output_pt->tetrahedronattributelist[n] =
            input_pt->tetrahedronattributelist[n];
        }
      }

      // Add in the neighbour list
      if (input_pt->neighborlist != (int*)NULL)
      {
        output_pt->neighborlist = new int[n_tetra * 4];
        for (int n = 0; n < (n_tetra * 4); ++n)
        {
          output_pt->neighborlist = input_pt->neighborlist;
        }
      }
    } // End of tetrahedron section

    // Now arrange the facet list
    output_pt->numberoffacets = input_pt->numberoffacets;
    const int n_facet = output_pt->numberoffacets;
    if (n_facet > 0)
    {
      output_pt->facetlist = new tetgenio::facet[n_facet];
      for (int n = 0; n < n_facet; ++n)
      {
        tetgenio::facet* input_f_pt = &input_pt->facetlist[n];
        tetgenio::facet* output_f_pt = &output_pt->facetlist[n];

        // Copy polygons and holes from the facets
        output_f_pt->numberofpolygons = input_f_pt->numberofpolygons;

        // Loop over polygons
        const int n_poly = output_f_pt->numberofpolygons;
        if (n_poly > 0)
        {
          output_f_pt->polygonlist = new tetgenio::polygon[n_poly];
          for (int p = 0; p < n_poly; ++p)
          {
            tetgenio::polygon* output_p_pt = &output_f_pt->polygonlist[p];
            tetgenio::polygon* input_p_pt = &input_f_pt->polygonlist[p];
            // Find out how many vertices each polygon has
            output_p_pt->numberofvertices = input_p_pt->numberofvertices;
            // Now copy of the vertices
            const int n_vertex = output_p_pt->numberofvertices;
            if (n_vertex > 0)
            {
              output_p_pt->vertexlist = new int[n_vertex];
              for (int v = 0; v < n_vertex; ++v)
              {
                output_p_pt->vertexlist[v] = input_p_pt->vertexlist[v];
              }
            }
          }
        }

        // Hole information
        output_f_pt->numberofholes = input_f_pt->numberofholes;
        const int n_hole = output_f_pt->numberofholes;
        if (n_hole > 0)
        {
          throw OomphLibError("Don't know how to handle holes yet\n",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      } // end of loop over facets

      // Add the facetmarkers if there are any
      if (input_pt->facetmarkerlist != (int*)NULL)
      {
        output_pt->facetmarkerlist = new int[n_facet];
        for (int n = 0; n < n_facet; ++n)
        {
          output_pt->facetmarkerlist[n] = input_pt->facetmarkerlist[n];
        }
      }
    }

    // Now the holes
    output_pt->numberofholes = input_pt->numberofholes;
    const int n_hole = output_pt->numberofholes;
    if (n_hole > 0)
    {
      output_pt->holelist = new double[n_hole * 3];
      for (int n = 0; n < (n_hole * 3); ++n)
      {
        output_pt->holelist[n] = input_pt->holelist[n];
      }
    }

    // Now the regions
    output_pt->numberofregions = input_pt->numberofregions;
    const int n_region = output_pt->numberofregions;
    if (n_region > 0)
    {
      output_pt->regionlist = new double[n_region * 5];
      for (int n = 0; n < (n_region * 5); ++n)
      {
        output_pt->regionlist[n] = input_pt->regionlist[n];
      }
    }

    // Now the facet constraints
    output_pt->numberoffacetconstraints = input_pt->numberoffacetconstraints;
    const int n_facet_const = output_pt->numberoffacetconstraints;
    if (n_facet_const > 0)
    {
      output_pt->facetconstraintlist = new double[n_facet_const * 2];
      for (int n = 0; n < (n_facet_const * 2); ++n)
      {
        output_pt->facetconstraintlist[n] = input_pt->facetconstraintlist[n];
      }
    }

    // Now the segment constraints
    output_pt->numberofsegmentconstraints =
      input_pt->numberofsegmentconstraints;
    const int n_seg_const = output_pt->numberofsegmentconstraints;
    if (n_seg_const > 0)
    {
      output_pt->segmentconstraintlist = new double[n_seg_const * 2];
      for (int n = 0; n < (n_seg_const * 2); ++n)
      {
        output_pt->segmentconstraintlist[n] =
          input_pt->segmentconstraintlist[n];
      }
    }

    // Now the face list
    output_pt->numberoftrifaces = input_pt->numberoftrifaces;
    const int n_tri_face = output_pt->numberoftrifaces;
    if (n_tri_face > 0)
    {
      output_pt->trifacelist = new int[n_tri_face * 3];
      for (int n = 0; n < (n_tri_face * 3); ++n)
      {
        output_pt->trifacelist[n] = input_pt->trifacelist[n];
      }

      output_pt->trifacemarkerlist = new int[n_tri_face];
      for (int n = 0; n < n_tri_face; ++n)
      {
        output_pt->trifacemarkerlist[n] = input_pt->trifacemarkerlist[n];
      }
    }

    // Now the edge list
    output_pt->numberofedges = input_pt->numberofedges;
    const int n_edge = output_pt->numberofedges;
    if (n_edge > 0)
    {
      output_pt->edgelist = new int[n_edge * 2];
      for (int n = 0; n < (n_edge * 2); ++n)
      {
        output_pt->edgelist[n] = input_pt->edgelist[n];
      }

      output_pt->edgemarkerlist = new int[n_edge];
      for (int n = 0; n < n_edge; ++n)
      {
        output_pt->edgemarkerlist[n] = input_pt->edgemarkerlist[n];
      }
    }
  }


} // namespace oomph
#endif
