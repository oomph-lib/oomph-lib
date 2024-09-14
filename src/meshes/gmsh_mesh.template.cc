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
#ifndef OOMPH_GEOMPACK_MESH_TEMPLATE_CC
#define OOMPH_GEOMPACK_MESH_TEMPLATE_CC

#include "gmsh_mesh.template.h"


namespace oomph
{
    /// ////////////////////////////////////////////////////////////////////
    /// ////////////////////////////////////////////////////////////////////
    /// ////////////////////////////////////////////////////////////////////


    //========================================================================
    /// Build unstructured tet mesh based on output from scaffold
    //========================================================================
    template<class ELEMENT>
    void GmshMesh<ELEMENT>::build_from_scaffold(TimeStepper* time_stepper_pt)
    {
        // Mesh can only be built with 3D Telements.
        MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

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
        /*if (use_attributes)
        {
            Boundary_region_element_pt.resize(nbound);
            Face_index_region_at_boundary.resize(nbound);
        }*/

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
                    if (boundaries_pt != nullptr)
                    {
                        // Create new boundary node
                        Node* new_node_pt = finite_element_pt(e)->construct_boundary_node(j, time_stepper_pt);

                        // Give it a number (not necessarily the global node
                        // number in the scaffold mesh -- we just need something
                        // to keep track...)
                        global_count++;
                        global_number[scaffold_node_pt] = global_count;

                        // Add to boundaries
                        for(auto it = boundaries_pt->begin(); it != boundaries_pt->end();
                             ++it)
                        {
                            auto b = *it;
                            add_boundary_node(b, new_node_pt);
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
            /*if (false)
            {
                element_attribute_map[Tmp_mesh_pt->element_attribute(e)].push_back(
                        finite_element_pt(e));
            }*/
        }

        // Now let's construct lists
        // Find the number of attributes
        /*if (use_attributes)
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
        }*/

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
         const unsigned faces_on_edge[6][2] = {{0, 1}, {0, 2},{1, 2},
                                               {0, 3}, {2, 3}, {1, 3}};

        // Loop over all elements
        for (unsigned e = 0; e < nelem; e++)
        {
            // Cache pointers to the elements
            FiniteElement* const elem_pt = this->finite_element_pt(e);
            FiniteElement* const tmp_elem_pt = Tmp_mesh_pt->finite_element_pt(e);

            // The number of edge nodes is 4 + 6*(n_node1d-2)
            unsigned n_edge_node = 8 + 6 * (n_node_1d - 2);

            // Now loop over the edge nodes
            // assuming that the numbering scheme is the same as the triangles
            // which puts edge nodes in ascending order.
            // We don't have any higher than quadratic at the moment, so it's
            // a bit academic really.

            // Only bother if we actually have some edge nodes
            if (n_edge_node > 8)
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
                        unsigned face_boundary_id = Tmp_mesh_pt->face_boundary(e, faces_on_edge[j][i]);
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
            for (unsigned j = 0; j < 6; ++j)
            {
                // Find the boundary id of the face
                boundary_id = Tmp_mesh_pt->face_boundary(e, j);


                if (boundary_id > 0)
                {
                    Boundary_element_pt[boundary_id - 1].push_back(elem_pt);

                    // This face indexing is not correct hope oomph-lib dev can help here!
                    Vector<unsigned > gmsh_face_map(6,0);
                    gmsh_face_map[0] = -2;
                    gmsh_face_map[1] = -1;
                    gmsh_face_map[2] = -3;
                    gmsh_face_map[3] =  1;
                    gmsh_face_map[4] =  2;
                    gmsh_face_map[5] =  3;

                    // Tetgen Face 0 is our Face 3
                    // Tetgen Face 1 is our Face 2
                    // Tetgen Face 2 is our Face 1
                    // Tetgen Face 3 is our Face 0

                    Face_index_at_boundary[boundary_id - 1].push_back(gmsh_face_map[j]);

                    // If using regions set up the boundary information
                    /*if (use_attributes)
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
                    }*/
                }
            } // End of loop over faces


            // Lookup scheme has now been setup
            Lookup_for_elements_next_boundary_is_setup = true;


        } // end for e


    } // end function


} // namespace oomph
#endif

