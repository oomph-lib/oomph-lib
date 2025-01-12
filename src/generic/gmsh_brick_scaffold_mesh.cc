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

#include "gmsh_brick_scaffold_mesh.h"


namespace oomph
{
    /// Hexaheadral constituents where each Hexa is consit of 8 nodes
    /// 6 faces/quads and 12 edges.
    enum Hexaheadral{
        nNodes = 8,
        nFaces = 6,
        nEdges = 12
    };
    //======================================================================
    /// Constructor: Pass the filename of the hexahedral *.msh file
    /// The assumptions are that the nodes have been assigned boundary
    /// information which is used in the nodal construction to make sure
    /// that BoundaryNodes are constructed. Any additional boundaries are
    /// determined from the face boundary information.
    //======================================================================
    GmshBrickScaffoldMesh::GmshBrickScaffoldMesh(const std::string& file_name, bool verbose)
    {
        // Process the element file
        // --------------------------
        GMSH::Gmsh rg(file_name, verbose);
        Boundaries = rg.boundaries;
        Ordered_boundary_conditions = rg.orderedBC;
        const std::vector<GMSH::Node>& nodes = rg.getNodes();
        std::vector<GMSH::Edge> edges = rg.getEdges();

        // dimension of the mesh
        const unsigned dim = rg.getMeshDim();
        std::vector<GMSH::Quad> faces = rg.getQuads();
        std::vector<GMSH::Hexa> bricks= rg.getHexas();


        unsigned n_local_face = rg.getFacesPerElement();
        unsigned n_local_edge = rg.getEdgesPerElement();


        if (dim ==2){
            std::ostringstream error_stream;
            error_stream << "Mesh is 2D!" <<std::endl;
            
            throw OomphLibError(
                    error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
        
        // Read in number of elements
        unsigned n_element = bricks.size();

        // Read in number of nodes per element // this should be constant = 8
        unsigned n_local_node = rg.getNodesPerElement();


        // Throw an error if we have anything other than [Hexas or Quads] but 
        // linear simplices
        if ((n_local_node != Hexaheadral::nNodes) && (n_local_node != 4))
        {
            std::ostringstream error_stream;
            error_stream
                    << "Gmsh Reader should only be used to generate 8-noded hexahedra!\n"
                    << "Your gmsh input file, contains data for " << n_local_node
                    << "-noded hexahedra" << std::endl;

            throw OomphLibError(
                    error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }

        // Element attributes may be used to distinguish internal regions
        // NOTE: This stores doubles because gmsh forces us to!
        // Element_attribute.resize(n_element, 0.0);


        // Resize storage for the global node numbers listed element-by-element
        Global_node.resize(n_element*n_local_node);

        // Initialise (global) node counter
        unsigned k = 0;

        // Read in global node numbers
        for (auto &brick: bricks)
        {
            for (auto node: brick.nodesTag) {
                Global_node[k] = node;
                k++;
            }
        }


        // Resize the Element vector
        Element_pt.resize(n_element);

        // Process node file
        //--------------------

        // Read in the number of nodes
        unsigned n_node = nodes.size();

        // Create a vector of boolean so as not to create the same node twice
        std::vector<bool> done(n_node, false);

        // Resize the Node vector
        Node_pt.resize(n_node);


        // Create storage for nodal positions and boundary markers
        Vector<double> x_node(n_node);
        Vector<double> y_node(n_node);
        Vector<double> z_node(n_node);
        Vector<Vector<int>>    bound(n_node);

        // Check if there are attributes
        // If not
        for (unsigned i = 0; i < n_node; i++)
        {
            x_node[i] = nodes[i].coord[0];
            y_node[i] = nodes[i].coord[1];
            z_node[i] = nodes[i].coord[2];

            bound[i].resize(nodes[i].boundaries.size());
            int c = 0;
            for (auto it = nodes[i].boundaries.begin(); 
                 it != nodes[i].boundaries.end(); ++it) 
	    {
                bound[i][c] = *it;
                c++;
            }

        }

        // Determine the highest boundary index
        //------------------------------------
        unsigned n_bound = rg.getBoundSize();

        // Process face to extract boundary faces
        //--------------------------------------------

        // Number of faces in faces
        unsigned n_face = faces.size();


        // Storage for the global node numbers (in the gmsh 1-based
        // numbering scheme!) of the first, second and third  node in
        // each segment
        Vector<unsigned> first_node(n_face);
        Vector<unsigned> second_node(n_face);
        Vector<unsigned> third_node(n_face);
        Vector<unsigned> forth_node(n_face);

        // Storage for the boundary marker for each face
        Vector<int> face_boundary(n_face);

        // Dummy for global face number
        unsigned dummy_face_number;

        // Storage for the (boundary) faces associated with each node.
        // Nodes are indexed using Reader's 1-based scheme, which is why
        // there is a +1 here
        Vector<std::set<unsigned>> node_on_faces(n_node);

        for (int i = 0; i < faces.size(); ++i)
        {
            first_node[i]  = faces[i].nodesTag[0];
            second_node[i] = faces[i].nodesTag[1];
            third_node[i]  = faces[i].nodesTag[2];
            forth_node[i]  = faces[i].nodesTag[3];

            // \todo: what if more than b.c.?
            auto it = faces[i].boundaries.begin();
            face_boundary[i] = *it;

            // Add the face index to each node
            node_on_faces[first_node[i]].insert(i);
            node_on_faces[second_node[i]].insert(i);
            node_on_faces[third_node[i]].insert(i);
            node_on_faces[forth_node[i]].insert(i);
        }

        // Set number of boundaries
        if (n_bound > 0)
        {
            this->set_nboundary(n_bound);
        }

        // Translate enumeration
        Vector<unsigned> gmsh_to_oomph_node(8);
        gmsh_to_oomph_node[0] = 0;
        gmsh_to_oomph_node[1] = 1;
        gmsh_to_oomph_node[2] = 3;
        gmsh_to_oomph_node[3] = 2;

        gmsh_to_oomph_node[4] = 4;
        gmsh_to_oomph_node[5] = 5;
        gmsh_to_oomph_node[6] = 7;
        gmsh_to_oomph_node[7] = 6;

        // Create the elements
        unsigned counter = 0;
        for (unsigned e = 0; e < n_element; e++)
        {
            Element_pt[e] = new QElement<3, 2>;

            // Loop over the other nodes
            for (unsigned j = 0; j < n_local_node; j++)
            {
                unsigned global_node_number = Global_node[counter];
                if (!done[global_node_number])
                {
                    // If we're on a boundary
                    // [gmsh convention 1 or greater is a boundary, 
                    // 0 not a boundary]
                    for (int i = 0; i < bound[global_node_number].size(); ++i)
                    {
                        if (bound[global_node_number][i] > 0)
                        {
                            // Construct the boundary node
                            Node_pt[global_node_number] = 
                            finite_element_pt(e)->construct_boundary_node(
                                                        gmsh_to_oomph_node[j]);

                            // Add to the boundary lookup scheme.
                            // bound[global_node_number]-1 because vector 
                            // indexing starts from 0
                            add_boundary_node(bound[global_node_number][i] - 1,
                                                  Node_pt[global_node_number]);
                        }
                        else
                        {
                            Node_pt[global_node_number] = 
                            finite_element_pt(e)->construct_node(
                                                        gmsh_to_oomph_node[j]);
                        }

                        done[global_node_number] = true;
                        Node_pt[global_node_number]->x(0) = 
                        			  x_node[global_node_number];
                        Node_pt[global_node_number]->x(1) = 
                                                  y_node[global_node_number];
                        Node_pt[global_node_number]->x(2) = 
                                                  z_node[global_node_number];
                    }

                }
                else
                {
                    // Otherwise copy the pointer over
                    finite_element_pt(e)->node_pt(gmsh_to_oomph_node[j]) = 
                    				Node_pt[global_node_number];
                }
                counter++;
            }
        }


        // Resize the "matrix" that stores the boundary id for each
        // face in each element.
        Face_boundary.resize(n_element);
        Face_index.resize(n_element);
        Edge_index.resize(n_element);


        // 0-based index scheme used to construct a global lookup for
        // each face that will be used to uniquely construct interior
        // face nodes.
        // The faces that lie on boundaries will have already been allocated so
        // we can start from there
        for (auto& face: faces) {
            if (*face.boundaries.begin() == 0) Nglobal_face++;
        }


        // Storage for the edges associated with each node.
        // Nodes are indexed using Reader's 1-based scheme, which is why
        // there is a +1 here
        Vector<std::set<unsigned>> node_on_edges(n_node);

        // 0-based index scheme used to construct a global lookup for each
        // edge that will be used to uniquely construct interior edge nodes
        Nglobal_edge = edges.size();

        // Conversion from the edge numbers to the nodes at the end
        // of each edge

        // Loop over the elements
        for (unsigned e = 0; e < n_element; e++)
        {
            // Each element has six faces
            Face_boundary[e].resize(n_local_face);
            Face_index[e].resize(n_local_face);
            // By default, each face is NOT on a boundary
            for (unsigned i = 0; i < n_local_face; i++)
            {
                Face_boundary[e][i] = 0;
            }

            // Read out the global node numbers of the nodes from
            // tetgen's 1-based numbering.
            // The offset is to match the offset used above
            const unsigned element_offset = e * n_local_node;
            unsigned glob_num[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            for (unsigned i = 0; i < n_local_node; ++i)
            {
                glob_num[i] = Global_node[element_offset + i];
            }

            // Now we know the global node numbers of the elements' four nodes
            // in tetgen's 1-based numbering.

            // Determine whether any of the element faces have already been allocated
            // an index. This may be because they are located on boundaries, or have
            // already been visited The global index for the i-th face will be stored
            // in Face_index[i]

            GMSH::Hexa& brick = bricks[e];
            unsigned n_local_faces = brick.quadsTag.size();
            
            // Loop over the local faces in the element
            for (unsigned i = 0; i < n_local_faces; ++i)
            {
                //std::vector<int> facesTags = brick.quadsTag;

                // Otherwise we already have a face
                GMSH::Quad& face = faces[brick.quadsTag[i]];

                auto it_b = face.boundaries.begin();
                if (*it_b > 0 )
                {
                    // Set the face index
                    Face_index[e][i] = face.quadTag;

                    // Allocate the boundary index, if it's a boundary
                    Face_boundary[e][i] = *it_b;

                    // Add the nodes to the boundary look-up scheme in
                    // oomph-lib (0-based) index
                    for (int nodeTag : face.nodesTag)
                    {
                        auto b = *it_b - 1;
                        // Don't add the omitted node
                        add_boundary_node(b, Node_pt[nodeTag]);
                    }
                }
                else
                {
                    // Allocate the next global index
                    Face_index[e][i] = face.quadTag;

                }
            } // end of loop over the faces


            Edge_index[e].resize(n_local_edge);

            // Loop over the element edges and assign global edge numbers 
            // we have 12 edges on 1 hexa element
            for (unsigned i = 0; i < brick.edgesTag.size(); ++i)
            {
                std::vector<int> edgesTags = brick.edgesTag;

                GMSH::Edge& edge = edges[edgesTags[i]];

                // Allocate the next global index
                Edge_index[e][i] = edge.edgeTag;

                // Otherwise we already have an edge

                // Set the edge index from the result of the intersection
                Edge_index[e][i] = edge.edgeTag;

            } // end for edge

        } // end for e

        // Now determine whether any edges lie on boundaries by using the
        // face boundary scheme and

        // Resize the storage
        Edge_boundary.resize(Nglobal_edge, false);

        // Now loop over all the boundary faces and mark that all edges
        // must also lie on the boundary
        for (unsigned i = 0; i < edges.size(); ++i)
        {
            if (*edges[i].boundaries.begin() > 0)
            {
                Edge_boundary[i] = true;
            }
        }
    } // end of constructor


} // namespace oomph
