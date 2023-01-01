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
#ifndef OOMPH_THIN_LAYER_BRICK_ON_TET_MESH_TEMPLATE_CC
#define OOMPH_THIN_LAYER_BRICK_ON_TET_MESH_TEMPLATE_CC


#include "../solid/solid_elements.h"
#include "thin_layer_brick_on_tet_mesh.template.h"


namespace oomph
{
  //=====================================================================
  /// Constructor: Specify (quadratic) tet mesh, boundary IDs of
  /// boundary on which the current mesh is to be erected (in an FSI context
  /// this boundary tends to be the FSI boundary of the fluid mesh. Also
  /// specify the uniform thickness of layer, and the number of element layers.
  /// The vectors stored in in_out_boundary_ids contain the boundary
  /// IDs of the other boundaries in the tet mesh. In an FSI context
  /// these typically identify the in/outflow boundaries in the fluid
  /// mesh. The boundary enumeration of the current mesh follows the
  /// one of the underlying fluid mesh: The enumeration of the FSI boundary
  /// matches (to enable the setup of the FSI matching); the "in/outflow"
  /// faces in this mesh inherit the same enumeration as the in/outflow
  /// faces in the underlying fluid mesh. Finally, the "outer" boundary
  /// gets its own boundary ID.
  /// Timestepper defaults to steady pseudo-timestepper.
  //=====================================================================
  template<class ELEMENT>
  ThinLayerBrickOnTetMesh<ELEMENT>::ThinLayerBrickOnTetMesh(
    Mesh* tet_mesh_pt,
    const Vector<unsigned>& boundary_ids,
    ThicknessFctPt thickness_fct_pt,
    const unsigned& nlayer,
    const Vector<Vector<unsigned>>& in_out_boundary_id,
    TimeStepper* time_stepper_pt)
    : Thickness_fct_pt(thickness_fct_pt)
  {
    // Mesh can only be built with 3D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

    // Figure out if the tet mesh is a solid mesh
    bool tet_mesh_is_solid_mesh = false;
    if (dynamic_cast<SolidFiniteElement*>(tet_mesh_pt->element_pt(0)) != 0)
    {
      tet_mesh_is_solid_mesh = true;
    }

    // Setup lookup scheme for local coordinates on triangular faces.
    // The local coordinates identify the points on the triangular
    // FaceElements on which we place the bottom layer of the
    // brick nodes.
    Vector<Vector<double>> s_face(19);
    for (unsigned i = 0; i < 19; i++)
    {
      s_face[i].resize(2);

      switch (i)
      {
          // Vertex nodes

        case 0:
          s_face[i][0] = 1.0;
          s_face[i][1] = 0.0;
          break;

        case 1:
          s_face[i][0] = 0.0;
          s_face[i][1] = 1.0;
          break;

        case 2:
          s_face[i][0] = 0.0;
          s_face[i][1] = 0.0;
          break;

          // Midside nodes

        case 3:
          s_face[i][0] = 0.5;
          s_face[i][1] = 0.5;
          break;

        case 4:
          s_face[i][0] = 0.0;
          s_face[i][1] = 0.5;
          break;


        case 5:
          s_face[i][0] = 0.5;
          s_face[i][1] = 0.0;
          break;


          // Quarter side nodes

        case 6:
          s_face[i][0] = 0.75;
          s_face[i][1] = 0.25;
          break;

        case 7:
          s_face[i][0] = 0.25;
          s_face[i][1] = 0.75;
          break;

        case 8:
          s_face[i][0] = 0.0;
          s_face[i][1] = 0.75;
          break;

        case 9:
          s_face[i][0] = 0.0;
          s_face[i][1] = 0.25;
          break;

        case 10:
          s_face[i][0] = 0.25;
          s_face[i][1] = 0.0;
          break;

        case 11:
          s_face[i][0] = 0.75;
          s_face[i][1] = 0.0;
          break;

          // Central node

        case 12:
          s_face[i][0] = 1.0 / 3.0;
          s_face[i][1] = 1.0 / 3.0;
          break;


          // Vertical internal midside nodes connecting 2 and 3

        case 13:
          s_face[i][0] = 5.0 / 24.0;
          s_face[i][1] = 5.0 / 24.0;
          break;

        case 14:
          s_face[i][0] = 5.0 / 12.0;
          s_face[i][1] = 5.0 / 12.0;
          break;

          // Internal midside nodes connecting nodes 0 and 4

        case 15:
          s_face[i][1] = 5.0 / 24.0;
          s_face[i][0] = 7.0 / 12.0; // 1.0-2.0*5.0/24.0;
          break;

        case 16:
          s_face[i][1] = 5.0 / 12.0;
          s_face[i][0] = 1.0 / 6.0; // 1.0-2.0*5.0/12.0;
          break;


          // Internal midside nodes connecting nodes 1 and 5

        case 17:
          s_face[i][0] = 5.0 / 24.0;
          s_face[i][1] = 7.0 / 12.0; // 1.0-2.0*5.0/24.0;
          break;

        case 18:
          s_face[i][0] = 5.0 / 12.0;
          s_face[i][1] = 1.0 / 6.0; // 1.0-2.0*5.0/12.0;
          break;
      }
    }


    // Translation scheme for inverted FaceElements
    MapMatrixMixed<int, unsigned, unsigned> translate;

    // Initialise with identify mapping
    for (unsigned i = 0; i < 19; i++)
    {
      translate(-1, i) = i;
      translate(1, i) = i;
    }
    translate(-1, 6) = 11;
    translate(-1, 11) = 6;
    translate(-1, 3) = 5;
    translate(-1, 5) = 3;
    translate(-1, 18) = 14;
    translate(-1, 14) = 18;
    translate(-1, 7) = 10;
    translate(-1, 10) = 7;
    translate(-1, 13) = 17;
    translate(-1, 17) = 13;
    translate(-1, 1) = 2;
    translate(-1, 2) = 1;
    translate(-1, 9) = 8;
    translate(-1, 8) = 9;

    // Lookup scheme relating "fluid" nodes to newly created "solid" nodes
    // (terminology for fsi problem)
    std::map<Node*, Node*> solid_node_pt;

    // Look up scheme for quarter edge nodes
    std::map<Edge, Node*> quarter_edge_node;

    // Map to store normal vectors for all surface nodes, labeled
    // by node on FSI surface
    std::map<Node*, Vector<Vector<double>>> normals;

    // Map of nodes connected to node on the tet surface, labeled by
    // node on FSI surface
    std::map<Node*, Vector<Node*>> connected_node_pt;

    // Number of elements in brick mesh
    Element_pt.reserve(3 * boundary_ids.size() * nlayer);

    // Get total number of distinct boundary IDs that we touch
    // in the fluid mesh
    std::set<unsigned> all_bnd;

    // Loop over all boundaries in tet mesh that make up the FSI interface
    unsigned nb = boundary_ids.size();
    for (unsigned ib = 0; ib < nb; ib++)
    {
      // Boundary number in "fluid" tet mesh
      unsigned b = boundary_ids[ib];

      // Loop over boundary nodes in the fluid mesh on that
      // boundary
      unsigned nnod = tet_mesh_pt->nboundary_node(b);
      for (unsigned j = 0; j < nnod; j++)
      {
        Node* nod_pt = tet_mesh_pt->boundary_node_pt(b, j);

        // Get pointer to set of boundaries this node is located on
        std::set<unsigned>* bnd_pt;
        nod_pt->get_boundaries_pt(bnd_pt);

        // Add
        for (std::set<unsigned>::iterator it = (*bnd_pt).begin();
             it != (*bnd_pt).end();
             it++)
        {
          all_bnd.insert(*it);
        }
      }
    }

    // Highest boundary ID
    unsigned highest_fluid_bound_id =
      *std::max_element(all_bnd.begin(), all_bnd.end());

    // Figure out which boundaries are actually on fsi boundary
    std::vector<bool> is_on_fsi_boundary(highest_fluid_bound_id + 1, false);
    for (unsigned ib = 0; ib < nb; ib++)
    {
      is_on_fsi_boundary[boundary_ids[ib]] = true;
    }


    // Figure out which boundaries are on the identified in/outflow boundaries
    unsigned n = in_out_boundary_id.size();
    Vector<std::vector<bool>> is_on_in_out_boundary(n);
    Vector<std::set<unsigned>> in_out_boundary_id_set(n);
    for (unsigned j = 0; j < n; j++)
    {
      is_on_in_out_boundary[j].resize(highest_fluid_bound_id + 1, false);
      unsigned nb = in_out_boundary_id[j].size();
      for (unsigned ib = 0; ib < nb; ib++)
      {
        is_on_in_out_boundary[j][in_out_boundary_id[j][ib]] = true;
      }
    }

    // Total number of boundaries: All boundaries that we touch
    // in the fluid mesh (the FSI boundary and the boundaries
    // on the in/outflow faces -- we flip these up and use
    // them for all the boundary faces in adjacent stacks
    // of solid elements) plus one additional boundary for
    // the outer boundary.
    unsigned maxb = highest_fluid_bound_id + 2;

    // Set number of boundaries
    set_nboundary(maxb);

    // Get ready for boundary lookup scheme
    Boundary_element_pt.resize(maxb);
    Face_index_at_boundary.resize(maxb);


    // Loop over all boundaries in tet mesh that make up the FSI interface
    nb = boundary_ids.size();
    for (unsigned ib = 0; ib < nb; ib++)
    {
      // Boundary number in "fluid" tet mesh
      unsigned b = boundary_ids[ib];


      // We'll setup boundary coordinates for this one
      Boundary_coordinate_exists[b] = true;

      // Remember for future reference
      FSI_boundary_id.push_back(b);

      // Loop over all elements on this boundary
      unsigned nel = tet_mesh_pt->nboundary_element(b);
      for (unsigned e = 0; e < nel; e++)
      {
        // Get pointer to the bulk fluid element that is adjacent to boundary b
        FiniteElement* bulk_elem_pt = tet_mesh_pt->boundary_element_pt(b, e);

        // Find the index of the face of element e along boundary b
        int face_index = tet_mesh_pt->face_index_at_boundary(b, e);

        // Create new face element
        FaceElement* face_el_pt = 0;
        if (tet_mesh_is_solid_mesh)
        {
#ifdef PARANOID
          if (dynamic_cast<SolidTElement<3, 3>*>(bulk_elem_pt) == 0)
          {
            std::ostringstream error_stream;
            error_stream
              << "Tet-element cannot be cast to SolidTElement<3,3>.\n"
              << "ThinBrickOnTetMesh can only be erected on mesh containing\n"
              << "quadratic tets." << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif
          face_el_pt =
            new DummyFaceElement<SolidTElement<3, 3>>(bulk_elem_pt, face_index);
        }
        else
        {
#ifdef PARANOID
          if (dynamic_cast<TElement<3, 3>*>(bulk_elem_pt) == 0)
          {
            std::ostringstream error_stream;
            error_stream
              << "Tet-element cannot be cast to TElement<3,3>.\n"
              << "ThinBrickOnTetMesh can only be erected on mesh containing\n"
              << "quadratic tets." << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif
          face_el_pt =
            new DummyFaceElement<TElement<3, 3>>(bulk_elem_pt, face_index);
        }


        // Specify boundary id in bulk mesh (needed to extract
        // boundary coordinate)
        face_el_pt->set_boundary_number_in_bulk_mesh(b);

        // Create storage for stack of brick elements
        Vector<Vector<FiniteElement*>> new_el_pt(3);

        // Sign of normal to detect inversion of FaceElement
        int normal_sign;

        // Loop over the three bricks that are built on the current
        //---------------------------------------------------------
        // triangular face
        //----------------
        for (unsigned j = 0; j < 3; j++)
        {
          // Build stack of bricks
          new_el_pt[j].resize(nlayer);
          for (unsigned ilayer = 0; ilayer < nlayer; ilayer++)
          {
            new_el_pt[j][ilayer] = new ELEMENT;
            Element_pt.push_back(new_el_pt[j][ilayer]);
          }

          Boundary_element_pt[b].push_back(new_el_pt[j][0]);
          Face_index_at_boundary[b].push_back(-3);

          // Associate zero-th node with vertex of triangular face
          //------------------------------------------------------
          unsigned j_local = 0;

          // Get normal sign
          normal_sign = face_el_pt->normal_sign();

          // Get coordinates etc of point from face: Vertex nodes enumerated
          // first....
          Vector<double> s = s_face[translate(normal_sign, j)];
          Vector<double> zeta(2);
          Vector<double> x(3);
          Vector<double> unit_normal(3);
          face_el_pt->interpolated_zeta(s, zeta);
          face_el_pt->interpolated_x(s, x);
          face_el_pt->outer_unit_normal(s, unit_normal);


          // Get node in the "fluid" mesh from face
          Node* fluid_node_pt = face_el_pt->node_pt(translate(normal_sign, j));

          // Has the corresponding "solid" node already been created?
          Node* existing_node_pt = solid_node_pt[fluid_node_pt];
          if (existing_node_pt == 0)
          {
            // Create new node
            Node* new_node_pt = new_el_pt[j][0]->construct_boundary_node(
              j_local, time_stepper_pt);
            Node_pt.push_back(new_node_pt);

            //...and remember it
            solid_node_pt[fluid_node_pt] = new_node_pt;

            // Set coordinates
            new_node_pt->x(0) = x[0];
            new_node_pt->x(1) = x[1];
            new_node_pt->x(2) = x[2];

            // Set boundary stuff -- boundary IDs copied from fluid
            bool only_on_fsi = true;
            std::set<unsigned>* bnd_pt;
            fluid_node_pt->get_boundaries_pt(bnd_pt);
            for (std::set<unsigned>::iterator it = (*bnd_pt).begin();
                 it != (*bnd_pt).end();
                 it++)
            {
              if (!is_on_fsi_boundary[(*it)]) only_on_fsi = false;
              add_boundary_node((*it), new_node_pt);
            }
            new_node_pt->set_coordinates_on_boundary(b, zeta);
            normals[new_node_pt].push_back(unit_normal);


            // If bottom node is only on FSI boundary, the nodes above
            // are not boundary nodes, apart from the last one!
            if (only_on_fsi)
            {
              // Create other nodes in bottom layer
              Node* new_nod_pt =
                new_el_pt[j][0]->construct_node(j_local + 9, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              // One layer thick?
              if (nlayer == 1)
              {
                Node* new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }
              else
              {
                Node* new_nod_pt = new_el_pt[j][0]->construct_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }

              // Now do other layers
              for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
              {
                // Copy bottom node from below
                new_el_pt[j][ilayer]->node_pt(j_local) =
                  connected_node_pt[new_node_pt][2 * ilayer - 1];

                // Create new nodes
                Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
                  j_local + 9, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);

                // Last node is boundary node
                if (ilayer != (nlayer - 1))
                {
                  Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
                    j_local + 18, time_stepper_pt);
                  connected_node_pt[new_node_pt].push_back(new_nod_pt);
                  Node_pt.push_back(new_nod_pt);
                }
                else
                {
                  Node* new_nod_pt =
                    new_el_pt[j][ilayer]->construct_boundary_node(
                      j_local + 18, time_stepper_pt);
                  connected_node_pt[new_node_pt].push_back(new_nod_pt);
                  Node_pt.push_back(new_nod_pt);
                }
              }
            }
            else
            {
              // Create other boundary nodes in bottom layer
              Node* new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                j_local + 9, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                j_local + 18, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              // Now do other layers
              for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
              {
                // Copy bottom node from below
                new_el_pt[j][ilayer]->node_pt(j_local) =
                  connected_node_pt[new_node_pt][2 * ilayer - 1];

                // Create new boundary nodes
                Node* new_nod_pt =
                  new_el_pt[j][ilayer]->construct_boundary_node(
                    j_local + 9, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
                new_nod_pt = new_el_pt[j][ilayer]->construct_boundary_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }
            }
          }
          else
          {
            // Add (repeated) bottom node to its other boundary and add
            // coordinates
            existing_node_pt->set_coordinates_on_boundary(b, zeta);
            normals[existing_node_pt].push_back(unit_normal);

            // Get pointer to nodes in bottom layer
            new_el_pt[j][0]->node_pt(j_local) = existing_node_pt;
            new_el_pt[j][0]->node_pt(j_local + 9) =
              connected_node_pt[existing_node_pt][0];
            new_el_pt[j][0]->node_pt(j_local + 18) =
              connected_node_pt[existing_node_pt][1];

            // Now do other layers
            for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
            {
              new_el_pt[j][ilayer]->node_pt(j_local) =
                connected_node_pt[existing_node_pt][2 * ilayer - 1];
              new_el_pt[j][ilayer]->node_pt(j_local + 9) =
                connected_node_pt[existing_node_pt][2 * ilayer];
              new_el_pt[j][ilayer]->node_pt(j_local + 18) =
                connected_node_pt[existing_node_pt][2 * ilayer + 1];
            }
          }


          // Second node with midside node in triangular face
          //-------------------------------------------------
          j_local = 2;

          // Get coordinates etc of point from face: Midside nodes enumerated
          // after vertex nodes
          s = s_face[translate(normal_sign, j + 3)];
          face_el_pt->interpolated_zeta(s, zeta);
          face_el_pt->interpolated_x(s, x);
          face_el_pt->outer_unit_normal(s, unit_normal);

          // Get node in the "fluid" mesh from face
          fluid_node_pt = face_el_pt->node_pt(translate(normal_sign, j + 3));

          // Has the corresponding "solid" node already been created?
          existing_node_pt = solid_node_pt[fluid_node_pt];
          if (existing_node_pt == 0)
          {
            // Create new node
            Node* new_node_pt = new_el_pt[j][0]->construct_boundary_node(
              j_local, time_stepper_pt);
            Node_pt.push_back(new_node_pt);

            // ...and remember it
            solid_node_pt[fluid_node_pt] = new_node_pt;

            // Set coordinates
            new_node_pt->x(0) = x[0];
            new_node_pt->x(1) = x[1];
            new_node_pt->x(2) = x[2];

            // Set boundary stuff -- boundary IDs copied from fluid
            bool only_on_fsi = true;
            std::set<unsigned>* bnd_pt;
            fluid_node_pt->get_boundaries_pt(bnd_pt);
            for (std::set<unsigned>::iterator it = (*bnd_pt).begin();
                 it != (*bnd_pt).end();
                 it++)
            {
              if (!is_on_fsi_boundary[(*it)]) only_on_fsi = false;
              add_boundary_node((*it), new_node_pt);
            }
            new_node_pt->set_coordinates_on_boundary(b, zeta);
            normals[new_node_pt].push_back(unit_normal);

            // If bottom node is only on FSI boundary, the nodes above
            // are not boundary nodes, apart from the last one!
            if (only_on_fsi)
            {
              // Create other nodes in bottom layer
              Node* new_nod_pt =
                new_el_pt[j][0]->construct_node(j_local + 9, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              // One layer thick?
              if (nlayer == 1)
              {
                Node* new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }
              else
              {
                Node* new_nod_pt = new_el_pt[j][0]->construct_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }

              // Now do other layers
              for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
              {
                // Copy bottom node from below
                new_el_pt[j][ilayer]->node_pt(j_local) =
                  connected_node_pt[new_node_pt][2 * ilayer - 1];

                // Create new nodes
                Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
                  j_local + 9, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);

                // Last node is boundary node
                if (ilayer != (nlayer - 1))
                {
                  Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
                    j_local + 18, time_stepper_pt);
                  connected_node_pt[new_node_pt].push_back(new_nod_pt);
                  Node_pt.push_back(new_nod_pt);
                }
                else
                {
                  Node* new_nod_pt =
                    new_el_pt[j][ilayer]->construct_boundary_node(
                      j_local + 18, time_stepper_pt);
                  connected_node_pt[new_node_pt].push_back(new_nod_pt);
                  Node_pt.push_back(new_nod_pt);
                }
              }
            }
            else
            {
              // Create other boundary nodes in bottom layer
              Node* new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                j_local + 9, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                j_local + 18, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              // Now do other layers
              for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
              {
                // Copy bottom node from below
                new_el_pt[j][ilayer]->node_pt(j_local) =
                  connected_node_pt[new_node_pt][2 * ilayer - 1];

                // Create new boundary nodes
                Node* new_nod_pt =
                  new_el_pt[j][ilayer]->construct_boundary_node(
                    j_local + 9, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
                new_nod_pt = new_el_pt[j][ilayer]->construct_boundary_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }
            }
          }
          else
          {
            // Add (repeated) bottom node to its other boundary and add
            // coordinates
            existing_node_pt->set_coordinates_on_boundary(b, zeta);
            normals[existing_node_pt].push_back(unit_normal);

            // Get pointer to nodes in bottom layer
            new_el_pt[j][0]->node_pt(j_local) = existing_node_pt;
            new_el_pt[j][0]->node_pt(j_local + 9) =
              connected_node_pt[existing_node_pt][0];
            new_el_pt[j][0]->node_pt(j_local + 18) =
              connected_node_pt[existing_node_pt][1];

            // Now do other layers
            for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
            {
              new_el_pt[j][ilayer]->node_pt(j_local) =
                connected_node_pt[existing_node_pt][2 * ilayer - 1];
              new_el_pt[j][ilayer]->node_pt(j_local + 9) =
                connected_node_pt[existing_node_pt][2 * ilayer];
              new_el_pt[j][ilayer]->node_pt(j_local + 18) =
                connected_node_pt[existing_node_pt][2 * ilayer + 1];
            }
          }


          // First node is quarter-edge node on triangular face
          //---------------------------------------------------
          j_local = 1;

          // Get coordinates of point from face: Quarter edge nodes enumerated
          // after midside nodes
          s = s_face[translate(normal_sign, 6 + 2 * j)];
          face_el_pt->interpolated_zeta(s, zeta);
          face_el_pt->interpolated_x(s, x);
          face_el_pt->outer_unit_normal(s, unit_normal);

          // Create Edge
          Edge edge(face_el_pt->node_pt(translate(normal_sign, j)),
                    face_el_pt->node_pt(translate(normal_sign, j + 3)));

          // Does node already exist?
          existing_node_pt = quarter_edge_node[edge];
          if (existing_node_pt == 0)
          {
            // Create new node
            Node* new_node_pt = new_el_pt[j][0]->construct_boundary_node(
              j_local, time_stepper_pt);
            Node_pt.push_back(new_node_pt);

            //...and remember it
            quarter_edge_node[edge] = new_node_pt;

            // Set coordinates
            new_node_pt->x(0) = x[0];
            new_node_pt->x(1) = x[1];
            new_node_pt->x(2) = x[2];

            // Set boundary stuff -- boundary IDs copied from fluid
            std::set<unsigned>* bnd1_pt;
            edge.node1_pt()->get_boundaries_pt(bnd1_pt);
            std::set<unsigned>* bnd2_pt;
            edge.node2_pt()->get_boundaries_pt(bnd2_pt);
            std::set<unsigned> bnd;
            set_intersection((*bnd1_pt).begin(),
                             (*bnd1_pt).end(),
                             (*bnd2_pt).begin(),
                             (*bnd2_pt).end(),
                             inserter(bnd, bnd.begin()));
            bool only_on_fsi = true;
            for (std::set<unsigned>::iterator it = bnd.begin(); it != bnd.end();
                 it++)
            {
              if (!is_on_fsi_boundary[(*it)]) only_on_fsi = false;
              add_boundary_node((*it), new_node_pt);
            }
            new_node_pt->set_coordinates_on_boundary(b, zeta);
            normals[new_node_pt].push_back(unit_normal);


            // If bottom node is only on FSI boundary, the nodes above
            // are not boundary nodes, apart from the last one!
            if (only_on_fsi)
            {
              // Create other nodes in bottom layer
              Node* new_nod_pt =
                new_el_pt[j][0]->construct_node(j_local + 9, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              // One layer thick?
              if (nlayer == 1)
              {
                Node* new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }
              else
              {
                Node* new_nod_pt = new_el_pt[j][0]->construct_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }

              // Now do other layers
              for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
              {
                // Copy bottom node from below
                new_el_pt[j][ilayer]->node_pt(j_local) =
                  connected_node_pt[new_node_pt][2 * ilayer - 1];

                // Create new nodes
                Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
                  j_local + 9, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);

                // Last node is boundary node
                if (ilayer != (nlayer - 1))
                {
                  Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
                    j_local + 18, time_stepper_pt);
                  connected_node_pt[new_node_pt].push_back(new_nod_pt);
                  Node_pt.push_back(new_nod_pt);
                }
                else
                {
                  Node* new_nod_pt =
                    new_el_pt[j][ilayer]->construct_boundary_node(
                      j_local + 18, time_stepper_pt);
                  connected_node_pt[new_node_pt].push_back(new_nod_pt);
                  Node_pt.push_back(new_nod_pt);
                }
              }
            }
            else
            {
              // Create other boundary nodes in bottom layer
              Node* new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                j_local + 9, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                j_local + 18, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              // Now do other layers
              for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
              {
                // Copy bottom node from below
                new_el_pt[j][ilayer]->node_pt(j_local) =
                  connected_node_pt[new_node_pt][2 * ilayer - 1];

                // Create new boundary nodes
                Node* new_nod_pt =
                  new_el_pt[j][ilayer]->construct_boundary_node(
                    j_local + 9, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
                new_nod_pt = new_el_pt[j][ilayer]->construct_boundary_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }
            }
          }
          else
          {
            // Add (repeated) bottom node to its other boundary and add
            // coordinates
            existing_node_pt->set_coordinates_on_boundary(b, zeta);
            normals[existing_node_pt].push_back(unit_normal);

            // Get pointer to nodes in bottom layer
            new_el_pt[j][0]->node_pt(j_local) = existing_node_pt;
            new_el_pt[j][0]->node_pt(j_local + 9) =
              connected_node_pt[existing_node_pt][0];
            new_el_pt[j][0]->node_pt(j_local + 18) =
              connected_node_pt[existing_node_pt][1];


            // Now do other layers
            for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
            {
              new_el_pt[j][ilayer]->node_pt(j_local) =
                connected_node_pt[existing_node_pt][2 * ilayer - 1];
              new_el_pt[j][ilayer]->node_pt(j_local + 9) =
                connected_node_pt[existing_node_pt][2 * ilayer];
              new_el_pt[j][ilayer]->node_pt(j_local + 18) =
                connected_node_pt[existing_node_pt][2 * ilayer + 1];
            }
          }


          // Third node is three-quarter-edge node on triangular face
          //---------------------------------------------------------
          j_local = 3;

          // Create Edge
          unsigned other_node = 0;
          unsigned jj = 0;
          switch (j)
          {
            case 0:
              other_node = 5;
              jj = 11;
              break;
            case 1:
              other_node = 3;
              jj = 7;
              break;
            case 2:
              other_node = 4;
              jj = 9;
              break;
          }
          Edge edge2(face_el_pt->node_pt(translate(normal_sign, j)),
                     face_el_pt->node_pt(translate(normal_sign, other_node)));

          // Get coordinates of point from face:
          s = s_face[translate(normal_sign, jj)];
          face_el_pt->interpolated_zeta(s, zeta);
          face_el_pt->interpolated_x(s, x);
          face_el_pt->outer_unit_normal(s, unit_normal);

          // Does node already exist?
          existing_node_pt = quarter_edge_node[edge2];
          if (existing_node_pt == 0)
          {
            // Create new node
            Node* new_node_pt = new_el_pt[j][0]->construct_boundary_node(
              j_local, time_stepper_pt);
            Node_pt.push_back(new_node_pt);

            //..and remember it
            quarter_edge_node[edge2] = new_node_pt;

            // Set coordinates
            new_node_pt->x(0) = x[0];
            new_node_pt->x(1) = x[1];
            new_node_pt->x(2) = x[2];

            // Set boundary stuff  -- boundary IDs copied from fluid
            std::set<unsigned>* bnd1_pt;
            edge2.node1_pt()->get_boundaries_pt(bnd1_pt);
            std::set<unsigned>* bnd2_pt;
            edge2.node2_pt()->get_boundaries_pt(bnd2_pt);
            std::set<unsigned> bnd;
            set_intersection((*bnd1_pt).begin(),
                             (*bnd1_pt).end(),
                             (*bnd2_pt).begin(),
                             (*bnd2_pt).end(),
                             inserter(bnd, bnd.begin()));
            bool only_on_fsi = true;
            for (std::set<unsigned>::iterator it = bnd.begin(); it != bnd.end();
                 it++)
            {
              if (!is_on_fsi_boundary[(*it)]) only_on_fsi = false;
              add_boundary_node((*it), new_node_pt);
            }
            new_node_pt->set_coordinates_on_boundary(b, zeta);
            normals[new_node_pt].push_back(unit_normal);

            // If bottom node is only on FSI boundary, the nodes above
            // are not boundary nodes, apart from the last one!
            if (only_on_fsi)
            {
              // Create other nodes in bottom layer
              Node* new_nod_pt =
                new_el_pt[j][0]->construct_node(j_local + 9, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              // One layer thick?
              if (nlayer == 1)
              {
                Node* new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }
              else
              {
                Node* new_nod_pt = new_el_pt[j][0]->construct_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }

              // Now do other layers
              for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
              {
                // Copy bottom node from below
                new_el_pt[j][ilayer]->node_pt(j_local) =
                  connected_node_pt[new_node_pt][2 * ilayer - 1];

                // Create new nodes
                Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
                  j_local + 9, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
                // Last node is boundary node
                if (ilayer != (nlayer - 1))
                {
                  Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
                    j_local + 18, time_stepper_pt);
                  connected_node_pt[new_node_pt].push_back(new_nod_pt);
                  Node_pt.push_back(new_nod_pt);
                }
                else
                {
                  Node* new_nod_pt =
                    new_el_pt[j][ilayer]->construct_boundary_node(
                      j_local + 18, time_stepper_pt);
                  connected_node_pt[new_node_pt].push_back(new_nod_pt);
                  Node_pt.push_back(new_nod_pt);
                }
              }
            }
            else
            {
              // Create other boundary nodes in bottom layer
              Node* new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                j_local + 9, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
                j_local + 18, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);

              // Now do other layers
              for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
              {
                // Copy bottom node from below
                new_el_pt[j][ilayer]->node_pt(j_local) =
                  connected_node_pt[new_node_pt][2 * ilayer - 1];

                // Create new boundary nodes
                Node* new_nod_pt =
                  new_el_pt[j][ilayer]->construct_boundary_node(
                    j_local + 9, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
                new_nod_pt = new_el_pt[j][ilayer]->construct_boundary_node(
                  j_local + 18, time_stepper_pt);
                connected_node_pt[new_node_pt].push_back(new_nod_pt);
                Node_pt.push_back(new_nod_pt);
              }
            }
          }
          else
          {
            // Add (repeated) bottom node to its other boundary and add
            // coordinates
            existing_node_pt->set_coordinates_on_boundary(b, zeta);
            normals[existing_node_pt].push_back(unit_normal);

            // Get pointer to nodes in bottom layer
            new_el_pt[j][0]->node_pt(j_local) = existing_node_pt;
            new_el_pt[j][0]->node_pt(j_local + 9) =
              connected_node_pt[existing_node_pt][0];
            new_el_pt[j][0]->node_pt(j_local + 18) =
              connected_node_pt[existing_node_pt][1];

            // Now do other layers
            for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
            {
              new_el_pt[j][ilayer]->node_pt(j_local) =
                connected_node_pt[existing_node_pt][2 * ilayer - 1];
              new_el_pt[j][ilayer]->node_pt(j_local + 9) =
                connected_node_pt[existing_node_pt][2 * ilayer];
              new_el_pt[j][ilayer]->node_pt(j_local + 18) =
                connected_node_pt[existing_node_pt][2 * ilayer + 1];
            }
          }


          // Fourth node is unique for all elements
          //--------------------------------------
          j_local = 4;

          // Create new node
          Node* new_node_pt =
            new_el_pt[j][0]->construct_boundary_node(j_local, time_stepper_pt);
          Node_pt.push_back(new_node_pt);

          jj = 0;
          switch (j)
          {
            case 0:
              jj = 15;
              break;
            case 1:
              jj = 17;
              break;
            case 2:
              jj = 13;
              break;
          }

          // Get coordinates etc of point from face:
          s = s_face[translate(normal_sign, jj)];
          face_el_pt->interpolated_zeta(s, zeta);
          face_el_pt->interpolated_x(s, x);
          face_el_pt->outer_unit_normal(s, unit_normal);

          // Set coordinates
          new_node_pt->x(0) = x[0];
          new_node_pt->x(1) = x[1];
          new_node_pt->x(2) = x[2];

          // Set boundary stuff
          add_boundary_node(b, new_node_pt);
          new_node_pt->set_coordinates_on_boundary(b, zeta);
          normals[new_node_pt].push_back(unit_normal);

          // Create other nodes in bottom layer
          Node* new_nod_pt =
            new_el_pt[j][0]->construct_node(j_local + 9, time_stepper_pt);
          connected_node_pt[new_node_pt].push_back(new_nod_pt);
          Node_pt.push_back(new_nod_pt);

          // One layer thick?
          if (nlayer == 1)
          {
            Node* new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
              j_local + 18, time_stepper_pt);
            connected_node_pt[new_node_pt].push_back(new_nod_pt);
            Node_pt.push_back(new_nod_pt);
          }
          else
          {
            Node* new_nod_pt =
              new_el_pt[j][0]->construct_node(j_local + 18, time_stepper_pt);
            connected_node_pt[new_node_pt].push_back(new_nod_pt);
            Node_pt.push_back(new_nod_pt);
          }

          // Now do other layers
          for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
          {
            // Copy bottom node from below
            new_el_pt[j][ilayer]->node_pt(j_local) =
              connected_node_pt[new_node_pt][2 * ilayer - 1];

            // Create new nodes
            Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
              j_local + 9, time_stepper_pt);
            connected_node_pt[new_node_pt].push_back(new_nod_pt);
            Node_pt.push_back(new_nod_pt);
            // Last node is boundary node
            if (ilayer != (nlayer - 1))
            {
              Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
                j_local + 18, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);
            }
            else
            {
              Node* new_nod_pt = new_el_pt[j][ilayer]->construct_boundary_node(
                j_local + 18, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);
            }
          }


          // Fifth node is created by all elements (internal to this
          //--------------------------------------------------------
          // this patch of elements to can't have been created elsewhere)
          //-------------------------------------------------------------

          j_local = 5;

          // Create new node
          new_node_pt =
            new_el_pt[j][0]->construct_boundary_node(j_local, time_stepper_pt);
          Node_pt.push_back(new_node_pt);

          // Get coordinates of point from face:
          jj = 0;
          switch (j)
          {
            case 0:
              jj = 14;
              break;
            case 1:
              jj = 16;
              break;
            case 2:
              jj = 18;
              break;
          }

          // Get coordinates etc from face
          s = s_face[translate(normal_sign, jj)];
          face_el_pt->interpolated_zeta(s, zeta);
          face_el_pt->interpolated_x(s, x);
          face_el_pt->outer_unit_normal(s, unit_normal);

          // Set coordinates
          new_node_pt->x(0) = x[0];
          new_node_pt->x(1) = x[1];
          new_node_pt->x(2) = x[2];

          // Set boundary stuff
          add_boundary_node(b, new_node_pt);
          new_node_pt->set_coordinates_on_boundary(b, zeta);
          normals[new_node_pt].push_back(unit_normal);

          // Create other nodes in bottom layer
          new_nod_pt =
            new_el_pt[j][0]->construct_node(j_local + 9, time_stepper_pt);
          connected_node_pt[new_node_pt].push_back(new_nod_pt);
          Node_pt.push_back(new_nod_pt);

          // One layer thick?
          if (nlayer == 1)
          {
            Node* new_nod_pt = new_el_pt[j][0]->construct_boundary_node(
              j_local + 18, time_stepper_pt);
            connected_node_pt[new_node_pt].push_back(new_nod_pt);
            Node_pt.push_back(new_nod_pt);
          }
          else
          {
            Node* new_nod_pt =
              new_el_pt[j][0]->construct_node(j_local + 18, time_stepper_pt);
            connected_node_pt[new_node_pt].push_back(new_nod_pt);
            Node_pt.push_back(new_nod_pt);
          }

          // Now do other layers
          for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
          {
            // Copy bottom node from below
            new_el_pt[j][ilayer]->node_pt(j_local) =
              connected_node_pt[new_node_pt][2 * ilayer - 1];

            // Create other nodes
            Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
              j_local + 9, time_stepper_pt);
            connected_node_pt[new_node_pt].push_back(new_nod_pt);
            Node_pt.push_back(new_nod_pt);

            // Last node is boundary node
            if (ilayer != (nlayer - 1))
            {
              Node* new_nod_pt = new_el_pt[j][ilayer]->construct_node(
                j_local + 18, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);
            }
            else
            {
              Node* new_nod_pt = new_el_pt[j][ilayer]->construct_boundary_node(
                j_local + 18, time_stepper_pt);
              connected_node_pt[new_node_pt].push_back(new_nod_pt);
              Node_pt.push_back(new_nod_pt);
            }
          }

        } // End over the three bricks erected on current triangular face


        // Last element builds central node as its node 8
        //-----------------------------------------------

        unsigned j_local = 8;

        // Create new node
        Node* new_node_pt =
          new_el_pt[2][0]->construct_boundary_node(j_local, time_stepper_pt);
        Node_pt.push_back(new_node_pt);

        // Get coordinates etc of point from face: Central node is
        // node 12 in face enumeration.
        Vector<double> s = s_face[12];
        Vector<double> zeta(2);
        Vector<double> x(3);
        Vector<double> unit_normal(3);
        face_el_pt->interpolated_zeta(s, zeta);
        face_el_pt->interpolated_x(s, x);
        face_el_pt->outer_unit_normal(s, unit_normal);

        // Set coordinates
        new_node_pt->x(0) = x[0];
        new_node_pt->x(1) = x[1];
        new_node_pt->x(2) = x[2];

        // Set boundary stuff
        add_boundary_node(b, new_node_pt);
        new_node_pt->set_coordinates_on_boundary(b, zeta);
        normals[new_node_pt].push_back(unit_normal);

        // Create other nodes in bottom layer
        Node* new_nod_pt =
          new_el_pt[2][0]->construct_node(j_local + 9, time_stepper_pt);
        connected_node_pt[new_node_pt].push_back(new_nod_pt);
        Node_pt.push_back(new_nod_pt);

        // One layer thick?
        if (nlayer == 1)
        {
          Node* new_nod_pt = new_el_pt[2][0]->construct_boundary_node(
            j_local + 18, time_stepper_pt);
          connected_node_pt[new_node_pt].push_back(new_nod_pt);
          Node_pt.push_back(new_nod_pt);
        }
        else
        {
          Node* new_nod_pt =
            new_el_pt[2][0]->construct_node(j_local + 18, time_stepper_pt);
          connected_node_pt[new_node_pt].push_back(new_nod_pt);
          Node_pt.push_back(new_nod_pt);
        }

        // Now do other layers
        for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
        {
          // Copy bottom node from below
          new_el_pt[2][ilayer]->node_pt(j_local) =
            connected_node_pt[new_node_pt][2 * ilayer - 1];

          // Create other nodes
          Node* new_nod_pt =
            new_el_pt[2][ilayer]->construct_node(j_local + 9, time_stepper_pt);
          connected_node_pt[new_node_pt].push_back(new_nod_pt);
          Node_pt.push_back(new_nod_pt);

          // Last node is boundary node
          if (ilayer != (nlayer - 1))
          {
            Node* new_nod_pt = new_el_pt[2][ilayer]->construct_node(
              j_local + 18, time_stepper_pt);
            connected_node_pt[new_node_pt].push_back(new_nod_pt);
            Node_pt.push_back(new_nod_pt);
          }
          else
          {
            Node* new_nod_pt = new_el_pt[2][ilayer]->construct_boundary_node(
              j_local + 18, time_stepper_pt);
            connected_node_pt[new_node_pt].push_back(new_nod_pt);
            Node_pt.push_back(new_nod_pt);
          }
        }

        // Other elements copy that node across
        new_el_pt[1][0]->node_pt(j_local) = new_node_pt;
        new_el_pt[0][0]->node_pt(j_local) = new_node_pt;

        new_el_pt[1][0]->node_pt(j_local + 9) =
          connected_node_pt[new_node_pt][0];
        new_el_pt[0][0]->node_pt(j_local + 9) =
          connected_node_pt[new_node_pt][0];

        new_el_pt[1][0]->node_pt(j_local + 18) =
          connected_node_pt[new_node_pt][1];
        new_el_pt[0][0]->node_pt(j_local + 18) =
          connected_node_pt[new_node_pt][1];

        // Now do layers
        for (unsigned ilayer = 1; ilayer < nlayer; ilayer++)
        {
          new_el_pt[1][ilayer]->node_pt(j_local) =
            connected_node_pt[new_node_pt][2 * ilayer - 1];
          new_el_pt[0][ilayer]->node_pt(j_local) =
            connected_node_pt[new_node_pt][2 * ilayer - 1];

          new_el_pt[1][ilayer]->node_pt(j_local + 9) =
            connected_node_pt[new_node_pt][2 * ilayer];
          new_el_pt[0][ilayer]->node_pt(j_local + 9) =
            connected_node_pt[new_node_pt][2 * ilayer];

          new_el_pt[1][ilayer]->node_pt(j_local + 18) =
            connected_node_pt[new_node_pt][2 * ilayer + 1];
          new_el_pt[0][ilayer]->node_pt(j_local + 18) =
            connected_node_pt[new_node_pt][2 * ilayer + 1];
        }


        // Nodes 6 and 7 in all elements are the same as nodes 2 and 5
        //------------------------------------------------------------
        // in previous element around the patch
        //-------------------------------------
        for (unsigned ilayer = 0; ilayer < nlayer; ilayer++)
        {
          for (unsigned j = 0; j < 3; j++)
          {
            unsigned offset = 9 * j;
            new_el_pt[2][ilayer]->node_pt(6 + offset) =
              new_el_pt[1][ilayer]->node_pt(2 + offset);

            new_el_pt[1][ilayer]->node_pt(6 + offset) =
              new_el_pt[0][ilayer]->node_pt(2 + offset);

            new_el_pt[0][ilayer]->node_pt(6 + offset) =
              new_el_pt[2][ilayer]->node_pt(2 + offset);

            new_el_pt[2][ilayer]->node_pt(7 + offset) =
              new_el_pt[1][ilayer]->node_pt(5 + offset);

            new_el_pt[1][ilayer]->node_pt(7 + offset) =
              new_el_pt[0][ilayer]->node_pt(5 + offset);

            new_el_pt[0][ilayer]->node_pt(7 + offset) =
              new_el_pt[2][ilayer]->node_pt(5 + offset);
          }
        }

        // Outer boundary is the last one
        Outer_boundary_id = maxb - 1;

        // Number of identified in/outflow domain boundaries
        // (remember they're broken up into separeate boundaries with oomph-lib)
        unsigned nb_in_out = is_on_in_out_boundary.size();
        In_out_boundary_id.resize(nb_in_out);

        // Now loop over the elements in the stacks again
        // and add all connected nodes to the appopriate non-FSI
        // boundary
        for (unsigned j_stack = 0; j_stack < 3; j_stack++)
        {
          // Bottom element
          FiniteElement* el_pt = new_el_pt[j_stack][0];

          // Loop over nodes in bottom layer
          for (unsigned j = 0; j < 9; j++)
          {
            // Get nodes above...
            Node* nod_pt = el_pt->node_pt(j);
            Vector<Node*> layer_node_pt = connected_node_pt[nod_pt];
            unsigned n = layer_node_pt.size();

            // Get boundary affiliation
            std::set<unsigned>* bnd_pt;
            nod_pt->get_boundaries_pt(bnd_pt);

            // Loop over boundaries
            for (std::set<unsigned>::iterator it = (*bnd_pt).begin();
                 it != (*bnd_pt).end();
                 it++)
            {
              // Ignore FSI surface!
              if (!is_on_fsi_boundary[(*it)])
              {
                // Loop over connnected nodes in layers above
                unsigned ilayer = 0;
                for (unsigned k = 0; k < n; k++)
                {
                  // Add to boundary
                  add_boundary_node((*it), layer_node_pt[k]);

                  int face_index = 0;

                  // Use edge node on bottom node layer to assess
                  // the element/s affiliation with boundary
                  if (j == 1) face_index = -2;
                  if (j == 3) face_index = -1;
                  if (j == 5) face_index = 1;
                  if (j == 7) face_index = 2;

                  if (face_index != 0)
                  {
                    // Use middle level in vertical direction
                    // to assess the element's affiliation with boundary
                    if (k % 2 == 1)
                    {
                      Boundary_element_pt[(*it)].push_back(
                        new_el_pt[j_stack][ilayer]);
                      Face_index_at_boundary[(*it)].push_back(face_index);
                      ilayer++;
                    }
                  }

                  // Add to lookup scheme that allows the nodes
                  // associated with an identified macroscopic in/outflow
                  // boundary to recovered collectively.

                  // Loop over macroscopic in/outflow boundaries
                  for (unsigned jj = 0; jj < nb_in_out; jj++)
                  {
                    if (is_on_in_out_boundary[jj][(*it)])
                    {
                      in_out_boundary_id_set[jj].insert((*it));
                    }
                  }
                }
              }
            }

            // Last connected node is on outer boundary
            add_boundary_node(Outer_boundary_id, layer_node_pt[n - 1]);

            // Use central node on bottom node layer to assess
            // the element/s affiliation with outer boundary
            if (j == 4)
            {
              Boundary_element_pt[Outer_boundary_id].push_back(
                new_el_pt[j_stack][nlayer - 1]);
              int face_index = 3;
              Face_index_at_boundary[Outer_boundary_id].push_back(face_index);
            }
          }
        }

        // Cleanup
        delete face_el_pt;
      }
    }


    // Copy boundary IDs across
    for (unsigned jj = 0; jj < n; jj++)
    {
      for (std::set<unsigned>::iterator it = in_out_boundary_id_set[jj].begin();
           it != in_out_boundary_id_set[jj].end();
           it++)
      {
        In_out_boundary_id[jj].push_back((*it));
      }
    }


#ifdef PARANOID
    // Check
    unsigned nel = Element_pt.size();
    for (unsigned e = 0; e < nel; e++)
    {
      FiniteElement* el_pt = finite_element_pt(e);
      for (unsigned j = 0; j < 27; j++)
      {
        if (el_pt->node_pt(j) == 0)
        {
          // Throw an error
          std::ostringstream error_stream;
          error_stream << "Null node in element " << e << " node " << j
                       << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
    }
#endif


    // Average unit normals
    std::ofstream outfile;
    bool doc_normals = false; // keep alive for future debugging
    if (doc_normals) outfile.open("normals.dat");
    for (std::map<Node*, Vector<Vector<double>>>::iterator it = normals.begin();
         it != normals.end();
         it++)
    {
      Vector<double> unit_normal(3, 0.0);
      unsigned nnormals = ((*it).second).size();
      for (unsigned j = 0; j < nnormals; j++)
      {
        for (unsigned i = 0; i < 3; i++)
        {
          unit_normal[i] += ((*it).second)[j][i];
        }
      }
      double norm = 0.0;
      for (unsigned i = 0; i < 3; i++)
      {
        norm += unit_normal[i] * unit_normal[i];
      }
      for (unsigned i = 0; i < 3; i++)
      {
        unit_normal[i] /= sqrt(norm);
      }

      Node* base_node_pt = (*it).first;
      Vector<double> base_pos(3);
      base_node_pt->position(base_pos);
      double h_thick;
      Thickness_fct_pt(base_pos, h_thick);
      Vector<Node*> layer_node_pt = connected_node_pt[base_node_pt];
      unsigned n = layer_node_pt.size();
      for (unsigned j = 0; j < n; j++)
      {
        for (unsigned i = 0; i < 3; i++)
        {
          layer_node_pt[j]->x(i) =
            base_pos[i] + h_thick * double(j + 1) / double(n) * unit_normal[i];
        }
      }
      if (doc_normals)
      {
        outfile << ((*it).first)->x(0) << " " << ((*it).first)->x(1) << " "
                << ((*it).first)->x(2) << " " << unit_normal[0] << " "
                << unit_normal[1] << " " << unit_normal[2] << "\n";
      }
    }
    if (doc_normals) outfile.close();


    // Lookup scheme has now been setup yet
    Lookup_for_elements_next_boundary_is_setup = true;
  }


} // namespace oomph

#endif
