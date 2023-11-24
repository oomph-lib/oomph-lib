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
#ifndef OOMPH_BRICK_FROM_TET_MESH_TEMPLATE_CC
#define OOMPH_BRICK_FROM_TET_MESH_TEMPLATE_CC


#include "brick_from_tet_mesh.template.h"

namespace oomph
{
  //=======================================================================
  /// Build fct: Pass pointer to existing tet mesh and timestepper
  /// Specialisation for XdaTetMesh<TElement<3,3> >
  //=======================================================================
  template<class ELEMENT>
  void BrickFromTetMesh<ELEMENT>::build_mesh(
    XdaTetMesh<TElement<3, 3>>* tet_mesh_pt, TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with 3D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3, 3);

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

    // Set number of boundaries
    unsigned nb = tet_mesh_pt->nboundary();
    set_nboundary(nb);

    // Get ready for boundary lookup scheme
    Boundary_element_pt.resize(nb);
    Face_index_at_boundary.resize(nb);

    // Maps to check which nodes have already been done

    // Map that stores the new brick node corresponding to an existing tet node
    std::map<Node*, Node*> tet_node_node_pt;

    // Map that stores node on an edge between two brick nodes
    std::map<Edge, Node*> brick_edge_node_pt;

    // Map that stores node on face spanned by three tet nodes
    std::map<TFace, Node*> tet_face_node_pt;

    // Create the four Dummy bricks:
    //------------------------------
    Vector<DummyBrickElement*> dummy_q_el_pt(4);
    for (unsigned e = 0; e < 4; e++)
    {
      dummy_q_el_pt[e] = new DummyBrickElement;
      for (unsigned j = 0; j < 8; j++)
      {
        dummy_q_el_pt[e]->construct_node(j);
      }
    }

    // Loop over the elements in the tet mesh
    unsigned n_el_tet = tet_mesh_pt->nelement();
    for (unsigned e_tet = 0; e_tet < n_el_tet; e_tet++)
    {
      // Cast to ten-noded tet
      TElement<3, 3>* tet_el_pt =
        dynamic_cast<TElement<3, 3>*>(tet_mesh_pt->element_pt(e_tet));

#ifdef PARANOID
      if (tet_el_pt == 0)
      {
        std::ostringstream error_stream;
        error_stream
          << "BrickFromTetMesh can only built from tet mesh containing\n"
          << "ten-noded tets.\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Storage for the centroid node for this tet
      Node* centroid_node_pt = 0;

      // Internal mid brick-face nodes
      Node* top_mid_face_node0_pt = 0;
      Node* right_mid_face_node0_pt = 0;
      Node* back_mid_face_node0_pt = 0;

      Node* top_mid_face_node1_pt = 0;
      Node* right_mid_face_node1_pt = 0;

      Node* top_mid_face_node2_pt = 0;

      // Newly created brick elements
      FiniteElement* brick_el0_pt = 0;
      FiniteElement* brick_el1_pt = 0;
      FiniteElement* brick_el2_pt = 0;
      FiniteElement* brick_el3_pt = 0;


      // First brick element is centred at node 0 of tet:
      //-------------------------------------------------
      {
        // Assign coordinates of dummy element
        for (unsigned j = 0; j < 8; j++)
        {
          Node* nod_pt = dummy_q_el_pt[0]->node_pt(j);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          switch (j)
          {
            case 0:
              tet_el_pt->local_coordinate_of_node(0, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 1:
              tet_el_pt->local_coordinate_of_node(4, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 2:
              tet_el_pt->local_coordinate_of_node(6, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 3:
              // label 13 in initial sketch: Mid face node on face spanned by
              // tet nodes 0,1,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 0.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 4:
              tet_el_pt->local_coordinate_of_node(5, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 5:
              // label 11 in initial sketch: Mid face node on face spanned
              // by tet nodes 0,1,2
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 6:
              // label 12 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,2,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 0.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 7:
              // label 14 in initial sketch: Centroid
              s_tet[0] = 0.25;
              s_tet[1] = 0.25;
              s_tet[2] = 0.25;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
          }
        }


        // Create actual zeroth brick element
        FiniteElement* el_pt = new ELEMENT;
        brick_el0_pt = el_pt;
        Element_pt.push_back(el_pt);

        TFace face0(
          tet_el_pt->node_pt(0), tet_el_pt->node_pt(1), tet_el_pt->node_pt(2));

        TFace face1(
          tet_el_pt->node_pt(0), tet_el_pt->node_pt(2), tet_el_pt->node_pt(3));

        TFace face2(
          tet_el_pt->node_pt(0), tet_el_pt->node_pt(1), tet_el_pt->node_pt(3));


        // Tet vertex nodes along edges emanating from node 0 in brick
        Vector<Vector<unsigned>> tet_edge_node(3);
        tet_edge_node[0].resize(2);
        tet_edge_node[0][0] = 4;
        tet_edge_node[0][1] = 1;
        tet_edge_node[1].resize(2);
        tet_edge_node[1][0] = 6;
        tet_edge_node[1][1] = 3;
        tet_edge_node[2].resize(2);
        tet_edge_node[2][0] = 5;
        tet_edge_node[2][1] = 2;

        // Node number of tet vertex that node 0 in brick is centred on
        unsigned central_tet_vertex = 0;

        Node* tet_node_pt = 0;
        Node* old_node_pt = 0;

        // Corner node
        {
          unsigned j = 0;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(central_tet_vertex);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-edge node on tet edge 0
        {
          unsigned j = 2;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[0][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid vertex node of tet edge 1
        {
          unsigned j = 6;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[1][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-vertex node of tet edge 2
        {
          unsigned j = 18;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[2][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in the middle of tet face0, spanned by
        // tet vertices 0, 1, 2. Enumerated "11" in initial sketch.
        {
          unsigned j = 20;

          // Need new node?
          old_node_pt = tet_face_node_pt[face0];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face0.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face0] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face1, spanned by
        // tet vertices 0, 2, 3. Enumerated "12" in initial sketch.
        {
          unsigned j = 24;

          // Need new node?
          old_node_pt = tet_face_node_pt[face1];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face1.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face1] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face2, spanned by
        // tet vertices 0, 1, 3. Enumerated "13" in initial sketch.
        {
          unsigned j = 8;

          // Need new node?
          old_node_pt = tet_face_node_pt[face2];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face2.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face2] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in centroid of tet. Only built for first element.
        // Enumerated "13" in initial sketch.
        {
          unsigned j = 26;

          // Always new
          {
            Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            centroid_node_pt = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
        }


        // Internal brick node -- always built
        {
          unsigned j = 13;

          // Always new
          {
            Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
        }

        // Brick edge node between brick nodes 0 and 2
        {
          unsigned j = 1;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(2));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 0 and 6
        {
          unsigned j = 3;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(6));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 2 and 8
        {
          unsigned j = 5;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 6 and 8
        {
          unsigned j = 7;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 18 and 20
        {
          unsigned j = 19;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 18 and 24
        {
          unsigned j = 21;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 20 and 26
        {
          unsigned j = 23;

          // Need new node?
          Edge edge(el_pt->node_pt(20), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 24 and 26
        {
          unsigned j = 25;

          // Need new node?
          Edge edge(el_pt->node_pt(24), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 0 and 18
        {
          unsigned j = 9;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(18));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 2 and 20
        {
          unsigned j = 11;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 6 and 24
        {
          unsigned j = 15;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 8 and 26
        {
          unsigned j = 17;

          // Need new node?
          Edge edge(el_pt->node_pt(8), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 2
        {
          unsigned j = 10;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 1 and 2
        {
          unsigned j = 12;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[1][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 1
        {
          unsigned j = 4;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[1][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Top mid brick-face node -- only built by first element
        {
          unsigned j = 22;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
          top_mid_face_node0_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }


        // Right mid brick-face node -- only built by first element
        {
          unsigned j = 14;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
          right_mid_face_node0_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }


        // Back mid brick-face node -- only built by first element
        {
          unsigned j = 16;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
          back_mid_face_node0_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }
      }


      // Second brick element is centred at node 1 of tet:
      //--------------------------------------------------
      {
        // Assign coordinates of dummy element
        for (unsigned j = 0; j < 8; j++)
        {
          Node* nod_pt = dummy_q_el_pt[1]->node_pt(j);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          switch (j)
          {
            case 0:
              tet_el_pt->local_coordinate_of_node(1, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 1:
              tet_el_pt->local_coordinate_of_node(9, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 2:
              tet_el_pt->local_coordinate_of_node(4, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 3:
              // label 13 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,1,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 0.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 4:
              tet_el_pt->local_coordinate_of_node(7, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 5:
              // label 10 in initial sketch: Mid face node on face
              // spanned by tet nodes 1,2,3
              s_tet[0] = 0.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 6:
              // label 11 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,1,2
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 7:
              // label 14 in initial sketch: Centroid
              s_tet[0] = 0.25;
              s_tet[1] = 0.25;
              s_tet[2] = 0.25;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
          }
        }


        // Create actual first brick element
        FiniteElement* el_pt = new ELEMENT;
        brick_el1_pt = el_pt;
        Element_pt.push_back(el_pt);

        TFace face0(
          tet_el_pt->node_pt(1), tet_el_pt->node_pt(3), tet_el_pt->node_pt(2));

        TFace face1(
          tet_el_pt->node_pt(1), tet_el_pt->node_pt(0), tet_el_pt->node_pt(2));

        TFace face2(
          tet_el_pt->node_pt(1), tet_el_pt->node_pt(0), tet_el_pt->node_pt(3));

        // Tet vertex nodes along edges emanating from node 0 in brick
        Vector<Vector<unsigned>> tet_edge_node(3);
        tet_edge_node[0].resize(2);
        tet_edge_node[0][0] = 9;
        tet_edge_node[0][1] = 3;
        tet_edge_node[1].resize(2);
        tet_edge_node[1][0] = 4;
        tet_edge_node[1][1] = 0;
        tet_edge_node[2].resize(2);
        tet_edge_node[2][0] = 7;
        tet_edge_node[2][1] = 2;

        // Node number of tet vertex that node 0 in brick is centred on
        unsigned central_tet_vertex = 1;

        Node* tet_node_pt = 0;
        Node* old_node_pt = 0;

        // Corner node
        {
          unsigned j = 0;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(central_tet_vertex);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-edge node on tet edge 0
        {
          unsigned j = 2;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[0][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid vertex node of tet edge 1
        {
          unsigned j = 6;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[1][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-vertex node of tet edge 2
        {
          unsigned j = 18;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[2][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in the middle of tet face0
        {
          unsigned j = 20;

          // Need new node?
          old_node_pt = tet_face_node_pt[face0];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face0.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face0] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face1
        {
          unsigned j = 24;

          // Need new node?
          old_node_pt = tet_face_node_pt[face1];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face1.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face1] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of face2
        {
          unsigned j = 8;

          // Need new node?
          old_node_pt = tet_face_node_pt[face2];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face2.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face2] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in centroid of tet. Only built for first element.
        // Enumerated "13" in initial sketch.
        {
          unsigned j = 26;

          // Always copied
          el_pt->node_pt(j) = centroid_node_pt;
        }


        // Internal brick node -- always built
        {
          unsigned j = 13;

          // Always new
          {
            Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
        }


        // Brick edge node between brick nodes 0 and 2
        {
          unsigned j = 1;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(2));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 0 and 6
        {
          unsigned j = 3;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(6));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 2 and 8
        {
          unsigned j = 5;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 6 and 8
        {
          unsigned j = 7;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 18 and 20
        {
          unsigned j = 19;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 18 and 24
        {
          unsigned j = 21;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 20 and 26
        {
          unsigned j = 23;

          // Need new node?
          Edge edge(el_pt->node_pt(20), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 24 and 26
        {
          unsigned j = 25;

          // Need new node?
          Edge edge(el_pt->node_pt(24), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 0 and 18
        {
          unsigned j = 9;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(18));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 2 and 20
        {
          unsigned j = 11;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 6 and 24
        {
          unsigned j = 15;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 8 and 26
        {
          unsigned j = 17;

          // Need new node?
          Edge edge(el_pt->node_pt(8), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 2
        {
          unsigned j = 10;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 1 and 2
        {
          unsigned j = 12;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[1][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 1
        {
          unsigned j = 4;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[1][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Top mid brick-face node -- only built by this element
        {
          unsigned j = 22;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
          top_mid_face_node1_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }


        // Right mid brick-face node -- only built by this element
        {
          unsigned j = 14;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
          right_mid_face_node1_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }


        // Back mid brick-face node copy from previous element
        {
          unsigned j = 16;

          // Always copied
          el_pt->node_pt(j) = right_mid_face_node0_pt;
        }
      }


      // Third brick element is centred at node 3 of tet:
      //-------------------------------------------------
      {
        // Assign coordinates of dummy element
        for (unsigned j = 0; j < 8; j++)
        {
          Node* nod_pt = dummy_q_el_pt[2]->node_pt(j);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          switch (j)
          {
            case 0:
              tet_el_pt->local_coordinate_of_node(3, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 1:
              tet_el_pt->local_coordinate_of_node(6, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 2:
              tet_el_pt->local_coordinate_of_node(9, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 3:
              // label 13 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,1,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 0.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 4:
              tet_el_pt->local_coordinate_of_node(8, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 5:
              // label 12 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,2,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 0.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 6:
              // label 10 in initial sketch: Mid face node on face
              // spanned by tet nodes 1,2,3
              s_tet[0] = 0.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 7:
              // label 14 in initial sketch: Centroid
              s_tet[0] = 0.25;
              s_tet[1] = 0.25;
              s_tet[2] = 0.25;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
          }
        }


        // Create actual second brick element
        FiniteElement* el_pt = new ELEMENT;
        brick_el2_pt = el_pt;
        Element_pt.push_back(el_pt);

        TFace face0(
          tet_el_pt->node_pt(3), tet_el_pt->node_pt(0), tet_el_pt->node_pt(2));

        TFace face1(
          tet_el_pt->node_pt(3), tet_el_pt->node_pt(1), tet_el_pt->node_pt(2));

        TFace face2(
          tet_el_pt->node_pt(3), tet_el_pt->node_pt(1), tet_el_pt->node_pt(0));

        // Tet vertex nodes along edges emanating from node 0 in brick
        Vector<Vector<unsigned>> tet_edge_node(3);
        tet_edge_node[0].resize(2);
        tet_edge_node[0][0] = 6;
        tet_edge_node[0][1] = 0;
        tet_edge_node[1].resize(2);
        tet_edge_node[1][0] = 9;
        tet_edge_node[1][1] = 1;
        tet_edge_node[2].resize(2);
        tet_edge_node[2][0] = 8;
        tet_edge_node[2][1] = 2;

        // Node number of tet vertex that node 0 in brick is centred on
        unsigned central_tet_vertex = 3;

        Node* tet_node_pt = 0;
        Node* old_node_pt = 0;

        // Corner node
        {
          unsigned j = 0;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(central_tet_vertex);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-edge node on tet edge 0
        {
          unsigned j = 2;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[0][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid vertex node of tet edge 1
        {
          unsigned j = 6;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[1][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-vertex node of tet edge 2
        {
          unsigned j = 18;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[2][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in the middle of tet face0
        {
          unsigned j = 20;

          // Need new node?
          old_node_pt = tet_face_node_pt[face0];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face0.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face0] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face1
        {
          unsigned j = 24;

          // Need new node?
          old_node_pt = tet_face_node_pt[face1];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face1.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face1] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face2
        {
          unsigned j = 8;

          // Need new node?
          old_node_pt = tet_face_node_pt[face2];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face2.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face2] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in centroid of tet. Only built for first element.
        // Enumerated "13" in initial sketch.
        {
          unsigned j = 26;

          // Always copied
          el_pt->node_pt(j) = centroid_node_pt;
        }


        // Internal brick node -- always built
        {
          unsigned j = 13;

          // Always new
          {
            Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
        }

        // Brick edge node between brick nodes 0 and 2
        {
          unsigned j = 1;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(2));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 0 and 6
        {
          unsigned j = 3;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(6));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 2 and 8
        {
          unsigned j = 5;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 6 and 8
        {
          unsigned j = 7;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 18 and 20
        {
          unsigned j = 19;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 18 and 24
        {
          unsigned j = 21;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 20 and 26
        {
          unsigned j = 23;

          // Need new node?
          Edge edge(el_pt->node_pt(20), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 24 and 26
        {
          unsigned j = 25;

          // Need new node?
          Edge edge(el_pt->node_pt(24), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 0 and 18
        {
          unsigned j = 9;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(18));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 2 and 20
        {
          unsigned j = 11;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 6 and 24
        {
          unsigned j = 15;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 8 and 26
        {
          unsigned j = 17;

          // Need new node?
          Edge edge(el_pt->node_pt(8), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 2
        {
          unsigned j = 10;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 1 and 2
        {
          unsigned j = 12;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[1][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));
          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 1
        {
          unsigned j = 4;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[1][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Top mid brick-face node -- only built by this element
        {
          unsigned j = 22;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
          top_mid_face_node2_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }


        // Right mid brick-face node copy from first element
        {
          unsigned j = 14;

          // Always copied
          el_pt->node_pt(j) = back_mid_face_node0_pt;
        }


        // Back mid brick-face node copy from previous element
        {
          unsigned j = 16;

          // Always copied
          el_pt->node_pt(j) = right_mid_face_node1_pt;
        }
      }


      // Fourth brick element is centred at node 2 of tet:
      //--------------------------------------------------
      {
        // Assign coordinates of dummy element
        for (unsigned j = 0; j < 8; j++)
        {
          Node* nod_pt = dummy_q_el_pt[3]->node_pt(j);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          switch (j)
          {
            case 0:
              tet_el_pt->local_coordinate_of_node(2, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 1:
              tet_el_pt->local_coordinate_of_node(7, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 2:
              tet_el_pt->local_coordinate_of_node(5, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 3:
              // label 11 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,1,2
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 4:
              tet_el_pt->local_coordinate_of_node(8, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 5:
              // label 10 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,2,3
              s_tet[0] = 0.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 6:
              // label 12 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,2,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 0.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 7:
              // label 14 in initial sketch: Centroid
              s_tet[0] = 0.25;
              s_tet[1] = 0.25;
              s_tet[2] = 0.25;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
          }
        }


        // Create actual third brick element
        FiniteElement* el_pt = new ELEMENT;
        brick_el3_pt = el_pt;
        Element_pt.push_back(el_pt);

        TFace face0(
          tet_el_pt->node_pt(1), tet_el_pt->node_pt(2), tet_el_pt->node_pt(3));

        TFace face1(
          tet_el_pt->node_pt(0), tet_el_pt->node_pt(2), tet_el_pt->node_pt(3));

        TFace face2(
          tet_el_pt->node_pt(0), tet_el_pt->node_pt(1), tet_el_pt->node_pt(2));

        // Tet vertex nodes along edges emanating from node 0 in brick
        Vector<Vector<unsigned>> tet_edge_node(3);
        tet_edge_node[0].resize(2);
        tet_edge_node[0][0] = 7;
        tet_edge_node[0][1] = 1;
        tet_edge_node[1].resize(2);
        tet_edge_node[1][0] = 5;
        tet_edge_node[1][1] = 0;
        tet_edge_node[2].resize(2);
        tet_edge_node[2][0] = 8;
        tet_edge_node[2][1] = 3;

        // Node number of tet vertex that node 0 in brick is centred on
        unsigned central_tet_vertex = 2;

        Node* tet_node_pt = 0;
        Node* old_node_pt = 0;

        // Corner node
        {
          unsigned j = 0;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(central_tet_vertex);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-edge node on tet edge 0
        {
          unsigned j = 2;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[0][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid vertex node of tet edge 1
        {
          unsigned j = 6;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[1][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-vertex node of tet edge 2
        {
          unsigned j = 18;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[2][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in the middle of tet face0
        {
          unsigned j = 20;

          // Need new node?
          old_node_pt = tet_face_node_pt[face0];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face0.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face0] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face1
        {
          unsigned j = 24;

          // Need new node?
          old_node_pt = tet_face_node_pt[face1];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face1.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face1] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face2
        {
          unsigned j = 8;

          // Need new node?
          old_node_pt = tet_face_node_pt[face2];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face2.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face2] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in centroid of tet. Only built for first element.
        // Enumerated "13" in initial sketch.
        {
          unsigned j = 26;

          // Always copied
          el_pt->node_pt(j) = centroid_node_pt;
        }


        // Internal brick node -- always built
        {
          unsigned j = 13;

          // Always new
          {
            Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
        }

        // Brick edge node between brick nodes 0 and 2
        {
          unsigned j = 1;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(2));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 0 and 6
        {
          unsigned j = 3;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(6));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 2 and 8
        {
          unsigned j = 5;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 6 and 8
        {
          unsigned j = 7;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 18 and 20
        {
          unsigned j = 19;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 18 and 24
        {
          unsigned j = 21;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 20 and 26
        {
          unsigned j = 23;

          // Need new node?
          Edge edge(el_pt->node_pt(20), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 24 and 26
        {
          unsigned j = 25;

          // Need new node?
          Edge edge(el_pt->node_pt(24), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 0 and 18
        {
          unsigned j = 9;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(18));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 2 and 20
        {
          unsigned j = 11;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 6 and 24
        {
          unsigned j = 15;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 8 and 26
        {
          unsigned j = 17;

          // Need new node?
          Edge edge(el_pt->node_pt(8), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 2
        {
          unsigned j = 10;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 1 and 2
        {
          unsigned j = 12;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[1][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 1
        {
          unsigned j = 4;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[1][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Top mid brick-face node copy from top of second element
        {
          unsigned j = 22;

          // Always copied
          el_pt->node_pt(j) = top_mid_face_node2_pt;
        }


        // Right mid brick-face node copy from top of first element
        {
          unsigned j = 14;

          // Always copied
          el_pt->node_pt(j) = top_mid_face_node1_pt;
        }


        // Back mid brick-face node  copy from top of zeroth element
        {
          unsigned j = 16;

          // Always copied
          el_pt->node_pt(j) = top_mid_face_node0_pt;
        }
      }


      // Check if the four faces of the tet are on a boundary.
      // If they are, add the nodes in the brick mesh to those
      // boundaries.

      // Loop over faces of tet (in its own face index enumeration)
      for (int face_index = 0; face_index < 4; face_index++)
      {
        // Identify face and coordinates in it
        TFace* face_pt = 0;
        switch (face_index)
        {
          case 0:
            // Face 0: s[0]=0; opposite node 0
            face_pt = new TFace(tet_el_pt->node_pt(1),
                                tet_el_pt->node_pt(2),
                                tet_el_pt->node_pt(3));
            break;

          case 1:
            // Face 1: s[1]=0; opposite node 1
            face_pt = new TFace(tet_el_pt->node_pt(0),
                                tet_el_pt->node_pt(2),
                                tet_el_pt->node_pt(3));

            break;


          case 2:
            // Face 2: s[2]=0; opposite node 2
            face_pt = new TFace(tet_el_pt->node_pt(0),
                                tet_el_pt->node_pt(1),
                                tet_el_pt->node_pt(3));

            break;

          case 3:
            // Face 3: s[0]+s[1]+s[2]=1; opposite node 3
            face_pt = new TFace(tet_el_pt->node_pt(0),
                                tet_el_pt->node_pt(1),
                                tet_el_pt->node_pt(2));
            break;
        }


        if (face_pt->is_boundary_face())
        {
          std::set<unsigned> bnd;
          std::set<unsigned>* bnd_pt = &bnd;
          face_pt->get_boundaries_pt(bnd_pt);

#ifdef PARANOID
          if ((*bnd_pt).size() > 1)
          {
            std::ostringstream error_stream;
            error_stream << "TFace should only be on one boundary.\n";
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

          if ((*bnd_pt).size() == 1)
          {
            // Create new face element
            FaceElement* face_el_pt = 0;
            if (tet_mesh_is_solid_mesh)
            {
#ifdef PARANOID
              if (dynamic_cast<SolidTElement<3, 3>*>(tet_el_pt) == 0)
              {
                std::ostringstream error_stream;
                error_stream
                  << "Tet-element cannot be cast to SolidTElement<3,3>.\n"
                  << "BrickFromTetMesh can only be built from\n"
                  << "mesh containing quadratic tets.\n"
                  << std::endl;
                throw OomphLibError(error_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif
              // Build the face element
              face_el_pt = new DummyFaceElement<SolidTElement<3, 3>>(
                tet_el_pt, face_index);
            }
            else
            {
#ifdef PARANOID
              if (dynamic_cast<TElement<3, 3>*>(tet_el_pt) == 0)
              {
                std::ostringstream error_stream;
                error_stream << "Tet-element cannot be cast to TElement<3,3>.\n"
                             << "BrickFromTetMesh can only be built from\n"
                             << "mesh containing quadratic tets.\n"
                             << std::endl;
                throw OomphLibError(error_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif

              // Build the face element
              face_el_pt =
                new DummyFaceElement<TElement<3, 3>>(tet_el_pt, face_index);
            }


            // Specify boundary id in bulk mesh (needed to extract
            // boundary coordinate)
            unsigned b = (*(*bnd_pt).begin());
            Boundary_coordinate_exists[b] = true;
            face_el_pt->set_boundary_number_in_bulk_mesh(b);

            // Now set up the brick nodes on this face, enumerated as
            // in s_face
            Vector<Node*> brick_face_node_pt(19);

            switch (face_index)
            {
              case 0:
                brick_face_node_pt[0] = brick_el1_pt->node_pt(0);
                brick_face_node_pt[1] = brick_el3_pt->node_pt(0);
                brick_face_node_pt[2] = brick_el2_pt->node_pt(0);

                brick_face_node_pt[3] = brick_el1_pt->node_pt(18);
                brick_face_node_pt[4] = brick_el2_pt->node_pt(18);
                brick_face_node_pt[5] = brick_el1_pt->node_pt(2);

                brick_face_node_pt[6] = brick_el1_pt->node_pt(9);
                brick_face_node_pt[7] = brick_el3_pt->node_pt(1);
                brick_face_node_pt[8] = brick_el3_pt->node_pt(9);

                brick_face_node_pt[9] = brick_el2_pt->node_pt(9);
                brick_face_node_pt[10] = brick_el2_pt->node_pt(3);
                brick_face_node_pt[11] = brick_el1_pt->node_pt(1);

                brick_face_node_pt[12] = brick_el1_pt->node_pt(20);

                brick_face_node_pt[13] = brick_el2_pt->node_pt(12);
                brick_face_node_pt[14] = brick_el1_pt->node_pt(19);

                brick_face_node_pt[15] = brick_el1_pt->node_pt(10);
                brick_face_node_pt[16] = brick_el2_pt->node_pt(21);

                brick_face_node_pt[17] = brick_el3_pt->node_pt(10);
                brick_face_node_pt[18] = brick_el1_pt->node_pt(11);
                break;

              case 1:
                brick_face_node_pt[0] = brick_el0_pt->node_pt(0);
                brick_face_node_pt[1] = brick_el3_pt->node_pt(0);
                brick_face_node_pt[2] = brick_el2_pt->node_pt(0);

                brick_face_node_pt[3] = brick_el0_pt->node_pt(18);
                brick_face_node_pt[4] = brick_el2_pt->node_pt(18);
                brick_face_node_pt[5] = brick_el0_pt->node_pt(6);

                brick_face_node_pt[6] = brick_el0_pt->node_pt(9);
                brick_face_node_pt[7] = brick_el3_pt->node_pt(3);
                brick_face_node_pt[8] = brick_el3_pt->node_pt(9);

                brick_face_node_pt[9] = brick_el2_pt->node_pt(9);
                brick_face_node_pt[10] = brick_el2_pt->node_pt(1);
                brick_face_node_pt[11] = brick_el0_pt->node_pt(3);

                brick_face_node_pt[12] = brick_el0_pt->node_pt(24);

                brick_face_node_pt[13] = brick_el2_pt->node_pt(10);
                brick_face_node_pt[14] = brick_el0_pt->node_pt(21);

                brick_face_node_pt[15] = brick_el0_pt->node_pt(12);
                brick_face_node_pt[16] = brick_el3_pt->node_pt(21);

                brick_face_node_pt[17] = brick_el3_pt->node_pt(12);
                brick_face_node_pt[18] = brick_el0_pt->node_pt(15);
                break;

              case 2:
                brick_face_node_pt[0] = brick_el0_pt->node_pt(0);
                brick_face_node_pt[1] = brick_el1_pt->node_pt(0);
                brick_face_node_pt[2] = brick_el2_pt->node_pt(0);

                brick_face_node_pt[3] = brick_el0_pt->node_pt(2);
                brick_face_node_pt[4] = brick_el1_pt->node_pt(2);
                brick_face_node_pt[5] = brick_el0_pt->node_pt(6);

                brick_face_node_pt[6] = brick_el0_pt->node_pt(1);
                brick_face_node_pt[7] = brick_el1_pt->node_pt(3);
                brick_face_node_pt[8] = brick_el1_pt->node_pt(1);

                brick_face_node_pt[9] = brick_el2_pt->node_pt(3);
                brick_face_node_pt[10] = brick_el2_pt->node_pt(1);
                brick_face_node_pt[11] = brick_el0_pt->node_pt(3);

                brick_face_node_pt[12] = brick_el0_pt->node_pt(8);

                brick_face_node_pt[13] = brick_el2_pt->node_pt(4);
                brick_face_node_pt[14] = brick_el0_pt->node_pt(5);

                brick_face_node_pt[15] = brick_el0_pt->node_pt(4);
                brick_face_node_pt[16] = brick_el1_pt->node_pt(5);

                brick_face_node_pt[17] = brick_el1_pt->node_pt(4);
                brick_face_node_pt[18] = brick_el0_pt->node_pt(7);
                break;

              case 3:
                brick_face_node_pt[0] = brick_el1_pt->node_pt(0);
                brick_face_node_pt[1] = brick_el3_pt->node_pt(0);
                brick_face_node_pt[2] = brick_el0_pt->node_pt(0);

                brick_face_node_pt[3] = brick_el1_pt->node_pt(18);
                brick_face_node_pt[4] = brick_el3_pt->node_pt(6);
                brick_face_node_pt[5] = brick_el1_pt->node_pt(6);

                brick_face_node_pt[6] = brick_el1_pt->node_pt(9);
                brick_face_node_pt[7] = brick_el3_pt->node_pt(1);
                brick_face_node_pt[8] = brick_el3_pt->node_pt(3);

                brick_face_node_pt[9] = brick_el0_pt->node_pt(9);
                brick_face_node_pt[10] = brick_el0_pt->node_pt(1);
                brick_face_node_pt[11] = brick_el1_pt->node_pt(3);

                brick_face_node_pt[12] = brick_el1_pt->node_pt(24);

                brick_face_node_pt[13] = brick_el0_pt->node_pt(10);
                brick_face_node_pt[14] = brick_el1_pt->node_pt(21);

                brick_face_node_pt[15] = brick_el1_pt->node_pt(12);
                brick_face_node_pt[16] = brick_el3_pt->node_pt(7);

                brick_face_node_pt[17] = brick_el3_pt->node_pt(4);
                brick_face_node_pt[18] = brick_el1_pt->node_pt(15);
                break;
            }

            // Provide possibility for translation -- may need to add
            // face index to this; see ThinLayerBrickOnTetMesh.
            Vector<unsigned> translate(19);

            // Initialise with identity mapping
            for (unsigned i = 0; i < 19; i++)
            {
              translate[i] = i;
            }

            // Visit all the nodes on that face
            for (unsigned j = 0; j < 19; j++)
            {
              // Which node is it?
              Node* brick_node_pt = brick_face_node_pt[translate[j]];

              // Get coordinates etc of point from face
              Vector<double> s = s_face[j];
              Vector<double> zeta(2);
              Vector<double> x(3);
              face_el_pt->interpolated_zeta(s, zeta);
              face_el_pt->interpolated_x(s, x);

#ifdef PARANOID
              // Check that the coordinates match (within tolerance)
              double dist = sqrt(pow(brick_node_pt->x(0) - x[0], 2) +
                                 pow(brick_node_pt->x(1) - x[1], 2) +
                                 pow(brick_node_pt->x(2) - x[2], 2));
              if (dist > BrickFromTetMeshHelper::Face_position_tolerance)
              {
                std::ofstream brick0;
                std::ofstream brick1;
                std::ofstream brick2;
                std::ofstream brick3;
                brick0.open("full_brick0.dat");
                brick1.open("full_brick1.dat");
                brick2.open("full_brick2.dat");
                brick3.open("full_brick3.dat");
                for (unsigned j = 0; j < 27; j++)
                {
                  brick0 << brick_el0_pt->node_pt(j)->x(0) << " "
                         << brick_el0_pt->node_pt(j)->x(1) << " "
                         << brick_el0_pt->node_pt(j)->x(2) << "\n";

                  brick1 << brick_el1_pt->node_pt(j)->x(0) << " "
                         << brick_el1_pt->node_pt(j)->x(1) << " "
                         << brick_el1_pt->node_pt(j)->x(2) << "\n";

                  brick2 << brick_el2_pt->node_pt(j)->x(0) << " "
                         << brick_el2_pt->node_pt(j)->x(1) << " "
                         << brick_el2_pt->node_pt(j)->x(2) << "\n";

                  brick3 << brick_el3_pt->node_pt(j)->x(0) << " "
                         << brick_el3_pt->node_pt(j)->x(1) << " "
                         << brick_el3_pt->node_pt(j)->x(2) << "\n";
                }
                brick0.close();
                brick1.close();
                brick2.close();
                brick3.close();

                std::ofstream full_face;
                full_face.open("full_face.dat");
                for (unsigned j = 0; j < 6; j++)
                {
                  full_face << face_el_pt->node_pt(j)->x(0) << " "
                            << face_el_pt->node_pt(j)->x(1) << " "
                            << face_el_pt->node_pt(j)->x(2) << "\n";
                }
                full_face.close();

                // Get normal sign
                int normal_sign = face_el_pt->normal_sign();

                std::ostringstream error_stream;
                error_stream
                  << "During assignment of boundary cordinates, the distance\n"
                  << "between brick node and reference point in \n"
                  << "triangular FaceElement is " << dist << " which \n"
                  << "is bigger than the tolerance defined in \n"
                  << "BrickFromTetMeshHelper::Face_position_tolerance="
                  << BrickFromTetMeshHelper::Face_position_tolerance << ".\n"
                  << "If this is tolerable, increase the tolerance \n"
                  << "(it's defined in a namespace and therefore publically\n"
                  << "accessible). If not, the Face may be inverted in which \n"
                  << "case you should re-implement the translation scheme,\n"
                  << "following the pattern used in the "
                     "ThinLayerBrickOnTetMesh."
                  << "\nThe required code fragements are already here but \n"
                  << "the translation is the unit map.\n"
                  << "To aid the diagnostics, the files full_brick[0-3].dat\n"
                  << "contain the coordinates of the 27 nodes in the four\n"
                  << "bricks associated with the current tet and "
                     "full_face.dat\n"
                  << "contains the coordinates of the 6 nodes in the "
                     "FaceElement"
                  << "\nfrom which the boundary coordinates are extracted.\n"
                  << "FYI: The normal_sign of the face is: " << normal_sign
                  << std::endl;
                throw OomphLibError(error_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif

              // Set boundary stuff
              add_boundary_node(b, brick_node_pt);
              brick_node_pt->set_coordinates_on_boundary(b, zeta);
            }

            // Add appropriate brick elements to boundary lookup scheme
            switch (face_index)
            {
              case 0:
                Boundary_element_pt[b].push_back(brick_el1_pt);
                Face_index_at_boundary[b].push_back(-2);
                Boundary_element_pt[b].push_back(brick_el2_pt);
                Face_index_at_boundary[b].push_back(-1);
                Boundary_element_pt[b].push_back(brick_el3_pt);
                Face_index_at_boundary[b].push_back(-2);
                break;

              case 1:
                Boundary_element_pt[b].push_back(brick_el0_pt);
                Face_index_at_boundary[b].push_back(-1);
                Boundary_element_pt[b].push_back(brick_el2_pt);
                Face_index_at_boundary[b].push_back(-2);
                Boundary_element_pt[b].push_back(brick_el3_pt);
                Face_index_at_boundary[b].push_back(-1);
                break;

              case 2:
                Boundary_element_pt[b].push_back(brick_el0_pt);
                Face_index_at_boundary[b].push_back(-3);
                Boundary_element_pt[b].push_back(brick_el1_pt);
                Face_index_at_boundary[b].push_back(-3);
                Boundary_element_pt[b].push_back(brick_el2_pt);
                Face_index_at_boundary[b].push_back(-3);
                break;

              case 3:
                Boundary_element_pt[b].push_back(brick_el0_pt);
                Face_index_at_boundary[b].push_back(-2);
                Boundary_element_pt[b].push_back(brick_el1_pt);
                Face_index_at_boundary[b].push_back(-1);
                Boundary_element_pt[b].push_back(brick_el3_pt);
                Face_index_at_boundary[b].push_back(-3);
                break;
            }
            // Cleanup
            delete face_el_pt;
          }
        }
        // Cleanup
        delete face_pt;
      }
    }

    // Lookup scheme has now been setup
    Lookup_for_elements_next_boundary_is_setup = true;

    // Get number of distinct boundaries specified
    // in the original xda enumeration.
    unsigned n_xda_boundaries = tet_mesh_pt->nxda_boundary();

    // Copy collective IDs across
    Boundary_id.resize(n_xda_boundaries);
    for (unsigned xda_b = 0; xda_b < n_xda_boundaries; xda_b++)
    {
      Boundary_id[xda_b] = tet_mesh_pt->oomph_lib_boundary_ids(xda_b);
    }


    // Cleanup
    for (unsigned e = 0; e < 4; e++)
    {
      for (unsigned j = 0; j < 8; j++)
      {
        delete dummy_q_el_pt[e]->node_pt(j);
      }
      delete dummy_q_el_pt[e];
    }
  }

  //=======================================================================
  /// Build fct: Pass pointer to existing tet mesh and timestepper
  /// Specialisation for TetgenMesh<TElement<3,3> >
  //=======================================================================
  template<class ELEMENT>
  void BrickFromTetMesh<ELEMENT>::build_mesh(
    TetgenMesh<TElement<3, 3>>* tet_mesh_pt, TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with 3D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3, 3);

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

    // Set number of boundaries
    unsigned nb = tet_mesh_pt->nboundary();
    set_nboundary(nb);

    // Get ready for boundary lookup scheme
    Boundary_element_pt.resize(nb);
    Face_index_at_boundary.resize(nb);

    // Maps to check which nodes have already been done

    // Map that stores the new brick node corresponding to an existing tet node
    std::map<Node*, Node*> tet_node_node_pt;

    // Map that stores node on an edge between two brick nodes
    std::map<Edge, Node*> brick_edge_node_pt;

    // Map that stores node on face spanned by three tet nodes
    std::map<TFace, Node*> tet_face_node_pt;

    // Create the four Dummy bricks:
    //------------------------------
    Vector<DummyBrickElement*> dummy_q_el_pt(4);
    for (unsigned e = 0; e < 4; e++)
    {
      dummy_q_el_pt[e] = new DummyBrickElement;
      for (unsigned j = 0; j < 8; j++)
      {
        dummy_q_el_pt[e]->construct_node(j);
      }
    }

    // Loop over the elements in the tet mesh
    unsigned n_el_tet = tet_mesh_pt->nelement();
    for (unsigned e_tet = 0; e_tet < n_el_tet; e_tet++)
    {
      // Cast to ten-noded tet
      TElement<3, 3>* tet_el_pt =
        dynamic_cast<TElement<3, 3>*>(tet_mesh_pt->element_pt(e_tet));

#ifdef PARANOID
      if (tet_el_pt == 0)
      {
        std::ostringstream error_stream;
        error_stream
          << "BrickFromTetMesh can only built from tet mesh containing\n"
          << "ten-noded tets.\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Storage for the centroid node for this tet
      Node* centroid_node_pt = 0;

      // Internal mid brick-face nodes
      Node* top_mid_face_node0_pt = 0;
      Node* right_mid_face_node0_pt = 0;
      Node* back_mid_face_node0_pt = 0;

      Node* top_mid_face_node1_pt = 0;
      Node* right_mid_face_node1_pt = 0;

      Node* top_mid_face_node2_pt = 0;

      // Newly created brick elements
      FiniteElement* brick_el0_pt = 0;
      FiniteElement* brick_el1_pt = 0;
      FiniteElement* brick_el2_pt = 0;
      FiniteElement* brick_el3_pt = 0;


      // First brick element is centred at node 0 of tet:
      //-------------------------------------------------
      {
        // Assign coordinates of dummy element
        for (unsigned j = 0; j < 8; j++)
        {
          Node* nod_pt = dummy_q_el_pt[0]->node_pt(j);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          switch (j)
          {
            case 0:
              tet_el_pt->local_coordinate_of_node(0, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 1:
              tet_el_pt->local_coordinate_of_node(4, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 2:
              tet_el_pt->local_coordinate_of_node(6, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 3:
              // label 13 in initial sketch: Mid face node on face spanned by
              // tet nodes 0,1,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 0.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 4:
              tet_el_pt->local_coordinate_of_node(5, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 5:
              // label 11 in initial sketch: Mid face node on face spanned
              // by tet nodes 0,1,2
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 6:
              // label 12 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,2,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 0.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 7:
              // label 14 in initial sketch: Centroid
              s_tet[0] = 0.25;
              s_tet[1] = 0.25;
              s_tet[2] = 0.25;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
          }
        }


        // Create actual zeroth brick element
        FiniteElement* el_pt = new ELEMENT;
        brick_el0_pt = el_pt;
        Element_pt.push_back(el_pt);

        TFace face0(
          tet_el_pt->node_pt(0), tet_el_pt->node_pt(1), tet_el_pt->node_pt(2));

        TFace face1(
          tet_el_pt->node_pt(0), tet_el_pt->node_pt(2), tet_el_pt->node_pt(3));

        TFace face2(
          tet_el_pt->node_pt(0), tet_el_pt->node_pt(1), tet_el_pt->node_pt(3));


        // Tet vertex nodes along edges emanating from node 0 in brick
        Vector<Vector<unsigned>> tet_edge_node(3);
        tet_edge_node[0].resize(2);
        tet_edge_node[0][0] = 4;
        tet_edge_node[0][1] = 1;
        tet_edge_node[1].resize(2);
        tet_edge_node[1][0] = 6;
        tet_edge_node[1][1] = 3;
        tet_edge_node[2].resize(2);
        tet_edge_node[2][0] = 5;
        tet_edge_node[2][1] = 2;

        // Node number of tet vertex that node 0 in brick is centred on
        unsigned central_tet_vertex = 0;

        Node* tet_node_pt = 0;
        Node* old_node_pt = 0;

        // Corner node
        {
          unsigned j = 0;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(central_tet_vertex);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-edge node on tet edge 0
        {
          unsigned j = 2;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[0][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid vertex node of tet edge 1
        {
          unsigned j = 6;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[1][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-vertex node of tet edge 2
        {
          unsigned j = 18;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[2][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in the middle of tet face0, spanned by
        // tet vertices 0, 1, 2. Enumerated "11" in initial sketch.
        {
          unsigned j = 20;

          // Need new node?
          old_node_pt = tet_face_node_pt[face0];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face0.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face0] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face1, spanned by
        // tet vertices 0, 2, 3. Enumerated "12" in initial sketch.
        {
          unsigned j = 24;

          // Need new node?
          old_node_pt = tet_face_node_pt[face1];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face1.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face1] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face2, spanned by
        // tet vertices 0, 1, 3. Enumerated "13" in initial sketch.
        {
          unsigned j = 8;

          // Need new node?
          old_node_pt = tet_face_node_pt[face2];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face2.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face2] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in centroid of tet. Only built for first element.
        // Enumerated "13" in initial sketch.
        {
          unsigned j = 26;

          // Always new
          {
            Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            centroid_node_pt = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
        }


        // Internal brick node -- always built
        {
          unsigned j = 13;

          // Always new
          {
            Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
        }

        // Brick edge node between brick nodes 0 and 2
        {
          unsigned j = 1;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(2));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 0 and 6
        {
          unsigned j = 3;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(6));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 2 and 8
        {
          unsigned j = 5;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 6 and 8
        {
          unsigned j = 7;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 18 and 20
        {
          unsigned j = 19;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 18 and 24
        {
          unsigned j = 21;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 20 and 26
        {
          unsigned j = 23;

          // Need new node?
          Edge edge(el_pt->node_pt(20), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 24 and 26
        {
          unsigned j = 25;

          // Need new node?
          Edge edge(el_pt->node_pt(24), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 0 and 18
        {
          unsigned j = 9;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(18));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 2 and 20
        {
          unsigned j = 11;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 6 and 24
        {
          unsigned j = 15;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 8 and 26
        {
          unsigned j = 17;

          // Need new node?
          Edge edge(el_pt->node_pt(8), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 2
        {
          unsigned j = 10;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 1 and 2
        {
          unsigned j = 12;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[1][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 1
        {
          unsigned j = 4;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[1][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Top mid brick-face node -- only built by first element
        {
          unsigned j = 22;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
          top_mid_face_node0_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }


        // Right mid brick-face node -- only built by first element
        {
          unsigned j = 14;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
          right_mid_face_node0_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }


        // Back mid brick-face node -- only built by first element
        {
          unsigned j = 16;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[0]->interpolated_s_tet(s, s_tet);
          back_mid_face_node0_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }
      }


      // Second brick element is centred at node 1 of tet:
      //--------------------------------------------------
      {
        // Assign coordinates of dummy element
        for (unsigned j = 0; j < 8; j++)
        {
          Node* nod_pt = dummy_q_el_pt[1]->node_pt(j);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          switch (j)
          {
            case 0:
              tet_el_pt->local_coordinate_of_node(1, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 1:
              tet_el_pt->local_coordinate_of_node(9, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 2:
              tet_el_pt->local_coordinate_of_node(4, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 3:
              // label 13 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,1,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 0.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 4:
              tet_el_pt->local_coordinate_of_node(7, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 5:
              // label 10 in initial sketch: Mid face node on face
              // spanned by tet nodes 1,2,3
              s_tet[0] = 0.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 6:
              // label 11 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,1,2
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 7:
              // label 14 in initial sketch: Centroid
              s_tet[0] = 0.25;
              s_tet[1] = 0.25;
              s_tet[2] = 0.25;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
          }
        }


        // Create actual first brick element
        FiniteElement* el_pt = new ELEMENT;
        brick_el1_pt = el_pt;
        Element_pt.push_back(el_pt);

        TFace face0(
          tet_el_pt->node_pt(1), tet_el_pt->node_pt(3), tet_el_pt->node_pt(2));

        TFace face1(
          tet_el_pt->node_pt(1), tet_el_pt->node_pt(0), tet_el_pt->node_pt(2));

        TFace face2(
          tet_el_pt->node_pt(1), tet_el_pt->node_pt(0), tet_el_pt->node_pt(3));

        // Tet vertex nodes along edges emanating from node 0 in brick
        Vector<Vector<unsigned>> tet_edge_node(3);
        tet_edge_node[0].resize(2);
        tet_edge_node[0][0] = 9;
        tet_edge_node[0][1] = 3;
        tet_edge_node[1].resize(2);
        tet_edge_node[1][0] = 4;
        tet_edge_node[1][1] = 0;
        tet_edge_node[2].resize(2);
        tet_edge_node[2][0] = 7;
        tet_edge_node[2][1] = 2;

        // Node number of tet vertex that node 0 in brick is centred on
        unsigned central_tet_vertex = 1;

        Node* tet_node_pt = 0;
        Node* old_node_pt = 0;

        // Corner node
        {
          unsigned j = 0;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(central_tet_vertex);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-edge node on tet edge 0
        {
          unsigned j = 2;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[0][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid vertex node of tet edge 1
        {
          unsigned j = 6;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[1][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-vertex node of tet edge 2
        {
          unsigned j = 18;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[2][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in the middle of tet face0
        {
          unsigned j = 20;

          // Need new node?
          old_node_pt = tet_face_node_pt[face0];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face0.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face0] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face1
        {
          unsigned j = 24;

          // Need new node?
          old_node_pt = tet_face_node_pt[face1];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face1.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face1] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of face2
        {
          unsigned j = 8;

          // Need new node?
          old_node_pt = tet_face_node_pt[face2];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face2.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face2] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in centroid of tet. Only built for first element.
        // Enumerated "13" in initial sketch.
        {
          unsigned j = 26;

          // Always copied
          el_pt->node_pt(j) = centroid_node_pt;
        }


        // Internal brick node -- always built
        {
          unsigned j = 13;

          // Always new
          {
            Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
        }


        // Brick edge node between brick nodes 0 and 2
        {
          unsigned j = 1;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(2));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 0 and 6
        {
          unsigned j = 3;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(6));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 2 and 8
        {
          unsigned j = 5;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 6 and 8
        {
          unsigned j = 7;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 18 and 20
        {
          unsigned j = 19;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 18 and 24
        {
          unsigned j = 21;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 20 and 26
        {
          unsigned j = 23;

          // Need new node?
          Edge edge(el_pt->node_pt(20), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 24 and 26
        {
          unsigned j = 25;

          // Need new node?
          Edge edge(el_pt->node_pt(24), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 0 and 18
        {
          unsigned j = 9;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(18));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 2 and 20
        {
          unsigned j = 11;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 6 and 24
        {
          unsigned j = 15;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 8 and 26
        {
          unsigned j = 17;

          // Need new node?
          Edge edge(el_pt->node_pt(8), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 2
        {
          unsigned j = 10;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 1 and 2
        {
          unsigned j = 12;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[1][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 1
        {
          unsigned j = 4;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[1][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Top mid brick-face node -- only built by this element
        {
          unsigned j = 22;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
          top_mid_face_node1_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }


        // Right mid brick-face node -- only built by this element
        {
          unsigned j = 14;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[1]->interpolated_s_tet(s, s_tet);
          right_mid_face_node1_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }


        // Back mid brick-face node copy from previous element
        {
          unsigned j = 16;

          // Always copied
          el_pt->node_pt(j) = right_mid_face_node0_pt;
        }
      }


      // Third brick element is centred at node 3 of tet:
      //-------------------------------------------------
      {
        // Assign coordinates of dummy element
        for (unsigned j = 0; j < 8; j++)
        {
          Node* nod_pt = dummy_q_el_pt[2]->node_pt(j);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          switch (j)
          {
            case 0:
              tet_el_pt->local_coordinate_of_node(3, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 1:
              tet_el_pt->local_coordinate_of_node(6, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 2:
              tet_el_pt->local_coordinate_of_node(9, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 3:
              // label 13 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,1,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 0.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 4:
              tet_el_pt->local_coordinate_of_node(8, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 5:
              // label 12 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,2,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 0.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 6:
              // label 10 in initial sketch: Mid face node on face
              // spanned by tet nodes 1,2,3
              s_tet[0] = 0.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 7:
              // label 14 in initial sketch: Centroid
              s_tet[0] = 0.25;
              s_tet[1] = 0.25;
              s_tet[2] = 0.25;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
          }
        }


        // Create actual second brick element
        FiniteElement* el_pt = new ELEMENT;
        brick_el2_pt = el_pt;
        Element_pt.push_back(el_pt);

        TFace face0(
          tet_el_pt->node_pt(3), tet_el_pt->node_pt(0), tet_el_pt->node_pt(2));

        TFace face1(
          tet_el_pt->node_pt(3), tet_el_pt->node_pt(1), tet_el_pt->node_pt(2));

        TFace face2(
          tet_el_pt->node_pt(3), tet_el_pt->node_pt(1), tet_el_pt->node_pt(0));

        // Tet vertex nodes along edges emanating from node 0 in brick
        Vector<Vector<unsigned>> tet_edge_node(3);
        tet_edge_node[0].resize(2);
        tet_edge_node[0][0] = 6;
        tet_edge_node[0][1] = 0;
        tet_edge_node[1].resize(2);
        tet_edge_node[1][0] = 9;
        tet_edge_node[1][1] = 1;
        tet_edge_node[2].resize(2);
        tet_edge_node[2][0] = 8;
        tet_edge_node[2][1] = 2;

        // Node number of tet vertex that node 0 in brick is centred on
        unsigned central_tet_vertex = 3;

        Node* tet_node_pt = 0;
        Node* old_node_pt = 0;

        // Corner node
        {
          unsigned j = 0;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(central_tet_vertex);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-edge node on tet edge 0
        {
          unsigned j = 2;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[0][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid vertex node of tet edge 1
        {
          unsigned j = 6;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[1][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-vertex node of tet edge 2
        {
          unsigned j = 18;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[2][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in the middle of tet face0
        {
          unsigned j = 20;

          // Need new node?
          old_node_pt = tet_face_node_pt[face0];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face0.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face0] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face1
        {
          unsigned j = 24;

          // Need new node?
          old_node_pt = tet_face_node_pt[face1];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face1.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face1] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face2
        {
          unsigned j = 8;

          // Need new node?
          old_node_pt = tet_face_node_pt[face2];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face2.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face2] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in centroid of tet. Only built for first element.
        // Enumerated "13" in initial sketch.
        {
          unsigned j = 26;

          // Always copied
          el_pt->node_pt(j) = centroid_node_pt;
        }


        // Internal brick node -- always built
        {
          unsigned j = 13;

          // Always new
          {
            Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
        }

        // Brick edge node between brick nodes 0 and 2
        {
          unsigned j = 1;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(2));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 0 and 6
        {
          unsigned j = 3;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(6));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 2 and 8
        {
          unsigned j = 5;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 6 and 8
        {
          unsigned j = 7;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 18 and 20
        {
          unsigned j = 19;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 18 and 24
        {
          unsigned j = 21;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 20 and 26
        {
          unsigned j = 23;

          // Need new node?
          Edge edge(el_pt->node_pt(20), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 24 and 26
        {
          unsigned j = 25;

          // Need new node?
          Edge edge(el_pt->node_pt(24), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 0 and 18
        {
          unsigned j = 9;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(18));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 2 and 20
        {
          unsigned j = 11;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 6 and 24
        {
          unsigned j = 15;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 8 and 26
        {
          unsigned j = 17;

          // Need new node?
          Edge edge(el_pt->node_pt(8), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 2
        {
          unsigned j = 10;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 1 and 2
        {
          unsigned j = 12;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[1][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));
          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 1
        {
          unsigned j = 4;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[1][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Top mid brick-face node -- only built by this element
        {
          unsigned j = 22;
          Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
          Node_pt.push_back(new_node_pt);
          Vector<double> s(3);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          el_pt->local_coordinate_of_node(j, s);
          dummy_q_el_pt[2]->interpolated_s_tet(s, s_tet);
          top_mid_face_node2_pt = new_node_pt;
          tet_el_pt->interpolated_x(s_tet, x_tet);
          new_node_pt->x(0) = x_tet[0];
          new_node_pt->x(1) = x_tet[1];
          new_node_pt->x(2) = x_tet[2];
        }


        // Right mid brick-face node copy from first element
        {
          unsigned j = 14;

          // Always copied
          el_pt->node_pt(j) = back_mid_face_node0_pt;
        }


        // Back mid brick-face node copy from previous element
        {
          unsigned j = 16;

          // Always copied
          el_pt->node_pt(j) = right_mid_face_node1_pt;
        }
      }


      // Fourth brick element is centred at node 2 of tet:
      //--------------------------------------------------
      {
        // Assign coordinates of dummy element
        for (unsigned j = 0; j < 8; j++)
        {
          Node* nod_pt = dummy_q_el_pt[3]->node_pt(j);
          Vector<double> s_tet(3);
          Vector<double> x_tet(3);
          switch (j)
          {
            case 0:
              tet_el_pt->local_coordinate_of_node(2, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 1:
              tet_el_pt->local_coordinate_of_node(7, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 2:
              tet_el_pt->local_coordinate_of_node(5, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 3:
              // label 11 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,1,2
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 4:
              tet_el_pt->local_coordinate_of_node(8, s_tet);
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 5:
              // label 10 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,2,3
              s_tet[0] = 0.0;
              s_tet[1] = 1.0 / 3.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 6:
              // label 12 in initial sketch: Mid face node on face
              // spanned by tet nodes 0,2,3
              s_tet[0] = 1.0 / 3.0;
              s_tet[1] = 0.0;
              s_tet[2] = 1.0 / 3.0;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
            case 7:
              // label 14 in initial sketch: Centroid
              s_tet[0] = 0.25;
              s_tet[1] = 0.25;
              s_tet[2] = 0.25;
              nod_pt->set_value(0, s_tet[0]);
              nod_pt->set_value(1, s_tet[1]);
              nod_pt->set_value(2, s_tet[2]);
              tet_el_pt->interpolated_x(s_tet, x_tet);
              nod_pt->x(0) = x_tet[0];
              nod_pt->x(1) = x_tet[1];
              nod_pt->x(2) = x_tet[2];
              break;
          }
        }


        // Create actual third brick element
        FiniteElement* el_pt = new ELEMENT;
        brick_el3_pt = el_pt;
        Element_pt.push_back(el_pt);

        TFace face0(
          tet_el_pt->node_pt(1), tet_el_pt->node_pt(2), tet_el_pt->node_pt(3));

        TFace face1(
          tet_el_pt->node_pt(0), tet_el_pt->node_pt(2), tet_el_pt->node_pt(3));

        TFace face2(
          tet_el_pt->node_pt(0), tet_el_pt->node_pt(1), tet_el_pt->node_pt(2));

        // Tet vertex nodes along edges emanating from node 0 in brick
        Vector<Vector<unsigned>> tet_edge_node(3);
        tet_edge_node[0].resize(2);
        tet_edge_node[0][0] = 7;
        tet_edge_node[0][1] = 1;
        tet_edge_node[1].resize(2);
        tet_edge_node[1][0] = 5;
        tet_edge_node[1][1] = 0;
        tet_edge_node[2].resize(2);
        tet_edge_node[2][0] = 8;
        tet_edge_node[2][1] = 3;

        // Node number of tet vertex that node 0 in brick is centred on
        unsigned central_tet_vertex = 2;

        Node* tet_node_pt = 0;
        Node* old_node_pt = 0;

        // Corner node
        {
          unsigned j = 0;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(central_tet_vertex);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-edge node on tet edge 0
        {
          unsigned j = 2;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[0][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid vertex node of tet edge 1
        {
          unsigned j = 6;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[1][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node coindides with mid-vertex node of tet edge 2
        {
          unsigned j = 18;

          // Need new node?
          tet_node_pt = tet_el_pt->node_pt(tet_edge_node[2][0]);
          old_node_pt = tet_node_node_pt[tet_node_pt];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (tet_node_pt->is_on_boundary())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_node_node_pt[tet_node_pt] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in the middle of tet face0
        {
          unsigned j = 20;

          // Need new node?
          old_node_pt = tet_face_node_pt[face0];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face0.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face0] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face1
        {
          unsigned j = 24;

          // Need new node?
          old_node_pt = tet_face_node_pt[face1];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face1.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face1] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick vertex node in the middle of tet face2
        {
          unsigned j = 8;

          // Need new node?
          old_node_pt = tet_face_node_pt[face2];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face2.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face2] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick vertex node in centroid of tet. Only built for first element.
        // Enumerated "13" in initial sketch.
        {
          unsigned j = 26;

          // Always copied
          el_pt->node_pt(j) = centroid_node_pt;
        }


        // Internal brick node -- always built
        {
          unsigned j = 13;

          // Always new
          {
            Node* new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
        }

        // Brick edge node between brick nodes 0 and 2
        {
          unsigned j = 1;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(2));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 0 and 6
        {
          unsigned j = 3;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(6));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 2 and 8
        {
          unsigned j = 5;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 6 and 8
        {
          unsigned j = 7;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(8));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 18 and 20
        {
          unsigned j = 19;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 18 and 24
        {
          unsigned j = 21;

          // Need new node?
          Edge edge(el_pt->node_pt(18), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 20 and 26
        {
          unsigned j = 23;

          // Need new node?
          Edge edge(el_pt->node_pt(20), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 24 and 26
        {
          unsigned j = 25;

          // Need new node?
          Edge edge(el_pt->node_pt(24), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }

        // Brick edge node between brick nodes 0 and 18
        {
          unsigned j = 9;

          // Need new node?
          Edge edge(el_pt->node_pt(0), el_pt->node_pt(18));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 2 and 20
        {
          unsigned j = 11;

          // Need new node?
          Edge edge(el_pt->node_pt(2), el_pt->node_pt(20));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 6 and 24
        {
          unsigned j = 15;

          // Need new node?
          Edge edge(el_pt->node_pt(6), el_pt->node_pt(24));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Brick edge node between brick nodes 8 and 26
        {
          unsigned j = 17;

          // Need new node?
          Edge edge(el_pt->node_pt(8), el_pt->node_pt(26));
          old_node_pt = brick_edge_node_pt[edge];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (edge.is_boundary_edge())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            brick_edge_node_pt[edge] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 2
        {
          unsigned j = 10;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 1 and 2
        {
          unsigned j = 12;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[1][0]),
                     tet_el_pt->node_pt(tet_edge_node[2][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Mid brick-face node associated with face
        // spanned by mid-vertex nodes associated with tet edges 0 and 1
        {
          unsigned j = 4;

          // Need new node?
          TFace face(tet_el_pt->node_pt(central_tet_vertex),
                     tet_el_pt->node_pt(tet_edge_node[0][0]),
                     tet_el_pt->node_pt(tet_edge_node[1][0]));

          old_node_pt = tet_face_node_pt[face];
          if (old_node_pt == 0)
          {
            Node* new_node_pt = 0;
            if (face.is_boundary_face())
            {
              new_node_pt = el_pt->construct_boundary_node(j, time_stepper_pt);
            }
            else
            {
              new_node_pt = el_pt->construct_node(j, time_stepper_pt);
            }
            tet_face_node_pt[face] = new_node_pt;
            Node_pt.push_back(new_node_pt);
            Vector<double> s(3);
            Vector<double> s_tet(3);
            Vector<double> x_tet(3);
            el_pt->local_coordinate_of_node(j, s);
            dummy_q_el_pt[3]->interpolated_s_tet(s, s_tet);
            tet_el_pt->interpolated_x(s_tet, x_tet);
            new_node_pt->x(0) = x_tet[0];
            new_node_pt->x(1) = x_tet[1];
            new_node_pt->x(2) = x_tet[2];
          }
          // Node already exists
          else
          {
            el_pt->node_pt(j) = old_node_pt;
          }
        }


        // Top mid brick-face node copy from top of second element
        {
          unsigned j = 22;

          // Always copied
          el_pt->node_pt(j) = top_mid_face_node2_pt;
        }


        // Right mid brick-face node copy from top of first element
        {
          unsigned j = 14;

          // Always copied
          el_pt->node_pt(j) = top_mid_face_node1_pt;
        }


        // Back mid brick-face node  copy from top of zeroth element
        {
          unsigned j = 16;

          // Always copied
          el_pt->node_pt(j) = top_mid_face_node0_pt;
        }
      }


      // Check if the four faces of the tet are on a boundary.
      // If they are, add the nodes in the brick mesh to those
      // boundaries.

      // Loop over faces of tet (in its own face index enumeration)
      for (int face_index = 0; face_index < 4; face_index++)
      {
        // Identify face and coordinates in it
        TFace* face_pt = 0;
        switch (face_index)
        {
          case 0:
            // Face 0: s[0]=0; opposite node 0
            face_pt = new TFace(tet_el_pt->node_pt(1),
                                tet_el_pt->node_pt(2),
                                tet_el_pt->node_pt(3));
            break;

          case 1:
            // Face 1: s[1]=0; opposite node 1
            face_pt = new TFace(tet_el_pt->node_pt(0),
                                tet_el_pt->node_pt(2),
                                tet_el_pt->node_pt(3));

            break;


          case 2:
            // Face 2: s[2]=0; opposite node 2
            face_pt = new TFace(tet_el_pt->node_pt(0),
                                tet_el_pt->node_pt(1),
                                tet_el_pt->node_pt(3));

            break;

          case 3:
            // Face 3: s[0]+s[1]+s[2]=1; opposite node 3
            face_pt = new TFace(tet_el_pt->node_pt(0),
                                tet_el_pt->node_pt(1),
                                tet_el_pt->node_pt(2));
            break;
        }


        if (face_pt->is_boundary_face())
        {
          std::set<unsigned> bnd;
          std::set<unsigned>* bnd_pt = &bnd;
          face_pt->get_boundaries_pt(bnd_pt);

#ifdef PARANOID
          if ((*bnd_pt).size() > 1)
          {
            std::ostringstream error_stream;
            error_stream << "TFace should only be on one boundary.\n";
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

          if ((*bnd_pt).size() == 1)
          {
            // Create new face element
            FaceElement* face_el_pt = 0;
            if (tet_mesh_is_solid_mesh)
            {
#ifdef PARANOID
              if (dynamic_cast<SolidTElement<3, 3>*>(tet_el_pt) == 0)
              {
                std::ostringstream error_stream;
                error_stream
                  << "Tet-element cannot be cast to SolidTElement<3,3>.\n"
                  << "BrickFromTetMesh can only be built from\n"
                  << "mesh containing quadratic tets.\n"
                  << std::endl;
                throw OomphLibError(error_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif
              // Build the face element
              face_el_pt = new DummyFaceElement<SolidTElement<3, 3>>(
                tet_el_pt, face_index);
            }
            else
            {
#ifdef PARANOID
              if (dynamic_cast<TElement<3, 3>*>(tet_el_pt) == 0)
              {
                std::ostringstream error_stream;
                error_stream << "Tet-element cannot be cast to TElement<3,3>.\n"
                             << "BrickFromTetMesh can only be built from\n"
                             << "mesh containing quadratic tets.\n"
                             << std::endl;
                throw OomphLibError(error_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif

              // Build the face element
              face_el_pt =
                new DummyFaceElement<TElement<3, 3>>(tet_el_pt, face_index);
            }


            // Specify boundary id in bulk mesh (needed to extract
            // boundary coordinate)
            unsigned b = (*(*bnd_pt).begin());
            Boundary_coordinate_exists[b] = true;
            face_el_pt->set_boundary_number_in_bulk_mesh(b);

            // Now set up the brick nodes on this face, enumerated as
            // in s_face
            Vector<Node*> brick_face_node_pt(19);

            switch (face_index)
            {
              case 0:
                brick_face_node_pt[0] = brick_el1_pt->node_pt(0);
                brick_face_node_pt[1] = brick_el3_pt->node_pt(0);
                brick_face_node_pt[2] = brick_el2_pt->node_pt(0);

                brick_face_node_pt[3] = brick_el1_pt->node_pt(18);
                brick_face_node_pt[4] = brick_el2_pt->node_pt(18);
                brick_face_node_pt[5] = brick_el1_pt->node_pt(2);

                brick_face_node_pt[6] = brick_el1_pt->node_pt(9);
                brick_face_node_pt[7] = brick_el3_pt->node_pt(1);
                brick_face_node_pt[8] = brick_el3_pt->node_pt(9);

                brick_face_node_pt[9] = brick_el2_pt->node_pt(9);
                brick_face_node_pt[10] = brick_el2_pt->node_pt(3);
                brick_face_node_pt[11] = brick_el1_pt->node_pt(1);

                brick_face_node_pt[12] = brick_el1_pt->node_pt(20);

                brick_face_node_pt[13] = brick_el2_pt->node_pt(12);
                brick_face_node_pt[14] = brick_el1_pt->node_pt(19);

                brick_face_node_pt[15] = brick_el1_pt->node_pt(10);
                brick_face_node_pt[16] = brick_el2_pt->node_pt(21);

                brick_face_node_pt[17] = brick_el3_pt->node_pt(10);
                brick_face_node_pt[18] = brick_el1_pt->node_pt(11);
                break;

              case 1:
                brick_face_node_pt[0] = brick_el0_pt->node_pt(0);
                brick_face_node_pt[1] = brick_el3_pt->node_pt(0);
                brick_face_node_pt[2] = brick_el2_pt->node_pt(0);

                brick_face_node_pt[3] = brick_el0_pt->node_pt(18);
                brick_face_node_pt[4] = brick_el2_pt->node_pt(18);
                brick_face_node_pt[5] = brick_el0_pt->node_pt(6);

                brick_face_node_pt[6] = brick_el0_pt->node_pt(9);
                brick_face_node_pt[7] = brick_el3_pt->node_pt(3);
                brick_face_node_pt[8] = brick_el3_pt->node_pt(9);

                brick_face_node_pt[9] = brick_el2_pt->node_pt(9);
                brick_face_node_pt[10] = brick_el2_pt->node_pt(1);
                brick_face_node_pt[11] = brick_el0_pt->node_pt(3);

                brick_face_node_pt[12] = brick_el0_pt->node_pt(24);

                brick_face_node_pt[13] = brick_el2_pt->node_pt(10);
                brick_face_node_pt[14] = brick_el0_pt->node_pt(21);

                brick_face_node_pt[15] = brick_el0_pt->node_pt(12);
                brick_face_node_pt[16] = brick_el3_pt->node_pt(21);

                brick_face_node_pt[17] = brick_el3_pt->node_pt(12);
                brick_face_node_pt[18] = brick_el0_pt->node_pt(15);
                break;

              case 2:
                brick_face_node_pt[0] = brick_el0_pt->node_pt(0);
                brick_face_node_pt[1] = brick_el1_pt->node_pt(0);
                brick_face_node_pt[2] = brick_el2_pt->node_pt(0);

                brick_face_node_pt[3] = brick_el0_pt->node_pt(2);
                brick_face_node_pt[4] = brick_el1_pt->node_pt(2);
                brick_face_node_pt[5] = brick_el0_pt->node_pt(6);

                brick_face_node_pt[6] = brick_el0_pt->node_pt(1);
                brick_face_node_pt[7] = brick_el1_pt->node_pt(3);
                brick_face_node_pt[8] = brick_el1_pt->node_pt(1);

                brick_face_node_pt[9] = brick_el2_pt->node_pt(3);
                brick_face_node_pt[10] = brick_el2_pt->node_pt(1);
                brick_face_node_pt[11] = brick_el0_pt->node_pt(3);

                brick_face_node_pt[12] = brick_el0_pt->node_pt(8);

                brick_face_node_pt[13] = brick_el2_pt->node_pt(4);
                brick_face_node_pt[14] = brick_el0_pt->node_pt(5);

                brick_face_node_pt[15] = brick_el0_pt->node_pt(4);
                brick_face_node_pt[16] = brick_el1_pt->node_pt(5);

                brick_face_node_pt[17] = brick_el1_pt->node_pt(4);
                brick_face_node_pt[18] = brick_el0_pt->node_pt(7);
                break;

              case 3:
                brick_face_node_pt[0] = brick_el1_pt->node_pt(0);
                brick_face_node_pt[1] = brick_el3_pt->node_pt(0);
                brick_face_node_pt[2] = brick_el0_pt->node_pt(0);

                brick_face_node_pt[3] = brick_el1_pt->node_pt(18);
                brick_face_node_pt[4] = brick_el3_pt->node_pt(6);
                brick_face_node_pt[5] = brick_el1_pt->node_pt(6);

                brick_face_node_pt[6] = brick_el1_pt->node_pt(9);
                brick_face_node_pt[7] = brick_el3_pt->node_pt(1);
                brick_face_node_pt[8] = brick_el3_pt->node_pt(3);

                brick_face_node_pt[9] = brick_el0_pt->node_pt(9);
                brick_face_node_pt[10] = brick_el0_pt->node_pt(1);
                brick_face_node_pt[11] = brick_el1_pt->node_pt(3);

                brick_face_node_pt[12] = brick_el1_pt->node_pt(24);

                brick_face_node_pt[13] = brick_el0_pt->node_pt(10);
                brick_face_node_pt[14] = brick_el1_pt->node_pt(21);

                brick_face_node_pt[15] = brick_el1_pt->node_pt(12);
                brick_face_node_pt[16] = brick_el3_pt->node_pt(7);

                brick_face_node_pt[17] = brick_el3_pt->node_pt(4);
                brick_face_node_pt[18] = brick_el1_pt->node_pt(15);
                break;
            }

            // Provide possibility for translation -- may need to add
            // face index to this; see ThinLayerBrickOnTetMesh.
            Vector<unsigned> translate(19);

            // Initialise with identity mapping
            for (unsigned i = 0; i < 19; i++)
            {
              translate[i] = i;
            }

            // Visit all the nodes on that face
            for (unsigned j = 0; j < 19; j++)
            {
              // Which node is it?
              Node* brick_node_pt = brick_face_node_pt[translate[j]];

              // Get coordinates etc of point from face
              Vector<double> s = s_face[j];
              Vector<double> zeta(2);
              Vector<double> x(3);
              face_el_pt->interpolated_zeta(s, zeta);
              face_el_pt->interpolated_x(s, x);

#ifdef PARANOID
              // Check that the coordinates match (within tolerance)
              double dist = sqrt(pow(brick_node_pt->x(0) - x[0], 2) +
                                 pow(brick_node_pt->x(1) - x[1], 2) +
                                 pow(brick_node_pt->x(2) - x[2], 2));
              if (dist > BrickFromTetMeshHelper::Face_position_tolerance)
              {
                std::ofstream brick0;
                std::ofstream brick1;
                std::ofstream brick2;
                std::ofstream brick3;
                brick0.open("full_brick0.dat");
                brick1.open("full_brick1.dat");
                brick2.open("full_brick2.dat");
                brick3.open("full_brick3.dat");
                for (unsigned j = 0; j < 27; j++)
                {
                  brick0 << brick_el0_pt->node_pt(j)->x(0) << " "
                         << brick_el0_pt->node_pt(j)->x(1) << " "
                         << brick_el0_pt->node_pt(j)->x(2) << "\n";

                  brick1 << brick_el1_pt->node_pt(j)->x(0) << " "
                         << brick_el1_pt->node_pt(j)->x(1) << " "
                         << brick_el1_pt->node_pt(j)->x(2) << "\n";

                  brick2 << brick_el2_pt->node_pt(j)->x(0) << " "
                         << brick_el2_pt->node_pt(j)->x(1) << " "
                         << brick_el2_pt->node_pt(j)->x(2) << "\n";

                  brick3 << brick_el3_pt->node_pt(j)->x(0) << " "
                         << brick_el3_pt->node_pt(j)->x(1) << " "
                         << brick_el3_pt->node_pt(j)->x(2) << "\n";
                }
                brick0.close();
                brick1.close();
                brick2.close();
                brick3.close();

                std::ofstream full_face;
                full_face.open("full_face.dat");
                for (unsigned j = 0; j < 6; j++)
                {
                  full_face << face_el_pt->node_pt(j)->x(0) << " "
                            << face_el_pt->node_pt(j)->x(1) << " "
                            << face_el_pt->node_pt(j)->x(2) << "\n";
                }
                full_face.close();

                // Get normal sign
                int normal_sign = face_el_pt->normal_sign();

                std::ostringstream error_stream;
                error_stream
                  << "During assignment of boundary cordinates, the distance\n"
                  << "between brick node and reference point in \n"
                  << "triangular FaceElement is " << dist << " which \n"
                  << "is bigger than the tolerance defined in \n"
                  << "BrickFromTetMeshHelper::Face_position_tolerance="
                  << BrickFromTetMeshHelper::Face_position_tolerance << ".\n"
                  << "If this is tolerable, increase the tolerance \n"
                  << "(it's defined in a namespace and therefore publically\n"
                  << "accessible). If not, the Face may be inverted in which \n"
                  << "case you should re-implement the translation scheme,\n"
                  << "following the pattern used in the "
                     "ThinLayerBrickOnTetMesh."
                  << "\nThe required code fragements are already here but \n"
                  << "the translation is the unit map.\n"
                  << "To aid the diagnostics, the files full_brick[0-3].dat\n"
                  << "contain the coordinates of the 27 nodes in the four\n"
                  << "bricks associated with the current tet and "
                     "full_face.dat\n"
                  << "contains the coordinates of the 6 nodes in the "
                     "FaceElement"
                  << "\nfrom which the boundary coordinates are extracted.\n"
                  << "FYI: The normal_sign of the face is: " << normal_sign
                  << std::endl;
                throw OomphLibError(error_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif

              // Set boundary stuff
              add_boundary_node(b, brick_node_pt);
              brick_node_pt->set_coordinates_on_boundary(b, zeta);
            }

            // Add appropriate brick elements to boundary lookup scheme
            switch (face_index)
            {
              case 0:
                Boundary_element_pt[b].push_back(brick_el1_pt);
                Face_index_at_boundary[b].push_back(-2);
                Boundary_element_pt[b].push_back(brick_el2_pt);
                Face_index_at_boundary[b].push_back(-1);
                Boundary_element_pt[b].push_back(brick_el3_pt);
                Face_index_at_boundary[b].push_back(-2);
                break;

              case 1:
                Boundary_element_pt[b].push_back(brick_el0_pt);
                Face_index_at_boundary[b].push_back(-1);
                Boundary_element_pt[b].push_back(brick_el2_pt);
                Face_index_at_boundary[b].push_back(-2);
                Boundary_element_pt[b].push_back(brick_el3_pt);
                Face_index_at_boundary[b].push_back(-1);
                break;

              case 2:
                Boundary_element_pt[b].push_back(brick_el0_pt);
                Face_index_at_boundary[b].push_back(-3);
                Boundary_element_pt[b].push_back(brick_el1_pt);
                Face_index_at_boundary[b].push_back(-3);
                Boundary_element_pt[b].push_back(brick_el2_pt);
                Face_index_at_boundary[b].push_back(-3);
                break;

              case 3:
                Boundary_element_pt[b].push_back(brick_el0_pt);
                Face_index_at_boundary[b].push_back(-2);
                Boundary_element_pt[b].push_back(brick_el1_pt);
                Face_index_at_boundary[b].push_back(-1);
                Boundary_element_pt[b].push_back(brick_el3_pt);
                Face_index_at_boundary[b].push_back(-3);
                break;
            }
            // Cleanup
            delete face_el_pt;
          }
        }
        // Cleanup
        delete face_pt;
      }
    }

    // Lookup scheme has now been setup
    Lookup_for_elements_next_boundary_is_setup = true;

    // Get number of distinct boundaries specified
    // in the original xda enumeration.
    //  unsigned n_xda_boundaries=tet_mesh_pt->nboundary();

    // Copy collective IDs across
    //  Boundary_id.resize(n_xda_boundaries);
    //  for (unsigned xda_b=0;xda_b<n_xda_boundaries;xda_b++)
    //   {
    //    //Boundary_id[xda_b]=tet_mesh_pt->oomph_lib_boundary_ids(xda_b);
    //    Boundary_id[xda_b]=xda_b;
    //   }


    // Cleanup
    for (unsigned e = 0; e < 4; e++)
    {
      for (unsigned j = 0; j < 8; j++)
      {
        delete dummy_q_el_pt[e]->node_pt(j);
      }
      delete dummy_q_el_pt[e];
    }
  }

} // namespace oomph

#endif
