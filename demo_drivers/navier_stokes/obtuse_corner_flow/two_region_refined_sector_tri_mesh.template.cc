// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
#ifndef TWO_REGION_REFINED_SECTOR_TRI_MESH_TEMPLATE_CC
#define TWO_REGION_REFINED_SECTOR_TRI_MESH_TEMPLATE_CC

#include <algorithm>

// Simple 2D triangle mesh class
#include "generic/Telements.h"
#include "two_region_refined_sector_tri_mesh.template.h"

namespace oomph
{
  //====================================================================
  /// Constructor for simple 2D triangular mesh class:
  ///
  /// n_radial    : number of elements in the radial direction;
  /// n_azimuthal : number of elements in the azimuthal direction;
  /// radius      : radius of the sector
  /// angle       : interior angle of the sector
  /// Ordering of elements: 'inner left' to 'inner right' then 'outward'
  ///
  //====================================================================
  template<class ELEMENT>
  TwoRegionRefinedSectorTriMesh<ELEMENT>::TwoRegionRefinedSectorTriMesh(
    const unsigned& n_radial,
    const unsigned& n_radial_region,
    const double& radial_geometric_sequence_base,
    const unsigned& n_azimuthal,
    const double& radius,
    const double& angle,
    TimeStepper* time_stepper_pt)
    : N_radial(n_radial),
      N_radial_region(n_radial_region),
      Radial_geometric_sequence_base(radial_geometric_sequence_base),
      N_azimuthal(n_azimuthal),
      Radius(radius),
      Angle(angle)
  {
    using namespace MathematicalConstants;

    // Mesh can only be built with 2D Telements.
    MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(2);

    // Set number of boundaries
    this->set_nboundary(6);

    // Allocate the store for the elements
    unsigned n_element = (N_radial - 1) * N_azimuthal * 2 + N_azimuthal;
    Element_pt.resize(n_element, 0);

    // Create first element
    Element_pt[0] = new ELEMENT;

    // Currently this mesh only works for 3 and 6 noded triangles
    if ((finite_element_pt(0)->nnode_1d() != 2) &&
        (finite_element_pt(0)->nnode_1d() != 3))
    {
      throw OomphLibError(
        "Currently this mesh only works for 3 & 6-noded triangles",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    unsigned n_node = 0;
    // Unless nnode_1d returned as !=2 default no extra nodes
    unsigned additional_nnode = 0;

    // Allocate storage if no extra nodes
    if (finite_element_pt(0)->nnode_1d() == 2)
    {
      n_node = (N_radial) * (N_azimuthal + 1) + 1;
    }

    if (finite_element_pt(0)->nnode_1d() == 3)
    {
      additional_nnode = 3;
      // Allocate storage for nodes (including extra)
      n_node = (N_radial) * (N_azimuthal + 1) + 1 +
               N_radial * (N_azimuthal + 1) + N_radial * N_azimuthal +
               (N_radial - 1) * (N_azimuthal);
    }

    Node_pt.resize(n_node);

    // Set up geometrical data
    //------------------------
    unsigned long node_count = 0;
    unsigned long element_count = 0;
    // Set the values of the increments
    double xstep = Radius / (N_radial);
    double ystep = Angle / (N_azimuthal);


    // FIRST NODE (lower left corner)
    //-------------------------------

    // Set the corner node

    // Allocate memory for the corner node of the first element
    //(which is not created loop as all subsequent bottom corners already exist)
    Node_pt[node_count] =
      finite_element_pt(0)->construct_node(0, time_stepper_pt);

    // Set the pointer from the element
    finite_element_pt(0)->node_pt(0) = Node_pt[node_count];

    // CURRENT (GLOBAL) NODE
    // Set the position
    Node_pt[node_count]->x(0) = map_to_sector_x(0, 0);
    Node_pt[node_count]->x(1) = map_to_sector_y(0, 0);

    // Add node to any relevant boundaries
    this->convert_to_boundary_node(Node_pt[node_count]);
    add_boundary_node(Inner_slip_boundary_id, Node_pt[node_count]);
    add_boundary_node(Inner_free_surface_boundary_id, Node_pt[node_count]);

    // Increment counter
    node_count++;

    // CREATE THE ELEMENTS (each has 3 local nodes)
    //--------------------------------------------------
    // Elements are created two at a time, the first being lower triangular
    // and the second upper triangular, so that they form a containing box.
    // Local nodes are numbered anti-clockwise with node_pt(1) being the
    // right-angle corner of the triangle.
    // The global node number, node_count, is considered to define
    // a containing box, with node_count defined as the node number
    // of the bottom left corner of the box.
    for (unsigned j = 0; j < N_azimuthal + 1; ++j)
    {
      // On the first row
      const unsigned i = 0;

      // CURRENT TRIANGLE
      if (j < N_azimuthal && i < N_radial)
      {
        // Create one elements
        //------------------------------
        if (element_count != 0) // 0th element already exists
        {
          Element_pt[element_count] = new ELEMENT;
        }

        // Allocate memory for nodes in the current triangle
        //--------------------------------------------
        // Allocate the top right corner node
        Node_pt[node_count + N_radial] =
          finite_element_pt(element_count)->construct_node(2, time_stepper_pt);

        // If on bottom row, construct the bottom right corner node
        if (j == 0)
        {
          Node_pt[node_count] = finite_element_pt(element_count)
                                  ->construct_node(1, time_stepper_pt);
        }

        // Bottom left corner node is already allocated

        // Set the pointers from the elements
        //----------------------------------
        finite_element_pt(element_count)->node_pt(0) = Node_pt[0];
        finite_element_pt(element_count)->node_pt(1) = Node_pt[node_count];
        finite_element_pt(element_count)->node_pt(2) =
          Node_pt[node_count + N_radial];

        element_count += 1;

        // Add node to any relevant boundaries
        if (j == 0)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          add_boundary_node(Inner_slip_boundary_id, Node_pt[node_count]);
        }
        if (j == N_azimuthal)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          add_boundary_node(Inner_free_surface_boundary_id,
                            Node_pt[node_count]);
        }
      }


      for (unsigned i = 1; i < N_radial + 1; ++i)
      {
        // CURRENT BOX
        //(nodes on RHS or top edge of domain do not define a box)
        if (j < N_azimuthal && i < N_radial)
        {
          // Create two elements
          //------------------------------
          Element_pt[element_count] = new ELEMENT;
          Element_pt[element_count + 1] = new ELEMENT;


          // Allocate memory for nodes in the current box
          //--------------------------------------------
          // If on bottom row, allocate the bottom right corner node
          if (j == 0)
          {
            Node_pt[node_count + 1] = finite_element_pt(element_count)
                                        ->construct_node(2, time_stepper_pt);
          }

          // Always need to allocate top right corner node
          Node_pt[node_count + N_radial + 1] =
            finite_element_pt(element_count + 1)
              ->construct_node(1, time_stepper_pt);

          // Bottom left corner node is already allocated


          // Set the pointers from the elements
          //----------------------------------
          finite_element_pt(element_count)->node_pt(0) =
            Node_pt[node_count + N_radial];
          finite_element_pt(element_count)->node_pt(1) = Node_pt[node_count];
          finite_element_pt(element_count)->node_pt(2) =
            Node_pt[node_count + 1];
          finite_element_pt(element_count + 1)->node_pt(0) =
            Node_pt[node_count + 1];
          finite_element_pt(element_count + 1)->node_pt(1) =
            Node_pt[node_count + N_radial + 1];
          finite_element_pt(element_count + 1)->node_pt(2) =
            Node_pt[node_count + N_radial];

          element_count += 2;
        }

        // CURRENT (GLOBAL) NODE
        // Set the position
        Node_pt[node_count]->x(0) = map_to_sector_x(i * xstep, j * ystep);
        Node_pt[node_count]->x(1) = map_to_sector_y(i * xstep, j * ystep);

        // Add node to any relevant boundaries
        if (j == 0)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          if (i <= N_radial_region)
          {
            add_boundary_node(Inner_slip_boundary_id, Node_pt[node_count]);
          }
          if (i >= N_radial_region)
          {
            add_boundary_node(Slip_boundary_id, Node_pt[node_count]);
          }
        }
        if (j == N_azimuthal)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          if (i <= N_radial_region)
          {
            add_boundary_node(Inner_free_surface_boundary_id,
                              Node_pt[node_count]);
          }
          if (i >= N_radial_region)
          {
            add_boundary_node(Free_surface_boundary_id, Node_pt[node_count]);
          }
        }
        if (i == N_radial_region)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          add_boundary_node(Inner_boundary_id, Node_pt[node_count]);
        }
        if (i == N_radial)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          add_boundary_node(Far_field_boundary_id, Node_pt[node_count]);
        }

        // Increment counter
        node_count++;
      }
    }


    if (additional_nnode == 3)
    {
      // Reset element counter to allow looping over existing elements
      // to add extra nodes.
      // Note that the node_count is not reset so no entries are overwritten
      element_count = 0;
      for (unsigned j = 0; j < N_azimuthal + 1; ++j)
      {
        const unsigned i = 0;

        // The additional nodes will be added in stages for ease of
        // application. Note: local numbering follows the anti-clockwise
        // pattern, node 3 halfway between  nodes 0-1, 4 between 1-2 and the
        // 5th local node between 2-0.

        // Restricted to stop creation of additional nodes outside the mesh
        if (j < N_azimuthal)
        {
          // If on the bottom row middle-node at bottom needs allocating
          if (j == 0)
          {
            Node_pt[node_count] = finite_element_pt(element_count)
                                    ->construct_node(3, time_stepper_pt);
          }
          // Always create the top node.
          Node_pt[node_count + N_radial] =
            finite_element_pt(element_count)
              ->construct_node(5, time_stepper_pt);

          // Set pointers to the top/bottom middle nodes
          // Bottom node
          finite_element_pt(element_count)->node_pt(3) = Node_pt[node_count];
          finite_element_pt(element_count)->node_pt(5) =
            Node_pt[node_count + N_radial];

          // Increase the element counter to add top/bottom nodes to
          // next pair of element on next pass
          element_count += 1;
        } // End extra top/bottom node construction

        // Set position of current (Global) Node
        Node_pt[node_count]->x(0) =
          map_to_sector_x(double(i + 0.5) * xstep, j * ystep);
        Node_pt[node_count]->x(1) =
          map_to_sector_y(double(i + 0.5) * xstep, j * ystep);

        // Add node to any applicable boundaries (node 4's can only be top
        // or bottom boundaries)
        if (j == 0)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          if (i < N_radial_region)
          {
            add_boundary_node(Inner_slip_boundary_id, Node_pt[node_count]);
          }
          if (i >= N_radial_region)
          {
            add_boundary_node(Slip_boundary_id, Node_pt[node_count]);
          }
        }
        if (j == N_azimuthal)
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
          if (i < N_radial_region)
          {
            add_boundary_node(Inner_free_surface_boundary_id,
                              Node_pt[node_count]);
          }
          if (i >= N_radial_region)
          {
            add_boundary_node(Free_surface_boundary_id, Node_pt[node_count]);
          }
        }

        // Update node_count
        node_count++;

        // Note: i counter reduced by 1 since i axis runs through middle of
        // elements on x-axis
        for (unsigned i = 1; i < N_radial; ++i)
        {
          // The additional nodes will be added in stages for ease of
          // application. Note: local numbering follows the anti-clockwise
          // pattern, node 3 halfway between  nodes 0-1, 4 between 1-2 and the
          // 5th local node between 2-0.

          // Restricted to stop creation of additional nodes outside the mesh
          if (j < N_azimuthal)
          {
            // If on the bottom row middle-node at bottom needs allocating
            if (j == 0)
            {
              Node_pt[node_count] = finite_element_pt(element_count)
                                      ->construct_node(4, time_stepper_pt);
            }

            // Due to directions of iteration node at middle of top box edge
            // (in next element) always needs allocating
            Node_pt[node_count + N_radial] =
              finite_element_pt(element_count + 1)
                ->construct_node(4, time_stepper_pt);

            // Set pointers to the top/bottom middle nodes
            // Bottom node
            finite_element_pt(element_count)->node_pt(4) = Node_pt[node_count];
            // Top node
            finite_element_pt(element_count + 1)->node_pt(4) =
              Node_pt[node_count + N_radial];

            // Increase the element counter to add top/bottom nodes to
            // next pair of element on next pass
            element_count += 2;
          } // End extra top/bottom node construction

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) =
            map_to_sector_x((double(i) + 0.5) * xstep, j * ystep);
          Node_pt[node_count]->x(1) =
            map_to_sector_y((double(i) + 0.5) * xstep, j * ystep);

          // Add node to any applicable boundaries (node 4's can only be top
          // or bottom boundaries)
          if (j == 0)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            if (i < N_radial_region)
            {
              add_boundary_node(Inner_slip_boundary_id, Node_pt[node_count]);
            }
            if (i >= N_radial_region)
            {
              add_boundary_node(Slip_boundary_id, Node_pt[node_count]);
            }
          }
          if (j == N_azimuthal)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            this->convert_to_boundary_node(Node_pt[node_count]);
            if (i < N_radial_region)
            {
              add_boundary_node(Inner_free_surface_boundary_id,
                                Node_pt[node_count]);
            }
            if (i >= N_radial_region)
            {
              add_boundary_node(Free_surface_boundary_id, Node_pt[node_count]);
            }
          }

          // Update node_count
          node_count++;
        }
      }

      // Next stage of additional node implementation involes the middle left
      // and right nodes (local number 3 on each triangle)

      // Element counter reset for second loop over existing elements
      element_count = 0;
      // Note: j counter reduced by 1 since j axis runs through middle of
      // elements on y-axis
      for (unsigned j = 0; j < N_azimuthal; ++j)
      {
        const unsigned i = 1;

        if (j < N_azimuthal && i < N_radial + 1)
        {
          // Make the node on the right hand side
          Node_pt[node_count] = finite_element_pt(element_count)
                                  ->construct_node(4, time_stepper_pt);

          // Set pointers from the elements to new nodes
          finite_element_pt(element_count)->node_pt(4) = Node_pt[node_count];

          // Increase element_count by 1
          element_count += 1;
        } // End extra left/right node construction

        // Set position of current (Global) Node
        Node_pt[node_count]->x(0) =
          map_to_sector_x(double(i) * xstep, (double(j) + 0.5) * ystep);
        Node_pt[node_count]->x(1) =
          map_to_sector_y(double(i) * xstep, (double(j) + 0.5) * ystep);

        // Update node_count
        node_count++;

        // Loop from left to right
        for (unsigned i = 2; i < N_radial + 1; ++i)
        {
          if (j < N_azimuthal && i < N_radial + 1)
          {
            // The mid node on the right hand side always needs constructing
            // within the bounds of the mesh
            Node_pt[node_count] = finite_element_pt(element_count + 1)
                                    ->construct_node(3, time_stepper_pt);

            // Set pointers from the elements to new nodes
            finite_element_pt(element_count)->node_pt(3) =
              Node_pt[node_count - 1];
            finite_element_pt(element_count + 1)->node_pt(3) =
              Node_pt[node_count];

            // Increase element_count by 2
            element_count += 2;
          } // End extra left/right node construction

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) =
            map_to_sector_x(double(i) * xstep, (double(j) + 0.5) * ystep);
          Node_pt[node_count]->x(1) =
            map_to_sector_y(double(i) * xstep, (double(j) + 0.5) * ystep);

          // Add node to any applicable boundaries again - only be left/right
          if (i == N_radial_region)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(Inner_boundary_id, Node_pt[node_count]);
          }
          if (i == N_radial)
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
            add_boundary_node(Far_field_boundary_id, Node_pt[node_count]);
          }

          // Update node_count
          node_count++;
        }
      }

      // Final stage of inserting extra nodes - inclusion of the local
      // number 5 (middle of hypot. edge)

      element_count = 0;
      // Note: both i,j counters reduced by 1 since j axis runs through middle
      // of elements in both x,y
      for (unsigned j = 0; j < N_azimuthal; ++j)
      {
        // const unsigned i = 0;

        // Increase element_count by 1
        element_count += 1;

        for (unsigned i = 1; i < N_radial; ++i)
        {
          // The middle node always needs constructing
          Node_pt[node_count] = finite_element_pt(element_count)
                                  ->construct_node(5, time_stepper_pt);

          // Set pointers from the elements to new nodes
          finite_element_pt(element_count)->node_pt(5) = Node_pt[node_count];
          finite_element_pt(element_count + 1)->node_pt(5) =
            Node_pt[node_count];

          // Increase element_count by 2
          element_count += 2;
          // End extra left/right node construction

          // Set position of current (Global) Node
          Node_pt[node_count]->x(0) = map_to_sector_x(
            (double(i) + 0.5) * xstep, (double(j) + 0.5) * ystep);
          Node_pt[node_count]->x(1) = map_to_sector_y(
            (double(i) + 0.5) * xstep, (double(j) + 0.5) * ystep);

          // All nodes are internal in this structure so no boundaries
          // applicable

          // Update node_count
          node_count++;
        }
      }
    }
    // End of extra nodes for 6 noded trianglur elements

    setup_boundary_element_info();
  }

  //================================================================
  /// Setup lookup schemes which establish which elements are located
  /// next to which boundaries (Doc to outfile if it's open).
  //================================================================
  template<class ELEMENT>
  void TwoRegionRefinedSectorTriMesh<ELEMENT>::setup_boundary_element_info(
    std::ostream& outfile)
  {
    // Should we document the output here
    bool doc = false;

    if (outfile) doc = true;

    // Number of boundaries
    unsigned nbound = nboundary();

    // Wipe/allocate storage for arrays
    Boundary_element_pt.clear();
    Face_index_at_boundary.clear();
    Boundary_element_pt.resize(nbound);
    Face_index_at_boundary.resize(nbound);

    // Temporary vector of vectors of pointers to elements on the boundaries:
    // This is a vector to ensure that order is strictly preserved
    Vector<Vector<FiniteElement*>> vector_of_boundary_element_pt;
    vector_of_boundary_element_pt.resize(nbound);

    // Matrix map for working out the fixed face for elements on boundary
    MapMatrixMixed<unsigned, FiniteElement*, int> face_identifier;

    // Loop over elements
    //-------------------
    unsigned nel = nelement();

    // Get pointer to vector of boundaries that the
    // node lives on
    Vector<std::set<unsigned>*> boundaries_pt(3, 0);

    // Data needed to deal with edges through the
    // interior of the domain
    std::map<Edge, unsigned> edge_count;
    std::map<Edge, TriangleBoundaryHelper::BCInfo> edge_bcinfo;
    std::map<Edge, TriangleBoundaryHelper::BCInfo> face_info;
    MapMatrixMixed<unsigned, FiniteElement*, int> face_count;
    Vector<unsigned> bonus(nbound);

    // When using internal boundaries, an edge can be related to more than
    // one element (because of both sides of the internal boundaries)
    std::map<Edge, Vector<TriangleBoundaryHelper::BCInfo>> edge_internal_bnd;

    for (unsigned e = 0; e < nel; e++)
    {
      // Get pointer to element
      FiniteElement* fe_pt = finite_element_pt(e);

      if (doc)
      {
        outfile << "Element: " << e << " " << fe_pt << std::endl;
      }

      // Only include 2D elements! Some meshes contain interface elements too.
      if (fe_pt->dim() == 2)
      {
        // Loop over the element's nodes and find out which boundaries they're
        // on
        // ----------------------------------------------------------------------

        // We need only loop over the corner nodes
        for (unsigned i = 0; i < 3; i++)
        {
          fe_pt->node_pt(i)->get_boundaries_pt(boundaries_pt[i]);
        }

        // Find the common boundaries of each edge
        Vector<std::set<unsigned>> edge_boundary(3);

        // Edge 0 connects points 1 and 2
        //-----------------------------

        if (boundaries_pt[1] && boundaries_pt[2])
        {
          // Create the corresponding edge
          Edge edge0(fe_pt->node_pt(1), fe_pt->node_pt(2));

          // Update infos about this edge
          TriangleBoundaryHelper::BCInfo info;
          info.Face_id = 0;
          info.FE_pt = fe_pt;

          std::set_intersection(boundaries_pt[1]->begin(),
                                boundaries_pt[1]->end(),
                                boundaries_pt[2]->begin(),
                                boundaries_pt[2]->end(),
                                std::insert_iterator<std::set<unsigned>>(
                                  edge_boundary[0], edge_boundary[0].begin()));
          std::set<unsigned>::iterator it0 = edge_boundary[0].begin();

          // Edge does exist:
          if (edge_boundary[0].size() > 0)
          {
            info.Boundary = *it0;

            // How many times this edge has been visited
            edge_count[edge0]++;

            // Update edge_bcinfo
            edge_bcinfo.insert(std::make_pair(edge0, info));

            // ... and also update the info associated with internal bnd
            edge_internal_bnd[edge0].push_back(info);
          }
        }

        // Edge 1 connects points 0 and 2
        //-----------------------------

        if (boundaries_pt[0] && boundaries_pt[2])
        {
          std::set_intersection(boundaries_pt[0]->begin(),
                                boundaries_pt[0]->end(),
                                boundaries_pt[2]->begin(),
                                boundaries_pt[2]->end(),
                                std::insert_iterator<std::set<unsigned>>(
                                  edge_boundary[1], edge_boundary[1].begin()));

          // Create the corresponding edge
          Edge edge1(fe_pt->node_pt(0), fe_pt->node_pt(2));

          // Update infos about this edge
          TriangleBoundaryHelper::BCInfo info;
          info.Face_id = 1;
          info.FE_pt = fe_pt;
          std::set<unsigned>::iterator it1 = edge_boundary[1].begin();

          // Edge does exist:
          if (edge_boundary[1].size() > 0)
          {
            info.Boundary = *it1;

            // How many times this edge has been visited
            edge_count[edge1]++;

            // Update edge_bcinfo
            edge_bcinfo.insert(std::make_pair(edge1, info));

            // ... and also update the info associated with internal bnd
            edge_internal_bnd[edge1].push_back(info);
          }
        }

        // Edge 2 connects points 0 and 1
        //-----------------------------

        if (boundaries_pt[0] && boundaries_pt[1])
        {
          std::set_intersection(boundaries_pt[0]->begin(),
                                boundaries_pt[0]->end(),
                                boundaries_pt[1]->begin(),
                                boundaries_pt[1]->end(),
                                std::insert_iterator<std::set<unsigned>>(
                                  edge_boundary[2], edge_boundary[2].begin()));

          // Create the corresponding edge
          Edge edge2(fe_pt->node_pt(0), fe_pt->node_pt(1));

          // Update infos about this edge
          TriangleBoundaryHelper::BCInfo info;
          info.Face_id = 2;
          info.FE_pt = fe_pt;
          std::set<unsigned>::iterator it2 = edge_boundary[2].begin();

          // Edge does exist:
          if (edge_boundary[2].size() > 0)
          {
            info.Boundary = *it2;

            // How many times this edge has been visited
            edge_count[edge2]++;

            // Update edge_bcinfo
            edge_bcinfo.insert(std::make_pair(edge2, info));

            // ... and also update the info associated with internal bnd
            edge_internal_bnd[edge2].push_back(info);
          }
        }


#ifdef PARANOID

        // Check if edge is associated with multiple boundaries

        // We now know whether any edges lay on the boundaries
        for (unsigned i = 0; i < 3; i++)
        {
          // How many boundaries are there
          unsigned count = 0;

          // Loop over all the members of the set and add to the count
          // and set the boundary
          for (std::set<unsigned>::iterator it = edge_boundary[i].begin();
               it != edge_boundary[i].end();
               ++it)
          {
            ++count;
          }

          // If we're on more than one boundary, this is weird, so die
          if (count > 1)
          {
            std::ostringstream error_stream;
            error_stream << "Edge " << i << " is located on " << count
                         << " boundaries.\n";
            error_stream << "This is rather strange, so I'm going to die\n";
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }

#endif

        // Now we set the pointers to the boundary sets to zero
        for (unsigned i = 0; i < 3; i++)
        {
          boundaries_pt[i] = 0;
        }
      }
    } // end of loop over all elements

    // Loop over all edges that are located on a boundary
    typedef std::map<Edge, TriangleBoundaryHelper::BCInfo>::iterator ITE;
    for (ITE it = edge_bcinfo.begin(); it != edge_bcinfo.end(); it++)
    {
      Edge current_edge = it->first;
      unsigned bound = it->second.Boundary;

      // If the edge has been visited only once
      if (edge_count[current_edge] == 1)
      {
        // Count the edges that are on the same element and on the same boundary
        face_count(static_cast<unsigned>(bound), it->second.FE_pt) =
          face_count(static_cast<unsigned>(bound), it->second.FE_pt) + 1;

        // If such edges exist, let store the corresponding element
        if (face_count(bound, it->second.FE_pt) > 1)
        {
          // Update edge's infos
          TriangleBoundaryHelper::BCInfo info;
          info.Face_id = it->second.Face_id;
          info.FE_pt = it->second.FE_pt;
          info.Boundary = it->second.Boundary;

          // Add it to FIinfo, that stores infos of problematic elements
          face_info.insert(std::make_pair(current_edge, info));

          // How many edges on which boundary have to be added
          bonus[bound]++;
        }
        else
        {
          // Add element and face to the appropriate vectors
          // Does the pointer already exits in the vector
          Vector<FiniteElement*>::iterator b_el_it = std::find(
            vector_of_boundary_element_pt[static_cast<unsigned>(bound)].begin(),
            vector_of_boundary_element_pt[static_cast<unsigned>(bound)].end(),
            it->second.FE_pt);

          // Only insert if we have not found it (i.e. got to the end)
          if (b_el_it ==
              vector_of_boundary_element_pt[static_cast<unsigned>(bound)].end())
          {
            vector_of_boundary_element_pt[static_cast<unsigned>(bound)]
              .push_back(it->second.FE_pt);
          }

          // set_of_boundary_element_pt[static_cast<unsigned>(bound)].insert(
          // it->second.FE_pt);
          face_identifier(static_cast<unsigned>(bound), it->second.FE_pt) =
            it->second.Face_id;
        }
      }

    } // End of "adding-boundaries"-loop


    // Now copy everything across into permanent arrays
    //-------------------------------------------------

    // Loop over boundaries
    for (unsigned i = 0; i < nbound; i++)
    {
      // Number of elements on this boundary that have to be added
      // in addition to other elements
      unsigned bonus1 = bonus[i];

      // Number of elements on this boundary
      unsigned nel = vector_of_boundary_element_pt[i].size() + bonus1;

      // Allocate storage for the coordinate identifiers
      Face_index_at_boundary[i].resize(nel);

      unsigned e_count = 0;
      typedef Vector<FiniteElement*>::iterator IT;
      for (IT it = vector_of_boundary_element_pt[i].begin();
           it != vector_of_boundary_element_pt[i].end();
           it++)
      {
        // Recover pointer to element
        FiniteElement* fe_pt = *it;

        // Add to permanent storage
        Boundary_element_pt[i].push_back(fe_pt);

        Face_index_at_boundary[i][e_count] = face_identifier(i, fe_pt);

        // Increment counter
        e_count++;
      }
      // We add the elements that have two or more edges on this boundary
      for (ITE itt = face_info.begin(); itt != face_info.end(); itt++)
      {
        if (itt->second.Boundary == i)
        {
          // Add to permanent storage
          Boundary_element_pt[i].push_back(itt->second.FE_pt);

          Face_index_at_boundary[i][e_count] = itt->second.Face_id;

          e_count++;
        }
      }

    } // End of loop over boundaries

    // Doc?
    //-----
    if (doc)
    {
      // Loop over boundaries
      for (unsigned i = 0; i < nbound; i++)
      {
        unsigned nel = Boundary_element_pt[i].size();
        outfile << "Boundary: " << i << " is adjacent to " << nel << " elements"
                << std::endl;

        // Loop over elements on given boundary
        for (unsigned e = 0; e < nel; e++)
        {
          FiniteElement* fe_pt = Boundary_element_pt[i][e];
          outfile << "Boundary element:" << fe_pt
                  << " Face index of boundary is "
                  << Face_index_at_boundary[i][e] << std::endl;
        }
      }
    }

    // Lookup scheme has now been setup yet
    Lookup_for_elements_next_boundary_is_setup = true;
  }

} // namespace oomph
#endif
