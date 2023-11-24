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
#ifndef OOMPH_QUARTER_CIRCLE_SECTOR_MESH_TEMPLATE_CC
#define OOMPH_QUARTER_CIRCLE_SECTOR_MESH_TEMPLATE_CC


#include "quarter_circle_sector_mesh.template.h"

namespace oomph
{
  //====================================================================
  /// Constructor for deformable 2D Ring mesh class. Pass pointer to
  /// geometric object that specifies the wall, start and end coordinates on the
  /// geometric object, and the fraction along
  /// which the dividing line is to be placed, and the timestepper
  /// (defaults to (Steady) default timestepper defined in Mesh).
  /// Nodal positions are determined via macro-element-based representation
  /// of the Domain (as a QuarterCircleSectorDomain).
  //====================================================================
  template<class ELEMENT>
  QuarterCircleSectorMesh<ELEMENT>::QuarterCircleSectorMesh(
    GeomObject* wall_pt,
    const double& xi_lo,
    const double& fract_mid,
    const double& xi_hi,
    TimeStepper* time_stepper_pt)
    : Wall_pt(wall_pt), Xi_lo(xi_lo), Fract_mid(fract_mid), Xi_hi(xi_hi)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Build macro element-based domain
    Domain_pt = new QuarterCircleSectorDomain(wall_pt, xi_lo, fract_mid, xi_hi);

    // Set the number of boundaries
    set_nboundary(3);

    // We have only bothered to parametrise boundary 1
    Boundary_coordinate_exists[1] = true;

    // Allocate the store for the elements
    Element_pt.resize(3);

    // Create first element
    Element_pt[0] = new ELEMENT;

    // Read out the number of linear points in the element
    unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

    // Can now allocate the store for the nodes
    Node_pt.resize(n_p * n_p + (n_p - 1) * n_p + (n_p - 1) * (n_p - 1));


    Vector<double> s(2);
    Vector<double> r(2);

    // Storage for the intrinsic boundary coordinate
    Vector<double> zeta(1);


    // Set up geometrical data
    //------------------------

    // Initialise node counter
    unsigned long node_count = 0;


    // Now assign the topology
    // Boundaries are numbered 0 1 2 from the bottom proceeding anticlockwise


    // FIRST ELEMENT (lower left corner)
    //

    // Set the corner node (on boundaries 0 and 2)
    //-------------------------------------------

    // Create the ll node
    Node_pt[node_count] =
      finite_element_pt(0)->construct_boundary_node(0, time_stepper_pt);

    // Set the pointer from the element to the node
    finite_element_pt(0)->node_pt(0) = Node_pt[node_count];

    // Set the position of the ll node
    s[0] = -1.0;
    s[1] = -1.0;
    Domain_pt->macro_element_pt(0)->macro_map(s, r);
    Node_pt[node_count]->x(0) = r[0];
    Node_pt[node_count]->x(1) = r[1];

    // Add the node to the boundaries
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(2, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // First row is on boundary 0:
    //---------------------------
    for (unsigned l1 = 1; l1 < n_p; l1++)
    {
      // Local node number
      unsigned jnod_local = l1;

      // Create the node
      Node_pt[node_count] = finite_element_pt(0)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(0)->node_pt(jnod_local) = Node_pt[node_count];

      // Set the position of the node
      s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
      s[1] = -1.0;
      Domain_pt->macro_element_pt(0)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(0, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }


    // Loop over the other rows of nodes
    //------------------------------------
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // First node in this row is on boundary 2:
      //-----------------------------------------

      // Local node number
      unsigned jnod_local = n_p * l2;

      // Create the node
      Node_pt[node_count] = finite_element_pt(0)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(0)->node_pt(jnod_local) = Node_pt[node_count];


      // Set the position of the node
      s[0] = -1.0;
      s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
      ;
      Domain_pt->macro_element_pt(0)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(2, Node_pt[node_count]);

      // Increment the node number
      node_count++;


      // The other nodes are in the interior
      //------------------------------------
      // Loop over the other node columns
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // Local node number
        unsigned jnod_local = l1 + n_p * l2;

        // Create the node
        Node_pt[node_count] =
          finite_element_pt(0)->construct_node(jnod_local, time_stepper_pt);

        // Set the pointer from the element to the node
        finite_element_pt(0)->node_pt(jnod_local) = Node_pt[node_count];

        // Set the position of the node
        s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
        s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
        Domain_pt->macro_element_pt(0)->macro_map(s, r);
        Node_pt[node_count]->x(0) = r[0];
        Node_pt[node_count]->x(1) = r[1];

        // Increment the node number
        node_count++;
      }
    }

    // SECOND ELEMENT (lower right corner)
    //
    // Create element
    Element_pt[1] = new ELEMENT;

    // Loop over the first column (already exists!)
    //---------------------------------------------
    for (unsigned l2 = 0; l2 < n_p; l2++)
    {
      // Node number in existing element
      unsigned jnod_local_old = (n_p - 1) + l2 * n_p;

      // Set the pointer from the element to the node
      finite_element_pt(1)->node_pt(l2 * n_p) =
        finite_element_pt(0)->node_pt(jnod_local_old);
    }

    // Loop over the other node columns (apart from last one)
    //------------------------------------------------------
    for (unsigned l1 = 1; l1 < n_p - 1; l1++)
    {
      // First node is at the bottom (on boundary 0)
      //--------------------------------------------

      // Local node number
      unsigned jnod_local = l1;

      // Create the node
      Node_pt[node_count] = finite_element_pt(1)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(1)->node_pt(jnod_local) = Node_pt[node_count];

      // Set the position of the node
      s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
      s[1] = -1.0;
      Domain_pt->macro_element_pt(1)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(0, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Now loop over the interior nodes in this column
      //-------------------------------------------------
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Local node number
        unsigned jnod_local = l1 + l2 * n_p;

        // Create the node
        Node_pt[node_count] =
          finite_element_pt(1)->construct_node(jnod_local, time_stepper_pt);

        // Set the pointer from the element to the node
        finite_element_pt(1)->node_pt(jnod_local) = Node_pt[node_count];

        // Set the position of the node
        s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
        s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
        Domain_pt->macro_element_pt(1)->macro_map(s, r);
        Node_pt[node_count]->x(0) = r[0];
        Node_pt[node_count]->x(1) = r[1];

        // Increment the node number
        node_count++;
      }
    }

    // Last column (on boundary 1)
    //----------------------------

    // First node is at the bottom (and hence also on boundary 0)
    //-----------------------------------------------------------

    // Local node number
    unsigned jnod_local = n_p - 1;

    // Create the node
    Node_pt[node_count] = finite_element_pt(1)->construct_boundary_node(
      jnod_local, time_stepper_pt);

    // Set the pointer from the element to the node
    finite_element_pt(1)->node_pt(jnod_local) = Node_pt[node_count];

    // Set the position of the node
    s[0] = 1.0;
    s[1] = -1.0;
    Domain_pt->macro_element_pt(1)->macro_map(s, r);
    Node_pt[node_count]->x(0) = r[0];
    Node_pt[node_count]->x(1) = r[1];

    // Add the node to the boundaries
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(1, Node_pt[node_count]);

    // Set the intrinsic coordinate on the boundary 1
    zeta[0] = Xi_lo + 0.5 * (1.0 + s[1]) * Fract_mid * (Xi_hi - Xi_lo);
    Node_pt[node_count]->set_coordinates_on_boundary(1, zeta);

    // Increment the node number
    node_count++;

    // Now do the remaining nodes in last column (only on boundary 1)
    //---------------------------------------------------------------
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // Local node number
      unsigned jnod_local = (n_p - 1) + l2 * n_p;

      // Create the node
      Node_pt[node_count] = finite_element_pt(1)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(1)->node_pt(jnod_local) = Node_pt[node_count];

      // Set the position of the node
      s[0] = 1.0;
      s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
      Domain_pt->macro_element_pt(1)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(1, Node_pt[node_count]);

      // Set the intrinsic coordinate on the boundary 1
      zeta[0] = Xi_lo + 0.5 * (1.0 + s[1]) * Fract_mid * (Xi_hi - Xi_lo);
      Node_pt[node_count]->set_coordinates_on_boundary(1, zeta);


      // Increment the node number
      node_count++;
    }


    // THIRD ELEMENT (upper left corner)
    //
    // Create element
    Element_pt[2] = new ELEMENT;

    // Loop over the first row (has already been created via element 0)
    //-----------------------------------------------------------------
    for (unsigned l1 = 0; l1 < n_p; l1++)
    {
      // Node number in existing element
      unsigned jnod_local_old = n_p * (n_p - 1) + l1;

      // Local node number here
      unsigned jnod_local = l1;

      // Set the pointer from the element to the node
      finite_element_pt(2)->node_pt(jnod_local) =
        finite_element_pt(0)->node_pt(jnod_local_old);
    }


    // Loop over the remaining nodes in the last column (has already
    //--------------------------------------------------------------
    // been created via element 1)
    //----------------------------
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // Node number in existing element
      unsigned jnod_local_old = n_p * (n_p - 1) + l2;

      // Local node number here
      unsigned jnod_local = (n_p - 1) + l2 * n_p;

      // Set the pointer from the element to the node
      finite_element_pt(2)->node_pt(jnod_local) =
        finite_element_pt(1)->node_pt(jnod_local_old);
    }


    // Loop over the nodes in rows (apart from last one which is on boundary 1)
    //-------------------------------------------------------------------------
    for (unsigned l2 = 1; l2 < n_p - 1; l2++)
    {
      // First node in this row is on boundary 2:
      //-----------------------------------------

      // Local node number
      unsigned jnod_local = n_p * l2;

      // Create the node
      Node_pt[node_count] = finite_element_pt(2)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(2)->node_pt(jnod_local) = Node_pt[node_count];

      // Set the position of the node
      s[0] = -1.0;
      s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
      ;
      Domain_pt->macro_element_pt(2)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(2, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // The other nodes are in the interior
      //------------------------------------
      // Loop over the other node columns
      for (unsigned l1 = 1; l1 < n_p - 1; l1++)
      {
        // Local node number
        unsigned jnod_local = l1 + n_p * l2;

        // Create the node
        Node_pt[node_count] =
          finite_element_pt(2)->construct_node(jnod_local, time_stepper_pt);

        // Set the pointer from the element to the node
        finite_element_pt(2)->node_pt(jnod_local) = Node_pt[node_count];

        // Set the position of the node
        s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
        s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
        Domain_pt->macro_element_pt(2)->macro_map(s, r);
        Node_pt[node_count]->x(0) = r[0];
        Node_pt[node_count]->x(1) = r[1];

        // Increment the node number
        node_count++;
      }
    }


    // Top left corner is on boundaries 1 and 2:
    //------------------------------------------

    // Local node number
    jnod_local = n_p * (n_p - 1);

    // Create the node
    Node_pt[node_count] = finite_element_pt(2)->construct_boundary_node(
      jnod_local, time_stepper_pt);

    // Set the pointer from the element to the node
    finite_element_pt(2)->node_pt(jnod_local) = Node_pt[node_count];

    // Set the position of the node
    s[0] = -1.0;
    s[1] = 1.0;
    Domain_pt->macro_element_pt(2)->macro_map(s, r);
    Node_pt[node_count]->x(0) = r[0];
    Node_pt[node_count]->x(1) = r[1];

    // Add the node to the boundaries
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(2, Node_pt[node_count]);

    // Set the intrinsic coordinate on the boundary 1
    zeta[0] = Xi_hi + 0.5 * (s[0] + 1.0) * (1.0 - Fract_mid) * (Xi_lo - Xi_hi);
    Node_pt[node_count]->set_coordinates_on_boundary(1, zeta);


    // Increment the node number
    node_count++;

    // Rest of top row is on boundary 1 only:
    //---------------------------------------
    for (unsigned l1 = 1; l1 < n_p - 1; l1++)
    {
      // Local node number
      unsigned jnod_local = n_p * (n_p - 1) + l1;

      // Create the node
      Node_pt[node_count] = finite_element_pt(2)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(2)->node_pt(jnod_local) = Node_pt[node_count];

      // Set the position of the node
      s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
      s[1] = 1.0;
      Domain_pt->macro_element_pt(2)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(1, Node_pt[node_count]);

      // Set the intrinsic coordinate on the boundary 1
      zeta[0] =
        Xi_hi + 0.5 * (s[0] + 1.0) * (1.0 - Fract_mid) * (Xi_lo - Xi_hi);
      Node_pt[node_count]->set_coordinates_on_boundary(1, zeta);

      // Increment the node number
      node_count++;
    }

    // Loop over all elements and set macro element pointer
    // to enable MacroElement-based node update
    unsigned n_element = this->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get pointer to full element type
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(this->element_pt(e));

      // Set pointer to macro element
      el_pt->set_macro_elem_pt(this->Domain_pt->macro_element_pt(e));
    }

    // Setup boundary element lookup schemes
    setup_boundary_element_info();
  }


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // Algebraic-mesh-version of RefineableQuarterCircleSectorMesh
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Setup algebraic update operation, based on individual
  /// blocks. Mesh is suspended from the `wall' GeomObject pointed to
  /// by wall_pt. The lower right corner of the mesh is located at the
  /// wall's coordinate xi_lo, the upper left corner at xi_hi;
  /// The dividing line between the two blocks at the outer ring
  /// is located at the fraction fract_mid between these two coordinates.
  /// Node updating strategy:
  /// - the shape of the central box remains rectangular; the position
  ///   of its top right corner is located at a fixed fraction
  ///   of the width and height of the domain. Nodes in this region
  ///   are located at fixed horizontal and vertical fractions of the box.
  /// - Nodes in the two outer "ring" elements (bottom right and top left)
  ///   are located on straight lines running from the edges of the
  ///   central box to the outer wall.
  //======================================================================
  template<class ELEMENT>
  void AlgebraicRefineableQuarterCircleSectorMesh<
    ELEMENT>::setup_algebraic_node_update()
  {
#ifdef PARANOID
    /// Pointer to algebraic element in central box
    AlgebraicElementBase* central_box_pt =
      dynamic_cast<AlgebraicElementBase*>(Mesh::element_pt(0));

    if (central_box_pt == 0)
    {
      std::ostringstream error_message;
      error_message
        << "Element in AlgebraicRefineableQuarterCircleSectorMesh must be\n ";
      error_message << "derived from AlgebraicElementBase\n";
      error_message << "but it is of type:  "
                    << typeid(Mesh::element_pt(0)).name() << std::endl;
      std::string function_name =
        "AlgebraicRefineableQuarterCircleSectorMesh::";
      function_name += "setup_algebraic_node_update()";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Number of nodes in elements
    unsigned nnodes = Mesh::finite_element_pt(0)->nnode();

    // Coordinates of central box
    double x_box = Mesh::finite_element_pt(0)->node_pt(nnodes - 1)->x(0);
    double y_box = Mesh::finite_element_pt(0)->node_pt(nnodes - 1)->x(1);

    // Get wall position in bottom right corner
    Vector<double> r_br(2);
    Vector<double> xi(1);
    xi[0] = QuarterCircleSectorMesh<ELEMENT>::Xi_lo;
    QuarterCircleSectorMesh<ELEMENT>::Wall_pt->position(xi, r_br);

    // Establish fractional width of central box
    Lambda_x = x_box / r_br[0];

    // Find corresponding wall element/local coordinate
    GeomObject* obj_br_pt;
    Vector<double> s_br(1);
    this->Wall_pt->locate_zeta(xi, obj_br_pt, s_br);

    obj_br_pt->position(s_br, r_br);

    // Get wall position in top left corner
    Vector<double> r_tl(2);
    xi[0] = QuarterCircleSectorMesh<ELEMENT>::Xi_hi;
    QuarterCircleSectorMesh<ELEMENT>::Wall_pt->position(xi, r_tl);

    // Establish fractional height of central box
    Lambda_y = y_box / r_tl[1];

    // Find corresponding wall element/local coordinate
    GeomObject* obj_tl_pt;
    Vector<double> s_tl(1);
    this->Wall_pt->locate_zeta(xi, obj_tl_pt, s_tl);


    // Element 0: central box
    //-----------------------
    {
      unsigned ielem = 0;
      FiniteElement* el_pt = Mesh::finite_element_pt(ielem);

      // Loop over all nodes in the element and give them the update
      // info appropriate for the current element
      for (unsigned jnod = 0; jnod < nnodes; jnod++)
      {
        // Nodal coordinates in undeformed mesh
        double x = Mesh::finite_element_pt(ielem)->node_pt(jnod)->x(0);
        double y = Mesh::finite_element_pt(ielem)->node_pt(jnod)->x(1);

        // The update function requires two geometric objects
        Vector<GeomObject*> geom_object_pt(2);

        // The update function requires four parameters:
        Vector<double> ref_value(4);

        // First reference value: fractional x-position inside box
        ref_value[0] = x / x_box;

        // Second reference value: fractional y-position inside box
        ref_value[1] = y / y_box;

        // Wall element at bottom right end of wall mesh:
        geom_object_pt[0] = obj_br_pt;

        // Local coordinate in this wall element (Note:
        // we'll have to recompute this reference
        // when the mesh is refined  as we might fall off the element otherwise)
        ref_value[2] = s_br[0];


        // Wall element at top left end of wall mesh:
        geom_object_pt[1] = obj_tl_pt;

        // Local coordinate in this wall element. Note:
        // we'll have to recompute this reference
        // when the mesh is refined  as we might fall off the element otherwise)
        ref_value[3] = s_tl[0];

        // Setup algebraic update for node: Pass update information
        dynamic_cast<AlgebraicNode*>(el_pt->node_pt(jnod))
          ->add_node_update_info(Central_box, // enumerated ID
                                 this, // mesh
                                 geom_object_pt, // vector of geom objects
                                 ref_value); // vector of ref. values
      }
    }

    // Element 1: Bottom right box
    //----------------------------
    {
      unsigned ielem = 1;
      FiniteElement* el_pt = Mesh::finite_element_pt(ielem);

      // Loop over all nodes in the element and give them the update
      // info appropriate for the current element

      // Double loop over nodes
      unsigned nnod_lin =
        dynamic_cast<ELEMENT*>(Mesh::finite_element_pt(ielem))->nnode_1d();
      for (unsigned i0 = 0; i0 < nnod_lin; i0++)
      {
        // Fraction in the s_0-direction
        double rho_0 = double(i0) / double(nnod_lin - 1);

        for (unsigned i1 = 0; i1 < nnod_lin; i1++)
        {
          // Fraction in the s_1-direction
          double rho_1 = double(i1) / double(nnod_lin - 1);

          // Node number
          unsigned jnod = i0 + i1 * nnod_lin;

          // The update function requires three geometric objects
          Vector<GeomObject*> geom_object_pt(3);

          // The update function requires five parameters:
          Vector<double> ref_value(5);

          // First reference value: fractional s0-position inside box
          ref_value[0] = rho_0;

          // Second reference value: fractional s1-position inside box
          ref_value[1] = rho_1;

          // Wall element at bottom right end of wall mesh:
          geom_object_pt[0] = obj_br_pt;

          // Local coordinate in this wall element. Note:
          // We'll have to recompute this reference
          // when the mesh is refined  as we might fall off the element
          // otherwise
          ref_value[2] = s_br[0];

          // Wall element at top left end of wall mesh:
          geom_object_pt[1] = obj_tl_pt;

          // Local coordinate in this wall element. Note:
          // We'll have to recompute this reference
          // when the mesh is refined  as we might fall off the element
          // otherwise
          ref_value[3] = s_tl[0];

          // Reference point on wall
          Vector<double> xi_wall(1);
          xi_wall[0] = QuarterCircleSectorMesh<ELEMENT>::Xi_lo +
                       rho_1 * QuarterCircleSectorMesh<ELEMENT>::Fract_mid *
                         (QuarterCircleSectorMesh<ELEMENT>::Xi_hi -
                          QuarterCircleSectorMesh<ELEMENT>::Xi_lo);

          // Identify wall element number and local coordinate of
          // reference point on wall
          GeomObject* obj_wall_pt;
          Vector<double> s_wall(1);
          this->Wall_pt->locate_zeta(xi_wall, obj_wall_pt, s_wall);

          // Wall element at that contians reference point:
          geom_object_pt[2] = obj_wall_pt;

          // Local coordinate in this wall element. Note:
          // We'll have to recompute this reference
          // when the mesh is refined  as we might fall off the element
          // otherwise
          ref_value[4] = s_wall[0];

          // Setup algebraic update for node: Pass update information
          dynamic_cast<AlgebraicNode*>(el_pt->node_pt(jnod))
            ->add_node_update_info(Lower_right_box, // enumerated ID
                                   this, // mesh
                                   geom_object_pt, // vector of geom objects
                                   ref_value); // vector of ref. vals
        }
      }
    }


    // Element 2: Top left box
    //---------------------------
    {
      unsigned ielem = 2;
      FiniteElement* el_pt = Mesh::finite_element_pt(ielem);

      // Double loop over nodes
      unsigned nnod_lin =
        dynamic_cast<ELEMENT*>(Mesh::finite_element_pt(ielem))->nnode_1d();

      for (unsigned i0 = 0; i0 < nnod_lin; i0++)
      {
        // Fraction in the s_0-direction
        double rho_0 = double(i0) / double(nnod_lin - 1);

        for (unsigned i1 = 0; i1 < nnod_lin; i1++)
        {
          // Fraction in the s_1-direction
          double rho_1 = double(i1) / double(nnod_lin - 1);

          // Node number
          unsigned jnod = i0 + i1 * nnod_lin;

          // The update function requires three geometric objects
          Vector<GeomObject*> geom_object_pt(3);

          // The update function requires five parameters:
          Vector<double> ref_value(5);

          // First reference value: fractional s0-position inside box
          ref_value[0] = rho_0;

          // Second  reference value: fractional s1-position inside box
          ref_value[1] = rho_1;

          // Wall element at bottom right end of wall mesh:
          geom_object_pt[0] = obj_br_pt;

          // Local coordinate in this wall element. Note:
          // We'll have to recompute this reference
          // when the mesh is refined  as we might fall off the element
          // otherwise
          ref_value[2] = s_br[0];

          // Wall element at top left end of wall mesh:
          geom_object_pt[1] = obj_tl_pt;

          // Local coordinate in this wall element. Note:
          // We'll have to recompute this reference
          // when the mesh is refined  as we might fall off the element
          // otherwise
          ref_value[3] = s_tl[0];

          // Reference point on wall
          Vector<double> xi_wall(1);
          xi_wall[0] = QuarterCircleSectorMesh<ELEMENT>::Xi_hi +
                       rho_0 *
                         (1.0 - QuarterCircleSectorMesh<ELEMENT>::Fract_mid) *
                         (QuarterCircleSectorMesh<ELEMENT>::Xi_lo -
                          QuarterCircleSectorMesh<ELEMENT>::Xi_hi);

          // Identify wall element number and local coordinate of
          // reference point on wall
          GeomObject* obj_wall_pt;
          Vector<double> s_wall(1);
          this->Wall_pt->locate_zeta(xi_wall, obj_wall_pt, s_wall);

          // Wall element at that contians reference point:
          geom_object_pt[2] = obj_wall_pt;

          // Local coordinate in this wall element. Note:
          // We'll have to recompute this reference
          // when the mesh is refined  as we might fall off the element
          // otherwise
          ref_value[4] = s_wall[0];

          // Setup algebraic update for node: Pass update information
          dynamic_cast<AlgebraicNode*>(el_pt->node_pt(jnod))
            ->add_node_update_info(Upper_left_box, // Enumerated ID
                                   this, // mesh
                                   geom_object_pt, // vector of geom objects
                                   ref_value); // vector of ref. vals
        }
      }
    }
  }


  //======================================================================
  /// Algebraic update function: Update in central box according
  /// to wall shape at time level t (t=0: present; t>0: previous)
  //======================================================================
  template<class ELEMENT>
  void AlgebraicRefineableQuarterCircleSectorMesh<
    ELEMENT>::node_update_in_central_box(const unsigned& t,
                                         AlgebraicNode*& node_pt)
  {
#ifdef PARANOID
    // We're updating the nodal positions (!) at time level t
    // and determine them by evaluating the wall GeomObject's
    // position at that gime level. I believe this only makes sense
    // if the t-th history value in the positional timestepper
    // actually represents previous values (rather than some
    // generalised quantity). Hence if this function is called with
    // t>nprev_values(), we issue a warning and terminate the execution.
    // It *might* be possible that the code still works correctly
    // even if this condition is violated (e.g. if the GeomObject's
    // position() function returns the appropriate "generalised"
    // position value that is required by the timestepping scheme but it's
    // probably worth flagging this up and forcing the user to manually switch
    // off this warning if he/she is 100% sure that this is kosher.
    if (t > node_pt->position_time_stepper_pt()->nprev_values())
    {
      std::string error_message =
        "Trying to  update the nodal position at a time level\n";
      error_message += "beyond the number of previous values in the nodes'\n";
      error_message += "position timestepper. This seems highly suspect!\n";
      error_message += "If you're sure the code behaves correctly\n";
      error_message += "in your application, remove this warning \n";
      error_message += "or recompile with PARNOID switched off.\n";

      std::string function_name =
        "AlgebraicRefineableQuarterCircleSectorMesh::";
      function_name += "node_update_in_central_box()",
        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Extract references for update in central box by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value(Central_box));

    // Extract geometric objects for update in central box by copy construction
    Vector<GeomObject*> geom_object_pt(
      node_pt->vector_geom_object_pt(Central_box));

    // First reference value: fractional x-position of node inside box
    double rho_x = ref_value[0];

    // Second reference value: fractional y-position of node inside box
    double rho_y = ref_value[1];


    // Wall position in bottom right corner:

    // Pointer to wall element:
    GeomObject* obj_br_pt = geom_object_pt[0];

    // Eulerian dimension
    unsigned n_dim = obj_br_pt->ndim();

    // Local coordinate:
    Vector<double> s_br(1);
    s_br[0] = ref_value[2];

    // Get wall position
    Vector<double> r_br(n_dim);
    obj_br_pt->position(t, s_br, r_br);

    // Wall position in top left corner:

    // Pointer to wall element:
    GeomObject* obj_tl_pt = geom_object_pt[1];

    // Local coordinate:
    Vector<double> s_tl(1);
    s_tl[0] = ref_value[3];

    Vector<double> r_tl(n_dim);
    obj_tl_pt->position(t, s_tl, r_tl);

    // Assign new nodal coordinate
    node_pt->x(t, 0) = r_br[0] * Lambda_x * rho_x;
    node_pt->x(t, 1) = r_tl[1] * Lambda_y * rho_y;
  }


  //====================================================================
  /// Algebraic update function: Update in lower right box according
  /// to wall shape at time level t (t=0: present; t>0: previous)
  //====================================================================
  template<class ELEMENT>
  void AlgebraicRefineableQuarterCircleSectorMesh<
    ELEMENT>::node_update_in_lower_right_box(const unsigned& t,
                                             AlgebraicNode*& node_pt)
  {
#ifdef PARANOID
    // We're updating the nodal positions (!) at time level t
    // and determine them by evaluating the wall GeomObject's
    // position at that gime level. I believe this only makes sense
    // if the t-th history value in the positional timestepper
    // actually represents previous values (rather than some
    // generalised quantity). Hence if this function is called with
    // t>nprev_values(), we issue a warning and terminate the execution.
    // It *might* be possible that the code still works correctly
    // even if this condition is violated (e.g. if the GeomObject's
    // position() function returns the appropriate "generalised"
    // position value that is required by the timestepping scheme but it's
    // probably worth flagging this up and forcing the user to manually switch
    // off this warning if he/she is 100% sure that this is kosher.
    if (t > node_pt->position_time_stepper_pt()->nprev_values())
    {
      std::string error_message =
        "Trying to update the nodal position at a time level";
      error_message += "beyond the number of previous values in the nodes'";
      error_message += "position timestepper. This seems highly suspect!";
      error_message += "If you're sure the code behaves correctly";
      error_message += "in your application, remove this warning ";
      error_message += "or recompile with PARNOID switched off.";

      std::string function_name =
        "AlgebraicRefineableQuarterCircleSectorMesh::";
      function_name += "node_update_in_lower_right_box()",
        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Extract references for update in central box by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value(Lower_right_box));

    // Extract geometric objects for update in central box by copy construction
    Vector<GeomObject*> geom_object_pt(
      node_pt->vector_geom_object_pt(Lower_right_box));

    // First reference value: fractional s0-position of node inside box
    double rho_0 = ref_value[0];

    // Second reference value: fractional s1-position of node inside box
    double rho_1 = ref_value[1];

    // Wall position in bottom right corner:

    // Pointer to wall element:
    GeomObject* obj_br_pt = geom_object_pt[0];

    // Eulerian dimension
    unsigned n_dim = obj_br_pt->ndim();

    // Local coordinate:
    Vector<double> s_br(1);
    s_br[0] = ref_value[2];

    // Get wall position
    Vector<double> r_br(n_dim);
    obj_br_pt->position(t, s_br, r_br);

    // Wall position in top left corner:

    // Pointer to wall element:
    GeomObject* obj_tl_pt = geom_object_pt[1];

    // Local coordinate:
    Vector<double> s_tl(1);
    s_tl[0] = ref_value[3];

    Vector<double> r_tl(n_dim);
    obj_tl_pt->position(t, s_tl, r_tl);

    // Position of the corner of the central box
    Vector<double> r_box(n_dim);
    r_box[0] = Lambda_x * r_br[0];
    r_box[1] = Lambda_y * r_tl[1];

    // Position Vector to left end of box
    Vector<double> r_left(n_dim);
    r_left[0] = Lambda_x * r_br[0] + rho_1 * (r_box[0] - Lambda_x * r_br[0]);
    r_left[1] = Lambda_x * r_br[1] + rho_1 * (r_box[1] - Lambda_x * r_br[1]);

    // Wall position

    // Pointer to wall element:
    GeomObject* obj_wall_pt = geom_object_pt[2];

    // Local coordinate:
    Vector<double> s_wall(1);
    s_wall[0] = ref_value[4];

    Vector<double> r_wall(n_dim);
    obj_wall_pt->position(t, s_wall, r_wall);

    // Assign new nodal coordinate
    node_pt->x(t, 0) = r_left[0] + rho_0 * (r_wall[0] - r_left[0]);
    node_pt->x(t, 1) = r_left[1] + rho_0 * (r_wall[1] - r_left[1]);
  }
  //====================================================================
  /// Algebraic update function: Update in upper left box according
  /// to wall shape at time level t (t=0: present; t>0: previous)
  //====================================================================
  template<class ELEMENT>
  void AlgebraicRefineableQuarterCircleSectorMesh<
    ELEMENT>::node_update_in_upper_left_box(const unsigned& t,
                                            AlgebraicNode*& node_pt)
  {
#ifdef PARANOID
    // We're updating the nodal positions (!) at time level t
    // and determine them by evaluating the wall GeomObject's
    // position at that gime level. I believe this only makes sense
    // if the t-th history value in the positional timestepper
    // actually represents previous values (rather than some
    // generalised quantity). Hence if this function is called with
    // t>nprev_values(), we issue a warning and terminate the execution.
    // It *might* be possible that the code still works correctly
    // even if this condition is violated (e.g. if the GeomObject's
    // position() function returns the appropriate "generalised"
    // position value that is required by the timestepping scheme but it's
    // probably worth flagging this up and forcing the user to manually switch
    // off this warning if he/she is 100% sure that this is kosher.
    if (t > node_pt->position_time_stepper_pt()->nprev_values())
    {
      std::string error_message =
        "Trying to update the nodal position at a time level";
      error_message += "beyond the number of previous values in the nodes'";
      error_message += "position timestepper. This seems highly suspect!";
      error_message += "If you're sure the code behaves correctly";
      error_message += "in your application, remove this warning ";
      error_message += "or recompile with PARNOID switched off.";

      std::string function_name =
        "AlgebraicRefineableQuarterCircleSectorMesh::";
      function_name += "node_update_in_upper_left_box()";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Extract references for update in central box by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value(Upper_left_box));

    // Extract geometric objects for update  in central box by copy construction
    Vector<GeomObject*> geom_object_pt(
      node_pt->vector_geom_object_pt(Upper_left_box));

    // First reference value: fractional s0-position of node inside box
    double rho_0 = ref_value[0];

    // Second  reference value: fractional s1-position of node inside box
    double rho_1 = ref_value[1];

    // Wall position in bottom right corner:

    // Pointer to wall element:
    GeomObject* obj_br_pt = geom_object_pt[0];

    // Eulerian dimension
    unsigned n_dim = obj_br_pt->ndim();

    // Local coordinate:
    Vector<double> s_br(1);
    s_br[0] = ref_value[2];

    // Get wall position
    Vector<double> r_br(n_dim);
    obj_br_pt->position(t, s_br, r_br);

    // Wall position in top left corner:

    // Pointer to wall element:
    GeomObject* obj_tl_pt = node_pt->geom_object_pt(1);

    // Local coordinate:
    Vector<double> s_tl(1);
    s_tl[0] = node_pt->ref_value(3);

    Vector<double> r_tl(n_dim);
    obj_tl_pt->position(t, s_tl, r_tl);

    // Position of the corner of the central box
    Vector<double> r_box(n_dim);
    r_box[0] = Lambda_x * r_br[0];
    r_box[1] = Lambda_y * r_tl[1];

    // Position Vector to top face of central box
    Vector<double> r_top(n_dim);
    r_top[0] = Lambda_y * r_tl[0] + rho_0 * (r_box[0] - Lambda_y * r_tl[0]);
    r_top[1] = Lambda_y * r_tl[1] + rho_0 * (r_box[1] - Lambda_y * r_tl[1]);

    // Wall position

    // Pointer to wall element:
    GeomObject* obj_wall_pt = node_pt->geom_object_pt(2);

    // Local coordinate:
    Vector<double> s_wall(1);
    s_wall[0] = node_pt->ref_value(4);

    Vector<double> r_wall(n_dim);
    obj_wall_pt->position(t, s_wall, r_wall);


    // Assign new nodal coordinate
    node_pt->x(t, 0) = r_top[0] + rho_1 * (r_wall[0] - r_top[0]);
    node_pt->x(t, 1) = r_top[1] + rho_1 * (r_wall[1] - r_top[1]);
  }


  //======================================================================
  /// Update algebraic update function for nodes in lower right box.
  //======================================================================
  template<class ELEMENT>
  void AlgebraicRefineableQuarterCircleSectorMesh<
    ELEMENT>::update_node_update_in_lower_right_box(AlgebraicNode*& node_pt)
  {
    // Extract references for update in central box by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value(Lower_right_box));

    // Extract geometric objects for updatein central box by copy construction
    Vector<GeomObject*> geom_object_pt(
      node_pt->vector_geom_object_pt(Lower_right_box));


    // Now remove the update  info to allow overwriting below
    node_pt->kill_node_update_info(Lower_right_box);

    // Fixed reference value: Start coordinate on wall
    double xi_lo = QuarterCircleSectorMesh<ELEMENT>::Xi_lo;

    // Fixed reference value: Fractional position of dividing line
    double fract_mid = QuarterCircleSectorMesh<ELEMENT>::Fract_mid;

    // Fixed reference value: End coordinate on wall
    double xi_hi = QuarterCircleSectorMesh<ELEMENT>::Xi_hi;

    // Second reference value: fractional s1-position of node inside box
    double rho_1 = ref_value[1];


    // Update reference to wall point in bottom right corner
    //------------------------------------------------------

    // Wall position in bottom right corner:
    Vector<double> xi(1);
    xi[0] = xi_lo;

    // Find corresponding wall element/local coordinate
    GeomObject* obj_br_pt;
    Vector<double> s_br(1);
    this->Wall_pt->locate_zeta(xi, obj_br_pt, s_br);

    // Wall element at bottom right end of wall mesh:
    geom_object_pt[0] = obj_br_pt;

    // Local coordinate in this wall element. Note:
    // We'll have to recompute this reference
    // when the mesh is refined  as we might fall off the element otherwise
    ref_value[2] = s_br[0];


    // Update reference to wall point in upper left corner
    //----------------------------------------------------

    // Wall position in top left corner
    xi[0] = xi_hi;

    // Find corresponding wall element/local coordinate
    GeomObject* obj_tl_pt;
    Vector<double> s_tl(1);
    this->Wall_pt->locate_zeta(xi, obj_tl_pt, s_tl);

    // Wall element at top left end of wall mesh:
    geom_object_pt[1] = obj_tl_pt;

    // Local coordinate in this wall element. Note:
    // We'll have to recompute this reference
    // when the mesh is refined  as we might fall off the element otherwise
    ref_value[3] = s_tl[0];


    // Update reference to reference point on wall
    //--------------------------------------------

    // Reference point on wall
    Vector<double> xi_wall(1);
    xi_wall[0] = xi_lo + rho_1 * fract_mid * (xi_hi - xi_lo);

    // Identify wall element number and local coordinate of
    // reference point on wall
    GeomObject* obj_wall_pt;
    Vector<double> s_wall(1);
    this->Wall_pt->locate_zeta(xi_wall, obj_wall_pt, s_wall);

    // Wall element at that contians reference point:
    geom_object_pt[2] = obj_wall_pt;

    // Local coordinate in this wall element. Note:
    // We'll have to recompute this reference
    // when the mesh is refined  as we might fall off the element otherwise
    ref_value[4] = s_wall[0];

    // Setup algebraic update for node: Pass update information
    node_pt->add_node_update_info(Lower_right_box, // Enumerated ID
                                  this, // mesh
                                  geom_object_pt, // vector of geom objects
                                  ref_value); // vector of ref. vals
  }

  //======================================================================
  /// Update algebraic update function for nodes in upper left box
  //======================================================================
  template<class ELEMENT>
  void AlgebraicRefineableQuarterCircleSectorMesh<
    ELEMENT>::update_node_update_in_upper_left_box(AlgebraicNode*& node_pt)
  {
    // Extract references for update in central box by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value(Upper_left_box));

    // Extract geometric objects for update in central box by copy construction
    Vector<GeomObject*> geom_object_pt(
      node_pt->vector_geom_object_pt(Upper_left_box));

    // Now remove the update info to allow overwriting below
    node_pt->kill_node_update_info(Upper_left_box);

    // Fixed reference value: Start coordinate on wall
    double xi_lo = QuarterCircleSectorMesh<ELEMENT>::Xi_lo;

    // Fixed reference value:  Fractional position of dividing line
    double fract_mid = QuarterCircleSectorMesh<ELEMENT>::Fract_mid;

    // Fixed reference value: End coordinate on wall
    double xi_hi = QuarterCircleSectorMesh<ELEMENT>::Xi_hi;

    // First reference value: fractional s0-position of node inside box
    double rho_0 = ref_value[0];

    // Update reference to wall point in bottom right corner
    //------------------------------------------------------

    // Wall position in bottom right corner:
    Vector<double> xi(1);
    xi[0] = xi_lo;

    // Find corresponding wall element/local coordinate
    GeomObject* obj_br_pt;
    Vector<double> s_br(1);
    this->Wall_pt->locate_zeta(xi, obj_br_pt, s_br);

    // Wall element at bottom right end of wall mesh:
    geom_object_pt[0] = obj_br_pt;

    // Local coordinate in this wall element. Note:
    // We'll have to recompute this reference
    // when the mesh is refined  as we might fall off the element otherwise
    ref_value[2] = s_br[0];

    // Update reference to wall point in upper left corner
    //----------------------------------------------------

    // Wall position in top left corner
    xi[0] = xi_hi;

    // Find corresponding wall element/local coordinate
    GeomObject* obj_tl_pt;
    Vector<double> s_tl(1);
    this->Wall_pt->locate_zeta(xi, obj_tl_pt, s_tl);

    // Wall element at top left end of wall mesh:
    geom_object_pt[1] = obj_tl_pt;

    // Local coordinate in this wall element. Note:
    // We'll have to recompute this reference
    // when the mesh is refined  as we might fall off the element otherwise
    ref_value[3] = s_tl[0];


    // Update reference to reference point on wall
    //--------------------------------------------

    // Reference point on wall
    Vector<double> xi_wall(1);
    xi_wall[0] = xi_hi + rho_0 * (1.0 - fract_mid) * (xi_lo - xi_hi);

    // Identify reference point on wall
    GeomObject* obj_wall_pt;
    Vector<double> s_wall(1);
    this->Wall_pt->locate_zeta(xi_wall, obj_wall_pt, s_wall);

    // Wall element at that contians reference point:
    geom_object_pt[2] = obj_wall_pt;

    // Local coordinate in this wall element. Note:
    // We'll have to recompute this reference
    // when the mesh is refined  as we might fall off the element otherwise
    ref_value[4] = s_wall[0];

    // Setup algebraic update for node: Pass update information
    node_pt->add_node_update_info(Upper_left_box, // Enumerated ID
                                  this, // mesh
                                  geom_object_pt, // vector of geom objects
                                  ref_value); // vector of ref. vals
  }


} // namespace oomph
#endif
