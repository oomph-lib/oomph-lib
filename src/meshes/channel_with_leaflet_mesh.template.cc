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

#ifndef OOMPH_CHANNEL_WITH_LEAFLET_MESH_TEMPLATE_CC
#define OOMPH_CHANNEL_WITH_LEAFLET_MESH_TEMPLATE_CC


// Include the headers file
#include "channel_with_leaflet_mesh.template.h"


namespace oomph
{
  //===================================================================
  /// Constructor: Pass pointer to GeomObject that represents the leaflet,
  /// the length of the domain to left and right of the leaflet, the
  /// height of the leaflet and the overall height of the channel,
  /// the number of element columns to the left and right of the leaflet,
  /// the number of rows of elements from the bottom of the channel to
  /// the end of the leaflet, the number of rows of elements above the
  /// end of the leaflet. Timestepper defaults to Steady default
  /// Timestepper defined in the Mesh base class
  //===================================================================
  template<class ELEMENT>
  ChannelWithLeafletMesh<ELEMENT>::ChannelWithLeafletMesh(
    GeomObject* leaflet_pt,
    const double& lleft,
    const double& lright,
    const double& hleaflet,
    const double& htot,
    const unsigned& nleft,
    const unsigned& nright,
    const unsigned& ny1,
    const unsigned& ny2,
    TimeStepper* time_stepper_pt)
    : SimpleRectangularQuadMesh<ELEMENT>(
        nright + nleft, ny1 + ny2, lright + lleft, htot, time_stepper_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Copy pointer to leaflet into private storage
    Leaflet_pt = leaflet_pt;

    // Create the ChannelWithLeafletDomain with the wall
    // represented by the geometric object
    Domain_pt = new ChannelWithLeafletDomain(
      leaflet_pt, lleft, lright, hleaflet, htot, nleft, nright, ny1, ny2);


    // Total number of (macro/finite)elements
    unsigned nmacro = (ny1 + ny2) * (nleft + nright);

    // Loop over all elements and set macro element pointer
    for (unsigned e = 0; e < nmacro; e++)
    {
      // Get pointer to finite element
      FiniteElement* el_pt = this->finite_element_pt(e);

      // Set pointer to macro element
      el_pt->set_macro_elem_pt(this->Domain_pt->macro_element_pt(e));
    }

    // Update the nodal positions corresponding to their
    // macro element representations.
    this->node_update();

    // Change the numbers of boundaries
    this->set_nboundary(7);

    // Remove the nodes of boundary 0
    this->remove_boundary_nodes(0);

    // Get the number of nodes along the element edge from first element
    unsigned nnode_1d = this->finite_element_pt(0)->nnode_1d();

    // Boundary 0 will be before the wall, boundary 6 after it
    for (unsigned e = 0; e < nleft; e++)
    {
      for (unsigned i = 0; i < nnode_1d; i++)
      {
        Node* nod_pt = this->finite_element_pt(e)->node_pt(i);
        // do not add the last node : it will go to boundary 6
        if (e != nleft - 1 || i != 2)
        {
          this->add_boundary_node(0, nod_pt);
        }
      }
    }

    for (unsigned e = nleft; e < nleft + nright; e++)
    {
      for (unsigned i = 0; i < nnode_1d; i++)
      {
        Node* nod_pt = this->finite_element_pt(e)->node_pt(i);
        this->add_boundary_node(6, nod_pt);
      }
    }

    // Vector of Lagrangian coordinates used as boundary coordinate
    Vector<double> zeta(1);

    // Set the wall as a boundary
    for (unsigned k = 0; k < ny1; k++)
    {
      unsigned e = (nleft + nright) * k + nleft - 1;
      for (unsigned i = 0; i < nnode_1d; i++)
      {
        Node* nod_pt =
          this->finite_element_pt(e)->node_pt((i + 1) * nnode_1d - 1);
        this->convert_to_boundary_node(nod_pt);
        this->add_boundary_node(4, nod_pt);

        // Set coordinates
        zeta[0] = double(k) * hleaflet / double(ny1) +
                  double(i) * hleaflet / double(ny1) / double(nnode_1d - 1);
        nod_pt->set_coordinates_on_boundary(4, zeta);
      }
    }


    // Duplicate the nodes of the wall and assign then as a boundary.
    // This will make one boundary for the east of the elements at the
    // left of the wall, and one for the west of the elements at the right
    // of the wall.
    // This is necessary to use TaylorHoodElements, because it will allow
    // a discontinuity of the pressure accross the wall.
    // We separate the lower element from the rest as we do not want to
    // add the same node twice to the boundary, and the upper element as its
    // upper node must be the same one than the node of the upper element
    // at the right of the wall

    // Lower element
    unsigned e = nleft - 1;

    // Allocate storage for newly created node outside
    // so we can refer to the most recently created one below.
    Node* nod_pt = 0;

    // duplicate the east nodes and add them to the 6th boundary
    // Add the first node to the 0th boundary (horizontal)
    for (unsigned i = 0; i < nnode_1d; i++)
    {
      nod_pt = this->finite_element_pt(e)->construct_boundary_node(
        (i + 1) * nnode_1d - 1, time_stepper_pt);
      this->add_boundary_node(5, nod_pt);
      if (i == 0)
      {
        this->add_boundary_node(0, nod_pt);
      }
      this->add_node_pt(nod_pt);
      // Set coordinates
      zeta[0] = i * hleaflet / double(ny1) / double(nnode_1d - 1);
      nod_pt->set_coordinates_on_boundary(5, zeta);
    }


    // Other elements just at the left of the wall
    for (unsigned k = 1; k < ny1; k++)
    {
      e = (nleft + nright) * k + nleft - 1;

      // add the upper node of the previous element
      this->finite_element_pt(e)->node_pt(nnode_1d - 1) = nod_pt;
      this->add_boundary_node(5, nod_pt);
      // Set coordinates
      zeta[0] = k * hleaflet / double(ny1);
      nod_pt->set_coordinates_on_boundary(5, zeta);

      // Loop over other nodes on element's eastern edge
      for (unsigned i = 1; i < nnode_1d; i++)
      {
        // Don't duplicate the node at the very top of the "obstacle"
        if ((k == ny1 - 1) && (i == nnode_1d - 1))
        {
          // Get the node but don't do anything else
          nod_pt = this->finite_element_pt(e)->node_pt(nnode_1d * nnode_1d - 1);
        }
        else
        {
          // Overwrite the node with a boundary node
          nod_pt = this->finite_element_pt(e)->construct_boundary_node(
            (i + 1) * nnode_1d - 1, time_stepper_pt);
          this->add_node_pt(nod_pt);
        }

        // Add node to boundary
        this->add_boundary_node(5, nod_pt);
        // Set coordinates
        zeta[0] = double(k) * hleaflet / double(ny1) +
                  double(i) * hleaflet / double(ny1) / double(nnode_1d - 1);
        nod_pt->set_coordinates_on_boundary(5, zeta);
      }
    }

    this->node_update();

    // Re-setup lookup scheme that establishes which elements are located
    // on the mesh boundaries (doesn't need to be wiped)
    this->setup_boundary_element_info();

    // We have parametrised boundary 4 and 5
    this->Boundary_coordinate_exists[4] = true;
    this->Boundary_coordinate_exists[5] = true;

  } // end of constructor


  /// ////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Setup algebraic node update.  Leaflet is
  /// assumed to be in its undeformed (straight vertical) position!
  //=================================================================
  template<class ELEMENT>
  void AlgebraicChannelWithLeafletMesh<ELEMENT>::setup_algebraic_node_update()
  {
    // Extract some reference lengths from the Domain.
    double hleaflet = this->domain_pt()->hleaflet();
    double htot = this->domain_pt()->htot();
    double lleft = this->domain_pt()->lleft();
    double lright = this->domain_pt()->lright();

    // Loop over all nodes in mesh
    unsigned nnod = this->nnode();
    for (unsigned j = 0; j < nnod; j++)
    {
      // Get pointer to node
      AlgebraicNode* nod_pt = node_pt(j);

      // Get coordinates
      double x = nod_pt->x(0);
      double y = nod_pt->x(1);

      // Quick check to know if the wall is in its undeformed position
      // It actually checks that the top of the wall is in (x_0,hleaflet)
      Vector<double> zeta(1);
      Vector<double> r(2);
      zeta[0] = Hleaflet;
      this->Leaflet_pt->position(zeta, r);
      if ((r[0] != X_0) || (r[1] != hleaflet))
      {
        std::ostringstream error_stream;
        error_stream << "Wall must be in its undeformed position when\n"
                     << "algebraic node update information is set up!\n "
                     << r[0] << " " << X_0 << " " << r[1] << " " << hleaflet
                     << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }


      // The update function requires four parameters:
      Vector<double> ref_value(4);

      // Part I
      if ((x <= X_0) && (y <= hleaflet))
      {
        // First reference value: y0
        ref_value[0] = y;

        // zeta coordinate on wall : fourth reference value
        //(needed for adaptativity)
        Vector<double> zeta(1);
        zeta[0] = y;
        ref_value[3] = zeta[0];

        // Sub-geomobject corresponding to the zeta coordinate on the wall
        GeomObject* geom_obj_pt;
        Vector<double> s(1);
        this->Leaflet_pt->locate_zeta(zeta, geom_obj_pt, s);

        // Create vector of geomobject for add_node_update_info()
        Vector<GeomObject*> geom_object_pt(1);
        geom_object_pt[0] = geom_obj_pt;

        // Second  reference value: Reference local coordinate
        // in the sub-geomobject
        ref_value[1] = s[0];

        // Third reference value:fraction of the horizontal line
        // between the edge and the wall
        ref_value[2] = (lleft + x - X_0) / lleft;


        // Setup algebraic update for node: Pass update information
        // to AlgebraicNode:
        nod_pt->add_node_update_info(1, // ID
                                     this, // mesh
                                     geom_object_pt, // vector of geom objects
                                     ref_value); // vector of ref. values
      }
      // Part II
      if ((x >= X_0) && (y <= hleaflet))
      {
        // First reference value: y0
        ref_value[0] = y;

        // zeta coordinate on wall: fourth reference value
        //(needed for adaptativity)
        Vector<double> zeta(1);
        zeta[0] = y;
        ref_value[3] = zeta[0];

        // Sub-geomobject corresponding to the zeta coordinate on the wall
        GeomObject* geom_obj_pt;
        Vector<double> s(1);
        this->Leaflet_pt->locate_zeta(zeta, geom_obj_pt, s);

        // Create vector of geomobject for add_node_update_info()
        Vector<GeomObject*> geom_object_pt(1);
        geom_object_pt[0] = geom_obj_pt;

        // Second  reference value: Reference local coordinate
        // in the sub-geomobject
        ref_value[1] = s[0];

        // Third reference value:fraction of the horizontal line
        // between the edge and the wall
        ref_value[2] = (x - X_0) / lright;

        // Setup algebraic update for node: Pass update information
        // to AlgebraicNode:
        nod_pt->add_node_update_info(2, // ID
                                     this, // mesh
                                     geom_object_pt, // vector of geom objects
                                     ref_value); // vector of ref. values
      }
      // Part III
      if ((x <= X_0) && (y >= hleaflet))
      {
        // First reference value: y0
        ref_value[0] = y;

        // Second reference value: zeta coordinate on the middle line
        ref_value[1] = (y - hleaflet) / (htot - hleaflet);

        // Third reference value:fraction of the horizontal line
        // between the edge and the middle line
        ref_value[2] = (lleft + x - X_0) / lleft;

        // geomobject
        Vector<GeomObject*> geom_object_pt(1);
        geom_object_pt[0] = this->Leaflet_pt;

        // Setup algebraic update for node: Pass update information
        // to AlgebraicNode:
        nod_pt->add_node_update_info(3, // ID
                                     this, // mesh
                                     geom_object_pt, // vector of geom objects
                                     ref_value); // vector of ref. values
      }
      // Part IV
      if ((x >= X_0) && (y >= hleaflet))
      {
        // First reference value: y0
        ref_value[0] = y;

        // Second reference value: zeta coordinate on wall
        ref_value[1] = (y - hleaflet) / (htot - hleaflet);

        // Third reference value:fraction of the horizontal line
        // between the edge and the wall
        ref_value[2] = (x - X_0) / lright;

        // geomobject
        Vector<GeomObject*> geom_object_pt(1);
        geom_object_pt[0] = this->Leaflet_pt;

        // Setup algebraic update for node: Pass update information
        // to AlgebraicNode:
        nod_pt->add_node_update_info(4, // ID
                                     this, // mesh
                                     geom_object_pt, // vector of geom objects
                                     ref_value); // vector of ref. values
      }
    }

  } // end of setup_algebraic_node_update


  //=================================================================
  /// Perform algebraic node update
  //=================================================================
  template<class ELEMENT>
  void AlgebraicChannelWithLeafletMesh<ELEMENT>::algebraic_node_update(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    unsigned id = node_pt->node_update_fct_id();

    switch (id)
    {
      case 1:
        node_update_I(t, node_pt);
        break;

      case 2:
        node_update_II(t, node_pt);
        break;

      case 3:
        node_update_III(t, node_pt);
        break;

      case 4:
        node_update_IV(t, node_pt);
        break;

      default:
        std::ostringstream error_message;
        error_message << "The node update fct id is " << id
                      << ", but it should only be one of " << 1 << ", " << 2
                      << ", " << 3 << " or " << 4 << std::endl;
        std::string function_name =
          "AlgebraicChannelWithLeafletMesh::algebraic_node_update()";

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }

  } // end of algebraic_node_update()


  //=================================================================
  /// Node update for region I
  //=================================================================
  template<class ELEMENT>
  void AlgebraicChannelWithLeafletMesh<ELEMENT>::node_update_I(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // relevant data of the domain for part I
    double lleft = this->domain_pt()->lleft();

    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to wall geom object
    GeomObject* leaflet_pt = geom_object_pt[0];

    // Coordinates of the steady node on the left boundary
    // corresponding to the current node
    double y0 = ref_value[0];
    double x0 = -lleft + X_0;

    // second reference value: local coordinate on wall
    Vector<double> s(1);
    s[0] = ref_value[1];


    // Get position vector to wall at timestep t
    Vector<double> r_wall(2);
    leaflet_pt->position(t, s, r_wall);


    // Third reference value : fraction of the horizontal line
    // between the edge and the wall
    double r = ref_value[2];


    // Assign new nodal coordinates
    node_pt->x(t, 0) = x0 + r * (r_wall[0] - x0);
    node_pt->x(t, 1) = y0 + r * (r_wall[1] - y0);
  }


  //=================================================================
  /// Node update for region II
  //=================================================================
  template<class ELEMENT>
  void AlgebraicChannelWithLeafletMesh<ELEMENT>::node_update_II(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // relevant data of the domain for part II
    double lright = this->domain_pt()->lright();

    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to wall geom object
    GeomObject* leaflet_pt = geom_object_pt[0];

    // Coordinates of the steady node on the right boundary
    // corresponding to the current node
    double y0 = ref_value[0];
    double x0 = X_0 + lright;

    // Second reference value: Zeta coordinate on wall
    Vector<double> s(1);
    s[0] = ref_value[1];

    // Get position vector to wall at timestep t
    Vector<double> r_wall(2);
    leaflet_pt->position(t, s, r_wall);

    // Third reference value : fraction of the horizontal line
    // between the wall and the right edge
    double r = ref_value[2];

    // Assign new nodal coordinates
    node_pt->x(t, 0) = r_wall[0] + r * (x0 - r_wall[0]);
    node_pt->x(t, 1) = r_wall[1] + r * (y0 - r_wall[1]);
  }


  //=================================================================
  /// Slanted bound : helper function
  //=================================================================
  template<class ELEMENT>
  void AlgebraicChannelWithLeafletMesh<ELEMENT>::slanted_bound_up(
    const unsigned& t, const Vector<double>& zeta, Vector<double>& r)
  {
    /// Coordinates of the point on the boundary beetween the upper
    /// and the lower part, in the same column, at the east.
    double htot = this->domain_pt()->htot();

    Vector<double> xi(1);
    xi[0] = Hleaflet;

    Vector<double> r_join(2);

    this->Leaflet_pt->position(t, xi, r_join);

    r[0] = r_join[0] + zeta[0] * (X_0 - r_join[0]);
    r[1] = r_join[1] + zeta[0] * (htot - r_join[1]);
  }

  //=================================================================
  /// Node update for region III
  //=================================================================
  template<class ELEMENT>
  void AlgebraicChannelWithLeafletMesh<ELEMENT>::node_update_III(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // relevant data of the domain for part I
    double lleft = this->domain_pt()->lleft();

    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Coordinates of the steady node on the left boundary
    // corresponding to the current node
    double y0 = ref_value[0];
    double x0 = -lleft + X_0;

    Vector<double> s(1);
    s[0] = ref_value[1];

    // Get position vector
    Vector<double> r_line(2);
    slanted_bound_up(t, s, r_line);

    // Third reference value : fraction of the horizontal line
    // between the edge and the middle line
    double r = ref_value[2];

    // Assign new nodal coordinates
    node_pt->x(t, 0) = x0 + r * (r_line[0] - x0);
    node_pt->x(t, 1) = y0 + r * (r_line[1] - y0);
  }

  //=================================================================
  /// Node update for region IV
  //=================================================================
  template<class ELEMENT>
  void AlgebraicChannelWithLeafletMesh<ELEMENT>::node_update_IV(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // relevant data of the domain for part I
    double lright = this->domain_pt()->lright();

    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Coordinates of the steady node on the left boundary
    // corresponding to the current node
    double y0 = ref_value[0];
    double x0 = X_0 + lright;

    // Second reference value: Zeta coordinate on the middle line
    Vector<double> s(1);
    s[0] = ref_value[1];

    // Get position vector  at timestep t
    Vector<double> r_line(2);
    slanted_bound_up(t, s, r_line);

    // Third reference value : fraction of the horizontal line
    // between the middle line and the right edge
    double r = ref_value[2];

    // Assign new nodal coordinates
    node_pt->x(t, 0) = r_line[0] + r * (x0 - r_line[0]);
    node_pt->x(t, 1) = r_line[1] + r * (y0 - r_line[1]);
  }


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


  //===================================================================
  ///  Update the node update functions
  //===================================================================
  template<class ELEMENT>
  void RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>::update_node_update(
    AlgebraicNode*& node_pt)
  {
    // Extract ID
    unsigned id = node_pt->node_update_fct_id();

    if ((id == 1) || (id == 2))
    {
      // Extract reference values for node update by copy construction
      Vector<double> ref_value(node_pt->vector_ref_value());

      // Get zeta coordinate on wall
      Vector<double> zeta_wall(1);
      zeta_wall[0] = ref_value[3];

      // Get the sub-geomobject and the local coordinate
      Vector<double> s(1);
      GeomObject* geom_obj_pt;
      this->Leaflet_pt->locate_zeta(zeta_wall, geom_obj_pt, s);

      // Extract geometric objects for update by copy construction
      Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

      // Update the pointer to the (sub-)GeomObject within which the
      // reference point is located. (If the wall is simple GeomObject
      // this is the same as Leaflet_pt; if it's a compound GeomObject
      // this points to the sub-object)
      geom_object_pt[0] = geom_obj_pt;

      // Update second reference value: Reference local coordinate
      // in wall sub-element
      ref_value[1] = s[0];


      if (id == 1)
      {
        // Kill the existing node update info
        node_pt->kill_node_update_info(1);

        // Setup algebraic update for node: Pass update information
        node_pt->add_node_update_info(1, // id
                                      this, // mesh
                                      geom_object_pt, // vector of geom objects
                                      ref_value); // vector of ref. values
      }
      else if (id == 2)
      {
        // Kill the existing node update info
        node_pt->kill_node_update_info(2);

        // Setup algebraic update for node: Pass update information
        node_pt->add_node_update_info(2, // id
                                      this, // mesh
                                      geom_object_pt, // vector of geom objects
                                      ref_value); // vector of ref. values
      }
    }
  }

} // namespace oomph

#endif
