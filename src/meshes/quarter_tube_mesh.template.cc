// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_QUARTER_TUBE_MESH_TEMPLATE_CC
#define OOMPH_QUARTER_TUBE_MESH_TEMPLATE_CC

#include "quarter_tube_mesh.template.h"


namespace oomph
{
  //====================================================================
  /// Constructor for deformable quarter tube mesh class.
  /// The domain is specified by the GeomObject that
  /// identifies boundary 3. Pass pointer to geometric object that
  /// specifies the wall, start and end coordinates on the
  /// geometric object, and the fraction along
  /// which the dividing line is to be placed, and the timestepper.
  /// Timestepper defaults to Static dummy timestepper.
  //====================================================================
  template<class ELEMENT>
  QuarterTubeMesh<ELEMENT>::QuarterTubeMesh(GeomObject* wall_pt,
                                            const Vector<double>& xi_lo,
                                            const double& fract_mid,
                                            const Vector<double>& xi_hi,
                                            const unsigned& nlayer,
                                            TimeStepper* time_stepper_pt)
    : Wall_pt(wall_pt), Xi_lo(xi_lo), Fract_mid(fract_mid), Xi_hi(xi_hi)
  {
    // Mesh can only be built with 3D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

    // Build macro element-based domain
    Domain_pt = new QuarterTubeDomain(wall_pt, xi_lo, fract_mid, xi_hi, nlayer);

    // Set the number of boundaries
    set_nboundary(5);

    // We have only bothered to parametrise boundary 3
    Boundary_coordinate_exists[3] = true;

    // Allocate the store for the elements
    unsigned nelem = 3 * nlayer;
    Element_pt.resize(nelem);

    // Create  dummy element so we can determine the number of nodes
    ELEMENT* dummy_el_pt = new ELEMENT;

    // Read out the number of linear points in the element
    unsigned n_p = dummy_el_pt->nnode_1d();

    // Kill the element
    delete dummy_el_pt;

    // Can now allocate the store for the nodes
    unsigned nnodes_total =
      (n_p * n_p + (n_p - 1) * n_p + (n_p - 1) * (n_p - 1)) *
      (1 + nlayer * (n_p - 1));
    Node_pt.resize(nnodes_total);


    Vector<double> s(3);
    Vector<double> r(3);

    // Storage for the intrinsic boundary coordinate
    Vector<double> zeta(2);

    // Loop over elements and create all nodes
    for (unsigned ielem = 0; ielem < nelem; ielem++)
    {
      // Create element
      Element_pt[ielem] = new ELEMENT;

      // Loop over rows in z/s_2-direction
      for (unsigned i2 = 0; i2 < n_p; i2++)
      {
        // Loop over rows in y/s_1-direction
        for (unsigned i1 = 0; i1 < n_p; i1++)
        {
          // Loop over rows in x/s_0-direction
          for (unsigned i0 = 0; i0 < n_p; i0++)
          {
            // Local node number
            unsigned jnod_local = i0 + i1 * n_p + i2 * n_p * n_p;

            // Create the node
            Node* node_pt = finite_element_pt(ielem)->construct_node(
              jnod_local, time_stepper_pt);

            // Set the position of the node from macro element mapping
            s[0] = -1.0 + 2.0 * double(i0) / double(n_p - 1);
            s[1] = -1.0 + 2.0 * double(i1) / double(n_p - 1);
            s[2] = -1.0 + 2.0 * double(i2) / double(n_p - 1);
            Domain_pt->macro_element_pt(ielem)->macro_map(s, r);

            node_pt->x(0) = r[0];
            node_pt->x(1) = r[1];
            node_pt->x(2) = r[2];
          }
        }
      }
    }

    // Initialise number of global nodes
    unsigned node_count = 0;

    // Tolerance for node killing:
    double node_kill_tol = 1.0e-12;

    // Check for error in node killing
    bool stopit = false;

    // Loop over elements
    for (unsigned ielem = 0; ielem < nelem; ielem++)
    {
      // Which layer?
      unsigned ilayer = unsigned(ielem / 3);

      // Which macro element?
      switch (ielem % 3)
      {
          // Macro element 0: Central box
          //-----------------------------
        case 0:


          // Loop over rows in z/s_2-direction
          for (unsigned i2 = 0; i2 < n_p; i2++)
          {
            // Loop over rows in y/s_1-direction
            for (unsigned i1 = 0; i1 < n_p; i1++)
            {
              // Loop over rows in x/s_0-direction
              for (unsigned i0 = 0; i0 < n_p; i0++)
              {
                // Local node number
                unsigned jnod_local = i0 + i1 * n_p + i2 * n_p * n_p;

                // Has the node been killed?
                bool killed = false;

                // First layer of all nodes in s_2 direction gets killed
                // and re-directed to nodes in previous element layer
                if ((i2 == 0) && (ilayer > 0))
                {
                  // Neighbour element
                  unsigned ielem_neigh = ielem - 3;

                  // Node in neighbour element
                  unsigned i0_neigh = i0;
                  unsigned i1_neigh = i1;
                  unsigned i2_neigh = n_p - 1;
                  unsigned jnod_local_neigh =
                    i0_neigh + i1_neigh * n_p + i2_neigh * n_p * n_p;

                  // Check:
                  for (unsigned i = 0; i < 3; i++)
                  {
                    double error = std::fabs(
                      finite_element_pt(ielem)->node_pt(jnod_local)->x(i) -
                      finite_element_pt(ielem_neigh)
                        ->node_pt(jnod_local_neigh)
                        ->x(i));
                    if (error > node_kill_tol)
                    {
                      oomph_info << "Error in node killing for i " << i << " "
                                 << error << std::endl;
                      stopit = true;
                    }
                  }

                  // Kill node
                  delete finite_element_pt(ielem)->node_pt(jnod_local);
                  killed = true;

                  // Set pointer to neighbour:
                  finite_element_pt(ielem)->node_pt(jnod_local) =
                    finite_element_pt(ielem_neigh)->node_pt(jnod_local_neigh);
                }

                // No duplicate node: Copy across to mesh
                if (!killed)
                {
                  Node_pt[node_count] =
                    finite_element_pt(ielem)->node_pt(jnod_local);

                  // Set boundaries:

                  // Back: Boundary 0
                  if ((i2 == 0) && (ilayer == 0))
                  {
                    this->convert_to_boundary_node(Node_pt[node_count]);
                    add_boundary_node(0, Node_pt[node_count]);
                  }

                  // Front: Boundary 4
                  if ((i2 == n_p - 1) && (ilayer == nlayer - 1))
                  {
                    this->convert_to_boundary_node(Node_pt[node_count]);
                    add_boundary_node(4, Node_pt[node_count]);
                  }

                  // Left symmetry plane: Boundary 1
                  if (i0 == 0)
                  {
                    this->convert_to_boundary_node(Node_pt[node_count]);
                    add_boundary_node(1, Node_pt[node_count]);
                  }

                  // Bottom symmetry plane: Boundary 2
                  if (i1 == 0)
                  {
                    this->convert_to_boundary_node(Node_pt[node_count]);
                    add_boundary_node(2, Node_pt[node_count]);
                  }

                  // Increment node counter
                  node_count++;
                }
              }
            }
          }


          break;

          // Macro element 1: Lower right box
          //---------------------------------
        case 1:


          // Loop over rows in z/s_2-direction
          for (unsigned i2 = 0; i2 < n_p; i2++)
          {
            // Loop over rows in y/s_1-direction
            for (unsigned i1 = 0; i1 < n_p; i1++)
            {
              // Loop over rows in x/s_0-direction
              for (unsigned i0 = 0; i0 < n_p; i0++)
              {
                // Local node number
                unsigned jnod_local = i0 + i1 * n_p + i2 * n_p * n_p;

                // Has the node been killed?
                bool killed = false;

                // First layer of all nodes in s_2 direction gets killed
                // and re-directed to nodes in previous element layer
                if ((i2 == 0) && (ilayer > 0))
                {
                  // Neighbour element
                  unsigned ielem_neigh = ielem - 3;

                  // Node in neighbour element
                  unsigned i0_neigh = i0;
                  unsigned i1_neigh = i1;
                  unsigned i2_neigh = n_p - 1;
                  unsigned jnod_local_neigh =
                    i0_neigh + i1_neigh * n_p + i2_neigh * n_p * n_p;


                  // Check:
                  for (unsigned i = 0; i < 3; i++)
                  {
                    double error = std::fabs(
                      finite_element_pt(ielem)->node_pt(jnod_local)->x(i) -
                      finite_element_pt(ielem_neigh)
                        ->node_pt(jnod_local_neigh)
                        ->x(i));
                    if (error > node_kill_tol)
                    {
                      oomph_info << "Error in node killing for i " << i << " "
                                 << error << std::endl;
                      stopit = true;
                    }
                  }

                  // Kill node
                  delete finite_element_pt(ielem)->node_pt(jnod_local);
                  killed = true;

                  // Set pointer to neighbour:
                  finite_element_pt(ielem)->node_pt(jnod_local) =
                    finite_element_pt(ielem_neigh)->node_pt(jnod_local_neigh);
                }
                // Not in first layer:
                else
                {
                  // Duplicate node: kill and set pointer to central element
                  if (i0 == 0)
                  {
                    // Neighbour element
                    unsigned ielem_neigh = ielem - 1;

                    // Node in neighbour element
                    unsigned i0_neigh = n_p - 1;
                    unsigned i1_neigh = i1;
                    unsigned i2_neigh = i2;
                    unsigned jnod_local_neigh =
                      i0_neigh + i1_neigh * n_p + i2_neigh * n_p * n_p;


                    // Check:
                    for (unsigned i = 0; i < 3; i++)
                    {
                      double error = std::fabs(
                        finite_element_pt(ielem)->node_pt(jnod_local)->x(i) -
                        finite_element_pt(ielem_neigh)
                          ->node_pt(jnod_local_neigh)
                          ->x(i));
                      if (error > node_kill_tol)
                      {
                        oomph_info << "Error in node killing for i " << i << " "
                                   << error << std::endl;
                        stopit = true;
                      }
                    }

                    // Kill node
                    delete finite_element_pt(ielem)->node_pt(jnod_local);
                    killed = true;

                    // Set pointer to neighbour:
                    finite_element_pt(ielem)->node_pt(jnod_local) =
                      finite_element_pt(ielem_neigh)->node_pt(jnod_local_neigh);
                  }
                }

                // No duplicate node: Copy across to mesh
                if (!killed)
                {
                  Node_pt[node_count] =
                    finite_element_pt(ielem)->node_pt(jnod_local);

                  // Set boundaries:

                  // Back: Boundary 0
                  if ((i2 == 0) && (ilayer == 0))
                  {
                    this->convert_to_boundary_node(Node_pt[node_count]);
                    add_boundary_node(0, Node_pt[node_count]);
                  }

                  // Front: Boundary 4
                  if ((i2 == n_p - 1) && (ilayer == nlayer - 1))
                  {
                    this->convert_to_boundary_node(Node_pt[node_count]);
                    add_boundary_node(4, Node_pt[node_count]);
                  }

                  // Bottom symmetry plane: Boundary 2
                  if (i1 == 0)
                  {
                    this->convert_to_boundary_node(Node_pt[node_count]);
                    add_boundary_node(2, Node_pt[node_count]);
                  }

                  // Tube wall: Boundary 3
                  if (i0 == n_p - 1)
                  {
                    this->convert_to_boundary_node(Node_pt[node_count]);
                    add_boundary_node(3, Node_pt[node_count]);


                    // Get axial boundary coordinate
                    zeta[0] = Xi_lo[0] +
                              (double(ilayer) + double(i2) / double(n_p - 1)) *
                                (Xi_hi[0] - Xi_lo[0]) / double(nlayer);

                    // Get azimuthal boundary coordinate
                    zeta[1] = Xi_lo[1] + double(i1) / double(n_p - 1) * 0.5 *
                                           (Xi_hi[1] - Xi_lo[1]);

                    Node_pt[node_count]->set_coordinates_on_boundary(3, zeta);
                  }

                  // Increment node counter
                  node_count++;
                }
              }
            }
          }

          break;


          // Macro element 2: Top left box
          //--------------------------------
        case 2:

          // Loop over rows in z/s_2-direction
          for (unsigned i2 = 0; i2 < n_p; i2++)
          {
            // Loop over rows in y/s_1-direction
            for (unsigned i1 = 0; i1 < n_p; i1++)
            {
              // Loop over rows in x/s_0-direction
              for (unsigned i0 = 0; i0 < n_p; i0++)
              {
                // Local node number
                unsigned jnod_local = i0 + i1 * n_p + i2 * n_p * n_p;

                // Has the node been killed?
                bool killed = false;

                // First layer of all nodes in s_2 direction gets killed
                // and re-directed to nodes in previous element layer
                if ((i2 == 0) && (ilayer > 0))
                {
                  // Neighbour element
                  unsigned ielem_neigh = ielem - 3;

                  // Node in neighbour element
                  unsigned i0_neigh = i0;
                  unsigned i1_neigh = i1;
                  unsigned i2_neigh = n_p - 1;
                  unsigned jnod_local_neigh =
                    i0_neigh + i1_neigh * n_p + i2_neigh * n_p * n_p;

                  // Check:
                  for (unsigned i = 0; i < 3; i++)
                  {
                    double error = std::fabs(
                      finite_element_pt(ielem)->node_pt(jnod_local)->x(i) -
                      finite_element_pt(ielem_neigh)
                        ->node_pt(jnod_local_neigh)
                        ->x(i));
                    if (error > node_kill_tol)
                    {
                      oomph_info << "Error in node killing for i " << i << " "
                                 << error << std::endl;
                      stopit = true;
                    }
                  }

                  // Kill node
                  delete finite_element_pt(ielem)->node_pt(jnod_local);
                  killed = true;

                  // Set pointer to neighbour:
                  finite_element_pt(ielem)->node_pt(jnod_local) =
                    finite_element_pt(ielem_neigh)->node_pt(jnod_local_neigh);
                }
                // Not in first layer:
                else
                {
                  // Duplicate node: kill and set pointer to node in bottom
                  // right element
                  if (i0 == n_p - 1)
                  {
                    // Neighbour element
                    unsigned ielem_neigh = ielem - 1;

                    // Node in neighbour element
                    unsigned i0_neigh = i1;
                    unsigned i1_neigh = n_p - 1;
                    unsigned i2_neigh = i2;
                    unsigned jnod_local_neigh =
                      i0_neigh + i1_neigh * n_p + i2_neigh * n_p * n_p;

                    // Check:
                    for (unsigned i = 0; i < 3; i++)
                    {
                      double error = std::fabs(
                        finite_element_pt(ielem)->node_pt(jnod_local)->x(i) -
                        finite_element_pt(ielem_neigh)
                          ->node_pt(jnod_local_neigh)
                          ->x(i));
                      if (error > node_kill_tol)
                      {
                        oomph_info << "Error in node killing for i " << i << " "
                                   << error << std::endl;
                        stopit = true;
                      }
                    }

                    // Kill node
                    delete finite_element_pt(ielem)->node_pt(jnod_local);
                    killed = true;

                    // Set pointer to neighbour:
                    finite_element_pt(ielem)->node_pt(jnod_local) =
                      finite_element_pt(ielem_neigh)->node_pt(jnod_local_neigh);
                  }


                  // Duplicate node: kill and set pointer to central element
                  if ((i1 == 0) && (i0 != n_p - 1))
                  {
                    // Neighbour element
                    unsigned ielem_neigh = ielem - 2;

                    // Node in neighbour element
                    unsigned i0_neigh = i0;
                    unsigned i1_neigh = n_p - 1;
                    unsigned i2_neigh = i2;
                    unsigned jnod_local_neigh =
                      i0_neigh + i1_neigh * n_p + i2_neigh * n_p * n_p;

                    // Check:
                    for (unsigned i = 0; i < 3; i++)
                    {
                      double error = std::fabs(
                        finite_element_pt(ielem)->node_pt(jnod_local)->x(i) -
                        finite_element_pt(ielem_neigh)
                          ->node_pt(jnod_local_neigh)
                          ->x(i));
                      if (error > node_kill_tol)
                      {
                        oomph_info << "Error in node killing for i " << i << " "
                                   << error << std::endl;
                        stopit = true;
                      }
                    }

                    // Kill node
                    delete finite_element_pt(ielem)->node_pt(jnod_local);
                    killed = true;

                    // Set pointer to neighbour:
                    finite_element_pt(ielem)->node_pt(jnod_local) =
                      finite_element_pt(ielem_neigh)->node_pt(jnod_local_neigh);
                  }

                  // No duplicate node: Copy across to mesh
                  if (!killed)
                  {
                    Node_pt[node_count] =
                      finite_element_pt(ielem)->node_pt(jnod_local);

                    // Set boundaries:

                    // Back: Boundary 0
                    if ((i2 == 0) && (ilayer == 0))
                    {
                      this->convert_to_boundary_node(Node_pt[node_count]);
                      add_boundary_node(0, Node_pt[node_count]);
                    }

                    // Front: Boundary 4
                    if ((i2 == n_p - 1) && (ilayer == nlayer - 1))
                    {
                      this->convert_to_boundary_node(Node_pt[node_count]);
                      add_boundary_node(4, Node_pt[node_count]);
                    }

                    // Left symmetry plane: Boundary 1
                    if (i0 == 0)
                    {
                      this->convert_to_boundary_node(Node_pt[node_count]);
                      add_boundary_node(1, Node_pt[node_count]);
                    }


                    // Tube wall: Boundary 3
                    if (i1 == n_p - 1)
                    {
                      this->convert_to_boundary_node(Node_pt[node_count]);
                      add_boundary_node(3, Node_pt[node_count]);


                      // Get axial boundary coordinate
                      zeta[0] =
                        Xi_lo[0] +
                        (double(ilayer) + double(i2) / double(n_p - 1)) *
                          (Xi_hi[0] - Xi_lo[0]) / double(nlayer);

                      // Get azimuthal boundary coordinate
                      zeta[1] = Xi_hi[1] - double(i0) / double(n_p - 1) * 0.5 *
                                             (Xi_hi[1] - Xi_lo[1]);

                      Node_pt[node_count]->set_coordinates_on_boundary(3, zeta);
                    }

                    // Increment node counter
                    node_count++;
                  }
                }
              }
            }
          }

          break;
      }
    }

    // Terminate if there's been an error
    if (stopit)
    {
      std::ostringstream error_stream;
      error_stream << "Error in killing nodes\n"
                   << "The most probable cause is that the domain is not\n"
                   << "compatible with the mesh.\n"
                   << "For the QuarterTubeMesh, the domain must be\n"
                   << "topologically consistent with a quarter tube with a\n"
                   << "non-curved centreline.\n";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Setup boundary element lookup schemes
    setup_boundary_element_info();
  }

  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // Algebraic-mesh-version of RefineableQuarterTubeMesh
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Setup algebraic node update data, based on 3 regions, each
  /// undergoing a different node update strategy. These regions are
  /// defined by the three MacroElements in each of the nlayer slices
  /// of the QuarterTubeDomain used to build the mesh.
  /// The Mesh is suspended from the `wall' GeomObject pointed to
  /// by wall_pt. The lower right edge of the mesh is located at the
  /// wall's coordinate xi[1]==xi_lo[1], the upper left edge at
  /// xi[1]=xi_hi[1], i.e. a view looking down the tube length.
  /// The dividing line between the two outer regions is located
  /// at the fraction fract_mid between these two coordinates.
  /// Node updating strategy:
  /// - the starting cross sectional shape along the tube length is
  ///   assumed to be uniform
  /// - the cross sectional shape of the central region remains
  ///   rectangular; the position of its top right corner is located
  ///   at a fixed fraction of the starting width and height of the
  ///   domain measured at xi[1]==xi_lo[1] and xi[1]==xi_hi[1]
  ///   respectively. Nodes in this region are located at fixed
  ///   horizontal and vertical fractions of the region.
  /// - Nodes in the two outer regions (bottom right and top left)
  ///   are located on straight lines running from the edges of the
  ///   central region to the outer wall.
  //======================================================================
  template<class ELEMENT>
  void AlgebraicRefineableQuarterTubeMesh<
    ELEMENT>::setup_algebraic_node_update()
  {
#ifdef PARANOID
    /// Pointer to first algebraic element in central region
    AlgebraicElementBase* algebraic_element_pt =
      dynamic_cast<AlgebraicElementBase*>(Mesh::element_pt(0));

    if (algebraic_element_pt == 0)
    {
      std::ostringstream error_message;
      error_message
        << "Element in AlgebraicRefineableQuarterTubeMesh must be\n ";
      error_message << "derived from AlgebraicElementBase\n";
      error_message << "but it is of type:  "
                    << typeid(Mesh::element_pt(0)).name() << std::endl;
      std::string function_name = "AlgebraicRefineableQuarterTubeMesh::";
      function_name += "setup_algebraic_node_update()";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find number of nodes in an element from the zeroth element
    unsigned nnodes_3d = Mesh::finite_element_pt(0)->nnode();

    // also find number of nodes in 1d line and 2d slice
    unsigned nnodes_1d = Mesh::finite_element_pt(0)->nnode_1d();
    unsigned nnodes_2d = nnodes_1d * nnodes_1d;

    // find node number of a top-left and a bottom-right node in an element
    // (orientation: looking down tube)
    unsigned tl_node = nnodes_2d - nnodes_1d;
    unsigned br_node = nnodes_1d - 1;

    // find x & y distances to top-right node in element 0 - this is the same
    // node as the top-left node of element 1
    double x_c_element = Mesh::finite_element_pt(1)->node_pt(tl_node)->x(0);
    double y_c_element = Mesh::finite_element_pt(1)->node_pt(tl_node)->x(1);

    // Get x-distance to bottom-right edge of wall, i.e. coord of node
    // at bottom-right of bottom-right of element 1
    double x_wall = Mesh::finite_element_pt(1)->node_pt(br_node)->x(0);

    // Get y-distance to top-left edge of wall, i.e. coord of node
    // at top-left of element 2
    double y_wall = Mesh::finite_element_pt(2)->node_pt(tl_node)->x(1);

    // Establish fractional widths in central region
    Lambda_x = Centre_box_size * x_c_element / x_wall;
    Lambda_y = Centre_box_size * y_c_element / y_wall;

    // how many elements are there?
    unsigned nelements = Mesh::nelement();

    // loop through the elements
    for (unsigned e = 0; e < nelements; e++)
    {
      // get pointer to element
      FiniteElement* el_pt = Mesh::finite_element_pt(e);

      // set region id
      unsigned region_id = e % 3;

      // find the first node for which to set up node update info - must
      // bypass the first nnodes_2d nodes after the first 3 elements
      unsigned nstart = nnodes_2d;
      if (e < 3)
      {
        nstart = 0;
      }

      // loop through the nodes,
      for (unsigned n = nstart; n < nnodes_3d; n++)
      {
        // find z coordinate of node
        // NOTE: to implement axial spacing replace z by z_spaced where
        // z_spaced=axial_spacing_fct(z) when finding the GeomObjects
        // and local coords below
        // BUT store z as the third reference value since this is the value
        // required by update_node_update()
        double z = el_pt->node_pt(n)->x(2);

        // Find wall GeomObject and the GeomObject local coordinates
        // at bottom-right edge of wall
        Vector<double> xi(2);
        xi[0] = z;
        xi[1] = QuarterTubeMesh<ELEMENT>::Xi_lo[1];
        GeomObject* obj_br_pt;
        Vector<double> s_br(2);
        this->Wall_pt->locate_zeta(xi, obj_br_pt, s_br);

        // Find wall GeomObject/local coordinate
        // at top-left edge of wall
        xi[1] = QuarterTubeMesh<ELEMENT>::Xi_hi[1];
        GeomObject* obj_tl_pt;
        Vector<double> s_tl(2);
        this->Wall_pt->locate_zeta(xi, obj_tl_pt, s_tl);

        // Element 0: central region
        //--------------------------
        if (region_id == 0)
        {
          // Nodal coordinates in x and y direction
          double x = el_pt->node_pt(n)->x(0);
          double y = el_pt->node_pt(n)->x(1);

          // The update function requires two geometric objects
          Vector<GeomObject*> geom_object_pt(2);

          // Wall GeomObject at bottom right end of wall mesh:
          geom_object_pt[0] = obj_br_pt;

          // Wall GeomObject at top left end of wall mesh:
          geom_object_pt[1] = obj_tl_pt;

          // The update function requires seven parameters:
          Vector<double> ref_value(7);

          // First reference value: fractional x-position inside region
          ref_value[0] = x / x_c_element;

          // Second cfractional y-position inside region
          ref_value[1] = y / y_c_element;

          // Third reference value: initial z coordinate of node
          ref_value[2] = z;

          // Fourth and fifth reference values:
          // local coordinate in GeomObject at bottom-right of wall.
          // Note: must recompute this reference for new nodes
          // since values interpolated from parent nodes will
          // be wrong (this is true for all local coordinates
          // within GeomObjects)
          ref_value[3] = s_br[0];
          ref_value[4] = s_br[1];

          // Sixth and seventh reference values:
          // local coordinate in GeomObject at top-left of wall.
          ref_value[5] = s_tl[0];
          ref_value[6] = s_tl[1];

          // Setup algebraic update for node: Pass update information
          static_cast<AlgebraicNode*>(el_pt->node_pt(n))
            ->add_node_update_info(Central_region, // enumerated ID
                                   this, // mesh
                                   geom_object_pt, // vector of geom objects
                                   ref_value); // vector of ref. values
        }

        // Element 1: Bottom right region
        //-------------------------------
        else if (region_id == 1)
        {
          // Fractional distance between nodes
          double fract = 1.0 / double(nnodes_1d - 1);

          // Fraction in the s_0-direction
          double rho_0 = fract * double(n % nnodes_1d);

          // Fraction in the s_1-direction
          double rho_1 = fract * double((n / nnodes_1d) % nnodes_1d);

          // Find coord of reference point on wall
          xi[1] = QuarterTubeMesh<ELEMENT>::Xi_lo[1] +
                  rho_1 * QuarterTubeMesh<ELEMENT>::Fract_mid *
                    (QuarterTubeMesh<ELEMENT>::Xi_hi[1] -
                     QuarterTubeMesh<ELEMENT>::Xi_lo[1]);

          // Identify wall GeomObject and local coordinate of
          // reference point in GeomObject
          GeomObject* obj_wall_pt;
          Vector<double> s_wall(2);
          this->Wall_pt->locate_zeta(xi, obj_wall_pt, s_wall);

          // The update function requires three geometric objects
          Vector<GeomObject*> geom_object_pt(3);

          // Wall element at bottom-right end of wall mesh:
          geom_object_pt[0] = obj_br_pt;

          // Wall element at top left end of wall mesh:
          geom_object_pt[1] = obj_tl_pt;

          // Wall element at that contians reference point:
          geom_object_pt[2] = obj_wall_pt;

          // The update function requires nine parameters:
          Vector<double> ref_value(9);

          // First reference value: fractional s0-position inside region
          ref_value[0] = rho_0;

          // Second reference value: fractional s1-position inside region
          ref_value[1] = rho_1;

          // Thrid reference value: initial z coord of node
          ref_value[2] = z;

          // Fourth and fifth reference values:
          //   local coordinates at bottom-right of wall.
          ref_value[3] = s_br[0];
          ref_value[4] = s_br[1];

          // Sixth and seventh reference values:
          //   local coordinates at top-left of wall.
          ref_value[5] = s_tl[0];
          ref_value[6] = s_tl[1];

          // Eighth and ninth reference values:
          //   local coordinate of reference point.
          ref_value[7] = s_wall[0];
          ref_value[8] = s_wall[1];

          // Setup algebraic update for node: Pass update information
          static_cast<AlgebraicNode*>(el_pt->node_pt(n))
            ->add_node_update_info(Lower_right_region, // enumerated ID
                                   this, // mesh
                                   geom_object_pt, // vector of geom objects
                                   ref_value); // vector of ref. vals
        }

        // Element 2: Top left region
        //---------------------------
        else if (region_id == 2)
        {
          // Fractional distance between nodes
          double fract = 1.0 / double(nnodes_1d - 1);

          // Fraction in the s_0-direction
          double rho_0 = fract * double(n % nnodes_1d);

          // Fraction in the s_1-direction
          double rho_1 = fract * double((n / nnodes_1d) % nnodes_1d);

          // Find coord of reference point on wall
          xi[1] = QuarterTubeMesh<ELEMENT>::Xi_hi[1] +
                  rho_0 * (1.0 - QuarterTubeMesh<ELEMENT>::Fract_mid) *
                    (QuarterTubeMesh<ELEMENT>::Xi_lo[1] -
                     QuarterTubeMesh<ELEMENT>::Xi_hi[1]);

          // Identify GeomObject and local coordinate at
          // reference point on wall
          GeomObject* obj_wall_pt;
          Vector<double> s_wall(2);
          this->Wall_pt->locate_zeta(xi, obj_wall_pt, s_wall);

          // The update function requires three geometric objects
          Vector<GeomObject*> geom_object_pt(3);

          // Wall element at bottom-right of wall mesh:
          geom_object_pt[0] = obj_br_pt;

          // Wall element at top-left of wall mesh:
          geom_object_pt[1] = obj_tl_pt;

          // Wall element contianing reference point:
          geom_object_pt[2] = obj_wall_pt;

          // The update function requires nine parameters:
          Vector<double> ref_value(9);

          // First reference value: fractional s0-position inside region
          ref_value[0] = rho_0;

          // Second reference value: fractional s1-position inside region
          ref_value[1] = rho_1;

          // Third reference value: initial z coord
          ref_value[2] = z;

          // Fourth and fifth reference values:
          //   local coordinates in bottom-right GeomObject
          ref_value[3] = s_br[0];
          ref_value[4] = s_br[1];

          // Sixth and seventh reference values:
          //   local coordinates in top-left GeomObject
          ref_value[5] = s_tl[0];
          ref_value[6] = s_tl[1];

          // Eighth and ninth reference values:
          //   local coordinates at reference point
          ref_value[7] = s_wall[0];
          ref_value[8] = s_wall[1];

          // Setup algebraic update for node: Pass update information
          static_cast<AlgebraicNode*>(el_pt->node_pt(n))
            ->add_node_update_info(Upper_left_region, // Enumerated ID
                                   this, // mesh
                                   geom_object_pt, // vector of geom objects
                                   ref_value); // vector of ref. vals
        }
      }
    }
  }

  //======================================================================
  /// Algebraic update function: Update in central region according
  /// to wall shape at time level t (t=0: present; t>0: previous)
  //======================================================================
  template<class ELEMENT>
  void AlgebraicRefineableQuarterTubeMesh<ELEMENT>::node_update_central_region(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
#ifdef PARANOID
    // We're updating the nodal positions (!) at time level t
    // and determine them by evaluating the wall GeomObject's
    // position at that time level. I believe this only makes sense
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
        "Trying to update the nodal position at a time level\n";
      error_message += "beyond the number of previous values in the nodes'\n";
      error_message += "position timestepper. This seems highly suspect!\n";
      error_message += "If you're sure the code behaves correctly\n";
      error_message += "in your application, remove this warning \n";
      error_message += "or recompile with PARNOID switched off.\n";

      std::string function_name = "AlgebraicRefineableQuarterTubeMesh::";
      function_name += "node_update_central_region()",
        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Extract references for update in central region by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value(Central_region));

    // Extract geometric objects for update in central region by copy
    // construction
    Vector<GeomObject*> geom_object_pt(
      node_pt->vector_geom_object_pt(Central_region));

    // First reference value: fractional x-position of node inside region
    double rho_x = ref_value[0];

    // Second reference value: fractional y-position of node inside region
    double rho_y = ref_value[1];

    // Wall position in bottom right corner:

    // Pointer to GeomObject
    GeomObject* obj_br_pt = geom_object_pt[0];

    // Local coordinate at bottom-right reference point:
    Vector<double> s_br(2);
    s_br[0] = ref_value[3];
    s_br[1] = ref_value[4];

    // Get wall position
    unsigned n_dim = obj_br_pt->ndim();
    Vector<double> r_br(n_dim);
    obj_br_pt->position(t, s_br, r_br);

    // Wall position in top left corner:

    // Pointer to GeomObject
    GeomObject* obj_tl_pt = geom_object_pt[1];

    // Local coordinate:
    Vector<double> s_tl(2);
    s_tl[0] = ref_value[5];
    s_tl[1] = ref_value[6];

    // Get wall position
    Vector<double> r_tl(n_dim);
    obj_tl_pt->position(t, s_tl, r_tl);

    // Assign new nodal coordinate
    node_pt->x(t, 0) = r_br[0] * Lambda_x * rho_x;
    node_pt->x(t, 1) = r_tl[1] * Lambda_y * rho_y;
    node_pt->x(t, 2) = (r_br[2] + r_tl[2]) * 0.5;
  }


  //====================================================================
  /// Algebraic update function: Update in lower-right region according
  /// to wall shape at time level t (t=0: present; t>0: previous)
  //====================================================================
  template<class ELEMENT>
  void AlgebraicRefineableQuarterTubeMesh<
    ELEMENT>::node_update_lower_right_region(const unsigned& t,
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

      std::string function_name = "AlgebraicRefineableQuarterTubeSectorMesh::";
      function_name += "node_update_lower_right_region()",
        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Extract references for update in lower-right region by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value(Lower_right_region));

    // Extract geometric objects for update in lower-right region
    // by copy construction
    Vector<GeomObject*> geom_object_pt(
      node_pt->vector_geom_object_pt(Lower_right_region));

    // First reference value: fractional s0-position of node inside region
    double rho_0 = ref_value[0];

    // Use s_squashed to modify fractional s0 position
    rho_0 = this->Domain_pt->s_squashed(rho_0);

    // Second reference value: fractional s1-position of node inside region
    double rho_1 = ref_value[1];

    // Wall position in bottom right corner:

    // Pointer to GeomObject
    GeomObject* obj_br_pt = geom_object_pt[0];

    // Local coordinate
    Vector<double> s_br(2);
    s_br[0] = ref_value[3];
    s_br[1] = ref_value[4];

    // Eulerian dimension
    unsigned n_dim = obj_br_pt->ndim();

    // Get wall position
    Vector<double> r_br(n_dim);
    obj_br_pt->position(t, s_br, r_br);

    // Wall position in top left corner:

    // Pointer to GeomObject
    GeomObject* obj_tl_pt = geom_object_pt[1];

    // Local coordinate
    Vector<double> s_tl(2);
    s_tl[0] = ref_value[5];
    s_tl[1] = ref_value[6];

    Vector<double> r_tl(n_dim);
    obj_tl_pt->position(t, s_tl, r_tl);

    // Wall position at reference point:

    // Pointer to GeomObject
    GeomObject* obj_wall_pt = geom_object_pt[2];

    // Local coordinate
    Vector<double> s_wall(2);
    s_wall[0] = ref_value[7];
    s_wall[1] = ref_value[8];

    Vector<double> r_wall(n_dim);
    obj_wall_pt->position(t, s_wall, r_wall);

    // Position Vector to corner of the central region
    Vector<double> r_corner(n_dim);
    r_corner[0] = Lambda_x * r_br[0];
    r_corner[1] = Lambda_y * r_tl[1];
    r_corner[2] = (r_br[2] + r_tl[2]) * 0.5;

    // Position Vector to left edge of region
    Vector<double> r_left(n_dim);
    r_left[0] = r_corner[0];
    r_left[1] = rho_1 * r_corner[1];
    r_left[2] = r_corner[2];

    // Assign new nodal coordinate
    node_pt->x(t, 0) = r_left[0] + rho_0 * (r_wall[0] - r_left[0]);
    node_pt->x(t, 1) = r_left[1] + rho_0 * (r_wall[1] - r_left[1]);
    node_pt->x(t, 2) = r_left[2] + rho_0 * (r_wall[2] - r_left[2]);
  }
  //====================================================================
  /// Algebraic update function: Update in upper left region according
  /// to wall shape at time level t (t=0: present; t>0: previous)
  //====================================================================
  template<class ELEMENT>
  void AlgebraicRefineableQuarterTubeMesh<
    ELEMENT>::node_update_upper_left_region(const unsigned& t,
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

      std::string function_name = "AlgebraicRefineableQuarterTubeMesh::";
      function_name += "node_update_upper_left_region()";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Extract references for update in upper-left region by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value(Upper_left_region));

    // Extract geometric objects for update in upper-left region
    // by copy construction
    Vector<GeomObject*> geom_object_pt(
      node_pt->vector_geom_object_pt(Upper_left_region));

    // First reference value: fractional s0-position of node inside region
    double rho_0 = ref_value[0];

    // Second  reference value: fractional s1-position of node inside region
    double rho_1 = ref_value[1];

    // Use s_squashed to modify fractional s1 position
    rho_1 = this->Domain_pt->s_squashed(rho_1);

    // Wall position in bottom right corner:

    // Pointer to GeomObject
    GeomObject* obj_br_pt = geom_object_pt[0];

    // Eulerian dimension
    unsigned n_dim = obj_br_pt->ndim();

    // Local coordinate
    Vector<double> s_br(2);
    s_br[0] = ref_value[3];
    s_br[1] = ref_value[4];

    // Get wall position
    Vector<double> r_br(n_dim);
    obj_br_pt->position(t, s_br, r_br);

    // Wall position in top left corner:

    // Pointer to GeomObject
    GeomObject* obj_tl_pt = node_pt->geom_object_pt(1);

    // Local coordinate
    Vector<double> s_tl(2);
    s_tl[0] = ref_value[5];
    s_tl[1] = ref_value[6];

    Vector<double> r_tl(n_dim);
    obj_tl_pt->position(t, s_tl, r_tl);

    // Wall position at reference point:

    // Pointer to GeomObject
    GeomObject* obj_wall_pt = node_pt->geom_object_pt(2);

    // Local coordinate
    Vector<double> s_wall(2);
    s_wall[0] = ref_value[7];
    s_wall[1] = ref_value[8];

    Vector<double> r_wall(n_dim);
    obj_wall_pt->position(t, s_wall, r_wall);

    // Position vector to corner of the central region
    Vector<double> r_corner(n_dim);
    r_corner[0] = Lambda_x * r_br[0];
    r_corner[1] = Lambda_y * r_tl[1];
    r_corner[2] = (r_br[2] + r_tl[2]) * 0.5;

    // Position Vector to top face of central region
    Vector<double> r_top(n_dim);
    r_top[0] = rho_0 * r_corner[0];
    r_top[1] = r_corner[1];
    r_top[2] = r_corner[2];

    // Assign new nodal coordinate
    node_pt->x(t, 0) = r_top[0] + rho_1 * (r_wall[0] - r_top[0]);
    node_pt->x(t, 1) = r_top[1] + rho_1 * (r_wall[1] - r_top[1]);
    node_pt->x(t, 2) = r_top[2] + rho_1 * (r_wall[2] - r_top[2]);
  }


  //======================================================================
  /// Update algebraic update function for nodes in region defined by
  /// region_id.
  //======================================================================
  template<class ELEMENT>
  void AlgebraicRefineableQuarterTubeMesh<
    ELEMENT>::update_node_update_in_region(AlgebraicNode*& node_pt,
                                           int& region_id)
  {
    // Extract references by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value(region_id));

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(
      node_pt->vector_geom_object_pt(region_id));

    // Now remove the update info to allow overwriting below
    node_pt->kill_node_update_info(region_id);

    // Fixed reference value: Start coordinate on wall
    double xi_lo = QuarterTubeMesh<ELEMENT>::Xi_lo[1];

    // Fixed reference value: Fractional position of dividing line
    double fract_mid = QuarterTubeMesh<ELEMENT>::Fract_mid;

    // Fixed reference value: End coordinate on wall
    double xi_hi = QuarterTubeMesh<ELEMENT>::Xi_hi[1];

    // get initial z-coordinate of node
    // NOTE: use modified z found using axial_spacing_fct()
    // to implement axial spacing
    double z = ref_value[2];

    // Update reference to wall point in bottom right corner
    //------------------------------------------------------

    // Wall position in bottom right corner:
    Vector<double> xi(2);
    xi[0] = z;
    xi[1] = xi_lo;

    // Find corresponding wall element/local coordinate
    GeomObject* obj_br_pt;
    Vector<double> s_br(2);
    this->Wall_pt->locate_zeta(xi, obj_br_pt, s_br);

    // Wall element at bottom right end of wall mesh:
    geom_object_pt[0] = obj_br_pt;

    // Local coordinate in this wall element.
    ref_value[3] = s_br[0];
    ref_value[4] = s_br[1];


    // Update reference to wall point in upper left corner
    //----------------------------------------------------

    // Wall position in top left corner
    xi[1] = xi_hi;

    // Find corresponding wall element/local coordinate
    GeomObject* obj_tl_pt;
    Vector<double> s_tl(2);
    this->Wall_pt->locate_zeta(xi, obj_tl_pt, s_tl);

    // Wall element at top left end of wall mesh:
    geom_object_pt[1] = obj_tl_pt;

    // Local coordinate in this wall element.
    ref_value[5] = s_tl[0];
    ref_value[6] = s_tl[1];

    if (region_id != Central_region)
    {
      // Update reference to reference point on wall
      //--------------------------------------------

      // Reference point on wall
      if (region_id == Lower_right_region)
      {
        // Second reference value: fractional s1-position of node inside region
        double rho_1 = ref_value[1];

        // position of reference point
        xi[1] = xi_lo + rho_1 * fract_mid * (xi_hi - xi_lo);
      }
      else // case of Upper_left region
      {
        // First reference value: fractional s0-position of node inside region
        double rho_0 = ref_value[0];

        // position of reference point
        xi[1] = xi_hi - rho_0 * (1.0 - fract_mid) * (xi_hi - xi_lo);
      }

      // Identify wall element number and local coordinate of
      // reference point on wall
      GeomObject* obj_wall_pt;
      Vector<double> s_wall(2);
      this->Wall_pt->locate_zeta(xi, obj_wall_pt, s_wall);

      // Wall element at that contians reference point:
      geom_object_pt[2] = obj_wall_pt;

      // Local coordinate in this wall element.
      ref_value[7] = s_wall[0];
      ref_value[8] = s_wall[1];
    }

    // Setup algebraic update for node: Pass update information
    node_pt->add_node_update_info(region_id, // Enumerated ID
                                  this, // mesh
                                  geom_object_pt, // vector of geom objects
                                  ref_value); // vector of ref. vals
  }


} // namespace oomph
#endif
