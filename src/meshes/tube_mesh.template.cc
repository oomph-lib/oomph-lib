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
#ifndef OOMPH_TUBE_MESH_TEMPLATE_CC
#define OOMPH_TUBE_MESH_TEMPLATE_CC

#include "tube_mesh.template.h"


namespace oomph
{
  //====================================================================
  /// \short Constructor for deformable quarter tube mesh class.
  /// The domain is specified by the GeomObject that
  /// identifies the entire volume.
  //====================================================================
  template<class ELEMENT>
  TubeMesh<ELEMENT>::TubeMesh(GeomObject* volume_pt,
                              const Vector<double>& centreline_limits,
                              const Vector<double>& theta_positions,
                              const Vector<double>& radius_box,
                              const unsigned& nlayer,
                              TimeStepper* time_stepper_pt)
    : Volume_pt(volume_pt)
  {
    // Mesh can only be built with 3D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

    // Build macro element-based domain
    Domain_pt = new TubeDomain(
      volume_pt, centreline_limits, theta_positions, radius_box, nlayer);

    // Set the number of boundaries
    set_nboundary(3);

    // We have only bothered to parametrise boundary 1
    Boundary_coordinate_exists[1] = true;

    // Allocate the store for the elements
    unsigned nelem = 5 * nlayer;
    Element_pt.resize(nelem);

    // Create  dummy element so we can determine the number of nodes
    ELEMENT* dummy_el_pt = new ELEMENT;

    // Read out the number of linear points in the element
    unsigned n_p = dummy_el_pt->nnode_1d();

    // Kill the element
    delete dummy_el_pt;

    // Can now allocate the store for the nodes
    unsigned nnodes_total =
      (n_p * n_p + 4 * (n_p - 1) * (n_p - 1)) * (1 + nlayer * (n_p - 1));
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

    // Define pine
    const double pi = MathematicalConstants::Pi;

    // Loop over elements
    for (unsigned ielem = 0; ielem < nelem; ielem++)
    {
      // Which layer?
      unsigned ilayer = unsigned(ielem / 5);

      // Which macro element?
      switch (ielem % 5)
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
                  unsigned ielem_neigh = ielem - 5;

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

                  // Front: Boundary 2
                  if ((i2 == n_p - 1) && (ilayer == nlayer - 1))
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

          // Macro element 1: Bottom
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
                  unsigned ielem_neigh = ielem - 5;

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
                  if (i1 == (n_p - 1))
                  {
                    // Neighbour element
                    unsigned ielem_neigh = ielem - 1;

                    // Node in neighbour element
                    unsigned i0_neigh = i0;
                    unsigned i1_neigh = 0;
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

                  // Front: Boundary 2
                  if ((i2 == n_p - 1) && (ilayer == nlayer - 1))
                  {
                    this->convert_to_boundary_node(Node_pt[node_count]);
                    add_boundary_node(2, Node_pt[node_count]);
                  }

                  // Bottom: outer boundary
                  if (i1 == 0)
                  {
                    this->convert_to_boundary_node(Node_pt[node_count]);
                    add_boundary_node(1, Node_pt[node_count]);

                    // Get axial boundary coordinate
                    zeta[0] = centreline_limits[0] +
                              (double(ilayer) + double(i2) / double(n_p - 1)) *
                                (centreline_limits[1] - centreline_limits[0]) /
                                double(nlayer);

                    // Get azimuthal boundary coordinate
                    zeta[1] = theta_positions[0] +
                              double(i1) / double(n_p - 1) * 0.5 *
                                (theta_positions[1] - theta_positions[0]);

                    Node_pt[node_count]->set_coordinates_on_boundary(1, zeta);
                  }
                  // Increment node counter
                  node_count++;
                }
              }
            }
          } // End of loop over nodes

          break;


          // Macro element 2: Right element
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
                  unsigned ielem_neigh = ielem - 5;

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
                  // Duplicate node: kill and set pointer to previous element
                  if (i1 == 0)
                  {
                    // Neighbour element
                    unsigned ielem_neigh = ielem - 1;

                    // Node in neighbour element
                    unsigned i0_neigh = n_p - 1;
                    unsigned i1_neigh = n_p - 1 - i0;
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
                  if ((i0 == 0) && (i1 != 0))
                  {
                    // Neighbour element
                    unsigned ielem_neigh = ielem - 2;

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

                    // Front: Boundary 2
                    if ((i2 == n_p - 1) && (ilayer == nlayer - 1))
                    {
                      this->convert_to_boundary_node(Node_pt[node_count]);
                      add_boundary_node(2, Node_pt[node_count]);
                    }

                    // Tube boundary
                    if (i0 == n_p - 1)
                    {
                      this->convert_to_boundary_node(Node_pt[node_count]);
                      add_boundary_node(1, Node_pt[node_count]);

                      // Get axial boundary coordinate
                      zeta[0] =
                        centreline_limits[0] +
                        (double(ilayer) + double(i2) / double(n_p - 1)) *
                          (centreline_limits[1] - centreline_limits[0]) /
                          double(nlayer);

                      // Get azimuthal boundary coordinate
                      zeta[1] = theta_positions[1] +
                                double(i1) / double(n_p - 1) * 0.5 *
                                  (theta_positions[2] - theta_positions[1]);

                      Node_pt[node_count]->set_coordinates_on_boundary(1, zeta);
                    }

                    // Increment node counter
                    node_count++;
                  }
                }
              }
            }
          }

          break;

          // Macro element 3: Top element
          //--------------------------------
        case 3:

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
                  unsigned ielem_neigh = ielem - 5;

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
                  // Duplicate node: kill and set pointer to previous element
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
                    unsigned ielem_neigh = ielem - 3;

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

                    // Front: Boundary 2
                    if ((i2 == n_p - 1) && (ilayer == nlayer - 1))
                    {
                      this->convert_to_boundary_node(Node_pt[node_count]);
                      add_boundary_node(2, Node_pt[node_count]);
                    }

                    // Tube boundary
                    if (i1 == n_p - 1)
                    {
                      this->convert_to_boundary_node(Node_pt[node_count]);
                      add_boundary_node(1, Node_pt[node_count]);

                      // Get axial boundary coordinate
                      zeta[0] =
                        centreline_limits[0] +
                        (double(ilayer) + double(i2) / double(n_p - 1)) *
                          (centreline_limits[1] - centreline_limits[0]) /
                          double(nlayer);

                      // Get azimuthal boundary coordinate
                      zeta[1] = theta_positions[3] +
                                double(i0) / double(n_p - 1) * 0.5 *
                                  (theta_positions[2] - theta_positions[3]);

                      Node_pt[node_count]->set_coordinates_on_boundary(1, zeta);
                    }

                    // Increment node counter
                    node_count++;
                  }
                }
              }
            }
          }
          break;


          // Macro element 4: Left element
          //--------------------------------
        case 4:

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
                  unsigned ielem_neigh = ielem - 5;

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
                  // Duplicate node: kill and set pointer to previous element
                  if (i1 == n_p - 1)
                  {
                    // Neighbour element
                    unsigned ielem_neigh = ielem - 1;

                    // Node in neighbour element
                    unsigned i0_neigh = 0;
                    unsigned i1_neigh = n_p - 1 - i0;
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
                  if ((i0 == n_p - 1) && (i1 != n_p - 1))
                  {
                    // Neighbour element
                    unsigned ielem_neigh = ielem - 4;

                    // Node in neighbour element
                    unsigned i0_neigh = 0;
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

                  // Duplicate node: kill and set pointer to other ring element
                  if ((i1 == 0) && (i0 != n_p - 1))
                  {
                    // Neighbour element
                    unsigned ielem_neigh = ielem - 3;

                    // Node in neighbour element
                    unsigned i0_neigh = 0;
                    unsigned i1_neigh = i0;
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

                    // Front: Boundary 2
                    if ((i2 == n_p - 1) && (ilayer == nlayer - 1))
                    {
                      this->convert_to_boundary_node(Node_pt[node_count]);
                      add_boundary_node(2, Node_pt[node_count]);
                    }

                    // Tube boundary
                    if (i0 == 0)
                    {
                      this->convert_to_boundary_node(Node_pt[node_count]);
                      add_boundary_node(1, Node_pt[node_count]);

                      // Get axial boundary coordinate
                      zeta[0] =
                        centreline_limits[0] +
                        (double(ilayer) + double(i2) / double(n_p - 1)) *
                          (centreline_limits[1] - centreline_limits[0]) /
                          double(nlayer);

                      // Get azimuthal boundary coordinate
                      zeta[1] =
                        theta_positions[0] + 2.0 * pi +
                        double(i1) / double(n_p - 1) * 0.5 *
                          (theta_positions[3] - theta_positions[0] + 2.0 * pi);

                      Node_pt[node_count]->set_coordinates_on_boundary(1, zeta);
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
                   << "For the TubeMesh, the domain must be\n"
                   << "topologically consistent with a quarter tube with a\n"
                   << "non-curved centreline.\n";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Setup boundary element lookup schemes
    setup_boundary_element_info();
  }

} // namespace oomph
#endif
