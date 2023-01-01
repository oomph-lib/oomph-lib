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
#ifndef OOMPH_EIGHTH_SPHERE_MESH_TEMPLATE_CC
#define OOMPH_EIGHTH_SPHERE_MESH_TEMPLATE_CC

#include "eighth_sphere_mesh.template.h"


namespace oomph
{
  //======================================================================
  /// Constructor for the eighth of a sphere mesh: Pass timestepper;
  /// defaults to static default timestepper.
  //======================================================================
  template<class ELEMENT>
  EighthSphereMesh<ELEMENT>::EighthSphereMesh(const double& radius,
                                              TimeStepper* time_stepper_pt)
    : Radius(radius)
  {
    // Mesh can only be built with 3D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

    // Set the number of boundaries
    this->set_nboundary(4);

    // Provide storage for the four elements
    this->Element_pt.resize(4);

    // Set the domain pointer: Pass the radius of the sphere
    Domain_pt = new EighthSphereDomain(Radius);


    Vector<double> s(3), s_fraction(3);
    Vector<double> r(3);

    // Create first element
    //--------------------
    this->Element_pt[0] = new ELEMENT;

    // Give it nodes

    // How many 1D nodes are there?
    unsigned nnode1d = this->finite_element_pt(0)->nnode_1d();

    // Loop over nodes
    for (unsigned i = 0; i < nnode1d; i++)
    {
      for (unsigned j = 0; j < nnode1d; j++)
      {
        for (unsigned k = 0; k < nnode1d; k++)
        {
          unsigned jnod = k * nnode1d * nnode1d + j * nnode1d + i;

          // If we're on a boundary, make a boundary node
          if ((i == 0) || (j == 0) || (k == 0))
          {
            // Allocate the node according to the element's wishes
            this->Node_pt.push_back(
              this->finite_element_pt(0)->construct_boundary_node(
                jnod, time_stepper_pt));
          }
          // Otherwise it's a normal node
          else
          {
            // Allocate the node according to the element's wishes
            this->Node_pt.push_back(this->finite_element_pt(0)->construct_node(
              jnod, time_stepper_pt));
          }

          // Work out the node's coordinates in the finite element's local
          // coordinate system
          this->finite_element_pt(0)->local_fraction_of_node(jnod, s_fraction);


          // Get the position of the node from macro element mapping
          s[0] = -1.0 + 2.0 * s_fraction[0];
          s[1] = -1.0 + 2.0 * s_fraction[1];
          s[2] = -1.0 + 2.0 * s_fraction[2];
          Domain_pt->macro_element_pt(0)->macro_map(s, r);

          // Assign coordinates
          this->finite_element_pt(0)->node_pt(jnod)->x(0) = r[0];
          this->finite_element_pt(0)->node_pt(jnod)->x(1) = r[1];
          this->finite_element_pt(0)->node_pt(jnod)->x(2) = r[2];

          // Add the node to the appropriate boundary
          if (i == 0)
            add_boundary_node(0, this->finite_element_pt(0)->node_pt(jnod));
          if (j == 0)
            add_boundary_node(1, this->finite_element_pt(0)->node_pt(jnod));
          if (k == 0)
            add_boundary_node(2, this->finite_element_pt(0)->node_pt(jnod));
        }
      }
    }


    // Create a second element
    //------------------------
    this->Element_pt[1] = new ELEMENT;
    ;

    // Give it nodes

    // Loop over nodes
    for (unsigned i = 0; i < nnode1d; i++)
    {
      for (unsigned j = 0; j < nnode1d; j++)
      {
        for (unsigned k = 0; k < nnode1d; k++)
        {
          unsigned jnod = k * nnode1d * nnode1d + j * nnode1d + i;

          // If node has not been yet created, create it
          if (i > 0)
          {
            // If the node is on a boundary, make a boundary node
            if ((i == nnode1d - 1) || (j == 0) || (k == 0))
            {
              this->Node_pt.push_back(
                this->finite_element_pt(1)->construct_boundary_node(
                  jnod, time_stepper_pt));
            }
            // Otherwise make a normal node
            else
            {
              this->Node_pt.push_back(
                this->finite_element_pt(1)->construct_node(jnod,
                                                           time_stepper_pt));
            }

            // Work out the node's coordinates in the finite element's local
            // coordinate system
            this->finite_element_pt(1)->local_fraction_of_node(jnod,
                                                               s_fraction);

            // Get the position of the node from macro element mapping
            s[0] = -1.0 + 2.0 * s_fraction[0];
            s[1] = -1.0 + 2.0 * s_fraction[1];
            s[2] = -1.0 + 2.0 * s_fraction[2];
            Domain_pt->macro_element_pt(1)->macro_map(s, r);

            // Assign coordinate
            this->finite_element_pt(1)->node_pt(jnod)->x(0) = r[0];
            this->finite_element_pt(1)->node_pt(jnod)->x(1) = r[1];
            this->finite_element_pt(1)->node_pt(jnod)->x(2) = r[2];

            // Add the node to the approriate boundary
            if (j == 0)
              add_boundary_node(1, this->finite_element_pt(1)->node_pt(jnod));
            if (k == 0)
              add_boundary_node(2, this->finite_element_pt(1)->node_pt(jnod));
            if (i == nnode1d - 1)
              add_boundary_node(3, this->finite_element_pt(1)->node_pt(jnod));
          }

          // ...else use the node already created
          else
          {
            this->finite_element_pt(1)->node_pt(jnod) =
              this->finite_element_pt(0)->node_pt(jnod + nnode1d - 1);
          }
        }
      }
    }


    // Create a third element
    //------------------------
    this->Element_pt[2] = new ELEMENT;

    // Give it nodes

    // Loop over nodes
    for (unsigned i = 0; i < nnode1d; i++)
    {
      for (unsigned j = 0; j < nnode1d; j++)
      {
        for (unsigned k = 0; k < nnode1d; k++)
        {
          unsigned jnod = k * nnode1d * nnode1d + j * nnode1d + i;

          // If the node has not been yet created, create it
          if ((i < nnode1d - 1) && (j > 0))
          {
            // If it's on a boundary, make a boundary node
            if ((i == 0) || (j == nnode1d - 1) || (k == 0))
            {
              // Allocate the node according to the element's wishes
              this->Node_pt.push_back(
                this->finite_element_pt(2)->construct_boundary_node(
                  jnod, time_stepper_pt));
            }
            // Otherwise allocate a normal node
            else
            {
              // Allocate the node according to the element's wishes
              this->Node_pt.push_back(
                this->finite_element_pt(2)->construct_node(jnod,
                                                           time_stepper_pt));
            }

            // Work out the node's coordinates in the finite element's local
            // coordinate system
            this->finite_element_pt(2)->local_fraction_of_node(jnod,
                                                               s_fraction);

            // Get the position of the node from macro element mapping
            s[0] = -1.0 + 2.0 * s_fraction[0];
            s[1] = -1.0 + 2.0 * s_fraction[1];
            s[2] = -1.0 + 2.0 * s_fraction[2];
            Domain_pt->macro_element_pt(2)->macro_map(s, r);

            // Assign coordinates
            this->finite_element_pt(2)->node_pt(jnod)->x(0) = r[0];
            this->finite_element_pt(2)->node_pt(jnod)->x(1) = r[1];
            this->finite_element_pt(2)->node_pt(jnod)->x(2) = r[2];

            // Add the node to the appropriate boundary
            if (i == 0)
              add_boundary_node(0, this->finite_element_pt(2)->node_pt(jnod));
            if (k == 0)
              add_boundary_node(2, this->finite_element_pt(2)->node_pt(jnod));
            if (j == nnode1d - 1)
              add_boundary_node(3, this->finite_element_pt(2)->node_pt(jnod));
          }

          // ...else use the nodes already created
          else
          {
            // If the node belongs to the element 0
            if (j == 0)
              this->finite_element_pt(2)->node_pt(jnod) =
                this->finite_element_pt(0)->node_pt(jnod +
                                                    nnode1d * (nnode1d - 1));

            // ...else it belongs to the element 1
            else if (i == nnode1d - 1)
              this->finite_element_pt(2)->node_pt(jnod) =
                this->finite_element_pt(1)->node_pt(nnode1d * nnode1d * k + j +
                                                    i * nnode1d);
          }
        }
      }
    }


    // Create the fourth element
    //-------------------------
    this->Element_pt[3] = new ELEMENT;

    // Give it nodes

    // Loop over nodes
    for (unsigned i = 0; i < nnode1d; i++)
    {
      for (unsigned j = 0; j < nnode1d; j++)
      {
        for (unsigned k = 0; k < nnode1d; k++)
        {
          unsigned jnod = k * nnode1d * nnode1d + j * nnode1d + i;

          // If the node has not been yet created, create it
          if ((k > 0) && (i < nnode1d - 1) && (j < nnode1d - 1))
          {
            // If it's on a boundary, allocate a boundary node
            if ((i == 0) || (j == 0) || (k == nnode1d - 1))
            {
              // Allocate the node according to the element's wishes
              this->Node_pt.push_back(
                this->finite_element_pt(3)->construct_boundary_node(
                  jnod, time_stepper_pt));
            }
            // Otherwise allocate a normal node
            else
            {
              // Allocate the node according to the element's wishes
              this->Node_pt.push_back(
                this->finite_element_pt(3)->construct_node(jnod,
                                                           time_stepper_pt));
            }

            // Work out the node's coordinates in the finite element's local
            // coordinate system
            this->finite_element_pt(3)->local_fraction_of_node(jnod,
                                                               s_fraction);

            // Get the position of the node from macro element mapping
            s[0] = -1.0 + 2.0 * s_fraction[0];
            s[1] = -1.0 + 2.0 * s_fraction[1];
            s[2] = -1.0 + 2.0 * s_fraction[2];
            Domain_pt->macro_element_pt(3)->macro_map(s, r);

            // Assign coordinates
            this->finite_element_pt(3)->node_pt(jnod)->x(0) = r[0];
            this->finite_element_pt(3)->node_pt(jnod)->x(1) = r[1];
            this->finite_element_pt(3)->node_pt(jnod)->x(2) = r[2];

            // Add the node to the appropriate boundary
            if (i == 0)
              add_boundary_node(0, this->finite_element_pt(3)->node_pt(jnod));
            if (j == 0)
              add_boundary_node(1, this->finite_element_pt(3)->node_pt(jnod));
            if (k == nnode1d - 1)
              add_boundary_node(3, this->finite_element_pt(3)->node_pt(jnod));
          }

          // ...otherwise the node was already created: use it.
          else
          {
            // if k=0 then the node belongs to element 0
            if (k == 0)
            {
              this->finite_element_pt(3)->node_pt(jnod) =
                this->finite_element_pt(0)->node_pt(jnod + nnode1d * nnode1d *
                                                             (nnode1d - 1));
            }
            else
            {
              // else if i==nnode1d-1 the node already exists in element 1
              if (i == nnode1d - 1)
              {
                this->finite_element_pt(3)->node_pt(jnod) =
                  this->finite_element_pt(1)->node_pt(
                    k + i * nnode1d * nnode1d + j * nnode1d);
              }
              else
              // else, the node exists in element 2
              {
                this->finite_element_pt(3)->node_pt(jnod) =
                  this->finite_element_pt(2)->node_pt(i + k * nnode1d +
                                                      j * nnode1d * nnode1d);
              }
            }
          }
        }
      }
    }

    // Setup boundary element lookup schemes
    setup_boundary_element_info();
  }

} // namespace oomph
#endif
