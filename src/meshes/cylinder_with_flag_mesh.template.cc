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
#ifndef OOMPH_CYLINDER_WITH_FLAG_MESH_TEMPLATE_CC
#define OOMPH_CYLINDER_WITH_FLAG_MESH_TEMPLATE_CC

// Include the headers file
#include "cylinder_with_flag_mesh.template.h"


namespace oomph
{
  //=============================================================
  /// Constructor. Pass the pointers to the GeomObjects that parametrise
  /// the cylinder, the three edges of the flag, the length and height of the
  /// domain, the length and height of the flag, the coordinates of the
  /// centre of the cylinder and its radius. Timestepper defaults to Steady
  /// default timestepper.
  //=============================================================
  template<class ELEMENT>
  CylinderWithFlagMesh<ELEMENT>::CylinderWithFlagMesh(
    Circle* cylinder_pt,
    GeomObject* top_flag_pt,
    GeomObject* bottom_flag_pt,
    GeomObject* tip_flag_pt,
    const double& length,
    const double& height,
    const double& flag_length,
    const double& flag_height,
    const double& centre_x,
    const double& centre_y,
    const double& a,
    TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Create the domain
    Domain_pt = new CylinderWithFlagDomain(cylinder_pt,
                                           top_flag_pt,
                                           bottom_flag_pt,
                                           tip_flag_pt,
                                           length,
                                           height,
                                           flag_length,
                                           flag_height,
                                           centre_x,
                                           centre_y,
                                           a);

    // Initialise the node counter
    unsigned long node_count = 0;

    // Vectors used to get data from domains
    Vector<double> s(2);
    Vector<double> r(2);

    // Setup temporary storage for the Node
    Vector<Node*> tmp_node_pt;

    // Now blindly loop over the macro elements and associate and finite
    // element with each
    unsigned nmacro_element = Domain_pt->nmacro_element();
    for (unsigned e = 0; e < nmacro_element; e++)
    {
      // Create the FiniteElement and add to the Element_pt Vector
      Element_pt.push_back(new ELEMENT);

      // Read out the number of linear points in the element
      unsigned np =
        dynamic_cast<ELEMENT*>(this->finite_element_pt(e))->nnode_1d();

      // Loop over nodes in the column
      for (unsigned l1 = 0; l1 < np; l1++)
      {
        // Loop over the nodes in the row
        for (unsigned l2 = 0; l2 < np; l2++)
        {
          // Allocate the memory for the node
          tmp_node_pt.push_back(
            this->finite_element_pt(e)->construct_boundary_node(
              l1 * np + l2, time_stepper_pt));

          // Read out the position of the node from the macro element
          s[0] = -1.0 + 2.0 * (double)l2 / (double)(np - 1);
          s[1] = -1.0 + 2.0 * (double)l1 / (double)(np - 1);
          Domain_pt->macro_element_pt(e)->macro_map(s, r);

          // Set the position of the node
          tmp_node_pt[node_count]->x(0) = r[0];
          tmp_node_pt[node_count]->x(1) = r[1];

          // Increment the node number
          node_count++;
        }
      }
    } // End of loop over macro elements

    // Now the elements have been created, but there will be nodes in
    // common, need to loop over the common edges and sort it, by reassigning
    // pointers and the deleting excess nodes

    // Read out the number of linear points in the element
    unsigned np =
      dynamic_cast<ELEMENT*>(this->finite_element_pt(0))->nnode_1d();

    // Edge between Elements 0 and 1
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 1 to be the same as in element 0
      this->finite_element_pt(1)->node_pt(n * np) =
        this->finite_element_pt(0)->node_pt((np - 1) * np + np - 1 - n);

      // Remove the nodes in element 1 from the temporaray node list
      delete tmp_node_pt[np * np + n * np];
      tmp_node_pt[np * np + n * np] = 0;
    }

    // Edge between Elements 1 and 2
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 2 to be the same as in element 1
      this->finite_element_pt(2)->node_pt(n * np) =
        this->finite_element_pt(1)->node_pt(np * n + np - 1);

      // Remove the nodes in element 2 from the temporaray node list
      delete tmp_node_pt[2 * np * np + n * np];
      tmp_node_pt[2 * np * np + n * np] = 0;
    }

    // Edge between Element 2 and 3
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 3 to be the same as in element 2
      this->finite_element_pt(3)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(2)->node_pt(np * n + np - 1);

      // Remove the nodes in element 3 from the temporaray node list
      delete tmp_node_pt[3 * np * np + np * (np - 1) + n];
      tmp_node_pt[3 * np * np + np * (np - 1) + n] = 0;
    }

    // Edge between Element 5 and 4
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 4 to be the same as in element 5
      this->finite_element_pt(4)->node_pt(n) =
        this->finite_element_pt(5)->node_pt(np * (np - n - 1) + np - 1);

      // Remove the nodes in element 4 from the temporaray node list
      delete tmp_node_pt[4 * np * np + n];
      tmp_node_pt[4 * np * np + n] = 0;
    }

    // Edge between Elements 6 and 5
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 5 to be the same as in element 6
      this->finite_element_pt(5)->node_pt(n * np) =
        this->finite_element_pt(6)->node_pt(np * n + np - 1);

      // Remove the nodes in element 5 from the temporaray node list
      delete tmp_node_pt[5 * np * np + n * np];
      tmp_node_pt[5 * np * np + n * np] = 0;
    }

    // Edge between Elements 0 and 6
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 6 to be the same as in element 0
      this->finite_element_pt(6)->node_pt(n * np) =
        this->finite_element_pt(0)->node_pt(n);

      // Remove the nodes in element 6 from the temporaray node list
      delete tmp_node_pt[6 * np * np + n * np];
      tmp_node_pt[6 * np * np + n * np] = 0;
    }

    // Edge between Elements 2 and 7
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 7 to be the same as in element 2
      this->finite_element_pt(7)->node_pt(n * np) =
        this->finite_element_pt(2)->node_pt((np - 1) * np + np - 1 - n);

      // Remove the nodes in element 7 from the temporaray node list
      delete tmp_node_pt[7 * np * np + n * np];
      tmp_node_pt[7 * np * np + n * np] = 0;
    }

    // Edge between Elements 3 and 8
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 8 to be the same as in element 3
      this->finite_element_pt(8)->node_pt(n * np) =
        this->finite_element_pt(3)->node_pt(np * n + np - 1);

      // Remove the nodes in element 8 from the temporaray node list
      delete tmp_node_pt[8 * np * np + n * np];
      tmp_node_pt[8 * np * np + n * np] = 0;
    }

    // Edge between Elements 4 and 9
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 9 to be the same as in element 4
      this->finite_element_pt(9)->node_pt(n * np) =
        this->finite_element_pt(4)->node_pt(np * n + np - 1);

      // Remove the nodes in element 9 from the temporaray node list
      delete tmp_node_pt[9 * np * np + n * np];
      tmp_node_pt[9 * np * np + n * np] = 0;
    }


    // Edge between Elements 5 and 10
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 10 to be the same as in element 5
      this->finite_element_pt(10)->node_pt(n * np) =
        this->finite_element_pt(5)->node_pt(n);

      // Remove the nodes in element 10 from the temporaray node list
      delete tmp_node_pt[10 * np * np + n * np];
      tmp_node_pt[10 * np * np + n * np] = 0;
    }


    // Edge between Elements 7 and 11
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 11 to be the same as in element 7
      this->finite_element_pt(11)->node_pt(n * np) =
        this->finite_element_pt(7)->node_pt(np * n + np - 1);

      // Remove the nodes in element 11 from the temporaray node list
      delete tmp_node_pt[11 * np * np + n * np];
      tmp_node_pt[11 * np * np + n * np] = 0;
    }

    // Edge between Elements 8 and 12
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 12 to be the same as in element 8
      this->finite_element_pt(12)->node_pt(n * np) =
        this->finite_element_pt(8)->node_pt(np * n + np - 1);

      // Remove the nodes in element 12 from the temporaray node list
      delete tmp_node_pt[12 * np * np + n * np];
      tmp_node_pt[12 * np * np + n * np] = 0;
    }

    // Edge between Elements 9 and 13
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 13 to be the same as in element 9
      this->finite_element_pt(13)->node_pt(n * np) =
        this->finite_element_pt(9)->node_pt(np * n + np - 1);

      // Remove the nodes in element 13 from the temporaray node list
      delete tmp_node_pt[13 * np * np + n * np];
      tmp_node_pt[13 * np * np + n * np] = 0;
    }

    // Edge between Elements 10 and 14
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 14 to be the same as in element 10
      this->finite_element_pt(14)->node_pt(n * np) =
        this->finite_element_pt(10)->node_pt(np * n + np - 1);

      // Remove the nodes in element 14 from the temporaray node list
      delete tmp_node_pt[14 * np * np + n * np];
      tmp_node_pt[14 * np * np + n * np] = 0;
    }


    // Edge between Elements 7 and 8
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 8 to be the same as in element 7
      this->finite_element_pt(8)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(7)->node_pt(n);

      // Remove the nodes in element 8 from the temporaray node list
      delete tmp_node_pt[8 * np * np + np * (np - 1) + n];
      tmp_node_pt[8 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Elements 9 and 10
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 10 to be the same as in element 9
      this->finite_element_pt(10)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(9)->node_pt(n);

      // Remove the nodes in element 10 from the temporaray node list
      delete tmp_node_pt[10 * np * np + np * (np - 1) + n];
      tmp_node_pt[10 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Elements 11 and 15
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 15 to be the same as in element 11
      this->finite_element_pt(15)->node_pt(n * np) =
        this->finite_element_pt(11)->node_pt(np * n + np - 1);

      // Remove the nodes in element 15 from the temporaray node list
      delete tmp_node_pt[15 * np * np + n * np];
      tmp_node_pt[15 * np * np + n * np] = 0;
    }

    // Edge between Elements 12 and 16
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 16 to be the same as in element 12
      this->finite_element_pt(16)->node_pt(n * np) =
        this->finite_element_pt(12)->node_pt(np * n + np - 1);

      // Remove the nodes in element 16 from the temporaray node list
      delete tmp_node_pt[16 * np * np + n * np];
      tmp_node_pt[16 * np * np + n * np] = 0;
    }

    // Edge between Elements 13 and 17
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 17 to be the same as in element 13
      this->finite_element_pt(17)->node_pt(n * np) =
        this->finite_element_pt(13)->node_pt(np * n + np - 1);

      // Remove the nodes in element 17 from the temporaray node list
      delete tmp_node_pt[17 * np * np + n * np];
      tmp_node_pt[17 * np * np + n * np] = 0;
    }

    // Edge between Elements 14 and 18
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 18 to be the same as in element 14
      this->finite_element_pt(18)->node_pt(n * np) =
        this->finite_element_pt(14)->node_pt(np * n + np - 1);

      // Remove the nodes in element 18 from the temporaray node list
      delete tmp_node_pt[18 * np * np + n * np];
      tmp_node_pt[18 * np * np + n * np] = 0;
    }


    // Edge between Elements 11 and 12
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 12 to be the same as in element 11
      this->finite_element_pt(12)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(11)->node_pt(n);

      // Remove the nodes in element 12 from the temporaray node list
      delete tmp_node_pt[12 * np * np + np * (np - 1) + n];
      tmp_node_pt[12 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Elements 13 and 14
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 14 to be the same as in element 13
      this->finite_element_pt(14)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(13)->node_pt(n);

      // Remove the nodes in element 14 from the temporaray node list
      delete tmp_node_pt[14 * np * np + np * (np - 1) + n];
      tmp_node_pt[14 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Element 15 and 19
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 19 to be the same as in element 15
      this->finite_element_pt(19)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(15)->node_pt(np * n + np - 1);

      // Remove the nodes in element 19 from the temporaray node list
      delete tmp_node_pt[19 * np * np + np * (np - 1) + n];
      tmp_node_pt[19 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Elements 19 and 16
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 19 to be the same as in element 16
      this->finite_element_pt(16)->node_pt(np * n + np - 1) =
        this->finite_element_pt(19)->node_pt(n * np);

      // Remove the nodes in element 16 from the temporaray node list
      delete tmp_node_pt[16 * np * np + np * (np - 1) + n];
      tmp_node_pt[16 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Elements 15 and 16
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 16 to be the same as in element 15
      this->finite_element_pt(16)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(15)->node_pt(n);

      // Remove the nodes in element 16 from the temporaray node list
      delete tmp_node_pt[16 * np * np + np * (np - 1) + n];
      tmp_node_pt[16 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Element 18 and 20
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 20 to be the same as in element 18
      this->finite_element_pt(20)->node_pt(n) =
        this->finite_element_pt(18)->node_pt(np * (np - n - 1) + np - 1);

      // Remove the nodes in element 20 from the temporaray node list
      delete tmp_node_pt[20 * np * np + n];
      tmp_node_pt[20 * np * np + n] = 0;
    }


    // Edge between Elements 17 and 20
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 20 to be the same as in element 17
      this->finite_element_pt(20)->node_pt(n * np) =
        this->finite_element_pt(17)->node_pt(np * n + np - 1);

      // Remove the nodes in element 20 from the temporaray node list
      delete tmp_node_pt[20 * np * np + n * np];
      tmp_node_pt[20 * np * np + n * np] = 0;
    }


    // Edge between Elements 17 and 18
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 18 to be the same as in element 17
      this->finite_element_pt(18)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(17)->node_pt(n);

      // Remove the nodes in element 18 from the temporaray node list
      delete tmp_node_pt[18 * np * np + np * (np - 1) + n];
      tmp_node_pt[18 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Elements 19 and 21
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 21 to be the same as in element 19
      this->finite_element_pt(21)->node_pt(n * np) =
        this->finite_element_pt(19)->node_pt(np * n + np - 1);

      // Remove the nodes in element 21 from the temporaray node list
      delete tmp_node_pt[21 * np * np + n * np];
      tmp_node_pt[21 * np * np + n * np] = 0;
    }


    // Edge between Elements 21 and 22
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 22 to be the same as in element 21
      this->finite_element_pt(22)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(21)->node_pt(n);

      // Remove the nodes in element 22 from the temporaray node list
      delete tmp_node_pt[22 * np * np + np * (np - 1) + n];
      tmp_node_pt[22 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Elements 20 and 23
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 23 to be the same as in element 20
      this->finite_element_pt(23)->node_pt(n * np) =
        this->finite_element_pt(20)->node_pt(np * n + np - 1);

      // Remove the nodes in element 23 from the temporaray node list
      delete tmp_node_pt[23 * np * np + n * np];
      tmp_node_pt[23 * np * np + n * np] = 0;
    }


    // Edge between Elements 23 and 22
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 22 to be the same as in element 23
      this->finite_element_pt(22)->node_pt(n) =
        this->finite_element_pt(23)->node_pt(np * (np - 1) + n);

      // Remove the nodes in element 22 from the temporaray node list
      delete tmp_node_pt[22 * np * np + np * (np - 1) + n];
      tmp_node_pt[22 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Elements 21 and 24
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 24 to be the same as in element 21
      this->finite_element_pt(24)->node_pt(n * np) =
        this->finite_element_pt(21)->node_pt(np * n + np - 1);

      // Remove the nodes in element 24 from the temporaray node list
      delete tmp_node_pt[24 * np * np + n * np];
      tmp_node_pt[24 * np * np + n * np] = 0;
    }


    // Edge between Elements 22 and 25
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 25 to be the same as in element 22
      this->finite_element_pt(25)->node_pt(n * np) =
        this->finite_element_pt(22)->node_pt(np * n + np - 1);

      // Remove the nodes in element 25 from the temporaray node list
      delete tmp_node_pt[25 * np * np + n * np];
      tmp_node_pt[25 * np * np + n * np] = 0;
    }


    // Edge between Elements 23 and 26
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 26 to be the same as in element 23
      this->finite_element_pt(26)->node_pt(n * np) =
        this->finite_element_pt(23)->node_pt(np * n + np - 1);

      // Remove the nodes in element 26 from the temporaray node list
      delete tmp_node_pt[26 * np * np + n * np];
      tmp_node_pt[26 * np * np + n * np] = 0;
    }


    // Edge between Elements 24 and 25
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 25 to be the same as in element 24
      this->finite_element_pt(25)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(24)->node_pt(n);

      // Remove the nodes in element 25 from the temporaray node list
      delete tmp_node_pt[25 * np * np + np * (np - 1) + n];
      tmp_node_pt[25 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Elements 25 and 26
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 26 to be the same as in element 25
      this->finite_element_pt(26)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(25)->node_pt(n);

      // Remove the nodes in element 26 from the temporaray node list
      delete tmp_node_pt[26 * np * np + np * (np - 1) + n];
      tmp_node_pt[26 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Element 24 and 27
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 27 to be the same as in element 24
      this->finite_element_pt(27)->node_pt(np * (np - 1) + n) =
        this->finite_element_pt(24)->node_pt(np * n + np - 1);

      // Remove the nodes in element 27 from the temporaray node list
      delete tmp_node_pt[27 * np * np + np * (np - 1) + n];
      tmp_node_pt[27 * np * np + np * (np - 1) + n] = 0;
    }


    // Edge between Elements 25 and 27
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 27 to be the same as in element 25
      this->finite_element_pt(27)->node_pt(n * np) =
        this->finite_element_pt(25)->node_pt(np * n + np - 1);

      // Remove the nodes in element 27 from the temporaray node list
      delete tmp_node_pt[27 * np * np + n * np];
      tmp_node_pt[27 * np * np + n * np] = 0;
    }


    // Edge between Element 26 and 27
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 27 to be the same as in element 26
      this->finite_element_pt(27)->node_pt(n) =
        this->finite_element_pt(26)->node_pt(np * (np - n - 1) + np - 1);

      // Remove the nodes in element 27 from the temporaray node list
      delete tmp_node_pt[27 * np * np + n];
      tmp_node_pt[27 * np * np + n] = 0;
    }


    // Edge between Elements 27 and 28
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 28 to be the same as in element 27
      this->finite_element_pt(28)->node_pt(n * np) =
        this->finite_element_pt(27)->node_pt(np * n + np - 1);

      // Remove the nodes in element 28 from the temporaray node list
      delete tmp_node_pt[28 * np * np + n * np];
      tmp_node_pt[28 * np * np + n * np] = 0;
    }


    // Edge between Elements 28 and 29
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 29 to be the same as in element 28
      this->finite_element_pt(29)->node_pt(n * np) =
        this->finite_element_pt(28)->node_pt(np * n + np - 1);

      // Remove the nodes in element 29 from the temporaray node list
      delete tmp_node_pt[29 * np * np + n * np];
      tmp_node_pt[29 * np * np + n * np] = 0;
    }


    // Edge between Elements 29 and 30
    for (unsigned n = 0; n < np; n++)
    {
      // Set the nodes in element 30 to be the same as in element 29
      this->finite_element_pt(30)->node_pt(n * np) =
        this->finite_element_pt(29)->node_pt(np * n + np - 1);

      // Remove the nodes in element 30 from the temporaray node list
      delete tmp_node_pt[30 * np * np + n * np];
      tmp_node_pt[30 * np * np + n * np] = 0;
    }

    // Now set the actual true nodes
    for (unsigned long n = 0; n < node_count; n++)
    {
      if (tmp_node_pt[n] != 0)
      {
        Node_pt.push_back(tmp_node_pt[n]);
      }
    }

    // Finally set the nodes on the boundaries
    this->set_nboundary(8);

    for (unsigned n = 0; n < np; n++)
    {
      // Left hand side
      this->add_boundary_node(3, this->finite_element_pt(0)->node_pt(n * np));
      // Right hand side
      this->add_boundary_node(
        1, this->finite_element_pt(30)->node_pt(n * np + np - 1));

      // First part of lower boundary
      this->add_boundary_node(0, this->finite_element_pt(6)->node_pt(n));

      // First part of upper boundary
      this->add_boundary_node(
        2, this->finite_element_pt(1)->node_pt(np * (np - 1) + n));

      // First part of hole boundary
      this->add_boundary_node(4, this->finite_element_pt(3)->node_pt(np * n));

      // First part of lower flag
      this->add_boundary_node(
        5, this->finite_element_pt(4)->node_pt(np * (np - 1) + n));

      // First part of upper flag
      this->add_boundary_node(6, this->finite_element_pt(3)->node_pt(n));

      // Right part of flag
      this->add_boundary_node(7, this->finite_element_pt(22)->node_pt(n * np));
    }

    for (unsigned n = 1; n < np; n++)
    {
      // Second part of lower boundary
      this->add_boundary_node(0, this->finite_element_pt(10)->node_pt(n));

      // Second part of upper boundary
      this->add_boundary_node(
        2, this->finite_element_pt(7)->node_pt(np * (np - 1) + n));

      // Next part of lower flag
      this->add_boundary_node(
        5, this->finite_element_pt(9)->node_pt(np * (np - 1) + n));

      // Next part of upper flag
      this->add_boundary_node(6, this->finite_element_pt(8)->node_pt(n));
    }

    for (unsigned n = np - 2; n > 0; n--)
    {
      // Next part of hole
      this->add_boundary_node(4, this->finite_element_pt(2)->node_pt(n));
    }

    for (unsigned n = 1; n < np; n++)
    {
      // Next part of lower boundary
      this->add_boundary_node(0, this->finite_element_pt(14)->node_pt(n));

      // Next part of upper boundary
      this->add_boundary_node(
        2, this->finite_element_pt(11)->node_pt(np * (np - 1) + n));

      // Next part of lower flag
      this->add_boundary_node(
        5, this->finite_element_pt(13)->node_pt(np * (np - 1) + n));

      // Next part of upper flag
      this->add_boundary_node(6, this->finite_element_pt(12)->node_pt(n));
    }

    for (unsigned n = np - 1; n > 0; n--)
    {
      // Next part of hole
      this->add_boundary_node(4, this->finite_element_pt(1)->node_pt(n));
    }

    for (unsigned n = 1; n < np; n++)
    {
      // Next part of lower boundary
      this->add_boundary_node(0, this->finite_element_pt(18)->node_pt(n));
      // Next part of upper boundary
      this->add_boundary_node(
        2, this->finite_element_pt(15)->node_pt(np * (np - 1) + n));

      // Next part of lower flag
      this->add_boundary_node(
        5, this->finite_element_pt(17)->node_pt(np * (np - 1) + n));

      // Next part of upper flag
      this->add_boundary_node(6, this->finite_element_pt(16)->node_pt(n));
    }

    for (unsigned n = np - 1; n > 0; n--)
    {
      // Next part of hole
      this->add_boundary_node(
        4, this->finite_element_pt(0)->node_pt(n * np + np - 1));
    }


    for (unsigned n = 1; n < np; n++)
    {
      // Next part of lower boundary
      this->add_boundary_node(0, this->finite_element_pt(23)->node_pt(n));
      // Next part of upper boundary
      this->add_boundary_node(
        2, this->finite_element_pt(21)->node_pt(np * (np - 1) + n));

      // Next part of hole
      this->add_boundary_node(
        4, this->finite_element_pt(6)->node_pt(np * (np - 1) + n));

      // Next part of lower flag
      this->add_boundary_node(
        5, this->finite_element_pt(20)->node_pt(np * (np - 1) + n));

      // Next part of upper flag
      this->add_boundary_node(6, this->finite_element_pt(19)->node_pt(n));
    }

    for (unsigned n = 0; n < np; n++)
    {
      // Next part of hole
      this->add_boundary_node(
        4, this->finite_element_pt(6)->node_pt(np * (np - 1) + n));
    }


    for (unsigned n = 1; n < np; n++)
    {
      // Next part of lower boundary
      this->add_boundary_node(0, this->finite_element_pt(26)->node_pt(n));
      // Next part of upper boundary
      this->add_boundary_node(
        2, this->finite_element_pt(24)->node_pt(np * (np - 1) + n));

      // Next part of hole
      this->add_boundary_node(
        4, this->finite_element_pt(5)->node_pt(np * (np - 1) + n));
    }


    for (unsigned n = 1; n < np; n++)
    {
      // Next part of lower boundary
      this->add_boundary_node(0, this->finite_element_pt(28)->node_pt(n));
      // Next part of upper boundary
      this->add_boundary_node(
        2, this->finite_element_pt(28)->node_pt(np * (np - 1) + n));

      // Next part of hole
      this->add_boundary_node(4, this->finite_element_pt(4)->node_pt(np * n));
    }

    for (unsigned n = 1; n < np; n++)
    {
      // Next part of lower boundary
      this->add_boundary_node(0, this->finite_element_pt(29)->node_pt(n));
      // Next part of upper boundary
      this->add_boundary_node(
        2, this->finite_element_pt(29)->node_pt(np * (np - 1) + n));
    }

    for (unsigned n = 1; n < np; n++)
    {
      // Next part of lower boundary
      this->add_boundary_node(0, this->finite_element_pt(30)->node_pt(n));
      // Next part of upper boundary
      this->add_boundary_node(
        2, this->finite_element_pt(30)->node_pt(np * (np - 1) + n));
    }


    this->node_update();
    setup_boundary_element_info();

    // Set boundary coordinates on the flag

    // Vector of Lagrangian coordinates used as boundary coordinate
    Vector<double> zeta(1);

    // loop on nodes of boundary 5
    unsigned nnode = this->nboundary_node(5);
    for (unsigned k = 0; k < nnode; k++)
    {
      Node* nod_pt = this->boundary_node_pt(5, k);
      zeta[0] = double(k) * flag_length / double(nnode - 1);
      nod_pt->set_coordinates_on_boundary(5, zeta);
    }

    // loop on nodes of boundary 6
    nnode = this->nboundary_node(6);
    for (unsigned k = 0; k < nnode; k++)
    {
      Node* nod_pt = this->boundary_node_pt(6, k);
      zeta[0] = double(k) * flag_length / double(nnode - 1);
      nod_pt->set_coordinates_on_boundary(6, zeta);
    }

    // loop on nodes of boundary 7
    nnode = this->nboundary_node(7);
    for (unsigned k = 0; k < nnode; k++)
    {
      Node* nod_pt = this->boundary_node_pt(7, k);
      zeta[0] = -flag_height / 2. + double(k) / double(nnode - 1) * flag_height;
      nod_pt->set_coordinates_on_boundary(7, zeta);
    }

    // We have parametrised boundary 5,6 and 7
    this->Boundary_coordinate_exists[5] = true;
    this->Boundary_coordinate_exists[6] = true;
    this->Boundary_coordinate_exists[7] = true;

    // Loop over all elements and set macro element pointer
    for (unsigned e = 0; e < 31; e++)
    {
      dynamic_cast<ELEMENT*>(this->element_pt(e))
        ->set_macro_elem_pt(this->Domain_pt->macro_element_pt(e));
    }
  }


  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Setup algebraic node update
  //=================================================================
  template<class ELEMENT>
  void AlgebraicCylinderWithFlagMesh<ELEMENT>::setup_algebraic_node_update()
  {
    // The update function requires six parameters in some cases:
    Vector<double> ref_value(6);
    for (unsigned i = 0; i < 5; i++)
    {
      ref_value[i] = 0.0;
    }

    // Part I : macro elements 8,12,16
    for (unsigned k = 0; k < 3; k++)
    {
      FiniteElement* el_pt = this->finite_element_pt(8 + k * 4);
      unsigned nnode = el_pt->nnode();
      for (unsigned i = 0; i < nnode; i++)
      {
        // Get local coordinates
        Vector<double> coord_loc(2);
        el_pt->local_coordinate_of_node(i, coord_loc);

        // First reference value : horizontal fraction
        ref_value[0] = 0.5 * (coord_loc[0] + 1.0);

        // Second reference value : vertical fraction
        ref_value[1] = 0.5 * (coord_loc[1] + 1.0);

        // Third reference value : zeta coordinate on flag
        ref_value[2] = double(k + 1) / 5. * Flag_length +
                       ref_value[0] * 1. / 5. * Flag_length;

        // Sub-geomobject corresponding to the zeta coordinate on the flag
        GeomObject* geom_obj_pt;
        Vector<double> s(1);
        Vector<double> zeta(1);
        zeta[0] = ref_value[2];
        Top_flag_pt->locate_zeta(zeta, geom_obj_pt, s);

        // Create vector of geomobject for add_node_update_info()
        Vector<GeomObject*> geom_object_pt(1);
        geom_object_pt[0] = geom_obj_pt;

        // Fourth reference value : local coordinate in the sub geomobject
        ref_value[3] = s[0];

        // Fifth reference value : x coordinate
        ref_value[4] = el_pt->node_pt(i)->x(0);

        // Setup algebraic update for node: Pass update information
        // to AlgebraicNode:
        dynamic_cast<AlgebraicNode*>(el_pt->node_pt(i))
          ->add_node_update_info(1, // ID
                                 this, // mesh
                                 geom_object_pt, // vector of geom objects
                                 ref_value); // vector of ref. values
      }
    }

    // Part II : macro elements 9,13,17
    for (unsigned k = 0; k < 3; k++)
    {
      FiniteElement* el_pt = this->finite_element_pt(9 + k * 4);
      unsigned nnode = el_pt->nnode();
      for (unsigned i = 0; i < nnode; i++)
      {
        // Get local coordinates
        Vector<double> coord_loc(2);
        el_pt->local_coordinate_of_node(i, coord_loc);

        // First reference value : horizontal fraction
        ref_value[0] = 0.5 * (coord_loc[0] + 1.0);

        // Second reference value : vertical fraction
        ref_value[1] = 0.5 * (coord_loc[1] + 1.0);

        // Third reference value : zeta coordinate on flag
        ref_value[2] = double(k + 1) / 5. * Flag_length +
                       ref_value[0] * 1. / 5. * Flag_length;

        // Sub-geomobject corresponding to the zeta coordinate on the flag
        GeomObject* geom_obj_pt;
        Vector<double> s(1);
        Vector<double> zeta(1);
        zeta[0] = ref_value[2];
        Bottom_flag_pt->locate_zeta(zeta, geom_obj_pt, s);

        // Create vector of geomobject for add_node_update_info()
        Vector<GeomObject*> geom_object_pt(1);
        geom_object_pt[0] = geom_obj_pt;

        // Fourth reference value : local coordinate in the sub geomobject
        ref_value[3] = s[0];

        // Fifth reference value : x coordinate
        ref_value[4] = el_pt->node_pt(i)->x(0);

        // Setup algebraic update for node: Pass update information
        // to AlgebraicNode:
        dynamic_cast<AlgebraicNode*>(el_pt->node_pt(i))
          ->add_node_update_info(2, // ID
                                 this, // mesh
                                 geom_object_pt, // vector of geom objects
                                 ref_value); // vector of ref. values
      }
    }

    // Part III : macro element 22
    FiniteElement* el_pt = this->finite_element_pt(22);
    unsigned nnode = el_pt->nnode();
    for (unsigned i = 0; i < nnode; i++)
    {
      // Get local coordinates
      Vector<double> coord_loc(2);
      el_pt->local_coordinate_of_node(i, coord_loc);

      // First reference value : horizontal fraction
      ref_value[0] = 0.5 * (coord_loc[0] + 1.0);

      // Second reference value : vertical fraction
      ref_value[1] = 0.5 * (coord_loc[1] + 1.0);

      // Third reference value : zeta coordinate on flag
      ref_value[2] = coord_loc[1] * Flag_height / 2.;

      // Sub-geomobject corresponding to the zeta coordinate on the flag
      GeomObject* geom_obj_pt;
      Vector<double> s(1);
      Vector<double> zeta(1);
      zeta[0] = ref_value[2];
      Tip_flag_pt->locate_zeta(zeta, geom_obj_pt, s);

      // Create vector of geomobject for add_node_update_info()
      Vector<GeomObject*> geom_object_pt(1);
      geom_object_pt[0] = geom_obj_pt;

      // Fourth reference value : local coordinate in the sub geomobject
      ref_value[3] = s[0];

      // Setup algebraic update for node: Pass update information
      // to AlgebraicNode:
      dynamic_cast<AlgebraicNode*>(el_pt->node_pt(i))
        ->add_node_update_info(3, // ID
                               this, // mesh
                               geom_object_pt, // vector of geom objects
                               ref_value); // vector of ref. values
    }

    // Part IV : macro element 21
    el_pt = this->finite_element_pt(21);
    nnode = el_pt->nnode();
    for (unsigned i = 0; i < nnode; i++)
    {
      // Get local coordinates
      Vector<double> coord_loc(2);
      el_pt->local_coordinate_of_node(i, coord_loc);

      // First reference value : horizontal fraction
      ref_value[0] = 0.5 * (coord_loc[0] + 1.0);

      // Second reference value : vertical fraction
      ref_value[1] = 0.5 * (coord_loc[1] + 1.0);

      // Sub-geomobject corresponding to the tip of the Tip_flag
      GeomObject* geom_obj_pt;
      Vector<double> s(1);
      Vector<double> zeta(1);
      zeta[0] = Flag_height / 2.;
      Tip_flag_pt->locate_zeta(zeta, geom_obj_pt, s);

      // Create vector of geomobject for add_node_update_info()
      Vector<GeomObject*> geom_object_pt(1);
      geom_object_pt[0] = geom_obj_pt;

      // Third reference value : local coordinate in the sub geomobject
      ref_value[2] = s[0];

      // Setup algebraic update for node: Pass update information
      // to AlgebraicNode:
      dynamic_cast<AlgebraicNode*>(el_pt->node_pt(i))
        ->add_node_update_info(4, // ID
                               this, // mesh
                               geom_object_pt, // vector of geom objects
                               ref_value); // vector of ref. values
    }

    // Part V : macro element 23
    el_pt = this->finite_element_pt(23);
    nnode = el_pt->nnode();
    for (unsigned i = 0; i < nnode; i++)
    {
      // Get local coordinates
      Vector<double> coord_loc(2);
      el_pt->local_coordinate_of_node(i, coord_loc);

      // First reference value : horizontal fraction
      ref_value[0] = 0.5 * (coord_loc[0] + 1.0);

      // Second reference value : vertical fraction
      ref_value[1] = 0.5 * (coord_loc[1] + 1.0);

      // Sub-geomobject corresponding to the tip of the Tip_flag
      GeomObject* geom_obj_pt;
      Vector<double> s(1);
      Vector<double> zeta(1);
      zeta[0] = -Flag_height / 2.;
      Tip_flag_pt->locate_zeta(zeta, geom_obj_pt, s);

      // Create vector of geomobject for add_node_update_info()
      Vector<GeomObject*> geom_object_pt(1);
      geom_object_pt[0] = geom_obj_pt;

      // Third reference value : local coordinate in the sub geomobject
      ref_value[2] = s[0];

      // Setup algebraic update for node: Pass update information
      // to AlgebraicNode:
      dynamic_cast<AlgebraicNode*>(el_pt->node_pt(i))
        ->add_node_update_info(5, // ID
                               this, // mesh
                               geom_object_pt, // vector of geom objects
                               ref_value); // vector of ref. values
    }

    // Part VI = macro element 19
    el_pt = this->finite_element_pt(19);
    nnode = el_pt->nnode();
    for (unsigned i = 0; i < nnode; i++)
    {
      // Get local coordinates
      Vector<double> coord_loc(2);
      el_pt->local_coordinate_of_node(i, coord_loc);

      // First reference value : horizontal fraction
      ref_value[0] = 0.5 * (coord_loc[0] + 1.0);

      // Second reference value : vertical fraction
      ref_value[1] = 0.5 * (coord_loc[1] + 1.0);

      // Third reference value : zeta coordinate on the flag
      ref_value[2] =
        4. / 5. * Flag_length + ref_value[0] * 1. / 5. * Flag_length;

      // Sub-geomobject
      GeomObject* geom_obj_pt;
      Vector<double> s(1);
      Vector<double> zeta(1);
      zeta[0] = ref_value[2];
      Top_flag_pt->locate_zeta(zeta, geom_obj_pt, s);

      // Create vector of geomobject for add_node_update_info()
      Vector<GeomObject*> geom_object_pt(1);
      geom_object_pt[0] = geom_obj_pt;

      // Third reference value : local coordinate in the sub geomobject
      ref_value[3] = s[0];

      // Setup algebraic update for node: Pass update information
      // to AlgebraicNode:
      dynamic_cast<AlgebraicNode*>(el_pt->node_pt(i))
        ->add_node_update_info(6, // ID
                               this, // mesh
                               geom_object_pt, // vector of geom objects
                               ref_value); // vector of ref. values
    }


    // Part VII = macro element 20
    el_pt = this->finite_element_pt(20);
    nnode = el_pt->nnode();
    for (unsigned i = 0; i < nnode; i++)
    {
      // Get local coordinates
      Vector<double> coord_loc(2);
      el_pt->local_coordinate_of_node(i, coord_loc);

      // First reference value : horizontal fraction
      ref_value[0] = 0.5 * (coord_loc[0] + 1.0);

      // Second reference value : vertical fraction
      ref_value[1] = 0.5 * (coord_loc[1] + 1.0);

      // Third reference value : zeta coordinate on the flag
      ref_value[2] =
        4. / 5. * Flag_length + ref_value[0] * 1. / 5. * Flag_length;

      // Sub-geomobject
      GeomObject* geom_obj_pt;
      Vector<double> s(1);
      Vector<double> zeta(1);
      zeta[0] = ref_value[2];
      Bottom_flag_pt->locate_zeta(zeta, geom_obj_pt, s);

      // Create vector of geomobject for add_node_update_info()
      Vector<GeomObject*> geom_object_pt(1);
      geom_object_pt[0] = geom_obj_pt;

      // Third reference value : local coordinate in the sub geomobject
      ref_value[3] = s[0];

      // Setup algebraic update for node: Pass update information
      // to AlgebraicNode:
      dynamic_cast<AlgebraicNode*>(el_pt->node_pt(i))
        ->add_node_update_info(7, // ID
                               this, // mesh
                               geom_object_pt, // vector of geom objects
                               ref_value); // vector of ref. values
    }

    // Part VIII : macro element 3
    el_pt = this->finite_element_pt(3);
    nnode = el_pt->nnode();
    for (unsigned i = 0; i < nnode; i++)
    {
      // Get local coordinates
      Vector<double> coord_loc(2);
      el_pt->local_coordinate_of_node(i, coord_loc);

      // First reference value : horizontal fraction
      ref_value[0] = 0.5 * (coord_loc[0] + 1.0);

      // Second reference value : vertical fraction
      ref_value[1] = 0.5 * (coord_loc[1] + 1.0);

      // Third reference value : zeta coordinate on flag at reference point
      ref_value[2] = ref_value[0] * 1. / 5. * Flag_length;

      // Sub-geomobject corresponding to the zeta coordinate on the flag
      // at the reference point
      GeomObject* geom_obj_pt;
      Vector<double> s(1);
      Vector<double> zeta(1);
      zeta[0] = ref_value[2];
      Top_flag_pt->locate_zeta(zeta, geom_obj_pt, s);

      // Fourth reference value : local coordinate in the sub geomobject
      ref_value[3] = s[0];

      // Create vector of geomobject for add_node_update_info()
      Vector<GeomObject*> geom_object_pt(2);
      geom_object_pt[0] = geom_obj_pt;

      // Fifth reference value : zeta coordinate on flag at end of macro element
      ref_value[4] = 1. / 5. * Flag_length;

      // Sub-geomobject corresponding to the zeta coordinate on the flag
      // at the end of the macro element
      zeta[0] = ref_value[4];
      Top_flag_pt->locate_zeta(zeta, geom_obj_pt, s);

      // Add geom object
      geom_object_pt[1] = geom_obj_pt;

      // Sixth reference value : local coordinate in the sub geomobject
      ref_value[5] = s[0];


      // Setup algebraic update for node: Pass update information
      // to AlgebraicNode:
      dynamic_cast<AlgebraicNode*>(el_pt->node_pt(i))
        ->add_node_update_info(8, // ID
                               this, // mesh
                               geom_object_pt, // vector of geom objects
                               ref_value); // vector of ref. values
    }


    // Part IX : macro element 4
    el_pt = this->finite_element_pt(4);
    nnode = el_pt->nnode();
    for (unsigned i = 0; i < nnode; i++)
    {
      // Get local coordinates
      Vector<double> coord_loc(2); /**set the size ??*/
      el_pt->local_coordinate_of_node(i, coord_loc);

      // First reference value : horizontal fraction
      ref_value[0] = 0.5 * (coord_loc[0] + 1.0);

      // Second reference value : vertical fraction
      ref_value[1] = 0.5 * (coord_loc[1] + 1.0);

      // Third reference value : zeta coordinate on flag
      ref_value[2] = ref_value[0] * 1. / 5. * Flag_length;

      // Sub-geomobject corresponding to the zeta coordinate on the flag
      // at the reference point
      GeomObject* geom_obj_pt;
      Vector<double> s(1);
      Vector<double> zeta(1);
      zeta[0] = ref_value[2];
      Bottom_flag_pt->locate_zeta(zeta, geom_obj_pt, s);

      // Fourth reference value : local coordinate in the sub geomobject
      ref_value[3] = s[0];

      // Create vector of geomobject for add_node_update_info()
      Vector<GeomObject*> geom_object_pt(2);
      geom_object_pt[0] = geom_obj_pt;

      // Fifth reference value : zeta coordinate on flag at end of macro element
      ref_value[4] = 1. / 5. * Flag_length;

      // Sub-geomobject corresponding to the zeta coordinate on the flag
      // at the end of the macro element
      zeta[0] = ref_value[4];
      Bottom_flag_pt->locate_zeta(zeta, geom_obj_pt, s);

      // Add geom object
      geom_object_pt[1] = geom_obj_pt;

      // Sixth reference value : local coordinate in the sub geomobject
      ref_value[5] = s[0];

      // Setup algebraic update for node: Pass update information
      // to AlgebraicNode:
      dynamic_cast<AlgebraicNode*>(el_pt->node_pt(i))
        ->add_node_update_info(9, // ID
                               this, // mesh
                               geom_object_pt, // vector of geom objects
                               ref_value); // vector of ref. values
    }


  } // end of setup_algebraic_node_update


  //=================================================================
  /// The algebraic node update function
  //=================================================================
  template<class ELEMENT>
  void AlgebraicCylinderWithFlagMesh<ELEMENT>::algebraic_node_update(
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

      case 5:
        node_update_V(t, node_pt);
        break;

      case 6:
        node_update_VI(t, node_pt);
        break;

      case 7:
        node_update_VII(t, node_pt);
        break;

      case 8:
        node_update_VIII(t, node_pt);
        break;

      case 9:
        node_update_IX(t, node_pt);
        break;

      default:
        std::ostringstream error_message;
        error_message << "Wrong id " << id << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }


  } // end of algebraic_node_update()

  //=================================================================
  /// Node update for region I
  //=================================================================
  template<class ELEMENT>
  void AlgebraicCylinderWithFlagMesh<ELEMENT>::node_update_I(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to geom object
    GeomObject* flag_pt = geom_object_pt[0];

    // Point on the line y=p11[1] corresponding to the initial x.
    Vector<double> ref_point(2);
    ref_point[0] = ref_value[4];
    ref_point[1] = 0.778024390 * Height;

    // Point on the flag
    Vector<double> flag_point(2);
    Vector<double> zeta(1);
    zeta[0] = ref_value[3];
    flag_pt->position(t, zeta, flag_point);

    // Third reference value : fraction of the vertical line
    // between the straight line y = p11[1] and the flag
    double r = ref_value[1];

    // Assign new nodal coordinates
    node_pt->x(t, 0) =
      ref_point[0] + (1.0 - r) * (flag_point[0] - ref_point[0]);
    node_pt->x(t, 1) =
      ref_point[1] + (1.0 - r) * (flag_point[1] - ref_point[1]);
  }


  //=================================================================
  /// Node update for region II
  //=================================================================
  template<class ELEMENT>
  void AlgebraicCylinderWithFlagMesh<ELEMENT>::node_update_II(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to geom object
    GeomObject* flag_pt = geom_object_pt[0];

    // Point on the line y=p37[1] corresponding to the initial x.
    Vector<double> ref_point(2);
    ref_point[0] = ref_value[4];
    ref_point[1] = 0.197585366 * Height;

    // Point on the flag
    Vector<double> flag_point(2);
    Vector<double> zeta(1);
    zeta[0] = ref_value[3];
    flag_pt->position(t, zeta, flag_point);

    // Third reference value : fraction of the vertical line
    // between the straight line y = p11[1] and the flag
    double r = ref_value[1];

    // Assign new nodal coordinates
    node_pt->x(t, 0) = ref_point[0] + r * (flag_point[0] - ref_point[0]);
    node_pt->x(t, 1) = ref_point[1] + r * (flag_point[1] - ref_point[1]);
  }

  //=================================================================
  /// Node update for region III
  //=================================================================
  template<class ELEMENT>
  void AlgebraicCylinderWithFlagMesh<ELEMENT>::node_update_III(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // useful points
    Vector<double> p15(2);
    Vector<double> p35(2);

    p15[0] = 0.285123967 * Length;
    p15[1] = 0.625 * Height;

    p35[0] = 0.285123967 * Length;
    p35[1] = 0.350609756 * Height;

    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to geom object
    GeomObject* flag_pt = geom_object_pt[0];

    // Point on the line x=p15[0]
    Vector<double> ref_point(2);
    ref_point[0] = p15[0];
    ref_point[1] = p35[1] + ref_value[1] * (p15[1] - p35[1]);

    // Point on the flag
    Vector<double> flag_point(2);
    Vector<double> zeta(1);
    zeta[0] = ref_value[3];
    flag_pt->position(t, zeta, flag_point);

    // Third reference value : fraction of the horizontal line
    // between the flag and the horizontal straight line in x=p15[0]
    double r = ref_value[0];

    // Assign new nodal coordinates
    node_pt->x(t, 0) = flag_point[0] + r * (ref_point[0] - flag_point[0]);
    node_pt->x(t, 1) = flag_point[1] + r * (ref_point[1] - flag_point[1]);
  }

  //=================================================================
  /// Node update for region IV
  //=================================================================
  template<class ELEMENT>
  void AlgebraicCylinderWithFlagMesh<ELEMENT>::node_update_IV(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // Useful points
    Vector<double> p15(2);
    Vector<double> p25(2);
    Vector<double> top_flag(2);

    p15[0] = 0.285123967 * Length;
    p15[1] = 0.625 * Height;

    p25[0] = Centre_x +
             A * sqrt(1.0 - Flag_height * Flag_height / (4.0 * A * A)) +
             Flag_length;
    p25[1] = Centre_y + Flag_height / 2.0;

    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to geom object
    GeomObject* flag_pt = geom_object_pt[0];

    Vector<double> zeta(1);
    zeta[0] = ref_value[2];
    flag_pt->position(t, zeta, top_flag);

    // point on the line linking p15 et top_flag
    Vector<double> p1(2);
    p1[0] = top_flag[0] + ref_value[0] * (p15[0] - top_flag[0]);
    p1[1] = top_flag[1] + ref_value[0] * (p15[1] - top_flag[1]);

    // Point on the line y = Height;
    Vector<double> p2(2);
    p2[0] = p25[0] + ref_value[0] * (p15[0] - p25[0]);
    p2[1] = Height;

    // Connect those points with the vertical fraction ref_value[1]
    node_pt->x(t, 0) = p1[0] + ref_value[1] * (p2[0] - p1[0]);
    node_pt->x(t, 1) = p1[1] + ref_value[1] * (p2[1] - p1[1]);
  }

  //=================================================================
  /// Node update for region V
  //=================================================================
  template<class ELEMENT>
  void AlgebraicCylinderWithFlagMesh<ELEMENT>::node_update_V(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // Useful points
    Vector<double> p31(2);
    Vector<double> p35(2);
    Vector<double> top_flag(2);

    p31[0] = Centre_x +
             A * sqrt(1.0 - Flag_height * Flag_height / (4.0 * A * A)) +
             Flag_length;
    p31[1] = Centre_y - Flag_height / 2.;

    p35[0] = 0.285123967 * Length;
    p35[1] = 0.350609756 * Height;

    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to geom object
    GeomObject* flag_pt = geom_object_pt[0];

    Vector<double> zeta(1);
    zeta[0] = ref_value[2];

    flag_pt->position(t, zeta, top_flag);

    // point on the line linking p35 et top_flag
    Vector<double> p1(2);
    p1[0] = top_flag[0] + ref_value[0] * (p35[0] - top_flag[0]);
    p1[1] = top_flag[1] + ref_value[0] * (p35[1] - top_flag[1]);

    // Point on the line y = 0.0;
    Vector<double> p2(2);
    p2[0] = p31[0] + ref_value[0] * (p35[0] - p31[0]);
    p2[1] = 0.;

    // Connect those points with the vertical fraction ref_value[1]
    node_pt->x(t, 0) = p2[0] + ref_value[1] * (p1[0] - p2[0]);
    node_pt->x(t, 1) = p2[1] + ref_value[1] * (p1[1] - p2[1]);
  }

  //=================================================================
  /// Node update for region VI
  //=================================================================
  template<class ELEMENT>
  void AlgebraicCylinderWithFlagMesh<ELEMENT>::node_update_VI(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Useful points
    Vector<double> p14(2);
    Vector<double> p5(2);
    Vector<double> point_flag(2);

    p14[0] = 0.211596 * Length;
    p14[1] = 0.778024390 * Height;

    p5[0] = 0.239596 * Length;
    p5[1] = Height;

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to geom object
    GeomObject* flag_pt = geom_object_pt[0];

    Vector<double> zeta(1);
    zeta[0] = ref_value[3];
    flag_pt->position(t, zeta, point_flag);

    // point on the line linking p14 et p5
    Vector<double> p1(2);
    p1[0] = p14[0] + ref_value[0] * (p5[0] - p14[0]);
    p1[1] = p14[1] + ref_value[0] * (p5[1] - p14[1]);


    // Connect those points with the vertical fraction ref_value[1]
    node_pt->x(t, 0) = point_flag[0] + ref_value[1] * (p1[0] - point_flag[0]);
    node_pt->x(t, 1) = point_flag[1] + ref_value[1] * (p1[1] - point_flag[1]);
  }

  //=================================================================
  /// Node update for region VII
  //=================================================================
  template<class ELEMENT>
  void AlgebraicCylinderWithFlagMesh<ELEMENT>::node_update_VII(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Useful points
    Vector<double> p40(2);
    Vector<double> p45(2);
    Vector<double> point_flag(2);

    p40[0] = 0.211596 * Length;
    p40[1] = 0.197585366 * Height;

    p45[0] = 0.239596 * Length;
    p45[1] = 0.0;

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to geom object
    GeomObject* flag_pt = geom_object_pt[0];

    Vector<double> zeta(1);
    zeta[0] = ref_value[3];
    flag_pt->position(t, zeta, point_flag);

    // point on the line linking p40 et p45
    Vector<double> p1(2);
    p1[0] = p40[0] + ref_value[0] * (p45[0] - p40[0]);
    p1[1] = p40[1] + ref_value[0] * (p45[1] - p40[1]);


    // Connect those points with the vertical fraction ref_value[1]
    node_pt->x(t, 0) =
      point_flag[0] + (1 - ref_value[1]) * (p1[0] - point_flag[0]);
    node_pt->x(t, 1) =
      point_flag[1] + (1 - ref_value[1]) * (p1[1] - point_flag[1]);
  }


  //=================================================================
  /// Node update for region VIII
  //=================================================================
  template<class ELEMENT>
  void AlgebraicCylinderWithFlagMesh<ELEMENT>::node_update_VIII(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // Useful point
    Vector<double> p11(2);
    p11[0] = 0.127596 * Length;
    p11[1] = 0.778024390 * Height;

    /// Extreme angles on circle
    double zeta_circle_top = atan(1.0);
    double zeta_circle_bot = asin(Flag_height / 2. / A);

    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to geom object containing the reference point
    GeomObject* flag_ref_pt = geom_object_pt[0];

    // Pointer to geom object containing the end of the macro element
    GeomObject* flag_end_pt = geom_object_pt[1];

    double omega_horizontal = ref_value[0];
    double omega_vertical = ref_value[1];

    // end of the macro element on the flag
    Vector<double> flag_end(2);
    Vector<double> zeta(1);
    zeta[0] = ref_value[5];
    flag_end_pt->position(t, zeta, flag_end);

    Vector<double> outer_point(p11);

    // Get reference point on circle
    Vector<double> ref_point_on_circle(2);
    zeta[0] =
      zeta_circle_bot + (zeta_circle_top - zeta_circle_bot) * omega_vertical;
    Cylinder_pt->position(zeta, ref_point_on_circle);

    // Get reference point on right line
    Vector<double> ref_point_on_right_line(2);
    ref_point_on_right_line[0] = // flag_end[0];
      outer_point[0] + (flag_end[0] - outer_point[0]) * (1.0 - omega_vertical);
    ref_point_on_right_line[1] =
      outer_point[1] + (flag_end[1] - outer_point[1]) * (1.0 - omega_vertical);

    // Get reference point on flag
    Vector<double> ref_point_on_flag(2);
    zeta[0] = ref_value[3];
    flag_ref_pt->position(t, zeta, ref_point_on_flag);

    // Get bottom-most point on circle
    Vector<double> circle_bot(2);
    zeta[0] = zeta_circle_bot;
    Cylinder_pt->position(zeta, circle_bot);

    // Get reference point on horizontal fraction of straight line
    // connecting the two bottom most reference points
    Vector<double> r_bot(2);
    r_bot[0] = circle_bot[0] + (flag_end[0] - circle_bot[0]) * omega_horizontal;
    r_bot[1] = circle_bot[1] + (flag_end[1] - circle_bot[1]) * omega_horizontal;

    // Place point on horizontal fraction of straight line
    // connecting reference points -- this won't match the
    // curved top boundary adjacent to the flag
    node_pt->x(t, 0) =
      ref_point_on_circle[0] +
      (ref_point_on_right_line[0] - ref_point_on_circle[0]) * omega_horizontal;
    node_pt->x(t, 1) =
      ref_point_on_circle[1] +
      (ref_point_on_right_line[1] - ref_point_on_circle[1]) * omega_horizontal;

    // Correct by scaled difference between bottom straight line
    // and bent flag
    node_pt->x(t, 0) +=
      (ref_point_on_flag[0] - r_bot[0]) * (1.0 - omega_vertical);
    node_pt->x(t, 1) +=
      (ref_point_on_flag[1] - r_bot[1]) * (1.0 - omega_vertical);
  }

  //=================================================================
  /// Node update for region IX
  //=================================================================
  template<class ELEMENT>
  void AlgebraicCylinderWithFlagMesh<ELEMENT>::node_update_IX(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
    // Useful point
    Vector<double> p37(2);
    p37[0] = 0.127596 * Length;
    p37[1] = 0.197585366 * Height;

    /// Extreme angles on circle
    double zeta_circle_top = -asin(Flag_height / 2. / A);
    double zeta_circle_bot = -atan(1.0);

    // Extract reference values for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to geom object containing the reference point
    GeomObject* flag_ref_pt = geom_object_pt[0];

    // Pointer to geom object containing the end of macro element
    GeomObject* flag_end_pt = geom_object_pt[1];

    double omega_horizontal = ref_value[0];
    double omega_vertical = ref_value[1];

    // end of the macro element on the flag
    Vector<double> flag_end(2);
    Vector<double> zeta(1);
    zeta[0] = ref_value[5];
    flag_end_pt->position(t, zeta, flag_end);

    Vector<double> outer_point(p37);

    // Get reference point on circle
    Vector<double> ref_point_on_circle(2);
    zeta[0] =
      zeta_circle_bot + (zeta_circle_top - zeta_circle_bot) * omega_vertical;
    Cylinder_pt->position(zeta, ref_point_on_circle);

    // Get reference point on right line
    Vector<double> ref_point_on_right_line(2);
    ref_point_on_right_line[0] = // flag_end[0];
      outer_point[0] + (flag_end[0] - outer_point[0]) * omega_vertical;
    ref_point_on_right_line[1] =
      outer_point[1] + (flag_end[1] - outer_point[1]) * omega_vertical;

    // Get reference point on flag
    Vector<double> ref_point_on_flag(2);
    zeta[0] = ref_value[3];
    flag_ref_pt->position(t, zeta, ref_point_on_flag);

    // Get top-most point on circle
    Vector<double> circle_top(2);
    zeta[0] = zeta_circle_top;
    Cylinder_pt->position(zeta, circle_top);

    // Get reference point on horizontal fraction of straight line
    // connecting the two top most reference points
    Vector<double> r_top(2);
    r_top[0] = circle_top[0] + (flag_end[0] - circle_top[0]) * omega_horizontal;
    r_top[1] = circle_top[1] + (flag_end[1] - circle_top[1]) * omega_horizontal;

    // Place point on horizontal fraction of straight line
    // connecting reference points -- this won't match the
    // curved top boundary adjacent to the flag
    node_pt->x(t, 0) =
      ref_point_on_circle[0] +
      (ref_point_on_right_line[0] - ref_point_on_circle[0]) * omega_horizontal;
    node_pt->x(t, 1) =
      ref_point_on_circle[1] +
      (ref_point_on_right_line[1] - ref_point_on_circle[1]) * omega_horizontal;

    // Correct by scaled difference between top straight line
    // and bent flag
    node_pt->x(t, 0) += (ref_point_on_flag[0] - r_top[0]) * omega_vertical;
    node_pt->x(t, 1) += (ref_point_on_flag[1] - r_top[1]) * omega_vertical;
  }


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


  //===================================================================
  ///  Update the node update functions
  //===================================================================
  template<class ELEMENT>
  void RefineableAlgebraicCylinderWithFlagMesh<ELEMENT>::update_node_update(
    AlgebraicNode*& node_pt)
  {
    // Extract ID
    unsigned id = node_pt->node_update_fct_id();


    if (id == 8)
    {
      // Extract geometric objects for update by copy construction
      Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

      // Extract reference values for node update by copy construction
      Vector<double> ref_value(node_pt->vector_ref_value());

      // Get zeta coordinate of reference point
      Vector<double> zeta_ref_flag(1);
      zeta_ref_flag[0] = ref_value[2];

      // Get the sub-geomobject and the local coordinate containing the
      // reference point
      Vector<double> s(1);
      GeomObject* geom_obj_pt;
      this->Top_flag_pt->locate_zeta(zeta_ref_flag, geom_obj_pt, s);


      // Update the pointer to the (sub-)GeomObject within which the
      // reference point is located.
      geom_object_pt[0] = geom_obj_pt;

      // Update second reference value: Reference local coordinate
      // in flag sub-element
      ref_value[3] = s[0];


      // Get zeta coordinate of point at end of macro element
      Vector<double> zeta_end_flag(1);
      zeta_end_flag[0] = ref_value[4];

      // Get the sub-geomobject and the local coordinate containing the
      // point at the end of the macro element
      this->Top_flag_pt->locate_zeta(zeta_end_flag, geom_obj_pt, s);


      // Update the pointer to the (sub-)GeomObject within which the
      // point at the end of the macro element is located
      geom_object_pt[1] = geom_obj_pt;

      // Update second reference value: Reference local coordinate
      // in flag sub-element
      ref_value[5] = s[0];

      // Kill the existing node update info
      node_pt->kill_node_update_info(8);

      // Setup algebraic update for node: Pass update information
      node_pt->add_node_update_info(8, // id
                                    this, // mesh
                                    geom_object_pt, // vector of geom objects
                                    ref_value); // vector of ref. values
    }
    else if (id == 9)
    {
      // Extract geometric objects for update by copy construction
      Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

      // Extract reference values for node update by copy construction
      Vector<double> ref_value(node_pt->vector_ref_value());

      // Get zeta coordinate of reference point
      Vector<double> zeta_ref_flag(1);
      zeta_ref_flag[0] = ref_value[2];

      // Get the sub-geomobject and the local coordinate containing the
      // reference point
      Vector<double> s(1);
      GeomObject* geom_obj_pt;
      this->Bottom_flag_pt->locate_zeta(zeta_ref_flag, geom_obj_pt, s);

      // Update the pointer to the (sub-)GeomObject within which the
      // reference point is located.
      geom_object_pt[0] = geom_obj_pt;

      // Update second reference value: Reference local coordinate
      // in flag sub-element
      ref_value[3] = s[0];


      // Get zeta coordinate of point at end of macro element
      Vector<double> zeta_end_flag(1);
      zeta_end_flag[0] = ref_value[4];

      // Get the sub-geomobject and the local coordinate containing the
      // point at the end of the macro element
      this->Bottom_flag_pt->locate_zeta(zeta_end_flag, geom_obj_pt, s);

      // Update the pointer to the (sub-)GeomObject within which the
      // point at the end of the macro element is located
      geom_object_pt[1] = geom_obj_pt;

      // Update second reference value: Reference local coordinate
      // in flag sub-element
      ref_value[5] = s[0];

      // Kill the existing node update info
      node_pt->kill_node_update_info(9);

      // Setup algebraic update for node: Pass update information
      node_pt->add_node_update_info(9, // id
                                    this, // mesh
                                    geom_object_pt, // vector of geom objects
                                    ref_value); // vector of ref. values
    }


    if ((id == 1) || (id == 6))
    {
      // Extract reference values for node update by copy construction
      Vector<double> ref_value(node_pt->vector_ref_value());

      // Get zeta coordinate on flag
      Vector<double> zeta_flag(1);
      zeta_flag[0] = ref_value[2];

      // Get the sub-geomobject and the local coordinate
      Vector<double> s(1);
      GeomObject* geom_obj_pt;
      this->Top_flag_pt->locate_zeta(zeta_flag, geom_obj_pt, s);

      // Extract geometric objects for update by copy construction
      Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

      // Update the pointer to the (sub-)GeomObject within which the
      // reference point is located. (If the flag is simple GeomObject
      // this is the same as Leaflet_pt; if it's a compound GeomObject
      // this points to the sub-object)
      geom_object_pt[0] = geom_obj_pt;

      // Update second reference value: Reference local coordinate
      // in flag sub-element
      ref_value[3] = s[0];

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
      else if (id == 6)
      {
        // Kill the existing node update info
        node_pt->kill_node_update_info(6);

        // Setup algebraic update for node: Pass update information
        node_pt->add_node_update_info(6, // id
                                      this, // mesh
                                      geom_object_pt, // vector of geom objects
                                      ref_value); // vector of ref. values
      }
    }


    if ((id == 2) || (id == 7))
    {
      // Extract reference values for node update by copy construction
      Vector<double> ref_value(node_pt->vector_ref_value());

      // Get zeta coordinate on flag
      Vector<double> zeta_flag(1);
      zeta_flag[0] = ref_value[2];

      // Get the sub-geomobject and the local coordinate
      Vector<double> s(1);
      GeomObject* geom_obj_pt;
      this->Bottom_flag_pt->locate_zeta(zeta_flag, geom_obj_pt, s);

      // Extract geometric objects for update by copy construction
      Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

      // Update the pointer to the (sub-)GeomObject within which the
      // reference point is located. (If the flag is simple GeomObject
      // this is the same as Leaflet_pt; if it's a compound GeomObject
      // this points to the sub-object)
      geom_object_pt[0] = geom_obj_pt;

      // Update second reference value: Reference local coordinate
      // in flag sub-element
      ref_value[3] = s[0];

      if (id == 2)
      {
        // Kill the existing node update info
        node_pt->kill_node_update_info(2);

        // Setup algebraic update for node: Pass update information
        node_pt->add_node_update_info(2, // id
                                      this, // mesh
                                      geom_object_pt, // vector of geom objects
                                      ref_value); // vector of ref. values
      }
      else if (id == 7)
      {
        // Kill the existing node update info
        node_pt->kill_node_update_info(7);

        // Setup algebraic update for node: Pass update information
        node_pt->add_node_update_info(7, // id
                                      this, // mesh
                                      geom_object_pt, // vector of geom objects
                                      ref_value); // vector of ref. values
      }
    }

    if ((id == 3) || (id == 4) || (id == 5))
    {
      // Extract reference values for node update by copy construction
      Vector<double> ref_value(node_pt->vector_ref_value());

      // Get zeta coordinate on flag
      Vector<double> zeta_flag(1);
      if (id == 3)
      {
        zeta_flag[0] = ref_value[2];
      }
      else if (id == 4)
      {
        zeta_flag[0] = this->Flag_height / 2.;
      }
      else if (id == 5)
      {
        zeta_flag[0] = -this->Flag_height / 2.;
      }

      // Get the sub-geomobject and the local coordinate
      Vector<double> s(1);
      GeomObject* geom_obj_pt;
      this->Tip_flag_pt->locate_zeta(zeta_flag, geom_obj_pt, s);

      // Extract geometric objects for update by copy construction
      Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

      // Update the pointer to the (sub-)GeomObject within which the
      // reference point is located. (If the flag is simple GeomObject
      // this is the same as Leaflet_pt; if it's a compound GeomObject
      // this points to the sub-object)
      geom_object_pt[0] = geom_obj_pt;


      if (id == 3)
      {
        // Update second reference value: Reference local coordinate
        // in flag sub-element
        ref_value[3] = s[0];

        // Kill the existing node update info
        node_pt->kill_node_update_info(3);

        // Setup algebraic update for node: Pass update information
        node_pt->add_node_update_info(3, // id
                                      this, // mesh
                                      geom_object_pt, // vector of geom objects
                                      ref_value); // vector of ref. values
      }
      else if (id == 4)
      {
        // Update second reference value: Reference local coordinate
        // in flag sub-element
        ref_value[2] = s[0];

        // Kill the existing node update info
        node_pt->kill_node_update_info(4);

        // Setup algebraic update for node: Pass update information
        node_pt->add_node_update_info(4, // id
                                      this, // mesh
                                      geom_object_pt, // vector of geom objects
                                      ref_value); // vector of ref. values
      }
      else if (id == 5)
      {
        // Update second reference value: Reference local coordinate
        // in flag sub-element
        ref_value[2] = s[0];

        // Kill the existing node update info
        node_pt->kill_node_update_info(5);

        // Setup algebraic update for node: Pass update information
        node_pt->add_node_update_info(5, // id
                                      this, // mesh
                                      geom_object_pt, // vector of geom objects
                                      ref_value); // vector of ref. values
      }
    }
  }


} // namespace oomph

#endif
