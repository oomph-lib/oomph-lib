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
#include "geompack_scaffold_mesh.h"

namespace oomph
{
  //=====================================================================
  /// \short Constructor: Pass the filename of the mesh files
  //=====================================================================
  GeompackQuadScaffoldMesh::GeompackQuadScaffoldMesh(
    const std::string& mesh_file_name, const std::string& curve_file_name)
  {
    // Read the output mesh file to find informations about the nodes
    // and elements of the mesh

    // Process mesh file
    //---------------------
    std::ifstream mesh_file(mesh_file_name.c_str(), std::ios_base::in);

    // Number of nodes
    unsigned n_node;
    mesh_file >> n_node;

    // Coordinates of nodes and extra information index
    Vector<double> x_node(n_node);
    Vector<double> y_node(n_node);
    Vector<int> vertinfo(n_node);
    for (unsigned i = 0; i < n_node; i++)
    {
      mesh_file >> x_node[i];
      mesh_file >> y_node[i];
      mesh_file >> vertinfo[i];
    }

    // Extra information (nodes lying on a curve)
    unsigned n_vx;
    mesh_file >> n_vx;

    // Dummy storage for the node code
    int dummy_node_code;

    // Storage for the curve indice
    Vector<int> icurv(n_vx); // it's negative if not used

    // Dummy storage for the location of the point on the curve
    double dummy_ucurv;

    for (unsigned i = 0; i < n_vx; i++)
    {
      mesh_file >> dummy_node_code;
      mesh_file >> icurv[i];
      mesh_file >> dummy_ucurv;
    }

    // Number of local nodes per element
    unsigned n_local_node;
    mesh_file >> n_local_node;

    // Number of elements
    unsigned n_element;
    mesh_file >> n_element;

    // Storage for global node numbers for all elements
    Vector<unsigned> global_node(n_local_node * n_element);

    // Storage for edge information
    // (needed for a possible construction of midside node
    // in the following build from scaffold function)
    Vector<int> edgeinfo(n_local_node * n_element);

    // Initialize counter
    unsigned k = 0;

    // Read global node numbers for all elements
    for (unsigned i = 0; i < n_element; i++)
    {
      for (unsigned j = 0; j < n_local_node; j++)
      {
        mesh_file >> global_node[k];
        k++;
      }
    }

    // Initialize counter
    unsigned l = 0;

    // Read the edge information
    for (unsigned i = 0; i < n_element; i++)
    {
      for (unsigned j = 0; j < n_local_node; j++)
      {
        mesh_file >> edgeinfo[l];
        l++;
      }
    }

    mesh_file.close();

    // Create a vector of boolean so as not to create the same node twice
    std::vector<bool> done(n_node);
    for (unsigned i = 0; i < n_node; i++)
    {
      done[i] = false;
    }

    // Resize the Node vector
    Node_pt.resize(n_node, 0);

    // Resize the Element vector
    Element_pt.resize(n_element);


    // Process curve file to extract information about boundaries
    // ----------------------------------------------------------

    // Important: the input file must NOT have NURBS curve
    std::ifstream curve_file(curve_file_name.c_str(), std::ios_base::in);

    // Number of curves
    unsigned n_curv;
    curve_file >> n_curv;

    // Storage of several information for each curve
    Vector<Vector<int>> curv;

    // Resize to n_curv rows
    curv.resize(n_curv);

    // Curve type
    unsigned type;

    // Loop over the curves to extract information
    for (unsigned i = 0; i < n_curv; i++)
    {
      curve_file >> type;
      if (type == 1)
      {
        curv[i].resize(4);
        curv[i][0] = type;
        for (unsigned j = 1; j < 4; j++)
        {
          curve_file >> curv[i][j];
        }
      }
      else if (type == 2)
      {
        curv[i].resize(5);
        curv[i][0] = type;
        for (unsigned j = 1; j < 5; j++)
        {
          curve_file >> curv[i][j];
        }
      }
      else
      {
        std::ostringstream error_stream;
        error_stream << "Current we can only process curves of\n"
                     << "type 1 (straight lines) and 2 (circular arcs\n"
                     << "You've specified: type " << type << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    curve_file.close();

    // Searching the number of boundaries
    int d = 0;
    for (unsigned i = 0; i < n_curv; i++)
    {
      if (d < curv[i][1])
      {
        d = curv[i][1]; // the boundary code is the 2nd value of each curve
      }
    }
    oomph_info << "The number of boundaries is " << d << std::endl;

    // Set number of boundaries
    if (d > 0)
    {
      set_nboundary(d);
    }

    // Search the boundary number of node located on a boundary
    // If after this, boundary_of_node[j][0] is -1 then node j
    // is not located on any boundary.
    // If boundary_of_node[j][0] is positive, the node is located
    // on the boundary indicated by that number.
    // If boundary_of_node[j][1] is also positive, the node is also
    // located on that boundary. Note: We're ignoring the (remote)
    // possibility that node is located on 3 or more boundaries
    // as one of them would have to be an internal boundary which
    // would be odd...
    Vector<Vector<int>> boundary_of_node;
    boundary_of_node.resize(n_node);
    unsigned n;
    for (unsigned i = 0; i < n_node; i++)
    {
      n = 0;
      boundary_of_node[i].resize(2);
      boundary_of_node[i][0] = -1;
      boundary_of_node[i][1] = -1;
      if (vertinfo[i] == 2) // it's an endpoint
      {
        for (unsigned j = 0; j < n_curv; j++)
        {
          for (unsigned m = 2; m < curv[j].size(); m++)
          {
            if (curv[j][m] ==
                static_cast<int>(i + 1)) // node number begins at 1
            { // in the mesh file !!!
              boundary_of_node[i][n] = curv[j][1];
              n++;
            }
          }
        }
      }
      if (vertinfo[i] > 20)
      {
        int a = 0;
        a = (vertinfo[i]) / 20;
        int b;
        b = icurv[a - 1]; // 1st value of vector at [0] !!
        boundary_of_node[i][0] =
          curv[b - 1][1]; // 1st value of vector at [0] !!
      }
    }


    // Create the elements
    //--------------------

    unsigned count = 0;
    unsigned c;
    for (unsigned e = 0; e < n_element; e++)
    {
      // Build simplex four node quad in the scaffold mesh
      Element_pt[e] = new QElement<2, 2>;

      // Construction of the two first nodes of the element
      for (unsigned j = 0; j < 2; j++)
      {
        c = global_node[count];
        if (done[c - 1] == false) // c-1 because node number begins
        // at 1 in the mesh file
        {
          // If the node is located on a boundary construct a boundary node
          if ((d > 0) && ((boundary_of_node[c - 1][0] > 0) ||
                          (boundary_of_node[c - 1][1] > 0)))
          {
            // Construct a boundary node
            Node_pt[c - 1] = finite_element_pt(e)->construct_boundary_node(j);
            // Add to the look=up schemes
            if (boundary_of_node[c - 1][0] > 0)
            {
              add_boundary_node(boundary_of_node[c - 1][0] - 1, Node_pt[c - 1]);
            }
            if (boundary_of_node[c - 1][1] > 0)
            {
              add_boundary_node(boundary_of_node[c - 1][1] - 1, Node_pt[c - 1]);
            }
          }
          // Otherwise construct a normal node
          else
          {
            Node_pt[c - 1] = finite_element_pt(e)->construct_node(j);
          }
          done[c - 1] = true;
          Node_pt[c - 1]->x(0) = x_node[c - 1];
          Node_pt[c - 1]->x(1) = y_node[c - 1];
        }
        else
        {
          finite_element_pt(e)->node_pt(j) = Node_pt[c - 1];
        }
        count++;
      }

      // Construction of the third node not in anticlockwise order
      c = global_node[count + 1];
      if (done[c - 1] ==
          false) // c-1 because node number begins at 1 in the mesh file
      {
        // If the node is on a boundary, construct a boundary node
        if ((d > 0) && ((boundary_of_node[c - 1][0] > 0) ||
                        (boundary_of_node[c - 1][1] > 0)))
        {
          // Construct the node
          Node_pt[c - 1] = finite_element_pt(e)->construct_boundary_node(2);
          // Add to the look-up schemes
          if (boundary_of_node[c - 1][0] > 0)
          {
            add_boundary_node(boundary_of_node[c - 1][0] - 1, Node_pt[c - 1]);
          }
          if (boundary_of_node[c - 1][1] > 0)
          {
            add_boundary_node(boundary_of_node[c - 1][1] - 1, Node_pt[c - 1]);
          }
        }
        // otherwise construct a normal node
        else
        {
          Node_pt[c - 1] = finite_element_pt(e)->construct_node(2);
        }
        done[c - 1] = true;
        Node_pt[c - 1]->x(0) = x_node[c - 1];
        Node_pt[c - 1]->x(1) = y_node[c - 1];
      }
      else
      {
        finite_element_pt(e)->node_pt(2) = Node_pt[c - 1];
      }

      count++;

      // Construction of the fourth node
      c = global_node[count - 1];
      if (done[c - 1] ==
          false) // c-1 because node number begins at 1 in the mesh file
      {
        // If the node is on a boundary, constuct a boundary node
        if ((d > 0) && ((boundary_of_node[c - 1][0] > 0) ||
                        (boundary_of_node[c - 1][1] > 0)))
        {
          // Construct the boundary node
          Node_pt[c - 1] = finite_element_pt(e)->construct_boundary_node(3);
          // Add to the look-up schemes
          if (boundary_of_node[c - 1][0] > 0)
          {
            add_boundary_node(boundary_of_node[c - 1][0] - 1, Node_pt[c - 1]);
          }
          if (boundary_of_node[c - 1][1] > 0)
          {
            add_boundary_node(boundary_of_node[c - 1][1] - 1, Node_pt[c - 1]);
          }
        }
        // Otherwise construct a normal node
        else
        {
          Node_pt[c - 1] = finite_element_pt(e)->construct_node(3);
        }
        done[c - 1] = true;
        Node_pt[c - 1]->x(0) = x_node[c - 1];
        Node_pt[c - 1]->x(1) = y_node[c - 1];
      }
      else
      {
        finite_element_pt(e)->node_pt(3) = Node_pt[c - 1];
      }

      count++;
    }
  }

} // namespace oomph
