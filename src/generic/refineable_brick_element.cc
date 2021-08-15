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
#include <algorithm>

#include "mesh.h"
#include "algebraic_elements.h"
#include "macro_element_node_update_element.h"
#include "refineable_brick_element.h"


namespace oomph
{
  //========================================================================
  ///  Print corner nodes, use colour (default "BLACK")
  ///  in right order so that tecplot can draw a cube without crossed lines
  //========================================================================
  void RefineableQElement<3>::output_corners(std::ostream& outfile,
                                             const std::string& colour) const
  {
    Vector<double> s(3);
    Vector<double> corner(3);

    outfile << "ZONE I=2, J=2, K=2 C=" << colour << std::endl;

    s[0] = -1.0;
    s[1] = -1.0;
    s[2] = -1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << corner[2] << " "
            << Number << std::endl;

    s[0] = -1.0;
    s[1] = -1.0;
    s[2] = 1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << corner[2] << " "
            << Number << std::endl;

    s[0] = -1.0;
    s[1] = 1.0;
    s[2] = -1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << corner[2] << " "
            << Number << std::endl;

    s[0] = -1.0;
    s[1] = 1.0;
    s[2] = 1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << corner[2] << " "
            << Number << std::endl;

    // next face


    s[0] = 1.0;
    s[1] = -1.0;
    s[2] = -1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << corner[2] << " "
            << Number << std::endl;

    s[0] = 1.0;
    s[1] = -1.0;
    s[2] = 1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << corner[2] << " "
            << Number << std::endl;

    s[0] = 1.0;
    s[1] = 1.0;
    s[2] = -1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << corner[2] << " "
            << Number << std::endl;

    s[0] = 1.0;
    s[1] = 1.0;
    s[2] = 1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << corner[2] << " "
            << Number << std::endl;


    //  outfile << "TEXT  CS = GRID, X = " << corner[0] <<
    //   ",Y = " << corner[1] << ",Z = " << corner[2] <<
    //   ", HU = GRID, H = 0.01, AN = MIDCENTER, T =\""
    //          << Number << "\"" << std::endl;
  }


  //==================================================================
  /// Setup static matrix for coincidence between son nodal points and
  /// father boundaries:
  ///
  /// Father_bound[nnode_1d](nnode_son,son_type)={RU/RF/RD/RB/.../ OMEGA}
  ///
  /// so that node nnode_son in element of type son_type lies
  /// on boundary/vertex Father_bound[nnode_1d](nnode_son,son_type) in its
  /// father element. If the node doesn't lie on a boundary
  /// the value is OMEGA.
  //==================================================================
  void RefineableQElement<3>::setup_father_bounds()
  {
    using namespace OcTreeNames;

    // Find the number of nodes along a 1D edge
    unsigned n_p = nnode_1d();

    // Allocate space for the boundary information
    Father_bound[n_p].resize(n_p * n_p * n_p, 8);

    // Initialise: By default points are not on the boundary
    for (unsigned n = 0; n < n_p * n_p * n_p; n++)
    {
      for (unsigned ison = 0; ison < 8; ison++)
      {
        Father_bound[n_p](n, ison) = Tree::OMEGA;
      }
    }

    for (int i_son = LDB; i_son <= RUF; i_son++)
    {
      // vector representing the son
      Vector<int> vect_son(3);
      // vector representing (at the end) the boundaries
      Vector<int> vect_bound(3);
      vect_son = OcTree::Direction_to_vector[i_son];
      for (unsigned i0 = 0; i0 < n_p; i0++)
      {
        for (unsigned i1 = 0; i1 < n_p; i1++)
        {
          for (unsigned i2 = 0; i2 < n_p; i2++)
          {
            // Initialisation to make it work
            for (unsigned i = 0; i < 3; i++)
            {
              vect_bound[i] = -vect_son[i];
            }
            // Seting up the boundaries coordinates as if the coordinates
            // of vect_bound had been initialised to 0 if it were so,
            // vect_bound would be the vector of the boundaries in the son
            // itself.

            if (i0 == 0)
            {
              vect_bound[0] = -1;
            }
            if (i0 == n_p - 1)
            {
              vect_bound[0] = 1;
            }
            if (i1 == 0)
            {
              vect_bound[1] = -1;
            }
            if (i1 == n_p - 1)
            {
              vect_bound[1] = 1;
            }
            if (i2 == 0)
            {
              vect_bound[2] = -1;
            }
            if (i2 == n_p - 1)
            {
              vect_bound[2] = 1;
            }

            // The effect of this is to filter the boundaries to keep only the
            // ones which are effectively father boundaries.
            // -- if the node is not on a "i0 boundary", we still
            //    have vect_bound[0]=-vect_son[0]
            //    and the result is vect_bound[0]=0
            //    (he is not on this boundary for his father)
            // -- if he is on a son's boundary  which is not one of
            //    the father -> same thing
            // -- if he is on a boundary which is one of his father's,
            //    vect_bound[i]=vect_son[i]
            //    and the new vect_bound[i] is the same as the old one
            for (int i = 0; i < 3; i++)
            {
              vect_bound[i] = (vect_bound[i] + vect_son[i]) / 2;
            }

            // Return the result as {U,R,D,...RDB,LUF,OMEGA}
            Father_bound[n_p](i0 + n_p * i1 + n_p * n_p * i2, i_son) =
              OcTree::Vector_to_direction[vect_bound];

          } // Loop over i2
        } // Loop over i1
      } // Loop over i0
    } // Loop over i_son

  } // setup_father_bounds()


  //==================================================================
  /// Determine Vector of boundary conditions along the element's boundary
  /// bound.
  ///
  /// This function assumes that the same boundary condition is applied
  /// along the entire area of an element's face (of course, the
  /// vertices combine the boundary conditions of their two adjacent edges
  /// in the most restrictive combination. Hence, if we're at a vertex,
  /// we apply the most restrictive boundary condition of the
  /// two adjacent edges. If we're on an edge (in its proper interior),
  /// we apply the least restrictive boundary condition of all nodes
  /// along the edge.
  ///
  /// Usual convention:
  ///   - bound_cons[ival]=0 if value ival on this boundary is free
  ///   - bound_cons[ival]=1 if value ival on this boundary is pinned
  //==================================================================
  void RefineableQElement<3>::get_bcs(int bound, Vector<int>& bound_cons) const
  {
    using namespace OcTreeNames;

    // Max. number of nodal data values in element
    unsigned nvalue = ncont_interpolated_values();
    // Set up temporary vectors to hold edge boundary conditions
    Vector<int> bound_cons1(nvalue), bound_cons2(nvalue);
    Vector<int> bound_cons3(nvalue);

    Vector<int> vect1(3), vect2(3), vect3(3);
    Vector<int> vect_elem;
    Vector<int> notzero;
    int n = 0;

    vect_elem = OcTree::Direction_to_vector[bound];

    // Just to see if bound is a face, an edge, or a vertex, n stores
    // the number of non-zero values in the vector reprensenting the bound
    // and the vector notzero stores the position of these values
    for (int i = 0; i < 3; i++)
    {
      if (vect_elem[i] != 0)
      {
        n++;
        notzero.push_back(i);
      }
    }

    switch (n)
    {
        // If there is only one non-zero value, bound is a face
      case 1:
        get_face_bcs(bound, bound_cons);
        break;

        // If there are two non-zero values, bound is an edge
      case 2:

        for (unsigned i = 0; i < 3; i++)
        {
          vect1[i] = 0;
          vect2[i] = 0;
        }
        // vect1 and vect2 are the vector of the two faces adjacent to bound
        vect1[notzero[0]] = vect_elem[notzero[0]];
        vect2[notzero[1]] = vect_elem[notzero[1]];

        get_face_bcs(OcTree::Vector_to_direction[vect1], bound_cons1);
        get_face_bcs(OcTree::Vector_to_direction[vect2], bound_cons2);
        // get the most restrictive bc
        for (unsigned k = 0; k < nvalue; k++)
        {
          bound_cons[k] = (bound_cons1[k] || bound_cons2[k]);
        }
        break;

        // If there are three non-zero value, bound is a vertex
      case 3:

        for (unsigned i = 0; i < 3; i++)
        {
          vect1[i] = 0;
          vect2[i] = 0;
          vect3[i] = 0;
        }
        // vectors to the three adjacent faces of the vertex
        vect1[0] = vect_elem[0];
        vect2[1] = vect_elem[1];
        vect3[2] = vect_elem[2];

        get_face_bcs(OcTree::Vector_to_direction[vect1], bound_cons1);
        get_face_bcs(OcTree::Vector_to_direction[vect2], bound_cons2);
        get_face_bcs(OcTree::Vector_to_direction[vect3], bound_cons3);


        // set the bcs to the most restrictive ones
        for (unsigned k = 0; k < nvalue; k++)
        {
          bound_cons[k] = (bound_cons1[k] || bound_cons2[k] || bound_cons3[k]);
        }
        break;

      default:
        throw OomphLibError("Make sure you are not giving OMEGA as bound",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }

  //==================================================================
  /// Determine Vector of boundary conditions along the element's
  /// face (R/L/U/D/B/F) -- BC is the least restrictive combination
  /// of all the nodes on this face
  ///
  /// Usual convention:
  ///   - bound_cons[ival]=0 if value ival on this boundary is free
  ///   - bound_cons[ival]=1 if value ival on this boundary is pinned
  //==================================================================
  void RefineableQElement<3>::get_face_bcs(const int& face,
                                           Vector<int>& bound_cons) const
  {
    using namespace OcTreeNames;

    // Number of nodes along 1D edge
    unsigned n_p = nnode_1d();
    // the four corner nodes on the boundary
    unsigned node1, node2, node3, node4;

    // Set the four corner nodes for the face
    switch (face)
    {
      case U:
        node1 = n_p * n_p * n_p - 1;
        node2 = n_p * n_p - 1;
        node3 = n_p * (n_p - 1);
        node4 = n_p * (n_p * n_p - 1);
        break;

      case D:
        node1 = 0;
        node2 = n_p - 1;
        node3 = (n_p * n_p + 1) * (n_p - 1);
        node4 = n_p * n_p * (n_p - 1);
        break;

      case R:
        node1 = n_p - 1;
        node2 = (n_p * n_p + 1) * (n_p - 1);
        node3 = n_p * n_p * n_p - 1;
        node4 = n_p * n_p - 1;
        break;

      case L:
        node1 = 0;
        node2 = n_p * (n_p - 1);
        node3 = n_p * (n_p * n_p - 1);
        node4 = n_p * n_p * (n_p - 1);
        break;

      case B:
        node1 = 0;
        node2 = n_p - 1;
        node3 = n_p * n_p - 1;
        node4 = n_p * (n_p - 1);
        break;

      case F:
        node1 = n_p * n_p * n_p - 1;
        node2 = n_p * (n_p * n_p - 1);
        node3 = n_p * n_p * (n_p - 1);
        node4 = (n_p - 1) * (n_p * n_p + 1);
        break;

      default:
        std::ostringstream error_stream;
        error_stream << "Wrong edge " << face << " passed\n";

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Max. number of nodal data values in element
    unsigned maxnvalue = ncont_interpolated_values();

    // Loop over all values, the least restrictive value is
    // the multiplication of the boundary conditions at the 4 nodes
    // Assuming that free is always zero and pinned is one
    for (unsigned k = 0; k < maxnvalue; k++)
    {
      bound_cons[k] =
        node_pt(node1)->is_pinned(k) * node_pt(node2)->is_pinned(k) *
        node_pt(node3)->is_pinned(k) * node_pt(node4)->is_pinned(k);
    }
  }


  //==================================================================
  /// Given an element edge/vertex, return a Vector which contains
  /// all the (mesh-)boundary numbers that this element edge/vertex
  /// lives on.
  ///
  /// For proper edges, the boundary is the one (if any) that is shared by
  /// both vertex nodes). For vertex nodes, we just return their
  /// boundaries.
  //==================================================================
  void RefineableQElement<3>::get_boundaries(const int& element,
                                             std::set<unsigned>& boundary) const
  {
    using namespace OcTreeNames;

    // Number of 1d nodes along an edge
    unsigned n_p = nnode_1d();
    // Left and right-hand nodes
    int node[4];
    int a = 0, b = 0;
    int n = 0;
    Vector<int> vect_face(3);
    vect_face = OcTree::Direction_to_vector[element];
    // Set the left (lower) and right (upper) nodes for the edge

    // this is to know what is the type of element (face, edge, vertex)
    // we just need to count the number of values equal to 0 in the
    // vector representing this element

    // Local storage for the node-numbers in given directions
    // Initialise to zero (assume at LH end of each domain)
    int one_d_node_number[3] = {0, 0, 0};
    // n stores the number of values equal to 0, a is the position of the
    // last 0-value ;b is the position of the last non0-value
    for (int i = 0; i < 3; i++)
    {
      if (vect_face[i] == 0)
      {
        a = i;
        n++;
      }
      else
      {
        b = i;
        // If we are at the extreme (RH) end of the face,
        // set the node number accordingly
        if (vect_face[i] == 1)
        {
          one_d_node_number[i] = n_p - 1;
        }
      }
    }

    switch (n)
    {
        // if n=0 element is a vertex, and need to look at only one node
      case 0:
        node[0] = one_d_node_number[0] + n_p * one_d_node_number[1] +
                  n_p * n_p * one_d_node_number[2];
        node[1] = node[0];
        node[2] = node[0];
        node[3] = node[0];
        break;

        // if n=1 element is an edge, and we need to look at two nodes
      case 1:
        if (a == 0)
        {
          node[0] = (n_p - 1) + n_p * one_d_node_number[1] +
                    n_p * n_p * one_d_node_number[2];
          node[1] =
            n_p * one_d_node_number[1] + n_p * n_p * one_d_node_number[2];
        }
        else if (a == 1)
        {
          node[0] = n_p * (n_p - 1) + one_d_node_number[0] +
                    n_p * n_p * one_d_node_number[2];
          node[1] = one_d_node_number[0] + n_p * n_p * one_d_node_number[2];
        }
        else if (a == 2)
        {
          node[0] = one_d_node_number[0] + n_p * one_d_node_number[1] +
                    n_p * n_p * (n_p - 1);
          node[1] = one_d_node_number[0] + n_p * one_d_node_number[1];
        }
        node[2] = node[1];
        node[3] = node[1];
        break;

        // if n=2 element is a face, and we need to look at its 4 nodes
      case 2:
        if (b == 0)
        {
          node[0] =
            one_d_node_number[0] + n_p * n_p * (n_p - 1) + n_p * (n_p - 1);
          node[1] = one_d_node_number[0] + n_p * (n_p - 1);
          node[2] = one_d_node_number[0] + n_p * n_p * (n_p - 1);
          node[3] = one_d_node_number[0];
        }
        else if (b == 1)
        {
          node[0] =
            n_p * one_d_node_number[1] + n_p * n_p * (n_p - 1) + (n_p - 1);
          node[1] = n_p * one_d_node_number[1] + (n_p - 1);
          node[2] = n_p * one_d_node_number[1] + n_p * n_p * (n_p - 1);
          node[3] = n_p * one_d_node_number[1];
        }
        else if (b == 2)
        {
          node[0] =
            n_p * n_p * one_d_node_number[2] + n_p * (n_p - 1) + (n_p - 1);
          node[1] = n_p * n_p * one_d_node_number[2] + (n_p - 1);
          node[2] = n_p * n_p * one_d_node_number[2] + n_p * (n_p - 1);
          node[3] = n_p * n_p * one_d_node_number[2];
        }
        break;
      default:
        throw OomphLibError("Make sure you are not giving OMEGA as boundary",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }


    // Empty boundary set: Edge does not live on any boundary
    boundary.clear();

    // Storage for the boundaries at the four nodes
    Vector<std::set<unsigned>*> node_boundaries_pt(4, 0);

    // Loop over the four nodes and get the boundary information
    for (unsigned i = 0; i < 4; i++)
    {
      node_pt(node[i])->get_boundaries_pt(node_boundaries_pt[i]);
    }


    // Now work out the intersections
    Vector<std::set<unsigned>> boundary_aux(2);

    for (unsigned i = 0; i < 2; i++)
    {
      // If the two nodes both lie on boundaries
      if ((node_boundaries_pt[2 * i] != 0) &&
          (node_boundaries_pt[2 * i + 1] != 0))
      {
        // Find the intersection (the common boundaries) of these nodes
        std::set_intersection(
          node_boundaries_pt[2 * i]->begin(),
          node_boundaries_pt[2 * i]->end(),
          node_boundaries_pt[2 * i + 1]->begin(),
          node_boundaries_pt[2 * i + 1]->end(),
          inserter(boundary_aux[i], boundary_aux[i].begin()));
      }
    }

    // Now calculate the total intersection
    set_intersection(boundary_aux[0].begin(),
                     boundary_aux[0].end(),
                     boundary_aux[1].begin(),
                     boundary_aux[1].end(),
                     inserter(boundary, boundary.begin()));
  }


  //===================================================================
  /// Return the value of the intrinsic boundary coordinate interpolated
  /// along the face
  //===================================================================
  void RefineableQElement<3>::interpolated_zeta_on_face(
    const unsigned& boundary,
    const int& face,
    const Vector<double>& s,
    Vector<double>& zeta)
  {
    using namespace OcTreeNames;

    // Number of nodes along an edge
    unsigned nnodes_1d = nnode_1d();

    // Number of nodes on a face
    unsigned nnodes_2d = nnodes_1d * nnodes_1d;

    // Total number of nodes
    unsigned nnodes_3d = nnode();

    // Storage for the shape functions
    Shape psi(nnodes_3d);

    // Get the shape functions at the passed position
    this->shape(s, psi);

    // Unsigned data that give starts and increments for the loop
    // over the nodes on the faces.
    unsigned start = 0, increment1 = 1, increment2 = 1;

    // Flag to record if actually on a face or an edge
    bool on_edge = true;

    // Which face?
    switch (face)
    {
      case L:
#ifdef PARANOID
        if (s[0] != -1.0)
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Left face\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Start is zero (bottom-left-back corner)
        increment1 = nnodes_1d;
        increment2 = 0;
        on_edge = false;
        break;

      case R:
#ifdef PARANOID
        if (s[0] != 1.0)
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Right face\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Start is bottom-right-back corner
        start = nnodes_1d - 1;
        increment1 = nnodes_1d;
        increment2 = 0;
        on_edge = false;
        break;

      case D:
#ifdef PARANOID
        if (s[1] != -1.0)
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Bottom face\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Start is zero and increments2 is nnode_2d-nnode_1d
        increment2 = nnodes_2d - nnodes_1d;
        on_edge = false;
        break;

      case U:
#ifdef PARANOID
        if (s[1] != 1.0)
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Upper face\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Start is top-left-back corner and increments2 is nnode_2d-nnode_1d
        start = nnodes_2d - nnodes_1d;
        increment2 = start;
        on_edge = false;
        break;

      case B:
#ifdef PARANOID
        if (s[2] != -1.0)
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Back face\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Start is zero and increments are 1 and 0
        increment2 = 0;
        on_edge = false;
        break;

      case F:
#ifdef PARANOID
        if (s[2] != 1.0)
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Front face\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Start is bottom-left-front corner
        start = nnodes_3d - nnodes_2d;
        increment2 = 0;
        on_edge = false;
        break;

      case LF:
#ifdef PARANOID
        if ((s[0] != -1.0) || (s[2] != 1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Front-Left edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        start = nnodes_3d - nnodes_2d;
        increment1 = nnodes_1d;
        break;

      case LD:
#ifdef PARANOID
        if ((s[0] != -1.0) || (s[1] != -1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Bottom-Left edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        increment1 = nnodes_2d;
        break;

      case LU:
#ifdef PARANOID
        if ((s[0] != -1.0) || (s[1] != 1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Upper-Left edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        start = nnodes_2d - nnodes_1d;
        increment1 = nnodes_2d;
        break;

      case LB:
#ifdef PARANOID
        if ((s[0] != -1.0) || (s[2] != -1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Back-Left edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        increment1 = nnodes_1d;
        break;

      case RF:
#ifdef PARANOID
        if ((s[0] != 1.0) || (s[2] != 1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Front-Right edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        start = nnodes_3d - 1;
        increment1 = -nnodes_1d;
        break;

      case RD:
#ifdef PARANOID
        if ((s[0] != 1.0) || (s[1] != -1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Bottom-Right edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        start = nnodes_1d - 1;
        increment1 = nnodes_2d;
        break;

      case RU:
#ifdef PARANOID
        if ((s[0] != 1.0) || (s[1] != 1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Upper-Right edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        start = nnodes_2d - 1;
        increment1 = nnodes_2d;
        break;

      case RB:
#ifdef PARANOID
        if ((s[0] != 1.0) || (s[2] != -1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Back-Right edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        start = nnodes_1d - 1;
        increment1 = nnodes_1d;
        break;

      case DB:
#ifdef PARANOID
        if ((s[1] != -1.0) || (s[2] != -1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Back-Bottom edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        break;

      case DF:
#ifdef PARANOID
        if ((s[1] != -1.0) || (s[2] != 1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Front-Bottom edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        start = nnodes_3d - nnodes_2d;
        break;

      case UB:
#ifdef PARANOID
        if ((s[1] != 1.0) || (s[2] != -1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Back-Upper edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        start = nnodes_2d - nnodes_1d;
        break;

      case UF:
#ifdef PARANOID
        if ((s[1] != 1.0) || (s[2] != 1.0))
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1] << " " << s[2]
                       << " is not on Upper-Front edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        start = nnodes_3d - nnodes_1d;
        break;

      default:

        std::ostringstream error_stream;
        error_stream << "Wrong face " << OcTree::Direct_string[face]
                     << " passed" << std::endl;
        error_stream << "Trouble at : s= [" << s[0] << " " << s[1] << " "
                     << s[2] << "]\n";
        Vector<double> x(3);
        this->interpolated_x(s, x);
        error_stream << "corresponding to : x= [" << x[0] << " " << x[1] << " "
                     << x[2] << "]\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Initialise the intrinsic coordinate
    zeta[0] = 0.0;
    zeta[1] = 0.0;

    // Set the start node number
    unsigned node = start;

    if (on_edge)
    {
      for (unsigned l1 = 0; l1 < nnodes_1d; l1++)
      {
        // Get the intrinsic coordinate
        Vector<double> node_zeta(2);
        node_pt(node)->get_coordinates_on_boundary(boundary, node_zeta);

        // Now multiply by the shape function
        zeta[0] += node_zeta[0] * psi(node);
        zeta[1] += node_zeta[1] * psi(node);

        // Update node
        node += increment1;
      }
    }
    else
    {
      for (unsigned l2 = 0; l2 < nnodes_1d; l2++)
      {
        for (unsigned l1 = 0; l1 < nnodes_1d; l1++)
        {
          // Get the intrinsic coordinate
          Vector<double> node_zeta(2);
          node_pt(node)->get_coordinates_on_boundary(boundary, node_zeta);

          // Now multiply by the shape function
          zeta[0] += node_zeta[0] * psi(node);
          zeta[1] += node_zeta[1] * psi(node);

          // Update node
          node += increment1;
        }
        // Update node
        node += increment2;
      }
    }
  }

  //===================================================================
  /// If a neighbouring element has already created a node at a
  /// position corresponding to the local fractional position within the
  /// present element, s_fraction, return
  /// a pointer to that node. If not, return NULL (0).
  //===================================================================
  Node* RefineableQElement<3>::node_created_by_neighbour(
    const Vector<double>& s_fraction, bool& is_periodic)
  {
    using namespace OcTreeNames;

    // Calculate the faces/edges on which the node lies
    Vector<int> faces;
    Vector<int> edges;

    if (s_fraction[0] == 0.0)
    {
      faces.push_back(L);
      if (s_fraction[1] == 0.0)
      {
        edges.push_back(LD);
      }
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(LB);
      }
      if (s_fraction[1] == 1.0)
      {
        edges.push_back(LU);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(LF);
      }
    }

    if (s_fraction[0] == 1.0)
    {
      faces.push_back(R);
      if (s_fraction[1] == 0.0)
      {
        edges.push_back(RD);
      }
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(RB);
      }
      if (s_fraction[1] == 1.0)
      {
        edges.push_back(RU);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(RF);
      }
    }

    if (s_fraction[1] == 0.0)
    {
      faces.push_back(D);
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(DB);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(DF);
      }
    }

    if (s_fraction[1] == 1.0)
    {
      faces.push_back(U);
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(UB);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(UF);
      }
    }

    if (s_fraction[2] == 0.0)
    {
      faces.push_back(B);
    }

    if (s_fraction[2] == 1.0)
    {
      faces.push_back(F);
    }

    // Find the number of faces
    unsigned n_face = faces.size();

    // Find the number of edges
    unsigned n_edge = edges.size();

    Vector<unsigned> translate_s(3);
    Vector<double> s_lo_neigh(3);
    Vector<double> s_hi_neigh(3);
    Vector<double> s(3);

    int neigh_face, diff_level;

    // Loop over the faces on which the node lies
    //-------------------------------------------
    for (unsigned j = 0; j < n_face; j++)
    {
      // Boolean to indicate whether or not the neighbour has a
      // different Tree root
      bool in_neighbouring_tree;

      // Find pointer to neighbouring element along face
      OcTree* neigh_pt;
      neigh_pt = octree_pt()->gteq_face_neighbour(faces[j],
                                                  translate_s,
                                                  s_lo_neigh,
                                                  s_hi_neigh,
                                                  neigh_face,
                                                  diff_level,
                                                  in_neighbouring_tree);

      // Neighbour exists
      if (neigh_pt != 0)
      {
        // Have its nodes been created yet?
        if (neigh_pt->object_pt()->nodes_built())
        {
          // We now need to translate the nodal location, defined in terms
          // of the fractional coordinates of the present element into
          // those of its neighbour. For this we use the information returned
          // to use from the octree function.

          // Calculate the local coordinate in the neighbour
          // Note that we need to use the translation scheme in case
          // the local coordinates are swapped in the neighbour.
          for (unsigned i = 0; i < 3; i++)
          {
            s[i] = s_lo_neigh[i] +
                   s_fraction[translate_s[i]] * (s_hi_neigh[i] - s_lo_neigh[i]);
          }

          // Find the node in the neighbour
          Node* neighbour_node_pt =
            neigh_pt->object_pt()->get_node_at_local_coordinate(s);

          // If there is a node, return it
          if (neighbour_node_pt != 0)
          {
            // Now work out whether it's a periodic boundary. This is
            // only possible if we have moved into a neighbouring tree
            if (in_neighbouring_tree)
            {
              // Return whether the neighbour is periodic
              is_periodic =
                octree_pt()->root_pt()->is_neighbour_periodic(faces[j]);
            }

            // Return the neighbour node pointer
            return neighbour_node_pt;
          }
        } // if (neigh_pt->object_pt()->nodes_built())
      } // if (neigh_pt!=0)
    } // for (unsigned j=0;j<n_face;j++)

    // Loop over the edges on which the node lies
    //------------------------------------------
    for (unsigned j = 0; j < n_edge; j++)
    {
      // Even if we restrict ourselves to true edge neighbours (i.e.
      // elements that are not also face neighbours) there may be multiple
      // edge neighbours across the edges between multiple root octrees.
      // When making the first call to OcTree::gteq_true_edge_neighbour(...)
      // we simply return the first one of these multiple edge neighbours
      // (if there are any at all, of course) and also return the total number
      // of true edge neighbours. If the node in question already exists
      // on the first edge neighbour we're done. If it doesn't it may exist
      // on other edge neighbours so we repeat the process over all
      // other edge neighbours (bailing out if a node is found, of course).

      // Initially return the zero-th true edge neighbour
      unsigned i_root_edge_neighbour = 0;

      // Initialise the total number of true edge neighbours
      unsigned nroot_edge_neighbour = 0;

      // Keep searching until we've found the node or until we've checked
      // all available edge neighbours
      bool keep_searching = true;
      while (keep_searching)
      {
        // Find pointer to neighbouring element along edge
        OcTree* neigh_pt;
        neigh_pt = octree_pt()->gteq_true_edge_neighbour(edges[j],
                                                         i_root_edge_neighbour,
                                                         nroot_edge_neighbour,
                                                         translate_s,
                                                         s_lo_neigh,
                                                         s_hi_neigh,
                                                         neigh_face,
                                                         diff_level);

        // Neighbour exists
        if (neigh_pt != 0)
        {
          // Have its nodes been created yet?
          if (neigh_pt->object_pt()->nodes_built())
          {
            // We now need to translate the nodal location, defined in terms
            // of the fractional coordinates of the present element into
            // those of its neighbour. For this we use the information returned
            // to use from the octree function.

            // Calculate the local coordinate in the neighbour
            // Note that we need to use the translation scheme in case
            // the local coordinates are swapped in the neighbour.
            for (unsigned i = 0; i < 3; i++)
            {
              s[i] = s_lo_neigh[i] + s_fraction[translate_s[i]] *
                                       (s_hi_neigh[i] - s_lo_neigh[i]);
            }

            // Find the node in the neighbour
            Node* neighbour_node_pt =
              neigh_pt->object_pt()->get_node_at_local_coordinate(s);

            // If there is a node, return it
            if (neighbour_node_pt != 0)
            {
              // Get the faces on which the edge lies
              Vector<int> faces_attached_to_edge =
                OcTree::faces_of_common_edge(edges[j]);

              // Get the number of entries in the vector
              unsigned n_faces_attached_to_edge = faces_attached_to_edge.size();

              // Loop over the faces
              for (unsigned i_face = 0; i_face < n_faces_attached_to_edge;
                   i_face++)
              {
                // Is the node periodic in the face direction?
                is_periodic = octree_pt()->root_pt()->is_neighbour_periodic(
                  faces_attached_to_edge[i_face]);

                // Check if the edge is periodic in the i_face-th face direction
                if (is_periodic)
                {
                  // We're done!
                  break;
                }
              } // for (unsigned
                // i_face=0;i_face<n_faces_attached_to_edge;i_face++)

              // Return the neighbour node pointer
              return neighbour_node_pt;
            } // if (neighbour_node_pt!=0)
          } // if (neigh_pt->object_pt()->nodes_built())
        } // if (neigh_pt!=0)

        // Keep searching, but only if there are further edge neighbours
        // Try next root edge neighbour
        i_root_edge_neighbour++;

        // Have we exhausted the search?
        if (i_root_edge_neighbour >= nroot_edge_neighbour)
        {
          // Stop searching
          keep_searching = false;
        }
      } // End of while keep searching over all true edge neighbours
    } // End of loop over edges

    // Node not found, return null
    return 0;
  }


  //==================================================================
  /// Build the element by doing the following:
  /// - Give it nodal positions (by establishing the pointers to its
  ///   nodes)
  /// - In the process create new nodes where required (i.e. if they
  ///   don't exist in father element or have already been created
  ///   while building new neighbour elements). Node building
  ///   involves the following steps:
  ///   - Get nodal position from father element.
  ///   - Establish the time-history of the newly created nodal point
  ///     (its coordinates and the previous values) consistent with
  ///     the father's history.
  ///   - Determine the boundary conditions of the nodes (newly
  ///     created nodes can only lie on the interior of any
  ///     edges of the father element -- this makes it possible to
  ///     to figure out what their bc should be...)
  ///   - Add node to the mesh's stoarge scheme for the boundary nodes.
  ///   - Add the new node to the mesh itself
  ///   - Doc newly created nodes in file "new_nodes.dat" in the directory
  ////    of the DocInfo object (only if it's open!)
  /// - Finally, excute the element-specific further_build()
  ///   (empty by default -- must be overloaded for specific elements).
  ///   This deals with any build operations that are not included
  ///   in the generic process outlined above. For instance, in
  ///   Crouzeix Raviart elements we need to initialise the internal
  ///   pressure values in manner consistent with the pressure
  ///   distribution in the father element.
  //==================================================================
  void RefineableQElement<3>::build(Mesh*& mesh_pt,
                                    Vector<Node*>& new_node_pt,
                                    bool& was_already_built,
                                    std::ofstream& new_nodes_file)
  {
    using namespace OcTreeNames;

    // Number of dimensions
    unsigned n_dim = 3;

    // Get the number of 1d nodes
    unsigned n_p = nnode_1d();

    // Check whether static father_bound needs to be created
    if (Father_bound[n_p].nrow() == 0)
    {
      setup_father_bounds();
    }

    // Pointer to my father (in octree impersonation)
    OcTree* father_pt = dynamic_cast<OcTree*>(octree_pt()->father_pt());

    // What type of son am I? Ask my octree representation...
    int son_type = octree_pt()->son_type();

    // Has somebody build me already? (If any nodes have not been built)
    if (!nodes_built())
    {
#ifdef PARANOID
      if (father_pt == 0)
      {
        std::string error_message =
          "Something fishy here: I have no father and yet \n";
        error_message += "I have no nodes. Who has created me then?!\n";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Indicate status:
      was_already_built = false;

      // Return pointer to father element
      RefineableQElement<3>* father_el_pt =
        dynamic_cast<RefineableQElement<3>*>(father_pt->object_pt());

      // Timestepper should be the same for all nodes in father
      // element -- use it create timesteppers for new nodes
      TimeStepper* time_stepper_pt =
        father_el_pt->node_pt(0)->time_stepper_pt();

      // Number of history values (incl. present)
      unsigned ntstorage = time_stepper_pt->ntstorage();

      // Currently we can't handle the case of generalised coordinates
      // since we haven't established how they should be interpolated
      // Buffer this case:
      if (father_el_pt->node_pt(0)->nposition_type() != 1)
      {
        throw OomphLibError("Can't handle generalised nodal positions (yet).",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      Vector<int> s_lo(n_dim);
      Vector<int> s_hi(n_dim);
      Vector<double> s(n_dim);
      Vector<double> x(n_dim);

      // Setup vertex coordinates in father element:
      //--------------------------------------------
      // find the s_lo coordinates
      s_lo = octree_pt()->Direction_to_vector[son_type];

      // Just scale them, because the Direction_to_vector
      // doesn't really gives s_lo;
      for (unsigned i = 0; i < n_dim; i++)
      {
        s_lo[i] = (s_lo[i] + 1) / 2 - 1;
      }

      // setup s_hi (Actually s_hi[i]=s_lo[i]+1)
      for (unsigned i = 0; i < n_dim; i++)
      {
        s_hi[i] = s_lo[i] + 1;
      }

      // Pass macro element pointer on to sons and
      // set coordinates in macro element
      if (father_el_pt->macro_elem_pt() != 0)
      {
        set_macro_elem_pt(father_el_pt->macro_elem_pt());
        for (unsigned i = 0; i < n_dim; i++)
        {
          s_macro_ll(i) =
            father_el_pt->s_macro_ll(i) +
            0.5 * (s_lo[i] + 1.0) *
              (father_el_pt->s_macro_ur(i) - father_el_pt->s_macro_ll(i));
          s_macro_ur(i) =
            father_el_pt->s_macro_ll(i) +
            0.5 * (s_hi[i] + 1.0) *
              (father_el_pt->s_macro_ur(i) - father_el_pt->s_macro_ll(i));
        }
      }


      // If the father element hasn't been generated yet, we're stuck...
      if (father_el_pt->node_pt(0) == 0)
      {
        throw OomphLibError(
          "Trouble: father_el_pt->node_pt(0)==0\n Can't build son element!\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        unsigned jnod = 0;
        Vector<double> s_fraction(n_dim);

        // Loop over nodes in element
        for (unsigned i0 = 0; i0 < n_p; i0++)
        {
          // Get the fractional position of the node in the direction of s[0]
          s_fraction[0] = local_one_d_fraction_of_node(i0, 0);
          // Local coordinate in father element
          s[0] = s_lo[0] + (s_hi[0] - s_lo[0]) * s_fraction[0];

          for (unsigned i1 = 0; i1 < n_p; i1++)
          {
            // Get the fractional position of the node in the direction of s[1]
            s_fraction[1] = local_one_d_fraction_of_node(i1, 1);
            // Local coordinate in father element
            s[1] = s_lo[1] + (s_hi[1] - s_lo[1]) * s_fraction[1];

            for (unsigned i2 = 0; i2 < n_p; i2++)
            {
              // Get the fractional position of the node in the direction of
              // s[2]
              s_fraction[2] = local_one_d_fraction_of_node(i2, 2);

              // Local coordinate in father element
              s[2] = s_lo[2] + (s_hi[2] - s_lo[2]) * s_fraction[2];

              // Local node number
              jnod = i0 + n_p * i1 + n_p * n_p * i2;

              // Initialise flag: So far, this node hasn't been built
              // or copied yet
              bool node_done = false;

              // Get the pointer to the node in the father; returns NULL
              // if there is not a node
              Node* created_node_pt =
                father_el_pt->get_node_at_local_coordinate(s);

              // Does this node already exist in father element?
              //------------------------------------------------
              bool node_exists_in_father = false;
              if (created_node_pt != 0)
              {
                // Remember this!
                node_exists_in_father = true;

                // Copy node across
                node_pt(jnod) = created_node_pt;

                // Make sure that we update the values at the node so that
                // they are consistent with the present representation.
                // This is only need for mixed interpolation where the value
                // at the father could now become active.

                // Loop over all history values
                for (unsigned t = 0; t < ntstorage; t++)
                {
                  // Get values from father element
                  // Note: get_interpolated_values() sets Vector size itself.
                  Vector<double> prev_values;
                  father_el_pt->get_interpolated_values(t, s, prev_values);
                  // Find the minimum number of values
                  //(either those stored at the node, or those returned by
                  // the function)
                  unsigned n_val_at_node = created_node_pt->nvalue();
                  unsigned n_val_from_function = prev_values.size();
                  // Use the ternary conditional operator here
                  unsigned n_var = n_val_at_node < n_val_from_function ?
                                     n_val_at_node :
                                     n_val_from_function;
                  // Assign the values that we can
                  for (unsigned k = 0; k < n_var; k++)
                  {
                    created_node_pt->set_value(t, k, prev_values[k]);
                  }
                }

                // Node has been created by copy
                node_done = true;
              }
              // Node does not exist in father element but might already
              //--------------------------------------------------------
              // have been created by neighbouring elements
              //-------------------------------------------
              else
              {
                // Boolean to check if the node is periodic
                bool is_periodic = false;

                // Was the node created by one of its neighbours
                // Whether or not the node lies on an edge can be determined
                // from the fractional position
                created_node_pt =
                  node_created_by_neighbour(s_fraction, is_periodic);

                // If so, then copy the pointer across
                if (created_node_pt != 0)
                {
                  // Now the node must be on a boundary, but we don't know which
                  // one. The returned created_node_pt is actually the
                  // neighbouring periodic node
                  Node* neighbour_node_pt = created_node_pt;

                  // Determine the edge on which the new node will live
                  int father_bound = Father_bound[n_p](jnod, son_type);

                  // Storage for the set of Mesh boundaries on which the
                  // appropriate father edge lives.
                  std::set<unsigned> boundaries;

                  // Only get the boundaries if we are at the edge of
                  // an element. Nodes in the centre of an element cannot be
                  // on Mesh boundaries
                  if (father_bound != Tree::OMEGA)
                  {
                    father_el_pt->get_boundaries(father_bound, boundaries);
                  }

#ifdef PARANOID
                  // Case where a new node lives on more than one boundary
                  // seems fishy enough to flag
                  if (boundaries.size() > 2)
                  {
                    throw OomphLibError(
                      "boundaries.size()>2 seems a bit strange..\n",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
                  }
#endif

                  // If the node is periodic and is definitely a boundary node.
                  // NOTE: the reason for this is that confusion can arise when
                  // a node is created on an edge that joins a periodic face.
                  if ((is_periodic) && (boundaries.size() > 0))
                  {
                    // Create node and set the pointer to it from the element
                    created_node_pt =
                      construct_boundary_node(jnod, time_stepper_pt);

                    // Make the node periodic from the neighbour
                    created_node_pt->make_periodic(neighbour_node_pt);

                    // Add to vector of new nodes
                    new_node_pt.push_back(created_node_pt);

                    // Loop over # of history values
                    for (unsigned t = 0; t < ntstorage; t++)
                    {
                      // Get position from father element -- this uses the macro
                      // element representation if appropriate. If the node
                      // turns out to be a hanging node later on, then
                      // its position gets adjusted in line with its
                      // hanging node interpolation.
                      Vector<double> x_prev(n_dim);
                      father_el_pt->get_x(t, s, x_prev);
                      // Set previous positions of the new node
                      for (unsigned i = 0; i < n_dim; i++)
                      {
                        created_node_pt->x(t, i) = x_prev[i];
                      }
                    }

                    // Next, we Update the boundary lookup schemes
                    // Loop over the boundaries stored in the set
                    for (std::set<unsigned>::iterator it = boundaries.begin();
                         it != boundaries.end();
                         ++it)
                    {
                      // Add the node to the boundary
                      mesh_pt->add_boundary_node(*it, created_node_pt);

                      // If we have set an intrinsic coordinate on this
                      // mesh boundary then it must also be interpolated on
                      // the new node
                      // Now interpolate the intrinsic boundary coordinate
                      if (mesh_pt->boundary_coordinate_exists(*it) == true)
                      {
                        Vector<double> zeta(2, 0.0);
                        father_el_pt->interpolated_zeta_on_face(
                          *it, father_bound, s, zeta);
                        created_node_pt->set_coordinates_on_boundary(*it, zeta);
                      }
                    }

                    // Make sure that we add the node to the mesh
                    mesh_pt->add_node_pt(created_node_pt);
                  } // End of periodic case
                  // Otherwise the node is not periodic, so just set the
                  // pointer to the neighbours node
                  else
                  {
                    node_pt(jnod) = created_node_pt;
                  }
                  node_done = true;
                }
                // Node does not exist in neighbour element but might already
                //-----------------------------------------------------------
                // have been created by a son of a neighbouring element
                //-----------------------------------------------------
                else
                {
                  // Was the node created by one of its neighbours' sons
                  // Whether or not the node lies on an edge can be calculated
                  // by from the fractional position
                  bool is_periodic = false;
                  created_node_pt =
                    node_created_by_son_of_neighbour(s_fraction, is_periodic);

                  // If the node was so created, assign the pointers
                  if (created_node_pt != 0)
                  {
                    // If the node is periodic
                    if (is_periodic)
                    {
                      // Now the node must be on a boundary, but we don't know
                      // which one The returned created_node_pt is actually the
                      // neighbouring periodic node
                      Node* neighbour_node_pt = created_node_pt;

                      // Determine the edge on which the new node will live
                      int father_bound = Father_bound[n_p](jnod, son_type);

                      // Storage for the set of Mesh boundaries on which the
                      // appropriate father edge lives.
                      std::set<unsigned> boundaries;

                      // Only get the boundaries if we are at the edge of
                      // an element. Nodes in the centre of an element cannot be
                      // on Mesh boundaries
                      if (father_bound != Tree::OMEGA)
                      {
                        father_el_pt->get_boundaries(father_bound, boundaries);
                      }

#ifdef PARANOID
                      // Case where a new node lives on more than one boundary
                      // seems fishy enough to flag
                      if (boundaries.size() > 2)
                      {
                        throw OomphLibError(
                          "boundaries.size()>2 seems a bit strange..\n",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
                      }

                      // Case when there are no boundaries, we are in big
                      // trouble
                      if (boundaries.size() == 0)
                      {
                        std::ostringstream error_stream;
                        error_stream
                          << "Periodic node is not on a boundary...\n"
                          << "Coordinates: " << created_node_pt->x(0) << " "
                          << created_node_pt->x(1) << "\n";
                        throw OomphLibError(error_stream.str(),
                                            OOMPH_CURRENT_FUNCTION,
                                            OOMPH_EXCEPTION_LOCATION);
                      }
#endif

                      // Create node and set the pointer to it from the element
                      created_node_pt =
                        construct_boundary_node(jnod, time_stepper_pt);

                      // Make the node periodic from the neighbour
                      created_node_pt->make_periodic(neighbour_node_pt);

                      // Add to vector of new nodes
                      new_node_pt.push_back(created_node_pt);

                      // Loop over # of history values
                      for (unsigned t = 0; t < ntstorage; t++)
                      {
                        // Get position from father element -- this uses the
                        // macro element representation if appropriate. If the
                        // node turns out to be a hanging node later on, then
                        // its position gets adjusted in line with its
                        // hanging node interpolation.
                        Vector<double> x_prev(n_dim, 0.0);
                        father_el_pt->get_x(t, s, x_prev);
                        // Set previous positions of the new node
                        for (unsigned i = 0; i < n_dim; i++)
                        {
                          created_node_pt->x(t, i) = x_prev[i];
                        }
                      }

                      // Next, we Update the boundary lookup schemes
                      // Loop over the boundaries stored in the set
                      for (std::set<unsigned>::iterator it = boundaries.begin();
                           it != boundaries.end();
                           ++it)
                      {
                        // Add the node to the boundary
                        mesh_pt->add_boundary_node(*it, created_node_pt);

                        // If we have set an intrinsic coordinate on this
                        // mesh boundary then it must also be interpolated on
                        // the new node
                        // Now interpolate the intrinsic boundary coordinate
                        if (mesh_pt->boundary_coordinate_exists(*it) == true)
                        {
                          Vector<double> zeta(2, 0.0);
                          father_el_pt->interpolated_zeta_on_face(
                            *it, father_bound, s, zeta);
                          created_node_pt->set_coordinates_on_boundary(*it,
                                                                       zeta);
                        }
                      }

                      // Make sure that we add the node to the mesh
                      mesh_pt->add_node_pt(created_node_pt);
                    } // End of periodic case
                    // Otherwise the node is not periodic, so just set the
                    // pointer to the neighbours node
                    else
                    {
                      node_pt(jnod) = created_node_pt;
                    }
                    // Node has been created
                    node_done = true;
                  } // Node does not exist in son of neighbouring element
                } // Node does not exist in neighbouring element
              } // Node does not exist in father element

              // Node already exists in father: No need to do anything else!
              // otherwise deal with boundary information etc.
              if (!node_exists_in_father)
              {
                // Check boundary status
                //----------------------

                // Firstly, we need to determine whether or not a node lies
                // on the boundary before building it, because
                // we actually assign a different type of node on boundaries.

                // If the new node lives on a face that is
                // shared with a face of its father element,
                // it needs to inherit the bounday conditions
                // from the father face
                int father_bound = Father_bound[n_p](jnod, son_type);

                // Storage for the set of Mesh boundaries on which the
                // appropriate father face lives.
                // [New nodes should always be mid-edge nodes in father
                // and therefore only live on one boundary but just to
                // play it safe...]
                std::set<unsigned> boundaries;

                // Only get the boundaries if we are at the edge of
                // an element. Nodes in the centre of an element cannot be
                // on Mesh boundaries
                if (father_bound != Tree::OMEGA)
                {
                  father_el_pt->get_boundaries(father_bound, boundaries);
                }

#ifdef PARANOID
                // Case where a new node lives on more than two boundaries
                // seems fishy enough to flag
                if (boundaries.size() > 2)
                {
                  throw OomphLibError(
                    "boundaries.size()>2 seems a bit strange..\n",
                    OOMPH_CURRENT_FUNCTION,
                    OOMPH_EXCEPTION_LOCATION);
                }
#endif

                // If the node lives on a mesh boundary,
                // then we need to create a boundary node
                if (boundaries.size() > 0)
                {
                  // Do we need a new node?
                  if (!node_done)
                  {
                    // Create node and set the internal pointer
                    created_node_pt =
                      construct_boundary_node(jnod, time_stepper_pt);
                    // Add to vector of new nodes
                    new_node_pt.push_back(created_node_pt);
                  }

                  // Now we need to work out whether to pin the values at
                  // the new node based on the boundary conditions applied at
                  // its Mesh boundary

                  // Get the boundary conditions from the father.
                  // Note: We can only deal with the values that are
                  //       continuously interpolated in the bulk element.
                  unsigned n_cont = ncont_interpolated_values();
                  Vector<int> bound_cons(n_cont);
                  father_el_pt->get_bcs(father_bound, bound_cons);

                  // Loop over the continuously interpolated values and pin,
                  // if necessary
                  for (unsigned k = 0; k < n_cont; k++)
                  {
                    if (bound_cons[k])
                    {
                      created_node_pt->pin(k);
                    }
                  }

                  // Solid node? If so, deal with the positional boundary
                  // conditions:
                  SolidNode* solid_node_pt =
                    dynamic_cast<SolidNode*>(created_node_pt);
                  if (solid_node_pt != 0)
                  {
                    // Get the positional boundary conditions from the father:
                    unsigned n_dim = created_node_pt->ndim();
                    Vector<int> solid_bound_cons(n_dim);
                    RefineableSolidQElement<3>* father_solid_el_pt =
                      dynamic_cast<RefineableSolidQElement<3>*>(father_el_pt);
#ifdef PARANOID
                    if (father_solid_el_pt == 0)
                    {
                      std::string error_message = "We have a SolidNode outside "
                                                  "a refineable SolidElement\n";
                      error_message +=
                        "during mesh refinement -- this doesn't make sense";

                      throw OomphLibError(error_message,
                                          OOMPH_CURRENT_FUNCTION,
                                          OOMPH_EXCEPTION_LOCATION);
                    }
#endif
                    father_solid_el_pt->get_solid_bcs(father_bound,
                                                      solid_bound_cons);

                    // Loop over the positions and pin, if necessary
                    for (unsigned k = 0; k < n_dim; k++)
                    {
                      if (solid_bound_cons[k])
                      {
                        solid_node_pt->pin_position(k);
                      }
                    }
                  } // End of if solid_node_pt

                  // Next update the boundary look-up schemes

                  // Loop over the boundaries in the set
                  for (std::set<unsigned>::iterator it = boundaries.begin();
                       it != boundaries.end();
                       ++it)
                  {
                    // Add the node to the bounadry
                    mesh_pt->add_boundary_node(*it, created_node_pt);

                    // If we have set an intrinsic coordinate on this
                    // mesh boundary then it must also be interpolated on
                    // the new node
                    // Now interpolate the intrinsic boundary coordinate
                    if (mesh_pt->boundary_coordinate_exists(*it) == true)
                    {
                      // Usually there will be two coordinates
                      Vector<double> zeta(2);
                      father_el_pt->interpolated_zeta_on_face(
                        *it, father_bound, s, zeta);

                      created_node_pt->set_coordinates_on_boundary(*it, zeta);
                    }
                  }
                }
                // Otherwise the node is not on a Mesh boundary and
                // we create a normal "bulk" node
                else
                {
                  // Do we need a new node?
                  if (!node_done)
                  {
                    // Create node and set the pointer to it from the element
                    created_node_pt = construct_node(jnod, time_stepper_pt);
                    // Add to vector of new nodes
                    new_node_pt.push_back(created_node_pt);
                  }
                }

                // In the first instance use macro element or FE representation
                // to create past and present nodal positions.
                // (THIS STEP SHOULD NOT BE SKIPPED FOR ALGEBRAIC
                // ELEMENTS AS NOT ALL OF THEM NECESSARILY IMPLEMENT
                // NONTRIVIAL NODE UPDATE FUNCTIONS. CALLING
                // THE NODE UPDATE FOR SUCH ELEMENTS/NODES WILL LEAVE
                // THEIR NODAL POSITIONS WHERE THEY WERE (THIS IS APPROPRIATE
                // ONCE THEY HAVE BEEN GIVEN POSITIONS) BUT WILL
                // NOT ASSIGN SENSIBLE INITIAL POSITONS!


                // Have we created a new node?
                if (!node_done)
                {
                  // Loop over # of history values
                  for (unsigned t = 0; t < ntstorage; t++)
                  {
                    // Get position from father element -- this uses the macro
                    // element representation if appropriate. If the node
                    // turns out to be a hanging node later on, then
                    // its position gets adjusted in line with its
                    // hanging node interpolation.
                    Vector<double> x_prev(n_dim);
                    father_el_pt->get_x(t, s, x_prev);

                    // Set previous positions of the new node
                    for (unsigned i = 0; i < n_dim; i++)
                    {
                      created_node_pt->x(t, i) = x_prev[i];
                    }
                  }

                  // Now set the values
                  // Loop over all history values
                  for (unsigned t = 0; t < ntstorage; t++)
                  {
                    // Get values from father element
                    // Note: get_interpolated_values() sets Vector size itself.
                    Vector<double> prev_values;
                    father_el_pt->get_interpolated_values(t, s, prev_values);

                    // Initialise the values at the new node
                    unsigned n_value = created_node_pt->nvalue();
                    for (unsigned k = 0; k < n_value; k++)
                    {
                      created_node_pt->set_value(t, k, prev_values[k]);
                    }
                  }

                  // Add new node to mesh
                  mesh_pt->add_node_pt(created_node_pt);

                } // End of whether the node has been created by us or not

              } // End of if node already existed in father in which case
              // everything above gets bypassed.

              // Check if the element is an algebraic element
              AlgebraicElementBase* alg_el_pt =
                dynamic_cast<AlgebraicElementBase*>(this);

              // If the element is an algebraic element, setup
              // node position (past and present) from algebraic update
              // function.
              // NOTE: YES, THIS NEEDS TO BE CALLED REPEATEDLY IF THE
              // NODE IS MEMBER OF MULTIPLE ELEMENTS: THEY ALL ASSIGN
              // THE SAME NODAL POSITIONS BUT WE NEED TO ADD THE REMESH
              // INFO FOR *ALL* ROOT ELEMENTS!
              if (alg_el_pt != 0)
              {
                // Build algebraic node update info for new node
                // This sets up the node update data for all node update
                // functions that are shared by all nodes in the father
                // element
                alg_el_pt->setup_algebraic_node_update(
                  node_pt(jnod), s, father_el_pt);
              }

              // We have built the node and we are documenting
              if ((!node_done) && (new_nodes_file.is_open()))
              {
                new_nodes_file << node_pt(jnod)->x(0) << " "
                               << node_pt(jnod)->x(1) << " "
                               << node_pt(jnod)->x(2) << std::endl;
              }

            } // End of Z loop over nodes in element

          } // End of vertical loop over nodes in element

        } // End of horizontal loop over nodes in element

        // If the element is a MacroElementNodeUpdateElement, set
        // the update parameters for the current element's nodes --
        // all this needs is the vector of (pointers to the)
        // geometric objects that affect the MacroElement-based
        // node update -- this is the same as that in the father element
        MacroElementNodeUpdateElementBase* father_m_el_pt =
          dynamic_cast<MacroElementNodeUpdateElementBase*>(father_el_pt);
        if (father_m_el_pt != 0)
        {
          // Get vector of geometric objects from father (construct vector
          // via copy operation)
          Vector<GeomObject*> geom_object_pt(father_m_el_pt->geom_object_pt());

          // Cast current element to MacroElementNodeUpdateElement:
          MacroElementNodeUpdateElementBase* m_el_pt =
            dynamic_cast<MacroElementNodeUpdateElementBase*>(this);

#ifdef PARANOID
          if (m_el_pt == 0)
          {
            std::string error_message =
              "Failed to cast to MacroElementNodeUpdateElementBase*\n";
            error_message +=
              "Strange -- if the father is a MacroElementNodeUpdateElement\n";
            error_message += "the son should be too....\n";

            throw OomphLibError(
              error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
          }
#endif
          // Build update info by passing vector of geometric objects:
          // This sets the current element to be the update element
          // for all of the element's nodes -- this is reversed
          // if the element is ever un-refined in the father element's
          // rebuild_from_sons() function which overwrites this
          // assignment to avoid nasty segmentation faults that occur
          // when a node tries to update itself via an element that no
          // longer exists...
          m_el_pt->set_node_update_info(geom_object_pt);
        }

#ifdef OOMPH_HAS_MPI
        // Is the new element a halo element?
        if (tree_pt()->father_pt()->object_pt()->is_halo())
        {
          Non_halo_proc_ID =
            tree_pt()->father_pt()->object_pt()->non_halo_proc_ID();
        }
#endif

        // Is it an ElementWithMovingNodes?
        ElementWithMovingNodes* aux_el_pt =
          dynamic_cast<ElementWithMovingNodes*>(this);

        // Pass down the information re the method for the evaluation
        // of the shape derivatives
        if (aux_el_pt != 0)
        {
          ElementWithMovingNodes* aux_father_el_pt =
            dynamic_cast<ElementWithMovingNodes*>(father_el_pt);

#ifdef PARANOID
          if (aux_father_el_pt == 0)
          {
            std::string error_message =
              "Failed to cast to ElementWithMovingNodes*\n";
            error_message +=
              "Strange -- if the son is a ElementWithMovingNodes\n";
            error_message += "the father should be too....\n";

            throw OomphLibError(
              error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
          }
#endif

          // If evaluating the residuals by finite differences in the father
          // continue to do so in the child
          if (aux_father_el_pt
                ->are_dresidual_dnodal_coordinates_always_evaluated_by_fd())
          {
            aux_el_pt
              ->enable_always_evaluate_dresidual_dnodal_coordinates_by_fd();
          }


          aux_el_pt->method_for_shape_derivs() =
            aux_father_el_pt->method_for_shape_derivs();

          // If bypassing the evaluation of fill_in_jacobian_from_geometric_data
          // continue to do so
          if (aux_father_el_pt
                ->is_fill_in_jacobian_from_geometric_data_bypassed())
          {
            aux_el_pt->enable_bypass_fill_in_jacobian_from_geometric_data();
          }
        }

        // Now do further build (if any)
        further_build();

      } // Sanity check: Father element has been generated

    } // End for element has not been built yet
    else
    {
      was_already_built = true;
    }
  }

  //====================================================================
  /// Set up all hanging nodes.
  ///
  /// (Wrapper to avoid specification of the unwanted output file).
  //====================================================================
  void RefineableQElement<3>::setup_hanging_nodes(
    Vector<std::ofstream*>& output_stream)
  {
#ifdef PARANOID
    if (output_stream.size() != 6)
    {
      throw OomphLibError("There must be six output streams",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    using namespace OcTreeNames;

    // Setup hanging nodes on each edge of the element
    oc_hang_helper(-1, U, *(output_stream[0]));
    oc_hang_helper(-1, D, *(output_stream[1]));
    oc_hang_helper(-1, L, *(output_stream[2]));
    oc_hang_helper(-1, R, *(output_stream[3]));
    oc_hang_helper(-1, B, *(output_stream[4]));
    oc_hang_helper(-1, F, *(output_stream[5]));
  }

  //================================================================
  /// Internal function that sets up the hanging node scheme for a
  /// particular value
  //=================================================================
  void RefineableQElement<3>::setup_hang_for_value(const int& value_id)
  {
    using namespace OcTreeNames;

    std::ofstream dummy_hangfile;
    // Setup hanging nodes on each edge of the element
    oc_hang_helper(value_id, U, dummy_hangfile);
    oc_hang_helper(value_id, D, dummy_hangfile);
    oc_hang_helper(value_id, R, dummy_hangfile);
    oc_hang_helper(value_id, L, dummy_hangfile);
    oc_hang_helper(value_id, B, dummy_hangfile);
    oc_hang_helper(value_id, F, dummy_hangfile);
  }

  //=================================================================
  /// Internal function to set up the hanging nodes on a particular
  /// face of the element
  //=================================================================
  void RefineableQElement<3>::oc_hang_helper(const int& value_id,
                                             const int& my_face,
                                             std::ofstream& output_hangfile)
  {
    using namespace OcTreeNames;

    // Number of dimensions
    unsigned n_dim = 3;

    Vector<unsigned> translate_s(n_dim);
    Vector<double> s_lo_neigh(n_dim);
    Vector<double> s_hi_neigh(n_dim);
    int neigh_face;
    int diff_level;
    bool in_neighbouring_tree;

    // Find pointer to neighbour in this direction
    OcTree* neigh_pt;
    neigh_pt = octree_pt()->gteq_face_neighbour(my_face,
                                                translate_s,
                                                s_lo_neigh,
                                                s_hi_neigh,
                                                neigh_face,
                                                diff_level,
                                                in_neighbouring_tree);

    // Neighbour exists
    if (neigh_pt != 0)
    {
      // Different sized element?
      if (diff_level != 0)
      {
        // Test for the periodic node case
        // Are we crossing a periodic boundary
        bool is_periodic = false;
        if (in_neighbouring_tree)
        {
          is_periodic = tree_pt()->root_pt()->is_neighbour_periodic(my_face);
        }

        // If it is periodic we actually need to get the node in
        // the neighbour of the neighbour (which will be a parent of
        // the present element) so that the "fixed" coordinate is
        // correctly calculated.
        // The idea is to replace the neigh_pt and associated data
        // with those of the neighbour of the neighbour
        if (is_periodic)
        {
          // Required data for the neighbour finding routine
          Vector<unsigned> translate_s_in_neigh(n_dim, 0.0);
          Vector<double> s_lo_neigh_of_neigh(n_dim, 0.0);
          Vector<double> s_hi_neigh_of_neigh(n_dim, 0.0);
          int neigh_face_of_neigh;
          int diff_level_of_neigh;
          bool in_neighbouring_tree_of_neigh;

          // Find pointer to neighbour of the neighbour on the edge
          // that we are currently considering
          OcTree* neigh_of_neigh_pt;
          neigh_of_neigh_pt =
            neigh_pt->gteq_face_neighbour(neigh_face,
                                          translate_s_in_neigh,
                                          s_lo_neigh_of_neigh,
                                          s_hi_neigh_of_neigh,
                                          neigh_face_of_neigh,
                                          diff_level_of_neigh,
                                          in_neighbouring_tree_of_neigh);

          // Set the value of the NEW neighbour and edge
          neigh_pt = neigh_of_neigh_pt;
          neigh_face = neigh_face_of_neigh;

          // Set the values of the translation schemes
          // Need to find the values of s_lo and s_hi
          // in the neighbour of the neighbour

          // Get the minimum and maximum values of the coordinate
          // in the neighbour (don't like this, but I think it's
          // necessary) Note that these values are hardcoded
          // in the quadtrees at some point!!
          double s_min = neigh_pt->object_pt()->s_min();
          double s_max = neigh_pt->object_pt()->s_max();
          Vector<double> s_lo_frac(n_dim), s_hi_frac(n_dim);
          // Work out the fractional position of the low and high points
          // of the original element
          for (unsigned i = 0; i < n_dim; i++)
          {
            s_lo_frac[i] = (s_lo_neigh[i] - s_min) / (s_max - s_min);
            s_hi_frac[i] = (s_hi_neigh[i] - s_min) / (s_max - s_min);
          }

          // We should now be able to construct the low and high points in
          // the neighbour of the neighbour
          for (unsigned i = 0; i < n_dim; i++)
          {
            s_lo_neigh[i] = s_lo_neigh_of_neigh[i] +
                            s_lo_frac[translate_s_in_neigh[i]] *
                              (s_hi_neigh_of_neigh[i] - s_lo_neigh_of_neigh[i]);
            s_hi_neigh[i] = s_lo_neigh_of_neigh[i] +
                            s_hi_frac[translate_s_in_neigh[i]] *
                              (s_hi_neigh_of_neigh[i] - s_lo_neigh_of_neigh[i]);
          }

          // Finally we must sort out the translation scheme
          Vector<unsigned> temp_translate(n_dim, 0.0);
          for (unsigned i = 0; i < n_dim; i++)
          {
            temp_translate[i] = translate_s_in_neigh[translate_s[i]];
          }
          for (unsigned i = 0; i < n_dim; i++)
          {
            translate_s[i] = temp_translate[i];
          }
        } // End of special treatment for periodic hanging nodes

        // Number of nodes in one dimension
        unsigned n_p = ninterpolating_node_1d(value_id);
        // Storage for the local nodes along the face of the element
        Node* local_node_pt = 0;

        // Loop over nodes along the face
        for (unsigned i0 = 0; i0 < n_p; i0++)
        {
          for (unsigned i1 = 0; i1 < n_p; i1++)
          {
            // Storage for the fractional position of the node in the element
            Vector<double> s_fraction(n_dim);

            // Local node number
            switch (my_face)
            {
              case U:
                s_fraction[0] =
                  local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
                s_fraction[1] = 1.0;
                s_fraction[2] =
                  local_one_d_fraction_of_interpolating_node(i1, 2, value_id);
                local_node_pt = interpolating_node_pt(
                  i0 + (n_p - 1) * n_p + n_p * n_p * i1, value_id);
                break;

              case D:
                s_fraction[0] =
                  local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
                s_fraction[1] = 0.0;
                s_fraction[2] =
                  local_one_d_fraction_of_interpolating_node(i1, 2, value_id);
                local_node_pt =
                  interpolating_node_pt(i0 + n_p * n_p * i1, value_id);
                break;

              case R:
                s_fraction[0] = 1.0;
                s_fraction[1] =
                  local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
                s_fraction[2] =
                  local_one_d_fraction_of_interpolating_node(i1, 2, value_id);
                local_node_pt = interpolating_node_pt(
                  n_p - 1 + i0 * n_p + i1 * n_p * n_p, value_id);
                break;

              case L:
                s_fraction[0] = 0.0;
                s_fraction[1] =
                  local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
                s_fraction[2] =
                  local_one_d_fraction_of_interpolating_node(i1, 2, value_id);
                local_node_pt =
                  interpolating_node_pt(n_p * i0 + i1 * n_p * n_p, value_id);
                break;

              case B:
                s_fraction[0] =
                  local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
                s_fraction[1] =
                  local_one_d_fraction_of_interpolating_node(i1, 1, value_id);
                s_fraction[2] = 0.0;
                local_node_pt = interpolating_node_pt(i0 + i1 * n_p, value_id);
                break;

              case F:
                s_fraction[0] =
                  local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
                s_fraction[1] =
                  local_one_d_fraction_of_interpolating_node(i1, 1, value_id);
                s_fraction[2] = 1.0;
                local_node_pt = interpolating_node_pt(
                  i0 + i1 * n_p + (n_p - 1) * n_p * n_p, value_id);
                break;

              default:
                throw OomphLibError("my_face not U, D, L, R, B, F\n",
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
            }

            // Set up the coordinate in the neighbour
            Vector<double> s_in_neighb(n_dim);
            for (unsigned i = 0; i < n_dim; i++)
            {
              s_in_neighb[i] =
                s_lo_neigh[i] +
                s_fraction[translate_s[i]] * (s_hi_neigh[i] - s_lo_neigh[i]);
            }

            // Find the Node in the neighbouring element
            Node* neighbouring_node_pt =
              neigh_pt->object_pt()->get_interpolating_node_at_local_coordinate(
                s_in_neighb, value_id);

            // If the neighbour does not have a node at this point
            if (0 == neighbouring_node_pt)
            {
              // Do we need to make a hanging node, we assume that we don't
              // initially
              bool make_hanging_node = false;

              // If the node is not hanging geometrically, then we must make it
              // hang
              if (!local_node_pt->is_hanging())
              {
                make_hanging_node = true;
              }
              // Otherwise is could be hanging geometrically, but still require
              // a different hanging scheme for this data value
              else
              {
                if ((value_id != -1) && (local_node_pt->hanging_pt(value_id) ==
                                         local_node_pt->hanging_pt()))
                {
                  make_hanging_node = true;
                }
              }

              if (make_hanging_node == true)
              {
                // Cache refineable element used here
                RefineableElement* const obj_pt = neigh_pt->object_pt();

                // Get shape functions in neighbour element
                Shape psi(obj_pt->ninterpolating_node(value_id));
                obj_pt->interpolating_basis(s_in_neighb, psi, value_id);

                // Allocate the storage for the Hang pointer
                // We know that it will hold n_p*n_p nodes
                HangInfo* hang_pt = new HangInfo(n_p * n_p);

                // Loop over nodes on edge in neighbour and mark them as nodes
                // that this node depends on
                unsigned n_neighbour;

                // Number of nodes along edge in neighbour element
                for (unsigned ii0 = 0; ii0 < n_p; ii0++)
                {
                  for (unsigned ii1 = 0; ii1 < n_p; ii1++)
                  {
                    switch (neigh_face)
                    {
                      case U:
                        n_neighbour = ii0 + n_p * (n_p - 1) + ii1 * n_p * n_p;
                        break;

                      case D:
                        n_neighbour = ii0 + ii1 * n_p * n_p;
                        break;

                      case L:
                        n_neighbour = ii0 * n_p + ii1 * n_p * n_p;
                        break;

                      case R:
                        n_neighbour = (n_p - 1) + ii0 * n_p + ii1 * n_p * n_p;
                        break;

                      case B:
                        n_neighbour = ii0 + ii1 * n_p;
                        break;

                      case F:
                        n_neighbour = ii0 + ii1 * n_p + n_p * n_p * (n_p - 1);
                        break;
                      default:
                        throw OomphLibError("neigh_face not U, L, R, B, F\n",
                                            OOMPH_CURRENT_FUNCTION,
                                            OOMPH_EXCEPTION_LOCATION);
                    }

                    // Push back neighbouring node and weight into
                    // Vector of (pointers to)
                    // master nodes and weights
                    // The weight is merely the value of the shape function
                    // corresponding to the node in the neighbour
                    hang_pt->set_master_node_pt(
                      ii0 * n_p + ii1,
                      obj_pt->interpolating_node_pt(n_neighbour, value_id),
                      psi[n_neighbour]);
                  }
                }
                // Now set the hanging data for the position
                // This also constrains the data values associated with the
                // value id
                local_node_pt->set_hanging_pt(hang_pt, value_id);
              }

              if (output_hangfile.is_open())
              {
                // output_hangfile
                output_hangfile << local_node_pt->x(0) << " "
                                << local_node_pt->x(1) << " "
                                << local_node_pt->x(2) << std::endl;
              }
            }
            else
            {
#ifdef PARANOID
              if (local_node_pt != neighbouring_node_pt)
              {
                std::ofstream reportage("dodgy.dat", std::ios_base::app);
                reportage << local_node_pt->x(0) << " " << local_node_pt->x(1)
                          << " " << local_node_pt->x(2) << std::endl;
                reportage.close();

                std::ostringstream warning_stream;
                warning_stream
                  << "SANITY CHECK in oc_hang_helper      \n"
                  << "Current node      " << local_node_pt << " at "
                  << "(" << local_node_pt->x(0) << ", " << local_node_pt->x(1)
                  << ", " << local_node_pt->x(2) << ")" << std::endl
                  << " is not hanging and has " << std::endl
                  << "Neighbour's node  " << neighbouring_node_pt << " at "
                  << "(" << neighbouring_node_pt->x(0) << ", "
                  << neighbouring_node_pt->x(1) << ", "
                  << neighbouring_node_pt->x(2) << ")" << std::endl
                  << "even though the two should be "
                  << "identical" << std::endl;
                OomphLibWarning(warning_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
              }
#endif
            }

            // If we are doing the position then
            if (value_id == -1)
            {
              // Get the nodal position from neighbour element
              Vector<double> x_in_neighb(n_dim);
              neigh_pt->object_pt()->interpolated_x(s_in_neighb, x_in_neighb);

              // Fine adjust the coordinates (macro map will pick up boundary
              // accurately but will lead to different element edges)
              local_node_pt->x(0) = x_in_neighb[0];
              local_node_pt->x(1) = x_in_neighb[1];
              local_node_pt->x(2) = x_in_neighb[2];
            }
          }
        }
      }
    }
  }


  //=================================================================
  /// Check inter-element continuity of
  /// - nodal positions
  /// - (nodally) interpolated function values
  //====================================================================
  void RefineableQElement<3>::check_integrity(double& max_error)
  {
    using namespace OcTreeNames;

    // std::ofstream error_file("errors.dat");

    // Number of nodes along edge
    unsigned n_p = nnode_1d();

    // Number of timesteps (incl. present) for which continuity is
    // to be checked.
    unsigned n_time = 1;

    // Initialise errors
    max_error = 0.0;
    Vector<double> max_error_x(3, 0.0);
    double max_error_val = 0.0;

    // Set the faces
    Vector<int> faces(6);
    faces[0] = D;
    faces[1] = U;
    faces[2] = L;
    faces[3] = R;
    faces[4] = B;
    faces[5] = F;

    // Loop over the edges
    for (unsigned face_counter = 0; face_counter < 6; face_counter++)
    {
      Vector<double> s(3), s_lo_neigh(3), s_hi_neigh(3);
      Vector<double> s_fraction(3);
      Vector<unsigned> translate_s(3);
      int neigh_face, diff_level;
      int my_face;
      bool in_neighbouring_tree;

      my_face = faces[face_counter];

      // Find pointer to neighbour in this direction
      OcTree* neigh_pt;
      neigh_pt = octree_pt()->gteq_face_neighbour(my_face,
                                                  translate_s,
                                                  s_lo_neigh,
                                                  s_hi_neigh,
                                                  neigh_face,
                                                  diff_level,
                                                  in_neighbouring_tree);

      // Neighbour exists and has existing nodes
      if ((neigh_pt != 0) && (neigh_pt->object_pt()->nodes_built()))
      {
        // if (diff_level!=0)
        {
          // Need to exclude periodic nodes from this check
          // There are only periodic nodes if we are in a neighbouring tree
          bool is_periodic = false;
          if (in_neighbouring_tree)
          {
            // Is it periodic
            is_periodic =
              tree_pt()->root_pt()->is_neighbour_periodic(faces[face_counter]);
          }

          // Loop over nodes along the edge
          for (unsigned i0 = 0; i0 < n_p; i0++)
          {
            for (unsigned i1 = 0; i1 < n_p; i1++)
            {
              // Local node number
              unsigned n = 0;
              switch (face_counter)
              {
                case 0:
                  // Fractions
                  s_fraction[0] = local_one_d_fraction_of_node(i0, 0);
                  s_fraction[1] = 0.0;
                  s_fraction[2] = local_one_d_fraction_of_node(i1, 2);
                  // Set local node number
                  n = i0 + i1 * n_p * n_p;
                  break;

                case 1:
                  s_fraction[0] = local_one_d_fraction_of_node(i0, 0);
                  s_fraction[1] = 1.0;
                  s_fraction[2] = local_one_d_fraction_of_node(i1, 2);
                  // Set local node number
                  n = i0 + n_p * (n_p - 1) + i1 * n_p * n_p;
                  break;

                case 2:
                  s_fraction[0] = 0.0;
                  s_fraction[1] = local_one_d_fraction_of_node(i0, 1);
                  s_fraction[2] = local_one_d_fraction_of_node(i1, 2);
                  // Set local node number
                  n = n_p * i0 + i1 * n_p * n_p;
                  break;

                case 3:
                  s_fraction[0] = 1.0;
                  s_fraction[1] = local_one_d_fraction_of_node(i0, 1);
                  s_fraction[2] = local_one_d_fraction_of_node(i1, 2);
                  // Set local node number
                  n = n_p - 1 + n_p * i0 + n_p * n_p * i1;
                  break;

                case 4:
                  s_fraction[0] = local_one_d_fraction_of_node(i0, 0);
                  s_fraction[1] = local_one_d_fraction_of_node(i1, 1);
                  s_fraction[2] = 0.0;
                  // Set local node number
                  n = i0 + n_p * i1;
                  break;

                case 5:
                  s_fraction[0] = local_one_d_fraction_of_node(i0, 0);
                  s_fraction[1] = local_one_d_fraction_of_node(i1, 1);
                  s_fraction[2] = 1.0;
                  // Set local node number
                  n = i0 + n_p * i1 + n_p * n_p * (n_p - 1);
                  break;
              }


              // We have to check if the hi and lo directions along the
              // face are inverted or not
              Vector<double> s_in_neighb(3);
              for (unsigned i = 0; i < 3; i++)
              {
                s[i] = -1.0 + 2.0 * s_fraction[i];
                s_in_neighb[i] =
                  s_lo_neigh[i] +
                  s_fraction[translate_s[i]] * (s_hi_neigh[i] - s_lo_neigh[i]);
              }


              // Loop over timesteps
              for (unsigned t = 0; t < n_time; t++)
              {
                // Get the nodal position from neighbour element
                Vector<double> x_in_neighb(3);
                neigh_pt->object_pt()->interpolated_x(
                  t, s_in_neighb, x_in_neighb);

                // Check error only if the node is NOT periodic
                if (is_periodic == false)
                {
                  // Check error
                  for (unsigned i = 0; i < 3; i++)
                  {
                    // Find the spatial error
                    double err =
                      std::fabs(node_pt(n)->x(t, i) - x_in_neighb[i]);

                    // If it's bigger than our tolerance, say so
                    if (err > 1e-9)
                    {
                      oomph_info << "errx[" << i << "], t x, x_neigh: " << err
                                 << " " << t << " " << node_pt(n)->x(t, i)
                                 << " " << x_in_neighb[i] << std::endl;
                      oomph_info << "at " << node_pt(n)->x(0) << " "
                                 << node_pt(n)->x(1) << " " << node_pt(n)->x(2)
                                 << " " << std::endl;
                    }

                    // If it's bigger than the previous max error, it is the
                    // new max error!
                    if (err > max_error_x[i])
                    {
                      max_error_x[i] = err;
                    }
                  }
                }

                // Get the values from neighbour element. Note: # of values
                // gets set by routine (because in general we don't know
                // how many interpolated values a certain element has
                Vector<double> values_in_neighb;
                neigh_pt->object_pt()->get_interpolated_values(
                  t, s_in_neighb, values_in_neighb);

                // Get the values in current element.
                Vector<double> values;
                get_interpolated_values(t, s, values);

                // Now figure out how many continuously interpolated
                // values there are
                unsigned num_val =
                  neigh_pt->object_pt()->ncont_interpolated_values();

                // Check error
                for (unsigned ival = 0; ival < num_val; ival++)
                {
                  double err = std::fabs(values[ival] - values_in_neighb[ival]);

                  if (err > 1.0e-10)
                  {
                    oomph_info << node_pt(n)->x(0) << " " << node_pt(n)->x(1)
                               << " " << node_pt(n)->x(2) << " \n# "
                               << "erru (S)" << err << " " << ival << " " << n
                               << " " << values[ival] << " "
                               << values_in_neighb[ival] << std::endl;
                    // error_file<<"ZONE"<<std::endl
                    //          <<  node_pt(n)->x(0) << " "
                    //          <<  node_pt(n)->x(1) << " "
                    //          <<  node_pt(n)->x(2) << std::endl;
                  }

                  if (err > max_error_val)
                  {
                    max_error_val = err;
                  }
                }
              }
            }
          }
        }
      }
    }

    max_error = max_error_x[0];
    if (max_error_x[1] > max_error) max_error = max_error_x[1];
    if (max_error_x[2] > max_error) max_error = max_error_x[2];
    if (max_error_val > max_error) max_error = max_error_val;

    if (max_error > 1e-9)
    {
      oomph_info << "\n#------------------------------------ \n#Max error ";
      oomph_info << max_error_x[0] << " " << max_error_x[1] << " "
                 << max_error_x[2] << " " << max_error_val << std::endl;
      oomph_info << "#------------------------------------ \n " << std::endl;
    }

    // error_file.close();
  }


  //==================================================================
  /// Determine vector of solid (positional) boundary conditions along
  /// the element's boundary (or vertex) bound (S/W/N/E/SW/SE/NW/NE).
  ///
  /// This function assumes that the same boundary condition is applied
  /// along the entire length of an element's edge (of course, the
  /// vertices combine the boundary conditions of their two adjacent edges
  /// in the most restrictive combination. Hence, if we're at a vertex,
  /// we apply the most restrictive boundary condition of the
  /// two adjacent edges. If we're on an edge (in its proper interior),
  /// we apply the least restrictive boundary condition of all nodes
  /// along the edge.
  ///
  /// Usual convention:
  ///   - solid_bound_cons[i]=0 if displacement in coordinate direction i
  ///                           on this boundary is free.
  ///   - solid_bound_cons[i]=1 if it's pinned.
  //==================================================================
  void RefineableSolidQElement<3>::get_solid_bcs(
    int bound, Vector<int>& solid_bound_cons) const
  {
    using namespace OcTreeNames;

    // Spatial dimension of all nodes
    unsigned n_dim = this->nodal_dimension();
    // Set up temporary vectors to hold edge boundary conditions
    Vector<int> bound_cons1(n_dim), bound_cons2(n_dim);
    Vector<int> bound_cons3(n_dim);

    Vector<int> vect1(3), vect2(3), vect3(3);
    Vector<int> vect_elem;
    Vector<int> notzero;
    int n = 0;

    vect_elem = OcTree::Direction_to_vector[bound];

    // Just to see if bound is a face, an edge, or a vertex, n stores
    // the number of non-zero values in the vector reprensenting the bound
    // and the vector notzero stores the position of these values
    for (int i = 0; i < 3; i++)
    {
      if (vect_elem[i] != 0)
      {
        n++;
        notzero.push_back(i);
      }
    }

    switch (n)
    {
        // If there is only one non-zero value, bound is a face
      case 1:
        get_face_solid_bcs(bound, solid_bound_cons);
        break;

        // If there are two non-zero values, bound is an edge
      case 2:

        for (unsigned i = 0; i < 3; i++)
        {
          vect1[i] = 0;
          vect2[i] = 0;
        }
        // vect1 and vect2 are the vector of the two faces adjacent to bound
        vect1[notzero[0]] = vect_elem[notzero[0]];
        vect2[notzero[1]] = vect_elem[notzero[1]];

        get_face_solid_bcs(OcTree::Vector_to_direction[vect1], bound_cons1);
        get_face_solid_bcs(OcTree::Vector_to_direction[vect2], bound_cons2);
        // get the most restrictive bc
        for (unsigned k = 0; k < n_dim; k++)
        {
          solid_bound_cons[k] = (bound_cons1[k] || bound_cons2[k]);
        }
        break;

        // If there are three non-zero value, bound is a vertex
      case 3:

        for (unsigned i = 0; i < 3; i++)
        {
          vect1[i] = 0;
          vect2[i] = 0;
          vect3[i] = 0;
        }
        // vectors to the three adjacent faces of the vertex
        vect1[0] = vect_elem[0];
        vect2[1] = vect_elem[1];
        vect3[2] = vect_elem[2];

        get_face_solid_bcs(OcTree::Vector_to_direction[vect1], bound_cons1);
        get_face_solid_bcs(OcTree::Vector_to_direction[vect2], bound_cons2);
        get_face_solid_bcs(OcTree::Vector_to_direction[vect3], bound_cons3);


        // set the bcs to the most restrictive ones
        for (unsigned k = 0; k < n_dim; k++)
        {
          solid_bound_cons[k] =
            (bound_cons1[k] || bound_cons2[k] || bound_cons3[k]);
        }
        break;

      default:
        throw OomphLibError("Make sure you are not giving  OMEGA as bound",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }

  //==================================================================
  /// Determine Vector of solid (positional) boundary conditions along
  /// the element's edge (S/N/W/E) -- BC is the least restrictive combination
  /// of all the nodes on this edge
  ///
  /// Usual convention:
  ///   - solid_bound_cons[i]=0 if displacement in coordinate direction i
  ///                           on this boundary is free
  ///   - solid_bound_cons[i]=1 if it's pinned
  //==================================================================
  void RefineableSolidQElement<3>::get_face_solid_bcs(
    const int& face, Vector<int>& solid_bound_cons) const
  {
    using namespace OcTreeNames;

    // Number of nodes along 1D edge
    unsigned n_p = nnode_1d();
    // the four corner nodes on the boundary
    unsigned node1, node2, node3, node4;

    // Set the four corner nodes for the face
    switch (face)
    {
      case U:
        node1 = n_p * n_p * n_p - 1;
        node2 = n_p * n_p - 1;
        node3 = n_p * (n_p - 1);
        node4 = n_p * (n_p * n_p - 1);
        break;

      case D:
        node1 = 0;
        node2 = n_p - 1;
        node3 = (n_p * n_p + 1) * (n_p - 1);
        node4 = n_p * n_p * (n_p - 1);
        break;

      case R:
        node1 = n_p - 1;
        node2 = (n_p * n_p + 1) * (n_p - 1);
        node3 = n_p * n_p * n_p - 1;
        node4 = n_p * n_p - 1;
        break;

      case L:
        node1 = 0;
        node2 = n_p * (n_p - 1);
        node3 = n_p * (n_p * n_p - 1);
        node4 = n_p * n_p * (n_p - 1);
        break;

      case B:
        node1 = 0;
        node2 = n_p - 1;
        node3 = n_p * n_p - 1;
        node4 = n_p * (n_p - 1);
        break;

      case F:
        node1 = n_p * n_p * n_p - 1;
        node2 = n_p * (n_p * n_p - 1);
        node3 = n_p * n_p * (n_p - 1);
        node4 = (n_p - 1) * (n_p * n_p + 1);
        break;

      default:
        std::ostringstream error_stream;
        error_stream << "Wrong edge " << face << " passed\n";

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Cast to solid nodes
    SolidNode* solid_node1_pt = dynamic_cast<SolidNode*>(node_pt(node1));
    SolidNode* solid_node2_pt = dynamic_cast<SolidNode*>(node_pt(node2));
    SolidNode* solid_node3_pt = dynamic_cast<SolidNode*>(node_pt(node3));
    SolidNode* solid_node4_pt = dynamic_cast<SolidNode*>(node_pt(node4));

    //#ifdef PARANOID
    if (solid_node1_pt == 0)
    {
      throw OomphLibError(
        "Corner node 1 cannot be cast to SolidNode --> something is wrong",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    if (solid_node2_pt == 0)
    {
      throw OomphLibError(
        "Corner node 2 cannot be cast to SolidNode --> something is wrong",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    if (solid_node3_pt == 0)
    {
      throw OomphLibError(
        "Corner node 3 cannot be cast to SolidNode --> something is wrong",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    if (solid_node4_pt == 0)
    {
      throw OomphLibError(
        "Corner node 4 cannot be cast to SolidNode --> something is wrong",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    //#endif

    // Number of coordinates
    unsigned n_dim = this->nodal_dimension();

    // Loop over all directions, the least restrictive value is
    // the multiplication of the boundary conditions at the 4 nodes
    // Assuming that free is always zero and pinned is one
    for (unsigned k = 0; k < n_dim; k++)
    {
      solid_bound_cons[k] = solid_node1_pt->position_is_pinned(k) *
                            solid_node2_pt->position_is_pinned(k) *
                            solid_node3_pt->position_is_pinned(k) *
                            solid_node4_pt->position_is_pinned(k);
    }
  }


  //========================================================================
  /// Static matrix for coincidence between son nodal points and
  /// father boundaries
  ///
  //========================================================================
  std::map<unsigned, DenseMatrix<int>> RefineableQElement<3>::Father_bound;


} // namespace oomph
