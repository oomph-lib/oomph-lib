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
#include <algorithm>

#include "mesh.h"
#include "algebraic_elements.h"
#include "macro_element_node_update_element.h"
#include "refineable_quad_element.h"

namespace oomph
{
  //==================================================================
  /// Setup static matrix for coincidence between son nodal points and
  /// father boundaries:
  ///
  /// Father_boundd[nnode_1d](nnode_son,son_type)={SW/SE/NW/NE/S/E/N/W/OMEGA}
  ///
  /// so that node nnode_son in element of type son_type lies
  /// on boundary/vertex Father_boundd[nnode_1d](nnode_son,son_type) in its
  /// father element. If the node doesn't lie on a boundary
  /// the value is OMEGA.
  //==================================================================
  void RefineableQElement<2>::setup_father_bounds()
  {
    using namespace QuadTreeNames;

    // Find the number of nodes along a 1D edge
    unsigned n_p = nnode_1d();
    // Allocate space for the boundary information
    Father_bound[n_p].resize(n_p * n_p, 4);

    // Initialise: By default points are not on the boundary
    for (unsigned n = 0; n < n_p * n_p; n++)
    {
      for (unsigned ison = 0; ison < 4; ison++)
      {
        Father_bound[n_p](n, ison) = Tree::OMEGA;
      }
    }

    // Southwest son
    //--------------
    // SW node (0) is the SW node of the parent
    Father_bound[n_p](0, SW) = SW;
    // Southern boundary is the southern boundary of the parent
    for (unsigned n = 1; n < n_p; n++)
    {
      Father_bound[n_p](n, SW) = S;
    }
    // Western boundary is the western boundary of the parent
    for (unsigned n = 1; n < n_p; n++)
    {
      Father_bound[n_p](n_p * n, SW) = W;
    }
    // Other boundaries are in the interior

    // Northwest son
    //--------------
    // NW node (n_p*(n_p-1))is the NW node of the parent
    Father_bound[n_p](n_p * (n_p - 1), NW) = NW;
    // Northern boundary is the northern boundary of the parent
    for (unsigned n = 1; n < n_p; n++)
    {
      Father_bound[n_p](n_p * (n_p - 1) + n, NW) = N;
    }
    // Western boundary is the western boundary of the parent
    for (unsigned n = 0; n < (n_p - 1); n++)
    {
      Father_bound[n_p](n_p * n, NW) = W;
    }
    // Other boundaries are in the interior

    // Northeast son
    //--------------
    // NE node (n_p*n_p-1) is the NE node of the parent
    Father_bound[n_p](n_p * n_p - 1, NE) = NE;
    // Northern boundary is the northern boundary of the parent
    for (unsigned n = 0; n < (n_p - 1); n++)
    {
      Father_bound[n_p](n_p * (n_p - 1) + n, NE) = N;
    }
    // Eastern boundary is the eastern boundary of the parent
    for (unsigned n = 0; n < (n_p - 1); n++)
    {
      Father_bound[n_p](n_p - 1 + n_p * n, NE) = E;
    }
    // Other boundaries are in the interior

    // Southeast son
    //--------------
    // SE node (n_p-1) is the SE node of the parent
    Father_bound[n_p](n_p - 1, SE) = SE;
    // Southern boundary is the southern boundary of the parent
    for (unsigned n = 0; n < (n_p - 1); n++)
    {
      Father_bound[n_p](n, SE) = S;
    }
    // Eastern boundary is the eastern boundary of the parent
    for (unsigned n = 1; n < n_p; n++)
    {
      Father_bound[n_p](n_p - 1 + n_p * n, SE) = E;
    }
  }

  //==================================================================
  /// Determine Vector of boundary conditions along the element's boundary
  /// (or vertex) bound (S/W/N/E/SW/SE/NW/NE).
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
  ///   - bound_cons[ival]=0 if value ival on this boundary is free
  ///   - bound_cons[ival]=1 if value ival on this boundary is pinned
  //==================================================================
  void RefineableQElement<2>::get_bcs(int bound, Vector<int>& bound_cons) const
  {
    using namespace QuadTreeNames;

    // Max. number of nodal data values in element
    unsigned nvalue = ncont_interpolated_values();
    // Set up temporary vectors to hold edge boundary conditions
    Vector<int> bound_cons1(nvalue), bound_cons2(nvalue);

    // Which boundary are we on?
    switch (bound)
    {
        // If on edge, just get the bcs
      case N:
      case S:
      case W:
      case E:
        get_edge_bcs(bound, bound_cons);
        break;

        // Most restrictive boundary at SE corner
      case SE:
        get_edge_bcs(S, bound_cons1);
        get_edge_bcs(E, bound_cons2);

        for (unsigned k = 0; k < nvalue; k++)
        {
          bound_cons[k] = (bound_cons1[k] || bound_cons2[k]);
        }
        break;

        // Most restrictive boundary at SW corner
      case SW:
        get_edge_bcs(S, bound_cons1);
        get_edge_bcs(W, bound_cons2);

        for (unsigned k = 0; k < nvalue; k++)
        {
          bound_cons[k] = (bound_cons1[k] || bound_cons2[k]);
        }
        break;

        // Most restrictive boundary at NW corner
      case NW:
        get_edge_bcs(N, bound_cons1);
        get_edge_bcs(W, bound_cons2);

        for (unsigned k = 0; k < nvalue; k++)
        {
          bound_cons[k] = (bound_cons1[k] || bound_cons2[k]);
        }
        break;

        // Most restrictive boundary at NE corner
      case NE:
        get_edge_bcs(N, bound_cons1);
        get_edge_bcs(E, bound_cons2);

        for (unsigned k = 0; k < nvalue; k++)
        {
          bound_cons[k] = (bound_cons1[k] || bound_cons2[k]);
        }
        break;

      default:
        throw OomphLibError(
          "Wrong boundary", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  //==================================================================
  /// Determine Vector of boundary conditions along the element's
  /// edge (S/N/W/E) -- BC is the least restrictive combination
  /// of all the nodes on this edge
  ///
  /// Usual convention:
  ///   - bound_cons[ival]=0 if value ival on this boundary is free
  ///   - bound_cons[ival]=1 if value ival on this boundary is pinned
  //==================================================================
  void RefineableQElement<2>::get_edge_bcs(const int& edge,
                                           Vector<int>& bound_cons) const
  {
    using namespace QuadTreeNames;

    // Number of nodes along 1D edge
    unsigned n_p = nnode_1d();
    // Left- and Right-hand nodes
    unsigned left_node, right_node;

    // Set the left (lower) and right (upper) nodes for the edge
    switch (edge)
    {
      case N:
        left_node = n_p * (n_p - 1);
        right_node = n_p * n_p - 1;
        break;

      case S:
        left_node = 0;
        right_node = n_p - 1;
        break;

      case W:
        left_node = 0;
        right_node = n_p * (n_p - 1);
        break;

      case E:
        left_node = n_p - 1;
        right_node = n_p * n_p - 1;
        break;

      default:
        std::ostringstream error_stream;
        error_stream << "Wrong edge " << edge << " passed to get_edge_bcs(..)"
                     << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Max. number of nodal data values in element
    unsigned maxnvalue = ncont_interpolated_values();

    // Loop over all values, the least restrictive value is
    // the boundary condition at the left multiplied by that at the right
    // Assuming that free is always zero and pinned is one
    for (unsigned k = 0; k < maxnvalue; k++)
    {
      bound_cons[k] =
        node_pt(left_node)->is_pinned(k) * node_pt(right_node)->is_pinned(k);
    }
  }

  //==================================================================
  /// Given an element edge/vertex, return a set that contains
  /// all the (mesh-)boundary numbers that this element edge/vertex
  /// lives on.
  ///
  /// For proper edges, the boundary is the one (if any) that is shared by
  /// both vertex nodes). For vertex nodes, we just return their
  /// boundaries.
  //==================================================================
  void RefineableQElement<2>::get_boundaries(const int& edge,
                                             std::set<unsigned>& boundary) const
  {
    using namespace QuadTreeNames;

    // Number of 1d nodes along an edge
    unsigned n_p = nnode_1d();
    // Left and right-hand nodes
    int left_node = -1, right_node = -1;

    // Set the left (lower) and right (upper) nodes for the edge
    switch (edge)
    {
      case N:
        left_node = n_p * (n_p - 1);
        right_node = n_p * n_p - 1;
        break;

      case S:
        left_node = 0;
        right_node = n_p - 1;
        break;

      case W:
        left_node = 0;
        right_node = n_p * (n_p - 1);
        break;

      case E:
        left_node = n_p - 1;
        right_node = n_p * n_p - 1;
        break;

        // Vertices do not have left nodes!
      case SE:
        right_node = n_p - 1;
        break;

      case SW:
        right_node = 0;
        break;

      case NE:
        right_node = n_p * n_p - 1;
        break;

      case NW:
        right_node = n_p * (n_p - 1);
        break;

      default:
        std::ostringstream error_stream;
        error_stream << "Wrong edge " << edge << " passed" << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Empty the boundary set: Edge does not live on any boundary
    boundary.clear();

    // Storage for the nodes that the right node lives on
    std::set<unsigned>* right_boundaries_pt = 0;
    // Get the boundaries that the right node lives on
    node_pt(right_node)->get_boundaries_pt(right_boundaries_pt);

    // If the right node lives on some boundaries
    if (right_boundaries_pt != 0)
    {
      // If the node is a vertex then add the boundaries at the node
      // into the vector boundary
      if (left_node < 0)
      {
        copy(right_boundaries_pt->begin(),
             right_boundaries_pt->end(),
             inserter(boundary, boundary.begin()));
      }
      // Otherwise only add if the boundary also exists at the left hand node
      else
      {
        std::set<unsigned>* left_boundaries_pt = 0;
        node_pt(left_node)->get_boundaries_pt(left_boundaries_pt);
        // If the left node is on some boundaries
        if (left_boundaries_pt != 0)
        {
          // Use the standard algorithms to compute the boundaries in
          // common between the left and right nodes
          std::set_intersection(right_boundaries_pt->begin(),
                                right_boundaries_pt->end(),
                                left_boundaries_pt->begin(),
                                left_boundaries_pt->end(),
                                inserter(boundary, boundary.begin()));
        }
      }
    }
  }


  //===================================================================
  /// Return the value of the intrinsic boundary coordinate interpolated
  /// along the edge (S/W/N/E)
  //===================================================================
  void RefineableQElement<2>::interpolated_zeta_on_edge(
    const unsigned& boundary,
    const int& edge,
    const Vector<double>& s,
    Vector<double>& zeta)
  {
    using namespace QuadTreeNames;

    // Number of 1D nodes along an edge
    unsigned n_p = nnode_1d();

    // Storage for the shape functions
    Shape psi(n_p * n_p);
    // Get the shape functions at the passed position
    this->shape(s, psi);

    // Unsigned data that give starts and multipliers for the loop
    // over the nodes on the edges.
    unsigned start = 0, multiplier = 1;

    // Which edge?
    switch (edge)
    {
      case S:
#ifdef PARANOID
        if (s[1] != -1.0)
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1]
                       << " is not on South edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Start is zero and multiplier is one
        break;

      case N:
#ifdef PARANOID
        if (s[1] != 1.0)
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1]
                       << " is not on North edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Start from the top left corner of the element, multiplier still one
        start = n_p * (n_p - 1);
        break;

      case W:
#ifdef PARANOID
        if (s[0] != -1.0)
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1]
                       << " is not on West edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Loop over left-hand edge of element (start from zero)
        multiplier = n_p;
        break;

      case E:
#ifdef PARANOID
        if (s[0] != 1.0)
        {
          std::ostringstream error_stream;
          error_stream << "Coordinate " << s[0] << " " << s[1]
                       << " is not on East edge\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Start from the bottom right-hand corner
        start = n_p - 1;
        // Loop over the right-hand edge of the element
        multiplier = n_p;
        break;


      default:
        std::ostringstream error_stream;
        error_stream << "Wrong edge " << edge << " passed" << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Initialise the intrinsic coordinate
    double inter_zeta = 0.0;
    // Loop over the nodes on the edge
    for (unsigned n = 0; n < n_p; n++)
    {
      // Get the node number
      unsigned node_number = start + multiplier * n;
      // Now get the intrinsic coordinate
      node_pt(node_number)->get_coordinates_on_boundary(boundary, zeta);
      // Now multiply by the shape function
      inter_zeta += zeta[0] * psi(node_number);
    }

    // Set the value of the intrinsic coordinate
    zeta[0] = inter_zeta;
  }


  //===================================================================
  /// If a neighbouring element has already created a node at
  /// a position corresponding to the local fractional position within the
  /// present element, s_fraction, return
  /// a pointer to that node. If not, return NULL (0). If the node is
  /// periodic the flag is_periodic will be true
  //===================================================================
  Node* RefineableQElement<2>::node_created_by_neighbour(
    const Vector<double>& s_fraction, bool& is_periodic)
  {
    using namespace QuadTreeNames;

    // Calculate the edges on which the node lies
    Vector<int> edges;
    if (s_fraction[0] == 0.0)
    {
      edges.push_back(W);
    }
    if (s_fraction[0] == 1.0)
    {
      edges.push_back(E);
    }
    if (s_fraction[1] == 0.0)
    {
      edges.push_back(S);
    }
    if (s_fraction[1] == 1.0)
    {
      edges.push_back(N);
    }

    // Find the number of edges
    unsigned n_size = edges.size();
    // If there are no edges, then there is no neighbour, return 0
    if (n_size == 0)
    {
      return 0;
    }

    Vector<unsigned> translate_s(2);
    Vector<double> s_lo_neigh(2);
    Vector<double> s_hi_neigh(2);
    Vector<double> s(2);

    int neigh_edge, diff_level;
    bool in_neighbouring_tree;

    // Loop over the edges
    for (unsigned j = 0; j < n_size; j++)
    {
      // Find pointer to neighbouring element along edge
      QuadTree* neigh_pt;
      neigh_pt = quadtree_pt()->gteq_edge_neighbour(edges[j],
                                                    translate_s,
                                                    s_lo_neigh,
                                                    s_hi_neigh,
                                                    neigh_edge,
                                                    diff_level,
                                                    in_neighbouring_tree);

      // Neighbour exists
      if (neigh_pt != 0)
      {
        // Have its nodes been created yet?
        // if(neigh_pt->object_pt()->node_pt(0)!=0)
        if (neigh_pt->object_pt()->nodes_built())
        {
          // We now need to translate the nodal location
          // as defined in terms of the fractional coordinates of the present
          // element into those of its neighbour

          // Calculate the local coordinate in the neighbour
          // Note that we need to use the translation scheme in case
          // the local coordinates are swapped in the neighbour.
          for (unsigned i = 0; i < 2; i++)
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
            // Now work out whether it's a periodic boundary
            // only possible if we have moved into a neighbouring tree
            if (in_neighbouring_tree)
            {
              // Return whether the neighbour is periodic
              is_periodic =
                quadtree_pt()->root_pt()->is_neighbour_periodic(edges[j]);
            }
            // Return the pointer to the neighbouring node
            return neighbour_node_pt;
          }
        }
      }
    }
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
  ///   - Doc newly created nodes in "new_nodes.dat" stored in the directory
  ///     of the DocInfo object (only if it's open!)
  /// - Finally, excute the element-specific further_build()
  ///   (empty by default -- must be overloaded for specific elements).
  ///   This deals with any build operations that are not included
  ///   in the generic process outlined above. For instance, in
  ///   Crouzeix Raviart elements we need to initialise the internal
  ///   pressure values in manner consistent with the pressure
  ///   distribution in the father element.
  //==================================================================
  void RefineableQElement<2>::build(Mesh*& mesh_pt,
                                    Vector<Node*>& new_node_pt,
                                    bool& was_already_built,
                                    std::ofstream& new_nodes_file)
  {
    using namespace QuadTreeNames;

    // Get the number of 1d nodes
    unsigned n_p = nnode_1d();

    // Check whether static father_bound needs to be created
    if (Father_bound[n_p].nrow() == 0)
    {
      setup_father_bounds();
    }

    // Pointer to my father (in quadtree impersonation)
    QuadTree* father_pt = dynamic_cast<QuadTree*>(quadtree_pt()->father_pt());

    // What type of son am I? Ask my quadtree representation...
    int son_type = Tree_pt->son_type();

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
      RefineableQElement<2>* father_el_pt =
        dynamic_cast<RefineableQElement<2>*>(father_pt->object_pt());

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

      Vector<double> s_lo(2);
      Vector<double> s_hi(2);
      Vector<double> s(2);
      Vector<double> x(2);

      // Setup vertex coordinates in father element:
      //--------------------------------------------
      switch (son_type)
      {
        case SW:
          s_lo[0] = -1.0;
          s_hi[0] = 0.0;
          s_lo[1] = -1.0;
          s_hi[1] = 0.0;
          break;

        case SE:
          s_lo[0] = 0.0;
          s_hi[0] = 1.0;
          s_lo[1] = -1.0;
          s_hi[1] = 0.0;
          break;

        case NE:
          s_lo[0] = 0.0;
          s_hi[0] = 1.0;
          s_lo[1] = 0.0;
          s_hi[1] = 1.0;
          break;

        case NW:
          s_lo[0] = -1.0;
          s_hi[0] = 0.0;
          s_lo[1] = 0.0;
          s_hi[1] = 1.0;
          break;
      }

      // Pass macro element pointer on to sons and
      // set coordinates in macro element
      // hierher why can I see this?
      if (father_el_pt->Macro_elem_pt != 0)
      {
        set_macro_elem_pt(father_el_pt->Macro_elem_pt);
        for (unsigned i = 0; i < 2; i++)
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
        Vector<double> x_small(2);
        Vector<double> x_large(2);

        Vector<double> s_fraction(2);
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

            // Local node number
            jnod = i0 + n_p * i1;

            // Check whether the father's node is periodic if so, complain
            /* {
             Node* father_node_pt = father_el_pt->node_pt(jnod);
             if((father_node_pt->is_a_copy()) ||
                (father_node_pt->position_is_a_copy()))
              {
               throw OomphLibError(
                "Can't handle periodic nodes (yet).",
                OOMPH_CURRENT_FUNCTION,
                OOMPH_EXCEPTION_LOCATION);
              }
              }*/

            // Initialise flag: So far, this node hasn't been built
            // or copied yet
            bool node_done = false;

            // Get the pointer to the node in the father, returns NULL
            // if there is not node
            Node* created_node_pt =
              father_el_pt->get_node_at_local_coordinate(s);

            // Does this node already exist in father element?
            //------------------------------------------------
            if (created_node_pt != 0)
            {
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
              // Was the node created by one of its neighbours
              // Whether or not the node lies on an edge can be calculated
              // by from the fractional position
              bool is_periodic = false;
              ;
              created_node_pt =
                node_created_by_neighbour(s_fraction, is_periodic);

              // If the node was so created, assign the pointers
              if (created_node_pt != 0)
              {
                // If the node is periodic
                if (is_periodic)
                {
                  // Now the node must be on a boundary, but we don't know which
                  // one
                  // The returned created_node_pt is actually the neighbouring
                  // periodic node
                  Node* neighbour_node_pt = created_node_pt;

                  // Determine the edge on which the new node will live
                  int father_bound = Father_bound[n_p](jnod, son_type);

                  // Storage for the set of Mesh boundaries on which the
                  // appropriate father edge lives.
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
                  // Case where a new node lives on more than one boundary
                  // seems fishy enough to flag
                  if (boundaries.size() > 1)
                  {
                    throw OomphLibError(
                      "boundaries.size()!=1 seems a bit strange..\n",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
                  }

                  // Case when there are no boundaries, we are in big trouble
                  if (boundaries.size() == 0)
                  {
                    std::ostringstream error_stream;
                    error_stream << "Periodic node is not on a boundary...\n"
                                 << "Coordinates: " << created_node_pt->x(0)
                                 << " " << created_node_pt->x(1) << "\n";
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
                    // Get position from father element -- this uses the macro
                    // element representation if appropriate. If the node
                    // turns out to be a hanging node later on, then
                    // its position gets adjusted in line with its
                    // hanging node interpolation.
                    Vector<double> x_prev(2);
                    father_el_pt->get_x(t, s, x_prev);
                    // Set previous positions of the new node
                    for (unsigned i = 0; i < 2; i++)
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
                      Vector<double> zeta(1);
                      father_el_pt->interpolated_zeta_on_edge(
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
                // Node has been created
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
                ;
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
                    // Case where a new node lives on more than one boundary
                    // seems fishy enough to flag
                    if (boundaries.size() > 1)
                    {
                      throw OomphLibError(
                        "boundaries.size()!=1 seems a bit strange..\n",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
                    }

                    // Case when there are no boundaries, we are in big trouble
                    if (boundaries.size() == 0)
                    {
                      std::ostringstream error_stream;
                      error_stream << "Periodic node is not on a boundary...\n"
                                   << "Coordinates: " << created_node_pt->x(0)
                                   << " " << created_node_pt->x(1) << "\n";
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
                      // Get position from father element -- this uses the macro
                      // element representation if appropriate. If the node
                      // turns out to be a hanging node later on, then
                      // its position gets adjusted in line with its
                      // hanging node interpolation.
                      Vector<double> x_prev(2);
                      father_el_pt->get_x(t, s, x_prev);
                      // Set previous positions of the new node
                      for (unsigned i = 0; i < 2; i++)
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
                        Vector<double> zeta(1);
                        father_el_pt->interpolated_zeta_on_edge(
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
                  // Node has been created
                  node_done = true;
                } // Node does not exist in son of neighbouring element
              } // Node does not exist in neighbouring element
            } // Node does not exist in father element

            // Node has not been built anywhere ---> build it here
            if (!node_done)
            {
              // Firstly, we need to determine whether or not a node lies
              // on the boundary before building it, because
              // we actually assign a different type of node on boundaries.

              // The node can only be on a Mesh boundary if it
              // lives on an edge that is shared with an edge of its
              // father element; i.e. it is not created inside the father
              // element Determine the edge on which the new node will live
              int father_bound = Father_bound[n_p](jnod, son_type);

              // Storage for the set of Mesh boundaries on which the
              // appropriate father edge lives.
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
              // Case where a new node lives on more than one boundary
              // seems fishy enough to flag
              if (boundaries.size() > 1)
              {
                throw OomphLibError(
                  "boundaries.size()!=1 seems a bit strange..\n",
                  OOMPH_CURRENT_FUNCTION,
                  OOMPH_EXCEPTION_LOCATION);
              }
#endif

              // If the node lives on a mesh boundary,
              // then we need to create a boundary node
              if (boundaries.size() > 0)
              {
                // Create node and set the pointer to it from the element
                created_node_pt =
                  construct_boundary_node(jnod, time_stepper_pt);
                // Add to vector of new nodes
                new_node_pt.push_back(created_node_pt);

                // Now we need to work out whether to pin the values at
                // the new node based on the boundary conditions applied at
                // its Mesh boundary

                // Get the boundary conditions from the father
                Vector<int> bound_cons(ncont_interpolated_values());
                father_el_pt->get_bcs(father_bound, bound_cons);

                // Loop over the values and pin, if necessary
                unsigned n_value = created_node_pt->nvalue();
                for (unsigned k = 0; k < n_value; k++)
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
                  RefineableSolidQElement<2>* father_solid_el_pt =
                    dynamic_cast<RefineableSolidQElement<2>*>(father_el_pt);
#ifdef PARANOID
                  if (father_solid_el_pt == 0)
                  {
                    std::string error_message =
                      "We have a SolidNode outside a refineable SolidElement\n";
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
                    Vector<double> zeta(1);
                    father_el_pt->interpolated_zeta_on_edge(
                      *it, father_bound, s, zeta);

                    created_node_pt->set_coordinates_on_boundary(*it, zeta);
                  }
                }
              }
              // Otherwise the node is not on a Mesh boundary and
              // we create a normal "bulk" node
              else
              {
                // Create node and set the pointer to it from the element
                created_node_pt = construct_node(jnod, time_stepper_pt);
                // Add to vector of new nodes
                new_node_pt.push_back(created_node_pt);
              }

              // Now we set the position and values at the newly created node

              // In the first instance use macro element or FE representation
              // to create past and present nodal positions.
              // (THIS STEP SHOULD NOT BE SKIPPED FOR ALGEBRAIC
              // ELEMENTS AS NOT ALL OF THEM NECESSARILY IMPLEMENT
              // NONTRIVIAL NODE UPDATE FUNCTIONS. CALLING
              // THE NODE UPDATE FOR SUCH ELEMENTS/NODES WILL LEAVE
              // THEIR NODAL POSITIONS WHERE THEY WERE (THIS IS APPROPRIATE
              // ONCE THEY HAVE BEEN GIVEN POSITIONS) BUT WILL
              // NOT ASSIGN SENSIBLE INITIAL POSITONS!

              // Loop over # of history values
              for (unsigned t = 0; t < ntstorage; t++)
              {
                // Get position from father element -- this uses the macro
                // element representation if appropriate. If the node
                // turns out to be a hanging node later on, then
                // its position gets adjusted in line with its
                // hanging node interpolation.
                Vector<double> x_prev(2);
                father_el_pt->get_x(t, s, x_prev);

                // Set previous positions of the new node
                for (unsigned i = 0; i < 2; i++)
                {
                  created_node_pt->x(t, i) = x_prev[i];
                }
              }

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

            } // End of case when we build the node ourselves

            // Check if the element is an algebraic element
            AlgebraicElementBase* alg_el_pt =
              dynamic_cast<AlgebraicElementBase*>(this);

            // If the element is an algebraic element, setup
            // node position (past and present) from algebraic node update
            // function. This over-writes previous assingments that
            // were made based on the macro-element/FE representation.
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

            // If we have built the node and we are documenting our progress
            // write the (hopefully consistent position) to  the outputfile
            if ((!node_done) && (new_nodes_file.is_open()))
            {
              new_nodes_file << node_pt(jnod)->x(0) << " "
                             << node_pt(jnod)->x(1) << std::endl;
            }

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
        // Pass on non-halo proc id
        Non_halo_proc_ID =
          tree_pt()->father_pt()->object_pt()->non_halo_proc_ID();
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
    }
    // Element has already been built
    else
    {
      was_already_built = true;
    }
  }

  //====================================================================
  ///  Print corner nodes, use colour (default "BLACK")
  //====================================================================
  void RefineableQElement<2>::output_corners(std::ostream& outfile,
                                             const std::string& colour) const
  {
    Vector<double> s(2);
    Vector<double> corner(2);

    outfile << "ZONE I=2,J=2, C=" << colour << std::endl;

    s[0] = -1.0;
    s[1] = -1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << Number << std::endl;

    s[0] = 1.0;
    s[1] = -1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << Number << std::endl;

    s[0] = -1.0;
    s[1] = 1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << Number << std::endl;

    s[0] = 1.0;
    s[1] = 1.0;
    get_x(s, corner);
    outfile << corner[0] << " " << corner[1] << " " << Number << std::endl;


    outfile << "TEXT  CS = GRID, X = " << corner[0] << ",Y = " << corner[1]
            << ", HU = GRID, H = 0.01, AN = MIDCENTER, T =\"" << Number << "\""
            << std::endl;
  }

  //====================================================================
  /// Set up all hanging nodes. If we are documenting the output then
  /// open the output files and pass the open files to the helper function
  //====================================================================
  void RefineableQElement<2>::setup_hanging_nodes(
    Vector<std::ofstream*>& output_stream)
  {
#ifdef PARANOID
    if (output_stream.size() != 4)
    {
      throw OomphLibError("There must be four output streams",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    using namespace QuadTreeNames;

    // Setup hanging nodes on each edge of the element
    quad_hang_helper(-1, S, *(output_stream[0]));
    quad_hang_helper(-1, N, *(output_stream[1]));
    quad_hang_helper(-1, W, *(output_stream[2]));
    quad_hang_helper(-1, E, *(output_stream[3]));
  }

  //================================================================
  /// Internal function that sets up the hanging node scheme for
  /// a particular continuously interpolated value
  //===============================================================
  void RefineableQElement<2>::setup_hang_for_value(const int& value_id)
  {
    using namespace QuadTreeNames;

    std::ofstream dummy_hangfile;
    quad_hang_helper(value_id, S, dummy_hangfile);
    quad_hang_helper(value_id, N, dummy_hangfile);
    quad_hang_helper(value_id, W, dummy_hangfile);
    quad_hang_helper(value_id, E, dummy_hangfile);
  }


  //=================================================================
  /// Internal function to set up the hanging nodes on a particular
  /// edge of the element
  //=================================================================
  void RefineableQElement<2>::quad_hang_helper(const int& value_id,
                                               const int& my_edge,
                                               std::ofstream& output_hangfile)
  {
    using namespace QuadTreeNames;

    Vector<unsigned> translate_s(2);
    Vector<double> s_lo_neigh(2);
    Vector<double> s_hi_neigh(2);
    int neigh_edge, diff_level;
    bool in_neighbouring_tree;

    // Find pointer to neighbour in this direction
    QuadTree* neigh_pt;
    neigh_pt = quadtree_pt()->gteq_edge_neighbour(my_edge,
                                                  translate_s,
                                                  s_lo_neigh,
                                                  s_hi_neigh,
                                                  neigh_edge,
                                                  diff_level,
                                                  in_neighbouring_tree);

    // Neighbour exists and all nodes have been created
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
          is_periodic = tree_pt()->root_pt()->is_neighbour_periodic(my_edge);
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
          Vector<unsigned> translate_s_in_neigh(2);
          Vector<double> s_lo_neigh_of_neigh(2);
          Vector<double> s_hi_neigh_of_neigh(2);
          int neigh_edge_of_neigh, diff_level_of_neigh;
          bool in_neighbouring_tree_of_neigh;

          // Find pointer to neighbour of the neighbour on the edge
          // that we are currently considering
          QuadTree* neigh_of_neigh_pt;
          neigh_of_neigh_pt =
            neigh_pt->gteq_edge_neighbour(neigh_edge,
                                          translate_s_in_neigh,
                                          s_lo_neigh_of_neigh,
                                          s_hi_neigh_of_neigh,
                                          neigh_edge_of_neigh,
                                          diff_level_of_neigh,
                                          in_neighbouring_tree_of_neigh);

          // Set the value of the NEW neighbour and edge
          neigh_pt = neigh_of_neigh_pt;
          neigh_edge = neigh_edge_of_neigh;

          // Set the values of the translation schemes
          // Need to find the values of s_lo and s_hi
          // in the neighbour of the neighbour

          // Get the minimum and maximum values of the coordinate
          // in the neighbour (don't like this, but I think it's
          // necessary) Note that these values are hardcoded
          // in the quadtrees at some point!!
          double s_min = neigh_pt->object_pt()->s_min();
          double s_max = neigh_pt->object_pt()->s_max();
          Vector<double> s_lo_frac(2), s_hi_frac(2);
          // Work out the fractional position of the low and high points
          // of the original element
          for (unsigned i = 0; i < 2; i++)
          {
            s_lo_frac[i] = (s_lo_neigh[i] - s_min) / (s_max - s_min);
            s_hi_frac[i] = (s_hi_neigh[i] - s_min) / (s_max - s_min);
          }

          // We should now be able to construct the low and high points in
          // the neighbour of the neighbour
          for (unsigned i = 0; i < 2; i++)
          {
            s_lo_neigh[i] = s_lo_neigh_of_neigh[i] +
                            s_lo_frac[translate_s_in_neigh[i]] *
                              (s_hi_neigh_of_neigh[i] - s_lo_neigh_of_neigh[i]);
            s_hi_neigh[i] = s_lo_neigh_of_neigh[i] +
                            s_hi_frac[translate_s_in_neigh[i]] *
                              (s_hi_neigh_of_neigh[i] - s_lo_neigh_of_neigh[i]);
          }

          // Finally we must sort out the translation scheme
          Vector<unsigned> temp_translate(2);
          for (unsigned i = 0; i < 2; i++)
          {
            temp_translate[i] = translate_s_in_neigh[translate_s[i]];
          }
          for (unsigned i = 0; i < 2; i++)
          {
            translate_s[i] = temp_translate[i];
          }
        } // End of special treatment for periodic hanging nodes

        // Number of nodes in one dimension
        unsigned n_p = ninterpolating_node_1d(value_id);
        // Storage for the local nodes along the edge of the quadtree
        Node* local_node_pt = 0;
        // Loop over nodes along the edge
        for (unsigned i0 = 0; i0 < n_p; i0++)
        {
          // Storage for the fractional position of the node in the element
          Vector<double> s_fraction(2);

          // Find the local node and the fractional position of the node
          // which depends on the edge, of course
          switch (my_edge)
          {
            case N:
              s_fraction[0] =
                local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
              s_fraction[1] = 1.0;
              local_node_pt =
                interpolating_node_pt(i0 + n_p * (n_p - 1), value_id);
              break;

            case S:
              s_fraction[0] =
                local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
              s_fraction[1] = 0.0;
              local_node_pt = interpolating_node_pt(i0, value_id);
              break;

            case E:
              s_fraction[0] = 1.0;
              s_fraction[1] =
                local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
              local_node_pt =
                interpolating_node_pt(n_p - 1 + n_p * i0, value_id);
              break;

            case W:
              s_fraction[0] = 0.0;
              s_fraction[1] =
                local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
              local_node_pt = interpolating_node_pt(n_p * i0, value_id);
              break;

            default:
              throw OomphLibError("my_edge not N, S, W, E\n",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
          }

          // Calculate the local coordinates of the node in the neighbour
          Vector<double> s_in_neighb(2);
          for (unsigned i = 0; i < 2; i++)
          {
            s_in_neighb[i] = s_lo_neigh[i] + s_fraction[translate_s[i]] *
                                               (s_hi_neigh[i] - s_lo_neigh[i]);
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

            // If the node is not hanging geometrically, then we must make
            // it hang
            if (!local_node_pt->is_hanging())
            {
              make_hanging_node = true;
            }
            // Otherwise, it could be hanging geometrically, but still
            // require a different hanging scheme for this data value
            else
            {
              if (local_node_pt->hanging_pt(value_id) ==
                  local_node_pt->hanging_pt())
              {
                make_hanging_node = true;
              }
            }

            // If we do need to make the hanging node, then let's do it
            if (make_hanging_node == true)
            {
              // Cache refineable element used here
              RefineableElement* const obj_pt = neigh_pt->object_pt();

              // Get shape functions in neighbour element
              Shape psi(obj_pt->ninterpolating_node(value_id));
              obj_pt->interpolating_basis(s_in_neighb, psi, value_id);

              // Allocate the storage for the Hang pointer
              // which contains n_p nodes
              HangInfo* hang_pt = new HangInfo(n_p);

              // Loop over nodes on edge in neighbour and mark them as nodes
              // that this node depends on
              unsigned n_neighbour;

              // Number of nodes along edge in neighbour element
              for (unsigned n_edge = 0; n_edge < n_p; n_edge++)
              {
                switch (neigh_edge)
                {
                  case N:
                    n_neighbour = n_p * (n_p - 1) + n_edge;
                    break;

                  case S:
                    n_neighbour = n_edge;
                    break;

                  case W:
                    n_neighbour = n_p * n_edge;
                    break;

                  case E:
                    n_neighbour = n_p * n_edge + (n_p - 1);
                    break;

                  default:
                    throw OomphLibError("neigh_edge not N, S, W, E\n",
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                }

                // Push back neighbouring node and weight into
                // Vector of (pointers to)
                // master nodes and weights
                // The weight is merely the value of the shape function
                // corresponding to the node in the neighbour
                hang_pt->set_master_node_pt(
                  n_edge,
                  obj_pt->interpolating_node_pt(n_neighbour, value_id),
                  psi[n_neighbour]);
              }

              // Now set the hanging data for the position
              // This also constrains the data values associated with the
              // value id
              local_node_pt->set_hanging_pt(hang_pt, value_id);
            }

            // Dump the output if the file has been openeed
            if (output_hangfile.is_open())
            {
              output_hangfile << local_node_pt->x(0) << " "
                              << local_node_pt->x(1) << std::endl;
            }
          }
          // Otherwise check that the nodes are the same
          else
          {
#ifdef PARANOID
            if (local_node_pt != neighbouring_node_pt)
            {
              std::ostringstream warning_stream;
              warning_stream << "SANITY CHECK in quad_hang_helper      \n"
                             << "Current node      " << local_node_pt << " at "
                             << "(" << local_node_pt->x(0) << ", "
                             << local_node_pt->x(1) << ")"
                             << " is not hanging and has " << std::endl
                             << "Neighbour's node  " << neighbouring_node_pt
                             << " at "
                             << "(" << neighbouring_node_pt->x(0) << ", "
                             << neighbouring_node_pt->x(1) << ")" << std::endl
                             << "even though the two should be "
                             << "identical" << std::endl;
              OomphLibWarning(warning_stream.str(),
                              "RefineableQElement<2>::quad_hang_helper()",
                              OOMPH_EXCEPTION_LOCATION);
            }
#endif
          }

          // If we are doing the position, then
          if (value_id == -1)
          {
            // Get the nodal position from neighbour element
            Vector<double> x_in_neighb(2);
            neigh_pt->object_pt()->interpolated_x(s_in_neighb, x_in_neighb);

            // Fine adjust the coordinates (macro map will pick up boundary
            // accurately but will lead to different element edges)
            local_node_pt->x(0) = x_in_neighb[0];
            local_node_pt->x(1) = x_in_neighb[1];
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
  // template<unsigned NNODE_1D>
  void RefineableQElement<2>::check_integrity(double& max_error)
  {
    using namespace QuadTreeNames;

    // Number of nodes along edge
    unsigned n_p = nnode_1d();

    // Number of timesteps (incl. present) for which continuity is
    // to be checked.
    unsigned n_time = 1;

    // Initialise errors
    max_error = 0.0;
    Vector<double> max_error_x(2, 0.0);
    double max_error_val = 0.0;

    Vector<int> edges(4);
    edges[0] = S;
    edges[1] = N;
    edges[2] = W;
    edges[3] = E;

    // Loop over the edges
    for (unsigned edge_counter = 0; edge_counter < 4; edge_counter++)
    {
      Vector<unsigned> translate_s(2);
      Vector<double> s(2), s_lo_neigh(2), s_hi_neigh(2), s_fraction(2);
      int neigh_edge, diff_level;
      bool in_neighbouring_tree;

      // Find pointer to neighbour in this direction
      QuadTree* neigh_pt;
      neigh_pt = quadtree_pt()->gteq_edge_neighbour(edges[edge_counter],
                                                    translate_s,
                                                    s_lo_neigh,
                                                    s_hi_neigh,
                                                    neigh_edge,
                                                    diff_level,
                                                    in_neighbouring_tree);

      // Neighbour exists and has existing nodes
      if ((neigh_pt != 0) && (neigh_pt->object_pt()->nodes_built()))
      {
        // Need to exclude periodic nodes from this check
        // There are only periodic nodes if we are in a neighbouring tree
        bool is_periodic = false;
        if (in_neighbouring_tree)
        {
          // Is it periodic
          is_periodic =
            tree_pt()->root_pt()->is_neighbour_periodic(edges[edge_counter]);
        }

        // Loop over nodes along the edge
        for (unsigned i0 = 0; i0 < n_p; i0++)
        {
          // Storage for pointer to the local node
          Node* local_node_pt = 0;

          switch (edge_counter)
          {
            case 0:
              // Local fraction of node
              s_fraction[0] = local_one_d_fraction_of_node(i0, 0);
              s_fraction[1] = 0.0;
              // Get pointer to local node
              local_node_pt = node_pt(i0);
              break;

            case 1:
              // Local fraction of node
              s_fraction[0] = local_one_d_fraction_of_node(i0, 0);
              s_fraction[1] = 1.0;
              // Get pointer to local node
              local_node_pt = node_pt(i0 + n_p * (n_p - 1));
              break;

            case 2:
              // Local fraction of node
              s_fraction[0] = 0.0;
              s_fraction[1] = local_one_d_fraction_of_node(i0, 1);
              // Get pointer to local node
              local_node_pt = node_pt(n_p * i0);
              break;

            case 3:
              // Local fraction of node
              s_fraction[0] = 1.0;
              s_fraction[1] = local_one_d_fraction_of_node(i0, 1);
              // Get pointer to local node
              local_node_pt = node_pt(n_p - 1 + n_p * i0);
              break;
          }

          // Calculate the local coordinate and the local coordinate as viewed
          // from the neighbour
          Vector<double> s_in_neighb(2);
          for (unsigned i = 0; i < 2; i++)
          {
            // Local coordinate in this element
            s[i] = -1.0 + 2.0 * s_fraction[i];
            // Local coordinate in the neighbour
            s_in_neighb[i] = s_lo_neigh[i] + s_fraction[translate_s[i]] *
                                               (s_hi_neigh[i] - s_lo_neigh[i]);
          }

          // Loop over timesteps
          for (unsigned t = 0; t < n_time; t++)
          {
            // Get the nodal position from neighbour element
            Vector<double> x_in_neighb(2);
            neigh_pt->object_pt()->interpolated_x(t, s_in_neighb, x_in_neighb);

            // Check error only if the node is NOT periodic
            if (is_periodic == false)
            {
              for (int i = 0; i < 2; i++)
              {
                // Find the spatial error
                double err = std::fabs(local_node_pt->x(t, i) - x_in_neighb[i]);

                // If it's bigger than our tolerance, say so
                if (err > 1e-9)
                {
                  oomph_info << "errx " << err << " " << t << " "
                             << local_node_pt->x(t, i) << " " << x_in_neighb[i]
                             << std::endl;

                  oomph_info << "at " << local_node_pt->x(0) << " "
                             << local_node_pt->x(1) << std::endl;
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

            // Now figure out how many continuously interpolated values there
            // are
            unsigned num_val =
              neigh_pt->object_pt()->ncont_interpolated_values();

            // Check error
            for (unsigned ival = 0; ival < num_val; ival++)
            {
              double err = std::fabs(values[ival] - values_in_neighb[ival]);

              if (err > 1.0e-10)
              {
                oomph_info << local_node_pt->x(0) << " " << local_node_pt->x(1)
                           << " \n# "
                           << "erru (S)" << err << " " << ival << " "
                           << get_node_number(local_node_pt) << " "
                           << values[ival] << " " << values_in_neighb[ival]
                           << std::endl;
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

    max_error = max_error_x[0];
    if (max_error_x[1] > max_error) max_error = max_error_x[1];
    if (max_error_val > max_error) max_error = max_error_val;

    if (max_error > 1e-9)
    {
      oomph_info << "\n#------------------------------------ \n#Max error ";
      oomph_info << max_error_x[0] << " " << max_error_x[1] << " "
                 << max_error_val << std::endl;
      oomph_info << "#------------------------------------ \n " << std::endl;
    }
  }


  //========================================================================
  /// Static matrix for coincidence between son nodal points and
  /// father boundaries
  ///
  //========================================================================
  std::map<unsigned, DenseMatrix<int>> RefineableQElement<2>::Father_bound;


  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


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
  void RefineableSolidQElement<2>::get_solid_bcs(
    int bound, Vector<int>& solid_bound_cons) const
  {
    using namespace QuadTreeNames;

    // Spatial dimension of all nodes
    unsigned n_dim = this->nodal_dimension();

    // Set up temporary vectors to hold edge boundary conditions
    Vector<int> bound_cons1(n_dim), bound_cons2(n_dim);

    // Which boundary are we on?
    switch (bound)
    {
        // If on edge, just get the bcs
      case N:
      case S:
      case W:
      case E:
        get_edge_solid_bcs(bound, solid_bound_cons);
        break;

        // Most restrictive boundary at SE corner
      case SE:
        get_edge_solid_bcs(S, bound_cons1);
        get_edge_solid_bcs(E, bound_cons2);

        for (unsigned k = 0; k < n_dim; k++)
        {
          solid_bound_cons[k] = (bound_cons1[k] || bound_cons2[k]);
        }
        break;

        // Most restrictive boundary at SW corner
      case SW:
        get_edge_solid_bcs(S, bound_cons1);
        get_edge_solid_bcs(W, bound_cons2);

        for (unsigned k = 0; k < n_dim; k++)
        {
          solid_bound_cons[k] = (bound_cons1[k] || bound_cons2[k]);
        }
        break;

        // Most restrictive boundary at NW corner
      case NW:
        get_edge_solid_bcs(N, bound_cons1);
        get_edge_solid_bcs(W, bound_cons2);

        for (unsigned k = 0; k < n_dim; k++)
        {
          solid_bound_cons[k] = (bound_cons1[k] || bound_cons2[k]);
        }
        break;

        // Most restrictive boundary at NE corner
      case NE:
        get_edge_solid_bcs(N, bound_cons1);
        get_edge_solid_bcs(E, bound_cons2);

        for (unsigned k = 0; k < n_dim; k++)
        {
          solid_bound_cons[k] = (bound_cons1[k] || bound_cons2[k]);
        }
        break;

      default:
        throw OomphLibError(
          "Wrong boundary", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
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
  void RefineableSolidQElement<2>::get_edge_solid_bcs(
    const int& edge, Vector<int>& solid_bound_cons) const
  {
    using namespace QuadTreeNames;

    // Number of nodes along 1D edge
    unsigned n_p = nnode_1d();

    // Left- and Right-hand nodes
    unsigned left_node, right_node;

    // Set the left (lower) and right (upper) nodes for the edge
    switch (edge)
    {
      case N:
        left_node = n_p * (n_p - 1);
        right_node = n_p * n_p - 1;
        break;

      case S:
        left_node = 0;
        right_node = n_p - 1;
        break;

      case W:
        left_node = 0;
        right_node = n_p * (n_p - 1);
        break;

      case E:
        left_node = n_p - 1;
        right_node = n_p * n_p - 1;
        break;

      default:
        std::ostringstream error_stream;
        error_stream << "Wrong edge " << edge
                     << " passed to get_solid_edge_bcs(..)" << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    // Cast to SolidNodes
    SolidNode* left_node_pt = dynamic_cast<SolidNode*>(node_pt(left_node));
    SolidNode* right_node_pt = dynamic_cast<SolidNode*>(node_pt(right_node));
#ifdef PARANOID
    if (left_node_pt == 0)
    {
      throw OomphLibError(
        "Left node cannot be cast to SolidNode --> something is wrong",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
    if (right_node_pt == 0)
    {
      throw OomphLibError(
        "Right node cannot be cast to SolidNode --> something is wrong",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Number of coordinate directions
    unsigned n_dim = this->nodal_dimension();

    // Loop over all directions, the least restrictive value is
    // the boundary condition at the left multiplied by that at the right
    // Assuming that free is always zero and pinned is one
    for (unsigned k = 0; k < n_dim; k++)
    {
      solid_bound_cons[k] = left_node_pt->position_is_pinned(k) *
                            right_node_pt->position_is_pinned(k);
    }
  }

  /// Build required templates
  // template class RefineableQElement<2>;

} // namespace oomph
