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
// Non-inline and non-templated functions for BinaryTree and BinaryTreeForest
// classes

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include <set>

// oomph-lib headers
#include "binary_tree.h"
#include "refineable_line_element.h"

namespace oomph
{
  //========================================================================
  /// Boolean indicating that static member data has been setup.
  //========================================================================
  bool BinaryTree::Static_data_has_been_setup = false;

  //========================================================================
  /// Colours for neighbours in various directions (static data).
  //========================================================================
  Vector<std::string> BinaryTree::Colour;

  //========================================================================
  /// S_base(direction): Initial value for coordinate s on the edge
  /// indicated by direction (L/R) (static data).
  //========================================================================
  Vector<double> BinaryTree::S_base;

  //========================================================================
  /// Translate (enumerated) directions into strings (static data).
  //========================================================================
  Vector<std::string> BinaryTree::Direct_string;

  //========================================================================
  /// Get opposite edge, e.g. Reflect_edge[N]=S (static data)
  //========================================================================
  Vector<int> BinaryTree::Reflect_edge;

  //========================================================================
  /// Array of direction/segment adjacency scheme:
  /// Is_adjacent(i_vertex,j_segment): Is vertex adjacent to segment?
  /// (static data)
  //========================================================================
  DenseMatrix<bool> BinaryTree::Is_adjacent;

  //========================================================================
  /// Reflection scheme: Reflect(direction,segment): Get mirror of segment
  /// in specified direction. E.g. Reflect(L,L)=R (static data)
  //========================================================================
  DenseMatrix<int> BinaryTree::Reflect;

  //========================================================================
  /// Set up the static data stored in the BinaryTree -- this needs to be
  /// called before BinaryTrees can be used. Automatically called by
  /// RefineableLineMesh constructor.
  //========================================================================
  void BinaryTree::setup_static_data()
  {
    using namespace BinaryTreeNames;

#ifdef PARANOID
    if (Tree::OMEGA != BinaryTree::OMEGA)
    {
      std::ostringstream error_stream;
      error_stream << "Inconsistent enumeration!  \n    Tree::OMEGA="
                   << Tree::OMEGA << "\nBinaryTree::OMEGA=" << BinaryTree::OMEGA
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Set flag to indicate that static data has been setup
    Static_data_has_been_setup = true;

    // Tecplot colours for neighbours in various directions
    Colour.resize(27);

    Colour[L] = "RED";
    Colour[R] = "CYAN";
    Colour[OMEGA] = "YELLOW";

    // S_base(direction): Initial value for coordinate s on the
    // edge indicated by direction (L/R)
    S_base.resize(27);

    S_base[L] = -1.0;
    S_base[R] = 1.0;

    // Translate (enumerated) directions into strings
    Direct_string.resize(27);

    Direct_string[L] = "L";
    Direct_string[R] = "R";
    Direct_string[OMEGA] = "OMEGA";

    // Build direction/segment adjacency scheme:
    // Is_adjacent(i_vertex,j_segment): Is vertex adjacent to segment?
    Is_adjacent.resize(27, 27);

    Is_adjacent(L, L) = true;
    Is_adjacent(R, L) = false;
    Is_adjacent(L, R) = false;
    Is_adjacent(R, R) = true;

    // Build reflection scheme: Reflect(direction,segment):
    // Get mirror of segment in direction
    Reflect.resize(27, 27);

    Reflect(L, L) = R;
    Reflect(R, L) = R;
    Reflect(L, R) = L;
    Reflect(R, R) = L;

    // Get opposite edge, e.g. Reflect_edge(L)=R
    Reflect_edge.resize(27);

    Reflect_edge[L] = R;
    Reflect_edge[R] = L;
  }

  //========================================================================
  /// Return pointer to greater or equal-sized edge neighbour in specified
  /// \c direction; also provide info regarding the relative size of the
  /// neighbour:
  /// - In the present binary tree, the left vertex is located at the local
  ///   coordinate s = -1. This point is located at the local coordinate
  ///   s = \c s_in_neighbour[0] in the neighbouring binary tree.
  /// - We're looking for a neighbour in the specified \c direction. When
  ///   viewed from the neighbouring binary tree, the edge that separates
  ///   the present binary tree from its neighbour is the neighbour's
  ///   \c edge edge. Since in 1D there can be no rotation between the two
  ///   binary trees, this is a simple reflection. For instance, if we're
  ///   looking for a neighhbour in the \c L [eft] \c direction, \c edge
  ///   will be \c R [ight].
  /// - \c diff_level <= 0 indicates the difference in refinement levels
  ///   between the two neighbours. If \c diff_level==0, the neighbour has
  ///   the same size as the current binary tree.
  /// - \c in_neighbouring_tree returns true if we have had to flip to a
  ///   different root, even if that root is actually the same (as it can
  ///   be in periodic problems).
  //========================================================================
  BinaryTree* BinaryTree::gteq_edge_neighbour(const int& direction,
                                              Vector<double>& s_in_neighbour,
                                              int& edge,
                                              int& diff_level,
                                              bool& in_neighbouring_tree) const
  {
    using namespace BinaryTreeNames;

#ifdef PARANOID
    if ((direction != L) && (direction != R))
    {
      std::ostringstream error_stream;
      error_stream << "Direction " << direction << " is not L or R"
                   << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Initialise in_neighbouring tree to false. It will be set to true
    // during the recursion if we do actually hop over into the neighbour
    in_neighbouring_tree = false;

    // Maximum level to which we're allowed to descend (we only want
    // greater-or-equal-sized neighbours)
    int max_level = Level;

    // Current element has the following root:
    BinaryTreeRoot* orig_root_pt = dynamic_cast<BinaryTreeRoot*>(Root_pt);

    // Initialise offset in local coordinate
    double s_diff = 0;

    // Initialise difference in level
    diff_level = 0;

    // Find neighbour
    BinaryTree* return_pt = gteq_edge_neighbour(direction,
                                                s_diff,
                                                diff_level,
                                                in_neighbouring_tree,
                                                max_level,
                                                orig_root_pt);

    BinaryTree* neighb_pt = return_pt;

    // If neighbour exists...
    if (neighb_pt != 0)
    {
      // What's the direction of the interfacial edge when viewed from within
      // the neighbour element?
      edge = Reflect_edge[direction];

      // What's the local coordinate s at which this edge is located in the
      // neighbouring element?
      s_in_neighbour[0] = S_base[edge];
    }
    return return_pt;
  }

  //========================================================================
  /// Find `greater-or-equal-sized edge neighbour' in given direction (L/R).
  ///
  /// This is an auxiliary routine which allows neighbour finding in
  /// adjacent binary trees. Needs to keep track of previous son types
  /// and the maximum level to which search is performed.
  ///
  /// Parameters:
  /// - direction (L/R): Direction in which neighbour has to be found.
  /// - s_diff: Offset of left vertex from corresponding vertex in
  ///   neighbour. Note that this is input/output as it needs to be
  ///   incremented/decremented during the recursive calls to this function.
  /// - edge: We're looking for the neighbour across our edge 'direction'
  ///   (L/R). When viewed from the neighbour, this edge is `edge' (L/R).
  ///   Since there is no relative rotation between neighbours this is a
  ///   mere reflection, e.g. direction=L  --> edge=R etc.
  /// - diff_level <= 0 indicates the difference in binary tree levels
  ///   between the current element and its neighbour.
  /// - max_level is the maximum level to which the neighbour search is
  ///   allowed to proceed. This is again necessary because in a forest,
  ///   the neighbour search isn't based on pure recursion.
  /// - orig_root_pt identifies the root node of the element whose
  ///   neighbour we're really trying to find by all these recursive calls.
  //========================================================================
  BinaryTree* BinaryTree::gteq_edge_neighbour(
    const int& direction,
    double& s_diff,
    int& diff_level,
    bool& in_neighbouring_tree,
    int max_level,
    BinaryTreeRoot* const& orig_root_pt) const
  {
    using namespace BinaryTreeNames;

#ifdef PARANOID
    if ((direction != L) && (direction != R))
    {
      std::ostringstream error_stream;
      error_stream << "Direction " << direction << " is not L or R"
                   << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    BinaryTree* next_el_pt;
    BinaryTree* return_el_pt;

    // STEP 1: Find the neighbour's father
    // -------

    // Does the element have a father?
    if (Father_pt != 0)
    {
      // If the present segment (whose location inside its father element
      // is specified by Son_type) is adjacent to the father's edge in the
      // required direction, then its neighbour has a different father
      // ---> we need to climb up the tree to the father and find his
      // neighbour in the required direction. Note that this is the cunning
      // recursive part. The returning may not stop until we hit the very
      // top of the tree, when the element does NOT have a father.
      if (Is_adjacent(direction, Son_type))
      {
        next_el_pt = dynamic_cast<BinaryTree*>(Father_pt)->gteq_edge_neighbour(
          direction,
          s_diff,
          diff_level,
          in_neighbouring_tree,
          max_level,
          orig_root_pt);
      }
      // If the present segment is not adjacent to the father's edge in the
      // required direction, then the neighbour has the same father and is
      // obtained by the appropriate reflection inside the father element.
      // This will only be called if we have not left the original tree.
      else
      {
        next_el_pt = dynamic_cast<BinaryTree*>(Father_pt);
      }

      // We're about to ascend one level
      diff_level -= 1;

      // Work out position of left corner of present edge in its father element
      s_diff += pow(0.5, -diff_level);

      // STEP 2: We have now located the neighbour's father and need to
      // ------- find the appropriate son.

      // Buffer cases where the neighbour (and hence its father) lie outside
      // the boundary
      if (next_el_pt != 0)
      {
        // If the father is a leaf then we can't descend to the same
        // level as the present node ---> simply return the father himself
        // as the (greater) neighbour. Same applies if we are about to
        // descend lower than the max_level (in a  neighbouring tree).
        if ((next_el_pt->Son_pt.size() == 0) ||
            (next_el_pt->Level > max_level - 1))
        {
          return_el_pt = next_el_pt;
        }
        // We have located the neighbour's father: The position of the
        // neighbour is obtained by `reflecting' the position of the
        // node itself. We know exactly how to reflect because we know which
        // son type we are and we have the pointer to the neighbour's father.
        else
        {
          int son_segment = Reflect(direction, Son_type);

          // The next element in the tree is the appropriate son of the
          // neighbour's father
          return_el_pt =
            dynamic_cast<BinaryTree*>(next_el_pt->Son_pt[son_segment]);

          // Work out position of left corner of present edge in next
          // higher element
          s_diff -= pow(0.5, -diff_level);

          // We have just descended one level
          diff_level += 1;
        }
      }
      // The neighbour's father lies outside the boundary --> the neighbour
      // itself does too --> return NULL
      else
      {
        return_el_pt = 0;
      }
    }
    // Element does not have a father --> check if it has a neighbouring
    // tree in the appropriate direction
    else
    {
      // Find neighbouring root
      if (Root_pt->neighbour_pt(direction) != 0)
      {
        // In this case we have moved to a neighbour, so set the flag
        in_neighbouring_tree = true;
        return_el_pt =
          dynamic_cast<BinaryTreeRoot*>(Root_pt->neighbour_pt(direction));
      }
      // No neighbouring tree, so there really is no neighbour --> return NULL
      else
      {
        return_el_pt = 0;
      }
    }

    return return_el_pt;
  }

  //========================================================================
  /// Self-test: Check neighbour finding routine. For each element in the
  /// tree and for each vertex, determine the distance between the vertex
  /// and its position in the neighbour. If the difference is less than
  /// Tree::Max_neighbour_finding_tolerance return success (0), otherwise
  /// failure (1).
  //========================================================================
  unsigned BinaryTree::self_test()
  {
    // Stick pointers to all nodes into Vector and number elements
    // in the process
    Vector<Tree*> all_nodes_pt;
    stick_all_tree_nodes_into_vector(all_nodes_pt);
    long int count = 0;
    unsigned long num_nodes = all_nodes_pt.size();
    for (unsigned long i = 0; i < num_nodes; i++)
    {
      all_nodes_pt[i]->object_pt()->set_number(++count);
    }

    // Check neighbours (distance between hanging nodes) -- don't print
    // (keep output streams closed)
    double max_error = 0.0;
    std::ofstream neighbours_file;
    std::ofstream neighbours_txt_file;
    BinaryTree::doc_neighbours(
      all_nodes_pt, neighbours_file, neighbours_txt_file, max_error);

    if (max_error > Max_neighbour_finding_tolerance)
    {
      oomph_info << "\n \n Failed self_test() for BinaryTree: Max. error "
                 << max_error << std::endl
                 << std::endl;
      return 1;
    }
    else
    {
      oomph_info << "\n \n Passed self_test() for BinaryTree: Max. error "
                 << max_error << std::endl
                 << std::endl;
      return 0;
    }
  }


  //========================================================================
  /// Constructor for BinaryTreeForest:
  ///
  /// Pass:
  ///  - trees_pt[], the Vector of pointers to the constituent trees
  ///    (BinaryTreeRoot objects)
  //========================================================================
  BinaryTreeForest::BinaryTreeForest(Vector<TreeRoot*>& trees_pt)
    : TreeForest(trees_pt)
  {
#ifdef LEAK_CHECK
    LeakCheckNames::BinaryTreeForest_build += 1;
#endif

    using namespace BinaryTreeNames;

    // Set up the neighbours
    find_neighbours();
  }

  //========================================================================
  /// Set up the neighbour lookup schemes for all constituent binary trees.
  //========================================================================
  void BinaryTreeForest::find_neighbours()
  {
    using namespace BinaryTreeNames;

    unsigned numtrees = ntree();
    unsigned n = 0; // to store nnode1d
    if (numtrees > 0)
    {
      n = Trees_pt[0]->object_pt()->nnode_1d();
    }
    else
    {
      throw OomphLibError(
        "Trying to setup the neighbour scheme for an empty forest\n",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    // Number of vertex nodes: 2
    unsigned n_vertex_node = 2;

    // Find connected trees by identifying those whose associated elements
    // share a common vertex node
    std::map<Node*, std::set<unsigned>> tree_assoc_with_vertex_node;

    // Loop over all trees
    for (unsigned i = 0; i < numtrees; i++)
    {
      // Loop over the vertex nodes of the associated element
      for (unsigned j = 0; j < n_vertex_node; j++)
      {
        Node* nod_pt = dynamic_cast<LineElementBase*>(Trees_pt[i]->object_pt())
                         ->vertex_node_pt(j);
        tree_assoc_with_vertex_node[nod_pt].insert(i);
      }
    }

    // For each tree we store a set of neighbouring trees
    // i.e. trees that share a node
    Vector<std::set<unsigned>> neighb_tree(numtrees);

    // Loop over vertex nodes
    for (std::map<Node*, std::set<unsigned>>::iterator it =
           tree_assoc_with_vertex_node.begin();
         it != tree_assoc_with_vertex_node.end();
         it++)
    {
      // Loop over connected elements twice
      for (std::set<unsigned>::iterator it_el1 = it->second.begin();
           it_el1 != it->second.end();
           it_el1++)
      {
        unsigned i = (*it_el1);
        for (std::set<unsigned>::iterator it_el2 = it->second.begin();
             it_el2 != it->second.end();
             it_el2++)
        {
          unsigned j = (*it_el2);
          // These two elements are connected
          if (i != j)
          {
            neighb_tree[i].insert(j);
          }
        }
      }
    }

    // Loop over all trees
    for (unsigned i = 0; i < numtrees; i++)
    {
      // Loop over their neighbours
      for (std::set<unsigned>::iterator it = neighb_tree[i].begin();
           it != neighb_tree[i].end();
           it++)
      {
        unsigned j = (*it);

        // is it the left-hand neighbour?
        bool is_L_neighbour = Trees_pt[j]->object_pt()->get_node_number(
                                Trees_pt[i]->object_pt()->node_pt(0)) != -1;

        // is it the right-hand neighbour?
        bool is_R_neighbour = Trees_pt[j]->object_pt()->get_node_number(
                                Trees_pt[i]->object_pt()->node_pt(n - 1)) != -1;

        if (is_L_neighbour) Trees_pt[i]->neighbour_pt(L) = Trees_pt[j];
        if (is_R_neighbour) Trees_pt[i]->neighbour_pt(R) = Trees_pt[j];
      }
    } // End of loop over all trees
  }

  //========================================================================
  /// Document and check all the neighbours in all the nodes in the forest.
  //========================================================================
  void BinaryTreeForest::check_all_neighbours(DocInfo& doc_info)
  {
    // Create Vector of elements
    Vector<Tree*> all_tree_nodes_pt;
    this->stick_all_tree_nodes_into_vector(all_tree_nodes_pt);

    // Create storage for information files
    std::ofstream neigh_file;
    std::ofstream neigh_txt_file;

    // If we are documenting the results, then open the files
    if (doc_info.is_doc_enabled())
    {
      std::ostringstream fullname;
      fullname << doc_info.directory() << "/neighbours" << doc_info.number()
               << ".dat";
      oomph_info << "opened " << fullname.str() << " to doc neighbours"
                 << std::endl;
      neigh_file.open(fullname.str().c_str());
      fullname.str("");
      fullname << doc_info.directory() << "/neighbours" << doc_info.number()
               << ".txt";
      oomph_info << "opened " << fullname.str() << " to doc neighbours"
                 << std::endl;
      neigh_txt_file.open(fullname.str().c_str());
    }

    // Call the standard documentation function
    double max_error = 0.0;
    BinaryTree::doc_neighbours(
      all_tree_nodes_pt, neigh_file, neigh_txt_file, max_error);

    // If the error is too large, complain
    if (max_error > Tree::max_neighbour_finding_tolerance())
    {
      std::ostringstream error_stream;
      error_stream << "Max. error in binary tree neighbour finding: "
                   << max_error << " is too big" << std::endl;
      error_stream
        << "i.e. bigger than Tree::max_neighbour_finding_tolerance()="
        << Tree::max_neighbour_finding_tolerance() << std::endl;

      // Close the files if they were opened
      if (doc_info.is_doc_enabled())
      {
        neigh_file.close();
        neigh_txt_file.close();
      }

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      oomph_info << "Max. error in binary tree neighbour finding: " << max_error
                 << " is OK" << std::endl;
      oomph_info
        << "i.e. less than BinaryTree::max_neighbour_finding_tolerance()="
        << BinaryTree::max_neighbour_finding_tolerance() << std::endl;
    }

    // Close the files if they were opened
    if (doc_info.is_doc_enabled())
    {
      neigh_file.close();
      neigh_txt_file.close();
    }
  }

  //========================================================================
  /// Self test: Check neighbour finding routine. For each element in the
  /// tree and for each vertex, determine the distance between the vertex
  /// and its position in the neighbour. If the difference is less than
  /// Tree::Max_neighbour_finding_tolerance return success (0), otherwise
  /// failure (1).
  //========================================================================
  unsigned BinaryTreeForest::self_test()
  {
    // Stick pointers to all nodes into Vector and number elements
    // in the process
    Vector<Tree*> all_forest_nodes_pt;
    stick_all_tree_nodes_into_vector(all_forest_nodes_pt);
    long int count = 0;
    unsigned long num_nodes = all_forest_nodes_pt.size();
    for (unsigned long i = 0; i < num_nodes; i++)
    {
      all_forest_nodes_pt[i]->object_pt()->set_number(++count);
    }

    // Check neighbours (distance between hanging nodes) -- don't print
    // (keep output streams closed)
    double max_error = 0.0;
    std::ofstream neighbours_file;
    std::ofstream neighbours_txt_file;
    BinaryTree::doc_neighbours(
      all_forest_nodes_pt, neighbours_file, neighbours_txt_file, max_error);
    if (max_error > BinaryTree::max_neighbour_finding_tolerance())
    {
      oomph_info << "\n \n Failed self_test() for BinaryTree: Max. error "
                 << max_error << std::endl
                 << std::endl;
      return 1;
    }
    else
    {
      oomph_info << "\n \n Passed self_test() for BinaryTree: Max. error "
                 << max_error << std::endl
                 << std::endl;
      return 0;
    }
  }

  //========================================================================
  /// Doc/check all neighbours of binary tree ("nodes") contained in the
  /// Vector forest_node_pt. Output into neighbours_file which can be viewed
  /// from tecplot with BinaryTreeNeighbours.mcr. Neighbour info and errors
  /// are displayed on neighbours_txt_file. Finally, compute the maximum
  /// error between vertices when viewed from neighbouring element. Output
  /// is suppressed if the output streams are closed.
  //========================================================================
  void BinaryTree::doc_neighbours(Vector<Tree*> forest_nodes_pt,
                                  std::ofstream& neighbours_file,
                                  std::ofstream& neighbours_txt_file,
                                  double& max_error)
  {
    using namespace BinaryTreeNames;

    int diff_level;
    double s_diff;
    bool in_neighbouring_tree;
    int edge = OMEGA;

    Vector<double> s(1);
    Vector<double> x(1);
    Vector<int> prev_son_type;

    Vector<double> s_in_neighbour(1);

    Vector<double> x_small(1);
    Vector<double> x_large(1);

    // Initialise error in vertex positions
    max_error = 0.0;

    // Loop over all elements to assign numbers for plotting
    // -----------------------------------------------------
    unsigned long num_nodes = forest_nodes_pt.size();
    for (unsigned long i = 0; i < num_nodes; i++)
    {
      // Set number
      forest_nodes_pt[i]->object_pt()->set_number(i);
    }

    // Loop over all elements for checks
    // ---------------------------------
    for (unsigned long i = 0; i < num_nodes; i++)
    {
      // Doc the element itself
      BinaryTree* el_pt = dynamic_cast<BinaryTree*>(forest_nodes_pt[i]);

      // If the object is incomplete, complain
      if (el_pt->object_pt()->nodes_built())
      {
        // Print it
        if (neighbours_file.is_open())
        {
          dynamic_cast<RefineableQElement<1>*>(el_pt->object_pt())
            ->output_corners(neighbours_file, "BLACK");
        }

        // Loop over directions to find neighbours
        // ---------------------------------------
        for (int direction = L; direction <= R; direction++)
        {
          // Initialise difference in levels and coordinate offset
          diff_level = 0;
          s_diff = 0.0;

          // Find greater-or-equal-sized neighbour...
          BinaryTree* neighb_pt = el_pt->gteq_edge_neighbour(
            direction, s_in_neighbour, edge, diff_level, in_neighbouring_tree);

          // If neighbour exist and nodes are created: Doc it
          if ((neighb_pt != 0) && (neighb_pt->object_pt()->nodes_built()))
          {
            // Doc neighbour stats
            if (neighbours_txt_file.is_open())
            {
              neighbours_txt_file
                << Direct_string[direction] << " neighbour of "
                << el_pt->object_pt()->number() << " is "
                << neighb_pt->object_pt()->number() << " diff_level "
                << diff_level << " s_diff " << s_diff
                << " inside neighbour the edge is " << Direct_string[edge]
                << std::endl
                << std::endl;
            }

            // Plot neighbour in the appropriate colour
            if (neighbours_file.is_open())
            {
              dynamic_cast<RefineableQElement<1>*>(neighb_pt->object_pt())
                ->output_corners(neighbours_file, Colour[direction]);
            }

            // Check that local coordinates in the larger element (obtained
            // via s_diff) lead to the same spatial point as the node vertices
            // in the current element
            {
              if (neighbours_file.is_open())
              {
                neighbours_file << "ZONE I=1 \n";
              }

              // Left vertex:
              // ------------

              // Get coordinates in large (neighbour) element
              s[0] = s_in_neighbour[0];
              neighb_pt->object_pt()->get_x(s, x_large);

              // Get coordinates in small element
              Vector<double> s(1);
              s[0] = S_base[direction];
              el_pt->object_pt()->get_x(s, x_small);

              // Need to exclude periodic nodes from this check. There can
              // only be periodic nodes if we have moved into the neighbour
              bool is_periodic = false;
              if (in_neighbouring_tree)
              {
                // Is the node periodic?
                is_periodic =
                  el_pt->root_pt()->is_neighbour_periodic(direction);
              }

              double error = 0.0;
              // Only bother to calculate the error if the node is NOT periodic
              if (is_periodic == false)
              {
                error += pow(x_small[0] - x_large[0], 2);
              }

              // Take the root of the square error
              error = sqrt(error);
              if (neighbours_txt_file.is_open())
              {
                neighbours_txt_file << "Error (1) " << error << std::endl;
              }

              if (std::fabs(error) > max_error)
              {
                max_error = std::fabs(error);
              }

              if (neighbours_file.is_open())
              {
                neighbours_file << x_large[0] << "  0 \n";
              }

              // Right vertex:
              // -------------

              // Get coordinates in large (neighbour) element
              s[0] = s_in_neighbour[0];
              neighb_pt->object_pt()->get_x(s, x_large);

              // Get coordinates in small element
              s[0] = S_base[direction];
              el_pt->object_pt()->get_x(s, x_small);

              error = 0.0;
              // Only do this if we are NOT periodic
              if (is_periodic == false)
              {
                error += pow(x_small[0] - x_large[0], 2);
              }
              // Take the root of the square error
              error = sqrt(error);

              // error =
              // sqrt(pow(x_small[0]-x_large[0],2)+pow(x_small[1]-x_large[1],2));
              if (neighbours_txt_file.is_open())
              {
                neighbours_txt_file << "Error (2) " << error << std::endl;
              }

              if (std::fabs(error) > max_error)
              {
                max_error = std::fabs(error);
              }

              if (neighbours_file.is_open())
              {
                neighbours_file << x_large[0] << "  0 \n";
              }
            }
            //        else
            //         {
            //          // No neighbour: write dummy zone so tecplot can find
            //          four
            //          // neighbours for every element
            //          if (neighbours_file.is_open())
            //            {
            //             neighbours_file << "ZONE I=1 \n";
            //             neighbours_file << "-0.05 -0.05   0 \n";
            //             neighbours_file << "-0.05 -0.05   0 \n";
            //            }
            //         }
          }
          // If neighbour does not exist: Insert blank zones into file
          // so that tecplot can find four neighbours for every element
          else
          {
            if (neighbours_file.is_open())
            {
              neighbours_file << "ZONE \n 0.00  0 \n";
              neighbours_file << "ZONE I=1 \n";
              neighbours_file << "-0.05  0 \n";
              neighbours_file << "-0.05  0 \n";
            }
          }
        }
      } // End of case when element can be documented
    }
  }

} // namespace oomph
