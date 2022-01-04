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
// Header file for binary tree and binary tree forest classes
#ifndef OOMPH_BINARY_TREE_HEADER
#define OOMPH_BINARY_TREE_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// oomph-lib headers
#include "tree.h"
#include "matrices.h"

namespace oomph
{
  //======================================================================
  /// Namespace for BinaryTree directions
  //======================================================================
  namespace BinaryTreeNames
  {
    /// Directions (L/R). OMEGA is used if a direction is undefined
    /// in a certain context
    enum
    {
      L,
      R,
      OMEGA = 26
    };
  }; // namespace BinaryTreeNames

  // Forward class definition for class representing the root of a BinaryTree
  class BinaryTreeRoot;

  //======================================================================
  /// BinaryTree class: Recursively defined, generalised binary tree.
  ///
  /// A BinaryTree has:
  /// - a pointer to the object (of type RefineableQElement<1>) that it
  ///   represents in a mesh refinement context.
  /// - a Vector of pointers to its two (L/R) sons (which are
  ///   themselves binary trees). If the Vector of pointers to the sons
  ///   has zero length, the BinaryTree is a "leaf node" in the overall
  ///   binary tree.
  /// - a pointer to its father. If this pointer is NULL, the BinaryTree
  ///   is the root node of the overall binary tree.
  /// This data is stored in the Tree base class.
  ///
  /// The tree can also be part of a forest. If that is the case, the root
  /// will have pointers to the roots of neighbouring binary trees.
  ///
  /// The objects contained in the binary tree are assumed to be
  /// line elements whose geometry is parametrised by local coordinates
  /// \f$ {\bf s} \in [-1,1] \f$.
  ///
  /// The tree can be traversed and actions performed at all its
  /// "nodes" or only at the leaf "nodes" ("nodes" without sons).
  ///
  /// Finally, the leaf "nodes" can be split depending on
  /// a criteria defined by the object.
  ///
  /// Note that BinaryTrees are only generated by splitting existing
  /// BinaryTrees. Therefore, the constructors are protected. The
  /// only BinaryTree that "Joe User" can create is the (derived) class
  /// BinaryTreeRoot.
  //======================================================================
  class BinaryTree : public virtual Tree
  {
  public:
    /// Destructor. Note: Deleting a binary tree also deletes the
    /// objects associated with all non-leaf nodes!
    virtual ~BinaryTree() {}

    /// Broken copy constructor
    BinaryTree(const BinaryTree& dummy) = delete;

    /// Broken assignment operator
    void operator=(const BinaryTree&) = delete;

    /// Overload the function construct_son to ensure that the son
    /// is a specific BinaryTree and not a general Tree.
    Tree* construct_son(RefineableElement* const& object_pt,
                        Tree* const& father_pt,
                        const int& son_type)
    {
      BinaryTree* temp_binary_pt =
        new BinaryTree(object_pt, father_pt, son_type);
      return temp_binary_pt;
    }

    /// Return pointer to greater or equal-sized edge neighbour
    /// in specified \c direction; also provide info regarding the relative
    /// size of the neighbour:
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
    ///   between the two neighbours. If \c diff_level==0, the neighbour
    ///   has the same size as the current binary tree.
    /// - \c in_neighbouring_tree indicates whether the neighbour is actually
    ///   in another tree in the forest. The introduction of this flag
    ///   was necessitated by periodic problems where a TreeRoot can be its
    ///   own neighbour.
    BinaryTree* gteq_edge_neighbour(const int& direction,
                                    Vector<double>& s_in_neighbour,
                                    int& edge,
                                    int& diff_level,
                                    bool& in_neighbouring_tree) const;

    /// Self-test: Check all neighbours. Return success (0)
    /// if the maximum distance between corresponding points in the
    /// neighbours is less than the tolerance specified in the
    /// static value BinaryTree::Max_neighbour_finding_tolerance.
    unsigned self_test();

    /// Set up the static data, reflection schemes, etc.
    static void setup_static_data();

    /// Doc/check all neighbours of binary tree (nodes) contained in the
    /// Vector forest_node_pt. Output into neighbours_file which can be viewed
    /// from tecplot with BinaryTreeNeighbours.mcr. Neighbour info and errors
    /// are displayed on neighbours_txt_file. Finally, compute the maximum
    /// error between vertices when viewed from the neighbouring element.
    /// If the two filestreams are closed, output is suppressed.
    static void doc_neighbours(Vector<Tree*> forest_nodes_pt,
                               std::ofstream& neighbours_file,
                               std::ofstream& neighbours_txt_file,
                               double& max_error);

    /// Translate (enumerated) directions into strings
    static Vector<std::string> Direct_string;

  protected:
    /// Default constructor (empty and broken)
    BinaryTree()
    {
      throw OomphLibError(
        "Don't call an empty constructor for a BinaryTree object",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Default constructor for empty (root) tree: no father, no sons;
    /// just pass a pointer to its object. Protected because BinaryTrees can
    /// only be created internally, during the split operation. Only
    /// BinaryTreeRoots can be created externally.
    BinaryTree(RefineableElement* const& object_pt) : Tree(object_pt) {}

    /// Constructor for tree that has a father: Pass it the pointer
    /// to its object, the pointer to its father and tell it what type of son
    /// (L/R) it is. Protected because BinaryTrees can only be created
    /// internally, during the split operation. Only BinaryTreeRoots can be
    /// created externally.
    BinaryTree(RefineableElement* const& object_pt,
               Tree* const& father_pt,
               const int& son_type)
      : Tree(object_pt, father_pt, son_type)
    {
    }

    /// Boolean indicating that static member data has been setup
    static bool Static_data_has_been_setup;

  private:
    /// Find greater or equal-sized edge neighbour in direction.
    /// Auxiliary internal routine which passes additional information around.
    BinaryTree* gteq_edge_neighbour(const int& direction,
                                    double& s_diff,
                                    int& diff_level,
                                    bool& in_neighbouring_tree,
                                    int max_level,
                                    BinaryTreeRoot* const& orig_root_pt) const;

    /// Colours for neighbours in various directions
    static Vector<std::string> Colour;

    /// S_base(direction): Initial value for coordinate s on the edge
    /// indicated by direction (L/R)
    static Vector<double> S_base;

    /// Get opposite edge, e.g. Reflect_edge[L]=R
    static Vector<int> Reflect_edge;

    /// Array of direction/segment adjacency scheme:
    /// Is_adjacent(i_vertex,j_segment): Is vertex adjacent to segment?
    static DenseMatrix<bool> Is_adjacent;

    /// Reflection scheme: Reflect(direction,segment): Get mirror
    /// of segment in specified direction. E.g. Reflect(L,L)=R.
    static DenseMatrix<int> Reflect;
  };


  //======================================================================
  /// BinaryTreeRoot is a BinaryTree that forms the root of a (recursive)
  /// binary tree. The "root node" is special as it holds additional
  /// information about its neighbours.
  //======================================================================
  class BinaryTreeRoot : public virtual BinaryTree, public virtual TreeRoot
  {
  public:
    /// Constructor for the (empty) root binary tree: Pass pointer to
    /// associated object, a RefineableQElement<1>.
    BinaryTreeRoot(RefineableElement* const& object_pt)
      : Tree(object_pt), BinaryTree(object_pt), TreeRoot(object_pt)
    {
#ifdef PARANOID
      // Check that static member data has been setup
      if (!Static_data_has_been_setup)
      {
        std::string error_message =
          "Static member data hasn't been setup yet.\n";
        error_message +=
          "Call BinaryTree::setup_static_data() before creating\n";
        error_message += "any BinaryTreeRoots\n";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    /// Broken copy constructor
    BinaryTreeRoot(const BinaryTreeRoot& dummy) = delete;

    /// Broken assignment operator
    void operator=(const BinaryTreeRoot&) = delete;

    /// If binary_tree_root_pt is a neighbour, return the direction
    /// (L/R) in which it is found, otherwise return OMEGA
    int direction_of_neighbour(BinaryTreeRoot* binary_tree_root_pt)
    {
      using namespace BinaryTreeNames;

      if (Neighbour_pt[L] == binary_tree_root_pt)
      {
        return L;
      }
      if (Neighbour_pt[R] == binary_tree_root_pt)
      {
        return R;
      }

      // If we get here, it's not a neighbour
      return OMEGA;
    }
  };


  //======================================================================
  /// A BinaryTreeForest consists of a collection of BinaryTreeRoots.
  /// Each member tree can have neighbours to its left and right.
  //======================================================================
  class BinaryTreeForest : public TreeForest
  {
  public:
    /// Default constructor (empty and broken)
    BinaryTreeForest()
    {
      // Throw an error
      throw OomphLibError(
        "Don't call an empty constructor for a BinaryTreeForest object",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Constructor: Pass vector of pointers to the roots of the
    /// constituent BinaryTrees
    BinaryTreeForest(Vector<TreeRoot*>& trees_pt);

    /// Broken copy constructor
    BinaryTreeForest(const BinaryTreeForest& dummy) = delete;

    /// Broken assignment operator
    void operator=(const BinaryTreeForest&) = delete;

    /// Destructor: Delete the constituent binary trees (and thus
    /// the objects associated with its non-leaf nodes!)
    virtual ~BinaryTreeForest() {}

    /// Document and check all the neighbours of all the nodes in
    /// the forest. DocInfo object specifies the output directory and file
    /// numbers for the various files. If \c doc_info.disable_doc() has been
    /// called, no output is created.
    void check_all_neighbours(DocInfo& doc_info);

    /// A line mesh cannot have hanging nodes so make this function empty
    void open_hanging_node_files(DocInfo& doc_info,
                                 Vector<std::ofstream*>& output_stream)
    {
    }

    /// Self-test: Check all neighbours. Return success (0) if the
    /// maximum distance between corresponding points in the neighbours is
    /// less than the tolerance specified in the static value
    /// BinaryTree::Max_neighbour_finding_tolerance.
    unsigned self_test();

  private:
    /// Construct the neighbour lookup scheme
    void find_neighbours();

    /// Return pointer to i-th root binary tree in this forest (performs
    /// a dynamic cast from the TreeRoot to a BinaryTreeRoot).
    BinaryTreeRoot* binary_tree_pt(const unsigned& i)
    {
      return dynamic_cast<BinaryTreeRoot*>(Trees_pt[i]);
    }

    /// Given the number i of the root binary tree in this forest,
    /// return a pointer to its neighbour in the specified direction.
    /// NULL if neighbour doesn't exist. (This does the dynamic cast
    /// from a TreeRoot to a BinaryTreeRoot internally).
    BinaryTreeRoot* binary_neigh_pt(const unsigned& i, const int& direction)
    {
      return dynamic_cast<BinaryTreeRoot*>(
        Trees_pt[i]->neighbour_pt(direction));
    }
  };

} // namespace oomph

#endif
