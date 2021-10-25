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
#ifndef OOMPH_GENERIC_REFINEABLE_QUAD_MESH_HEADER
#define OOMPH_GENERIC_REFINEABLE_QUAD_MESH_HEADER

#include <limits.h>

#include "quad_mesh.h"
#include "refineable_mesh.h"
#include "refineable_quad_element.h"
// Include to fill in additional_setup_shared_node_scheme() function
#include "refineable_mesh.template.cc"

namespace oomph
{
  //=======================================================================
  /// Intermediate mesh class that implements the mesh adaptation functions
  /// specified in the TreeBasedRefineableMesh class for meshes that contain the
  /// refineable variant of QElement s [The class ELEMENT provided
  /// as the template parameter must be of type
  /// RefineableQElement<2>].
  ///
  /// Mesh adaptation/refinement is implemented by QuadTree
  /// procedures and any concrete implementation of this class needs to
  /// provide a QuadTreeForest representation of the initial (coarse) mesh.
  //=======================================================================
  template<class ELEMENT>
  class RefineableQuadMesh : public virtual TreeBasedRefineableMesh<ELEMENT>,
                             public virtual QuadMeshBase
  {
  public:
    /// Constructor: Setup static quadtree data
    RefineableQuadMesh()
    {
      // QuadTree static data needs to be setup before quadtree-based mesh
      // refinement works
      QuadTree::setup_static_data();
    }

    /// Broken copy constructor
    RefineableQuadMesh(const RefineableQuadMesh& dummy) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const RefineableQuadMesh&) = delete;*/

    /// Destructor:
    virtual ~RefineableQuadMesh() {}


    /// Set up the tree forest associated with the Mesh.
    /// Forwards call to setup_quadtree_forest()
    virtual void setup_tree_forest()
    {
      setup_quadtree_forest();
    }

    /// Set up QuadTreeForest. Wipes any existing tree structure below
    /// the minimum refinement level and regards the elements at that level
    /// as the root trees in the forest.
    void setup_quadtree_forest()
    {
      // A forst pointer is setup at least once when the mesh is initially
      // created in serial and stays around
      if (this->Forest_pt != 0)
      {
        // Get all the tree nodes
        Vector<Tree*> all_tree_nodes_pt;
        this->Forest_pt->stick_all_tree_nodes_into_vector(all_tree_nodes_pt);

        // Get min and max refinement level from the tree
        unsigned local_min_ref = 0;
        unsigned local_max_ref = 0;
        this->get_refinement_levels(local_min_ref, local_max_ref);

#ifdef OOMPH_HAS_MPI

        // Reconcile between processors: If (e.g. following
        // distribution/pruning) the mesh has no elements on this processor)
        // then ignore its contribution to the poll of max/min refinement levels
        int int_local_min_ref = local_min_ref;

        if (this->nelement() == 0)
        {
          int_local_min_ref = INT_MAX;
        }

        int int_min_ref = 0;
        MPI_Allreduce(&int_local_min_ref,
                      &int_min_ref,
                      1,
                      MPI_INT,
                      MPI_MIN,
                      Comm_pt->mpi_comm());

        unsigned min_ref = unsigned(int_min_ref);

#else

        unsigned min_ref = local_min_ref;

#endif

        // If we have no elements there's nothing more to be done --
        // we only came in here to participate in the communication
        if (this->nelement() == 0)
        {
          // Flush the Forest's current trees
          this->Forest_pt->flush_trees();

          // Delete the old Forest
          delete this->Forest_pt;

          // Empty dummy vector to build empty forest
          Vector<TreeRoot*> trees_pt;

          // Make a new (empty) Forest
          this->Forest_pt = new QuadTreeForest(trees_pt);
          return;
        }

        // Vector to store trees for new Forest
        Vector<TreeRoot*> trees_pt;

        // Loop over tree nodes (e.g. elements)
        unsigned n_tree_nodes = all_tree_nodes_pt.size();
        for (unsigned e = 0; e < n_tree_nodes; e++)
        {
          Tree* tree_pt = all_tree_nodes_pt[e];

          // If the object_pt has been flushed then we don't want to keep
          // this tree
          if (tree_pt->object_pt() != 0)
          {
            // Get the refinement level of the current tree node
            RefineableElement* el_pt =
              dynamic_cast<RefineableElement*>(tree_pt->object_pt());
            unsigned level = el_pt->refinement_level();

            // If we are below the minimum refinement level, remove tree
            if (level < min_ref)
            {
              // Flush sons for this tree
              tree_pt->flush_sons();

              // Delete the tree (no recursion)
              delete tree_pt;

              // Delete the element
              delete el_pt;
            }
            else if (level == min_ref)
            {
              // Get the sons (if there are any) and store them
              unsigned n_sons = tree_pt->nsons();
              Vector<Tree*> backed_up_sons(n_sons);
              for (unsigned i_son = 0; i_son < n_sons; i_son++)
              {
                backed_up_sons[i_son] = tree_pt->son_pt(i_son);
              }

              // Make the element into a new tree-root
              QuadTreeRoot* tree_root_pt = new QuadTreeRoot(el_pt);

              // Pass sons
              tree_root_pt->set_son_pt(backed_up_sons);

              // Loop over sons and make the new treeroot their father
              for (unsigned i_son = 0; i_son < n_sons; i_son++)
              {
                Tree* son_pt = backed_up_sons[i_son];

                // Tell the son about its new father (which is also the root)
                son_pt->set_father_pt(tree_root_pt);
                son_pt->root_pt() = tree_root_pt;

                // ...and then tell all the descendants too
                Vector<Tree*> all_sons_pt;
                son_pt->stick_all_tree_nodes_into_vector(all_sons_pt);
                unsigned n = all_sons_pt.size();
                for (unsigned i = 0; i < n; i++)
                {
                  all_sons_pt[i]->root_pt() = tree_root_pt;
                }
              }

              // Add tree-root to the trees_pt vector
              trees_pt.push_back(tree_root_pt);

              // Now kill the original (non-root) tree: First
              // flush sons for this tree
              tree_pt->flush_sons();

              // ...then delete the tree (no recursion)
              delete tree_pt;
            }
          }
          else // tree_pt->object_pt() is null, so delete tree
          {
            // Flush sons for this tree
            tree_pt->flush_sons();

            // Delete the tree (no recursion)
            delete tree_pt;
          }
        }

        // Flush the Forest's current trees
        this->Forest_pt->flush_trees();

        // Delete the old Forest
        delete this->Forest_pt;

        // Make a new Forest with the trees_pt roots created earlier
        this->Forest_pt = new QuadTreeForest(trees_pt);
      }
      else // Create a new Forest from scratch in the "usual" uniform way
      {
        // Each finite element in the coarse base mesh gets associated
        // with (the root of) a QuadTree. Store QuadTreeRoots in vector:
        Vector<TreeRoot*> trees_pt;

        // Loop over all elements, build corresponding QuadTree
        // and store QuadTreeRoots in vector:
        unsigned n_element = this->nelement();
        for (unsigned e = 0; e < n_element; e++)
        {
          // Get pointer to full element type
          ELEMENT* el_pt = dynamic_cast<ELEMENT*>(this->element_pt(e));

          // Build associated quadtree(root) -- pass pointer to corresponding
          // finite element and add the pointer to vector of quadtree (roots):
          trees_pt.push_back(new QuadTreeRoot(el_pt));
        }

        // Plant QuadTreeRoots in QuadTreeForest
        this->Forest_pt = new QuadTreeForest(trees_pt);
      }
    }
  };

} // namespace oomph

#endif
