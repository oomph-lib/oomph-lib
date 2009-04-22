//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#ifndef QMESH2OOMPH_D_HEADER
#define QMESH2OOMPH_D_HEADER

#include "quad_mesh.h"
#include "refineable_mesh.h"
#include "refineable_quad_element.h"

namespace oomph
{

//=======================================================================
/// Intermediate mesh class that implements the mesh adaptation functions
/// specified in the RefineableMesh class for meshes that contain the
/// refineable variant of QElement s [The class ELEMENT provided
/// as the template parameter must be of type 
/// RefineableQElement<2>].
/// 
/// Mesh adaptation/refinement is implemented by QuadTree 
/// procedures and any concrete implementation of this class needs to
/// provide a QuadTreeForest representation of the initial (coarse) mesh.
//=======================================================================
template <class ELEMENT>
class RefineableQuadMesh : public virtual RefineableMesh<ELEMENT>, 
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
 RefineableQuadMesh(const RefineableQuadMesh& dummy) 
  { 
   BrokenCopy::broken_copy("RefineableQuadMesh");
  } 
 
 /// Broken assignment operator
 void operator=(const RefineableQuadMesh&) 
  {
   BrokenCopy::broken_assign("RefineableQuadMesh");
  }

 /// Destructor:
 virtual ~RefineableQuadMesh() {}


 /// \short Set up the tree forest associated with the Mesh. 
 /// Forwards call to setup_quadtree_forest()
 virtual void setup_tree_forest()
  {
   setup_quadtree_forest();
  }

 /// \short Set up QuadTreeForest. Wipes any existing tree structure below
 /// the minimum refinement level and regards the elements at that level
 /// as the root trees in the forest.
 void setup_quadtree_forest()
  {
   if (this->Forest_pt!=0)
    {
     // Get all the tree nodes
     Vector<Tree*> all_tree_nodes_pt;
     this->Forest_pt->stick_all_tree_nodes_into_vector(all_tree_nodes_pt);

     // Get min and max refinement level from the tree
     unsigned min_ref;
     unsigned max_ref;
     this->get_refinement_levels(min_ref,max_ref);

     // Vector to store trees for new Forest
     Vector<TreeRoot*> trees_pt;

     // Loop over tree nodes (e.g. elements)
     unsigned n_tree_nodes=all_tree_nodes_pt.size();
     for (unsigned e=0;e<n_tree_nodes;e++)
      {
       Tree* tree_pt=all_tree_nodes_pt[e];

       // If the object_pt has been flushed then we don't want to keep
       // this tree
       if (tree_pt->object_pt()!=0)
        {
         // Get the refinement level of the current tree node
         RefineableElement* el_pt=dynamic_cast<RefineableElement*>
          (tree_pt->object_pt());
         unsigned level=el_pt->refinement_level();

         // If we are below the minimum refinement level, remove tree
         if (level<min_ref)
          {
           // Flush sons for this tree
           tree_pt->flush_sons();

           // Delete the tree (no recursion)
           delete tree_pt;

           // Delete the element
           delete el_pt;
          }
         else if (level==min_ref)
          {
           // Get the sons (if there are any) and store them
           unsigned n_sons=tree_pt->nsons();
           Vector<Tree*> backed_up_sons(n_sons);
           for (unsigned i_son=0;i_son<n_sons;i_son++)
            {
             backed_up_sons[i_son]=tree_pt->son_pt(i_son);
            }

           // Make the element into a new treeroot
           QuadTreeRoot* tree_root_pt=new QuadTreeRoot(el_pt);

           // Loop over sons and make the new treeroot their father
           for (unsigned i_son=0;i_son<n_sons;i_son++)
            {
             Tree* son_pt=backed_up_sons[i_son];
             son_pt->set_father_pt(tree_root_pt);
            }

           // Add treeroot to the trees_pt vector
           trees_pt.push_back(tree_root_pt);
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
     unsigned n_element=this->nelement();
     for (unsigned e=0;e<n_element;e++)
      {
       // Get pointer to full element type 
       ELEMENT* el_pt=dynamic_cast<ELEMENT*>(this->element_pt(e));
     
       // Build associated quadtree(root) -- pass pointer to corresponding
       // finite element and add the pointer to vector of quadtree (roots):
       trees_pt.push_back(new QuadTreeRoot(el_pt));
      } 
   
     // Plant QuadTreeRoots in QuadTreeForest
     this->Forest_pt = new QuadTreeForest(trees_pt);
    }
  }

};

}

#endif
