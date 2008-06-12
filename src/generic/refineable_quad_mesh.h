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

 /// Set up QuadTreeForest. Wipes any existing tree structure and
 /// regards the currently active elements as the root trees in the forest.
 void setup_quadtree_forest()
  {
   // This wipes all elements/quadtrees in the tree representation
   // but leaves the leaf elements alone.
   if (this->Forest_pt!=0) delete this->Forest_pt;
      
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

};

}

#endif
