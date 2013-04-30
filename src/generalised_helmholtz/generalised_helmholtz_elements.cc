//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
//Non-inline functions for gen Helmholtz elements
#include "generalised_helmholtz_elements.h"



namespace oomph
{



//===================================================================
/// \short Namespace with functions that allow the construction of
/// PML layers on axis aligned boundaries
//===================================================================
namespace PMLHelper
{
 /// helper function for sorting the right boundary nodes
 bool sorter_right_boundary(Node* nod_i_pt, Node* nod_j_pt)
 {
  return (nod_i_pt -> x(1) < nod_j_pt -> x(1));
 }
 
 /// helper function for sorting the top boundary nodes
 bool sorter_top_boundary(Node* nod_i_pt, Node* nod_j_pt)
 {  
  return (nod_i_pt -> x(0) < nod_j_pt -> x(0));
 }
 
 /// helper function for sorting the left boundary nodes
 bool sorter_left_boundary(Node* nod_i_pt, Node* nod_j_pt)
 {  
  return (nod_i_pt -> x(1) > nod_j_pt -> x(1));
 }
 
 /// helper function for sorting the bottom boundary nodes
 bool sorter_bottom_boundary(Node* nod_i_pt, Node* nod_j_pt)
 {
  return ( nod_i_pt -> x(0) > nod_j_pt -> x(0));
 }


 //==========================================================================
 /// "Constructor" for PML mesh,aligned with the right physical domain boundary
 //==========================================================================
 Mesh* create_right_pml_mesh(Mesh* bulk_mesh_pt,
                             const unsigned& right_boundary_id,
                             const unsigned& n_x_right_pml,
                             const double& width_x_right_pml)
 {  
  // Look at the right boundary of the triangular mesh
  unsigned n_right_boundary_node = 
   bulk_mesh_pt -> nboundary_node(right_boundary_id);
  
  // Create a vector of ordered boundary nodes
  Vector<Node*> ordered_right_boundary_node_pt(n_right_boundary_node);
  
  // Fill the vector with the nodes on the respective boundary
  for (unsigned n=0; n<n_right_boundary_node; n++)
   {
    ordered_right_boundary_node_pt[n] = 
     bulk_mesh_pt->boundary_node_pt(right_boundary_id,n);
   }
  
  // Sort them from lowest to highest (in y coordinate)
  std::sort(ordered_right_boundary_node_pt.begin(),
            ordered_right_boundary_node_pt.end(),
            sorter_right_boundary);
  
  // The number of elements in y is taken from the triangular mesh
  unsigned n_y_right_pml = 
   bulk_mesh_pt -> nboundary_element(right_boundary_id);
  
  // Specific PML sizes needed, taken directly from physical domain
  double l_pml_right_x_start = 
   ordered_right_boundary_node_pt[0] -> x(0);
  /// \short PML layer with added to the bulk mesh coordinate
  double l_pml_right_x_end   = 
   width_x_right_pml 
   + ordered_right_boundary_node_pt[0] -> x(0);
  double l_pml_right_y_start = 
   ordered_right_boundary_node_pt[0] -> x(1);
  double l_pml_right_y_end   = 
   ordered_right_boundary_node_pt[n_right_boundary_node-1] -> x(1);
  
  // Rectangular boundary id to be merged with triangular mesh
  unsigned right_quadPML_boundary_id = 3;
  
  /// \short Get the number of nodes along the (bulk) element edges
  /// from the first element
  unsigned nnode_1d=bulk_mesh_pt->finite_element_pt(0)->nnode_1d();

  // Create the mesh to be designated to the PML
  Mesh* pml_right_mesh_pt = 0;

  // Build the right one
  if (nnode_1d==3)
   {
    pml_right_mesh_pt=
     new PMLQuadMesh<QGeneralisedHelmholtzElement<2,3> >
     (bulk_mesh_pt, right_boundary_id, right_quadPML_boundary_id,  
      n_x_right_pml, n_y_right_pml, 
      l_pml_right_x_start, l_pml_right_x_end, 
      l_pml_right_y_start, l_pml_right_y_end);
   }
  else if (nnode_1d==4)
   {
    pml_right_mesh_pt=
     new PMLQuadMesh<QGeneralisedHelmholtzElement<2,4> >
     (bulk_mesh_pt, right_boundary_id, right_quadPML_boundary_id,  
      n_x_right_pml, n_y_right_pml, 
      l_pml_right_x_start, l_pml_right_x_end, 
      l_pml_right_y_start, l_pml_right_y_end);  
   }
  else
   {
    std::ostringstream error_stream;
    error_stream << "Currently I can't build PML meshes for   "
                 << std::endl
                 <<" bulk meshes with nnode_1d = " << nnode_1d
                 << std::endl
                 << "nodes along the element edge. Please fix this...  "
                 << std::endl;
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  // Enable PML damping on the entire mesh
  unsigned n_element_pml_right = pml_right_mesh_pt->nelement();
  for(unsigned e=0;e<n_element_pml_right;e++)
   {
    // Upcast
    GeneralisedHelmholtzBase<2>* el_pt = 
     dynamic_cast<GeneralisedHelmholtzBase<2>*>
     (pml_right_mesh_pt->element_pt(e));
    el_pt -> enable_pml(0, l_pml_right_x_start, l_pml_right_x_end);
   }
  
  /// \short Exterior boundary needs to be set to Dirichlet 0
  /// in both real and imaginary parts
  unsigned n_bound_pml_right = pml_right_mesh_pt->nboundary();
  for(unsigned b=0;b<n_bound_pml_right;b++)
   {
    unsigned n_node = pml_right_mesh_pt -> nboundary_node(b);
    for (unsigned n=0;n<n_node;n++)
     {
      if (b==1) {
       pml_right_mesh_pt -> boundary_node_pt(b,n)->pin(0);
       pml_right_mesh_pt -> boundary_node_pt(b,n)->pin(1); 
      }
     }
   }
  
  for(unsigned i=0;i<n_bound_pml_right;i++)
   {
    // How many nodes are there on this boundary?
    unsigned n_node = pml_right_mesh_pt->nboundary_node(i);
    
    // Loop over the nodes on boundary
    for (unsigned n=0;n<n_node;n++)
     {
      // Get pointer to node
      Node* nod_pt=pml_right_mesh_pt->boundary_node_pt(i,n);
      if (i==1)
       { 
        nod_pt->set_value(0,0.0);
        nod_pt->set_value(1,0.0);
       }
     }
   }
  
  /// \short Return the finalized mesh, with PML enabled 
  /// and boundary conditions added
  return pml_right_mesh_pt;
 }

 //==========================================================================
 /// "Constructor" for PML mesh, aligned with the top physical domain boundary
 //==========================================================================
 Mesh* create_top_pml_mesh(Mesh* bulk_mesh_pt,
                           const unsigned& top_boundary_id,
                           const unsigned& n_y_top_pml,
                           const double& width_y_top_pml)
 {  
  // Look at the top boundary of the triangular mesh
  unsigned n_top_boundary_node = 
   bulk_mesh_pt -> nboundary_node(top_boundary_id);
  
  // Create a vector of ordered boundary nodes
  Vector<Node*> ordered_top_boundary_node_pt(n_top_boundary_node);
  
  // Fill the vector with the nodes on the respective boundary
  for (unsigned n=0; n<n_top_boundary_node; n++)
   {
    ordered_top_boundary_node_pt[n] = 
     bulk_mesh_pt->boundary_node_pt(top_boundary_id,n);
   }
  
  // Sort them from lowest to highest (in y coordinate)
  std::sort(ordered_top_boundary_node_pt.begin(),
            ordered_top_boundary_node_pt.end(),
            sorter_top_boundary);
 
  // The number of elements in x is taken from the triangular mesh
  unsigned n_x_top_pml = bulk_mesh_pt -> nboundary_element(top_boundary_id);
  
  // Specific PML sizes needed, taken directly from physical domain
  double l_pml_top_x_start = 
   ordered_top_boundary_node_pt[0] -> x(0);
  double l_pml_top_x_end   = 
   ordered_top_boundary_node_pt[n_top_boundary_node-1] -> x(0);
  double l_pml_top_y_start = 
   ordered_top_boundary_node_pt[0] -> x(1);
  /// \short PML layer width added to the bulk mesh coordinate
  double l_pml_top_y_end   = 
   width_y_top_pml 
   + ordered_top_boundary_node_pt[0] -> x(1);
  
  unsigned top_quadPML_boundary_id = 0;

  /// \short Get the number of nodes along the (bulk) element edges
  /// from the first element
  unsigned nnode_1d=bulk_mesh_pt->finite_element_pt(0)->nnode_1d();

  // Create the mesh to be designated to the PML
  Mesh* pml_top_mesh_pt = 0;

  // Build the top PML mesh
  if (nnode_1d==3)
   {
    pml_top_mesh_pt=
     new PMLQuadMesh<QGeneralisedHelmholtzElement<2,3> >
     (bulk_mesh_pt, top_boundary_id, top_quadPML_boundary_id,  
      n_x_top_pml, n_y_top_pml, 
      l_pml_top_x_start, l_pml_top_x_end, 
      l_pml_top_y_start, l_pml_top_y_end);
   }
  else if (nnode_1d==4)
   {
    pml_top_mesh_pt=
     new PMLQuadMesh<QGeneralisedHelmholtzElement<2,4> >
     (bulk_mesh_pt, top_boundary_id, top_quadPML_boundary_id,  
      n_x_top_pml, n_y_top_pml, 
      l_pml_top_x_start, l_pml_top_x_end, 
      l_pml_top_y_start, l_pml_top_y_end); 
   }
  else {
   std::ostringstream error_stream;
   error_stream << "Currently I can't build PML meshes for   "
                << std::endl
                <<" bulk meshes with nnode_1d = " << nnode_1d
                << std::endl
                << "nodes along the element edge. Please fix this...  "
                << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
  
  // Enable PML damping on the entire mesh
  unsigned n_element_pml_top = pml_top_mesh_pt->nelement();
  for(unsigned e=0;e<n_element_pml_top;e++)
   {
    // Upcast
    GeneralisedHelmholtzBase<2>* el_pt = 
     dynamic_cast<GeneralisedHelmholtzBase<2>*>
     (pml_top_mesh_pt->element_pt(e));
    el_pt -> enable_pml(1, l_pml_top_y_start, l_pml_top_y_end);
   }
  
  /// \short Exterior boundary needs to be set to Dirichlet 0
  /// in both real and imaginary parts
  unsigned n_bound_pml_top = pml_top_mesh_pt->nboundary();
  for(unsigned b=0;b<n_bound_pml_top;b++)
   {
    unsigned n_node = pml_top_mesh_pt -> nboundary_node(b);
    for (unsigned n=0;n<n_node;n++)
     {
      if (b==2) {
       pml_top_mesh_pt -> boundary_node_pt(b,n)->pin(0);
       pml_top_mesh_pt -> boundary_node_pt(b,n)->pin(1); 
      }
     }
   }
  
  for(unsigned i=0;i<n_bound_pml_top;i++)
   {
    // How many nodes are there on this boundary?
    unsigned n_node = pml_top_mesh_pt->nboundary_node(i);
    
    // Loop over the nodes on boundary
    for (unsigned n=0;n<n_node;n++)
     {
      // Get pointer to node
      Node* nod_pt=pml_top_mesh_pt->boundary_node_pt(i,n);
      if (i==2)
       { 
        nod_pt->set_value(0, 0.0);
        nod_pt->set_value(1, 0.0);
       }
     }
   }

  /// \short Return the finalized mesh, with PML enabled 
  /// and boundary conditions added
  return pml_top_mesh_pt;
 }
 
 //==========================================================================
 /// "Constructor" for PML mesh, aligned with the left physical domain boundary
 //==========================================================================
 Mesh* create_left_pml_mesh(Mesh* bulk_mesh_pt,
                            const unsigned& left_boundary_id,
                            const unsigned& n_x_left_pml,
                            const double& width_x_left_pml)
 {
  // Look at the left boundary of the triangular mesh
  unsigned n_left_boundary_node = 
   bulk_mesh_pt -> nboundary_node(left_boundary_id);
  
  // Create a vector of ordered boundary nodes
  Vector<Node*> ordered_left_boundary_node_pt(n_left_boundary_node);
  
  // Fill the vector with the nodes on the respective boundary
  for (unsigned n=0; n<n_left_boundary_node; n++)
   {
    ordered_left_boundary_node_pt[n] = 
     bulk_mesh_pt->boundary_node_pt(left_boundary_id,n);
   }
  
  // Sort them from lowest to highest (in y coordinate)
  std::sort(ordered_left_boundary_node_pt.begin(),
            ordered_left_boundary_node_pt.end(),
            sorter_left_boundary);
  
  // The number of elements in y is taken from the triangular mesh
  unsigned n_y_left_pml = bulk_mesh_pt -> nboundary_element(left_boundary_id);
 
  // Specific PML sizes needed, taken directly from physical domain
  /// \short PML layer width subtracted from left bulk mesh coordinate
  double l_pml_left_x_start = 
   - width_x_left_pml 
   + ordered_left_boundary_node_pt[n_left_boundary_node-1] -> x(0);
  double l_pml_left_x_end   = 
   ordered_left_boundary_node_pt[n_left_boundary_node-1] -> x(0);
  double l_pml_left_y_start = 
   ordered_left_boundary_node_pt[n_left_boundary_node-1] -> x(1);
  double l_pml_left_y_end   = 
   ordered_left_boundary_node_pt[0] -> x(1);

  unsigned left_quadPML_boundary_id = 1;

  /// \short Get the number of nodes along the (bulk) element edges
  /// from the first element
  unsigned nnode_1d=bulk_mesh_pt->finite_element_pt(0)->nnode_1d();

  // Create the mesh to be designated to the PML
  Mesh* pml_left_mesh_pt = 0;

  // Build the left PML mesh
  if (nnode_1d==3)
   {
    pml_left_mesh_pt=
     new PMLQuadMesh<QGeneralisedHelmholtzElement<2,3> >
     (bulk_mesh_pt, left_boundary_id, left_quadPML_boundary_id,  
      n_x_left_pml, n_y_left_pml, 
      l_pml_left_x_start, l_pml_left_x_end, 
      l_pml_left_y_start, l_pml_left_y_end);
   }
  else if (nnode_1d==4)
   {
    pml_left_mesh_pt=
     new PMLQuadMesh<QGeneralisedHelmholtzElement<2,4> >
     (bulk_mesh_pt, left_boundary_id, left_quadPML_boundary_id,  
      n_x_left_pml, n_y_left_pml, 
      l_pml_left_x_start, l_pml_left_x_end, 
      l_pml_left_y_start, l_pml_left_y_end);  
   }
  else{
   std::ostringstream error_stream;
   error_stream << "Currently I can't build PML meshes for   "
                << std::endl
                <<" bulk meshes with nnode_1d = " << nnode_1d
                 << std::endl
                << "nodes along the element edge. Please fix this...  "
                << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
  
  // Enable PML damping on the entire mesh
  unsigned n_element_pml_left = pml_left_mesh_pt->nelement();
  for(unsigned e=0;e<n_element_pml_left;e++)
   {
    // Upcast
    GeneralisedHelmholtzBase<2>* el_pt = 
     dynamic_cast<GeneralisedHelmholtzBase<2>*>
     (pml_left_mesh_pt->element_pt(e));
    el_pt -> enable_pml(0, l_pml_left_x_end, l_pml_left_x_start);
   }

  /// \short Exterior boundary needs to be set to Dirichlet 0
  /// in both real and imaginary parts
  unsigned n_bound_pml_left = pml_left_mesh_pt->nboundary();
  for(unsigned b=0;b<n_bound_pml_left;b++)
   {
    unsigned n_node = pml_left_mesh_pt -> nboundary_node(b);
    for (unsigned n=0;n<n_node;n++)
     {
      if (b==3) {
       pml_left_mesh_pt -> boundary_node_pt(b,n)->pin(0);
       pml_left_mesh_pt -> boundary_node_pt(b,n)->pin(1); 
      }
     }
   }

  for(unsigned i=0;i<n_bound_pml_left;i++)
   {
    // How many nodes are there on this boundary?
    unsigned n_node = pml_left_mesh_pt->nboundary_node(i);
    
    // Loop over the nodes on boundary
    for (unsigned n=0;n<n_node;n++)
     {
      // Get pointer to node
      Node* nod_pt=pml_left_mesh_pt->boundary_node_pt(i,n);
      if (i==3)
       { 
        nod_pt->set_value(0,0.0);
        nod_pt->set_value(1,0.0);
       }
     }
   }

  /// \short Return the finalized mesh, with PML enabled 
  /// and boundary conditions added
  return pml_left_mesh_pt;
 }

 //==========================================================================
 ///"Constructor" for PML mesh,aligned with the bottom physical domain boundary
 //==========================================================================
 Mesh* create_bottom_pml_mesh(Mesh* bulk_mesh_pt,
                              const unsigned& bottom_boundary_id,
                              const unsigned& n_y_bottom_pml,
                              const double& width_y_bottom_pml)
 {
  // Look at the bottom boundary of the triangular mesh
  unsigned n_bottom_boundary_node = 
   bulk_mesh_pt -> nboundary_node(bottom_boundary_id);
 
  // Create a vector of ordered boundary nodes
  Vector<Node*> ordered_bottom_boundary_node_pt(n_bottom_boundary_node);

  // Fill the vector with the nodes on the respective boundary
  for (unsigned n=0; n<n_bottom_boundary_node; n++)
   {
    ordered_bottom_boundary_node_pt[n] = 
     bulk_mesh_pt->boundary_node_pt(bottom_boundary_id,n);
   }

  // Sort them from highest to lowest (in x coordinate)
  std::sort(ordered_bottom_boundary_node_pt.begin(),
            ordered_bottom_boundary_node_pt.end(),
            sorter_bottom_boundary);
 
  // The number of elements in y is taken from the triangular mesh
  unsigned n_x_bottom_pml = 
   bulk_mesh_pt -> nboundary_element(bottom_boundary_id);
 
  // Specific PML sizes needed, taken directly from physical domain
  double l_pml_bottom_x_start = 
   ordered_bottom_boundary_node_pt[n_bottom_boundary_node-1] -> x(0);
  double l_pml_bottom_x_end   = 
   ordered_bottom_boundary_node_pt[0] -> x(0);
  /// \short PML layer width subtracted from the bulk mesh lower 
  /// boundary coordinate
  double l_pml_bottom_y_start = 
   - width_y_bottom_pml 
   + ordered_bottom_boundary_node_pt[0] -> x(1);
  double l_pml_bottom_y_end   = 
   ordered_bottom_boundary_node_pt[0] -> x(1);

  unsigned bottom_quadPML_boundary_id = 2;

  /// \short Get the number of nodes along the (bulk) element edges
  /// from the first element
  unsigned nnode_1d=bulk_mesh_pt->finite_element_pt(0)->nnode_1d();

  // Create the mesh to be designated to the PML
  Mesh* pml_bottom_mesh_pt = 0;

  // Build the bottom PML mesh
  if (nnode_1d==3)
   {
    pml_bottom_mesh_pt=
     new PMLQuadMesh<QGeneralisedHelmholtzElement<2,3> > 
     (bulk_mesh_pt, bottom_boundary_id, bottom_quadPML_boundary_id,  
      n_x_bottom_pml, n_y_bottom_pml, 
      l_pml_bottom_x_start, l_pml_bottom_x_end, 
      l_pml_bottom_y_start, l_pml_bottom_y_end);
   }
  else if (nnode_1d==4)
   {
    pml_bottom_mesh_pt=
     new PMLQuadMesh<QGeneralisedHelmholtzElement<2,4> > 
     (bulk_mesh_pt, bottom_boundary_id, bottom_quadPML_boundary_id,  
      n_x_bottom_pml, n_y_bottom_pml, 
      l_pml_bottom_x_start, l_pml_bottom_x_end, 
      l_pml_bottom_y_start, l_pml_bottom_y_end); 
   }
  else {
   std::ostringstream error_stream;
   error_stream << "Currently I can't build PML meshes for   "
                << std::endl
                <<" bulk meshes with nnode_1d = " << nnode_1d
                << std::endl
                << "nodes along the element edge. Please fix this...  "
                << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
  
  // Enable PML damping on the entire mesh
  unsigned n_element_pml_bottom = pml_bottom_mesh_pt->nelement();
  for(unsigned e=0;e<n_element_pml_bottom;e++)
   {
    // Upcast
    GeneralisedHelmholtzBase<2>* el_pt = 
     dynamic_cast<GeneralisedHelmholtzBase<2>*>
     (pml_bottom_mesh_pt->element_pt(e));
    el_pt -> enable_pml(1, l_pml_bottom_y_end, l_pml_bottom_y_start);
   }

  /// \short Exterior boundary needs to be set to Dirichlet 0
  /// in both real and imaginary parts
  unsigned n_bound_pml_bottom = pml_bottom_mesh_pt->nboundary();
  for(unsigned b=0;b<n_bound_pml_bottom;b++)
   {
    unsigned n_node = pml_bottom_mesh_pt -> nboundary_node(b);
    for (unsigned n=0;n<n_node;n++)
     {
      if (b==0) {
       pml_bottom_mesh_pt -> boundary_node_pt(b,n)->pin(0);
       pml_bottom_mesh_pt -> boundary_node_pt(b,n)->pin(1); 
      }
     }
   }
  
  for(unsigned i=0;i<n_bound_pml_bottom;i++)
   {
    // How many nodes are there on this boundary?
    unsigned n_node = pml_bottom_mesh_pt->nboundary_node(i);
    
    // Loop over the nodes on boundary
    for (unsigned n=0;n<n_node;n++)
     {
      // Get pointer to node
      Node* nod_pt=pml_bottom_mesh_pt->boundary_node_pt(i,n);
      if (i==0)
       { 
        nod_pt->set_value(0,0.0);
        nod_pt->set_value(1,0.0);
       }
     }
   }

  /// \short Return the finalized mesh, with PML enabled 
  /// and boundary conditions added
  return pml_bottom_mesh_pt;
 }

 //==========================================================================
 /// \short "Constructor" for PML top right corner mesh,
 /// aligned with the existing PML meshes
 //==========================================================================
Mesh* create_top_right_pml_mesh(Mesh* pml_right_mesh_pt, 
                                Mesh* pml_top_mesh_pt, 
                                Mesh* bulk_mesh_pt,
                                const unsigned& right_boundary_id)
 {

  /// \short Relevant boundary id's to be used in construction
  /// Parent id refers to already existing PML meshes
  unsigned parent_boundary_x_id = 2;
  unsigned parent_boundary_y_id = 1;
  // Current id refers to the mesh that is to be constructed
  unsigned current_boundary_x_id = 0;
  unsigned current_boundary_y_id = 3;

  // Look at the right boundary of the triangular mesh
  unsigned n_right_boundary_node = 
   bulk_mesh_pt -> nboundary_node(right_boundary_id);
 
  // Create a vector of ordered boundary nodes
  Vector<Node*> ordered_right_boundary_node_pt(n_right_boundary_node);

  // Fill the vector with the nodes on the respective boundary
  for (unsigned n=0; n<n_right_boundary_node; n++)
   {
    ordered_right_boundary_node_pt[n] = 
     bulk_mesh_pt->boundary_node_pt(right_boundary_id,n);
   }

  // Sort them from lowest to highest (in y coordinate)
  std::sort(ordered_right_boundary_node_pt.begin(),
            ordered_right_boundary_node_pt.end(),
            sorter_right_boundary);

  /// \short Number of elements and boundary nodes to be acted upon during
  /// construction are extracted from the 'parent' PML meshes
  unsigned n_x_right_pml = 
   pml_right_mesh_pt -> nboundary_element(parent_boundary_x_id);
  unsigned n_y_top_pml = 
   pml_top_mesh_pt -> nboundary_element(parent_boundary_y_id);
  unsigned n_x_boundary_nodes = 
   pml_right_mesh_pt -> nboundary_node(parent_boundary_x_id);
  unsigned n_y_boundary_nodes = 
   pml_top_mesh_pt -> nboundary_node(parent_boundary_y_id);

  /// \short Specific PML sizes needed, taken directly from physical domain
  /// and existing PML meshes
  double l_pml_right_x_start = 
   ordered_right_boundary_node_pt[n_right_boundary_node-1] -> x(0);
  double l_pml_right_x_end =
   pml_right_mesh_pt -> 
   boundary_node_pt(parent_boundary_x_id, n_x_boundary_nodes-1)-> x(0);
  double l_pml_top_y_start = 
   ordered_right_boundary_node_pt[n_right_boundary_node-1] -> x(1);
  double l_pml_top_y_end = 
   pml_top_mesh_pt -> 
   boundary_node_pt(parent_boundary_y_id, n_y_boundary_nodes-1)-> x(1);

  /// \short Get the number of nodes along the (bulk) element edges
  /// from the first element
  unsigned nnode_1d=bulk_mesh_pt->finite_element_pt(0)->nnode_1d();

  // Create the mesh to be designated to the PML
  Mesh* pml_top_right_mesh_pt = 0;

  // Build the top right corner PML mesh
  if (nnode_1d==3)
   {
    pml_top_right_mesh_pt = 
     new PMLCornerQuadMesh<QGeneralisedHelmholtzElement<2,3> >
     (pml_right_mesh_pt, pml_top_mesh_pt, bulk_mesh_pt, 
      ordered_right_boundary_node_pt[n_right_boundary_node-1],
      parent_boundary_x_id, parent_boundary_y_id, 
      current_boundary_x_id, current_boundary_y_id,  
      n_x_right_pml, n_y_top_pml, 
      l_pml_right_x_start, l_pml_right_x_end, 
      l_pml_top_y_start, l_pml_top_y_end);
   }
  else if (nnode_1d==4) 
   {
    pml_top_right_mesh_pt = 
     new PMLCornerQuadMesh<QGeneralisedHelmholtzElement<2,4> >
     (pml_right_mesh_pt, pml_top_mesh_pt, bulk_mesh_pt, 
      ordered_right_boundary_node_pt[n_right_boundary_node-1],
      parent_boundary_x_id, parent_boundary_y_id, 
      current_boundary_x_id, current_boundary_y_id,  
      n_x_right_pml, n_y_top_pml, 
      l_pml_right_x_start, l_pml_right_x_end, 
      l_pml_top_y_start, l_pml_top_y_end); 
   }
  else {
   std::ostringstream error_stream;
   error_stream << "Currently I can't build PML meshes for   "
                << std::endl
                <<" bulk meshes with nnode_1d = " << nnode_1d
                << std::endl
                << "nodes along the element edge. Please fix this...  "
                << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

  // Enable PML damping on the entire mesh
  /// \short The enabling must be perfromed in both x- and y-directions
  /// as this is a corner PML mesh
  unsigned n_element_pml_top_right = pml_top_right_mesh_pt->nelement();
  for(unsigned e=0;e<n_element_pml_top_right;e++)
   {
    // Upcast
    GeneralisedHelmholtzBase<2>* el_pt = 
     dynamic_cast<GeneralisedHelmholtzBase<2>*>     
     (pml_top_right_mesh_pt->element_pt(e));
    el_pt -> enable_pml(0, l_pml_right_x_start, l_pml_right_x_end);
    el_pt -> enable_pml(1, l_pml_top_y_start, l_pml_top_y_end);
   }
  
  /// \short Exterior boundary needs to be set to Dirichlet 0
  /// in both real and imaginary parts
  unsigned n_bound_pml_top_right = pml_top_right_mesh_pt->nboundary();
  for(unsigned b=0;b<n_bound_pml_top_right;b++)
   {
    unsigned n_node = pml_top_right_mesh_pt -> nboundary_node(b);
    for (unsigned n=0;n<n_node;n++)
     {
      if ((b==1)||(b==2)) {
       pml_top_right_mesh_pt -> boundary_node_pt(b,n)->pin(0);
       pml_top_right_mesh_pt -> boundary_node_pt(b,n)->pin(1); 
      }
     }
   }

  for(unsigned i=0;i<n_bound_pml_top_right;i++)
   {
    // How many nodes are there on this boundary?
    unsigned n_node = pml_top_right_mesh_pt->nboundary_node(i);
    
    // Loop over the nodes on boundary
    for (unsigned n=0;n<n_node;n++)
     {
      // Get pointer to node
      Node* nod_pt=pml_top_right_mesh_pt->boundary_node_pt(i,n);
      if ((i==1)||(i==2))
       { 
        nod_pt->set_value(0,0.0);
        nod_pt->set_value(1,0.0);
       }
     }
   }
  
  /// \short Return the finalized mesh, with PML enabled 
  /// and boundary conditions added
  return pml_top_right_mesh_pt;
 }

 //==========================================================================
 /// \short  "Constructor" for PML bottom right corner mesh, 
 /// aligned with the existing PML meshes
 //==========================================================================
 Mesh* create_bottom_right_pml_mesh(Mesh* pml_right_mesh_pt, 
                                    Mesh* pml_bottom_mesh_pt, 
                                    Mesh* bulk_mesh_pt,
                                    const unsigned& right_boundary_id)
 {

  /// \short Relevant boundary id's to be used in construction
  /// Parent id refers to already existing PML meshes
  unsigned parent_boundary_x_id = 0;
  unsigned parent_boundary_y_id = 1;
  // Current id refers to the mesh that is to be constructed
  unsigned current_boundary_x_id = 2;
  unsigned current_boundary_y_id = 3;

  // Look at the right boundary of the triangular mesh
  unsigned n_right_boundary_node = 
   bulk_mesh_pt -> nboundary_node(right_boundary_id);
 
  // Create a vector of ordered boundary nodes
  Vector<Node*> ordered_right_boundary_node_pt(n_right_boundary_node);

  // Fill the vector with the nodes on the respective boundary
  for (unsigned n=0; n<n_right_boundary_node; n++)
   {
    ordered_right_boundary_node_pt[n] = 
     bulk_mesh_pt->boundary_node_pt(right_boundary_id,n);
   }

  // Sort them from lowest to highest (in y coordinate)
  std::sort(ordered_right_boundary_node_pt.begin(),
            ordered_right_boundary_node_pt.end(),
            sorter_right_boundary);

  /// \short Number of elements and boundary nodes to be acted upon during
  /// construction are extracted from the 'parent' PML meshes
  unsigned n_x_right_pml = 
   pml_right_mesh_pt -> nboundary_element(parent_boundary_x_id);
  unsigned n_y_bottom_pml = 
   pml_bottom_mesh_pt -> nboundary_element(parent_boundary_y_id);
  unsigned n_x_boundary_nodes = 
   pml_right_mesh_pt -> nboundary_node(parent_boundary_x_id);

  /// \short Specific PML sizes needed, taken directly from physical domain
  /// and existing PML meshes
  double l_pml_right_x_start = 
   ordered_right_boundary_node_pt[0] -> x(0);
  double l_pml_right_x_end = 
   pml_right_mesh_pt -> 
   boundary_node_pt(parent_boundary_x_id, n_x_boundary_nodes-1)-> x(0);
  double l_pml_bottom_y_start = 
   pml_bottom_mesh_pt -> boundary_node_pt(parent_boundary_y_id, 0)-> x(1);
  double l_pml_bottom_y_end = 
   ordered_right_boundary_node_pt[0] -> x(1);

  /// \short Get the number of nodes along the (bulk) element edges
  /// from the first element
  unsigned nnode_1d=bulk_mesh_pt->finite_element_pt(0)->nnode_1d();

  // Create the mesh to be designated to the PML
  Mesh* pml_bottom_right_mesh_pt = 0;

  // Build the bottom right corner PML mesh
  if (nnode_1d==3)
   {
    pml_bottom_right_mesh_pt=
     new PMLCornerQuadMesh<QGeneralisedHelmholtzElement<2,3> >
     (pml_right_mesh_pt, pml_bottom_mesh_pt, bulk_mesh_pt, 
      ordered_right_boundary_node_pt[0],
      parent_boundary_x_id, parent_boundary_y_id, 
      current_boundary_x_id, current_boundary_y_id,  
      n_x_right_pml, n_y_bottom_pml, 
      l_pml_right_x_start, l_pml_right_x_end, 
      l_pml_bottom_y_start, l_pml_bottom_y_end);    
   }
  else if (nnode_1d==4)
   {
    pml_bottom_right_mesh_pt=
     new PMLCornerQuadMesh<QGeneralisedHelmholtzElement<2,4> >
     (pml_right_mesh_pt, pml_bottom_mesh_pt, bulk_mesh_pt, 
      ordered_right_boundary_node_pt[0],
      parent_boundary_x_id, parent_boundary_y_id, 
      current_boundary_x_id, current_boundary_y_id,  
      n_x_right_pml, n_y_bottom_pml, 
      l_pml_right_x_start, l_pml_right_x_end, 
      l_pml_bottom_y_start, l_pml_bottom_y_end);  
   }
  else {
   std::ostringstream error_stream;
   error_stream << "Currently I can't build PML meshes for   "
                << std::endl
                <<" bulk meshes with nnode_1d = " << nnode_1d
                << std::endl
                << "nodes along the element edge. Please fix this...  "
                << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
  
  // Enable PML damping on the entire mesh
  /// \short The enabling must be perfromed in both x- and y-directions
  /// as this is a corner PML mesh
  unsigned n_element_pml_bottom_right = 
   pml_bottom_right_mesh_pt->nelement();

  for(unsigned e=0;e<n_element_pml_bottom_right;e++)
   {
    // Upcast
    GeneralisedHelmholtzBase<2>* el_pt = 
     dynamic_cast<GeneralisedHelmholtzBase<2>*>
     (pml_bottom_right_mesh_pt->element_pt(e));
    el_pt -> enable_pml(0, l_pml_right_x_start, l_pml_right_x_end);
    el_pt -> enable_pml(1, l_pml_bottom_y_end, l_pml_bottom_y_start);
   }

  /// \short Exterior boundary needs to be set to Dirichlet 0
  /// in both real and imaginary parts
  unsigned n_bound_pml_bottom_right = 
   pml_bottom_right_mesh_pt->nboundary();
  
  for(unsigned b=0;b<n_bound_pml_bottom_right;b++)
   {
    unsigned n_node = pml_bottom_right_mesh_pt -> nboundary_node(b);
    for (unsigned n=0;n<n_node;n++)
     {
      if ((b==0)||(b==1)) {
       pml_bottom_right_mesh_pt -> boundary_node_pt(b,n)->pin(0);
       pml_bottom_right_mesh_pt -> boundary_node_pt(b,n)->pin(1); 
      }
     }
   }

  for(unsigned i=0;i<n_bound_pml_bottom_right;i++)
   {
    // How many nodes are there on this boundary?
    unsigned n_node = pml_bottom_right_mesh_pt->nboundary_node(i);
    
    // Loop over the nodes on boundary
    for (unsigned n=0;n<n_node;n++)
     {
      // Get pointer to node
      Node* nod_pt=pml_bottom_right_mesh_pt->boundary_node_pt(i,n);
      if ((i==0)||(i==1))
       { 
        nod_pt->set_value(0,0.0);
        nod_pt->set_value(1,0.0);
       }
     }
   }
  
  /// \short Return the finalized mesh, with PML enabled 
  /// and boundary conditions added
  return pml_bottom_right_mesh_pt;
 }

 //==========================================================================
 /// \short "Constructor" for PML top left corner mesh, 
 /// aligned with the existing PML meshes
 //==========================================================================
 Mesh* create_top_left_pml_mesh(Mesh* pml_left_mesh_pt, 
                                Mesh* pml_top_mesh_pt, 
                                Mesh* bulk_mesh_pt,
                                const unsigned& left_boundary_id)
 {

  /// \short Relevant boundary id's to be used in construction
  /// Parent id refers to already existing PML meshes
  unsigned parent_boundary_x_id = 2;
  unsigned parent_boundary_y_id = 3;
  // Current id refers to the mesh that is to be constructed
  unsigned current_boundary_x_id = 0;
  unsigned current_boundary_y_id = 1;

  // Look at the left boundary of the triangular mesh
  unsigned n_left_boundary_node = 
   bulk_mesh_pt -> nboundary_node(left_boundary_id);
 
  // Create a vector of ordered boundary nodes
  Vector<Node*> ordered_left_boundary_node_pt(n_left_boundary_node);

  // Fill the vector with the nodes on the respective boundary
  for (unsigned n=0; n<n_left_boundary_node; n++)
   {
    ordered_left_boundary_node_pt[n] = 
     bulk_mesh_pt->boundary_node_pt(left_boundary_id,n);
   }

  /// \short Sort them from lowest to highest (in y coordinate)
  /// sorter_right_boundary is still functional, as the sorting
  /// is performed by the same criterion
  std::sort(ordered_left_boundary_node_pt.begin(),
            ordered_left_boundary_node_pt.end(),
            sorter_right_boundary);

  /// \short Number of elements and boundary nodes to be acted upon during
  /// construction are extracted from the 'parent' PML meshes
  unsigned n_x_left_pml = 
   pml_left_mesh_pt -> nboundary_element(parent_boundary_x_id);
  unsigned n_y_top_pml = 
   pml_top_mesh_pt -> nboundary_element(parent_boundary_y_id);
  unsigned n_y_boundary_nodes = 
   pml_top_mesh_pt -> nboundary_node(parent_boundary_y_id);

  /// \short Specific PML sizes needed, taken directly from physical domain
  /// and existing PML meshes
  double l_pml_left_x_start =  
   pml_left_mesh_pt -> boundary_node_pt(parent_boundary_x_id, 0)-> x(0);
  double l_pml_left_x_end = 
   ordered_left_boundary_node_pt[n_left_boundary_node-1] -> x(0);
  double l_pml_top_y_start = 
   ordered_left_boundary_node_pt[n_left_boundary_node-1] -> x(1);
  double l_pml_top_y_end = 
   pml_top_mesh_pt -> 
   boundary_node_pt(parent_boundary_y_id, n_y_boundary_nodes-1)-> x(1);

  /// \short Get the number of nodes along the (bulk) element edges
  /// from the first element
  unsigned nnode_1d=bulk_mesh_pt->finite_element_pt(0)->nnode_1d();

  // Create the mesh to be designated to the PML
  Mesh* pml_top_left_mesh_pt = 0;

  // Build the top left corner PML mesh
  if (nnode_1d==3)
   {
    pml_top_left_mesh_pt=
     new PMLCornerQuadMesh<QGeneralisedHelmholtzElement<2,3> >
     (pml_left_mesh_pt, pml_top_mesh_pt, bulk_mesh_pt, 
      ordered_left_boundary_node_pt[n_left_boundary_node-1], 
      parent_boundary_x_id, parent_boundary_y_id, 
      current_boundary_x_id, current_boundary_y_id,  
      n_x_left_pml, n_y_top_pml, 
      l_pml_left_x_start, l_pml_left_x_end, 
      l_pml_top_y_start, l_pml_top_y_end);
   }
  else if (nnode_1d==4)
   {
    pml_top_left_mesh_pt=
     new PMLCornerQuadMesh<QGeneralisedHelmholtzElement<2,4> >
     (pml_left_mesh_pt, pml_top_mesh_pt, bulk_mesh_pt, 
      ordered_left_boundary_node_pt[n_left_boundary_node-1], 
      parent_boundary_x_id, parent_boundary_y_id, 
      current_boundary_x_id, current_boundary_y_id,  
      n_x_left_pml, n_y_top_pml, 
      l_pml_left_x_start, l_pml_left_x_end, 
      l_pml_top_y_start, l_pml_top_y_end);
   }
  else {
    std::ostringstream error_stream;
    error_stream << "Currently I can't build PML meshes for   "
                 << std::endl
                 <<" bulk meshes with nnode_1d = " << nnode_1d
                 << std::endl
                 << "nodes along the element edge. Please fix this...  "
                 << std::endl;
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  // Enable PML damping on the entire mesh
  /// \short The enabling must be perfromed in both x- and y-directions
  /// as this is a corner PML mesh
  unsigned n_element_pml_top_left = pml_top_left_mesh_pt->nelement();

  for(unsigned e=0;e<n_element_pml_top_left;e++)
   {
    // Upcast
    GeneralisedHelmholtzBase<2>* el_pt = 
     dynamic_cast<GeneralisedHelmholtzBase<2>*>
     (pml_top_left_mesh_pt->element_pt(e));
    el_pt -> enable_pml(0, l_pml_left_x_end, l_pml_left_x_start);
    el_pt -> enable_pml(1, l_pml_top_y_start, l_pml_top_y_end);
   }

  /// \short Exterior boundary needs to be set to Dirichlet 0
  /// in both real and imaginary parts
  unsigned n_bound_pml_top_left = pml_top_left_mesh_pt->nboundary();
  for(unsigned b=0;b<n_bound_pml_top_left;b++)
   {
    unsigned n_node = pml_top_left_mesh_pt -> nboundary_node(b);
    for (unsigned n=0;n<n_node;n++)
     {
     if ((b==2)||(b==3)) {
      pml_top_left_mesh_pt -> boundary_node_pt(b,n)->pin(0);
      pml_top_left_mesh_pt -> boundary_node_pt(b,n)->pin(1); 
     }
     }
   }
  
  for(unsigned i=0;i<n_bound_pml_top_left;i++)
   {
    // How many nodes are there on this boundary?
    unsigned n_node = pml_top_left_mesh_pt->nboundary_node(i);
    
    // Loop over the nodes on boundary
    for (unsigned n=0;n<n_node;n++)
     {
      // Get pointer to node
      Node* nod_pt=pml_top_left_mesh_pt->boundary_node_pt(i,n);
      if ((i==2)||(i==3))
       { 
        nod_pt->set_value(0,0.0);
        nod_pt->set_value(1,0.0);
       }
     }
   }

  /// \short Return the finalized mesh, with PML enabled 
  /// and boundary conditions added
  return pml_top_left_mesh_pt;
 }

 //==========================================================================
 /// \short "Constructor" for PML bottom left corner mesh, 
 /// aligned with the existing PML meshes
 //==========================================================================
 Mesh* create_bottom_left_pml_mesh(Mesh* pml_left_mesh_pt, 
                                   Mesh* pml_bottom_mesh_pt, 
                                   Mesh* bulk_mesh_pt,
                                   const unsigned& left_boundary_id)
 {

  /// \short Relevant boundary id's to be used in construction
  /// Parent id refers to already existing PML meshes
  unsigned parent_boundary_x_id = 0;
  unsigned parent_boundary_y_id = 3;
  // Current id refers to the mesh that is to be constructed
  unsigned current_boundary_x_id = 2;
  unsigned current_boundary_y_id = 1;

  // Look at the left boundary of the triangular mesh
  unsigned n_left_boundary_node = 
   bulk_mesh_pt -> nboundary_node(left_boundary_id);
 
  // Create a vector of ordered boundary nodes
  Vector<Node*> ordered_left_boundary_node_pt(n_left_boundary_node);

  // Fill the vector with the nodes on the respective boundary
  for (unsigned n=0; n<n_left_boundary_node; n++)
   {
    ordered_left_boundary_node_pt[n] = 
     bulk_mesh_pt->boundary_node_pt(left_boundary_id,n);
   }

  /// \short Sort them from lowest to highest (in y coordinate)
  /// sorter_right_boundary is still functional, as the sorting
  /// is performed by the same criterion
  std::sort(ordered_left_boundary_node_pt.begin(),
            ordered_left_boundary_node_pt.end(),
            sorter_right_boundary);

  /// \short Number of elements and boundary nodes to be acted upon during
  /// construction are extracted from the 'parent' PML meshes
  unsigned n_x_left_pml = 
   pml_left_mesh_pt -> nboundary_element(parent_boundary_x_id);
  unsigned n_y_bottom_pml = 
   pml_bottom_mesh_pt -> nboundary_element(parent_boundary_y_id);

  /// \short Specific PML sizes needed, taken directly from physical domain
  /// and existing PML meshes
  double l_pml_left_x_start =  
   pml_left_mesh_pt -> boundary_node_pt(parent_boundary_x_id, 0)-> x(0);
  double l_pml_left_x_end = 
   ordered_left_boundary_node_pt[n_left_boundary_node-1] -> x(0);
  double l_pml_bottom_y_start =
   pml_bottom_mesh_pt -> boundary_node_pt(parent_boundary_y_id, 0)-> x(1);
  double l_pml_bottom_y_end = 
   ordered_left_boundary_node_pt[0] -> x(1);

  /// \short Get the number of nodes along the (bulk) element edges
  /// from the first element
  unsigned nnode_1d=bulk_mesh_pt->finite_element_pt(0)->nnode_1d();

  // Create the mesh to be designated to the PML
  Mesh* pml_bottom_left_mesh_pt = 0;

  // Build the bottom left corner PML mesh
  if (nnode_1d==3)
   {
    pml_bottom_left_mesh_pt=
     new PMLCornerQuadMesh<QGeneralisedHelmholtzElement<2,3> >
     (pml_left_mesh_pt, pml_bottom_mesh_pt, bulk_mesh_pt, 
      ordered_left_boundary_node_pt[0], 
      parent_boundary_x_id, parent_boundary_y_id, 
      current_boundary_x_id, current_boundary_y_id,  
      n_x_left_pml, n_y_bottom_pml, 
      l_pml_left_x_start, l_pml_left_x_end, 
      l_pml_bottom_y_start, l_pml_bottom_y_end);
   }
  else if (nnode_1d==4) 
   {
    pml_bottom_left_mesh_pt=
     new PMLCornerQuadMesh<QGeneralisedHelmholtzElement<2,4> >
     (pml_left_mesh_pt, pml_bottom_mesh_pt, bulk_mesh_pt, 
      ordered_left_boundary_node_pt[0], 
      parent_boundary_x_id, parent_boundary_y_id, 
      current_boundary_x_id, current_boundary_y_id,  
      n_x_left_pml, n_y_bottom_pml, 
      l_pml_left_x_start, l_pml_left_x_end, 
      l_pml_bottom_y_start, l_pml_bottom_y_end);
   }
  else {
   std::ostringstream error_stream;
   error_stream << "Currently I can't build PML meshes for   "
                << std::endl
                <<" bulk meshes with nnode_1d = " << nnode_1d
                << std::endl
                << "nodes along the element edge. Please fix this...  "
                << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
  
  //Enable PML damping on the entire mesh
  /// \short The enabling must be perfromed in both x- and y-directions
  /// as this is a corner PML mesh
  unsigned n_element_pml_bottom_left = pml_bottom_left_mesh_pt->nelement();
  for(unsigned e=0;e<n_element_pml_bottom_left;e++)
   {
    // Upcast
    GeneralisedHelmholtzBase<2>* el_pt = 
     dynamic_cast<GeneralisedHelmholtzBase<2>*>
     (pml_bottom_left_mesh_pt->element_pt(e));
    el_pt -> enable_pml(0, l_pml_left_x_end, l_pml_left_x_start);
    el_pt -> enable_pml(1, l_pml_bottom_y_end, l_pml_bottom_y_start);
   }

  /// \short Exterior boundary needs to be set to Dirichlet 0
  /// in both real and imaginary parts
  unsigned n_bound_pml_bottom_left = pml_bottom_left_mesh_pt->nboundary();
  for(unsigned b=0;b<n_bound_pml_bottom_left;b++)
   {
    unsigned n_node = pml_bottom_left_mesh_pt -> nboundary_node(b);
    for (unsigned n=0;n<n_node;n++)
     {
      if ((b==0)||(b==3)) {
       pml_bottom_left_mesh_pt -> boundary_node_pt(b,n)->pin(0);
       pml_bottom_left_mesh_pt -> boundary_node_pt(b,n)->pin(1); 
      }
     }
   }
  
  for(unsigned i=0;i<n_bound_pml_bottom_left;i++)
   {
    // How many nodes are there on this boundary?
    unsigned n_node = pml_bottom_left_mesh_pt->nboundary_node(i);
    
    // Loop over the nodes on boundary
    for (unsigned n=0;n<n_node;n++)
     {
      // Get pointer to node
      Node* nod_pt=pml_bottom_left_mesh_pt->boundary_node_pt(i,n);
      if ((i==0)||(i==3))
       { 
        nod_pt->set_value(0,0.0);
        nod_pt->set_value(1,0.0);
       }
     }
   }
  
  /// \short Return the finalized mesh, with PML enabled 
  /// and boundary conditions added
  return pml_bottom_left_mesh_pt;
 }
 
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//======================================================================
/// Set the data for the number of Variables at each node, always two
/// (real and imag part) in every case
//======================================================================
 template<unsigned DIM, unsigned NNODE_1D>
 const unsigned QGeneralisedHelmholtzElement<DIM,NNODE_1D>::Initial_Nvalue = 2;


//======================================================================
/// Compute element residual Vector and/or element Jacobian matrix 
/// 
/// flag=1: compute both
/// flag=0: compute only residual Vector
///
/// Pure version without hanging nodes
//======================================================================
template <unsigned DIM>
void  GeneralisedHelmholtzEquations<DIM>::
fill_in_generic_residual_contribution_helmholtz(Vector<double> &residuals, 
                                                DenseMatrix<double> &jacobian, 
                                                const unsigned& flag) 
{
 //Find out how many nodes there are
 const unsigned n_node = nnode();

 //Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);
 
 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();
 
 //Integers to store the local equation and unknown numbers
 int local_eqn_real=0, local_unknown_real=0;
 int local_eqn_imag=0, local_unknown_imag=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   //Call the derivatives of the shape and test functions
   double J = dshape_and_dtest_eulerian_at_knot_helmholtz(ipt,psi,dpsidx,
                                                          test,dtestdx);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   //Calculate local values of unknown
   //Allocate and initialise to zero
   std::complex<double> interpolated_u(0.0,0.0);
   Vector<double> interpolated_x(DIM,0.0);
   Vector< std::complex<double> > interpolated_dudx(DIM);
   
   //Calculate function value and derivatives:
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     // Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       interpolated_x[j] += raw_nodal_position(l,j)*psi(l);
      }
     
     //Get the nodal value of the helmholtz unknown
     const std::complex<double> 
      u_value(raw_nodal_value(l,u_index_helmholtz().real()),
              raw_nodal_value(l,u_index_helmholtz().imag()));

     //Add to the interpolated value
     interpolated_u += u_value*psi(l);
     
     // Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       interpolated_dudx[j] += u_value*dpsidx(l,j);
      }
    }
   
   //Get source function
   //-------------------
   std::complex<double> source(0.0,0.0);
   get_source_helmholtz(ipt,interpolated_x,source);

   // Get the value of the soundspeed at the respective point
   double c(0.0); 
   get_c_helmholtz(ipt,interpolated_x,c);
   
   // Get the value of the absorption at the respective point
   double alpha(0.0); 
   get_alpha_helmholtz(ipt,interpolated_x,alpha);

   /// \short Use soundspeed, absorption and frequency 
   /// to precompute expressions to be subsequently used
   double c_squared = c * c;
   double alpha_squared = alpha * alpha;
   double omega_squared = omega() * omega();
   // First expression containing c, alpha and omega
   double cao_expr1 = 2.0 * c * alpha * omega();
   // Second expression containing c and alpha 
   double ca_expr2 = alpha_squared * c_squared;

   // Declare a vector of complex numbers for pml weights on the Laplace bit
   Vector< std::complex<double> > pml_stiffness_weight(DIM);
   // Declare a complex number for pml weights on the mass matrix bit
   std::complex<double> pml_mass_weight = std::complex<double>(1.0,0.0);

   /// \short All the PML weights that participate in the assemby process 
   /// are computed here. pml_stiffness_weight will contain the entries
   /// for the Laplace bit, while pml_mass_weight contains the contributions
   /// to the Helmholtz bit. Both default to 1.0, should the PML not be 
   /// enabled via enable_pml.
   compute_pml_coefficients(ipt, interpolated_x, 
                            pml_stiffness_weight, 
                            pml_mass_weight); 

   // Assemble residuals and Jacobian
   //--------------------------------
   // Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {

     // first, compute the real part contribution 
     //-------------------------------------------  
    
     //Get the local equation
     local_eqn_real = nodal_local_eqn(l,u_index_helmholtz().real());
     local_eqn_imag = nodal_local_eqn(l,u_index_helmholtz().imag());

     /*IF it's not a boundary condition*/
     if(local_eqn_real >= 0)
      {
       // Add body force/source term and Helmholtz bit    
       residuals[local_eqn_real] += 
        ( source.real() - ( 
           pml_mass_weight.real() * ( 
            (omega_squared-ca_expr2)*interpolated_u.real() - 
             cao_expr1 * interpolated_u.imag() ) 
          -pml_mass_weight.imag() * ( 
            (omega_squared-ca_expr2)*interpolated_u.imag() + 
             cao_expr1 * interpolated_u.real() ) 
        ) )*test(l)*W;
       
       // The Laplace bit
       for(unsigned k=0;k<DIM;k++)
        {
         residuals[local_eqn_real] += c_squared * (
          pml_stiffness_weight[k].real() * interpolated_dudx[k].real()
         -pml_stiffness_weight[k].imag() * interpolated_dudx[k].imag()
         )*dtestdx(l,k)*W;
        }
       
       // Calculate the jacobian
       //-----------------------
       if(flag)
        {
         //Loop over the velocity shape functions again
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           local_unknown_real = nodal_local_eqn(l2,u_index_helmholtz().real());
           local_unknown_imag = nodal_local_eqn(l2,u_index_helmholtz().imag());

           //If at a non-zero degree of freedom add in the entry
           if(local_unknown_real >= 0)
            {
             //Add contribution to Elemental Matrix
             for(unsigned i=0;i<DIM;i++)
              {
               jacobian(local_eqn_real,local_unknown_real) 
                += c_squared * pml_stiffness_weight[i].real() 
                             * dpsidx(l2,i)*dtestdx(l,i)*W;
              }
             // Add the helmholtz contribution
             jacobian(local_eqn_real,local_unknown_real) 
              += (pml_mass_weight.imag() * cao_expr1 - 
                  pml_mass_weight.real() * (omega_squared -  ca_expr2)) 
                * psi(l2)*test(l)*W;
            }
           //If at a non-zero degree of freedom add in the entry
           if(local_unknown_imag >= 0)
            {
             //Add contribution to Elemental Matrix
             for(unsigned i=0;i<DIM;i++)
              {
               jacobian(local_eqn_real,local_unknown_imag) 
                -= c_squared * pml_stiffness_weight[i].imag() 
                             * dpsidx(l2,i)*dtestdx(l,i)*W;

              }
             // Add the helmholtz contribution
             jacobian(local_eqn_real,local_unknown_imag) 
              += (pml_mass_weight.real() * cao_expr1 + 
                  pml_mass_weight.imag() * (omega_squared -  ca_expr2)) 
                * psi(l2)*test(l)*W;
            }
          }
        }
      }
     
     // Second, compute the imaginary part contribution 
     //------------------------------------------------

     //Get the local equation
     local_eqn_imag = nodal_local_eqn(l,u_index_helmholtz().imag());
     local_eqn_real = nodal_local_eqn(l,u_index_helmholtz().real());

     /*IF it's not a boundary condition*/
     if(local_eqn_imag >= 0)
      {
       // Add body force/source term and Helmholtz bit
       residuals[local_eqn_imag] += 
        ( source.imag() - ( 
           pml_mass_weight.imag() * ( 
            (omega_squared-ca_expr2)*interpolated_u.real() - 
             cao_expr1 * interpolated_u.imag() ) 
         + pml_mass_weight.real() * ( 
            (omega_squared-ca_expr2)*interpolated_u.imag() + 
            cao_expr1 * interpolated_u.real() ) 
         ) )*test(l)*W;
       
       // The Laplace bit
       for(unsigned k=0;k<DIM;k++)
        {
         residuals[local_eqn_imag] += c_squared * (
               pml_stiffness_weight[k].imag() * interpolated_dudx[k].real()
              +pml_stiffness_weight[k].real() * interpolated_dudx[k].imag()   
         )*dtestdx(l,k)*W;
        }
       
       // Calculate the jacobian
       //-----------------------
       if(flag)
        {
         //Loop over the velocity shape functions again
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           local_unknown_imag = nodal_local_eqn(l2,u_index_helmholtz().imag());
           local_unknown_real = nodal_local_eqn(l2,u_index_helmholtz().real());

           //If at a non-zero degree of freedom add in the entry
           if(local_unknown_imag >= 0)
            {
             //Add contribution to Elemental Matrix
             for(unsigned i=0;i<DIM;i++)
              {
               jacobian(local_eqn_imag,local_unknown_imag) 
                += c_squared * pml_stiffness_weight[i].real() 
                             * dpsidx(l2,i)*dtestdx(l,i)*W;
              }
             // Add the helmholtz contribution
             jacobian(local_eqn_imag,local_unknown_imag)
              += (pml_mass_weight.imag() * cao_expr1 - 
                  pml_mass_weight.real() * (omega_squared -  ca_expr2)) 
                * psi(l2)*test(l)*W;
            }
           if(local_unknown_real >= 0)
            {
             //Add contribution to Elemental Matrix
             for(unsigned i=0;i<DIM;i++)
              {
               jacobian(local_eqn_imag,local_unknown_real) 
                += c_squared * pml_stiffness_weight[i].imag() 
                             * dpsidx(l2,i)*dtestdx(l,i)*W;
              }
             // Add the helmholtz contribution
             jacobian(local_eqn_imag,local_unknown_real)
              += (-pml_mass_weight.real() * cao_expr1 - 
                   pml_mass_weight.imag() * (omega_squared -  ca_expr2)) 
                 * psi(l2)*test(l)*W;

            }
          }
        }
      }   
    }
  } // End of loop over integration points
}   

 
//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
template <unsigned DIM>
unsigned  GeneralisedHelmholtzEquations<DIM>::self_test()
{

 bool passed=true;

 // Check lower-level stuff
 if (FiniteElement::self_test()!=0)
  {
   passed=false;
  }

 // Return verdict
 if (passed)
  {
   return 0;
  }
 else
  {
   return 1;
  }
   
}


//======================================================================
/// Output function:
///
///   x,y,u_re,u_imag   or    x,y,z,u_re,u_imag
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void  GeneralisedHelmholtzEquations<DIM>::output(std::ostream &outfile, 
                                    const unsigned &nplot)
{

 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   std::complex<double> u(interpolated_u_generalised_helmholtz(s));
   for(unsigned i=0;i<DIM;i++) 
    {
     outfile << interpolated_x(s,i) << " ";
    }
   outfile << u.real() << " " << u.imag() << std::endl;   
   
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);

}
 



//======================================================================
/// Output function for real part of full time-dependent solution
///
///  u = Re( (u_r +i u_i) exp(-i omega t)
///
/// at phase angle omega t = phi.
///
///   x,y,u   or    x,y,z,u
///
/// Output at nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void  GeneralisedHelmholtzEquations<DIM>::output_real(std::ostream &outfile, 
                                           const double& phi,
                                           const unsigned &nplot)
{
 
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   std::complex<double> u(interpolated_u_generalised_helmholtz(s));
   for(unsigned i=0;i<DIM;i++) 
    {
     outfile << interpolated_x(s,i) << " ";
    }
   outfile << u.real()*cos(phi)+u.imag()*sin(phi) << std::endl;   
   
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);

}

//======================================================================
/// Output function for imaginary part of full time-dependent solution
///
///  u = Im( (u_r +i u_i) exp(-i omega t))
///
/// at phase angle omega t = phi.
///
///   x,y,u   or    x,y,z,u
///
/// Output at nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void  GeneralisedHelmholtzEquations<DIM>::output_imag(std::ostream &outfile, 
                                           const double& phi,
                                           const unsigned &nplot)
{
 
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   std::complex<double> u(interpolated_u_generalised_helmholtz(s));
   for(unsigned i=0;i<DIM;i++) 
    {
     outfile << interpolated_x(s,i) << " ";
    }
   outfile << u.imag()*cos(phi)-u.real()*sin(phi) << std::endl;   
   
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);

}
 
 
//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void  GeneralisedHelmholtzEquations<DIM>::output(FILE* file_pt,
                                    const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Tecplot header info
 fprintf(file_pt,"%s",tecplot_zone_string(nplot).c_str());

 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   std::complex<double> u(interpolated_u_generalised_helmholtz(s));

   for(unsigned i=0;i<DIM;i++) 
    {
     fprintf(file_pt,"%g ",interpolated_x(s,i));
    }

   for(unsigned i=0;i<DIM;i++) 
    {
     fprintf(file_pt,"%g ",interpolated_x(s,i));
    }
   fprintf(file_pt,"%g ",u.real());
   fprintf(file_pt,"%g \n",u.imag());

  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(file_pt,nplot);
}



//======================================================================
 /// Output exact solution
 /// 
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points.
 ///
 ///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
template <unsigned DIM>
void GeneralisedHelmholtzEquations<DIM>::output_fct(std::ostream &outfile, 
                                       const unsigned &nplot, 
                  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
  
  // Vector for coordintes
  Vector<double> x(DIM);
  
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Exact solution Vector 
 Vector<double> exact_soln(2);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,u_exact
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   outfile << exact_soln[0] << " " <<  exact_soln[1] << std::endl;  
  }
 
 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);
}



//======================================================================
/// Output function for real part of full time-dependent fct
///
///  u = Re( (u_r +i u_i) exp(-i omega t)
///
/// at phase angle omega t = phi.
///
///   x,y,u   or    x,y,z,u
///
/// Output at nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void GeneralisedHelmholtzEquations<DIM>::output_real_fct(
 std::ostream &outfile, 
 const double& phi,
 const unsigned &nplot, 
 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
  
  // Vector for coordintes
  Vector<double> x(DIM);
  
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Exact solution Vector 
 Vector<double> exact_soln(2);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,u_exact
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   outfile << exact_soln[0]*cos(phi)+ exact_soln[1]*sin(phi) << std::endl;  
  }
 
 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);
}

//======================================================================
/// Output function for imaginary part of full time-dependent fct
///
///  u = Im( (u_r +i u_i) exp(-i omega t))
///
/// at phase angle omega t = phi.
///
///   x,y,u   or    x,y,z,u
///
/// Output at nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void GeneralisedHelmholtzEquations<DIM>::output_imag_fct(
 std::ostream &outfile, 
 const double& phi,
 const unsigned &nplot, 
 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
  
  // Vector for coordintes
  Vector<double> x(DIM);
  
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Exact solution Vector 
 Vector<double> exact_soln(2);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,u_exact
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   outfile << exact_soln[1]*cos(phi) - exact_soln[0]*sin(phi) << std::endl;  
  }
 
 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);
}




//======================================================================
 /// Validate against exact solution
 /// 
 /// Solution is provided via function pointer.
 /// Plot error at a given number of plot points.
 ///
//======================================================================
template <unsigned DIM>
void GeneralisedHelmholtzEquations<DIM>::compute_error(std::ostream &outfile, 
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double& error, double& norm)
{ 
 
 // Initialise
 error=0.0;
 norm=0.0;
 
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Vector for coordintes
 Vector<double> x(DIM);
 
 //Find out how many nodes there are in the element
 unsigned n_node = nnode();
 
 Shape psi(n_node);
 
 //Set the value of n_intpt
 unsigned n_intpt = integral_pt()->nweight();
  
 // Tecplot 
 outfile << "ZONE" << std::endl;
 
 // Exact solution Vector 
 Vector<double> exact_soln(2);
 
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   
   //Assign values of s
   for(unsigned i=0;i<DIM;i++)
    {
     s[i] = integral_pt()->knot(ipt,i);
    }
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   // Get jacobian of mapping
   double J=J_eulerian(s);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   // Get FE function value
   std::complex<double> u_fe=interpolated_u_generalised_helmholtz(s);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,error
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   outfile << exact_soln[0] << " "  << exact_soln[1] << " " 
           << exact_soln[0]-u_fe.real() << " " << exact_soln[1]-u_fe.imag() 
           << std::endl;  
   
   // Add to error and norm
   norm+=(exact_soln[0]*exact_soln[0]+exact_soln[1]*exact_soln[1])*W;
   error+=( (exact_soln[0]-u_fe.real())*(exact_soln[0]-u_fe.real())+
            (exact_soln[1]-u_fe.imag())*(exact_soln[1]-u_fe.imag()) )*W;
   
  }
}




//======================================================================
 /// Compute norm of fe solution
//======================================================================
template <unsigned DIM>
void GeneralisedHelmholtzEquations<DIM>::compute_norm(double& norm)
{ 
 
 // Initialise
 norm=0.0;
 
 //Vector of local coordinates
 Vector<double> s(2);
 
 // Vector for coordintes
 Vector<double> x(2);
 
 //Find out how many nodes there are in the element
  unsigned n_node = nnode();
  
  Shape psi(n_node);
  
  //Set the value of n_intpt
  unsigned n_intpt = integral_pt()->nweight();
  
  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    
    //Assign values of s
    for(unsigned i=0;i<2;i++)
     {
      s[i] = integral_pt()->knot(ipt,i);
     }
    
    //Get the integral weight
    double w = integral_pt()->weight(ipt);
    
    // Get jacobian of mapping
    double J=J_eulerian(s);
    
    //Premultiply the weights and the Jacobian
    double W = w*J;
    
    // Get FE function value
    std::complex<double> u_fe=interpolated_u_generalised_helmholtz(s);
    
    // Add to  norm
    norm+=(u_fe.real()*u_fe.real()+u_fe.imag()*u_fe.imag())*W;
    
   }
 }
 
//====================================================================
// Force build of templates
//====================================================================
template class GeneralisedHelmholtzEquations<1>;
template class GeneralisedHelmholtzEquations<2>;
template class GeneralisedHelmholtzEquations<3>;

template class QGeneralisedHelmholtzElement<1,2>;
template class QGeneralisedHelmholtzElement<1,3>;
template class QGeneralisedHelmholtzElement<1,4>;

template class QGeneralisedHelmholtzElement<2,2>;
template class QGeneralisedHelmholtzElement<2,3>;
template class QGeneralisedHelmholtzElement<2,4>;

template class QGeneralisedHelmholtzElement<3,2>;
template class QGeneralisedHelmholtzElement<3,3>;
template class QGeneralisedHelmholtzElement<3,4>;

}
