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
#ifndef OOMPH_TRIANGLE_MESH_TEMPLATE_CC
#define OOMPH_TRIANGLE_MESH_TEMPLATE_CC

#include <iostream>

#include "triangle_mesh.template.h"
#include "../generic/map_matrix.h"
#include "../generic/multi_domain.h"
#include "../generic/projection.h"
#include "../generic/face_element_as_geometric_object.h"

namespace oomph
{

 //======================================================================
 /// Build with the help of the scaffold mesh coming
 /// from the triangle mesh generator Triangle.
 //======================================================================
 template<class ELEMENT>
  void
  TriangleMesh<ELEMENT>::build_from_scaffold(TimeStepper* time_stepper_pt,
    const bool &use_attributes)
 {
  // Mesh can only be built with 2D Telements.
  MeshChecker::assert_geometric_element<TElementGeometricBase,ELEMENT>(2);
  
  // Create space for elements
  unsigned nelem = Tmp_mesh_pt->nelement();
  Element_pt.resize(nelem);
  
  // Create space for nodes
  unsigned nnode_scaffold = Tmp_mesh_pt->nnode();
  
  // Create a map storing the node_id of the mesh used to update the
  // node position in the update_triangulateio function
  std::map<Node*, unsigned> old_global_number;
  
  // Store the TriangulateIO node id
  for (unsigned inod = 0; inod < nnode_scaffold; inod++)
   {
    Node* old_node_pt = Tmp_mesh_pt->node_pt(inod);
    old_global_number[old_node_pt] = inod;
   }
  
  // Initialize the old node id vector
  Oomph_vertex_nodes_id.resize(nnode_scaffold);
  
  // Create space for nodes
  Node_pt.resize(nnode_scaffold, 0);
  
  // Set number of boundaries
  unsigned nbound = Tmp_mesh_pt->nboundary();
  
  // Resize the boundary information
  set_nboundary(nbound);
  Boundary_element_pt.resize(nbound);
  Face_index_at_boundary.resize(nbound);
  
  //If we have different regions, then resize the region
  //information
  if (use_attributes)
   {
    Boundary_region_element_pt.resize(nbound);
    Face_index_region_at_boundary.resize(nbound);
   }

  // Loop over elements in scaffold mesh, visit their nodes
  for (unsigned e = 0; e < nelem; e++)
   {
    Element_pt[e] = new ELEMENT;
   }
  
  //Number of nodes per element from the scaffold mesh
  unsigned nnod_el = Tmp_mesh_pt->finite_element_pt(0)->nnode();
  
  // Setup map to check the (pseudo-)global node number
  // Nodes whose number is zero haven't been copied across
  // into the mesh yet.
  std::map<Node*, unsigned> global_number;
  unsigned global_count = 0;
  
  // Map of Element attribute pairs
  std::map<double, Vector<FiniteElement*> > element_attribute_map;
  
  // Loop over elements in scaffold mesh, visit their nodes
  for (unsigned e = 0; e < nelem; e++)
   {
    // Loop over all nodes in element
    for (unsigned j = 0; j < nnod_el; j++)
     {
      // Pointer to node in the scaffold mesh
      Node* scaffold_node_pt = Tmp_mesh_pt->finite_element_pt(e)->node_pt(j);

      // Get the (pseudo-)global node number in scaffold mesh
      // (It's zero [=default] if not visited this one yet)
      unsigned j_global = global_number[scaffold_node_pt];

      // Haven't done this one yet
      if (j_global == 0)
       {
        // Find and store the node_id in the old nodes map
        Oomph_vertex_nodes_id[global_count] =
         old_global_number[scaffold_node_pt];

        // Get pointer to set of mesh boundaries that this
        // scaffold node occupies; NULL if the node is not on any boundary
        std::set<unsigned>* boundaries_pt;
        scaffold_node_pt->get_boundaries_pt(boundaries_pt);

        //Storage for the new node
        Node* new_node_pt = 0;

        //Is it on boundaries
        if (boundaries_pt != 0)
         {
          //Create new boundary node
          new_node_pt = finite_element_pt(e)-> construct_boundary_node(j,
                                                                       time_stepper_pt);

          // Add to boundaries
          for (std::set<unsigned>::iterator it = boundaries_pt->begin(); it
                != boundaries_pt->end(); ++it)
           {
            add_boundary_node(*it, new_node_pt);
           }
         }
        //Build normal node
        else
         {
          //Create new normal node
          new_node_pt = finite_element_pt(e)-> construct_node(j,
                                                              time_stepper_pt);
         }

        // Give it a number (not necessarily the global node
        // number in the scaffold mesh -- we just need something
        // to keep track...)
        global_count++;
        global_number[scaffold_node_pt] = global_count;

        // Copy new node, created using the NEW element's construct_node
        // function into global storage, using the same global
        // node number that we've just associated with the
        // corresponding node in the scaffold mesh
        Node_pt[global_count - 1] = new_node_pt;

        // Assign coordinates
        for (unsigned i = 0; i < finite_element_pt(e)->dim(); i++)
         {
          new_node_pt->x(i) = scaffold_node_pt->x(i);
         }
       }
      // This one has already been done: Copy accross
      else
       {
        finite_element_pt(e)->node_pt(j) = Node_pt[j_global - 1];
       }
     }

    if (use_attributes)
     {
      element_attribute_map[Tmp_mesh_pt->element_attribute(e)].push_back(
       finite_element_pt(e));
     }
   }

  //Now let's construct lists
  //Find the number of attributes
  if (use_attributes)
   {
    unsigned n_attribute = element_attribute_map.size();
    //There are n_attribute different regions
    Region_attribute.resize(n_attribute);
    //Copy the vectors in the map over to our internal storage
    unsigned count = 0;
    for (std::map<double, Vector<FiniteElement*> >::iterator it =
          element_attribute_map.begin(); it != element_attribute_map.end(); ++it)
     {
      Region_attribute[count] = it->first;
      Region_element_pt[static_cast<unsigned>(Region_attribute[count])] = 
       it->second;
      ++count;
     }
   }

  // At this point we've created all the elements and
  // created their vertex nodes. Now we need to create
  // the additional (midside and internal) nodes!

  unsigned boundary_id;

  // Get number of nodes along element edge and dimension of element (2)
  // from first element
  unsigned n_node_1d = finite_element_pt(0)->nnode_1d();
  unsigned dim = finite_element_pt(0)->dim();

  // Storage for the local coordinate of the new node
  Vector<double> s(dim);

  // Get number of nodes in the element from first element
  unsigned n_node = finite_element_pt(0)->nnode();

  //Storage for each global edge of the mesh
  unsigned n_global_edge = Tmp_mesh_pt->nglobal_edge();
  Vector < Vector<Node*> > nodes_on_global_edge(n_global_edge);

  // Loop over elements
  for (unsigned e = 0; e < nelem; e++)
   {
    //Cache pointers to the elements
    FiniteElement* const elem_pt = finite_element_pt(e);
    FiniteElement* const tmp_elem_pt = Tmp_mesh_pt->finite_element_pt(e);

    //The number of edge nodes is  3*(nnode_1d-1)
    unsigned n_edge_node = 3 * (n_node_1d - 1);

    //If there are any more nodes, these are internal and can be
    //constructed and added directly to the mesh
    for (unsigned n = n_edge_node; n < n_node; ++n)
     {
      // Create new node (it can never be a boundary node)
      Node* new_node_pt = elem_pt->construct_node(n, time_stepper_pt);

      // What are the node's local coordinates?
      elem_pt->local_coordinate_of_node(n, s);

      // Find the coordinates of the new node from the existing
      // and fully-functional element in the scaffold mesh
      for (unsigned i = 0; i < dim; i++)
       {
        new_node_pt->x(i) = tmp_elem_pt->interpolated_x(s, i);
       }

      //Add the node to the mesh's global look-up scheme
      Node_pt.push_back(new_node_pt);
     }

    //Now loop over the mid-side edge nodes
    //Start from node number 3
    unsigned n = 3;

    // Loop over edges
    for (unsigned j = 0; j < 3; j++)
     {
      //Find the boundary id of the edge
      boundary_id = Tmp_mesh_pt->edge_boundary(e, j);

      //Find the global edge index
      unsigned edge_index = Tmp_mesh_pt->edge_index(e, j);

      //If the nodes on the edge have not been allocated, construct them
      if (nodes_on_global_edge[edge_index].size() == 0)
       {
        //Loop over the nodes on the edge excluding the ends
        for (unsigned j2 = 0; j2 < n_node_1d - 2; ++j2)
         {
          //Storage for the new node
          Node* new_node_pt = 0;

          //If the edge is on a boundary, construct a boundary node
          if (boundary_id > 0)
           {
            new_node_pt = elem_pt->construct_boundary_node(n, time_stepper_pt);
            //Add it to the boundary
            this->add_boundary_node(boundary_id - 1, new_node_pt);
           }
          //Otherwise construct a normal node
          else
           {
            new_node_pt = elem_pt->construct_node(n, time_stepper_pt);
           }

          // What are the node's local coordinates?
          elem_pt->local_coordinate_of_node(n, s);

          // Find the coordinates of the new node from the existing
          // and fully-functional element in the scaffold mesh
          for (unsigned i = 0; i < dim; i++)
           {
            new_node_pt->x(i) = tmp_elem_pt->interpolated_x(s, i);
           }

          //Add to the global node list
          Node_pt.push_back(new_node_pt);

          //Add to the edge index
          nodes_on_global_edge[edge_index].push_back(new_node_pt);
          //Increment the node number
          ++n;
         }
       }
      //Otherwise just set the pointers
      //using the fact that the next time the edge is visited
      //the nodes must be arranged in the other order because all
      //triangles have the same orientation
      else
       {
        //Loop over the nodes on the edge excluding the ends
        for (unsigned j2 = 0; j2 < n_node_1d - 2; ++j2)
         {
          //Set the local node from the edge but indexed the other
          //way around
          elem_pt->node_pt(n) = nodes_on_global_edge[edge_index][n_node_1d - 3
                                                                 - j2];
          ++n;
         }
       }

      //Set the elements adjacent to the boundary from the
      //boundary id information
      if (boundary_id > 0)
       {
        Boundary_element_pt[boundary_id - 1].push_back(elem_pt);
        //Need to put a shift in here because of an inconsistent naming
        //convention between triangle and face elements
        Face_index_at_boundary[boundary_id - 1].push_back((j + 2) % 3);

        //If using regions set up the boundary information
        if (use_attributes)
         {
          unsigned tmp_region =
           static_cast<unsigned> (Tmp_mesh_pt->element_attribute(e));
          //Element adjacent to boundary
          Boundary_region_element_pt[boundary_id - 1]
           [tmp_region].push_back(elem_pt);
          //Need to put a shift in here because of an inconsistent naming
          //convention between triangle and face elements
          Face_index_region_at_boundary[boundary_id - 1]
           [tmp_region].push_back((j + 2) % 3);
         }
       }

     } //end of loop over edges
   } //end of loop over elements

  // Lookup scheme has now been setup
  Lookup_for_elements_next_boundary_is_setup = true;

 }

 //======================================================================
 /// Setup boundary coordinate on boundary b. Doc Faces
 /// in outfile. Boundary coordinate increases continously along
 /// polygonal boundary. It's zero at the lexicographically
 /// smallest node on the boundary.
 //======================================================================
 template<class ELEMENT>
void TriangleMesh<ELEMENT>::setup_boundary_coordinates(const unsigned& b,
                                                       std::ofstream& outfile)
 {
  // Temporary storage for face elements
  Vector<FiniteElement*> face_el_pt;
  
  // Temporary storage for number of elements adjacent to the boundary
  unsigned nel = 0;
  
  // Temporary storage for elements adjacent to the boundary that
  // have an common edge
  unsigned n_repeated_ele = 0;
  
  unsigned n_regions = this->nregion();
  
  // Temporary storage for already done nodes
  Vector < std::pair<Node*, Node *> > done_nodes_pt;
  
  // Temporary created node pair
  std::pair<Node*, Node*> tmp_pair;
  std::pair<Node*, Node*> tmp_pair_inverse;
  
  // Re-used variable for storing the number of nodes of the
  // current element
  unsigned n_nodes;
  
  //If there is more than one region then only use
  //boundary coordinates from the bulk side (region 0)
  if (n_regions > 1)
   {
    for (unsigned rr = 0 ; rr < n_regions; rr++)
     {
      unsigned region_id = static_cast<unsigned>(this->Region_attribute[rr]);
      
      // Loop over all elements on boundaries in region i_r
      unsigned nel_in_region = this->nboundary_element_in_region(b, region_id);
      unsigned nel_repetead_in_region = 0;
      
// #ifdef PARANOID
//       if (nel_in_region==0)
//        {
//         std::ostringstream warning_message;
//         warning_message
//          << "There are no elements associated with boundary (" << b << ")\n"
//          << "in region (" << region_id << "). This could happen because:\n"
//          << "1) You did not specify boundaries with this boundary id.\n"
//          << "---- Review carefully the indexing of your boundaries.\n"
//          << "2) The boundary (" << b << ") is not associated with region ("
//          << region_id << ").\n"
//          << "---- The boundary does not touch the region\n.";
//         OomphLibWarning(warning_message.str(),
//                         "TriangleMesh::setup_boundary_coordinates()",
//                         OOMPH_EXCEPTION_LOCATION);
        
//        }
// #endif
      
      // Only bother to do anything else, if there are elements
      // associated with the boundary and the current region
      if (nel_in_region > 0)
       {
        
        // ***********************************************************
        // All the elements -- first inclusive -- (check for
        // repeated ones)
        bool repeated = false;
        
        // Loop over the bulk elements adjacent to boundary b
        for (unsigned e = 0; e < nel_in_region; e++)
         {
          // Get pointer to the bulk element that is adjacent to boundary b
          FiniteElement* bulk_elem_pt =
           this->boundary_element_in_region_pt(b, region_id, e);
          
          //Find the index of the face of element e along boundary b
          int face_index=this->face_index_at_boundary_in_region(b,region_id,e);
          
          // Before adding the new element we need to be sure that
          // the edge that this element represent has not been
          // already added
          FiniteElement* tmp_ele_pt = new DummyFaceElement<ELEMENT> (
           bulk_elem_pt, face_index);
          
          // Current number of added elements
          //unsigned c_ele_added = face_el_pt.size();
          
          n_nodes = tmp_ele_pt->nnode();
          
          tmp_pair = std::make_pair(tmp_ele_pt->node_pt(0),
                                    tmp_ele_pt->node_pt(n_nodes - 1));
          
          tmp_pair_inverse = std::make_pair(tmp_ele_pt->node_pt(n_nodes - 1),
                                            tmp_ele_pt->node_pt(0));
          
          // Search for repeated nodes
          unsigned repeated_nodes_size = done_nodes_pt.size();
          for (unsigned l = 0; l < repeated_nodes_size; l++)
           {
            if (tmp_pair == done_nodes_pt[l] || tmp_pair_inverse
                == done_nodes_pt[l])
             {
              nel_repetead_in_region++;
              repeated = true;
             }
           }
          
          done_nodes_pt.push_back(tmp_pair);
          
          // Create new face element
          if (!repeated)
           {
            face_el_pt.push_back(tmp_ele_pt);
           }
          else
           {
            // Clean up
             delete tmp_ele_pt;
           }
          
          // Re-start
          repeated = false;
          
          // Output faces?
          if (outfile.is_open())
           {
            face_el_pt[face_el_pt.size() - 1]->output(outfile);
           }
         } // for nel
        
        nel += nel_in_region;
        
        n_repeated_ele += nel_repetead_in_region;
        
       } // nel > 0
     } // for n_regions
   } // n_regions > 1
  //Otherwise it's just the normal boundary functions
  else
   {
    // Loop over all elements on boundaries
    nel = this->nboundary_element(b);
    
// #ifdef PARANOID
//     if (nel==0)
//      {
//       std::ostringstream warning_message;
//       warning_message
//        << "There are no elements associated with boundary (" << b << ").\n"
//        << "This could happen because you did not specify boundaries with\n"
//        << "this boundary id. Review carefully the indexing of your\n"
//        << "boundaries.";
//       OomphLibWarning(warning_message.str(),
//                       "TriangleMesh::setup_boundary_coordinates()",
//                       OOMPH_EXCEPTION_LOCATION);
//      }
// #endif
    
    //Only bother to do anything else, if there are elements
    if (nel > 0)
     {
      // Check for repeated ones
      bool repeated = false;
      
      // Loop over the bulk elements adjacent to boundary b
      for (unsigned e = 0; e < nel; e++)
       {
        // Get pointer to the bulk element that is adjacent to boundary b
        FiniteElement* bulk_elem_pt = this->boundary_element_pt(b, e);
        
        //Find the index of the face of element e along boundary b
        int face_index = this->face_index_at_boundary(b, e);
        
        // Before adding the new element we need to be sure that
        // the edge that this element represent has not been
        // already added
        FiniteElement* tmp_ele_pt = new DummyFaceElement<ELEMENT> (
         bulk_elem_pt, face_index);
        
        n_nodes = tmp_ele_pt->nnode();
        
        tmp_pair = std::make_pair(tmp_ele_pt->node_pt(0),
                                  tmp_ele_pt->node_pt(n_nodes - 1));
        
        tmp_pair_inverse = std::make_pair(tmp_ele_pt->node_pt(n_nodes - 1),
                                          tmp_ele_pt->node_pt(0));
        
        // Search for repeated nodes
        unsigned repeated_nodes_size = done_nodes_pt.size();
        for (unsigned l = 0; l < repeated_nodes_size; l++)
         {
          if (tmp_pair == done_nodes_pt[l] || tmp_pair_inverse
              == done_nodes_pt[l])
           {
            n_repeated_ele++;
            repeated = true;
           }
         }
        
        done_nodes_pt.push_back(tmp_pair);
        
        // Create new face element
        if (!repeated)
         {
          face_el_pt.push_back(tmp_ele_pt);
         }
        else
         {
          // Free the repeated bulk element!!
          delete tmp_ele_pt;
          tmp_ele_pt = 0;
         }
        
        // Re-start
        repeated = false;
        
        // Output faces?
        if (outfile.is_open())
         {
          face_el_pt[face_el_pt.size() - 1]->output(outfile);
         }
       }
     }
    
   }
  nel -= n_repeated_ele;
  
  //Only bother to do anything else, if there are elements
  if (nel > 0)
   {
    // Put first element into ordered list
    std::list<FiniteElement*> ordered_el_pt;
    FiniteElement* el_pt = face_el_pt[0];
    ordered_el_pt.push_back(el_pt);
    
    // Count elements that have been done
    unsigned count_done = 0;
    
    // Keep track of who's done
    std::map<FiniteElement*, bool> done_el;
    
    // Keep track of which element is inverted
    std::map<FiniteElement*, bool> is_inverted;
    
    // Number of nodes
    unsigned nnod = el_pt->nnode();
    
    // Fit in the other elements in at most nel^2 loops
    for (unsigned ee = 1; ee < nel; ee++)
     {
      // Loop over all elements to check if they fit to the right
      // or the left of the current one
      for (unsigned e = 1; e < nel; e++)
       {
        // Candidate element
        el_pt = face_el_pt[e];
        
        // Is it done yet?
        if (!done_el[el_pt])
         {
          
          // Left and rightmost elements
          FiniteElement* first_el_pt = (*ordered_el_pt.begin());
          std::list<FiniteElement*>::iterator it = ordered_el_pt.end();
          it--;
          FiniteElement* last_el_pt = *it;
          
          // Left and rightmost nodes
          Node* left_node_pt = first_el_pt->node_pt(0);
          if (is_inverted[first_el_pt])
           {
            left_node_pt = first_el_pt->node_pt(nnod - 1);
           }
          Node* right_node_pt = last_el_pt->node_pt(nnod - 1);
          if (is_inverted[last_el_pt])
           {
            right_node_pt = last_el_pt->node_pt(0);
           }
          
          // New element fits at the left of first element and is not inverted
          if (left_node_pt == el_pt->node_pt(nnod - 1))
           {
            ordered_el_pt.push_front(el_pt);
            done_el[el_pt] = true;
            count_done++;
            is_inverted[el_pt] = false;
           }
          // New element fits at the left of first element and is inverted
          else if (left_node_pt == el_pt->node_pt(0))
           {
            ordered_el_pt.push_front(el_pt);
            done_el[el_pt] = true;
            count_done++;
            is_inverted[el_pt] = true;
           }
          // New element fits on the right of last element and is not inverted
          else if (right_node_pt == el_pt->node_pt(0))
           {
            ordered_el_pt.push_back(el_pt);
            done_el[el_pt] = true;
            count_done++;
            is_inverted[el_pt] = false;
           }
          // New element fits on the right of last element and is inverted
          else if (right_node_pt == el_pt->node_pt(nnod - 1))
           {
            ordered_el_pt.push_back(el_pt);
            done_el[el_pt] = true;
            count_done++;
            is_inverted[el_pt] = true;
           }
          
          if (done_el[el_pt])
           {
            break;
           }
         }
       }
     }
    
    // Are we done?
    if (count_done != (nel - 1))
     {
      std::ostringstream error_message;
      error_message 
       << "Was only able to setup boundary coordinate on "
       << "boundary " << b << "\nfor " << count_done << " of " << nel
       << " face elements. This usually means\n"
       << "that the boundary is not simply connected.\n\n"
       << "Re-run the setup_boundary_coordinates() function\n"
       << "with an output file specified " << "as the second argument.\n"
       << "This will file will contain FaceElements that\n"
       << "oomph-lib believes to be located on the boundary.\n" << std::endl;
      throw OomphLibError(
       error_message.str(),
       OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
     }
    
    // First node
    FiniteElement* first_el_pt = *ordered_el_pt.begin();
    Node* first_node_pt = first_el_pt->node_pt(0);
    if (is_inverted[first_el_pt])
     {
      first_node_pt = first_el_pt->node_pt(nnod - 1);
     }
    
    //Last node
    FiniteElement* last_el_pt = ordered_el_pt.back();
    Node* last_node_pt = last_el_pt->node_pt(nnod - 1);
    if (is_inverted[last_el_pt])
     {
      last_node_pt = last_el_pt->node_pt(0);
     }
    
    // Coordinates of left node
    double x_left = first_node_pt->x(0);
    double y_left = first_node_pt->x(1);
    
    // Initialise boundary coordinate
    Vector<double> zeta(1, 0.0);
    
    // Set boundary coordinate
    first_node_pt->set_coordinates_on_boundary(b, zeta);
    
    // Lexicographically bottom left node
    std::set<Node*> all_nodes_pt;
    all_nodes_pt.insert(first_node_pt);
    
    // Now loop over nodes in order
    for (std::list<FiniteElement*>::iterator it = ordered_el_pt.begin(); it
          != ordered_el_pt.end(); it++)
     {
      // Get element
      FiniteElement* el_pt = *it;
      
      // Start node and increment
      unsigned k_nod = 1;
      int nod_diff = 1;
      if (is_inverted[el_pt])
       {
        k_nod = nnod - 2;
        nod_diff = -1;
       }
      
      // Loop over nodes
      for (unsigned j = 1; j < nnod; j++)
       {
        Node* nod_pt = el_pt->node_pt(k_nod);
        k_nod += nod_diff;
        
        // Coordinates of right node
        double x_right = nod_pt->x(0);
        double y_right = nod_pt->x(1);
        
        // Increment boundary coordinate
        zeta[0] += sqrt(
         (x_right - x_left) * (x_right - x_left) + (y_right - y_left)
         * (y_right - y_left));
        
        // Set boundary coordinate
        nod_pt->set_coordinates_on_boundary(b, zeta);
        
        // Increment reference coordinate
        x_left = x_right;
        y_left = y_right;
        
        // Get lexicographically bottom left node but only
        // use vertex nodes as candidates
        all_nodes_pt.insert(nod_pt);
       }
     }
    
    //The nodes have been assigned arc-length coordinates
    //from one end or the other of the connected segment.
    
    //If the boundary has a geometric object representation then
    //scale the coordinates to match those of the geometric object
    GeomObject* const geom_object_pt = this->boundary_geom_object_pt(b);
    if (geom_object_pt != 0)
     {
      Vector<double> bound_coord_limits = this->boundary_coordinate_limits(b);
      
      //Get the position of the ends of the geometric object
      Vector<double> zeta(1);
      Vector<double> first_geom_object_location(2);
      Vector<double> last_geom_object_location(2);
      zeta[0] = bound_coord_limits[0];
      geom_object_pt->position(zeta, first_geom_object_location);
      zeta[0] = bound_coord_limits[1];
      geom_object_pt->position(zeta, last_geom_object_location);
      
      //Calculate the errors in position between the first and last nodes
      //and the endpoints of the geometric object
      double error = 0.0;
      double tmp_error = 0.0;
      for (unsigned i = 0; i < 2; i++)
       {
        const double dist = first_geom_object_location[i]
         - first_node_pt->x(i);
        tmp_error += dist * dist;
       }
      error += sqrt(tmp_error);
      tmp_error = 0.0;
      for (unsigned i = 0; i < 2; i++)
       {
        const double dist =
         last_geom_object_location[i] - last_node_pt->x(i);
        tmp_error += dist * dist;
       }
      error += sqrt(tmp_error);
      
      //Calculate the errors in position between the first and last nodes
      //and the endpoints of the geometric object if reversed
      double rev_error = 0.0;
      tmp_error = 0.0;
      for (unsigned i = 0; i < 2; i++)
       {
        const double dist =
         first_geom_object_location[i] - last_node_pt->x(i);
        tmp_error += dist * dist;
       }
      rev_error += sqrt(tmp_error);
      tmp_error = 0.0;
      for (unsigned i = 0; i < 2; i++)
       {
        const double dist =
         last_geom_object_location[i] - first_node_pt->x(i);
        tmp_error += dist * dist;
       }
      rev_error += sqrt(tmp_error);
      
      // Number of polyline vertices along this boundary
      unsigned n_vertex = Polygonal_vertex_arclength_info[b].size();
      
      // Get polygonal vertex data
      Vector < std::pair<double, double> > polygonal_vertex_arclength
       = Polygonal_vertex_arclength_info[b];
      
      //If the (normal) error is small than reversed then we have the
      //coordinate direction correct.
      //If not then we must reverse it
      //bool reversed = false;
      if (error < rev_error)
       {
        // Coordinates are aligned (btw: don't delete this block -- there's
        // a final else below to catch errors!)
        //reversed = false;
       }
      else if (error > rev_error)
       {
        //reversed = true;
        
        //Reverse the limits of the boundary coordinates along the
        //geometric object
        double temp = bound_coord_limits[0];
        bound_coord_limits[0] = bound_coord_limits[1];
        bound_coord_limits[1] = temp;
        for (unsigned v = 0; v < n_vertex; v++)
         {
          polygonal_vertex_arclength[v].first
           = Polygonal_vertex_arclength_info[b][v].first;
          
          polygonal_vertex_arclength[v].second
           = Polygonal_vertex_arclength_info[b][n_vertex - v - 1].second;
         }
       }
      else
       {
        std::ostringstream error_stream;
        error_stream 
         << "Something very strange has happened.\n"
         << "The error between the endpoints of the geometric object\n"
         << "and the first and last nodes on the boundary is the same\n"
         << "irrespective of the direction of the coordinate.\n"
         << "This probably means that things are way soff.\n"
         << "The errors are " << error << " and " << rev_error << "\n";
        throw OomphLibError(error_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
       }
      
      //Get the total arclength of the edge
      last_node_pt->get_coordinates_on_boundary(b, zeta);
      
      double zeta_old_range = zeta[0];
      double zeta_new_range = bound_coord_limits[1] - bound_coord_limits[0];
      
      // Re-assign boundary coordinate for the case where boundary
      // is represented by polygon
      bool use_old = false;
      if (n_vertex==0) use_old = true;
      
      //Now scale the coordinates accordingly
      for (std::set<Node*>::iterator it = all_nodes_pt.begin(); it
            != all_nodes_pt.end(); it++)
       {
        Node* nod_pt = (*it);
        
        // Get coordinate based on arclength along polygonal repesentation
        nod_pt->get_coordinates_on_boundary(b, zeta);
                
        if (use_old)
         {
          // Boundary is actually a polygon -- simply rescale
          zeta[0] = bound_coord_limits[0] + (zeta_new_range / zeta_old_range)
           * zeta[0];
         }
        else
         {
          // Scale such that vertex nodes stay where they were
          bool found = false;
          
          // Loop over vertex nodes
          for (unsigned v = 1; v < n_vertex; v++)
           {
            if ((zeta[0] >= polygonal_vertex_arclength[v - 1].first)
                && (zeta[0] <= polygonal_vertex_arclength[v].first))
             {              
              // Increment in intrinsic coordinate along geom object
              double delta_zeta = (polygonal_vertex_arclength[v].second
                                   -polygonal_vertex_arclength[v - 1].second);
              // Increment in arclength along segment
              double delta_polyarc = (polygonal_vertex_arclength[v].first
                                      -polygonal_vertex_arclength[v - 1].first);
              
              // Mapped arclength coordinate
              double zeta_new = polygonal_vertex_arclength[v - 1].second
               + delta_zeta * (zeta[0]
                               - polygonal_vertex_arclength[v - 1].first) / 
               delta_polyarc;
              zeta[0] = zeta_new;
              
              // Success!
              found = true;
              
              // Bail out
              break;
             }
           }
          
          // If we still haven't found it's probably the last point along
          if (!found)
           {
#ifdef PARANOID
            double diff=std::fabs(zeta[0]-
                                  polygonal_vertex_arclength[n_vertex-1].first);
            if (diff>ToleranceForVertexMismatchInPolygons::Tolerable_error)
             {
              std::ostringstream error_stream;
              error_stream
               << "Wasn't able to locate the polygonal arclength exactly\n"
               << "during re-setup of boundary coordinates and have\n"
               << "assumed that we're dealing with the final point along\n"
               << "the curvilinear segment and encountered some roundoff\n"
               << "However, the difference in the polygonal zeta coordinates\n"
               << "between zeta[0] " << zeta[0]
               << " and the originally stored value "
               << polygonal_vertex_arclength[n_vertex-1].first << "\n"
               << "is " << diff << " which exceeds the threshold specified\n"
               << "in the publically modifiable variable\n"
               << "ToleranceForVertexMismatchInPolygons::Tolerable_error\n"
               << "whose current value is: "
               << ToleranceForVertexMismatchInPolygons::Tolerable_error
               << "\nPlease check your mesh carefully and increase the\n"
               << "threshold if you're sure this is appropriate\n";
              throw OomphLibError(error_stream.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
             }
#endif
            zeta[0] = polygonal_vertex_arclength[n_vertex - 1].second;
            }
         }
        
        // Assign updated coordinate
        nod_pt->set_coordinates_on_boundary(b, zeta);
       }
     }
    else
     {
      //Only use end points of the whole segment and pick the bottom left
      //node
      Node* bottom_left_node_pt = first_node_pt;
      if (last_node_pt->x(1) < bottom_left_node_pt->x(1))
       {
        bottom_left_node_pt = last_node_pt;
       }
      else if (last_node_pt->x(1) == bottom_left_node_pt->x(1))
       {
        if (last_node_pt->x(0) < bottom_left_node_pt->x(0))
         {
          bottom_left_node_pt = last_node_pt;
         }
       }
      
      // Now adjust boundary coordinate so that the bottom left node
      // has a boundary coordinate of zero and that zeta increases
      // away from that point
      bottom_left_node_pt->get_coordinates_on_boundary(b, zeta);
      double zeta_ref = zeta[0];
      double zeta_max = 0.0;
      for (std::set<Node*>::iterator it = all_nodes_pt.begin(); it
            != all_nodes_pt.end(); it++)
       {
        Node* nod_pt = (*it);
        nod_pt->get_coordinates_on_boundary(b, zeta);
        zeta[0] -= zeta_ref;
        //If direction is reversed, then take absolute value
        if (zeta[0] < 0.0)
         {
          zeta[0] = std::fabs(zeta[0]);
         }
        if (zeta[0] > zeta_max)
         {
          zeta_max = zeta[0];
         }
        nod_pt->set_coordinates_on_boundary(b, zeta);
       }
      
      //Scale all surface coordinates by the max
      for (std::set<Node*>::iterator it = all_nodes_pt.begin(); it
            != all_nodes_pt.end(); it++)
       {
        Node* nod_pt = (*it);
        nod_pt->get_coordinates_on_boundary(b, zeta);
        zeta[0] /= zeta_max;
        nod_pt->set_coordinates_on_boundary(b, zeta);
       }
     }
    
    // Cleanup
    for (unsigned e = 0; e < nel; e++)
     {
      delete face_el_pt[e];
      face_el_pt[e] = 0;
     }
    
    //Nodes_on_boundary_pt[b] = all_nodes_pt;
    
   }
  
  // Indicate that boundary coordinate has been set up
  Boundary_coordinate_exists[b] = true;
  
 }
 
#ifdef OOMPH_HAS_TRIANGLE_LIB

 //======================================================================
 /// Create TriangulateIO object from TriangleMeshPolygon (outer_boundary)
 /// TriangleMeshPolygon (holes) and TriangleMeshPolyLines (curves in the
 /// domain)
 //======================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::build_triangulateio(
   TriangleMeshPolygon* &outer_boundary_pt,
   Vector<TriangleMeshPolygon*> &internal_polygon_pt,
   Vector<TriangleMeshOpenCurve*> &internal_polylines_pt,
   Vector<Vector<double> > &extra_holes_coordinates,
   std::map<unsigned, Vector<double> > &regions_coordinates,
   std::map<unsigned, double> &regions_areas,
   TriangulateIO& triangulate_io)
   {
  // triangulate_io initialization
  TriangleHelper::initialise_triangulateio(triangulate_io);

  // Count the global number of boundary polylines (used for
  // sanity check)
  unsigned n_boundglobalseg=0;

  // Total number of vertices and segments for outer boundary,
  // internal closed boundaries and internal open boundaries
  unsigned n_outer_vertices = 0;
  unsigned n_internal_closed_vertices = 0;
  unsigned n_internal_open_vertices = 0;

  unsigned n_outer_segments = 0;
  unsigned n_internal_closed_segments = 0;
  unsigned n_internal_open_segments = 0;

  // Total number of vertices and segments
  unsigned n_global_vertices=0;
  unsigned n_global_segments=0;

  // General purpose counters
  unsigned count_tri=0;
  int edge_segment=1;

  // Reusable variable
  unsigned n_polyline_vertices=0;

  // We need to know the dimension of the global number of vertices
  // and segments to initialise the pointlist, segmentlist and
  // segmentmarkerlist.

  // **************************************************************
  // Build the TriangulateIO object (1st stage)
  // Counting vertices and segments in the outer boundary(1.1)
  // **************************************************************

  // Get number of polyline
  unsigned n_outer_boundary_polylines=0;
  n_outer_boundary_polylines = outer_boundary_pt->npolyline();

  // Loop over the boundary polylines
  for(unsigned count_seg=0;count_seg<n_outer_boundary_polylines;count_seg++)
   {
    n_outer_vertices +=
      outer_boundary_pt->curve_section_pt(count_seg)->nvertex()-1;
    n_outer_segments +=
      outer_boundary_pt->curve_section_pt(count_seg)->nsegment();
   }

  // **************************************************************
  // Build the TriangulateIO object (1st stage)
  // Counting vertices and segments in the internal closed
  // boundaries (1.2)
  // **************************************************************

  // Count the global closed boundaries boundary polyline
  unsigned n_holes = internal_polygon_pt.size();

  // Counter initialisation
  unsigned n_holepolyline=0;

  // Loop over the holes boundary
  for(unsigned count_hole=0;count_hole<n_holes;count_hole++)
   {
    n_holepolyline=internal_polygon_pt[count_hole]->npolyline();
    for(unsigned count_seg=0;count_seg<n_holepolyline;count_seg++)
     {
      n_internal_closed_vertices +=
        internal_polygon_pt[count_hole]->
        polyline_pt(count_seg)->nvertex()-1;
      n_internal_closed_segments +=
        internal_polygon_pt[count_hole]->
        polyline_pt(count_seg)->nsegment();
     }
   }

  // **************************************************************
  // Build the TriangulateIO object (1st stage)
  // Counting vertices and segments in the outer boundary (1.3)
  // Finalising the counting of outer boundary vertices
  // **************************************************************

  // If there's just one boundary. All the vertices should be counted
  // for the polygon...
  if(n_outer_boundary_polylines==1)
   {
    n_outer_vertices+=1;
   }

  // **************************************************************
  // Build the TriangulateIO object (1st stage)
  // Counting vertices and segments in the internal closed
  // boundaries (1.4)
  // Finalising the counting of internal closed boundaries vertices
  // **************************************************************

  // ...and for the hole.
  for(unsigned count_hole=0;count_hole<n_holes;count_hole++)
   {
    n_holepolyline=internal_polygon_pt[count_hole]->npolyline();
    if(n_holepolyline==1)
     {
      n_internal_closed_vertices+=1;
     }
   }

  // **************************************************************
  // Build the TriangulateIO object (1st stage)
  // Counting vertices and segments in the internal open boundaries
  // (1.5)
  // **************************************************************

  // Number of open polylines
  unsigned n_open_polylines = internal_polylines_pt.size();

  // Loop over the open polylines
  for (unsigned b =0; b < n_open_polylines; b++)
   {
    // Number of polylines on the current open polyline
    unsigned n_polylines = internal_polylines_pt[b]->ncurve_section();
    // Loop over the constituent polylines
    for(unsigned p = 0; p < n_polylines; p++)
     {
      // We are not double counting internal vertices
      // We have a special case if there is just one polyline per open curve
      if (p == 0)
	{
	  n_internal_open_vertices+=
	    internal_polylines_pt[b]->curve_section_pt(p)->nvertex();
	}
      else
	{
	  n_internal_open_vertices+=
	    internal_polylines_pt[b]->curve_section_pt(p)->nvertex()-1;
	}
	
      // When two lines are connected (via a vertex) we need to avoid
      // counting the connection vertex twice
      // If the polyline is connected to an outer or another internal
      // boundary then we need to avoid double counting it
      // We do not have to store the vertices coordinates twice
      n_internal_open_vertices-=
        internal_polylines_pt[b]->
        polyline_pt(p)->is_initial_vertex_connected();

      n_internal_open_vertices-=
        internal_polylines_pt[b]->
        polyline_pt(p)->is_final_vertex_connected();
      
      n_internal_open_segments+=
        internal_polylines_pt[b]->curve_section_pt(p)->nsegment();

     }
    
   }

#ifdef PARANOID
  if (n_outer_vertices!=n_outer_segments)
   {
    std::ostringstream error_message;
    error_message
    << "The number of vertices and segments in the\n"
    << "outer boundary is different. The number of\n"
    << "vertices is " << n_outer_vertices << " and\n"
    << "the number of segments is " << n_outer_segments
    << std::endl;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
   }

  if (n_internal_closed_vertices!=n_internal_closed_segments)
   {
    std::ostringstream error_message;
    error_message
    << "The number of vertices and segments in the\n"
    << "internal closed boundaries is different. The "
    << "number of vertices is " << n_internal_closed_vertices
    << " and the number of segments is "
    << n_internal_closed_segments << std::endl;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
   }
#endif

  n_global_vertices =
    n_outer_vertices +
    n_internal_closed_vertices +
    n_internal_open_vertices;
  n_global_segments =
    n_outer_segments +
    n_internal_closed_segments +
    n_internal_open_segments;

  // Store the global number of vertices and segments
  // in the list.
  triangulate_io.numberofpoints = n_global_vertices;
  triangulate_io.numberofsegments = n_global_segments;

  // Preparing the TriangulateIO objects to store the values
  triangulate_io.pointlist =
    (double *) malloc(triangulate_io.numberofpoints * 2 * sizeof(double));
  triangulate_io.segmentlist =
    (int *) malloc(triangulate_io.numberofsegments * 2 * sizeof(int));
  triangulate_io.segmentmarkerlist =
    (int *) malloc(triangulate_io.numberofsegments * sizeof(int));

  // Map for storing connections information (used for recovering the
  // 'edge_segment' number when performing the connections)
  // The map uses the boundary id and the vertex number as a 'key-pair'
  // for storing the 'edge_segment' when it is connected
  std::map<std::pair<unsigned, unsigned>, unsigned> vertices_connections;

  // Temporal storage when using the 'key' on the above map
  std::pair<unsigned, unsigned> pair_id_vertex;

  // **************************************************************
  // Store values in the TriangulateIO object (2nd stage)
  // **************************************************************

  // **************************************************************
  // Store values in the TriangulateIO object (2nd stage)
  // Outer boundary first (2.1)
  // **************************************************************
  // Storing all the values in the list
  for(unsigned count_seg=0;count_seg<n_outer_boundary_polylines;count_seg++)
   {
    // Storing the number of the vertices
    n_polyline_vertices =
      outer_boundary_pt->curve_section_pt(count_seg)->nvertex()-1;

    // If there's just one boundary. All the vertices have to be counted
    if(n_outer_boundary_polylines==1)
     {
      n_polyline_vertices+=1;
     }

    // Store the segmen_id if given by the user
    unsigned idpolyline=outer_boundary_pt->curve_section_pt(count_seg)
       ->boundary_id();

    //Storing the coordinates for each points
    for(unsigned count_vertices=0;count_vertices<n_polyline_vertices;
      count_vertices++)
     {
      triangulate_io.pointlist[count_tri]=
        outer_boundary_pt->polyline_pt(count_seg)->
        vertex_coordinate(count_vertices)[0];
      triangulate_io.pointlist[count_tri+1]= outer_boundary_pt->
        polyline_pt(count_seg)->vertex_coordinate(count_vertices)[1];

      // Store the segment values
      // If the segment is not the last one, take the next node
      if(count_seg==(n_outer_boundary_polylines-1) &&
        count_vertices==(n_polyline_vertices-1))
       {
        // Adding the information to the vertices_connections map
        // The vertex "count_vertices" of the boundary with
        // "boundary_id" = idpolyline is stored on the "edge_segment"
        pair_id_vertex = std::make_pair(idpolyline, count_vertices);
        vertices_connections[pair_id_vertex] = edge_segment;
        
        // Added line to allow connect to a vertex by using any reference
        // of it (as last vertex of previous polyline or as first vertex
        // of the current polyline)
        if (count_seg>0 && count_vertices==0)
         {
          unsigned id_previous_polyline = 
           outer_boundary_pt->curve_section_pt(count_seg-1)->boundary_id();

          unsigned n_previous_polyline_vertices =
           outer_boundary_pt->curve_section_pt(count_seg-1)->nvertex()-1;

          pair_id_vertex = std::make_pair(
           id_previous_polyline, n_previous_polyline_vertices);
          vertices_connections[pair_id_vertex] = edge_segment;    

         }
        
        // By doing this we close the boundary
        // By setting that the segment is composed by the
        // vertex with number "edge_segment" and the vertex
        // with number 1 it means that this is closing the boundary
        triangulate_io.segmentlist[count_tri]=edge_segment;
        triangulate_io.segmentlist[count_tri+1]=1;
       }
      else
       {
        // Adding the information to the vertices_connections map
        // The vertex "count_vertices" of the boundary with
        // "boundary_id" = idpolyline is stored on the "edge_segment"
        pair_id_vertex = std::make_pair(idpolyline, count_vertices);
        vertices_connections[pair_id_vertex] = edge_segment;

        // Added line to allow connect to a vertex by using any reference
        // of it (as last vertex of previous polyline or as first vertex
        // of the current polyline)
        if (count_seg>0 && count_vertices==0)
         {
          unsigned id_previous_polyline = 
           outer_boundary_pt->curve_section_pt(count_seg-1)->boundary_id();

          unsigned n_previous_polyline_vertices =
           outer_boundary_pt->curve_section_pt(count_seg-1)->nvertex()-1;

          pair_id_vertex = std::make_pair(
           id_previous_polyline, n_previous_polyline_vertices);
          vertices_connections[pair_id_vertex] = edge_segment;    

         }

        triangulate_io.segmentlist[count_tri]=edge_segment;
        triangulate_io.segmentlist[count_tri+1]=edge_segment+1;
        edge_segment++;
       }

      // Store the marker list of the segments
      triangulate_io.segmentmarkerlist[count_tri/2]=idpolyline+1;

      // Increment counter
      count_tri+=2;
      n_boundglobalseg++;
     }
   }

  // **************************************************************
  // Store values in the TriangulateIO object (2nd stage)
  // Internal closed boundaries (2.2)
  // **************************************************************

  // Add the hole's vertices coordinates and segments if provided
  // Storing all the values in the list
  for(unsigned count_hole=0;count_hole<n_holes;count_hole++)
   {
    // Store the value of the first starting counting node for the hole
    // The hole's boundary doesn't share node with the external polygon
    edge_segment++;
    int hole_vertex_start = edge_segment;
    n_holepolyline=internal_polygon_pt[count_hole]->npolyline();

    for(unsigned count_seg=0;count_seg<n_holepolyline;count_seg++)
     {
      // Storing the number of the vertices
      n_polyline_vertices = internal_polygon_pt[count_hole]->
        polyline_pt(count_seg)->nvertex()-1;

      // If there's just one boundary. All the vertices should be counted
      if(n_holepolyline==1)
       {
        n_polyline_vertices +=1;
       }

      // Store the segment_id if given by the user
      unsigned idpolyline=internal_polygon_pt[count_hole]->
        polyline_pt(count_seg)->boundary_id();

      // Store the coordinates for each points
      for(unsigned count_vertices=0;count_vertices<n_polyline_vertices;
        count_vertices++)
       {
        triangulate_io.pointlist[count_tri]=
          internal_polygon_pt[count_hole]->polyline_pt(count_seg)->
          vertex_coordinate(count_vertices)[0];
        triangulate_io.pointlist[count_tri+1]=
          internal_polygon_pt[count_hole]->
          polyline_pt(count_seg)->
          vertex_coordinate(count_vertices)[1];

        // Store the segments values
        // If the segment is not the last one, take the next node
        if(count_seg==(n_holepolyline-1) &&
          count_vertices==(n_polyline_vertices-1))
         {
          // Adding the information to the vertices_connections map
          // The vertex "count_vertices" of the boundary with
          // "boundary_id" = idpolyline is stored on the "edge_segment"
          pair_id_vertex = std::make_pair(idpolyline, count_vertices);
          vertices_connections[pair_id_vertex] = edge_segment;
          
          // Added line to allow connect to a vertex by using any reference
          // of it (as last vertex of previous polyline or as first vertex
          // o f the current polyline)
          if (count_seg>0 && count_vertices==0)
           {
            unsigned id_previous_polyline = 
             internal_polygon_pt[count_hole]->curve_section_pt(count_seg-1)->
             boundary_id();
            
            unsigned n_previous_polyline_vertices =
             internal_polygon_pt[count_hole]->curve_section_pt(count_seg-1)->
             nvertex()-1;
            
            pair_id_vertex = std::make_pair(
             id_previous_polyline, n_previous_polyline_vertices);
            vertices_connections[pair_id_vertex] = edge_segment;    
            
           }
          
          // Closing the boundary by setting the last vertex
          // same as the first vertex of the internal closed boundary
          // (hole_vertex_start)
          triangulate_io.segmentlist[count_tri]=edge_segment;
          triangulate_io.segmentlist[count_tri+1]=hole_vertex_start;
         }
        else
         {
          // Adding the information to the vertices_connections map
          // The vertex "count_vertices" of the boundary with
          // "boundary_id" = idpolyline is stored on the "edge_segment"
          pair_id_vertex = std::make_pair(idpolyline, count_vertices);
          vertices_connections[pair_id_vertex] = edge_segment;

          // Added line to allow connect to a vertex by using any reference
          // of it (as last vertex of previous polyline or as first vertex
          // o f the current polyline)
          if (count_seg>0 && count_vertices==0)
           {
            unsigned id_previous_polyline = 
             internal_polygon_pt[count_hole]->curve_section_pt(count_seg-1)->
             boundary_id();
            
            unsigned n_previous_polyline_vertices =
             internal_polygon_pt[count_hole]->curve_section_pt(count_seg-1)->
             nvertex()-1;
            
            pair_id_vertex = std::make_pair(
             id_previous_polyline, n_previous_polyline_vertices);
            vertices_connections[pair_id_vertex] = edge_segment;    
            
           }
          
          triangulate_io.segmentlist[count_tri]=edge_segment;
          triangulate_io.segmentlist[count_tri+1]=edge_segment+1;
          edge_segment++;
         }

        // Store the marker list of the segments
        triangulate_io.segmentmarkerlist[count_tri/2]= idpolyline+1;

        count_tri+=2;
        n_boundglobalseg++;
       }
     }
   }

  // **************************************************************
  // Store values in the TriangulateIO object (2nd stage)
  // Internal open boundaries (2.3)
  // **************************************************************

  // For managing the counting of the number of vertices and
  // segments separately
  unsigned count_tri_seg = count_tri;

  // Index for segment list
  unsigned seg_lst_idx = count_tri_seg;

  // Increment edge segment counter
  edge_segment++;

  for (unsigned b = 0; b < n_open_polylines; b++)
   {
    // Number of constituent curve section on the open curve
    unsigned n_polyline = internal_polylines_pt[b]->ncurve_section();

    // Temporal storage for current open curve
    TriangleMeshOpenCurve *open_polyline_pt =
      internal_polylines_pt[b];

    // Loop over all the polylines on the open curve
    for(unsigned count_seg=0;count_seg<n_polyline;count_seg++)
     {
      // Temporal storage for the current polyline
      TriangleMeshPolyLine *polyline_pt =
        open_polyline_pt->polyline_pt(count_seg);

      // Storing the number of the vertices
      n_polyline_vertices = polyline_pt->nvertex();

      // Store the boundary id if given by the user
      unsigned idpolyline = polyline_pt->boundary_id();

      // Do the initial and final vertex separately
      // *********************************************************************
      // Initial vertex on first polyline
      // *********************************************************************
      if (count_seg == 0 && polyline_pt->is_initial_vertex_connected())
       {
        // If it is connected then we already know its coordinates,
        // therefore we do not need to add them to the list but we
        // need to search for the corresponding 'edge_segment' that
        // stores them

        // For this purpose we need the 'boudary_id' of the boundary
        // to which the current boundary is connected and the
        // associated vertex number
        unsigned connected_id_polyline =
          polyline_pt->initial_vertex_connected_bnd_id();

        unsigned connected_n_vertex =
          polyline_pt->initial_vertex_connected_n_vertex();

        // Create the pair
        pair_id_vertex =
          std::make_pair(connected_id_polyline, connected_n_vertex);

        // Here is where we "recover" the 'edge_segment'
        unsigned recovered_edge_segment =
          vertices_connections[pair_id_vertex];

        // Adding the information to the vertices_connections map.
        // The vertex "0" of "boundary_id" = idpolyline is connected
        // to the "recovered_edge_segment"
        pair_id_vertex = std::make_pair(idpolyline, 0);
        vertices_connections[pair_id_vertex] = recovered_edge_segment;

        // Here is where we use the corresponding "recovered_edge_segment"
        // number and create the connection in the "triangulate_io"
        // structure
        triangulate_io.segmentlist[seg_lst_idx]=recovered_edge_segment;

        // Increment the segment list index
        seg_lst_idx++;
       }
      else if (count_seg == 0)
       {
        // If it is not connected then we do not know its
        // coordinates and then we need to add them to the list
        triangulate_io.pointlist[count_tri]=
          polyline_pt->vertex_coordinate(0)[0];
        triangulate_io.pointlist[count_tri+1]=
          polyline_pt->vertex_coordinate(0)[1];

        // Increment counter
        count_tri+=2;

        // Adding the information to the vertices_connections map.
        // The vertex "0" of "boundary_id" = idpolyline is connected
        // to the "edge_segment"
        pair_id_vertex = std::make_pair(idpolyline, 0);
        vertices_connections[pair_id_vertex] = edge_segment;

        // It is not connected to the edges so there is no need to look
        // for the corresponding edge_segment
        triangulate_io.segmentlist[seg_lst_idx] = edge_segment;

        // Increment the edge segment counter since we have used it
        edge_segment++;

        // Increment the segment list index
        seg_lst_idx++;

       }

      // *********************************************************************
      // Medium vertices (no the very first and last one)
      // *********************************************************************
      // Variables for having control of the 'for' loop limits.
      // 1) When dealing with the first polyline we have already processed
      //    the first vertex therefore we should start on the next vertex
      // 2) When dealing with the last polyline we should process the last
      //    vertex out of the 'for' loop
      unsigned lower_limit_index = 0;
      unsigned upper_limit_index = n_polyline_vertices;

      // The first vertex is always treated by the previous polyline, therefore
      // we always start on the second vertex
      lower_limit_index = 1;

      // If we are dealing with the last segment
      if (count_seg == n_polyline-1)
       {upper_limit_index = n_polyline_vertices-1;}

      // Storing the coordinates for each point in this polyline
      for(unsigned count_vertices=lower_limit_index;
        count_vertices < upper_limit_index; count_vertices++)
       {
        // It is not connected then we do not know its coordinates.
        // We need to add them to the list
        triangulate_io.pointlist[count_tri]=
          polyline_pt->vertex_coordinate(count_vertices)[0];
        triangulate_io.pointlist[count_tri+1]=
          polyline_pt->vertex_coordinate(count_vertices)[1];

        // Increment counter
        count_tri+=2;

        // Adding the information to the vertices_connections map.
        // The vertex "count_vertices" of "boudary_id" = idpolyline
        // is connected to the "edge_segment"
        pair_id_vertex = std::make_pair(idpolyline, count_vertices);
        vertices_connections[pair_id_vertex] = edge_segment;

        // Added line to allow connect to a vertex by using any reference
        // of it (as last vertex of previous polyline or as first vertex
        // o f the current polyline)
        if (count_seg<n_polyline-1&&count_vertices==upper_limit_index-1)
         {
          unsigned id_next_polyline = 
           open_polyline_pt->polyline_pt(count_seg+1)->boundary_id();

          unsigned n_next_polyline_vertex = 0;

          pair_id_vertex = std::make_pair(
           id_next_polyline, n_next_polyline_vertex);
          vertices_connections[pair_id_vertex] = edge_segment;    

         }

        // Since it is not an initial or final vertex neither is
        // connected there is no need to look for the corresponding
        // 'edge_segment'
        triangulate_io.segmentlist[seg_lst_idx]=edge_segment;
        triangulate_io.segmentlist[seg_lst_idx+1]=edge_segment;

        // Increment the edge segment counter since we have used it
        edge_segment++;

        // Increment the segment list index
        seg_lst_idx+=2;

        // Store the marker list of the segments
        triangulate_io.segmentmarkerlist[count_tri_seg/2]=idpolyline+1;

        // Increment counter
        count_tri_seg+=2;

        // Increase the total number of segments
        n_boundglobalseg++;

       }

      // *********************************************************************
      // Final vertex on last polyline
      // *********************************************************************
      // Last vertex (Is it connected?)
      if (count_seg == n_polyline - 1 &&
        polyline_pt->is_final_vertex_connected())
       {
        // If it is connected then we already know its coordinates,
        // therefore we do not need to add them to the list but we
        // need to search for the corresponding edge_segment that
        // stores them

        // For this purpose we need the boudary_id of the boundary
        // to which the current boundary is connected and the
        // associated vertex number
        unsigned connected_id_polyline =
          polyline_pt->final_vertex_connected_bnd_id();

        unsigned connected_n_vertex =
          polyline_pt->final_vertex_connected_n_vertex();

        // Create the pair
        pair_id_vertex =
          std::make_pair(connected_id_polyline, connected_n_vertex);

        // Here is where we "recover" the 'edge_segment'
        unsigned recovered_edge_segment =
          vertices_connections[pair_id_vertex];

        // Compute the total number of vertices
        unsigned n_vertices = polyline_pt->nvertex();

        // Adding the information to the vertices_connections map.
        // The vertex "n_vertices-1" (last one in this special case) of
        // "boudary_id" = idpolyline is connected to the
        // "recovered_edge_segment"
        pair_id_vertex = std::make_pair(idpolyline, n_vertices-1);
        vertices_connections[pair_id_vertex] = recovered_edge_segment;

        // Here is where we use the corresponding "recovered_edge_segment"
        // number and create the connection in the "triangulate_io"
        // structure
        triangulate_io.segmentlist[seg_lst_idx]=recovered_edge_segment;

        // Increment the segment list index
        seg_lst_idx++;

	// Store the marker list of the segments
	triangulate_io.segmentmarkerlist[count_tri_seg/2]=idpolyline+1;

	// Increment counter
	count_tri_seg+=2;

	// Increase the total number of segments
	n_boundglobalseg++;

       }
      else if (count_seg == n_polyline - 1)
       {
        // Compute the total number of vertices
        unsigned n_vertices = polyline_pt->nvertex();

        // If it is not connected then we do not know its coordinates,
        // therefore we need to add them to the list (the last vertex
        // is "n_vertices-1"
        triangulate_io.pointlist[count_tri]=
          polyline_pt->vertex_coordinate(n_vertices-1)[0];
        triangulate_io.pointlist[count_tri+1]=
          polyline_pt->vertex_coordinate(n_vertices-1)[1];

        // Increment counter
        count_tri+=2;

        // Adding the information to the vertices_connections map.
        // The vertex "n_vertices-1" of the "boudary_id" = idpolyline is
        // connected to the "edge_segment"
        pair_id_vertex = std::make_pair(idpolyline, n_vertices-1);
        vertices_connections[pair_id_vertex] = edge_segment;

        // It is not connected to the edges so there is no need
        // to look for the corresponding edge_segment
        triangulate_io.segmentlist[seg_lst_idx] = edge_segment;

        // Increment the edge segment counter since we have used it
        edge_segment++;

        // Increment the segment list index
        seg_lst_idx++;

	// Store the marker list of the segments
	triangulate_io.segmentmarkerlist[count_tri_seg/2]=idpolyline+1;

	// Increment counter
	count_tri_seg+=2;

	// Increase the total number of segments
	n_boundglobalseg++;

       }

     }

   }

  // ***************************************************************
  // Check if the number of segments and previous number of segments
  // are the same
  // ***************************************************************
  if(n_boundglobalseg!=n_global_segments)
   {
    std::ostringstream error_stream;
    error_stream << "Error building TriangulateIO object.\n"
      << "There is a difference between the number of\n"
      << "global segments and the ones created in the\n"
      << "TriangulateIO object.\n"
      << "The numbers are: " << n_global_segments << " for "
      << "the total number of segments in PolyObjects and\n"
      << n_boundglobalseg << " for the total counted in the "
      << "TriangulateIO object.\n"
      << "Please, check TriangleMeshPolygon, TriangleMeshPolygon\n"
      << "list (holes) and TriangleMeshPolyLine (open curves)"
      <<std::endl;
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
   }

  // ****************************************************************
  // Store the region values in the TriangulateIO object (3th stage)
  // ****************************************************************

  // Number of extra regions
  unsigned n_regions = regions_coordinates.size();

  // Check for any defined region
  if(n_regions > 0)
   {
    triangulate_io.numberofregions = n_regions;
    triangulate_io.regionlist = (double*)
         malloc(triangulate_io.numberofregions * 4 * sizeof(double));

    std::map<unsigned, Vector<double> >::iterator it_regions;

    //Loop over the regions map
    unsigned p = 1;
    for(it_regions = regions_coordinates.begin(); 
        it_regions != regions_coordinates.end(); 
        it_regions++)
     {
      unsigned region_id = (*it_regions).first;
      triangulate_io.regionlist[4*p-4] =
        ((*it_regions).second)[0];
      triangulate_io.regionlist[4*p-3] =
        ((*it_regions).second)[1];
      triangulate_io.regionlist[4*p-2] =
        static_cast<double>(region_id);
      triangulate_io.regionlist[4*p-1] = regions_areas[region_id];
      p++;
     }

   }

  // ******************************************************************
  // Store the hole coordinates in the TriangulateIO object (4th stage)
  // ******************************************************************

  unsigned n_extra_holes = extra_holes_coordinates.size();
  unsigned n_real_holes = 0;
  Vector<unsigned> index_holes;

  // Count how many closed boundaries should be treated as holes
  for(unsigned c_hole=0;c_hole<n_holes;c_hole++)
   {
    if (!internal_polygon_pt[c_hole]->internal_point().empty())
     {
      n_real_holes++;
      index_holes.push_back(c_hole);
     }
   }

  // Storing the hole centre coordinates
  triangulate_io.numberofholes = n_extra_holes + n_real_holes;
  triangulate_io.holelist =
    (double*) malloc(triangulate_io.numberofholes * 2 * sizeof(double));

  unsigned count_hole = 0;
  for(;count_hole<n_real_holes*2;count_hole+=2)
   {
    unsigned index_hole = index_holes[count_hole/2];
    triangulate_io.holelist[count_hole]=
      internal_polygon_pt[index_hole]->
      internal_point()[0];

    triangulate_io.holelist[count_hole+1]=
      internal_polygon_pt[index_hole]->
      internal_point()[1];
   }

  unsigned c_extra_hole = 0;
  for(;count_hole<2*(n_extra_holes + n_real_holes);
    count_hole+=2,c_extra_hole++)
   {
    triangulate_io.holelist[count_hole] =
      extra_holes_coordinates[c_extra_hole][0];
    triangulate_io.holelist[count_hole+1] =
      extra_holes_coordinates[c_extra_hole][1];
   }

 }

 //========================================================================
 /// Create TriangulateIO object via the .poly file
 //========================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::build_triangulateio(
   const std::string& poly_file_name,
   TriangulateIO& triangulate_io)
  {

   // Process poly file
   // -----------------
   std::ifstream poly_file(poly_file_name.c_str(),std::ios_base::in);
   if(!poly_file)
    {
     throw OomphLibError("Error opening .poly file\n",
                         OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
    }

   // Initialize triangulateio structure
   TriangleHelper::initialise_triangulateio(triangulate_io);

   // Ignore the first line with structure description
   poly_file.ignore(80,'\n');

   // Read and store number of nodes
   unsigned invertices;
   poly_file>>invertices;
   triangulate_io.numberofpoints=invertices;

   // Initialisation of the point list
   triangulate_io.pointlist =
   (double *) malloc(triangulate_io.numberofpoints * 2 * sizeof(double));

   // Read and store spatial dimension of nodes
   unsigned mesh_dim;
   poly_file>>mesh_dim;

   if(mesh_dim == 0)
    {
     mesh_dim=2;
    }

#ifdef PARANOID
   if(mesh_dim!=2)
    {
     throw OomphLibError("The dimension must be 2\n",
                         OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // Read and check the flag for attributes
   unsigned nextras;
   poly_file>> nextras;

   triangulate_io.numberofpointattributes = 0;
   triangulate_io.pointattributelist = (double *) NULL;

   // Read and check the flag for boundary markers
   unsigned nodemarkers;
   poly_file>>nodemarkers;
   triangulate_io.pointmarkerlist = (int *) NULL;

#ifdef PARANOID
   // Reading the .poly with the oomph.lib we need
   // to set the point attribute and markers to 0
   if(nextras!=0 || nodemarkers!=0)
    {
     oomph_info << "===================================================="
     << std::endl<<std::endl;
     oomph_info <<"Reading the .poly file via oomph_lib \n"
     <<"point's attribute and point's markers \n"
     <<"are automatically set to 0"<<std::endl;
     oomph_info << "===================================================="
     <<std::endl;
    }
#endif

   // Dummy for node number (and attribute or markers if included)
   unsigned dummy_value;
   unsigned count_point=0;
   std::string test_string;

   // Skip line with commentary
   getline(poly_file,test_string,'#');
   poly_file.ignore(80,'\n');

   // Read and store all the nodes coordinates
   // (hole's vertices as well)
   for(unsigned count=0;count<invertices;count++)
    {
     poly_file>>dummy_value;
     poly_file>>triangulate_io.pointlist[count_point];
     poly_file>>triangulate_io.pointlist[count_point+1];
     if(nextras!=0 || nodemarkers!=0)
      {
       for(unsigned j=0;j<nextras;j++)
        {
         poly_file>>dummy_value;
        }
      }
     else if(nextras!=0 && nodemarkers!=0)
      {
       for(unsigned j=0;j<nextras;j++)
        {
         poly_file>>dummy_value;
         poly_file>>dummy_value;
        }
      }
     // Read the next line
     poly_file.ignore(80,'\n');

     // Skip line with commentary for internal box whether found
     if(poly_file.get() == '#')
      {
       poly_file.ignore(80,'\n');
      }
     // If read the char should be put back in the string

     else
      {
       poly_file.unget();
      }
     count_point+=2;
    }

   // The line with the segment's commentary has been skipped
   // by the command of the last loop

   // Read and store the number of segments
   unsigned dummy_seg;
   unsigned inelements;
   poly_file>>inelements;

   unsigned segment_markers;
   poly_file>>segment_markers;

   // Marker list should be provided by the user to assign
   // each segment to a boundary
#ifdef PARANOID
   if(segment_markers!=1)
    {

     std::ostringstream error_stream;
     error_stream
     <<"The segment marker should be provided \n"
     <<"In order to assign each segment to a boundary \n "<< std::endl;

     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
    }
#endif

   triangulate_io.numberofsegments = inelements;
   triangulate_io.segmentlist =
   (int *) malloc(triangulate_io.numberofsegments * 2 * sizeof(int));
   triangulate_io.segmentmarkerlist =
   (int *) malloc(triangulate_io.numberofsegments * sizeof(int));

   // Read all the segments edges and markers
   for(unsigned i=0;i<2*inelements;i+=2)
    {
     poly_file>>dummy_seg;
     poly_file>>triangulate_io.segmentlist[i];
     poly_file>>triangulate_io.segmentlist[i+1];
     if(segment_markers!=0)
      {
       poly_file>>triangulate_io.segmentmarkerlist[i/2];
      }

     //Skip line with commentary
     poly_file.ignore(80,'\n');
    }

   // Read and store the number of holes if given
   // Skip line with commentary
   if(getline(poly_file,test_string,'#'))
    {
     poly_file.ignore(80,'\n');

     unsigned dummy_hole;
     unsigned nhole;
     poly_file>>nhole;

     triangulate_io.numberofholes = nhole;
     triangulate_io.holelist =
     (double *) malloc(triangulate_io.numberofholes * 2 * sizeof(double));

     // Loop over the holes to get centre coords and store value onto the
     // TriangulateIO object
     for(unsigned i=0;i<2*nhole;i+=2)
      {
       poly_file>>dummy_hole;
       poly_file>>triangulate_io.holelist[i];
       poly_file>>triangulate_io.holelist[i+1];
      }
    }

   // Read and store the number of regions if given
   // Skip line with commentary
   if(getline(poly_file,test_string,'#'))
    {
     poly_file.ignore(80,'\n');

     unsigned dummy_region;
     unsigned nregion;
     poly_file>>nregion;
     std::cerr << "Regions: "<< nregion << std::endl;
     getchar();

     triangulate_io.numberofregions = nregion;
     triangulate_io.regionlist =
     (double *) malloc(triangulate_io.numberofregions * 4 * sizeof(double));

     // Loop over the regions to get coords and store value onto the
     // TriangulateIO object
     for(unsigned i=0;i<nregion;i++)
      {
       poly_file>>dummy_region;
       poly_file>>triangulate_io.regionlist[4*i];
       poly_file>>triangulate_io.regionlist[4*i+1];
       poly_file>>triangulate_io.regionlist[4*i+2];
       triangulate_io.regionlist[4*i+3] = 0.0;
      }
    }

  }

#endif

 //======================================================================
 /// Move the nodes on boundaries with associated Geometric Objects so
 /// that the exactly coincide with the geometric object. This requires
 /// that the boundary coordinates are set up consistently
 //======================================================================
 template<class ELEMENT>
  void
  TriangleMesh<ELEMENT>::snap_nodes_onto_geometric_objects()
  {
   //Loop over all boundaries
   const unsigned n_bound = this->nboundary();
   for (unsigned b = 0; b < n_bound; b++)
    {
     //Find the geometric object
     GeomObject* const geom_object_pt = this->boundary_geom_object_pt(b);

     //If there is one
     if (geom_object_pt != 0)
      {
       Vector<double> b_coord(1);
       Vector<double> new_x(2);
       const unsigned n_boundary_node = this->nboundary_node(b);
       for (unsigned n = 0; n < n_boundary_node; ++n)
        {
         //Get the boundary node and coordinates
         Node* const nod_pt = this->boundary_node_pt(b, n);
         nod_pt->get_coordinates_on_boundary(b, b_coord);

         //Get the position and time history according to the underlying
         //geometric object, assuming that it has the same timestepper
         //as the nodes....
         unsigned n_tvalues = 1 + nod_pt->position_time_stepper_pt()
          ->nprev_values();
         for (unsigned t = 0; t < n_tvalues; ++t)
          {
           //Get the position according to the underlying geometric object
           geom_object_pt->position(t, b_coord, new_x);

           //Move the node
           for (unsigned i = 0; i < 2; i++)
            {
             nod_pt->x(t, i) = new_x[i];
            }
          }
        }
      }
    }
  }

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

#ifdef OOMPH_HAS_TRIANGLE_LIB

//========================================================================
/// Build a new TriangulateIO object based on target areas specified
//========================================================================
template <class ELEMENT>
void RefineableTriangleMesh<ELEMENT>::refine_triangulateio(
  TriangulateIO& triangulate_io,
  const Vector<double>& target_area,
  struct TriangulateIO& triangle_refine)
 {

  //  Initialize
  TriangleHelper::initialise_triangulateio(triangle_refine);

  // Store the global number of vertices and segments
  // in the list
  unsigned n_points = triangulate_io.numberofpoints;
  triangle_refine.numberofpoints=n_points;

  unsigned n_segments=triangulate_io.numberofsegments;
  triangle_refine.numberofsegments=n_segments;

  // Initialization of the TriangulateIO objects to store the values
  triangle_refine.pointlist =
  (double *) malloc(triangulate_io.numberofpoints * 2 * sizeof(double));
  triangle_refine.pointmarkerlist =
  (int *) malloc(triangulate_io.numberofpoints * sizeof(int));
  triangle_refine.segmentlist =
  (int *) malloc(triangulate_io.numberofsegments * 2 * sizeof(int));
  triangle_refine.segmentmarkerlist =
  (int *) malloc(triangulate_io.numberofsegments * sizeof(int));

  // Storing the point's coordinates in the list
  // and in two vectors with x and y coordinates
  Vector<double> x_coord (n_points);
  Vector<double> y_coord (n_points);

  for(unsigned count_point=0;count_point<n_points*2;count_point++)
   {
    triangle_refine.pointlist[count_point]=
    triangulate_io.pointlist[count_point];

    // Even vaules represent the x coordinate
    // Odd values represent the y coordinate
    if (count_point%2==0)
     {
      x_coord[count_point/2] = triangulate_io.pointlist[count_point];
     }
    else
     {
      y_coord[(count_point-1)/2] = triangulate_io.pointlist[count_point];
     }
   }

  // Store the point's markers in the list
  for(unsigned count_marker=0;count_marker<n_points;count_marker++)
   {
    triangle_refine.pointmarkerlist[count_marker]=
    triangulate_io.pointmarkerlist[count_marker];
   }

  // Storing the segment's edges in the list
  for(unsigned count_seg=0;count_seg<n_segments*2;count_seg++)
   {
    triangle_refine.segmentlist[count_seg]=
    triangulate_io.segmentlist[count_seg];
   }

  // Store the segment's markers in the list
  for(unsigned count_markers=0;count_markers<n_segments;count_markers++)
   {
    triangle_refine.segmentmarkerlist[count_markers]=
    triangulate_io.segmentmarkerlist[count_markers];
   }

  // Store the hole's center coordinates
  unsigned n_holes = triangulate_io.numberofholes;
  triangle_refine.numberofholes = n_holes;

  triangle_refine.holelist =
  (double*) malloc(triangulate_io.numberofholes * 2 * sizeof(double));

  // Loop over the holes to get centre coords
  for(unsigned count_hole=0;count_hole<n_holes*2;count_hole++)
   {
    triangle_refine.holelist[count_hole] = triangulate_io.holelist[count_hole];
   }

  // Store the triangles values
  unsigned n_triangles = triangulate_io.numberoftriangles;
  triangle_refine.numberoftriangles = n_triangles;

#ifdef PARANOID
  if (n_triangles!=target_area.size())
   {
    std::stringstream err;
    err << "Number of triangles in triangulate_io="
    << n_triangles << " doesn't match\n"
    << "size of target area vector ("
    << target_area.size() << ")\n";
    throw OomphLibError(
      err.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
   }
#endif

  unsigned n_corners = triangulate_io.numberofcorners;
  triangle_refine.numberofcorners = n_corners;

  triangle_refine.trianglelist =
  (int *) malloc(triangulate_io.numberoftriangles * 3 * sizeof(int));

  // Store the triangle's corners in the list and get element sizes
  for(unsigned count_tri=0;count_tri<n_triangles*3;count_tri++)
   {
    triangle_refine.trianglelist[count_tri]=
    triangulate_io.trianglelist[count_tri];
   }

  // Store the triangle's area in the list
  triangle_refine.trianglearealist =
  (double *) malloc(triangulate_io.numberoftriangles * sizeof(double));
  for(unsigned count_area=0;count_area<n_triangles;count_area++)
   {
    triangle_refine.trianglearealist[count_area]=target_area[count_area];
   }

  // Store the triangles attributes in the list
  triangle_refine.numberoftriangleattributes =
  triangulate_io.numberoftriangleattributes;

  triangle_refine.triangleattributelist =
  (double *) malloc(
    triangulate_io.numberoftriangles *
    triangulate_io.numberoftriangleattributes * sizeof(double));
  for(unsigned count_attribute=0;
    count_attribute<(n_triangles*triangulate_io.numberoftriangleattributes);
    count_attribute++)
   {
    triangle_refine.triangleattributelist[count_attribute] =
    triangulate_io.triangleattributelist[count_attribute];
   }

 }

#endif

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#ifdef OOMPH_HAS_TRIANGLE_LIB

//======================================================================
/// Adapt problem based on specified elemental error estimates
//======================================================================
 template <class ELEMENT>
 void RefineableTriangleMesh<ELEMENT>::adapt(
  const Vector<double>& elem_error)
 {

  double t_start_overall=TimingHelpers::timer();

  // Get refinement targets
  Vector<double> target_area(elem_error.size());
  double min_angle=compute_area_target(elem_error,
                                       target_area);


  // Post-process to allow only quantised target areas
  // in an attempt to more closely mimick the structured
  // case and limit the diffusion of small elements.
  bool quantised_areas=true; 
  if (quantised_areas)
   {
    unsigned n=target_area.size();
    double total_area=0;
    for (unsigned e=0;e<n;e++)
     {
      total_area+=this->finite_element_pt(e)->size();
     }
    for (unsigned e=0;e<n;e++)
     {
      unsigned level=unsigned(ceil(log(target_area[e]/total_area)/log(1.0/3.0)))-1;
      double new_target_area=total_area*pow(1.0/3.0,int(level));
      target_area[e]=new_target_area;
     }
   }

  // std::ofstream hierher;
  // hierher.open((Global_string_for_annotation:: String[0]+"overall_target_areas"+
  //               StringConversion::to_string(Global_unsigned::Number)+".dat").c_str());
  
  // Get maximum target area
  unsigned n=target_area.size();
  double max_area=0.0;
  double min_area=DBL_MAX;
  for (unsigned e=0;e<n;e++)
   {
    if (target_area[e]>max_area) max_area=target_area[e];
    if (target_area[e]<min_area) min_area=target_area[e];

    // hierher << (finite_element_pt(e)->node_pt(0)->x(0)+
    //            finite_element_pt(e)->node_pt(1)->x(0)+
    //           finite_element_pt(e)->node_pt(2)->x(0))/3.0 << " "
    //        << (finite_element_pt(e)->node_pt(0)->x(1)+
    //            finite_element_pt(e)->node_pt(1)->x(1)+
    //            finite_element_pt(e)->node_pt(2)->x(1))/3.0 << " "
    //        << target_area[e] << " " 
    //        << finite_element_pt(e)->size() << " "
    //        << elem_error[e]  << " " << std::endl;
   }

  //hierher.close();


  oomph_info << "Maximum target area: " << max_area << std::endl;
  oomph_info << "Minimum target area: " << min_area << std::endl;
  oomph_info << "Number of elements to be refined "
             << this->Nrefined << std::endl;
  oomph_info << "Number of elements to be unrefined "
             << this->Nunrefined << std::endl;
  oomph_info << "Min angle "<< min_angle << std::endl;

  double orig_max_area, orig_min_area;
  this->max_and_min_element_size(orig_max_area, orig_min_area);
  oomph_info << "Max/min element size in original mesh: "
             << orig_max_area << " "
             << orig_min_area << std::endl;

  // Check if boundaries need to be updated (regardless of
  // requirements of bulk error estimator) but don't do anything!
  bool check_only=true;
  bool outer_boundary_update_necessary= false;
  bool inner_boundary_update_necessary= false;
  bool inner_open_boundary_update_necessary=false;

  outer_boundary_update_necessary=
   this->update_polygon_using_face_mesh(this->Outer_boundary_pt,check_only);

  // Do not waste time if we already know that it is necessary an update
  // on the boundary representation
  if (!outer_boundary_update_necessary)
   {
    // Check if we need to generate a new 1D mesh representation of
    // the inner hole boundaries
    unsigned nhole=this->Internal_polygon_pt.size();
    Vector<Vector<double> > internal_point_coord(nhole);
    inner_boundary_update_necessary=
     this->surface_remesh_for_inner_hole_boundaries(internal_point_coord,
                                                    check_only);

    // If there was not necessary a change even on the internal closed
    // curve then finally check for the open curves as well
    if (!inner_boundary_update_necessary)
     {
      unsigned n_open_polyline = this->Internal_open_curve_pt.size();
      for (unsigned i = 0; i < n_open_polyline; i++)
       {
        inner_open_boundary_update_necessary=
         this->update_open_curve_using_face_mesh(
          this->Internal_open_curve_pt[i], check_only);

        // If at least one needs modification then break the for loop
        if (inner_open_boundary_update_necessary) break;
       }
     }
   }

  // Should we bother to adapt?
  if ( (Nrefined > 0) || (Nunrefined > max_keep_unrefined()) ||
       (min_angle < min_permitted_angle()) || (outer_boundary_update_necessary)
       || (inner_boundary_update_necessary)
       || (inner_open_boundary_update_necessary) )
   {

    if (! ( (Nrefined > 0) || (Nunrefined > max_keep_unrefined()) ) )
     {

      if ( (outer_boundary_update_necessary)
           || (inner_boundary_update_necessary)
           || (inner_open_boundary_update_necessary) )
       {
        oomph_info
         << "Mesh regeneration triggered by inaccurate interface/surface\n"
         << "representation; setting Nrefined to number of elements.\n"
         << "outer_boundary_update_necessary     : " 
         << outer_boundary_update_necessary << "\n"
         << "inner_boundary_update_necessary     : "
         << inner_boundary_update_necessary << "\n"
         << "inner_open_boundary_update_necessary: "
         << inner_open_boundary_update_necessary << "\n"; 
        Nrefined=nelement();
       }
      else
       {
        oomph_info
         << "Mesh regeneration triggered by min angle criterion;\n"
         << "setting Nrefined to number of elements.\n";
        Nrefined=nelement();
       }
     }

    //Generate a new 1D mesh representation of the inner hole boundaries
    unsigned nhole=this->Internal_polygon_pt.size();
    Vector<Vector<double> > internal_point_coord(nhole);
    this->surface_remesh_for_inner_hole_boundaries(internal_point_coord);

    //Update the representation of the outer boundary
    this->update_polygon_using_face_mesh(this->Outer_boundary_pt);

    // After updating outer and internal closed boundaries it is also
    // necessary to update internal boundaries.
    unsigned n_open_polyline = this->Internal_open_curve_pt.size();
    for (unsigned i = 0; i < n_open_polyline; i++)
     {
      this->update_open_curve_using_face_mesh(
       this->Internal_open_curve_pt[i]);
     }

    //If there is not a geometric object associated with the boundary
    //then reset the boundary coordinates so that the lengths are consistent
    //in the new mesh and the old mesh.
    const unsigned n_boundary = this->nboundary();
    for(unsigned b=0;b<n_boundary;++b)
     {
      if(this->boundary_geom_object_pt(b)==0)
       {
        this->setup_boundary_coordinates(b);
       }
     }

    // Update the region information by setting the coordinates from the
    // centroid of the first element in each region (which should allow
    // automatic updates when the regions deform)
    {
     unsigned n_region = this->nregion();
     if(n_region > 1)
      {
       for(std::map<unsigned, Vector<double> >::iterator it =
            this->Regions_coordinates.begin();
           it!=this->Regions_coordinates.end(); ++it)
        {
         //Storage for the approximate centroid
         Vector<double> centroid(2,0.0);

         //Get the region id
         unsigned region_id = it->first;

         //Report information
         oomph_info << "Region " << region_id << ": "
                    << it->second[0] << " " << it->second[1] << " ";
         
         //Check that there is at least one element in the region
         unsigned n_region_element = this->nregion_element(region_id);
         if(n_region_element > 0)
          {
           //Cache pointer to the first element
           FiniteElement* const elem_pt = this->region_element_pt(region_id,0);

           //Loop over the corners of the triangle and average
           for(unsigned n=0;n<3;n++)
            {
             Node* const nod_pt = elem_pt->node_pt(n);
             for(unsigned i=0;i<2;i++) {centroid[i] += nod_pt->x(i);}
            }
           for(unsigned i=0;i<2;i++) {centroid[i] /= 3;} 
           //Now we have the centroid set it
           it->second = centroid;

           oomph_info << "   ,    " << 
            it->second[0] << " " << it->second[1] << std::endl;
          } //end of case when there is at least one element
        }
      }
    }
              

    // Are we dealing with a solid mesh?
    SolidMesh* solid_mesh_pt=dynamic_cast<SolidMesh*>(this);

    // Build temporary uniform background mesh
    //----------------------------------------
    // with area set by maximum required area
    //---------------------------------------
    RefineableTriangleMesh<ELEMENT>* tmp_new_mesh_pt=0;

    TriangleMeshClosedCurve* closed_curve_pt=this->Outer_boundary_pt;
    unsigned nh=this->Internal_polygon_pt.size();
    Vector<TriangleMeshClosedCurve*> hole_pt(nh);
    for (unsigned i=0;i<nh;i++)
     {
      hole_pt[i]=this->Internal_polygon_pt[i];
     }
    unsigned nib=this->Internal_open_curve_pt.size();
    Vector<TriangleMeshOpenCurve*> open_curves_pt(nib);
    for (unsigned i=0;i<nib;i++)
     {
      open_curves_pt[i]=this->Internal_open_curve_pt[i];
     }

    // *****************************************************************
    // Gather all the information and use the TriangleMeshParameters
    // object which help us on the manage of all TriangleMesh object's
    // information
    
    // Create the TriangleMeshParameters objects with the outer boundary
    // as the only one parameter
    TriangleMeshParameters triangle_mesh_parameters(closed_curve_pt);
    
    // Pass information about the holes
    triangle_mesh_parameters.internal_closed_curve_pt() = hole_pt;
    
    // Pass information about the internal open boundaries
    triangle_mesh_parameters.internal_open_curves_pt() = open_curves_pt;
    
    // Pass information about max allowed area for elements
    triangle_mesh_parameters.element_area() = max_area;
    
    // Pass information about the extra holes (not defined with closed
    // boundaries)
    triangle_mesh_parameters.extra_holes_coordinates() =
     this->Extra_holes_coordinates;
    
    //Pass information about regions
    triangle_mesh_parameters.regions_coordinates() =
     this->Regions_coordinates;
    
    // *****************************************************************

    if (solid_mesh_pt!=0)
     {
      tmp_new_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>
       (triangle_mesh_parameters, this->Time_stepper_pt);
     }
    else
     {
      tmp_new_mesh_pt=new RefineableTriangleMesh<ELEMENT>
       (triangle_mesh_parameters, this->Time_stepper_pt);
     }

    // Snap to curvilinear boundaries (some code duplication as this
    // is repeated below but helper function would take so many
    // arguments that it's nearly as messy...

    //Pass the boundary geometric objects to the new mesh
    tmp_new_mesh_pt->boundary_geom_object_pt() =
     this->boundary_geom_object_pt();

    //Reset the boundary coordinates if there is
    //a geometric object associated with the boundary
    tmp_new_mesh_pt->boundary_coordinate_limits() =
     this->boundary_coordinate_limits();
    for (unsigned b=0;b<n_boundary;b++)
     {
      if(tmp_new_mesh_pt->boundary_geom_object_pt(b)!=0)
       {
        tmp_new_mesh_pt->setup_boundary_coordinates(b);
       }
     }

    //Output the mesh before any snapping takes place
    //   tmp_new_mesh_pt->output("pre_mesh_nodes_snapped_0.dat");

    //Move the nodes on the new boundary onto the
    //old curvilinear boundary. If the boundary is straight this
    //will do precisely nothing but will be somewhat inefficient
    for(unsigned b=0;b<n_boundary;b++)
     {
      this->snap_nodes_onto_boundary(tmp_new_mesh_pt,b);
     }

    // Update mesh further?
    if (Mesh_update_fct_pt!=0)
     {
      Mesh_update_fct_pt(tmp_new_mesh_pt);
     }

    //If we have a continuation problem 
    //any problem in which the timestepper is a "generalisedtimestepper",
    //which will have been set by the problem, then ensure
    //all data in the new mesh has the appropriate timestepper
    /*if(dynamic_cast<GeneralisedTimeStepper*>(this->Time_stepper_pt))
      {
      tmp_new_mesh_pt->set_nodal_and_elemental_time_stepper(
      this->Time_stepper_pt);
      tmp_new_mesh_pt->set_mesh_level_time_stepper(this->Time_stepper_pt);
      }*/

    //Output the mesh after the snapping has taken place
    //tmp_new_mesh_pt->output("mesh_nodes_snapped_0.dat"); 
    //this->output("existing_mesh.dat"); 

    // Get the TriangulateIO object associated with that mesh
    TriangulateIO tmp_new_triangulateio=
     tmp_new_mesh_pt->triangulateio_representation();
    RefineableTriangleMesh<ELEMENT>* new_mesh_pt=0;

    // Map storing target areas for elements in temporary
    // TriangulateIO mesh
    std::map<GeneralisedElement*,double> target_area_map;


    
    // Adjust size of bins
    unsigned backup_bin_x=Multi_domain_functions::Nx_bin;
    unsigned backup_bin_y=Multi_domain_functions::Ny_bin;
    Multi_domain_functions::Nx_bin=Nbin_x_for_area_transfer;
    Multi_domain_functions::Ny_bin=Nbin_y_for_area_transfer;
    
    // Make a mesh as geom object representation of the temporary
    // mesh -- this also builds up the internal bin structure 
    // from which we'll recover the target areas
    double t0_geom_obj=TimingHelpers::timer();
    
    // If the mesh is a solid mesh then do the mapping based on the
    // Eulerian coordinates
    bool backup= MeshAsGeomObject::Use_eulerian_coordinates_during_setup;
    if (solid_mesh_pt!=0)
     {
      MeshAsGeomObject::Use_eulerian_coordinates_during_setup=true;
     }
    MeshAsGeomObject* mesh_geom_obj_pt =
     new MeshAsGeomObject(this);
    if (solid_mesh_pt!=0)
     {
      MeshAsGeomObject::Use_eulerian_coordinates_during_setup=backup;
     }

    // Reset
    Multi_domain_functions::Nx_bin=backup_bin_x;
    Multi_domain_functions::Ny_bin=backup_bin_y;
    
    oomph_info << "time for setup of mesh as geom obj: "
               << TimingHelpers::timer()-t0_geom_obj
               << std::endl;
    
    // Do some stats
    {
     unsigned max_entry=0;
     unsigned min_entry=UINT_MAX;
     unsigned tot_entry=0;
     unsigned nempty=0;
     Vector<Vector<std::pair<FiniteElement*,Vector<double> > > > 
      bin_content=mesh_geom_obj_pt->bin_content();
     mesh_geom_obj_pt->bin_content();
     unsigned nbin=bin_content.size();
     for (unsigned b=0;b<nbin;b++)
      {
       unsigned nentry=bin_content[b].size();
       if (nentry>max_entry) max_entry=nentry;
       if (nentry<min_entry) min_entry=nentry;
       if (nentry==0) nempty++;
       tot_entry+=nentry;
      }
     oomph_info 
      << "Before bin diffusion: nbin, nempty, min, max, av entries: "
      << nbin << " " 
      << nempty << " " 
      << min_entry << " " 
      << max_entry << " " 
      << double(tot_entry)/double(nbin) << " "
      << std::endl;
    }
    
    // Fill bin by diffusion
    double t0_bin_diff=TimingHelpers::timer();
    oomph_info << "going into diffusion bit...\n";
    mesh_geom_obj_pt->fill_bin_by_diffusion();
    oomph_info << "back from diffusion bit...\n";
    oomph_info << "time for bin diffusion: "
               << TimingHelpers::timer()-t0_bin_diff
               << std::endl;
    
    

    // Do some stats
    {
     unsigned max_entry=0;
     unsigned min_entry=UINT_MAX;
     unsigned tot_entry=0;
     unsigned nempty=0;
     Vector<Vector<std::pair<FiniteElement*,Vector<double> > > > 
      bin_content=mesh_geom_obj_pt->bin_content();
     mesh_geom_obj_pt->bin_content();
     unsigned nbin=bin_content.size();
     for (unsigned b=0;b<nbin;b++)
      {
       unsigned nentry=bin_content[b].size();
       if (nentry>max_entry) max_entry=nentry;
       if (nentry<min_entry) min_entry=nentry;
       if (nentry==0) nempty++;
       tot_entry+=nentry;
      }
     oomph_info 
      << "After bin diffusion: nbin, nempty, min, max, av entries: "
      << nbin << " " 
      << nempty << " " 
      << min_entry << " " 
      << max_entry << " " 
      << double(tot_entry)/double(nbin) << " "
      << std::endl;
    }
    

    // Set up a map from pointer to element to its number
    // in the mesh
    std::map<GeneralisedElement*,unsigned> element_number;
    unsigned nelem=this->nelement();
    for (unsigned e=0;e<nelem;e++)
     {
      element_number[this->element_pt(e)]=e;
     }
    

    // Now start iterating to refine mesh recursively
    //-----------------------------------------------
    bool done=false;
    unsigned iter=0;
    double t_iter=TimingHelpers::timer();
    while (!done)
     {

      // Accept by default but overwrite if things go wrong below
      done=true;
      double t_start=TimingHelpers::timer();




      double t0_loop_int_pts=TimingHelpers::timer();


      // Get ready for next assignment of target areas
      target_area_map.clear();
          
      // Loop over elements in new (tmp) mesh and visit all
      // its integration points. Check where it's located in the bin
      // structure of the current mesh and pass the target area
      // to the new element
      nelem=tmp_new_mesh_pt->nelement();

      for (unsigned e=0;e<nelem;e++)
       { // start loop el
        ELEMENT* el_pt=dynamic_cast<ELEMENT*>(tmp_new_mesh_pt->element_pt(e));
        unsigned nint=el_pt->integral_pt()->nweight();
        for (unsigned ipt=0;ipt<nint;ipt++)
         {
          // Get the coordinate of current point
          Vector<double> s(2);
          for(unsigned i=0;i<2;i++)
           {
            s[i] = el_pt->integral_pt()->knot(ipt,i);
           }

          Vector<double> x(2);
          el_pt->interpolated_x(s,x);

          // Find the bin that contains that point and its contents
          int bin_number=0;
          Vector<std::pair<FiniteElement*,Vector<double> > >
           sample_point_pairs;
          mesh_geom_obj_pt->get_bin(x,bin_number,sample_point_pairs);
              
          // Did we find it?
          if (bin_number<0)
           {
            // Not even within bin boundaries... odd
            std::stringstream error_message;
            error_message
             << "Very odd -- we're looking for a point[ "
             << x[0] << " " << x[1] << " ] that's not even \n"
             << "located within the bin boundaries.\n";
            throw OomphLibError(error_message.str(),
                                "RefineableTriangleMesh::adapt()",
                                OOMPH_EXCEPTION_LOCATION);
           }
          else
           {
            // Pass target area to new element
            unsigned n=sample_point_pairs.size();
            if (n>0)
             {
              for (unsigned ee=0;ee<n;ee++)
               {
                // Get ee-th element (in currrent mesh) in bin
                GeneralisedElement* current_el_pt=sample_point_pairs[ee].first;
                unsigned e_current=element_number[current_el_pt];                    

                // Go for smallest target area of any element in this bin
                // to force "one level" of refinement (the one-level-ness is
                // enforced below by limiting the actual reduction in area
                if (target_area_map[el_pt]!=0)
                 {
                  target_area_map[el_pt]=
                   std::min(target_area_map[el_pt], 
                            target_area[e_current]);
                 }
                else
                 {
                  target_area_map[el_pt]=target_area[e_current];
                 }
               }
             }
            else
             {
              std::stringstream error_message;
              error_message
               << "Point not found within bin structure\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
             }
           }
         }
       }


      oomph_info << "time for loop over int pts in new mesh: "
                 << TimingHelpers::timer()-t0_loop_int_pts  
                 << std::endl;


      // //hierher
      // {
      // hierher.open((Global_string_for_annotation:: String[0]+"binned_target_areas"+
      //               StringConversion::to_string(Global_unsigned::Number)+".dat").c_str());
  
      // Vector<Vector<std::pair<FiniteElement*,Vector<double> > > > bin_content=
      //  mesh_geom_obj_pt->bin_content();
      // unsigned nbin=bin_content.size();
      // for (unsigned b=0;b<nbin;b++)
      //  {
      //   unsigned nentry=bin_content[b].size();
      //   for (unsigned entry=0;entry<nentry;entry++)
      //    {
      //     FiniteElement* el_pt=bin_content[b][entry].first;
      //     GeneralisedElement* gen_el_pt=bin_content[b][entry].first;
      //     Vector<double> s=bin_content[b][entry].second;
      //     Vector<double> x(2);
      //     el_pt->interpolated_x(s,x);
      //     unsigned e_current=element_number[gen_el_pt];
      //     hierher << x[0] << " " << x[1] << " " 
      //             << target_area[e_current] << " " 
      //             << el_pt->size() << " " 
      //             << std::endl;
      //    }
      //  }
      // hierher.close();
      // }
                    

      oomph_info << "CPU for transfer of target areas (iter "
                 << iter << ") " << TimingHelpers::timer()-t_start
                 << std::endl;


      // // Output mesh
      // tmp_new_mesh_pt->output(("intermediate_mesh"+
      //                         StringConversion::to_string(iter)+".dat").c_str()); 
      
      // hierher.open((Global_string_for_annotation:: String[0]+"target_areas_intermediate_mesh_iter"+
      //               StringConversion::to_string(iter)+"_"+
      //               StringConversion::to_string(Global_unsigned::Number)+".dat").c_str());
     
  
      // Now copy into target area for temporary mesh but limit to
      // the equivalent of one sub-division per iteration
      unsigned nel_new=tmp_new_mesh_pt->nelement();
      Vector<double> new_target_area(nel_new);
      for (unsigned e=0;e<nel_new;e++)
       {
        // No target area found for this element -- keep its size
        // by setting target area to -1 for triangle
        double new_area=target_area_map[tmp_new_mesh_pt->element_pt(e)];
        if (new_area<=0.0)
         {
          std::ostringstream error_stream;
          error_stream << "This shouldn't happen! Element whose centroid is at"
                       <<  (tmp_new_mesh_pt->finite_element_pt(e)->node_pt(0)->x(0)+
                            tmp_new_mesh_pt->finite_element_pt(e)->node_pt(1)->x(0)+
                            tmp_new_mesh_pt->finite_element_pt(e)->node_pt(2)->x(0))/3.0 << " "
                       << (tmp_new_mesh_pt->finite_element_pt(e)->node_pt(0)->x(1)+
                           tmp_new_mesh_pt->finite_element_pt(e)->node_pt(1)->x(1)+
                           tmp_new_mesh_pt->finite_element_pt(e)->node_pt(2)->x(1))/3.0 << " "
                       << " has no target area assigned\n"; 
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
         }
        else
         {
          // Limit target area to the equivalent of uniform
          // refinement during this stage of the iteration
          new_target_area[e]=new_area;
          if (new_target_area[e]<
              tmp_new_mesh_pt->finite_element_pt(e)->size()/3.0)
           {
            new_target_area[e]=
             tmp_new_mesh_pt->finite_element_pt(e)->size()/3.0;

            // We'll need to give it another go later
            done=false;
           }
         }               
       
      
        // hierher 
        // << (tmp_new_mesh_pt->finite_element_pt(e)->node_pt(0)->x(0)+
        //     tmp_new_mesh_pt->finite_element_pt(e)->node_pt(1)->x(0)+
        //     tmp_new_mesh_pt->finite_element_pt(e)->node_pt(2)->x(0))/3.0 << " "
        // << (tmp_new_mesh_pt->finite_element_pt(e)->node_pt(0)->x(1)+
        //     tmp_new_mesh_pt->finite_element_pt(e)->node_pt(1)->x(1)+
        //     tmp_new_mesh_pt->finite_element_pt(e)->node_pt(2)->x(1))/3.0 << " "
        // << new_target_area[e] << " " 
        // << tmp_new_mesh_pt->finite_element_pt(e)->size()<< std::endl;
        
       }	
      
      if (done) 
       {
        oomph_info << "All area adjustments accomodated by max. permitted area reduction \n";
       }
      else
       {
        oomph_info << "NOT all area adjustments accomodated by max. permitted area reduction \n";
       }
      
      //hierher.close();
      //pause("doced binned_target_areas.dat and interemdiate mesh targets");

      // Now create the new mesh from TriangulateIO structure
      //-----------------------------------------------------
      // associated with uniform background mesh and the
      //------------------------------------------------
      // associated target element sizes.
      //---------------------------------

      // Solid mesh?
      if (solid_mesh_pt!=0)
       {
        new_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>
         (new_target_area,
          tmp_new_triangulateio,
          this->Time_stepper_pt,
          this->Use_attributes);
       }
      // No solid mesh
      else
       {
        new_mesh_pt=new RefineableTriangleMesh<ELEMENT>
         (new_target_area,
          tmp_new_triangulateio,
          this->Time_stepper_pt,
          this->Use_attributes);
       }

      // Snap to curvilinear boundaries (some code duplication as this
      // is repeated below but helper function would take so many
      // arguments that it's nearly as messy...

      //Pass the boundary  geometric objects to the new mesh 
      new_mesh_pt->boundary_geom_object_pt() =
       this->boundary_geom_object_pt();

      // Reset the boundary coordinates if there is
      // a geometric object associated with the boundary
      new_mesh_pt->boundary_coordinate_limits() =
       this->boundary_coordinate_limits();
      for (unsigned b=0;b<n_boundary;b++)
       {
        if(new_mesh_pt->boundary_geom_object_pt(b)!=0)
         {
          new_mesh_pt->setup_boundary_coordinates(b);
         }
       }

      //Output the mesh before any snapping takes place
      //   new_mesh_pt->output("pre_mesh_nodes_snapped_1.dat"); 

      //Move the nodes on the new boundary onto the 
      //old curvilinear boundary. If the boundary is straight this 
      //will do precisely nothing but will be somewhat inefficient
      for(unsigned b=0;b<n_boundary;b++)
       {
        this->snap_nodes_onto_boundary(new_mesh_pt,b);
       }

      // Update mesh further?
      if (Mesh_update_fct_pt!=0)
       {
        Mesh_update_fct_pt(new_mesh_pt);
       }

      //If we have a continuation problem 
      //any problem in which the timestepper is a "generalisedtimestepper",
      //which will have been set by the problem, then ensure
      //all data in the new mesh has the appropriate timestepper
      if(dynamic_cast<GeneralisedTimeStepper*>(this->Time_stepper_pt))
       {
        new_mesh_pt->set_nodal_and_elemental_time_stepper(
         this->Time_stepper_pt,false);
        new_mesh_pt->set_mesh_level_time_stepper(this->Time_stepper_pt,false);
       }

      //Output the mesh after the snapping has taken place
      //   new_mesh_pt->output("mesh_nodes_snapped_1.dat"); 

      // Not done: get ready for another iteration
      iter++;
      delete tmp_new_mesh_pt;
      if (!done)
       {
        tmp_new_mesh_pt=new_mesh_pt;
        tmp_new_triangulateio=new_mesh_pt->triangulateio_representation();
       }

     } // end of iteration

    //Delete the temporary geometric object representation of the current mesh
    delete mesh_geom_obj_pt;
    
    oomph_info << "CPU for iterative generation of new mesh " 
               << TimingHelpers::timer()-t_iter
               << std::endl;

    double t_proj=TimingHelpers::timer();
    oomph_info << "about to do projection\n";

    // Project current solution onto new mesh
    //---------------------------------------
    ProjectionProblem<ELEMENT>* project_problem_pt=
     new ProjectionProblem<ELEMENT>;
    project_problem_pt->mesh_pt()=new_mesh_pt;
    //project_problem_pt->disable_suppress_output_during_projection();
    project_problem_pt->project(this);

    oomph_info << "CPU for projection of solution onto new mesh " 
               << TimingHelpers::timer()-t_proj
               << std::endl;

    double t_rest=TimingHelpers::timer();

    //Flush the old mesh 
    unsigned nnod=nnode();
    for(unsigned j=nnod;j>0;j--)
     {
      delete Node_pt[j-1];
      Node_pt[j-1] = 0;
     }
    unsigned nel=nelement();
    for(unsigned e=nel;e>0;e--)
     {
      delete Element_pt[e-1];
      Element_pt[e-1] = 0;
     }

    // Now copy back to current mesh
    //------------------------------
    nnod=new_mesh_pt->nnode();
    Node_pt.resize(nnod);
    nel=new_mesh_pt->nelement();
    Element_pt.resize(nel);
    for(unsigned j=0;j<nnod;j++)
     {
      Node_pt[j] = new_mesh_pt->node_pt(j);
     }
    for(unsigned e=0;e<nel;e++)
     {
      Element_pt[e] = new_mesh_pt->element_pt(e);
     }

    //Copy the boundary schemes
    unsigned nbound=new_mesh_pt->nboundary();
    Boundary_element_pt.resize(nbound);
    Face_index_at_boundary.resize(nbound);
    Boundary_node_pt.resize(nbound);
    for (unsigned b=0;b<nbound;b++)
     {
      unsigned nel=new_mesh_pt->nboundary_element(b);
      Boundary_element_pt[b].resize(nel);
      Face_index_at_boundary[b].resize(nel);
      for (unsigned e=0;e<nel;e++)
       {
        Boundary_element_pt[b][e]=new_mesh_pt->boundary_element_pt(b,e);
        Face_index_at_boundary[b][e]=new_mesh_pt->face_index_at_boundary(b,e);
       }
      unsigned nnod=new_mesh_pt->nboundary_node(b);
      Boundary_node_pt[b].resize(nnod);
      for (unsigned j=0;j<nnod;j++)
       {
        Boundary_node_pt[b][j]=new_mesh_pt->boundary_node_pt(b,j);
       }
     }

    //Also copy over the new boundary and region information
    unsigned n_region = new_mesh_pt->nregion();
    //Only bother if we have regions
    if(n_region > 1)
     {
      //Deal with the region information first
      this->Region_attribute.resize(n_region);
      for(unsigned r=0;r<n_region;r++)
       {
        this->Region_attribute[r] = new_mesh_pt->region_attribute(r);
        // Get the region id
        unsigned r_id = static_cast<unsigned>(this->Region_attribute[r]);
        //Find the number of elements in the region
        unsigned n_region_element = new_mesh_pt->nregion_element(r_id);
        this->Region_element_pt[r_id].resize(n_region_element);
        for(unsigned e=0;e<n_region_element;e++)
         {
          this->Region_element_pt[r_id][e] = 
           new_mesh_pt->region_element_pt(r_id,e);
         }
       }

      //Now the boundary region information
      this->Boundary_region_element_pt.resize(nbound);
      this->Face_index_region_at_boundary.resize(nbound);

      //Now loop over the boundaries
      for(unsigned b=0;b<nbound;++b)
       {
        for (unsigned rr = 0 ; rr < n_region; rr++)
         {
          // The region id
          unsigned r = static_cast<unsigned>(this->Region_attribute[rr]);
          
          unsigned n_boundary_el_in_region =
           new_mesh_pt->nboundary_element_in_region(b,r);

          if(n_boundary_el_in_region > 0)
           {
            //Allocate storage in the map
            this->Boundary_region_element_pt[b][r].
             resize(n_boundary_el_in_region);
            this->Face_index_region_at_boundary[b][r].
             resize(n_boundary_el_in_region);

            //Copy over the information
            for(unsigned e=0;e<n_boundary_el_in_region;++e)
             {
              this->Boundary_region_element_pt[b][r][e]
               = new_mesh_pt->boundary_element_in_region_pt(b,r,e);
              this->Face_index_region_at_boundary[b][r][e]
               = new_mesh_pt->face_index_at_boundary_in_region(b,r,e);
             }
           }
         }
       } //End of loop over boundaries

     } //End of case when more than one region

    // Copy and create new associated boundary

    //Snap the newly created nodes onto any geometric objects
    this->snap_nodes_onto_geometric_objects();

    // Copy the IDs of the vertex nodes
    this->Oomph_vertex_nodes_id=new_mesh_pt->oomph_vertex_nodes_id();

    // Copy TriangulateIO representation
    TriangleHelper::clear_triangulateio(this->Triangulateio);
    bool quiet=true;
    this->Triangulateio=
     TriangleHelper::deep_copy_of_triangulateio_representation(
      new_mesh_pt->triangulateio_representation(),quiet);

    // Flush the mesh
    new_mesh_pt->flush_element_and_node_storage();

    // Delete the mesh and the problem
    delete new_mesh_pt;
    delete project_problem_pt;

    double max_area;
    double min_area;
    this->max_and_min_element_size(max_area, min_area);
    oomph_info << "Max/min element size in adapted mesh: "
               << max_area << " "
               << min_area << std::endl;

    oomph_info << "CPU for final bits... " 
               << TimingHelpers::timer()-t_rest
               << std::endl;


   }
  else
   {
    oomph_info << "Not enough benefit in adaptation.\n";
    Nrefined=0;
    Nunrefined=0;
   }

  oomph_info <<"CPU for adaptation: "
             << TimingHelpers::timer()-t_start_overall << std::endl;

 }

//=========================================================================
/// \short Helper function that updates the input polygon's PSLG
/// by using the end-points of elements from FaceMesh(es) that are
/// constructed for the boundaries associated with the segments of the
/// polygon. Optional boolean is used to run it as test only (if
/// true is specified as input) in which case polygon isn't actually
/// modified. Returned boolean indicates if polygon was (or would have
/// been -- if called with check_only=false) changed.
//=========================================================================
template<class ELEMENT>
bool RefineableTriangleMesh<ELEMENT>::
update_polygon_using_face_mesh(TriangleMeshPolygon* polygon_pt,
  const bool& check_only)
  {
 // Boolean that indicates whether an actual update of the polygon
 // was performed or not
 bool unrefinement_was_performed=false;
 bool refinement_was_performed=false;
 bool max_length_applied = false;

 //Loop over the number of polylines
 unsigned n_polyline = polygon_pt->npolyline();

 // Get face mesh representation of all polylines, possibly
 // with segments re-distributed to maintain an approximately
 // even sub-division of the polygon
 Vector<Mesh*> face_mesh_pt;
 get_face_mesh_representation(polygon_pt,face_mesh_pt);

 // Create vertices for the polylines by using the vertices
 // of the FaceElements
 Vector<double> vertex_coord(3); // zeta,x,y
 Vector<double> bound_left(1);
 Vector<double> bound_right(1);

 for(unsigned p=0;p<n_polyline;p++)
  {
   // Set of coordinates that will be placed on the boundary
   // Set entries are ordered on first  entry in vector which stores
   // the boundary coordinate so the vertices come out in order!
   std::set<Vector<double> > vertex_nodes;

   //Get the boundary id
   unsigned bound=polygon_pt->curve_section_pt(p)->boundary_id();

   // Loop over the face elements (ordered) and add their vertices
   unsigned n_face_element = face_mesh_pt[p]->nelement();

   //n_count = 0;
   for(unsigned e=0;e<n_face_element;++e)
    {
     FiniteElement* el_pt = face_mesh_pt[p]->finite_element_pt(e);
     unsigned n_node = el_pt->nnode();

     //Add the left-hand node to the set:

     // Boundary coordinate
     el_pt->node_pt(0)->get_coordinates_on_boundary(bound,bound_left);
     vertex_coord[0] = bound_left[0];

     // Actual coordinates
     for(unsigned i=0;i<2;i++)
      {
       vertex_coord[i+1] = el_pt->node_pt(0)->x(i);
      }
     vertex_nodes.insert(vertex_coord);

     //Add the right-hand nodes to the set:

     //Boundary coordinate
     el_pt->node_pt(n_node-1)->get_coordinates_on_boundary(bound,bound_right);
     vertex_coord[0] = bound_right[0];

     // Actual coordinates
     for(unsigned i=0;i<2;i++)
      {
       vertex_coord[i+1] = el_pt->node_pt(n_node-1)->x(i);
      }
     vertex_nodes.insert(vertex_coord);

    }

   // Now turn into vector for ease of handling...
   unsigned n_poly_vertex = vertex_nodes.size();
   Vector<Vector<double> > tmp_vector_vertex_node(n_poly_vertex);
   unsigned count=0;
   for(std::set<Vector<double> >::iterator it = vertex_nodes.begin();
     it!=vertex_nodes.end();++it)
    {
     tmp_vector_vertex_node[count].resize(3);
     tmp_vector_vertex_node[count][0] = (*it)[0];
     tmp_vector_vertex_node[count][1] = (*it)[1];
     tmp_vector_vertex_node[count][2] = (*it)[2];
     ++count;
    }

   // Size of the vector
   unsigned n_vertex=tmp_vector_vertex_node.size();

   // Tolerance below which the middle point can be deleted
   // (ratio of deflection to element length)
   double unrefinement_tolerance=
    polygon_pt->polyline_pt(p)->unrefinement_tolerance();

   //------------------------------------------------------
   // Unrefinement
   //------------------------------------------------------
   if (unrefinement_tolerance>0.0 && n_vertex>=3)
    {

     unrefinement_was_performed =
       unrefine_boundary(bound, tmp_vector_vertex_node,
			 unrefinement_tolerance, check_only);

     // In this case the "unrefinement_was_performed" variable
     // tell us if the update had been performed when calling
     // with check_oly=false
     if (check_only && unrefinement_was_performed)
      {
       // Cleanup (but only the elements -- the nodes still exist in
       // the bulk mesh!
       for(unsigned p=0;p<n_polyline;p++)
        {
         face_mesh_pt[p]->flush_node_storage();
         delete face_mesh_pt[p];
        }
       return true;
      }

    } // end of unrefinement

   //------------------------------------------------
   // Refinement
   //------------------------------------------------
   double refinement_tolerance=
    polygon_pt->polyline_pt(p)->refinement_tolerance();
   if (refinement_tolerance>0.0)
    {
     refinement_was_performed =
       refine_boundary(face_mesh_pt[p], tmp_vector_vertex_node,
         refinement_tolerance, check_only);

     // In this case the "refinement_was_performed" variable
     // tell us if the update had been performed when calling
     // with check_oly=false
     if (check_only && refinement_was_performed)
      {
       // Cleanup (but only the elements -- the nodes still exist in
       // the bulk mesh!
       for(unsigned p=0;p<n_polyline;p++)
        {
         face_mesh_pt[p]->flush_node_storage();
         delete face_mesh_pt[p];
        }
       return true;
      }

    } // end refinement

   //------------------------------------------------
   // Maximum length constrait
   //-----------------------------------------------
   double maximum_length = polygon_pt->polyline_pt(p)->maximum_length();
   if (maximum_length > 0.0)
    {
     max_length_applied = 
      apply_max_length_constraint(face_mesh_pt[p], 
                                  tmp_vector_vertex_node,
                                  maximum_length);
     
     // In this case the max length criteria was applied, check if 
     // check_only=false
     if (check_only && max_length_applied)
      {
       // Cleanup (but only the elements -- the nodes still exist in
       // the bulk mesh!
       for(unsigned p=0;p<n_polyline;p++)
        {
         face_mesh_pt[p]->flush_node_storage();
         delete face_mesh_pt[p];
        }
       return true;
      }

    }

   // For further processing the three-dimensional vector
   // has to be reduced to a two-dimensional vector
   n_vertex=tmp_vector_vertex_node.size();
   Vector<Vector<double> > vector_vertex_node(n_vertex);

   for(unsigned i=0;i<n_vertex;i++)
    {
     vector_vertex_node[i].resize(2);
     vector_vertex_node[i][0]=tmp_vector_vertex_node[i][1];
     vector_vertex_node[i][1]=tmp_vector_vertex_node[i][2];
    }

   if ( (p > 0) && !check_only )
    {
     //Final end point of previous line
     Vector<double> final_vertex_of_previous_segment;
     unsigned n_prev_vertex =
       polygon_pt->curve_section_pt(p-1)->nvertex();
     final_vertex_of_previous_segment =
       polygon_pt->polyline_pt(p-1)->
       vertex_coordinate(n_prev_vertex-1);

     unsigned prev_seg_boundary_id =
       polygon_pt->curve_section_pt(p-1)->boundary_id();

     //Find the error between the final vertex of the previous
     //line and the first vertex of the current line
     double error = 0.0;
     for(unsigned i=0;i<2;i++)
      {
       const double dist =
         final_vertex_of_previous_segment[i] -
         (*vector_vertex_node.begin())[i];
       error += dist*dist;
      }
     error = sqrt(error);

     //If the error is bigger than the tolerance then
     //we probably need to reverse, but better check
     if(error > ToleranceForVertexMismatchInPolygons::Tolerable_error)
      {
       //Find the error between the final vertex of the previous
       //line and the last vertex of the current line
       double rev_error = 0.0;
       for(unsigned i=0;i<2;i++)
        {
         const double dist =
           final_vertex_of_previous_segment[i] -
           (*--vector_vertex_node.end())[i];
         rev_error += dist*dist;
        }
       rev_error = sqrt(rev_error);

       if(rev_error >
         ToleranceForVertexMismatchInPolygons::Tolerable_error)
        {
         // It could be possible that the first segment be reversed and we
         // did not notice it because this check does not apply for the
         // first segment. We can verify if the first segment is reversed
         // by using the vertex number 1
         if (p == 1)
          {

           //Initial end point of previous line
           Vector<double> initial_vertex_of_previous_segment;

           initial_vertex_of_previous_segment =
             polygon_pt->polyline_pt(p-1)->
             vertex_coordinate(0);

           unsigned prev_seg_boundary_id =
             polygon_pt->curve_section_pt(p-1)->boundary_id();

           //Find the error between the initial vertex of the previous
           //line and the first vertex of the current line
           double error = 0.0;
           for(unsigned i=0;i<2;i++)
            {
             const double dist =
               initial_vertex_of_previous_segment[i] -
               (*vector_vertex_node.begin())[i];
             error += dist*dist;
            }
           error = sqrt(error); // Reversed only the previous one

           //If the error is bigger than the tolerance then
           //we probably need to reverse, but better check
           if(error > ToleranceForVertexMismatchInPolygons::Tolerable_error)
            {
             //Find the error between the final vertex of the previous
             //line and the last vertex of the current line
             double rev_error = 0.0;
             for(unsigned i=0;i<2;i++)
              {
               const double dist =
                 initial_vertex_of_previous_segment[i] -
                 (*--vector_vertex_node.end())[i];
               rev_error += dist*dist;
              }
             rev_error = sqrt(rev_error); // Reversed both the current one and
             // the previous one

             if (rev_error >
               ToleranceForVertexMismatchInPolygons::Tolerable_error)
              {

               std::ostringstream error_stream;
               error_stream <<
                 "The distance between the first node of the current\n" <<
                 "line segment (boundary " << bound << ") and either end of "
                 << "the previous line segment\n"
                 << "(boundary " << prev_seg_boundary_id << ") is bigger than "
                 << "the desired tolerance " <<
                 ToleranceForVertexMismatchInPolygons::Tolerable_error << ".\n"
                 << "This suggests that the polylines defining the polygonal\n"
                 << "representation are not properly ordered.\n"
                 << "Fail on last vertex of polyline: ("
                 << prev_seg_boundary_id<< ") and\nfirst vertex of polyline ("
                 << bound << ").\nThis should have failed when first trying to "
                 << "construct the\npolygon.\n";
               throw OomphLibError(error_stream.str(),
                                   OOMPH_CURRENT_FUNCTION,
                 OOMPH_EXCEPTION_LOCATION);

              }
             else
              {
               // Reverse both
               //Reverse the current vector to line up with the previous one
               std::reverse(vector_vertex_node.begin(),vector_vertex_node.end());
               polygon_pt->polyline_pt(p-1)->reverse();
              }

            }
           else
            {
             // Reverse the previous one
             polygon_pt->polyline_pt(p-1)->reverse();
            }

          } // if p == 1
         else
          {

           std::ostringstream error_stream;
           error_stream <<
             "The distance between the first node of the current\n" <<
             "line segment (boundary " << bound << ") and either end of "
             << "the previous line segment\n"
             << "(boundary " << prev_seg_boundary_id << ") is bigger than the "
             << "desired tolerance " <<
             ToleranceForVertexMismatchInPolygons::Tolerable_error << ".\n"
             <<
             "This suggests that the polylines defining the polygonal\n" <<
             "representation are not properly ordered.\n"
             << "Fail on last vertex of polyline: (" << prev_seg_boundary_id
             << ") and\nfirst vertex of polyline (" << bound << ").\n"
             << "This should have failed when first trying to construct the\n"
             << "polygon.\n";
           throw OomphLibError(error_stream.str(),
                               OOMPH_CURRENT_FUNCTION,
             OOMPH_EXCEPTION_LOCATION);
          }

        }
       else
        {
         //Reverse the current vector to line up with the previous one
         std::reverse(vector_vertex_node.begin(),vector_vertex_node.end());
        }

      }
    }

   if(!check_only)
    {
     //Now update the polyline according to the new vertices but
     //first check if the object is allowed to delete the representation
     //or if it should be done by other object

     bool delete_it_on_destructor = false;

     std::set<TriangleMeshCurveSection*>::iterator it =
       this->Free_curve_section_pt.find(polygon_pt->curve_section_pt(p));

     if (it!=this->Free_curve_section_pt.end())
      {
       this->Free_curve_section_pt.erase(it);
       delete polygon_pt->curve_section_pt(p);
       delete_it_on_destructor = true;
      }

     polygon_pt->curve_section_pt(p) =
      new TriangleMeshPolyLine(vector_vertex_node,bound);
     
     // Establish refinement and unrefinement tolerance
     polygon_pt->curve_section_pt(p)->set_unrefinement_tolerance(
      unrefinement_tolerance);
     polygon_pt->curve_section_pt(p)->set_refinement_tolerance(
      refinement_tolerance);
     
     // Establish the maximum length constraint
     polygon_pt->curve_section_pt(p)->set_maximum_length(maximum_length);
     
     // Update the Boundary - Polyline map
     this->Boundary_curve_section_pt[bound] = polygon_pt->curve_section_pt(p);
     
     if (delete_it_on_destructor)
      {
       this->Free_curve_section_pt.insert(polygon_pt->curve_section_pt(p));
      }

    }
  } // n_polylines

 // Cleanup (but only the elements -- the nodes still exist in
 // the bulk mesh!
 for(unsigned p=0;p<n_polyline;p++)
  {
   face_mesh_pt[p]->flush_node_storage();
   delete face_mesh_pt[p];
  }

 if(check_only)
  {
   // if we end up all the way down here, no update of the internal boundaries
   // is necessary (in case we only check)
   return false;
  }
 else
  {
   // if we not only check, but actually perform the update and end up
   // all the way down here then we indicate whether an update was performed
   // or not
   return (unrefinement_was_performed || refinement_was_performed || max_length_applied);
  }

  }

//=========================================================================
/// \short Helper function that updates the input open curve by using
/// end-points of elements from FaceMesh(es) that are constructed for the
/// boundaries associated with the polylines. Optional boolean is used to
/// run it as test only (if true is specified as input) in which case the
/// polylines are not actually modified. Returned boolean indicates if
/// polylines were (or would have been -- if called with check_only=false)
/// changed.
//=========================================================================
template<class ELEMENT>
bool RefineableTriangleMesh<ELEMENT>::update_open_curve_using_face_mesh(
  TriangleMeshOpenCurve* open_polyline_pt,
  const bool& check_only)
  {
  // Boolean that indicates whether an actual update of the polylines
  // were performed or not
  bool unrefinement_was_performed=false;
  bool refinement_was_performed=false;
  bool max_length_applied = false;

  //Loop over the number of polylines
  unsigned n_polyline = open_polyline_pt->ncurve_section();

  // Get face mesh representation of all polylines, possibly
  // with segments re-distributed to maintain an approximately
  // even sub-division of the polygon
  Vector<Mesh*> face_mesh_pt;
  get_face_mesh_representation(open_polyline_pt, face_mesh_pt);

  // Create vertices for the polylines by using the vertices
  // of the FaceElements
  Vector<double> vertex_coord(3); // zeta,x,y
  Vector<double> bound_left(1);
  Vector<double> bound_right(1);

  for(unsigned p=0;p<n_polyline;p++)
   {
    // Set of coordinates that will be placed on the boundary
    // Set entries are ordered on first  entry in vector which stores
    // the boundary coordinate so the vertices come out in order!
    std::set<Vector<double> > vertex_nodes;

    //Get the boundary id
    unsigned bound=open_polyline_pt->curve_section_pt(p)->boundary_id();

    // Loop over the face elements (ordered) and add their vertices
    unsigned n_face_element = face_mesh_pt[p]->nelement();

    //n_count = 0;
    for(unsigned e=0;e<n_face_element;++e)
     {
      FiniteElement* el_pt = face_mesh_pt[p]->finite_element_pt(e);
      unsigned n_node = el_pt->nnode();

      //Add the left-hand node to the set:

      // Boundary coordinate
      el_pt->node_pt(0)->get_coordinates_on_boundary(bound,bound_left);
      vertex_coord[0] = bound_left[0];

      // Actual coordinates
      for(unsigned i=0;i<2;i++)
       {
        vertex_coord[i+1] = el_pt->node_pt(0)->x(i);
       }
      vertex_nodes.insert(vertex_coord);

      //Add the right-hand nodes to the set:

      //Boundary coordinate
      el_pt->node_pt(n_node-1)->get_coordinates_on_boundary(bound,bound_right);
      vertex_coord[0] = bound_right[0];

      // Actual coordinates
      for(unsigned i=0;i<2;i++)
       {
        vertex_coord[i+1] = el_pt->node_pt(n_node-1)->x(i);
       }
      vertex_nodes.insert(vertex_coord);

     }

    // Now turn into vector for ease of handling...
    unsigned n_poly_vertex = vertex_nodes.size();
    Vector<Vector<double> > tmp_vector_vertex_node(n_poly_vertex);
    unsigned count=0;
    for(std::set<Vector<double> >::iterator it = vertex_nodes.begin();
      it!=vertex_nodes.end();++it)
     {
      tmp_vector_vertex_node[count].resize(3);
      tmp_vector_vertex_node[count][0] = (*it)[0];
      tmp_vector_vertex_node[count][1] = (*it)[1];
      tmp_vector_vertex_node[count][2] = (*it)[2];
      ++count;
     }

    // Size of the vector
    unsigned n_vertex=tmp_vector_vertex_node.size();

    // Tolerance below which the middle point can be deleted
    // (ratio of deflection to element length)
    double unrefinement_tolerance=
      open_polyline_pt->polyline_pt(p)->unrefinement_tolerance();

    //------------------------------------------------------
    // Unrefinement
    //------------------------------------------------------
    if (unrefinement_tolerance>0.0 && n_vertex>=3)
     {
      unrefinement_was_performed =
        unrefine_boundary(bound, tmp_vector_vertex_node,
			  unrefinement_tolerance, check_only);

      // In this case the unrefinement_was_performed variable actually
      // tell us if the update had been performed when calling
      // with check_only=false
      if (check_only && unrefinement_was_performed)
       {
        // Cleanup (but only the elements -- the nodes still exist in
        // the bulk mesh!
        for(unsigned p=0;p<n_polyline;p++)
         {
          face_mesh_pt[p]->flush_node_storage();
          delete face_mesh_pt[p];
         }
        return true;
       }

     } // end of unrefinement

    //------------------------------------------------
    /// Refinement
    //------------------------------------------------
    double refinement_tolerance=
      open_polyline_pt->polyline_pt(p)->refinement_tolerance();
    if (refinement_tolerance>0.0)
     {
      refinement_was_performed =
        refine_boundary(face_mesh_pt[p], tmp_vector_vertex_node,
          refinement_tolerance, check_only);

      // In this case the unrefinement_was_performed variable actually
      // tell us if the update had been performed when calling
      // with check_only=false
      if (check_only && refinement_was_performed)
       {
        // Cleanup (but only the elements -- the nodes still exist in
        // the bulk mesh!
        for(unsigned p=0;p<n_polyline;p++)
         {
          face_mesh_pt[p]->flush_node_storage();
          delete face_mesh_pt[p];
         }
        return true;
       }

     } // end refinement

    //------------------------------------------------
    // Maximum length constrait
    //-----------------------------------------------
    double maximum_length = open_polyline_pt->polyline_pt(p)->maximum_length();
    if (maximum_length > 0.0)
      {
       bool max_length_applied = false;
       max_length_applied = 
        apply_max_length_constraint(face_mesh_pt[p], 
                                    tmp_vector_vertex_node,
                                    maximum_length);
       
       // In this case the max length criteria was applied, check if 
       // check_only=false
       if (check_only && max_length_applied)
        {
         // Cleanup (but only the elements -- the nodes still exist in
         // the bulk mesh!
         for(unsigned p=0;p<n_polyline;p++)
          {
           face_mesh_pt[p]->flush_node_storage();
           delete face_mesh_pt[p];
          }
         return true;
        }
       
      }
    
    // For further processing the three-dimensional vector
    // has to be reduced to a two-dimensional vector
    n_vertex=tmp_vector_vertex_node.size();
    Vector<Vector<double> > vector_vertex_node(n_vertex);

    for(unsigned i=0;i<n_vertex;i++)
     {
      vector_vertex_node[i].resize(2);
      vector_vertex_node[i][0]=tmp_vector_vertex_node[i][1];
      vector_vertex_node[i][1]=tmp_vector_vertex_node[i][2];
     }

    //Check whether the segments are continguous (first vertex of this
    //segment is equal to last vertex of previous segment).
    //If not, we should reverse the order of the current segment.
    //This check only applies for segments other than the first.
    //We only bother with this check, if we actually perform an update
    //of the polyline, i.e. if it's not only a check
    if( (p > 0) && !check_only )
     {
      //Final end point of previous line
      Vector<double> final_vertex_of_previous_segment;
      open_polyline_pt->polyline_pt(p-1)->
       final_vertex_coordinate(final_vertex_of_previous_segment);
      
      unsigned prev_seg_boundary_id =
       open_polyline_pt->curve_section_pt(p-1)->boundary_id();
      
      //Find the error between the final vertex of the previous
      //line and the first vertex of the current line
      double error = 0.0;
      for(unsigned i=0;i<2;i++)
       {
        const double dist =
         final_vertex_of_previous_segment[i] -
         (*vector_vertex_node.begin())[i];
        error += dist*dist;
       }
      error = sqrt(error);
      
      //If the error is bigger than the tolerance then
      //we probably need to reverse, but better check
      if(error > ToleranceForVertexMismatchInPolygons::Tolerable_error)
       {
        //Find the error between the final vertex of the previous
        //line and the first vertex of the current line
        error = 0.0;
        for(unsigned i=0;i<2;i++)
         {
          const double dist =
           final_vertex_of_previous_segment[i] -
           (*--vector_vertex_node.end())[i];
          error += dist*dist;
         }
        error = sqrt(error);
        
        if (error > ToleranceForVertexMismatchInPolygons::Tolerable_error)
         {
          // If no found it is possible that the previous polyline be reversed
          // Check for that case
          //Initial point of previous line
          Vector<double> initial_vertex_of_previous_segment;
          open_polyline_pt->polyline_pt(p-1)->
           initial_vertex_coordinate(initial_vertex_of_previous_segment);
          
          //Find the error between the initial vertex of the previous
          //line and the first vertex of the current line
          error = 0.0;
          for(unsigned i=0;i<2;i++)
           {
            const double dist =
             initial_vertex_of_previous_segment[i] -
             (*vector_vertex_node.begin())[i];
            error += dist*dist;
           }
          error = sqrt(error);
          
          //If the error is bigger than the tolerance then
          //we probably need to reverse, but better check
          if(error > ToleranceForVertexMismatchInPolygons::Tolerable_error)
           {
            //Find the error between the final vertex of the previous
            //line and the first vertex of the current line
            error = 0.0;
            for(unsigned i=0;i<2;i++)
             {
              const double dist =
               initial_vertex_of_previous_segment[i] -
               (*--vector_vertex_node.end())[i];
              error += dist*dist;
             }
            error = sqrt(error);
            
            if(error >
               ToleranceForVertexMismatchInPolygons::Tolerable_error)
             {
              std::ostringstream error_stream;
              error_stream 
               <<"The distance between the first node of the current\n"
               <<"line segment (boundary " << bound
               <<") and either end of the previous line segment\n"
               <<"(boundary " << prev_seg_boundary_id << ") is bigger than "
               <<"the desired tolerance " <<
               ToleranceForVertexMismatchInPolygons::Tolerable_error << ".\n"
               <<"This suggests that the polylines defining the open curve\n"
               <<"representation are not properly ordered.\n"
               <<"Fail on last vertex of polyline: (" << prev_seg_boundary_id
               <<") and\nfirst vertex of polyline (" << bound << ").\n"
               <<"This should have failed when first trying to construct\n"
               <<"the open curve.\n";
              throw OomphLibError(error_stream.str(),
                                  OOMPH_CURRENT_FUNCTION,
               OOMPH_EXCEPTION_LOCATION);
             }
            else // We have to reverse both
             {
              // First reverse the previous polyline
              open_polyline_pt->polyline_pt(p-1)->reverse();
              // Then reverse the current polyline
              std::reverse(vector_vertex_node.begin(),
                           vector_vertex_node.end());
             }
           }
          else
           {
            // Reverse the previous polyline only
            open_polyline_pt->polyline_pt(p-1)->reverse();
           }
         }
        else
         {
          //Reverse the current vector to line up with the previous one
          std::reverse(vector_vertex_node.begin(),vector_vertex_node.end());
         }

       }

     } // if p > 0

    if(!check_only)
     {
      // Now update the polyline according to the new vertices
      // The new one representation
      TriangleMeshPolyLine *tmp_polyline =
        new TriangleMeshPolyLine(vector_vertex_node, bound);

      // Create a temporal "curve section" version of the recently created
      // polyline
      TriangleMeshCurveSection *tmp_curve_section = tmp_polyline;
      
      // Copy the unrefinement and refinement information
      tmp_polyline->set_unrefinement_tolerance(
       unrefinement_tolerance);
      tmp_polyline->set_refinement_tolerance(
       refinement_tolerance);
      
      // Establish the maximum length constraint
      tmp_polyline->set_maximum_length(maximum_length);
      
      // *****************************************************************
      // We need to pass the connection information from the
      // old polyline to the new one (if the polyline is connected
      // or not. After this we need to restore the vertices number used for
      // connections.
      this->compute_connection_information(
       open_polyline_pt->polyline_pt(p), tmp_curve_section);
      
      // Re-store the connection vertices
      // Since we have added or eliminated nodes (vertices) on the
      // mesh/domain then we need to restore the connections (it means,
      // compute the new vertices numbers for connections since they could
      // be moved). Call it with the "first_try" option enabled
      const bool restored_connections = 
       this->restore_connections_on_internal_boundary(tmp_polyline);
      
      // Check if it was possible to restore the connections, if that
      // is not the case then invert the connection information of the
      // polyline and try to restore the connections again
      if (!restored_connections)
       {
        // Invert the connection information of the polyline
        const bool invert_connection_information = true;
        this->compute_connection_information(
         open_polyline_pt->polyline_pt(p), tmp_curve_section, 
         invert_connection_information);
        
        // Last chage to restore the connections
        const bool first_try = false;
        this->restore_connections_on_internal_boundary(tmp_polyline,first_try);
       }
      
      std::set<TriangleMeshCurveSection*>::iterator it =
        this->Free_curve_section_pt.find(open_polyline_pt->curve_section_pt(p));
      
      bool delete_it_on_destructor = false;

      if (it!=this->Free_curve_section_pt.end())
       {
        // Free previous representation only if you created
        this->Free_curve_section_pt.erase(it);
        delete open_polyline_pt->curve_section_pt(p);
        delete_it_on_destructor = true;
       }

      // *****************************************************************
      // Copying the new representation
      open_polyline_pt->curve_section_pt(p) = tmp_polyline;

      // Update the Boundary <--> PolyLine map
      this->Boundary_curve_section_pt[bound] =
        open_polyline_pt->curve_section_pt(p);

      if (delete_it_on_destructor)
       {
        this->Free_curve_section_pt.insert(
          open_polyline_pt->curve_section_pt(p));
       }

     }

   } // n_polylines

  // Cleanup (but only the elements -- the nodes still exist in
  // the bulk mesh!
  for(unsigned p=0;p<n_polyline;p++)
   {
    face_mesh_pt[p]->flush_node_storage();
    delete face_mesh_pt[p];
   }

  if(check_only)
   {
    // if we end up all the way down here, no update of the internal boundaries
    // is necessary (in case we only check)
    return false;
   }
  else
   {
    // if we not only check, but actually perform the update and end up
    // all the way down here then we indicate whether an update was performed
    // or not
    return (unrefinement_was_performed || refinement_was_performed || max_length_applied);
   }

  }

//=========================================================================
/// \short Helper function that performs the unrefinement process
/// on the specified boundary by using the provided vertices
/// representation. Optional boolean is used to run it as test only (if
/// true is specified as input) in which case vertex coordinates aren't
/// actually modified. Returned boolean indicates if polyline was (or
/// would have been -- if called with check_only=false) changed.
//=========================================================================
template<class ELEMENT>
bool RefineableTriangleMesh<ELEMENT>::
unrefine_boundary(const unsigned &b,
  Vector<Vector<double> > &vector_bnd_vertices,
  double &unrefinement_tolerance,
  const bool &check_only)
 {
  // *************************************************************
  // Compute the vector of vertices not allowed for deletion!!!
  Vector<Vector<double> > no_del_vertex;

  // Get the set of nodes on this boundary that should not be
  // deleted
  if (!this->Boundary_connections_pt[b].empty())
   {
    // For easy to use
    std::set<Vector<double> > tmp_set =
    this->Boundary_connections_pt[b];

    unsigned counter = 0;
    no_del_vertex.resize(tmp_set.size());

    // Loop over all the elements on the set and get the
    // vertex coordinates of the nodes
    for (std::set<Vector<double> >::iterator it = tmp_set.begin();
      it!=tmp_set.end();++it)
     {
      no_del_vertex[counter].resize(2);
      no_del_vertex[counter][0] = (*it)[0];
      no_del_vertex[counter][1] = (*it)[1];
      counter++;
     }
   }
  // *************************************************************

  // Boolean that indicates whether an actual update of the vertex
  // coordinates was performed or not
  bool unrefinement_was_performed=false;

  unsigned n_vertex = vector_bnd_vertices.size();

  // Initialise counter that indicates at which vertex we're currently
  // considering for deletion
  unsigned counter=1;

  // Loop over the nodes; start with the second one and increment by two
  // this way a "pack" of three nodes will be considered for calculation:
  // the middle-node (which is to be deleted or not) and the adjacent
  // nodes
  for(unsigned i=1;i<=n_vertex-2;i+=2)
   {
    // Maths from http://www.cgafaq.info/wiki/Circle_Through_Three_Points
    double a_x=vector_bnd_vertices[i-1][1];
    double a_y=vector_bnd_vertices[i-1][2];
    double b_x=vector_bnd_vertices[i][1];
    double b_y=vector_bnd_vertices[i][2];
    double c_x=vector_bnd_vertices[i+1][1];
    double c_y=vector_bnd_vertices[i+1][2];

    double a=b_x-a_x;
    double b=b_y-a_y;
    double c=c_x-a_x;
    double d=c_y-a_y;

    double e=a*(a_x+b_x)+b*(a_y+b_y);
    double f=c*(a_x+c_x)+d*(a_y+c_y);

    double g=2.0*(a*(c_y-b_y)-b*(c_x-b_x));

    bool do_it=false;
    if (std::fabs(g)<1.0e-14)
     {
      do_it=true;
      if(check_only) {return true;}
     }
    else
     {
      double p_x=(d*e-b*f)/g;
      double p_y=(a*f-c*e)/g;

      double r=sqrt(pow((a_x-p_x),2)+pow((a_y-p_y),2));

      double rhalfca_x=0.5*(a_x-c_x);
      double rhalfca_y=0.5*(a_y-c_y);

      double halfca_squared=pow(rhalfca_x,2)+pow(rhalfca_y,2);

      double sticky_out_bit=r-sqrt(std::fabs((r*r) - halfca_squared));

      // If sticky out bit divided by distance between end nodes
      // is less than tolerance the boundary is so flat that we
      // can safely kill the node
      if ((sticky_out_bit/(2.0*sqrt(halfca_squared)))<
        unrefinement_tolerance)
       {
        do_it=true;
        if(check_only) {return true;}
       }
     }

    if (do_it)
     {
      // Last chance for not delete the node
      // Verify if it is not a node that should not be delete
      unsigned size = no_del_vertex.size();
      for (unsigned i = 0; i < size; i++)
       {
        // Compute the distance between the proposed node to delete
        // and the ones that should not be deleted
        double x = no_del_vertex[i][0];
        double y = no_del_vertex[i][1];
        double error = (b_x - x)*(b_x - x) + (b_y - y)*(b_y - y);
        error = sqrt(error);

        if(error < ToleranceForVertexMismatchInPolygons::Tolerable_error)
         {
          do_it = false;
          break;
         }

       }
     }

    // Remove node?
    if (do_it)
     {
      vector_bnd_vertices[i].resize(0);
     }

    // Increase the counter, that indicates the number of the
    // next middle node
    counter+=2;
   }

  // coming out of here the value of counter is the index of the
  // last node on the polyline counter=n_vertex-1 (in case of an
  // even number of nodes) or counter has the value of the number
  // of nodes on the polyline counter=n_vertex (in case of an odd
  // number of nodes

  // Special treatment for the end of the polyline:
  // If the number of nodes is even, then the previous loop stopped
  // at the last but second node, i.e. the current value of counter
  // is the index of the last node. If that's the case, the last but
  // one node needs to be treated separately
  if( (counter)==(n_vertex-1) )
   {
    // Set the last but one node as middle node
    unsigned i=vector_bnd_vertices.size()-2;

    // Index of the current! last but second node (considering any
    // previous deletion)
    unsigned n=0;

    if(vector_bnd_vertices[counter-2].size()!=0)
     {
      // if the initial last but second node does still exist then
      // this one is obviously also the current last but second one
      n=counter-2;
     }
    else
     {
      // if the initial last but second node was deleted then the
      // initial last but third node is the current last but second
      // node
      n=counter-3;
     }

    // CODE DUPLICATION -- CAN'T BE BOTHERED TO WRITE A SEPARATE
    // FUNCTION FOR THIS; PROBABLY WORTH DOING IF/WHEN THERE'S
    // A MISTAKE IN ANY OF THIS AND IT NEEDS TO BE FIXED...

    // Maths from http://www.cgafaq.info/wiki/Circle_Through_Three_Points
    double a_x=vector_bnd_vertices[n][1];
    double a_y=vector_bnd_vertices[n][2];
    double b_x=vector_bnd_vertices[i][1];
    double b_y=vector_bnd_vertices[i][2];
    double c_x=vector_bnd_vertices[i+1][1];
    double c_y=vector_bnd_vertices[i+1][2];

    double a=b_x-a_x;
    double b=b_y-a_y;
    double c=c_x-a_x;
    double d=c_y-a_y;

    double e=a*(a_x+b_x)+b*(a_y+b_y);
    double f=c*(a_x+c_x)+d*(a_y+c_y);

    double g=2.0*(a*(c_y-b_y)-b*(c_x-b_x));

    bool do_it=false;
    if (std::fabs(g)<1.0e-14)
     {
      do_it=true;
      if(check_only) {return true;}
     }
    else
     {
      double p_x=(d*e-b*f)/g;
      double p_y=(a*f-c*e)/g;

      double r=sqrt(pow((a_x-p_x),2)+pow((a_y-p_y),2));

      double rhalfca_x=0.5*(a_x-c_x);
      double rhalfca_y=0.5*(a_y-c_y);

      double halfca_squared=pow(rhalfca_x,2)+pow(rhalfca_y,2);

      double sticky_out_bit=r-sqrt(std::fabs((r*r) - halfca_squared));

      // If sticky out bit divided by distance between end nodes
      // is less than tolerance the boundary is so flat that we
      // can safely kill the node
      if ((sticky_out_bit/(2.0*sqrt(halfca_squared)))<
        unrefinement_tolerance)
       {
        do_it=true;
        if(check_only) {return true;}
       }
     }

    // Last chance for not delete the node
    // Verify if it is not a node that should not be delete
    if (do_it)
     {
      unsigned size = no_del_vertex.size();
      for (unsigned i = 0; i < size; i++)
       {
        // Compute the distance between the proposed node to delete
        // and the ones that should not be deleted
        double x = no_del_vertex[i][0];
        double y = no_del_vertex[i][1];
        double error = (b_x - x)*(b_x - x) + (b_y - y)*(b_y - y);
        error = sqrt(error);

        if(error <
          ToleranceForVertexMismatchInPolygons::Tolerable_error)
         {
          do_it = false;
          break;
         }

       }
     }

    // Remove node?
    if (do_it)
     {
      vector_bnd_vertices[i].resize(0);
     }
   }

  // Create another vector, which will only contain entries of
  // nodes that still exist
  Vector<Vector<double> > compact_vector;
  compact_vector.reserve(n_vertex);
  for (unsigned i=0;i<n_vertex;i++)
   {
    // If the entry was not deleted include it in the new vector
    if (vector_bnd_vertices[i].size()!=0)
     {
      compact_vector.push_back(vector_bnd_vertices[i]);
     }
   }

  /// Get the size of the vector that now includes all remaining nodes
  n_vertex =compact_vector.size();

  // If the size of the vector containing the remaining nodes is
  // different from the size of the vector before the unrefinement
  // routine (with the original nodes)
  // then the polyline was obviously updated
  if( n_vertex != vector_bnd_vertices.size() )
   {
    unrefinement_was_performed=true;
   }

  /// Copy back
  vector_bnd_vertices.resize(n_vertex);
  for(unsigned i=0;i<n_vertex;i++)
   {
    vector_bnd_vertices[i].resize(3);
    vector_bnd_vertices[i][0]=compact_vector[i][0];
    vector_bnd_vertices[i][1]=compact_vector[i][1];
    vector_bnd_vertices[i][2]=compact_vector[i][2];
   }

  return unrefinement_was_performed;

 }

//=========================================================================
/// \short Helper function that performs the refinement process
/// on the specified boundary by using the provided vertices
/// representation. Optional boolean is used to run it as test only (if
/// true is specified as input) in which case vertex coordinates aren't
/// actually modified. Returned boolean indicates if polyline was (or
/// would have been -- if called with check_only=false) changed.
//=========================================================================
template<class ELEMENT>
bool RefineableTriangleMesh<ELEMENT>::
refine_boundary(Mesh* face_mesh_pt,
  Vector<Vector<double> > &vector_bnd_vertices,
  double &refinement_tolerance,
  const bool &check_only)
 {

  // Boolean that indicates whether an actual update of the vertex
  // coordinates was performed or not
  bool refinement_was_performed=false;

  // Create a geometric object from the mesh to represent
  //the curvilinear boundary

  MeshAsGeomObject* mesh_geom_obj_pt =
    new MeshAsGeomObject(face_mesh_pt);

  // Get the total number of current vertices
  unsigned n_vertex=vector_bnd_vertices.size();

  // Create a new (temporary) vector for the nodes, so
  // that new nodes can be stored
  Vector<Vector<double> > extended_vector;

  // Reserve memory space for twice the number of already
  // existing nodes (worst case)
  extended_vector.reserve(2*n_vertex);

  // Loop over the nodes until the last but one node
  for(unsigned inod=0;inod<n_vertex-1;inod++)
   {
    // Get local coordinate of "left" node
    double zeta_left=vector_bnd_vertices[inod][0];

    // Get position vector of "left" node
    Vector<double> R_left(2);
    for(unsigned i=0;i<2;i++)
     {
      R_left[i]=vector_bnd_vertices[inod][i+1];
     }

    // Get local coordinate of "right" node
    double zeta_right=vector_bnd_vertices[inod+1][0];

    // Get position vector of "right" node
    Vector<double> R_right(2);
    for(unsigned i=0;i<2;i++)
     {
      R_right[i]=vector_bnd_vertices[inod+1][i+1];
     }

    // Get the boundary coordinate of the midpoint
    Vector<double> zeta_mid(1);
    zeta_mid[0]=0.5*(zeta_left+zeta_right);

    // Get the position vector of the midpoint on the
    // curvilinear boundary
    Vector<double> R_mid(2);
    mesh_geom_obj_pt->position(zeta_mid,R_mid);

    // Get the position vector of the midpoint on the straight
    // line connecting "left" and "right" node
    Vector<double> R_mid_polygon(2);
    for(unsigned i=0;i<2;i++)
     {
      R_mid_polygon[i]=0.5*(R_right[i]+R_left[i]);
     }

    // Calculate the distance between the midpoint on the curvilinear
    // boundary and the midpoint on the straight line
    double distance=sqrt((R_mid[0]-R_mid_polygon[0])*
      (R_mid[0]-R_mid_polygon[0])+
      (R_mid[1]-R_mid_polygon[1])*
      (R_mid[1]-R_mid_polygon[1]));

    // Calculating the length of the straight line
    double length=sqrt((R_right[0]-R_left[0])*(R_right[0]-R_left[0])+
      (R_right[1]-R_left[1])*(R_right[1]-R_left[1]));

    // If the ratio of distance between the midpoints to the length
    // of the straight line is larger than the tolerance
    // specified for the criterion when points can be deleted,
    // create a new node and add it to the (temporary) vector
    if((distance/length) > refinement_tolerance)
     {

      if(check_only)
       {
        // Delete the allocated memory for the geometric object
        // that represents the curvilinear boundary
        delete mesh_geom_obj_pt;
        return true;
       }

      Vector<double> new_node(3);
      new_node[0]=zeta_mid[0];
      new_node[1]=R_mid[0];
      new_node[2]=R_mid[1];

      // Include the "left" node in the new "temporary" vector
      extended_vector.push_back(vector_bnd_vertices[inod]);

      // Include the new node as well
      extended_vector.push_back(new_node);

     }
    else
     {
      // Include the "left" node in the new "temporary" vector
      // and move on to the next node
      extended_vector.push_back(vector_bnd_vertices[inod]);
     }
   } // end of loop over nodes

  // Add the last node to the vector
  extended_vector.push_back(vector_bnd_vertices[n_vertex-1]);

  /// Get the size of the vector that now includes all added nodes
  n_vertex=extended_vector.size();

  // If the size of the vector including the added nodes is
  // different from the size of the vector before the refinement
  // routine then the polyline was obviously updated
  if( n_vertex != vector_bnd_vertices.size() )
   {
    refinement_was_performed=true;
   }

  // Copy across
  vector_bnd_vertices.resize(n_vertex);
  for(unsigned i=0;i<n_vertex;i++)
   {
    vector_bnd_vertices[i].resize(3);
    vector_bnd_vertices[i][0]=extended_vector[i][0];
    vector_bnd_vertices[i][1]=extended_vector[i][1];
    vector_bnd_vertices[i][2]=extended_vector[i][2];
   }

  // Delete the allocated memory for the geometric object
  // that represents the curvilinear boundary
  delete mesh_geom_obj_pt;

  return refinement_was_performed;

 }

 template<class ELEMENT>
 bool RefineableTriangleMesh<ELEMENT>::
 apply_max_length_constraint(Mesh* face_mesh_pt,
                             Vector<Vector<double> > &vector_bnd_vertices,
                             double &max_length_constraint)
 {

  // Boolean that indicates whether an actual update of the vertex
  // coordinates was performed or not
  bool max_length_applied=false;

  // Create a geometric object from the mesh to represent
  //the curvilinear boundary
  MeshAsGeomObject* mesh_geom_obj_pt =
   new MeshAsGeomObject(face_mesh_pt);

  // Get the total number of current vertices
  unsigned n_vertex=vector_bnd_vertices.size();

  // Create a new (temporary) vector for the nodes, so
  // that new nodes can be stored
  Vector<Vector<double> > extended_vector;

  // Reserve memory space for twice the number of already
  // existing nodes (worst case)
  //extended_vector.reserve(2*n_vertex);

  // Loop over the nodes until the last but one node
  for(unsigned inod=0;inod<n_vertex-1;inod++)
   {
    // Get local coordinate of "left" node
    double zeta_left=vector_bnd_vertices[inod][0];

    // Get position vector of "left" node
    Vector<double> R_left(2);
    for(unsigned i=0;i<2;i++)
     {
      R_left[i]=vector_bnd_vertices[inod][i+1];
     }

    // Get local coordinate of "right" node
    double zeta_right=vector_bnd_vertices[inod+1][0];

    // Get position vector of "right" node
    Vector<double> R_right(2);
    for(unsigned i=0;i<2;i++)
     {
      R_right[i]=vector_bnd_vertices[inod+1][i+1];
     }

    // Include the "left" node in the new "temporary" vector
    extended_vector.push_back(vector_bnd_vertices[inod]);

    // Check whether the current distance between the left and right node 
    // is longer than the specified constraint or not
    double length=std::fabs(zeta_right-zeta_left);

    // Do we need to introduce new nodes?
    if (length > max_length_constraint)
     {
      double n_pts = length/max_length_constraint;
      // We only want the integer part
      unsigned n_points = static_cast<unsigned>(n_pts);
      double zeta_increment = (zeta_right-zeta_left)/((double)n_points+1);
       
      Vector<double> zeta(1);
      // Create the n_points+1 points inside the segment
      for(unsigned s=1;s<n_points+1;s++)
       {
        // Get the coordinates
        zeta[0]= zeta_left + zeta_increment*double(s);
        Vector<double> vertex(2);
        mesh_geom_obj_pt->position(zeta, vertex);

        // Create the new node
        Vector<double> new_node(3);
        new_node[0]=zeta[0];
        new_node[1]=vertex[0];
        new_node[2]=vertex[1];

        // Include the new node
        extended_vector.push_back(new_node);
       }
     }
   }

  // Add the last node to the vector
  extended_vector.push_back(vector_bnd_vertices[n_vertex-1]);

  /// Get the size of the vector that now includes all added nodes
  n_vertex=extended_vector.size();

  // If the size of the vector including the added nodes is
  // different from the size of the vector before applying the maximum length
  // constraint then the polyline was obviously updated
  if( n_vertex != vector_bnd_vertices.size() )
   {
    max_length_applied = true;
   }

  // Copy across
  vector_bnd_vertices.resize(n_vertex);
  for(unsigned i=0;i<n_vertex;i++)
   {
    vector_bnd_vertices[i].resize(3);
    vector_bnd_vertices[i][0]=extended_vector[i][0];
    vector_bnd_vertices[i][1]=extended_vector[i][1];
    vector_bnd_vertices[i][2]=extended_vector[i][2];
   }

  // Delete the allocated memory for the geometric object
  // that represents the curvilinear boundary
  delete mesh_geom_obj_pt;

  return max_length_applied;

 }

//=========================================================================
/// \short Gets the associated vertex number according to the vertex
/// coordinates on the destination boundary
//=========================================================================
template<class ELEMENT>
bool RefineableTriangleMesh<ELEMENT>::
get_connected_vertex_number_on_dst_boundary(
  Vector<double> &vertex_coordinates,
  const unsigned &dst_bnd_id,
  unsigned &vertex_number)
 {

  bool found_associated_vertex_number = false;

  // Get the pointer to the associated polyline by using the boundary id
  TriangleMeshPolyLine *dst_polyline =
    this->boundary_polyline_pt(dst_bnd_id);

  unsigned n_vertices = dst_polyline->nvertex();

  // Loop over the vertices and return the closest vertex
  // to the given vertex coordinates
  for (unsigned i = 0; i < n_vertices; i++)
   {
    Vector<double> current_vertex =
    dst_polyline->vertex_coordinate(i);

    double error =
    (vertex_coordinates[0] - current_vertex[0])*
    (vertex_coordinates[0] - current_vertex[0])
    +
    (vertex_coordinates[1] - current_vertex[1])*
    (vertex_coordinates[1] - current_vertex[1]);

    error = sqrt(error);

    if(error <
      ToleranceForVertexMismatchInPolygons::Tolerable_error)
     {
      vertex_number = i;
      found_associated_vertex_number = true;
      break;
     }

   }

  return found_associated_vertex_number;

 }

//=========================================================================
/// \short Restore the connections on the specific internal boundary since
/// there could be changes on the vertices numbering when adding or
/// erasing nodes (vertices)
//=========================================================================
template<class ELEMENT>
bool RefineableTriangleMesh<ELEMENT>::
restore_connections_on_internal_boundary(
 TriangleMeshPolyLine* polyline_pt, const bool first_try)
{
 
#ifdef PARANOID
 // Get the associated boundary id of the current polyline
 unsigned bnd_id = polyline_pt->boundary_id();
#endif
 
 // Verify if this polyline is connected to another polyline
 
 // ****************************************************************
 // 1) If it is connected then get the proper vertex number on the
 //    destination polyline
 // 2) Assign the proper vertex number
 // ****************************************************************
 
 // ****************************************************************
 // First check the initial end
 if (polyline_pt->is_initial_vertex_connected())
  {
   // We need to get the boundary id of the destination/connected
   // boundary
   unsigned dst_bnd_id_initial = 
    polyline_pt->initial_vertex_connected_bnd_id();
   
   // Get the vertex number according to the vertex coordinates
   // on the source boundary
   Vector<double> src_vertex_coordinates_initial =
    polyline_pt->vertex_coordinate(0);
   
   unsigned n_vertex_connection_initial;
   
   bool found_vertex_on_dst_boundary_initial =
    get_connected_vertex_number_on_dst_boundary(
     src_vertex_coordinates_initial, 
     dst_bnd_id_initial, 
     n_vertex_connection_initial);
   
   // If no found it is because the original polyline could be
   // reversed, try again after reversing!!!
   if (!found_vertex_on_dst_boundary_initial)
    {
     if (first_try)
      {
       return false;
      }
     else
      {
       // If no found again then there is a problem with the vertices
#ifdef PARANOID
       std::ostringstream error_message;
       error_message 
        << "It was not possible to find the associated "
        << "vertex number on the\ndestination boundary ("
        << dst_bnd_id_initial
        << ").\nThe source boundary is (" << bnd_id << ") and "
        << "the vertex trying to find on\nthe destination boundary "
        << "is (" << src_vertex_coordinates_initial[0] << ","
        << src_vertex_coordinates_initial[1]<< ")\n"
        << "Initial vertex connection\n";
       throw OomphLibError(
        error_message.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
#endif
      } // else if (first_try)
     
    } // if (!found_vertex_on_dst_boundary_initial)
   
   polyline_pt->initial_vertex_connected_n_vertex() = 
    n_vertex_connection_initial;
   
  }
 
 // ****************************************************************
 // and now the final end
 if (polyline_pt->is_final_vertex_connected())
  {
   // We need to get the boundary id of the destination/connected
   // boundary
   unsigned dst_bnd_id_final = polyline_pt->final_vertex_connected_bnd_id();
   
   // Get the vertex number according to the vertex coordinates
   // on the source boundary
   unsigned tmp_n_vertices = polyline_pt->nvertex();
   Vector<double> src_vertex_coordinates_final =
    polyline_pt->vertex_coordinate(tmp_n_vertices-1);
   
   unsigned n_vertex_connection_final;
   
   bool found_vertex_on_dst_boundary_final =
    get_connected_vertex_number_on_dst_boundary(
     src_vertex_coordinates_final, 
     dst_bnd_id_final, 
     n_vertex_connection_final);
   
   // If no found it is because the original polyline could be
   // reversed, try again after reversing!!!
   if (!found_vertex_on_dst_boundary_final)
    {
     if (first_try)
      {
       return false;
      }
     else
      {
       // If no found again then there is a problem with the vertices
#ifdef PARANOID
       std::ostringstream error_message;
       error_message 
        << "It was not possible to find the associated "
        << "vertex number on the\ndestination boundary ("
        << dst_bnd_id_final 
        << ").\nThe source boundary is (" << bnd_id << ") and "
        << "the vertex trying to find on\nthe destination boundary "
        << "is (" << src_vertex_coordinates_final[0] << ","
        << src_vertex_coordinates_final[1]<< ")\n"
        << "Initial vertex connection\n";
       throw OomphLibError(
        error_message.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
#endif
      }
     
    } // if (!found_vertex_on_dst_boundary_final)
   
   polyline_pt->final_vertex_connected_n_vertex() = 
    n_vertex_connection_final;
   
  }
 
 return true;
 
}


 
//=========================================================================
/// \short Helper function
/// Creates an unsorted face mesh representation from the specified
/// boundary id. It means that the elements are not sorted along the
/// boundary
//=========================================================================
template<class ELEMENT>
void RefineableTriangleMesh<ELEMENT>::
create_unsorted_face_mesh_representation(
  const unsigned &boundary_id,
  Mesh* face_mesh_pt)
 {

  // Create a face mesh adjacent to specified boundary.
  // The face mesh consists of FaceElements that may also be
  // interpreted as GeomObjects

  // Build the face mesh
  this->template build_face_mesh
  <ELEMENT,FaceElementAsGeomObject>(boundary_id, face_mesh_pt);

  // Find the total number of added elements
  unsigned n_element = face_mesh_pt->nelement();
  // Loop over the elements
  for(unsigned e=0;e<n_element;e++)
   {

    //Cast the element pointer to the correct thing!
    FaceElementAsGeomObject<ELEMENT>* el_pt=
    dynamic_cast<FaceElementAsGeomObject<ELEMENT>*>
    (face_mesh_pt->element_pt(e));

    // Set bulk boundary number
    el_pt->set_boundary_number_in_bulk_mesh(boundary_id);

   }

 }

//=========================================================================
/// \short Helper function
/// Creates a sorted face mesh representation of the specified PolyLine
/// It means that the elements are sorted along the boundary
//=========================================================================
template<class ELEMENT>
void RefineableTriangleMesh<ELEMENT>::
create_sorted_face_mesh_representation(
  const unsigned &boundary_id,
  Mesh* face_mesh_pt,
  std::map<FiniteElement*, bool> &is_inverted,
  bool &inverted_face_mesh)
 {

  Mesh *tmp_unsorted_face_mesh_pt = new Mesh();

  // First step we get the unsorted version of the face mesh
  create_unsorted_face_mesh_representation(
    boundary_id, tmp_unsorted_face_mesh_pt);

  // Once with the unsorted version of the face mesh
  // only left to sort it out!!!

  // Put all face elements in order
  //-------------------------------

  // Put first element into ordered list
  // Temporal list for sorting the elements
  std::list<FiniteElement*> ordered_el_pt;
  FiniteElement* el_pt = tmp_unsorted_face_mesh_pt->finite_element_pt(0);
  ordered_el_pt.push_back(el_pt);

  // Number of nodes
  unsigned nnod=el_pt->nnode();

  // Count elements that have been done
  unsigned count_done=0;

  // How many face elements are there?
  unsigned n_face_element = tmp_unsorted_face_mesh_pt->nelement();

  // Keep track of who's done
  std::map<FiniteElement*,bool> done_el;

  is_inverted.clear();

  // Fit in the other elements in at most nel^2 loops
  for (unsigned ee=1;ee<n_face_element;ee++)
   {
    // Loop over all elements to check if they fit to the right
    // or the left of the current one
    for (unsigned e=1;e<n_face_element;e++)
     {
      // Candidate element
      el_pt=tmp_unsorted_face_mesh_pt->finite_element_pt(e);

      // Is it done yet?
      if (!done_el[el_pt])
       {
        // Left and rightmost elements
        FiniteElement* first_el_pt=(*ordered_el_pt.begin());
        std::list<FiniteElement*>::iterator it=ordered_el_pt.end();
        it--;
        FiniteElement* last_el_pt=*it;

        // Left and rightmost nodes
        Node* left_node_pt=first_el_pt->node_pt(0);
        if (is_inverted[first_el_pt])
         {
          left_node_pt=first_el_pt->node_pt(nnod-1);
         }
        Node* right_node_pt=last_el_pt->node_pt(nnod-1);
        if (is_inverted[last_el_pt])
         {
          right_node_pt=last_el_pt->node_pt(0);
         }

        // New element fits at the left of first element and is not inverted
        if (left_node_pt==el_pt->node_pt(nnod-1))
         {
          ordered_el_pt.push_front(el_pt);
          done_el[el_pt]=true;
          count_done++;
          is_inverted[el_pt]=false;
         }
        // New element fits at the left of first element and is inverted

        else if (left_node_pt==el_pt->node_pt(0))
         {
          ordered_el_pt.push_front(el_pt);
          done_el[el_pt]=true;
          count_done++;
          is_inverted[el_pt]=true;
         }
        // New element fits on the right of last element and is not inverted

        else if(right_node_pt==el_pt->node_pt(0))
         {
          ordered_el_pt.push_back(el_pt);
          done_el[el_pt]=true;
          count_done++;
          is_inverted[el_pt]=false;
         }
        // New element fits on the right of last element and is inverted

        else if (right_node_pt==el_pt->node_pt(nnod-1))
         {
          ordered_el_pt.push_back(el_pt);
          done_el[el_pt]=true;
          count_done++;
          is_inverted[el_pt]=true;
         }

        if (done_el[el_pt])
         {
          break;
         }
       }
     }
   }

  // Are we done?
  if (count_done!=(n_face_element-1))
   {
    std::ostringstream error_message;
    error_message
    << "When ordering FaceElements on  "
    << "boundary " << boundary_id << " only managed to order \n" << count_done
    << " of " << n_face_element << " face elements.\n"
    << std::endl;
    throw OomphLibError(
      error_message.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
   }

  // Now make a mesh that contains the FaceElements in order
  // Remember that we currently have a list, not a mesh of sorted elements

  // Fill it
  for (std::list<FiniteElement*>::iterator it=ordered_el_pt.begin();
    it!=ordered_el_pt.end();it++)
   {
    // Get element
    FiniteElement* el_pt=*it;

    // add this face element to the order original mesh
    face_mesh_pt->add_element_pt(el_pt);
   }

  // Verify if face mesh representation is not inverted according to the
  // polyline specified by the user, it means that the initial and the
  // final vertex does really correspond to the first and last vertex
  // respectively, if not, state that the face mesh representation is
  // inverted

  // Get the associated polyline representation to the boundary
  TriangleMeshPolyLine *bnd_polyline =
  this->Boundary_curve_section_pt[boundary_id];

  // Get the really first vertex
  Vector<double> first_vertex =
  bnd_polyline->vertex_coordinate(0);

  // Now get the first node based on the face mesh representation
  // First get access to the first element
  FiniteElement* first_el_pt =
  face_mesh_pt->finite_element_pt(0);

  // Now get access to the first node
  unsigned n_node = first_el_pt->nnode();
  // Get the very first node (taking into account if it is
  // inverted or not!!)
  Node* first_node_pt = first_el_pt->node_pt(0);
  if (is_inverted[first_el_pt])
   {
    first_node_pt = first_el_pt->node_pt(n_node-1);
   }

  double error = (first_node_pt->x(0) - first_vertex[0])*
  (first_node_pt->x(0) - first_vertex[0]) +
  (first_node_pt->x(1) - first_vertex[1])*
  (first_node_pt->x(1) - first_vertex[1]);

  error = sqrt(error);

  if(error <
    ToleranceForVertexMismatchInPolygons::Tolerable_error)
   {
    inverted_face_mesh = false;
   }
  else
   {
    inverted_face_mesh = true;
   }

 }

//=========================================================================
/// Helper function to construct face mesh representation of all polylines,
/// possibly with segments re-distributed between polylines
/// to maintain an approximately even sub-division of the polygon
//=========================================================================
template<class ELEMENT>
void RefineableTriangleMesh<ELEMENT>::
get_face_mesh_representation(TriangleMeshPolygon* polygon_pt,
  Vector<Mesh*>& face_mesh_pt)
 {

  // Number of polylines
  unsigned n_polyline = polygon_pt->npolyline();
  face_mesh_pt.resize(n_polyline);

  // Are we eligible for re-distributing polyline segments between
  // polylines? We're not if any of the boundaries are associated
  // with a GeomObject because we're then tied to the start and
  // end coordinates along it.
  bool eligible_for_segment_redistribution=true;

  // Loop over constituent polylines
  for(unsigned p=0;p<n_polyline;p++)
   {

    //Get the boundary id of the polyline
    unsigned bound =
    polygon_pt->polyline_pt(p)->boundary_id();

    //If the boundary has a geometric object representation then
    //we can't redistribute
    GeomObject* const geom_object_pt =
    this->boundary_geom_object_pt(bound);
    if(geom_object_pt!=0)
     {
      eligible_for_segment_redistribution=false;
     }

    face_mesh_pt[p] = new Mesh();
    create_unsorted_face_mesh_representation(
      bound, face_mesh_pt[p]);

   }

  if (!polygon_pt->is_redistribution_of_segments_between_polylines_enabled())
   {
    return;
   }

  //If there is more than one region we have to think... Die for now.
  if(this->nregion() > 1)
   {
    std::ostringstream warn_message;
    warn_message
    << "Can't currently re-distribute segments between polylines if there\n"
    << "are multiple regions; returning..." << std::endl;
    OomphLibWarning(warn_message.str(),
      "RefineableTriangleMesh::get_face_mesh_representation()",
      OOMPH_EXCEPTION_LOCATION);
    return;
   }

  // Redistribution overruled
  if (!eligible_for_segment_redistribution)
   {
    std::ostringstream warn_message;
    warn_message
    << "Over-ruling re-distribution of segments between polylines\n"
    << "because at least one boundary is associated with a GeomObject."
    << "Returning..." << std::endl;
    OomphLibWarning(warn_message.str(),
      "RefineableTriangleMesh::get_face_mesh_representation()",
      OOMPH_EXCEPTION_LOCATION);
    return;
   }

  // Create a vector for ordered face mesh
  Vector<Mesh*> ordered_face_mesh_pt(n_polyline);

  // Storage for the total arclength of polygon
  double s_total=0.0;

  // Storage for first and last nodes on polylines so we can figure
  // out if they are inverted relative to each other
  Vector<Node*> first_polyline_node_pt(n_polyline);
  Vector<Node*> last_polyline_node_pt(n_polyline);
  std::vector<bool> is_reversed(n_polyline,false);

  // Loop over constituent polylines 
  for(unsigned p=0;p<n_polyline;p++)
   {

    // Put all face elements in order
    //-------------------------------

    // Put first element into ordered list
    std::list<FiniteElement*> ordered_el_pt;
    FiniteElement* el_pt=face_mesh_pt[p]->finite_element_pt(0);
    ordered_el_pt.push_back(el_pt);

    // Number of nodes
    unsigned nnod=el_pt->nnode();

    // Default for first and last node on polyline
    first_polyline_node_pt[p]=el_pt->node_pt(0);
    last_polyline_node_pt[p]=el_pt->node_pt(nnod-1);

    // Count elements that have been done
    unsigned count_done=0;

    // How many face elements are there?
    unsigned n_face_element = face_mesh_pt[p]->nelement();

    //Get the boundary id of the polyline
    unsigned bound =
    polygon_pt->polyline_pt(p)->boundary_id();

    // Keep track of who's done
    std::map<FiniteElement*,bool> done_el;

    // Keep track of which element is inverted
    std::map<FiniteElement*,bool> is_inverted;

    // Fit in the other elements in at most nel^2 loops
    for (unsigned ee=1;ee<n_face_element;ee++)
     {
      // Loop over all elements to check if they fit to the right
      // or the left of the current one
      for (unsigned e=1;e<n_face_element;e++)
       {
        // Candidate element
        el_pt=face_mesh_pt[p]->finite_element_pt(e);

        // Is it done yet?
        if (!done_el[el_pt])
         {
          // Left and rightmost elements 
          FiniteElement* first_el_pt=(*ordered_el_pt.begin());
          std::list<FiniteElement*>::iterator it=ordered_el_pt.end();
          it--;
          FiniteElement* last_el_pt=*it;

          // Left and rightmost nodes
          Node* left_node_pt=first_el_pt->node_pt(0);
          if (is_inverted[first_el_pt])
           {
            left_node_pt=first_el_pt->node_pt(nnod-1);
           }
          Node* right_node_pt=last_el_pt->node_pt(nnod-1);
          if (is_inverted[last_el_pt])
           {
            right_node_pt=last_el_pt->node_pt(0);
           }

          // New element fits at the left of first element and is not inverted
          if (left_node_pt==el_pt->node_pt(nnod-1))
           {
            ordered_el_pt.push_front(el_pt);
            done_el[el_pt]=true;
            count_done++;
            is_inverted[el_pt]=false;
            first_polyline_node_pt[p]=el_pt->node_pt(0);
           }
          // New element fits at the left of first element and is inverted

          else if (left_node_pt==el_pt->node_pt(0))
           {
            ordered_el_pt.push_front(el_pt);
            done_el[el_pt]=true;
            count_done++;
            is_inverted[el_pt]=true;
            first_polyline_node_pt[p]=el_pt->node_pt(nnod-1);
           }
          // New element fits on the right of last element and is not inverted

          else if(right_node_pt==el_pt->node_pt(0))
           {
            ordered_el_pt.push_back(el_pt);
            done_el[el_pt]=true;
            count_done++;
            is_inverted[el_pt]=false;
            last_polyline_node_pt[p]=el_pt->node_pt(nnod-1);
           }
          // New element fits on the right of last element and is inverted

          else if (right_node_pt==el_pt->node_pt(nnod-1))
           {
            ordered_el_pt.push_back(el_pt);
            done_el[el_pt]=true;
            count_done++;
            is_inverted[el_pt]=true;
            last_polyline_node_pt[p]=el_pt->node_pt(0);
           }

          if (done_el[el_pt])
           {
            break;
           }
         }
       }
     }

    // Are we done?
    if (count_done!=(n_face_element-1))
     {
      std::ostringstream error_message;
      error_message
      << "When ordering FaceElements on  "
      << "boundary " << bound << " only managed to order \n" << count_done
      << " of " << n_face_element << " face elements.\n"
      << std::endl;
      throw OomphLibError(
        error_message.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
     }

    // Now make a mesh that contains the FaceElements in order
    ordered_face_mesh_pt[p] = new Mesh;

    // Fill it
    for (std::list<FiniteElement*>::iterator it=ordered_el_pt.begin();
      it!=ordered_el_pt.end();it++)
     {
      // Get element
      FiniteElement* el_pt=*it;

      // add this face element to the order original mesh
      ordered_face_mesh_pt[p]->add_element_pt(el_pt);
     }

    //Get the arclength along the polygon
    for(unsigned e=0;e<n_face_element;++e)
     {
      FiniteElement* el_pt=ordered_face_mesh_pt[p]->finite_element_pt(e);
      unsigned n_node=el_pt->nnode();
      double element_length_squared=0.0;
      for(unsigned i=0;i<2;i++)
       {
        element_length_squared += pow(el_pt->node_pt(n_node-1)->x(i)-
          el_pt->node_pt(0)->x(i),2);
       }

      // Determine element length
      double element_length=sqrt(element_length_squared);

      // Add this length to the total arclength
      s_total += element_length;
     }

    // Empty the original meshes
    face_mesh_pt[p]->flush_element_and_node_storage();
   }

  // Is first one reversed?
  if ((last_polyline_node_pt[0]==first_polyline_node_pt[1])||
    (last_polyline_node_pt[0]==last_polyline_node_pt[1]))
   {
    is_reversed[0]=false;
   }
  else if ((first_polyline_node_pt[0]==first_polyline_node_pt[1])||
    (first_polyline_node_pt[0]==last_polyline_node_pt[1]))
   {
    is_reversed[0]=true;
   }

  // Reorder the face meshes so that they are contiguous
  Vector<Mesh*> tmp_face_mesh_pt(n_polyline);
  std::vector<bool> mesh_done(n_polyline,false);
  Vector<unsigned> old_polyline_number(n_polyline);

  // Initial entry
  tmp_face_mesh_pt[0]=ordered_face_mesh_pt[0];
  unsigned current=0;
  old_polyline_number[0]=0;
  unsigned count_found=0;

  // Fill in the next entries
  for(unsigned p=1;p<n_polyline;p++)
   {
    Node* end_node_pt=last_polyline_node_pt[current];
    if (is_reversed[current])
     {
      end_node_pt=first_polyline_node_pt[current];
     }

    // Loop over all remaining face meshes to see which one fits
    for(unsigned pp=1;pp<n_polyline;pp++)
     {
      if (!mesh_done[pp])
       {
        // Current one is not reversed, candidate is not reversed
        if ((!is_reversed[current])&&
          (end_node_pt==first_polyline_node_pt[pp]))
         {
          tmp_face_mesh_pt[p]=ordered_face_mesh_pt[pp];
          mesh_done[pp]=true;
          is_reversed[pp]=false;
          old_polyline_number[p]=pp;
          current=pp;
          count_found++;
          break;
         }
        // Current one is not reversed, candidate is reversed

        else if ((!is_reversed[current])&&
          (end_node_pt==last_polyline_node_pt[pp]))
         {
          tmp_face_mesh_pt[p]=ordered_face_mesh_pt[pp];
          mesh_done[pp]=true;
          is_reversed[pp]=true;
          old_polyline_number[p]=pp;
          current=pp;
          count_found++;
          break;
         }
        // Current one is reversed, candidate is not reversed

        else if ((is_reversed[current])&&
          (end_node_pt==first_polyline_node_pt[pp]))
         {
          tmp_face_mesh_pt[p]=ordered_face_mesh_pt[pp];
          mesh_done[pp]=true;
          is_reversed[pp]=false;
          old_polyline_number[p]=pp;
          current=pp;
          count_found++;
          break;
         }
        // Current one is reversed, candidate is reversed

        else if ((is_reversed[current])&&
          (end_node_pt==last_polyline_node_pt[pp]))
         {
          tmp_face_mesh_pt[p]=ordered_face_mesh_pt[pp];
          mesh_done[pp]=true;
          is_reversed[pp]=true;
          old_polyline_number[p]=pp;
          current=pp;
          count_found++;
          break;
         }
       }
     }
   }

#ifdef PARANOID
  if (count_found!=n_polyline-1)
   {
    std::ostringstream error_message;
    error_message << "Only found " << count_found
    << " out of " << n_polyline-1
    << " polylines to be fitted in.\n";
    throw OomphLibError(
      error_message.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Now overwrite the re-ordered data
  for (unsigned i=0;i<n_polyline;i++)
   {
    ordered_face_mesh_pt[i]=tmp_face_mesh_pt[i];
   }

  // Now do an approximate equidistribution of polylines
  //----------------------------------------------------
  double s=0.0;
  unsigned new_face_id=0;

  // Matrix map to indicate if node must not be removed from specified
  // boundary (!=0) or not (=0). Initialises itself to zero
  std::map<Node*,std::map<unsigned,unsigned> >
  node_must_not_be_removed_from_boundary_flag;

  // Loop over the old face mesh
  for(unsigned p=0;p<n_polyline;p++)
   {
    // Loop over the face elements
    unsigned n_face_element = ordered_face_mesh_pt[p]->nelement();
    for (unsigned e=0;e<n_face_element;e++)
     {
      unsigned el_number=e;
      if (is_reversed[p])
       {
        el_number=n_face_element-e-1;
       }

      FiniteElement* el_pt=
      ordered_face_mesh_pt[p]->finite_element_pt(el_number);
      unsigned n_node = el_pt->nnode();

      // Determine element length
      double element_length_squared=0.0;
      for(unsigned i=0;i<2;i++)
       {
        element_length_squared += pow(el_pt->node_pt(n_node-1)->x(i)-
          el_pt->node_pt(0)->x(i),2);
       }
      double element_length=sqrt(element_length_squared);

      // Add this length to the total arclength
      s += element_length;

      // Check if the current 'arclength' is less than the
      // whole 'arclength' divided by the number of polylines
      if(s < s_total/double(n_polyline)+1e-6)
       {
        // If so add this face element to the new face mesh
        face_mesh_pt[new_face_id]->add_element_pt(el_pt);

        unsigned bound_old =
        polygon_pt->polyline_pt(old_polyline_number[p])->boundary_id();

        unsigned bound_new =
        polygon_pt->polyline_pt(new_face_id)->boundary_id();

        // Loop over the nodes in the element
        for(unsigned i=0;i<n_node;i++)
         {
          // Get the pointer to the node
          Node* nod_pt=el_pt->node_pt(i);

          // If the two boundary id's are different, the face element's nodes
          // have to be added to the new boundary
          if(bound_new != bound_old)
           {
            // Add it to the new boundary
            add_boundary_node(bound_new,nod_pt);

            // We are happy for this node to be removed from the
            // old boundary? 
            node_must_not_be_removed_from_boundary_flag[nod_pt][bound_old]+=0;
           }

          // If the face element hasn't moved, its nodes MUST remain
          // on that boundary (incl. any nodes that ar shared by
          // FaceElements that have moved (see above)

          else
           {
            node_must_not_be_removed_from_boundary_flag[nod_pt][bound_old]+=1;
           }
         }
       }

      // If not, reset the current 'arclength' to zero,
      // increase the new face id by one and go one element
      // back by decreasing e by one to make sure the current
      // element gets added to the next face mesh      

      else
       {
        if(new_face_id!=n_polyline-1)
         {
          s=0.0;
          new_face_id++;
          --e;
         }
        else
         {
          s=0.0;
          --e;
         }
       }
     }
   } // end of loop over all polylines -- they are now re-distributed


  // Loop over all nodes on the boundaries of the polygon to remove
  // nodes from boundaries they are no longer on
  unsigned move_count=0;
  for (std::map<Node*,std::map<unsigned,unsigned> >::iterator
    it=node_must_not_be_removed_from_boundary_flag.begin();
    it!=node_must_not_be_removed_from_boundary_flag.end();it++)
   {
    // Get the node
    Node* nod_pt=(*it).first;

    // Now we loop over the boundaries that this node is on
    for (std::map<unsigned,unsigned>::iterator
      it_2=(*it).second.begin();it_2!=(*it).second.end();it_2++)
     {
      // Get the boundary id
      unsigned bound=(*it_2).first;

      // Remove it from that boundary?
      if((*it_2).second==0)
       {
        remove_boundary_node(bound,nod_pt);
        move_count++;
       }
     }
   }

  // Loop over the new face mesh to assign new boundary IDs
  for(unsigned p=0;p<n_polyline;p++)
   {
    //Get the boundary id of the polyline
    unsigned bound =
    polygon_pt->polyline_pt(p)->boundary_id();

    // Loop over the face elements
    unsigned n_face_element = face_mesh_pt[p]->nelement();
    for(unsigned e=0;e<n_face_element;e++)
     {
      //Cast the element pointer to the correct thing!
      FaceElementAsGeomObject<ELEMENT>* el_pt=
      dynamic_cast<FaceElementAsGeomObject<ELEMENT>*>
      (face_mesh_pt[p]->element_pt(e));

      // Set bulk boundary number
      el_pt->set_boundary_number_in_bulk_mesh(bound);
     }
   }

  // Update look-up for elements next to boundary
  setup_boundary_element_info();

  // Now re-create the boundary coordinates
  for(unsigned p=0;p<n_polyline;p++)
   {
    //Get the boundary id of the polyline
    unsigned bound =
    polygon_pt->polyline_pt(p)->boundary_id();

    // Do it
    this->setup_boundary_coordinates(bound);
   }

  // Clean up
  for(unsigned p=0;p<n_polyline;p++)
   {
    // Flush the nodes from the face mesh to make sure we
    // don't delete them (the face mesh that we're returning from here
    // still needs them!)
    ordered_face_mesh_pt[p]->flush_element_and_node_storage();
    delete ordered_face_mesh_pt[p];
   }

 }

//=========================================================================
/// Helper function to construct face mesh representation of all polylines
//=========================================================================
template<class ELEMENT>
void RefineableTriangleMesh<ELEMENT>::
get_face_mesh_representation(
  TriangleMeshOpenCurve* open_polyline_pt,
  Vector<Mesh*>& face_mesh_pt)
{

 // Number of polylines
 unsigned n_polyline = open_polyline_pt->ncurve_section();
 face_mesh_pt.resize(n_polyline);

 // Loop over constituent polylines
 for(unsigned p=0;p<n_polyline;p++)
  {

   //Get the boundary id of the polyline
   unsigned bound =
    open_polyline_pt->curve_section_pt(p)->boundary_id();

   face_mesh_pt[p] = new Mesh();
   create_unsorted_face_mesh_representation(
    bound, face_mesh_pt[p]);

  }

}

//======================================================================
/// Update the PSLG that define the inner boundaries of the mesh.
///Optional boolean is used to run it as test only (if 
/// true is specified as input) in which case PSLG isn't actually
/// modified. Returned boolean indicates if PSLG was (or would have
/// been -- if called with check_only=false) changed.
//======================================================================
template <class ELEMENT>
bool RefineableTriangleMesh<ELEMENT>::
surface_remesh_for_inner_hole_boundaries(Vector<Vector<double> >
  &internal_point_coord,
  const bool& check_only)
 {
  //Boolean to indicate whether an actual update of the internal
  // holes was performed
  bool update_was_performed=false;
  //Loop over the number of internal boundaries
  unsigned n_hole = internal_point_coord.size();
  for(unsigned ihole=0;ihole<n_hole;ihole++)
   {
    //Cache the pointer to the polygon representation
    TriangleMeshPolygon* const poly_pt
    = this->Internal_polygon_pt[ihole];

    //Can the polygon update its own configuration, in which case this
    //is easy
    if(poly_pt->can_update_reference_configuration())
     {
      poly_pt->reset_reference_configuration();

      // Initialize Vector hole_coordinates
      internal_point_coord[ihole].resize(2);

      // Get the vector of hole coordinates
      internal_point_coord[ihole]=poly_pt->internal_point();
     }
    //Otherwise we have to work much harder

    else
     {
      //if we only want to check whether an update of the inner
      //hole is necessary
      if(check_only)
       {
        //is it necessary?
        bool update_necessary=
        this->update_polygon_using_face_mesh(poly_pt,check_only);

        //Yes?
        if(update_necessary)
         {
          //then we have to adapt and return 'true'
          return true;
         }
       }
      //if we not only want to check, then we actually perform
      //the update

      else
       {
        update_was_performed=
        this->update_polygon_using_face_mesh(poly_pt);
       }

      if (!poly_pt->internal_point().empty())
       {

        //Now sort out the hole coordinates
        Vector<double> vertex_coord;
        unsigned n_polyline = poly_pt->npolyline();

        // Initialize Vector hole_coordinates
        vertex_coord.resize(2);
        internal_point_coord[ihole].resize(2);

        //Hole centre will be found by averaging the position of
        //all vertex nodes
        internal_point_coord[ihole][0] = 0.0;
        internal_point_coord[ihole][1] = 0.0;

        for(unsigned p=0;p<n_polyline;p++)
         {
          Vector<double> poly_ave(2,0.0);
          //How many vertices are there in the segment
          unsigned n_vertex = poly_pt->polyline_pt(p)->nvertex();
          for(unsigned v=0;v<n_vertex;v++)
           {
            vertex_coord = poly_pt->polyline_pt(p)->vertex_coordinate(v);
            for(unsigned i=0;i<2;i++)
             {
              poly_ave[i] += vertex_coord[i];
             }
           }

          //Add the average polyline coordinate to the hole centre
          for(unsigned i=0;i<2;i++)
           {
            internal_point_coord[ihole][i] += poly_ave[i]/n_vertex;
           }
         }

        //Now average out the hole centre
        for(unsigned i=0;i<2;i++)
         {
          internal_point_coord[ihole][i] /= n_polyline;
         }

        //Set the new hole centre
        poly_pt->internal_point() = internal_point_coord[ihole];
        
      }

     }

   } //End of the action (n_hole for)

  if(check_only)
   {
    // If we make it up to here and we only check then no update is required
    return false;
   }
  else
   {
    // otherwise indicate whether an actual update was performed
    return update_was_performed;
   }

 } //End of the loop of internal boundaries

 //======================================================================
 /// Create the polylines and fill associate data structures, used when
 /// creating from a mesh from polyfiles
 //======================================================================
 template<class ELEMENT>
 void RefineableTriangleMesh<ELEMENT>::
 create_polylines_from_polyfiles(const std::string& node_file_name,
                                 const std::string& poly_file_name)
 {
  // Get the nodes coordinates (the index of the nodes to build the
  // polylines is the one used in the node_file_name file)
  // Process node file
  // -----------------
  std::ifstream node_file(node_file_name.c_str(),std::ios_base::in);

  // Check that the file actually opened correctly
  if(!node_file.is_open())
   {
    std::string error_msg("Failed to open node file: ");
    error_msg += "\"" + node_file_name + "\".";
    throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  // Read number of nodes
  unsigned nnodes;
  node_file >> nnodes;
   
  // Spatial dimension of nodes
  unsigned dimension;
  node_file >> dimension;
  
#ifdef PARANOID
  if(dimension!=2)
   {
    throw OomphLibError("The dimension must be 2\n",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
   
  // Storage the nodes vertices
  Vector<double> x_node(nnodes);
  Vector<double> y_node(nnodes);

  // Number of attributes
  unsigned npoint_attributes;
  node_file >> npoint_attributes;;
  
  // Flag for boundary markers
  unsigned boundary_markers_flag;
  node_file >> boundary_markers_flag;

  // Dummy for node number
  unsigned dummy_node_number;
  // Dummy for node attribute
  unsigned dummy_node_attribute;
  // Dummy for node boundary
  unsigned dummy_node_boundary;
   
  // Load in nodal posititions, point attributes
  // and boundary markers
  for(unsigned i=0;i<nnodes;i++)
   {
    node_file>>dummy_node_number;
    node_file>>x_node[i];
    node_file>>y_node[i];
    for(unsigned j=0;j<npoint_attributes;++j)
     {
      node_file>>dummy_node_attribute;
     }
    if(boundary_markers_flag)
     {
      node_file>>dummy_node_boundary;
     }
   }
  node_file.close();
  
  // Get the segments information and use that info. to create the
  // polylines
  
  // A map to store the segments associated to a boundary, non sorted
  std::map<unsigned,Vector<std::pair<unsigned,unsigned> > > 
   unsorted_boundary_segments;
  
  // Independent storage for the boundaries ids found in the segments so that
  // the polylines, and therefore polygons be created in the order they appear
  // in the polyfile
  Vector<unsigned> sorted_boundaries_ids;
   
  // Process poly file to extract edges
  //-----------------------------------
   
  // Open poly file
  std::ifstream poly_file(poly_file_name.c_str(),std::ios_base::in);

  // Check that the file actually opened correctly
  if(!poly_file.is_open())
   {
    std::string error_msg("Failed to open poly file: ");
    error_msg += "\"" + poly_file_name + "\".";
    throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  // Number of nodes in poly file --- these will be ignore
  unsigned n_node_poly;
  poly_file >> n_node_poly;

  // Dimension
  poly_file >> dimension;

  // Attribute flag
  unsigned attribute_flag;
  poly_file >> attribute_flag;

  // Flag for boundary markers
  poly_file >> boundary_markers_flag;
  
  // Ignore node information: Note: No, we can't extract the
  // actual nodes themselves from here!
  unsigned dummy;
  for(unsigned i=0;i<n_node_poly;i++)
   {
    //Read in (and discard) node number and x and y coordinates
    poly_file>>dummy;
    poly_file>>dummy;
    poly_file>>dummy;
    //read in the attributes
    for(unsigned j=0;j<attribute_flag;++j)
     {
      poly_file >> dummy;
     }
    //read in the boundary marker
    if(boundary_markers_flag==1)
     {
      poly_file>>dummy;
     }
   }
  
  // Variable used to read the values from the input file
  unsigned read_value;
  
  // Number of segments
  poly_file >> read_value;
  const unsigned nglobal_segments = read_value;
  
  // Boundary marker flag
  poly_file >> boundary_markers_flag;
  
  // Global segment number
  unsigned global_segment_number;

  // Node identifier set (used to identify possible internal boundaries)
  std::set<unsigned> nodes_ids;

  // Extract information for each segment
  for(unsigned i=0;i<nglobal_segments;i++)
   {
    // Node id on the edge of the segment
    unsigned lnode_id; // left node
    unsigned rnode_id; // right node
    unsigned bnd_id;   // boundary id associated to the current segment
    poly_file >> global_segment_number;
    poly_file >> lnode_id;
    poly_file >> rnode_id;
    nodes_ids.insert(lnode_id);
    nodes_ids.insert(rnode_id);
    if(boundary_markers_flag)
     {
      poly_file >> bnd_id;
     }
    
    // Store the segments info. (use bnd_id - 1 because the nodes and
    // elements associated the bnd_id have been associated by external
    // methods to bnd_id - 1)
    unsorted_boundary_segments[bnd_id-1].push_back(
     std::make_pair(lnode_id, rnode_id));
    
    // Add the boundary id to the vector of boundaries ids only if it
    // has not been added, the polylines will be created using this
    // order
    
    // Get the number of boundaries ids currently sorted
    const unsigned nsorted_boundaries_ids = 
     sorted_boundaries_ids.size();
    // Flag to know if the boundary id was found
    bool boundary_id_found = false;
    for (unsigned ib = 0; ib < nsorted_boundaries_ids; ib++)
     {
      if (sorted_boundaries_ids[ib] == bnd_id - 1)
       {
        boundary_id_found = true;
        break;
       } // if (sorted_boundaries_ids[ib] == bnd_id - 1)
     } // for (ib < nsorted_boundaries_ids)
    
    // If th boundary id has not been added, then add it!!!
    if (!boundary_id_found)
     {
      sorted_boundaries_ids.push_back(bnd_id - 1);
     } // if (!boundary_id_found)
    
   }
  
  // Verify if there are internal boundaries defined, if that is the
  // case we can not continue since we are not yet supporting internal
  // boundaries defined in polyfiles to created a mesh that may be
  // adapted
#ifdef PARANOID
  if (nglobal_segments != nodes_ids.size())
   {
    std::ostringstream error_message;
    error_message
     << "The number of nodes (" << nodes_ids.size() << ") and segments (" 
     << nglobal_segments << ") is different.\nThis may mean that there  "
     << "are internal non-closed boundaries defined in\nthe polyfile. "
     << "If you need this feature please use the TriangleMeshPoyLine\n"
     << "and TriangleMeshCurviLine objects to define your domain.\n\n";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // Now sort the segments associated to a boundary to create a contiguous
  // polyline, but first check that the number of found boundaries be the
  // same as the current number of boundaries in the mesh
  const unsigned nboundary = unsorted_boundary_segments.size();
  
#ifdef PARANOID
  if (nboundary != this->nboundary())
   {
    std::ostringstream error_message;
    error_message
     << "The number of boundaries on the mesh (" << this->nboundary() 
     << ") is different from the number of\nboundaries read from the "
     << "polyfiles (" << unsorted_boundary_segments.size() << ")!!!\n\n\n";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // Get the number of sorted boundaries ids and check that it matches
  // with the total number of boundaries
  const unsigned nsorted_boundaries_ids = 
   sorted_boundaries_ids.size();
#ifdef PARANOID
  if (nsorted_boundaries_ids != this->nboundary())
   {
    std::ostringstream error_message;
    error_message
     << "The number of boundaries on the mesh (" << this->nboundary() 
     << ") is different from the number of\nsorted boundaries ids read "
     << "from the polyfiles (" << nsorted_boundaries_ids << ")!!!\n\n\n";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // Sorted segments (to create a polyline -- boundary)
  std::map<unsigned, std::list<unsigned> > sorted_boundary_segments;
  
  // Go through all the found boundaries
  std::map<unsigned,Vector<std::pair<unsigned,unsigned> > >::iterator it;
  
  for (it = unsorted_boundary_segments.begin(); 
       it != unsorted_boundary_segments.end(); 
       it++)
   {
    // Get the current boundary id, only look for the segments
    // associated with this boundary
    const unsigned bnd_id = (*it).first;
    Vector<std::pair<unsigned, unsigned> > segments_edges = (*it).second;

    // Now sort the segments associated to this boundary
    std::map<std::pair<unsigned, unsigned>, bool> segment_done;
    const unsigned nsegments = segments_edges.size();
    
    // Sorted nodes for the current segment
    std::list<unsigned> sorted_segments;
    
    // Get the left and right node of the zero segment
    unsigned left_node_id = segments_edges[0].first;
    unsigned right_node_id = segments_edges[0].second;
    
    // ...  and add it to the sorted segments structure
    sorted_segments.push_back(left_node_id);
    sorted_segments.push_back(right_node_id);

    // Mark the current segment as done
    segment_done[segments_edges[0]] = true;

    // Set the number of sorted segments
    unsigned nsorted_segments = 1;

    while(nsorted_segments < nsegments)
     {
      for (unsigned i = 1; i < nsegments; i++)
       {
        // Check if the i-th segments has been done
        if (!segment_done[segments_edges[i]])
         {
          // Get the left and right node id
          unsigned current_left_node_id = segments_edges[i].first;
          unsigned current_right_node_id = segments_edges[i].second;
          
          // Now check if the current segment can be added to the left
          // or right side of the sorted segments
          if (current_left_node_id == right_node_id)
           {
            // Add the current_right_node_id to the right of the sorted
            // segments
            sorted_segments.push_back(current_right_node_id);
            // Increase the number of sorted segments
            nsorted_segments++;
            // Mark the segment as done
            segment_done[segments_edges[i]] = true;
            // Update the right most node
            right_node_id = current_right_node_id;
            // Break the for loop
            break;
           }
          else if (current_right_node_id == left_node_id)
           {
            // Add the current_left_node_id to the left of the sorted
            // segments
            sorted_segments.push_front(current_left_node_id);
            // Increase the number of sorted segments
            nsorted_segments++;
            // Mark the segment as done
            segment_done[segments_edges[i]] = true;
            // Update the left most node
            left_node_id = current_left_node_id;
            // Break the for loop
            break;
           }
          else if (current_left_node_id == left_node_id)
           {
            // Add the current_right_node_id to the left of the sorted
            // segments
            sorted_segments.push_front(current_right_node_id);
            // Increase the number of sorted segments
            nsorted_segments++;
            // Mark the segment as done
            segment_done[segments_edges[i]] = true;
            // Update the left most node
            left_node_id = current_right_node_id;
            // Break the for loop
            break;
           }
          else if (current_right_node_id == right_node_id)
           {
            // Add the current_left_node_id to the right of the sorted
            // segments
            sorted_segments.push_back(current_left_node_id);
            // Increase the number of sorted segments
            nsorted_segments++;
            // Mark the segment as done
            segment_done[segments_edges[i]] = true;
            // Update the left most node
            right_node_id = current_left_node_id;
            // Break the for loop
            break;
           }
         } // if (!segment_done[segments_edges[i]])
       } // for (i < nsegments)
     } // while(nsorted_segments < nsegments)
    
    sorted_boundary_segments[bnd_id] = sorted_segments;
    
   } // for (unsorted_boundary_segments.begin(); 
     //      unsorted_boundary_segments.end())

#ifdef PARANOID
  if (sorted_boundary_segments.size() != this->nboundary())
   {
    std::ostringstream error_message;
    error_message
     << "The number of boundaries on the mesh (" << this->nboundary() 
     << ") is different from the number\nof sorted boundaries to create the "
     << "polylines (" << sorted_boundary_segments.size() << ")\n\n";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Now we have the sorted nodes, we can create the polylines by
  // getting the vertices of the nodes
  Vector<TriangleMeshPolyLine*> polylines_pt(nboundary);
  unsigned current_polyline = 0;
  
  // Go through the sorted boundaries using the sorted boundaries ids
  for (unsigned ib = 0; ib < nsorted_boundaries_ids; ib++)
   {
    // Get the boundary id from the vector of sorted boundaries ids
    const unsigned bnd_id = sorted_boundaries_ids[ib];
    
    // Create a vector representation for ease to use
    // Get the vertices of the nodes that create the boundary / polyline
    Vector<unsigned> nodes_ids;
    for (std::list<unsigned>::iterator it_list = 
          sorted_boundary_segments[bnd_id].begin();
         it_list != sorted_boundary_segments[bnd_id].end(); 
         it_list++)
     {nodes_ids.push_back((*it_list));}
    
    // Get the number of vertices for the polyline
    const unsigned nvertices = nodes_ids.size();

    // The storage for the vertices
    Vector<Vector<double> > vertices(nvertices);
    
    // Now get the vertices of the nodes of the current boundary
    for (unsigned i = 0; i < nvertices; i++)
     {
      // Get the vertices
      vertices[i].resize(2);
      vertices[i][0] = x_node[nodes_ids[i]-1];
      vertices[i][1] = y_node[nodes_ids[i]-1];
     }
    
    // Now create the polyline 
     
    // Note: The bnd_id is the real bnd_id (from the input file) - 1
    // since nodes and elements of the current boundary have been
    // associated to bnd_id - 1)
    polylines_pt[current_polyline] = 
     new TriangleMeshPolyLine(vertices, bnd_id);
    
    // Updates bnd_id<--->curve section map
    this->Boundary_curve_section_pt[bnd_id] = 
     dynamic_cast<TriangleMeshCurveSection*>(polylines_pt[current_polyline]);
    
    // Increase the index for the polyline storage
    current_polyline++;
    
   } // for (it_sorted = sorted_boundary_segments.begin();
     //      it_sorted != sorted_boundary_segments.end())
  
  // Now create the polygons or closed curves
  // Sort the polylines to create polygons
  unsigned nsorted_polylines = 0;

  // Number of created polygons
  unsigned npolygons = 0;

  // Storage for the polygons
  Vector<TriangleMeshPolygon*> polygons_pt;
  
  // Mark the already done polylines
  std::map<unsigned, bool> polyline_done;
  while(nsorted_polylines < nboundary)
   {
    // Storage for the curve sections that create a polygon
    std::list<TriangleMeshCurveSection*> sorted_curve_sections_pt;
    
    unsigned init_poly = 0;
#ifdef PARANOID
    bool found_root_polyline = false;
#endif
    // Get the left and right node of the current polyline
    for (unsigned i = 0; i < nboundary; i++)
     {
      if (!polyline_done[i])
       {
        init_poly = i;
        // Increase the number of sorted polylines
        nsorted_polylines++;
#ifdef PARANOID
        // Mark as found the root polyline
        found_root_polyline = true;
#endif
        // Mark the polyline as done
        polyline_done[i] = true;
        // Add the polyline to the curve sections storage
        sorted_curve_sections_pt.push_back(polylines_pt[i]);
        // Break the loop to set we have found a root polyline
        break;
       }
     }
    
#ifdef PARANOID
    if (!found_root_polyline)
     {
      std::ostringstream error_message;
      error_message
       << "Was not possible to found the root polyline to create polygons\n\n";
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    // Get the associated boundary to the current polyline
    const unsigned bnd_id = polylines_pt[init_poly]->boundary_id();
    // Get the initial and final node id of the current polyline
    unsigned left_node_id = sorted_boundary_segments[bnd_id].front();
    unsigned right_node_id = sorted_boundary_segments[bnd_id].back();
    
    // Flag to know that we already have a closed polygon
    bool closed_polygon = false;
    
    do
     {
      // Go through all the polylines
      for (unsigned i = init_poly; i < nboundary; i++)
       {
        // Check that the polyline has not been currently done
        if (!polyline_done[i])
         {
          // Get the initial and final nodes id of the current polyline
          
          // Get the associated boundary to the current polyline
          const unsigned cbnd_id = polylines_pt[i]->boundary_id();
          // Get the initial and final node id of the current polyline
          unsigned cleft_node_id = sorted_boundary_segments[cbnd_id].front();
          unsigned cright_node_id = sorted_boundary_segments[cbnd_id].back();
          
          // Check if the polyline goes to the left or right of the
          // current sorted polylines
          if (cleft_node_id == right_node_id)
           {
            // Add the polyline to the curve section storage
            sorted_curve_sections_pt.push_back(polylines_pt[i]);
            // Mark the polyline as done
            polyline_done[i] = true;
            // Update the right node
            right_node_id = cright_node_id;
            // Increase the number of done polyines
            nsorted_polylines++;
            // Break the for loop
            break;
           }
          else if (cright_node_id == left_node_id)
           {
            // Add the polyline to the curve section storage
            sorted_curve_sections_pt.push_front(polylines_pt[i]);
            // Mark the polyline as done
            polyline_done[i] = true;
            // Update the right node
            left_node_id = cleft_node_id;
            // Increase the number of done polyines
            nsorted_polylines++;
            // Break the for loop
            break;
           }
          else if (cleft_node_id == left_node_id)
           {
            // First reverse the polyline
            polylines_pt[i]->reverse();
            // Add the polyline to the curve section storage
            sorted_curve_sections_pt.push_front(polylines_pt[i]);
            // Mark the polyline as done
            polyline_done[i] = true;
            // Update the right node
            left_node_id = cright_node_id;
            // Increase the number of done polyines
            nsorted_polylines++;
            // Break the for loop
            break;
           }
          else if (cright_node_id == right_node_id)
           {
            // First reverse the polyline
            polylines_pt[i]->reverse();
            // Add the polyline to the curve section storage
            sorted_curve_sections_pt.push_back(polylines_pt[i]);
            // Mark the polyline as done
            polyline_done[i] = true;
            // Update the right node
            right_node_id = cleft_node_id;
            // Increase the number of done polyines
            nsorted_polylines++;
            // Break the for loop
            break;
           }
         } // if (!polyline_done[i])
        
       } // for (i < nboundary)
      
      // We have created a polygon
      if (left_node_id == right_node_id)
       {
        // Set the flag as true
        closed_polygon = true;  
       }
      
     }while(nsorted_polylines < nboundary && !closed_polygon);
    
#ifdef PARANOID
    if (!closed_polygon)
     {
      std::ostringstream error_message;
      error_message
       << "It was not possible to create a closed curve, these are the "
       << "vertices of the already sorted polylines\n\n";
      unsigned cpolyline = 0;
      for (std::list<TriangleMeshCurveSection*>::iterator it_list = 
            sorted_curve_sections_pt.begin(); 
           it_list != sorted_curve_sections_pt.end(); 
           it_list++)
       {
        error_message << "Polyline (" << cpolyline << ")\n";
        TriangleMeshPolyLine *tmp_poly_pt =
         dynamic_cast<TriangleMeshPolyLine*>((*it_list));
        const unsigned nvertex = tmp_poly_pt->nvertex();
        for (unsigned v = 0; v < nvertex; v++)
         {
          error_message <<"("<<tmp_poly_pt->vertex_coordinate(v)[0] 
                        <<", "<<tmp_poly_pt->vertex_coordinate(v)[1]<<")\n";
         }
        error_message << "\n";
        cpolyline++;
       }
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    // Create a vector version to create the polygon from the sorted
    // polyines
    Vector<TriangleMeshCurveSection*> tmp_sorted_curve_sections_pt;
    for (std::list<TriangleMeshCurveSection*>::iterator it_list = 
          sorted_curve_sections_pt.begin(); 
         it_list != sorted_curve_sections_pt.end();
         it_list++)
     {tmp_sorted_curve_sections_pt.push_back((*it_list));}
    
    // Create a new polygon by using the new created polylines
    TriangleMeshPolygon *polygon_pt =
     new TriangleMeshPolygon(tmp_sorted_curve_sections_pt);
    
    // Keep track of new created polygons that need to be deleted!!!
    this->Free_polygon_pt.insert(polygon_pt);
    
    // Store the polygon in the polygons storages
    polygons_pt.push_back(polygon_pt);
    
    npolygons++;
    
   } // while(nsorted_polylines < nboundary)
  
  // ------------------------------------------------------------------
  // Before filling the data structures we need to identify the outer
  // closed boundary and the inner closed boundaries.
  // If the nodes are not in order we throw a warning message
  
  // Index for the polygon that is currently considered as the outer
  // boundary
  unsigned index_outer = 0;
  
  for (unsigned idx_outer = 0; idx_outer < npolygons; idx_outer++)
   {
    // Get the vertices of the outer boundary
    Vector<Vector<double> > outer_vertex_coordinates;
    
    // Flag to know if ALL the inner closed boundaries are inside the
    // outer closed boundary
    bool all_inner_inside = true;
    
    // Number of polylines of the outer boundary
    const unsigned nouter_polylines = polygons_pt[idx_outer]->npolyline();
    for (unsigned p = 0; p < nouter_polylines; p++)
     {
      TriangleMeshPolyLine* tmp_poly_pt = 
       polygons_pt[idx_outer]->polyline_pt(p);
      const unsigned nvertex = tmp_poly_pt->nvertex();
      for (unsigned v = 0; v < nvertex; v++)
       {
        Vector<double> current_vertex = tmp_poly_pt->vertex_coordinate(v);
        outer_vertex_coordinates.push_back(current_vertex);
       } // for (v < nvertex)
     } // for (p < nouter_polylines)
    
    // Now get the vertices for the inner boundaries 
    
    // First get the number of inner closed boundaries (polygons size
    // minus one because one of the polygons is considered to be the
    // outer closed boundary
    const unsigned ninner_polygons = polygons_pt.size() - 1;
    
    // Store the vertices of the inner closed boundaries
    Vector<Vector<Vector<double> > > inner_vertex_coordinates(ninner_polygons);
    // Get all the vertices of the inner closed boundaries
    for (unsigned i = 0; i <= ninner_polygons; i++)
     {    
      if (i != idx_outer)
       {
        // Number of polylines of the current internal closed boundary
        const unsigned ninner_polylines = polygons_pt[i]->npolyline();
        for (unsigned p = 0; p < ninner_polylines; p++)
         {
          TriangleMeshPolyLine* tmp_poly_pt = polygons_pt[i]->polyline_pt(p);
          const unsigned nvertex = tmp_poly_pt->nvertex();
          for (unsigned v = 0; v < nvertex; v++)
           {
            Vector<double> current_vertex = tmp_poly_pt->vertex_coordinate(v);
            if (i < idx_outer)
             {
              inner_vertex_coordinates[i].push_back(current_vertex);
             }
            else if (i > idx_outer)
             {
              inner_vertex_coordinates[i-1].push_back(current_vertex);
             }
           } // for (v < nvertex)
          
         } // for (p < ninner_polylines)
        
       } // if (i != index_outer)
      
     } // for (i <= ninner_polygons)
    
    // Now check that ALL the vertices of ALL the internal closed
    // boundaries are inside the outer closed boundary
    for (unsigned i = 0; i < ninner_polygons; i++)
     {
      // Get the number of vertices in the current internal closed
      // boundary
      const unsigned nvertex_internal = inner_vertex_coordinates[i].size();
      for (unsigned v = 0; v < nvertex_internal; v++)
       {
        // Get a vertex in the current internal closed boundary
        Vector<double> current_point = inner_vertex_coordinates[i][v];
        all_inner_inside &= 
         this->is_point_inside_polygon(outer_vertex_coordinates,
                                       current_point);
        
        // Check if we should continue checking for more points inside
        // the current proposed outer boundary
        if (!all_inner_inside)
         {
          // Break the "for" for the vertices
          break;
         }
        
       } // for (v < nvertex_internal)
      
      // Check if we should continue checking for more inner closed
      // boundaries inside the current proposed outer boundary
      if (!all_inner_inside)
       {
        // Break the "for" for the inner boundaries
        break;
       }
      
     } // for (i < ninner_polygons)
    
    // Check if all the vertices of all the polygones are inside the
    // current proposed outer boundary
    if (all_inner_inside)
     {
      index_outer = idx_outer;
      break;
     }
    
   } // for (idx_outer < npolygons)
  
#ifdef PARANOID
  // Check if the first nodes listed in the polyfiles correspond to
  // the outer boundary, if that is not the case then throw a warning
  // message
  if (index_outer != 0)
   {
    std::ostringstream warning_message;
    warning_message
     << "The first set of nodes listed in the input polyfiles does not\n"
     << "correspond to the outer closed boundary. This may lead to\n"
     << "problems at the adaptation stage if the holes coordinates\n"
     << "are no correctly associated to the inner closed boundaries.\n"
     << "You can check the generated mesh by calling the output() method\n"
     << "from the mesh object '(problem.mesh_pt()->output(string))'\n\n";
    OomphLibWarning(warning_message.str(),
                    OOMPH_CURRENT_FUNCTION,
                    OOMPH_EXCEPTION_LOCATION);
   } // if (index_outer != 0)
#endif
  
  // ------------------------------------------------------------------
  // Now fill the data structures
  
  // Store outer polygon
  this->Outer_boundary_pt = polygons_pt[index_outer];
  
  this->Internal_polygon_pt.resize(npolygons-1);
  for (unsigned i = 0; i < npolygons; i++)
   {
    if (i != index_outer)
     {
      if (i < index_outer)
       {
        // Store internal polygons by copy constructor
        this->Internal_polygon_pt[i] = polygons_pt[i];
       }
      else if (i > index_outer)
       {
        // Store internal polygons by copy constructor
        this->Internal_polygon_pt[i-1] = polygons_pt[i];
       }
     } // if (i != index_outer)
   } // for (i < npolygons)
  
  // Before assigning the hole vertex coordinate to the inner closed
  // boundaries check that the holes are listed in orderm if that is
  // not the case the associate each hole vertex coordinate to the
  // inner closed boundaries
  
  // Store the vertices of the inner closed boundaries
  Vector<Vector<Vector<double> > > inner_vertex_coordinates(npolygons-1);
  // Get all the vertices of the inner closed boundaries
  for (unsigned i = 0; i < npolygons-1; i++)
   {    
    // Number of polylines of the current internal closed boundary
    const unsigned ninner_polylines = 
     this->Internal_polygon_pt[i]->npolyline();
    for (unsigned p = 0; p < ninner_polylines; p++)
     {
      TriangleMeshPolyLine* tmp_poly_pt = 
       this->Internal_polygon_pt[i]->polyline_pt(p);
      // Number of vertices of the current polyline in the current
      // internal closed polygon
      const unsigned nvertex = tmp_poly_pt->nvertex();
      for (unsigned v = 0; v < nvertex; v++)
       {
        Vector<double> current_vertex = tmp_poly_pt->vertex_coordinate(v);
        inner_vertex_coordinates[i].push_back(current_vertex);
       } // for (v < nvertex)
      
     } // for (p < ninner_polylines)
    
   } // for (i <= ninner_polygons)
  
  // Holes information
  unsigned nholes;
  poly_file >> nholes;
  
#ifdef PARANOID
  if (npolygons > 1 && (npolygons - 1) != nholes)
   {
    std::ostringstream error_message;
    error_message
     << "The number of holes (" << nholes << ") does not correspond "
     << "with the number\nof internal polygons (" 
     << npolygons - 1 <<")\n\n"
     << "Using polyfiles as input does not currently allows the\n"
     << "definition of more than one outer polygon\n\n";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // Storage for the holes
  Vector<Vector<double> > hole_coordinates(nholes);
  
  // Dummy for hole number
  unsigned dummy_hole;
  // Loop over the holes to get centre coords
  for(unsigned ihole=0;ihole<nholes;ihole++)
   {
    hole_coordinates[ihole].resize(2);
    // Read the centre value
    poly_file >> dummy_hole;
    poly_file >> hole_coordinates[ihole][0];
    poly_file >> hole_coordinates[ihole][1];
   }
  
  // Vector that store the index of the hole coordinate that
  // correspond to each internal closed polygon
  Vector<unsigned> index_hole_of_internal_polygon(npolygons-1);
  std::map<unsigned, bool> hole_done;
  
  // Now associate each hole vertex to a corresponding internal closed
  // polygon
  for (unsigned i = 0; i < npolygons-1; i++)
   {
    // Find which hole is associated to each internal closed boundary
    for (unsigned h = 0; h < nholes; h++)
     {
      // If the hole has not been previously associated
      if (!hole_done[h])
       {
        // Get the hole coordinate
        Vector<double> current_point = hole_coordinates[h];
        
        const bool hole_in_polygon =
         this->is_point_inside_polygon(inner_vertex_coordinates[i],
                                       current_point);
        
        // If the hole is inside the polygon
        if (hole_in_polygon)
         {
          // Mark the hole as done
          hole_done[h] = true;
          // Associate the current hole with the current inner closed
          // boundary
          index_hole_of_internal_polygon[i] = h;
          // Break the search
          break;
         }
        
       } // if (!hole_done[h])
      
     } // for (h < nholes)
    
   } // for (i < npolygons-1)
  
#ifdef PARANOID
  if (hole_done.size() != npolygons-1)
   {
    std::ostringstream error_message;
    error_message
     << "Not all the holes were associated to an internal closed boundary\n"
     << "Only ("<<hole_done.size()<<") holes were assigned for a total of\n"
     << "(" << npolygons-1 << ") internal closed boundaries.\n"
     << "You can check the generated mesh by calling the output() method\n"
     << "from the mesh object '(problem.mesh_pt()->output(string))'\n\n";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   } // if (index_hole != ihole)
#endif
  
  // Assign the holes coordinates to the internal polygons
  for (unsigned ihole = 0; ihole < nholes; ihole++)
   {
    // Get the index hole of the current internal closed polygon
    const unsigned index_hole = index_hole_of_internal_polygon[ihole];
#ifdef PARANOID
    // Check if the hole index is the same as the internal closed
    // boundary, it means that the holes were listed in the same order
    // as the nodes of the internal closed boundaries
    if (index_hole != ihole)
     {
      std::ostringstream error_message;
      error_message
       << "The hole vertices coordinates are not listed in the same order\n"
       << "as the nodes that define the internal closed boundaries.\n"
       << "This may lead to problems in case that the holes coordinates\n"
       << "were no properly assigned to the internal closed boundaries.\n"
       << "You can check the generated mesh by calling the output() method\n"
       << "from the mesh object '(problem.mesh_pt()->output(string))'\n\n";
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     } // if (index_hole != ihole)
#endif
    
    // Set the hole coordinate for the internal polygon
    this->Internal_polygon_pt[ihole]->internal_point() = 
     hole_coordinates[index_hole];
   }
  
  // Ignore the first line with structure description
  poly_file.ignore(80,'\n');
  
  // Regions information
  unsigned nregions;
  
  // Extract regions information
  // But first check if there are regions or not
  std::string regions_info_string;
  
  // Read line up to termination sign
  getline(poly_file, regions_info_string);
  
  // Check if the read string is a number or a comment wrote by triangle,
  // if it is a number then that is the number of regions 
  if (isdigit(regions_info_string.c_str()[0]))
   {
    nregions = std::atoi(regions_info_string.c_str());
   }
  else
   {
    nregions = 0;
   }
  
  // The regions coordinates
  std::map<unsigned, Vector<double> > regions_coordinates;

  // Dummy for regions number
  unsigned dummy_region;

  unsigned region_id;

  // Loop over the regions to get their coords
  for(unsigned iregion=0;iregion<nregions;iregion++)
   {
    Vector<double> tmp_region_coordinates(2);
    // Read the regions coordinates
    poly_file >> dummy_region;
    poly_file >> tmp_region_coordinates[0];
    poly_file >> tmp_region_coordinates[1];
    poly_file >> region_id;
    regions_coordinates[region_id].resize(2);
    regions_coordinates[region_id][0] = tmp_region_coordinates[0];
    regions_coordinates[region_id][1] = tmp_region_coordinates[1];
    
    // Ignore the first line with structure description
    poly_file.ignore(80,'\n');
    
    // Verify if not using the default region number (zero)
    if (region_id == 0) 
     {
      std::ostringstream error_message;
      error_message << "Please use another region id different from zero.\n"
                    << "It is internally used as the default region number.\n";
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
    
   }

  // Store the extra regions coordinates
  this->Regions_coordinates = regions_coordinates;

  poly_file.close();
  
 }

//======================================================================
/// Helper function that checks if a given point is inside a polygon
//======================================================================
template <class ELEMENT>
bool RefineableTriangleMesh<ELEMENT>::is_point_inside_polygon(
 Vector<Vector<double> > polygon_vertices, Vector<double> point)
{
 // Total number of vertices (the first and last vertex should be the
 // same)
 const unsigned nvertex = polygon_vertices.size();
 
#ifdef PARANOID
 // Check that the first and last vertex of the given polygon are the
 // same
 const double dist_first_last = 
  sqrt(((polygon_vertices[0][0] - polygon_vertices[nvertex-1][0]) *
        (polygon_vertices[0][0] - polygon_vertices[nvertex-1][0])) +
       ((polygon_vertices[0][1] - polygon_vertices[nvertex-1][1]) *
        (polygon_vertices[0][1] - polygon_vertices[nvertex-1][1])));
 
 if (dist_first_last > 
     ToleranceForVertexMismatchInPolygons::Tolerable_error)
  {
   std::ostringstream error_stream;
   error_stream
    << "The start ("
    << polygon_vertices[0][0]<<", "<<polygon_vertices[0][1] 
    << ") and end ("
    << polygon_vertices[nvertex-1][0]<<", "<<polygon_vertices[nvertex-1][1]
    << ") points of the polygon don't match when judged\n"
    << "with the tolerance ("
    << ToleranceForVertexMismatchInPolygons::Tolerable_error
    << ") which is specified in the namespace \nvariable "
    << "ToleranceForVertexMismatchInPolygons::"
    << "Tolerable_error.\n\n"
    << "Feel free to adjust this or to recompile the "
    << "code without\n"
    << "paranoia if you think this is OK...\n"
    << std::endl;
   throw OomphLibError(
    error_stream.str(),
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 // Counter for number of intersections
 unsigned intersect_counter=0;
 
 //Get first vertex
 Vector<double> p1=polygon_vertices[0];
 for (unsigned i=1;i<=nvertex;i++)
  {
   // Get second vertex by wrap-around
   Vector<double> p2 = polygon_vertices[i%nvertex];

   if (point[1] > std::min(p1[1],p2[1]))
    {
     if (point[1] <= std::max(p1[1],p2[1]))
      {
       if (point[0] <= std::max(p1[0],p2[0]))
        {
         if (p1[1] != p2[1])
          {
           double xintersect =
            (point[1]-p1[1])*(p2[0]-p1[0])/
            (p2[1]-p1[1])+p1[0];
           if ( (p1[0] == p2[0]) ||
                (point[0] <= xintersect) )
            {
             intersect_counter++;
            }
          }
        }
      }
    }
   p1 = p2;
  }
 
 // Even number of intersections: outside
 if (intersect_counter%2==0)
  {
   return false;
  } // if (intersect_counter%2==0)
 else
  {
   return true;
  }
 
}

//======================================================================
/// Move the boundary nodes onto the boundary defined by the old mesh
//======================================================================
template <class ELEMENT>
void RefineableTriangleMesh<ELEMENT>::snap_nodes_onto_boundary(
  RefineableTriangleMesh<ELEMENT>* &new_mesh_pt, const unsigned &b)
 {

  // Quick return
  if (!Boundary_coordinate_exists[b])
   {
    return;
   }

  //Firstly we set the boundary coordinates of the new nodes
  //In case the mapping between the geometric object's intrinsic coordinate
  //and the arc-length coordinate is nonlinear. This is only an approximation, 
  //but it will ensure that the nodes that were input to triangle will
  //retain exactly the same boundary coordinates and then linear interpolation
  //is used between those values for any newly created nodes.

  //Create a vector of existing boundary nodes with their boundary
  //coordinate as the first entry so that we can use standard sort algorithms
  Vector<double> node_coord(3);
  Vector<double> b_coord(1);
  const unsigned n_boundary_node = this->nboundary_node(b);

  // Quick return if there are no nodes
  if (n_boundary_node==0)
   {
    return;
   }

  Vector<Vector<double> > old_boundary_node(n_boundary_node);
  for(unsigned n=0;n<n_boundary_node;n++)
   {
    Node* nod_pt = this->boundary_node_pt(b,n);
    nod_pt->get_coordinates_on_boundary(b,b_coord);
    node_coord[0] = b_coord[0];
    node_coord[1] = nod_pt->x(0);
    node_coord[2] = nod_pt->x(1);
    old_boundary_node[n] = node_coord;
   }

  //Sort the vector
  std::sort(old_boundary_node.begin(),old_boundary_node.end());

  //Set up an equivalent ordered vector for the new nodes, based on the
  //current coordinate which is the scaled arc-length.
  //Also provide storage for the original node index,
  //the mapped coordinate and a flag to indicate whether the mapped
  //coordinate has been assigned.
  const unsigned n_new_boundary_node = new_mesh_pt->nboundary_node(b);
  Vector<Vector<double> > new_boundary_node(n_new_boundary_node);
  //There will be six data associated with each node
  node_coord.resize(6,0.0);
  for(unsigned n=0;n<n_new_boundary_node;n++)
   {
    Node* nod_pt = new_mesh_pt->boundary_node_pt(b,n);
    nod_pt->get_coordinates_on_boundary(b,b_coord);
    node_coord[0] = b_coord[0];
    node_coord[1] = nod_pt->x(0);
    node_coord[2] = nod_pt->x(1);
    node_coord[3] = n;
    new_boundary_node[n] = node_coord;
   }

  //Sort the new boundary nodes based on their arc-length coordinate
  std::sort(new_boundary_node.begin(),new_boundary_node.end());

  //We now have two sets of nodes ordered by a coordinate that acts in the
  //same direction and has the same limits.

  //Loop over the vector of new nodes and allocate exactly the same coordinate
  //as the old nodes at points of coincidence
  unsigned old_index = 0;
  for(unsigned n=0;n<n_new_boundary_node;++n)
   {
    //Loop over the set of old nodes and if the x and y coordinates
    //coincide with the new node copy accross the new boundary coordinate
    for(unsigned m=old_index;m<n_boundary_node;++m)
     {
      if(
        (std::fabs(old_boundary_node[m][1] - new_boundary_node[n][1])<1.0e-14)
        &&
        (std::fabs(old_boundary_node[m][2] - new_boundary_node[n][2])<1.0e-14))
       {
        //Store the boundary coordinate from the old mesh
        new_boundary_node[n][4] = old_boundary_node[m][0];
        //Say that it has been stored
        new_boundary_node[n][5] = 1.0;
        //For efficiency, we can start the iteration from here next time round
        //because both vectors are ordered
        old_index = m;
        break;
       }
     }
   }

  //Check that the end-points have new boundary coordinates allocated
#ifdef PARANOID
  if((new_boundary_node[0][5]==0.0) ||
    (new_boundary_node[n_new_boundary_node-1][5] == 0.0))
   {
    std::ostringstream error_stream;
    error_stream <<
    "New boundary coordinates not found for the first and/or last nodes\n"
    << "on the boundary " << b
    << ". This should not happen because these\nlimits should have been setup"
    << "in the constructor\n";
    error_stream <<
    "The distance between the new and old nodes is probably outside\n"
    <<"our tolerance.\n";
    error_stream.precision(20);
    error_stream << "Old boundaries: \n";
    error_stream <<
    old_boundary_node[0][1] << " " << old_boundary_node[0][2]
    << " : " <<
    old_boundary_node[n_boundary_node-1][1] << " " <<
    old_boundary_node[n_boundary_node-1][2] << "\n";
    error_stream << "New boundaries: \n" <<
    new_boundary_node[0][1] << " " << new_boundary_node[0][2] << " : " <<
    new_boundary_node[n_new_boundary_node-1][1] << " " <<
    new_boundary_node[n_new_boundary_node-1][2] << "\n";

    OomphLibWarning(error_stream.str(),
      "RefineableTriangleMesh::snap_nodes_onto_boundary()",
      OOMPH_EXCEPTION_LOCATION);
   }
#endif

  //The end points should always be present, so we
  //can (and must) always add them in exactly
  new_boundary_node[0][4] = new_boundary_node[0][0];

  /// Correct!? Because assigned again below
  new_boundary_node[n_new_boundary_node-1][4] =
  new_boundary_node[0][5] = 1.0;

  new_boundary_node[n_new_boundary_node-1][4] =
  new_boundary_node[n_new_boundary_node-1][0];
  new_boundary_node[n_new_boundary_node-1][5] = 1.0;

  //Create a list of boundary nodes that must be moved
  std::list<unsigned> nodes_to_be_snapped;

  //Now loop over the interior nodes again and 
  //use linear interpolation to fill in any unassigned coordiantes
  for(unsigned n=1;n<n_new_boundary_node-1;++n)
   {
    //If the new boundary coordinate has NOT been allocated
    if(new_boundary_node[n][5]==0.0)
     {
      //Add its (unsorted) node number to the list
      nodes_to_be_snapped.push_back(
        static_cast<unsigned>(new_boundary_node[n][3]));

      //We assume that the previous nodal value has been assigned
      //and read out the old and new boundary coordinates
      double zeta_old_low = new_boundary_node[n-1][0];
      double zeta_new_low = new_boundary_node[n-1][4];

      //Loop over the nodes above the current node until 
      //we find the next one that has been allocated
      for(unsigned m=n+1;m<n_new_boundary_node;++m)
       {
        if(new_boundary_node[m][5]==1.0)
         {
          //Read out the old boundary coordinate
          double zeta_old_high = new_boundary_node[m][0];
          double zeta_new_high = new_boundary_node[m][4];
          //Use linear interpolation to assign the new boundary coordinate
          double frac = (new_boundary_node[n][0] - zeta_old_low)/
          (zeta_old_high - zeta_old_low);
          new_boundary_node[n][4] = zeta_new_low
          + frac*(zeta_new_high - zeta_new_low);
          new_boundary_node[n][5] = 1.0;
          break;
         }
       }
     }
   }

  //Loop over all the nodes and set the new boundary coordinate
  for(unsigned n=0;n<n_new_boundary_node;++n)
   {
    if(new_boundary_node[n][5]==0)
     {
      throw OomphLibError(
        "New boundary coordinate not assigned\n",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
     }

    //get the old coordinate
    new_mesh_pt->boundary_node_pt(
      b,static_cast<unsigned>(new_boundary_node[n][3]))
    ->get_coordinates_on_boundary(b,b_coord);
    //Set the new coordinate
    b_coord[0] = new_boundary_node[n][4];
    new_mesh_pt->boundary_node_pt(
      b,static_cast<unsigned>(new_boundary_node[n][3]))
    ->set_coordinates_on_boundary(b,b_coord);
   }

  Mesh* face_mesh_pt = new Mesh();
  create_unsorted_face_mesh_representation(b, face_mesh_pt);

  //Now that the coordinates have been set up we can do the snapping

  //Now make the mesh as geometric object
  MeshAsGeomObject* mesh_geom_obj_pt = new MeshAsGeomObject(face_mesh_pt);

  //Now assign the new nodes positions based on the old meshes
  //potentially curvilinear boundary (its geom object incarnation)
  Vector<double> new_x(2);

  //Loop over the nodes that need to be snapped
  for(std::list<unsigned>::iterator it=nodes_to_be_snapped.begin();
    it!=nodes_to_be_snapped.end();++it)
   {
    //Read out the boundary node number
    unsigned n = *it;
    //Get the boundary coordinate of all new nodes
    Node* const nod_pt = new_mesh_pt->boundary_node_pt(b,n);
    nod_pt->get_coordinates_on_boundary(b,b_coord);
    //Let's find boundary coordinates of the new node
    mesh_geom_obj_pt->position(b_coord,new_x);
    //Now snap to the boundary
    for(unsigned i=0;i<2;i++)
     {
      nod_pt->x(i) = new_x[i];
     }
   }

  //Delete the allocated memory for the geometric object and face mesh
  delete mesh_geom_obj_pt;
  // Flush the nodes from the face mesh to make sure we
  // don't delete them (the bulk mesh still needs them!)
  face_mesh_pt->flush_node_storage();
  delete face_mesh_pt;

  //Fix up the elements adjacent to the boundary

  // Dummy six node element for sorting out bubble node for
  // seven node enriched quadratic triangles
  TElement<2,3> dummy_six_node_element;
  for (unsigned j=0;j<6;j++)
   {
    dummy_six_node_element.construct_node(j);
   }

  //This should definitely become a triangular element member function
  //Loop over elements 
  unsigned n_bound_el = new_mesh_pt->nboundary_element(b);
  for(unsigned e=0;e<n_bound_el;e++)
   {
    FiniteElement* el_pt = new_mesh_pt->boundary_element_pt(b,e);

    // Deal with different numbers of nodes separately
    unsigned nnod=el_pt->nnode();

// #ifdef PARANOID
//     // Flag to indicate if we successully classified/dealt with the element
//     bool success=false;
// #endif

    // Simplex element: Nothing to be done other than error checking
    if (nnod==3)
     {
#ifdef PARANOID
      // Try to cast to a simplex element
      TElement<2,2>* t_el_pt=dynamic_cast<TElement<2,2>*>(el_pt);
      if (t_el_pt==0)
       {
        throw OomphLibError(
          "Have a three-noded element that's not a TElement<2,2>",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
       }
      // If I get there I must not have thrown :)
      //success=true;
#endif
     }
    // Quadratic element (or enriched quadratic)

    else if ((nnod==6)||(nnod==7))
     {

#ifdef PARANOID
      // Try to cast to a quadratic element
      TElement<2,3>* t_el_pt=dynamic_cast<TElement<2,3>*>(el_pt);
      if (t_el_pt==0)
       {
        if (nnod==6)
         {
          throw OomphLibError(
            "Have a six-noded element that's not a TElement<2,3>",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
         }
        else
         {
          throw OomphLibError(
            "Have a seven-noded element that's not a TElement<2,3>",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
         }
       }
      // If I get there I must not have thrown :)
      //success=true;
#endif  
      // Deal with six noded stuff for all (normal and enriched) elements 

      ///----------------------------------------------------------------
      ///         Repositioning of mid-side nodes
      ///----------------------------------------------------------------

      //Side between 0 and 1
      if(el_pt->node_pt(3)->is_on_boundary(b))
       {
        //Make sure that the node I'm about to move is NOT on
        //a boundary
        if(!el_pt->node_pt(5)->is_on_boundary())
         {
          //Reset the internal nodes
          for(unsigned i=0;i<2;i++)
           {
            el_pt->node_pt(5)->x(i) =
            0.5*(el_pt->node_pt(0)->x(i) + el_pt->node_pt(2)->x(i));
           }
         }
        //Make sure that the node I'm about to move is NOT on
        //a boundary
        if(!el_pt->node_pt(4)->is_on_boundary())
         {
          //Reset the internal nodes
          for(unsigned i=0;i<2;i++)
           {
            el_pt->node_pt(4)->x(i) =
            0.5*(el_pt->node_pt(1)->x(i) + el_pt->node_pt(2)->x(i));
           }
         }
       }

      //Side between 1 and 2
      if(el_pt->node_pt(4)->is_on_boundary(b))
       {
        //Make sure that the node I'm about to move is NOT on
        //a boundary
        if(!el_pt->node_pt(5)->is_on_boundary())
         {
          //Reset the internal nodes
          for(unsigned i=0;i<2;i++)
           {
            el_pt->node_pt(5)->x(i) =
            0.5*(el_pt->node_pt(0)->x(i) + el_pt->node_pt(2)->x(i));
           }
         }
        //Make sure that the node I'm about to move is NOT on
        //a boundary
        if(!el_pt->node_pt(3)->is_on_boundary())
         {
          //Reset the internal nodes
          for(unsigned i=0;i<2;i++)
           {
            el_pt->node_pt(3)->x(i) =
            0.5*(el_pt->node_pt(0)->x(i) + el_pt->node_pt(1)->x(i));
           }
         }
       }

      //Side between 0 and 2
      if(el_pt->node_pt(5)->is_on_boundary(b))
       {
        //Make sure that the node I'm about to move is NOT on
        //a boundary
        if(!el_pt->node_pt(4)->is_on_boundary())
         {
          //Reset the internal nodes
          for(unsigned i=0;i<2;i++)
           {
            el_pt->node_pt(4)->x(i) =
            0.5*(el_pt->node_pt(1)->x(i) + el_pt->node_pt(2)->x(i));
           }
         }
        //Make sure that the node I'm about to move is NOT on
        //a boundary
        if(!el_pt->node_pt(3)->is_on_boundary())
         {
          //Reset the internal nodes
          for(unsigned i=0;i<2;i++)
           {
            el_pt->node_pt(3)->x(i) =
            0.5*(el_pt->node_pt(0)->x(i) + el_pt->node_pt(1)->x(i));
           }
         }
       }

      // If it's seven noded it's likely to be an enriched one: Deal with
      // the central (bubble) node
      if (nnod==7)
       {
        // Try to cast to an enriched quadratic element
        TBubbleEnrichedElement<2,3>* t_el_pt=
        dynamic_cast<TBubbleEnrichedElement<2,3>*>(el_pt);
        if (t_el_pt==0)
         {
          throw OomphLibError(
            "Have seven-noded element that's not a TBubbleEnrichedElement<2,3>",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
         }

        // Assign the new non-bubble coordinates to the six noded dummy element
        for (unsigned j=0;j<6;j++)
         {
          for (unsigned i=0;i<2;i++)
           {
            dummy_six_node_element.node_pt(j)->x(i)=el_pt->node_pt(j)->x(i);
           }
         }

        // Local coordinate of enriched node
        unsigned j_enriched=6;
        Vector<double> s(2);
        el_pt->local_coordinate_of_node(j_enriched,s);

        // Get its position from non-enriched element
        Vector<double> x(2);
        dummy_six_node_element.interpolated_x(s,x);
        el_pt->node_pt(j_enriched)->x(0) = x[0];
        el_pt->node_pt(j_enriched)->x(1) = x[1];
       }
     }
    // Any other case cannot be dealt with at the moment

    else
     {
      std::ostringstream error_stream;
      error_stream
      << "Cannot deal with this particular " << nnod
      << "-noded element yet.\n"
      << "Please implement this yourself.\n";
      throw OomphLibError(error_stream.str(),
OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
     }
   }

  // Cleanup
  for (unsigned j=0;j<6;j++)
   {
    delete dummy_six_node_element.node_pt(j);
   }

 }



#endif // #ifdef OOMPH_HAS_TRIANGLE_LIB

}
#endif
