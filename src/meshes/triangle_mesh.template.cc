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
#include "../rigid_body/rigid_body_elements.h"


namespace oomph
{
 //============================================================
 /// Static empty vector for use as a default argument
 /// to to constructor, which will specify that all 
 /// "holes" remain unfilled by default
 //============================================================
 template<class ELEMENT>
 std::set<unsigned> TriangleMesh<ELEMENT>::Empty_fill_index;
 
//======================================================================
/// Build with the help of the scaffold mesh coming  
/// from the triangle mesh generator Triangle.
//======================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::build_from_scaffold(TimeStepper* time_stepper_pt,
                                                 const bool &use_attributes)
 {   
  // Create space for elements
  unsigned nelem=Tmp_mesh_pt->nelement();
  Element_pt.resize(nelem);
  
   // Create space for nodes
  unsigned nnode_scaffold=Tmp_mesh_pt->nnode();
  
  // Create a map storing the node_id of the mesh used to update the 
  // node position in the update_triangulateio function
  std::map<Node*,unsigned> old_global_number;
  
  // Store the TriangulateIO node id 
  for(unsigned inod=0;inod<nnode_scaffold;inod++)
   {
    Node* old_node_pt=Tmp_mesh_pt->node_pt(inod);
    old_global_number[old_node_pt]=inod;
   }
  
  // Initialize the old node id vector
  Oomph_vertex_nodes_id.resize(nnode_scaffold);
 
  // Create space for nodes     
  Node_pt.resize(nnode_scaffold);

  // Set number of boundaries
  unsigned nbound=Tmp_mesh_pt->nboundary();

  // Resize the boundary information
  set_nboundary(nbound);
  Boundary_element_pt.resize(nbound);
  Face_index_at_boundary.resize(nbound);

  //If we have different regions, then resize the region
  //information
  if(use_attributes)
   {
    Boundary_region_element_pt.resize(nbound);
    Face_index_region_at_boundary.resize(nbound);
   }

  // Loop over elements in scaffold mesh, visit their nodes
  for (unsigned e=0;e<nelem;e++)
   {
    Element_pt[e]=new ELEMENT;
   }
   
  // In the first instance build all nodes from within all the elements
  unsigned nnod_el=Tmp_mesh_pt->finite_element_pt(0)->nnode();
  // Loop over elements in scaffold mesh, visit their nodes
  for (unsigned e=0;e<nelem;e++)
   {
    // Loop over all nodes in element
    for (unsigned j=0;j<nnod_el;j++)
     {
      // Create new node, using the NEW element's construct_node
      // member function
      finite_element_pt(e)->construct_node(j,time_stepper_pt); 
     }
   }

  // Map of Element attribute pairs
  std::map<double,Vector<FiniteElement*> > element_attribute_map;
   
   
  // Setup map to check the (pseudo-)global node number 
  // Nodes whose number is zero haven't been copied across
  // into the mesh yet.
  std::map<Node*,unsigned> global_number;
  unsigned global_count=0;

  // Loop over elements in scaffold mesh, visit their nodes
  for (unsigned e=0;e<nelem;e++)
   {
    // Loop over all nodes in element
    for (unsigned j=0;j<nnod_el;j++)
     {
      // Pointer to node in the scaffold mesh
      Node* scaffold_node_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(j);

      // Get the (pseudo-)global node number in scaffold mesh
      // (It's zero [=default] if not visited this one yet)
      unsigned j_global=global_number[scaffold_node_pt];

      // Haven't done this one yet
      if (j_global==0)
       {
        // Find and store the node_id in the old nodes map
        Oomph_vertex_nodes_id[global_count]=
         old_global_number[scaffold_node_pt];
        
        // Give it a number (not necessarily the global node 
        // number in the scaffold mesh -- we just need something
        // to keep track...)
        global_count++;
        global_number[scaffold_node_pt]=global_count;

        // Copy new node, created using the NEW element's construct_node
        // function into global storage, using the same global
        // node number that we've just associated with the 
        // corresponding node in the scaffold mesh
        Node_pt[global_count-1]=finite_element_pt(e)->node_pt(j);

        // Assign coordinates
        for(unsigned i=0;i<finite_element_pt(e)->dim();i++)
         {
          Node_pt[global_count-1]->x(i)=scaffold_node_pt->x(i);
         }

         
        // Get pointer to set of mesh boundaries that this 
        // scaffold node occupies; NULL if the node is not on any boundary
        std::set<unsigned>* boundaries_pt;
        scaffold_node_pt->get_boundaries_pt(boundaries_pt);

        // Loop over the mesh boundaries that the node in the scaffold mesh
        // occupies and assign new node to the same ones.
        if (boundaries_pt!=0)
         {
          this->convert_to_boundary_node(Node_pt[global_count-1]);
          for(std::set<unsigned>::iterator it=boundaries_pt->begin();
              it!=boundaries_pt->end();++it)
           {
            add_boundary_node(*it,Node_pt[global_count-1]);
           }
         }
       }
      // This one has already been done: Kill it
      else
       {
        // Kill it
        delete finite_element_pt(e)->node_pt(j);

        // Overwrite the element's pointer to local node by
        // pointer to the already existing node -- identified by
        // the number kept in the map
        finite_element_pt(e)->node_pt(j)=Node_pt[j_global-1];         
       }
     }

    if(use_attributes)
     {
      element_attribute_map[Tmp_mesh_pt->element_attribute(e)].push_back(
       finite_element_pt(e));
     }

   }

  //Now let's construct lists
  //Find the number of attributes
  if(use_attributes)
   {
    unsigned n_attribute = element_attribute_map.size();
    //There are n_attribute different regions
    Region_element_pt.resize(n_attribute);
    Region_attribute.resize(n_attribute);
    //Copy the vectors in the map over to our internal storage
    unsigned count = 0;
    for(std::map<double,Vector<FiniteElement*> >::iterator it =
         element_attribute_map.begin(); it != element_attribute_map.end();++it)
     {
      Region_attribute[count] = it->first;
      Region_element_pt[count] = it->second;
      ++count;
     }
   }

  // At this point we've created all the elements and 
  // created their vertex nodes. Now we need to create
  // the additional (midside and internal) nodes!


  // We'll first create all local nodes for all elements
  // and then delete the superfluous ones that have
  // a matching node in an adjacent element.

  unsigned boundary_id;   

  // Get number of nodes along element edge and dimension of element (2)
  // from first element
  unsigned nnode_1d=finite_element_pt(0)->nnode_1d();
  unsigned dim=finite_element_pt(0)->dim();

  // Storage for the local coordinate of the new node
  Vector<double> s(dim);

  // Get number of nodes in the element from first element
  unsigned nnode=finite_element_pt(0)->nnode();

  // Loop over all elements
  for (unsigned e=0;e<nelem;e++)
   {
    FiniteElement* const elem_pt = finite_element_pt(e);
    FiniteElement* const tmp_elem_pt = Tmp_mesh_pt->finite_element_pt(e);
    // Loop over the new nodes in the element and create them.
    for(unsigned j=3;j<nnode;j++)
     {

      // Create new node
      Node* new_node_pt=elem_pt->construct_node(j,time_stepper_pt);

      // What are the node's local coordinates?
      elem_pt->local_coordinate_of_node(j,s);

      // Find the coordinates of the new node from the existing
      // and fully-functional element in the scaffold mesh
      for(unsigned i=0;i<dim;i++)
       {
        new_node_pt->x(i)= tmp_elem_pt->interpolated_x(s,i);
       }
   
      // Searching if the new node is on a boundary,
      // and if it is, add the zero-based boundary id to it

      // The general convention is that midside nodes are numbered
      // consecutively, in anti-clockwise direction; they are numbered 
      // after the corner nodes and before the internal nodes. 


      // These are the edge nodes on the element's edge 0:
      if((2<j)&&(j<nnode_1d+1))
       {
        boundary_id=Tmp_mesh_pt->edge_boundary(e,0);
        if(boundary_id>0)
         {
          this->convert_to_boundary_node(new_node_pt);
          add_boundary_node(boundary_id-1,new_node_pt);
         }          
       }
      // These are the edge nodes on the element's edge 1:
      if((nnode_1d<j)&&(j<2*nnode_1d-1))
       {
        boundary_id=Tmp_mesh_pt->edge_boundary(e,1); 
        if(boundary_id>0)
         {
          this->convert_to_boundary_node(new_node_pt);
          add_boundary_node(boundary_id-1,new_node_pt);
         } 
       }
      // These are the edge nodes on the element's edge 2:
      if((2*(nnode_1d-1)<j) && (j<3*(nnode_1d-1)))
       {
        boundary_id=Tmp_mesh_pt->edge_boundary(e,2); 
        if(boundary_id>0)
         {
          this->convert_to_boundary_node(new_node_pt);
          add_boundary_node(boundary_id-1,new_node_pt);
         } 
       }

     } // end of loop over new nodes

    //Set up the boundary element information
    //Loop over the edges
    for(unsigned j=0;j<3;j++)
     {
      boundary_id = Tmp_mesh_pt->edge_boundary(e,j);
      if(boundary_id > 0)
       {
        Boundary_element_pt[boundary_id-1].push_back(elem_pt);
        //Need to put a shift in here because of an inconsistent naming 
        //convention between triangle and face elements
        Face_index_at_boundary[boundary_id-1].push_back((j+2)%3);

        //If using regions set up the boundary information
        if(use_attributes)
         {
          //Element adjacent to boundary
          Boundary_region_element_pt[boundary_id-1]
           [static_cast<unsigned>(Tmp_mesh_pt->element_attribute(e))].
           push_back(elem_pt);
          //Need to put a shift in here because of an inconsistent naming 
          //convention between triangle and face elements
          Face_index_region_at_boundary[boundary_id-1]
           [static_cast<unsigned>(Tmp_mesh_pt->element_attribute(e))].
           push_back((j+2)%3);
         }
       }
     }
     
   } //end of loop over elements
    
   
  // Lookup scheme has now been setup yet
  Lookup_for_elements_next_boundary_is_setup=true;

  // If condition to check the segment local node number (all located on edges)
        
  if(nnode_1d >4)
   {
    oomph_info << "===================================================="
               << std::endl<<std::endl;
    oomph_info << "Warning -- using a terribly inefficient scheme to " 
               << std::endl;
    oomph_info << "determine duplicated edge nodes" << std::endl;
    oomph_info << std::endl << std::endl;
   
  
    // Tolerance for detecting repeated nodes
    double tolerance=1.0e-10;
     
    // Loop over elements
    for (unsigned e=0;e<nelem;e++)
     {
     
      // Loop over new local nodes
      for(unsigned j=3;j<nnode;j++)
       {
        // Pointer to the element's local node
        Node* node_pt=finite_element_pt(e)->node_pt(j);
       
        // By default, we assume the node is new
        bool is_new=true;
       
        // Now loop over all nodes already stored in the
        // Mesh's node list and check if the local node considered
        // at the moment is located at the same position (within
        // the tolerance specified above). Note: We start searching
        // at the number of nodes that created before
        // we added any new (midside and other) nodes. 
        unsigned nnod=Node_pt.size();
        for(unsigned k=nnode_scaffold;k<nnod;k++)
         { 
          double a=node_pt->x(0);
          double b=node_pt->x(1);
          double c=Node_pt[k]->x(0);
          double d=Node_pt[k]->x(1);
          // Nodes coincide
          if(sqrt((c-a)*(c-a)+(d-b)*(d-b))<tolerance)
           {
            //If the node is on a boundary                    
            //Need to remove the node from the boundary list in the mesh
            //before deletion  
            if(node_pt->is_on_boundary()) 
             {     
              std::set<unsigned>* boundaries_pt;                            
              //Find the boundaries on which the nodes live                   
              node_pt->get_boundaries_pt(boundaries_pt);     

              //Now copy the boundary indices into local storage because the
              //set will change as the node is removed from the boundaries 
              const unsigned size = boundaries_pt->size();       
              unsigned boundary_index[size];                     
              unsigned counter=0;                        
              for(std::set<unsigned>::iterator it = boundaries_pt->begin();
                  it != boundaries_pt->end();++it)         
               {        
                boundary_index[counter] = *it;              
                ++counter;                               
               }      
              
              //Loop over these boundaries and remove the node from them
              for(unsigned i=0;i<size;i++)               
               {                                             
                this->remove_boundary_node(boundary_index[i],node_pt);    
               }   
             }

            // Delete local node in element...
            delete finite_element_pt(e)->node_pt(j);
            // ... and reset pointer to the existing node
            finite_element_pt(e)->node_pt(j)=Node_pt[k];
            // Node is not new!
            is_new=false;
            break;
           }
         }
        // If the local node considered at the moment is new, add it to 
        // the mesh's vector of nodes
        if(is_new)
         {
          Node_pt.push_back(node_pt);
         }
       }
     }

    oomph_info << "Done! If this took too long for you, reimplement it... " 
               << std::endl;
    oomph_info << "===================================================="
               << std::endl<<std::endl;    

   }
  // Use efficient scheme for 3 nodes for segment
  else if(nnode_1d==3)
   {     
    // Map storing the mid-side of an edge; edge identified by
    // pointers to vertex nodes
    MapMatrix<Node*,Node*> central_edge_node_pt;
    Node* edge_node1_pt=0;
    Node* edge_node2_pt=0;
    
    // Loop over elements
    for (unsigned e=0;e<nelem;e++)
     {
      //For standard quadratic elements, all nodes are edge nodes
      unsigned n_edge_node = nnode;
      //If we have enriched elements, we need to treat the central
      //node slightly differently
      if(nnode==7) 
       {
        // Pointer to the element's central node. This 
        // definitely gets added because it cannot be
        // duplicated.
        Node* node_pt=finite_element_pt(e)->node_pt(6);
        Node_pt.push_back(node_pt);
        //There are only n-1 edge nodes
        n_edge_node -= 1;
       }

      // Loop over new local nodes
      for(unsigned j=3;j<n_edge_node;j++)
       {
        // Pointer to the element's local node
        Node* node_pt=finite_element_pt(e)->node_pt(j);
       
       
        // Switch on local node number (all located on edges)
        switch (j)
         {        
          // Node 3 is located between nodes 0 and 1
         case 3:
         
          edge_node1_pt=finite_element_pt(e)->node_pt(0);
          edge_node2_pt=finite_element_pt(e)->node_pt(1);
          break;
          // Node 4 is located between nodes 1 and 2
         case 4:
         
          edge_node1_pt=finite_element_pt(e)->node_pt(1);
          edge_node2_pt=finite_element_pt(e)->node_pt(2);
          break;
          // Node 6 is located between nodes 2 and 0
         case 5:
         
          edge_node1_pt=finite_element_pt(e)->node_pt(2);
          edge_node2_pt=finite_element_pt(e)->node_pt(0);
          break;

         default:
          //Error
          throw OomphLibError("More than six or seven nodes in TriangleMesh",
                              "TriangleMesh::build_from_scaffold()",
                              OOMPH_EXCEPTION_LOCATION);
         }

        if (central_edge_node_pt(edge_node1_pt,edge_node2_pt)==0)
         {
          // This is superfluous: central_edge_node_pt(edge_node1_pt,
          // edge_node2_pt)=node_pt;
          central_edge_node_pt(edge_node2_pt,edge_node1_pt)=node_pt;
          Node_pt.push_back(node_pt);
         }    
        else
         {
          //If the node is on a boundary                    
          //Need to remove the node from the boundary list in the mesh
          //before deletion  
          if(node_pt->is_on_boundary()) 
           {     
            std::set<unsigned>* boundaries_pt;                            
            //Find the boundaries on which the nodes live                   
            node_pt->get_boundaries_pt(boundaries_pt);     
            
            //Now copy the boundary indices into local storage because the
            //set will change as the node is removed from the boundaries 
            const unsigned size = boundaries_pt->size();       
            unsigned boundary_index[size];                     
            unsigned counter=0;                        
            for(std::set<unsigned>::iterator it = boundaries_pt->begin();
                it != boundaries_pt->end();++it)         
             {        
              boundary_index[counter] = *it;              
              ++counter;                               
             }      
            
            //Loop over these boundaries and remove the node from them
            for(unsigned i=0;i<size;i++)               
             {                                             
              this->remove_boundary_node(boundary_index[i],node_pt);    
             }   
           }
          
          // Delete local node in element...
          delete finite_element_pt(e)->node_pt(j);
          
          // ... and reset pointer to the existing node
          finite_element_pt(e)->node_pt(j)=
           central_edge_node_pt(edge_node1_pt,edge_node2_pt);
         }
       
       }
     }   
   }

  // Use efficient scheme for 4 nodes for segment
  else if(nnode_1d==4)
   {
   
    // Map storing the mid-side of an edge; edge identified by
    // pointers to vertex nodes
    MapMatrix<Node*,Node*> first_edge_node_pt ;
    MapMatrix<Node*,Node*> second_edge_node_pt;
    Node* edge_node1_4d_pt=0;
    Node* edge_node2_4d_pt=0;
   
   
    // Loop over elements
    for (unsigned e=0;e<nelem;e++)
     {
     
      // Pointer to the element's central node. This 
      // definitely gets added because it cannot be
      // duplicated.
      Node* node_pt=finite_element_pt(e)->node_pt(9);
      Node_pt.push_back(node_pt);
     
      // Loop over edge nodes to check if they are duplicated
      for(unsigned j=3;j<nnode-1;j+=2)
       {
        // Pointer to the element's local node
        Node* node_pt=finite_element_pt(e)->node_pt(j);
       
        // Switch on local node number (all located on edges)
        switch (j)
         {
         
          // Node 3 and 4 are located between nodes 0 and 1
         case 3:
         
          edge_node1_4d_pt=finite_element_pt(e)->node_pt(0);
          edge_node2_4d_pt=finite_element_pt(e)->node_pt(1);        
          break;
          // Node 5 and 6 are located between nodes 1 and 2
         case 5:
         
          edge_node1_4d_pt=finite_element_pt(e)->node_pt(1);
          edge_node2_4d_pt=finite_element_pt(e)->node_pt(2);        
          break;
          // Node 7 and 8 are located between nodes 2 and 0
         case 7:
         
          edge_node1_4d_pt=finite_element_pt(e)->node_pt(2);
          edge_node2_4d_pt=finite_element_pt(e)->node_pt(0);          
          break;

         default:
          //Error 
          throw OomphLibError("More than ten nodes in TriangleMesh, now loop",
                              "TriangleMesh::build_from_scaffold()",
                              OOMPH_EXCEPTION_LOCATION);
         }
       
        if (second_edge_node_pt(edge_node1_4d_pt,edge_node2_4d_pt)==0)
         {
          first_edge_node_pt(edge_node2_4d_pt,edge_node1_4d_pt)=node_pt;
          Node_pt.push_back(node_pt);
         
          node_pt=finite_element_pt(e)->node_pt(j+1);
          second_edge_node_pt(edge_node2_4d_pt,edge_node1_4d_pt)=node_pt;
          Node_pt.push_back(node_pt);
         }
        else
         {

          // Delete local node in element...
          delete finite_element_pt(e)->node_pt(j);
          delete finite_element_pt(e)->node_pt(j+1);
         
          // ... and reset pointer to the existing nodes
          finite_element_pt(e)->node_pt(j)=
           second_edge_node_pt(edge_node1_4d_pt,edge_node2_4d_pt);
          finite_element_pt(e)->node_pt(j+1)=
           first_edge_node_pt(edge_node1_4d_pt,edge_node2_4d_pt);
         
         } 
       
       }      
     }  
   }             
 }



   
 
//======================================================================
/// Setup boundary coordinate on boundary b. Doc Faces
/// in outfile. Boundary coordinate increases continously along
/// polygonal boundary. It's zero at the lexicographically
/// smallest node on the boundary.
//======================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::setup_boundary_coordinates(const unsigned& b,
                                                        std::ofstream& outfile)
 {
  // Temporary storage for face elements
  Vector<FiniteElement*> face_el_pt;

  // Temporary storage for number of elements adjacent to the boundary
  unsigned nel = 0;

  //If there is more than one region then only use 
  //boundary coordinates from the bulk side (region 0)
  if(this->nregion() > 1)
   {
    // Loop over all elements on boundaries from region 0
    nel=this->nboundary_element_in_region(b,0);
    
#ifdef PARANOID
    if (nel==0)
     {
      std::ostringstream error_message;
      error_message 
       << "No boundary elements in region 0. You're in trouble.\n" 
       << "Ask Andrew to sort this out!\n";
      throw OomphLibError(error_message.str(),
                          "TriangleMesh::setup_boundary_coordinates()",
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    // Loop over the bulk elements adjacent to boundary b
    for(unsigned e=0;e<nel;e++)
     {
      // Get pointer to the bulk element that is adjacent to boundary b
      // in region 0
      FiniteElement* bulk_elem_pt = 
       this->boundary_element_pt_in_region(b,0,e);
      
      //Find the index of the face of element e along boundary b
      // in region 0
      int face_index = this->face_index_at_boundary_in_region(b,0,e);
      
      // Create new face element 
      face_el_pt.push_back(new DummyFaceElement<ELEMENT>(
                            bulk_elem_pt,face_index));   
      
      // Output faces?
      if (outfile.is_open()) 
       {
        face_el_pt[face_el_pt.size()-1]->output(outfile); 
       }
     }
   }
  //Otherwise it's just the normal boundary functions
  else
   {
    // Loop over all elements on boundaries
    nel=this->nboundary_element(b);
    
    // Loop over the bulk elements adjacent to boundary b
    for(unsigned e=0;e<nel;e++)
     {
      // Get pointer to the bulk element that is adjacent to boundary b
      FiniteElement* bulk_elem_pt = this->boundary_element_pt(b,e);
      
      //Find the index of the face of element e along boundary b
      int face_index = this->face_index_at_boundary(b,e);
      
      // Create new face element 
      face_el_pt.push_back(new DummyFaceElement<ELEMENT>(
                            bulk_elem_pt,face_index));   
      
      // Output faces?
      if (outfile.is_open()) 
       {
        face_el_pt[face_el_pt.size()-1]->output(outfile); 
       }
     }  
   }
  
  //Only bother to do anything else, if there are elements
  if(nel > 0)
   {
    // Put first element into ordered list
    std::list<FiniteElement*> ordered_el_pt;
    FiniteElement* el_pt=face_el_pt[0];
    ordered_el_pt.push_back(el_pt);
    
    // Count elements that have been done
    unsigned count_done=0;
    
    // Keep track of who's done
    std::map<FiniteElement*,bool> done_el;
    
    // Keep track of which element is inverted
    std::map<FiniteElement*,bool> is_inverted;
    
    // Number of nodes
    unsigned nnod=el_pt->nnode();

    // Fit in the other elements in at most nel^2 loops
    for (unsigned ee=1;ee<nel;ee++)
     {
      // Loop over all elements to check if they fit to the right
      // or the left of the current one
      for (unsigned e=1;e<nel;e++)
       {
        // Candidate element
        el_pt=face_el_pt[e];

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
    if (count_done!=(nel-1))
     {
      std::ostringstream error_message;
      error_message << "Was only able to setup boundary coordinate on " 
                    << "boundary " << b << "\nfor " << count_done 
                    << " of " << nel << " face elements. This usually means\n"
                    << "that the boundary is not simply connected.\n\n"
                    << "Re-run the setup_boundary_coordintes() function\n"
                    << "with an output file specified "
                    << "as the second argument.\n"
                    << "This will file will contain FaceElements that\n"
                    << "oomph-lib believes to be located on the boundary.\n"
                    << std::endl;
      throw OomphLibError(error_message.str(),
                          "TriangleMesh::setup_boundary_coordinates()",
                          OOMPH_EXCEPTION_LOCATION);
     }

    // First node
    FiniteElement* first_el_pt=*ordered_el_pt.begin();
    Node* first_node_pt=first_el_pt->node_pt(0);
    if (is_inverted[first_el_pt]) first_node_pt=first_el_pt->node_pt(nnod-1);  

    //Last node
    FiniteElement* last_el_pt = ordered_el_pt.back();
    Node* last_node_pt = last_el_pt->node_pt(nnod-1);
    if (is_inverted[last_el_pt]) last_node_pt=last_el_pt->node_pt(0);  

     
    // Coordinates of left node
    double x_left=first_node_pt->x(0);
    double y_left=first_node_pt->x(1);
     
    // Initialise boundary coordinate
    Vector<double> zeta(1,0.0);
         
    // Set boundary coordinate
    first_node_pt->set_coordinates_on_boundary(b,zeta);

    // Lexicographically bottom left node
    std::set<Node*> all_nodes_pt;
    all_nodes_pt.insert(first_node_pt);

    // Now loop over nodes in order
    for (std::list<FiniteElement*>::iterator it=ordered_el_pt.begin();
         it!=ordered_el_pt.end();it++)
     {
      // Get element
      FiniteElement* el_pt=*it;
      
      // Start node and increment
      unsigned k_nod=1;
      int nod_diff=1;
      if (is_inverted[el_pt])
       {
        k_nod=nnod-2;
        nod_diff=-1;
       }
       
      // Loop over nodes
      for (unsigned j=1;j<nnod;j++)
       {
        Node* nod_pt=el_pt->node_pt(k_nod);
        k_nod+=nod_diff;
         
        // Coordinates of right node
        double x_right=nod_pt->x(0);
        double y_right=nod_pt->x(1);
         
        // Increment boundary coordinate
        zeta[0]+=sqrt((x_right-x_left)*(x_right-x_left)+
                      (y_right-y_left)*(y_right-y_left));
                  
        // Set boundary coordinate
        nod_pt->set_coordinates_on_boundary(b,zeta);
         
        // Increment reference coordinate
        x_left=x_right;
        y_left=y_right;
         
        // Get lexicographically bottom left node but only 
        // use vertex nodes as candidates
        all_nodes_pt.insert(nod_pt);
       }
     }

    //The nodes have been assigned arc-length coordinates
    //from one end or the other of the connected segement.

    //If the boundary has a geometric object representation then
    //scale the coordinates to match those of the geometric object
    GeomObject* const geom_object_pt = this->boundary_geom_object_pt(b);
    if(geom_object_pt!=0)
     {
      oomph_info << "A geometric_object is set up for boundary " << b << "\n";
      Vector<double> bound_coord_limits = this->boundary_coordinate_limits(b);
      
      //Get the position of the ends of the geometric object
      Vector<double> zeta(1);
      Vector<double> first_geom_object_location(2);
      Vector<double> last_geom_object_location(2);
      zeta[0] = bound_coord_limits[0];
      geom_object_pt->position(zeta,first_geom_object_location);
      zeta[0] = bound_coord_limits[1];
      geom_object_pt->position(zeta,last_geom_object_location);

      //Calculate the errors in position between the first and last nodes
      //and the endpoints of the geometric object
      double error=0.0;
      double tmp_error = 0.0;
      for(unsigned i=0;i<2;i++)
       {
        const double dist =
         first_geom_object_location[i] - first_node_pt->x(i); 
        tmp_error += dist*dist;
       }
      error += sqrt(tmp_error);
      tmp_error=0.0;
      for(unsigned i=0;i<2;i++)
       {
        const double dist  =
         last_geom_object_location[i] - last_node_pt->x(i);
        tmp_error += dist*dist;
       }
      error += sqrt(tmp_error);
      
      //Calculate the errors in position between the first and last nodes
      //and the endpoints of the geometric object if reversed
      double rev_error=0.0;
      tmp_error = 0.0;
      for(unsigned i=0;i<2;i++)
       {
        const double dist  =
         first_geom_object_location[i] - last_node_pt->x(i);
        tmp_error += dist*dist;
       }
      rev_error += sqrt(tmp_error);
      tmp_error=0.0;
      for(unsigned i=0;i<2;i++)
       {
        const double dist  =
         last_geom_object_location[i] - first_node_pt->x(i);
        tmp_error += dist*dist;
       }
      rev_error += sqrt(tmp_error);
      
      //If the (normal) error is small than reversed then we have the 
      //coordinate direction correct.
      //If not then we must reverse it
      if(error < rev_error)
       {
        oomph_info << "Coordinates are aligned\n";
       }
      else if(error > rev_error)
       {
        oomph_info << "Coordinates are reversed\n";
        //Reverse the boundary coordinates
        double temp = bound_coord_limits[0];
        bound_coord_limits[0] = bound_coord_limits[1];
        bound_coord_limits[1] = temp;
       }
      else
       {
        std::ostringstream error_stream;
        error_stream << "Something very strange has happened.\n" << 
         "The error between the endpoints of the geometric object\n" <<
         "and the first and last nodes on the boundary is the same\n" <<
         "irrespective of the direction of the coordinate.\n" << 
         "This probably means that things are way off.\n" << 
         "The errors are " << error  << " and " << rev_error << "\n";
        throw OomphLibError(error_stream.str(),
                            "TriangleMesh::setup_boundary_coordinates()",
                            OOMPH_EXCEPTION_LOCATION);
       }

      //Get the total arclength of the edge
      last_node_pt->get_coordinates_on_boundary(b,zeta);
      double zeta_old_range = zeta[0]; 
      double zeta_new_range = bound_coord_limits[1] - bound_coord_limits[0];

      //Now scale the coordinates accordingly
      for (std::set<Node*>::iterator it=all_nodes_pt.begin();
           it!=all_nodes_pt.end();it++)
       {    
        Node* nod_pt=(*it);
        nod_pt->get_coordinates_on_boundary(b,zeta);
        zeta[0] = bound_coord_limits[0] + 
         (zeta_new_range/zeta_old_range)*zeta[0];
        nod_pt->set_coordinates_on_boundary(b,zeta);
       }
     }
    else
     {
      oomph_info << "No geometric object for boundary " << b << "\n";
      
      //Only use end points of the whole segment and pick the bottom left
      //node
      Node* bottom_left_node_pt=first_node_pt;
      if (last_node_pt->x(1)<bottom_left_node_pt->x(1))
       {
        bottom_left_node_pt=last_node_pt;
       }
      else if (last_node_pt->x(1)==bottom_left_node_pt->x(1))
       {
        if (last_node_pt->x(0)<bottom_left_node_pt->x(0))            
         {
          bottom_left_node_pt=last_node_pt;
         }
       }
      
      
      // Now adjust boundary coordinate so that the bottom left node
      // has a boundary coordinate of zero and that zeta increases
      // away from that point
      bottom_left_node_pt->get_coordinates_on_boundary(b,zeta);
      double zeta_ref = zeta[0];
      double zeta_max = 0.0;
      for (std::set<Node*>::iterator it=all_nodes_pt.begin();
           it!=all_nodes_pt.end();it++)
       {    
        Node* nod_pt=(*it);
        nod_pt->get_coordinates_on_boundary(b,zeta);
        zeta[0]-=zeta_ref;
        //If direction is reversed, then take absolute value
        if(zeta[0] < 0.0)
         {
          zeta[0] = std::abs(zeta[0]);
         }
        if(zeta[0] > zeta_max) {zeta_max = zeta[0];}
        nod_pt->set_coordinates_on_boundary(b,zeta);
       }

      //Scale all surface coordinates by the max
      for (std::set<Node*>::iterator it=all_nodes_pt.begin();
           it!=all_nodes_pt.end();it++)
       {    
        Node* nod_pt=(*it);
        nod_pt->get_coordinates_on_boundary(b,zeta);
        zeta[0] /= zeta_max;
        nod_pt->set_coordinates_on_boundary(b,zeta);
       }
     }
      
    // Cleanup
    for(unsigned e=0;e<nel;e++)
     {
      delete face_el_pt[e];
      face_el_pt[e]=0;
     }  
    
   }
  
  // Indicate that boundary coordinate has been set up
  Boundary_coordinate_exists[b]=true;
 
 }

#ifdef OOMPH_HAS_TRIANGLE_LIB  

//======================================================================
/// Create TriangulateIO object from TriangleMeshPolygon and 
/// TriangleMeshInternalPolygon
//======================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::build_triangulateio(
  TriangleMeshPolygon* &outer_boundary_pt,
  Vector<TriangleMeshInternalPolygon*> &internal_polygon_pt,
  std::set<unsigned> &fill_index,
  TriangulateIO& triangulate_io)
 {
  // triangulate_io initialization
  TriangleHelper::initialise_triangulateio(triangulate_io);
   
  // Build the TriangulateIO in object

  // Get number of polyline
  unsigned n_boundline=0;
  n_boundline = outer_boundary_pt->npolyline();
  
  // Count the global polygon boundary polyline
  unsigned n_boundglobalseg=0;

  // Count the global holes noundary polyline
  unsigned n_holes = internal_polygon_pt.size();
  
  // Initialize global hole boundary counter
  unsigned multi_hole_polyline=0;

  // Deal with holes
  for(unsigned count_hole=0;count_hole<n_holes;count_hole++)
   {
    unsigned inner_hole_polyline=internal_polygon_pt[count_hole]->npolyline();
    
    // Add value to the global polyline counter
    multi_hole_polyline+=inner_hole_polyline;
   }
  
  // Counter initialization
  unsigned n_holepolyline=0;
  unsigned n_globalvertices=0;
  unsigned n_polylinevertices=0;
  unsigned count_tri=0;
  int edge_segment=1;
     
  // We need to know the dimension of the global number of vertices
  // and segments (basically the same, if the input values aren't wrong)
  // to initialize the pointlist, segmentlist and segmentmarkerlis.
     
  // Loop over the boundary polylines 
  for(unsigned count_seg=0;count_seg<n_boundline;count_seg++)
   {
    n_globalvertices += outer_boundary_pt->polyline_pt(count_seg)->
     nvertex()-1;
   }
  // Loop over the holes boundary
  for(unsigned count_hole=0;count_hole<n_holes;count_hole++)
   {
    n_holepolyline=internal_polygon_pt[count_hole]->npolyline();
    for(unsigned count_seg=0;count_seg<n_holepolyline;count_seg++)
     {
      n_globalvertices += internal_polygon_pt[count_hole]->
       polyline_pt(count_seg)->nvertex()-1;
     }
   }
  // If there's just one boundary. All the vertices should be counted
  // for the polygon...
  if(n_boundline==1)
   {
    n_globalvertices +=1;
   }

  // ...and for the hole.
  for(unsigned count_hole=0;count_hole<n_holes;count_hole++)
   {
    n_holepolyline=internal_polygon_pt[count_hole]->npolyline();
    if(n_holepolyline==1)
     {
      n_globalvertices +=1;
     }
   }

  // Store the global number of vertices and segments
  // in the list. They are supposed to be the same!!
  triangulate_io.numberofpoints=n_globalvertices;
  triangulate_io.numberofsegments=n_globalvertices;
      
  // Prepearing the TriangulateIO objects to store the values     
  triangulate_io.pointlist = 
   (double *) malloc(triangulate_io.numberofpoints * 2 * sizeof(double));
  triangulate_io.segmentlist = 
   (int *) malloc(triangulate_io.numberofsegments * 2 * sizeof(int));
  triangulate_io.segmentmarkerlist = 
   (int *) malloc(triangulate_io.numberofsegments * sizeof(int));
  
  // Storing all the values in the list
  for(unsigned count_seg=0;count_seg<n_boundline;count_seg++)
   {
    // Storing the number of the vertices
    n_polylinevertices = 
     outer_boundary_pt->polyline_pt(count_seg)->nvertex()-1;

    // If there's just one boundary. All the vertices have to be counted   
    if(n_boundline==1)
     {
      n_polylinevertices +=1;
     }
    
    // Store the segmen_id if given by the user
    unsigned idpolyline=outer_boundary_pt->polyline_pt(count_seg)
     ->boundary_id();
    
    // Initialize to zero the segment counter
    unsigned count_bound_segment=0;
    
    //Storing the coordinates for each points
    for(unsigned count_vertices=0;count_vertices<n_polylinevertices;
        count_vertices++)
     {

      triangulate_io.pointlist[count_tri]= outer_boundary_pt->polyline_pt
       (count_seg)
       ->vertex_coordinate(count_vertices)[0];
      triangulate_io.pointlist[count_tri+1]= outer_boundary_pt->
       polyline_pt(count_seg)->vertex_coordinate(count_vertices)[1];

      // Store the segment values
      // If the segment is not the last one, take the next node
      if(count_seg==(n_boundline-1) && 
         count_vertices==(n_polylinevertices-1))
       {
        triangulate_io.segmentlist[count_tri]=edge_segment;
        triangulate_io.segmentlist[count_tri+1]=1;
       }
      else
       {
        triangulate_io.segmentlist[count_tri]=edge_segment;
        triangulate_io.segmentlist[count_tri+1]=edge_segment+1;
        edge_segment++;
       }

      // Store the marker list of the segments
      triangulate_io.segmentmarkerlist[count_tri/2]=idpolyline+1;

      // Increment counter
      count_bound_segment++;
      count_tri+=2;
      n_boundglobalseg++;
     }
   }

  // Initialize  global hole boundary counter
  unsigned nmulti_hole_bound=0;
  
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
      n_polylinevertices = internal_polygon_pt[count_hole]->
       polyline_pt(count_seg)->nvertex()-1;
      
      // If there's just one boundary. All the vertices should be counted   
      if(n_holepolyline==1)
       {
        n_polylinevertices +=1;
       }

      // Store the segment_id if given by the user
      unsigned idpolyline=internal_polygon_pt[count_hole]->
       polyline_pt(count_seg)->boundary_id();
      
      // Initialize 
      unsigned count_bound_segment=0;

      // Store the coordinates for each points
      for(unsigned count_vertices=0;count_vertices<n_polylinevertices;
          count_vertices++)
       {
        triangulate_io.pointlist[count_tri]= internal_polygon_pt[count_hole]->
         polyline_pt(count_seg)->
         vertex_coordinate(count_vertices)[0];
        triangulate_io.pointlist[count_tri+1]=internal_polygon_pt[count_hole]->
         polyline_pt(count_seg)->
         vertex_coordinate(count_vertices)[1];

        // Store the segments values
        // If the segment is not the last one, take the next node
        if(count_seg==(n_holepolyline-1) && 
           count_vertices==(n_polylinevertices-1))
         {
          triangulate_io.segmentlist[count_tri]=edge_segment;
          triangulate_io.segmentlist[count_tri+1]=hole_vertex_start;
         }
        else
         {
          triangulate_io.segmentlist[count_tri]=edge_segment;
          triangulate_io.segmentlist[count_tri+1]=edge_segment+1;
          edge_segment++;
         }

        // Store the marker list of the segments
        triangulate_io.segmentmarkerlist[count_tri/2]= idpolyline+1;
        
        // Increment counter
        count_bound_segment++;
        count_tri+=2;
        n_boundglobalseg++;
       }
     }
    // Add the previous number of boundary
    nmulti_hole_bound+=n_holepolyline;
   }

  // Check if the number of segments and vertices are the the same
  if(n_boundglobalseg!=n_globalvertices)
   {
    std::ostringstream error_stream;
    error_stream  << "Error building TriangulateIO object\n"
                  << "Please, check TriangleMeshPolyLine and\n"
                  << "TriangleMeshPolygon provided"
                  <<std::endl;      
    throw OomphLibError(error_stream.str(),
                        "TriangleMeshBoundary::TriangleMeshBoundary()",
                        OOMPH_EXCEPTION_LOCATION);
   }


  //Find the number of filled holes
  unsigned n_filled = fill_index.size();

  //Fill in the internal regions if needed
  if(n_filled > 0)
   {
    triangulate_io.numberofregions = n_filled;
    triangulate_io.regionlist = 
     (double*) malloc(triangulate_io.numberofregions * 4 * sizeof(double));
    
    //Use the "filled" hole centre coordinate to define the region
    //And the attribute will be indexed from 2
    unsigned region_count = 1;
    //Loop over the filled regions
    for(std::set<unsigned>::iterator it = fill_index.begin();
        it!=fill_index.end();++it)
     {
      triangulate_io.regionlist[4*region_count-4] = 
       internal_polygon_pt[*it]->internal_point()[0];
      triangulate_io.regionlist[4*region_count-3] = 
       internal_polygon_pt[*it]->internal_point()[1];
      triangulate_io.regionlist[4*region_count-2] = 
       static_cast<double>(region_count);
      triangulate_io.regionlist[4*region_count-1] = 0.0;
      //Increase the number of regions
      ++region_count;
     }
   }
 
  // Storing the hole center coordinates
  triangulate_io.numberofholes = n_holes - n_filled;
  triangulate_io.holelist =
   (double*) malloc(triangulate_io.numberofholes * 2 * sizeof(double));

  for(unsigned count_hole=0;count_hole<n_holes*2;count_hole+=2)
   {
    //Only add the hole, if it is not filled in
    if(fill_index.find(count_hole)==fill_index.end())
     {
      triangulate_io.holelist[count_hole]=internal_polygon_pt[count_hole/2]->
       internal_point()[0];
      triangulate_io.holelist[count_hole+1]=internal_polygon_pt[count_hole/2]->
       internal_point()[1];
     }
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
                        "TriangleMesh<ELEMENT>::build_triangulateio()",
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
                        "TriangleMesh<ELEMENT>::create_triangle_input_data()",
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
  string test_string; 
 
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
                        "TriangleMesh<ELEMENT>::create_triangle_input_data",
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
 }

#endif

//======================================================================
/// Move the nodes on boundaries with associated Geometric Objects so 
/// that the exactly coincide with the geometric object. This requires
/// that the boundary coordinates are set up consistently
//======================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::snap_nodes_onto_geometric_objects()
 {
  //Loop over all boundaries
  const unsigned n_bound = this->nboundary();
  for(unsigned b=0;b<n_bound;b++)
   {
    //Find the geometric object
    GeomObject* const geom_object_pt = this->boundary_geom_object_pt(b);
    
    //If there is one
    if(geom_object_pt!=0)
     {
      Vector<double> b_coord(1);
      Vector<double> new_x(2);
      const unsigned n_boundary_node = this->nboundary_node(b);
      for(unsigned n=0;n<n_boundary_node;++n)
       {
        //Get the boundary node and coordinates
        Node* const nod_pt = this->boundary_node_pt(b,n);
        nod_pt->get_coordinates_on_boundary(b,b_coord);

        //Get the position and time history according to the underlying
        //geometric object, assuming that it has the same timestepper
        //as the nodes....
        unsigned n_tstorage = nod_pt->ntstorage();
        for(unsigned t=0;t<n_tstorage;++t)
         {
          //Get the position according to the underlying geometric object
          geom_object_pt->position(t,b_coord,new_x);

          //Move the node
          for(unsigned i=0;i<2;i++)
           {
            nod_pt->x(t,i) = new_x[i];
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
     "RefineableTriangleMesh::refine_triangulateio(",
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
 void RefineableTriangleMesh<ELEMENT>::adapt(OomphCommunicator* comm_pt,
                                             const Vector<double>& elem_error)
 {    
  // Get refinement targets
  Vector<double> target_area(elem_error.size());
  double min_angle=compute_area_target(elem_error,
                                       target_area);
  // Get maximum target area
  unsigned n=target_area.size();
  double max_area=0.0;
  double min_area=DBL_MAX;
  for (unsigned e=0;e<n;e++)
   {
    if (target_area[e]>max_area) max_area=target_area[e];
    if (target_area[e]<min_area) min_area=target_area[e];
   }
  
  oomph_info << "Maximum target area: " << max_area << std::endl;
  oomph_info << "Minimum target area: " << min_area << std::endl;
  oomph_info << "Number of elements to be refined " 
             << this->Nrefined << std::endl;
  oomph_info << "Number of elements to be unrefined "
             << this->Nunrefined << std::endl;
  oomph_info << "Min angle "<< min_angle << std::endl;

  double orig_max_area, orig_min_area;
  this->max_and_min_area(orig_max_area, orig_min_area);
  oomph_info << "Max/min area in original mesh: " 
             << orig_max_area  << " "
             << orig_min_area << std::endl;    

  // Should we bother to adapt?
  if ( (Nrefined > 0) || (Nunrefined > max_keep_unrefined()) ||
       (min_angle < min_permitted_angle()) )
   {

    if (! ( (Nrefined > 0) || (Nunrefined > max_keep_unrefined()) ) )
     {
      oomph_info 
       << "Mesh regeneration triggered by min angle criterion\n";
     }


    //Generate a new 1D mesh representation of the inner hole boundaries
    unsigned nhole=this->Internal_polygon_pt.size();
    Vector<Vector<double> > internal_point_coord(nhole);
    this->surface_remesh_for_inner_hole_boundaries(internal_point_coord);

    //Update the representation of the outer boundary
    this->surface_remesh_for_outer_boundary();

    //If there is not a geometric object associated with the boundary
    //the reset the boundary coordinates so that the lengths are consistent
    //in the new mesh and the old mesh.
    const  unsigned n_boundary = this->nboundary();
    for(unsigned b=0;b<n_boundary;++b)
     {
      if(this->boundary_geom_object_pt(b)==0)
       {
        this->setup_boundary_coordinates(b);
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
    Vector<TriangleMeshInternalClosedCurve*> hole_pt(nh);
    for (unsigned i=0;i<nh;i++)
     {
      hole_pt[i]=this->Internal_polygon_pt[i];
     }
    if (solid_mesh_pt!=0)
     {
      tmp_new_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>
       (closed_curve_pt,
        hole_pt,
        max_area,
        this->Time_stepper_pt,
        this->Fill_index,
        this->Use_attributes);
     }
    else
     {
      tmp_new_mesh_pt=new RefineableTriangleMesh<ELEMENT>
       (closed_curve_pt,
        hole_pt,
        max_area,
        this->Time_stepper_pt,
        this->Fill_index,
        this->Use_attributes);
     }

    // Snap to curvilinear boundaries (some code duplication as this
    // is repeated below but helper function would take so many
    // arguments that it's nearly as messy...
    
    //Pass the boundary  geometric objects to the new mesh
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
    //tmp_new_mesh_pt->output("pre_mesh_nodes_snapped_0.dat");
    
    //Move the nodes on the new boundary onto the 
    //old curvilinear boundary
    //If the boundary is straight this will do precisely nothing
    //but will be somewhat inefficient
    for(unsigned b=0;b<n_boundary;b++)
     {
      this->snap_nodes_onto_boundary(tmp_new_mesh_pt,b);
     }
    
    //Output the mesh after the snapping has taken place
    //tmp_new_mesh_pt->output("mesh_nodes_snapped_0.dat"); 
    
    // Get the TriangulateIO object associated with that mesh
    TriangulateIO tmp_new_triangulateio=
     tmp_new_mesh_pt->triangulateio_representation();
    
#ifdef PARANOID
    if (this->Problem_pt==0) 
     {
      throw OomphLibError("Problem pointer must be set with problem_pt()",
                          "TriangleMesh::adapt()",
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    RefineableTriangleMesh<ELEMENT>* new_mesh_pt=0;

    // Map storing target areas for elements in temporary 
    // TriangulateIO mesh
    std::map<GeneralisedElement*,double> target_area_map;


    //////////////////////////////////////////////////////////////
    // NOTE: Repeated setup of multidomain interaction could
    // be avoided by setting up a sufficiently fine bin
    // for the original mesh and reading out the target
    // area information from there
    //////////////////////////////////////////////////////////////

    // Now start iterating to refine mesh recursively
    //-----------------------------------------------
    bool done=false;
    unsigned iter=0;
    while (!done)
     {
      
      // "Project" target areas from current mesh onto uniform
      //------------------------------------------------------
      // background mesh
      //----------------
      
      // Temporarily switch on projection capabilities to allow
      // storage of pointer to external element.
      // Need to do this for both meshes to ensure that 
      // matching is done based on Eulerian coordinates for both
      // (in case we're dealing with solid meshes where the
      // locate_zeta would otherwise use the Lagrangian coordintes).
      unsigned nelem=this->nelement();
      for (unsigned e=0;e<nelem;e++)
       {
        dynamic_cast<ELEMENT*>(this->element_pt(e))->enable_projection();
       }
      unsigned nelem2=tmp_new_mesh_pt->nelement();
      for (unsigned e=0;e<nelem2;e++)
       {
        dynamic_cast<ELEMENT*>(tmp_new_mesh_pt->element_pt(e))->
         enable_projection();
       }

      // Set up multi domain interactions so we can figure out
      // which element in the intermediate uniform mesh is co-located
      // with given element in current mesh (which is to be refined)
      Multi_domain_functions::setup_multi_domain_interaction
       <ELEMENT>(this->Problem_pt,this,tmp_new_mesh_pt);
      
      target_area_map.clear();
      for (unsigned e=0;e<nelem;e++)
       {
        ELEMENT* el_pt=dynamic_cast<ELEMENT*>(this->element_pt(e));
        unsigned nint=el_pt->integral_pt()->nweight();
        for (unsigned ipt=0;ipt<nint;ipt++)
         {
          GeneralisedElement* ext_el_pt=el_pt->external_element_pt(0,ipt);

          // Use max. rather than min area of any element overlapping the
          // the current element, otherwise we get a rapid outward diffusion
          // of small elements
          target_area_map[ext_el_pt]=std::max(target_area_map[ext_el_pt],
                                              target_area[e]);
         }

        // Switch off projection capability          
        dynamic_cast<ELEMENT*>(this->element_pt(e))->disable_projection();
       }
      for (unsigned e=0;e<nelem2;e++)
       {
        dynamic_cast<ELEMENT*>(tmp_new_mesh_pt->element_pt(e))->
         disable_projection();
       }      

      // Now copy into target area for temporary mesh but limit to
      // the equivalent of one sub-division per iteration
      done=true;
      unsigned nel_new=tmp_new_mesh_pt->nelement();
      Vector<double> new_target_area(nel_new);
      for (unsigned e=0;e<nel_new;e++)
       {
        // No target area found for this element -- keep its size
        // by setting target area to -1 for triangle
        double new_area=target_area_map[tmp_new_mesh_pt->element_pt(e)];
        if (new_area<=0.0) 
         {
          new_target_area[e]=-1.0; 
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
       }
      
      

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
      //new_mesh_pt->output("pre_mesh_nodes_snapped_1.dat"); 
      
      //Move the nodes on the new boundary onto the 
      //old curvilinear boundary
      //If the boundary is straight this will do precisely nothing
      //but will be somewhat inefficient
      for(unsigned b=0;b<n_boundary;b++)
       {
        this->snap_nodes_onto_boundary(new_mesh_pt,b);
       }
      
      //Output the mesh after the snapping has taken place
      //new_mesh_pt->output("mesh_nodes_snapped_1.dat"); 

      
      // Not done: get ready for another iteration
      iter++;
      delete tmp_new_mesh_pt;
      if (!done)
       {
        tmp_new_mesh_pt=new_mesh_pt;
        tmp_new_triangulateio=new_mesh_pt->triangulateio_representation();
       }
      
     } // end of iteration
    
    // Project current solution onto new mesh
    //---------------------------------------
    ProjectionProblem<ELEMENT>* project_problem_pt=
     new ProjectionProblem<ELEMENT>;
    project_problem_pt->mesh_pt()=new_mesh_pt;
    project_problem_pt->project(this);
    
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
      this->Region_element_pt.resize(n_region);
      this->Region_attribute.resize(n_region);
      for(unsigned r=0;r<n_region;r++)
       {
        this->Region_attribute[r] = new_mesh_pt->region_attribute(r);
        //Find the number of elements in the region
        unsigned n_region_element = new_mesh_pt->nregion_element(r);
        this->Region_element_pt[r].resize(n_region_element);
        for(unsigned e=0;e<n_region_element;e++)
         {
          this->Region_element_pt[r][e] = new_mesh_pt->region_element_pt(r,e);
         }
       }

      //Now the boundary region information
      this->Boundary_region_element_pt.resize(nbound);
      this->Face_index_region_at_boundary.resize(nbound);
      
      //Now loop over the boundaries
      for(unsigned b=0;b<nbound;++b)
       {
        //Loop over the regions
        for(unsigned r=0;r<n_region;++r)
         {
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
               = new_mesh_pt->boundary_element_pt_in_region(b,r,e);
              this->Face_index_region_at_boundary[b][r][e] 
               = new_mesh_pt->face_index_at_boundary_in_region(b,r,e);
             }
           }
         }
       } //End of loop over boundaries

     } //End of case when more than one region

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

    // Solid mesh?
    if (solid_mesh_pt!=0)
     {
      // Warning
      std::stringstream error_message;
      error_message 
       << "Lagrangian coordinates are currently not projected but are\n"
       << "are re-set during adaptation. This is not appropriate for\n"
       << "real solid mechanics problems!\n";
      OomphLibWarning(error_message.str(),
                      "RefineableTriangleMesh::adapt()",
                      OOMPH_EXCEPTION_LOCATION);
      
      // Reset Lagrangian coordinates
      dynamic_cast<SolidMesh*>(this)->set_lagrangian_nodal_coordinates();
     }
    
    double max_area;
    double min_area;
    this->max_and_min_area(max_area, min_area);
    oomph_info << "Max/min area in adapted mesh: " 
               << max_area  << " "
               << min_area << std::endl;    
   }
  else
   {
    oomph_info << "Not enough benefit in adaptation.\n";
    Nrefined=0;
    Nunrefined=0;
   }
 }

 //======================================================================
 /// Helper function that updates the input polygon's PSLG
 /// by using the end-points of elements from FaceMesh(es) that are
 /// constructed for the boundaries associated with the segments of the
 /// polygon.
 //======================================================================
 template<class ELEMENT>
 void RefineableTriangleMesh<ELEMENT>::
 update_polygon_using_face_mesh(TriangleMeshPolygon* polygon_pt)
 {
  Vector<double> vertex_coord(3);
  Vector<double> bound_left(1);
  Vector<double> bound_right(1);
  
  //Loop over the number of polylines
  unsigned n_polyline = polygon_pt->npolyline();
     for(unsigned p=0;p<n_polyline;p++)
      {
       //Set of coordinates that will be placed on the boundary
       std::set<Vector<double> > vertex_nodes;
       
       //Get the boundary id of each polyline
       unsigned bound = 
        polygon_pt->polyline_pt(p)->boundary_id();
       
       // Create a face mesh adjacent to the fluid mesh's b-th boundary. 
       // The face mesh consists of FaceElements that may also be 
       // interpreted as GeomObjects
       Mesh* face_mesh_pt = new Mesh;
       this->template build_face_mesh
        <ELEMENT,FaceElementAsGeomObject>(bound, face_mesh_pt);
       
       // Loop over these new face elements and tell them the boundary number
       // from the bulk fluid mesh -- this is required to they can
       // get access to the boundary coordinates!
       unsigned n_face_element = face_mesh_pt->nelement();
       for(unsigned e=0;e<n_face_element;e++)
        {
         //Cast the element pointer to the correct thing!
         FaceElementAsGeomObject<ELEMENT>* el_pt=
          dynamic_cast<FaceElementAsGeomObject<ELEMENT>*>
          (face_mesh_pt->element_pt(e));
         
         // Set bulk boundary number
         el_pt->set_boundary_number_in_bulk_mesh(bound);
        }
       
       //Now we have the face mesh loop over the face elements and 
       //print out the end points
       for(unsigned e=0;e<n_face_element;++e)
        {
         unsigned n_node = face_mesh_pt->finite_element_pt(e)->nnode();
         FiniteElement* el_pt = face_mesh_pt->finite_element_pt(e);
         
         //Add the left-hand node to the list
         el_pt->node_pt(0)
          ->get_coordinates_on_boundary(bound,bound_left);
         vertex_coord[0] = bound_left[0];
         for(unsigned i=0;i<2;i++)
          {
           vertex_coord[i+1] = el_pt->node_pt(0)->x(i);
          }
         
         vertex_nodes.insert(vertex_coord);
         
         //Add the right-hand nodes to the list
         el_pt->node_pt(n_node-1)
          ->get_coordinates_on_boundary(bound,bound_right);
         vertex_coord[0] = bound_right[0];
         for(unsigned i=0;i<2;i++)
          {
           vertex_coord[i+1] = el_pt->node_pt(n_node-1)->x(i);
          }
         
         vertex_nodes.insert(vertex_coord);
         
         //Worry about bulk refinement here?
        }
        
       //Delete the allocated memory for theface mesh
       face_mesh_pt->flush_element_and_node_storage();
       delete face_mesh_pt;
       
       //Turn the set into a vector
       //Firstly trim any elements below a minimum size
       {
        double min_length = 0.01;
        std::set<Vector<double> > ::iterator it = vertex_nodes.begin();
        //Get the location of the "leftmost" (first) vertex
        Vector<double> left_vertex = *it;
        //Loop over all other vertices starting from the first+1
        //and stopping before the final vertex
        for(++it;it!=--vertex_nodes.end();++it)
         {
          //What is the actual length between the "left" vertex
          //and the current vertex
          double length = 0.0;
          for(unsigned i=0;i<2;i++)
           {
            length += pow(((*it)[i+1] - left_vertex[i+1]),2.0);
           }
          //If smaller than minimum delete the current vertex
          //Need to be careful when deleting entries in a set that
          //is being iterated over.
          if(sqrt(length) < min_length)
           {
            //Say so
            oomph_info << "Surface element too small: ";
            oomph_info << "Removing node at " << (*it)[1] << " " << (*it)[2]
                       << "\n";
            //Store the current value of the iterator
            std::set<Vector<double> >::iterator tmp_it = it;
            //Go back to the previous entry with the loop iterator
            --it;
            //Erase the offending entry
            vertex_nodes.erase(tmp_it);
           }
          //Otherwise, the "left" vertex becomes the current vertex
          else
           {
            left_vertex = *it;
           }
         }

        //If the final element is too small then remove the interior node
        //Get the last node
        it = --vertex_nodes.end();
        Vector<double> right_vertex = *it;
        //Now decrease iterator
        --it;
        
        //What is the actual length between the "right" vertex
        //and the current vertex
        double length = 0.0;
        for(unsigned i=0;i<2;i++)
         {
          length += pow(((*it)[i+1] - right_vertex[i+1]),2.0);
         }
        //If smaller than minimum delete the current vertex
        //Need to be careful when deleting entries in a set that
        //is being iterated over.
        if(sqrt(length) < min_length)
         {
          //Say so
          oomph_info << "Surface element too small: ";
          oomph_info << "Removing node at " << (*it)[1] << " " << (*it)[2]
                     << "\n";
          //Store the current value of the iterator
          std::set<Vector<double> >::iterator tmp_it = it;
          //Erase the offending entry
          vertex_nodes.erase(tmp_it);
         }
       }
        
       unsigned n_poly_vertex = vertex_nodes.size();
       Vector<Vector<double> > vector_vertex_node(n_poly_vertex);
       unsigned count=0;
       for(std::set<Vector<double> >::iterator it = vertex_nodes.begin();
           it!=vertex_nodes.end();++it)
        {
         vector_vertex_node[count].resize(2);
         vector_vertex_node[count][0] = (*it)[1];
         vector_vertex_node[count][1] = (*it)[2];
         ++count;
        }
       
       
       //Check whether the segments are continguous (first vertex of this
       //segment is equal to last vertex of previous segment).
       //If not, we should reverse the order of the current segment.
       //This check only applies for segments other than the first.
       if(p > 0)
        {
         //Final end point of previous line
         Vector<double> final_vertex_of_previous_segment;
         unsigned n_prev_vertex = 
          polygon_pt->polyline_pt(p-1)->nvertex();
         final_vertex_of_previous_segment = 
          polygon_pt->polyline_pt(p-1)->
          vertex_coordinate(n_prev_vertex-1);
         
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
             std::ostringstream error_stream;
             error_stream << 
              "The distance between the first node of the current\n" << 
              "line segment and either and of the previous line segment\n"
                          <<
              "is bigger than the desired tolerance " <<
              ToleranceForVertexMismatchInPolygons::Tolerable_error << ".\n"
                          << 
              "This suggests that the polylines defining the polygonal\n" <<
              "representation of the hole are not properly ordered.\n" <<
              "This should have failed when first trying to construct the\n"
                          << "polygon.\n";
             throw OomphLibError(error_stream.str(),
                                 "TriangleMesh::adapt()",
                                 OOMPH_EXCEPTION_LOCATION);
            }
           
           //Reverse the current vector to line up with the previous one
           std::reverse(vector_vertex_node.begin(),vector_vertex_node.end());
          }
        }
       
       //Now update the polyline according to the new vertices
       //Need offset because of stupid triangle
       delete polygon_pt->polyline_pt(p);
       polygon_pt->polyline_pt(p) = 
        new TriangleMeshPolyLine(vector_vertex_node,bound); //+1);
      }
    }



//======================================================================
/// Update the PSLG that define the inner boundaries of the mesh.
///
//======================================================================
 template <class ELEMENT>
 void RefineableTriangleMesh<ELEMENT>::
 surface_remesh_for_inner_hole_boundaries(Vector<Vector<double> >
                                          &internal_point_coord)
 {
  unsigned n_hole = internal_point_coord.size();
  for(unsigned ihole=0;ihole<n_hole;ihole++)
   {
    //Can I cast it
    RigidBodyTriangleMeshInternalPolygon* poly_pt
     = dynamic_cast<RigidBodyTriangleMeshInternalPolygon*>
     (this->Internal_polygon_pt[ihole]);

    //If it's a polygon, this is easy!
    if(poly_pt!=0)
     {
      poly_pt->reset_reference_configuration();
      
      // Initialize Vector hole_coordinates
      internal_point_coord[ihole].resize(2);
      
      // Get the vector of hole coordinates
      internal_point_coord[ihole]=
       this->Internal_polygon_pt[ihole]->internal_point();
     }
    //Otherwise we have to work much harder
    else
     {
      //Update the polygon associated with the ihole-th hole
      this->update_polygon_using_face_mesh(this->Internal_polygon_pt[ihole]);

      //Now sort out the hole coordinates
      Vector<double> vertex_coord;
      unsigned n_polyline = this->Internal_polygon_pt[ihole]->npolyline();
      
      vertex_coord.resize(2);
      // Initialize Vector hole_coordinates
      internal_point_coord[ihole].resize(2);
      
      //Hole centre will be found by averaging the position of 
      //all vertex nodes
      internal_point_coord[ihole][0] = 0.0;
      internal_point_coord[ihole][1] = 0.0;
      
      for(unsigned p=0;p<n_polyline;p++)
       {
        Vector<double> poly_ave(2,0.0);
        //How many vertices are there in the segment
        unsigned n_vertex =
         this->Internal_polygon_pt[ihole]->polyline_pt(p)->nvertex();
        for(unsigned v=0;v<n_vertex;v++)
         {
          vertex_coord = 
           this->Internal_polygon_pt[ihole]->polyline_pt(p)->
           vertex_coordinate(v);
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
      this->Internal_polygon_pt[ihole]->internal_point() =
       internal_point_coord[ihole];
     }
   }
 }


//======================================================================
/// Update the PSLG that define the outer boundary of the mesh.
///
//======================================================================
 template <class ELEMENT>
 void RefineableTriangleMesh<ELEMENT>::
 surface_remesh_for_outer_boundary()
 {
  //Only if there is a geometric object associated with the first boundary
  // (HACK)
  if(this->boundary_geom_object_pt(0)!=0)
   {
    //Update the polygon associated with the outer boundary
    this->update_polygon_using_face_mesh(this->Outer_boundary_pt);
   }
 }


//======================================================================
/// Move the boundary nodes onto the boundary defined by the old mesh
//======================================================================
 template <class ELEMENT>
 void RefineableTriangleMesh<ELEMENT>:: snap_nodes_onto_boundary(
  RefineableTriangleMesh<ELEMENT>* &new_mesh_pt, const unsigned &b)
 {

  // Quick return
  if (!Boundary_coordinate_exists[b])
   {
    oomph_info << "Not snapping nodes on boundary " << b 
               << " because no boundary coordinate has been set up.\n";
    return;
   }

  //Firstly we set the boundary coordinates of the new nodes
  //In case the mapping between the geometric object's intrinsic coordiante
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
    oomph_info << "Not snapping nodes on boundary " << b 
               << " because there aren't any.\n";
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
       (std::abs(old_boundary_node[m][1] - new_boundary_node[n][1]) < 1.0e-14)
       && 
       (std::abs(old_boundary_node[m][2] - new_boundary_node[n][2]) < 1.0e-14))
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
                 << 
     "on the boundary " << b 
                 << 
     ". This should not happen because these limits should have been setup\n"
                 <<
     "in the constructor\n";
    error_stream << 
     "The distance between the new and old nodes is probably outside\n"
                 << 
     "our tolerance.\n";
    error_stream.precision(20);
    error_stream << "Old boundaries: \n";
    error_stream << 
     old_boundary_node[0][1] << " " << old_boundary_node[0][2] 
                 << " : " <<
     old_boundary_node[n_boundary_node-1][1] << " "  <<
     old_boundary_node[n_boundary_node-1][2] << "\n";
    error_stream << "New boundaries: \n" <<
     new_boundary_node[0][1] << " " << new_boundary_node[0][2] << " : " <<
     new_boundary_node[n_new_boundary_node-1][1] << " "  <<
     new_boundary_node[n_new_boundary_node-1][2] << "\n";
    
    OomphLibWarning(error_stream.str(),
                    "RefineableTriangleMesh::snap_nodes_onto_boundary()",
                    OOMPH_EXCEPTION_LOCATION);
   }
#endif

  //The end points should always be present, so we
  //can (and must) always add them in exactly
  new_boundary_node[0][4] = new_boundary_node[0][0];
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
       "RefineableTriangleMesh::snap_nodes_onto_boundary()",
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

  //Now that the coordinates have been set up we can do the snapping

  // Create a face mesh adjacent to the fluid mesh's b-th boundary. 
  // The face mesh consists of FaceElements that may also be 
  // interpreted as GeomObjects
  Mesh* face_mesh_pt = new Mesh;
  this->template build_face_mesh
   <ELEMENT,FaceElementAsGeomObject>(b, face_mesh_pt);
    
  // Loop over these new face elements and tell them the boundary number
  // from the bulk fluid mesh -- this is required to they can
  // get access to the boundary coordinates!
  unsigned n_face_element = face_mesh_pt->nelement();
  for(unsigned e=0;e<n_face_element;e++)
   {
    //Cast the element pointer to the correct thing!
    FaceElementAsGeomObject<ELEMENT>* el_pt=
     dynamic_cast<FaceElementAsGeomObject<ELEMENT>*>
     (face_mesh_pt->element_pt(e));
    
    // Set bulk boundary number
    el_pt->set_boundary_number_in_bulk_mesh(b);
   }   
  
  //Now make the mesh as geometric object
  MeshAsGeomObject* mesh_geom_obj_pt = new MeshAsGeomObject(face_mesh_pt);
    
  //Now assign the new nodes positions based on the old meshes
  //potentially curvilinear boundary (its geom object incarnation)
  Vector<double> new_x(2);
  //for(unsigned n=0;n<n_new_boundary_node;n++)
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
  face_mesh_pt->flush_element_and_node_storage();
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

#ifdef PARANOID
    // Flag to indicate if we successully classified/dealt with the element
    bool success=false;
#endif

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
         "RefineableTriangleMesh::snap_nodes_onto_boundary()",
         OOMPH_EXCEPTION_LOCATION);
       }
      // If I get there I must not have thrown :)
      success=true;
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
           "RefineableTriangleMesh::snap_nodes_onto_boundary()",
           OOMPH_EXCEPTION_LOCATION);
         }
        else
         {
          throw OomphLibError(
           "Have a seven-noded element that's not a TElement<2,3>",
           "RefineableTriangleMesh::snap_nodes_onto_boundary()",
           OOMPH_EXCEPTION_LOCATION);
         }
       }
      // If I get there I must not have thrown :)
      success=true;
#endif  
      // Deal with six noded stuff for all (normal and enriched) elements 
      
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
           "RefineableTriangleMesh::snap_nodes_onto_boundary()",
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
                          "RefineableTriangleMesh::snap_nodes_onto_boundary()",
                          OOMPH_EXCEPTION_LOCATION);
     }
   }

  // Cleanup
  for (unsigned j=0;j<6;j++)
   {
    delete dummy_six_node_element.node_pt(j);
   }

 }

#endif
 
}
#endif
