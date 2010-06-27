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


#include "triangle_mesh.template.h"
#include "../generic/map_matrix.h"


#include <iostream>
using namespace std;

namespace oomph
{

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
  
  // Store the triangulateio node id 
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
      // Loop over new local nodes
      for(unsigned j=3;j<nnode;j++)
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
          throw OomphLibError("More than six nodes in TriangleMesh",
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
 
  // Loop over all elements on boundaries
  unsigned nel=this->nboundary_element(b);
  if (nel>0)
   {
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
    Node* bottom_left_node_pt=first_node_pt;

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
         
        // Get lexicographically bottom left node
        all_nodes_pt.insert(nod_pt);
        if (nod_pt->x(1)<bottom_left_node_pt->x(1))
         {
          bottom_left_node_pt=nod_pt;
         }
        else if (nod_pt->x(1)==bottom_left_node_pt->x(1))
         {
          if (nod_pt->x(0)<bottom_left_node_pt->x(0))            
           {
            bottom_left_node_pt=nod_pt;
           }
         }
       }
     }

    // Now adjust boundary coordinate so that the bottom left node
    // has a boundary coordinate of zero and that zeta increases
    // away from that point
    bottom_left_node_pt->get_coordinates_on_boundary(b,zeta);
    double zeta_ref=zeta[0];
    for (std::set<Node*>::iterator it=all_nodes_pt.begin();
         it!=all_nodes_pt.end();it++)
     {    
      Node* nod_pt=(*it);
      nod_pt->get_coordinates_on_boundary(b,zeta);

      // hierher Check way the self-test failed here!
      zeta[0]-=zeta_ref;
      if (zeta[0]<0.0)
       {
        zeta[0]=abs(zeta[0]);
       }
      nod_pt->set_coordinates_on_boundary(b,zeta);
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
//======================================================================
/// Clear the Triangulateio object
//======================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::clear_triangulateio()
 {
       
  // triangulateio object clear
  // Commentary starting with "/*" has been kept
  // from the original triangulate function 
  // (external_src/oomph_triangle/triangle.h or tricall.c)

  // Clear the point,attribute and marker list 
  free(Triangulateio.pointlist); 
  Triangulateio.numberofpoints = 0;
  free(Triangulateio.pointattributelist);
  free(Triangulateio.pointmarkerlist); 
  Triangulateio.numberofpointattributes = 0;

  // Clear the triangle, attribute,neighbor and area list 
  free(Triangulateio.trianglelist);    
  free(Triangulateio.triangleattributelist);
  free(Triangulateio.trianglearealist);
  free(Triangulateio.neighborlist);
  Triangulateio.numberoftriangles = 0; 
  Triangulateio.numberofcorners = 0;
  Triangulateio.numberoftriangleattributes = 0;

  // Clear the segment and marker list 
  free(Triangulateio.segmentlist);
  free(Triangulateio.segmentmarkerlist);

  // Clear edge list 
  free(Triangulateio.holelist);

  // Clear edge, maerke and norm list 
  free(Triangulateio.edgelist);
  free(Triangulateio.edgemarkerlist);
  free(Triangulateio.normlist);
  Triangulateio.numberofedges = 0;
  
  // Clear region list 
  free(Triangulateio.regionlist);
  Triangulateio.numberofedges = 0;
  Triangulateio.numberofregions = 0;

  // Do we need to clear the Triangulateio itself as well?
  // Delete the Triangulateio
  // free(Triangulateio);
 }

//======================================================================
/// Initialize the triagulateio structure
//======================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::initialize_triangulateio(
  triangulateio &triangle_out)
 {
       
  // triangulateio object initialization
  // Commentary starting with "/*" has been kept
  // from the original triangulate function 
  // (external_src/oomph_triangle/triangle.h or tricall.c)

  // Initialize the point list 
  /* Not needed if -N switch used */
  triangle_out.pointlist = (double *) NULL; 
  triangle_out.numberofpoints = 0;
  triangle_out.pointattributelist = (double *) NULL;  
  triangle_out.numberofpointattributes = 0;

  /* Not needed if -N or -B switch used */
  triangle_out.pointmarkerlist = (int *) NULL; 

  // Initialize the triangle list 
  /* Not needed if -E switch used */
  triangle_out.trianglelist = (int *) NULL;    
  
  /* Not needed if -E switch used or number of triangle attributes is zero:*/
  triangle_out.triangleattributelist = (double *) NULL;
  triangle_out.trianglearealist = (double *) NULL;
  triangle_out.numberoftriangles = 0; 
  triangle_out.numberofcorners = 0;
  triangle_out.numberoftriangleattributes = 0;

  // Initialize the segment list 
  /* Needed only if segments are output (-p or -c) and -P not used:*/
  triangle_out.segmentlist = (int *) NULL;
  
  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  triangle_out.segmentmarkerlist = (int *) NULL;
  
  // Initialize edge list 
  triangle_out.edgelist = (int *) NULL;
  triangle_out.edgemarkerlist = (int *) NULL;
  triangle_out.numberofedges = 0;
  
  // Initialize region list 
  triangle_out.regionlist = (double *) NULL;
  triangle_out.numberofedges = 0;
  triangle_out.numberofregions = 0;
 }

//======================================================================
/// Build TriangleMeshPolygon and TriangleMeshHolePolygon
/// to triangulateio object
//======================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::build_triangulateio(
  TriangleMeshPolygon* &outer_boundary_pt,
  Vector<TriangleMeshHolePolygon*> &inner_hole_pt,
  triangulateio &triangle_in)
 {

  // triangle_in initialization
  this->initialize_triangulateio(triangle_in);
   
  // Build the triangulateio in object

  // Get number of polyline
  unsigned n_boundline=0;
  n_boundline = outer_boundary_pt->npolyline();
  
  // Count the global polygon boundary polyline
  unsigned n_boundglobalseg=0;

  // Count the global holes noundary polyline
  unsigned n_holes=0;
  n_holes = inner_hole_pt.size();
  
  // Initialize global hole boundary counter
  unsigned multi_hole_polyline=0;

  // Resize the inner box sub_boundary
  for(unsigned count_hole=0;count_hole<n_holes;count_hole++)
   {
    unsigned inner_hole_polyline=inner_hole_pt[count_hole]->npolyline();
    
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
    n_holepolyline=inner_hole_pt[count_hole]->npolyline();
    for(unsigned count_seg=0;count_seg<n_holepolyline;count_seg++)
     {
      n_globalvertices += inner_hole_pt[count_hole]->
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
    n_holepolyline=inner_hole_pt[count_hole]->npolyline();
    if(n_holepolyline==1)
     {
      n_globalvertices +=1;
     }
   }

  // Store the global number of vertices and segments
  // in the list. They are supposed to be the same!!
  triangle_in.numberofpoints=n_globalvertices;
  triangle_in.numberofsegments=n_globalvertices;
      
  // Prepearing the triangulateio objects to store the values     
  triangle_in.pointlist = 
   (double *) malloc(triangle_in.numberofpoints * 2 * sizeof(double));
  triangle_in.segmentlist = 
   (int *) malloc(triangle_in.numberofsegments * 2 * sizeof(int));
  triangle_in.segmentmarkerlist = 
   (int *) malloc(triangle_in.numberofsegments * sizeof(int));

  // Initialisation sub_bound counter
  unsigned count_sub_bound=1;     
  
  // Storing all the values in the list
  for(unsigned count_seg=0;count_seg<n_boundline;count_seg++)
   {

    // Storing the number of the vertices
    n_polylinevertices = outer_boundary_pt->polyline_pt(count_seg)->nvertex()-1;

    // If there's just one boundary. All the vertices have to be counted   
    if(n_boundline==1)
     {
      n_polylinevertices +=1;
     }
    
    // Store the segmen_id if given by the user
    unsigned idpolyline=outer_boundary_pt->polyline_pt(count_seg)
     ->boundary_id();
    
    // Initialize to zero the sub_bound counter
    unsigned count_bound_segment=0;

    // Initialize the vector of sub boundary
    Vector<unsigned>sub_bound_id(n_polylinevertices);
    
    //Storing the coordinates for each points
    for(unsigned count_vertices=0;count_vertices<n_polylinevertices;
        count_vertices++)
     {

      triangle_in.pointlist[count_tri]= outer_boundary_pt->polyline_pt
       (count_seg)
       ->vertex_coordinate(count_vertices)[0];
      triangle_in.pointlist[count_tri+1]= outer_boundary_pt->
       polyline_pt(count_seg)->vertex_coordinate(count_vertices)[1];

      // Store the segment values
      // If the segment is not the last one, take the next node
      if(count_seg==(n_boundline-1) && 
         count_vertices==(n_polylinevertices-1))
       {
        triangle_in.segmentlist[count_tri]=edge_segment;
        triangle_in.segmentlist[count_tri+1]=1;
       }
      else
       {
        triangle_in.segmentlist[count_tri]=edge_segment;
        triangle_in.segmentlist[count_tri+1]=edge_segment+1;
        edge_segment++;
       }

      // Store the marker list of the segments
      // The count_sub_bound is used, instead of the idpolyline
      triangle_in.segmentmarkerlist[count_tri/2]=count_sub_bound;
      
      // -1 because of the different enumeration between oomph_lib mesh
      // and the triangulateio structure!
            
      // Build the vector os sub boundary id for the boundary "idpolyline"
      sub_bound_id[count_bound_segment]=count_sub_bound-1;

      // Increment counter
      count_sub_bound++;
      count_bound_segment++;
      count_tri+=2;
      n_boundglobalseg++;
     }
    
    // Add sub boundary vector to the map
    Sub_boundary_id[idpolyline]=sub_bound_id;
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
    n_holepolyline=inner_hole_pt[count_hole]->npolyline();
    
    for(unsigned count_seg=0;count_seg<n_holepolyline;count_seg++)
     {
      // Storing the number of the vertices
      n_polylinevertices = inner_hole_pt[count_hole]->polyline_pt(count_seg)->
       nvertex()-1;
      
      // If there's just one boundary. All the vertices should be counted   
      if(n_holepolyline==1)
       {
        n_polylinevertices +=1;
       }

      // Store the segment_id if given by the user
      unsigned idpolyline=inner_hole_pt[count_hole]->
       polyline_pt(count_seg)->boundary_id();
      
      // Initialize to zero the sub_bound counter
      unsigned count_bound_segment=0;
      
      // Initialize the vector of sub boundary
      Vector<unsigned>sub_bound_id(n_polylinevertices);

      // Store the coordinates for each points
      for(unsigned count_vertices=0;count_vertices<n_polylinevertices;
          count_vertices++)
       {
        triangle_in.pointlist[count_tri]= inner_hole_pt[count_hole]->
         polyline_pt(count_seg)->
         vertex_coordinate(count_vertices)[0];
        triangle_in.pointlist[count_tri+1]= inner_hole_pt[count_hole]->
         polyline_pt(count_seg)->
         vertex_coordinate(count_vertices)[1];

        // Store the segments values
        // If the segment is not the last one, take the next node
        if(count_seg==(n_holepolyline-1) && 
           count_vertices==(n_polylinevertices-1))
         {
          triangle_in.segmentlist[count_tri]=edge_segment;
          triangle_in.segmentlist[count_tri+1]=hole_vertex_start;
         }
        else
         {
          triangle_in.segmentlist[count_tri]=edge_segment;
          triangle_in.segmentlist[count_tri+1]=edge_segment+1;
          edge_segment++;
         }

        // Store the marker list of the segments
        // Check if the boundary id has been provided
        triangle_in.segmentmarkerlist[count_tri/2]=count_sub_bound; 
        
        // -1 because of the different enumeration between oomph_lib mesh
        // and the triangulateio structure!
        // Build the vector of sub boundary id for the boundary "idpolyline"
        sub_bound_id[count_bound_segment]=count_sub_bound-1;
        
        // Increment counter
        count_bound_segment++;
        count_sub_bound++;
        count_tri+=2;
        n_boundglobalseg++;
       }
      
      // Add sub boundary vector to the map
      Sub_boundary_id[idpolyline]=sub_bound_id;
     }
    // Add the previous number of boundary
    nmulti_hole_bound+=n_holepolyline;
   }

  // Check if the number of segments and vertices are the the same
  if(n_boundglobalseg!=n_globalvertices)
   {
    std::ostringstream error_stream;
    error_stream  << "Error building triangulateio object\n"
                  << "Please, check TriangleMeshPolyLine and\n"
                  << "TriangleMeshPolygon provided"
                  <<std::endl;      
    throw OomphLibError(error_stream.str(),
                        "TriangleMeshBoundary::TriangleMeshBoundary()",
                        OOMPH_EXCEPTION_LOCATION);
   }
     
  // Storing the hole center coordinates
  triangle_in.numberofholes = n_holes;
  triangle_in.holelist =
   (double*) malloc(triangle_in.numberofholes * 2 * sizeof(double));

  for(unsigned count_hole=0;count_hole<n_holes*2;count_hole+=2)
   {
    triangle_in.holelist[count_hole] = inner_hole_pt[count_hole/2]->
     hole_coordinate()[0];
    triangle_in.holelist[count_hole+1] = inner_hole_pt[count_hole/2]->
     hole_coordinate()[1];
   }

  // triangulateio in object built
 }

//========================================================================
/// Create triangulateio object via the .poly file
//========================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::build_triangulateio(
  const std::string& poly_file_name,
  triangulateio &triangle_data)
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
  initialize_triangulateio(triangle_data);
  
  // Ignore the first line with structure description
  poly_file.ignore(80,'\n');
 
  // Read and store number of nodes
  unsigned invertices;
  poly_file>>invertices;
  triangle_data.numberofpoints=invertices; 

  // Initialisation of the point list
  triangle_data.pointlist = 
   (double *) malloc(triangle_data.numberofpoints * 2 * sizeof(double));

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

  triangle_data.numberofpointattributes = 0;
  triangle_data.pointattributelist = (double *) NULL;
   
  // Read and check the flag for boundary markers
  unsigned nodemarkers;
  poly_file>>nodemarkers;
  triangle_data.pointmarkerlist = (int *) NULL;

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
    poly_file>>triangle_data.pointlist[count_point];
    poly_file>>triangle_data.pointlist[count_point+1];
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

  triangle_data.numberofsegments = inelements;
  triangle_data.segmentlist = 
   (int *) malloc(triangle_data.numberofsegments * 2 * sizeof(int));
  triangle_data.segmentmarkerlist = 
   (int *) malloc(triangle_data.numberofsegments * sizeof(int));
 
  
  // Initialisation sub_bound counter
  unsigned count_sub_bound=0;
  
  // Read all the segments edges and markers
  for(unsigned i=0;i<2*inelements;i+=2)
   {
    poly_file>>dummy_seg;
    poly_file>>triangle_data.segmentlist[i];
    poly_file>>triangle_data.segmentlist[i+1];
    if(segment_markers!=0)
     {
      poly_file>>triangle_data.segmentmarkerlist[i/2];
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
 
    triangle_data.numberofholes = nhole;
    triangle_data.holelist = 
     (double *) malloc(triangle_data.numberofholes * 2 * sizeof(double));

    // Loop over the holes to get centre coords and store value onto the 
    // triangulateio object
    for(unsigned i=0;i<2*nhole;i+=2)
     {
      poly_file>>dummy_hole;
      poly_file>>triangle_data.holelist[i];
      poly_file>>triangle_data.holelist[i+1];
     }
   }  
 }

//========================================================================
/// Build a new triangulateio object, copying the previous triangulateio
/// and updating the maximum area for each element, driven by the
/// estimate computed
//========================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::refine_triangulateio(
  struct triangulateio &triangle_in,
  Vector<double> &error_elem,
  const double &error_target,
  struct triangulateio &triangle_refine)
 {
  
  //  Initialize triangulateio structure
  this->initialize_triangulateio(triangle_refine);

  // Store the global number of vertices and segments
  // in the list  
  unsigned n_points = triangle_in.numberofpoints;
  triangle_refine.numberofpoints=n_points;
  
  unsigned n_segments=triangle_in.numberofsegments;
  triangle_refine.numberofsegments=n_segments;
      
  // Initialization of the triangulateio objects to store the values     
  triangle_refine.pointlist = 
   (double *) malloc(triangle_in.numberofpoints * 2 * sizeof(double));
  triangle_refine.pointmarkerlist = 
   (int *) malloc(triangle_in.numberofpoints * sizeof(int));
  triangle_refine.segmentlist = 
   (int *) malloc(triangle_in.numberofsegments * 2 * sizeof(int));
  triangle_refine.segmentmarkerlist = 
   (int *) malloc(triangle_in.numberofsegments * sizeof(int));
     
  // Storing the point's coordinates in the list
  // and in two vectors with x and y coordinates
  Vector<double> x_coord (n_points);
  Vector<double> y_coord (n_points);
  
  for(unsigned count_point=0;count_point<n_points*2;count_point++)
   {
    triangle_refine.pointlist[count_point]=triangle_in.pointlist[count_point];
    
    // Even vaules represent the x coordinate
    // Odd values represent the y coordinate
    if (count_point%2==0)
     {
      x_coord[count_point*0.5] = triangle_in.pointlist[count_point];
     }
    else
     {
      y_coord[(count_point-1)*0.5] = triangle_in.pointlist[count_point];
     }
   }

  // Store the point's markers in the list
  for(unsigned count_marker=0;count_marker<n_points;count_marker++)
   {
    triangle_refine.pointmarkerlist[count_marker]=
     triangle_in.pointmarkerlist[count_marker];
   }

  // Storing the segment's edges in the list
  for(unsigned count_seg=0;count_seg<n_segments*2;count_seg++)
   {
    triangle_refine.segmentlist[count_seg]=triangle_in.segmentlist[count_seg];
   }

  // Store the segment's markers in the list
  for(unsigned count_markers=0;count_markers<n_segments;count_markers++)
   {
    triangle_refine.segmentmarkerlist[count_markers]=
     triangle_in.segmentmarkerlist[count_markers];
   }

  // Store the hole's center coordinates
  unsigned n_holes = triangle_in.numberofholes;
  triangle_refine.numberofholes = n_holes;

  triangle_refine.holelist =
   (double*) malloc(triangle_in.numberofholes * 2 * sizeof(double));

  // Loop over the holes to get centre coords  
  for(unsigned count_hole=0;count_hole<n_holes*2;count_hole++)
   {
    triangle_refine.holelist[count_hole] = triangle_in.holelist[count_hole];
   }
  
  // Store the triangles values  
  unsigned n_triangles = triangle_in.numberoftriangles;
  triangle_refine.numberoftriangles = n_triangles;
  
  unsigned n_corners = triangle_in.numberofcorners;
  triangle_refine.numberofcorners = n_corners;
  
  triangle_refine.trianglelist = 
   (int *) malloc(triangle_in.numberoftriangles * 3 * sizeof(int)); 
 
  // Store the triangle's corners in the list
  // Fill in the vector of the element's area list
  Vector<double> elem_area(n_triangles);
   
  for(unsigned count_tri=0;count_tri<n_triangles*3;count_tri++)
   {
    triangle_refine.trianglelist[count_tri]=
     triangle_in.trianglelist[count_tri];

    // In order to refine the meshing we need to know the area value 
    // of each element (triangle). The area value is obtained by the formula:
    // area = |(Ax(By-Cy)+Bx(Cy-Ay)+Cx(Ay-By))/2|
    
    // We compute one element at a time, hence every 3 list's values
    if (count_tri%3==0)
     {
      unsigned count_tri_a = triangle_in.trianglelist[count_tri];
      double ax= x_coord[count_tri_a-1];
      double ay= y_coord[count_tri_a-1];
      
      unsigned count_tri_b = triangle_in.trianglelist[count_tri+1];
      double bx= x_coord[count_tri_b-1];
      double by= y_coord[count_tri_b-1];
      
      unsigned count_tri_c = triangle_in.trianglelist[count_tri+2];
      double cx= x_coord[count_tri_c-1];
      double cy= y_coord[count_tri_c-1];

      // Area value
      elem_area[count_tri/3] = (ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))*0.5;
      
     }
   }

  // Store the triangle's area in the list
  // This list give the refine criteria. Area parameter
  // is the maximum area of the index triangle 
  triangle_refine.trianglearealist = 
   (double *) malloc(triangle_in.numberoftriangles * sizeof(double)); 
 
  // The new area value is computed according to
  // the error estimate of each element
  double area = 0.0;
  double no_refine = -1;

  for(unsigned count_area=0;count_area<n_triangles;count_area++)
   {
    // Initialize error ratio
    double error_ratio=0;
    
    // Find the best element's area value  
    // according to its error estimate
    if(error_elem[count_area]!=0.0)
     {
      error_ratio=error_target/error_elem[count_area];
     }
    else
     {
      error_ratio=Max_error_ratio;
     }
    
    // Condition to save a segmation fault for error_ratio
    // values over the double range
    if(error_ratio>=Min_error_ratio && error_ratio<=Max_error_ratio )
     {
      // With unrefinement        
      area = abs(elem_area[count_area]*error_ratio); 
      triangle_refine.trianglearealist[count_area]=area;
     }
    else
     {
      // Triangulateio function doesn't refine if the 
      // are value is -1
      triangle_refine.trianglearealist[count_area]=no_refine;
     }
   }
  
  // triangulateio in object built
 }
//==============================================================
/// Write a Triangulateio_object file of the triangulateio object
/// String s is add to assign a different value for
/// input and/or output structure.
/// The function give the same result of the "report" function
/// included in the tricall.c, esternal_src.
//==============================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::write_triangulateio(
  struct triangulateio &triangle, 
  std::string &s)
 {

  std::ofstream outfile;
  char filename[100];

  sprintf(filename,"Triangulateio_object_%s.dat",s.c_str());
  outfile.open(filename);
  outfile <<"# Triangulateio object values:\n\n"<<std::endl;

  // Write points coordinates
  if(triangle.numberofpoints!=0)
   {
    outfile <<"# Triangulateio number of points is:"
            <<triangle.numberofpoints <<std::endl;
   }
  if(triangle.pointlist != NULL)
   {
    outfile <<"# Vertex coordinates are:"<<std::endl;
    for(int k=0;k<triangle.numberofpoints*2;k+=2)
     {
      outfile << (k*0.5)+1 << " "
              << triangle.pointlist[k] << " "
              << triangle.pointlist[k+1] << std::endl;
     }
   }

  // Write points attribute list
  if(triangle.numberofpointattributes!=0)
   {
    outfile <<"# Triangulateio number of points attributelist is:"
            <<triangle.numberofpointattributes <<std::endl;
   }
  if(triangle.pointattributelist != NULL)
   {
       
    outfile <<"# Vertex attribute are:"<<std::endl;
    for(int k=0;k<triangle.numberofpointattributes;k++)
     {
      outfile << triangle.pointattributelist[k] << std::endl;
     }
   }

  // Write point markers list
  if(triangle.pointmarkerlist != NULL)
   {
    outfile <<"# Vertex Markers are:"<<std::endl;
    for(int k=0;k<triangle.numberofpoints;k++)
     {
      outfile << triangle.pointmarkerlist[k] << std::endl;
     }
   }
  
  // Write the 1.node file used by the showme function  
  std::ofstream nodefile;
  char nodename[100];

  sprintf(nodename,"file_%s.1.node",s.c_str());
  nodefile.open(nodename);
  nodefile <<triangle.numberofpoints<<" 2 "
           <<triangle.numberofpointattributes
           <<" 0"<<std::endl;
  for(int j=0;j<triangle.numberofpoints*2;j+=2)
   {
    nodefile <<(j/2)+1<<" " << triangle.pointlist[j] << " "
            << triangle.pointlist[j+1] << std::endl;
   }
  nodefile.close();
  

  // Write segments edge elements
  if(triangle.numberofsegments!=0)
   {
    outfile <<"# The number of segments is:"
            <<triangle.numberofsegments<<std::endl;
   }
  if(triangle.segmentlist != NULL)
   {
    outfile <<"# Segments are:"<<std::endl;
    for(int k=0;k<triangle.numberofsegments*2;k+=2)
     {
      outfile << triangle.segmentlist[k] << "  "
              << triangle.segmentlist[k+1] <<std::endl;
     }
   }

  // Write segments markers list
  if(triangle.segmentmarkerlist != NULL)
   {
    outfile <<"# Segments Markers are:"<<std::endl;
    for(int k=0;k<triangle.numberofsegments;k++)
     {
      outfile << triangle.segmentmarkerlist[k] << std::endl;
     }
   }

  // Write regions
  if(triangle.numberofregions!=0)
   {
    outfile <<"# The number of region is:"
            <<triangle.numberofregions<<std::endl;
   }

  // Write holes
  if(triangle.numberofholes!=0)
   {
    outfile <<"# The number of holes is:"
            <<triangle.numberofholes<<std::endl;
   }
  if(triangle.holelist != NULL)
   {
    outfile <<"#  Holes are:"<<std::endl;
    for(int k=0;k<triangle.numberofholes*2;k+=2)
     {
      outfile << triangle.holelist[k] << "  "
              << triangle.holelist[k+1] <<std::endl;
     }
   }

  // Write triangles
  if(triangle.numberoftriangles!=0)
   {
    outfile <<"# Triangulateio number of triangles:"
            <<triangle.numberoftriangles <<std::endl;
   }
  if(triangle.numberofcorners!=0)
   {
    outfile <<"# Triangulateio number of corners:"
            <<triangle.numberofcorners<<std::endl;
   }
  if(triangle.numberoftriangleattributes!=0)
   {
    outfile <<"# Triangulateio number of triangles attributes:"
            <<triangle.numberoftriangleattributes <<std::endl;
   }
  if(triangle.trianglelist != NULL)
   {
       
    outfile <<"# Traingles are:"<<std::endl;
    for(int k=0;k<triangle.numberoftriangles*3;k+=3)
     {
      outfile << triangle.trianglelist[k] << " "
              << triangle.trianglelist[k+1] << " "
              << triangle.trianglelist[k+2] << std::endl;
     }
   }

  if(triangle.trianglearealist != NULL)
   {
    outfile <<"# Triangle's areas are:"<<std::endl;
     for(int k=0;k<triangle.numberoftriangles;k++)
     {
      outfile << triangle.trianglearealist[k] << std::endl;
     }
   }

  if(triangle.trianglelist != NULL)
   {
    
    // Write the 1.ele file used by the showme function  
    std::ofstream elefile;
    char elename[100];
  
    sprintf(elename,"file_%s.1.ele",s.c_str());
    elefile.open(elename);
    elefile <<triangle.numberoftriangles<<" 3 0"<<std::endl;
    for(int j=0;j<triangle.numberoftriangles*3;j+=3)
     {
      elefile <<(j/3)+1<<" "<< triangle.trianglelist[j] << " "
              << triangle.trianglelist[j+1] << " "
              << triangle.trianglelist[j+2] << std::endl;
     }
    elefile.close();
   }
  
  outfile.close();
 }

}
#endif
