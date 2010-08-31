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
#include "../generic/multi_domain.h"
#include "../generic/projection.h"


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
         
        // Get lexicographically bottom left node but only 
        // use vertex nodes as candidates
        all_nodes_pt.insert(nod_pt);
        if (j==(nnod-1))
         {
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
/// Create TriangulateIO object from TriangleMeshPolygon and 
/// TriangleMeshHolePolygon
//======================================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::build_triangulateio(
  TriangleMeshPolygon* &outer_boundary_pt,
  Vector<TriangleMeshHolePolygon*> &inner_hole_pt,
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
  triangulate_io.numberofpoints=n_globalvertices;
  triangulate_io.numberofsegments=n_globalvertices;
      
  // Prepearing the TriangulateIO objects to store the values     
  triangulate_io.pointlist = 
   (double *) malloc(triangulate_io.numberofpoints * 2 * sizeof(double));
  triangulate_io.segmentlist = 
   (int *) malloc(triangulate_io.numberofsegments * 2 * sizeof(int));
  triangulate_io.segmentmarkerlist = 
   (int *) malloc(triangulate_io.numberofsegments * sizeof(int));

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
      // The count_sub_bound is used, instead of the idpolyline
      triangulate_io.segmentmarkerlist[count_tri/2]=count_sub_bound;
      
      // -1 because of the different enumeration between oomph_lib mesh
      // and the TriangulateIO structure!
            
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
        triangulate_io.pointlist[count_tri]= inner_hole_pt[count_hole]->
         polyline_pt(count_seg)->
         vertex_coordinate(count_vertices)[0];
        triangulate_io.pointlist[count_tri+1]= inner_hole_pt[count_hole]->
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
        // Check if the boundary id has been provided
        triangulate_io.segmentmarkerlist[count_tri/2]=count_sub_bound; 
        
        // -1 because of the different enumeration between oomph_lib mesh
        // and the TriangulateIO structure!
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
    error_stream  << "Error building TriangulateIO object\n"
                  << "Please, check TriangleMeshPolyLine and\n"
                  << "TriangleMeshPolygon provided"
                  <<std::endl;      
    throw OomphLibError(error_stream.str(),
                        "TriangleMeshBoundary::TriangleMeshBoundary()",
                        OOMPH_EXCEPTION_LOCATION);
   }
     
  // Storing the hole center coordinates
  triangulate_io.numberofholes = n_holes;
  triangulate_io.holelist =
   (double*) malloc(triangulate_io.numberofholes * 2 * sizeof(double));

  for(unsigned count_hole=0;count_hole<n_holes*2;count_hole+=2)
   {
    triangulate_io.holelist[count_hole] = inner_hole_pt[count_hole/2]->
     hole_coordinate()[0];
    triangulate_io.holelist[count_hole+1] = inner_hole_pt[count_hole/2]->
     hole_coordinate()[1];
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
 
  
  // Initialisation sub_bound counter
  unsigned count_sub_bound=0;
  
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


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

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
 }


//==============================================================
/// Write a Triangulateio_object file of the TriangulateIO object
/// String s is add to assign a different value for
/// input and/or output structure.
/// The function give the same result of the "report" function
/// included in the tricall.c, esternal_src.
//==============================================================
 template <class ELEMENT>
 void TriangleMesh<ELEMENT>::write_triangulateio(TriangulateIO& triangle, 
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


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

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


    // Update the holes' reference configuration (vertices and position
    // of the point in the hole)
    unsigned nhole=this->Inner_hole_pt.size();
    Vector<Vector<double> > hole_centre_coord(nhole);
    for(unsigned ihole=0;ihole<nhole;ihole++)
     {
      // Reset inner hole reference configuration
      this->Inner_hole_pt[ihole]->reset_reference_configuration();
      
      // Initialize Vector hole_coordinates
      hole_centre_coord[ihole].resize(2);
      
      // Get the vector of hole coordinates
      hole_centre_coord[ihole]=this->Inner_hole_pt[ihole]->hole_coordinate();
     }
    
    // Update the TriangulateIO structure according to the new nodes elastic 
    // displacement.
    this->update_triangulateio(hole_centre_coord);
    
    // Get the updated TriangulateIO object
    TriangulateIO tmp_triangulateio=this->triangulateio_representation();

    // Are we dealing with a solid mesh?
    SolidMesh* solid_mesh_pt=dynamic_cast<SolidMesh*>(this);

    // Build temporary uniform background mesh
    //----------------------------------------
    // with area set by maximum required area
    //---------------------------------------
    RefineableTriangleMesh<ELEMENT>* tmp_new_mesh_pt=0;
    if (solid_mesh_pt==0)
     {
      tmp_new_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>
       (this->Outer_boundary_pt,
        this->Inner_hole_pt,
        max_area,
        this->Time_stepper_pt);
     }
    else
     {
      tmp_new_mesh_pt=new RefineableTriangleMesh<ELEMENT>
       (this->Outer_boundary_pt,
        this->Inner_hole_pt,
        max_area,
        this->Time_stepper_pt);
     }

//     oomph_info << "Built background mesh with area: " << max_area << std::endl;
//     tmp_new_mesh_pt->output("background_mesh.dat");
//     this->output("actual_mesh_before_adapt.dat");

//     pause("done");

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
          this->Time_stepper_pt);
       }      
      // No solid mesh
      else
       { 
        new_mesh_pt=new RefineableTriangleMesh<ELEMENT>
         (new_target_area,
          tmp_new_triangulateio,
          this->Time_stepper_pt);
       }    
      
      
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

 
}
#endif
