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
#ifndef OOMPH_TRIANGLE_MESH_TEMPLATE_CC
#define OOMPH_TRIANGLE_MESH_TEMPLATE_CC


#include "triangle_mesh.template.h"
#include "../generic/map_matrix.h"


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
 Node_pt.resize(nnode_scaffold);

  // Set number of boundaries
 unsigned nbound=Tmp_mesh_pt->nboundary();

 //Resize the boundary information
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

 //Map of Element attribute pairs
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


 if (nnode_1d!=3)
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
 // Use efficient scheme
 else
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
       
       // By default, we assume the node is not new
       bool is_new=false;
       
       // Switch on local node number (all located on edges)
       switch (j)
        {
         
         // Node 3 is located between nodes 0 and 1
        case 3:
         
         edge_node1_pt=finite_element_pt(e)->node_pt(0);
         edge_node2_pt=finite_element_pt(e)->node_pt(1);
         if (central_edge_node_pt(edge_node1_pt,edge_node2_pt)==0)
          {
           is_new=true;
           central_edge_node_pt(edge_node1_pt,edge_node2_pt)=node_pt;
           central_edge_node_pt(edge_node2_pt,edge_node1_pt)=node_pt;
          }
         break;
         
         
         // Node 4 is located between nodes 1 and 2
        case 4:
         
         edge_node1_pt=finite_element_pt(e)->node_pt(1);
         edge_node2_pt=finite_element_pt(e)->node_pt(2);
         if (central_edge_node_pt(edge_node1_pt,edge_node2_pt)==0)
          {
           is_new=true;
           central_edge_node_pt(edge_node1_pt,edge_node2_pt)=node_pt;
           central_edge_node_pt(edge_node2_pt,edge_node1_pt)=node_pt;
          }
         break;
         
         
         // Node 6 is located between nodes 0 and 2
        case 5:
         
         edge_node1_pt=finite_element_pt(e)->node_pt(0);
         edge_node2_pt=finite_element_pt(e)->node_pt(2);
         if (central_edge_node_pt(edge_node1_pt,edge_node2_pt)==0)
          {
           is_new=true;
           central_edge_node_pt(edge_node1_pt,edge_node2_pt)=node_pt;
           central_edge_node_pt(edge_node2_pt,edge_node1_pt)=node_pt;
          }
         break;
         
        default:
         //Error
         throw OomphLibError("More than six nodes in TriangleMesh",
                             "TriangleMesh::build_from_scaffold()",
                             OOMPH_EXCEPTION_LOCATION);
        }
       
       if (is_new)
        {
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
                   << "with an output file specified as the second argument.\n"
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
     zeta[0]-=zeta_ref;
     if (zeta[0]<0.0) zeta[0]=-zeta[0];
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

}
#endif
