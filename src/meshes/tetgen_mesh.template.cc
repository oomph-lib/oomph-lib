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
#ifndef OOMPH_TETGEN_MESH_TEMPLATE_CC
#define OOMPH_TETGEN_MESH_TEMPLATE_CC


#include<algorithm>

#include "tetgen_mesh.template.h"
#include "../generic/map_matrix.h"



namespace oomph
{

//========================================================================
/// Build unstructured tet mesh based on output from scaffold
//========================================================================
template <class ELEMENT>
void TetgenMesh<ELEMENT>::build_from_scaffold(TimeStepper* time_stepper_pt)
{   
 // Create space for elements
 unsigned nelem=Tmp_mesh_pt->nelement();
 Element_pt.resize(nelem);
   
 // Create space for nodes
 unsigned nnode_scaffold=Tmp_mesh_pt->nnode();
 Node_pt.resize(nnode_scaffold);
   
 // Set number of boundaries
 unsigned nbound=Tmp_mesh_pt->nboundary();
 set_nboundary(nbound);
   
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
       Node_pt[global_count-1]->x(0)=scaffold_node_pt->x(0);
       Node_pt[global_count-1]->x(1)=scaffold_node_pt->x(1);
       Node_pt[global_count-1]->x(2)=scaffold_node_pt->x(2);
         
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
  }

 // At this point we've created all the elements and 
 // created their vertex nodes. Now we need to create
 // the additional (midside and internal) nodes!


 // We'll first create all local nodes for all elements
 // and then delete the superfluous ones that have
 // a matching node in an adjacent element.

 unsigned boundary_id;   

 // Get number of nodes along element edge and dimension of element (3)
 // from first element
 unsigned nnode_1d=finite_element_pt(0)->nnode_1d();

 // At the moment we're only able to deal with nnode_1d=2 or 3.
 if ((nnode_1d!=2)&&(nnode_1d!=3))
  {
   std::ostringstream error_message;
   error_message << "Mesh generation from tetgen currently only works\n";
   error_message << "for nnode_1d = 2 or 3. You're trying to use it\n";
   error_message << "for nnode_1d=" << nnode_1d << std::endl;
   
   throw OomphLibError(error_message.str(),
                       "TetgenMesh::build_from_scaffold()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Spatial dimension of element = number of local coordinates
 unsigned dim=finite_element_pt(0)->dim();
  
 // Storage for the local coordinate of the new node
 Vector<double> s(dim);

 // Get number of nodes in the element from first element
 unsigned nnode=finite_element_pt(0)->nnode();

 // Loop over all elements
 for (unsigned e=0;e<nelem;e++)
  {
   // Loop over the new nodes in the element and create them.
   for(unsigned j=4;j<nnode;j++)
    {

     // Create new node
     Node* new_node_pt=finite_element_pt(e)->
      construct_node(j,time_stepper_pt);

     // What are the node's local coordinates?
     finite_element_pt(e)->local_coordinate_of_node(j,s);

     // Find the coordinates of the new node from the existing
     // and fully-functional element in the scaffold mesh
     for(unsigned i=0;i<dim;i++)
      {
       new_node_pt->x(i)=
        Tmp_mesh_pt->finite_element_pt(e)->interpolated_x(s,i);
      }
   
    } // end of loop over new nodes
  } //end of loop over elements
   
 
 // Map storing the mid-side of an edge; edge identified by
 // pointers to vertex nodes
 MapMatrix<Node*,Node*> central_edge_node_pt;
 Node* edge_node1_pt=0;
 Node* edge_node2_pt=0;

 // Loop over elements
 for (unsigned e=0;e<nelem;e++)
  {
   // Loop over new local nodes
   for(unsigned j=4;j<nnode;j++)
    {
     // Pointer to the element's local node
     Node* node_pt=finite_element_pt(e)->node_pt(j);
    
     // By default, we assume the node is not new
     bool is_new=false;
     
     // This will have to be changed for higher-order elements
     //=======================================================
     
     // Switch on local node number (all located on edges)
     switch (j)
      {
       
       // Node 4 is located between nodes 0 and 1
      case 4:
       
       edge_node1_pt=finite_element_pt(e)->node_pt(0);
       edge_node2_pt=finite_element_pt(e)->node_pt(1);
       if (central_edge_node_pt(edge_node1_pt,edge_node2_pt)==0)
        {
         is_new=true;
         central_edge_node_pt(edge_node1_pt,edge_node2_pt)=node_pt;
         central_edge_node_pt(edge_node2_pt,edge_node1_pt)=node_pt;
        }
       break;
           
           
       // Node 5 is located between nodes 0 and 2
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
           

       // Node 6 is located between nodes 0 and 3
      case 6:

       edge_node1_pt=finite_element_pt(e)->node_pt(0);
       edge_node2_pt=finite_element_pt(e)->node_pt(3);
       if (central_edge_node_pt(edge_node1_pt,edge_node2_pt)==0)
        {
         is_new=true;
         central_edge_node_pt(edge_node1_pt,edge_node2_pt)=node_pt;
         central_edge_node_pt(edge_node2_pt,edge_node1_pt)=node_pt;
        }
       break;

       
       // Node 7 is located between nodes 1 and 2
      case 7:
       
       edge_node1_pt=finite_element_pt(e)->node_pt(1);
       edge_node2_pt=finite_element_pt(e)->node_pt(2);
       if (central_edge_node_pt(edge_node1_pt,edge_node2_pt)==0)
        {
         is_new=true;
         central_edge_node_pt(edge_node1_pt,edge_node2_pt)=node_pt;
         central_edge_node_pt(edge_node2_pt,edge_node1_pt)=node_pt;
        }
       break;


       // Node 8 is located between nodes 2 and 3
      case 8:
      
       edge_node1_pt=finite_element_pt(e)->node_pt(2);
       edge_node2_pt=finite_element_pt(e)->node_pt(3);
       if (central_edge_node_pt(edge_node1_pt,edge_node2_pt)==0)
        {
         is_new=true;
         central_edge_node_pt(edge_node1_pt,edge_node2_pt)=node_pt;
         central_edge_node_pt(edge_node2_pt,edge_node1_pt)=node_pt;
        }
       break;
     
           
       // Node 9 is located between nodes 1 and 3
      case 9:
           
       edge_node1_pt=finite_element_pt(e)->node_pt(1);
       edge_node2_pt=finite_element_pt(e)->node_pt(3);
       if (central_edge_node_pt(edge_node1_pt,edge_node2_pt)==0)
        {
         is_new=true;
         central_edge_node_pt(edge_node1_pt,edge_node2_pt)=node_pt;
         central_edge_node_pt(edge_node2_pt,edge_node1_pt)=node_pt;
        }
       break;
       
      default:
       //Error
       throw OomphLibError("More than ten nodes in Tet Element",
                           "TetgenMesh::build_from_scaffold()",
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

 //Boundary conditions
  
 // Matrix map to check if a node has already been added to 
 // the boundary number b 
 MapMatrixMixed<Node*,unsigned,bool> bound_node_done;

 // Loop over elements
 for (unsigned e=0;e<nelem;e++)
  {
   // Loop over new local nodes
   for(unsigned j=4;j<nnode;j++)
    {
     // Loop over the boundaries
     for(unsigned bo=0;bo<nbound;bo++)
      {
       // Pointer to the element's local node
       Node* loc_node_pt=finite_element_pt(e)->node_pt(j);

       // value of the map for the node and boundary specified
       bool bound_test=bound_node_done(loc_node_pt,bo);
   

       // This will have to be changed for higher-order elements
       //=======================================================

       //  These are the face nodes on the element's face 0:
       if(  (bound_test==false)  &&   ( (j==4) || (j==5) || (j==7) ))
        {
         boundary_id=Tmp_mesh_pt->face_boundary(e,0);
         if((boundary_id-1)==bo)
          {
           this->convert_to_boundary_node(loc_node_pt);
           add_boundary_node(boundary_id-1,loc_node_pt);
           bound_test=true;
          }          
        }  

       // These are the face nodes on the element's face 1:
       if (  (bound_test==false)  &&   ( (j==4) || (j==6) || (j==9) ))
        {
         boundary_id=Tmp_mesh_pt->face_boundary(e,1);
         if((boundary_id-1)==bo)
          {
           this->convert_to_boundary_node(loc_node_pt);
           add_boundary_node(boundary_id-1,loc_node_pt);
           bound_test=true;
          }          
        }

       // These are the face nodes on the element's face 2:
       if (  (bound_test==false)  &&   ( (j==5) || (j==6) || (j==8) ))
        {
         boundary_id=Tmp_mesh_pt->face_boundary(e,2);
         if((boundary_id-1)==bo)
          {
           this->convert_to_boundary_node(loc_node_pt);
           add_boundary_node(boundary_id-1,loc_node_pt);
           bound_test=true;
          }          
        }

       // These are the face nodes on the element's face 3:
       if  (  (bound_test==false)  && ( (j==7) || (j==8) || (j==9) ))
        {
         boundary_id=Tmp_mesh_pt->face_boundary(e,3);
         if((boundary_id-1)==bo)
          {
           this->convert_to_boundary_node(loc_node_pt);
           add_boundary_node(boundary_id-1,loc_node_pt);
           bound_test=true;
          }          
        }


      } // end for bo
    } //end for j
  } //end for e
 

}   // end function
 



//======================================================================
///  Setup boundary coordinate on boundary b which is
/// assumed to be planar. Boundary coordinates are the
/// x-y coordinates in the plane of that boundary with the
/// x-axis along the line from the (lexicographically)
/// "lower left" to the "upper right" node. The y axis
/// is obtained by taking the cross-product of the positive
/// x direction with the outer unit normal computed by
/// the face elements (or its negative if switch_normal is set
/// to true). Doc faces in output file.
//======================================================================
template <class ELEMENT>
void TetgenMesh<ELEMENT>::setup_boundary_coordinates(const unsigned& b,
                                                     const bool& switch_normal,
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

      
   // Loop over all nodes to find the lower left and upper right ones
   Node* lower_left_node_pt=this->boundary_node_pt(b,0);
   Node* upper_right_node_pt=this->boundary_node_pt(b,0);

   // Loop over all nodes on the boundary
   unsigned nnod=this->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {

     //Get node
     Node* nod_pt=this->boundary_node_pt(b,j);

     // Primary criterion for lower left: Does it have a lower z-coordinate?
     if (nod_pt->x(2)<lower_left_node_pt->x(2))
      {
       lower_left_node_pt=nod_pt;
      }
     // ...or is it a draw?
     else if (nod_pt->x(2)==lower_left_node_pt->x(2))
      {
       // If it's a draw: Does it have a lower y-coordinate?
       if (nod_pt->x(1)<lower_left_node_pt->x(1))
        {
         lower_left_node_pt=nod_pt;
        }
       // ...or is it a draw?
       else if (nod_pt->x(1)==lower_left_node_pt->x(1))
        {

         // If it's a draw: Does it have a lower x-coordinate?
         if (nod_pt->x(0)<lower_left_node_pt->x(0))
          {
           lower_left_node_pt=nod_pt;
          }
        }
      }

     // Primary criterion for upper right: Does it have a higher z-coordinate?
     if (nod_pt->x(2)>upper_right_node_pt->x(2))
      {
       upper_right_node_pt=nod_pt;
      }
     // ...or is it a draw?
     else if (nod_pt->x(2)==upper_right_node_pt->x(2))
      {
       // If it's a draw: Does it have a higher y-coordinate?
       if (nod_pt->x(1)>upper_right_node_pt->x(1))
        {
         upper_right_node_pt=nod_pt;
        }
       // ...or is it a draw?
       else if (nod_pt->x(1)==upper_right_node_pt->x(1))
        {

         // If it's a draw: Does it have a higher x-coordinate?
         if (nod_pt->x(0)>upper_right_node_pt->x(0))
          {
           upper_right_node_pt=nod_pt;
          }
        }
      }
    }

   // Prepare storage for boundary coordinates
   Vector<double> zeta(2);

   // Unit vector connecting lower left and upper right nodes
   Vector<double> b0(3);
   b0[0]=upper_right_node_pt->x(0)-lower_left_node_pt->x(0);
   b0[1]=upper_right_node_pt->x(1)-lower_left_node_pt->x(1);
   b0[2]=upper_right_node_pt->x(2)-lower_left_node_pt->x(2);
   
   // Normalise
   double inv_length=1.0/sqrt(b0[0]*b0[0]+b0[1]*b0[1]+b0[2]*b0[2]);
   b0[0]*=inv_length;
   b0[1]*=inv_length;
   b0[2]*=inv_length;

   // Get (outer) unit normal to first face element
   Vector<double> normal(3);
   Vector<double> s(2,0.0);
   dynamic_cast<DummyFaceElement<ELEMENT>*>(face_el_pt[0])->
    outer_unit_normal(s,normal);

   if (switch_normal)
    {
     normal[0]=-normal[0];
     normal[1]=-normal[1];
     normal[2]=-normal[2];
    }

#ifdef PARANOID
 
   // Check that all elements have the same normal
   for(unsigned e=0;e<nel;e++)
    {
     // Get (outer) unit normal to face element
     Vector<double> my_normal(3);
     dynamic_cast<DummyFaceElement<ELEMENT>*>(face_el_pt[0])->
      outer_unit_normal(s,my_normal);
     
     // Dot product should be one!
     double error=      
      normal[0]*my_normal[0]+
      normal[1]*my_normal[1]+
      normal[2]*my_normal[2];

     error-=1.0;
     if (switch_normal) error+=1.0;
     
     if (error>Tolerance_for_boundary_finding)
      {
       std::ostringstream error_message;
       error_message 
        << "Error in alingment of normals (dot product-1) " 
        << error << " for element " << e << std::endl
        << "exceeds tolerance specified by the static member data\n"
        << "TetMeshBase::Tolerance_for_boundary_finding = " 
        << Tolerance_for_boundary_finding << std::endl
        << "This usually means that the boundary is not planar.\n\n"
        << "You can suppress this error message by recompiling \n"
        << "recompiling without PARANOID or by changing the tolerance.\n";
       throw OomphLibError(error_message.str(),
                           "TetgenMesh::setup_boundary_coordinates()",
                           OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif
     
      

   

   // Cross-product to get second in-plane vector, normal to b0
   Vector<double> b1(3);
   b1[0]=b0[1]*normal[2]-b0[2]*normal[1];
   b1[1]=b0[2]*normal[0]-b0[0]*normal[2];
   b1[2]=b0[0]*normal[1]-b0[1]*normal[0];



   // Assign boundary coordinates: projection onto the axes
   for (unsigned j=0;j<nnod;j++)
    {
     //Get node
     Node* nod_pt=this->boundary_node_pt(b,j);
     
     // Difference vector to lower left corner
     Vector<double> delta(3);
     delta[0]=nod_pt->x(0)-lower_left_node_pt->x(0);
     delta[1]=nod_pt->x(1)-lower_left_node_pt->x(1);
     delta[2]=nod_pt->x(2)-lower_left_node_pt->x(2);

     // Project
     zeta[0]=delta[0]*b0[0]+delta[1]*b0[1]+delta[2]*b0[2];
     zeta[1]=delta[0]*b1[0]+delta[1]*b1[1]+delta[2]*b1[2];

#ifdef PARANOID

     // Check:
     double error=
      pow(lower_left_node_pt->x(0)+zeta[0]*b0[0]+zeta[1]*b1[0]-nod_pt->x(0),2)+
      pow(lower_left_node_pt->x(1)+zeta[0]*b0[1]+zeta[1]*b1[1]-nod_pt->x(1),2)+
      pow(lower_left_node_pt->x(2)+zeta[0]*b0[2]+zeta[1]*b1[2]-nod_pt->x(2),2);

     if (sqrt(error)>Tolerance_for_boundary_finding)
      {
       std::ostringstream error_message;
       error_message 
        << "Error in setup of boundary coordinate " 
        << sqrt(error) << std::endl
        << "exceeds tolerance specified by the static member data\n"
        << "TetMeshBase::Tolerance_for_boundary_finding = " 
        << Tolerance_for_boundary_finding << std::endl
        << "This usually means that the boundary is not planar.\n\n"
        << "You can suppress this error message by recompiling \n"
        << "recompiling without PARANOID or by changing the tolerance.\n";
       throw OomphLibError(error_message.str(),
                           "TetgenMesh::setup_boundary_coordinates()",
                           OOMPH_EXCEPTION_LOCATION);
      }
#endif

     // Set boundary coordinate
     nod_pt->set_coordinates_on_boundary(b,zeta);
    }

  }

 // Indicate that boundary coordinate has been set up
 Boundary_coordinate_exists[b]=true;

}



//========================================================================
/// Non-delaunay split elements that have three faces on a boundary
/// into sons.
//========================================================================
template <class ELEMENT>
void TetgenMesh<ELEMENT>::split_elements_in_corners(
 TimeStepper* time_stepper_pt)
{
 
 // Setup boundary lookup scheme if required
 if (!Lookup_for_elements_next_boundary_is_setup)
  {
   setup_boundary_element_info();
  }
 
 // Find out how many nodes we have along each element edge
 unsigned nnode_1d=finite_element_pt(0)->nnode_1d();
 
 // At the moment we're only able to deal with nnode_1d=2 or 3.
 if ((nnode_1d!=2)&&(nnode_1d!=3))
  {
   std::ostringstream error_message;
   error_message << "Mesh generation from tetgen currently only works\n";
   error_message << "for nnode_1d = 2 or 3. You're trying to use it\n";
   error_message << "for nnode_1d=" << nnode_1d << std::endl;
   
   throw OomphLibError(error_message.str(),
                       "TetgenMesh::split_elements_in_corners()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Map to store how many boundaries elements are located on
 std::map<FiniteElement*,unsigned> count;
 
 // Loop over all boundaries
 unsigned nb=this->nboundary();
 for (unsigned b=0;b<nb;b++)
  {
   // Loop over elements next to that boundary
   unsigned nel=this->nboundary_element(b);
   for (unsigned e=0;e<nel;e++)
    {
     // Get pointer to element
     FiniteElement* el_pt=boundary_element_pt(b,e);
     
     // Bump up counter
     count[el_pt]++;
    }
  }
 
 // To avoid having to check the map for all elements (which will
 // inflate it to the size of all elements!), move offending elements
 // into set
 std::set<FiniteElement*> elements_to_be_split;
 for (std::map<FiniteElement*,unsigned>::iterator it=count.begin();
      it!=count.end();it++)
  {
   // Get pointer to element: Key
   FiniteElement* el_pt=it->first;
   
   // Does it have to be split?
   if (it->second>2)
    {
     elements_to_be_split.insert(el_pt);
    }
  }

 // Vector for retained or newly built elements
 unsigned nel=this->nelement();
 Vector<FiniteElement*> new_el_pt;
 new_el_pt.reserve(nel);
 
 // Now loop over all elements
 for (unsigned e=0;e<nel;e++)
  {
   
   // Get pointer to element
   FiniteElement* el_pt=this->finite_element_pt(e);
   
   // Does it have to be split? I.e. is it in the set?
    std::set<FiniteElement*>::iterator it=
     std::find(elements_to_be_split.begin(),elements_to_be_split.end(),el_pt);

    // It's not in the set, so iterator has reached end
   if (it==elements_to_be_split.end())
    {
     // Carry it across
     new_el_pt.push_back(el_pt);
    }
   // It's in the set of elements to be split
   else
    {
          // Storage for new nodes -- Note: All new nodes are interior and 
     // therefore cannot be boundary nodes!
     Node* node0_pt=0;
     Node* node1_pt=0;
     Node* node2_pt=0;
     Node* node3_pt=0;
     Node* node4_pt=0;

     // Create first new element      
     FiniteElement* el1_pt=new ELEMENT;
     
     // Copy existing nodes
     el1_pt->node_pt(0)=el_pt->node_pt(0);
     el1_pt->node_pt(1)=el_pt->node_pt(1);
     el1_pt->node_pt(3)=el_pt->node_pt(3);
     if (nnode_1d==3)
      {
       el1_pt->node_pt(4)=el_pt->node_pt(4);
       el1_pt->node_pt(6)=el_pt->node_pt(6);
       el1_pt->node_pt(9)=el_pt->node_pt(9);
      }

     // Create new nodes
     node0_pt=el1_pt->construct_node(2,time_stepper_pt);
     if (nnode_1d==3)
      {
       node1_pt=el1_pt->construct_node(5,time_stepper_pt); 
       node2_pt=el1_pt->construct_node(7,time_stepper_pt); 
       node4_pt=el1_pt->construct_node(8,time_stepper_pt); 
      }
     
     
     // Create second new element
     FiniteElement* el2_pt=new ELEMENT;
     
     // Copy existing nodes
     el2_pt->node_pt(0)=el_pt->node_pt(0);
     el2_pt->node_pt(1)=el_pt->node_pt(1);
     el2_pt->node_pt(2)=el_pt->node_pt(2);
     if (nnode_1d==3)
      {
       el2_pt->node_pt(4)=el_pt->node_pt(4);
       el2_pt->node_pt(5)=el_pt->node_pt(5);
       el2_pt->node_pt(7)=el_pt->node_pt(7);
      }

     // Create new node
     if (nnode_1d==3)
      {
       node3_pt=el2_pt->construct_node(8,time_stepper_pt); 
      }
     
     // Copy existing new nodes
     el2_pt->node_pt(3)=node0_pt;
     if (nnode_1d==3)
      {
       el2_pt->node_pt(6)=node1_pt;
       el2_pt->node_pt(9)=node2_pt;
      }
     
     // Create third new element
     FiniteElement* el3_pt=new ELEMENT;
     
     // Copy existing nodes
     el3_pt->node_pt(1)=el_pt->node_pt(1);
     el3_pt->node_pt(2)=el_pt->node_pt(2);
     el3_pt->node_pt(3)=el_pt->node_pt(3);
     if (nnode_1d==3)
      {
       el3_pt->node_pt(7)=el_pt->node_pt(7);
       el3_pt->node_pt(8)=el_pt->node_pt(8);
       el3_pt->node_pt(9)=el_pt->node_pt(9);
      }

     // Copy existing new nodes
     el3_pt->node_pt(0)=node0_pt;
     if (nnode_1d==3)
      {
       el3_pt->node_pt(4)=node2_pt;
       el3_pt->node_pt(5)=node3_pt;
       el3_pt->node_pt(6)=node4_pt;
      }
     
     
     
     // Create fourth new element
     FiniteElement* el4_pt=new ELEMENT;
     
     // Copy existing nodes
     el4_pt->node_pt(0)=el_pt->node_pt(0);
     el4_pt->node_pt(2)=el_pt->node_pt(2);
     el4_pt->node_pt(3)=el_pt->node_pt(3);
     if (nnode_1d==3)
      {
       el4_pt->node_pt(5)=el_pt->node_pt(5);
       el4_pt->node_pt(6)=el_pt->node_pt(6);
       el4_pt->node_pt(8)=el_pt->node_pt(8);
      }

     // Copy existing new nodes
     el4_pt->node_pt(1)=node0_pt;
     if (nnode_1d==3)
      {
       el4_pt->node_pt(4)=node1_pt;
       el4_pt->node_pt(7)=node3_pt;
       el4_pt->node_pt(9)=node4_pt;
      }


     // Add new elements and nodes
     new_el_pt.push_back(el1_pt);
     new_el_pt.push_back(el2_pt);
     new_el_pt.push_back(el3_pt);
     new_el_pt.push_back(el4_pt);
       
     this->add_node_pt(node0_pt);
     this->add_node_pt(node1_pt);
     this->add_node_pt(node2_pt);
     this->add_node_pt(node3_pt);
     this->add_node_pt(node4_pt);


     // Set nodal positions
     for (unsigned i=0;i<3;i++)
      {
       node0_pt->x(i)=0.25*(el_pt->node_pt(0)->x(i)+
                            el_pt->node_pt(1)->x(i)+
                            el_pt->node_pt(2)->x(i)+
                            el_pt->node_pt(3)->x(i));


       if (nnode_1d==3)
        {
         node1_pt->x(i)=0.5*(el_pt->node_pt(0)->x(i)+node0_pt->x(i));
         node2_pt->x(i)=0.5*(el_pt->node_pt(1)->x(i)+node0_pt->x(i));
         node3_pt->x(i)=0.5*(el_pt->node_pt(2)->x(i)+node0_pt->x(i));
         node4_pt->x(i)=0.5*(el_pt->node_pt(3)->x(i)+node0_pt->x(i));
        }
      }

//      ofstream junk("junk.dat");
//      el_pt->output(junk);
//      el1_pt->output(junk);
//      el2_pt->output(junk);
//      el3_pt->output(junk);
//      el4_pt->output(junk);
//      junk.close();
 
     // Kill old element
     delete el_pt;
       
    }
  }
   
 // Flush element storage
 Element_pt.clear();
   
 // Copy across
 nel=new_el_pt.size();
 Element_pt.resize(nel);
 for (unsigned e=0;e<nel;e++)
  {
   Element_pt[e]=new_el_pt[e];
  }

 // Setup boundary lookup scheme again
 setup_boundary_element_info();
   
}

}
#endif
