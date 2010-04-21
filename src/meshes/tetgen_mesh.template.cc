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
#include "../generic/Telements.h"
#include "../generic/map_matrix.h"



namespace oomph
{


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



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
   
 // Build elements
 for (unsigned e=0;e<nelem;e++)
  {
   Element_pt[e]=new ELEMENT;
  }
   
 // Number of nodes per element
 unsigned nnod_el=Tmp_mesh_pt->finite_element_pt(0)->nnode();

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
       // Get pointer to set of mesh boundaries that this 
       // scaffold node occupies; NULL if the node is not on any boundary
       std::set<unsigned>* boundaries_pt;
       scaffold_node_pt->get_boundaries_pt(boundaries_pt);

       // Is it on boundaries?
       if (boundaries_pt!=0)
        {
         // Create new boundary node
         Node* new_node_pt=finite_element_pt(e)->
          construct_boundary_node(j,time_stepper_pt);
         
         // Give it a number (not necessarily the global node 
         // number in the scaffold mesh -- we just need something
         // to keep track...) 
         global_count++;
         global_number[scaffold_node_pt]=global_count;

         // Add to boundaries 
         for(std::set<unsigned>::iterator it=boundaries_pt->begin();
             it!=boundaries_pt->end();++it)
          {
           add_boundary_node(*it,new_node_pt);
          }
        }
       // Build normal node
       else
        {
         // Create new normal node
         finite_element_pt(e)->construct_node(j,time_stepper_pt); 

         // Give it a number (not necessarily the global node 
         // number in the scaffold mesh -- we just need something
         // to keep track...)
         global_count++;
         global_number[scaffold_node_pt]=global_count;
        }

       // Copy new node, created using the NEW element's construct_node
       // function into global storage, using the same global
       // node number that we've just associated with the 
       // corresponding node in the scaffold mesh
       Node_pt[global_count-1]=finite_element_pt(e)->node_pt(j);

       // Assign coordinates
       Node_pt[global_count-1]->x(0)=scaffold_node_pt->x(0);
       Node_pt[global_count-1]->x(1)=scaffold_node_pt->x(1);
       Node_pt[global_count-1]->x(2)=scaffold_node_pt->x(2);
         
      }
     // This one has already been done: Copy across
     else
      {
       finite_element_pt(e)->node_pt(j)=Node_pt[j_global-1];         
      }
    }
  }

 // At this point we've created all the elements and 
 // created their vertex nodes. Now we need to create
 // the additional (midside and internal) nodes!
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

 // Map storing the mid-side of an edge; edge identified by
 // pointers to vertex nodes in scaffold mesh
 MapMatrix<Node*,Node*> central_edge_node_pt;
 Node* edge_node1_pt=0;
 Node* edge_node2_pt=0;

 // Loop over all elements
 for (unsigned e=0;e<nelem;e++)
  {
   // Loop over the new nodes in the element and create them.
   for(unsigned j=4;j<nnode;j++)
    {
     
     // Figure out edge nodes
     switch(j)
      {

       // Node 4 is between nodes 0 and 1
      case 4:

       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(0);
       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(1);
       break;


       // Node 5 is between nodes 0 and 2
      case 5:
              
       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(0);
       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(2);
       break;

       // Node 6 is between nodes 0 and 3
      case 6:
       
       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(0);
       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(3);
       break;

       // Node 7 is between nodes 1 and 2
      case 7:
       
       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(1);
       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(2);
       break;

       // Node 8 is between nodes 2 and 3
      case 8:
       
       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(2);
       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(3);
       break;

       // Node 9 is between nodes 1 and 3
      case 9:
       
       edge_node1_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(1);
       edge_node2_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(3);
       break;

       break;

       //Error       
      default:

       throw OomphLibError("More than ten nodes in Tet Element",
                           "TetgenMesh::build_from_scaffold()",
                           OOMPH_EXCEPTION_LOCATION);
      }
     

     

     // Do we need a boundary node?
     bool need_boundary_node=false;
     // hierher Fine tune via intersection of boundary sets
     if (edge_node1_pt->is_on_boundary()&&edge_node2_pt->is_on_boundary())
      {
       need_boundary_node=true;
      }
     
     // Do we need a new node?
     if (central_edge_node_pt(edge_node1_pt,edge_node2_pt)==0)
      {       
       Node* new_node_pt=0;

       // Create new  boundary node       
       if (need_boundary_node)
        {
         new_node_pt=finite_element_pt(e)->
          construct_boundary_node(j,time_stepper_pt);
        }
       // Create new normal node
       else
        {
         new_node_pt=finite_element_pt(e)->
          construct_node(j,time_stepper_pt);
        }
       Node_pt.push_back(new_node_pt);
       
       // Now indicate existence of newly created mideside node in map 
       central_edge_node_pt(edge_node1_pt,edge_node2_pt)=new_node_pt;
       central_edge_node_pt(edge_node2_pt,edge_node1_pt)=new_node_pt;
       
       // What are the node's local coordinates?
       finite_element_pt(e)->local_coordinate_of_node(j,s);
       
       // Find the coordinates of the new node from the existing
       // and fully-functional element in the scaffold mesh
       for(unsigned i=0;i<dim;i++)
        {
         new_node_pt->x(i)=
          Tmp_mesh_pt->finite_element_pt(e)->interpolated_x(s,i);
        }
      }
     else
      {
       // Set pointer to the existing node
       finite_element_pt(e)->node_pt(j)=
        central_edge_node_pt(edge_node1_pt,edge_node2_pt);
      }

    } // end of loop over new nodes
  } //end of loop over elements
   
 
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
     // Pointer to the element's local node
     Node* loc_node_pt=finite_element_pt(e)->node_pt(j);
     
     // This will have to be changed for higher-order elements
     //=======================================================
     
     //  These are the face nodes on the element's face 0:
     if ( (j==4) || (j==5) || (j==7) )
      {
       boundary_id=Tmp_mesh_pt->face_boundary(e,0);
       if (boundary_id!=0)
        {
         if (!bound_node_done(loc_node_pt,boundary_id-1))
          {
           add_boundary_node(boundary_id-1,loc_node_pt);
           bound_node_done(loc_node_pt,boundary_id-1)=true;
          }
        }
      }          
     
     
     // These are the face nodes on the element's face 1:
     if ( (j==4) || (j==6) || (j==9) )
      {
       boundary_id=Tmp_mesh_pt->face_boundary(e,1);
       if (boundary_id!=0)
        {
         if (!bound_node_done(loc_node_pt,boundary_id-1))
          {
           add_boundary_node(boundary_id-1,loc_node_pt);       
           bound_node_done(loc_node_pt,boundary_id-1)=true;
          }
        }          
      }
     
     // These are the face nodes on the element's face 2:
     if ( (j==5) || (j==6) || (j==8) )
      {
       boundary_id=Tmp_mesh_pt->face_boundary(e,2);
       if (boundary_id!=0)
        {
         if (!bound_node_done(loc_node_pt,boundary_id-1))
          {
           add_boundary_node(boundary_id-1,loc_node_pt);
           bound_node_done(loc_node_pt,boundary_id-1)=true;
          }
        }          
      }


     // These are the face nodes on the element's face 3:
     if  ( (j==7) || (j==8) || (j==9) )
      {
       boundary_id=Tmp_mesh_pt->face_boundary(e,3);
       if (boundary_id!=0)
        {
         if (!bound_node_done(loc_node_pt,boundary_id-1))
          {
           add_boundary_node(boundary_id-1,loc_node_pt);
           bound_node_done(loc_node_pt,boundary_id-1)=true;
          }
        }          
      }

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

 // Cleanup
 unsigned n=face_el_pt.size();
 for (unsigned e=0;e<n;e++)
  {
   delete face_el_pt[e];
  }

}








//======================================================================
/// Snap boundaries specified by the IDs listed in boundary_id to
/// a quadratric surface, specified in the file 
/// quadratic_surface_file_name. This is usually used with vmtk-based
/// meshes for which oomph-lib's xda to poly conversion code produces the files
/// "quadratic_fsi_boundary.dat" and "quadratic_outer_solid_boundary.dat"
/// which specify the quadratic FSI boundary (for the fluid and the solid)
/// and the quadratic representation of the outer boundary of the solid. 
/// When used with these files, the flag switch_normal should be
/// set to true when calling the function for the outer boundary of the
/// solid. The DocInfo object can be used to label optional output
/// files. (Uses directory and label).
//======================================================================
template <class ELEMENT>
void TetgenMesh<ELEMENT>::snap_to_quadratic_surface(
 const Vector<unsigned>& boundary_id,
 const std::string& quadratic_surface_file_name,
 const bool& switch_normal,
 DocInfo& doc_info)
{

 // Aux storage for processing input
 char dummy[100];
 
 // Prepare to doc nodes that couldn't be snapped 
 std::ofstream the_file_non_snapped;
 std::string non_snapped_filename="non_snapped_nodes_"+doc_info.label()+".dat";
 
 // Read the number of nodes and elements (quadratic facets)
 std::ifstream infile(quadratic_surface_file_name.c_str(),std::ios_base::in);
 unsigned  n_node;
 infile>>n_node;

 // Ignore rest of line
 infile.getline(dummy, 100);

 // Number of quadratic facets
 unsigned  nel;
 infile>> nel;

 // Ignore rest of line
 infile.getline(dummy, 100);

 // Check that the number of elements matches
 if (nel!=boundary_id.size())
  {
   std::ostringstream error_message;
   error_message 
    << "Number of quadratic facets specified in  "
    << quadratic_surface_file_name << "is: " << nel
    << "\nThis doesn't match the number of planar boundaries \n" 
    << "specified in boundary_id which is " << boundary_id.size() 
    << std::endl;
   throw OomphLibError(error_message.str(),
                       "TetgenMesh::snap_to_quadratic_surface()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Temporary storage for face elements
 Vector<FreeStandingFaceElement<ELEMENT>*> face_el_pt;
 
 // Loop over quadratic face elements -- one for each facet
 for(unsigned e=0;e<nel;e++)
  {
   // hierher clean up 
   face_el_pt.push_back(new FreeStandingFaceElement<ELEMENT>);   
  }


 // Now build nodes
 unsigned n_dim=3;
 unsigned n_position_type=1;
 unsigned initial_n_value=0; 
 Vector<Node*> nod_pt(n_node);
 unsigned node_nmbr;
 std::map<unsigned,unsigned> node_number;
 std::ofstream node_file;
 for (unsigned j=0;j<n_node;j++)
  {
   nod_pt[j]=new BoundaryNode<Node>(n_dim,n_position_type,initial_n_value);
   infile >> nod_pt[j]->x(0);
   infile >> nod_pt[j]->x(1);
   infile >> nod_pt[j]->x(2);
   infile >> node_nmbr;
   node_number[node_nmbr]=j;
  }
 
 
 // Now assign nodes to elements -- each element represents
 // distinct boundary; assign enumeration as specified by
 // boundary_id.
 for(unsigned e=0;e<nel;e++)
  {
   unsigned nnod_el=face_el_pt[e]->nnode();
   unsigned j_global;
   for (unsigned j=0;j<nnod_el;j++)
    {
     infile >> j_global;
     face_el_pt[e]->node_pt(j)=nod_pt[node_number[j_global]]; 
     face_el_pt[e]->node_pt(j)->add_to_boundary(boundary_id[e]);
    }
   face_el_pt[e]->set_boundary_number_in_bulk_mesh(boundary_id[e]);
   face_el_pt[e]->set_nodal_dimension(3);
  }
 

 // Setup boundary coordinates for each facet, using
 // the same strategy as for the simplex boundaries 
 // (there's some code duplication here but it doesn't
 // seem worth it to rationlise this as the interface would
 // be pretty clunky).
 for (unsigned e=0;e<nel;e++)
  {
   Vector<Vector<double> >vertex_boundary_coord(3);
   
   // Loop over all nodes to find the lower left and upper right ones
   Node* lower_left_node_pt=face_el_pt[e]->node_pt(0);
   Node* upper_right_node_pt=face_el_pt[e]->node_pt(0);
   
   // Loop over all vertex nodes
   for (unsigned j=0;j<3;j++)
    {
     //Get node
     Node* nod_pt=face_el_pt[e]->node_pt(j);
     
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

   // Get (outer) unit normal to face element -- note that 
   // with the enumeration chosen in oomph-lib's xda to poly
   // conversion code the sign below is correct for the outer
   // unit normal on the FSI interface.
   Vector<double> tang1(3);
   Vector<double> tang2(3);
   Vector<double> normal(3);
   tang1[0]=face_el_pt[e]->node_pt(1)->x(0)-face_el_pt[e]->node_pt(0)->x(0);
   tang1[1]=face_el_pt[e]->node_pt(1)->x(1)-face_el_pt[e]->node_pt(0)->x(1);
   tang1[2]=face_el_pt[e]->node_pt(1)->x(2)-face_el_pt[e]->node_pt(0)->x(2);
   tang2[0]=face_el_pt[e]->node_pt(2)->x(0)-face_el_pt[e]->node_pt(0)->x(0);
   tang2[1]=face_el_pt[e]->node_pt(2)->x(1)-face_el_pt[e]->node_pt(0)->x(1);
   tang2[2]=face_el_pt[e]->node_pt(2)->x(2)-face_el_pt[e]->node_pt(0)->x(2);
   normal[0]=-(tang1[1]*tang2[2]-tang1[2]*tang2[1]);
   normal[1]=-(tang1[2]*tang2[0]-tang1[0]*tang2[2]);
   normal[2]=-(tang1[0]*tang2[1]-tang1[1]*tang2[0]);

   // Normalise
   inv_length=
    1.0/sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
   normal[0]*=inv_length;
   normal[1]*=inv_length;
   normal[2]*=inv_length;

   // Change direction -- usually for outer boundary of solid
   if (switch_normal)
    {
     normal[0]=-normal[0];
     normal[1]=-normal[1];
     normal[2]=-normal[2];
    }
   
   // Cross-product to get second in-plane vector, normal to b0
   Vector<double> b1(3);
   b1[0]=b0[1]*normal[2]-b0[2]*normal[1];
   b1[1]=b0[2]*normal[0]-b0[0]*normal[2];
   b1[2]=b0[0]*normal[1]-b0[1]*normal[0];

   // Assign boundary coordinates for vertex nodes: projection onto the axes
   for (unsigned j=0;j<3;j++)
    {
     //Get node
     Node* nod_pt=face_el_pt[e]->node_pt(j);
     
     // Difference vector to lower left corner
     Vector<double> delta(3);
     delta[0]=nod_pt->x(0)-lower_left_node_pt->x(0);
     delta[1]=nod_pt->x(1)-lower_left_node_pt->x(1);
     delta[2]=nod_pt->x(2)-lower_left_node_pt->x(2);
     
     // Project
     zeta[0]=delta[0]*b0[0]+delta[1]*b0[1]+delta[2]*b0[2];
     zeta[1]=delta[0]*b1[0]+delta[1]*b1[1]+delta[2]*b1[2];
     
     vertex_boundary_coord[j].resize(2);
     vertex_boundary_coord[j][0]=zeta[0];
     vertex_boundary_coord[j][1]=zeta[1];

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
                           "TetgenMesh::snap_to_quadratic_surface()",
                           OOMPH_EXCEPTION_LOCATION);
      }
#endif
     
     // Set boundary coordinate
     nod_pt->set_coordinates_on_boundary(boundary_id[e],zeta);
     
    }

   // Assign boundary coordinates of three midside nodes by linear 
   // interpolation (corresponding to a flat facet)

   // Node 3 is between 0 and 1
   zeta[0]=0.5*(vertex_boundary_coord[0][0]+vertex_boundary_coord[1][0]);
   zeta[1]=0.5*(vertex_boundary_coord[0][1]+vertex_boundary_coord[1][1]);
   face_el_pt[e]->node_pt(3)->set_coordinates_on_boundary(boundary_id[e],zeta);

   // Node 4 is between 1 and 2
   zeta[0]=0.5*(vertex_boundary_coord[1][0]+vertex_boundary_coord[2][0]);
   zeta[1]=0.5*(vertex_boundary_coord[1][1]+vertex_boundary_coord[2][1]);
   face_el_pt[e]->node_pt(4)->set_coordinates_on_boundary(boundary_id[e],zeta);

   // Node 5 is between 2 and 0
   zeta[0]=0.5*(vertex_boundary_coord[2][0]+vertex_boundary_coord[0][0]);
   zeta[1]=0.5*(vertex_boundary_coord[2][1]+vertex_boundary_coord[0][1]);
   face_el_pt[e]->node_pt(5)->set_coordinates_on_boundary(boundary_id[e],zeta);

  }

 
 // Loop over elements/facets = boundaries to snap
 bool success=true;
 for(unsigned b=0;b<nel;b++)
  {

   //Doc boundary coordinates on quadratic patches
   std::ofstream the_file;
   std::ofstream the_file_before;
   std::ofstream the_file_after;
   if (doc_info.doc_flag())
    {
     std::ostringstream filename;
     filename << doc_info.directory() << "/quadratic_coordinates_" 
              << doc_info.label() << b << ".dat";
     the_file.open(filename.str().c_str());
    
     std::ostringstream filename1;
     filename1 << doc_info.directory() << "/quadratic_nodes_before_"
               << doc_info.label() << b << ".dat";
     the_file_before.open(filename1.str().c_str());
     
     std::ostringstream filename2;
     filename2 << doc_info.directory() << "/quadratic_nodes_after_" 
               << doc_info.label() << b << ".dat";
     the_file_after.open(filename2.str().c_str());
    }     

   //Cast the element pointer
   FreeStandingFaceElement<ELEMENT>* el_pt=face_el_pt[b];
   
   // Doc boundary coordinate on quadratic facet representation
   if (doc_info.doc_flag())
    {
     Vector<double> s(2);
     Vector<double> zeta(2);
     Vector<double> x(3);
     unsigned n_plot=3;
     the_file << el_pt->tecplot_zone_string(n_plot);
     
     // Loop over plot points
     unsigned num_plot_points=el_pt->nplot_points(n_plot);
     for (unsigned iplot=0;iplot<num_plot_points;iplot++)
      {         
       // Get local coordinates of plot point
       el_pt->get_s_plot(iplot,n_plot,s);         
       el_pt->interpolated_zeta(s,zeta);
       el_pt->interpolated_x(s,x);
       for (unsigned i=0;i<3;i++)
        {
         the_file << x[i] << " ";
        }
       for (unsigned i=0;i<2;i++)
        {
         the_file << zeta[i] << " ";
        }
       the_file << std::endl;
      }
     el_pt->write_tecplot_zone_footer(the_file,n_plot);
   
//      std::cout << "Facet " << b << " corresponds to quadratic boundary " 
//                << boundary_id[b] << " which contains " 
//                << this->nboundary_element(boundary_id[b]) 
//                << " element[s] " << std::endl;
    }
   
   
   // Loop over bulk elements that are adjacent to quadratic boundary
   Vector<double> boundary_zeta(2);
   Vector<double> quadratic_coordinates(3);
   GeomObject* geom_obj_pt=0;
   Vector<double> s_geom_obj(2);
   unsigned nel=this->nboundary_element(boundary_id[b]);
   for (unsigned e=0;e<nel;e++)
    {    
     // Get pointer to the bulk element that is adjacent to boundary
     FiniteElement* bulk_elem_pt=this->boundary_element_pt(boundary_id[b],e);
     
     // Loop over nodes
     unsigned nnod=bulk_elem_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=bulk_elem_pt->node_pt(j);       
       if (nod_pt->is_on_boundary(boundary_id[b]))
        {
         nod_pt->get_coordinates_on_boundary(boundary_id[b],boundary_zeta);
         
         // Doc it?
         if (doc_info.doc_flag())
          {
           the_file_before << nod_pt->x(0) << " " 
                           << nod_pt->x(1) << " " 
                           << nod_pt->x(2) << " " 
                           << boundary_zeta[0]  << " " 
                           << boundary_zeta[1]  << " " 
                           << std::endl;
          }
         
         // Find local coordinate in quadratic facet
         el_pt->locate_zeta(boundary_zeta,geom_obj_pt,s_geom_obj);
         if (geom_obj_pt!=0)
          {
           geom_obj_pt->position(s_geom_obj,quadratic_coordinates);
           nod_pt->x(0)=quadratic_coordinates[0];
           nod_pt->x(1)=quadratic_coordinates[1];
           nod_pt->x(2)=quadratic_coordinates[2];
          }
         else
          {
           // Get ready to bail out below...
           success=false;

           std::ostringstream error_message;
           error_message 
            << "Warning: Couldn't find GeomObject during snapping to\n"
            << "quadratic surface for boundary " << boundary_id[b] 
            << ". I'm leaving the node where it was. Will bail out later.\n";   
           OomphLibWarning(
            error_message.str(),
            "TetgenMesh::snap_to_quadratic_surface()",
            OOMPH_EXCEPTION_LOCATION);
           if (!the_file_non_snapped.is_open())
            {
             the_file_non_snapped.open(non_snapped_filename.c_str());
            }
           the_file_non_snapped << nod_pt->x(0) << " " 
                                << nod_pt->x(1) << " " 
                                << nod_pt->x(2) << " " 
                                << boundary_zeta[0]  << " " 
                                << boundary_zeta[1]  << " " 
                                << std::endl;
          }
         
         // Doc it?
         if (doc_info.doc_flag())
          {
           the_file_after << nod_pt->x(0) << " " 
                          << nod_pt->x(1) << " " 
                          << nod_pt->x(2) << " " 
                          << boundary_zeta[0]  << " " 
                          << boundary_zeta[1]  << " " 
                          << std::endl;
          }
        }
      }
    }
  
   // Close doc file
   the_file.close();
   the_file_before.close();
   the_file_after.close();
  }

 // Bail out?
 if (!success)
  {
   std::ostringstream error_message;
   error_message 
    << "Warning: Couldn't find GeomObject during snapping to\n"
    << "quadratic surface. Bailing out.\n"
    << "Nodes that couldn't be snapped are contained in \n"
    << "file: " << non_snapped_filename << ".\n"
    << "This problem may arise because the switch_normal flag was \n"
    << "set wrongly.\n";
   throw OomphLibError(
    error_message.str(),
    "TetgenMesh::snap_to_quadratic_surface()",
    OOMPH_EXCEPTION_LOCATION);
  }
 
 // Close
 if (!the_file_non_snapped.is_open())
  {
   the_file_non_snapped.close();
  }

 // Kill auxiliary FreeStandingFaceElements
 for(unsigned e=0;e<nel;e++)
  {
   delete face_el_pt[e];
   face_el_pt[e]=0;
  }

 // Kill boundary nodes
 unsigned nn=nod_pt.size();
 for (unsigned j=0;j<nn;j++)
  {
   delete nod_pt[j];
  }

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
     node0_pt=el1_pt->construct_boundary_node(2,time_stepper_pt);
     if (nnode_1d==3)
      {
       node1_pt=el1_pt->construct_boundary_node(5,time_stepper_pt); 
       node2_pt=el1_pt->construct_boundary_node(7,time_stepper_pt); 
       node4_pt=el1_pt->construct_boundary_node(8,time_stepper_pt); 
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
       node3_pt=el2_pt->construct_boundary_node(8,time_stepper_pt); 
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
