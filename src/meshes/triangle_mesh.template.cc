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
  Node_pt.resize(nnode_scaffold,0);

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
   
  //Number of nodes per element from the scaffold mesh
  unsigned nnod_el=Tmp_mesh_pt->finite_element_pt(0)->nnode();
   
  // Setup map to check the (pseudo-)global node number 
  // Nodes whose number is zero haven't been copied across
  // into the mesh yet.
  std::map<Node*,unsigned> global_number;
  unsigned global_count=0;
 
  // Map of Element attribute pairs
  std::map<double,Vector<FiniteElement*> > element_attribute_map;

  // Loop over elements in scaffold mesh, visit their nodes
  for (unsigned e=0;e<nelem;e++)
   {
    // Loop over all nodes in element
    for (unsigned j=0;j<nnod_el;j++)
     {
      // Pointer to node in the scaffold mesh
      Node* scaffold_node_pt=Tmp_mesh_pt->finite_element_pt(e)->node_pt(j);

      //unsigned j_global = Tmp_mesh_pt->global_node_number(e*nnod_el+j) -1;

      // Get the (pseudo-)global node number in scaffold mesh
      // (It's zero [=default] if not visited this one yet)
      unsigned j_global=global_number[scaffold_node_pt];

      // Haven't done this one yet
      if (j_global==0)
      //if(Node_pt[j_global]==0)
       {
        // Find and store the node_id in the old nodes map
        Oomph_vertex_nodes_id[global_count]=
        //Oomph_vertex_nodes_id[j_global] = 
         old_global_number[scaffold_node_pt];
        
        // Get pointer to set of mesh boundaries that this 
        // scaffold node occupies; NULL if the node is not on any boundary
        std::set<unsigned>* boundaries_pt;
        scaffold_node_pt->get_boundaries_pt(boundaries_pt);

        //Storage for the new node
        Node* new_node_pt = 0;

        //Is it on boundaries
        if(boundaries_pt!=0)
         {
          //Create new boundary node
          new_node_pt = finite_element_pt(e)->
           construct_boundary_node(j,time_stepper_pt);
          
          // Add to boundaries
          for(std::set<unsigned>::iterator it=boundaries_pt->begin();
              it!=boundaries_pt->end();++it)
           {
            add_boundary_node(*it,new_node_pt);
           }
         }
        //Build normal node
        else
         {
          //Create new normal node
          new_node_pt =  finite_element_pt(e)->
           construct_node(j,time_stepper_pt);
         }
        
        // Give it a number (not necessarily the global node 
        // number in the scaffold mesh -- we just need something
        // to keep track...)
        global_count++;
        global_number[scaffold_node_pt]=global_count;
        
        // Copy new node, created using the NEW element's construct_node
        // function into global storage, using the same global
        // node number that we've just associated with the 
        // corresponding node in the scaffold mesh
        Node_pt[global_count-1] = new_node_pt;
        //Node_pt[j_global] = new_node_pt;

        // Assign coordinates
        for(unsigned i=0;i<finite_element_pt(e)->dim();i++)
         {
          new_node_pt->x(i)=scaffold_node_pt->x(i);
         }
       }
      // This one has already been done: Copy accross
      else
       {
        //finite_element_pt(e)->node_pt(j)=Node_pt[j_global];
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

  unsigned boundary_id;   

  // Get number of nodes along element edge and dimension of element (2)
  // from first element
  unsigned n_node_1d=finite_element_pt(0)->nnode_1d();
  unsigned dim=finite_element_pt(0)->dim();

  // Storage for the local coordinate of the new node
  Vector<double> s(dim);

  // Get number of nodes in the element from first element
  unsigned n_node=finite_element_pt(0)->nnode();
  
  //Storage for each global edge of the mesh
  unsigned n_global_edge = Tmp_mesh_pt->nglobal_edge();
  Vector<Vector<Node*> > nodes_on_global_edge(n_global_edge);

  // Loop over elements
  for (unsigned e=0;e<nelem;e++)
   {
    //Cache pointers to the elements
    FiniteElement* const elem_pt = finite_element_pt(e);
    FiniteElement* const tmp_elem_pt = Tmp_mesh_pt->finite_element_pt(e);
    
    //The number of edge nodes is  3*(nnode_1d-1)
    unsigned n_edge_node = 3*(n_node_1d-1);
    
    //If there are any more nodes, these are internal and can be 
    //constructed and added directly to the mesh
    for(unsigned n=n_edge_node;n<n_node;++n)
     {
      // Create new node (it can never be a boundary node)
      Node* new_node_pt=elem_pt->construct_node(n,time_stepper_pt);
      
      // What are the node's local coordinates?
      elem_pt->local_coordinate_of_node(n,s);
      
      // Find the coordinates of the new node from the existing
      // and fully-functional element in the scaffold mesh
      for(unsigned i=0;i<dim;i++)
       {
        new_node_pt->x(i)= tmp_elem_pt->interpolated_x(s,i);
       }
      
      //Add the node to the mesh's global look-up scheme
      Node_pt.push_back(new_node_pt);
     }

    //Now loop over the mid-side edge nodes
    //Start from node number 3
    unsigned n = 3;
    
    // Loop over edges
    for(unsigned j=0;j<3;j++)
     {
      //Find the boundary id of the edge
      boundary_id = Tmp_mesh_pt->edge_boundary(e,j);

      //Find the global edge index
      unsigned edge_index = Tmp_mesh_pt->edge_index(e,j);
      
      //If the nodes on the edge have not been allocated, construct them
      if(nodes_on_global_edge[edge_index].size()==0)
       {
        //Loop over the nodes on the edge excluding the ends
        for(unsigned j2=0;j2<n_node_1d-2;++j2)
         {
          //Storage for the new node
          Node* new_node_pt = 0;
          
          //If the edge is on a boundary, construct a boundary node
          if(boundary_id > 0)
           {
            new_node_pt = 
             elem_pt->construct_boundary_node(n,time_stepper_pt);
            //Add it to the boundary
            this->add_boundary_node(boundary_id-1,new_node_pt);
           }
          //Otherwise construct a normal node
          else
           {
            new_node_pt = 
             elem_pt->construct_node(n,time_stepper_pt);
           }
          
          // What are the node's local coordinates?
          elem_pt->local_coordinate_of_node(n,s);
          
          // Find the coordinates of the new node from the existing
          // and fully-functional element in the scaffold mesh
          for(unsigned i=0;i<dim;i++)
           {
            new_node_pt->x(i)= tmp_elem_pt->interpolated_x(s,i);
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
        for(unsigned j2=0;j2<n_node_1d-2;++j2)
         {
          //Set the local node from the edge but indexed the other 
          //way around
          elem_pt->node_pt(n) = 
           nodes_on_global_edge[edge_index][n_node_1d-3-j2];
          ++n;
         }
       }
      
      //Set the elements adjacent to the boundary from the 
      //boundary id information
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

     } //end of loop over edges
   } //end of loop over elements
   
    // Lookup scheme has now been setup 
  Lookup_for_elements_next_boundary_is_setup=true;
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
      OomphLibError(error_message.str(), // hierher re-enable throw
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
      OomphLibError(error_message.str(), // hierher re-enable throw
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
    //from one end or the other of the connected segment.

    //If the boundary has a geometric object representation then
    //scale the coordinates to match those of the geometric object
    GeomObject* const geom_object_pt = this->boundary_geom_object_pt(b);
    if(geom_object_pt!=0)
     {
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
      
      // Number of polyline vertices along this boundary
      unsigned n_vertex=Polygonal_vertex_arclength_info[b].size();

      // Get polygonal vertex data
      Vector<std::pair<double,double> > polygonal_vertex_arclength=
       Polygonal_vertex_arclength_info[b];

      //If the (normal) error is small than reversed then we have the 
      //coordinate direction correct.
      //If not then we must reverse it
      bool reversed=false;
      if (error< rev_error)
       {
        // Coordinates are aligned (btw: don't delete this block -- there's
        // a final else below to catch errors!)
        reversed=false;
       }
      else if (error > rev_error)
       {
        reversed=true;

        //Reverse the limits of the boundary coordinates along the
        //geometric object
        double temp = bound_coord_limits[0];
        bound_coord_limits[0] = bound_coord_limits[1];
        bound_coord_limits[1] = temp;
        for (unsigned v=0;v<n_vertex;v++)
         {
          polygonal_vertex_arclength[v].first=
           Polygonal_vertex_arclength_info[b][v].first;

          polygonal_vertex_arclength[v].second=
           Polygonal_vertex_arclength_info[b][n_vertex-v-1].second;
         }
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

      // Re-assign boundary coordinate for the case where boundary
      // is represented by polygon
      unsigned use_old=false; 
      if (n_vertex==0) use_old=true;
      
      //Now scale the coordinates accordingly
      for (std::set<Node*>::iterator it=all_nodes_pt.begin();
           it!=all_nodes_pt.end();it++)
       {            
        Node* nod_pt=(*it);

        // Get coordinate based on arclength along polygonal repesentation
        nod_pt->get_coordinates_on_boundary(b,zeta);
        
        if (use_old)
         {
          // Boundary is actually a polygon -- simply rescale
          zeta[0] = bound_coord_limits[0] + 
           (zeta_new_range/zeta_old_range)*zeta[0];
         }
        else
         {
          // Scale such that vertex nodes stay where they were
          bool found=false;

          // Loop over vertex nodes
          for (unsigned v=1;v<n_vertex;v++)
           {
            if ((zeta[0]>=polygonal_vertex_arclength[v-1].first)&&
                (zeta[0]<=polygonal_vertex_arclength[v].first))
             {
              
              // Increment in intrinsic coordinate along geom object
              double delta_zeta=(polygonal_vertex_arclength[v].second-
                                 polygonal_vertex_arclength[v-1].second);
              // Increment in arclength along segment
              double delta_polyarc=(polygonal_vertex_arclength[v].first-
                                    polygonal_vertex_arclength[v-1].first);

              // Mapped arclength coordinate
              double zeta_new=polygonal_vertex_arclength[v-1].second+
               delta_zeta*(zeta[0]-polygonal_vertex_arclength[v-1].first)/
               delta_polyarc;
              zeta[0]=zeta_new;

              // Success!
              found=true;

              // Bail out
              break;
             }
           }
          
          // If we still haven't found it's probably the last point along
          if (!found)
           {
#ifdef PARANOID
            double diff=abs(zeta[0]-
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
               << " and the originallly stored value " 
               << polygonal_vertex_arclength[n_vertex-1].first << "\n"
               << "is " << diff << " which exceeds the threshold specified\n"
               << "in the publically modifiable variable\n"
               << "ToleranceForVertexMismatchInPolygons::Tolerable_error\n"
               << "whose current value is: "
               << ToleranceForVertexMismatchInPolygons::Tolerable_error
               << "\nPlease check your mesh carefully and increase the\n"
               << "threshold if you're sure this is appropriate\n";
              throw OomphLibError(error_stream.str(),
                                  "TriangleMesh::setup_boundary_coordinates()",
                                  OOMPH_EXCEPTION_LOCATION);
             }
#endif
            zeta[0]=polygonal_vertex_arclength[n_vertex-1].second;
           }
         }

        // Assign updated coordinate
        nod_pt->set_coordinates_on_boundary(b,zeta);
       }
     }
    else
     {
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


  double t_start_overall=TimingHelpers::timer();

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
  this->max_and_min_element_size(orig_max_area, orig_min_area);
  oomph_info << "Max/min element size in original mesh: " 
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
    this->update_polygon_using_face_mesh(this->Outer_boundary_pt);

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

    //Output the mesh after the snapping has taken place
    //tmp_new_mesh_pt->output("mesh_nodes_snapped_0.dat"); 

    // Get the TriangulateIO object associated with that mesh
    TriangulateIO tmp_new_triangulateio=
     tmp_new_mesh_pt->triangulateio_representation();
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
      double t_start=TimingHelpers::timer();
      if (Do_area_transfer_by_projection)
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
        double t_start=TimingHelpers::timer();
        Problem* tmp_problem_pt=new Problem; 
        Multi_domain_functions::setup_multi_domain_interaction
         <ELEMENT>(tmp_problem_pt,this,tmp_new_mesh_pt);      
        delete tmp_problem_pt;
        tmp_problem_pt=0;
        oomph_info  
         <<"CPU for multi-domain setup for projection of areas (iter "
         << iter << ") " << TimingHelpers::timer()-t_start 
         << std::endl;
        
        
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
       }
      // Do by direct diffused bin lookup
      else
       {
        
        
        // Adjust size of bins
        unsigned backup_bin_x=Multi_domain_functions::Nx_bin;
        unsigned backup_bin_y=Multi_domain_functions::Ny_bin;
        Multi_domain_functions::Nx_bin=Nbin_x_for_area_transfer;
        Multi_domain_functions::Ny_bin=Nbin_y_for_area_transfer;
        
        // Make a mesh as geom object representation of the temporary
        // mesh -- this also builds up the internal bin structure 
        // from which we'll recover the target areas
        MeshAsGeomObject* tmp_mesh_geom_obj_pt = 
         new MeshAsGeomObject(tmp_new_mesh_pt);
        
        // Reset
        Multi_domain_functions::Nx_bin=backup_bin_x;
        Multi_domain_functions::Ny_bin=backup_bin_y;
        
        // Fill bin by diffusion
        tmp_mesh_geom_obj_pt->fill_bin_by_diffusion();

        // Get ready for next assignment of target areas
        target_area_map.clear();
        
        // Loop over elements in current (fine) mesh and visit all
        // its integration points. Check where it's located in the bin
        // structure of the temporary mesh and pass the target area
        // to the associated element(s)
        unsigned nelem=this->nelement();
        for (unsigned e=0;e<nelem;e++)
         {
          ELEMENT* el_pt=dynamic_cast<ELEMENT*>(this->element_pt(e));
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
            tmp_mesh_geom_obj_pt->get_bin(x,bin_number,sample_point_pairs);
            
            // Did we find it?
            if (bin_number<0)
             {
              // Not even within bin boundaries... odd
#ifdef PARANOID
              std::stringstream error_message;
              error_message 
               << "Very odd -- we're looking for a point that's not even \n"
               << "located within the bin boundaries.\n";
              OomphLibWarning(error_message.str(),
                              "RefineableTriangleMesh::adapt()",
                              OOMPH_EXCEPTION_LOCATION);
#endif
             }
            else
             {
              // Pass target area to all new elements associated
              // with this bin
              unsigned n=sample_point_pairs.size();
              if (n>0)
               {
                for (unsigned ee=0;ee<n;ee++)
                 {
                  // Get ee-th new element in bin
                  GeneralisedElement* ext_el_pt=sample_point_pairs[ee].first;

                  // Use max. rather than min area of any element overlapping 
                  // the current element, otherwise we get a rapid outward 
                  // diffusion of small elements
                  target_area_map[ext_el_pt]=
                   std::max(target_area_map[ext_el_pt],
                            target_area[e]);
                }
               }
              else
               {
                std::stringstream error_message;
                error_message 
                 << "Point not found within bin structure\n";
                throw OomphLibError(error_message.str(),
                                    "RefineableTriangleMesh::adapt()",
                                    OOMPH_EXCEPTION_LOCATION);
               }
             }
           }
         }

       } // endif for projection/diffused bin
      
      oomph_info << "CPU for transfer of target areas (iter "
                 << iter << ") " << TimingHelpers::timer()-t_start 
                 << std::endl;
      
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
    
    double max_area;
    double min_area;
    this->max_and_min_element_size(max_area, min_area);
    oomph_info << "Max/min element size in adapted mesh: " 
               << max_area  << " "
               << min_area << std::endl;    
   }
  else
   {
    oomph_info << "Not enough benefit in adaptation.\n";
    Nrefined=0;
    Nunrefined=0;
   }

   oomph_info  <<"CPU for adaptation: "
               << TimingHelpers::timer()-t_start_overall << std::endl;

 }


//=========================================================================
 /// Helper function that updates the input polygon's PSLG
 /// by using the end-points of elements from FaceMesh(es) that are
 /// constructed for the boundaries associated with the segments of the
 /// polygon.
//=========================================================================
 template<class ELEMENT>
 void RefineableTriangleMesh<ELEMENT>::
 update_polygon_using_face_mesh(TriangleMeshPolygon* polygon_pt)
 {
  
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
    unsigned bound=polygon_pt->polyline_pt(p)->boundary_id();
    
    // Loop over the face elements (ordered) and add their vertices
    unsigned n_face_element = face_mesh_pt[p]->nelement();
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
    double unrefinement_tolerance=polygon_pt->polyline_unrefinement_tolerance();

    //------------------------------------------------------
    // Unrefinement 
    //------------------------------------------------------
    if (unrefinement_tolerance>0.0)
     {

      // Initialise counter that indicates at which vertex we're currently
      // considering for deletion
      unsigned counter=1;
      
      // Loop over the nodes; start with the second one and increment by two
      // this way a "pack" of three nodes will be considered for calculation:
      // the middle-node (which is to be deleted or not) and the adjacent 
      // nodes
      if (n_vertex>=3)
       {
        for(unsigned i=1;i<n_vertex-2;i+=2)
         {
          // Maths from http://www.cgafaq.info/wiki/Circle_Through_Three_Points
          double a_x=tmp_vector_vertex_node[i-1][0];
          double a_y=tmp_vector_vertex_node[i-1][1];
          double b_x=tmp_vector_vertex_node[i][0];
          double b_y=tmp_vector_vertex_node[i][1];
          double c_x=tmp_vector_vertex_node[i+1][0];
          double c_y=tmp_vector_vertex_node[i+1][1];

          double a=b_x-a_x;
          double b=b_y-a_y;
          double c=c_x-a_x;
          double d=c_y-a_y;

          double e=a*(a_x+b_x)+b*(a_y+b_y);
          double f=a*(a_x+c_x)+b*(a_y+c_y);

          double g=2.0*(a*(c_y-b_y)-b*(c_x-b_x));

          bool do_it=false;
          if (g<1.0e-14)
           {
            do_it=true;
           }
          else
           {
            double p_x=(d*e-b*f)/g;
            double p_y=(a*f-c*e)/g;
            
            double r=sqrt(pow((a_x-p_x),2)+pow((a_y-p_y),2));
            
            double rpc_x=c_x-p_x;
            double rpc_y=c_y-p_y;

            double rhalfca_x=0.5*(a_x-c_x);
            double rhalfca_y=0.5*(a_y-c_y);
            
            double halfca_squared=pow(rhalfca_x,2)+pow(rhalfca_y,2);
            double pc_squared=pow(rpc_x,2)+pow(rpc_y,2);

            double sticky_out_bit=r-sqrt(pc_squared-halfca_squared);

            // If sticky out bit divided by distance between end nodes 
            // is less than tolerance the boundary is so flat that we 
            // can safely kill the node
            if ((sticky_out_bit/(2.0*sqrt(halfca_squared)))<
                 unrefinement_tolerance)
             {
              do_it=true;
             } 
           }


          // Remove node?
          if (do_it)
           {
            tmp_vector_vertex_node[i].resize(0);
           }
          
          // Increase the counter, that indicates the number of the 
          // current middle node
          counter+=2;
         }


        
        // Special treatment for the end of the polyline:
        // If the for loop ended at the last but second node and this node
        // was not deleted, then check if the last but one node can be deleted
        if((counter==n_vertex-3)&&
           (tmp_vector_vertex_node[counter].size()!=0))
         {
          // Set the last but one node as middle node
          unsigned i=tmp_vector_vertex_node.size()-2;


          // CODE DUPLICATION -- CAN'T BE BOTHERED TO WRITE A SEPARATE
          // FUNCTION FOR THIS; PROBABLY WORTH DOING IF/WHEN THERE'S
          // A MISTAKE IN ANY OF THIS AND IT NEEDS TO BE FIXED...

          // Maths from http://www.cgafaq.info/wiki/Circle_Through_Three_Points
          double a_x=tmp_vector_vertex_node[i-1][0];
          double a_y=tmp_vector_vertex_node[i-1][1];
          double b_x=tmp_vector_vertex_node[i][0];
          double b_y=tmp_vector_vertex_node[i][1];
          double c_x=tmp_vector_vertex_node[i+1][0];
          double c_y=tmp_vector_vertex_node[i+1][1];

          double a=b_x-a_x;
          double b=b_y-a_y;
          double c=c_x-a_x;
          double d=c_y-a_y;

          double e=a*(a_x+b_x)+b*(a_y+b_y);
          double f=a*(a_x+c_x)+b*(a_y+c_y);

          double g=2.0*(a*(c_y-b_y)-b*(c_x-b_x));

          bool do_it=false;
          if (g<1.0e-14)
           {
            do_it=true;
           }
          else
           {
            double p_x=(d*e-b*f)/g;
            double p_y=(a*f-c*e)/g;
            
            double r=sqrt(pow((a_x-p_x),2)+pow((a_y-p_y),2));
            
            double rpc_x=c_x-p_x;
            double rpc_y=c_y-p_y;

            double rhalfca_x=0.5*(a_x-c_x);
            double rhalfca_y=0.5*(a_y-c_y);
            
            double halfca_squared=pow(rhalfca_x,2)+pow(rhalfca_y,2);
            double pc_squared=pow(rpc_x,2)+pow(rpc_y,2);

            double sticky_out_bit=r-sqrt(pc_squared-halfca_squared);

            // If sticky out bit divided by distance between end nodes 
            // is less than tolerance the boundary is so flat that we 
            // can safely kill the node
            if ((sticky_out_bit/(2.0*sqrt(halfca_squared)))<
                 unrefinement_tolerance)
             {
              do_it=true;
             } 
           }
          
          // Remove node?
          if (do_it)
           {
            tmp_vector_vertex_node[i].resize(0);
           }
         }
        
        // Create another vector, which will only contain entries of 
        //nodes that still exist
        Vector<Vector<double> > compact_vector;
        compact_vector.reserve(n_vertex);
        for (unsigned i=0;i<n_vertex;i++)
         {
          // If the entry was not deleted include it in the new vector
          if (tmp_vector_vertex_node[i].size()!=0)
           {
            compact_vector.push_back(tmp_vector_vertex_node[i]);
           }
         }
        
        /// Get the size of the vector that now includes all remaining nodes
        n_vertex =compact_vector.size();
        
        /// Copy back
        tmp_vector_vertex_node.resize(n_vertex);
        for(unsigned i=0;i<n_vertex;i++)
         {
          tmp_vector_vertex_node[i].resize(3);
          tmp_vector_vertex_node[i][0]=compact_vector[i][0];
          tmp_vector_vertex_node[i][1]=compact_vector[i][1];
          tmp_vector_vertex_node[i][2]=compact_vector[i][2];
         }
       } // end of if for at least three points
     } // end of unrefinement
    
    
    
    //------------------------------------------------
    /// Refinement
    //------------------------------------------------
    double refinement_tolerance=polygon_pt->polyline_refinement_tolerance();
    if (refinement_tolerance>0)
     {
      // Create a geometric object from the mesh to represent
      //the curvilinear boundary
      MeshAsGeomObject* mesh_geom_obj_pt = 
       new MeshAsGeomObject(face_mesh_pt[p]);
      
      // Get the total number of current vertices
      unsigned n_vertex=tmp_vector_vertex_node.size();
      
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
        double zeta_left=tmp_vector_vertex_node[inod][0];
        
        // Get position vector of "left" node
        Vector<double> R_left(2);
        for(unsigned i=0;i<2;i++)
         {
          R_left[i]=tmp_vector_vertex_node[inod][i+1];
         }
        
        // Get local coordinate of "right" node
        double zeta_right=tmp_vector_vertex_node[inod+1][0];
        
        // Get position vector of "right" node
        Vector<double> R_right(2);
        for(unsigned i=0;i<2;i++)
         {
          R_right[i]=tmp_vector_vertex_node[inod+1][i+1];
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
          Vector<double> new_node(3);
          new_node[0]=zeta_mid[0];
          new_node[1]=R_mid[0];
          new_node[2]=R_mid[1];
          
          // Include the "left" node in the new "temporary" vector
          extended_vector.push_back(tmp_vector_vertex_node[inod]);
          
          // Include the new node as well
          extended_vector.push_back(new_node);
         }
        else
         {
          // Include the "left" node in the new "temporary" vector
          // and move on to the next node
          extended_vector.push_back(tmp_vector_vertex_node[inod]);
         }
       } // end of loop over nodes
      
      // Add the last node to the vector
      extended_vector.push_back(tmp_vector_vertex_node[n_vertex-1]);
      
      /// Get the size of the vector that now includes all added nodes
      n_vertex=extended_vector.size();
      
      // Copy across
      tmp_vector_vertex_node.resize(n_vertex);
      for(unsigned i=0;i<n_vertex;i++)
       {
        tmp_vector_vertex_node[i].resize(3);
        tmp_vector_vertex_node[i][0]=extended_vector[i][0];
        tmp_vector_vertex_node[i][1]=extended_vector[i][1];
        tmp_vector_vertex_node[i][2]=extended_vector[i][2];
       }
      
      // Delete the allocated memory for the geometric object
      // that represents the curvilinear boundary
      delete mesh_geom_obj_pt;
      
     } // end refinement
    
    
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
    delete polygon_pt->polyline_pt(p);
    polygon_pt->polyline_pt(p) = 
     new TriangleMeshPolyLine(vector_vertex_node,bound); 
   }
  
  // Cleanup (but only the elements -- the nodes still exist in
  // the bulk mesh!
  for(unsigned p=0;p<n_polyline;p++)
   {
    face_mesh_pt[p]->flush_node_storage();
    delete face_mesh_pt[p];
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
  
  // Are we eligigible for re-distributing polyline segments between
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
    GeomObject* const geom_object_pt = this->boundary_geom_object_pt(bound);
    if(geom_object_pt!=0)
     {
      eligible_for_segment_redistribution=false;
     }

    // Create a face mesh adjacent to the fluid mesh's b-th boundary. 
    // The face mesh consists of FaceElements that may also be 
    // interpreted as GeomObjects
    face_mesh_pt[p] = new Mesh;
    this->template build_face_mesh
     <ELEMENT,FaceElementAsGeomObject>(bound, face_mesh_pt[p]);
    
    // Loop over these new face elements and tell them the boundary number
    // from the bulk fluid mesh -- this is required to they can
    // get access to the boundary coordinates!
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
       "RefineableTriangleMesh::get_face_mesh_representation()",
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
     "RefineableTriangleMesh::get_face_mesh_representation()",
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
 


//======================================================================
/// Update the PSLG that define the inner boundaries of the mesh.
///
//======================================================================
 template <class ELEMENT>
 void RefineableTriangleMesh<ELEMENT>::
 surface_remesh_for_inner_hole_boundaries(Vector<Vector<double> >
                                          &internal_point_coord)
 {
  //Loop over the number of internal boundaries
  unsigned n_hole = internal_point_coord.size();
  for(unsigned ihole=0;ihole<n_hole;ihole++)
   {
    //Cache the pointer to the polygon representation
    TriangleMeshInternalPolygon* const poly_pt 
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
      //Update the polygon associated with the ihole-th hole
      this->update_polygon_using_face_mesh(poly_pt);

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
   } //End of the action
 } //End of the loop of internal boundaries


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
