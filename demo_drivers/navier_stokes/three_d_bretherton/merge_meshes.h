//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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

namespace oomph
{
 namespace MeshHelper
  {
   template<class MESH1, class MESH2> 
    void merge_spine_meshes(MESH1* mesh_pt, unsigned bound1,
                            MESH2* addmesh_pt, 
                            int *addmesh_map_boundary, 
                            int total_boundaries, unsigned spine_flag)
    {   
     int bound2 = -1;
//Map of the commom nodes from the nodes in the added mesh to the old mesh (It is needed for the added elements to point to the old nodes and not duplicate this nodes) 
     std::map<Node*,Node*> map_bound_node;
     
     std::map<Spine*, Spine*> map_spines;
     
// We have to look for the shared boundary of the second mwsh. We identify it giving a value of -1 in the mappring of the old to the new boundary conditions
     for(unsigned i = 0; i<addmesh_pt->nboundary(); i++)
      {
       if(addmesh_map_boundary[i] == -1)
        bound2 = i;
      }
     
#ifdef PARANOID
     
     if(bound2 == -1) 
      {
       throw OomphLibError("Error setting the shared boundary conditions",
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }
     
     
     
// Another control. Test that the two shared boundaries have the same number of nodes
     if(addmesh_pt->nboundary_node(bound2) != mesh_pt->nboundary_node(bound1) ) 
      {
       std::ostringstream error_stream;
       error_stream 
        << "Error: different number of nodes in the shared boundaries:" 
        << "Boundary "<<bound1<<" in the original mesh has "
        <<mesh_pt->nboundary_node(bound1)<<" nodes "<<std::endl<<
        "and Boundary "<<bound2<<" in the added mesh has "
        <<addmesh_pt->nboundary_node(bound2)<<" nodes."<<std::endl;
       
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
      }
#endif

//We create the mappping
   
  //square distance between two nodes
  double d;
// Minimun distance between two nodes (only useful in case of error)
  double dmin;
  //Dimession of the space
  unsigned dim = addmesh_pt->node_pt(0)->ndim();
  
  unsigned long nodecounter;  


 
  for( unsigned long i = 0; i< addmesh_pt->nboundary_node(bound2);i++)
   {
     nodecounter =0;
     dmin = 100.0;
     do
     {
      d = 0.0;
      for(unsigned int k=0;k<dim;k++)
       {
         d += (addmesh_pt->boundary_node_pt(bound2,i)->x(k) -  mesh_pt->boundary_node_pt(bound1,nodecounter)->x(k) )* 
              (addmesh_pt->boundary_node_pt(bound2,i)->x(k) -  mesh_pt->boundary_node_pt(bound1,nodecounter)->x(k) );        
       }
      

      nodecounter++;
      if(dmin>d) dmin = d;
      
     }while((nodecounter < mesh_pt->nboundary_node(bound1)) && (d>1E-10));
     
     if((nodecounter == mesh_pt->nboundary_node(bound1)) && (d>1E-10))
     {
      std::ostringstream error_stream;
      error_stream 
       <<"Error doing the mapping between shared boundaries:\n"<< 
       "it could not be found minimum distance in node "<<i<<"\n"
       <<"Minimum found distance = "<<sqrt(dmin)<<"\n"
       <<"Position in the added mesh  = "
       <<addmesh_pt->boundary_node_pt(bound2,i)->x(0)
       <<"  "<<addmesh_pt->boundary_node_pt(bound2,i)->x(1)<<"   "
       <<addmesh_pt->boundary_node_pt(bound2,i)->x(2)<<"\n";

      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
     else
      { 
 
       // Creating the nodes map
        map_bound_node[addmesh_pt->boundary_node_pt(bound2,i)] =  mesh_pt->boundary_node_pt(bound1,nodecounter-1);       
         // Once we have indentified the corresponding nodes we map the spines
        SpineNode *add_spinenode = dynamic_cast<SpineNode*>(addmesh_pt->boundary_node_pt(bound2,i));
        SpineNode *org_spinenode = dynamic_cast<SpineNode*>( mesh_pt->boundary_node_pt(bound1,nodecounter-1) );
       
#ifdef PARANOID
        if(add_spinenode ==0)
        {
         throw OomphLibError(
          "There are nodes on the shared boundary of the added mesh which are not Spine Nodes\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
         }
      
         if(org_spinenode ==0)
         {
          throw OomphLibError(
           "There are nodes on the shared boundary of the original mesh which are not Spine Nodes\n",
           OOMPH_CURRENT_FUNCTION,
           OOMPH_EXCEPTION_LOCATION);
         } //mesh3_pt
#endif

         Spine *add_spine = add_spinenode->spine_pt();
         Spine *org_spine = org_spinenode->spine_pt();
         

        //we check the map from this  spine has not been before asigned  (one spine is pointed by many nodes)
         if(map_spines[add_spine] == 0)
          {
             map_spines[add_spine] = org_spine;
          }
         else
          {
//We assure that the nodes that belong to  one spine in one mesh are as well belonging to only one spine in the other mesh
#ifdef PARANOID
             if( map_spines[add_spine] != org_spine )
              {
               throw OomphLibError(
                "The spines map could not be performed because the shared spines are pointing to different nodes\n",
                OOMPH_CURRENT_FUNCTION,
                OOMPH_EXCEPTION_LOCATION);
              }
#endif
        
          }
      }
   
  }
  
  
  //Now we loop over the elements of the added mesh on the shared boundary and reasign the node and spine pointers to the ones 
  // of the original mesh. We need to do a distition between ELMENTS and SURFACE ELEMENT

// REMARCK!!!!!  the boundary_element scheme was not implemented in the mesh design, so that we have to loop over all the elements


   for(unsigned int i = 0; i< addmesh_pt->nbulk(); i++ )
   {
    //Pointer to the element
    FiniteElement* el_pt = addmesh_pt->bulk_element_pt(i);
    
    // Loop over the nodes in the element
    for(unsigned j =0;j< el_pt->nnode();j++)
     {
      //In case the node is in the boundary (and therefore in the map), 
      //we reasign the node
      if(map_bound_node[ el_pt->node_pt(j) ] != 0)
       { 
        el_pt->node_pt(j) = map_bound_node[ el_pt->node_pt(j)];   
       }
     }
   }


/// Now the surface elements

  for(unsigned int i = 0; i< addmesh_pt->ninterface_element(); i++ )
   {    
//Pointer to the element
    FiniteElement* surfel_pt = addmesh_pt->interface_element_pt(i);
   
    // Loop over the nodes in the element
     for(unsigned j =0;j< surfel_pt->nnode();j++)
     {
      //In case the node is in the boundary (and therefore in the map), 
      //we reasign the node
        if(map_bound_node[ surfel_pt->node_pt(j) ] != 0)
       { 
        surfel_pt->node_pt(j) = map_bound_node[ surfel_pt->node_pt(j)];   
       }
     }
   }

  
 
//We add the nodes which are not in the shared boundary
// and set the spine flag and the mesh for a later update

  unsigned long n_nodes_addmesh = addmesh_pt->nnode();
  unsigned int zaehler = 0;
   for(unsigned j=0;j<n_nodes_addmesh;j++)
  {
   if( map_bound_node[ addmesh_pt->node_pt(j)] == 0)
    { 
        addmesh_pt->node_pt(j)->node_update_fct_id() = spine_flag;
        addmesh_pt->node_pt(j)->spine_mesh_pt() = mesh_pt;
        mesh_pt->add_node_pt( addmesh_pt->node_pt(j) );       
    }
   else
    {
       zaehler++;
    }
  }

//We add the spines which are not in the shared boundary
//The ones in the shared boundary will be deleted
 unsigned long n_spines_addmesh = addmesh_pt->nspine();
  for(unsigned j=0;j<n_spines_addmesh;j++)
  {
   if( map_spines[ addmesh_pt->spine_pt(j)] == 0)
       mesh_pt->add_spine_pt(addmesh_pt->spine_pt(j));   
  }


#ifdef PARANOID

// Another control
  if(zaehler != addmesh_pt->nboundary_node(bound2))
   {
    std::ostringstream error_stream;
    error_stream
     <<"Error: you have added "
     <<(zaehler-addmesh_pt->nboundary_node(bound2))
     <<" nodes too much to the mesh.\n"
     <<"(This control should be removed in case we do not want to copy all the nodes of the shared boundaries)\n";

    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);                        
   }


#endif

 //We add the elements to the mesh

   
 
// add all the elements to the mesh
for(unsigned i = 0; i< addmesh_pt->nelement(); i++)
  {
    mesh_pt->add_element_pt(addmesh_pt->element_pt(i));
  }

 


//We reset the boundary conditions of the old mesh (the shared boundary is not more in a boundary)
    mesh_pt->remove_boundary_nodes(bound1); 


 // We reset the number of boundaries
  mesh_pt->set_nboundary(total_boundaries);



  //We reset the boundary conditions of the new nodes
  for(int i=0;i<(int)(addmesh_pt->nboundary());i++)
   {
    //Loop over the boundary nodes
     for(unsigned long j =0;j<addmesh_pt->nboundary_node(i);j++)
     {
    
         //We do not reset the commom boundary,whose nodes will be deleted  at the end
       if(i!=bound2)
        {
          //Create a pointer to the node
          Node* node_pt =  addmesh_pt->boundary_node_pt(i,j);
          
          // We remove the old boundaries in case it is not yet included in a new boundary
           bool alr_included = 0;
           for(unsigned k =0;k<mesh_pt->nboundary_node(i);k++)
           {
             if(node_pt == mesh_pt->boundary_node_pt(i,k) )
                   alr_included = 1;
           }
           if(!alr_included)
             node_pt->remove_from_boundary(i);
          
          // We test again not to include the nodes which will be deleted
           if( map_bound_node[node_pt] == 0)
           { 
              mesh_pt->add_boundary_node( addmesh_map_boundary[i], node_pt ); 
           }
            // if not we have to include to the boundarie (in the case that they were not included before)the maped nodes
           else
            {
             Node* map_node_pt =  map_bound_node[node_pt];
             if(!map_node_pt->is_on_boundary(addmesh_map_boundary[i]) )
                  mesh_pt->add_boundary_node(addmesh_map_boundary[i], map_node_pt ); 
            }
         }
       }
   }


  

// We remove the shared spines of the second mesh
 for(unsigned j=0;j<n_spines_addmesh;j++)
  {
   if( map_spines[ addmesh_pt->spine_pt(j)] != 0)
      delete addmesh_pt->spine_pt(j);
  }

 

// At the end we remove the shared nodes of the second mesh for avoiding problems with memory leaking
  
  for(unsigned long j =0;j<addmesh_pt->nboundary_node(bound2);j++)
     delete addmesh_pt->boundary_node_pt(bound2,j);



// Small control
#ifdef PARANOID
  for(unsigned i = 0; i<mesh_pt->nboundary(); i ++)
   {
    std::cout<<"Boundary "<<i<<" has "<<mesh_pt->nboundary_node(i)<<" nodes."<<std::endl;
   }
#endif
    }
  }
}
