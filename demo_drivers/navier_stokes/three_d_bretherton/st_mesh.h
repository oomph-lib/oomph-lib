//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
#ifndef ST_MESH_HEADER
#define  ST_MESH_HEADER

// oomph-lib includes
#include "generic/spines.h"

// Include the Meshes which form the Saffman-Taylor mesh 
#include "meshes/simple_cubic_mesh.h" 
#include "comb_can_spine_mesh.h"
#include "comb_tip_spine_mesh.h"


using namespace oomph;




//======================================================================
/// Mesh for the Saffman-Taylor problem in a cell of arbitrary aspect ratio
/// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
/// e.g  SpineElement<CrozierRaviartElement<2>) in the boundary 5
///
/// Remarck, we set the flag of every mesh twice. it is not logical but saves programming time
//======================================================================


template <class ELEMENT, class INTERFACE_ELEMENT>
class STSpineMesh : public SpineMesh
{

public:

 /// Constructor: Pass number of elements in x-direction, number of
 /// The composed mesh is too complicated for giving xmin,xmax etc.. Nevertheless we keep nx, ny, nz making reference 
//   to the elements in each direction of each cubic mesh


 STSpineMesh(const unsigned int &nel_can, const unsigned int &nel, const unsigned int &nel_liq ,const double &alpha, 
             const double &length_can, const double &length_tip,const double &length_liq,const double &height, const double &radius, 
              TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper );


 /// Access functions for pointers to interface elements
 FiniteElement* &interface_element_pt(const unsigned long &i) 
  {return Interface_element_pt[i];}

 /// Number of elements on interface
 unsigned long ninterface_element() const {return Interface_element_pt.size();}

 /// Access functions for pointers to interface line elements
 FiniteElement* &interface_line_element_pt(const unsigned long &i) 
  {return Interface_line_element_pt[i];}


 /// Number of elements on the line interface
 unsigned long ninterfaceline() const {return Interface_line_element_pt.size();}
 
 ///Access functions for pointers to elements in bulk
 FiniteElement* &bulk_element_pt(const unsigned long &i) 
  {return Bulk_element_pt[i];}

 ///Number of elements in bulk 
 unsigned long nbulk() const {return Bulk_element_pt.size();}

// /Access functions for pointers to elements in the outlet (here identified as boundary 3)
 FiniteElement* &bulk_outlet_element_pt(const unsigned long &i) 
  {return Bulk_outlet_element_pt[i];}

//Number of outlet elements
 unsigned long nbulkoutlet() const {return Bulk_outlet_element_pt.size();} 



// /Access functions for pointers to elements in the inlet (here identified as boundary 1)
 FiniteElement* &bulk_inlet_element_pt(const unsigned long &i) 
  {return Bulk_inlet_element_pt[i];}

//Number of intlet elements
 unsigned long nbulkinlet() const {return Bulk_inlet_element_pt.size();} 


//Acces functions for the fixed coordinates of the bulk outlet and inlet elements
 
 //Face index for the outlet elements
 int face_index_outlet() {return  Face_index_outlet;}
 
 //Face index for the inlet elements
 int face_index_inlet() {return  Face_index_inlet;}

 /// General node update function implements pure virtual function 
 /// defined in SpineMesh base class and performs specific node update
 /// actions:  along vertical spines

/* virtual void spine_node_update(SpineNode* spine_node_pt)
  {
   //Get fraction along the spine
   double W = spine_node_pt->fraction();
   //Get spine height
   double H = spine_node_pt->h();
   
   //Set the value of z 
   spine_node_pt->x(2) = 0.0 + W*H;
   }*/


void full_output(std::ostream &outfile, const unsigned &n_plot)
 {
 //Loop over the elements and call their output functions
 //Assign Element_pt_range
 unsigned long Element_pt_range = Bulk_element_pt.size();
 for(unsigned long e=0;e<Element_pt_range;e++)
  {
   // Try to cast to Navier Stokes Element
   NavierStokesEquations<3>* el_pt=dynamic_cast<NavierStokesEquations<3>* >(Bulk_element_pt[e]);
   if (el_pt==0)
    {
     oomph_info << "Can't execute output(...) for non Navier Stokes Elements" 
               << std::endl;
    }
   else
    {
     el_pt->full_output(outfile,n_plot);
    }
  }
 }


void surface_output(std::ostream &outfile, const unsigned &n_plot)
 {
 //Loop over the elements and call their output functions
 //Assign Element_pt_range
 unsigned long Element_pt_range = Interface_element_pt.size();
 for(unsigned long e=0;e<Element_pt_range;e++)
  {
   // Try to cast to Navier Stokes Element
   FiniteElement* el_pt=dynamic_cast<FiniteElement* >(Interface_element_pt[e]);
   if (el_pt==0)
    {
     oomph_info << "Can't execute surface output(...) for non Finite Element elements" 
               << std::endl;
    }
   else
    {
     el_pt->output(outfile,n_plot);
    }
  }
 }



 virtual void spine_node_update(SpineNode* spine_node_pt)
  {
  
 //Get de flag
   unsigned flag = spine_node_pt->node_update_fct_id();
 
// Lets first chech whether it is a real SpinNode (not in the liquid mesh)
   if(flag != 2)
    {
     //Get fraction along the spine
     double W = spine_node_pt->fraction();
     //Get spine height
     double H = spine_node_pt->h();
     //Get the position of the node at the origin
      Vector<double> origin_pt(3);
     origin_pt[0] = spine_node_pt->spine_pt()->geom_parameter(0); 
     origin_pt[1] = spine_node_pt->spine_pt()->geom_parameter(1); 
     origin_pt[2] = spine_node_pt->spine_pt()->geom_parameter(2); 
 
  

     if(flag ==0) // Canyon mesh
      {
        //set scale norm vector
        double norm = sqrt( (Xsim -  origin_pt[0] ) *  (Xsim -  origin_pt[0] )+(Zsim -  origin_pt[2] ) *  (Zsim -  origin_pt[2])  );

        //Set th x coordinate
        spine_node_pt->x(0) = origin_pt[0] + H*W* (Xsim -  origin_pt[0] )/norm;
   
        //Set the value of z
        spine_node_pt->x(2) = origin_pt[2] + H*W* (Zsim -  origin_pt[2] )/norm;    


      }

     if(flag ==1) //Tip mesh
      {
     
        //set scale norm vector
         double norm = sqrt( (Xsim -  origin_pt[0] ) *  (Xsim -  origin_pt[0] )+ (Ysim -  origin_pt[1] ) *  (Ysim -  origin_pt[1] ) + 
                             (Zsim-  origin_pt[2] ) *  (Zsim -  origin_pt[2])  );

        //Set the x coordinate
         spine_node_pt->x(0) = origin_pt[0] + H*W* (Xsim -  origin_pt[0] )/norm;
     
         //Set the y coordinate
         spine_node_pt->x(1) = origin_pt[1] + H*W* (Ysim -  origin_pt[1] )/norm;
   
         //Set the value of z
          spine_node_pt->x(2) = origin_pt[2] + H*W* (Zsim -  origin_pt[2] )/norm;
      }
    }
  }

protected:


//Number of elements 
 unsigned int Nel_can;
 unsigned int Nel;
 unsigned int Nel_liq;

 /// Aspect ratio, length, height and radius
 double Alpha;
 double Length_can;
 double Length_tip;
 double Length_liq;
 double Height;
 double Radius; 

// Axis of simmetry of the canyon  mesh
 double Xsim;
 double Zsim;

//Point of simety of the tip mesh (in the axis of canyon mesh)
 double Ysim;


 //Face index for the outlet elements
 int Face_index_outlet;

 //Face index for the inlet elements
 int Face_index_inlet;



 /// Vector of pointers to element in the fluid layer
 Vector <FiniteElement *> Bulk_element_pt;

 /// Vector of pointers to interface elements
 Vector<FiniteElement *> Interface_element_pt;

/// Vector of pointers to the bulk outlet elements in the fluid layer
 Vector <FiniteElement *> Bulk_outlet_element_pt;

/// Vector of pointers to the bulk inlet elements in the fluid layer
 Vector <FiniteElement *> Bulk_inlet_element_pt;

// Vector of pointers to the surface elements which will generate the LinContElement
 Vector <FiniteElement *> Interface_line_element_pt;

 /// Helper function to actually build the single-layer spine mesh 
 /// (called from various constructors)
 virtual void build_single_layer_mesh(TimeStepper* time_stepper_pt);

// add side_mesh to the problem mesh
 void add_side_spinemesh(unsigned int bound1,  CombTipSpineMesh<ELEMENT,INTERFACE_ELEMENT > *addmesh_pt,int *addmesh_map_boundary, int total_boundaries, unsigned flag);

 void add_mesh(unsigned int bound1, Mesh *addmesh_pt, int *addmesh_map_boundary, int total_boundaries, unsigned  spine_flag);

};

 


//===========================================================================
/// Constructor for spine 3D mesh: 
///
/// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
/// e.g  SpineElement<CrozierRaviartElement<2>)

//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
STSpineMesh<ELEMENT, INTERFACE_ELEMENT>::STSpineMesh(const unsigned int &nel_can, const unsigned int &nel, const unsigned int &nel_liq , 
                                                     const double &alpha, const double &length_can, const double &length_tip,const double &length_liq,
                                                             const double &height, const double &radius, 
                                                              TimeStepper* time_stepper_pt): 
 Nel_can(nel_can), Nel(nel), Nel_liq(nel_liq), Alpha(alpha), Length_can(length_can), Length_tip(length_tip), Length_liq(length_liq), Height(height), Radius(radius)
{
 // We've called the "generic" constructor for the RectangularQuadMesh
 // which doesn't do much...
 // Now set up the parameters that characterise the mesh geometry
 // then call the constructor

 
 // Now build the mesh: 
 

 build_single_layer_mesh(time_stepper_pt);

}



//===========================================================================
/// Helper function that actually builds the single-layer spine mesh
/// based on the parameters set in the various constructors
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
void STSpineMesh<ELEMENT, INTERFACE_ELEMENT>::build_single_layer_mesh(
 TimeStepper* time_stepper_pt) 
{


// Axis pf simmetry
 Xsim = 0.0;
 Ysim = 0.0;
 Zsim = 1.0 * Height/2;


 // Build first canion mesh with flag 0 
 CombCanSpineMesh<ELEMENT,INTERFACE_ELEMENT >* mesh1_pt = new CombCanSpineMesh<ELEMENT,INTERFACE_ELEMENT > (Nel,Nel_can,Alpha,Length_can,Height,Radius,0,time_stepper_pt);

 mesh1_pt->node_update();

// We copy this mesh to our problem mesh 
 

// add the nodes to the mesh
//we reasign the spine_nodes to point to the ptoblem mesh
 for(unsigned i = 0; i< mesh1_pt->nnode(); i++)
  {
   //Set spine update flag equal to 0 for the nodes of this function
    mesh1_pt->node_pt(i)->node_update_fct_id() = 0;
    mesh1_pt->node_pt(i)->spine_mesh_pt() = this;
    this->add_node_pt(mesh1_pt->node_pt(i));
  }

   // add the bulk elements to the mesh
 for(unsigned i = 0; i< mesh1_pt->nbulk(); i++)
  {
    Bulk_element_pt.push_back(mesh1_pt->bulk_element_pt(i));
  }

   // add the interface elements to the mesh
for(unsigned i = 0; i< mesh1_pt->ninterface_element(); i++)
  {
    Interface_element_pt.push_back(mesh1_pt->interface_element_pt(i));
  }

// add all the elements to the mesh
for(unsigned i = 0; i< mesh1_pt->nelement(); i++)
  {
    Element_pt.push_back(mesh1_pt->element_pt(i));
  }


   // add the spines  to the mesh
for(unsigned i = 0; i< mesh1_pt->nspine(); i++)
  {
    Spine_pt.push_back(mesh1_pt->spine_pt(i));
  }

//add the boundaries
 set_nboundary( mesh1_pt->nboundary()); 

for(unsigned b =0; b< mesh1_pt->nboundary(); b++)
  {
   for(unsigned i = 0; i<mesh1_pt->nboundary_node(b); i++)
    {
     // We remove the node from the old boundary mesh ( This is not necessary and it is mainly wrotten for avoiding the warning message) 
     Node* node_pt =  mesh1_pt->boundary_node_pt(b,i );
     node_pt->remove_from_boundary(b);
     this->add_boundary_node( b , node_pt);
    }
  }

// Set the outlet bulk elements

 for(unsigned long e =0;e< mesh1_pt->nbulkoutlet();e++)
  Bulk_outlet_element_pt.push_back(mesh1_pt->bulk_outlet_element_pt(e));

 this->Face_index_outlet = mesh1_pt->face_index_outlet();

 for(unsigned long e =0;e< mesh1_pt->ninterfaceline();e++)
  Interface_line_element_pt.push_back(mesh1_pt->interface_line_element_pt(e));



//We create the tip mesh with flag 1
 CombTipSpineMesh<ELEMENT,INTERFACE_ELEMENT >* mesh2_pt = new  CombTipSpineMesh<ELEMENT,INTERFACE_ELEMENT > (Nel,Alpha,Length_tip,Height,Radius,1,time_stepper_pt);

//Resignation pointer for the boundary conditions of the second mesh
  int addmesh_map_boundary[7];  

  for(unsigned int i =0;i<7;i++)
   {
    addmesh_map_boundary[i] = i;
   } 
   addmesh_map_boundary[3] = -1; 

// Add the tip mesh
    add_side_spinemesh(1, mesh2_pt ,addmesh_map_boundary, 7,1);

    std::cout<<"Tip mesh added"<<std::endl;

//Creating the cubic mesh
   double lycubicmin = -(Length_tip + Length_liq); 
    double lycubicmax = -Length_tip; 
   unsigned Nzcubic =  2*Nel;
   double heightcubic = Height;
   double alphacubic = Alpha/2;
  
 SimpleCubicMesh<ELEMENT >* mesh3_pt = new  SimpleCubicMesh<ELEMENT >(Nel,Nel_liq,Nzcubic,0,alphacubic,lycubicmin,lycubicmax,0, heightcubic,time_stepper_pt);


// Set the inlet elementes
 //  int nxin = mesh3_pt->nx();  //It is not working in this way (ask Andrew)
 //  int nyin = mesh3_pt->ny();
 //  int nzin = mesh3_pt->nz(); 
 int nxin = Nel;
 int nyin = Nel_liq;
 int nzin = Nzcubic;  

   int jin = 0; 
    for(int k =0; k<nzin;k++)
    {
     for(int i =0; i<nxin;i++)
      {
       ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh3_pt->element_pt(i+jin*nxin+k*nxin*nyin)); 
       Bulk_inlet_element_pt.push_back(el_pt);
      }
    }

  // Constant coordinate for the face inlet elementes
   Face_index_inlet = -2;

//Resignation pointer for the boundary conditions of the third mesh   
   for(unsigned int i =0;i<7;i++)
   {
    addmesh_map_boundary[i] = i;
   } 
    addmesh_map_boundary[3] = -1;
    addmesh_map_boundary[5] = 6; 



// Add the cubic mesh with flag 2
  add_mesh(1, mesh3_pt, addmesh_map_boundary,7,2);

//   cout<<"Cubic mesh added"<<endl;

   std::cout<<"The complete mesh consists of  "<<nbulk()<< " bulk elements and "<<ninterface_element()<<" interface elements"<<std::endl;

// At the end we destroy the block meshes in two steps
 
// 1. We disconect the meshesh from the nodes and elements (if not we delete them at the time we delete the meshes)
  mesh1_pt->flush_spine_element_and_node_storage();
  mesh2_pt->flush_spine_element_and_node_storage();
  mesh3_pt->flush_element_and_node_storage();
 
//2. Simply delete (this function can give a segmentation fault error when there are spines still attached to the old update functions written in the constructor)
  delete mesh1_pt;
  delete mesh2_pt;
  delete mesh3_pt;


 
}




//Agregate a MySpineMesh to the problem mesh mesh_pt
//The function neeeds the pointer to the new added mes;, the sahred boundary of the problem mesh;
//a pointer which is a map between the boundaries in the addeed mesh and their new values when the mesh 
//is attached (we will set as convection boundary = -1 for the shared boundary), the new total value of boundaries and a flag 
// which will be used by the function update_node ( so that all the sub meshes with the same flag will be updated using the same function)

//This function can not be extended to adding every spine_mesh because the class SpineMesh lacks the bulk and interface element pointers

template<class ELEMENT, class INTERFACE_ELEMENT> void  STSpineMesh<ELEMENT, INTERFACE_ELEMENT>::add_side_spinemesh(unsigned int bound1, 
                                 CombTipSpineMesh<ELEMENT,INTERFACE_ELEMENT > *addmesh_pt, int *addmesh_map_boundary, int total_boundaries, unsigned spine_flag)
{   

 
 //Call the generic function
 MeshHelper::merge_spine_meshes(this,bound1,addmesh_pt,addmesh_map_boundary,
                                total_boundaries,spine_flag);
   
// add the bulk elements to the mesh
  for(unsigned i = 0; i< addmesh_pt->nbulk(); i++)
  {
    this->Bulk_element_pt.push_back(addmesh_pt->bulk_element_pt(i));
  }

   // add the interface elements to the mesh
for(unsigned i = 0; i< addmesh_pt->ninterface_element(); i++)
  {
    this->Interface_element_pt.push_back(addmesh_pt->interface_element_pt(i));
  }

}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Add a mesh without spines to our mesh

template<class ELEMENT, class INTERFACE_ELEMENT> void  STSpineMesh<ELEMENT, INTERFACE_ELEMENT>::add_mesh(unsigned int bound1, 
                                 Mesh *addmesh_pt, int *addmesh_map_boundary, int total_boundaries, unsigned  spine_flag)
{   
 
   int bound2 = -1;

   Mesh *mesh_pt = this;

//Map of the commom nodes from the nodes in the added mesh to the old mesh (It is needed for the added elements to point to the old nodes and not duplicate this nodes) 
  std::map<Node*,Node*> map_bound_node;

// We have to look for the shared boundary of the second mwsh. We identify it giving a value of -1 in the mappring of the old to the new boundary conditions
   for(unsigned i = 0; i<addmesh_pt->nboundary(); i++)
  {
    if(addmesh_map_boundary[i] == -1)
       bound2 = i;
   }
 
//#ifdef PARANOID
 
  if(bound2 == -1) {
   throw OomphLibError("Error setting the shared boundary conditions\n",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }


// Another control. Test that the two shared boundaries have the same number of nodes
  if(addmesh_pt->nboundary_node(bound2) != mesh_pt->nboundary_node(bound1) ) 
   {
    std::ostringstream error_stream;
    error_stream
     <<"Different number of nodes in the shared boundaries:\n" <<
     "Boundary "<<bound1<<" in the original mesh has "
     <<mesh_pt->nboundary_node(bound1)<<" nodes \n"<<
     "and Boundary "<<bound2<<" in the added mesh has "
     <<addmesh_pt->nboundary_node(bound2)<<" nodes.\n";
    
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);                        
   }

//#endif

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
       <<"Error doing the mapping between shared boundaries:\n" <<
       "it could not be found minimum distance in node "<<i<<"\n" 
       <<"Minimum found distance = "<<sqrt(dmin)<<"\n"
       <<"Position in the added mesh  = "
       <<addmesh_pt->boundary_node_pt(bound2,i)->x(0)<<"  "
       <<addmesh_pt->boundary_node_pt(bound2,i)->x(1)<<"   "
       <<addmesh_pt->boundary_node_pt(bound2,i)->x(2)<<"\n";

      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);                       
     }
     else
      { 
 
       // Creating the nodes map
        map_bound_node[addmesh_pt->boundary_node_pt(bound2,i)] =  mesh_pt->boundary_node_pt(bound1,nodecounter-1);       
          
      }
      
   
  }
   //Now we loop over the elements of the added mesh on the shared boundary and reasign the node  pointers to the ones 
  // of the original mesh. 

// REMARCK!!!!!  the boundary_element scheme was not implemented in the mesh design, so that we have to loop over all the elements


   for(unsigned int i = 0; i< addmesh_pt->nelement(); i++ )
   {
    
//Pointer to the element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>( addmesh_pt->element_pt(i));
       
    // Loop over the nodes in the element
     for(unsigned j =0;j< el_pt->nnode();j++)
     {
      //In case the node is in the boundary (and therefore in the map), we reasign the node
        if(map_bound_node[ el_pt->node_pt(j) ] != 0)
       { 
        el_pt->node_pt(j) = map_bound_node[ el_pt->node_pt(j)];   
       }
     }
  }

  
//We add the nodes which are not in the shared boundary


  unsigned long n_nodes_addmesh = addmesh_pt->nnode();
  unsigned int zaehler = 0;
   for(unsigned j=0;j<n_nodes_addmesh;j++)
  {
   if( map_bound_node[ addmesh_pt->node_pt(j)] == 0)
    { 
       SpineNode* nod_pt =  dynamic_cast<SpineNode*>(addmesh_pt->node_pt(j));
       nod_pt->spine_mesh_pt() = this;
       nod_pt->node_update_fct_id() = spine_flag;
       mesh_pt->add_node_pt( addmesh_pt->node_pt(j) );       
    }
   else
    {
       zaehler++;
    }
  }


//#ifdef PARANOID

// Another control
  if(zaehler != addmesh_pt->nboundary_node(bound2))
   {
    std::ostringstream error_stream;
    error_stream
     <<"You have added "
     <<(zaehler-addmesh_pt->nboundary_node(bound2))
     <<" nodes too much to the mesh.\n"
     <<"(This control should be removed in case we do not want to copy all the nodes of the shared boundaries)\n";
    
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);                       
   }

 
//#endif

 //We add the elements to the mesh
   
  for(unsigned i = 0; i< addmesh_pt->nelement(); i++)
  {
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>( addmesh_pt->element_pt(i));
   //  FiniteElement* el_pt =    addmesh_pt->finite_element_pt(i);
    Element_pt.push_back(el_pt);
    Bulk_element_pt.push_back(el_pt );
  }


  // We reset the number of boundaries
  mesh_pt->set_nboundary(total_boundaries);

 
 //We reset the boundary conditions of the old mesh (the shared boundary is not more in a boundary)
    mesh_pt->remove_boundary_nodes(bound1); 

 

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
           for(unsigned k =0;k<this->nboundary_node(i);k++)
           {
             if(node_pt == this->boundary_node_pt(i,k) )
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


// At the end we remove the shared nodes of the second mesh for avoiding problems with memory leaking
  
  for(unsigned long j =0;j<addmesh_pt->nboundary_node(bound2);j++)
     delete addmesh_pt->boundary_node_pt(bound2,j);

 // We reset the number of boundaries
  mesh_pt->set_nboundary(total_boundaries);

// Small control
//#ifdef PARANOID
  for(unsigned i = 0; i<mesh_pt->nboundary(); i ++)
   {
    std::cout<<"Boundary "<<i
             <<" has "<<mesh_pt->nboundary_node(i)<<" nodes."<<std::endl;
   }
//#endif

  }


#endif


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


