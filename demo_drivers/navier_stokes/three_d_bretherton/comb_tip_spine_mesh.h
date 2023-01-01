//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// ?????? ASK ANDREW
//#ifndef SINGLE_LAYER_SPINE_MESH_HEADER
//#define SINGLE_LAYER_SPINE_MESH_HEADER

#ifndef COMB_TIP_MESH_HEADER
#define  COMB_TIP_MESH_HEADER

// oomph-lib includes
#include "generic/spines.h"
// This has been changed to make reference to my cubic mesh
//#include "meshes/Simple_Cubic_Mesh.h" 

#include "tip_spine_mesh.h"

//#include "my_spine_mesh.h"

using namespace oomph;

//======================================================================
/// Single-layer spine mesh class derived from standard cubic 3D mesh.
/// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
/// e.g  SpineElement<CrozierRaviartElement<2>)
///
// *****We do not introduce periodic meshes  at this point beacause we should change the complete bulk mesh*******
///
// and a surface layer of corresponding Spine interface elements, 
/// of type INTERFACE_ELEMENT, e.g. 
/// SpineSurfaceNavierStokesInterfaceElement<ELEMENT> for 3D planar problems.
///
/// This mesh has been carefully desihned so that the numeration of the nodes on the boundaries 0 and 5 (bottom and top)
/// coincides with the numeration of the spines
//======================================================================


template <class ELEMENT, class INTERFACE_ELEMENT>
class CombTipSpineMesh : public SpineMesh
{

public:

 /// Constructor: Pass number of elements in x-direction, number of
 /// The composed mesh is too complicated for giving xmin,xmax etc.. Nevertheless we keep nx, ny, nz making reference 
//   to the elements in each direction of each cubic mesh


 CombTipSpineMesh(const unsigned int  &nel, const double &alpha, const double &length, const double &height, const double &radius, const unsigned &flag,
                 TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);


 /// Access functions for pointers to interface elements
 FiniteElement* &interface_element_pt(const unsigned long &i) 
  {return Interface_element_pt[i];}

 /// Number of elements on interface

 unsigned long ninterface_element() const {return Interface_element_pt.size();}
 
 /// Access functions for pointers to elements in bulk
 FiniteElement* &bulk_element_pt(const unsigned long &i) 
  {return Bulk_element_pt[i];}


 /// Number of elements in bulk 
 unsigned long nbulk() const {return Bulk_element_pt.size();}



 void flush_spine_element_and_node_storage()
  { 
    //Clear vectors of pointers to the nodes and elements
    Node_pt.clear();
    Bulk_element_pt.clear();
    Interface_element_pt.clear();
    Element_pt.clear();
    Spine_pt.clear();
  }





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

 virtual void spine_node_update(SpineNode* spine_node_pt)
  {
 //Get fraction along the spine
   double W = spine_node_pt->fraction();
   //Get spine height
   double H = spine_node_pt->h();
   //Get the position of the node at the origin
   // SpineNode* origin_pt = dynamic_cast<SpineNode*>(spine_node_pt->spine_pt()->geom_data_pt(0));
   Vector<double> origin_pt(3);
   origin_pt[0] = spine_node_pt->spine_pt()->geom_parameter(0); 
   origin_pt[1] = spine_node_pt->spine_pt()->geom_parameter(1); 
   origin_pt[2] = spine_node_pt->spine_pt()->geom_parameter(2); 

   //set scale norm vector
   double norm = sqrt( (Xax -  origin_pt[0] ) *  (Xax -  origin_pt[0] )+ (Yax -  origin_pt[1] ) *  (Yax -  origin_pt[1] ) + 
                       (Zax -  origin_pt[2] ) *  (Zax -  origin_pt[2])  );

   //Set the x coordinate
   spine_node_pt->x(0) = origin_pt[0] + H*W* (Xax -  origin_pt[0] )/norm;
   
   //Set the y coordinate
   spine_node_pt->x(1) = origin_pt[1] + H*W* (Yax -  origin_pt[1] )/norm;
   
   //Set the value of z
   spine_node_pt->x(2) = origin_pt[2] + H*W* (Zax -  origin_pt[2] )/norm;
 
   
  }

protected:

 //Number of elements (Nx = Ny = Nel, Nz = 2*Nel) 
 unsigned int Nel;

 /// Aspect ratio, lengt, heighth and radius
 double Alpha;
 double Length;
 double Height;
 double Radius;
 unsigned  Flag;

 ///  Axis of Symmetry
 double Xax;
 double Yax;
 double Zax;


 /// Vector of pointers to element in the fluid layer
 Vector <FiniteElement *> Bulk_element_pt;

 /// Vector of pointers to interface elements
 Vector<FiniteElement *> Interface_element_pt;


// Rotation functions
 void rotate_90_zeq1_xeq0( SpineMesh* rot_mesh_pt)
  {
   unsigned long nnode = rot_mesh_pt->nnode();
   for(unsigned long i = 0; i<nnode;i++)
    {
      Node* node_pt = rot_mesh_pt->node_pt(i);
      double x_node = node_pt->x(0);
      double z_node = node_pt->x(2);
      node_pt->x(0) = 1.0 - z_node;
      node_pt->x(2) = 1.0 + x_node;
    }
   

// We rotate as well the geometric parameters which were originally set as the coordinate of the node at the origin of the spine
 unsigned long nspines = rot_mesh_pt->nspine();
  for(unsigned long i = 0; i<nspines;i++)
    {
      Spine* spine_pt = rot_mesh_pt->spine_pt(i);
      double x_spine = spine_pt->geom_parameter(0);
      double z_spine = spine_pt->geom_parameter(2);
      spine_pt->geom_parameter(0) = 1.0 - z_spine;
      spine_pt->geom_parameter(2) = 1.0 + x_spine;
    }

   //we update again the nodes
    rot_mesh_pt->node_update();
  } 



 void rotate_90_zeq1_yeq0( SpineMesh* rot_mesh_pt)
  {
   unsigned long nnode = rot_mesh_pt->nnode();
   for(unsigned long i = 0; i<nnode;i++)
    {
      Node* node_pt = rot_mesh_pt->node_pt(i);
      double y_node = node_pt->x(1);
      double z_node = node_pt->x(2);
      node_pt->x(1) =  -1.0 + z_node;
      node_pt->x(2) = 1.0 - y_node;
    }

  // We rotate as well the geometric parameters which were originally set as the coordinate of the node at the origin of the spine
 
 unsigned long nspines = rot_mesh_pt->nspine();
 for(unsigned long i = 0; i<nspines;i++)
    {
      Spine* spine_pt = rot_mesh_pt->spine_pt(i);
      double y_spine = spine_pt->geom_parameter(1);
      double z_spine = spine_pt->geom_parameter(2);
      spine_pt->geom_parameter(1) = -1.0 + z_spine;
      spine_pt->geom_parameter(2) = 1.0 - y_spine;
    }

   //we update again the nodes
    rot_mesh_pt->node_update();
  } 



 /// Helper function to actually build the single-layer spine mesh 
 /// (called from various constructors)
 virtual void build_single_layer_mesh(TimeStepper* time_stepper_pt);


 void add_side_spinemesh(unsigned int bound1, MyTipMesh<ELEMENT,INTERFACE_ELEMENT > *addmesh_pt,int *addmesh_map_boundary, int total_boundaries, unsigned flag);

 void change_boundaries(Mesh *pt_mesh,unsigned int oldbound, unsigned int newbound, unsigned int total_boundaries);
};
 


//===========================================================================
/// Constructor for spine 3D mesh: 
///
/// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
/// e.g  SpineElement<CrozierRaviartElement<2>)

//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
CombTipSpineMesh<ELEMENT, INTERFACE_ELEMENT>::CombTipSpineMesh(const unsigned int &nel, const double &alpha, const double &length, const double &height, const double &radius, 
                                                               const unsigned &flag,TimeStepper* time_stepper_pt): 
 Nel(nel), Alpha(alpha), Length(length), Height(height), Radius(radius), Flag(flag)
{
 // We've called the "generic" constructor for the RectangularQuadMesh
 // which doesn't do much...
 // Now set up the parameters that characterise the mesh geometry
 // then call the constructor


 // Periodic?
//  MySimpleCubicQuadMesh<ELEMENT >::Xperiodic = false;
 
 // Now build the mesh: 
 

 build_single_layer_mesh(time_stepper_pt);

}



//===========================================================================
/// Helper function that actually builds the single-layer spine mesh
/// based on the parameters set in the various constructors
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
void CombTipSpineMesh<ELEMENT, INTERFACE_ELEMENT>::build_single_layer_mesh(
 TimeStepper* time_stepper_pt) 
{


// Axis pf simmetry
 Xax = 0.0;
 Yax = 0.0;
 Zax = 1.0;

 unsigned n_x = Nel;
 unsigned n_y = Nel;
 unsigned n_z = Nel;

 // Build the first canion mesh: 
MyTipMesh<ELEMENT,INTERFACE_ELEMENT >* mesh1_pt = new MyTipMesh<ELEMENT,INTERFACE_ELEMENT > (n_x,n_y,n_z,0.0,1.0,-1.0,0.0,0.0,1.0,Radius,0,time_stepper_pt);

 mesh1_pt->node_update();

 //mesh1_pt->output("tip1.dat");
// We copy this mesh to our problem mesh 
 

// add the nodes to the mesh
//we reasign the spine_nodes to point to the ptoblem mesh
 for(unsigned i = 0; i< mesh1_pt->nnode(); i++)
  {
   //Set spine update flag equal to Flag for the nodes of this function
    mesh1_pt->node_pt(i)->node_update_fct_id() = Flag;
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



//We create a second mesh
 MyTipMesh<ELEMENT,INTERFACE_ELEMENT >* mesh2_pt = new MyTipMesh<ELEMENT,INTERFACE_ELEMENT > (n_x,n_y,n_z,0.0,1.0,-1.0,0.0,0.0,1.0,Radius,1,time_stepper_pt);

 //mesh2_pt->output("tip2.dat");

//Resignation pointer for the boundary conditions of the second mesh
  int addmesh_map_boundary[6];  

  for(unsigned int i =0;i<6;i++)
   {
    addmesh_map_boundary[i] = i;
   }

   addmesh_map_boundary[2] = 6; 
   addmesh_map_boundary[0] = 2; 
   addmesh_map_boundary[4] = -1; 

  add_side_spinemesh(2, mesh2_pt ,addmesh_map_boundary, 7,Flag);
//  cout<<"After adding 2 the mesh has "<<nbulk()<<" bulk elements"<<endl;



//We create a third mesh
 MyTipMesh<ELEMENT,INTERFACE_ELEMENT >* mesh3_pt = new MyTipMesh<ELEMENT,INTERFACE_ELEMENT > (n_x,n_y,n_z,0.0,1.0,-1.0,0.0,0.0,1.0,Radius,2,time_stepper_pt);

 //mesh3_pt->output("tip3.dat");

 change_boundaries( mesh3_pt,2,3,6);

 //std::cout<<"After change the boundary 2 have "
 //        <<mesh3_pt->nboundary_node(2)<<" nodes"<<std::endl;

 
//Resignation pointer for the boundary conditions of the third mesh

  for(unsigned int i =0;i<6;i++)
   {
    addmesh_map_boundary[i] = i;
   }
   addmesh_map_boundary[1] = 7;  
   addmesh_map_boundary[0] = 1;
   addmesh_map_boundary[3] = -1; 

  add_side_spinemesh(1, mesh3_pt ,addmesh_map_boundary, 8,Flag);

// cout<<"After adding 3 the mesh has "<<nbulk()<<" bulk elements"<<endl;


//We create a fourth  mesh
 MyTipMesh<ELEMENT,INTERFACE_ELEMENT >* mesh4_pt = new MyTipMesh<ELEMENT,INTERFACE_ELEMENT > (n_x,n_y,n_z,0.0,1.0,-1.0,0.0,0.0,1.0,Radius,0,time_stepper_pt);


// We miror the mesh twice
  rotate_90_zeq1_xeq0(mesh4_pt);

  // mesh4_pt->output("tip4.dat");
//Resignation pointer for the boundary conditions of the fifth mesh

  for(unsigned int i =0;i<6;i++)
   {
    addmesh_map_boundary[i] = i;
   }
 addmesh_map_boundary[0] = 2;
   addmesh_map_boundary[1] =7 ;
   addmesh_map_boundary[2] = 6; 
  addmesh_map_boundary[4] = -1; 
//cout<<"Before adding 4 the mesh has "<<nnode()<<" nodes "<<endl;

 add_side_spinemesh(6, mesh4_pt ,addmesh_map_boundary, 8,Flag);

 


//cout<<"After adding 4 the mesh has "<<nnode()<<" nodes "<<endl;
//We create a fourth mesh
  MyTipMesh<ELEMENT,INTERFACE_ELEMENT >* mesh5_pt = new MyTipMesh<ELEMENT,INTERFACE_ELEMENT > (n_x,n_y,n_z,0.0,1.0,-1.0,0.0,0.0,1.0,Radius,0,time_stepper_pt);

// We miror the mesh twice
   rotate_90_zeq1_yeq0(mesh5_pt);

   //mesh5_pt->output("tip5.dat");

  change_boundaries( mesh5_pt,2,3,6);

//Resignation pointer for the boundary conditions of the fourth mesh

  for(unsigned int i =0;i<6;i++)
   {
    addmesh_map_boundary[i] = i;
   }
 addmesh_map_boundary[0] = 1;  
   addmesh_map_boundary[1] = 6; 
  addmesh_map_boundary[3] = -1; 

// Add to the general mesh
 add_side_spinemesh(7, mesh5_pt ,addmesh_map_boundary, 8,Flag);

// cout<<"After adding 5 in  the tip mesh the boundary 1 has"<<this->nboundary_node(1)<<"nodes"<<endl;

// cout<<"After adding 5 the mesh has "<<nbulk()<<" bulk elements"<<endl;




// We create the last mesh
 MyTipMesh<ELEMENT,INTERFACE_ELEMENT >* mesh6_pt = new MyTipMesh<ELEMENT,INTERFACE_ELEMENT > (n_x,n_y,n_z,0.0,1.0,-1.0,0.0,0.0,1.0,Radius,1,time_stepper_pt);
 
   rotate_90_zeq1_xeq0(mesh6_pt);

   // mesh6_pt->output("tip6.dat");

 change_boundaries( mesh6_pt,4,1,6);

  for(unsigned int i =0;i<6;i++)
   {
    addmesh_map_boundary[i] = i;
   }
   addmesh_map_boundary[0] = 6;  
   addmesh_map_boundary[1] = -1; 
   addmesh_map_boundary[2] = 4; 

// Add to the general mesh
  add_side_spinemesh(6, mesh6_pt ,addmesh_map_boundary, 7,Flag);

  std::cout<<"The tip mesh consists of  "
           <<nbulk()<<" bulk elements and "<<ninterface_element()
           << " interface elements."<<std::endl;

// At the end we destroy the block meshes in two steps
 
// 1. We disconect the meshesh from the nodes and elements (if not we delete them at the time we delete the meshes)
  mesh1_pt->flush_spine_element_and_node_storage();
  mesh2_pt->flush_spine_element_and_node_storage();
  mesh3_pt->flush_spine_element_and_node_storage();
  mesh4_pt->flush_spine_element_and_node_storage();
  mesh5_pt->flush_spine_element_and_node_storage();
  mesh6_pt->flush_spine_element_and_node_storage();

//2. Simply delete (this function can give a segmentation fault error when there are spines still attached to the old update functions written in the constructor)
  delete mesh1_pt;
  delete mesh2_pt;
  delete mesh3_pt;
  delete mesh4_pt;
  delete mesh5_pt;
  delete mesh6_pt; 




// until now we have just worked with a square channel
// introduce the aspect ratio and different lenght (there is a factor 2 becasue 1 = halve height= b/2)
  int ncount;
  ncount = this->nnode();
  for(int i = 0; i< ncount; i++)
   {
     Node* node_pt =  this->node_pt(i);
     node_pt->x(0) =  node_pt->x(0) * Alpha/2;
     node_pt->x(1) =  node_pt->x(1) * Length;
     node_pt->x(2) =   node_pt->x(2) * Height/2;
   }
// Update also the origin of the spines
 unsigned long nspines = this->nspine();
 for(unsigned long i = 0; i<nspines;i++)
    {
      Spine* spine_pt = this->spine_pt(i);
      double x_spine = spine_pt->geom_parameter(0);
      double y_spine = spine_pt->geom_parameter(1);
      double z_spine = spine_pt->geom_parameter(2);
      spine_pt->geom_parameter(0) = x_spine *Alpha/2.0;
      spine_pt->geom_parameter(1) = y_spine * Length;
      spine_pt->geom_parameter(2) = z_spine * Height/2.0;
    }


// Update the axis of simmetry
  Zax = Zax *  Height/2.0;

// Update the spines heigth according to the new geometry
  ncount = this->nboundary_node(5);
  
  for(int i = 0; i< ncount; i++)
   {
     SpineNode* spine_node_pt =  dynamic_cast<SpineNode*>(this->boundary_node_pt(5,i));
     Spine* spine_pt = spine_node_pt->spine_pt();
      double x_spine = spine_pt->geom_parameter(0);
      double y_spine = spine_pt->geom_parameter(1);
      double z_spine = spine_pt->geom_parameter(2);
    
      double  dorg = sqrt( (spine_node_pt->x(0) -  x_spine )*(spine_node_pt->x(0) -  x_spine )+
                 (spine_node_pt->x(1) -  y_spine)*(spine_node_pt->x(1) -  y_spine)+ 
               (spine_node_pt->x(2) -  z_spine)*(spine_node_pt->x(2) -  z_spine) );
     spine_node_pt->h() = dorg;
     
   }

 
}










//Agregate a MySpineMesh to the problem mesh mesh_pt
//The function neeeds the pointer to the new added mes;, the sahred boundary of the problem mesh;
//a pointer which is a map between the boundaries in the addeed mesh and their new values when the mesh 
//is attached (we will set as convection boundary = -1 for the shared boundary), the new total value of boundaries and a flag 
// which will be used by the function update_node ( so that all the sub meshes with the same flag will be updated using the same function)

//This function can not be extended to adding every spine_mesh because the class SpineMesh lacks the bulk and interface element pointers

template<class ELEMENT, class INTERFACE_ELEMENT> void  CombTipSpineMesh<ELEMENT, INTERFACE_ELEMENT>::add_side_spinemesh(unsigned int bound1, 
                                MyTipMesh<ELEMENT,INTERFACE_ELEMENT > *addmesh_pt, int *addmesh_map_boundary, int total_boundaries, unsigned spine_flag)
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




/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//Merges this boundaries of the problem mesh. Both boundaries must contain nodes exactly in the same positions and equivalent spines.
//The nodes of one of the boundaries will be completly deleted and the ones of the other will keep on the mesh but not any more in the former boundary
//(they will remain in other boundaries different from bound1 and bound2)
//The function neeeds the boundary numbers to be merged and the end total number of boundaries (the merged boundaries can still exists but they 
// would contain  0 nodes after calling this function)

//This function can not be extended to adding every spine_mesh because the class SpineMesh lacks the bulk and interface element pointers

template<class ELEMENT, class INTERFACE_ELEMENT> void  CombTipSpineMesh<ELEMENT, INTERFACE_ELEMENT>::change_boundaries(Mesh *pt_mesh,  
                                                       unsigned int oldbound, unsigned int newbound, unsigned int total_boundaries)
{  

 for(unsigned int i = 0; i<pt_mesh->nboundary_node(oldbound);i++)
  {
   if( !(pt_mesh->boundary_node_pt(oldbound,i)->is_on_boundary(newbound)) )     //Do not add the nodes which were before on the boundary
     pt_mesh->add_boundary_node(newbound, pt_mesh->boundary_node_pt(oldbound,i) );
  }
   pt_mesh->remove_boundary_nodes(oldbound);
   pt_mesh->set_nboundary(total_boundaries);
 
}


#endif

  


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

