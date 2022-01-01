//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
#ifndef COMB_CAN_MESH_HEADER
#define  COMB_CAN_MESH_HEADER

// Standard includes
#include <algorithm>
#include <iostream>
#include <cmath>
#include <map>


// oomph-lib includes
#include "generic/spines.h"
// This has been changed to make reference to my cubic mesh
//#include "meshes/Simple_Cubic_Mesh.h" 

#include "canyon_spine_mesh.h"



//#include "my_spine_mesh.h"
using namespace oomph;


//======================================================================
/// This mesh combines a number of CanyonMeshes into a Canyon around 
/// an entire quarter tube.
//======================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
class CombCanSpineMesh : public SpineMesh
{

public:

 /// Constructor: Pass number of elements in x-direction, number of
 /// The composed mesh is too complicated for giving xmin,xmax etc.. Nevertheless we keep nx, ny, nz making reference 
//   to the elements in each direction of each cubic mesh
 CombCanSpineMesh(const unsigned int &nel_xz, const unsigned int &nel_y, const double &alpha, const double &length, const double &height, 
                  const double &radius, const unsigned int &flag, TimeStepper* time_stepper_pt= &Mesh::Default_TimeStepper);
 
 
 /// Access functions for pointers to interface elements
 FiniteElement* &interface_element_pt(const unsigned long &i) 
  {return Interface_element_pt[i];}

 /// Access functions for pointers to elements in bulk
 FiniteElement* &bulk_element_pt(const unsigned long &i) 
  {return Bulk_element_pt[i];}

// /Access functions for pointers to elements in the outlet (here identified as boundary 3)
 FiniteElement* &bulk_outlet_element_pt(const unsigned long &i) 
  {return Bulk_outlet_element_pt[i];}


// /Access functions for pointers to interface elements in the outlet (here identified as boundary 3)
 FiniteElement* & interface_line_element_pt(const unsigned long &i) 
  {return  Interface_line_element_pt[i];}

//Acces functions for the fixed coordinates of the bulk outlet elements
 
 //Face index for the outlet elements
 int face_index_outlet() {return  Face_index_outlet;}

 /// Number of elements in bulk 
 unsigned long nbulk() const {return Bulk_element_pt.size();}

 /// Number of elements on interface
 unsigned long ninterface_element() const {return Interface_element_pt.size();}

//Number of outlet elements
 unsigned long nbulkoutlet() const {return Bulk_outlet_element_pt.size();} 

// Number of interface line elements
 unsigned long ninterfaceline() const {return  Interface_line_element_pt.size();} 

 void flush_spine_element_and_node_storage()
  { 
    //Clear vectors of pointers to the nodes and elements
    Node_pt.clear();
    Bulk_element_pt.clear();
    Interface_element_pt.clear();
    Element_pt.clear();
    Spine_pt.clear();
  }


 virtual void spine_node_update(SpineNode* spine_node_pt)
  {
   //Get fraction along the spine
   double W = spine_node_pt->fraction();
   //Get spine height
   double H = spine_node_pt->h();
  
  //Get the position of the node at the origin
   Vector<double> origin_pt(3);
   origin_pt[0] = spine_node_pt->spine_pt()->geom_parameter(0); 
   origin_pt[2] = spine_node_pt->spine_pt()->geom_parameter(2); 

   //set scale norm vector
   double norm = 
    sqrt( (Xsim -  origin_pt[0] ) *  (Xsim -  origin_pt[0] ) +
          (Zsim -  origin_pt[2] ) *  (Zsim -  origin_pt[2])  );
   
   //Set th x coordinate
   spine_node_pt->x(0) = origin_pt[0] + H*W* (Xsim -  origin_pt[0] )/norm;
   
   //Set the value of z
   spine_node_pt->x(2) = origin_pt[2] + H*W* (Zsim -  origin_pt[2] )/norm; 
 
  }

protected:


//Number of elements in x and z direction
 unsigned int Nel_xz; //elements in x direction = N, elements in z direction = 2*N
 unsigned int Nel_y; //Number of elements in y direction

 /// Aspect ratio, length, height and radius
 double Alpha;
 double Length;
 double Height;
 double Radius;
 unsigned int Flag;

 ///  Axis of Symmetry
 double Xsim;
 double Zsim;

 //Face index for the outlet elements
 int Face_index_outlet;



/// Vector of pointers to element in the fluid layer
 Vector <FiniteElement *> Bulk_element_pt;

 /// Vector of pointers to interface elements
 Vector<FiniteElement *> Interface_element_pt;

/// Vector of pointers to the bulk outlet elements in the fluid layer
 Vector <FiniteElement *> Bulk_outlet_element_pt;

// Vector of pointers to the surface elements which will generate the LinContElement
 Vector <FiniteElement *> Interface_line_element_pt;

 void rotate_90( SpineMesh* rot_mesh_pt)
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

// Add the outlet elements from a quadratic mesh
 void add_outlet_bulk_elements(SimpleCubicMesh<ELEMENT > *add_mesh_pt)
  {
   int nx = add_mesh_pt->nx();
   int ny = add_mesh_pt->ny();
   int nz = add_mesh_pt->nz(); 
   int j = ny -1; 
    for(int k =0; k<nz;k++)
    {
     for(int i =0; i<nx;i++)
      {
       ELEMENT *el_pt = dynamic_cast<ELEMENT*>(add_mesh_pt->element_pt(i+j*nx+k*nx*ny)); 
       Bulk_outlet_element_pt.push_back(el_pt);
      }
    }

  }

 //Add the Interface Elements which will generate the LineCont Elements
 void add_line_interface_elements( MyCanyonMesh<ELEMENT,INTERFACE_ELEMENT > *add_mesh_pt )
  {
     int nx = add_mesh_pt->nx();
     int ny = add_mesh_pt->ny(); 
     int j = ny -1; 
     for(int i =0; i<nx;i++)
      {
        INTERFACE_ELEMENT *el_pt = dynamic_cast<INTERFACE_ELEMENT*>(add_mesh_pt->interface_element_pt(i+j*nx)); 
        Interface_line_element_pt.push_back(el_pt); 
      }
  }

 /// Helper function to actually build the single-layer spine mesh 
 /// (called from various constructors)
 virtual void build_single_layer_mesh(TimeStepper* time_stepper_pt);

// add side_mesh to the problem mesh
 void add_side_spinemesh(unsigned int bound1, MyCanyonMesh<ELEMENT,INTERFACE_ELEMENT > *addmesh_pt,int *addmesh_map_boundary, int total_boundaries, unsigned flag);


};

 


//===========================================================================
/// Constructor for spine 3D mesh: 
///
/// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
/// e.g  SpineElement<CrozierRaviartElement<2>)

//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
CombCanSpineMesh<ELEMENT, INTERFACE_ELEMENT>::CombCanSpineMesh(const unsigned int &nel_xz, const unsigned int &nel_y, 
                                                              const double &alpha, const double &length, const double &height, const double &radius, 
                                                              const unsigned int &flag, TimeStepper* time_stepper_pt): 
 Nel_xz(nel_xz), Nel_y(nel_y), Alpha(alpha), Length(length), Height(height), Radius(radius), Flag(flag)
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
void CombCanSpineMesh<ELEMENT, INTERFACE_ELEMENT>::build_single_layer_mesh(
 TimeStepper* time_stepper_pt) 
{


// Axis pf simmetry
 Xsim = 0.0;
 Zsim = 1.0;

// Fixed coordinate for the outlet bulk elements

 Face_index_outlet = 2;

 unsigned n_x = Nel_xz;
 unsigned n_y = Nel_y;
 unsigned n_z = Nel_xz;

 // Build the first canion mesh: 
MyCanyonMesh<ELEMENT,INTERFACE_ELEMENT >* mesh1_pt = new MyCanyonMesh<ELEMENT,INTERFACE_ELEMENT > (n_x,n_y,n_z,0.0,1.0,0.0,1.0,0.0,1.0,Radius,0,time_stepper_pt);


//mesh1_pt->output("can1.dat"); 

// We copy this mesh to our problem mesh 
 

// add the nodes to the mesh
//we reasign the spine_nodes to point to the ptoblem mesh
 for(unsigned i = 0; i< mesh1_pt->nnode(); i++)
  {
   //Set spine update flag equal to 0 for the nodes of this function
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

// add the outlet bulk and interface line elements
 add_outlet_bulk_elements(mesh1_pt);
 add_line_interface_elements(mesh1_pt);

//We create a second mesh
 MyCanyonMesh<ELEMENT,INTERFACE_ELEMENT >* mesh2_pt = new MyCanyonMesh<ELEMENT,INTERFACE_ELEMENT > (n_x,n_y,n_z,0.0,1.0,0.0,1.0,0.0,1.0,Radius,1,time_stepper_pt);


 //mesh2_pt->output("can2.dat"); 

//Resignation pointer for the boundary conditions of the second mesh
  int addmesh_map_boundary[6];  

  for(unsigned int i =0;i<6;i++)
   {
    addmesh_map_boundary[i] = i;
   }

   addmesh_map_boundary[0] = 2; 
   addmesh_map_boundary[2] = 6; 
   addmesh_map_boundary[4] = -1; 

// add the outlet bulk and interface line elements
 add_outlet_bulk_elements(mesh2_pt);
 add_line_interface_elements(mesh2_pt);

  add_side_spinemesh(2, mesh2_pt ,addmesh_map_boundary, 7,Flag);


//We create a third mesh
 MyCanyonMesh<ELEMENT,INTERFACE_ELEMENT >* mesh3_pt = new MyCanyonMesh<ELEMENT,INTERFACE_ELEMENT > (n_x,n_y,n_z,0.0,1.0,0.0,1.0,0.0,1.0,Radius,0,time_stepper_pt);

//We rotate the mesh
 rotate_90(mesh3_pt);
 
 //mesh3_pt->output("can3.dat"); 

//Resignation pointer for the boundary conditions of the third mesh

  for(unsigned int i =0;i<6;i++)
   {
    addmesh_map_boundary[i] = i;
   }

   addmesh_map_boundary[2] = 6; 
   addmesh_map_boundary[0] = 2; 
   addmesh_map_boundary[4] = -1; 

  add_side_spinemesh(6, mesh3_pt ,addmesh_map_boundary, 7,Flag);

// add the outlet bulk and interface line elements
 add_outlet_bulk_elements(mesh3_pt);
 add_line_interface_elements(mesh3_pt);

//We create a fourth mesh
 MyCanyonMesh<ELEMENT,INTERFACE_ELEMENT >* mesh4_pt = new MyCanyonMesh<ELEMENT,INTERFACE_ELEMENT > (n_x,n_y,n_z,0.0,1.0,0.0,1.0,0.0,1.0,Radius,1,time_stepper_pt);

//We rotate the mesh
  rotate_90(mesh4_pt);
 
  //mesh4_pt->output("can4.dat"); 

//Resignation pointer for the boundary conditions of the third mesh

  for(unsigned int i =0;i<6;i++)
   {
    addmesh_map_boundary[i] = i;
   }

   addmesh_map_boundary[2] = 4; 
   addmesh_map_boundary[0] = 6;
   addmesh_map_boundary[4] = -1;  

// add the outlet bulk and interface line elements
 add_outlet_bulk_elements(mesh4_pt);
 add_line_interface_elements(mesh4_pt);

  add_side_spinemesh(6, mesh4_pt ,addmesh_map_boundary, 7,Flag);


  std::cout
   <<"The canyon mesh consists of  "
   <<nbulk()<<" bulk elements and "
   <<ninterface_element()<< " interface elements."<<std::endl;


// At the end we destroy the block meshes in two steps
 
// 1. We disconect the meshesh from the nodes and elements (if not we delete them at the time we delete the meshes)
  mesh1_pt->flush_spine_element_and_node_storage();
  mesh2_pt->flush_spine_element_and_node_storage();
  mesh3_pt->flush_spine_element_and_node_storage();
  mesh4_pt->flush_spine_element_and_node_storage();

//2. Simply delete (this function can give a segmentation fault error when there are spines still attached to the old update functions written in the constructor)
  delete mesh1_pt;
  delete mesh2_pt;
  delete mesh3_pt;
  delete mesh4_pt;


// until now we have just worked with a square channel
// introduce the aspect ratio and different lenght (there is a factor 2 becasue 1 = halve height= b/2)
  int ncount;
  ncount = this->nnode();
  for(int i = 0; i< ncount; i++)
   {
     Node* node_pt =  this->node_pt(i);
     node_pt->x(0) =  node_pt->x(0) * Alpha/2;
     node_pt->x(1) =  node_pt->x(1) * Length;
     node_pt->x(2) =  node_pt->x(2) * Height/2;
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


  Zsim = Zsim * Height/2.0;

// Update the spines heigth according to the new geometry
  ncount = this->nboundary_node(5);
  
  for(int i = 0; i< ncount; i++)
   {
      SpineNode* spine_node_pt =  dynamic_cast<SpineNode*>(this->boundary_node_pt(5,i));
     Spine* spine_pt = spine_node_pt->spine_pt();
      double x_spine = spine_pt->geom_parameter(0);
      double z_spine = spine_pt->geom_parameter(2);

     double dorg = sqrt( (spine_node_pt->x(0) -  x_spine)*(spine_node_pt->x(0) -  x_spine)+ 
               (spine_node_pt->x(2) -  z_spine)*(spine_node_pt->x(2) -  z_spine) );
     spine_node_pt->h() = dorg;
     
   }

   //Get the position of the node at the origin
//   SpineNode* origin_pt = dynamic_cast<SpineNode*>(spine_node_pt->spine_pt()->geom_data_pt(0));


}




//Agregate a MySpineMesh to the problem mesh mesh_pt
//The function neeeds the pointer to the new added mes;, the sahred boundary of the problem mesh;
//a pointer which is a map between the boundaries in the addeed mesh and their new values when the mesh 
//is attached (we will set as convection boundary = -1 for the shared boundary), the new total value of boundaries and a flag 
// which will be used by the function update_node ( so that all the sub meshes with the same flag will be updated using the same function)

//This function can not be extended to adding every spine_mesh because the class SpineMesh lacks the bulk and interface element pointers

template<class ELEMENT, class INTERFACE_ELEMENT> void  CombCanSpineMesh<ELEMENT, INTERFACE_ELEMENT>::add_side_spinemesh(unsigned int bound1, 
                                MyCanyonMesh<ELEMENT,INTERFACE_ELEMENT > *addmesh_pt, int *addmesh_map_boundary, int total_boundaries, unsigned spine_flag)
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



#endif


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


