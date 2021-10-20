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
#ifndef MY_CANYON_SPINE_MESH_HEADER
#define  MY_CANYON_SPINE_MESH_HEADER

// oomph-lib includes
#include "generic/spines.h"
// This has been changed to make reference to my cubic mesh
//#include "meshes/simple_cubic_mesh.h" 


using namespace oomph;
//======================================================================
///Quarter of a canyon mesh
///Mesh with outer part assumed to be a straight line
///and inner part a cylinder with given radius, but axially uniform.
///The outer face is given by the plane z = zmin 
///The cylinder axis is given by the line x = xmin, z = zmax
///The interface elements are constructed on the inner (cylindrical)
///boundary.
//======================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
class MyCanyonMesh : public  SimpleCubicMesh<ELEMENT >, 
                          public SpineMesh
{

public:

 ///  Constructor: Pass number of elements in x-direction, number of
 /// elements in y-direction, number of elements in z-direciton,
 /// the maximum and minimum values of x, y and z. The radius of the 
 /// inner cylinder is R and the rotation flag specifies whether the mesh
 /// is aligned with the lower wall (z=zmin) or the side wall (x=max)
 /// The later is accomplished by a simple rotation of all coordinates 
 /// through a right-angle.
 MyCanyonMesh(const unsigned &nx, 
              const unsigned &ny,
              const unsigned &nz,
              const double &x_min,
              const double &x_max,
              const double &y_min,
              const double &y_max,
              const double &z_min,
              const double &z_max,
              const double &R,
              const unsigned &rotation_flag,
              TimeStepper* time_stepper_pt=
              &Mesh::Default_TimeStepper);

 
 /// Access functions for pointers to interface elements
 FiniteElement* &interface_element_pt(const unsigned long &i) 
  {return Interface_element_pt[i];}

 /// Number of elements on interface
 unsigned long ninterface_element() const 
  {return Interface_element_pt.size();}
 
 ///Access functions for pointers to elements in bulk
 FiniteElement* &bulk_element_pt(const unsigned long &i) 
  {return Bulk_element_pt[i];}
 
 ///Number of elements in bulk 
 unsigned long nbulk() const {return Bulk_element_pt.size();}
 
// Reurn radius of the canyon
 double radius() {return Radius;}
 
// Change radius
 void change_radius(double R)
  {
   double scaled_R = R/(this->Zmax - this->Zmin);
   if ( (scaled_R<0.95) && (scaled_R>0.05) )
    Radius = R;
   else
    {
     std::ostringstream error_stream;
     error_stream 
      <<
      "We can not constract so small elements so close to the origin\n"
      << "Choose a different radius (actual scaled radius = "<<scaled_R<<
       std::endl;
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
  }


 ///  General node update function implements pure virtual function 
 /// defined in SpineMesh base class and performs specific node update
 /// actions:  along vertical spines
 virtual void spine_node_update(SpineNode* spine_node_pt)
  {

  //Get fraction along the spine
   double W = spine_node_pt->fraction();
   //Get spine height
   double H = spine_node_pt->h();
  
   //Get the position of the node at the origin
   //The geometric data is is the x,y,z position of the wall
   //to which the spine is affixed
   Vector<double> origin_pt(3);
   origin_pt[0] = spine_node_pt->spine_pt()->geom_parameter(0); 
   origin_pt[2] = spine_node_pt->spine_pt()->geom_parameter(2); 

   //set scale norm vector
   double norm = 
    sqrt( (Xsim -  origin_pt[0] ) *  (Xsim -  origin_pt[0] )+
          (Zsim -  origin_pt[2] ) *  (Zsim -  origin_pt[2])  );

   //Set th x coordinate
   spine_node_pt->x(0) = origin_pt[0] + H*W* (Xsim -  origin_pt[0] )/norm;
   
   //Set the value of z
   spine_node_pt->x(2) = origin_pt[2] + H*W* (Zsim -  origin_pt[2] )/norm;
   
  }


//This function wull be called only at the begining for setting the initial spine appropiated for makeing the canyon shape
//We will send the node wich is in the inmovile face (origin of the spines)
 double init_spine_height(SpineNode* spine_node_pt)
  {

   //set distance to the axis
   double daxis = sqrt( (this->Xsim -  spine_node_pt->x(0) ) *(this->Xsim - spine_node_pt->x(0) )+(this->Zsim -  spine_node_pt->x(2) ) *  (this->Zsim -  spine_node_pt->x(2)));
   double init_h = ( daxis - radius() ) ; 
   return init_h;
  
  }

 //========================================================
 /// Flush storage for spines, elements and nodes by emptying the
 /// vectors that store the pointers to them. This is
 /// useful if a particular mesh is only built to generate
 /// a small part of a bigger mesh. Once the elements and
 /// nodes have been created, they are typically copied
 /// into the new mesh and the auxiliary mesh can be
 /// deleted. However, if we simply call the destructor
 /// of the auxiliary mesh, it will also wipe out
 /// the nodes and elements, because it still "thinks"
 /// it's in charge of these...
 //========================================================


  void flush_spine_element_and_node_storage()
  { 
    //Clear vectors of pointers to the nodes and elements
    Node_pt.clear();
    Bulk_element_pt.clear();
    Interface_element_pt.clear();
    Element_pt.clear();
    Spine_pt.clear();
  }




protected:

 /// Vector of pointers to element in the fluid layer
 Vector <FiniteElement *> Bulk_element_pt;

 /// Vector of pointers to interface elements
 Vector<FiniteElement *> Interface_element_pt;

/// Sacled radius of the canyon (must be between 0 and one
 double Radius;

 ///  Helper function to actually build the single-layer spine mesh 
 /// (called from various constructors)
 virtual void build_single_layer_mesh(TimeStepper* time_stepper_pt);
 

//Axis of the cilinder
 double Xsim;
 
 double Zsim;

// Rotation flag

 unsigned Rotation_flag;

};

 

//#endif


//===========================================================================
/// Constructor for spine 3D mesh: Pass number of elements in x-direction, 
/// number of elements in y-direction, number elements in z-direction, 
/// axial length, deep and height of layer, 
/// and pointer to timestepper (defaults to Static timestepper).
///
/// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
/// e.g  SpineElement<CrozierRaviartElement<2>)

/// and a surface layer of corresponding Spine interface elements
/// of type INTERFACE_ELEMENT, e.g.
/// SpineLineNavierStokesInterfaceElement<ELEMENT> for 2D planar
/// problems.
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
MyCanyonMesh<ELEMENT, INTERFACE_ELEMENT>::MyCanyonMesh(
 const unsigned &nx, const unsigned &ny,  const unsigned &nz,
 const double &x_min, const double &x_max, const double &y_min, 
 const double &y_max, const double &z_min, const double &z_max, 
 const double &R, const unsigned &rotation_flag, TimeStepper* time_stepper_pt) :
  SimpleCubicMesh<ELEMENT >(nx,ny,nz,x_min,x_max,y_min,y_max,z_min,z_max,
                            time_stepper_pt)
{
 // We've called the "generic" constructor for the RectangularQuadMesh
 // which doesn't do much...
 // Now set up the parameters that characterise the mesh geometry
 // then call the constructor
 
// Axis of the cilinder
  this->Xsim = 0.0;
 
  this->Zsim = 1.0;
  
  //Set rotation
  this->Rotation_flag = rotation_flag;
  
  change_radius(R);
  
  build_single_layer_mesh(time_stepper_pt);
}

//===========================================================================
/// Helper function that actually builds the single-layer spine mesh
/// based on the parameters set in the various constructors
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
void MyCanyonMesh<ELEMENT, INTERFACE_ELEMENT>::build_single_layer_mesh(
 TimeStepper* time_stepper_pt) 
{
 
// Set rotation
 if(Rotation_flag == 1)  //rotation around y--axis
   for(unsigned i = 0;i< this->nnode();i++)
    { 
      double x_node = this->node_pt(i)->x(0);
      double z_node = this->node_pt(i)->x(2);
      this->node_pt(i)->x(0) = -z_node +1.0;
      this->node_pt(i)->x(2) = x_node;
    }



 //Read out the number of elements in the x-direction
 unsigned n_x = this->Nx;
 unsigned n_y = this->Ny;
 unsigned n_z = this->Nz;



 //Set up the pointers to elements in the bulk
 unsigned nbulk=nelement();
 Bulk_element_pt.reserve(nbulk);
 for(unsigned e=0;e<nbulk;e++)
  {
//   Bulk_element_pt.push_back(this->finite_element_pt(e));
     Bulk_element_pt.push_back(dynamic_cast<ELEMENT*>(this->finite_element_pt(e)));
  }

 //Allocate memory for the spines and fractions along spines
 //---------------------------------------------------------

 //Read out number of linear points in the element
 unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

 //Allocate store for the spines: (different in the case of periodic meshes !!)
   Spine_pt.reserve(((n_p-1)*n_x+1)*((n_p-1)*n_y+1));
  

// NOW WE GO THROUH ALL THE ELEMENTS OF THE LOWER LAYERS AND ATTACHING THE SPINES

// FIRST ELEMENT: Element 0

  

for(unsigned l1=0;l1<n_p;l1++) //y loop over the nodes
{
 for(unsigned l2=0;l2<n_p;l2++)  //x loop over the nodes
  {

   //Node j + i*np
// Vector to the origin node, it will be passed as argument for the a later update of the spine
    Vector<double> origin_node(3);
   origin_node[0] = element_node_pt(0,l2+l1*n_p)->x(0);
   origin_node[1] = element_node_pt(0,l2+l1*n_p)->x(1);
   origin_node[2] = element_node_pt(0,l2+l1*n_p)->x(2);

   double hinit  = init_spine_height( element_node_pt(0,l2+l1*n_p) );
   //Assign the new spine within the cilinder
   Spine* new_spine_pt=new Spine( hinit);
   new_spine_pt->set_geom_parameter(origin_node);
   Spine_pt.push_back(new_spine_pt);

   // Get pointer to node
   SpineNode* nod_pt=element_node_pt(0,l2+l1*n_p); // Element 0; node j + i*n_p
   //Set the pointer to the spine
   nod_pt->spine_pt() = new_spine_pt;
   //Set the fraction
   nod_pt->fraction() = 0.0;
   // Pointer to the mesh that implements the update fct
   nod_pt->spine_mesh_pt() = this; 

  //Loop vertically along the spine
   //Loop over the elements 
   for(unsigned long k=0;k<n_z;k++)
   {
    //Loop over the vertical nodes, apart from the first
    for(unsigned l3=1;l3<n_p;l3++)
     {
      // Get pointer to node
      SpineNode* nod_pt=element_node_pt(k*n_x*n_y,l3*n_p*n_p+l2+l1*n_p);
      //Set the pointer to the spine
      nod_pt->spine_pt() = new_spine_pt;
      //Set the fraction
      nod_pt->fraction()=(double(k)+double(l3)/double(n_p-1))/double(n_z);
      // Pointer to the mesh that implements the update fct
       nod_pt->spine_mesh_pt() = this; 
     }
   }
  }
}



 //LOOP OVER OTHER ELEMENTS IN THE FIRST ROW
 //-----------------------------------------

// The procedure is the same but we have to identify the before defined spines for not defining them two times

for(unsigned j=1;j<n_x;j++) //loop over the elements in the first row
{

  for(unsigned l1=0;l1<n_p;l1++) //y loop over the nodes
  {
   
   for(unsigned l2=1;l2<n_p;l2++)  //x loop over the nodes
    {

     //Node j + i*np
    Vector<double> origin_node(3);
    origin_node[0] = element_node_pt(j,l2+l1*n_p)->x(0);
    origin_node[1] = element_node_pt(j,l2+l1*n_p)->x(1);
    origin_node[2] = element_node_pt(j,l2+l1*n_p)->x(2);
    double hinit  = init_spine_height( element_node_pt(j,l2+l1*n_p) ); 
   //Assign the new spine
    Spine* new_spine_pt=new Spine( hinit);   
    new_spine_pt->set_geom_parameter(origin_node);

     Spine_pt.push_back(new_spine_pt);


     // Get pointer to node
     SpineNode* nod_pt=element_node_pt(j,l2+l1*n_p); // Element j; node l2 + l1*n_p
     //Set the pointer to the spine
     nod_pt->spine_pt() = new_spine_pt;
     //Set the fraction
     nod_pt->fraction() = 0.0;
     // Pointer to the mesh that implements the update fct
     nod_pt->spine_mesh_pt() = this; 

    //Loop vertically along the spine
     //Loop over the elements 
     for(unsigned long k=0;k<n_z;k++)
     {
      //Loop over the vertical nodes, apart from the first
      for(unsigned l3=1;l3<n_p;l3++)
       {
        // Get pointer to node
        SpineNode* nod_pt=element_node_pt(j+k*n_x*n_y,l3*n_p*n_p+l2+l1*n_p);
        //Set the pointer to the spine
        nod_pt->spine_pt() = new_spine_pt;
        //Set the fraction
        nod_pt->fraction()=(double(k)+double(l3)/double(n_p-1))/double(n_z);
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this; 
       }
     }
     
     
    }
  }
}

//REST OF THE ELEMENTS
// Now we loop over the rest of the elements. We will separe the first of each row being al the rest equal



for(unsigned long i=1;i<n_y;i++)
{
//FIRST ELEMENT OF THE ROW
   
 //First line of nodes is copied from the element of the bottom
 
  for(unsigned l1=1;l1<n_p;l1++) //y loop over the nodes
  {
   
   for(unsigned l2=0;l2<n_p;l2++)  //x loop over the nodes
    {

     //Node j + i*np
    Vector<double> origin_node(3); 
    origin_node[0] = element_node_pt(i*n_x,l2+l1*n_p)->x(0);
    origin_node[1] = element_node_pt(i*n_x,l2+l1*n_p)->x(1);
    origin_node[2] = element_node_pt(i*n_x,l2+l1*n_p)->x(2);
     double hinit  = init_spine_height( element_node_pt(i*n_x,l2+l1*n_p) );   
    //Assign the new spine 
    Spine* new_spine_pt= new Spine( hinit);
    new_spine_pt->set_geom_parameter(origin_node);

     Spine_pt.push_back(new_spine_pt);


     // Get pointer to node
     SpineNode* nod_pt=element_node_pt(i*n_x,l2+l1*n_p); // Element i*n_x; node l2 + l1*n_p
     //Set the pointer to the spine
     nod_pt->spine_pt() = new_spine_pt;
     //Set the fraction
     nod_pt->fraction() = 0.0;
     // Pointer to the mesh that implements the update fct
     nod_pt->spine_mesh_pt() = this; 

    //Loop vertically along the spine
     //Loop over the elements 
     for(unsigned long k=0;k<n_z;k++)
     {
      //Loop over the vertical nodes, apart from the first
      for(unsigned l3=1;l3<n_p;l3++)
       {
        // Get pointer to node
        SpineNode* nod_pt=element_node_pt(i*n_x+k*n_x*n_y,l3*n_p*n_p+l2+l1*n_p);
        //Set the pointer to the spine
        nod_pt->spine_pt() = new_spine_pt;
        //Set the fraction
        nod_pt->fraction()=(double(k)+double(l3)/double(n_p-1))/double(n_z);
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this; 
       }
     }
     
     
    }
  }


 


  // REST OF THE ELEMENTS OF THE ROW
 for(unsigned j=1;j<n_x;j++)
  {


//First line of nodes is copied from the element of the bottom
   
   for(unsigned l1=1;l1<n_p;l1++) //y loop over the nodes
    {
     
     for(unsigned l2=1;l2<n_p;l2++)  //x loop over the nodes
      {
 
       //Node j + i*np
        Vector<double> origin_node(3);
        origin_node[0] = element_node_pt(j+i*n_x,l2+l1*n_p)->x(0);
        origin_node[1] = element_node_pt(j+i*n_x,l2+l1*n_p)->x(1);
        origin_node[2] = element_node_pt(j+i*n_x,l2+l1*n_p)->x(2);
         double hinit  = init_spine_height( element_node_pt(j+i*n_x,l2+l1*n_p) );   
        //Assign the new spine with radius as unit length
        Spine* new_spine_pt=new Spine( hinit);
       new_spine_pt->set_geom_parameter(origin_node);


       Spine_pt.push_back(new_spine_pt);


       // Get pointer to node
       SpineNode* nod_pt=element_node_pt(j+i*n_x,l2+l1*n_p); // Element j + i*n_x; node l2 + l1*n_p
       //Set the pointer to the spine
       nod_pt->spine_pt() = new_spine_pt;
       //Set the fraction
       nod_pt->fraction() = 0.0;
       // Pointer to the mesh that implements the update fct
       nod_pt->spine_mesh_pt() = this; 

      //Loop vertically along the spine
       //Loop over the elements 
       for(unsigned long k=0;k<n_z;k++)
       {
        //Loop over the vertical nodes, apart from the first
        for(unsigned l3=1;l3<n_p;l3++)
         {
          // Get pointer to node
          SpineNode* nod_pt=element_node_pt(j+i*n_x+k*n_x*n_y,l3*n_p*n_p+l2+l1*n_p);
          //Set the pointer to the spine
          nod_pt->spine_pt() = new_spine_pt;
          //Set the fraction
          nod_pt->fraction()=(double(k)+double(l3)/double(n_p-1))/double(n_z);
          // Pointer to the mesh that implements the update fct
          nod_pt->spine_mesh_pt() = this; 
         }
       }
       
       
      }
    }

  }


}






 //Assign the 2D Line elements
 //---------------------------

 //Get the present number of elements
 unsigned long element_count = Element_pt.size();

 //Loop over the elements on the plane
 for(unsigned l1=0;l1<n_y;l1++)
   for(unsigned l2=0;l2<n_x;l2++)
 {
  {
  //Construct a new 2D surface element on the face on which the local
  //coordinate 2 is fixed at its max. value (1)
   FiniteElement *interface_element_element_pt =
    new INTERFACE_ELEMENT(finite_element_pt(n_x*n_y*(n_z-1)+l2+l1*n_x),3);


   //Push it back onto the stack
   Element_pt.push_back(interface_element_element_pt); 

   //Push it back onto the stack of interface elements
   Interface_element_pt.push_back(interface_element_element_pt);

   element_count++;
   }
 }


//Update the node positions in the mesh
  this->node_update();

}



#endif
