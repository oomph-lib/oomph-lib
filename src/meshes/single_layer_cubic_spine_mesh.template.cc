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
#ifndef OOMPH_SINGLE_LAYER_CUBIC_SPINE_MESH_TEMPLATE_CC
#define OOMPH_SINGLE_LAYER_CUBIC_SPINE_MESH_TEMPLATE_CC

#include "single_layer_cubic_spine_mesh.template.h"
#include "simple_cubic_mesh.template.cc"

namespace oomph
{

//===========================================================================
/// Constructor for spine 3D mesh: Pass number of elements in x-direction, 
/// number of elements in y-direction, number elements in z-direction, 
/// length, width and height of layer, 
/// and pointer to timestepper (defaults to Static timestepper).
///
/// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
/// e.g  SpineElement<QCrouzeixRaviartElement<3>)
/// and a surface layer of corresponding Spine interface elements
/// of type INTERFACE_ELEMENT, e.g.
/// SpineSurfaceFluidInterfaceElement<ELEMENT> for 2D planar
/// problems.
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
SingleLayerCubicSpineMesh<ELEMENT, INTERFACE_ELEMENT>::
SingleLayerCubicSpineMesh(
 const unsigned &nx, const unsigned &ny,  const unsigned &nz, 
 const double &lx, const double &ly, const double &h, 
 TimeStepper* time_stepper_pt) :
 SimpleCubicMesh<ELEMENT >(nx,ny,nz,lx,ly,h,time_stepper_pt)
{
 // Mesh can only be built with 3D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(3);

 // Now build the single layer mesh on top of the existing mesh
 build_single_layer_mesh(time_stepper_pt);
}

//===========================================================================
/// Helper function that actually builds the single-layer spine mesh
/// based on the parameters set in the various constructors
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
void SingleLayerCubicSpineMesh<ELEMENT, INTERFACE_ELEMENT>::
build_single_layer_mesh(
 TimeStepper* time_stepper_pt) 
{
 // Mesh can only be built with 3D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(3);

 //Read out the number of elements in the x-direction
 unsigned n_x = this->Nx;
 //Read out the number of elements in the y-direction
 unsigned n_y = this->Ny;
 //Read out the number of elements in the z-direction
 unsigned n_z = this->Nz;

 //Set up the pointers to elements in the bulk
 unsigned nbulk=nelement();
 Bulk_element_pt.reserve(nbulk);
 for(unsigned e=0;e<nbulk;e++)
  {
   Bulk_element_pt.push_back(this->finite_element_pt(e));
  }

 //Allocate memory for the spines and fractions along spines
 
 //Read out number of linear points in the element
 unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();
 
 //Allocate store for the spines: (different in the case of periodic meshes !!)
 Spine_pt.reserve(((n_p-1)*n_x+1)*((n_p-1)*n_y+1));
 
 //Now we loop over all the elements and attach the spines

 // FIRST ELEMENT: Element 0
 //Loop over the nodes on the base of the element
 for(unsigned l1=0;l1<n_p;l1++) //y loop over the nodes
  {
   for(unsigned l2=0;l2<n_p;l2++)  //x loop over the nodes
    {
     //Assign the new spine with length h
     Spine* new_spine_pt=new Spine(1.0);
     Spine_pt.push_back(new_spine_pt);
     
     // Get pointer to node
     SpineNode* nod_pt = element_node_pt(0, l2 + l1*n_p); 
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
 
// The procedure is the same but we have to identify the 
// before defined spines for not defining them two times
 
 for(unsigned j=1;j<n_x;j++) //loop over the elements in the first row
  {
   for(unsigned l1=0;l1<n_p;l1++) //y loop over the nodes
    {
     // First we copy the last row of nodes into the first row of the new element (and extend to the third dimension)
   for(unsigned l2=1;l2<n_p;l2++)  //x loop over the nodes
    {
     //Node j + i*np
     //Assign the new spine with unit length
     Spine* new_spine_pt=new Spine(1.0);
     Spine_pt.push_back(new_spine_pt);
     
     // Get pointer to node
     SpineNode* nod_pt=element_node_pt(j,l2+l1*n_p); 

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
     //Assign the new spine with unit length
     Spine* new_spine_pt=new Spine(1.0);
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
       //Assign the new spine with unit length
       Spine* new_spine_pt=new Spine(1.0);
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
  //coordinate 2 is fixed at its max. value (1) -- Face 3
   FiniteElement *interface_element_element_pt =
    new INTERFACE_ELEMENT(finite_element_pt(n_x*n_y*(n_z-1)+l2+l1*n_x),3);

   //Push it back onto the stack
   Element_pt.push_back(interface_element_element_pt); 

   //Push it back onto the stack of interface elements
   Interface_element_pt.push_back(interface_element_element_pt);

   //Now assign the spines to the elements
   element_count++;
   }
  }

}

}
#endif
