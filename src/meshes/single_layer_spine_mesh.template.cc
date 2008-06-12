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
#ifndef OOMPH_SINGLE_LAYER_SPINE_MESH_TEMPLATE_CC
#define OOMPH_SINGLE_LAYER_SPINE_MESH_TEMPLATE_CC

#include "single_layer_spine_mesh.template.h"
#include "rectangular_quadmesh.template.cc"


namespace oomph
{



//===========================================================================
/// Constructor for spine 2D mesh: Pass number of elements in x-direction, 
/// number of elements in y-direction, axial length and height of layer, 
/// and pointer to timestepper (defaults to Static timestepper).
///
/// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
/// e.g  SpineElement<QCrouzeixRaviartElement<2>)
/// and a surface layer of corresponding Spine interface elements
/// of type INTERFACE_ELEMENT, e.g.
/// SpineLineFluidInterfaceElement<ELEMENT> for 2D planar
/// problems.
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
SingleLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::SingleLayerSpineMesh(
 const unsigned &nx, const unsigned &ny,
 const double &lx, const double &h, TimeStepper* time_stepper_pt) :
 RectangularQuadMesh<ELEMENT >(nx,ny,0.0,lx,0.0,h,false,false,
                               time_stepper_pt)
{
 // We've called the "generic" constructor for the RectangularQuadMesh
 // which doesn't do much...
 
 // Now build the mesh: 
 build_single_layer_mesh(time_stepper_pt);
}



//===========================================================================
/// Constuctor for spine 2D mesh: Pass number of elements in x-direction, 
/// number of elements in y-direction, axial length and height of layer, 
/// a boolean flag to make the mesh periodic in the x-direction,
/// and pointer to timestepper (defaults to Static timestepper).
///
/// The mesh contains a layer of elements (of type ELEMENT)
/// and a surface layer of corresponding Spine interface elements (of type
/// SpineLineFluidInterfaceElement<ELEMENT>).
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
SingleLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::SingleLayerSpineMesh(
 const unsigned &nx, const unsigned &ny,
 const double &lx, const double &h,  const bool& periodic_in_x,
 TimeStepper* time_stepper_pt) :
 RectangularQuadMesh<ELEMENT >(nx,ny,0.0,lx,0.0,h,periodic_in_x,false,
                               time_stepper_pt)
{
 // We've called the "generic" constructor for the RectangularQuadMesh
 // which doesn't do much...
 
 // Now build the mesh: 
 build_single_layer_mesh(time_stepper_pt);
}





//===========================================================================
/// Helper function that actually builds the single-layer spine mesh
/// based on the parameters set in the various constructors
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
void SingleLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::build_single_layer_mesh(
 TimeStepper* time_stepper_pt) 
{

 // Build the underlying quad mesh: 
 RectangularQuadMesh<ELEMENT >::build_mesh(time_stepper_pt);
 //Read out the number of elements in the x-direction
 unsigned n_x = this->Nx;
 unsigned n_y = this->Ny;

 //Set up the pointers to elements in the bulk
 unsigned nbulk=nelement();
 Bulk_element_pt.reserve(nbulk);
 for(unsigned e=0;e<nbulk;e++)
  {
   Bulk_element_pt.push_back(this->finite_element_pt(e));
  }

 //Allocate memory for the spines and fractions along spines
 //---------------------------------------------------------

 //Read out number of linear points in the element
 unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

 //Allocate store for the spines:
 if (this->Xperiodic)
  {
   Spine_pt.reserve((n_p-1)*n_x);
  }
 else
  {
   Spine_pt.reserve((n_p-1)*n_x+1);
  }


 //FIRST SPINE
 //-----------

 //Element 0
 //Node 0
 //Assign the new spine with unit length
 Spine* new_spine_pt=new Spine(1.0);
 Spine_pt.push_back(new_spine_pt);


 // Get pointer to node
 SpineNode* nod_pt=element_node_pt(0,0);
 //Set the pointer to the spine
 nod_pt->spine_pt() = new_spine_pt;
 //Set the fraction
 nod_pt->fraction() = 0.0;
 // Pointer to the mesh that implements the update fct
 nod_pt->spine_mesh_pt() = this; 

 //Loop vertically along the spine
 //Loop over the elements 
 for(unsigned long i=0;i<n_y;i++)
  {
   //Loop over the vertical nodes, apart from the first
   for(unsigned l1=1;l1<n_p;l1++)
    {
     // Get pointer to node
     SpineNode* nod_pt=element_node_pt(i*n_x,l1*n_p);
     //Set the pointer to the spine
     nod_pt->spine_pt() = new_spine_pt;
     //Set the fraction
     nod_pt->fraction()=(double(i)+double(l1)/double(n_p-1))/double(n_y);
     // Pointer to the mesh that implements the update fct
     nod_pt->spine_mesh_pt() = this; 
    }
  }


 //LOOP OVER OTHER SPINES
 //----------------------

 //Now loop over the elements horizontally
 for(unsigned long j=0;j<n_x;j++)
  {
   //Loop over the nodes in the elements horizontally, ignoring 
   //the first column

   // Last spine needs special treatment in x-periodic meshes:
   unsigned n_pmax=n_p;
   if ((this->Xperiodic)&&(j==n_x-1)) n_pmax=n_p-1;
    
   for(unsigned l2=1;l2<n_pmax;l2++)
    {
     //Assign the new spine with unit height
     new_spine_pt=new Spine(1.0);
     Spine_pt.push_back(new_spine_pt);

     // Get the node
     SpineNode* nod_pt=element_node_pt(j,l2);
     //Set the pointer to spine
     nod_pt->spine_pt() = new_spine_pt;
     //Set the fraction
     nod_pt->fraction() = 0.0;
     // Pointer to the mesh that implements the update fct
     nod_pt->spine_mesh_pt() = this; 

     //Loop vertically along the spine
     //Loop over the elements 
     for(unsigned long i=0;i<n_y;i++)
      {
       //Loop over the vertical nodes, apart from the first
       for(unsigned l1=1;l1<n_p;l1++)
        {
         // Get the node
         SpineNode* nod_pt=element_node_pt(i*n_x+j,l1*n_p+l2);
         //Set the pointer to the spine
         nod_pt->spine_pt() = new_spine_pt;
         //Set the fraction
         nod_pt->fraction()=(double(i)+double(l1)/double(n_p-1))/double(n_y);
         // Pointer to the mesh that implements the update fct
         nod_pt->spine_mesh_pt() = this; 
        }  
      }
    }
  }


 // Last spine needs special treatment for periodic meshes
 // because it's the same as the first one...
 if (this->Xperiodic)
  {
   // Last spine is the same as first one...
   Spine* final_spine_pt=Spine_pt[0];  

   // Get the node
   SpineNode* nod_pt=element_node_pt((n_x-1),(n_p-1));

   //Set the pointer for the first node
   nod_pt->spine_pt() = final_spine_pt;
   //Set the fraction to be the same as for the nodes on the first row
   nod_pt->fraction() = element_node_pt(0,0)->fraction();
   // Pointer to the mesh that implements the update fct
   nod_pt->spine_mesh_pt() = element_node_pt(0,0)->spine_mesh_pt();
   
   //Now loop vertically along the spine
   for(unsigned i=0;i<n_y;i++)
    {
     //Loop over the vertical nodes, apart from the first
     for(unsigned l1=1;l1<n_p;l1++)
      {
       // Get the node
       SpineNode* nod_pt=element_node_pt(i*n_x+(n_x-1),l1*n_p+(n_p-1));

       //Set the pointer to the spine
       nod_pt->spine_pt() = final_spine_pt;
       //Set the fraction to be the same as in first row
       nod_pt->fraction() = element_node_pt(i*n_x,l1*n_p)->fraction();
       // Pointer to the mesh that implements the update fct
       nod_pt->spine_mesh_pt() = 
        element_node_pt(i*n_x,l1*n_p)->spine_mesh_pt();
      }
    }
   
  }
 
 //Assign the 1D Line elements
 //---------------------------

 //Get the present number of elements
 unsigned long element_count = Element_pt.size();

 //Loop over the horizontal elements
 for(unsigned i=0;i<n_x;i++)
  {
  //Construct a new 1D line element on the face on which the local
  //coordinate 1 is fixed at its max. value (1) --- This is face 2
   FiniteElement *interface_element_element_pt =
    new INTERFACE_ELEMENT(finite_element_pt(n_x*(n_y-1)+i),2);

   //Push it back onto the stack
   Element_pt.push_back(interface_element_element_pt); 

   //Push it back onto the stack of interface elements
   Interface_element_pt.push_back(interface_element_element_pt);

   element_count++;
  }

}


//======================================================================
/// Reorder the elements, so we loop over them vertically first
/// (advantageous in "wide" domains if a frontal solver is used).
//======================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
void SingleLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::element_reorder()
{

 unsigned n_x = this->Nx;
 unsigned n_y = this->Ny;
 //Find out how many elements there are
 unsigned long Nelement = nelement();
 //Find out how many fluid elements there are
 unsigned long Nfluid = n_x*n_y;
 //Create a dummy array of elements
 Vector<FiniteElement *> dummy;

 //Loop over the elements in horizontal order
 for(unsigned long j=0;j<n_x;j++)
  {
   //Loop over the elements in lower layer vertically
   for(unsigned long i=0;i<n_y;i++)
    {
     //Push back onto the new stack
     dummy.push_back(finite_element_pt(n_x*i + j));
    }

   //Push back the line element onto the stack
   dummy.push_back(finite_element_pt(Nfluid+j));
  }

 //Now copy the reordered elements into the element_pt
 for(unsigned long e=0;e<Nelement;e++)
  {
   Element_pt[e] = dummy[e];
  }
}

}
#endif
