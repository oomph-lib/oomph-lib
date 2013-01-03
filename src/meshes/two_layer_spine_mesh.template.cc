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
#ifndef OOMPH_TWO_LAYER_SPINE_MESH_TEMPLATE_CC
#define OOMPH_TWO_LAYER_SPINE_MESH_TEMPLATE_CC

#include "two_layer_spine_mesh.template.h"
#include "rectangular_quadmesh.template.cc"


namespace oomph
{


//===========================================================================
/// Constuctor for spine 2D mesh: Pass number of elements in x-direction, 
/// number of elements in y-direction in bottom and top layer, respectively,
/// axial length and height of top and bottom layers, and pointer to 
/// timestepper (defaults to Static timestepper).
///
/// The mesh contains two layers of elements (of type ELEMENT;
/// e.g  SpineElement<QCrouzeixRaviartElement<2>)
/// and an interfacial layer of corresponding Spine interface elements 
/// of type INTERFACE_ELEMENT, e.g.
/// SpineLineFluidInterfaceElement<ELEMENT> for 2D planar
/// problems.
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::TwoLayerSpineMesh(
 const unsigned &nx, const unsigned &ny1, const unsigned &ny2,
 const double &lx, const double &h1, const double &h2,
 TimeStepper* time_stepper_pt) :
 RectangularQuadMesh<ELEMENT >(nx,ny1+ny2,0.0,lx,0.0,h1+h2,false,false,
                               time_stepper_pt)
{
 // Mesh can only be built with 2D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

 // We've called the "generic" constructor for the RectangularQuadMesh
 // which doesn't do much...
 // Now set up the parameters that characterise the mesh geometry
 // then call the build function

 // Number of elements in bottom and top layers
 Ny1 = ny1;
 Ny2 = ny2; 
 
 // Set height of upper and lower layers
 H1 = h1;
 H2 = h2;
 
 // Now build the mesh: 
 build_two_layer_mesh(time_stepper_pt);

}




//===========================================================================
/// Constuctor for spine 2D mesh: Pass number of elements in x-direction, 
/// number of elements in y-direction in bottom and top layer, respectively,
/// axial length and height of top and bottom layers, a boolean
/// flag to make the mesh periodic in the x-direction, and pointer to 
/// timestepper (defaults to Static timestepper).
///
/// The mesh contains two layers of elements (of type ELEMENT;
/// e.g  SpineElement<QCrouzeixRaviartElement<2>)
/// and an interfacial layer of corresponding Spine interface elements 
/// of type INTERFACE_ELEMENT, e.g.
/// SpineLineFluidInterfaceElement<ELEMENT> for 2D planar
/// problems.
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::TwoLayerSpineMesh(
 const unsigned &nx, const unsigned &ny1, const unsigned &ny2,
 const double &lx, const double &h1, const double &h2,
 const bool& periodic_in_x, TimeStepper* time_stepper_pt) :
 RectangularQuadMesh<ELEMENT >(nx,ny1+ny2,0.0,lx,0.0,h1+h2,periodic_in_x,
                               false,time_stepper_pt)
{ 
 // Mesh can only be built with 2D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

 // We've called the "generic" constructor for the RectangularQuadMesh
 // which doesn't do much...
 // Now set up the parameters that characterise the mesh geometry
 // then call the constructor

 // Number of elements in bottom and top layers
 Ny1 = ny1;
 Ny2 = ny2; 
 
 // Set height of upper and lower layers
 H1 = h1;
 H2 = h2;
 
 // Now build the mesh: 
 build_two_layer_mesh(time_stepper_pt);

}




//===========================================================================
/// Constuctor for spine 2D mesh: Pass number of elements in x-direction, 
/// number of elements in y-direction in bottom and top layer, respectively,
/// axial length and height of top and bottom layers, a boolean
/// flag to make the mesh periodic in the x-direction, a boolean flag to
/// specify whether or not to call the "build_two_layer_mesh" function,
/// and pointer to timestepper (defaults to Static timestepper).
///
/// The mesh contains two layers of elements (of type ELEMENT;
/// e.g  SpineElement<QCrouzeixRaviartElement<2>)
/// and an interfacial layer of corresponding Spine interface elements 
/// of type INTERFACE_ELEMENT, e.g.
/// SpineLineFluidInterfaceElement<ELEMENT> for 2D planar
/// problems.
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::TwoLayerSpineMesh(
 const unsigned &nx, const unsigned &ny1, const unsigned &ny2,
 const double &lx, const double &h1, const double &h2,
 const bool& periodic_in_x, const bool& build_mesh,
 TimeStepper* time_stepper_pt) :
 RectangularQuadMesh<ELEMENT >(nx,ny1+ny2,0.0,lx,0.0,h1+h2,periodic_in_x,
                               false,time_stepper_pt)
{
 // Mesh can only be built with 2D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

 // We've called the "generic" constructor for the RectangularQuadMesh
 // which doesn't do much...
 // Now set up the parameters that characterise the mesh geometry
 // then call the constructor

 // Number of elements in bottom and top layers
 Ny1 = ny1;
 Ny2 = ny2; 
 
 // Set height of upper and lower layers
 H1 = h1;
 H2 = h2;
 
 // Only build the mesh here if build_mesh=true
 // This is useful when calling this constructor from a derived class
 // (such as Axisym2x6TwoLayerSpineMesh) where the mesh building
 // needs to be called from *its* constructor and this constructor is
 // only used to pass arguments to the RectangularQuadMesh constructor.
 if(build_mesh) { build_two_layer_mesh(time_stepper_pt); }

}


//==================================================================
/// \short The spacing function for the x co-ordinate, which is the
/// same as the default function.
//==================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
double TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
x_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
{
 //Calculate the values of equal increments in nodal values
 double xstep = (this->Xmax-this->Xmin)/((this->Np-1)*this->Nx);
 //Return the appropriate value
 return (this->Xmin + xstep*((this->Np-1)*xelement + xnode));
}

//==================================================================
/// \short The spacing function for the y co-ordinates, which splits
/// the region into two regions (1 and 2), according to the 
/// heights H1 and H2, with Ny1 and Ny2 elements respectively.
template<class ELEMENT, class INTERFACE_ELEMENT>
double TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
y_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
{
 //Set up some spacing parameters
 //The lower region a starts at Ymin
 double Ymin = RectangularQuadMesh<ELEMENT >::Ymin;
 //The upper region a ends at Ymax
 double Ymax = RectangularQuadMesh<ELEMENT >::Ymax;
 //Number of nodes per element
 unsigned n_p = RectangularQuadMesh<ELEMENT >::Np;
 //The lower region starts at Ymin
 double y1init = Ymin;
 //The upper region starts at H1 - Ymin
 double y2init = H1 - Ymin;
 //Calculate the space between each node in each region,
 //Assumming uniform spacing
 //Lower region has a length (H1-Ymin)
 double y1step = (H1-Ymin)/((n_p-1)*Ny1);
 //Upper region has a length (Ymax-H1)
 double y2step = (Ymax-H1)/((n_p-1)*Ny2);
 
 //Now return the actual node position, it's different in the two
 //regions, of course
 if(yelement < Ny1) 
  {
   return (y1init + y1step*((n_p-1)*yelement + ynode));
  }
 else
  {
   return (y2init + y2step*((n_p-1)*(yelement-Ny1) + ynode));
  }
}

//===========================================================================
/// Helper function that actually builds the two-layer spine mesh
/// based on the parameters set in the various constructors
//===========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
void TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::build_two_layer_mesh(
 TimeStepper* time_stepper_pt) 
{
 // Mesh can only be built with 2D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

 // Build the underlying quad mesh: 
 RectangularQuadMesh<ELEMENT >::build_mesh(time_stepper_pt);

 //Set up the pointers to elements in the upper and lower fluid
 //Calculate numbers of nodes in upper and lower regions
 unsigned long n_lower = this->Nx*Ny1;
 unsigned long n_upper = this->Nx*Ny2;
 //Loop over lower elements and push back
 Lower_layer_element_pt.reserve(n_lower);
 for(unsigned e=0;e<n_lower;e++)
  {
   Lower_layer_element_pt.push_back(this->finite_element_pt(e));
  }
 //Loop over upper elements and push back
 Upper_layer_element_pt.reserve(n_upper);
 for(unsigned e=n_lower;e<(n_lower+n_upper);e++)
  {
   Upper_layer_element_pt.push_back(this->finite_element_pt(e));
  }
 


 //Allocate memory for the spines and fractions along spines
 //---------------------------------------------------------

 //Read out number of linear points in the element
 unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

 //Allocate store for the spines:
 if (this->Xperiodic)
  {
   Spine_pt.reserve((n_p-1)*this->Nx);
  }
 else
  {
   Spine_pt.reserve((n_p-1)*this->Nx+1);
  }


 //FIRST SPINE
 //-----------

 //Element 0
 //Node 0
 //Assign the new spine of height of the lower layer
 Spine* new_spine_pt=new Spine(H1);
 Spine_pt.push_back(new_spine_pt);

 // Get pointer to node
 SpineNode* nod_pt=element_node_pt(0,0);

 //Set the pointer to the spine
 nod_pt->spine_pt() = new_spine_pt;
 //Set the fraction
 nod_pt->fraction() = 0.0;
 // Pointer to the mesh that implements the update fct
 nod_pt->spine_mesh_pt() = this; 
 // ID of update function within this mesh: 0 = lower; 1 = upper
 nod_pt->node_update_fct_id() = 0;

 //Loop vertically along the spine
 //Loop over the elements in fluid 1
 for(unsigned long i=0;i<Ny1;i++)
  {
   //Loop over the vertical nodes, apart from the first
   for(unsigned l1=1;l1<n_p;l1++)
    {
     // Get pointer to node
     SpineNode* nod_pt=element_node_pt(i*this->Nx,l1*n_p);
     //Set the pointer to the spine
     nod_pt->spine_pt() = new_spine_pt;
     //Set the fraction
     nod_pt->fraction() = (nod_pt->x(1)-this->Ymin)/(H1);
     // Pointer to the mesh that implements the update fct
     nod_pt->spine_mesh_pt() = this; 
     // ID of update function within this mesh: 0 = lower; 1 = upper
     nod_pt->node_update_fct_id() = 0;
    }
  }

 //Loop over the elements in fluid 2
 for(unsigned long i=0;i<Ny2;i++)
  {
   //Loop over vertical nodes, apart from the first
   for(unsigned l1=1;l1<n_p;l1++)
    {
     // Get pointer to node
     SpineNode* nod_pt=element_node_pt((Ny1+i)*this->Nx,l1*n_p);

     //Set the pointer to the spine
     nod_pt->spine_pt() = new_spine_pt;
     //Set the fraction
     nod_pt->fraction() =(nod_pt->x(1)-(this->Ymin+H1))/(H2);
     // Pointer to the mesh that implements the update fct
     nod_pt->spine_mesh_pt() = this; 
     // ID of update function within this mesh: 0 = lower; 1 = upper
     nod_pt->node_update_fct_id() = 1;
    }
  }


 //LOOP OVER OTHER SPINES
 //----------------------

 //Now loop over the elements horizontally
 for(unsigned long j=0;j<this->Nx;j++)
  {
   //Loop over the nodes in the elements horizontally, ignoring 
   //the first column

   // Last spine needs special treatment in x-periodic meshes:
   unsigned n_pmax=n_p;
   if ((this->Xperiodic)&&(j==this->Nx-1)) n_pmax=n_p-1;
    
   for(unsigned l2=1;l2<n_pmax;l2++)
    {
     //Assign the new spine with length the height of the lower layer
     new_spine_pt=new Spine(H1);
     Spine_pt.push_back(new_spine_pt);

     // Get pointer to node
     SpineNode* nod_pt=element_node_pt(j,l2);

     //Set the pointer to the spine
     nod_pt->spine_pt() = new_spine_pt;
     //Set the fraction
     nod_pt->fraction() = 0.0;
     // Pointer to the mesh that implements the update fct
     nod_pt->spine_mesh_pt() = this; 
     // ID of update function within this mesh: 0 = lower; 1 = upper
     nod_pt->node_update_fct_id() = 0;

     //Loop vertically along the spine
     //Loop over the elements in fluid 1
     for(unsigned long i=0;i<Ny1;i++)
      {
       //Loop over the vertical nodes, apart from the first
       for(unsigned l1=1;l1<n_p;l1++)
        {
         // Get pointer to node
         SpineNode* nod_pt=element_node_pt(i*this->Nx+j,l1*n_p+l2);
         //Set the pointer to the spine
         nod_pt->spine_pt() = new_spine_pt;
         //Set the fraction
         nod_pt->fraction() = (nod_pt->x(1)-this->Ymin)/H1;
         // Pointer to the mesh that implements the update fct
         nod_pt->spine_mesh_pt() = this; 
         // ID of update function within this mesh: 0 = lower; 1 = upper
         nod_pt->node_update_fct_id() = 0;
        }  
      }

     //Loop over the elements in fluid 2
     for(unsigned long i=0;i<Ny2;i++)
      {
       //Loop over vertical nodes, apart from the first
       for(unsigned l1=1;l1<n_p;l1++)
        {
         // Get pointer to node
         SpineNode* nod_pt=element_node_pt((Ny1+i)*this->Nx+j,l1*n_p+l2);

         //Set the pointer to the spine
         nod_pt->spine_pt() = new_spine_pt;
         //Set the fraction
         nod_pt->fraction() = (nod_pt->x(1)-(this->Ymin+H1))/H2;
         // Pointer to the mesh that implements the update fct
         nod_pt->spine_mesh_pt() = this; 
         // ID of update function within this mesh: 0 = lower; 1 = upper
         nod_pt->node_update_fct_id() = 1;
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

   // Get pointer to node
   SpineNode* nod_pt=element_node_pt((this->Nx-1),(n_p-1));

   //Set the pointer to the spine
   nod_pt->spine_pt() = final_spine_pt;
   //Set the fraction to be the same as for the nodes on the first row
   nod_pt->fraction() = element_node_pt(0,0)->fraction();
   // Pointer to the mesh that implements the update fct
   nod_pt->spine_mesh_pt() = element_node_pt(0,0)->spine_mesh_pt();
   // ID of update function within this mesh: 0 = lower; 1 = upper
   nod_pt->node_update_fct_id() = element_node_pt(0,0)->node_update_fct_id();
   
   //Now loop vertically along the spine
   for(unsigned i=0;i<(Ny1+Ny2);i++)
    {
     //Loop over the vertical nodes, apart from the first
     for(unsigned l1=1;l1<n_p;l1++)
      {
       // Get pointer to node
       SpineNode* nod_pt = 
        element_node_pt(i*this->Nx+(this->Nx-1),l1*n_p+(n_p-1));

       //Set the pointer to the spine
       nod_pt->spine_pt() = final_spine_pt;
       //Set the fraction to be the same as in first row
       nod_pt->fraction() = element_node_pt(i*this->Nx,l1*n_p)->fraction();
       // ID of update function within this mesh: 0 = lower; 1 = upper
       nod_pt->node_update_fct_id() = 
        element_node_pt(i*this->Nx,l1*n_p)->node_update_fct_id();
       // Pointer to the mesh that implements the update fct
       nod_pt->spine_mesh_pt() = element_node_pt(i*this->Nx,l1*n_p)
        ->spine_mesh_pt();
      }
    }
  }
 

 //Assign the 1D Line elements
 //---------------------------

 //Get the present number of elements
 unsigned long element_count = Element_pt.size();

 //Loop over the horizontal elements
 for(unsigned i=0;i<this->Nx;i++)
  {
  //Construct a new 1D line element on the face on which the local
  //coordinate 1 is fixed at its max. value (1) -- Face 2
   FiniteElement *interface_element_element_pt =
    new INTERFACE_ELEMENT(finite_element_pt(this->Nx*(Ny1-1)+i),2);

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
void TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::element_reorder()
{

 unsigned Nx = this->Nx;
 //Find out how many elements there are
 unsigned long Nelement = nelement();
 //Find out how many fluid elements there are
 unsigned long Nfluid = Nx*(Ny1+Ny2);
 //Create a dummy array of elements
 Vector<FiniteElement *> dummy;

 //Loop over the elements in horizontal order
 for(unsigned long j=0;j<Nx;j++)
  {
   //Loop over the elements in lower layer vertically
   for(unsigned long i=0;i<Ny1;i++)
    {
     //Push back onto the new stack
     dummy.push_back(finite_element_pt(Nx*i + j));
    }

   //Push back the line element onto the stack
   dummy.push_back(finite_element_pt(Nfluid+j));

   //Loop over the elements in upper layer vertically
   for(unsigned long i=Ny1;i<(Ny1+Ny2);i++)
    {
     //Push back onto the new stack
     dummy.push_back(finite_element_pt(Nx*i + j));
    }
  }

 //Now copy the reordered elements into the element_pt
 for(unsigned long e=0;e<Nelement;e++)
  {
   Element_pt[e] = dummy[e];
  }

}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// //=========================================================================
// /// \short Constructor for a non-uniform two layer spine mesh, with 
// /// element layout in the lower fluid reflected in the upper. Three
// /// distinct y regions need numbers of element specified and two 
// /// x regions.
// //=========================================================================
// template <class ELEMENT, class INTERFACE_ELEMENT >
// Axisym2x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
// Axisym2x6TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb, 
//                            const unsigned &nya1, const unsigned &nyb1,
//                            const unsigned &nyc1, const unsigned &nya2, 
//                            const unsigned &nyb2, const unsigned &nyc2,
//                            const double &lx, 
//                            const double &h1, const double &h2,
//                            TimeStepper* time_stepper_pt) :
//  TwoLayerSpineMesh<ELEMENT, INTERFACE_ELEMENT >()
// {
//  // We've called the "generic" constructor for the RectangularQuadMesh
//  // which doesn't do much...
//  // Now set up the parameters that characterise the mesh geometry
//  // then call the constructor

//  Nxa = nxa; Nxb = nxb;
//  Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1;
//  Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2;
//  Xfraction = 0.8;
//  Yfraction1 = 0.2;
//  Yfraction2 = 0.8;
 
//  // Check validity of Xfraction, Yfraction1 and Yfraction2.
//  if (Yfraction1 < 0.0 || Yfraction1>Yfraction2)
//   {
//    throw OomphLibError("Invalid Yfraction1",
//                        "AxiSym2x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//  if (Yfraction2 < Yfraction1 || Yfraction2>1.0)
//   {
//    throw OomphLibError("Invalid Yfraction2",
//                        "AxiSym2x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//  if (Xfraction<0.0 || Xfraction>1.0)
//   {
//    throw OomphLibError("Invalid xfraction",
//                        "AxiSym2x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }


//  // Number of elements in x direction
//  RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb; 
 
//  // Number of elements in bottom and top layers
//  this->Ny1 = Nya1+Nyb1+Nyc1;
//  this->Ny2 = Nya2+Nyb2+Nyc2;
 
//  // Number of elements in y direction
//  RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
//  // Min. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
//  // Max. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
//  // Min. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
//  // Max. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymax = h1+h2;
 
//  // Periodic?
//  RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
//  // Set height of upper and lower layers
//  this->H1 = h1;
//  this->H2 = h2;
 
//  // Now build the mesh: 
//  this->build_two_layer_mesh(time_stepper_pt);
 
// }


//=========================================================================
/// \short Constructor for a non-uniform two layer spine mesh, with 
/// element layout in the lower fluid reflected in the upper. Three
/// distinct y regions need numbers of element specified and two 
/// x regions. The fractions of these regions are specified in this 
/// constructor.
// PATRICKFLAG this is the one I've been editing...
//=========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT >
Axisym2x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
Axisym2x6TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb, 
                           const double &x_frac, 
                           const unsigned &nya1, const unsigned &nyb1, 
                           const unsigned &nyc1, const unsigned &nya2,
                           const unsigned &nyb2, const unsigned &nyc2, 
                           const double &y1_frac1, const double &y1_frac2,
                           const double &y2_frac1, const double &y2_frac2,
                           const double &lx, 
                           const double &h1, const double &h2,
                           TimeStepper* time_stepper_pt) :
 TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>(nxa+nxb,nya1+nyb1+nyc1,
                                              nya2+nyb2+nyc2,lx,h1,h2,
                                              false,false,time_stepper_pt)
{
 // Mesh can only be built with 2D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

 // We've called the "generic" constructor for the TwoLayerSpineMesh and
 // set the "build_mesh" flag to false. This is done so that we do not
 // call "build_two_layer_mesh(...)" prematurely. We now set up the
 // parameters that characterise this particular mesh's geometry before
 // calling "build_two_layer_mesh(...)".
 
 Nxa = nxa; Nxb = nxb;
 Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1;
 Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2;
 Xfraction = x_frac;
 Y1fraction1 = y1_frac1; Y1fraction2 = y1_frac2;
 Y2fraction1 = y2_frac1; Y2fraction2 = y2_frac2;

 // Check validaty of Xfraction
 if (Xfraction<0.0 || Xfraction>1.0)
  {
   throw OomphLibError("Invalid Xfraction",
                       "Axisym2x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Check validaty of Y1fraction1 and Y2fraction1
 if (Y1fraction1<0.0 || Y1fraction1>Y1fraction2 || Y1fraction1>1.0)
  {
   throw OomphLibError("Invalid Y1fraction1",
                       "Axisym2x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Y2fraction1<0.0 || Y2fraction1>Y2fraction2 || Y2fraction1>1.0)
  {
   throw OomphLibError("Invalid Y2fraction1",
                       "Axisym2x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Check validaty of Y1fraction2 and Y2fraction2
 if (Y1fraction2<0.0 || Y1fraction2<Y1fraction1 || Y1fraction2>1.0)
  {
   throw OomphLibError("Invalid Y1fraction2",
                       "Axisym2x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Y2fraction2<0.0 || Y2fraction2<Y2fraction1 || Y2fraction2>1.0)
  {
   throw OomphLibError("Invalid Y2fraction2",
                       "Axisym2x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 
 // Number of elements in x direction
// RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb; 
 
 // // Store number of elements in bottom and top layers
//  this->Ny1 = Nya1+Nyb1+Nyc1;
//  this->Ny2 = Nya2+Nyb2+Nyc2;
 
// //  // Number of elements in y direction
// //  RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
// //  // Min. x coordinate
// //  RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
// //  // Max. x coordinate
// //  RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
// //  // Min. y coordinate
// //  RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
// //  // Max. y coordinate
// //  RectangularQuadMesh<ELEMENT >::Ymax = h1+h2;
 
// //  // Periodic?
// //  RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
//  // Store height of upper and lower layers
//  this->H1 = h1;
//  this->H2 = h2;

 // Now build the mesh: 
 this->build_two_layer_mesh(time_stepper_pt);
 
}


/// \short The spacing function for the x co-ordinates with two 
/// regions.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym2x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
x_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
 {
  // Set up some spacing parameters

  // Left region starts at Xmin
  double Xmin =  RectangularQuadMesh<ELEMENT >::Xmin;
  
  // Right region ends at Xmax
  double Xmax =  RectangularQuadMesh<ELEMENT >::Xmax;

  // Number of nodes per element (in one direction)
  unsigned n_p =  RectangularQuadMesh<ELEMENT >::Np;

  // Left region starts at Xmin
  double x1init = Xmin;

  // Right region starts at Xmin + Xfraction(Xmax-Xmin)
  double x2init = Xmin + Xfraction*(Xmax-Xmin);

  // Assuming uniform spacing, calculate the spacing between
  // the nodes in each region

  // Left region has a length Xfraction*(Xmax-Xmin)
  double x1step = Xfraction*(Xmax-Xmin)/((n_p-1)*Nxa);

  // Right region has a length (1.0-Xfraction)*(Xmax-Xmin)
  double x2step = (1.0-Xfraction)*(Xmax-Xmin)/((n_p-1)*Nxb);
  
  // Now set up the particular spacing in each region
  if(xelement < Nxa)
   {
    return (x1init + x1step*((n_p-1)*xelement + xnode));
   }
  else
   {
    return (x2init + x2step*((n_p-1)*(xelement-Nxa) + xnode));
   }
 }

/// \short The spacing function for the y co-ordinates with three
/// regions in each fluid.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym2x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
y_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
{
 //Set up some spacing parameters
 //The lower region a starts at Ymin
 double Ymin = RectangularQuadMesh<ELEMENT >::Ymin;
 //The interface is at Ymid
 double Ymid = this->H1;
 //The upper region a ends at Ymax
 double Ymax = RectangularQuadMesh<ELEMENT >::Ymax;
 //Number of nodes per element
 unsigned n_p = RectangularQuadMesh<ELEMENT >::Np;
 //The lower region a starts at Ymin
 double y1init = Ymin;
 //The lower region b starts at Ymin + Y1fraction1*(Ymid-Ymin)
 double y2init = Ymin + Y1fraction1*(Ymid-Ymin);
 //The lower region c starts at Ymin + Y1fraction2*(Ymid-Ymin)
 double y3init = Ymin + Y1fraction2*(Ymid-Ymin);
 //The upper region c starts at Ymid
 double y4init = Ymid;
 //The upper region b starts at Ymax - Y2fraction2*(Ymax-Ymid)
 double y5init = Ymax - Y2fraction2*(Ymax-Ymid);
 //The upper region a starts at Ymax - Y2fraction1*(Ymax-Ymid)
 double y6init = Ymax - Y2fraction1*(Ymax-Ymid);
 //Calculate the space between each node in each region,
 //Assumming uniform spacing
 //Lower region a has a length Y1fraction1(Ymid-Ymin)
 double y1step = Y1fraction1*(Ymid-Ymin)/((n_p-1)*Nya1);
 //Lower region b has a length (Y1fraction2-Y1fraction1)*(Ymid-Ymin)
 double y2step = (Y1fraction2-Y1fraction1)*(Ymid-Ymin)/((n_p-1)*Nyb1);
 //Lower region c has a length (1.0-Y1fraction2)*(Ymid-Ymin)
 double y3step = (1.0-Y1fraction2)*(Ymid-Ymin)/((n_p-1)*Nyc1);
 //Upper region c has a length (1.0-Y1fraction2)*(Ymax-Ymid)
 double y4step = (1.0-Y2fraction2)*(Ymax-Ymid)/((n_p-1)*Nyc2);
 //Upper region b has a length (Y1fraction2-Y1fraction1)*(Ymax-Ymid)
 double y5step = (Y2fraction2-Y2fraction1)*(Ymax-Ymid)/((n_p-1)*Nyb2);
 //Upper region a has a length Y1fraction1(Ymax-Ymid)
 double y6step = Y2fraction1*(Ymax-Ymid)/((n_p-1)*Nya2);

 //Now return the actual node position, it's different in the two
 //regions, of course
 if(yelement < Nya1) 
  {
   return (y1init + y1step*((n_p-1)*yelement + ynode));
  }
 else if (yelement < Nya1+Nyb1)
  {
   return (y2init + y2step*((n_p-1)*(yelement-Nya1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1)
  {
   return (y3init + y3step*((n_p-1)*(yelement-Nya1-Nyb1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyc2)
  {
   return (y4init + y4step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyc2+Nyb2)
  {
   return (y5init + y5step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyc2) + ynode));
  }
 else
  {
   return (y6init + 
           y6step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyc2-Nyb2) + ynode));
  }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


// //=========================================================================
// /// \short Constructor for a non-uniform two layer spine mesh, with 
// /// element layout in the lower fluid reflected in the upper. Three
// /// distinct y regions need numbers of element specified and two 
// /// x regions.
// //=========================================================================
// template <class ELEMENT, class INTERFACE_ELEMENT >
// Axisym3x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
// Axisym3x6TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb, 
//                            const unsigned &nxc,
//                            const unsigned &nya1, const unsigned &nyb1, 
//                            const unsigned &nyc1, const unsigned &nya2, 
//                            const unsigned &nyb2, const unsigned &nyc2, 
//                            const double &lx, 
//                            const double &h1, const double &h2,
//                            TimeStepper* time_stepper_pt) :
//  TwoLayerSpineMesh<ELEMENT, INTERFACE_ELEMENT >()
// {
//  // We've called the "generic" constructor for the RectangularQuadMesh
//  // which doesn't do much...
//  // Now set up the parameters that characterise the mesh geometry
//  // then call the constructor

//  Nxa = nxa; Nxb = nxb; Nxc = nxc; 
//  Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1;
//  Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2;
//  Xfraction1 = 0.2;
//  Xfraction2 = 0.8;
//  Yfraction1 = 0.2;
//  Yfraction2 = 0.8;
 
//  // Check validity of Xfraction, Yfraction1 and Yfraction2.
//  if (Yfraction1 < 0.0 || Yfraction1>Yfraction2)
//   {
//    throw OomphLibError("Invalid Yfraction1",
//                        "AxiSymm3x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//  if (Yfraction2 < Yfraction1 || Yfraction2>1.0)
//   {
//    throw OomphLibError("Invalid Yfraction2",
//                        "AxiSymm3x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//  if (Xfraction1<0.0 || Xfraction1>Xfraction2)
//   {
//    throw OomphLibError("Invalid Xfraction1",
//                        "AxiSymm3x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//  if (Xfraction2 < Xfraction1 || Xfraction2>1.0)
//   {
//     throw OomphLibError("Invalid Xfraction2",
//                         "AxiSymm3x6TwoLayerSpineMesh",
//                         OOMPH_EXCEPTION_LOCATION);
//   }

//  // Number of elements in x direction
//  RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb+Nxc; 
 
//  // Number of elements in bottom and top layers
//  this->Ny1 = Nya1+Nyb1+Nyc1;
//  this->Ny2 = Nya2+Nyb2+Nyc2;
 
//  // Number of elements in y direction
//  RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
//  // Min. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
//  // Max. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
//  // Min. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
//  // Max. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymax = h1;
 
//  // Periodic?
//  RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
//  // Set height of upper and lower layers
//  this->H1 = h1;
//  this->H2 = h2;
 
//  // Now build the mesh: 
//  this->build_two_layer_mesh(time_stepper_pt);
 
// }


//=========================================================================
/// \short Constructor for a non-uniform two layer spine mesh, with 
/// element layout in the lower fluid reflected in the upper. Three
/// distinct y regions need numbers of element specified and two 
/// x regions. The fractions of these regions are specified in this 
/// constructor.
// PATRICKFLAG I've also been editing this one...
//=========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT >
Axisym3x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
Axisym3x6TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb, 
                           const unsigned &nxc, 
                           const double &x_frac1, const double &x_frac2,
                           const unsigned &nya1, const unsigned &nyb1, 
                           const unsigned &nyc1, const unsigned &nya2, 
                           const unsigned &nyb2, const unsigned &nyc2,
                           const double &y1_frac1, const double &y1_frac2,
                           const double &y2_frac1, const double &y2_frac2,
                           const double &lx, 
                           const double &h1, const double &h2,
                           TimeStepper* time_stepper_pt) :
 TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>(nxa+nxb+nxc,nya1+nyb1+nyc1,
                                              nya2+nyb2+nyc2,lx,h1,h2,
                                              false,false,time_stepper_pt)
{
 // Mesh can only be built with 2D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

 // We've called the "generic" constructor for the TwoLayerSpineMesh and
 // set the "build_mesh" flag to false. This is done so that we do not
 // call "build_two_layer_mesh(...)" prematurely. We now set up the
 // parameters that characterise this particular mesh's geometry before
 // calling "build_two_layer_mesh(...)".

 Nxa = nxa; Nxb = nxb; Nxc = nxc; 
 Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1;
 Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2;
 Xfraction1 = x_frac1;  Xfraction2 = x_frac2;
 Y1fraction1 = y1_frac1; Y1fraction2 = y1_frac2;
 Y2fraction1 = y2_frac1; Y2fraction2 = y2_frac2;

 // Check validaty of Xfraction1
 if (Xfraction1<0.0 || Xfraction1>Xfraction2 || Xfraction1>1.0)
  {
   throw OomphLibError("Invalid Xfraction1",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Check validaty of Xfraction2
 if (Xfraction2<0.0 || Xfraction2<Xfraction1 || Xfraction2>1.0)
  {
   throw OomphLibError("Invalid Xfraction2",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Check validaty of Y1fraction1 and Y2fraction1
 if (Y1fraction1<0.0 || Y1fraction1>Y1fraction2 || Y1fraction1>1.0)
  {
   throw OomphLibError("Invalid Y1fraction1",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Y2fraction1<0.0 || Y2fraction1>Y2fraction2 || Y2fraction1>1.0)
  {
   throw OomphLibError("Invalid Y2fraction1",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Check validaty of Y1fraction2 and Y2fraction2
 if (Y1fraction2<0.0 || Y1fraction2<Y1fraction1 || Y1fraction2>1.0)
  {
   throw OomphLibError("Invalid Y1fraction2",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Y2fraction2<0.0 || Y2fraction2<Y2fraction1 || Y2fraction2>1.0)
  {
   throw OomphLibError("Invalid Y2fraction2",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
//  // Number of elements in x direction
//  RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb+Nxc; 
 
//  // Number of elements in bottom and top layers
//  this->Ny1 = Nya1+Nyb1+Nyc1;
//  this->Ny2 = Nya2+Nyb2+Nyc2;
 
//  // Number of elements in y direction
//  RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
//  // Min. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
//  // Max. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
//  // Min. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
//  // Max. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymax = h1+h2;
 
//  // Periodic?
//  RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
//  // Set height of upper and lower layers
//  this->H1 = h1;
//  this->H2 = h2;

 // Now build the mesh: 
 this->build_two_layer_mesh(time_stepper_pt);
 
}


/// \short The spacing function for the x co-ordinates with two 
/// regions.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym3x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
x_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
 {
  //Set up some spacing parameters
  //region a starts at Xmin
  double Xmin =  RectangularQuadMesh<ELEMENT >::Xmin;
  //region c ends at Xmax
  double Xmax =  RectangularQuadMesh<ELEMENT >::Xmax;
  //Number of nodes per element
  unsigned n_p =  RectangularQuadMesh<ELEMENT >::Np;
  //region a starts at Xmin
  double x1init = Xmin;
  //region b starts at Xmin + Xfraction1*(Xmax-Xmin)
  double x2init = Xmin + Xfraction1*(Xmax-Xmin);
  //region c starts at Xmin + Xfraction2*(Xmax-Xmin)
  double x3init = Xmin + Xfraction2*(Xmax-Xmin);
  //Calculate the spacing between the nodes in each region
  //Assuming uniform spacing
  //region a has a length Xfraction1*(Xmax-Xmin)
  double x1step = Xfraction1*(Xmax-Xmin)/((n_p-1)*Nxa);
  //region b has a length (Xfraction2-Xfraction1)*(Xmax-Xmin)
  double x2step = (Xfraction2-Xfraction1)*(Xmax-Xmin)/((n_p-1)*Nxb);
  //region c has a length (1.0-Xfraction2)*(Xmax-Xmin)
  double x3step = (1.0-Xfraction2)*(Xmax-Xmin)/((n_p-1)*Nxc);
  
  //Now set up the particular spacing 
  //(it's different in the two different regions)
  if(xelement < Nxa)
   {
    return (x1init + x1step*((n_p-1)*xelement + xnode));
   }
  else
   {
    if(xelement < Nxa+Nxb)
     {
      return (x2init + x2step*((n_p-1)*(xelement-Nxa) + xnode));
     }
    else
     {
      return (x3init + x3step*((n_p-1)*(xelement-Nxa-Nxb) + xnode));
     }
   }
 }

/// \short The spacing function for the y co-ordinates with three
/// regions in each fluid.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym3x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
y_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
{
 //Set up some spacing parameters
 //The lower region a starts at Ymin
 double Ymin = RectangularQuadMesh<ELEMENT >::Ymin;
 //The interface is at Ymid
 double Ymid = this->H1;
 //The upper region a ends at Ymax
 double Ymax = RectangularQuadMesh<ELEMENT >::Ymax;
 //Number of nodes per element
 unsigned n_p = RectangularQuadMesh<ELEMENT >::Np;
 //The lower region a starts at Ymin
 double y1init = Ymin;
 //The lower region b starts at Ymin + Y1fraction1*(Ymid-Ymin)
 double y2init = Ymin + Y1fraction1*(Ymid-Ymin);
 //The lower region c starts at Ymin + Y1fraction2*(Ymid-Ymin)
 double y3init = Ymin + Y1fraction2*(Ymid-Ymin);
 //The upper region c starts at Ymid
 double y4init = Ymid;
 //The upper region b starts at Ymax - Y2fraction2*(Ymax-Ymid)
 double y5init = Ymax - Y2fraction2*(Ymax-Ymid);
 //The upper region a starts at Ymax - Y2fraction1*(Ymax-Ymid)
 double y6init = Ymax - Y2fraction1*(Ymax-Ymid);
 //Calculate the space between each node in each region,
 //Assumming uniform spacing
 //Lower region a has a length Y1fraction1(Ymid-Ymin)
 double y1step = Y1fraction1*(Ymid-Ymin)/((n_p-1)*Nya1);
 //Lower region b has a length (Y1fraction2-Y1fraction1)*(Ymid-Ymin)
 double y2step = (Y1fraction2-Y1fraction1)*(Ymid-Ymin)/((n_p-1)*Nyb1);
 //Lower region c has a length (1.0-Y1fraction2)*(Ymid-Ymin)
 double y3step = (1.0-Y1fraction2)*(Ymid-Ymin)/((n_p-1)*Nyc1);
 //Upper region c has a length (1.0-Y1fraction2)*(Ymax-Ymid)
 double y4step = (1.0-Y2fraction2)*(Ymax-Ymid)/((n_p-1)*Nyc2);
 //Upper region b has a length (Y1fraction2-Y1fraction1)*(Ymax-Ymid)
 double y5step = (Y2fraction2-Y2fraction1)*(Ymax-Ymid)/((n_p-1)*Nyb2);
 //Upper region a has a length Y1fraction1(Ymax-Ymid)
 double y6step = Y2fraction1*(Ymax-Ymid)/((n_p-1)*Nya2);

 //Now return the actual node position, it's different in the two
 //regions, of course
 if(yelement < Nya1) 
  {
   return (y1init + y1step*((n_p-1)*yelement + ynode));
  }
 else if (yelement < Nya1+Nyb1)
  {
   return (y2init + y2step*((n_p-1)*(yelement-Nya1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1)
  {
   return (y3init + y3step*((n_p-1)*(yelement-Nya1-Nyb1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyc2)
  {
   return (y4init + y4step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyc2+Nyb2)
  {
   return (y5init + y5step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyc2) + ynode));
  }
 else
  {
   return (y6init + 
           y6step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyc2-Nyb2) + ynode));
  }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//=========================================================================
/// \short Constructor for a non-uniform two layer spine mesh, with 
/// element layout in the lower fluid reflected in the upper. Three
/// distinct y regions need numbers of element specified and two 
/// x regions.
//=========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT >
Axisym3x8TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
Axisym3x8TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb, 
                           const unsigned &nxc, 
                           const unsigned &nya1, const unsigned &nyb1, 
                           const unsigned &nyc1, const unsigned &nyd1, 
                           const unsigned &nya2, const unsigned &nyb2, 
                           const unsigned &nyc2, const unsigned &nyd2, 
                           const double &lx, 
                           const double &h1, const double &h2,
                           TimeStepper* time_stepper_pt) :
 TwoLayerSpineMesh<ELEMENT, INTERFACE_ELEMENT >()
{
 // Mesh can only be built with 2D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

 // We've called the "generic" constructor for the RectangularQuadMesh
 // which doesn't do much...
 // Now set up the parameters that characterise the mesh geometry
 // then call the constructor

 Nxa = nxa; Nxb = nxb; Nxc = nxc; 
 Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1; Nyd1 = nyd1;
 Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2; Nyd2 = nyd2;
 Xfraction1 = 0.2;
 Xfraction2 = 0.8;
 Yfraction1 = 0.2;
 Yfraction2 = 0.8;
 Yfraction3 = 0.93;
 
 // Check validity of Xfraction, Yfraction1 and Yfraction2.
 if (Yfraction1 < 0.0 || Yfraction1>Yfraction2)
  {
   throw OomphLibError("Invalid Yfraction1",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Yfraction2 < Yfraction1 || Yfraction2>Yfraction3)
  {
    throw OomphLibError("Invalid Yfraction2",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Yfraction3 < Yfraction2 || Yfraction3>1.0)
  {
   throw OomphLibError("Invalid Yfraction3",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Xfraction1<0.0 || Xfraction1>Xfraction2)
  {
   throw OomphLibError("Invalid Xfraction1",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Xfraction2 < Xfraction1 || Xfraction2>1.0)
  {
   throw OomphLibError("Invalid Xfraction2",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Number of elements in x direction
 RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb+Nxc; 
 
 // Number of elements in bottom and top layers
 this->Ny1 = Nya1+Nyb1+Nyc1+Nyd1;
 this->Ny2 = Nya2+Nyb2+Nyc2+Nyd2;
 
 // Number of elements in y direction
 RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
 // Min. x coordinate
 RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
 // Max. x coordinate
 RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
 // Min. y coordinate
 RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
 // Max. y coordinate
 RectangularQuadMesh<ELEMENT >::Ymax = h1+h2;
 
 // Periodic?
 RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
 // Set height of upper and lower layers
 this->H1 = h1;
 this->H2 = h2;
 
 // Now build the mesh: 
 this->build_two_layer_mesh(time_stepper_pt);
 
}


//=========================================================================
/// \short Constructor for a non-uniform two layer spine mesh, with 
/// element layout in the lower fluid reflected in the upper. Three
/// distinct y regions need numbers of element specified and two 
/// x regions. The fractions of these regions are specified in this 
/// constructor.
//=========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT >
Axisym3x8TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
Axisym3x8TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb,
                           const unsigned &nxc,
                           const double &Xfrac1, const double &Xfrac2,
                           const unsigned &nya1, const unsigned &nyb1,
                           const unsigned &nyc1, const unsigned &nyd1, 
                           const unsigned &nya2, const unsigned &nyb2,
                           const unsigned &nyc2, const unsigned &nyd2, 
                           const double &Yfrac1, const double &Yfrac2,
                           const double &Yfrac3, 
                           const double &lx, 
                           const double &h1, const double &h2,
                           TimeStepper* time_stepper_pt) :
 TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>()
{
 // We've called the "generic" constructor for the TwoLayerSpineMesh
 // which doesn't do much...
 // Now set up the parameters that characterise the mesh geometry
 // then call the constructor

 Nxa = nxa; Nxb = nxb; Nxc = nxc; 
 Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1; Nyd1 = nyd1;
 Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2; Nyd2 = nyd2;
 Xfraction1 = Xfrac1;
 Xfraction2 = Xfrac2;
 Yfraction1 = Yfrac1;
 Yfraction2 = Yfrac2;
 Yfraction3 = Yfrac3;

 // Check validity of Xfraction, Yfraction1 and Yfraction2.
 if (Yfraction1 < 0.0 || Yfraction1>Yfraction2)
  {
   throw OomphLibError("Invalid Yfraction1",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Yfraction2 < Yfraction1 || Yfraction2>Yfraction3)
  {
   throw OomphLibError("Invalid Yfraction2",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Yfraction3 < Yfraction2 || Yfraction3>1.0)
  {
   throw OomphLibError("Invalid Yfraction3",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Xfraction1<0.0 || Xfraction1>Xfraction2)
  {
    throw OomphLibError("Invalid Xfraction1",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Xfraction2 < Xfraction1 || Xfraction2>1.0)
  {
   throw OomphLibError("Invalid Xfraction2",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Number of elements in x direction
 RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb+Nxc; 
 
 // Number of elements in bottom and top layers
 this->Ny1 = Nya1+Nyb1+Nyc1+Nyd1;
 this->Ny2 = Nya1+Nyb1+Nyc1+Nyd1;
 
 // Number of elements in y direction
 RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
 // Min. x coordinate
 RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
 // Max. x coordinate
 RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
 // Min. y coordinate
 RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
 // Max. y coordinate
 RectangularQuadMesh<ELEMENT >::Ymax = h1+h2;
 
 // Periodic?
 RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
 // Set height of upper and lower layers
 this->H1 = h1;
 this->H2 = h2;

 // Now build the mesh: 
 this->build_two_layer_mesh(time_stepper_pt);
 
}


/// \short The spacing function for the x co-ordinates with two 
/// regions.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym3x8TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
x_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
 {
  //Set up some spacing parameters
  //region a starts at Xmin
  double Xmin =  RectangularQuadMesh<ELEMENT >::Xmin;
  //region c ends at Xmax
  double Xmax =  RectangularQuadMesh<ELEMENT >::Xmax;
  //Number of nodes per element
  unsigned n_p =  RectangularQuadMesh<ELEMENT >::Np;
  //region a starts at Xmin
  double x1init = Xmin;
  //region b starts at Xmin + Xfraction1*(Xmax-Xmin)
  double x2init = Xmin + Xfraction1*(Xmax-Xmin);
  //region c starts at Xmin + Xfraction2*(Xmax-Xmin)
  double x3init = Xmin + Xfraction2*(Xmax-Xmin);
  //Calculate the spacing between the nodes in each region
  //Assuming uniform spacing
  //region a has a length Xfraction1*(Xmax-Xmin)
  double x1step = Xfraction1*(Xmax-Xmin)/((n_p-1)*Nxa);
  //region b has a length (Xfraction2-Xfraction1)*(Xmax-Xmin)
  double x2step = (Xfraction2-Xfraction1)*(Xmax-Xmin)/((n_p-1)*Nxb);
  //region c has a length (1.0-Xfraction2)*(Xmax-Xmin)
  double x3step = (1.0-Xfraction2)*(Xmax-Xmin)/((n_p-1)*Nxc);
  
  //Now set up the particular spacing 
  //(it's different in the two different regions)
  if(xelement < Nxa)
   {
    return (x1init + x1step*((n_p-1)*xelement + xnode));
   }
  else
   {
    if(xelement < Nxa+Nxb)
     {
      return (x2init + x2step*((n_p-1)*(xelement-Nxa) + xnode));
     }
    else
     {
      return (x3init + x3step*((n_p-1)*(xelement-Nxa-Nxb) + xnode));
     }
   }
 }

/// \short The spacing function for the y co-ordinates with three
/// regions in each fluid.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym3x8TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
y_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
{
 //Set up some spacing parameters
 //The lower region a starts at Ymin
 double Ymin = RectangularQuadMesh<ELEMENT >::Ymin;
 //The interface is at Ymid
 double Ymid = this->H1;
 //The upper region a ends at Ymax
 double Ymax = RectangularQuadMesh<ELEMENT >::Ymax;
 //Number of nodes per element
 unsigned n_p = RectangularQuadMesh<ELEMENT >::Np;
 //The lower region a starts at Ymin
 double y1init = Ymin;
 //The lower region b starts at Ymin + Yfraction1*(Ymid-Ymin)
 double y2init = Ymin + Yfraction1*(Ymid-Ymin);
 //The lower region c starts at Ymin + Yfraction2*(Ymid-Ymin)
 double y3init = Ymin + Yfraction2*(Ymid-Ymin);
 //The lower region d starts at Ymin + Yfraction3*(Ymid-Ymin)
 double y4init = Ymin + Yfraction3*(Ymid-Ymin);
 //The upper region d starts at Ymid
 double y5init = Ymid;
 //The upper region c starts at Ymax - Yfraction3*(Ymax-Ymid)
 double y6init = Ymax - Yfraction3*(Ymax-Ymid);
 //The upper region b starts at Ymax - Yfraction2*(Ymax-Ymid)
 double y7init = Ymax - Yfraction2*(Ymax-Ymid);
 //The upper region a starts at Ymax - Yfraction1*(Ymax-Ymid)
 double y8init = Ymax - Yfraction1*(Ymax-Ymid);
 //Calculate the space between each node in each region,
 //Assumming uniform spacing
 //Lower region a has a length Yfraction1*(Ymid-Ymin)
 double y1step = Yfraction1*(Ymid-Ymin)/((n_p-1)*Nya1);
 //Lower region b has a length (Yfraction2-Yfraction1)*(Ymid-Ymin)
 double y2step = (Yfraction2-Yfraction1)*(Ymid-Ymin)/((n_p-1)*Nyb1);
 //Lower region c has a length (Yfraction3-Yfraction2)*(Ymin-Ymin)
 double y3step = (Yfraction3-Yfraction2)*(Ymid-Ymin)/((n_p-1)*Nyc1);
 //Lower region d has length (1.0-Yfraction3)*(Ymin-Ymin)
 double y4step = (1.0-Yfraction3)*(Ymid-Ymin)/((n_p-1)*Nyd1);
 //Upper region d has a length = Lower region d
 double y5step = (1.0-Yfraction3)*(Ymax-Ymid)/((n_p-1)*Nyd2);;
 //Upper region c has a length = lower region c
 double y6step = (Yfraction3-Yfraction2)*(Ymax-Ymid)/((n_p-1)*Nyc2);
 //Upper region b has a length = lower region b
 double y7step = (Yfraction2-Yfraction1)*(Ymax-Ymid)/((n_p-1)*Nyb2);
 //Upper region a has a length = lower region a
 double y8step = Yfraction1*(Ymax-Ymid)/((n_p-1)*Nya2);
 
 //Now return the actual node position, it's different in the two
 //regions, of course
 if(yelement < Nya1) 
  {
   return (y1init + y1step*((n_p-1)*yelement + ynode));
  }
 else if (yelement < Nya1+Nyb1)
  {
   return (y2init + y2step*((n_p-1)*(yelement-Nya1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1)
  {
   return (y3init + y3step*((n_p-1)*(yelement-Nya1-Nyb1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyd1)
  {
   return (y4init + y4step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyd1+Nyd2)
  {
   return (y5init + y5step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyd1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyd1+Nyd2+Nyc2)
  {
   return (y6init + 
           y6step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyd1-Nyd2) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyd1+Nyd2+Nyc2+Nyb2)
  {
   return (y7init + 
           y7step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyd1-Nyd2-Nyc2) + ynode));
  }
 else
  {
   return (y8init + 
           y8step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyd1-Nyd2-Nyc2-Nyb2) + 
                   ynode));
  }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

}
#endif
