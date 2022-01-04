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
#ifndef OOMPH_TWO_LAYER_PERTURBED_SPINE_MESH_TEMPLATE_CC
#define OOMPH_TWO_LAYER_PERTURBED_SPINE_MESH_TEMPLATE_CC

#include "two_layer_perturbed_spine_mesh.template.h"

namespace oomph
{
 /// "Magic" negative number that indicates that the indices at
 /// which the perturbations to the nodal y-positions are stored have
 /// not been set up
 template<class ELEMENT>
 int TwoLayerPerturbedSpineMesh<ELEMENT>::
 Perturbation_to_nodal_position_indices_not_set_up = -64;

 //=======================================================================
 /// Constructor for 2D perturbed spine mesh: Pass in:
 ///   - number of elements in x-direction
 ///   - number of elements in y-direction in bottom layer
 ///   - number of elements in y-direction in top layer
 ///   - axial length
 ///   - height of bottom layer
 ///   - height of top layer
 ///   - pointer to the base mesh (where we compute the base flow)
 ///   - pointer to timestepper (defaults to Steady timestepper)
 ///
 /// The mesh contains two layers of elements (of type ELEMENT;
 /// e.g  PerturbedSpineElement<QCrouzeixRaviartElement<2>)
 /// and an interfacial layer of corresponding PerturbedSpine interface
 /// elements of type INTERFACE_ELEMENT, e.g.
 /// PerturbedSpineLineFluidInterfaceElement<ELEMENT> for 2D planar
 /// problems.
 ///
 /// We set YC_index and YS_index, which store the indices at which the
 /// cosine and sine parts (respectively) of the perturbation to the
 /// nodal y-position are stored in the bulk element, to a "magic"
 /// number which indicates that the indices have not been set up yet.
 //=======================================================================
 template<class ELEMENT>
 TwoLayerPerturbedSpineMesh<ELEMENT>::
 TwoLayerPerturbedSpineMesh(const unsigned &nx,
                            const unsigned &ny1,
                            const unsigned &ny2,
                            const double &lx,
                            const double &h1,
                            const double &h2,
                            SpineMesh* &base_mesh_pt,
                            TimeStepper* time_stepper_pt) :
  RectangularQuadMesh<ELEMENT>(nx,ny1+ny2,0.0,lx,0.0,h1+h2,false,false,
                               time_stepper_pt),
  Base_mesh_pt(base_mesh_pt),
  YC_index(Perturbation_to_nodal_position_indices_not_set_up),
  YS_index(Perturbation_to_nodal_position_indices_not_set_up)
 {
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
 
 
 
 //=======================================================================
 /// Constructor for 2D perturbed spine mesh: Pass in:
 ///   - number of elements in x-direction
 ///   - number of elements in y-direction in bottom layer
 ///   - number of elements in y-direction in top layer
 ///   - axial length
 ///   - height of bottom layer
 ///   - height of top layer
 ///   - boolean flag to make the mesh periodic in the x-direction
 ///   - pointer to the base mesh (where we compute the base flow)
 ///   - pointer to timestepper (defaults to Steady timestepper)
 ///
 /// The mesh contains two layers of elements (of type ELEMENT;
 /// e.g  PerturbedSpineElement<QCrouzeixRaviartElement<2>)
 /// and an interfacial layer of corresponding PerturbedSpine interface
 /// elements of type INTERFACE_ELEMENT, e.g.
 /// PerturbedSpineLineFluidInterfaceElement<ELEMENT> for 2D planar
 /// problems.
 ///
 /// We set YC_index and YS_index, which store the indices at which the
 /// cosine and sine parts (respectively) of the perturbation to the
 /// nodal y-position are stored in the bulk element, to a "magic"
 /// number which indicates that the indices have not been set up yet.
 //=======================================================================
 template<class ELEMENT>
 TwoLayerPerturbedSpineMesh<ELEMENT>::
 TwoLayerPerturbedSpineMesh(const unsigned &nx,
                            const unsigned &ny1,
                            const unsigned &ny2,
                            const double &lx,
                            const double &h1,
                            const double &h2,
                            const bool& periodic_in_x,
                            SpineMesh* &base_mesh_pt,
                            TimeStepper* time_stepper_pt) :
  RectangularQuadMesh<ELEMENT >(nx,ny1+ny2,0.0,lx,0.0,h1+h2,
                                periodic_in_x,false,time_stepper_pt),
  Base_mesh_pt(base_mesh_pt),
  YC_index(Perturbation_to_nodal_position_indices_not_set_up),
  YS_index(Perturbation_to_nodal_position_indices_not_set_up)
 {
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

 

 //=======================================================================
 /// Constructor for 2D perturbed spine mesh: Pass in:
 ///   - number of elements in x-direction
 ///   - number of elements in y-direction in bottom layer
 ///   - number of elements in y-direction in top layer
 ///   - axial length
 ///   - height of bottom layer
 ///   - height of top layer
 ///   - boolean flag to make the mesh periodic in the x-direction
 ///   - boolean flag to specify whether or not to call the
 ///     "build_two_layer_mesh" function
 ///   - pointer to the base mesh (where we compute the base flow)
 ///   - pointer to timestepper (defaults to Steady timestepper)
 ///
 /// The mesh contains two layers of elements (of type ELEMENT;
 /// e.g  PerturbedSpineElement<QCrouzeixRaviartElement<2>)
 /// and an interfacial layer of corresponding PerturbedSpine interface
 /// elements of type INTERFACE_ELEMENT, e.g.
 /// PerturbedSpineLineFluidInterfaceElement<ELEMENT> for 2D planar
 /// problems.
 ///
 /// We set YC_index and YS_index, which store the indices at which the
 /// cosine and sine parts (respectively) of the perturbation to the
 /// nodal y-position are stored in the bulk element, to a "magic"
 /// number which indicates that the indices have not been set up yet.
 //=======================================================================
 template<class ELEMENT>
 TwoLayerPerturbedSpineMesh<ELEMENT>::
 TwoLayerPerturbedSpineMesh(const unsigned &nx,
                            const unsigned &ny1,
                            const unsigned &ny2,
                            const double &lx,
                            const double &h1,
                            const double &h2,
                            const bool& periodic_in_x,
                            const bool& build_mesh,
                            SpineMesh* &base_mesh_pt,
                            TimeStepper* time_stepper_pt) :
  RectangularQuadMesh<ELEMENT >(nx,ny1+ny2,0.0,lx,0.0,h1+h2,
                                periodic_in_x,false,time_stepper_pt),
  Base_mesh_pt(base_mesh_pt),
  YC_index(Perturbation_to_nodal_position_indices_not_set_up),
  YS_index(Perturbation_to_nodal_position_indices_not_set_up)
 {
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
  // (such as PatrickAxisym2x6TwoLayerPerturbedSpineMesh) where the mesh
  // building needs to be called from *its* constructor and this constructor
  // is only used to pass arguments to the RectangularQuadMesh constructor.
  if(build_mesh) { build_two_layer_mesh(time_stepper_pt); }
  
}


//==================================================================
/// The spacing function for the x co-ordinate, which is the
/// same as the default function.
//==================================================================
 template<class ELEMENT>
 double TwoLayerPerturbedSpineMesh<ELEMENT>::
 x_spacing_function(unsigned xelement, unsigned xnode,
                    unsigned yelement, unsigned ynode)
 {
  // Calculate the values of equal increments in nodal values
  double xstep = (this->Xmax-this->Xmin)/((this->Np-1)*this->Nx);

  // Return the appropriate value
  return (this->Xmin + xstep*((this->Np-1)*xelement + xnode));
 }

//==================================================================
/// The spacing function for the y co-ordinates, which splits
/// the region into two regions (1 and 2), according to the 
/// heights H1 and H2, with Ny1 and Ny2 elements respectively.
//==================================================================
 template<class ELEMENT>
 double TwoLayerPerturbedSpineMesh<ELEMENT>::
 y_spacing_function(unsigned xelement, unsigned xnode,
                    unsigned yelement, unsigned ynode)
 {
  // Set up some spacing parameters
  
  // The lower region a starts at Ymin
  double Ymin = RectangularQuadMesh<ELEMENT >::Ymin;

  // The upper region a ends at Ymax
  double Ymax = RectangularQuadMesh<ELEMENT >::Ymax;

  // Number of nodes per element
  unsigned n_p = RectangularQuadMesh<ELEMENT >::Np;

  // The lower region starts at Ymin
  double y1init = Ymin;

  // The upper region starts at H1 - Ymin
  double y2init = H1 - Ymin;

  // Calculate the space between each node in each region, assuming
  // uniform spacing

  // Lower region has a length (H1-Ymin)
  double y1step = (H1-Ymin)/((n_p-1)*Ny1);

  // Upper region has a length (Ymax-H1)
  double y2step = (Ymax-H1)/((n_p-1)*Ny2);
  
  // Now return the actual node position (different in the two regions)
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
 template<class ELEMENT>
 void TwoLayerPerturbedSpineMesh<ELEMENT>::
 build_two_layer_mesh(TimeStepper* time_stepper_pt) 
 {
  
  // Build the underlying quad mesh: 
  RectangularQuadMesh<ELEMENT >::build_mesh(time_stepper_pt);
  
  // Set up the pointers to elements in the upper and lower fluid
  // ------------------------------------------------------------

  // Calculate numbers of nodes in upper and lower regions
  const unsigned long n_lower = this->Nx*Ny1;
  const unsigned long n_upper = this->Nx*Ny2;

  // Loop over lower elements and push back
  Lower_layer_element_pt.reserve(n_lower);
  for(unsigned e=0;e<n_lower;e++)
   {
    Lower_layer_element_pt.push_back(this->finite_element_pt(e));
   }
  // Loop over upper elements and push back
  Upper_layer_element_pt.reserve(n_upper);
  for(unsigned e=n_lower;e<(n_lower+n_upper);e++)
   {
    Upper_layer_element_pt.push_back(this->finite_element_pt(e));
   }
 
  // Allocate memory for the spines and fractions along spines
  // ---------------------------------------------------------

  // Read out number of linear points in the element
  const unsigned n_p =
   dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

  // Allocate store for the spines:
  if(this->Xperiodic) { PerturbedSpine_pt.reserve((n_p-1)*this->Nx); }
  else { PerturbedSpine_pt.reserve((n_p-1)*this->Nx+1); }

  // FIRST SPINE
  // -----------

  // Element 0
  // Node 0

  // Create a perturbed spine and pass it a pointer to the corresponding
  // base spine
  PerturbedSpine* new_perturbed_spine_pt =
   new PerturbedSpine(time_stepper_pt,Base_mesh_pt->spine_pt(0));
  PerturbedSpine_pt.push_back(new_perturbed_spine_pt);

  // Get pointer to node
  PerturbedSpineNode* nod_pt = element_node_pt(0,0);

  // Set the pointer to the spine
  nod_pt->perturbed_spine_pt() = new_perturbed_spine_pt;
  // Set the fraction
  nod_pt->fraction() = 0.0;
  // Pointer to the mesh that implements the update fct
  nod_pt->spine_mesh_pt() = this; 
  // ID of update function within this mesh: 0 = lower; 1 = upper
  nod_pt->node_update_fct_id() = 0;

  // Pass the spine a pointer to this node
  new_perturbed_spine_pt->node_at_bottom_of_spine_pt() = nod_pt;
  
  // Loop vertically along the spine
  // Loop over the elements in fluid 1
  for(unsigned long i=0;i<Ny1;i++)
   {
    // Loop over the vertical nodes, apart from the first
    for(unsigned l1=1;l1<n_p;l1++)
     {
      // Get pointer to node
      PerturbedSpineNode* nod_pt=element_node_pt(i*this->Nx,l1*n_p);
      // Set the pointer to the spine
      nod_pt->perturbed_spine_pt() = new_perturbed_spine_pt;
      // Set the fraction
      nod_pt->fraction() = (nod_pt->x(1)-this->Ymin)/(H1);
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = this; 
      // ID of update function within this mesh: 0 = lower; 1 = upper
      nod_pt->node_update_fct_id() = 0;
     }
   }
  
  // Loop over the elements in fluid 2
  for(unsigned long i=0;i<Ny2;i++)
   {
    // Loop over vertical nodes, apart from the first
    for(unsigned l1=1;l1<n_p;l1++)
     { 
      // Get pointer to node
      PerturbedSpineNode* nod_pt=element_node_pt((Ny1+i)*this->Nx,l1*n_p);
      
      //Set the pointer to the spine
      nod_pt->perturbed_spine_pt() = new_perturbed_spine_pt;
      //Set the fraction
      nod_pt->fraction() =(nod_pt->x(1)-(this->Ymin+H1))/(H2);
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = this; 
      // ID of update function within this mesh: 0 = lower; 1 = upper
      nod_pt->node_update_fct_id() = 1;
      
      // Pass the spine a pointer to the "top" node
      if(i==(Ny2-1) && l1==(n_p-1))
       {
        new_perturbed_spine_pt->node_at_top_of_spine_pt() = nod_pt;
       }
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
      // Create a perturbed spine and pass it a pointer to the
      // corresponding base spine
      new_perturbed_spine_pt =
       new PerturbedSpine(time_stepper_pt,
                          Base_mesh_pt->spine_pt(((n_pmax-1)*j)+l2));
      PerturbedSpine_pt.push_back(new_perturbed_spine_pt);
      
      // Get pointer to node
      PerturbedSpineNode* nod_pt=element_node_pt(j,l2);
      
      //Set the pointer to the spine
      nod_pt->perturbed_spine_pt() = new_perturbed_spine_pt;
      //Set the fraction
      nod_pt->fraction() = 0.0;
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = this; 
      // ID of update function within this mesh: 0 = lower; 1 = upper
      nod_pt->node_update_fct_id() = 0;

     // Pass the spine a pointer to this node
     new_perturbed_spine_pt->node_at_bottom_of_spine_pt() = nod_pt;

     //Loop vertically along the spine
     //Loop over the elements in fluid 1
     for(unsigned long i=0;i<Ny1;i++)
      {
       //Loop over the vertical nodes, apart from the first
       for(unsigned l1=1;l1<n_p;l1++)
        {
         // Get pointer to node
         PerturbedSpineNode* nod_pt=element_node_pt(i*this->Nx+j,l1*n_p+l2);
         //Set the pointer to the spine
         nod_pt->perturbed_spine_pt() = new_perturbed_spine_pt;
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
         PerturbedSpineNode* nod_pt=element_node_pt((Ny1+i)*this->Nx+j,l1*n_p+l2);

         //Set the pointer to the spine
         nod_pt->perturbed_spine_pt() = new_perturbed_spine_pt;
         //Set the fraction
         nod_pt->fraction() = (nod_pt->x(1)-(this->Ymin+H1))/H2;
         // Pointer to the mesh that implements the update fct
         nod_pt->spine_mesh_pt() = this; 
         // ID of update function within this mesh: 0 = lower; 1 = upper
         nod_pt->node_update_fct_id() = 1;

         // Pass the spine a pointer to the "top" node
         if(i==(Ny2-1) && l1==(n_p-1))
          {
           new_perturbed_spine_pt->node_at_top_of_spine_pt() = nod_pt;
          }
        }
      }
    }
  }


 // Last spine needs special treatment for periodic meshes
 // because it's the same as the first one...
 if (this->Xperiodic)
  {
   // Last spine is the same as first one...
   PerturbedSpine* final_perturbed_spine_pt=PerturbedSpine_pt[0];   

   // Get pointer to node
   PerturbedSpineNode* nod_pt=element_node_pt((this->Nx-1),(n_p-1));

   //Set the pointer to the spine
   nod_pt->perturbed_spine_pt() = final_perturbed_spine_pt;
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
       PerturbedSpineNode* nod_pt = 
        element_node_pt(i*this->Nx+(this->Nx-1),l1*n_p+(n_p-1));

       //Set the pointer to the spine
       nod_pt->perturbed_spine_pt() = final_perturbed_spine_pt;
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
 
}


//======================================================================
/// Reorder the elements, so we loop over them vertically first
/// (advantageous in "wide" domains if a frontal solver is used).
//======================================================================
template <class ELEMENT>
void TwoLayerPerturbedSpineMesh<ELEMENT>::element_reorder()
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


}
#endif
