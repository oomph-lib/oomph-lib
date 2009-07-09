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
#ifndef OOMPH_ONE_D_MESH_TEMPLATE_CC
#define OOMPH_ONE_D_MESH_TEMPLATE_CC


// oomph-lib headers
#include "one_d_mesh.template.h"



namespace oomph
{


//===========================================================================
/// The generic mesh construction routine --- this contains all the hard
/// work and is called by all constructors
//===========================================================================
template <class ELEMENT>
void OneDMesh<ELEMENT>::build_mesh(TimeStepper* time_stepper_pt) 
{
 // Set the length of the domain
 Length = Xmax - Xmin;

 // Set the number of boundaries; There are 2 boundaries.
 set_nboundary(2);

 //Allocate the storage for the pointers to the elements
 Element_pt.resize(N);

 //Allocate the memory for the first element
 Element_pt[0] = new ELEMENT;
 
 // Read out the number of nodes in the element (the member function
 // nnode_1d() is implemented in QElement)
 unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

 // We can now allocate the storarage for the pointers to the 
 // nodes in the mesh 
 Node_pt.resize(1+(n_p-1)*N);

 // Initialise the node counter
 unsigned node_count=0;


 double xinit = Xmin;

 //Set the value of the increment between nodes
 //double xstep = Length/double((n_p-1)*N);

 //Set the length of the element
 double el_length = Length/double(N);
 
 //Storage for local coordinate in element
 Vector<double> s_fraction;

 //FIRST ELEMENT 
 //-------------

 //Set the first node

 // Allocate memory for the node, using the element's own construct_node
 // function -- only the element knows what type of nodes it needs!
 Node_pt[node_count] = 
  finite_element_pt(0)->construct_boundary_node(0,time_stepper_pt);
 
 //Set the position of the node
 node_pt(node_count)->x(0) = xinit;

 //Add the node to the boundary 0
 add_boundary_node(0,Node_pt[node_count]);
  
 //Increment the counter for the nodes
 node_count++;

 //Loop over the other nodes in the first element
 for(unsigned jnod=1;jnod<n_p;jnod++)
  {
   //Allocate memory for the nodes
   Node_pt[node_count] = 
    finite_element_pt(0)->construct_node(jnod,time_stepper_pt);
   
   finite_element_pt(0)->local_fraction_of_node(jnod,s_fraction);
   
   //Set the position of the node, linear mapping
   node_pt(node_count)->x(0) = xinit + el_length*s_fraction[0];
   
   //Increment the node number
   node_count++;
  }
  



 //CENTRAL ELEMENTS
 //----------------
 for(unsigned ielem=1;ielem<(N-1);ielem++)
  {
   //Allocate memory for new element
   Element_pt[ielem] = new ELEMENT;

   //First node is same as last node in neighbouring element on the left
   finite_element_pt(ielem)->node_pt(0) = 
    finite_element_pt(ielem-1)->node_pt((n_p-1));

   //New nodes 
   for(unsigned jnod=1;jnod<n_p;jnod++)
    {
     //Allocate memory for the nodes, as before
     Node_pt[node_count] = 
      finite_element_pt(ielem)->construct_node(jnod,time_stepper_pt);
     
     //Get the local coordinate
     finite_element_pt(ielem)->local_fraction_of_node(jnod,s_fraction);

     //Set the position of the node
     node_pt(node_count)->x(0) = xinit + el_length*(ielem + s_fraction[0]);

     //Increment the node number
     node_count++;
    }

  }


 // FINAL ELEMENT
 //--------------
 
 // Allocate memory for element
 Element_pt[N-1] = new ELEMENT;

 // First node is same as in last node in neighbouring element on the left
 finite_element_pt(N-1)->node_pt(0) = finite_element_pt(N-2)->node_pt(n_p-1);
 
 // New central nodes (ignore last one which needs special treatment
 // because it's on the boundary)
 for(unsigned jnod=1;jnod<(n_p-1);jnod++)
  {
   //Allocate memory for nodes, as before
   Node_pt[node_count] = 
    finite_element_pt(N-1)->construct_node(jnod,time_stepper_pt);

   //Get the local coordinate of the node
   finite_element_pt(N-1)->local_fraction_of_node(jnod,s_fraction);
   
   //Set the position of the node
   node_pt(node_count)->x(0) = xinit + el_length*(N-1 + s_fraction[0]);

   //Increment the node number
   node_count++;
  }


 //New final node

 //Allocate memory for the node, as before
 Node_pt[node_count] = 
  finite_element_pt(N-1)->construct_boundary_node(n_p-1,time_stepper_pt);
 
 //Set the position of the node
 node_pt(node_count)->x(0) = xinit + Length;

 //Add the node to the boundary 1
 add_boundary_node(1,Node_pt[node_count]);
 
 //Increment the node number
 node_count++;

}


}
#endif
