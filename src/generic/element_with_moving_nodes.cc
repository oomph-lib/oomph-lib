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
//Functions for the ElementWithMovingNode class
#include "element_with_moving_nodes.h"
#include "geom_objects.h"

namespace oomph
{

 /////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////
 //Functions for the ElementWithMovingNodes class
 ///////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////

 //====================================================================
 /// Return a set of all geometric data associated with the element's node
 /// update function
 //======================================================================
 void ElementWithMovingNodes::assemble_set_of_all_geometric_data(
  std::set<Data*> &unique_geom_data_pt)
 {
  //First clear the set (just in case)
  unique_geom_data_pt.clear();
  
  //Get number of nodes
  const unsigned n_node = this->nnode();
  
  // Loop over all nodes
  for(unsigned n=0;n<n_node;n++)
   {
    //Cache pointer to the Node
    Node* const nod_pt = this->node_pt(n);
    
    //Is the node hanging
    const bool node_is_hanging = nod_pt->is_hanging();
    
    // Default number of master nodes
    unsigned nmaster=1;
    
    // Default: Node isn't hanging so it's its own master node
    Node* master_node_pt=nod_pt;
    
    // Cache the hanging point
    HangInfo* hang_info_pt = 0;
    
    // Find tghe number of master nodes if the node is hanging
    if (node_is_hanging) 
     {
      hang_info_pt = nod_pt->hanging_pt();
      nmaster = hang_info_pt->nmaster();
     }
    
    // Loop over all master nodes
    for (unsigned imaster=0;imaster<nmaster;imaster++)
     {
      //Get the master node
      if (node_is_hanging)
       {master_node_pt = hang_info_pt->master_node_pt(imaster);}
      
      //Find the number of data
      const unsigned n_geom_data = master_node_pt->ngeom_data();
      //If there are geometric data add them to the set
      if(n_geom_data > 0)
       {
        //Get vector of geometric data involved in the geometric 
        //change of this node
        Data** node_geom_data_pt = master_node_pt->all_geom_data_pt();
        
        for(unsigned i=0;i<n_geom_data;i++)
         {
          unique_geom_data_pt.insert(node_geom_data_pt[i]);
         }
       }
      
          
      // Find the number of geometric objects
      unsigned n_geom_obj = master_node_pt->ngeom_object();
      //If there are geometric objects, add them to the set
      if(n_geom_obj > 0)
       {
        // Get vector of geometric objects involved in the default 
        // update function for this (master) node.
        // Vector is constructed by copy operation.
        GeomObject** geom_object_pt = master_node_pt->all_geom_object_pt();
        
        //Loop over the geometric objects
        for(unsigned i=0;i<n_geom_obj;i++)
         {
          //Get the next geometric object
          GeomObject* geom_obj_pt=geom_object_pt[i];
          
          // Number of items of geometric data that the geometric
          // object depends on
          unsigned n_geom_data=geom_obj_pt->ngeom_data();
          
          // Loop over geometric data and add to set (use set to ensure 
          // that each one is only counted once)
          for (unsigned idata=0;idata<n_geom_data;idata++)
           {
            unique_geom_data_pt.insert(geom_obj_pt->geom_data_pt(idata));
           }         
         }
       } //End of add geom object loop
     }
   }
 }
 
 //=================================================================
 /// Construct the vector of (unique) geometric data
 //=================================================================
 void ElementWithMovingNodes::complete_setup_of_dependencies()
 {
  // This set will hold the pointers to all the unique (geometric) Data that 
  // affects the shape of this element
  std::set<Data*> unique_geom_data_pt;
  //Assemble that data
  this->assemble_set_of_all_geometric_data(unique_geom_data_pt);
  
  // Resize storage for the pointers to the Data items that are
  //involved in the element's node update operation.
  Geom_data_pt.clear();
  
  // Loop over all the unique remaining Data items involved in the
  // node update operations
  typedef std::set<Data*>::iterator IT;
  for (IT it=unique_geom_data_pt.begin();it!=unique_geom_data_pt.end();it++)
   {
    Geom_data_pt.push_back(*it);
   }
 }


 //==================================================================
 /// Assign local equation numbers for the geometric data associated
 ///with the element.
 //==================================================================
 void ElementWithMovingNodes::assign_all_generic_local_eqn_numbers()
 {
  //Get local number of dofs so far
  unsigned local_eqn_number = this->ndof();
  
  //Set the number of data
  const unsigned n_geom_data = ngeom_data();
  
  //If we have any geometric data
  if(n_geom_data > 0)
   {
    // Work out total number of values involved
    // Initialise from the first object
    unsigned n_total_values = Geom_data_pt[0]->nvalue();
    //Add the values from the other data
    for(unsigned i=1;i<n_geom_data;i++)
     {n_total_values += Geom_data_pt[i]->nvalue();}
    
    //If allocated delete the old storage
    if(Geometric_data_local_eqn)
     {
      delete[] Geometric_data_local_eqn[0];
      delete[] Geometric_data_local_eqn;
     }
    
    //If there are no values, we are done, null out the storage and
    //return
    if(n_total_values==0) {Geometric_data_local_eqn=0; return;}
    
    //Resize the storage for the spine local equation numbers
    //Firstly allocate the row pointers
    Geometric_data_local_eqn = new int*[n_geom_data];
    //Now allocate storage for all the equation numbers
    Geometric_data_local_eqn[0] = new int[n_total_values];
    //Initially all local equations are unclassified
    for(unsigned i=0;i<n_total_values;i++)
     {Geometric_data_local_eqn[0][i] = Data::Is_unclassified;}
    
    //Loop over the remaining rows and set their pointers
    for(unsigned i=1;i<n_geom_data;++i)
     {
      //Initially set the pointer to the i-th row to the pointer
      //to the i-1th row
      Geometric_data_local_eqn[i] = Geometric_data_local_eqn[i-1];
      //Now increase the row pointer by the number of values 
      //stored at the i-1th geometric data
      Geometric_data_local_eqn[i] += Geom_data_pt[i-1]->nvalue();
     }
    
    
    //A local queue to store the global equation numbers
    std::deque<unsigned long> global_eqn_number_queue;
    
    //Loop over the node update data
    for(unsigned i=0;i<n_geom_data;i++)
     {
      // Pointer to spine height data
      Data* data_pt=Geom_data_pt[i];
      
      // Loop over values at this Data item
      unsigned n_value= data_pt->nvalue();
      for (unsigned j=0;j<n_value;j++)
       {
        // Get global equation number
        long eqn_number=data_pt->eqn_number(j);
        //If equation number positive 
        if (eqn_number>=0)
         {
          //Add the global equation number to our queue
          global_eqn_number_queue.push_back(eqn_number);
          //Add to local value
          Geometric_data_local_eqn[i][j] = local_eqn_number;
          local_eqn_number++;
         }
        else
         {
          //Set the local scheme to be pinned
          Geometric_data_local_eqn[i][j] = Data::Is_pinned;
         }
       }
     }
    
    //Now add our global equations numbers to the internal element storage
    this->add_global_eqn_numbers(global_eqn_number_queue);
   }
 }
 
 

 //==================================================================
 /// Calculate the node-update--related entries in the
 /// Jacobian, using finite differencing. The vector passed
 /// in residuals has to contain the nonlinear reasiduals
 /// evaluated for the current values of the unknowns.
 //==================================================================
 void ElementWithMovingNodes::fill_in_jacobian_from_geometric_data_by_fd(
  Vector<double> &residuals, DenseMatrix<double> &jacobian)
 {
  //Get number of Data items involved in node update operations
  const unsigned n_geometric_data = ngeom_data();
  //If there is nothing to be done, then leave
  if(n_geometric_data == 0) {return;}
  
  //Create a new residuals Vector
  const unsigned n_dof = this->ndof();
  
  //Create newres vector
  Vector<double> newres(n_dof);
  
  //Integer storage for the local unknown
  int local_unknown=0;
  
  //Use the default finite difference step
  const double fd_step = GeneralisedElement::Default_fd_jacobian_step;
  
  //Loop over the Data items that affect the node update operations
  for(unsigned i=0;i<n_geometric_data;i++)
   {
    //Loop over values
    unsigned n_value = Geom_data_pt[i]->nvalue();
    for(unsigned j=0;j<n_value;j++)
     {
      local_unknown = geometric_data_local_eqn(i,j);
      //If the value is free
      if(local_unknown >= 0)
       {
        //Get a pointer to the spine value
        double *value_pt = Geom_data_pt[i]->value_pt(j);
        //Save the old value
        double old_var = *value_pt;
        
        //Increment the variable
        *value_pt += fd_step;
        
        //Update the whole element (Bit inefficient)
        this->node_update();
        
        //Calculate the new residuals
        this->get_residuals(newres);
        
        //Now do finite differences
        for(unsigned m=0;m<n_dof;m++)
         {
          double sum = (newres[m] - residuals[m])/fd_step;
          //Stick the entry into the Jacobian matrix
          jacobian(m,local_unknown) = sum;
         }
        
        //Reset the variable
        *value_pt = old_var;
        
        //We're relying on the total node update in the next loop
       }
     }
   }
  
  // Node update the element one final time to get things back to
  // the original state
  this->node_update();
  
 }
 
}
