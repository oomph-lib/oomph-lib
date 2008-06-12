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
#include "fsi.h"
#include "integral.h"
#include "algebraic_elements.h"

namespace oomph
{

//======================================================================
/// Namespace for "global" FSI functions
//======================================================================
namespace FSI_functions
{

 /// Strouhal number St = a/(UT) for application of no slip condition.
 /// Initialised to 1.0.
 double Strouhal_for_no_slip=1.0;

 /// Apply no-slip condition for N.St. on a moving wall node 
 /// u = St dR/dt, where the Strouhal number St = a/(UT) is defined by
 /// FSI_functions::Strouhal_for_no_slip and is initialised to 1.0.
 /// Note: This requires the x,y,[z] velocity components to be stored
 /// in nodal values 0,1,[2]. This is the default for all currently
 /// existing Navier-Stokes elements. If you use any others, 
 /// use this function at your own risk. 
 void apply_no_slip_on_moving_wall(Node* node_pt)
  {
   // Find out spatial dimension of node
   unsigned ndim=node_pt->ndim();
   // Assign velocity
   for (unsigned i=0;i<ndim;i++)
    {
     node_pt->set_value(i,Strouhal_for_no_slip*(node_pt->dposition_dt(i))); 
    }
  }

}

 //================================================================
 /// Constructor: Assign storage -- pass the Eulerian 
 /// dimension of the "adjacent" fluid elements and the
 /// number of local coordinates required to parametrise
 /// the wall element. E.g. for a FSIKirchhoffLoveBeam,
 /// bounding a 2D fluid domain ndim_fluid=2 and nlagr_solid=1
 /// NOTE: GeomObject constructor must be called explicitly
 /// as GeomObject(nlagr_solid, ndim_fluid) from the constructor of the
 /// derived class because of virtual inheritance!
 /// Note that, by default, we're only providing storage for fluid
 /// loading on the "front" of the element. Call 
 /// FSIWallElement::enable_fluid_loading_on_both_sides()
 /// to to enable fluid loading on the back, too.
 //====================================================================
 FSIWallElement::FSIWallElement(const unsigned& nlagr_solid, 
                                const unsigned& ndim_fluid) : 
  GeomObject(nlagr_solid, ndim_fluid) 
 {
  // Number of Gauss integration points
   unsigned n_intpt=integral_pt()->nweight();

   // By default only the "front" face is loaded
   Only_front_is_loaded_by_fluid=true;

   // Only need lookup schemes for fluid elements that are 
   // adjacent to the element's "front"
   Adjacent_fluid_element_pt.resize(1);
   Adjacent_fluid_local_coord.resize(1);
   
   // Storage for FSI fluid elements that are adjacent to Gauss point
   Adjacent_fluid_element_pt[0].resize(n_intpt);
   
   // Storage for local coordinates in FSI fluid elements that are 
   // adjacent to wall Gauss point
   Adjacent_fluid_local_coord[0].resize(n_intpt);

   // Resize of vectors of local coordinates & initialise element pointers
   for (unsigned i=0;i<n_intpt;i++)
    {
     Adjacent_fluid_element_pt[0][i]=0;
     Adjacent_fluid_local_coord[0][i].resize(ndim_fluid);
    }
   
   // Set pointer to default stress ratio
   Q_pt = &Default_Q_Value;
   
   // By default we do include the influence of the external load
   // data on the jacobian and residuals. This may have to be switched
   // off, either because the "user" deems the coupling to be
   // irrelevant or because we're solving an auxiliary solids-only
   // problem, e.g. during the assignment of initial conditions
   // for a time-dependent FSI problem.
   Add_external_load_data = true;
 }
 
 
 
 //==================================================================
 /// \short Allow element to be loaded by fluid on both
 /// sides. (Resizes containers for lookup schemes and initialises
 /// data associated with elements at the "back" of the FSIWallElement
 /// to NULL.
 //==================================================================
 void FSIWallElement::enable_fluid_loading_on_both_sides()
 {
  

  // Dimension of fluid elements has been stored in Eulerian dimension
  // of the FSIWallElement's GeomObject representation
  unsigned ndim_fluid=ndim();

  // Number of Gauss integration points
  unsigned n_intpt=integral_pt()->nweight();
  
  // Both faces are loaded
  Only_front_is_loaded_by_fluid=false;
  
  // Add storage for lookup schemes for fluid elements that are 
  // adjacent to the element's "back"
  Adjacent_fluid_element_pt.resize(2);
  Adjacent_fluid_local_coord.resize(2);

  // Storage for FSI fluid elements that are adjacent to Gauss point
  Adjacent_fluid_element_pt[1].resize(n_intpt);
  
  // Storage for local coordinates in FSI fluid elements that are 
   // adjacent to wall Gauss point
   Adjacent_fluid_local_coord[1].resize(n_intpt);

   // Resize of vectors of local coordinates & initialise element pointers
   for (unsigned i=0;i<n_intpt;i++)
    {
     Adjacent_fluid_element_pt[1][i]=0;
     Adjacent_fluid_local_coord[1][i].resize(ndim_fluid);
    }
 }


 //==================================================================
 /// Get the contribution to the load vector provided by 
 /// the adjacent fluid element: Pass number of integration point
 /// in solid element,
 /// and the unit normal vector (pointing into the fluid on the "front"
 /// of the FSIWallElement) and return the load vector.\n
 /// Note that the load is non-dimensionalised on the wall-stress scale, 
 /// i.e. it is obtained by computing the traction (on the fluid stress-scale)
 /// from the adjacent fluid element and then multiplying it by 
 /// the stress-scale-ratio \f$ Q. \f$.  \n
 /// If element is loaded on both sides, the fluid load computed
 /// from the reverse element is \b subtracted, because it's
 /// computed with the same normal vector.
 //==================================================================
 void FSIWallElement::fluid_load_vector(const unsigned& intpt,
                                        const Vector<double>& N,
                                        Vector<double>& load)
  {

   // Size of load vector
   unsigned n_load=load.size();

   // Initialise
   for (unsigned i=0;i<n_load;i++) load[i]=0.0;
   
   // Create vector for fluid load
   Vector<double> fluid_load(n_load);

   // Loop over front and back if required: Get number of fluid-loaded faces
   unsigned n_loaded_face=2;
   if (Only_front_is_loaded_by_fluid) n_loaded_face=1;
 
   for (unsigned face=0;face<n_loaded_face;face++)
    {
     //Get the local coordinate in the fluid element (copy 
     // operation for Vector)
     Vector<double> s_adjacent(adjacent_fluid_local_coord(face,intpt));
 
     //Call the load function for adjacent element
     if (Adjacent_fluid_element_pt[face][intpt]!=0)
      {
       Adjacent_fluid_element_pt[face][intpt]->
        get_load(s_adjacent,N,fluid_load);
      }
     else
      {
#ifdef PARANOID
       OomphLibWarning("Info: No adjacent element set in FSIWallElement.\n",
                       "FSIWallElement::fluid_load_vector()",
                       OOMPH_EXCEPTION_LOCATION);
#endif
      }

     // Adjust sign of fluid traction. If normal on the front
     // points into the fluid, the normal at the back points away
     // from it.
     double sign=1.0;
     if (face==1) sign=-1.0;
     
     //The load is scaled by the stress-scale ratio Q
     for(unsigned i=0;i<n_load;i++)
      {
       load[i] += fluid_load[i]*sign*q();
      }
     
    } // end of loop over faces

  }


 //==================================================================
 /// Update the nodal positions in all fluid elements that affect 
 /// the traction on this FSIWallElement
 //==================================================================
 void FSIWallElement::node_update_adjacent_fluid_elements()
 {
  // Don't update elements repeatedly
  std::map<FSIFluidElement*,bool> done;

  // Number of integration points
  unsigned n_intpt=integral_pt()->nweight();
  
  // Loop over front and back if required: Get number of fluid-loaded faces
  unsigned n_loaded_face=2;
  if (Only_front_is_loaded_by_fluid) n_loaded_face=1;
  for (unsigned face=0;face<n_loaded_face;face++)
   {
    
    // Loop over all integration points in wall element
    for (unsigned iint=0;iint<n_intpt;iint++)
     {
      // Get fluid element that affects this integration point
      FSIFluidElement* el_f_pt=Adjacent_fluid_element_pt[face][iint];
      
      // Is there an adjacent fluid element?
      if (el_f_pt!=0)
       {
        // Have we update its positions yet?
        if (!done[el_f_pt])
         {
          // Update nodal positions
          el_f_pt->node_update();
          done[el_f_pt]=true;
         }
       }
     }
   }
 }


//=================================================================
/// Static default value for the ratio of stress scales
/// used in the fluid and solid equations (default is 1.0)
//=================================================================
double FSIWallElement::Default_Q_Value=1.0;


//==========================================================================
/// This function determines the all Data that affects the load on the
/// element, and adds its global equation numbers to the 
/// local-to-global look-up scheme. Note that we only include
/// Data items into the element's Load_data if they are not already
/// included in the element's nodal positional Data, its internal
/// or external Data.
//==========================================================================
void FSIWallElement::assign_load_data_local_eqn_numbers()
{
 //Clear all the data storage
 Load_data_pt.clear(); 
 Load_data_index.clear();

 //If desired, determine the Data that affects the load but is not
 //already included in elements other generic Data.
 //The conditional test is done here, so that the vectors are cleared
 //and not re-filled if Add_external_load_data is false
 if(Add_external_load_data)
  {
   // Number of integration points
   unsigned n_intpt = integral_pt()->nweight();

   //A set of all FSIFluidElements that affect the load
   std::set<FSIFluidElement*> all_load_elements_pt;
   
   // Loop over front and back if required: Get number of fluid-loaded faces
   unsigned n_loaded_face=2;
   if (Only_front_is_loaded_by_fluid) n_loaded_face=1;
   for (unsigned face=0;face<n_loaded_face;face++)
    {  
     //Loop over the integration points and the adjacent element at
     //each integration point to the set
     for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
       //Add the element adjacent to the element into the set
       all_load_elements_pt.insert(Adjacent_fluid_element_pt[face][ipt]);
      }
    }
  
   //For safety erase any null pointers
   all_load_elements_pt.erase(0);
   
   // Storage for a pairs of load data (pointer to Data and the index
   // of the load value within this Data object) affecting the element
   std::set<std::pair<Data*,unsigned> > paired_load_data;
   
   //Storage for a set of external geometric Data affecting the element
   std::set<Data*> external_geometric_data_pt;
   
   //Now, loop over all the adjacent FSIFluidElements and get the data 
   //that affects the load.
   for(std::set<FSIFluidElement*>::iterator it=all_load_elements_pt.begin();
       it != all_load_elements_pt.end(); it++)
    {
     // Add the "direct" load data (usually velocities and pressures)
     // to the set
     (*it)->identify_load_data(paired_load_data);
     
     // Add the "indirect" load data to the set: All geometric Data that 
     // affects the nodal positions in the FSIFluidElement and thus indirectly
     // affects the load.
     (*it)->identify_geometric_data(external_geometric_data_pt);
    }
   
   // Now loop over the geometric data of the FSIWallElement itself 
   // and erase them from the external_geometric_data_pt 
   // because these data are actually intrinsic 
   // data of this element and are counted and numbered elsewhere
   unsigned n_geom_data = ngeom_data();
   for(unsigned j=0;j<n_geom_data;j++)
    {
     external_geometric_data_pt.erase(geom_data_pt(j));
    }
   
   //It is possible that the geometric data may have already been added as
   //external data. We should erase any common entries from the Load_data
   //but not touch the external data that has been set up by a "knowledgeable"
   //user
   unsigned n_external = nexternal_data();
   for(unsigned j=0;j<n_external;j++)
    {
     external_geometric_data_pt.erase(external_data_pt(j));
    }
   
   //Add the pairs of data to the load vector
   for(std::set<std::pair<Data*,unsigned> >::iterator it=
        paired_load_data.begin();
       it != paired_load_data.end();it++)
    {
     Load_data_pt.push_back(it->first);
     Load_data_index.push_back(it->second);
    }
   
   // Now we can add all the remaining geometric data to the load vector
   for(std::set<Data*>::iterator it=external_geometric_data_pt.begin();
       it != external_geometric_data_pt.end(); it++)
    {
     //Find the number of values stored in the geometric data
     unsigned n_value = (*it)->nvalue();
     //Loop over the values
     for(unsigned j=0;j<n_value;j++)
      {
       //Add data to the Load data
       Load_data_pt.push_back((*it));
       Load_data_index.push_back(j);
      }
    }
  }
 
 //External load data has now been specified
 
 //Find the number of load data
 unsigned n_load_data = nload_data();
 
 //Resize the storage for the load data local equation numbers
 //initially all local equations are unclassified
 Load_data_local_eqn.resize(n_load_data, Data::Is_unclassified);
 
 //If there are load data fill in the internal storage
 if(n_load_data > 0)
  {
   //Find the number of local equations assigned so far
   unsigned local_eqn_number = ndof();
   
   //A local queue to store the global equation numbers
   std::deque<unsigned long> global_eqn_number_queue;
 
   //Now loop over the load data again assign local equation numbers
   for(unsigned i=0;i<n_load_data;i++)
    {
     //Get the GLOBAL equation number
     long eqn_number = Load_data_pt[i]->eqn_number(Load_data_index[i]);

     //If the GLOBAL equation number is positive (a free variable)
     if(eqn_number >= 0)
      {
       //Add the GLOBAL equation number to the local queue
       global_eqn_number_queue.push_back(eqn_number);
       //Add the local equation number to the local scheme
       Load_data_local_eqn[i] = local_eqn_number;
       //Increase the local number
       local_eqn_number++;
      }
     else
      {
       //Set the local scheme to be pinned
       Load_data_local_eqn[i] = Data::Is_pinned;
      }
    }
   //Now add our global equations numbers to the internal element storage
   add_global_eqn_numbers(global_eqn_number_queue);
  }
}


//======================================================================
/// Get FE Jacobian by systematic finite differencing w.r.t.
/// nodal positition Data, internal and external Data and any
/// load Data that is not included in the previous categories.
/// This is a re-implementation of the generic FD routines with 
/// they key difference being that any updates of values are followed
/// by a node update in the adjacent fluid elements since their
/// position (and hence the shear stresses they exert onto the solid)
/// may be indirectly affected by these. For greater efficiency
/// this may be overloaded in derived classes, e.g. if it is known
/// that for a specific FSIWallElement, the internal Data does not
/// affect the nodal positions in adjacent fluid elements.
/// 
/// \todo Note: skipping the node update doesn't
/// seem to do much harm --> investigate properly and maybe allow it
/// to be switched off.
//======================================================================
void  FSIWallElement::fill_in_jacobian_from_solid_position_and_external_by_fd(
 DenseMatrix<double>& jacobian)
{
// bool use_first_order_fd=false;

 //Should probably give this a more global scope
 double fd_step = SolidFiniteElement::Default_fd_jacobian_step;
 
 //Create a new residuals Vector
 unsigned n_dof = ndof();
 Vector<double> residuals(n_dof), newres(n_dof); // , newres_minus(n_dof);
 
 //Get the residuals. Note that this computes the full residuals, so in 
 //coupled problems all solid dof entries will be computed by FDs
 get_residuals(residuals);
 
 //Get number of nodes
 unsigned n_node = nnode();
 
 //Get the number of position dofs and dimensions for the nodes
 unsigned n_position_type = this->nnodal_position_type();
 unsigned n_dim = this->nodal_dimension();
 
 //Integer for the local unknown
 int local_unknown=0;


 // FD w.r.t. nodal position Data:
 //-------------------------------

 //Loop over the nodes
 for(unsigned l=0;l<n_node;l++)
  {
   //Loop over position dofs
   for(unsigned k=0;k<n_position_type;k++)
    {
     //Loop over dimension
     for(unsigned i=0;i<n_dim;i++)
      {
       //If the variable is free
       local_unknown = position_local_eqn(l,k,i);
       if(local_unknown >= 0)
        {
         //Save old value of the (generalised) Eulerian nodal position
         double old_var = node_pt(l)->x_gen(k,i);
         
         //Increment the  (generalised) Eulerian nodal position
         node_pt(l)->x_gen(k,i) += fd_step;
         
         //I am FSI element: Need to update adjacent fluid nodes/elements
         node_update_adjacent_fluid_elements();
         
         //Calculate the new residuals
         get_residuals(newres);
         
//         if (use_first_order_fd)
          {
           //Do forward finite differences
           for(unsigned m=0;m<n_dof;m++)
            {
             //Stick the entry into the Jacobian matrix
             jacobian(m,local_unknown) = (newres[m] - residuals[m])/fd_step;
            }
          }
//          else
//           {
//            //Take backwards step for the  (generalised) Eulerian nodal
//            // position
//            node_pt(l)->x_gen(k,i) = old_var-fd_step;

//            //I am FSI element: Need to update adjacent fluid nodes/elements
//            node_update_adjacent_fluid_elements();
         
//            //Calculate the new residuals at backward position
//            get_residuals(newres_minus);

//            //Do central finite differences
//            for(unsigned m=0;m<n_dof;m++)
//             {
//              //Stick the entry into the Jacobian matrix
//              jacobian(m,local_unknown) =
//               (newres[m] - newres_minus[m])/(2.0*fd_step);
//             }
//           }

         // Reset the (generalised) Eulerian nodal position
         // (no node-update in adjacent fluid elements req'd
         // as we'll reset it next time we go around this loop)
         node_pt(l)->x_gen(k,i) = old_var;
        }
      }
    }
  }
  


 // FD w.r.t. load Data that affects the residual vector but is not
 //----------------------------------------------------------------
 // included in the internal, external or the element's own 
 //--------------------------------------------------------
 // nodal position Data
 //--------------------
 
 unsigned n_load_data = nload_data();
 for(unsigned i=0;i<n_load_data;i++)
  {
   // Free or BC?
   local_unknown = Load_data_local_eqn[i];
   if (local_unknown >= 0)
    {
     //Store a pointer to the load value
     double *value_pt = Load_data_pt[i]->value_pt(Load_data_index[i]);

     //Save the old value of the load value
     double old_var = *value_pt;
     
     //Increment the value
     *value_pt += fd_step;
     
     // I am FSI element: need to update affected fluid nodes/elements
     // because the load Data may include values that affect the nodal
     // position in adjacent fluid elements (via shear stresses!)
     node_update_adjacent_fluid_elements();
     
     //Calculate the new residuals
     get_residuals(newres);

//     if (use_first_order_fd)
      {
       //Do forward finite differences
       for(unsigned m=0;m<n_dof;m++)
        {
         //Stick the entry into the Jacobian matrix
         jacobian(m,local_unknown) = (newres[m] - residuals[m])/fd_step;
        }
      }
//      else
//       {
//        //Take backwards step
//        *value_pt = old_var-fd_step;

//        // I am FSI element: need to update affected fluid nodes/elements
//        // because the load Data may include values that affect the nodal
//        // position in adjacent fluid elements (via shear stresses!)
//        node_update_adjacent_fluid_elements();
     
//        //Calculate the new residuals at backward position
//        get_residuals(newres_minus);

//        //Do central finite differences
//        for(unsigned m=0;m<n_dof;m++)
//         {
//          //Stick the entry into the Jacobian matrix
//          jacobian(m,local_unknown) =
//           (newres[m] - newres_minus[m])/(2.0*fd_step);
//         }
//       }
          
     // Reset the variables 
     // (no node-update in adjacent fluid elements req'd
     // as we'll reset it next time we go around this loop)
     *value_pt = old_var;
    }
  }
 
 // FD w.r.t. external Data:
 //-------------------------
 
 unsigned n_external=nexternal_data();
 for (unsigned idata=0;idata<n_external;idata++)
  {
   
   Data* ext_data_pt=external_data_pt(idata);
   
   // Loop over values in Data
   unsigned initial_nvalue=ext_data_pt->nvalue();
   for (unsigned ival=0;ival<initial_nvalue;ival++)
    {
     // Free or BC?
     local_unknown = external_local_eqn(idata,ival);
     if (local_unknown >= 0)
      {
       //Store a pointer to the external data's value pointers
       double *value_pt = ext_data_pt->value_pt(ival);
       
       //Save the old value
       double old_var= *value_pt;
       
       //Increment the value
       *value_pt += fd_step;
       
       // I am FSI element: need to update affected fluid nodes/elements
       node_update_adjacent_fluid_elements();
       
       //Calculate the new residuals
       get_residuals(newres);
       

//       if (use_first_order_fd)
        {
         //Do forward finite differences
         for(unsigned m=0;m<n_dof;m++)
          {
           //Stick the entry into the Jacobian matrix
           jacobian(m,local_unknown) = (newres[m] - residuals[m])/fd_step;
          }
        }
//        else
//         {
//          //Take backwards step 
//          *value_pt = old_var-fd_step;
       
//          // I am FSI element: need to update affected fluid nodes/elements
//          node_update_adjacent_fluid_elements();

//          //Calculate the new residuals at backward position
//          get_residuals(newres_minus);

//          //Do central finite differences
//          for(unsigned m=0;m<n_dof;m++)
//           {
//            //Stick the entry into the Jacobian matrix
//            jacobian(m,local_unknown) =
//             (newres[m] - newres_minus[m])/(2.0*fd_step);
//           }
//         }

       // Reset the variables
       // (no node-update in adjacent fluid elements req'd
       // as we'll reset it next time we go around this loop)
       *value_pt = old_var;
      }
    }
  }


 // FD w.r.t. internal Data:
 //-------------------------

 unsigned n_internal=ninternal_data();
 for (unsigned idata=0;idata<n_internal;idata++)
  {
   
   Data* int_data_pt=internal_data_pt(idata);
   
   // Loop over values in Data
   unsigned initial_nvalue=int_data_pt->nvalue();
   for (unsigned ival=0;ival<initial_nvalue;ival++)
    {
     // Free or BC?
     local_unknown = internal_local_eqn(idata,ival);
     if (local_unknown >= 0)
      {
       //Store a pointer to the internal data pointer
       double *value_pt = int_data_pt->value_pt(ival);

       //Save the old value
       double old_var = *value_pt;
       
       //Increment the value
       *value_pt += fd_step;
       
       // I am FSI element: need to update affected fluid nodes/elements
       node_update_adjacent_fluid_elements();
       
       //Calculate the new residuals
       get_residuals(newres);


//       if (use_first_order_fd)
        {
         //Do forward finite differences
         for(unsigned m=0;m<n_dof;m++)
          {
           //Stick the entry into the Jacobian matrix
           jacobian(m,local_unknown) = (newres[m] - residuals[m])/fd_step;
          }
        }
//        else
//         {
//          //Take backwards step 
//          *value_pt=old_var-fd_step;
       
//          // I am FSI element: need to update affected fluid nodes/elements
//          node_update_adjacent_fluid_elements();

//          //Calculate the new residuals at backward position
//          get_residuals(newres_minus);

//          //Do central finite differences
//          for(unsigned m=0;m<n_dof;m++)
//           {
//            //Stick the entry into the Jacobian matrix
//            jacobian(m,local_unknown) =
//             (newres[m] - newres_minus[m])/(2.0*fd_step);
//           }
//         }
              
       // Reset the variables
       // (no node-update in adjacent fluid elements req'd
       // as we'll reset it next time we go around this loop)
       *value_pt = old_var;
      }
    }
  }
 
 // Cleanup and reset the lot (only needed at the end)
 node_update_adjacent_fluid_elements();
 
}



}
