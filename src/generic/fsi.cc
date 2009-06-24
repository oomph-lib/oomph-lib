///LIC// ====================================================================
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

 /// Boolean flag to specify whether to use external storage in the
 /// setup of fluid load info for solid elements.  Default value is true.
 bool Use_external_storage=true;

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


 /// Output file to document the boundary coordinate 
 /// along the FSI interface of the fluid mesh during call to
 /// setup_fluid_load_info_for_solid_elements()
 std::ofstream Doc_boundary_coordinate_file;


}

 
//===================================================================
/// Flag that allows the suppression of warning messages
//====================================================================
bool FSIWallElement::Dont_warn_about_missing_adjacent_fluid_elements=false;




 //================================================================
 /// Compete build of element: Assign storage -- pass the Eulerian 
 /// dimension of the "adjacent" fluid elements and the
 /// number of local coordinates required to parametrise
 /// the wall element. E.g. for a FSIKirchhoffLoveBeam,
 /// bounding a 2D fluid domain ndim_fluid=2 and nlagr_solid=1.
 /// Note that, by default, we're only providing storage for fluid
 /// loading on the "front" of the element. Call 
 /// FSIWallElement::enable_fluid_loading_on_both_sides()
 /// to enable fluid loading on the back, too.
 //====================================================================
 void FSIWallElement::setup_fsi_wall_element(const unsigned& nlagr_solid, 
                                             const unsigned& ndim_fluid) 
 {
  
  // Set storage for underlying GeomObject
  set_nlagrangian_and_ndim(nlagr_solid, ndim_fluid);
  
  // Set source element storage - one interaction
  this->set_ninteraction(1);
 }
 
 
 //==================================================================
 /// \short Allow element to be loaded by fluid on both
 /// sides. (Resizes containers for lookup schemes and initialises
 /// data associated with elements at the "back" of the FSIWallElement
 /// to NULL.
 //==================================================================
 void FSIWallElement::enable_fluid_loading_on_both_sides()
 {
  // Both faces are loaded
  Only_front_is_loaded_by_fluid=false;
  
  // Set source element storage - two interactions
  this->set_ninteraction(2);
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
     Vector<double> s_adjacent(external_element_local_coord(face,intpt));
 
     //Call the load function for adjacent element
     FSIFluidElement* el_f_pt=dynamic_cast<FSIFluidElement*>
      (external_element_pt(face,intpt));
     if (el_f_pt!=0)
      {
       el_f_pt->get_load(s_adjacent,N,fluid_load);
      }
     else
      {
#ifdef PARANOID
       if (!Dont_warn_about_missing_adjacent_fluid_elements)
        {
         std::ostringstream warning_stream;
         warning_stream 
          << "Info: No adjacent element set in FSIWallElement.\n\n"
          << "Note: you can disable this message by setting \n      "
          << "FSIWallElement::Dont_warn_about_missing_adjacent_fluid_elements"
          << "\n      to true or recompiling without PARANOID.\n";
         OomphLibWarning(warning_stream.str(),
                         "FSIWallElement::fluid_load_vector()",
                         OOMPH_EXCEPTION_LOCATION);
        }
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
      FSIFluidElement* el_f_pt=dynamic_cast<FSIFluidElement*>
       (external_element_pt(face,iint));
      
      // Is there an adjacent fluid element?
      if (el_f_pt!=0)
       {
        // Have we updated its positions yet?
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

 //===================================================================
/// This function identifies all data that is required for the 
/// interaction
//=======================================================================
 void FSIWallElement::identify_all_field_data_for_external_interaction(
  Vector<std::set<FiniteElement*> > const &external_elements_pt,
  std::set<std::pair<Data*,unsigned> > &paired_interaction_data) 
 {
  //Loop over each inteaction
  const unsigned n_interaction = this->ninteraction();
  for(unsigned i=0;i<n_interaction;i++)
   {
    //Loop over each element in the set 
    for(std::set<FiniteElement*>::iterator it=external_elements_pt[i].begin();
       it != external_elements_pt[i].end(); it++)
     {
      //Cast the element  the specific fluid element
      FSIFluidElement *external_fluid_el_pt = 
        dynamic_cast<FSIFluidElement*>(*it);
  
      //Only add pressure load if so desired
      if (Ignore_shear_stress_in_jacobian)
       {
        // Add the "pressure" load data
        external_fluid_el_pt->identify_pressure_data(paired_interaction_data);
       }
     else
   {
    // Add the "direct" load data (usually velocities and pressures)
    // to the set
    external_fluid_el_pt->identify_load_data(paired_interaction_data);
   }
     }
   } //End of loop over interactions
 }

 /// \short Function that must return all geometric data involved 
 /// in the desired interactions from the external element
 void FSIWallElement::identify_all_geometric_data_for_external_interaction(
  Vector<std::set<FiniteElement*> > const &external_elements_pt,
  std::set<Data*> &external_geometric_data_pt) 
 {
  //If we are ignoring the shear stress, then we can ignore all
  //external geometric data
  if(Ignore_shear_stress_in_jacobian) {return;}
  else
   {
  //Loop over each inteaction
  const unsigned n_interaction = this->ninteraction();
  for(unsigned i=0;i<n_interaction;i++)
   {
    //Loop over each element in the set 
    for(std::set<FiniteElement*>::iterator it=external_elements_pt[i].begin();
       it != external_elements_pt[i].end(); it++)
     {
     (*it)->identify_geometric_data(external_geometric_data_pt);
     }
   }
   }
 }


}
