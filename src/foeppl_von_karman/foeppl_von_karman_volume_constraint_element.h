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
//Header file for FoepplvonKarman elements
#ifndef OOMPH_FVK_VOLUME_CONSTRAINT_ELEMENT_HEADER
#define OOMPH_FVK_VOLUME_CONSTRAINT_ELEMENT_HEADER

#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

#include "../generic/elements.h"
#include "../meshes/triangle_mesh.h"
#include "Tfoeppl_von_karman_elements.h"

#include <fstream>

namespace oomph
{


//=============================================================
/// A class which allows the user to specifiy a prescribed
/// volume (as opposed to a prescribed pressure) for in the region
/// bounded by the membrane.
/// Effectively adds an equation to the system for pressure.
/// There would usually only be a single instance of this
/// element in a problem.
//=============================================================
 template <class ELEMENT, template<class> class MESH>
  class FoepplvonKarmanVolumeConstraintElement : 
 public virtual GeneralisedElement
 {
   public:

  /// \short Constructor. Takes pointer to mesh of Foeppl von Karman
  /// elements and a vector of unsigneds which identifies the
  /// regions within it that contribute to the controlled volume
  /// defined as \int w dA (i.e. the "volume under the membrane").
  /// The optional final argument specifies the initial value
  /// for the pressure that is "traded" in exchange for the
  /// volume constraint.
  FoepplvonKarmanVolumeConstraintElement(
   MESH<ELEMENT> *bounding_mesh_pt,
   const Vector<unsigned>& contributing_region,
   const double& pressure = 0.0) :
   Bounding_mesh_pt(bounding_mesh_pt),
   Contributing_region(contributing_region),
   Prescribed_volume_pt(0)
    {
     // Create instance of pressure that is traded for volume constraint
     Volume_control_pressure_pt = new Data(1);
     Volume_control_pressure_pt->set_value(0,pressure);

     // Add the data object as internal data for this element
     Pressure_data_index = add_internal_data(Volume_control_pressure_pt);
          
     // Loop over the regions that contribute to volume
     unsigned n_contributing_regions = Contributing_region.size();
     for(unsigned r = 0; r < n_contributing_regions; r++)
      {
       unsigned n_inner_el = 
        Bounding_mesh_pt->nregion_element(Contributing_region[r]);
       
       /// Loop over the elements in the contributing regions
       for(unsigned e = 0; e < n_inner_el; e++)
        {
         ELEMENT *el_pt = dynamic_cast<ELEMENT *>
          (Bounding_mesh_pt->region_element_pt(Contributing_region[r],e));
         
         // Add their nodes as external data to volume constraint element
         unsigned n_node = el_pt->nnode();
         for(unsigned n = 0; n < n_node; n++)
          {
           add_external_data(el_pt->node_pt(n));
           
           //Also add the nodal positions as external data in the 
           // case that the nodes are SolidNodes
           SolidNode *solid_node_pt = 
            dynamic_cast<SolidNode*>(el_pt->node_pt(n));
           if(solid_node_pt != 0)
            {
             add_external_data(solid_node_pt->variable_position_pt());
            }
          }
        } /// end of loop over elements
      } /// end of loop over contributing regions
    }
   
   // Destructor (empty)
   ~FoepplvonKarmanVolumeConstraintElement()
    {}
   
   /// Broken copy constructor
   FoepplvonKarmanVolumeConstraintElement(
    const FoepplvonKarmanVolumeConstraintElement&)
    { 
     BrokenCopy::broken_copy("FoepplvonKarmanVolumeConstraintElement");
    } 
   
   /// Broken assignment operator
   void operator=(const FoepplvonKarmanVolumeConstraintElement&) 
    {
     BrokenCopy::broken_assign("FoepplvonKarmanVolumeConstraintElement");
    }
   
   /// \short Returns the volume "under the elements" in the constrained
   /// regions
   double get_bounded_volume()
   {
    double bounded_volume = 0.0;
    
    // Loop over the bounded regions
    unsigned n_contributing_regions = Contributing_region.size();
    for(unsigned r = 0; r < n_contributing_regions; r++)
     {
      // Loop over the elements in the bounded regions
      unsigned n_inner_el = 
       Bounding_mesh_pt->nregion_element(Contributing_region[r]);
      for(unsigned e = 0; e < n_inner_el; e++)
       {
        // Add the contribution
        ELEMENT *el_pt = dynamic_cast<ELEMENT*>
         (Bounding_mesh_pt->region_element_pt(Contributing_region[r],e));
        bounded_volume += el_pt->get_bounded_volume();
       }
     }
    return bounded_volume;
   }
   
   /// Assign the equation number for the new equation
   void assign_additional_local_eqn_numbers()
   {
    Volume_control_local_eqn = internal_local_eqn(Pressure_data_index,0);
   }
   
   /// Add the contribution to the appropriate elemtent of the 
   /// residual vector which effectively minimises the difference
   /// between the actual volume and the prescribed volume
   void fill_in_contribution_to_residuals(Vector<double> &residuals)
   {
    if(Volume_control_local_eqn >= 0)
     {
      residuals[Volume_control_local_eqn] += get_bounded_volume()
       - (*Prescribed_volume_pt);
     }
   }
   
   /// Set pointer to target volume
   void set_prescribed_volume(double* volume_pt)
   {
    Prescribed_volume_pt = volume_pt;
   }
   
   /// \short Access to Data object whose single value contains the pressure
   /// that has been "traded" for the volume constraint.
   Data *pressure_data_pt() const
    {
     return Volume_control_pressure_pt;
    }
   
   protected:

   /// \short Data object whose single value contains the pressure
   /// that has been "traded" for the volume constraint.
   Data *Volume_control_pressure_pt;


   /// \short Unsigned indicating which internal Data object 
   /// stores the pressure
   unsigned Pressure_data_index;

   /// Local equation number of volume constraint
   int Volume_control_local_eqn;
   
   /// \short Pointer to mesh of Foeppl von Karman elements that bound
   /// the prescribed volume
   MESH<ELEMENT>* Bounding_mesh_pt;
   
   /// \short Region IDs in the bounding mesh that contribute to the 
   /// prescribed/controlled volume
   Vector<unsigned> Contributing_region;

   /// Pointer to target volume
   double* Prescribed_volume_pt;
 };
 
}

#endif
