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

namespace oomph
{


//=============================================================
/// A class which allows the user to specifiy a prescribed
/// volume (as opposed to a prescribed pressure) for the bubble.
/// Effectively adds an equation to the system for pressure
/// which minimises the difference between the actual bubble
/// volume and the prescribed volume.
/// There would usually only be a single instance of this
/// element in a problem.
//=============================================================

class FoepplvonKarmanVolumeConstraintElement : 
 public virtual GeneralisedElement
{
public:
 /// Constructor
 /// Takes an optional pointer to the existing value of the pressure
 /// from a previous instance of this element class
 FoepplvonKarmanVolumeConstraintElement(
   RefineableTriangleMesh<
    ProjectableFoepplvonKarmanElement<TFoepplvonKarmanElement<3> > >
     *my_mesh_pt,
   double *inherited_pressure = 0) :
  My_mesh_pt(my_mesh_pt)
  {
   Volume_control_pressure_pt = new Data(1);
   if(inherited_pressure != 0)
    {
     Volume_control_pressure_pt->set_value(0,*inherited_pressure);
    }
   else
    {
     // Initial guess of pressure
     Volume_control_pressure_pt->set_value(0,0.1);
    }

   /// Add the data object as internal data for this element
   Pressure_data_index = add_internal_data(Volume_control_pressure_pt);

   if(My_mesh_pt->nregion() > 1)
    {
     unsigned n_inner_el;
     unsigned n_node;
     /// Loop over the bubble regions
     for(unsigned r = 1; r <= 2; r++)
      {
       n_inner_el = My_mesh_pt->nregion_element(r);
       /// Loop over the elements in the bubble regions
       for(unsigned e = 0; e < n_inner_el; e++)
        {
         ProjectableFoepplvonKarmanElement<TFoepplvonKarmanElement<3> > *el_pt
          = dynamic_cast<
           ProjectableFoepplvonKarmanElement<TFoepplvonKarmanElement<3> > *>
            (My_mesh_pt->region_element_pt(r,e));
         n_node = el_pt->nnode();

         for(unsigned n = 0; n < n_node; n++)
          {
           /// Add data for the bubble region as external data for this element
           add_external_data(el_pt->node_pt(n));
          }
        } /// end of loop over elements
      } /// end of loop over bubble regions
    }
  }

 //~FoepplvonKarmanVolumeConstraintElement()
 // {
 // }

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

 /// Returns the volume in the bubble region from the current state
 /// of the solution by looping over the elements and summing their
 /// contributions to the integral
 double get_bubble_volume()
  {
   double bubble_volume = 0;
   if(My_mesh_pt->nregion() > 1)
    {
     unsigned n_inner_el;
     /// Loop over the bubble regions
     for(unsigned r = 1; r <= 2; r++)
      {
       n_inner_el = My_mesh_pt->nregion_element(r);
       /// Loop over the elements in the bubble regions
       for(unsigned e = 0; e < n_inner_el; e++)
        {
         ProjectableFoepplvonKarmanElement<TFoepplvonKarmanElement<3> > *el_pt
          = dynamic_cast<
           ProjectableFoepplvonKarmanElement<TFoepplvonKarmanElement<3> > *>
            (My_mesh_pt->region_element_pt(r,e));
         /// Add the contribution to the total volume
         bubble_volume += el_pt->get_integral_w();
        }
      }
    }
   return bubble_volume;
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
     residuals[Volume_control_local_eqn] += get_bubble_volume()
      - Prescribed_volume;
    }
  }

 void set_prescribed_volume(const double& volume)
  {
   Prescribed_volume = volume;
  }

 Data *pressure_data_pt() const
  {
   return Volume_control_pressure_pt;
  }

protected:
 Data *Volume_control_pressure_pt;
 unsigned Pressure_data_index;
 int Volume_control_local_eqn;
 double Prescribed_volume;
 RefineableTriangleMesh<
  ProjectableFoepplvonKarmanElement<TFoepplvonKarmanElement<3> > > *My_mesh_pt;
};

}

#endif
