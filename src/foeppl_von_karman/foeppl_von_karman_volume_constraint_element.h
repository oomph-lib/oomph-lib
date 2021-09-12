// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Header file for FoepplvonKarman elements
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
  /// A class which allows the user to specify a prescribed
  /// volume (as opposed to a prescribed pressure) for in the region
  /// bounded by the membrane.
  /// Effectively adds an equation to the system for pressure.
  /// There would usually only be a single instance of this
  /// element in a problem.
  //=============================================================
  template<class ELEMENT, template<class> class MESH>
  class FoepplvonKarmanVolumeConstraintElement
    : public virtual GeneralisedElement
  {
  public:
    /// \short Constructor. Takes pointer to mesh of Foeppl von Karman
    /// elements and a vector of unsigneds which identifies the
    /// regions within it that contribute to the controlled volume
    /// defined as int w dA (i.e. the "volume under the membrane").
    /// The optional final argument specifies the initial value
    /// for the pressure that is "traded" in exchange for the
    /// volume constraint.
    FoepplvonKarmanVolumeConstraintElement(
      MESH<ELEMENT>* bounding_mesh_pt,
      const Vector<unsigned>& contributing_region,
      const double& pressure = 0.0)
      : Bounding_mesh_pt(bounding_mesh_pt),
        Contributing_region(contributing_region),
        Prescribed_volume_pt(0)
    {
      // Create instance of pressure that is traded for volume constraint
      Volume_control_pressure_pt = new Data(1);
      Volume_control_pressure_pt->set_value(0, pressure);

      // Add the data object as internal data for this element
      Pressure_data_index = add_internal_data(Volume_control_pressure_pt);
    }

    // Destructor (empty)
    ~FoepplvonKarmanVolumeConstraintElement() {}

    /// Broken copy constructor
    FoepplvonKarmanVolumeConstraintElement(
      const FoepplvonKarmanVolumeConstraintElement&) = delete;

    /// Broken assignment operator
    void operator=(const FoepplvonKarmanVolumeConstraintElement&) = delete;

    /// \short Returns the volume "under the elements" in the constrained
    /// regions
    double get_bounded_volume()
    {
      double bounded_volume = 0.0;

      // Loop over the bounded regions
      unsigned n_contributing_regions = Contributing_region.size();
      for (unsigned r = 0; r < n_contributing_regions; r++)
      {
        // Loop over the elements in the bounded regions
        unsigned n_inner_el =
          Bounding_mesh_pt->nregion_element(Contributing_region[r]);
        for (unsigned e = 0; e < n_inner_el; e++)
        {
          // Add the contribution
          ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
            Bounding_mesh_pt->region_element_pt(Contributing_region[r], e));
          bounded_volume += el_pt->get_bounded_volume();
        }
      }
      return bounded_volume;
    }

    /// Assign the equation number for the new equation
    void assign_additional_local_eqn_numbers()
    {
      Volume_control_local_eqn = internal_local_eqn(Pressure_data_index, 0);
    }

    /// \short Fill in residual: Difference between actual and
    /// prescribed bounded volume.
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      if (Volume_control_local_eqn >= 0)
      {
        residuals[Volume_control_local_eqn] += -(*Prescribed_volume_pt);
      }
    }


    /// Fill in contribution to elemental residual and Jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      if (Volume_control_local_eqn >= 0)
      {
        // Only add contribution to residual; Jacobian (derivs w.r.t. to
        // this element's external data is handled by FvK elements)
        residuals[Volume_control_local_eqn] -= (*Prescribed_volume_pt);


        /*  // We are in charge... */
        /*  else */
        /*   { */
        /*    // NOTE: This is very slow but keep code alive -- can be */
        /*    // recycled if/when we ever make the Jacobians in the */
        /*    // Foeppl von Karman elements analytical. */
        /*    double t_start=TimingHelpers::timer(); */

        /*    // Setup lookup scheme between local and global data */
        /*    std::map<unsigned,unsigned> local_external_eqn; */
        /*    unsigned next=nexternal_data(); */
        /*    for (unsigned j=0;j<next;j++) */
        /*     { */
        /*      Data* data_pt=external_data_pt(j); */
        /*      unsigned nval=data_pt->nvalue(); */
        /*      for (unsigned k=0;k<nval;k++) */
        /*       { */
        /*        int local_eqn=external_local_eqn(j,k); */
        /*        if (local_eqn>=0) */
        /*         { */
        /*          int global_eqn=data_pt->eqn_number(k); */
        /*          local_external_eqn[global_eqn]=local_eqn; */
        /*         } */
        /*       } */
        /*     } */


        /*    double t_end=TimingHelpers::timer(); */
        /*    oomph_info << "Time for local_external_eqn setup: "  */
        /*               << t_end-t_start << std::endl; */
        /*    t_start=TimingHelpers::timer(); */


        /*    // Add initial contribution */
        /*    residuals[Volume_control_local_eqn] -= (*Prescribed_volume_pt); */

        /*    // Initialise total bounded volume */
        /*    double bounded_volume = 0.0; */

        /*    // Loop over the bounded regions */
        /*    unsigned n_contributing_regions = Contributing_region.size(); */
        /*    for(unsigned r = 0; r < n_contributing_regions; r++) */
        /*     { */
        /*      // Loop over the elements in the bounded regions */
        /*      unsigned n_inner_el =  */
        /*       Bounding_mesh_pt->nregion_element(Contributing_region[r]); */
        /*      for(unsigned e = 0; e < n_inner_el; e++) */
        /*       { */
        /*        // Get element */
        /*        ELEMENT *el_pt = dynamic_cast<ELEMENT*> */
        /*         (Bounding_mesh_pt->region_element_pt(Contributing_region[r],e));
         */

        /*        // Get element's contribution to bounded volume and */
        /*        // derivative w.r.t. its unknowns */
        /*        double el_bounded_volume=0.0; */

        /*        std::map<unsigned,double> d_bounded_volume_d_unknown; */
        /*        el_pt->fill_in_d_bounded_volume_d_unknown(el_bounded_volume,
         */
        /*                                                  d_bounded_volume_d_unknown);
         */

        /*        // Add contribution to bounded volume */
        /*        bounded_volume += el_bounded_volume; */


        /*        // Add contribution to Jacobian */
        /*        for (std::map<unsigned,double>::iterator it= */
        /*              d_bounded_volume_d_unknown.begin();it!= */
        /*              d_bounded_volume_d_unknown.end();it++) */
        /*         { */
        /*          unsigned global_eqn=(*it).first; */
        /*          unsigned local_eqn=local_external_eqn[global_eqn]; */
        /*          jacobian(Volume_control_local_eqn,local_eqn)+=(*it).second;
         */
        /*         } */
        /*       } */
        /*     } */
        /*    // Add contribution to residuals */
        /*    residuals[Volume_control_local_eqn]+=bounded_volume; */


        /*    t_end=TimingHelpers::timer(); */
        /*    oomph_info << "Time for second part (actual work) of
         * fill...jacobian: "  */
        /*               << t_end-t_start << std::endl; */

        /*   } */
      }
    }

    /// Set pointer to target volume
    void set_prescribed_volume(double* volume_pt)
    {
      Prescribed_volume_pt = volume_pt;
    }

    /// \short Access to Data object whose single value contains the pressure
    /// that has been "traded" for the volume constraint.
    Data* pressure_data_pt() const
    {
      return Volume_control_pressure_pt;
    }

  protected:
    /// \short Data object whose single value contains the pressure
    /// that has been "traded" for the volume constraint.
    Data* Volume_control_pressure_pt;

    /// \short Unsigned indicating which internal Data object
    /// stores the pressure
    unsigned Pressure_data_index;

    /// Local equation number of volume constraint
    int Volume_control_local_eqn;

    /// \short Pointer to mesh of Foeppl von Karman elements that bound
    /// the prescribed volume; NULL if the FvK elements that contribute
    /// to the bounding volume are in charge of adding their own
    /// contribution to the volume constraint equation (and the
    /// associated Jacobian)
    MESH<ELEMENT>* Bounding_mesh_pt;

    /// \short Region IDs in the bounding mesh that contribute to the
    /// prescribed/controlled volume
    Vector<unsigned> Contributing_region;

    /// Pointer to target volume
    double* Prescribed_volume_pt;
  };

} // namespace oomph

#endif
