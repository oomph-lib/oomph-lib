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
#ifndef OOMPH_BRETHERTON_SPINE_MESH_HEADER
#define OOMPH_BRETHERTON_SPINE_MESH_HEADER


// The mesh
#include "single_layer_spine_mesh.template.h"
#include "simple_rectangular_quadmesh.template.h"

namespace oomph
{
  //============================================================================
  /// Mesh for 2D Bretherton problem -- based on single layer mesh. Templated
  /// by spine-ified Navier-Stokes element type (e.g.
  /// SpineElement<QCrouzeixRaviartElement<2> > and the corresponding
  /// interface element (e.g.
  /// SpineLineFluidInterfaceElement<SpineElement<QCrouzeixRaviartElement<2> > >
  ///
  //============================================================================
  template<class ELEMENT, class INTERFACE_ELEMENT>
  class BrethertonSpineMesh : public SingleLayerSpineMesh<ELEMENT>
  {
  public:
    /// \short Constructor. Arguments:
    /// - nx1:   Number of elements along wall in deposited film region
    /// - nx2:   Number of elements along wall in horizontal transition region
    /// - nx3:   Number of elements along wall in channel region
    /// - nhalf: Number of elements in vertical transition region (there are
    ///          twice as many elements across the channel)
    /// - nh:    Number of elements across the deposited film
    /// - h:     Thickness of deposited film
    /// - zeta0:   Start coordinate on wall
    /// - zeta1:   Wall coordinate of start of transition region
    /// - zeta2:   Wall coordinate of end of liquid filled region (inflow)
    /// - lower_wall_pt: Pointer to geometric object that represents the lower
    /// wall
    /// - upper_wall_pt: Pointer to geometric object that represents the upper
    /// wall
    /// - time_stepper_pt: Pointer to timestepper; defaults to Steady
    BrethertonSpineMesh(
      const unsigned& nx1,
      const unsigned& nx2,
      const unsigned& nx3,
      const unsigned& nh,
      const unsigned& nhalf,
      const double& h,
      GeomObject* lower_wall_pt,
      GeomObject* upper_wall_pt,
      const double& zeta_start,
      const double& zeta_transition_start,
      const double& zeta_transition_end,
      const double& zeta_end,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);


    /// Access functions for pointers to interface elements
    FiniteElement*& interface_element_pt(const unsigned long& i)
    {
      return Interface_element_pt[i];
    }

    /// Number of elements on interface
    unsigned long ninterface_element() const
    {
      return Interface_element_pt.size();
    }

    /// Access functions for pointers to elements in bulk
    FiniteElement*& bulk_element_pt(const unsigned long& i)
    {
      return Bulk_element_pt[i];
    }

    /// Number of elements in bulk
    unsigned long nbulk() const
    {
      return Bulk_element_pt.size();
    }


    /// Number of free-surface spines (i.e. excluding the dummy spines
    /// in the channel region)
    unsigned nfree_surface_spines()
    {
      unsigned np = this->finite_element_pt(0)->nnode_1d();
      return 2 * (Nx1 + Nx2 + Nhalf) * (np - 1);
    }


    /// Recalculate the spine lengths after repositioning
    double find_distance_to_free_surface(GeomObject* const& fs_geom_object_pt,
                                         Vector<double>& initial_zeta,
                                         const Vector<double>& spine_base,
                                         const Vector<double>& spine_end);

    /// Reposition the spines in response to changes in geometry
    void reposition_spines(const double& zeta_lo_transition_start,
                           const double& zeta_lo_transition_end,
                           const double& zeta_up_transition_start,
                           const double& zeta_up_transition_end);

    /// Pin all spines so the mesh can be used for computation
    /// without free surfaces
    void pin_all_spines()
    {
      unsigned n_spine = this->nspine();
      for (unsigned i = 0; i < n_spine; i++)
      {
        this->spine_pt(i)->spine_height_pt()->pin(0);
      }
    }

    /// \short General node update function implements pure virtual function
    /// defined in SpineMesh base class and performs specific update
    /// actions, depending on the node update fct id stored for each node.
    void spine_node_update(SpineNode* spine_node_pt)
    {
      unsigned id = spine_node_pt->node_update_fct_id();
      switch (id)
      {
        case 0:
          spine_node_update_film_lower(spine_node_pt);
          break;

        case 1:
          spine_node_update_horizontal_transition_lower(spine_node_pt);
          break;

        case 2:
          spine_node_update_vertical_transition_lower(spine_node_pt);
          break;

        case 3:
          spine_node_update_vertical_transition_upper(spine_node_pt);
          break;

        case 4:
          spine_node_update_horizontal_transition_upper(spine_node_pt);
          break;

        case 5:
          spine_node_update_film_upper(spine_node_pt);
          break;

        case 6:
          spine_node_update_channel(spine_node_pt);
          break;

        default:
          std::ostringstream error_message;
          error_message << "Incorrect spine update id " << id << std::endl;

          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }


    /// \short Pointer to control element (just under the symmetry line
    /// near the bubble tip, so the bubble tip is located at
    /// s=[1.0,1.0] in this element.
    ELEMENT* control_element_pt()
    {
      return Control_element_pt;
    }


    /// \short Read the value of the spine centre's vertical fraction
    double spine_centre_fraction() const
    {
      return *Spine_centre_fraction_pt;
    }

    /// \short Set the pointer to the spine centre's vertial fraction
    void set_spine_centre_fraction_pt(double* const& fraction_pt)
    {
      Spine_centre_fraction_pt = fraction_pt;
    }

  protected:
    /// \short Update function for the deposited film region in the
    /// lower part of the domain: Vertical spines
    void spine_node_update_film_lower(SpineNode* spine_node_pt)
    {
      // Get fraction along the spine
      double w = spine_node_pt->fraction();
      // Get spine height
      double h = spine_node_pt->h();

      // Get wall coordinate
      Vector<double> s_lo(1);
      s_lo[0] = spine_node_pt->spine_pt()->geom_parameter(0);

      // Get position vector to wall
      Vector<double> r_wall_lo(2);
      spine_node_pt->spine_pt()->geom_object_pt(0)->position(s_lo, r_wall_lo);

      // Update the nodal position
      spine_node_pt->x(0) = r_wall_lo[0];
      spine_node_pt->x(1) = r_wall_lo[1] + w * h;
    }


    /// \short Update function for the horizontal transitition region in the
    /// lower part of the domain: Spine points from wall to origin
    void spine_node_update_horizontal_transition_lower(SpineNode* spine_node_pt)
    {
      // Get fraction along the spine
      double w = spine_node_pt->fraction();
      // Get spine height
      double h = spine_node_pt->h();

      // Get wall coordinate
      Vector<double> s_lo(1);
      s_lo[0] = spine_node_pt->spine_pt()->geom_parameter(0);

      // Get position vector to wall
      Vector<double> r_wall_lo(2);
      spine_node_pt->spine_pt()->geom_object_pt(0)->position(s_lo, r_wall_lo);

      // Get the spine centre
      Vector<double> s_transition_lo(1), s_transition_up(1);
      s_transition_lo[0] = spine_node_pt->spine_pt()->geom_parameter(1);
      s_transition_up[0] = spine_node_pt->spine_pt()->geom_parameter(2);
      Vector<double> r_transition_lo(2), r_transition_up(2);
      spine_node_pt->spine_pt()->geom_object_pt(1)->position(s_transition_lo,
                                                             r_transition_lo);
      spine_node_pt->spine_pt()->geom_object_pt(2)->position(s_transition_up,
                                                             r_transition_up);

      Vector<double> spine_centre(2);
      // Horizontal position is always halfway between the walls
      spine_centre[0] = 0.5 * (r_transition_lo[0] + r_transition_up[0]);
      // Vertical spine centre is given by fraction between the walls
      spine_centre[1] =
        r_transition_lo[1] +
        spine_centre_fraction() * (r_transition_up[1] - r_transition_lo[1]);

      // Get vector twoards spine origin
      Vector<double> N(2);
      N[0] = spine_centre[0] - r_wall_lo[0];
      N[1] = spine_centre[1] - r_wall_lo[1];
      double inv_length = 1.0 / sqrt(N[0] * N[0] + N[1] * N[1]);
      // Update the nodal position
      spine_node_pt->x(0) = r_wall_lo[0] + w * h * N[0] * inv_length;
      spine_node_pt->x(1) = r_wall_lo[1] + w * h * N[1] * inv_length;
    }


    /// \short Update function for the vertical transitition region in the
    /// lower part of the domain: Spine points to origin
    void spine_node_update_vertical_transition_lower(SpineNode* spine_node_pt)
    {
      // Get fraction along the spine
      double w = spine_node_pt->fraction();
      // Get spine height
      double h = spine_node_pt->h();

      // Get local coordinates
      Vector<double> s_lo(1), s_up(1);
      s_lo[0] = spine_node_pt->spine_pt()->geom_parameter(0);
      s_up[0] = spine_node_pt->spine_pt()->geom_parameter(1);

      // Get position vector to wall
      Vector<double> r_lo(2), r_up(2);
      spine_node_pt->spine_pt()->geom_object_pt(0)->position(s_lo, r_lo);
      spine_node_pt->spine_pt()->geom_object_pt(1)->position(s_up, r_up);

      // Get fraction along vertical line across
      double vertical_fraction = spine_node_pt->spine_pt()->geom_parameter(2);

      // Origin of spine
      Vector<double> S0(2);
      S0[0] = r_lo[0] + vertical_fraction * (r_up[0] - r_lo[0]);
      S0[1] = r_lo[1] + vertical_fraction * (r_up[1] - r_lo[1]);


      // Get the spine centre
      Vector<double> s_transition_lo(1), s_transition_up(1);
      s_transition_lo[0] = spine_node_pt->spine_pt()->geom_parameter(3);
      s_transition_up[0] = spine_node_pt->spine_pt()->geom_parameter(4);
      Vector<double> r_transition_lo(2), r_transition_up(2);
      spine_node_pt->spine_pt()->geom_object_pt(2)->position(s_transition_lo,
                                                             r_transition_lo);
      spine_node_pt->spine_pt()->geom_object_pt(3)->position(s_transition_up,
                                                             r_transition_up);

      Vector<double> spine_centre(2);
      // Horizontal position is always halfway between the walls
      spine_centre[0] = 0.5 * (r_transition_lo[0] + r_transition_up[0]);
      // Vertical spine centre is given by fraction between the walls
      spine_centre[1] =
        r_transition_lo[1] +
        spine_centre_fraction() * (r_transition_up[1] - r_transition_lo[1]);

      // Get vector towards origin
      Vector<double> N(2);
      N[0] = spine_centre[0] - S0[0];
      N[1] = spine_centre[1] - S0[1];

      double inv_length = 1.0 / sqrt(N[0] * N[0] + N[1] * N[1]);
      // Update the nodal position
      spine_node_pt->x(0) = S0[0] + w * h * N[0] * inv_length;
      spine_node_pt->x(1) = S0[1] + w * h * N[1] * inv_length;
    }


    /// \short Update function for the vertical transitition region in the
    /// upper part of the domain: Spine points to origin
    void spine_node_update_vertical_transition_upper(SpineNode* spine_node_pt)
    {
      // Get fraction along the spine
      double w = spine_node_pt->fraction();
      // Get spine height
      double h = spine_node_pt->h();

      // Get local coordinates
      Vector<double> s_lo(1), s_up(1);
      s_lo[0] = spine_node_pt->spine_pt()->geom_parameter(0);
      s_up[0] = spine_node_pt->spine_pt()->geom_parameter(1);

      // Get position vector to wall
      Vector<double> r_lo(2), r_up(2);
      spine_node_pt->spine_pt()->geom_object_pt(0)->position(s_lo, r_lo);
      spine_node_pt->spine_pt()->geom_object_pt(1)->position(s_up, r_up);

      // Get fraction along vertical line across
      double vertical_fraction = spine_node_pt->spine_pt()->geom_parameter(2);

      // Origin of spine
      Vector<double> S0(2);
      S0[0] = r_lo[0] + vertical_fraction * (r_up[0] - r_lo[0]);
      S0[1] = r_lo[1] + vertical_fraction * (r_up[1] - r_lo[1]);


      // Get the spine centre
      Vector<double> s_transition_lo(1), s_transition_up(1);
      s_transition_lo[0] = spine_node_pt->spine_pt()->geom_parameter(3);
      s_transition_up[0] = spine_node_pt->spine_pt()->geom_parameter(4);
      Vector<double> r_transition_lo(2), r_transition_up(2);
      spine_node_pt->spine_pt()->geom_object_pt(2)->position(s_transition_lo,
                                                             r_transition_lo);
      spine_node_pt->spine_pt()->geom_object_pt(3)->position(s_transition_up,
                                                             r_transition_up);

      Vector<double> spine_centre(2);
      // Horizontal position is always halfway between the walls
      spine_centre[0] = 0.5 * (r_transition_lo[0] + r_transition_up[0]);
      // Vertical spine centre is given by fraction between the walls
      spine_centre[1] =
        r_transition_lo[1] +
        spine_centre_fraction() * (r_transition_up[1] - r_transition_lo[1]);

      // Get vector towards origin
      Vector<double> N(2);
      N[0] = spine_centre[0] - S0[0];
      N[1] = spine_centre[1] - S0[1];

      double inv_length = 1.0 / sqrt(N[0] * N[0] + N[1] * N[1]);
      // Update the nodal position
      spine_node_pt->x(0) = S0[0] + w * h * N[0] * inv_length;
      spine_node_pt->x(1) = S0[1] + w * h * N[1] * inv_length;
    }


    /// \short Update function for the horizontal transitition region in the
    /// upper part of the domain: Spine points towards origin
    void spine_node_update_horizontal_transition_upper(SpineNode* spine_node_pt)
    {
      // Get fraction along the spine
      double w = spine_node_pt->fraction();
      // Get spine height
      double h = spine_node_pt->h();

      // Get wall coordinate
      Vector<double> s_up(1);
      s_up[0] = spine_node_pt->spine_pt()->geom_parameter(0);

      // Get position vector to wall
      Vector<double> r_wall_up(2);
      spine_node_pt->spine_pt()->geom_object_pt(0)->position(s_up, r_wall_up);

      // Get the spine centre
      Vector<double> s_transition_lo(1), s_transition_up(1);
      s_transition_lo[0] = spine_node_pt->spine_pt()->geom_parameter(1);
      s_transition_up[0] = spine_node_pt->spine_pt()->geom_parameter(2);
      Vector<double> r_transition_lo(2), r_transition_up(2);
      spine_node_pt->spine_pt()->geom_object_pt(1)->position(s_transition_lo,
                                                             r_transition_lo);
      spine_node_pt->spine_pt()->geom_object_pt(2)->position(s_transition_up,
                                                             r_transition_up);

      Vector<double> spine_centre(2);
      // Horizontal position is always halfway between the walls
      spine_centre[0] = 0.5 * (r_transition_lo[0] + r_transition_up[0]);
      // Vertical spine centre is given by fraction between the walls
      spine_centre[1] =
        r_transition_lo[1] +
        spine_centre_fraction() * (r_transition_up[1] - r_transition_lo[1]);

      // Get vector towards origin
      Vector<double> N(2);
      N[0] = spine_centre[0] - r_wall_up[0];
      N[1] = spine_centre[1] - r_wall_up[1];
      double inv_length = 1.0 / sqrt(N[0] * N[0] + N[1] * N[1]);

      // Update the nodal position
      spine_node_pt->x(0) = r_wall_up[0] + w * h * N[0] * inv_length;
      spine_node_pt->x(1) = r_wall_up[1] + w * h * N[1] * inv_length;
    }


    /// \short Update function for the deposited film region in the
    /// upper part of the domain: Vertical spines
    void spine_node_update_film_upper(SpineNode* spine_node_pt)
    {
      // Get fraction along the spine
      double w = spine_node_pt->fraction();
      // Get spine height
      double h = spine_node_pt->h();

      // Get wall coordinate
      Vector<double> s_up(1);
      s_up[0] = spine_node_pt->spine_pt()->geom_parameter(0);

      // Get position vector to wall
      Vector<double> r_wall_up(2);
      spine_node_pt->spine_pt()->geom_object_pt(0)->position(s_up, r_wall_up);

      // Update the nodal position
      spine_node_pt->x(0) = r_wall_up[0];
      spine_node_pt->x(1) = r_wall_up[1] - w * h;
    }


    /// \short Update function for the nodes in the channel region ahead
    /// of the finger tip: Nodes are evenly distributed along vertical
    /// lines between the top and bottom walls
    void spine_node_update_channel(SpineNode* spine_node_pt)
    {
      // Get fraction along the spine
      double w = spine_node_pt->fraction();

      // Get upper and lower local coordinates
      Vector<double> s_lo(1), s_up(1);
      s_lo[0] = spine_node_pt->spine_pt()->geom_parameter(0);
      s_up[0] = spine_node_pt->spine_pt()->geom_parameter(1);

      // Get position vector to lower wall
      Vector<double> r_lo(2), r_up(2);
      spine_node_pt->spine_pt()->geom_object_pt(0)->position(s_lo, r_lo);
      spine_node_pt->spine_pt()->geom_object_pt(1)->position(s_up, r_up);

      // Update the nodal position
      spine_node_pt->x(0) = r_lo[0] + w * (r_up[0] - r_lo[0]);
      spine_node_pt->x(1) = r_lo[1] + w * (r_up[1] - r_lo[1]);
    }


    /// \short Initial reordering elements of the elements -- before
    /// the channel mesh is added. Vertical stacks of elements, each topped by
    /// their interface element -- this is (currently) identical to the
    /// version in the SingleLayerSpineMesh but it's important
    /// that the element reordering is maintained in exactly this form
    /// for the rest of the mesh generation process to work properly.
    /// Therefore we keep a copy of the function in here.
    void initial_element_reorder();

    /// Number of elements along wall in deposited film region
    unsigned Nx1;

    /// Number of elements along wall in horizontal transition region
    unsigned Nx2;

    /// Number of elements along wall in channel region
    unsigned Nx3;

    /// \short Number of elements in vertical transition region (there are
    /// twice as many elements across the channel)
    unsigned Nhalf;

    /// Number of elements across the deposited film
    unsigned Nh;

    /// Thickness of deposited film
    double H;

    /// Pointer to geometric object that represents the upper wall
    GeomObject* Upper_wall_pt;

    /// Pointer to geometric object that represents the lower wall
    GeomObject* Lower_wall_pt;

    /// Start coordinate on wall
    double Zeta_start;

    /// Wall coordinate of end of liquid filled region (inflow)
    double Zeta_end;

    /// Wall coordinate of start of the transition region
    double Zeta_transition_start;

    /// Wall coordinate of end of transition region
    double Zeta_transition_end;

    /// Pointer to vertical fraction of the spine centre
    double* Spine_centre_fraction_pt;

    /// Default spine fraction
    double Default_spine_centre_fraction;

    /// \short Pointer to control element (just under the symmetry line
    /// near the bubble tip; the bubble tip is located at s=[1.0,1.0]
    /// in this element
    ELEMENT* Control_element_pt;

    /// Vector of pointers to element in the fluid layer
    Vector<FiniteElement*> Bulk_element_pt;

    /// Vector of pointers to interface elements
    Vector<FiniteElement*> Interface_element_pt;
  };

} // namespace oomph

#endif
