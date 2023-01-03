// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_BRETHERTON_SPINE_MESH_TEMPLATE_CC
#define OOMPH_BRETHERTON_SPINE_MESH_TEMPLATE_CC

#include "bretherton_spine_mesh.template.h"
#include "../generic/mesh_as_geometric_object.h"
#include "../generic/face_element_as_geometric_object.h"

#include "single_layer_spine_mesh.template.cc"
#include "simple_rectangular_quadmesh.template.cc"


namespace oomph
{
  //======================================================================
  /// Constructor. Arguments:
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
  /// - time_stepper_pt: Pointer to timestepper; defaults to Static
  //======================================================================
  template<class ELEMENT, class INTERFACE_ELEMENT>
  BrethertonSpineMesh<ELEMENT, INTERFACE_ELEMENT>::BrethertonSpineMesh(
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
    TimeStepper* time_stepper_pt)
    : SingleLayerSpineMesh<ELEMENT>(
        2 * (nx1 + nx2 + nhalf), nh, 1.0, h, time_stepper_pt),
      Nx1(nx1),
      Nx2(nx2),
      Nx3(nx3),
      Nhalf(nhalf),
      Nh(nh),
      H(h),
      Upper_wall_pt(upper_wall_pt),
      Lower_wall_pt(lower_wall_pt),
      Zeta_start(zeta_start),
      Zeta_end(zeta_end),
      Zeta_transition_start(zeta_transition_start),
      Zeta_transition_end(zeta_transition_end),
      Spine_centre_fraction_pt(&Default_spine_centre_fraction),
      Default_spine_centre_fraction(0.5)
  {
    // Add the created elements to the bulk domain
    unsigned n_bulk = this->nelement();
    Bulk_element_pt.resize(n_bulk);
    for (unsigned e = 0; e < n_bulk; e++)
    {
      Bulk_element_pt[e] = this->finite_element_pt(e);
    }

    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Mesh can only be built with spine elements
    MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(2);

    // Initialise start of transition region to zero
    // Zeta_transition_start = 0.0;

    // Length of deposited film region
    double llayer = Zeta_transition_start - Zeta_start;

    // Length of transition region
    double d = Zeta_transition_end - Zeta_transition_start;

    // Work out radius of circular cap from lower and upper wall
    Vector<double> r_wall_lo(2), r_wall_up(2);
    Vector<double> zeta(1), s_lo(1), s_up(1);
    GeomObject *lower_sub_geom_object_pt = 0, *upper_sub_geom_object_pt = 0;

    GeomObject* lower_transition_geom_object_pt = 0;
    GeomObject* upper_transition_geom_object_pt = 0;
    Vector<double> s_transition_lo(1), s_transition_up(1);
    Vector<double> spine_centre(2);

    // Find position of lower and upper walls at start of transition region
    zeta[0] = Zeta_transition_start;
    Lower_wall_pt->position(zeta, r_wall_lo);
    Upper_wall_pt->position(zeta, r_wall_up);
    // Radius is the difference between the film thickness and the distance
    // to the lower wall
    double radius = -r_wall_lo[1] - H;

    // Check to non-symmetric mesh
    if (std::fabs(r_wall_lo[1] + r_wall_up[1]) > 1.0e-4)
    {
      oomph_info << "\n\n=================================================== "
                 << std::endl;
      oomph_info << "Warning: " << std::endl;
      oomph_info << "-------- " << std::endl;
      oomph_info << " " << std::endl;
      oomph_info << "Upper and lower walls are not symmetric at zeta=0"
                 << std::endl;
      oomph_info << "Your initial mesh will look very odd." << std::endl;
      oomph_info << "y-coordinates of walls at zeta=0 are: " << r_wall_lo[1]
                 << " and " << r_wall_up[1] << std::endl
                 << std::endl;
      oomph_info << "===================================================\n\n "
                 << std::endl;
    }


    // Add the interface elements
    // Loop over the horizontal elements
    {
      unsigned n_x = 2 * (nx1 + nx2 + nhalf);
      unsigned n_y = nh;
      for (unsigned i = 0; i < n_x; i++)
      {
        // Construct a new 1D line element on the face on which the local
        // coordinate 1 is fixed at its max. value (1) --- This is face 2
        FiniteElement* interface_element_element_pt = new INTERFACE_ELEMENT(
          this->finite_element_pt(n_x * (n_y - 1) + i), 2);

        // Push it back onto the stack
        this->Element_pt.push_back(interface_element_element_pt);

        // Push it back onto the stack of interface elements
        this->Interface_element_pt.push_back(interface_element_element_pt);
      }
    }

    // Reorder elements: Vertical stacks of elements, each topped by
    // their interface element -- this is (currently) identical to the
    // version in the SingleLayerSpineMesh but it's important
    // that element reordering is maintained in exactly this form
    // so to be on the safe side, we move the function in here.
    initial_element_reorder();

    // Store pointer to control element
    Control_element_pt = dynamic_cast<ELEMENT*>(
      this->element_pt((nx1 + nx2 + nhalf) * (nh + 1) - 2));

    // Temporary storage for boundary lookup scheme
    Vector<std::set<Node*>> set_boundary_node_pt(6);

    // Boundary 1 -> 3; 2 -> 4; 3 -> 5
    for (unsigned ibound = 1; ibound <= 3; ibound++)
    {
      unsigned numnod = this->nboundary_node(ibound);
      for (unsigned j = 0; j < numnod; j++)
      {
        set_boundary_node_pt[ibound + 2].insert(
          this->boundary_node_pt(ibound, j));
      }
    }

    // Get number of nodes per element
    unsigned nnod = this->finite_element_pt(0)->nnode();

    // Get number of nodes along element edge
    unsigned np = this->finite_element_pt(0)->nnode_1d();

    // Initialise number of elements in previous regions:
    unsigned n_prev_elements = 0;

    // We've now built the straight single-layer mesh and need to change the
    // the update functions for all nodes so that the domain
    // deforms into the Bretherton shape

    // Loop over elements in lower deposited film region
    // -------------------------------------------------
    {
      // Increments in wall coordinate
      double dzeta_el = llayer / double(nx1);
      double dzeta_node = llayer / double(nx1 * (np - 1));

      // Loop over elements horizontally
      for (unsigned i = 0; i < nx1; i++)
      {
        // Start of wall coordinate
        double zeta_lo = Zeta_start + double(i) * dzeta_el;

        // Loop over elements vertically
        for (unsigned j = 0; j < nh; j++)
        {
          // Work out element number in overall mesh
          unsigned e = n_prev_elements + i * (nh + 1) + j;

          // Get pointer to element
          FiniteElement* el_pt = this->finite_element_pt(e);

          // Loop over its nodes
          for (unsigned i0 = 0; i0 < np; i0++)
          {
            for (unsigned i1 = 0; i1 < np; i1++)
            {
              // Node number:
              unsigned n = i0 * np + i1;

              // Get spine node
              SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(n));

              // Set update fct id
              nod_pt->node_update_fct_id() = 0;

              // Provide spine with additional storage for wall coordinate
              // and wall geom object:
              if (i0 == 0)
              {
                // Get the Lagrangian coordinate in the Lower Wall
                zeta[0] = zeta_lo + double(i1) * dzeta_node;
                // Get the geometric object and local coordinate
                Lower_wall_pt->locate_zeta(
                  zeta, lower_sub_geom_object_pt, s_lo);

                // The local coordinate is a geometric parameter
                // This needs to be set (rather than added) because the
                // same spine may be visited more than once
                Vector<double> parameters(1, s_lo[0]);
                nod_pt->spine_pt()->set_geom_parameter(parameters);

                // Adjust spine height
                nod_pt->spine_pt()->height() = H;

                // The sub geom object is one (and only) geom object
                // for spine:
                Vector<GeomObject*> geom_object_pt(1);
                geom_object_pt[0] = lower_sub_geom_object_pt;

                // Pass geom object(s) to spine
                nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

                // Push the node back onto boundaries
                if (j == 0) set_boundary_node_pt[0].insert(nod_pt);
              }
            }
          }
        }
      }

      // Increment number of previous elements
      n_prev_elements += nx1 * (nh + 1);
    }

    {
      // Calculate the centre for the spine nodes in the transition region
      zeta[0] = Zeta_transition_start;
      // Get the geometric objects on the walls at the start of the transition
      // region
      Lower_wall_pt->locate_zeta(
        zeta, lower_transition_geom_object_pt, s_transition_lo);
      Upper_wall_pt->locate_zeta(
        zeta, upper_transition_geom_object_pt, s_transition_up);

      // Find the Eulerian coordinates of the walls at the transition region
      lower_transition_geom_object_pt->position(s_transition_lo, r_wall_lo);
      upper_transition_geom_object_pt->position(s_transition_up, r_wall_up);

      // Take the average of these positions to define the origin of the spines
      // in the transition region Horizontal position is always halfway
      spine_centre[0] = 0.5 * (r_wall_lo[0] + r_wall_up[0]);

      // Vertical Position is given by a specified fraction
      // between the upper and lower walls
      spine_centre[1] =
        r_wall_lo[1] + spine_centre_fraction() * (r_wall_up[1] - r_wall_lo[1]);
    }


    // Loop over elements in lower horizontal transition region
    // --------------------------------------------------------
    {
      // Increments in wall coordinate
      double dzeta_el = d / double(nx2);
      double dzeta_node = d / double(nx2 * (np - 1));

      // Loop over elements horizontally
      for (unsigned i = 0; i < nx2; i++)
      {
        // Start of wall coordinate
        double zeta_lo = Zeta_transition_start + double(i) * dzeta_el;

        // Loop over elements vertically
        for (unsigned j = 0; j < nh; j++)
        {
          // Work out element number in overall mesh
          unsigned e = n_prev_elements + i * (nh + 1) + j;

          // Get pointer to element
          FiniteElement* el_pt = this->finite_element_pt(e);

          // Loop over its nodes
          for (unsigned i0 = 0; i0 < np; i0++)
          {
            for (unsigned i1 = 0; i1 < np; i1++)
            {
              // Node number:
              unsigned n = i0 * np + i1;

              // Get spine node
              SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(n));

              // Set update id
              nod_pt->node_update_fct_id() = 1;

              // Provide spine with additional storage for wall coordinate
              if (i0 == 0)
              {
                // Get the Lagrangian coordinate in the Lower Wall
                zeta[0] = zeta_lo + double(i1) * dzeta_node;
                // Get the sub geometric object and local coordinate
                Lower_wall_pt->locate_zeta(
                  zeta, lower_sub_geom_object_pt, s_lo);

                // Pass geometric parameter to the spine
                Vector<double> parameters(3);
                parameters[0] = s_lo[0];
                parameters[1] = s_transition_lo[0];
                parameters[2] = s_transition_up[0];
                nod_pt->spine_pt()->set_geom_parameter(parameters);

                // Get position vector to wall
                lower_sub_geom_object_pt->position(s_lo, r_wall_lo);

                // Get normal vector towards origin
                Vector<double> N(2);
                N[0] = spine_centre[0] - r_wall_lo[0];
                N[1] = spine_centre[1] - r_wall_lo[1];
                double length = sqrt(N[0] * N[0] + N[1] * N[1]);
                nod_pt->spine_pt()->height() = length - radius;

                // Lower sub geom object is one (and only) geom object
                // for spine:
                Vector<GeomObject*> geom_object_pt(3);
                geom_object_pt[0] = lower_sub_geom_object_pt;
                geom_object_pt[1] = lower_transition_geom_object_pt;
                geom_object_pt[2] = upper_transition_geom_object_pt;

                // Pass geom object(s) to spine
                nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

                // Push the node back onto boundaries
                if (j == 0) set_boundary_node_pt[0].insert(nod_pt);
              }
            }
          }
        }
      }

      // Increment number of previous elements
      n_prev_elements += nx2 * (nh + 1);
    }

    // Loop over elements in lower vertical transition region
    // --------------------------------------------------------
    {
      for (unsigned i = 0; i < nhalf; i++)
      {
        // Loop over elements vertically
        for (unsigned j = 0; j < nh; j++)
        {
          // Work out element number in overall mesh
          unsigned e = n_prev_elements + i * (nh + 1) + j;

          // Get pointer to element
          FiniteElement* el_pt = this->finite_element_pt(e);

          // Loop over its nodes
          for (unsigned i0 = 0; i0 < np; i0++)
          {
            for (unsigned i1 = 0; i1 < np; i1++)
            {
              // Node number:
              unsigned n = i0 * np + i1;

              // Get spine node
              SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(n));

              // Set update id
              nod_pt->node_update_fct_id() = 2;

              // Provide spine with additional storage for fraction along
              // vertical line
              if (i0 == 0)
              {
                // Get position vectors to wall
                zeta[0] = Zeta_transition_end;
                // Get the sub geometric objects
                Lower_wall_pt->locate_zeta(
                  zeta, lower_sub_geom_object_pt, s_lo);
                Upper_wall_pt->locate_zeta(
                  zeta, upper_sub_geom_object_pt, s_up);

                lower_sub_geom_object_pt->position(s_lo, r_wall_lo);
                upper_sub_geom_object_pt->position(s_up, r_wall_up);

                // Set vertical fraction
                double vertical_fraction =
                  (double(i) + double(i1) / double(np - 1)) /
                  double(2.0 * nhalf);

                // Add the geometric parameters in order
                Vector<double> parameters(5);
                parameters[0] = s_lo[0];
                parameters[1] = s_up[0];
                parameters[2] = vertical_fraction;
                parameters[3] = s_transition_lo[0];
                parameters[4] = s_transition_up[0];
                nod_pt->spine_pt()->set_geom_parameter(parameters);

                // Origin of spine
                Vector<double> S0(2);
                S0[0] = r_wall_lo[0] +
                        vertical_fraction * (r_wall_up[0] - r_wall_lo[0]);
                S0[1] = r_wall_lo[1] +
                        vertical_fraction * (r_wall_up[1] - r_wall_lo[1]);

                // Get normal vector towards origin
                Vector<double> N(2);
                N[0] = spine_centre[0] - S0[0];
                N[1] = spine_centre[1] - S0[1];

                double length = sqrt(N[0] * N[0] + N[1] * N[1]);
                nod_pt->spine_pt()->height() = length - radius;

                // Lower and Upper wall sub geom objects affect  spine:
                Vector<GeomObject*> geom_object_pt(4);
                geom_object_pt[0] = lower_sub_geom_object_pt;
                geom_object_pt[1] = upper_sub_geom_object_pt;
                geom_object_pt[2] = lower_transition_geom_object_pt;
                geom_object_pt[3] = upper_transition_geom_object_pt;

                // Pass geom object(s) to spine
                nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

                // Push the node back onto boundaries
                if (j == 0) set_boundary_node_pt[1].insert(nod_pt);
              }
            }
          }
        }
      }

      // Increment number of previous elements
      n_prev_elements += nhalf * (nh + 1);
    }


    // Loop over elements in upper vertical transition region
    // ------------------------------------------------------
    {
      for (unsigned i = 0; i < nhalf; i++)
      {
        // Loop over elements vertically
        for (unsigned j = 0; j < nh; j++)
        {
          // Work out element number in overall mesh
          unsigned e = n_prev_elements + i * (nh + 1) + j;

          // Get pointer to element
          FiniteElement* el_pt = this->finite_element_pt(e);

          // Loop over its nodes
          for (unsigned i0 = 0; i0 < np; i0++)
          {
            for (unsigned i1 = 0; i1 < np; i1++)
            {
              // Node number:
              unsigned n = i0 * np + i1;

              // Get spine node
              SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(n));

              // Set update id
              nod_pt->node_update_fct_id() = 3;

              // Provide spine with additional storage for fraction along
              // vertical line
              if (i0 == 0)
              {
                // Get position vectors to wall
                zeta[0] = Zeta_transition_end;
                // Get the sub geometric objects
                Lower_wall_pt->locate_zeta(
                  zeta, lower_sub_geom_object_pt, s_lo);
                Upper_wall_pt->locate_zeta(
                  zeta, upper_sub_geom_object_pt, s_up);

                lower_sub_geom_object_pt->position(s_lo, r_wall_lo);
                upper_sub_geom_object_pt->position(s_up, r_wall_up);

                // Set vertical fraction
                double vertical_fraction =
                  0.5 + (double(i) + double(i1) / double(np - 1)) /
                          double(2.0 * nhalf);

                // Add the geometric parameters in order
                Vector<double> parameters(5);
                parameters[0] = s_lo[0];
                parameters[1] = s_up[0];
                parameters[2] = vertical_fraction;
                parameters[3] = s_transition_lo[0];
                parameters[4] = s_transition_up[0];
                nod_pt->spine_pt()->set_geom_parameter(parameters);

                // Origin of spine
                Vector<double> S0(2);
                S0[0] = r_wall_lo[0] +
                        vertical_fraction * (r_wall_up[0] - r_wall_lo[0]);
                S0[1] = r_wall_lo[1] +
                        vertical_fraction * (r_wall_up[1] - r_wall_lo[1]);

                // Get normal vector towards origin
                Vector<double> N(2);
                N[0] = spine_centre[0] - S0[0];
                N[1] = spine_centre[1] - S0[1];

                double length = sqrt(N[0] * N[0] + N[1] * N[1]);
                nod_pt->spine_pt()->height() = length - radius;

                // Upper and Lower wall geom objects affect spine
                Vector<GeomObject*> geom_object_pt(4);
                geom_object_pt[0] = lower_sub_geom_object_pt;
                geom_object_pt[1] = upper_sub_geom_object_pt;
                geom_object_pt[2] = lower_transition_geom_object_pt;
                geom_object_pt[3] = upper_transition_geom_object_pt;

                // Pass geom object(s) to spine
                nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

                // Push the node back onto boundaries
                if (j == 0) set_boundary_node_pt[1].insert(nod_pt);
              }
            }
          }
        }
      }
      // Increment number of previous elements
      n_prev_elements += nhalf * (nh + 1);
    }


    // Loop over elements in upper horizontal transition region
    // --------------------------------------------------------
    {
      // Increments in wall coordinate
      double dzeta_el = d / double(nx2);
      double dzeta_node = d / double(nx2 * (np - 1));

      // Loop over elements horizontally
      for (unsigned i = 0; i < nx2; i++)
      {
        // Start of wall coordinate
        double zeta_lo = Zeta_transition_end - double(i) * dzeta_el;

        // Loop over elements vertically
        for (unsigned j = 0; j < nh; j++)
        {
          // Work out element number in overall mesh
          unsigned e = n_prev_elements + i * (nh + 1) + j;

          // Get pointer to element
          FiniteElement* el_pt = this->finite_element_pt(e);

          // Loop over its nodes
          for (unsigned i0 = 0; i0 < np; i0++)
          {
            for (unsigned i1 = 0; i1 < np; i1++)
            {
              // Node number:
              unsigned n = i0 * np + i1;

              // Get spine node
              SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(n));

              // Set update id
              nod_pt->node_update_fct_id() = 4;

              // Provide spine with additional storage for wall coordinate
              if (i0 == 0)
              {
                // Get the Lagrangian coordinate in the Upper wall
                zeta[0] = zeta_lo - double(i1) * dzeta_node;
                // Get the sub geometric object and local coordinate
                Upper_wall_pt->locate_zeta(
                  zeta, upper_sub_geom_object_pt, s_up);

                // Pass geometric parameter to spine
                Vector<double> parameters(3);
                parameters[0] = s_up[0];
                parameters[1] = s_transition_lo[0];
                parameters[2] = s_transition_up[0];
                nod_pt->spine_pt()->set_geom_parameter(parameters);

                // Get position vector to wall
                upper_sub_geom_object_pt->position(s_up, r_wall_up);

                // Get normal vector towards origin
                Vector<double> N(2);
                N[0] = spine_centre[0] - r_wall_up[0];
                N[1] = spine_centre[1] - r_wall_up[1];
                double length = sqrt(N[0] * N[0] + N[1] * N[1]);
                nod_pt->spine_pt()->height() = length - radius;

                // upper wall sub geom object is one (and only) geom object
                // for spine:
                Vector<GeomObject*> geom_object_pt(3);
                geom_object_pt[0] = upper_sub_geom_object_pt;
                geom_object_pt[1] = lower_transition_geom_object_pt;
                geom_object_pt[2] = upper_transition_geom_object_pt;

                // Pass geom object(s) to spine
                nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

                // Push the node back onto boundaries
                if (j == 0) set_boundary_node_pt[2].insert(nod_pt);
              }
            }
          }
        }
      }

      // Increment number of previous elements
      n_prev_elements += nx2 * (nh + 1);
    }


    // Loop over elements in upper deposited film region
    // -------------------------------------------------
    {
      // Increments in wall coordinate
      double dzeta_el = llayer / double(nx1);
      double dzeta_node = llayer / double(nx1 * (np - 1));

      // Loop over elements horizontally
      for (unsigned i = 0; i < nx1; i++)
      {
        // Start of wall coordinate
        double zeta_lo = Zeta_transition_start - double(i) * dzeta_el;

        // Loop over elements vertically
        for (unsigned j = 0; j < nh; j++)
        {
          // Work out element number in overall mesh
          unsigned e = n_prev_elements + i * (nh + 1) + j;

          // Get pointer to element
          FiniteElement* el_pt = this->finite_element_pt(e);

          // Loop over its nodes
          for (unsigned i0 = 0; i0 < np; i0++)
          {
            for (unsigned i1 = 0; i1 < np; i1++)
            {
              // Node number:
              unsigned n = i0 * np + i1;

              // Get spine node
              SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(n));

              // Set update id
              nod_pt->node_update_fct_id() = 5;

              // Provide spine with additional storage for wall coordinate
              if (i0 == 0)
              {
                // Get the Lagrangian coordinate in the Upper wall
                zeta[0] = zeta_lo - double(i1) * dzeta_node;
                // Get the geometric object and local coordinate
                Upper_wall_pt->locate_zeta(
                  zeta, upper_sub_geom_object_pt, s_up);

                // The local coordinate is a geometric parameter
                Vector<double> parameters(1, s_up[0]);
                nod_pt->spine_pt()->set_geom_parameter(parameters);

                // Adjust spine height
                nod_pt->spine_pt()->height() = H;

                // upper sub geom object is one (and only) geom object
                // for spine:
                Vector<GeomObject*> geom_object_pt(1);
                geom_object_pt[0] = upper_sub_geom_object_pt;

                // Pass geom object(s) to spine
                nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

                // Push the node back onto boundaries
                if (j == 0) set_boundary_node_pt[2].insert(nod_pt);
              }
            }
          }
        }
      }
    }

    // Wipe the boundary lookup schemes
    this->remove_boundary_nodes();
    this->set_nboundary(6);

    // Copy from sets to vectors
    for (unsigned ibound = 0; ibound < 6; ibound++)
    {
      typedef std::set<Node*>::iterator IT;
      for (IT it = set_boundary_node_pt[ibound].begin();
           it != set_boundary_node_pt[ibound].end();
           it++)
      {
        this->add_boundary_node(ibound, *it);
      }
    }


    // Loop over the free surface boundary (4) and set a boundary coordinate
    {
      // Boundary coordinate
      Vector<double> zeta(1);
      unsigned n_boundary_node = this->nboundary_node(4);
      for (unsigned n = 0; n < n_boundary_node; ++n)
      {
        zeta[0] = this->boundary_node_pt(4, n)->x(0);
        this->boundary_node_pt(4, n)->set_coordinates_on_boundary(4, zeta);
      }
    }


    // Now add the rectangular mesh in the channel ahead of the finger tip
    //--------------------------------------------------------------------

    // Build a temporary version of a SimpleRectangularQuadMesh as
    // a unit square
    SimpleRectangularQuadMesh<ELEMENT>* aux_mesh_pt =
      new SimpleRectangularQuadMesh<ELEMENT>(
        Nx3, 2 * nhalf, 1.0, 1.0, time_stepper_pt);

    // Wipe the boundary information in the auxilliary mesh
    aux_mesh_pt->remove_boundary_nodes();

    // Copy all nodes in the auxiliary mesh into a set (from where they
    // can easily be removed)
    nnod = aux_mesh_pt->nnode();
    std::set<Node*> aux_node_pt;
    for (unsigned i = 0; i < nnod; i++)
    {
      aux_node_pt.insert(aux_mesh_pt->node_pt(i));
    }

    // Loop over elements in first column and set pointers
    // to the nodes in their first column to the ones that already exist

    // Boundary node number for first node in element
    unsigned first_bound_node = 0;

    // Loop over elements
    for (unsigned e = 0; e < 2 * nhalf; e++)
    {
      FiniteElement* el_pt = aux_mesh_pt->finite_element_pt(e * Nx3);
      // Loop over first column of nodes
      for (unsigned i = 0; i < np; i++)
      {
        // Node number in element
        unsigned n = i * np;
        // Remove node from set and kill it
        if ((e < 1) || (i > 0))
        {
          aux_node_pt.erase(el_pt->node_pt(n));
          delete el_pt->node_pt(n);
        }
        // Set pointer to existing node
        el_pt->node_pt(n) = this->boundary_node_pt(1, first_bound_node + i);
      }

      // Increment first node number
      first_bound_node += np - 1;
    }

    // Now add all the remaining nodes in the auxiliary mesh
    // to the actual one
    typedef std::set<Node*>::iterator IT;
    for (IT it = aux_node_pt.begin(); it != aux_node_pt.end(); it++)
    {
      this->Node_pt.push_back(*it);
    }

    // Store number of elements before the new bit gets added:
    unsigned nelement_orig = this->Element_pt.size();

    // Add the elements to the actual mesh
    unsigned nelem = aux_mesh_pt->nelement();
    for (unsigned e = 0; e < nelem; e++)
    {
      this->Element_pt.push_back(aux_mesh_pt->element_pt(e));
      // Don't forget to add them to the bulk
      this->Bulk_element_pt.push_back(aux_mesh_pt->finite_element_pt(e));
    }

    // Now wipe the boundary node storage scheme for the
    // nodes that used to be at the end of the domain:
    this->remove_boundary_nodes(1);


    // FIRST SPINE
    //-----------

    // Element 0
    // Node 0
    // Assign the new spine with unit height (pinned dummy -- never used)
    Spine* new_spine_pt = new Spine(1.0);
    this->Spine_pt.push_back(new_spine_pt);
    new_spine_pt->spine_height_pt()->pin(0);

    // Get pointer to node
    SpineNode* nod_pt = this->element_node_pt(nelement_orig + 0, 0);
    // Set the pointer to node
    nod_pt->spine_pt() = new_spine_pt;
    // Set the fraction
    nod_pt->fraction() = 0.0;
    // Pointer to the mesh that implements the update fct
    nod_pt->spine_mesh_pt() = this;
    // ID for the update function
    nod_pt->node_update_fct_id() = 6;

    // Need to get the local coordinates for the upper and lower wall
    zeta[0] = Zeta_transition_end;
    // Get the sub geometric objects
    Lower_wall_pt->locate_zeta(zeta, lower_sub_geom_object_pt, s_lo);
    Upper_wall_pt->locate_zeta(zeta, upper_sub_geom_object_pt, s_up);

    lower_sub_geom_object_pt->position(s_lo, r_wall_lo);
    upper_sub_geom_object_pt->position(s_up, r_wall_up);

    // Pass additional data to spine
    Vector<double> parameters(2);
    parameters[0] = s_lo[0];
    parameters[1] = s_up[0];
    new_spine_pt->set_geom_parameter(parameters);

    // Lower and upper wall sub geom objects affect update
    // for spine:
    Vector<GeomObject*> geom_object_pt(2);
    geom_object_pt[0] = lower_sub_geom_object_pt;
    geom_object_pt[1] = upper_sub_geom_object_pt;

    // Pass geom object(s) to spine
    new_spine_pt->set_geom_object_pt(geom_object_pt);

    // Loop vertically along the spine
    // Loop over the elements
    for (unsigned long i = 0; i < 2 * nhalf; i++)
    {
      // Loop over the vertical nodes, apart from the first
      for (unsigned l1 = 1; l1 < np; l1++)
      {
        // Get pointer to node
        SpineNode* nod_pt =
          this->element_node_pt(nelement_orig + i * Nx3, l1 * np);
        // Set the pointer to the spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() =
          (double(i) + double(l1) / double(np - 1)) / double(2 * nhalf);
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;
        // ID for the update function
        nod_pt->node_update_fct_id() = 6;
      }
    }


    // LOOP OVER OTHER SPINES
    //----------------------

    // Now loop over the elements horizontally
    for (unsigned long j = 0; j < Nx3; j++)
    {
      // Loop over the nodes in the elements horizontally, ignoring
      // the first column
      for (unsigned l2 = 1; l2 < np; l2++)
      {
        // Assign the new spine with unit height (pinned dummy; never used)
        new_spine_pt = new Spine(1.0);
        this->Spine_pt.push_back(new_spine_pt);
        new_spine_pt->spine_height_pt()->pin(0);

        // Get the node
        SpineNode* nod_pt = this->element_node_pt(nelement_orig + j, l2);
        // Set the pointer to the spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() = 0.0;
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;
        // ID for the update function
        nod_pt->node_update_fct_id() = 6;

        // Add to boundary lookup scheme
        this->add_boundary_node(0, nod_pt);
        if ((j == Nx3 - 1) && (l2 == np - 1))
        {
          this->add_boundary_node(1, nod_pt);
        }

        // Increment in wall coordinate
        double dzeta_el = (Zeta_end - Zeta_transition_end) / double(Nx3);
        double dzeta_nod = dzeta_el / double(np - 1);

        // Get wall coordinate
        zeta[0] = Zeta_transition_end + j * dzeta_el + l2 * dzeta_nod;

        // Get the sub geometric objects
        Lower_wall_pt->locate_zeta(zeta, lower_sub_geom_object_pt, s_lo);
        Upper_wall_pt->locate_zeta(zeta, upper_sub_geom_object_pt, s_up);

        lower_sub_geom_object_pt->position(s_lo, r_wall_lo);
        upper_sub_geom_object_pt->position(s_up, r_wall_up);

        // Add geometric parameters to spine
        Vector<double> parameters(2);
        parameters[0] = s_lo[0];
        parameters[1] = s_up[0];
        new_spine_pt->set_geom_parameter(parameters);

        // Lower and upper sub geom objects affect update
        // for spine:
        Vector<GeomObject*> geom_object_pt(2);
        geom_object_pt[0] = lower_sub_geom_object_pt;
        geom_object_pt[1] = upper_sub_geom_object_pt;

        // Pass geom object(s) to spine
        new_spine_pt->set_geom_object_pt(geom_object_pt);


        // Loop vertically along the spine
        // Loop over the elements
        for (unsigned long i = 0; i < 2 * nhalf; i++)
        {
          // Loop over the vertical nodes, apart from the first
          for (unsigned l1 = 1; l1 < np; l1++)
          {
            // Get node
            SpineNode* nod_pt =
              this->element_node_pt(nelement_orig + i * Nx3 + j, l1 * np + l2);
            // Set the pointer to the spine
            nod_pt->spine_pt() = new_spine_pt;
            // Set the fraction
            nod_pt->fraction() =
              (double(i) + double(l1) / double(np - 1)) / double(2 * nhalf);
            // Pointer to the mesh that implements the update fct
            nod_pt->spine_mesh_pt() = this;
            // ID for the update function
            nod_pt->node_update_fct_id() = 6;

            // Add to boundary lookup scheme
            if ((j == Nx3 - 1) && (l2 == np - 1))
            {
              this->add_boundary_node(1, nod_pt);
            }

            // Add to boundary lookup scheme
            if ((i == 2 * nhalf - 1) && (l1 == np - 1))
            {
              this->add_boundary_node(2, nod_pt);
            }
          }
        }
      }
    }

    // (Re-)setup lookup scheme that establishes which elements are located
    // on the mesh boundaries
    this->setup_boundary_element_info();

    // Flush the storage for elements and nodes in the auxiliary mesh
    // so it can be safely deleted
    aux_mesh_pt->flush_element_and_node_storage();

    // Kill the auxiliary mesh
    delete aux_mesh_pt;
  }


  //======================================================================
  /// Reorder elements: Vertical stacks of elements, each topped by
  /// their interface element -- this is (currently) identical to the
  /// version in the SingleLayerSpineMesh but it's important
  /// that element reordering is maintained in exactly this form
  /// so to be on the safe side, we move the function in here.
  //======================================================================
  template<class ELEMENT, class INTERFACE_ELEMENT>
  void BrethertonSpineMesh<ELEMENT,
                           INTERFACE_ELEMENT>::initial_element_reorder()
  {
    unsigned Nx = this->Nx;
    unsigned Ny = this->Ny;
    // Find out how many elements there are
    unsigned long Nelement = this->nelement();
    // Find out how many fluid elements there are
    unsigned long Nfluid = Nx * Ny;
    // Create a dummy array of elements
    Vector<FiniteElement*> dummy;

    // Loop over the elements in horizontal order
    for (unsigned long j = 0; j < Nx; j++)
    {
      // Loop over the elements in lower layer vertically
      for (unsigned long i = 0; i < Ny; i++)
      {
        // Push back onto the new stack
        dummy.push_back(this->finite_element_pt(Nx * i + j));
      }

      // Push back the line element onto the stack
      dummy.push_back(this->finite_element_pt(Nfluid + j));
    }

    // Now copy the reordered elements into the element_pt
    for (unsigned long e = 0; e < Nelement; e++)
    {
      this->Element_pt[e] = dummy[e];
    }
  }

  //=======================================================================
  /// Calculate the distance from the spine base to the free surface, i.e.
  /// the spine height.
  //=======================================================================
  template<class ELEMENT, class INTERFACE_ELEMENT>
  double BrethertonSpineMesh<ELEMENT, INTERFACE_ELEMENT>::
    find_distance_to_free_surface(GeomObject* const& fs_geom_object_pt,
                                  Vector<double>& initial_zeta,
                                  const Vector<double>& spine_base,
                                  const Vector<double>& spine_end)
  {
    // Now we need to find the intersection points
    double lambda = 0.5;

    Vector<double> dx(2);
    Vector<double> new_free_x(2);
    DenseDoubleMatrix jacobian(2, 2, 0.0);
    Vector<double> spine_x(2);
    Vector<double> free_x(2);

    double maxres = 100.0;

    unsigned count = 0;
    // Let's Newton method it
    do
    {
      count++;
      // Find the spine's location
      for (unsigned i = 0; i < 2; ++i)
      {
        spine_x[i] = spine_base[i] + lambda * (spine_end[i] - spine_base[i]);
      }

      // Find the free_surface location
      fs_geom_object_pt->position(initial_zeta, free_x);

      for (unsigned i = 0; i < 2; ++i)
      {
        dx[i] = spine_x[i] - free_x[i];
      }

      // Calculate the entries of the jacobian matrix
      jacobian(0, 0) = (spine_end[0] - spine_base[0]);
      jacobian(1, 0) = (spine_end[1] - spine_base[1]);

      // Calculate the others by finite differences
      double FD_Jstep = 1.0e-6;
      double old_zeta = initial_zeta[0];
      initial_zeta[0] += FD_Jstep;
      fs_geom_object_pt->position(initial_zeta, new_free_x);

      for (unsigned i = 0; i < 2; ++i)
      {
        jacobian(i, 1) = (free_x[i] - new_free_x[i]) / FD_Jstep;
      }

      // Now let's solve the damn thing
      jacobian.solve(dx);

      lambda -= dx[0];
      initial_zeta[0] = old_zeta - dx[1];

      // How are we doing

      for (unsigned i = 0; i < 2; ++i)
      {
        spine_x[i] = spine_base[i] + lambda * (spine_end[i] - spine_base[i]);
      }
      fs_geom_object_pt->position(initial_zeta, free_x);

      for (unsigned i = 0; i < 2; ++i)
      {
        dx[i] = spine_x[i] - free_x[i];
      }
      maxres = std::fabs(dx[0]) > std::fabs(dx[1]) ? std::fabs(dx[0]) :
                                                     std::fabs(dx[1]);

      if (count > 100)
      {
        throw OomphLibError("Failed to converge after 100 steps",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

    } while (maxres > 1.0e-8);


    oomph_info << "DONE " << initial_zeta[0] << std::endl;
    double spine_length = sqrt(pow((spine_base[0] - spine_end[0]), 2.0) +
                               pow((spine_base[1] - spine_end[1]), 2.0));

    return (lambda * spine_length); // Divided by spine length
  }

  //=======================================================================
  /// Reposition the spines that emenate from the lower wall
  //=======================================================================
  template<class ELEMENT, class INTERFACE_ELEMENT>
  void BrethertonSpineMesh<ELEMENT, INTERFACE_ELEMENT>::reposition_spines(
    const double& zeta_lo_transition_start,
    const double& zeta_lo_transition_end,
    const double& zeta_up_transition_start,
    const double& zeta_up_transition_end)
  {
    // Firstly create a geometric object of the free surface
    Mesh* fs_mesh_pt = new Mesh;
    this->template build_face_mesh<ELEMENT, FaceElementAsGeomObject>(
      4, fs_mesh_pt);

    // Loop over these new face elements and set the boundary number
    // of the bulk mesh
    unsigned n_face_element = fs_mesh_pt->nelement();
    // Loop over the elements
    for (unsigned e = 0; e < n_face_element; e++)
    {
      // Cast the element pointer to the correct thing!
      dynamic_cast<FaceElementAsGeomObject<ELEMENT>*>(fs_mesh_pt->element_pt(e))
        ->set_boundary_number_in_bulk_mesh(4);
    }

    // Now make a single geometric object that represents the collection of
    // geometric objects that form the boundary of the bulk mesh. Two
    // Eulerian coordinates, one intrinsic coordinate.
    MeshAsGeomObject* fs_geom_object_pt = new MeshAsGeomObject(fs_mesh_pt);


    // Length of deposited film region
    double llayer_lower = zeta_lo_transition_start - Zeta_start;
    double llayer_upper = zeta_up_transition_start - Zeta_start;

    // Length of transition region
    double d_lower = zeta_lo_transition_end - zeta_lo_transition_start;
    double d_upper = zeta_up_transition_end - zeta_up_transition_start;

    // Work out radius of circular cap from lower and upper wall
    Vector<double> r_wall_lo(2), r_wall_up(2);
    Vector<double> zeta(1), s_lo(1), s_up(1);
    GeomObject *lower_sub_geom_object_pt = 0, *upper_sub_geom_object_pt = 0;

    GeomObject* lower_transition_geom_object_pt = 0;
    GeomObject* upper_transition_geom_object_pt = 0;
    Vector<double> s_transition_lo(1), s_transition_up(1);
    Vector<double> spine_centre(2);

    // Get number of nodes along element edge
    unsigned np = this->finite_element_pt(0)->nnode_1d();

    // Calculate the centre for the spine nodes in the transition region
    {
      // Get the geometric objects on the walls at the start of the transition
      // region
      // Lower wall
      zeta[0] = zeta_lo_transition_start;
      Lower_wall_pt->locate_zeta(
        zeta, lower_transition_geom_object_pt, s_transition_lo);
      // Upper wall
      zeta[0] = zeta_up_transition_start;
      Upper_wall_pt->locate_zeta(
        zeta, upper_transition_geom_object_pt, s_transition_up);

      // Find the Eulerian coordinates of the walls at the transition region
      lower_transition_geom_object_pt->position(s_transition_lo, r_wall_lo);
      upper_transition_geom_object_pt->position(s_transition_up, r_wall_up);

      // Take the average of these positions to define the origin of the spines
      // in the transition region Horizontal position is always halfway
      spine_centre[0] = 0.5 * (r_wall_lo[0] + r_wall_up[0]);

      // Vertical Position is given by a specified fraction
      // between the upper and lower walls
      spine_centre[1] =
        r_wall_lo[1] + spine_centre_fraction() * (r_wall_up[1] - r_wall_lo[1]);
    }


    // Initialise number of elements in previous regions:
    unsigned n_prev_elements = 0;

    // Storage for the end of the spines
    Vector<double> spine_end(2);
    Vector<double> fs_zeta(1, 0.0);

    // Loop over elements in lower deposited film region
    // -------------------------------------------------
    {
      oomph_info << "LOWER FILM " << std::endl;
      // Increments in wall coordinate
      double dzeta_el = llayer_lower / double(Nx1);
      double dzeta_node = llayer_lower / double(Nx1 * (np - 1));

      // Loop over elements horizontally
      for (unsigned i = 0; i < Nx1; i++)
      {
        // Start of wall coordinate
        double zeta_lo = Zeta_start + double(i) * dzeta_el;

        // Work out element number in overall mesh
        unsigned e = n_prev_elements + i * (Nh + 1);

        // Get pointer to lower element
        FiniteElement* el_pt = this->finite_element_pt(e);

        // Loop over its nodes "horizontally"
        for (unsigned i1 = 0; i1 < (np - 1); i1++)
        {
          // Get spine node
          SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(i1));

          // Get the Lagrangian coordinate in the Lower Wall
          zeta[0] = zeta_lo + double(i1) * dzeta_node;
          // Reset the boundary coordinate
          nod_pt->set_coordinates_on_boundary(0, zeta);
          // Get the geometric object and local coordinate
          Lower_wall_pt->locate_zeta(zeta, lower_sub_geom_object_pt, s_lo);

          // The local coordinate is a geometric parameter
          // This needs to be set (rather than added) because the
          // same spine may be visited more than once
          Vector<double> parameters(1, s_lo[0]);
          nod_pt->spine_pt()->set_geom_parameter(parameters);

          // The sub geom object is one (and only) geom object
          // for spine:
          Vector<GeomObject*> geom_object_pt(1);
          geom_object_pt[0] = lower_sub_geom_object_pt;

          // Pass geom object(s) to spine
          nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

          // Get the wall position at the bottom of the spine
          lower_sub_geom_object_pt->position(s_lo, r_wall_lo);
          // The end of the spine is vertically above the base
          spine_end[0] = r_wall_lo[0];
          spine_end[1] = spine_centre[1];
          nod_pt->spine_pt()->height() = find_distance_to_free_surface(
            fs_geom_object_pt, fs_zeta, r_wall_lo, spine_end);
        }
      }

      // Increment number of previous elements
      n_prev_elements += Nx1 * (Nh + 1);
    }


    // Loop over elements in lower horizontal transition region
    // --------------------------------------------------------
    {
      oomph_info << "LOWER HORIZONTAL TRANSITION " << std::endl;
      // Increments in wall coordinate
      double dzeta_el = d_lower / double(Nx2);
      double dzeta_node = d_lower / double(Nx2 * (np - 1));

      // Loop over elements horizontally
      for (unsigned i = 0; i < Nx2; i++)
      {
        // Start of wall coordinate
        double zeta_lo = zeta_lo_transition_start + double(i) * dzeta_el;

        // Work out element number in overall mesh
        unsigned e = n_prev_elements + i * (Nh + 1);

        // Get pointer to element
        FiniteElement* el_pt = this->finite_element_pt(e);

        // Loop over its nodes
        for (unsigned i1 = 0; i1 < (np - 1); i1++)
        {
          // Get spine node
          SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(i1));

          // Get the Lagrangian coordinate in the Lower Wall
          zeta[0] = zeta_lo + double(i1) * dzeta_node;
          // Reset the boundary coordinate
          nod_pt->set_coordinates_on_boundary(0, zeta);
          // Get the sub geometric object and local coordinate
          Lower_wall_pt->locate_zeta(zeta, lower_sub_geom_object_pt, s_lo);

          // Pass geometric parameter to the spine
          Vector<double> parameters(3);
          parameters[0] = s_lo[0];
          parameters[1] = s_transition_lo[0];
          parameters[2] = s_transition_up[0];
          nod_pt->spine_pt()->set_geom_parameter(parameters);

          // Lower sub geom object is one (and only) geom object
          // for spine:
          Vector<GeomObject*> geom_object_pt(3);
          geom_object_pt[0] = lower_sub_geom_object_pt;
          geom_object_pt[1] = lower_transition_geom_object_pt;
          geom_object_pt[2] = upper_transition_geom_object_pt;

          // Pass geom object(s) to spine
          nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

          // Get position vector to wall
          lower_sub_geom_object_pt->position(s_lo, r_wall_lo);
          // The end of the spine is the spine centre,so the height is easy(ish)
          nod_pt->spine_pt()->height() = find_distance_to_free_surface(
            fs_geom_object_pt, fs_zeta, r_wall_lo, spine_centre);
        }
      }

      // Increment number of previous elements
      n_prev_elements += Nx2 * (Nh + 1);
    }

    // Loop over elements in lower vertical transition region
    // --------------------------------------------------------
    {
      oomph_info << "LOWER VERTICAL TRANSITION " << std::endl;
      for (unsigned i = 0; i < Nhalf; i++)
      {
        // Work out element number in overall mesh
        unsigned e = n_prev_elements + i * (Nh + 1);

        // Get pointer to element
        FiniteElement* el_pt = this->finite_element_pt(e);

        // Loop over its nodes
        for (unsigned i1 = 0; i1 < (np - 1); i1++)
        {
          // Get spine node
          // Note that I have to loop over the second row of nodes
          // because the first row are updated in region 6 and so
          // you get the wrong spines from them (doh!)
          SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(np + i1));

          // Get position vectors to wall
          zeta[0] = zeta_lo_transition_end;
          // Get the sub geometric objects
          Lower_wall_pt->locate_zeta(zeta, lower_sub_geom_object_pt, s_lo);
          zeta[0] = zeta_up_transition_end;
          Upper_wall_pt->locate_zeta(zeta, upper_sub_geom_object_pt, s_up);

          lower_sub_geom_object_pt->position(s_lo, r_wall_lo);
          upper_sub_geom_object_pt->position(s_up, r_wall_up);

          // Set vertical fraction
          double vertical_fraction =
            (double(i) + double(i1) / double(np - 1)) / double(2.0 * Nhalf);

          // Add the geometric parameters in order
          Vector<double> parameters(5);
          parameters[0] = s_lo[0];
          parameters[1] = s_up[0];
          parameters[2] = vertical_fraction;
          parameters[3] = s_transition_lo[0];
          parameters[4] = s_transition_up[0];
          nod_pt->spine_pt()->set_geom_parameter(parameters);

          // Origin of spine
          Vector<double> S0(2);
          S0[0] =
            r_wall_lo[0] + vertical_fraction * (r_wall_up[0] - r_wall_lo[0]);
          S0[1] =
            r_wall_lo[1] + vertical_fraction * (r_wall_up[1] - r_wall_lo[1]);

          // Lower and Upper wall sub geom objects affect  spine:
          Vector<GeomObject*> geom_object_pt(4);
          geom_object_pt[0] = lower_sub_geom_object_pt;
          geom_object_pt[1] = upper_sub_geom_object_pt;
          geom_object_pt[2] = lower_transition_geom_object_pt;
          geom_object_pt[3] = upper_transition_geom_object_pt;

          // Pass geom object(s) to spine
          nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

          // Calculate the spine height
          nod_pt->spine_pt()->height() = find_distance_to_free_surface(
            fs_geom_object_pt, fs_zeta, S0, spine_centre);
        }
      }

      // Increment number of previous elements
      n_prev_elements += Nhalf * (Nh + 1);
    }


    // Loop over elements in upper vertical transition region
    // --------------------------------------------------------
    {
      oomph_info << "UPPER VERTICAL TRANSITION" << std::endl;
      for (unsigned i = 0; i < Nhalf; i++)
      {
        // Work out element number in overall mesh
        unsigned e = n_prev_elements + i * (Nh + 1);

        // Get pointer to element
        FiniteElement* el_pt = this->finite_element_pt(e);

        // Loop over its nodes
        for (unsigned i1 = 0; i1 < (np - 1); i1++)
        {
          // Get spine node
          // Note that I have to loop over the second row of nodes
          // because the first row are updated in region 6 and so
          // you get the wrong spines from them (doh!)
          SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(np + i1));

          // Get position vectors to wall
          zeta[0] = zeta_lo_transition_end;
          // Get the sub geometric objects
          Lower_wall_pt->locate_zeta(zeta, lower_sub_geom_object_pt, s_lo);
          zeta[0] = zeta_up_transition_end;
          Upper_wall_pt->locate_zeta(zeta, upper_sub_geom_object_pt, s_up);

          lower_sub_geom_object_pt->position(s_lo, r_wall_lo);
          upper_sub_geom_object_pt->position(s_up, r_wall_up);


          // Set vertical fraction
          double vertical_fraction =
            0.5 +
            (double(i) + double(i1) / double(np - 1)) / double(2.0 * Nhalf);

          // Add the geometric parameters in order
          Vector<double> parameters(5);
          parameters[0] = s_lo[0];
          parameters[1] = s_up[0];
          parameters[2] = vertical_fraction;
          parameters[3] = s_transition_lo[0];
          parameters[4] = s_transition_up[0];
          nod_pt->spine_pt()->set_geom_parameter(parameters);

          // Origin of spine
          Vector<double> S0(2);
          S0[0] =
            r_wall_lo[0] + vertical_fraction * (r_wall_up[0] - r_wall_lo[0]);
          S0[1] =
            r_wall_lo[1] + vertical_fraction * (r_wall_up[1] - r_wall_lo[1]);

          // Lower and Upper wall sub geom objects affect  spine:
          Vector<GeomObject*> geom_object_pt(4);
          geom_object_pt[0] = lower_sub_geom_object_pt;
          geom_object_pt[1] = upper_sub_geom_object_pt;
          geom_object_pt[2] = lower_transition_geom_object_pt;
          geom_object_pt[3] = upper_transition_geom_object_pt;

          // Pass geom object(s) to spine
          nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

          // Calculate the spine height
          nod_pt->spine_pt()->height() = find_distance_to_free_surface(
            fs_geom_object_pt, fs_zeta, S0, spine_centre);
        }
      }

      // Increment number of previous elements
      n_prev_elements += Nhalf * (Nh + 1);
    }


    // Loop over elements in upper horizontal transition region
    // --------------------------------------------------------
    {
      oomph_info << "UPPER HORIZONTAL TRANSITION " << std::endl;
      // Increments in wall coordinate
      double dzeta_el = d_upper / double(Nx2);
      double dzeta_node = d_upper / double(Nx2 * (np - 1));

      // Loop over elements horizontally
      for (unsigned i = 0; i < Nx2; i++)
      {
        // Start of wall coordinate
        double zeta_lo = zeta_up_transition_end - double(i) * dzeta_el;

        // Work out element number in overall mesh
        unsigned e = n_prev_elements + i * (Nh + 1);

        // Get pointer to element
        FiniteElement* el_pt = this->finite_element_pt(e);

        // Loop over its nodes
        for (unsigned i1 = 0; i1 < (np - 1); i1++)
        {
          // Get spine node (same comment)
          SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(np + i1));

          // Get the Lagrangian coordinate in the Lower Wall
          zeta[0] = zeta_lo - double(i1) * dzeta_node;
          // Reset the boundary coordinate
          el_pt->node_pt(i1)->set_coordinates_on_boundary(2, zeta);
          // Get the sub geometric object and local coordinate
          Upper_wall_pt->locate_zeta(zeta, upper_sub_geom_object_pt, s_up);

          // Pass geometric parameter to the spine
          Vector<double> parameters(3);
          parameters[0] = s_up[0];
          parameters[1] = s_transition_lo[0];
          parameters[2] = s_transition_up[0];
          nod_pt->spine_pt()->set_geom_parameter(parameters);

          // Lower sub geom object is one (and only) geom object
          // for spine:
          Vector<GeomObject*> geom_object_pt(3);
          geom_object_pt[0] = upper_sub_geom_object_pt;
          geom_object_pt[1] = lower_transition_geom_object_pt;
          geom_object_pt[2] = upper_transition_geom_object_pt;

          // Pass geom object(s) to spine
          nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);


          // Get position vector to wall
          upper_sub_geom_object_pt->position(s_up, r_wall_up);
          // Find spine height
          nod_pt->spine_pt()->height() = find_distance_to_free_surface(
            fs_geom_object_pt, fs_zeta, r_wall_up, spine_centre);
        }
      }

      // Increment number of previous elements
      n_prev_elements += Nx2 * (Nh + 1);
    }


    // Loop over elements in upper deposited film region
    // -------------------------------------------------
    {
      oomph_info << "UPPER THIN FILM" << std::endl;
      // Increments in wall coordinate
      double dzeta_el = llayer_upper / double(Nx1);
      double dzeta_node = llayer_upper / double(Nx1 * (np - 1));

      // Loop over elements horizontally
      for (unsigned i = 0; i < Nx1; i++)
      {
        // Start of wall coordinate
        double zeta_lo = zeta_up_transition_start - double(i) * dzeta_el;

        // Work out element number in overall mesh
        unsigned e = n_prev_elements + i * (Nh + 1);

        // Get pointer to element
        FiniteElement* el_pt = this->finite_element_pt(e);

        // Loop over its nodes "horizontally"
        for (unsigned i1 = 0; i1 < (np - 1); i1++)
        {
          // Get spine node
          SpineNode* nod_pt = dynamic_cast<SpineNode*>(el_pt->node_pt(i1));

          // Get the Lagrangian coordinate in the Upper wall
          zeta[0] = zeta_lo - double(i1) * dzeta_node;
          // Reset coordinate on boundary
          nod_pt->set_coordinates_on_boundary(2, zeta);
          // Get the geometric object and local coordinate
          Upper_wall_pt->locate_zeta(zeta, upper_sub_geom_object_pt, s_up);

          // The local coordinate is a geometric parameter
          Vector<double> parameters(1, s_up[0]);
          nod_pt->spine_pt()->set_geom_parameter(parameters);

          // upper sub geom object is one (and only) geom object
          // for spine:
          Vector<GeomObject*> geom_object_pt(1);
          geom_object_pt[0] = upper_sub_geom_object_pt;

          // Pass geom object(s) to spine
          nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

          // Get the wall position at the bottom of the spine
          upper_sub_geom_object_pt->position(s_up, r_wall_up);
          spine_end[0] = r_wall_up[0];
          spine_end[1] = spine_centre[1];
          // Find the new spine height
          nod_pt->spine_pt()->height() = find_distance_to_free_surface(
            fs_geom_object_pt, fs_zeta, r_wall_up, spine_end);
        }
      }


      // Increment number of previous elements
      n_prev_elements += Nx1 * (Nh + 1);
    }


    // Additional mesh
    {
      unsigned e = n_prev_elements;

      // Get pointer to node
      SpineNode* nod_pt = this->element_node_pt(e, 0);

      // Need to get the local coordinates for the upper and lower wall
      zeta[0] = zeta_lo_transition_end;
      // Get the sub geometric objects
      Lower_wall_pt->locate_zeta(zeta, lower_sub_geom_object_pt, s_lo);
      zeta[0] = zeta_up_transition_end;
      Upper_wall_pt->locate_zeta(zeta, upper_sub_geom_object_pt, s_up);

      lower_sub_geom_object_pt->position(s_lo, r_wall_lo);
      upper_sub_geom_object_pt->position(s_up, r_wall_up);

      // Pass additional data to spine
      Vector<double> parameters(2);
      parameters[0] = s_lo[0];
      parameters[1] = s_up[0];
      nod_pt->spine_pt()->set_geom_parameter(parameters);

      // Lower and upper wall sub geom objects affect update
      // for spine:
      Vector<GeomObject*> geom_object_pt(2);
      geom_object_pt[0] = lower_sub_geom_object_pt;
      geom_object_pt[1] = upper_sub_geom_object_pt;

      // Pass geom object(s) to spine
      nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);

      // LOOP OVER OTHER SPINES
      //----------------------

      // Now loop over the elements horizontally
      for (unsigned long j = 0; j < Nx3; j++)
      {
        unsigned e = n_prev_elements + j;

        // Loop over the nodes in the elements horizontally, ignoring
        // the first column
        for (unsigned l2 = 0; l2 < np; l2++)
        {
          // Get pointer to node at the base of the spine
          SpineNode* nod_pt = this->element_node_pt(e, l2);

          // Increment in wall coordinate
          double dzeta_el_lower =
            (Zeta_end - zeta_lo_transition_end) / double(Nx3);
          double dzeta_nod_lower = dzeta_el_lower / double(np - 1);

          double dzeta_el_upper =
            (Zeta_end - zeta_up_transition_end) / double(Nx3);
          double dzeta_nod_upper = dzeta_el_upper / double(np - 1);

          // Get wall coordinate
          zeta[0] =
            zeta_lo_transition_end + j * dzeta_el_lower + l2 * dzeta_nod_lower;
          // Reset the boundary coordinate
          nod_pt->set_coordinates_on_boundary(0, zeta);

          // Get the sub geometric objects
          Lower_wall_pt->locate_zeta(zeta, lower_sub_geom_object_pt, s_lo);

          zeta[0] =
            zeta_up_transition_end + j * dzeta_el_upper + l2 * dzeta_nod_upper;
          // Reset the upper boundary coordinate
          this->element_node_pt(e + Nx3 * (2 * Nhalf - 1), np * (np - 1) + l2)
            ->set_coordinates_on_boundary(2, zeta);
          Upper_wall_pt->locate_zeta(zeta, upper_sub_geom_object_pt, s_up);

          lower_sub_geom_object_pt->position(s_lo, r_wall_lo);
          upper_sub_geom_object_pt->position(s_up, r_wall_up);

          // Add geometric parameters to spine
          Vector<double> parameters(2);
          parameters[0] = s_lo[0];
          parameters[1] = s_up[0];
          nod_pt->spine_pt()->set_geom_parameter(parameters);

          // Lower and upper sub geom objects affect update
          // for spine:
          Vector<GeomObject*> geom_object_pt(2);
          geom_object_pt[0] = lower_sub_geom_object_pt;
          geom_object_pt[1] = upper_sub_geom_object_pt;

          // Pass geom object(s) to spine
          nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);
        }
      }
    }


    // Clean up all the memory
    // Can delete the Geometric object
    delete fs_geom_object_pt;
    // Need to be careful with the FaceMesh, because we can't delete the nodes
    // Loop
    for (unsigned e = n_face_element; e > 0; e--)
    {
      delete fs_mesh_pt->element_pt(e - 1);
      fs_mesh_pt->element_pt(e - 1) = 0;
    }
    // Now clear all element and node storage (should maybe fine-grain this)
    fs_mesh_pt->flush_element_and_node_storage();
    // Finally delete the mesh
    delete fs_mesh_pt;
  }

} // namespace oomph

#endif
