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
// Non-inline functions for Rigid Body Elements
#include "immersed_rigid_body_elements.h"

namespace oomph
{
  /// Static default value for physical constants
  /// Zero gives instantaneous force and torque balances --- no solid intertia
  double ImmersedRigidBodyElement::Default_Physical_Constant_Value = 0.0;

  /// Static default value for physical ratio
  double ImmersedRigidBodyElement::Default_Physical_Ratio_Value = 1.0;

  /// Static default gravity direction vector
  Vector<double> ImmersedRigidBodyElement::Default_Gravity_vector(2, 0.0);

  //=======================================================================
  /// Work out the position derivative taking into account the movement
  /// relative to the original centre of mass
  //======================================================================
  void ImmersedRigidBodyElement::dposition_dt(const Vector<double>& zeta,
                                              const unsigned& j,
                                              Vector<double>& drdt)
  {
    // Switch on time level
    switch (j)
    {
        // Current time, just return the position
      case 0:
        return position(zeta, drdt);
        break;

        // Time derivative
      case 1:
      {
        // Get the initial position of the underlying geometric object
        Vector<double> initial_x(2);
        Geom_object_pt->position(zeta, initial_x);
        // Scale relative to the centre of mass
        double X = initial_x[0] - Initial_centre_of_mass[0];
        double Y = initial_x[1] - Initial_centre_of_mass[1];

        // Now calculate the original angle and radius
        double phi_orig = atan2(Y, X);
        double r_orig = sqrt(X * X + Y * Y);

        // Get first time derivatives of all displacement data
        Vector<double> veloc(3);
        // Get the velocity of the data
        this->Centre_displacement_data_pt->time_stepper_pt()->time_derivative(
          1, this->Centre_displacement_data_pt, veloc);

        // Now use the chain rule to specify the boundary velocities
        drdt[0] = veloc[0];
        drdt[1] = veloc[1];

        if (Include_geometric_rotation)
        {
          drdt[0] +=
            -r_orig *
            sin(phi_orig + this->Centre_displacement_data_pt->value(2)) *
            veloc[2];
          drdt[1] +=
            r_orig *
            cos(phi_orig + this->Centre_displacement_data_pt->value(2)) *
            veloc[2];
        }
      }
        // Done
        return;
        break;

      default:
        std::ostringstream warning_stream;
        warning_stream
          << "Using default (static) assignment " << j
          << "-th time derivative in GeomObject::dposition_dt(...) is zero\n"
          << "Overload for your specific geometric object if this is not \n"
          << "appropriate. \n";
        OomphLibWarning(warning_stream.str(),
                        "GeomObject::dposition_dt()",
                        OOMPH_EXCEPTION_LOCATION);

        unsigned n = drdt.size();
        for (unsigned i = 0; i < n; i++)
        {
          drdt[i] = 0.0;
        }
        break;
    }
  }

  //=======================================================================
  /// Output the position of the centre of gravity including velocities
  /// and accelerations
  //======================================================================
  void ImmersedRigidBodyElement::output_centre_of_gravity(std::ostream& outfile)
  {
    // Get timestepper
    TimeStepper* time_stepper_pt =
      this->Centre_displacement_data_pt->time_stepper_pt();

    // Get first time derivatives of all displacement data
    Vector<double> veloc(3);
    time_stepper_pt->time_derivative(
      1, this->Centre_displacement_data_pt, veloc);

    // Get second time derivatives of all displacement data
    Vector<double> accel(3);
    time_stepper_pt->time_derivative(
      2, this->Centre_displacement_data_pt, accel);

    outfile << time_stepper_pt->time() << " "
            << Initial_centre_of_mass[0] +
                 this->Centre_displacement_data_pt->value(0)
            << " "
            << Initial_centre_of_mass[1] +
                 this->Centre_displacement_data_pt->value(1)
            << " " << Initial_Phi + this->Centre_displacement_data_pt->value(2)
            << " " << veloc[0] << " " << veloc[1] << " " << veloc[2] << " "
            << accel[0] << " " << accel[1] << " " << accel[2] << std::endl;
  }


  //======================================================================
  /// Obtain the external force and torque on the body from specified
  /// function pointers and also from a drag mesh, if there is one
  //=====================================================================
  void ImmersedRigidBodyElement::get_force_and_torque(const double& time,
                                                      Vector<double>& force,
                                                      double& torque)
  {
    // Get external force
    if (External_force_fct_pt == 0)
    {
      force[0] = 0.0;
      force[1] = 0.0;
    }
    else
    {
      External_force_fct_pt(time, force);
    }

    // Get external torque
    if (External_torque_fct_pt == 0)
    {
      torque = 0.0;
    }
    else
    {
      External_torque_fct_pt(time, torque);
    }

    // Add drag from any (fluid) mesh attached to surface of body
    Vector<double> element_drag_force(2);
    Vector<double> element_drag_torque(1);
    if (Drag_mesh_pt == 0)
    {
      return;
    }
    else
    {
      unsigned nel = Drag_mesh_pt->nelement();

      for (unsigned e = 0; e < nel; e++)
      {
        dynamic_cast<ElementWithDragFunction*>(Drag_mesh_pt->element_pt(e))
          ->get_drag_and_torque(element_drag_force, element_drag_torque);
        force[0] += element_drag_force[0];
        force[1] += element_drag_force[1];
        torque += element_drag_torque[0];
      }
    }
  }

  //=======================================================================
  /// Set the external drag mesh, which should consist of
  /// NavierStokesSurfaceDragTorqueElements and then read out the
  /// appropriate load and geometric data from those elements and set
  /// as external data for this element.
  //=======================================================================
  void ImmersedRigidBodyElement::set_drag_mesh(Mesh* const& drag_mesh_pt)
  {
    // Delete the external hijacked data
    this->delete_external_hijacked_data();
    // Flush any existing external data
    this->flush_external_data();
    // Set the pointer
    Drag_mesh_pt = drag_mesh_pt;

    // Allocate storage for all geometric data in the mesh
    std::set<Data*> bulk_geometric_data_pt;
    // Allocate storage for all load data in the mesh
    std::set<std::pair<Data*, unsigned>> bulk_load_data_pt;

    // Loop over all elements in the drag mesh
    const unsigned n_element = drag_mesh_pt->nelement();
    for (unsigned e = 0; e < n_element; ++e)
    {
      // Cast the bulk element associated with each FaceElement to
      // an FSIFluidElement
      FSIFluidElement* bulk_elem_pt = dynamic_cast<FSIFluidElement*>(
        dynamic_cast<FaceElement*>(drag_mesh_pt->element_pt(e))
          ->bulk_element_pt());
      // Check that the cast worked
      if (bulk_elem_pt == 0)
      {
        throw OomphLibError("Drag mesh must consist of FSIFluidElements\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Add the geometric and load data from the bulk element to the
      // set allocated above
      bulk_elem_pt->identify_geometric_data(bulk_geometric_data_pt);
      bulk_elem_pt->identify_load_data(bulk_load_data_pt);
    }

    // Need to add all these data as external data to this (Rigid Body) object
    for (std::set<Data*>::iterator it = bulk_geometric_data_pt.begin();
         it != bulk_geometric_data_pt.end();
         ++it)
    {
      this->add_external_data(*it);
    }

    // Now do the same but make custom data for the load data
    for (std::set<std::pair<Data*, unsigned>>::iterator it =
           bulk_load_data_pt.begin();
         it != bulk_load_data_pt.end();
         ++it)
    {
      Data* temp_data_pt = new HijackedData(it->second, it->first);
      List_of_external_hijacked_data.push_back(
        this->add_external_data(temp_data_pt));
    }
  }

  //======================================================================
  /// Initialise the internal data
  //=====================================================================
  void ImmersedRigidBodyElement::initialise(TimeStepper* const& time_stepper_pt)
  {
    // This could be calculated by an integral around the boundary
    Initial_centre_of_mass.resize(2, 0.0);

    // Temporary hack
    if (time_stepper_pt == 0)
    {
      return;
    }

    // Provide Data for centre-of-mass displacement internally
    if (Centre_displacement_data_pt == 0)
    {
      Centre_displacement_data_pt = new Data(time_stepper_pt, 3);

      // I've created it so I have to tidy up too!
      Displacement_data_is_internal = true;

      // Centre displacement is internal Data for this element
      Index_for_centre_displacement =
        add_internal_data(Centre_displacement_data_pt);
    }
    // Data created externally, so somebody else will clean up
    else
    {
      Displacement_data_is_internal = false;

      // Centre displacement is external Data for this element
      Index_for_centre_displacement =
        add_external_data(Centre_displacement_data_pt);
    }
  }


  //=======================================================================
  /// Calculate the contributions to the residuals and the jacobian
  //======================================================================
  void ImmersedRigidBodyElement::get_residuals_rigid_body_generic(
    Vector<double>& residuals, DenseMatrix<double>& jacobian, const bool& flag)
  {
    // Get timestepper and time
    TimeStepper* timestepper_pt =
      this->Centre_displacement_data_pt->time_stepper_pt();
    double time = timestepper_pt->time();

    // Get second time derivatives of all displacement data
    Vector<double> accel(3);
    timestepper_pt->time_derivative(
      2, this->Centre_displacement_data_pt, accel);

    // Get force and torque
    Vector<double> external_force(2);
    double external_torque;
    get_force_and_torque(time, external_force, external_torque);

    // Get the timescale ratio product of Reynolds number
    // and Strouhal number squared
    const double Lambda_sq =
      this->re() * this->st() * this->st() * this->density_ratio();

    // Get the effective ReInvFr which must be multiplied by the
    // density ratio to compute the gravitational load on the rigid body
    const double scaled_re_inv_fr = this->re_invfr() * this->density_ratio();
    // Get the gravitational load
    Vector<double> G = g();

    // Newton's law
    int local_eqn = 0;
    local_eqn = this->centre_displacement_local_eqn(0);
    if (local_eqn >= 0)
    {
      residuals[local_eqn] = Lambda_sq * Mass * accel[0] - external_force[0] -
                             Mass * scaled_re_inv_fr * G[0];

      // Get Jacobian too?
      if (flag)
      {
        jacobian(local_eqn, local_eqn) =
          Lambda_sq * Mass * timestepper_pt->weight(2, 0);
      }
    }

    local_eqn = this->centre_displacement_local_eqn(1);
    if (local_eqn >= 0)
    {
      residuals[local_eqn] = Lambda_sq * Mass * accel[1] - external_force[1] -
                             Mass * scaled_re_inv_fr * G[1];
      // Get Jacobian too?
      if (flag)
      {
        jacobian(local_eqn, local_eqn) =
          Lambda_sq * Mass * timestepper_pt->weight(2, 0);
      }
    }

    local_eqn = this->centre_displacement_local_eqn(2);
    if (local_eqn >= 0)
    {
      residuals[local_eqn] =
        Lambda_sq * Moment_of_inertia * accel[2] - external_torque;
      // Get Jacobian too?
      if (flag)
      {
        jacobian(local_eqn, local_eqn) =
          Lambda_sq * Moment_of_inertia * timestepper_pt->weight(2, 0);
      }
    }
  }


  //=======================================================================
  /// Constructor: Specify coordinates of a point inside the hole
  /// and a vector of pointers to TriangleMeshPolyLines
  /// that define the boundary segments of the polygon.
  /// Each TriangleMeshPolyLine has its own boundary ID and can contain
  /// multiple (straight-line) segments. The optional final argument
  /// is a pointer to a Data object whose three values represent
  /// the two displacements of and the rotation angle about the polygon's
  /// centre of mass.
  //=======================================================================
  ImmersedRigidBodyTriangleMeshPolygon::ImmersedRigidBodyTriangleMeshPolygon(
    const Vector<double>& hole_center,
    const Vector<TriangleMeshCurveSection*>& boundary_polyline_pt,
    TimeStepper* const& time_stepper_pt,
    Data* const& centre_displacement_data_pt)
    : TriangleMeshCurve(boundary_polyline_pt),
      TriangleMeshClosedCurve(boundary_polyline_pt, hole_center),
      TriangleMeshPolygon(boundary_polyline_pt, hole_center),
      ImmersedRigidBodyElement(time_stepper_pt, centre_displacement_data_pt)
  {
    // The underlying geometric object can be used to update the configuration
    // internally before a remesh
    this->Can_update_configuration = true;

    // Original rotation angle is zero
    Initial_Phi = 0.0;

    // Compute coordinates of centre of gravity etc
    Vector<double> r_left(2);
    Vector<double> r_right(2);
    Mass = 0.0;
    Initial_centre_of_mass[0] = 0.0;
    Initial_centre_of_mass[1] = 0.0;
    double inertia_x = 0.0;
    double inertia_y = 0.0;

    // Loop over polylines
    unsigned nboundary = boundary_polyline_pt.size();
    for (unsigned i = 0; i < nboundary; i++)
    {
      // Loop over the segments to get the vertex coordinates
      unsigned nseg = boundary_polyline_pt[i]->nsegment();
      for (unsigned j = 0; j < nseg; j++)
      {
        // Get the vertex coordinates
        r_left = this->polyline_pt(i)->vertex_coordinate(j);
        r_right = this->polyline_pt(i)->vertex_coordinate(j + 1);

        // Mass (area)
        Mass += 0.5 * (r_left[0] * r_right[1] - r_right[0] * r_left[1]);

        // Centroid
        Initial_centre_of_mass[0] +=
          (r_left[0] + r_right[0]) *
          (r_left[0] * r_right[1] - r_right[0] * r_left[1]);
        Initial_centre_of_mass[1] +=
          (r_left[1] + r_right[1]) *
          (r_left[0] * r_right[1] - r_right[0] * r_left[1]);
      }
      if (nboundary == 1)
      {
        // Get the vertex coordinates
        r_left = this->polyline_pt(0)->vertex_coordinate(nseg);
        r_right = this->polyline_pt(0)->vertex_coordinate(0);

        // Mass (area)
        Mass += 0.5 * (r_left[0] * r_right[1] - r_right[0] * r_left[1]);

        // Centroid
        Initial_centre_of_mass[0] +=
          (r_left[0] + r_right[0]) *
          (r_left[0] * r_right[1] - r_right[0] * r_left[1]);
        Initial_centre_of_mass[1] +=
          (r_left[1] + r_right[1]) *
          (r_left[0] * r_right[1] - r_right[0] * r_left[1]);
      }
    }

    // Normalise
    Initial_centre_of_mass[0] /= (6.0 * Mass);
    Initial_centre_of_mass[1] /= (6.0 * Mass);

    // Another loop over polylines for moment of inertia
    for (unsigned i = 0; i < nboundary; i++)
    {
      // Loop over the segments to get the vertex coordinates
      unsigned nseg = boundary_polyline_pt[i]->nsegment();
      for (unsigned j = 0; j < nseg; j++)
      {
        // Get the vertex coordinates
        r_left = this->polyline_pt(i)->vertex_coordinate(j);
        r_right = this->polyline_pt(i)->vertex_coordinate(j + 1);

        // Get moment about centroid
        r_left[0] -= Initial_centre_of_mass[0];
        r_left[1] -= Initial_centre_of_mass[1];
        r_right[0] -= Initial_centre_of_mass[0];
        r_right[1] -= Initial_centre_of_mass[1];

        // Moment of inertia
        inertia_x += 1.0 / 12.0 *
                     (r_left[1] * r_left[1] + r_left[1] * r_right[1] +
                      r_right[1] * r_right[1]) *
                     (r_left[0] * r_right[1] - r_right[0] * r_left[1]);

        inertia_y += 1.0 / 12.0 *
                     (r_left[0] * r_left[0] + r_left[0] * r_right[0] +
                      r_right[0] * r_right[0]) *
                     (r_left[0] * r_right[1] - r_right[0] * r_left[1]);
      }

      if (nboundary == 1)
      {
        // Get the vertex coordinates
        r_left = this->polyline_pt(0)->vertex_coordinate(nseg);
        r_right = this->polyline_pt(0)->vertex_coordinate(0);

        // Get moment about centroid
        r_left[0] -= Initial_centre_of_mass[0];
        r_left[1] -= Initial_centre_of_mass[1];
        r_right[0] -= Initial_centre_of_mass[0];
        r_right[1] -= Initial_centre_of_mass[1];

        // Moment of inertia
        inertia_x += 1.0 / 12.0 *
                     (r_left[1] * r_left[1] + r_left[1] * r_right[1] +
                      r_right[1] * r_right[1]) *
                     (r_left[0] * r_right[1] - r_right[0] * r_left[1]);

        inertia_y += 1.0 / 12.0 *
                     (r_left[0] * r_left[0] + r_left[0] * r_right[0] +
                      r_right[0] * r_right[0]) *
                     (r_left[0] * r_right[1] - r_right[0] * r_left[1]);
      }
    }

    // Polar moment of inertia is sum of two orthogonal planar moments
    Moment_of_inertia = inertia_x + inertia_y;

    //    // Tested for circular and elliptical cross section
    //    cout << "Mass             : " << Mass << std::endl;
    //    cout << "Moment of inertia: " << Moment_of_inertia << std::endl;
    //    cout << "X_c              : " << Initial_centre_of_mass[0] <<
    //    std::endl; cout << "Y_c              : " << Initial_centre_of_mass[1]
    //    << std::endl; pause("done");


    // Assign the intrinsic coordinate
    this->assign_zeta();

    //  {
    //   unsigned n_poly = this->npolyline();
    //   for(unsigned p=0;p<n_poly;++p)
    //    {
    //     std::cout << "Polyline " << p << "\n";
    //     std::cout << "-----------------------\n";
    //     unsigned n_vertex = Zeta_vertex[p].size();
    //     for(unsigned v=0;v<n_vertex;v++)
    //      {
    //       std::cout << v << " " << Zeta_vertex[p][v] << "\n";
    //      }
    //    }
    // }
  }


  //===============================================================
  /// Update the reference configuration by re-setting the original
  /// position of the vertices to their current ones, re-set the
  /// original position of the centre of mass, and the displacements
  /// and rotations relative to it
  //===============================================================
  void ImmersedRigidBodyTriangleMeshPolygon::reset_reference_configuration()
  {
    Vector<double> x_orig(2);
    Vector<double> r(2);

    // Loop over the polylines and update their vertex positions
    unsigned npoly = this->ncurve_section();
    for (unsigned i = 0; i < npoly; i++)
    {
      TriangleMeshPolyLine* poly_line_pt = this->polyline_pt(i);
      unsigned nvertex = poly_line_pt->nvertex();
      for (unsigned j = 0; j < nvertex; j++)
      {
        x_orig = poly_line_pt->vertex_coordinate(j);
        this->apply_rigid_body_motion(0, x_orig, r);
        poly_line_pt->vertex_coordinate(j) = r;
      }
    }

    // Update coordinates of hole
    Vector<double> orig_hole_coord(this->internal_point());
    this->apply_rigid_body_motion(0, orig_hole_coord, this->internal_point());

    // Update centre of gravity
    double x_displ = Centre_displacement_data_pt->value(0);
    double y_displ = Centre_displacement_data_pt->value(1);
    double phi_displ = Centre_displacement_data_pt->value(2);
    Initial_centre_of_mass[0] += x_displ;
    Initial_centre_of_mass[1] += y_displ;
    Initial_Phi += phi_displ;

    // Reset displacement and rotation ("non-previous-value"
    // history values stay)
    TimeStepper* timestepper_pt =
      Centre_displacement_data_pt->time_stepper_pt();
    unsigned nprev = timestepper_pt->nprev_values();
    for (unsigned t = 0; t < nprev; t++)
    {
      Centre_displacement_data_pt->set_value(
        t, 0, Centre_displacement_data_pt->value(t, 0) - x_displ);

      Centre_displacement_data_pt->set_value(
        t, 1, Centre_displacement_data_pt->value(t, 1) - y_displ);

      Centre_displacement_data_pt->set_value(
        t, 2, Centre_displacement_data_pt->value(t, 2) - phi_displ);
    }
  }

} // namespace oomph
