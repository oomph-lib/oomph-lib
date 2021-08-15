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
// Add in the header file
#include "rectangle_with_moving_cylinder_mesh.template.h"

// PML mesh helpers
#include "../generic/pml_meshes.h"

// Namespace extension
namespace oomph
{
  //=============================================================================
  /// Rectangular domain with circular whole
  //=============================================================================
  void RectangleWithHoleAndAnnularRegionDomain::
    project_point_on_cylinder_to_annular_boundary(const unsigned& time,
                                                  const Vector<double>& xi,
                                                  Vector<double>& r)
  {
    // Lagrangian coordinate corresponds to theta=pi on the circle
    Vector<double> xi_left(1, 4.0 * atan(1.0));

    // Point on circle corresponding to theta=pi
    Vector<double> point_left(2, 0.0);

    // Get the position of the point on the circle
    Cylinder_pt->position(time, xi_left, point_left);

    // Lagrangian coordinate corresponds to theta=0 on the circle
    Vector<double> xi_right(1, 0.0);

    // Point on circle corresponding to theta=0
    Vector<double> point_right(2, 0.0);

    // Get the position of the point on the circle
    Cylinder_pt->position(time, xi_right, point_right);

    // Storage for the coordinates of the centre
    Vector<double> point_centre(2, 0.0);

    // Loop over the coordinates
    for (unsigned i = 0; i < 2; i++)
    {
      // Calculate the i-th coordinate of the central point
      point_centre[i] = 0.5 * (point_left[i] + point_right[i]);
    }

    // Calculate the x-coordinate of the projectd point
    r[0] = point_centre[0] + Annular_region_radius * cos(xi[0]);

    // Calculate the y-coordinate of the projectd point
    r[1] = point_centre[1] + Annular_region_radius * sin(xi[0]);
  } // End of project_point_on_cylinder_to_annular_boundary


  /// \short Helper function that, given the Lagrangian coordinate, xi,
  /// (associated with a point on the cylinder), returns the corresponding
  /// point on the outer boundary of the annular region (where the inner
  /// boundary is prescribed by the boundary of the cylinder)
  void RectangleWithHoleAndAnnularRegionDomain::
    project_point_on_cylinder_to_annular_boundary(const double& time,
                                                  const Vector<double>& xi,
                                                  Vector<double>& r)
  {
    // Lagrangian coordinate corresponds to theta=pi on the circle
    Vector<double> xi_left(1, 4.0 * atan(1.0));

    // Point on circle corresponding to theta=pi
    Vector<double> point_left(2, 0.0);

    // Get the position of the point on the circle
    Cylinder_pt->position(time, xi_left, point_left);

    // Lagrangian coordinate corresponds to theta=0 on the circle
    Vector<double> xi_right(1, 0.0);

    // Point on circle corresponding to theta=0
    Vector<double> point_right(2, 0.0);

    // Get the position of the point on the circle
    Cylinder_pt->position(time, xi_right, point_right);

    // Storage for the coordinates of the centre
    Vector<double> point_centre(2, 0.0);

    // Loop over the coordinates
    for (unsigned i = 0; i < 2; i++)
    {
      // Calculate the i-th coordinate of the central point
      point_centre[i] = 0.5 * (point_left[i] + point_right[i]);
    }

    // Calculate the x-coordinate of the projectd point
    r[0] = point_centre[0] + Annular_region_radius * cos(xi[0]);

    // Calculate the y-coordinate of the projectd point
    r[1] = point_centre[1] + Annular_region_radius * sin(xi[0]);
  } // End of project_point_on_cylinder_to_annular_boundary


  /// \short Parametrisation of macro element boundaries: f(s) is the position
  /// vector to macro-element m's boundary in the specified direction [N/S/E/W]
  /// at the specified discrete time level (time=0: present; time>0: previous)
  void RectangleWithHoleAndAnnularRegionDomain::macro_element_boundary(
    const double& time,
    const unsigned& m,
    const unsigned& direction,
    const Vector<double>& s,
    Vector<double>& f)
  {
#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Create warning about time argument being moved to the front
    std::string error_message_string =
      "Order of function arguments has changed ";

    // Finish the string off
    error_message_string += "between versions 0.8 and 0.85";

    // Output a warning
    OomphLibWarning(
      error_message_string, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif

    // Lagrangian coordinate along surface of cylinder
    Vector<double> xi(1, 0.0);

    // Point on circle
    Vector<double> point_on_circle(2, 0.0);

    // Point on the outer boundary of the annular region
    Vector<double> point_on_annular_ring(2, 0.0);

    // Use the QuadTree enumeration for face directions
    using namespace QuadTreeNames;

    // Switch on the macro element
    switch (m)
    {
      // Macro element 0 is immediately to the left of the cylinder, outside
      // the inner annular region
      case 0:

        // Switch on the direction
        switch (direction)
        {
          case N:
            xi[0] = 3.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(Upper_mid_left, point_on_annular_ring, s[0], f);
            break;

          case S:
            xi[0] = -3.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(Lower_mid_left, point_on_annular_ring, s[0], f);
            break;

          case W:
            linear_interpolate(Lower_mid_left, Upper_mid_left, s[0], f);
            break;

          case E:
            xi[0] = 5.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 1 is immediately above the cylinder, outside
      // the inner annular region
      case 1:

        switch (direction)
        {
          case N:
            linear_interpolate(Upper_mid_left, Upper_mid_right, s[0], f);
            break;

          case S:
            xi[0] = 3.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case W:
            xi[0] = 3.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, Upper_mid_left, s[0], f);
            break;

          case E:
            xi[0] = 1.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, Upper_mid_right, s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect:  " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 2 is immediately to the right of the cylinder, outside
      // the inner annular region
      case 2:

        switch (direction)
        {
          case N:
            xi[0] = 1.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, Upper_mid_right, s[0], f);
            break;

          case S:
            xi[0] = -1.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, Lower_mid_right, s[0], f);
            break;

          case W:
            xi[0] = -atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case E:
            linear_interpolate(Lower_mid_right, Upper_mid_right, s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect:  " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 3 is immediately below cylinder, outside
      // the inner annular region
      case 3:

        switch (direction)
        {
          case N:
            xi[0] = -3.0 * atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case S:
            linear_interpolate(Lower_mid_left, Lower_mid_right, s[0], f);
            break;

          case W:
            xi[0] = -3.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(Lower_mid_left, point_on_annular_ring, s[0], f);
            break;

          case E:
            xi[0] = -1.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(Lower_mid_right, point_on_annular_ring, s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 4, is immediately to the left of the cylinder
      // lying in the inner annular region
      case 4:

        // Switch on the face of the m-th macro element
        switch (direction)
        {
          case N:
            xi[0] = 3.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, point_on_circle, s[0], f);
            break;

          case S:
            xi[0] = -3.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, point_on_circle, s[0], f);
            break;

          case W:
            xi[0] = 5.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case E:
            xi[0] = 5.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 5, is immediately above the cylinder in the
      // inner annular region
      case 5:

        switch (direction)
        {
          case N:
            xi[0] = 3.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case S:
            xi[0] = 3.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case W:
            xi[0] = 3.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_circle, point_on_annular_ring, s[0], f);
            break;

          case E:
            xi[0] = 1.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_circle, point_on_annular_ring, s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 6, is immediately to the right of the cylinder
      // lying in the inner annular region
      case 6:

        switch (direction)
        {
          case N:
            xi[0] = 1.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_circle, point_on_annular_ring, s[0], f);
            break;

          case S:
            xi[0] = -1.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_circle, point_on_annular_ring, s[0], f);
            break;

          case W:
            xi[0] = -atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case E:
            xi[0] = -atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 7, is immediately below the cylinder in the
      // inner annular region
      case 7:

        switch (direction)
        {
          case N:
            xi[0] = -3.0 * atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case S:
            xi[0] = -3.0 * atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case W:
            xi[0] = -3.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, point_on_circle, s[0], f);
            break;

          case E:
            xi[0] = -1.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, point_on_circle, s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      default:

        std::ostringstream error_stream;
        error_stream << "Wrong macro element number" << m << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  } // End of macro_element_boundary


  /// \short Parametrisation of macro element boundaries: f(s) is the position
  /// vector to macro-element m's boundary in the specified direction [N/S/E/W]
  /// at the specified discrete time level (time=0: present; time>0: previous)
  void RectangleWithHoleAndAnnularRegionDomain::macro_element_boundary(
    const unsigned& time,
    const unsigned& m,
    const unsigned& direction,
    const Vector<double>& s,
    Vector<double>& f)
  {
#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "RectangleWithHoleAndAnnularRegionDomain::macro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    // Lagrangian coordinate along surface of cylinder
    Vector<double> xi(1, 0.0);

    // Point on circle
    Vector<double> point_on_circle(2, 0.0);

    // Point on the outer boundary of the annular region
    Vector<double> point_on_annular_ring(2, 0.0);

    // Use the QuadTree enumeration for face directions
    using namespace QuadTreeNames;

    // Switch on the macro element
    switch (m)
    {
      // Macro element 0 is immediately to the left of the cylinder, outside
      // the inner annular region
      case 0:

        // Switch on the direction
        switch (direction)
        {
          case N:
            xi[0] = 3.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(Upper_mid_left, point_on_annular_ring, s[0], f);
            break;

          case S:
            xi[0] = -3.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(Lower_mid_left, point_on_annular_ring, s[0], f);
            break;

          case W:
            linear_interpolate(Lower_mid_left, Upper_mid_left, s[0], f);
            break;

          case E:
            xi[0] = 5.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 1 is immediately above the cylinder, outside
      // the inner annular region
      case 1:

        switch (direction)
        {
          case N:
            linear_interpolate(Upper_mid_left, Upper_mid_right, s[0], f);
            break;

          case S:
            xi[0] = 3.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case W:
            xi[0] = 3.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, Upper_mid_left, s[0], f);
            break;

          case E:
            xi[0] = 1.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, Upper_mid_right, s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect:  " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 2 is immediately to the right of the cylinder, outside
      // the inner annular region
      case 2:

        switch (direction)
        {
          case N:
            xi[0] = 1.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, Upper_mid_right, s[0], f);
            break;

          case S:
            xi[0] = -1.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, Lower_mid_right, s[0], f);
            break;

          case W:
            xi[0] = -atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case E:
            linear_interpolate(Lower_mid_right, Upper_mid_right, s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect:  " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 3 is immediately below cylinder, outside
      // the inner annular region
      case 3:

        switch (direction)
        {
          case N:
            xi[0] = -3.0 * atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case S:
            linear_interpolate(Lower_mid_left, Lower_mid_right, s[0], f);
            break;

          case W:
            xi[0] = -3.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(Lower_mid_left, point_on_annular_ring, s[0], f);
            break;

          case E:
            xi[0] = -1.0 * atan(1.0);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(Lower_mid_right, point_on_annular_ring, s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 4, is immediately to the left of the cylinder
      // lying in the inner annular region
      case 4:

        // Switch on the face of the m-th macro element
        switch (direction)
        {
          case N:
            xi[0] = 3.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, point_on_circle, s[0], f);
            break;

          case S:
            xi[0] = -3.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, point_on_circle, s[0], f);
            break;

          case W:
            xi[0] = 5.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case E:
            xi[0] = 5.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 5, is immediately above the cylinder in the
      // inner annular region
      case 5:

        switch (direction)
        {
          case N:
            xi[0] = 3.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case S:
            xi[0] = 3.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case W:
            xi[0] = 3.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_circle, point_on_annular_ring, s[0], f);
            break;

          case E:
            xi[0] = 1.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_circle, point_on_annular_ring, s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 6, is immediately to the right of the cylinder
      // lying in the inner annular region
      case 6:

        switch (direction)
        {
          case N:
            xi[0] = 1.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_circle, point_on_annular_ring, s[0], f);
            break;

          case S:
            xi[0] = -1.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_circle, point_on_annular_ring, s[0], f);
            break;

          case W:
            xi[0] = -atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case E:
            xi[0] = -atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      // Macro element 7, is immediately below the cylinder in the
      // inner annular region
      case 7:

        switch (direction)
        {
          case N:
            xi[0] = -3.0 * atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case S:
            xi[0] = -3.0 * atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            project_point_on_cylinder_to_annular_boundary(time, xi, f);
            break;

          case W:
            xi[0] = -3.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, point_on_circle, s[0], f);
            break;

          case E:
            xi[0] = -1.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            project_point_on_cylinder_to_annular_boundary(
              time, xi, point_on_annular_ring);
            linear_interpolate(point_on_annular_ring, point_on_circle, s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "Direction is incorrect: " << direction
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      default:

        std::ostringstream error_stream;
        error_stream << "Wrong macro element number" << m << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  } // End of macro_element_boundary

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=============================================================================
  /// Domain-based mesh for rectangular mesh with circular hole
  //=============================================================================
  /// Constructor: Pass pointer to geometric object that
  /// represents the cylinder, the length and height of the domain.
  /// The GeomObject must be parametrised such that
  /// \f$\zeta \in [0,2\pi]\f$ sweeps around the circumference
  /// in anticlockwise direction. Timestepper defaults to Steady
  /// default timestepper.
  template<class ELEMENT>
  RectangleWithHoleAndAnnularRegionMesh<ELEMENT>::
    RectangleWithHoleAndAnnularRegionMesh(GeomObject* cylinder_pt,
                                          const double& annular_region_radius,
                                          const double& length,
                                          TimeStepper* time_stepper_pt)
  {
    // Create the domain
    Domain_pt = new RectangleWithHoleAndAnnularRegionDomain(
      cylinder_pt, annular_region_radius, length);

    // Initialise the node counter
    unsigned long node_count = 0;

    // Vectors used to get data from domains
    Vector<double> s(2), r(2);

    // Setup temporary storage for the Nodes
    Vector<Node*> tmp_node_pt;

    // Number of macro elements in the domain
    unsigned nmacro_element = Domain_pt->nmacro_element();

    // Now blindly loop over the macro elements and associate a finite
    // element with each
    for (unsigned e = 0; e < nmacro_element; e++)
    {
      // Create the FiniteElement and add to the Element_pt Vector
      Element_pt.push_back(new ELEMENT);

      // Read out the number of linear points in the element
      unsigned np = dynamic_cast<ELEMENT*>(finite_element_pt(e))->nnode_1d();

      // Loop over nodes in the column
      for (unsigned l1 = 0; l1 < np; l1++)
      {
        // Loop over the nodes in the row
        for (unsigned l2 = 0; l2 < np; l2++)
        {
          // Allocate the memory for the node
          tmp_node_pt.push_back(finite_element_pt(e)->construct_node(
            l1 * np + l2, time_stepper_pt));

          // Read out the position of the node from the macro element
          s[0] = -1.0 + 2.0 * double(l2) / double(np - 1);
          s[1] = -1.0 + 2.0 * double(l1) / double(np - 1);

          // Use the macro element map from local coordinates to global
          // coordinates
          Domain_pt->macro_element_pt(e)->macro_map(s, r);

          // Set the position of the node
          tmp_node_pt[node_count]->x(0) = r[0];
          tmp_node_pt[node_count]->x(1) = r[1];

          // Increment the node number
          node_count++;
        }
      } // for (unsigned l1=0;l1<np;l1++)
    } // for (unsigned e=0;e<nmacro_element;e++)

    //----------------------------------------------------------------
    // Now the elements have been created, but there will be nodes in
    // common, need to loop over the common edges and sort it, by
    // reassigning pointers and the deleting excess nodes.
    //----------------------------------------------------------------
    // Read out the number of linear points in the element
    unsigned np = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

    // Edge between Elements 0 and 1 (miss out the node which coincides
    // with 4 elements; this will be dealt with separately later)
    for (unsigned n = 1; n < np; n++)
    {
      // Set the nodes in element 1 to be the same as in element 0
      finite_element_pt(1)->node_pt(n * np) =
        finite_element_pt(0)->node_pt((np * np - 1) - n);

      // Remove the nodes in element 1 from the temporary node list
      delete tmp_node_pt[np * np + n * np];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[np * np + n * np] = 0;
    }

    // Edge between Elements 0 and 3 (miss out the node which coincides
    // with 4 elements; this will be dealt with separately later)
    for (unsigned n = 0; n < np - 1; n++)
    {
      // Set the nodes in element 3 to be the same as in element 0
      finite_element_pt(3)->node_pt(n * np) = finite_element_pt(0)->node_pt(n);

      // Remove the nodes in element 3 from the temporary node list
      delete tmp_node_pt[3 * np * np + n * np];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[3 * np * np + n * np] = 0;
    }

    // Edge between Element 1 and 2 (miss out the node which coincides
    // with 4 elements; this will be dealt with separately later)
    for (unsigned n = 1; n < np; n++)
    {
      // Set the nodes in element 2 to be the same as in element 1
      finite_element_pt(2)->node_pt(np * (np - 1) + n) =
        finite_element_pt(1)->node_pt((np - 1) + n * np);

      // Remove the nodes in element 2 from the temporary node list
      delete tmp_node_pt[2 * np * np + np * (np - 1) + n];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[2 * np * np + np * (np - 1) + n] = 0;
    }

    // Edge between Element 3 and 2 (miss out the node which coincides
    // with 4 elements; this will be dealt with separately later)
    for (unsigned n = 1; n < np; n++)
    {
      // Set the nodes in element 2 to be the same as in element 3
      finite_element_pt(2)->node_pt(n) =
        finite_element_pt(3)->node_pt((np * np - 1) - n * np);

      // Remove the nodes in element 2 from the temporary node list
      delete tmp_node_pt[2 * np * np + n];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[2 * np * np + n] = 0;
    }

    // Edge between Elements 4 and 5 (miss out the node which coincides
    // with 4 elements; this will be dealt with separately later)
    for (unsigned n = 0; n < np - 1; n++)
    {
      // Set the nodes in element 1 to be the same as in element 0
      finite_element_pt(5)->node_pt(n * np) =
        finite_element_pt(4)->node_pt((np * np - 1) - n);

      // Remove the nodes in element 1 from the temporary node list
      delete tmp_node_pt[5 * np * np + n * np];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[5 * np * np + n * np] = 0;
    }

    // Edge between Elements 4 and 7 (miss out the node which coincides
    // with 4 elements; this will be dealt with separately later)
    for (unsigned n = 1; n < np; n++)
    {
      // Set the nodes in element 3 to be the same as in element 0
      finite_element_pt(7)->node_pt(n * np) = finite_element_pt(4)->node_pt(n);

      // Remove the nodes in element 3 from the temporary node list
      delete tmp_node_pt[7 * np * np + n * np];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[7 * np * np + n * np] = 0;
    }

    // Edge between Element 5 and 6 (miss out the node which coincides
    // with 4 elements; this will be dealt with separately later)
    for (unsigned n = 0; n < np - 1; n++)
    {
      // Set the nodes in element 2 to be the same as in element 1
      finite_element_pt(6)->node_pt(np * (np - 1) + n) =
        finite_element_pt(5)->node_pt((n + 1) * np - 1);

      // Remove the nodes in element 2 from the temporary node list
      delete tmp_node_pt[6 * np * np + np * (np - 1) + n];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[6 * np * np + np * (np - 1) + n] = 0;
    }

    // Edge between Element 7 and 6 (miss out the node which coincides
    // with 4 elements; this will be dealt with separately later)
    for (unsigned n = 0; n < np - 1; n++)
    {
      // Set the nodes in element 2 to be the same as in element 3
      finite_element_pt(6)->node_pt(n) =
        finite_element_pt(7)->node_pt((np * np - 1) - n * np);

      // Remove the nodes in element 2 from the temporary node list
      delete tmp_node_pt[6 * np * np + n];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[6 * np * np + n] = 0;
    }

    // Edge between Elements 0 and 4 (miss out the corner nodes as they have
    // coincide with 4 elements rather than just 2)
    for (unsigned n = 1; n < np - 1; n++)
    {
      // Set the nodes in element 1 to be the same as in element 0
      finite_element_pt(0)->node_pt((np - 1) + n * np) =
        finite_element_pt(4)->node_pt(n * np);

      // Remove the nodes in element 1 from the temporary node list
      delete tmp_node_pt[(n + 1) * np - 1];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[(n + 1) * np - 1] = 0;
    }

    // Edge between Elements 1 and 5 (miss out the corner nodes as they have
    // coincide with 4 elements rather than just 2)
    for (unsigned n = 1; n < np - 1; n++)
    {
      // Set the nodes in element 3 to be the same as in element 0
      finite_element_pt(1)->node_pt(n) =
        finite_element_pt(5)->node_pt(np * (np - 1) + n);

      // Remove the nodes in element 1 from the temporary node list
      delete tmp_node_pt[np * np + n];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[np * np + n] = 0;
    }

    // Edge between Element 2 and 6 (miss out the corner nodes as they have
    // coincide with 4 elements rather than just 2)
    for (unsigned n = 1; n < np - 1; n++)
    {
      // Set the nodes in element 2 to be the same as in element 1
      finite_element_pt(2)->node_pt(np * (np - 1) - n * np) =
        finite_element_pt(6)->node_pt((np * np - 1) - n * np);

      // Remove the nodes in element 1 from the temporary node list
      delete tmp_node_pt[2 * np * np + np * (np - 1) - n * np];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[2 * np * np + np * (np - 1) - n * np] = 0;
    }

    // Edge between Element 3 and 7 (miss out the corner nodes as they have
    // coincide with 4 elements rather than just 2)
    for (unsigned n = 1; n < np - 1; n++)
    {
      // Set the nodes in element 2 to be the same as in element 3
      finite_element_pt(3)->node_pt((np * np - 1) - n) =
        finite_element_pt(7)->node_pt((np - 1) - n);

      // Remove the nodes in element 1 from the temporary node list
      delete tmp_node_pt[3 * np * np + (np * np - 1) - n];

      // Make the corresponding pointer a null pointer
      tmp_node_pt[3 * np * np + (np * np - 1) - n] = 0;
    }

    //--------------------------------------------------------------------
    // Now we'll deal with the corner nodes which lie in 4 elements rather
    // than just two elements. There are four of these to deal with.
    // Specifically:
    //     Node 0: Lies in elements 0, 3, 4 and 7 (copy from element 0)
    //     Node 1: Lies in elements 0, 1, 4 and 5 (copy from element 0)
    //     Node 2: Lies in elements 1, 2, 5 and 6 (copy from element 1)
    //     Node 3: Lies in elements 2, 3, 6 and 7 (copy from element 2)
    //--------------------------------------------------------------------
    //---------------
    // Corner node 0:
    //---------------
    // Set the corner node in element 3 to be the same as in element 0
    finite_element_pt(3)->node_pt(np * (np - 1)) =
      finite_element_pt(0)->node_pt(np - 1);

    // Set the corner node in element 4 to be the same as in element 0
    finite_element_pt(4)->node_pt(0) = finite_element_pt(0)->node_pt(np - 1);

    // Set the corner node in element 7 to be the same as in element 0
    finite_element_pt(7)->node_pt(0) = finite_element_pt(0)->node_pt(np - 1);

    // Remove the overwritten node in element 3 from the temporary node list
    delete tmp_node_pt[3 * np * np + np * (np - 1)];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[3 * np * np + np * (np - 1)] = 0;

    // Remove the overwritten node in element 4 from the temporary node list
    delete tmp_node_pt[4 * np * np];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[4 * np * np] = 0;

    // Remove the overwritten node in element 7 from the temporary node list
    delete tmp_node_pt[7 * np * np];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[7 * np * np] = 0;

    //---------------
    // Corner node 1:
    //---------------
    // Set the corner node in element 3 to be the same as in element 0
    finite_element_pt(1)->node_pt(0) =
      finite_element_pt(0)->node_pt(np * np - 1);

    // Set the corner node in element 4 to be the same as in element 0
    finite_element_pt(4)->node_pt(np * (np - 1)) =
      finite_element_pt(0)->node_pt(np * np - 1);

    // Set the corner node in element 7 to be the same as in element 0
    finite_element_pt(5)->node_pt(np * (np - 1)) =
      finite_element_pt(0)->node_pt(np * np - 1);

    // Remove the overwritten node in element 1 from the temporary node list
    delete tmp_node_pt[np * np];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[np * np] = 0;

    // Remove the overwritten node in element 4 from the temporary node list
    delete tmp_node_pt[4 * np * np + np * (np - 1)];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[4 * np * np + np * (np - 1)] = 0;

    // Remove the overwritten node in element 5 from the temporary node list
    delete tmp_node_pt[5 * np * np + np * (np - 1)];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[5 * np * np + np * (np - 1)] = 0;

    //---------------
    // Corner node 2:
    //---------------
    // Set the corner node in element 2 to be the same as in element 1
    finite_element_pt(2)->node_pt(np * (np - 1)) =
      finite_element_pt(1)->node_pt(np - 1);

    // Set the corner node in element 5 to be the same as in element 1
    finite_element_pt(5)->node_pt(np * np - 1) =
      finite_element_pt(1)->node_pt(np - 1);

    // Set the corner node in element 6 to be the same as in element 1
    finite_element_pt(6)->node_pt(np * np - 1) =
      finite_element_pt(1)->node_pt(np - 1);

    // Remove the overwritten node in element 2 from the temporary node list
    delete tmp_node_pt[2 * np * np + np * (np - 1)];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[2 * np * np + np * (np - 1)] = 0;

    // Remove the overwritten node in element 5 from the temporary node list
    delete tmp_node_pt[5 * np * np + np * np - 1];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[5 * np * np + np * np - 1] = 0;

    // Remove the overwritten node in element 6 from the temporary node list
    delete tmp_node_pt[6 * np * np + np * np - 1];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[6 * np * np + np * np - 1] = 0;

    //---------------
    // Corner node 3:
    //---------------
    // Set the corner node in element 3 to be the same as in element 2
    finite_element_pt(3)->node_pt(np * np - 1) =
      finite_element_pt(2)->node_pt(0);

    // Set the corner node in element 4 to be the same as in element 2
    finite_element_pt(6)->node_pt(np - 1) = finite_element_pt(2)->node_pt(0);

    // Set the corner node in element 7 to be the same as in element 2
    finite_element_pt(7)->node_pt(np - 1) = finite_element_pt(2)->node_pt(0);

    // Remove the overwritten node in element 2 from the temporary node list
    delete tmp_node_pt[3 * np * np + np * np - 1];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[3 * np * np + np * np - 1] = 0;

    // Remove the overwritten node in element 5 from the temporary node list
    delete tmp_node_pt[6 * np * np + np - 1];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[6 * np * np + np - 1] = 0;

    // Remove the overwritten node in element 6 from the temporary node list
    delete tmp_node_pt[7 * np * np + np - 1];

    // Make the corresponding pointer a null pointer
    tmp_node_pt[7 * np * np + np - 1] = 0;

    //----------------------------------------------------------------------
    // Now all the nodes have been set up and the nodes which coincided with
    // each other have been deleted and the correct pointers have been set.
    // All that remains is to copy the nodes in the temporary storage into
    // the global storage vector Node_pt.
    //----------------------------------------------------------------------
    // Now set the actual true nodes
    for (unsigned long n = 0; n < node_count; n++)
    {
      // If it's not a null pointer
      if (tmp_node_pt[n] != 0)
      {
        // Add it to the global Node list
        Node_pt.push_back(tmp_node_pt[n]);
      }
    } // for (unsigned long n=0;n<node_count;n++)

    //-------------------------------------------------------------
    // Now sort out which boundary each node lies on. The outer box
    // has boundary numbering 0 for the South face, 1 for the East
    // face, 2 for the North face and 3 for the West face. Finally,
    // the cylinder boundary corresponds to boundary 4.
    //-------------------------------------------------------------
    // Finally set the nodes on the boundaries
    set_nboundary(5);

    // Loop over the nodes along an edge
    for (unsigned n = 0; n < np; n++)
    {
      // West face
      Node* nod_pt = finite_element_pt(0)->node_pt(n * np);
      convert_to_boundary_node(nod_pt);
      add_boundary_node(3, nod_pt);

      // North face
      nod_pt = finite_element_pt(1)->node_pt(np * (np - 1) + n);
      convert_to_boundary_node(nod_pt);
      add_boundary_node(2, nod_pt);

      // East face
      nod_pt = finite_element_pt(2)->node_pt(n * np + np - 1);
      convert_to_boundary_node(nod_pt);
      add_boundary_node(1, nod_pt);

      // South face
      nod_pt = finite_element_pt(3)->node_pt(n);
      convert_to_boundary_node(nod_pt);
      add_boundary_node(0, nod_pt);
    }

    // Loop over all but one node on an edge
    for (unsigned n = 0; n < np - 1; n++)
    {
      // First part of hole
      Node* nod_pt = finite_element_pt(4)->node_pt((np - 1) + n * np);
      convert_to_boundary_node(nod_pt);
      add_boundary_node(4, nod_pt);

      // Next part of hole
      nod_pt = finite_element_pt(5)->node_pt(n);
      convert_to_boundary_node(nod_pt);
      add_boundary_node(4, nod_pt);

      // Next part of hole
      nod_pt = finite_element_pt(6)->node_pt(np * (np - 1) - n * np);
      convert_to_boundary_node(nod_pt);
      add_boundary_node(4, nod_pt);

      // Final part of hole boundary
      nod_pt = finite_element_pt(7)->node_pt((np * np - 1) - n);
      convert_to_boundary_node(nod_pt);
      add_boundary_node(4, nod_pt);
    }

#ifdef PARANOID
    // Make sure there are no duplicate nodes (i.e. two or more nodes that
    // occupy the same Eulerian position)
    if (check_for_repeated_nodes())
    {
      // Throw a warning
      OomphLibWarning("The mesh contains repeated nodes!\n",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif
  } // End of RectangleWithHoleAndAnnularRegionMesh

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=============================================================================
  /// Constructor. Pass pointer to geometric object that
  /// represents the cylinder; hierher the length and height of the domain.
  /// The GeomObject must be parametrised such that
  /// \f$\zeta \in [0,2\pi]\f$ sweeps around the circumference
  /// in anticlockwise direction. Timestepper defaults to Steady
  /// default timestepper.
  //=============================================================================
  template<class ELEMENT>
  RefineableQuadMeshWithMovingCylinder<ELEMENT>::
    RefineableQuadMeshWithMovingCylinder(GeomObject* cylinder_pt,
                                         const double& annular_region_radius,
                                         const double& length_of_central_box,
                                         const double& x_left,
                                         const double& x_right,
                                         const double& height,
                                         TimeStepper* time_stepper_pt)
  {
    // Do you want a coarse problem?
    Coarse_problem = false; // true; // hierher

    // Vector for pointers to all the temporary meshes
    Vector<Mesh*> all_temp_mesh_pt;

    // Build central mesh
    Central_mesh_pt =
      new RefineableRectangleWithHoleAndAnnularRegionMesh<ELEMENT>(
        cylinder_pt,
        annular_region_radius,
        length_of_central_box,
        time_stepper_pt);

    // Add this mesh to the list of temproary meshes
    all_temp_mesh_pt.push_back(Central_mesh_pt);

    // Bulk mesh left boundary id
    unsigned int left_boundary_id = 3;

    // Bulk mesh top boundary id
    unsigned int top_boundary_id = 2;

    // Bulk mesh right boundary id
    unsigned int right_boundary_id = 1;

    // Bulk mesh bottom boundary id
    unsigned int bottom_boundary_id = 0;

    // Build the surrounding meshes

    // Right mesh
    double l_right = x_right - 0.5 * length_of_central_box;
    unsigned n_right =
      std::max(unsigned(1), unsigned(l_right / length_of_central_box));
    if (Coarse_problem)
    {
      n_right = 1;
    }

    Mesh* surrounding_right_mesh_pt =
      TwoDimensionalPMLHelper::create_right_pml_mesh<PMLLayerElement<ELEMENT>>(
        Central_mesh_pt, right_boundary_id, n_right, l_right, time_stepper_pt);

    all_temp_mesh_pt.push_back(surrounding_right_mesh_pt);

    // Remove nodes from top and bottom boundaries (0 and 2)
    for (unsigned ibound_rm = 0; ibound_rm <= 2; ibound_rm += 2)
    {
      unsigned num_nod = surrounding_right_mesh_pt->nboundary_node(ibound_rm);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        Node* nod_pt =
          surrounding_right_mesh_pt->boundary_node_pt(ibound_rm, inod);
        if (nod_pt->is_on_boundary(ibound_rm))
        {
          nod_pt->remove_from_boundary(ibound_rm);
        }
      }
    }

    // Top mesh
    double l_top = 0.5 * height - 0.5 * length_of_central_box;
    unsigned n_top =
      std::max(unsigned(1), unsigned(l_top / length_of_central_box));
    if (Coarse_problem)
    {
      n_top = 1;
    }
    Mesh* surrounding_top_mesh_pt =
      TwoDimensionalPMLHelper::create_top_pml_mesh<PMLLayerElement<ELEMENT>>(
        Central_mesh_pt, top_boundary_id, n_top, l_top, time_stepper_pt);
    all_temp_mesh_pt.push_back(surrounding_top_mesh_pt);

    // Remove nodes from right and left boundaries (1 and 3)
    for (unsigned ibound_rm = 1; ibound_rm <= 3; ibound_rm += 2)
    {
      unsigned num_nod = surrounding_top_mesh_pt->nboundary_node(ibound_rm);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        Node* nod_pt =
          surrounding_top_mesh_pt->boundary_node_pt(ibound_rm, inod);
        if (nod_pt->is_on_boundary(ibound_rm))
        {
          nod_pt->remove_from_boundary(ibound_rm);
        }
      }
    }


    // Left mesh
    double l_left = -0.5 * length_of_central_box - x_left;
    unsigned n_left =
      std::max(unsigned(1), unsigned(l_left / length_of_central_box));
    if (Coarse_problem)
    {
      n_left = 1;
    }
    Mesh* surrounding_left_mesh_pt =
      TwoDimensionalPMLHelper::create_left_pml_mesh<PMLLayerElement<ELEMENT>>(
        Central_mesh_pt, left_boundary_id, n_left, l_left, time_stepper_pt);
    all_temp_mesh_pt.push_back(surrounding_left_mesh_pt);

    // Remove nodes from top and bottom boundaries (0 and 2)
    for (unsigned ibound_rm = 0; ibound_rm <= 2; ibound_rm += 2)
    {
      unsigned num_nod = surrounding_left_mesh_pt->nboundary_node(ibound_rm);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        Node* nod_pt =
          surrounding_left_mesh_pt->boundary_node_pt(ibound_rm, inod);
        if (nod_pt->is_on_boundary(ibound_rm))
        {
          nod_pt->remove_from_boundary(ibound_rm);
        }
      }
    }

    // Bottom mesh
    double l_bottom = 0.5 * height - 0.5 * length_of_central_box;
    unsigned n_bottom =
      std::max(unsigned(1), unsigned(l_bottom / length_of_central_box));
    if (Coarse_problem)
    {
      n_bottom = 1;
    }
    Mesh* surrounding_bottom_mesh_pt =
      TwoDimensionalPMLHelper::create_bottom_pml_mesh<PMLLayerElement<ELEMENT>>(
        Central_mesh_pt,
        bottom_boundary_id,
        n_bottom,
        l_bottom,
        time_stepper_pt);
    all_temp_mesh_pt.push_back(surrounding_bottom_mesh_pt);


    // Remove nodes from right and left boundaries (1 and 3)
    for (unsigned ibound_rm = 1; ibound_rm <= 3; ibound_rm += 2)
    {
      unsigned num_nod = surrounding_bottom_mesh_pt->nboundary_node(ibound_rm);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        Node* nod_pt =
          surrounding_bottom_mesh_pt->boundary_node_pt(ibound_rm, inod);
        if (nod_pt->is_on_boundary(ibound_rm))
        {
          nod_pt->remove_from_boundary(ibound_rm);
        }
      }
    }

    // Build corner Surrounding meshes
    Mesh* surrounding_top_right_mesh_pt =
      TwoDimensionalPMLHelper::create_top_right_pml_mesh<
        PMLLayerElement<ELEMENT>>(surrounding_right_mesh_pt,
                                  surrounding_top_mesh_pt,
                                  Central_mesh_pt,
                                  right_boundary_id,
                                  time_stepper_pt);
    all_temp_mesh_pt.push_back(surrounding_top_right_mesh_pt);

    Mesh* surrounding_bottom_right_mesh_pt =
      TwoDimensionalPMLHelper::create_bottom_right_pml_mesh<
        PMLLayerElement<ELEMENT>>(surrounding_right_mesh_pt,
                                  surrounding_bottom_mesh_pt,
                                  Central_mesh_pt,
                                  right_boundary_id,
                                  time_stepper_pt);
    all_temp_mesh_pt.push_back(surrounding_bottom_right_mesh_pt);

    Mesh* surrounding_top_left_mesh_pt =
      TwoDimensionalPMLHelper::create_top_left_pml_mesh<
        PMLLayerElement<ELEMENT>>(surrounding_left_mesh_pt,
                                  surrounding_top_mesh_pt,
                                  Central_mesh_pt,
                                  left_boundary_id,
                                  time_stepper_pt);
    all_temp_mesh_pt.push_back(surrounding_top_left_mesh_pt);

    Mesh* surrounding_bottom_left_mesh_pt =
      TwoDimensionalPMLHelper::create_bottom_left_pml_mesh<
        PMLLayerElement<ELEMENT>>(surrounding_left_mesh_pt,
                                  surrounding_bottom_mesh_pt,
                                  Central_mesh_pt,
                                  left_boundary_id,
                                  time_stepper_pt);
    all_temp_mesh_pt.push_back(surrounding_bottom_left_mesh_pt);


    // Copy all elements across
    unsigned n_mesh = all_temp_mesh_pt.size();
    unsigned nel = 0;
    unsigned nnod = 0;
    for (unsigned i = 0; i < n_mesh; i++)
    {
      nel += all_temp_mesh_pt[i]->nelement();
      nnod += all_temp_mesh_pt[i]->nnode();
    }
    this->Element_pt.resize(nel);
    this->Node_pt.resize(nnod);

    unsigned count_el = 0;
    unsigned count_nod = 0;
    for (unsigned i = 0; i < n_mesh; i++)
    {
      unsigned nel = all_temp_mesh_pt[i]->nelement();
      for (unsigned e = 0; e < nel; e++)
      {
        this->Element_pt[count_el] = all_temp_mesh_pt[i]->element_pt(e);
        count_el++;
      }
      unsigned nnod = all_temp_mesh_pt[i]->nnode();
      for (unsigned j = 0; j < nnod; j++)
      {
        this->Node_pt[count_nod] = all_temp_mesh_pt[i]->node_pt(j);
        count_nod++;
      }
    }

    // Remove nodes on outer boundaries (0 to 3) of central mesh
    // from their respective boundaries (otherwise they get added
    // to the mesh boundaries during adaptation)
    for (unsigned b = 0; b < 4; b++)
    {
      unsigned nnod = Central_mesh_pt->nboundary_node(b);
      for (unsigned j = 0; j < nnod; j++)
      {
        Node* nod_pt = Central_mesh_pt->boundary_node_pt(b, j);
        nod_pt->remove_from_boundary(b);
      }
    }

    // Setup boundaries:
    //==================
    this->set_nboundary(5);


    Vector<Mesh*> pml_mesh_pt;

    // Top meshes:
    //------------
    {
      // Upper boundary in top mesh (same number in this mesh)
      unsigned ibound = 2;

      // Loop over relevant meshes
      Vector<Mesh*> all_mesh_pt;
      all_mesh_pt.push_back(surrounding_top_left_mesh_pt);
      all_mesh_pt.push_back(surrounding_top_mesh_pt);
      all_mesh_pt.push_back(surrounding_top_right_mesh_pt);
      unsigned n_mesh = all_mesh_pt.size();
      for (unsigned m = 0; m < n_mesh; m++)
      {
        pml_mesh_pt.push_back(all_mesh_pt[m]);
        unsigned num_nod = all_mesh_pt[m]->nboundary_node(ibound);
        for (unsigned inod = 0; inod < num_nod; inod++)
        {
          Node* nod_pt = all_mesh_pt[m]->boundary_node_pt(ibound, inod);
          this->add_boundary_node(ibound, nod_pt);
        }
      }
    }

    // Left meshes:
    //------------
    {
      // Left boundary (same number in this mesh)
      unsigned ibound = 3;

      pml_mesh_pt.push_back(surrounding_left_mesh_pt);

      // Loop over relevant meshes
      Vector<Mesh*> all_mesh_pt;
      all_mesh_pt.push_back(surrounding_top_left_mesh_pt);
      all_mesh_pt.push_back(surrounding_left_mesh_pt);
      all_mesh_pt.push_back(surrounding_bottom_left_mesh_pt);
      unsigned n_mesh = all_mesh_pt.size();
      for (unsigned m = 0; m < n_mesh; m++)
      {
        unsigned num_nod = all_mesh_pt[m]->nboundary_node(ibound);
        for (unsigned inod = 0; inod < num_nod; inod++)
        {
          Node* nod_pt = all_mesh_pt[m]->boundary_node_pt(ibound, inod);
          this->add_boundary_node(ibound, nod_pt);
        }
      }
    }

    // Bottom meshes:
    //---------------
    {
      // Bottom boundary (same number in this mesh)
      unsigned ibound = 0;

      // Loop over relevant meshes
      Vector<Mesh*> all_mesh_pt;
      all_mesh_pt.push_back(surrounding_bottom_left_mesh_pt);
      all_mesh_pt.push_back(surrounding_bottom_mesh_pt);
      all_mesh_pt.push_back(surrounding_bottom_right_mesh_pt);
      unsigned n_mesh = all_mesh_pt.size();
      for (unsigned m = 0; m < n_mesh; m++)
      {
        pml_mesh_pt.push_back(all_mesh_pt[m]);
        unsigned num_nod = all_mesh_pt[m]->nboundary_node(ibound);
        for (unsigned inod = 0; inod < num_nod; inod++)
        {
          Node* nod_pt = all_mesh_pt[m]->boundary_node_pt(ibound, inod);
          this->add_boundary_node(ibound, nod_pt);
        }
      }
    }

    // Right meshes:
    //--------------
    {
      // Right boundary  (same number in this mesh)
      unsigned ibound = 1;

      pml_mesh_pt.push_back(surrounding_right_mesh_pt);

      // Loop over relevant meshes
      Vector<Mesh*> all_mesh_pt;
      all_mesh_pt.push_back(surrounding_bottom_right_mesh_pt);
      all_mesh_pt.push_back(surrounding_right_mesh_pt);
      all_mesh_pt.push_back(surrounding_top_right_mesh_pt);
      unsigned n_mesh = all_mesh_pt.size();
      for (unsigned m = 0; m < n_mesh; m++)
      {
        unsigned num_nod = all_mesh_pt[m]->nboundary_node(ibound);
        for (unsigned inod = 0; inod < num_nod; inod++)
        {
          Node* nod_pt = all_mesh_pt[m]->boundary_node_pt(ibound, inod);
          this->add_boundary_node(ibound, nod_pt);
        }
      }
    }

    // No slip on boundary of inner mesh
    //----------------------------------
    {
      // Inner boundary
      unsigned ibound = 4;

      // Pin both velocities
      unsigned num_nod = Central_mesh_pt->nboundary_node(ibound);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        Node* nod_pt = Central_mesh_pt->boundary_node_pt(ibound, inod);
        this->add_boundary_node(ibound, nod_pt);
      }
    }

    // Setup quadtree forest for mesh refinement
    this->setup_quadtree_forest();

    // Setup boundary element lookup schemes
    this->setup_boundary_element_info();

    // Cleanup. NOTE: Can't delete Central_mesh_pt as it's responsible for
    // deleting Domain_pt but
    unsigned n_pml_mesh = pml_mesh_pt.size();
    for (unsigned j = 0; j < n_pml_mesh; j++)
    {
      pml_mesh_pt[j]->flush_element_and_node_storage();
      delete pml_mesh_pt[j];
      pml_mesh_pt[j] = 0;
    }
    Central_mesh_pt->flush_element_and_node_storage();
  } // End of RefineableQuadMeshWithMovingCylinder
} // namespace oomph
