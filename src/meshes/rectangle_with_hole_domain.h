// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_RECTANGLE_WITH_HOLE_DOMAIN_HEADER
#define OOMPH_RECTANGLE_WITH_HOLE_DOMAIN_HEADER


// Generic includes
#include "../generic/quadtree.h"
#include "../generic/geom_objects.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"


namespace oomph
{
  //===========================================================
  /// Rectangular domain with circular whole
  //===========================================================
  class RectangleWithHoleDomain : public Domain
  {
  public:
    /// Constructor. Pass pointer to geometric object that
    /// represents the cylinder, the length of the (square) domain.
    /// The GeomObject must be parametrised such that
    /// \f$\zeta \in [0,2\pi]\f$ sweeps around the circumference
    /// in anticlockwise direction.
    RectangleWithHoleDomain(GeomObject* cylinder_pt, const double& length)
      : Cylinder_pt(cylinder_pt)
    {
      // Vertices of rectangle
      Lower_left.resize(2);
      Lower_left[0] = -0.5 * length;
      Lower_left[1] = -0.5 * length;

      Upper_left.resize(2);
      Upper_left[0] = -0.5 * length;
      Upper_left[1] = 0.5 * length;

      Lower_right.resize(2);
      Lower_right[0] = 0.5 * length;
      Lower_right[1] = -0.5 * length;

      Upper_right.resize(2);
      Upper_right[0] = 0.5 * length;
      Upper_right[1] = 0.5 * length;


      // Coordinates of points where the "radial" lines from central
      // cylinder meet the upper and lower boundaries
      Lower_mid_left.resize(2);
      Lower_mid_left[0] = -0.5 * length;
      Lower_mid_left[1] = -0.5 * length;

      Upper_mid_left.resize(2);
      Upper_mid_left[0] = -0.5 * length;
      Upper_mid_left[1] = 0.5 * length;

      Lower_mid_right.resize(2);
      Lower_mid_right[0] = 0.5 * length;
      Lower_mid_right[1] = -0.5 * length;

      Upper_mid_right.resize(2);
      Upper_mid_right[0] = 0.5 * length;
      Upper_mid_right[1] = 0.5 * length;


      // There are four macro elements
      Macro_element_pt.resize(4);

      // Build the 2D macro elements
      for (unsigned i = 0; i < 4; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    }


    /// Destructor: Empty; cleanup done in base class
    ~RectangleWithHoleDomain() {}

    /// Helper function to interpolate linearly between the
    /// "right" and "left" points; \f$ s \in [-1,1] \f$
    void linear_interpolate(Vector<double> left,
                            Vector<double> right,
                            const double& s,
                            Vector<double>& f)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        f[i] = left[i] + (right[i] - left[i]) * 0.5 * (s + 1.0);
      }
    }


    /// Parametrisation of macro element boundaries: f(s) is the position
    /// vector to macro-element m's boundary in the specified direction
    /// [N/S/E/W] at the specfied discrete time level (time=0: present; time>0:
    /// previous)
    void macro_element_boundary(const unsigned& time,
                                const unsigned& m,
                                const unsigned& direction,
                                const Vector<double>& s,
                                Vector<double>& f)
    {
#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
      // Warn about time argument being moved to the front
      OomphLibWarning(
        "Order of function arguments has changed between versions 0.8 and 0.85",
        "RectangleWithHoleDomain::macro_element_boundary(...)",
        OOMPH_EXCEPTION_LOCATION);
#endif

      // Lagrangian coordinate along surface of cylinder
      Vector<double> xi(1);

      // Point on circle
      Vector<double> point_on_circle(2);

      using namespace QuadTreeNames;

      // Switch on the macro element
      switch (m)
      {
          // Macro element 0, is is immediately left of the cylinder
        case 0:

          switch (direction)
          {
            case N:
              xi[0] = 3.0 * atan(1.0);
              Cylinder_pt->position(time, xi, point_on_circle);
              linear_interpolate(Upper_mid_left, point_on_circle, s[0], f);
              break;

            case S:
              xi[0] = -3.0 * atan(1.0);
              Cylinder_pt->position(time, xi, point_on_circle);
              linear_interpolate(Lower_mid_left, point_on_circle, s[0], f);
              break;

            case W:
              linear_interpolate(Lower_mid_left, Upper_mid_left, s[0], f);
              break;

            case E:
              xi[0] = 5.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
              Cylinder_pt->position(time, xi, f);
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

        // Macro element 1, is immediately above the cylinder
        case 1:

          switch (direction)
          {
            case N:
              linear_interpolate(Upper_mid_left, Upper_mid_right, s[0], f);
              break;

            case S:
              xi[0] = 3.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
              Cylinder_pt->position(time, xi, f);
              break;

            case W:
              xi[0] = 3.0 * atan(1.0);
              Cylinder_pt->position(time, xi, point_on_circle);
              linear_interpolate(point_on_circle, Upper_mid_left, s[0], f);
              break;

            case E:
              xi[0] = 1.0 * atan(1.0);
              Cylinder_pt->position(time, xi, point_on_circle);
              linear_interpolate(point_on_circle, Upper_mid_right, s[0], f);
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

          // Macro element 2, is immediately right of the cylinder
        case 2:

          switch (direction)
          {
            case N:
              xi[0] = 1.0 * atan(1.0);
              Cylinder_pt->position(time, xi, point_on_circle);
              linear_interpolate(point_on_circle, Upper_mid_right, s[0], f);
              break;

            case S:
              xi[0] = -1.0 * atan(1.0);
              Cylinder_pt->position(time, xi, point_on_circle);
              linear_interpolate(point_on_circle, Lower_mid_right, s[0], f);
              break;

            case W:
              xi[0] = -atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
              Cylinder_pt->position(time, xi, f);
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

          // Macro element 3, is immediately below cylinder
        case 3:

          switch (direction)
          {
            case N:
              xi[0] = -3.0 * atan(1.0) + 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
              Cylinder_pt->position(time, xi, f);
              break;

            case S:
              linear_interpolate(Lower_mid_left, Lower_mid_right, s[0], f);
              break;

            case W:
              xi[0] = -3.0 * atan(1.0);
              Cylinder_pt->position(time, xi, point_on_circle);
              linear_interpolate(Lower_mid_left, point_on_circle, s[0], f);
              break;

            case E:
              xi[0] = -1.0 * atan(1.0);
              Cylinder_pt->position(time, xi, point_on_circle);
              linear_interpolate(Lower_mid_right, point_on_circle, s[0], f);
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

        default:

          std::ostringstream error_stream;
          error_stream << "Wrong macro element number" << m << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }


  private:
    /// Lower left corner of rectangle
    Vector<double> Lower_left;

    /// Lower right corner of rectangle
    Vector<double> Lower_right;

    /// Where the "radial" line from circle meets lower boundary on left
    Vector<double> Lower_mid_left;

    /// Where the "radial" line from circle meets lower boundary on right
    Vector<double> Lower_mid_right;

    /// Upper left corner of rectangle
    Vector<double> Upper_left;

    /// Upper right corner of rectangle
    Vector<double> Upper_right;

    /// Where the "radial" line from circle meets upper boundary on left
    Vector<double> Upper_mid_left;

    /// Where the "radial" line from circle meets upper boundary on right
    Vector<double> Upper_mid_right;

    /// Pointer to geometric object that represents the central cylinder
    GeomObject* Cylinder_pt;
  };


} // namespace oomph
#endif
