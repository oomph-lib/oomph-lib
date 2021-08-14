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
#include "cylinder_with_flag_domain.h"


namespace oomph
{
  //=======================================================================
  /// Constructor, Pass the pointers to the GeomObjects that parametrise
  /// the cylinder, the three edges of the flag, the length and height of the
  ///  domain, the length and height of the flag, the coordinates of the
  /// centre of the cylinder and its radius.
  //=======================================================================
  CylinderWithFlagDomain::CylinderWithFlagDomain(Circle* cylinder_pt,
                                                 GeomObject* top_flag_pt,
                                                 GeomObject* bottom_flag_pt,
                                                 GeomObject* tip_flag_pt,
                                                 const double& length,
                                                 const double& height,
                                                 const double& flag_length,
                                                 const double& flag_height,
                                                 const double& centre_x,
                                                 const double& centre_y,
                                                 const double& a)
    : Cylinder_pt(cylinder_pt),
      Top_flag_pt(top_flag_pt),
      Bottom_flag_pt(bottom_flag_pt),
      Tip_flag_pt(tip_flag_pt),
      Lx(flag_length),
      Ly(flag_height),
      Centre_x(centre_x),
      Centre_y(centre_y),
      A(a)
  {
    // Vertices of rectangle
    // Those are points of references of the domain
    // to help create the macro_element_boundary sub_functions
    p1.resize(2);
    p1[0] = 0.0;
    p1[1] = height;

    p2.resize(2);
    p2[0] = Centre_x;
    p2[1] = height;

    p3.resize(2);
    p3[0] = 0.155596 * length;
    p3[1] = height;

    p4.resize(2);
    p4[0] = 0.183596 * length;
    p4[1] = height;

    p5.resize(2);
    p5[0] = 0.239596 * length;
    p5[1] = height;

    p6.resize(2);
    p6[0] = 0.285123967 * length;
    p6[1] = height;

    p7.resize(2);
    p7[0] = 0.433884298 * length;
    p7[1] = height;

    p8.resize(2);
    p8[0] = 0.578512397 * length;
    p8[1] = height;

    p9.resize(2);
    p9[0] = 0.789256198 * length;
    p9[1] = height;

    p10.resize(2);
    p10[0] = length;
    p10[1] = height;

    p11.resize(2);
    p11[0] = 0.127596 * length;
    p11[1] = 0.778024390 * height;

    p12.resize(2);
    p12[0] = 0.155596 * length;
    p12[1] = 0.778024390 * height;

    p13.resize(2);
    p13[0] = 0.183596 * length;
    p13[1] = 0.778024390 * height;

    p14.resize(2);
    p14[0] = 0.211596 * length;
    p14[1] = 0.778024390 * height;

    p15.resize(2);
    p15[0] = 0.285123967 * length;
    p15[1] = 0.625 * height;

    p16.resize(2);
    p16[0] = 0.351239669 * length;
    p16[1] = 0.625 * height;

    p18.resize(2);
    p18[0] = Centre_x;
    p18[1] = Centre_y + A;

    p33.resize(2);
    p33[0] = Centre_x;
    p33[1] = Centre_y - A;

    p35.resize(2);
    p35[0] = 0.285123967 * length;
    p35[1] = 0.350609756 * height;

    p36.resize(2);
    p36[0] = 0.351239669 * length;
    p36[1] = 0.350609756 * height;

    p37.resize(2);
    p37[0] = 0.127596 * length;
    p37[1] = 0.197585366 * height;

    p38.resize(2);
    p38[0] = 0.155596 * length;
    p38[1] = 0.197585366 * height;

    p39.resize(2);
    p39[0] = 0.183596 * length;
    p39[1] = 0.197585366 * height;

    p40.resize(2);
    p40[0] = 0.211596 * length;
    p40[1] = 0.197585366 * height;

    p41.resize(2);
    p41[0] = 0.0;
    p41[1] = 0.;

    p42.resize(2);
    p42[0] = Centre_x;
    p42[1] = 0.;

    p43.resize(2);
    p43[0] = 0.155596 * length;
    p43[1] = 0.;

    p44.resize(2);
    p44[0] = 0.183596 * length;
    p44[1] = 0.;

    p45.resize(2);
    p45[0] = 0.239596 * length;
    p45[1] = 0.;

    p46.resize(2);
    p46[0] = 0.285123967 * length;
    p46[1] = 0.;

    p47.resize(2);
    p47[0] = 0.433884298 * length;
    p47[1] = 0.;

    p48.resize(2);
    p48[0] = 0.578512397 * length;
    p48[1] = 0.;

    p49.resize(2);
    p49[0] = 0.789256198 * length;
    p49[1] = 0.;

    p50.resize(2);
    p50[0] = length;
    p50[1] = 0.;


    // Allocate storage for variable points
    p21.resize(2);
    p22.resize(2);
    p23.resize(2);
    p24.resize(2);
    p25.resize(2);
    p27.resize(2);
    p28.resize(2);
    p29.resize(2);
    p30.resize(2);
    p31.resize(2);

    // There are 31 macro elements
    Macro_element_pt.resize(31);

    // Build the 2D macro elements
    for (unsigned i = 0; i < 31; i++)
    {
      Macro_element_pt[i] = new QMacroElement<2>(this, i);
    }

  } // end of constructor


  //============================================================================
  /// Parametrisation of macro element boundaries: f(s) is the position
  /// vector to macro-element m's boundary in the specified direction [N/S/E/W]
  /// at the specfied discrete time level (time=0: present; time>0: previous)
  //============================================================================
  void CylinderWithFlagDomain::macro_element_boundary(const unsigned& time,
                                                      const unsigned& m,
                                                      const unsigned& direction,
                                                      const Vector<double>& s,
                                                      Vector<double>& f)
  {
    // Use Quadtree names for directions
    using namespace QuadTreeNames;

#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "CylinderWithFlagDomain::macro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif


    // Lagrangian coordinate along surface of cylinder
    Vector<double> xi(1);

    // Point on circle
    Vector<double> point_on_circle(2);

    // Lagrangian coordinates on the flag
    Vector<double> zeta(1);


    // Definition of the points that depend on the shape of the flags
    zeta[0] = 1. / 5. * Lx;
    Top_flag_pt->position(time, zeta, p21);

    zeta[0] = 2. / 5. * Lx;
    Top_flag_pt->position(time, zeta, p22);

    zeta[0] = 3. / 5. * Lx;
    Top_flag_pt->position(time, zeta, p23);

    zeta[0] = 4. / 5. * Lx;
    Top_flag_pt->position(time, zeta, p24);

    zeta[0] = Ly / 2.;
    Tip_flag_pt->position(time, zeta, p25);

    zeta[0] = 1. / 5. * Lx;
    Bottom_flag_pt->position(time, zeta, p27);

    zeta[0] = 2. / 5. * Lx;
    Bottom_flag_pt->position(time, zeta, p28);

    zeta[0] = 3. / 5. * Lx;
    Bottom_flag_pt->position(time, zeta, p29);

    zeta[0] = 4. / 5. * Lx;
    Bottom_flag_pt->position(time, zeta, p30);

    zeta[0] = -Ly / 2.;
    Tip_flag_pt->position(time, zeta, p31);


    std::ostringstream error_message;

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
            linear_interpolate(p1, point_on_circle, s[0], f);
            break;

          case S:
            xi[0] = -3.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            linear_interpolate(p41, point_on_circle, s[0], f);
            break;

          case W:
            linear_interpolate(p41, p1, s[0], f);
            break;

          case E:
            xi[0] = 5.0 * atan(1.0) - 2.0 * atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 1, is immediately above the cylinder
      case 1:

        switch (direction)
        {
          case N:
            linear_interpolate(p1, p2, s[0], f);
            break;

          case S:
            xi[0] = 3.0 * atan(1.0) - atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case W:
            xi[0] = 3.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            linear_interpolate(point_on_circle, p1, s[0], f);
            break;

          case E:
            linear_interpolate(p18, p2, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 2, is immediately right of the cylinder
      case 2:

        switch (direction)
        {
          case N:
            linear_interpolate(p2, p11, s[0], f);
            break;

          case S:
            xi[0] = 2.0 * atan(1.0) - atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case W:
            linear_interpolate(p18, p2, s[0], f);
            break;

          case E:
            xi[0] = atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            linear_interpolate(point_on_circle, p11, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 3, is immediately below cylinder
      case 3:

        switch (direction)
        {
          case N:
            xi[0] = atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            linear_interpolate(point_on_circle, p11, s[0], f);
            break;

          case S:
            xi[0] = (1. + s[0]) / 2. * 1. / 5. * Lx;
            Top_flag_pt->position(time, xi, f);
            break;

          case W:
            xi[0] = asin(Ly / A / 2.) +
                    (atan(1.0) - asin(Ly / A / 2.)) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case E:
            linear_interpolate(p21, p11, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 4, is right hand block 1
      case 4:

        switch (direction)
        {
          case N:
            xi[0] = (1. + s[0]) / 2. * 1. / 5. * Lx;
            Bottom_flag_pt->position(time, xi, f);
            break;

          case S:
            xi[0] = -atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            linear_interpolate(point_on_circle, p37, s[0], f);
            break;

          case W:
            xi[0] =
              -atan(1.0) + (atan(1.0) - asin(Ly / A / 2.)) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case E:
            linear_interpolate(p37, p27, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 5, is right hand block 2
      case 5:

        switch (direction)
        {
          case N:
            xi[0] = 6 * atan(1.0) + atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case S:
            linear_interpolate(p42, p37, s[0], f);
            break;

          case W:
            linear_interpolate(p42, p33, s[0], f);
            break;

          case E:
            xi[0] = -atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            linear_interpolate(p37, point_on_circle, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 6, is right hand block 3
      case 6:

        switch (direction)
        {
          case N:
            xi[0] = 5.0 * atan(1.0) + atan(1.0) * 0.5 * (1.0 + s[0]);
            Cylinder_pt->position(time, xi, f);
            break;

          case S:
            linear_interpolate(p41, p42, s[0], f);
            break;

          case W:
            xi[0] = 5.0 * atan(1.0);
            Cylinder_pt->position(time, xi, point_on_circle);
            linear_interpolate(p41, point_on_circle, s[0], f);
            break;

          case E:
            linear_interpolate(p42, p33, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 7, is right hand block 4
      case 7:

        switch (direction)
        {
          case N:
            linear_interpolate(p2, p3, s[0], f);
            break;

          case S:
            linear_interpolate(p11, p12, s[0], f);
            break;

          case W:
            linear_interpolate(p11, p2, s[0], f);
            break;

          case E:
            linear_interpolate(p12, p3, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 8:

        switch (direction)
        {
          case N:
            linear_interpolate(p11, p12, s[0], f);
            break;

          case S:
            xi[0] = 1. / 5. * Lx + (1. + s[0]) / 2. * 1. / 5. * Lx;
            Top_flag_pt->position(time, xi, f);
            break;

          case W:
            linear_interpolate(p21, p11, s[0], f);
            break;

          case E:
            linear_interpolate(p22, p12, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;
      case 9:

        switch (direction)
        {
          case N:
            xi[0] = 1. / 5. * Lx + (1. + s[0]) / 2. * 1. / 5. * Lx;
            Bottom_flag_pt->position(time, xi, f);
            break;

          case S:
            linear_interpolate(p37, p38, s[0], f);
            break;

          case W:
            linear_interpolate(p37, p27, s[0], f);
            break;

          case E:
            linear_interpolate(p38, p28, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 10:

        switch (direction)
        {
          case N:
            linear_interpolate(p37, p38, s[0], f);
            break;

          case S:
            linear_interpolate(p42, p43, s[0], f);
            break;

          case W:
            linear_interpolate(p42, p37, s[0], f);
            break;

          case E:
            linear_interpolate(p43, p38, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;
      case 11:

        switch (direction)
        {
          case N:
            linear_interpolate(p3, p4, s[0], f);
            break;

          case S:
            linear_interpolate(p12, p13, s[0], f);
            break;

          case W:
            linear_interpolate(p12, p3, s[0], f);
            break;

          case E:
            linear_interpolate(p13, p4, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 12:

        switch (direction)
        {
          case N:
            linear_interpolate(p12, p13, s[0], f);
            break;

          case S:
            //    linear_interpolate(p22,p23,s[0],f);
            xi[0] = 2. / 5. * Lx + (1. + s[0]) / 2. * 1. / 5. * Lx;
            Top_flag_pt->position(time, xi, f);
            break;

          case W:
            linear_interpolate(p22, p12, s[0], f);
            break;

          case E:
            linear_interpolate(p23, p13, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 13:

        switch (direction)
        {
          case N:
            xi[0] = 2. / 5. * Lx + (1. + s[0]) / 2. * 1. / 5. * Lx;
            Bottom_flag_pt->position(time, xi, f);
            break;

          case S:
            linear_interpolate(p38, p39, s[0], f);
            break;

          case W:
            linear_interpolate(p38, p28, s[0], f);
            break;

          case E:
            linear_interpolate(p39, p29, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 14:

        switch (direction)
        {
          case N:
            linear_interpolate(p38, p39, s[0], f);
            break;

          case S:
            linear_interpolate(p43, p44, s[0], f);
            break;

          case W:
            linear_interpolate(p43, p38, s[0], f);
            break;

          case E:
            linear_interpolate(p44, p39, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 15:

        switch (direction)
        {
          case N:
            linear_interpolate(p4, p5, s[0], f);
            break;

          case S:
            linear_interpolate(p13, p14, s[0], f);
            break;

          case W:
            linear_interpolate(p13, p4, s[0], f);
            break;

          case E:
            linear_interpolate(p14, p5, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 16:

        switch (direction)
        {
          case N:
            linear_interpolate(p13, p14, s[0], f);
            break;

          case S:
            xi[0] = 3. / 5. * Lx + (1. + s[0]) / 2. * 1. / 5. * Lx;
            Top_flag_pt->position(time, xi, f);
            break;

          case W:
            linear_interpolate(p23, p13, s[0], f);
            break;

          case E:
            linear_interpolate(p24, p14, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 17:

        switch (direction)
        {
          case N:
            xi[0] = 3. / 5. * Lx + (1. + s[0]) / 2. * 1. / 5. * Lx;
            Bottom_flag_pt->position(time, xi, f);
            break;

          case S:
            linear_interpolate(p39, p40, s[0], f);
            break;

          case W:
            linear_interpolate(p39, p29, s[0], f);
            break;

          case E:
            linear_interpolate(p40, p30, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 18:

        switch (direction)
        {
          case N:
            linear_interpolate(p39, p40, s[0], f);
            break;

          case S:
            linear_interpolate(p44, p45, s[0], f);
            break;

          case W:
            linear_interpolate(p44, p39, s[0], f);
            break;

          case E:
            linear_interpolate(p45, p40, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 19:

        switch (direction)
        {
          case N:
            linear_interpolate(p14, p5, s[0], f);
            break;

          case S:
            xi[0] = 4. / 5. * Lx + (1. + s[0]) / 2. * 1. / 5. * Lx;
            Top_flag_pt->position(time, xi, f);
            break;

          case W:
            linear_interpolate(p24, p14, s[0], f);
            break;

          case E:
            linear_interpolate(p25, p5, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 20:

        switch (direction)
        {
          case N:
            xi[0] = 4. / 5. * Lx + (1. + s[0]) / 2. * 1. / 5. * Lx;
            Bottom_flag_pt->position(time, xi, f);
            break;

          case S:
            linear_interpolate(p40, p45, s[0], f);
            break;

          case W:
            linear_interpolate(p40, p30, s[0], f);
            break;

          case E:
            linear_interpolate(p45, p31, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 21:

        switch (direction)
        {
          case N:
            linear_interpolate(p5, p6, s[0], f);
            break;

          case S:
            linear_interpolate(p25, p15, s[0], f);
            break;

          case W:
            linear_interpolate(p25, p5, s[0], f);
            break;

          case E:
            linear_interpolate(p15, p6, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 22:

        switch (direction)
        {
          case N:
            linear_interpolate(p25, p15, s[0], f);
            break;

          case S:
            linear_interpolate(p31, p35, s[0], f);
            break;

          case W:
            xi[0] = s[0] * Ly / 2.;
            Tip_flag_pt->position(time, xi, f);
            break;

          case E:
            linear_interpolate(p35, p15, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 23:

        switch (direction)
        {
          case N:
            linear_interpolate(p31, p35, s[0], f);
            break;

          case S:
            linear_interpolate(p45, p46, s[0], f);
            break;

          case W:
            linear_interpolate(p45, p31, s[0], f);
            break;

          case E:
            linear_interpolate(p46, p35, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 24:

        switch (direction)
        {
          case N:
            linear_interpolate(p6, p7, s[0], f);
            break;

          case S:
            linear_interpolate(p15, p16, s[0], f);
            break;

          case W:
            linear_interpolate(p15, p6, s[0], f);
            break;

          case E:
            linear_interpolate(p16, p7, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 25:

        switch (direction)
        {
          case N:
            linear_interpolate(p15, p16, s[0], f);
            break;

          case S:
            linear_interpolate(p35, p36, s[0], f);
            break;

          case W:
            linear_interpolate(p35, p15, s[0], f);
            break;

          case E:
            linear_interpolate(p36, p16, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 26:

        switch (direction)
        {
          case N:
            linear_interpolate(p35, p36, s[0], f);
            break;

          case S:
            linear_interpolate(p46, p47, s[0], f);
            break;

          case W:
            linear_interpolate(p46, p35, s[0], f);
            break;

          case E:
            linear_interpolate(p47, p36, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 27:

        switch (direction)
        {
          case N:
            linear_interpolate(p16, p7, s[0], f);
            break;

          case S:
            linear_interpolate(p36, p47, s[0], f);
            break;

          case W:
            linear_interpolate(p36, p16, s[0], f);
            break;

          case E:
            linear_interpolate(p47, p7, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 28:

        switch (direction)
        {
          case N:
            linear_interpolate(p7, p8, s[0], f);
            break;

          case S:
            linear_interpolate(p47, p48, s[0], f);
            break;

          case W:
            linear_interpolate(p47, p7, s[0], f);
            break;

          case E:
            linear_interpolate(p48, p8, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 29:

        switch (direction)
        {
          case N:
            linear_interpolate(p8, p9, s[0], f);
            break;

          case S:
            linear_interpolate(p48, p49, s[0], f);
            break;

          case W:
            linear_interpolate(p48, p8, s[0], f);
            break;

          case E:
            linear_interpolate(p49, p9, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      case 30:

        switch (direction)
        {
          case N:
            linear_interpolate(p9, p10, s[0], f);
            break;

          case S:
            linear_interpolate(p49, p50, s[0], f);
            break;

          case W:
            linear_interpolate(p49, p9, s[0], f);
            break;

          case E:
            linear_interpolate(p50, p10, s[0], f);
            break;

          default:
            error_message << "Direction is incorrect: " << direction
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

      default:

        error_message << "Wrong macro element number" << m << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }

  } // end of macro element boundary

} // namespace oomph
