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
#ifndef OOMPH_CYLINDER_WITH_FLAG_DOMAIN_HEADER
#define OOMPH_CYLINDER_WITH_FLAG_DOMAIN_HEADER


// Generic includes
#include "../generic/geom_objects.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"


namespace oomph
{
  //===========================================================
  /// Domain for cylinder with flag as in Turek benchmark.
  //===========================================================
  class CylinderWithFlagDomain : public Domain
  {
  public:
    /// Constructor. Pass the pointers to the GeomObjects that parametrise
    /// the cylinder, the three edges of the flag, the length and height of the
    /// domain, the length and height of the flag, the coordinates of the
    /// centre of the cylinder and its radius.
    CylinderWithFlagDomain(Circle* cylinder_pt,
                           GeomObject* top_flag_pt,
                           GeomObject* bottom_flag_pt,
                           GeomObject* tip_flag_pt,
                           const double& length,
                           const double& height,
                           const double& flag_length,
                           const double& flag_height,
                           const double& centre_x,
                           const double& centre_y,
                           const double& a);


    /// Destructor: Emtpy because clean up happens in base class
    /// as a service to the user!
    ~CylinderWithFlagDomain() {}


    /// Parametrisation of macro element boundaries: f(s) is the position
    /// vector to macro-element m's boundary in the specified direction
    /// [N/S/E/W] at the specfied discrete time level (time=0: present; time>0:
    /// previous)
    void macro_element_boundary(const unsigned& time,
                                const unsigned& m,
                                const unsigned& direction,
                                const Vector<double>& s,
                                Vector<double>& f);

    /// Access fct to GeomObject (of type Circle)
    /// that represents the cylinder
    Circle* cylinder_pt()
    {
      return Cylinder_pt;
    }

    /// Access fct to GeomObjects for top, bottom and tip
    GeomObject*& bottom_flag_pt()
    {
      return Bottom_flag_pt;
    }
    GeomObject*& top_flag_pt()
    {
      return Top_flag_pt;
    }
    GeomObject*& tip_flag_pt()
    {
      return Tip_flag_pt;
    }

  private:
    /// Helper function to interpolate linearly between the
    /// "right" and "left" points; \f$ s \in [-1,1] \f$
    void linear_interpolate(const Vector<double>& left,
                            const Vector<double>& right,
                            const double& s,
                            Vector<double>& f)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        f[i] = left[i] + (right[i] - left[i]) * 0.5 * (s + 1.0);
      }
    }

    // Helper points
    Vector<double> p1;
    Vector<double> p2;
    Vector<double> p3;
    Vector<double> p4;
    Vector<double> p5;
    Vector<double> p6;
    Vector<double> p7;
    Vector<double> p8;
    Vector<double> p9;
    Vector<double> p10;
    Vector<double> p11;
    Vector<double> p12;
    Vector<double> p13;
    Vector<double> p14;
    Vector<double> p15;
    Vector<double> p16;
    Vector<double> p17;
    Vector<double> p18;
    Vector<double> p19;
    Vector<double> p20;
    Vector<double> p21;
    Vector<double> p22;
    Vector<double> p23;
    Vector<double> p24;
    Vector<double> p25;
    Vector<double> p26;
    Vector<double> p27;
    Vector<double> p28;
    Vector<double> p29;
    Vector<double> p30;
    Vector<double> p31;
    Vector<double> p32;
    Vector<double> p33;
    Vector<double> p34;
    Vector<double> p35;
    Vector<double> p36;
    Vector<double> p37;
    Vector<double> p38;
    Vector<double> p39;
    Vector<double> p40;
    Vector<double> p41;
    Vector<double> p42;
    Vector<double> p43;
    Vector<double> p44;
    Vector<double> p45;
    Vector<double> p46;
    Vector<double> p47;
    Vector<double> p48;
    Vector<double> p49;
    Vector<double> p50;


    /// Pointer to geometric object that represents the central cylinder
    Circle* Cylinder_pt;

    /// Pointer to geometric object that represents the top of the flag
    GeomObject* Top_flag_pt;

    /// Pointer to geometric object that represents the bottom of the flag
    GeomObject* Bottom_flag_pt;

    /// Pointer to geometric object that represents the tip of the flag
    GeomObject* Tip_flag_pt;

    // Length of the flag
    double Lx;

    // Thickness of the flag
    double Ly;

    // Centre of the cylinder : x coordinate
    double Centre_x;

    // Centre of the cylinder : y coordinate
    double Centre_y;

    // Radius of the cylinder
    double A;


  }; // end of domain

} // namespace oomph

#endif
