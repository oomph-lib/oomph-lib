// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
// Include guards
#ifndef OOMPH_QUARTER_PIPE_DOMAIN_HEADER
#define OOMPH_QUARTER_PIPE_DOMAIN_HEADER

// Generic oomph-lib includes
#include "../generic/octree.h"
#include "../generic/domain.h"
#include "../generic/geom_objects.h"

namespace oomph
{
  //================================================================
  /// Domain representing a quarter pipe
  //================================================================
  class QuarterPipeDomain : public Domain
  {
  public:
    /// Constructor: Pass number of elements in various directions,
    /// the inner and outer radius and the length of the tube
    QuarterPipeDomain(const unsigned& ntheta,
                      const unsigned& nr,
                      const unsigned& nz,
                      const double& rmin,
                      const double& rmax,
                      const double& length)
      : Ntheta(ntheta),
        Nr(nr),
        Nz(nz),
        Rmin(rmin),
        Rmax(rmax),
        Length(length),
        Axial_spacing_fct_pt(&default_axial_spacing_fct)
    {
      // Number of macroelements
      unsigned nmacro = nr * ntheta * nz;

      // Resize
      Macro_element_pt.resize(nmacro);

      // Create macro elements
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<3>(this, i);
      }

      // Make geom object representing the outer and inner boundaries of
      // the cross section
      Inner_boundary_cross_section_pt = new Ellipse(rmin, rmin);
      Outer_boundary_cross_section_pt = new Ellipse(rmax, rmax);
    }

    /// Broken copy constructor
    QuarterPipeDomain(const QuarterPipeDomain&) = delete;

    /// Broken assignment operator
    void operator=(const QuarterPipeDomain&) = delete;

    /// Destructor:
    ~QuarterPipeDomain()
    {
      // Note: macro elements are cleaned up in base class...
      delete Outer_boundary_cross_section_pt;
      delete Inner_boundary_cross_section_pt;
    }

    /// Typedef for function pointer for function that implements
    /// axial spacing of macro elements
    typedef double (*AxialSpacingFctPt)(const double& xi);

    /// Function pointer for function that  implements
    /// axial spacing of macro elements
    AxialSpacingFctPt& axial_spacing_fct_pt()
    {
      return Axial_spacing_fct_pt;
    }

    /// Function that implements
    /// axial spacing of macro elements
    double axial_spacing_fct(const double& xi)
    {
      return Axial_spacing_fct_pt(xi);
    }

    /// Vector representation of the i_macro-th macro element
    /// boundary i_direct (U/D/L/R/F/B) at time level t
    /// (t=0: present; t>0: previous): f(s).
    void macro_element_boundary(const unsigned& t,
                                const unsigned& i_macro,
                                const unsigned& i_direct,
                                const Vector<double>& s,
                                Vector<double>& f);

  private:
    /// Number of elements azimuthal direction
    unsigned Ntheta;

    /// Number of elements radial direction
    unsigned Nr;

    /// Number of elements axial direction
    unsigned Nz;

    /// Inner radius
    double Rmin;

    /// Outer radius
    double Rmax;

    /// Length
    double Length;

    /// Geom object representing the outer boundary of
    /// the cross section
    GeomObject* Outer_boundary_cross_section_pt;

    /// Geom object representing the inner boundary of
    /// the cross section
    GeomObject* Inner_boundary_cross_section_pt;

    /// Function pointer for function that implements
    /// axial spacing of macro elements
    AxialSpacingFctPt Axial_spacing_fct_pt;

    /// Default for function that  implements
    /// axial spacing of macro elements
    static double default_axial_spacing_fct(const double& xi)
    {
      return xi;
    }

    /// Boundary of macro element zeta \f$ \in [-1,1]x[-1,1] \f$
    void r_U(const unsigned& t,
             const Vector<double>& zeta,
             Vector<double>& f,
             const double& rmin,
             const double& rmax,
             const double& thetamin,
             const double& thetamax,
             const double& zmin,
             const double& zmax);

    /// Boundary of macro element zeta \f$ \in [-1,1]x[-1,1] \f$
    void r_L(const unsigned& t,
             const Vector<double>& zeta,
             Vector<double>& f,
             const double& rmin,
             const double& rmax,
             const double& thetamin,
             const double& thetamax,
             const double& zmin,
             const double& zmax);

    /// Boundary of macro element zeta \f$ \in [-1,1]x[-1,1] \f$
    void r_D(const unsigned& t,
             const Vector<double>& zeta,
             Vector<double>& f,
             const double& rmin,
             const double& rmax,
             const double& thetamin,
             const double& thetamax,
             const double& zmin,
             const double& zmax);

    /// Boundary of macro element zeta \f$ \in [-1,1]x[-1,1] \f$
    void r_R(const unsigned& t,
             const Vector<double>& zeta,
             Vector<double>& f,
             const double& rmin,
             const double& rmax,
             const double& thetamin,
             const double& thetamax,
             const double& zmin,
             const double& zmax);

    /// Boundary of macro element zeta \f$ \in [-1,1]x[-1,1] \f$
    void r_F(const unsigned& t,
             const Vector<double>& zeta,
             Vector<double>& f,
             const double& rmin,
             const double& rmax,
             const double& thetamin,
             const double& thetamax,
             const double& zmin,
             const double& zmax);

    /// Boundary of macro element zeta \f$ \in [-1,1]x[-1,1] \f$
    void r_B(const unsigned& t,
             const Vector<double>& zeta,
             Vector<double>& f,
             const double& rmin,
             const double& rmax,
             const double& thetamin,
             const double& thetamax,
             const double& zmin,
             const double& zmax);

  }; // endofclass


  //=================================================================
  /// Vector representation of the  imacro-th macro element
  /// boundary idirect (U/D/L/R/F/B) at time level t:
  /// f(s)
  //=================================================================
  void QuarterPipeDomain::macro_element_boundary(const unsigned& t,
                                                 const unsigned& imacro,
                                                 const unsigned& idirect,
                                                 const Vector<double>& s,
                                                 Vector<double>& f)
  {
    using namespace OcTreeNames;

    const double pi = MathematicalConstants::Pi;

    // Match the elements number with the position
    unsigned num_z = imacro / (Nr * Ntheta);
    unsigned num_y = (imacro % (Nr * Ntheta)) / Ntheta;
    unsigned num_x = imacro % Ntheta;

    // Define the extreme coordinates

    // radial direction
    double rmin = Rmin + (Rmax - Rmin) * double(num_y) / double(Nr);
    double rmax = Rmin + (Rmax - Rmin) * double(num_y + 1) / double(Nr);

    // theta direction
    double thetamin = (pi / 2.0) * (1.0 - double(num_x + 1) / double(Ntheta));
    double thetamax = (pi / 2.0) * (1.0 - double(num_x) / double(Ntheta));

    // zdirection (tube)
    double zmin = double(num_z) * Length / double(Nz);
    double zmax = double(num_z + 1) * Length / double(Nz);


    // Which direction?
    if (idirect == U)
    {
      r_U(t, s, f, rmin, rmax, thetamin, thetamax, zmin, zmax);
    }
    else if (idirect == D)
    {
      r_D(t, s, f, rmin, rmax, thetamin, thetamax, zmin, zmax);
    }
    else if (idirect == L)
    {
      r_L(t, s, f, rmin, rmax, thetamin, thetamax, zmin, zmax);
    }
    else if (idirect == R)
    {
      r_R(t, s, f, rmin, rmax, thetamin, thetamax, zmin, zmax);
    }
    else if (idirect == F)
    {
      r_F(t, s, f, rmin, rmax, thetamin, thetamax, zmin, zmax);
    }
    else if (idirect == B)
    {
      r_B(t, s, f, rmin, rmax, thetamin, thetamax, zmin, zmax);
    }
    else
    {
      std::ostringstream error_stream;
      error_stream << "idirect is " << idirect << " not one of U, D, L, R, F, B"
                   << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Now redistribute points in the axial direction
    double z_frac = f[2] / Length;
    f[2] = Length * axial_spacing_fct(z_frac);
  }


  //=================================================================
  /// Left face of a macro element \f$ s \in [-1,1]*[-1,1] \f$
  //=================================================================
  void QuarterPipeDomain::r_L(const unsigned& t,
                              const Vector<double>& s,
                              Vector<double>& f,
                              const double& rmin,
                              const double& rmax,
                              const double& thetamin,
                              const double& thetamax,
                              const double& zmin,
                              const double& zmax)
  {
    Vector<double> x(1);
    x[0] = thetamax;

    // Point on outer wall
    Vector<double> r_outer(2);
    Outer_boundary_cross_section_pt->position(t, x, r_outer);

    // Point on inner wall
    Vector<double> r_inner(2);
    Inner_boundary_cross_section_pt->position(t, x, r_inner);

    // Get layer boundaries
    Vector<double> r_top(2);
    Vector<double> r_bot(2);
    for (unsigned i = 0; i < 2; i++)
    {
      r_top[i] =
        r_inner[i] + (r_outer[i] - r_inner[i]) * (rmax - Rmin) / (Rmax - Rmin);
      r_bot[i] =
        r_inner[i] + (r_outer[i] - r_inner[i]) * (rmin - Rmin) / (Rmax - Rmin);
    }

    // Compute coordinates
    f[0] = r_bot[0] + (0.5 * (s[0] + 1.0)) * (r_top[0] - r_bot[0]);
    f[1] = r_bot[1] + (0.5 * (s[0] + 1.0)) * (r_top[1] - r_bot[1]);
    f[2] = zmin + (zmax - zmin) * (0.5 * (s[1] + 1.0));
  }

  //=================================================================
  /// Right face of a macro element \f$ s \in [-1,1]*[-1,1] \f$
  //=================================================================
  void QuarterPipeDomain::r_R(const unsigned& t,
                              const Vector<double>& s,
                              Vector<double>& f,
                              const double& rmin,
                              const double& rmax,
                              const double& thetamin,
                              const double& thetamax,
                              const double& zmin,
                              const double& zmax)
  {
    Vector<double> x(1);
    x[0] = thetamin;

    // Point on outer wall
    Vector<double> r_outer(2);
    Outer_boundary_cross_section_pt->position(t, x, r_outer);

    // Point on inner wall
    Vector<double> r_inner(2);
    Inner_boundary_cross_section_pt->position(t, x, r_inner);

    // Get layer boundaries
    Vector<double> r_top(2);
    Vector<double> r_bot(2);
    for (unsigned i = 0; i < 2; i++)
    {
      r_top[i] =
        r_inner[i] + (r_outer[i] - r_inner[i]) * (rmax - Rmin) / (Rmax - Rmin);
      r_bot[i] =
        r_inner[i] + (r_outer[i] - r_inner[i]) * (rmin - Rmin) / (Rmax - Rmin);
    }

    // Compute coordinates
    f[0] = r_bot[0] + (0.5 * (s[0] + 1.0)) * (r_top[0] - r_bot[0]);
    f[1] = r_bot[1] + (0.5 * (s[0] + 1.0)) * (r_top[1] - r_bot[1]);
    f[2] = zmin + (zmax - zmin) * (0.5 * (s[1] + 1.0));
  }


  //=================================================================
  /// Left face of a macro element \f$s \in [-1,1]*[-1,1] \f$
  //=================================================================
  void QuarterPipeDomain::r_D(const unsigned& t,
                              const Vector<double>& s,
                              Vector<double>& f,
                              const double& rmin,
                              const double& rmax,
                              const double& thetamin,
                              const double& thetamax,
                              const double& zmin,
                              const double& zmax)
  {
    Vector<double> x(1);
    x[0] = thetamax + (0.5 * (s[0] + 1.0)) * (thetamin - thetamax);

    // Point on outer wall
    Vector<double> r_outer(2);
    Outer_boundary_cross_section_pt->position(t, x, r_outer);

    // Point on inner wall
    Vector<double> r_inner(2);
    Inner_boundary_cross_section_pt->position(t, x, r_inner);

    // Get layer
    for (unsigned i = 0; i < 2; i++)
    {
      f[i] =
        r_inner[i] + (r_outer[i] - r_inner[i]) * (rmin - Rmin) / (Rmax - Rmin);
    }
    f[2] = zmin + (zmax - zmin) * (0.5 * (s[1] + 1.0));
  }


  //=================================================================
  /// Right face of a macro element \f$ s \in [-1,1]*[-1,1] \f$
  //=================================================================
  void QuarterPipeDomain::r_U(const unsigned& t,
                              const Vector<double>& s,
                              Vector<double>& f,
                              const double& rmin,
                              const double& rmax,
                              const double& thetamin,
                              const double& thetamax,
                              const double& zmin,
                              const double& zmax)
  {
    Vector<double> x(1);
    x[0] = thetamax + (0.5 * (s[0] + 1.0)) * (thetamin - thetamax);

    // Point on outer wall
    Vector<double> r_outer(2);
    Outer_boundary_cross_section_pt->position(t, x, r_outer);

    // Point on inner wall
    Vector<double> r_inner(2);
    Inner_boundary_cross_section_pt->position(t, x, r_inner);

    // Get layer
    for (unsigned i = 0; i < 2; i++)
    {
      f[i] =
        r_inner[i] + (r_outer[i] - r_inner[i]) * (rmax - Rmin) / (Rmax - Rmin);
    }
    f[2] = zmin + (zmax - zmin) * (0.5 * (s[1] + 1.0));
  }


  //=================================================================
  /// Front face of a macro element \f$ s \in [-1,1]*[-1,1] \f$
  //=================================================================
  void QuarterPipeDomain::r_F(const unsigned& t,
                              const Vector<double>& s,
                              Vector<double>& f,
                              const double& rmin,
                              const double& rmax,
                              const double& thetamin,
                              const double& thetamax,
                              const double& zmin,
                              const double& zmax)
  {
    Vector<double> x(1);
    x[0] = thetamax + (0.5 * (s[0] + 1.0)) * (thetamin - thetamax);

    // Point on outer wall
    Vector<double> r_outer(2);
    Outer_boundary_cross_section_pt->position(t, x, r_outer);

    // Point on inner wall
    Vector<double> r_inner(2);
    Inner_boundary_cross_section_pt->position(t, x, r_inner);

    // Get layer
    double rad = rmin + (0.5 * (s[1] + 1.0)) * (rmax - rmin);
    for (unsigned i = 0; i < 2; i++)
    {
      f[i] =
        r_inner[i] + (r_outer[i] - r_inner[i]) * (rad - Rmin) / (Rmax - Rmin);
    }
    f[2] = zmax;
  }


  //=================================================================
  /// Back face of a macro element \f$ s \in [-1,1]*[-1,1] \f$
  //=================================================================
  void QuarterPipeDomain::r_B(const unsigned& t,
                              const Vector<double>& s,
                              Vector<double>& f,
                              const double& rmin,
                              const double& rmax,
                              const double& thetamin,
                              const double& thetamax,
                              const double& zmin,
                              const double& zmax)
  {
    Vector<double> x(1);
    x[0] = thetamax + (0.5 * (s[0] + 1.0)) * (thetamin - thetamax);

    // Point on outer wall
    Vector<double> r_outer(2);
    Outer_boundary_cross_section_pt->position(t, x, r_outer);

    // Point on inner wall
    Vector<double> r_inner(2);
    Inner_boundary_cross_section_pt->position(t, x, r_inner);

    // Get layer
    double rad = rmin + (0.5 * (s[1] + 1.0)) * (rmax - rmin);
    for (unsigned i = 0; i < 2; i++)
    {
      f[i] =
        r_inner[i] + (r_outer[i] - r_inner[i]) * (rad - Rmin) / (Rmax - Rmin);
    }
    f[2] = zmin;
  }


} // namespace oomph

#endif
