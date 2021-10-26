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
#ifndef OOMPH_ANNULAR_DOMAIN_HEADER
#define OOMPH_ANNULAR_DOMAIN_HEADER

// Generic oomph-lib includes
#include "../generic/quadtree.h"
#include "../generic/domain.h"
#include "../generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  /// Annular domain
  //=================================================================
  class AnnularDomain : public Domain
  {
  public:
    /// Constructor: Specify azimuthal fraction (1.0 is 360 degrees)
    /// number of macro elements in azimuthal and radial direction,
    /// inner radius and thickness. Rotate mesh by angle phi.
    AnnularDomain(const double& azimuthal_fraction,
                  const unsigned& ntheta,
                  const unsigned& nr,
                  const double& a,
                  const double& h,
                  const double& phi)
      : Azimuthal_fraction(azimuthal_fraction),
        Inner_radius(a),
        Thickness(h),
        Ntheta(ntheta),
        Nr(nr),
        Phi(phi)
    {
      const unsigned n_macro = ntheta * nr;
      Macro_element_pt.resize(n_macro);

      // Create the macro elements
      for (unsigned i = 0; i < n_macro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    }

    /// Broken copy constructor
    AnnularDomain(const AnnularDomain&) = delete;

    /// Broken assignment operator
    void operator=(const AnnularDomain&) = delete;

    /// Destructor: Empty; cleanup done in base class
    ~AnnularDomain() {}


    /// Vector representation of the  i_macro-th macro element
    /// boundary i_direct (N/S/W/E) at time level t
    /// (t=0: present; t>0: previous):
    /// f(s).
    void macro_element_boundary(const unsigned& t,
                                const unsigned& i_macro,
                                const unsigned& i_direct,
                                const Vector<double>& s,
                                Vector<double>& f);

  private:
    /// Azimuthal fraction
    double Azimuthal_fraction;

    /// Inner radius
    double Inner_radius;

    /// Thickness
    double Thickness;

    /// Number of macro elements in azimuthal direction
    unsigned Ntheta;

    /// Number of macro elements in radial direction
    unsigned Nr;

    /// Rotation angle
    double Phi;
  };


  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Vector representation of the  imacro-th macro element
  /// boundary idirect (N/S/W/E) at time level t
  /// (t=0: present; t>0: previous): f(s)
  //=================================================================
  void AnnularDomain::macro_element_boundary(const unsigned& t,
                                             const unsigned& imacro,
                                             const unsigned& idirect,
                                             const Vector<double>& s,
                                             Vector<double>& f)
  {
    using namespace QuadTreeNames;

    // Get coordinates of macro element
    unsigned i_theta = imacro % Ntheta;
    unsigned i_r = (imacro - i_theta) / Ntheta;

    // Angle and radius limits
    double theta_lo = Azimuthal_fraction * 2.0 * MathematicalConstants::Pi *
                      double(i_theta) / double(Ntheta);

    double theta_hi = Azimuthal_fraction * 2.0 * MathematicalConstants::Pi *
                      double(i_theta + 1) / double(Ntheta);

    // Revert direction (convoluted -- don't ask. It mirrors what happens
    // in the mesh...
    theta_lo = -MathematicalConstants::Pi +
               Azimuthal_fraction * 2.0 * MathematicalConstants::Pi - theta_lo;
    theta_hi = -MathematicalConstants::Pi +
               Azimuthal_fraction * 2.0 * MathematicalConstants::Pi - theta_hi;

    double r_lo = Inner_radius + Thickness * double(i_r) / double(Nr);
    double r_hi = Inner_radius + Thickness * double(i_r + 1) / double(Nr);

    // Actual radius and angle
    double r = 0.0;
    double theta = 0.0;

    // Which direction?
    switch (idirect)
    {
      case N:

        theta = theta_lo + 0.5 * (s[0] + 1.0) * (theta_hi - theta_lo);
        r = r_hi;

        break;

      case S:

        theta = theta_lo + 0.5 * (s[0] + 1.0) * (theta_hi - theta_lo);
        r = r_lo;

        break;

      case W:

        theta = theta_lo;
        r = r_lo + 0.5 * (s[0] + 1.0) * (r_hi - r_lo);

        break;

      case E:

        theta = theta_hi;
        r = r_lo + 0.5 * (s[0] + 1.0) * (r_hi - r_lo);

        break;

      default:
        std::ostringstream error_stream;
        error_stream << "idirect is " << idirect << " not one of N, S, W, E"
                     << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    f[0] = r * cos(theta + Phi);
    f[1] = r * sin(theta + Phi);
  }

} // namespace oomph

#endif
