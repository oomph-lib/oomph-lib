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
#include "quadtree.h"


#include "macro_element.h"
#include "geom_objects.h"
#include "domain.h"

namespace oomph
{
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  // Warped cube domain
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Vector representation of the  imacro-th macro element
  /// boundary idirect (L/R/D/B/F) at time level t (t=0: present; t>0:
  /// previous): f(s)
  //=================================================================
  void WarpedCubeDomain::macro_element_boundary(const unsigned& t,
                                                const unsigned& imacro,
                                                const unsigned& idirect,
                                                const Vector<double>& s,
                                                Vector<double>& f)
  {
    using namespace OcTreeNames;

#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "WarpedCubeDomain::macro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    // Which direction?
    if (idirect == L)
    {
      r_L(t, s, f);
    }
    else if (idirect == R)
    {
      r_R(t, s, f);
    }
    else if (idirect == D)
    {
      r_D(t, s, f);
    }
    else if (idirect == U)
    {
      r_U(t, s, f);
    }
    else if (idirect == B)
    {
      r_B(t, s, f);
    }
    else if (idirect == F)
    {
      r_F(t, s, f);
    }
    else
    {
      std::ostringstream error_stream;
      error_stream << "idirect is " << idirect << " not one of U, D, L, R, B, F"
                   << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //#####################################################################


  //=======================================================================
  /// Left boundary face
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void WarpedCubeDomain::r_L(const unsigned& t,
                             const Vector<double>& zeta,
                             Vector<double>& f)
  {
    f[0] = -1.0;
    f[1] = zeta[0];
    f[2] = zeta[1];


    // Warp it
    warp_it(f);
  }


  //=======================================================================
  /// Right boundary face
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void WarpedCubeDomain::r_R(const unsigned& t,
                             const Vector<double>& zeta,
                             Vector<double>& f)
  {
    f[0] = 1.0;
    f[1] = zeta[0];
    f[2] = zeta[1];

    // Warp it
    warp_it(f);
  }


  //=======================================================================
  /// Down boundary face
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void WarpedCubeDomain::r_D(const unsigned& t,
                             const Vector<double>& zeta,
                             Vector<double>& f)
  {
    f[0] = zeta[0];
    f[1] = -1.0;
    f[2] = zeta[1];

    // Warp it
    warp_it(f);
  }


  //=======================================================================
  /// Up boundary face
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void WarpedCubeDomain::r_U(const unsigned& t,
                             const Vector<double>& zeta,
                             Vector<double>& f)
  {
    f[0] = zeta[0];
    f[1] = 1.0;
    f[2] = zeta[1];


    // Warp it
    warp_it(f);
  }


  //=======================================================================
  /// Back boundary face
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void WarpedCubeDomain::r_B(const unsigned& t,
                             const Vector<double>& zeta,
                             Vector<double>& f)
  {
    f[0] = zeta[0];
    f[1] = zeta[1];
    f[2] = -1.0;

    // Warp it
    warp_it(f);
  }


  //=======================================================================
  /// Front boundary face
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void WarpedCubeDomain::r_F(const unsigned& t,
                             const Vector<double>& zeta,
                             Vector<double>& f)
  {
    f[0] = zeta[0];
    f[1] = zeta[1];
    f[2] = 1.0;

    // Warp it
    warp_it(f);
  }


  //=======================================================================
  /// Warp the unit cube
  //=======================================================================
  void WarpedCubeDomain::warp_it(Vector<double>& f)
  {
    Vector<double> f_aux(f);
    double x = 0.5 * (1.0 + f_aux[0]);
    double y = 0.5 * (1.0 + f_aux[1]);
    double z = 0.5 * (1.0 + f_aux[2]);
    f[0] = (1.0 + x) * cos(y + 0.5 * z);
    f[1] = (2.0 + 3 * x) * sin(y + 0.5 * z);
    f[2] = sin(z) + 0.1 * (x * x + y * y);
  }


} // namespace oomph
