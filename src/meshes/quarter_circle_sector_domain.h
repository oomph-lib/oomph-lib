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
// Include guards
#ifndef OOMPH_QUARTER_CIRCLE_SECTOR_DOMAIN_HEADER
#define OOMPH_QUARTER_CIRCLE_SECTOR_DOMAIN_HEADER

// Generic oomph-lib includes
#include "../generic/quadtree.h"
#include "../generic/domain.h"
#include "../generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  ///  Circular sector as domain. Domain is bounded by
  /// curved boundary which is represented by a GeomObject. Domain is
  /// parametrised by three macro elements.
  //=================================================================
  class QuarterCircleSectorDomain : public Domain
  {
  public:
    ///  Constructor: Pass boundary object and start and end coordinates
    /// and fraction along boundary object where outer ring is divided.
    QuarterCircleSectorDomain(GeomObject* boundary_geom_object_pt,
                              const double& xi_lo,
                              const double& fract_mid,
                              const double& xi_hi)
      : Xi_lo(xi_lo),
        Fract_mid(fract_mid),
        Xi_hi(xi_hi),
        Wall_pt(boundary_geom_object_pt),
        BL_squash_fct_pt(&default_BL_squash_fct)
    {
      // There are three macro elements
      unsigned nmacro = 3;

      // Resize
      Macro_element_pt.resize(nmacro);

      // Create macro elements
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    }


    /// Broken copy constructor
    QuarterCircleSectorDomain(const QuarterCircleSectorDomain&) = delete;

    /// Broken assignment operator
    void operator=(const QuarterCircleSectorDomain&) = delete;

    /// Destructor: empty; cleanup done in base class
    ~QuarterCircleSectorDomain() {}

    ///  Typedef for function pointer for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value
    typedef double (*BLSquashFctPt)(const double& s);


    ///  Function pointer for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value
    BLSquashFctPt& bl_squash_fct_pt()
    {
      return BL_squash_fct_pt;
    }


    ///  Function that squashes the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value
    double s_squashed(const double& s)
    {
      return BL_squash_fct_pt(s);
    }

    ///  Vector representation of the  i_macro-th macro element
    /// boundary i_direct (N/S/W/E) at time level t
    /// (t=0: present; t>0: previous):
    /// f(s). Note that the local coordinate \b s is a 1D
    /// Vector rather than a scalar -- this is unavoidable because
    /// this function implements the pure virtual function in the
    /// Domain base class.
    void macro_element_boundary(const unsigned& t,
                                const unsigned& i_macro,
                                const unsigned& i_direct,
                                const Vector<double>& s,
                                Vector<double>& f);


  private:
    /// Lower limit for the (1D) coordinates along the wall
    double Xi_lo;

    /// Fraction along wall where outer ring is to be divided
    double Fract_mid;

    /// Upper limit for the (1D) coordinates along the wall
    double Xi_hi;

    /// Pointer to geometric object that represents the curved wall
    GeomObject* Wall_pt;

    ///  Function pointer for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value
    BLSquashFctPt BL_squash_fct_pt;

    ///  Default for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value: Identity.
    static double default_BL_squash_fct(const double& s)
    {
      return s;
    }

    ///  Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_left_N(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    ///  Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_left_W(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    ///  Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_left_S(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    ///  Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_left_E(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    ///  Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_N(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    ///  Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_W(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    ///  Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_S(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    ///  Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_E(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    ///  Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_N(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_E(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_S(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_W(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);
  };


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Vector representation of the  imacro-th macro element
  /// boundary idirect (N/S/W/E) at time level t (t=0: present; t>0: previous):
  /// f(s)
  //=================================================================
  void QuarterCircleSectorDomain::macro_element_boundary(
    const unsigned& t,
    const unsigned& imacro,
    const unsigned& idirect,
    const Vector<double>& s,
    Vector<double>& f)
  {
#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "QuarterCircleSectorDomain::macro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    // Which macro element?
    // --------------------
    switch (imacro)
    {
      using namespace QuadTreeNames;

        // Macro element 0: Central box
      case 0:

        // Which direction?
        if (idirect == N)
        {
          r_centr_N(t, s, f);
        }
        else if (idirect == S)
        {
          r_centr_S(t, s, f);
        }
        else if (idirect == W)
        {
          r_centr_W(t, s, f);
        }
        else if (idirect == E)
        {
          r_centr_E(t, s, f);
        }
        else
        {
          std::ostringstream error_stream;
          error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                       << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 1: Bottom right
      case 1:

        // Which direction?
        if (idirect == N)
        {
          r_bot_right_N(t, s, f);
        }
        else if (idirect == S)
        {
          r_bot_right_S(t, s, f);
        }
        else if (idirect == W)
        {
          r_bot_right_W(t, s, f);
        }
        else if (idirect == E)
        {
          r_bot_right_E(t, s, f);
        }
        else
        {
          std::ostringstream error_stream;
          error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                       << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 2:Top left
      case 2:

        // Which direction?
        if (idirect == N)
        {
          r_top_left_N(t, s, f);
        }
        else if (idirect == S)
        {
          r_top_left_S(t, s, f);
        }
        else if (idirect == W)
        {
          r_top_left_W(t, s, f);
        }
        else if (idirect == E)
        {
          r_top_left_E(t, s, f);
        }
        else
        {
          std::ostringstream error_stream;
          error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                       << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        break;

      default:

        // Error
        std::ostringstream error_stream;
        error_stream << "Wrong imacro " << imacro << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //=================================================================
  /// Northern edge of top left macro element \f$ s \in [-1,1] \f$
  //=================================================================
  void QuarterCircleSectorDomain::r_top_left_N(const unsigned& t,
                                               const Vector<double>& s,
                                               Vector<double>& f)
  {
    Vector<double> x(1);

    // Coordinate along wall
    x[0] = Xi_hi + (Fract_mid * (Xi_hi - Xi_lo) - Xi_hi) * 0.5 * (s[0] + 1.0);

    Wall_pt->position(t, x, f);
  }


  //=================================================================
  /// Western edge of top left macro element \f$s \in [-1,1] \f$
  //=================================================================
  void QuarterCircleSectorDomain::r_top_left_W(const unsigned& t,
                                               const Vector<double>& s,
                                               Vector<double>& f)
  {
    Vector<double> x(1);

    // Top left corner
    Vector<double> r_top(2);
    x[0] = Xi_hi;

    Wall_pt->position(t, x, r_top);

    f[0] = 0.0;
    f[1] = 0.5 * r_top[1] * (1.0 + s_squashed(0.5 * (s[0] + 1.0)));
  }


  //=================================================================
  /// Southern edge of top left macro element \f$ s \in [-1,1] \f$
  //=================================================================
  void QuarterCircleSectorDomain::r_top_left_S(const unsigned& t,
                                               const Vector<double>& s,
                                               Vector<double>& f)
  {
    Vector<double> x(1);

    // Top left corner
    Vector<double> r_top(2);
    x[0] = Xi_hi;

    Wall_pt->position(t, x, r_top);


    // Bottom right corner
    Vector<double> r_bot(2);
    x[0] = 0.0;

    Wall_pt->position(t, x, r_bot);

    f[0] = 0.5 * r_bot[0] * 0.5 * (s[0] + 1.0);
    f[1] = 0.5 * r_top[1];
  }


  //=================================================================
  /// Eastern edge of top left macro element \f$ s \in [-1,1] \f$
  //=================================================================
  void QuarterCircleSectorDomain::r_top_left_E(const unsigned& t,
                                               const Vector<double>& s,
                                               Vector<double>& f)
  {
    Vector<double> x(1);

    // Top left corner
    Vector<double> r_top(2);
    x[0] = Xi_hi;

    Wall_pt->position(t, x, r_top);

    // Bottom right corner
    Vector<double> r_bot(2);
    x[0] = Xi_lo;

    Wall_pt->position(t, x, r_bot);

    // Halfway along wall
    Vector<double> r_half(2);
    x[0] = Xi_lo + Fract_mid * (Xi_hi - Xi_lo);

    Wall_pt->position(t, x, r_half);

    f[0] = 0.5 * (r_bot[0] + s_squashed(0.5 * (s[0] + 1.0)) *
                               (2.0 * r_half[0] - r_bot[0]));
    f[1] = 0.5 * (r_top[1] + s_squashed(0.5 * (s[0] + 1.0)) *
                               (2.0 * r_half[1] - r_top[1]));
  }


  //=================================================================
  /// Northern edge of bottom right macro element
  //=================================================================
  void QuarterCircleSectorDomain::r_bot_right_N(const unsigned& t,
                                                const Vector<double>& s,
                                                Vector<double>& f)
  {
    r_top_left_E(t, s, f);
  }

  //=================================================================
  /// Western edge of bottom right macro element
  //=================================================================
  void QuarterCircleSectorDomain::r_bot_right_W(const unsigned& t,
                                                const Vector<double>& s,
                                                Vector<double>& f)
  {
    Vector<double> x(1);

    // Top left corner
    Vector<double> r_top(2);
    x[0] = Xi_hi;

    Wall_pt->position(t, x, r_top);

    // Bottom right corner
    Vector<double> r_bot(2);
    x[0] = Xi_lo;

    Wall_pt->position(t, x, r_bot);

    f[0] = 0.5 * r_bot[0];
    f[1] = 0.5 * r_top[1] * 0.5 * (s[0] + 1.0);
  }

  //=================================================================
  /// Southern edge of bottom right macro element
  //=================================================================
  void QuarterCircleSectorDomain::r_bot_right_S(const unsigned& t,
                                                const Vector<double>& s,
                                                Vector<double>& f)
  {
    Vector<double> x(1);

    // Bottom right corner
    Vector<double> r_bot(2);
    x[0] = Xi_lo;
    Wall_pt->position(t, x, r_bot);


    f[0] = 0.5 * r_bot[0] * (1.0 + s_squashed(0.5 * (s[0] + 1.0)));
    f[1] = 0.0;
  }

  //=================================================================
  /// Eastern edge of bottom right macro element
  //=================================================================
  void QuarterCircleSectorDomain::r_bot_right_E(const unsigned& t,
                                                const Vector<double>& s,
                                                Vector<double>& f)
  {
    Vector<double> x(1);

    // Coordinate along wall
    x[0] = Xi_lo + (Fract_mid * (Xi_hi - Xi_lo) - Xi_lo) * (s[0] + 1.0) * 0.5;

    Wall_pt->position(t, x, f);
  }


  //=================================================================
  /// Northern edge of central box
  //=================================================================
  void QuarterCircleSectorDomain::r_centr_N(const unsigned& t,
                                            const Vector<double>& s,
                                            Vector<double>& f)
  {
    r_top_left_S(t, s, f);
  }


  //=================================================================
  /// Eastern edge of central box
  //=================================================================
  void QuarterCircleSectorDomain::r_centr_E(const unsigned& t,
                                            const Vector<double>& s,
                                            Vector<double>& f)
  {
    r_bot_right_W(t, s, f);
  }


  //=================================================================
  /// Southern edge of central box
  //=================================================================
  void QuarterCircleSectorDomain::r_centr_S(const unsigned& t,
                                            const Vector<double>& s,
                                            Vector<double>& f)
  {
    Vector<double> x(1);

    // Bottom right corner
    Vector<double> r_bot(2);
    x[0] = Xi_lo;
    Wall_pt->position(t, x, r_bot);

    f[0] = 0.5 * r_bot[0] * 0.5 * (s[0] + 1.0);
    f[1] = 0.0;
  }


  //=================================================================
  /// Western  edge of central box
  //=================================================================
  void QuarterCircleSectorDomain::r_centr_W(const unsigned& t,
                                            const Vector<double>& s,
                                            Vector<double>& f)
  {
    Vector<double> x(1);

    // Top left corner
    Vector<double> r_top(2);
    x[0] = Xi_hi;
    Wall_pt->position(t, x, r_top);

    f[0] = 0.0;
    f[1] = 0.5 * r_top[1] * 0.5 * (s[0] + 1.0);
  }


} // namespace oomph

#endif
