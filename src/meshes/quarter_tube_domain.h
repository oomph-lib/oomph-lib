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
// Include guards
#ifndef OOMPH_QUARTER_TUBE_DOMAIN_HEADER
#define OOMPH_QUARTER_TUBE_DOMAIN_HEADER

// Generic oomph-lib includes
#include "../generic/quadtree.h"
#include "../generic/domain.h"
#include "../generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  /// Quarter tube as domain. Domain is bounded by
  /// curved boundary which is represented by a GeomObject. Domain is
  /// parametrised by three macro elements in each of the nlayer slices
  //=================================================================
  class QuarterTubeDomain : public Domain
  {
  public:
    /// Constructor: Pass boundary object and start and end coordinates
    /// and fraction along boundary object where outer ring is divided.
    /// We form nlayer axial slices.
    QuarterTubeDomain(GeomObject* boundary_geom_object_pt,
                      const Vector<double>& xi_lo,
                      const double& fract_mid,
                      const Vector<double>& xi_hi,
                      const unsigned& nlayer)
      : Xi_lo(xi_lo),
        Fract_mid(fract_mid),
        Xi_hi(xi_hi),
        Nlayer(nlayer),
        Wall_pt(boundary_geom_object_pt),
        BL_squash_fct_pt(&default_BL_squash_fct),
        Axial_spacing_fct_pt(&default_axial_spacing_fct)
    {
      // There are three macro elements
      unsigned nmacro = 3 * nlayer;

      // Resize
      Macro_element_pt.resize(nmacro);

      // Create macro elements
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<3>(this, i);
      }
    }

    /// Broken copy constructor
    QuarterTubeDomain(const QuarterTubeDomain&) = delete;

    /// Broken assignment operator
    void operator=(const QuarterTubeDomain&) = delete;

    /// Destructor: empty; cleanup done in base class
    ~QuarterTubeDomain() {}

    /// Typedef for function pointer for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value
    typedef double (*BLSquashFctPt)(const double& s);


    /// Function pointer for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value
    BLSquashFctPt& bl_squash_fct_pt()
    {
      return BL_squash_fct_pt;
    }


    /// Function that squashes the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value.
    double s_squashed(const double& s)
    {
      return BL_squash_fct_pt(s);
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


    /// Vector representation of the  i_macro-th macro element
    /// boundary i_direct (L/R/D/U/B/F) at time level t
    /// (t=0: present; t>0: previous):
    /// f(s).
    void macro_element_boundary(const unsigned& t,
                                const unsigned& i_macro,
                                const unsigned& i_direct,
                                const Vector<double>& s,
                                Vector<double>& f);

  private:
    /// Lower limit for the coordinates along the wall
    Vector<double> Xi_lo;

    /// Fraction along wall where outer ring is to be divided
    double Fract_mid;

    /// Upper limit for the coordinates along the wall
    Vector<double> Xi_hi;

    /// Number of layers
    unsigned Nlayer;

    /// Pointer to geometric object that represents the curved wall
    GeomObject* Wall_pt;


    /// Function pointer for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value
    BLSquashFctPt BL_squash_fct_pt;


    /// Default for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value: Identity.
    static double default_BL_squash_fct(const double& s)
    {
      return s;
    }


    /// Function pointer for function that implements
    /// axial spacing of macro elements
    AxialSpacingFctPt Axial_spacing_fct_pt;


    /// Default for function that  implements
    /// axial spacing of macro elements
    static double default_axial_spacing_fct(const double& xi)
    {
      return xi;
    }


    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_L(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);

    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_R(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);

    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_D(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);

    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_U(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);

    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_B(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);

    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_F(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);


    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_L(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);

    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_R(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);

    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_D(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);

    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_U(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);

    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_B(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);

    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_F(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);


    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_L(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);

    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_R(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);

    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_D(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);

    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_U(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);

    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_B(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);

    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_F(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);
  };


  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Vector representation of the  imacro-th macro element
  /// boundary idirect (L/R/D/U/B/F) at time level t
  /// (t=0: present; t>0: previous): f(s)
  //=================================================================
  void QuarterTubeDomain::macro_element_boundary(const unsigned& t,
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
      "QuarterTubeDomain::macro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif


    unsigned ilayer = unsigned(imacro / 3);

    // Which macro element?
    // --------------------
    switch (imacro % 3)
    {
        // Macro element 0: Central box
      case 0:

        // Which direction?
        if (idirect == L)
        {
          r_centr_L(t, s, ilayer, f);
        }
        else if (idirect == R)
        {
          r_centr_R(t, s, ilayer, f);
        }
        else if (idirect == D)
        {
          r_centr_D(t, s, ilayer, f);
        }
        else if (idirect == U)
        {
          r_centr_U(t, s, ilayer, f);
        }
        else if (idirect == B)
        {
          r_centr_B(t, s, ilayer, f);
        }
        else if (idirect == F)
        {
          r_centr_F(t, s, ilayer, f);
        }
        else
        {
          std::ostringstream error_stream;
          error_stream << "idirect is " << idirect
                       << " not one of L, R, D, U, B, F" << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        break;


        // Macro element 1: Bottom right
      case 1:

        // Which direction?
        if (idirect == L)
        {
          r_bot_right_L(t, s, ilayer, f);
        }
        else if (idirect == R)
        {
          r_bot_right_R(t, s, ilayer, f);
        }
        else if (idirect == D)
        {
          r_bot_right_D(t, s, ilayer, f);
        }
        else if (idirect == U)
        {
          r_bot_right_U(t, s, ilayer, f);
        }
        else if (idirect == B)
        {
          r_bot_right_B(t, s, ilayer, f);
        }
        else if (idirect == F)
        {
          r_bot_right_F(t, s, ilayer, f);
        }
        else
        {
          std::ostringstream error_stream;
          error_stream << "idirect is " << idirect
                       << " not one of L, R, D, U, B, F" << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }


        break;

      // Macro element 2:Top left
      case 2:

        // Which direction?
        if (idirect == L)
        {
          r_top_left_L(t, s, ilayer, f);
        }
        else if (idirect == R)
        {
          r_top_left_R(t, s, ilayer, f);
        }
        else if (idirect == D)
        {
          r_top_left_D(t, s, ilayer, f);
        }
        else if (idirect == U)
        {
          r_top_left_U(t, s, ilayer, f);
        }
        else if (idirect == B)
        {
          r_top_left_B(t, s, ilayer, f);
        }
        else if (idirect == F)
        {
          r_top_left_F(t, s, ilayer, f);
        }
        else
        {
          std::ostringstream error_stream;
          error_stream << "idirect is " << idirect
                       << " not one of L, R, D, U, B, F" << std::endl;

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


  //=======================================================================
  /// Boundary of central box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_centr_L(const unsigned& t,
                                    const Vector<double>& zeta,
                                    const unsigned& i_layer,
                                    Vector<double>& f)
  {
    // Wall coordinates along top edge of wall
    Vector<double> x(2);
    x[0] =
      Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                   axial_spacing_fct((0.5 * (1.0 + zeta[1]) + double(i_layer)) /
                                     double(Nlayer));
    x[1] = Xi_hi[1];

    // Get position vector to upper edge of wall
    Vector<double> r_top(3);
    Wall_pt->position(t, x, r_top);

    // Scale it down to half the height
    f[0] = r_top[0];
    f[1] = r_top[1] * 0.25 * (1.0 + zeta[0]);
    // Warp it:
    double rho = 0.0; // 0.25*(1.0+zeta[0]);
    f[2] = x[0] + rho * (r_top[2] - x[0]);

    // f[2]=r_top[2];
  }


  //=======================================================================
  /// Boundary of central box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_centr_R(const unsigned& t,
                                    const Vector<double>& zeta,
                                    const unsigned& i_layer,
                                    Vector<double>& f)
  {
    // Note the repetition in the calculation, there is some scope
    // for optimisation

    // Wall coordinates along top edge of wall
    Vector<double> x(2);
    x[0] =
      Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                   axial_spacing_fct((0.5 * (1.0 + zeta[1]) + double(i_layer)) /
                                     double(Nlayer));
    x[1] = Xi_hi[1];

    // Get position vector to upper edge of wall
    Vector<double> r_top(3);
    Wall_pt->position(t, x, r_top);

    // Wall coordinates along bottom edge of wall
    x[0] =
      Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                   axial_spacing_fct((0.5 * (1.0 + zeta[1]) + double(i_layer)) /
                                     double(Nlayer));
    x[1] = Xi_lo[1];

    // Get position vector to bottom edge of wall
    Vector<double> r_bottom(3);
    Wall_pt->position(t, x, r_bottom);

    // Scale it down to half the height, halfway across width
    f[0] = 0.5 * r_bottom[0];
    f[1] = r_top[1] * 0.25 * (1.0 + zeta[0]);

    // Warp it:
    double rho = 0.0; // 0.25*(1.0+zeta[0]);
    f[2] = x[0] + rho * (r_top[2] - x[0]);
    // f[2]=r_top[2];
  }


  //=======================================================================
  /// Boundary of central box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_centr_D(const unsigned& t,
                                    const Vector<double>& zeta,
                                    const unsigned& i_layer,
                                    Vector<double>& f)
  {
    // Wall coordinates along bottom edge of wall
    Vector<double> x(2);
    x[0] =
      Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                   axial_spacing_fct((0.5 * (1.0 + zeta[1]) + double(i_layer)) /
                                     double(Nlayer));
    x[1] = Xi_lo[1];

    // Get position vector to bottom edge of wall
    Vector<double> r_bottom(3);
    Wall_pt->position(t, x, r_bottom);

    // Scale it down to half the width
    f[0] = r_bottom[0] * 0.25 * (1.0 + zeta[0]);
    f[1] = r_bottom[1];

    // Warp it:
    double rho = 0.0; // 0.25*(1.0+zeta[0]);
    f[2] = x[0] + rho * (r_bottom[2] - x[0]);
    // f[2]=r_bottom[2];
  }


  //=======================================================================
  /// Boundary of central box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_centr_U(const unsigned& t,
                                    const Vector<double>& zeta,
                                    const unsigned& i_layer,
                                    Vector<double>& f)
  {
    // Wall coordinates along top edge of wall
    Vector<double> x(2);
    x[0] =
      Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                   axial_spacing_fct((0.5 * (1.0 + zeta[1]) + double(i_layer)) /
                                     double(Nlayer));
    x[1] = Xi_hi[1];

    // Get position vector to upper edge of wall
    Vector<double> r_top(3);
    Wall_pt->position(t, x, r_top);

    // Wall coordinates along bottom edge of wall
    x[0] =
      Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                   axial_spacing_fct((0.5 * (1.0 + zeta[1]) + double(i_layer)) /
                                     double(Nlayer));
    x[1] = Xi_lo[1];

    // Get position vector to bottom edge of wall
    Vector<double> r_bottom(3);
    Wall_pt->position(t, x, r_bottom);

    // Scale it down to half the width
    f[0] = r_bottom[0] * 0.25 * (1.0 + zeta[0]);
    f[1] = 0.5 * r_top[1];

    // Warp it:
    double rho = 0.0; // 0.25*(1.0+zeta[0]);
    f[2] = x[0] + rho * (r_bottom[2] - x[0]);
    // f[2]=r_bottom[2];
  }


  //=======================================================================
  /// Boundary of central box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_centr_B(const unsigned& t,
                                    const Vector<double>& zeta,
                                    const unsigned& i_layer,
                                    Vector<double>& f)
  {
    // Wall coordinates along bottom edge of wall
    Vector<double> x(2);
    x[0] = Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                        axial_spacing_fct(double(i_layer) / double(Nlayer));
    x[1] = Xi_lo[1];

    // Get position vector to bottom edge of wall
    Vector<double> r_bottom(3);
    Wall_pt->position(t, x, r_bottom);


    // Wall coordinates along top edge of wall
    x[0] = Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                        axial_spacing_fct(double(i_layer) / double(Nlayer));
    x[1] = Xi_hi[1];

    // Get position vector to top edge of wall
    Vector<double> r_top(3);
    Wall_pt->position(t, x, r_top);

    // Map it
    f[0] = r_bottom[0] * 0.25 * (1.0 + zeta[0]);
    f[1] = r_top[1] * 0.25 * (1.0 + zeta[1]);

    // Warp it:
    double rho = 0.0; // 0.25*(1.0+zeta[1]);
    f[2] = x[0] + rho * (r_top[2] - x[0]);
    // f[2]=r_top[2];
  }


  //=======================================================================
  /// Boundary of central box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_centr_F(const unsigned& t,
                                    const Vector<double>& zeta,
                                    const unsigned& i_layer,
                                    Vector<double>& f)
  {
    // Wall coordinates along bottom edge of wall
    Vector<double> x(2);
    x[0] = Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                        axial_spacing_fct(double(1 + i_layer) / double(Nlayer));
    x[1] = Xi_lo[1];

    // Get position vector to bottom edge of wall
    Vector<double> r_bottom(3);
    Wall_pt->position(t, x, r_bottom);


    // Wall coordinates along top edge of wall
    x[0] = Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                        axial_spacing_fct(double(1 + i_layer) / double(Nlayer));
    x[1] = Xi_hi[1];

    // Get position vector to top edge of wall
    Vector<double> r_top(3);
    Wall_pt->position(t, x, r_top);

    // Map it
    f[0] = r_bottom[0] * 0.25 * (1.0 + zeta[0]);
    f[1] = r_top[1] * 0.25 * (1.0 + zeta[1]);

    // Warp it:
    double rho = 0.0; // 0.25*(1.0+zeta[1]);
    f[2] = x[0] + rho * (r_top[2] - x[0]);
    // f[2]=r_top[2];
  }


  //#####################################################################


  //=======================================================================
  /// Boundary of bottom right box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_bot_right_L(const unsigned& t,
                                        const Vector<double>& zeta,
                                        const unsigned& i_layer,
                                        Vector<double>& f)
  {
    r_centr_R(t, zeta, i_layer, f);
  }


  //=======================================================================
  /// Boundary of bottom right box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_bot_right_R(const unsigned& t,
                                        const Vector<double>& zeta,
                                        const unsigned& i_layer,
                                        Vector<double>& f)
  {
    // Wall coordinates
    Vector<double> x(2);
    x[0] =
      Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                   axial_spacing_fct((0.5 * (1.0 + zeta[1]) + double(i_layer)) /
                                     double(Nlayer));
    x[1] = Xi_lo[1] + Fract_mid * (Xi_hi[1] - Xi_lo[1]) * 0.5 * (1.0 + zeta[0]);

    // Get position vector on wall
    Wall_pt->position(t, x, f);
  }


  //=======================================================================
  /// Boundary of bottom right box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_bot_right_D(const unsigned& t,
                                        const Vector<double>& zeta,
                                        const unsigned& i_layer,
                                        Vector<double>& f)
  {
    // Wall coordinates along bottom edge of wall
    Vector<double> x(2);
    x[0] =
      Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                   axial_spacing_fct((0.5 * (1.0 + zeta[1]) + double(i_layer)) /
                                     double(Nlayer));
    x[1] = Xi_lo[1];

    // Get position vector to bottom edge of wall
    Vector<double> r_bottom(3);
    Wall_pt->position(t, x, r_bottom);

    // Scale it down to half the width
    f[0] = 0.5 * r_bottom[0] * (1.0 + s_squashed(0.5 * (1.0 + zeta[0])));
    f[1] = r_bottom[1];

    // Warp it:
    double rho = s_squashed(0.5 * (1.0 + zeta[0]));
    f[2] = x[0] + rho * (r_bottom[2] - x[0]);
    // f[2]=r_bottom[2];
  }


  //=======================================================================
  /// Boundary of bottom right box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_bot_right_U(const unsigned& t,
                                        const Vector<double>& zeta,
                                        const unsigned& i_layer,
                                        Vector<double>& f)
  {
    // Wall coordinates of dividing line
    Vector<double> x(2);
    x[0] =
      Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                   axial_spacing_fct((0.5 * (1.0 + zeta[1]) + double(i_layer)) /
                                     double(Nlayer));
    x[1] = Xi_lo[1] + Fract_mid * (Xi_hi[1] - Xi_lo[1]);

    // Get position vector on dividing line
    Vector<double> r_div(3);
    Wall_pt->position(t, x, r_div);


    // Position vector to corner of central box
    Vector<double> zeta_central(2);
    Vector<double> r_central(3);
    zeta_central[0] = 1.0;
    zeta_central[1] = zeta[1];
    r_centr_R(t, zeta_central, i_layer, r_central);


    // Straight line across
    f[0] = r_central[0] +
           (r_div[0] - r_central[0]) * s_squashed(0.5 * (1.0 + zeta[0]));
    f[1] = r_central[1] +
           (r_div[1] - r_central[1]) * s_squashed(0.5 * (1.0 + zeta[0]));
    f[2] = r_central[2] +
           (r_div[2] - r_central[2]) * s_squashed(0.5 * (1.0 + zeta[0]));
  }


  //=======================================================================
  /// Boundary of bottom right box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_bot_right_B(const unsigned& t,
                                        const Vector<double>& zeta,
                                        const unsigned& i_layer,
                                        Vector<double>& f)
  {
    // Wall coordinates
    Vector<double> x(2);
    x[0] = Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                        axial_spacing_fct(double(i_layer) / double(Nlayer));
    x[1] = Xi_lo[1] + Fract_mid * (Xi_hi[1] - Xi_lo[1]) * 0.5 * (1.0 + zeta[1]);

    // Get position vector to wall
    Vector<double> r_wall(3);
    Wall_pt->position(t, x, r_wall);

    // Get position vector on central box
    Vector<double> zeta_central(2);
    Vector<double> r_central(3);
    zeta_central[0] = zeta[1];
    zeta_central[1] = -1.0;
    r_centr_R(t, zeta_central, i_layer, r_central);


    // Straight line across
    f[0] = r_central[0] +
           (r_wall[0] - r_central[0]) * s_squashed(0.5 * (1.0 + zeta[0]));
    f[1] = r_central[1] +
           (r_wall[1] - r_central[1]) * s_squashed(0.5 * (1.0 + zeta[0]));
    f[2] = r_central[2] +
           (r_wall[2] - r_central[2]) * s_squashed(0.5 * (1.0 + zeta[0]));
  }


  //=======================================================================
  /// Boundary of bottom right box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_bot_right_F(const unsigned& t,
                                        const Vector<double>& zeta,
                                        const unsigned& i_layer,
                                        Vector<double>& f)
  {
    // Wall coordinates
    Vector<double> x(2);
    x[0] = Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                        axial_spacing_fct(double(1 + i_layer) / double(Nlayer));
    x[1] = Xi_lo[1] + Fract_mid * (Xi_hi[1] - Xi_lo[1]) * 0.5 * (1.0 + zeta[1]);

    // Get position vector to wall
    Vector<double> r_wall(3);
    Wall_pt->position(t, x, r_wall);

    // Get position vector on central box
    Vector<double> zeta_central(2);
    Vector<double> r_central(3);
    zeta_central[0] = zeta[1];
    zeta_central[1] = 1.0;
    r_centr_R(t, zeta_central, i_layer, r_central);


    // Straight line across
    f[0] = r_central[0] +
           (r_wall[0] - r_central[0]) * s_squashed(0.5 * (1.0 + zeta[0]));
    f[1] = r_central[1] +
           (r_wall[1] - r_central[1]) * s_squashed(0.5 * (1.0 + zeta[0]));
    f[2] = r_central[2] +
           (r_wall[2] - r_central[2]) * s_squashed(0.5 * (1.0 + zeta[0]));
  }


  //#####################################################################


  //=======================================================================
  /// Boundary of top left box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_top_left_L(const unsigned& t,
                                       const Vector<double>& zeta,
                                       const unsigned& i_layer,
                                       Vector<double>& f)
  {
    // Wall coordinates along top edge of wall
    Vector<double> x(2);
    x[0] =
      Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                   axial_spacing_fct((0.5 * (1.0 + zeta[1]) + double(i_layer)) /
                                     double(Nlayer));
    x[1] = Xi_hi[1];

    // Get position vector to upper edge of wall
    Vector<double> r_top(3);
    Wall_pt->position(t, x, r_top);

    // Scale it down to half the height
    f[0] = r_top[0];
    f[1] = 0.5 * r_top[1] * (1.0 + s_squashed(0.5 * (1.0 + zeta[0])));

    // Warp it:
    double rho = s_squashed(0.5 * (1.0 + zeta[0]));
    f[2] = x[0] + rho * (r_top[2] - x[0]);
    // f[2]=r_top[2];
  }


  //=======================================================================
  /// Boundary of top left box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_top_left_R(const unsigned& t,
                                       const Vector<double>& zeta,
                                       const unsigned& i_layer,
                                       Vector<double>& f)
  {
    // Swap coordinates
    Vector<double> zeta_br(2);
    zeta_br[0] = zeta[0];
    zeta_br[1] = zeta[1];
    r_bot_right_U(t, zeta_br, i_layer, f);
  }


  //=======================================================================
  /// Boundary of top left box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_top_left_D(const unsigned& t,
                                       const Vector<double>& zeta,
                                       const unsigned& i_layer,
                                       Vector<double>& f)
  {
    r_centr_U(t, zeta, i_layer, f);
  }


  //=======================================================================
  /// Boundary of top left box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_top_left_U(const unsigned& t,
                                       const Vector<double>& zeta,
                                       const unsigned& i_layer,
                                       Vector<double>& f)
  {
    // Wall coordinates
    Vector<double> x(2);
    x[0] =
      Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                   axial_spacing_fct((0.5 * (1.0 + zeta[1]) + double(i_layer)) /
                                     double(Nlayer));
    x[1] = Xi_hi[1] +
           (Xi_lo[1] - Xi_hi[1]) * (1 - Fract_mid) * 0.5 * (1.0 + zeta[0]);

    // Get position vector on wall
    Wall_pt->position(t, x, f);
  }


  //=======================================================================
  /// Boundary of top left box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_top_left_B(const unsigned& t,
                                       const Vector<double>& zeta,
                                       const unsigned& i_layer,
                                       Vector<double>& f)
  {
    // Wall coordinates
    Vector<double> x(2);
    x[0] = Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                        axial_spacing_fct(double(i_layer) / double(Nlayer));
    x[1] = Xi_hi[1] +
           (Xi_lo[1] - Xi_hi[1]) * (1.0 - Fract_mid) * 0.5 * (1.0 + zeta[0]);


    // Get position vector to wall
    Vector<double> r_wall(3);
    Wall_pt->position(t, x, r_wall);


    // Get position vector on central box
    Vector<double> zeta_central(2);
    Vector<double> r_central(3);
    zeta_central[0] = zeta[0];
    zeta_central[1] = -1.0;
    r_centr_U(t, zeta_central, i_layer, r_central);

    // Straight line across
    f[0] = r_central[0] +
           (r_wall[0] - r_central[0]) * s_squashed(0.5 * (1.0 + zeta[1]));
    f[1] = r_central[1] +
           (r_wall[1] - r_central[1]) * s_squashed(0.5 * (1.0 + zeta[1]));
    f[2] = r_central[2] +
           (r_wall[2] - r_central[2]) * s_squashed(0.5 * (1.0 + zeta[1]));
  }


  //=======================================================================
  /// Boundary of top left box macro element in layer i_layer
  /// zeta \f$ \in [-1,1]^2 \f$
  //=======================================================================
  void QuarterTubeDomain::r_top_left_F(const unsigned& t,
                                       const Vector<double>& zeta,
                                       const unsigned& i_layer,
                                       Vector<double>& f)
  {
    // Wall coordinates
    Vector<double> x(2);
    x[0] = Xi_lo[0] + (Xi_hi[0] - Xi_lo[0]) *
                        axial_spacing_fct(double(1 + i_layer) / double(Nlayer));
    x[1] = Xi_hi[1] +
           (Xi_lo[1] - Xi_hi[1]) * (1.0 - Fract_mid) * 0.5 * (1.0 + zeta[0]);


    // Get position vector to wall
    Vector<double> r_wall(3);
    Wall_pt->position(t, x, r_wall);


    // Get position vector on central box
    Vector<double> zeta_central(2);
    Vector<double> r_central(3);
    zeta_central[0] = zeta[0];
    zeta_central[1] = 1.0;
    r_centr_U(t, zeta_central, i_layer, r_central);

    // Straight line across
    f[0] = r_central[0] +
           (r_wall[0] - r_central[0]) * s_squashed(0.5 * (1.0 + zeta[1]));
    f[1] = r_central[1] +
           (r_wall[1] - r_central[1]) * s_squashed(0.5 * (1.0 + zeta[1]));
    f[2] = r_central[2] +
           (r_wall[2] - r_central[2]) * s_squashed(0.5 * (1.0 + zeta[1]));
  }


} // namespace oomph

#endif
