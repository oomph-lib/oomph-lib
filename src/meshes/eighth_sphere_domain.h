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
#ifndef OOMPH_EIGHTH_SPHERE_DOMAIN_HEADER
#define OOMPH_EIGHTH_SPHERE_DOMAIN_HEADER

// Generic oomph-lib includes
#include "../generic/quadtree.h"
#include "../generic/domain.h"
#include "../generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  ///  Eighth sphere as domain. Domain is
  /// parametrised by four macro elements
  //=================================================================
  class EighthSphereDomain : public Domain
  {
  public:
    ///  Constructor: Pass the radius of the sphere.
    EighthSphereDomain(const double& radius) : Radius(radius)
    {
      // There are four macro elements
      unsigned nmacro = 4;

      // Resize
      Macro_element_pt.resize(nmacro);

      // Create the 3D Q macro elements
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<3>(this, i);
      }
    }

    /// Broken copy constructor
    EighthSphereDomain(const EighthSphereDomain&) = delete;

    /// Broken assignment operator
    void operator=(const EighthSphereDomain&) = delete;

    /// Destructor: Empty; cleanup done in base class
    ~EighthSphereDomain() {}


    ///  Vector representation of the  imacro-th macro element
    /// boundary idirect (L/R/D/U/B/F) at time level t
    /// (t=0: present; t>0: previous):
    /// f(s).
    void macro_element_boundary(const unsigned& t,
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
        "EighthSphereDomain::macro_element_boundary(...)",
        OOMPH_EXCEPTION_LOCATION);
#endif

      // Which macro element?
      // --------------------
      switch (imacro)
      {
          // Macro element 0: Central box
        case 0:

          if (idirect == L)
          {
            r_centr_L(t, s, f);
          }
          else if (idirect == R)
          {
            r_centr_R(t, s, f);
          }
          else if (idirect == D)
          {
            r_centr_D(t, s, f);
          }
          else if (idirect == U)
          {
            r_centr_U(t, s, f);
          }
          else if (idirect == B)
          {
            r_centr_B(t, s, f);
          }
          else if (idirect == F)
          {
            r_centr_F(t, s, f);
          }
          else
          {
            std::ostringstream error_message;
            error_message << "idirect is " << OcTree::Direct_string[idirect]
                          << "not one of L, R, U, D, B, F" << std::endl;

            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          break;


          // Macro element 1:right
        case 1:

          // Which direction?
          if (idirect == L)
          {
            r_right_L(t, s, f);
          }
          else if (idirect == R)
          {
            r_right_R(t, s, f);
          }
          else if (idirect == D)
          {
            r_right_D(t, s, f);
          }
          else if (idirect == U)
          {
            r_right_U(t, s, f);
          }
          else if (idirect == B)
          {
            r_right_B(t, s, f);
          }
          else if (idirect == F)
          {
            r_right_F(t, s, f);
          }
          else
          {
            std::ostringstream error_message;
            error_message << "idirect is " << OcTree::Direct_string[idirect]
                          << "not one of L, R, U, D, B, F" << std::endl;

            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          break;

          // Macro element 2: Up
        case 2:

          // Which direction?
          if (idirect == L)
          {
            r_up_L(t, s, f);
          }
          else if (idirect == R)
          {
            r_up_R(t, s, f);
          }
          else if (idirect == D)
          {
            r_up_D(t, s, f);
          }
          else if (idirect == U)
          {
            r_up_U(t, s, f);
          }
          else if (idirect == B)
          {
            r_up_B(t, s, f);
          }
          else if (idirect == F)
          {
            r_up_F(t, s, f);
          }
          else
          {
            std::ostringstream error_message;
            error_message << "idirect is " << OcTree::Direct_string[idirect]
                          << "not one of L, R, U, D, B, F" << std::endl;

            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          break;

          // Macro element 3: Front
        case 3:
          // Which direction?
          if (idirect == L)
          {
            r_front_L(t, s, f);
          }
          else if (idirect == R)
          {
            r_front_R(t, s, f);
          }
          else if (idirect == D)
          {
            r_front_D(t, s, f);
          }
          else if (idirect == U)
          {
            r_front_U(t, s, f);
          }
          else if (idirect == B)
          {
            r_front_B(t, s, f);
          }
          else if (idirect == F)
          {
            r_front_F(t, s, f);
          }
          else
          {
            std::ostringstream error_message;
            error_message << "idirect is " << OcTree::Direct_string[idirect]
                          << "not one of L, R, U, D, B, F" << std::endl;

            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          break;

        default:

          std::ostringstream error_message;
          error_message << "imacro is " << OcTree::Direct_string[idirect]
                        << ", but should be in the range 0-3" << std::endl;

          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }


  private:
    // Radius of the sphere
    double Radius;

    ///  Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_L(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_R(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_D(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_U(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_B(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_F(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_L(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_R(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of  right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_D(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_U(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of  right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_B(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of  right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_F(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_L(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_R(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_D(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_U(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_B(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_F(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_L(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_R(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_D(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_U(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_B(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    ///  Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_F(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);
  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Boundary of central box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_centr_L(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    f[0] = 0;
    f[1] = Radius * 0.25 * (1.0 + zeta[0]);
    f[2] = Radius * 0.25 * (1.0 + zeta[1]);
  }


  //=================================================================
  ///  Boundary of central box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_centr_R(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    f[0] = Radius * 0.5;
    f[1] = Radius * 0.25 * (1.0 + zeta[0]);
    f[2] = Radius * 0.25 * (1.0 + zeta[1]);
  }


  //=================================================================
  ///  Boundary of central box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_centr_D(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    f[0] = Radius * 0.25 * (1.0 + zeta[0]);
    f[1] = 0;
    f[2] = Radius * 0.25 * (1.0 + zeta[1]);
  }


  //=================================================================
  ///  Boundary of central box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_centr_U(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    f[0] = Radius * 0.25 * (1.0 + zeta[0]);
    f[1] = Radius * 0.5;
    f[2] = Radius * 0.25 * (1.0 + zeta[1]);
  }


  //=================================================================
  ///  Boundary of central box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_centr_B(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    f[0] = Radius * 0.25 * (1.0 + zeta[0]);
    f[1] = Radius * 0.25 * (1.0 + zeta[1]);
    f[2] = 0.0;
  }

  //=================================================================
  ///  Boundary of central box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_centr_F(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    f[0] = Radius * 0.25 * (1.0 + zeta[0]);
    f[1] = Radius * 0.25 * (1.0 + zeta[1]);
    f[2] = Radius * 0.5;
  }


  //=================================================================
  ///  Boundary of right box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_right_L(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    r_centr_R(t, zeta, f);
  }


  //=================================================================
  ///  Boundary of right box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_right_R(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    double k0 = 0.5 * (1.0 + zeta[0]);
    double k1 = 0.5 * (1.0 + zeta[1]);
    Vector<double> p(3);
    Vector<double> point1(3);
    Vector<double> point2(3);

    point1[0] = Radius - 0.5 * Radius * k0;
    point1[1] = 0.5 * Radius * k0;
    point1[2] = 0.0;
    point2[0] = 0.5 * Radius - k0 * Radius / 6.0;
    point2[1] = k0 * Radius / 3.0;
    point2[2] = 0.5 * Radius - k0 * Radius / 6.0;

    for (unsigned i = 0; i < 3; i++)
    {
      p[i] = point1[i] + k1 * (point2[i] - point1[i]);
    }
    double alpha = Radius / std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);

    for (unsigned i = 0; i < 3; i++)
    {
      f[i] = alpha * p[i];
    }
  }

  //=================================================================
  ///  Boundary of  right box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_right_D(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    // position vector on sphere
    Vector<double> on_sphere(3);
    Vector<double> temp_zeta(2);
    temp_zeta[0] = -1.0;
    temp_zeta[1] = zeta[1];
    r_right_R(t, temp_zeta, on_sphere);

    // position vector on center box
    Vector<double> on_center(3);
    on_center[0] = 0.5 * Radius;
    on_center[1] = 0.0;
    on_center[2] = 0.5 * Radius * 0.5 * (zeta[1] + 1.0);
    // strait line across
    for (unsigned i = 0; i < 3; i++)
    {
      f[i] =
        on_center[i] + 0.5 * (zeta[0] + 1.0) * (on_sphere[i] - on_center[i]);
    }
  }


  //=================================================================
  ///  Boundary of right box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_right_U(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    // position vector on sphere
    Vector<double> on_sphere(3);
    Vector<double> temp_zeta(2);
    temp_zeta[0] = 1.0;
    temp_zeta[1] = zeta[1];
    r_right_R(t, temp_zeta, on_sphere);

    // position vector on center box
    Vector<double> on_center(3);
    on_center[0] = 0.5 * Radius;
    on_center[1] = 0.5 * Radius;
    on_center[2] = 0.5 * Radius * 0.5 * (zeta[1] + 1.0);
    // strait line across
    for (unsigned i = 0; i < 3; i++)
    {
      f[i] =
        on_center[i] + 0.5 * (zeta[0] + 1.0) * (on_sphere[i] - on_center[i]);
    }
  }

  //=================================================================
  ///  Boundary of  right box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_right_B(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    // position vector on sphere
    Vector<double> on_sphere(3);
    Vector<double> temp_zeta(2);
    temp_zeta[0] = zeta[1];
    temp_zeta[1] = -1.0;
    r_right_R(t, temp_zeta, on_sphere);

    // position vector on center box
    Vector<double> on_center(3);
    on_center[0] = 0.5 * Radius;
    on_center[1] = 0.5 * Radius * 0.5 * (zeta[1] + 1.0);
    on_center[2] = 0.0;
    // strait line across
    for (unsigned i = 0; i < 3; i++)
    {
      f[i] =
        on_center[i] + 0.5 * (zeta[0] + 1.0) * (on_sphere[i] - on_center[i]);
    }
  }


  //=================================================================
  ///  Boundary of  right box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_right_F(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    // position vector on sphere
    Vector<double> on_sphere(3);
    Vector<double> temp_zeta(2);
    temp_zeta[0] = zeta[1];
    temp_zeta[1] = 1.0;
    r_right_R(t, temp_zeta, on_sphere);

    // position vector on center box
    Vector<double> on_center(3);
    on_center[0] = 0.5 * Radius;
    on_center[1] = 0.5 * Radius * 0.5 * (zeta[1] + 1.0);
    on_center[2] = 0.5 * Radius;
    // strait line across
    for (unsigned i = 0; i < 3; i++)
    {
      f[i] =
        on_center[i] + 0.5 * (zeta[0] + 1.0) * (on_sphere[i] - on_center[i]);
    }
  }

  //=================================================================
  ///  Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_up_L(const unsigned& t,
                                  const Vector<double>& zeta,
                                  Vector<double>& f)
  {
    // position vector on sphere
    Vector<double> on_sphere(3);
    Vector<double> temp_zeta(2);
    temp_zeta[0] = -1.0;
    temp_zeta[1] = zeta[1];
    r_up_U(t, temp_zeta, on_sphere);

    // position vector on center box
    Vector<double> on_center(3);
    on_center[0] = 0.0;
    on_center[1] = 0.5 * Radius;
    on_center[2] = 0.5 * Radius * 0.5 * (zeta[1] + 1.0);
    // strait line across
    for (unsigned i = 0; i < 3; i++)
    {
      f[i] =
        on_center[i] + 0.5 * (zeta[0] + 1.0) * (on_sphere[i] - on_center[i]);
    }
  }


  //=================================================================
  ///  Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_up_R(const unsigned& t,
                                  const Vector<double>& zeta,
                                  Vector<double>& f)
  {
    r_right_U(t, zeta, f);
  }


  //=================================================================
  ///  Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_up_D(const unsigned& t,
                                  const Vector<double>& zeta,
                                  Vector<double>& f)
  {
    r_centr_U(t, zeta, f);
  }

  //=================================================================
  ///  Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_up_U(const unsigned& t,
                                  const Vector<double>& zeta,
                                  Vector<double>& f)
  {
    double k0 = 0.5 * (1.0 + zeta[0]);
    double k1 = 0.5 * (1.0 + zeta[1]);
    Vector<double> p(3);
    Vector<double> point1(3);
    Vector<double> point2(3);

    point1[0] = 0.5 * Radius * k0;
    point1[1] = Radius - 0.5 * Radius * k0;
    point1[2] = 0;
    point2[0] = k0 * Radius / 3.0;
    point2[1] = 0.5 * Radius - k0 * Radius / 6.0;
    point2[2] = 0.5 * Radius - k0 * Radius / 6.0;

    for (unsigned i = 0; i < 3; i++)
    {
      p[i] = point1[i] + k1 * (point2[i] - point1[i]);
    }
    double alpha = Radius / std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);

    for (unsigned i = 0; i < 3; i++)
    {
      f[i] = alpha * p[i];
    }
  }

  //=================================================================
  ///  Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_up_B(const unsigned& t,
                                  const Vector<double>& zeta,
                                  Vector<double>& f)
  {
    // position vector on sphere
    Vector<double> on_sphere(3);
    Vector<double> temp_zeta(2);
    temp_zeta[0] = zeta[0];
    temp_zeta[1] = -1.0;
    r_up_U(t, temp_zeta, on_sphere);

    // position vector on center box
    Vector<double> on_center(3);
    on_center[0] = 0.5 * Radius * 0.5 * (zeta[0] + 1.0);
    on_center[1] = 0.5 * Radius;
    on_center[2] = 0.0;
    // strait line across
    for (unsigned i = 0; i < 3; i++)
    {
      f[i] =
        on_center[i] + 0.5 * (zeta[1] + 1.0) * (on_sphere[i] - on_center[i]);
    }
  }


  //=================================================================
  ///  Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_up_F(const unsigned& t,
                                  const Vector<double>& zeta,
                                  Vector<double>& f)
  {
    // position vector on sphere
    Vector<double> on_sphere(3);
    Vector<double> temp_zeta(2);
    temp_zeta[0] = zeta[0];
    temp_zeta[1] = 1.0;
    r_up_U(t, temp_zeta, on_sphere);

    // position vector on center box
    Vector<double> on_center(3);
    on_center[0] = 0.5 * Radius * 0.5 * (zeta[0] + 1.0);
    on_center[1] = 0.5 * Radius;
    on_center[2] = 0.5 * Radius;
    // straight line across
    for (unsigned i = 0; i < 3; i++)
    {
      f[i] =
        on_center[i] + 0.5 * (zeta[1] + 1.0) * (on_sphere[i] - on_center[i]);
    }
  }


  //=================================================================
  ///  Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_front_L(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    // position vector on sphere
    Vector<double> on_sphere(3);
    Vector<double> temp_zeta(2);
    temp_zeta[0] = -1.0;
    temp_zeta[1] = zeta[0];
    r_front_F(t, temp_zeta, on_sphere);

    // position vector on center box
    Vector<double> on_center(3);
    on_center[0] = 0.0;
    on_center[1] = 0.5 * Radius * 0.5 * (zeta[0] + 1.0);
    on_center[2] = 0.5 * Radius;
    // straight line across
    for (unsigned i = 0; i < 3; i++)
    {
      f[i] =
        on_center[i] + 0.5 * (zeta[1] + 1.0) * (on_sphere[i] - on_center[i]);
    }
  }


  //=================================================================
  ///  Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_front_R(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    Vector<double> zeta2(2);
    zeta2[0] = zeta[1];
    zeta2[1] = zeta[0];
    r_right_F(t, zeta2, f);
  }


  //=================================================================
  ///  Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_front_D(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    // position vector on sphere
    Vector<double> on_sphere(3);
    Vector<double> temp_zeta(2);
    temp_zeta[0] = zeta[0];
    temp_zeta[1] = -1.0;
    r_front_F(t, temp_zeta, on_sphere);

    // position vector on center box
    Vector<double> on_center(3);
    on_center[0] = 0.5 * Radius * 0.5 * (zeta[0] + 1.0);
    on_center[1] = 0.0;
    on_center[2] = 0.5 * Radius;
    // straight line across
    for (unsigned i = 0; i < 3; i++)
    {
      f[i] =
        on_center[i] + 0.5 * (zeta[1] + 1.0) * (on_sphere[i] - on_center[i]);
    }
  }


  //=================================================================
  ///  Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_front_U(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    r_up_F(t, zeta, f);
  }


  //=================================================================
  ///  Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_front_B(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    r_centr_F(t, zeta, f);
  }


  //=================================================================
  /// Boundary of top left box macro element
  /// zeta \f$ \in [-1,1]^2 \f$
  //=================================================================
  void EighthSphereDomain::r_front_F(const unsigned& t,
                                     const Vector<double>& zeta,
                                     Vector<double>& f)
  {
    double k0 = 0.5 * (1.0 + zeta[0]);
    double k1 = 0.5 * (1.0 + zeta[1]);
    Vector<double> p(3);
    Vector<double> point1(3);
    Vector<double> point2(3);

    point1[0] = 0.5 * Radius * k0;
    point1[1] = 0.0;
    point1[2] = Radius - k0 * 0.5 * Radius;
    point2[0] = k0 * Radius / 3.0;
    point2[1] = 0.5 * Radius - k0 * Radius / 6.0;
    point2[2] = 0.5 * Radius - k0 * Radius / 6.0;

    for (unsigned i = 0; i < 3; i++)
    {
      p[i] = point1[i] + k1 * (point2[i] - point1[i]);
    }
    double alpha = Radius / std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);

    for (unsigned i = 0; i < 3; i++)
    {
      f[i] = alpha * p[i];
    }
  }


} // namespace oomph

#endif
