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

#ifndef OOMPH_FULL_CIRCLE_DOMAIN_HEADER
#define OOMPH_FULL_CIRCLE_DOMAIN_HEADER

// Generic oomph-lib includes
#include "../generic/quadtree.h"
#include "../generic/domain.h"
#include "../generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  ///  Topologically circular domain, e.g. a tube cross section.
  /// The entire domain must be defined by a GeomObject with the
  /// following convention: zeta[0] is the radial coordinate and
  /// zeta[1] is the theta coordinate around the cross-sectin.
  /// The outer boundary must lie at zeta[0] = 1.
  ///
  /// The domain is
  /// parametrised by five macro elements (a central box surrounded by
  /// four curved elements). The labelling of the macro elements is shown
  /// below.
  ///
  ///                       ----------------------------
  ///                       |\                        /|
  ///                       | \       Macro          / |
  ///                       |  3      Element 3     2  |
  ///                       |   \                  /   |
  ///                       |    \----------------/    |
  ///                       |     |               |    |
  ///                       | 4   |    Macro      |    |
  ///                       |     |    Element 0  | 2  |
  ///                       |     |               |    |
  ///                       |     -----------------    |
  ///                       |    /                 \   |
  ///                       |   0     Macro         1  |
  ///                       |  /      Element 1      \ |
  ///                       | /                       \|
  ///                       |/-------------------------|
  ///
  ///
  //=================================================================
  class FullCircleDomain : public Domain
  {
  public:
    ///  Constructor: Pass geometric  object; the theta locations
    /// marking the division between
    /// the elements of the outer ring, labelled from the lower left to the
    /// upper left in order, theta should be in the range \f$-\pi\f$ to
    /// \f$\pi\f$; and  the corresponding fractions of the
    /// radius at which the central box is to be placed.
    FullCircleDomain(GeomObject* area_geom_object_pt,
                     const Vector<double>& theta_positions,
                     const Vector<double>& radius_box)
      : Theta_positions(theta_positions),
        Radius_box(radius_box),
        Area_pt(area_geom_object_pt)
    {
      // There are five macro elements
      const unsigned n_macro = 5;
      Macro_element_pt.resize(n_macro);

      // Create the macro elements
      for (unsigned i = 0; i < n_macro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    }

    /// Broken copy constructor
    FullCircleDomain(const FullCircleDomain&) = delete;

    /// Broken assignment operator
    void operator=(const FullCircleDomain&) = delete;


    /// Destructor: Empty; cleanup done in base class
    ~FullCircleDomain() {}


    ///  Vector representation of the  i_macro-th macro element
    /// boundary i_direct (N/S/W/E) at time level t
    /// (t=0: present; t>0: previous):
    /// f(s).
    void macro_element_boundary(const unsigned& t,
                                const unsigned& i_macro,
                                const unsigned& i_direct,
                                const Vector<double>& s,
                                Vector<double>& f);

  private:
    ///  Storage for the dividing lines on the boundary
    /// starting from the lower left and proceeding anticlockwise to
    /// the upper left
    Vector<double> Theta_positions;

    /// Storage for the fraction of the radius at which the central box
    /// should be located corresponding to the chosen values of theta.
    Vector<double> Radius_box;

    /// Pointer to geometric object that represents the domain
    GeomObject* Area_pt;

    ///   A very little linear interpolation helper.
    /// Interpolate from the low
    /// point to the high point using the coordinate s, which is
    /// assumed to run from -1 to 1.
    void lin_interpolate(const Vector<double>& low,
                         const Vector<double>& high,
                         const double& s,
                         Vector<double>& f)
    {
      // Loop over all coordinates (we are working in 2D)
      for (unsigned i = 0; i < 2; i++)
      {
        f[i] = low[i] + (high[i] - low[i]) * 0.5 * (s + 1.0);
      }
    }
  };


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


  //=================================================================
  ///  Vector representation of the  imacro-th macro element
  /// boundary idirect (N/S/W/E) at time level t
  /// (t=0: present; t>0: previous): f(s)
  //=================================================================
  void FullCircleDomain::macro_element_boundary(const unsigned& t,
                                                const unsigned& imacro,
                                                const unsigned& idirect,
                                                const Vector<double>& s,
                                                Vector<double>& f)
  {
    using namespace QuadTreeNames;

    // Get the coordinates of the corners of the box
    Vector<Vector<double>> Box(4);
    // Get the corresponding coordinates on the boundary
    Vector<Vector<double>> Wall(4);

    // Storage for position in the area
    Vector<double> zeta(2);

    // Loop over the angles
    for (unsigned j = 0; j < 4; j++)
    {
      // Set up the storage
      Box[j].resize(2);
      Wall[j].resize(2);

      // Set the other values of zeta
      zeta[0] = Radius_box[j];
      zeta[1] = Theta_positions[j];
      // Now get the values
      Area_pt->position(t, zeta, Box[j]);

      // Now get the position on the boundary
      zeta[0] = 1.0;
      Area_pt->position(t, zeta, Wall[j]);
    }

    // Define pi
    const double pi = MathematicalConstants::Pi;

    // Which macro element?
    // --------------------
    switch (imacro)
    {
        // Macro element 0: Central box
      case 0:

        // Choose a direction
        switch (idirect)
        {
          case N:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[3], Box[2], s[0], f);
            break;

          case S:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[0], Box[1], s[0], f);

          case W:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[0], Box[3], s[0], f);
            break;

          case E:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[1], Box[2], s[0], f);
            break;

          default:
            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                         << std::endl;

            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
            break;
        }

        break;

        // Macro element 1: Bottom
      case 1:

        // Choose a direction
        switch (idirect)
        {
          case N:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[0], Box[1], s[0], f);
            break;

          case S:
            // Get the position on the wall by linearly interpolating in theta
            zeta[0] = 1.0;
            zeta[1] =
              Theta_positions[0] +
              (Theta_positions[1] - Theta_positions[0]) * 0.5 * (s[0] + 1.0);

            Area_pt->position(t, zeta, f);
            break;

          case W:
            // Now linearly interpolate between the wall and the box
            lin_interpolate(Wall[0], Box[0], s[0], f);
            break;

          case E:
            // Linearly interpolate between the wall and the box
            lin_interpolate(Wall[1], Box[1], s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                         << std::endl;

            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
            break;
        }


        break;

      // Macro element 2:Right
      case 2:

        // Which direction?
        switch (idirect)
        {
          case N:
            // Linearly interpolate between the box and the wall
            lin_interpolate(Box[2], Wall[2], s[0], f);
            break;

          case S:
            // Linearly interpolate between the box and the wall
            lin_interpolate(Box[1], Wall[1], s[0], f);
            break;

          case W:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[1], Box[2], s[0], f);
            break;

          case E:
            // Get the position on the wall by linearly interpolating in theta
            zeta[0] = 1.0;
            zeta[1] =
              Theta_positions[1] +
              (Theta_positions[2] - Theta_positions[1]) * 0.5 * (s[0] + 1.0);

            Area_pt->position(t, zeta, f);
            break;

          default:
            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect << " not one of N, S, W, E"
                         << std::endl;

            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 3: Top
      case 3:

        // Which direction?
        switch (idirect)
        {
          case N:
            // Get the position on the wall by linearly interpolating in theta
            zeta[0] = 1.0;
            zeta[1] =
              Theta_positions[3] +
              (Theta_positions[2] - Theta_positions[3]) * 0.5 * (s[0] + 1.0);

            Area_pt->position(t, zeta, f);
            break;

          case S:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[3], Box[1], s[0], f);
            break;

          case W:
            // Linearly interpolate between the box and the wall
            lin_interpolate(Box[3], Wall[3], s[0], f);
            break;

          case E:
            // Linearly interpolate between the box and the wall
            lin_interpolate(Box[2], Wall[2], s[0], f);
            break;

          default:
            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                         << std::endl;

            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;


        // Macro element 4: Left
      case 4:

        // Which direction?
        switch (idirect)
        {
          case N:
            // Linearly interpolate between the wall and the box
            lin_interpolate(Wall[3], Box[3], s[0], f);
            break;

          case S:
            // Linearly interpolate between the wall and the box
            lin_interpolate(Wall[0], Box[0], s[0], f);
            break;

          case W:
            // Entirely on the wall, Need to be careful
            // because this is the point at which theta passes through the
            //"branch cut"
            zeta[0] = 1.0;
            zeta[1] = Theta_positions[0] + 2.0 * pi +
                      (Theta_positions[3] - (Theta_positions[0] + 2.0 * pi)) *
                        0.5 * (s[0] + 1.0);

            Area_pt->position(t, zeta, f);
            break;

          case E:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[0], Box[3], s[0], f);
            break;

          default:
            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect << " not one of N, S, W, E"
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

} // namespace oomph

#endif
