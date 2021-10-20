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
#ifndef OOMPH_TUBE_DOMAIN_HEADER
#define OOMPH_TUBE_DOMAIN_HEADER

// Generic oomph-lib includes
#include "../generic/quadtree.h"
#include "../generic/domain.h"
#include "../generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  ///  Tube as a domain. The entire domain must be defined by
  /// a GeomObject with the following convention: zeta[0] is the coordinate
  /// along the centreline, zeta[1] is the theta coordinate around the tube
  /// wall and zeta[2] is the radial coordinate.
  /// The outer boundary must lie at zeta[2] = 1.
  ///
  /// The domain is
  /// parametrised by five macro elements (a central box surrounded by
  /// four curved elements) in each of the nlayer slices. The labelling
  /// of the macro elements is as follows with the zeta[0] coordinate
  /// coming out of the page.
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
  class TubeDomain : public Domain
  {
  public:
    ///  Constructor: Pass geometric  object; start and end limit of the
    /// centreline coordinate; the theta locations marking the division between
    /// the elements of the outer ring, labelled from the lower left to the
    /// upper left in order, theta should be in the range \f$-\pi\f$ to
    /// \f$\pi\f$; the corresponding fractions of the
    /// radius at which the central box is to be placed; and the number of
    /// layers in the domain
    TubeDomain(GeomObject* volume_geom_object_pt,
               const Vector<double>& centreline_limits,
               const Vector<double>& theta_positions,
               const Vector<double>& radius_box,
               const unsigned& nlayer)
      : Centreline_limits(centreline_limits),
        Theta_positions(theta_positions),
        Radius_box(radius_box),
        Nlayer(nlayer),
        Volume_pt(volume_geom_object_pt)
    {
      // There are five macro elements
      const unsigned n_macro = 5 * nlayer;
      Macro_element_pt.resize(n_macro);

      // Create the macro elements
      for (unsigned i = 0; i < n_macro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<3>(this, i);
      }
    }

    /// Broken copy constructor
    TubeDomain(const TubeDomain&) = delete;

    /// Broken assignment operator
    void operator=(const TubeDomain&) = delete;

    /// Destructor: Empty; cleanup done in base class
    ~TubeDomain() {}

    ///  Vector representation of the  i_macro-th macro element
    /// boundary i_direct (L/R/D/U/B/F) at time level t
    /// (t=0: present; t>0: previous):
    /// f(s).
    void macro_element_boundary(const unsigned& t,
                                const unsigned& i_macro,
                                const unsigned& i_direct,
                                const Vector<double>& s,
                                Vector<double>& f);

  private:
    /// Storage for the limits of the centreline coordinate
    Vector<double> Centreline_limits;

    ///  Storage for the dividing lines on the boundary
    /// starting from the lower left and proceeding anticlockwise to
    /// the upper left
    Vector<double> Theta_positions;

    /// Storage for the fraction of the radius at which the central box
    /// should be located corresponding to the chosen values of theta.
    Vector<double> Radius_box;

    /// Number of axial layers
    unsigned Nlayer;

    /// Pointer to geometric object that represents the domain
    GeomObject* Volume_pt;

    ///   A very little linear interpolation helper.
    /// Interpolate from the low
    /// point to the high point using the coordinate s which is
    /// assumed to run from -1 to 1.
    void lin_interpolate(const Vector<double>& low,
                         const Vector<double>& high,
                         const double& s,
                         Vector<double>& f)
    {
      // Loop over all coordinates
      for (unsigned i = 0; i < 3; i++)
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
  /// boundary idirect (L/R/D/U/B/F) at time level t
  /// (t=0: present; t>0: previous): f(s)
  //=================================================================
  void TubeDomain::macro_element_boundary(const unsigned& t,
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
      "TubeDomain::macro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    // Get the number of the layer
    unsigned ilayer = unsigned(imacro / 5);

    // Get all required coordinates of the corners of the box
    // at each end of the layer
    Vector<Vector<Vector<double>>> Box(2);

    // Get the centreline coordinates at the ends of the layer
    Vector<double> zeta_centre(2);
    // Storage for position in the volume
    Vector<double> zeta(3);

    // Loop over the layers
    for (unsigned i = 0; i < 2; i++)
    {
      // Resize the storage
      Box[i].resize(4);

      // Get the centreline coordinate
      zeta_centre[i] =
        Centreline_limits[0] + (ilayer + i) *
                                 (Centreline_limits[1] - Centreline_limits[0]) /
                                 (double)(Nlayer);

      // Now get the coordinates of each corner
      zeta[0] = zeta_centre[i];

      // Loop over the angles
      for (unsigned j = 0; j < 4; j++)
      {
        // Set up the storage
        Box[i][j].resize(3);

        // Set the other values of zeta
        zeta[1] = Theta_positions[j];
        zeta[2] = Radius_box[j];

        // Now get the values
        Volume_pt->position(t, zeta, Box[i][j]);
      }
    }

    // Local storage for positions on the boundaries
    Vector<double> pos_1(3), pos_2(3);

    const double pi = MathematicalConstants::Pi;

    // Which macro element?
    // --------------------
    switch (imacro % 5)
    {
        // Macro element 0: Central box
      case 0:

        // Choose a direction
        switch (idirect)
        {
          case L:
            // Need to get the position from the domain
            // Get the centreline position
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            // Get the lower corner
            zeta[1] = Theta_positions[0];
            zeta[2] = Radius_box[0];
            Volume_pt->position(t, zeta, pos_1);

            // Get the upper corner
            zeta[1] = Theta_positions[3];
            zeta[2] = Radius_box[3];
            Volume_pt->position(t, zeta, pos_2);

            // Now interpolate between the two corner positions
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;

          case R:
            // Need to get the position from the domain
            // Get the centreline position
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            // Get the lower corner
            zeta[1] = Theta_positions[1];
            zeta[2] = Radius_box[1];
            Volume_pt->position(t, zeta, pos_1);

            // Get the upper corner
            zeta[1] = Theta_positions[2];
            zeta[2] = Radius_box[2];
            Volume_pt->position(t, zeta, pos_2);

            // Now interpolate between the positions
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;

          case D:
            // Need to get the position from the domain
            // Get the centreline position
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            // Get the left corner
            zeta[1] = Theta_positions[0];
            zeta[2] = Radius_box[0];
            Volume_pt->position(t, zeta, pos_1);

            // Get the right corner
            zeta[1] = Theta_positions[1];
            zeta[2] = Radius_box[1];
            Volume_pt->position(t, zeta, pos_2);
            // Now interpolate between the positions
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;

          case U:
            // Need to get the position from the domain
            // Get the centreline position
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            // Get the left corner
            zeta[1] = Theta_positions[3];
            zeta[2] = Radius_box[3];
            Volume_pt->position(t, zeta, pos_1);

            // Get the right corner
            zeta[1] = Theta_positions[2];
            zeta[2] = Radius_box[2];
            Volume_pt->position(t, zeta, pos_2);

            // Now interpolate between the positions
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;

          case B:
            // Linearly interpolate between lower left and lower right corners
            lin_interpolate(Box[0][0], Box[0][1], s[0], pos_1);
            // Linearly interpolate between upper left and upper right corners
            lin_interpolate(Box[0][3], Box[0][2], s[0], pos_2);
            // Finally, linearly interpolate between the upper and lower
            // positions
            lin_interpolate(pos_1, pos_2, s[1], f);
            break;

          case F:
            // Linearly interpolate between lower left and lower right corners
            lin_interpolate(Box[1][0], Box[1][1], s[0], pos_1);
            // Linearly interpolate between upper left and upper right corners
            lin_interpolate(Box[1][3], Box[1][2], s[0], pos_2);
            // Finally, linearly interpolate between the upper and lower
            // positions
            lin_interpolate(pos_1, pos_2, s[1], f);
            break;

          default:
            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect
                         << " not one of L, R, D, U, B, F" << std::endl;

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
          case L:
            // Get the position on the wall
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] = Theta_positions[0];
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Get the position on the box
            zeta[2] = Radius_box[0];
            Volume_pt->position(t, zeta, pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;

          case R:
            // Get the position on the wall
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] = Theta_positions[1];
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Get the position from the box
            zeta[2] = Radius_box[1];
            Volume_pt->position(t, zeta, pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;

          case D:
            // This is entrirly on the wall
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] =
              Theta_positions[0] +
              (Theta_positions[1] - Theta_positions[0]) * 0.5 * (s[0] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, f);
            break;

          case U:
            // Need to get the position from the domain
            // Get the centreline position
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            // Get the left corner
            zeta[1] = Theta_positions[0];
            zeta[2] = Radius_box[0];
            Volume_pt->position(t, zeta, pos_1);

            // Get the right corner
            zeta[1] = Theta_positions[1];
            zeta[2] = Radius_box[1];
            Volume_pt->position(t, zeta, pos_2);
            // Now interpolate between the positions
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;

          case B:
            // Get the position on the wall
            zeta[0] = zeta_centre[0];
            zeta[1] =
              Theta_positions[0] +
              (Theta_positions[1] - Theta_positions[0]) * 0.5 * (s[0] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Now linearly interpolate the position on the box
            lin_interpolate(Box[0][0], Box[0][1], s[0], pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_1, pos_2, s[1], f);
            break;

          case F:
            // Get the position on the wall
            zeta[0] = zeta_centre[1];
            zeta[1] =
              Theta_positions[0] +
              (Theta_positions[1] - Theta_positions[0]) * 0.5 * (s[0] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Now linearly interpolate the position on the box
            lin_interpolate(Box[1][0], Box[1][1], s[0], pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_1, pos_2, s[1], f);
            break;


          default:

            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect
                         << " not one of L, R, D, U, B, F" << std::endl;

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
          case L:
            // Need to get the position from the domain
            // Get the centreline position
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            // Get the lower corner
            zeta[1] = Theta_positions[1];
            zeta[2] = Radius_box[1];
            Volume_pt->position(t, zeta, pos_1);

            // Get the upper corner
            zeta[1] = Theta_positions[2];
            zeta[2] = Radius_box[2];
            Volume_pt->position(t, zeta, pos_2);
            // Now interpolate between the positions
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;

          case R:
            // Entirely on the wall
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] =
              Theta_positions[1] +
              (Theta_positions[2] - Theta_positions[1]) * 0.5 * (s[0] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, f);
            break;

          case U:
            // Get the position on the wall
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] = Theta_positions[2];
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Get the position of the box
            zeta[2] = Radius_box[2];
            Volume_pt->position(t, zeta, pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_2, pos_1, s[0], f);
            break;

          case D:
            // Get the position on the wall
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] = Theta_positions[1];
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Get the position of the box
            zeta[2] = Radius_box[1];
            Volume_pt->position(t, zeta, pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_2, pos_1, s[0], f);
            break;

          case B:
            // Get the position on the wall
            zeta[0] = zeta_centre[0];
            zeta[1] =
              Theta_positions[1] +
              (Theta_positions[2] - Theta_positions[1]) * 0.5 * (s[1] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Now linearly interpolate the position on the box
            lin_interpolate(Box[0][1], Box[0][2], s[1], pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_2, pos_1, s[0], f);
            break;

          case F:
            // Get the position on the wall
            zeta[0] = zeta_centre[1];
            zeta[1] =
              Theta_positions[1] +
              (Theta_positions[2] - Theta_positions[1]) * 0.5 * (s[1] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Now linearly interpolate the position on the box
            lin_interpolate(Box[1][1], Box[1][2], s[1], pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_2, pos_1, s[0], f);
            break;


          default:
            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect
                         << " not one of L, R, D, U, B, F" << std::endl;

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
          case L:
            // Get the position on the wall
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] = Theta_positions[3];
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Get the position on the box
            zeta[2] = Radius_box[3];
            Volume_pt->position(t, zeta, pos_2);


            // Now linearly interpolate between the two
            lin_interpolate(pos_2, pos_1, s[0], f);
            break;

          case R:
            // Get the position on the wall
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] = Theta_positions[2];
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Get the position on the box
            zeta[2] = Radius_box[2];
            Volume_pt->position(t, zeta, pos_2);


            // Now linearly interpolate between the two
            lin_interpolate(pos_2, pos_1, s[0], f);
            break;

          case D:
            // This is entirely on the box
            // Need to get the position from the domain
            // Get the centreline position
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] = Theta_positions[3];
            zeta[2] = Radius_box[3];
            // Get the lower corner
            Volume_pt->position(t, zeta, pos_2);

            // Get the upper corner
            zeta[1] = Theta_positions[2];
            zeta[2] = Radius_box[2];
            // Get the upper corner
            Volume_pt->position(t, zeta, pos_1);
            // Now interpolate between the positions
            lin_interpolate(pos_2, pos_1, s[0], f);
            break;

          case U:
            // This is entirely on the wall
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] =
              Theta_positions[3] +
              (Theta_positions[2] - Theta_positions[3]) * 0.5 * (s[0] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, f);
            break;

          case B:
            // Get the position on the wall
            zeta[0] = zeta_centre[0];
            zeta[1] =
              Theta_positions[3] +
              (Theta_positions[2] - Theta_positions[3]) * 0.5 * (s[0] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Now linearly interpolate the position on the box
            lin_interpolate(Box[0][3], Box[0][2], s[0], pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_2, pos_1, s[1], f);
            break;

          case F:
            // Get the position on the wall
            zeta[0] = zeta_centre[1];
            zeta[1] =
              Theta_positions[3] +
              (Theta_positions[2] - Theta_positions[3]) * 0.5 * (s[0] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Now linearly interpolate the position on the box
            lin_interpolate(Box[1][3], Box[1][2], s[0], pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_2, pos_1, s[1], f);
            break;


          default:
            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect
                         << " not one of L, R, D, U, B, F" << std::endl;

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
          case L:
            // Entirely on the wall, Need to be careful
            // because this is the point at which theta passes through the
            // branch cut
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] = Theta_positions[0] + 2.0 * pi +
                      (Theta_positions[3] - (Theta_positions[0] + 2.0 * pi)) *
                        0.5 * (s[0] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, f);
            break;

          case R:
            // Entirely on the box
            // Need to get the position from the domain
            // Get the centreline position
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] = Theta_positions[0];
            zeta[2] = Radius_box[0];
            // Get the lower corner
            Volume_pt->position(t, zeta, pos_2);

            // Get the upper corner
            zeta[1] = Theta_positions[3];
            zeta[2] = Radius_box[3];
            // Get the upper corner
            Volume_pt->position(t, zeta, pos_1);

            lin_interpolate(pos_2, pos_1, s[0], f);
            break;

          case D:
            // Get the position on the wall
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] = Theta_positions[0];
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);


            // Get the position on the box
            zeta[2] = Radius_box[0];
            Volume_pt->position(t, zeta, pos_2);


            // Now linearly interpolate between the two
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;

          case U:
            // Get the position on the wall
            zeta[0] = zeta_centre[0] +
                      (zeta_centre[1] - zeta_centre[0]) * 0.5 * (s[1] + 1.0);
            zeta[1] = Theta_positions[3];
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Get the position on the box
            zeta[2] = Radius_box[3];
            Volume_pt->position(t, zeta, pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;

          case B:
            // Get the position on the wall
            // Again be careful of the branch cut
            zeta[0] = zeta_centre[0];
            zeta[1] = Theta_positions[0] + 2.0 * pi +
                      (Theta_positions[3] - (Theta_positions[0] + 2.0 * pi)) *
                        0.5 * (s[1] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Now linearly interpolate the position on the box
            lin_interpolate(Box[0][0], Box[0][3], s[1], pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;


          case F:
            // Get the position on the wall
            // Again be careful of the branch cut
            zeta[0] = zeta_centre[1];
            zeta[1] = Theta_positions[0] + 2.0 * pi +
                      (Theta_positions[3] - (Theta_positions[0] + 2.0 * pi)) *
                        0.5 * (s[1] + 1.0);
            zeta[2] = 1.0;
            Volume_pt->position(t, zeta, pos_1);

            // Now linearly interpolate the position on the box
            lin_interpolate(Box[1][0], Box[1][3], s[1], pos_2);

            // Now linearly interpolate between the two
            lin_interpolate(pos_1, pos_2, s[0], f);
            break;


          default:
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

} // namespace oomph

#endif
