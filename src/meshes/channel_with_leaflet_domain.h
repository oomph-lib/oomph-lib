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
#ifndef OOMPH_CHANNEL_WITH_LEAFLET_DOMAIN_HEADER
#define OOMPH_CHANNEL_WITH_LEAFLET_DOMAIN_HEADER


// Generic includes
#include "../generic/geom_objects.h"
#include "../generic/domain.h"
#include "../generic/macro_element.h"

namespace oomph
{
  //===========================================================
  /// Rectangular domain with a leaflet blocking the lower
  /// half
  //===========================================================
  class ChannelWithLeafletDomain : public Domain
  {
  public:
    /// Constructor: Pass pointer to GeomObject that represents the
    /// leaflet,
    /// the length of the domain to left and right of the leaflet, the
    /// height of the leaflet and the overall height of the channel,
    /// the number of element columns to the left and right of the leaflet,
    /// the number of rows of elements from the bottom of the channel to
    /// the end of the leaflet, the number of rows of elements above the
    /// end of the leaflet.
    ChannelWithLeafletDomain(GeomObject* leaflet_pt,
                             const double& lleft,
                             const double& lright,
                             const double& hleaflet,
                             const double& htot,
                             const unsigned& nleft,
                             const unsigned& nright,
                             const unsigned& ny1,
                             const unsigned& ny2)
    {
      // Copy assignments into private storage
      Lleft = lleft;
      Lright = lright;
      Hleaflet = hleaflet;
      Htot = htot;
      Nleft = nleft;
      Nright = nright;
      Ny1 = ny1;
      Ny2 = ny2;
      Leaflet_pt = leaflet_pt;

      // Store origin of leaflet for fast reference
      Vector<double> zeta(1);
      zeta[0] = 0.0;
      Vector<double> r(2);
      Leaflet_pt->position(zeta, r);
      X_0 = r[0];

      // Number of macro elements
      unsigned nmacro = (Ny1 + Ny2) * (Nleft + Nright);
      Macro_element_pt.resize(nmacro);

      // Build the 2D macro elements
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    }

    /// Destructor: Empty; cleanup done in base class
    ~ChannelWithLeafletDomain() {}

    /// Total height of domain (width of channel)
    double htot()
    {
      return Htot;
    }

    /// Height of leaflet
    double hleaflet()
    {
      return Hleaflet;
    }

    /// Length of domain to the left of leaflet
    double lleft()
    {
      return Lleft;
    }

    /// Length of domain to the right of leaflet
    double lright()
    {
      return Lright;
    }

    /// Pointer to the wall
    GeomObject*& leaflet_pt()
    {
      return Leaflet_pt;
    };

    /// Parametrisation of macro element boundaries
    void macro_element_boundary(const unsigned& t,
                                const unsigned& imacro,
                                const unsigned& idirect,
                                const Vector<double>& zeta,
                                Vector<double>& r);

  protected:
    /// Helper function
    void macro_bound_I_N(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& i,
                         const unsigned& j);

    /// Helper function
    void macro_bound_I_S(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& i,
                         const unsigned& j);

    /// Helper function
    void macro_bound_I_W(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& i,
                         const unsigned& j);

    /// Helper function
    void macro_bound_I_E(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& i,
                         const unsigned& j);

    /// Helper function
    void macro_bound_II_N(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_II_S(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_II_W(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_II_E(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_III_N(const unsigned& t,
                           const Vector<double>& zeta,
                           Vector<double>& r,
                           const unsigned& i,
                           const unsigned& j);

    /// Helper function
    void macro_bound_III_S(const unsigned& t,
                           const Vector<double>& zeta,
                           Vector<double>& r,
                           const unsigned& i,
                           const unsigned& j);

    /// Helper function
    void macro_bound_III_W(const unsigned& t,
                           const Vector<double>& zeta,
                           Vector<double>& r,
                           const unsigned& i,
                           const unsigned& j);

    /// Helper function
    void macro_bound_III_E(const unsigned& t,
                           const Vector<double>& zeta,
                           Vector<double>& r,
                           const unsigned& i,
                           const unsigned& j);

    /// Helper function
    void macro_bound_IV_N(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_IV_S(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_IV_W(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_IV_E(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);
    /// Helper function
    void slanted_bound_up(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r);


    /// Length of the domain to the right of the leaflet
    double Lright;

    /// Length of the domain to the left of the leaflet
    double Lleft;

    /// Lagrangian coordinate at end of leaflet
    double Hleaflet;

    /// Total width of the channel
    double Htot;

    /// Number of macro element columnns to the right of the leaflet
    unsigned Nright;

    /// Number of macro element columns to the left of the leaflet
    unsigned Nleft;

    /// Number of macro element rows up to the end of the leaflet
    unsigned Ny1;

    /// Number of macro element rows above the leaflet
    unsigned Ny2;

    /// Center of the domain : origin of the leaflet, extracted
    /// from GeomObject and stored for fast access.
    double X_0;

    /// Pointer to leaflet
    GeomObject* Leaflet_pt;
  };


  //===================================================================
  /// Parametrisation of macro element boundaries
  //===================================================================
  void ChannelWithLeafletDomain::macro_element_boundary(
    const unsigned& t,
    const unsigned& imacro,
    const unsigned& idirect,
    const Vector<double>& zeta,
    Vector<double>& r)
  {
#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "ChannelWithLeafletDomain::macro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    using namespace QuadTreeNames;

    // Number of the line and of the colum in the whole domain
    // Beware : iline starts on zero
    unsigned iline = int(floor(double(imacro) / double(Nleft + Nright)));
    unsigned icol = imacro % (Nright + Nleft);

    // Number of the line and of the colum in the part considered
    unsigned i, j;

    // Left low part of the domain : Part I
    //------------------------------------
    if ((iline < Ny1) && (icol < Nleft))
    {
      i = iline;
      j = icol;

      switch (idirect)
      {
        case N:
          macro_bound_I_N(t, zeta, r, i, j);
          break;
        case S:
          macro_bound_I_S(t, zeta, r, i, j);
          break;
        case W:
          macro_bound_I_W(t, zeta, r, i, j);
          break;
        case E:
          macro_bound_I_E(t, zeta, r, i, j);
          break;
      }
    }
    // Right low part of the domain : Part II
    //--------------------------------------
    else if ((iline < Ny1) && (icol >= Nleft))
    {
      i = iline;
      j = icol - Nleft;

      switch (idirect)
      {
        case N:
          macro_bound_II_N(t, zeta, r, i, j);
          break;
        case S:
          macro_bound_II_S(t, zeta, r, i, j);
          break;
        case W:
          macro_bound_II_W(t, zeta, r, i, j);
          break;
        case E:
          macro_bound_II_E(t, zeta, r, i, j);
          break;
      }
    }
    // Left upper part of the domain : Part III
    //----------------------------------------
    else if ((iline >= Ny1) && (icol < Nleft))
    {
      i = iline - Ny1;
      j = icol;

      switch (idirect)
      {
        case N:
          macro_bound_III_N(t, zeta, r, i, j);
          break;
        case S:
          macro_bound_III_S(t, zeta, r, i, j);
          break;
        case W:
          macro_bound_III_W(t, zeta, r, i, j);
          break;
        case E:
          macro_bound_III_E(t, zeta, r, i, j);
          break;
      }
    }
    // Right upper part of the domain : Part IV
    //-----------------------------------------
    else if ((iline >= Ny1) && (icol >= Nleft))
    {
      i = iline - Ny1;
      j = icol - Nleft;

      switch (idirect)
      {
        case N:
          macro_bound_IV_N(t, zeta, r, i, j);
          break;
        case S:
          macro_bound_IV_S(t, zeta, r, i, j);
          break;
        case W:
          macro_bound_IV_W(t, zeta, r, i, j);
          break;
        case E:
          macro_bound_IV_E(t, zeta, r, i, j);
          break;
      }
    }
  }


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  // Helper functions for region I (lower left region)
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //=====================================================================
  /// Helper function for eastern boundary in lower left region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_I_E(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& i,
                                                 const unsigned& j)
  {
    // Find x,y on the wall corresponding to the position of the macro element

    // xi_wall varies from xi0 to xi1 on the wall
    double xi0, xi1;
    xi0 = double(i) * Hleaflet / double(Ny1);
    xi1 = double(i + 1) * Hleaflet / double(Ny1);

    Vector<double> xi_wall(1);
    xi_wall[0] = xi0 + (1.0 + zeta[0]) / 2.0 * (xi1 - xi0);

    Vector<double> r_wall(2);
    Leaflet_pt->position(t, xi_wall, r_wall);

    // Find x,y on a vertical line corresponding
    // to the position of the macro element

    // the vertical line goes from y0 to y1
    double y0, y1;
    y0 = double(i) * Hleaflet / double(Ny1);
    y1 = double(i + 1) * Hleaflet / double(Ny1);

    Vector<double> r_vert(2);
    r_vert[0] = -Lleft + X_0;
    r_vert[1] = y0 + (1.0 + zeta[0]) / 2.0 * (y1 - y0);

    // Parameter with value 0 in -Lleft and value 1 on the wall.
    double s = double(j + 1) / double(Nleft);

    /// Final expression of r
    r[0] = r_vert[0] + s * (r_wall[0] - r_vert[0]);
    r[1] = r_vert[1] + s * (r_wall[1] - r_vert[1]);
  }


  //=====================================================================
  /// Helper function for western boundary in lower left region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_I_W(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& i,
                                                 const unsigned& j)
  {
    // Find x,y on the wall corresponding to the position of the macro element

    // xi_wall varies from xi0 to xi1 on the wall
    double xi0, xi1;
    xi0 = double(i) * Hleaflet / double(Ny1);
    xi1 = double(i + 1) * Hleaflet / double(Ny1);

    Vector<double> xi_wall(1);
    xi_wall[0] = xi0 + (1.0 + zeta[0]) / 2.0 * (xi1 - xi0);

    Vector<double> r_wall(2);
    Leaflet_pt->position(t, xi_wall, r_wall);

    // Find x,y on a vertical line corresponding
    // to the position of the macro element

    // the vertical line goes from y0 to y1
    double y0, y1;
    y0 = double(i) * Hleaflet / double(Ny1);
    y1 = double(i + 1) * Hleaflet / double(Ny1);

    Vector<double> r_vert(2);
    r_vert[0] = -Lleft + X_0;
    r_vert[1] = y0 + (1.0 + zeta[0]) / 2.0 * (y1 - y0);

    // Parameter with value 0 in -Lleft and value 1 on the wall.
    double s = double(j) / double(Nleft);

    /// Final expression of r
    r[0] = r_vert[0] + s * (r_wall[0] - r_vert[0]);
    r[1] = r_vert[1] + s * (r_wall[1] - r_vert[1]);
  }


  //=====================================================================
  /// Helper function for northern boundary in lower left region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_I_N(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& i,
                                                 const unsigned& j)
  {
    // Find the coordinates of the two corners of the north boundary
    Vector<double> xi(1);
    Vector<double> r_left(2);
    Vector<double> r_right(2);
    xi[0] = 1;
    macro_bound_I_W(t, xi, r_left, i, j);
    macro_bound_I_E(t, xi, r_right, i, j);

    // Connect those two points with a straight line
    r[0] = r_left[0] + (1.0 + zeta[0]) / 2.0 * (r_right[0] - r_left[0]);
    r[1] = r_left[1] + (1.0 + zeta[0]) / 2.0 * (r_right[1] - r_left[1]);
  }


  //=====================================================================
  /// Helper function for southern boundary in lower left region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_I_S(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& i,
                                                 const unsigned& j)
  {
    /// Find the coordinates of the two corners of the south boundary
    Vector<double> xi(1);
    Vector<double> r_left(2);
    Vector<double> r_right(2);
    xi[0] = -1.0;
    macro_bound_I_W(t, xi, r_left, i, j);
    macro_bound_I_E(t, xi, r_right, i, j);

    // Connect those two points with a straight line
    r[0] = r_left[0] + (1.0 + zeta[0]) / 2.0 * (r_right[0] - r_left[0]);
    r[1] = r_left[1] + (1.0 + zeta[0]) / 2.0 * (r_right[1] - r_left[1]);
  }


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  // Helper functions for region II (lower right region)
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //=====================================================================
  /// Helper function for eastern boundary in lower right region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_II_E(const unsigned& t,
                                                  const Vector<double>& zeta,
                                                  Vector<double>& r,
                                                  const unsigned& i,
                                                  const unsigned& j)

  {
    // Find x,y on the wall corresponding to the position of the macro element

    // xi_wall varies from xi0 to xi1 on the wall
    double xi0, xi1;
    xi0 = double(i) * Hleaflet / double(Ny1);
    xi1 = double(i + 1) * Hleaflet / double(Ny1);

    Vector<double> xi_wall(1);
    xi_wall[0] = xi0 + (1.0 + zeta[0]) / 2.0 * (xi1 - xi0);

    Vector<double> r_wall(2);
    Leaflet_pt->position(t, xi_wall, r_wall);

    // Find x,y on a vertical line corresponding
    // to the position of the macro element

    // the vertical line goes from y0 to y1
    double y0, y1;
    y0 = double(i) * Hleaflet / double(Ny1);
    y1 = double(i + 1) * Hleaflet / double(Ny1);

    Vector<double> r_vert(2);
    r_vert[0] = Lright + X_0;
    r_vert[1] = y0 + (1.0 + zeta[0]) / 2.0 * (y1 - y0);

    // Parameter with value 0 in Lright and value 1 on the wall.
    double s = double(Nright - j - 1) / double(Nright); /***Change****/

    /// Final expression of r
    r[0] = r_vert[0] + s * (r_wall[0] - r_vert[0]);
    r[1] = r_vert[1] + s * (r_wall[1] - r_vert[1]);
  }


  //=====================================================================
  /// Helper function for western boundary in lower right region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_II_W(const unsigned& t,
                                                  const Vector<double>& zeta,
                                                  Vector<double>& r,
                                                  const unsigned& i,
                                                  const unsigned& j)
  {
    // Abscissa of the origin of the boudary east

    // Find x,y on the wall corresponding to the position of the macro element

    // xi_wall varies from xi0 to xi1 on the wall
    double xi0, xi1;
    xi0 = double(i) * Hleaflet / double(Ny1);
    xi1 = double(i + 1) * Hleaflet / double(Ny1);

    Vector<double> xi_wall(1);
    xi_wall[0] = xi0 + (1.0 + zeta[0]) / 2.0 * (xi1 - xi0);

    Vector<double> r_wall(2);
    Leaflet_pt->position(t, xi_wall, r_wall);

    // Find x,y on a vertical line corresponding
    // to the position of the macro element

    // the vertical line goes from y0 to y1
    double y0, y1;
    y0 = double(i) * Hleaflet / double(Ny1);
    y1 = double(i + 1) * Hleaflet / double(Ny1);

    Vector<double> r_vert(2);
    r_vert[0] = Lright + X_0;
    r_vert[1] = y0 + (1.0 + zeta[0]) / 2.0 * (y1 - y0);

    // Parameter with value 0 in -Lleft and value 1 on the wall.
    double s = double(Nright - j) / double(Nright); /***Change****/

    // Final expression of r
    r[0] = r_vert[0] + s * (r_wall[0] - r_vert[0]);
    r[1] = r_vert[1] + s * (r_wall[1] - r_vert[1]);
  }


  //=====================================================================
  /// Helper function for northern boundary in lower right region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_II_N(const unsigned& t,
                                                  const Vector<double>& zeta,
                                                  Vector<double>& r,
                                                  const unsigned& i,
                                                  const unsigned& j)
  {
    // Find the coordinates of the two corners of the north boundary
    Vector<double> xi(1);
    Vector<double> r_left(2);
    Vector<double> r_right(2);
    xi[0] = 1;
    macro_bound_II_W(t, xi, r_left, i, j);
    macro_bound_II_E(t, xi, r_right, i, j);

    // Connect those two points with a straight line
    r[0] = r_left[0] + (1.0 + zeta[0]) / 2.0 * (r_right[0] - r_left[0]);
    r[1] = r_left[1] + (1.0 + zeta[0]) / 2.0 * (r_right[1] - r_left[1]);
  }


  //=====================================================================
  /// Helper function for southern boundary in lower right region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_II_S(const unsigned& t,
                                                  const Vector<double>& zeta,
                                                  Vector<double>& r,
                                                  const unsigned& i,
                                                  const unsigned& j)
  {
    // Find the coordinates of the two corners of the south boundary
    Vector<double> xi(1);
    Vector<double> r_left(2);
    Vector<double> r_right(2);
    xi[0] = -1.0;
    macro_bound_II_W(t, xi, r_left, i, j);
    macro_bound_II_E(t, xi, r_right, i, j);

    // Connect those two points with a straight line
    r[0] = r_left[0] + (1.0 + zeta[0]) / 2.0 * (r_right[0] - r_left[0]);
    r[1] = r_left[1] + (1.0 + zeta[0]) / 2.0 * (r_right[1] - r_left[1]);
  }


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  // Helper functions for region III (upper left region)
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //=====================================================================
  /// Describe the line between the boundary north of the domain (at x=X_0)
  /// and the top of the wall, when zeta goes from 0 to 1.
  //=====================================================================
  void ChannelWithLeafletDomain::slanted_bound_up(const unsigned& t,
                                                  const Vector<double>& zeta,
                                                  Vector<double>& r)
  {
    // Coordinates of the point on the boundary beetween the upper
    // and the lower part, in the same column, at the east.
    Vector<double> xi(1);
    xi[0] = Hleaflet;

    Vector<double> r_join(2);

    Leaflet_pt->position(t, xi, r_join);

    r[0] = r_join[0] + zeta[0] * (X_0 - r_join[0]);
    r[1] = r_join[1] + zeta[0] * (Htot - r_join[1]);
  }


  //=====================================================================
  /// Helper function for eastern boundary in upper left region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_III_E(const unsigned& t,
                                                   const Vector<double>& zeta,
                                                   Vector<double>& r,
                                                   const unsigned& i,
                                                   const unsigned& j)
  {
    // Find x,y on the slanted straight line (SSL) corresponding to
    // the position of the macro element

    // xi_line varies from xi0 to xi1 on the SSL
    double xi0, xi1;
    xi0 = double(i) / double(Ny2);
    xi1 = double(i + 1) / double(Ny2);

    Vector<double> xi_line(1);
    xi_line[0] = xi0 + (1.0 + zeta[0]) / 2.0 * (xi1 - xi0);

    Vector<double> r_line(2);
    slanted_bound_up(t, xi_line, r_line);

    // Find x,y on a vertical line corresponding
    // to the position of the macro element

    // the vertical line goes from y0 to y1
    double y0, y1;
    y0 = double(i) * (Htot - Hleaflet) / double(Ny2) + Hleaflet;
    y1 = double(i + 1) * (Htot - Hleaflet) / double(Ny2) + Hleaflet;
    ;

    Vector<double> r_vert(2);
    r_vert[0] = -Lleft + X_0;
    r_vert[1] = y0 + (1.0 + zeta[0]) / 2.0 * (y1 - y0);

    // Parameter with value 0 in Lright and value 1 on the wall.
    double s = double(j + 1) / double(Nleft); /***Change****/

    /// Final expression of r
    r[0] = r_vert[0] + s * (r_line[0] - r_vert[0]);
    r[1] = r_vert[1] + s * (r_line[1] - r_vert[1]);
  }


  //=====================================================================
  /// Helper function for western boundary in upper left region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_III_W(const unsigned& t,
                                                   const Vector<double>& zeta,
                                                   Vector<double>& r,
                                                   const unsigned& i,
                                                   const unsigned& j)
  {
    // Find x,y on the slanted straight line (SSL) corresponding to
    // the position of the macro element

    // xi_line varies from xi0 to xi1 on the SSL
    double xi0, xi1;
    xi0 = double(i) / double(Ny2);
    xi1 = double(i + 1) / double(Ny2);

    Vector<double> xi_line(1);
    xi_line[0] = xi0 + (1.0 + zeta[0]) / 2.0 * (xi1 - xi0);

    Vector<double> r_line(2);
    slanted_bound_up(t, xi_line, r_line);

    // Find x,y on a vertical line corresponding
    // to the position of the macro element

    // the vertical line goes from y0 to y1
    double y0, y1;
    y0 = double(i) * (Htot - Hleaflet) / double(Ny2) + Hleaflet;
    y1 = double(i + 1) * (Htot - Hleaflet) / double(Ny2) + Hleaflet;
    ;

    Vector<double> r_vert(2);
    r_vert[0] = -Lleft + X_0;
    r_vert[1] = y0 + (1.0 + zeta[0]) / 2.0 * (y1 - y0);

    // Parameter with value 0 in Lright and value 1 on the wall.
    double s = double(j) / double(Nleft); /***Change****/

    // Final expression of r
    r[0] = r_vert[0] + s * (r_line[0] - r_vert[0]);
    r[1] = r_vert[1] + s * (r_line[1] - r_vert[1]);
  }


  //=====================================================================
  /// Helper function for northern boundary in upper left region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_III_N(const unsigned& t,
                                                   const Vector<double>& zeta,
                                                   Vector<double>& r,
                                                   const unsigned& i,
                                                   const unsigned& j)
  {
    // Find the coordinates of the two corners of the north boundary
    Vector<double> xi(1);
    Vector<double> r_left(2);
    Vector<double> r_right(2);
    xi[0] = 1;
    macro_bound_III_W(t, xi, r_left, i, j);
    macro_bound_III_E(t, xi, r_right, i, j);

    // Connect those two points with a straight line
    r[0] = r_left[0] + (1.0 + zeta[0]) / 2.0 * (r_right[0] - r_left[0]);
    r[1] = r_left[1] + (1.0 + zeta[0]) / 2.0 * (r_right[1] - r_left[1]);
  }


  //=====================================================================
  /// Helper function for southern boundary in upper left region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_III_S(const unsigned& t,
                                                   const Vector<double>& zeta,
                                                   Vector<double>& r,
                                                   const unsigned& i,
                                                   const unsigned& j)
  {
    // Find the coordinates of the two corners of the south boundary
    Vector<double> xi(1);
    Vector<double> r_left(2);
    Vector<double> r_right(2);
    xi[0] = -1;
    macro_bound_III_W(t, xi, r_left, i, j);
    macro_bound_III_E(t, xi, r_right, i, j);

    // Connect those two points with a straight line
    r[0] = r_left[0] + (1.0 + zeta[0]) / 2.0 * (r_right[0] - r_left[0]);
    r[1] = r_left[1] + (1.0 + zeta[0]) / 2.0 * (r_right[1] - r_left[1]);
  }


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  // Helper functions for region IV (upper right region)
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //=====================================================================
  /// Helper function for eastern boundary in upper right region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_IV_E(const unsigned& t,
                                                  const Vector<double>& zeta,
                                                  Vector<double>& r,
                                                  const unsigned& i,
                                                  const unsigned& j)
  {
    // Find x,y on the slanted straight line (SSL) corresponding to
    // the position of the macro element

    // xi_line varies from xi0 to xi1 on the SSL
    double xi0, xi1;
    xi0 = double(i) / double(Ny2);
    xi1 = double(i + 1) / double(Ny2);

    Vector<double> xi_line(1);
    xi_line[0] = xi0 + (1.0 + zeta[0]) / 2.0 * (xi1 - xi0);

    Vector<double> r_line(2);
    slanted_bound_up(t, xi_line, r_line);

    // Find x,y on a vertical line corresponding
    // to the position of the macro element

    // the vertical line goes from y0 to y1
    double y0, y1;
    y0 = double(i) * (Htot - Hleaflet) / double(Ny2) + Hleaflet;
    y1 = double(i + 1) * (Htot - Hleaflet) / double(Ny2) + Hleaflet;
    ;

    Vector<double> r_vert(2);
    r_vert[0] = Lright + X_0;
    r_vert[1] = y0 + (1.0 + zeta[0]) / 2.0 * (y1 - y0);

    // Parameter with value 0 in Lright and value 1 on the wall.
    double s = double(Nright - j - 1) / double(Nright); /***Change****/

    // Final expression of r
    r[0] = r_vert[0] + s * (r_line[0] - r_vert[0]);
    r[1] = r_vert[1] + s * (r_line[1] - r_vert[1]);
  }


  //=====================================================================
  /// Helper function for western boundary in upper right region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_IV_W(const unsigned& t,
                                                  const Vector<double>& zeta,
                                                  Vector<double>& r,
                                                  const unsigned& i,
                                                  const unsigned& j)
  {
    // Find x,y on the slanted straight line (SSL) corresponding to
    // the position of the macro element

    // xi_line varies from xi0 to xi1 on the SSL
    double xi0, xi1;
    xi0 = double(i) / double(Ny2);
    xi1 = double(i + 1) / double(Ny2);

    Vector<double> xi_line(1);
    xi_line[0] = xi0 + (1.0 + zeta[0]) / 2.0 * (xi1 - xi0);

    Vector<double> r_line(2);
    slanted_bound_up(t, xi_line, r_line);

    // Find x,y on a vertical line corresponding
    // to the position of the macro element

    // The vertical line goes from y0 to y1
    double y0, y1;
    y0 = double(i) * (Htot - Hleaflet) / double(Ny2) + Hleaflet;
    y1 = double(i + 1) * (Htot - Hleaflet) / double(Ny2) + Hleaflet;
    ;

    Vector<double> r_vert(2);
    r_vert[0] = Lright + X_0;
    r_vert[1] = y0 + (1.0 + zeta[0]) / 2.0 * (y1 - y0);

    // Parameter with value 0 in Lright and value 1 on the wall.
    double s = double(Nright - j) / double(Nright); /***Change****/

    // Final expression of r
    r[0] = r_vert[0] + s * (r_line[0] - r_vert[0]);
    r[1] = r_vert[1] + s * (r_line[1] - r_vert[1]);
  }


  //=====================================================================
  /// Helper function for northern boundary in upper right region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_IV_N(const unsigned& t,
                                                  const Vector<double>& zeta,
                                                  Vector<double>& r,
                                                  const unsigned& i,
                                                  const unsigned& j)
  {
    // Find the coordinates of the two corners of the north boundary
    Vector<double> xi(1);
    Vector<double> r_left(2);
    Vector<double> r_right(2);
    xi[0] = 1;
    macro_bound_IV_W(t, xi, r_left, i, j);
    macro_bound_IV_E(t, xi, r_right, i, j);

    // Connect those two points with a straight line
    r[0] = r_left[0] + (1.0 + zeta[0]) / 2.0 * (r_right[0] - r_left[0]);
    r[1] = r_left[1] + (1.0 + zeta[0]) / 2.0 * (r_right[1] - r_left[1]);
  }


  //=====================================================================
  /// Helper function for southern boundary in upper right region
  //=====================================================================
  void ChannelWithLeafletDomain::macro_bound_IV_S(const unsigned& t,
                                                  const Vector<double>& zeta,
                                                  Vector<double>& r,
                                                  const unsigned& i,
                                                  const unsigned& j)
  {
    // Find the coordinates of the two corners of the south boundary
    Vector<double> xi(1);
    Vector<double> r_left(2);
    Vector<double> r_right(2);
    xi[0] = -1;
    macro_bound_IV_W(t, xi, r_left, i, j);
    macro_bound_IV_E(t, xi, r_right, i, j);

    // Connect those two points with a straight line
    r[0] = r_left[0] + (1.0 + zeta[0]) / 2.0 * (r_right[0] - r_left[0]);
    r[1] = r_left[1] + (1.0 + zeta[0]) / 2.0 * (r_right[1] - r_left[1]);
  }


} // namespace oomph


#endif
