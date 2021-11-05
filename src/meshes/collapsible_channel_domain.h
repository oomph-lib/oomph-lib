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
#ifndef OOMPH_COLLAPSIBLE_CHANNEL_DOMAIN_HEADER
#define OOMPH_COLLAPSIBLE_CHANNEL_DOMAIN_HEADER

// Generic oomph-lib includes
#include "../generic/quadtree.h"
#include "../generic/domain.h"
#include "../generic/geom_objects.h"

namespace oomph
{
  //==================================================================
  /// Collapsible channel domain
  //==================================================================
  class CollapsibleChannelDomain : public Domain
  {
  public:
    /// Constructor: Pass the number of (macro-)elements,
    /// the domain lengths in the x- and y-direction
    /// and the pointer to the geometric object that specifies
    /// the shape of the "collapsible" segment.
    CollapsibleChannelDomain(const unsigned& nup,
                             const unsigned& ncollapsible,
                             const unsigned& ndown,
                             const unsigned& ny,
                             const double& lup,
                             const double& lcollapsible,
                             const double& ldown,
                             const double& ly,
                             GeomObject* wall_pt)
      : BL_squash_fct_pt(&default_BL_squash_fct),
        Axial_spacing_fct_pt(&default_axial_spacing_fct),
        Rotate_domain(false)
    {
      Nup = nup;
      Ncollapsible = ncollapsible;
      Ndown = ndown;
      Ny = ny;
      Lup = lup;
      Lcollapsible = lcollapsible;
      Ldown = ldown;
      Ly = ly;
      Wall_pt = wall_pt;

      // Total number of macro elements
      unsigned nmacro = (Nup + Ncollapsible + Ndown) * Ny;

      // Build the macro elements
      Macro_element_pt.resize(nmacro);
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    }


    /// Destructor: emtpy; cleanup done in base class
    ~CollapsibleChannelDomain() {}


    /// Number of vertical columns of macro elements the upstream section
    unsigned nup()
    {
      return Nup;
    }

    /// Number of vertical clumns of macro elements in the "collapsible" segment
    unsigned ncollapsible()
    {
      return Ncollapsible;
    }

    /// Number of vertical columns of macro elements in the downstream section
    unsigned ndown()
    {
      return Ndown;
    }

    /// Number of macro-elements across the channel
    unsigned ny()
    {
      return Ny;
    }

    /// Length of upstream section
    double l_up()
    {
      return Lup;
    }

    /// Length of collapsible segment
    double l_collapsible()
    {
      return Lcollapsible;
    }

    /// Length of downstream section
    double l_down()
    {
      return Ldown;
    }

    /// Width of channel
    double l_y()
    {
      return Ly;
    }

    /// Access to pointer to the geometric object that parametrises
    /// the collapsible wall
    GeomObject*& wall_pt()
    {
      return Wall_pt;
    }


    /// Access to pointer to the geometric object that parametrises
    /// the collapsible wall (const version)
    GeomObject* wall_pt() const
    {
      return Wall_pt;
    }

    /// Typedef for function pointer for function that squashes
    /// the macro elements near the wall to help resolution of any
    /// wall boundary layers.
    typedef double (*BLSquashFctPt)(const double& s);


    /// Default for function that squashes
    /// the macro elements near the walls. Identity.
    static double default_BL_squash_fct(const double& s)
    {
      return s;
    }

    /// Function pointer for function that squashes
    /// the macro elements near wall. Default mapping (identity)
    /// leaves the y-coordinate of the nodal points unchanged.
    BLSquashFctPt& bl_squash_fct_pt()
    {
      return BL_squash_fct_pt;
    }

    /// Function that squashes the macro elements near the wall.
    /// Input argument should vary between 0 and 1; function should return
    /// stretched/squashed coordinate in the same range. Default implementation
    /// is the identity; can be overloaded by specifying a different
    /// function pointer with bl_squash_fct_pt().
    double s_squash(const double& s)
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


    /// Vector representation of the  imacro-th macro element
    /// boundary idirect (N/S/W/E) at time level t
    /// (t=0: present; t>0: previous): \f$ {\bf r}({\bf zeta}) \f$
    /// Note that the local coordinate \b zeta is a 1D
    /// Vector rather than a scalar -- this is unavoidable because
    /// this function implements the pure virtual function in the
    /// Domain base class.
    void macro_element_boundary(const unsigned& t,
                                const unsigned& imacro,
                                const unsigned& idirect,
                                const Vector<double>& zeta,
                                Vector<double>& r);


    /// Rotate the domain (for axisymmetric problems)
    void enable_rotate_domain()
    {
      Rotate_domain = true;
    }

    /// Undo rotation of the domain (for axisymmetric problems)
    void disable_rotate_domain()
    {
      Rotate_domain = false;
    }


  private:
    /// Northern boundary of the macro element imacro in the
    /// upstream (part=0) or downstream (part=1) sections
    void r_N_straight(const Vector<double>& zeta,
                      Vector<double>& r,
                      const unsigned& imacro,
                      const unsigned& part);

    /// Western boundary of the macro element imacro in the
    /// upstream (part=0) or downstream (part=1) sections
    void r_W_straight(const Vector<double>& zeta,
                      Vector<double>& r,
                      const unsigned& imacro,
                      const unsigned& part);

    /// Southern boundary of the macro element imacro in the
    ///  upstream (part=0) or downstream (part=1) sections
    void r_S_straight(const Vector<double>& zeta,
                      Vector<double>& r,
                      const unsigned& imacro,
                      const unsigned& part);

    /// Eastern boundary of the macro element imacro in the
    /// upstream (part=0) or downstream (part=1) sections
    void r_E_straight(const Vector<double>& zeta,
                      Vector<double>& r,
                      const unsigned& imacro,
                      const unsigned& part);

    /// Northern boundary of the macro element imacro in the collapsible
    /// section
    void r_N_collapsible(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& imacro);

    /// Western boundary of the macro element imacro in the collapsible
    /// section
    void r_W_collapsible(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& imacro);

    /// Southern boundary of the macro element imacro in the collapsible
    /// section
    void r_S_collapsible(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& imacro);

    /// Eastern boundary of the macro element imacro in the collapsible
    /// section
    void r_E_collapsible(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& imacro);


    /// Function pointer for function that squashes
    /// the macro elements near the walls
    BLSquashFctPt BL_squash_fct_pt;

    /// Function pointer for function that implements
    /// axial spacing of macro elements
    AxialSpacingFctPt Axial_spacing_fct_pt;

    /// Default for function that  implements
    /// axial spacing of macro elements
    static double default_axial_spacing_fct(const double& xi)
    {
      return xi;
    }


    /// Number of vertical element columns in upstream section
    unsigned Nup;

    /// Number of vertical element columns in "collapsible" section
    unsigned Ncollapsible;

    /// Number of vertical element columns in downstream section
    unsigned Ndown;

    /// Number of macro elements across channel
    unsigned Ny;

    /// x-length in the upstream part of the channel
    double Lup;

    /// x-length in the "collapsible" part of the channel
    double Lcollapsible;

    /// x-length in the downstream part of the channel
    double Ldown;

    /// Width
    double Ly;

    /// Pointer to the geometric object that parametrises the collapsible wall
    GeomObject* Wall_pt;

    /// Rotate domain (for axisymmetric problems, say)
    bool Rotate_domain;
  };


  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////

  //===================================================================
  /// Vector representation of the  imacro-th macro element
  /// boundary idirect (N/S/W/E) at time level t
  /// (t=0: present; t>0: previous): \f$ {\bf r}({\bf zeta}) \f$
  /// Note that the local coordinate \b zeta is a 1D
  /// Vector rather than a scalar -- this is unavoidable because
  /// this function implements the pure virtual function in the
  /// Domain base class.
  //=================================================================
  void CollapsibleChannelDomain::macro_element_boundary(
    const unsigned& t,
    const unsigned& imacro,
    const unsigned& idirect,
    const Vector<double>& zeta,
    Vector<double>& r)
  {
    using namespace QuadTreeNames;

#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "CollapsibleChannelDomain::macro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    // Total number of vertical columns of (macro-)elements
    unsigned n_x = Nup + Ncollapsible + Ndown;

    // Which direction?
    if (idirect == N)
    {
      // Upstream part of the channel
      if ((imacro % n_x) < Nup)
      {
        r_N_straight(zeta, r, imacro, 0);
      }
      // Downstream part of channel
      else if ((imacro % n_x) >= Nup + Ncollapsible)
      {
        r_N_straight(zeta, r, imacro, 1);
      }
      // Collapsible segment
      else if (((imacro % n_x) < Nup + Ncollapsible) && ((imacro % n_x) >= Nup))
      {
        r_N_collapsible(t, zeta, r, imacro);
      }
      else
      {
        std::ostringstream error_stream;
        error_stream << "Never get here! imacro, idirect: " << imacro << " "
                     << idirect << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
    else if (idirect == S)
    {
      // Upstream part
      if ((imacro % n_x) < Nup)
      {
        r_S_straight(zeta, r, imacro, 0);
      }
      // Downstream part
      else if ((imacro % n_x) >= Nup + Ncollapsible)
      {
        r_S_straight(zeta, r, imacro, 1);
      }
      // "Collapsible" bit
      else if (((imacro % n_x) < Nup + Ncollapsible) && ((imacro % n_x) >= Nup))
      {
        r_S_collapsible(t, zeta, r, imacro);
      }
      else
      {
        std::ostringstream error_stream;
        error_stream << "Never get here! imacro, idirect: " << imacro << " "
                     << idirect << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
    else if (idirect == E)
    {
      // Upstream bit
      if ((imacro % n_x) < Nup)
      {
        r_E_straight(zeta, r, imacro, 0);
      }
      // Downstream bit
      else if ((imacro % n_x) >= Nup + Ncollapsible)
      {
        r_E_straight(zeta, r, imacro, 1);
      }
      // "Collapsible" bit
      else if (((imacro % n_x) < Nup + Ncollapsible) && ((imacro % n_x) >= Nup))
      {
        r_E_collapsible(t, zeta, r, imacro);
      }
      else
      {
        std::ostringstream error_stream;
        error_stream << "Never get here! imacro, idirect: " << imacro << " "
                     << idirect << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    else if (idirect == W)
    {
      // Upstream bit
      if ((imacro % n_x) < Nup)
      {
        r_W_straight(zeta, r, imacro, 0);
      }
      // Downstream bit
      else if ((imacro % n_x) >= Nup + Ncollapsible)
      {
        r_W_straight(zeta, r, imacro, 1);
      }
      // "Collapsible" bit
      else if (((imacro % n_x) < Nup + Ncollapsible) && ((imacro % n_x) >= Nup))
      {
        r_W_collapsible(t, zeta, r, imacro);
      }
      else
      {
        std::ostringstream error_stream;
        error_stream << "Never get here! imacro, idirect: " << imacro << " "
                     << idirect << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Rotate?
    if (Rotate_domain)
    {
      double radius = r[1];
      double z = r[0];

      r[0] = radius;
      r[1] = -z;
    }
  }


  //===========================================================================
  /// Western edge of the  macro element in the upstream (part=0)
  /// or downstream (part=1) parts of the channel; \f$ \zeta \in [-1,1] \f$
  //===========================================================================
  void CollapsibleChannelDomain::r_W_straight(const Vector<double>& zeta,
                                              Vector<double>& r,
                                              const unsigned& imacro,
                                              const unsigned& part)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Where are we?
    switch (part)
    {
      case 0: // in the upstream part of the channel

        // Parametrize the boundary
        r[0] = axial_spacing_fct(double(x) * (Lup / double(Nup)));
        r[1] = (double(y) + (0.5 * (1.0 + zeta[0]))) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      case 1: // in the downstream part of the channel

        // Parametrizes the boundary
        r[0] = axial_spacing_fct(double(x - Nup - Ncollapsible) *
                                   (Ldown / double(Ndown)) +
                                 Lup + Lcollapsible);
        r[1] = (double(y) + (0.5 * (1.0 + zeta[0]))) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      default:

        std::ostringstream error_stream;
        error_stream << "Never get here! part=" << part << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //===========================================================================
  /// Eastern  edge of the  macro element in the straight parts
  /// of the channel; \f$ \zeta \in [-1,1] \f$
  /// part=0 in the upstream part, part=1 in the downstream part
  //===========================================================================
  void CollapsibleChannelDomain::r_E_straight(const Vector<double>& zeta,
                                              Vector<double>& r,
                                              const unsigned& imacro,
                                              const unsigned& part)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Where are we?
    switch (part)
    {
      case 0: // in the upstream part of the channel

        // Parametrizes the boundary
        r[0] = axial_spacing_fct((double(x) + 1.0) * (Lup / double(Nup)));
        r[1] = (double(y) + (0.5 * (1.0 + zeta[0]))) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      case 1: // in the downstream part of the channel

        // Parametrizes the boundary
        r[0] = axial_spacing_fct((double(x - Nup - Ncollapsible) + 1.0) *
                                   (Ldown / double(Ndown)) +
                                 Lup + Lcollapsible);
        r[1] = (double(y) + (0.5 * (1.0 + zeta[0]))) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      default:


        std::ostringstream error_stream;
        error_stream << "Never get here! part=" << part << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  //==========================================================================
  /// Northern edge of the  macro element in the straight parts of
  /// the channel; \f$ \zeta \in [-1,1] \f$
  /// part=0 in the left part, part=1 in the right part
  //==========================================================================
  void CollapsibleChannelDomain::r_N_straight(const Vector<double>& zeta,
                                              Vector<double>& r,
                                              const unsigned& imacro,
                                              const unsigned& part)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));


    // Where are we?
    switch (part)
    {
      case 0: // in the upstream part of the channel

        // Parametrizes the boundary
        r[0] = axial_spacing_fct((double(x) + (0.5 * (1.0 + zeta[0]))) *
                                 (Lup / double(Nup)));
        r[1] = (double(y) + 1.0) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      case 1: // in the downstream part of the channel

        // Parametrizes the boundary
        r[0] = axial_spacing_fct(
          (double(x - Nup - Ncollapsible) + (0.5 * (1.0 + zeta[0]))) *
            (Ldown / double(Ndown)) +
          Lup + Lcollapsible);
        r[1] = (double(y) + 1.0) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      default:


        std::ostringstream error_stream;
        error_stream << "Never get here! part=" << part << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //=========================================================================
  /// Southern  edge of the  macro element in the straight parts of
  /// the channel; \f$ \zeta \in [-1,1] \f$
  /// part=0 in the left part, part=1 in the right part
  //=========================================================================
  void CollapsibleChannelDomain::r_S_straight(const Vector<double>& zeta,
                                              Vector<double>& r,
                                              const unsigned& imacro,
                                              const unsigned& part)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Where are we?
    switch (part)
    {
      case 0: // in the upstream bit

        // Parametrizes the boundary
        r[0] = axial_spacing_fct((double(x) + (0.5 * (1 + zeta[0]))) *
                                 (Lup / double(Nup)));
        r[1] = double(y) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      case 1: // in the downstream bit

        // Parametrizes the boundary
        r[0] = axial_spacing_fct(
          (double(x - Nup - Ncollapsible) + (0.5 * (1 + zeta[0]))) *
            (Ldown / double(Ndown)) +
          Lup + Lcollapsible);
        r[1] = double(y) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      default:


        std::ostringstream error_stream;
        error_stream << "Never get here! part=" << part << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //========================================================================
  /// Western edge of the  macro element in the collapsible part of the
  /// channel; \f$ \zeta \in [-1,1] \f$.
  //========================================================================
  void CollapsibleChannelDomain::r_W_collapsible(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& imacro)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Vector of Lagrangian coordinates
    Vector<double> xi(1);
    xi[0] = double(x - Nup) * (Lcollapsible / double(Ncollapsible));

    // Position vector on upper wall:
    Vector<double> r_wall(2);
    Wall_pt->position(t, xi, r_wall);

    // Point will be located on straight line from bottom to top wall
    double fract = (double(y) + (0.5 * (1.0 + zeta[0]))) / double(Ny);

    // Map it via squash fct
    fract = s_squash(fract);

    // x-cooordinate -- straight line from fixed position on the bottom
    // wall to moving position on the top wall
    r[0] = Lup + xi[0] + (r_wall[0] - (xi[0] + Lup)) * fract;

    // y-coordinate
    r[1] = r_wall[1] * fract;
  }


  //=========================================================================
  /// Eastern  edge of the  macro element in the collapsible part of the
  /// channel; \f$ \zeta \in [-1,1] \f$
  //=========================================================================
  void CollapsibleChannelDomain::r_E_collapsible(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& imacro)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Vector of Lagrangian coordinates
    Vector<double> xi(1);
    xi[0] = (double(x - Nup) + 1.0) * (Lcollapsible / double(Ncollapsible));

    // Position vector on upper wall:
    Vector<double> r_wall(2);
    Wall_pt->position(t, xi, r_wall);

    // Point will be located on straight line from bottom to top wall
    double fract = (double(y) + (0.5 * (1.0 + zeta[0]))) / double(Ny);

    // Map it via squash fct
    fract = s_squash(fract);

    // x-cooordinate -- straight line from fixed position on the bottom
    // wall to moving position on the top wall
    r[0] = Lup + xi[0] + (r_wall[0] - (xi[0] + Lup)) * fract;

    // y-coordinate
    r[1] = r_wall[1] * fract;
  }


  //==========================================================================
  /// Northern edge of the  macro element in the collapsible part of the
  /// channel; \f$ \zeta \in [-1,1] \f$
  //==========================================================================
  void CollapsibleChannelDomain::r_N_collapsible(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& imacro)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Vector of Lagrangian coordinates
    Vector<double> xi(1);
    xi[0] = (double(x - Nup) + (0.5 * (1.0 + zeta[0]))) *
            (Lcollapsible / double(Ncollapsible));

    // Position vector on upper wall:
    Vector<double> r_wall(2);
    Wall_pt->position(t, xi, r_wall);

    // Point will be located on straight line from bottom to top wall
    double fract = (double(y) + 1.0) / double(Ny);

    // Map it via squash fct
    fract = s_squash(fract);

    // x-cooordinate -- straight line from fixed position on the bottom
    // wall to moving position on the top wall
    r[0] = Lup + xi[0] + (r_wall[0] - (xi[0] + Lup)) * fract;

    // y-coordinate
    r[1] = r_wall[1] * fract;
  }


  //========================================================================
  /// Southern  edge of the  macro element in the collapsible part of the
  /// channel; \f$ \zeta \in [-1,1] \f$
  //========================================================================
  void CollapsibleChannelDomain::r_S_collapsible(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& imacro)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Vector of Lagrangian coordinates
    Vector<double> xi(1);
    xi[0] = (double(x - Nup) + (0.5 * (1.0 + zeta[0]))) *
            (Lcollapsible / double(Ncollapsible));

    // Position vector on upper wall:
    Vector<double> r_wall(2);
    Wall_pt->position(t, xi, r_wall);

    // Point will be located on straight line from bottom to top wall
    double fract = double(y) / double(Ny);

    // Map it via squash fct
    fract = s_squash(fract);

    // x-cooordinate -- straight line from fixed position on the bottom
    // wall to moving position on the top wall
    r[0] = Lup + xi[0] + (r_wall[0] - (xi[0] + Lup)) * fract;

    // y-coordinate
    r[1] = r_wall[1] * fract;
  }


} // namespace oomph

#endif
