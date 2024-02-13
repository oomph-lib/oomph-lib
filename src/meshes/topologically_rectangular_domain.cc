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
#include "topologically_rectangular_domain.h"


namespace oomph
{
  //=============================================================================
  /// Constructor - domain boundaries are described with four boundary
  /// function pointers describing the topology of the north, east, south, and
  /// west boundaries
  //=============================================================================
  TopologicallyRectangularDomain::TopologicallyRectangularDomain(
    BoundaryFctPt north_pt,
    BoundaryFctPt east_pt,
    BoundaryFctPt south_pt,
    BoundaryFctPt west_pt)
  {
    // domain comprises one macro element
    Macro_element_pt.resize(1);

    // Create the macro element
    Macro_element_pt[0] = new QMacroElement<2>(this, 0);

    // set boundary function pointers
    North_boundary_fn_pt = north_pt;
    East_boundary_fn_pt = east_pt;
    South_boundary_fn_pt = south_pt;
    West_boundary_fn_pt = west_pt;

    // by default derivative boundary function pointers are null
    dNorth_boundary_fn_pt = 0;
    dEast_boundary_fn_pt = 0;
    dSouth_boundary_fn_pt = 0;
    dWest_boundary_fn_pt = 0;

    // also by default second derivate boundary function pointers are null
    d2North_boundary_fn_pt = 0;
    d2East_boundary_fn_pt = 0;
    d2South_boundary_fn_pt = 0;
    d2West_boundary_fn_pt = 0;

    // paranoid self check to ensure that ends of boundaries meet
#ifdef PARANOID
    // error message stream
    std::ostringstream error_message;
    bool error_flag = false;
    // check NE corner
    {
      Vector<double> x_N(2);
      (*North_boundary_fn_pt)(1.0, x_N);
      Vector<double> x_E(2);
      (*East_boundary_fn_pt)(1.0, x_E);
      if (x_N[0] != x_E[0] || x_N[1] != x_E[1])
      {
        error_message << "North and East Boundaries do not meet at the "
                      << "North East Corner.\n"
                      << "North Boundary : x[0] = " << x_N[0] << "\n"
                      << "                 x[1] = " << x_N[1] << "\n"
                      << "East Boundary : x[0] = " << x_E[0] << "\n"
                      << "                x[1] = " << x_E[1] << "\n\n";
        error_flag = true;
      }
    }
    // check SE corner
    {
      Vector<double> x_S(2);
      (*South_boundary_fn_pt)(1.0, x_S);
      Vector<double> x_E(2);
      (*East_boundary_fn_pt)(-1.0, x_E);
      if (x_S[0] != x_E[0] || x_S[1] != x_E[1])
      {
        error_message << "South and East Boundaries do not meet at the "
                      << "South East Corner.\n"
                      << "South Boundary : x[0] = " << x_S[0] << "\n"
                      << "                 x[1] = " << x_S[1] << "\n"
                      << "East Boundary : x[0] = " << x_E[0] << "\n"
                      << "                x[1] = " << x_E[1] << "\n\n";
        error_flag = true;
      }
    }
    // check SW corner
    {
      Vector<double> x_S(2);
      (*South_boundary_fn_pt)(-1.0, x_S);
      Vector<double> x_W(2);
      (*West_boundary_fn_pt)(-1.0, x_W);
      if (x_S[0] != x_W[0] || x_S[1] != x_W[1])
      {
        error_message << "South and West Boundaries do not meet at the "
                      << "South West Corner.\n"
                      << "South Boundary : x[0] = " << x_S[0] << "\n"
                      << "                 x[1] = " << x_S[1] << "\n"
                      << "West Boundary : x[0] = " << x_W[0] << "\n"
                      << "                x[1] = " << x_W[1] << "\n\n";
        error_flag = true;
      }
    }
    // check NW corner
    {
      Vector<double> x_N(2);
      (*North_boundary_fn_pt)(-1.0, x_N);
      Vector<double> x_W(2);
      (*West_boundary_fn_pt)(1.0, x_W);
      if (x_N[0] != x_W[0] || x_N[1] != x_W[1])
      {
        error_message << "North and West Boundaries do not meet at the "
                      << "North West Corner.\n"
                      << "North Boundary : x[0] = " << x_N[0] << "\n"
                      << "                 x[1] = " << x_N[1] << "\n"
                      << "West Boundary : x[0] = " << x_W[0] << "\n"
                      << "                x[1] = " << x_W[1] << "\n\n";
        error_flag = true;
      }
    }
    if (error_flag)
    {
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif
  }

  //=============================================================================
  /// Constructor - takes length of domain in x and y direction as
  /// arguements. Assumes domain is rectangular, and the south west (lower
  /// left) corner is at 0,0.
  //=============================================================================
  TopologicallyRectangularDomain::TopologicallyRectangularDomain(
    const double& l_x, const double& l_y)
  {
    // domain comprises one macro element
    Macro_element_pt.resize(1);

    // Create the macro element
    Macro_element_pt[0] = new QMacroElement<2>(this, 0);

    // set the boundary function pointers to zero
    North_boundary_fn_pt = 0;
    East_boundary_fn_pt = 0;
    South_boundary_fn_pt = 0;
    West_boundary_fn_pt = 0;

    // resize x vectors
    x_south_west.resize(2);
    x_north_east.resize(2);

    // set south west corner to 0,0
    x_south_west[0] = 0.0;
    x_south_west[1] = 0.0;

    // set north east corner
    x_north_east[0] = l_x;
    x_north_east[1] = l_y;

    // by default derivative boundary function pointers are null
    dNorth_boundary_fn_pt = 0;
    dEast_boundary_fn_pt = 0;
    dSouth_boundary_fn_pt = 0;
    dWest_boundary_fn_pt = 0;

    // also by default second derivate boundary function pointers are null
    d2North_boundary_fn_pt = 0;
    d2East_boundary_fn_pt = 0;
    d2South_boundary_fn_pt = 0;
    d2West_boundary_fn_pt = 0;
  }

  //=============================================================================
  /// Constructor - takes the minimum and maximum coordinates of the
  /// of an assumed rectanguler domain in the x and y direction
  //=============================================================================
  TopologicallyRectangularDomain::TopologicallyRectangularDomain(
    const double& x_min,
    const double& x_max,
    const double& y_min,
    const double& y_max)
  {
    // domain comprises one macro element
    Macro_element_pt.resize(1);

    // Create the macro element
    Macro_element_pt[0] = new QMacroElement<2>(this, 0);

    // set the boundary function pointers to zero
    North_boundary_fn_pt = 0;
    East_boundary_fn_pt = 0;
    South_boundary_fn_pt = 0;
    West_boundary_fn_pt = 0;

    // resize x vectors
    x_south_west.resize(2);
    x_north_east.resize(2);

    // set vector values
    x_south_west[0] = x_min;
    x_south_west[1] = y_min;
    x_north_east[0] = x_max;
    x_north_east[1] = y_max;

    // by default derivative boundary function pointers are null
    dNorth_boundary_fn_pt = 0;
    dEast_boundary_fn_pt = 0;
    dSouth_boundary_fn_pt = 0;
    dWest_boundary_fn_pt = 0;

    // also by default second derivate boundary function pointers are null
    d2North_boundary_fn_pt = 0;
    d2East_boundary_fn_pt = 0;
    d2South_boundary_fn_pt = 0;
    d2West_boundary_fn_pt = 0;
  }

  //=============================================================================
  /// allows the boundary derivate function pointers to be set. To
  /// compute the derivatives of the problem domain global coordinates (x_i) wrt
  /// the macro element coordinates (m_i), dx_i/dm_t is required along the
  /// domain boundaries (where dm_t is the macro element coordinate tangential
  /// to the domain boundary). The derivatives dx_i/dm_t can either be
  /// prescribed with function pointers, or if the function pointers are not
  /// provided then dx_i/dm_t is computed with finite differencing.
  /// Note - these functions are only required for domains contructed with
  /// boundary function pointers
  //=============================================================================
  void TopologicallyRectangularDomain::set_boundary_derivative_functions(
    BoundaryFctPt d_north_pt,
    BoundaryFctPt d_east_pt,
    BoundaryFctPt d_south_pt,
    BoundaryFctPt d_west_pt)
  {
    // set the boundary derivate function pointers
    dNorth_boundary_fn_pt = d_north_pt;
    dEast_boundary_fn_pt = d_east_pt;
    dSouth_boundary_fn_pt = d_south_pt;
    dWest_boundary_fn_pt = d_west_pt;
  }


  //=============================================================================
  /// allows the boundary second derivate function pointers to be set.
  /// To compute the second derivatives of the problem domain global
  /// coordinates (x_i) wrt the macro element coordinates (m_i), d2x_i/dm_t^2
  /// is required along the domain boundaries (where dm_t is the macro element
  /// coordinate tangential to the domain boundary). The derivatives
  /// d2x_i/dm_t^2 can either be prescribed with function pointers, or if the
  /// function pointers are not provided then dx_i/dm_t is computed with finite
  /// differencing.
  /// Note - these functions are only required for domains contructed with
  /// boundary function pointers
  //=============================================================================
  void TopologicallyRectangularDomain::set_boundary_second_derivative_functions(
    BoundaryFctPt d2_north_pt,
    BoundaryFctPt d2_east_pt,
    BoundaryFctPt d2_south_pt,
    BoundaryFctPt d2_west_pt)
  {
    // set the boundary derivate function pointers
    d2North_boundary_fn_pt = d2_north_pt;
    d2East_boundary_fn_pt = d2_east_pt;
    d2South_boundary_fn_pt = d2_south_pt;
    d2West_boundary_fn_pt = d2_west_pt;
  }


  //=============================================================================
  /// returns the global coordinate position (f) of macro element position s
  /// on boundary i_direct (e.g. N/S/W/E in 2D) at time t (no time dependence)
  //=============================================================================
  void TopologicallyRectangularDomain::macro_element_boundary(
    const unsigned& t,
    const unsigned& i_macro,
    const unsigned& i_direct,
    const Vector<double>& s,
    Vector<double>& f)
  {
    // use quad tree edge names to label edge of domain
    using namespace QuadTreeNames;

#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "TopologicallyRectangularDomain::macro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    // north boundary
    if (i_direct == N) r_N(s, f);
    // east boundary
    else if (i_direct == E)
      r_E(s, f);
    // south boundary
    else if (i_direct == S)
      r_S(s, f);
    // west boundary
    else if (i_direct == W)
      r_W(s, f);
  }

  //=============================================================================
  /// returns the derivates of the global coordinate position (f) wrt to the
  /// macro element coordinate at macro macro element position s on boundary
  /// i_direct (e.g. N/S/W/E in 2D) at time t (no time dependence)
  //=============================================================================
  void TopologicallyRectangularDomain::dmacro_element_boundary(
    const unsigned& t,
    const unsigned& i_macro,
    const unsigned& i_direct,
    const Vector<double>& s,
    Vector<double>& f)
  {
    // use quad tree edge names to label edge of domain
    using namespace QuadTreeNames;

#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "TopologicallyRectangularDomain::dmacro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    // north boundary
    if (i_direct == N) dr_N(s, f);
    // east boundary
    else if (i_direct == E)
      dr_E(s, f);
    // south boundary
    else if (i_direct == S)
      dr_S(s, f);
    // west boundary
    else if (i_direct == W)
      dr_W(s, f);
  }

  //=============================================================================
  /// returns the second derivates of the global coordinate position (f) wrt to
  /// the macro element coordinate at macro macro element position s on boundary
  /// i_direct (e.g. N/S/W/E in 2D) at time t (no time dependence)
  //=============================================================================
  void TopologicallyRectangularDomain::d2macro_element_boundary(
    const unsigned& t,
    const unsigned& i_macro,
    const unsigned& i_direct,
    const Vector<double>& s,
    Vector<double>& f)
  {
    // use quad tree edge names to label edge of domain
    using namespace QuadTreeNames;


#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "TopologicallyRectangularDomain::d2macro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    // north boundary
    if (i_direct == N) d2r_N(s, f);
    // east boundary
    else if (i_direct == E)
      d2r_E(s, f);
    // south boundary
    else if (i_direct == S)
      d2r_S(s, f);
    // west boundary
    else if (i_direct == W)
      d2r_W(s, f);
  }

  //=============================================================================
  /// takes the macro element coordinate position along the north
  /// boundary and returns the global coordinate position along that boundary
  //=============================================================================
  void TopologicallyRectangularDomain::r_N(const Vector<double>& s,
                                           Vector<double>& f)
  {
    if (North_boundary_fn_pt != 0)
    {
      (*North_boundary_fn_pt)(s[0], f);
    }
    else
    {
      f[0] =
        x_south_west[0] + (s[0] + 1) / 2 * (x_north_east[0] - x_south_west[0]);
      f[1] = x_north_east[1];
    }
  }

  //=============================================================================
  /// takes the macro element coordinate position along the east
  /// boundary and returns the global coordinate position along that boundary
  //=============================================================================
  void TopologicallyRectangularDomain::r_E(const Vector<double>& s,
                                           Vector<double>& f)
  {
    if (East_boundary_fn_pt != 0)
    {
      (*East_boundary_fn_pt)(s[0], f);
    }
    else
    {
      f[0] = x_north_east[0];
      f[1] =
        x_south_west[1] + (s[0] + 1) / 2 * (x_north_east[1] - x_south_west[1]);
    }
  }

  //=============================================================================
  /// takes the macro element coordinate position along the south
  /// boundary and returns the global coordinate position along that boundary
  //=============================================================================
  void TopologicallyRectangularDomain::r_S(const Vector<double>& s,
                                           Vector<double>& f)
  {
    if (South_boundary_fn_pt != 0)
    {
      (*South_boundary_fn_pt)(s[0], f);
    }
    else
    {
      f[0] =
        x_south_west[0] + (s[0] + 1) / 2 * (x_north_east[0] - x_south_west[0]);
      f[1] = x_south_west[1];
    }
  }

  //=============================================================================
  /// takes the macro element coordinate position along the west
  /// boundary and returns the global coordinate position along that boundary
  /// access down boundary function pointer
  //=============================================================================
  void TopologicallyRectangularDomain::r_W(const Vector<double>& s,
                                           Vector<double>& f)
  {
    if (West_boundary_fn_pt != 0)
    {
      (*West_boundary_fn_pt)(s[0], f);
    }
    else
    {
      f[0] = x_south_west[0];
      f[1] =
        x_south_west[1] + (s[0] + 1) / 2 * (x_north_east[1] - x_south_west[1]);
    }
  }

  //=============================================================================
  /// takes the macro element coordinate position along the north
  /// boundary and returns the derivates of the global coordinates with respect
  /// to the boundary
  //=============================================================================
  void TopologicallyRectangularDomain::dr_N(const Vector<double>& s,
                                            Vector<double>& dr)
  {
    // if N boundary fn provided
    if (North_boundary_fn_pt != 0)
    {
      // if dN boundary fn provided
      if (dNorth_boundary_fn_pt != 0)
      {
        (*dNorth_boundary_fn_pt)(s[0], dr);
      }
      // else compute by finite differencing
      else
      {
        const double h = 10e-8;
        Vector<double> x_N_left(2);
        (*North_boundary_fn_pt)(s[0] - 0.5 * h, x_N_left);
        Vector<double> x_N_right(2);
        (*North_boundary_fn_pt)(s[0] + 0.5 * h, x_N_right);
        dr[0] = (x_N_right[0] - x_N_left[0]) / h;
        dr[1] = (x_N_right[1] - x_N_left[1]) / h;
      }
    }
    // if N boundary fn not provided then mesh is rectangular
    else
    {
      dr[0] = (x_north_east[0] - x_south_west[0]) / 2;
      dr[1] = 0;
    }
  }

  //=============================================================================
  /// takes the macro element coordinate position along the E
  /// boundary and returns the derivates of the global coordinates with respect
  /// to the boundary
  //=============================================================================
  void TopologicallyRectangularDomain::dr_E(const Vector<double>& s,
                                            Vector<double>& dr)
  {
    // if E boundary fn provided
    if (East_boundary_fn_pt != 0)
    {
      // if dE boundary fn provided
      if (dEast_boundary_fn_pt != 0)
      {
        (*dEast_boundary_fn_pt)(s[0], dr);
      }
      // else compute by finite differencing
      else
      {
        const double h = 10e-8;
        Vector<double> x_E_down(2);
        (*East_boundary_fn_pt)(s[0] - 0.5 * h, x_E_down);
        Vector<double> x_E_up(2);
        (*East_boundary_fn_pt)(s[0] + 0.5 * h, x_E_up);
        dr[0] = (x_E_up[0] - x_E_down[0]) / h;
        dr[1] = (x_E_up[1] - x_E_down[1]) / h;
      }
    }
    // if E boundary fn not provided then mesh is rectangular
    else
    {
      dr[0] = 0;
      dr[1] = (x_north_east[1] - x_south_west[1]) / 2;
    }
  }

  //=============================================================================
  /// takes the macro element coordinate position along the south
  /// boundary and returns the derivates of the global coordinates with respect
  /// to the boundary
  //=============================================================================
  void TopologicallyRectangularDomain::dr_S(const Vector<double>& s,
                                            Vector<double>& dr)
  {
    // if S boundary fn provided
    if (South_boundary_fn_pt != 0)
    {
      // if dS boundary fn provided
      if (dSouth_boundary_fn_pt != 0)
      {
        (*dSouth_boundary_fn_pt)(s[0], dr);
      }
      // else compute by finite differencing
      else
      {
        const double h = 10e-8;
        Vector<double> x_N_left(2);
        (*South_boundary_fn_pt)(s[0] - 0.5 * h, x_N_left);
        Vector<double> x_N_right(2);
        (*South_boundary_fn_pt)(s[0] + 0.5 * h, x_N_right);
        dr[0] = (x_N_right[0] - x_N_left[0]) / h;
        dr[1] = (x_N_right[1] - x_N_left[1]) / h;
      }
    }
    // if S boundary fn not provided then mesh is rectangular
    else
    {
      dr[0] = (x_north_east[0] - x_south_west[0]) / 2;
      dr[1] = 0;
    }
  }

  //=============================================================================
  /// takes the macro element coordinate position along the W
  /// boundary and returns the derivates of the global coordinates with respect
  /// to the boundary
  //=============================================================================
  void TopologicallyRectangularDomain::dr_W(const Vector<double>& s,
                                            Vector<double>& dr)
  {
    // if W boundary fn provided
    if (West_boundary_fn_pt != 0)
    {
      // if dW boundary fn provided
      if (dWest_boundary_fn_pt != 0)
      {
        (*dWest_boundary_fn_pt)(s[0], dr);
      }
      // else compute by finite differencing
      else
      {
        const double h = 10e-8;
        Vector<double> x_W_down(2);
        (*West_boundary_fn_pt)(s[0] - 0.5 * h, x_W_down);
        Vector<double> x_W_up(2);
        (*West_boundary_fn_pt)(s[0] + 0.5 * h, x_W_up);
        dr[0] = (x_W_up[0] - x_W_down[0]) / h;
        dr[1] = (x_W_up[1] - x_W_down[1]) / h;
      }
    }
    // if E boundary fn not provided then mesh is rectangular
    else
    {
      dr[0] = 0;
      dr[1] = (x_north_east[1] - x_south_west[1]) / 2;
    }
  }

  //=============================================================================
  /// takes the macro element coordinate position along the north
  /// boundary and returns the second derivates of the global coordinates with
  /// respect to the boundary
  //=============================================================================
  void TopologicallyRectangularDomain::d2r_N(const Vector<double>& s,
                                             Vector<double>& d2r)
  {
    // if N boundary fn provided
    if (North_boundary_fn_pt != 0)
    {
      // if d2N boundary fn provided
      if (d2North_boundary_fn_pt != 0)
      {
        (*d2North_boundary_fn_pt)(s[0], d2r);
      }
      // else compute by finite differencing
      else
      {
        // if dN boundary fn provided finite difference d2N from it
        if (dNorth_boundary_fn_pt != 0)
        {
          const double h = 10e-8;
          Vector<double> dx_N_left(2);
          (*dNorth_boundary_fn_pt)(s[0] - 0.5 * h, dx_N_left);
          Vector<double> dx_N_right(2);
          (*dNorth_boundary_fn_pt)(s[0] + 0.5 * h, dx_N_right);
          d2r[0] = (dx_N_right[0] - dx_N_left[0]) / h;
          d2r[1] = (dx_N_right[1] - dx_N_left[1]) / h;
        }
        // else finite difference from N boundary fn
        else
        {
          const double h = 10e-8;
          Vector<double> N_left(2);
          (*North_boundary_fn_pt)(s[0] - h, N_left);
          Vector<double> N_right(2);
          (*North_boundary_fn_pt)(s[0] + h, N_right);
          Vector<double> N_centre(2);
          (*North_boundary_fn_pt)(s[0], N_centre);
          d2r[0] = (N_right[0] + N_left[0] - 2 * N_centre[0]) / (h * h);
          d2r[1] = (N_right[1] + N_left[1] - 2 * N_centre[1]) / (h * h);
        }
      }
    }
    // if N boundary fn not provided then mesh is rectangular
    else
    {
      d2r[0] = 0;
      d2r[1] = 0;
    }
  }

  //=============================================================================
  /// takes the macro element coordinate position along the east
  /// boundary and returns the second derivates of the global coordinates with
  /// respect to the boundary
  //=============================================================================
  void TopologicallyRectangularDomain::d2r_E(const Vector<double>& s,
                                             Vector<double>& d2r)
  {
    // if E boundary fn provided
    if (East_boundary_fn_pt != 0)
    {
      // if d2E boundary fn provided
      if (d2East_boundary_fn_pt != 0)
      {
        (*d2East_boundary_fn_pt)(s[0], d2r);
      }
      // else compute by finite differencing
      else
      {
        // if dE boundary fn provided finite difference d2E from it
        if (dEast_boundary_fn_pt != 0)
        {
          const double h = 10e-8;
          Vector<double> dx_E_lower(2);
          (*dEast_boundary_fn_pt)(s[0] - 0.5 * h, dx_E_lower);
          Vector<double> dx_E_upper(2);
          (*dEast_boundary_fn_pt)(s[0] + 0.5 * h, dx_E_upper);
          d2r[0] = (dx_E_upper[0] - dx_E_lower[0]) / h;
          d2r[1] = (dx_E_upper[1] - dx_E_lower[1]) / h;
        }
        // else finite difference from E boundary fn
        else
        {
          const double h = 10e-8;
          Vector<double> E_left(2);
          (*East_boundary_fn_pt)(s[0] - h, E_left);
          Vector<double> E_right(2);
          (*East_boundary_fn_pt)(s[0] + h, E_right);
          Vector<double> E_centre(2);
          (*East_boundary_fn_pt)(s[0], E_centre);
          d2r[0] = (E_right[0] + E_left[0] - 2 * E_centre[0]) / (h * h);
          d2r[1] = (E_right[1] + E_left[1] - 2 * E_centre[1]) / (h * h);
        }
      }
    }
    // if E boundary fn not provided then mesh is rectangular
    else
    {
      d2r[0] = 0;
      d2r[1] = 0;
    }
  }

  //=============================================================================
  /// takes the macro element coordinate position along the south
  /// boundary and returns the second derivates of the global coordinates with
  /// respect to the boundary
  //=============================================================================
  void TopologicallyRectangularDomain::d2r_S(const Vector<double>& s,
                                             Vector<double>& d2r)
  {
    // if S boundary fn provided
    if (South_boundary_fn_pt != 0)
    {
      // if d2S boundary fn provided
      if (d2South_boundary_fn_pt != 0)
      {
        (*d2South_boundary_fn_pt)(s[0], d2r);
      }
      // else compute by finite differencing
      else
      {
        // if dS boundary fn provided finite difference d2S from it
        if (dSouth_boundary_fn_pt != 0)
        {
          const double h = 10e-8;
          Vector<double> dx_S_left(2);
          (*dSouth_boundary_fn_pt)(s[0] - 0.5 * h, dx_S_left);
          Vector<double> dx_S_right(2);
          (*dSouth_boundary_fn_pt)(s[0] + 0.5 * h, dx_S_right);
          d2r[0] = (dx_S_right[0] - dx_S_left[0]) / h;
          d2r[1] = (dx_S_right[1] - dx_S_left[1]) / h;
        }
        // else finite difference from S boundary fn
        else
        {
          const double h = 10e-8;
          Vector<double> S_left(2);
          (*South_boundary_fn_pt)(s[0] - h, S_left);
          Vector<double> S_right(2);
          (*South_boundary_fn_pt)(s[0] + h, S_right);
          Vector<double> S_centre(2);
          (*South_boundary_fn_pt)(s[0], S_centre);
          d2r[0] = (S_right[0] + S_left[0] - 2 * S_centre[0]) / (h * h);
          d2r[1] = (S_right[1] + S_left[1] - 2 * S_centre[1]) / (h * h);
        }
      }
    }
    // if S boundary fn not provided then mesh is rectangular
    else
    {
      d2r[0] = 0;
      d2r[1] = 0;
    }
  }

  //=============================================================================
  /// takes the macro element coordinate position along the west
  /// boundary and returns the second derivates of the global coordinates with
  /// respect to the boundary
  //=============================================================================
  void TopologicallyRectangularDomain::d2r_W(const Vector<double>& s,
                                             Vector<double>& d2r)
  {
    // if W boundary fn provided
    if (West_boundary_fn_pt != 0)
    {
      // if d2W boundary fn provided
      if (d2West_boundary_fn_pt != 0)
      {
        (*d2West_boundary_fn_pt)(s[0], d2r);
      }
      // else compute by finite differencing
      else
      {
        // if dW boundary fn provided finite difference d2W from it
        if (dWest_boundary_fn_pt != 0)
        {
          const double h = 10e-8;
          Vector<double> dx_W_lower(2);
          (*dWest_boundary_fn_pt)(s[0] - 0.5 * h, dx_W_lower);
          Vector<double> dx_W_upper(2);
          (*dWest_boundary_fn_pt)(s[0] + 0.5 * h, dx_W_upper);
          d2r[0] = (dx_W_upper[0] - dx_W_lower[0]) / h;
          d2r[1] = (dx_W_upper[1] - dx_W_lower[1]) / h;
        }
        // else finite difference from W boundary fn
        else
        {
          const double h = 10e-8;
          Vector<double> W_left(2);
          (*West_boundary_fn_pt)(s[0] - h, W_left);
          Vector<double> W_right(2);
          (*West_boundary_fn_pt)(s[0] + h, W_right);
          Vector<double> W_centre(2);
          (*West_boundary_fn_pt)(s[0], W_centre);
          d2r[0] = (W_right[0] + W_left[0] - 2 * W_centre[0]) / (h * h);
          d2r[1] = (W_right[1] + W_left[1] - 2 * W_centre[1]) / (h * h);
        }
      }
    }
    // if W boundary fn not provided then mesh is rectangular
    else
    {
      d2r[0] = 0;
      d2r[1] = 0;
    }
  }
} // namespace oomph
