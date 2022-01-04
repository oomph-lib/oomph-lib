// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
#include "macro_element.h"
#include "domain.h"

namespace oomph
{
  //=================================================================
  /// Get global position r(S) at discrete time level t.
  /// t=0: Present time; t>0: previous timestep.
  //=================================================================
  void QMacroElement<2>::macro_map(const double& t,
                                   const Vector<double>& s,
                                   Vector<double>& r)
  {
    using namespace QuadTreeNames;

    Vector<double> bound_N(2);
    Vector<double> bound_S(2);
    Vector<double> bound_W(2);
    Vector<double> bound_E(2);

    Vector<double> diff_N(2);
    Vector<double> diff_S(2);
    Vector<double> diff_W(2);
    Vector<double> diff_E(2);

    Vector<double> f_rect(2);

    Vector<double> corner_SE(2);
    Vector<double> corner_SW(2);
    Vector<double> corner_NE(2);
    Vector<double> corner_NW(2);

    Vector<double> edge_N(2);
    Vector<double> edge_S(2);

    Vector<double> zeta(1);

    // Determine Vectors to corners
    zeta[0] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, corner_SE);
    zeta[0] = -1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, corner_SW);
    zeta[0] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, corner_NE);
    zeta[0] = -1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, corner_NW);


    // Get the position on the N/S/W/E edges
    zeta[0] = s[0];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, bound_N);
    zeta[0] = s[0];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, bound_S);
    zeta[0] = s[1];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::W, zeta, bound_W);
    zeta[0] = s[1];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::E, zeta, bound_E);


    for (int i = 0; i < 2; i++)
    {
      // Position on the straight S/W edges of the rectangle formed
      // by the corner points
      edge_S[i] =
        corner_SW[i] + (corner_SE[i] - corner_SW[i]) * 0.5 * (s[0] + 1.0);
      edge_N[i] =
        corner_NW[i] + (corner_NE[i] - corner_NW[i]) * 0.5 * (s[0] + 1.0);

      // Position inside rectangle
      f_rect[i] = edge_S[i] + (edge_N[i] - edge_S[i]) * 0.5 * (s[1] + 1.0);

      // Get difference between curved edge and point in rectangle
      diff_N[i] = bound_N[i] - f_rect[i];
      diff_S[i] = bound_S[i] - f_rect[i];
      diff_E[i] = bound_E[i] - f_rect[i];
      diff_W[i] = bound_W[i] - f_rect[i];

      // Map it...
      r[i] = f_rect[i] + diff_S[i] * (1.0 - 0.5 * (s[1] + 1.0)) +
             diff_N[i] * 0.5 * (s[1] + 1.0) +
             diff_W[i] * (1.0 - 0.5 * (s[0] + 1.0)) +
             diff_E[i] * 0.5 * (s[0] + 1.0);
    }
  }

  //=================================================================
  /// Get global position r(S) at discrete time level t.
  /// t=0: Present time; t>0: previous timestep.
  //=================================================================
  void QMacroElement<2>::macro_map(const unsigned& t,
                                   const Vector<double>& S,
                                   Vector<double>& r)
  {
    using namespace QuadTreeNames;

    Vector<double> bound_N(2);
    Vector<double> bound_S(2);
    Vector<double> bound_W(2);
    Vector<double> bound_E(2);

    Vector<double> diff_N(2);
    Vector<double> diff_S(2);
    Vector<double> diff_W(2);
    Vector<double> diff_E(2);

    Vector<double> f_rect(2);

    Vector<double> corner_SE(2);
    Vector<double> corner_SW(2);
    Vector<double> corner_NE(2);
    Vector<double> corner_NW(2);

    Vector<double> edge_N(2);
    Vector<double> edge_S(2);

    Vector<double> zeta(1);

    // Determine Vectors to corners
    zeta[0] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, corner_SE);
    zeta[0] = -1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, corner_SW);
    zeta[0] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, corner_NE);
    zeta[0] = -1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, corner_NW);


    // Get the position on the N/S/W/E edges
    zeta[0] = S[0];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, bound_N);
    zeta[0] = S[0];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, bound_S);
    zeta[0] = S[1];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::W, zeta, bound_W);
    zeta[0] = S[1];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::E, zeta, bound_E);


    for (int i = 0; i < 2; i++)
    {
      // Position on the straight S/W edges of the rectangle formed
      // by the corner points
      edge_S[i] =
        corner_SW[i] + (corner_SE[i] - corner_SW[i]) * 0.5 * (S[0] + 1.0);
      edge_N[i] =
        corner_NW[i] + (corner_NE[i] - corner_NW[i]) * 0.5 * (S[0] + 1.0);

      // Position inside rectangle
      f_rect[i] = edge_S[i] + (edge_N[i] - edge_S[i]) * 0.5 * (S[1] + 1.0);

      // Get difference between curved edge and point in rectangle
      diff_N[i] = bound_N[i] - f_rect[i];
      diff_S[i] = bound_S[i] - f_rect[i];
      diff_E[i] = bound_E[i] - f_rect[i];
      diff_W[i] = bound_W[i] - f_rect[i];

      // Map it...
      r[i] = f_rect[i] + diff_S[i] * (1.0 - 0.5 * (S[1] + 1.0)) +
             diff_N[i] * 0.5 * (S[1] + 1.0) +
             diff_W[i] * (1.0 - 0.5 * (S[0] + 1.0)) +
             diff_E[i] * 0.5 * (S[0] + 1.0);
    }
  }


  //=================================================================
  /// Output all macro element boundaries as tecplot zones
  //=================================================================
  void QMacroElement<2>::output_macro_element_boundaries(std::ostream& outfile,
                                                         const unsigned& nplot)
  {
    using namespace QuadTreeNames;

    Vector<double> s(1);
    Vector<double> f(2);
    // Dump at present time
    unsigned t = 0;
    for (unsigned idirect = N; idirect <= W; idirect++)
    {
      outfile << "ZONE I=" << nplot << std::endl;
      for (unsigned j = 0; j < nplot; j++)
      {
        s[0] = -1.0 + 2.0 * double(j) / double(nplot - 1);
        Domain_pt->macro_element_boundary(
          t, Macro_element_number, idirect, s, f);
        outfile << f[0] << " " << f[1] << std::endl;
      }
    }
  }


  //=============================================================================
  /// Assembles the jacobian of the mapping from the macro coordinates to
  /// the global coordinates
  //=============================================================================
  void QMacroElement<2>::assemble_macro_to_eulerian_jacobian(
    const unsigned& t, const Vector<double>& S, DenseMatrix<double>& jacobian)
  {
    using namespace QuadTreeNames;

    Vector<double> bound_N(2);
    Vector<double> bound_S(2);
    Vector<double> bound_W(2);
    Vector<double> bound_E(2);

    Vector<double> dbound_N(2);
    Vector<double> dbound_S(2);
    Vector<double> dbound_E(2);
    Vector<double> dbound_W(2);

    Vector<double> corner_SE(2);
    Vector<double> corner_SW(2);
    Vector<double> corner_NE(2);
    Vector<double> corner_NW(2);

    Vector<double> zeta(1);


    // Determine Vectors to corners
    zeta[0] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, corner_SE);
    zeta[0] = -1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, corner_SW);
    zeta[0] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, corner_NE);
    zeta[0] = -1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, corner_NW);


    // Get the position and first derivativeson the N/S/W/E edges
    zeta[0] = S[0];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, bound_N);
    Domain_pt->dmacro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, dbound_N);
    zeta[0] = S[0];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, bound_S);
    Domain_pt->dmacro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, dbound_S);
    zeta[0] = S[1];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::W, zeta, bound_W);
    Domain_pt->dmacro_element_boundary(
      t, Macro_element_number, QuadTreeNames::W, zeta, dbound_W);
    zeta[0] = S[1];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::E, zeta, bound_E);
    Domain_pt->dmacro_element_boundary(
      t, Macro_element_number, QuadTreeNames::E, zeta, dbound_E);


    // dr0/dm0
    jacobian(0, 0) =
      0.25 * (corner_SW[0] - corner_SE[0] + corner_NW[0] - corner_NE[0] -
              corner_NE[0] * S[1] + corner_NW[0] * S[1] + corner_SE[0] * S[1] -
              corner_SW[0] * S[1]) +
      0.5 * (dbound_S[0] + dbound_N[0] - bound_W[0] + bound_E[0] -
             dbound_S[0] * S[1] + dbound_N[0] * S[1]);
    // dr1/dm0
    jacobian(0, 1) =
      0.25 * (corner_SW[1] - corner_SE[1] + corner_NW[1] - corner_NE[1] -
              corner_NE[1] * S[1] + corner_NW[1] * S[1] + corner_SE[1] * S[1] -
              corner_SW[1] * S[1]) +
      0.5 * (dbound_S[1] + dbound_N[1] - bound_W[1] + bound_E[1] -
             dbound_S[1] * S[1] + dbound_N[1] * S[1]);
    // dr0/dm1
    jacobian(1, 0) =
      0.25 * (corner_SW[0] + corner_SE[0] - corner_NW[0] - corner_NE[0] +
              corner_SE[0] * S[0] - corner_SW[0] * S[0] - corner_NE[0] * S[0] +
              corner_NW[0] * S[0]) +
      0.5 * (-bound_S[0] + bound_N[0] + dbound_W[0] + dbound_E[0] -
             dbound_W[0] * S[0] + dbound_E[0] * S[0]);
    // dr1/dm1
    jacobian(1, 1) =
      0.25 * (corner_SW[1] + corner_SE[1] - corner_NW[1] - corner_NE[1] +
              corner_SE[1] * S[0] - corner_SW[1] * S[0] - corner_NE[1] * S[0] +
              corner_NW[1] * S[0]) +
      0.5 * (-bound_S[1] + bound_N[1] + dbound_W[1] + dbound_E[1] -
             dbound_W[1] * S[0] + dbound_E[1] * S[0]);
  }


  //=============================================================================
  /// Assembles the second derivative jacobian of the mapping from the
  /// macro coordinates to global coordinates x
  //=============================================================================
  void QMacroElement<2>::assemble_macro_to_eulerian_jacobian2(
    const unsigned& t, const Vector<double>& S, DenseMatrix<double>& jacobian2)
  {
    using namespace QuadTreeNames;

    Vector<double> bound_N(2);
    Vector<double> bound_S(2);
    Vector<double> bound_W(2);
    Vector<double> bound_E(2);

    Vector<double> dbound_N(2);
    Vector<double> dbound_S(2);
    Vector<double> dbound_E(2);
    Vector<double> dbound_W(2);

    Vector<double> d2bound_N(2);
    Vector<double> d2bound_S(2);
    Vector<double> d2bound_E(2);
    Vector<double> d2bound_W(2);

    Vector<double> corner_SE(2);
    Vector<double> corner_SW(2);
    Vector<double> corner_NE(2);
    Vector<double> corner_NW(2);

    Vector<double> zeta(1);


    // Determine Vectors to corners
    zeta[0] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, corner_SE);
    zeta[0] = -1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, corner_SW);
    zeta[0] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, corner_NE);
    zeta[0] = -1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, corner_NW);


    // Get the position and first derivativeson the N/S/W/E edges
    zeta[0] = S[0];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, bound_N);
    Domain_pt->dmacro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, dbound_N);
    Domain_pt->d2macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::N, zeta, d2bound_N);
    zeta[0] = S[0];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, bound_S);
    Domain_pt->dmacro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, dbound_S);
    Domain_pt->d2macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::S, zeta, d2bound_S);
    zeta[0] = S[1];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::W, zeta, bound_W);
    Domain_pt->dmacro_element_boundary(
      t, Macro_element_number, QuadTreeNames::W, zeta, dbound_W);
    Domain_pt->d2macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::W, zeta, d2bound_W);
    zeta[0] = S[1];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::E, zeta, bound_E);
    Domain_pt->dmacro_element_boundary(
      t, Macro_element_number, QuadTreeNames::E, zeta, dbound_E);
    Domain_pt->d2macro_element_boundary(
      t, Macro_element_number, QuadTreeNames::E, zeta, d2bound_E);


    // d2x0/dm0^2
    jacobian2(0, 0) = 0.5 * (d2bound_S[0] + d2bound_N[0] - d2bound_S[0] * S[1] +
                             d2bound_N[0] * S[1]);
    // d2x0/dm1^2
    jacobian2(1, 0) = 0.5 * (d2bound_W[0] + d2bound_E[0] - d2bound_W[0] * S[0] +
                             d2bound_E[0] * S[0]);
    // d2x0/dm0dm1
    jacobian2(2, 0) =
      0.25 * (-corner_NE[0] + corner_NW[0] + corner_SE[0] - corner_SW[0]) +
      0.5 * (-dbound_W[0] + dbound_E[0] - dbound_S[0] + dbound_N[0]);
    // d2x1/dm0^2
    jacobian2(0, 1) = 0.5 * (d2bound_S[1] + d2bound_N[1] - d2bound_S[1] * S[1] +
                             d2bound_N[1] * S[1]);
    // d2x1/dm1^2
    jacobian2(1, 1) = 0.5 * (d2bound_W[1] + d2bound_E[1] - d2bound_W[1] * S[0] +
                             d2bound_E[1] * S[0]);
    // d2x1/dm0dm1
    jacobian2(2, 1) =
      0.25 * (-corner_NE[1] + corner_NW[1] + corner_SE[1] - corner_SW[1]) +
      0.5 * (-dbound_W[1] + dbound_E[1] - dbound_S[1] + dbound_N[1]);
  }


  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //=================================================================
  /// Get global position r(S) at discrete time level t.
  /// t=0: Present time; t>0: previous timestep.
  //=================================================================
  void QMacroElement<3>::macro_map(const unsigned& t,
                                   const Vector<double>& S,
                                   Vector<double>& r)
  {
    // get the eight corners
    Vector<double> corner_LDB(3);
    Vector<double> corner_RDB(3);
    Vector<double> corner_LUB(3);
    Vector<double> corner_RUB(3);
    Vector<double> corner_LDF(3);
    Vector<double> corner_RDF(3);
    Vector<double> corner_LUF(3);
    Vector<double> corner_RUF(3);

    Vector<double> zeta(2);
    zeta[0] = -1.0;
    zeta[1] = -1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::B, zeta, corner_LDB);
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, corner_LUB);
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::F, zeta, corner_LDF);
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::R, zeta, corner_RDB);
    zeta[0] = 1.0;
    zeta[1] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::B, zeta, corner_RUB);
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, corner_RDF);
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::L, zeta, corner_LUF);
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::R, zeta, corner_RUF);


    // get the position of the 4 corners of the center slice
    Vector<double> corner_LD(3);
    Vector<double> corner_RD(3);
    Vector<double> corner_LU(3);
    Vector<double> corner_RU(3);

    zeta[0] = -1.0;
    zeta[1] = S[2];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, corner_LD);
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, corner_LU);
    zeta[0] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, corner_RD);
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, corner_RU);

    // get position on the B,F faces;
    Vector<double> face_B(3);
    Vector<double> face_F(3);

    zeta[0] = S[0];
    zeta[1] = S[1];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::B, zeta, face_B);
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::F, zeta, face_F);


    // get position on the edges of the middle slice
    Vector<double> edge_mid_L(3);
    Vector<double> edge_mid_R(3);
    Vector<double> edge_mid_D(3);
    Vector<double> edge_mid_U(3);
    zeta[0] = S[0];
    zeta[1] = S[2];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, edge_mid_U);

    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, edge_mid_D);
    zeta[0] = S[1];
    zeta[1] = S[2];
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::L, zeta, edge_mid_L);

    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::R, zeta, edge_mid_R);

    // get position on the edges of the back slice
    Vector<double> edge_back_L(3);
    Vector<double> edge_back_R(3);
    Vector<double> edge_back_D(3);
    Vector<double> edge_back_U(3);
    zeta[0] = S[0];
    zeta[1] = -1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, edge_back_U);

    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, edge_back_D);
    zeta[0] = S[1];
    zeta[1] = -1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::L, zeta, edge_back_L);

    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::R, zeta, edge_back_R);

    // get position on the edges of the front slice
    Vector<double> edge_front_L(3);
    Vector<double> edge_front_R(3);
    Vector<double> edge_front_D(3);
    Vector<double> edge_front_U(3);
    zeta[0] = S[0];
    zeta[1] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, edge_front_U);

    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, edge_front_D);
    zeta[0] = S[1];
    zeta[1] = 1.0;
    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::L, zeta, edge_front_L);

    Domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::R, zeta, edge_front_R);


    for (unsigned i = 0; i < 3; i++)
    {
      // Position on the middle slice
      // ============================
      double slice_mid;

      // the points on the up and down edges of the middle "rectangular slice"
      double ptUp, ptDo;
      ptUp = corner_LU[i] + (corner_RU[i] - corner_LU[i]) * 0.5 * (S[0] + 1.0);
      ptDo = corner_LD[i] + (corner_RD[i] - corner_LD[i]) * 0.5 * (S[0] + 1.0);
      // position in the rectangular middle slice
      slice_mid = ptDo + 0.5 * (1.0 + S[1]) * (ptUp - ptDo);

      // get the differences to the edges of the middle slice
      double diff_L, diff_R, diff_D, diff_U;
      diff_L = edge_mid_L[i] - slice_mid;
      diff_R = edge_mid_R[i] - slice_mid;
      diff_D = edge_mid_D[i] - slice_mid;
      diff_U = edge_mid_U[i] - slice_mid;

      // Map it to get the position in the middle slice
      slice_mid +=
        diff_L * (1.0 - 0.5 * (S[0] + 1.0)) + diff_R * 0.5 * (S[0] + 1.0) +
        diff_D * (1.0 - 0.5 * (S[1] + 1.0)) + diff_U * 0.5 * (S[1] + 1.0);


      // Position on the back slice
      //===========================
      double slice_back;

      // the points on the up and down edges of the back "rectangular slice"
      ptUp =
        corner_LUB[i] + (corner_RUB[i] - corner_LUB[i]) * 0.5 * (S[0] + 1.0);
      ptDo =
        corner_LDB[i] + (corner_RDB[i] - corner_LDB[i]) * 0.5 * (S[0] + 1.0);
      // position in the rectangular back slice
      slice_back = ptDo + 0.5 * (1.0 + S[1]) * (ptUp - ptDo);

      // get the differences to the edges of the middle slice
      diff_L = edge_back_L[i] - slice_back;
      diff_R = edge_back_R[i] - slice_back;
      diff_D = edge_back_D[i] - slice_back;
      diff_U = edge_back_U[i] - slice_back;

      // Map it to get the position in the back slice
      slice_back +=
        diff_L * (1.0 - 0.5 * (S[0] + 1.0)) + diff_R * 0.5 * (S[0] + 1.0) +
        diff_D * (1.0 - 0.5 * (S[1] + 1.0)) + diff_U * 0.5 * (S[1] + 1.0);

      // Position on the front slice
      //============================
      double slice_front;

      // the points on the up and down edges of the back "rectangular slice"
      ptUp =
        corner_LUF[i] + (corner_RUF[i] - corner_LUF[i]) * 0.5 * (S[0] + 1.0);
      ptDo =
        corner_LDF[i] + (corner_RDF[i] - corner_LDF[i]) * 0.5 * (S[0] + 1.0);
      // position in the rectangular back slice
      slice_front = ptDo + 0.5 * (1.0 + S[1]) * (ptUp - ptDo);

      // get the differences to the edges of the middle slice
      diff_L = edge_front_L[i] - slice_front;
      diff_R = edge_front_R[i] - slice_front;
      diff_D = edge_front_D[i] - slice_front;
      diff_U = edge_front_U[i] - slice_front;

      // Map it to get the position in the back slice
      slice_front +=
        diff_L * (1.0 - 0.5 * (S[0] + 1.0)) + diff_R * 0.5 * (S[0] + 1.0) +
        diff_D * (1.0 - 0.5 * (S[1] + 1.0)) + diff_U * 0.5 * (S[1] + 1.0);

      // Get difference between the back and front slices and the actual
      // boundary
      // ========================================================================

      double diff_back = face_B[i] - slice_back;
      double diff_front = face_F[i] - slice_front;

      // final map
      //==========

      r[i] = slice_mid + 0.5 * (1 + S[2]) * diff_front +
             0.5 * (1 - S[2]) * diff_back;
    }
  }


  //=================================================================
  /// Output all macro element boundaries as tecplot zones
  //=================================================================
  void QMacroElement<3>::output_macro_element_boundaries(std::ostream& outfile,
                                                         const unsigned& nplot)
  {
    using namespace OcTreeNames;

    Vector<double> s(2);
    Vector<double> f(3);
    // Dump at present time
    unsigned t = 0;
    for (unsigned idirect = L; idirect <= F; idirect++)
    {
      outfile << "ZONE I=" << nplot << ", J=" << nplot << std::endl;
      for (unsigned i = 0; i < nplot; i++)
      {
        s[1] = -1.0 + 2.0 * double(i) / double(nplot - 1);
        for (unsigned j = 0; j < nplot; j++)
        {
          s[0] = -1.0 + 2.0 * double(j) / double(nplot - 1);
          Domain_pt->macro_element_boundary(
            t, Macro_element_number, idirect, s, f);
          outfile << f[0] << " " << f[1] << " " << f[2] << std::endl;
        }
      }
    }
  }

} // namespace oomph
