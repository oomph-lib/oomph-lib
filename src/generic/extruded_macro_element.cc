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
// Necessary oomph-lib headers
#include "extruded_macro_element.h"
#include "extruded_domain.h"

namespace oomph
{
  //=================================================================
  /// Get global position r(s) at the continuous time, t
  //=================================================================
  void QExtrudedMacroElement<3>::macro_map(const unsigned& t,
                                           const Vector<double>& s,
                                           Vector<double>& r)
  {
    // Make sure that t=0 otherwise this doesn't make sense
    if (t != 0)
    {
      // Create an output stream
      std::ostringstream error_message_stream;

      // Create an error message
      error_message_stream << "This output function outputs a space-time\n"
                           << "representation of the solution. As such, it\n"
                           << "does not make sense to output the solution\n"
                           << "at a previous time level!" << std::endl;

      // Throw an error
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Storage for the eight corners
    Vector<double> corner_LDB(3, 0.0);
    Vector<double> corner_RDB(3, 0.0);
    Vector<double> corner_LUB(3, 0.0);
    Vector<double> corner_RUB(3, 0.0);
    Vector<double> corner_LDF(3, 0.0);
    Vector<double> corner_RDF(3, 0.0);
    Vector<double> corner_LUF(3, 0.0);
    Vector<double> corner_RUF(3, 0.0);

    // Lagrangian coordinates of a point on a 2D surface in 3D space
    Vector<double> zeta(2, 0.0);

    // First coordinate of the bottom-left coordinates of a surface
    zeta[0] = -1.0;

    // Second coordinate of the bottom-left coordinates of a surface
    zeta[1] = -1.0;

    // Calculate the space-time coordinates of the LDB corner
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::B, zeta, corner_LDB);

    // Calculate the space-time coordinates of the LUB corner
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, corner_LUB);

    // Calculate the space-time coordinates of the LDF corner
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::F, zeta, corner_LDF);

    // Calculate the space-time coordinates of the RDB corner
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::R, zeta, corner_RDB);

    // First coordinate of the the top-right coordinates of a surface
    zeta[0] = 1.0;

    // Second coordinate of the the top-right coordinates of a surface
    zeta[1] = 1.0;

    // Calculate the space-time coordinates of the RUB corner
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::B, zeta, corner_RUB);

    // Calculate the space-time coordinates of the RDF corner
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, corner_RDF);

    // Calculate the space-time coordinates of the LUF corner
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::L, zeta, corner_LUF);

    // Calculate the space-time coordinates of the RUF corner
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::R, zeta, corner_RUF);

    // Get the position of the 4 corners of the center slice
    Vector<double> corner_LD(3, 0.0);
    Vector<double> corner_RD(3, 0.0);
    Vector<double> corner_LU(3, 0.0);
    Vector<double> corner_RU(3, 0.0);

    // Set the coordinates of the point on a surface
    zeta[0] = -1.0;
    zeta[1] = s[2];

    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, corner_LD);
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, corner_LU);
    zeta[0] = 1.0;
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, corner_RD);
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, corner_RU);

    // Get position on the B,F faces;
    Vector<double> face_B(3);
    Vector<double> face_F(3);

    zeta[0] = s[0];
    zeta[1] = s[1];
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::B, zeta, face_B);
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::F, zeta, face_F);


    // Get position on the edges of the middle slice
    Vector<double> edge_mid_L(3);
    Vector<double> edge_mid_R(3);
    Vector<double> edge_mid_D(3);
    Vector<double> edge_mid_U(3);
    zeta[0] = s[0];
    zeta[1] = s[2];
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, edge_mid_U);

    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, edge_mid_D);
    zeta[0] = s[1];
    zeta[1] = s[2];
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::L, zeta, edge_mid_L);

    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::R, zeta, edge_mid_R);

    // Get position on the edges of the back slice
    Vector<double> edge_back_L(3);
    Vector<double> edge_back_R(3);
    Vector<double> edge_back_D(3);
    Vector<double> edge_back_U(3);
    zeta[0] = s[0];
    zeta[1] = -1.0;
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, edge_back_U);

    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, edge_back_D);
    zeta[0] = s[1];
    zeta[1] = -1.0;
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::L, zeta, edge_back_L);

    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::R, zeta, edge_back_R);

    // Get position on the edges of the front slice
    Vector<double> edge_front_L(3);
    Vector<double> edge_front_R(3);
    Vector<double> edge_front_D(3);
    Vector<double> edge_front_U(3);
    zeta[0] = s[0];
    zeta[1] = 1.0;
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::U, zeta, edge_front_U);

    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::D, zeta, edge_front_D);
    zeta[0] = s[1];
    zeta[1] = 1.0;
    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::L, zeta, edge_front_L);

    Extruded_domain_pt->macro_element_boundary(
      t, Macro_element_number, OcTreeNames::R, zeta, edge_front_R);

    // The number of dimensions (=space+time)
    unsigned n_dim = 3;

    // Loop over the coordinate directions
    for (unsigned i = 0; i < n_dim; i++)
    {
      //-----------------------------
      // Position on the middle slice
      //-----------------------------
      double slice_mid;

      // The points on the up and down edges of the middle "rectangular slice"
      double point_up = 0.0;
      double point_down = 0.0;

      // Calculate the point on the upper edge
      point_up =
        corner_LU[i] + (corner_RU[i] - corner_LU[i]) * 0.5 * (s[0] + 1.0);

      // Calculate the point on the lower edge
      point_down =
        corner_LD[i] + (corner_RD[i] - corner_LD[i]) * 0.5 * (s[0] + 1.0);

      // Position in the rectangular middle slice
      slice_mid = point_down + 0.5 * (1.0 + s[1]) * (point_up - point_down);

      // Get the differences to the edges of the middle slice
      double diff_L, diff_R, diff_D, diff_U;
      diff_L = edge_mid_L[i] - slice_mid;
      diff_R = edge_mid_R[i] - slice_mid;
      diff_D = edge_mid_D[i] - slice_mid;
      diff_U = edge_mid_U[i] - slice_mid;

      // Map it to get the position in the middle slice
      slice_mid +=
        (diff_L * (1.0 - 0.5 * (s[0] + 1.0)) + diff_R * 0.5 * (s[0] + 1.0) +
         diff_D * (1.0 - 0.5 * (s[1] + 1.0)) + diff_U * 0.5 * (s[1] + 1.0));

      //---------------------------
      // Position on the back slice
      //---------------------------
      double slice_back;

      // Calculate the point on the upper edge of the back "rectangular slice"
      point_up =
        corner_LUB[i] + (corner_RUB[i] - corner_LUB[i]) * 0.5 * (s[0] + 1.0);

      // Calculate the point on the lower edge of the back "rectangular slice"
      point_down =
        corner_LDB[i] + (corner_RDB[i] - corner_LDB[i]) * 0.5 * (s[0] + 1.0);

      // Position in the rectangular back slice
      slice_back = point_down + 0.5 * (1.0 + s[1]) * (point_up - point_down);

      // Get the differences to the edges of the middle slice
      diff_L = edge_back_L[i] - slice_back;
      diff_R = edge_back_R[i] - slice_back;
      diff_D = edge_back_D[i] - slice_back;
      diff_U = edge_back_U[i] - slice_back;

      // Map it to get the position in the back slice
      slice_back +=
        (diff_L * (1.0 - 0.5 * (s[0] + 1.0)) + diff_R * 0.5 * (s[0] + 1.0) +
         diff_D * (1.0 - 0.5 * (s[1] + 1.0)) + diff_U * 0.5 * (s[1] + 1.0));

      //----------------------------
      // Position on the front slice
      //----------------------------
      double slice_front;

      // Calculate the point on the upper edge of the back "rectangular slice"
      point_up =
        corner_LUF[i] + (corner_RUF[i] - corner_LUF[i]) * 0.5 * (s[0] + 1.0);

      // Calculate the point on the lower edge of the back "rectangular slice"
      point_down =
        corner_LDF[i] + (corner_RDF[i] - corner_LDF[i]) * 0.5 * (s[0] + 1.0);

      // Position in the rectangular back slice
      slice_front = point_down + 0.5 * (1.0 + s[1]) * (point_up - point_down);

      // Get the differences to the edges of the middle slice
      diff_L = edge_front_L[i] - slice_front;
      diff_R = edge_front_R[i] - slice_front;
      diff_D = edge_front_D[i] - slice_front;
      diff_U = edge_front_U[i] - slice_front;

      // Map it to get the position in the back slice
      slice_front +=
        (diff_L * (1.0 - 0.5 * (s[0] + 1.0)) + diff_R * 0.5 * (s[0] + 1.0) +
         diff_D * (1.0 - 0.5 * (s[1] + 1.0)) + diff_U * 0.5 * (s[1] + 1.0));

      //-------------------------------------------------------------------------
      // Get difference between the back and front slices and the actual
      // boundary
      //-------------------------------------------------------------------------
      double diff_back = face_B[i] - slice_back;
      double diff_front = face_F[i] - slice_front;

      //----------
      // Final map
      //----------
      r[i] = slice_mid + 0.5 * (1.0 + s[2]) * diff_front +
             0.5 * (1.0 - s[2]) * diff_back;
    }
  } // End of macro_map

  //=================================================================
  /// Output all macro element boundaries as tecplot zones
  //=================================================================
  void QExtrudedMacroElement<3>::output_macro_element_boundaries(
    std::ostream& outfile, const unsigned& nplot)
  {
    // Use the OcTree enumeration for corners, edges and faces
    using namespace OcTreeNames;

    // Storage for the local coordinates (of a point on a face)
    Vector<double> s(2);

    // Storage for the global coordinates
    Vector<double> x(3);

    // Loop over the faces
    for (unsigned idirect = L; idirect <= F; idirect++)
    {
      // Output the header
      outfile << "ZONE I=" << nplot << ", J=" << nplot << std::endl;

      // Loop over the plot points in the second surface direction
      for (unsigned i = 0; i < nplot; i++)
      {
        // Calculate the second local coordinate associated with this plot point
        s[1] = -1.0 + 2.0 * double(i) / double(nplot - 1);

        // Loop over the plot points in the first surface direction
        for (unsigned j = 0; j < nplot; j++)
        {
          // Calculate the first local coordinate associated with this plot
          // point
          s[0] = -1.0 + 2.0 * double(j) / double(nplot - 1);

          // To make the extrusion machinery consistent with the Domain
          // machinery a time level argument has to be provided. For our
          // purposes we set this to zero to ensure that the appropriate output
          // is provided.
          unsigned t = 0;

          // Get the global coordinates associated with these local coordinates
          Extruded_domain_pt->macro_element_boundary(
            t, Macro_element_number, idirect, s, x);

          // Output the global (space-time) coordinates
          outfile << x[0] << " " << x[1] << " " << x[2] << std::endl;
        }
      } // for (unsigned i=0;i<nplot;i++)
    } // for (unsigned idirect=L;idirect<=F;idirect++)
  } // End of output_macro_element_boundaries
} // End of namespace oomph
