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
#ifndef OOMPH_ANNULAR_MESH_TEMPLATE_CC
#define OOMPH_ANNULAR_MESH_TEMPLATE_CC

#include "annular_mesh.template.h"

namespace oomph
{
  //===================================================================
  /// Helper function to Wrap mesh into annular shape
  //==================================================================
  template<class ELEMENT>
  void TwoDAnnularMesh<ELEMENT>::wrap_into_annular_shape(
    const double& a,
    const double& h,
    const double& azimuthal_fraction,
    const double& phi)
  {
    // Create the hole
    Ellipse ellipse(a, a);

    // Set all the initial positions of the nodes
    Vector<double> xi(1);
    Vector<double> base(2);
    Vector<double> N(2);
    const unsigned n_node = this->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      // Pointer to node
      Node* nod_pt = this->node_pt(n);

      // Get the angle of the node -- rotate such that jump in angle
      // appears at periodic boundary. Shrink domain slightly
      // to keep angle unique
      xi[0] = (1.0 - 1.0e-10) * (-azimuthal_fraction * nod_pt->x(0)) * 2.0 *
                MathematicalConstants::Pi +
              MathematicalConstants::Pi -
              (1.0 - azimuthal_fraction) * 2.0 * MathematicalConstants::Pi;

      // Rotate
      xi[0] += phi;

      // Get the node's fraction in the radial direction
      double w = nod_pt->x(1);

      // Get the position on the ellipse base
      ellipse.position(xi, base);

      // Get the unit normal, if it were a circle , by normalising the base
      double norm = sqrt(base[0] * base[0] + base[1] * base[1]);
      N[0] = base[0] / norm;
      N[1] = base[1] / norm;

      // Set circular film from the ellipse
      nod_pt->x(0) = base[0] + w * (h + a - norm) * N[0];
      nod_pt->x(1) = base[1] + w * (h + a - norm) * N[1];

      // Set boundary coordinates
      Vector<double> xi_bound(1);

      // Polar angle for boundary coordinate on boundary 0
      if (nod_pt->is_on_boundary(0))
      {
        xi_bound[0] = atan2(nod_pt->x(1), nod_pt->x(0));
        nod_pt->set_coordinates_on_boundary(0, xi_bound);
      }

      // Radius for boundary coordinate on boundary 1
      if (nod_pt->is_on_boundary(1))
      {
        xi_bound[0] = sqrt(pow(nod_pt->x(0), 2) + pow(nod_pt->x(1), 2));
        nod_pt->set_coordinates_on_boundary(1, xi_bound);
      }

      // Polar angle for boundary coordinate on boundary 2
      if (nod_pt->is_on_boundary(2))
      {
        xi_bound[0] = atan2(nod_pt->x(1), nod_pt->x(0));
        nod_pt->set_coordinates_on_boundary(2, xi_bound);
      }

      // Radius for boundary coordinate on boundary 3
      if (nod_pt->is_on_boundary(3))
      {
        xi_bound[0] = sqrt(pow(nod_pt->x(0), 2) + pow(nod_pt->x(1), 2));
        nod_pt->set_coordinates_on_boundary(3, xi_bound);
      }
    }

    this->Boundary_coordinate_exists[0] = true;
    this->Boundary_coordinate_exists[1] = true;
    this->Boundary_coordinate_exists[2] = true;
    this->Boundary_coordinate_exists[3] = true;
  }


} // namespace oomph


#endif
