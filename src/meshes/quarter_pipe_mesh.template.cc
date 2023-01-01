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
#ifndef OOMPH_QUARTER_PIPE_MESH_TEMPLATE_CC
#define OOMPH_QUARTER_PIPE_MESH_TEMPLATE_CC

#include "quarter_pipe_mesh.template.h"


namespace oomph
{
  //====================================================================
  /// Constructor: Pass number of elements in various directions,
  /// the inner and outer radius and the length of the tube
  //====================================================================
  template<class ELEMENT>
  QuarterPipeMesh<ELEMENT>::QuarterPipeMesh(const unsigned& ntheta,
                                            const unsigned& nr,
                                            const unsigned& nz,
                                            const double& rmin,
                                            const double& rmax,
                                            const double& length,
                                            TimeStepper* time_stepper_pt)
    : SimpleCubicMesh<ELEMENT>(
        ntheta, nr, nz, 1.0, 1.0, length, time_stepper_pt)
  {
    // Mesh can only be built with 3D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

    // Variables declaration
    Ntheta = ntheta;
    Nr = nr;
    Nz = nz;
    Rmin = rmin;
    Rmax = rmax;
    Length = length;

    // Build macro element-based domain
    Domain_pt = new QuarterPipeDomain(ntheta, nr, nz, rmin, rmax, length);

    // Loop over all elements
    unsigned nel = this->nelement();
    for (unsigned e = 0; e < nel; e++)
    {
      // Try to cast to FiniteElement
      FiniteElement* el_pt = dynamic_cast<FiniteElement*>(this->element_pt(e));

      // Set macro element pointer
      el_pt->set_macro_elem_pt(Domain_pt->macro_element_pt(e));
    }

    // Update node coordinates with macroelement coordinates,
    // updating solid coordinates too.
    this->node_update(true);

    // Setup boundary coordinates on inner boundary (boundary 1)
    unsigned b = 1;
    unsigned nnod = this->nboundary_node(b);
    for (unsigned j = 0; j < nnod; j++)
    {
      // Pointer to node
      Node* nod_pt = this->boundary_node_pt(b, j);

      // Get the Eulerian coordinates
      double x = nod_pt->x(0);
      double y = nod_pt->x(1);
      double z = nod_pt->x(2);

      // Polar angle
      double phi = atan2(y, x);

      // Set boundary coordinates
      Vector<double> zeta(2);
      zeta[0] = z;
      zeta[1] = phi;
      nod_pt->set_coordinates_on_boundary(b, zeta);
    }
    this->Boundary_coordinate_exists[b] = true;
  }

} // namespace oomph


#endif
