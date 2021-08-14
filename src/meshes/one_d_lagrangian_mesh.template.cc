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
#ifndef OOMPH_ONE_D_LAGRANGIAN_MESH_TEMPLATE_CC
#define OOMPH_ONE_D_LAGRANGIAN_MESH_TEMPLATE_CC

// The templated member functions of OneDLagrangianMesh
#include "one_d_lagrangian_mesh.template.h"

// Include the templated member functions of the OneDMesh
#include "one_d_mesh.template.cc"


namespace oomph
{
  //=======================================================================
  /// Constructor for 1D mesh:
  /// n_element   : number of elements
  /// length  : length of domain
  /// undef_eulerian_posn_pt: pointer to geom object that describes
  ///                         the initial Eulerian position of the mesh.
  /// time_stepper_pt  : timestepper
  //=======================================================================
  template<class ELEMENT>
  OneDLagrangianMesh<ELEMENT>::OneDLagrangianMesh(
    const unsigned& n_element,
    const double& length,
    GeomObject* undef_eulerian_posn_pt,
    TimeStepper* time_stepper_pt)
    : OneDMesh<ELEMENT>(n_element, length, time_stepper_pt),
      Undef_eulerian_posn_pt(undef_eulerian_posn_pt)
  {
    // Mesh can only be built with 1D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(1);

    // Now set the lagrangian coordinates of the nodes
    set_lagrangian_nodal_coordinates();

    // Set the default slopes
    assign_default_element_gradients();

    // Now set up the Eulerian position of the nodal points
    assign_undeformed_positions();
  }

  //==================================================
  /// Constructor for 1D mesh:
  /// n_element   : number of elements
  /// xmin        : minimum coordinate value (LH end)
  /// xmax        : maximum coordinate value (RH end)
  /// undef_eulerian_posn_pt: pointer to geom object that describes
  ///                         the initial Eulerian position of the mesh.
  /// time_stepper_pt  : timestepper
  //==================================================
  template<class ELEMENT>
  OneDLagrangianMesh<ELEMENT>::OneDLagrangianMesh(
    const unsigned& n_element,
    const double& xmin,
    const double& xmax,
    GeomObject* undef_eulerian_posn_pt,
    TimeStepper* time_stepper_pt)
    : OneDMesh<ELEMENT>(n_element, xmin, xmax, time_stepper_pt),
      Undef_eulerian_posn_pt(undef_eulerian_posn_pt)
  {
    // Mesh can only be built with 1D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(1);

    // Now set the lagrangian coordinates of the nodes
    set_lagrangian_nodal_coordinates();

    // Set the default slopes
    assign_default_element_gradients();

    // Now set up the Eulerian position of the nodal points
    assign_undeformed_positions();
  }


  //=====================================================================
  /// Set the default (initial) gradients within each element, which
  /// are merely the distances between the nodes, scaled by 0.5 because
  /// the elements have length 2 in local coordinates.
  /// N.B. This only works for QHermiteElements at the moment
  //=====================================================================
  template<class ELEMENT>
  void OneDLagrangianMesh<ELEMENT>::assign_default_element_gradients()
  {
    // Find cast pointer to first element
    ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(finite_element_pt(0));
    // Do we need to worry about the slopes
    // Read out number of position dofs
    unsigned n_lagrangian_type = cast_element_pt->nnodal_lagrangian_type();

    // If this is greater than 1 set the slopes, which are the distances between
    // nodes
    if (n_lagrangian_type > 1)
    {
      // Read out the number of linear points in the element
      unsigned n_p = cast_element_pt->nnode_1d();
      // Set the values of the distance between each node
      double xstep =
        OneDMesh<ELEMENT>::Length / double((n_p - 1) * OneDMesh<ELEMENT>::N);
      // Loop over the nodes and set the slopes
      unsigned long n_node = nnode();
      for (unsigned long n = 0; n < n_node; n++)
      {
        node_pt(n)->xi_gen(1, 0) = 0.5 * xstep;
      }
    }
  }

  //======================================================================
  /// Set the initial (2D Eulerian!) positions of the nodes
  //======================================================================
  template<class ELEMENT>
  void OneDLagrangianMesh<ELEMENT>::assign_undeformed_positions()
  {
    // Lagrangian coordinate
    Vector<double> xi(1);

    // Find the number Eulerian coordinates, assume same for all nodes
    unsigned n_dim = node_pt(0)->ndim();

    // Find cast pointer to first element
    ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(finite_element_pt(0));

    // Do we need to worry about the slopes
    // Read out number of position dofs
    unsigned n_lagrangian_type = cast_element_pt->nnodal_lagrangian_type();

    // Setup position Vector and derivatives (they *should* dimension themselves
    // in call to position() etc)
    Vector<double> R(n_dim);
    // Derivative wrt Lagrangian coordinate
    DenseMatrix<double> a(1, n_dim);
    // Second derivative wrt Lagrangian coordinate
    RankThreeTensor<double> dadxi(1, 1, n_dim);

    // Find out how many nodes there are
    unsigned long n_node = nnode();

    // Loop over all the nodes
    for (unsigned long n = 0; n < n_node; n++)
    {
      // Lagrangian coordinate of node
      xi[0] = node_pt(n)->xi(0);

      // Get the undeformed midplane
      Undef_eulerian_posn_pt->d2position(xi, R, a, dadxi);

      // Loop over coordinate directions
      for (unsigned i = 0; i < n_dim; i++)
      {
        // Set the position
        node_pt(n)->x_gen(0, i) = R[i];

        if (n_lagrangian_type > 1)
        {
          // Set the derivative wrt Lagrangian coordinates
          // Note that we need to scale by the length of each element here!!
          // and the 0.5 comes from the fact that our reference element has
          // length 2.0
          node_pt(n)->x_gen(1, i) =
            0.5 * a(0, i) *
            (OneDMesh<ELEMENT>::Length / double(OneDMesh<ELEMENT>::N));
        }
      }
    }
  }


} // namespace oomph
#endif
