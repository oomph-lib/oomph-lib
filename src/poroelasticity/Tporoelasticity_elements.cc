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
#include "Tporoelasticity_elements.h"

namespace oomph
{
  /// Lowest-order Raviart-Thomas based Darcy equation element

  /// Constructor. Order 0 elements have 1 pressure dof and no internal
  /// velocity dofs
  template<>
  TPoroelasticityElement<0>::TPoroelasticityElement()
    : TElement<2, 3>(), PoroelasticityEquations<2>(), Sign_edge(3, 1)
  {
    P_internal_data_index = this->add_internal_data(new Data(1));
  }

  /// Destructor
  template<>
  TPoroelasticityElement<0>::~TPoroelasticityElement()
  {
  }

  /// Return the total number of computational basis functions for u
  template<>
  unsigned TPoroelasticityElement<0>::nq_basis() const
  {
    return 3;
  }

  /// Return the number of edge basis functions for u
  template<>
  unsigned TPoroelasticityElement<0>::nq_basis_edge() const
  {
    return 3;
  }

  /// Returns the local form of the q basis at local coordinate s
  template<>
  void TPoroelasticityElement<0>::get_q_basis_local(const Vector<double>& s,
                                                    Shape& q_basis) const
  {
    // RT_0 basis functions
    q_basis(0, 0) = Sign_edge[0] * std::sqrt(2) * s[0];
    q_basis(0, 1) = Sign_edge[0] * std::sqrt(2) * s[1];

    q_basis(1, 0) = Sign_edge[1] * (s[0] - 1);
    q_basis(1, 1) = Sign_edge[1] * s[1];

    q_basis(2, 0) = Sign_edge[2] * s[0];
    q_basis(2, 1) = Sign_edge[2] * (s[1] - 1);
  }

  /// Returns the local form of the q basis and dbasis/ds at local coordinate s
  template<>
  void TPoroelasticityElement<0>::get_div_q_basis_local(
    const Vector<double>& s, Shape& div_q_basis_ds) const
  {
    div_q_basis_ds(0) = Sign_edge[0] * 2 * std::sqrt(2);
    div_q_basis_ds(1) = Sign_edge[1] * 2;
    div_q_basis_ds(2) = Sign_edge[2] * 2;

    // Scale the basis by the ratio of the length of the edge to the length of
    // the corresponding edge on the reference element
    scale_basis(div_q_basis_ds);
  }

  /// Returns the number of gauss points along each edge of the element
  template<>
  unsigned TPoroelasticityElement<0>::nedge_gauss_point() const
  {
    return 1;
  }

  /// Return the total number of pressure basis functions
  template<>
  unsigned TPoroelasticityElement<0>::np_basis() const
  {
    return 1;
  }

  /// Return the pressure basis
  template<>
  void TPoroelasticityElement<0>::get_p_basis(const Vector<double>& s,
                                              Shape& p_basis) const
  {
    p_basis(0) = 1.0;
  }

  /// The number of values stored at each node
  template<>
  const unsigned TPoroelasticityElement<0>::Initial_Nvalue[6] =
    //{0,0,0,1,1,1};
    {2, 2, 2, 3, 3, 3};

  /// Conversion scheme from an edge degree of freedom to the node it's stored
  /// at
  template<>
  const unsigned TPoroelasticityElement<0>::Q_edge_conv[3] = {3, 4, 5};

  /// The points along each edge where the fluxes are taken to be
  template<>
  const double TPoroelasticityElement<0>::Gauss_point[1] = {0.5};

  /// Second-order Raviart-Thomas based Darcy equation element

  /// Constructor. Order 1 elements have 3 pressure dofs and 2 internal
  /// velocity dofs
  template<>
  TPoroelasticityElement<1>::TPoroelasticityElement()
    : TElement<2, 3>(), PoroelasticityEquations<2>(), Sign_edge(3, 1)
  {
    // RT_1 elements have 2 internal degrees of freedom for u, and 3 for p
    Q_internal_data_index = this->add_internal_data(new Data(2));
    P_internal_data_index = this->add_internal_data(new Data(3));
  }

  /// Destructor
  template<>
  TPoroelasticityElement<1>::~TPoroelasticityElement()
  {
  }

  /// Return the total number of computational basis functions for u
  template<>
  unsigned TPoroelasticityElement<1>::nq_basis() const
  {
    return 8;
  }

  /// Return the number of edge basis functions for u
  template<>
  unsigned TPoroelasticityElement<1>::nq_basis_edge() const
  {
    return 6;
  }

  /// Returns the local form of the q basis at local coordinate s
  template<>
  void TPoroelasticityElement<1>::get_q_basis_local(const Vector<double>& s,
                                                    Shape& q_basis) const
  {
    // RT_1 basis functions

    double g1, g2;

    g1 = edge_gauss_point(0, 0);
    g2 = edge_gauss_point(0, 1);
    q_basis(0, 0) =
      Sign_edge[0] * std::sqrt(2) * s[0] * (s[1] - g2) / (g1 - g2);
    q_basis(0, 1) =
      Sign_edge[0] * std::sqrt(2) * s[1] * (s[1] - g2) / (g1 - g2);

    q_basis(1, 0) =
      Sign_edge[0] * std::sqrt(2) * s[0] * (s[1] - g1) / (g2 - g1);
    q_basis(1, 1) =
      Sign_edge[0] * std::sqrt(2) * s[1] * (s[1] - g1) / (g2 - g1);

    g1 = edge_gauss_point(1, 0);
    g2 = edge_gauss_point(1, 1);
    q_basis(2, 0) = Sign_edge[1] * (s[0] - 1) * (s[1] - g1) / (g2 - g1);
    q_basis(2, 1) = Sign_edge[1] * s[1] * (s[1] - g1) / (g2 - g1);

    q_basis(3, 0) = Sign_edge[1] * (s[0] - 1) * (s[1] - g2) / (g1 - g2);
    q_basis(3, 1) = Sign_edge[1] * s[1] * (s[1] - g2) / (g1 - g2);

    g1 = edge_gauss_point(2, 0);
    g2 = edge_gauss_point(2, 1);
    q_basis(4, 0) = Sign_edge[2] * s[0] * (s[0] - g2) / (g1 - g2);
    q_basis(4, 1) = Sign_edge[2] * (s[1] - 1) * (s[0] - g2) / (g1 - g2);

    q_basis(5, 0) = Sign_edge[2] * s[0] * (s[0] - g1) / (g2 - g1);
    q_basis(5, 1) = Sign_edge[2] * (s[1] - 1) * (s[0] - g1) / (g2 - g1);

    q_basis(6, 0) = s[1] * s[0];
    q_basis(6, 1) = s[1] * (s[1] - 1);

    q_basis(7, 0) = s[0] * (s[0] - 1);
    q_basis(7, 1) = s[0] * s[1];
  }

  /// Returns the local form of the q basis and dbasis/ds at local coordinate s
  template<>
  void TPoroelasticityElement<1>::get_div_q_basis_local(
    const Vector<double>& s, Shape& div_q_basis_ds) const
  {
    double g1, g2;

    g1 = edge_gauss_point(0, 0);
    g2 = edge_gauss_point(0, 1);
    div_q_basis_ds(0) =
      Sign_edge[0] * std::sqrt(2) * (3 * s[1] - 2 * g2) / (g1 - g2);
    div_q_basis_ds(1) =
      Sign_edge[0] * std::sqrt(2) * (2 * g1 - 3 * s[1]) / (g1 - g2);

    g1 = edge_gauss_point(1, 0);
    g2 = edge_gauss_point(1, 1);
    div_q_basis_ds(2) = Sign_edge[1] * (2 * g1 - 3 * s[1]) / (g1 - g2);
    div_q_basis_ds(3) = Sign_edge[1] * (3 * s[1] - 2 * g2) / (g1 - g2);

    g1 = edge_gauss_point(2, 0);
    g2 = edge_gauss_point(2, 1);
    div_q_basis_ds(4) = Sign_edge[2] * (3 * s[0] - 2 * g2) / (g1 - g2);
    div_q_basis_ds(5) = Sign_edge[2] * (2 * g1 - 3 * s[0]) / (g1 - g2);

    div_q_basis_ds(6) = 3 * s[1] - 1;
    div_q_basis_ds(7) = 3 * s[0] - 1;

    // Scale the basis by the ratio of the length of the edge to the length of
    // the corresponding edge on the reference element
    scale_basis(div_q_basis_ds);
  }

  /// Returns the number of gauss points along each edge of the element
  template<>
  unsigned TPoroelasticityElement<1>::nedge_gauss_point() const
  {
    return 2;
  }

  /// Return the total number of pressure basis functions
  template<>
  unsigned TPoroelasticityElement<1>::np_basis() const
  {
    return 3;
  }

  /// Return the pressure basis
  template<>
  void TPoroelasticityElement<1>::get_p_basis(const Vector<double>& s,
                                              Shape& p_basis) const
  {
    p_basis(0) = 1.0;
    p_basis(1) = s[0];
    p_basis(2) = s[1];
  }

  /// The number of values stored at each node
  template<>
  const unsigned TPoroelasticityElement<1>::Initial_Nvalue[6] = {
    2, 2, 2, 4, 4, 4};

  /// Conversion scheme from an edge degree of freedom to the node it's stored
  /// at
  template<>
  const unsigned TPoroelasticityElement<1>::Q_edge_conv[3] = {3, 4, 5};

  /// The points along each edge where the fluxes are taken to be
  template<>
  const double TPoroelasticityElement<1>::Gauss_point[2] = {
    0.5 - std::sqrt(3.0) / 6.0, 0.5 + std::sqrt(3.0) / 6.0};

  // Force building of templates
  template class TPoroelasticityElement<0>;
  template class TPoroelasticityElement<1>;

} // namespace oomph
