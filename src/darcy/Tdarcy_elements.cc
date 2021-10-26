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
#include "Tdarcy_elements.h"

namespace oomph
{
  /// ///////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////
  /// Lowest-order Raviart-Thomas based Darcy equation element
  /// ///////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////


  //===========================================================================
  /// Constructor. Order 0 elements have 1 pressure dof and no internal
  /// velocity dofs
  //===========================================================================
  template<>
  TRaviartThomasDarcyElement<0>::TRaviartThomasDarcyElement()
    : TElement<2, 3>(), DarcyEquations<2>(), Sign_edge(3, 1)
  {
    P_internal_data_index = this->add_internal_data(new Data(1));
  }

  //===========================================================================
  /// Destructor
  //===========================================================================
  template<>
  TRaviartThomasDarcyElement<0>::~TRaviartThomasDarcyElement()
  {
  }


  //===========================================================================
  /// Return the number of edge basis functions for q
  //===========================================================================
  template<>
  unsigned TRaviartThomasDarcyElement<0>::nq_basis_edge() const
  {
    return 3;
  }


  //===========================================================================
  /// Return the number of internal basis functions for q
  //===========================================================================
  template<>
  unsigned TRaviartThomasDarcyElement<0>::nq_basis_internal() const
  {
    return 0;
  }

  //===========================================================================
  /// Compute the local form of the q basis at local coordinate s
  //===========================================================================
  template<>
  void TRaviartThomasDarcyElement<0>::get_q_basis_local(const Vector<double>& s,
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

  //===========================================================================
  /// Compute the local form of the q basis and dbasis/ds at local coordinate s
  //===========================================================================
  template<>
  void TRaviartThomasDarcyElement<0>::get_div_q_basis_local(
    const Vector<double>& s, Shape& div_q_basis_ds) const
  {
    div_q_basis_ds(0) = Sign_edge[0] * 2 * std::sqrt(2);
    div_q_basis_ds(1) = Sign_edge[1] * 2;
    div_q_basis_ds(2) = Sign_edge[2] * 2;

    // Scale the basis by the ratio of the length of the edge to the length of
    // the corresponding edge on the reference element
    // hierher explain please
    scale_basis(div_q_basis_ds);
  }

  //===========================================================================
  /// Return the number of flux interpolation points along each edge
  /// of the element
  //===========================================================================
  template<>
  unsigned TRaviartThomasDarcyElement<0>::nedge_flux_interpolation_point() const
  {
    return 1;
  }

  //===========================================================================
  /// Return the equation number of the n-th pressure degree of freedom
  //===========================================================================
  template<>
  int TRaviartThomasDarcyElement<0>::p_local_eqn(const unsigned& n) const
  {
    return this->internal_local_eqn(P_internal_data_index, n);
  }

  //===========================================================================
  /// Return the nth pressure value
  //===========================================================================
  template<>
  double TRaviartThomasDarcyElement<0>::p_value(const unsigned& n) const
  {
    return this->internal_data_pt(P_internal_data_index)->value(n);
  }

  //===========================================================================
  /// Return the total number of pressure basis functions
  //===========================================================================
  template<>
  unsigned TRaviartThomasDarcyElement<0>::np_basis() const
  {
    return 1;
  }

  //===========================================================================
  /// Compute the pressure basis
  //===========================================================================
  template<>
  void TRaviartThomasDarcyElement<0>::get_p_basis(const Vector<double>& s,
                                                  Shape& p_basis) const
  {
    p_basis(0) = 1.0;
  }


  //===========================================================================
  /// Recovery order for Z2 error estimator
  //===========================================================================
  template<>
  unsigned TRaviartThomasDarcyElement<0>::nrecovery_order()
  {
    return 2; // Need to experiment with this...
  }

  //===========================================================================
  /// The number of values stored at each node: A single flux dof is
  /// stored at each midside node
  //===========================================================================
  template<>
  const unsigned TRaviartThomasDarcyElement<0>::Initial_Nvalue[6] = {
    0, 0, 0, 1, 1, 1};


  //===========================================================================
  ///  Face index associated with edge flux degree of freedom
  //===========================================================================
  template<>
  const unsigned TRaviartThomasDarcyElement<0>::Face_index_of_edge_flux[3] = {
    2, 0, 1};


  //===========================================================================
  /// Conversion scheme from an edge degree of freedom to the node it's stored
  /// at
  /// Fluxes are stored at mid-side nodes
  //===========================================================================
  template<>
  const unsigned TRaviartThomasDarcyElement<0>::Q_edge_conv[3] = {3, 4, 5};


  //===========================================================================
  /// The points along each edge where the fluxes are taken to be
  // hierher explain please
  //===========================================================================
  template<>
  const double TRaviartThomasDarcyElement<0>::Flux_interpolation_point[1] = {
    0.5};


  /// ///////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////
  /// Second-order Raviart-Thomas based Darcy equation element
  /// ///////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////


  //===========================================================================
  /// Constructor. Order 1 elements have 3 pressure dofs and 2 internal
  /// velocity dofs
  //===========================================================================
  template<>
  TRaviartThomasDarcyElement<1>::TRaviartThomasDarcyElement()
    : TElement<2, 3>(), DarcyEquations<2>(), Sign_edge(3, 1)
  {
    // RT_1 elements have 2 internal degrees of freedom for q, and 3 for p
    Q_internal_data_index = this->add_internal_data(new Data(2));
    P_internal_data_index = this->add_internal_data(new Data(3));
  }

  //===========================================================================
  /// Destructor
  //===========================================================================
  template<>
  TRaviartThomasDarcyElement<1>::~TRaviartThomasDarcyElement()
  {
  }


  //===========================================================================
  /// Return the number of edge basis functions for q
  //===========================================================================
  template<>
  unsigned TRaviartThomasDarcyElement<1>::nq_basis_edge() const
  {
    return 6;
  }


  //===========================================================================
  /// Return the number of internal basis functions for q
  //===========================================================================
  template<>
  unsigned TRaviartThomasDarcyElement<1>::nq_basis_internal() const
  {
    return 2;
  }

  //===========================================================================
  /// Compute the local form of the q basis at local coordinate s
  //===========================================================================
  template<>
  void TRaviartThomasDarcyElement<1>::get_q_basis_local(const Vector<double>& s,
                                                        Shape& q_basis) const
  {
    // RT_1 basis functions

    Vector<double> g1_vect = edge_flux_interpolation_point(0, 0);
    Vector<double> g2_vect = edge_flux_interpolation_point(0, 1);
    double g1 = g1_vect[0];
    double g2 = g2_vect[0];
    q_basis(0, 0) =
      Sign_edge[0] * std::sqrt(2.0) * s[0] * (s[1] - g2) / (g1 - g2);
    q_basis(0, 1) =
      Sign_edge[0] * std::sqrt(2.0) * s[1] * (s[1] - g2) / (g1 - g2);

    q_basis(1, 0) =
      Sign_edge[0] * std::sqrt(2.0) * s[0] * (s[1] - g1) / (g2 - g1);
    q_basis(1, 1) =
      Sign_edge[0] * std::sqrt(2.0) * s[1] * (s[1] - g1) / (g2 - g1);


    g1_vect = edge_flux_interpolation_point(1, 0);
    g2_vect = edge_flux_interpolation_point(1, 1);
    g1 = g1_vect[0];
    g2 = g2_vect[0];
    q_basis(2, 0) = Sign_edge[1] * (s[0] - 1.0) * (s[1] - g1) / (g2 - g1);
    q_basis(2, 1) = Sign_edge[1] * s[1] * (s[1] - g1) / (g2 - g1);

    q_basis(3, 0) = Sign_edge[1] * (s[0] - 1.0) * (s[1] - g2) / (g1 - g2);
    q_basis(3, 1) = Sign_edge[1] * s[1] * (s[1] - g2) / (g1 - g2);


    g1_vect = edge_flux_interpolation_point(2, 0);
    g2_vect = edge_flux_interpolation_point(2, 1);
    g1 = g1_vect[0];
    g2 = g2_vect[0];
    q_basis(4, 0) = Sign_edge[2] * s[0] * (s[0] - g2) / (g1 - g2);
    q_basis(4, 1) = Sign_edge[2] * (s[1] - 1.0) * (s[0] - g2) / (g1 - g2);

    q_basis(5, 0) = Sign_edge[2] * s[0] * (s[0] - g1) / (g2 - g1);
    q_basis(5, 1) = Sign_edge[2] * (s[1] - 1.0) * (s[0] - g1) / (g2 - g1);

    q_basis(6, 0) = s[1] * s[0];
    q_basis(6, 1) = s[1] * (s[1] - 1.0);

    q_basis(7, 0) = s[0] * (s[0] - 1.0);
    q_basis(7, 1) = s[0] * s[1];
  }


  //===========================================================================
  /// Compute the local form of the q basis and dbasis/ds at local coordinate s
  //===========================================================================
  template<>
  void TRaviartThomasDarcyElement<1>::get_div_q_basis_local(
    const Vector<double>& s, Shape& div_q_basis_ds) const
  {
    Vector<double> g1_vect = edge_flux_interpolation_point(0, 0);
    Vector<double> g2_vect = edge_flux_interpolation_point(0, 1);
    double g1 = g1_vect[0];
    double g2 = g2_vect[0];
    div_q_basis_ds(0) =
      Sign_edge[0] * std::sqrt(2.0) * (3.0 * s[1] - 2.0 * g2) / (g1 - g2);
    div_q_basis_ds(1) =
      Sign_edge[0] * std::sqrt(2.0) * (2.0 * g1 - 3.0 * s[1]) / (g1 - g2);


    g1_vect = edge_flux_interpolation_point(1, 0);
    g2_vect = edge_flux_interpolation_point(1, 1);
    g1 = g1_vect[0];
    g2 = g2_vect[0];
    div_q_basis_ds(2) = Sign_edge[1] * (2.0 * g1 - 3.0 * s[1]) / (g1 - g2);
    div_q_basis_ds(3) = Sign_edge[1] * (3.0 * s[1] - 2.0 * g2) / (g1 - g2);


    g1_vect = edge_flux_interpolation_point(2, 0);
    g2_vect = edge_flux_interpolation_point(2, 1);
    g1 = g1_vect[0];
    g2 = g2_vect[0];
    div_q_basis_ds(4) = Sign_edge[2] * (3.0 * s[0] - 2.0 * g2) / (g1 - g2);
    div_q_basis_ds(5) = Sign_edge[2] * (2.0 * g1 - 3.0 * s[0]) / (g1 - g2);

    div_q_basis_ds(6) = 3.0 * s[1] - 1.0;
    div_q_basis_ds(7) = 3.0 * s[0] - 1.0;

    // Scale the basis by the ratio of the length of the edge to the length of
    // the corresponding edge on the reference element
    scale_basis(div_q_basis_ds);
  }

  //===========================================================================
  /// Return the number of flux_interpolation points along each edge of
  /// the element
  //===========================================================================
  template<>
  unsigned TRaviartThomasDarcyElement<1>::nedge_flux_interpolation_point() const
  {
    return 2;
  }

  //===========================================================================
  /// Return the equation number of the n-th pressure degree of freedom
  //===========================================================================
  template<>
  int TRaviartThomasDarcyElement<1>::p_local_eqn(const unsigned& n) const
  {
    return this->internal_local_eqn(P_internal_data_index, n);
  }

  //===========================================================================
  /// Return the nth pressure value
  //===========================================================================
  template<>
  double TRaviartThomasDarcyElement<1>::p_value(const unsigned& n) const
  {
    return this->internal_data_pt(P_internal_data_index)->value(n);
  }

  //===========================================================================
  /// Return the total number of pressure basis functions
  //===========================================================================
  template<>
  unsigned TRaviartThomasDarcyElement<1>::np_basis() const
  {
    return 3;
  }

  //===========================================================================
  /// Compute the pressure basis
  //===========================================================================
  template<>
  void TRaviartThomasDarcyElement<1>::get_p_basis(const Vector<double>& s,
                                                  Shape& p_basis) const
  {
    p_basis(0) = 1.0;
    p_basis(1) = s[0];
    p_basis(2) = s[1];
  }


  //===========================================================================
  /// Recovery order for Z2 error estimator
  //===========================================================================
  template<>
  unsigned TRaviartThomasDarcyElement<1>::nrecovery_order()
  {
    return 2; // Need to experiment with this...
  }

  //===========================================================================
  /// The number of values stored at each node. Two flux values are stored
  /// at each midside node.
  //===========================================================================
  template<>
  const unsigned TRaviartThomasDarcyElement<1>::Initial_Nvalue[6] = {
    0, 0, 0, 2, 2, 2};


  //===========================================================================
  ///  Face index associated with edge flux degree of freedom
  //===========================================================================
  template<>
  const unsigned TRaviartThomasDarcyElement<1>::Face_index_of_edge_flux[6] = {
    2, 2, 0, 0, 1, 1};

  //===========================================================================
  /// Conversion scheme from an edge degree of freedom to the node it's stored
  /// at
  /// Fluxes are stored at midside nodes
  //===========================================================================
  template<>
  const unsigned TRaviartThomasDarcyElement<1>::Q_edge_conv[3] = {3, 4, 5};


  //===========================================================================
  /// The points along each edge where the fluxes are taken to be
  //===========================================================================
  template<>
  const double TRaviartThomasDarcyElement<1>::Flux_interpolation_point[2] = {
    0.5 - std::sqrt(3.0) / 6.0, 0.5 + std::sqrt(3.0) / 6.0};


  //===========================================================================
  // Force building of templates
  //===========================================================================
  template class TRaviartThomasDarcyElement<0>;
  template class TRaviartThomasDarcyElement<1>;

} // namespace oomph
