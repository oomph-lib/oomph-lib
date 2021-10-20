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
// Header for an element that couples a linearised axisymmetric
// Navier-Stokes element to a non-linear axisymmetric Navier-Stokes
// element via a multi domain approach

// oomph-lib headers
#include "generic.h"
#include "navier_stokes.h"
#include "linearised_navier_stokes_elements.h"
#include "refineable_linearised_navier_stokes_elements.h"

// Use the oomph namespace
using namespace oomph;


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


//======================================================================
/// Build a LinearisedQTaylorHood element that inherits from
/// ElementWithExternalElement so that it can "communicate" with an
/// axisymmetric Navier-Stokes element that provides the base flow
/// velocities and their derivatives w.r.t. global coordinates (r and z)
//======================================================================
class LinearisedQTaylorHoodMultiDomainElement
  : public virtual LinearisedQTaylorHoodElement,
    public virtual ElementWithExternalElement
{
public:
  /// Constructor: call the underlying constructors
  LinearisedQTaylorHoodMultiDomainElement()
    : LinearisedQTaylorHoodElement(), ElementWithExternalElement()
  {
    // There are two interactions: the base flow velocities and their
    // derivatives w.r.t. global coordinates
    this->set_ninteraction(2);

    // Do not include any external interaction data when computing
    // the element's Jacobian
    // ElementWithExternalElement::ignore_external_interaction_data();

    /// Do not include any external geometric data when computing
    /// the element's Jacobian.
    ElementWithExternalElement::ignore_external_geometric_data();
  }

  ///  Overload get_base_flow_u(...) to return the external
  /// element's velocity components at the integration point
  virtual void get_base_flow_u(const double& time,
                               const unsigned& ipt,
                               const Vector<double>& x,
                               Vector<double>& result) const
  {
    // Set interaction index to 0
    const unsigned interaction = 0;

    // Get a pointer to the external element that computes the base flow.
    // We know that it's an axisymmetric Navier-Stokes element.
    const QTaylorHoodElement<DIM>* base_flow_el_pt =
      dynamic_cast<QTaylorHoodElement<DIM>*>(
        external_element_pt(interaction, ipt));

    // Provide storage for local coordinates in the external element
    // which correspond to the integration point ipt
    Vector<double> s_external(2);

    // Determine local coordinates in the external element which correspond
    // to the integration point ipt
    s_external = external_element_local_coord(interaction, ipt);

    // Get the DIM velocity components interpolated from the external element
    for (unsigned i = 0; i < DIM; i++)
    {
      result[i] = base_flow_el_pt->interpolated_u_nst(s_external, i);
    }

  } // End of overloaded get_base_flow_u function

  ///  Overload get_base_flow_dudx(...) to return the derivatives of
  /// the external element's velocity components w.r.t. global coordinates
  /// at the integration point
  virtual void get_base_flow_dudx(const double& time,
                                  const unsigned& ipt,
                                  const Vector<double>& x,
                                  DenseMatrix<double>& result) const
  {
    // Set interaction index to 1
    const unsigned interaction = 1;

    // Get a pointer to the external element that computes the base flow.
    // We know that it's an axisymmetric Navier-Stokes element.
    const QTaylorHoodElement<DIM>* base_flow_el_pt =
      dynamic_cast<QTaylorHoodElement<DIM>*>(
        external_element_pt(interaction, ipt));

    // Provide storage for local coordinates in the external element
    // which correspond to the integration point ipt
    Vector<double> s_external(2);

    // Determine local coordinates in the external element which correspond
    // to the integration point ipt
    s_external = external_element_local_coord(interaction, ipt);

    // Loop over velocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Loop over coordinate directions and get derivatives of velocity
      // components from the external element
      for (unsigned j = 0; j < DIM; j++)
      {
        result(i, j) = base_flow_el_pt->interpolated_dudx_nst(s_external, i, j);
      }
    }

  } // End of overloaded get_base_flow_dudx function


  ///  Compute the element's residual vector and the Jacobian matrix
  void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                        DenseMatrix<double>& jacobian)
  {
    // Get the analytical contribution from the basic linearised element
    LinearisedQTaylorHoodElement::fill_in_contribution_to_jacobian(residuals,
                                                                   jacobian);

    // Get the off-diagonal terms by finite differencing
    this->fill_in_jacobian_from_external_interaction_by_fd(residuals, jacobian);
  }
};


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


//======================================================================
/// Build a LinearisedQCrouzeixRaviart element that inherits
/// from ElementWithExternalElement so that it can "communicate" with an
/// axisymmetric Navier-Stokes element that provides the base flow
/// velocities and their derivatives w.r.t. global coordinates (r and z)
//======================================================================
class LinearisedQCrouzeixRaviartMultiDomainElement
  : public virtual LinearisedQCrouzeixRaviartElement,
    public virtual ElementWithExternalElement
{
public:
  /// Constructor: call the underlying constructors
  LinearisedQCrouzeixRaviartMultiDomainElement()
    : LinearisedQCrouzeixRaviartElement(), ElementWithExternalElement()
  {
    // There are two interactions: the base flow velocities and their
    // derivatives w.r.t. global coordinates
    this->set_ninteraction(2);

    // Do not include any external interaction data when computing
    // the element's Jacobian
    // ElementWithExternalElement::ignore_external_interaction_data();

    /// Do not include any external geometric data when computing
    /// the element's Jacobian.
    ElementWithExternalElement::ignore_external_geometric_data();
  }

  ///  Overload get_base_flow_u(...) to return the external
  /// element's velocity components at the integration point
  virtual void get_base_flow_u(const double& time,
                               const unsigned& ipt,
                               const Vector<double>& x,
                               Vector<double>& result) const
  {
    // Set interaction index to 0
    const unsigned interaction = 0;

    // Get a pointer to the external element that computes the base flow.
    // We know that it's an axisymmetric Navier-Stokes element.
    const QCrouzeixRaviartElement<DIM>* base_flow_el_pt =
      dynamic_cast<QCrouzeixRaviartElement<DIM>*>(
        external_element_pt(interaction, ipt));

    // Provide storage for local coordinates in the external element
    // which correspond to the integration point ipt
    Vector<double> s_external(2);

    // Determine local coordinates in the external element which correspond
    // to the integration point ipt
    s_external = external_element_local_coord(interaction, ipt);

    // Get the DIM velocity components interpolated from the external element
    for (unsigned i = 0; i < DIM; i++)
    {
      result[i] = base_flow_el_pt->interpolated_u_nst(s_external, i);
    }

  } // End of overloaded get_base_flow_u function

  ///  Overload get_base_flow_dudx(...) to return the derivatives of
  /// the external element's velocity components w.r.t. global coordinates
  /// at the integration point
  virtual void get_base_flow_dudx(const double& time,
                                  const unsigned& ipt,
                                  const Vector<double>& x,
                                  DenseMatrix<double>& result) const
  {
    // Set interaction index to 1
    const unsigned interaction = 1;

    // Get a pointer to the external element that computes the base flow.
    // We know that it's an axisymmetric Navier-Stokes element.
    const QCrouzeixRaviartElement<DIM>* base_flow_el_pt =
      dynamic_cast<QCrouzeixRaviartElement<DIM>*>(
        external_element_pt(interaction, ipt));

    // Provide storage for local coordinates in the external element
    // which correspond to the integration point ipt
    Vector<double> s_external(2);

    // Determine local coordinates in the external element which correspond
    // to the integration point ipt
    s_external = external_element_local_coord(interaction, ipt);

    // Loop over velocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Loop over coordinate directions and get derivatives of velocity
      // components from the external element
      for (unsigned j = 0; j < DIM; j++)
      {
        result(i, j) = base_flow_el_pt->interpolated_dudx_nst(s_external, i, j);
      }
    }

  } // End of overloaded get_base_flow_dudx function


  ///  Compute the element's residual vector and the Jacobian matrix
  void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                        DenseMatrix<double>& jacobian)
  {
    // Get the analytical contribution from the basic linearised element
    LinearisedQCrouzeixRaviartElement::fill_in_contribution_to_jacobian(
      residuals, jacobian);

    // Get the off-diagonal terms by finite differencing
    this->fill_in_jacobian_from_external_interaction_by_fd(residuals, jacobian);
  }
};


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


//======================================================================
///  Build a RefineableLinearisedQTaylorHood element
/// that inherits from ElementWithExternalElement so that it can
/// "communicate" with an axisymmetric Navier-Stokes element that
/// provides the base flow velocities and their derivatives w.r.t.
/// global coordinates (r and z)
//======================================================================
class RefineableLinearisedQTaylorHoodMultiDomainElement
  : public virtual RefineableLinearisedQTaylorHoodElement,
    public virtual ElementWithExternalElement
{
public:
  /// Constructor: call the underlying constructors
  RefineableLinearisedQTaylorHoodMultiDomainElement()
    : RefineableLinearisedQTaylorHoodElement(), ElementWithExternalElement()
  {
    // There are two interactions: the base flow velocities and their
    // derivatives w.r.t. global coordinates
    this->set_ninteraction(2);

    // Do not include any external interaction data when computing
    // the element's Jacobian
    // ElementWithExternalElement::ignore_external_interaction_data();

    /// Do not include any external geometric data when computing
    /// the element's Jacobian.
    ElementWithExternalElement::ignore_external_geometric_data();
  }

  ///  Overload get_base_flow_u(...) to return the external
  /// element's velocity components at the integration point
  virtual void get_base_flow_u(const double& time,
                               const unsigned& ipt,
                               const Vector<double>& x,
                               Vector<double>& result) const
  {
    // Set interaction index to 0
    const unsigned interaction = 0;

    // Get a pointer to the external element that computes the base flow.
    // We know that it's an axisymmetric Navier-Stokes element.
    const QTaylorHoodElement<DIM>* base_flow_el_pt =
      dynamic_cast<QTaylorHoodElement<DIM>*>(
        external_element_pt(interaction, ipt));

    // Provide storage for local coordinates in the external element
    // which correspond to the integration point ipt
    Vector<double> s_external(2);

    // Determine local coordinates in the external element which correspond
    // to the integration point ipt
    s_external = external_element_local_coord(interaction, ipt);

    // Get the three velocity components interpolated from the external element
    for (unsigned i = 0; i < DIM; i++)
    {
      result[i] = base_flow_el_pt->interpolated_u_nst(s_external, i);
    }

  } // End of overloaded get_base_flow_u function

  ///  Overload get_base_flow_dudx(...) to return the derivatives of
  /// the external element's velocity components w.r.t. global coordinates
  /// at the integration point
  virtual void get_base_flow_dudx(const double& time,
                                  const unsigned& ipt,
                                  const Vector<double>& x,
                                  DenseMatrix<double>& result) const
  {
    // Set interaction index to 1
    const unsigned interaction = 1;

    // Get a pointer to the external element that computes the base flow.
    // We know that it's an axisymmetric Navier-Stokes element.
    const QTaylorHoodElement<DIM>* base_flow_el_pt =
      dynamic_cast<QTaylorHoodElement<DIM>*>(
        external_element_pt(interaction, ipt));

    // Provide storage for local coordinates in the external element
    // which correspond to the integration point ipt
    Vector<double> s_external(2);

    // Determine local coordinates in the external element which correspond
    // to the integration point ipt
    s_external = external_element_local_coord(interaction, ipt);

    // Loop over velocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Loop over coordinate directions and get derivatives of velocity
      // components from the external element
      for (unsigned j = 0; j < DIM; j++)
      {
        result(i, j) = base_flow_el_pt->interpolated_dudx_nst(s_external, i, j);
      }
    }

  } // End of overloaded get_base_flow_dudx function


  ///  Compute the element's residual vector and the Jacobian matrix
  void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                        DenseMatrix<double>& jacobian)
  {
    // Get the analytical contribution from the basic patricklinearised element
    RefineableLinearisedQTaylorHoodElement::fill_in_contribution_to_jacobian(
      residuals, jacobian);

    // Get the off-diagonal terms by finite differencing
    this->fill_in_jacobian_from_external_interaction_by_fd(residuals, jacobian);
  }
};


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


//======================================================================
///  Build a RefineableLinearisedQCrouzeixRaviart element
/// that inherits from ElementWithExternalElement so that it can
/// "communicate" with an axisymmetric Navier-Stokes element that
/// provides the base flow velocities and their derivatives w.r.t.
/// global coordinates (r and z)
//======================================================================
class RefineableLinearisedQCrouzeixRaviartMultiDomainElement
  : public virtual RefineableLinearisedQCrouzeixRaviartElement,
    public virtual ElementWithExternalElement
{
public:
  /// Constructor: call the underlying constructors
  RefineableLinearisedQCrouzeixRaviartMultiDomainElement()
    : RefineableLinearisedQCrouzeixRaviartElement(),
      ElementWithExternalElement()
  {
    // There are two interactions: the base flow velocities and their
    // derivatives w.r.t. global coordinates
    this->set_ninteraction(2);

    // Do not include any external interaction data when computing
    // the element's Jacobian
    // ElementWithExternalElement::ignore_external_interaction_data();

    /// Do not include any external geometric data when computing
    /// the element's Jacobian.
    ElementWithExternalElement::ignore_external_geometric_data();
  }

  ///  Overload get_base_flow_u(...) to return the external
  /// element's velocity components at the integration point
  virtual void get_base_flow_u(const double& time,
                               const unsigned& ipt,
                               const Vector<double>& x,
                               Vector<double>& result) const
  {
    // Set interaction index to 0
    const unsigned interaction = 0;

    // Get a pointer to the external element that computes the base flow.
    // We know that it's an axisymmetric Navier-Stokes element.
    const QCrouzeixRaviartElement<DIM>* base_flow_el_pt =
      dynamic_cast<QCrouzeixRaviartElement<DIM>*>(
        external_element_pt(interaction, ipt));

    // Provide storage for local coordinates in the external element
    // which correspond to the integration point ipt
    Vector<double> s_external(2);

    // Determine local coordinates in the external element which correspond
    // to the integration point ipt
    s_external = external_element_local_coord(interaction, ipt);

    // Get the three velocity components interpolated from the external element
    for (unsigned i = 0; i < DIM; i++)
    {
      result[i] = base_flow_el_pt->interpolated_u_nst(s_external, i);
    }

  } // End of overloaded get_base_flow_u function

  ///  Overload get_base_flow_dudx(...) to return the derivatives of
  /// the external element's velocity components w.r.t. global coordinates
  /// at the integration point
  virtual void get_base_flow_dudx(const double& time,
                                  const unsigned& ipt,
                                  const Vector<double>& x,
                                  DenseMatrix<double>& result) const
  {
    // Set interaction index to 1
    const unsigned interaction = 1;

    // Get a pointer to the external element that computes the base flow.
    // We know that it's an axisymmetric Navier-Stokes element.
    const QCrouzeixRaviartElement<DIM>* base_flow_el_pt =
      dynamic_cast<QCrouzeixRaviartElement<DIM>*>(
        external_element_pt(interaction, ipt));

    // Provide storage for local coordinates in the external element
    // which correspond to the integration point ipt
    Vector<double> s_external(2);

    // Determine local coordinates in the external element which correspond
    // to the integration point ipt
    s_external = external_element_local_coord(interaction, ipt);

    // Loop over velocity components
    for (unsigned i = 0; i < DIM; i++)
    {
      // Loop over coordinate directions and get derivatives of velocity
      // components from the external element
      for (unsigned j = 0; j < DIM; j++)
      {
        result(i, j) = base_flow_el_pt->interpolated_dudx_nst(s_external, i, j);
      }
    }

  } // End of overloaded get_base_flow_dudx function


  ///  Compute the element's residual vector and the Jacobian matrix
  void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                        DenseMatrix<double>& jacobian)
  {
    // Get the analytical contribution from the basic patricklinearised element
    RefineableLinearisedQCrouzeixRaviartElement::
      fill_in_contribution_to_jacobian(residuals, jacobian);

    // Get the off-diagonal terms by finite differencing
    this->fill_in_jacobian_from_external_interaction_by_fd(residuals, jacobian);
  }
};
