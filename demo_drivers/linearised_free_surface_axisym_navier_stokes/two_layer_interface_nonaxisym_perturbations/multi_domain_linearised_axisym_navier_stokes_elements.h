//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
// Header for an element that couples a linearised axisymmetric
// Navier-Stokes element to a non-linear axisymmetric Navier-Stokes
// element via a multi domain approach

// oomph-lib headers
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "linearised_axisym_navier_stokes_elements.h"
#include "perturbed_spines.h"

// Use the oomph namespace
using namespace oomph;



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//======================================================================
/// Build a LinearisedAxisymmetricQCrouzeixRaviart element that inherits
/// from ElementWithExternalElement so that it can "communicate" with an
/// axisymmetric Navier-Stokes element that provides the base flow
/// velocities and their derivatives w.r.t. global coordinates (r and z)
//======================================================================
class LinearisedAxisymmetricQCrouzeixRaviartMultiDomainElement : 
public virtual LinearisedAxisymmetricQCrouzeixRaviartElement,
 public virtual ElementWithExternalElement
{
  public:
 
 /// Constructor: call the underlying constructors
 LinearisedAxisymmetricQCrouzeixRaviartMultiDomainElement() :
  LinearisedAxisymmetricQCrouzeixRaviartElement(),
  ElementWithExternalElement()
   {
    // There are three interactions: the base flow velocities, their 
    // derivatives w.r.t. global coordinates and the base flow pressure
    this->set_ninteraction(5);
    
    // Do not include any external interaction data when computing
    // the element's Jacobian
    ElementWithExternalElement::ignore_external_interaction_data();

    /// Do not include any external geometric data when computing
    /// the element's Jacobian.
    ElementWithExternalElement::ignore_external_geometric_data();
   }
  
  /// Overload get_base_flow_u(...) to return the external
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
   const AxisymmetricQCrouzeixRaviartElement* base_flow_el_pt =
    dynamic_cast<AxisymmetricQCrouzeixRaviartElement*>(
     external_element_pt(interaction,ipt));

   // Provide storage for local coordinates in the external element
   // which correspond to the integration point ipt
   Vector<double> s_external(2);

   // Determine local coordinates in the external element which correspond
   // to the integration point ipt
   s_external = external_element_local_coord(interaction,ipt);

   // Get the three velocity components interpolated from the external element
   for(unsigned i=0;i<3;i++)
    {
     result[i] = base_flow_el_pt->interpolated_u_axi_nst(s_external,i);
    }

  } // End of overloaded get_base_flow_u function

  /// Overload get_base_flow_dudx(...) to return the derivatives of
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
   const AxisymmetricQCrouzeixRaviartElement* base_flow_el_pt =
    dynamic_cast<AxisymmetricQCrouzeixRaviartElement*>(
     external_element_pt(interaction,ipt));
   
   // Provide storage for local coordinates in the external element
   // which correspond to the integration point ipt
   Vector<double> s_external(2);
   
   // Determine local coordinates in the external element which correspond
   // to the integration point ipt
   s_external = external_element_local_coord(interaction,ipt);
   
   // Loop over velocity components
   for(unsigned i=0;i<3;i++)
    {
     // Loop over coordinate directions and get derivatives of velocity
     // components from the external element
     for(unsigned j=0;j<2;j++)
      {
       result(i,j)
        = base_flow_el_pt->interpolated_dudx_axi_nst(s_external,i,j);
      }
    }

  } // End of overloaded get_base_flow_dudx function

  /// Overload get_base_flow_p(...) to return the external
  /// element's pressure at the integration point
  virtual void get_base_flow_p(const double& time,
                               const unsigned& ipt,
                               const Vector<double>& x,
                               double& result) const
  {
   // Set interaction index to 2
   const unsigned interaction = 2;
   
   // Get a pointer to the external element that computes the base flow.
   // We know that it's an axisymmetric Navier-Stokes element.
   const AxisymmetricQCrouzeixRaviartElement* base_flow_el_pt =
    dynamic_cast<AxisymmetricQCrouzeixRaviartElement*>(
     external_element_pt(interaction,ipt));

   // Provide storage for local coordinates in the external element
   // which correspond to the integration point ipt
   Vector<double> s_external(2);

   // Determine local coordinates in the external element which correspond
   // to the integration point ipt
   s_external = external_element_local_coord(interaction,ipt);

   // Get the pressure interpolated from the external element
   result = base_flow_el_pt->interpolated_p_axi_nst(s_external);

  } // End of overloaded get_base_flow_p function

  /// Overload get_base_flow_dudt(...) to return the derivatives of
  /// the external element's velocity components w.r.t. time
  /// at the integration point
  virtual void get_base_flow_dudt(const double& time,
                                  const unsigned& ipt,
                                  const Vector<double>& x,
                                  Vector<double>& result) const
  {
   // Set interaction index to 3
   const unsigned interaction = 3;
   
   // Get a pointer to the external element that computes the base flow.
   // We know that it's an axisymmetric Navier-Stokes element.
   const AxisymmetricQCrouzeixRaviartElement* base_flow_el_pt =
    dynamic_cast<AxisymmetricQCrouzeixRaviartElement*>(
     external_element_pt(interaction,ipt));
   
   // Provide storage for local coordinates in the external element
   // which correspond to the integration point ipt
   Vector<double> s_external(2);
   
   // Determine local coordinates in the external element which correspond
   // to the integration point ipt
   s_external = external_element_local_coord(interaction,ipt);
   
   // Loop over velocity components
   for(unsigned i=0;i<3;i++)
    {
     // Get dudt components from the external element
     result[i] = base_flow_el_pt->interpolated_dudt_axi_nst(s_external,i);
    }

  } // End of overloaded get_base_flow_dudt function

  /// Overload get_base_flow_duds(...) to return the derivatives of
  /// the external element's velocity components w.r.t. local coordinates
  /// at the integration point
  virtual void get_base_flow_duds(const double& time,
                                  const unsigned& ipt,
                                  const Vector<double>& x,
                                  DenseMatrix<double>& result) const
  {
   // Set interaction index to 4
   const unsigned interaction = 4;
   
   // Get a pointer to the external element that computes the base flow.
   // We know that it's an axisymmetric Navier-Stokes element.
   const AxisymmetricQCrouzeixRaviartElement* base_flow_el_pt =
    dynamic_cast<AxisymmetricQCrouzeixRaviartElement*>(
     external_element_pt(interaction,ipt));
   
   // Provide storage for local coordinates in the external element
   // which correspond to the integration point ipt
   Vector<double> s_external(2);
   
   // Determine local coordinates in the external element which correspond
   // to the integration point ipt
   s_external = external_element_local_coord(interaction,ipt);
   
   // Loop over velocity components
   for(unsigned i=0;i<3;i++)
    {
     // Loop over coordinate directions and get derivatives of velocity
     // components from the external element
     for(unsigned j=0;j<2;j++)
      {
       result(i,j)
        = base_flow_el_pt->interpolated_duds_axi_nst(s_external,i,j);
      }
    }

  } // End of overloaded get_base_flow_duds function

  /// Compute the element's residual vector and the Jacobian matrix
  void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                        DenseMatrix<double> &jacobian)
  {
   // Get the analytical contribution from the basic patricklinearised element
   LinearisedAxisymmetricQCrouzeixRaviartElement::
    fill_in_contribution_to_jacobian(residuals,jacobian);
  }
  
};



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



namespace oomph
{

 //==========================================================================
 /// Face geometry of the linearised axisymmetric CrouzeixRaviart elements
 //==========================================================================
 template<>
  class FaceGeometry<LinearisedAxisymmetricQCrouzeixRaviartMultiDomainElement>
  : public virtual QElement<1,3>
  {
    public:
   FaceGeometry() : QElement<1,3>() {}
  };

}
