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
// The element-independent guts for imposition of "constant volume"
// constraints in free surface/interface problems.


#include "constrained_volume_elements.h"

namespace oomph
{
  //=====================================================================
  /// Fill in the residuals for the volume constraint
  //====================================================================
  void VolumeConstraintElement::
    fill_in_generic_contribution_to_residuals_volume_constraint(
      Vector<double>& residuals)
  {
    // Note: This element can only be used with the associated
    // VolumeConstraintBoundingElement elements which compute the actual
    // enclosed volume; here we only add the contribution to the
    // residual; everything else, incl. the derivatives of this
    // residual w.r.t. the nodal positions of the
    // VolumeConstraintBoundingElements
    // is handled by them
    const int local_eqn = this->ptraded_local_eqn();
    if (local_eqn >= 0)
    {
      residuals[local_eqn] -= *Prescribed_volume_pt;
    }
  }

  //===========================================================================
  /// Constructor: Pass pointer to target volume. "Pressure" value that
  /// "traded" for the volume contraint is created internally (as a Data
  /// item with a single pressure value)
  //===========================================================================
  VolumeConstraintElement::VolumeConstraintElement(double* prescribed_volume_pt)
  {
    // Store pointer to prescribed volume
    Prescribed_volume_pt = prescribed_volume_pt;

    // Create data, add as internal data and record the index
    // (gets deleted automatically in destructor of GeneralisedElement)
    External_or_internal_data_index_of_traded_pressure =
      add_internal_data(new Data(1));

    // Record that we've created it as internal data
    Traded_pressure_stored_as_internal_data = true;

    // ...and stored the "traded pressure" value as first value
    Index_of_traded_pressure_value = 0;
  }

  //======================================================================
  /// Constructor: Pass pointer to target volume, pointer to Data
  /// item whose value specified by index_of_traded_pressure represents
  /// the "Pressure" value that "traded" for the volume contraint.
  /// The Data is stored as external Data for this element.
  //======================================================================
  VolumeConstraintElement::VolumeConstraintElement(
    double* prescribed_volume_pt,
    Data* p_traded_data_pt,
    const unsigned& index_of_traded_pressure)
  {
    // Store pointer to prescribed volume
    Prescribed_volume_pt = prescribed_volume_pt;

    // Add as external data and record the index
    External_or_internal_data_index_of_traded_pressure =
      add_external_data(p_traded_data_pt);

    // Record that it is external data
    Traded_pressure_stored_as_internal_data = false;

    // Record index
    Index_of_traded_pressure_value = index_of_traded_pressure;
  }


  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////


  //====================================================================
  /// Helper function to fill in contributions to residuals
  /// (remember that part of the residual is added by the
  /// the associated VolumeConstraintElement). This is specific for
  /// 1D line elements that bound 2D cartesian fluid elements.
  //====================================================================
  void LineVolumeConstraintBoundingElement::
    fill_in_generic_residual_contribution_volume_constraint(
      Vector<double>& residuals)
  {
    // Add in the volume constraint term if required
    const int local_eqn = this->ptraded_local_eqn();
    if (local_eqn >= 0)
    {
      // Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // Set up memeory for the shape functions
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);

      // Set the value of n_intpt
      const unsigned n_intpt = this->integral_pt()->nweight();

      // Storage for the local coordinate
      Vector<double> s(1);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the local coordinate at the integration point
        s[0] = this->integral_pt()->knot(ipt, 0);

        // Get the integral weight
        double W = this->integral_pt()->weight(ipt);

        // Call the derivatives of the shape function at the knot point
        this->dshape_local_at_knot(ipt, psif, dpsifds);

        // Get position and tangent vector
        Vector<double> interpolated_t1(2, 0.0);
        Vector<double> interpolated_x(2, 0.0);
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directional components
          for (unsigned i = 0; i < 2; i++)
          {
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
            interpolated_t1[i] += this->nodal_position(l, i) * dpsifds(l, 0);
          }
        }

        // Calculate the length of the tangent Vector
        double tlength = interpolated_t1[0] * interpolated_t1[0] +
                         interpolated_t1[1] * interpolated_t1[1];

        // Set the Jacobian of the line element
        double J = sqrt(tlength);

        // Now calculate the normal Vector
        Vector<double> interpolated_n(2);
        this->outer_unit_normal(ipt, interpolated_n);

        // Assemble dot product
        double dot = 0.0;
        for (unsigned k = 0; k < 2; k++)
        {
          dot += interpolated_x[k] * interpolated_n[k];
        }

        // Add to residual with sign chosen so that the volume is
        // positive when the elements bound the fluid
        residuals[local_eqn] += 0.5 * dot * W * J;
      }
    }
  }


  //====================================================================
  /// Return this element's contribution to the total volume enclosed
  //====================================================================
  double LineVolumeConstraintBoundingElement::contribution_to_enclosed_volume()
  {
    // Initialise
    double vol = 0.0;

    // Find out how many nodes there are
    const unsigned n_node = this->nnode();

    // Set up memeory for the shape functions
    Shape psif(n_node);
    DShape dpsifds(n_node, 1);

    // Set the value of n_intpt
    const unsigned n_intpt = this->integral_pt()->nweight();

    // Storage for the local cooridinate
    Vector<double> s(1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the local coordinate at the integration point
      s[0] = this->integral_pt()->knot(ipt, 0);

      // Get the integral weight
      double W = this->integral_pt()->weight(ipt);

      // Call the derivatives of the shape function at the knot point
      this->dshape_local_at_knot(ipt, psif, dpsifds);

      // Get position and tangent vector
      Vector<double> interpolated_t1(2, 0.0);
      Vector<double> interpolated_x(2, 0.0);
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directional components
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += this->nodal_position(l, i) * psif(l);
          interpolated_t1[i] += this->nodal_position(l, i) * dpsifds(l, 0);
        }
      }

      // Calculate the length of the tangent Vector
      double tlength = interpolated_t1[0] * interpolated_t1[0] +
                       interpolated_t1[1] * interpolated_t1[1];

      // Set the Jacobian of the line element
      double J = sqrt(tlength);

      // Now calculate the normal Vector
      Vector<double> interpolated_n(2);
      this->outer_unit_normal(ipt, interpolated_n);

      // Assemble dot product
      double dot = 0.0;
      for (unsigned k = 0; k < 2; k++)
      {
        dot += interpolated_x[k] * interpolated_n[k];
      }

      // Add to volume with sign chosen so that the volume is
      // positive when the elements bound the fluid
      vol += 0.5 * dot * W * J;
    }

    return vol;
  }


  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////

  //====================================================================
  /// Return this element's contribution to the total volume enclosed
  //====================================================================
  double AxisymmetricVolumeConstraintBoundingElement::
    contribution_to_enclosed_volume()
  {
    // Initialise
    double vol = 0.0;

    // Find out how many nodes there are
    const unsigned n_node = this->nnode();

    // Set up memeory for the shape functions
    Shape psif(n_node);
    DShape dpsifds(n_node, 1);

    // Set the value of n_intpt
    const unsigned n_intpt = this->integral_pt()->nweight();

    // Storage for the local cooridinate
    Vector<double> s(1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the local coordinate at the integration point
      s[0] = this->integral_pt()->knot(ipt, 0);

      // Get the integral weight
      double W = this->integral_pt()->weight(ipt);

      // Call the derivatives of the shape function at the knot point
      this->dshape_local_at_knot(ipt, psif, dpsifds);

      // Get position and tangent vector
      Vector<double> interpolated_t1(2, 0.0);
      Vector<double> interpolated_x(2, 0.0);
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directional components
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += this->nodal_position(l, i) * psif(l);
          interpolated_t1[i] += this->nodal_position(l, i) * dpsifds(l, 0);
        }
      }

      // Calculate the length of the tangent Vector
      double tlength = interpolated_t1[0] * interpolated_t1[0] +
                       interpolated_t1[1] * interpolated_t1[1];

      // Set the Jacobian of the line element
      double J = sqrt(tlength) * interpolated_x[0];

      // Now calculate the normal Vector
      Vector<double> interpolated_n(2);
      this->outer_unit_normal(ipt, interpolated_n);

      // Assemble dot product
      double dot = 0.0;
      for (unsigned k = 0; k < 2; k++)
      {
        dot += interpolated_x[k] * interpolated_n[k];
      }

      // Add to volume with sign chosen so that the volume is
      // positive when the elements bound the fluid
      vol += dot * W * J / 3.0;
    }

    return vol;
  }


  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //====================================================================
  /// Helper function to fill in contributions to residuals
  /// (remember that part of the residual is added by the
  /// the associated VolumeConstraintElement). This is specific for
  /// axisymmetric line elements that bound 2D axisymmetric fluid elements.
  /// The only difference from the 1D case is the addition of the
  /// weighting factor in the integrand and division of the dot product
  /// by three.
  //====================================================================
  void AxisymmetricVolumeConstraintBoundingElement::
    fill_in_generic_residual_contribution_volume_constraint(
      Vector<double>& residuals)
  {
    // Add in the volume constraint term if required
    const int local_eqn = this->ptraded_local_eqn();
    if (local_eqn >= 0)
    {
      // Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // Set up memeory for the shape functions
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);

      // Set the value of n_intpt
      const unsigned n_intpt = this->integral_pt()->nweight();

      // Storage for the local cooridinate
      Vector<double> s(1);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the local coordinate at the integration point
        s[0] = this->integral_pt()->knot(ipt, 0);

        // Get the integral weight
        double W = this->integral_pt()->weight(ipt);

        // Call the derivatives of the shape function at the knot point
        this->dshape_local_at_knot(ipt, psif, dpsifds);

        // Get position and tangent vector
        Vector<double> interpolated_t1(2, 0.0);
        Vector<double> interpolated_x(2, 0.0);
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directional components
          for (unsigned i = 0; i < 2; i++)
          {
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
            interpolated_t1[i] += this->nodal_position(l, i) * dpsifds(l, 0);
          }
        }

        // Calculate the length of the tangent Vector
        double tlength = interpolated_t1[0] * interpolated_t1[0] +
                         interpolated_t1[1] * interpolated_t1[1];

        // Set the Jacobian of the line element
        // multiplied by r (x[0])
        double J = sqrt(tlength) * interpolated_x[0];

        // Now calculate the normal Vector
        Vector<double> interpolated_n(2);
        this->outer_unit_normal(ipt, interpolated_n);

        // Assemble dot product
        double dot = 0.0;
        for (unsigned k = 0; k < 2; k++)
        {
          dot += interpolated_x[k] * interpolated_n[k];
        }

        // Add to residual with sign chosen so that the volume is
        // positive when the elements bound the fluid
        residuals[local_eqn] += dot * W * J / 3.0;
      }
    }
  }

  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //=================================================================
  /// Helper function to fill in contributions to residuals
  /// (remember that part of the residual is added by the
  /// the associated VolumeConstraintElement). This is specific for
  /// 2D surface elements that bound 3D cartesian fluid elements.
  //=================================================================
  void SurfaceVolumeConstraintBoundingElement::
    fill_in_generic_residual_contribution_volume_constraint(
      Vector<double>& residuals)
  {
    // Add in the volume constraint term if required
    const int local_eqn = this->ptraded_local_eqn();
    if (local_eqn >= 0)
    {
      // Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // Set up memeory for the shape functions
      Shape psif(n_node);
      DShape dpsifds(n_node, 2);

      // Set the value of n_intpt
      const unsigned n_intpt = this->integral_pt()->nweight();

      // Storage for the local coordinate
      Vector<double> s(2);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the local coordinate at the integration point
        for (unsigned i = 0; i < 2; i++)
        {
          s[i] = this->integral_pt()->knot(ipt, i);
        }

        // Get the integral weight
        double W = this->integral_pt()->weight(ipt);

        // Call the derivatives of the shape function at the knot point
        this->dshape_local_at_knot(ipt, psif, dpsifds);

        // Get position and tangent vector
        DenseMatrix<double> interpolated_g(2, 3, 0.0);
        Vector<double> interpolated_x(3, 0.0);
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directional components
          for (unsigned i = 0; i < 3; i++)
          {
            // Cache the nodal position
            const double x_ = this->nodal_position(l, i);
            // Calculate the position
            interpolated_x[i] += x_ * psif(l);
            // Calculate the two tangent vecors
            for (unsigned j = 0; j < 2; j++)
            {
              interpolated_g(j, i) += x_ * dpsifds(l, j);
            }
          }
        }

        // Define the (non-unit) normal vector (cross product of tangent
        // vectors)
        Vector<double> interpolated_n(3);
        interpolated_n[0] = interpolated_g(0, 1) * interpolated_g(1, 2) -
                            interpolated_g(0, 2) * interpolated_g(1, 1);
        interpolated_n[1] = interpolated_g(0, 2) * interpolated_g(1, 0) -
                            interpolated_g(0, 0) * interpolated_g(1, 2);
        interpolated_n[2] = interpolated_g(0, 0) * interpolated_g(1, 1) -
                            interpolated_g(0, 1) * interpolated_g(1, 0);

        // The determinant of the local metric tensor is
        // equal to the length of the normal vector, so
        // rather than dividing and multipling, we just leave it out
        // completely, which is why there is no J in the sum below

        // We can now find the sign to get the OUTER UNIT normal
        for (unsigned i = 0; i < 3; i++)
        {
          interpolated_n[i] *= this->normal_sign();
        }

        // Assemble dot product
        double dot = 0.0;
        for (unsigned k = 0; k < 3; k++)
        {
          dot += interpolated_x[k] * interpolated_n[k];
        }

        // Add to residual with the sign chosen such that the volume is
        // positive when the faces are assembled such that the outer unit
        // normal is directed out of the bulk fluid
        residuals[local_eqn] += dot * W / 3.0;
      }
    }
  }


} // namespace oomph
