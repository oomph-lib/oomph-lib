// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
#include "young_laplace_contact_angle_elements.h"
#include "young_laplace_elements.h"
#include "refineable_young_laplace_elements.h"


namespace oomph
{
  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element  and the
  /// index of the face to which the element is attached.
  //===========================================================================
  template<class ELEMENT>
  YoungLaplaceContactAngleElement<ELEMENT>::YoungLaplaceContactAngleElement(
    FiniteElement* const& bulk_el_pt, const int& face_index)
    : FaceGeometry<ELEMENT>(), FaceElement()
  {
    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index, this);

    // Initialise the prescribed contact angle pointer to zero
    Prescribed_cos_gamma_pt = 0;

#ifdef PARANOID
    // Extract the dimension of the problem from the dimension of
    // the first node
    unsigned dim = node_pt(0)->ndim();
    if (dim != 2)
    {
      throw OomphLibError(
        "YoungLaplaceContactAngleElement only work with 2D nodes",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif
  }


  //================================================================
  /// Compute the element's contribution to the residual vector
  //================================================================
  template<class ELEMENT>
  void YoungLaplaceContactAngleElement<
    ELEMENT>::fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Set up memory for the shape functions
    Shape psi(n_node);

    // Number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    // Note: We need the coordinate itself below for the evaluation
    // of the contact line vectors even though we use the *at_knot
    // version for the various shape-function-related functions
    Vector<double> s(1);

    // Integers to hold the local equation and unknown numbers
    int local_eqn = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign value of s
      s[0] = integral_pt()->knot(ipt, 0);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Find the shape functions
      shape_at_knot(ipt, psi);

      // Get the prescribed cos_gamma
      double cos_gamma = prescribed_cos_gamma();

      // Get the various contact line vectors
      Vector<double> tangent(3);
      Vector<double> normal(3);
      Vector<double> spine(3);
      double norm_of_drds;
      contact_line_vectors(s, tangent, normal, spine, norm_of_drds);

      // Get beta factor:

      // Cross product of spine and tangent to contact line is
      // the wall normal
      Vector<double> wall_normal(3);
      ELEMENT::cross_product(spine, tangent, wall_normal);

      // Normalise
      double norm = ELEMENT::two_norm(wall_normal);
      for (unsigned i = 0; i < 3; i++) wall_normal[i] /= norm;

      // Take cross product with tangent to get the normal to the
      // contact line parallel to wall
      Vector<double> normal_to_contact_line_parallel_to_wall(3);
      ELEMENT::cross_product(
        tangent, wall_normal, normal_to_contact_line_parallel_to_wall);

      double beta =
        ELEMENT::scalar_product(spine, normal_to_contact_line_parallel_to_wall);


      // Now add to the appropriate equations

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        local_eqn = u_local_eqn(l);
        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          // Add to residual:
          residuals[local_eqn] -= beta * cos_gamma * psi[l] * norm_of_drds * w;
        }
      }
    }
  }


  //========================================================================
  /// Get the actual contact angle
  //========================================================================
  template<class ELEMENT>
  double YoungLaplaceContactAngleElement<ELEMENT>::actual_cos_contact_angle(
    const Vector<double>& s)
  {
    // Get pointer to bulk element
    ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

    double cos_gamma = 0.0;

    // Spine case
    if (bulk_elem_pt->use_spines())
    {
      // Get various contact line vectors
      Vector<double> tangent(3);
      Vector<double> normal(3);
      Vector<double> spine(3);
      double norm_of_drds;
      contact_line_vectors(s, tangent, normal, spine, norm_of_drds);

      // Get the wall normal: Both the tangent to the contact line
      // and the spine vector are tangential to the wall:
      Vector<double> wall_normal(3);
      Vector<double> axe_ez(3, 0.0);
      axe_ez[2] = 1.0;
      ELEMENT::cross_product(axe_ez, tangent, wall_normal);

      // Normalise
      double norm = 0.0;
      for (unsigned i = 0; i < 3; i++) norm += wall_normal[i] * wall_normal[i];
      for (unsigned i = 0; i < 3; i++)
      {
        wall_normal[i] /= sqrt(norm);
      }

      // Get the auxiliary unit vector that's normal to
      // the contact line and tangent to the wall
      Vector<double> aux(3);
      ELEMENT::cross_product(tangent, wall_normal, aux);

      // Cos of contact angle is dot product with wall normal
      cos_gamma = ELEMENT::scalar_product(aux, normal);
    }

    // Cartesian case
    else
    {
      // Get local coordinates in bulk element by copy construction
      Vector<double> s_bulk(local_coordinate_in_bulk(s));

      // Number of nodes in bulk element
      unsigned nnode_bulk = bulk_elem_pt->nnode();


#ifdef PARANOID
      // Dimension of (= number of local coordinates in) bulk element
      unsigned dim_bulk = bulk_elem_pt->dim();

      if (dim_bulk != 2)
      {
        throw OomphLibError(
          "YoungLaplaceContactAngleElements only work with 2D bulk elements",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Set up memory for the shape functions
      Shape psi(nnode_bulk);
      DShape dpsidzeta(nnode_bulk, 2);

      // Call the derivatives of the shape and test functions
      // in the bulk -- must do this via s because this point
      // is not an integration point the bulk element!
      (void)bulk_elem_pt->dshape_eulerian(s_bulk, psi, dpsidzeta);

      // Get the gradient at s
      Vector<double> gradient_u(2, 0.0);

      // Calculate function value and derivatives:
      //-----------------------------------------
      // Loop over nodes
      for (unsigned l = 0; l < nnode_bulk; l++)
      {
        // Loop over directions
        for (unsigned j = 0; j < 2; j++)
        {
          gradient_u[j] += bulk_elem_pt->u(l) * dpsidzeta(l, j);
        }
      }

      // Get the outer unit normal to boundary
      Vector<double> outer_normal(2, 0.0);
      outer_unit_normal(s, outer_normal);

      // Compute the cosinus of the angle
      double gradient_norm_2 =
        ELEMENT::two_norm(gradient_u) * ELEMENT::two_norm(gradient_u);
      cos_gamma = ELEMENT::scalar_product(gradient_u, outer_normal) /
                  sqrt(1 + gradient_norm_2);
    }

    return cos_gamma;
  }


  //========================================================================
  /// Get unit tangent and normal to contact line and the spine itself (this
  /// allows the wall normal to be worked out by a cross product). Final
  /// argument gives the norm of dR/ds where R is the vector to the
  /// contact line and s the local coordinate in the element.
  //========================================================================
  template<class ELEMENT>
  void YoungLaplaceContactAngleElement<ELEMENT>::contact_line_vectors(
    const Vector<double>& s,
    Vector<double>& tangent,
    Vector<double>& normal,
    Vector<double>& spine,
    double& norm_of_drds)
  {
    // Get pointer to bulk element
    ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

    // Dimension of (= number of local coordinates in) bulk element
    unsigned dim_bulk = bulk_elem_pt->dim();

#ifdef PARANOID
    if (dim_bulk != 2)
    {
      throw OomphLibError(
        "YoungLaplaceContactAngleElements only work with 2D bulk elements",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Which local coordinate is const in bulk element?
    unsigned s_fixed_index_in_bulk;

    DenseMatrix<double> ds_bulk_ds_face(2, 1);
    this->get_ds_bulk_ds_face(s, ds_bulk_ds_face, s_fixed_index_in_bulk);

    // Get local coordinates in bulk element by copy construction
    Vector<double> s_bulk(local_coordinate_in_bulk(s));

    // Number of nodes in bulk element
    unsigned nnode_bulk = bulk_elem_pt->nnode();

    // Get the shape functions and their derivatives w.r.t.
    // to the local coordinates in the bulk element
    Shape psi_bulk(nnode_bulk);
    DShape dpsi_bulk(nnode_bulk, dim_bulk);
    bulk_elem_pt->dshape_local(s_bulk, psi_bulk, dpsi_bulk);

    // Displacement along spine
    double interpolated_u = 0.0;

    // Derivative of u w.r.t. tangential and pseudo-normal bulk coordinate
    double interpolated_du_ds_tangent = 0.0;
    double interpolated_du_ds_pseudo_normal = 0.0;

    // Global intrinsic coordinates of current point for evaluation of
    // spine and spine base
    Vector<double> interpolated_zeta(dim_bulk, 0.0);


    // Derivatives of zeta (in bulk) w.r.t. to tangent and pseudo-normal
    // local coordinate
    Vector<double> interpolated_dzeta_ds_tangent(dim_bulk);
    Vector<double> interpolated_dzeta_ds_pseudo_normal(dim_bulk);

    // Loop over nodes
    for (unsigned l = 0; l < nnode_bulk; l++)
    {
      interpolated_u += bulk_elem_pt->u(l) * psi_bulk(l);

      // Chain rule for tangent derivative
      double aux = 0.0;
      for (unsigned i = 0; i < dim_bulk; i++)
      {
        aux += dpsi_bulk(l, i) * ds_bulk_ds_face(i, 0);
      }
      interpolated_du_ds_tangent += bulk_elem_pt->u(l) * aux;

      // Straight derivative w.r.t. one of the bulk coordinates
      // that's fixed along the face
      interpolated_du_ds_pseudo_normal +=
        bulk_elem_pt->u(l) * dpsi_bulk(l, s_fixed_index_in_bulk);

      // Loop over directions
      for (unsigned j = 0; j < dim_bulk; j++)
      {
        interpolated_zeta[j] +=
          bulk_elem_pt->nodal_position(l, j) * psi_bulk(l);

        // Chain rule for tangent derivative
        double aux = 0.0;
        for (unsigned i = 0; i < dim_bulk; i++)
        {
          aux += dpsi_bulk(l, i) * ds_bulk_ds_face(i, 0);
        }
        interpolated_dzeta_ds_tangent[j] +=
          bulk_elem_pt->nodal_position(l, j) * aux;

        // Straight derivative w.r.t. one of the bulk coordinates
        // that's fixed along the face
        interpolated_dzeta_ds_pseudo_normal[j] +=
          bulk_elem_pt->nodal_position(l, j) *
          dpsi_bulk(l, s_fixed_index_in_bulk);
      }
    }

    // Auxiliary vector (tangent to non-fixed bulk coordinate but
    // not necessarily normal to contact line)
    Vector<double> aux_vector(3);
    double tang_norm = 0.0;
    double aux_norm = 0.0;

    if (bulk_elem_pt->use_spines())
    {
      // Get the spine and spine base vector at this point from the bulk element
      Vector<double> spine_base(3, 0.0);
      Vector<Vector<double>> dspine;
      ELEMENT::allocate_vector_of_vectors(2, 3, dspine);
      Vector<Vector<double>> dspine_base;
      ELEMENT::allocate_vector_of_vectors(2, 3, dspine_base);
      bulk_elem_pt->get_spine(interpolated_zeta, spine, dspine);
      bulk_elem_pt->get_spine_base(interpolated_zeta, spine_base, dspine_base);

      // Derivative of spine and spine base w.r.t. local coordinate in
      // FaceElement:
      Vector<double> dspine_ds_tangent(3, 0.0);
      Vector<double> dspine_base_ds_tangent(3, 0.0);
      Vector<double> dspine_ds_pseudo_normal(3, 0.0);
      Vector<double> dspine_base_ds_pseudo_normal(3, 0.0);
      for (unsigned i = 0; i < 3; i++)
      {
        dspine_ds_tangent[i] +=
          dspine[0][i] * interpolated_dzeta_ds_tangent[0] +
          dspine[1][i] * interpolated_dzeta_ds_tangent[1];

        dspine_base_ds_tangent[i] +=
          dspine_base[0][i] * interpolated_dzeta_ds_tangent[0] +
          dspine_base[1][i] * interpolated_dzeta_ds_tangent[1];

        dspine_ds_pseudo_normal[i] +=
          dspine[0][i] * interpolated_dzeta_ds_pseudo_normal[0] +
          dspine[1][i] * interpolated_dzeta_ds_pseudo_normal[1];

        dspine_base_ds_pseudo_normal[i] +=
          dspine_base[0][i] * interpolated_dzeta_ds_pseudo_normal[0] +
          dspine_base[1][i] * interpolated_dzeta_ds_pseudo_normal[1];
      }

      // Auxiliary vector (tangent to non-fixed bulk coordinate but
      // not necessarily normal to contact line)
      for (unsigned i = 0; i < 3; i++)
      {
        tangent[i] = dspine_base_ds_tangent[i] +
                     interpolated_du_ds_tangent * spine[i] +
                     interpolated_u * dspine_ds_tangent[i];
        tang_norm += tangent[i] * tangent[i];


        aux_vector[i] = dspine_base_ds_pseudo_normal[i] +
                        interpolated_du_ds_pseudo_normal * spine[i] +
                        interpolated_u * dspine_ds_pseudo_normal[i];
        aux_norm += aux_vector[i] * aux_vector[i];
      }
    }

    // Cartesian case
    else
    {
      for (unsigned i = 0; i < 2; i++)
      {
        tangent[i] = interpolated_dzeta_ds_tangent[i];
        tang_norm += tangent[i] * tangent[i];
        aux_vector[i] = interpolated_dzeta_ds_pseudo_normal[i];
        aux_norm += aux_vector[i] * aux_vector[i];
      }

      tangent[2] = interpolated_du_ds_tangent;
      tang_norm += tangent[2] * tangent[2];
      aux_vector[2] = interpolated_du_ds_pseudo_normal;
      aux_norm += aux_vector[2] * aux_vector[2];
    }

    // Norm of the rate of change
    norm_of_drds = sqrt(tang_norm);

    // Normalise
    double tang_norm_fact = 1.0 / sqrt(tang_norm);
    double aux_norm_fact = 1.0 / sqrt(aux_norm);
    for (unsigned i = 0; i < 3; i++)
    {
      tangent[i] *= tang_norm_fact;
      aux_vector[i] *= aux_norm_fact;
    }


    // Normal to meniscus is the cross product between the
    // two contact line vectors:
    Vector<double> meniscus_normal(3);
    ELEMENT::cross_product(tangent, aux_vector, meniscus_normal);

    // Calculate the norm
    double men_norm_fact = 0.0;
    for (unsigned i = 0; i < 3; i++)
    {
      men_norm_fact += meniscus_normal[i] * meniscus_normal[i];
    }

    // Normalise and adjust direction
    double sign = -double(normal_sign());
    for (unsigned i = 0; i < 3; i++)
    {
      meniscus_normal[i] *= sign / sqrt(men_norm_fact);
    }

    // The actual (bi) normal to the contact line is given
    // by another cross product
    ELEMENT::cross_product(meniscus_normal, tangent, normal);
  }


  //============================================================
  // Build the required elements
  //============================================================
  template class YoungLaplaceContactAngleElement<QYoungLaplaceElement<2>>;
  template class YoungLaplaceContactAngleElement<QYoungLaplaceElement<3>>;
  template class YoungLaplaceContactAngleElement<QYoungLaplaceElement<4>>;

  template class YoungLaplaceContactAngleElement<
    RefineableQYoungLaplaceElement<2>>;
  template class YoungLaplaceContactAngleElement<
    RefineableQYoungLaplaceElement<3>>;
  template class YoungLaplaceContactAngleElement<
    RefineableQYoungLaplaceElement<4>>;

} // namespace oomph
