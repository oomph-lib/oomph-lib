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
// Non-inline member functions for generic elements

#include "stored_shape_function_elements.h"
#include "shape.h"
#include "integral.h"


namespace oomph
{
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  //  Functions for finite elements
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  //=========================================================================
  /// Delete all the objects stored in the vectors:
  /// Shape_stored, DShape_local_stored, and D2Shape_local_stored
  //========================================================================
  void StorableShapeElementBase::delete_all_shape_local_stored()
  {
    delete_shape_local_stored();
    delete_dshape_local_stored();
    delete_d2shape_local_stored();
  }


  //=========================================================================
  /// Delete stored shape functions
  //========================================================================
  void StorableShapeElementBase::delete_shape_local_stored()
  {
    // Can this element delete the stored data?
    if ((Can_delete_shape_local_stored) && (Shape_stored_pt))
    {
      // Find the number of entries in the shape function storage vector
      // N.B. This *should* be the same for all three vectors
      unsigned n_intpt = Shape_stored_pt->size();
      // Loop over the entries of each vector, in reverse and delete them
      for (unsigned ipt = n_intpt; ipt > 0; ipt--)
      {
        delete (*Shape_stored_pt)[ipt - 1];
        (*Shape_stored_pt)[ipt - 1] = 0;
      }
      // Delete the vector itself
      delete Shape_stored_pt;
    }
    // Reset the pointer to zero {Must do this even if copied}
    Shape_stored_pt = 0;
  }


  //=========================================================================
  /// Delete stored derivatives of shape functions w.r.t. to local
  /// coordinates
  //========================================================================
  void StorableShapeElementBase::delete_dshape_local_stored()
  {
    // Can this element delete the stored data?
    if ((Can_delete_shape_local_stored) && (DShape_local_stored_pt))
    {
      unsigned n_intpt = DShape_local_stored_pt->size();
      for (unsigned ipt = n_intpt; ipt > 0; ipt--)
      {
        delete (*DShape_local_stored_pt)[ipt - 1];
        (*DShape_local_stored_pt)[ipt - 1] = 0;
      }
      // Delete the vector itself
      delete DShape_local_stored_pt;
    }
    // Reset the pointer to zero, must do this, even if copied
    DShape_local_stored_pt = 0;
  }


  //=========================================================================
  /// Delete stored second derivatives of shape functions w.r.t. to local
  /// coordinates
  //========================================================================
  void StorableShapeElementBase::delete_d2shape_local_stored()
  {
    // Can this element delete the stored data?
    if ((Can_delete_shape_local_stored) && (D2Shape_local_stored_pt))
    {
      unsigned n_intpt = D2Shape_local_stored_pt->size();
      for (unsigned ipt = n_intpt; ipt > 0; ipt--)
      {
        delete (*D2Shape_local_stored_pt)[ipt - 1];
        (*D2Shape_local_stored_pt)[ipt - 1] = 0;
      }
      // Delete the vector itself
      delete D2Shape_local_stored_pt;
    }
    // Reset the pointer to zero, must do it, even if copied
    D2Shape_local_stored_pt = 0;
  }


  //=========================================================================
  /// Delete all stored quantities related to derivatives of shape
  /// fcts w.r.t. to global Eulerian coordinates
  //========================================================================
  void StorableShapeElementBase::delete_all_dshape_eulerian_stored()
  {
    delete_dshape_eulerian_stored();
    delete_d2shape_eulerian_stored();
    delete_J_eulerian_stored();
  }

  //=========================================================================
  /// Delete stored derivatives w.r.t. Eulerian coordinates
  //========================================================================
  void StorableShapeElementBase::delete_dshape_eulerian_stored()
  {
    // If the storage has been allocated and we can delete it
    if ((Can_delete_dshape_eulerian_stored) && (DShape_eulerian_stored_pt))
    {
      // Find the number of entries in the first vector
      unsigned n_intpt = DShape_eulerian_stored_pt->size();
      // Loop over the entries of the vectors, in reverse and delete them
      for (unsigned ipt = n_intpt; ipt > 0; ipt--)
      {
        delete (*DShape_eulerian_stored_pt)[ipt - 1];
        (*DShape_eulerian_stored_pt)[ipt - 1] = 0;
      }
      // Delete the vector itself
      delete DShape_eulerian_stored_pt;
    }
    // Reset the pointer to zero, even if copied
    DShape_eulerian_stored_pt = 0;
  }


  //=========================================================================
  /// Delete stored 2nd derivatives w.r.t. Eulerian coordinates
  //========================================================================
  void StorableShapeElementBase::delete_d2shape_eulerian_stored()
  {
    // If the storage has been allocated
    if ((Can_delete_dshape_eulerian_stored) && (D2Shape_eulerian_stored_pt))
    {
      // Find the number of entries in the second vector
      unsigned n_intpt = D2Shape_eulerian_stored_pt->size();
      // Loop over the entries in reverse and delete them
      for (unsigned ipt = n_intpt; ipt > 0; ipt--)
      {
        delete (*D2Shape_eulerian_stored_pt)[ipt - 1];
        (*D2Shape_eulerian_stored_pt)[ipt - 1] = 0;
      }
      // Delete the vector itself
      delete D2Shape_eulerian_stored_pt;
    }
    // Reset the pointer to zero, even if copied
    D2Shape_eulerian_stored_pt = 0;
  }


  //=========================================================================
  /// Delete stored Jacobian of mapping between local and global Eulerian
  /// coordinates
  //========================================================================
  void StorableShapeElementBase::delete_J_eulerian_stored()
  {
    // If the element originally allocated the storage, delete it
    if (Can_delete_dshape_eulerian_stored)
    {
      // Delete the stored Jacobians
      delete Jacobian_eulerian_stored_pt;
    }
    // Reset the pointer to zero, even if copied
    Jacobian_eulerian_stored_pt = 0;
  }


  //======================================================================
  /// \short The destructor cleans up the memory allocated
  /// for shape function storage.
  //=======================================================================
  StorableShapeElementBase::~StorableShapeElementBase()
  {
    // Merely need to call the private functions to clean up storages
    delete_all_shape_local_stored();
    delete_all_dshape_eulerian_stored();
  }

  //=======================================================================
  /// Set the spatial integration scheme and also calculate the values of the
  /// shape functions and their derivatives w.r.t. the local coordinates,
  /// placing the values into storage so that they may be re-used,
  /// without recalculation
  //=======================================================================
  void StorableShapeElementBase::set_integration_scheme(
    Integral* const& integral_pt)
  {
    // Assign the integration scheme
    FiniteElement::set_integration_scheme(integral_pt);

    // If we are storing the shape functions and first and second derivatives
    if (D2Shape_local_stored_pt != 0)
    {
      pre_compute_d2shape_local_at_knots();
    }
    // If we are storing the shape functions and first derivatives
    else if (DShape_local_stored_pt != 0)
    {
      pre_compute_dshape_local_at_knots();
    }
    // If we are storing the shape functions
    else if (Shape_stored_pt != 0)
    {
      pre_compute_shape_at_knots();
    }

    // If we are storing Eulerian first and second derivatives, recompute them
    if (D2Shape_eulerian_stored_pt != 0)
    {
      pre_compute_d2shape_eulerian_at_knots();
    }
    // If we are storing Eulerian first derivatives, recompute them
    else if (DShape_eulerian_stored_pt != 0)
    {
      pre_compute_dshape_eulerian_at_knots();
    }
    // If we are storing the Jacobian of the mapping from local to Eulerian
    // coordinates, recompute it
    else if (Jacobian_eulerian_stored_pt != 0)
    {
      pre_compute_J_eulerian_at_knots();
    }
  }

  //========================================================================
  /// Calculate the shape functions at the integration points and store in
  /// internal storage of the element
  //========================================================================
  void StorableShapeElementBase::pre_compute_shape_at_knots()
  {
    // Find the number of nodes in the element
    unsigned n_node = nnode();
#ifdef PARANOID
    if (n_node == 0)
    {
      std::string error_message =
        "FiniteElement::Node_pt must be resized to a value greater than\n";
      error_message += "zero before calling pre_compute_shape_at_knots()";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif
    // Find number of interpolated position dofs
    unsigned n_position_type = nnodal_position_type();

    // In case we have an exisiting integration scheme, wipe the stored data
    delete_shape_local_stored();
    // Element is now in charge of deleting its own stored data again
    Can_delete_shape_local_stored = true;

    // Allocate internal storage for the shape functions
    Shape_stored_pt = new Vector<Shape*>;

    // Storage for the shape functions and their local derivatives
    Shape psi(n_node, n_position_type);

    // Loop over the integration points
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the shape functions
      FiniteElement::shape_at_knot(ipt, psi);

      // Set up local storage for the shape functions and derivatives
      Shape* psi_pt = new Shape(n_node, n_position_type);

      // Copy the values of the shape functions and their local derivatives
      // into the element storage
      for (unsigned n = 0; n < n_node; n++)
      {
        for (unsigned k = 0; k < n_position_type; k++)
        {
          (*psi_pt)(n, k) = psi(n, k);
        }
      }

      // Add the pointers to the shape functions to the internal storage
      Shape_stored_pt->push_back(psi_pt);
    } // End of loop over integration points
  }

  //========================================================================
  /// Calculate the shape functions and thir first derivatives
  /// w.r.t. local coordinates at the integration points and store in
  /// internal storage of the element
  //========================================================================
  void StorableShapeElementBase::pre_compute_dshape_local_at_knots()
  {
    // Find the number of nodes in the element
    unsigned n_node = nnode();
#ifdef PARANOID
    if (n_node == 0)
    {
      std::string error_message =
        "FiniteElement::Node_pt must be resized to a value greater than\n";
      error_message +=
        "zero before calling pre_compute_dshape_local_at_knots()";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif
    // Find number of interpolated position dofs
    unsigned n_position_type = nnodal_position_type();
    // Find spatial dimension of the element
    unsigned Dim = dim();

    // In case we have an exisiting integration scheme, wipe the stored data
    delete_shape_local_stored();
    delete_dshape_local_stored();
    // Element is now in charge of deleting its own stored data again
    Can_delete_shape_local_stored = true;

    // Allocate internal storage for the shape functions
    Shape_stored_pt = new Vector<Shape*>;
    DShape_local_stored_pt = new Vector<DShape*>;

    // Storage for the shape functions and their local derivatives
    Shape psi(n_node, n_position_type);
    DShape dpsids(n_node, n_position_type, Dim);

    // Loop over the integration points
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the shape functions and local derivatives at the integration point
      FiniteElement::dshape_local_at_knot(ipt, psi, dpsids);

      // Set up local storage for the shape functions and derivatives
      Shape* psi_pt = new Shape(n_node, n_position_type);
      DShape* dpsids_pt = new DShape(n_node, n_position_type, Dim);

      // Copy the values of the shape functions and their local derivatives
      // into the element storage
      for (unsigned n = 0; n < n_node; n++)
      {
        for (unsigned k = 0; k < n_position_type; k++)
        {
          (*psi_pt)(n, k) = psi(n, k);

          for (unsigned i = 0; i < Dim; i++)
          {
            (*dpsids_pt)(n, k, i) = dpsids(n, k, i);
          }
        }
      }

      // Add the pointers to the shape functions and derivatives to the internal
      // storage
      Shape_stored_pt->push_back(psi_pt);
      DShape_local_stored_pt->push_back(dpsids_pt);
    } // End of loop over integration points
  }

  //========================================================================
  /// Calculate the shape functions and thir first and second derivatives
  /// w.r.t. local coordinates at the integration points and store in
  /// internal storage of the element
  //========================================================================
  void StorableShapeElementBase::pre_compute_d2shape_local_at_knots()
  {
    // Find the number of nodes in the element
    unsigned n_node = nnode();
#ifdef PARANOID
    if (n_node == 0)
    {
      std::string error_message =
        "FiniteElement::Node_pt must be resized to a value greater than\n";
      error_message +=
        "zero before calling pre_compute_d2shape_local_at_knots()";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif
    // Find number of interpolated position dofs
    unsigned n_position_type = nnodal_position_type();
    // Find spatial dimension of the element
    unsigned Dim = dim();

    // Find the number of second derivatives required
    // N.B. We are assuming that the mixed derivatives are symmetric here
    unsigned n_deriv = 0;
    switch (Dim)
    {
      case 1:
        n_deriv = 1;
        break;
      case 2:
        n_deriv = 3;
        break;
      case 3:
        n_deriv = 6;
        break;
      default:
        oomph_info << "Really more than 3 dimensions?" << std::endl;
        break;
    }

    // In case we have an exisiting integration scheme, wipe the stored data
    delete_shape_local_stored();
    delete_dshape_local_stored();
    delete_d2shape_local_stored();
    // Element is now in charge of deleting its own stored data again
    Can_delete_shape_local_stored = true;

    // Allocate internal storage for the shape functions
    Shape_stored_pt = new Vector<Shape*>;
    DShape_local_stored_pt = new Vector<DShape*>;
    D2Shape_local_stored_pt = new Vector<DShape*>;

    // Storage for the shape functions and their local derivatives
    Shape psi(n_node, n_position_type);
    DShape dpsids(n_node, n_position_type, Dim);
    DShape d2psids(n_node, n_position_type, n_deriv);

    // Loop over the integration points
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the shape functions and local derivatives at the integration point
      FiniteElement::d2shape_local_at_knot(ipt, psi, dpsids, d2psids);

      // Set up local storage for the shape functions and derivatives
      Shape* psi_pt = new Shape(n_node, n_position_type);
      DShape* dpsids_pt = new DShape(n_node, n_position_type, Dim);
      DShape* d2psids_pt = new DShape(n_node, n_position_type, n_deriv);

      // Copy the values of the shape functions and their local derivatives
      // into the element storage
      for (unsigned n = 0; n < n_node; n++)
      {
        for (unsigned k = 0; k < n_position_type; k++)
        {
          (*psi_pt)(n, k) = psi(n, k);

          for (unsigned i = 0; i < Dim; i++)
          {
            (*dpsids_pt)(n, k, i) = dpsids(n, k, i);
          }

          for (unsigned i = 0; i < n_deriv; i++)
          {
            (*d2psids_pt)(n, k, i) = d2psids(n, k, i);
          }
        }
      }

      // Add the pointers to the shape functions and derivatives to the internal
      // storage
      Shape_stored_pt->push_back(psi_pt);
      DShape_local_stored_pt->push_back(dpsids_pt);
      D2Shape_local_stored_pt->push_back(d2psids_pt);
    } // End of loop over integration points
  }

  //=======================================================================
  /// Calculate the value of the Jacobian of the mapping from local to
  /// global coordinates at the integration points and store internally
  //=======================================================================
  void StorableShapeElementBase::pre_compute_J_eulerian_at_knots()
  {
    // Delete previously existing storage
    delete_J_eulerian_stored();
    // Now we're in change of deletion again
    Can_delete_dshape_eulerian_stored = true;

    // Allocate storage for the stored Jacobian values
    Jacobian_eulerian_stored_pt = new Vector<double>;

    // Loop over the integration points
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Add the value of the Jacobian to the internally stored vector
      Jacobian_eulerian_stored_pt->push_back(
        FiniteElement::J_eulerian_at_knot(ipt));
    }
  }

  //========================================================================
  /// Calculate the values of the derivatives of the shape functions at the
  /// integration points and store in the internal storage of the element
  //========================================================================
  void StorableShapeElementBase::pre_compute_dshape_eulerian_at_knots()
  {
    // Pre-compute the basic shape functions
    pre_compute_shape_at_knots();

    // Find the number of nodes
    unsigned n_node = nnode();
    // Get the number of position types and the dimension from the element
    unsigned n_position_type = this->nnodal_position_type();
    unsigned n_dim = this->nodal_dimension();

    // Delete the exisiting stored objects
    delete_J_eulerian_stored();
    delete_dshape_eulerian_stored();
    // Now we're in change of deletion again
    Can_delete_dshape_eulerian_stored = true;

    // Allocate storage for the stored shape function derivatives
    DShape_eulerian_stored_pt = new Vector<DShape*>;
    Jacobian_eulerian_stored_pt = new Vector<double>;

    // Assign local variables for the shape function and derivatives
    Shape psi(n_node, n_position_type);
    DShape dpsidx(n_node, n_position_type, n_dim);

    // Loop over the integration points
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the values of the shape function and derivatives at the
      // integration point and add to the value of the Jacobian to the
      // internally stored vector
      Jacobian_eulerian_stored_pt->push_back(
        FiniteElement::dshape_eulerian_at_knot(ipt, psi, dpsidx));

      // Set up local storage for the shape function derivatives
      DShape* dpsidx_pt = new DShape(n_node, n_position_type, n_dim);

      // Now copy the values over
      for (unsigned l = 0; l < n_node; l++)
      {
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // First derivatives
          for (unsigned i = 0; i < n_dim; i++)
          {
            (*dpsidx_pt)(l, k, i) = dpsidx(l, k, i);
          }
        }
      }

      // Add the pointer to the vector of stored DShape objects
      DShape_eulerian_stored_pt->push_back(dpsidx_pt);
    } // End of loop over integration points
  }

  //========================================================================
  /// Calculate the values of the first and second derivatives of the shape
  /// functions at the integration points and store in the internal storage
  /// of the element
  //========================================================================
  void StorableShapeElementBase::pre_compute_d2shape_eulerian_at_knots()
  {
    // Pre-compute the basic shape functions
    pre_compute_shape_at_knots();

    // Find the number of nodes
    unsigned n_node = nnode();
    // Get the number of position types and the dimension from element
    unsigned n_position_type = this->nnodal_position_type();
    unsigned n_dim = this->nodal_dimension();

    // Find the number of second derivatives required
    // N.B. We are assuming that the mixed derivatives are symmetric here
    unsigned n_deriv = 0;
    switch (n_dim)
    {
      case 1:
        n_deriv = 1;
        break;
      case 2:
        n_deriv = 3;
        break;
      case 3:
        n_deriv = 6;
        break;
      default:
        oomph_info << "Really more than 3 dimensions?" << std::endl;
        break;
    }

    // Delete the existing objects, if there are any
    delete_J_eulerian_stored();
    delete_dshape_eulerian_stored();
    delete_d2shape_eulerian_stored();
    // Now we're in change of deletion again
    Can_delete_dshape_eulerian_stored = true;

    // Allocate storage for the stored shape function derivatives
    DShape_eulerian_stored_pt = new Vector<DShape*>;
    D2Shape_eulerian_stored_pt = new Vector<DShape*>;
    Jacobian_eulerian_stored_pt = new Vector<double>;

    // Assign local variables for the shape function and derivatives
    Shape psi(n_node, n_position_type);
    DShape dpsidx(n_node, n_position_type, n_dim);
    DShape d2psidx(n_node, n_position_type, n_deriv);

    // Loop over the integration points
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the values of the shape function and derivatives at the
      // integration point and assign the value of the Jacobian to the
      // internally stored vector
      Jacobian_eulerian_stored_pt->push_back(
        FiniteElement::d2shape_eulerian_at_knot(ipt, psi, dpsidx, d2psidx));

      // Set up local storage for the shape function derivatives
      DShape* dpsidx_pt = new DShape(n_node, n_position_type, n_dim);
      DShape* d2psidx_pt = new DShape(n_node, n_position_type, n_deriv);

      // Now copy the values over
      for (unsigned l = 0; l < n_node; l++)
      {
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // First derivatives
          for (unsigned i = 0; i < n_dim; i++)
          {
            (*dpsidx_pt)(l, k, i) = dpsidx(l, k, i);
          }

          // Second derivatives
          for (unsigned i = 0; i < n_deriv; i++)
          {
            (*d2psidx_pt)(l, k, i) = d2psidx(l, k, i);
          }
        }
      }

      // Add the pointers to the shape function derivatives to the internal
      // storage
      DShape_eulerian_stored_pt->push_back(dpsidx_pt);
      D2Shape_eulerian_stored_pt->push_back(d2psidx_pt);
    } // End of loop over the shape functions
  }

  //=========================================================================
  /// \short Return the shape function stored at the ipt-th integration
  /// point.
  //=========================================================================
  void StorableShapeElementBase::shape_at_knot(const unsigned& ipt,
                                               Shape& psi) const
  {
    // If we are not storing the shape functions, calculate the values
    if (Shape_stored_pt == 0)
    {
      FiniteElement::shape_at_knot(ipt, psi);
    }
    else
    {
      // Read out the stored shape functions
      // Note this will copy the values by pointer (fast)
      psi = (*Shape_stored_pt)[ipt];
    }
  }


  //=========================================================================
  /// \short Return the shape function and its derivatives w.r.t. the local
  /// coordinates at the ipt-th integration point.
  //=========================================================================
  void StorableShapeElementBase::dshape_local_at_knot(const unsigned& ipt,
                                                      Shape& psi,
                                                      DShape& dpsids) const
  {
    // If we are not storing the first derivatives, calculate them
    if (DShape_local_stored_pt == 0)
    {
      FiniteElement::dshape_local_at_knot(ipt, psi, dpsids);
    }
    else
    {
      // Read out the stored shape functions
      // Set the internal pointers in psi and dpsids
      psi = (*Shape_stored_pt)[ipt];
      dpsids = (*DShape_local_stored_pt)[ipt];
    }
  }

  //=========================================================================
  /// \short Return the shape function and its first and second derivatives
  /// w.r.t. the local coordinates at the ipt-th integration point.
  //=========================================================================
  void StorableShapeElementBase::d2shape_local_at_knot(const unsigned& ipt,
                                                       Shape& psi,
                                                       DShape& dpsids,
                                                       DShape& d2psids) const
  {
    // If we are not storing the second derivatives, calculate them on the fly
    if (D2Shape_local_stored_pt == 0)
    {
      FiniteElement::d2shape_local_at_knot(ipt, psi, dpsids, d2psids);
    }
    else
    {
      // Read out the stored shape functions
      // Set the internal pointers in psi, dpsids, and d2psids
      psi = (*Shape_stored_pt)[ipt];
      dpsids = (*DShape_local_stored_pt)[ipt];
      d2psids = (*D2Shape_local_stored_pt)[ipt];
    }
  }


  //==========================================================================
  /// \short Compute the geometric shape functions, and
  /// derivatives w.r.t eulerian coordinates at the ipt-th integration point.
  /// If the values have already been computed, return the stored values.
  //==========================================================================
  double StorableShapeElementBase::dshape_eulerian_at_knot(const unsigned& ipt,
                                                           Shape& psi,
                                                           DShape& dpsidx) const
  {
    // If we are not storing the values, return the calculated values
    if (DShape_eulerian_stored_pt == 0)
    {
      return FiniteElement::dshape_eulerian_at_knot(ipt, psi, dpsidx);
    }
    else
    {
      // Set internal pointers in the shape functions.
      psi = (*Shape_stored_pt)[ipt];
      dpsidx = (*DShape_eulerian_stored_pt)[ipt];

      // Return the stored value of the jacobian
      return ((*Jacobian_eulerian_stored_pt)[ipt]);
    }
  }

  //==========================================================================
  /// \short Return the geometric shape functions, first and second
  /// derivatives w.r.t eulerian coordinates at the ipt-th integration point.
  /// If the values have already been computed, return the stored values.
  //==========================================================================
  double StorableShapeElementBase::d2shape_eulerian_at_knot(
    const unsigned& ipt, Shape& psi, DShape& dpsidx, DShape& d2psidx) const
  {
    // If we are not storing the values return the calculated values
    if (D2Shape_eulerian_stored_pt == 0)
    {
      return FiniteElement::d2shape_eulerian_at_knot(ipt, psi, dpsidx, d2psidx);
    }
    else
    {
      // Set internal pointers in the shape functions
      psi = (*Shape_stored_pt)[ipt];
      dpsidx = (*DShape_eulerian_stored_pt)[ipt];
      d2psidx = (*D2Shape_eulerian_stored_pt)[ipt];

      // Return the stored value of the jacobian
      return ((*Jacobian_eulerian_stored_pt)[ipt]);
    }
  }

  //==========================================================================
  /// \short Return the Jacobian of the mapping between the local and global
  /// coordinates. If the value has been precomputed return that
  //==========================================================================
  double StorableShapeElementBase::J_eulerian_at_knot(const unsigned& ipt) const
  {
    // If we are not storing the values, return the calculated values
    if (Jacobian_eulerian_stored_pt == 0)
    {
      return FiniteElement::J_eulerian_at_knot(ipt);
    }
    else
    {
      // Return the stored value of the jacobian
      return ((*Jacobian_eulerian_stored_pt)[ipt]);
    }
  }

  //========================================================================
  /// Set the shape functions referenced in the internal vectors to be
  /// those stored in the StorableShapeElementBase pointed to by element_pt.
  /// Using this function will allow a saving in the storage required
  /// for integration schemes in the (most common)
  /// case when a large number of elements have the same integration scheme
  //=======================================================================
  void StorableShapeElementBase::set_shape_local_stored_from_element(
    StorableShapeElementBase* const& element_pt)
  {
#ifdef PARANOID
    // Check that we aren't nulling out
    if (element_pt->shape_stored_pt() == 0)
    {
      std::string error_message =
        "Element does not have stored shape functions\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Only do this if the referenced Shape objects are not already the same
    // Assume that if the stored shape functions are the same, the rest will be
    if (Shape_stored_pt != element_pt->shape_stored_pt())
    {
      // Delete the existing data
      delete_all_shape_local_stored();
      // Now this element can no longer delete the data pointed at by
      // the internal vectors
      Can_delete_shape_local_stored = false;

      // Assign the pointers
      Shape_stored_pt = element_pt->shape_stored_pt();
      DShape_local_stored_pt = element_pt->dshape_local_stored_pt();
      D2Shape_local_stored_pt = element_pt->d2shape_local_stored_pt();
    }
  }

  //========================================================================
  /// Set the stored derivatives of shape functions w.r.t Eulerian coordinates
  /// to be
  /// those stored in the StorableShapeElementBase pointed to by element_pt.
  /// Using this function will allow a saving in the storage required
  /// for integration schemes in the (most common)
  /// case when a large number of elements have the same integration scheme
  //=======================================================================
  void StorableShapeElementBase::set_dshape_eulerian_stored_from_element(
    StorableShapeElementBase* const& element_pt)
  {
    set_shape_local_stored_from_element(element_pt);

    // Only do this if the referenced Shape objects are not already the same
    // Assume that if the stored shape functions are the same, the rest will be
    if (DShape_eulerian_stored_pt != element_pt->dshape_eulerian_stored_pt())
    {
      // Delete the existing data
      delete_all_dshape_eulerian_stored();
      // Now this element can no longer delete the data pointed at by
      // the internal vectors
      Can_delete_dshape_eulerian_stored = false;

      // Assign the pointers
      DShape_eulerian_stored_pt = element_pt->dshape_eulerian_stored_pt();
      D2Shape_eulerian_stored_pt = element_pt->d2shape_eulerian_stored_pt();
      Jacobian_eulerian_stored_pt = element_pt->jacobian_eulerian_stored_pt();
    }
  }


  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  //  Functions for solid elements with stored shape functions
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  //=========================================================================
  /// Delete all the objects and storage associated with stored derivatives
  /// of shape functions with respect to Lagrangian coordinates
  //=========================================================================
  void StorableShapeSolidElementBase::delete_all_dshape_lagrangian_stored()
  {
    delete_dshape_lagrangian_stored();
    delete_d2shape_lagrangian_stored();
    delete_J_lagrangian_stored();
  }


  //=========================================================================
  /// Delete all the objects stored in the vectors:
  /// DShape_lagrangian_stored
  //========================================================================
  void StorableShapeSolidElementBase::delete_dshape_lagrangian_stored()
  {
    // If storage has been allocated
    if ((Can_delete_dshape_lagrangian_stored) && (DShape_lagrangian_stored_pt))
    {
      // Find the number of entries in the shape function storage vector
      unsigned n_intpt = DShape_lagrangian_stored_pt->size();
      // Loop over the entries of the vectors, in reverse and delete them
      for (unsigned ipt = n_intpt; ipt > 0; ipt--)
      {
        delete (*DShape_lagrangian_stored_pt)[ipt - 1];
        (*DShape_lagrangian_stored_pt)[ipt - 1] = 0;
      }
      // Delete the actual vector
      delete DShape_lagrangian_stored_pt;
    }
    // Reset the pointer to zero
    DShape_lagrangian_stored_pt = 0;
  }


  //=========================================================================
  /// Delete all the objects stored in the vectors:
  /// D2Shape_lagrangian_stored
  //========================================================================
  void StorableShapeSolidElementBase::delete_d2shape_lagrangian_stored()
  {
    // If storage has been allocated
    if ((Can_delete_dshape_lagrangian_stored) && (D2Shape_lagrangian_stored_pt))
    {
      // Now find the number of entries in the second vector
      unsigned n_intpt = D2Shape_lagrangian_stored_pt->size();
      // Loop over the entries in reverse and delete them
      for (unsigned ipt = n_intpt; ipt > 0; ipt--)
      {
        delete (*D2Shape_lagrangian_stored_pt)[ipt - 1];
        (*D2Shape_lagrangian_stored_pt)[ipt - 1] = 0;
      }
      // Delete the actual vector
      delete D2Shape_lagrangian_stored_pt;
    }
    // Reset the pointer to zero
    D2Shape_lagrangian_stored_pt = 0;
  }

  //=======================================================================
  /// Delete the stored Jacobian of the mapping between the Lagrangian and
  /// local coordinates.
  //=======================================================================
  void StorableShapeSolidElementBase::delete_J_lagrangian_stored()
  {
    // If we allocated the storage, delete it
    if (Can_delete_dshape_lagrangian_stored)
    {
      // Delete the stored Jacobian
      delete Jacobian_lagrangian_stored_pt;
    }
    // Reset the pointer to zero
    Jacobian_lagrangian_stored_pt = 0;
  }

  //========================================================================
  /// Calculate the values of the derivatives of the shape functions at the
  /// integration points and store in the internal storage of the element
  //========================================================================
  void StorableShapeSolidElementBase::pre_compute_dshape_lagrangian_at_knots()
  {
    // Pre-compute the basic shape functions
    pre_compute_shape_at_knots();

    // Find the number of nodes
    unsigned n_node = nnode();
    // Get the number of position types and the dimension from first node
    // N.B. Assume that it is the same for all nodes
    unsigned n_lagrangian_type =
      static_cast<SolidNode*>(node_pt(0))->nlagrangian_type();
    unsigned n_lagrangian = static_cast<SolidNode*>(node_pt(0))->nlagrangian();

    // Delete the exisiting stored objects
    delete_J_lagrangian_stored();
    delete_dshape_lagrangian_stored();
    // Now we're in charge of deletion again
    Can_delete_dshape_lagrangian_stored = true;

    // Allocate storage for the stored shape function derivatives
    DShape_lagrangian_stored_pt = new Vector<DShape*>;
    Jacobian_lagrangian_stored_pt = new Vector<double>;

    // Assign local variables for the shape function and derivatives
    Shape psi(n_node, n_lagrangian_type);
    DShape dpsidxi(n_node, n_lagrangian_type, n_lagrangian);

    // Loop over the integration points
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the values of the shape function and derivatives at the
      // integration point and add to the value of the Jacobian to the
      // internally stored vector
      Jacobian_lagrangian_stored_pt->push_back(
        SolidFiniteElement::dshape_lagrangian_at_knot(ipt, psi, dpsidxi));

      // Set up local storage for the shape function derivatives
      DShape* dpsidxi_pt = new DShape(n_node, n_lagrangian_type, n_lagrangian);

      // Now copy the values over
      for (unsigned l = 0; l < n_node; l++)
      {
        for (unsigned k = 0; k < n_lagrangian_type; k++)
        {
          // First derivatives
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            (*dpsidxi_pt)(l, k, i) = dpsidxi(l, k, i);
          }
        }
      }

      // Add the pointer to the vector of stored DShape objects
      DShape_lagrangian_stored_pt->push_back(dpsidxi_pt);
    } // End of loop over integration points
  }

  //========================================================================
  /// Calculate the values of the first and second derivatives of the shape
  /// functions at the integration points and store in the internal storage
  /// of the element
  //========================================================================
  void StorableShapeSolidElementBase::pre_compute_d2shape_lagrangian_at_knots()
  {
    // Pre-compute the basic shape functions
    pre_compute_shape_at_knots();

    // Find the number of nodes
    unsigned n_node = nnode();
    // Get the number of position types and the dimension from first node
    // N.B. Assume that it is the same for all nodes
    unsigned n_lagrangian_type =
      static_cast<SolidNode*>(node_pt(0))->nlagrangian_type();
    unsigned n_lagrangian = static_cast<SolidNode*>(node_pt(0))->nlagrangian();

    // Find the number of second derivatives required
    // N.B. We are assuming that the mixed derivatives are symmetric here
    unsigned n_deriv = 0;
    switch (n_lagrangian)
    {
      case 1:
        n_deriv = 1;
        break;
      case 2:
        n_deriv = 3;
        break;
      case 3:
        n_deriv = 6;
        break;
      default:
        oomph_info << "Really more than 3 dimensions?" << std::endl;
        break;
    }

    // Delete the existing objects, if there are any
    delete_J_lagrangian_stored();
    delete_dshape_lagrangian_stored();
    delete_d2shape_lagrangian_stored();
    // We are in charge of deleting again
    Can_delete_dshape_lagrangian_stored = true;

    // Allocate storage for the stored shape function derivatives
    DShape_lagrangian_stored_pt = new Vector<DShape*>;
    D2Shape_lagrangian_stored_pt = new Vector<DShape*>;
    Jacobian_lagrangian_stored_pt = new Vector<double>;

    // Assign local variables for the shape function and derivatives
    Shape psi(n_node, n_lagrangian_type);
    DShape dpsidxi(n_node, n_lagrangian_type, n_lagrangian);
    DShape d2psidxi(n_node, n_lagrangian_type, n_deriv);

    // Loop over the integration points
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the values of the shape function and derivatives at the
      // integration point and assign the value of the Jacobian to the
      // internally stored vector
      Jacobian_lagrangian_stored_pt->push_back(
        SolidFiniteElement::d2shape_lagrangian_at_knot(
          ipt, psi, dpsidxi, d2psidxi));

      // Set up local storage for the shape function derivatives
      DShape* dpsidxi_pt = new DShape(n_node, n_lagrangian_type, n_lagrangian);
      DShape* d2psidxi_pt = new DShape(n_node, n_lagrangian_type, n_deriv);

      // Now copy the values over
      for (unsigned l = 0; l < n_node; l++)
      {
        for (unsigned k = 0; k < n_lagrangian_type; k++)
        {
          // First derivatives
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            (*dpsidxi_pt)(l, k, i) = dpsidxi(l, k, i);
          }

          // Second derivatives
          for (unsigned i = 0; i < n_deriv; i++)
          {
            (*d2psidxi_pt)(l, k, i) = d2psidxi(l, k, i);
          }
        }
      }

      // Add the pointers to the shape function derivatives to the internal
      // storage
      DShape_lagrangian_stored_pt->push_back(dpsidxi_pt);
      D2Shape_lagrangian_stored_pt->push_back(d2psidxi_pt);
    } // End of loop over the shape functions
  }


  //==========================================================================
  /// \short Compute the geometric shape functions, and
  /// derivatives w.r.t Lagrangian coordinates at the ipt-th integration point.
  /// If the values have already been computed, return the stored values.
  //==========================================================================
  double StorableShapeSolidElementBase::dshape_lagrangian_at_knot(
    const unsigned& ipt, Shape& psi, DShape& dpsidxi) const
  {
    // If we are not storing the values, return the calculated values
    if (DShape_lagrangian_stored_pt == 0)
    {
      return SolidFiniteElement::dshape_lagrangian_at_knot(ipt, psi, dpsidxi);
    }
    else
    {
      // Set the internal pointers in the shape functions
      psi = shape_stored_pt(ipt);
      dpsidxi = (*DShape_lagrangian_stored_pt)[ipt];

      // Return the stored value of the jacobian
      return ((*Jacobian_lagrangian_stored_pt)[ipt]);
    }
  }

  //==========================================================================
  /// \short Compute the geometric shape functions, first and second
  /// derivatives w.r.t Lagrangian coordinates at the ipt-th integration point.
  /// If the values have already been computed, return the stored values.
  //==========================================================================
  double StorableShapeSolidElementBase::d2shape_lagrangian_at_knot(
    const unsigned& ipt, Shape& psi, DShape& dpsidxi, DShape& d2psidxi) const
  {
    // If we are not storing the values return the calculated values
    if (D2Shape_lagrangian_stored_pt == 0)
    {
      return SolidFiniteElement::d2shape_lagrangian_at_knot(
        ipt, psi, dpsidxi, d2psidxi);
    }
    else
    {
      // Set the internal values of the pointers in the Shape objects
      psi = shape_stored_pt(ipt);
      dpsidxi = (*DShape_lagrangian_stored_pt)[ipt];
      d2psidxi = (*D2Shape_lagrangian_stored_pt)[ipt];

      // Return the stored value of the jacobian
      return ((*Jacobian_lagrangian_stored_pt)[ipt]);
    }
  }


  //========================================================================
  /// Set the stored derivatives of shape functions w.r.t Lagrangian coordinates
  /// to be
  /// those stored in the StorableShapeElementBase pointed to by element_pt.
  /// Using this function will allow a saving in the storage required
  /// for integration schemes in the (most common)
  /// case when a large number of elements have the same integration scheme
  //=======================================================================
  void StorableShapeSolidElementBase::set_dshape_lagrangian_stored_from_element(
    StorableShapeSolidElementBase* const& element_pt)
  {
    set_shape_local_stored_from_element(element_pt);

    // Only do this if the referenced Shape objects are not already the same
    // Assume that if the stored shape functions are the same, the rest will be
    if (DShape_lagrangian_stored_pt !=
        element_pt->dshape_lagrangian_stored_pt())
    {
      // Delete the existing data
      delete_all_dshape_lagrangian_stored();
      // Now this element can no longer delete the data pointed at by
      // the internal vectors
      Can_delete_dshape_lagrangian_stored = false;

      // Assign the pointers
      DShape_lagrangian_stored_pt = element_pt->dshape_lagrangian_stored_pt();
      D2Shape_lagrangian_stored_pt = element_pt->d2shape_lagrangian_stored_pt();
      Jacobian_lagrangian_stored_pt =
        element_pt->jacobian_lagrangian_stored_pt();
    }
  }


} // namespace oomph
