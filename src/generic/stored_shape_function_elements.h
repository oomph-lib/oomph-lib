// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
// Header file for classes that overload elements and implement
// stored shape functions.

// Include guard to prevent multiple inclusions of the header
#ifndef OOMPH_STORED_SHAPE_FUNCTION_ELEMENTS_HEADER
#define OOMPH_STORED_SHAPE_FUNCTION_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include "elements.h"


namespace oomph
{
// Debugging flag
#define OOMPH_STORED_SHAPE_FUNCTIONS_VERBOSE
#undef OOMPH_STORED_SHAPE_FUNCTIONS_VERBOSE


  //==========================================================================
  /// Base class for elements that allow storage of precomputed shape functions
  /// and their derivatives w.r.t to the local and global (Eulerian)
  /// coordinates at the element's integration points.
  //==========================================================================
  class StorableShapeElementBase : public virtual FiniteElement
  {
  private:
    /// Pointer to storage for the pointers to the nodal shape functions
    /// at the integration points (knots)
    // ALH: Note that the vector must be Shape* because we do not know,
    // a priori, how much storage to allocate in each Shape object.
    // N.B. This could in 99% of cases be static, but that would not
    // permit individual elements to have different integration schemes.
    // If this proves to be a problem, one can use the the function
    // set_local_shape_stored_from_element(), which sets this pointer to
    // point to the pointer allocated by another element.
    Vector<Shape*>* Shape_stored_pt;

    /// Pointer to storage for the pointers to the derivatives of the
    /// nodal shape functions w.r.t. the local coordinates at integration points
    Vector<DShape*>* DShape_local_stored_pt;

    /// Pointer to storage for the pointers to the second derivatives of
    /// the nodal shape functions w.r.t. the local coordinates at integration
    /// points
    Vector<DShape*>* D2Shape_local_stored_pt;

    /// Boolean to determine whether the element can delete the stored
    /// local shape functions
    bool Can_delete_shape_local_stored;

    /// Pointer to storage for the derivatives of the
    /// shape functions w.r.t. global coordinates at integration points
    Vector<DShape*>* DShape_eulerian_stored_pt;

    /// Pointer to storage for the second derivatives of the
    /// shape functions w.r.t. global coordinates at integration points
    Vector<DShape*>* D2Shape_eulerian_stored_pt;

    /// Pointer to storage for the Jacobian of the element w.r.t
    /// global coordinates
    Vector<double>* Jacobian_eulerian_stored_pt;

    /// Boolean to determine whether the element can delete the stored
    /// derivatives of shape functions w.r.t. global coordinates
    bool Can_delete_dshape_eulerian_stored;

  public:
    /// Constructor, set most storage pointers to NULL.
    // By default the element can delete its own stored shape functions
    StorableShapeElementBase()
      : FiniteElement(),
        Shape_stored_pt(0),
        DShape_local_stored_pt(0),
        D2Shape_local_stored_pt(0),
        Can_delete_shape_local_stored(true),
        DShape_eulerian_stored_pt(0),
        D2Shape_eulerian_stored_pt(0),
        Jacobian_eulerian_stored_pt(0),
        Can_delete_dshape_eulerian_stored(true)
    {
    }

    /// The destructor cleans up the static memory allocated
    /// for shape function storage. Internal and external data get
    /// wiped by the GeneralisedElement destructor; nodes get
    /// killed in mesh destructor.
    virtual ~StorableShapeElementBase();

    /// Broken copy constructor
    StorableShapeElementBase(const StorableShapeElementBase&) = delete;

    /// Broken assignment operator
    void operator=(const StorableShapeElementBase&) = delete;

    /// Delete all the objects stored in the [...]_local_stored_pt
    /// vectors and delete the vectors themselves
    void delete_all_shape_local_stored();

    /// Delete stored shape functions
    void delete_shape_local_stored();

    /// Delete stored derivatives of shape functions w.r.t. to
    /// local coordinates
    void delete_dshape_local_stored();

    /// Delete stored 2nd derivatives of shape functions w.r.t. to
    /// local coordinates
    void delete_d2shape_local_stored();

    /// Delete all storage related to deriatives of shape fcts
    /// w.r.t. to global Eulerian coords
    void delete_all_dshape_eulerian_stored();

    /// Delete stored deriatives of shape fcts w.r.t. to global
    /// Eulerian coords
    void delete_dshape_eulerian_stored();

    /// Delete stored second derivatives of shape functions w.r.t. to
    /// global Eulerian coordinates
    void delete_d2shape_eulerian_stored();

    /// Delete stored Jacobian of mapping between local and global
    /// (Eulerian) coordinates
    void delete_J_eulerian_stored();

    /// Set the spatial integration scheme -- overloaded from the
    /// finite element base class since a change in the integration scheme
    /// forces a recomputation of the shape fcts at the integration points.
    virtual void set_integration_scheme(Integral* const& integral_pt);

    /// Return a pointer to the vector of pointers to the
    /// stored shape functions
    inline Vector<Shape*>*& shape_stored_pt()
    {
      return Shape_stored_pt;
    }

    /// Return a pointer to the vector of pointers to the
    /// stored shape functions (const version)
    inline Vector<Shape*>* const& shape_stored_pt() const
    {
      return Shape_stored_pt;
    }

    /// Return a pointer to the stored shape function at the ipt-th
    /// integration point
    inline Shape*& shape_stored_pt(const unsigned& ipt)
    {
      return (*Shape_stored_pt)[ipt];
    }

    /// Return a pointer to the stored shape function at the ipt-th
    /// integration point (const version)
    inline Shape* const& shape_stored_pt(const unsigned& ipt) const
    {
      return (*Shape_stored_pt)[ipt];
    }

    /// Return a pointer to the vector of pointers to the stored first
    /// derivatives of the shape functions w.r.t the local coordinates
    inline Vector<DShape*>*& dshape_local_stored_pt()
    {
      return DShape_local_stored_pt;
    }

    /// Return a pointer to the vector of pointers to the stored first
    /// derivatives of the shape functions w.r.t the local coordinates
    /// (const version)
    inline Vector<DShape*>* const& dshape_local_stored_pt() const
    {
      return DShape_local_stored_pt;
    }

    /// Return a pointer to the vector of pointers to the stored second
    /// derivatives of the shape functions w.r.t the local coordinates
    inline Vector<DShape*>*& d2shape_local_stored_pt()
    {
      return D2Shape_local_stored_pt;
    }

    /// Return a pointer to the vector of pointers to the stored second
    /// derivatives of the shape functions w.r.t the local coordinates
    /// (const version)
    inline Vector<DShape*>* const& d2shape_local_stored_pt() const
    {
      return D2Shape_local_stored_pt;
    }

    /// Return a pointer to the vector of pointers to the stored first
    /// derivatives of the shape functions w.r.t the global (eulerian)
    /// coordinates
    inline Vector<DShape*>*& dshape_eulerian_stored_pt()
    {
      return DShape_eulerian_stored_pt;
    }

    /// Return a pointer to the vector of pointers to the stored first
    /// derivatives of the shape functions w.r.t the global (eulerian)
    /// coordinates (const version)
    inline Vector<DShape*>* const& dshape_eulerian_stored_pt() const
    {
      return DShape_eulerian_stored_pt;
    }

    /// Return a pointer to the vector of pointers to the stored second
    /// derivatives of the shape functions w.r.t the global (eulerian)
    /// coordinates
    inline Vector<DShape*>*& d2shape_eulerian_stored_pt()
    {
      return D2Shape_eulerian_stored_pt;
    }

    /// Return a pointer to the vector of pointers to the stored second
    /// derivatives of the shape functions w.r.t the global (eulerian)
    /// coordinates (const version)
    inline Vector<DShape*>* const& d2shape_eulerian_stored_pt() const
    {
      return D2Shape_eulerian_stored_pt;
    }

    /// Return a pointer to the vector of Jacobians of
    /// the mapping between the local and global (eulerian) coordinates
    inline Vector<double>*& jacobian_eulerian_stored_pt()
    {
      return Jacobian_eulerian_stored_pt;
    }

    /// Return a pointer to the vector of Jacobians of
    /// the mapping between the local and global (eulerian) coordinates
    /// (const version)
    inline Vector<double>* const& jacobian_eulerian_stored_pt() const
    {
      return Jacobian_eulerian_stored_pt;
    }

    /// Set the shape functions pointed to internally to be
    /// those pointed to by the FiniteElement element_pt (In most
    /// cases all elements of the same type have the same number of
    /// integration points so the shape function values and their
    /// local derivatives are the same --> They only need to be
    /// stored by one element). Calling this function deletes the locally
    /// created storage and re-directs the pointers to the stored
    /// shape function of the specified element.
    void set_shape_local_stored_from_element(
      StorableShapeElementBase* const& element_pt);

    /// Set the derivatives of stored shape functions with respect
    /// to the global coordinates to be the same as
    /// those pointed to by the FiniteElement element_pt. Note that this
    /// function also calls set_shape_local_stored_from_element(element_pt).
    /// so that the local derivatives are also stored.
    /// Calling this function only makes sense for uniformly-spaced meshes with
    /// elements of equal sizes.
    void set_dshape_eulerian_stored_from_element(
      StorableShapeElementBase* const& element_pt);

    /// Calculate the shape functions at the integration points
    /// and store the results internally
    void pre_compute_shape_at_knots();

    /// Return the geometric shape function at the ipt-th integration
    /// point
    void shape_at_knot(const unsigned& ipt, Shape& psi) const;

    /// Calculate the shape functions and first derivatives w.r.t. local
    /// coordinatess at the integration points and store the results internally
    void pre_compute_dshape_local_at_knots();

    /// Return the geometric shape function and its derivative w.r.t.
    /// the local coordinates at the ipt-th integration point. If pre-computed
    /// values have been stored, they will be used.
    void dshape_local_at_knot(const unsigned& ipt,
                              Shape& psi,
                              DShape& dpsids) const;

    /// Calculate the second derivatives of the shape functions
    /// w.r.t. local coordinates at the integration points and store the
    /// results internally
    void pre_compute_d2shape_local_at_knots();

    /// Return the geometric shape function and its first and
    /// second derivatives w.r.t.
    /// the local coordinates at the ipt-th integration point. If pre-computed
    /// values have been stored, they will be used.
    void d2shape_local_at_knot(const unsigned& ipt,
                               Shape& psi,
                               DShape& dpsids,
                               DShape& d2psids) const;

    /// Calculate the Jacobian of the mapping from local to global
    /// coordinates at the integration points and store the results
    /// internally.
    void pre_compute_J_eulerian_at_knots();

    /// Return the Jacobian of the mapping from local to global
    /// coordinates at the ipt-th integration point
    double J_eulerian_at_knot(const unsigned& ipt) const;

    /// Calculate the first derivatives of the shape functions
    /// w.r.t the global coordinates at the integration points and store
    /// the results internally
    void pre_compute_dshape_eulerian_at_knots();

    /// Return the geometric shape functions and also first
    /// derivatives w.r.t. global coordinates at the ipt-th integration point.
    /// If the values of the shape functions and derivatives have been
    /// pre-computed, these will be used
    double dshape_eulerian_at_knot(const unsigned& ipt,
                                   Shape& psi,
                                   DShape& dpsidx) const;


    /// Calculate the first and second derivatives of the shape
    /// functions w.r.t global coordinates at the integration points and
    /// store the results internally.
    void pre_compute_d2shape_eulerian_at_knots();

    /// Return the geometric shape functions and also first
    /// and second derivatives w.r.t. global coordinates at ipt-th integration
    /// point. If the values of the shape functions and derivatives have been
    /// pre-computred, these will be used
    double d2shape_eulerian_at_knot(const unsigned& ipt,
                                    Shape& psi,
                                    DShape& dpsidx,
                                    DShape& d2psidx) const;


    /*  /// Diagnostic */
    /*  void tell_me() */
    /*   { */
    /*    oomph_info << "Diagnostic" << std::endl; */
    /*    oomph_info << Shape_stored_pt << " "; */
    /*    oomph_info << DShape_local_stored_pt << " "; */
    /*    oomph_info << D2Shape_local_stored_pt << " "; */
    /*    oomph_info << DShape_eulerian_stored_pt << " "; */
    /*    oomph_info << D2Shape_eulerian_stored_pt << " "; */
    /*    oomph_info << Jacobian_eulerian_stored_pt << std::endl; */
    /*   } */
  };


  //============================================================================
  /// Base class for solid elements that allow storage of precomputed
  /// shape functions and their derivatives w.r.t to the local and global
  /// (Lagrangian) coordinates at the element's integration points.
  //============================================================================
  class StorableShapeSolidElementBase : public virtual StorableShapeElementBase,
                                        public virtual SolidFiniteElement
  {
  private:
    /// Pointer to storage for the pointers to the derivatives of the
    /// shape functions w.r.t. Lagrangian coordinates at integration points
    Vector<DShape*>* DShape_lagrangian_stored_pt;

    /// Pointer to storage for the pointers to the second derivatives of
    /// the shape functions w.r.t. Lagrangian coordinates at integration points
    Vector<DShape*>* D2Shape_lagrangian_stored_pt;

    /// Pointer to storage for the Jacobian of the mapping between
    /// the local and the global Lagrangian coordinates
    Vector<double>* Jacobian_lagrangian_stored_pt;

    /// Boolean to determine whether the element can delete the stored
    /// shape function derivatives w.r.t. the Lagrangian coordinate
    bool Can_delete_dshape_lagrangian_stored;

  public:
    /// Constructor: Set defaults: Nothing is stored
    StorableShapeSolidElementBase()
      : StorableShapeElementBase(),
        SolidFiniteElement(),
        DShape_lagrangian_stored_pt(0),
        D2Shape_lagrangian_stored_pt(0),
        Jacobian_lagrangian_stored_pt(0),
        Can_delete_dshape_lagrangian_stored(true)
    {
    }

    /// Destructor to clean up any allocated memory
    virtual ~StorableShapeSolidElementBase()
    {
      delete_all_dshape_lagrangian_stored();
    }

    /// Broken copy constructor
    StorableShapeSolidElementBase(const StorableShapeSolidElementBase&) =
      delete;

    /// Broken assignment operator
    void operator=(const StorableShapeSolidElementBase&) = delete;

    /// Delete all the objects stored in the [...]_lagrangian_stored_pt
    /// vectors and delete the vectors themselves
    void delete_all_dshape_lagrangian_stored();

    /// Delete all the objects stored in the
    /// Lagrangian_stored vectors
    void delete_dshape_lagrangian_stored();

    /// Delete stored second derivatives of shape functions w.r.t.
    /// Lagrangian coordinates
    void delete_d2shape_lagrangian_stored();

    /// Delete stored Jaocbian of mapping between local and Lagrangian
    /// coordinates
    void delete_J_lagrangian_stored();

    /// Overload the set_integration_scheme to recompute any stored
    /// derivatives w.r.t. Lagrangian coordinates
    void set_integration_scheme(Integral* const& integral_pt)
    {
      StorableShapeElementBase::set_integration_scheme(integral_pt);

      // If we are storing Lagrangian first and second derivatives, recompute
      // them
      if (D2Shape_lagrangian_stored_pt != 0)
      {
        pre_compute_d2shape_lagrangian_at_knots();
      }
      // If we are storing Lagrangian first derivatives, recompute them
      else if (DShape_lagrangian_stored_pt != 0)
      {
        pre_compute_dshape_lagrangian_at_knots();
      }
    }

    /// Return a pointer to the vector of pointers to the stored first
    /// derivatives of the shape functions w.r.t the global (eulerian)
    /// coordinates
    inline Vector<DShape*>*& dshape_lagrangian_stored_pt()
    {
      return DShape_lagrangian_stored_pt;
    }

    /// Return a pointer to the vector of pointers to the stored first
    /// derivatives of the shape functions w.r.t the global (eulerian)
    /// coordinates (const version)
    inline Vector<DShape*>* const& dshape_lagrangian_stored_pt() const
    {
      return DShape_lagrangian_stored_pt;
    }

    /// Return a pointer to the vector of pointers to the stored second
    /// derivatives of the shape functions w.r.t the global (eulerian)
    /// coordinates
    inline Vector<DShape*>*& d2shape_lagrangian_stored_pt()
    {
      return D2Shape_lagrangian_stored_pt;
    }

    /// Return a pointer to the vector of pointers to the stored second
    /// derivatives of the shape functions w.r.t the global (eulerian)
    /// coordinates (const version)
    inline Vector<DShape*>* const& d2shape_lagrangian_stored_pt() const
    {
      return D2Shape_lagrangian_stored_pt;
    }

    /// Return a pointer to the vector of Jacobians of
    /// the mapping between the local and global (eulerian) coordinates
    inline Vector<double>*& jacobian_lagrangian_stored_pt()
    {
      return Jacobian_lagrangian_stored_pt;
    }

    /// Return a pointer to the vector of Jacobians of
    /// the mapping between the local and global (eulerian) coordinates
    /// (const version)
    inline Vector<double>* const& jacobian_lagrangian_stored_pt() const
    {
      return Jacobian_lagrangian_stored_pt;
    }


    /// Calculate the first derivatives of the shape functions w.r.t
    /// Lagrangian coordinates at the integration points and store the
    /// results internally.
    void pre_compute_dshape_lagrangian_at_knots();

    /// Return the geometric shape functions and also first
    /// derivatives w.r.t. Lagrangian coordinates at ipt-th integration point.
    /// If the values of the shape function and derivatives have been
    /// pre-computed, they will be used.
    double dshape_lagrangian_at_knot(const unsigned& ipt,
                                     Shape& psi,
                                     DShape& dpsidxi) const;

    /// Calculate the first and second derivatives of the
    /// shape functions w.r.t
    /// Lagrangian coordinates at the integration points and store the
    /// results internally
    void pre_compute_d2shape_lagrangian_at_knots();

    /// Return  the geometric shape functions and also first
    /// and second derivatives w.r.t. Lagrangian coordinates at
    /// the ipt-th integration point. If the values have been pre-computed,
    /// they will be used.
    /// Returns Jacobian of mapping from Lagrangian to local coordinates.
    double d2shape_lagrangian_at_knot(const unsigned& ipt,
                                      Shape& psi,
                                      DShape& dpsidxi,
                                      DShape& d2psidxi) const;


    /// Set the derivatives of stored shape functions with respect
    /// to the lagrangian coordinates to be the same as
    /// those pointed to by the FiniteElement element_pt. Note that this
    /// function also calls set_shape_local_stored_from_element(element_pt).
    /// so that the local derivatives are also stored.
    /// Calling this function only makes sense for uniformly-spaced meshes with
    /// elements of equal sizes.
    void set_dshape_lagrangian_stored_from_element(
      StorableShapeSolidElementBase* const& element_pt);
  };


  //========================================================================
  /// Templated wrapper that attaches the ability to store the shape
  /// functions and their derivatives w.r.t. to the local and global
  /// (Eulerian) coordinates at the integration points to the
  /// element specified by the template parameter.
  //========================================================================
  template<class ELEMENT>
  class StorableShapeElement : public virtual StorableShapeElementBase,
                               public virtual ELEMENT
  {
  public:
    /// Constructor, set most storage pointers to zero
    // By default the element can delete its own stored shape functions
    StorableShapeElement() : StorableShapeElementBase(), ELEMENT()
    {
      // Reset the integration scheme and force a (re)compute of the
      // the stored shape functions and their derivatives w.r.t. to the
      // local coordinates.
      this->set_integration_scheme(this->integral_pt());
    }

    /// Empty virtual destructor
    virtual ~StorableShapeElement() {}

    /// Broken copy constructor
    StorableShapeElement(const StorableShapeElement&) = delete;

    /// Broken assignment operator
    void operator=(const StorableShapeElement&) = delete;
  };


  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////


  //============================================================================
  /// Templated wrapper that attaches the ability to store the shape
  /// functions and their derivatives w.r.t. to the local and global
  /// (Eulerian) coordinates at the integration points to the
  /// SolidFiniteElement specified by the template parameter.
  //============================================================================
  template<class ELEMENT>
  class StorableShapeSolidElement
    : public virtual StorableShapeSolidElementBase,
      public virtual ELEMENT
  {
  public:
    /// Constructor: Set defaults
    StorableShapeSolidElement() : StorableShapeSolidElementBase(), ELEMENT()
    {
      // Reset the integration scheme and force a (re-)computation
      // of the shape functions and their derivatives w.r.t. to the
      // local coordinates at the integration points.
      this->set_integration_scheme(this->integral_pt());
    }

    /// Destructor to clean up any allocated memory
    virtual ~StorableShapeSolidElement() {}


    /// Broken copy constructor
    StorableShapeSolidElement(const StorableShapeSolidElement&) = delete;

    /// Broken assignment operator
    void operator=(const StorableShapeSolidElement&) = delete;
  };

} // namespace oomph

#endif
