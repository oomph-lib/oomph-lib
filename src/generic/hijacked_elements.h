// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
// Header file for hijacked elements

// Include guards to prevent multiple inclusion of the header
#ifndef OOMPH_HIJACKED_ELEMENTS_HEADER
#define OOMPH_HIJACKED_ELEMENTS_HEADER


// oomph-lib header
#include "elements.h"
#include "spines.h"

namespace oomph
{
  //========================================================================
  /// HijackedElement base class that provides storage and access funcitons
  /// for pointers to the global equation numbers that are hijacked by
  /// the HijackedElement. A default residuals multiplier is also provided.
  //========================================================================
  class HijackedElementBase
  {
  public:
    /// Constructor, initialise the pointer to the equation numbers
    /// for the storage to zero
    HijackedElementBase()
      : Hijacked_global_eqn_number_pt(0), Hijacked_local_eqn_number_pt(0)
    {
      // Set the default value of the pointer to the residuals multiplier.
      // The default value is zero, so that the element is "completely hijacked"
      // and the exisiting contributions to the residuals and jacobian are
      // totally wiped
      Residual_multiplier_pt = &Default_residual_multiplier;
    }

    /// Destructor, destroy the storage for the equation numbers
    virtual ~HijackedElementBase();

    /// Reset the hijacked data pt, so that none of the equations in the
    /// element are hijacked.
    void unhijack_all_data()
    {
      Hijacked_global_eqn_number_pt->clear();
    }

    /// Return the value of the residual multiplier
    inline const double& residual_multiplier() const
    {
      return *Residual_multiplier_pt;
    }

    /// Return the pointer to the residual multiplier
    inline double*& residual_multiplier_pt()
    {
      return Residual_multiplier_pt;
    }

  protected:
    /// Pointer to a Set of pointers to the equation numbers that will be
    /// hijacked by this element. Note that these MUST be pointers because
    /// hijacking information is set BEFORE global equation numbers have
    /// been assigned.
    std::set<long*>* Hijacked_global_eqn_number_pt;

    /// Pointer to  a vector of integers
    /// containing the local equation numbers of
    /// any hijacked variables in the element.
    Vector<int>* Hijacked_local_eqn_number_pt;

    /// Pointer to a double that multiplies the contribution to
    /// the residuals from the original element.
    /// This is usually used as a homotopy parameter to permit a smooth
    /// transition between different types of boundary conditions, rather than
    /// switching them on or off abruptly
    double* Residual_multiplier_pt;

    /// Static default value for the double that multiplies the original
    /// residuals
    static double Default_residual_multiplier;

    /// Mark the global equation, addressed by global_eqn_pt,
    /// as hijacked by this element.
    void hijack_global_eqn(long* const& global_eqn_pt);

    /// The global equation, addressed by global_eqn_pt,
    /// is no longer hijacked by this element.
    void unhijack_global_eqn(long* const& global_eqn_pt);
  };


  //========================================================================
  /// Hijacked elements are elements in which one or more
  /// Data values that affect the element's residuals, are determined
  /// by another element -- the data values are then said to have
  /// been hijacked by another element. The main functionality
  /// added by the Hijacked element class is that it wipes out
  /// those entries in the element's residual vector and those rows in
  /// the element's Jacobian matrix that are determined by the "other"
  /// elements that have hijacked the values. Note that for continuation
  /// in homotopy parameters, it may be desriable to multiply the residuals
  /// and corresponding jacobian entries by a "homotopy parameter". The
  /// value of this parameter can be set by assigning residual_multiplier_pt()
  /// which has a default value of zero. Note: it would be possible to extend
  /// the functionality so that different residuals are multiplied by
  /// different values, but will this ever be required?
  //========================================================================
  template<class ELEMENT>
  class Hijacked : public virtual ELEMENT, public virtual HijackedElementBase
  {
  public:
    /// Constructor, call the constructors of the base elements
    Hijacked() : ELEMENT(), HijackedElementBase() {}

    /// Constructor used for hijacking face elements
    Hijacked(FiniteElement* const& element_pt, const int& face_index)
      : ELEMENT(element_pt, face_index), HijackedElementBase()
    {
    }


    /// Constructor used for hijacking face elements with specification
    /// of ID of additional variables
    Hijacked(FiniteElement* const& element_pt,
             const int& face_index,
             const unsigned& id = 0)
      : ELEMENT(element_pt, face_index, id), HijackedElementBase()
    {
    }


    /// Hijack the i-th value stored at internal data n.
    /// Optionally return a custom-made (copied) data object that contains only
    /// the hijacked value. This can be used as the input to other elements.
    /// Note that the calling program assumes responsibility for this
    /// data object and must clean it up. The default is that
    /// the data object is returned
    Data* hijack_internal_value(const unsigned& n,
                                const unsigned& i,
                                const bool& return_data = true)
    {
      // Initialise pointer to zero
      Data* temp_data_pt = 0;

      // If desired,
      // Create a new Data object containing only the value that is to be
      // hijacked
      if (return_data)
      {
        temp_data_pt = new HijackedData(i, this->internal_data_pt(n));
      }

      // Mark the value as hijacked
      hijack_global_eqn(this->internal_data_pt(n)->eqn_number_pt(i));

      // Return the pointer to the data
      return temp_data_pt;
    }


    /// Hijack the i-th value stored at external data n.
    /// Optionally return a custom-made (copied) data object that contains only
    /// the hijacked value. Note that the calling program assumes
    /// responsibility for this data object and must clean it up.
    /// The default is that the data object is returned
    Data* hijack_external_value(const unsigned& n,
                                const unsigned& i,
                                const bool& return_data = true)
    {
      // Initialise pointer to zero
      Data* temp_data_pt = 0;
      // If desired
      // create a new Data object containing only the value that is to be
      // hijacked
      if (return_data)
      {
        temp_data_pt = new HijackedData(i, this->external_data_pt(n));
      }

      // Mark the value as hijacked
      hijack_global_eqn(this->external_data_pt(n)->eqn_number_pt(i));

      // Return the pointer to the data
      return temp_data_pt;
    }

    /// Hijack the i-th value stored at node n.
    /// Optionally return a custom-made (copied) data object that contains only
    /// the hijacked value. Once again, the calling program must
    /// clean up the allocated Data object.
    /// The default is that the data object is returned
    Data* hijack_nodal_value(const unsigned& n,
                             const unsigned& i,
                             const bool& return_data = true)
    {
      // Initialise pointer to zero
      Data* temp_data_pt = 0;
      // If desired
      // create a new Data object containing only the value that is to be
      // hijacked
      if (return_data)
      {
        temp_data_pt = new HijackedData(i, this->node_pt(n));
      }

      // Mark the value as hijacked
      hijack_global_eqn(this->node_pt(n)->eqn_number_pt(i));

      // Return the pointer to the data, which may be null
      return temp_data_pt;
    }

    /// Hijack the i-th positional value stored at node n.
    /// Optionaly return a custom-made (copied) data object that contains only
    /// the hijacked value. Again, responsibility for the memory allocated
    /// lies with the calling function.
    /// The default is that the data object is returned
    Data* hijack_nodal_position_value(const unsigned& n,
                                      const unsigned& i,
                                      const bool& return_data = true)
    {
      // Can we do the casting?
      SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(this->node_pt(n));
      if (solid_node_pt == 0)
      {
        std::string error_message = "Failed to cast to SolidNode\n ";
        error_message += "You may be trying to hijack a non-elastic element\n";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Initialise pointer to zero
      Data* temp_data_pt = 0;
      // If desired
      // create a new Data object containing only the value that is to be
      // hijacked
      if (return_data)
      {
        temp_data_pt =
          new HijackedData(i, solid_node_pt->variable_position_pt());
      }

      // Mark the value as hijacked
      hijack_global_eqn(
        solid_node_pt->variable_position_pt()->eqn_number_pt(i));

      // Return the pointer to the data
      return temp_data_pt;
    }

    /// Hijack the i-th value stored at the spine that affects
    /// local node n.
    /// Optionally return a custom-made (copied) data object that contains only
    /// the hijacked value. Deletion must be handled at the higher level
    /// The default is that the data object is returned
    Data* hijack_nodal_spine_value(const unsigned& n,
                                   const unsigned& i,
                                   const bool& return_data = true)
    {
      // Can we actually do this casting
      SpineNode* spine_node_pt = dynamic_cast<SpineNode*>(this->node_pt(n));
      if (spine_node_pt == 0)
      {
        std::string error_message = "Failed to cast to SpineNode\n ";
        error_message += "You may be trying to hijack a non-spine element\n";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Initialise pointer to zero
      Data* temp_data_pt = 0;
      // If desired
      // create a new Data object containing only the value that is to be
      // hijacked
      if (return_data)
      {
        temp_data_pt =
          new HijackedData(i, spine_node_pt->spine_pt()->spine_height_pt());
      }

      // Mark the data as hijacked
      hijack_global_eqn(
        spine_node_pt->spine_pt()->spine_height_pt()->eqn_number_pt(i));

      // Return the pointer to the data
      return temp_data_pt;
    }


    /// Set up the local equation numbers for the underlying element,
    /// then set up the local arrays to hold the hijacked variables.
    /// If the boolean argument is true then pointers to the associated degrees
    /// of freedom are stored in the array Dof_pt
    void assign_local_eqn_numbers(const bool& store_local_dof_pt)
    {
      // If things have already been allocated,
      // clear the local hijacked array, so that if the equation numbers
      // are reassigned after changes in the boundary conditions, the
      // correct terms will be in the array.
      if (Hijacked_local_eqn_number_pt != 0)
      {
        Hijacked_local_eqn_number_pt->clear();
      }
      // Otherwise allocate it
      else
      {
        Hijacked_local_eqn_number_pt = new Vector<int>;
      }


      // Call the hijacked element's assign_local_eqn_numbers
      ELEMENT::assign_local_eqn_numbers(store_local_dof_pt);

      // If any values have been hijacked, we need to find the corresponding
      // local equation numbers
      if (Hijacked_global_eqn_number_pt != 0)
      {
        // Now loop over the local array that stores GLOBAL equation numbers
        // to check if any of the local values have in fact been hijacked
        // by "somebody"
        unsigned n_dof = this->ndof();
        for (unsigned i = 0; i < n_dof; i++)
        {
          // Loop over *all* hijacked data
          for (std::set<long*>::iterator it =
                 Hijacked_global_eqn_number_pt->begin();
               it != Hijacked_global_eqn_number_pt->end();
               ++it)
          {
            // If the hijacked (global!) equation is not pinned
            if (*(*it) >= 0)
            {
              // Get the (unsigned) global equation number:
              // Note that it is only unsigned AFTER the test above.
              unsigned long hijacked_eqn_number = *(*it);

              // If a GLOBAL variable used by this element is hijacked add the
              // variable's LOCAL index to the array
              if (hijacked_eqn_number == this->eqn_number(i))
              {
                Hijacked_local_eqn_number_pt->push_back(i);
                break;
              }
            }
          }
        }
      }
    }

    /// Get the residuals from the underlying element, but then wipe the
    /// entries in the residual vector that correspond to hijacked
    /// values -- they will be computed by other elements.
    void get_residuals(Vector<double>& residuals)
    {
      // Get parent redisuals
      ELEMENT::get_residuals(residuals);
      // Multiply any hijacked dofs by the residual multiplier
      //(default value is zero)
      unsigned n_hijacked = Hijacked_local_eqn_number_pt->size();
      for (unsigned i = 0; i < n_hijacked; i++)
      {
        residuals[(*Hijacked_local_eqn_number_pt)[i]] *= residual_multiplier();
      }
    }


    /// Get the residuals and Jacobian matrix from the underlying
    /// element, but then wipe the entries in the residual vector and the
    /// rows in the Jacobian matrix that correspond to hijacked
    /// values -- they will be computed by other elements.
    void get_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian)
    {
      // Call the element's get jacobian function
      ELEMENT::get_jacobian(residuals, jacobian);
      // Wipe any hijacked dofs
      unsigned n_hijacked = Hijacked_local_eqn_number_pt->size();
      unsigned n_dof = this->ndof();
      for (unsigned i = 0; i < n_hijacked; i++)
      {
        // Multiply any hijacked dofs by the residual multiplier
        //(default value is zero)
        residuals[(*Hijacked_local_eqn_number_pt)[i]] *= residual_multiplier();
        // Multiply the row in the Jacobian matrix by the residual
        // multiplier for consistency
        for (unsigned j = 0; j < n_dof; j++)
        {
          jacobian((*Hijacked_local_eqn_number_pt)[i], j) *=
            residual_multiplier();
        }
      }
    }
  };


  //============================================================================
  /// Explicit definition of the face geometry of hijacked elements:
  /// the same as the face geometry of the underlying element
  //============================================================================
  template<class ELEMENT>
  class FaceGeometry<Hijacked<ELEMENT>> : public virtual FaceGeometry<ELEMENT>
  {
  public:
    /// Constructor
    FaceGeometry() : FaceGeometry<ELEMENT>() {}

  protected:
  };

  //============================================================================
  /// Explicit definition of the face geometry of hijacked elements:
  /// the same as the face geometry of the underlying element
  //============================================================================
  template<class ELEMENT>
  class FaceGeometry<FaceGeometry<Hijacked<ELEMENT>>>
    : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    /// Constructor
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}

  protected:
  };

  //============================================================================
  /// Explicit definition of the face geometry of hijacked elements:
  /// the same as the face geometry of the underlying element
  //============================================================================
  template<class ELEMENT>
  class FaceGeometry<Hijacked<FaceGeometry<ELEMENT>>>
    : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    /// Constructor
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}

  protected:
  };


} // namespace oomph

#endif
