// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2026 Matthias Heil and Andrew Hazel
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

// Include guard to prevent multiple inclusions of this header
#ifndef OOMPH_REFINEABLE_AXISYM_CYLINDRICAL_ELASTICITY_WITH_PRESSURE_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_AXISYM_CYLINDRICAL_ELASTICITY_WITH_PRESSURE_ELEMENTS_HEADER

// oomph-lib headers
#include "axisym_cylindrical_solid_with_pressure_elements.h"
#include "../generic/refineable_quad_element.h"
#include "../generic/refineable_brick_element.h"
#include "../generic/error_estimator.h"

namespace oomph
{
  //========================================================================
  /// Class for Refineable axisymmetric PVD equations in cylindrical coords
  //========================================================================
  class RefineableAxisymmetricCylindricalPVDWithPressureEquations
    : public virtual AxisymmetricCylindricalPVDWithPressureEquations,
      public virtual RefineableSolidElement,
      public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// Constructor
    RefineableAxisymmetricCylindricalPVDWithPressureEquations()
      : AxisymmetricCylindricalPVDWithPressureEquations(),
        RefineableElement(),
        RefineableSolidElement(),
        ElementWithZ2ErrorEstimator()
    {
    }

    /// Call the residuals including hanging node cases
    void fill_in_contribution_to_residuals_axisym_pvd_with_pressure(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag);

    /// No values are interpolated in this element (pure solid)
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      values.clear();
    }

    /// No values are interpolated in this element (pure solid)
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      values.clear();
    }

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      return 4;
    }

    /// Get 'flux' for Z2 error recovery:   Upper triangular entries
    /// in strain tensor.
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
#ifdef PARANOID
      unsigned num_entries = 4;
      if (flux.size() != num_entries)
      {
        std::ostringstream error_message;
        error_message << "The flux vector has the wrong number of entries, "
                      << flux.size() << ", whereas it should be " << num_entries
                      << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get strain matrix
      DenseMatrix<double> strain(3);
      this->get_strain(s, strain);

      // Pack into flux Vector
      unsigned icount = 0;

      // Start with diagonal terms
      for (unsigned i = 0; i < 3; i++)
      {
        flux[icount] = strain(i, i);
        icount++;
      }

      // Only non-zero off-diagonal
      flux[icount] = strain(0, 1);
      icount++;
    }

    /// Number of continuously interpolated values: 0 (pure solid problem)
    unsigned ncont_interpolated_values() const
    {
      return 0;
    }

    /// Return a pointer to the solid node at which pressure dof l2 is stored
    virtual Node* solid_pressure_node_pt(const unsigned& l)
    {
      return 0;
    }

    /// Further build function, pass the pointers down to the sons
    void further_build()
    {
      RefineableAxisymmetricCylindricalPVDWithPressureEquations*
        cast_father_element_pt = dynamic_cast<
          RefineableAxisymmetricCylindricalPVDWithPressureEquations*>(
          this->father_element_pt());

      // Do whatever needs to be done in the base class
      RefineableSolidElement::further_build();

      // Set pointer to body force function
      this->Body_force_fct_pt = cast_father_element_pt->body_force_fct_pt();

      // Set pointer to the contitutive law
      this->Constitutive_law_pt = cast_father_element_pt->constitutive_law_pt();

      // Set the timescale ratio (non-dim. density)
      this->Lambda_sq_pt = cast_father_element_pt->lambda_sq_pt();

      // Set the mass damping parameter
      this->Eta_mass_pt = cast_father_element_pt->eta_mass_pt();

      // Set the compressibility/incompressibility status
      this->Incompressible = cast_father_element_pt->is_incompressible();
    }
  };

  //========================================================================
  /// Class for refineable QPVDElement elements
  //========================================================================
  class RefineableQAxisymmetricCylindricalPVDWithPressureElement
    : public virtual QAxisymmetricCylindricalPVDWithPressureElement,
      public virtual RefineableAxisymmetricCylindricalPVDWithPressureEquations,
      public virtual RefineableSolidQElement<2>
  {
  private:
    /// Unpin all solid pressure dofs
    void unpin_elemental_solid_pressure_dofs()
    {
      // find the index at which the pressure is stored
      int solid_p_index = this->solid_p_nodal_index();
      unsigned n_node = this->nnode();
      // loop over nodes
      for (unsigned n = 0; n < n_node; n++)
      {
        this->node_pt(n)->unpin(solid_p_index);
      }
    }

  public:
    /// Constructor:
    RefineableQAxisymmetricCylindricalPVDWithPressureElement()
      : QAxisymmetricCylindricalPVDWithPressureElement(),
        RefineableElement(),
        RefineableSolidElement(),
        RefineableAxisymmetricCylindricalPVDWithPressureEquations(),
        RefineableSolidQElement<2>()
    {
    }

    /// Overload the number of additional solid dofs at each node, we
    /// shall always assign 1, otherwise it's a real pain
    unsigned required_nvalue(const unsigned& n) const
    {
      return 1;
    }

    /// Number of continuously interpolated values (1) solid pressure
    unsigned ncont_interpolated_values() const
    {
      return 1;
    }

    /// Interpolate the solid pressures
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      // There is only one solid pressure, initialise to zero
      values.resize(1);

      // Get the interpolated value
      values[0] = this->interpolated_solid_p(s);
    }

    /// The time-dependent verion
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      // There is only one solid pressure, initialise to zero
      values.resize(1);
      // The solid pressure does not depend on time!
      values[0] = this->interpolated_solid_p(s);
    }

    /// Empty rebuild from sons, no need to reconstruct anything here
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// Pin the redundant solid pressure
    void pin_elemental_redundant_nodal_solid_pressures()
    {
      // Find the index of the solid pressure
      int solid_p_index = this->solid_p_nodal_index();
      // Let's pin all pressure nodes
      unsigned n_node = this->nnode();
      for (unsigned l = 0; l < n_node; l++)
      {
        // Pin the solid pressure
        this->node_pt(l)->pin(solid_p_index);
      }

      // Now loop over the pressure nodes and unpin the solid pressures
      unsigned n_solid_pres = this->nsolid_pres();
      // Loop over these nodes and unpin the solid pressures
      for (unsigned l = 0; l < n_solid_pres; l++)
      {
        Node* nod_pt = this->solid_pressure_node_pt(l);
        if (!nod_pt->is_hanging(solid_p_index))
        {
          nod_pt->unpin(solid_p_index);
        }
      }
    }

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QAxisymmetricCylindricalPVDWithPressureElement::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QAxisymmetricCylindricalPVDWithPressureElement::vertex_node_pt(j);
    }

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return 2;
    }

    /// The pressure "nodes" are a
    /// subset of the nodes, so when value_id==0, the n-th pressure
    /// node is returned.
    Node* interpolating_node_pt(const unsigned& n, const int& value_id)

    {
#ifdef PARANOID
      RefineableElement::check_value_id(1, value_id);
#endif
      // If we are at the value return the solid pressure node
      if (value_id == 0)
      {
        return this->solid_pressure_node_pt(n);
      }
      // Otherwise return the nodal values
      else
      {
        return this->node_pt(n);
      }
    }

    /// The pressure nodes are the corner nodes, so when value_id==0,
    /// the fraction is the same as the 1d node number, 0 or 1.
    double local_one_d_fraction_of_interpolating_node(const unsigned& n1d,
                                                      const unsigned& i,
                                                      const int& value_id)
    {
#ifdef PARANOID
      RefineableElement::check_value_id(1, value_id);
#endif
      // If it's the only value, we have the pressure
      if (value_id == 0)
      {
        // The pressure nodes are just located on the boundaries at 0 or 1
        return double(n1d);
      }
      // Otherwise we have the geometric nodes
      else
      {
        return this->local_one_d_fraction_of_node(n1d, i);
      }
    }

    /// The velocity nodes are the same as the geometric nodes. The
    /// pressure nodes must be calculated by using the same methods as
    /// the geometric nodes, but by recalling that there are only two pressure
    /// nodes per edge.
    Node* get_interpolating_node_at_local_coordinate(const Vector<double>& s,
                                                     const int& value_id)
    {
#ifdef PARANOID
      RefineableElement::check_value_id(1, value_id);
#endif

      // If we are calculating solid pressure nodes
      if (value_id == 0)
      {
        // Storage for the index of the pressure node
        unsigned total_index = 0;
        // The number of nodes along each 1d edge is 2.
        unsigned NNODE_1D = 2;
        // Storage for the index along each boundary
        Vector<int> index(2);
        // Loop over the coordinates
        for (unsigned i = 0; i < 2; i++)
        {
          // If we are at the lower limit, the index is zero
          if (s[i] == -1.0)
          {
            index[i] = 0;
          }
          // If we are at the upper limit, the index is the number of nodes
          // minus 1
          else if (s[i] == 1.0)
          {
            index[i] = NNODE_1D - 1;
          }
          // Otherwise, we have to calculate the index in general
          else
          {
            // For uniformly spaced nodes the 0th node number would be
            double float_index = 0.5 * (1.0 + s[i]) * (NNODE_1D - 1);
            index[i] = int(float_index);
            // What is the excess. This should be safe because the
            // taking the integer part rounds down
            double excess = float_index - index[i];
            // If the excess is bigger than our tolerance there is no node,
            // return null
            if ((excess > FiniteElement::Node_location_tolerance) &&
                ((1.0 - excess) > FiniteElement::Node_location_tolerance))
            {
              return 0;
            }
          }
          /// Construct the general pressure index from the components.
          total_index +=
            index[i] * static_cast<unsigned>(pow(static_cast<float>(NNODE_1D),
                                                 static_cast<int>(i)));
        }
        // If we've got here we have a node, so let's return a pointer to it
        return this->solid_pressure_node_pt(total_index);
      }
      // Otherwise velocity nodes are the same as pressure nodes
      else
      {
        return this->get_node_at_local_coordinate(s);
      }
    }


    /// The number of 1d pressure nodes is 2, otherwise we have
    /// the positional nodes
    unsigned ninterpolating_node_1d(const int& value_id)
    {
#ifdef PARANOID
      RefineableElement::check_value_id(1, value_id);
#endif

      if (value_id == 0)
      {
        return 2;
      }
      else
      {
        return this->nnode_1d();
      }
    }

    /// The number of pressure nodes is 4. The number of
    /// velocity nodes is the same as the number of geometric nodes.
    unsigned ninterpolating_node(const int& value_id)
    {
#ifdef PARANOID
      RefineableElement::check_value_id(1, value_id);
#endif

      if (value_id == 0)
      {
        return 4;
      }
      else
      {
        return this->nnode();
      }
    }

    /// The basis interpolating the pressure is given by pshape().
    /// / The basis interpolating the velocity is shape().
    void interpolating_basis(const Vector<double>& s,
                             Shape& psi,
                             const int& value_id) const
    {
#ifdef PARANOID
      RefineableElement::check_value_id(1, value_id);
#endif

      if (value_id == 0)
      {
        return this->solid_pshape(s, psi);
      }
      else
      {
        return this->shape(s, psi);
      }
    }

    ///  Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes.
    void further_setup_hanging_nodes()
    {
      this->setup_hang_for_value(this->solid_p_nodal_index());
    }

    // Return a pointer to the solid node at which pressure dof l2 is stored
    Node* solid_pressure_node_pt(const unsigned& l)
    {
      return this->node_pt(this->Pconv[l]);
    }
  };

  //==============================================================
  /// FaceGeometry of the
  /// RefineableQAxisymmetricCylindricalPVDWithPressureElement
  //==============================================================
  template<>
  class FaceGeometry<RefineableQAxisymmetricCylindricalPVDWithPressureElement>
    : public virtual SolidQElement<1, 3>
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<1, 3>() {}
  };

  //==============================================================
  /// FaceGeometry of the FaceGeometry of the
  /// RefineableQAxisymmetricCylindricalPVDWithPressureElement
  //==============================================================
  template<>
  class FaceGeometry<
    FaceGeometry<RefineableQAxisymmetricCylindricalPVDWithPressureElement>>
    : public virtual PointElement
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : PointElement() {}
  };

} // namespace oomph

#endif
