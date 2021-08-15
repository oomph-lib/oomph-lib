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
// Header file for refineable linearised axisymmetric Navier-Stokes elements

#ifndef OOMPH_REFINEABLE_LINEARISED_AXISYM_NAVIER_STOKES_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_LINEARISED_AXISYM_NAVIER_STOKES_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// oomph-lib includes
#include "../generic/refineable_quad_element.h"
#include "../generic/error_estimator.h"
#include "linearised_axisym_navier_stokes_elements.h"

namespace oomph
{
  //=======================================================================
  /// \short Refineable version of the linearised axisymmetric
  /// Navier--Stokes equations
  //=======================================================================
  class RefineableLinearisedAxisymmetricNavierStokesEquations
    : public virtual LinearisedAxisymmetricNavierStokesEquations,
      public virtual RefineableElement,
      public virtual ElementWithZ2ErrorEstimator
  {
  protected:
    /// \short Pointer to n_p-th pressure node (Default: NULL,
    /// indicating that pressure is not based on nodal interpolation).
    virtual Node* pressure_node_pt(const unsigned& n_p)
    {
      return NULL;
    }

    /// \short Unpin all pressure dofs in the element
    virtual void unpin_elemental_pressure_dofs() = 0;

    /// \short Pin unused nodal pressure dofs (empty by default, because
    /// by default pressure dofs are not associated with nodes)
    virtual void pin_elemental_redundant_nodal_pressure_dofs() {}

  public:
    /// \short Empty Constructor
    RefineableLinearisedAxisymmetricNavierStokesEquations()
      : LinearisedAxisymmetricNavierStokesEquations(),
        RefineableElement(),
        ElementWithZ2ErrorEstimator()
    {
    }

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      // 3 diagonal strain rates, 3 off diagonal
      return 6;
    }

    /// \short Get 'flux' for Z2 error recovery:   Upper triangular entries
    /// in strain rate tensor.
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      // Specify the number of velocity dimensions
      unsigned DIM = 3;

#ifdef PARANOID
      unsigned num_entries = DIM + ((DIM * DIM) - DIM) / 2;
      if (flux.size() != num_entries)
      {
        std::ostringstream error_message;
        error_message << "The flux vector has the wrong number of entries, "
                      << flux.size() << ", whereas it should be " << num_entries
                      << std::endl;
        throw OomphLibError(error_message.str(),
                            "RefineableLinearisedAxisymmetricNavierStokesEquati"
                            "ons::get_Z2_flux()",
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get strain rate matrix
      DenseMatrix<double> strainrate(DIM);
      this->strain_rate(s, strainrate);

      // Pack into flux Vector
      unsigned icount = 0;

      // Start with diagonal terms
      for (unsigned i = 0; i < DIM; i++)
      {
        flux[icount] = strainrate(i, i);
        icount++;
      }

      // Off diagonals row by row
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = i + 1; j < DIM; j++)
        {
          flux[icount] = strainrate(i, j);
          icount++;
        }
      }
    }

    /// Fill in the geometric Jacobian, which in this case is r
    double geometric_jacobian(const Vector<double>& x)
    {
      return x[0];
    }

    ///  Further build: pass the pointers down to the sons
    void further_build()
    {
      // Find the father element
      RefineableLinearisedAxisymmetricNavierStokesEquations*
        cast_father_element_pt =
          dynamic_cast<RefineableLinearisedAxisymmetricNavierStokesEquations*>(
            this->father_element_pt());

      // Set the viscosity ratio pointer
      this->Viscosity_Ratio_pt = cast_father_element_pt->viscosity_ratio_pt();

      // Set the density ratio pointer
      this->Density_Ratio_pt = cast_father_element_pt->density_ratio_pt();

      // Set pointer to global Reynolds number
      this->Re_pt = cast_father_element_pt->re_pt();

      // Set pointer to global Reynolds number x Strouhal number (=Womersley)
      this->ReSt_pt = cast_father_element_pt->re_st_pt();

      // Set pointer to azimuthal mode number
      this->Azimuthal_Mode_Number_pt =
        cast_father_element_pt->azimuthal_mode_number_pt();

      // Set the ALE_is_disabled flag
      this->ALE_is_disabled = cast_father_element_pt->ALE_is_disabled;
    }

    /// \short Loop over all elements in Vector (which typically contains
    /// all the elements in a fluid mesh) and pin the nodal pressure degrees
    /// of freedom that are not being used. Function uses
    /// the member function
    /// - \c RefineableLinearisedAxisymmetricNavierStokesEquations::
    ///       pin_all_nodal_pressure_dofs()
    /// .
    /// which is empty by default and should be implemented for
    /// elements with nodal pressure degrees of freedom
    /// (e.g. for refineable Taylor-Hood.)
    static void pin_redundant_nodal_pressures(
      const Vector<GeneralisedElement*>& element_pt)
    {
      // Loop over all elements to brutally pin all nodal pressure degrees of
      // freedom
      const unsigned n_element = element_pt.size();
      for (unsigned e = 0; e < n_element; e++)
      {
        dynamic_cast<RefineableLinearisedAxisymmetricNavierStokesEquations*>(
          element_pt[e])
          ->pin_elemental_redundant_nodal_pressure_dofs();
      }
    }

    /// Unpin all pressure dofs in elements listed in Vector
    static void unpin_all_pressure_dofs(
      const Vector<GeneralisedElement*>& element_pt)
    {
      // Loop over all elements to brutally unpin all nodal pressure degrees of
      // freedom and internal pressure degrees of freedom
      const unsigned n_element = element_pt.size();
      for (unsigned e = 0; e < n_element; e++)
      {
        dynamic_cast<RefineableLinearisedAxisymmetricNavierStokesEquations*>(
          element_pt[e])
          ->unpin_elemental_pressure_dofs();
      }
    }


  private:
    /// \short Add element's contribution to the elemental residual vector
    /// and/or Jacobian matrix
    /// flag=1: compute both
    /// flag=0: compute only residual vector
    void fill_in_generic_residual_contribution_linearised_axi_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag);

  }; // End of RefineableLinearisedAxisymmetricNavierStokesEquations class defn


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// \short Refineable version of linearised axisymmetric quadratic
  /// Crouzeix-Raviart elements
  //=======================================================================
  class RefineableLinearisedAxisymmetricQCrouzeixRaviartElement
    : public LinearisedAxisymmetricQCrouzeixRaviartElement,
      public virtual RefineableLinearisedAxisymmetricNavierStokesEquations,
      public virtual RefineableQElement<2>
  {
  private:
    /// Unpin all the internal pressure freedoms
    void unpin_elemental_pressure_dofs()
    {
      const unsigned n_pres = this->npres_linearised_axi_nst();
      // Loop over pressure dofs and unpin
      for (unsigned l = 0; l < n_pres; l++)
      {
        // There are two pressure components
        for (unsigned i = 0; i < 2; i++)
        {
          this->internal_data_pt(P_linearised_axi_nst_internal_index[i])
            ->unpin(l);
        }
      }
    }

  public:
    /// Constructor
    RefineableLinearisedAxisymmetricQCrouzeixRaviartElement()
      : RefineableElement(),
        RefineableLinearisedAxisymmetricNavierStokesEquations(),
        RefineableQElement<2>(),
        LinearisedAxisymmetricQCrouzeixRaviartElement()
    {
    }

    /// Number of continuously interpolated values: 6 (velocities)
    unsigned ncont_interpolated_values() const
    {
      return 6;
    }

    /// Rebuild from sons: Reconstruct pressure from the (merged) sons
    void rebuild_from_sons(Mesh*& mesh_pt)
    {
      using namespace QuadTreeNames;

      // Central pressure value:
      // -----------------------

      // Use average of the sons central pressure values
      // Other options: Take average of the four (discontinuous)
      // pressure values at the father's midpoint]

      Vector<double> av_press(2, 0.0);

      // Loop over the sons
      for (unsigned ison = 0; ison < 4; ison++)
      {
        // Loop over the two pressure components
        for (unsigned i = 0; i < 2; i++)
        {
          // Add the sons midnode pressure
          av_press[i] +=
            quadtree_pt()
              ->son_pt(ison)
              ->object_pt()
              ->internal_data_pt(P_linearised_axi_nst_internal_index[i])
              ->value(0);
        }
      }

      // Use the average
      for (unsigned i = 0; i < 2; i++)
      {
        internal_data_pt(P_linearised_axi_nst_internal_index[i])
          ->set_value(0, 0.25 * av_press[i]);
      }

      Vector<double> slope1(2, 0.0), slope2(2, 0.0);

      // Loop over pressure components
      for (unsigned i = 0; i < 2; i++)
      {
        // Slope in s_0 direction
        // ----------------------

        // Use average of the 2 FD approximations based on the
        // elements central pressure values
        // [Other options: Take average of the four
        // pressure derivatives]

        slope1[i] = quadtree_pt()
                      ->son_pt(SE)
                      ->object_pt()
                      ->internal_data_pt(P_linearised_axi_nst_internal_index[i])
                      ->value(0) -
                    quadtree_pt()
                      ->son_pt(SW)
                      ->object_pt()
                      ->internal_data_pt(P_linearised_axi_nst_internal_index[i])
                      ->value(0);

        slope2[i] = quadtree_pt()
                      ->son_pt(NE)
                      ->object_pt()
                      ->internal_data_pt(P_linearised_axi_nst_internal_index[i])
                      ->value(0) -
                    quadtree_pt()
                      ->son_pt(NW)
                      ->object_pt()
                      ->internal_data_pt(P_linearised_axi_nst_internal_index[i])
                      ->value(0);

        // Use the average
        internal_data_pt(P_linearised_axi_nst_internal_index[i])
          ->set_value(1, 0.5 * (slope1[i] + slope2[i]));

        // Slope in s_1 direction
        // ----------------------

        // Use average of the 2 FD approximations based on the
        // elements central pressure values
        // [Other options: Take average of the four
        // pressure derivatives]

        slope1[i] = quadtree_pt()
                      ->son_pt(NE)
                      ->object_pt()
                      ->internal_data_pt(P_linearised_axi_nst_internal_index[i])
                      ->value(0) -
                    quadtree_pt()
                      ->son_pt(SE)
                      ->object_pt()
                      ->internal_data_pt(P_linearised_axi_nst_internal_index[i])
                      ->value(0);

        slope2[i] = quadtree_pt()
                      ->son_pt(NW)
                      ->object_pt()
                      ->internal_data_pt(P_linearised_axi_nst_internal_index[i])
                      ->value(0) -
                    quadtree_pt()
                      ->son_pt(SW)
                      ->object_pt()
                      ->internal_data_pt(P_linearised_axi_nst_internal_index[i])
                      ->value(0);

        // Use the average
        internal_data_pt(P_linearised_axi_nst_internal_index[i])
          ->set_value(2, 0.5 * (slope1[i] + slope2[i]));
      } // End of loop over pressure components
    }

    /// \short Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return 2;
    }

    /// \short Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return LinearisedAxisymmetricQCrouzeixRaviartElement::nvertex_node();
    }

    /// \short Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return LinearisedAxisymmetricQCrouzeixRaviartElement::vertex_node_pt(j);
    }

    /// \short Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Set the velocity dimension of the element
      const unsigned DIM = 3;

      // Determine size of values Vector: U^C, U^S, W^C, W^S, V^C, V^S
      const unsigned n_values = 2 * DIM;

      // Set size of values Vector and initialise to zero
      values.resize(n_values, 0.0);

      // Calculate velocities: values[0],...
      for (unsigned i = 0; i < n_values; i++)
      {
        values[i] = interpolated_u_linearised_axi_nst(s, i);
      }
    }

    /// \short Get all function values [U^C,U^S,...,P^S] at previous timestep t
    /// (t=0: present; t>0: previous timestep).
    /// \n
    /// Note: Given the generality of the interface (this function is
    /// usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    /// \n
    /// Note: No pressure history is kept, so pressure is always
    /// the current value.
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      const unsigned DIM = 3;

      // Set size of Vector: U^C, U^S, W^C, W^S, V^C, V^S
      values.resize(2 * DIM);

      // Initialise to zero
      for (unsigned i = 0; i < 2 * DIM; i++)
      {
        values[i] = 0.0;
      }

      // Determine number of nodes in element
      const unsigned n_node = nnode();

      // Shape functions
      Shape psif(n_node);
      shape(s, psif);

      // Calculate velocities: values[0],...
      for (unsigned i = 0; i < (2 * DIM); i++)
      {
        // Get the local index at which the i-th velocity is stored
        const unsigned u_local_index = u_index_linearised_axi_nst(i);
        for (unsigned l = 0; l < n_node; l++)
        {
          values[i] += nodal_value(t, l, u_local_index) * psif[l];
        }
      }
    }

    /// \short Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes. Empty
    void further_setup_hanging_nodes() {}

    /// Further build for Crouzeix_Raviart interpolates the internal
    /// pressure dofs from father element: Make sure pressure values and
    /// dp/ds agree between fathers and sons at the midpoints of the son
    /// elements.
    void further_build()
    {
      // Call the generic further build
      RefineableLinearisedAxisymmetricNavierStokesEquations::further_build();

      using namespace QuadTreeNames;

      // What type of son am I? Ask my quadtree representation...
      int son_type = quadtree_pt()->son_type();

      // Pointer to my father (in element impersonation)
      RefineableElement* father_el_pt = quadtree_pt()->father_pt()->object_pt();

      Vector<double> s_father(2);

      // Son midpoint is located at the following coordinates in father element:

      // South west son
      if (son_type == SW)
      {
        s_father[0] = -0.5;
        s_father[1] = -0.5;
      }
      // South east son
      else if (son_type == SE)
      {
        s_father[0] = 0.5;
        s_father[1] = -0.5;
      }
      // North east son
      else if (son_type == NE)
      {
        s_father[0] = 0.5;
        s_father[1] = 0.5;
      }

      // North west son
      else if (son_type == NW)
      {
        s_father[0] = -0.5;
        s_father[1] = 0.5;
      }

      // Pressure values in father element
      // ---------------------------------

      // Find pointer to father element
      RefineableLinearisedAxisymmetricQCrouzeixRaviartElement*
        cast_father_el_pt = dynamic_cast<
          RefineableLinearisedAxisymmetricQCrouzeixRaviartElement*>(
          father_el_pt);

      // Set up storage for pressure in father element
      Vector<double> press(2, 0.0);

      // Loop over pressure components
      for (unsigned i = 0; i < 2; i++)
      {
        // Get pressure from father element
        press[i] =
          cast_father_el_pt->interpolated_p_linearised_axi_nst(s_father, i);

        // Pressure value gets copied straight into internal dof:
        internal_data_pt(P_linearised_axi_nst_internal_index[i])
          ->set_value(0, press[i]);

        // The slopes get copied from father
        internal_data_pt(P_linearised_axi_nst_internal_index[i])
          ->set_value(
            1,
            0.5 * cast_father_el_pt
                    ->internal_data_pt(P_linearised_axi_nst_internal_index[i])
                    ->value(1));

        internal_data_pt(P_linearised_axi_nst_internal_index[i])
          ->set_value(
            2,
            0.5 * cast_father_el_pt
                    ->internal_data_pt(P_linearised_axi_nst_internal_index[i])
                    ->value(2));
      } // End of loop over pressure components
    }

  }; // End of RefineableLinearisedAxisymmetricQCrouzeixRaviartElement defn


  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// \short Face geometry of the refineable linearised axisym
  /// Crouzeix-Raviart elements
  //=======================================================================
  template<>
  class FaceGeometry<RefineableLinearisedAxisymmetricQCrouzeixRaviartElement>
    : public virtual FaceGeometry<LinearisedAxisymmetricQCrouzeixRaviartElement>
  {
  public:
    FaceGeometry()
      : FaceGeometry<LinearisedAxisymmetricQCrouzeixRaviartElement>()
    {
    }
  };


  //=======================================================================
  /// \short Face geometry of face geometric of the refineable linearised
  /// axisym Crouzeix-Raviart elements
  //=======================================================================
  template<>
  class FaceGeometry<
    FaceGeometry<RefineableLinearisedAxisymmetricQCrouzeixRaviartElement>>
    : public virtual FaceGeometry<
        FaceGeometry<LinearisedAxisymmetricQCrouzeixRaviartElement>>
  {
  public:
    FaceGeometry()
      : FaceGeometry<
          FaceGeometry<LinearisedAxisymmetricQCrouzeixRaviartElement>>()
    {
    }
  };


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// \short Refineable version of linearised axisymmetric quadratic
  /// Taylor-Hood elements
  //=======================================================================
  class RefineableLinearisedAxisymmetricQTaylorHoodElement
    : public LinearisedAxisymmetricQTaylorHoodElement,
      public virtual RefineableLinearisedAxisymmetricNavierStokesEquations,
      public virtual RefineableQElement<2>
  {
  private:
    /// Pointer to n_p-th pressure node
    Node* pressure_node_pt(const unsigned& n_p)
    {
      return this->node_pt(this->Pconv[n_p]);
    }

    /// Unpin all pressure dofs
    void unpin_elemental_pressure_dofs()
    {
      // Determine number of nodes in element
      const unsigned n_node = this->nnode();

      // Get nodal indeces of the two pressure components
      Vector<int> p_index(2);
      for (unsigned i = 0; i < 2; i++)
      {
        p_index[i] = this->p_index_linearised_axi_nst(i);
      }

      // Loop over nodes and unpin both pressure components
      for (unsigned n = 0; n < n_node; n++)
      {
        for (unsigned i = 0; i < 2; i++)
        {
          this->node_pt(n)->unpin(p_index[i]);
        }
      }
    }

    ///  Unpin the proper nodal pressure dofs
    void pin_elemental_redundant_nodal_pressure_dofs()
    {
      // Determine number of nodes in element
      const unsigned n_node = this->nnode();

      // Get nodal indeces of the two pressure components
      Vector<int> p_index(2);
      for (unsigned i = 0; i < 2; i++)
      {
        p_index[i] = this->p_index_linearised_axi_nst(i);
      }

      // Loop over all nodes and pin all the nodal pressures
      for (unsigned n = 0; n < n_node; n++)
      {
        for (unsigned i = 0; i < 2; i++)
        {
          this->node_pt(n)->pin(p_index[i]);
        }
      }

      // Loop over all actual pressure nodes and unpin if they're not hanging
      const unsigned n_pres = this->npres_linearised_axi_nst();
      for (unsigned l = 0; l < n_pres; l++)
      {
        Node* nod_pt = this->node_pt(this->Pconv[l]);
        for (unsigned i = 0; i < 2; i++)
        {
          if (!nod_pt->is_hanging(p_index[i]))
          {
            nod_pt->unpin(p_index[i]);
          }
        }
      }
    }

  public:
    /// \short Constructor:
    RefineableLinearisedAxisymmetricQTaylorHoodElement()
      : RefineableElement(),
        RefineableLinearisedAxisymmetricNavierStokesEquations(),
        RefineableQElement<2>(),
        LinearisedAxisymmetricQTaylorHoodElement()
    {
    }

    /// \short Number of values (pinned or dofs) required at node n.
    /// Bumped up to 8 so we don't have to worry if a hanging mid-side node
    /// gets shared by a corner node (which has extra degrees of freedom)
    unsigned required_nvalue(const unsigned& n) const
    {
      return 8;
    }

    /// \short Number of continuously interpolated values: 8
    /// (6 velocities + 2 pressures)
    unsigned ncont_interpolated_values() const
    {
      return 8;
    }

    /// Rebuild from sons: empty
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// \short Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return 2;
    }

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return LinearisedAxisymmetricQTaylorHoodElement::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return LinearisedAxisymmetricQTaylorHoodElement::vertex_node_pt(j);
    }

    /// \short Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Set the velocity dimension of the element
      const unsigned DIM = 3;

      // Determine size of values Vector:
      // U^C, U^S, W^C, W^S, V^C, V^S, P^C, P^S
      const unsigned n_values = 2 * (DIM + 1);

      // Set size of values Vector and initialise to zero
      values.resize(n_values, 0.0);

      // Calculate velocities: values[0],...
      for (unsigned i = 0; i < (2 * DIM); i++)
      {
        values[i] = interpolated_u_linearised_axi_nst(s, i);
      }

      // Calculate pressure: values[DIM], values[DIM+1]
      for (unsigned i = 0; i < 2; i++)
      {
        values[DIM + i] = interpolated_p_linearised_axi_nst(s, i);
      }
    }

    /// \short Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Set the velocity dimension of the element
      const unsigned DIM = 3;

      // Set size of values Vector: U^C, U^S, W^C, W^S, V^C, V^S, P^C, P^S
      values.resize(2 * (DIM + 1));

      // Initialise all entries to zero
      for (unsigned i = 0; i < 2 * (DIM + 1); i++)
      {
        values[i] = 0.0;
      }

      // Determine number of nodes in element
      const unsigned n_node = nnode();

      // Shape functions
      Shape psif(n_node);
      shape(s, psif);

      // Calculate velocities: values[0],...
      for (unsigned i = 0; i < (2 * DIM); i++)
      {
        // Get the nodal coordinate of the velocity
        const unsigned u_nodal_index = u_index_linearised_axi_nst(i);
        for (unsigned l = 0; l < n_node; l++)
        {
          values[i] += nodal_value(t, l, u_nodal_index) * psif[l];
        }
      }

      // Calculate pressure: values[DIM], values[DIM+1]
      // (no history is carried in the pressure)
      for (unsigned i = 0; i < 2; i++)
      {
        values[DIM + i] = interpolated_p_linearised_axi_nst(s, i);
      }
    }

    ///  \short Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes. The two pressure components
    /// are stored at the 6th and 7th location in each node
    void further_setup_hanging_nodes()
    {
      const unsigned DIM = 3;
      for (unsigned i = 0; i < 2; i++)
      {
        this->setup_hang_for_value((2 * DIM) + i);
      }
    }

    /// \short The velocities are isoparametric and so the "nodes"
    /// interpolating the velocities are the geometric nodes. The
    /// pressure "nodes" are a subset of the nodes, so when n_value==6
    /// or 7, the n-th pressure node is returned.
    Node* interpolating_node_pt(const unsigned& n, const int& n_value)

    {
      const int DIM = 3;

      // The only different nodes are the pressure nodes
      if (n_value == (2 * DIM) || n_value == ((2 * DIM) + 1))
      {
        return this->node_pt(this->Pconv[n]);
      }

      // The other variables are interpolated via the usual nodes
      else
      {
        return this->node_pt(n);
      }
    }

    /// \short The pressure nodes are the corner nodes, so when n_value==6
    /// or 7, the fraction is the same as the 1d node number, 0 or 1.
    double local_one_d_fraction_of_interpolating_node(const unsigned& n1d,
                                                      const unsigned& i,
                                                      const int& n_value)
    {
      const int DIM = 3;
      if (n_value == (2 * DIM) || n_value == ((2 * DIM) + 1))
      {
        // The pressure nodes are just located on the boundaries at 0 or 1
        return double(n1d);
      }
      // Otherwise the velocity nodes are the same as the geometric ones
      else
      {
        return this->local_one_d_fraction_of_node(n1d, i);
      }
    }

    /// \short The velocity nodes are the same as the geometric nodes.
    /// The pressure nodes must be calculated by using the same methods
    /// as the geometric nodes, but by recalling that there are only two
    /// pressure nodes per edge.
    Node* get_interpolating_node_at_local_coordinate(const Vector<double>& s,
                                                     const int& n_value)
    {
      const int DIM = 3;

      // If we are calculating pressure nodes
      if (n_value == static_cast<int>(2 * DIM) ||
          n_value == static_cast<int>((2 * DIM) + 1))
      {
        // Storage for the index of the pressure node
        unsigned total_index = 0;
        // The number of nodes along each 1d edge is 2.
        const unsigned NNODE_1D = 2;
        // Storage for the index along each boundary
        // Note that it's only a 2D spatial element
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
        return this->node_pt(this->Pconv[total_index]);
      }
      // Otherwise velocity nodes are the same as pressure nodes
      else
      {
        return this->get_node_at_local_coordinate(s);
      }
    }


    /// \short The number of 1d pressure nodes is 2, the number of 1d
    /// velocity nodes is the same as the number of 1d geometric nodes.
    unsigned ninterpolating_node_1d(const int& n_value)
    {
      const int DIM = 3;
      if (n_value == (2 * DIM) || n_value == ((2 * DIM) + 1))
      {
        return 2;
      }
      else
      {
        return this->nnode_1d();
      }
    }

    /// \short The number of pressure nodes is 4. The number of
    /// velocity nodes is the same as the number of geometric nodes.
    unsigned ninterpolating_node(const int& n_value)
    {
      const int DIM = 3;
      if (n_value == (2 * DIM) || n_value == ((2 * DIM) + 1))
      {
        return 4;
      }
      else
      {
        return this->nnode();
      }
    }

    /// \short The basis interpolating the pressure is given by pshape().
    //// The basis interpolating the velocity is shape().
    void interpolating_basis(const Vector<double>& s,
                             Shape& psi,
                             const int& n_value) const
    {
      const int DIM = 3;
      if (n_value == (2 * DIM) || n_value == ((2 * DIM) + 1))
      {
        return this->pshape_linearised_axi_nst(s, psi);
      }
      else
      {
        return this->shape(s, psi);
      }
    }

  }; // End of RefineableLinearisedAxisymmetricQTaylorHoodElement class defn


  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// \short Face geometry of the refineable linearised axisym
  /// Taylor-Hood elements
  //=======================================================================
  template<>
  class FaceGeometry<RefineableLinearisedAxisymmetricQTaylorHoodElement>
    : public virtual FaceGeometry<LinearisedAxisymmetricQTaylorHoodElement>
  {
  public:
    FaceGeometry() : FaceGeometry<LinearisedAxisymmetricQTaylorHoodElement>() {}
  };


  //=======================================================================
  /// \short Face geometry of face geometric of the refineable linearised
  /// axisym Taylor-Hood elements
  //=======================================================================
  template<>
  class FaceGeometry<
    FaceGeometry<RefineableLinearisedAxisymmetricQTaylorHoodElement>>
    : public virtual FaceGeometry<
        FaceGeometry<LinearisedAxisymmetricQTaylorHoodElement>>
  {
  public:
    FaceGeometry()
      : FaceGeometry<FaceGeometry<LinearisedAxisymmetricQTaylorHoodElement>>()
    {
    }
  };


} // End of namespace oomph

#endif
