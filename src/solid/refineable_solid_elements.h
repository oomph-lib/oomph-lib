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
// Header file for refineable solid mechanics elements

// Include guard to prevent multiple inclusions of this header
#ifndef OOMPH_REFINEABLE_ELASTICITY_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_ELASTICITY_ELEMENTS_HEADER

// oomph-lib headers
#include "solid_elements.h"
#include "../generic/refineable_quad_element.h"
#include "../generic/refineable_brick_element.h"
#include "../generic/error_estimator.h"

namespace oomph
{
  //========================================================================
  /// Class for Refineable PVD equations
  //========================================================================
  template<unsigned DIM>
  class RefineablePVDEquations : public virtual PVDEquations<DIM>,
                                 public virtual RefineableSolidElement,
                                 public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// Constructor
    RefineablePVDEquations()
      : PVDEquations<DIM>(),
        RefineableElement(),
        RefineableSolidElement(),
        ElementWithZ2ErrorEstimator()
    {
    }

    /// Call the residuals including hanging node cases
    void fill_in_generic_contribution_to_residuals_pvd(
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
      // DIM Diagonal strain rates and DIM*(DIM-1)/2 off diagonal terms
      return DIM + DIM * (DIM - 1) / 2;
    }

    /// Get 'flux' for Z2 error recovery:   Upper triangular entries
    /// in strain tensor.
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
#ifdef PARANOID
      unsigned num_entries = DIM + ((DIM * DIM) - DIM) / 2;
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
      DenseMatrix<double> strain(DIM);
      this->get_strain(s, strain);

      // Pack into flux Vector
      unsigned icount = 0;

      // Start with diagonal terms
      for (unsigned i = 0; i < DIM; i++)
      {
        flux[icount] = strain(i, i);
        icount++;
      }

      // Off diagonals row by row
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = i + 1; j < DIM; j++)
        {
          flux[icount] = strain(i, j);
          icount++;
        }
      }
    }

    /// Number of continuously interpolated values: 0 (pure solid problem)
    unsigned ncont_interpolated_values() const
    {
      return 0;
    }

    // Return a pointer to the solid node at which pressure dof l2 is stored
    // This is only required so that the generic templating in
    // PseudoSolidNodeUpdateElements works OK
    virtual Node* solid_pressure_node_pt(const unsigned& l)
    {
      return 0;
    }

    /// Further build function, pass the pointers down to the sons
    void further_build()
    {
      RefineablePVDEquations<DIM>* cast_father_element_pt =
        dynamic_cast<RefineablePVDEquations<DIM>*>(this->father_element_pt());

      // Do whatever needs to be done in the base class
      RefineableSolidElement::further_build();

      // Set pointer to isotropic growth function
      this->Isotropic_growth_fct_pt =
        cast_father_element_pt->isotropic_growth_fct_pt();

      // Set pointer to body force function
      this->Body_force_fct_pt = cast_father_element_pt->body_force_fct_pt();

      // Set pointer to the contitutive law
      this->Constitutive_law_pt = cast_father_element_pt->constitutive_law_pt();

      // Set the timescale ratio (non-dim. density)
      this->Lambda_sq_pt = cast_father_element_pt->lambda_sq_pt();

      // Set the flag that switches inertia on/off
      this->Unsteady = cast_father_element_pt->is_inertia_enabled();

      // Evaluation of Jacobian by same method as father
      this->Evaluate_jacobian_by_fd =
        cast_father_element_pt->is_jacobian_evaluated_by_fd();
    }
  };

  //========================================================================
  /// Class for refineable QPVDElement elements
  //========================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class RefineableQPVDElement : public virtual QPVDElement<DIM, NNODE_1D>,
                                public virtual RefineablePVDEquations<DIM>,
                                public virtual RefineableSolidQElement<DIM>
  {
  public:
    /// Constructor:
    RefineableQPVDElement()
      : QPVDElement<DIM, NNODE_1D>(),
        RefineableElement(),
        RefineableSolidElement(),
        RefineablePVDEquations<DIM>(),
        RefineableSolidQElement<DIM>()
    {
    }

    /// Empty rebuild from sons, no need to reconstruct anything here
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QPVDElement<DIM, NNODE_1D>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QPVDElement<DIM, NNODE_1D>::vertex_node_pt(j);
    }

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return NNODE_1D - 1;
    }

    ///  No additional hanging node procedures are required
    /// for the solid elements.
    void further_setup_hanging_nodes() {}
  };

  //==============================================================
  /// FaceGeometry of the 2D RefineableQPVDElement elements
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<RefineableQPVDElement<2, NNODE_1D>>
    : public virtual SolidQElement<1, NNODE_1D>
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<1, NNODE_1D>() {}
  };

  //==============================================================
  /// FaceGeometry of the FaceGeometry of the 2D RefineableQPVDElement
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<FaceGeometry<RefineableQPVDElement<2, NNODE_1D>>>
    : public virtual PointElement
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : PointElement() {}
  };


  //==============================================================
  /// FaceGeometry of the 3D RefineableQPVDElement elements
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<RefineableQPVDElement<3, NNODE_1D>>
    : public virtual SolidQElement<2, NNODE_1D>
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<2, NNODE_1D>() {}
  };

  //==============================================================
  /// FaceGeometry of the FaceGeometry of the 3D RefineableQPVDElement
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<FaceGeometry<RefineableQPVDElement<3, NNODE_1D>>>
    : public virtual SolidQElement<1, NNODE_1D>
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<1, NNODE_1D>() {}
  };


  //===========================================================================
  /// Class for Refineable solid mechanics elements in near-incompressible/
  /// incompressible formulation, so a pressure is included! In this case,
  /// the pressure interpolation is discontinuous, a la Crouzeix Raviart
  //===========================================================================
  template<unsigned DIM>
  class RefineablePVDEquationsWithPressure
    : public virtual PVDEquationsWithPressure<DIM>,
      public virtual RefineableSolidElement,
      public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// Constructor:
    RefineablePVDEquationsWithPressure()
      : PVDEquationsWithPressure<DIM>(),
        RefineableElement(),
        RefineableSolidElement(),
        ElementWithZ2ErrorEstimator()
    {
    }

    /// Add element's contribution to elemental residual vector and/or
    /// Jacobian matrix
    /// flag=1: compute both
    /// flag=0: compute only residual vector
    void fill_in_generic_residual_contribution_pvd_with_pressure(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
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
      // DIM Diagonal strain rates and DIM*(DIM-1)/2 off diagonal terms
      return DIM + DIM * (DIM - 1) / 2;
    }

    // Get 'flux' for Z2 error recovery:   Upper triangular entries
    /// in strain tensor.
    //----------------------------------------------------------------
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      // Find the dimension of the problem
#ifdef PARANOID
      unsigned num_entries = DIM + ((DIM * DIM) - DIM) / 2;
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
      DenseMatrix<double> strain(DIM);
      this->get_strain(s, strain);

      // Pack into flux Vector
      unsigned icount = 0;

      // Start with diagonal terms
      for (unsigned i = 0; i < DIM; i++)
      {
        flux[icount] = strain(i, i);
        icount++;
      }

      // Off diagonals row by row
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = i + 1; j < DIM; j++)
        {
          flux[icount] = strain(i, j);
          icount++;
        }
      }
    }

    /// Number of continuously interpolated values: 0 (pure solid problem)
    unsigned ncont_interpolated_values() const
    {
      return 0;
    }

    // Return a pointer to the solid node at which pressure dof l2 is stored
    virtual Node* solid_pressure_node_pt(const unsigned& l)
    {
      return 0;
    }


    /// Pass the generic stuff down to the sons
    void further_build()
    {
      RefineablePVDEquationsWithPressure<DIM>* cast_father_element_pt =
        dynamic_cast<RefineablePVDEquationsWithPressure<DIM>*>(
          this->father_element_pt());

      // Do whatever needs to be done in the base class
      RefineableSolidElement::further_build();

      // Set pointer to isotropic growth function
      this->Isotropic_growth_fct_pt =
        cast_father_element_pt->isotropic_growth_fct_pt();

      // Set pointer to body force function
      this->Body_force_fct_pt = cast_father_element_pt->body_force_fct_pt();

      // Set pointer to the contitutive law
      this->Constitutive_law_pt = cast_father_element_pt->constitutive_law_pt();

      // Set the timescale ratio (non-dim. density)
      this->Lambda_sq_pt = cast_father_element_pt->lambda_sq_pt();

      // Set the flag that switches inertia on/off
      this->Unsteady = cast_father_element_pt->is_inertia_enabled();

      // Set the incompressibility flag
      this->Incompressible = cast_father_element_pt->is_incompressible();

      // Evaluation of Jacobian by same method as father
      this->Evaluate_jacobian_by_fd =
        cast_father_element_pt->is_jacobian_evaluated_by_fd();
    }


    /// Compute the diagonal of the displacement mass matrix for
    /// LSC preconditioner
    void get_mass_matrix_diagonal(Vector<double>& mass_diag);
  };

  //===========================================================================
  /// Class for refineable solid mechanics elements in near-incompressible/
  /// incompressible formulation, so a pressure is included! In this case,
  /// the pressure interpolation is discontinuous, a la Crouzeix Raviart,
  /// and the displacement is always quadratic.
  //===========================================================================
  template<unsigned DIM>
  class RefineableQPVDElementWithPressure
    : public virtual QPVDElementWithPressure<DIM>,
      public virtual RefineablePVDEquationsWithPressure<DIM>,
      public virtual RefineableSolidQElement<DIM>
  {
  private:
    /// Unpin all solid pressure dofs
    void unpin_elemental_solid_pressure_dofs()
    {
      unsigned n_pres = this->npres_solid();
      // loop over pressure dofs and unpin them
      for (unsigned l = 0; l < n_pres; l++)
      {
        this->internal_data_pt(this->P_solid_internal_index)->unpin(l);
      }
    }

  public:
    /// Constructor:
    RefineableQPVDElementWithPressure()
      : QPVDElementWithPressure<DIM>(),
        RefineableElement(),
        RefineableSolidElement(),
        RefineablePVDEquationsWithPressure<DIM>(),
        RefineableSolidQElement<DIM>()
    {
    }


    /// Reconstruct the pressure from the sons
    /// Must be specialized for each dimension
    inline void rebuild_from_sons(Mesh*& mesh_pt);

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QPVDElementWithPressure<DIM>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QPVDElementWithPressure<DIM>::vertex_node_pt(j);
    }

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return 2;
    } // NNODE_1D-1;}

    ///  No additional hanging node procedures are required for
    /// discontinuous solid pressures.
    void further_setup_hanging_nodes() {}

    ///  Further build: Interpolate the solid pressure values
    /// Again this must be specialised for each dimension
    inline void further_build();


    /// Number of continuously interpolated values: 0 (pure solid problem)
    unsigned ncont_interpolated_values() const
    {
      return 0;
    }
  };

  //======================================================================
  /// FaceGeometry of the 2D RefineableQPVDElementWithPressure
  //=======================================================================
  template<>
  class FaceGeometry<RefineableQPVDElementWithPressure<2>>
    : public virtual SolidQElement<1, 3>
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<1, 3>() {}
  };


  //===========================================================================
  /// FaceGeometry of the FaceGeometry of the 2D
  /// RefineableQPVDElementWithPressure
  //============================================================================
  template<>
  class FaceGeometry<FaceGeometry<RefineableQPVDElementWithPressure<2>>>
    : public virtual PointElement
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : PointElement() {}
  };

  //====================================================================
  /// 2D rebuild from sons: reconstruct the solid pressure
  //===================================================================
  template<>
  inline void RefineableQPVDElementWithPressure<2>::rebuild_from_sons(
    Mesh*& mesh_pt)
  {
    using namespace QuadTreeNames;

    // Storage for the solid pressure of each son
    double centre_solid_p[4] = {0.0, 0.0, 0.0, 0.0};

    // Loop over the sons and assign the central solid pressure
    for (unsigned ison = 0; ison < 4; ison++)
    {
      // Add the sons midnode pressure
      centre_solid_p[ison] =
        dynamic_cast<RefineableQPVDElementWithPressure<2>*>(
          quadtree_pt()->son_pt(ison)->object_pt())
          ->solid_p(0);
    }

    // Use the average for the central solid pressure
    double p_value = 0.25 * (centre_solid_p[0] + centre_solid_p[1] +
                             centre_solid_p[2] + centre_solid_p[3]);
    // Actually set the pressure
    set_solid_p(0, p_value);


    // Slope in s_0 direction
    //----------------------

    // Use average of the 2 FD approximations based on the
    // elements central pressure values
    // [Other options: Take average of the four
    // pressure derivatives]
    double slope1 = centre_solid_p[SE] - centre_solid_p[SW];
    double slope2 = centre_solid_p[NE] - centre_solid_p[NW];


    // Use the average value
    p_value = 0.5 * (slope1 + slope2);
    set_solid_p(1, p_value);

    // Slope in s_1 direction
    //----------------------

    // Use average of the 2 FD approximations based on the
    // elements central pressure values
    // [Other options: Take average of the four
    // pressure derivatives]

    slope1 = centre_solid_p[NE] - centre_solid_p[SE];
    slope2 = centre_solid_p[NW] - centre_solid_p[SW];

    // Use the average
    p_value = 0.5 * (slope1 + slope2);
    set_solid_p(2, p_value);
  }


  //==============================================================
  /// 2D further build interpolates the internal solid pressure
  /// from the father.
  //=============================================================
  template<>
  inline void RefineableQPVDElementWithPressure<2>::further_build()
  {
    RefineablePVDEquationsWithPressure<2>::further_build();

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

    // Pressure value in father element
    RefineableQPVDElementWithPressure<2>* cast_father_element_pt =
      dynamic_cast<RefineableQPVDElementWithPressure<2>*>(father_el_pt);
    double press = cast_father_element_pt->interpolated_solid_p(s_father);

    // Pressure  value gets copied straight into internal dof:
    set_solid_p(0, press);

    // The slopes get copied from father and halved
    set_solid_p(1, 0.5 * cast_father_element_pt->solid_p(1));
    set_solid_p(2, 0.5 * cast_father_element_pt->solid_p(2));
  }


  //===============================================================
  /// 3D rebuild from sons: reconstruct the pressure
  //===============================================================
  template<>
  inline void RefineableQPVDElementWithPressure<3>::rebuild_from_sons(
    Mesh*& mesh_pt)
  {
    using namespace OcTreeNames;

    // Storage for the central solid pressure of each son
    double centre_solid_p[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Loop over the sons and assign the central solid pressure
    for (unsigned ison = 0; ison < 8; ison++)
    {
      centre_solid_p[ison] = octree_pt()
                               ->son_pt(ison)
                               ->object_pt()
                               ->internal_data_pt(this->P_solid_internal_index)
                               ->value(0);
    }

    // Central pressure value:
    //-----------------------

    // Use average of the sons central pressure values
    // Other options: Take average of the four (discontinuous)
    // pressure values at the father's midpoint]
    double av_press = 0.0;

    // Loop over the sons and sum the centre pressures
    for (unsigned ison = 0; ison < 8; ison++)
    {
      av_press += centre_solid_p[ison];
    }

    // Use the average
    internal_data_pt(this->P_solid_internal_index)
      ->set_value(0, 0.125 * av_press);


    // Slope in s_0 direction
    //----------------------

    // Use average of the 4 FD approximations based on the
    // elements central pressure values
    // [Other options: Take average of the four
    // pressure derivatives]

    double slope1 = centre_solid_p[RDF] - centre_solid_p[LDF];
    double slope2 = centre_solid_p[RUF] - centre_solid_p[LUF];
    double slope3 = centre_solid_p[RDB] - centre_solid_p[LDB];
    double slope4 = centre_solid_p[RUB] - centre_solid_p[LUB];

    // Use the average
    internal_data_pt(this->P_solid_internal_index)
      ->set_value(1, 0.25 * (slope1 + slope2 + slope3 + slope4));


    // Slope in s_1 direction
    //----------------------

    // Use average of the 4 FD approximations based on the
    // elements central pressure values
    // [Other options: Take average of the four
    // pressure derivatives]
    slope1 = centre_solid_p[LUB] - centre_solid_p[LDB];
    slope2 = centre_solid_p[RUB] - centre_solid_p[RDB];
    slope3 = centre_solid_p[LUF] - centre_solid_p[LDF];
    slope4 = centre_solid_p[RUF] - centre_solid_p[RDF];

    // Use the average
    internal_data_pt(this->P_solid_internal_index)
      ->set_value(2, 0.25 * (slope1 + slope2 + slope3 + slope4));


    // Slope in s_2 direction
    //----------------------

    // Use average of the 4 FD approximations based on the
    // elements central pressure values
    // [Other options: Take average of the four
    // pressure derivatives]
    slope1 = centre_solid_p[LUF] - centre_solid_p[LUB];
    slope2 = centre_solid_p[RUF] - centre_solid_p[RUB];
    slope3 = centre_solid_p[LDF] - centre_solid_p[LDB];
    slope4 = centre_solid_p[RDF] - centre_solid_p[RDB];

    // Use the average
    internal_data_pt(this->P_solid_internal_index)
      ->set_value(3, 0.25 * (slope1 + slope2 + slope3 + slope4));
  }


  //================================================================
  ///  3D Further build: Interpolate the solid pressure values
  //=================================================================
  template<>
  inline void RefineableQPVDElementWithPressure<3>::further_build()

  {
    RefineablePVDEquationsWithPressure<3>::further_build();

    using namespace OcTreeNames;

    // What type of son am I? Ask my quadtree representation...
    int son_type = octree_pt()->son_type();

    // Pointer to my father (in element impersonation)
    RefineableQElement<3>* father_el_pt = dynamic_cast<RefineableQElement<3>*>(
      octree_pt()->father_pt()->object_pt());

    Vector<double> s_father(3);

    // Son midpoint is located at the following coordinates in father element:
    for (unsigned i = 0; i < 3; i++)
    {
      s_father[i] = 0.5 * OcTree::Direction_to_vector[son_type][i];
    }

    // Pressure value in father element
    RefineableQPVDElementWithPressure<3>* cast_father_element_pt =
      dynamic_cast<RefineableQPVDElementWithPressure<3>*>(father_el_pt);

    double press = cast_father_element_pt->interpolated_solid_p(s_father);


    // Pressure  value gets copied straight into internal dof:
    set_solid_p(0, press);

    // The slopes get copied from father and halved
    for (unsigned i = 1; i < 4; i++)
    {
      set_solid_p(i, 0.5 * cast_father_element_pt->solid_p(i));
    }
  }

  //=========================================================================
  /// FaceGeometry of the 3D RefineableQPVDElementWithPressure
  //========================================================================
  template<>
  class FaceGeometry<RefineableQPVDElementWithPressure<3>>
    : public virtual SolidQElement<2, 3>
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<2, 3>() {}
  };


  //========================================================================
  /// FaceGeometry of the FaceGeometry of the 3D
  /// RefineableQPVDElementWithPressure
  //==========================================================================
  template<>
  class FaceGeometry<FaceGeometry<RefineableQPVDElementWithPressure<3>>>
    : public virtual SolidQElement<1, 3>
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<1, 3>() {}
  };


  //===========================================================================
  /// Class for refineable solid mechanics elements in near-incompressible/
  /// incompressible formulation, so a pressure is included! These elements
  /// include a continuously interpolated pressure a la Taylor Hood/
  //===========================================================================
  template<unsigned DIM>
  class RefineableQPVDElementWithContinuousPressure
    : public virtual QPVDElementWithContinuousPressure<DIM>,
      public virtual RefineablePVDEquationsWithPressure<DIM>,
      public virtual RefineableSolidQElement<DIM>
  {
  public:
    /// Constructor:
    RefineableQPVDElementWithContinuousPressure()
      : QPVDElementWithContinuousPressure<DIM>(),
        RefineableElement(),
        RefineableSolidElement(),
        RefineablePVDEquationsWithPressure<DIM>(),
        RefineableSolidQElement<DIM>()
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

    /// Empty rebuild from sons, empty
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// OK, interpolate the solid pressures
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      // There is only one solid pressure, initialise to zero
      values.resize(1);

      // Get the interpolated value
      values[0] = this->interpolated_solid_p(s);
    }

    /// OK get the time-dependent verion
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      // There is only one solid pressure, initialise to zero
      values.resize(1);
      // The solid pressure does not depend on time!
      values[0] = this->interpolated_solid_p(s);
    }


    /// Unpin all pressure dofs
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
      unsigned n_solid_pres = this->npres_solid();
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
      return QPVDElementWithContinuousPressure<DIM>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QPVDElementWithContinuousPressure<DIM>::vertex_node_pt(j);
    }

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return 2;
    } // NNODE_1D-1;}


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
        Vector<int> index(DIM);
        // Loop over the coordinates
        for (unsigned i = 0; i < DIM; i++)
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

    /// The number of pressure nodes is 2^DIM. The number of
    /// velocity nodes is the same as the number of geometric nodes.
    unsigned ninterpolating_node(const int& value_id)
    {
#ifdef PARANOID
      RefineableElement::check_value_id(1, value_id);
#endif

      if (value_id == 0)
      {
        return static_cast<unsigned>(pow(2.0, static_cast<int>(DIM)));
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


  //=========================================================================
  /// FaceGeometry of the 2D RefineableQPVDElementWithContinuousPressure
  /// elements
  //=========================================================================
  template<>
  class FaceGeometry<RefineableQPVDElementWithContinuousPressure<2>>
    : public virtual SolidQElement<1, 3>
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<1, 3>() {}
  };

  //=========================================================================
  /// FaceGeometry of the FaceGeometry of the 2D
  /// RefineableQPVDElementWithContinuousPressure
  //=========================================================================
  template<>
  class FaceGeometry<
    FaceGeometry<RefineableQPVDElementWithContinuousPressure<2>>>
    : public virtual PointElement
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : PointElement() {}
  };

  //=========================================================================
  /// FaceGeometry of the 3D RefineableQPVDElementWithContinuousPressure
  //=========================================================================
  template<>
  class FaceGeometry<RefineableQPVDElementWithContinuousPressure<3>>
    : public virtual SolidQElement<2, 3>
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<2, 3>() {}
  };

  //=========================================================================
  /// FaceGeometry of the FaceGeometry of the 3D
  /// RefineableQPVDElementWithContinuousPressue
  //=========================================================================
  template<>
  class FaceGeometry<
    FaceGeometry<RefineableQPVDElementWithContinuousPressure<3>>>
    : public virtual SolidQElement<1, 3>
  {
  public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<1, 3>() {}
  };

} // namespace oomph

#endif
