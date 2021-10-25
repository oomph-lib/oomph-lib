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
// Header file for refineable linear elasticity elements

// Include guard to prevent multiple inclusions of this header
#ifndef OOMPH_REFINEABLE_LINEAR_ELASTICITY_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_LINEAR_ELASTICITY_ELEMENTS_HEADER

// oomph-lib headers
#include "linear_elasticity_elements.h"
#include "../generic/refineable_quad_element.h"
#include "../generic/refineable_brick_element.h"
#include "../generic/error_estimator.h"

namespace oomph
{
  //========================================================================
  /// Class for Refineable LinearElasticity equations
  //========================================================================
  template<unsigned DIM>
  class RefineableLinearElasticityEquations
    : public virtual LinearElasticityEquations<DIM>,
      public virtual RefineableElement,
      public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// Constructor
    RefineableLinearElasticityEquations()
      : LinearElasticityEquations<DIM>(),
        RefineableElement(),
        ElementWithZ2ErrorEstimator()
    {
    }


    /// Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Create enough initialised storage
      values.resize(DIM, 0.0);

      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Shape functions
      Shape psi(n_node);
      this->shape(s, psi);

      // Calculate displacements
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get the index at which the i-th velocity is stored
        unsigned u_nodal_index = this->u_index_linear_elasticity(i);
        for (unsigned l = 0; l < n_node; l++)
        {
          values[i] += this->nodal_value(t, l, u_nodal_index) * psi(l);
        }
      }
    }

    /// Get the current interpolated values (displacements).
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines),the values Vector sets its own size in here.
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      values.resize(DIM);
      this->interpolated_u_linear_elasticity(s, values);
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

    /// Number of continuously interpolated values: DIM
    unsigned ncont_interpolated_values() const
    {
      return DIM;
    }

    /// Further build function, pass the pointers down to the sons
    void further_build()
    {
      RefineableLinearElasticityEquations<DIM>* cast_father_element_pt =
        dynamic_cast<RefineableLinearElasticityEquations<DIM>*>(
          this->father_element_pt());

      // Set pointer to body force function
      this->Body_force_fct_pt = cast_father_element_pt->body_force_fct_pt();

      // Set pointer to the contitutive law
      this->Elasticity_tensor_pt =
        cast_father_element_pt->elasticity_tensor_pt();

      // Set the timescale ratio (non-dim. density)
      this->Lambda_sq_pt = cast_father_element_pt->lambda_sq_pt();

      /// Set the flag that switches inertia on/off
      this->Unsteady = cast_father_element_pt->is_inertia_enabled();
    }


  private:
    /// Overloaded helper function to take hanging nodes into account
    void fill_in_generic_contribution_to_residuals_linear_elasticity(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag);
  };


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  //========================================================================
  /// Class for refineable QLinearElasticityElement elements
  //========================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class RefineableQLinearElasticityElement
    : public virtual QLinearElasticityElement<DIM, NNODE_1D>,
      public virtual RefineableLinearElasticityEquations<DIM>,
      public virtual RefineableQElement<DIM>
  {
  public:
    /// Constructor:
    RefineableQLinearElasticityElement()
      : QLinearElasticityElement<DIM, NNODE_1D>(),
        RefineableElement(),
        RefineableLinearElasticityEquations<DIM>(),
        RefineableQElement<DIM>()
    {
    }

    /// Empty rebuild from sons, no need to reconstruct anything here
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QLinearElasticityElement<DIM, NNODE_1D>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QLinearElasticityElement<DIM, NNODE_1D>::vertex_node_pt(j);
    }

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return NNODE_1D - 1;
    }

    ///  No additional hanging node procedures are required
    void further_setup_hanging_nodes() {}
  };


  //======================================================================
  /// p-refineable version of 2D QLinearElasticityElement elements
  //======================================================================
  template<unsigned DIM>
  class PRefineableQLinearElasticityElement
    : public QLinearElasticityElement<DIM, 2>,
      public virtual RefineableLinearElasticityEquations<DIM>,
      public virtual PRefineableQElement<DIM>
  {
  public:
    /// Constructor, simply call the other constructors
    PRefineableQLinearElasticityElement()
      : RefineableElement(),
        RefineableLinearElasticityEquations<DIM>(),
        PRefineableQElement<DIM>(),
        QLinearElasticityElement<DIM, 2>()
    {
      // Set integration scheme
      // (To avoid memory leaks in pre-build and p-refine where new
      // integration schemes are created)
      this->set_integration_scheme(new GaussLobattoLegendre<DIM, 2>);
    }

    /// Destructor (to avoid memory leaks)
    ~PRefineableQLinearElasticityElement()
    {
      delete this->integral_pt();
    }


    /// Broken copy constructor
    PRefineableQLinearElasticityElement(
      const PRefineableQLinearElasticityElement<DIM>& dummy) = delete;

    /// Broken assignment operator
    /* void operator=(const PRefineableQLinearElasticityElement<DIM>&) =
     * delete;*/

    void further_build();

    /// Number of continuously interpolated values: 1
    unsigned ncont_interpolated_values() const
    {
      return 1;
    }

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QLinearElasticityElement<DIM, 2>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QLinearElasticityElement<DIM, 2>::vertex_node_pt(j);
    }

    /// Order of recovery shape functions for Z2 error estimation:
    /// - Same order as shape functions.
    // unsigned nrecovery_order()
    // {
    //  if(this->nnode_1d() < 4) {return (this->nnode_1d()-1);}
    //  else {return 3;}
    // }
    /// - Constant recovery order, since recovery order of the first element
    ///   is used for the whole mesh.
    unsigned nrecovery_order()
    {
      return 3;
    }

    void compute_energy_error(
      std::ostream& outfile,
      FiniteElement::SteadyExactSolutionFctPt exact_grad_pt,
      double& error,
      double& norm);
  };


  //==============================================================
  /// FaceGeometry of the 2D RefineableQLinearElasticityElement elements
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<RefineableQLinearElasticityElement<2, NNODE_1D>>
    : public virtual QElement<1, NNODE_1D>
  {
  public:
    // Make sure that we call the constructor of the QElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : QElement<1, NNODE_1D>() {}
  };

  //==============================================================
  /// FaceGeometry of the FaceGeometry of the 2D
  /// RefineableQLinearElasticityElement
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<
    FaceGeometry<RefineableQLinearElasticityElement<2, NNODE_1D>>>
    : public virtual PointElement
  {
  public:
    // Make sure that we call the constructor of the QElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : PointElement() {}
  };


  //==============================================================
  /// FaceGeometry of the 3D RefineableQLinearElasticityElement elements
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<RefineableQLinearElasticityElement<3, NNODE_1D>>
    : public virtual QElement<2, NNODE_1D>
  {
  public:
    // Make sure that we call the constructor of the QElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : QElement<2, NNODE_1D>() {}
  };

  //==============================================================
  /// FaceGeometry of the FaceGeometry of the 3D
  /// RefineableQLinearElasticityElement
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<
    FaceGeometry<RefineableQLinearElasticityElement<3, NNODE_1D>>>
    : public virtual QElement<1, NNODE_1D>
  {
  public:
    // Make sure that we call the constructor of the QElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : QElement<1, NNODE_1D>() {}
  };


} // namespace oomph

#endif
