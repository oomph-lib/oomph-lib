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
// Header file for refineable QPoissonElement elements

#ifndef OOMPH_REFINEABLE_POISSON_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_POISSON_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// oomph-lib headers
#include "../generic/refineable_quad_element.h"
#include "../generic/refineable_brick_element.h"
#include "../generic/hp_refineable_elements.h"
#include "../generic/error_estimator.h"
#include "poisson_elements.h"

namespace oomph
{
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Refineable version of Poisson equations
  ///
  ///
  //======================================================================
  template<unsigned DIM>
  class RefineablePoissonEquations : public virtual PoissonEquations<DIM>,
                                     public virtual RefineableElement,
                                     public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// \short Constructor, simply call other constructors
    RefineablePoissonEquations()
      : PoissonEquations<DIM>(),
        RefineableElement(),
        ElementWithZ2ErrorEstimator()
    {
    }

    /// Broken copy constructor
    RefineablePoissonEquations(const RefineablePoissonEquations<DIM>& dummy)
    {
      BrokenCopy::broken_copy("RefineablePoissonEquations");
    }

    /// Broken assignment operator
    void operator=(const RefineablePoissonEquations<DIM>&)
    {
      BrokenCopy::broken_assign("RefineablePoissonEquations");
    }

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      return DIM;
    }

    /// Get 'flux' for Z2 error recovery:  Standard flux.from Poisson equations
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      this->get_flux(s, flux);
    }

    /// Get error against and norm of exact flux
    void compute_exact_Z2_error(
      std::ostream& outfile,
      FiniteElement::SteadyExactSolutionFctPt exact_flux_pt,
      double& error,
      double& norm);

    /// \short Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Set size of Vector: u
      values.resize(1);

      // Find number of nodes
      unsigned n_node = nnode();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      values[0] = 0.0;

      // Find the index at which the poisson unknown is stored
      unsigned u_nodal_index = this->u_index_poisson();

      // Loop over the local nodes and sum up the values
      for (unsigned l = 0; l < n_node; l++)
      {
        values[0] += this->nodal_value(l, u_nodal_index) * psi[l];
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
      if (t != 0)
      {
        std::string error_message =
          "Time-dependent version of get_interpolated_values() ";
        error_message += "not implemented for this element \n";
        throw OomphLibError(
          error_message,
          "RefineablePoissonEquations::get_interpolated_values()",
          OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        // Make sure that we call this particular object's steady
        // get_interpolated_values (it could get overloaded lower down)
        RefineablePoissonEquations<DIM>::get_interpolated_values(s, values);
      }
    }


    ///  Further build: Copy source function pointer from father element
    void further_build()
    {
      this->Source_fct_pt = dynamic_cast<RefineablePoissonEquations<DIM>*>(
                              this->father_element_pt())
                              ->source_fct_pt();
    }


  private:
    /// \short Add element's contribution to elemental residual vector and/or
    /// Jacobian matrix
    /// flag=1: compute both
    /// flag=0: compute only residual vector
    void fill_in_generic_residual_contribution_poisson(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag);

    /// \short Compute derivatives of elemental residual vector with respect
    /// to nodal coordinates. Overwrites default implementation in
    /// FiniteElement base class.
    /// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
    virtual void get_dresidual_dnodal_coordinates(
      RankThreeTensor<double>& dresidual_dnodal_coordinates);
  };


  //======================================================================
  /// Refineable version of 2D QPoissonElement elements
  ///
  ///
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class RefineableQPoissonElement
    : public QPoissonElement<DIM, NNODE_1D>,
      public virtual RefineablePoissonEquations<DIM>,
      public virtual RefineableQElement<DIM>
  {
  public:
    /// \short Constructor, simply call the other constructors
    RefineableQPoissonElement()
      : RefineableElement(),
        RefineablePoissonEquations<DIM>(),
        RefineableQElement<DIM>(),
        QPoissonElement<DIM, NNODE_1D>()
    {
    }


    /// Broken copy constructor
    RefineableQPoissonElement(
      const RefineableQPoissonElement<DIM, NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("RefineableQuadPoissonElement");
    }

    /// Broken assignment operator
    void operator=(const RefineableQPoissonElement<DIM, NNODE_1D>&)
    {
      BrokenCopy::broken_assign("RefineableQuadPoissonElement");
    }

    /// Number of continuously interpolated values: 1
    unsigned ncont_interpolated_values() const
    {
      return 1;
    }

    /// \short Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QPoissonElement<DIM, NNODE_1D>::nvertex_node();
    }

    /// \short Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QPoissonElement<DIM, NNODE_1D>::vertex_node_pt(j);
    }

    /// Rebuild from sons: empty
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// \short Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return (NNODE_1D - 1);
    }

    ///  \short Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes. Empty.
    void further_setup_hanging_nodes() {}
  };


  //======================================================================
  /// p-refineable version of 2D QPoissonElement elements
  //======================================================================
  template<unsigned DIM>
  class PRefineableQPoissonElement
    : public QPoissonElement<DIM, 2>,
      public virtual RefineablePoissonEquations<DIM>,
      public virtual PRefineableQElement<DIM>
  {
  public:
    /// \short Constructor, simply call the other constructors
    PRefineableQPoissonElement()
      : RefineableElement(),
        RefineablePoissonEquations<DIM>(),
        PRefineableQElement<DIM>(),
        QPoissonElement<DIM, 2>()
    {
      // Set integration scheme
      // (To avoid memory leaks in pre-build and p-refine where new
      // integration schemes are created)
      this->set_integration_scheme(new GaussLobattoLegendre<DIM, 2>);
    }

    /// Destructor (to avoid memory leaks)
    ~PRefineableQPoissonElement()
    {
      delete this->integral_pt();
    }


    /// Broken copy constructor
    PRefineableQPoissonElement(const PRefineableQPoissonElement<DIM>& dummy)
    {
      BrokenCopy::broken_copy("PRefineableQPoissonElement");
    }

    /// Broken assignment operator
    void operator=(const PRefineableQPoissonElement<DIM>&)
    {
      BrokenCopy::broken_assign("PRefineableQPoissonElement");
    }

    void further_build();

    /// Number of continuously interpolated values: 1
    unsigned ncont_interpolated_values() const
    {
      return 1;
    }

    /// \short Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QPoissonElement<DIM, 2>::nvertex_node();
    }

    /// \short Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QPoissonElement<DIM, 2>::vertex_node_pt(j);
    }

    /// \short Order of recovery shape functions for Z2 error estimation:
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


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the RefineableQuadPoissonElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<RefineableQPoissonElement<DIM, NNODE_1D>>
    : public virtual QElement<DIM - 1, NNODE_1D>
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QElement<DIM - 1, NNODE_1D>() {}
  };

} // namespace oomph

#endif
