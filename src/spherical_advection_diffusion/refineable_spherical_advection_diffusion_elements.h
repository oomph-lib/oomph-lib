// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision$
// LIC//
// LIC// $LastChangedDate$
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
// Header file for elements that solve the advection diffusion equation
// and that can be refined.

#ifndef OOMPH_REFINEABLE_SPHERICAL_ADVECTION_DIFFUSION_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_SPHERICAL_ADVECTION_DIFFUSION_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// oomph-lib headers
#include "../generic/refineable_quad_element.h"
#include "../generic/refineable_brick_element.h"
#include "../generic/error_estimator.h"
#include "spherical_advection_diffusion_elements.h"

namespace oomph
{
  //======================================================================
  /// \short A version of the Advection Diffusion in spherical
  /// coordinates equations that can be
  /// used with non-uniform mesh refinement. In essence, the class overloads
  /// the fill_in_generic_residual_contribution_spherical_adv_diff()
  /// function so that contributions
  /// from hanging nodes (or alternatively in-compatible function values)
  /// are taken into account.
  //======================================================================
  class RefineableSphericalAdvectionDiffusionEquations :
    public virtual SphericalAdvectionDiffusionEquations,
    public virtual RefineableElement,
    public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// \short Empty Constructor
    RefineableSphericalAdvectionDiffusionEquations() :
      SphericalAdvectionDiffusionEquations(),
      RefineableElement(),
      ElementWithZ2ErrorEstimator()
    {
    }

    /// Broken copy constructor
    RefineableSphericalAdvectionDiffusionEquations(
      const RefineableSphericalAdvectionDiffusionEquations &dummy)
    {
      BrokenCopy::broken_copy("RefineableSphericalAdvectionDiffusionEquations");
    }

    /// Broken assignment operator
    void operator=(const RefineableSphericalAdvectionDiffusionEquations &)
    {
      BrokenCopy::broken_assign(
        "RefineableSphericalAdvectionDiffusionEquations");
    }

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      return 2;
    }

    /// \short Get 'flux' for Z2 error recovery:
    /// Standard flux.from AdvectionDiffusion equations
    void get_Z2_flux(const Vector<double> &s, Vector<double> &flux)
    {
      this->get_flux(s, flux);
    }

    /// \short Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const Vector<double> &s,
                                 Vector<double> &values)
    {
      // Set size of Vector: u
      values.resize(1);

      // Find number of nodes
      unsigned n_node = nnode();

      // Find the index at which the unknown is stored
      unsigned u_nodal_index = this->u_index_spherical_adv_diff();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      values[0] = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        values[0] += this->nodal_value(l, u_nodal_index) * psi[l];
      }
    }

    /// \short Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const unsigned &t,
                                 const Vector<double> &s,
                                 Vector<double> &values)
    {
      // Set size of Vector: u
      values.resize(1);

      // Find number of nodes
      const unsigned n_node = nnode();

      // Find the index at which the unknown is stored
      const unsigned u_nodal_index = this->u_index_spherical_adv_diff();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      values[0] = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        values[0] += this->nodal_value(t, l, u_nodal_index) * psi[l];
      }
      // }

      if (t != 0)
      {
        std::string error_message =
          "Time-dependent version of get_interpolated_values() ";
        error_message += "not implemented for this element \n";
        throw OomphLibError(error_message,
                            "RefineableSphericalAdvectionDiffusionEquations::"
                            "get_interpolated_values()",
                            OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        // Make sure that we call the appropriate steady version
        //(the entire function might be overloaded lower down)
        RefineableSphericalAdvectionDiffusionEquations::get_interpolated_values(
          s, values);
      }
    }

    /// Fill in the geometric Jacobian, which in this case is r*r*sin(theta)
    double geometric_jacobian(const Vector<double> &x)
    {
      return x[0] * x[0] * sin(x[1]);
    }

    ///  Further build: Copy source function pointer from father element
    void further_build()
    {
      RefineableSphericalAdvectionDiffusionEquations *cast_father_element_pt =
        dynamic_cast<RefineableSphericalAdvectionDiffusionEquations *>(
          this->father_element_pt());

      // Set the values of the pointers from the father
      this->Source_fct_pt = cast_father_element_pt->source_fct_pt();
      this->Wind_fct_pt = cast_father_element_pt->wind_fct_pt();
      this->Pe_pt = cast_father_element_pt->pe_pt();
      this->PeSt_pt = cast_father_element_pt->pe_st_pt();

      // Set the ALE status
      // this->ALE_is_disabled = cast_father_element_pt->ALE_is_disabled;
    }

    /// \short Compute the derivatives of the i-th component of
    /// velocity at point s with respect
    /// to all data that can affect its value. In addition, return the global
    /// equation numbers corresponding to the data.
    /// Overload the non-refineable version to take account of hanging node
    /// information
    void dinterpolated_u_adv_diff_ddata(const Vector<double> &s,
                                        Vector<double> &du_ddata,
                                        Vector<unsigned> &global_eqn_number)
    {
      // Find number of nodes
      unsigned n_node = this->nnode();
      // Local shape function
      Shape psi(n_node);
      // Find values of shape function at the given local coordinate
      this->shape(s, psi);

      // Find the index at which the velocity component is stored
      const unsigned u_nodal_index = this->u_index_spherical_adv_diff();

      // Storage for hang info pointer
      HangInfo *hang_info_pt = 0;
      // Storage for global equation
      int global_eqn = 0;

      // Find the number of dofs associated with interpolated u
      unsigned n_u_dof = 0;
      for (unsigned l = 0; l < n_node; l++)
      {
        unsigned n_master = 1;

        // Local bool (is the node hanging)
        bool is_node_hanging = this->node_pt(l)->is_hanging();

        // If the node is hanging, get the number of master nodes
        if (is_node_hanging)
        {
          hang_info_pt = this->node_pt(l)->hanging_pt();
          n_master = hang_info_pt->nmaster();
        }
        // Otherwise there is just one master node, the node itself
        else
        {
          n_master = 1;
        }

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Get the equation number
          if (is_node_hanging)
          {
            // Get the equation number from the master node
            global_eqn =
              hang_info_pt->master_node_pt(m)->eqn_number(u_nodal_index);
          }
          else
          {
            // Global equation number
            global_eqn = this->node_pt(l)->eqn_number(u_nodal_index);
          }

          // If it's positive add to the count
          if (global_eqn >= 0)
          {
            ++n_u_dof;
          }
        }
      }

      // Now resize the storage schemes
      du_ddata.resize(n_u_dof, 0.0);
      global_eqn_number.resize(n_u_dof, 0);

      // Loop over th nodes again and set the derivatives
      unsigned count = 0;
      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        unsigned n_master = 1;
        double hang_weight = 1.0;

        // Local bool (is the node hanging)
        bool is_node_hanging = this->node_pt(l)->is_hanging();

        // If the node is hanging, get the number of master nodes
        if (is_node_hanging)
        {
          hang_info_pt = this->node_pt(l)->hanging_pt();
          n_master = hang_info_pt->nmaster();
        }
        // Otherwise there is just one master node, the node itself
        else
        {
          n_master = 1;
        }

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // If the node is hanging get weight from master node
          if (is_node_hanging)
          {
            // Get the hang weight from the master node
            hang_weight = hang_info_pt->master_weight(m);
          }
          else
          {
            // Node contributes with full weight
            hang_weight = 1.0;
          }

          // Get the equation number
          if (is_node_hanging)
          {
            // Get the equation number from the master node
            global_eqn =
              hang_info_pt->master_node_pt(m)->eqn_number(u_nodal_index);
          }
          else
          {
            // Global equation number
            global_eqn = this->node_pt(l)->eqn_number(u_nodal_index);
          }

          if (global_eqn >= 0)
          {
            // Set the global equation number
            global_eqn_number[count] = global_eqn;
            // Set the derivative with respect to the unknown
            du_ddata[count] = psi[l] * hang_weight;
            // Increase the counter
            ++count;
          }
        }
      }
    }

  protected:
    /// \short Add the element's contribution to the elemental residual vector
    /// and/or Jacobian matrix
    /// flag=1: compute both
    /// flag=0: compute only residual vector
    void fill_in_generic_residual_contribution_spherical_adv_diff(
      Vector<double> &residuals,
      DenseMatrix<double> &jacobian,
      DenseMatrix<double> &mass_matrix,
      unsigned flag);
  };

  //======================================================================
  /// \short Refineable version of QSphericalAdvectionDiffusionElement.
  /// Inherit from the standard QSphericalAdvectionDiffusionElement and the
  /// appropriate refineable geometric element and the refineable equations.
  //======================================================================
  template<unsigned NNODE_1D>
  class RefineableQSphericalAdvectionDiffusionElement :
    public QSphericalAdvectionDiffusionElement<NNODE_1D>,
    public virtual RefineableSphericalAdvectionDiffusionEquations,
    public virtual RefineableQElement<2>
  {
  public:
    /// \short Empty Constructor:
    RefineableQSphericalAdvectionDiffusionElement() :
      RefineableElement(),
      RefineableSphericalAdvectionDiffusionEquations(),
      RefineableQElement<2>(),
      QSphericalAdvectionDiffusionElement<NNODE_1D>()
    {
    }

    /// Broken copy constructor
    RefineableQSphericalAdvectionDiffusionElement(
      const RefineableQSphericalAdvectionDiffusionElement<NNODE_1D> &dummy)
    {
      BrokenCopy::broken_copy("RefineableQSphericalAdvectionDiffusionElement");
    }

    /// Broken assignment operator
    void operator=(
      const RefineableQSphericalAdvectionDiffusionElement<NNODE_1D> &)
    {
      BrokenCopy::broken_assign(
        "RefineableQSphericalAdvectionDiffusionElement");
    }

    /// Number of continuously interpolated values: 1
    unsigned ncont_interpolated_values() const
    {
      return 1;
    }

    /// \short Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QSphericalAdvectionDiffusionElement<NNODE_1D>::nvertex_node();
    }

    /// \short Pointer to the j-th vertex node in the element
    Node *vertex_node_pt(const unsigned &j) const
    {
      return QSphericalAdvectionDiffusionElement<NNODE_1D>::vertex_node_pt(j);
    }

    /// Rebuild from sons: empty
    void rebuild_from_sons(Mesh *&mesh_pt) {}

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

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Face geometry for the RefineableQSphericalAdvectionDiffusionElement
  /// elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<RefineableQSphericalAdvectionDiffusionElement<NNODE_1D>> :
    public virtual QElement<1, NNODE_1D>
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QElement<1, NNODE_1D>() {}
  };

} // namespace oomph

#endif
