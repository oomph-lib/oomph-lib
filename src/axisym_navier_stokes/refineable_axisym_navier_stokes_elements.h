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
// Header file for refineable axisymmetric quad Navier Stokes elements
#ifndef OOMPH_REFINEABLE_AXISYMMETRIC_NAVIER_STOKES_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_AXISYMMETRIC_NAVIER_STOKES_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// Oomph lib includes
#include "../generic/refineable_quad_element.h"
#include "../generic/error_estimator.h"
#include "axisym_navier_stokes_elements.h"

namespace oomph
{
  //======================================================================
  /// Refineable version of the Axisymmetric Navier--Stokes equations
  //======================================================================
  class RefineableAxisymmetricNavierStokesEquations
    : public virtual AxisymmetricNavierStokesEquations,
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
    RefineableAxisymmetricNavierStokesEquations()
      : AxisymmetricNavierStokesEquations(),
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
      if (flux.size() < num_entries)
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
      RefineableAxisymmetricNavierStokesEquations* cast_father_element_pt =
        dynamic_cast<RefineableAxisymmetricNavierStokesEquations*>(
          this->father_element_pt());

      // Set the viscosity ratio pointer
      this->Viscosity_Ratio_pt = cast_father_element_pt->viscosity_ratio_pt();
      // Set the density ratio pointer
      this->Density_Ratio_pt = cast_father_element_pt->density_ratio_pt();
      // Set pointer to global Reynolds number
      this->Re_pt = cast_father_element_pt->re_pt();
      // Set pointer to global Reynolds number x Strouhal number (=Womersley)
      this->ReSt_pt = cast_father_element_pt->re_st_pt();
      // Set pointer to global Reynolds number x inverse Froude number
      this->ReInvFr_pt = cast_father_element_pt->re_invfr_pt();
      // Set pointer to the global Reynolds number x inverse Rossby number
      this->ReInvRo_pt = cast_father_element_pt->re_invro_pt();
      // Set pointer to global gravity Vector
      this->G_pt = cast_father_element_pt->g_pt();

      // Set pointer to body force function
      this->Body_force_fct_pt =
        cast_father_element_pt->axi_nst_body_force_fct_pt();

      // Set pointer to volumetric source function
      this->Source_fct_pt = cast_father_element_pt->source_fct_pt();

      // Set the ALE_is_disabled flag
      this->ALE_is_disabled = cast_father_element_pt->ALE_is_disabled;
    }

    /// \short Compute the derivatives of the i-th component of
    /// velocity at point s with respect
    /// to all data that can affect its value. In addition, return the global
    /// equation numbers corresponding to the data.
    /// Overload the non-refineable version to take account of hanging node
    /// information
    void dinterpolated_u_axi_nst_ddata(const Vector<double>& s,
                                       const unsigned& i,
                                       Vector<double>& du_ddata,
                                       Vector<unsigned>& global_eqn_number)
    {
      // Find number of nodes
      unsigned n_node = this->nnode();
      // Local shape function
      Shape psi(n_node);
      // Find values of shape function at the given local coordinate
      this->shape(s, psi);

      // Find the index at which the velocity component is stored
      const unsigned u_nodal_index = this->u_index_axi_nst(i);

      // Storage for hang info pointer
      HangInfo* hang_info_pt = 0;
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
            // Local equation number
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


    /// \short  Loop over all elements in Vector (which typically contains
    /// all the elements in a fluid mesh) and pin the nodal pressure degrees
    /// of freedom that are not being used. Function uses
    /// the member function
    /// - \c RefineableAxisymmetricNavierStokesEquations::
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
      unsigned n_element = element_pt.size();
      for (unsigned e = 0; e < n_element; e++)
      {
        dynamic_cast<RefineableAxisymmetricNavierStokesEquations*>(
          element_pt[e])
          ->pin_elemental_redundant_nodal_pressure_dofs();
      }
    }

    /// \short Unpin all pressure dofs in elements listed in vector.
    static void unpin_all_pressure_dofs(
      const Vector<GeneralisedElement*>& element_pt)
    {
      // Loop over all elements to brutally unpin all nodal pressure degrees of
      // freedom and internal pressure degrees of freedom
      unsigned n_element = element_pt.size();
      for (unsigned e = 0; e < n_element; e++)
      {
        dynamic_cast<RefineableAxisymmetricNavierStokesEquations*>(
          element_pt[e])
          ->unpin_elemental_pressure_dofs();
      }
    }

    /// \short Compute derivatives of elemental residual vector with respect to
    /// nodal coordinates. This function computes these terms analytically and
    /// overwrites the default implementation in the FiniteElement base class.
    /// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
    virtual void get_dresidual_dnodal_coordinates(
      RankThreeTensor<double>& dresidual_dnodal_coordinates);

  private:
    /// \short Add element's contribution to the elemental residual vector
    /// and/or Jacobian matrix and mass matrix
    /// flag=2: compute all
    /// flag=1: compute both residual and Jacobian
    /// flag=0: compute only residual vector
    void fill_in_generic_residual_contribution_axi_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag);

    /// \short Add element's contribution to the derivative of the
    /// elemental residual vector
    /// and/or Jacobian matrix and/or mass matrix
    /// flag=2: compute all
    /// flag=1: compute both residual and Jacobian
    /// flag=0: compute only residual vector
    void fill_in_generic_dresidual_contribution_axi_nst(
      double* const& parameter_pt,
      Vector<double>& dres_dparam,
      DenseMatrix<double>& djac_dparam,
      DenseMatrix<double>& dmass_matrix_dparam,
      unsigned flag)
    {
      throw OomphLibError("Not yet implemented\n",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    /// \short Compute the hessian tensor vector products required to
    /// perform continuation of bifurcations analytically
    void fill_in_contribution_to_hessian_vector_products(
      Vector<double> const& Y,
      DenseMatrix<double> const& C,
      DenseMatrix<double>& product)
    {
      throw OomphLibError("Not yet implemented\n",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  };

  //======================================================================
  /// Refineable version of Axisymmetric Quad Taylor Hood elements.
  /// (note that unlike the cartesian version this is not scale-able
  /// to higher dimensions!)
  //======================================================================
  class RefineableAxisymmetricQTaylorHoodElement
    : public AxisymmetricQTaylorHoodElement,
      public virtual RefineableAxisymmetricNavierStokesEquations,
      public virtual RefineableQElement<2>
  {
  private:
    /// \short Pointer to n_p-th pressure node
    Node* pressure_node_pt(const unsigned& n_p)
    {
      return this->node_pt(this->Pconv[n_p]);
    }

    /// Unpin all pressure dofs
    void unpin_elemental_pressure_dofs()
    {
      unsigned n_node = this->nnode();
      int p_index = this->p_nodal_index_axi_nst();
      // loop over nodes
      for (unsigned n = 0; n < n_node; n++)
      {
        this->node_pt(n)->unpin(p_index);
      }
    }

    ///  Unpin the proper nodal pressure dofs
    void pin_elemental_redundant_nodal_pressure_dofs()
    {
      // Loop over all nodes
      unsigned n_node = this->nnode();
      int p_index = this->p_nodal_index_axi_nst();
      // loop over all nodes and pin all  the nodal pressures
      for (unsigned n = 0; n < n_node; n++)
      {
        this->node_pt(n)->pin(p_index);
      }

      // Loop over all actual pressure nodes and unpin if they're not hanging
      unsigned n_pres = this->npres_axi_nst();
      for (unsigned l = 0; l < n_pres; l++)
      {
        Node* nod_pt = this->node_pt(this->Pconv[l]);
        if (!nod_pt->is_hanging(p_index))
        {
          nod_pt->unpin(p_index);
        }
      }
    }

  public:
    /// \short Constructor:
    RefineableAxisymmetricQTaylorHoodElement()
      : RefineableElement(),
        RefineableAxisymmetricNavierStokesEquations(),
        RefineableQElement<2>(),
        AxisymmetricQTaylorHoodElement()
    {
    }

    /// Number of values (pinned or dofs) required at node n.
    // Bumped up to 4 so we don't have to worry if a hanging mid-side node
    // gets shared by a corner node (which has extra degrees of freedom)
    unsigned required_nvalue(const unsigned& n) const
    {
      return 4;
    }

    /// Number of continuously interpolated values: 4 (3 velocities + 1
    /// pressure)
    unsigned ncont_interpolated_values() const
    {
      return 4;
    }

    /// Rebuild from sons: empty
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// \short Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return 2;
    }

    /// \short Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return AxisymmetricQTaylorHoodElement::nvertex_node();
    }

    /// \short Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return AxisymmetricQTaylorHoodElement::vertex_node_pt(j);
    }

    /// \short Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Set the velocity dimension of the element
      unsigned DIM = 3;

      // Set size of values Vector: u,w,v,p and initialise to zero
      values.resize(DIM + 1, 0.0);

      // Calculate velocities: values[0],...
      for (unsigned i = 0; i < DIM; i++)
      {
        values[i] = interpolated_u_axi_nst(s, i);
      }

      // Calculate pressure: values[DIM]
      values[DIM] = interpolated_p_axi_nst(s);
    }

    /// \short Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      unsigned DIM = 3;

      // Set size of Vector: u,w,v,p
      values.resize(DIM + 1);

      // Initialise
      for (unsigned i = 0; i < DIM + 1; i++)
      {
        values[i] = 0.0;
      }

      // Find out how many nodes there are
      unsigned n_node = nnode();

      // Shape functions
      Shape psif(n_node);
      shape(s, psif);

      // Calculate velocities: values[0],...
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get the nodal coordinate of the velocity
        unsigned u_nodal_index = u_index_axi_nst(i);
        for (unsigned l = 0; l < n_node; l++)
        {
          values[i] += nodal_value(t, l, u_nodal_index) * psif[l];
        }
      }

      // Calculate pressure: values[DIM]
      //(no history is carried in the pressure)
      values[DIM] = interpolated_p_axi_nst(s);
    }

    ///  \short Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes. The pressures are stored
    /// at the 3rd location in each node
    void further_setup_hanging_nodes()
    {
      int DIM = 3;
      this->setup_hang_for_value(DIM);
    }

    /// \short The velocities are isoparametric and so the "nodes" interpolating
    /// the velocities are the geometric nodes. The pressure "nodes" are a
    /// subset of the nodes, so when n_value==DIM, the n-th pressure
    /// node is returned.
    Node* interpolating_node_pt(const unsigned& n, const int& n_value)

    {
      int DIM = 3;
      // The only different nodes are the pressure nodes
      if (n_value == DIM)
      {
        return this->node_pt(this->Pconv[n]);
      }
      // The other variables are interpolated via the usual nodes
      else
      {
        return this->node_pt(n);
      }
    }

    /// \short The pressure nodes are the corner nodes, so when n_value==DIM,
    /// the fraction is the same as the 1d node number, 0 or 1.
    double local_one_d_fraction_of_interpolating_node(const unsigned& n1d,
                                                      const unsigned& i,
                                                      const int& n_value)
    {
      int DIM = 3;
      if (n_value == DIM)
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

    /// \short The velocity nodes are the same as the geometric nodes. The
    /// pressure nodes must be calculated by using the same methods as
    /// the geometric nodes, but by recalling that there are only two pressure
    /// nodes per edge.
    Node* get_interpolating_node_at_local_coordinate(const Vector<double>& s,
                                                     const int& n_value)
    {
      int DIM = 3;
      // If we are calculating pressure nodes
      if (n_value == static_cast<int>(DIM))
      {
        // Storage for the index of the pressure node
        unsigned total_index = 0;
        // The number of nodes along each 1d edge is 2.
        unsigned NNODE_1D = 2;
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


    /// \short The number of 1d pressure nodes is 2, the number of 1d velocity
    /// nodes is the same as the number of 1d geometric nodes.
    unsigned ninterpolating_node_1d(const int& n_value)
    {
      int DIM = 3;
      if (n_value == DIM)
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
      int DIM = 3;
      if (n_value == DIM)
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
      int DIM = 3;
      if (n_value == DIM)
      {
        return this->pshape_axi_nst(s, psi);
      }
      else
      {
        return this->shape(s, psi);
      }
    }
  };


  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Face geometry of the RefineableQuadQTaylorHoodElements
  //=======================================================================
  template<>
  class FaceGeometry<RefineableAxisymmetricQTaylorHoodElement>
    : public virtual FaceGeometry<AxisymmetricQTaylorHoodElement>
  {
  public:
    FaceGeometry() : FaceGeometry<AxisymmetricQTaylorHoodElement>() {}
  };


  //=======================================================================
  /// Face geometry of the RefineableQuadQTaylorHoodElements
  //=======================================================================
  template<>
  class FaceGeometry<FaceGeometry<RefineableAxisymmetricQTaylorHoodElement>>
    : public virtual FaceGeometry<FaceGeometry<AxisymmetricQTaylorHoodElement>>
  {
  public:
    FaceGeometry()
      : FaceGeometry<FaceGeometry<AxisymmetricQTaylorHoodElement>>()
    {
    }
  };


  //======================================================================
  /// Refineable version of Axisymmetric Quad Crouzeix Raviart elements
  /// (note that unlike the cartesian version this is not scale-able
  /// to higher dimensions!)
  //======================================================================
  class RefineableAxisymmetricQCrouzeixRaviartElement
    : public AxisymmetricQCrouzeixRaviartElement,
      public virtual RefineableAxisymmetricNavierStokesEquations,
      public virtual RefineableQElement<2>
  {
  private:
    /// Unpin all the internal pressure freedoms
    void unpin_elemental_pressure_dofs()
    {
      unsigned n_pres = this->npres_axi_nst();
      // loop over pressure dofs and unpin
      for (unsigned l = 0; l < n_pres; l++)
      {
        this->internal_data_pt(P_axi_nst_internal_index)->unpin(l);
      }
    }

  public:
    /// \short Constructor:
    RefineableAxisymmetricQCrouzeixRaviartElement()
      : RefineableElement(),
        RefineableAxisymmetricNavierStokesEquations(),
        RefineableQElement<2>(),
        AxisymmetricQCrouzeixRaviartElement()
    {
    }

    /// Number of continuously interpolated values: 3 (velocities)
    unsigned ncont_interpolated_values() const
    {
      return 3;
    }

    /// Rebuild from sons: Reconstruct pressure from the (merged) sons
    void rebuild_from_sons(Mesh*& mesh_pt)
    {
      using namespace QuadTreeNames;

      // Central pressure value:
      //-----------------------

      // Use average of the sons central pressure values
      // Other options: Take average of the four (discontinuous)
      // pressure values at the father's midpoint]

      double av_press = 0.0;

      // Loop over the sons
      for (unsigned ison = 0; ison < 4; ison++)
      {
        // Add the sons midnode pressure
        av_press += quadtree_pt()
                      ->son_pt(ison)
                      ->object_pt()
                      ->internal_data_pt(P_axi_nst_internal_index)
                      ->value(0);
      }

      // Use the average
      internal_data_pt(P_axi_nst_internal_index)->set_value(0, 0.25 * av_press);


      // Slope in s_0 direction
      //----------------------

      // Use average of the 2 FD approximations based on the
      // elements central pressure values
      // [Other options: Take average of the four
      // pressure derivatives]

      double slope1 = quadtree_pt()
                        ->son_pt(SE)
                        ->object_pt()
                        ->internal_data_pt(P_axi_nst_internal_index)
                        ->value(0) -
                      quadtree_pt()
                        ->son_pt(SW)
                        ->object_pt()
                        ->internal_data_pt(P_axi_nst_internal_index)
                        ->value(0);

      double slope2 = quadtree_pt()
                        ->son_pt(NE)
                        ->object_pt()
                        ->internal_data_pt(P_axi_nst_internal_index)
                        ->value(0) -
                      quadtree_pt()
                        ->son_pt(NW)
                        ->object_pt()
                        ->internal_data_pt(P_axi_nst_internal_index)
                        ->value(0);


      // Use the average
      internal_data_pt(P_axi_nst_internal_index)
        ->set_value(1, 0.5 * (slope1 + slope2));


      // Slope in s_1 direction
      //----------------------

      // Use average of the 2 FD approximations based on the
      // elements central pressure values
      // [Other options: Take average of the four
      // pressure derivatives]

      slope1 = quadtree_pt()
                 ->son_pt(NE)
                 ->object_pt()
                 ->internal_data_pt(P_axi_nst_internal_index)
                 ->value(0) -
               quadtree_pt()
                 ->son_pt(SE)
                 ->object_pt()
                 ->internal_data_pt(P_axi_nst_internal_index)
                 ->value(0);

      slope2 = quadtree_pt()
                 ->son_pt(NW)
                 ->object_pt()
                 ->internal_data_pt(P_axi_nst_internal_index)
                 ->value(0) -
               quadtree_pt()
                 ->son_pt(SW)
                 ->object_pt()
                 ->internal_data_pt(P_axi_nst_internal_index)
                 ->value(0);


      // Use the average
      internal_data_pt(P_axi_nst_internal_index)
        ->set_value(2, 0.5 * (slope1 + slope2));
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
      return AxisymmetricQCrouzeixRaviartElement::nvertex_node();
    }

    /// \short Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return AxisymmetricQCrouzeixRaviartElement::vertex_node_pt(j);
    }

    /// \short Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      unsigned DIM = 3;

      // Set size of Vector: u,w,v and initialise to zero
      values.resize(DIM, 0.0);

      // Calculate velocities: values[0],...
      for (unsigned i = 0; i < DIM; i++)
      {
        values[i] = interpolated_u_axi_nst(s, i);
      }
    }

    /// \short Get all function values [u,v..,p] at previous timestep t
    /// (t=0: present; t>0: previous timestep).
    /// \n
    /// Note: Given the generality of the interface (this function is
    /// usually called from black-box documentation or interpolation routines),
    /// the values Vector sets its own size in here.
    /// \n
    /// Note: No pressure history is kept, so pressure is always
    /// the current value.
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      unsigned DIM = 3;

      // Set size of Vector: u,w,v
      values.resize(DIM);

      // Initialise
      for (unsigned i = 0; i < DIM; i++)
      {
        values[i] = 0.0;
      }

      // Find out how many nodes there are
      unsigned n_node = nnode();

      // Shape functions
      Shape psif(n_node);
      shape(s, psif);

      // Calculate velocities: values[0],...
      for (unsigned i = 0; i < DIM; i++)
      {
        // Get the local index at which the i-th velocity is stored
        unsigned u_local_index = u_index_axi_nst(i);
        for (unsigned l = 0; l < n_node; l++)
        {
          values[i] += nodal_value(t, l, u_local_index) * psif[l];
        }
      }
    }

    ///  \short Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes. Empty
    void further_setup_hanging_nodes() {}

    /// Further build for Crouzeix_Raviart interpolates the internal
    /// pressure dofs from father element: Make sure pressure values and
    /// dp/ds agree between fathers and sons at the midpoints of the son
    /// elements.
    void further_build()
    {
      // Call the generic further build
      RefineableAxisymmetricNavierStokesEquations::further_build();

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
      RefineableAxisymmetricQCrouzeixRaviartElement* cast_father_el_pt =
        dynamic_cast<RefineableAxisymmetricQCrouzeixRaviartElement*>(
          father_el_pt);
      double press = cast_father_el_pt->interpolated_p_axi_nst(s_father);

      // Pressure value gets copied straight into internal dof:
      internal_data_pt(P_axi_nst_internal_index)->set_value(0, press);


      // The slopes get copied from father
      internal_data_pt(P_axi_nst_internal_index)
        ->set_value(
          1,
          0.5 * cast_father_el_pt->internal_data_pt(P_axi_nst_internal_index)
                  ->value(1));

      internal_data_pt(P_axi_nst_internal_index)
        ->set_value(
          2,
          0.5 * cast_father_el_pt->internal_data_pt(P_axi_nst_internal_index)
                  ->value(2));
    }
  };


  //=======================================================================
  /// Face geometry of the RefineableQuadQCrouzeixRaviartElements
  //=======================================================================
  template<>
  class FaceGeometry<RefineableAxisymmetricQCrouzeixRaviartElement>
    : public virtual FaceGeometry<AxisymmetricQCrouzeixRaviartElement>
  {
  public:
    FaceGeometry() : FaceGeometry<AxisymmetricQCrouzeixRaviartElement>() {}
  };


  //=======================================================================
  /// Face geometry of the RefineableQuadQCrouzeixRaviartElements
  //=======================================================================
  template<>
  class FaceGeometry<
    FaceGeometry<RefineableAxisymmetricQCrouzeixRaviartElement>>
    : public virtual FaceGeometry<
        FaceGeometry<AxisymmetricQCrouzeixRaviartElement>>
  {
  public:
    FaceGeometry()
      : FaceGeometry<FaceGeometry<AxisymmetricQCrouzeixRaviartElement>>()
    {
    }
  };


} // namespace oomph

#endif
