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
// Header file for Navier Stokes elements with singularity
#ifndef SINGULAR_AXISYM_NAVIER_STOKES_ELEMENTS_HEADER
#define SINGULAR_AXISYM_NAVIER_STOKES_ELEMENTS_HEADER

// oomph-lib includes
#include "navier_stokes.h"
#include "axisym_navier_stokes.h"

// local includes
#include "singular_navier_stokes_solution_elements.h"

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

namespace oomph
{
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Templated wrapper to add handling of singularities to
  /// the underlying Navier-Stokes element (specified via the template
  /// argument). Slightly inefficient because the integration loop
  /// is repeated here. The only alternative is to add the additional
  /// functionality into the Navier-Stokes elements. Not pretty either!
  /// NOTE Element ia assumed to be of Taylor Hood type with pressures
  /// stored at nodes.
  ///
  /// Dirichlet BCs are applied in hijacking manner and must be imposed
  /// from each element that shares a boundary node that is subject to
  /// a Dirichlet condition.
  //=======================================================================
  template<class BASIC_AXISYM_NAVIER_STOKES_ELEMENT>
  class SingularAxisymNavierStokesElement
    : public virtual BASIC_AXISYM_NAVIER_STOKES_ELEMENT
  {
  private:
    double Error;

    bool IsAugmented;
    /// Vector of pointers to SingularNavierStokesSolutionElement objects
    Vector<SingularNavierStokesSolutionElement<
      SingularAxisymNavierStokesElement<BASIC_AXISYM_NAVIER_STOKES_ELEMENT>>*>
      C_equation_elements_pt;

    /// Vector indicating which velocity component of
    /// which node is subject to Dirichlet BC
    /// [size = number of nodes; initialised to false]
    Vector<std::vector<bool>> Node_is_subject_to_velocity_dirichlet_bcs;

    /// Imposed values of velocity component at nodes
    /// that are subject to Dirichlet BC
    /// [size = number of nodes; initialised to zero]
    Vector<Vector<double>> Imposed_velocity_values_at_node;

    /// Vector indicating which pressure dof is subject to Dirichlet BC
    /// [size = number of pressure dofs; initialised to false]
    std::vector<bool> Pressure_dof_is_subject_to_dirichlet_bc;

    /// Imposed value at pressure dofs
    /// that are subject to Dirichlet BC
    /// [size = number of pressure dof; initialised to zero]
    Vector<double> Imposed_value_at_pressure_dof;

  public:
    /// Constructor
    SingularAxisymNavierStokesElement()
      : BASIC_AXISYM_NAVIER_STOKES_ELEMENT(), Error(0.0), IsAugmented(false)
    {
      // Find the number of nodes in the element
      const unsigned n_node = this->nnode();

      // Find the number of velocity components
      const unsigned cached_n_u_nst = this->n_u_nst();

      // Initialise the vector indicating which node is subject to velocity
      // Dirichlet BCs. The size of the vector is equal to the number of nodes
      // in the element. Each component of the vector is a vector of booleans
      // indicating if the velocity components at the corresponding node are
      // subject to Dirichlet BC. By default, no node is subject to Dirichlet
      // BC, so the vector is full of false.
      Node_is_subject_to_velocity_dirichlet_bcs.resize(n_node);
      for (unsigned j = 0; j < n_node; j++)
      {
        Node_is_subject_to_velocity_dirichlet_bcs[j].resize(cached_n_u_nst);
        for (unsigned d = 0; d < cached_n_u_nst; d++)
        {
          Node_is_subject_to_velocity_dirichlet_bcs[j][d] = false;
        }
      }

      // Initialise the vector of imposed velocity values on the nodes
      // subject to Dirichlet BC. The size of the vector is equal to the
      // number of nodes in the element. Each component of the vector is
      // a vector of length the dimension of the problem. This vector contains
      // the imposed values of the velocity vector at the corresponding node.
      // If a node is not subject to Dirichlet BC, its imposed values are zero.
      // By default, no node is subject to Dirichlet BC so the vector is full
      // of zeros
      Imposed_velocity_values_at_node.resize(n_node);
      for (unsigned j = 0; j < n_node; j++)
      {
        Imposed_velocity_values_at_node[j].resize(cached_n_u_nst);
        for (unsigned d = 0; d < cached_n_u_nst; d++)
        {
          Imposed_velocity_values_at_node[j][d] = 0.0;
        }
      }

      // Find the number of pressure dofs
      unsigned n_pres = this->npres_nst();

      // Initialise the vector indicating which pressure dof is subject to
      // Dirichlet BCs. The size of the vector is equal to the number of
      // pressure dofs in the element. Each component of the vector is a boolean
      // indicating if the corresponding pressure unknown is subject to
      // Dirichlet BC. By default, no pressure dof is subject to Dirichlet BC,
      // so the vector is full of false.
      Pressure_dof_is_subject_to_dirichlet_bc.resize(n_pres);
      for (unsigned l = 0; l < n_pres; l++)
      {
        Pressure_dof_is_subject_to_dirichlet_bc[l] = false;
      }

      // Initialise the vector of imposed values on the pressure dofs
      // subject to Dirichlet BC. The size of the vector is equal to the
      // number of pressure dofs in the element. Each component of
      // the vector contains the imposed value at the corresponding pressure
      // dof. If a pressure dof is not subject to Dirichlet BC, its imposed
      // value is zero. By default, no pressure dof is subject to Dirichlet BC
      // so the vector is full of zeros
      Imposed_value_at_pressure_dof.resize(n_node);
      for (unsigned l = 0; l < n_pres; l++)
      {
        Imposed_value_at_pressure_dof[l] = 0.0;
      }
    } // End of constructor

    // Set error value for post-processing
    void set_error(const double& error)
    {
      Error = error;
    }

    void get_error(double& error)
    {
      error = Error;
    }

    //==============================================================
    /// Get strain-rate tensor: \f$ e_{ij} \f$  where
    /// \f$ i,j = r,z,\theta \f$ (in that order)
    //==============================================================
    void strain_rate(const Vector<double>& s,
                     DenseMatrix<double>& strainrate) const
    {
#ifdef PARANOID
      if ((strainrate.ncol() != 3) || (strainrate.nrow() != 3))
      {
        std::ostringstream error_message;
        error_message << "The strain rate has incorrect dimensions "
                      << strainrate.ncol() << " x " << strainrate.nrow()
                      << " Not 3" << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Find out how many nodes there are in the element
      unsigned n_node = this->nnode();

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node, 2);

      // Call the derivatives of the shape functions
      this->dshape_eulerian(s, psi, dpsidx);

      // Radius
      double interpolated_r = 0.0;

      // Velocity components and their derivatives
      double ur = 0.0;
      double durdr = 0.0;
      double durdz = 0.0;
      double uz = 0.0;
      double duzdr = 0.0;
      double duzdz = 0.0;
      double uphi = 0.0;
      double duphidr = 0.0;
      double duphidz = 0.0;

      // Loop over nodes to assemble velocities and their derivatives
      // w.r.t. to r and z (x_0 and x_1)
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_r += this->nodal_position(l, 0) * psi[l];

        ur += u_reconstructed(l, 0) * psi[l];
        uz += u_reconstructed(l, 1) * psi[l];
        uphi += u_reconstructed(l, 2) * psi[l];

        durdr += u_reconstructed(l, 0) * dpsidx(l, 0);
        durdz += u_reconstructed(l, 0) * dpsidx(l, 1);

        duzdr += u_reconstructed(l, 1) * dpsidx(l, 0);
        duzdz += u_reconstructed(l, 1) * dpsidx(l, 1);

        duphidr += u_reconstructed(l, 2) * dpsidx(l, 0);
        duphidz += u_reconstructed(l, 2) * dpsidx(l, 1);
      }


      // Assign strain rates without negative powers of the radius
      // and zero those with:
      strainrate(0, 0) = durdr;
      strainrate(0, 1) = 0.5 * (durdz + duzdr);
      strainrate(1, 0) = strainrate(0, 1);
      strainrate(0, 2) = 0.0;
      strainrate(2, 0) = strainrate(0, 2);
      strainrate(1, 1) = duzdz;
      strainrate(1, 2) = 0.5 * duphidz;
      strainrate(2, 1) = strainrate(1, 2);
      strainrate(2, 2) = 0.0;


      // Overwrite the strain rates with negative powers of the radius
      // unless we're at the origin
      if (std::fabs(interpolated_r) > 1.0e-16)
      {
        double inverse_radius = 1.0 / interpolated_r;
        strainrate(0, 2) = 0.5 * (duphidr - inverse_radius * uphi);
        strainrate(2, 0) = strainrate(0, 2);
        strainrate(2, 2) = inverse_radius * ur;
      }
    }

    /// Number of values (pinned or dofs) required at node n. Can
    /// be overwritten for hanging node version
    inline virtual unsigned required_nvalue(const unsigned& n) const
    {
      return BASIC_AXISYM_NAVIER_STOKES_ELEMENT::required_nvalue(n);
      //      // Basic element
      //      // dim * [u] + 1 * [p]
      //
      //      // This element
      //      // dim * [u] + 1 * [p_tilde] + dim * [u_tilde]
      //
      //      return TTaylorHoodElement<2>::required_nvalue(n) + this->dim();
    }

    // Make element augmented
    void augment()
    {
      // Add extra data
      const unsigned n_node = this->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        unsigned n_value = this->node_pt(n)->nvalue();
        if (this->is_pressure_node(n))
        {
          if (n_value != this->n_u_nst() * 2 + 1)
          {
            n_value += this->n_u_nst();
            this->node_pt(n)->resize(n_value);
          }
        }
        else
        {
          if (n_value != this->n_u_nst() * 2)
          {
            n_value += this->n_u_nst();
            this->node_pt(n)->resize(n_value);
          }
        }
      }

      IsAugmented = true;
    }

    void setup_new_data()
    {
      // Copy over data
      const unsigned n_node = this->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        Vector<double> local_u_bar = u_bar(this->node_pt(n)->position());
        const unsigned n_dim = 3;
        for (unsigned i = 0; i < n_dim; i++)
        {
          this->node_pt(n)->set_value(u_index_axi_nst(n, i),
                                      u_reconstructed(n, i) - local_u_bar[i]);
        }
      }
    }


    // Make element un-augmented
    void unaugment()
    {
      IsAugmented = false;
    }

    // Access function for IsAugmented, return (constant bool)
    const bool is_augmented()
    {
      return IsAugmented;
    }

    /// Return the index at which the i-th unknown velocity component
    /// is stored. The default value, i, is appropriate for single-physics
    /// problems.
    /// In derived multi-physics elements, this function should be
    /// overloaded to reflect the chosen storage scheme.
    virtual inline unsigned u_index_axi_nst(const unsigned& n,
                                            const unsigned& i) const
    {
      unsigned index = -1;
      if (IsAugmented)
      {
        if (this->is_pressure_node(n))
        {
          // dim * [u] + 1 * [p] + i * [u_tilde]
          index = this->n_u_nst() + 1 + i;
        }
        else
        {
          // dim * [u] + i * [u_tilde]
          index = this->n_u_nst() + i;
        }
      }
      else
      {
        index = BASIC_AXISYM_NAVIER_STOKES_ELEMENT::u_index_axi_nst(n, i);
      }
      return index;
    }

    virtual inline unsigned u_reconstructed_index(const unsigned& n,
                                                  const unsigned& i) const
    {
      // unsigned index = -1;
      // if (this->is_pressure_node(n))
      //{
      //   // dim * [u] + 1 * [p] + i * [u_tilde]
      //   index = this->n_u_nst() + 1 + i;
      // }
      // else
      //{
      //   // dim * [u] + i * [u_tilde]
      //   index = this->n_u_nst() + i;
      // }
      // return index;

      // i * [u]
      return i;
    }

    double u_reconstructed(const unsigned& n, const unsigned& i) const
    {
      // i * [u]
      return this->nodal_value(n, this->u_reconstructed_index(n, i));
    }

    virtual inline unsigned axi_momentum_index_nst(const unsigned& n,
                                                   const unsigned& i) const
    {
      return i;
    }

    virtual inline int total_velocity_eqn_index(const unsigned& n,
                                                const unsigned& i) const
    {
      // return this->u_reconstructed_index(n, i);
      // return this->u_index_axi_nst(n, i);

      unsigned index = -1;
      if (this->is_pressure_node(n))
      {
        // dim * [u] + 1 * [p] + i * [u_tilde]
        index = this->n_u_nst() + 1 + i;
      }
      else
      {
        // dim * [u] + i * [u_tilde]
        index = this->n_u_nst() + i;
      }
      return index;
    }

    virtual inline int total_velocity_local_eqn(const unsigned& n,
                                                const unsigned& i) const
    {
      // Note: u_index_axi_nst references u_tilde_index

      // return nodal_local_eqn(n, u_index_axi_nst(n, i));
      return this->nodal_local_eqn(n, total_velocity_eqn_index(n, i));
    }

    virtual inline int u_reconstructed_local_unknown(const unsigned& n,
                                                     const unsigned& i) const
    {
      // Note: u_index_axi_nst references u_tilde_index

      // return nodal_local_eqn(n, u_index_axi_nst(n, i));
      return this->nodal_local_eqn(n, this->u_reconstructed_index(n, i));
    }

  public:
    void pin()
    {
      const unsigned n_node = this->nnode();
      for (unsigned i = 0; i < n_node; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          this->node_pt(i)->pin(this->axi_momentum_index_nst(i, j));
        }
        if (this->IsAugmented)
        {
          for (unsigned j = 0; j < 2; j++)
          {
            this->node_pt(i)->pin(this->total_velocity_eqn_index(i, j));
          }
        }
      }
      const unsigned& n_pres = 3;
      for (unsigned i = 0; i < n_pres; i++)
      {
        this->node_pt(i)->pin(this->p_index_axi_nst());
      }
    }

    void pin_total_velocity_eqn(const unsigned& i_node, const unsigned& i_dim)
    {
      if (this->IsAugmented)
      {
        this->node_pt(i_node)->pin(
          this->total_velocity_eqn_index(i_node, i_dim));
      }
    }

    void set_total_velocity_value(const unsigned& i_node,
                                  const unsigned& i_dim,
                                  const double& value)
    {
      if (this->IsAugmented)
      {
        this->node_pt(i_node)->set_value(
          this->u_reconstructed_index(i_node, i_dim), value);
      }
    }

    void set_velocity_value(const unsigned& i_node,
                            const unsigned& i_dim,
                            const double& value)
    {
      this->node_pt(i_node)->set_value(this->u_index_axi_nst(i_node, i_dim),
                                       value);
    }

    void pin_momentum_eqn(const unsigned& i_node, const unsigned& i_dim)
    {
      this->node_pt(i_node)->pin(this->axi_momentum_index_nst(i_node, i_dim));
    }

    /// Impose Dirichlet BC on the d-th component of the velocity
    /// (including the singular contribution) at the j-th node
    void impose_velocity_dirichlet_bc_on_node(const unsigned& j,
                                              const unsigned& d)
    {
      Node_is_subject_to_velocity_dirichlet_bcs[j][d] = true;
    }

    /// Undo Dirichlet BC on the d-th velocity component (including the
    /// singular contribution) of the jth node
    void undo_velocity_dirichlet_bc_on_node(const unsigned& j,
                                            const unsigned& d)
    {
      Node_is_subject_to_velocity_dirichlet_bcs[j][d] = false;
    }

    /// Specify Dirichlet boundary value for the d-th velocity component
    /// (including the singular contribution) at the j-th local node
    void set_velocity_dirichlet_value_on_node(const unsigned& j,
                                              const unsigned& d,
                                              const double& value)
    {
      Imposed_velocity_values_at_node[j][d] = value;
    }


    /// Impose Dirichlet BC at the j-th pressure dof
    void impose_dirichlet_bc_on_pressure_dof(const unsigned& j)
    {
      Pressure_dof_is_subject_to_dirichlet_bc[j] = true;
    }

    /// Undo Dirichlet BC at the j-th pressure dof
    void undo_dirichlet_bc_on_pressure_dof(const unsigned& j)
    {
      Pressure_dof_is_subject_to_dirichlet_bc[j] = false;
    }

    /// Specify Dirichlet boundary value for the j-th pressure dof
    void set_dirichlet_value_on_pressure_dof(const unsigned& j,
                                             const double& value)
    {
      Imposed_value_at_pressure_dof[j] = value;
    }

    /// Access function to vector of pointers to
    /// SingularNavierStokesSolutionElements
    Vector<SingularNavierStokesSolutionElement<
      SingularAxisymNavierStokesElement<BASIC_AXISYM_NAVIER_STOKES_ELEMENT>>*>
    c_equation_elements_pt()
    {
      return C_equation_elements_pt;
    }

    /// Add pointer to associated SingularNavierStokesSolutionElement that
    /// determines the value of the amplitude of the singular functions (and
    /// gives access to the singular functions). The unknown amplitude
    /// becomes external Data for this element so assign_eqn_numbers() must
    /// be called after this function has been called.
    void add_c_equation_element_pt(
      SingularNavierStokesSolutionElement<
        SingularAxisymNavierStokesElement<BASIC_AXISYM_NAVIER_STOKES_ELEMENT>>*
        c_pt)
    {
      // Add the element
      C_equation_elements_pt.push_back(c_pt);

      // Add the additional unknown of this object as external data in the
      // Navier-Stokes element
      this->add_external_data(c_pt->internal_data_pt(0));
    }

    /// Assuming Taylor-Hood elements
    virtual double p_axi_nst(const unsigned& n) const
    {
      return this->nodal_value(n, this->p_nodal_index_axi_nst());
    }

    /// Derivative of pressure in direction indicated by pointer to unsigned
    double dpdx_fe_only(Vector<double> s, const unsigned* direction_pt)
    {
      // Find the number of pressure dofs in the wrapped Navier-Stokes
      // element pointed by
      // the SingularNavierStokesSolutionElement class
      unsigned n_pres = this->npres_nst();

      // Find the dimension of the problem
      unsigned cached_dim = this->dim();

      // Set up memory for the pressure shape functions and their
      // derivatives
      Shape psip(n_pres), testp(n_pres);
      DShape dpsipdx(n_pres, cached_dim), dtestpdx(n_pres, cached_dim);

      // Compute the pressure shape functions and their derivatives
      // at the local coordinate S_in_wrapped_navier_stokes_element
      // (Test fcts not really needed but nobody's got around to writing
      // a fct that only picks out the basis fcts.
      this->dpshape_and_dptest_eulerian_nst(s, psip, dpsipdx, testp, dtestpdx);

      // Initialise the derivative used for the residual
      double interpolated_dpdx_fe_only = 0.0;

      // Compute the derivative used for the residual. The direction of the
      // derivative is given by *Direction_pt
      for (unsigned j = 0; j < n_pres; j++)
      {
        interpolated_dpdx_fe_only += p_axi_nst(j) * dpsipdx(j, *direction_pt);
      }

      return interpolated_dpdx_fe_only;
    }


    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      BASIC_AXISYM_NAVIER_STOKES_ELEMENT::fill_in_contribution_to_residuals(
        residuals);

      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      if (this->IsAugmented)
      {
        fill_in_generic_residual_contribution_additional_terms(
          residuals, GeneralisedElement::Dummy_matrix, 0);
      }
    }

    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1 and dummy mass
      // matrix
      if (!this->IsAugmented)
      {
        BASIC_AXISYM_NAVIER_STOKES_ELEMENT::fill_in_contribution_to_jacobian(
          residuals, jacobian);

        this->fill_in_jacobian_from_external_by_fd(residuals, jacobian);
      }

      else
      {
        SolidFiniteElement::fill_in_contribution_to_jacobian(residuals,
                                                             jacobian);
        // if (this->IsAugmented)
        //{
        //   fill_in_generic_residual_contribution_additional_terms(
        //     residuals, jacobian, 1);
        // }
      }
    }

    /// Overload the output function
    /// x, y, [z,] u, v, [w], p, u_fe, v_fe, [w_fe], p_fe,
    // u_sing, v_sing, [w_sing], p_sing
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      const bool IsStressOutputIncluded = false;
      const bool IsSingularOutputInclude = false;

      // Find the dimension of the problem
      unsigned cached_dim = this->dim();

      // Vector of local coordinates
      Vector<double> s(cached_dim);
      Vector<double> xi(cached_dim);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);
        this->interpolated_xi(s, xi);

        Vector<double> x(cached_dim);
        for (unsigned i = 0; i < cached_dim; i++)
        {
          outfile << this->interpolated_x(s, i) << " ";
          x[i] = this->interpolated_x(s, i);
        }

        // Total velocity
        Vector<double> velocity;

        if (IsAugmented)
        {
          velocity = interpolated_u_nst_new(s);
        }
        else
        {
          velocity = this->interpolated_u_nst_fe_only(s);
        }
        for (unsigned i = 0; i < this->n_u_nst(); i++)
        {
          outfile << velocity[i] << " ";
        }

        // "Total" pressure
        outfile << this->interpolated_p_nst_fe_only(s) << " ";

        // Lagrange coordinates
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << this->interpolated_x(s, i) - xi[i] << " ";
        }

        // Error
        outfile << Error << " ";

        // Size
        outfile << this->size() << " ";

        if (IsStressOutputIncluded)
        {
          // Total stress
          DenseMatrix<double> stress_tensor = interpolated_stress_tensor(s);
          for (unsigned i = 0; i < cached_dim; i++)
          {
            for (unsigned j = 0; j < cached_dim; j++)
            {
              outfile << stress_tensor(i, j) << " ";
            }
          }
        }

        if (IsSingularOutputInclude)
        {
          // Finite element Velocity
          Vector<double> velocity_fe_only(this->n_u_nst(), 0.0);
          if (this->IsAugmented)
          {
            velocity_fe_only = interpolated_u_nst_fe_only(s);
          }

          for (unsigned i = 0; i < this->n_u_nst(); i++)
          {
            outfile << velocity_fe_only[i] << " ";
          }

          // Singular Velocity
          Vector<double> velocity_bar(this->n_u_nst(), 0.0);
          if (this->IsAugmented)
          {
            velocity_bar = u_bar(x);
          }
          for (unsigned i = 0; i < this->n_u_nst(); i++)
          {
            outfile << velocity_bar[i] << " ";
          }

          if (IsStressOutputIncluded)
          {
            // Singular stress
            DenseMatrix<double> stress_tensor_bar(cached_dim, cached_dim, 0.0);
            if (this->IsAugmented)
            {
              stress_tensor_bar = interpolated_stress_tensor_bar(s);
            }
            for (unsigned i = 0; i < cached_dim; i++)
            {
              for (unsigned j = 0; j < cached_dim; j++)
              {
                outfile << stress_tensor_bar(i, j) << " ";
              }
            }
          }
        } // End of augmented

        // Error
        double local_error = 0.0;
        this->get_error(local_error);
        outfile << local_error << " ";

        // Size
        outfile << this->size() << " ";

        outfile << std::endl;
      }

      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);
    }


    /// Overloaded compute error function; uses FE+singular parts
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       const bool& include_pressure,
                       double& error,
                       double& norm)
    {
      unsigned cached_dim = this->dim();

      error = 0.0;
      norm = 0.0;

      // Vector of local coordinates
      Vector<double> s(cached_dim);

      // Vector for coordintes
      Vector<double> x(cached_dim);

      // Set the value of n_intpt
      unsigned n_intpt = this->integral_pt()->nweight();


      outfile << "ZONE" << std::endl;

      // Exact solution Vector (u,v,[w],p)
      Vector<double> exact_soln(this->n_u_nst() + 1);
      Vector<double> computed_soln(this->n_u_nst() + 1);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s
        for (unsigned i = 0; i < cached_dim; i++)
        {
          s[i] = this->integral_pt()->knot(ipt, i);
        }

        // Get the integral weight
        double w = this->integral_pt()->weight(ipt);

        // Get jacobian of mapping
        double J = this->J_eulerian(s);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get x position as Vector
        this->interpolated_x(s, x);

        // Get exact solution at this point
        (*exact_soln_pt)(x, exact_soln);
        Vector<double> u_comp = interpolated_u_nst_sum(s);
        for (unsigned i = 0; i < this->n_u_nst(); i++)
        {
          computed_soln[i] = u_comp[i];
        }
        unsigned hi_limit = this->n_u_nst();
        if (include_pressure)
        {
          computed_soln[this->n_u_nst()] = my_interpolated_p_nst(s);
          hi_limit = this->n_u_nst() + 1;
        }


        // Velocity error
        for (unsigned i = 0; i < hi_limit; i++)
        {
          norm += exact_soln[i] * exact_soln[i] * W;
          error += (exact_soln[i] - computed_soln[i]) *
                   (exact_soln[i] - computed_soln[i]) * W;
        }

        // Output x,y,...,u_exact,...]
        for (unsigned i = 0; i < cached_dim; i++)
        {
          outfile << x[i] << " ";
        }

        // Output x,y,[z],u_error,v_error,[w_error], [p_error]
        for (unsigned i = 0; i < hi_limit; i++)
        {
          outfile << exact_soln[i] - computed_soln[i] << " ";
        }
        outfile << std::endl;
      }
    }

    /// Get 'flux' for Z2 error recovery:   Upper triangular entries
    /// in strain rate tensor.
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
#ifdef PARANOID
      unsigned num_entries = 6;
      if (flux.size() < num_entries)
      {
        std::ostringstream error_message;
        error_message << "The flux vector has the wrong number of entries, "
                      << flux.size() << ", whereas it should be at least "
                      << num_entries << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get strain rate matrix
      DenseMatrix<double> strainrate(3);
      this->strain_rate(s, strainrate);

      // Pack into flux Vector
      unsigned icount = 0;

      // Start with diagonal terms
      for (unsigned i = 0; i < 3; i++)
      {
        flux[icount] = strainrate(i, i);
        icount++;
      }

      // Off diagonals row by row
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = i + 1; j < 3; j++)
        {
          flux[icount] = strainrate(i, j);
          icount++;
        }
      }
    }


    /// Overloaded version of the interpolated velocity solution including
    /// the singular contributions
    inline Vector<double> interpolated_u_nst_sum(const Vector<double>& s) const
    {
      // Find the dimension of the problem
      unsigned cached_dim = this->dim();

      // Find number of nodes
      const unsigned n_node = this->nnode();

      // Initialise value of u
      Vector<double> interpolated_u(this->n_u_nst(), 0.0);

      // Calculate the global coordinate associated with s
      Vector<double> x(cached_dim, 0.0);

      // Local shape function
      Shape psif(n_node);

      // Find values of shape function
      this->shape(s, psif);

      // Loop over the local nodes and sum
      for (unsigned j = 0; j < n_node; j++)
      {
        // Loop over the spatial directions
        for (unsigned d = 0; d < cached_dim; d++)
        {
          x[d] += this->nodal_position(j, d) * psif(j);
        }
        for (unsigned d = 0; d < this->n_u_nst(); d++)
        {
          interpolated_u[d] +=
            this->nodal_value(j, u_index_axi_nst(j, d)) * psif[j];
        }
      }

      // Add the contribution of the singularities (all of them, summed)
      Vector<double> u_bar_local = u_bar(x);
      for (unsigned d = 0; d < this->n_u_nst(); d++)
      {
        interpolated_u[d] += u_bar_local[d];
      }
      return (interpolated_u);
    }

    /// Overloaded version of the interpolated velocity solution including
    /// the singular contributions
    Vector<double> interpolated_u_nst_new(const Vector<double>& s)
    {
      // Find number of nodes
      const unsigned n_node = this->nnode();

      // Initialise value of u
      Vector<double> interpolated_u(this->n_u_nst(), 0.0);

      // Local shape function
      Shape psif(n_node);

      // Find values of shape function
      this->shape(s, psif);

      // Loop over the spatial directions
      for (unsigned d = 0; d < this->n_u_nst(); d++)
      {
        // Loop over the local nodes and sum
        for (unsigned j = 0; j < n_node; j++)
        {
          interpolated_u[d] += u_reconstructed(j, d) * psif[j];
        }
      }
      return interpolated_u;
    }

    void interpolated_u_axi_nst(const Vector<double>& s, Vector<double>& veloc)
    {
      if (IsAugmented)
      {
        veloc = interpolated_u_nst_new(s);
      }
      else
      {
        veloc = this->interpolated_u_nst_fe_only(s);
      }
    }

    /// Return FE interpolated derivatives of velocity component u[i]
    /// w.r.t spatial global coordinate direction x[j] at local coordinate s
    double interpolated_dudx_axi_nst(const Vector<double>& s,
                                     const unsigned& i,
                                     const unsigned& j)
    {
      // Find number of nodes
      const unsigned n_node = this->nnode();

      // Find the dimension of the problem
      const unsigned cached_dim = this->dim();

      // Local shape function
      Shape psif(n_node);
      DShape dpsifdx(n_node, cached_dim);

      // Find values of shape function
      this->dshape_eulerian(s, psif, dpsifdx);

      // Get the velocity gradient
      double interpolated_dudx(0.0);
      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        double u_value(0.0);
        if (this->is_augmented())
        {
          u_value = u_reconstructed(l, i);
        }
        else
        {
          u_value = this->nodal_value(l, u_index_axi_nst(l, i));
        }
        interpolated_dudx += u_value * dpsifdx(l, j);
      }

      return interpolated_dudx;
    }

    /// Return FE interpolated derivatives of velocity component u[i]
    /// w.r.t time at local coordinate s
    double interpolated_dudt_axi_nst(const Vector<double>& s,
                                     const unsigned& i) const
    {
      // Determine number of nodes
      const unsigned n_node = this->nnode();

      // Allocate storage for local shape function
      Shape psif(n_node);

      // Find values of shape function
      this->shape(s, psif);

      // Initialise value of dudt
      double interpolated_dudt = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_dudt += du_dt_axi_nst(l, i) * psif[l];
      }

      return (interpolated_dudt);
    }

    /// i-th component of du/dt at local node n.
    /// Uses suitably interpolated value for hanging nodes.
    double du_dt_axi_nst(const unsigned& n, const unsigned& i) const
    {
      // Get the data's timestepper
      TimeStepper* time_stepper_pt = this->node_pt(n)->time_stepper_pt();

      // Initialise dudt
      double dudt = 0.0;
      // Loop over the timesteps, if there is a non Steady timestepper
      if (!time_stepper_pt->is_steady())
      {
        // Number of timsteps (past & present)
        const unsigned n_time = time_stepper_pt->ntstorage();

        // Add the contributions to the time derivative
        for (unsigned t = 0; t < n_time; t++)
        {
          dudt += time_stepper_pt->weight(1, t) *
                  this->nodal_value(t, n, this->u_reconstructed_index(n, i));
        }
      }

      return dudt;
    }

    /// Return FE interpolated derivatives of velocity component u[i]
    /// w.r.t spatial local coordinate direction s[j] at local coordinate s
    double interpolated_duds_axi_nst(const Vector<double>& s,
                                     const unsigned& i,
                                     const unsigned& j) const
    {
      // Determine number of nodes
      const unsigned n_node = this->nnode();

      // Allocate storage for local shape function and its derivatives
      // with respect to space
      Shape psif(n_node);
      DShape dpsifds(n_node, 2);

      // Find values of shape function (ignore jacobian)
      (void)this->dshape_local(s, psif, dpsifds);

      // Initialise value of duds
      double interpolated_duds = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_duds +=
          this->nodal_value(l, this->u_reconstructed_index(l, i)) *
          dpsifds(l, j);
      }

      return (interpolated_duds);
    }

    /// Version of interpolated pressure including the singular
    /// contributions
    inline double my_interpolated_p_nst(const Vector<double>& s) const
    {
      // Initialise pressure value with fe part
      double interpolated_p = interpolated_p_nst_fe_only(s);

      // Calculate the global coordinate associated with s
      unsigned cached_dim = this->dim();
      Vector<double> x(cached_dim, 0.0);

      // Find number of nodes
      const unsigned n_node = this->nnode();

      // Local shape function
      Shape psif(n_node);

      // Find values of shape function
      this->shape(s, psif);

      // Loop over the local nodes and sum
      for (unsigned j = 0; j < n_node; j++)
      {
        for (unsigned d = 0; d < cached_dim; d++)
        {
          x[d] += this->raw_nodal_position(j, d) * psif(j);
        }
      }

      // Add the singularities contribution (all of them, summed)
      interpolated_p += p_bar(x);

      return interpolated_p;
    }

    DenseMatrix<double> interpolated_stress_tensor(const Vector<double>& s)
    {
      // Find number of nodes
      const unsigned n_node = this->nnode();

      // Find the dimension of the problem
      const unsigned cached_dim = this->dim();

      // Local shape function
      Shape psif(n_node);
      DShape dpsifdx(n_node, cached_dim);

      // Find values of shape function
      this->dshape_eulerian(s, psif, dpsifdx);

      // Get the velocity gradient
      DenseMatrix<double> interpolated_dudx(cached_dim, cached_dim, 0.0);
      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the spatial directions
        for (unsigned i = 0; i < cached_dim; i++)
        {
          const double u_value = this->nodal_value(l, u_index_axi_nst(l, i));
          for (unsigned j = 0; j < cached_dim; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      // Construct the stress tensor
      // Loop over the spatial directions
      DenseMatrix<double> stress_tensor(cached_dim, cached_dim, 0.0);
      const double visc_ratio = this->viscosity_ratio();

      for (unsigned i = 0; i < cached_dim; i++)
      {
        for (unsigned j = 0; j < cached_dim; j++)
        {
          stress_tensor(i, j) -=
            visc_ratio * (interpolated_dudx(i, j) +
                          this->Gamma[i] * interpolated_dudx(j, i));
        }
      }

      return stress_tensor;
    }

    DenseMatrix<double> interpolated_stress_tensor_bar(const Vector<double>& s)
    {
      // Find the dimension of the problem
      const unsigned cached_dim = this->dim();
      DenseMatrix<double> stress_tensor(cached_dim, cached_dim, 0.0);
      if (this->IsAugmented)
      {
        Vector<double> x(cached_dim);
        for (unsigned i = 0; i < cached_dim; i++)
        {
          x[i] = this->interpolated_x(s, i);
        }
        Vector<Vector<double>> grad_u_bar_local = this->grad_u_bar(x);

        // Construct the stress tensor
        // Loop over the spatial directions
        const double visc_ratio = this->viscosity_ratio();
        for (unsigned i = 0; i < cached_dim; i++)
        {
          for (unsigned j = 0; j < cached_dim; j++)
          {
            stress_tensor(i, j) -=
              visc_ratio * (grad_u_bar_local[i][j] +
                            this->Gamma[i] * grad_u_bar_local[j][i]);
          }
        }
      }

      return stress_tensor;
    }


  private:
    /// Evaluate i-th u_bar (i-th velocity singular function
    /// incl. the amplitudes) function at Eulerian position x
    Vector<double> u_bar(const unsigned& i, const Vector<double>& x) const
    {
      return C_equation_elements_pt[i]->u_bar(x);
    }

    /// Evaluate sum of all velocity singular fcts
    /// (incl. the amplitude) at Eulerian position x
    Vector<double> u_bar(const Vector<double>& x) const
    {
      // Find the number of singularities
      unsigned n_sing = C_equation_elements_pt.size();

      // Find the dimension of the problem
      Vector<double> sum(this->n_u_nst(), 0.0);
      for (unsigned s = 0; s < n_sing; s++)
      {
        Vector<double> u_bar_local = C_equation_elements_pt[s]->u_bar(x);
        for (unsigned i = 0; i < this->n_u_nst(); i++)
        {
          sum[i] += u_bar_local[i];
        }
      }
      return sum;
    }

    /// Evaluate i-th grad_u_bar (i-th gradient of velocity singular
    /// fct incl. the amplitude) function at Eulerian position x;
    /// grad[i][j] = du_i/dx_j
    Vector<Vector<double>> grad_u_bar(const unsigned& i,
                                      const Vector<double>& x) const
    {
      return C_equation_elements_pt[i]->grad_u_bar(x);
    }

    /// Evaluate gradient of sum of all velocity singular fcts
    /// (incl. the amplitudes) at Eulerian position x: grad[i][j] =
    /// du_i/dx_j
    Vector<Vector<double>> grad_u_bar(const Vector<double>& x) const
    {
      // Find the number of singularities
      unsigned n_sing = C_equation_elements_pt.size();

      // Find the dimension of the problem
      unsigned cached_dim = this->dim();
      Vector<Vector<double>> sum(cached_dim);
      for (unsigned i = 0; i < cached_dim; i++)
      {
        sum[i].resize(cached_dim, 0.0);
      }
      for (unsigned s = 0; s < n_sing; s++)
      {
        Vector<Vector<double>> grad_u_bar_local =
          C_equation_elements_pt[s]->grad_u_bar(x);
        for (unsigned i = 0; i < cached_dim; i++)
        {
          for (unsigned j = 0; j < cached_dim; j++)
          {
            sum[i][j] += grad_u_bar_local[i][j];
          }
        }
      }
      return sum;
    }


    /// Evaluate i-th "raw" velocity singular function at Eulerian position
    /// x
    Vector<double> velocity_singular_function(const unsigned& i,
                                              const Vector<double>& x) const
    {
      return C_equation_elements_pt[i]->velocity_singular_function(x);
    }


    /// Evaluate gradient of i-th "raw" velocity singular function at
    /// Eulerian position x
    Vector<Vector<double>> grad_velocity_singular_function(
      const unsigned& i, const Vector<double>& x) const
    {
      return C_equation_elements_pt[i]->grad_velocity_singular_function(x);
    }

    /// Evaluate i-th pressure singular fct (without the
    /// amplitude) at Eulerian position x
    double pressure_singular_function(const unsigned& i,
                                      const Vector<double>& x) const
    {
      return C_equation_elements_pt[i]->pressure_singular_function(x);
    }


    /// Evaluate i-th p_bar (i-th pressure singular fct (incl. the
    /// amplitude) at Eulerian position x
    double p_bar(const unsigned& i, const Vector<double>& x) const
    {
      return C_equation_elements_pt[i]->p_bar(x);
    }

    /// Evaluate sum of all pressure singular fcts
    /// (incl. the amplitudes) at Eulerian position x
    double p_bar(const Vector<double>& x) const
    {
      // Find the number of singularities
      unsigned n_sing = C_equation_elements_pt.size();

      double sum = 0.0;
      for (unsigned i = 0; i < n_sing; i++)
      {
        sum += C_equation_elements_pt[i]->p_bar(x);
      }
      return sum;
    }

    /// Return FE representation of velocity solution WITHOUT singular
    /// contributions at local coordinate s
    Vector<double> interpolated_u_nst_fe_only(const Vector<double>& s) const
    {
      // Find number of nodes
      const unsigned n_node = this->nnode();

      // Initialise value of u
      Vector<double> interpolated_u(this->n_u_nst(), 0.0);

      // Local shape function
      Shape psif(n_node);

      // Find values of shape function
      this->shape(s, psif);

      // Loop over the velocity components
      for (unsigned d = 0; d < this->n_u_nst(); d++)
      {
        // Add the contribution of u_FE
        // Loop over the local nodes and sum
        for (unsigned j = 0; j < n_node; j++)
        {
          interpolated_u[d] +=
            this->nodal_value(j, u_index_axi_nst(j, d)) * psif[j];
        }
      }

      return interpolated_u;
    }


  public:
    /// Return FE representation of pressure solution WITHOUT singular
    /// contributions at local coordinate s
    double interpolated_p_nst_fe_only(const Vector<double>& s) const
    {
      // Find number of pressure degrees of freedom
      unsigned n_pres = this->npres_nst();

      // Set up memory for pressure shape and test functions
      Shape psip(n_pres), testp(n_pres);

      // Compute the values of the pressure shape and test functions at the
      // local coordinate s
      this->pshape_axi_nst(s, psip, testp);

      // Initialise pressure value
      double interpolated_p = 0.0;

      // Add the contribution of p_FE
      // Loop over the pressure dof and sum
      for (unsigned l = 0; l < n_pres; l++)
      {
        interpolated_p += p_axi_nst(l) * psip[l];
      }

      return interpolated_p;
    }

  private:
    void fill_in_generic_residual_contribution_additional_terms(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      // Find the dimension of the problem
      unsigned cached_dim = this->dim();

      // Find the number of singularities
      unsigned n_sing = C_equation_elements_pt.size();

      // Find the local equation number of the additional unknowns
      Vector<int> local_equation_number_C(n_sing);
      for (unsigned i = 0; i < n_sing; i++)
      {
        local_equation_number_C[i] = this->external_local_eqn(i, 0);
      }

      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Find out how many pressure dofs there are
      unsigned n_pres = this->npres_nst();

      // integer to store the local equations
      int local_eqn = 0, local_unknown = 0;

      // Get Physical Variables from Element
      // Reynolds number must be multiplied by the density ratio
      double scaled_re = this->re() * this->density_ratio();

      // Set up memory for the velocity shape and test functions
      Shape psif(n_node), testf(n_node);
      DShape dpsifdx(n_node, cached_dim), dtestfdx(n_node, cached_dim);

      // Set up memory for pressure shape and test functions
      Shape psip(n_pres), testp(n_pres);

      // Number of integration points
      unsigned n_intpt = this->integral_pt()->nweight();

      // Set the vector to hold local coordinates
      Vector<double> s(cached_dim);

      // Cachec viscosity ratio
      double visc_ratio = this->viscosity_ratio();

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s
        for (unsigned d = 0; d < cached_dim; d++)
        {
          s[d] = this->integral_pt()->knot(ipt, d);
        }

        // Get the integral weight
        double w = this->integral_pt()->weight(ipt);

        // Call the derivatives of the velocity shape and test functions
        double J = this->dshape_and_dtest_eulerian_at_knot_axi_nst(
          ipt, psif, dpsifdx, testf, dtestfdx);

        // Call the pressure shape and test functions
        this->pshape_axi_nst(s, psip, testp);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Initialise the global coordinate and velocity
        Vector<double> interpolated_x(cached_dim, 0.0);
        Vector<double> interpolated_u_tilde(this->n_u_nst(), 0.0);
        Vector<double> interpolated_u(this->n_u_nst(), 0.0);
        DenseMatrix<double> interpolated_dudx(cached_dim, cached_dim, 0.0);

        // Calculate the global coordinate associated with s
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < cached_dim; i++)
          {
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
          }
          // Loop over directions
          for (unsigned i = 0; i < this->n_u_nst(); i++)
          {
            double u_value = this->nodal_value(l, u_index_axi_nst(l, i));
            interpolated_u_tilde[i] += u_value * psif[l];
            interpolated_u[i] += u_reconstructed(l, i) * psif[l];
          }
          for (unsigned i = 0; i < cached_dim; i++)
          {
            double u_value = this->nodal_value(l, u_index_axi_nst(l, i));
            for (unsigned j = 0; j < cached_dim; j++)
            {
              interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
            }
          }
        }

        // Get sum of singular functions
        Vector<double> u_bar_local = this->u_bar(interpolated_x);
        Vector<Vector<double>> grad_u_bar_local =
          this->grad_u_bar(interpolated_x);
        double p_bar_local = this->p_bar(interpolated_x);

        // Singular functions
        Vector<Vector<double>> u_hat_local(n_sing);
        for (unsigned s = 0; s < n_sing; s++)
        {
          u_hat_local[s] = this->velocity_singular_function(s, interpolated_x);
        }
        Vector<Vector<Vector<double>>> grad_u_hat_local(n_sing);
        for (unsigned s = 0; s < n_sing; s++)
        {
          grad_u_hat_local[s] =
            this->grad_velocity_singular_function(s, interpolated_x);
        }
        Vector<double> p_hat_local(n_sing);
        for (unsigned s = 0; s < n_sing; s++)
        {
          p_hat_local[s] = this->pressure_singular_function(s, interpolated_x);
        }

        const double r = interpolated_x[0];

        // MOMENTUM EQUATIONS
        //-------------------

        // Loop over the velocity test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // --------------------------------
          // FIRST (RADIAL) MOMENTUM EQUATION
          // --------------------------------
          local_eqn = this->axi_momentum_local_eqn(l, 0);

          // If it is not pinned
          if (local_eqn >= 0)
          {
            // If it is not a Dirichlet BC
            if (not(Node_is_subject_to_velocity_dirichlet_bcs[l][0]))
            {
              // Stress contribution
              // -------------------
              // Pressure
              residuals[local_eqn] +=
                p_bar_local * (testf[l] + r * dtestfdx(l, 0)) * W;

              // Shear stress
              residuals[local_eqn] -= visc_ratio * r * (1.0 + this->Gamma[0]) *
                                      grad_u_bar_local[0][0] * dtestfdx(l, 0) *
                                      W;

              residuals[local_eqn] -=
                visc_ratio * r *
                (grad_u_bar_local[0][1] +
                 this->Gamma[0] * grad_u_bar_local[1][0]) *
                dtestfdx(l, 1) * W;

              residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) *
                                      u_bar_local[0] * testf[l] * W / r;

              // Add singular convective terms
              // -----------------------------
              // Radial
              residuals[local_eqn] -= scaled_re * r * u_bar_local[0] *
                                      grad_u_bar_local[0][0] * testf[l] * W;

              residuals[local_eqn] -=
                scaled_re *
                (r * u_bar_local[0] * interpolated_dudx(0, 0) +
                 r * interpolated_u[0] * grad_u_bar_local[0][0]) *
                testf[l] * W;

              // Axial
              residuals[local_eqn] -= scaled_re * r * u_bar_local[1] *
                                      grad_u_bar_local[0][1] * testf[l] * W;

              residuals[local_eqn] -=
                scaled_re *
                (r * u_bar_local[1] * interpolated_dudx(0, 1) +
                 r * interpolated_u[1] * grad_u_bar_local[0][1]) *
                testf[l] * W;

              // Azimuthal
              residuals[local_eqn] +=
                scaled_re * u_bar_local[2] * u_bar_local[2] * testf[l] * W;

              residuals[local_eqn] += scaled_re *
                                      (u_bar_local[2] * interpolated_u[2] +
                                       interpolated_u[2] * u_bar_local[2]) *
                                      testf[l] * W;

              // Jacobian
              if (flag)
              {
                for (unsigned ss = 0; ss < n_sing; ss++)
                {
                  const int local_unknown = local_equation_number_C[ss];
                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      p_hat_local[ss] * (testf[l] + r * dtestfdx(l, 0)) * W;

                    jacobian(local_eqn, local_unknown) -=
                      visc_ratio * r * (1.0 + this->Gamma[0]) *
                      grad_u_hat_local[ss][0][0] * dtestfdx(l, 0) * W;

                    jacobian(local_eqn, local_unknown) -=
                      visc_ratio * r *
                      (grad_u_hat_local[ss][0][1] +
                       this->Gamma[0] * grad_u_hat_local[ss][1][0]) *
                      dtestfdx(l, 1) * W;

                    jacobian(local_eqn, local_unknown) -=
                      visc_ratio * (1.0 + this->Gamma[0]) * u_hat_local[ss][0] *
                      testf[l] * W / r;


                    // Radial
                    jacobian(local_eqn, local_unknown) -=
                      scaled_re * r *
                      (u_hat_local[ss][0] * grad_u_bar_local[0][0] +
                       u_bar_local[0] * grad_u_hat_local[ss][0][0]) *
                      testf[l] * W;

                    jacobian(local_eqn, local_unknown) -=
                      scaled_re *
                      (r * u_hat_local[ss][0] * interpolated_dudx(0, 0) +
                       r * interpolated_u[0] * grad_u_hat_local[ss][0][0]) *
                      testf[l] * W;

                    // Axial
                    jacobian(local_eqn, local_unknown) -=
                      scaled_re *
                      (r * u_hat_local[ss][1] * grad_u_bar_local[0][1] +
                       r * u_bar_local[1] * grad_u_hat_local[ss][0][1]) *
                      testf[l] * W;

                    jacobian(local_eqn, local_unknown) -=
                      scaled_re *
                      (r * u_hat_local[ss][1] * interpolated_dudx(0, 1) +
                       r * interpolated_u[1] * grad_u_hat_local[ss][0][1]) *
                      testf[l] * W;

                    // Azimuthal
                    jacobian(local_eqn, local_unknown) +=
                      scaled_re *
                      (u_hat_local[ss][2] * u_bar_local[2] +
                       u_bar_local[2] * u_hat_local[ss][2]) *
                      testf[l] * W;

                    jacobian(local_eqn, local_unknown) +=
                      scaled_re *
                      (u_hat_local[ss][2] * interpolated_u[2] +
                       interpolated_u[2] * u_hat_local[ss][2]) *
                      testf[l] * W;
                  }
                }

                // Loop over the velocity shape functions again
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  local_unknown = this->u_axi_nst_local_unknown(l2, 0);
                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) -=
                      scaled_re *
                      (r * u_bar_local[0] * dpsifdx(0, 0) +
                       r * psif[0] * grad_u_bar_local[0][0]) *
                      testf[l] * W;
                  }


                  local_unknown = this->u_axi_nst_local_unknown(l2, 1);
                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) -=
                      scaled_re *
                      (r * u_bar_local[1] * dpsifdx(0, 1) +
                       r * psif[1] * grad_u_bar_local[0][1]) *
                      testf[l] * W;
                  }

                  local_unknown = this->u_axi_nst_local_unknown(l2, 2);
                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      scaled_re *
                      (u_bar_local[2] * psif[2] + psif[2] * u_bar_local[2]) *
                      testf[l] * W;
                  }
                }


                // NOTE: Pressure contribution not include as for my case it
                // is zero

              } // End of Jacobian
            } // End of check of the Dirichlet status
          } // End of check of the pin status

          // --------------------------------
          // SECOND (AXIAL) MOMENTUM EQUATION
          // --------------------------------
          local_eqn = this->axi_momentum_local_eqn(l, 1);

          // If it is not pinned
          if (local_eqn >= 0)
          {
            // If it is not a Dirichlet BC
            if (not(Node_is_subject_to_velocity_dirichlet_bcs[l][0]))
            {
              // Stress contribution
              // -------------------
              // Pressure
              residuals[local_eqn] += p_bar_local * r * dtestfdx(l, 1) * W;

              // Shear stress
              residuals[local_eqn] -=
                visc_ratio * r *
                (grad_u_bar_local[1][0] +
                 this->Gamma[1] * grad_u_bar_local[0][1]) *
                dtestfdx(l, 0) * W;

              residuals[local_eqn] -= visc_ratio * r * (1.0 + this->Gamma[1]) *
                                      grad_u_bar_local[1][1] * dtestfdx(l, 1) *
                                      W;


              // Add singular convective terms
              // -----------------------------
              // Radial
              residuals[local_eqn] -= scaled_re * r * u_bar_local[0] *
                                      grad_u_bar_local[1][0] * testf[l] * W;

              residuals[local_eqn] -=
                scaled_re *
                (r * u_bar_local[0] * interpolated_dudx(1, 0) +
                 r * interpolated_u[0] * grad_u_bar_local[1][0]) *
                testf[l] * W;

              // Axial
              residuals[local_eqn] -= scaled_re * r * u_bar_local[1] *
                                      grad_u_bar_local[1][1] * testf[l] * W;

              residuals[local_eqn] -=
                scaled_re *
                (r * u_bar_local[1] * interpolated_dudx(1, 1) +
                 r * interpolated_u[1] * grad_u_bar_local[1][1]) *
                testf[l] * W;

              // Jacobian
              if (flag)
              {
                for (unsigned ss = 0; ss < n_sing; ss++)
                {
                  const int local_unknown = local_equation_number_C[ss];
                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      p_hat_local[ss] * r * dtestfdx(l, 1) * W;

                    jacobian(local_eqn, local_unknown) -=
                      visc_ratio * r *
                      (grad_u_hat_local[ss][1][0] +
                       this->Gamma[1] * grad_u_hat_local[ss][0][0]) *
                      dtestfdx(l, 0) * W;

                    jacobian(local_eqn, local_unknown) -=
                      visc_ratio * r * (1.0 + this->Gamma[1]) *
                      grad_u_hat_local[ss][1][1] * dtestfdx(l, 1) * W;

                    // Radial
                    jacobian(local_eqn, local_unknown) -=
                      scaled_re *
                      (r * u_hat_local[ss][0] * grad_u_bar_local[1][0] +
                       r * u_bar_local[0] * grad_u_hat_local[ss][1][0]) *
                      testf[l] * W;

                    jacobian(local_eqn, local_unknown) -=
                      scaled_re *
                      (r * u_hat_local[ss][0] * interpolated_dudx(1, 0) +
                       r * interpolated_u[0] * grad_u_hat_local[ss][1][0]) *
                      testf[l] * W;

                    // Axial
                    jacobian(local_eqn, local_unknown) -=
                      scaled_re *
                      (r * u_hat_local[ss][1] * grad_u_bar_local[1][1] +
                       r * u_bar_local[1] * grad_u_hat_local[ss][1][1]) *
                      testf[l] * W;

                    jacobian(local_eqn, local_unknown) -=
                      scaled_re *
                      (r * u_hat_local[ss][1] * interpolated_dudx(1, 1) +
                       r * interpolated_u[1] * grad_u_hat_local[ss][1][1]) *
                      testf[l] * W;
                  }
                }

                // Loop over the velocity shape functions again
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Radial
                  local_unknown = this->u_axi_nst_local_unknown(l2, 0);

                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) -=
                      scaled_re *
                      (r * u_bar_local[0] * dpsifdx(1, 0) +
                       r * psif[0] * grad_u_bar_local[1][0]) *
                      testf[l] * W;
                  }

                  // Axial
                  local_unknown = this->u_axi_nst_local_unknown(l2, 1);

                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) -=
                      scaled_re *
                      (r * u_bar_local[1] * dpsifdx(1, 1) +
                       r * psif[1] * grad_u_bar_local[1][1]) *
                      testf[l] * W;
                  }
                }
              } // End of Jacobian
            } // End of check of the Dirichlet status
          } // End of check of the pin status


          // ADDITIONAL VELOCITY DATA EQUATIONS
          // ----------------------------------

          for (unsigned i = 0; i < this->n_u_nst(); i++)
          {
            // Find its local equation number
            local_eqn = this->total_velocity_local_eqn(l, i);

            // If it is not pinned
            if (local_eqn >= 0)
            {
              residuals[local_eqn] +=
                (interpolated_u[i] -
                 (interpolated_u_tilde[i] + u_bar_local[i])) *
                testf[l] * W;

              // Jacobian
              if (flag)
              {
                // Singular contribution
                for (unsigned ss = 0; ss < n_sing; ss++)
                {
                  local_unknown = local_equation_number_C[ss];
                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) -=
                      u_hat_local[ss][i] * testf[l] * W;
                  }
                }

                // Velocity contribution
                // If at a non-zero degree of freedom add in the
                // entry
                // Loop over the velocity test functions
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  local_unknown = this->u_axi_nst_local_unknown(l2, i);
                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) -=
                      psif(l2) * testf[l] * W;
                  }
                }

                // Tilde velocity contribution
                // If at a non-zero degree of freedom add in the
                // entry
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  local_unknown = u_reconstructed_local_unknown(l2, i);
                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      psif(l2) * testf[l] * W;
                  }
                }
              }
            }
          }
        } // End of loop over test functions

        // CONTINUITY EQUATION
        //-------------------

        // Loop over the shape functions
        for (unsigned l = 0; l < n_pres; l++)
        {
          local_eqn = this->p_local_eqn(l);

          // If not a boundary conditions
          if (local_eqn >= 0)
          {
            // If not subject to Dirichlet BC
            if (not(Pressure_dof_is_subject_to_dirichlet_bc[l]))
            {
              residuals[local_eqn] +=
                (u_bar_local[0] + r * grad_u_bar_local[0][0] +
                 r * grad_u_bar_local[1][1]) *
                testp[l] * W;

              /*CALCULATE THE JACOBIAN*/
              if (flag)
              {
                // Loop over the singularities and add the
                // contributions of the additional
                // unknowns associated with them to the
                // jacobian if they are not pinned
                for (unsigned ss = 0; ss < n_sing; ss++)
                {
                  local_unknown = local_equation_number_C[ss];
                  if (local_unknown >= 0)
                  {
                    double sum = u_hat_local[ss][0];
                    for (unsigned k = 0; k < cached_dim; k++)
                    {
                      sum += r * grad_u_hat_local[ss][k][k];
                    }
                    jacobian(local_eqn, local_unknown) += sum * testp[l] * W;
                  }
                }
              } /*End of Jacobian calculation*/
            }
          } // End of if not boundary condition
        } // End of loop over l
      } // End of loop over integration points

      // VELOCITY DIRICHLET BCS
      //-----------------------
      Vector<double> u_bar_at_node(this->n_u_nst());
      Vector<Vector<double>> u_hat_at_node(n_sing);
      for (unsigned i = 0; i < n_sing; i++)
      {
        u_hat_at_node[i].resize(this->n_u_nst());
      }

      // Loop over the nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Find the global coordinate of the node
        Vector<double> global_coordinate(cached_dim);
        for (unsigned d = 0; d < cached_dim; d++)
        {
          global_coordinate[d] = this->raw_nodal_position(l, d);
        }

        // Get singular velocity at node
        u_bar_at_node = this->u_bar(global_coordinate);
        for (unsigned i = 0; i < n_sing; i++)
        {
          u_hat_at_node[i] = velocity_singular_function(i, global_coordinate);
        }

        // Loop over the velocity components
        for (unsigned d = 0; d < this->n_u_nst(); d++)
        {
          // Find its local equation number
          local_eqn = this->axi_momentum_local_eqn(l, d);
          // local_eqn = this->total_velocity_local_eqn(l, d);

          // If it is not pinned
          if (local_eqn >= 0)
          {
            // If it is a Dirichlet boundary condition
            if (Node_is_subject_to_velocity_dirichlet_bcs[l][d])
            {
              // Initialise the residual
              residuals[local_eqn] = 0.0;

              // Add the contribution of the nodal value
              residuals[local_eqn] +=
                this->nodal_value(l, u_index_axi_nst(l, d));

              // Add the contribution of the singularities (all of them,
              // summed)
              residuals[local_eqn] += u_bar_at_node[d];

              // Substract the imposed Dirichlet value
              residuals[local_eqn] -= Imposed_velocity_values_at_node[l][d];

              if (flag)
              {
                // Wipe the existing entries
                unsigned n_dof = this->ndof();
                for (unsigned j = 0; j < n_dof; j++)
                {
                  jacobian(local_eqn, j) = 0.0;
                }

                // Add diagonal entry
                jacobian(local_eqn, local_eqn) += 1.0;


                // Add derivative w.r.t. the Cs
                for (unsigned i = 0; i < n_sing; i++)
                {
                  // Find the contribution of the additional unknowns to
                  // the jacobian
                  local_unknown = local_equation_number_C[i];
                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) += u_hat_at_node[i][d];
                  }
                }
              }
            }
          }
        }
      }


      // PRESSURE DIRICHLET BCS
      //-----------------------

      // Loop over the pressure dofs
      for (unsigned l = 0; l < n_pres; l++)
      {
        // Find its local equation number
        local_eqn = this->p_local_eqn(l);

        // If it is not pinned
        if (local_eqn >= 0)
        {
          // If it is subject to a Dirichlet BC
          if (Pressure_dof_is_subject_to_dirichlet_bc[l])
          {
            // Find its global coordinate
            // This conversionly works for Taylor Hood type elements
            // but there's not much point assigning pressure dofs
            Node* p_nod_pt = this->node_pt(this->Pconv[l]);

            oomph_info << "Constrained pressure node: " << this->Pconv[l]
                       << " at: ";

            Vector<double> global_coordinate(cached_dim, 0.0);
            for (unsigned d = 0; d < cached_dim; d++)
            {
              global_coordinate[d] = p_nod_pt->x(d);
              oomph_info << global_coordinate[d] << " ";
            }
            oomph_info << std::endl;

            // Initialise its residual component
            residuals[local_eqn] = 0.0;

            // Add the contribution of the pressure unknown
            residuals[local_eqn] += p_axi_nst(l);

            // Add singular contributions
            residuals[local_eqn] += p_bar(global_coordinate);

            // Substract the imposed pressure value
            residuals[local_eqn] -= Imposed_value_at_pressure_dof[l];

            if (flag)
            {
              // Wipe the existing entries
              unsigned n_dof = this->ndof();
              for (unsigned j = 0; j < n_dof; j++)
              {
                jacobian(local_eqn, j) = 0.0;
              }

              // Add diagonal entry
              jacobian(local_eqn, local_eqn) += 1.0;

              // Add derivative w.r.t. the Cs
              for (unsigned i = 0; i < n_sing; i++)
              {
                // Find the contribution of the additional unknowns to
                // the jacobian
                local_unknown = local_equation_number_C[i];
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    pressure_singular_function(i, global_coordinate);
                }
              }
            }
          }
        }
      }
    }
  };


  //=======================================================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<SingularAxisymNavierStokesElement<ELEMENT>>
    : public virtual FaceGeometry<ELEMENT>
  {
  public:
    FaceGeometry() : FaceGeometry<ELEMENT>() {}
  };


  //=======================================================================
  /// Face geometry of the Face Geometry for element is the same as
  /// that for the underlying wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<FaceGeometry<SingularAxisymNavierStokesElement<ELEMENT>>>
    : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}
  };

} // namespace oomph

#endif
