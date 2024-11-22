#ifndef SINGULAR_NAVIER_STOKES_ELEMENTS_HEADER
#define SINGULAR_NAVIER_STOKES_ELEMENTS_HEADER

// oomph-lib includes
#include "navier_stokes.h"

// local includes
#include "singular_navier_stokes_solution_elements.h"

namespace oomph
{
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
  template<class BASIC_NAVIER_STOKES_ELEMENT>
  class SingularNavierStokesElement
    : public virtual BASIC_NAVIER_STOKES_ELEMENT,
      public virtual FiniteElement

  {
  private:
    bool IsAugmented;
    /// Vector of pointers to SingularNavierStokesSolutionElement objects
    Vector<SingularNavierStokesSolutionElement<
      SingularNavierStokesElement<BASIC_NAVIER_STOKES_ELEMENT>>*>
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
    SingularNavierStokesElement() : IsAugmented(false)
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

    /// Get 'flux' for Z2 error recovery:   Upper triangular entries
    /// in strain rate tensor.
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
#ifdef PARANOID
      unsigned num_entries =
        this->dim() + (this->dim() * (this->dim() - 1)) / 2;
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
      DenseMatrix<double> strainrate(this->dim());
      this->strain_rate(s, strainrate);

      // Pack into flux Vector
      unsigned icount = 0;

      // Start with diagonal terms
      for (unsigned i = 0; i < this->dim(); i++)
      {
        flux[icount] = strainrate(i, i);
        icount++;
      }

      // Off diagonals row by row
      for (unsigned i = 0; i < this->dim(); i++)
      {
        for (unsigned j = i + 1; j < this->dim(); j++)
        {
          flux[icount] = strainrate(i, j);
          icount++;
        }
      }
    }

    void strain_rate(const Vector<double>& s,
                     DenseMatrix<double>& strainrate) const
    {
      // Velocity gradient matrix
      DenseMatrix<double> dudx(this->dim());

      // Find out how many nodes there are in the element
      unsigned n_node = nnode();

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node, this->dim());

      // Call the derivatives of the shape functions
      dshape_eulerian(s, psi, dpsidx);

      // Initialise to zero
      for (unsigned i = 0; i < this->dim(); i++)
      {
        for (unsigned j = 0; j < this->dim(); j++)
        {
          dudx(i, j) = 0.0;
        }
      }

      // Loop over veclocity components
      for (unsigned i = 0; i < this->dim(); i++)
      {
        // Loop over derivative directions
        for (unsigned j = 0; j < this->dim(); j++)
        {
          // Loop over nodes
          for (unsigned l = 0; l < n_node; l++)
          {
            // Get the index at which the i-th velocity is stored
            unsigned u_nodal_index = u_reconstructed_index(l, i);
            dudx(i, j) += nodal_value(l, u_nodal_index) * dpsidx(l, j);
          }
        }
      }

      // Loop over veclocity components
      for (unsigned i = 0; i < this->dim(); i++)
      {
        // Loop over derivative directions
        for (unsigned j = 0; j < this->dim(); j++)
        {
          strainrate(i, j) = 0.5 * (dudx(i, j) + dudx(j, i));
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
    /// In derived multi-physics elements, this function should be overloaded
    /// to reflect the chosen storage scheme.
    virtual inline unsigned u_index_nst(const unsigned& n,
                                        const unsigned& i) const override
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
        index = BASIC_NAVIER_STOKES_ELEMENT::u_index_nst(n, i);
      }
      return index;
    }

    virtual inline unsigned u_reconstructed_index(const unsigned& n,
                                                  const unsigned& i) const
    {
      // i * [u]
      return i;
    }

    virtual inline unsigned momentum_index_nst(const unsigned& n,
                                               const unsigned& i) const
    {
      return u_index_nst(n, i);
      // return u_reconstructed_index(n, i);
    }

    virtual inline unsigned total_velocity_eqn_index(const unsigned& n,
                                                     const unsigned& i) const
    {
      return u_reconstructed_index(n, i);
      // return u_index_nst(n, i);
    }

    virtual inline int total_velocity_local_eqn(const unsigned& n,
                                                const unsigned& i) const
    {
      return nodal_local_eqn(n, total_velocity_eqn_index(n, i));
    }

    virtual inline int u_reconstructed_local_unknown(const unsigned& n,
                                                     const unsigned& i) const
    {
      return nodal_local_eqn(n, this->u_reconstructed_index(n, i));
    }

    double u_reconstructed(const unsigned& n, const unsigned& i) const
    {
      return nodal_value(n, this->u_reconstructed_index(n, i));
    }

    // /// Idea: Everyone outside the element used this function which relates
    // to the total velocity, and the internals swap the total velocity with
    // the tilde part. This means that only this element needs to know that it
    // is augmented and Dirchlet and Neumann boundary conditions would be
    // imposed automatically.
    // /// Compute vector of FE interpolated velocity u at local coordinate s
    // void interpolated_u_nst(const Vector<double>& s,
    //                         Vector<double>& veloc) const
    // {
    //   // Find number of nodes
    //   unsigned n_node = nnode();
    //   // Local shape function
    //   Shape psi(n_node);
    //   // Find values of shape function
    //   shape(s, psi);

    //   for (unsigned i = 0; i < DIM; i++)
    //   {
    //     // Initialise value of u
    //     veloc[i] = 0.0;
    //     // Loop over the local nodes and sum
    //     for (unsigned l = 0; l < n_node; l++)
    //     {
    //       // Index at which the nodal value is stored
    //       unsigned u_nodal_index = u_reconstructed_index(l, i);
    //       veloc[i] += nodal_value(l, u_nodal_index) * psi[l];
    //     }
    //   }
    // }
    // /// Return FE interpolated velocity u[i] at local coordinate s
    // double interpolated_u_nst(const Vector<double>& s, const unsigned& i)
    // const
    // {
    //   // Find number of nodes
    //   unsigned n_node = nnode();
    //   // Local shape function
    //   Shape psi(n_node);
    //   // Find values of shape function
    //   shape(s, psi);

    //   // Initialise value of u
    //   double interpolated_u = 0.0;
    //   // Loop over the local nodes and sum
    //   for (unsigned l = 0; l < n_node; l++)
    //   {
    //     // Get nodal index at which i-th velocity is stored
    //     unsigned u_nodal_index = u_reconstructed_index(l, i);
    //     interpolated_u += nodal_value(l, u_nodal_index) * psi[l];
    //   }

    //   return (interpolated_u);
    // }

    // virtual inline int momentum_local_eqn(const unsigned& n,
    //                                       const unsigned& i) const
    // {
    //   // Note: u_reconstructed_index references to total u

    //   return nodal_local_eqn(n, u_index_nst(n, i));
    //   // return nodal_local_eqn(n, this->u_reconstructed_index(n, i));
    // }

    /// Return the index at which the unknown pressure component
    /// is stored.
    // virtual inline int p_nodal_index_nst() const
    // {
    //   // dim * [u]
    //   return this->dim();
    // }


  protected:
    // Increase the required number of values by dim
    //    virtual unsigned required_nvalue(const unsigned& n) const
    //    {
    //      // Basic element
    //      // dim * [u] + 1 * [p]
    //
    //      // This element
    //      // dim * [u] + 1 * [p_tilde] + dim * [u_tilde]
    //
    //      return TTaylorHoodElement<2>::required_nvalue(n) + this->dim();
    //    }


  public:
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
      this->node_pt(i_node)->set_value(this->u_index_nst(i_node, i_dim), value);
    }

    void pin_momentum_eqn(const unsigned& i_node, const unsigned& i_dim)
    {
      this->node_pt(i_node)->pin(this->momentum_index_nst(i_node, i_dim));
    }

    void pin_fluid()
    {
      for (unsigned i_node = 0; i_node < this->nnode(); i_node++)
      {
        for (unsigned i_dim = 0; i_dim < this->dim(); i_dim++)
        {
          this->node_pt(i_node)->pin(this->momentum_index_nst(i_node, i_dim));
        }
        if (i_node < 3)
        {
          this->node_pt(i_node)->pin(this->p_nodal_index_nst());
        }
      }
    }

    void pin_total_velocity()
    {
      if (this->IsAugmented)
      {
        for (unsigned i_node = 0; i_node < this->nnode(); i_node++)
        {
          for (unsigned i_dim = 0; i_dim < this->dim(); i_dim++)
          {
            this->node_pt(i_node)->pin(
              this->total_velocity_eqn_index(i_node, i_dim));
          }
        }
      }
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
      SingularNavierStokesElement<BASIC_NAVIER_STOKES_ELEMENT>>*>
    c_equation_elements_pt()
    {
      return C_equation_elements_pt;
    }

    /// Add pointer to associated SingularNavierStokesSolutionElement that
    /// determines the value of the amplitude of the singular functions (and
    /// gives access to the singular functions). The unknown amplitude becomes
    /// external Data for this element so assign_eqn_numbers() must be
    /// called after this function has been called.
    void add_c_equation_element_pt(
      SingularNavierStokesSolutionElement<
        SingularNavierStokesElement<BASIC_NAVIER_STOKES_ELEMENT>>* c_pt)
    {
      // Add the element
      C_equation_elements_pt.push_back(c_pt);

      // Add the additional unknown of this object as external data in the
      // Navier-Stokes element
      this->add_external_data(c_pt->internal_data_pt(0));
    }

    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      this->fill_in_generic_residual_contribution_wrapped_nst(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1 and dummy mass
      // matrix
      this->fill_in_generic_residual_contribution_wrapped_nst(
        residuals, jacobian, 1);
      // FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    /// Overload the output function
    /// x, y, [z,] u, v, [w], p, u_fe, v_fe, [w_fe], p_fe,
    // u_sing, v_sing, [w_sing], p_sing
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      const bool IsStressOutputIncluded = true;
      const bool IsSingularOutputInclude = true;

      // Find the dimension of the problem
      unsigned cached_dim = this->dim();

      // Vector of local coordinates
      Vector<double> s(cached_dim);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);

        Vector<double> x(cached_dim);
        for (unsigned i = 0; i < cached_dim; i++)
        {
          outfile << this->interpolated_x(s, i) << " ";
          x[i] = this->interpolated_x(s, i);
        }

        // Total velocity
        Vector<double> velocity(cached_dim, 0.0);

        if (IsAugmented)
        {
          velocity = interpolated_u_reconstructed(s);
        }
        else
        {
          this->interpolated_u_nst(s, velocity);
        }
        for (unsigned i = 0; i < cached_dim; i++)
        {
          outfile << velocity[i] << " ";
        }

        // "Total" pressure
        outfile << this->interpolated_p_nst(s) << " ";

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

        // Finite element Velocity
        Vector<double> velocity_fe_only(cached_dim, 0.0);
        if (this->IsAugmented)
        {
          this->interpolated_u_nst(s, velocity_fe_only);
        }

        for (unsigned i = 0; i < cached_dim; i++)
        {
          outfile << velocity_fe_only[i] << " ";
        }

        // Singular Velocity
        if (IsSingularOutputInclude)
        {
          Vector<double> velocity_bar(cached_dim, 0.0);
          if (this->IsAugmented)
          {
            velocity_bar = u_bar(x);
          }
          for (unsigned i = 0; i < cached_dim; i++)
          {
            outfile << velocity_bar[i] << " ";
          }

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
        } // End of augmented

        // Output the stored error
        outfile << this->get_error() << " ";

        // Output the element size
        outfile << this->size() << " ";

        // Output the continuity residual
        outfile << continuity_residual(s) << " ";

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
      Vector<double> exact_soln(cached_dim + 1);
      Vector<double> computed_soln(cached_dim + 1);

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
        Vector<double> u_comp = interpolated_u_reconstructed(s);
        for (unsigned i = 0; i < cached_dim; i++)
        {
          computed_soln[i] = u_comp[i];
        }
        unsigned hi_limit = cached_dim;
        if (include_pressure)
        {
          computed_soln[cached_dim] = this->interpolated_p_nst(s);
          hi_limit = cached_dim + 1;
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


    /// Overloaded version of the interpolated velocity solution including
    /// the singular contributions
    Vector<double> interpolated_u_reconstructed(const Vector<double>& s)
    {
      // Find the dimension of the problem
      unsigned cached_dim = this->dim();

      // Find number of nodes
      const unsigned n_node = this->nnode();

      // Initialise value of u
      Vector<double> interpolated_u(cached_dim, 0.0);

      // Local shape function
      Shape psif(n_node);

      // Find values of shape function
      this->shape(s, psif);

      // Loop over the spatial directions
      for (unsigned d = 0; d < cached_dim; d++)
      {
        // Loop over the local nodes and sum
        for (unsigned j = 0; j < n_node; j++)
        {
          interpolated_u[d] += u_reconstructed(j, d) * psif[j];
        }
      }
      return interpolated_u;
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
      double interpolated_p = this->interpolated_p_nst(s);
      DenseMatrix<double> interpolated_dudx(cached_dim, cached_dim, 0.0);
      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the spatial directions
        for (unsigned i = 0; i < cached_dim; i++)
        {
          for (unsigned j = 0; j < cached_dim; j++)
          {
            const double u_value = this->u_reconstructed(l, i);
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
        stress_tensor(i, i) += interpolated_p;
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

    double continuity_residual(const Vector<double>& s)
    {
      // Find the dimension of the problem
      const unsigned cached_dim = this->dim();
      DenseMatrix<double> stress_tensor(cached_dim, cached_dim, 0.0);

      // Find number of nodes
      const unsigned n_node = this->nnode();

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
          for (unsigned j = 0; j < cached_dim; j++)
          {
            const double u_value = this->u_nst(l, i);
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }
      if (this->IsAugmented)
      {
        Vector<double> x(cached_dim);
        for (unsigned i = 0; i < cached_dim; i++)
        {
          x[i] = this->interpolated_x(s, i);
        }
        Vector<Vector<double>> grad_u_bar_local = this->grad_u_bar(x);
        // Loop over the spatial directions
        for (unsigned i = 0; i < cached_dim; i++)
        {
          for (unsigned j = 0; j < cached_dim; j++)
          {
            interpolated_dudx(i, j) += grad_u_bar_local[i][j];
          }
        }
      }

      return interpolated_dudx(0, 0) + interpolated_dudx(1, 1);
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
      unsigned cached_dim = this->dim();
      Vector<double> sum(cached_dim, 0.0);
      for (unsigned s = 0; s < n_sing; s++)
      {
        Vector<double> u_bar_local = C_equation_elements_pt[s]->u_bar(x);
        for (unsigned i = 0; i < cached_dim; i++)
        {
          sum[i] += u_bar_local[i];
        }
      }
      return sum;
    }

    /// Evaluate sum of all velocity singular fcts
    /// (incl. the amplitude) at Eulerian position x
    double u_bar(const Vector<double>& x, const unsigned& i) const
    {
      // Find the number of singularities
      unsigned n_sing = C_equation_elements_pt.size();

      // Find the dimension of the problem
      double sum = 0.0;
      for (unsigned s = 0; s < n_sing; s++)
      {
        Vector<double> u_bar_local = C_equation_elements_pt[s]->u_bar(x);
        sum += u_bar_local[i];
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
    /// (incl. the amplitudes) at Eulerian position x: grad[i][j] = du_i/dx_j
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


    /// Evaluate i-th "raw" velocity singular function at Eulerian position x
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

  private:
    void fill_in_generic_residual_contribution_additional_terms(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      // Find the dimension of the problem
      unsigned cached_dim = this->dim();

      // Do all the singular functions satisfy the Stokes eqn?
      bool all_singular_functions_satisfy_stokes_equation = true;

      // Find the number of singularities
      unsigned n_sing = C_equation_elements_pt.size();

      // Find the local equation number of the additional unknowns
      Vector<int> local_equation_number_C(n_sing);
      for (unsigned i = 0; i < n_sing; i++)
      {
        local_equation_number_C[i] = this->external_local_eqn(i, 0);
        if (!(C_equation_elements_pt[i]
                ->singular_function_satisfies_stokes_equation()))
        {
          all_singular_functions_satisfy_stokes_equation = false;
        }
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
        double J = this->dshape_and_dtest_eulerian_at_knot_nst(
          ipt, psif, dpsifdx, testf, dtestfdx);

        // Call the pressure shape and test functions
        this->pshape_nst(s, psip, testp);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Initialise the global coordinate and velocity
        Vector<double> interpolated_x(cached_dim, 0.0);
        Vector<double> interpolated_u_tilde(cached_dim, 0.0);
        Vector<double> interpolated_u_reconstructed(cached_dim, 0.0);
        DenseMatrix<double> interpolated_dudx(cached_dim, cached_dim, 0.0);

        // Calculate the global coordinate associated with s
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < cached_dim; i++)
          {
            const double u_value = nodal_value(l, u_index_nst(l, i));
            interpolated_u_tilde[i] += u_value * psif[l];
            for (unsigned j = 0; j < cached_dim; j++)
            {
              interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
            }
            interpolated_u_reconstructed[i] += u_reconstructed(l, i) * psif[l];
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
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

        // MOMENTUM EQUATIONS
        //-------------------

        // Loop over the velocity test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the velocity components
          for (unsigned i = 0; i < cached_dim; i++)
          {
            // Find its local equation number
            local_eqn = this->momentum_local_eqn(l, i);

            // If it is not pinned
            if (local_eqn >= 0)
            {
              // If it is not a Dirichlet BC
              if (not(Node_is_subject_to_velocity_dirichlet_bcs[l][i]))
              {
                // Stress contribution
                // -------------------
                // if (!all_singular_functions_satisfy_stokes_equation)
                {
                  residuals[local_eqn] += p_bar_local * dtestfdx(l, i) * W;
                  for (unsigned k = 0; k < cached_dim; k++)
                  {
                    residuals[local_eqn] -=
                      visc_ratio *
                      (grad_u_bar_local[i][k] +
                       this->Gamma[i] * grad_u_bar_local[k][i]) *
                      dtestfdx(l, k) * W;
                  }

                  // Jacobian
                  if (flag)
                  {
                    for (unsigned ss = 0; ss < n_sing; ss++)
                    {
                      const int local_unknown = local_equation_number_C[ss];
                      if (local_unknown >= 0)
                      {
                        jacobian(local_eqn, local_unknown) +=
                          p_hat_local[ss] * dtestfdx(l, i) * W;
                        for (unsigned k = 0; k < cached_dim; k++)
                        {
                          jacobian(local_eqn, local_unknown) -=
                            visc_ratio *
                            (grad_u_hat_local[ss][i][k] +
                             this->Gamma[i] * grad_u_hat_local[ss][k][i]) *
                            dtestfdx(l, k) * W;
                        }
                      }
                    }
                  }
                }

                // Nonlinear term. Always add (unless Re=0)
                //-----------------------------------------
                if (scaled_re > 0.0)
                {
                  // Add singular convective terms
                  double sum = 0.0;
                  for (unsigned k = 0; k < cached_dim; k++)
                  {
                    sum += u_bar_local[k] * (grad_u_bar_local[i][k] +
                                             interpolated_dudx(i, k)) +
                           interpolated_u_tilde[k] * grad_u_bar_local[i][k];
                  }
                  residuals[local_eqn] -= scaled_re * sum * testf[l] * W;

                  // Calculate the jacobian
                  //-----------------------
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
                        double sum = 0.0;
                        for (unsigned k = 0; k < cached_dim; k++)
                        {
                          sum +=
                            u_hat_local[ss][k] * (grad_u_bar_local[i][k] +
                                                  interpolated_dudx(i, k)) +
                            (u_bar_local[k] + interpolated_u_tilde[k]) *
                              grad_u_hat_local[ss][i][k];
                        }
                        jacobian(local_eqn, local_unknown) -=
                          scaled_re * sum * testf[l] * W;
                      }
                    }


                    // Loop over the velocity shape functions again
                    for (unsigned l2 = 0; l2 < n_node; l2++)
                    {
                      // Loop over the velocity components again
                      for (unsigned i2 = 0; i2 < cached_dim; i2++)
                      {
                        // If at a non-zero degree of freedom add in the
                        // entry
                        local_unknown = this->u_local_unknown(l2, i2);
                        if (local_unknown >= 0)
                        {
                          double sum = 0.0;
                          if (i == i2)
                          {
                            for (unsigned k = 0; k < cached_dim; k++)
                            {
                              sum += u_bar_local[k] * dpsifdx(l2, k);
                            }
                          }
                          sum += psif(l2) * grad_u_bar_local[i][i2];

                          // Add contribution to Elemental Matrix
                          jacobian(local_eqn, local_unknown) -=
                            scaled_re * sum * testf[l] * W;
                        }
                      }
                    }
                  }
                }

              } // End of check of the Dirichlet status
            } // End of check of the pin status
          } // End of loop over velocity components
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
              // Source term
              double aux = 0.0;

              // Loop over velocity components
              for (unsigned k = 0; k < cached_dim; k++)
              {
                aux += grad_u_bar_local[k][k];
              }

              residuals[local_eqn] += aux * testp[l] * W;

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
                    double sum = 0.0;
                    for (unsigned k = 0; k < cached_dim; k++)
                    {
                      sum += grad_u_hat_local[ss][k][k];
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
      Vector<double> u_bar_at_node(cached_dim);
      Vector<Vector<double>> u_hat_at_node(n_sing);
      for (unsigned i = 0; i < n_sing; i++)
      {
        u_hat_at_node[i].resize(cached_dim);
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
        for (unsigned d = 0; d < cached_dim; d++)
        {
          // Find its local equation number
          local_eqn = this->momentum_local_eqn(l, d);
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
              residuals[local_eqn] += this->u_nst(l, d);

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
            residuals[local_eqn] += this->p_nst(l);

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

    void fill_in_generic_residual_contribution_total_velocity(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      // Find the dimension of the problem
      unsigned cached_dim = this->dim();

      // Do all the singular functions satisfy the Stokes eqn?
      bool all_singular_functions_satisfy_stokes_equation = true;

      // Find the number of singularities
      unsigned n_sing = C_equation_elements_pt.size();

      // Find the local equation number of the additional unknowns
      Vector<int> local_equation_number_C(n_sing);
      for (unsigned i = 0; i < n_sing; i++)
      {
        local_equation_number_C[i] = this->external_local_eqn(i, 0);
        if (!(C_equation_elements_pt[i]
                ->singular_function_satisfies_stokes_equation()))
        {
          all_singular_functions_satisfy_stokes_equation = false;
        }
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
        double J = this->dshape_and_dtest_eulerian_at_knot_nst(
          ipt, psif, dpsifdx, testf, dtestfdx);

        // Call the pressure shape and test functions
        this->pshape_nst(s, psip, testp);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Initialise the global coordinate and velocity
        Vector<double> interpolated_x(cached_dim, 0.0);
        Vector<double> interpolated_u_tilde(cached_dim, 0.0);
        Vector<double> interpolated_u_reconstructed(cached_dim, 0.0);
        DenseMatrix<double> interpolated_dudx(cached_dim, cached_dim, 0.0);

        // Calculate the global coordinate associated with s
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < cached_dim; i++)
          {
            const double u_value = nodal_value(l, u_index_nst(l, i));
            interpolated_u_tilde[i] += u_value * psif[l];
            for (unsigned j = 0; j < cached_dim; j++)
            {
              interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
            }
            interpolated_u_reconstructed[i] += u_reconstructed(l, i) * psif[l];
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
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

        // MOMENTUM EQUATIONS
        //-------------------

        // Loop over the velocity test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the velocity components
          for (unsigned i = 0; i < cached_dim; i++)
          {
            // Additional velocity data
            // ------------------------

            // Find its local equation number
            local_eqn = this->total_velocity_local_eqn(l, i);

            // If it is not pinned
            if (local_eqn >= 0)
            {
              // residuals[local_eqn] +=
              //   (interpolated_u_reconstructed[i] -
              //    (interpolated_u_tilde[i] + u_bar_local[i])) *
              //   testf[l] * W;
              Vector<double> pos_n(2, 0.0);
              for (unsigned k = 0; k < 2; k++)
              {
                pos_n[k] = this->nodal_position(l, k);
              }

              residuals[local_eqn] +=
                (u_reconstructed(l, i) -
                 (nodal_value(l, u_index_nst(l, i)) + u_bar(pos_n, i)));

              // Jacobian
              if (flag)
              {
                Vector<Vector<double>> u_hat_local2(n_sing);
                for (unsigned s = 0; s < n_sing; s++)
                {
                  u_hat_local2[s] = this->velocity_singular_function(s, pos_n);
                }

                // Singular contribution
                for (unsigned ss = 0; ss < n_sing; ss++)
                {
                  local_unknown = local_equation_number_C[ss];
                  if (local_unknown >= 0)
                  {
                    jacobian(local_eqn, local_unknown) -= u_hat_local2[ss][i];
                  }
                }

                // Velocity contribution
                // If at a non-zero degree of freedom add in the
                // entry
                // Loop over the velocity test functions
                local_unknown = this->u_local_unknown(l, i);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -= 1;
                }

                // Tilde velocity contribution
                // If at a non-zero degree of freedom add in the
                // entry
                local_unknown = u_reconstructed_local_unknown(l, i);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) += 1;
                }
              }
            }
          } // End of loop over velocity components
        } // End of loop over test functions

      } // End of loop over integration points
    }

    /// Overloaded fill-in function
    void fill_in_generic_residual_contribution_wrapped_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      // Get the contribution from the underlying wrapped element first
      BASIC_NAVIER_STOKES_ELEMENT::fill_in_generic_residual_contribution_nst(
        residuals, jacobian, GeneralisedElement::Dummy_matrix, flag);

      if (this->IsAugmented)
      {
        fill_in_generic_residual_contribution_additional_terms(
          residuals, jacobian, flag);
        fill_in_generic_residual_contribution_total_velocity(
          residuals, jacobian, flag);
      }
    }
  };


  //=======================================================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<SingularNavierStokesElement<ELEMENT>>
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
  class FaceGeometry<FaceGeometry<SingularNavierStokesElement<ELEMENT>>>
    : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}
  };
} // namespace oomph
#endif
