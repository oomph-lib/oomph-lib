#ifndef LINEARISED_ELASTIC_AXISYM_FLUID_INTERFACE_ELEMENT_HEADER
#define LINEARISED_ELASTIC_AXISYM_FLUID_INTERFACE_ELEMENT_HEADER

#include <algorithm>

//#include "overlaying_my_linear_element.h"
#include "overlaying_Tlinear_axisym_ns_pvd_elements.h"
#include "linearised_axisymmetric_fluid_interface_elements.h"
#include "debug_jacobian_elements.h"

namespace oomph
{
  template<class BULK_ELEMENT>
  class LinearisedElasticAxisymmetricFluidInterfaceElement
    : public virtual LinearisedAxisymmetricFluidInterfaceElement<BULK_ELEMENT>,
      public virtual TElement<1, 3>,
      public virtual DebugJacobianFiniteElement
  {
  private:
    TimeStepper* Time_stepper_pt;

    /// Storage for the location of the Lagrange multiplier
    /// (If other additional values have been added we need
    /// to add the Lagrange multiplier at the end)
    Vector<Vector<unsigned>> Lagrange_index;

    /// Storage for the location of the Lagrange multiplier
    /// (If other additional values have been added we need
    /// to add the Lagrange multiplier at the end)
    Vector<unsigned> First_index;

    /// Return the index at which the lagrange multiplier is
    /// stored at the n-th node
    inline unsigned lagrange_index(const unsigned& n, const unsigned& i)
    {
      return this->Lagrange_index[n][i];
    }

    virtual double lagrange_multiplier(const unsigned& n, const unsigned& i)
    {
      return this->nodal_value(n, this->Lagrange_index[n][i]);
    }

    /// Equation number of the kinematic BC associated with node j.
    /// (This is the equation for the Lagrange multiplier)
    int kinematic_local_eqn(const unsigned& n, const unsigned& i)
    {
      // Get the index of the nodal value associated with Lagrange multiplier
      return this->nodal_local_eqn(n, this->lagrange_index(n, i));
    }

    /// Access function that returns the local equation number
    /// for the i in {0 = R, 1 = Z}, j in {0 = C, 1 = S} displacement equation
    /// that corresponds to the n-th local node. This must be overloaded by
    /// specific interface elements and depends on the method for handing the
    /// free-surface deformation.
    int displacement_local_eqn(const unsigned& n,
                               const unsigned& i,
                               const unsigned& j)
    {
      return this->nodal_local_eqn(n, u_index_linear_elasticity(n, i, j));
    }

    // /// Hijacking the kinematic condition corresponds to hijacking the
    // /// variables associated with the Lagrange multipliers that are assigned
    // /// on construction of this element.
    // void hijack_kinematic_conditions(const Vector<unsigned>&
    // bulk_node_number)
    // {
    //   // Loop over all the nodes that are passed in
    //   for (Vector<unsigned>::const_iterator it = bulk_node_number.begin();
    //        it != bulk_node_number.end();
    //        ++it)
    //   {
    //     // Hijack the appropriate value and delete the returned Node
    //     delete this->hijack_nodal_value(*it, this->lagrange_index(*it));
    //   }
    // }


  public:
    /// Return the index at which the i-th unknown displacement
    /// component is stored. The default value, i, is appropriate for
    /// single-physics problems.
    virtual inline unsigned u_index_linear_elasticity(const unsigned& n,
                                                      const unsigned& i,
                                                      const unsigned& j) const
    {
      return i * 2 + j;
    }

    // Pin displacement values: n in {0,...,nnode()}, i in {R, Z}, j in {C, S}
    void pin_Xhat(const unsigned& n, const unsigned& i, const unsigned& j)
    {
      this->node_pt(n)->pin(u_index_linear_elasticity(n, i, j));
    }

    // Set displacement values: n in {0,...,nnode()}, i in {R, Z}, j in {C, S}
    void set_value_Xhat(const unsigned& n,
                        const unsigned& i,
                        const unsigned& j,
                        const double& value)
    {
      this->node_pt(n)->set_value(u_index_linear_elasticity(n, i, j), value);
    }

    void pin_lagrange_multiplier(const unsigned& n, const unsigned& i)
    {
      this->node_pt(n)->pin(this->lagrange_index(n, i));
    }

    void pin_lagrange_multiplier(const unsigned& n,
                                 const unsigned& i,
                                 const double& value)
    {
      this->node_pt(n)->pin(this->lagrange_index(n, i));
      this->node_pt(n)->set_value(this->lagrange_index(n, i), value);
    }

    // Get normal vector to the end point corresponding to face_index (-1,1),
    // i = 0 for the cosine component and i = 1 for the sine component
    Vector<double> displacement_gradient(const int& face_index,
                                         const unsigned i)
    {
      Vector<double> gradient(2, 0.0);
      Vector<double> s(1, 0.0);
      if (face_index == -1)
      {
        s[0] = 0;
      }
      else if (face_index == 1)
      {
        s[0] = 1;
      }
      else
      {
        throw OomphLibError("Specified face index is incorrect\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      const unsigned n_node = this->nnode();
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);
      this->dshape_local(s, psif, dpsifds);

      Vector<double> interpolated_t1(2, 0.0);
      double interpolated_dRds = 0.0;
      double interpolated_dZds = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        double dpsifds_ = dpsifds(l, 0);

        // Calculate the tangent vector
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_t1[i] += this->nodal_position(l, i) * dpsifds_;
        }

        interpolated_dRds += this->Xhat(l, 0, i) * dpsifds_;
        interpolated_dZds += this->Xhat(l, 1, i) * dpsifds_;
      }

      const double tlength = interpolated_t1[0] * interpolated_t1[0] +
                             interpolated_t1[1] * interpolated_t1[1];
      // Set the Jacobian of the line element
      const double J = sqrt(tlength);

      gradient[0] = interpolated_dRds / J;
      gradient[1] = interpolated_dZds / J;

      return gradient;
    }

    // Get normal vector to the end point corresponding to face_index (-1,1),
    // i = 0 for the cosine component and i = 1 for the sine component
    Vector<double> normal(const int& face_index)
    {
      Vector<double> gradient(2, 0.0);
      Vector<double> s(1, 0.0);
      if (face_index == -1)
      {
        s[0] = 0;
      }
      else if (face_index == 1)
      {
        s[0] = 1;
      }
      else
      {
        throw OomphLibError("Specified face index is incorrect\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }


      Vector<double> unit_normal(2, 0.0);
      Vector<Vector<double>> tang_vec(1);
      tang_vec[0].resize(2, 0.0);
      this->continuous_tangent_and_outer_unit_normal(s, tang_vec, unit_normal);

      return unit_normal;
    }

    // Get normal vector to the end point corresponding to face_index (-1,1),
    // i = 0 for the cosine component and i = 1 for the sine component
    Vector<double> tangent(const int& face_index)
    {
      Vector<double> gradient(2, 0.0);
      Vector<double> s(1, 0.0);
      if (face_index == -1)
      {
        s[0] = 0;
      }
      else if (face_index == 1)
      {
        s[0] = 1;
      }
      else
      {
        throw OomphLibError("Specified face index is incorrect\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      Vector<double> unit_normal(2, 0.0);
      Vector<Vector<double>> tang_vec(1);
      tang_vec[0].resize(2, 0.0);
      this->continuous_tangent_and_outer_unit_normal(s, tang_vec, unit_normal);
      Vector<double> continuous_t1 = tang_vec[0];

      const unsigned n_node = nnode();
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);
      this->dshape_local(s, psif, dpsifds);
      Vector<double> interpolated_t1(2, 0.0);

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        double dpsifds_ = dpsifds(l, 0);

        // Calculate the tangent vector
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_t1[i] += this->nodal_position(l, i) * dpsifds_;
        }
      }

      bool is_element_reversed = (continuous_t1[0] * interpolated_t1[0] +
                                  continuous_t1[1] * interpolated_t1[1]) < 0;
      if (is_element_reversed)
      {
        continuous_t1[0] = -continuous_t1[0];
        continuous_t1[1] = -continuous_t1[1];
      }

      return continuous_t1;
    }

  public:
    bool is_element_reversed()
    {
      Vector<double> s(1, 0.0);
      const unsigned n_node = nnode();
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);
      this->dshape_local(s, psif, dpsifds);
      Vector<double> interpolated_t1(2, 0.0);
      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        double dpsifds_ = dpsifds(l, 0);

        // Calculate the tangent vector
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_t1[i] += this->nodal_position(l, i) * dpsifds_;
        }
      }

      Vector<double> unit_normal(2, 0.0);
      Vector<Vector<double>> tang_vec(1);
      tang_vec[0].resize(2, 0.0);
      this->continuous_tangent_and_outer_unit_normal(s, tang_vec, unit_normal);
      Vector<double> continuous_t1 = tang_vec[0];

      return (continuous_t1[0] * interpolated_t1[0] +
              continuous_t1[1] * interpolated_t1[1]) < 0;
    }

    // Displacement values: n in {0,...,nnode()}, i in {R, Z}, j in {C, S}
    virtual double Xhat(const unsigned& n, const unsigned& i, const unsigned& j)
    {
      return this->nodal_value(n, u_index_linear_elasticity(n, i, j));
    }

    // Time derivative of displacement values: n in {0,...,nnode()}, i in {R,
    // Z}, j in {C, S}
    virtual double dXhatdt(const unsigned& n,
                           const unsigned& i,
                           const unsigned& j)
    {
      // Get the data's positional timestepper
      TimeStepper* time_stepper_pt = this->node_pt(n)->time_stepper_pt();

      // Initialise dXhat/dt
      double dXhatdt = 0.0;

      // Loop over the timesteps, if there is a non Steady timestepper
      if (!time_stepper_pt->is_steady())
      {
        // Get the nodal index
        const unsigned u_nodal_index = u_index_linear_elasticity(n, i, j);

        // Number of timsteps (past & present)
        const unsigned n_time = time_stepper_pt->ntstorage();

        // Add the contributions to the time derivative
        for (unsigned t = 0; t < n_time; t++)
        {
          dXhatdt +=
            time_stepper_pt->weight(1, t) * nodal_value(t, n, u_nodal_index);
        }
      }

      return dXhatdt;
    }

  public:
    LinearisedElasticAxisymmetricFluidInterfaceElement(
      FiniteElement* const& element_pt,
      const int& face_index,
      TimeStepper* const& time_stepper_pt,
      const unsigned& id = 0)
      : LinearisedAxisymmetricFluidInterfaceElement<BULK_ELEMENT>(),
        TElement<1, 3>(),
        Time_stepper_pt(time_stepper_pt)
    {
      // Attach the geometrical information to the element
      // This function also assigned nbulk_value from required_nvalue of the
      // bulk element
      element_pt->build_face_element(face_index, this);

#ifdef PARANOID
      // Is it refineable
      RefineableElement* ref_el_pt =
        dynamic_cast<RefineableElement*>(element_pt);
      if (ref_el_pt != 0)
      {
        if (this->has_hanging_nodes())
        {
          throw OomphLibError(
            "This flux element will not work correctly if nodes are hanging\n",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif

      // Find the index at which the velocity unknowns are stored
      // from the bulk element and resize the local storage scheme
      BULK_ELEMENT* cast_element_pt = dynamic_cast<BULK_ELEMENT*>(element_pt);
      const unsigned n_u_index = cast_element_pt->n_u_lin_axi_nst();
      this->U_index_interface.resize(n_u_index);
      for (unsigned i = 0; i < n_u_index; i++)
      {
        this->U_index_interface[i] = cast_element_pt->u_index_lin_axi_nst(i);
      }

      // Read out the number of nodes on the face
      unsigned n_node_face = this->nnode();

      // Set the additional data values in the face
      // There is always also one additional values at each node --- the
      // Lagrange multiplier
      Vector<unsigned> additional_data_values(n_node_face);
      for (unsigned n = 0; n < n_node_face; n++)
      {
        // Now add one to the addtional values at every single node
        additional_data_values[n] = 2;
      }

      // Now add storage for Lagrange multipliers and set the map containing
      // the position of the first entry of this face element's
      // additional values.
      this->add_additional_values(additional_data_values, id);

      // Now I can just store the lagrange index offset to give the storage
      // location of the nodes
      First_index.resize(n_node_face);
      Lagrange_index.resize(n_node_face);
      for (unsigned n = 0; n < n_node_face; ++n)
      {
        Lagrange_index[n].resize(2);
        const unsigned first_index =
          dynamic_cast<BoundaryNodeBase*>(this->node_pt(n))
            ->index_of_first_value_assigned_by_face_element(id);
        First_index[n] = first_index;
        for (unsigned j = 0; j < 2; j++)
        {
          Lagrange_index[n][j] = first_index + j;
        }
      }

      // Loop over bulk nodes
      const unsigned n_node_bulk = element_pt->nnode();
      for (unsigned n = 0; n < n_node_bulk; n++)
      {
        // Test if node is in element
        if (std::find(this->Bulk_node_number.begin(),
                      this->Bulk_node_number.end(),
                      n) == this->Bulk_node_number.end())
        {
          // Add node
          add_external_data(element_pt->node_pt(n));
        }
      }
    }

    /// Calculate the residuals by calling the generic residual contribution.
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Add the residual contributions
      this->fill_in_generic_residual_contribution_interface(
        residuals,
        GeneralisedElement::Dummy_matrix,
        GeneralisedElement::Dummy_matrix,
        0);
    }

    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      this->fill_in_generic_residual_contribution_interface(
        residuals, jacobian, GeneralisedElement::Dummy_matrix, 1);
      // FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    // void fill_in_contribution_to_jacobian_and_mass_matrix(
    //   Vector<double>& residuals,
    //   DenseMatrix<double>& jacobian,
    //   DenseMatrix<double>& mass_matrix)
    //{
    //   FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    //   const unsigned n_dofs = this->ndof();
    //   Vector<double> dummy_residuals(n_dofs);
    //   DenseMatrix<double> dummy_jacobian(n_dofs, n_dofs);
    //   fill_in_generic_residual_contribution_interface(
    //     dummy_residuals, dummy_jacobian, mass_matrix, 2);
    // }
    /// Add the elemental contribution to the jacobian and mass matrices
    /// and the residuals vector. Note that
    /// this function will NOT initialise the residuals vector or the jacobian
    /// matrix. It must be called after the residuals vector and
    /// jacobian matrix have been initialised to zero. The default
    /// is to use finite differences to calculate the jacobian
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Add the contribution to the residuals
      this->fill_in_contribution_to_residuals(residuals);

      // Allocate storage for the full residuals (residuals of entire
      // element)
      unsigned n_dof = this->ndof();
      Vector<double> full_residuals(n_dof);
      // Get the residuals for the entire element
      this->get_residuals(full_residuals);

      // Make the timestepper steady
      Time_stepper_pt->make_steady();

      // Use the fill in contribution to jacobian
      fill_in_contribution_to_jacobian(full_residuals, jacobian);

      Time_stepper_pt->undo_make_steady();

      DenseMatrix<double> unsteady_jacobian(n_dof, n_dof);
      fill_in_contribution_to_jacobian(full_residuals, unsteady_jacobian);

      const double dt = Time_stepper_pt->time_pt()->dt();
      for (unsigned i = 0; i < n_dof; i++)
      {
        for (unsigned j = 0; j < n_dof; j++)
        {
          /// The 2/3 is due to the BDF<2> scheme.
          mass_matrix(j, i) +=
            (2.0 / 3.0) * dt * (-unsteady_jacobian(j, i) + jacobian(j, i));
        }
      }
    }

    // virtual void add_additional_residual_contributions(
    //   Vector<double>& residuals,
    //   DenseMatrix<double>& jacobian,
    //   const unsigned& flag,
    //   const Shape& psif,
    //   const DShape& dpsifds,
    //   const Vector<double>& interpolated_n,
    //   const double& r,
    //   const double& W,
    //   const double& J)
    // {
    //   const unsigned n_node = this->nnode();

    //   Vector<double> interpolated_lagrange_multiplier(2, 0.0);
    //   for (unsigned l = 0; l < n_node; l++)
    //   {
    //     const double psif_ = psif(l);
    //     for (unsigned i = 0; i < 2; i++)
    //     {
    //       interpolated_lagrange_multiplier[i] +=
    //         nodal_value(l, lagrange_index(l, i)) * psif_;
    //     }
    //   }

    //   int local_eqn = 0, local_unknown = 0;
    //   // Loop over the test functions
    //   for (unsigned l = 0; l < n_node; l++)
    //   {
    //     // Get local equation number of first velocity value at this node
    //     // RC
    //     local_eqn = this->nodal_local_eqn(
    //       l, First_index[l] + u_index_linear_elasticity(0, 0));

    //     if (local_eqn >= 0)
    //     {
    //       residuals[local_eqn] += interpolated_lagrange_multiplier[0]*
    //     }
    //   }
    // }


    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      const bool has_output_derivatives = true;
      const unsigned el_dim = this->dim();
      const unsigned n_dim = this->nodal_dimension();
      const unsigned n_node = this->nnode();
      // Set output Vector
      Vector<double> s(el_dim);

      // Tecplot header info
      // outfile << tecplot_zone_string(n_plot);
      // bool is_reversed = this->is_element_reversed();

      // Loop over plot points
      unsigned num_plot_points = nplot_points(n_plot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        get_s_plot(iplot, n_plot, s);

        Shape psif(n_node);
        DShape dpsifds(n_node, 1);
        this->dshape_local(s, psif, dpsifds);
        Vector<double> interpolated_t1(2, 0.0);
        Vector<double> unit_normal(2, 0.0);
        Vector<double> continuous_t1(2, 0.0);

        // Initialise value of u
        double interpolated_UC = 0.0;
        double interpolated_US = 0.0;
        double interpolated_VC = 0.0;
        double interpolated_VS = 0.0;
        double interpolated_WC = 0.0;
        double interpolated_WS = 0.0;
        Vector<double> interpolated_dUds(6, 0.0);

        double interpolated_RC = 0.0;
        double interpolated_RS = 0.0;
        double interpolated_ZC = 0.0;
        double interpolated_ZS = 0.0;

        double interpolated_dRCds = 0.0;
        double interpolated_dRSds = 0.0;
        double interpolated_dZCds = 0.0;
        double interpolated_dZSds = 0.0;

        double interpolated_dRCdt = 0.0;
        double interpolated_dRSdt = 0.0;
        double interpolated_dZCdt = 0.0;
        double interpolated_dZSdt = 0.0;

        double interpolated_lagrange_multiplierC = 0.0;
        double interpolated_lagrange_multiplierS = 0.0;

        double base_lagrange_multiplier = 0.0;

        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          const double psif_ = psif(l);
          double dpsifds_ = dpsifds(l, 0);

          interpolated_UC += this->u(l, 0) * psif_;
          interpolated_US += this->u(l, 1) * psif_;
          interpolated_VC += this->u(l, 2) * psif_;
          interpolated_VS += this->u(l, 3) * psif_;
          interpolated_WC += this->u(l, 4) * psif_;
          interpolated_WS += this->u(l, 5) * psif_;
          for (unsigned k = 0; k < 6; k++)
          {
            interpolated_dUds[k] += this->u(l, k) * dpsifds_;
          }

          interpolated_RC += this->Xhat(l, 0, 0) * psif_;
          interpolated_RS += this->Xhat(l, 0, 1) * psif_;
          interpolated_ZC += this->Xhat(l, 1, 0) * psif_;
          interpolated_ZS += this->Xhat(l, 1, 1) * psif_;

          interpolated_dRCds += this->Xhat(l, 0, 0) * dpsifds_;
          interpolated_dRSds += this->Xhat(l, 0, 1) * dpsifds_;
          interpolated_dZCds += this->Xhat(l, 1, 0) * dpsifds_;
          interpolated_dZSds += this->Xhat(l, 1, 1) * dpsifds_;

          interpolated_dRCdt += this->dXhatdt(l, 0, 0) * psif_;
          interpolated_dRSdt += this->dXhatdt(l, 0, 1) * psif_;
          interpolated_dZCdt += this->dXhatdt(l, 1, 0) * psif_;
          interpolated_dZSdt += this->dXhatdt(l, 1, 1) * psif_;

          // Calculate the tangent vector
          for (unsigned i = 0; i < 2; i++)
          {
            interpolated_t1[i] += this->nodal_position(l, i) * dpsifds_;
          }

          interpolated_lagrange_multiplierC +=
            nodal_value(l, lagrange_index(l, 0)) * psif_;
          interpolated_lagrange_multiplierS +=
            nodal_value(l, lagrange_index(l, 1)) * psif_;
        }

        Vector<double> x(3, 0.0);
        this->interpolated_x(s, x);
        this->get_base_lagrange_multiplier(
          node_pt(0)->time_stepper_pt()->time(),
          iplot,
          x,
          base_lagrange_multiplier);

        Vector<Vector<double>> tang_vec(1);
        tang_vec[0].resize(2, 0.0);
        this->continuous_tangent_and_outer_unit_normal(
          s, tang_vec, unit_normal);
        continuous_t1 = tang_vec[0];

        // Compare tangent vectors
        bool is_element_reversed = false;
        is_element_reversed = (continuous_t1[0] * interpolated_t1[0] +
                               continuous_t1[1] * interpolated_t1[1]) < 0;

        const double tlength = interpolated_t1[0] * interpolated_t1[0] +
                               interpolated_t1[1] * interpolated_t1[1];

        // Set the Jacobian of the line element
        const double J = sqrt(tlength);

        if (is_element_reversed)
        {
          interpolated_dRCds = -interpolated_dRCds;
          interpolated_dRSds = -interpolated_dRSds;
          interpolated_dZCds = -interpolated_dZCds;
          interpolated_dZSds = -interpolated_dZSds;

          for (unsigned k = 0; k < 6; k++)
          {
            interpolated_dUds[k] = -interpolated_dUds[k];
          }
        }

        // Output the x,y,u,v
        for (unsigned i = 0; i < n_dim; i++)
          outfile << this->interpolated_x(s, i) << " ";

        outfile << this->get_base_u(s, 0) << " ";
        outfile << this->get_base_u(s, 1) << " ";

        outfile << interpolated_UC << " ";
        outfile << interpolated_US << " ";
        outfile << interpolated_VC << " ";
        outfile << interpolated_VS << " ";
        outfile << interpolated_WC << " ";
        outfile << interpolated_WS << " ";

        outfile << interpolated_RC << " ";
        outfile << interpolated_RS << " ";
        outfile << interpolated_ZC << " ";
        outfile << interpolated_ZS << " ";

        outfile << interpolated_lagrange_multiplierC << " ";
        outfile << interpolated_lagrange_multiplierS << " ";

        outfile << unit_normal[0] << " ";
        outfile << unit_normal[1];

        if (has_output_derivatives)
        {
          outfile << " ";
          outfile << interpolated_dRCds / J << " ";
          outfile << interpolated_dRSds / J << " ";
          outfile << interpolated_dZCds / J << " ";
          outfile << interpolated_dZSds / J << " ";

          outfile << interpolated_dRCdt << " ";
          outfile << interpolated_dRSdt << " ";
          outfile << interpolated_dZCdt << " ";
          outfile << interpolated_dZSdt << " ";

          outfile << base_lagrange_multiplier;
        }
        outfile << " ";
        outfile << continuous_t1[0] << " ";
        outfile << continuous_t1[1] << " ";

        for (unsigned k = 0; k < 6; k++)
        {
          outfile << interpolated_dUds[k] << " ";
        }

        outfile << std::endl;
      }
      // write_tecplot_zone_footer(outfile, n_plot);
    }

    void output(std::ostream& outfile)
    {
      output(outfile, 3);
    }
  };
} // namespace oomph

#endif
