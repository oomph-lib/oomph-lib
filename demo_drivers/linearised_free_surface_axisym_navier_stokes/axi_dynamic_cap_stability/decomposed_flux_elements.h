#ifndef DECOMPOSED_FLUX_ELEMENTS_HEADER
#define DECOMPOSED_FLUX_ELEMENTS_HEADER

namespace oomph
{
  template<class ELEMENT>
  class DecomposedFluxElements : public virtual FaceElement,
                                 public virtual FaceGeometry<ELEMENT>
  {
    friend class LinearisedElasticAxisymmetricFluidInterfaceElement<ELEMENT>;

  private:
    unsigned Data_number_of_lagrange_multiplier;
    unsigned Index_of_flux;

  public:
    DecomposedFluxElements(FiniteElement* const& element_pt,
                           const int& face_index,
                           Data* lagrange_data_pt,
                           unsigned data_index)
      : FaceElement(), FaceGeometry<ELEMENT>()
    {
      // Attach the geometrical information to the element, by
      // making the face element from the bulk element
      element_pt->build_face_element(face_index, this);

      // Setup data pointer
      Data_number_of_lagrange_multiplier = add_external_data(lagrange_data_pt);
      Index_of_flux = data_index;
    }

    void fill_in_generic_contribution_to_residuals_decomposed_integral(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      // Add in the volume constraint term if required
      const int local_eqn = this->external_local_eqn(
        Data_number_of_lagrange_multiplier, Index_of_flux);
      if (local_eqn >= 0)
      {
        const unsigned wc_index = 2;

        // Find out how many nodes there are
        const unsigned n_node = this->nnode();

        // Set up memeory for the shape functions
        Shape psif(n_node);
        DShape dpsifds(n_node, 1);

        // Set the value of n_intpt
        const unsigned n_intpt = this->integral_pt()->nweight();

        // Storage for the local cooridinate
        Vector<double> s(1);

        // Get a pointer to the parent element
        ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

        const double lagrange_multiplier =
          this->external_data_pt(Data_number_of_lagrange_multiplier)
            ->value(Index_of_flux);

        // Loop over the integration points
        for (unsigned ipt = 0; ipt < n_intpt; ipt++)
        {
          // Get the local coordinate at the integration point
          s[0] = this->integral_pt()->knot(ipt, 0);

          // Get the integral weight
          double W = this->integral_pt()->weight(ipt);

          // Call the derivatives of the shape function at the knot point
          this->dshape_local_at_knot(ipt, psif, dpsifds);

          // Get position and tangent vector
          Vector<double> interpolated_t1(2, 0.0);
          Vector<double> interpolated_x(2, 0.0);
          Vector<double> interpolated_vertical_velocity(1, 0.0);
          for (unsigned l = 0; l < n_node; l++)
          {
            // Loop over directional components
            for (unsigned i = 0; i < 2; i++)
            {
              interpolated_x[i] += this->nodal_position(l, i) * psif(l);
              interpolated_t1[i] += this->nodal_position(l, i) * dpsifds(l, 0);
            }
            // Cosine component of the vertical velocity
            interpolated_vertical_velocity[0] +=
              this->nodal_value(l, bulk_el_pt->u_index_lin_axi_nst(wc_index)) *
              psif(l);
          }

          // Calculate the length of the tangent Vector
          double tlength = interpolated_t1[0] * interpolated_t1[0] +
                           interpolated_t1[1] * interpolated_t1[1];

          // Set the Jacobian of the line element
          // multiplied by r (x[0])
          double J = sqrt(tlength) * interpolated_x[0];

          for (unsigned l = 0; l < n_node; l++)
          {
            // Add to residual with sign chosen so that the volume is
            // positive when the elements bound the fluid
            residuals[local_eqn] +=
              interpolated_vertical_velocity[0] * psif(l) * W * J;

            int local_eqn_2 =
              nodal_local_eqn(l, bulk_el_pt->u_index_lin_axi_nst(wc_index));
            if (local_eqn_2 >= 0)
            {
              residuals[local_eqn_2] += lagrange_multiplier * psif(l) * W * J;
            }
          }
        }
      }
    }

    /// Fill in the residuals for the volume constraint
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      this->fill_in_generic_contribution_to_residuals_decomposed_integral(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    // /// Fill in the residuals and jacobian for the volume constraint
    // void fill_in_contribution_to_jacobian(Vector<double>& residuals,
    //                                       DenseMatrix<double>& jacobian)
    // {
    //   this->fill_in_generic_contribution_to_residuals_decomposed_integral(
    //     residuals, jacobian, 1);

    //   FiniteElement::fill_in_contribution_to_jacobian_bf
    // }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    void output(std::ostream& outfile)
    {
      const unsigned nplot = 5;
      const unsigned n_node = nnode();
      // Get a pointer to the parent element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());
      const unsigned wc_index = 4 + 2;

      // Loop over plot points
      for (unsigned iplot = 0; iplot < nplot; iplot++)
      {
        // Get local coordinates of plot point
        Vector<double> s(1);
        get_s_plot(iplot, nplot, s);
        Shape psi(n_node);
        this->shape(s, psi);

        Vector<double> interpolated_vertical_velocity(1, 0.0);
        for (unsigned l = 0; l < n_node; l++)
        {
          // Cosine component of the vertical velocity
          interpolated_vertical_velocity[0] +=
            this->nodal_value(l, bulk_el_pt->u_index_lin_axi_nst(wc_index)) *
            psi(l);
        }

        // Output global coordinates to file
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << interpolated_x(s, i) << " ";
        }

        outfile << interpolated_vertical_velocity[0] << std::endl;
      }
    }
  };

} // namespace oomph

#endif
