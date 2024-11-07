#ifndef FULL_DECOMPOSED_INTEGRAL_ELEMENTS_HEADER
#define FULL_DECOMPOSED_INTEGRAL_ELEMENTS_HEADER

namespace oomph
{
  template<class ELEMENT>
  class FullDecomposedIntegralElements : public virtual FaceElement,
                                     public virtual FaceGeometry<ELEMENT>
  {
    friend class LinearisedElasticAxisymmetricFluidInterfaceElement<ELEMENT>;

  private:
    unsigned Data_number_of_integrand;
    unsigned Index_of_integrand;

  public:
    FullDecomposedIntegralElements(FiniteElement* const& element_pt,
                               const int& face_index,
                               Data* data_pt,
                               unsigned data_index)
      : FaceElement(), FaceGeometry<ELEMENT>()
    {
      // Attach the geometrical information to the element, by
      // making the face element from the bulk element
      element_pt->build_face_element(face_index, this);

      // Setup data pointer
      Data_number_of_integrand = add_external_data(data_pt);
      Index_of_integrand = data_index;
    }

    // Displacement values: n in {0,...,nnode()}, i in {R, Z}, j in {C, S}
    virtual double Xhat(const unsigned& n, const unsigned& i, const unsigned& j)
    {
      return this->nodal_value(n, u_index_linear_elasticity(n, i, j));
    }

    /// Return the index at which the i-th unknown displacement
    /// component is stored. The default value, i, is appropriate for
    /// single-physics problems.
    virtual inline unsigned u_index_linear_elasticity(const unsigned& n,
                                                      const unsigned& i,
                                                      const unsigned& j) const
    {
      return i * 2 + j;
    }


    void fill_in_generic_contribution_to_residuals_decomposed_integral(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      // Add in the volume constraint term if required
      const int local_eqn =
        this->external_local_eqn(Data_number_of_integrand, Index_of_integrand);
      if (local_eqn >= 0)
      {
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
        OverlayingMyLinearElement* bulk_el_pt =
          dynamic_cast<OverlayingMyLinearElement*>(this->bulk_element_pt());

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
          Vector<double> continuous_t1(2, 0.0);
          Vector<double> interpolated_x(2, 0.0);

          double interpolated_RC = 0;
          double interpolated_RS = 0;
          double interpolated_ZC = 0;
          double interpolated_ZS = 0;

          double interpolated_dRCds = 0;
          double interpolated_dRSds = 0;
          double interpolated_dZCds = 0;
          double interpolated_dZSds = 0;

          for (unsigned l = 0; l < n_node; l++)
          {
            const double psif_ = psif(l);
            double dpsifds_ = dpsifds(l, 0);

            // Loop over directional components
            for (unsigned i = 0; i < 2; i++)
            {
              interpolated_x[i] += this->nodal_position(l, i) * psif(l);
              interpolated_t1[i] += this->nodal_position(l, i) * dpsifds(l, 0);
            }


            interpolated_RC += this->Xhat(l, 0, 0) * psif_;
            interpolated_RS += this->Xhat(l, 0, 1) * psif_;
            interpolated_ZC += this->Xhat(l, 1, 0) * psif_;
            interpolated_ZS += this->Xhat(l, 1, 1) * psif_;

            interpolated_dRCds += this->Xhat(l, 0, 0) * dpsifds_;
            interpolated_dRSds += this->Xhat(l, 0, 1) * dpsifds_;
            interpolated_dZCds += this->Xhat(l, 1, 0) * dpsifds_;
            interpolated_dZSds += this->Xhat(l, 1, 1) * dpsifds_;
          }

          Vector<double> unit_normal(2, 0.0);
          Vector<Vector<double>> tang_vec(1);
          tang_vec[0].resize(2, 0.0);
          continuous_tangent_and_outer_unit_normal(s, tang_vec, unit_normal);
          continuous_t1 = tang_vec[0];

          bool is_element_reversed = false;
          is_element_reversed = (continuous_t1[0] * interpolated_t1[0] +
                                 continuous_t1[1] * interpolated_t1[1]) < 0;

          if (is_element_reversed)
          {
            interpolated_dRCds = -interpolated_dRCds;
            interpolated_dRSds = -interpolated_dRSds;
            interpolated_dZCds = -interpolated_dZCds;
            interpolated_dZSds = -interpolated_dZSds;
          }


          // Calculate the length of the tangent Vector
          double tlength = interpolated_t1[0] * interpolated_t1[0] +
                           interpolated_t1[1] * interpolated_t1[1];

          // Set the Jacobian of the line element
          // multiplied by r (x[0])
          double J = sqrt(tlength) * interpolated_x[0];


          // Add to residual with sign chosen so that the volume is
          // positive when the elements bound the fluid
          residuals[local_eqn] += (continuous_t1[1] * interpolated_RC * J +
                                   interpolated_x[0] * interpolated_dZCds) *
                                    interpolated_x[0] * W +
                                  (continuous_t1[0] * interpolated_RC * J +
                                   interpolated_x[0] * interpolated_dRCds) *
                                    interpolated_x[1] * W +
                                  (-continuous_t1[1] * interpolated_RC +
                                   continuous_t1[0] * interpolated_ZC) *
                                    J * W;
        }
      }
    }

    /// Fill in the residuals for the volume constraint
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      this->fill_in_generic_contribution_to_residuals_decomposed_integral(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    // /// Fill in the residuals and jacobian for the volume constraint
    // void fill_in_contribution_to_jacobian(Vector<double>& residuals,
    //                                       DenseMatrix<double>& jacobian)
    // {
    //   this->fill_in_generic_contribution_to_residuals_decomposed_integral(
    //     residuals, jacobian, 1);

    //   FiniteElement::fill_in_contribution_to_jacobian_bf
    // }
  };

  class IntegralElement : public GeneralisedElement
  {
  private:
    unsigned Data_number;

  public:
    IntegralElement() : GeneralisedElement()
    {
      Data_number = this->add_internal_data(new Data(1), true);
    }

    Data* const& get_data_pt()
    {
      return this->internal_data_pt(Data_number);
    }


    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      fill_in_contribution_to_generic_residuals(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      fill_in_contribution_to_generic_residuals(residuals, jacobian, 1);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      fill_in_contribution_to_generic_residuals(residuals, jacobian, 1);
    }

    void fill_in_contribution_to_generic_residuals(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      unsigned value_index = 0;
      const int local_eqn = this->internal_local_eqn(Data_number, value_index);
      if (local_eqn >= 0)
      {
        residuals[local_eqn] -=
          this->internal_data_pt(Data_number)->value(value_index);

        if (flag)
        {
          const int local_unknown =
            this->internal_local_eqn(Data_number, value_index);
          jacobian(local_eqn, local_unknown) -= 1;
        }
      }
    }
  };

  class DataTradingElement : public GeneralisedElement
  {
  private:
    unsigned Value_data_number;
    unsigned Traded_data_number;

  public:
    DataTradingElement(Data* const& value_data_pt, Data* const& traded_data_pt)
    {
      Value_data_number = this->add_external_data(value_data_pt, true);
      Traded_data_number = this->add_external_data(traded_data_pt, true);
    }

    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      fill_in_contribution_to_generic_residuals(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      fill_in_contribution_to_generic_residuals(residuals, jacobian, 1);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      fill_in_contribution_to_generic_residuals(residuals, jacobian, 1);
    }

    void fill_in_contribution_to_generic_residuals(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      unsigned value_index = 0;
      const int local_eqn =
        this->external_local_eqn(Traded_data_number, value_index);
      if (local_eqn >= 0)
      {
        residuals[local_eqn] +=
          this->external_data_pt(Value_data_number)->value(value_index);
        if (flag)
        {
          const int local_unknown =
            this->external_local_eqn(Value_data_number, value_index);
          jacobian(local_eqn, local_unknown) += 1;
        }
      }
    }
  };

} // namespace oomph

#endif

