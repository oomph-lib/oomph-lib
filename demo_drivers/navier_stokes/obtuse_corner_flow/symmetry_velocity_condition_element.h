#ifndef SYMMETRY_VELOCITY_CONDITION_ELEMENT_HEADER
#define SYMMETRY_VELOCITY_CONDITION_ELEMENT_HEADER

namespace oomph
{
  class SymmetryVelocityConditionElement : public GeneralisedElement
  {
  public:
    SymmetryVelocityConditionElement(Node* const& node) : GeneralisedElement()
    {
      add_external_data(dynamic_cast<Data*>(node));
    }

    void fill_in_contribution_to_residuals(Vector<double>& residuals)

    {
      fill_in_contribution_to_residuals_symmetry_condition(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      GeneralisedElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    /// Fill in contribution from Jacobian
    // void fill_in_contribution_to_jacobian(Vector<double>& residuals,
    //                                       DenseMatrix<double>& jacobian)
    // {
    //   // Fill in analytic contribution of internal equations
    //   fill_in_contribution_to_residuals_symmetry_condition(
    //     residuals, jacobian, 1);
    // }

    void fill_in_contribution_to_residuals_symmetry_condition(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
    {
      const unsigned uc_index = 4;
      const unsigned us_index = 5;
      const unsigned vc_index = 8;
      const unsigned vs_index = 9;
      const double uc = this->external_data_pt(0)->value(uc_index);
      const double us = this->external_data_pt(0)->value(us_index);
      const double vc = this->external_data_pt(0)->value(vc_index);
      const double vs = this->external_data_pt(0)->value(vs_index);
      int local_eqn = this->external_local_eqn(0, uc_index);
      if (local_eqn >= 0)
      {
        residuals[local_eqn] += uc + vs;
      }
      local_eqn = this->external_local_eqn(0, us_index);
      if (local_eqn >= 0)
      {
        residuals[local_eqn] += us - vc;
      }
    }
  };

} // namespace oomph

#endif
