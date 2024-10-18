#ifndef HEIGHT_ELEMENT_HEADER
#define HEIGHT_ELEMENT_HEADER

#include <iostream>
#include "generic.h"

namespace oomph
{
  class HeightElement : public GeneralisedElement
  {
  private:
    int Traded_data_number;

  public:
    HeightElement(SolidNode* const& inner_node, SolidNode* const& outer_node)
      : Traded_data_number(-1)
    {
      add_external_data(inner_node->variable_position_pt());
      add_external_data(outer_node->variable_position_pt());
      // Storage for the height
      add_internal_data(new Data(1));
    }

    void set_parameter_data_pt(Data* const& parameter_data_pt)
    {
      // Only if not in external data
      if (Traded_data_number < 0)
      {
        Traded_data_number = this->add_external_data(parameter_data_pt);
      }
      else
      {
        throw("Can't set height step parameter a second time at the moment.");
      }
    }

    void pin_height()
    {
      this->internal_data_pt(0)->pin(0);
    }

    void unpin_height()
    {
      this->internal_data_pt(0)->unpin(0);
    }

    virtual void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      fill_in_generic_contribution_to_residuals_and_jacobian(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    virtual void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    // virtual void fill_in_contribution_to_jacobian(Vector<double>& residuals,
    //                                               DenseMatrix<double>&
    //                                               jacobian)
    //{
    //   fill_in_generic_contribution_to_residuals_and_jacobian(residuals,
    //                                                          jacobian);
    // }

    void fill_in_generic_contribution_to_residuals_and_jacobian(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      int local_eqn = this->internal_local_eqn(0, 0);
      // If the height is pinned then we trade for the external data.
      if (local_eqn < 0)
      {
        local_eqn = this->external_local_eqn(Traded_data_number, 0);
      }

      if (local_eqn >= 0)
      {
        residuals[local_eqn] += -this->internal_data_pt(0)->value(0) +
                                this->external_data_pt(1)->value(1) -
                                this->external_data_pt(0)->value(1);
      }
    }

    void step_height(const double& ds)
    {
      this->internal_data_pt(0)->set_value(
        0, this->internal_data_pt(0)->value(0) + ds);
    }

    void output(std::ostream& outfile)
    {
      outfile << this->internal_data_pt(0)->value(0) << std::endl;
    }
  };

} // namespace oomph

#endif
