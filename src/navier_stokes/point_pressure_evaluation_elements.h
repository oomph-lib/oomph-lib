#ifndef POINT_PRESSURE_EVALUATION_ELEMENTS_HEADER
#define POINT_PRESSURE_EVALUATION_ELEMENTS_HEADER

#include "generic.h"

namespace oomph
{
  //==================CLASS FOR THE PRESSURE CONTRIBUTION================
  /// This class adds the finite element pressure at the evaluation point
  /// to the residual for the singular eigensolution function.
  ///
  /// R_C += +- p_FE (Evaluation_point)
  ///
  /// and thus regularises the FE solution by matching the pressure at
  /// two locations. If the amplitude of the  singular solution is known,
  /// pin C.
  //=====================================================================
  class PointPressureEvaluationElement : public GeneralisedElement
  {
  private:
    // Storage for the bulk element
    int Pressure_index;
    int Pressure_data_index;
    int Scaling_index;
    int Scaling_data_index;
    bool Is_adding_to_residuals;
    Node* Node_pt;

  public:
    // Constructor
    PointPressureEvaluationElement(Node* const& node_pt,
                                   const unsigned& pressure_value_index)
      : GeneralisedElement(),
        Pressure_index(pressure_value_index),
        Pressure_data_index(-1),
        Scaling_index(0),
        Scaling_data_index(-1),
        Is_adding_to_residuals(true),
        Node_pt(node_pt)
    {
      Pressure_data_index = add_external_data(Node_pt);
    }

    void set_add_to_residuals()
    {
      Is_adding_to_residuals = true;
    }
    void set_subtract_from_residuals()
    {
      Is_adding_to_residuals = false;
    }

    // Set and add the pressure data as external data
    void set_pressure_data_pt(Data* const& scaling_data_pt)
    {
      Scaling_data_index = add_external_data(scaling_data_pt);
    }

    // Calculate the element's residual vector
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_pressure_contribution(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      this->fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    void fill_in_contribution_to_dresiduals_dparameter(
      double* const& parameter_pt, Vector<double>& dres_dparam)
    {
    }

  protected:
    // Generic residual and Jacobian routine
    void fill_in_generic_residual_contribution_pressure_contribution(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
    {
      // Set the residual multiplier, dependent on whether we are adding or
      // subtracting to the residuals
      double multiplier = 1.0;
      if (!Is_adding_to_residuals)
      {
        multiplier = -1.0;
      }

      // Set the local equation
      int local_eqn = 0;

      // Add to singular function scaling residual
      local_eqn = this->external_local_eqn(Scaling_data_index, Scaling_index);

      // If the equation is not pinned
      if (local_eqn >= 0)
      {
        // Add (or subtract) the pressure at the evaluation point
        residuals[local_eqn] += Node_pt->value(Pressure_index) * multiplier;

        // If the Jacobian flag is on, add to the Jacobian
        if (flag)
        {
          // Initialise a variable for the local_unknown
          int local_unknown = 0;

          // The residual depends on the pressure at each of the bulk
          // elements nodes, which are stored here as external data.
          local_unknown =
            this->external_local_eqn(Pressure_data_index, Pressure_index);

          // If not pinned
          if (local_unknown > 0)
          {
            // Add the contribution of the node to the local jacobian
            jacobian(local_eqn, local_unknown) += multiplier;
          }
        }
      }
    }

    // Overwrite the output function
    void output(std::ostream& outfile)
    {
      // Vector of local coordinates
      const unsigned n_dim = 2;

      // Spatial coordinates are one higher
      for (unsigned i = 0; i < n_dim + 1; i++)
      {
        outfile << Node_pt->x(i) << " ";
      }

      // Output the pressure
      outfile << Node_pt->value(Pressure_index) << " ";

      // End of line
      outfile << std::endl;
    }
  };
} // namespace oomph
#endif
