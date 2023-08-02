#include "vorticity_smoother.h"

namespace oomph
{
  namespace VorticityRecoveryHelpers
  {
    /// Constructor
    RecoveryHelper::RecoveryHelper()
    {
      // Set the default (max) order of vorticity derivative to recover
      Maximum_order_of_vorticity_derivative = -1;

      // Set the default (max) order of velocity derivative to recover
      Maximum_order_of_velocity_derivative = 0;
    }

    /// The maximum order of derivatives calculated in the vorticity recovery
    int RecoveryHelper::maximum_order_of_vorticity_derivative() const
    {
      // Return the appropriate value
      return Maximum_order_of_vorticity_derivative;
    } // End of get_maximum_order_of_vorticity_derivative

    /// The maximum order of derivatives calculated in the velocity recovery
    int RecoveryHelper::maximum_order_of_velocity_derivative() const
    {
      // Return the appropriate value
      return Maximum_order_of_velocity_derivative;
    } // End of get_maximum_order_of_velocity_derivative

    /// The maximum order of derivatives calculated in the vorticity recovery
    void RecoveryHelper::set_maximum_order_of_vorticity_derivative(
      const int& max_deriv)
    {
      // Make sure the user has supplied a valid input value
      if ((max_deriv < -1) || (max_deriv > 3))
      {
        // Throw an error
        throw OomphLibError("Invalid input value! Should be between -1 and 3!",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Return the appropriate value
      Maximum_order_of_vorticity_derivative = max_deriv;

      // Calculate the value of Number_of_values_per_field with the updated
      // value of Maximum_order_of_vorticity_derivative
      calculate_number_of_values_per_field();
    } // End of set_maximum_order_of_vorticity_derivative

    /// The maximum order of derivatives calculated in the velocity recovery
    void RecoveryHelper::set_maximum_order_of_velocity_derivative(
      const int& max_deriv)
    {
      // Make sure the user has supplied a valid input value. Note, unlike the
      // vorticity, we always output the zeroth derivative of the velocity
      // so we don't use -1 as an input
      if ((max_deriv < 0) || (max_deriv > 1))
      {
        // Throw an error
        throw OomphLibError("Invalid input value! Should be between 0 and 3!",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Return the appropriate value
      Maximum_order_of_velocity_derivative = max_deriv;

      // Calculate the value of Number_of_values_per_field with the updated
      // value of Maximum_order_of_velocity_derivative
      calculate_number_of_values_per_field();
    } // End of set_maximum_order_of_vorticity_derivative

    /// Calculates the number of values per field given the number of
    /// vorticity and velocity derivatives to recover (stored as private data)
    void RecoveryHelper::calculate_number_of_values_per_field()
    {
      // Output: u,v,p
      Number_of_values_per_field = 3;

      // Loop over the vorticity derivatives
      for (unsigned i = 0;
           i < unsigned(Maximum_order_of_vorticity_derivative + 1);
           i++)
      {
        // Update the number of values per field
        Number_of_values_per_field += npartial_derivative(i);
      }

      // Loop over the velocity derivatives
      for (unsigned i = 1;
           i < unsigned(Maximum_order_of_velocity_derivative + 1);
           i++)
      {
        // Update the number of values per field
        Number_of_values_per_field += 2 * npartial_derivative(i);
      }
    } // End of calculate_number_of_values_per_field

    /// Helper function that determines the number of n-th order
    /// partial derivatives in d-dimensions. Specifically there are
    /// (n+d-1)(choose)(d-1) possible n-th order partial derivatives in
    /// d-dimensions. Implementation makes use of the code found at:
    ///    www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
    unsigned RecoveryHelper::npartial_derivative(const unsigned& n) const
    {
      // This will only work in 2D so n_dim is always 2
      unsigned n_dim = 2;

      // Calculate m
      unsigned n_bins = n + n_dim - 1;

      // Calculate k
      unsigned k = n_dim - 1;

      // Initialise the result
      unsigned value = 1;

      // Since C(n_bins,k)=C(n_bins,n_bins-k)
      if (k > n_bins - k)
      {
        // Replace k
        k = n_bins - k;
      }

      // Calculate [n_bins*(n_bins-1)*...*(n_bins-k+1)]/[k*(k-1)*...*1]
      for (unsigned i = 0; i < k; ++i)
      {
        // First update
        value *= (n_bins - i);

        // Second update
        value /= (i + 1);
      }

      // Return the result
      return value;
    } // End of npartial_derivative

    /// Number of continuously interpolated values:
    unsigned RecoveryHelper::ncont_interpolated_values() const
    {
      // Return the number of values used per field
      return Number_of_values_per_field;
    } // End of ncont_interpolated_values

    // Create an instance of the RecoveryHelper
    RecoveryHelper Recovery_helper;
  } // namespace VorticityRecoveryHelpers
} // namespace oomph
