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
#ifndef OOMPH_TRAPEZOID_RULE_H
#define OOMPH_TRAPEZOID_RULE_H

#include "problem.h"
#include "timesteppers.h"

namespace oomph
{
  /// Trapezoid rule time stepping scheme.
  ///
  /// This method requires a value of dy/dt at the initial time. The
  /// implementation of this calculation is exactly the same as is used for
  /// explicit time stepping.
  ///
  /// The function setup_initial_derivative(Problem* problem_pt) should be
  /// called after the initial conditions have been set, but before beginning
  /// time stepping, to compute this initial value of dy/dt.
  ///
  /// Warning: moving nodes not implemented (I have no test case).
  class TR : public TimeStepper
  {
    // The standard trapezoid rule is:

    // y_{n+1} = y_n + dt_n (f(t_n, y_n) + f(t_{n+1}, y_{n+1}))/2

    // where y is the vector of nodal values and f is the derivative function
    // (as in standard ODE theory). However we want to calculate time steps
    // using only one function evaluation (i.e. f) per time step because
    // function evaluations are expensive. So we require f(t_0, y_0) as
    // initial input and at each step store f(t_{n+1}, y_{n+1}) as a history
    // value for use in the calculatio of the next time step.


  public:
    /// Constructor, storage for two history derivatives (one for TR and
    /// one for the predictor step), one history value, present value and
    /// predicted value.
    TR(const bool& adaptive = false) : TimeStepper(2 + 2 + 1, 1)
    {
      // Set the weight for the zero-th derivative
      Weight(0, 0) = 1.0;

      Initial_derivative_set = false;

      Adaptive_Flag = adaptive;

      // Initialise adaptivity stuff
      Predictor_weight.assign(4, 0.0);
      Error_weight = 0.0;

      Shift_f = true;
    }

    /// Virtual destructor
    virtual ~TR() {}

    /// Return the actual order of the scheme
    unsigned order() const
    {
      return 2;
    }

    /// Set the weights
    void set_weights()
    {
      double dt = Time_pt->dt(0);
      Weight(1, 0) = 2.0 / dt;
      Weight(1, 1) = -2.0 / dt;
      Weight(1, derivative_index(0)) = -1.0; // old derivative
    }

    void set_error_weights()
    {
      double dt = Time_pt->dt(0);
      double dtprev = Time_pt->dt(1);
      Error_weight = 1 / (3 * (1 + (dtprev / dt)));
    }

    /// Function to set the predictor weights
    void set_predictor_weights()
    {
      // Read the value of the time steps
      double dt = Time_pt->dt(0);
      double dtprev = Time_pt->dt(1);
      double dtr = dt / dtprev;

      // y weights
      Predictor_weight[0] = 0.0;
      Predictor_weight[1] = 1.0;

      // dy weights
      Predictor_weight[derivative_index(0)] = (dt / 2) * (2 + dtr);
      Predictor_weight[derivative_index(1)] = -(dt / 2) * dtr;
    }

    /// Number of previous values available.
    unsigned nprev_values() const
    {
      return 1;
    }

    /// Location in data of derivatives
    unsigned derivative_index(const unsigned& t) const
    {
#ifdef PARANOID
      if (t >= 2)
      {
        std::string err = "Only storing two derivatives.";
        throw OomphLibError(
          err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
      }
#endif
      return t + 2;
    }

    /// Location of predicted value
    unsigned predicted_value_index() const
    {
      return derivative_index(1) + 1;
    }

    /// Number of timestep increments that need to be stored by the scheme
    unsigned ndt() const
    {
      return 2;
    }

    /// Initialise the time-history for the Data values, corresponding
    /// to an impulsive start.
    void assign_initial_values_impulsive(Data* const& data_pt)
    {
      throw OomphLibError("Function not yet implemented",
                          OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
    }

    /// Initialise the time-history for the nodal positions
    /// corresponding to an impulsive start.
    void assign_initial_positions_impulsive(Node* const& node_pt)
    {
      throw OomphLibError("Function not yet implemented",
                          OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
    }

    void actions_after_timestep(Problem* problem_pt)
    {
      // only ever want this to be true for one step
      Shift_f = true;
    }

    void actions_before_timestep(Problem* problem_pt)
    {
#ifdef PARANOID
      if (!Initial_derivative_set)
      {
        std::string err = "Initial derivative of TR has not been set";
        throw OomphLibError(
          err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    void setup_initial_derivative(Problem* problem_pt)
    {
      oomph_info
        << "Solving for derivative at initial time."
        << " Warning: if residual is not in the correct form this may fail."
        << std::endl;

      // Get the derivative at initial time
      DoubleVector f0;
      problem_pt->get_dvaluesdt(f0);

      // store in derivatives slot ready for use in timestepping.
      problem_pt->set_dofs(this->derivative_index(0), f0);

      Initial_derivative_set = true;
      Shift_f = false;
    }


    /// This function updates the Data's time history so that
    /// we can advance to the next timestep.
    void shift_time_values(Data* const& data_pt)
    {
      // Find number of values stored
      unsigned n_value = data_pt->nvalue();

      // Get previous step dt, time_pt has already been shifted so it's in
      // slot 1.
      double dtn = time_pt()->dt(1);

      // Loop over the values
      for (unsigned j = 0; j < n_value; j++)
      {
        // Set previous values to the previous value, if not a copy
        if (data_pt->is_a_copy(j) == false)
        {
          if (Shift_f)
          {
            // Calculate time derivative at step n by rearranging the TR
            // formula.
            double fnm1 = data_pt->value(derivative_index(0), j);
            double ynm1 = data_pt->value(1, j);
            double yn = data_pt->value(0, j);
            double fn = (2 / dtn) * (yn - ynm1) - fnm1;

            data_pt->set_value(derivative_index(0), j, fn);

            // Shift time derivative
            data_pt->set_value(derivative_index(1), j, fnm1);
          }
          else
          {
            std::cout << "didn't shift derivatives" << std::endl;
          }

          // Shift value
          data_pt->set_value(1, j, data_pt->value(0, j));
        }
      }
    }


    bool Initial_derivative_set;

    bool Shift_f;

    /// This function advances the time history of the positions at a
    /// node.
    void shift_time_positions(Node* const& node_pt)
    {
      // do nothing, not implemented for moving nodes
    }


    /// Function to calculate predicted positions at a node
    void calculate_predicted_positions(Node* const& node_pt)
    {
      // Loop over the dimensions
      unsigned n_dim = node_pt->ndim();
      for (unsigned j = 0; j < n_dim; j++)
      {
        // If the node is not copied
        if (!node_pt->position_is_a_copy(j))
        {
          // Initialise the predictor to zero
          double predicted_value = 0.0;

          // Loop over all the stored data and add appropriate values
          // to the predictor
          for (unsigned i = 1; i < 4; i++)
          {
            predicted_value += node_pt->x(i, j) * Predictor_weight[i];
          }

          // Store the predicted value
          node_pt->x(predicted_value_index(), j) = predicted_value;
        }
      }
    }

    /// Function to calculate predicted data values in a Data object
    void calculate_predicted_values(Data* const& data_pt)
    {
      // Loop over the values
      unsigned n_value = data_pt->nvalue();
      for (unsigned j = 0; j < n_value; j++)
      {
        // If the value is not copied
        if (!data_pt->is_a_copy(j))
        {
          // Loop over all the stored data and add appropriate
          // values to the predictor
          double predicted_value = 0.0;
          for (unsigned i = 1; i < 4; i++)
          {
            predicted_value += data_pt->value(i, j) * Predictor_weight[i];
          }

          // Store the predicted value
          data_pt->set_value(predicted_value_index(), j, predicted_value);
        }
      }
    }


    /// Compute the error in the position i at a node
    double temporal_error_in_position(Node* const& node_pt, const unsigned& i)
    {
      return Error_weight *
             (node_pt->x(i) - node_pt->x(predicted_value_index(), i));
    }

    /// Compute the error in the value i in a Data structure
    double temporal_error_in_value(Data* const& data_pt, const unsigned& i)
    {
      return Error_weight *
             (data_pt->value(i) - data_pt->value(predicted_value_index(), i));
    }


  private:
    /// Private data for the predictor weights
    Vector<double> Predictor_weight;

    /// Private data for the error weight
    double Error_weight;

    /// Broken copy constructor
    TR(const TR& dummy) = delete;

    /// Broken assignment operator
    void operator=(const TR& dummy) = delete;
  };

} // namespace oomph

#endif
