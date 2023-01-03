// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_IMPLICIT_MIDPOINT_RULE_H
#define OOMPH_IMPLICIT_MIDPOINT_RULE_H

// oomph-lib headers
#include "Vector.h"
#include "nodes.h"
#include "matrices.h"
#include "timesteppers.h"
#include "explicit_timesteppers.h"

namespace oomph
{
  // Forward decl. so that we can have function of Problem*
  class Problem;


  /// Implicit midpoint rule base class for the two implementations.
  class IMRBase : public TimeStepper
  {
  public:
    /// Constructor with initialisation
    IMRBase(const bool& adaptive = false)
      : TimeStepper(2, 1) // initialise weights later
    {
      Adaptive_Flag = adaptive;
      Is_steady = false;
      Type = "Midpoint method";
      Predict_by_explicit_step = true;

      // If adaptive then we are storing predicted values in slot
      // 4. Otherwise we aren't storing them so leave it as -1.
      if (adaptive)
      {
        Predictor_storage_index = 4;
      }


      // Storage for weights needs to be 2x(2 + 0/2) (the +2 is extra history
      // for ebdf3 predictor if adaptive). This means we provide ways to
      // calculate the zeroth and first derivatives and in calculations we
      // use 2 + 0/2 time values.

      // But actually (for some reason...) the amount of data storage
      // allocated for the timestepper in each node is determined by the
      // size of weight. So we need to store an extra dummy entry in order
      // to have storage space for the predicted value at t_{n+1}.

      // Just leave space for predictor values anyway, it's not expensive.
      Weight.resize(2, 5, 0.0);

      // Assume that set_weights() or make_steady() is called before use to
      // set up the values stored in Weight.
    }

    /// Destructor
    virtual ~IMRBase() {}

    /// Setup weights for time derivative calculations.
    virtual void set_weights() = 0;

    /// Number of history values to interpolate over to get the "current"
    /// value.
    virtual unsigned nprev_values_for_value_at_evaluation_time() const = 0;

    /// Actual order (accuracy) of the scheme
    unsigned order() const
    {
      return 2;
    }

    /// Number of timestep increments that are required by the scheme
    unsigned ndt() const
    {
      return nprev_values();
    }

    /// ??ds
    unsigned nprev_values() const
    {
      return 4;
    }

    /// This function advances the Data's time history so that
    /// we can move on to the next timestep
    void shift_time_values(Data* const& data_pt);

    /// This function advances the time history of the positions
    /// at a node.
    void shift_time_positions(Node* const& node_pt);

    /// Set the weights for the error computation. This is not used
    /// by midpoint rule.
    void set_error_weights() {}

    /// Set the weights for the predictor previous timestep. This is not
    /// used by midpint rule.
    void set_predictor_weights() {}

    /// not implemented (??ds TODO)
    void assign_initial_values_impulsive(Data* const& data_pt)
    {
      std::string err = "Not implemented";
      throw OomphLibError(
        err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    void assign_initial_positions_impulsive(Node* const& node_pt)
    {
      std::string err = "Not implemented";
      throw OomphLibError(
        err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    void calculate_predicted_positions(Node* const& node_pt)
    {
      std::string err = "Not implemented";
      throw OomphLibError(
        err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    double temporal_error_in_position(Node* const& node_pt, const unsigned& i)
    {
      std::string err = "Not implemented";
      throw OomphLibError(
        err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Adaptivity
    void calculate_predicted_values(Data* const& data_pt);
    double temporal_error_in_value(Data* const& data_pt, const unsigned& i);
  };

  /// The "real" implementation of the implicit midpoint rule. Implemented
  /// by calculation of residuals etc. at half step. This requires
  /// non-trivial modifications to the element's residual and Jacobian
  /// calculation functions to interpolate values to the midpoint. As such
  /// IMRByBDF should be preferred.
  ///
  /// However this class must be used when multiple different time steppers
  /// are being used simultaneously for different parts of the problem.
  class IMR : public IMRBase
  {
  public:
    /// Common mistakes when using this implementation of midpoint:
    /// * Didn't include the d_u_evaltime_by_du_np1 term in the Jacobian.
    /// * Didn't properly interpolate time/values/x/derivatives in
    ///   jacobian/residual (all of these MUST be evaluated at the midpoint).


    /// Constructor with initialisation
    IMR(const bool& adaptive = false) : IMRBase(adaptive) {}

    /// Destructor, predictor_pt handled by base
    virtual ~IMR() {}

    /// Setup weights for time derivative calculations.
    void set_weights()
    {
      // Set the weights for zero-th derivative (i.e. the value to use in
      // newton solver calculations, implicit midpoint method uses the
      // average of previous and current values).
      Weight(0, 0) = 0.5;
      Weight(0, 1) = 0.5;

      // Set the weights for 1st time derivative
      double dt = Time_pt->dt(0);
      Weight(1, 0) = 1.0 / dt;
      Weight(1, 1) = -1.0 / dt;
    }

    /// Number of history values to interpolate over to get the
    /// "current" value.
    unsigned nprev_values_for_value_at_evaluation_time() const
    {
      return 2;
    }

  private:
    /// Inaccessible copy constructor.
    IMR(const IMR& dummy) {}

    /// Inaccessible assignment operator.
    void operator=(const IMR& dummy) {}
  };


  /// Implementation of implicit midpoint rule by taking half a step of bdf1
  /// then applying an update to all dofs. This implementation *should* work
  /// with any existing problem for which the BDF methods work.
  ///
  /// The exception is when multiple different time steppers are being used
  /// simultaneously for different parts of the problem. In this case the
  /// IMR class should be used.
  class IMRByBDF : public IMRBase
  {
  public:
    /// Constructor with initialisation
    IMRByBDF(const bool& adaptive = false) : IMRBase(adaptive)
    {
      Update_pinned = true;
    }

    /// Destructor
    virtual ~IMRByBDF() {}

    /// Setup weights for time derivative calculations.
    void set_weights()
    {
      // Use weights from bdf1
      double dt = Time_pt->dt(0);
      Weight(1, 0) = 1.0 / dt;
      Weight(1, 1) = -1.0 / dt;
    }

    /// Number of history values to interpolate over to get the
    /// "current" value. Evaluation time is the end of the bdf1 "half-step",
    /// so only need one value as normal.
    unsigned nprev_values_for_value_at_evaluation_time() const
    {
      return 1;
    }

    /// Half the timestep before starting solve
    void actions_before_timestep(Problem* problem_pt);

    /// Take problem from t={n+1/2} to t=n+1 by algebraic update and restore
    /// time step.
    void actions_after_timestep(Problem* problem_pt);

    /// Should we update pinned variables after the half-step?
    bool Update_pinned;

  private:
    /// Inaccessible copy constructor.
    IMRByBDF(const IMRByBDF& dummy) {}

    /// Inaccessible assignment operator.
    void operator=(const IMRByBDF& dummy) {}
  };

} // namespace oomph

#endif
