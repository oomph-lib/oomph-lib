// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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

#include "implicit_midpoint_rule.h"
#include "problem.h"

// Needed for mipoint update...
#include "mesh.h"
#include "elements.h"
#include "nodes.h"

namespace oomph
{
  /// \short This function advances the Data's time history so that
  /// we can move on to the next timestep
  void IMRBase::shift_time_values(Data* const& data_pt)
  {
    // Loop over the values, set previous values to the previous value, if
    // not a copy.
    for (unsigned j = 0, nj = data_pt->nvalue(); j < nj; j++)
    {
      if (!data_pt->is_a_copy(j))
      {
        for (unsigned t = ndt(); t > 0; t--)
        {
          data_pt->set_value(t, j, data_pt->value(t - 1, j));
        }
      }
    }
  }

  ///\short This function advances the time history of the positions
  /// at a node. ??ds Untested: I have no problems with moving nodes.
  void IMRBase::shift_time_positions(Node* const& node_pt)
  {
    // Find the number of coordinates
    unsigned n_dim = node_pt->ndim();
    // Find the number of position types
    unsigned n_position_type = node_pt->nposition_type();

    // Find number of stored timesteps
    unsigned n_tstorage = ntstorage();

    // Storage for the velocity
    double velocity[n_position_type][n_dim];

    // If adaptive, find the velocities
    if (adaptive_flag())
    {
      // Loop over the variables
      for (unsigned i = 0; i < n_dim; i++)
      {
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Initialise velocity to zero
          velocity[k][i] = 0.0;
          // Loop over all history data
          for (unsigned t = 0; t < n_tstorage; t++)
          {
            velocity[k][i] += Weight(1, t) * node_pt->x_gen(t, k, i);
          }
        }
      }
    }

    // Loop over the positions
    for (unsigned i = 0; i < n_dim; i++)
    {
      // If the position is not a copy
      if (node_pt->position_is_a_copy(i) == false)
      {
        // Loop over the position types
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Loop over stored times, and set values to previous values
          for (unsigned t = ndt(); t > 0; t--)
          {
            node_pt->x_gen(t, k, i) = node_pt->x_gen(t - 1, k, i);
          }
        }
      }
    }
  }


  /// Dummy - just check that the values that
  /// problem::calculate_predicted_values() has been called right.
  void IMRBase::calculate_predicted_values(Data* const& data_pt)
  {
    if (adaptive_flag())
    {
      // Can't do it here, but we can check that the predicted values have
      // been updated.
      check_predicted_values_up_to_date();
    }
  }


  double IMRBase::temporal_error_in_value(Data* const& data_pt,
                                          const unsigned& i)
  {
    if (adaptive_flag())
    {
      // predicted value is more accurate so just compare with that
      return data_pt->value(i) - data_pt->value(Predictor_storage_index, i);
    }
    else
    {
      std::string err("Tried to get the temporal error from a non-adaptive");
      err += " time stepper.";
      throw OomphLibError(
        err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  /// Half the timestep before starting solve
  void IMRByBDF::actions_before_timestep(Problem* problem_pt)
  {
    // Check that this is the only time stepper
#ifdef PARANOID
    if (problem_pt->ntime_stepper() != 1)
    {
      std::string err = "This implementation of midpoint can only work with a ";
      err += "single time stepper, try using IMR instead (but check ";
      err += "your Jacobian and residual calculations very carefully for "
             "compatability).";
      throw OomphLibError(
        err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
    }
#endif

    time_pt()->dt() /= 2;
    time_pt()->time() -= time_pt()->dt();

    // Set the weights again because we just changed the step size
    set_weights();


    if (problem_pt->use_predictor_values_as_initial_guess())
    {
      // Shift the initial guess to the midpoint so that it is an initial
      // guess for the newton step that we actually take.

      // Optimisation possiblity: here we update all values three time
      // (initial prediction, copy into initial guess, interpolate to
      // midpoint), could probably avoid this with more fancy code if
      // needed.

      // Get dofs at previous time step
      DoubleVector dof_n;
      problem_pt->get_dofs(1, dof_n);

      // Update dofs at current step to be the average of current and previous
      for (unsigned i = 0; i < problem_pt->ndof(); i++)
      {
        problem_pt->dof(i) = (problem_pt->dof(i) + dof_n[i]) / 2;
      }
    }
  }

  /// Local (not exported in header) helper function to handle midpoint
  /// update on a data object.
  void post_midpoint_update(Data* dat_pt, const bool& update_pinned)
  {
    if (!dat_pt->is_a_copy())
    {
      for (unsigned j = 0, nj = dat_pt->nvalue(); j < nj; j++)
      {
        int eqn = dat_pt->eqn_number(j);
        if (update_pinned || eqn >= 0)
        {
          double ynp1 = 2 * dat_pt->value(0, j) - dat_pt->value(1, j);
          dat_pt->set_value(0, j, ynp1);
        }
      }
    }
  }

  /// Take problem from t={n+1/2} to t=n+1 by algebraic update and restore
  /// time step.
  void IMRByBDF::actions_after_timestep(Problem* problem_pt)
  {
#ifdef PARANOID
    // Do it as dofs too to compare
    const unsigned ndof = problem_pt->ndof();
    DoubleVector dof_n, dof_np1;
    problem_pt->get_dofs(1, dof_n);
    problem_pt->get_dofs(dof_np1);

    for (unsigned i = 0; i < ndof; i++)
    {
      dof_np1[i] += dof_np1[i] - dof_n[i];
    }
#endif

    // First deal with global data
    for (unsigned i = 0, ni = problem_pt->nglobal_data(); i < ni; i++)
    {
      post_midpoint_update(problem_pt->global_data_pt(i), Update_pinned);
    }

    // Next element internal data
    for (unsigned i = 0, ni = problem_pt->mesh_pt()->nelement(); i < ni; i++)
    {
      GeneralisedElement* ele_pt = problem_pt->mesh_pt()->element_pt(i);
      for (unsigned j = 0, nj = ele_pt->ninternal_data(); j < nj; j++)
      {
        post_midpoint_update(ele_pt->internal_data_pt(j), Update_pinned);
      }
    }

    // Now the nodes
    for (unsigned i = 0, ni = problem_pt->mesh_pt()->nnode(); i < ni; i++)
    {
      post_midpoint_update(problem_pt->mesh_pt()->node_pt(i), Update_pinned);
    }

    // Update time
    problem_pt->time_pt()->time() += problem_pt->time_pt()->dt();

    // Return step size to its full value
    problem_pt->time_pt()->dt() *= 2;

#ifdef PARANOID
    using namespace StringConversion;
    DoubleVector actual_dof_np1;
    problem_pt->get_dofs(actual_dof_np1);

    for (unsigned j = 0; j < actual_dof_np1.nrow(); j++)
    {
      if (std::abs(actual_dof_np1[j] - dof_np1[j]) > 1e-8)
      {
        std::string err = "Got different values doing midpoint update via "
                          "extracted dofs than doing it in place!";
        err += to_string(actual_dof_np1[j]) + " vs " + to_string(dof_np1[j]);
        throw OomphLibError(
          err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
      }
    }

#endif
  }

} // namespace oomph
