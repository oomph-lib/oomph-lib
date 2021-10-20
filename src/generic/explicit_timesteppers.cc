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
// Non-inline member functions for the explicit timesteppers
#include "explicit_timesteppers.h"
#include "timesteppers.h"

namespace oomph
{
  //===================================================================
  /// Dummy value of time always set to zero
  //==================================================================
  double ExplicitTimeSteppableObject::Dummy_time_value = 0.0;

  //==================================================================
  /// A single virtual function that returns the residuals
  /// vector multiplied by the inverse mass matrix
  //=================================================================
  void ExplicitTimeSteppableObject::get_dvaluesdt(DoubleVector& minv_res)
  {
    std::ostringstream error_stream;
    error_stream
      << "Empty default function called.\n"
      << "The function must return the solution x of the linear system\n"
      << "                    M x = R\n"
      << "in order for the object to be used by an ExplicitTimeStepper.\n"
      << "NOTE: It is the responsibility of the object to set the size \n"
      << "      of the vector x\n";


    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  //=======================================================================
  /// Function that should get the values of the dofs in the object
  //=======================================================================
  void ExplicitTimeSteppableObject::get_dofs(DoubleVector& dofs) const
  {
    std::ostringstream error_stream;
    error_stream
      << "Empty default function called.\n"
      << "The function must return the current values of the degrees of \n"
      << "freedom in the object.\n"
      << "Note: It is the responsibility of the object to set the size\n"
      << "of the vector\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  //=======================================================================
  /// Function that should get the values of the dofs in the object
  //=======================================================================
  void ExplicitTimeSteppableObject::get_dofs(const unsigned& t,
                                             DoubleVector& dofs) const
  {
    std::ostringstream error_stream;
    error_stream
      << "Empty default function called.\n"
      << "The function must return the t'th history values of the degrees of \n"
      << "freedom in the object.\n"
      << "Note: It is the responsibility of the object to set the size\n"
      << "of the vector\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  //=====================================================================
  /// Function that sets the values of the dofs in the object
  //====================================================================
  void ExplicitTimeSteppableObject::set_dofs(const DoubleVector& dofs)
  {
    std::ostringstream error_stream;
    error_stream
      << "Empty default function called.\n"
      << "The function must set the current values of the degrees of \n"
      << "freedom in the object.\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


  //====================================================================
  /// Function that adds the lambda multiplied by the increment_dofs
  /// vector to the dofs in the object
  //====================================================================
  void ExplicitTimeSteppableObject::add_to_dofs(
    const double& lambda, const DoubleVector& increment_dofs)
  {
    std::ostringstream error_stream;
    error_stream
      << "Empty default function called.\n"
      << "The function must add lambda multiplied by the contents of the\n"
      << "input vector to the degrees of freedom in the object.\n"
      << "Note: It is the responsibility of the object to ensure that the\n"
      << "      the degrees of freedom are in the same order as those \n"
      << "      returned by get_dvaluesdt()\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  //==================================================================
  /// Virtual function that should be overloaded to return access
  /// to the local time in the object
  //=================================================================
  double& ExplicitTimeSteppableObject::time()
  {
    std::ostringstream error_stream;
    error_stream << "Empty default function called.\n"
                 << "The function must return a reference to the local time in "
                    "the object\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);

    return Dummy_time_value;
  }

  ///  Virtual function that should be overloaded to return a pointer to a
  /// Time object.
  Time* ExplicitTimeSteppableObject::time_pt() const
  {
    std::ostringstream error_stream;
    error_stream
      << "Empty default function called.\n"
      << "The function must return a pointer to an oomph-lib Time object.\n";
    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);

    return 0;
  }


  //================================================================
  /// Euler timestepping x^{t+1} = x^{t} + dt M^{-1} L(x^{t})
  //=================================================================
  void Euler::timestep(ExplicitTimeSteppableObject* const& object_pt,
                       const double& dt)
  {
    object_pt->actions_before_explicit_timestep();
    object_pt->actions_before_explicit_stage();

    // Vector that will hold the inverse mass matrix multiplied by the
    // residuals
    DoubleVector minv_res;
    // Get M^{-1} R
    object_pt->get_dvaluesdt(minv_res);

    // Add the result to the unknowns
    object_pt->add_to_dofs(dt, minv_res);
    // Increment the time by the appropriate amount
    object_pt->time() += dt;

    // Call any actions required after the change in the unknowns
    object_pt->actions_after_explicit_stage();
    object_pt->actions_after_explicit_timestep();
  }


  //====================================================================
  // Broken default timestep function for RungeKutta schemes
  //====================================================================
  template<unsigned ORDER>
  void RungeKutta<ORDER>::timestep(
    ExplicitTimeSteppableObject* const& object_pt, const double& dt)
  {
    std::ostringstream error_stream;
    error_stream << "Timestep not implemented for order " << ORDER << "\n";
    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  //===================================================================
  /// Explicit specialisation for fourth-order RK scheme
  //==================================================================
  template<>
  void RungeKutta<4>::timestep(ExplicitTimeSteppableObject* const& object_pt,
                               const double& dt)
  {
    object_pt->actions_before_explicit_timestep();

    // Stage 1
    // ============================================================
    object_pt->actions_before_explicit_stage();

    // Store the initial values and initial time
    DoubleVector u;
    object_pt->get_dofs(u);

    // Now get the first unknowns
    DoubleVector k1;
    object_pt->get_dvaluesdt(k1);

    // Add to the residuals
    object_pt->add_to_dofs(0.5 * dt, k1);
    // Increment the time
    object_pt->time() += 0.5 * dt;
    object_pt->actions_after_explicit_stage();


    // Stage 2
    // ============================================================
    object_pt->actions_before_explicit_stage();

    // Get the next unknowns
    DoubleVector k2;
    object_pt->get_dvaluesdt(k2);

    // Now reset the residuals
    object_pt->set_dofs(u);
    object_pt->add_to_dofs(0.5 * dt, k2);
    // Time remains the same
    object_pt->actions_after_explicit_stage();

    // Stage 3
    // ============================================================
    object_pt->actions_before_explicit_stage();

    // Get the next unknowns
    DoubleVector k3;
    object_pt->get_dvaluesdt(k3);

    // Now reset the residuals
    object_pt->set_dofs(u);
    object_pt->add_to_dofs(dt, k3);
    // Increment the time (now at initial_time + dt)
    object_pt->time() += 0.5 * dt;
    object_pt->actions_after_explicit_stage();

    // Stage 4
    // ============================================================
    object_pt->actions_before_explicit_stage();

    // Get the final unknowns
    DoubleVector k4;
    object_pt->get_dvaluesdt(k4);

    // Set the final values of the unknowns
    object_pt->set_dofs(u);
    object_pt->add_to_dofs((dt / 6.0), k1);
    object_pt->add_to_dofs((dt / 3.0), k2);
    object_pt->add_to_dofs((dt / 3.0), k3);
    object_pt->add_to_dofs((dt / 6.0), k4);
    object_pt->actions_after_explicit_stage();

    object_pt->actions_after_explicit_timestep();
  }

  //===================================================================
  /// Explicit specialisation for second-order RK scheme
  //==================================================================
  template<>
  void RungeKutta<2>::timestep(ExplicitTimeSteppableObject* const& object_pt,
                               const double& dt)
  {
    object_pt->actions_before_explicit_timestep();

    // Stage 1
    // ============================================================
    object_pt->actions_before_explicit_stage();

    // Store the initial values
    DoubleVector u;
    object_pt->get_dofs(u);

    // Get f1 (time derivative at t0, y0) and add to dofs
    DoubleVector f1;
    object_pt->get_dvaluesdt(f1);
    object_pt->add_to_dofs(dt, f1);

    // Advance time to t1 = t0 + dt
    object_pt->time() += dt;

    object_pt->actions_after_explicit_stage();


    // Stage 2
    // ============================================================
    object_pt->actions_before_explicit_stage();

    // get f2 (with t=t1, y = y0 + h f0)
    DoubleVector f2;
    object_pt->get_dvaluesdt(f2);

    // Final answer is starting dofs + h/2 * (f1 + f2)
    object_pt->set_dofs(u);
    object_pt->add_to_dofs(0.5 * dt, f1);
    object_pt->add_to_dofs(0.5 * dt, f2);

    object_pt->actions_after_explicit_stage();

    // Done, do actions
    object_pt->actions_after_explicit_timestep();
  }


  //=================================================================
  // General constructor for LowOrder RK schemes
  //==================================================================
  template<unsigned ORDER>
  LowStorageRungeKutta<ORDER>::LowStorageRungeKutta()
  {
    Type = "LowStorageRungeKutta";
  }

  //======================================================================
  /// Broken default timestep for LowStorageRungeKutta
  //======================================================================
  template<unsigned ORDER>
  void LowStorageRungeKutta<ORDER>::timestep(
    ExplicitTimeSteppableObject* const& object_pt, const double& dt)
  {
    std::ostringstream error_stream;
    error_stream << "Timestep not implemented for order " << ORDER << "\n";
    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  //================================================================
  // Specialised constructor for fourth-order RK scheme
  //================================================================
  template<>
  LowStorageRungeKutta<4>::LowStorageRungeKutta()
  {
    Type = "LowStorageRungeKutta";

    A.resize(5);
    A[0] = 0.0;
    A[1] = -567301805773.0 / 1357537059087.0;
    A[2] = -2404267990393.0 / 2016746695238.0;
    A[3] = -3550918686646.0 / 2091501179385.0;
    A[4] = -1275806237668.0 / 842570457699.0;

    B.resize(5);
    B[0] = 1432997174477.0 / 9575080441755.0;
    B[1] = 5161836677717.0 / 13612068292357.0;
    B[2] = 1720146321549.0 / 2090206949498.0;
    B[3] = 3134564353537.0 / 4481467310338.0;
    B[4] = 2277821191437.0 / 14882151754819.0;

    C.resize(5);
    C[0] = B[0];
    C[1] = 2526269341429.0 / 6820363962896.0;
    C[2] = 2006345519317.0 / 3224310063776.0;
    C[3] = 2802321613138.0 / 2924317926251.0;
    C[4] = 1.0;
  }


  // Explicit specialisation for fourth-order RK scheme
  template<>
  void LowStorageRungeKutta<4>::timestep(
    ExplicitTimeSteppableObject* const& object_pt, const double& dt)
  {
    object_pt->actions_before_explicit_timestep();

    // Store the initial time
    const double initial_time = object_pt->time();
    // Temporary storage
    DoubleVector k;

    // Temporary storage for the inverse mass matrix multiplied by the residuals
    DoubleVector minv_res;

    // Loop over the number of steps in the scheme
    for (unsigned i = 0; i < 5; i++)
    {
      object_pt->actions_before_explicit_stage();

      // Get the inverse mass matrix multiplied by the current value
      // of the residuals
      object_pt->get_dvaluesdt(minv_res);
      // Get the values of k
      const unsigned n_dof = minv_res.nrow();

      // First time round resize k and initialise to zero
      if (i == 0)
      {
        k.build(minv_res.distribution_pt(), 0.0);
      }
      // Now construct the next value of k
      for (unsigned n = 0; n < n_dof; n++)
      {
        k[n] *= A[i];
        k[n] += dt * minv_res[n];
      }

      // Add to the residuals
      object_pt->add_to_dofs(B[i], k);
      // Set the new time
      object_pt->time() = initial_time + C[i] * dt;

      object_pt->actions_after_explicit_stage();
    }

    object_pt->actions_after_explicit_timestep();
  }

  // ??ds this could be heavily optimised if needed. Keeping it simple for
  // now
  void EBDF3::timestep(ExplicitTimeSteppableObject* const& object_pt,
                       const double& dt)
  {
    using namespace StringConversion;

    object_pt->actions_before_explicit_timestep();
    object_pt->actions_before_explicit_stage();

    // Storage indicies for the history values that we need
    unsigned tn = 1;
    unsigned tnm1 = tn + 1;
    unsigned tnm2 = tnm1 + 1;

    // Check dts are the same, this will need to be removed if ebdf3 is being
    // used as something other than a predictor... But seeing as it isn't
    // stable that isn't likely.
#ifdef PARANOID
    if (std::abs(dt - object_pt->time_pt()->dt(0)) > 1e-15)
    {
      std::string err =
        "Inconsistent dts! Predictor is stepping by " + to_string(dt);
      err += " but dt(0) = " + to_string(object_pt->time_pt()->dt(0));
      throw OomphLibError(
        err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
    }
#endif

    // Get older dt values
    double dtnm1 = object_pt->time_pt()->dt(1);
    double dtnm2 = object_pt->time_pt()->dt(2);

    // Calculate weights for these dts
    set_weights(dt, dtnm1, dtnm2);

    // Get derivative value at step n (even though this uses values from t=0 =
    // step n+1, it's ok because we haven't changed the values in that slot
    // yet).
    DoubleVector fn;
    object_pt->get_dvaluesdt(fn);
    fn *= Fn_weight;

    // Extract history values and multiply by their weights
    DoubleVector ynp1, yn, ynm1, ynm2;
    object_pt->get_dofs(tn, yn);
    yn *= Yn_weight;
    object_pt->get_dofs(tnm1, ynm1);
    ynm1 *= Ynm1_weight;
    object_pt->get_dofs(tnm2, ynm2);
    ynm2 *= Ynm2_weight;


    // Add everything together
    ynp1 = yn;
    ynp1 += ynm1;
    ynp1 += ynm2;
    ynp1 += fn;

    // Done, update things in the object
    object_pt->set_dofs(ynp1);
    object_pt->time() += dt;

    // Actions functions
    object_pt->actions_after_explicit_stage();
    object_pt->actions_after_explicit_timestep();
  }

  /// Calculate the weights for this set of step sizes.
  void EBDF3::set_weights(const double& dtn,
                          const double& dtnm1,
                          const double& dtnm2)
  {
    using namespace std;

    // If this is slow we can probably optimise by doing direct
    // multiplication instead of using pow.

    // Generated using sympy from my python ode code.

    double denom = pow(dtnm1, 4) * dtnm2 + 2 * pow(dtnm1, 3) * pow(dtnm2, 2) +
                   pow(dtnm1, 2) * pow(dtnm2, 3);

    Yn_weight =
      -(2 * pow(dtn, 3) * dtnm1 * dtnm2 + pow(dtn, 3) * pow(dtnm2, 2) +
        3 * pow(dtn, 2) * pow(dtnm1, 2) * dtnm2 +
        3 * pow(dtn, 2) * dtnm1 * pow(dtnm2, 2) + pow(dtn, 2) * pow(dtnm2, 3) -
        pow(dtnm1, 4) * dtnm2 - 2 * pow(dtnm1, 3) * pow(dtnm2, 2) -
        pow(dtnm1, 2) * pow(dtnm2, 3)) /
      denom;

    Ynm1_weight =
      -(-pow(dtn, 3) * pow(dtnm1, 2) - 2 * pow(dtn, 3) * dtnm1 * dtnm2 -
        pow(dtn, 3) * pow(dtnm2, 2) - pow(dtn, 2) * pow(dtnm1, 3) -
        3 * pow(dtn, 2) * pow(dtnm1, 2) * dtnm2 -
        3 * pow(dtn, 2) * dtnm1 * pow(dtnm2, 2) - pow(dtn, 2) * pow(dtnm2, 3)) /
      denom;

    Ynm2_weight =
      -(pow(dtn, 3) * pow(dtnm1, 2) + pow(dtn, 2) * pow(dtnm1, 3)) / denom;

    Fn_weight =
      -(-pow(dtn, 3) * pow(dtnm1, 2) * dtnm2 -
        pow(dtn, 3) * dtnm1 * pow(dtnm2, 2) -
        2 * pow(dtn, 2) * pow(dtnm1, 3) * dtnm2 -
        3 * pow(dtn, 2) * pow(dtnm1, 2) * pow(dtnm2, 2) -
        pow(dtn, 2) * dtnm1 * pow(dtnm2, 3) - dtn * pow(dtnm1, 4) * dtnm2 -
        2 * dtn * pow(dtnm1, 3) * pow(dtnm2, 2) -
        dtn * pow(dtnm1, 2) * pow(dtnm2, 3)) /
      denom;
  }


  // Force build of templates
  template class RungeKutta<4>;
  template class LowStorageRungeKutta<4>;
} // namespace oomph
