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
// Include the header
#include "timesteppers.h"

// Required so that we can delete the predictor.
#include "explicit_timesteppers.h"


namespace oomph
{
  /// Destructor for time stepper, in .cc file so that we don't have to
  /// include explicit_timesteppers.h in header.
  TimeStepper::~TimeStepper()
  {
    // Make sure explicit predictor gets deleted if it was set
    delete Explicit_predictor_pt;
    Explicit_predictor_pt = 0;
  }


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  // Static variables for steady timestepper
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////

  template<unsigned NSTEPS>
  double Steady<NSTEPS>::One = 1.0;

  template<unsigned NSTEPS>
  double Steady<NSTEPS>::Zero = 0.0;

  template<unsigned NSTEPS>
  Time Steady<NSTEPS>::Dummy_time(NSTEPS);

  // Force the instantiation for all NSTEPS values that may be
  // required in other timesteppers
  template class Steady<0>;
  template class Steady<1>;
  template class Steady<2>;
  template class Steady<3>;
  template class Steady<4>;


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  // Weights for specific bdf schemes
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Assign the values of the weights
  //=======================================================================
  template<>
  void BDF<1>::set_weights()
  {
    double dt = Time_pt->dt(0);
    Weight(1, 0) = 1.0 / dt;
    Weight(1, 1) = -1.0 / dt;
  }


  //======================================================================
  /// Set the predictor weights (Gresho and Sani pg. 270)
  //======================================================================
  template<>
  void BDF<1>::set_predictor_weights()
  {
    // If it's adaptive set the predictor weights
    if (adaptive_flag())
    {
      double dt = Time_pt->dt(0);

      // Set the predictor weights
      Predictor_weight[0] = 0.0;
      Predictor_weight[1] = 1.0;

      // Velocity weight
      Predictor_weight[2] = dt;
    }
  }


  //=======================================================================
  /// Calculate the predicted values and store them at the appropriate
  /// location in the data structure
  /// This function must be called after the time-values have been shifted!
  //=======================================================================
  template<>
  void BDF<1>::calculate_predicted_values(Data* const& data_pt)
  {
    // If it's adaptive calculate the values
    if (adaptive_flag())
    {
      // Find number of values
      unsigned n_value = data_pt->nvalue();
      // Loop over the values
      for (unsigned j = 0; j < n_value; j++)
      {
        // If the value is not copied
        if (data_pt->is_a_copy(j) == false)
        {
          // Now Initialise the predictor to zero
          double predicted_value = 0.0;
          // Now loop over all the stored data and add appropriate values
          // to the predictor
          for (unsigned i = 1; i < 3; i++)
          {
            predicted_value += data_pt->value(i, j) * Predictor_weight[i];
          }
          // Store the predicted value
          data_pt->set_value(Predictor_storage_index, j, predicted_value);
        }
      }
    }
  }


  //=======================================================================
  /// Calculate predictions for the positions
  //=======================================================================
  template<>
  void BDF<1>::calculate_predicted_positions(Node* const& node_pt)
  {
    // Only do this if adaptive
    if (adaptive_flag())
    {
      // Find number of dimensions of the problem
      unsigned n_dim = node_pt->ndim();
      // Loop over the dimensions
      for (unsigned j = 0; j < n_dim; j++)
      {
        // If the node is not copied
        if (node_pt->position_is_a_copy(j) == false)
        {
          // Initialise the predictor to zero
          double predicted_value = 0.0;
          // Now loop over all the stored data and add appropriate values
          // to the predictor
          for (unsigned i = 1; i < 3; i++)
          {
            predicted_value += node_pt->x(i, j) * Predictor_weight[i];
          }
          // Store the predicted value
          node_pt->x(Predictor_storage_index, j) = predicted_value;
        }
      }
    }
  }


  //=======================================================================
  /// Set the error weights  (Gresho and Sani pg. 270)
  //=======================================================================
  template<>
  void BDF<1>::set_error_weights()
  {
    Error_weight = 0.5;
  }

  //===================================================================
  /// Function to compute the error in position i at node
  //===================================================================
  template<>
  double BDF<1>::temporal_error_in_position(Node* const& node_pt,
                                            const unsigned& i)
  {
    if (adaptive_flag())
    {
      // Just return the error
      return Error_weight *
             (node_pt->x(i) - node_pt->x(Predictor_storage_index, i));
    }
    else
    {
      return 0.0;
    }
  }


  //=========================================================================
  /// Function to calculate the error in the data value i
  //=========================================================================
  template<>
  double BDF<1>::temporal_error_in_value(Data* const& data_pt,
                                         const unsigned& i)
  {
    if (adaptive_flag())
    {
      // Just return the error
      return Error_weight *
             (data_pt->value(i) - data_pt->value(Predictor_storage_index, i));
    }
    else
    {
      return 0.0;
    }
  }

  //=======================================================================
  /// Assign the values of the weights. The scheme used is from Gresho &
  /// Sani, Incompressible Flow and the Finite Element Method
  /// (advection-diffusion and isothermal laminar flow), 1998, pg. 715 (with
  /// some algebraic rearrangement).
  //=======================================================================
  template<>
  void BDF<2>::set_weights()
  {
    double dt = Time_pt->dt(0);
    double dtprev = Time_pt->dt(1);
    Weight(1, 0) = 1.0 / dt + 1.0 / (dt + dtprev);
    Weight(1, 1) = -(dt + dtprev) / (dt * dtprev);
    Weight(1, 2) = dt / ((dt + dtprev) * dtprev);

    if (adaptive_flag())
    {
      Weight(1, 3) = 0.0;
      Weight(1, 4) = 0.0;
    }
  }

  //========================================================================
  /// Calculate the predictor weights. The scheme used is from Gresho &
  /// Sani, Incompressible Flow and the Finite Element Method
  /// (advection-diffusion and isothermal laminar flow), 1998, pg. 715.
  //========================================================================
  template<>
  void BDF<2>::set_predictor_weights()
  {
    // If it's adaptive set the predictor weights
    if (adaptive_flag())
    {
      // Read the value of the previous timesteps
      double dt = Time_pt->dt(0);
      double dtprev = Time_pt->dt(1);

      // Set the predictor weights
      Predictor_weight[0] = 0.0;
      Predictor_weight[1] = 1.0 - (dt * dt) / (dtprev * dtprev);
      Predictor_weight[2] = (dt * dt) / (dtprev * dtprev);
      // Acceleration term
      Predictor_weight[3] = (1.0 + dt / dtprev) * dt;
    }
  }

  //=======================================================================
  /// Calculate the predicted values and store them at the appropriate
  /// location in the data structure
  /// This function must be called after the time-values have been shifted!
  //=======================================================================
  template<>
  void BDF<2>::calculate_predicted_values(Data* const& data_pt)
  {
    // If it's adaptive calculate the values
    if (adaptive_flag())
    {
      // Find number of values
      unsigned n_value = data_pt->nvalue();
      // Loop over the values
      for (unsigned j = 0; j < n_value; j++)
      {
        // If the value is not copied
        if (data_pt->is_a_copy(j) == false)
        {
          // Now Initialise the predictor to zero
          double predicted_value = 0.0;
          // Now loop over all the stored data and add appropriate values
          // to the predictor
          for (unsigned i = 1; i < 4; i++)
          {
            predicted_value += data_pt->value(i, j) * Predictor_weight[i];
          }
          // Store the predicted value
          // Note that predictor is stored as the FIFTH entry in this scheme
          data_pt->set_value(Predictor_storage_index, j, predicted_value);
        }
      }
    }
  }

  //=======================================================================
  /// Calculate predictions for the positions
  //=======================================================================
  template<>
  void BDF<2>::calculate_predicted_positions(Node* const& node_pt)
  {
    // Only do this if adaptive
    if (adaptive_flag())
    {
      // Find number of dimensions of the problem
      unsigned n_dim = node_pt->ndim();
      // Loop over the dimensions
      for (unsigned j = 0; j < n_dim; j++)
      {
        // If the node is not copied
        if (node_pt->position_is_a_copy(j) == false)
        {
          // Initialise the predictor to zero
          double predicted_value = 0.0;
          // Now loop over all the stored data and add appropriate values
          // to the predictor
          for (unsigned i = 1; i < 4; i++)
          {
            predicted_value += node_pt->x(i, j) * Predictor_weight[i];
          }
          // Store the predicted value
          node_pt->x(Predictor_storage_index, j) = predicted_value;
        }
      }
    }
  }


  //=======================================================================
  /// Function that sets the error weights
  //=======================================================================
  template<>
  void BDF<2>::set_error_weights()
  {
    if (adaptive_flag())
    {
      double dt = Time_pt->dt(0);
      double dtprev = Time_pt->dt(1);
      // Calculate the error weight
      Error_weight =
        pow((1.0 + dtprev / dt), 2.0) /
        (1.0 + 3.0 * (dtprev / dt) + 4.0 * pow((dtprev / dt), 2.0) +
         2.0 * pow((dtprev / dt), 3.0));
    }
  }


  //===================================================================
  /// Function to compute the error in position i at node
  //===================================================================
  template<>
  double BDF<2>::temporal_error_in_position(Node* const& node_pt,
                                            const unsigned& i)
  {
    if (adaptive_flag())
    {
      // Just return the error
      return Error_weight *
             (node_pt->x(i) - node_pt->x(Predictor_storage_index, i));
    }
    else
    {
      return 0.0;
    }
  }


  //=========================================================================
  /// Function to calculate the error in the data value i
  //=========================================================================
  template<>
  double BDF<2>::temporal_error_in_value(Data* const& data_pt,
                                         const unsigned& i)
  {
    if (adaptive_flag())
    {
      // Just return the error
      return Error_weight *
             (data_pt->value(i) - data_pt->value(Predictor_storage_index, i));
    }
    else
    {
      return 0.0;
    }
  }


  //====================================================================
  /// Assign the values of the weights; pass the value of the timestep
  //====================================================================
  template<>
  void BDF<4>::set_weights()
  {
#ifdef PARANOID
    double dt0 = Time_pt->dt(0);
    for (unsigned i = 0; i < Time_pt->ndt(); i++)
    {
      if (dt0 != Time_pt->dt(i))
      {
        throw OomphLibError("BDF4 currently only works for fixed timesteps \n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif
    double dt = Time_pt->dt(0);
    Weight(1, 0) = 25.0 / 12.0 / dt;
    Weight(1, 1) = -48.0 / 12.0 / dt;
    Weight(1, 2) = 36.0 / 12.0 / dt;
    Weight(1, 3) = -16.0 / 12.0 / dt;
    Weight(1, 4) = 3.0 / 12.0 / dt;
  }


  //======================================================================
  /// Calculate the predictor weights
  //======================================================================
  template<>
  void BDF<4>::set_predictor_weights()
  {
    // Read the value of the previous timesteps
    // double dt=Time_pt->dt(0);
    // double dtprev=Time_pt->dt(1);

    throw OomphLibError(
      "Not implemented yet", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  //=======================================================================
  /// Calculate the predicted values and store them at the appropriate
  /// location in the data structure
  /// This function must be called after the time-values have been shifted!
  //=======================================================================
  template<>
  void BDF<4>::calculate_predicted_values(Data* const& data_pt)
  {
    throw OomphLibError(
      "Not implemented yet", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  /// Calculate predictions for the positions
  template<>
  void BDF<4>::calculate_predicted_positions(Node* const& node_pt)
  {
    throw OomphLibError(
      "Not implemented yet", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  /// Function that sets the error weights
  template<>
  void BDF<4>::set_error_weights()
  {
    throw OomphLibError(
      "Not implemented yet", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  //===================================================================
  /// Function to compute the error in position i at node
  //===================================================================
  template<>
  double BDF<4>::temporal_error_in_position(Node* const& node_pt,
                                            const unsigned& i)
  {
    throw OomphLibError(
      "Not implemented yet", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    return 0.0;
  }


  //=========================================================================
  /// Function to calculate the error in the data value i
  //=========================================================================
  template<>
  double BDF<4>::temporal_error_in_value(Data* const& data_pt,
                                         const unsigned& i)
  {
    throw OomphLibError(
      "Not implemented yet", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    return 0.0;
  }


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //=========================================================================
  /// Initialise the time-history for the values,
  /// corresponding to an impulsive start.
  //=========================================================================
  template<unsigned NSTEPS>
  void Newmark<NSTEPS>::assign_initial_values_impulsive(Data* const& data_pt)
  {
    // Find number of values stored in Data object
    unsigned n_value = data_pt->nvalue();

    // Loop over values
    for (unsigned j = 0; j < n_value; j++)
    {
      // Set previous values to the initial value, if not a copy
      if (data_pt->is_a_copy(j) == false)
      {
        for (unsigned t = 1; t <= NSTEPS; t++)
        {
          data_pt->set_value(t, j, data_pt->value(j));
        }
      }

      // Initial velocity and initial acceleration are zero
      data_pt->set_value(NSTEPS + 1, j, 0.0);
      data_pt->set_value(NSTEPS + 2, j, 0.0);
    }
  }


  //=========================================================================
  /// Initialise the time-history for the values,
  /// corresponding to an impulsive start.
  //=========================================================================
  template<unsigned NSTEPS>
  void Newmark<NSTEPS>::assign_initial_positions_impulsive(Node* const& node_pt)
  {
    // Find the number of coordinates
    unsigned n_dim = node_pt->ndim();
    // Find the number of position types
    unsigned n_position_type = node_pt->nposition_type();

    // Loop over values
    for (unsigned i = 0; i < n_dim; i++)
    {
      // If not a copy
      if (node_pt->position_is_a_copy(i) == false)
      {
        // Loop over the position types
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Set previous value to the initial value, if not a copy
          for (unsigned t = 1; t <= NSTEPS; t++)
          {
            node_pt->x_gen(t, k, i) = node_pt->x_gen(k, i);
          }

          // Initial velocity and initial acceleration are zero
          node_pt->x_gen(NSTEPS + 1, k, i) = 0.0;
          node_pt->x_gen(NSTEPS + 2, k, i) = 0.0;
        }
      }
    }
  }


  //=========================================================================
  ///  Initialise the time-history for the Data values,
  /// so that the Newmark representations for current veloc and
  /// acceleration are exact.
  //=========================================================================
  template<unsigned NSTEPS>
  void Newmark<NSTEPS>::assign_initial_data_values(
    Data* const& data_pt,
    Vector<InitialConditionFctPt> initial_value_fct,
    Vector<InitialConditionFctPt> initial_veloc_fct,
    Vector<InitialConditionFctPt> initial_accel_fct)
  {
    // Set weights
    set_weights();

#ifdef PARANOID
    // Check if the 3 vectors of functions have the same size
    if (initial_value_fct.size() != initial_veloc_fct.size() ||
        initial_value_fct.size() != initial_accel_fct.size())
    {
      throw OomphLibError("The Vectors of fcts for values, velocities and "
                          "acceleration must be the same size",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // The number of functions in each vector (they should be the same)
    unsigned n_fcts = initial_value_fct.size();

#ifdef PARANOID

    // If there are more data values at the node than functions, issue a warning
    unsigned n_value = data_pt->nvalue();
    if (n_value > n_fcts && !Shut_up_in_assign_initial_data_values)
    {
      std::stringstream message;
      message << "There are more data values than initial value fcts.\n";
      message << "Only the first " << n_fcts << " data values will be set\n";
      OomphLibWarning(
        message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Loop over elements of the function vectors
    for (unsigned j = 0; j < n_fcts; j++)
    {
      if (initial_value_fct[j] == 0)
      {
#ifdef PARANOID
        if (!Shut_up_in_assign_initial_data_values)
        {
          std::stringstream message;
          message << "Ignoring value " << j << " in assignment of ic.\n";
          OomphLibWarning(message.str(),
                          "Newmark<NSTEPS>::assign_initial_data_values",
                          OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }
      else
      {
        // Value itself at current and previous times
        for (unsigned t = 0; t <= NSTEPS; t++)
        {
          double time_local = Time_pt->time(t);
          data_pt->set_value(t, j, initial_value_fct[j](time_local));
        }

        // Now, rather than assigning the values for veloc and accel
        // in the Newmark scheme directly, we solve a linear system
        // to determine the values required to make the Newmark
        // representations of the  veloc and accel at the current (!)
        // time are exact.

        // Initial time: The value itself
        double time_ = Time_pt->time();
        double U = initial_value_fct[j](time_);

        // Value itself at previous time
        time_ = Time_pt->time(1);
        double U0 = initial_value_fct[j](time_);

        // Veloc (time deriv) at present time -- this is what the
        // Newmark scheme should return!
        time_ = Time_pt->time(0);
        double Udot = initial_veloc_fct[j](time_);

        // Acccel (2nd time deriv) at present time -- this is what the
        // Newmark scheme should return!
        time_ = Time_pt->time(0);
        double Udotdot = initial_accel_fct[j](time_);

        Vector<double> vect(2);
        vect[0] = Udotdot - Weight(2, 0) * U - Weight(2, 1) * U0;
        vect[1] = Udot - Weight(1, 0) * U - Weight(1, 1) * U0;

        DenseDoubleMatrix matrix(2, 2);
        matrix(0, 0) = Weight(2, NSTEPS + 1);
        matrix(0, 1) = Weight(2, NSTEPS + 2);
        matrix(1, 0) = Weight(1, NSTEPS + 1);
        matrix(1, 1) = Weight(1, NSTEPS + 2);

        matrix.solve(vect);

        // Discrete veloc (time deriv) at previous time , to be entered into the
        // Newmark scheme  so that its prediction for the *current* veloc
        // and accel is correct.
        data_pt->set_value(NSTEPS + 1, j, vect[0]);

        // Discrete veloc (2nd time deriv) at previous time, to be entered into
        // the Newmark scheme  so that its prediction for the *current* veloc
        // and accel is correct.
        data_pt->set_value(NSTEPS + 2, j, vect[1]);
      }
    }
  }


  //=========================================================================
  ///  Initialise the time-history for the Data values,
  /// so that the Newmark representations for current veloc and
  /// acceleration are exact.
  //=========================================================================
  template<unsigned NSTEPS>
  void Newmark<NSTEPS>::assign_initial_data_values(
    Node* const& node_pt,
    Vector<NodeInitialConditionFctPt> initial_value_fct,
    Vector<NodeInitialConditionFctPt> initial_veloc_fct,
    Vector<NodeInitialConditionFctPt> initial_accel_fct)
  {
    // Set weights
    set_weights();


#ifdef PARANOID
    // Check if the 3 vectors of functions have the same size
    if (initial_value_fct.size() != initial_veloc_fct.size() ||
        initial_value_fct.size() != initial_accel_fct.size())
    {
      throw OomphLibError("The Vectors of fcts for values, velocities and "
                          "acceleration must be the same size",
                          "Newmark<NSTEPS>::assign_initial_data_values",
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // The number of functions in each vector (they should be the same)
    unsigned n_fcts = initial_value_fct.size();

#ifdef PARANOID

    // If there are more data values at the node than functions, issue a warning
    unsigned n_value = node_pt->nvalue();
    if (n_value > n_fcts && !Shut_up_in_assign_initial_data_values)
    {
      std::stringstream message;
      message << "There are more nodal values than initial value fcts.\n";
      message << "Only the first " << n_fcts << " nodal values will be set\n";
      OomphLibWarning(
        message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Get current nodal coordinates
    unsigned n_dim = node_pt->ndim();
    Vector<double> x(n_dim);
    for (unsigned i = 0; i < n_dim; i++) x[i] = node_pt->x(i);

    // Loop over values
    for (unsigned j = 0; j < n_fcts; j++)
    {
      if (initial_value_fct[j] == 0)
      {
#ifdef PARANOID
        if (!Shut_up_in_assign_initial_data_values)
        {
          std::stringstream message;
          message << "Ignoring value " << j << " in assignment of ic.\n";
          OomphLibWarning(
            message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }
      else
      {
        // Value itself at current and previous times
        for (unsigned t = 0; t <= NSTEPS; t++)
        {
          double time_local = Time_pt->time(t);
          node_pt->set_value(t, j, initial_value_fct[j](time_local, x));
        }

        // Now, rather than assigning the values for veloc and accel
        // in the Newmark scheme directly, we solve a linear system
        // to determine the values required to make the Newmark
        // representations of the  veloc and accel at the current (!)
        // time are exact.

        // Initial time: The value itself
        double time_ = Time_pt->time();
        double U = initial_value_fct[j](time_, x);

        // Value itself at previous time
        time_ = Time_pt->time(1);
        double U0 = initial_value_fct[j](time_, x);

        // Veloc (time deriv) at present time -- this is what the
        // Newmark scheme should return!
        time_ = Time_pt->time(0);
        double Udot = initial_veloc_fct[j](time_, x);

        // Acccel (2nd time deriv) at present time -- this is what the
        // Newmark scheme should return!
        time_ = Time_pt->time(0);
        double Udotdot = initial_accel_fct[j](time_, x);

        Vector<double> vect(2);
        vect[0] = Udotdot - Weight(2, 0) * U - Weight(2, 1) * U0;
        vect[1] = Udot - Weight(1, 0) * U - Weight(1, 1) * U0;

        DenseDoubleMatrix matrix(2, 2);
        matrix(0, 0) = Weight(2, NSTEPS + 1);
        matrix(0, 1) = Weight(2, NSTEPS + 2);
        matrix(1, 0) = Weight(1, NSTEPS + 1);
        matrix(1, 1) = Weight(1, NSTEPS + 2);

        matrix.solve(vect);

        // Discrete veloc (time deriv) at previous time , to be entered into the
        // Newmark scheme  so that its prediction for the *current* veloc
        // and accel is correct.
        node_pt->set_value(NSTEPS + 1, j, vect[0]);

        // Discrete veloc (2nd time deriv) at previous time, to be entered into
        // the Newmark scheme  so that its prediction for the *current* veloc
        // and accel is correct.
        node_pt->set_value(NSTEPS + 2, j, vect[1]);
      }
    }
  }


  //=========================================================================
  ///  First step in a two-stage procedure to assign
  /// the history values for the Newmark scheme so
  /// that the veloc and accel that are computed by the scheme
  /// are correct at the current time.
  ///
  /// Call this function for t_deriv=0,1,2,3.
  /// When calling with
  /// - t_deriv=0 :  data_pt(0,*) should contain the values at the
  ///                previous timestep
  /// - t_deriv=1 :  data_pt(0,*) should contain the current values;
  ///                they get stored (temporarily) in data_pt(1,*).
  /// - t_deriv=2 :  data_pt(0,*) should contain the current velocities
  ///                (first time derivs); they get stored (temporarily)
  ///                in data_pt(NSTEPS+1,*).
  /// - t_deriv=3 :  data_pt(0,*) should contain the current accelerations
  ///                (second time derivs); they get stored (temporarily)
  ///                in data_pt(NSTEPS+2,*).
  /// .
  /// Follow this by calls to
  /// \code
  /// assign_initial_data_values_stage2(...)
  /// \endcode
  //=========================================================================
  template<unsigned NSTEPS>
  void Newmark<NSTEPS>::assign_initial_data_values_stage1(
    const unsigned t_deriv, Data* const& data_pt)
  {
    // Find number of values stored
    unsigned n_value = data_pt->nvalue();

    // Loop over values
    for (unsigned j = 0; j < n_value; j++)
    {
      // Which deriv are we dealing with?
      switch (t_deriv)
      {
          // 0-th deriv: the value itself at present time
        case 0:

          data_pt->set_value(1, j, data_pt->value(j));
          break;

          // 1st deriv: the veloc at present time
        case 1:

          data_pt->set_value(NSTEPS + 1, j, data_pt->value(j));
          break;

          // 2nd deriv: the accel at present time
        case 2:

          data_pt->set_value(NSTEPS + 2, j, data_pt->value(j));
          break;


          // None other are possible!
        default:

          std::ostringstream error_message_stream;
          error_message_stream << "t_deriv " << t_deriv
                               << " is not possible with a Newmark scheme "
                               << std::endl;

          throw OomphLibError(error_message_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }
  }

  //=========================================================================
  /// Second step in a two-stage procedure to assign
  /// the history values for the Newmark scheme so
  /// that the veloc and accel that are computed by the scheme
  /// are correct at the current time.
  ///
  /// This assigns appropriate values for the "previous
  /// velocities and accelerations" so that their current
  /// values, which were defined in assign_initial_data_values_stage1(...),
  /// are represented exactly by the Newmark scheme.
  //=========================================================================
  template<unsigned NSTEPS>
  void Newmark<NSTEPS>::assign_initial_data_values_stage2(Data* const& data_pt)
  {
    // Find number of values stored
    unsigned n_value = data_pt->nvalue();

    // Loop over values
    for (unsigned j = 0; j < n_value; j++)
    {
      // Rather than assigning the previous veloc and accel
      // directly, solve linear system to determine their values
      // so that the Newmark representations are exact.

      // The exact value at present time (stored in stage 1)
      double U = data_pt->value(1, j);

      // Exact value at previous time
      double U0 = data_pt->value(0, j);

      // Veloc (1st time deriv) at present time (stored in stage 1)
      double Udot = data_pt->value(NSTEPS + 1, j);

      // Acccel (2nd time deriv) at present time (stored in stage 1)
      double Udotdot = data_pt->value(NSTEPS + 2, j);

      Vector<double> vect(2);
      vect[0] = Udotdot - Weight(2, 0) * U - Weight(2, 1) * U0;
      vect[1] = Udot - Weight(1, 0) * U - Weight(1, 1) * U0;

      DenseDoubleMatrix matrix(2, 2);
      matrix(0, 0) = Weight(2, NSTEPS + 1);
      matrix(0, 1) = Weight(2, NSTEPS + 2);
      matrix(1, 0) = Weight(1, NSTEPS + 1);
      matrix(1, 1) = Weight(1, NSTEPS + 2);

      matrix.solve(vect);


      // Now assign the discrete coefficients
      // so that the Newmark scheme returns the correct
      // results for the present values, and 1st and 2nd derivatives:

      // Value at present time
      data_pt->set_value(0, j, U);

      // Value at previous time
      data_pt->set_value(1, j, U0);

      // Veloc (time deriv) at previous time
      data_pt->set_value(NSTEPS + 1, j, vect[0]);

      // Acccel (2nd time deriv) at previous time
      data_pt->set_value(NSTEPS + 2, j, vect[1]);
    }
  }


  //=========================================================================
  /// This function updates the Data's time history so that
  /// we can advance to the next timestep.
  //=========================================================================
  template<unsigned NSTEPS>
  void Newmark<NSTEPS>::shift_time_values(Data* const& data_pt)
  {
    // Find number of values stored
    unsigned n_value = data_pt->nvalue();

    Vector<double> veloc(n_value);
    time_derivative(1, data_pt, veloc);

    Vector<double> accel(n_value);
    time_derivative(2, data_pt, accel);

    // Loop over values
    for (unsigned j = 0; j < n_value; j++)
    {
      // Set previous values/veloc/accel to present values/veloc/accel,
      // if not a copy
      if (data_pt->is_a_copy(j) == false)
      {
        for (unsigned t = NSTEPS; t > 0; t--)
        {
          data_pt->set_value(t, j, data_pt->value(t - 1, j));
        }
        data_pt->set_value(NSTEPS + 1, j, veloc[j]);
        data_pt->set_value(NSTEPS + 2, j, accel[j]);
      }
    }
  }

  //=========================================================================
  /// This function updates a nodal time history so that
  /// we can advance to the next timestep.
  //=========================================================================
  template<unsigned NSTEPS>
  void Newmark<NSTEPS>::shift_time_positions(Node* const& node_pt)
  {
    // Find the number of coordinates
    unsigned n_dim = node_pt->ndim();
    // Find the number of position types
    unsigned n_position_type = node_pt->nposition_type();

    // Storage for the velocity and acceleration
    double veloc[n_position_type][n_dim];
    double accel[n_position_type][n_dim];

    // Find number of stored time values
    unsigned n_tstorage = ntstorage();

    // Loop over the variables
    for (unsigned i = 0; i < n_dim; i++)
    {
      for (unsigned k = 0; k < n_position_type; k++)
      {
        veloc[k][i] = accel[k][i] = 0.0;
        // Loop over all the history data
        for (unsigned t = 0; t < n_tstorage; t++)
        {
          veloc[k][i] += weight(1, t) * node_pt->x_gen(t, k, i);
          accel[k][i] += weight(2, t) * node_pt->x_gen(t, k, i);
        }
      }
    }

    // Loop over the position variables
    for (unsigned i = 0; i < n_dim; i++)
    {
      // If not a copy
      if (node_pt->position_is_a_copy(i) == false)
      {
        // Loop over position types
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Set previous values/veloc/accel to present values/veloc/accel,
          // if not a copy
          for (unsigned t = NSTEPS; t > 0; t--)
          {
            node_pt->x_gen(t, k, i) = node_pt->x_gen(t - 1, k, i);
          }

          node_pt->x_gen(NSTEPS + 1, k, i) = veloc[k][i];
          node_pt->x_gen(NSTEPS + 2, k, i) = accel[k][i];
        }
      }
    }
  }

  //=========================================================================
  /// Set weights
  //=========================================================================
  template<unsigned NSTEPS>
  void Newmark<NSTEPS>::set_weights()
  {
    double dt = Time_pt->dt(0);

    // Weights for second derivs
    Weight(2, 0) = 2.0 / (Beta2 * dt * dt);
    Weight(2, 1) = -2.0 / (Beta2 * dt * dt);
    for (unsigned t = 2; t <= NSTEPS; t++)
    {
      Weight(2, t) = 0.0;
    }
    Weight(2, NSTEPS + 1) = -2.0 / (dt * Beta2);
    Weight(2, NSTEPS + 2) = (Beta2 - 1.0) / Beta2;

    // Weights for first derivs.
    Weight(1, 0) = Beta1 * dt * Weight(2, 0);
    Weight(1, 1) = Beta1 * dt * Weight(2, 1);
    for (unsigned t = 2; t <= NSTEPS; t++)
    {
      Weight(1, t) = 0.0;
    }
    Weight(1, NSTEPS + 1) = 1.0 + Beta1 * dt * Weight(2, NSTEPS + 1);
    Weight(1, NSTEPS + 2) =
      dt * (1.0 - Beta1) + Beta1 * dt * Weight(2, NSTEPS + 2);
  }


  //===================================================================
  // Force build of templates: These are all the ones that might be
  // needed if Newmark is used together with BDF schemes (they require
  // up to 4 previous timesteps)
  //===================================================================
  template class Newmark<1>;
  template class Newmark<2>;
  template class Newmark<3>;
  template class Newmark<4>;


  //=========================================================================
  /// This function updates the Data's time history so that
  /// we can advance to the next timestep.
  //=========================================================================
  template<unsigned NSTEPS>
  void NewmarkBDF<NSTEPS>::shift_time_values(Data* const& data_pt)
  {
    // Find number of values stored
    unsigned n_value = data_pt->nvalue();

    // Storage for the velocity and acceleration
    Vector<double> veloc(n_value, 0.0);
    Vector<double> accel(n_value, 0.0);

    // Find number of stored time values
    unsigned n_tstorage = this->ntstorage();

    // Loop over the variables
    for (unsigned i = 0; i < n_value; i++)
    {
      // Loop over all the history data
      for (unsigned t = 0; t < n_tstorage; t++)
      {
        veloc[i] += Newmark_veloc_weight[t] * data_pt->value(t, i);
        accel[i] += this->weight(2, t) * data_pt->value(t, i);
      }
    }


    // Loop over values
    for (unsigned j = 0; j < n_value; j++)
    {
      // Set previous values/veloc/accel to present values/veloc/accel,
      // if not a copy
      if (data_pt->is_a_copy(j) == false)
      {
        for (unsigned t = NSTEPS; t > 0; t--)
        {
          data_pt->set_value(t, j, data_pt->value(t - 1, j));
        }
        data_pt->set_value(NSTEPS + 1, j, veloc[j]);
        data_pt->set_value(NSTEPS + 2, j, accel[j]);
      }
    }
  }

  //=========================================================================
  /// This function updates a nodal time history so that
  /// we can advance to the next timestep.
  //=========================================================================
  template<unsigned NSTEPS>
  void NewmarkBDF<NSTEPS>::shift_time_positions(Node* const& node_pt)
  {
    // Find the number of coordinates
    unsigned n_dim = node_pt->ndim();
    // Find the number of position types
    unsigned n_position_type = node_pt->nposition_type();

    // Storage for the velocity and acceleration
    double veloc[n_position_type][n_dim];
    double accel[n_position_type][n_dim];

    // Find number of stored time values
    unsigned n_tstorage = this->ntstorage();

    // Loop over the variables
    for (unsigned i = 0; i < n_dim; i++)
    {
      for (unsigned k = 0; k < n_position_type; k++)
      {
        veloc[k][i] = accel[k][i] = 0.0;
        // Loop over all the history data
        for (unsigned t = 0; t < n_tstorage; t++)
        {
          veloc[k][i] += Newmark_veloc_weight[t] * node_pt->x_gen(t, k, i);
          accel[k][i] += this->weight(2, t) * node_pt->x_gen(t, k, i);
        }
      }
    }

    // Loop over the position variables
    for (unsigned i = 0; i < n_dim; i++)
    {
      // If not a copy
      if (node_pt->position_is_a_copy(i) == false)
      {
        // Loop over position types
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Set previous values/veloc/accel to present values/veloc/accel,
          // if not a copy
          for (unsigned t = NSTEPS; t > 0; t--)
          {
            node_pt->x_gen(t, k, i) = node_pt->x_gen(t - 1, k, i);
          }

          node_pt->x_gen(NSTEPS + 1, k, i) = veloc[k][i];
          node_pt->x_gen(NSTEPS + 2, k, i) = accel[k][i];
        }
      }
    }
  }


  //=========================================================================
  /// Set weights
  //=========================================================================
  template<>
  void NewmarkBDF<1>::set_weights()
  {
    // Set Newmark weights for second derivatives
    double dt = Time_pt->dt(0);
    Weight(2, 0) = 2.0 / (Beta2 * dt * dt);
    Weight(2, 1) = -2.0 / (Beta2 * dt * dt);
    Weight(2, 2) = -2.0 / (dt * Beta2);
    Weight(2, 3) = (Beta2 - 1.0) / Beta2;

    // Set BDF weights for first derivatives
    Weight(1, 0) = 1.0 / dt;
    Weight(1, 1) = -1.0 / dt;
    Weight(1, 2) = 0.0;
    Weight(1, 3) = 0.0;

    // Orig Newmark weights for first derivs.
    set_newmark_veloc_weights(dt);
  }

  //=========================================================================
  /// Set weights
  //=========================================================================
  template<>
  void NewmarkBDF<2>::set_weights()
  {
    // Set Newmark weights for second derivatives
    double dt = Time_pt->dt(0);
    Weight(2, 0) = 2.0 / (Beta2 * dt * dt);
    Weight(2, 1) = -2.0 / (Beta2 * dt * dt);
    Weight(2, 2) = 0.0;
    Weight(2, 3) = -2.0 / (dt * Beta2);
    Weight(2, 4) = (Beta2 - 1.0) / Beta2;

    // Set BDF weights for first derivatives
    if (Degrade_to_bdf1_for_first_derivs)
    {
      this->Weight(1, 0) = 1.0 / dt;
      this->Weight(1, 1) = -1.0 / dt;
      unsigned nweights = this->Weight.ncol();
      for (unsigned i = 2; i < nweights; i++)
      {
        this->Weight(1, i) = 0.0;
      }
    }
    else
    {
      double dtprev = Time_pt->dt(1);
      Weight(1, 0) = 1.0 / dt + 1.0 / (dt + dtprev);
      Weight(1, 1) = -(dt + dtprev) / (dt * dtprev);
      Weight(1, 2) = dt / ((dt + dtprev) * dtprev);
      Weight(1, 3) = 0.0;
      Weight(1, 4) = 0.0;
    }

    // Orig Newmark weights for first derivs.
    set_newmark_veloc_weights(dt);
  }

  //=========================================================================
  /// Set weights
  //=========================================================================
  template<>
  void NewmarkBDF<4>::set_weights()
  {
    // Set Newmark weights for second derivatives
    double dt = Time_pt->dt(0);
    Weight(2, 0) = 2.0 / (Beta2 * dt * dt);
    Weight(2, 1) = -2.0 / (Beta2 * dt * dt);
    Weight(2, 2) = 0.0;
    Weight(2, 3) = 0.0;
    Weight(2, 4) = 0.0;
    Weight(2, 5) = -2.0 / (dt * Beta2);
    Weight(2, 6) = (Beta2 - 1.0) / Beta2;

    // Set BDF weights for first derivatives
#ifdef PARANOID
    for (unsigned i = 0; i < Time_pt->ndt(); i++)
    {
      if (dt != Time_pt->dt(i))
      {
        throw OomphLibError("BDF4 currently only works for fixed timesteps \n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Set BDF weights for first derivatives
    if (Degrade_to_bdf1_for_first_derivs)
    {
      this->Weight(1, 0) = 1.0 / dt;
      this->Weight(1, 1) = -1.0 / dt;
      unsigned nweights = this->Weight.ncol();
      for (unsigned i = 2; i < nweights; i++)
      {
        this->Weight(1, i) = 0.0;
      }
    }
    else
    {
      Weight(1, 0) = 25.0 / 12.0 / dt;
      Weight(1, 1) = -48.0 / 12.0 / dt;
      Weight(1, 2) = 36.0 / 12.0 / dt;
      Weight(1, 3) = -16.0 / 12.0 / dt;
      Weight(1, 4) = 3.0 / 12.0 / dt;
      Weight(1, 5) = 0.0;
      Weight(1, 6) = 0.0;
    }

    // Orig Newmark weights for first derivs.
    set_newmark_veloc_weights(dt);
  }

  //===================================================================
  // Force build of templates.
  //===================================================================
  template class NewmarkBDF<1>;
  template class NewmarkBDF<2>;
  template class NewmarkBDF<4>;


} // namespace oomph
