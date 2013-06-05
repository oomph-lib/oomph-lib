#ifndef OOMPH_MIDPOINT_METHOD_H
#define OOMPH_MIDPOINT_METHOD_H

/*
  description of file goes here
*/

// //oomph-lib headers
// #include "Vector.h"
// #include "nodes.h"
// #include "matrices.h"
// #include "timesteppers.h"

#include "generic.h"

#include "./poly_interp.h"

using namespace oomph;

namespace oomph
{

  // =======================================================================
  /// Implicit midpoint method, implemented via a BDF1 step followed by an
  /// update.

  /// Common mistakes when using midpoint:

  /// * Didn't include the d_u_evaltime_by_du_np1 term in the Jacobian.

  /// * Didn't properly interpolate time/values/x/derivatives in
  /// jacobian/residual (all of these MUST be evaluated at the midpoint).


  //========================================================================
  class MidpointMethod : public TimeStepper
  {
  public:

    /// Constructor with initialisation
    MidpointMethod(bool adaptive=false, unsigned n_interpolation_points=2,
                   double fudge_factor=1.0) :
      TimeStepper(2,1), // initialise weights later
      Fudge_factor(fudge_factor),
      N_interp(n_interpolation_points),
      Predictor_storage_index(nprev_values()+1),
      Dy_tnph_storage_index(nprev_values()+2)
    {
      Adaptive_Flag = adaptive;
      Is_steady = false;
      Type = "Midpoint method";

      // Storage for weights needs to be 2x(n_interp +1) (or at least
      // 2x2). This means we provide ways to calculate the zeroth and first
      // derivatives and in calculations we use n_interp+1 previous time
      // values.

      // But actually (for some reason...) the amount of data storage
      // allocated for the timestepper in each node is determined by the
      // size of weight. So we need to store two extra dummy entries in
      // order to have storage space for the auxillary values used in
      // adaptive timestepping.
      Weight.resize(2, nprev_values() + 1 + 2, 0.0);

      // Set the weights for zero-th derivative (i.e. the value to use in
      // newton solver calculations, implicit midpoint method uses the
      // average of previous and current values).
      Weight(0, 0) = 0.5;
      Weight(0, 1) = 0.5;
    }

    /// Destructor
    virtual ~MidpointMethod() {}

    /// Setup weights for time derivative calculations.
    void set_weights()
    {
      double dt=Time_pt->dt(0);
      Weight(1,0) = 1.0/dt;
      Weight(1,1) = -1.0/dt;
    }

    /// Actual order (accuracy) of the scheme
    unsigned order() {return 2;}

    /// Number of timestep increments that are required by the scheme
    unsigned ndt() {return nprev_values();}

    /// Number of previous values that actually hold previous values
    unsigned nprev_values() {return 1 + N_interp;}
    // unsigned nprev_values() {return 2;}

    /// \short ??ds
    unsigned nprev_values_for_value_at_evaluation_time() {return 2;}

    /// \short This function advances the Data's time history so that
    /// we can move on to the next timestep
    void shift_time_values(Data* const &data_pt);

    /// \short This function advances the time history of the positions
    /// at a node.
    void shift_time_positions(Node* const &node_pt);

    /// \short Set the weights for the error computation. This is not used
    /// by the midpoint method since interpolation is needed.
    void set_error_weights() {}

    /// \short Set the weights for the predictor previous timestep. This is not
    /// used by the midpoint method since interpolation is needed.
    void set_predictor_weights() {}

    /// not implemented (TODO)
    void assign_initial_values_impulsive(Data* const &data_pt) {abort();}
    void assign_initial_positions_impulsive(Node* const &node_pt) {abort();}
    void calculate_predicted_positions(Node* const &node_pt) {}
    double temporal_error_in_position(Node* const &node_pt, const unsigned &i)
    { abort(); return 0.0;}
    void undo_make_steady(){abort();}

    // Adaptivity
    void calculate_predicted_values(Data* const &data_pt);
    double temporal_error_in_value(Data* const &data_pt, const unsigned &i);

    double Fudge_factor;

  private:

    /// Number of interpolation points to use.
    unsigned N_interp;

    /// Dummy time index of predictor values in data.
    unsigned Predictor_storage_index;

    /// Dummy time index of y'(t_nph) in data.
    unsigned Dy_tnph_storage_index;

    /// Inaccessible copy constructor.
    MidpointMethod(const MidpointMethod &dummy) {}

    /// Inaccessible assignment operator.
    void operator=(const MidpointMethod &dummy) {}
  };


  /// \short This function advances the Data's time history so that
  /// we can move on to the next timestep
  void MidpointMethod::shift_time_values(Data* const &data_pt)
  {
    //Loop over the values, set previous values to the previous value, if
    //not a copy.
    for(unsigned j=0, nj=data_pt->nvalue(); j<nj; j++)
      {
        if(!data_pt->is_a_copy(j))
          {
            for(unsigned t=ndt(); t>0; t--)
              {
                data_pt->set_value(t,j,data_pt->value(t-1,j));
              }
          }
      }
  }

  ///\short This function advances the time history of the positions
  ///at a node. ??ds I have no idea what I'm doing here!
  void MidpointMethod::shift_time_positions(Node* const &node_pt)
  {
   //Find the number of coordinates
   unsigned n_dim = node_pt->ndim();
   //Find the number of position types
   unsigned n_position_type = node_pt->nposition_type();

   //Find number of stored timesteps
   unsigned n_tstorage = ntstorage();

   //Storage for the velocity
   double velocity[n_position_type][n_dim];

   //If adaptive, find the velocities
   if(adaptive_flag())
    {
     //Loop over the variables
     for(unsigned i=0;i<n_dim;i++)
      {
       for(unsigned k=0;k<n_position_type;k++)
        {
         //Initialise velocity to zero
         velocity[k][i] =0.0;
         //Loop over all history data
         for(unsigned t=0;t<n_tstorage;t++)
          {
           velocity[k][i] += Weight(1,t)*node_pt->x_gen(t,k,i);
          }
        }
      }
    }

   //Loop over the positions
   for(unsigned i=0;i<n_dim;i++)
    {
     //If the position is not a copy
     if(node_pt->position_is_a_copy(i) == false)
      {
       //Loop over the position types
       for(unsigned k=0;k<n_position_type;k++)
        {
         //Loop over stored times, and set values to previous values
         for(unsigned t=ndt();t>0;t--)
          {
           node_pt->x_gen(t,k,i) = node_pt->x_gen(t-1,k,i);
          }

         //If we are using the adaptive scheme, set the velocity
         if(adaptive_flag())
          {
           node_pt->x_gen(ndt()+1,k,i) = velocity[k][i];
          }
        }
      }
    }
  }

  ///
  void MidpointMethod::calculate_predicted_values(Data* const &data_pt)
  {
    if(adaptive_flag())
      {
        // Cache some values
        double tnph = (time_pt()->time(0) + time_pt()->time(1))/2;
        double dt_n = time_pt()->dt(); //??ds is this right?

        // Collect the data needed for interpolation
        Vector<double> known_times(N_interp);
        Vector<Vector<double> > known_values(N_interp);
        unsigned nvalue = data_pt->nvalue();
        for(unsigned i_t=0; i_t<N_interp; i_t++)
          {
            known_times[i_t] = time_pt()->time(i_t);

            for(unsigned v=0; v<nvalue; v++)
              {
                known_values[i_t].push_back(data_pt->value(i_t, v));
              }
          }

        // Interpolate y(t_n+1/2), y'_n and y'(t_n+1/2)
        BarycentricLagrangeInterpolator interpolator(known_times, known_values);

        Vector<double> y_tnph, dy_tnph, dy_tn;
        interpolator.eval(tnph, y_tnph);
        interpolator.eval_derivative(tnph, 1, dy_tnph);
        interpolator.eval_derivative(time_pt()->time(1), 1, dy_tn);

        // Use the interpolated values to calculate an estimate to
        // y_np1. Basically just using AB2.
        for(unsigned v=0, nv=data_pt->nvalue(); v<nv; v++)
          {
            data_pt->set_value(Predictor_storage_index, v,
                               y_tnph[v] + (dt_n/4) * ( 3*dy_tnph[v] - dy_tn[v]));

            // Also store the interpolated values of y'(t_nph) for later
            data_pt->set_value(Dy_tnph_storage_index, v, dy_tnph[v]);
          }
      }
  }

  /// \short
  double MidpointMethod::temporal_error_in_value(Data* const &data_pt,
                                                 const unsigned &i)
  {
    if(adaptive_flag())
      {
        // Calculate error. Hopefully see my paper ??ds
        double y_np1_MP = data_pt->value(0, i);
        double y_n = data_pt->value(1, i);
        double y_np1_pred = data_pt->value(Predictor_storage_index, i);
        double dy_tnph = data_pt->value(Dy_tnph_storage_index, i);;
        double dt_n = time_pt()->dt(); //??ds is this right?
        double a_n = dt_n*dy_tnph + y_n - y_np1_MP;

        return (4*(y_np1_MP - y_np1_pred) + 5*a_n) * Fudge_factor;
      }
    else
      {
        std::string err("Tried to get the temporal error from a non-adaptive");
        err += " time stepper.";
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
  }

} // End of oomph namespace

#endif
