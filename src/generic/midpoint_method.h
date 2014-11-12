#ifndef OOMPH_MIDPOINT_METHOD_H
#define OOMPH_MIDPOINT_METHOD_H

//oomph-lib headers
#include "Vector.h"
#include "nodes.h"
#include "matrices.h"
#include "timesteppers.h"

#include "./poly_interp.h"

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
    unsigned nprev_values_for_value_at_evaluation_time() const {return 2;}

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

    /// not implemented (??ds TODO)
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



} // End of oomph namespace

#endif
