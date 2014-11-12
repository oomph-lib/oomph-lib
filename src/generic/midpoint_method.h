#ifndef OOMPH_MIDPOINT_METHOD_H
#define OOMPH_MIDPOINT_METHOD_H

//oomph-lib headers
#include "Vector.h"
#include "nodes.h"
#include "matrices.h"
#include "timesteppers.h"
#include "explicit_timesteppers.h"

namespace oomph
{

 // Forward decl. so that we can have function of Problem*
 class Problem;
 

  // =======================================================================
  /// Implicit midpoint rule base class for the two implementations
  //========================================================================
  class MidpointMethodBase : public TimeStepper
  {
  public:

   /// Constructor with initialisation
   MidpointMethodBase(bool adaptive=false) :
    TimeStepper(2,1) // initialise weights later
   {
    Adaptive_Flag = adaptive;
    Is_steady = false;
    Type = "Midpoint method";
    Predict_by_explicit_step=true;
    Predictor_pt = 0;

    // If adaptive then we are storing predicted values in slot
    // 4. Otherwise we aren't storing them so leave it as -1.
    if(adaptive) 
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
    virtual ~MidpointMethodBase() 
    {
     delete Predictor_pt; Predictor_pt = 0;
    }

    /// Setup weights for time derivative calculations.
    virtual void set_weights()=0;

   /// \short ??ds
   virtual unsigned nprev_values_for_value_at_evaluation_time() const=0;

    /// Actual order (accuracy) of the scheme
    unsigned order() {return 2;}

    /// Number of timestep increments that are required by the scheme
    unsigned ndt() {return nprev_values();}

    /// ??ds
    unsigned nprev_values() {return 4;}

    /// Get the pointer to the timestepper to use as a predictor in
    /// adaptivity.
    ExplicitTimeStepper* predictor_pt() {return Predictor_pt;}

    /// Check that the predicted values are the ones we want.
    void check_predicted_values_up_to_date() const
     {
#ifdef PARANOID
      if(std::abs(time_pt()->time() - Predicted_time) > 1e-15)
       {
        throw OomphLibError("Predicted values are not from the correct time step",
                            OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
       }
#endif
     }

    /// \short This function advances the Data's time history so that
    /// we can move on to the next timestep
    void shift_time_values(Data* const &data_pt);

    /// \short This function advances the time history of the positions
    /// at a node.
    void shift_time_positions(Node* const &node_pt);

    /// \short Set the weights for the error computation. This is not used
    /// by midpoint rule.
    void set_error_weights() {}

    /// \short Set the weights for the predictor previous timestep. This is not
    /// used by midpint rule.
    void set_predictor_weights() {}

    /// not implemented (??ds TODO)
    void assign_initial_values_impulsive(Data* const &data_pt) {abort();}
    void assign_initial_positions_impulsive(Node* const &node_pt) {abort();}
    void calculate_predicted_positions(Node* const &node_pt) {}
    double temporal_error_in_position(Node* const &node_pt, const unsigned &i)
    { abort(); return 0.0;}

    // Adaptivity
    void calculate_predicted_values(Data* const &data_pt);
    double temporal_error_in_value(Data* const &data_pt, const unsigned &i);


   /// Store the time that the predicted values currently stored are at,
   /// to compare for paranoid checks.
   double Predicted_time;

   void set_predictor_pt(ExplicitTimeStepper* _pred_pt)
    {
     Predictor_pt = _pred_pt;
    }

  private:

    /// Time stepper to use to calculate predictor value
    ExplicitTimeStepper* Predictor_pt;
  };

 /// Implicit midpoint rule implemented by calculation of residuals etc. at
 /// half step.
class MidpointMethod : public MidpointMethodBase
{
public:
 /// Common mistakes when using this implementation of midpoint:
 /// * Didn't include the d_u_evaltime_by_du_np1 term in the Jacobian.
 /// * Didn't properly interpolate time/values/x/derivatives in
 ///   jacobian/residual (all of these MUST be evaluated at the midpoint).


 /// Constructor with initialisation
 MidpointMethod(bool adaptive=false) : MidpointMethodBase(adaptive) {}

 /// Destructor, predictor_pt handled by base
 virtual ~MidpointMethod() {}

 /// Setup weights for time derivative calculations.
 void set_weights()
 {
  // Set the weights for zero-th derivative (i.e. the value to use in
  // newton solver calculations, implicit midpoint method uses the
  // average of previous and current values).
  Weight(0, 0) = 0.5;
  Weight(0, 1) = 0.5;

  // Set the weights for 1st time derivative
  double dt=Time_pt->dt(0);
  Weight(1,0) = 1.0/dt;
  Weight(1,1) = -1.0/dt;
 }

 /// \short ??ds
 unsigned nprev_values_for_value_at_evaluation_time() const {return 2;}

private:

 /// Inaccessible copy constructor.
 MidpointMethod(const MidpointMethod &dummy) {}

 /// Inaccessible assignment operator.
 void operator=(const MidpointMethod &dummy) {}
};


/// Implementation of implicit midpoint rule by taking half a step of bdf1
/// then applying an update to all dofs.
class MidpointMethodByBDF : public MidpointMethodBase
{
public:
 /// Constructor with initialisation
 MidpointMethodByBDF(bool adaptive=false) : MidpointMethodBase(adaptive) {}

 /// Destructor
 virtual ~MidpointMethodByBDF() {}

 /// Setup weights for time derivative calculations.
 void set_weights()
 {
  // Use weights from bdf1
  double dt=Time_pt->dt(0);
  Weight(1,0) = 1.0/dt;
  Weight(1,1) = -1.0/dt;
 }

 /// \short Evaluation time is the end of the bdf1 "half-step", so only need one
 /// value.
 unsigned nprev_values_for_value_at_evaluation_time() const {return 1;}

 /// Half the timestep before starting solve
 void actions_before_timestep(Problem* problem_pt);

 /// Take problem from t={n+1/2} to t=n+1 by algebraic update and restore
 /// time step.
 void actions_after_timestep(Problem* problem_pt);

private:

 /// Time stepper to use to calculate predictor value
 ExplicitTimeStepper* Predictor_pt;

 /// Inaccessible copy constructor.
 MidpointMethodByBDF(const MidpointMethodByBDF &dummy) {}

 /// Inaccessible assignment operator.
 void operator=(const MidpointMethodByBDF &dummy) {}

};

} // End of oomph namespace

#endif
