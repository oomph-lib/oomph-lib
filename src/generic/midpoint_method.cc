
#include "midpoint_method.h"
#include "problem.h"

namespace oomph
{

 /// \short This function advances the Data's time history so that
 /// we can move on to the next timestep
 void MidpointMethodBase::shift_time_values(Data* const &data_pt)
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
 void MidpointMethodBase::shift_time_positions(Node* const &node_pt)
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
       }
     }
   }
 }


 /// Dummy - just check that the values that
 /// problem::calculate_predicted_values() has been called right.
 void MidpointMethodBase::calculate_predicted_values(Data* const &data_pt)
 {
  if(adaptive_flag())
   {
    // Can't do it here, but we can check that the predicted values have
    // been updated.
    check_predicted_values_up_to_date();
   }
 }


 double MidpointMethodBase::temporal_error_in_value(Data* const &data_pt,
                                                    const unsigned &i)
 {
  if(adaptive_flag())
   {
    // predicted value is more accurate so just compare with that
    //??ds is sign important? probably not...
    return data_pt->value(i) - data_pt->value(Predictor_storage_index, i);
   }
  else
   {
    std::string err("Tried to get the temporal error from a non-adaptive");
    err += " time stepper.";
    throw OomphLibError(err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
   }
 }

 /// Half the timestep before starting solve
 void MidpointMethodByBDF::actions_before_timestep(Problem* problem_pt)
  {

   // Check that this is the only time stepper
#ifdef PARANOID
   if(problem_pt->ntime_stepper() != 1)
    {
     std::string err = "This implementation of midpoint can only work with a ";
     err += "single time stepper, try using MidpointMethod instead (but check ";
     err += "your Jacobian and residual calculations very carefully for compatability).";
     throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                         OOMPH_CURRENT_FUNCTION);
    }
#endif

   time_pt()->dt() /= 2;
   time_pt()->time() -= time_pt()->dt();

   // Set the weights again because we just changed the step size
   set_weights();
  }

 /// Take problem from t={n+1/2} to t=n+1 by algebraic update and restore
 /// time step.
 void MidpointMethodByBDF::actions_after_timestep(Problem* problem_pt)
  {
   const unsigned ndof = problem_pt->ndof();

   // Get dofs at previous time step
   DoubleVector dof_n;
   problem_pt->get_dofs(1, dof_n);

   // Update dofs at current step
   for(unsigned i=0; i<ndof; i++)
    {
     problem_pt->dof(i) += problem_pt->dof(i) - dof_n[i];
    }

   // Update time
   problem_pt->time_pt()->time() += problem_pt->time_pt()->dt();

   // Return step size to its full value
   problem_pt->time_pt()->dt() *= 2;

  }

} // End of oomph namespace
