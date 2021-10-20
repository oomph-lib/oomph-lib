//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================

// Include oomph-lib headers
#include "generic.h"

namespace oomph
{
 
 //====================================================================
 /// Self-starting BDF2 timestepper class
 //====================================================================
 class SelfStartingBDF2 : public TimeStepper
 {
   private:
  
  /// Private data for the predictor weights
  Vector<double> Predictor_weight;
  
  /// Private data for the error weight 
  double Error_weight;
  
  ///  Bool to indicate if the timestepper is in `BDF1' mode
  /// i.e. it behaves just as if it were a BDF1 timestepper.
  /// This status may be achieved temporarily by calling
  /// turn_on_bdf1_mode(). It's reset to the default BDF2 behaviour
  /// by calling turn_off_bdf1_mode().
  bool BDF1_mode;
  
   public:
  
  /// Constructor for the case when we allow adaptive timestepping.
  /// Set BDF1_mode flag to false (default is BDF2 mode)
  SelfStartingBDF2(const bool& adaptive=false) : TimeStepper(3,1), BDF1_mode(false)
   {
    Type="BDF";
    
    // If it's adaptive, we need to allocate additional space to 
    // carry along a prediction and an acceleration
    if(adaptive)
     {
      // Set the adaptive flag to be true
      Adaptive_Flag = true;
      
      // Set the size of the Predictor_Weights Vector
      // N.B. The size is correct for BDF<2>, but may be wrong for others
      Predictor_weight.resize(4);
      
      // The order of the scheme is 1, i.e. Weights has two entries
      
      //Resize the weights to the appropriate size
      Weight.resize(2,5);
      
      // Initialise
      for(unsigned i=0;i<2;i++)
       {
        // Initialise each weight to zero
        for(unsigned j=0;j<5;j++)
         {
          Weight(i,j) = 0.0;
         }
       }
      // Set the weight for the zero-th derivative
      Weight(0,0) = 1.0;
     }
   }
   
   /// Broken copy constructor
   SelfStartingBDF2(const SelfStartingBDF2&) 
    { 
     BrokenCopy::broken_copy("SelfStartingBDF2");
    } 
   
   /// Broken assignment operator
   void operator=(const SelfStartingBDF2&) 
    {
     BrokenCopy::broken_assign("SelfStartingBDF2");
    }
   
   /// Return the actual order of the scheme
   unsigned order() const { return 2; }
   
   ///  Function to make the time stepper temporarily work
   /// as if it is BDF1. This is trivially achieved by resetting the
   /// weights to those of BDF1, before setting the weight corresponding
   /// to history value 2 to zero.
   void turn_on_bdf1_mode()
   {
    // Update flag
    BDF1_mode = true;
    
    // Reset weights
    set_weights_bdf1();
   }
   
   /// Flag to indicate if the timestepper is working in BDF1 mode
   bool bdf1_mode() { return BDF1_mode; }
 
   ///  Reset the timestepper to BDF2: Set the "BDF1_mode"
   /// flag to false and re-compute the weights
   void turn_off_bdf1_mode() 
   {
    // Update flag
    BDF1_mode = false;
    
    // Reset_weights
    set_weights_bdf2();
   }

   ///  Reset the is_steady status to its default (false) and
   /// re-compute the weights
   void undo_make_steady() 
   {
    // Update flag
    Is_steady = false;
    
    // Reset weights
    if(BDF1_mode)
     {
      set_weights_bdf1();
     }
    else
     {
      set_weights_bdf2();
     }
   }

   ///  Initialise the time-history for the Data values,
   /// corresponding to an impulsive start.
   void assign_initial_values_impulsive(Data* const &data_pt)
   {
    // Find number of values stored
    const unsigned n_value = data_pt->nvalue();
    
    // Loop over values
    for(unsigned j=0;j<n_value;j++)
     {
      // Set previous values to the initial value, if not a copy
      if(data_pt->is_a_copy(j) == false)
       {
        for(unsigned t=1;t<=2;t++)
         {
          data_pt->set_value(t,j,data_pt->value(j));
         }
        
        // If it's adaptive 
        if(adaptive_flag())
         {
          // Initial velocity is zero
          data_pt->set_value(3,j,0.0);
          // Initial prediction is the value
          data_pt->set_value(4,j,data_pt->value(j));
         }
       }
     }
   }
   
   
   ///  Initialise the time-history for the nodal positions
   /// corresponding to an impulsive start.
   void assign_initial_positions_impulsive(Node* const &node_pt)
   {
    // Find the dimension of the node
    const unsigned n_dim = node_pt->ndim();
    // Find the number of position types at the node
    const unsigned n_position_type = node_pt->nposition_type();
    
    // Loop over the position variables
    for(unsigned i=0;i<n_dim;i++)
     {
      // If the position is not copied 
      // We copy entire coordinates at once
      if(node_pt->position_is_a_copy(i) == false)
       {
        // Loop over the position types
        for(unsigned k=0;k<n_position_type;k++)
         {
          // Set previous values to the initial value, if not a copy
          for(unsigned t=1;t<=2;t++)
           {
            node_pt->x_gen(t,k,i) = node_pt->x_gen(k,i);
           }
          
          // If it's adaptive 
          if(adaptive_flag())
           {
            // Initial mesh velocity is zero
            node_pt->x_gen(3,k,i) = 0.0;
            // Initial prediction is the value
            node_pt->x_gen(4,k,i) = node_pt->x_gen(k,i);
           }
         }
       }
     }
   }
   
   
   ///  Typedef for function that returns the (scalar) initial
   /// value at a given value of the continuous time t.
   typedef double (*InitialConditionFctPt)(const double& t);
   
   ///  Initialise the time-history for the Data values,
   /// corresponding to given time history, specified by
   /// Vector of function pointers.
   void assign_initial_data_values(Data* const &data_pt, 
                                   Vector<InitialConditionFctPt> 
                                   initial_value_fct)
   {
    // The time history stores the previous function values
    const unsigned n_time_value = ntstorage();
    
    // Find number of values stored
    const unsigned n_value = data_pt->nvalue();
    
    // Loop over current and stored timesteps
    for(unsigned t=0;t<n_time_value;t++)
     {
      // Get corresponding continous time
      double time = Time_pt->time(t);
      
      // Loop over values
      for(unsigned j=0;j<n_value;j++)
       {
        data_pt->set_value(t,j,initial_value_fct[j](time));
       }
     }
   }
   
   ///  This function updates the Data's time history so that
   /// we can advance to the next timestep. For BDF schemes,
   /// we simply push the values backwards...
   void shift_time_values(Data* const &data_pt)
   {
    // Find number of values stored
    const unsigned n_value = data_pt->nvalue();
    // Storage for velocity need to be here to be in scope
    Vector<double> velocity(n_value);
    
    // If adaptive, find the velocities
    if(adaptive_flag()) { time_derivative(1,data_pt,velocity); }
    
    // Loop over the values
    for(unsigned j=0;j<n_value;j++)
     {
      // Set previous values to the previous value, if not a copy
      if(data_pt->is_a_copy(j) == false)
       {
        // Loop over times, in reverse order
        for(unsigned t=2;t>0;t--)
         {
          data_pt->set_value(t,j,data_pt->value(t-1,j));
         }
        
        // If we are using the adaptive scheme
        if(adaptive_flag())
         {
          // Set the velocity
          data_pt->set_value(3,j,velocity[j]);
         }
       } 
     }
   }
   
   ///  This function advances the time history of the positions
   /// at a node. 
   void shift_time_positions(Node* const &node_pt)
   {
    // Find the number of coordinates
    const unsigned n_dim = node_pt->ndim();
    // Find the number of position types
    const unsigned n_position_type = node_pt->nposition_type();
    
    // Find number of stored timesteps
    const unsigned n_tstorage = ntstorage();
    
    // Storage for the velocity
    double velocity[n_position_type][n_dim];
    
    // If adaptive, find the velocities
    if(adaptive_flag())
     {
      // Loop over the variables
      for(unsigned i=0;i<n_dim;i++)
       {
        for(unsigned k=0;k<n_position_type;k++)
         {
          // Initialise velocity to zero
          velocity[k][i] =0.0;
          // Loop over all history data
          for(unsigned t=0;t<n_tstorage;t++)
           {
            velocity[k][i] += Weight(1,t)*node_pt->x_gen(t,k,i);
           }
         }
       }
     }
    
    // Loop over the positions
    for(unsigned i=0;i<n_dim;i++)
     {
      // If the position is not a copy
      if(node_pt->position_is_a_copy(i) == false)
       {
        // Loop over the position types
        for(unsigned k=0;k<n_position_type;k++)
         {
          // Loop over stored times, and set values to previous values
          for(unsigned t=2;t>0;t--)
           {
            node_pt->x_gen(t,k,i) = node_pt->x_gen(t-1,k,i);
           }
          
          // If we are using the adaptive scheme, set the velocity
          if(adaptive_flag())
           {
            node_pt->x_gen(3,k,i) = velocity[k][i];
           }
         }
       }
     }
   }

   ///  Implementation of pure virtual function in base class.
   /// Set weights corresponding to either a BDF1 scheme or a BDF2
   /// scheme, depending on the flag.
   void set_weights()
   {
    if(BDF1_mode) { set_weights_bdf1(); }
    else { set_weights_bdf2(); }
   }

   /// Set the weights to those corresponding to BDF1
   void set_weights_bdf1();
  
   /// Set the weights to those corresponding to BDF2
   void set_weights_bdf2();
   
   /// Number of previous values available.
   unsigned nprev_values() const { return 2; }
   
   /// Number of timestep increments that need to be stored by the scheme
   unsigned ndt() const { return 2; }
   
   /// Function to set the predictor weights corresponding to BDF1
   void set_predictor_weights_bdf1();
   
   /// Function to set the predictor weights corresponding to BDF2
   void set_predictor_weights_bdf2();

   /// Function to calculate predicted positions at a node (BDF1)
   void calculate_predicted_positions_bdf1(Node* const &node_pt);
   
   /// Function to calculate predicted positions at a node (BDF2)
   void calculate_predicted_positions_bdf2(Node* const &node_pt);

   /// Function to calculate predicted data values in a Data object (BDF1)
   void calculate_predicted_values_bdf1(Data* const &data_pt);
   
   /// Function to calculate predicted data values in a Data object (BDF2)
   void calculate_predicted_values_bdf2(Data* const &data_pt);

   /// Function to set the error weights corresponding to BDF1
   void set_error_weights_bdf1();

   /// Function to set the error weights corresponding to BDF2
   void set_error_weights_bdf2();

   /// Compute the error in the position i at a node (BDF1)
   double temporal_error_in_position_bdf1(Node* const &node_pt, 
                                          const unsigned &i);

   /// Compute the error in the position i at a node (BDF2)
   double temporal_error_in_position_bdf2(Node* const &node_pt, 
                                          const unsigned &i);
   
   /// Compute the error in the value i in a Data structure (BDF1)
   double temporal_error_in_value_bdf1(Data* const &data_pt, 
                                       const unsigned &i);

   /// Compute the error in the value i in a Data structure (BDF2)
   double temporal_error_in_value_bdf2(Data* const &data_pt, 
                                       const unsigned &i);
   
 };



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

 
 // BDF1 functions
 // --------------

 
 //=======================================================================
 /// Assign the values of the weights (BDF1)
 //=======================================================================
 void  SelfStartingBDF2::set_weights_bdf1()
 {
  // Set BDF1 weights as normal
  double dt=Time_pt->dt(0);
  Weight(1,0) = 1.0/dt;
  Weight(1,1) = -1.0/dt;
  
  // Set weight associated with history value 2 to zero so that it has
  // no effect on the dudt calculation
  Weight(1,2) = 0.0;
 }
 
 
 //======================================================================
 /// Calculate the predictor weights (BDF1)
 //======================================================================
 void SelfStartingBDF2::set_predictor_weights_bdf1()
 {
  //Read the value of the previous timesteps
  //double dt=Time_pt->dt(0);
  //double dtprev=Time_pt->dt(1);
  
  throw OomphLibError("Not implemented yet",
                      "SelfStartingBDF2::set_predictor_weights_bdf1()",
                      OOMPH_EXCEPTION_LOCATION);
 }


 //=======================================================================
 /// Calculate the predicted values and store them at the appropriate
 /// location in the data structure (BDF1)
 /// This function must be called after the time-values have been shifted!
 //=======================================================================
 void SelfStartingBDF2::calculate_predicted_values_bdf1(Data* const &data_pt)
 {
  throw OomphLibError("Not implemented yet",
                      "SelfStartingBDF2::calculate_predicted_weights_bdf1()",
                      OOMPH_EXCEPTION_LOCATION);
 }
 
 
 //=======================================================================
 /// Calculate predictions for the positions (BDF1)
 //=======================================================================
 void SelfStartingBDF2::calculate_predicted_positions_bdf1(
  Node* const &node_pt)
 {
  throw OomphLibError("Not implemented yet",
                      "SelfStartingBDF2::calculate_predicted_positions_bdf1()",
                      OOMPH_EXCEPTION_LOCATION);
 }
 

 //=======================================================================
 /// Function that sets the error weights (BDF1)
 //=======================================================================
 void SelfStartingBDF2::set_error_weights_bdf1()
 {
  throw OomphLibError("Not implemented yet",
                      "SelfStartingBDF2::set_error_weights_bdf1()",
                      OOMPH_EXCEPTION_LOCATION);
 }
 
 //===================================================================
 /// Function to compute the error in position i at node (BDF1)
 //===================================================================
 double SelfStartingBDF2::temporal_error_in_position_bdf1(
  Node* const &node_pt,const unsigned &i)
 {
  throw OomphLibError("Not implemented yet",
                      "SelfStartingBDF2::temporal_error_in_position_bdf1()",
                      OOMPH_EXCEPTION_LOCATION);
  return 0.0;
 }
 
 
 //=========================================================================
 /// Function to calculate the error in the data value i (BDF1)
 //=========================================================================
 double SelfStartingBDF2::temporal_error_in_value_bdf1(Data* const &data_pt,
                                                       const unsigned &i)
 {
  throw OomphLibError("Not implemented yet",
                      "SelfStartingBDF2::temporal_error_in_value_bdf1()",
                      OOMPH_EXCEPTION_LOCATION);
  return 0.0;
 }
 
 
 // BDF2 functions
 // --------------


 //========================================================================
 /// Assign the values of the weights; pass the value of the timestep (BDF2)
 //========================================================================
 void SelfStartingBDF2::set_weights_bdf2()
 {
  double dt=Time_pt->dt(0);
  double dtprev=Time_pt->dt(1);
  Weight(1,0) =  1.0/dt + 1.0/(dt + dtprev);
  Weight(1,1) = -(dt + dtprev)/(dt*dtprev);
  Weight(1,2) = dt/((dt+dtprev)*dtprev);
  
  if(adaptive_flag())
   {
    Weight(1,3) = 0.0;
    Weight(1,4) = 0.0;
   }
 }
 
 //======================================================================
 /// Calculate the predictor weights (BDF2)
 //======================================================================
 void SelfStartingBDF2::set_predictor_weights_bdf2()
 {
  //If it's adaptive set the predictor weights
  if(adaptive_flag())
   {
    //Read the value of the previous timesteps
    double dt=Time_pt->dt(0);
    double dtprev=Time_pt->dt(1);
    
    //Set the predictor weights
    Predictor_weight[0] = 0.0;
    Predictor_weight[1] = 1.0 - (dt*dt)/(dtprev*dtprev);
    Predictor_weight[2] = (dt*dt)/(dtprev*dtprev);
    //Acceleration term
    Predictor_weight[3] = (1.0 + dt/dtprev)*dt;
   }
 }
 
 //=======================================================================
 /// Calculate the predicted values and store them at the appropriate
 /// location in the data structure (BDF2)
 /// This function must be called after the time-values have been shifted!
 //=======================================================================
 void SelfStartingBDF2::calculate_predicted_values_bdf2(Data* const &data_pt)
 {
  //If it's adaptive calculate the values
  if(adaptive_flag())
   {
    //Find number of values
    unsigned n_value = data_pt->nvalue();
    //Loop over the values
    for(unsigned j=0;j<n_value;j++)
     {
      //If the value is not copied
      if(data_pt->is_a_copy(j) == false)
       {
        //Now Initialise the predictor to zero
        double predicted_value = 0.0;
        //Now loop over all the stored data and add appropriate values
        //to the predictor
        for(unsigned i=1;i<4;i++)
         {
          predicted_value += data_pt->value(i,j)*Predictor_weight[i];
         }
        //Store the predicted value
        //Note that predictor is stored as the FIFTH entry in this scheme
        data_pt->set_value(4,j,predicted_value);
       }
     }
   }
 }
 
 //=======================================================================
 /// Calculate predictions for the positions (BDF2)
 //=======================================================================
 void SelfStartingBDF2::calculate_predicted_positions_bdf2(
  Node* const &node_pt)
 {
  //Only do this if adaptive
  if(adaptive_flag())
   {
    //Find number of dimensions of the problem
    unsigned n_dim = node_pt->ndim();
    //Loop over the dimensions
    for(unsigned j=0;j<n_dim;j++)
     {
      //If the node is not copied
      if(node_pt->position_is_a_copy(j) == false)
       {
        //Initialise the predictor to zero
        double predicted_value = 0.0;
        //Now loop over all the stored data and add appropriate values
        //to the predictor
        for(unsigned i=1;i<4;i++)
         {
          predicted_value += node_pt->x(i,j)*Predictor_weight[i];
         }
        //Store the predicted value
        //Note that predictor is stored as the FIFTH entry in this scheme
        node_pt->x(4,j) = predicted_value;
       }
     }
   }
 }
 
 
 //=======================================================================
 /// Function that sets the error weights (BDF2)
 //=======================================================================
 void SelfStartingBDF2::set_error_weights_bdf2()
 {
  if(adaptive_flag())
   {
    double dt=Time_pt->dt(0);
    double dtprev=Time_pt->dt(1);
    //Calculate the error weight
    Error_weight = pow((1.0 + dtprev/dt),2.0)/
     (1.0 + 3.0*(dtprev/dt) + 4.0*pow((dtprev/dt),2.0) + 
      2.0*pow((dtprev/dt),3.0));
   }
 }
 
 
 //===================================================================
 /// Function to compute the error in position i at node (BDF2)
 //===================================================================
 double SelfStartingBDF2::temporal_error_in_position_bdf2(
  Node* const &node_pt,const unsigned &i)
 {
  if(adaptive_flag())
   {
    //Just return the error
    return Error_weight*(node_pt->x(i) - node_pt->x(4,i));
   }
  else
   {
    return 0.0;
   }
 }
 
 
 //=========================================================================
 /// Function to calculate the error in the data value i (BDF2)
 //=========================================================================
 double SelfStartingBDF2::temporal_error_in_value_bdf2(Data* const &data_pt,
                                                       const unsigned &i)
 {
  if(adaptive_flag())
   {
    //Just return the error
    return Error_weight*(data_pt->value(i) - data_pt->value(4,i));
   }
  else
   {
    return 0.0;
   }
 }
 
 
} // End of namespace oomph
