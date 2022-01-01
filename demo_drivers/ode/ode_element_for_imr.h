//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_ODE_ELEMENT_FOR_IMR_H
#define OOMPH_ODE_ELEMENT_FOR_IMR_H


#include "generic.h"
#include "ode.h"

namespace oomph
{
 /// Class for ODE elements when using the implicit midpoint rule time
 /// integrator IMR. The residual and Jacobian functions are
 /// modified to interpolate in time. Note that the normal element can be
 /// used for the IMR implementation in the class IMRByBDF.
 class IMRODEElement : public ODEElement 
 {
 public:
  /// Constructor
  IMRODEElement(TimeStepper* time_stepper_pt,
                SolutionFunctorBase* exact_solution_pt)
   : ODEElement(time_stepper_pt, exact_solution_pt) {}
  
  /// Virtual destructor
  virtual ~IMRODEElement() {}

  /// For IMR implemented using IMR we have to interpolate
  /// everything in time (to get the midpoint values). This function does
  /// that for the u values.
  Vector<double> time_interpolate_u() const
  {

   // Get pointer to one-and-only internal data object
   Data* dat_pt = internal_data_pt(0);

   // Get time stepper
   TimeStepper* ts_pt = dat_pt->time_stepper_pt();
   
   // Number of history values to interpolate over for value/time calculation
   const unsigned nj = ts_pt->nprev_values_for_value_at_evaluation_time();

   Vector<double> u(nvalue(), 0.0);
   for(unsigned j=0; j<nj; j++)
    {
     for(unsigned i=0; i<nvalue(); i++)
      {
       u[i] += dat_pt->value(j, i) * ts_pt->weight(0, j);
      }
    }

   return u;   
  }

  /// For IMR implemented using IMR we have to interpolate
  /// everything in time (to get the midpoint values). This function does
  /// that for the continuous time.
  double time_interpolate_time() const
  {

   // Get pointer to one-and-only internal data object
   Data* dat_pt = internal_data_pt(0);

   // Get time stepper
   TimeStepper* ts_pt = dat_pt->time_stepper_pt();

   // Number of history values to interpolate over for value/time calculation
   const unsigned nj = ts_pt->nprev_values_for_value_at_evaluation_time();

   double t = 0;
   for(unsigned j=0; j<nj; j++)
    {
     t += ts_pt->time_pt()->time(j) * ts_pt->weight(0, j);
    }

   return t;
  }
  
  /// Modified get residuals: have to interpolate in time for IMR when
  /// using the IMR implementation so we can't use the normal
  /// ODE element implementation.
  void fill_in_contribution_to_residuals(Vector<double>& residuals)
   {

    // Get pointer to one-and-only internal data object
    Data* dat_pt = internal_data_pt(0);

    // Get time stepper
    TimeStepper* ts_pt = dat_pt->time_stepper_pt();

    // Get values **interpolated in time**
    Vector<double> u = time_interpolate_u();

    // Get **interpolated** continuous time
    double t = time_interpolate_time();

    Vector<double> deriv = derivative_function(t, u);
    for(unsigned j=0, nj=deriv.size(); j<nj; j++)
     {
      // Get dudt approximation from time stepper
      double dudt = ts_pt->time_derivative(1, dat_pt, j);

      // Residual is difference between the exact derivative and our
      // time stepper's derivative estimate.
      residuals[j] = deriv[j] - dudt;
     }
   }

  /// Modified get jacobian: two differences from the normal verson here
  /// which are required when using the IMR implementation: 
  /// 1) We have to interpolate in time to calculate the analytical
  /// contribution. 
  /// 2) The interpolation in time within the residual calculation modifies
  /// the Jacobian (because the jacobian is the derivative w.r.t. the
  /// n+1-th values).
  void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                        DenseMatrix<double>& jacobian)
  {
   // Get residuals
   fill_in_contribution_to_residuals(residuals);

   if(Exact_solution_pt->have_jacobian() && !Use_fd_jacobian)
    {
     const unsigned n = nvalue();
     TimeStepper* ts_pt = internal_data_pt(0)->time_stepper_pt();

     // Get values **interpolated in time**
     Vector<double> u = time_interpolate_u();

     // Get **interpolated** continuous time
     double t = time_interpolate_time();

     // get df/du jacobian analytically
     Vector<double> dummy;
     Exact_solution_pt->jacobian(t, dummy, u, jacobian);

     // Multiply by the weight of the history value 0 (i.e. at time n+1) in
     // the calculation of the time-interpolated "current" value (0th
     // derivative). Need to do this because of the chain rule, just write
     // out the derivation of dr/du for nvalue=1 to see it.
     const double w00 = ts_pt->weight(0,0);
     for(unsigned i=0; i<n; i++)
      {
       for(unsigned j=0; j<n; j++)
        {
         jacobian(i, j) *= w00;
        }
      }

     // We actually need jacobian of residual = f(imr_t, imr_u) - dudt so
     // subtract diagonal (dudt)/du term now.
     const double w10 = ts_pt->weight(1,0);
     for(unsigned i=0; i<n; i++)
      {
       jacobian(i, i) -= w10;
      }
    }
   else
    {
     // Use FD for jacobian
     GeneralisedElement::fill_in_jacobian_from_internal_by_fd
      (residuals, jacobian, true);
    }
   }

  // Mass matrix calculation should be unaffected by the use of IMR.
  
 private:
  /// Broken copy constructor
  IMRODEElement(const IMRODEElement& dummy)
  {BrokenCopy::broken_copy("IMRODEElement");}
  
 };

} // End of oomph namespace

#endif
