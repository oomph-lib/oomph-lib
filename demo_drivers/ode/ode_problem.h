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
#ifndef OOMPH_ODE_PROBLEM_H
#define OOMPH_ODE_PROBLEM_H


#include "generic.h"
#include "ode.h"

#include "my_problem.h"


namespace oomph
{

 namespace VectorHelpers
 {
  inline double two_norm(const Vector<double>& a)
  {
   return std::sqrt(dot(a,a));
  }

  inline Vector<double> vector_diff(const Vector<double>& a, const Vector<double>& b)
  {
   check_lengths_match(a,b);
   unsigned ni = a.size();
   Vector<double> diff(ni, 0.0);
   for(unsigned i=0; i<ni; i++) {diff[i] = a[i] - b[i];}
   return diff;
  }
 }

 class ODEProblem : public MyProblem
 {
 public:

  /// constructor
  ODEProblem()
  {
   // Don't output to trace file every step, often too many steps
   Always_write_trace = false;

   Use_fd_jacobian = false;
  }

  virtual ~ODEProblem() {}

  virtual void build(Vector<Mesh*>& bulk_mesh_pt)
  {
   // Call the underlying build
   MyProblem::build(bulk_mesh_pt);

   // Set up the global mesh
   build_global_mesh();

   element_pt()->Use_fd_jacobian = Use_fd_jacobian;

   // assign equation numbers
   this->assign_eqn_numbers();
   oomph_info << "Number of equations: " << ndof() << std::endl;
  }

  void my_set_initial_condition(const SolutionFunctorBase& ic)
  {
   // Loop over current & previous timesteps
   const unsigned nprev_values = time_stepper_pt()->nprev_values();
   for(unsigned t=0; t<nprev_values+1; t++)
   {
    double time = time_pt()->time(t);

    std::cout << "setting IC at time =" << time << std::endl;

    // Get + set the (only) value
    Vector<double> dummy(nvalue(), 1.0);
    Vector<double> values = ic(time, dummy);

    for(unsigned j=0, nj=nvalue(); j<nj; j++)
    {
     mesh_pt()->element_pt(0)->internal_data_pt(0)
      ->set_value(t, j, values[j]);
    }
   }

   actions_after_set_initial_condition();
  }

  virtual void write_additional_trace_headers(std::ofstream& trace_file) const
  {
   trace_file << Trace_seperator << "exact";
  }

  virtual void write_additional_trace_data(const unsigned& t_hist,
                                           std::ofstream& trace_file) const
  {
   // Create a temporary vector to store the exact solution
   Vector<double> output_vector=exact_solution(time_pt()->time(t_hist));

   // Output the vector
   unsigned output_vector_length=output_vector.size();
   if (output_vector_length==0)
   {
    trace_file << Trace_seperator << "[]";
   }
   else
   {
    trace_file << Trace_seperator << "[" << output_vector[0];
    if (output_vector_length>1)
    {
     for (unsigned i=1;i<output_vector_length;i++)
     {
      trace_file << ", " << output_vector[i];
     }
    }
    trace_file << "]";
   }
  }

  virtual double get_error_norm(const unsigned& t_hist=0) const
  {
   Vector<double> val = trace_values(t_hist);
   Vector<double> exact = exact_solution(time_pt()->time(t_hist));

   return VectorHelpers::two_norm(VectorHelpers::vector_diff(val, exact));
  }

  Vector<double> exact_solution(const double& time) const
  {
   ODEElement* el_pt = checked_dynamic_cast<ODEElement*>
    (mesh_pt()->element_pt(0));
   return el_pt->exact_solution(time);
  }

  /// Error norm: use abs(error in data).
  double global_temporal_error_norm()
  {
   Data* dat_pt=mesh_pt()->element_pt(0)->internal_data_pt(0);

   return std::abs(ts_pt()->temporal_error_in_value(dat_pt, 0));
  }

  Vector<double> trace_values(const unsigned& t_hist=0) const
  {
   return solution(t_hist);
  }

  ODEElement* element_pt()
  {return checked_dynamic_cast<ODEElement*>(mesh_pt()->element_pt(0));}

  const ODEElement* element_pt() const
  {return checked_dynamic_cast<ODEElement*>(mesh_pt()->element_pt(0));}

  TimeStepper* ts_pt() const
  {
   return element_pt()->internal_data_pt(0)->time_stepper_pt();
  }

  unsigned nvalue() const
  {
   return element_pt()->nvalue();
  }

  // Output solution
  void output_solution(const unsigned& t, std::ostream& outstream,
		       const unsigned& npoints=2) const
  {
   // Create a temporary vector to store the exact solution
   Vector<double> output_vector=solution(t);

   // Vector length
   unsigned output_vector_length=output_vector.size();
   
   // Output the vector (to the command window)
   if (output_vector_length==0)
   {
    std::cout << Trace_seperator << "[]";
   }
   else
   {
    std::cout << Trace_seperator << "[" << output_vector[0];
    if (output_vector_length>1)
    {
     for (unsigned i=1;i<output_vector_length;i++)
     {
      std::cout << ", " << output_vector[i];
     }
    }
    std::cout << "]";
   }

   // Output the vector (to outstream)
   if (output_vector_length==0)
   {
    outstream << Trace_seperator << "[]";
   }
   else
   {
    outstream << Trace_seperator << "[" << output_vector[0];
    if (output_vector_length>1)
    {
     for (unsigned i=1;i<output_vector_length;i++)
     {
      outstream << ", " << output_vector[i];
     }
    }
    outstream << "]";
   }
      
   // npoints is ignored
  }

  Vector<double> solution(const unsigned& timestep=0) const
  {
   Data* dat_pt=mesh_pt()->element_pt(0)->internal_data_pt(0);
   Vector<double> solution(nvalue(), 0.0);
   dat_pt->value(timestep, solution);
   return solution;
  }

  bool Use_fd_jacobian;

 };

 namespace Factories
 {
  ///  Make a timestepper from an input argument using
  /// new.
  TimeStepper* time_stepper_factory(const std::string& ts_name)
  {
 
   // Always make timestepper adaptive, we can control adaptivity by
   // calling adaptive or non adaptive newton solve.
   bool adaptive_flag = true;

   if(ts_name == "bdf1")
   {
    return new BDF<1>(adaptive_flag);
   }
   else if(ts_name == "bdf2")
   {
    return new BDF<2>(adaptive_flag);
   }
   else if(ts_name == "real-imr")
   {
    IMR* mp_pt = new IMR(adaptive_flag);
    ExplicitTimeStepper* pred_pt = new EBDF3;
    mp_pt->set_predictor_pt(pred_pt);
    return mp_pt;
   }
   else if(ts_name == "imr")
   {
    IMRByBDF* mp_pt = new IMRByBDF(adaptive_flag);
    ExplicitTimeStepper* pred_pt = new EBDF3;
    mp_pt->set_predictor_pt(pred_pt);
    return mp_pt;
   }
   else if(ts_name == "steady")
   {
    // 2 steps so that we have enough space to do reasonable time
    // derivative estimates in e.g. energy derivatives.
    return new Steady<3>;
   }
   else if(ts_name == "tr")
   {
    // 2 steps so that we have enough space to do reasonable time
    // derivative estimates in e.g. energy derivatives.
    return new TR(adaptive_flag);
   }
   else
   {
    std::string err = "Unrecognised time stepper name";
    throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
			OOMPH_EXCEPTION_LOCATION);
   }
  }
 }


} // End of oomph namespace

#endif
