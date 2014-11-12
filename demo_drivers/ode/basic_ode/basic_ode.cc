
#include "generic.h"
#include "ode.h"

#include "../ode_problem.h"

#include "../ode_example_functions.h"

using namespace oomph;

namespace Factories
{
 /// \short Make a timestepper from an input argument using
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
  else if((ts_name == "midpoint") || (ts_name == "old-imr"))
   {
    MidpointMethod* mp_pt = new MidpointMethod(adaptive_flag);
    ExplicitTimeStepper* pred_pt = new EBDF3;
    mp_pt->set_predictor_pt(pred_pt);
    return mp_pt;
   }
  else if((ts_name == "midpoint-bdf") || (ts_name == "imr"))
   {
    MidpointMethodByBDF* mp_pt = new MidpointMethodByBDF(adaptive_flag);
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
  // else if(ts_name == "tr")
  //  {
  //   // 2 steps so that we have enough space to do reasonable time
  //   // derivative estimates in e.g. energy derivatives.
  //   return new TR(adaptive_flag);
  //  }
  else
   {
    std::string err = "Unrecognised time stepper name";
    throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
 }
}

int main()
{

 ODEProblem problem;


 problem.Exact_solution_pt = ODEFactories::exact_solutions_factory("damped_oscillation");


 TimeStepper* time_stepper_pt = Factories::time_stepper_factory("bdf2");
 double tmax = 5.0;
 double tol = 0;
 double dt = 0.1;


 Vector<Mesh*> mesh_pts;
 mesh_pts.push_back(new Mesh);
 mesh_pts[0]->add_element_pt(new ODEElement(time_stepper_pt, problem.Exact_solution_pt));
 problem.build(mesh_pts);

 problem.Doc_info.set_directory("Validation");



 // Set all dts to the value given in args
 problem.initialise_dt(dt);

 // Set values useing the initial condition function
 problem.set_initial_condition(*problem.Exact_solution_pt);

 problem.initial_doc();

 // Time step to end or to max number of steps
 while(problem.time() < tmax)
  {
   // Output some basic info
   oomph_info
    << std::endl
    << std::endl
    << "Time step " << problem.N_steps_taken
    << ", time = " << problem.time()
    << ", dt = " << dt << std::endl
    << "=============================================" << std::endl
    << std::endl;

   // Do the newton solve (different ones depending flags set)
   dt = problem.smart_time_step(dt, tol);

   // Output
   problem.doc_solution();
  }

 problem.final_doc(); 
 
 return 0;
}
