
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

int main(int argc, char *argv[])
{
 using namespace CommandLineArgs;

 // Parse cli arguments
 double dt = 0.1;
 specify_command_line_flag("-dt", &dt, "Time step size");

 double tmax = 5.0;
 specify_command_line_flag("-tmax", &tmax, "Time at which to stop");

 double tol = 0.0;
 specify_command_line_flag("-tol", &tol,
                           "Adaptive time step tolerance (default 0 = fixed step)");

 double dt_initial = 1e-5;
 specify_command_line_flag("-dt-initial", &dt_initial,
                           "Initial dt for adaptive time step selection.");

 std::string ts_name = "bdf2";
 specify_command_line_flag("-ts", &ts_name, "The time stepper to use.");

 std::string ode_name = "sin";
 specify_command_line_flag("-ode", &ode_name, "The ODE to solve.");

 std::string outdir = "results";
 specify_command_line_flag("-outdir", &outdir, "Directory to write output to.");

 setup(argc, argv);
 parse_and_assign(argc, argv, true);
 doc_all_flags();


 // Build problem
 ODEProblem problem;

 problem.Exact_solution_pt = ODEFactories::exact_solutions_factory(ode_name);

 TimeStepper* time_stepper_pt = Factories::time_stepper_factory(ts_name);


 Vector<Mesh*> mesh_pts;
 mesh_pts.push_back(new Mesh);
 mesh_pts[0]->add_element_pt(new ODEElement(time_stepper_pt, problem.Exact_solution_pt));
 problem.build(mesh_pts);

 problem.Doc_info.set_directory(outdir);



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
