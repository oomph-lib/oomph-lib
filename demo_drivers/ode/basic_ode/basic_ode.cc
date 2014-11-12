
#include "generic.h"
#include "ode.h"

#include "../ode_problem.h"
#include "../ode_element_for_imr.h"
#include "../ode_example_functions.h"

using namespace oomph;

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

 std::string element_type = "normal";
 specify_command_line_flag("-element-type", &element_type, "Element to use.");

 setup(argc, argv);
 parse_and_assign(argc, argv, true);
 doc_all_flags();


 // Build problem
 ODEProblem problem;

 problem.Exact_solution_pt = ODEFactories::exact_solutions_factory(ode_name);

 TimeStepper* time_stepper_pt = Factories::time_stepper_factory(ts_name);

 Vector<Mesh*> mesh_pts;
 mesh_pts.push_back(new Mesh);
 if(element_type == "normal")
  {
   mesh_pts[0]->add_element_pt(new ODEElement(time_stepper_pt, problem.Exact_solution_pt));
  }
 else if(element_type == "imr_element")
  {
   mesh_pts[0]->add_element_pt(new IMRODEElement(time_stepper_pt, problem.Exact_solution_pt));
  }
 else 
  {
   std::string err = "Unrecognised element type.";
   throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 problem.build(mesh_pts);

 problem.Doc_info.set_directory(outdir);



 // Set all dts to the value given in args
 problem.initialise_dt(dt);

 // Set values using the initial condition function (initialisation for
 // trapezoid rule is automatically handled here by the
 // actions_after_set_initial_conditions() function in MyProblem).
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

   // Do the newton solve (adaptive or not depending on tol)
   dt = problem.smart_time_step(dt, tol);

   // Output
   problem.doc_solution();
  }

 problem.final_doc(); 
 
 return 0;
}
