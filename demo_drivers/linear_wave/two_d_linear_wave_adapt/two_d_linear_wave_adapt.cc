//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Driver for the adaptive solution of the 2D wave equation 
// in a rectanglular domain 

//Generic includes
#include "generic.h"

//Linear wave includes
#include "linear_wave.h"

// Mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;

//==start_of_tanh_solution============================================
/// Namespace for exact solution for LinearWave equation 
/// with sharp step 
//====================================================================
namespace TanhSolnForLinearWave
{

 /// Parameter for steepness of step
 double Alpha;

 /// Orientation of step wave
 double Phi;

 /// Exact solution
 double exact_u(const double& time, const Vector<double>& x)
 {
  double zeta=cos(Phi)*x[0]+sin(Phi)*x[1];
  return tanh(1.0-Alpha*(zeta-time));
 }
 
 /// 1st time-deriv of exact solution
 double exact_dudt(const double& time, const Vector<double>& x)
 {
  double zeta=cos(Phi)*x[0]+sin(Phi)*x[1];
  return Alpha/(cosh(1.0-Alpha*(zeta-time))*
              cosh(1.0-Alpha*(zeta-time)));
 }

 /// 2nd time-deriv of exact solution
 double exact_d2udt2(const double& time, const Vector<double>& x)
 {
  double zeta=cos(Phi)*x[0]+sin(Phi)*x[1];
  return -2.0*Alpha*Alpha*tanh(1.0-Alpha*(zeta-time))/
   (cosh(1.0-Alpha*(zeta-time))*cosh(1.0-Alpha*(zeta-time)));
 }


 /// Exact solution as a vector
 void get_exact_u(const double& time, const Vector<double>& x, 
                  Vector<double>& u)
 {
  u[0]=exact_u(time,x);
  u[1]=exact_dudt(time,x);
  u[2]=exact_d2udt2(time,x);
 }

 /// Source function to make it an exact solution 
 void get_source(const double& time, const Vector<double>& x, double& source)
 {
  source=0.0;
 }

} // end of tanh solution



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//===start_of_problem_class===========================================
/// LinearWave problem in rectanglular domain
//====================================================================
template<class ELEMENT, class TIMESTEPPER>
class LinearWaveProblem : public Problem
{

public:

 /// Constructor: pass number of elements in x and y directions
 /// bool indicating impulsive or "smooth" start 
 /// and pointer to source function
 LinearWaveProblem(const unsigned& nx, const unsigned& ny, 
                   const bool& impulsive_start,
                   LinearWaveEquations<2>::LinearWaveSourceFctPt source_fct_pt);

 /// Destructor (empty) 
 ~LinearWaveProblem() {}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem specs after solve (empty)
 void actions_after_implicit_timestep() {}

 /// Update the problem specs before next timestep: 
 /// Set time-dependent Dirchlet boundary from exact solution.
 void actions_before_implicit_timestep()
  {
   Vector<typename TIMESTEPPER::NodeInitialConditionFctPt> 
    initial_value_fct(1);
   Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
    initial_veloc_fct(1);
   Vector<typename TIMESTEPPER::NodeInitialConditionFctPt> 
    initial_accel_fct(1);
   
   // Assign values for analytical value, veloc and accel:
   initial_value_fct[0]=&TanhSolnForLinearWave::exact_u;
   initial_veloc_fct[0]=&TanhSolnForLinearWave::exact_dudt;
   initial_accel_fct[0]=&TanhSolnForLinearWave::exact_d2udt2;
   
   // Loop over boundaries
   unsigned num_bound=mesh_pt()->nboundary();
   for (unsigned ibound=0;ibound<num_bound;ibound++)
    {
     // Loop over boundary nodes
     unsigned num_nod=mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Set the boundary condition from the exact solution
       Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

       bool use_direct_assignment=false;
       if (use_direct_assignment)
        {
         // Set nodal coordinates for evaluation of BC:
         Vector<double> x(2);
         x[0]=nod_pt->x(0);
         x[1]=nod_pt->x(1);
         
         // Set exact solution at current time
         nod_pt->
          set_value(0,
                    TanhSolnForLinearWave::exact_u(time_pt()->time(),x));
        }
       else
        {  
         // Get timestepper
         TIMESTEPPER* timestepper_pt=dynamic_cast<TIMESTEPPER*>
          (time_stepper_pt());
         
         // Assign the history values
         timestepper_pt->assign_initial_data_values(nod_pt, 
                                                    initial_value_fct,
                                                    initial_veloc_fct,
                                                    initial_accel_fct);
        }
      }
    }
  } // end of actions before timestep


 /// Access function for the specific mesh
 RefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Set initial condition (incl history values)
 void set_initial_condition();

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Do unsteady run 
 void unsteady_run();

private:

 // Trace file
 ofstream Trace_file;

 // Impulsive start?
 bool Impulsive_start;

}; // end of problem class



//===start_of_constructor=================================================
/// Constructor for LinearWave problem 
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
LinearWaveProblem<ELEMENT,TIMESTEPPER>::LinearWaveProblem(
 const unsigned& nx, const unsigned& ny, const bool& impulsive_start, 
 LinearWaveEquations<2>::LinearWaveSourceFctPt source_fct_pt) :
 Impulsive_start(impulsive_start)
{ 

 //Allocate the timestepper -- this constructs the time object as well
 add_time_stepper_pt(new TIMESTEPPER());

 // Set up parameters for exact solution
 //-------------------------------------

 // Steepness of tanh profile
 TanhSolnForLinearWave::Alpha=10.0;

 // Orientation of step wave
 TanhSolnForLinearWave::Phi=MathematicalConstants::Pi/180.0*30.0;

 // Set up mesh
 //------------

 // # of elements in x-direction
 unsigned Nx=nx;

 // # of elements in y-direction
 unsigned Ny=ny;

 // Domain length in x-direction
 double Lx=1.0;

 // Domain length in y-direction
 double Ly=2.0;

 // Build and assign mesh
 Problem::mesh_pt()=new RefineableRectangularQuadMesh<ELEMENT>(
  Nx,Ny,Lx,Ly,time_stepper_pt());

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
    }
  } //end of boundary conditions

 // Complete build of elements
 // --------------------------

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor 
unsigned n_element = mesh_pt()->nelement();
for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = source_fct_pt;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor




//===start_of_set_initial_condition=======================================
/// Set initial condition.
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void LinearWaveProblem<ELEMENT,TIMESTEPPER>::set_initial_condition()
{ 

 // Get timestepper
 TIMESTEPPER* timestepper_pt=dynamic_cast<TIMESTEPPER*>(time_stepper_pt());


 // Impulsive start
 //----------------
 if (Impulsive_start)
  {
   // Loop over the nodes to set initial conditions everywhere
   unsigned num_nod=mesh_pt()->nnode();
   for (unsigned jnod=0;jnod<num_nod;jnod++)
    {
     // Pointer to node
     Node* nod_pt=mesh_pt()->node_pt(jnod);

     // Get nodal coordinates
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     // Assign initial value from exact solution
     nod_pt->set_value(0,TanhSolnForLinearWave::exact_u(time_pt()->time(),x));

     // Set history values so that they are consistent with an impulsive
     // start from this value
     timestepper_pt->assign_initial_values_impulsive(nod_pt);
    }
  } // end impulsive start

 // "Smooth" start from analytical time history
 //--------------------------------------------
 else
  {

   // Vector of function pointers to functions that specify the
   // value, and the first and second time-derivatives of the
   // function used as the initial condition
   Vector<typename TIMESTEPPER::NodeInitialConditionFctPt> 
    initial_value_fct(1);
   Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
    initial_veloc_fct(1);
   Vector<typename TIMESTEPPER::NodeInitialConditionFctPt> 
    initial_accel_fct(1);
   
   // Assign values for analytical value, veloc and accel:
   initial_value_fct[0]=&TanhSolnForLinearWave::exact_u;
   initial_veloc_fct[0]=&TanhSolnForLinearWave::exact_dudt;
   initial_accel_fct[0]=&TanhSolnForLinearWave::exact_d2udt2;
   
   // Assign Newmark history values so that Newmark approximations
   // for velocity and accel are correct at initial time:

   // Loop over the nodes to set initial conditions everywhere
   unsigned num_nod=mesh_pt()->nnode();
   for (unsigned jnod=0;jnod<num_nod;jnod++)
    {
     // Pointer to node
     Node* nod_pt=mesh_pt()->node_pt(jnod);
    
     // Assign the history values
     timestepper_pt->assign_initial_data_values(nod_pt, 
                                                initial_value_fct,
                                                initial_veloc_fct,
                                                initial_accel_fct);
    } // end of smooth start


   // Paranoia: Check that the initial values were assigned correctly
   double err_max=0.0;
   for (unsigned jnod=0;jnod<num_nod;jnod++)
    {
     // Pointer to node
     Node* nod_pt=mesh_pt()->node_pt(jnod);

     // Get nodal coordinates
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     // Get exact value and first and second time-derivatives
     double u_exact=
      TanhSolnForLinearWave::exact_u(time_pt()->time(),x);
     double dudt_exact=
      TanhSolnForLinearWave::exact_dudt(time_pt()->time(),x);
     double d2udt2_exact=
      TanhSolnForLinearWave::exact_d2udt2(time_pt()->time(),x);
    
     // Get Newmark approximations for zero-th, first and second 
     // time-derivatives of the nodal values. 
     double u_fe=timestepper_pt->time_derivative(0,nod_pt,0);
     double dudt_fe=timestepper_pt->time_derivative(1,nod_pt,0);
     double d2udt2_fe=timestepper_pt->time_derivative(2,nod_pt,0);
     
     // Error
     double error=sqrt(pow(u_exact-u_fe,2)+
                       pow(dudt_exact-dudt_fe,2)+
                       pow(d2udt2_exact-d2udt2_fe,2));
     if (error>err_max) err_max=error;
    }
   cout << "Max. error in assignment of initial condition " 
        << err_max << std::endl;   
  }


} // end of set initial condition



//===start_of_doc_solution================================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void LinearWaveProblem<ELEMENT,TIMESTEPPER>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 cout << std::endl;
 cout << "=================================================" << std::endl;
 cout << "Docing solution for t=" << time_pt()->time() << std::endl;
 cout << "=================================================" << std::endl;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
           << time_pt()->time() << "\"";
 some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << std::endl;
 some_file << "1" << std::endl;
 some_file << "2" << std::endl;
 some_file << " 0 0" << std::endl;
 some_file << time_pt()->time()*20.0 << " 0" << std::endl;

 some_file.close();

 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,time_pt()->time(),
                       TanhSolnForLinearWave::get_exact_u); 
 some_file.close();

 // Doc error
 //----------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,
                          TanhSolnForLinearWave::get_exact_u,
                          time_pt()->time(),
                          error,norm); 
 some_file.close();
 cout << "error: " << error << std::endl; 
 cout << "norm : " << norm << std::endl << std::endl;

 // Write trace file
 Trace_file << time_pt()->time() << " " << time_pt()->dt()
            << " " << mesh_pt()->nelement() << " " 
            << error << " " << norm << std::endl;

} // end of doc solution

 


//===start_of_unsteady_run================================================
/// Unsteady run.
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void LinearWaveProblem<ELEMENT,TIMESTEPPER>::unsteady_run()
{

 // Setup labels for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 // Initialise time
 double time0=0.0;
 time_pt()->time()=time0;

 // Set initial timestep
 double dt=0.005;
 time_pt()->initialise_dt(dt);

 // Set IC
 set_initial_condition();

 //Output initial condition
 doc_solution(doc_info);
 
 //Increment counter for solutions 
 doc_info.number()++;

 // Maximum time
 double t_max=4.0;

 // Number of steps
 unsigned nstep=unsigned(t_max/dt);

 // If validation run only do 2 timesteps
 if (CommandLineArgs::Argc>1)
  { 
   nstep=2; 
   cout << "Validation run -- only doing two timesteps." << std::endl;
  }
 

 // First timestep: Larger number of maximum adaptations
 unsigned max_adapt=10;
 bool first=true;
 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   //Take fixed timestep with spatial adaptivity
   unsteady_newton_solve(dt,max_adapt,first);
      
   // Now we've done the first step
   max_adapt=1;
   first=false;

   //Output solution
   doc_solution(doc_info);
     
   //Increment counter for solutions 
   doc_info.number()++;
  }

 // Close trace file
 Trace_file.close();

} // end of unsteady run




/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//===start_of_main========================================================
/// Demonstrate how to solve LinearWave problem.
//========================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments: If a command line argument is specied
 // we regard this as validation run.
 CommandLineArgs::setup(argc,argv);
 
 // Number of elements in x direction
 unsigned n_x=10;

 // Number of elements in y direction
 unsigned n_y=20;

 // Impulsive start?
 bool impulsive_start=false;
 
 // Build problem
 LinearWaveProblem<RefineableQLinearWaveElement<2,3>, Newmark<1> >
  problem(n_x,n_y,impulsive_start,&TanhSolnForLinearWave::get_source);
 
 // Run it
 problem.unsteady_run();


}; // end of main
