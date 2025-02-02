//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
// Driver for Rayleigh-type problem: 2D channel whose upper
// wall oscillates periodically.

// The oomphlib headers
#include "generic.h"
#include "navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//===start_of_namespace=================================================
/// Namespace for global parameters
//======================================================================
namespace Global_Parameters
{
 /// Reynolds number
 double Re;

 /// Womersley = Reynolds times Strouhal
 double ReSt;

 /// Flag for long/short run: Default =  perform long run
 unsigned Long_run_flag=1;

 /// Flag for impulsive start: Default = start from exact
 /// time-periodic solution. 
 unsigned Impulsive_start_flag=0;

} // end of namespace


//==start_of_exact_solution=============================================
/// Namespace for exact solution
//======================================================================
namespace ExactSoln
{
 
 /// Exact solution of the problem as a vector
 void get_exact_u(const double& t, const Vector<double>& x, Vector<double>& u)
 {
  double y=x[1];
  // I=sqrt(-1)
  complex<double> I(0.0,1.0);
  // omega
  double omega=2.0*MathematicalConstants::Pi;
  // lambda
  complex<double> lambda(0.0,omega*Global_Parameters::ReSt);
  lambda = I*sqrt(lambda);

  // Solution
  complex<double> sol(
   exp(complex<double>(0.0,omega*t)) * 
   (exp(lambda*complex<double>(0.0,y))-exp(lambda*complex<double>(0.0,-y)))
   /(exp(I*lambda)-exp(-I*lambda)) );
  
  // Assign real solution
  u.resize(2);
  u[0]=real(sol);
  u[1]=0.0;
 }

 /// Exact solution of the problem as a scalar
 void get_exact_u(const double& t, const double& y,double& u)
 {
  // I=sqrt(-1)
  complex<double> I(0.0,1.0);
  // omega
  double omega=2.0*MathematicalConstants::Pi;
  // lambda
  complex<double> lambda(0.0,omega*Global_Parameters::ReSt);
  lambda = I*sqrt(lambda);
  // Solution
  complex<double> sol(
   exp(complex<double>(0.0,omega*t)) * 
   (exp(lambda*complex<double>(0.0,y))-exp(lambda*complex<double>(0.0,-y)))
   /(exp(I*lambda)-exp(-I*lambda)) );

  // Assign real solution
  u=real(sol);

 }

} // end of exact_solution


//===start_of_problem_class=============================================
/// Rayleigh-type problem: 2D channel whose upper
/// wall oscillates periodically.
//======================================================================
template<class ELEMENT, class TIMESTEPPER>
class RayleighProblem : public Problem
{
public:

 /// Constructor: Pass number of elements in x and y directions and 
 /// lengths and ALE flag
 RayleighProblem(const unsigned &nx, const unsigned &ny, 
                 const double &lx, const double &ly,
                 const bool& use_ale);

 //Update before solve is empty
 void actions_before_newton_solve() {}

 /// Update after solve is empty
 void actions_after_newton_solve() {}

 //Actions before timestep: Update no slip on upper oscillating wall
 void actions_before_implicit_timestep()
  {
   // No slip on upper boundary
   unsigned ibound=2;
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get exact solution
     double y=mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
     double time=time_pt()->time();
     double soln;
     ExactSoln::get_exact_u(time,y,soln);
     
     // Assign exact solution to boundary
     mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,soln);
     mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
    }

  } // end of actions_before_implicit_timestep
   
 /// Run an unsteady simulation
 void unsteady_run(DocInfo& doc_info); 
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Set initial condition (incl previous timesteps) according
 /// to specified function. 
 void set_initial_condition();

private:

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to proper element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  }

 /// Trace file
 ofstream Trace_file;

 /// Flag for ALE
 bool Use_ALE;

}; // end of problem_class


//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT,class TIMESTEPPER>
RayleighProblem<ELEMENT,TIMESTEPPER>::RayleighProblem
(const unsigned &nx, const unsigned &ny,
 const double &lx, const double& ly,
 const bool& use_ale) : Use_ALE(use_ale)
{
 //Allocate the timestepper
 add_time_stepper_pt(new TIMESTEPPER); 

 //Now create the mesh with periodic boundary conditions in x direction
 bool periodic_in_x=true;
 Problem::mesh_pt() = 
  new RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,periodic_in_x,
                                   time_stepper_pt());

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here
 unsigned num_bound=mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // No slip on top and bottom
     if ((ibound==0)||(ibound==2))
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
     // Horizontal outflow on the left (and right -- right bc not
     // strictly necessary because of symmetry)
     else if ((ibound==1)||(ibound==3))
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
    }
  } // end loop over boundaries
 
 if (Use_ALE)
  {
   std::cout << "Enabling ALE" << std::endl;
  }
 else
  {
   std::cout << "Disabling ALE" << std::endl;
  }
 

 //Complete the problem setup to make the elements fully functional

 //Loop over the elements
 unsigned n_el = mesh_pt()->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Parameters::Re;
   el_pt->re_st_pt() = &Global_Parameters::ReSt;

   // Enable/disable ALE?
   if (Use_ALE)
    {
     el_pt->enable_ALE();
    }
   else
    {
     el_pt->disable_ALE();
    }
  }

 // Now pin the pressure in first element at value 0 to 0.0
 fix_pressure(0,0,0.0);

 //Assgn equation numbers
 cout << assign_eqn_numbers() << std::endl; 
} // end of constructor




//======================start_of_set_initial_condition====================
/// Set initial condition: Assign previous and current values
/// from exact solution.
//========================================================================
template<class ELEMENT,class TIMESTEPPER>
void RayleighProblem<ELEMENT,TIMESTEPPER>::set_initial_condition()
{ 

 // Check that timestepper is from the BDF family
 if (time_stepper_pt()->type()!="BDF")
  {
   std::ostringstream error_stream;
   error_stream << "Timestepper has to be from the BDF family!\n"
                << "You have specified a timestepper from the "
                << time_stepper_pt()->type() << " family" << std::endl;

   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Backup time in global Time object
 double backed_up_time=time_pt()->time();
         
 // Past history needs to be established for t=time0-deltat, ...
 // Then provide current values (at t=time0) which will also form
 // the initial guess for the first solve at t=time0+deltat
 
 // Vector of exact solution value
 Vector<double> soln(2);
 Vector<double> x(2);

 //Find number of nodes in mesh
 unsigned num_nod = mesh_pt()->nnode();

 // Set continuous times at previous timesteps:
 // How many previous timesteps does the timestepper use?
 int nprev_steps=time_stepper_pt()->nprev_values();
 Vector<double> prev_time(nprev_steps+1);
 for (int t=nprev_steps;t>=0;t--)
  {
   prev_time[t]=time_pt()->time(unsigned(t));
  } 

 // Loop over current & previous timesteps
 for (int t=nprev_steps;t>=0;t--)
  {
   // Continuous time
   double time=prev_time[t];
   cout << "setting IC at time =" << time << std::endl;
   
   // Loop over the nodes to set initial guess everywhere
   for (unsigned n=0;n<num_nod;n++)
    {
     // Get nodal coordinates
     x[0]=mesh_pt()->node_pt(n)->x(0);
     x[1]=mesh_pt()->node_pt(n)->x(1);

     // Get exact solution at previous time
     ExactSoln::get_exact_u(time,x,soln);
     
     // Assign solution
     mesh_pt()->node_pt(n)->set_value(t,0,soln[0]);
     mesh_pt()->node_pt(n)->set_value(t,1,soln[1]);
     
     // Loop over coordinate directions: Mesh doesn't move, so 
     // previous position = present position
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->node_pt(n)->x(t,i)=x[i];
      }
    } 
  }

 // Reset backed up time for global timestepper
 time_pt()->time()=backed_up_time;

} // end of set_initial_condition


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT,class TIMESTEPPER>
void RayleighProblem<ELEMENT,TIMESTEPPER>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);

 // Write file as a tecplot text object
 some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
           << time_pt()->time() << "\"";
 // ...and draw a horizontal line whose length is proportional
 // to the elapsed time
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
                       ExactSoln::get_exact_u); 
 some_file.close();

 // Doc error
 //----------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,
                          ExactSoln::get_exact_u,
                          time_pt()->time(),
                          error,norm); 
 some_file.close();

 // Doc solution and error
 //-----------------------
 cout << "error: " << error << std::endl; 
 cout << "norm : " << norm << std::endl << std::endl;

 // Get time, position and exact soln at control node 
 unsigned n_control=37;
 Vector<double> x(2), u(2);
 double time=time_pt()->time();
 Node* node_pt=
  dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(n_control))->node_pt(1);
 x[0] = node_pt->x(0);
 x[1] = node_pt->x(1);
 ExactSoln::get_exact_u(time,x,u);

 // Write trace file
 Trace_file << time << " " 
            << x[0] << " "
            << x[1] << " "
            << node_pt->value(0) << " "
            << node_pt->value(1) << " "
            << u[0] << " "
            << u[1] << " " 
            << abs(u[0]-node_pt->value(0)) << " "
            << abs(u[1]-node_pt->value(1)) << " " 
            << error << " "
            << norm << " "
            << std::endl;


} // end_of_doc_solution   


//===start_of_unsteady_run=====================================================
/// Unsteady run...
//=============================================================================
template<class ELEMENT,class TIMESTEPPER>
void RayleighProblem<ELEMENT,TIMESTEPPER>::unsteady_run(DocInfo& doc_info)
{

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 // Write tecplot header for trace file
 Trace_file << "time" << ",      " 
            << "x" << ",         "
            << "y" << ",         "
            << "u_1" << ",       "
            << "u_2" << ",       "
            << "u_exact_1" << ", "
            << "u_exact_2" << ", "
            << "error_1" << ",   "
            << "error_2" << ",   " 
            << "L2 error" << ",  "
            << "L2 norm" << ",   " << std::endl;

 //Set value of dt
 double dt = 0.025;

 if (Global_Parameters::Impulsive_start_flag==1)
  {
   // Initialise all history values for an impulsive start
   assign_initial_values_impulsive(dt);
   cout << "IC = impulsive start" << std::endl;
  }
 else
  {
   // Initialise timestep
   initialise_dt(dt);
   // Set initial conditions.
   set_initial_condition();
   cout << "IC = exact solution" << std::endl;
  } 

 //Now do many timesteps
 unsigned ntsteps=80;

 // If validation run only do 5 timesteps
 if (Global_Parameters::Long_run_flag==0)
  { 
   ntsteps=5; 
   cout << "validation run" << std::endl;
  }

 // Doc initial condition
 doc_solution(doc_info);
 
 // increment counter
 doc_info.number()++;
                                                           
 // Initialise timer
 clock_t t_start = clock();

 //Loop over the timesteps
 for(unsigned t=1;t<=ntsteps;t++)
  {
   cout << "TIMESTEP " << t << std::endl;
   
   //Take one fixed timestep
   unsteady_newton_solve(dt);

   //Output the time
   //cout << "Time is now " << time_pt()->time() << std::endl;
   
   // Doc solution
   // doc_solution(doc_info);

   // increment counter
   // doc_info.number()++;
  }

 // Stop timer
 clock_t t_end = clock();
 
 std::cout << "Total time for solution: "
           << double(t_end-t_start)/CLOCKS_PER_SEC
           << std::endl; 
 
 //Output final solution for check
 doc_solution(doc_info);

} // end of unsteady run


//===start_of_main======================================================
/// Driver code for Rayleigh channel problem
//======================================================================
int main(int argc, char* argv[]) 
{

 /// Convert command line arguments (if any) into flags:
 if (argc==1)
  {
   cout << "No command line arguments specified -- using defaults." << std::endl;
  }
 else if (argc==3)
  {
   cout << "Two command line arguments specified:" << std::endl;
   // Flag for long run
   Global_Parameters::Long_run_flag=atoi(argv[1]);
   /// Flag for impulsive start
   Global_Parameters::Impulsive_start_flag=atoi(argv[2]);
  }
 else
  {
   std::string error_message = 
    "Wrong number of command line arguments. Specify none or two.\n";
   error_message +=
    "Arg1: Long_run_flag [0/1]\n";
   error_message +=
    "Arg2: Impulsive_start_flag [0/1]\n";

   throw OomphLibError(error_message,
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 cout << "Long run flag: " 
      <<  Global_Parameters::Long_run_flag << std::endl;
 cout << "Impulsive start flag: " 
      <<  Global_Parameters::Impulsive_start_flag << std::endl;


 // Set physical parameters:

 // Womersley number = Reynolds number (St = 1)
 Global_Parameters::ReSt = 10.0;
 Global_Parameters::Re = Global_Parameters::ReSt;

 //Horizontal length of domain
 double lx = 1.0;

 //Vertical length of domain
 double ly = 1.0;

 // Number of elements in x-direction
 unsigned nx=5;

 // Number of elements in y-direction
 unsigned ny=10;

 // Loop over run with/without ALE
 for (unsigned ale_flag=0;ale_flag<2;ale_flag++)
  {

   bool use_ale=true;
   if (ale_flag==1) use_ale=false;

   // Solve with Crouzeix-Raviart elements
   {
    // Set up doc info
    DocInfo doc_info;
    doc_info.number()=0;
    if (ale_flag==0)
     {
      doc_info.set_directory("RESLT_CR");
     }
    else
     {
      doc_info.set_directory("RESLT_CR_ALE");
     }

    //Set up problem
    RayleighProblem<QCrouzeixRaviartElement<2>,BDF<2> > 
     problem(nx,ny,lx,ly,use_ale);
    
    // Run the unsteady simulation
    problem.unsteady_run(doc_info);
   }
   
   // Solve with Taylor-Hood elements
   {
    // Set up doc info
    DocInfo doc_info;
    doc_info.number()=0;
    if (ale_flag==0)
     {
      doc_info.set_directory("RESLT_TH");
     }
    else
     {
      doc_info.set_directory("RESLT_TH_ALE");
     }
    
    //Set up problem
    RayleighProblem<QTaylorHoodElement<2>,BDF<2> >
     problem(nx,ny,lx,ly,use_ale);
    
    // Run the unsteady simulation
    problem.unsteady_run(doc_info);
   }
  }

} // end of main
