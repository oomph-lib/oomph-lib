//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
// Driver for 2D unsteady heat problem with temporal adaptivity and
// optional restart from disk

//Generic routines
#include "generic.h"

// The unsteady heat equations
#include "unsteady_heat.h"

// Mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;


/////////////////////////////////////////////////////////////////////// 
///////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// 

//======start_of_ExactSolnForUnsteadyHeat=====================
/// Namespace for forced exact solution for UnsteadyHeat equation 
//====================================================================
namespace ExactSolnForUnsteadyHeat
{

 /// Factor controlling the rate of change
 double Gamma=10.0;

 /// Wavenumber
 double K=3.0;

 /// Angle of bump
 double Phi=1.0; 
 
 /// Exact solution as a Vector
 void get_exact_u(const double& time, const Vector<double>& x, 
                  Vector<double>& u)
 {
  double zeta=cos(Phi)*x[0]+sin(Phi)*x[1];
  u[0]=sin(K*zeta)*
       0.5*(1.0+tanh(Gamma*cos(2.0*MathematicalConstants::Pi*time)));
 }
 
 /// Exact solution as a scalar
 void get_exact_u(const double& time, const Vector<double>& x, double& u)
 {
  double zeta=cos(Phi)*x[0]+sin(Phi)*x[1];
  u=sin(K*zeta)*
       0.5*(1.0+tanh(Gamma*cos(2.0*MathematicalConstants::Pi*time)));
 }
 
 /// Source function to make it an exact solution 
 void get_source(const double& time, const Vector<double>& x, double& source)
 {
  source=
   -0.5*sin(K*(cos(Phi)*x[0]+sin(Phi)*x[1]))*K*K*pow(cos(Phi),2.0)*(
    0.1E1+tanh(Gamma*cos(0.2E1*0.3141592653589793E1*time)))-
   0.5*sin(K*(cos(Phi)*x[0]+sin(Phi)*x[1]))*K*K*pow(sin(Phi),2.0)*
   (0.1E1+tanh(Gamma*cos(0.2E1*0.3141592653589793E1*time)))+
   0.1E1*sin(K*(cos(Phi)*x[0]+sin(Phi)*x[1]))*
   (1.0-pow(tanh(Gamma*cos(0.2E1*0.3141592653589793E1*time)),2.0))*
   Gamma*sin(0.2E1*0.3141592653589793E1*time)*0.3141592653589793E1;
 }

} // end of ExactSolnForUnsteadyHeat

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//=====start_of_problem_class=========================================
/// UnsteadyHeat problem
//====================================================================
template<class ELEMENT>
class UnsteadyHeatProblem : public Problem
{

public:

 /// Constructor
 UnsteadyHeatProblem(UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt 
 source_fct_pt);

 /// Destructor (empty)
 ~UnsteadyHeatProblem(){}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// \short Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem specs after solve (empty)
 void actions_after_implicit_timestep() {}

 /// \short Update the problem specs before next timestep: 
 /// Set Dirchlet boundary conditions from exact solution.
 void actions_before_implicit_timestep();

 /// \short Set initial condition (incl previous timesteps) according
 /// to specified function. 
 void set_initial_condition();

 /// Doc the solution
 void doc_solution(DocInfo& doc_info, ofstream& trace_file);
 
 /// \short Dump problem to disk to allow for restart.
 void dump_it(ofstream& dump_file);

 /// \short Read problem for restart from specified restart file.
 void restart(ifstream& restart_file);

 /// Global error norm for adaptive time-stepping
 double global_temporal_error_norm();

private:

 /// Pointer to source function
 UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt Source_fct_pt;
 
 /// Pointer to control node at which the solution is documented
 Node* Control_node_pt;

}; // end of problem class




//========start_of_constructor============================================
/// Constructor for UnsteadyHeat problem in square domain
//========================================================================
template<class ELEMENT>
UnsteadyHeatProblem<ELEMENT>::UnsteadyHeatProblem(
 UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt source_fct_pt) : 
 Source_fct_pt(source_fct_pt)
{ 
 
 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 add_time_stepper_pt(new BDF<2>(true));

 // Setup parameters for exact solution
 // -----------------------------------

 // Parameter controlling the rate of change
 ExactSolnForUnsteadyHeat::Gamma=5.0;

 // Wavenumber
 ExactSolnForUnsteadyHeat::K=2.0;

 // Setup mesh
 //-----------

 // Number of elements in x and y directions
 unsigned nx=5;
 unsigned ny=5;

 // Lengths in x and y directions
 double lx=1.0;
 double ly=1.0;

 // Build mesh
 mesh_pt() = new RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,time_stepper_pt());

 // Choose a control node at which the solution is documented
 //----------------------------------------------------------
 // Total number of elements
 unsigned n_el=mesh_pt()->nelement();

 // Choose an element in the middle
 unsigned control_el=unsigned(n_el/2);

 // Choose its first node as the control node
 Control_node_pt=mesh_pt()->finite_element_pt(control_el)->node_pt(0);

 oomph_info << "Recording trace of the solution at: " 
            << Control_node_pt->x(0) << " "
            << Control_node_pt->x(1) << std::endl;
 

 // Set the boundary conditions for this problem: 
 // ---------------------------------------------
 // All nodes are free by default -- just pin the ones that have 
 // Dirichlet conditions here. 
 unsigned n_bound = mesh_pt()->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   unsigned n_node = mesh_pt()->nboundary_node(b);
   for (unsigned n=0;n<n_node;n++)
    {
     mesh_pt()->boundary_node_pt(b,n)->pin(0); 
    }
  } // end of set boundary conditions


 // Complete the build of all elements so they are fully functional
 //----------------------------------------------------------------

 // Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;

   // Set pointer to continous time
   el_pt->time_pt()=time_pt();
  }

 // Do equation numbering
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//=========start of actions_before_implicit_timestep===============================
/// \short Actions before timestep: update the domain, then reset the 
/// boundary conditions for the current time.
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::actions_before_implicit_timestep()
{
 // Get current time
 double time=time_pt()->time();
 
 //Loop over the boundaries
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
     double u;
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);
     // Get current values of the boundary conditions from the
     // exact solution
     ExactSolnForUnsteadyHeat::get_exact_u(time,x,u);
     nod_pt->set_value(0,u);
    }
  }
} // end of actions_before_implicit_timestep



//======================start_of_set_initial_condition====================
/// \short Set initial condition: Assign previous and current values
/// from exact solution or from restart file.
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::set_initial_condition()
{ 

 // Pointer to restart file
 ifstream* restart_file_pt=0;

 // Restart?
 //---------
 // Restart file specified via command line [all programs have at least
 // a single command line argument: their name. Ignore this here.]
 if ((CommandLineArgs::Argc==1)||(CommandLineArgs::Argc==2))
  {
   oomph_info << "No restart -- setting IC from exact solution" << std::endl;
  }
 else if (CommandLineArgs::Argc==3)
  {
   // Open restart file from stem
   char filename[100];
   sprintf(filename,"%s_on_proc%i.dat",CommandLineArgs::Argv[1],
           this->communicator_pt()->my_rank());
   restart_file_pt= new ifstream(filename,ios_base::in);
   if (restart_file_pt!=0)
    {
     oomph_info << "Have opened " << filename << 
      " for restart. " << std::endl;
    }
   else
    {
     std::ostringstream error_stream;
     error_stream 
      << "ERROR while trying to open " << filename << 
      " for restart." << std::endl;
     
     throw OomphLibError(
      error_stream.str(),
      "UnsteadyHeatProblem<ELEMENT>::set_initial_condition()",
      OOMPH_EXCEPTION_LOCATION);
    }
  }
 // Wrong number of command line arguments?
 else 
  {
   std::ostringstream error_stream;
   error_stream << "Can only specify one or two input file(s).\n" 
                << "You specified the following command line arguments:\n";
   //Fix this
   CommandLineArgs::output();
   
   throw OomphLibError( 
    error_stream.str(),
    "UnsteadyHeatProblem<ELEMENT>::set_initial_condition()",
    OOMPH_EXCEPTION_LOCATION);
  }


 // Read restart data:
 //-------------------
 if (restart_file_pt!=0)
  {
   // Read the data from restart file. This also initialises
   // the previous timesteps and sets up the weights for the timestepper(s)
   // in the problem
   restart(*restart_file_pt);
  }
 // Assign initial condition from exact solution
 //---------------------------------------------
 else
  {
   
   // Choose initial timestep
   double dt0=0.05;

   // Initialise timestep -- also sets the weights for all timesteppers
   // in the problem.
   initialise_dt(dt0);

   // Backup time in global Time object
   double backed_up_time=time_pt()->time();
   
   // Past history fo needs to be established for t=time0-deltat, ...
   // Then provide current values (at t=time0) which will also form
   // the initial guess for the first solve at t=time0+deltat
   
   // Vector of exact solution value
   Vector<double> soln(1);
   Vector<double> x(2);
   
   //Find number of nodes in mesh
   unsigned num_nod = mesh_pt()->nnode();
   
   // Set continuous times at previous timesteps
   int nprev_steps=time_stepper_pt()->nprev_values();
   Vector<double> prev_time(nprev_steps+1);
   for (int itime=nprev_steps;itime>=0;itime--)
    {
     prev_time[itime]=time_pt()->time(unsigned(itime));
    } 
   
   // Loop over current & previous timesteps
   for (int itime=nprev_steps;itime>=0;itime--)
    {
     // Continuous time
     double time=prev_time[itime];
     oomph_info << "setting IC at time =" << time << std::endl;
     
     // Loop over the nodes to set initial guess everywhere
     for (unsigned n=0;n<num_nod;n++)
      {
       // Get nodal coordinates
       x[0]=mesh_pt()->node_pt(n)->x(0);
       x[1]=mesh_pt()->node_pt(n)->x(1);
       
       // Get intial solution
       ExactSolnForUnsteadyHeat::get_exact_u(time,x,soln);
       
       // Assign solution
       mesh_pt()->node_pt(n)->set_value(itime,0,soln[0]);
       
       // Loop over coordinate directions: Previous position = present position
       for (unsigned i=0;i<2;i++)
        {
         mesh_pt()->node_pt(n)->x(itime,i)=x[i];
        }
      } 
    }
   
   // Reset backed up time for global timestepper
   time_pt()->time()=backed_up_time;

  }


} // end of set_initial_condition



//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::
doc_solution(DocInfo& doc_info,ofstream& trace_file)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;


 oomph_info << std::endl;
 oomph_info << "=================================================" << std::endl;
 oomph_info << "Docing solution for t=" << time_pt()->time() << std::endl;
 oomph_info << "=================================================" << std::endl;


 oomph_info << " Timestep: " << doc_info.number() << std::endl;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),this->communicator_pt()->my_rank());
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

 // Write dummy zones
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "0.0 0.0 0.0" << std::endl;
 some_file << "1.0 0.0 0.0" << std::endl;
 some_file << "0.0 1.0 0.0" << std::endl;
 some_file << "1.0 1.0 0.0" << std::endl;
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "0.0 0.0 1.0" << std::endl;
 some_file << "1.0 0.0 1.0" << std::endl;
 some_file << "0.0 1.0 1.0" << std::endl;
 some_file << "1.0 1.0 1.0" << std::endl;
 some_file.close();


 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,time_pt()->time(),
                       ExactSolnForUnsteadyHeat::get_exact_u); 
 // Write dummy zones
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "0.0 0.0 0.0" << std::endl;
 some_file << "1.0 0.0 0.0" << std::endl;
 some_file << "0.0 1.0 0.0" << std::endl;
 some_file << "1.0 1.0 0.0" << std::endl;
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "0.0 0.0 1.0" << std::endl;
 some_file << "1.0 0.0 1.0" << std::endl;
 some_file << "0.0 1.0 1.0" << std::endl;
 some_file << "1.0 1.0 1.0" << std::endl;
 some_file.close();

 // Doc error
 //----------
 double error,norm;
 sprintf(filename,"%s/error%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,
                          ExactSolnForUnsteadyHeat::get_exact_u,
                          time_pt()->time(),
                          error,norm); 
 some_file.close();

 // Collect error and norm from all processors
 Vector<double> local_stuff(2);
 local_stuff[0]=error;
 local_stuff[1]=norm;
 Vector<double> global_stuff(2);
#ifdef OOMPH_HAS_MPI
 MPI_Allreduce(&local_stuff[0],&global_stuff[0],2,MPI_DOUBLE,MPI_SUM,
               communicator_pt()->mpi_comm());
#else
 global_stuff[0]=error;
 global_stuff[1]=norm;
#endif

 // Doc solution and error
 //-----------------------
 oomph_info << "error: " << global_stuff[0] << std::endl; 
 oomph_info << "norm : " << global_stuff[1] << std::endl << std::endl;

 // Get exact solution at control node
 Vector<double> x_ctrl(2);
 x_ctrl[0]=Control_node_pt->x(0);
 x_ctrl[1]=Control_node_pt->x(1);
 double u_exact;
 ExactSolnForUnsteadyHeat::get_exact_u(time_pt()->time(),x_ctrl,u_exact);

 // Write trace file on root only
 if (Communicator_pt->my_rank()==0)
  {
   trace_file << time_pt()->time() << " " 
              << time_pt()->dt() << " " 
              << Control_node_pt->value(0) << " " 
              << u_exact << " " 
              << error   << " " 
              << norm    << std::endl;
  }

 // Write restart file
 sprintf(filename,"%s/restart%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 dump_it(some_file);
 some_file.close();

} // end of doc_solution




//=====start_of_dump_it===================================================
/// Dump the solution to disk to allow for restart
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::dump_it(ofstream& dump_file)
{

 // Call generic dump()
 Problem::dump(dump_file); 

} // end of dump_it



//=================start_of_restart=======================================
/// Read solution from disk for restart
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::restart(ifstream& restart_file)
{

 // Read the generic problem data from restart file
 Problem::read(restart_file);

} // end of restart




//========start_of_global_temporal_error_norm==============================
/// Global error norm for adaptive timestepping: RMS error, based on
/// difference between predicted and actual value at all nodes.
//=========================================================================
template<class ELEMENT>
double UnsteadyHeatProblem<ELEMENT>::global_temporal_error_norm()
{

#ifdef OOMPH_HAS_MPI

 double global_error = 0.0;
   
 //Find out how many nodes there are in the problem
 unsigned n_node = mesh_pt()->nnode();

 // Loop over the nodes and calculate the estimated error in the values
 // for non-haloes
 for(unsigned i=0;i<n_node;i++)
  {
   Node* nod_pt=mesh_pt()->node_pt(i);

   if (!(nod_pt->is_halo()))
    {
     // Get error in solution: Difference between predicted and actual
     // value for nodal value 0
     double error = nod_pt->time_stepper_pt()->
      temporal_error_in_value(nod_pt,0);
     
     //Add the square of the individual error to the global error
     global_error += error*error;
    }
  }

 // Accumulate
 int n_node_local=n_node;
 int n_node_total=0;
 MPI_Allreduce(&n_node_local,&n_node_total,1,MPI_INT,MPI_SUM,
               communicator_pt()->mpi_comm());

 double global_error_total=0.0;
 MPI_Allreduce(&global_error,&global_error_total,1,MPI_DOUBLE,MPI_SUM,
               communicator_pt()->mpi_comm());
 
 // Divide by the number of nodes
 global_error_total /= double(n_node_total);

 // Return square root...
 return sqrt(global_error_total);

#else

 double global_error = 0.0;
   
 //Find out how many nodes there are in the problem
 unsigned n_node = mesh_pt()->nnode();

 //Loop over the nodes and calculate the estimated error in the values
 for(unsigned i=0;i<n_node;i++)
  {
   // Get error in solution: Difference between predicted and actual
   // value for nodal value 0
   double error = mesh_pt()->node_pt(i)->time_stepper_pt()->
    temporal_error_in_value(mesh_pt()->node_pt(i),0);
   
   //Add the square of the individual error to the global error
   global_error += error*error;
  }
    
 // Divide by the number of nodes
 global_error /= double(n_node);

 // Return square root...
 return sqrt(global_error);

#endif

} // end of global_temporal_error_norm

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======start_of_main====================================================
/// \short Driver code for the adaptive solution of an
/// unsteady heat equation with option for restart from disk: 
/// Only a single command line argument is allowed.
/// If specified it is interpreted as the name of the restart file.
//========================================================================
int main(int argc, char* argv[])
{

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Build problem
 UnsteadyHeatProblem<QUnsteadyHeatElement<2,4> >
  problem(&ExactSolnForUnsteadyHeat::get_source);
 
 // Setup labels for output
 DocInfo doc_info;

 // Output directory
 doc_info.set_directory("RESLT");
 
 // Output number
 doc_info.number()=0;
 
 // Open a trace file (on root only)
 ofstream trace_file;
 if (problem.communicator_pt()->my_rank()==0)
  {
   char filename[100];   
   sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
   trace_file.open(filename);
   trace_file << "VARIABLES=\"time\",\"dt\",\"u<SUB>FE</SUB>\","
              << "\"u<SUB>exact</SUB>\",\"norm of error\",\"norm of solution\""
              << std::endl;
  }

 // Choose simulation interval
 double t_max=0.18; // 2.0;

 // Create storage for partitioning
 const unsigned n_element=problem.mesh_pt()->nelement();
 Vector<unsigned> element_partition(n_element);
 
 //No command line arguments -- simply distribute problem and
 // record distribution in file to allow restart
 if (CommandLineArgs::Argc==1)
  {
   // Distribute
   element_partition=problem.distribute();

   // Write partition to disk
   std::ofstream output_file;
   char filename[100];
   sprintf(filename,"%s/partitioning.dat",doc_info.directory().c_str());
   output_file.open(filename);
   for (unsigned e=0;e<n_element;e++)
    {
     output_file << element_partition[e] << std::endl;
    }
  }
 // Read in partitioning from disk using the specified file (during
 // self test)
 else if (CommandLineArgs::Argc==2)
  {
   // Read in partitioning from disk: Name of partitioning file specified
   // as second command line argument
   std::ifstream input_file;
   input_file.open(CommandLineArgs::Argv[1]);
   std::string input_string;
   for (unsigned e=0;e<n_element;e++)
    {
     getline(input_file,input_string,'\n');
     element_partition[e]=atoi(input_string.c_str());
    }
   
   // Now perform the distribution 
   problem.distribute(element_partition);
  }
 // Read in partitioning from disk using the specified file (during
 // restart)
 else if (CommandLineArgs::Argc==3)
  {
   // Read in partitioning from disk: Name of partitioning file specified
   // as second command line argument
   std::ifstream input_file;
   input_file.open(CommandLineArgs::Argv[2]);
   std::string input_string;
   for (unsigned e=0;e<n_element;e++)
    {
     getline(input_file,input_string,'\n');
     element_partition[e]=atoi(input_string.c_str());
    }
   
   // Now perform the distribution 
   problem.distribute(element_partition);
  }

 // Set IC (either analytical or restart)
 problem.set_initial_condition();

 // Initial timestep: Use the one used when setting up the initial
 // condition
 double dt=problem.time_pt()->dt();
 
 //Output initial condition
 problem.doc_solution(doc_info,trace_file);
 
 //Increment counter for solutions 
 doc_info.number()++;

 // Target error for adaptive timestepping
 double epsilon_t=1.0e-4;

 // Timestepping loop: Don't know how many steps we're going to take
 // in advance
 while (problem.time_pt()->time()<t_max)
  {
  
   // Take an adaptive timestep -- the input dt is the suggested timestep.
   // The adaptive timestepper will adjust dt until the required error
   // tolerance is satisfied. The function returns a suggestion
   // for the timestep that should be taken for the next step. This
   // is either the actual timestep taken this time or a larger
   // value if the solution was found to be "too accurate". 
   double dt_next=problem.adaptive_unsteady_newton_solve(dt,epsilon_t);

   // Use dt_next as suggestion for the next timestep
   dt=dt_next;
 
   //Output solution
   problem.doc_solution(doc_info,trace_file);
   
   //Increment counter for solutions 
   doc_info.number()++;

  } // end of timestepping loop
 
 // Close trace file
 if (problem.communicator_pt()->my_rank()==0)
  {
   trace_file.close();
  }

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

}; // end of main
