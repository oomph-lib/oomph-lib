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
//Driver code for a large-amplitude ring oscillation

//Generic stuff
#include "generic.h"

// The beam equations
#include "beam.h"

// The mesh
#include "meshes/one_d_lagrangian_mesh.h"

using namespace std;

using namespace oomph;


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//====start_of_namespace============================
/// Namespace for global parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Perturbation pressure
 double Pcos;

 /// Duration of transient load
 double T_kick;

 /// Load function: Perturbation pressure to force non-axisymmetric deformation
 void press_load(const Vector<double>& xi,
                 const Vector<double> &x,
                 const Vector<double>& N,
                 Vector<double>& load)
 {
  for(unsigned i=0;i<2;i++) 
   {
    load[i] = -Pcos*cos(2.0*xi[0])*N[i];
   }
 } //end of load


 /// Scaling factor for wall thickness (to be used in an exercise)
 double Alpha=1.0;

 /// Wall thickness -- 1/20 for default value of scaling factor
 double H=Alpha*1.0/20.0;

 /// Square of timescale ratio (i.e. non-dimensional density)  
 /// -- 1.0 for default value of scaling factor
 double Lambda_sq=pow(Alpha,2);

} // end of namespace





/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////



//===start_of_problem_class=============================================
/// Ring problem
//======================================================================
template<class ELEMENT, class TIMESTEPPER>
class ElasticRingProblem : public Problem
{

public:

 /// Constructor: Number of elements
 ElasticRingProblem(const unsigned& n_element);

 /// Access function for the specific mesh
 OneDLagrangianMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<OneDLagrangianMesh<ELEMENT>*>(Problem::mesh_pt());
  }

 /// Update function is empty 
 void actions_after_newton_solve() {}

 /// Update function is empty 
 void actions_before_newton_solve() {}
 
 /// Setup initial conditions
 void set_initial_conditions();

 /// Doc solution
 void doc_solution(DocInfo& doc_info);

 /// Do unsteady run
 void unsteady_run();

 /// Dump problem-specific parameter values, then dump
 /// generic problem data.
 void dump_it(ofstream& dump_file);

 /// Read problem-specific parameter values, then recover
 /// generic problem data.
 void restart(ifstream& restart_file);


private:

 /// Trace file for recording control data
 ofstream Trace_file;

 /// Flag for validation run: Default: 0 = no validation run
 unsigned Validation_run_flag;

 /// Restart flag specified via command line?
 bool Restart_flag;

}; // end of problem class





//===start_of_constructor===============================================
/// Constructor for elastic ring problem
//======================================================================
template<class ELEMENT,class TIMESTEPPER>
ElasticRingProblem<ELEMENT,TIMESTEPPER>::ElasticRingProblem
(const unsigned& n_element) :
 Validation_run_flag(0), //default: false
 Restart_flag(false)
{

 // Create the timestepper and add it to the Problem's collection of
 // timesteppers -- this creates the Problem's Time object.
 add_time_stepper_pt(new TIMESTEPPER());

 // Undeformed beam is an ellipse with unit axes
 GeomObject* undef_geom_pt=new Ellipse(1.0,1.0); 

 //Length of domain
 double length = MathematicalConstants::Pi/2.0;
 
 //Now create the (Lagrangian!) mesh
 Problem::mesh_pt() = new OneDLagrangianMesh<ELEMENT>(
  n_element,length,undef_geom_pt,Problem::time_stepper_pt()); 

 // Boundary condition: 

 // Bottom: 
 unsigned ibound=0;
 // No vertical displacement
 mesh_pt()->boundary_node_pt(ibound,0)->pin_position(1); 
 // Zero slope: Pin type 1 (slope) dof for displacement direction 0 
 mesh_pt()->boundary_node_pt(ibound,0)->pin_position(1,0);

 // Top: 
 ibound=1;
 // No horizontal displacement
 mesh_pt()->boundary_node_pt(ibound,0)->pin_position(0); 
 // Zero slope: Pin type 1 (slope) dof for displacement direction 1
 mesh_pt()->boundary_node_pt(ibound,0)->pin_position(1,1); 

 // Complete build of all elements so they are fully functional
 // -----------------------------------------------------------

 //Loop over the elements to set physical parameters etc.
 for(unsigned i=0;i<n_element;i++)
  {
   // Cast to proper element type
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   // Pass pointer to square of timescale ratio (non-dimensional density)
   elem_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;

   // Pass pointer to non-dimensional wall thickness
   elem_pt->h_pt() = &Global_Physical_Variables::H;

   // Function that specifies load vector
   elem_pt->load_vector_fct_pt() = &Global_Physical_Variables::press_load;
   
   // Assign the undeformed surface
   elem_pt->undeformed_beam_pt() = undef_geom_pt;
  }
 
 // Do equation numbering
 cout << "# of dofs " << assign_eqn_numbers() << std::endl;
 
} // end of constructor




//====start_of_doc_solution===============================================
/// Document solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void ElasticRingProblem<ELEMENT, TIMESTEPPER>::doc_solution(
 DocInfo& doc_info)
{ 
 
 cout << "Doc-ing step " <<  doc_info.number()
      << " for time " << time_stepper_pt()->time_pt()->time() << std::endl;
 
 
 // Loop over all elements to get global kinetic and potential energy
 unsigned n_elem=mesh_pt()->nelement();
 double global_kin=0;
 double global_pot=0;
 double pot,kin;
 for (unsigned ielem=0;ielem<n_elem;ielem++)
  {
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(ielem))->get_energy(pot,kin);
   global_kin+=kin;
   global_pot+=pot;
  }
  

 // Get pointer to last element to document displacement
 FiniteElement* trace_elem_pt=mesh_pt()->finite_element_pt(n_elem-1);
 
 // Vector of local coordinates at control point
 Vector<double> s_trace(1);
 s_trace[0]=1.0;
 
 // Write trace file: Time, control position, energies
 Trace_file << time_pt()->time()  << " " 
            << trace_elem_pt->interpolated_x(s_trace,1) 
            << " " << global_pot  << " " << global_kin
            << " " << global_pot + global_kin
            << std::endl; // end of output to trace file
  
 ofstream some_file;
 char filename[100];
  
 // Number of plot points
 unsigned npts=5;
  
 // Output solution 
 sprintf(filename,"%s/ring%i.dat",doc_info.directory().c_str(),
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
 
 
 // Loop over all elements do dump out previous solutions
 unsigned nsteps=time_stepper_pt()->nprev_values();
 for (unsigned t=0;t<=nsteps;t++)
  {     
   sprintf(filename,"%s/ring%i-%i.dat",doc_info.directory().c_str(),
           doc_info.number(),t);
   some_file.open(filename);
   unsigned n_elem=mesh_pt()->nelement();
   for (unsigned ielem=0;ielem<n_elem;ielem++)
    {
     dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(ielem))->
      output(t,some_file,npts);
    }
   some_file.close();
  } // end of output of previous solutions
  
 
 // Write restart file
 sprintf(filename,"%s/ring_restart%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 dump_it(some_file);
 some_file.close();

       
} // end of doc solution






//===start_of_dump_it====================================================
/// Dump problem-specific parameter values, then dump
/// generic problem data.
//=======================================================================
template<class ELEMENT,class TIMESTEPPER>
void ElasticRingProblem<ELEMENT,TIMESTEPPER>::dump_it(ofstream& dump_file)

{
  
 // Write Pcos
 dump_file << Global_Physical_Variables::Pcos << " # Pcos" << std::endl;

 // Write validation run flag
 dump_file << Validation_run_flag 
           << " # Validation run flag" << std::endl;

 // Call generic dump()
 Problem::dump(dump_file);

} // end of dump it


//==start_of_restart=====================================================
/// Read problem-specific parameter values, then recover
/// generic problem data.
//=======================================================================
template<class ELEMENT,class TIMESTEPPER>
void ElasticRingProblem<ELEMENT,TIMESTEPPER>::restart(ifstream& restart_file)
{
  
 string input_string;

 // Read line up to termination sign
 getline(restart_file,input_string,'#');
 // Ignore rest of line
 restart_file.ignore(80,'\n');
 // Read in pcos
 Global_Physical_Variables::Pcos=atof(input_string.c_str());
 
 // Read line up to termination sign
 getline(restart_file,input_string,'#');
 // Ignore rest of line
 restart_file.ignore(80,'\n');
 // Read in Long run flag
 Validation_run_flag=
  unsigned(atof(input_string.c_str()));

 // Call generic read()
 Problem::read(restart_file);

} // end of restart




//=====start_of_set_ic=====================================================
/// Setup initial conditions -- either restart from solution 
/// specified via command line or impulsive start. 
//=========================================================================
template<class ELEMENT,class TIMESTEPPER>
void ElasticRingProblem<ELEMENT,TIMESTEPPER>::set_initial_conditions()
{

 // No restart file --> impulsive start from initial configuration
 // assigned in the Lagrangian mesh.
 if (!Restart_flag)
  {
   // Set initial timestep
   double dt=1.0; 

   // Assign impulsive start
   assign_initial_values_impulsive(dt);
  }
 // Restart file specified via command line
 else
  {
   // Try to open restart file
   ifstream* restart_file_pt= 
    new ifstream(CommandLineArgs::Argv[2],ios_base::in);
   if (restart_file_pt!=0)
    {
     cout << "Have opened " << CommandLineArgs::Argv[2] << 
      " for restart. " << std::endl;
    }
   else
    {
     std::ostringstream error_stream;
     error_stream << "ERROR while trying to open " 
                  << CommandLineArgs::Argv[2] 
                  << " for restart." << std::endl;

     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }

   // Read restart data:
   restart(*restart_file_pt);

  }

} // end of set ic



//===start_of_unsteady_run=================================================
/// Solver loop to perform unsteady run.
//=========================================================================
template<class ELEMENT,class TIMESTEPPER>
void ElasticRingProblem<ELEMENT,TIMESTEPPER>::unsteady_run()
{
 
 // Convert command line arguments (if any) into flags:
 //----------------------------------------------------

 if (CommandLineArgs::Argc<2)
  {
   cout << "No command line arguments -- using defaults." 
        << std::endl;
  }
 else if (CommandLineArgs::Argc==2)
  {
   // Flag for validation run
   Validation_run_flag=atoi(CommandLineArgs::Argv[1]);
  }
 else if (CommandLineArgs::Argc==3)
  {
   // Flag for validation run
   Validation_run_flag=atoi(CommandLineArgs::Argv[1]);

   // Second argument is restart file. If it's there we're performing
   // a restart
   Restart_flag=true;
  }
 else
  {
   std::string error_message =
    "Wrong number of command line arguments. Specify two or fewer.\n";
   error_message += "Arg1: Validation_run_flag [0/1] for [false/true]\n";
   error_message += "Arg2: Name of restart_file (optional)\n";
   error_message += "No arguments specified --> default run\n";
   
   throw OomphLibError(error_message,
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }


 // Label for output
 DocInfo doc_info;

 // Output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;

 // Set up trace file
 char filename[100];
 sprintf(filename,"%s/trace_ring.dat",doc_info.directory().c_str());
 Trace_file.open(filename);
 
 Trace_file <<  "VARIABLES=\"time\",\"R<sub>ctrl</sub>\",\"E<sub>pot</sub>\"";
 Trace_file << ",\"E<sub>kin</sub>\",\"E<sub>kin</sub>+E<sub>pot</sub>\"" 
            << std::endl;

 // Perturbation pressure -- incl. scaling factor (Alpha=1.0 by default)
 Global_Physical_Variables::Pcos=1.0e-4
  *pow(Global_Physical_Variables::Alpha,3); 

 // Duration of transient load
 Global_Physical_Variables::T_kick=15.0;

 // Number of steps
 unsigned nstep=600; 
 if (Validation_run_flag==1) {nstep=10;}

 // Setup initial condition (either restart or impulsive start)
 set_initial_conditions();

 // Extract initial timestep as set up in set_initial_conditions()
 double dt=time_pt()->dt();

 // Output initial data
 doc_solution(doc_info);
 
 // If the run is restarted, dt() contains the size of the timestep
 // that was taken to reach the dumped solution. In a non-restarted run
 // the next step would have been taken with a slightly smaller 
 // timestep (see below). For a restarted run we therefore adjust 
 // the next timestep accordingly.
 if (Restart_flag&&(time_pt()->time()>10.0)&&(time_pt()->time()<100.0))
  {
   dt=0.995*dt;
  }

 // Time integration loop
 for(unsigned i=1;i<=nstep;i++)
  {
   // Switch off perturbation pressure 
   if (time_pt()->time()>Global_Physical_Variables::T_kick) 
     {
      /// Perturbation pressure
      Global_Physical_Variables::Pcos=0.0; 
     }
   
   // Solve
   unsteady_newton_solve(dt);
   
   // Doc solution
   doc_info.number()++;
   doc_solution(doc_info);
   
   // Reduce timestep for the next solve
   if ((time_pt()->time()>10.0)&&(time_pt()->time()<100.0))
    {
     dt=0.995*dt;
    }

  }
 
} // end of unsteady run




//===start_of_main=====================================================
/// Driver for oscillating ring problem 
//=====================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Number of elements
 unsigned nelem = 13;

 //Set up the problem
 ElasticRingProblem<HermiteBeamElement,Newmark<3> > problem(nelem);

 // Do unsteady run
 problem.unsteady_run();

} // end of main








