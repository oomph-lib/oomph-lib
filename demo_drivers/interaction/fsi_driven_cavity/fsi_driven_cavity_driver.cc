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
#include "fsi_driven_cavity_problem.h"


using namespace std;
using namespace oomph;


//============start_of_main====================================================
/// Driver code for a collapsible channel problem with FSI.
/// Presence of command line arguments indicates validation run with 
/// coarse resolution and small number of timesteps.
//=============================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // Number of elements in the domain (W. Wall standard resolution)
 unsigned nx=16; //32;
 unsigned ny=16; //32;

 // Length of the domain
 double lx=1.0;
 double ly=1.0;

 // Fractional height of the gap next to the moving lid
 double gap_fraction=1.0/16.0;
   
 // Period of lid oscillation
 double period=5.0;

 // Build the problem
 FSIDrivenCavityProblem
  <AlgebraicElement<QTaylorHoodElement<2> > > 
  problem(nx, ny, lx, ly, gap_fraction, period);


 // Timestep. 
 double dt=period/40.0; 
 
 // Maximum time for simulation
 double t_max=50.0; 

 // Initialise timestep 
 problem.initialise_dt(dt);

 // Apply initial condition
 problem.set_initial_condition();

 //Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 
 // Open a trace file 
 ofstream trace_file;
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);
 
 // Output the initial solution
 problem.doc_solution(doc_info, trace_file);

 // Increment step number
 doc_info.number()++;

 // Find number of timesteps (reduced for validation)
 unsigned nstep = unsigned(t_max/dt);
 if (CommandLineArgs::Argc>1)
  {
   nstep=3;
  }
 
 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   // Solve the problem
   problem.unsteady_newton_solve(dt);
   
   // Outpt the solution
   problem.doc_solution(doc_info, trace_file);
   
   // Step number
   doc_info.number()++;
  }


 // Close trace file.
 trace_file.close();

}//end of main
  
