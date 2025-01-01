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
// Generic oomph-lib header
#include "generic.h"

#include <fenv.h>
using namespace std;

using namespace oomph;

//==========================================================
// Global namespace for spring with contact
//==========================================================
namespace GlobalFct
{
 /// Force
 double Force=0.0;

 /// Spring stiffness
 double Stiffness=1.0;

 /// Max. displacement
 double U_max=5.5;

 /// Number of primary variables. Equal to number of raw "displacement"
 /// variables for contact problem; total number of degrees of freedom
 /// in that case is twice as big because of the unknown contact forces.
 unsigned N_primary=1;

 /// Reset unknowns to some nontrivial garbage
 void reset(Vector<double>& unknowns)
 {
  unsigned n=N_primary*2;
  unknowns.resize(n);
  for (unsigned j=0;j<n;j++)
   {
    unknowns[j]=double(j)+20.5;
   }
 }


 /// Global residual fct
 void get_residuals(const Vector<double>& param, 
                    const Vector<double>& unknowns, 
                    Vector<double>& residuals)
 {
  residuals[0]=Stiffness*(U_max-unknowns[0]*unknowns[0])+
   unknowns[1]*unknowns[1]-Force;

  residuals[1]=unknowns[0]*unknowns[1];
 }


 
 

 /// Plot "landscape" of residuals (only for 2D problems!)
 void plot_it(const std::string filename)
 {
  if (N_primary!=1)
   {
    oomph_info << "Skipping plot_it()\n";
    return;
   }

  ofstream outfile;
  outfile.open(filename.c_str());
  Vector<double> x(2);
  Vector<double> f(2);
  Vector<double> params(0);
  double x_min=0.0;
  double x_max=10.0;
  double y_min=-5.0;
  double y_max=5.0;
 
  unsigned nplot=100;
  outfile << "ZONE I=" << nplot << ", J=" << nplot << std::endl;
  for (unsigned i=0;i<nplot;i++)
   {
    x[0]=x_min+(x_max-x_min)*double(i)/double(nplot-1);
    for (unsigned j=0;j<nplot;j++)
     {
      x[1]=y_min+(y_max-y_min)*double(j)/double(nplot-1);
      get_residuals(params,x,f);
      outfile << x[0] << " " 
              << x[1] << " " 
              << f[0] << " " 
              << f[1] << "\n" ;
     }
   }
  outfile.close();

 }
}

/// ///////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////


//==start_of_main======================================================
/// 
//=====================================================================
int main(int argc, char **argv)
{
 //feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // CommandLineArgs::specify_command_line_flag("--bilinear");
  
 // CommandLineArgs::specify_command_line_flag("--completely_smooth");

 // CommandLineArgs::specify_command_line_flag("--single_kink");

 // CommandLineArgs::specify_command_line_flag("--old_version");

 // CommandLineArgs::specify_command_line_flag("--kuhn_tucker");

 // CommandLineArgs::specify_command_line_flag("--n_primary",
 //                                            &GlobalFct::N_primary);

 // Doc available command line flags
 CommandLineArgs::doc_available_flags();

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Number of steps for load incrementation
 unsigned nstep=10;
 unsigned nstep_total=2*nstep;

 BlackBoxFDNewtonSolver::Max_iter=30;
 
 ofstream outfile2;
 outfile2.open("trace.dat");

 BlackBoxFDNewtonSolver::Doc_Progress=true;
 Vector<double> param;
 Vector<double> unknowns;
 GlobalFct::reset(unknowns);
 
 GlobalFct::Force=0.0;
 for (unsigned j=0;j<nstep_total;j++)
  {
   
   oomph_info << "\n\nSolving for force: " << GlobalFct::Force
              << std::endl << std::endl;
   
   char name[100];   
   sprintf(name,"landscape%i.dat",j);
   std::string filename(name);
   GlobalFct::plot_it(filename);
   try
    {
     //BlackBoxFDNewtonSolver::Use_step_length_control=true;

     BlackBoxFDNewtonSolver::black_box_fd_newton_solve(
      GlobalFct::get_residuals,param,unknowns);
     
     oomph_info << "Number of Newton iterations to convergence: "
                << BlackBoxFDNewtonSolver::N_iter_taken++ << std::endl;
     
     char name[100];   
     sprintf(name,"soln%i.dat",j);
     ofstream outfile;
     outfile.open(name);
     outfile <<  unknowns[0] << " " 
             <<  unknowns[1] << "\n";
     outfile.close();
     
     oomph_info <<  "Displ: " << GlobalFct::U_max-unknowns[0]*unknowns[0]<< " " 
                <<  "Contact force: " << unknowns[1]*unknowns[1]<< " " 
                << std::endl;


     unsigned n=unknowns.size();
     outfile2 << GlobalFct::Force << " ";
     for (unsigned k=0;k<n;k++)
      {
       outfile2 << unknowns[k] << " ";
      }
     outfile2 << std::endl;
    }
   catch(OomphLibError) {}
   
   // Sweep up and then down
   if (j<nstep)
    {
     GlobalFct::Force+=2.1;
    }
   else
    {
     GlobalFct::Force-=2.1;
    }
  }
 
 outfile2.close();


 // //----------------------------------------

 // oomph_info << "\n\n\n\nTRY POSSIBLY SPURIOUS SOLUTION\n"
 //            <<         "==============================\n\n";

 // // Try spurious solution
 // GlobalFct::set_spurious_solution(unknowns);

 // BlackBoxFDNewtonSolver::black_box_fd_newton_solve(
 //  GlobalFct::get_residuals,param,unknowns);
 
 // oomph_info << "Number of Newton iterations to convergence: "
 //            << BlackBoxFDNewtonSolver::N_iter_taken++ << std::endl;
 
 // oomph_info << "\nPossibly spurious solution: "
 //            << unknowns[0] << " " <<  unknowns[1] << "\n\n";

 
 
} // end_of_main











