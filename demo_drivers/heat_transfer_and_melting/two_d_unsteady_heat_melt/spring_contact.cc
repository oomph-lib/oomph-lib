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
// Generic oomph-lib header
#include "generic.h"

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
  // Loop over primary unknowns
  for (unsigned j=0;j<N_primary;j++)
   {    
    
    // Mash 'em up
    double alpha=1.0/double(5*N_primary);
    double u=unknowns[j];
    double u_combined=u;
    for (unsigned k=0;k<N_primary;k++)
     {    
      if (k!=j)
       {
        u_combined+=unknowns[k]*alpha;
       }
     }
    
    // Contact force
    double fc=unknowns[N_primary+j];
    
    // Eqn numbers
    unsigned displ_eqn=j;
    unsigned contact_eqn=N_primary+j;
    
    // Equilibrium equation
    //=====================
    double N=double(N_primary);
    residuals[displ_eqn]=Stiffness*u_combined-
     Force*(2.0*j+2.0+alpha*N+alpha*N*N-2.0*alpha*j-2.0*alpha)/2.0
     +fc;
    
    // Contact equation
    //=================

    // Old version
    //------------
    if (CommandLineArgs::command_line_flag_has_been_set
        ("--old_version"))
     {
      
      // Equation for "melt" rate
      if (u<U_max)
       {
        // "Lower left" quadrant of (u-U_max,fc) space:
        // Linear variation with m forces things back to m=0
        residuals[contact_eqn]=fc;
       }
      else
       {
        // Lower right quadrant of (u-U_max,fc) space
        if (fc<0.0)
         {
          residuals[contact_eqn]=fc;
         }
        // Upper right quadrant of (u-U_max,fc) space
        else
         {
          // Linear variation forces things back to U_max and maintains
          // continuity with upper left and lower right quadrants
          residuals[contact_eqn]=u-U_max; 
         }
       }
     }
    // Completely smooth surface for melt rate residual
    //-------------------------------------------------
    else if (CommandLineArgs::command_line_flag_has_been_set
             ("--completely_smooth"))
     {
      double phi=atan2(fc,(U_max-u));
      if (phi<0) phi+=2.0*MathematicalConstants::Pi;
      
      // Quadratic fit to periodic polynomial that's only
      // above zero between phi=90 and 180 degrees.
      double t0 = 768.0/77.0/0.3141592653589793E1*phi-2176.0/77.0/(
       0.3141592653589793E1*0.3141592653589793E1)*phi*phi+128.0/7.0/(
        0.3141592653589793E1*0.3141592653589793E1*0.3141592653589793E1)*
       phi*phi*phi
       -256.0/77.0/(0.3141592653589793E1*0.3141592653589793E1*
                    0.3141592653589793E1*
                    0.3141592653589793E1)*phi*phi*phi*phi;
      
      residuals[contact_eqn]=(fc*fc+(U_max-u)*(U_max-u))*t0;
     }
    // Single kink
    //------------
    else if (CommandLineArgs::command_line_flag_has_been_set
             ("--single_kink"))
     {
      // Piecewise linear variation with a single kink
      if ((U_max-u)>fc)
       {
        residuals[contact_eqn]=fc;
       }
      else
       {
        residuals[contact_eqn]=(U_max-u);
       }
     }

    // Kuhn Tucker
    //------------
    else if (CommandLineArgs::command_line_flag_has_been_set
             ("--kuhn_tucker"))
     {
      residuals[contact_eqn]=fc*(U_max-u);
     }
    else
     {
      oomph_info << "Please specify which method you want to use for\n"
                 << "the enforcement of the contact condition: \n";
      CommandLineArgs::doc_available_flags();
      abort();
     }
   }
 }


 
 /// Set unknowns to spurious solution
 void set_spurious_solution(Vector<double>& unknowns)
 {
  Force=100.0;
  N_primary=1;
  
  Vector<double> param; 
  unknowns.resize(2); 
  Vector<double> residuals(2);
  
  // Displacement
  unknowns[0]=Force/Stiffness+0.5*BlackBoxFDNewtonSolver::Tol;
  
  // Contact force 
  unknowns[1]=-0.5*BlackBoxFDNewtonSolver::Tol;
  
  // Get residuals
  get_residuals(param,unknowns, residuals);
  oomph_info << "Residual for spurious solution: " 
             << residuals[0] << " " 
             << residuals[1] << "\n"; 
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
 
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 CommandLineArgs::specify_command_line_flag("--bilinear");
  
 CommandLineArgs::specify_command_line_flag("--completely_smooth");

 CommandLineArgs::specify_command_line_flag("--single_kink");

 CommandLineArgs::specify_command_line_flag("--old_version");

 CommandLineArgs::specify_command_line_flag("--kuhn_tucker");

 CommandLineArgs::specify_command_line_flag("--n_primary",
                                            &GlobalFct::N_primary);

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
     GlobalFct::Force+=1.0;
    }
   else
    {
     GlobalFct::Force-=1.0;
    }
  }
 
 outfile2.close();


 //----------------------------------------

 oomph_info << "\n\n\n\nTRY POSSIBLY SPURIOUS SOLUTION\n"
            <<         "==============================\n\n";

 // Try spurious solution
 GlobalFct::set_spurious_solution(unknowns);

 BlackBoxFDNewtonSolver::black_box_fd_newton_solve(
  GlobalFct::get_residuals,param,unknowns);
 
 oomph_info << "Number of Newton iterations to convergence: "
            << BlackBoxFDNewtonSolver::N_iter_taken++ << std::endl;
 
 oomph_info << "\nPossibly spurious solution: "
            << unknowns[0] << " " <<  unknowns[1] << "\n\n";

 
 
} // end_of_main











