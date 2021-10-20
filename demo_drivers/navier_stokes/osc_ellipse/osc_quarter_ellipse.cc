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
//Driver for 2D Navier-Stokes problem in moving domain

//Generic routines
#include "generic.h"

// The Navier Stokes equations
#include "navier_stokes.h"

// Mesh
#include "meshes/quarter_circle_sector_mesh.h"

using namespace std;

using namespace oomph;

//============start_of_MyEllipse===========================================
///  Oscillating ellipse
/// \f[ x = (A + \widehat{A} \cos(2\pi t/T)) \cos(\xi)  \f]
/// \f[ y = \frac{\sin(\xi)}{A + \widehat{A} \cos(2\pi t/T)}   \f]
/// Note that cross-sectional area is conserved.
//=========================================================================
class MyEllipse : public GeomObject
{

public:

 ///  Constructor:  Pass initial x-half axis, amplitude of x-variation, 
 /// period of oscillation and pointer to time object.
 MyEllipse(const double& a, const double& a_hat,
           const double& period, Time* time_pt) : 
  GeomObject(1,2), A(a), A_hat(a_hat), T(period), Time_pt(time_pt) {}

 /// Destructor: Empty
 virtual ~MyEllipse() {}

 ///  Current position vector to material point at 
 /// Lagrangian coordinate xi 
 void position(const Vector<double>& xi, Vector<double>& r) const
  {
   // Get current time:
   double time=Time_pt->time();

   // Position vector
   double axis=A+A_hat*cos(2.0*MathematicalConstants::Pi*time/T);
   r[0] = axis*cos(xi[0]);
   r[1] = (1.0/axis)*sin(xi[0]);
  } 

 ///  Parametrised position on object: r(xi). Evaluated at
 /// previous time level. t=0: current time; t>0: previous
 /// time level.
 void position(const unsigned& t, const Vector<double>& xi,
               Vector<double>& r) const
  {
   // Get current time:
   double time=Time_pt->time(t);
   
   // Position vector
   double axis=A+A_hat*cos(2.0*MathematicalConstants::Pi*time/T);
   r[0] = axis*cos(xi[0]);
   r[1] = (1.0/axis)*sin(xi[0]);
  } 

private:

 /// x-half axis
 double A;

 /// Amplitude of variation in x-half axis
 double A_hat;

 /// Period of oscillation
 double T;

 /// Pointer to time object
 Time* Time_pt;

}; // end of MyEllipse



/////////////////////////////////////////////////////////////////////// 
///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//===start_of_namespace=================================================
/// Namepspace for global parameters
//======================================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=100.0;

 /// Womersley = Reynolds times Strouhal
 double ReSt=100.0;

 /// x-Half axis length
 double A=1.0;

 /// x-Half axis amplitude
 double A_hat=0.1;

 /// Period of oscillations
 double T=1.0;

 /// Exact solution of the problem as a vector containing u,v,p
 void get_exact_u(const double& t, const Vector<double>& x, Vector<double>& u)
 {
  using namespace MathematicalConstants;

  // Strouhal number
  double St = ReSt/Re;

  // Half axis
  double a=A+A_hat*cos(2.0*Pi*t/T);
  double adot=-2.0*A_hat*Pi*sin(2.0*Pi*t/T)/T; 
  u.resize(3);

  // Velocity solution
  u[0]=adot*x[0]/a;
  u[1]=-adot*x[1]/a;

  // Pressure solution
  u[2]=(2.0*A_hat*Pi*Pi*Re*(x[0]*x[0]*St*cos(2.0*Pi*t/T)*A + 
                            x[0]*x[0]*St*A_hat - x[0]*x[0]*A_hat +
                            x[0]*x[0]*A_hat*cos(2.0*Pi*t/T)*cos(2.0*Pi*t/T) -
                            x[1]*x[1]*St*cos(2.0*Pi*t/T)*A - 
                            x[1]*x[1]*St*A_hat - x[1]*x[1]*A_hat +
                            x[1]*x[1]*A_hat*cos(2.0*Pi*t/T)*cos(2.0*Pi*t/T) ))
   /(T*T*(A*A + 2.0*A*A_hat*cos(2.0*Pi*t/T) + 
          A_hat*A_hat*cos(2.0*Pi*t/T)*cos(2.0*Pi*t/T) ));
 }


} // end of namespace




/////////////////////////////////////////////////////////////////////// 
///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=====start_of_problem_class=========================================
/// Navier-Stokes problem in an oscillating ellipse domain.
//====================================================================
template<class ELEMENT, class TIMESTEPPER>
class OscEllipseProblem : public Problem
{

public:

 /// Constructor
 OscEllipseProblem();

 /// Destructor (empty)
 ~OscEllipseProblem() {}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 ///  Update problem specs before solve (empty)
 void actions_before_newton_solve() {} 
 
 /// Actions before adapt (empty)
 void actions_before_adapt(){}

 /// Actions after adaptation, pin relevant pressures
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());

   // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now set the first pressure dof in the first element to 0.0
   fix_pressure(0,0,0.0);

  } // end of actions_after_adapt


 ///  Update the problem specs before next timestep
 void actions_before_implicit_timestep()
  {
   // Update the domain shape
   mesh_pt()->node_update();

   // Ring boundary: No slip; this implies that the velocity needs
   // to be updated in response to wall motion
   unsigned ibound=1;
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Which node are we dealing with?
     Node* node_pt=mesh_pt()->boundary_node_pt(ibound,inod);
     
     // Apply no slip
     FSI_functions::apply_no_slip_on_moving_wall(node_pt);
    }
  } 

 /// Update the problem specs after timestep (empty)
 void actions_after_implicit_timestep(){}
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Timestepping loop
 void unsteady_run(DocInfo& doc_info);

 ///  Set initial condition
 void set_initial_condition();

private:

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to proper element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

 /// Pointer to GeomObject that specifies the domain bondary
 GeomObject* Wall_pt;

}; // end of problem_class



//========start_of_constructor============================================
/// Constructor for Navier-Stokes problem on an oscillating ellipse domain.
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
OscEllipseProblem<ELEMENT,TIMESTEPPER>::OscEllipseProblem()
{ 

 //Create the timestepper and add it to the problem
 add_time_stepper_pt(new TIMESTEPPER);

 // Setup mesh
 //-----------

 // Build geometric object that forms the curvilinear domain boundary:
 // an oscillating ellipse

 // Half axes
 double a=Global_Physical_Variables::A;

 // Variations of half axes
 double a_hat=Global_Physical_Variables::A_hat;

 // Period of the oscillation
 double period=Global_Physical_Variables::T;

 // Create GeomObject that specifies the domain bondary
 Wall_pt=new MyEllipse(a,a_hat,period,Problem::time_pt()); 


 // Start and end coordinates of curvilinear domain boundary on ellipse
 double xi_lo=0.0;
 double xi_hi=MathematicalConstants::Pi/2.0;

 // Now create the mesh. Separating line between the two 
 // elements next to the curvilinear boundary is located half-way
 // along the boundary.
 double fract_mid=0.5;
 Problem::mesh_pt() = new RefineableQuarterCircleSectorMesh<ELEMENT>(
  Wall_pt,xi_lo,fract_mid,xi_hi,time_stepper_pt());

 // Set error estimator NOT NEEDED IN CURRENT PROBLEM SINCE 
 // WE'RE ONLY REFINING THE MESH UNIFORMLY
 //Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 //mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;


 // Fluid boundary conditions
 //--------------------------
 // Ring boundary: No slip; this also implies that the velocity needs
 // to be updated in response to wall motion
 unsigned ibound=1;
 {
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
     
    // Pin both velocities
    for (unsigned i=0;i<2;i++)
     {
      mesh_pt()->boundary_node_pt(ibound,inod)->pin(i);
     }
   }
 } // end boundary 1
 
 // Bottom boundary: 
 ibound=0;
 {
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Pin vertical velocity
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
    }
   }
 } // end boundary 0

 // Left boundary:
 ibound=2;
 {
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Pin horizontal velocity
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
    }
   }
 } // end boundary 2

 
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

   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
  }

 // Pin redundant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
 // Now set the first pressure dof in the first element to 0.0
 fix_pressure(0,0,0.0);

 // Do equation numbering
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor


//======================start_of_set_initial_condition====================
///  Set initial condition: Assign previous and current values
/// from exact solution.
//========================================================================
template<class ELEMENT,class TIMESTEPPER>
void OscEllipseProblem<ELEMENT,TIMESTEPPER>::set_initial_condition()
{ 
 // Backup time in global timestepper
 double backed_up_time=time_pt()->time();
   
 // Past history for velocities must be established for t=time0-deltat, ...
 // Then provide current values (at t=time0) which will also form
 // the initial guess for first solve at t=time0+deltat
 
 // Vector of exact solution value
 Vector<double> soln(3);
 Vector<double> x(2);
 
 //Find number of nodes in mesh
 unsigned num_nod = mesh_pt()->nnode();
 
 // Get continuous times at previous timesteps
 int nprev_steps=time_stepper_pt()->nprev_values();
 Vector<double> prev_time(nprev_steps+1);
 for (int itime=nprev_steps;itime>=0;itime--)
  {
   prev_time[itime]=time_pt()->time(unsigned(itime));
  }
 
 // Loop over current & previous timesteps (in outer loop because
 // the mesh also moves!)
 for (int itime=nprev_steps;itime>=0;itime--)
  {
   double time=prev_time[itime];
   
   // Set global time (because this is how the geometric object refers 
   // to continous time )
   time_pt()->time()=time;
   
   cout << "setting IC at time =" << time << std::endl;
   
   // Update the mesh for this value of the continuous time
   // (The wall object reads the continous time from the 
   // global time object)
   mesh_pt()->node_update(); 
   
   // Loop over the nodes to set initial guess everywhere
   for (unsigned jnod=0;jnod<num_nod;jnod++)
    {
     // Get nodal coordinates
     x[0]=mesh_pt()->node_pt(jnod)->x(0);
     x[1]=mesh_pt()->node_pt(jnod)->x(1);
     
     // Get exact solution (unsteady stagnation point flow)
     Global_Physical_Variables::get_exact_u(time,x,soln);
     
     // Assign solution
     mesh_pt()->node_pt(jnod)->set_value(itime,0,soln[0]);
     mesh_pt()->node_pt(jnod)->set_value(itime,1,soln[1]);
       
     // Loop over coordinate directions
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->node_pt(jnod)->x(itime,i)=x[i];
      }
    } 
  } // end of loop over previous timesteps
 
 // Reset backed up time for global timestepper
 time_pt()->time()=backed_up_time;
 
} // end of set_initial_condition



//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void OscEllipseProblem<ELEMENT,TIMESTEPPER>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;

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

 // Write dummy zones that force tecplot to keep the axis limits constant
 // while the domain is moving.
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "0.0 0.0 -0.65 -0.65 -200.0" << std::endl;
 some_file << "1.15 0.0 -0.65 -0.65 -200.0" << std::endl;
 some_file << "0.0 1.15 -0.65 -0.65 -200.0" << std::endl;
 some_file << "1.15 1.15 -0.65 -0.65 -200.0" << std::endl;
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "0.0 0.0 0.65 0.65 300.0" << std::endl;
 some_file << "1.15 0.0 0.65 0.65 300.0" << std::endl;
 some_file << "0.0 1.15 0.65 0.65 300.0" << std::endl;
 some_file << "1.15 1.15 0.65 0.65 300.0" << std::endl;

 some_file.close();

 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,time_pt()->time(),
                       Global_Physical_Variables::get_exact_u); 
 some_file.close();

 // Doc error
 //----------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,
                          Global_Physical_Variables::get_exact_u,
                          time_pt()->time(),
                          error,norm); 
 some_file.close();


 // Doc solution and error
 //-----------------------
 cout << "error: " << error << std::endl; 
 cout << "norm : " << norm << std::endl << std::endl;


 // Plot wall posn
 //---------------
 sprintf(filename,"%s/Wall%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 
 unsigned nplot=100;
 for (unsigned iplot=0;iplot<nplot;iplot++)
  {
   Vector<double> xi_wall(1), r_wall(2);
   xi_wall[0]=0.5*MathematicalConstants::Pi*double(iplot)/double(nplot-1);
   Wall_pt->position(xi_wall,r_wall);
   some_file << r_wall[0] << " " << r_wall[1] << std::endl;
  }
 some_file.close();
 
 // Increment number of doc
 doc_info.number()++;

} // end of doc_solution


//=======start_of_unsteady_run============================================
/// Unsteady run
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void OscEllipseProblem<ELEMENT,TIMESTEPPER>::unsteady_run(DocInfo& doc_info)
{

 // Specify duration of the simulation
 double t_max=3.0;

 // Initial timestep
 double dt=0.025;

 // Initialise timestep
 initialise_dt(dt);

 // Set initial conditions.
 set_initial_condition();

 // Alternative initial conditions: impulsive start; see exercise.
 //assign_initial_values_impulsive(); 

 // find number of steps
 unsigned nstep = unsigned(t_max/dt);

 // If validation: Reduce number of timesteps performed and 
 // use coarse-ish mesh
 if (CommandLineArgs::Argc>1)
  {
   nstep=2;
   refine_uniformly();
   cout << "validation run" << std::endl;
  }
 else
  {
   // Refine the mesh three times, to resolve the pressure distribution
   // (the velocities could be represented accurately on a much coarser mesh).
   refine_uniformly();
   refine_uniformly();
   refine_uniformly();
  }

 // Output solution initial 
 doc_solution(doc_info);

 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   cout << "TIMESTEP " << istep << std::endl;
   cout << "Time is now " << time_pt()->time() << std::endl;

   // Take timestep 
   unsteady_newton_solve(dt);
     
   //Output solution
   doc_solution(doc_info);
  }

} // end of unsteady_run


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//======start_of_main=================================================
/// Driver code for unsteady Navier-Stokes flow, driven by
/// oscillating ellipse. If the code is executed with command line
/// arguments, a validation run is performed. 
//====================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);


 // Solve with Crouzeix-Raviart elements
 {
  // Create DocInfo object with suitable directory name for output
  DocInfo doc_info;
  doc_info.set_directory("RESLT_CR");
  
  //Set up problem
  OscEllipseProblem<RefineableQCrouzeixRaviartElement<2>,BDF<2> > problem;
  
  // Run the unsteady simulation
  problem.unsteady_run(doc_info);
 }

 // Solve with Taylor-Hood elements
 {
  // Create DocInfo object with suitable directory name for output
  DocInfo doc_info;
  doc_info.set_directory("RESLT_TH");

  //Set up problem
  OscEllipseProblem<RefineableQTaylorHoodElement<2>,BDF<2> > problem;
  
  // Run the unsteady simulation
  problem.unsteady_run(doc_info);
 }
 

}; // end of main



