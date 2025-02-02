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
// Driver for 2D Navier Stokes flow, driven by oscillating ring
// Leave "left" symmetry boundary condition unconstrained to 
// allow outflow.

//Generic includes
#include "generic.h"
#include "navier_stokes.h"

//Need to instantiate templated mesh
#include "meshes/quarter_circle_sector_mesh.h"


using namespace std;

using namespace oomph;
  

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
  

//==================================================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=100.0;  // ADJUST_PRIORITY 

 /// Reynolds x Strouhal number
 double ReSt=100.0; // ADJUST_PRIORITY 

}



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//=====================================================================
/// Sarah's boundary layer solution for flow in oscillating ring
//=====================================================================
namespace SarahBL
{

 double epsilon,alpha,A,Omega,N;



/* The options were    : operatorarrow */
#include <math.h>
double Diss_sarah(double rho,double zeta,double t)
{
  double t1;
  double t10;
  double t11;
  double t13;
  double t14;
  double t19;
  double t2;
  double t20;
  double t21;
  double t22;
  double t25;
  double t26;
  double t39;
  double t43;
  double t5;
  double t6;
  double t8;
  {
    t1 = alpha*alpha;
    t2 = epsilon*epsilon;
    t5 = sin(N*zeta);
    t6 = t5*t5;
    t8 = Omega*Omega;
    t10 = 1.0+A;
    t11 = t10*t10;
    t13 = sqrt(2.0);
    t14 = sqrt(Omega);
    t19 = t13*t14*alpha*(1.0-rho)/2.0;
    t20 = exp(-t19);
    t21 = t20*t20;
    t22 = Omega*t;
    t25 = cos(-t19+t22+0.3141592653589793E1/4.0);
    t26 = t25*t25;
    t39 = cos(-t19+t22);
    t43 = cos(t22);
    return(t1*t2*t6*t8*Omega*t11*t21*t26+4.0*alpha*t2*t6*t14*Omega*t10*t20*t25*
(N-1.0)*(Omega*t10*t20*t39/2.0-Omega*t43));
  }
}
/* The options were    : operatorarrow */
#include <math.h>
double Kin_energy_sarah(double t)
{
  double t1;
  double t11;
  double t13;
  double t16;
  double t18;
  double t2;
  double t20;
  double t27;
  double t28;
  double t31;
  double t33;
  double t5;
  double t6;
  double t7;
  {
    t1 = epsilon*epsilon;
    t2 = Omega*Omega;
    t5 = Omega*t;
    t6 = cos(t5);
    t7 = t6*t6;
    t11 = 1/alpha;
    t13 = sqrt(Omega);
    t16 = 1.0+A;
    t18 = 0.3141592653589793E1/4.0;
    t20 = sin(t5+t18);
    t27 = t16*t16;
    t28 = 1/t13;
    t31 = sin(2.0*t5+t18);
    t33 = sqrt(2.0);
    return(t1*(0.3141592653589793E1*t2/N*t7/2.0+t11*0.3141592653589793E1*t13*
Omega*t16*t6*t20)+0.3141592653589793E1*t1*t2*t11*(t27*(t28*t31/4.0+t33*t28/4.0
)-2.0*t28*t16*t6*t20)/2.0);
  }
}


/* The options were    : operatorarrow */
#include <math.h>
double P_sarah(double rho,double zeta,double t)
{
  double t1;
  double t12;
  double t19;
  double t4;
  double t6;
  double t8;
  double t9;
  {
    t1 = Omega*Omega;
    t4 = pow(rho,1.0*N);
    t6 = cos(N*zeta);
    t8 = Omega*t;
    t9 = sin(t8);
    t12 = sqrt(Omega);
    t19 = cos(t8+0.3141592653589793E1/4.0);
    return(epsilon*(t1/N*t4*t6*t9-t12*Omega*(1.0+A)*t4*t6*t19/alpha));
  }
}

/* The options were    : operatorarrow */
#include <math.h>
double Total_Diss_lead_sarah(double t)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t13;
  double t16;
  double t4;
  double t5;
  double t6;
  {
    t1 = epsilon*epsilon;
    t4 = sqrt(2.0);
    t5 = Omega*Omega;
    t6 = sqrt(Omega);
    t10 = pow(1.0+A,2.0);
    t11 = Omega*t;
    t12 = sin(t11);
    t13 = cos(t11);
    t16 = t13*t13;
    return(-alpha*t1*0.3141592653589793E1*t4*t6*t5*t10*(-1.0+2.0*t12*t13-2.0*
t16)/8.0);
  }
}

/* The options were    : operatorarrow */
#include <math.h>
double Total_Diss_sarah(double t)
{
  double t1;
  double t14;
  double t15;
  double t18;
  double t19;
  double t2;
  double t20;
  double t4;
  double t6;
  double t7;
  double t8;
  {
    t1 = Omega*Omega;
    t2 = 0.3141592653589793E1*t1;
    t4 = epsilon*epsilon;
    t6 = Omega*t;
    t7 = cos(t6);
    t8 = t7*t7;
    t14 = sqrt(2.0);
    t15 = sqrt(Omega);
    t18 = 1.0+A;
    t19 = t18*t18;
    t20 = sin(t6);
    return(4.0*t2*(N-1.0)*t4*t8-alpha*t4*0.3141592653589793E1*t14*t15*t1*t19*(
-1.0+2.0*t20*t7-2.0*t8)/8.0+t4*t18*t2*(-10.0*A*t8-A+8.0*A*N*t8+22.0*t8-1.0-24.0
*N*t8)/8.0);
  }
}

/* The options were    : operatorarrow */
#include <math.h>
double U_sarah(double rho,double zeta,double t)
{
  double t1;
  double t12;
  double t14;
  double t19;
  double t21;
  double t22;
  double t23;
  double t25;
  double t3;
  double t35;
  double t4;
  double t41;
  double t43;
  double t46;
  double t47;
  double t49;
  double t5;
  double t52;
  double t55;
  double t57;
  double t59;
  double t6;
  double t62;
  double t64;
  double t8;
  double t87;
  double t89;
  double t9;
  {
    t1 = cos(zeta);
    t3 = N-1.0;
    t4 = pow(rho,1.0*t3);
    t5 = N*zeta;
    t6 = cos(t5);
    t8 = Omega*t;
    t9 = cos(t8);
    t12 = sin(zeta);
    t14 = sin(t5);
    t19 = sqrt(Omega);
    t21 = 1.0+A;
    t22 = t21*t4;
    t23 = 0.3141592653589793E1/4.0;
    t25 = sin(t8+t23);
    t35 = 1/alpha;
    t41 = sqrt(2.0);
    t43 = 1.0-rho;
    t46 = t41*t19*alpha*t43/2.0;
    t47 = exp(-t46);
    t49 = cos(-t46+t8);
    t52 = Omega*t9;
    t55 = N*t19;
    t57 = sin(-t46+t8+t23);
    t59 = t47*t57-t25;
    t62 = Omega*alpha;
    t64 = -t43*t3;
    t87 = t9*(1.0+t64);
    t89 = N*t35;
    return(epsilon*(t1*Omega*t4*t6*t9+t12*Omega*t4*t14*t9+(t1*N*t19*t22*t6*t25+
t12*N*t19*t22*t14*t25)*t35)-epsilon*t12*(t14*(Omega*t21*t47*t49-t52)+t14*(t55*
t21*t59-t62*t64*t9)*t35)+epsilon*t1*(t52*t6+t6*(-t55*t21*t59-t62*t43*t3*t9)*t35
)-(-t12*(-Omega*t14*t87-t89*t19*t21*t14*t25)+t1*(Omega*t6*t87+t89*t6*t19*t21*
t25))*epsilon);
  }
}

/* The options were    : operatorarrow */
#include <math.h>
double V_sarah(double rho,double zeta,double t)
{
  double t1;
  double t12;
  double t14;
  double t19;
  double t21;
  double t22;
  double t23;
  double t25;
  double t3;
  double t35;
  double t4;
  double t41;
  double t43;
  double t46;
  double t47;
  double t49;
  double t5;
  double t52;
  double t55;
  double t57;
  double t59;
  double t6;
  double t62;
  double t64;
  double t8;
  double t87;
  double t89;
  double t9;
  {
    t1 = sin(zeta);
    t3 = N-1.0;
    t4 = pow(rho,1.0*t3);
    t5 = N*zeta;
    t6 = cos(t5);
    t8 = Omega*t;
    t9 = cos(t8);
    t12 = cos(zeta);
    t14 = sin(t5);
    t19 = sqrt(Omega);
    t21 = 1.0+A;
    t22 = t21*t4;
    t23 = 0.3141592653589793E1/4.0;
    t25 = sin(t8+t23);
    t35 = 1/alpha;
    t41 = sqrt(2.0);
    t43 = 1.0-rho;
    t46 = t41*t19*alpha*t43/2.0;
    t47 = exp(-t46);
    t49 = cos(-t46+t8);
    t52 = Omega*t9;
    t55 = N*t19;
    t57 = sin(-t46+t8+t23);
    t59 = t47*t57-t25;
    t62 = Omega*alpha;
    t64 = -t43*t3;
    t87 = t9*(1.0+t64);
    t89 = N*t35;
    return(epsilon*(t1*Omega*t4*t6*t9-t12*Omega*t4*t14*t9+(t1*N*t19*t22*t6*t25-
t12*N*t19*t22*t14*t25)*t35)+epsilon*t12*(t14*(Omega*t21*t47*t49-t52)+t14*(t55*
t21*t59-t62*t64*t9)*t35)+epsilon*t1*(t52*t6+t6*(-t55*t21*t59-t62*t43*t3*t9)*t35
)-(t12*(-Omega*t14*t87-t89*t19*t21*t14*t25)+t1*(Omega*t6*t87+t89*t6*t19*t21*t25
))*epsilon);
  }
}

/* The options were    : operatorarrow */
#include <math.h>
double X_sarah(double rho,double zeta,double t)
{
  double t1;
  double t10;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t1 = cos(zeta);
    t3 = sin(Omega*t);
    t5 = N*zeta;
    t6 = cos(t5);
    t8 = sin(t5);
    t10 = sin(zeta);
    return(rho*(t1+epsilon*t3*(t6*t1-A*t8*t10)));
  }
}

/* The options were    : operatorarrow */
#include <math.h>
double Y_sarah(double rho,double zeta,double t)
{
  double t1;
  double t10;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t1 = sin(zeta);
    t3 = sin(Omega*t);
    t5 = N*zeta;
    t6 = cos(t5);
    t8 = sin(t5);
    t10 = cos(zeta);
    return(rho*(t1+epsilon*t3*(t6*t1+A*t8*t10)));
  }
}



 /// Residual function for buckled ring
  void buckled_ring_residual(const Vector<double>& params,
                             const Vector<double>& unknowns,
                             Vector<double>& residuals)
 {
  // Parameters
  double t=params[0];
  double xx=params[1];
  double yy=params[2];

  // Unknowns
  double rho=unknowns[0];
  double zeta=unknowns[1];

  // Residuals
  residuals[0]=X_sarah(rho,zeta,t)-xx;
  residuals[1]=Y_sarah(rho,zeta,t)-yy;
 }



 /// Exact solution: x,y,u,v,p
 void exact_soln(const double& time,
                 const Vector<double>& x, 
                 Vector<double>& soln)
 {
 
  // Primary variables only
  soln.resize(3);

  // Get rho and zeta with black-box Newton method:
 
  // Parameters: time, x, y
  Vector<double> params(3);
  params[0]=time;
  params[1]=x[0];
  params[2]=x[1];
  

  // Unknowns rho, zeta
  Vector<double> unknowns(2);

  // Initial guess:
  double rho=sqrt(x[0]*x[0]+x[1]*x[1]);
  double zeta=0.5*MathematicalConstants::Pi;
  if (abs(x[0])>1e-4) zeta=atan(x[1]/x[0]);

  // Copy across
  unknowns[0]=rho;
  unknowns[1]=zeta;

  // Call Newton solver
  BlackBoxFDNewtonSolver::black_box_fd_newton_solve(
   buckled_ring_residual,params,unknowns);
  
  // Extract rho, zeta

  // Copy across
  rho=unknowns[0];
  zeta=unknowns[1];


  // Double check the position
  double dx=X_sarah(rho,zeta,time)-x[0]; 
  double dy=Y_sarah(rho,zeta,time)-x[1]; 
  double error=sqrt(dx*dx+dy*dy);
  if (error>1.0e-6)
   {
    std::ostringstream warning_stream;
    warning_stream << "Trafo error " << error << std::endl;

    OomphLibWarning(warning_stream.str(),
                    "Sarah_BL::exact_soln()",
                    OOMPH_EXCEPTION_LOCATION);
   }

  // Velocities
  soln[0]=U_sarah(rho,zeta,time); 
  soln[1]=V_sarah(rho,zeta,time); 
  
  // Pressure on viscous scale
  soln[2]=P_sarah(rho,zeta,time)*(alpha*alpha);
  
 }
 


 
 /// Full exact solution: x,y,u,v,p,du/dt,dv/dt,diss
  void full_exact_soln(const double& time,
                      const Vector<double>& x, 
                      Vector<double>& soln)
 {
 
  // Full solution
  soln.resize(6);

  // Get rho and zeta with black-box Newton method:
  
  // Parameters: time, x, y
  Vector<double> params(3);
  params[0]=time;
  params[1]=x[0];
  params[2]=x[1];
  

  // Unknowns rho, zeta
  Vector<double> unknowns(2);

  // Initial guess:
  double rho=sqrt(x[0]*x[0]+x[1]*x[1]);
  double zeta=0.5*MathematicalConstants::Pi;
  if (abs(x[0])>1e-4) zeta=atan(x[1]/x[0]);

  // Copy across
  unknowns[0]=rho;
  unknowns[1]=zeta;

  // Call Newton solver
  BlackBoxFDNewtonSolver::black_box_fd_newton_solve(
   buckled_ring_residual,params,unknowns);
  
  // Extract rho, zeta

  // Copy across
  rho=unknowns[0];
  zeta=unknowns[1];


  // Double check the position
  double dx=X_sarah(rho,zeta,time)-x[0]; 
  double dy=Y_sarah(rho,zeta,time)-x[1]; 
  double error=sqrt(dx*dx+dy*dy);
  if (error>1.0e-6)
   {
    std::ostringstream warning_stream;
    warning_stream << "Trafo error " << error << std::endl;
    
    OomphLibWarning(warning_stream.str(),
                    "Sarah_BL::full_exact_soln()",
                    OOMPH_EXCEPTION_LOCATION);
   }

  // Velocities
  soln[0]=U_sarah(rho,zeta,time); 
  soln[1]=V_sarah(rho,zeta,time); 
  
  // Pressure on viscous scale
  soln[2]=P_sarah(rho,zeta,time)*(alpha*alpha);
 
  // du/dt dv/dt not coded up yet...
  soln[3]=0.0; //DUDT(xx,y,time);
  soln[4]=0.0; //DVDT(xx,y,time);

  if (rho<1.0e-6)
   {
    rho=1e-6;
   }
  soln[5]=Diss_sarah(rho,zeta,time);


 }
 
}






/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//====================================================================
/// Driver for oscillating ring problem: Wall performs oscillations
/// that resemble eigenmodes of freely oscillating ring and drives
/// viscous fluid flow -- we allow for outflow through the
/// "left" symmetry plane to accomodate mass conservation...
//====================================================================
template<class ELEMENT, class TIMESTEPPER>
class OscRingNStProblem : public Problem
{

public:

 /// Constructor: Pass timestep and function 
 /// pointer to the solution that provides the initial conditions for the fluid
 OscRingNStProblem(const double& dt,
                   FiniteElement::UnsteadyExactSolutionFctPt IC_fct_pt);

 /// Destructor (empty) 
 ~OscRingNStProblem(){}

 /// Get pointer to wall as geometric object
 GeomObject* wall_pt()
  {
   return Wall_pt;
  }

 /// Update after solve (empty)
 void actions_after_newton_solve(){}

 /// Update before solve (empty)
 void actions_before_newton_solve(){}


 /// Update the problem specs before next timestep:
 /// Update the fluid mesh and re-set velocity boundary conditions. 
 void actions_before_implicit_timestep()
 {
  // Update the fluid mesh -- auxiliary update function for algebraic
  // nodes automatically updates no slip condition.
  mesh_pt()->node_update(); 

 // Ring boundary: No slip; this also implies that the velocity needs
 // to be updated in response to wall motion
 unsigned ibound=1;
  {
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Which node are we dealing with?
     Node* node_pt=mesh_pt()->boundary_node_pt(ibound,inod);

     // Apply no slip
     FSI_functions::apply_no_slip_on_moving_wall(node_pt);
    }
  }
 }

 /// Run the time integration for ntsteps steps
 void unsteady_run(const unsigned &ntsteps, const bool& restarted,
                   DocInfo& doc_info);

 /// Set initial condition (incl previous timesteps) according
 /// to specified function. 
 void set_initial_condition();
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Access function for the fluid mesh
 RefineableQuarterCircleSectorMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableQuarterCircleSectorMesh<ELEMENT>*>(
    Problem::mesh_pt()); 
  }
 
 /// Dump generic problem data.
 void dump_it(ofstream& dump_file);

 /// Read generic problem data.
 void restart(ifstream& restart_file);

private:

 /// Write header for trace file
 void write_trace_file_header();

 /// Function pointer to set the intial condition
 FiniteElement::UnsteadyExactSolutionFctPt IC_Fct_pt;

 /// Pointer to wall
 GeomObject* Wall_pt;

 /// Trace file
 ofstream Trace_file;

 /// Pointer to node on coarsest mesh on which velocity is traced
 Node* Veloc_trace_node_pt;

 /// Pointer to node in symmetry plane on coarsest mesh at 
 /// which velocity is traced
 Node* Sarah_veloc_trace_node_pt;

};



//========================================================================
/// Pass initial time, initial timestep and function 
/// pointer to the solution that provides the initial conditions for
/// the fluid.
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
OscRingNStProblem<ELEMENT, TIMESTEPPER>::OscRingNStProblem(
 const double& dt, FiniteElement::UnsteadyExactSolutionFctPt IC_fct_pt) : 
 IC_Fct_pt(IC_fct_pt)
{ 

 // Period of oscillation
 double T=1.0;

 //Allocate the timestepper
 add_time_stepper_pt(new TIMESTEPPER());

 // Initialise timestep
 time_pt()->initialise_dt(dt);

 // Initialise time
 time_pt()->time()=0.0;

 // Parameters for pseudo-buckling ring 
 double eps_buckl=0.1;  // ADJUST_PRIORITY 
 double ampl_ratio=-0.5;  // ADJUST_PRIORITY 
 unsigned n_buckl=2; // ADJUST_PRIORITY 
 double r_0=1.0;
 
 // Build wall geometric element
 Wall_pt=new PseudoBucklingRing(eps_buckl,ampl_ratio,n_buckl,r_0,T,
                                time_stepper_pt());

 // Fluid mesh is suspended from wall between these two Lagrangian
 // coordinates:
 double xi_lo=0.0;
 double xi_hi=2.0*atan(1.0);

 // Fractional position of dividing line for two outer blocks in fluid mesh
 double fract_mid=0.5;


 // Build fluid mesh...
 Problem::mesh_pt()=new RefineableQuarterCircleSectorMesh<ELEMENT>(
  Wall_pt,xi_lo,fract_mid,xi_hi,time_stepper_pt());

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
  
 // Extract pointer to node at center of mesh (this node exists
 // at all refinement levels and can be used to doc continuous timetrace
 // of velocities)
 unsigned nnode=mesh_pt()->finite_element_pt(0)->nnode();
 Veloc_trace_node_pt=mesh_pt()->finite_element_pt(0)->node_pt(nnode-1);

 //  Extract pointer to node in symmetry plane (this node exists
 // at all refinement levels and can be used to doc continuous timetrace
 // of velocities)
 unsigned nnode_1d=dynamic_cast<ELEMENT*>(
  mesh_pt()->finite_element_pt(0))->nnode_1d();
 Sarah_veloc_trace_node_pt=mesh_pt()->
  finite_element_pt(0)->node_pt(nnode_1d-1);


 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 
 // Bottom boundary: 
 unsigned ibound=0;
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Pin vertical velocity
     {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1); 
     }
    }
  }
  
 // Ring boundary: No slip; this also implies that the velocity needs
 // to be updated in response to wall motion
 ibound=1;
  {
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Which node are we dealing with?
     Node* node_pt=mesh_pt()->boundary_node_pt(ibound,inod);

     // Pin both velocities
     for (unsigned i=0;i<2;i++)
      {
       node_pt->pin(i); 
      }
    }
  }

  // Left boundary is left free to allow for mass conservation
                        
 //Find out how many timesteppers there are
 unsigned Ntime_steppers = ntime_stepper();

 //Loop over them all and set the weights
 for(unsigned i=0;i<Ntime_steppers;i++)
  {
   time_stepper_pt(i)->set_weights();
  }


 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned Nelement = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<Nelement;i++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
  }


 //Attach the boundary conditions to the mesh
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
 // Set parameters for Sarah's asymptotic solution
 SarahBL::epsilon=static_cast<PseudoBucklingRing*>(Wall_pt)->eps_buckl();
 SarahBL::alpha=sqrt(Global_Physical_Variables::Re);
 SarahBL::A=static_cast<PseudoBucklingRing*>(Wall_pt)->ampl_ratio();
 SarahBL::N=static_cast<PseudoBucklingRing*>(Wall_pt)->n_buckl_float();
 SarahBL::Omega=2.0*MathematicalConstants::Pi;

}




//========================================================================
/// Set initial condition: Assign previous and current values
/// of the velocity from the velocity field specified via
/// the function pointer.
///
/// Values are assigned so that the velocities and accelerations
/// are correct for current time.
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void OscRingNStProblem<ELEMENT, TIMESTEPPER>::set_initial_condition()
{ 

 // Backup time in global timestepper
 double backed_up_time=time_pt()->time();
         
 // Past history for velocities needs to be established for t=time0-deltat, ...
 // Then provide current values (at t=time0) which will also form
 // the initial guess for first solve at t=time0+deltat

 // Vector of exact solution values (includes pressure)
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
 // the mesh might also move)
 for (int itime=nprev_steps;itime>=0;itime--)
  {
   double time=prev_time[itime];
   
   // Set global time (because this is how the geometric object refers 
   // to continous time 
   time_pt()->time()=time;
   
   cout << "setting IC at time =" << time << std::endl;
   
   // Update the fluid mesh for this value of the continuous time
   // (The wall object reads the continous time from the same
   // global time object)
   mesh_pt()->node_update(); 
   
   // Loop over the nodes to set initial guess everywhere
   for (unsigned jnod=0;jnod<num_nod;jnod++)
    {
     
     // Get nodal coordinates
     x[0]=mesh_pt()->node_pt(jnod)->x(0);
     x[1]=mesh_pt()->node_pt(jnod)->x(1);

     // Get intial solution
     (*IC_Fct_pt)(time,x,soln);
     
     // Loop over velocity directions (velocities are in soln[0] and soln[1]).
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->node_pt(jnod)->set_value(itime,i,soln[i]);
      }
     
     // Loop over coordinate directions
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->node_pt(jnod)->x(itime,i)=x[i];
      }
    } 
  }

 // Reset backed up time for global timestepper
 time_pt()->time()=backed_up_time;

}





//========================================================================
/// Doc the solution
///
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void OscRingNStProblem<ELEMENT, TIMESTEPPER>::doc_solution(DocInfo& doc_info)
{ 

 cout << "Doc-ing step " <<  doc_info.number()
      << " for time " << time_stepper_pt()->time_pt()->time() << std::endl;


 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 


 // Output solution on fluid mesh
 //-------------------------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 //some_file.precision(20);
 some_file.open(filename);
 unsigned nelem=mesh_pt()->nelement();
 for (unsigned ielem=0;ielem<nelem;ielem++)
  {
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(ielem))->
    full_output(some_file,npts);
  }
 some_file.close();
  

 // Plot wall posn
 //---------------
 sprintf(filename,"%s/Wall%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 
 unsigned nplot=100;
 Vector<double > xi_wall(1);
 Vector<double > r_wall(2);
  for (unsigned iplot=0;iplot<nplot;iplot++)
   {
    xi_wall[0]=0.5*MathematicalConstants::Pi*double(iplot)/double(nplot-1);
    Wall_pt->position(xi_wall,r_wall);
    some_file << r_wall[0] << " " << r_wall[1] << std::endl;
   }
  some_file.close();


 // Doc Sarah's asymptotic solution 
 //--------------------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,
                       time_stepper_pt()->time_pt()->time(), 
                       SarahBL::full_exact_soln);
 some_file.close();




 // Get position of control point
 //------------------------------
 Vector<double> r(2);
 Vector<double> xi(1);
 xi[0]=MathematicalConstants::Pi/2.0;
 wall_pt()->position(xi,r);



 // Get total volume (area) of computational domain, energies and average
 //----------------------------------------------------------------------
 // pressure
 //---------
 double area=0.0;
 double press_int=0.0;
 double diss=0.0;
 double kin_en=0.0;
 nelem=mesh_pt()->nelement();

 for (unsigned ielem=0;ielem<nelem;ielem++)
  {
   area+=mesh_pt()->finite_element_pt(ielem)->size();
   press_int+=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(ielem))
    ->pressure_integral();
   diss+=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(ielem))->
    dissipation();
   kin_en+=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(ielem))->
    kin_energy();
  }

 // Total kinetic energy in entire domain
 double global_kin=4.0*kin_en;

 // Max/min refinement level
 unsigned min_level;
 unsigned max_level;
 mesh_pt()->get_refinement_levels(min_level,max_level);


 // Total dissipation for quarter domain
 double time=time_pt()->time();
 double diss_asympt=SarahBL::Total_Diss_sarah(time)/4.0;

 // Asymptotic predictions for velocities at control point
 Vector<double> x_sarah(2);
 Vector<double> soln_sarah(3);
 x_sarah[0]=Sarah_veloc_trace_node_pt->x(0);
 x_sarah[1]=Sarah_veloc_trace_node_pt->x(1);
 SarahBL::exact_soln(time,x_sarah,soln_sarah);

 // Doc
 Trace_file << time_pt()->time() 
            << " " << r[1]
            << " " << global_kin 
            << " " << SarahBL::Kin_energy_sarah(time_pt()->time()) 
            << " " << static_cast<PseudoBucklingRing*>(Wall_pt)->r_0()
            << " " << area  
            << " " << press_int/area 
            << " " << diss  
            << " " << diss_asympt
            << " " << Veloc_trace_node_pt->x(0)
            << " " << Veloc_trace_node_pt->x(1) 
            << " " << Veloc_trace_node_pt->value(0)
            << " " << Veloc_trace_node_pt->value(1) 
            << " " << mesh_pt()->nelement() 
            << " " << ndof()
            << " " << min_level 
            << " " << max_level
            << " " << mesh_pt()->nrefinement_overruled() 
            << " " << mesh_pt()->max_error()  
            << " " << mesh_pt()->min_error()  
            << " " << mesh_pt()->max_permitted_error() 
            << " " << mesh_pt()->min_permitted_error() 
            << " " << mesh_pt()->max_keep_unrefined()
            << " " << doc_info.number()
            << " " << Sarah_veloc_trace_node_pt->x(0) 
            << " " << Sarah_veloc_trace_node_pt->x(1) 
            << " " << Sarah_veloc_trace_node_pt->value(0) 
            << " " << Sarah_veloc_trace_node_pt->value(1)
            << " " << x_sarah[0]
            << " " << x_sarah[1]
            << " " << soln_sarah[0]
            << " " << soln_sarah[1]
            << " " 
            << static_cast<PseudoBucklingRing*>(Wall_pt)->r_0()-1.0
            << std::endl;


 // Output fluid solution on coarse-ish mesh 
 //-----------------------------------------

 // Extract all elements from quadtree representation
 Vector<Tree*> all_element_pt;
 mesh_pt()->forest_pt()->
  stick_all_tree_nodes_into_vector(all_element_pt);

 // Build a coarse mesh
 Mesh* coarse_mesh_pt = new Mesh();

 //Loop over all elements and check if their refinement level is
 //equal to the mesh's minimum refinement level
 nelem=all_element_pt.size();
 for (unsigned ielem=0;ielem<nelem;ielem++)
  {
   Tree* el_pt=all_element_pt[ielem];
   if (el_pt->level()==min_level)
    {
     coarse_mesh_pt->add_element_pt(el_pt->object_pt());
    }
  }

 // Output fluid solution on coarse mesh
 sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 nelem=coarse_mesh_pt->nelement();
 for (unsigned ielem=0;ielem<nelem;ielem++)
  {
   dynamic_cast<ELEMENT*>(coarse_mesh_pt->element_pt(ielem))->
    full_output(some_file,npts);
  }
 some_file.close();

 // Write restart file
 sprintf(filename,"%s/restart%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 dump_it(some_file);
 some_file.close();

}



//========================================================================
/// Dump the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void OscRingNStProblem<ELEMENT, TIMESTEPPER>::dump_it(ofstream& dump_file)
{

 // Dump refinement status of mesh
 //mesh_pt()->dump_refinement(dump_file);

 // Call generic dump()
 Problem::dump(dump_file);
  
}



//========================================================================
/// Read solution from disk
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void OscRingNStProblem<ELEMENT, TIMESTEPPER>::restart(ifstream& restart_file)
{
 // Refine fluid mesh according to the instructions in restart file
 //mesh_pt()->refine(restart_file);

 // Read generic problem data
 Problem::read(restart_file);

 //Assign equation numbers
 //cout <<"\nNumber of equations: " << assign_eqn_numbers() 
 //     << std::endl<< std::endl; 
}



//========================================================================
/// Driver for timestepping the problem: Fixed timestep but 
/// guaranteed spatial accuracy. Beautiful, innit?
///
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void OscRingNStProblem<ELEMENT, TIMESTEPPER>::unsteady_run(
                                               const unsigned& ntsteps,
                                               const bool& restarted,
                                               DocInfo& doc_info)
{ 

 // Change max. residual for Newton iteration
 Max_residuals=1000;

 // Step number
 doc_info.number()=0;

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 //Use constant timestep for entire simulation
 double dt = time_pt()->dt();

 // Max. number of adaptations per solve
 unsigned max_adapt;

 // Max. number of adaptations per solve
 if (restarted)
  {
   max_adapt=0;
  }
 else
  {
   max_adapt=1;
  }
 
 // Max. and min. error for adaptive refinement/unrefinement
 mesh_pt()->max_permitted_error()=  0.5e-2;   
 mesh_pt()->min_permitted_error()=  0.5e-3;  

 // Don't allow refinement to drop under given level
 mesh_pt()->min_refinement_level()=2;

 // Don't allow refinement beyond given level 
 mesh_pt()->max_refinement_level()=6; 

 // Don't bother adapting the mesh if no refinement is required
 // and if less than ... elements are to be merged.
 mesh_pt()->max_keep_unrefined()=20;


 // Get max/min refinement levels in mesh
 unsigned min_refinement_level;
 unsigned max_refinement_level;
 mesh_pt()->get_refinement_levels(min_refinement_level,
                                        max_refinement_level);

 cout << "\n Initial mesh: min/max refinement levels: " 
      << min_refinement_level << " " << max_refinement_level << std::endl << std::endl;

 // Doc refinement targets
 mesh_pt()->doc_adaptivity_targets(cout);

 // Write header for trace file
 write_trace_file_header();

 // Doc initial condition
 doc_solution(doc_info);
 doc_info.number()++;

 // Switch off doc during solve
 doc_info.disable_doc();

 // If we set first to true, then initial guess will be re-assigned
 // after mesh adaptation. Don't want this if we've done a restart.
 bool first;
 bool shift;
 if (restarted)
  {
   first=false;
   shift=false;
   // Move time back by dt to make sure we're re-solving the read-in solution
   time_pt()->time()-=time_pt()->dt();
  }
 else
  {
   first=true;
   shift=true;
  }
 
 //Take the first fixed timestep with specified tolerance for Newton solver
 unsteady_newton_solve(dt,max_adapt,first,shift);

 // Switch doc back on
 doc_info.enable_doc();

 // Doc initial solution
 doc_solution(doc_info);
 doc_info.number()++;

 // Now do normal run; allow for one mesh adaptation per timestep
 max_adapt=1;

 //Loop over the remaining timesteps
 for(unsigned t=2;t<=ntsteps;t++)
  {

   // Switch off doc during solve
   doc_info.disable_doc();

   //Take fixed timestep
   first=false;
   unsteady_newton_solve(dt,max_adapt,first);

   // Switch doc back on
   doc_info.enable_doc();

   // Doc solution
   //if (icount%10==0)
    {
     doc_solution(doc_info);
     doc_info.number()++;
    }
  }

 // Close trace file
 Trace_file.close();

}




//========================================================================
/// Write trace file header
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void OscRingNStProblem<ELEMENT,TIMESTEPPER>::write_trace_file_header()
{

 // Doc parameters in trace file
 Trace_file << "# err_max " << mesh_pt()->max_permitted_error() << std::endl;
 Trace_file << "# err_min " << mesh_pt()->min_permitted_error() << std::endl;
 Trace_file << "# Re " << Global_Physical_Variables::Re << std::endl;
 Trace_file << "# St " << Global_Physical_Variables::ReSt/
                          Global_Physical_Variables::Re << std::endl;
 Trace_file << "# dt " << time_stepper_pt()->time_pt()->dt() << std::endl;
 Trace_file << "# initial # elements " << mesh_pt()->nelement() << std::endl;
 Trace_file << "# min_refinement_level " 
            << mesh_pt()->min_refinement_level() << std::endl;
 Trace_file << "# max_refinement_level " 
            << mesh_pt()->max_refinement_level() << std::endl;



 Trace_file << "VARIABLES=\"time\",\"V_c_t_r_l\",\"e_k_i_n\"";
 Trace_file << ",\"e_k_i_n_(_a_s_y_m_p_t_)\",\"R_0\",\"Area\"" ;
 Trace_file << ",\"Average pressure\",\"Total dissipation (quarter domain)\"";
 Trace_file << ",\"Asymptotic dissipation (quarter domain)\"";
 Trace_file << ",\"x<sub>1</sub><sup>(track)</sup>\"";
 Trace_file << ",\"x<sub>2</sub><sup>(track)</sup>\"";
 Trace_file << ",\"u<sub>1</sub><sup>(track)</sup>\"";
 Trace_file << ",\"u<sub>2</sub><sup>(track)</sup>\"";
 Trace_file << ",\"N<sub>element</sub>\"";
 Trace_file << ",\"N<sub>dof</sub>\"";
 Trace_file << ",\"max. refinement level\"";
 Trace_file << ",\"min. refinement level\"";
 Trace_file << ",\"# of elements whose refinement was over-ruled\"";
 Trace_file << ",\"max. error\"";
 Trace_file << ",\"min. error\"";
 Trace_file << ",\"max. permitted error\"";
 Trace_file << ",\"min. permitted error\"";
 Trace_file << ",\"max. permitted # of unrefined elements\"";
 Trace_file << ",\"doc number\"";
 Trace_file << ",\"x<sub>1</sub><sup>(track2 FE)</sup>\"";
 Trace_file << ",\"x<sub>2</sub><sup>(track2 FE)</sup>\"";
 Trace_file << ",\"u<sub>1</sub><sup>(track2 FE)</sup>\"";
 Trace_file << ",\"u<sub>2</sub><sup>(track2 FE)</sup>\"";
 Trace_file << ",\"x<sub>1</sub><sup>(track2 Sarah)</sup>\"";
 Trace_file << ",\"x<sub>2</sub><sup>(track2 Sarah)</sup>\"";
 Trace_file << ",\"u<sub>1</sub><sup>(track2 Sarah)</sup>\"";
 Trace_file << ",\"u<sub>2</sub><sup>(track2 Sarah)</sup>\"";
 Trace_file << ",\"R<sub>0</sub>-1\"";
 Trace_file << std::endl;

}






/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////




//========================================================================
/// Demonstrate how to solve OscRingNSt problem in deformable domain
/// with mesh adaptation
//========================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 //Do a certain number of timesteps per period
 unsigned nstep_per_period=40; // 80;  // ADJUST_PRIORITY
 unsigned nperiod=3;

 // Work out total number of steps and timestep
 unsigned nstep=nstep_per_period*nperiod;  
 double dt=1.0/double(nstep_per_period);

 // Set up the problem: Pass timestep and Sarah's asymptotic solution for
 // generation of initial condition
 OscRingNStProblem<RefineableQCrouzeixRaviartElement<2>,
  BDF<4> > problem(dt,&SarahBL::full_exact_soln);


 // Restart?
 //---------
 bool restarted=false;

 // Pointer to restart file
 ifstream* restart_file_pt=0;

 // No restart
 //-----------
 if (CommandLineArgs::Argc!=2)
  {
   cout << "No restart" << std::endl;
   restarted=false;

   // Refine uniformly
   problem.refine_uniformly();
   problem.refine_uniformly();
   problem.refine_uniformly();

   // Set initial condition on uniformly refined domain (if this isn't done
   // then the mesh contains the interpolation of the initial guess
   // on the original coarse mesh -- if the nodal values were zero in
   // the interior and nonzero on the boundary, then the the waviness of
   // of the interpolated i.g. between the nodes on the coarse mesh
   // gets transferred onto the fine mesh where we can do better
   problem.set_initial_condition();
  }
 
 // Restart (single command line argument)
 //---------------------------------------
 else if (CommandLineArgs::Argc==2)
  {
   restarted=true;

   // Open restart file
   restart_file_pt=new ifstream(CommandLineArgs::Argv[1],ios_base::in);
   if (restart_file_pt!=0)
    {
     cout << "Have opened " << CommandLineArgs::Argv[1] << 
      " for restart. " << std::endl;
    }
   else
    {
     std::ostringstream error_stream;
     error_stream 
      << "ERROR while trying to open " << CommandLineArgs::Argv[1] << 
      " for restart." << std::endl;

     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   // Do the actual restart
   problem.restart(*restart_file_pt);
  }


 // Two command line arguments: do validation run
 //----------------------------------------------
 if (CommandLineArgs::Argc==3)
  {
    nstep=3;
    cout << "Only doing nstep steps for validation: " << nstep << std::endl;
  }


 // Setup labels for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT"); 


 // Do unsteady run of the problem for nstep steps
 //-----------------------------------------------
 problem.unsteady_run(nstep,restarted,doc_info);



 // Validate the restart procedure
 //-------------------------------
 if (CommandLineArgs::Argc==3)
  {
   
   // Build coarser problem
   
   // Set up the problem: Pass timestep and Sarah's asymptotic solution for
   // generation of initial condition
   OscRingNStProblem<RefineableQCrouzeixRaviartElement<2>,
    BDF<4> > restarted_problem(dt,&SarahBL::full_exact_soln);
   
   // Setup labels for output
   DocInfo restarted_doc_info;
   
   // Output directory
   restarted_doc_info.set_directory("RESLT_restarted"); 
   
   // Step number
   restarted_doc_info.number()=0;
   
   // Copy by performing a restart from old problem
   restart_file_pt=new ifstream("RESLT/restart2.dat");
   
   // Do the actual restart
   restarted_problem.restart(*restart_file_pt);
   
   // Do unsteady run of the problem for one step
   unsigned nstep=2;
   bool restarted=true;
   restarted_problem.unsteady_run(nstep,restarted,restarted_doc_info);

  }
}










