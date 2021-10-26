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


namespace oomph
{

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
  if (std::abs(x[0])>1e-4) zeta=atan(x[1]/x[0]);

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
    std::cout << "Trafo error " << error << std::endl;
    //assert(false);
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
  if (std::abs(x[0])>1e-4) zeta=atan(x[1]/x[0]);

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
    std::cout << "Trafo error " << error << std::endl;
    //assert(false);
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


}

