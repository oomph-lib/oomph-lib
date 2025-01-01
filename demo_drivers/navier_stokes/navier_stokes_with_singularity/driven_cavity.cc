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
//Driver for 2D rectangular driven cavity


//Generic includes
#include "generic.h"
#include "navier_stokes.h"

#include "meshes/simple_rectangular_quadmesh.h"
#include "navier_stokes_elements_with_singularity.h"

using namespace std;

using namespace oomph;

/// ///////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////


//===============================================================
/// Overloaded Taylor Hood element with compute_error fct that 
/// includes the pressure
//===============================================================
class MyTaylorHoodElement : public virtual QTaylorHoodElement<2>
{

public:



 /// Overloaded compute error function; uses FE+singular parts
void compute_error(std::ostream &outfile,
                   FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                   const bool& include_pressure,
                   double& error, double& norm)
{

 unsigned cached_dim=this->dim();  

 error=0.0;
 norm=0.0;

 //Vector of local coordinates
 Vector<double> s(cached_dim);

 // Vector for coordintes
 Vector<double> x(cached_dim);

 //Set the value of n_intpt
 unsigned n_intpt = this->integral_pt()->nweight();
   

 outfile << "ZONE" << std::endl;
 
 // Exact solution Vector (u,v,[w],p)
 Vector<double> exact_soln(cached_dim+1);
 Vector<double> computed_soln(cached_dim+1);
   
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<cached_dim;i++)
    {
     s[i] = this->integral_pt()->knot(ipt,i);
    }

   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);

   // Get jacobian of mapping
   double J=this->J_eulerian(s);

   //Premultiply the weights and the Jacobian
   double W = w*J;

   // Get x position as Vector
   this->interpolated_x(s,x);

   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   for(unsigned i=0;i<cached_dim;i++)
    {
     computed_soln[i]=interpolated_u_nst(s,i);
    }
   unsigned hi_limit=cached_dim;
   if (include_pressure)
    {
     computed_soln[cached_dim]=interpolated_p_nst(s);
     hi_limit=cached_dim+1;
    }



   // Velocity error
   for(unsigned i=0;i<hi_limit;i++)
    {
     norm+=exact_soln[i]*exact_soln[i]*W;
     error+=(exact_soln[i]-computed_soln[i])*
      (exact_soln[i]-computed_soln[i])*W;
    }

   //Output x,y,...,u_exact
   for(unsigned i=0;i<cached_dim;i++)
    {
     outfile << x[i] << " ";
    }

   //Output x,y,[z],u_error,v_error,[w_error],[p_error[
   for(unsigned i=0;i<hi_limit;i++)
    {
     outfile << exact_soln[i]-computed_soln[i] << " ";
    }
   outfile << std::endl;   
  }
}


};





/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////




//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=100.0;

} // end_of_namespace


#include "finite_re_perturbation.h"



//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Parameters
{
 // Multiplier for the number of elements
 unsigned Element_multiplier = 1;
 
 // Basic # of elements in x-direction
 unsigned N_x=10;

 // Basic # of elements in y-direction
 unsigned N_y=10;
 
 // Domain length in x-direction
 double L_x=1.0;

 // Domain length in y-direction
 double L_y=1.0;
 
 // Pointer to the direction at which the derivative used in the residuals
 // of C1 and C2 will be computed
 unsigned Direction = 1;

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

 /// Inner blending radius
 double R_blend_inner=0.2; // 0.1;

 /// Outer blending radius
 double R_blend_outer=0.8; // 0.2;

 /// Blend function, centred at x_c
 double blend(const Vector<double>& x, const Vector<double>& x_c)
 {
  // shifted origin
  double dx=x[0]-x_c[0];
  double dy=x[1]-x_c[1];
  double r=sqrt(dx*dx+dy*dy);
  double result=0.0;
  if (r<=R_blend_inner)
   {
    result=1.0;
   }
  else if (r<=R_blend_outer)
   {
    result=0.5*(1.0+cos(MathematicalConstants::Pi*
                        (r-R_blend_inner)/(R_blend_outer-R_blend_inner)));
   }
  return result;
 }


 /// Gradient of blend function, centred at x_c
 Vector<double> grad_blend(const Vector<double>& x, const Vector<double>& x_c)
 {
  // shifted origin
  double dx=x[0]-x_c[0];
  double dy=x[1]-x_c[1];
  double r=sqrt(dx*dx+dy*dy);
  Vector<double> result(2,0.0);
  if ((r>R_blend_inner)&&(r<R_blend_outer))
   {
    result[0] = -sin(MathematicalConstants::Pi*
                     (sqrt(dx*dx+dy*dy)-R_blend_inner)/(
                      R_blend_outer-R_blend_inner))*
     MathematicalConstants::Pi/sqrt(dx*dx+dy*dy)*dx/(
      R_blend_outer-R_blend_inner)/2.0;

    result[1] = -sin(MathematicalConstants::Pi*
                      (sqrt(dx*dx+dy*dy)-R_blend_inner)/(
                       R_blend_outer-R_blend_inner))*
     MathematicalConstants::Pi/sqrt(dx*dx+dy*dy)*dy/(
      R_blend_outer-R_blend_inner)/2.0;
   }
  return result;
 }



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


 // Circular Couette flow
 //----------------------

 /// Amplitude of linearly varying part
 double A_couette=1.0;


 /// Amplitude of singular part
 double B_couette=1.0;


 /// Couette flow
 Vector<double> velocity_couette(const Vector<double>& x)
 {
  // Initialise velocity vector to return
  Vector<double> u(2,0.0);

  // Shifted origin
  double x_shifted=1.0+x[0];
  double y_shifted=1.0+x[1];


  // Polar coordinates
  double theta = atan2(y_shifted,x_shifted);
  double r=sqrt(x_shifted*x_shifted+y_shifted*y_shifted);

  double magnitude=A_couette*r+B_couette/r;
  u[0]=-magnitude*sin(theta);
  u[1]= magnitude*cos(theta);

  return u;
 }

 /// Couette Pressure
 double pressure_couette(const Vector<double>& x)
 {
  // Shifted origin
  double x_shifted=1.0+x[0];
  double y_shifted=1.0+x[1];

  // Polar coordinates
  double r=sqrt(x_shifted*x_shifted+y_shifted*y_shifted);

  // Pressure
  double press=Global_Physical_Variables::Re*
   (A_couette*A_couette*r*r/2.0+
    2.0*A_couette*B_couette*log(r)-
    B_couette*B_couette/(2.0*r*r));

  return press;
 }


 /// Combined exact solution (u,v,p)
 void exact_couette(const Vector<double>& x, Vector<double>& soln)
 {
  soln=velocity_couette(x);
  soln.resize(3);
  soln[2]=pressure_couette(x);            
 }


 /// Pseudo_singularity for computation of Couette flow
 Vector<double> velocity_pseudo_singularity_for_couette(const Vector<double>& x)
 {
  // Initialise velocity vector to return
  Vector<double> u(2,0.0);

  // Poiseuille flow
  u[0]=0.0;
  u[1]=x[0]*(1.0-x[0]);


  return u;
 }


 /// Pseudo_singularity for computation of Couette flow: grad[i][j] = du_i/dx_j
 Vector<Vector<double> > grad_velocity_pseudo_singularity_for_couette(
  const Vector<double>& x)
 {
  // Initialise grad matrix
  Vector<Vector<double> > grad_u(2);
  for (unsigned d=0;d<2;d++)
   {
    grad_u[d].resize(2);
    for (unsigned k=0;k<2;k++)
     {
      grad_u[d][k] = 0.0;
     }
   }

  // Compute \frac{\partial u}{\partial x} 
  grad_u[0][0] = 0.0;

  // Compute \frac{\partial u}{\partial y}
  grad_u[0][1] = 0.0; // 1.0-2.0*x[1];

  // Compute \frac{\partial v}{\partial x}
  grad_u[1][0] = 1.0-2.0*x[0]; // 0.0;

  // Compute \frac{\partial v}{\partial y}
  grad_u[1][1] = 0.0;

  return grad_u;
 }


 /// Pseudo_singularity for computation of Couette flow
 double pressure_pseudo_singularity_for_couette(const Vector<double>& x)
 {
  return -2.0*x[1]; //-2.0*x[0];
 }


 //----------------------------------------------------------------------


 /// Random center for blended couette pseudo singularity
 Vector<double> X_c_couette(2,0.0);

 /// Blended Pseudo_singularity for computation of Couette flow
 Vector<double> blended_velocity_pseudo_singularity_for_couette(const Vector<double>& x)
 {
  // Initialise velocity vector to return
  Vector<double> u=velocity_pseudo_singularity_for_couette(x);

  // Now blend
  double b=blend(x,X_c_couette);
  for (unsigned i=0;i<2;i++)
   {
    u[i]*=b;
   }
  return u;
 }


 /// Blended Pseudo_singularity for computation of Couette flow: 
 /// grad[i][j] = du_i/dx_j
 Vector<Vector<double> > blended_grad_velocity_pseudo_singularity_for_couette(
  const Vector<double>& x)
 {
  // unblended velocity
  Vector<double> u=velocity_pseudo_singularity_for_couette(x);
  Vector<Vector<double> > grad_u=
   grad_velocity_pseudo_singularity_for_couette(x);

  Vector<Vector<double> > blended_grad_u(2);
  for (unsigned d=0;d<2;d++)
   {
    blended_grad_u[d].resize(2);
   }

  // Now blend
  double b=blend(x,X_c_couette);
  Vector<double> grad_b=grad_blend(x,X_c_couette);
  for (unsigned i=0;i<2;i++)
   {
     for (unsigned j=0;j<2;j++)
      {
       blended_grad_u[i][j]=grad_u[i][j]*b+u[i]*grad_b[j];
      }
   }

  return blended_grad_u;
 }


 /// Blended Pseudo_singularity for computation of Couette flow
 double blended_pressure_pseudo_singularity_for_couette(const Vector<double>& x)
 {
  double p=pressure_pseudo_singularity_for_couette(x);

  // Now blend
  double b=blend(x,X_c_couette);

  return b*p;
 }






//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


 // Corner singularities at x=0 and x=1
 //------------------------------------

 /// Function that computes the fitting velocity solution near the corner (0,0)
 Vector<double> velocity_singularity1(const Vector<double>& x)
 {
  using namespace MathematicalConstants;

  // Initialise velocity vector to return
  Vector<double> u(2,0.0);
  
  // Find angle of polar coordinates centered at the point (0,0)
  double theta = atan2(x[1],x[0]);
  
  u[0] = -(0.5*Pi*theta*sin(theta) + theta*cos(theta) - 0.25*Pi*Pi*sin(theta))*sin(theta)/(-1.0 + 0.25*Pi*Pi) - 1.0*(-theta*sin(theta) + 0.5*Pi*theta*cos(theta) + 0.5*Pi*sin(theta) - 0.25*Pi*Pi*cos(theta) + cos(theta))*cos(theta)/(-1.0 + 0.25*Pi*Pi);

  u[1] = (0.5*Pi*theta*sin(theta) + theta*cos(theta) - 0.25*Pi*Pi*sin(theta))*cos(theta)/(-1.0 + 0.25*Pi*Pi) - 1.0*(-theta*sin(theta) + 0.5*Pi*theta*cos(theta) + 0.5*Pi*sin(theta) - 0.25*Pi*Pi*cos(theta) + cos(theta))*sin(theta)/(-1.0 + 0.25*Pi*Pi);

  // return the velocity vector
  return u;
 }

 /// Function that computes the fitting velocity solution near the corner (1,0)
 Vector<double> velocity_singularity2(const Vector<double>& x)
 {
  using namespace MathematicalConstants;

  // Variable change
  Vector<double> x1(2);
  x1[0] = 1.0-x[0];
  x1[1] = x[1];

  /// Polar coordinates centered at the point (1,0)
  //double r = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
  double theta = atan2(x1[1],x1[0]);

  // Initialise the velocity vector
  Vector<double> u(2,0.0);

  /// Value of the U component at the angle theta
  u[0] = (0.5*Pi*theta*sin(theta) + theta*cos(theta) - 0.25*Pi*Pi*sin(theta))*sin(theta)/(-1.0 + 0.25*Pi*Pi) + 1.0*(-theta*sin(theta) + 0.5*Pi*theta*cos(theta) + 0.5*Pi*sin(theta) - 0.25*Pi*Pi*cos(theta) + cos(theta))*cos(theta)/(-1.0 + 0.25*Pi*Pi);
  u[0] = - u[0];

  /// Value of the V component at the angle theta
  u[1] = -(0.5*Pi*theta*sin(theta) + theta*cos(theta) - 0.25*Pi*Pi*sin(theta))*cos(theta)/(-1.0 + 0.25*Pi*Pi) + 1.0*(-theta*sin(theta) + 0.5*Pi*theta*cos(theta) + 0.5*Pi*sin(theta) - 0.25*Pi*Pi*cos(theta) + cos(theta))*sin(theta)/(-1.0 + 0.25*Pi*Pi);

  return u;
 }

  
/// Function that computes the gradient of the fitting velocity solution near
/// the corner (0,0): grad[i][j] = du_i/dx_j
 Vector<Vector<double> > grad_velocity_singularity1(const Vector<double>& coord)
 {  
  using namespace MathematicalConstants;

  // Initialise the gradient matrix to return
  Vector<Vector<double> > grad_u(2);
  for (unsigned d=0;d<2;d++)
   {
    grad_u[d].resize(2);
    for (unsigned k=0;k<2;k++)
     {
      grad_u[d][k] = 0.0;
     }
   }

  // Cartesian coordinates
  double x = coord[0];
  double y = coord[1];
 
  // Value of \frac{\partial u}{\partial x}
  grad_u[0][0] = pow(x, 2)*(0.5*Pi*x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi, 2)*x/sqrt(pow(x, 2) + pow(y, 2)) + 1.0*x/sqrt(pow(x, 2) + pow(y, 2)) - 1.0*y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.5*Pi*y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi, 2))*pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L)) - x*y*(-x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.5*Pi*y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.25*pow(Pi, 2)*y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi, 2))*pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L)) - x*(-0.5*Pi*pow(x, 2)*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 1.0*pow(x, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.25*pow(Pi, 2)*pow(x, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 1.0*x*y*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 1.0*Pi*x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 1.0*pow(y, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.5*Pi*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi, 2)/sqrt(pow(x, 2) + pow(y, 2)) + 1.0/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(x, 2) + pow(y, 2))) + y*(pow(x, 2)*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.5*Pi*x*y*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 0.25*pow(Pi, 2)*x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.5*Pi*pow(y, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(x, 2) + pow(y, 2))) - (0.5*Pi*x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi, 2)*x/sqrt(pow(x, 2) + pow(y, 2)) + 1.0*x/sqrt(pow(x, 2) + pow(y, 2)) - 1.0*y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.5*Pi*y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(x, 2) + pow(y, 2)));

  // Value of \frac{\partial u}{\partial y}
  grad_u[0][1] = x*y*(0.5*Pi  *x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi  , 2)*x/sqrt(pow(x, 2) + pow(y, 2)) + 1.0*x/sqrt(pow(x, 2) + pow(y, 2)) - 1.0*y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.5*Pi  *y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L)) - x*(0.5*Pi  *pow(x, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 0.5*Pi  *x*y*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 2.0*x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.25*pow(Pi  , 2)*x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 1.0*pow(y, 2)*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 0.5*Pi  *pow(y, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 1.0*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.5*Pi  /sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*sqrt(pow(x, 2) + pow(y, 2))) - pow(y, 2)*(-x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.5*Pi  *y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.25*pow(Pi  , 2)*y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L)) + y*(-pow(x, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + x*y*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 0.5*Pi  *x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.5*Pi  *pow(y, 2)*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 0.25*pow(Pi  , 2)*pow(y, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 0.5*Pi  *atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.25*pow(Pi  , 2)/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*sqrt(pow(x, 2) + pow(y, 2))) + (-x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.5*Pi  *y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.25*pow(Pi  , 2)*y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*sqrt(pow(x, 2) + pow(y, 2)));

  // Value of \frac{\partial v}{\partial x}
  grad_u[1][0] = -pow(x, 2)*(x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.5*Pi  *y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi  , 2)*y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L)) + x*y*(0.5*Pi  *x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi  , 2)*x/sqrt(pow(x, 2) + pow(y, 2)) + 1.0*x/sqrt(pow(x, 2) + pow(y, 2)) - 1.0*y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.5*Pi  *y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L)) + x*(-pow(x, 2)*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 0.5*Pi  *x*y*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.25*pow(Pi  , 2)*x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 0.5*Pi  *pow(y, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*sqrt(pow(x, 2) + pow(y, 2))) - y*(-0.5*Pi  *pow(x, 2)*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 1.0*pow(x, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.25*pow(Pi  , 2)*pow(x, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 1.0*x*y*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 1.0*Pi  *x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 1.0*pow(y, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.5*Pi  *atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi  , 2)/sqrt(pow(x, 2) + pow(y, 2)) + 1.0/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*sqrt(pow(x, 2) + pow(y, 2))) + (x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.5*Pi  *y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi  , 2)*y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*sqrt(pow(x, 2) + pow(y, 2)));

  // Value of \frac{\partial v}{\partial y}
  grad_u[1][1] = -x*y*(x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.5*Pi  *y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi  , 2)*y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L)) + x*(pow(x, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - x*y*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.5*Pi  *x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 0.5*Pi  *pow(y, 2)*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.25*pow(Pi  , 2)*pow(y, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.5*Pi  *atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi  , 2)/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*sqrt(pow(x, 2) + pow(y, 2))) + pow(y, 2)*(0.5*Pi  *x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi  , 2)*x/sqrt(pow(x, 2) + pow(y, 2)) + 1.0*x/sqrt(pow(x, 2) + pow(y, 2)) - 1.0*y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.5*Pi  *y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L)) - y*(0.5*Pi  *pow(x, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 0.5*Pi  *x*y*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 2.0*x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 0.25*pow(Pi  , 2)*x*y/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) + 1.0*pow(y, 2)*atan2(y, x)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 0.5*Pi  *pow(y, 2)/pow(pow(x, 2) + pow(y, 2), 3.0L/2.0L) - 1.0*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.5*Pi  /sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*sqrt(pow(x, 2) + pow(y, 2))) - (0.5*Pi  *x*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) - 0.25*pow(Pi  , 2)*x/sqrt(pow(x, 2) + pow(y, 2)) + 1.0*x/sqrt(pow(x, 2) + pow(y, 2)) - 1.0*y*atan2(y, x)/sqrt(pow(x, 2) + pow(y, 2)) + 0.5*Pi  *y/sqrt(pow(x, 2) + pow(y, 2)))/((-1 + 0.25*pow(Pi  , 2))*sqrt(pow(x, 2) + pow(y, 2)));
 
  return grad_u;
 
 } // End of function


 /// Function that computes the gradient of the fitting velocity solution near
 /// the corner (1,0): grad[i][j] = du_i/dx_j
 Vector<Vector<double> > grad_velocity_singularity2(const Vector<double>& coord)
 { 
  using namespace MathematicalConstants;

  // Cartesian coordinates
  double x = coord[0];
  double y = coord[1];
  
  // Initialise the gradient matrix to return
  Vector<Vector<double> > grad_u(2);
  for (unsigned d=0;d<2;d++)
   {
    grad_u[d].resize(2);
    for (unsigned k=0;k<2;k++)
     {
      grad_u[d][k] = 0.0;
     }
   }

  // Value of \frac{\partial u}{\partial x}
  grad_u[0][0] = -y*(-x + 1)*(-0.5*Pi*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.25*pow(Pi, 2)*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) - (-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L)) - y*(-0.5*Pi*pow(y, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.5*Pi*y*(-x + 1)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 0.25*pow(Pi, 2)*y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - pow(-x + 1, 2)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2))) + pow(-x + 1, 2)*(-1.0*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*(-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 0.25*pow(Pi, 2)*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 1.0*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L)) + (-x + 1)*(-1.0*pow(y, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 1.0*y*(-x + 1)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 1.0*Pi*y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 0.5*Pi*pow(-x + 1, 2)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.25*pow(Pi, 2)*pow(-x + 1, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 1.0*pow(-x + 1, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.5*Pi*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 1.0/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.25*pow(Pi, 2)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2))) - (-1.0*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*(-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 0.25*pow(Pi, 2)*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 1.0*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2)));

  // Value of \frac{\partial u}{\partial y}
  grad_u[0][1] = pow(y, 2)*(-0.5*Pi*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.25*pow(Pi, 2)*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) - (-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L)) - y*(-x + 1)*(-1.0*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*(-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 0.25*pow(Pi, 2)*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 1.0*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L)) - y*(0.5*Pi*pow(y, 2)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.25*pow(Pi, 2)*pow(y, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + y*(-x + 1)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.5*Pi*y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - pow(-x + 1, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.5*Pi*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.25*pow(Pi, 2)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2))) + (-x + 1)*(1.0*pow(y, 2)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.5*Pi*pow(y, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.5*Pi*y*(-x + 1)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 2.0*y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 0.25*pow(Pi, 2)*y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 0.5*Pi*pow(-x + 1, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 1.0*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2))) - (-0.5*Pi*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.25*pow(Pi, 2)*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) - (-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2)));

  // Value of \frac{\partial v}{\partial x}
  grad_u[1][0] = -y*(-x + 1)*(-1.0*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*(-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 0.25*pow(Pi, 2)*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 1.0*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L)) - y*(-1.0*pow(y, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 1.0*y*(-x + 1)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 1.0*Pi*y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 0.5*Pi*pow(-x + 1, 2)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.25*pow(Pi, 2)*pow(-x + 1, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 1.0*pow(-x + 1, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.5*Pi*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 1.0/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.25*pow(Pi, 2)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2))) + pow(-x + 1, 2)*(0.5*Pi*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 0.25*pow(Pi, 2)*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) + (-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L)) + (-x + 1)*(0.5*Pi*pow(y, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 0.5*Pi*y*(-x + 1)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.25*pow(Pi, 2)*y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + pow(-x + 1, 2)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2))) - (0.5*Pi*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 0.25*pow(Pi, 2)*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) + (-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2)));

  // Value of \frac{\partial v}{\partial y}
  grad_u[1][1] = pow(y, 2)*(-1.0*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*(-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 0.25*pow(Pi, 2)*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 1.0*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L)) - y*(-x + 1)*(0.5*Pi*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 0.25*pow(Pi, 2)*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) + (-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L)) - y*(1.0*pow(y, 2)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.5*Pi*pow(y, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 0.5*Pi*y*(-x + 1)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 2.0*y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 0.25*pow(Pi, 2)*y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 0.5*Pi*pow(-x + 1, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - 1.0*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2))) + (-x + 1)*(-0.5*Pi*pow(y, 2)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 0.25*pow(Pi, 2)*pow(y, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) - y*(-x + 1)*atan2(y, -x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 0.5*Pi*y*(-x + 1)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + pow(-x + 1, 2)/pow(pow(y, 2) + pow(-x + 1, 2), 3.0L/2.0L) + 0.5*Pi*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 0.25*pow(Pi, 2)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2))) - (-1.0*y*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*y/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 0.5*Pi*(-x + 1)*atan2(y, -x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) - 0.25*pow(Pi, 2)*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)) + 1.0*(-x + 1)/sqrt(pow(y, 2) + pow(-x + 1, 2)))/((-1 + 0.25*pow(Pi, 2))*sqrt(pow(y, 2) + pow(-x + 1, 2)));


  grad_u[0][0]*=-1.0;
  grad_u[0][1]*=-1.0;
  grad_u[1][0]*=-1.0;
  grad_u[1][1]*=-1.0;

  return grad_u;
 }
 
 /// Function that computes the fitting pressure solution near the corner (0,0)
 double pressure_singularity1(const Vector<double>& x)
 {  
  using namespace MathematicalConstants;

  // Polar coordinates centered at the point (0,0)
  double r = sqrt(x[0]*x[0]+x[1]*x[1]);
  double theta = atan2(x[1],x[0]);
  
  double p = -(-0.25*Pi*Pi*Pi*sin(theta) + 1.0*Pi*sin(theta) - 0.5*Pi*Pi*cos(theta) + 2.0*cos(theta))/(r*(-0.0625*Pi*Pi*Pi*Pi - 1.0 + 0.5*Pi*Pi));
  
  //hierher
  if (std::isinf(p))
   {
    p=200.0;
   }
  
  return p;
 }

 /// Function that computes the fitting pressure solution near the corner (1,0)
 double pressure_singularity2(const Vector<double>& x)
 {
  using namespace MathematicalConstants;

  // Variable change
  Vector<double> x1(2);
  x1[0] = 1.0-x[0];
  x1[1] = x[1];

  /// Polar coordinates centered at the point (1,0)
  double r = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
  double theta = atan2(x1[1],x1[0]);
  
  double p = -(-1.0*Pi*sin(theta) + 0.25*Pi*Pi*Pi*sin(theta) - 2.0*cos(theta) + 0.5*Pi*Pi*cos(theta))/(r*(-0.0625*Pi*Pi*Pi*Pi - 1.0 + 0.5*Pi*Pi));

  // hierher 
  if (std::isinf(p))
   {
    p=200.0;
   }

  return p;
 }

 /// Function that computes the gradient of the fitting pressure solution near
 /// the corner (0,0)
 Vector<double> grad_pressure_singularity1(const Vector<double>& coord)
 {
  using namespace MathematicalConstants;

  // Cartesian coordinates
  double x = coord[0];
  double y = coord[1];
  
  // Initialise the gradient vector
  Vector<double> grad_p(2,0.0);

  // Value of \frac{\partial p}{\partial x}
  grad_p[0] = -x*(-0.5*Pi*Pi*x/sqrt(x*x + y*y) + 2.0*x/sqrt(x*x + y*y)
                  - 0.25*Pi*Pi*Pi*y/sqrt(x*x + y*y) + 1.0*Pi*y/sqrt(x*x + y*y))
   /((x*x + y*y)*sqrt(x*x + y*y)*(-0.0625*Pi*Pi*Pi*Pi - 1.0 + 0.5*Pi*Pi))
   + (-2.0*x*x/(x*x + y*y)*sqrt(x*x + y*y) + 0.5*Pi*Pi*x*x/(x*x + y*y)
      *sqrt(x*x + y*y) - 1.0*Pi*x*y/(x*x + y*y)*sqrt(x*x + y*y) + 0.25*Pi
      *Pi*Pi*x*y/(x*x + y*y)*sqrt(x*x + y*y) - 0.5*Pi*Pi/sqrt(x*x + y*y)
      + 2.0/sqrt(x*x + y*y))
   /(sqrt(x*x + y*y)*(-0.0625*Pi*Pi*Pi*Pi - 1.0 + 0.5*Pi*Pi));

  // Value of \frac{\partial p}{\partial y}
  grad_p[1] = -y*(-0.5*Pi*Pi*x/sqrt(x*x + y*y) + 2.0*x/sqrt(x*x + y*y)
                  - 0.25*Pi*Pi*Pi*y/sqrt(x*x + y*y) + 1.0*Pi*y/sqrt(x*x + y*y))
   /((x*x + y*y)*sqrt(x*x + y*y)*(-0.0625*Pi*Pi*Pi*Pi - 1.0 + 0.5*Pi*Pi))
   + (-2.0*x*y/(x*x + y*y)*sqrt(x*x + y*y) + 0.5*Pi*Pi*x*y/(x*x + y*y)
      *sqrt(x*x + y*y) - 1.0*Pi*y*y/(x*x + y*y)*sqrt(x*x + y*y) + 0.25*Pi
      *Pi*Pi*y*y/(x*x + y*y)*sqrt(x*x + y*y) - 0.25*Pi*Pi*Pi
      /sqrt(x*x + y*y) + 1.0*Pi/sqrt(x*x + y*y))
   /(sqrt(x*x + y*y)*(-0.0625*Pi*Pi*Pi*Pi - 1.0 + 0.5*Pi*Pi));

  return grad_p;
  
 } // End of function
   
   

 
 /// Function that computes the fitting solutions near the corner (0,0)
 void analytic_solution1(const Vector<double>& x, Vector<double>& solution)
 {
  using namespace MathematicalConstants;

  /// Polar coordinates centered at the point (0,0)
  double r = sqrt(x[0]*x[0] + x[1]*x[1]);
  double theta = atan2(x[1],x[0]);

  /// Value of the U component at the angle theta
  solution[0] = -(0.5*Pi*theta*sin(theta) + theta*cos(theta) - 0.25*Pi*Pi*sin(theta))*sin(theta)/(-1.0 + 0.25*Pi*Pi) - 1.0*(-theta*sin(theta) + 0.5*Pi*theta*cos(theta) + 0.5*Pi*sin(theta) - 0.25*Pi*Pi*cos(theta) + cos(theta))*cos(theta)/(-1.0 + 0.25*Pi*Pi);

  /// Value of the V component at the angle theta
  solution[1] = (0.5*Pi*theta*sin(theta) + theta*cos(theta) - 0.25*Pi*Pi*sin(theta))*cos(theta)/(-1.0 + 0.25*Pi*Pi) - 1.0*(-theta*sin(theta) + 0.5*Pi*theta*cos(theta) + 0.5*Pi*sin(theta) - 0.25*Pi*Pi*cos(theta) + cos(theta))*sin(theta)/(-1.0 + 0.25*Pi*Pi);

  /// Value of P at the angle theta
  if (x[0]==0.0 && x[1]==0.0)
   {
    solution[2] = -200.0;
   }
  else
   {
    solution[2] = -(-0.25*Pi*Pi*Pi*sin(theta) + 1.0*Pi*sin(theta) - 0.5*Pi*Pi*cos(theta) + 2.0*cos(theta))/(r*(-0.0625*Pi*Pi*Pi*Pi - 1.0 + 0.5*Pi*Pi));
   }
 }

 /// Function that computes the fitting solutions near the corner (1,0)
 void analytic_solution2(const Vector<double>& x, Vector<double>& solution)
 {
  using namespace MathematicalConstants;

  // Variable change
  Vector<double> x1(2);
  x1[0] = 1.0-x[0];
  x1[1] = x[1];

  /// Polar coordinates centered at the point (1,0)
  double r = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
  double theta = atan2(x1[1],x1[0]);

  /// Value of the U component at the angle theta
  solution[0] = (0.5*Pi*theta*sin(theta) + theta*cos(theta) - 0.25*Pi*Pi*sin(theta))*sin(theta)/(-1.0 + 0.25*Pi*Pi) + 1.0*(-theta*sin(theta) + 0.5*Pi*theta*cos(theta) + 0.5*Pi*sin(theta) - 0.25*Pi*Pi*cos(theta) + cos(theta))*cos(theta)/(-1.0 + 0.25*Pi*Pi);
  solution[0] = - solution[0];

  /// Value of the V component at the angle theta
  solution[1] = -(0.5*Pi*theta*sin(theta) + theta*cos(theta) - 0.25*Pi*Pi*sin(theta))*cos(theta)/(-1.0 + 0.25*Pi*Pi) + 1.0*(-theta*sin(theta) + 0.5*Pi*theta*cos(theta) + 0.5*Pi*sin(theta) - 0.25*Pi*Pi*cos(theta) + cos(theta))*sin(theta)/(-1.0 + 0.25*Pi*Pi);

  /// Value of P at the angle theta
  if (x1[0]==0.0 && x1[1]==0.0)
   {
    solution[2] = 200.0;
   }
  else
   {
    solution[2] = -(-1.0*Pi*sin(theta) + 0.25*Pi*Pi*Pi*sin(theta) - 2.0*cos(theta) + 0.5*Pi*Pi*cos(theta))/(r*(-0.0625*Pi*Pi*Pi*Pi - 1.0 + 0.5*Pi*Pi));
   }
 }




//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


 // Blended corner singularities
 //-----------------------------



 // Vector<double> blended_velocity_singularity1(const Vector<double>& x)
 // Vector<double> blended_velocity_singularity2(const Vector<double>& x)
 // Vector<Vector<double> > blended_grad_velocity_singularity1(const Vector<double>& coord)
 // Vector<Vector<double> > blended_grad_velocity_singularity2(const Vector<double>& coord)
 // double blended_pressure_singularity1(const Vector<double>& x)
 // double blended_pressure_singularity2(const Vector<double>& x)


 // Vector<double> grad_pressure_singularity1(const Vector<double>& coord)



 
} // end of namespace



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain
//====================================================================
template<class ELEMENT>
class RectangularDrivenCavityProblem : public Problem
{

public:


 /// Constructor
 RectangularDrivenCavityProblem();

 /// Destructor (empty)
 ~RectangularDrivenCavityProblem(){}

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure


 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the before solve (empty)
 void actions_before_newton_solve(){}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 // What it says...
 void actions_before_newton_step()
 {
  // CRDoubleMatrix jac;
  // DoubleVector res;
  // get_jacobian(res,jac);
  // jac.sparse_indexed_output("jac.dat");
  
  // unsigned n_dof=ndof();
  // DenseDoubleMatrix fd_jac(n_dof);
  // get_fd_jacobian(res,fd_jac);
  // fd_jac.sparse_indexed_output("fd_jac.dat");
  // pause("done");
 }

 /// Build line plot points
 void setup_line_plot_points()
  {
   // Get geom object
   Mesh_as_geom_object_pt=new MeshAsGeomObject(Navier_stokes_mesh_pt);

   // Number of radial lines
   unsigned n_phi=5;
   Line_plot_point.resize(n_phi);

   // How many points to you want?
   unsigned npt=100;
   double r_max=0.1;
   Vector<double> zeta(2);   
   Vector<double> s(2);
   GeomObject* geom_object_pt=0;
   for (unsigned i=0;i<n_phi;i++)
    {
     Line_plot_point[i].resize(npt);
     double phi=0.5*MathematicalConstants::Pi*double(i)/double(n_phi-1);
     oomph_info << "setup at phi : " << phi/(0.5*MathematicalConstants::Pi) 
                << " pi/2" << std::endl;
     for (unsigned j=0;j<npt;j++)
      {
       double r=r_max*double(j+1)/double(npt);
       zeta[0]=r*cos(phi);
       zeta[1]=r*sin(phi);
       Mesh_as_geom_object_pt->locate_zeta(zeta,geom_object_pt,s);
        if (geom_object_pt==0)
         {
          oomph_info << "Point : " 
                     << zeta[0] << " " 
                     << zeta[1] << " " 
                     << " not found in setup of line plots" 
                     << std::endl;
         }   
        Line_plot_point[i][j]=std::make_pair(geom_object_pt,s);
      }
    }
  }

private:
 
 /// Pointer to the element defining the unknown C1 (Note: eqn element is
 /// templated by the wrapped element!)
 SingularNavierStokesSolutionElement<ELEMENT>* Singular_navier_stokes_solution1_element_pt;

 /// Pointer to the element defining the unknown C2 (Note: eqn element is
 /// templated by the wrapped element!)
 SingularNavierStokesSolutionElement<ELEMENT>* Singular_navier_stokes_solution2_element_pt;
 
 /// Pointer to the C mesh associated with the elements defining the
 /// unknowns C1 and C2
 Mesh* Singular_navier_stokes_solution_mesh_pt;

 /// Pointer to mesh with Poisson elements
 Mesh* Navier_stokes_mesh_pt;

 /// Mesh as geom object representation of mesh
 MeshAsGeomObject* Mesh_as_geom_object_pt;

 /// Line_plot_point[i_phi][i_rho]
 Vector<Vector<std::pair<GeomObject*,Vector<double> > > > Line_plot_point;

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for RectangularDrivenCavity problem 
//========================================================================
template<class ELEMENT>
RectangularDrivenCavityProblem<ELEMENT>::RectangularDrivenCavityProblem()
{ 

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=Global_Parameters::N_x * Global_Parameters::Element_multiplier; 

 // # of elements in y-direction
 unsigned n_y=Global_Parameters::N_y * Global_Parameters::Element_multiplier; 

 // Domain length in x-direction
 double l_x=Global_Parameters::L_x;

 // Domain length in y-direction
 double l_y=Global_Parameters::L_y;

 // Build and assign mesh
 Navier_stokes_mesh_pt =
  new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);


 // Now we have the mesh, let's set up the plot points
 setup_line_plot_points();

 // Add Navier-Stokes mesh to global mesh
 add_sub_mesh(Navier_stokes_mesh_pt);

#ifdef DONT_USE_SINGULARITY

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = Navier_stokes_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= Navier_stokes_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       Navier_stokes_mesh_pt->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries


#ifdef COUETTE

 // Couette
 //--------
 {
  // Set the boundary conditions for this problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here. 
  Vector<double> x(2);
  Vector<double> u_couette(2);
  unsigned num_bound = Navier_stokes_mesh_pt->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
   {
    unsigned num_nod= Navier_stokes_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      // Get coordinates
      for (unsigned i=0;i<2;i++)
       {
        x[i]=Navier_stokes_mesh_pt->boundary_node_pt(ibound,inod)->x(i); 
       }
      
      // Get Couette velocity
      u_couette=Global_Parameters::velocity_couette(x);
      
      // Set veloc
      for (unsigned i=0;i<2;i++)
       {
        Navier_stokes_mesh_pt->boundary_node_pt(ibound,inod)->
         set_value(i,u_couette[i]); 
       }
     }
   } // end loop over boundaries
 }

#else

  // Setup tangential flow along boundary 0:
  unsigned ibound=0; 
  unsigned num_nod= Navier_stokes_mesh_pt->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Tangential flow
    unsigned i=0;
    Navier_stokes_mesh_pt->boundary_node_pt(ibound,inod)->set_value(i,1.0);
    // No penetration
    i=1;
    Navier_stokes_mesh_pt->boundary_node_pt(ibound,inod)->set_value(i,0.0);
   }

#endif




#else

 // Find number of nodes on one side of an element
 // Upcast from GeneralsedElement to the present element
 ELEMENT *spec_el_pt =
  dynamic_cast<ELEMENT*>(Navier_stokes_mesh_pt->element_pt(0));
 unsigned Nnode_1d = sqrt(spec_el_pt->nnode());


 // Boundary conditions
 //---------------------
 

 // LEFT BOUNDARY
 //----------------

 // Loop over the boundary elements
 for (unsigned i=0;i<n_y;i++)
  {
   // Number of the element
   double e = i*n_x;

   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Navier_stokes_mesh_pt->element_pt(e));

   // Loop over the boundary nodes
   for (unsigned k=0;k<Nnode_1d;k++)
    {

     // Local node number
     double j = k*Nnode_1d;

     // The x velocity component at node j is subject to Dirichlet BC
     // Its imposed value is 0.0 which is already the value by default
     // in the vector Imposed_velocity_values_at_node. So there is no
     // need to call the function set_velocity_dirichlet_value_on_node
     el_pt->impose_velocity_dirichlet_bc_on_node(j,0);
     
     // The y velocity component at node j is subject to Dirichlet BC
     // Its imposed value is 0.0 which is already the value by default
     // in the vector Imposed_velocity_values_at_node. So there is no
     // need to call the function set_velocity_dirichlet_value_on_node
     el_pt->impose_velocity_dirichlet_bc_on_node(j,1);


#ifdef COUETTE

     // Get coordinates
     Vector<double> x(2);
     for (unsigned i=0;i<2;i++)
      {
       x[i]=el_pt->node_pt(j)->x(i); 
      }
     
     // Get Couette velocity
     Vector<double> u_couette=Global_Parameters::velocity_couette(x);
     
     // Impose it
     el_pt->set_velocity_dirichlet_value_on_node(j,0,u_couette[0]);
     el_pt->set_velocity_dirichlet_value_on_node(j,1,u_couette[1]);
     
#endif



    }
  }

 
 // RIGHT BOUNDARY
 //----------------

 // Loop over the boundary elements
 for (unsigned i=0;i<n_y;i++)
  {
   // Number of the element
   double e = i*n_x + n_x-1;

   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Navier_stokes_mesh_pt->element_pt(e));

   // Loop over the boundary nodes
   for (unsigned k=0;k<Nnode_1d;k++)
    {
     // Local node number
     double j = k*Nnode_1d + Nnode_1d-1;

     // The x velocity component at node j is subject to Dirichlet BC
     // Its imposed value is 0.0 which is already the value by default
     // in the vector Imposed_velocity_values_at_node. So there is no
     // need to call the function set_velocity_dirichlet_value_on_node
     el_pt->impose_velocity_dirichlet_bc_on_node(j,0);
     
     // The y velocity component at node j is subject to Dirichlet BC
     // Its imposed value is 0.0 which is already the value by default
     // in the vector Imposed_velocity_values_at_node. So there is no
     // need to call the function set_velocity_dirichlet_value_on_node
     el_pt->impose_velocity_dirichlet_bc_on_node(j,1);


#ifdef COUETTE

     // Get coordinates
     Vector<double> x(2);
     for (unsigned i=0;i<2;i++)
      {
       x[i]=el_pt->node_pt(j)->x(i); 
      }
     
     // Get Couette velocity
     Vector<double> u_couette=Global_Parameters::velocity_couette(x);
     
     // Impose it
     el_pt->set_velocity_dirichlet_value_on_node(j,0,u_couette[0]);
     el_pt->set_velocity_dirichlet_value_on_node(j,1,u_couette[1]);
     
#endif


    }
  }


 // UPPER BOUNDARY
 //----------------

 // Loop over the boundary elements
 for (unsigned i=0;i<n_x;i++)
  {
   // Number of the element
   double e = (n_y-1)*n_x + i;

   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Navier_stokes_mesh_pt->element_pt(e));

   // Loop over the boundary nodes
   for (unsigned k=0;k<Nnode_1d;k++)
    {
     // Local node number
     double j = (Nnode_1d-1)*Nnode_1d + k;

     // The x velocity component at node j is subject to Dirichlet BC
     // Its imposed value is 0.0 which is already the value by default
     // in the vector Imposed_velocity_values_at_node. So there is no
     // need to call the function set_velocity_dirichlet_value_on_node
     el_pt->impose_velocity_dirichlet_bc_on_node(j,0);
     
     // The y velocity component at node j is subject to Dirichlet BC
     // Its imposed value is 0.0 which is already the value by default
     // in the vector Imposed_velocity_values_at_node. So there is no
     // need to call the function set_velocity_dirichlet_value_on_node
     el_pt->impose_velocity_dirichlet_bc_on_node(j,1);

#ifdef COUETTE

     // Get coordinates
     Vector<double> x(2);
     for (unsigned i=0;i<2;i++)
      {
       x[i]=el_pt->node_pt(j)->x(i); 
      }
     
     // Get Couette velocity
     Vector<double> u_couette=Global_Parameters::velocity_couette(x);
     
     // Impose it
     el_pt->set_velocity_dirichlet_value_on_node(j,0,u_couette[0]);
     el_pt->set_velocity_dirichlet_value_on_node(j,1,u_couette[1]);
     
     
#endif


    }
  }

 
 // BOTTOM BOUNDARY
 //----------------
 
 // Loop over the boundary elements
 for (unsigned e=0;e<n_x;e++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Navier_stokes_mesh_pt->element_pt(e));

   // Loop over the boundary nodes
   for (unsigned j=0;j<Nnode_1d;j++)
    {
     // The x velocity component at node j is subject to Dirichlet BC
     el_pt->impose_velocity_dirichlet_bc_on_node(j,0);
     
     // The y velocity component at node j is subject to Dirichlet BC
     // Its imposed value is 0.0 which is already the value by default
     // in the vector Imposed_velocity_values_at_node. So there is no
     // need to call the function set_velocity_dirichlet_value_on_node
     el_pt->impose_velocity_dirichlet_bc_on_node(j,1);



#ifdef COUETTE

     // Get coordinates
     Vector<double> x(2);
     for (unsigned i=0;i<2;i++)
      {
       x[i]=el_pt->node_pt(j)->x(i); 
      }
     
     // Get Couette velocity
     Vector<double> u_couette=Global_Parameters::velocity_couette(x);
     
     // Impose it
     el_pt->set_velocity_dirichlet_value_on_node(j,0,u_couette[0]);
     el_pt->set_velocity_dirichlet_value_on_node(j,1,u_couette[1]);
     
#else

     // // To break the symetry
     // double x = el_pt->nodal_position(j,0);
     // el_pt->set_velocity_dirichlet_value_on_node(j,0,1.0+sin(2.0*Global_Parameters::Pi*x)/2.0);

     // Its imposed value is 1.0
     el_pt->set_velocity_dirichlet_value_on_node(j,0,1.0);

     // // linear singularity + poiseuille
     // double y = el_pt->nodal_position(j,1);
     // el_pt->set_velocity_dirichlet_value_on_node(j,0,y*(1-y)+y);
     
#endif




    }
  }
 

 // PRESSURE AT (0.0,0.0)
 //----------------------
 #ifdef COUETTE

 // Upcast from GeneralsedElement to the present element
 ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Navier_stokes_mesh_pt->element_pt(0));

 // The pressure at (0.0,0.0) is subject to Dirichlet BC
 // Its imposed value is 0.0 which is already the value by default
 // in the vector Imposed_value_at_pressure_dof. So there is no
 // need to call the function set_dirichlet_value_on_pressure_dof
 el_pt->impose_dirichlet_bc_on_pressure_dof(0); 

 Vector<double> x_fixed_pressure(2,0.0);
 el_pt->set_dirichlet_value_on_pressure_dof
  (0,Global_Parameters::pressure_couette(x_fixed_pressure)); 

#else

 // Upcast from GeneralsedElement to the present element
 unsigned nel=Navier_stokes_mesh_pt->nelement();
 ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Navier_stokes_mesh_pt->element_pt(nel-1));

 // The pressure at (0.0,0.0) is subject to Dirichlet BC
 // Its imposed value is 0.0 which is already the value by default
 // in the vector Imposed_value_at_pressure_dof. So there is no
 // need to call the function set_dirichlet_value_on_pressure_dof
 el_pt->impose_dirichlet_bc_on_pressure_dof(3); 
 el_pt->set_dirichlet_value_on_pressure_dof(3,0.0);

#endif



 // Implementation of the singularities
 //------------------------------------

 // Create a C mesh 
 Singular_navier_stokes_solution_mesh_pt=new Mesh;


#ifdef COUETTE

 // Create a pointer to a C1EquationElement object
 Singular_navier_stokes_solution1_element_pt = new SingularNavierStokesSolutionElement<ELEMENT>;

 
 if (CommandLineArgs::command_line_flag_has_been_set("--blend"))
  {
   // Set the pointer to the first velocity singularity function for this
   // element, defined in Global_parameters namespace
   Singular_navier_stokes_solution1_element_pt->velocity_singular_fct_pt() =
    &Global_Parameters::blended_velocity_pseudo_singularity_for_couette;
   
   // Set the pointer to the gradient of the first velocity singularity function
   // for this element, defined in Global_parameters namespace
   Singular_navier_stokes_solution1_element_pt->grad_velocity_singular_fct_pt() =
    &Global_Parameters::blended_grad_velocity_pseudo_singularity_for_couette;
   
   // Set the pointer to the first pressure singularity function for this 
   // element, defined in Global_parameters namespace
   Singular_navier_stokes_solution1_element_pt->pressure_singular_fct_pt() =
    &Global_Parameters::blended_pressure_pseudo_singularity_for_couette;
  
   // The singular function satisfies the Stokes equation
   Singular_navier_stokes_solution1_element_pt->singular_function_satisfies_stokes_equation()=false;
  }
 else
  {
   // Set the pointer to the first velocity singularity function for this 
   // element, defined in Global_parameters namespace
   Singular_navier_stokes_solution1_element_pt->velocity_singular_fct_pt() =
    &Global_Parameters::velocity_pseudo_singularity_for_couette;
   
   // Set the pointer to the gradient of the first velocity singularity function
   // for this element, defined in Global_parameters namespace
   Singular_navier_stokes_solution1_element_pt->grad_velocity_singular_fct_pt() =
    &Global_Parameters::grad_velocity_pseudo_singularity_for_couette;
   
   // Set the pointer to the first pressure singularity function for this 
   // element, defined in Global_parameters namespace
   Singular_navier_stokes_solution1_element_pt->pressure_singular_fct_pt() =
    &Global_Parameters::pressure_pseudo_singularity_for_couette;
  
   // The singular function satisfies the Stokes equation
   Singular_navier_stokes_solution1_element_pt->singular_function_satisfies_stokes_equation()=true;
  }

 // Tell the C1EquationElement object about its associated Navier-Stokes
 // element, the component of the velocity whose derivative will be computed
 // in the residual, the local coordinate in the Navier_Stokes element at which
 // the residual will be computed, and the direction of the derivative
 Vector<double> s1_pt(2);
 s1_pt[0] = -1.0;
 s1_pt[1] = -1.0;
 ELEMENT *el1_pt = dynamic_cast<ELEMENT*>
  (Navier_stokes_mesh_pt->element_pt(0));
 Singular_navier_stokes_solution1_element_pt->
  set_wrapped_navier_stokes_element_pt
  (el1_pt,s1_pt,&Global_Parameters::Direction);

 // Add element to the C mesh
 Singular_navier_stokes_solution_mesh_pt->
  add_element_pt(Singular_navier_stokes_solution1_element_pt);

#else

 // Create a pointer to a C1EquationElement object
 Singular_navier_stokes_solution1_element_pt = new SingularNavierStokesSolutionElement<ELEMENT>;

 // Set the pointer to the first velocity singularity function for this element,
 // defined in Global_parameters namespace
 Singular_navier_stokes_solution1_element_pt->velocity_singular_fct_pt() =
  &Global_Parameters::velocity_singularity1;

 // Set the pointer to the gradient of the first velocity singularity function
 // for this element, defined in Global_parameters namespace
 Singular_navier_stokes_solution1_element_pt->grad_velocity_singular_fct_pt() =
  &Global_Parameters::grad_velocity_singularity1;

 // Set the pointer to the first pressure singularity function for this element,
 // defined in Global_parameters namespace
 Singular_navier_stokes_solution1_element_pt->pressure_singular_fct_pt() =
  &Global_Parameters::pressure_singularity1;

 // The singular function satisfies the Stokes equation
 Singular_navier_stokes_solution1_element_pt->singular_function_satisfies_stokes_equation()=true;

 // Tell the C1EquationElement object about its associated Navier-Stokes
 // element, the component of the velocity whose derivative will be computed
 // in the residual, the local coordinate in the Navier_Stokes element at which
 // the residual will be computed, and the direction of the derivative
 Vector<double> s1_pt(2);
 s1_pt[0] = -1.0;
 s1_pt[1] = -1.0;
 ELEMENT *el1_pt = dynamic_cast<ELEMENT*>
  (Navier_stokes_mesh_pt->element_pt(0));
 Singular_navier_stokes_solution1_element_pt->set_wrapped_navier_stokes_element_pt
  (el1_pt,s1_pt,&Global_Parameters::Direction);

 // Add element to the C mesh
 Singular_navier_stokes_solution_mesh_pt->
  add_element_pt(Singular_navier_stokes_solution1_element_pt);

 // Pin value of C
 Singular_navier_stokes_solution1_element_pt->pin_c();
 Singular_navier_stokes_solution1_element_pt->set_c(1.0);


 // Create a pointer to a C2EquationElement object
 Singular_navier_stokes_solution2_element_pt = new SingularNavierStokesSolutionElement<ELEMENT>;
 
 // Set the pointer to the second velocity singularity function for
 // this element, defined in Global_parameters namespace
 Singular_navier_stokes_solution2_element_pt->velocity_singular_fct_pt() =
  &Global_Parameters::velocity_singularity2;

 // Set the pointer to the gradient of the second velocity singularity function
 // for this element, defined in Global_parameters namespace
 Singular_navier_stokes_solution2_element_pt->grad_velocity_singular_fct_pt() =
  &Global_Parameters::grad_velocity_singularity2;

 // Set the pointer to the second pressure singularity function for 
 // this element, defined in Global_parameters namespace
 Singular_navier_stokes_solution2_element_pt->pressure_singular_fct_pt() =
  &Global_Parameters::pressure_singularity2;

 // The singular function satisfies the Stokes equation
 Singular_navier_stokes_solution2_element_pt->singular_function_satisfies_stokes_equation()=true;

 // Tell the C2EquationElement object about its associated Navier-Stokes
 // element, the component of the velocity whose derivative will be computed
 // in the residual, the local coordinate in the Navier_Stokes element at which
 // the residual will be computed, and the direction of the derivative
 Vector<double> s2_pt(2);
 s2_pt[0] = 1.0;
 s2_pt[1] = -1.0;
 ELEMENT *el2_pt = dynamic_cast<ELEMENT*>
  (Navier_stokes_mesh_pt->element_pt(n_x-1));
 Singular_navier_stokes_solution2_element_pt->set_wrapped_navier_stokes_element_pt
  (el2_pt,s2_pt,&Global_Parameters::Direction);
 
 // Add element to the C mesh
 Singular_navier_stokes_solution_mesh_pt->
  add_element_pt(Singular_navier_stokes_solution2_element_pt);

 // Pin value of C
 Singular_navier_stokes_solution2_element_pt->pin_c();
 Singular_navier_stokes_solution2_element_pt->set_c(1.0);

#endif

 // Add C mesh to global mesh
 add_sub_mesh(Singular_navier_stokes_solution_mesh_pt);

#endif

 // Build global mesh
 build_global_mesh();

 // Make better integration scheme
 const unsigned NPTS_1D=20;   
 GaussLegendre<2,NPTS_1D>* gauss_legendre_pt=                            
  new GaussLegendre<2,NPTS_1D>;
 
 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = Navier_stokes_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Navier_stokes_mesh_pt->
                                           element_pt(e));

   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Assert that we're on a fixed mesh
   el_pt->disable_ALE();

   // Set integration scheme?
   if (CommandLineArgs::command_line_flag_has_been_set
       ("--use_high_order_integration_scheme"))
    {
     el_pt->set_integration_scheme(gauss_legendre_pt);
    }

#ifndef DONT_USE_SINGULARITY

   // Set the pointer to the element that determines the amplitude
   // of the first singular fct
   el_pt->add_c_equation_element_pt(Singular_navier_stokes_solution1_element_pt);

#ifndef COUETTE

   // Set the pointer to the element that determines the amplitude
   // of the second singular fct
   el_pt->add_c_equation_element_pt(Singular_navier_stokes_solution2_element_pt);

#endif

#endif
   
  } // end loop over elements

#ifdef COUETTE

#ifdef DONT_USE_SINGULARITY
 
 // Now set the first FE pressure value in first element 
 Vector<double> x_fixed_pressure(2,0.0);
 fix_pressure(0,0,Global_Parameters::pressure_couette(x_fixed_pressure));

#endif
 
#else

#ifdef DONT_USE_SINGULARITY

 // Now set the last FE pressure value in last element
 unsigned nell=Navier_stokes_mesh_pt->nelement();
 fix_pressure(nell-1,3,0.0);

#endif

#endif
 

#ifdef OOMPH_HAS_MUMPS

 // Change solver
 linear_solver_pt()=new MumpsSolver;

 // Shut up
 dynamic_cast<MumpsSolver*>(linear_solver_pt())->
  enable_suppress_warning_about_MPI_COMM_WORLD();

#endif

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RectangularDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file,some_file2;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=10; // 100; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Navier_stokes_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output coarse solution 
 sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned npts_coarse=2;
 Navier_stokes_mesh_pt->output(some_file,npts_coarse);
 some_file.close();


 // Output perturbation solution 
 //-----------------------------
 sprintf(filename,"%s/perturbation_soln_two_term%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 PerturbationSolution::N_terms_for_plot=2;
 some_file.open(filename);
 Navier_stokes_mesh_pt->output_fct(
  some_file,npts,
  PerturbationSolution::perturbation_soln_for_plot);
 some_file.close();

 // Output perturbation solution 
 //-----------------------------
 sprintf(filename,"%s/perturbation_soln_one_term%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 PerturbationSolution::N_terms_for_plot=1;
 Navier_stokes_mesh_pt->output_fct(
  some_file,npts,
  PerturbationSolution::perturbation_soln_for_plot);
 some_file.close();


 // Output first-order perturbation solution (Stokes)
 //--------------------------------------------------
 sprintf(filename,"%s/first_order_perturbation%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Navier_stokes_mesh_pt->output_fct(
  some_file,npts,
  PerturbationSolution::first_order_perturbation_for_plot);
 some_file.close();



 // Output second-order perturbation solution (Stokes)
 //--------------------------------------------------
 sprintf(filename,"%s/second_order_perturbation%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Navier_stokes_mesh_pt->output_fct(
  some_file,npts,
  PerturbationSolution::second_order_perturbation_for_plot);
 some_file.close();



 // Do line plots
 sprintf(filename,"%s/line_plot%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 sprintf(filename,"%s/second_order_perturbation_line_plot%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file2.open(filename);
 Vector<double> s(2); 
 Vector<double> x(2); 
 Vector<double> u(2);
 double p=0.0;
 Vector<double> soln(3);
 unsigned nphi=Line_plot_point.size();
 for (unsigned i=0;i<nphi;i++)
  {
   unsigned nr=Line_plot_point[i].size();
   some_file << "ZONE T=\"phi = " << i << "/" << nphi-1 << " pi/2 \"" 
             << std::endl;  
   some_file2 << "ZONE T=\"phi = " << i << "/" << nphi << " pi/2 \"" 
             << std::endl;  
   for (unsigned j=0;j<nr;j++)
    {
     ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Line_plot_point[i][j].first);
     s=Line_plot_point[i][j].second;
     el_pt->interpolated_x(s,x);
#ifdef DONT_USE_SINGULARITY
     u[0]=el_pt->interpolated_u_nst(s,0);
     u[1]=el_pt->interpolated_u_nst(s,1);
#else
     u=el_pt->interpolated_u_nst(s);
#endif
     p=el_pt->interpolated_p_nst(s);
     some_file << x[0] << " " 
               << x[1] << " " 
               << sqrt(x[0]*x[0]+x[1]*x[1]) << " " 
               << u[0] << " " 
               << u[1] << " " 
               << p << " " 
               << std::endl;

     PerturbationSolution::second_order_perturbation_for_plot(x,soln);
     some_file2 << x[0] << " " 
                << x[1] << " " 
                << sqrt(x[0]*x[0]+x[1]*x[1]) << " " 
                << soln[0] << " " 
                << soln[1] << " " 
                << soln[2] << " " 
                << std::endl;
    }
  }
 some_file.close();
 some_file2.close();


#ifndef DONT_USE_SINGULARITY

 oomph_info << "Value of C1: " 
            << Singular_navier_stokes_solution1_element_pt->c() << std::endl;

#ifndef COUETTE
 oomph_info << "Value of C2: " 
            << Singular_navier_stokes_solution2_element_pt->c() << std::endl;
#endif

#endif


#ifdef COUETTE

 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_couette_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Navier_stokes_mesh_pt->output_fct(some_file,npts,
                                   Global_Parameters::exact_couette);
 some_file.close();

// #ifdef DONT_USE_SINGULARITY

//  // Doc error
//  //----------
//  double error,norm;
//  sprintf(filename,"%s/couette_error%i.dat",
//          doc_info.directory().c_str(),
//          doc_info.number());
//  some_file.open(filename);
//  Navier_stokes_mesh_pt->compute_error(some_file,
//                                       Global_Parameters::exact_couette,
//                                       error,norm);   
//  some_file.close();
//  cout << "error (without pressure): " << sqrt(error) << std::endl; 
//  cout << "norm  (without pressure): " << sqrt(norm) << std::endl << std::endl;

// #else

 // Doc error
 //----------
 double error=0.0;
 double norm=0.0;
 sprintf(filename,"%s/couette_error_with_pressure%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 bool include_pressure=true;
 unsigned nel=Navier_stokes_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   double el_error=0.0;
   double el_norm=0.0;
   dynamic_cast<ELEMENT*>(Navier_stokes_mesh_pt->element_pt(e))->
    compute_error(some_file,
                  Global_Parameters::exact_couette,
                  include_pressure,el_error,el_norm); 
   norm+=el_norm;
   error+=el_error;
  }
 some_file.close();
 cout << "error (with pressure): " << sqrt(error) << std::endl; 
 cout << "norm  (with pressure): " << sqrt(norm) << std::endl << std::endl;



 error=0.0;
 norm=0.0;
 sprintf(filename,"%s/couette_error_with_pressure%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 include_pressure=false;
 for (unsigned e=0;e<nel;e++)
  {
   double el_error=0.0;
   double el_norm=0.0;
   dynamic_cast<ELEMENT*>(Navier_stokes_mesh_pt->element_pt(e))->
    compute_error(some_file,
                  Global_Parameters::exact_couette,
                  include_pressure,el_error,el_norm); 
   norm+=el_norm;
   error+=el_error;
  }
 some_file.close();
 cout << "error (without pressure): " << sqrt(error) << std::endl; 
 cout << "norm  (without pressure): " << sqrt(norm) << std::endl << std::endl;

// #endif

#endif

} // end_of_doc_solution





/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////







//==start_of_main======================================================
/// Driver for RectangularDrivenCavity test problem -- test drive
/// with two different types of element.
//=====================================================================
int main(int argc, char **argv)
{

#ifdef OOMPH_HAS_MPI

 MPI_Helpers::init(argc,argv);

#endif


 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Use blending?
 CommandLineArgs::specify_command_line_flag("--blend");

 // Multiplier for elements
 CommandLineArgs::specify_command_line_flag(
  "--element_multiplier",
  &Global_Parameters::Element_multiplier);

 // Multiplier for elements
 CommandLineArgs::specify_command_line_flag(
  "--re",
  &Global_Physical_Variables::Re);

 // Use higher order integration scheme?
 CommandLineArgs::specify_command_line_flag(
  "--use_high_order_integration_scheme");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();


 // Check singularities and their gradients
 {
  ofstream some_file;
  some_file.open("u_sing.dat");
  unsigned nplot=100;
  some_file << "ZONE I=" << nplot << ", J=" << nplot << std::endl;
  Vector<double> x(2);
  for (unsigned i=0;i<nplot;i++)
   {
    x[0]=double(i)/double(nplot-1);
    for (unsigned j=0;j<nplot;j++)
     {
      x[1]=double(j)/double(nplot-1);
      
      // Get veloc
      //Vector<double> u=Global_Parameters::velocity_singularity1(x);
      //Vector<Vector<double> > grad_u=Global_Parameters::grad_velocity_singularity1(x);

      //Vector<double> u=Global_Parameters::velocity_singularity2(x);
      //Vector<Vector<double> > grad_u=Global_Parameters::grad_velocity_singularity2(x);

      Vector<double> u=Global_Parameters::blended_velocity_pseudo_singularity_for_couette(x);
      Vector<Vector<double> > grad_u=Global_Parameters::blended_grad_velocity_pseudo_singularity_for_couette(x);

      some_file << x[0] << " " << x[1] << " " 
                << u[0] << " " << u[1] << " " 
                << grad_u[0][0] << " " << grad_u[0][1] << " "
                << grad_u[1][0] << " " << grad_u[1][1] << " "
                << std::endl;
     }
   }
  some_file.close();
 }

 
 // Set up doc info
 // ---------------

 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 

#ifdef DONT_USE_SINGULARITY

 // Build the problem with QTaylorHoodElements
 RectangularDrivenCavityProblem<MyTaylorHoodElement> problem;

#else

 // Build the problem with singular QTaylorHoodElements
 RectangularDrivenCavityProblem<NavierStokesElementWithSingularity<
  MyTaylorHoodElement> > problem;

#endif

 // Solve the problem
 problem.newton_solve();
 
 // Outpt the solution
 problem.doc_solution(doc_info);
 
 // Step number
 doc_info.number()++;
 

} // end_of_main










