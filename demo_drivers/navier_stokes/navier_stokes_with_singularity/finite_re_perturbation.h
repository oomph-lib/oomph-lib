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

//====================================================================
/// Christian's perturbation solution
//====================================================================
namespace PerturbationSolution
{

  using namespace MathematicalConstants;
  using namespace Global_Physical_Variables;

  /// function to compute radial distance of this point
  double get_radial_distance(const Vector<double>& x)
  {
    return sqrt(x[0]*x[0] + x[1]*x[1]);
  }

  /// function to return angle of radial line from positive x-axis 
  double get_azimuthal_angle(const Vector<double>& x)
  {
    return atan2(x[1], x[0]);
  }

  /// function to convert velocity vectors from polar to Cartesian coords
  void polar_to_cartesian_velocity(const Vector<double> u_polar,
				   const double phi, 
				   Vector<double>& u_cart)
  {
    u_cart[0] = u_polar[0] * cos(phi) - u_polar[1] * sin(phi);
    u_cart[1] = u_polar[0] * sin(phi) + u_polar[1] * cos(phi);
  }
  
  /// Leading order solution (u and p are in polar coordinates)
  void leading_order_veloc_and_pressure(const double r, const double phi, 
					 Vector<double>& u_0, double& p_0)
  {
    // leading order velocity term, r component
    u_0[0] = 
      (((4 + 2*phi*Pi - pow(Pi,2))*cos(phi) + 2*(-2*phi + Pi)*sin(phi)))/(-4 + pow(Pi,2));

    // phi component
    u_0[1] = ((-4*phi*cos(phi) + Pi*(-2*phi + Pi)*sin(phi)))/(-4 + pow(Pi,2));
    
    // leading order pressure term
    p_0 = (4*(2*cos(phi) + Pi*sin(phi)))/((-4 + pow(Pi,2))*r);
  }

  /// First order correction (u and p are in polar coordinates)
  void first_order_veloc_and_pressure(const double r, const double phi, 
				      Vector<double>& u_1, double& p_1)
  {
    // first order velocity term, r component
    u_1[0] = (r*(-2*pow(Pi,2)*(8 + pow(Pi,2)) + 
       2*(8*pow(phi,2)*(-4 + pow(Pi,2)) - 8*phi*Pi*(2 + pow(Pi,2)) + pow(Pi,2)*(8 + pow(Pi,2)))*
        cos(2*phi) + (32*phi - 72*Pi - 64*pow(phi,2)*Pi + 24*phi*pow(Pi,2) + 6*pow(Pi,3) + pow(Pi,5))*
        sin(2*phi)))/(32.*pow(-4 + pow(Pi,2),2));

    // phi component
    u_1[1] = (r*(Pi*(24 - 14*pow(Pi,2) - pow(Pi,4) + 4*phi*Pi*(8 + pow(Pi,2))) + 
       (-64*pow(phi,2)*Pi + 8*phi*(12 + pow(Pi,2)) + Pi*(-24 + 14*pow(Pi,2) + pow(Pi,4)))*cos(2*phi) - 
       2*(24 + 10*pow(Pi,2) + pow(Pi,4) + 8*pow(phi,2)*(-4 + pow(Pi,2)) - 8*phi*Pi*(6 + pow(Pi,2)))*
        sin(2*phi)))/(32.*pow(-4 + pow(Pi,2),2));

    // first order pressure term
    p_1 = -(pow(Pi,2)*(-8 + pow(Pi,2))*log(r))/(4.*pow(-4 + pow(Pi,2),2));
  }

  /// Second order correction (u and p are in polar coordinates)
  void second_order_veloc_and_pressure(const double r, const double phi, 
				       Vector<double>& u_2, double& p_2)
  {
    // second order velocity term, r component
    u_2[0] = (pow(r,2)*(2*(-3952 + 3672*pow(Pi,2) - 51*pow(Pi,4) - 21*pow(Pi,6) + 
          192*pow(phi,3)*Pi*(4 + pow(Pi,2)) + 72*phi*Pi*(2 + pow(Pi,2))*(-8 + 3*pow(Pi,2)) + 
          144*pow(phi,2)*(20 + 5*pow(Pi,2) - 2*pow(Pi,4)))*cos(phi) + 
       2*(3952 - 3672*pow(Pi,2) + 51*pow(Pi,4) + 21*pow(Pi,6) + 288*pow(phi,3)*Pi*(-12 + pow(Pi,2)) - 
          432*pow(phi,2)*(-12 + 5*pow(Pi,2) + pow(Pi,4)) + 24*phi*Pi*(276 + 85*pow(Pi,2) + 9*pow(Pi,4)))*
        cos(3*phi) + (-768*pow(phi,3)*(4 + pow(Pi,2)) + 144*pow(phi,2)*Pi*(20 + 13*pow(Pi,2)) - 
          36*phi*(-32 + 32*pow(Pi,2) + 42*pow(Pi,4) + pow(Pi,6)) + 
          Pi*(-3648 + 2576*pow(Pi,2) + 564*pow(Pi,4) + 9*pow(Pi,6)))*sin(phi) + 
       (3072*Pi + 1152*pow(phi,3)*(4 - 3*pow(Pi,2)) + 432*pow(phi,2)*Pi*(36 + 5*pow(Pi,2)) - 
          pow(Pi,3)*(2560 + 252*pow(Pi,2) + 27*pow(Pi,4)) + 
          12*phi*(-1024 + 72*pow(Pi,2) + 54*pow(Pi,4) + 3*pow(Pi,6)))*sin(3*phi)))/
   (4608.*pow(-4 + pow(Pi,2),3));

    // phi component
    u_2[1] = (pow(r,2)*((-768*pow(phi,3)*(4 + pow(Pi,2)) + 144*pow(phi,2)*Pi*(-12 + 5*pow(Pi,2)) - 
          36*phi*(-224 - 16*pow(Pi,2) + 10*pow(Pi,4) + pow(Pi,6)) + 
          Pi*(2112 + 1424*pow(Pi,2) + 132*pow(Pi,4) + 9*pow(Pi,6)))*cos(phi) + 
       (-384*pow(phi,3)*(-4 + 3*pow(Pi,2)) + 48*pow(phi,2)*Pi*(156 + 11*pow(Pi,2)) - 
          Pi*(2112 + 1424*pow(Pi,2) + 132*pow(Pi,4) + 9*pow(Pi,6)) + 
          4*phi*(-1856 + 3*pow(Pi,2)*(6 + pow(Pi,2))*(28 + pow(Pi,2))))*cos(3*phi) + 
       2*(-80 - 3960*pow(Pi,2) + 231*pow(Pi,4) + 39*pow(Pi,6) - 192*pow(phi,3)*Pi*(4 + pow(Pi,2)) + 
          144*pow(phi,2)*(12 + 3*pow(Pi,2) + 2*pow(Pi,4)) - 72*phi*Pi*(-40 + 8*pow(Pi,2) + 3*pow(Pi,4)))*
        sin(phi) - 2*(80 - 888*pow(Pi,2) + 85*pow(Pi,4) + 9*pow(Pi,6) + 
          96*pow(phi,3)*Pi*(-12 + pow(Pi,2)) + 8*phi*Pi*(588 + 107*pow(Pi,2) + 9*pow(Pi,4)) - 
          48*pow(phi,2)*(-52 + 3*pow(Pi,2)*(9 + pow(Pi,2))))*sin(3*phi)))/(1536.*pow(-4 + pow(Pi,2),3));

    // second order pressure term
    p_2 = (r*(-((9632 - 6264*pow(Pi,2) + 552*pow(Pi,4) + 33*pow(Pi,6) + 144*phi*Pi*(28 - 13*pow(Pi,2)) + 
            192*pow(phi,3)*Pi*(4 + pow(Pi,2)) - 288*pow(phi,2)*(20 + pow(Pi,2)))*cos(phi)) + 
       9*(-352 + 80*pow(Pi,2) + 22*pow(Pi,4) + pow(Pi,6) + 32*pow(phi,2)*(4 - 3*pow(Pi,2)) + 
          16*phi*Pi*(36 + pow(Pi,2)))*cos(3*phi) + 
       (384*pow(phi,3)*(4 + pow(Pi,2)) - 144*pow(phi,2)*Pi*(-20 + 3*pow(Pi,2) + pow(Pi,4)) + 
          36*phi*(-256 - 72*pow(Pi,2) + 10*pow(Pi,4) + pow(Pi,6)) + 
          Pi*(7152 - 124*pow(Pi,2) + 24*pow(Pi,4) + 9*pow(Pi,6)))*sin(phi) - 
       36*(96*phi - 48*(-2 + pow(phi,2))*Pi - 56*phi*pow(Pi,2) + 4*(4 + pow(phi,2))*pow(Pi,3) - 
          4*phi*pow(Pi,4) + pow(Pi,5))*sin(3*phi)))/(576.*pow(-4 + pow(Pi,2),3));
  }

  /// array of function pointers for each term in expansion
  void (* veloc_and_pressure_term_pt [3])(const double,
					  const double,
					  Vector<double>&,
					  double&) =
  { & leading_order_veloc_and_pressure,
    & first_order_veloc_and_pressure,
    & second_order_veloc_and_pressure
  };
  


  /// Full solution. Final optional argument specifies 
  /// number of terms in expansion (Stokes or Stokes plus 
  /// leading order correction in Re)
  void full_soln(const Vector<double>& x,
		 Vector<double>& u, double& p, 
                 const unsigned& nterms=2)
  {
    // get radial distance to this point from origin
    double r = get_radial_distance(x);
    
    // get azimuthal angle of this point
    double phi = get_azimuthal_angle(x);
    
    // vector to hold each velocity term (Cartesian)
    Vector<double> u_i(2);

    // vector to hold velocity term (polar coordinates)
    Vector<double> u_polar(2);
    
    // variable to hold each pressure term
    double p_i = 0;

    // initialise solutions
    u[0] = 0.0;
    u[1] = 0.0;
    p    = 0.0;

    // -------------------------------------------
    // compute each term and add it to the solution
    for (unsigned i = 0; i < nterms; i++)
    {
      // compute the term (in polar coordinates)
      veloc_and_pressure_term_pt[i](r, phi, u_polar, p_i);

      // convert coordinates to Cartesians
      polar_to_cartesian_velocity(u_polar, phi, u_i);
	
      // multiply by appropriate power of Re and add it to the solution
      u[0] -= pow(Re, i) * u_i[0];
      u[1] -= pow(Re, i) * u_i[1];
      p    -= pow(Re, i) * p_i;      
    }
  }


  /// Number of terms in perturbation expansion for plot
  unsigned N_terms_for_plot=2;
  
  /// Combined solution (u,v,p) for plotting
  void perturbation_soln_for_plot(const Vector<double>& x, 
                                 Vector<double>& soln)
  {
   soln.resize(3);
   double p=0.0;
   full_soln(x,soln,p,N_terms_for_plot);
   if (std::isinf(p)) p=200.0;
   if (std::isnan(p)) p=200.0;
   soln[2]=p;
  }

  /// First order-perturbation for plotting
  void first_order_perturbation_for_plot(const Vector<double>& x, 
                                        Vector<double>& soln)
  { 
   soln.resize(3);

   // get radial distance to this point from origin
   double r = get_radial_distance(x);
   
   // get azimuthal angle of this point
   double phi = get_azimuthal_angle(x);
   
   // vector to hold each velocity term (Cartesian)
   Vector<double> u_i(2);
   
   // vector to hold velocity term (polar coordinates)
   Vector<double> u_polar(2);
   
   // variable to hold each pressure term
   double p_i = 0;
   
   // compute the term (in polar coordinates)
   veloc_and_pressure_term_pt[0](r, phi, u_polar, p_i);
   
   // convert coordinates to Cartesians
   polar_to_cartesian_velocity(u_polar, phi, u_i);
	
   // multiply by appropriate power of Re and add it to the solution
   soln[0]=-u_i[0];
   soln[1]=-u_i[1];
   soln[2]=-p_i;      
 
   for (unsigned i=0;i<3;i++)
    {
     if (std::isinf(soln[i])) soln[i]=200.0;
     if (std::isnan(soln[i])) soln[i]=200.0;
    }
   
  }

  /// Second order-perturbation for plotting
  void second_order_perturbation_for_plot(const Vector<double>& x, 
                                         Vector<double>& soln)
  {
   soln.resize(3);

   // get radial distance to this point from origin
   double r = get_radial_distance(x);
   
   // get azimuthal angle of this point
   double phi = get_azimuthal_angle(x);
   
   // vector to hold each velocity term (Cartesian)
   Vector<double> u_i(2);
   
   // vector to hold velocity term (polar coordinates)
   Vector<double> u_polar(2);
   
   // variable to hold each pressure term
   double p_i = 0;
   
   // compute the term (in polar coordinates)
   veloc_and_pressure_term_pt[1](r, phi, u_polar, p_i);
   
   // convert coordinates to Cartesians
   polar_to_cartesian_velocity(u_polar, phi, u_i);
	
   // multiply by appropriate power of Re and add it to the solution
   soln[0]=-u_i[0];
   soln[1]=-u_i[1];
   soln[2]=-p_i;      
 
   for (unsigned i=0;i<3;i++)
    {
     if (std::isinf(soln[i])) soln[i]=200.0;
     if (std::isnan(soln[i])) soln[i]=200.0;
    }
   
  }


  
} // end of namespace
