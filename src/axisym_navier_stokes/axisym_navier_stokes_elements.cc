//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
//Non-inline functions for NS elements
#include "axisym_navier_stokes_elements.h"


namespace oomph
{

/// Navier--Stokes equations static data
Vector<double> AxisymmetricNavierStokesEquations::Gamma(2,1.0);

//=================================================================
/// "Magic" negative number that indicates that the pressure is
/// not stored at a node. This cannot be -1 because that represents
/// the positional hanging scheme in the hanging_pt object of nodes
//=================================================================
 int AxisymmetricNavierStokesEquations::Pressure_not_stored_at_node = -100;

/// Navier--Stokes equations static data
double AxisymmetricNavierStokesEquations::
Default_Physical_Constant_Value = 0.0;

// Navier--Stokes equations static data
double AxisymmetricNavierStokesEquations::Default_Physical_Ratio_Value = 1.0;

/// Navier-Stokes equations default gravity vector
Vector<double> AxisymmetricNavierStokesEquations::
Default_Gravity_vector(3,0.0);

//======================================================================
/// Validate against exact velocity solution at given time.
/// Solution is provided via function pointer.
/// Plot at a given number of plot points and compute L2 error
/// and L2 norm of velocity solution over element.
//=======================================================================
void AxisymmetricNavierStokesEquations::
compute_error(std::ostream &outfile,
              FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
              const double& time,
              double& error, double& norm)
{
 error=0.0;
 norm=0.0;

 //Vector of local coordinates
 Vector<double> s(2);

 // Vector for coordintes
 Vector<double> x(2);

 //Set the value of Nintpt
 unsigned Nintpt = integral_pt()->nweight();
   
 outfile << "ZONE" << std::endl;

 // Exact solution Vector (u,v,w,p)
 Vector<double> exact_soln(4);
   
 //Loop over the integration points
 for(unsigned ipt=0;ipt<Nintpt;ipt++)
  {
   //Assign values of s
   for(unsigned i=0;i<2;i++) {s[i] = integral_pt()->knot(ipt,i);}

   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   // Get jacobian of mapping
   double J = J_eulerian(s);

   // Get x position as Vector
   interpolated_x(s,x);

   //Premultiply the weights and the Jacobian and r
   double W = w*J*x[0];

   // Get exact solution at this point
   (*exact_soln_pt)(time,x,exact_soln);

   // Velocity error
   for(unsigned i=0;i<3;i++)
    {
     norm+=exact_soln[i]*exact_soln[i]*W;
     error+=(exact_soln[i]-interpolated_u_axi_nst(s,i))*
      (exact_soln[i]-interpolated_u_axi_nst(s,i))*W;
    }

   //Output x,y,...,u_exact
   for(unsigned i=0;i<2;i++) {outfile << x[i] << " ";}

   //Output x,y,[z],u_error,v_error,[w_error]
   for(unsigned i=0;i<3;i++)
    {outfile << exact_soln[i]-interpolated_u_axi_nst(s,i) << " ";}

   outfile << std::endl;
  }
}

//======================================================================
/// Validate against exact velocity solution
/// Solution is provided via function pointer.
/// Plot at a given number of plot points and compute L2 error
/// and L2 norm of velocity solution over element.
//=======================================================================
void AxisymmetricNavierStokesEquations::
compute_error(std::ostream &outfile,
              FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
              double& error, double& norm)
{
 error=0.0;
 norm=0.0;
 
 //Vector of local coordinates
 Vector<double> s(2);
 
 // Vector for coordintes
 Vector<double> x(2);
 
 //Set the value of Nintpt
 unsigned Nintpt = integral_pt()->nweight();
 
 outfile << "ZONE" << std::endl;
 
 // Exact solution Vector (u,v,w,p)
 Vector<double> exact_soln(4);
 
 //Loop over the integration points
 for(unsigned ipt=0;ipt<Nintpt;ipt++)
  {
   
   //Assign values of s
   for(unsigned i=0;i<2;i++) {s[i] = integral_pt()->knot(ipt,i);}
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   // Get jacobian of mapping
   double J=J_eulerian(s);
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   //Premultiply the weights and the Jacobian and r
   double W = w*J*x[0];
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   // Velocity error
   for(unsigned i=0;i<3;i++)
    {
     norm+=exact_soln[i]*exact_soln[i]*W;
     error+=(exact_soln[i]-interpolated_u_axi_nst(s,i))*
      (exact_soln[i]-interpolated_u_axi_nst(s,i))*W;
    }

   //Output x,y,...,u_exact
   for(unsigned i=0;i<2;i++) {outfile << x[i] << " ";}

     //Output x,y,u_error,v_error,w_error
     for(unsigned i=0;i<3;i++)
      {outfile << exact_soln[i]-interpolated_u_axi_nst(s,i) << " ";}

     outfile << std::endl;
  }
}

//======================================================================
/// Output "exact" solution
/// Solution is provided via function pointer.
/// Plot at a given number of plot points.
/// Function prints as many components as are returned in solution Vector.
//=======================================================================
void AxisymmetricNavierStokesEquations::
output_fct(std::ostream &outfile, 
           const unsigned &nplot, 
           FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 
 //Vector of local coordinates
 Vector<double> s(2);
 
 // Vector for coordintes
 Vector<double> x(2);
 
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Exact solution Vector
 Vector<double> exact_soln;
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...
   for(unsigned i=0;i<2;i++)
    {
     outfile << x[i] << " ";
    }
   
   //Output "exact solution"
   for(unsigned i=0;i<exact_soln.size();i++)
    {
     outfile << exact_soln[i] << " ";
    }
   
   outfile << std::endl;
   
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);

}

//======================================================================
/// Output "exact" solution at a given time
/// Solution is provided via function pointer.
/// Plot at a given number of plot points.
/// Function prints as many components as are returned in solution Vector.
//=======================================================================
void AxisymmetricNavierStokesEquations::
output_fct(std::ostream &outfile,
           const unsigned &nplot, 
           const double& time,
           FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
{
 
 //Vector of local coordinates
 Vector<double> s(2);
 
 // Vector for coordintes
 Vector<double> x(2);
  
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Exact solution Vector
 Vector<double> exact_soln;
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(time,x,exact_soln);
   
   //Output x,y,...
   for(unsigned i=0;i<2;i++)
    {
     outfile << x[i] << " ";
    }
   
   //Output "exact solution"
   for(unsigned i=0;i<exact_soln.size();i++)
    {
     outfile << exact_soln[i] << " ";
    }
   
   outfile << std::endl;
   
  }
 
 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);

}

//==============================================================
/// Output function: Velocities only  
/// x,y,[z],u,v,[w]
/// in tecplot format at specified previous timestep (t=0: present;
/// t>0: previous timestep). Specified number of plot points in each
/// coordinate direction.
//==============================================================
void AxisymmetricNavierStokesEquations::
output_veloc(std::ostream &outfile, 
             const unsigned &nplot, 
             const unsigned &t) 
{
 //Find number of nodes
 unsigned n_node = nnode();
 
 //Local shape function
 Shape psi(n_node);

 //Vectors of local coordinates and coords and velocities
 Vector<double> s(2);
 Vector<double> interpolated_x(2);
 Vector<double> interpolated_u(3);


 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
  
   // Get shape functions
   shape(s,psi);

   // Loop over coordinate directions
   for(unsigned i=0;i<2;i++) 
    {
     interpolated_x[i]=0.0;
     //Loop over the local nodes and sum
     for(unsigned l=0;l<n_node;l++) 
      {interpolated_x[i] += nodal_position(t,l,i)*psi[l];}
    }
   
   //Loop over the velocity components
   for(unsigned i=0;i<3;i++) 
    {
     //Get the index at which the velocity is stored
     unsigned u_nodal_index = u_index_axi_nst(i);
     interpolated_u[i]=0.0;
     //Loop over the local nodes and sum
     for(unsigned l=0;l<n_node;l++) 
      {interpolated_u[i] += nodal_value(t,l,u_nodal_index)*psi[l];}
    }

   // Coordinates
   for(unsigned i=0;i<2;i++) {outfile << interpolated_x[i] << " ";}
   
   // Velocities
   for(unsigned i=0;i<3;i++) {outfile << interpolated_u[i] << " ";}
   
   outfile << std::endl;   
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);

}

//==============================================================
/// Output function: 
/// r,z,u,v,w,p
/// in tecplot format. Specified number of plot points in each
/// coordinate direction.
//==============================================================
void AxisymmetricNavierStokesEquations::output(std::ostream &outfile, 
                                               const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(2);
 
  // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
  
   // Coordinates
   for(unsigned i=0;i<2;i++) {outfile << interpolated_x(s,i) << " ";}
   
   // Velocities
   for(unsigned i=0;i<3;i++) {outfile << interpolated_u_axi_nst(s,i) << " ";}
   
   // Pressure
   outfile << interpolated_p_axi_nst(s)  << " ";
   
   outfile << std::endl;   
  }
 outfile << std::endl;

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);

}


//==============================================================
/// Output function: 
/// r,z,u,v,w,p
/// in tecplot format. Specified number of plot points in each
/// coordinate direction.
//==============================================================
void AxisymmetricNavierStokesEquations::output(FILE* file_pt,
                                               const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(2);
 
 // Tecplot header info
  fprintf(file_pt,"%s ",tecplot_zone_string(nplot).c_str());

 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
  
   // Coordinates
   for(unsigned i=0;i<2;i++) 
    {
     //outfile << interpolated_x(s,i) << " ";
     fprintf(file_pt,"%g ",interpolated_x(s,i));
    }
   
   // Velocities
   for(unsigned i=0;i<3;i++) 
    {
     //outfile << interpolated_u(s,i) << " ";
     fprintf(file_pt,"%g ",interpolated_u_axi_nst(s,i));
    }
   
   // Pressure
   //outfile << interpolated_p(s)  << " ";
   fprintf(file_pt,"%g ",interpolated_p_axi_nst(s));

   //outfile << std::endl;   
   fprintf(file_pt,"\n");
  }
 //outfile << std::endl;
 fprintf(file_pt,"\n");

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(file_pt,nplot);

}




//==============================================================
/// Return integral of dissipation over element
//==============================================================
double AxisymmetricNavierStokesEquations::dissipation() const
{  
 throw OomphLibError(
  "Check the dissipation calculation for axisymmetric NSt",
  "AxisymmetricNavierStokesEquations::dissipation()",
  OOMPH_EXCEPTION_LOCATION);
 
 // Initialise
 double diss=0.0;

 //Set the value of Nintpt
 unsigned Nintpt = integral_pt()->nweight();
 
 //Set the Vector to hold local coordinates
 Vector<double> s(2);
 
 //Loop over the integration points
 for(unsigned ipt=0;ipt<Nintpt;ipt++)
  {
   
   //Assign values of s
   for(unsigned i=0;i<2;i++)
    {
     s[i] = integral_pt()->knot(ipt,i);
    }
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   // Get Jacobian of mapping
   double J = J_eulerian(s);
   
   // Get strain rate matrix
   DenseMatrix<double> strainrate(3,3);
   strain_rate(s,strainrate);

   // Initialise
   double local_diss=0.0;
   for(unsigned i=0;i<3;i++) 
    {
     for(unsigned j=0;j<3;j++) 
      {    
       local_diss+=2.0*strainrate(i,j)*strainrate(i,j);
      }
    }
   
   diss+=local_diss*w*J;
  }

 return diss;

}

//==============================================================
/// \short Compute traction (on the viscous scale) at local
/// coordinate s for outer unit normal N
//==============================================================
void AxisymmetricNavierStokesEquations::traction(const Vector<double>& s,
                                   const Vector<double>& N, 
                                   Vector<double>& traction)
{
 throw OomphLibError(
  "Check the traction calculation for axisymmetric NSt",
  "AxisymmetricNavierStokesEquations::traction()",
  OOMPH_EXCEPTION_LOCATION);

 // Get velocity gradients
 DenseMatrix<double> strainrate(3,3);
 strain_rate(s,strainrate);
 
 // Get pressure
 double press=interpolated_p_axi_nst(s);
 
 // Loop over traction components
 for (unsigned i=0;i<3;i++)
  {
   traction[i]=-press*N[i];
   for (unsigned j=0;j<3;j++)
    {
     traction[i]+=2.0*strainrate(i,j)*N[j];
    }
  }
}

//==============================================================
/// Return dissipation at local coordinate s
//==============================================================
double AxisymmetricNavierStokesEquations::
dissipation(const Vector<double>& s) const
{  
 throw OomphLibError(
  "Check the dissipation calculation for axisymmetric NSt",
  "AxisymmetricNavierStokesEquations::dissipation()",
  OOMPH_EXCEPTION_LOCATION);

 // Get strain rate matrix
 DenseMatrix<double> strainrate(3,3);
 strain_rate(s,strainrate);
 
 // Initialise
 double local_diss=0.0;
 for(unsigned i=0;i<3;i++) 
  {
   for(unsigned j=0;j<3;j++) 
    {    
     local_diss+=2.0*strainrate(i,j)*strainrate(i,j);
    }
  }
 
 return local_diss;
}

//==============================================================
/// Get strain-rate tensor: \f$ e_{ij} \f$  where 
/// \f$ i,j = r,z,\theta \f$ (in that order)
//==============================================================
void AxisymmetricNavierStokesEquations::
strain_rate(const Vector<double>& s, DenseMatrix<double>& strainrate) const
{
 
#ifdef PARANOID
 if ((strainrate.ncol()!=3)||(strainrate.nrow()!=3))
  {
   std::ostringstream error_message;
   error_message  << "The strain rate has incorrect dimensions " 
                  << strainrate.ncol() << " x " 
                  << strainrate.nrow() << " Not 3" << std::endl;
   
   throw OomphLibError(error_message.str(),
                       "AxisymmetricNavierStokeEquations::strain_rate()",
                         OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 //Find out how many nodes there are in the element
 unsigned n_node = nnode();
 
 //Set up memory for the shape and test functions
 Shape psi(n_node);
 DShape dpsidx(n_node,2);
 
 //Call the derivatives of the shape functions
 dshape_eulerian(s,psi,dpsidx);
 
 // Radius
 double interpolated_r = 0.0;

 // Velocity components and their derivatives
 double ur=0.0;
 double durdr=0.0;
 double durdz=0.0;
 double uz=0.0;
 double duzdr=0.0;
 double duzdz=0.0;
 double uphi=0.0;
 double duphidr=0.0;
 double duphidz=0.0;

 //Get the local storage for the indices
 unsigned u_nodal_index[3];
 for(unsigned i=0;i<3;++i) {u_nodal_index[i] = u_index_axi_nst(i);}

 // Loop over nodes to assemble velocities and their derivatives
 // w.r.t. to r and z (x_0 and x_1)
 for(unsigned l=0;l<n_node;l++) 
  {   
   interpolated_r += nodal_position(l,0)*psi[l];

   ur += nodal_value(l,u_nodal_index[0])*psi[l];
   uz += nodal_value(l,u_nodal_index[1])*psi[l];
   uphi += nodal_value(l,u_nodal_index[2])*psi[l];

   durdr += nodal_value(l,u_nodal_index[0])*dpsidx(l,0);
   durdz += nodal_value(l,u_nodal_index[0])*dpsidx(l,1);

   duzdr += nodal_value(l,u_nodal_index[1])*dpsidx(l,0);
   duzdz += nodal_value(l,u_nodal_index[1])*dpsidx(l,1);

   duphidr += nodal_value(l,u_nodal_index[2])*dpsidx(l,0);
   duphidz += nodal_value(l,u_nodal_index[2])*dpsidx(l,1);
  }

 
 // Assign strain rates without negative powers of the radius
 // and zero those with:
 strainrate(0,0)=durdr;
 strainrate(0,1)=0.5*(durdz+duzdr);
 strainrate(1,0)=strainrate(0,1);
 strainrate(0,2)=0.0;
 strainrate(2,0)=strainrate(0,2);
 strainrate(1,1)=duzdz;
 strainrate(1,2)=0.5*duphidz;
 strainrate(2,1)=strainrate(1,2);
 strainrate(2,2)=0.0;
 
 
 // Overwrite the strain rates with negative powers of the radius
 // unless we're at the origin
 if (std::abs(interpolated_r)>1.0e-16)
  {
   double inverse_radius=1.0/interpolated_r;
   strainrate(0,2)=0.5*(duphidr-inverse_radius*uphi);
   strainrate(2,0)=strainrate(0,2);
   strainrate(2,2)=inverse_radius*ur;
  }
 

}


//==============================================================
///  \short Get integral of kinetic energy over element:
//==============================================================
double AxisymmetricNavierStokesEquations::kin_energy() const
{  

 throw OomphLibError(
  "Check the kinetic energy calculation for axisymmetric NSt",
  "AxisymmetricNavierStokesEquations::kin_energy()",
  OOMPH_EXCEPTION_LOCATION);

 // Initialise
 double kin_en=0.0;

 //Set the value of Nintpt
 unsigned Nintpt = integral_pt()->nweight();
 
 //Set the Vector to hold local coordinates
 Vector<double> s(2);
 
 //Loop over the integration points
 for(unsigned ipt=0;ipt<Nintpt;ipt++)
  {   
   //Assign values of s
   for(unsigned i=0;i<2;i++) {s[i] = integral_pt()->knot(ipt,i);}
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   //Get Jacobian of mapping
   double J = J_eulerian(s);
   
   // Loop over directions
   double veloc_squared=0.0;
   for(unsigned i=0;i<3;i++) 
    {
     veloc_squared+=interpolated_u_axi_nst(s,i)*interpolated_u_axi_nst(s,i);
    }
   
   kin_en+=0.5*veloc_squared*w*J*interpolated_x(s,0);
  }

 return kin_en;

}

//==============================================================
/// Return pressure integrated over the element
//==============================================================
double AxisymmetricNavierStokesEquations::pressure_integral() const
{

 // Initialise 
 double press_int=0;

 //Set the value of Nintpt
 unsigned Nintpt = integral_pt()->nweight();
 
 //Set the Vector to hold local coordinates
 Vector<double> s(2);
 
 //Loop over the integration points
 for(unsigned ipt=0;ipt<Nintpt;ipt++)
  {
   
   //Assign values of s
   for(unsigned i=0;i<2;i++)
    {
     s[i] = integral_pt()->knot(ipt,i);
    }
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   //Get Jacobian of mapping
   double J = J_eulerian(s);
   
   //Premultiply the weights and the Jacobian
   double W = w*J*interpolated_x(s,0);
   
   // Get pressure
   double press=interpolated_p_axi_nst(s);
   
   // Add
   press_int+=press*W;
   
  }
 
 return press_int;

}

//==============================================================
///  Compute the residuals for the Navier--Stokes 
///  equations; flag=1(or 0): do (or don't) compute the 
///  Jacobian as well. 
//==============================================================
void AxisymmetricNavierStokesEquations::
fill_in_generic_residual_contribution_axi_nst(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              DenseMatrix<double> &mass_matrix,
                                              unsigned flag)
{
 //Find out how many nodes there are
 unsigned n_node = nnode();
 
 //Find out how many pressure dofs there are
 unsigned n_pres = npres_axi_nst();
 
 //Get the nodal indices at which the velocity is stored
 unsigned u_nodal_index[3];
 for(unsigned i=0;i<3;++i) {u_nodal_index[i] = u_index_axi_nst(i);}

 //Set up memory for the shape and test functions
 //Note that there are only two dimensions, r and z in this problem
 Shape psif(n_node), testf(n_node);
 DShape dpsifdx(n_node,2), dtestfdx(n_node,2);
   
 //Set up memory for pressure shape and test functions
 Shape psip(n_pres), testp(n_pres);

 //Number of integration points
 unsigned Nintpt = integral_pt()->nweight();
   
 //Set the Vector to hold local coordinates (two dimensions)
 Vector<double> s(2);

 //Get Physical Variables from Element
 //Reynolds number must be multiplied by the density ratio
 double scaled_re = re()*density_ratio();
 double scaled_re_st = re_st()*density_ratio();
 double scaled_re_inv_fr = re_invfr()*density_ratio();
 double scaled_re_inv_ro = re_invro()*density_ratio();
 double visc_ratio = viscosity_ratio();
 double dens_ratio = density_ratio();
 Vector<double> G = g();

 //Integers used to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<Nintpt;ipt++)
  {
   //Assign values of s
   for(unsigned i=0;i<2;i++) s[i] = integral_pt()->knot(ipt,i);
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   //Call the derivatives of the shape and test functions
   double J = 
    dshape_and_dtest_eulerian_at_knot_axi_nst(ipt,psif,dpsifdx,testf,dtestfdx);
   
   //Call the pressure shape and test functions
   pshape_axi_nst(s,psip,testp);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Allocate storage for the position and the derivative of the 
   //mesh positions wrt time
   Vector<double> interpolated_x(2,0.0);
   Vector<double> mesh_velocity(2,0.0);
   //Allocate storage for the pressure, velocity components and their
   //time and spatial derivatives
   double interpolated_p=0.0;
   Vector<double> interpolated_u(3,0.0);
   Vector<double> dudt(3,0.0);
   DenseMatrix<double> interpolated_dudx(3,2,0.0);
   
   //Calculate pressure at integration point
   for(unsigned l=0;l<n_pres;l++) {interpolated_p += p_axi_nst(l)*psip[l];}
   
   //Calculate velocities and derivatives at integration point

   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Cache the shape function
     const double psif_ = psif(l);
     //Loop over the two coordinate directions
     for(unsigned i=0;i<2;i++)
      {
       interpolated_x[i] += this->raw_nodal_position(l,i)*psif_;
      }
       //mesh_velocity[i]  += dnodal_position_dt(l,i)*psif[l];

     //Loop over the three velocity directions
     for(unsigned i=0;i<3;i++)
      {
       //Get the u_value
       const double u_value = this->raw_nodal_value(l,u_nodal_index[i]);
       interpolated_u[i] += u_value*psif_;
       dudt[i]+= du_dt_axi_nst(l,i)*psif_;
       //Loop over derivative directions
       for(unsigned j=0;j<2;j++)
        {interpolated_dudx(i,j) += u_value*dpsifdx(l,j);}
      }
    }

   //Get the mesh velocity if ALE is enabled
   if(!ALE_is_disabled)
    {
     // Loop over nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       //Loop over the two coordinate directions
       for(unsigned i=0;i<2;i++)
        {
         mesh_velocity[i]  += this->raw_dnodal_position_dt(l,i)*psif(l);
        }
      }
    }
       
   
   //Get the user-defined body force terms
   Vector<double> body_force(3);
   get_body_force(time(),interpolated_x,body_force);
   
   //Get the user-defined source function
   double source = get_source_fct(time(),interpolated_x);

   //r is the first position component
   double r = interpolated_x[0];


   //MOMENTUM EQUATIONS
   //------------------
   
   //Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {

     //FIRST (RADIAL) MOMENTUM EQUATION
     local_eqn = nodal_local_eqn(l,u_nodal_index[0]);
     //If it's not a boundary condition
     if(local_eqn >= 0)
      {
       //Add the user-defined body force terms
       residuals[local_eqn] += 
        r*body_force[0]*testf[l]*W;

       //Add the gravitational body force term
       residuals[local_eqn] += r*scaled_re_inv_fr*testf[l]*G[0]*W;
       
       //Add the pressure gradient term
       residuals[local_eqn]  += 
        interpolated_p*(testf[l] + r*dtestfdx(l,0))*W;
       
       //Add in the stress tensor terms
       //The viscosity ratio needs to go in here to ensure
       //continuity of normal stress is satisfied even in flows
       //with zero pressure gradient!
       residuals[local_eqn] -= visc_ratio*
        r*(1.0+Gamma[0])*interpolated_dudx(0,0)*dtestfdx(l,0)*W;
       
       residuals[local_eqn] -= visc_ratio*r*
        (interpolated_dudx(0,1) + Gamma[0]*interpolated_dudx(1,0))*
        dtestfdx(l,1)*W;

       residuals[local_eqn] -= 
        visc_ratio*(1.0 + Gamma[0])*interpolated_u[0]*testf[l]*W/r;
       
       //Add in the inertial terms
       //du/dt term
       residuals[local_eqn] -= scaled_re_st*r*dudt[0]*testf[l]*W;

       //Convective terms
       residuals[local_eqn] -= 
        scaled_re*(r*interpolated_u[0]*interpolated_dudx(0,0) 
            - interpolated_u[2]*interpolated_u[2] 
            + r*interpolated_u[1]*interpolated_dudx(0,1))*testf[l]*W;

       //Mesh velocity terms
       if(!ALE_is_disabled)
        {
         for(unsigned k=0;k<2;k++)
          {
           residuals[local_eqn] += 
            scaled_re_st*r*mesh_velocity[k]*interpolated_dudx(0,k)*testf[l]*W;
          }
        }

       //Add the Coriolis term
       residuals[local_eqn] += 
        2.0*r*scaled_re_inv_ro*interpolated_u[2]*testf[l]*W;

       //CALCULATE THE JACOBIAN
       if(flag)
        {
         //Loop over the velocity shape functions again
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           local_unknown = nodal_local_eqn(l2,u_nodal_index[0]);
           //Radial velocity component
           if(local_unknown >= 0)
            {
             if(flag==2)
              {
               //Add the mass matrix
               mass_matrix(local_eqn,local_unknown) +=
                scaled_re_st*r*psif[l2]*testf[l]*W;
              }

             //Add contribution to the Jacobian matrix
             jacobian(local_eqn,local_unknown)
              -= visc_ratio*r*(1.0+Gamma[0])
              *dpsifdx(l2,0)*dtestfdx(l,0)*W;

             jacobian(local_eqn,local_unknown)
              -= visc_ratio*r*dpsifdx(l2,1)*dtestfdx(l,1)*W;
       
             jacobian(local_eqn,local_unknown)
              -= visc_ratio*(1.0 + Gamma[0])*psif[l2]*testf[l]*W/r;
       
             //Add in the inertial terms
             //du/dt term
             jacobian(local_eqn,local_unknown) 
              -= scaled_re_st*r*node_pt(l2)->time_stepper_pt()->weight(1,0)*
              psif[l2]*testf[l]*W;

             //Convective terms
             jacobian(local_eqn,local_unknown) -=
              scaled_re*(r*psif[l2]*interpolated_dudx(0,0) 
                  + r*interpolated_u[0]*dpsifdx(l2,0)
                  + r*interpolated_u[1]*dpsifdx(l2,1))*testf[l]*W;

             //Mesh velocity terms
             if(!ALE_is_disabled)
              {
               for(unsigned k=0;k<2;k++)
                {
                 jacobian(local_eqn,local_unknown) += 
                  scaled_re_st*r*mesh_velocity[k]*dpsifdx(l2,k)*testf[l]*W;
                }
              }
            }


           //Axial velocity component
           local_unknown = nodal_local_eqn(l2,u_nodal_index[1]);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*r*Gamma[0]*dpsifdx(l2,0)*dtestfdx(l,1)*W;
             
             //Convective terms
             jacobian(local_eqn,local_unknown) -= 
              scaled_re*r*psif[l2]*interpolated_dudx(0,1)*testf[l]*W;
            }

           //Azimuthal velocity component
           local_unknown = nodal_local_eqn(l2,u_nodal_index[2]);
           if(local_unknown >= 0)
            {
             //Convective terms
             jacobian(local_eqn,local_unknown) -= 
              - scaled_re*2.0*interpolated_u[2]*psif[l2]*testf[l]*W;

             //Coriolis terms
             jacobian(local_eqn,local_unknown) +=
              2.0*r*scaled_re_inv_ro*psif[l2]*testf[l]*W;
            }
          }
         
         /*Now loop over pressure shape functions*/
         /*This is the contribution from pressure gradient*/
         for(unsigned l2=0;l2<n_pres;l2++)
          {
           local_unknown = p_local_eqn(l2);
           /*If we are at a non-zero degree of freedom in the entry*/
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown)
              += psip[l2]*(testf[l] + r*dtestfdx(l,0))*W;
            }
          }
        } /*End of Jacobian calculation*/
       
      } //End of if not boundary condition statement
       
     //SECOND (AXIAL) MOMENTUM EQUATION
     local_eqn = nodal_local_eqn(l,u_nodal_index[1]);
     //If it's not a boundary condition
     if(local_eqn >= 0)
      {
       //Add the user-defined body force terms
       //Remember to multiply by the density ratio!
       residuals[local_eqn] += 
        r*dens_ratio*body_force[1]*testf[l]*W;
       
       //Add the gravitational body force term
       residuals[local_eqn] += r*scaled_re_inv_fr*testf[l]*G[1]*W;
       
       //Add the pressure gradient term
       residuals[local_eqn]  += r*interpolated_p*dtestfdx(l,1)*W;
       
       //Add in the stress tensor terms
       //The viscosity ratio needs to go in here to ensure
       //continuity of normal stress is satisfied even in flows
       //with zero pressure gradient!
       residuals[local_eqn] -= visc_ratio*
        r*(interpolated_dudx(1,0) + Gamma[1]*interpolated_dudx(0,1))
        *dtestfdx(l,0)*W;
       
       residuals[local_eqn] -= visc_ratio*r*
        (1.0 + Gamma[1])*interpolated_dudx(1,1)*dtestfdx(l,1)*W;
       
       //Add in the inertial terms
       //du/dt term
       residuals[local_eqn] -= scaled_re_st*r*dudt[1]*testf[l]*W;

       //Convective terms
       residuals[local_eqn] -= 
        scaled_re*(r*interpolated_u[0]*interpolated_dudx(1,0) 
            + r*interpolated_u[1]*interpolated_dudx(1,1))*testf[l]*W;

       //Mesh velocity terms
       if(!ALE_is_disabled)
        {
         for(unsigned k=0;k<2;k++)
          {
           residuals[local_eqn] += 
            scaled_re_st*r*mesh_velocity[k]*interpolated_dudx(1,k)*testf[l]*W;
          }
        }
       
       //CALCULATE THE JACOBIAN
       if(flag)
        {
         //Loop over the velocity shape functions again
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           local_unknown = nodal_local_eqn(l2,u_nodal_index[0]);
           //Radial velocity component
           if(local_unknown >= 0)
            {
             //Add in the stress tensor terms
             //The viscosity ratio needs to go in here to ensure
             //continuity of normal stress is satisfied even in flows
             //with zero pressure gradient!
             jacobian(local_eqn,local_unknown) -= 
              visc_ratio*r*Gamma[1]*dpsifdx(l2,1)*dtestfdx(l,0)*W;
       
             //Convective terms
             jacobian(local_eqn,local_unknown) -= 
              scaled_re*r*psif[l2]*interpolated_dudx(1,0)*testf[l]*W;
            }

           //Axial velocity component
           local_unknown = nodal_local_eqn(l2,u_nodal_index[1]);
           if(local_unknown >= 0)
            {
             if(flag==2)
              {
               //Add the mass matrix
               mass_matrix(local_eqn,local_unknown) +=
                scaled_re_st*r*psif[l2]*testf[l]*W;
              }


             //Add in the stress tensor terms
             //The viscosity ratio needs to go in here to ensure
             //continuity of normal stress is satisfied even in flows
             //with zero pressure gradient!
             jacobian(local_eqn,local_unknown) -= 
              visc_ratio*r*dpsifdx(l2,0)*dtestfdx(l,0)*W;
       
             jacobian(local_eqn,local_unknown) -= 
              visc_ratio*r*(1.0 + Gamma[1])*dpsifdx(l2,1)*
              dtestfdx(l,1)*W;
       
             //Add in the inertial terms
             //du/dt term
             jacobian(local_eqn,local_unknown) -= 
              scaled_re_st*r*node_pt(l2)->time_stepper_pt()->weight(1,0)*
              psif[l2]*testf[l]*W;

             //Convective terms
             jacobian(local_eqn,local_unknown) -= 
              scaled_re*(r*interpolated_u[0]*dpsifdx(l2,0) 
                  + r*psif[l2]*interpolated_dudx(1,1)
                  + r*interpolated_u[1]*dpsifdx(l2,1))*testf[l]*W;
       
             
             //Mesh velocity terms
             if(!ALE_is_disabled)
              {
               for(unsigned k=0;k<2;k++)
                {
                 jacobian(local_eqn,local_unknown) += 
                  scaled_re_st*r*mesh_velocity[k]*dpsifdx(l2,k)*testf[l]*W;
                }
              }
            }
           
           //There are no azimithal terms in the axial momentum equation
          } //End of loop over velocity shape functions
       
         /*Now loop over pressure shape functions*/
         /*This is the contribution from pressure gradient*/
         for(unsigned l2=0;l2<n_pres;l2++)
          {
           local_unknown = p_local_eqn(l2);
           /*If we are at a non-zero degree of freedom in the entry*/
           if(local_unknown >= 0)
              {
               jacobian(local_eqn,local_unknown)
                += r*psip[l2]*dtestfdx(l,1)*W;
              }
          }
        } /*End of Jacobian calculation*/
       
      } //End of AXIAL MOMENTUM EQUATION

     //THIRD (AZIMUTHAL) MOMENTUM EQUATION
     local_eqn = nodal_local_eqn(l,u_nodal_index[2]);
     if(local_eqn >= 0)
      {
       //Add the user-defined body force terms
       //Remember to multiply by the density ratio!
       residuals[local_eqn] += 
        r*dens_ratio*body_force[2]*testf[l]*W;
       
       //Add the gravitational body force term
       residuals[local_eqn] += r*scaled_re_inv_fr*testf[l]*G[2]*W;
       
       //There is NO pressure gradient term
       
       //Add in the stress tensor terms
       //The viscosity ratio needs to go in here to ensure
       //continuity of normal stress is satisfied even in flows
       //with zero pressure gradient!
       residuals[local_eqn] -= visc_ratio*
        (r*interpolated_dudx(2,0) - 
         Gamma[0]*interpolated_u[2])*dtestfdx(l,0)*W;
       
       residuals[local_eqn] -= visc_ratio*r*
        interpolated_dudx(2,1)*dtestfdx(l,1)*W;

       residuals[local_eqn] -= visc_ratio*
        ((interpolated_u[2]/r) - Gamma[0]*interpolated_dudx(2,0))*testf[l]*W;


       //Add in the inertial terms
       //du/dt term
       residuals[local_eqn] -= scaled_re_st*r*dudt[2]*testf[l]*W;

       //Convective terms
       residuals[local_eqn] -= 
        scaled_re*(r*interpolated_u[0]*interpolated_dudx(2,0)
            + interpolated_u[0]*interpolated_u[2]
            + r*interpolated_u[1]*interpolated_dudx(2,1))*testf[l]*W;
       
       //Mesh velocity terms
       if(!ALE_is_disabled)
        {
         for(unsigned k=0;k<2;k++)
          {
           residuals[local_eqn] += 
            scaled_re_st*r*mesh_velocity[k]*interpolated_dudx(2,k)*testf[l]*W;
          }
        }

       //Add the Coriolis term
       residuals[local_eqn] -= 
        2.0*r*scaled_re_inv_ro*interpolated_u[0]*testf[l]*W;
       
       //CALCULATE THE JACOBIAN
       if(flag)
        {
         //Loop over the velocity shape functions again
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           //Radial velocity component
           local_unknown = nodal_local_eqn(l2,u_nodal_index[0]);
           if(local_unknown >= 0)
            {
             //Convective terms
             jacobian(local_eqn,local_unknown) -= 
              scaled_re*(r*psif[l2]*interpolated_dudx(2,0)
                  + psif[l2]*interpolated_u[2])*testf[l]*W;

             //Coriolis term
             jacobian(local_eqn,local_unknown) -=
              2.0*r*scaled_re_inv_ro*psif[l2]*testf[l]*W;
            }

           //Axial velocity component
           local_unknown = nodal_local_eqn(l2,u_nodal_index[1]);
           if(local_unknown >= 0)
            {
             //Convective terms
             jacobian(local_eqn,local_unknown) -= 
              scaled_re*r*psif[l2]*interpolated_dudx(2,1)*testf[l]*W;
            }

           //Azimuthal velocity component
           local_unknown = nodal_local_eqn(l2,u_nodal_index[2]);
           if(local_unknown >= 0)
            {
             if(flag==2)
              {
               //Add the mass matrix
               mass_matrix(local_eqn,local_unknown) +=
                scaled_re_st*r*psif[l2]*testf[l]*W;
              }

             //Add in the stress tensor terms
             //The viscosity ratio needs to go in here to ensure
             //continuity of normal stress is satisfied even in flows
             //with zero pressure gradient!
             jacobian(local_eqn,local_unknown) -= 
              visc_ratio*(r*dpsifdx(l2,0) - 
                               Gamma[0]*psif[l2])*dtestfdx(l,0)*W;
       
             jacobian(local_eqn,local_unknown) -= 
              visc_ratio*r*dpsifdx(l2,1)*dtestfdx(l,1)*W;

             jacobian(local_eqn,local_unknown) -= 
              visc_ratio*((psif[l2]/r) - Gamma[0]*dpsifdx(l2,0))
              *testf[l]*W;
             
             //Add in the inertial terms
             //du/dt term
             jacobian(local_eqn,local_unknown) -= 
              scaled_re_st*r*node_pt(l2)->time_stepper_pt()->weight(1,0)*
              psif[l2]*testf[l]*W;

             //Convective terms
             jacobian(local_eqn,local_unknown) -= 
              scaled_re*(r*interpolated_u[0]*dpsifdx(l2,0)
                  + interpolated_u[0]*psif[l2]
                  + r*interpolated_u[1]*dpsifdx(l2,1))*testf[l]*W;
             
             //Mesh velocity terms
             if(!ALE_is_disabled)
              {
               for(unsigned k=0;k<2;k++)
                {
                 jacobian(local_eqn,local_unknown) += 
                  scaled_re_st*r*mesh_velocity[k]*dpsifdx(l2,k)*testf[l]*W;
                }
              }
            }
          }

         //There are no pressure terms
        } //End of Jacobian

      } //End of AZIMUTHAL EQUATION

    } //End of loop over shape functions
      
   
   //CONTINUITY EQUATION
   //-------------------
   
   //Loop over the shape functions
   for(unsigned l=0;l<n_pres;l++)
    {
     local_eqn = p_local_eqn(l);
     //If not a boundary conditions
     if(local_eqn >= 0)
      {
       // Source term
       residuals[local_eqn] -= r*source*testp[l]*W;

       //Gradient terms
       residuals[local_eqn] += 
        (interpolated_u[0] + r*interpolated_dudx(0,0) 
         + r*interpolated_dudx(1,1))*testp[l]*W;

       
       /*CALCULATE THE JACOBIAN*/
       if(flag)
        {
         /*Loop over the velocity shape functions*/
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           //Radial velocity component
           local_unknown = nodal_local_eqn(l2,u_nodal_index[0]);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              (psif[l2] + r*dpsifdx(l2,0))*testp[l]*W;
           
            }

           //Axial velocity component
           local_unknown = nodal_local_eqn(l2,u_nodal_index[1]);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              r*dpsifdx(l2,1)*testp[l]*W;
            }
          } /*End of loop over l2*/
        } /*End of Jacobian calculation*/
       
      } //End of if not boundary condition
    } //End of loop over l
  }
 
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

///Axisymmetric Crouzeix-Raviart elements
//Set the data for the number of Variables at each node
const unsigned AxisymmetricQCrouzeixRaviartElement::
Initial_Nvalue[9]={3,3,3,3,3,3,3,3,3};

//========================================================================
/// Number of values (pinned or dofs) required at node n.
//========================================================================
unsigned AxisymmetricQCrouzeixRaviartElement::
required_nvalue(const unsigned &n) const {return Initial_Nvalue[n];}

//========================================================================
/// Compute traction at local coordinate s for outer unit normal N
//========================================================================
void AxisymmetricQCrouzeixRaviartElement::get_traction(const Vector<double>& s, 
                                                     const Vector<double>& N, 
                                                     Vector<double>& traction)
{
 AxisymmetricNavierStokesEquations::traction(s,N,traction);
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

//Axisymmetric Taylor--Hood
//Set the data for the number of Variables at each node
const unsigned AxisymmetricQTaylorHoodElement::
Initial_Nvalue[9]={4,3,4,3,3,3,4,3,4};

//Set the data for the pressure conversion array
const unsigned AxisymmetricQTaylorHoodElement::Pconv[4]={0,2,6,8};


//========================================================================
/// Compute traction at local coordinate s for outer unit normal N
//========================================================================
void AxisymmetricQTaylorHoodElement::get_traction(const Vector<double>& s, 
                                                 const Vector<double>& N, 
                                                 Vector<double>& traction)
{
 AxisymmetricNavierStokesEquations::traction(s,N,traction);
}

}
