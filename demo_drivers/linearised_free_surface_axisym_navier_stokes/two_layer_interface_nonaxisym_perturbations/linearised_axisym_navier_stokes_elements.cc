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
// Non-inline functions for linearised axisymmetric Navier-Stokes elements

// oomph-lib includes
#include "linearised_axisym_navier_stokes_elements.h"
#include "perturbed_spines.h"

namespace oomph
{
 
 //=======================================================================
 /// Linearised axisymmetric Navier--Stokes equations static data
 //=======================================================================

 // Use the stress-divergence form by default (Gamma=1)
 Vector<double> LinearisedAxisymmetricNavierStokesEquations::Gamma(2,1.0);

 // "Magic" number to indicate pressure is not stored at node
 int LinearisedAxisymmetricNavierStokesEquations::
 Pressure_not_stored_at_node = -100;

 // Physical constants default to zero
 double LinearisedAxisymmetricNavierStokesEquations::
 Default_Physical_Constant_Value = 0.0;

 // Azimuthal mode number defaults to zero
 int LinearisedAxisymmetricNavierStokesEquations::
 Default_Azimuthal_Mode_Number_Value = 0;

 // Density/viscosity ratios default to one
 double LinearisedAxisymmetricNavierStokesEquations::
 Default_Physical_Ratio_Value = 1.0;

 // Gravity vector defaults to zero
 Vector<double> LinearisedAxisymmetricNavierStokesEquations::
 Default_Gravity_Vector(3,0.0);
 

 
 //=======================================================================
 /// Output function in tecplot format: Velocities only  
 /// r, z, U^C, U^S, V^C, V^S, W^C, W^S
 /// at specified previous timestep (t=0 present; t>0 previous timestep).
 /// Specified number of plot points in each coordinate direction.
 //=======================================================================
 void LinearisedAxisymmetricNavierStokesEquations::output_veloc(
  std::ostream &outfile, const unsigned &nplot, const unsigned &t)
 {
  // Determine number of nodes in element
  const unsigned n_node = nnode();
  
  // Provide storage for local shape functions
  Shape psi(n_node);
  
  // Provide storage for vectors of local coordinates and
  // global coordinates and velocities
  Vector<double> s(2);
  Vector<double> interpolated_x(2);
  Vector<double> interpolated_u(6);
    
  // Tecplot header info
  outfile << tecplot_zone_string(nplot);
  
  // Determine number of plot points
  const unsigned n_plot_points = nplot_points(nplot);

  // Loop over plot points
  for(unsigned iplot=0;iplot<n_plot_points;iplot++)
   {
    // Get local coordinates of plot point
    get_s_plot(iplot,nplot,s);
    
    // Get shape functions
    shape(s,psi);
    
    // Loop over coordinate directions
    for(unsigned i=0;i<2;i++) 
     {
      // Initialise global coordinate
      interpolated_x[i]=0.0;
      
      // Loop over the local nodes and sum
      for(unsigned l=0;l<n_node;l++)
       {
        interpolated_x[i] += nodal_position(t,l,i)*psi[l];
       }
     }
    
    // Loop over the velocity components
    for(unsigned i=0;i<6;i++) 
     {
      // Get the index at which the velocity is stored
      const unsigned u_nodal_index = u_index_lin_axi_nst(i);

      // Initialise velocity
      interpolated_u[i]=0.0;
      
      // Loop over the local nodes and sum
      for(unsigned l=0;l<n_node;l++)
       {
        interpolated_u[i] += nodal_value(t,l,u_nodal_index)*psi[l];
       }
     }
    
    // Output global coordinates to file
    for(unsigned i=0;i<2;i++) { outfile << interpolated_x[i] << " "; }
    
    // Output velocities to file
    for(unsigned i=0;i<6;i++) { outfile << interpolated_u[i] << " "; }
    
    outfile << std::endl;   
   }
  
  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,nplot);
  
 } // End of output_veloc



 //=======================================================================
 /// Output function in tecplot format:
 /// r, z, R^C, R^S, Z^C, Z^S, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S
 /// Specified number of plot points in each coordinate direction.
 //=======================================================================
 void LinearisedAxisymmetricNavierStokesEquations::output(
  std::ostream &outfile, const unsigned &nplot)
 {
  // Provide storage for vector of local coordinates
  Vector<double> s(2);
  
  // Tecplot header info
  outfile << tecplot_zone_string(nplot);
 
  // Determine number of plot points
  const unsigned n_plot_points = nplot_points(nplot);

  // Loop over plot points
  for(unsigned iplot=0;iplot<n_plot_points;iplot++)
   {
    // Get local coordinates of plot point
    get_s_plot(iplot,nplot,s);

    // Output global coordinates to file
    for(unsigned i=0;i<2;i++) { outfile << interpolated_x(s,i) << " "; }
  
    // Output perturbations to nodal positions to file
    for(unsigned i=0;i<4;i++)
     {
      outfile <<
       interpolated_nodal_position_perturbation_lin_axi_nst(s,i) << " ";
     }

    //  Output velocities to file
    for(unsigned i=0;i<6;i++)
     {
      outfile << interpolated_u_lin_axi_nst(s,i) << " ";
     }
   
    // Output pressure to file
    for(unsigned i=0;i<2;i++)
     {
      outfile << interpolated_p_lin_axi_nst(s,i)  << " ";
     }

    outfile << std::endl;   
   }
  outfile << std::endl;
  
  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,nplot);
  
 } // End of output



 //=======================================================================
 /// Output function in tecplot format: 
 /// r, z, R^C, R^S, Z^C, Z^S, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S
 /// Specified number of plot points in each coordinate direction.
 //=======================================================================
 void LinearisedAxisymmetricNavierStokesEquations::output(
  FILE* file_pt, const unsigned &nplot)
 {
  // Provide storage for vector of local coordinates
  Vector<double> s(2);
  
  // Tecplot header info
  fprintf(file_pt,"%s ",tecplot_zone_string(nplot).c_str());
  
  // Determine number of plot points
  const unsigned n_plot_points = nplot_points(nplot);  
  
  // Loop over plot points
  for(unsigned iplot=0;iplot<n_plot_points;iplot++)
   {
    // Get local coordinates of plot point
    get_s_plot(iplot,nplot,s);
  
    // Output global coordinates to file
    for(unsigned i=0;i<2;i++) { fprintf(file_pt,"%g ",interpolated_x(s,i)); }

    // Output perturbations to nodal positions to file
    for(unsigned i=0;i<4;i++)
     {
      fprintf(file_pt,"%g ",
              interpolated_nodal_position_perturbation_lin_axi_nst
              (s,i));
     }

    //  Output velocities to file
    for(unsigned i=0;i<6;i++)
     {
      fprintf(file_pt,"%g ",interpolated_u_lin_axi_nst(s,i));
     }
   
    // Output pressure to file
    for(unsigned i=0;i<2;i++)
     {
      fprintf(file_pt,"%g ",interpolated_p_lin_axi_nst(s,i));
     }

    fprintf(file_pt,"\n");
   }

  fprintf(file_pt,"\n");

  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(file_pt,nplot);
  
 } // End of output(...)



 
 //=======================================================================
 /// Get strain-rate tensor: \f$ e_{ij} \f$  where 
 /// \f$ i,j = r,z,\theta \f$ (in that order). We evaluate this tensor
 /// at a value of theta such that the product of theta and the azimuthal
 /// mode number (k) gives \f$ \pi/4 \f$. Therefore
 /// \f$ \cos(k \theta) = \sin(k \theta) = 1/\sqrt{2} \f$.
 //=======================================================================
 void LinearisedAxisymmetricNavierStokesEquations::
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
  
  // Determine number of nodes in element
  const unsigned n_node = nnode();

  //Set up memory for the shape and test functions
  Shape psi(n_node);
  DShape dpsidx(n_node,2);
  
  //Call the derivatives of the shape functions
  dshape_eulerian(s,psi,dpsidx);
  
  // Radius
  double interpolated_r = 0.0;
  
  // Velocity components and their derivatives
  double UC = 0.0, US = 0.0;
  double dUCdr = 0.0, dUSdr = 0.0;
  double dUCdz = 0.0, dUSdz = 0.0;
  double WC = 0.0, WS = 0.0;
  double dWCdr = 0.0, dWSdr = 0.0;
  double dWCdz = 0.0, dWSdz = 0.0;
  double VC = 0.0, VS = 0.0;
  double dVCdr = 0.0, dVSdr = 0.0;
  double dVCdz = 0.0, dVSdz = 0.0;
  
  //Get the local storage for the indices
  unsigned u_nodal_index[6];
  for(unsigned i=0;i<6;++i) {u_nodal_index[i] = u_index_lin_axi_nst(i);}
  
  // Loop over nodes to assemble velocities and their derivatives
  // w.r.t. r and z (x_0 and x_1)
  for(unsigned l=0;l<n_node;l++) 
   {   
    interpolated_r += nodal_position(l,0)*psi[l];
    
    UC += nodal_value(l,u_nodal_index[0])*psi[l];
    US += nodal_value(l,u_nodal_index[1])*psi[l];
    WC += nodal_value(l,u_nodal_index[2])*psi[l];
    WS += nodal_value(l,u_nodal_index[3])*psi[l];
    VC += nodal_value(l,u_nodal_index[4])*psi[l];
    VS += nodal_value(l,u_nodal_index[4])*psi[l];
    
    dUCdr += nodal_value(l,u_nodal_index[0])*dpsidx(l,0);
    dUSdr += nodal_value(l,u_nodal_index[1])*dpsidx(l,0);
    dWCdr += nodal_value(l,u_nodal_index[2])*dpsidx(l,0);
    dWSdr += nodal_value(l,u_nodal_index[3])*dpsidx(l,0);
    dVCdr += nodal_value(l,u_nodal_index[4])*dpsidx(l,0);
    dVSdr += nodal_value(l,u_nodal_index[5])*dpsidx(l,0);

    dUCdz += nodal_value(l,u_nodal_index[0])*dpsidx(l,1);
    dUSdz += nodal_value(l,u_nodal_index[1])*dpsidx(l,1);
    dWCdz += nodal_value(l,u_nodal_index[2])*dpsidx(l,1);
    dWSdz += nodal_value(l,u_nodal_index[3])*dpsidx(l,1);
    dVCdz += nodal_value(l,u_nodal_index[4])*dpsidx(l,1);
    dVSdz += nodal_value(l,u_nodal_index[5])*dpsidx(l,1);
   }

  // Cache azimuthal mode number
  const int k = this->azimuthal_mode_number();

  // We wish to evaluate the strain-rate tensor at a value of theta
  // such that k*theta = pi/4 radians. That way we pick up equal
  // contributions from the real and imaginary parts of the velocities.
  // sin(pi/4) = cos(pi/4) = 1/sqrt(2)
  const double cosktheta = 1.0/sqrt(2);
  const double sinktheta = cosktheta;
  
  // Assemble velocities and their derivatives w.r.t. r, z and theta
  // from real and imaginary parts
  const double ur = UC*cosktheta + US*sinktheta;
  const double utheta = VC*cosktheta + VS*sinktheta;

  const double durdr = dUCdr*cosktheta + dUSdr*sinktheta;
  const double durdz = dUCdz*cosktheta + dUSdz*sinktheta;
  const double durdtheta = k*US*cosktheta - k*UC*sinktheta;

  const double duzdr = dWCdr*cosktheta + dWSdr*sinktheta;
  const double duzdz = dWCdz*cosktheta + dWSdz*sinktheta;
  const double duzdtheta = k*WS*cosktheta - k*WC*sinktheta;

  const double duthetadr = dVCdr*cosktheta + dVSdr*sinktheta;
  const double duthetadz = dVCdz*cosktheta + dVSdz*sinktheta;
  const double duthetadtheta = k*VS*cosktheta - k*VC*sinktheta;

  // Assign strain rates without negative powers of the radius
  // and zero those with:
  strainrate(0,0)=durdr;
  strainrate(0,1)=0.5*(durdz+duzdr);
  strainrate(1,0)=strainrate(0,1);
  strainrate(0,2)=0.5*duthetadr;
  strainrate(2,0)=strainrate(0,2);
  strainrate(1,1)=duzdz;
  strainrate(1,2)=0.5*duthetadz;
  strainrate(2,1)=strainrate(1,2);
  strainrate(2,2)=0.0;
  
  
  // Overwrite the strain rates with negative powers of the radius
  // unless we're at the origin
  if (std::abs(interpolated_r)>1.0e-16)
   {
    double inverse_radius=1.0/interpolated_r;
    strainrate(0,2)=0.5*(duthetadr + inverse_radius*(durdtheta - utheta));
    strainrate(2,0)=strainrate(0,2);
    strainrate(2,2)=inverse_radius*(ur + duthetadtheta);
    strainrate(1,2)=0.5*(duthetadz + inverse_radius*duzdtheta);
    strainrate(2,1)=strainrate(1,2);
   }
  
 } // End of strain_rate
  


 //======================================================================= 
 /// Get integral of kinetic energy over element. Also return the
 /// derivative of the integral of the kinetic energy w.r.t. time.
 /// Note that this is the "raw" kinetic energy in the sense
 /// that the density ratio has not been included. In problems
 /// with two or more fluids the user will have to remember to
 /// premultiply certain elements by the appropriate density ratio.
 /// ratio.
 //==============================================================
 void LinearisedAxisymmetricNavierStokesEquations::dkin_energy_dt
 (double& dkin_en_dt, double& kin_en) const
 {  
  // Initialise kinetic energy and its deriv w.r.t. time
  kin_en = 0.0; dkin_en_dt = 0.0;

  // Determine number of nodes in the element
  const unsigned n_node = nnode();

  // Set up memory for the fluid (f) shape functions and their derivatives
  // w.r.t. local coordinates (s).
  Shape psif(n_node);
  DShape dpsifds(n_node,2);
  
  // Set the Vector to hold local coordinates
  Vector<double> s(2);
  
  // Get the nodal indices at which the velocity is stored
  unsigned u_nodal_index[6];
  for(unsigned i=0;i<6;++i) { u_nodal_index[i] = u_index_lin_axi_nst(i); }

  // Determine number of integration points
  const unsigned n_intpt = integral_pt()->nweight();
  
  // Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {   
    // Compute the local coordinates s of the integration point
    for(unsigned i=0;i<2;i++) { s[i] = integral_pt()->knot(ipt,i); }
    
    // Get the integral weight
    const double w = integral_pt()->weight(ipt);
    
    // Calculate the derivatives of the fluid shape functions w.r.t. the
    // local coordinates
    this->dshape_local_at_knot(ipt,psif,dpsifds);

    // Get Jacobian of mapping
    const double J = J_eulerian(s);
    
    // Allocate storage for the r-position and its deriv w.r.t. time
    double interpolated_r = 0.0;
    double interpolated_drdt = 0.0;

    // Allocate storage for the derivatives of the positions
    // w.r.t. local coordinates (s_1 and s_2)
    DenseMatrix<double> interpolated_dxds(2,2,0.0);

    // Allocate storage for the derivative w.r.t. time of the spatial
    // derivatives of the positions
    DenseMatrix<double> interpolated_d_dxds_dt(2,2,0.0);

    // Allocate storage for the velocity components (six of these)
    // and their derivatives w.r.t. time
    Vector<double> interpolated_u(6,0.0);
    Vector<double> dudt(6,0.0);

    // Loop over the element's nodes
    for(unsigned l=0;l<n_node;l++) 
     {
      // Cache the shape function
      const double psif_ = psif(l);

      // Loop over the two coordinate directions
      for(unsigned i=0;i<2;i++)
       {
        // Loop over the two coordinate directions (for derivatives)
        for(unsigned j=0;j<2;j++)
         {
          interpolated_dxds(i,j) += this->nodal_position(l,i)*dpsifds(l,j);

          interpolated_d_dxds_dt(i,j) += 
           this->dnodal_position_dt(l,i)*dpsifds(l,j);
         }
       }

      // Calculate the r-position
      interpolated_r += this->nodal_position(l,0)*psif_;

      // Calculate the derivative of the r-position w.r.t. time
      interpolated_drdt += this->dnodal_position_dt(l,0)*psif_;
      
      // Loop over the six velocity components
      for(unsigned i=0;i<6;i++)
       {
        // Get the value
        const double u_value = this->nodal_value(l,u_nodal_index[i]);
        
        // Add contribution
        interpolated_u[i] += u_value*psif_;
        
        // Add contribution to dudt
        dudt[i] += du_dt_lin_axi_nst(l,i)*psif_;
       }
     } // End of loop over the element's nodes

    // Compute sum of square of velocity components and sum of product of
    // velocity components and their time derivatives
    double veloc_squared = 0.0;
    double veloc_multiplied_by_time_deriv = 0.0;
    for(unsigned i=0;i<6;i++) 
     {
      veloc_squared += interpolated_u[i]*interpolated_u[i];
      veloc_multiplied_by_time_deriv += interpolated_u[i]*dudt[i];
     }

    // Compute the derivative of the Jacobian of the mapping w.r.t. time
    const double dJdt = interpolated_d_dxds_dt(0,0)*interpolated_dxds(1,1)
     + interpolated_d_dxds_dt(1,1)*interpolated_dxds(0,0)
     + interpolated_d_dxds_dt(1,0)*interpolated_dxds(0,1)
     + interpolated_d_dxds_dt(0,1)*interpolated_dxds(1,0);

    // Add contribution to kinetic energy (no density ratio)
    kin_en += veloc_squared*interpolated_r*J*w;

    // Add contributions to deriv of kinetic energy w.r.t. time
    dkin_en_dt += 2.0*veloc_multiplied_by_time_deriv*interpolated_r*J*w;
    dkin_en_dt += veloc_squared*interpolated_drdt*J*w;
    dkin_en_dt += veloc_squared*interpolated_r*dJdt*w;
   }
 } // End of dkin_energy_dt
 



 //======================================================================= 
 /// Compute the residuals for the linearised axisymmetric Navier--Stokes 
 /// equations; flag=1(or 0): do (or don't) compute the Jacobian as well.
 //======================================================================= 
 void LinearisedAxisymmetricNavierStokesEquations::
 fill_in_generic_residual_contribution_lin_axi_nst(
  Vector<double> &residuals, 
  DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix,
  unsigned flag)
 {
  // Get the time from the first node in the element
  const double time = this->node_pt(0)->time_stepper_pt()->time();

  // Determine number of nodes in the element
  const unsigned n_node = nnode();

  // Determine number of nodes along one edge
  const unsigned n_p = nnode_1d();
 
  // Determine how many pressure values there are associated with
  // a single pressure component
  const unsigned n_pres = npres_lin_axi_nst();
  
  // Get the nodal indices at which the velocity is stored
  unsigned u_nodal_index[6];
  for(unsigned i=0;i<6;++i)
   {
    u_nodal_index[i] = u_index_lin_axi_nst(i);
   }

  // Get the nodal indices at which the perturbations to the nodal
  // positions are stored
  unsigned xhat_nodal_index[4];
  for(unsigned i=0;i<4;++i)
   {
    xhat_nodal_index[i] = xhat_index_lin_axi_nst(i);
   }
  
  // Set up memory for the fluid shape and test functions
  // Note that there are two spatial dimensions, r and z, in this problem
  Shape psif(n_node), testf(n_node);
  DShape dpsifds(n_node,2), dpsifdx(n_node,2);
  DShape dtestfds(n_node,2), dtestfdx(n_node,2);
  
  // Set up memory for the pressure (p) shape and test functions
  Shape psip(n_pres), testp(n_pres);
  
  // Determine number of integration points
  const unsigned n_intpt = integral_pt()->nweight();
  
  // Set up memory for the vector to hold local coordinates (two dimensions)
  Vector<double> s(2);

  // Get physical variables from the element
  // (Reynolds number must be multiplied by the density ratio)
  const double scaled_re = re()*density_ratio();
  const double scaled_re_st = re_st()*density_ratio();
  const double scaled_re_inv_fr = re_invfr()*density_ratio();
  const double visc_ratio = viscosity_ratio();
  const Vector<double> G = g();
  const int k = azimuthal_mode_number();

  // Integers used to store the local equation and unknown numbers
  int local_eqn = 0, local_unknown = 0;

  // Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    // Assign values of the local coordinates s
    for(unsigned i=0;i<2;i++) { s[i] = integral_pt()->knot(ipt,i); }

    // Get the integral weight
    const double w = integral_pt()->weight(ipt);
    
    // Calculate the derivatives of the fluid shape functions w.r.t. the
    // local coordinates
    this->dshape_local_at_knot(ipt,psif,dpsifds);

    // Set the derivatives of the test functions w.r.t. local coords
    // equal to that of the shape functions
    dtestfds = dpsifds;

    // Calculate the fluid shape and test functions, and their derivatives
    // w.r.t. the global coordinates
    const double Jbar = 
     dshape_and_dtest_eulerian_at_knot_lin_axi_nst(ipt,psif,dpsifdx,
                                                          testf,dtestfdx);
    
    // Calculate the pressure shape and test functions
    pshape_lin_axi_nst(s,psip,testp);
    
    // Allocate storage for the position and the derivative of the 
    // mesh positions w.r.t. time
    Vector<double> interpolated_x(2,0.0);
    Vector<double> mesh_velocity(2,0.0);

    // Allocate storage for the derivatives of the unperturbed positions
    // w.r.t. local coordinates (s_1 and s_2)
    DenseMatrix<double> interpolated_dxbar_ds(2,2,0.0);

    // Allocate storage for the perturbed position xhat
    Vector<double> interpolated_xhat(4,0.0);

    // Allocate storage for the derivative of the perturbed position
    // w.r.t. time
    Vector<double> dxhat_dt(4,0.0);

    // Allocate storage for the derivatives of the perturbed positions
    // w.r.t. local coordinates (s_1 and s_2)
    DenseMatrix<double> interpolated_dxhat_ds(4,2,0.0);

    // Allocate storage for the velocity components (six of these)
    // and their derivatives w.r.t. time
    Vector<double> interpolated_u(6,0.0);
    Vector<double> dudt(6,0.0);
    
    // Allocate storage for the pressure components (two of these)
    Vector<double> interpolated_p(2,0.0);

    // Allocate storage for the derivatives of the velocity components
    // w.r.t. local coordinates (s_1 and s_2)
    DenseMatrix<double> interpolated_duds(6,2,0.0);
   
    // Calculate pressure at the integration point
    // -------------------------------------------

    // Loop over pressure degrees of freedom (associated with a single
    // pressure component) in the element
    for(unsigned l=0;l<n_pres;l++)
     {
      // Cache the shape function
      const double psip_ = psip(l);

      // Loop over the two pressure components
      for(unsigned i=0;i<2;i++)
       {
        // Get the value
        const double p_value = this->p_lin_axi_nst(l,i);

        // Add contribution
        interpolated_p[i] += p_value*psip_;
       }
     } // End of loop over the pressure degrees of freedom in the element
   
    // Calculate eulerian positions, perturbations to these positions,
    // ---------------------------------------------------------------
    // velocities and their derivatives at the integration point
    // ---------------------------------------------------------

    // Loop over the element's nodes
    for(unsigned l=0;l<n_node;l++) 
     {
      // Cache the shape function
      const double psif_ = psif(l);

      // Loop over the two coordinate directions
      for(unsigned i=0;i<2;i++)
       {
        // Calculate the unperturbed position xbar
        interpolated_x[i] += this->raw_nodal_position(l,i)*psif_;

        // Loop over the two coordinate directions (for derivatives)
        for(unsigned j=0;j<2;j++)
         {
          interpolated_dxbar_ds(i,j) +=
           this->raw_nodal_position(l,i)*dpsifds(l,j);
         }
       }
      
      // Loop over the four position perturbations R_k^C, R_k^S, Z_k^C, Z_k^S
      for(unsigned i=0;i<4;i++)
       {
        // Get the value
        const double xhat_value = this->raw_nodal_value(l,xhat_nodal_index[i]);

        // Add contribution
        interpolated_xhat[i] += xhat_value*psif_;

        // Loop over the two coordinate directions (for derivatives)
        for(unsigned j=0;j<2;j++)
         {
          interpolated_dxhat_ds(i,j) += xhat_value*dpsifds(l,j);
         }
       }
      
      // Loop over the six velocity components
      for(unsigned i=0;i<6;i++)
       {
        // Get the value
        const double u_value = this->raw_nodal_value(l,u_nodal_index[i]);
        
        // Add contribution
        interpolated_u[i] += u_value*psif_;
        
        // Add contribution to dudt
        dudt[i] += du_dt_lin_axi_nst(l,i)*psif_;
        
        // Loop over the two coordinate directions (for derivatives)
        for(unsigned j=0;j<2;j++)
         {
          interpolated_duds(i,j) += u_value*dpsifds(l,j);
         }
       }
     } // End of loop over the element's nodes

    // Get the mesh velocity if ALE is enabled
    if(!ALE_is_disabled)
     {
      // Loop over the element's nodes
      for(unsigned l=0;l<n_node;l++) 
       {
        // Loop over the two coordinate directions
        for(unsigned i=0;i<2;i++)
         {
          // Calculate the derivative of the unperturbed position xbar
          // w.r.t. time
          mesh_velocity[i] += this->raw_dnodal_position_dt(l,i)*psif(l);
         }

        // Loop over the four perturbed positions R_k^C, R_k^S, Z_k^C, Z_k^S
        for(unsigned i=0;i<4;i++)
         {
          // Calculate the derivative of the perturbed position xhat
          // w.r.t. time
          dxhat_dt[i] +=
           this->dnodal_position_perturbation_dt_lin_axi_nst(l,i)
           *psif(l);
         }
       }
     }

    // Compute the cosine part of the perturbed jacobian at this ipt
    const double JhatC =
     interpolated_dxbar_ds(0,0)*interpolated_dxhat_ds(2,1)
     + interpolated_dxbar_ds(1,1)*interpolated_dxhat_ds(0,0)
     - interpolated_dxbar_ds(0,1)*interpolated_dxhat_ds(2,0)
     - interpolated_dxbar_ds(1,0)*interpolated_dxhat_ds(0,1);

    // Compute the sine part of the perturbed jacobian at this ipt
    const double JhatS =
     interpolated_dxbar_ds(0,0)*interpolated_dxhat_ds(3,1)
     + interpolated_dxbar_ds(1,1)*interpolated_dxhat_ds(1,0)
     - interpolated_dxbar_ds(0,1)*interpolated_dxhat_ds(3,0)
     - interpolated_dxbar_ds(1,0)*interpolated_dxhat_ds(1,1);

    // Determine the inverse of the unperturbed jacobian
    const double invJbar = 1.0/Jbar;

    // Get the user-defined (base flow) body force terms
    Vector<double> body_force(3);
    get_body_force_base_flow(time,ipt,interpolated_x,body_force);
   
    // Get the user-defined (base flow) source function
    const double source = get_source_base_flow(time,ipt,interpolated_x);

    // Get velocities and their derivatives from base flow problem
    // -----------------------------------------------------------

    // Allocate storage for the velocity components of the base state
    // solution (initialise to zero)
    Vector<double> base_flow_u(3,0.0);

    // Get the base state solution velocity components
    get_base_flow_u(time,ipt,interpolated_x,base_flow_u);

    // Allocate storage for the derivatives of the base state solution's
    // velocity components w.r.t. local coordinates (s_1 and s_2) and
    // global coordinates (r and z)
    // N.B. the derivatives of the base flow components w.r.t. the
    // azimuthal coordinate direction (theta) are always zero since the
    // base flow is axisymmetric
    DenseMatrix<double> base_flow_duds(3,2,0.0);
    DenseMatrix<double> base_flow_dudx(3,2,0.0);

    // Get the derivatives of the base state solution
    // velocity components w.r.t. local and global coordinates
    get_base_flow_duds(time,ipt,interpolated_x,base_flow_duds);
    get_base_flow_dudx(time,ipt,interpolated_x,base_flow_dudx);

    // Allocate storage for the base state pressure at the current
    // integration point
    double base_flow_p = 0.0;

    // Allocate storage for the derivatives of the base state solution
    // velocity components w.r.t. time
    Vector<double> base_flow_dudt(3,0.0);

    // If ALE is enabled, get the base state pressure and the derivatives
    // of the base state velocity w.r.t. time (only needed in this case)
    if(!ALE_is_disabled)
     {
      get_base_flow_p(time,ipt,interpolated_x,base_flow_p);
      get_base_flow_dudt(time,ipt,interpolated_x,base_flow_dudt);
     }

    // Compute the following quantities
    const double interpolated_dUdRC =
     invJbar*(interpolated_dxbar_ds(1,1)*interpolated_duds(0,0)
              + base_flow_duds(0,0)*interpolated_dxhat_ds(2,1)
              - interpolated_dxbar_ds(1,0)*interpolated_duds(0,1)
              - base_flow_duds(0,1)*interpolated_dxhat_ds(2,0));
    const double interpolated_dUdRS =
     invJbar*(interpolated_dxbar_ds(1,1)*interpolated_duds(1,0)
              + base_flow_duds(0,0)*interpolated_dxhat_ds(3,1)
              - interpolated_dxbar_ds(1,0)*interpolated_duds(1,1)
              - base_flow_duds(0,1)*interpolated_dxhat_ds(3,0));
    const double interpolated_dWdRC =
     invJbar*(interpolated_dxbar_ds(1,1)*interpolated_duds(2,0)
              + base_flow_duds(1,0)*interpolated_dxhat_ds(2,1)
              - interpolated_dxbar_ds(1,0)*interpolated_duds(2,1)
              - base_flow_duds(1,1)*interpolated_dxhat_ds(2,0));
    const double interpolated_dWdRS =
     invJbar*(interpolated_dxbar_ds(1,1)*interpolated_duds(3,0)
              + base_flow_duds(1,0)*interpolated_dxhat_ds(3,1)
              - interpolated_dxbar_ds(1,0)*interpolated_duds(3,1)
              - base_flow_duds(1,1)*interpolated_dxhat_ds(3,0));
    const double interpolated_dVdRC =
     invJbar*(interpolated_dxbar_ds(1,1)*interpolated_duds(4,0)
              + base_flow_duds(2,0)*interpolated_dxhat_ds(2,1)
              - interpolated_dxbar_ds(1,0)*interpolated_duds(4,1)
              - base_flow_duds(2,1)*interpolated_dxhat_ds(2,0));
    const double interpolated_dVdRS =
     invJbar*(interpolated_dxbar_ds(1,1)*interpolated_duds(5,0)
              + base_flow_duds(2,0)*interpolated_dxhat_ds(3,1)
              - interpolated_dxbar_ds(1,0)*interpolated_duds(5,1)
              - base_flow_duds(2,1)*interpolated_dxhat_ds(3,0));
    const double interpolated_dUdZC =
     invJbar*(interpolated_dxbar_ds(0,0)*interpolated_duds(0,1)
              + base_flow_duds(0,1)*interpolated_dxhat_ds(0,0)
              - interpolated_dxbar_ds(0,1)*interpolated_duds(0,0)
              - base_flow_duds(0,0)*interpolated_dxhat_ds(0,1));
    const double interpolated_dUdZS =
     invJbar*(interpolated_dxbar_ds(0,0)*interpolated_duds(1,1)
              + base_flow_duds(0,1)*interpolated_dxhat_ds(1,0)
              - interpolated_dxbar_ds(0,1)*interpolated_duds(1,0)
              - base_flow_duds(0,0)*interpolated_dxhat_ds(1,1));
    const double interpolated_dWdZC =
     invJbar*(interpolated_dxbar_ds(0,0)*interpolated_duds(2,1)
              + base_flow_duds(1,1)*interpolated_dxhat_ds(0,0)
              - interpolated_dxbar_ds(0,1)*interpolated_duds(2,0)
              - base_flow_duds(1,0)*interpolated_dxhat_ds(0,1));
    const double interpolated_dWdZS =
     invJbar*(interpolated_dxbar_ds(0,0)*interpolated_duds(3,1)
              + base_flow_duds(1,1)*interpolated_dxhat_ds(1,0)
              - interpolated_dxbar_ds(0,1)*interpolated_duds(3,0)
              - base_flow_duds(1,0)*interpolated_dxhat_ds(1,1));
    const double interpolated_dVdZC =
     invJbar*(interpolated_dxbar_ds(0,0)*interpolated_duds(4,1)
              + base_flow_duds(2,1)*interpolated_dxhat_ds(0,0)
              - interpolated_dxbar_ds(0,1)*interpolated_duds(4,0)
              - base_flow_duds(2,0)*interpolated_dxhat_ds(0,1));
    const double interpolated_dVdZS =
     invJbar*(interpolated_dxbar_ds(0,0)*interpolated_duds(5,1)
              + base_flow_duds(2,1)*interpolated_dxhat_ds(1,0)
              - interpolated_dxbar_ds(0,1)*interpolated_duds(5,0)
              - base_flow_duds(2,0)*interpolated_dxhat_ds(1,1));

    // Define the following useful quantities...
    Vector<double> group_A(n_node,0.0);
    Vector<double> group_B(n_node,0.0);
    Vector<double> group_C(n_node,0.0);
    Vector<double> group_D(n_node,0.0);
    Vector<double> group_E(n_node,0.0);
    DenseMatrix<double> group_F(n_node,n_node,0.0);

    // Loop over the element's nodes
    for(unsigned l=0;l<n_node;l++) 
     {
      group_A[l] = interpolated_dxbar_ds(0,0)*dpsifds(l,1)
       - interpolated_dxbar_ds(0,1)*dpsifds(l,0);
      group_B[l] = interpolated_dxbar_ds(1,1)*dpsifds(l,0)
       - interpolated_dxbar_ds(1,0)*dpsifds(l,1);
      group_C[l] = base_flow_duds(0,0)*dpsifds(l,1)
       - base_flow_duds(0,1)*dpsifds(l,0);
      group_D[l] = base_flow_duds(1,0)*dpsifds(l,1)
       - base_flow_duds(1,1)*dpsifds(l,0);
      group_E[l] = base_flow_duds(2,0)*dpsifds(l,1)
       - base_flow_duds(2,1)*dpsifds(l,0); 

      // Loop over the element's nodes again
      for(unsigned l2=0;l2<n_node;l2++)
       {
        group_F(l,l2) = dtestfds(l,0)*dpsifds(l2,1)
         - dtestfds(l,1)*dpsifds(l2,0);
       }
     }

    // Allocate storage for the "effective" perturbed spine node fraction
    Vector<double> effective_perturbed_spine_node_fraction(n_node,0.0);

    // Loop over the element's nodes
    for(unsigned l=0;l<n_node;l++) 
     {
      // Upcast to a PerturbedSpineNode
      PerturbedSpineNode* nod_pt =
       dynamic_cast<PerturbedSpineNode*>(this->node_pt(l));

      // Cache node update fct id
      const unsigned cached_node_update_fct_id = nod_pt->node_update_fct_id();

      // Set effective perturbed spine node fraction -- this is different
      // depending on whether node l uses the lower (0) or upper (1) node
      // update function
      if(cached_node_update_fct_id==0)
       {
        effective_perturbed_spine_node_fraction[l] = nod_pt->fraction();
       }
      else if(cached_node_update_fct_id==1)
       {
        effective_perturbed_spine_node_fraction[l] = 1.0 - nod_pt->fraction();
       }
      else
       {
        throw OomphLibError("Incorrect node_update_fct_id",
                            "LinearisedAxisymmetricNavierStokesEquations::fill_in_generic_residual_contribution_lin_axi_nst()",
                            OOMPH_EXCEPTION_LOCATION);
       }
     }

    // Cache base flow velocities and their derivatives
    const double base_flow_ur = base_flow_u[0];
    const double base_flow_uz = base_flow_u[1];
    const double base_flow_utheta = base_flow_u[2];
    const double base_flow_durdr = base_flow_dudx(0,0);
    const double base_flow_durdz = base_flow_dudx(0,1);
    const double base_flow_duzdr = base_flow_dudx(1,0);
    const double base_flow_duzdz = base_flow_dudx(1,1);
    const double base_flow_duthetadr = base_flow_dudx(2,0);
    const double base_flow_duthetadz = base_flow_dudx(2,1);
    const double base_flow_durdt = base_flow_dudt[0];
    const double base_flow_duzdt = base_flow_dudt[1];
    const double base_flow_duthetadt = base_flow_dudt[2];

    // Cache r-component of position
    const double r = interpolated_x[0];

    // Cache perturbations to nodal positions
    const double interpolated_RC = interpolated_xhat[0];
    const double interpolated_RS = interpolated_xhat[1];
    const double interpolated_ZC = interpolated_xhat[2];
    const double interpolated_ZS = interpolated_xhat[3];

    // Cache temporal derivatives of the perturbations to nodal positions
    const double dRCdt = dxhat_dt[0];
    const double dRSdt = dxhat_dt[1];
    const double dZCdt = dxhat_dt[2];
    const double dZSdt = dxhat_dt[3];

    // Cache unknowns
    const double interpolated_UC = interpolated_u[0];
    const double interpolated_US = interpolated_u[1];
    const double interpolated_WC = interpolated_u[2];
    const double interpolated_WS = interpolated_u[3];
    const double interpolated_VC = interpolated_u[4];
    const double interpolated_VS = interpolated_u[5];
    const double interpolated_PC = interpolated_p[0];
    const double interpolated_PS = interpolated_p[1];
    
    // Cache temporal derivatives of the unknowns
    const double dUCdt = dudt[0];
    const double dUSdt = dudt[1];
    const double dWCdt = dudt[2];
    const double dWSdt = dudt[3];
    const double dVCdt = dudt[4];
    const double dVSdt = dudt[5];

    // ==================
    // MOMENTUM EQUATIONS
    // ==================

    // Loop over the fluid test functions
    for(unsigned l=0;l<n_node;l++)
     {
      // Cache test functions and their derivatives
      const double testf_ = testf(l);
      const double dtestfdr = dtestfdx(l,0);
      const double dtestfdz = dtestfdx(l,1);

      // Compute the following useful quantities...
      const double dtestfdRC = invJbar*
       (dtestfds(l,0)*interpolated_dxhat_ds(2,1)
        - dtestfds(l,1)*interpolated_dxhat_ds(2,0));
      const double dtestfdRS = invJbar*
       (dtestfds(l,0)*interpolated_dxhat_ds(3,1)
        - dtestfds(l,1)*interpolated_dxhat_ds(3,0));
      const double dtestfdZC = invJbar*
       (dtestfds(l,1)*interpolated_dxhat_ds(0,0)
        - dtestfds(l,0)*interpolated_dxhat_ds(0,1));
      const double dtestfdZS = invJbar*
       (dtestfds(l,1)*interpolated_dxhat_ds(1,0)
        - dtestfds(l,0)*interpolated_dxhat_ds(1,1));

      // ---------------------------------------------
      // FIRST (RADIAL) MOMENTUM EQUATION: COSINE PART
      // ---------------------------------------------

      // Get local equation number of first velocity value at this node
      local_eqn = nodal_local_eqn(l,u_nodal_index[0]);

      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] -= scaled_re_st*r*dUCdt*testf_*Jbar*w;
        residuals[local_eqn] -= 
         scaled_re_st*interpolated_RC*base_flow_durdt*testf_*Jbar*w;
        residuals[local_eqn] -= 
         scaled_re*r*base_flow_ur*interpolated_dUdRC*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[0]*interpolated_dUdRC*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_UC*base_flow_durdr*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dRCdt*base_flow_durdr*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RC*base_flow_ur*base_flow_durdr*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RC*mesh_velocity[0]*base_flow_durdr
         *testf_*Jbar*w;
        residuals[local_eqn] -=
         k*scaled_re*base_flow_utheta*interpolated_US*testf_*Jbar*w;
        residuals[local_eqn] +=
         k*scaled_re*base_flow_utheta*base_flow_durdr*interpolated_RS
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         k*scaled_re*base_flow_utheta*base_flow_durdz*interpolated_ZS
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         2*scaled_re*base_flow_utheta*interpolated_VC*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*base_flow_uz*interpolated_dUdZC*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[1]*interpolated_dUdZC*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_WC*base_flow_durdz*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dZCdt*base_flow_durdz*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RC*base_flow_uz*base_flow_durdz*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RC*mesh_velocity[1]*base_flow_durdz*testf_
         *Jbar*w;
        residuals[local_eqn] += interpolated_RC*body_force[0]*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_inv_fr*interpolated_RC*G[0]*testf_*Jbar*w;
        residuals[local_eqn] += r*base_flow_p*dtestfdRC*Jbar*w;
        residuals[local_eqn] += r*interpolated_PC*dtestfdr*Jbar*w;
        residuals[local_eqn] += interpolated_RC*base_flow_p*dtestfdr*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr*dtestfdRC*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*r*interpolated_dUdRC*dtestfdr*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr*JhatC*dtestfdr*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*interpolated_RC*base_flow_durdr*dtestfdr
         *Jbar*w;
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[0]*base_flow_duthetadr*dtestfdr*interpolated_RS
         *Jbar*w;
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[0]*base_flow_duthetadr*dtestfdz*interpolated_ZS
         *Jbar*w;
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[0]*interpolated_dVdRS*testf_*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*k*k*interpolated_UC*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*k*base_flow_durdr*interpolated_RC*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*k*base_flow_durdz*interpolated_ZC*testf_*Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*k*base_flow_utheta*dtestfdr*interpolated_RS*Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*k*base_flow_utheta*dtestfdz*interpolated_ZS*Jbar*w/r;
        residuals[local_eqn] -= visc_ratio*k*interpolated_VS*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*interpolated_RS*base_flow_utheta*testf_*Jbar*w/(r*r);
        residuals[local_eqn] -=
         visc_ratio*Gamma[0]*r*base_flow_duzdr*dtestfdZC*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*Gamma[0]*r*interpolated_dWdRC*dtestfdz*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*Gamma[0]*r*base_flow_duzdr*JhatC*dtestfdz*w;
        residuals[local_eqn] -=
         visc_ratio*Gamma[0]*interpolated_RC*base_flow_duzdr*dtestfdz*Jbar*w;
        residuals[local_eqn] -= visc_ratio*r*base_flow_durdz*dtestfdZC*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*interpolated_dUdZC*dtestfdz*Jbar*w;
        residuals[local_eqn] += visc_ratio*r*base_flow_durdz*JhatC*dtestfdz*w;
        residuals[local_eqn] -=
         visc_ratio*interpolated_RC*base_flow_durdz*dtestfdz*Jbar*w;
        residuals[local_eqn] += interpolated_PC*testf_*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*k*interpolated_VS*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*k*base_flow_duthetadr*interpolated_RS
         *testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*k*base_flow_duthetadz*interpolated_ZS
         *testf_*Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*interpolated_UC*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*base_flow_ur*interpolated_RC
         *testf_*Jbar*w/(r*r);
        residuals[local_eqn] -= scaled_re_st*r*base_flow_durdt*testf_*JhatC*w;
        residuals[local_eqn] +=
         scaled_re*base_flow_utheta*base_flow_utheta*testf_*JhatC*w;
        residuals[local_eqn] += r*body_force[0]*testf_*JhatC*w;
        residuals[local_eqn] += scaled_re_inv_fr*r*G[0]*testf_*JhatC*w;
        residuals[local_eqn] -= k*visc_ratio*base_flow_utheta*testf_*JhatS*w/r;
        residuals[local_eqn] += base_flow_p*testf_*JhatC*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_*JhatC*w/r;
     
        // Calculate the Jacobian
        // ----------------------

        if(flag)
         {
          // Loop over the velocity shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Radial velocity component (cosine part) U_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[0]);
            if(local_unknown >= 0)
             {
              if(flag==2)
               {
                // Add the mass matrix
                mass_matrix(local_eqn,local_unknown) +=
                 scaled_re_st*r*psif[l2]*testf_*Jbar*w;
               }
              
              // Add contributions to the Jacobian matrix
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->time_stepper_pt()->weight(1,0)*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_ur*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[0]*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_durdr*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_uz*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[1]*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*r*group_B[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*k*k*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*group_A[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*psif[l2]*testf_*Jbar*w/r;
             }
            
            // Radial velocity component (sine part) U_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[1]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*scaled_re*base_flow_utheta*psif[l2]*testf_*Jbar*w;
             }

            // Axial velocity component (cosine part) W_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[2]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_durdz*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*Gamma[0]*r*group_B[l2]*dtestfdz*w;
             }

            // Axial velocity component (sine part) W_k^S
            // has no contribution

            // Azimuthal velocity component (cosine part) V_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[4]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               2.0*scaled_re*base_flow_utheta*psif[l2]*testf_*Jbar*w;
             }

            // Azimuthal velocity component (sine part) V_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[5]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               k*visc_ratio*Gamma[0]*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*k*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*k*psif[l2]*testf_*Jbar*w/r;
             }
            
            // Perturbation to radial nodal coord (cosine part) R_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[0]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*psif[l2]*base_flow_durdt*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*r*psif[l2]
              *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
              *base_flow_durdr*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re*psif[l2]*base_flow_ur*base_flow_durdr*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*psif[l2]*mesh_velocity[0]*base_flow_durdr
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re*r*base_flow_uz*group_C[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*r*mesh_velocity[1]*group_C[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re*psif[l2]*base_flow_uz*base_flow_durdz*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*psif[l2]*mesh_velocity[1]*base_flow_durdz
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              psif[l2]*body_force[0]*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_inv_fr*psif[l2]*G[0]*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              psif[l2]*base_flow_p*dtestfdr*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr*group_B[l2]
              *dtestfdr*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*(1.0+Gamma[0])*psif[l2]*base_flow_durdr
              *dtestfdr*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*k*k*base_flow_durdr*psif[l2]*testf_*Jbar*w/r;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*Gamma[0]*r*base_flow_duzdr*group_F(l,l2)*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*Gamma[0]*r*base_flow_duzdr*group_B[l2]*dtestfdz*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*Gamma[0]*psif[l2]*base_flow_duzdr*dtestfdz*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*base_flow_durdz*group_F(l,l2)*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*group_C[l2]*dtestfdz*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*base_flow_durdz*group_B[l2]*dtestfdz*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*psif[l2]*base_flow_durdz*dtestfdz*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[0])*base_flow_ur*psif[l2]*testf_
              *Jbar*w/(r*r);
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*r*base_flow_durdt*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re*base_flow_utheta*base_flow_utheta*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              r*body_force[0]*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_inv_fr*r*G[0]*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              base_flow_p*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_*group_B[l2]*w/r;
            }

            // Perturbation to radial nodal coord (sine part) R_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[1]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              k*scaled_re*base_flow_utheta*base_flow_durdr*psif[l2]
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              k*visc_ratio*Gamma[0]*base_flow_duthetadr*dtestfdr*psif[l2]
              *Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*k*base_flow_utheta*dtestfdr*psif[l2]*Jbar*w/r;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*k*psif[l2]*base_flow_utheta*testf_*Jbar*w/(r*r);
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[0])*k*base_flow_duthetadr*psif[l2]
              *testf_*Jbar*w/r;
             jacobian(local_eqn,local_unknown) -=
              k*visc_ratio*base_flow_utheta*testf_*group_B[l2]*w/r;
            }
            
            // Perturbation to axial nodal coord (cosine part) Z_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[2]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_ur*group_C[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[0]*group_C[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
               *base_flow_durdz*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               r*base_flow_p*group_F(l,l2)*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr*group_F(l,l2)*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*r*group_C[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr*group_A[l2]
               *dtestfdr*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*k*k*base_flow_durdz*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*Gamma[0]*r*group_D[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*Gamma[0]*r*base_flow_duzdr*group_A[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*base_flow_durdz*group_A[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*base_flow_durdt*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re*base_flow_utheta*base_flow_utheta*testf_
               *group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               r*body_force[0]*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_inv_fr*r*G[0]*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               base_flow_p*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_*group_A[l2]*w/r;
             }
            
            // Perturbation to axial nodal coord (sine part) Z_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[3]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               k*scaled_re*base_flow_utheta*base_flow_durdz*psif[l2]
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               k*visc_ratio*Gamma[0]*base_flow_duthetadr*dtestfdz*psif[l2]
               *Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               k*visc_ratio*Gamma[0]*group_E[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*k*base_flow_utheta*dtestfdz*psif[l2]*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*(1.0+Gamma[0])*k*base_flow_duthetadz*psif[l2]
               *testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*base_flow_utheta*testf_*group_A[l2]*w/r;
             }

           } // End of loop over velocity shape functions
          
          // Now loop over pressure shape functions
          // (This is the contribution from pressure gradient)
          for(unsigned l2=0;l2<n_pres;l2++)
           {
            // Cosine part P_k^C
            local_unknown = p_local_eqn(l2,0);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) += r*psip[l2]*dtestfdr*Jbar*w;
              jacobian(local_eqn,local_unknown) += psip[l2]*testf_*Jbar*w;
             }

            // Sine part P_k^S has no contribution

           } // End of loop over pressure shape functions

          // Now loop over the spines in the element
          for(unsigned l2=0;l2<n_p;l2++)
           {
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown =
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,0);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum -= scaled_re*r*base_flow_ur*group_C[l2+(3*a)]*testf_*w;
                sum += scaled_re_st*r*mesh_velocity[0]*group_C[l2+(3*a)]
                 *testf_*w;
                sum += scaled_re_st*r*psif[l2+(3*a)]
                 *node_pt(l2+(3*a))->position_time_stepper_pt()->weight(1,0)
                 *base_flow_durdz*testf_*Jbar*w;
                sum += r*base_flow_p*group_F(l,(l2+(3*a)))*w;
                sum -= visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr
                 *group_F(l,(l2+(3*a)))*w;
                sum -= visc_ratio*(1.0+Gamma[0])*r*group_C[l2+(3*a)]
                 *dtestfdr*w;
                sum += visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr
                 *group_A[l2+(3*a)]*dtestfdr*w;
                sum += visc_ratio*k*k*base_flow_durdz*psif[l2+(3*a)]
                 *testf_*Jbar*w/r;
                sum -= visc_ratio*Gamma[0]*r*group_D[l2+(3*a)]*dtestfdz*w;
                sum += visc_ratio*Gamma[0]*r*base_flow_duzdr
                 *group_A[l2+(3*a)]*dtestfdz*w;
                sum += visc_ratio*r*base_flow_durdz*group_A[l2+(3*a)]
                 *dtestfdz*w;
                sum -= scaled_re_st*r*base_flow_durdt*testf_
                 *group_A[l2+(3*a)]*w;
                sum += scaled_re*base_flow_utheta*base_flow_utheta
                 *testf_*group_A[l2+(3*a)]*w;
                sum += r*body_force[0]*testf_*group_A[l2+(3*a)]*w;
                sum += scaled_re_inv_fr*r*G[0]*testf_*group_A[l2+(3*a)]*w;
                sum += base_flow_p*testf_*group_A[l2+(3*a)]*w;
                sum -= visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_
                 *group_A[l2+(3*a)]*w/r;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
            
            // Perturbed spine "height" (sine part) H_k^S
            local_unknown =
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,1);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum += k*scaled_re*base_flow_utheta*base_flow_durdz
                 *psif[l2+(3*a)]*testf_*Jbar*w;
                sum += k*visc_ratio*Gamma[0]*base_flow_duthetadr*dtestfdz
                 *psif[l2+(3*a)]*Jbar*w;
                sum += k*visc_ratio*Gamma[0]*group_E[l2+(3*a)]*testf_*w;
                sum -= visc_ratio*k*base_flow_utheta*dtestfdz
                 *psif[l2+(3*a)]*Jbar*w/r;
                sum += visc_ratio*(1.0+Gamma[0])*k*base_flow_duthetadz
                 *psif[l2+(3*a)]*testf_*Jbar*w/r;
                sum -= k*visc_ratio*base_flow_utheta*testf_
                 *group_A[l2+(3*a)]*w/r;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
           } // End of loop over spines in the element
         } // End of Jacobian calculation
        
       } // End of if not boundary condition statement
      
      // --------------------------------------------
      // SECOND (RADIAL) MOMENTUM EQUATION: SINE PART
      // --------------------------------------------
      
      // Get local equation number of second velocity value at this node
      local_eqn = nodal_local_eqn(l,u_nodal_index[1]);

      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] -= scaled_re_st*r*dUSdt*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re_st*interpolated_RS*base_flow_durdt*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*base_flow_ur*interpolated_dUdRS*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[0]*interpolated_dUdRS*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_US*base_flow_durdr*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dRSdt*base_flow_durdr*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RS*base_flow_ur*base_flow_durdr
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RS*mesh_velocity[0]*base_flow_durdr
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         k*scaled_re*base_flow_utheta*interpolated_UC*testf_*Jbar*w;
        residuals[local_eqn] -=
         k*scaled_re*base_flow_utheta*base_flow_durdr*interpolated_RC
         *testf_*Jbar*w;
        residuals[local_eqn] -=
         k*scaled_re*base_flow_utheta*base_flow_durdz*interpolated_ZC
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         2*scaled_re*base_flow_utheta*interpolated_VS*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*base_flow_uz*interpolated_dUdZS*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[1]*interpolated_dUdZS*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_WS*base_flow_durdz*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dZSdt*base_flow_durdz*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RS*base_flow_uz*base_flow_durdz*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RS*mesh_velocity[1]*base_flow_durdz
         *testf_*Jbar*w;
        residuals[local_eqn] += interpolated_RS*body_force[0]*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_inv_fr*interpolated_RS*G[0]*testf_*Jbar*w;
        residuals[local_eqn] += r*base_flow_p*dtestfdRS*Jbar*w;
        residuals[local_eqn] += r*interpolated_PS*dtestfdr*Jbar*w;
        residuals[local_eqn] += interpolated_RS*base_flow_p*dtestfdr*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr*dtestfdRS*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*r*interpolated_dUdRS*dtestfdr*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr*JhatS*dtestfdr*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*interpolated_RS*base_flow_durdr*dtestfdr
         *Jbar*w;
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[0]*base_flow_duthetadr*dtestfdr*interpolated_RC
         *Jbar*w;
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[0]*base_flow_duthetadr*dtestfdz*interpolated_ZC
         *Jbar*w;
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[0]*interpolated_dVdRC*testf_*Jbar*w;
        residuals[local_eqn] -= visc_ratio*k*k*interpolated_US*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*k*base_flow_durdr*interpolated_RS*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*k*base_flow_durdz*interpolated_ZS*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*base_flow_utheta*dtestfdr*interpolated_RC*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*base_flow_utheta*dtestfdz*interpolated_ZC*Jbar*w/r;
        residuals[local_eqn] += visc_ratio*k*interpolated_VC*testf_*Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*k*interpolated_RC*base_flow_utheta*testf_*Jbar*w/(r*r);
        residuals[local_eqn] -=
         visc_ratio*Gamma[0]*r*base_flow_duzdr*dtestfdZS*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*Gamma[0]*r*interpolated_dWdRS*dtestfdz*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*Gamma[0]*r*base_flow_duzdr*JhatS*dtestfdz*w;
        residuals[local_eqn] -=
         visc_ratio*Gamma[0]*interpolated_RS*base_flow_duzdr*dtestfdz*Jbar*w;
        residuals[local_eqn] -= visc_ratio*r*base_flow_durdz*dtestfdZS*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*interpolated_dUdZS*dtestfdz*Jbar*w;
        residuals[local_eqn] += visc_ratio*r*base_flow_durdz*JhatS*dtestfdz*w;
        residuals[local_eqn] -=
         visc_ratio*interpolated_RS*base_flow_durdz*dtestfdz*Jbar*w;
        residuals[local_eqn] += interpolated_PS*testf_*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*k*interpolated_VC*testf_*Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*k*base_flow_duthetadr*interpolated_RC
         *testf_*Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*k*base_flow_duthetadz*interpolated_ZC
         *testf_*Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*interpolated_US*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*base_flow_ur*interpolated_RS*testf_
         *Jbar*w/(r*r);
        residuals[local_eqn] -= scaled_re_st*r*base_flow_durdt*testf_*JhatS*w;
        residuals[local_eqn] +=
         scaled_re*base_flow_utheta*base_flow_utheta*testf_*JhatS*w;
        residuals[local_eqn] += r*body_force[0]*testf_*JhatS*w;
        residuals[local_eqn] += scaled_re_inv_fr*r*G[0]*testf_*JhatS*w;
        residuals[local_eqn] += k*visc_ratio*base_flow_utheta*testf_*JhatC*w/r;
        residuals[local_eqn] += base_flow_p*testf_*JhatS*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_*JhatS*w/r;
        
        
        // Calculate the Jacobian
        // ----------------------

        if(flag)
         {
          // Loop over the velocity shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Radial velocity component (cosine part) U_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[0]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               k*scaled_re*base_flow_utheta*psif[l2]*testf_*Jbar*w;
             }

            // Radial velocity component (sine part) U_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[1]);
            if(local_unknown >= 0)
             {
              if(flag==2)
               {
                // Add the mass matrix
                mass_matrix(local_eqn,local_unknown) +=
                 scaled_re_st*r*psif[l2]*testf_*Jbar*w;
               }
              
              // Add contributions to the Jacobian matrix
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->time_stepper_pt()->weight(1,0)*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_ur*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[0]*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_durdr*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_uz*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[1]*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*r*group_B[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*k*k*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*group_A[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*psif[l2]*testf_*Jbar*w/r;
             }

            // Axial velocity component (cosine part) W_k^C
            // has no contribution

            // Axial velocity component (sine part) W_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[3]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_durdz*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*Gamma[0]*r*group_B[l2]*dtestfdz*w;
             }

            // Azimuthal velocity component (cosine part) V_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[4]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[0]*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*k*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*(1.0+Gamma[0])*k*psif[l2]*testf_*Jbar*w/r;
             }
            
            // Azimuthal velocity component (sine part) V_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[5]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               2.0*scaled_re*base_flow_utheta*psif[l2]*testf_*Jbar*w;
             }
            
            // Perturbation to radial nodal coord (cosine part) R_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[0]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*scaled_re*base_flow_utheta*base_flow_durdr*psif[l2]
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[0]*base_flow_duthetadr*dtestfdr
               *psif[l2]*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*k*base_flow_utheta*dtestfdr*psif[l2]*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*k*psif[l2]*base_flow_utheta*testf_*Jbar*w/(r*r);
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*k*base_flow_duthetadr*psif[l2]
               *testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               k*visc_ratio*base_flow_utheta*testf_*group_B[l2]*w/r;
             }
            
            // Perturbation to radial nodal coord (sine part) R_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[1]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*psif[l2]*base_flow_durdt*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
               *base_flow_durdr*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*psif[l2]*base_flow_ur*base_flow_durdr*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*psif[l2]*mesh_velocity[0]*base_flow_durdr
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re*r*base_flow_uz*group_C[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*mesh_velocity[1]*group_C[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*psif[l2]*base_flow_uz*base_flow_durdz*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*psif[l2]*mesh_velocity[1]*base_flow_durdz
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               psif[l2]*body_force[0]*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_inv_fr*psif[l2]*G[0]*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               psif[l2]*base_flow_p*dtestfdr*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr*group_B[l2]
               *dtestfdr*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*psif[l2]*base_flow_durdr
               *dtestfdr*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*k*k*base_flow_durdr*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*Gamma[0]*r*base_flow_duzdr*group_F(l,l2)*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*Gamma[0]*r*base_flow_duzdr*group_B[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*Gamma[0]*psif[l2]*base_flow_duzdr*dtestfdz*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*base_flow_durdz*group_F(l,l2)*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*group_C[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*base_flow_durdz*group_B[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*psif[l2]*base_flow_durdz*dtestfdz*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*(1.0+Gamma[0])*base_flow_ur*psif[l2]
               *testf_*Jbar*w/(r*r);
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*base_flow_durdt*testf_*group_B[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re*base_flow_utheta*base_flow_utheta
               *testf_*group_B[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               r*body_force[0]*testf_*group_B[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_inv_fr*r*G[0]*testf_*group_B[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               base_flow_p*testf_*group_B[l2]*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_*group_B[l2]*w/r;
             }
            
            // Perturbation to axial nodal coord (cosine part) Z_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[2]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*scaled_re*base_flow_utheta*base_flow_durdz*psif[l2]
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[0]*base_flow_duthetadr*dtestfdz
               *psif[l2]*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[0]*group_E[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*k*base_flow_utheta*dtestfdz*psif[l2]*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*k*base_flow_duthetadz*psif[l2]
               *testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               k*visc_ratio*base_flow_utheta*testf_*group_A[l2]*w/r;
             }
            
            // Perturbation to axial nodal coord (sine part) Z_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[3]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_ur*group_C[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[0]*group_C[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
               *base_flow_durdz*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               r*base_flow_p*group_F(l,l2)*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr*group_F(l,l2)*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*r*group_C[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr*group_A[l2]
               *dtestfdr*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*k*k*base_flow_durdz*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*Gamma[0]*r*base_flow_duzdr*dtestfdz*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*Gamma[0]*r*group_D[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*base_flow_durdz*group_A[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*base_flow_durdt*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re*base_flow_utheta*base_flow_utheta*testf_
               *group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               r*body_force[0]*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_inv_fr*r*G[0]*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               base_flow_p*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_*group_A[l2]*w/r;
             }
            
           } // End of loop over velocity shape functions
          
          // Now loop over pressure shape functions
          // (This is the contribution from pressure gradient)
          for(unsigned l2=0;l2<n_pres;l2++)
           {
            // Cosine part P_k^C has no contribution
            
            // Sine part P_k^S
            local_unknown = p_local_eqn(l2,1);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) += r*psip[l2]*dtestfdr*Jbar*w;
              jacobian(local_eqn,local_unknown) += psip[l2]*testf_*Jbar*w;
             }
           } // End of loop over pressure shape functions

          // Now loop over the spines in the element
          for(unsigned l2=0;l2<n_p;l2++)
           {
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown =
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,0);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum -=
                 k*scaled_re*base_flow_utheta*base_flow_durdz
                 *psif[l2+(3*a)]*testf_*Jbar*w;
                sum -=
                 k*visc_ratio*Gamma[0]*base_flow_duthetadr*dtestfdz
                 *psif[l2+(3*a)]*Jbar*w;
                sum -= k*visc_ratio*Gamma[0]*group_E[l2+(3*a)]*testf_*w;
                sum +=
                 visc_ratio*k*base_flow_utheta*dtestfdz*psif[l2+(3*a)]
                 *Jbar*w/r;
                sum -=
                 visc_ratio*(1.0+Gamma[0])*k*base_flow_duthetadz
                 *psif[l2+(3*a)]*testf_*Jbar*w/r;
                sum +=
                 k*visc_ratio*base_flow_utheta*testf_*group_A[l2+(3*a)]*w/r;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
            
            // Perturbed spine "height" (sine part) H_k^S
            local_unknown =
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,1);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum -= scaled_re*r*base_flow_ur*group_C[l2+(3*a)]*testf_*w;
                sum +=
                 scaled_re_st*r*mesh_velocity[0]*group_C[l2+(3*a)]*testf_*w;
                sum +=
                 scaled_re_st*r*psif[l2+(3*a)]
                 *node_pt(l2+(3*a))->position_time_stepper_pt()->weight(1,0)
                 *base_flow_durdz*testf_*Jbar*w;
                sum += r*base_flow_p*group_F(l,(l2+(3*a)))*w;
                sum -=
                 visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr
                 *group_F(l,(l2+(3*a)))*w;
                sum -=
                 visc_ratio*(1.0+Gamma[0])*r*group_C[l2+(3*a)]*dtestfdr*w;
                sum +=
                 visc_ratio*(1.0+Gamma[0])*r*base_flow_durdr
                 *group_A[l2+(3*a)]*dtestfdr*w;
                sum +=
                 visc_ratio*k*k*base_flow_durdz*psif[l2+(3*a)]*testf_*Jbar*w/r;
                sum -= visc_ratio*Gamma[0]*r*group_D[l2+(3*a)]*dtestfdz*w;
                sum +=
                 visc_ratio*Gamma[0]*r*base_flow_duzdr*group_A[l2+(3*a)]
                 *dtestfdz*w;
                sum +=
                 visc_ratio*r*base_flow_durdz*group_A[l2+(3*a)]*dtestfdz*w;
                sum -=
                 scaled_re_st*r*base_flow_durdt*testf_*group_A[l2+(3*a)]*w;
                sum +=
                 scaled_re*base_flow_utheta*base_flow_utheta*testf_
                 *group_A[l2+(3*a)]*w;
                sum += r*body_force[0]*testf_*group_A[l2+(3*a)]*w;
                sum += scaled_re_inv_fr*r*G[0]*testf_*group_A[l2+(3*a)]*w;
                sum += base_flow_p*testf_*group_A[l2+(3*a)]*w;
                sum -=
                 visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_
                 *group_A[l2+(3*a)]*w/r;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
           } // End of loop over spines in the element
         } // End of Jacobian calculation
        
       } // End of if not boundary condition statement

      // --------------------------------------------
      // THIRD (AXIAL) MOMENTUM EQUATION: COSINE PART
      // --------------------------------------------

      // Get local equation number of third velocity value at this node
      local_eqn = nodal_local_eqn(l,u_nodal_index[2]);

      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] -= scaled_re_st*r*dWCdt*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re_st*interpolated_RC*base_flow_duzdt*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*base_flow_ur*interpolated_dWdRC*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[0]*interpolated_dWdRC*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_UC*base_flow_duzdr*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dRCdt*base_flow_duzdr*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RC*base_flow_ur*base_flow_duzdr*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RC*mesh_velocity[0]*base_flow_duzdr
         *testf_*Jbar*w;
        residuals[local_eqn] -=
         k*scaled_re*base_flow_utheta*interpolated_WS*testf_*Jbar*w;
        residuals[local_eqn] +=
         k*scaled_re*base_flow_utheta*base_flow_duzdr*interpolated_RS
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         k*scaled_re*base_flow_utheta*base_flow_duzdz*interpolated_ZS
         *testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*base_flow_uz*interpolated_dWdZC*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[1]*interpolated_dWdZC*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_WC*base_flow_duzdz*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dZCdt*base_flow_duzdz*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RC*base_flow_uz*base_flow_duzdz*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RC*mesh_velocity[1]*base_flow_duzdz
         *testf_*Jbar*w;
        residuals[local_eqn] += interpolated_RC*body_force[1]*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_inv_fr*interpolated_RC*G[1]*testf_*Jbar*w;
        residuals[local_eqn] -= visc_ratio*r*base_flow_duzdr*dtestfdRC*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*interpolated_dWdRC*dtestfdr*Jbar*w;
        residuals[local_eqn] += visc_ratio*r*base_flow_duzdr*JhatC*dtestfdr*w;
        residuals[local_eqn] -=
         visc_ratio*interpolated_RC*base_flow_duzdr*dtestfdr*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*Gamma[1]*r*base_flow_durdz*dtestfdRC*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*Gamma[1]*r*interpolated_dUdZC*dtestfdr*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*Gamma[1]*r*base_flow_durdz*JhatC*dtestfdr*w;
        residuals[local_eqn] -=
         visc_ratio*Gamma[1]*interpolated_RC*base_flow_durdz*dtestfdr*Jbar*w;
        residuals[local_eqn] -= visc_ratio*k*k*interpolated_WC*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*k*base_flow_duzdr*interpolated_RC*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*k*base_flow_duzdz*interpolated_ZC*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[1]*base_flow_duthetadz*dtestfdr*interpolated_RS
         *Jbar*w;
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[1]*base_flow_duthetadz*dtestfdz*interpolated_ZS
         *Jbar*w;
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[1]*interpolated_dVdZS*testf_*Jbar*w;
        residuals[local_eqn] += r*base_flow_p*dtestfdZC*Jbar*w;
        residuals[local_eqn] += r*interpolated_PC*dtestfdz*Jbar*w;
        residuals[local_eqn] += interpolated_RC*base_flow_p*dtestfdz*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz*dtestfdZC*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[1])*r*interpolated_dWdZC*dtestfdz*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz*JhatC*dtestfdz*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[1])*interpolated_RC*base_flow_duzdz
         *dtestfdz*Jbar*w;
        residuals[local_eqn] -= scaled_re_st*r*base_flow_duzdt*testf_*JhatC*w; 
        residuals[local_eqn] += r*body_force[1]*testf_*JhatC*w;
        residuals[local_eqn] += scaled_re_inv_fr*r*G[1]*testf_*JhatC*w;
        
        
        // Calculate the Jacobian
        // ----------------------

        if(flag)
         {
          // Loop over the velocity shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Radial velocity component (cosine part) U_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[0]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_duzdr*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*Gamma[1]*r*group_A[l2]*dtestfdr*w;
             }

            // Radial velocity component (sine part) U_k^S
            // has no contribution
            
            // Axial velocity component (cosine part) W_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[2]);
            if(local_unknown >= 0)
             {
              if(flag==2)
               {
                // Add the mass matrix
                mass_matrix(local_eqn,local_unknown) +=
                 scaled_re_st*r*psif[l2]*testf_*Jbar*w;
               }
              
              // Add contributions to the Jacobian matrix
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->time_stepper_pt()->weight(1,0)*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_ur*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[0]*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_uz*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[1]*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_duzdz*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*group_B[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*k*k*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[1])*r*group_A[l2]*dtestfdz*w;
             }
            
            // Axial velocity component (sine part) W_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[3]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*scaled_re*base_flow_utheta*psif[l2]*testf_*Jbar*w;
             }
           
            // Azimuthal velocity component (cosine part) V_k^C
            // has no contribution

            // Azimuthal velocity component (sine part) V_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[5]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               k*visc_ratio*Gamma[1]*group_A[l2]*testf_*w;
             }

            // Perturbation to radial nodal coord (cosine part) R_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[0]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*psif[l2]*base_flow_duzdt*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*r*psif[l2]
              *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
              *base_flow_duzdr*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re*psif[l2]*base_flow_ur*base_flow_duzdr*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*psif[l2]*mesh_velocity[0]*base_flow_duzdr
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re*r*base_flow_uz*group_D[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*r*mesh_velocity[1]*group_D[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re*psif[l2]*base_flow_uz*base_flow_duzdz*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*psif[l2]*mesh_velocity[1]*base_flow_duzdz
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              psif[l2]*body_force[1]*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_inv_fr*psif[l2]*G[1]*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*base_flow_duzdr*group_B[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*psif[l2]*base_flow_duzdr*dtestfdr*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*Gamma[1]*r*group_C[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*Gamma[1]*r*base_flow_durdz*group_B[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*Gamma[1]*psif[l2]*base_flow_durdz*dtestfdr*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*k*k*base_flow_duzdr*psif[l2]*testf_*Jbar*w/r;
             jacobian(local_eqn,local_unknown) -=
              r*base_flow_p*group_F(l,l2)*w;
             jacobian(local_eqn,local_unknown) +=
              psif[l2]*base_flow_p*dtestfdz*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz*group_F(l,l2)*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[1])*r*group_D[l2]*dtestfdz*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz*group_B[l2]
              *dtestfdz*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*(1.0+Gamma[1])*psif[l2]*base_flow_duzdz*dtestfdz
              *Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*r*base_flow_duzdt*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              r*body_force[1]*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_inv_fr*r*G[1]*testf_*group_B[l2]*w;
            }
            
            // Perturbation to radial nodal coord (sine part) R_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[1]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              k*scaled_re*base_flow_utheta*base_flow_duzdr*psif[l2]*testf_
              *Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              k*visc_ratio*Gamma[1]*base_flow_duthetadz*dtestfdr*psif[l2]
              *Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              k*visc_ratio*Gamma[1]*group_E[l2]*testf_*w;
            }

            // Perturbation to axial nodal coord (cosine part) Z_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[2]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=
              scaled_re*r*base_flow_ur*group_D[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*r*mesh_velocity[0]*group_D[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*r*psif[l2]
              *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
              *base_flow_duzdz*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*r*base_flow_duzdr*group_F(l,l2)*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*r*group_D[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*base_flow_duzdr*group_A[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*Gamma[1]*r*base_flow_durdz*group_F(l,l2)*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*Gamma[1]*r*base_flow_durdz*group_A[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*k*k*base_flow_duzdz*psif[l2]*testf_*Jbar*w/r;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz*group_A[l2]
              *dtestfdz*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*r*base_flow_duzdt*testf_*group_A[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              r*body_force[1]*testf_*group_A[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_inv_fr*r*G[1]*testf_*group_A[l2]*w;
            }
            
            // Perturbation to axial nodal coord (sine part) Z_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[3]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               k*scaled_re*base_flow_utheta*base_flow_duzdz*psif[l2]
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               k*visc_ratio*Gamma[1]*base_flow_duthetadz*dtestfdz
               *psif[l2]*Jbar*w;
             }
            
           } // End of loop over velocity shape functions
          
          // Now loop over pressure shape functions
          // (This is the contribution from pressure gradient)
          for(unsigned l2=0;l2<n_pres;l2++)
           {
            // Cosine part P_k^C
            local_unknown = p_local_eqn(l2,0);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) += r*psip[l2]*dtestfdz*Jbar*w;
             }

            // Sine part P_k^S has no contribution

           } // End of loop over pressure shape functions

          // Now loop over the spines in the element
          for(unsigned l2=0;l2<n_p;l2++)
           {
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown = 
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,0);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum -= scaled_re*r*base_flow_ur*group_D[l2+(3*a)]*testf_*w;
                sum +=
                 scaled_re_st*r*mesh_velocity[0]*group_D[l2+(3*a)]*testf_*w;
                sum +=
                 scaled_re_st*r*psif[l2+(3*a)]
                 *node_pt(l2+(3*a))->position_time_stepper_pt()->weight(1,0)
                 *base_flow_duzdz*testf_*Jbar*w;
                sum -= visc_ratio*r*base_flow_duzdr*group_F(l,(l2+(3*a)))*w;
                sum -= visc_ratio*r*group_D[l2+(3*a)]*dtestfdr*w;
                sum +=
                 visc_ratio*r*base_flow_duzdr*group_A[l2+(3*a)]*dtestfdr*w;
                sum -=
                 visc_ratio*Gamma[1]*r*base_flow_durdz*group_F(l,(l2+(3*a)))*w;
                sum +=
                 visc_ratio*Gamma[1]*r*base_flow_durdz*group_A[l2+(3*a)]
                 *dtestfdr*w;
                sum +=
                 visc_ratio*k*k*base_flow_duzdz*psif[l2+(3*a)]*testf_*Jbar*w/r;
                sum +=
                 visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz
                 *group_A[l2+(3*a)]*dtestfdz*w;
                sum -=
                 scaled_re_st*r*base_flow_duzdt*testf_*group_A[l2+(3*a)]*w;
                sum += r*body_force[1]*testf_*group_A[l2+(3*a)]*w;
                sum += scaled_re_inv_fr*r*G[1]*testf_*group_A[l2+(3*a)]*w;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
            
            // Perturbed spine "height" (sine part) H_k^S
            local_unknown = 
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,1);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum += k*scaled_re*base_flow_utheta*base_flow_duzdz
                 *psif[l2+(3*a)]*testf_*Jbar*w;
                sum += k*visc_ratio*Gamma[1]*base_flow_duthetadz*dtestfdz
                 *psif[l2+(3*a)]*Jbar*w;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
           } // End of loop over spines in the element

         } // End of Jacobian calculation
        
       } // End of if not boundary condition statement

      // -------------------------------------------
      // FOURTH (AXIAL) MOMENTUM EQUATION: SINE PART
      // -------------------------------------------

      // Get local equation number of fourth velocity value at this node
      local_eqn = nodal_local_eqn(l,u_nodal_index[3]);

      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] -= scaled_re_st*r*dWSdt*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re_st*interpolated_RS*base_flow_duzdt*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*base_flow_ur*interpolated_dWdRS*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[0]*interpolated_dWdRS*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_US*base_flow_duzdr*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dRSdt*base_flow_duzdr*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RS*base_flow_ur*base_flow_duzdr*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RS*mesh_velocity[0]*base_flow_duzdr
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         k*scaled_re*base_flow_utheta*interpolated_WC*testf_*Jbar*w;
        residuals[local_eqn] -=
         k*scaled_re*base_flow_utheta*base_flow_duzdr*interpolated_RC
         *testf_*Jbar*w;
        residuals[local_eqn] -=
         k*scaled_re*base_flow_utheta*base_flow_duzdz*interpolated_ZC
         *testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*base_flow_uz*interpolated_dWdZS*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[1]*interpolated_dWdZS*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_WS*base_flow_duzdz*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dZSdt*base_flow_duzdz*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RS*base_flow_uz*base_flow_duzdz*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RS*mesh_velocity[1]*base_flow_duzdz
         *testf_*Jbar*w;
        residuals[local_eqn] += interpolated_RS*body_force[1]*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_inv_fr*interpolated_RS*G[1]*testf_*Jbar*w;
        residuals[local_eqn] -= visc_ratio*r*base_flow_duzdr*dtestfdRS*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*interpolated_dWdRS*dtestfdr*Jbar*w;
        residuals[local_eqn] += visc_ratio*r*base_flow_duzdr*JhatS*dtestfdr*w;
        residuals[local_eqn] -=
         visc_ratio*interpolated_RS*base_flow_duzdr*dtestfdr*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*Gamma[1]*r*base_flow_durdz*dtestfdRS*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*Gamma[1]*r*interpolated_dUdZS*dtestfdr*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*Gamma[1]*r*base_flow_durdz*JhatS*dtestfdr*w;
        residuals[local_eqn] -=
         visc_ratio*Gamma[1]*interpolated_RS*base_flow_durdz*dtestfdr*Jbar*w;
        residuals[local_eqn] -= visc_ratio*k*k*interpolated_WS*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*k*base_flow_duzdr*interpolated_RS*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*k*base_flow_duzdz*interpolated_ZS*testf_*Jbar*w/r;
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[1]*base_flow_duthetadz*dtestfdr*interpolated_RC
         *Jbar*w;
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[1]*base_flow_duthetadz*dtestfdz*interpolated_ZC
         *Jbar*w;
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[1]*interpolated_dVdZC*testf_*Jbar*w;
        residuals[local_eqn] += r*base_flow_p*dtestfdZS*Jbar*w;
        residuals[local_eqn] += r*interpolated_PS*dtestfdz*Jbar*w;
        residuals[local_eqn] += interpolated_RS*base_flow_p*dtestfdz*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz*dtestfdZS*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[1])*r*interpolated_dWdZS*dtestfdz*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz*JhatS*dtestfdz*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[1])*interpolated_RS*base_flow_duzdz*dtestfdz
         *Jbar*w;
        residuals[local_eqn] -= scaled_re_st*r*base_flow_duzdt*testf_*JhatS*w; 
        residuals[local_eqn] += r*body_force[1]*testf_*JhatS*w;
        residuals[local_eqn] += scaled_re_inv_fr*r*G[1]*testf_*JhatS*w;

          
        // Calculate the Jacobian
        // ----------------------

        if(flag)
         {
          // Loop over the velocity shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Radial velocity component (cosine part) U_k^C
            // has no contribution

            // Radial velocity component (sine part) U_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[1]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_duzdr*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*Gamma[1]*r*group_A[l2]*dtestfdr*w;
             }
            
            // Axial velocity component (cosine part) W_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[2]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               k*scaled_re*base_flow_utheta*psif[l2]*testf_*Jbar*w;
             }

            // Axial velocity component (sine part) W_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[3]);
            if(local_unknown >= 0)
             {
              if(flag==2)
               {
                // Add the mass matrix
                mass_matrix(local_eqn,local_unknown) +=
                 scaled_re_st*r*psif[l2]*testf_*Jbar*w;
               }

              // Add contributions to the Jacobian matrix
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->time_stepper_pt()->weight(1,0)*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_ur*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[0]*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_uz*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[1]*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_duzdz*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*group_B[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*k*k*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[1])*r*group_A[l2]*dtestfdz*w;
             }
            
            // Azimuthal velocity component (cosine part) V_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[4]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[1]*group_A[l2]*testf_*w;
             }
            
            // Azimuthal velocity component (sine part) V_k^S
            // has no contribution
            
            // Perturbation to radial nodal coord (cosine part) R_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[0]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*scaled_re*base_flow_utheta*base_flow_duzdr*psif[l2]
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[1]*base_flow_duthetadz*dtestfdr
               *psif[l2]*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               k*visc_ratio*Gamma[1]*group_E[l2]*testf_*w;
             }
            
            // Perturbation to radial nodal coord (sine part) R_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[1]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*psif[l2]*base_flow_duzdt*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*r*psif[l2]
              *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
              *base_flow_duzdr*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re*psif[l2]*base_flow_ur*base_flow_duzdr*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*psif[l2]*mesh_velocity[0]*base_flow_duzdr
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re*r*base_flow_uz*group_D[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*r*mesh_velocity[1]*group_D[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re*psif[l2]*base_flow_uz*base_flow_duzdz*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*psif[l2]*mesh_velocity[1]*base_flow_duzdz
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              psif[l2]*body_force[1]*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_inv_fr*psif[l2]*G[1]*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*base_flow_duzdr*group_B[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*psif[l2]*base_flow_duzdr*dtestfdr*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*Gamma[1]*r*group_C[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*Gamma[1]*r*base_flow_durdz*group_B[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*Gamma[1]*psif[l2]*base_flow_durdz*dtestfdr*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*k*k*base_flow_duzdr*psif[l2]*testf_*Jbar*w/r;
             jacobian(local_eqn,local_unknown) -=
              r*base_flow_p*group_F(l,l2)*w;
             jacobian(local_eqn,local_unknown) +=
              psif[l2]*base_flow_p*dtestfdz*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz*group_F(l,l2)*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[1])*r*group_D[l2]*dtestfdz*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz*group_B[l2]
              *dtestfdz*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*(1.0+Gamma[1])*psif[l2]*base_flow_duzdz
              *dtestfdz*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*r*base_flow_duzdt*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              r*body_force[1]*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_inv_fr*r*G[1]*testf_*group_B[l2]*w;
            }
            
            // Perturbation to axial nodal coord (cosine part) Z_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[2]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*scaled_re*base_flow_utheta*base_flow_duzdz*psif[l2]
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[1]*base_flow_duthetadz*dtestfdz
               *psif[l2]*Jbar*w;
             }
            
            // Perturbation to axial nodal coord (sine part) Z_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[3]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_ur*group_D[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[0]*group_D[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
               *base_flow_duzdz*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*base_flow_duzdr*group_F(l,l2)*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*group_D[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*base_flow_duzdr*group_A[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*Gamma[1]*r*base_flow_durdz*group_F(l,l2)*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*Gamma[1]*r*base_flow_durdz*group_A[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*k*k*base_flow_duzdz*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz*group_A[l2]
               *dtestfdz*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*base_flow_duzdt*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               r*body_force[1]*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_inv_fr*r*G[1]*testf_*group_A[l2]*w;
             }
            
           } // End of loop over velocity shape functions
          
          // Now loop over pressure shape functions
          // (This is the contribution from pressure gradient)
          for(unsigned l2=0;l2<n_pres;l2++)
           {
            // Cosine part P_k^C has no contribution
            
            // Sine part P_k^S
            local_unknown = p_local_eqn(l2,1);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) += r*psip[l2]*dtestfdz*Jbar*w;
             }
           } // End of loop over pressure shape functions

          // Now loop over the spines in the element
          for(unsigned l2=0;l2<n_p;l2++)
           {
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown = 
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,0);
 
           if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum -=
                 k*scaled_re*base_flow_utheta*base_flow_duzdz
                 *psif[l2+(3*a)]*testf_*Jbar*w;
                sum -=
                 k*visc_ratio*Gamma[1]*base_flow_duthetadz*dtestfdz
                 *psif[l2+(3*a)]*Jbar*w;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
            
            // Perturbed spine "height" (sine part) H_k^S
            local_unknown = 
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,1);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum -= scaled_re*r*base_flow_ur*group_D[l2+(3*a)]*testf_*w;
                sum +=
                 scaled_re_st*r*mesh_velocity[0]*group_D[l2+(3*a)]*testf_*w;
                sum +=
                 scaled_re_st*r*psif[l2+(3*a)]
                 *node_pt(l2+(3*a))->position_time_stepper_pt()->weight(1,0)
                 *base_flow_duzdz*testf_*Jbar*w;
                sum -= visc_ratio*r*base_flow_duzdr*group_F(l,(l2+(3*a)))*w;
                sum -= visc_ratio*r*group_D[l2+(3*a)]*dtestfdr*w;
                sum +=
                 visc_ratio*r*base_flow_duzdr*group_A[l2+(3*a)]*dtestfdr*w;
                sum -=
                 visc_ratio*Gamma[1]*r*base_flow_durdz*group_F(l,(l2+(3*a)))*w;
                sum +=
                 visc_ratio*Gamma[1]*r*base_flow_durdz*group_A[l2+(3*a)]
                 *dtestfdr*w;
                sum +=
                 visc_ratio*k*k*base_flow_duzdz*psif[l2+(3*a)]*testf_*Jbar*w/r;
                sum +=
                 visc_ratio*(1.0+Gamma[1])*r*base_flow_duzdz
                 *group_A[l2+(3*a)]*dtestfdz*w;
                sum -=
                 scaled_re_st*r*base_flow_duzdt*testf_*group_A[l2+(3*a)]*w;
                sum += r*body_force[1]*testf_*group_A[l2+(3*a)]*w;
                sum += scaled_re_inv_fr*r*G[1]*testf_*group_A[l2+(3*a)]*w;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
           } // End of loop over spines in the element
          
         } // End of Jacobian calculation
        
       } // End of if not boundary condition statement

      // ------------------------------------------------
      // FIFTH (AZIMUTHAL) MOMENTUM EQUATION: COSINE PART
      // ------------------------------------------------

      //if(offensive_eqn>=0) cout<<"residuals[offensive_eqn] = "<<residuals[offensive_eqn]<<endl;

      // Get local equation number of fifth velocity value at this node
      local_eqn = nodal_local_eqn(l,u_nodal_index[4]);

      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] -= scaled_re_st*r*dVCdt*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re_st*interpolated_RC*base_flow_duthetadt*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*base_flow_ur*interpolated_dVdRC*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[0]*interpolated_dVdRC*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_UC*base_flow_duthetadr*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dRCdt*base_flow_duthetadr*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RC*base_flow_ur*base_flow_duthetadr
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RC*mesh_velocity[0]*base_flow_duthetadr
         *testf_*Jbar*w;
        residuals[local_eqn] -=
         k*scaled_re*base_flow_utheta*interpolated_VS*testf_*Jbar*w;
        residuals[local_eqn] +=
         k*scaled_re*base_flow_utheta*base_flow_duthetadr*interpolated_RS
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         k*scaled_re*base_flow_utheta*base_flow_duthetadz*interpolated_ZS
         *testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*base_flow_utheta*interpolated_UC*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_VC*base_flow_ur*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*base_flow_uz*interpolated_dVdZC*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[1]*interpolated_dVdZC*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_WC*base_flow_duthetadz*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dZCdt*base_flow_duthetadz*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RC*base_flow_uz*base_flow_duthetadz
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RC*mesh_velocity[1]*base_flow_duthetadz
         *testf_*Jbar*w;
        residuals[local_eqn] += interpolated_RC*body_force[2]*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_inv_fr*interpolated_RC*G[2]*testf_*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*base_flow_duthetadr*dtestfdRC*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*interpolated_dVdRC*dtestfdr*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*r*base_flow_duthetadr*JhatC*dtestfdr*w;
        residuals[local_eqn] -=
         visc_ratio*interpolated_RC*base_flow_duthetadr*dtestfdr*Jbar*w;
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[0]*interpolated_US*dtestfdr*Jbar*w;
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[0]*base_flow_durdr*interpolated_RS*dtestfdr*Jbar*w;
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[0]*base_flow_durdz*interpolated_ZS*dtestfdr*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*Gamma[0]*base_flow_utheta*dtestfdRC*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*Gamma[0]*interpolated_VC*dtestfdr*Jbar*w;
        residuals[local_eqn] -= k*base_flow_p*dtestfdr*interpolated_RS*Jbar*w;
        residuals[local_eqn] -= k*base_flow_p*dtestfdz*interpolated_ZS*Jbar*w;
        residuals[local_eqn] -= k*interpolated_PS*testf_*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*k*k*interpolated_VC*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*k*k*base_flow_duthetadr*interpolated_RC
         *testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*k*k*base_flow_duthetadz*interpolated_ZC
         *testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*k*base_flow_ur*dtestfdr*interpolated_RS
         *Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*k*base_flow_ur*dtestfdz*interpolated_ZS
         *Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*k*interpolated_US*testf_*Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*k*interpolated_RS*base_flow_ur*testf_
         *Jbar*w/(r*r);
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[0]*interpolated_WS*dtestfdz*Jbar*w;
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[0]*base_flow_duzdr*interpolated_RS*dtestfdz*Jbar*w;
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[0]*base_flow_duzdz*interpolated_ZS*dtestfdz*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*base_flow_duthetadz*dtestfdZC*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*interpolated_dVdZC*dtestfdz*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*r*base_flow_duthetadz*JhatC*dtestfdz*w;
        residuals[local_eqn] -=
         visc_ratio*interpolated_RC*base_flow_duthetadz*dtestfdz*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*Gamma[0]*interpolated_dVdRC*testf_*Jbar*w;
        residuals[local_eqn] += visc_ratio*k*interpolated_US*testf_*Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*k*base_flow_durdr*interpolated_RS*testf_*Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*k*base_flow_durdz*interpolated_ZS*testf_*Jbar*w/r;
        residuals[local_eqn] -= visc_ratio*interpolated_VC*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*base_flow_utheta*interpolated_RC*testf_*Jbar*w/(r*r);
        residuals[local_eqn] -=
         scaled_re_st*r*base_flow_duthetadt*testf_*JhatC*w; 
        residuals[local_eqn] -=
         scaled_re*base_flow_utheta*base_flow_ur*testf_*JhatC*w; 
        residuals[local_eqn] += r*body_force[2]*testf_*JhatC*w;
        residuals[local_eqn] += scaled_re_inv_fr*r*G[2]*testf_*JhatC*w;
        residuals[local_eqn] -= k*base_flow_p*testf_*JhatS*w;
        residuals[local_eqn] +=
         k*visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_*JhatS*w/r;
        residuals[local_eqn] -= visc_ratio*base_flow_utheta*testf_*JhatC*w/r;
        
        
        // Calculate the Jacobian
        // ----------------------

        if(flag)
         {
          // Loop over the velocity shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Radial velocity component (cosine part) U_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[0]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_duthetadr*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*base_flow_utheta*psif[l2]*testf_*Jbar*w;
             }

            // Radial velocity component (sine part) U_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[1]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[0]*psif[l2]*dtestfdr*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*(1.0+Gamma[0])*k*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*k*psif[l2]*testf_*Jbar*w/r;
             }

            // Axial velocity component (cosine part) W_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[2]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_duthetadz*testf_*Jbar*w;
             }

            // Axial velocity component (sine part) W_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[3]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -= 
               k*visc_ratio*Gamma[0]*psif[l2]*dtestfdz*Jbar*w;
             }
            
            // Azimuthal velocity component (cosine part) V_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[4]);
            if(local_unknown >= 0)
             {
              if(flag==2)
               {
                // Add the mass matrix
                mass_matrix(local_eqn,local_unknown) +=
                 scaled_re_st*r*psif[l2]*testf_*Jbar*w;
               }
              
              // Add contributions to the Jacobian matrix
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->time_stepper_pt()->weight(1,0)*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_ur*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[0]*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*psif[l2]*base_flow_ur*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_uz*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[1]*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*group_B[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*Gamma[0]*psif[l2]*dtestfdr*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*k*k*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*group_A[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*Gamma[0]*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*psif[l2]*testf_*Jbar*w/r;
             }
            
            // Azimuthal velocity component (sine part) V_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[5]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*scaled_re*base_flow_utheta*psif[l2]*testf_*Jbar*w;
             }

            // Perturbation to radial nodal coord (cosine part) R_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[0]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*psif[l2]*base_flow_duthetadt*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*r*psif[l2]
              *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
              *base_flow_duthetadr*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re*psif[l2]*base_flow_ur*base_flow_duthetadr*testf_
              *Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*psif[l2]*mesh_velocity[0]*base_flow_duthetadr
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re*r*base_flow_uz*group_E[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*r*mesh_velocity[1]*group_E[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re*psif[l2]*base_flow_uz*base_flow_duthetadz
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*psif[l2]*mesh_velocity[1]*base_flow_duthetadz
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              psif[l2]*body_force[2]*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_inv_fr*psif[l2]*G[2]*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*base_flow_duthetadr*group_B[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*psif[l2]*base_flow_duthetadr*dtestfdr*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[0])*k*k*base_flow_duthetadr
              *psif[l2]*testf_*Jbar*w/r;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*base_flow_duthetadz*group_F(l,l2)*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*group_E[l2]*dtestfdz*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*base_flow_duthetadz*group_B[l2]*dtestfdz*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*psif[l2]*base_flow_duthetadz*dtestfdz*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*base_flow_utheta*psif[l2]*testf_*Jbar*w/(r*r);
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*r*base_flow_duthetadt*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re*base_flow_utheta*base_flow_ur*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              r*body_force[2]*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_inv_fr*r*G[2]*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*base_flow_utheta*testf_*group_B[l2]*w/r;
            }

            // Perturbation to radial nodal coord (sine part) R_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[1]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              k*scaled_re*base_flow_utheta*base_flow_duthetadr*psif[l2]
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              k*visc_ratio*Gamma[0]*base_flow_durdr*psif[l2]*dtestfdr*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              k*base_flow_p*dtestfdr*psif[l2]*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[0])*k*base_flow_ur*dtestfdr*psif[l2]
              *Jbar*w/r;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*(1.0+Gamma[0])*k*psif[l2]*base_flow_ur*testf_
              *Jbar*w/(r*r);
             jacobian(local_eqn,local_unknown) +=
              k*visc_ratio*Gamma[0]*base_flow_duzdr*psif[l2]*dtestfdz*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*k*base_flow_durdr*psif[l2]*testf_*Jbar*w/r;
             jacobian(local_eqn,local_unknown) -=
              k*base_flow_p*testf_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              k*visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_*group_B[l2]*w/r;
            }

            // Perturbation to axial nodal coord (cosine part) Z_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[2]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=
              scaled_re*r*base_flow_ur*group_E[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*r*mesh_velocity[0]*group_E[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_st*r*psif[l2]
              *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
              *base_flow_duthetadz*testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*r*base_flow_duthetadr*group_F(l,l2)*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*r*group_E[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*base_flow_duthetadr*group_A[l2]*dtestfdr*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*Gamma[0]*base_flow_utheta*group_F(l,l2)*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[0])*k*k*base_flow_duthetadz*psif[l2]
              *testf_*Jbar*w/r;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*r*base_flow_duthetadz*group_A[l2]*dtestfdz*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*Gamma[0]*group_E[l2]*testf_*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re_st*r*base_flow_duthetadt*testf_*group_A[l2]*w;
             jacobian(local_eqn,local_unknown) -=
              scaled_re*base_flow_utheta*base_flow_ur*testf_*group_A[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              r*body_force[2]*testf_*group_A[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              scaled_re_inv_fr*r*G[2]*testf_*group_A[l2]*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*base_flow_utheta*testf_*group_A[l2]*w/r;
            }

            // Perturbation to axial nodal coord (sine part) Z_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[3]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              k*scaled_re*base_flow_utheta*base_flow_duthetadz*psif[l2]
              *testf_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              k*visc_ratio*Gamma[0]*base_flow_durdz*psif[l2]*dtestfdr*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              k*base_flow_p*dtestfdz*psif[l2]*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              visc_ratio*(1.0+Gamma[0])*k*base_flow_ur*dtestfdz*psif[l2]
              *Jbar*w/r;
             jacobian(local_eqn,local_unknown) +=
              k*visc_ratio*Gamma[0]*base_flow_duzdz*psif[l2]*dtestfdz*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              visc_ratio*k*base_flow_durdz*psif[l2]*testf_*Jbar*w/r;
             jacobian(local_eqn,local_unknown) -=
              k*base_flow_p*testf_*group_A[l2]*w;
             jacobian(local_eqn,local_unknown) +=
              k*visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_*group_A[l2]*w/r;
            }

           } // End of loop over velocity shape functions
          
          // Now loop over pressure shape functions
          // (This is the contribution from pressure gradient)
          for(unsigned l2=0;l2<n_pres;l2++)
           {
            // Cosine part P_k^C has no contribution
            
            // Sine part P_k^S
            local_unknown = p_local_eqn(l2,1);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -= k*psip[l2]*testf_*Jbar*w;
             }
           } // End of loop over pressure shape functions

          // Now loop over the spines in the element
          for(unsigned l2=0;l2<n_p;l2++)
           {
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown = 
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,0);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum -= scaled_re*r*base_flow_ur*group_E[l2+(3*a)]*testf_*w;
                sum +=
                 scaled_re_st*r*mesh_velocity[0]*group_E[l2+(3*a)]*testf_*w;
                sum +=
                 scaled_re_st*r*psif[l2+(3*a)]
                 *node_pt(l2+(3*a))->position_time_stepper_pt()->weight(1,0)
                 *base_flow_duthetadz*testf_*Jbar*w;
                sum -=
                 visc_ratio*r*base_flow_duthetadr*group_F(l,(l2+(3*a)))*w;
                sum -= visc_ratio*r*group_E[l2+(3*a)]*dtestfdr*w;
                sum +=
                 visc_ratio*r*base_flow_duthetadr*group_A[l2+(3*a)]*dtestfdr*w;
                sum +=
                 visc_ratio*Gamma[0]*base_flow_utheta*group_F(l,(l2+(3*a)))*w;
                sum +=
                 visc_ratio*(1.0+Gamma[0])*k*k*base_flow_duthetadz
                 *psif[l2+(3*a)]*testf_*Jbar*w/r;
                sum +=
                 visc_ratio*r*base_flow_duthetadz*group_A[l2+(3*a)]*dtestfdz*w;
                sum += visc_ratio*Gamma[0]*group_E[l2+(3*a)]*testf_*w;
                sum -=
                 scaled_re_st*r*base_flow_duthetadt*testf_*group_A[l2+(3*a)]*w;
                sum -=
                 scaled_re*base_flow_utheta*base_flow_ur*testf_
                 *group_A[l2+(3*a)]*w;
                sum += r*body_force[2]*testf_*group_A[l2+(3*a)]*w;
                sum += scaled_re_inv_fr*r*G[2]*testf_*group_A[l2+(3*a)]*w;
                sum -=
                 visc_ratio*base_flow_utheta*testf_*group_A[l2+(3*a)]*w/r;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
            
            // Perturbed spine "height" (sine part) H_k^S
            local_unknown = 
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,1);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum +=
                 k*scaled_re*base_flow_utheta*base_flow_duthetadz
                 *psif[l2+(3*a)]*testf_*Jbar*w;
                sum +=
                 k*visc_ratio*Gamma[0]*base_flow_durdz*psif[l2+(3*a)]
                 *dtestfdr*Jbar*w;
                sum -= k*base_flow_p*dtestfdz*psif[l2+(3*a)]*Jbar*w;
                sum +=
                 visc_ratio*(1.0+Gamma[0])*k*base_flow_ur*dtestfdz
                 *psif[l2+(3*a)]*Jbar*w/r;
                sum +=
                 k*visc_ratio*Gamma[0]*base_flow_duzdz*psif[l2+(3*a)]
                 *dtestfdz*Jbar*w;
                sum -=
                 visc_ratio*k*base_flow_durdz*psif[l2+(3*a)]*testf_*Jbar*w/r;
                sum -= k*base_flow_p*testf_*group_A[l2+(3*a)]*w;
                sum +=
                 k*visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_
                 *group_A[l2+(3*a)]*w/r;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
           } // End of loop over spines in the element
          
         } // End of Jacobian calculation
        
       } // End of if not boundary condition statement
      
      // ----------------------------------------------
      // SIXTH (AZIMUTHAL) MOMENTUM EQUATION: SINE PART
      // ----------------------------------------------

      //if(offensive_eqn>=0) cout<<"residuals[offensive_eqn] = "<<residuals[offensive_eqn]<<endl;

      // Get local equation number of sixth velocity value at this node
      local_eqn = nodal_local_eqn(l,u_nodal_index[5]);

      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] -= scaled_re_st*r*dVSdt*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re_st*interpolated_RS*base_flow_duthetadt*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*base_flow_ur*interpolated_dVdRS*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[0]*interpolated_dVdRS*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_US*base_flow_duthetadr*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dRSdt*base_flow_duthetadr*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RS*base_flow_ur*base_flow_duthetadr
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RS*mesh_velocity[0]*base_flow_duthetadr
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         k*scaled_re*base_flow_utheta*interpolated_VC*testf_*Jbar*w;
        residuals[local_eqn] -=
         k*scaled_re*base_flow_utheta*base_flow_duthetadr*interpolated_RC
         *testf_*Jbar*w;
        residuals[local_eqn] -=
         k*scaled_re*base_flow_utheta*base_flow_duthetadz*interpolated_ZC
         *testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*base_flow_utheta*interpolated_US*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_VS*base_flow_ur*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*base_flow_uz*interpolated_dVdZS*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*mesh_velocity[1]*interpolated_dVdZS*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*r*interpolated_WS*base_flow_duthetadz*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*r*dZSdt*base_flow_duthetadz*testf_*Jbar*w;
        residuals[local_eqn] -=
         scaled_re*interpolated_RS*base_flow_uz*base_flow_duthetadz
         *testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_st*interpolated_RS*mesh_velocity[1]*base_flow_duthetadz
         *testf_*Jbar*w;
        residuals[local_eqn] += interpolated_RS*body_force[2]*testf_*Jbar*w;
        residuals[local_eqn] +=
         scaled_re_inv_fr*interpolated_RS*G[2]*testf_*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*base_flow_duthetadr*dtestfdRS*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*interpolated_dVdRS*dtestfdr*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*r*base_flow_duthetadr*JhatS*dtestfdr*w;
        residuals[local_eqn] -=
         visc_ratio*interpolated_RS*base_flow_duthetadr*dtestfdr*Jbar*w;
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[0]*interpolated_UC*dtestfdr*Jbar*w;
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[0]*base_flow_durdr*interpolated_RC*dtestfdr*Jbar*w;
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[0]*base_flow_durdz*interpolated_ZC*dtestfdr*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*Gamma[0]*base_flow_utheta*dtestfdRS*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*Gamma[0]*interpolated_VS*dtestfdr*Jbar*w;
        residuals[local_eqn] += k*base_flow_p*dtestfdr*interpolated_RC*Jbar*w;
        residuals[local_eqn] += k*base_flow_p*dtestfdz*interpolated_ZC*Jbar*w;
        residuals[local_eqn] += k*interpolated_PC*testf_*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*k*k*interpolated_VS*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*k*k*base_flow_duthetadr*interpolated_RS
         *testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*k*k*base_flow_duthetadz*interpolated_ZS
         *testf_*Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*k*base_flow_ur*dtestfdr*interpolated_RC
         *Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*k*base_flow_ur*dtestfdz*interpolated_ZC
         *Jbar*w/r;
        residuals[local_eqn] -=
         visc_ratio*(1.0+Gamma[0])*k*interpolated_UC*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*(1.0+Gamma[0])*k*interpolated_RC*base_flow_ur
         *testf_*Jbar*w/(r*r);
        residuals[local_eqn] +=
         k*visc_ratio*Gamma[0]*interpolated_WC*dtestfdz*Jbar*w;
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[0]*base_flow_duzdr*interpolated_RC*dtestfdz
         *Jbar*w;
        residuals[local_eqn] -=
         k*visc_ratio*Gamma[0]*base_flow_duzdz*interpolated_ZC*dtestfdz
         *Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*base_flow_duthetadz*dtestfdZS*Jbar*w;
        residuals[local_eqn] -=
         visc_ratio*r*interpolated_dVdZS*dtestfdz*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*r*base_flow_duthetadz*JhatS*dtestfdz*w;
        residuals[local_eqn] -=
         visc_ratio*interpolated_RS*base_flow_duthetadz*dtestfdz*Jbar*w;
        residuals[local_eqn] +=
         visc_ratio*Gamma[0]*interpolated_dVdRS*testf_*Jbar*w;
        residuals[local_eqn] -= visc_ratio*k*interpolated_UC*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*base_flow_durdr*interpolated_RC*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*k*base_flow_durdz*interpolated_ZC*testf_*Jbar*w/r;
        residuals[local_eqn] -= visc_ratio*interpolated_VS*testf_*Jbar*w/r;
        residuals[local_eqn] +=
         visc_ratio*base_flow_utheta*interpolated_RS*testf_*Jbar*w/(r*r);
        residuals[local_eqn] -=
         scaled_re_st*r*base_flow_duthetadt*testf_*JhatS*w; 
        residuals[local_eqn] -=
         scaled_re*base_flow_utheta*base_flow_ur*testf_*JhatS*w; 
        residuals[local_eqn] += r*body_force[2]*testf_*JhatS*w;
        residuals[local_eqn] += scaled_re_inv_fr*r*G[2]*testf_*JhatS*w;
        residuals[local_eqn] += k*base_flow_p*testf_*JhatC*w;
        residuals[local_eqn] -=
         k*visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_*JhatC*w/r;
        residuals[local_eqn] -= visc_ratio*base_flow_utheta*testf_*JhatS*w/r;
        

        // Calculate the Jacobian
        // ----------------------

        if(flag)
         {
          // Loop over the velocity shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Radial velocity component (cosine part) U_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[0]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               k*visc_ratio*Gamma[0]*psif[l2]*dtestfdr*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*k*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*k*psif[l2]*testf_*Jbar*w/r;
             }

            // Radial velocity component (sine part) U_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[1]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_duthetadr*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*base_flow_utheta*psif[l2]*testf_*Jbar*w;
             }

            // Axial velocity component (cosine part) W_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[2]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               k*visc_ratio*Gamma[0]*psif[l2]*dtestfdz*Jbar*w;
             }

            // Axial velocity component (sine part) W_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[3]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*psif[l2]*base_flow_duthetadz*testf_*Jbar*w;
             }
            
            // Azimuthal velocity component (cosine part) V_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[4]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               k*scaled_re*base_flow_utheta*psif[l2]*testf_*Jbar*w;
             }

            // Azimuthal velocity component (sine part) V_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[5]);
            if(local_unknown >= 0)
             {
              if(flag==2)
               {
                // Add the mass matrix
                mass_matrix(local_eqn,local_unknown) +=
                 scaled_re_st*r*psif[l2]*testf_*Jbar*w;
               }
              
              // Add contributions to the Jacobian matrix
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->time_stepper_pt()->weight(1,0)*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_ur*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[0]*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*psif[l2]*base_flow_ur*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_uz*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[1]*group_A[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*group_B[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*Gamma[0]*psif[l2]*dtestfdr*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*k*k*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*group_A[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*Gamma[0]*group_B[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*psif[l2]*testf_*Jbar*w/r;
             }
            
            // Perturbation to radial nodal coord (cosine part) R_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[0]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*scaled_re*base_flow_utheta*base_flow_duthetadr
               *psif[l2]*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[0]*base_flow_durdr*psif[l2]*dtestfdr*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               k*base_flow_p*dtestfdr*psif[l2]*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*k*base_flow_ur*dtestfdr
               *psif[l2]*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*(1.0+Gamma[0])*k*psif[l2]*base_flow_ur
               *testf_*Jbar*w/(r*r);
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[0]*base_flow_duzdr*psif[l2]*dtestfdz*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*k*base_flow_durdr*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               k*base_flow_p*testf_*group_B[l2]*w;
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_
               *group_B[l2]*w/r;
             }
            
            // Perturbation to radial nodal coord (sine part) R_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[1]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*psif[l2]*base_flow_duthetadt*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
               *base_flow_duthetadr*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*psif[l2]*base_flow_ur*base_flow_duthetadr
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*psif[l2]*mesh_velocity[0]*base_flow_duthetadr
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re*r*base_flow_uz*group_E[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*mesh_velocity[1]*group_E[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*psif[l2]*base_flow_uz*base_flow_duthetadz
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*psif[l2]*mesh_velocity[1]*base_flow_duthetadz
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               psif[l2]*body_force[2]*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_inv_fr*psif[l2]*G[2]*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*base_flow_duthetadr*group_B[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*psif[l2]*base_flow_duthetadr*dtestfdr*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*(1.0+Gamma[0])*k*k*base_flow_duthetadr*psif[l2]
               *testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*base_flow_duthetadz*group_F(l,l2)*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*group_E[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*base_flow_duthetadz*group_B[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*psif[l2]*base_flow_duthetadz*dtestfdz*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*base_flow_utheta*psif[l2]*testf_*Jbar*w/(r*r);
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*base_flow_duthetadt*testf_*group_B[l2]*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*base_flow_utheta*base_flow_ur*testf_*group_B[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               r*body_force[2]*testf_*group_B[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_inv_fr*r*G[2]*testf_*group_B[l2]*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*base_flow_utheta*testf_*group_B[l2]*w/r;
             }
            
            // Perturbation to axial nodal coord (cosine part) Z_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[2]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*scaled_re*base_flow_utheta*base_flow_duthetadz*psif[l2]
               *testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[0]*base_flow_durdz*psif[l2]*dtestfdr*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               k*base_flow_p*dtestfdz*psif[l2]*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*(1.0+Gamma[0])*k*base_flow_ur*dtestfdz*psif[l2]
               *Jbar*w/r;
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*Gamma[0]*base_flow_duzdz*psif[l2]*dtestfdz*Jbar*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*k*base_flow_durdz*psif[l2]*testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               k*base_flow_p*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) -=
               k*visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_*group_A[l2]*w/r;
             }
            
            // Perturbation to axial nodal coord (sine part) Z_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[3]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               scaled_re*r*base_flow_ur*group_E[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*mesh_velocity[0]*group_E[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_st*r*psif[l2]
               *node_pt(l2)->position_time_stepper_pt()->weight(1,0)
               *base_flow_duthetadz*testf_*Jbar*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*base_flow_duthetadr*group_F(l,l2)*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*r*group_E[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*base_flow_duthetadr*group_A[l2]*dtestfdr*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*Gamma[0]*base_flow_utheta*group_F(l,l2)*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*(1.0+Gamma[0])*k*k*base_flow_duthetadz*psif[l2]
               *testf_*Jbar*w/r;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*r*base_flow_duthetadz*group_A[l2]*dtestfdz*w;
              jacobian(local_eqn,local_unknown) +=
               visc_ratio*Gamma[0]*group_E[l2]*testf_*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re_st*r*base_flow_duthetadt*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) -=
               scaled_re*base_flow_utheta*base_flow_ur*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               r*body_force[2]*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) +=
               scaled_re_inv_fr*r*G[2]*testf_*group_A[l2]*w;
              jacobian(local_eqn,local_unknown) -=
               visc_ratio*base_flow_utheta*testf_*group_A[l2]*w/r;
             }
            
           } // End of loop over velocity shape functions
          
          // Now loop over pressure shape functions
          // (This is the contribution from pressure gradient)
          for(unsigned l2=0;l2<n_pres;l2++)
           {
            // Cosine part P_k^C
            local_unknown = p_local_eqn(l2,0);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) += k*psip[l2]*testf_*Jbar*w;
             }

            // Sine part P_k^S has no contribution

           } // End of loop over pressure shape functions

          // Now loop over the spines in the element
          for(unsigned l2=0;l2<n_p;l2++)
           {
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown = 
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,0);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum -=
                 k*scaled_re*base_flow_utheta*base_flow_duthetadz
                 *psif[l2+(3*a)]*testf_*Jbar*w;
                sum -=
                 k*visc_ratio*Gamma[0]*base_flow_durdz*psif[l2+(3*a)]
                 *dtestfdr*Jbar*w;
                sum += k*base_flow_p*dtestfdz*psif[l2+(3*a)]*Jbar*w;
                sum -=
                 visc_ratio*(1.0+Gamma[0])*k*base_flow_ur*dtestfdz
                 *psif[l2+(3*a)]*Jbar*w/r;
                sum -=
                 k*visc_ratio*Gamma[0]*base_flow_duzdz*psif[l2+(3*a)]
                 *dtestfdz*Jbar*w;
                sum +=
                 visc_ratio*k*base_flow_durdz*psif[l2+(3*a)]*testf_*Jbar*w/r;
                sum += k*base_flow_p*testf_*group_A[l2+(3*a)]*w;
                sum -=
                 k*visc_ratio*(1.0+Gamma[0])*base_flow_ur*testf_
                 *group_A[l2+(3*a)]*w/r;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
            
            // Perturbed spine "height" (sine part) H_k^S
            local_unknown = 
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,1);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum -= scaled_re*r*base_flow_ur*group_E[l2+(3*a)]*testf_*w;
                sum +=
                 scaled_re_st*r*mesh_velocity[0]*group_E[l2+(3*a)]*testf_*w;
                sum +=
                 scaled_re_st*r*psif[l2+(3*a)]
                 *node_pt(l2+(3*a))->position_time_stepper_pt()->weight(1,0)
                 *base_flow_duthetadz*testf_*Jbar*w;
                sum -=
                 visc_ratio*r*base_flow_duthetadr*group_F(l,(l2+(3*a)))*w;
                sum -= visc_ratio*r*group_E[l2+(3*a)]*dtestfdr*w;
                sum +=
                 visc_ratio*r*base_flow_duthetadr*group_A[l2+(3*a)]
                 *dtestfdr*w;
                sum +=
                 visc_ratio*Gamma[0]*base_flow_utheta*group_F(l,(l2+(3*a)))*w;
                sum +=
                 visc_ratio*(1.0+Gamma[0])*k*k*base_flow_duthetadz
                 *psif[l2+(3*a)]*testf_*Jbar*w/r;
                sum +=
                 visc_ratio*r*base_flow_duthetadz*group_A[l2+(3*a)]*dtestfdz*w;
                sum += visc_ratio*Gamma[0]*group_E[l2+(3*a)]*testf_*w;
                sum -=
                 scaled_re_st*r*base_flow_duthetadt*testf_*group_A[l2+(3*a)]*w;
                sum -=
                 scaled_re*base_flow_utheta*base_flow_ur*testf_
                 *group_A[l2+(3*a)]*w;
                sum += r*body_force[2]*testf_*group_A[l2+(3*a)]*w;
                sum += scaled_re_inv_fr*r*G[2]*testf_*group_A[l2+(3*a)]*w;
                sum -=
                 visc_ratio*base_flow_utheta*testf_*group_A[l2+(3*a)]*w/r;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
           } // End of loop over spines in the element
          
         } // End of Jacobian calculation
        
       } // End of if not boundary condition statement

     } // End of loop over fluid test functions
    

    // ====================
    // CONTINUITY EQUATIONS
    // ====================

    // Loop over the pressure test functions
    for(unsigned l=0;l<n_pres;l++)
     {
      // Cache the test function
      const double testp_ = testp[l];

      // --------------------------------------
      // FIRST CONTINUITY EQUATION: COSINE PART
      // --------------------------------------

      // Get local equation number of first pressure value at this node
      local_eqn = p_local_eqn(l,0);

      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] += r*interpolated_dUdRC*testp_*Jbar*w;
        residuals[local_eqn] += interpolated_RC*base_flow_durdr*testp_*Jbar*w;
        residuals[local_eqn] += interpolated_UC*testp_*Jbar*w;
        residuals[local_eqn] += k*interpolated_VS*testp_*Jbar*w;
        residuals[local_eqn] -=
         k*base_flow_duthetadr*interpolated_RS*testp_*Jbar*w;
        residuals[local_eqn] -=
         k*base_flow_duthetadz*interpolated_ZS*testp_*Jbar*w;
        residuals[local_eqn] += r*interpolated_dWdZC*testp_*Jbar*w;
        residuals[local_eqn] += interpolated_RC*base_flow_duzdz*testp_*Jbar*w;
        residuals[local_eqn] -= interpolated_RC*source*testp_*Jbar*w;
        residuals[local_eqn] += base_flow_ur*testp_*JhatC*w;
        residuals[local_eqn] -= source*r*testp_*JhatC*w;


        // Calculate the Jacobian
        // ----------------------

        if(flag)
         {
          // Loop over the velocity shape functions
          for(unsigned l2=0;l2<n_node;l2++)
           { 
            // Radial velocity component (cosine part) U_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[0]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) += r*group_B[l2]*testp_*w;
              jacobian(local_eqn,local_unknown) += psif[l2]*testp_*Jbar*w;
             }

            // Radial velocity component (sine part) U_k^S
            // has no contribution

            // Axial velocity component (cosine part) W_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[2]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) += r*group_A[l2]*testp_*w;
             }

            // Axial velocity component (sine part) W_k^S
            // has no contribution

            // Azimuthal velocity component (cosine part) V_k^C
            // has no contribution

            // Azimuthal velocity component (sine part) V_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[5]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) += k*psif[l2]*testp_*Jbar*w;
             }

            // Perturbation to radial nodal coord (cosine part) R_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[0]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              psif[l2]*base_flow_durdr*testp_*Jbar*w;
             jacobian(local_eqn,local_unknown) -= r*group_D[l2]*testp_*w;
             jacobian(local_eqn,local_unknown) +=
              psif[l2]*base_flow_duzdz*testp_*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              psif[l2]*source*testp_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              base_flow_ur*testp_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) -=
              source*r*testp_*group_B[l2]*w;
            }

            // Perturbation to radial nodal coord (sine part) R_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[1]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=
              k*base_flow_duthetadr*psif[l2]*testp_*Jbar*w;
            }

            // Perturbation to axial nodal coord (cosine part) Z_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[2]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += r*group_C[l2]*testp_*w;
             jacobian(local_eqn,local_unknown) +=
              base_flow_ur*testp_*group_A[l2]*w;
             jacobian(local_eqn,local_unknown) -=
              source*r*testp_*group_A[l2]*w;
            }

            // Perturbation to axial nodal coord (sine part) Z_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[3]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=
              k*base_flow_duthetadz*psif[l2]*testp_*Jbar*w;
            }

           } // End of loop over velocity shape functions

          // Real and imaginary pressure components, P_k^C and P_k^S,
          // have no contribution

          // Now loop over the spines in the element
          for(unsigned l2=0;l2<n_p;l2++)
           {
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown =
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,0);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum += r*group_C[l2+(3*a)]*testp_*w;
                sum += base_flow_ur*testp_*group_A[l2+(3*a)]*w;
                sum -= source*r*testp_*group_A[l2+(3*a)]*w;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }

            // Perturbed spine "height" (sine part) H_k^S
            local_unknown = 
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,1);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum -= k*base_flow_duthetadz*psif[l2+(3*a)]*testp_*Jbar*w;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
           } // End of loop over spines in the element

         } // End of Jacobian calculation
        
       } // End of if not boundary condition statement

      // -------------------------------------
      // SECOND CONTINUITY EQUATION: SINE PART
      // -------------------------------------

      // Get local equation number of second pressure value at this node
      local_eqn = p_local_eqn(l,1);

      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] += r*interpolated_dUdRS*testp_*Jbar*w;
        residuals[local_eqn] += interpolated_RS*base_flow_durdr*testp_*Jbar*w;
        residuals[local_eqn] += interpolated_US*testp_*Jbar*w;
        residuals[local_eqn] -= k*interpolated_VC*testp_*Jbar*w;
        residuals[local_eqn] +=
         k*base_flow_duthetadr*interpolated_RC*testp_*Jbar*w;
        residuals[local_eqn] +=
         k*base_flow_duthetadz*interpolated_ZC*testp_*Jbar*w;
        residuals[local_eqn] += r*interpolated_dWdZS*testp_*Jbar*w;
        residuals[local_eqn] += interpolated_RS*base_flow_duzdz*testp_*Jbar*w;
        residuals[local_eqn] -= interpolated_RS*source*testp_*Jbar*w;
        residuals[local_eqn] += base_flow_ur*testp_*JhatS*w;
        residuals[local_eqn] -= source*r*testp_*JhatS*w;


        // Calculate the Jacobian
        // ----------------------

        if(flag)
         {
          // Loop over the velocity shape functions
          for(unsigned l2=0;l2<n_node;l2++)
           { 
            // Radial velocity component (cosine part) U_k^C
            // has no contribution

            // Radial velocity component (sine part) U_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[1]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) += r*group_B[l2]*testp_*w;
              jacobian(local_eqn,local_unknown) += psif[l2]*testp_*Jbar*w;
             }

            // Axial velocity component (cosine part) W_k^C
            // has no contribution

            // Axial velocity component (sine part) W_k^S
            local_unknown = nodal_local_eqn(l2,u_nodal_index[3]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) += r*group_A[l2]*testp_*w;
             }

            // Azimuthal velocity component (cosine part) V_k^C
            local_unknown = nodal_local_eqn(l2,u_nodal_index[4]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -= k*psif[l2]*testp_*Jbar*w;
             }

            // Azimuthal velocity component (sine part) V_k^S
            // has no contribution

            // Perturbation to radial nodal coord (cosine part) R_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[0]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              k*base_flow_duthetadr*psif[l2]*testp_*Jbar*w;
            }

            // Perturbation to radial nodal coord (sine part) R_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[1]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              psif[l2]*base_flow_durdr*testp_*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              r*group_D[l2]*testp_*w;
             jacobian(local_eqn,local_unknown) +=
              psif[l2]*base_flow_duzdz*testp_*Jbar*w;
             jacobian(local_eqn,local_unknown) -=
              psif[l2]*source*testp_*Jbar*w;
             jacobian(local_eqn,local_unknown) +=
              base_flow_ur*testp_*group_B[l2]*w;
             jacobian(local_eqn,local_unknown) -=
              source*r*testp_*group_B[l2]*w;
            }

            // Perturbation to axial nodal coord (cosine part) Z_k^C
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[2]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              k*base_flow_duthetadz*psif[l2]*testp_*Jbar*w;
            }

            // Perturbation to axial nodal coord (sine part) Z_k^S
            local_unknown = nodal_local_eqn(l2,xhat_nodal_index[3]);
            if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += r*group_C[l2]*testp_*w;
             jacobian(local_eqn,local_unknown) +=
              base_flow_ur*testp_*group_A[l2]*w;
             jacobian(local_eqn,local_unknown) -=
              source*r*testp_*group_A[l2]*w;
            }

           } // End of loop over velocity shape functions

          // Real and imaginary pressure components, P_k^C and P_k^S,
          // have no contribution

          // Now loop over the spines in the element
          for(unsigned l2=0;l2<n_p;l2++)
           {
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown = 
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,0);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum += k*base_flow_duthetadz*psif[l2+(3*a)]*testp_*Jbar*w;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
            
            // Perturbed spine "height" (sine part) H_k^S
            local_unknown = 
             get_local_eqn_number_corresponding_to_geometric_dofs(l2,1);

            if(local_unknown >= 0)
             {
              for(unsigned a=0;a<3;a++)
               {
                double sum = 0.0;
                
                sum += r*group_C[l2+(3*a)]*testp_*w;
                sum += base_flow_ur*testp_*group_A[l2+(3*a)]*w;
                sum -= source*r*testp_*group_A[l2+(3*a)]*w;
                
                jacobian(local_eqn,local_unknown) +=
                 sum*effective_perturbed_spine_node_fraction[l2+(3*a)];
               }
             }
           } // End of loop over spines in the element

         } // End of Jacobian calculation
        
       } // End of if not boundary condition statement

     } // End of loop over pressure test functions

   } // End of loop over the integration points
  
 } // End of fill_in_generic_residual_contribution_lin_axi_nst
 


/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////


 /// Linearised axisymmetric Crouzeix-Raviart elements
 /// -------------------------------------------------


 //=======================================================================
 /// Set the data for the number of variables at each node.
 /// At each node we have six velocity components and then four
 /// components of perturbations to nodal positions.
 //=======================================================================
 const unsigned LinearisedAxisymmetricQCrouzeixRaviartElement::
 Initial_Nvalue[9]={10,10,10,10,10,10,10,10,10};
 
 
 
 //========================================================================
 /// Number of values (pinned or dofs) required at node n
 //========================================================================
 unsigned LinearisedAxisymmetricQCrouzeixRaviartElement::
 required_nvalue(const unsigned &n) const { return Initial_Nvalue[n]; }
 
 
 
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
 

 /// Linearised axisymmetric Taylor-Hood elements
 /// --------------------------------------------


 //=======================================================================
 /// Set the data for the number of variables at each node.
 /// At each node we have six velocity components and then four
 /// components of perturbations to nodal positions. We have an
 /// additional two pressure components at each "corner node".
 //=======================================================================
 const unsigned LinearisedAxisymmetricQTaylorHoodElement::
 Initial_Nvalue[9]={12,10,12,10,10,10,12,10,12};



 //=======================================================================
 /// Set the data for the pressure conversion array
 //=======================================================================
 const unsigned LinearisedAxisymmetricQTaylorHoodElement::
 Pconv[4]={0,2,6,8};



} // End of oomph namespace
