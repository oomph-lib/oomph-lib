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
// Non-inline functions for linearised axisymmetric Navier-Stokes
// interface elements

// oomph-lib includes
#include "linearised_axisymmetric_fluid_interface_elements.h"

namespace oomph
{

 //=======================================================================
 /// Physical constants default to zero
 //=======================================================================
 double LinearisedAxisymmetricFluidInterfaceElement::
 Default_Physical_Constant_Value = 1.0;

 //=======================================================================
 /// Azimuthal mode number defaults to zero
 //=======================================================================
 int LinearisedAxisymmetricFluidInterfaceElement::
 Default_Azimuthal_Mode_Number_Value = 0;

 //=======================================================================
 ///  Calculate the residuals for the linearised axisymmetric
 /// interface element
 //=======================================================================
 void LinearisedAxisymmetricFluidInterfaceElement::
 fill_in_generic_residual_contribution_interface(
  Vector<double> &residuals, 
  DenseMatrix<double> &jacobian, 
  unsigned flag)
 {
  std::ostringstream error_stream;
  
  error_stream << "This function has been deliberately broken, since it\n"
               << "has been overloaded with one in the (derived) spine\n"
               << "version of this element, which makes certain assumptions:\n"
               << "namely that we are using vertical spines to perform our\n"
               << "node-updating procedures. At some point a more general\n"
               << "framework will want to be written, and this function\n"
               << "(here!) replaced with one that implements the generic\n"
               << "mathematics, which is independent of the node-update\n"
               << "scheme which is employed (and supplied by derived classes)."
               << std::endl;
  
  throw OomphLibError(error_stream.str(),
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }



 //=======================================================================
 /// Overload the output function
 //=======================================================================
 void LinearisedAxisymmetricFluidInterfaceElement::
 output(std::ostream &outfile, const unsigned &n_plot) {}



 //=======================================================================
 /// Overload the output function
 //=======================================================================
 void LinearisedAxisymmetricFluidInterfaceElement::
 output(FILE* file_pt, const unsigned &n_plot)
 {
  // Set output Vector
  Vector<double> s(1);
 
  // Tecplot header info 
  fprintf(file_pt,"ZONE I=%i \n",n_plot);
  
  // Loop over plot points
  for(unsigned l=0;l<n_plot;l++)
   {
    s[0] = -1.0 + l*2.0/(n_plot-1);
    
    // Output the x,y, and the 6 velocity components
    for(unsigned i=0;i<2;i++) fprintf(file_pt,"%g ",this->interpolated_x(s,i));
    for(unsigned i=0;i<6;i++) fprintf(file_pt,"%g ",this->interpolated_u(s,i));

    // Output two dummy pressures
    fprintf(file_pt,"0.0 0.0 \n");
   }
  fprintf(file_pt,"\n");
 }
 


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


 
 //=======================================================================
 /// Output the perturbation to the base interface position (computed by
 /// these interface elements) in tecplot format, using nplot points in
 /// each coordinate direction.
 /// WARNING: THIS FUNCTION ASSUMES THAT THE PERTURBATION TO THE
 /// INTERFACE IS IN A VERTICAL DIRECTION ONLY
 //=======================================================================
 template<class ELEMENT>
 void PerturbedSpineLinearisedAxisymmetricFluidInterfaceElement<ELEMENT>::
 output_perturbation_to_interface(std::ostream &outfile,const unsigned &nplot)
 {
  // Provide storage for vector of local coordinates
  Vector<double> s(1);
  
  // Tecplot header info
  outfile << tecplot_zone_string(nplot);
  
  // Determine number of plot points
  const unsigned n_plot_points = nplot_points(nplot);
  
  // Loop over plot points
  for(unsigned iplot=0;iplot<n_plot_points;iplot++)
   {
    // Get local coordinates of plot point
    get_s_plot(iplot,nplot,s);
    
    // Output the global r coordinate
    outfile << interpolated_x(s,0) << " ";
    
    // Output perturbation to interface
    outfile << interpolated_H(s,0) << " " << interpolated_H(s,1) << " ";
    
    outfile << std::endl;   
   }
  outfile << std::endl;
  
  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,nplot);
 }
 


 //=======================================================================
 /// Output the *combined* interface position
 /// (base position + perturbation) in tecplot format, using nplot points
 /// in each coordinate direction.
 /// WARNING: THIS FUNCTION ASSUMES THAT THE PERTURBATION TO THE
 /// INTERFACE IS IN A VERTICAL DIRECTION ONLY
 //=======================================================================
 template<class ELEMENT>
 void PerturbedSpineLinearisedAxisymmetricFluidInterfaceElement<ELEMENT>::
 output_interface_position(std::ostream &outfile, const unsigned &nplot)
 {
  // Provide storage for vector of local coordinates
  Vector<double> s(1);
  
  // Tecplot header info
  outfile << tecplot_zone_string(nplot);
  
  // Determine number of plot points
  const unsigned n_plot_points = nplot_points(nplot);
  
  // Loop over plot points
  for(unsigned iplot=0;iplot<n_plot_points;iplot++)
   {
    // Get local coordinates of plot point
    get_s_plot(iplot,nplot,s);
    
    // Output the global r coordinate
    outfile << interpolated_x(s,0) << " ";
    
    // Output the actual interface position (unperturbed+perturbed)
    // We must ensure that we increase the precision for this!
    outfile.precision(14);
    outfile << (interpolated_x(s,1) + interpolated_H(s,0)) << " "
            << (interpolated_x(s,1) + interpolated_H(s,1)) << " ";
    outfile.precision(6);
    
    outfile << std::endl;   
   }
  outfile << std::endl;
  
  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,nplot);
 }



 //=======================================================================
 /// Calculate the residuals for the linearised axisymmetric interface
 /// element.
 /// WARNING: THIS FUNCTION ASSUMES THAT THE PERTURBATION TO THE
 /// INTERFACE IS IN A VERTICAL DIRECTION ONLY
 //=======================================================================
 template<class ELEMENT>
 void PerturbedSpineLinearisedAxisymmetricFluidInterfaceElement<ELEMENT>::
 fill_in_generic_residual_contribution_interface(
  Vector<double> &residuals, 
  DenseMatrix<double> &jacobian, 
  unsigned flag)
 {
  // Determine number of nodes in the element
  const unsigned n_node = this->nnode();
  
  // Set up memory for the shape functions and their derivative w.r.t. the
  // local coordinate s
  Shape psif(n_node);
  DShape dpsifds(n_node,1);
  
  // Set up memory for the test functions and their derivative w.r.t. the
  // local coordinate s
  Shape testf(n_node);
  DShape dtestfds(n_node,1);
  
  // Storage for the local coordinate
  Vector<double> s(1);
  
  // Storage for the local coordinate in the parent element which
  // corresponds to the local coordinate in this element
  Vector<double> s_parent(2);
  
  // Get a pointer to the parent element
  LinearisedAxisymmetricNavierStokesEquations* bulk_el_pt =
   dynamic_cast<LinearisedAxisymmetricNavierStokesEquations*>
   (bulk_element_pt());
  
  // Determine the number of nodes in the parent element
  const unsigned n_node_parent = bulk_el_pt->nnode();
  
  // Set up memory for the parent test functions and their derivatives
  Shape testf_parent(n_node_parent);
  DShape dtestfdx_parent(n_node_parent,2);
  
  // Determine number of integration points
  const unsigned n_intpt = this->integral_pt()->nweight();
  
  // Get physical variables from the element
  const double Ca = ca();
  const double St = st();
  const int k = azimuthal_mode_number();
  
  // Integers used to store the local equation and unknown numbers
  int local_eqn = 0, local_unknown = 0;
  
  // Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    // Get the local coordinate at the integration point
    s[0] = integral_pt()->knot(ipt,0);
    
    // Get the corresponding local coordinate in the parent (bulk) element
    this->get_local_coordinate_in_bulk(s,s_parent);
    
    // Get the integral weight
    const double w = this->integral_pt()->weight(ipt);
    
    // Call the derivatives of the shape function
    this->dshape_local_at_knot(ipt,psif,dpsifds);
    
    // Set test functions equal to shape functions
    testf = psif;
    dtestfds = dpsifds;
    
    // Find the test functions and derivatives of the parent
    (void)bulk_el_pt->dshape_eulerian(s_parent,testf_parent,dtestfdx_parent);
    
    // Define and zero the tangent Vectors
    double interpolated_t1[2] = {0.0,0.0};
    Vector<double> interpolated_x(2,0.0);
    double interpolated_dx_dt[2] = {0.0,0.0};
    
    // Provide storage for velocity unknowns
    // (only need four of them in these equations)
    double interpolated_UC = 0.0;
    double interpolated_US = 0.0;
    double interpolated_WC = 0.0;
    double interpolated_WS = 0.0;
    
    // Provide storage for perturbed spine "heights"
    double interpolated_HC = 0.0;
    double interpolated_HS = 0.0;
    
    // Provide storage for derivatives of perturbed spine "heights" w.r.t.
    // the local coordinate
    double interpolated_dHCds = 0.0;
    double interpolated_dHSds = 0.0;
    
    // Provide storage for derivatives of perturbed spine "heights" w.r.t.
    // time
    double interpolated_dHCdt = 0.0;
    double interpolated_dHSdt = 0.0;
    
    // Loop over the shape functions
    for(unsigned l=0;l<n_node;l++)
     {
      // Upcast from general node to PerturbedSpineNode
      PerturbedSpineNode* perturbed_spine_node_pt =
       dynamic_cast<PerturbedSpineNode*>(this->node_pt(l));
      
      // Cache the shape function and its derivative
      const double psif_ = psif(l);
      const double dpsifds_ = dpsifds(l,0);
      
      // Calculate interpolated velocity components
      interpolated_UC += u(l,0)*psif_;
      interpolated_US += u(l,1)*psif_;
      interpolated_WC += u(l,2)*psif_;
      interpolated_WS += u(l,3)*psif_;
      
      // Calculate interpolated perturbed spine "heights"
      interpolated_HC +=
       perturbed_spine_node_pt->perturbed_spine_pt()->height(0)*psif_;
      interpolated_HS +=
       perturbed_spine_node_pt->perturbed_spine_pt()->height(1)*psif_;
      
      // Calculate derivatives of interpolated perturbed spine "heights"
      interpolated_dHCds +=
       perturbed_spine_node_pt->perturbed_spine_pt()->height(0)*dpsifds_;
      interpolated_dHSds +=
       perturbed_spine_node_pt->perturbed_spine_pt()->height(1)*dpsifds_;
      interpolated_dHCdt += this->dH_dt(l,0)*psif_;
      interpolated_dHSdt += this->dH_dt(l,1)*psif_;
      
      // Loop over directional components
      for(unsigned i=0;i<2;i++)
       {
        interpolated_x[i] += this->nodal_position(l,i)*psif_;
        interpolated_dx_dt[i] += this->dnodal_position_dt(l,i)*psif_;
        
        // Calculate the tangent vector
        interpolated_t1[i] += this->nodal_position(l,i)*dpsifds_;
       }
     }
    
    
    // Get velocities from base flow problem
    // -------------------------------------
    
    // Allocate storage for the velocity components of the base state
    // solution (initialise to zero)
    Vector<double> base_flow_u(3,0.0);
    
    // Get the user-defined base state solution velocity components
    bulk_el_pt->get_base_flow_u(
     bulk_el_pt->node_pt(0)->time_stepper_pt()->time(),
     ipt,
     interpolated_x,
     base_flow_u);
    
    // Cache base flow velocities
    const double interpolated_ur = base_flow_u[0];
    const double interpolated_utheta = base_flow_u[2];
    
    // The first positional coordinate is the radial coordinate
    const double r = interpolated_x[0];
    
    // Calculate the length of the tangent Vector
    const double tlength = interpolated_t1[0]*interpolated_t1[0] + 
     interpolated_t1[1]*interpolated_t1[1];
    
    // Set the Jacobian of the line element
    const double J = sqrt(tlength);
    
    // Normalise the tangent Vector
    interpolated_t1[0] /= J; interpolated_t1[1] /= J;
    
    // Now calculate the normal Vector
    Vector<double> interpolated_n(2);
    outer_unit_normal(ipt,interpolated_n);
    
    // Also get the (possibly variable) surface tension
    const double Sigma = this->sigma(s);
    
    // Loop over the test functions
    for(unsigned l=0;l<n_node;l++)
     {
      // Cache test function and its derivative
      const double testf_ = testf(l);
      const double dtestfds_ = dtestfds(l,0);
      
      // =================================================
      // START OF DYNAMIC BOUNDARY CONDITION CONTRIBUTIONS
      // =================================================
      
      // -------------------------------------------------------------
      // Contribution to first (radial) momentum equation: cosine part
      // -------------------------------------------------------------
      
      // Get local equation number of first velocity value at this node
      local_eqn = this->nodal_local_eqn(l,this->U_index_interface[0]);
      
      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] -=
         (Sigma/Ca)*interpolated_t1[1]*interpolated_dHCds*testf_*w;
        residuals[local_eqn] +=
         (Sigma/Ca)*r*interpolated_t1[1]*interpolated_t1[0]
         *interpolated_dHCds*dtestfds_*w/J;
        residuals[local_eqn] +=
         (Sigma/Ca)*k*k*J*interpolated_t1[1]*interpolated_t1[0]
         *interpolated_HC*testf_*w/r;
        
        // Calculate the Jacobian
        // ----------------------
        if(flag)
         {
          // Loop over the shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Cache the shape function and its derivative
            const double psif_ = psif(l2);
            const double dpsifds_ = dpsifds(l2,0);
            
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown = kinematic_local_eqn(l2,0);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               (Sigma/Ca)*interpolated_t1[1]*dpsifds_*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               (Sigma/Ca)*r*interpolated_t1[1]*interpolated_t1[0]
               *dpsifds_*dtestfds_*w/J;
              jacobian(local_eqn,local_unknown) +=
               (Sigma/Ca)*k*k*J*interpolated_t1[1]*interpolated_t1[0]
               *psif_*testf_*w/r;
             }
            
            // Perturbed spine "height" (sine part) H_K^S
            // has no contribution
            
           } // End of loop over shape functions
         } // End of Jacobian calculation         
       } // End of if not boundary condition statement
      
      // ------------------------------------------------------------
      // Contribution to second (radial) momentum equation: sine part
      // ------------------------------------------------------------
      
      // Get local equation number of second velocity value at this node
      local_eqn = this->nodal_local_eqn(l,this->U_index_interface[1]);
      
      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] -=
         (Sigma/Ca)*interpolated_t1[1]*interpolated_dHSds*testf_*w;
        residuals[local_eqn] +=
         (Sigma/Ca)*r*interpolated_t1[1]*interpolated_t1[0]
         *interpolated_dHSds*dtestfds_*w/J;
        residuals[local_eqn] += (Sigma/Ca)*k*k*J*interpolated_t1[1]
         *interpolated_t1[0]*interpolated_HS*testf_*w/r;
        
        // Calculate the Jacobian
        // ----------------------
        if(flag)
         {
          // Loop over the shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Cache the shape function and its derivative
            const double psif_ = psif(l2);
            const double dpsifds_ = dpsifds(l2,0);
            
            // Perturbed spine "height" (cosine part) H_k^C
            // has no contribution
            
            // Perturbed spine "height" (sine part) H_K^S
            local_unknown = kinematic_local_eqn(l2,1);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               (Sigma/Ca)*interpolated_t1[1]*dpsifds_*testf_*w;
              jacobian(local_eqn,local_unknown) +=
               (Sigma/Ca)*r*interpolated_t1[1]*interpolated_t1[0]
               *dpsifds_*dtestfds_*w/J;
              jacobian(local_eqn,local_unknown) +=
               (Sigma/Ca)*k*k*J*interpolated_t1[1]*interpolated_t1[0]
               *psif_*testf_*w/r;
             }
           } // End of loop over shape functions
         } // End of Jacobian calculation         
       } // End of if not boundary condition statement
      
      // ------------------------------------------------------------
      // Contribution to third (axial) momentum equation: cosine part
      // ------------------------------------------------------------
      
      // Get local equation number of third velocity value at this node
      local_eqn = this->nodal_local_eqn(l,this->U_index_interface[2]);
      
      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] -=
         (Sigma/Ca)*r*interpolated_t1[0]*interpolated_t1[0]
         *interpolated_dHCds*dtestfds_*w/J;
        residuals[local_eqn] -= (Sigma/Ca)*k*k*J*interpolated_t1[0]
         *interpolated_t1[0]*interpolated_HC*testf_*w/r;
        
        // Calculate the Jacobian
        // ----------------------
        if(flag)
         {
          // Loop over the shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Cache the shape function and its derivative
            const double psif_ = psif(l2);
            const double dpsifds_ = dpsifds(l2,0);
            
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown = kinematic_local_eqn(l2,0);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               (Sigma/Ca)*r*interpolated_t1[0]*interpolated_t1[0]
               *dpsifds_*dtestfds_*w/J;
              jacobian(local_eqn,local_unknown) -=
               (Sigma/Ca)*k*k*J*interpolated_t1[0]*interpolated_t1[0]
               *psif_*testf_*w/r;
             }
            
            // Perturbed spine "height" (sine part) H_K^S
            // has no contribution
            
           } // End of loop over shape functions
         } // End of Jacobian calculation         
       } // End of if not boundary condition statement
      
      // -----------------------------------------------------------
      // Contribution to fourth (axial) momentum equation: sine part
      // -----------------------------------------------------------
      
      // Get local equation number of fourth velocity value at this node
      local_eqn = this->nodal_local_eqn(l,this->U_index_interface[3]);
      
      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] -=
         (Sigma/Ca)*r*interpolated_t1[0]*interpolated_t1[0]
         *interpolated_dHSds*dtestfds_*w/J;
        residuals[local_eqn] -=
         (Sigma/Ca)*k*k*J*interpolated_t1[0]*interpolated_t1[0]
         *interpolated_HS*testf_*w/r;
        
        // Calculate the Jacobian
        // ----------------------
        if(flag)
         {
          // Loop over the shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Cache the shape function and its derivative
            const double psif_ = psif(l2);
            const double dpsifds_ = dpsifds(l2,0);
            
            // Perturbed spine "height" (cosine part) H_k^C
            // has no contribution
            
            // Perturbed spine "height" (sine part) H_K^S
            local_unknown = kinematic_local_eqn(l2,1);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               (Sigma/Ca)*r*interpolated_t1[0]*interpolated_t1[0]*dpsifds_
               *dtestfds_*w/J;
              jacobian(local_eqn,local_unknown) -=
               (Sigma/Ca)*k*k*J*interpolated_t1[0]*interpolated_t1[0]*psif_
               *testf_*w/r;
             }
           } // End of loop over shape functions
         } // End of Jacobian calculation         
       } // End of if not boundary condition statement
      
      // ----------------------------------------------------------------
      // Contribution to fifth (azimuthal) momentum equation: cosine part
      // ----------------------------------------------------------------
      
      // Get local equation number of fifth velocity value at this node
      local_eqn = this->nodal_local_eqn(l,this->U_index_interface[4]);
      
      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] -=
         (Sigma/Ca)*2.0*k*J*interpolated_t1[1]*interpolated_t1[0]
         *interpolated_HS*testf_*w/(r*r);
        residuals[local_eqn] +=
         (Sigma/Ca)*k*interpolated_t1[1]*interpolated_HS*dtestfds_*w/r;
        residuals[local_eqn] +=
         (Sigma/Ca)*k*interpolated_t1[1]*interpolated_dHSds*testf_*w/r;
        
        // Calculate the Jacobian
        // ----------------------
        if(flag)
         {
          // Loop over the shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Cache the shape function and its derivative
            const double psif_ = psif(l2);
            const double dpsifds_ = dpsifds(l2,0);
            
            // Perturbed spine "height" (cosine part) H_k^C
            // has no contribution
            
            // Perturbed spine "height" (sine part) H_K^S
            local_unknown = kinematic_local_eqn(l2,1);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               (Sigma/Ca)*2.0*k*J*interpolated_t1[1]*interpolated_t1[0]
               *psif_*testf_*w/(r*r);
              jacobian(local_eqn,local_unknown) +=
               (Sigma/Ca)*k*interpolated_t1[1]*psif_*dtestfds_*w/r;
              jacobian(local_eqn,local_unknown) +=
               (Sigma/Ca)*k*interpolated_t1[1]*dpsifds_*testf_*w/r;
             }
           } // End of loop over shape functions
         } // End of Jacobian calculation         
       } // End of if not boundary condition statement
      
      // --------------------------------------------------------------
      // Contribution to sixth (azimuthal) momentum equation: sine part
      // --------------------------------------------------------------
      
      // Get local equation number of sixth velocity value at this node
      local_eqn = this->nodal_local_eqn(l,this->U_index_interface[5]);
      
      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        residuals[local_eqn] +=
         (Sigma/Ca)*2.0*k*J*interpolated_t1[1]*interpolated_t1[0]
         *interpolated_HC*testf_*w/(r*r);
        residuals[local_eqn] -=
         (Sigma/Ca)*k*interpolated_t1[1]*interpolated_HC*dtestfds_*w/r;
        residuals[local_eqn] -=
         (Sigma/Ca)*k*interpolated_t1[1]*interpolated_dHCds*testf_*w/r;
        
        // Calculate the Jacobian
        // ----------------------
        if(flag)
         {
          // Loop over the shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Cache the shape function and its derivative
            const double psif_ = psif(l2);
            const double dpsifds_ = dpsifds(l2,0);
            
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown = kinematic_local_eqn(l2,0);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               (Sigma/Ca)*2.0*k*J*interpolated_t1[1]*interpolated_t1[0]
               *psif_*testf_*w/(r*r);
              jacobian(local_eqn,local_unknown) -=
               (Sigma/Ca)*k*interpolated_t1[1]*psif_*dtestfds_*w/r;
              jacobian(local_eqn,local_unknown) -=
               (Sigma/Ca)*k*interpolated_t1[1]*dpsifds_*testf_*w/r;
             }
            
            // Perturbed spine "height" (sine part) H_K^S
            // has no contribution
            
           } // End of loop over shape functions
         } // End of Jacobian calculation         
       } // End of if not boundary condition statement
      
      // ===================================================
      // START OF KINEMATIC BOUNDARY CONDITION CONTRIBUTIONS
      // ===================================================
      
      // Using the same shape functions for the spines, so can stay inside
      // the loop over test functions
      
      // ------------------------------------------------
      // First kinematic boundary condition (cosine part)
      // ------------------------------------------------
      
      // Get local equation number of first perturbed spine "height"
      // (this is the cosine part H_k^C)
      local_eqn = kinematic_local_eqn(l,0);
      
      // If the spine is not a boundary condition
      if(local_eqn >= 0) 
       {
        residuals[local_eqn] -=
         r*interpolated_ur*interpolated_dHCds*testf_*w;
        
        residuals[local_eqn] -=
         k*interpolated_utheta*interpolated_t1[0]*interpolated_HS*testf_*w*J;
        
        residuals[local_eqn] -=
         r*interpolated_t1[1]*interpolated_UC*testf_*w*J;
        
        residuals[local_eqn] +=
         r*interpolated_t1[0]*interpolated_WC*testf_*w*J;
        
        residuals[local_eqn] -=
         St*r*interpolated_t1[0]*interpolated_dHCdt*testf_*w*J;
        
        // Calculate the Jacobian
        // ----------------------
        if(flag)
         {
          // Loop over the shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Cache the shape function and its derivative
            const double psif_ = psif(l2);
            const double dpsifds_ = dpsifds(l2,0);
            
            // Radial velocity component (cosine part) U_k^C
            local_unknown = this->nodal_local_eqn(l2,U_index_interface[0]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               r*interpolated_t1[1]*psif_*testf_*w*J;
             }
            
            // Axial velocity component (cosine part) W_k^C
            local_unknown = nodal_local_eqn(l2,U_index_interface[2]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               r*interpolated_t1[0]*psif_*testf_*w*J;
             }
            
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown = kinematic_local_eqn(l2,0);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               r*interpolated_ur*dpsifds_*testf_*w;
              
              jacobian(local_eqn,local_unknown) -= 
               St*r*interpolated_t1[0]
               *node_pt(l2)->time_stepper_pt()->weight(1,0)*psif_*testf_*w*J;
             }
            
            // Perturbed spine "height" (sine part) H_k^S
            local_unknown = kinematic_local_eqn(l2,1);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               k*interpolated_utheta*interpolated_t1[0]*psif_*testf_*w*J;
             }
            
           } // End of loop over shape functions
         } // End of Jacobian contribution
       } // End of if not boundary condition statement
      
      // -----------------------------------------------
      // Second kinematic boundary condition (sine part)
      // -----------------------------------------------
      
      // Get local equation number of second perturbed spine "height"
      // (this is the cosine part H_k^S)
      local_eqn = kinematic_local_eqn(l,1);
      
      // If the spine is not a boundary condition
      if(local_eqn >= 0) 
       {
        residuals[local_eqn] -=
         r*interpolated_ur*interpolated_dHSds*testf_*w;
        
        residuals[local_eqn] +=
         k*interpolated_utheta*interpolated_t1[0]*interpolated_HC*testf_*w*J;
        
        residuals[local_eqn] -=
         r*interpolated_t1[1]*interpolated_US*testf_*w*J;
        
        residuals[local_eqn] +=
         r*interpolated_t1[0]*interpolated_WS*testf_*w*J;
        
        residuals[local_eqn] -=
         St*r*interpolated_t1[0]*interpolated_dHSdt*testf_*w*J;
        
        // Add in the jacobian
        if(flag)
         {
          // Loop over velocity shape functions
          for(unsigned l2=0;l2<n_node;l2++)
           {
            // Cache the shape function and its derivative
            const double psif_ = psif(l2);
            const double dpsifds_ = dpsifds(l2,0);
            
            // Radial velocity component (sine part) U_k^S
            local_unknown = this->nodal_local_eqn(l2,U_index_interface[1]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               r*interpolated_t1[1]*psif_*testf_*w*J;
             }
            
            // Axial velocity component (sine part) W_k^S
            local_unknown = nodal_local_eqn(l2,U_index_interface[3]);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               r*interpolated_t1[0]*psif_*testf_*w*J;
             }
            
            // Perturbed spine "height" (cosine part) H_k^C
            local_unknown = kinematic_local_eqn(l2,0);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) +=
               k*interpolated_utheta*interpolated_t1[0]*psif_*testf_*w*J;
             }
            
            // Perturbed spine "height" (sine part) H_k^S
            local_unknown = kinematic_local_eqn(l2,1);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               r*interpolated_ur*dpsifds_*testf_*w;
              
              jacobian(local_eqn,local_unknown) -=
               St*r*interpolated_t1[0]
               *node_pt(l2)->time_stepper_pt()->weight(1,0)*psif_*testf_*w*J;
             }
            
           } // End of loop over shape functions
         } // End of Jacobian contribution
       } // End of if not boundary condition statement
      
     } // End of loop over test functions
    
   } // End of loop over integration points
  
 } // End of fill_in_generic_residual_contribution_interface
 
 
 
} // End of oomph namespace
