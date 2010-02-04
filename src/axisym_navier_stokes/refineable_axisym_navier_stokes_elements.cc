//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
#include "refineable_axisym_navier_stokes_elements.h"


namespace oomph
{


//========================================================================
/// Add element's contribution to the elemental 
/// residual vector and/or Jacobian matrix.
/// flag=1: compute both
/// flag=0: compute only residual vector
//=======================================================================
void RefineableAxisymmetricNavierStokesEquations::
fill_in_generic_residual_contribution_axi_nst(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              DenseMatrix<double> &mass_matrix,
                                              unsigned flag)
{
 //The dimension is actually two
 unsigned DIM=2;
 
//Find out how many nodes there are
unsigned n_node = nnode();
   
//Find out how many pressure dofs there are
unsigned n_pres = npres_axi_nst();
 
// Get the local indices of the nodal coordinates
 unsigned u_nodal_index[3];
 for(unsigned i=0;i<3;++i) {u_nodal_index[i] = u_index_axi_nst(i);}

// Which nodal value represents the pressure? (Negative if pressure
// is not based on nodal interpolation).
int p_index = this->p_nodal_index_axi_nst();

// Local array of booleans that are true if the l-th pressure value is
// hanging (avoid repeated virtual function calls)
 bool pressure_dof_is_hanging[n_pres];
 //If the pressure is stored at a node
 if(p_index >= 0)
  {
   //Read out whether the pressure is hanging
   for(unsigned l=0;l<n_pres;++l)
    {
     pressure_dof_is_hanging[l] = 
      pressure_node_pt(l)->is_hanging(p_index);
    }
  }
 //Otherwise the pressure is not stored at a node and so cannot hang
 else
  {
   for(unsigned l=0;l<n_pres;++l)
    {pressure_dof_is_hanging[l] = false;}
  }
   
    
//Set up memory for the shape and test functions
Shape psif(n_node), testf(n_node);
DShape dpsifdx(n_node,DIM), dtestfdx(n_node,DIM);
    
    
//Set up memory for pressure shape and test functions
Shape psip(n_pres), testp(n_pres);
    
//Set the value of Nintpt
unsigned Nintpt = integral_pt()->nweight();
    
//Set the Vector to hold local coordinates
Vector<double> s(DIM);
    
//Get Physical Variables from Element
//Reynolds number must be multiplied by the density ratio
double scaled_re  = re()*density_ratio();
double scaled_re_st = re_st()*density_ratio();
double scaled_re_inv_fr = re_invfr()*density_ratio();
double scaled_re_inv_ro = re_invro()*density_ratio();
double visc_ratio = viscosity_ratio(); // hierher -- rewrite and
                                       // make consistent with 
                                       //non-refineable version
Vector<double> G = g();
    
//Integers to store the local equation and unknown numbers
int local_eqn=0, local_unknown=0;

//Local storage for pointers to hang info objects
 HangInfo *hang_info_pt=0, *hang_info2_pt=0; 
    
//Loop over the integration points
for(unsigned ipt=0;ipt<Nintpt;ipt++)
{
      
 //Assign values of s
 for(unsigned i=0;i<DIM;i++) {s[i] = integral_pt()->knot(ipt,i);}
     
 //Get the integral weight
 double w = integral_pt()->weight(ipt);
      
 //Call the derivatives of the shape and test functions
 double J = 
  dshape_and_dtest_eulerian_at_knot_axi_nst(ipt,psif,dpsifdx,testf,dtestfdx);
   
 //Call the pressure shape and test functions
 pshape_axi_nst(s,psip,testp);
      
 //Premultiply the weights and the Jacobian
 double W = w*J;
      
 //Calculate local values of the pressure and velocity components
 //--------------------------------------------------------------
 double interpolated_p=0.0;
 Vector<double> interpolated_u(DIM+1,0.0);
 Vector<double> interpolated_x(DIM,0.0);
 Vector<double> mesh_veloc(DIM,0.0);
 Vector<double> dudt(DIM+1,0.0);
 DenseMatrix<double> interpolated_dudx(DIM+1,DIM,0.0);

 //Calculate pressure
 for(unsigned l=0;l<n_pres;l++) {interpolated_p += p_axi_nst(l)*psip[l];}
      
      
 //Calculate velocities and derivatives
      
 // Loop over nodes
 for(unsigned l=0;l<n_node;l++) 
  {
   //Cache the shape function
   const double psif_ = psif(l);
   //Loop over directions
   for(unsigned i=0;i<DIM;i++)
    {
     interpolated_x[i] += nodal_position(l,i)*psif_;
     //mesh_veloc[i] +=dnodal_position_dt(l,i)*psif(l);
    }

   for(unsigned i=0;i<DIM+1;i++)
    {
     const double u_value = nodal_value(l,u_nodal_index[i]);
     interpolated_u[i] += u_value*psif_;
     dudt[i]+=du_dt_axi_nst(l,i)*psif_;
     //Loop over derivative directions for gradients
     for(unsigned j=0;j<DIM;j++)
      {                               
       interpolated_dudx(i,j) += u_value*dpsifdx(l,j);
      }
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
       mesh_veloc[i]  += this->dnodal_position_dt(l,i)*psif(l);
      }
    }
  }
       

 //Get the user-defined body force terms
 Vector<double> body_force(DIM+1);
 get_body_force(time(),ipt,interpolated_x,body_force);
      
 //Get the user-defined source function
 double source=get_source_fct(time(),ipt,interpolated_x);

 // r is the first postition component
 double r = interpolated_x[0];
  
 //MOMENTUM EQUATIONS
 //==================
 //Number of master nodes and storage for the weight of the shape function
 unsigned n_master=1; double hang_weight=1.0;
     
 //Loop over the nodes for the test functions/equations
 //----------------------------------------------------
 for(unsigned l=0;l<n_node;l++)
  {
   //Local boolean that indicates the hanging status of the node
   bool is_node_hanging = node_pt(l)->is_hanging();

   //If the node is hanging
   if(is_node_hanging)
    {
     //Get the hanging pointer
     hang_info_pt = node_pt(l)->hanging_pt();
     //Read out number of master nodes from hanging data
     n_master = hang_info_pt->nmaster();
    }
   //Otherwise the node is its own master
   else
    {
     n_master = 1;
    }
       
   //Loop over the master nodes
   for(unsigned m=0;m<n_master;m++)
    {
     // Loop over velocity components for equations
     for(unsigned i=0;i<DIM+1;i++)
      {
       //Get the equation number
       //If the node is hanging
       if(is_node_hanging)
        {
         //Get the equation number from the master node
         local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                          u_nodal_index[i]);
         //Get the hang weight from the master node
         hang_weight = hang_info_pt->master_weight(m);
        }
       //If the node is not hanging
       else
        {
         // Local equation number
         local_eqn = this->nodal_local_eqn(l,u_nodal_index[i]);
         // Node contributes with full weight
         hang_weight = 1.0;
        }
       
       //If it's not a boundary condition...
       if(local_eqn>= 0)
        {
         // initialise for residual calculation
         double sum=0.0;
             
         switch (i)
          {
           // RADIAL MOMENTUM EQUATION
          case 0:
           //Add the user-defined body force terms
           sum += 
            r*body_force[0]*testf[l]*W*hang_weight;
               
           //Add the gravitational body force term
           sum += r*scaled_re_inv_fr*testf[l]*G[0]*W*hang_weight;
               
           //Add the pressure gradient term
           sum  += 
            interpolated_p*(testf[l] + r*dtestfdx(l,0))*W*hang_weight;
               
           //Add in the stress tensor terms
           //The viscosity ratio needs to go in here to ensure
           //continuity of normal stress is satisfied even in flows
           //with zero pressure gradient!
           sum -= visc_ratio*
            r*(1.0+Gamma[0])*interpolated_dudx(0,0)*dtestfdx(l,0)*W
            *hang_weight;
               
           sum -= visc_ratio*r*
            (interpolated_dudx(0,1) + Gamma[0]*interpolated_dudx(1,0))*
            dtestfdx(l,1)*W*hang_weight;
               
           sum -= 
            visc_ratio*(1.0 + Gamma[0])*interpolated_u[0]
            *testf[l]*W*hang_weight/r;
               
           //Add in the inertial terms
           //du/dt term
           sum -= scaled_re_st*r*dudt[0]*testf[l]*W*hang_weight;
               
           //Convective terms
           sum -= 
            scaled_re*(r*interpolated_u[0]*interpolated_dudx(0,0) 
                - interpolated_u[2]*interpolated_u[2] 
                + r*interpolated_u[1]*interpolated_dudx(0,1))*testf[l]*W*
            hang_weight;
               
           //Mesh velocity terms
           if(!ALE_is_disabled)
            {
             for(unsigned k=0;k<2;k++)
              {
               sum += 
                scaled_re_st*r*mesh_veloc[k]*interpolated_dudx(0,k)*testf[l]*W*
                hang_weight;
              }
            }

           //Coriolis term
           sum +=  
            2.0*r*scaled_re_inv_ro*interpolated_u[2]*testf[l]*W*hang_weight;

           break;
                
           // AXIAL MOMENTUM EQUATION
          case 1:
           //If it's not a boundary condition
           //Add the user-defined body force terms
           //Remember to multiply by the density ratio!
           sum += r*body_force[1]*testf[l]*W*hang_weight;
                
           //Add the gravitational body force term
           sum += r*scaled_re_inv_fr*testf[l]*G[1]*W*hang_weight;
                
           //Add the pressure gradient term
           sum  += r*interpolated_p*dtestfdx(l,1)*W*hang_weight;
                
           //Add in the stress tensor terms
           //The viscosity ratio needs to go in here to ensure
           //continuity of normal stress is satisfied even in flows
           //with zero pressure gradient!
           sum -= visc_ratio*
            r*(interpolated_dudx(1,0) + 
               Gamma[1]*interpolated_dudx(0,1))*dtestfdx(l,0)*W*
            hang_weight;
                
           sum -= visc_ratio*r*
            (1.0 + Gamma[1])*interpolated_dudx(1,1)*dtestfdx(l,1)*W*
            hang_weight;
                
           //Add in the inertial terms
           //du/dt term
           sum -= scaled_re_st*r*dudt[1]*testf[l]*W*hang_weight;
                
           //Convective terms
           sum -= 
            scaled_re*(r*interpolated_u[0]*interpolated_dudx(1,0) 
                + r*interpolated_u[1]*interpolated_dudx(1,1))*
            testf[l]*W*hang_weight;
                
           //Mesh velocity terms
           if(!ALE_is_disabled)
            {
             for(unsigned k=0;k<2;k++)
              {
               sum += 
                scaled_re_st*r*mesh_veloc[k]*interpolated_dudx(1,k)*testf[l]*W
                *hang_weight;
              }
            }
           break;
                
           // AZIMUTHAL MOMENTUM EQUATION
          case 2:
           //Add the user-defined body force terms
           //Remember to multiply by the density ratio!
           sum += r*body_force[2]*testf[l]*W*hang_weight;
                
           //Add the gravitational body force term
           sum += r*scaled_re_inv_fr*testf[l]*G[2]*W*hang_weight;
                
           //There is NO pressure gradient term
                
           //Add in the stress tensor terms
           //The viscosity ratio needs to go in here to ensure
           //continuity of normal stress is satisfied even in flows
           //with zero pressure gradient!
           sum -= visc_ratio*
            (r*interpolated_dudx(2,0) - Gamma[0]*interpolated_u[2])
            *dtestfdx(l,0)*W*hang_weight;
                
           sum -= visc_ratio*r*interpolated_dudx(2,1)*
            dtestfdx(l,1)*W*hang_weight;
                
           sum -= visc_ratio*
            ((interpolated_u[2]/r) - Gamma[0]*interpolated_dudx(2,0))*
            testf[l]*W*hang_weight;
                
                
           //Add in the inertial terms
           //du/dt term
           sum -= scaled_re_st*r*dudt[2]*testf[l]*W*hang_weight;
                
           //Convective terms
           sum -= 
            scaled_re*(r*interpolated_u[0]*interpolated_dudx(2,0)
                + interpolated_u[0]*interpolated_u[2]
                + r*interpolated_u[1]*interpolated_dudx(2,1))*testf[l]*W
            *hang_weight;
                
           //Mesh velocity terms
           if(!ALE_is_disabled)
            {
             for(unsigned k=0;k<2;k++)
              {
               sum += scaled_re_st*r*mesh_veloc[k]*
                interpolated_dudx(2,k)*testf[l]*W*hang_weight;
              }
            }

           //Coriolis term
           sum -= 
            2.0*r*scaled_re_inv_ro*interpolated_u[0]*testf[l]*W*hang_weight;

           break;
          }

         // Add contribution
         residuals[local_eqn] += sum;
             
         //CALCULATE THE JACOBIAN
         if (flag)
          {
           //Number of master nodes and weights
           unsigned n_master2=1; double hang_weight2=1.0;
           //Loop over the velocity nodes for columns
           for(unsigned l2=0;l2<n_node;l2++)
            {
             //Local boolean for hanging status
             bool is_node2_hanging = node_pt(l2)->is_hanging();
                 
             //If the node is hanging
             if(is_node2_hanging)
              {
               hang_info2_pt = node_pt(l2)->hanging_pt();
               //Read out number of master nodes from hanging data
               n_master2 = hang_info2_pt->nmaster();
              }
             //Otherwise the node is its own master
             else
              {
               n_master2 = 1;
              }
                 
             //Loop over the master nodes
             for(unsigned m2=0;m2<n_master2;m2++)
              {
               //Loop over the velocity components
               for(unsigned i2=0;i2<DIM+1;i2++)
                {
                 //Get the number of the unknown
                 //If the node is hanging
                 if(is_node2_hanging)
                  {
                   //Get the equation number from the master node
                   local_unknown = 
                    this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                         u_nodal_index[i2]);
                   hang_weight2 = hang_info2_pt->master_weight(m2);
                  }
                 else
                  {
                   local_unknown = this->nodal_local_eqn(l2,u_nodal_index[i2]);
                   hang_weight2 = 1.0;
                  }
                     
                 // If the unknown is non-zero
                 if(local_unknown >= 0)
                  {
                   //Different results for i and i2
                   switch (i)
                    {
                     // RADIAL MOMENTUM EQUATION
                    case 0:
                     switch (i2)
                      {
                       // radial component
                      case 0:
                                              
                       //Add the mass matrix entries
                       if(flag==2)
                        {
                         mass_matrix(local_eqn,local_unknown) +=
                          scaled_re_st*r*psif[l2]*testf[l]*W
                          *hang_weight*hang_weight2;
                        }

                       //Add contribution to the Jacobian matrix
                       jacobian(local_eqn,local_unknown)
                        -= visc_ratio*r*(1.0+Gamma[0])
                        *dpsifdx(l2,0)*dtestfdx(l,0)*W*
                        hang_weight*hang_weight2;
                           
                       jacobian(local_eqn,local_unknown)
                        -= visc_ratio*r*dpsifdx(l2,1)*
                        dtestfdx(l,1)*W*hang_weight*hang_weight2;
                           
                       jacobian(local_eqn,local_unknown)
                        -= visc_ratio*(1.0 + Gamma[0])*psif[l2]*
                        testf[l]*W*hang_weight*hang_weight2/r;

                       //Add in the inertial terms
                       //du/dt term
                       jacobian(local_eqn,local_unknown)
                        -= scaled_re_st*r*node_pt(l2)->time_stepper_pt()->
                        weight(1,0)*psif[l2]*testf[l]*W
                        *hang_weight*hang_weight2;
                           
                       //Convective terms
                       jacobian(local_eqn,local_unknown) -=
                        scaled_re*(r*psif[l2]*interpolated_dudx(0,0) 
                            + r*interpolated_u[0]*dpsifdx(l2,0)
                            + r*interpolated_u[1]*dpsifdx(l2,1))
                        *testf[l]*W*hang_weight*hang_weight2;
                           
                       //Mesh velocity terms
                       if(!ALE_is_disabled)
                        {
                         for(unsigned k=0;k<2;k++)
                          {
                           jacobian(local_eqn,local_unknown)
                            += scaled_re_st*r*mesh_veloc[k]*dpsifdx(l2,k)
                            *testf[l]*W*hang_weight*hang_weight2;
                          }
                        }
                       break;
                            
                       // axial component
                      case 1:
                       jacobian(local_eqn,local_unknown) -=
                        visc_ratio*r*Gamma[0]*dpsifdx(l2,0)
                        *dtestfdx(l,1)*W*hang_weight*hang_weight2;
                           
                       //Convective terms
                       jacobian(local_eqn,local_unknown) -= 
                        scaled_re*r*psif[l2]*interpolated_dudx(0,1)*
                        testf[l]*W*hang_weight*hang_weight2;
                       break;
                           
                       // azimuthal component
                      case 2:
                       //Convective terms
                       jacobian(local_eqn,local_unknown) -= 
                        - scaled_re*2.0*interpolated_u[2]*psif[l2]
                        *testf[l]*W*hang_weight*hang_weight2;

                       //Coriolis terms
                       jacobian(local_eqn,local_unknown) +=
                        2.0*r*scaled_re_inv_ro*psif[l2]*testf[l]*W
                        *hang_weight*hang_weight2;

                       break;
                      } /*End of contribution radial momentum eqn*/
                     break;
                         
                     // AXIAL MOMENTUM EQUATION
                    case 1:
                     switch (i2)
                      {
                       // radial component
                      case 0:
                       //Add in the stress tensor terms
                       //The viscosity ratio needs to go in here to ensure
                       //continuity of normal stress is satisfied even 
                       //in flows with zero pressure gradient!
                       jacobian(local_eqn,local_unknown) -= 
                        visc_ratio*r*Gamma[1]*dpsifdx(l2,1)*
                        dtestfdx(l,0)*W*hang_weight*hang_weight2;
                            
                       //Convective terms
                       jacobian(local_eqn,local_unknown) -= 
                        scaled_re*r*psif[l2]*interpolated_dudx(1,0)*testf[l]*W
                        *hang_weight*hang_weight2;
                       break;
                           
                       // axial component
                      case 1:
                                                  
                       //Add the mass matrix terms
                       if(flag==2)
                        {
                         mass_matrix(local_eqn,local_unknown) +=
                          scaled_re_st*r*psif[l2]*testf[l]*W
                          *hang_weight*hang_weight2;
                        }

                       //Add in the stress tensor terms
                       //The viscosity ratio needs to go in here to ensure
                       //continuity of normal stress is satisfied even in 
                       //flows with zero pressure gradient!
                       jacobian(local_eqn,local_unknown) -= 
                        visc_ratio*r*dpsifdx(l2,0)*dtestfdx(l,0)*W
                        *hang_weight*hang_weight2;
                           
                       jacobian(local_eqn,local_unknown) -= 
                        visc_ratio*r*(1.0 + Gamma[1])*dpsifdx(l2,1)*
                        dtestfdx(l,1)*W*hang_weight*hang_weight2;

                       //Add in the inertial terms
                       //du/dt term
                       jacobian(local_eqn,local_unknown) -= 
                        scaled_re_st*r*node_pt(l2)->time_stepper_pt()->
                        weight(1,0)*psif[l2]*testf[l]*W
                        *hang_weight*hang_weight2;
                            
                       //Convective terms
                       jacobian(local_eqn,local_unknown) -= 
                        scaled_re*(r*interpolated_u[0]*dpsifdx(l2,0) 
                            + r*psif[l2]*interpolated_dudx(1,1)
                            + r*interpolated_u[1]*dpsifdx(l2,1))
                        *testf[l]*W*hang_weight*hang_weight2;
                           
                       //Mesh velocity terms
                       if(!ALE_is_disabled)
                        {
                         for(unsigned k=0;k<2;k++)
                          {
                           jacobian(local_eqn,local_unknown) 
                            += scaled_re_st*r*mesh_veloc[k]*dpsifdx(l2,k)
                            *testf[l]*W*hang_weight*hang_weight2;
                          }
                        }
                       break;
                            
                       // azimuthal component
                      case 2:
                       //There are no azimithal terms in the axial
                       //momentum equation
                       break;
                      }
                     break;
                          
                     // AZIMUTHAL MOMENTUM EQUATION
                    case 2:
                     switch (i2)
                      {
                       // radial component
                      case 0:
                       //Convective terms
                       jacobian(local_eqn,local_unknown) -= 
                        scaled_re*(r*psif[l2]*interpolated_dudx(2,0)
                            + psif[l2]*interpolated_u[2])*testf[l]
                        *W*hang_weight*hang_weight2;

                       //Coriolis term
                       jacobian(local_eqn,local_unknown) -=
                        2.0*r*scaled_re_inv_ro*psif[l2]*testf[l]*W
                        *hang_weight*hang_weight2;

                       break;
                            
                       // axial component
                      case 1:
                       //Convective terms
                       jacobian(local_eqn,local_unknown) -= 
                        scaled_re*r*psif[l2]*interpolated_dudx(2,1)*testf[l]
                        *W*hang_weight*hang_weight2;
                       break;
                            
                       // azimuthal component
                      case 2:
                                                   
                       //Add the mass matrix terms
                       if(flag==2)
                        {
                         mass_matrix(local_eqn,local_unknown) +=
                          scaled_re_st*r*psif[l2]*testf[l]*W
                          *hang_weight*hang_weight2;
                        }
                       
                       //Add in the stress tensor terms
                       //The viscosity ratio needs to go in here to ensure
                       //continuity of normal stress is satisfied even in 
                       //flows with zero pressure gradient!
                       jacobian(local_eqn,local_unknown) -= 
                        visc_ratio*
                        (r*dpsifdx(l2,0) - Gamma[0]*psif[l2])
                        *dtestfdx(l,0)*W*hang_weight*hang_weight2;
                            
                       jacobian(local_eqn,local_unknown) -= 
                        visc_ratio*r*dpsifdx(l2,1)*dtestfdx(l,1)
                        *W*hang_weight*hang_weight2;
                            
                       jacobian(local_eqn,local_unknown) -= 
                        visc_ratio*
                        ((psif[l2]/r) - Gamma[0]*dpsifdx(l2,0))
                        *testf[l]*W*hang_weight*hang_weight2;

                       //Add in the inertial terms
                       //du/dt term
                       jacobian(local_eqn,local_unknown) -= 
                        scaled_re_st*r*node_pt(l2)
                        ->time_stepper_pt()->weight(1,0)*
                        psif[l2]*testf[l]*W*hang_weight*hang_weight2;
                            
                       //Convective terms
                       jacobian(local_eqn,local_unknown) -= 
                        scaled_re*(r*interpolated_u[0]*dpsifdx(l2,0)
                            + interpolated_u[0]*psif[l2]
                            + r*interpolated_u[1]*dpsifdx(l2,1))
                        *testf[l]*W*hang_weight*hang_weight2;
                            
                       //Mesh velocity terms
                       if(!ALE_is_disabled)
                        {
                         for(unsigned k=0;k<2;k++)
                          {
                           jacobian(local_eqn,local_unknown) 
                            += scaled_re_st*r*mesh_veloc[k]*dpsifdx(l2,k)
                            *testf[l]*W*hang_weight*hang_weight2;
                          }
                        }
                       break;
                      }
                     break;
                    }
                  }
                }
              }
            } //End of loop over the nodes
                
                
           //Loop over the pressure shape functions
           for(unsigned l2=0;l2<n_pres;l2++)
            {
             //If the pressure dof is hanging
             if(pressure_dof_is_hanging[l2])
              {
               // Pressure dof is hanging so it must be nodal-based
               hang_info2_pt = 
                pressure_node_pt(l2)->hanging_pt(p_index);

               //Get the number of master nodes from the pressure node
               n_master2 = hang_info2_pt->nmaster();
              }
             //Otherwise the node is its own master
             else
              {
               n_master2 = 1;
              }
                 
             //Loop over the master nodes
             for(unsigned m2=0;m2<n_master2;m2++)
              {
               //Get the number of the unknown
               //If the pressure dof is hanging
               if(pressure_dof_is_hanging[l2])
                {
                 //Get the unknown from the node
                 local_unknown = 
                  local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                 p_index);
                 //Get the weight from the hanging object
                 hang_weight2 = hang_info2_pt->master_weight(m2);
                }
               else
                {
                 local_unknown = p_local_eqn(l2);
                 hang_weight2 = 1.0;
                }
                   
               //If the unknown is not pinned
               if(local_unknown >= 0)
                {
                 //Add in contributions to different equations
                 switch (i)
                  {
                   // RADIAL MOMENTUM EQUATION
                  case 0:
                   jacobian(local_eqn,local_unknown)
                    += psip[l2]*(testf[l] + r*dtestfdx(l,0))*W*
                    hang_weight*hang_weight2;
                   break;
                       
                   // AXIAL MOMENTUM EQUATION
                  case 1:
                   jacobian(local_eqn,local_unknown)
                    += r*psip[l2]*dtestfdx(l,1)*W*hang_weight*hang_weight2;
                   break;
                       
                   // AZIMUTHAL MOMENTUM EQUATION
                  case 2:
                   break;
                  }
                }
              }
            } //End of loop over pressure dofs
          }// End of Jacobian calculation
        } //End of if not boundary condition statement
      } //End of loop over velocity components 
    } //End of loop over master nodes
  } //End of loop over nodes
                          
        
 //CONTINUITY EQUATION
 //===================
        
 //Loop over the pressure shape functions
 for(unsigned l=0;l<n_pres;l++)
  {
   //If the pressure dof is hanging
   if(pressure_dof_is_hanging[l])
    {
     // Pressure dof is hanging so it must be nodal-based
     hang_info_pt = pressure_node_pt(l)->hanging_pt(p_index);
     //Get the number of master nodes from the pressure node
     n_master = hang_info_pt->nmaster();
    }
   //Otherwise the node is its own master
   else
    {
     n_master = 1;
    }
         
   //Loop over the master nodes
   for(unsigned m=0;m<n_master;m++)
    {
     //Get the number of the unknown
     //If the pressure dof is hanging
     if(pressure_dof_is_hanging[l])
      {
       local_eqn = local_hang_eqn(hang_info_pt->master_node_pt(m),p_index);
       hang_weight = hang_info_pt->master_weight(m);
      }
     else
      {
       local_eqn = p_local_eqn(l);
       hang_weight = 1.0;
      }
           
     //If the equation is not pinned
     if(local_eqn >= 0)
      {
       // Source term
       residuals[local_eqn] -= r*source*testp[l]*W*hang_weight;
              
       //Gradient terms
       residuals[local_eqn] += 
        (interpolated_u[0] + r*interpolated_dudx(0,0) 
         + r*interpolated_dudx(1,1))*testp[l]*W*hang_weight;
              
       //CALCULATE THE JACOBIAN
       //======================
       if (flag)
        {
         unsigned n_master2=1; double hang_weight2=1.0;
         //Loop over the velocity nodes for columns
         for(unsigned l2=0;l2<n_node;l2++)
          {
           //Local boolean to indicate whether the node is hanging
           bool is_node2_hanging = node_pt(l2)->is_hanging();
                 
           //If the node is hanging
           if(is_node2_hanging)
            {
             hang_info2_pt = node_pt(l2)->hanging_pt();
             //Read out number of master nodes from hanging data
             n_master2 = hang_info2_pt->nmaster();
            }
           //Otherwise the node is its own master
           else
            {
             n_master2 = 1;
            }
                 
           //Loop over the master nodes
           for(unsigned m2=0;m2<n_master2;m2++)
            {
             //Loop over the velocity components
             for(unsigned i2=0;i2<DIM+1;i2++)
              {
               //Get the number of the unknown
               //If the node is hanging
               if(is_node2_hanging)
                {
                 //Get the equation number from the master node
                 local_unknown = 
                  local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                 u_nodal_index[i2]);
                 hang_weight2 = hang_info2_pt->master_weight(m2);
                }
               else
                {
                 local_unknown = this->nodal_local_eqn(l2,u_nodal_index[i2]);
                 hang_weight2 = 1.0;
                }
                     
               //If the unknown is not pinned
               if(local_unknown >= 0)
                {
                 switch (i2)
                  {
                   // radial component
                  case 0:
                   jacobian(local_eqn,local_unknown) +=
                    (psif[l2] + r*dpsifdx(l2,0))*testp[l]*W*
                    hang_weight*hang_weight2;
                   break;

                   // axial component
                  case 1:
                   jacobian(local_eqn,local_unknown) +=
                    r*dpsifdx(l2,1)*testp[l]*W*hang_weight*hang_weight2;
                   break;

                   // azimuthal component
                  case 2:
                   break;
                  }
                }
              }
            }
          } //End of loop over nodes

         //NO PRESSURE CONTRIBUTIONS TO CONTINUITY EQUATION
        } //End of Jacobian calculation
      }
    }
  } //End of loop over pressure nodes
       
} // end of loop over integration points


}

}
