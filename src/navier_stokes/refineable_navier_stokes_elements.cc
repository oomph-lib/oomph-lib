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
 #include "refineable_navier_stokes_elements.h"



namespace oomph
{


//========================================================================
/// Add element's contribution to the elemental 
/// residual vector and/or Jacobian matrix.
/// flag=1: compute both
/// flag=0: compute only residual vector
//========================================================================
template<unsigned DIM>
 void RefineableNavierStokesEquations<DIM>::
fill_in_generic_residual_contribution_nst(Vector<double> &residuals, 
                                          DenseMatrix<double> &jacobian, 
                                          DenseMatrix<double> &mass_matrix,
                                          unsigned flag)
{
//Find out how many nodes there are
unsigned n_node = nnode();

//Find out how many pressure dofs there are
unsigned n_pres = this->npres_nst();

// Get the indices at which the velocity components are stored
 unsigned u_nodal_index[DIM];
 for(unsigned i=0;i<DIM;i++) {u_nodal_index[i] = this->u_index_nst(i);}
 
// Which nodal value represents the pressure? (Negative if pressure
// is not based on nodal interpolation).
int p_index = this->p_nodal_index_nst();

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

//Set the value of n_intpt
unsigned n_intpt = integral_pt()->nweight();

//Set the Vector to hold local coordinates
Vector<double> s(DIM);

//Get Physical Variables from Element
//Reynolds number must be multiplied by the density ratio
double scaled_re = this->re()*this->density_ratio();
double scaled_re_st = this->re_st()*this->density_ratio();
double scaled_re_inv_fr = this->re_invfr()*this->density_ratio();
double visc_ratio = this->viscosity_ratio();
Vector<double> G = this->g();

//Integers that store the local equations and unknowns
int local_eqn=0, local_unknown=0;

//Pointers to hang info objects
 HangInfo *hang_info_pt=0, *hang_info2_pt=0;

//Local boolean for ALE (or not)
 bool ALE_is_disabled_flag = this->ALE_is_disabled;

//Loop over the integration points
for(unsigned ipt=0;ipt<n_intpt;ipt++)
{
 
 //Assign values of s
 for(unsigned i=0;i<DIM;i++) {s[i] = integral_pt()->knot(ipt,i);}
 
 //Get the integral weight
 double w = integral_pt()->weight(ipt);
 
 //Call the derivatives of the shape and test functions
 double J = 
  this->dshape_and_dtest_eulerian_at_knot_nst(ipt,psif,dpsifdx,testf,dtestfdx);
 
 //Call the pressure shape and test functions
 this->pshape_nst(s,psip,testp);
 
 //Premultiply the weights and the Jacobian
 double W = w*J;
 
 //Calculate local values of the pressure and velocity components
 //--------------------------------------------------------------
 double interpolated_p=0.0;
 Vector<double> interpolated_u(DIM,0.0);
 Vector<double> interpolated_x(DIM,0.0);
 Vector<double> mesh_veloc(DIM,0.0);
 Vector<double> dudt(DIM,0.0);
 DenseMatrix<double> interpolated_dudx(DIM,DIM,0.0);
 
 //Calculate pressure
 for(unsigned l=0;l<n_pres;l++) {interpolated_p += this->p_nst(l)*psip[l];}
 
 
 //Calculate velocities and derivatives
 
 // Loop over nodes
 for(unsigned l=0;l<n_node;l++) 
  {   
   //Loop over directions
   for(unsigned i=0;i<DIM;i++)
    {
     //Get the nodal value
     double u_value = this->nodal_value(l,u_nodal_index[i]);
     interpolated_u[i] += u_value*psif[l];
     interpolated_x[i] += this->nodal_position(l,i)*psif[l];
     dudt[i] += this->du_dt_nst(l,i)*psif[l];
     
     //Loop over derivative directions for velocity gradients
     for(unsigned j=0;j<DIM;j++)
      {                               
       interpolated_dudx(i,j) += u_value*dpsifdx(l,j);
      }
    }
  }       
 
 if (!ALE_is_disabled_flag)
  {
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++)
    {
     //Loop over directions
     for(unsigned i=0;i<DIM;i++)
      {
       mesh_veloc[i] += this->dnodal_position_dt(l,i)*psif[l];
      }
    }       
  }

 //Get the user-defined body force terms
 Vector<double> body_force(DIM);
 this->get_body_force_nst(time(),ipt,s,interpolated_x,body_force);
 
 //Get the user-defined source function
 double source = this->get_source_nst(time(),ipt,interpolated_x);
 
 //MOMENTUM EQUATIONS
 //==================

 //Number of master nodes and storage for the weight of the shape function
 unsigned n_master=1; double hang_weight=1.0;

 //Loop over the nodes for the test functions/equations
 //----------------------------------------------------
 for(unsigned l=0;l<n_node;l++)
  {
   //Local boolean to indicate whether the node is hanging
   bool is_node_hanging = node_pt(l)->is_hanging();

   //If the node is hanging
   if(is_node_hanging)
    {
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
     for(unsigned i=0;i<DIM;i++)
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
         //Temporary variable to hold the residuals
         double sum=0.0;
         
         //Add the user-defined body force terms
         sum += body_force[i];
         
         //Add the gravitational body force term
         sum += scaled_re_inv_fr*G[i];
         
         //Add in the inertial term
         sum -= scaled_re_st*dudt[i];
         
         //Convective terms, including mesh velocity
         for(unsigned k=0;k<DIM;k++)
          {
           double tmp=scaled_re*interpolated_u[k];
           if (!ALE_is_disabled_flag) 
            {tmp -= scaled_re_st*mesh_veloc[k];}
           sum -= tmp*interpolated_dudx(i,k);
          }
         
         //Add the pressure gradient term
         sum = (sum*testf[l] + interpolated_p*dtestfdx(l,i))*W*hang_weight;
         
         //Add in the stress tensor terms
         //The viscosity ratio needs to go in here to ensure
         //continuity of normal stress is satisfied even in flows
         //with zero pressure gradient!
         for(unsigned k=0;k<DIM;k++)
          {
           sum -= visc_ratio*
            (interpolated_dudx(i,k) + this->Gamma[i]*interpolated_dudx(k,i))
            *dtestfdx(l,k)*W*hang_weight;
          }
         
         // Add contribution
         residuals[local_eqn] += sum;
         
         //CALCULATE THE JACOBIAN
         if(flag)
          {
           //Number of master nodes and weights
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
               for(unsigned i2=0;i2<DIM;i2++)
                {
                 //Get the number of the unknown
                 //If the node is hanging
                 if(is_node2_hanging)
                  {
                   //Get the equation number from the master node
                   local_unknown = 
                    this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                         u_nodal_index[i2]);
                   //Get the hang weights
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
                   //Add contribution to Elemental Matrix
                   jacobian(local_eqn,local_unknown)
                    -= visc_ratio*this->Gamma[i]*dpsifdx(l2,i)*
                    dtestfdx(l,i2)*W*hang_weight*hang_weight2;
                   
                   //Now add in the inertial terms
                   jacobian(local_eqn,local_unknown)
                    -= scaled_re*psif[l2]*interpolated_dudx(i,i2)*testf[l]*W*
                    hang_weight*hang_weight2;
                   
                   //Extra diagonal components if i2=i
                   if(i2 == i)
                    {
                     //Mass matrix entries
                     //Again note the positive sign because the mass
                     //matrix is taken on the other side of the equation
                     if(flag==2)
                      {
                       mass_matrix(local_eqn,local_unknown) +=
                        scaled_re_st*psif[l2]*testf[l]*W*
                        hang_weight*hang_weight2;
                      }
                        
                     // du/dt term
                     jacobian(local_eqn,local_unknown)
                      -= scaled_re_st*
                      node_pt(l2)->time_stepper_pt()->weight(1,0)*
                      psif[l2]*testf[l]*W*hang_weight*hang_weight2;

                     //Extra advective terms
                     for(unsigned k=0;k<DIM;k++)
                      {
                       double tmp=scaled_re*interpolated_u[k];
                       if (!ALE_is_disabled_flag) 
                        {tmp -= scaled_re_st*mesh_veloc[k];}
                       
                       jacobian(local_eqn,local_unknown) -=
                        tmp*dpsifdx(l2,k)*testf[l]*W*hang_weight*hang_weight2;
                      }

                     // Extra viscous terms
                     for(unsigned k=0;k<DIM;k++)
                      {
                       jacobian(local_eqn,local_unknown)
                        -= visc_ratio*dpsifdx(l2,k)*
                        dtestfdx(l,k)*W*hang_weight*hang_weight2;
                      }
                     
                    }
                  }
                }
              }
            }
           
           //Loop over the pressure shape functions
           for(unsigned l2=0;l2<n_pres;l2++)
            {
             //If the pressure dof is hanging
             if(pressure_dof_is_hanging[l2])
              {
               hang_info2_pt = 
                this->pressure_node_pt(l2)->hanging_pt(p_index);
               // Pressure dof is hanging so it must be nodal-based
               //Get the number of master nodes from the pressure node
               n_master2 =  hang_info2_pt->nmaster();
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
                 //Get the unknown from the master node
                 local_unknown = 
                  this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                       p_index);
                 //Get the weight from the hanging object
                 hang_weight2 = hang_info2_pt->master_weight(m2);
                }
               else
                {
                 local_unknown = this->p_local_eqn(l2);
                 hang_weight2 = 1.0;
                }
               
               //If the unknown is not pinned
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown)
                  += psip[l2]*dtestfdx(l,i)*W*hang_weight*hang_weight2;
                }
              }
            }
           
          }// End of Jacobian calculation
         
        } //End of if not boundary condition statement
       
      } //End of loop over components of non-hanging node
     
    } //End of loop over master nodes

  } // End of loop over nodes for equations
 
 
 
 //CONTINUITY EQUATION
 //===================
 
 //Loop over the pressure shape functions
 for(unsigned l=0;l<n_pres;l++)
  {
   //If the pressure dof is hanging
   if(pressure_dof_is_hanging[l])
    {
     // Pressure dof is hanging so it must be nodal-based
     // Get the hang info object
     hang_info_pt = this->pressure_node_pt(l)->hanging_pt(p_index);

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
       //Get the local equation from the master node
       local_eqn = 
        this->local_hang_eqn(hang_info_pt->master_node_pt(m),p_index);
       //Get the weight for the node
       hang_weight = hang_info_pt->master_weight(m);
      }
     else
      {
       local_eqn = this->p_local_eqn(l);
       hang_weight = 1.0;
      }

     //If the equation is not pinned
     if(local_eqn >= 0)
      {
       // Source term
       residuals[local_eqn] -= source*testp[l]*W*hang_weight;
       
       // Loop over velocity components
       for(unsigned k=0;k<DIM;k++)
        {
         residuals[local_eqn] += interpolated_dudx(k,k)*testp[l]*W*hang_weight;
        }
       
       //CALCULATE THE JACOBIAN
       if(flag)
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
             for(unsigned i2=0;i2<DIM;i2++)
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
                 
               //If the unknown is not pinned
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown)
                  += dpsifdx(l2,i2)*testp[l]*W*hang_weight*hang_weight2;
                }
              }
            }
          }

         // NO PRESSURE CONTRIBUTION TO THE JACOBIAN

        } //End of jacobian calculation
      }
    }
  } //End of loop over pressure variables

} //End of loop over integration points
}




//======================================================================
/// Compute derivatives of elemental residual vector with respect
/// to nodal coordinates. 
/// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
/// Overloads the FD-based version in the FE base class.
//======================================================================
template <unsigned DIM>
void RefineableNavierStokesEquations<DIM>::get_dresidual_dnodal_coordinates(
 RankThreeTensor<double>&
 dresidual_dnodal_coordinates)
{
 // Return immediately if there are no dofs
 if (ndof()==0) return;

 //Find out how many nodes there are
 unsigned n_node = nnode();
 
 //Find out how many pressure dofs there are
 unsigned n_pres = this->npres_nst();

 //Find the indices at which the local velocities are stored
 unsigned u_nodal_index[DIM];
 for(unsigned i=0;i<DIM;i++) {u_nodal_index[i] = this->u_index_nst(i);}
 
 // Which nodal value represents the pressure? (Negative if pressure
 // is not based on nodal interpolation).
 int p_index = this->p_nodal_index_nst();
 
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
 DShape dpsifdx_pls(n_node,DIM), dtestfdx_pls(n_node,DIM);

 //Set up memory for pressure shape and test functions
 Shape psip(n_pres), testp(n_pres);

 // Get number of shape controlling nodes 
 unsigned n_shape_controlling_node=nshape_controlling_nodes();

 // Deriatives of shape fct derivatives w.r.t. nodal coords
 RankFourTensor<double> d_dpsifdx_dX(DIM,n_shape_controlling_node,n_node,DIM);
 RankFourTensor<double> d_dtestfdx_dX(DIM,n_shape_controlling_node,n_node,DIM);

 // Derivative of Jacobian of mapping w.r.t. to nodal coords
 DenseMatrix<double> dJ_dX(DIM,n_shape_controlling_node);

 // Derivatives of derivative of u w.r.t. nodal coords
 RankFourTensor<double> d_dudx_dX(DIM,n_shape_controlling_node,DIM,DIM);

 // Derivatives of nodal velocities w.r.t. nodal coords:
 // Assumption: Interaction only local via no-slip so 
 // X_ij only affects U_ij.
 DenseMatrix<double> d_U_dX(DIM,n_shape_controlling_node,0.0);

 //Number of integration points
 unsigned n_intpt = integral_pt()->nweight();
   
 //Set the Vector to hold local coordinates
 Vector<double> s(DIM);

 //Get Physical Variables from Element
 //Reynolds number must be multiplied by the density ratio
 double scaled_re = this->re()*this->density_ratio();
 double scaled_re_st = this->re_st()*this->density_ratio();
 double scaled_re_inv_fr = this->re_invfr()*this->density_ratio();
 double visc_ratio = this->viscosity_ratio();
 Vector<double> G = this->g();
 
 // FD step 
 double eps_fd=GeneralisedElement::Default_fd_jacobian_step;
   
 // Pre-compute derivatives of nodal velocities w.r.t. nodal coords:
 // Assumption: Interaction only local via no-slip so 
 // X_ij only affects U_ij.
 bool element_has_node_with_aux_node_update_fct=false;


 std::map<Node*,unsigned> local_shape_controlling_node_lookup= 
  shape_controlling_node_lookup();
 
 // FD loop over shape-controlling nodes
 for (std::map<Node*,unsigned>::iterator it=
       local_shape_controlling_node_lookup.begin();
      it!=local_shape_controlling_node_lookup.end();
      it++)
  {  
   // Get node
   Node* nod_pt=it->first;
   
   // Get its number
   unsigned jj=it->second;
   
   // Only compute if there's a node-update fct involved
   if (nod_pt->has_auxiliary_node_update_fct_pt())
    {
     element_has_node_with_aux_node_update_fct=true;
     
     // Current nodal velocity
     Vector<double> u_ref(DIM);
     for (unsigned i=0;i<DIM;i++)
      {
       u_ref[i]=*(nod_pt->value_pt(u_nodal_index[i]));
      }
     
     // FD
     for (unsigned ii=0;ii<DIM;ii++)
      {
       // Make backup
       double backup=nod_pt->x(ii);
       
       // Do FD step. No node update required as we're
       // attacking the coordinate directly...
       nod_pt->x(ii)+=eps_fd;
       
       // Do auxiliary node update (to apply no slip)
       nod_pt->perform_auxiliary_node_update_fct();
       
       // Evaluate
       d_U_dX(ii,jj)=(*(nod_pt->value_pt(u_nodal_index[ii]))-u_ref[ii])/eps_fd;
       
       // Reset 
       nod_pt->x(ii)=backup;
       
       // Do auxiliary node update (to apply no slip)
       nod_pt->perform_auxiliary_node_update_fct();
      }
    }
  }

 //Integers to store the local equation numbers
 int local_eqn=0;

 //Pointers to hang info object
 HangInfo *hang_info_pt=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Assign values of s
   for(unsigned i=0;i<DIM;i++) s[i] = integral_pt()->knot(ipt,i);

   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   //Call the derivatives of the shape and test functions
   double J = 
    this->dshape_and_dtest_eulerian_at_knot_nst(ipt,psif,dpsifdx,
                                                testf,dtestfdx);
   
   //Call the pressure shape and test functions
   this->pshape_nst(s,psip,testp);
   
   //Calculate local values of the pressure and velocity components
   //Allocate
   double interpolated_p=0.0;
   Vector<double> interpolated_u(DIM,0.0);
   Vector<double> interpolated_x(DIM,0.0);
   Vector<double> mesh_velocity(DIM,0.0);
   Vector<double> dudt(DIM,0.0);
   DenseMatrix<double> interpolated_dudx(DIM,DIM,0.0);    

   //Calculate pressure
   for(unsigned l=0;l<n_pres;l++) interpolated_p += this->p_nst(l)*psip[l];
   
   //Calculate velocities and derivatives:

   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Loop over directions
     for(unsigned i=0;i<DIM;i++)
      {
       //Get the nodal value
       double u_value = nodal_value(l,u_nodal_index[i]);
       interpolated_u[i] += u_value*psif[l];
       interpolated_x[i] += nodal_position(l,i)*psif[l];
       dudt[i] += this->du_dt_nst(l,i)*psif[l];
       
       //Loop over derivative directions
       for(unsigned j=0;j<DIM;j++)
        {                               
         interpolated_dudx(i,j) += u_value*dpsifdx(l,j);
        }
      }
    }

   if (!this->ALE_is_disabled)
    {
     // Loop over nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       //Loop over directions
       for(unsigned i=0;i<DIM;i++)
        {
         mesh_velocity[i] += this->dnodal_position_dt(l,i)*psif[l];
        }
      }
    }

   // Get weight of actual nodal position/value in computation of mesh
   // velocity from positional/value time stepper
   double pos_time_weight=node_pt(0)->position_time_stepper_pt()->weight(1,0);
   double val_time_weight=node_pt(0)->time_stepper_pt()->weight(1,0);

   //Get the user-defined body force terms
   Vector<double> body_force(DIM);
   this->get_body_force_nst(time(),ipt,s,interpolated_x,body_force);
   
   //Get the user-defined source function
   double source = this->get_source_nst(time(),ipt,interpolated_x);

   
   // FD loop over shape-controlling nodes
   for (std::map<Node*,unsigned>::iterator it=
         local_shape_controlling_node_lookup.begin();
        it!=local_shape_controlling_node_lookup.end();
        it++)
    {  
     // Get node
     Node* nod_pt=it->first;
     
     // Get its number
     unsigned jj=it->second;
     
     // Loop over coordinate directions
     for (unsigned ii=0;ii<DIM;ii++)
      {
       // Make backup
       double backup=nod_pt->x(ii);
       
       // Do FD step. No node update required as we're
       // attacking the coordinate directly...
       nod_pt->x(ii)+=eps_fd;
       
       //Call the derivatives of the shape and test functions
       //at advanced level
       double J_pls = 
        this->dshape_and_dtest_eulerian_at_knot_nst(ipt,psif,dpsifdx_pls,
                                                    testf,dtestfdx_pls);
       
       // Assign
       dJ_dX(ii,jj)=(J_pls-J)/eps_fd;
       for (unsigned i=0;i<DIM;i++)
        {
         for (unsigned j=0;j<n_node;j++)
          {
           d_dpsifdx_dX(ii,jj,j,i)=(dpsifdx_pls(j,i)-dpsifdx(j,i))/eps_fd;
           d_dtestfdx_dX(ii,jj,j,i)=(dtestfdx_pls(j,i)-dtestfdx(j,i))/eps_fd;
          }
        }

       // Shape deriv of du_i/dx_j
       for (unsigned i=0;i<DIM;i++)
        {
         for (unsigned j=0;j<DIM;j++)
          {
           double aux=0.0;
           for (unsigned j_nod=0;j_nod<n_node;j_nod++)
            {
             aux+=nodal_value(j_nod,u_nodal_index[i])*
              d_dpsifdx_dX(ii,jj,j_nod,j);
            }
           d_dudx_dX(ii,jj,i,j)=aux;
          }
        }

       // Reset coordinate. No node update required as we're
       // attacking the coordinate directly...
       nod_pt->x(ii)=backup;
      }
    }

   // Get gradient of body force function
   DenseMatrix<double> d_body_force_dx(DIM,DIM,0.0);
   this->get_body_force_gradient_nst(time(),ipt,s,
                                     interpolated_x, d_body_force_dx);

   // Get gradient of source function
   Vector<double> source_gradient(DIM,0.0);
   this->get_source_gradient_nst(time(),ipt,interpolated_x, source_gradient);


   // Assemble shape derivatives
   //---------------------------

   //MOMENTUM EQUATIONS
   //------------------
   
   //Number of master nodes and storage for the weight of the shape function
   unsigned n_master=1; double hang_weight=1.0;

   // Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {

     //Local boolean to indicate whether the node is hanging
     bool is_node_hanging = node_pt(l)->is_hanging();
     
     //If the node is hanging
     if(is_node_hanging)
      {
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
       
       // Loop over coordinate directions
       for(unsigned i=0;i<DIM;i++)
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
         
         /*IF it's not a boundary condition*/
         if(local_eqn >= 0)
          {
           // Loop over coordinate directions
           for (unsigned ii=0;ii<DIM;ii++)
            {              
             // Loop over nodes
             for (unsigned jj=0;jj<n_shape_controlling_node;jj++)
              {       
               // Residual x deriv of Jacobian
               //-----------------------------
               
               //Add the user-defined body force terms
               double sum = body_force[i]*testf[l];
               
               //Add the gravitational body force term
               sum += scaled_re_inv_fr*testf[l]*G[i];
               
               //Add the pressure gradient term
               sum  += interpolated_p*dtestfdx(l,i);
               
               //Add in the stress tensor terms
               //The viscosity ratio needs to go in here to ensure
               //continuity of normal stress is satisfied even in flows
               //with zero pressure gradient!
               for(unsigned k=0;k<DIM;k++)
                {
                 sum -= visc_ratio*
                  (interpolated_dudx(i,k) + 
                   this->Gamma[i]*interpolated_dudx(k,i))*dtestfdx(l,k);
                }
               
               //Add in the inertial terms
               //du/dt term
               sum -= scaled_re_st*dudt[i]*testf[l];
               
               
               //Convective terms, including mesh velocity
               for(unsigned k=0;k<DIM;k++)
                {
                 double tmp=scaled_re*interpolated_u[k];
                 if (!this->ALE_is_disabled)
                  {
                   tmp-=scaled_re_st*mesh_velocity[k];
                  }
                 sum -= tmp*interpolated_dudx(i,k)*testf[l];
                }
               
               // Multiply through by deriv of Jacobian and integration weight
               dresidual_dnodal_coordinates(local_eqn,ii,jj)+=
                sum*dJ_dX(ii,jj)*w*hang_weight;
             
               // Derivative of residual x Jacobian
               //----------------------------------
               
               // Body force
               sum=d_body_force_dx(i,ii)*psif(jj)*testf(l);
               
               // Pressure gradient term
               sum += interpolated_p*d_dtestfdx_dX(ii,jj,l,i);
               
               // Viscous term
               for (unsigned k=0;k<DIM;k++)
                {
                 sum -= visc_ratio*(
                  (interpolated_dudx(i,k)+
                   this->Gamma[i]*interpolated_dudx(k,i))
                  *d_dtestfdx_dX(ii,jj,l,k)+                
                  (d_dudx_dX(ii,jj,i,k) + 
                   this->Gamma[i]*d_dudx_dX(ii,jj,k,i))
                  *dtestfdx(l,k));
                }
               
               //Convective terms, including mesh velocity
               for(unsigned k=0;k<DIM;k++)
                {
                 double tmp=scaled_re*interpolated_u[k];
                 if (!this->ALE_is_disabled)
                  {
                   tmp-=scaled_re_st*mesh_velocity[k];
                  }
                 sum -= tmp*d_dudx_dX(ii,jj,i,k)*testf(l);
                }
               sum+=scaled_re_st*pos_time_weight*
                psif(jj)*interpolated_dudx(i,ii)*testf(l);
               
               // Multiply through by Jacobian and integration weight
               dresidual_dnodal_coordinates(local_eqn,ii,jj)+=
                sum*J*w*hang_weight;
              }
            }
             
           

           // Derivs w.r.t. to nodal velocities
           //----------------------------------
           if (element_has_node_with_aux_node_update_fct)
            {
             // Loop over local nodes
             for (unsigned jj_local=0;jj_local<n_node;jj_local++)
              {
               
               //Number of master nodes and storage for the weight of 
               // the shape function
               unsigned n_master2=1; 
               double hang_weight2=1.0;
               HangInfo* hang_info2_pt=0;
               
               //Local boolean to indicate whether the node is hanging
               bool is_node_hanging2 = node_pt(jj_local)->is_hanging();
               
               Node* actual_shape_controlling_node_pt=node_pt(jj_local);

               //If the node is hanging
               if(is_node_hanging2)
                {
                 hang_info2_pt = node_pt(jj_local)->hanging_pt();
                 
                 //Read out number of master nodes from hanging data
                 n_master2 = hang_info2_pt->nmaster();
                }
               //Otherwise the node is its own master
               else
                {
                 n_master2 = 1;
                }
               
               //Loop over the master nodes
               for(unsigned mm=0;mm<n_master2;mm++)
                {                 

                 if(is_node_hanging2)
                  {
                   actual_shape_controlling_node_pt=
                    hang_info2_pt->master_node_pt(mm);
                   hang_weight2 = hang_info2_pt->master_weight(mm);
                  }

                 // Find the corresponding number
                 unsigned jj=local_shape_controlling_node_lookup[
                  actual_shape_controlling_node_pt];

                 // Loop over coordinate directions
                 for(unsigned ii=0;ii<DIM;ii++)
                  {
                   double sum=
                    -visc_ratio*this->Gamma[i]*dpsifdx(jj_local,i)*
                    dtestfdx(l,ii)
                    -scaled_re*psif(jj_local)*interpolated_dudx(i,ii)*testf(l);
                   if (i==ii)
                    {
                     sum-=scaled_re_st*val_time_weight*psif(jj_local)*testf(l);
                     for (unsigned k=0;k<DIM;k++)
                      {
                       sum-=visc_ratio*dpsifdx(jj_local,k)*dtestfdx(l,k);
                       double tmp=scaled_re*interpolated_u[k];
                       if (!this->ALE_is_disabled)
                        {
                         tmp-=scaled_re_st*mesh_velocity[k];
                        }
                       sum-=tmp*dpsifdx(jj_local,k)*testf(l); 
                      }
                    }

                   dresidual_dnodal_coordinates(local_eqn,ii,jj)+=
                    sum*d_U_dX(ii,jj)*J*w*hang_weight*hang_weight2; 
                  }
                }
              }
            }


          } // local_eqn>=0
        }
      }
    }
   
   
   //CONTINUITY EQUATION
   //-------------------
   
   //Loop over the shape functions
   for(unsigned l=0;l<n_pres;l++)
    {
     
     //If the pressure dof is hanging
     if(pressure_dof_is_hanging[l])
      {
       // Pressure dof is hanging so it must be nodal-based
       // Get the hang info object
       hang_info_pt = this->pressure_node_pt(l)->hanging_pt(p_index);
       
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
         //Get the local equation from the master node
         local_eqn = 
          this->local_hang_eqn(hang_info_pt->master_node_pt(m),p_index);
         //Get the weight for the node
         hang_weight = hang_info_pt->master_weight(m);
        }
       else
        {
         local_eqn = this->p_local_eqn(l);
         hang_weight = 1.0;
        }
       
       //If not a boundary conditions
       if(local_eqn >= 0)
        {
         // Loop over coordinate directions
         for (unsigned ii=0;ii<DIM;ii++)
          {              
           // Loop over nodes
           for (unsigned jj=0;jj<n_shape_controlling_node;jj++)
            {       
             
             // Residual x deriv of Jacobian
             //-----------------------------
             
             // Source term
             double aux=-source;
             
             //Loop over velocity components
             for(unsigned k=0;k<DIM;k++)
              {
               aux += interpolated_dudx(k,k);
              }
             
             // Multiply through by deriv of Jacobian and integration weight
             dresidual_dnodal_coordinates(local_eqn,ii,jj)+=
              aux*dJ_dX(ii,jj)*testp[l]*w*hang_weight;


             // Derivative of residual x Jacobian
             //----------------------------------          
             
             //Loop over velocity components
             aux=-source_gradient[ii]*psif(jj);
             for(unsigned k=0;k<DIM;k++)
              {
               aux += d_dudx_dX(ii,jj,k,k);
              }
             // Multiply through by Jacobian and integration weight
             dresidual_dnodal_coordinates(local_eqn,ii,jj)+=
              aux*testp[l]*J*w*hang_weight;
            }
          }
        
         
         
         // Derivs w.r.t. to nodal velocities
         //----------------------------------
         if (element_has_node_with_aux_node_update_fct)
          {
           // Loop over local nodes
           for (unsigned jj_local=0;jj_local<n_node;jj_local++)
            {
             
             //Number of master nodes and storage for the weight of 
             // the shape function
             unsigned n_master2=1; 
             double hang_weight2=1.0;
             HangInfo* hang_info2_pt=0;
             
             //Local boolean to indicate whether the node is hanging
             bool is_node_hanging2 = node_pt(jj_local)->is_hanging();
             
             Node* actual_shape_controlling_node_pt=node_pt(jj_local);
             
             //If the node is hanging
             if(is_node_hanging2)
              {
               hang_info2_pt = node_pt(jj_local)->hanging_pt();
               
               //Read out number of master nodes from hanging data
               n_master2 = hang_info2_pt->nmaster();
              }
             //Otherwise the node is its own master
             else
              {
               n_master2 = 1;
              }
             
             //Loop over the master nodes
             for(unsigned mm=0;mm<n_master2;mm++)
              {                 
               
               if(is_node_hanging2)
                {
                 actual_shape_controlling_node_pt=
                  hang_info2_pt->master_node_pt(mm);
                 hang_weight2 = hang_info2_pt->master_weight(mm);
                }
               
               // Find the corresponding number
               unsigned jj=local_shape_controlling_node_lookup[
                actual_shape_controlling_node_pt];
               
               // Loop over coordinate directions
               for(unsigned ii=0;ii<DIM;ii++)
                {
                 double aux=dpsifdx(jj_local,ii)*testp(l);
                 dresidual_dnodal_coordinates(local_eqn,ii,jj)+=
                  aux*d_U_dX(ii,jj)*J*w*hang_weight*hang_weight2;
                }       
              }     
            }
          }
        }
      }
    }

 }// End of loop over integration points
}   



     

//====================================================================
//// Force build of templates
//====================================================================
template class RefineableNavierStokesEquations<2>;
template class RefineableNavierStokesEquations<3>;


}
