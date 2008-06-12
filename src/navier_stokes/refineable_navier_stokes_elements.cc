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
 
//  //Initialise everything to zero
//  for(unsigned i=0;i<DIM;i++)
//   {
//    for(unsigned j=0;j<DIM;j++) {interpolated_dudx(i,j) = 0.0;}
//   }
 
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
 this->get_body_force_nst(time(),s,interpolated_x,body_force);
 
 //Get the user-defined source function
 double source = this->get_source_nst(time(),interpolated_x);
 
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
//            sum -= (scaled_re*interpolated_u[k] - scaled_re_st*mesh_veloc[k])*
//             interpolated_dudx(i,k);
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

//                        jacobian(local_eqn,local_unknown) -=
//                         (scaled_re*interpolated_u[k] -
//                          scaled_re_st*mesh_veloc[k])*
//                         dpsifdx(l2,k)*testf[l]*W*hang_weight*hang_weight2;
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
     //Get the hang info object
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

//====================================================================
//// Force build of templates
//====================================================================
template class RefineableNavierStokesEquations<2>;
template class RefineableNavierStokesEquations<3>;


}
