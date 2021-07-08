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
#include "refineable_unsteady_heat_elements.h"


namespace oomph
{


//========================================================================
/// Add element's contribution to the elemental 
/// residual vector and/or Jacobian matrix.
/// flag=1: compute both
/// flag=0: compute only residual vector
//========================================================================
template<unsigned DIM>
void RefineableUnsteadyHeatEquations<DIM>::
fill_in_generic_residual_contribution_ust_heat(Vector<double> &residuals, 
                                               DenseMatrix<double> &jacobian, 
                                               unsigned flag)
{

//Find out how many nodes there are in the element
unsigned n_node = nnode();

// Get continuous time from timestepper of first node
double time=node_pt(0)->time_stepper_pt()->time_pt()->time();

//Find the index at which the unknown is stored
 unsigned u_nodal_index = this->u_index_ust_heat(); 

//Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

//Set the value of n_intpt
unsigned n_intpt = integral_pt()->nweight();

//Set the Vector to hold local coordinates
Vector<double> s(DIM);

 //Get Alpha and beta parameters number
 double alpha_local = this->alpha();
 double beta_local = this->beta();

//Integers used to store the local equation number and local unknown
//indices for the residuals and jacobians
int local_eqn=0, local_unknown=0;

//Local storage for pointers to hang info objects
 HangInfo *hang_info_pt=0, *hang_info2_pt=0;

//Local variable to determine the ALE stuff
 bool ALE_is_disabled_flag = this->ALE_is_disabled; 

//Loop over the integration points
for(unsigned ipt=0;ipt<n_intpt;ipt++)
{
 
 //Assign values of s
 for(unsigned i=0;i<DIM;i++) s[i] = integral_pt()->knot(ipt,i);
 
 //Get the integral weight
 double w = integral_pt()->weight(ipt);
 
 //Call the derivatives of the shape and test functions
 double J = this->
  dshape_and_dtest_eulerian_at_knot_ust_heat(ipt,psi,dpsidx,test,dtestdx);
 
 //Premultiply the weights and the Jacobian
 double W = w*J;
 
 //Calculate local values of the function
 double dudt=0.0;
 double interpolated_u=0.0;
 
 //This needs to be a Vector to be ANSI C++, initialise to zero
 Vector<double> interpolated_x(DIM,0.0);
 Vector<double> interpolated_dudx(DIM,0.0);
 Vector<double> mesh_velocity(DIM,0.0);
 
 
 //Calculate function value and derivatives:
 //-----------------------------------------
 
 // Loop over nodes
 for(unsigned l=0;l<n_node;l++) 
  {
   //Get the value of the unknown at the node
   double u_value = this->nodal_value(l,u_nodal_index);
   interpolated_u += u_value*psi(l);
   dudt += this->du_dt_ust_heat(l)*psi(l);
   // Loop over directions
   for(unsigned j=0;j<DIM;j++)
    {
     interpolated_x[j] += nodal_position(l,j)*psi(l);
     interpolated_dudx[j] += u_value*dpsidx(l,j);
    }
  }

 if (!ALE_is_disabled_flag)
  {
   for(unsigned l=0;l<n_node;l++) 
    {
     // Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       mesh_velocity[j] += dnodal_position_dt(l,j)*psi(l);
      }
    }
  }

 //Get body force
 double source;
 this->get_source_ust_heat(time,ipt,interpolated_x,source);


 // Assemble residuals and Jacobian
 //================================
   
 // Loop over the nodes for the test functions 
 //-------------------------------------------
 for(unsigned l=0;l<n_node;l++)
  {
   //Local variables to store the number of master nodes and 
   //the weight associated with the shape function if the node is hanging
   unsigned n_master=1; double hang_weight=1.0;
   //Local bool (is the node hanging)
   bool is_node_hanging = this->node_pt(l)->is_hanging();


   //If the node is hanging, get the number of master nodes
   if(is_node_hanging)
    {
     hang_info_pt = this->node_pt(l)->hanging_pt();
     n_master = hang_info_pt->nmaster();
    }
   //Otherwise there is just one master node, the node itself
   else 
    {
     n_master = 1;
    }
     
   //Loop over the number of master nodes
   for(unsigned m=0;m<n_master;m++)
    {
     //Get the local equation number and hang_weight
     //If the node is hanging
     if(is_node_hanging)
      {
       //Read out the local equation from the master node
       local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                        u_nodal_index);
       //Read out the weight from the master node
       hang_weight = hang_info_pt->master_weight(m);
      }
     //If the node is not hanging
     else
      {
       //The local equation number comes from the node itself
       local_eqn = this->nodal_local_eqn(l,u_nodal_index);
       //The hang weight is one
       hang_weight = 1.0;
      }
       
     //If the nodal equation is not a boundary condition
     if(local_eqn >= 0)
      {
       // Add body force/source term and time derivative
       residuals[local_eqn]+=(alpha_local*dudt + source)*test(l)*W*hang_weight;
         
       // Mesh velocity and  Laplace operator itself
       for(unsigned k=0;k<DIM;k++)
        {
         double tmp=dtestdx(l,k)*beta_local;
         if (!ALE_is_disabled_flag) tmp -= alpha_local*mesh_velocity[k]*test(l);
         residuals[local_eqn] += interpolated_dudx[k]*tmp*W*hang_weight;
        }
         
       // Calculate the Jacobian
       if(flag)
        {
         //Local variables to store the number of master nodes
         //and the weights associated with each hanging node
         unsigned n_master2=1; double hang_weight2=1.0;
         //Loop over the nodes for the variables
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           //Local bool (is the node hanging)
           bool is_node2_hanging = this->node_pt(l2)->is_hanging();
           //If the node is hanging, get the number of master nodes
           if(is_node2_hanging)
            {
             hang_info2_pt = this->node_pt(l2)->hanging_pt();
             n_master2 = hang_info2_pt->nmaster();
            }
           //Otherwise there is one master node, the node itself
           else
            {
             n_master2 = 1;
            }
             
           //Loop over the master nodes
           for(unsigned m2=0;m2<n_master2;m2++)
            {
             //Get the local unknown and weight
             //If the node is hanging
             if(is_node2_hanging)
              {
               //Read out the local unknown from the master node
               local_unknown = 
                this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                     u_nodal_index);
               //Read out the hanging weight from the master node
               hang_weight2 = hang_info2_pt->master_weight(m2);
              }
             //If the node is not hanging
             else
              {
               //The local unknown number comes from the node
               local_unknown = this->nodal_local_eqn(l2,u_nodal_index);
               //The hang weight is one
               hang_weight2 = 1.0;
              }
               
             //If the unknown is not pinned
             if(local_unknown >= 0)
              {
               //Add contribution to Elemental Matrix
               // Mass matrix
               jacobian(local_eqn,local_unknown) 
                += alpha_local*test(l)*psi(l2)*
                this->node_pt(l2)->time_stepper_pt()->weight(1,0)
                *W*hang_weight*hang_weight2;

               // Laplace operator and mesh veloc bit
               for(unsigned k=0;k<DIM;k++)
                {
                 double tmp=dtestdx(l,k)*beta_local;
                 if(!ALE_is_disabled_flag)
                  {
                   tmp -= 
                    alpha_local*mesh_velocity[k]*test(l);
                  }
                 jacobian(local_eqn,local_unknown) += 
                  dpsidx(l2,k)*tmp*W*hang_weight*hang_weight2;
                }
              }
            } //End of loop over master nodes
          } 
        } //End of Jacobian calculation
      }
    } //End of loop over master nodes for residuals
  } //End of loop over nodes

} // End of loop over integration points


}

//====================================================================
// Force build of templates
//====================================================================
template class RefineableQUnsteadyHeatElement<2,2>;
template class RefineableQUnsteadyHeatElement<2,3>;
template class RefineableQUnsteadyHeatElement<2,4>;

template class RefineableQUnsteadyHeatElement<3,2>;
template class RefineableQUnsteadyHeatElement<3,3>;
template class RefineableQUnsteadyHeatElement<3,4>;

}
