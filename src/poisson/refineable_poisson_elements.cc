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
#include "refineable_poisson_elements.h"


namespace oomph
{


//========================================================================
/// Add element's contribution to the elemental 
/// residual vector and/or Jacobian matrix.
/// flag=1: compute both
/// flag=0: compute only residual vector
//========================================================================
template<unsigned DIM>
void RefineablePoissonEquations<DIM>::
fill_in_generic_residual_contribution_poisson(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              const unsigned& flag)
{

//Find out how many nodes there are in the element
unsigned n_node = nnode();

//Set up memory for the shape and test functions
Shape psi(n_node), test(n_node);
DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

//Set the value of n_intpt
unsigned n_intpt = integral_pt()->nweight();

//The local index at which the poisson variable is stored
 unsigned u_nodal_index = this->u_index_poisson(); 

//Integers to store the local equation and unknown numbers
int local_eqn=0, local_unknown=0;

// Local storage for pointers to hang_info objects
 HangInfo *hang_info_pt=0, *hang_info2_pt=0;

//Loop over the integration points
for(unsigned ipt=0;ipt<n_intpt;ipt++)
{
 //Get the integral weight
 double w = integral_pt()->weight(ipt);
 
 //Call the derivatives of the shape and test functions
 double J = 
  this->dshape_and_dtest_eulerian_at_knot_poisson(ipt,psi,dpsidx,test,dtestdx);
 
 //Premultiply the weights and the Jacobian
 double W = w*J;
 
 // Position and gradient
 Vector<double> interpolated_x(DIM,0.0);
 Vector<double> interpolated_dudx(DIM,0.0);
 
 //Calculate function value and derivatives:
 //-----------------------------------------
 
 // Loop over nodes
 for(unsigned l=0;l<n_node;l++) 
  {
   //Get the poisson value from the node 
   //(hanging-ness will be taken into account
   double u_value = this->nodal_value(l,u_nodal_index);

   // Loop over directions
   for(unsigned j=0;j<DIM;j++)
    {
     interpolated_x[j] += nodal_position(l,j)*psi(l);
     interpolated_dudx[j] += u_value*dpsidx(l,j);
    }
  }
 
 //Get body force
 double source;
 this->get_source_poisson(ipt,interpolated_x,source);
 
 
 // Assemble residuals and Jacobian
 
 // Loop over the nodes for the test functions 
 for(unsigned l=0;l<n_node;l++)
  {
   //Local variables used to store the number of master nodes and the
   //weight associated with the shape function if the node is hanging
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
   
   //Loop over the master nodes
   for(unsigned m=0;m<n_master;m++)
    {
     //Get the local equation number and hang_weight
     //If the node is hanging
     if(is_node_hanging)
      {
       //Read out the local equation number from the m-th master node
       local_eqn =  this->local_hang_eqn(hang_info_pt->master_node_pt(m),
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
       // Add body force/source term here 
       residuals[local_eqn] += source*test(l)*W*hang_weight;
       
       // The Poisson bit itself
       for(unsigned k=0;k<DIM;k++)
        {
         residuals[local_eqn] += 
          interpolated_dudx[k]*dtestdx(l,k)*W*hang_weight;
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
               for(unsigned i=0;i<DIM;i++)
                {
                 jacobian(local_eqn,local_unknown) += 
                  dpsidx(l2,i)*dtestdx(l,i)*W*hang_weight*hang_weight2;
                }
              }
            } //End of loop over master nodes
          } //End of loop over nodes
        } //End of Jacobian calculation
       
      } //End of case when residual equation is not pinned
    } //End of loop over master nodes for residual
  } //End of loop over nodes
 
} // End of loop over integration points
}


//======================================================================
/// Compute derivatives of elemental residual vector with respect
/// to nodal coordinates. 
/// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
/// Overloads the FD-based version in the FE base class.
//======================================================================
template <unsigned DIM>
void RefineablePoissonEquations<DIM>::get_dresidual_dnodal_coordinates(
 RankThreeTensor<double>&
 dresidual_dnodal_coordinates)
{

 //Find out how many nodes there are
 const unsigned n_node = nnode();

 //Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);
 DShape dpsidx_pls(n_node,DIM), dtestdx_pls(n_node,DIM);

 // Get number of shape controlling nodes 
 unsigned n_shape_controlling_node=nshape_controlling_nodes();

 // Deriatives of shape fct derivatives w.r.t. nodal coords
 RankFourTensor<double> d_dpsidx_dX(DIM,n_shape_controlling_node,n_node,DIM);
 RankFourTensor<double> d_dtestdx_dX(DIM,n_shape_controlling_node,n_node,DIM);

 // Derivative of Jacobian of mapping w.r.t. to nodal coords
 DenseMatrix<double> dJ_dX(DIM,n_shape_controlling_node);

 // Derivatives of derivative of u w.r.t. nodal coords
 RankThreeTensor<double> d_dudx_dX(DIM,n_shape_controlling_node,DIM);

 // Gradient of source fct
 Vector<double> d_source_dx(DIM);

 //Index at which the poisson unknown is stored
 const unsigned u_nodal_index = this->u_index_poisson();
 
 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();

 //Integers to store the local equation number
 int local_eqn=0;

 // Local storage for pointers to hang_info object
 HangInfo *hang_info_pt=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Call the derivatives of the shape and test functions
   double J = this->dshape_and_dtest_eulerian_at_knot_poisson(ipt,psi,dpsidx,
                                                              test,dtestdx);
       
   //Calculate local values 
   //Allocate and initialise to zero
   Vector<double> interpolated_x(DIM,0.0);
   Vector<double> interpolated_dudx(DIM,0.0);
   
   //Calculate function value and derivatives:
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Get the nodal value of the Poisson unknown
     double u_value = nodal_value(l,u_nodal_index);
     // Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       interpolated_x[j] += nodal_position(l,j)*psi(l);
       interpolated_dudx[j] += u_value*dpsidx(l,j);
      }
    }

   //Get source function
   //-------------------
   double source;
   this->get_source_poisson(ipt,interpolated_x,source);

   // FD step 
   double eps_fd=GeneralisedElement::Default_fd_jacobian_step;
   
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
        this->dshape_and_dtest_eulerian_at_knot_poisson(ipt,psi,dpsidx_pls,
                                                        test,dtestdx_pls);
       
       // Assign
       dJ_dX(ii,jj)=(J_pls-J)/eps_fd;
       for (unsigned i=0;i<DIM;i++)
        {
         for (unsigned j=0;j<n_node;j++)
          {
           d_dpsidx_dX(ii,jj,j,i)=(dpsidx_pls(j,i)-dpsidx(j,i))/eps_fd;
           d_dtestdx_dX(ii,jj,j,i)=(dtestdx_pls(j,i)-dtestdx(j,i))/eps_fd;
          }
        }

       // Shape deriv of du/dx_i
       for (unsigned i=0;i<DIM;i++)
        {
         double aux=0.0;
         for (unsigned j_nod=0;j_nod<n_node;j_nod++)
          {
           aux+=nodal_value(j_nod,u_nodal_index)*
            d_dpsidx_dX(ii,jj,j_nod,i);
          }
         d_dudx_dX(ii,jj,i)=aux;
        }
  
       // Reset coordinate. No node update required as we're
       // attacking the coordinate directly...
       nod_pt->x(ii)=backup;
      }
    }
   
   // Get gradient of source function
   this->get_source_gradient_poisson(ipt, interpolated_x, d_source_dx);

   
   // Assemble shape derivatives
   //---------------------------
       
   // Loop over the nodes for the test functions 
   for(unsigned l=0;l<n_node;l++)
    {
     //Local variables used to store the number of master nodes and the
     //weight associated with the shape function if the node is hanging
     unsigned n_master=1; 
     double hang_weight=1.0;
     
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
     
     //Loop over the master nodes
     for(unsigned m=0;m<n_master;m++)
      {
       //Get the local equation number and hang_weight
       //If the node is hanging
       if(is_node_hanging)
        {
         //Read out the local equation number from the m-th master node
         local_eqn =  this->local_hang_eqn(hang_info_pt->master_node_pt(m),
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
         // Loop over coordinate directions
         for (unsigned ii=0;ii<DIM;ii++)
          {              
           // Loop over shape controlling nodes
           for (unsigned jj=0;jj<n_shape_controlling_node;jj++)
            {       
             double sum=source*psi(l)*dJ_dX(ii,jj)+
              d_source_dx[ii]*psi(l)*psi(jj)*J;
             
             for (unsigned k=0;k<DIM;k++)
              {
               sum+=interpolated_dudx[k]*(dtestdx(l,k)*dJ_dX(ii,jj)+
                                          d_dtestdx_dX(ii,jj,l,k)*J)
                + d_dudx_dX(ii,jj,k)*dtestdx(l,k)*J;
              }
             
             // Multiply through by integration weight
             dresidual_dnodal_coordinates(local_eqn,ii,jj)+=sum*w*hang_weight;
            }
          }
        }
      }
    }
   
  } // End of loop over integration points
 


}   



//====================================================================
// Force build of templates
//====================================================================
template class RefineableQPoissonElement<2,2>;
template class RefineableQPoissonElement<2,3>;
template class RefineableQPoissonElement<2,4>;

template class RefineableQPoissonElement<3,2>;
template class RefineableQPoissonElement<3,3>;
template class RefineableQPoissonElement<3,4>;

}
