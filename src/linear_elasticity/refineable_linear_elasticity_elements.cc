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
//Non-inline member functions and static member data for refineable solid 
//mechanics elements

#include "refineable_linear_elasticity_elements.h"

namespace oomph
{

//====================================================================
/// Residuals for Refineable QLinearElasticityElements
//====================================================================
template<unsigned DIM>
void RefineableLinearElasticityEquations<DIM>::
fill_in_generic_contribution_to_residuals_linear_elasticity(
 Vector<double> &residuals,DenseMatrix<double> &jacobian,unsigned flag)
{

 //Find out how many nodes there are
 unsigned n_node = this->nnode();

 //Find out how many positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();

 if(n_position_type != 1)
  {
   throw OomphLibError(
    "LinearElasticity is not yet implemented for more than one position type",
    "LinearElasticityEquationsBase<DIM>::fill_in_generic_contribution_to_residuals_lin_elast()",
    OOMPH_EXCEPTION_LOCATION);
  }

 //Throw and error if an elasticity tensor has not been set
 if(this->Elasticity_tensor_pt==0)
  {
   throw OomphLibError(
    "No elasticity tensor set",
    "LinearElasticityEquationsBase<DIM>::fill_in_generic_contribution_to_residuals_lin_elast()",
    OOMPH_EXCEPTION_LOCATION);
  }
 
 //Find the indices at which the local velocities are stored
 unsigned u_nodal_index[DIM];
 for(unsigned i=0;i<DIM;i++) 
  {u_nodal_index[i] = this->u_index_linear_elasticity(i);}
 
 // Timescale ratio (non-dim density)
// double Lambda_sq = this->lambda_sq();
  
 //Set up memory for the shape functions
 Shape psi(n_node);
 DShape dpsidx(n_node,DIM);

 //Set the value of Nintpt -- the number of integration points
 unsigned n_intpt = this->integral_pt()->nweight();
   
  //Set the vector to hold the local coordinates in the element
 Vector<double> s(DIM);

 //Integer to store the local equation number
 int local_eqn=0, local_unknown=0;

 // Local storage for pointers to hang_info objects
 HangInfo *hang_info_pt=0, *hang_info2_pt=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Assign the values of s
   for(unsigned i=0;i<DIM;++i) {s[i] = this->integral_pt()->knot(ipt,i);}

   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);

   //Call the derivatives of the shape functions (and get Jacobian)
   double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);
   
   //Storage for Eulerian coordinates (initialised to zero)
   Vector<double> interpolated_x(DIM,0.0);

   //Calculate interpolated values of the derivative of global position
   //wrt lagrangian coordinates
   DenseMatrix<double> interpolated_dudx(DIM,DIM,0.0);

   // Setup memory for accelerations (initialised to zero)
   //Vector<double> accel(DIM,0.0);


   //Calculate displacements and derivatives and lagrangian coordinates
   for(unsigned l=0;l<n_node;l++)
    {
     //Loop over displacement components (deformed position)
     for(unsigned i=0;i<DIM;i++)
      {
       //Calculate the Lagrangian coordinates and the accelerations
       interpolated_x[i] += this->nodal_position(l,i)*psi(l);

       // Only compute accelerations if inertia is switched on
       // otherwise the timestepper might not be able to 
       // work out dx_gen_dt(2,...)
       //if (this->unsteady())
       // {
       //  accel[i] += this->dnodal_position_dt(2,l,i)*psi(l);
       // }
       
       //Get the nodal displacements
       const double u_value = nodal_value(l,u_nodal_index[i]);

       //Loop over derivative directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_dudx(i,j) += u_value*dpsidx(l,j);
        }
      }
    }

   //Get body force at current time
   Vector<double> b(DIM);
   this->body_force(interpolated_x,b);

   //Premultiply the weights and the Jacobian
   double W = w*J; 

//=====EQUATIONS OF ELASTICITY FROM PRINCIPLE OF VIRTUAL DISPLACEMENTS========
   //Number of master nodes and storage for the weight of the shape function
   unsigned n_master=1; double hang_weight=1.0;
   
   //Loop over the test functions, nodes of the element
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
       //Loop over the displacement components
       for(unsigned a=0;a<DIM;a++)
        {
         //Get the equation number
         if(is_node_hanging)
          {
           //Get the equation number from the master node
           local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                            u_nodal_index[a]);
           //Get the hang weight from the master node
           hang_weight = hang_info_pt->master_weight(m);
          }
         //Otherwise the node is not hanging
         else
          {
           local_eqn = this->nodal_local_eqn(l,u_nodal_index[a]);
           hang_weight = 1.0;
          }
         
         /*IF it's not a boundary condition*/
         if(local_eqn >= 0)
          {
           // Acceleration and body force
           residuals[local_eqn] += 
            (/*Lambda_sq*accel[a]*/-b[a])*psi(l)*W*hang_weight;
           
           // Stress term
           for(unsigned b=0;b<DIM;b++)
            {
             for(unsigned c=0;c<DIM;c++)
              {
               for(unsigned d=0;d<DIM;d++)
                {
                 //Add the stress terms to the residuals
                 residuals[local_eqn] +=
                  this->E(a,b,c,d)*interpolated_dudx(c,d)*dpsidx(l,b)*W*
                  hang_weight;
                }
              }
            }
           
           //Jacobian entries
           if(flag)
            {
             //Number of master nodes and weights
             unsigned n_master2=1; double hang_weight2=1.0;
             //Loop over the displacement basis functions again
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
                 //Loop over the displacement components again
                 for(unsigned c2=0;c2<DIM;c2++)
                  {
                   
                   //Get the number of the unknown
                   //If the node is hanging
                   if(is_node2_hanging)
                    {
                     //Get the equation number from the master node
                     local_unknown = 
                      this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                           u_nodal_index[c2]);
                     //Get the hang weights
                     hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                   else
                    {
                     local_unknown = 
                      this->nodal_local_eqn(l2,u_nodal_index[c2]);
                     hang_weight2 = 1.0;
                    }
                   
                   //If it's not pinned
                   if(local_unknown >= 0)
                    {
                     for(unsigned b=0;b<DIM;b++)
                      {
                       for(unsigned d=0;d<DIM;d++)
                        {
                         //Add the contribution to the Jacobian matrix
                         jacobian(local_eqn,local_unknown) +=
                          this->E(a,b,c2,d)*dpsidx(l2,d)*dpsidx(l,b)*W
                          *hang_weight*hang_weight2;
                        }
                      }
                    } //End of if not boundary condition
                  }
                }
              }
            } //End of jacobian calculation
           
          } //End of if not boundary condition
        } //End of loop over coordinate directions
      }
    } //End of loop over shape functions
  } //End of loop over integration points
}


/*
//===========================================================================
/// Add element's contribution to the elemental 
/// residual vector and/or Jacobian matrix.
/// flag=1: compute both
/// flag=0: compute only residual vector
//==========================================================================
template<unsigned DIM>
void RefineableLinearElasticityEquationsWithPressure<DIM>::
fill_in_generic_residual_contribution_pvd_with_pressure(
 Vector<double> &residuals, DenseMatrix<double> &jacobian, 
 unsigned flag)
{
  // Simply set up initial condition?
 if (Solid_ic_pt!=0)
  {
   get_residuals_for_ic(residuals);
   return;
  }

 //Find out how many nodes there are
 unsigned n_node = nnode();

 //Find out how many position types of dof there are
 unsigned n_position_type = this->nnodal_position_type();

 //Find out how many pressure dofs there are
 unsigned n_solid_pres = this->nsolid_pres();

 //Find out the index of the solid dof
 int solid_p_index = this->solid_p_nodal_index();
 //Local array of booleans that is true if the l-th pressure value is hanging
 //This is an optimization because it avoids repeated virtual function calls 
 bool solid_pressure_dof_is_hanging[n_solid_pres];
 //If the solid pressure is stored at a node 
 if(solid_p_index >= 0)
  {
   //Read out whether the solid pressure is hanging
   for(unsigned l=0;l<n_solid_pres;++l)
    {
     solid_pressure_dof_is_hanging[l] = 
      solid_pressure_node_pt(l)->is_hanging(solid_p_index);
    }
  }
 //Otherwise the pressure is not stored at a node and so 
 //it cannot hang
 else
  {
   for(unsigned l=0;l<n_solid_pres;++l)
    {solid_pressure_dof_is_hanging[l] = false;}
  }

 //Integer for storage of local equation and unknown numbers
 int local_eqn=0, local_unknown=0;

 // Timescale ratio (non-dim density)
 double Lambda_sq = this->lambda_sq();

 //Set up memory for the shape functions
 Shape psi(n_node,n_position_type);
 DShape dpsidxi(n_node,n_position_type,DIM);
 
 //Set up memory for the pressure shape functions
 Shape psisp(n_solid_pres);

 //Set the value of n_intpt
 unsigned n_intpt = integral_pt()->nweight();

  //Set the vector to hold the local coordinates in the element
 Vector<double> s(DIM);

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Assign the values of s
   for(unsigned i=0;i<DIM;++i) {s[i] = integral_pt()->knot(ipt,i);}

   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Call the derivatives of the shape functions
   double J = dshape_lagrangian_at_knot(ipt,psi,dpsidxi);

   //Call the pressure shape functions
   this->solid_pshape_at_knot(ipt,psisp);

   //Storage for Lagrangian coordinates (initialised to zero)
   Vector<double> interpolated_xi(DIM,0.0);

   // Deformed tangent vectors
   DenseMatrix<double> interpolated_G(DIM,DIM,0.0);

   // Setup memory for accelerations
   Vector<double> accel(DIM,0.0);
      
   //Calculate displacements and derivatives and lagrangian coordinates
   for(unsigned l=0;l<n_node;l++)
    {
     //Loop over positional dofs
     for(unsigned k=0;k<n_position_type;k++)
      {
       //Loop over displacement components (deformed position)
       for(unsigned i=0;i<DIM;i++)
        {
         //Calculate the lagrangian coordinates and the accelerations
         interpolated_xi[i] += lagrangian_position_gen(l,k,i)*psi(l,k);

         // Only compute accelerations if inertia is switched on
         // otherwise the timestepper might not be able to 
         // work out dx_gen_dt(2,...)
         if (this->unsteady())
          {
           accel[i] += dnodal_position_gen_dt(2,l,k,i)*psi(l,k);
          }

         //Loop over derivative directions
         for(unsigned j=0;j<DIM;j++)
          {
           interpolated_G(j,i) += nodal_position_gen(l,k,i)*dpsidxi(l,k,j);
          }
        }
      }
    }

   //Get isotropic growth factor
   double gamma=1.0;
   this->get_isotropic_growth(s,interpolated_xi,gamma);

   //Get body force at current time
   Vector<double> b(DIM);
   this->body_force(interpolated_xi,b);

   // We use Cartesian coordinates as the reference coordinate
   // system. In this case the undeformed metric tensor is always
   // the identity matrix -- stretched by the isotropic growth
   double diag_entry=pow(gamma,2.0/double(DIM)); 
   DenseMatrix<double> g(DIM);
   for(unsigned i=0;i<DIM;i++)
    {
     for(unsigned j=0;j<DIM;j++)
      {
       if(i==j) {g(i,j) = diag_entry;}
       else {g(i,j) = 0.0;}
      }
    }
     
   //Premultiply the undeformed volume ratio (from the isotropic
   // growth), the weights and the Jacobian
   double W = gamma*w*J; 
     
   //Calculate the interpolated solid pressure
   double interpolated_solid_p=0.0;
   for(unsigned l=0;l<n_solid_pres;l++)
    {
     interpolated_solid_p += this->solid_p(l)*psisp[l];
    }

   //Declare and calculate the deformed metric tensor
   DenseMatrix<double> G(DIM);

   //Assign values of G 
   for(unsigned i=0;i<DIM;i++)
    {
     //Do upper half of matrix
     //Note that j must be signed here for the comparison test to work
     //Also i must be cast to an int
     for(int j=(DIM-1);j>=static_cast<int>(i);j--)
      {
       //Initialise G(i,j) to zero
       G(i,j) = 0.0;
       //Now calculate the dot product
       for(unsigned k=0;k<DIM;k++)
        {
         G(i,j) += interpolated_G(i,k)*interpolated_G(j,k);
        }
      }
     //Matrix is symmetric so just copy lower half
     for(int j=(i-1);j>=0;j--)
      {
       G(i,j) = G(j,i);
      }
    }

   //Now calculate the deviatoric stress and all pressure-related
   //quantitites
   DenseMatrix<double> sigma_dev(DIM), Gup(DIM);
   double detG = 0.0;
   double gen_dil=0.0;
   double inv_kappa=0.0;

   // Incompressible: Compute the deviatoric part of the stress tensor, the
   // contravariant deformed metric tensor and the determinant
   // of the deformed covariant metric tensor.
   if(this->Incompressible)
    {
     this->get_stress(g,G,sigma_dev,Gup,detG);
    }
   // Nearly incompressible: Compute the deviatoric part of the 
   // stress tensor, the contravariant deformed metric tensor,
   // the generalised dilatation and the inverse bulk modulus.
   else
    {
     this->get_stress(g,G,sigma_dev,Gup,gen_dil,inv_kappa);
    }


//=====EQUATIONS OF ELASTICITY FROM PRINCIPLE OF VIRTUAL DISPLACEMENTS========
       
   //Loop over the test functions, nodes of the element
   for(unsigned l=0;l<n_node;l++)
    {
     //Get pointer to local node l
     Node* local_node_pt = node_pt(l);

     //If the node is NOT a hanging node
     if(local_node_pt->is_hanging()==false)
      {
       //Loop over the types of dof
       for(unsigned k=0;k<n_position_type;k++)
        {
         //Loop over the displacement components
         for(unsigned i=0;i<DIM;i++)
          {
           local_eqn = position_local_eqn(l,k,i);
           //IF it's not a boundary condition
           if(local_eqn >= 0)
            {
             // Acceleration and body force
             residuals[local_eqn] += 
              (Lambda_sq*accel[i]-b[i])*psi(l,k)*W;
             
             // Stress term
             for(unsigned a=0;a<DIM;a++)
              {
               for(unsigned b=0;b<DIM;b++)
                {
                 //Add the "stress" terms to the residuals
                 residuals[local_eqn] += 
                  (sigma_dev(a,b) - interpolated_solid_p*Gup(a,b))
                  *interpolated_G(a,i)*dpsidxi(l,k,b)*W;
                }
              }
             
             //Can add in the pressure jacobian terms
             if(flag)
              {
               //Loop over the pressure nodes
               for(unsigned l2=0;l2<n_solid_pres;l2++)
                {
                 //If the pressure dof is not hanging
                 if(solid_pressure_dof_is_hanging[l2]==false)
                  {
                   //If it's not a boundary condition
                   local_unknown = this->solid_p_local_eqn(l2);
                   if(local_unknown >= 0)
                    {
                     //Add the pressure terms to the jacobian
                     for(unsigned a=0;a<DIM;a++)
                      {
                       for(unsigned b=0;b<DIM;b++)
                        {
                         jacobian(local_eqn,local_unknown) -=
                          psisp[l2]*Gup(a,b)*
                          interpolated_G(a,i)*dpsidxi(l,k,b)*W;
                        }
                      }
                    }
                  }
                 //Otherwise the pressure dof is hanging
                 else
                  {
                   //Get the HangInfo object associated with
                   //the hanging solid pressure
                   HangInfo* hang_info2_pt = 
                    solid_pressure_node_pt(l2)->hanging_pt(solid_p_index);
                   
                   //Loop over all the master nodes
                   unsigned n_master2 = hang_info2_pt->nmaster();
                   for(unsigned m2=0;m2<n_master2;m2++)
                    {
                     //Get the equation numbers at the master node
                     local_unknown = 
                      local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                     solid_p_index);
                     
                     //Find the hanging node weight at the node
                     double hang_weight2
                      = hang_info2_pt->master_weight(m2);
                     
                     //If it's not a boundary condition
                     if(local_unknown >= 0)
                      {
                       for(unsigned a=0;a<DIM;a++)
                        {
                         for(unsigned b=0;b<DIM;b++)
                          {
                           jacobian(local_eqn,local_unknown) -=
                            psisp[l2]*Gup(a,b)*
                            interpolated_G(a,i)*dpsidxi(l,k,b)*W*
                            hang_weight2;
                          }
                        }
                      }
                    } //End of loop over master weights
                  }
                } //End of loop over pressure dofs
              } //End of Jacobian terms

            } //End of if not boundary condition
          } //End of loop over coordinate directions
        } //End of loop over types of dof
       
      } //End of if not hanging
     //Hanging case
     else
      {
       //Load the local hang info object
       HangInfo* hang_info_pt = local_node_pt->hanging_pt();
       //Loop over the master nodes
       unsigned n_master = hang_info_pt->nmaster();
       for(unsigned m=0;m<n_master;m++)
        {
         //Find the equation numbers
         DenseMatrix<int> Position_local_eqn_at_node;
         Position_local_eqn_at_node = 
          local_position_hang_eqn(hang_info_pt->master_node_pt(m));
         //Find the hanging node weight
         double hang_weight = hang_info_pt->master_weight(m);
         
         //Loop of types of dofs
         for(unsigned k=0;k<n_position_type;k++)
          {
           //Loop over the displacement components
           for(unsigned i=0;i<DIM;i++)
            {
             local_eqn = Position_local_eqn_at_node(k,i);
             //IF it's not a boundary condition
             if(local_eqn >= 0)
              {
               // Acceleration and body force
               residuals[local_eqn] += 
                (Lambda_sq*accel[i]-b[i])*psi(l,k)*W*hang_weight;
               
               // Stress term
               for(unsigned a=0;a<DIM;a++)
                {
                 for(unsigned b=0;b<DIM;b++)
                  {
                   //Add the stress terms to the residuals
                   residuals[local_eqn] +=
                    (sigma_dev(a,b) - interpolated_solid_p*Gup(a,b))
                    *interpolated_G(a,i)*dpsidxi(l,k,b)*W*hang_weight;
                  }
                }
               
               //Can add in the pressure jacobian terms
               if(flag)
                {
                 //Loop over the pressure nodes
                 for(unsigned l2=0;l2<n_solid_pres;l2++)
                  {
                   //If the pressure dof is not hanging
                   if(solid_pressure_dof_is_hanging[l2]==false)
                    {
                     local_unknown = this->solid_p_local_eqn(l2);
                     //If it's not a boundary condition
                     if(local_unknown >= 0)
                      {
                       //Add the pressure terms to the jacobian
                       for(unsigned a=0;a<DIM;a++)
                        {
                         for(unsigned b=0;b<DIM;b++)
                          {
                           jacobian(local_eqn,local_unknown) -=
                            psisp[l2]*Gup(a,b)*
                            interpolated_G(a,i)*dpsidxi(l,k,b)*W
                            *hang_weight;
                          }
                        }
                      }
                    }
                   //Otherwise the pressure dof is hanging
                   else
                    {
                     //Get the HangInfo object associated with the
                     //hanging solid pressure
                     HangInfo* hang_info2_pt =
                      solid_pressure_node_pt(l2)->hanging_pt(solid_p_index);

                     //Loop over all the master nodes
                     unsigned n_master2 = hang_info2_pt->nmaster();
                     for(unsigned m2=0;m2<n_master2;m2++)
                      {
                       //Get the equation numbers at the master node
                       local_unknown = 
                        local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                       solid_p_index);
                       
                       //Find the hanging node weight at the node
                       double hang_weight2
                        = hang_info2_pt->master_weight(m2);
                       
                       //If it's not a boundary condition
                       if(local_unknown >= 0)
                        {
                         //Add the pressure terms to the jacobian
                         for(unsigned a=0;a<DIM;a++)
                          {
                           for(unsigned b=0;b<DIM;b++)
                            {
                             jacobian(local_eqn,local_unknown) -=
                              psisp[l2]*Gup(a,b)*
                              interpolated_G(a,i)*dpsidxi(l,k,b)*W
                              *hang_weight*hang_weight2;
                            }
                          }
                        }
                      } //End of loop over master weights
                    }
                  } //End of loop over pressure dofs
                } //End of Jacobian terms
               
              } //End of if not boundary condition
            }
          }
        }
      } //End of hanging node case
     
    } //End of loop over shape functions
     
   //==============CONSTRAINT EQUATIONS FOR PRESSURE=====================
       
   //Now loop over the pressure degrees of freedom
   for(unsigned l=0;l<n_solid_pres;l++)
    {
     //If the pressure dof is NOT hanging
     if(solid_pressure_dof_is_hanging[l]==false)
      {
       local_eqn = this->solid_p_local_eqn(l);
       // Pinned (unlikely, actually) or real dof?
       if(local_eqn >= 0)
        {
         //For true incompressibility we need to conserve volume
         //so the determinant of the deformed metric tensor
         //needs to be equal to that of the undeformed one, which
         //is equal to the volumetric growth factor
         if(this->Incompressible)
          {
           residuals[local_eqn] += (detG - gamma)*psisp[l]*W;
           
           //No Jacobian terms since the pressure does not feature
           //in the incompressibility constraint
          }
         //Nearly incompressible: (Neg.) pressure given by product of
         //bulk modulus and generalised dilatation
         else
          {
           residuals[local_eqn] += 
            (inv_kappa*interpolated_solid_p + gen_dil)*psisp[l]*W;
           
           //Add in the jacobian terms
           if(flag)
            {
             //Loop over the pressure nodes again
             for(unsigned l2=0;l2<n_solid_pres;l2++)
              {
               //If the pressure is NOT hanging
               if(solid_pressure_dof_is_hanging[l2]==false)
                {
                 local_unknown = this->solid_p_local_eqn(l2);
                 //If not pinnned 
                 if(local_unknown >= 0)
                  {
                   jacobian(local_eqn,local_unknown)
                    += inv_kappa*psisp[l2]*psisp[l]*W;
                  }
                }
               //Otherwise, if it's hanging
               else
                {
                 //Get the HangInfo object associated with the
                 //hanging solid pressure
                 //Note that the pressure is stored at
                 //the index solid_p_index at the node
                 HangInfo* hang_info2_pt = 
                  solid_pressure_node_pt(l2)->hanging_pt(solid_p_index);
 
                 //Loop over all the master nodes
                 unsigned n_master2 = hang_info2_pt->nmaster();
                 for(unsigned m2=0;m2<n_master2;m2++)
                  {
                   //Get the equation numbers at the master node
                   local_unknown = 
                    local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                   solid_p_index);
                                      
                   //Find the hanging node weight at the node
                   double hang_weight2
                    = hang_info2_pt->master_weight(m2);
 
                   //If it's not a boundary condition
                   if(local_unknown >= 0)
                    {
                     jacobian(local_eqn,local_unknown)
                      += inv_kappa*psisp[l2]*psisp[l]*W*hang_weight2;
                    }
                  }
                } //End of hanging case
              } //End of loop over pressure dofs
            } //End of Jacobian
          } //End of nearly incompressible case
        } //End of if not boundary condition
      }
     //Otherwise the pressure is hanging
     else
      {
       //Get a pointer to the HangInfo object associated with the
       //solid pressure (stored at solid_p_index)
       HangInfo* hang_info_pt = 
        solid_pressure_node_pt(l)->hanging_pt(solid_p_index);

       //Loop over all the master nodes
       //Note that the pressure is stored at the inded solid_p_index 
       unsigned n_master = hang_info_pt->nmaster();
       for(unsigned m=0;m<n_master;m++)
        {
         //Get the equation numbers at the master node
         local_eqn = local_hang_eqn(hang_info_pt->master_node_pt(m),
                                    solid_p_index);
         
         
         //Find the hanging node weight at the node
         double hang_weight = hang_info_pt->master_weight(m);
         
         // Pinned (unlikely, actually) or real dof?
         if(local_eqn >= 0)
          {
           //For true incompressibility we need to conserve volume
           //so the determinant of the deformed metric tensor
           //needs to be equal to that of the undeformed one, which
           //is equal to the volumetric growth factor
           if(this->Incompressible)
            {
             residuals[local_eqn] += 
              (detG - gamma)*psisp[l]*W*hang_weight;
             
             //No Jacobian terms since the pressure does not feature
             //in the incompressibility constraint
            }
           //Nearly incompressible: (Neg.) pressure given by product of
           //bulk modulus and generalised dilatation
           else
            {
             residuals[local_eqn] += 
              (inv_kappa*interpolated_solid_p + gen_dil)
              *psisp[l]*W*hang_weight;
             
             //Add in the jacobian terms
             if(flag)
              {
               //Loop over the pressure nodes again
               for(unsigned l2=0;l2<n_solid_pres;l2++)
                {
                 //If the pressure is NOT hanging
                 if(solid_pressure_dof_is_hanging[l2]==false)
                  {
                   //If not pinnned
                   local_unknown = this->solid_p_local_eqn(l2);
                   if(local_unknown >= 0)
                    {
                     jacobian(local_eqn,local_unknown)
                      += inv_kappa*psisp[l2]*psisp[l]*W*hang_weight;
                    }
                  }
                 //Otherwise, if it's hanging
                 else
                  {
                   //Get pointer to hang info object
                   //Note that the pressure is stored at
                   //the index solid_p_index
                   HangInfo* hang_info2_pt = 
                    solid_pressure_node_pt(l2)->hanging_pt(solid_p_index);
                   
                   //Loop over all the master nodes
                   unsigned n_master2 = hang_info2_pt->nmaster();
                   for(unsigned m2=0;m2<n_master2;m2++)
                    {
                     //Get the equation numbers at the master node
                     local_unknown = 
                      local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                     solid_p_index);
                     
                     //Find the hanging node weight at the node
                     double hang_weight2
                      = hang_info2_pt->master_weight(m2);
                    
                     //If it's not a boundary condition
                     if(local_unknown >= 0)
                      {
                       jacobian(local_eqn,local_unknown)
                        += inv_kappa*psisp[l2]*psisp[l]*W*
                        hang_weight*hang_weight2;
                      }
                    }
                  } //End of hanging case
                } //End of loop over pressure dofs
              }//End of Jacobian
            } //End of nearly incompressible case
          } //End of if not boundary condition
        } //End of loop over master nodes
      } //End of hanging case

    } //End of loop over pressure dofs

  } //End of loop over integration points
}
*/   



//====================================================================
/// Forcing building of required templates
//====================================================================
template class RefineableLinearElasticityEquations<2>;
template class RefineableLinearElasticityEquations<3>;

//template class RefineableLinearElasticityEquationsWithPressure<2>;
//template class RefineableLinearElasticityEquationsWithPressure<3>;

}
