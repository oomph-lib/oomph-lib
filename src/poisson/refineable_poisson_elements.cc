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
/// Validate against exact flux.
/// Flux is provided via function pointer.
/// Plot error at a given number of plot points.
//========================================================================
template <unsigned DIM>
void RefineablePoissonEquations<DIM>::compute_exact_Z2_error(
 std::ostream &outfile, 
 FiniteElement::SteadyExactSolutionFctPt exact_flux_pt,
 double& error, double& norm)
{ 
 // Initialise
 error=0.0;
 norm=0.0;
 
 // Vector of local coordinates
 Vector<double> s(DIM);
 
 // Vector of global coordinates
 Vector<double> x(DIM);
 
 // Find out how many nodes there are in the element
 const unsigned n_node = nnode();
 
 // Allocate storage for shape functions
 Shape psi(n_node);
 
 // Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();
  
 // Tecplot
 outfile << "ZONE" << std::endl;
 
 // Exact flux Vector (size is DIM)
 Vector<double> exact_flux(DIM);
 
 // Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   
   // Assign values of s
   for(unsigned i=0;i<DIM;i++)
    {
     s[i] = integral_pt()->knot(ipt,i);
    }
   
   // Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   // Get jacobian of mapping
   double J = J_eulerian(s);
   
   // Premultiply the weights and the Jacobian
   double W = w*J;
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   // Allocate storage for FE flux Vector
   Vector<double> fe_flux(DIM);

   // Get FE flux as Vector
   get_Z2_flux(s,fe_flux);
   
   // Get exact flux at this point
   (*exact_flux_pt)(x,exact_flux);
   
   // Output x,y,...
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }

   // Output exact flux
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << exact_flux[i] << " ";
    }

   // Output error
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << exact_flux[i]-fe_flux[i] << " ";
    }
   outfile << std::endl;

   // Add to RMS error:
   double sum=0.0;
   double sum2=0.0;
   for(unsigned i=0;i<DIM;i++)
    {
     sum += (fe_flux[i]-exact_flux[i])*(fe_flux[i]-exact_flux[i]);
     sum2 += exact_flux[i]*exact_flux[i];
    }
   error += sum*W;
   norm += sum2*W;

  } // End of loop over the integration points
}



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
/// to nodal coordinates (fully analytical). 
/// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
/// Overloads the FD-based version in the FE base class.
//======================================================================
template <unsigned DIM>
void RefineablePoissonEquations<DIM>::get_dresidual_dnodal_coordinates(
 RankThreeTensor<double>&
 dresidual_dnodal_coordinates)
{
 // Determine number of nodes in element
 const unsigned n_node = nnode();

 // Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

 // Get number of shape controlling nodes 
 const unsigned n_shape_controlling_node = nshape_controlling_nodes();

 // Deriatives of shape fct derivatives w.r.t. nodal coords
 RankFourTensor<double> d_dpsidx_dX(DIM,n_shape_controlling_node,n_node,DIM);
 RankFourTensor<double> d_dtestdx_dX(DIM,n_shape_controlling_node,n_node,DIM);

 // Derivative of Jacobian of mapping w.r.t. to nodal coords
 DenseMatrix<double> dJ_dX(DIM,n_shape_controlling_node);

 // Derivatives of derivative of u w.r.t. nodal coords
 RankThreeTensor<double> d_dudx_dX(DIM,n_shape_controlling_node,DIM);

 // Source function and its gradient
 double source;
 Vector<double> d_source_dx(DIM);

 // Index at which the poisson unknown is stored
 const unsigned u_nodal_index = this->u_index_poisson();
 
 // Determine the number of integration points
 const unsigned n_intpt = integral_pt()->nweight();

 // Integer to store the local equation number
 int local_eqn=0;

 // Local storage for pointers to hang_info object
 HangInfo *hang_info_pt=0;

 // Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   // Get the integral weight
   double w = integral_pt()->weight(ipt);

   // Call the derivatives of the shape/test functions, as well as the
   // derivatives of these w.r.t. nodal coordinates and the derivative
   // of the jacobian of the mapping w.r.t. nodal coordinates
   const double J = this->dshape_and_dtest_eulerian_at_knot_poisson(
    ipt,psi,dpsidx,d_dpsidx_dX,test,dtestdx,d_dtestdx_dX,dJ_dX);
   
   // Calculate local values 
   // Allocate and initialise to zero
   Vector<double> interpolated_x(DIM,0.0);
   Vector<double> interpolated_dudx(DIM,0.0);

   // Calculate function value and derivatives:
   // -----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++)
    {
     // Get the nodal value of the Poisson unknown
     double u_value = nodal_value(l,u_nodal_index);

     // Loop over directions
     for(unsigned i=0;i<DIM;i++)
      {
       interpolated_x[i] += nodal_position(l,i)*psi(l);
       interpolated_dudx[i] += u_value*dpsidx(l,i);
      }
    }

   // Calculate derivative of du/dx_i w.r.t. nodal positions X_{pq}

   // Loop over shape controlling nodes
   for(unsigned q=0;q<n_shape_controlling_node;q++)
    {       
     // Loop over coordinate directions
     for(unsigned p=0;p<DIM;p++)
      {
       for(unsigned i=0;i<DIM;i++)
        {
         double aux=0.0;
         for(unsigned j=0;j<n_node;j++)
          {
           aux += nodal_value(j,u_nodal_index)*d_dpsidx_dX(p,q,j,i);
          }
         d_dudx_dX(p,q,i) = aux;
        }
      }
    }

   // Get source function
   this->get_source_poisson(ipt,interpolated_x,source);

   // Get gradient of source function
   this->get_source_gradient_poisson(ipt,interpolated_x,d_source_dx);

//   std::map<Node*,unsigned> local_shape_controlling_node_lookup=shape_controlling_node_lookup();

   // Assemble d res_{local_eqn} / d X_{pq}
   // -------------------------------------
       
   // Loop over the nodes for the test functions 
   for(unsigned l=0;l<n_node;l++)
    {
     // Local variables used to store the number of master nodes and the
     // weight associated with the shape function if the node is hanging
     unsigned n_master=1; 
     double hang_weight=1.0;
     
     // Local bool (is the node hanging)
     bool is_node_hanging = this->node_pt(l)->is_hanging();
     
     // If the node is hanging, get the number of master nodes
     if(is_node_hanging)
      {
       hang_info_pt = this->node_pt(l)->hanging_pt();
       n_master = hang_info_pt->nmaster();
      }
     // Otherwise there is just one master node, the node itself
     else
      {
       n_master = 1;
      }
     
     // Loop over the master nodes
     for(unsigned m=0;m<n_master;m++)
      {
       // Get the local equation number and hang_weight
       // If the node is hanging
       if(is_node_hanging)
        {
         // Read out the local equation number from the m-th master node
         local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                          u_nodal_index);
         
         // Read out the weight from the master node
         hang_weight = hang_info_pt->master_weight(m);
        }
       // If the node is not hanging
       else
        {
         // The local equation number comes from the node itself
         local_eqn = this->nodal_local_eqn(l,u_nodal_index);
         // The hang weight is one
         hang_weight = 1.0;
        }
       
       // If the nodal equation is not a boundary condition
       if(local_eqn >= 0)
        {
         // Loop over coordinate directions
         for(unsigned p=0;p<DIM;p++)
          {              
           // Loop over shape controlling nodes
           for(unsigned q=0;q<n_shape_controlling_node;q++)
            {       
             double sum = source*test(l)*dJ_dX(p,q)
              + d_source_dx[p]*test(l)*psi(q)*J;
             
             for(unsigned i=0;i<DIM;i++)
              {
               sum += interpolated_dudx[i]*(dtestdx(l,i)*dJ_dX(p,q)+
                                               d_dtestdx_dX(p,q,l,i)*J)
                + d_dudx_dX(p,q,i)*dtestdx(l,i)*J;
              }
             
             // Multiply through by integration weight
             dresidual_dnodal_coordinates(local_eqn,p,q) += sum*w*hang_weight;
            }
          }
        }
      }
    }
  } // End of loop over integration points
}

//===============================================================
/// \short Rebuild the element from its sons
/// Adjusts its p-order to be the maximum of its sons' p-orders
//===============================================================
template<>
void PRefineableQPoissonElement<1>::rebuild_from_sons(Mesh* &mesh_pt) 
 {
  // Get p-orders of sons
  unsigned n_sons = this->tree_pt()->nsons();
  Vector<unsigned> son_p_order(n_sons);
  unsigned max_son_p_order = 0;
  for (unsigned ison=0;ison<n_sons;ison++)
   {
    PRefineableElement* el_pt = dynamic_cast<PRefineableElement*>(this->tree_pt()->son_pt(ison)->object_pt());
    son_p_order[ison] = el_pt->p_order();
    if (son_p_order[ison] > max_son_p_order) max_son_p_order = son_p_order[ison];
   }
  
  unsigned old_Nnode = this->nnode();
  unsigned old_P_order = this->p_order();
  // Set p-order of the element
  this->p_order() = max_son_p_order;
  
  // Change integration scheme
  delete this->integral_pt();
  switch(this->p_order())
  {
  case 2:
   this->set_integration_scheme(new GaussLobattoLegendre<1,2>);
   break;
  case 3:
   this->set_integration_scheme(new GaussLobattoLegendre<1,3>);
   break;
  case 4:
   this->set_integration_scheme(new GaussLobattoLegendre<1,4>);
   break;
  case 5:
   this->set_integration_scheme(new GaussLobattoLegendre<1,5>);
   break;
  case 6:
   this->set_integration_scheme(new GaussLobattoLegendre<1,6>);
   break;
  case 7:
   this->set_integration_scheme(new GaussLobattoLegendre<1,7>);
   break;
  default:
   std::ostringstream error_message;
   error_message <<"\nERROR: Exceeded maximum polynomial order for";
   error_message <<"\n       integration scheme.\n";
   throw OomphLibError(error_message.str(),
                       "PRefineableQPoissonElement<1>::rebuild_from_sons()",
                       OOMPH_EXCEPTION_LOCATION);
  }

  // Back up pointers to old nodes before they are lost
  Vector<Node*> old_node_pt(old_Nnode);
  for (unsigned n=0; n<old_Nnode; n++)
   {
    old_node_pt[n] = this->node_pt(n);
   }
  
  // Allocate new space for Nodes (at the element level)
  this->set_n_node(this->p_order());
  
  // Copy vertex nodes and create new edge and internal nodes
  //---------------------------------------------------------
  
  // Copy vertex nodes
  this->node_pt(0) = old_node_pt[0];
  this->node_pt(this->p_order()-1) = old_node_pt[old_P_order-1];
  
  //The timestepper should be the same for all nodes and node 0 should
  //never be deleted.
  if(this->node_pt(0)==0)
   {
    throw OomphLibError("The Corner node (0) does not exist",
                        "PRefineableQPoissonElement<1>::rebuild_from_sons()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  
  TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();
  unsigned ntstorage = time_stepper_pt->ntstorage();

  unsigned jnod=0;
  Vector<double> s_fraction(1), s(1);
  //Loop over the nodes in the element
  unsigned n_p = this->nnode_1d();
  for(unsigned i0=0;i0<n_p;i0++)
   {
    //Get the fractional position of the node
    s_fraction[0] = this->local_one_d_fraction_of_node(i0,0);
    //Local coordinate
    s[0] = -1.0 + 2.0*s_fraction[0];

    //Set the local node number
    jnod = i0;
    
    //If the node has not been built
    if(this->node_pt(jnod)==0)
     {
      //First, find the son element in which the node should live
      
      //Find coordinates in the sons
      Vector<double> s_in_son(2);
      using namespace BinaryTreeNames;
      int son=-10;
      //If negative on the west side
      if(s_fraction[0] < 0.5)
       {
        //It's the southwest son
        son = L;
        s_in_son[0] =  -1.0 + 4.0*s_fraction[0];
       }
      else
       {
        //On the south side
        son = R;
        s_in_son[0] = -1.0 + 4.0*(s_fraction[0]-0.5);
       }
      
      //Get the pointer to the son element in which the new node
      //would live
      PRefineableQElement<1>* son_el_pt = 
       dynamic_cast<PRefineableQElement<1>*>(
        this->tree_pt()->son_pt(son)->object_pt());
      
      //Construct the new node
      this->node_pt(jnod) = this->construct_node(jnod,time_stepper_pt);
      
      //Now we set the position and values at the newly created node
      
      // In the first instance use macro element or FE representation
      // to create past and present nodal positions.
      // (THIS STEP SHOULD NOT BE SKIPPED FOR ALGEBRAIC
      // ELEMENTS AS NOT ALL OF THEM NECESSARILY IMPLEMENT
      // NONTRIVIAL NODE UPDATE FUNCTIONS. CALLING
      // THE NODE UPDATE FOR SUCH ELEMENTS/NODES WILL LEAVE
      // THEIR NODAL POSITIONS WHERE THEY WERE (THIS IS APPROPRIATE
      // ONCE THEY HAVE BEEN GIVEN POSITIONS) BUT WILL
      // NOT ASSIGN SENSIBLE INITIAL POSITONS!
      
      // Loop over # of history values
      for(unsigned t=0;t<ntstorage;t++)
       {
        using namespace BinaryTreeNames;
        //Get the position from the son
        Vector<double> x_prev(1);
        
        //Now let's fill in the value
        son_el_pt->get_x(t,s_in_son,x_prev);
        for(unsigned i=0;i<1;i++)
         {
          this->node_pt(jnod)->x(t,i) = x_prev[i];
         }
       }
      
      // Now set up the values
      // Loop over all history values
      for(unsigned t=0;t<ntstorage;t++)
       {
        // Get values from father element
        // Note: get_interpolated_values() sets Vector size itself.
        Vector<double> prev_values;
        son_el_pt->get_interpolated_values(t,s_in_son,prev_values);
        
        //Initialise the values at the new node
        for(unsigned k=0;k<this->node_pt(jnod)->nvalue();k++)
         {
          this->node_pt(jnod)->set_value(t,k,prev_values[k]);
         }
       }
      
      //Add the node to the mesh
      //(Must be added because if it wasn't needed, it wouldn't have
      //been created)
      mesh_pt->add_node_pt(this->node_pt(jnod));
      
     } //End of the case when we build the node ourselves
    
    //Algebraic stuff here
    //Check whether the element is an algebraic element
    AlgebraicElementBase* alg_el_pt =
     dynamic_cast<AlgebraicElementBase*>(this);
    
    //If we do have an algebraic element
    if(alg_el_pt!=0)
     {
      std::string error_message =
       "Have not implemented rebuilding from sons for";
      error_message +=
       "Algebraic Spectral elements yet\n";
      
      throw 
       OomphLibError(error_message,
                     "PRefineableQPoissonElement<1>::rebuild_from_sons()",
                     OOMPH_EXCEPTION_LOCATION);
     }
    
   }
 }

//===============================================================
/// \short Rebuild the element from its sons
/// Adjusts its p-order to be the maximum of its sons' p-orders
//===============================================================
template<>
void PRefineableQPoissonElement<2>::rebuild_from_sons(Mesh* &mesh_pt) 
 {
  // Get p-orders of sons
  unsigned n_sons = this->tree_pt()->nsons();
  Vector<unsigned> son_p_order(n_sons);
  unsigned max_son_p_order = 0;
  for (unsigned ison=0;ison<n_sons;ison++)
   {
    PRefineableElement* el_pt = dynamic_cast<PRefineableElement*>(this->tree_pt()->son_pt(ison)->object_pt());
    son_p_order[ison] = el_pt->p_order();
    if (son_p_order[ison] > max_son_p_order) max_son_p_order = son_p_order[ison];
   }
  
  unsigned old_Nnode = this->nnode();
  unsigned old_P_order = this->p_order();
  // Set p-order of the element
  this->p_order() = max_son_p_order;
  
  // Change integration scheme
  delete this->integral_pt();
  switch(this->p_order())
  {
  case 2:
   this->set_integration_scheme(new GaussLobattoLegendre<2,2>);
   break;
  case 3:
   this->set_integration_scheme(new GaussLobattoLegendre<2,3>);
   break;
  case 4:
   this->set_integration_scheme(new GaussLobattoLegendre<2,4>);
   break;
  case 5:
   this->set_integration_scheme(new GaussLobattoLegendre<2,5>);
   break;
  case 6:
   this->set_integration_scheme(new GaussLobattoLegendre<2,6>);
   break;
  case 7:
   this->set_integration_scheme(new GaussLobattoLegendre<2,7>);
   break;
  default:
   std::ostringstream error_message;
   error_message <<"\nERROR: Exceeded maximum polynomial order for";
   error_message <<"\n       integration scheme.\n";
   throw OomphLibError(error_message.str(),
                       "PRefineableQPoissonElement<2>::rebuild_from_sons()",
                       OOMPH_EXCEPTION_LOCATION);
  }

  // Back up pointers to old nodes before they are lost
  Vector<Node*> old_node_pt(old_Nnode);
  for (unsigned n=0; n<old_Nnode; n++)
   {
    old_node_pt[n] = this->node_pt(n);
   }
  
  // Allocate new space for Nodes (at the element level)
  this->set_n_node(this->p_order()*this->p_order());
  
  // Copy vertex nodes and create new edge and internal nodes
  //---------------------------------------------------------
  
  // Copy vertex nodes
  this->node_pt(0) = old_node_pt[0];
  this->node_pt(this->p_order()-1) = old_node_pt[old_P_order-1];
  this->node_pt(this->p_order()*(this->p_order()-1))
   = old_node_pt[(old_P_order)*(old_P_order-1)];
  this->node_pt(this->p_order()*this->p_order()-1)
   = old_node_pt[(old_P_order)*(old_P_order)-1];
  
  //The timestepper should be the same for all nodes and node 0 should
  //never be deleted.
  if(this->node_pt(0)==0)
   {
    throw OomphLibError("The Corner node (0) does not exist",
                        "PRefineableQPoissonElement<2>::rebuild_from_sons()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  
  TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();
  unsigned ntstorage = time_stepper_pt->ntstorage();

  unsigned jnod=0;
  Vector<double> s_fraction(2), s(2);
  //Loop over the nodes in the element
  unsigned n_p = this->nnode_1d();
  for(unsigned i0=0;i0<n_p;i0++)
   {
    //Get the fractional position of the node
    s_fraction[0] = this->local_one_d_fraction_of_node(i0,0);
    //Local coordinate
    s[0] = -1.0 + 2.0*s_fraction[0];

    for(unsigned i1=0;i1<n_p;i1++)
     {
      //Get the fractional position of the node in the direction of s[1]
      s_fraction[1] = this->local_one_d_fraction_of_node(i1,1);
      // Local coordinate in father element
      s[1] = -1.0 + 2.0*s_fraction[1];

      //Set the local node number
      jnod = i0 + n_p*i1;

      //If the node has not been built
      if(this->node_pt(jnod)==0)
       {
        //Has the node been created by one of its neighbours
        bool is_periodic = false;
        Node* created_node_pt = 
         this->node_created_by_neighbour(s_fraction,is_periodic);
        
        //If it has set the pointer
        if(created_node_pt!=0)
         {
          //If the node is periodic
          if(is_periodic)
           {
            throw OomphLibError(
             "Cannot handle periodic nodes in refineable spectral elements yet",
             "PRefineableQPoissonElement<2>::rebuild_from_sons()",
             OOMPH_EXCEPTION_LOCATION);
           }
          //Non-periodic case, just set the pointer
          else
           {
            this->node_pt(jnod) = created_node_pt;
           }
         }
        //Otherwise, we need to build it
        else
         {
          //First, find the son element in which the node should live
          
          //Find coordinates in the sons
          Vector<double> s_in_son(2);
          using namespace QuadTreeNames;
          int son=-10;
          //If negative on the west side
          if(s_fraction[0] < 0.5)
           {
            //On the south side
            if(s_fraction[1] < 0.5)
             {
              //It's the southwest son
              son = SW;
              s_in_son[0] =  -1.0 + 4.0*s_fraction[0];
              s_in_son[1] =  -1.0 + 4.0*s_fraction[1];
             }
            //On the north side
            else
             {
              //It's the northwest son
              son = NW;
              s_in_son[0] = -1.0 + 4.0*s_fraction[0];
              s_in_son[1] = -1.0 + 4.0*(s_fraction[1]-0.5);
             }
           }
          else
           {
            //On the south side
            if(s_fraction[1] < 0.5)
             {
              //It's the southeast son
              son = SE;
              s_in_son[0] =  -1.0 + 4.0*(s_fraction[0]-0.5);
              s_in_son[1] =  -1.0 + 4.0*s_fraction[1];
             }
            //On the north side
            else
             {
              //It's the northeast son
              son = NE;
              s_in_son[0] = -1.0 + 4.0*(s_fraction[0]-0.5);
              s_in_son[1] = -1.0 + 4.0*(s_fraction[1]-0.5);
             }
           }

          //Get the pointer to the son element in which the new node
          //would live
          PRefineableQElement<2>* son_el_pt = 
           dynamic_cast<PRefineableQElement<2>*>(
            this->tree_pt()->son_pt(son)->object_pt());
          
          //If we are rebuilding, then worry about the boundary conditions
          //Find the boundary of the node
          //Initially none
          int boundary=Tree::OMEGA;
          //If we are on the western boundary
          if(i0==0) {boundary = W;}
          //If we are on the eastern boundary
          else if(i0==n_p-1) {boundary = E;}
          
          //If we are on the southern boundary
          if(i1==0)
           {
            //If we already have already set the boundary, we're on a corner
            switch(boundary)
             {
             case W:
              boundary = SW;
              break;
             case E:
              boundary = SE;
              break;
              //Boundary not set
             default:
              boundary = S;
              break;
             }
           }
          //If we are the northern bounadry
          else if(i1==n_p-1)
           {
            //If we already have a boundary
            switch(boundary)
             {
             case W:
              boundary = NW;
              break;
             case E:
              boundary = NE;
              break;
             default:
              boundary = N;
              break;
             }
           }
          
          // set of boundaries that this edge in the son lives on
          std::set<unsigned> boundaries;

          //Now get the boundary conditions from the son
          //The boundaries will be common to the son because there can be
          //no rotations here
          if(boundary!=Tree::OMEGA)
           {
            son_el_pt->get_boundaries(boundary,boundaries);
           }
          
          // If the node lives on a boundary: 
          // Construct a boundary node, 
          // Get boundary conditions and
          // update all lookup schemes
          if(boundaries.size()>0)
           {
            //Construct the new node
            this->node_pt(jnod) = 
             this->construct_boundary_node(jnod,time_stepper_pt);
            
            //Get the boundary conditions from the son
            Vector<int> bound_cons(this->ncont_interpolated_values());
            son_el_pt->get_bcs(boundary,bound_cons);
            
            //Loop over the values and pin if necessary
            unsigned nval = this->node_pt(jnod)->nvalue();
            for(unsigned k=0;k<nval;k++)
             {
              if(bound_cons[k]) {this->node_pt(jnod)->pin(k);}
             }
            
            // Solid node? If so, deal with the positional boundary
            // conditions:
            SolidNode* solid_node_pt =
             dynamic_cast<SolidNode*>(this->node_pt(jnod));
            if (solid_node_pt!=0)
             {
              //Get the positional boundary conditions from the father:
              unsigned n_dim = this->node_pt(jnod)->ndim();
              Vector<int> solid_bound_cons(n_dim);
              RefineableSolidQElement<2>* son_solid_el_pt=
               dynamic_cast<RefineableSolidQElement<2>*>(son_el_pt);
#ifdef PARANOID
              if (son_solid_el_pt==0)
               {
                std::string error_message =
                 "We have a SolidNode outside a refineable SolidElement\n";
                error_message +=
                 "during mesh refinement -- this doesn't make sense\n";

                throw OomphLibError(
                 error_message,
                 "PRefineableQPoissonElement<2>::rebuild_from_sons()",
                 OOMPH_EXCEPTION_LOCATION);
               }
#endif
              son_solid_el_pt->get_solid_bcs(boundary,solid_bound_cons);
              
              //Loop over the positions and pin, if necessary
              for(unsigned k=0;k<n_dim;k++)
               {
                if (solid_bound_cons[k]) {solid_node_pt->pin_position(k);}
               }
             }

            //Next we update the boundary look-up schemes
            //Loop over the boundaries stored in the set
            for(std::set<unsigned>::iterator it = boundaries.begin();
                it != boundaries.end(); ++it)
             {
              //Add the node to the boundary
              mesh_pt->add_boundary_node(*it,this->node_pt(jnod));
              
              //If we have set an intrinsic coordinate on this
              //mesh boundary then it must also be interpolated on
              //the new node
              //Now interpolate the intrinsic boundary coordinate
              if(mesh_pt->boundary_coordinate_exists(*it)==true)
               {
                Vector<double> zeta(1);
                son_el_pt->interpolated_zeta_on_edge(*it,boundary,
                                                     s_in_son,zeta);
                
                this->node_pt(jnod)
                 ->set_coordinates_on_boundary(*it,zeta);
               }
             }
           }
          //Otherwise the node is not on a Mesh boundary 
          //and we create a normal "bulk" node
          else
           {
            //Construct the new node
            this->node_pt(jnod) = this->construct_node(jnod,time_stepper_pt);
           }

          //Now we set the position and values at the newly created node
          
          // In the first instance use macro element or FE representation
          // to create past and present nodal positions.
          // (THIS STEP SHOULD NOT BE SKIPPED FOR ALGEBRAIC
          // ELEMENTS AS NOT ALL OF THEM NECESSARILY IMPLEMENT
          // NONTRIVIAL NODE UPDATE FUNCTIONS. CALLING
          // THE NODE UPDATE FOR SUCH ELEMENTS/NODES WILL LEAVE
          // THEIR NODAL POSITIONS WHERE THEY WERE (THIS IS APPROPRIATE
          // ONCE THEY HAVE BEEN GIVEN POSITIONS) BUT WILL
          // NOT ASSIGN SENSIBLE INITIAL POSITONS!
            
          // Loop over # of history values
          for(unsigned t=0;t<ntstorage;t++)
           {
            using namespace QuadTreeNames;
            //Get the position from the son
            Vector<double> x_prev(2);
            
            //Now let's fill in the value
            son_el_pt->get_x(t,s_in_son,x_prev);
            for(unsigned i=0;i<2;i++)
             {
              this->node_pt(jnod)->x(t,i) = x_prev[i];
             }
           }

          // Now set up the values
          // Loop over all history values
          for(unsigned t=0;t<ntstorage;t++)
           {
            // Get values from father element
            // Note: get_interpolated_values() sets Vector size itself.
            Vector<double> prev_values;
            son_el_pt->get_interpolated_values(t,s_in_son,prev_values);
            
            //Initialise the values at the new node
            for(unsigned k=0;k<this->node_pt(jnod)->nvalue();k++)
             {
              this->node_pt(jnod)->set_value(t,k,prev_values[k]);
             }
           }
          
          //Add the node to the mesh
          //(Must be added because if it wasn't needed, it wouldn't have
          //been created)
          mesh_pt->add_node_pt(this->node_pt(jnod));

         } //End of the case when we build the node ourselves

        //Algebraic stuff here
        //Check whether the element is an algebraic element
        AlgebraicElementBase* alg_el_pt =
         dynamic_cast<AlgebraicElementBase*>(this);
        
        //If we do have an algebraic element
        if(alg_el_pt!=0)
         {
          std::string error_message =
           "Have not implemented rebuilding from sons for";
          error_message +=
             "Algebraic Spectral elements yet\n";

            throw 
             OomphLibError(error_message,
                           "PRefineableQPoissonElement<2>::rebuild_from_sons()",
                           OOMPH_EXCEPTION_LOCATION);
           }
          
       }
     }
   }
 }

//===============================================================
/// \short Rebuild the element from its sons
/// Adjusts its p-order to be the maximum of its sons' p-orders
//===============================================================
template<>
void PRefineableQPoissonElement<3>::rebuild_from_sons(Mesh* &mesh_pt) 
 {
  throw
   OomphLibError("This function is not yet implemented.",
                 "PRefineableQPoissonElement<3>::rebuild_from_sons()",
                 OOMPH_EXCEPTION_LOCATION);
 }

/// Get error against and norm of exact solution
template<unsigned DIM>
void PRefineableQPoissonElement<DIM>::compute_energy_error(
 std::ostream &outfile, FiniteElement::SteadyExactSolutionFctPt exact_grad_pt,
 double& error, double& norm)
 {
  // Initialise
  error=0.0;
  norm=0.0;

  //Vector of local coordinates
  Vector<double> s(DIM);

  // Vector for coordintes
  Vector<double> x(DIM);
  
  //Set the value of n_intpt
  unsigned n_intpt = this->integral_pt()->nweight();
  
  // Setup output structure: Conversion is fishy but it's only output...
  unsigned nplot;
  if (DIM==1)
   {
    nplot=n_intpt;
   }
  else 
   {
    nplot=unsigned(pow(n_intpt,1.0/double(DIM)));
   }
  
  // Tecplot header info
  outfile << this->tecplot_zone_string(nplot);
  
  // Exact gradient Vector
  Vector<double> exact_grad(DIM);
  
  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    
    //Assign values of s
    for(unsigned i=0;i<DIM;i++)
     {
      s[i] = this->integral_pt()->knot(ipt,i);
     }
    
    //Get the integral weight
    double w = this->integral_pt()->weight(ipt);
    
    // Get jacobian of mapping
    double J=this->J_eulerian(s);
    
    //Premultiply the weights and the Jacobian
    double W = w*J;
    
    // Get x position as Vector
    this->interpolated_x(s,x);
    
    // Get FE du/dx
    Vector<double> dudx_fe(DIM);
    PoissonEquations<DIM>::get_flux(s,dudx_fe);
    
    // Get exact gradient at this point
    (*exact_grad_pt)(x,exact_grad);
    
    //Output x,y,...,error
    for(unsigned i=0;i<DIM;i++)
     {
      outfile << x[i] << " ";
     }
    for(unsigned i=0;i<DIM;i++)
     {
      outfile << exact_grad[i] << " " << exact_grad[i]-dudx_fe[i] << std::endl;
     }
    
    // Add to error and norm
    for(unsigned i=0;i<DIM;i++)
     {
      norm+=exact_grad[i]*exact_grad[i]*W;
      error+=(exact_grad[i]-dudx_fe[i])*(exact_grad[i]-dudx_fe[i])*W;
     }
    
   }
 }

template<unsigned DIM>
void PRefineableQPoissonElement<DIM>::further_build()
 {
  if(this->tree_pt()->father_pt()!=0)
   {
    // Needed to set the source function pointer (if there is a father)
    RefineablePoissonEquations<DIM>::further_build();
   }
  // Now do the PRefineableQElement further_build()
  PRefineableQElement<DIM>::further_build();
 }



//====================================================================
// Force build of templates
//====================================================================
template class RefineableQPoissonElement<1,2>;
template class RefineableQPoissonElement<1,3>;
template class RefineableQPoissonElement<1,4>;

template class RefineableQPoissonElement<2,2>;
template class RefineableQPoissonElement<2,3>;
template class RefineableQPoissonElement<2,4>;

template class RefineableQPoissonElement<3,2>;
template class RefineableQPoissonElement<3,3>;
template class RefineableQPoissonElement<3,4>;

template class PRefineableQPoissonElement<1>;
template class PRefineableQPoissonElement<2>;
template class PRefineableQPoissonElement<3>;

}
