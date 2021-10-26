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
//Header for a multi-physics problem that couples a Navier--Stokes
//mesh to an advection diffusion mesh, giving Boussinesq convection

//Oomph-lib headers, we require the generic, advection-diffusion
//and navier-stokes elements.
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;

//======================class definitions==============================
/// Build RefineableQCrouzeixRaviartElementWithExternalElement 
/// that inherits from 
/// ElementWithExternalElement so that it can "communicate" with 
/// RefineableQAdvectionDiffusionElementWithExternalElement
//=====================================================================
template<unsigned DIM>
class RefineableQCrouzeixRaviartElementWithTwoExternalElement : 
 public virtual RefineableQCrouzeixRaviartElement<DIM>,
 public virtual ElementWithExternalElement
{

private:

 /// Pointer to a private data member, the thermal Rayleigh number
 double* Ra_T_pt;

 /// Pointer to the private data member, the solutal Rayleigh number
 double* Ra_S_pt;

 /// The static default value of the Rayleigh number
 static double Default_Physical_Constant_Value;

public: 

 /// Constructor: call the underlying constructors and 
 /// initialise the pointer to the Rayleigh number to point
 /// to the default value of 0.0.
 RefineableQCrouzeixRaviartElementWithTwoExternalElement() : 
  RefineableQCrouzeixRaviartElement<DIM>(),
  ElementWithExternalElement()
  {
   Ra_T_pt = &Default_Physical_Constant_Value;

   Ra_S_pt = &Default_Physical_Constant_Value;

   //There are two interactions
   this->set_ninteraction(2);

   //We do not need to add any external geometric data because the
   //element is fixed
   this->ignore_external_geometric_data();
  } 

 /// Access function for the Rayleigh number (const version)
 const double &ra_t() const {return *Ra_T_pt;}

 /// Access function for the pointer to the Rayleigh number
 double* &ra_t_pt() {return Ra_T_pt;}


  /// Access function for the solutal Rayleigh number (const version)
 const double &ra_s() const {return *Ra_S_pt;}

 /// Access function for the pointer to the solutal Rayleigh number
 double* &ra_s_pt() {return Ra_S_pt;}

 /// Call the underlying single-physics element's further_build()
 /// functions and make sure that the pointer to the Rayleigh number
 /// is passed to the sons
 void further_build()
  {
   RefineableQCrouzeixRaviartElement<DIM>::further_build();

   //Cast the pointer to the father element to the specific
   //element type
   RefineableQCrouzeixRaviartElementWithTwoExternalElement<DIM>* 
    cast_father_element_pt
    = dynamic_cast<
    RefineableQCrouzeixRaviartElementWithTwoExternalElement<DIM>*>(
     this->father_element_pt());

   //Set the pointer to the Rayleigh numbers to be the same as that in
   //the father
   this->Ra_T_pt = cast_father_element_pt->ra_t_pt();

   this->Ra_S_pt = cast_father_element_pt->ra_s_pt();
  }

 // Overload get_body_force_nst to get the temperature "body force"
 // from the "source" AdvectionDiffusion element via current integration point
 void get_body_force_nst(const double& time, const unsigned& ipt, 
                         const Vector<double> &s, const Vector<double> &x, 
                         Vector<double> &result);

/// Fill in the derivatives of the body force with respect to the
 /// external unknowns
 void get_dbody_force_nst_dexternal_element_data(
  const unsigned& ipt, 
  DenseMatrix<double> &result, Vector<unsigned> &global_eqn_number);


 /// Fill in the constituent elements' contribution to the residual vector.
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the residuals of the Navier-Stokes equations
   RefineableNavierStokesEquations<DIM>::fill_in_contribution_to_residuals(
    residuals);
  }

 /// Compute the element's residual vector and the Jacobian matrix.
 /// Jacobian is computed by finite-differencing.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
#ifdef USE_FD_JACOBIAN_FOR_MY_NAVIER_STOKES_ELEMENT   
   // This function computes the Jacobian by finite-differencing
   FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
#else
   //Get the contribution from the basic Navier--Stokes element
   RefineableQCrouzeixRaviartElement<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);
   //Get the off-diagonal terms analytically
   this->fill_in_off_diagonal_block_analytic(residuals,jacobian);
#endif

  }

 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   //Call the standard (Broken) function
   //which will prevent these elements from being used
   //in eigenproblems until replaced.
   FiniteElement::fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);
  }

 /// Compute the contribution of the external
 /// degrees of freedom (temperatures) on the Navier-Stokes equations
 void fill_in_off_diagonal_block_analytic(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
  {
   //Local storage for the index in the nodes at which the
   //Navier-Stokes velocities are stored (we know that this should be 0,1,2)
   unsigned u_nodal_nst[DIM];
   for(unsigned i=0;i<DIM;i++) 
    {u_nodal_nst[i] = this->u_index_nst(i);}

   //Find out how many nodes there are
   const unsigned n_node = this->nnode();
   
   //Set up memory for the shape and test functions and their derivatives
   Shape psif(n_node), testf(n_node);
   DShape dpsifdx(n_node,DIM), dtestfdx(n_node,DIM);
   
   //Number of integration points
   const unsigned n_intpt = this->integral_pt()->nweight();
   
   //Integers to store the local equations and unknowns
   int local_eqn=0, local_unknown=0;
   
   // Local storage for pointers to hang_info objects
   HangInfo *hang_info_pt=0;   

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = this->integral_pt()->weight(ipt);
     
     //Call the derivatives of the shape and test functions
     double J = 
      this->dshape_and_dtest_eulerian_at_knot_nst(ipt,psif,dpsifdx,
                                                  testf,dtestfdx);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;
     
     //Assemble the jacobian terms
     
     //Get the derivatives of the body force wrt the unknowns
     //of the external element
     DenseMatrix<double> dbody_dexternal_element_data;
     //Vector of global equation number corresponding to the external
     //element's data
     Vector<unsigned> global_eqn_number_of_external_element_data;
     //Get the appropriate derivatives
     this->get_dbody_force_nst_dexternal_element_data(
      ipt,dbody_dexternal_element_data,
      global_eqn_number_of_external_element_data);
     //Find out how many external data there are
     const unsigned n_external_element_data = 
      global_eqn_number_of_external_element_data.size();

     //Loop over the test functions
     for(unsigned l=0;l<n_node;l++)
      {
       //Assemble the contributions of the temperature to 
       //the Navier--Stokes equations (which arise through the buoyancy
       //body-force term)
       unsigned n_master = 1;
       double hang_weight = 1.0;
       
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
         //If the node is hanging get weight from master node
         if(is_node_hanging)
          {
           //Get the hang weight from the master node
           hang_weight = hang_info_pt->master_weight(m);
          }
         else
          {
           // Node contributes with full weight
           hang_weight = 1.0;
          }
         
         
         //Loop over the velocity components in the Navier--Stokes equtions
         for(unsigned i=0;i<DIM;i++)
          {
           //Get the equation number
           if(is_node_hanging)
            {
             //Get the equation number from the master node
             local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                              u_nodal_nst[i]);
            }
           else
            {
             // Local equation number
             local_eqn = this->nodal_local_eqn(l,u_nodal_nst[i]);
            }
           
           if(local_eqn >= 0)
            {
             //Loop over the external data
             for(unsigned l2=0;l2<n_external_element_data;l2++)
              { 
               //Find the local equation number corresponding to the global
               //unknown
               local_unknown = 
                this->local_eqn_number(
                 global_eqn_number_of_external_element_data[l2]);
               if(local_unknown >= 0)
                {
                 //Add contribution to jacobian matrix
                 jacobian(local_eqn,local_unknown) 
                  += dbody_dexternal_element_data(i,l2)*testf(l)*hang_weight*W;
                }
              }
            }
          }
        }
      }
    }
  } 

};

//======================class definitions==============================
/// Build MyAdvectionDiffusionElement that inherits from 
/// ElementWithExternalElement so that it can "communicate" with the 
/// MyNavierStokesElement
//=====================================================================
template<unsigned DIM>
class RefineableQAdvectionDiffusionElementWithExternalElement : 
 public virtual RefineableQAdvectionDiffusionElement<DIM,3>,
 public virtual ElementWithExternalElement
{

public:

 /// Constructor: call the underlying constructors
 RefineableQAdvectionDiffusionElementWithExternalElement() : 
  RefineableQAdvectionDiffusionElement<DIM,3>(), ElementWithExternalElement()
  { 
   //There is one interaction
   this->set_ninteraction(1);

   //We do not need to add any external geometric data because the
   //element is fixed
   this->ignore_external_geometric_data();
  }

 /// Overload the wind function in the advection-diffusion equations.
 /// This provides the coupling from the Navier--Stokes equations to the
 /// advection-diffusion equations because the wind is the fluid velocity,
 /// obtained from the source element in the other mesh
 void get_wind_adv_diff(const unsigned& ipt, const Vector<double> &s, 
                        const Vector<double>& x, Vector<double>& wind) const;

 /// Fill in the derivatives of the wind with respect to the
 /// external unknowns
 void get_dwind_adv_diff_dexternal_element_data(
  const unsigned& ipt, const unsigned &i,
  Vector<double> &result, Vector<unsigned> &global_eqn_number);


 /// Just call the fill_in_residuals for AdvDiff
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   RefineableAdvectionDiffusionEquations<DIM>::
    fill_in_contribution_to_residuals(residuals);
  }

 /// Compute the element's residual vector and the Jacobian matrix.
 /// Jacobian is computed by finite-differencing.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
#ifdef USE_FD_JACOBIAN_FOR_MY_ADVECTION_DIFFUSION_ELEMENT   
   // This function computes the Jacobian by finite-differencing
   FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
#else
   //Get the contribution from the basic Navier--Stokes element
   RefineableQAdvectionDiffusionElement<DIM,3>::
    fill_in_contribution_to_jacobian(residuals,jacobian);
   //Get the off-diagonal terms analytically
   this->fill_in_off_diagonal_block_analytic(residuals,jacobian);
#endif
  }

 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   //Call the standard (Broken) function
   //which will prevent these elements from being used
   //in eigenproblems until replaced.
   FiniteElement::fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);
  }


 /// Compute the contribution of the external
 /// degrees of freedom (velocities) on the advection-diffusion equations
 void fill_in_off_diagonal_block_analytic(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
  {
   //Local storage for the index in the nodes at which the temperature 
   //is stored
   const unsigned u_nodal_adv_diff = this->u_index_adv_diff();

   //Find out how many nodes there are
   const unsigned n_node = this->nnode();
   
   //Set up memory for the shape and test functions and their derivatives
   Shape psi(n_node), test(n_node);
   DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);
   
   //Number of integration points
   const unsigned n_intpt = this->integral_pt()->nweight();
   
   //Integers to store the local equations and unknowns
   int local_eqn=0, local_unknown=0;
   
   // Local storage for pointers to hang_info objects
   HangInfo *hang_info_pt=0;   

   //Get the peclet number
   const double peclet = this->pe();

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = this->integral_pt()->weight(ipt);
     
     //Call the derivatives of the shape and test functions
     double J = 
      this->dshape_and_dtest_eulerian_at_knot_adv_diff(ipt,psi,dpsidx,
                                                       test,dtestdx);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;
     
     //Calculate local values of the derivatives of the solution
     Vector<double> interpolated_dudx(DIM,0.0);
     // Loop over nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       // Loop over directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_dudx[j] += 
          this->nodal_value(l,u_nodal_adv_diff)*dpsidx(l,j);
        }
      }
     
     //Get the derivatives of the wind wrt the unknowns
     //of the external element
     Vector<double> dwind_dexternal_element_data;
     //Vector of global equation number corresponding to the external
     //element's data
     Vector<unsigned> global_eqn_number_of_external_element_data;

     //Loop over the wind directions
     for(unsigned i2=0;i2<DIM;i2++)
      {
       //Get the appropriate derivatives
       this->get_dwind_adv_diff_dexternal_element_data(
        ipt,i2,dwind_dexternal_element_data,
        global_eqn_number_of_external_element_data);
       

       //Find out how many external data there are
       const unsigned n_external_element_data = 
        global_eqn_number_of_external_element_data.size();
       
       //Loop over the test functions
       for(unsigned l=0;l<n_node;l++)
        {
         //Assemble the contributions of the velocities to 
         //the advection-diffusion equations
         unsigned n_master = 1;
         double hang_weight = 1.0;
         
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
           //If the node is hanging get weight from master node
           if(is_node_hanging)
            {
             //Get the hang weight from the master node
             hang_weight = hang_info_pt->master_weight(m);
            }
           else
            {
             // Node contributes with full weight
             hang_weight = 1.0;
            }
           
           //Get the equation number
           if(is_node_hanging)
            {
             //Get the equation number from the master node
             local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                              u_nodal_adv_diff);
            }
           else
            {
             // Local equation number
             local_eqn = this->nodal_local_eqn(l,u_nodal_adv_diff);
            }
           
           if(local_eqn >= 0)
            {
             //Loop over the external data
             for(unsigned l2=0;l2<n_external_element_data;l2++)
              { 
               //Find the local equation number corresponding to the global
               //unknown
               local_unknown = 
                this->local_eqn_number(
                 global_eqn_number_of_external_element_data[l2]);
               if(local_unknown >= 0)
                {
                 //Add contribution to jacobian matrix
                 jacobian(local_eqn,local_unknown) 
                  -= peclet*dwind_dexternal_element_data[l2]*
                  interpolated_dudx[i2]*test(l)*hang_weight*W;
                }
              }
            }
          }
        }
      }
    }
  }

};

//======================start_of_get_body_force_nst============================
// Overload get_body_force_nst to get the temperature "body force"
// from the "source" AdvectionDiffusion element via current integration point
//=============================================================================
template<unsigned DIM>
void RefineableQCrouzeixRaviartElementWithTwoExternalElement<DIM>::
get_body_force_nst
(const double& time,const unsigned& ipt,const Vector<double> &s,
 const Vector<double> &x,Vector<double> &result)
{
 //Let's choose the thermal interaction to be interaction zeo
 const unsigned t_interaction=0;
 //The solutal interaction is interaction 1
 const unsigned c_interaction=1;

 //Get the temperature field from the external element
 const double interpolated_t = 
  dynamic_cast<AdvectionDiffusionEquations<DIM>*>
  (external_element_pt(t_interaction,ipt))->
  interpolated_u_adv_diff(external_element_local_coord(t_interaction,ipt));

 
 //Get the interpolated concentration field from the external element
 const double interpolated_c =  dynamic_cast<AdvectionDiffusionEquations<DIM>*>
  (external_element_pt(c_interaction,ipt))->
  interpolated_u_adv_diff(external_element_local_coord(c_interaction,ipt));

 //Combine into the buoyancy force
 const double buoyancy = interpolated_t*ra_t() + interpolated_c*ra_s();


 // Get vector that indicates the direction of gravity from
 // the Navier-Stokes equations
 Vector<double> gravity(NavierStokesEquations<DIM>::g());
   
 // Temperature-dependent body force:
 for (unsigned i=0;i<DIM;i++)
  {
   result[i] = -gravity[i]*buoyancy;
  }

} //end of get_body_force_nst

//==========start_of_get_dbody_force=========== ===========================
/// Fill in the derivatives of the body force with respect to the external
/// unknowns in the Navier--Stokes equations
//=========================================================================
template<unsigned DIM>
void RefineableQCrouzeixRaviartElementWithTwoExternalElement<DIM>::
get_dbody_force_nst_dexternal_element_data(const unsigned &ipt,
                                           DenseMatrix<double> &result,
                                           Vector<unsigned> &global_eqn_number)
{
 // The temperature interaction index is 0 in this case
 unsigned t_interaction=0;
 // The concentration interaction index is 1 in this case
 unsigned c_interaction = 1;

 //Get the temperature interactions
 Vector<double> du_temp_ddata;
 Vector<unsigned> global_eqn_number_temp;

 //Get the interation data from the temperature element
 dynamic_cast<AdvectionDiffusionEquations<DIM>*>
  (external_element_pt(t_interaction,ipt))->
  dinterpolated_u_adv_diff_ddata(
  external_element_local_coord(t_interaction,ipt),du_temp_ddata,
  global_eqn_number_temp);

 //Get the concentration interactions
 Vector<double> du_conc_ddata;
 Vector<unsigned> global_eqn_number_conc;

 //Get the interation data from the temperature element
 dynamic_cast<AdvectionDiffusionEquations<DIM>*>
  (external_element_pt(c_interaction,ipt))->
  dinterpolated_u_adv_diff_ddata(
  external_element_local_coord(c_interaction,ipt),du_conc_ddata,
  global_eqn_number_conc);
 
 // Get vector that indicates the direction of gravity from
 // the Navier-Stokes equations
 Vector<double> gravity(NavierStokesEquations<DIM>::g());
  
 //Find the number of external data
 //Assuming that the temperature and concentration elements are separate
 //which they are!
 unsigned n_external_temp_data = du_temp_ddata.size();
 unsigned n_external_conc_data = du_conc_ddata.size();

 //Set the size of the matrix to be returned
 result.resize(DIM,n_external_temp_data+n_external_conc_data);

 // Temperature-dependent body force:
 for (unsigned i=0;i<DIM;i++)
  {
   //Loop over the temperature external data
   for(unsigned n=0;n<n_external_temp_data;n++)
    {
     result(i,n) = -gravity[i]*du_temp_ddata[n]*ra_t();
    }

   //Loop over the concentration external data
   //Loop over the temperature external data
   for(unsigned n=0;n<n_external_conc_data;n++)
    {
     result(i,n_external_temp_data+n) = -gravity[i]*du_conc_ddata[n]*ra_s();
    }
  }

 //Set the size of the global equation numbers
 global_eqn_number.resize(n_external_temp_data+n_external_conc_data);
 //Fill in the entries temperature first
 for(unsigned n=0;n<n_external_temp_data;n++)
  {
   global_eqn_number[n] = global_eqn_number_temp[n];
  }
 //Concentration second
 for(unsigned n=0;n<n_external_conc_data;n++)
  {
   global_eqn_number[n_external_temp_data+n] = global_eqn_number_conc[n];
  }
 

} // end_of_get_dbody_force



//==========================start_of_get_wind_adv_diff====================
/// Overload the wind function in the advection-diffusion equations.
/// This provides the coupling from the Navier--Stokes equations to the
/// advection-diffusion equations because the wind is the fluid velocity,
/// obtained from the source elements in the other mesh
//==========================================================================
template<unsigned DIM>
void RefineableQAdvectionDiffusionElementWithExternalElement<DIM>::
get_wind_adv_diff
(const unsigned& ipt,const Vector<double> &s,const Vector<double>& x, 
 Vector<double>& wind) const
{
 // The interatction is stored at index 0 of the NST element
 unsigned interaction=0;

 // Dynamic cast "other" element to correct type
 NavierStokesEquations<DIM>* source_el_pt=
  dynamic_cast<NavierStokesEquations<DIM>*>
  (external_element_pt(interaction,ipt));
 
 //The wind function is simply the velocity at the points of the source element
 source_el_pt->interpolated_u_nst
  (external_element_local_coord(interaction,ipt),wind);
}  //end of get_wind_adv_diff



//=============start_of_get_dwind==========================================
/// Fill in the derivatives of the wind with respect to the external
/// unknowns in the advection-diffusion equations
//=========================================================================
template<unsigned DIM>
void RefineableQAdvectionDiffusionElementWithExternalElement<DIM>::
get_dwind_adv_diff_dexternal_element_data(const unsigned &ipt,
                                          const unsigned &i,
                                          Vector<double> &result,
                                          Vector<unsigned> &global_eqn_number)
{
 // The interaction index is 0 in this case
 unsigned interaction=0;
 
 // Dynamic cast "other" element to correct type
 NavierStokesEquations<DIM>* source_el_pt=
  dynamic_cast<NavierStokesEquations<DIM>*>
  (external_element_pt(interaction,ipt));
  
 // Get the external element's derivatives of the velocity with respect
 // to the data. The wind is just the Navier--Stokes velocity, so this
 // is all that's required
 source_el_pt->dinterpolated_u_nst_ddata(
  external_element_local_coord(interaction,ipt),i,result,
  global_eqn_number);
} // end_of_get_dwind



//=========================================================================
/// Set the default physical value to be zero in 2D and 3D
//=========================================================================
template<>
double RefineableQCrouzeixRaviartElementWithTwoExternalElement<2>::
Default_Physical_Constant_Value=0.0;

template<>
double RefineableQCrouzeixRaviartElementWithTwoExternalElement<3>::
Default_Physical_Constant_Value=0.0;

/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////


