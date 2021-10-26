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
//
//Driver for a multi-physics problem that couples a Navier--Stokes
//mesh to two advection diffusion meshes, giving double-diffusive
//Boussinesq convection

//Oomph-lib headers, we require the generic, advection-diffusion,
//and navier-stokes elements.
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"

// Both meshes are the standard rectangular quadmesh
#include "meshes/rectangular_quadmesh.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;


//======================class definitions==============================
/// Build QCrouzeixRaviartElementWithTwoExternalElement that 
/// inherits from ElementWithExternalElement
/// so that it can "communicate" with two
/// QAdvectionDiffusionElementWithExternalElement
//=====================================================================
template<unsigned DIM>
class QCrouzeixRaviartElementWithTwoExternalElement : 
 public virtual QCrouzeixRaviartElement<DIM>,
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
 QCrouzeixRaviartElementWithTwoExternalElement() : 
  QCrouzeixRaviartElement<DIM>(),
  ElementWithExternalElement()
  {
   //Set the deafult values of the Rayleigh numbers
   Ra_T_pt = &Default_Physical_Constant_Value;

   Ra_S_pt = &Default_Physical_Constant_Value;

   //There are two interactions
   this->set_ninteraction(2);

   //We do not need to add any external geometric data because the
   //element is fixed
   this->ignore_external_geometric_data();
  } 

 /// Access function for the thermal Rayleigh number (const version)
 const double &ra_t() const {return *Ra_T_pt;}

 /// Access function for the pointer to the thermal Rayleigh number
 double* &ra_t_pt() {return Ra_T_pt;}

  /// Access function for the solutal Rayleigh number (const version)
 const double &ra_s() const {return *Ra_S_pt;}

 /// Access function for the pointer to the solutal Rayleigh number
 double* &ra_s_pt() {return Ra_S_pt;}

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


 /// Compute the element's residual vector and the Jacobian matrix.
 /// Jacobian is computed by finite-differencing or analytically
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
#ifdef USE_FD_JACOBIAN_FOR_MY_NAVIER_STOKES_ELEMENT   
   // This function computes the Jacobian by finite-differencing
   ElementWithExternalElement::
    fill_in_contribution_to_jacobian(residuals,jacobian);
#else
   //Get the contribution from the basic Navier--Stokes element
   QCrouzeixRaviartElement<DIM>::
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
       
       //Loop over the velocity components in the Navier--Stokes equtions
       for(unsigned i=0;i<DIM;i++)
        {
         //If it's not a boundary condition
         local_eqn = this->nodal_local_eqn(l,u_nodal_nst[i]);
         if(local_eqn >= 0)
          {
           //Loop over the external data
           for(unsigned l2=0;l2<n_external_element_data;l2++)
            { 
             //Find the local equation number corresponding to the global
             //unknown
             local_unknown = this->local_eqn_number(
              global_eqn_number_of_external_element_data[l2]);
             if(local_unknown >= 0)
              {
               //Add contribution to jacobian matrix
               jacobian(local_eqn,local_unknown) 
                += dbody_dexternal_element_data(i,l2)*testf(l)*W;
              }
            }
          }
        }
      }
    }
  }

};

//======================class definitions==============================
/// Build QAdvectionDiffusionElementWithExternalElement that inherits from 
/// ElementWithExternalElement
/// so that it can "communicate" with the 
/// QCrouzeixRaviartElementWithExternalElement
//=====================================================================
template<unsigned DIM>
class QAdvectionDiffusionElementWithExternalElement : 
 public virtual QAdvectionDiffusionElement<DIM,3>,
 public virtual ElementWithExternalElement
{

public:

 /// Constructor: call the underlying constructors
 QAdvectionDiffusionElementWithExternalElement() : 
  QAdvectionDiffusionElement<DIM,3>(),
  ElementWithExternalElement()
  { 
   //There is only one interaction
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


 /// Compute the element's residual vector and the Jacobian matrix.
 /// Jacobian is computed by finite-differencing.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
#ifdef USE_FD_JACOBIAN_FOR_MY_ADVECTION_DIFFUSION_ELEMENT   
   // This function computes the Jacobian by finite-differencing
   ElementWithExternalElement::
    fill_in_contribution_to_jacobian(residuals,jacobian);
#else
   //Get the contribution from the basic advection-diffusion element
   QAdvectionDiffusionElement<DIM,3>::
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
 /// degrees of freedom (velocities) on the AdvectionDiffusion equations
 void fill_in_off_diagonal_block_analytic(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
  {
   //Local storage for the  index at which the temperature is stored
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
          this->raw_nodal_value(l,u_nodal_adv_diff)*dpsidx(l,j);
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
         //Assemble the contributions of the velocities 
         //the Advection-Diffusion equations
         
         //If it's not a boundary condition
         local_eqn = this->nodal_local_eqn(l,u_nodal_adv_diff);
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
                interpolated_dudx[i2]*test(l)*W;
              }
            }
          }
        }
      }
    }
  }

};

//============================================================
/// Overload get_body_force_nst to get the temperature "body force"
/// from the "source" AdvectionDiffusions element via 
/// current integration point
//========================================================
template<unsigned DIM>
void QCrouzeixRaviartElementWithTwoExternalElement<DIM>::
get_body_force_nst(const double& time, 
                   const unsigned& ipt,
                   const Vector<double> &s,
                   const Vector<double> &x,
                   Vector<double> &result)
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
}


//=========================================================================
/// Fill in the derivatives of the body force with respect to the external
/// unknowns in the Navier--Stokes equations
//=========================================================================
template<unsigned DIM>
void QCrouzeixRaviartElementWithTwoExternalElement<DIM>::
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

}


//==========================================================================
/// Overload the wind function in the advection-diffusion equations.
/// This provides the coupling from the Navier--Stokes equations to the
/// advection-diffusion equations because the wind is the fluid velocity,
/// obtained from the source elements in the other mesh
//==========================================================================
template<unsigned DIM>
void QAdvectionDiffusionElementWithExternalElement<DIM>::get_wind_adv_diff
(const unsigned& ipt,const Vector<double> &s,const Vector<double>& x, 
 Vector<double>& wind) const
{
 // The interaction index is 0 in this case
 unsigned interaction=0;

 // Dynamic cast "other" element to correct type
 NavierStokesEquations<DIM>* source_el_pt=
  dynamic_cast<NavierStokesEquations<DIM>*>
  (external_element_pt(interaction,ipt));

 //The wind function is simply the velocity at the points of the "other" el
 source_el_pt->interpolated_u_nst
  (external_element_local_coord(interaction,ipt),wind);
}  

//=========================================================================
/// Fill in the derivatives of the wind with respect to the external
/// unknowns in the advection-diffusion equations
//=========================================================================
template<unsigned DIM>
void QAdvectionDiffusionElementWithExternalElement<DIM>::
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
}



//=========================================================================
/// Set the default physical value to be zero
//=========================================================================
template<>
double QCrouzeixRaviartElementWithTwoExternalElement<2>::
Default_Physical_Constant_Value=0.0;


//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{
 /// The Lewis number
 double Lewis = 10.0;

 /// Peclet number (identically one from our non-dimensionalisation)
 /// in both cases
 double Peclet=1.0;

 /// 1/Prandtl number
 double Inverse_Prandtl=1.0;

 /// Thermal Rayleigh number, set to be greater than 
 /// the threshold for linear instability
 double Rayleigh_T = 1800.0;

 /// Solutal Rayleigh number
 double Rayleigh_S = -1000;

 /// Length of domain
 double Lambda = 1.414;

 /// Gravity vector
 Vector<double> Direction_of_gravity(2);
  
} // end_of_namespace

/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////

//====== start_of_problem_class=======================================
/// 2D double-diffusive Convection  problem on three rectangular domains, 
/// discretised 
/// with Navier-Stokes and Advection-Diffusion elements. The specific type
/// of elements is specified via the template parameters.
//====================================================================
template<class NST_ELEMENT,class AD_ELEMENT> 
class DDConvectionProblem : public Problem
{

public:

 /// Constructor
 DDConvectionProblem();

 /// Destructor. Empty
 ~DDConvectionProblem() 
  {
   //Delete the meshes
   delete Conc_mesh_pt;
   delete Temp_mesh_pt;
   delete Nst_mesh_pt;
   //Delete the timestepper
   delete this->time_stepper_pt();
  }

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Actions before adapt:(empty)
 void actions_before_adapt(){}

 /// Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep() 
  {
   set_boundary_conditions(time_pt()->time());
  }

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to specific element and fix pressure (only for NST element...?)
   dynamic_cast<NST_ELEMENT*>(nst_mesh_pt()->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

 /// Doc the solution.
 void doc_solution();

 /// Set the boundary conditions
 void set_boundary_conditions(const double &time);

 /// Access function to the Navier-Stokes mesh
 RectangularQuadMesh<NST_ELEMENT>* nst_mesh_pt() 
  {
   return dynamic_cast<RectangularQuadMesh<NST_ELEMENT>*>(Nst_mesh_pt);
  }

 /// Access function to the Advection-Diffusion mesh
 RectangularQuadMesh<AD_ELEMENT>* temp_mesh_pt() 
  {
   return dynamic_cast<RectangularQuadMesh<AD_ELEMENT>*>(Temp_mesh_pt);
  }
 
 /// Access function to the Advection-Diffusion concentration mesh
 RectangularQuadMesh<AD_ELEMENT>* conc_mesh_pt() 
  {
   return dynamic_cast<RectangularQuadMesh<AD_ELEMENT>*>(Conc_mesh_pt);
  }

/// Get kinetic energy and kinetic energy flux
 void get_kinetic_energy(double &E, double &Edot)
  {
   //Reset values to zero
   E = 0.0; Edot=0.0;
   
   //Loop over the elements
   unsigned n_element = nst_mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     NST_ELEMENT* elem_pt = 
      dynamic_cast<NST_ELEMENT*>(nst_mesh_pt()->element_pt(e));
     
     E += elem_pt->kin_energy();
     Edot += elem_pt->d_kin_energy_dt();
    }
  }


private:
 
 /// DocInfo object
 DocInfo Doc_info;

protected:

 RectangularQuadMesh<NST_ELEMENT>* Nst_mesh_pt;
 RectangularQuadMesh<AD_ELEMENT>* Temp_mesh_pt;
 RectangularQuadMesh<AD_ELEMENT>* Conc_mesh_pt;

}; // end of problem class

//===========start_of_constructor=========================================
/// Constructor for convection problem
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
DDConvectionProblem<NST_ELEMENT,AD_ELEMENT>::DDConvectionProblem()
{

 //Allocate a timestepper
 add_time_stepper_pt(new BDF<2>);

 // Set output directory
 Doc_info.set_directory("RESLT");
 
 // # of elements in x-direction
 unsigned n_x=8;

 // # of elements in y-direction
 unsigned n_y=8;

 // Domain length in x-direction
 double l_x= Global_Physical_Variables::Lambda;

 // Domain length in y-direction
 double l_y=1.0;

 // Build two standard rectangular quadmesh
 Nst_mesh_pt = 
  new RectangularQuadMesh<NST_ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());
 Temp_mesh_pt = 
  new RectangularQuadMesh<AD_ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());
 Conc_mesh_pt = 
  new RectangularQuadMesh<AD_ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here

 //Loop over the boundaries
 unsigned num_bound = nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Loop over the number of nodes on the boundry
   unsigned num_nod= nst_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //If we are on the side-walls, the v-velocity 
     //satisfies natural boundary conditions, so we only pin the
     //first value
     if ((ibound==1) || (ibound==3)) 
      {
       nst_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
      }
     else // On the top and bottom walls, we have "stress-free" conditions
      //which actually corresponds to transverse stress free and normal 
      //zero velocity (symmetry)
      //Thus we pin the second value
      {
       nst_mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
    }
  }

 //Pin the zero-th pressure dof in element 0 and set its value to
 //zero:
 fix_pressure(0,0,0.0);

 //Loop over the boundaries of the AD mesh
 num_bound = temp_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Loop over the number of nodes on the boundry
   unsigned num_nod= temp_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //If we are on the side-walls, the temperature
     //satisfies natural boundary conditions, so we don't pin anything
     // in this mesh
     if ((ibound==1) || (ibound==3)) 
      {
      
      }
     //Otherwise pin the temperature
     else // pin all values
      {
       temp_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
      }
    }
  }


 //Loop over the boundaries of the AD mesh
 num_bound = conc_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Loop over the number of nodes on the boundry
   unsigned num_nod= conc_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //If we are on the side-walls, the concentration
     //satisfies natural boundary conditions, so we don't pin anything
     // in this mesh
     if ((ibound==1) || (ibound==3)) 
      {

      }
     //Otherwiwse pin the concentration
     else // pin all values
      {
       conc_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
      }
    }
  }


 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructors. 
 unsigned n_nst_element = nst_mesh_pt()->nelement();
 for(unsigned i=0;i<n_nst_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   NST_ELEMENT *el_pt = dynamic_cast<NST_ELEMENT*>
    (nst_mesh_pt()->element_pt(i));

   // Set the Reynolds number (1/Pr in our non-dimensionalisation)
   el_pt->re_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set the Thermal Rayleigh number
   el_pt->ra_t_pt() = &Global_Physical_Variables::Rayleigh_T;

   // Set the Solutal Rayleigh number
   el_pt->ra_s_pt() = &Global_Physical_Variables::Rayleigh_S;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;

   //The mesh is fixed, so we can disable ALE
   el_pt->disable_ALE();
  }

 unsigned n_temp_element = temp_mesh_pt()->nelement();
 for(unsigned i=0;i<n_temp_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   AD_ELEMENT *el_pt = dynamic_cast<AD_ELEMENT*>
    (temp_mesh_pt()->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;

   // Set the timescale to be the same as the Navier--Stokes
   // equations (1.0)
   el_pt->pe_st_pt() =&Global_Physical_Variables::Peclet;

   //The mesh is fixed, so we can disable ALE
   el_pt->disable_ALE();
  }

 unsigned n_conc_element = conc_mesh_pt()->nelement();
 for(unsigned i=0;i<n_conc_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   AD_ELEMENT *el_pt = dynamic_cast<AD_ELEMENT*>
    (conc_mesh_pt()->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Lewis;

   // Set the Peclet number multiplied by the Strouhal number
   el_pt->pe_st_pt() =&Global_Physical_Variables::Lewis;

   //The mesh is fixed, so we can disable ALE
   el_pt->disable_ALE();
  }


 // Set sources for temperature
 Multi_domain_functions::setup_multi_domain_interactions
  <NST_ELEMENT,AD_ELEMENT>(this,nst_mesh_pt(),temp_mesh_pt());

 // Set sources for concentrations
 Multi_domain_functions::setup_multi_domain_interactions
  <NST_ELEMENT,AD_ELEMENT>(this,nst_mesh_pt(),conc_mesh_pt(),1,0);
 
 // combine the submeshes
 add_sub_mesh(Nst_mesh_pt);
 add_sub_mesh(Temp_mesh_pt);
 add_sub_mesh(Conc_mesh_pt);
 build_global_mesh();

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor



//===========start_of_set_boundary_conditions================
/// Set the boundary conditions as a function of continuous 
/// time
//===========================================================
template<class NST_ELEMENT,class AD_ELEMENT>
void DDConvectionProblem<NST_ELEMENT,AD_ELEMENT>::set_boundary_conditions(
 const double &time)
{
 // Loop over all the boundaries on the NST mesh
 unsigned num_bound=nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=nst_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=nst_mesh_pt()->boundary_node_pt(ibound,inod);

     //Set the number of velocity components
     unsigned vel_max=2;

     //If we are on the side walls we only set the x-velocity.
     if((ibound==1) || (ibound==3)) {vel_max = 1;}

     //Set the pinned velocities to zero on NST mesh
     for(unsigned j=0;j<vel_max;j++) {nod_pt->set_value(j,0.0);}
    }
  }

 // Loop over all the boundaries on the AD mesh
 num_bound=temp_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=temp_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=temp_mesh_pt()->boundary_node_pt(ibound,inod);
     
     //If we are on the top boundary, set the temperature 
     //to -0.5 (cooled)
     if(ibound==2) {nod_pt->set_value(0,-0.5);}

     //If we are on the bottom boundary, set the temperature
     //to 0.5 (heated)
     if(ibound==0) {nod_pt->set_value(0,0.5);}
    }
  }



 // Loop over all the boundaries on the AD mesh
 num_bound=conc_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=conc_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=conc_mesh_pt()->boundary_node_pt(ibound,inod);

     //If we are on the top boundary, set the concentration to be low
     //to -0.5 (cooled)
     if(ibound==2) 
      {
       nod_pt->set_value(0,-0.5);
       
       //Add small concentration imperfection if desired
       double epsilon = 0.01;
       
       //Read out the x position
       double x = nod_pt->x(0);
       
       //Set a sinusoidal perturbation in the concentration
       double value = sin(2.0*MathematicalConstants::Pi*x/1.5)*
        epsilon*time*exp(-10.0*time);
       nod_pt->set_value(0, -0.5 + value);
      }

     //If we are on the bottom boundary, set the concentration to be high
     //to 0.5 (heated)
     if(ibound==0) {nod_pt->set_value(0,0.5);}
    }
  }


} // end_of_set_boundary_conditions

//===============start_doc_solution=======================================
/// Doc the solution
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
void DDConvectionProblem<NST_ELEMENT,AD_ELEMENT>::doc_solution()
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output whole solution (this will output elements from one mesh
 //----------------------  followed by the other mesh at the moment...?)
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

 Doc_info.number()++;
} // end of doc


//=======start_of_main================================================
/// Driver code for 2D Boussinesq convection problem
//====================================================================
int main(int argc, char **argv)
{
 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;

 //Construct our problem
 DDConvectionProblem<QCrouzeixRaviartElementWithTwoExternalElement<2>,
  QAdvectionDiffusionElementWithExternalElement<2> > 
  problem;

 // Apply the boundary condition at time zero
 problem.set_boundary_conditions(0.0);

 //Perform a single steady Newton solve
 problem.steady_newton_solve();

 //Document the solution
 problem.doc_solution();

  //Start a trace file
 ofstream trace("RESLT/trace.dat");
 //Local variables for the kinetic energy and its rate of change
 double E=0.0, Edot = 0.0;
 
 //Output to the trace file
 problem.get_kinetic_energy(E,Edot);
 trace << problem.time_pt()->time() << " " 
       << problem.nst_mesh_pt()->boundary_node_pt(1,8)->value(1) << " " 
       << E << " " << Edot << std::endl;
 
 //Set the timestep
 double dt = 0.01;

 //Initialise the value of the timestep and set an impulsive start
 problem.assign_initial_values_impulsive(dt);

 //Set the number of timesteps to our default value
 unsigned n_steps = 2000;

 //If we have a command line argument, perform fewer steps 
 //(used for self-test runs)
 if(argc > 1) {n_steps = 5;}

 //Perform n_steps timesteps
 for(unsigned i=0;i<n_steps;++i)
  {
   problem.unsteady_newton_solve(dt);
   problem.doc_solution();

 //Output to the trace file
   problem.get_kinetic_energy(E,Edot);
   trace << problem.time_pt()->time() << " " 
         << problem.nst_mesh_pt()->boundary_node_pt(1,8)->value(1) << " "
         << E << " " << Edot << std::endl;
  }

 trace.close();
} // end of main









