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
//Driver for a multi-physics problem that couples the Navier--Stokes
//equations to the advection diffusion equations and adds
//a temperature-dependent surface tension to give
//Bernard-Marangoni convection. The effects of a soluble surfactant
//are also included

//N.B. Surfactant calculations require validation!
//The diffusivity ratio is chosen to be large so that something different
//happens, i.e. surfactant diffuses back into the bulk and this
//prevents the oscillatory instability (maybe!)

//Oomph-lib headers, 
//We require the generic header
#include "generic.h"
//Our custom coupling of advection-diffusion and Navier--Stokes
#include "double_buoyant_navier_stokes_elements.h"
//The fluid interface elements
#include "fluid_interface.h"

// The mesh is our standard rectangular quadmesh
#include "meshes/single_layer_spine_mesh.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;

namespace oomph
{

//==================================================================
/// Spine-based Marangoni surface tension elements that add
/// a linear dependence on temperature and concentration
/// of a surface chemical to the surface tension, 
/// which decreases with increasing temperature. In addition, 
/// this element adds a flux contribution to the advection-diffusion
/// equation to represent heat loss at the free surface. This
/// introduces the Biot number.
//=================================================================
template<class ELEMENT>
class SpineLineMarangoniSurfactantFluidInterfaceElement :
public SpineLineFluidInterfaceElement<ELEMENT>
{
private:

 /// Pointer to a Biot number
 double *Bi_pt;

 /// Pointer to a Marangoni number
 double *Ma_pt;

 /// Pointer to an Elasticity number
 double *Beta_pt;

 /// Pointer to Surface Peclet number
 double *Peclet_S_pt;

 /// Pointer to the surface Peclect Strouhal number
 double *Peclet_Strouhal_S_pt;

 /// Pointer to the diffusion ratios
 double *D_pt;

 /// Pointer to the reaction ratios
 double *K_pt;

 /// Index at which the temperature is stored at the nodes
 unsigned T_index;

 /// Index at which the bulk concentration is stored at the nodes
 unsigned C_bulk_index;

 /// Index at which the surfactant concentration is stored at the
 /// nodes
 unsigned C_index;

 /// Default value of the physical constants
 static double Default_Physical_Constant_Value;

protected:

 /// Get the temperature
 double interpolated_T(const Vector<double> &s)
  {
     //Find number of nodes
   unsigned n_node = this->nnode();

   //Get the nodal index at which the unknown is stored
   const unsigned t_index = T_index;

   //Local shape function
   Shape psi(n_node);

   //Find values of shape function
   this->shape(s,psi);

   //Initialise value of t
   double T = 0.0;

   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_node;l++) 
    {
     T += this->nodal_value(l,t_index)*psi(l);
    }

   return(T);
  }


 /// Get the surfactant concentration
 double interpolated_C(const Vector<double> &s)
  {
     //Find number of nodes
   unsigned n_node = this->nnode();

   //Get the nodal index at which the unknown is stored
   const unsigned c_index = C_index;

   //Local shape function
   Shape psi(n_node);

   //Find values of shape function
   this->shape(s,psi);

   //Initialise value of C
   double C = 0.0;

   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_node;l++) 
    {
     C += this->nodal_value(l,c_index)*psi(l);
    }

   return(C);
  }


 /// The time derivative of the surface concentration
 double dcdt_surface(const unsigned &l) const
  {
   // Get the data's timestepper
   TimeStepper* time_stepper_pt= this->node_pt(l)->time_stepper_pt();
   
   //Initialise dudt
   double dcdt=0.0;
   //Loop over the timesteps, if there is a non Steady timestepper
   if (time_stepper_pt->type()!="Steady")
    {
     //Find the index at which the variable is stored
     const unsigned c_index = C_index;

     // Number of timsteps (past & present)
     const unsigned n_time = time_stepper_pt->ntstorage();
     
     for(unsigned t=0;t<n_time;t++)
      {
       dcdt += time_stepper_pt->weight(1,t)*this->nodal_value(t,l,c_index);
      }
    }
   return dcdt;
  }

 /// The surface tension function is linear in the
 /// temperature with constant of proportionality equal
 /// to the Marangoni number.
 double sigma(const Vector<double> &s)
  {
   //Find the number of shape functions
   const unsigned n_node = this->nnode();
   //Now get the shape fuctions at the local coordinate
   Shape psi(n_node);
   this->shape(s,psi);

   //Now interpolate the temperature and surfactant concentration
   double T = 0.0, C=0.0;
   for(unsigned l=0;l<n_node;l++)
    {
     T += this->nodal_value(l,T_index)*psi(l);
     C += this->nodal_value(l,C_index)*psi(l);
    }
   
   //Get the Marangoni, Capillary and Elasticity numbers
   double Ma = this->ma();
   double Ca = this->ca();
   double Beta = this->beta();
   //Return the variable surface tension
   //The additional multiplication by Ca will cancel with the 1/Ca
   //in the underlying equations
   return (1.0 - Ca*Ma*T - Ca*Beta*(C-1.0));
  }

 /// Fill in the contribution to the residuals
  /// Calculate the contribution to the jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                                       DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   this->fill_in_generic_residual_contribution_interface(residuals,jacobian,1);
   {
    //Use finite differences to handle temperature variations
    const unsigned n_node = this->nnode();
    //Find the number of dofs in the element
    const unsigned n_dof = this->ndof();
    //Create newres vector
    Vector<double> newres(n_dof);
    
    //Integer storage for local unknown
    int local_unknown=0;
    
    //Use the default finite difference step
    const double fd_step = this->Default_fd_jacobian_step;
    
    //Loop over the nodes
    for(unsigned n=0;n<n_node;n++)
     {
      //Get the number of values stored at the node
      unsigned t_index = this->T_index;

      //Get the local equation number
      local_unknown = this->nodal_local_eqn(n,t_index);
      //If it's not pinned
      if(local_unknown >= 0)
       {
        //Store a pointer to the nodal data value
        double *value_pt = this->node_pt(n)->value_pt(t_index);
        
        //Save the old value of the Nodal data
        double old_var = *value_pt;
       
        //Increment the value of the Nodal data
        *value_pt += fd_step;
       
        //Calculate the new residuals
        this->get_residuals(newres);
       
        //Do finite differences
        for(unsigned m=0;m<n_dof;m++)
         {
          double sum = (newres[m] - residuals[m])/fd_step;
          //Stick the entry into the Jacobian matrix
          jacobian(m,local_unknown) = sum;
         }
       
        //Reset the Nodal data
        *value_pt = old_var;
       }
     }
    
    //Use finite differences to handle the concentration variations
    //Loop over the nodes again
    for(unsigned n=0;n<n_node;n++)
     {
      //Get the number of values stored at the node
      unsigned c_index = this->C_index;

      //Get the local equation number
      local_unknown = this->nodal_local_eqn(n,c_index);
      //If it's not pinned
      if(local_unknown >= 0)
       {
        //Store a pointer to the nodal data value
        double *value_pt = this->node_pt(n)->value_pt(c_index);
        
        //Save the old value of the Nodal data
        double old_var = *value_pt;
       
        //Increment the value of the Nodal data
        *value_pt += fd_step;
       
        //Calculate the new residuals
        this->get_residuals(newres);
       
        //Do finite differences
        for(unsigned m=0;m<n_dof;m++)
         {
          double sum = (newres[m] - residuals[m])/fd_step;
          //Stick the entry into the Jacobian matrix
          jacobian(m,local_unknown) = sum;
         }
       
        //Reset the Nodal data
        *value_pt = old_var;
       }
     }
   }

   //Call the generic routine to handle the spine variables
   SpineElement<FaceGeometry<ELEMENT> >::
    fill_in_jacobian_from_geometric_data(jacobian);
  }

 
 /// Overload the Helper function to calculate the residuals and 
 /// jacobian entries. This particular function ensures that the
 /// additional entries are calculated inside the integration loop
 void add_additional_residual_contributions_interface(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,
  const unsigned &flag,
  const Shape &psif, const DShape &dpsifds,
  const DShape &dpsifdS, const DShape &dpsifds_div,
  const Vector<double> &s,
  const Vector<double> &interpolated_x, 
  const Vector<double> &interpolated_n, 
  const double &W,
  const double &J)
  {
   //Find the index at which the temperature is stored
   unsigned t_index = this->T_index;
   
   //Find out how many nodes there are
   unsigned n_node = this->nnode();
   
   //Read out the Bi number
   const double Bi = this->bi();

   //Now calculate the temperature at this point
   //Assuming the same shape functions are used (which they are)
   double T = 0.0;
   for(unsigned l=0;l<n_node;l++)
    {
     T += this->nodal_value(l,t_index)*psif(l);
    }

   //Storage for the local equation numbers and unknowns
   int local_eqn = 0, local_unknown = 0;

   //Now we add the flux term to the appropriate residuals
   for(unsigned l=0;l<n_node;l++)
    {
     //Read out the apprporiate local equation
     local_eqn = this->nodal_local_eqn(l,t_index);

     //If not a boundary condition
     if(local_eqn >= 0)
      {
       residuals[local_eqn] -= (1.0 + Bi*(T-1.0))*psif(l)*W*J;

       //We also need to worry about the jacobian terms
       if(flag)
        {
         //Loop over the nodes again
         for(unsigned l2=0;l2<n_node;l2++)
          {
           //Get the unknown
           local_unknown = this->nodal_local_eqn(l2,t_index);
           //If not a boundary condition
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -= Bi*psif(l2)*psif(l)*W*J;
            }
          }
        }
      }
    } //End of loop over the nodes


   //Bulk flux condition
   //Find the index at which the bulk is stored
   unsigned c_bulk_index = this->C_bulk_index;

   //Find the index at which the surface surfactant
   unsigned c_index = this->C_index;

   
   //Now calculate the bulk concentration at this point
   //Assuming the same shape functions are used (which they are)
   double C_bulk = 0.0;
   double C = 0.0;
   for(unsigned l=0;l<n_node;l++)
    {
     C_bulk += this->nodal_value(l,c_bulk_index)*psif(l);
     C += this->nodal_value(l,c_index)*psif(l);
    }

   //Get the reaction ratio
   const double K = this->k();

   //The transport between the two layers is given by the flux
   double flux = K*C_bulk - C;

   //Now we add the flux term to the appropriate residuals
   for(unsigned l=0;l<n_node;l++)
    {
     //Read out the apprporiate local equation
     local_eqn = this->nodal_local_eqn(l,c_bulk_index);

     //If not a boundary condition
     if(local_eqn >= 0)
      {
       residuals[local_eqn] -= flux*psif(l)*W*J;

       //We also need to worry about the jacobian terms
       if(flag)
        {
         //Loop over the nodes again
         for(unsigned l2=0;l2<n_node;l2++)
          {
           //Get the unknown
           local_unknown = this->nodal_local_eqn(l2,c_bulk_index);
           //If not a boundary condition
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -= K*psif(l2)*psif(l)*W*J;
            }
           
           local_unknown = this->nodal_local_eqn(l2,c_index);
           //If not a boundary condition
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += psif(l2)*psif(l)*W*J;
            }
          }
        }
      }
    } //End of loop over the nodes



   //Surface advection-diffusion equation
   //Get the diffusion ratio
   const double D = d();
   
    //Find the index at which the temperature is stored
   //unsigned c_index = this->C_index;
   Vector<unsigned> u_index = this->U_index_interface;

   //Read out the surface peclect number
   const double Pe_s = this->peclet_s();
   const double PeSt_s = this->peclet_strouhal_s();

   //Now calculate the concentration at this point
   //Assuming the same shape functions are used (which they are)
   double interpolated_C = 0.0;
   double interpolated_dCds = 0.0;
   double dCdt = 0.0;
   //The tangent vector
   const unsigned ndim = this->node_pt(0)->ndim();
   Vector<double> interpolated_tangent(ndim,0.0);
   Vector<double> interpolated_u(ndim,0.0);
   Vector<double> interpolated_duds(ndim,0.0);
   Vector<double> mesh_velocity(ndim,0.0);

   if(ndim != u_index.size())
    {
     throw OomphLibError("Dimension Incompatibility",
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }

   for(unsigned l=0;l<n_node;l++)
    {
     const double psi = psif(l);
     const double dpsi = dpsifds(l,0);
     interpolated_C += this->nodal_value(l,c_index)*psi;
     interpolated_dCds += this->nodal_value(l,c_index)*dpsi;
     dCdt += dcdt_surface(l)*psi;
     for(unsigned i=0;i<ndim;i++)
      {
       interpolated_tangent[i] += this->nodal_position(l,i)*dpsi;
       interpolated_u[i] += this->nodal_value(l,u_index[i])*psi;
       interpolated_duds[i] += this->nodal_value(l,u_index[i])*dpsi;
       mesh_velocity[i] += this->dnodal_position_dt(l,i)*psif(l);
      }
    }

   double u_tangent = 0.0, t_length = 0.0;
   for(unsigned i=0;i<ndim;i++) 
    {
     u_tangent += 
      (interpolated_u[i] - mesh_velocity[i])*interpolated_tangent[i];
     t_length  += interpolated_tangent[i]*interpolated_tangent[i];
    }

   //Now we add the flux term to the appropriate residuals
   for(unsigned l=0;l<n_node;l++)
    {
     //Read out the apprporiate local equation
     local_eqn = this->nodal_local_eqn(l,c_index);

     //If not a boundary condition
     if(local_eqn >= 0)
      {
       //Time derivative term
       residuals[local_eqn] -= PeSt_s*dCdt*psif(l)*W*J;

       //Diffusion term
       residuals[local_eqn] -= interpolated_dCds*dpsifds(l,0)*W/J;
       
       //Advection term in new formulation
       residuals[local_eqn] -= 
        Pe_s*(u_tangent*interpolated_dCds + 
              interpolated_C*
              (interpolated_tangent[0]*interpolated_duds[0]
               + interpolated_tangent[1]*interpolated_duds[1]))*psif(l)*W/J;

       //Now add the flux term
       residuals[local_eqn] += D*flux*psif(l)*W*J;
       //We also need to worry about the jacobian terms
       if(flag)
        {
         //Loop over the nodes again
         for(unsigned l2=0;l2<n_node;l2++)
          {
           //Loop over the velocity components
           for(unsigned i2=0;i2<ndim;i2++)
            {
             //Get the unknown
             local_unknown = this->nodal_local_eqn(l2,u_index[i2]);
             //If not a boundary condition
             if(local_unknown >= 0)
              {
               jacobian(local_eqn,local_unknown) 
                -= Pe_s*interpolated_tangent[i2]*
                (psif(l2)*interpolated_dCds + 
                 interpolated_C*dpsifds(l2,0))*psif(l)*W/J;
              }
            }
           
           //Get the unknown
           local_unknown = this->nodal_local_eqn(l2,c_bulk_index);
           //If not a boundary condition
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += D*K*psif(l2)*psif(l)*W*J;
            }
           
           local_unknown = this->nodal_local_eqn(l2,c_index);
           //If not a boundary condition
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -= D*psif(l2)*psif(l)*W*J;
            }
          }
        }
      }
    } //End of loop over the nodes

  }  

 
 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   //Add the contribution to the jacobian
   this->fill_in_contribution_to_jacobian(residuals,jacobian);
   //No mass matrix terms, but should probably do kinematic bit here
  }

public:
 /// Constructor that passes the bulk element and face index
 /// down to the underlying 
 SpineLineMarangoniSurfactantFluidInterfaceElement(
 FiniteElement* const &element_pt, const int &face_index) : 
  SpineLineFluidInterfaceElement<ELEMENT>
  (element_pt,face_index)
  {
   //Initialise the values
   Bi_pt = &Default_Physical_Constant_Value;
   Ma_pt = &Default_Physical_Constant_Value;
   Beta_pt = &Default_Physical_Constant_Value;
   Peclet_S_pt = &Default_Physical_Constant_Value;
   Peclet_Strouhal_S_pt = &Default_Physical_Constant_Value;
   K_pt = &Default_Physical_Constant_Value;
   D_pt = &Default_Physical_Constant_Value;

   //Cast the bulk element 
   ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(element_pt);
   //Now find the index at which the temperature is stored from the 
   //advection-diffusion part of the bulk element
   T_index = cast_element_pt->c_index_adv_diff_react(0);
   
   //Find the index at which the bulk surfactant concentration is stored
   C_bulk_index = cast_element_pt->c_index_adv_diff_react(1);

   //Add the additional surfactant terms to these surface elements
   
   //Read out the number of nodes on the face
   //For some reason I need to specify the this pointer here(!)
   unsigned n_node_face = this->nnode();
   //Set the additional data values in the face
   //There is one additional values at each node --- the lagrange multiplier
   Vector<unsigned> additional_data_values(n_node_face);
   for(unsigned i=0;i<n_node_face;i++) additional_data_values[i] = 1;
   //Resize the data arrays accordingly 
   this->resize_nodes(additional_data_values);

   //The C_index is the new final value
   //HACK HERE
   C_index = this->node_pt(0)->nvalue()-1;
  }

 /// Return the Biot number
 double bi() {return *Bi_pt;}
 
 /// Return the Marangoni number
 double ma() {return *Ma_pt;}

 /// Return the Elasticity number
 double beta() {return *Beta_pt;}

 /// Return the surface peclect number
 double peclet_s() {return *Peclet_S_pt;}

 /// Return the surface peclect strouhal number
 double peclet_strouhal_s() {return *Peclet_Strouhal_S_pt;}

 /// Return the diffusion ratio
 double d() {return *D_pt;}

 /// Return the reaction ratio
 double k() {return *K_pt;}

 /// Access function for pointer to the Marangoni number
 double* &ma_pt() {return Ma_pt;}

 /// Access function for pointer to the Biot number
 double* &bi_pt() {return Bi_pt;}

 /// Access function for pointer to the Elasticity number
 double* &beta_pt() {return Beta_pt;}

 /// Access function for pointer to the surface Peclet number
 double* &peclet_s_pt() {return Peclet_S_pt;}

 /// Access function for pointer to the surface Peclet x Strouhal number
 double* &peclet_strouhal_s_pt() {return Peclet_Strouhal_S_pt;}

 /// Access function for pointer to the diffusion ratios
 double* &d_pt() {return D_pt;}

 /// Access function for pointer to the reaction ratios
 double* &k_pt() {return K_pt;}


void output(std::ostream &outfile, const unsigned &n_plot)
{
 //Set output Vector
 Vector<double> s(1);
 
 //Tecplot header info 
 outfile << "ZONE I=" << n_plot << std::endl;
 
 //Loop over plot points
 for(unsigned l=0;l<n_plot;l++)
  {
   s[0] = -1.0 + l*2.0/(n_plot-1);
   
   //Output the x,y,u,v 
   for(unsigned i=0;i<2;i++) outfile << this->interpolated_x(s,i) << " ";
   for(unsigned i=0;i<2;i++) outfile << this->interpolated_u(s,i) << " ";      
   //Output a dummy pressure
   outfile << 0.0 << " ";
   //Output the temperature
   outfile << interpolated_T(s) << " " 
           << interpolated_C(s) << std::endl;
  }
 outfile << std::endl;
}



};


//Define the default physical value to be one
template<class ELEMENT>
double SpineLineMarangoniSurfactantFluidInterfaceElement<ELEMENT>::
Default_Physical_Constant_Value = 1.0;

}

//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{
 /// Diffusivity  (identically one from our non-dimensionalisation)
 Vector<double> D(2,1.0);

 /// Timescales for transport equations (identically one from our
 /// non-dimensionalisation)
 Vector<double> Tau(2,1.0);

 /// 1/Prandtl number
 double Inverse_Prandtl=0.1;//1.0;

 /// Rayleigh number, set to be zero so that
 /// there are no gravitational effects
 double Rayleigh = 0.0;

 /// Scaled Bond number (Bo/Ca), set to be zero
 /// so that there are no gravitational effects
 double Scaled_Bond = 0.0;
 
 /// Biot number
 double Biot = 0.0; //1.0;

 /// Marangoni number (just above the threshold for 
 /// linear instability)
 double Marangoni = 125.0;

 /// Capillary number (of which the results are independent
 /// for a pinned surface)
 double Capillary = 0.0045;

 /// Surface Elasticity number
 double Beta = 100.0;//0.1;//100.0;//1.0;//1.0e5;

 /// Surface Peclet number
 double Peclet_S = 1.0;

 /// \shorT Sufrace Peclet number multiplied by Strouhal number
 double Peclet_St_S = 100.0;

 /// The ratio of adsorption-desorption times
 double K = 1.0; //1.0;

 /// The ratio of bulk diffusion to surface diffusion
 double DD_s = 1000.0;

 /// Gravity vector
 Vector<double> Direction_of_gravity(2);
  
} // end_of_namespace

/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////

//====== start_of_problem_class=======================================
/// 2D Convection  problem on rectangular domain, discretised 
/// with refineable elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT,class INTERFACE_ELEMENT> 
class ConvectionProblem : public Problem
{

public:

 /// Constructor. The boolean indicates whether the free surface
 //should be pinned or not in the first instance
 ConvectionProblem(const bool &pin=true);

 /// Destructor. Empty
 ~ConvectionProblem() {}

 /// Release the free surface so that it can move
 void unpin_surface()
  {
   //Only bother if the surface is pinned
   if(Surface_pinned)
    {
     Surface_pinned = false;
     
     //Unpin the heights of all the spines in the middle
     unsigned n_spine = Bulk_mesh_pt->nspine();
     for(unsigned n=0;n<n_spine;n++)
      {
       Bulk_mesh_pt->spine_pt(n)->spine_height_pt()->unpin(0);
      }
     
     //If we on the top wall, v velocity is no longer pinned
     unsigned ibound=2;
     //Loop over the number of nodes on the boundary
     unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //Loop over the desired values stored at the nodes and unpin
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->unpin(1);
      }
     
     //Unfix the pressure
     unfix_pressure(0,0);

     // Loop over the elements to set up element-specific 
     // and re-enable ALE
     unsigned n_element = Bulk_mesh_pt->nelement();
     for(unsigned i=0;i<n_element;i++)
      {
       // Upcast from GeneralsedElement to the present element
       ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

       el_pt->enable_ALE();
      }

     //Reassign the equation number
     std::cout << "Surface unpinned to give " 
               << assign_eqn_numbers() << " equation numbers\n";
    }
  }


 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Remember to update the nodes if the surface is not pinned
 void actions_before_newton_convergence_check()
  {if(!Surface_pinned) {Bulk_mesh_pt->node_update();}}

 /// Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep() 
  {set_boundary_conditions(time_pt()->time());}

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure


 /// UnFix pressure in element e at pressure dof pdof and set to pvalue
 void unfix_pressure(const unsigned &e, const unsigned &pdof)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    unfix_pressure(pdof);
  } // end_of_unfix_pressure


 /// Doc the solution.
 void doc_solution(std::ofstream &trace);

 /// Set the boundary conditions
 void set_boundary_conditions(const double &time);

 /// Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 SingleLayerSpineMesh<ELEMENT>* Bulk_mesh_pt;

 Mesh* Surface_mesh_pt;
 
private:
 
 /// DocInfo object
 DocInfo Doc_info;

 /// Boolean to indicate whether the surface is pinned
 bool Surface_pinned;

}; // end of problem class

//===========start_of_constructor=========================================
/// Constructor for convection problem
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
ConvectionProblem<ELEMENT,INTERFACE_ELEMENT>::
ConvectionProblem(const bool &pin) : Surface_pinned(pin)
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
 double l_x=3.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build a standard rectangular quadmesh
 Bulk_mesh_pt = 
  new SingleLayerSpineMesh<ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());



 //Create the surface mesh that will contain the interface elements
 //First create storage, but with no elements or nodes
 Surface_mesh_pt = new Mesh;

 //Loop over the horizontal elements
 for(unsigned i=0;i<n_x;i++)
  {
   //Construct a new 1D line element on the face on which the local
   //coordinate 1 is fixed at its max. value (1) --- This is face 2
   FiniteElement *interface_element_pt =
    new INTERFACE_ELEMENT(
     Bulk_mesh_pt->finite_element_pt(n_x*(n_y-1)+i),2);
   
   //Push it back onto the stack
   this->Surface_mesh_pt->add_element_pt(interface_element_pt); 
  }
 // Add the two sub-meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all sub-meshes into a single mesh
 build_global_mesh();

 
 //Pin the heights of all the spines if the surface is pinned
 if(Surface_pinned)
  {
   unsigned n_spine = Bulk_mesh_pt->nspine();
   for(unsigned n=0;n<n_spine;n++)
    {
     Bulk_mesh_pt->spine_pt(n)->spine_height_pt()->pin(0);
    }
  }

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here

 //Loop over the boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the minimum index to be pinned (all values by default)
   unsigned val_min=0;
   //Set the maximum index to be pinned (all values by default)
   unsigned val_max=4;
   //If we are on the side-walls, the v-velocity, temperature
   //and concentration
   //satisfy natural boundary conditions, so we only pin the
   //first value
   if((ibound==1) || (ibound==3)) {val_max=1;}

   //If we on the top wall, v velocity is pinned
   if(ibound==2) 
    {
     //If the surface is pinned, pin the v velocity
     if(Surface_pinned) {val_min=1; val_max=2;}
     //Otherwise pin nothing
     else {val_min=0; val_max=0;}
    }

   //If we are on the lower wall, concentration is free (no flux)
   if(ibound==0) {val_max = 3;}

   //Loop over the number of nodes on the boundary
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //Loop over the desired values stored at the nodes and pin
     for(unsigned j=val_min;j<val_max;j++)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(j);
      }
    }
  }

 //Pin the zero-th pressure dof in element 0 and set its value to
 //zero:
 if(Surface_pinned) {fix_pressure(0,0,0.0);}

 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor.
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

   // Set the diffusivities number
   el_pt->diff_pt() = &Global_Physical_Variables::D;

   // Set the timescales
   el_pt->tau_pt() =&Global_Physical_Variables::Tau;

   // Set the Reynolds number (1/Pr in our non-dimensionalisation)
   el_pt->re_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set the Re/Fr equal to Bo/Ca, the scaled Bond number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::Scaled_Bond;

   // Set the Rayleigh number
   el_pt->ra_t_pt() = &Global_Physical_Variables::Rayleigh;

   el_pt->ra_s_pt() = &Global_Physical_Variables::Rayleigh;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;

   //If the mesh is fixed, we can disable ALE
   if(Surface_pinned) {el_pt->disable_ALE();}
  }


  // Loop over the interface elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor.
 unsigned n_interface  = Surface_mesh_pt->nelement();
 for(unsigned i=0;i<n_interface;i++)
  {
   // Upcast from GeneralsedElement to the present element
   INTERFACE_ELEMENT *el_pt = dynamic_cast<INTERFACE_ELEMENT*>(
    Surface_mesh_pt->element_pt(i));
   
   // Set the Biot number
   el_pt->bi_pt() = &Global_Physical_Variables::Biot;

   // Set the Marangoni number
   el_pt->ma_pt() =&Global_Physical_Variables::Marangoni;

   // Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Capillary;

   // Set the surface elasticity number
   el_pt->beta_pt() = &Global_Physical_Variables::Beta;

   // Set the surface peclect number
   el_pt->peclet_s_pt() = &Global_Physical_Variables::Peclet_S;

   // Set the surface peclect number multiplied by strouhal number
   el_pt->peclet_strouhal_s_pt() = &Global_Physical_Variables::Peclet_St_S;

   // Set the reaction ratio
   el_pt->k_pt() = &Global_Physical_Variables::K;

   // Set the diffustion ratio
   el_pt->d_pt() = &Global_Physical_Variables::DD_s;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor



//===========start_of_set_boundary_conditions================
/// Set the boundary conditions as a function of continuous 
/// time
//===========================================================
template<class ELEMENT,class INTERFACE_ELEMENT>
void ConvectionProblem<ELEMENT,INTERFACE_ELEMENT>::set_boundary_conditions(
 const double &time)
{
 //Set initial temperature profile
 if(time <= 0.0)
  {
   unsigned n_node = Bulk_mesh_pt->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt = Bulk_mesh_pt->node_pt(n);
     //Set linear variation
     nod_pt->set_value(2,2.0 - nod_pt->x(1));
     //And in uniformly distributed surfactant
     nod_pt->set_value(3,1.0/Global_Physical_Variables::K);
    }

   //Set the initial surface concentration to be one
   unsigned ibound = 2;
   n_node = Bulk_mesh_pt->nboundary_node(ibound);
   for(unsigned n=0;n<n_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(ibound,n)->set_value(4,1.0);
    }
  }

 // Loop over the boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

     //Set the number of velocity components
     unsigned vel_max=2;

     //If we are on the side walls we only set the x-velocity.
     if((ibound==1) || (ibound==3)) {vel_max = 1;}

     //If we are on the top boundary, do not set the velocities
     //(yet)
     if(ibound==2) {vel_max = 0;}

     //Set the pinned velocities to zero
     for(unsigned j=0;j<vel_max;j++) {nod_pt->set_value(j,0.0);}

     //If we are on the bottom boundary peturb the velocity
     if(ibound==0) //2 
      {
       //Add small velocity imperfection if desired
       double epsilon = 0.01;

       //Read out the x position
       double x = nod_pt->x(0);

       //Set a sinusoidal perturbation in the vertical velocity
       //This perturbation is mass conserving
       double value = sin(2.0*MathematicalConstants::Pi*x/3.0)*
        epsilon*time*exp(-time);
       nod_pt->set_value(1,value);
      }

            
       //If we are on the bottom boundary, set the temperature
       //to 2 (heated)
     if(ibound==0) {nod_pt->set_value(2,2.0);}
    }
  }
} // end_of_set_boundary_conditions

//===============start_doc_solution=======================================
/// Doc the solution
//========================================================================
template<class ELEMENT,class INTERFACE_ELEMENT>
void ConvectionProblem<ELEMENT,INTERFACE_ELEMENT>::doc_solution(
 ofstream &trace)
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   Bulk_mesh_pt->finite_element_pt(e)->output(some_file,npts);
  }
 some_file.close();

 //Output the interface
 sprintf(filename,"%s/int%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);

 unsigned n_interface  = Surface_mesh_pt->nelement();
 for(unsigned i=0;i<n_interface;i++)
  {
   Surface_mesh_pt->finite_element_pt(i)->output(some_file,npts);
  }
 some_file.close();

 trace << time_pt()->time() << " " 
       << Bulk_mesh_pt->boundary_node_pt(2,8)->value(2) << " "
       << Bulk_mesh_pt->boundary_node_pt(2,8)->value(3) << " "
       << Bulk_mesh_pt->boundary_node_pt(2,8)->value(4) << std::endl;


 Doc_info.number()++;
} // end of doc


//=======start_of_main================================================
/// Driver code for 2D Boussinesq convection problem
//====================================================================
int main(int argc, char **argv)
{
 ofstream trace("RESLT/trace.dat");

 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;

 Global_Physical_Variables::D[1] = 0.01;

 //Construct our problem
 ConvectionProblem<SpineElement<DoubleBuoyantQCrouzeixRaviartElement<2> >,
  SpineLineMarangoniSurfactantFluidInterfaceElement<DoubleBuoyantQCrouzeixRaviartElement<2> > > 
problem;

 // Apply the boundary condition at time zero
 problem.set_boundary_conditions(0.0);
 
 //Perform a single steady Newton solve
 problem.steady_newton_solve();
 
 //Document the solution
 problem.doc_solution(trace);

 //Now release the interface for real fun
 //problem.unpin_surface();

 //Set the timestep
 double dt = 0.1;

 //Initialise the value of the timestep and set initial values 
 //of previous time levels assuming an impulsive start.
 problem.assign_initial_values_impulsive(dt);

 //Set the number of timesteps to our default value
 unsigned n_steps = 4000;

 //If we have a command line argument, perform fewer steps 
 //(used for self-test runs)
 if(argc > 1) {n_steps = 5;}

 //Perform n_steps timesteps
 for(unsigned i=0;i<n_steps;++i)
  {
   problem.unsteady_newton_solve(dt);
   problem.doc_solution(trace);
  }

} // end of main









