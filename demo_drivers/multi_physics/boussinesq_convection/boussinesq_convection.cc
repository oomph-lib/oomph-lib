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
//Driver for a multi-physics problem that couples the Navier--Stokes
//equations to the advection diffusion equations, 
//giving Boussinesq convection

//Oomph-lib headers, we require the generic, advection-diffusion
//and navier-stokes elements.
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"

// The mesh is our standard rectangular quadmesh
#include "meshes/rectangular_quadmesh.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;


//======================class definition==============================
///A class that solves the Boussinesq approximation of the Navier--Stokes
///and energy equations by coupling two pre-existing classes. 
///The QAdvectionDiffusionElement with bi-quadratic interpolation for the
///scalar variable (temperature) and
///QCrouzeixRaviartElement which solves the Navier--Stokes equations
///using bi-quadratic interpolation for the velocities and a discontinuous
///bi-linear interpolation for the pressure. Note that we are free to 
///choose the order in which we store the variables at the nodes. In this
///case we choose to store the variables in the order fluid velocities
///followed by temperature. We must, therefore, overload the function
///AdvectionDiffusionEquations<DIM>::u_index_adv_diff() to indicate that
///the temperature is stored at the DIM-th position not the 0-th. We do not
///need to overload the corresponding function in the 
///NavierStokesEquations<DIM> class because the velocities are stored
///first.
//=========================================================================
template<unsigned DIM>
class BuoyantQCrouzeixRaviartElement :
 public virtual QAdvectionDiffusionElement<DIM,3>,
 public virtual QCrouzeixRaviartElement<DIM>
{

private:

 /// Pointer to a private data member, the Rayleigh number
 double* Ra_pt;

 /// The static default value of the Rayleigh number
 static double Default_Physical_Constant_Value;

public:

 /// \short Constructor: call the underlying constructors and 
 /// initialise the pointer to the Rayleigh number to point
 /// to the default value of 0.0.
 BuoyantQCrouzeixRaviartElement() : QAdvectionDiffusionElement<DIM,3>(),
                                   QCrouzeixRaviartElement<DIM>() 
  {
   Ra_pt = &Default_Physical_Constant_Value;
  }

 ///\short The required number of values stored at the nodes is the sum of the
 ///required values of the two single-physics  elements. Note that this step is
 ///generic for any multi-physics element of this type.
 unsigned required_nvalue(const unsigned &n) const
  {return (QAdvectionDiffusionElement<DIM,3>::required_nvalue(n) +
           QCrouzeixRaviartElement<DIM>::required_nvalue(n));}

 ///Access function for the Rayleigh number (const version)
 const double &ra() const {return *Ra_pt;}

 ///Access function for the pointer to the Rayleigh number
 double* &ra_pt() {return Ra_pt;}

 /// Final override for disable ALE
 void disable_ALE() 
  {
   //Disable ALE in both sets of equations
   NavierStokesEquations<DIM>::disable_ALE();
   AdvectionDiffusionEquations<DIM>::disable_ALE();
  }

 /// Final override for enable ALE
 void enable_ALE() 
  {
   //Enable ALE in both sets of equations
   NavierStokesEquations<DIM>::enable_ALE();
   AdvectionDiffusionEquations<DIM>::enable_ALE();
  }


 ///  Overload the standard output function with the broken default
 void output(ostream &outfile) {FiniteElement::output(outfile);}

 /// \short Output function:  
 ///  Output x, y, u, v, p, theta at Nplot^DIM plot points
 // Start of output function
 void output(ostream &outfile, const unsigned &nplot)
  {
   //vector of local coordinates
   Vector<double> s(DIM);
   
   // Tecplot header info
   outfile << this->tecplot_zone_string(nplot);
   
   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,nplot,s);
     
     // Output the position of the plot point
     for(unsigned i=0;i<DIM;i++) 
      {outfile << this->interpolated_x(s,i) << " ";}
     
     // Output the fluid velocities at the plot point
     for(unsigned i=0;i<DIM;i++) 
      {outfile << this->interpolated_u_nst(s,i) << " ";}
     
     // Output the fluid pressure at the plot point
     outfile << this->interpolated_p_nst(s)  << " ";
  
     // Output the temperature (the advected variable) at the plot point
     outfile << this->interpolated_u_adv_diff(s) << std::endl;   
    }
   outfile << std::endl;
   
   // Write tecplot footer (e.g. FE connectivity lists)
   this->write_tecplot_zone_footer(outfile,nplot);
  } //End of output function


 /// \short C-style output function: Broken default
 void output(FILE* file_pt)
  {FiniteElement::output(file_pt);}

 ///  \short C-style output function: Broken default
 void output(FILE* file_pt, const unsigned &n_plot)
  {FiniteElement::output(file_pt,n_plot);}

 /// \short Output function for an exact solution: Broken default
 void output_fct(ostream &outfile, const unsigned &Nplot,
                 FiniteElement::SteadyExactSolutionFctPt 
                 exact_soln_pt)
  {FiniteElement::output_fct(outfile,Nplot,exact_soln_pt);}


 /// \short Output function for a time-dependent exact solution:
 /// Broken default.
 void output_fct(ostream &outfile, const unsigned &Nplot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt 
                 exact_soln_pt)
  {
   FiniteElement::
    output_fct(outfile,Nplot,time,exact_soln_pt);
  }

 ///\short Overload the index at which the temperature 
 ///variable is stored. We choose to store is after the fluid velocities.
 inline unsigned u_index_adv_diff() const {return DIM;}
  
 /// \short Validate against exact solution at given time
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points and compute L2 error
 /// and L2 norm of velocity solution over element
 /// Call the broken default
 void compute_error(ostream &outfile,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time,
                    double& error, double& norm)
  {FiniteElement::compute_error(outfile,exact_soln_pt,
                                time,error,norm);}
 
 /// \short Validate against exact solution.
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points and compute L2 error
 /// and L2 norm of velocity solution over element
 /// Call the broken default
 void compute_error(ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm)
  {FiniteElement::compute_error(outfile,exact_soln_pt,error,norm);}

 /// \short Overload the wind function in the advection-diffusion equations.
 /// This provides the coupling from the Navier--Stokes equations to the
 /// advection-diffusion equations because the wind is the fluid velocity.
 void get_wind_adv_diff(const unsigned& ipt,
                        const Vector<double> &s, const Vector<double>& x,
                        Vector<double>& wind) const
 {
   //The wind function is simply the velocity at the points
   this->interpolated_u_nst(s,wind);
 }


 /// \short Overload the body force in the Navier-Stokes equations
 /// This provides the coupling from the advection-diffusion equations
 /// to the Navier--Stokes equations, the body force is the
 /// temperature multiplied by the Rayleigh number acting in the
 /// direction opposite to gravity. 
 void get_body_force_nst(const double& time, const unsigned& ipt,
                         const Vector<double> &s, const Vector<double> &x,
                         Vector<double> &result)
  {
   //Get the temperature
   const double interpolated_t = this->interpolated_u_adv_diff(s);

   // Get vector that indicates the direction of gravity from
   // the Navier-Stokes equations
   Vector<double> gravity(NavierStokesEquations<DIM>::g());
   
   // Temperature-dependent body force:
   for (unsigned i=0;i<DIM;i++)
    {
     result[i] = -gravity[i]*interpolated_t*ra();
    }
  }

 /// \short Calculate the element's contribution to the residual vector.
 /// Recall that fill_in_* functions MUST NOT initialise the entries 
 /// in the vector to zero. This allows us to call the 
 /// fill_in_* functions of the constituent single-physics elements
 /// sequentially, without wiping out any previously computed entries.
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Fill in the residuals of the Navier-Stokes equations
   NavierStokesEquations<DIM>::fill_in_contribution_to_residuals(residuals);

   //Fill in the residuals of the advection-diffusion eqautions
   AdvectionDiffusionEquations<DIM>::fill_in_contribution_to_residuals(
    residuals);
  }


//-----------Finite-difference the entire jacobian-----------------------
//-----------------------------------------------------------------------
#ifdef USE_FD_JACOBIAN_FOR_BUOYANT_Q_ELEMENT


 ///\short Compute the element's residual vector and the Jacobian matrix.
 /// Jacobian is computed by finite-differencing.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   // This function computes the Jacobian by finite-differencing
   FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
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

#else
//--------Finite-difference off-diagonal-blocks in jacobian--------------
//-----------------------------------------------------------------------
#ifdef USE_OFF_DIAGONAL_FD_JACOBIAN_FOR_BUOYANT_Q_ELEMENT

 ///\short Helper function to get the off-diagonal blocks of the Jacobian
 ///matrix by finite differences
 void fill_in_off_diagonal_jacobian_blocks_by_fd(Vector<double> &residuals,
                                                 DenseMatrix<double> &jacobian)
  {
   //Local storage for the index in the nodes at which the
   //Navier-Stokes velocities are stored (we know that this should be 0,1,2)
   unsigned u_nodal_nst[DIM];
   for(unsigned i=0;i<DIM;i++) 
    {u_nodal_nst[i] = this->u_index_nst(i);}

   //Local storage for the  index at which the temperature is stored
   unsigned u_nodal_adv_diff = this->u_index_adv_diff();

   //Find the total number of unknowns in the elements
   unsigned n_dof = this->ndof();

   //Temporary storage for residuals
   Vector<double> newres(n_dof);

   //Storage for local unknown and local equation
   int local_unknown =0, local_eqn = 0;
   
   //Set the finite difference step
   double fd_step = FiniteElement::Default_fd_jacobian_step;

   //Find the number of nodes
   unsigned n_node = this->nnode();
   
   //Calculate the contribution of the Navier--Stokes velocities to the
   //advection-diffusion equations
   
   //Loop over the nodes
   for(unsigned n=0;n<n_node;n++)
    {
     //There are DIM values of the  velocities
     for(unsigned i=0;i<DIM;i++)
      {
       //Get the local velocity equation number
       local_unknown = nodal_local_eqn(n,u_nodal_nst[i]);

       //If it's not pinned
       if(local_unknown >= 0)
        {
         //Get a pointer to the velocity value
         double *value_pt = this->node_pt(n)->value_pt(u_nodal_nst[i]); 

         //Save the old value
         double old_var = *value_pt;

         //Increment the value 
         *value_pt += fd_step;

         //Get the altered advection-diffusion residuals.
         AdvectionDiffusionEquations<DIM>::get_residuals(newres);

         //Now fill in the Advection-Diffusion sub-block
         //of the jacobian
         for(unsigned m=0;m<n_node;m++)
          {
           //Local equation for temperature
           local_eqn = this->nodal_local_eqn(m,u_nodal_adv_diff);

           //If it's not a boundary condition
           if(local_eqn >= 0)
           {
            double sum = (newres[local_eqn] - residuals[local_eqn])/fd_step;
            jacobian(local_eqn,local_unknown) = sum;
           }
          }
         
         //Reset the nodal data
         *value_pt = old_var;
        }
      }


     //Calculate the contribution of the temperature to the Navier--Stokes
     //equations
     {
      //Get the local equation number
      local_unknown = this->nodal_local_eqn(n,u_nodal_adv_diff);
      
      //If it's not pinned
      if(local_unknown >= 0)
       {
        //Get a pointer to the concentration value
        double *value_pt = this->node_pt(n)->value_pt(u_nodal_adv_diff); 

        //Save the old value
        double old_var = *value_pt;

        //Increment the value (Really need access)
        *value_pt += fd_step;
        
        //Get the altered Navier--Stokes residuals
        NavierStokesEquations<DIM>::get_residuals(newres);
         
        //Now fill in the Navier-Stokes sub-block
        for(unsigned m=0;m<n_node;m++)
         {
          //Loop over the fluid velocities
          for(unsigned j=0;j<DIM;j++)
           {
            //Local fluid equation
            local_eqn = this->nodal_local_eqn(m,u_nodal_nst[j]);
            if(local_eqn >= 0)
             {
              double sum = (newres[local_eqn] - residuals[local_eqn])/fd_step;
              jacobian(local_eqn,local_unknown) = sum;
             }
           }
         }
        
        //Reset the nodal data
        *value_pt = old_var;
       }
     }
     
    } //End of loop over nodes
  }


 ///\short Compute the element's residual Vector and the Jacobian matrix.
 /// Use finite-differencing only for the off-diagonal blocks.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {

   //Calculate the Navier-Stokes contributions (diagonal block and residuals)
   NavierStokesEquations<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //Calculate the advection-diffusion contributions 
   //(diagonal block and residuals)
   AdvectionDiffusionEquations<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //We now fill in the off-diagonal (interaction) blocks by finite
   //differences.
   fill_in_off_diagonal_jacobian_blocks_by_fd(residuals,jacobian);
  } //End of jacobian calculation


 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   //Get the analytic diagonal terms
   NavierStokesEquations<DIM>::
    fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);
   
   AdvectionDiffusionEquations<DIM>::
    fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);

   //Now fill in the off-diagonal blocks
   fill_in_off_diagonal_jacobian_blocks_by_fd(residuals,jacobian);
  }


 //--------------------Analytic jacobian---------------------------------
//-----------------------------------------------------------------------
#else

 ///\short Helper function to get the off-diagonal blocks of the Jacobian
 ///matrix analytically
 void fill_in_off_diagonal_jacobian_blocks_analytic(
  Vector<double> &residuals, DenseMatrix<double> &jacobian)
  {
   //We now fill in the off-diagonal (interaction) blocks analytically
   //This requires knowledge of exactly how the residuals are assembled
   //within the parent elements and involves yet another loop over
   //the integration points!
   
   //Local storage for the index in the nodes at which the
   //Navier-Stokes velocities are stored (we know that this should be 0,1,2)
   unsigned u_nodal_nst[DIM];
   for(unsigned i=0;i<DIM;i++) {u_nodal_nst[i] = this->u_index_nst(i);}
   
   //Local storage for the  index at which the temperature is stored
   const unsigned u_nodal_adv_diff = this->u_index_adv_diff();

   //Find out how many nodes there are
   const unsigned n_node = this->nnode();
   
   //Set up memory for the shape and test functions and their derivatives
   Shape psif(n_node), testf(n_node);
   DShape dpsifdx(n_node,DIM), dtestfdx(n_node,DIM);
   
   //Number of integration points
   const unsigned n_intpt = this->integral_pt()->nweight();
   
   //Get Physical Variables from Element
   double Ra = this->ra();
   double Pe = this->pe();
   Vector<double> gravity = this->g();
   
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
     
     //Calculate local values of temperature derivatives
     //Allocate
     Vector<double> interpolated_du_adv_diff_dx(DIM,0.0);
     
     // Loop over nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       //Get the nodal value
       double u_value = this->raw_nodal_value(l,u_nodal_adv_diff);
       //Loop over the derivative directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_du_adv_diff_dx[j] += u_value*dpsifdx(l,j);
        }
      }
     
     //Assemble the jacobian terms
     
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
           //Loop over the velocity shape functions again
           for(unsigned l2=0;l2<n_node;l2++)
            { 
             //We have only the temperature degree of freedom to consider
             //If it's non-zero add in the contribution
             local_unknown = this->nodal_local_eqn(l2,u_nodal_adv_diff);
             if(local_unknown >= 0)
              {
               //Add contribution to jacobian matrix
               jacobian(local_eqn,local_unknown) 
                += -gravity[i]*psif(l2)*Ra*testf(l)*W;
              }
            }
          }
        }
       
       //Assemble the contributions of the fluid velocities to the
       //advection-diffusion equation for the temperature
       {
        local_eqn = this->nodal_local_eqn(l,u_nodal_adv_diff);
        //If it's not pinned
        if(local_eqn >= 0)
         {
          //Loop over the shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            //Loop over the velocity degrees of freedom
            for(unsigned i2=0;i2<DIM;i2++)
             {
              //Get the local unknown
              local_unknown = this->nodal_local_eqn(l2,u_nodal_nst[i2]);
              //If it's not pinned
              if(local_unknown >= 0)
               {
                //Add the contribution to the jacobian matrix
                jacobian(local_eqn,local_unknown)
                 -= Pe*psif(l2)*interpolated_du_adv_diff_dx[i2]*testf(l)*W;
               }
             }
           }
         }
       }
       
      }
    }
  }


 ///\short Compute the element's residual Vector and the Jacobian matrix.
 /// Use analytic expressions for the off-diagonal blocks
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
  {
   
   //Calculate the Navier-Stokes contributions (diagonal block and residuals)
   NavierStokesEquations<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //Calculate the advection-diffusion contributions 
   //(diagonal block and residuals)
   AdvectionDiffusionEquations<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //Fill in the off diagonal blocks analytically
   fill_in_off_diagonal_jacobian_blocks_analytic(residuals,jacobian);
  }


 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   //Get the analytic diagonal terms for the jacobian and mass matrix
   NavierStokesEquations<DIM>::
    fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);
   
   AdvectionDiffusionEquations<DIM>::
    fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);

   //Now fill in the off-diagonal blocks in the jacobian matrix
   fill_in_off_diagonal_jacobian_blocks_analytic(residuals,jacobian);
  }
 
#endif
#endif

};

//=========================================================================
/// Set the default physical value to be zero
//=========================================================================
template<>
double BuoyantQCrouzeixRaviartElement<2>::Default_Physical_Constant_Value=0.0;


//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{
 /// Peclet number (identically one from our non-dimensionalisation)
 double Peclet=1.0;

 /// 1/Prandtl number
 double Inverse_Prandtl=1.0;

 /// \short Rayleigh number, set to be greater than 
 /// the threshold for linear instability
 double Rayleigh = 1800.0;

 /// Gravity vector
 Vector<double> Direction_of_gravity(2);
  
} // end_of_namespace

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//====== start_of_problem_class=======================================
/// 2D Convection  problem on rectangular domain, discretised 
/// with refineable elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class ConvectionProblem : public Problem
{

public:

 ///Constructor
 ConvectionProblem();

 /// Destructor. Empty
 ~ConvectionProblem() {}

 /// \short Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Actions before adapt:(empty)
 void actions_before_adapt(){}

 /// \short Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep() 
  {
   set_boundary_conditions(time_pt()->time());
  }

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

 /// \short Doc the solution.
 void doc_solution();

 /// \short Set the boundary conditions
 void set_boundary_conditions(const double &time);

 /// \short Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }
 
private:
 
 /// DocInfo object
 DocInfo Doc_info;

}; // end of problem class

//===========start_of_constructor=========================================
/// \short Constructor for convection problem
//========================================================================
template<class ELEMENT>
ConvectionProblem<ELEMENT>::ConvectionProblem()
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
 Problem::mesh_pt() = 
  new RectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here

 //Loop over the boundaries
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the maximum index to be pinned (all values by default)
   unsigned val_max=3;
   //If we are on the side-walls, the v-velocity and temperature
   //satisfy natural boundary conditions, so we only pin the
   //first value
   if((ibound==1) || (ibound==3)) {val_max=1;}

   //Loop over the number of nodes on the boundry
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //Loop over the desired values stored at the nodes and pin
     for(unsigned j=0;j<val_max;j++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(j);
      }
    }
  }

 //Pin the zero-th pressure dof in element 0 and set its value to
 //zero:
 fix_pressure(0,0,0.0);

 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor.
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;

   // Set the Peclet number multiplied by the Strouhal number
   el_pt->pe_st_pt() =&Global_Physical_Variables::Peclet;

   // Set the Reynolds number (1/Pr in our non-dimensionalisation)
   el_pt->re_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set the Rayleigh number
   el_pt->ra_pt() = &Global_Physical_Variables::Rayleigh;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;

   //The mesh is fixed, so we can disable ALE
   el_pt->disable_ALE();

   // Set pointer to the continuous time
   el_pt->time_pt() = time_pt();
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor



//===========start_of_set_boundary_conditions================
/// Set the boundary conditions as a function of continuous 
/// time
//===========================================================
template<class ELEMENT>
void ConvectionProblem<ELEMENT>::set_boundary_conditions(
 const double &time)
{
 // Loop over the boundaries
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

     //Set the number of velocity components
     unsigned vel_max=2;

     //If we are on the side walls we only set the x-velocity.
     if((ibound==1) || (ibound==3)) {vel_max = 1;}

     //Set the pinned velocities to zero
     for(unsigned j=0;j<vel_max;j++) {nod_pt->set_value(j,0.0);}

     //If we are on the top boundary
     if(ibound==2) 
      {
       //Set the temperature to -0.5 (cooled)
       nod_pt->set_value(2,-0.5);

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
     //to 0.5 (heated)
     if(ibound==0) {nod_pt->set_value(2,0.5);}
    }
  }
} // end_of_set_boundary_conditions

//===============start_doc_solution=======================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ConvectionProblem<ELEMENT>::doc_solution()
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
 ConvectionProblem<BuoyantQCrouzeixRaviartElement<2> > problem;

 // Apply the boundary condition at time zero
 problem.set_boundary_conditions(0.0);

 //Perform a single steady Newton solve
 problem.steady_newton_solve();

 //Document the solution
 problem.doc_solution();

 //Set the timestep
 double dt = 0.1;

 //Initialise the value of the timestep and set an impulsive start
 problem.assign_initial_values_impulsive(dt);

 //Set the number of timesteps to our default value
 unsigned n_steps = 200;

 //If we have a command line argument, perform fewer steps 
 //(used for self-test runs)
 if(argc > 1) {n_steps = 5;}

 //Perform n_steps timesteps
 for(unsigned i=0;i<n_steps;++i)
  {
   problem.unsteady_newton_solve(dt);
   problem.doc_solution();
  }

} // end of main









