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
#ifndef BUOYANT_NAVIER_STOKES_ELEMENTS_HEADER
#define BUOYANT_NAVIER_STOKES_ELEMENTS_HEADER

//Oomph-lib headers, we require the generic, advection-diffusion-reaction,
//navier-stokes elements and fluid interface elements
#include "generic.h"
#include "advection_diffusion_reaction.h"
#include "navier_stokes.h"

namespace oomph
{

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
 public virtual QAdvectionDiffusionReactionElement<1,DIM,3>,
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
 BuoyantQCrouzeixRaviartElement() : 
  QAdvectionDiffusionReactionElement<1,DIM,3>(),
  QCrouzeixRaviartElement<DIM>() 
  {
   Ra_pt = &Default_Physical_Constant_Value;
  }
 
 /// UnPin p_dof-th pressure dof 
 void unfix_pressure(const unsigned &p_dof)
  {
   this->internal_data_pt(this->P_nst_internal_index)->unpin(p_dof);
  }


 ///\short The required number of values stored at the nodes is the sum of the
 ///required values of the two single-physics  elements. Note that this step is
 ///generic for any multi-physics element of this type.
 unsigned required_nvalue(const unsigned &n) const
  {return (QAdvectionDiffusionReactionElement<1,DIM,3>::required_nvalue(n) +
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
   AdvectionDiffusionReactionEquations<1,DIM>::disable_ALE();
  }

 /// Final override for enable ALE
 void enable_ALE() 
  {
   //Enable ALE in both sets of equations
   NavierStokesEquations<DIM>::enable_ALE();
   AdvectionDiffusionReactionEquations<1,DIM>::enable_ALE();
  }


 ///  Overload the standard output function with the broken default
 void output(std::ostream &outfile) {FiniteElement::output(outfile);}

 /// \short Output function:  
 ///  Output x, y, u, v, p, theta at Nplot^DIM plot points
 // Start of output function
 void output(std::ostream &outfile, const unsigned &nplot)
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
     outfile << this->interpolated_c_adv_diff_react(s,0) << std::endl;   
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
 void output_fct(std::ostream &outfile, const unsigned &Nplot,
                 FiniteElement::SteadyExactSolutionFctPt 
                 exact_soln_pt)
  {FiniteElement::output_fct(outfile,Nplot,exact_soln_pt);}


 /// \short Output function for a time-dependent exact solution:
 /// Broken default.
 void output_fct(std::ostream &outfile, const unsigned &Nplot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt 
                 exact_soln_pt)
  {
   FiniteElement::
    output_fct(outfile,Nplot,time,exact_soln_pt);
  }

 ///\short Overload the index at which the temperature 
 ///variable is stored. We choose to store is after the fluid velocities.
 inline unsigned c_index_adv_diff_react(const unsigned &i) const {return DIM;}
  
 /// \short Validate against exact solution at given time
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points and compute L2 error
 /// and L2 norm of velocity solution over element
 /// Call the broken default
 void compute_error(std::ostream &outfile,
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
 void compute_error(std::ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm)
  {FiniteElement::compute_error(outfile,exact_soln_pt,error,norm);}

 /// \short Overload the wind function in the advection-diffusion equations.
 /// This provides the coupling from the Navier--Stokes equations to the
 /// advection-diffusion equations because the wind is the fluid velocity.
 void get_wind_adv_diff_react(const Vector<double> &s, const Vector<double>& x,
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
 void get_body_force_nst(const double& time, const Vector<double> &s,
                         const Vector<double> &x, 
                         Vector<double> &result)
  {
   // Get vector that indicates the direction of gravity from
   // the Navier-Stokes equations
   Vector<double> gravity(NavierStokesEquations<DIM>::g());
   
   // Temperature-dependent body force:
   for (unsigned i=0;i<DIM;i++)
    {
     result[i] = -gravity[i]*this->interpolated_c_adv_diff_react(s,0)*ra();
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
   AdvectionDiffusionReactionEquations<1,DIM>::
    fill_in_contribution_to_residuals(residuals);
  }


#ifdef USE_FD_JACOBIAN_FOR_BUOYANT_Q_CROZIER_RAVIART_ELEMENT


 ///\short Compute the element's residual vector and the Jacobian matrix.
 /// Jacobian is computed by finite-differencing.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   // This function computes the Jacobian by finite-differencing
   FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
  }

#else

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
   unsigned T_nodal_adv_diff_react = this->c_index_adv_diff_react(0);

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
         AdvectionDiffusionReactionEquations<1,DIM>::get_residuals(newres);

         //Now fill in the Advection-DiffusionReaction sub-block
         //of the jacobian
         for(unsigned m=0;m<n_node;m++)
          {
           //Local equation for temperature
           local_eqn = this->nodal_local_eqn(m,T_nodal_adv_diff_react);

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
      local_unknown = this->nodal_local_eqn(n,T_nodal_adv_diff_react);

      //If it's not pinned
      if(local_unknown >= 0)
       {
        //Get a pointer to the concentration value
        double *value_pt = this->node_pt(n)->value_pt(T_nodal_adv_diff_react); 

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
   AdvectionDiffusionReactionEquations<1,DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //Add in the off-diagonal blocks
   fill_in_off_diagonal_jacobian_blocks_by_fd(residuals,jacobian);
  } //End of jacobian calculation

#endif


 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   NavierStokesEquations<DIM>::
    fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);
   
   AdvectionDiffusionReactionEquations<1,DIM>::
    fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);

   fill_in_off_diagonal_jacobian_blocks_by_fd(residuals,jacobian);
  }


};

//=========================================================================
/// Set the default physical value to be zero
//=========================================================================
template<>
double BuoyantQCrouzeixRaviartElement<2>::Default_Physical_Constant_Value=0.0;


//=======================================================================
/// Face geometry of the 2D Buoyant Crouzeix_Raviart elements
//=======================================================================
template<unsigned int DIM>
class FaceGeometry<BuoyantQCrouzeixRaviartElement<DIM> >: 
public virtual QElement<DIM-1,3>
{
  public:
 FaceGeometry() : QElement<DIM-1,3>() {}
};

//=======================================================================
/// Face geometry of the 2D Buoyant Crouzeix_Raviart elements
//=======================================================================
template<>
class FaceGeometry<FaceGeometry<BuoyantQCrouzeixRaviartElement<2> > >: 
public virtual PointElement
{
  public:
 FaceGeometry() : PointElement() {}
};

} //End of the oomph namespace

#endif
