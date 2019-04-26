//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
#ifndef DOUBLE_BUOYANT_NAVIER_STOKES_ELEMENTS_HEADER
#define DOUBLE_BUOYANT_NAVIER_STOKES_ELEMENTS_HEADER

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
///The QAdvectionDiffusionReactionElement with 
///bi-quadratic interpolation for the
///scalar variables (temperature and concentration) and
///QCrouzeixRaviartElement which solves the Navier--Stokes equations
///using bi-quadratic interpolation for the velocities and a discontinuous
///bi-linear interpolation for the pressure. Note that we are free to 
///choose the order in which we store the variables at the nodes. In this
///case we choose to store the variables in the order fluid velocities
///followed by bulk concentration followed by micelle concentration.
///We must, therefore, overload the function
///AdvectionDiffusionReactionEquations<2,DIM>::u_index_adv_diff_react() 
///to indicate that
///the concentraion is stored at the DIM-th position not the 0-th. We do not
///need to overload the corresponding function in the 
///NavierStokesEquations<DIM> class because the velocities are stored
///first.
//=========================================================================
template<unsigned DIM>
class DoubleBuoyantQCrouzeixRaviartElement :
 public virtual QAdvectionDiffusionReactionElement<2,DIM,3>,
 public virtual QCrouzeixRaviartElement<DIM>
{
  
private:

  ///Pointer to private data. The value of Km
  double *Km_pt;

  //Pointer to private data. The value of N
  double *N_pt;
  
 /// The static default value of the Rayleigh number
 static double Default_Physical_Constant_Value;

public:

 /// \short Constructor: call the underlying constructors and 
 /// initialise the pointer to the Rayleigh number to point
 /// to the default value of 0.0.
 DoubleBuoyantQCrouzeixRaviartElement() : 
  QAdvectionDiffusionReactionElement<2,DIM,3>(),
  QCrouzeixRaviartElement<DIM>() 
  {
    Km_pt = &Default_Physical_Constant_Value;
    N_pt = &Default_Physical_Constant_Value;
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
  {return (QAdvectionDiffusionReactionElement<2,DIM,3>::required_nvalue(n) +
           QCrouzeixRaviartElement<DIM>::required_nvalue(n));}


 ///Access function for the transfer constant
 const double &km() const {return *Km_pt;}

 ///Access function for the pointer to transfer constant
 double* &km_pt() {return Km_pt;}

 ///Access function for the number of monomers in the micelle
 const double &n() const {return *N_pt;}

 ///Access function for the pointer to the number of monomers in the micelle
 double* &n_pt() {return N_pt;}
  
 /// Final override for disable ALE
 void disable_ALE() 
  {
   //Disable ALE in both sets of equations
   NavierStokesEquations<DIM>::disable_ALE();
   AdvectionDiffusionReactionEquations<2,DIM>::disable_ALE();
  }

 /// Final override for enable ALE
 void enable_ALE() 
  {
   //Enable ALE in both sets of equations
   NavierStokesEquations<DIM>::enable_ALE();
   AdvectionDiffusionReactionEquations<2,DIM>::enable_ALE();
  }

  ///Overload the reaction terms to couple the concentration and
  ///micelle terms
  void get_reaction_adv_diff_react(const unsigned &ipt,
				   const Vector<double> &C,
				   Vector<double> &R) const
  {
    //Compute the flux between equations (equation (2.22))
    const double J_m = this->km()*(pow(C[0],this->n()) - C[1]);
    //Return the reaction terms
    R[0] = J_m;
    R[1] = -J_m;
  }

  //Fill in the derivatives of the reaction terms
  void get_reaction_deriv_adv_diff_react(const unsigned& ipt,
					 const Vector<double> &C,
					 DenseMatrix<double> &dRdC)
  const
  {
    const double Km = this->km(); const double N = this->n();
    //Fill in the derivative terms
    dRdC(0,0) = Km*N*pow(C[0],N-1); dRdC(0,1) = -Km;
    dRdC(1,0) = -dRdC(0,0); dRdC(1,1) = -dRdC(0,1);
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

     //Output the temperature and the solute concentration
     for(unsigned i=0;i<2;i++)
      {
       outfile << this->interpolated_c_adv_diff_react(s,i) << " ";
      }
     outfile << "\n";
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

 ///\short Overload the index at which the temperature and solute
 ///concentration variables are stored. 
 // We choose to store them after the fluid velocities.
 inline unsigned c_index_adv_diff_react(const unsigned &i) const 
  {return DIM+i;}
  
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
 void get_wind_adv_diff_react(const unsigned& ipt,
                              const Vector<double> &s, const Vector<double>& x,
                              Vector<double>& wind) const
 {
   //The wind function is simply the velocity at the points
   this->interpolated_u_nst(s,wind);
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
   AdvectionDiffusionReactionEquations<2,DIM>::
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

   //Local storage for the  index at which the temperature and
   //solute are  stored
   unsigned C_nodal_adv_diff_react[2];
   for(unsigned r=0;r<2;r++)
    {C_nodal_adv_diff_react[r] = this->c_index_adv_diff_react(r);}
   

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
       local_unknown = this->nodal_local_eqn(n,u_nodal_nst[i]);

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
         //Do this using fill_in because get_residuals has never been 
         //overloaded, and will actually compute all residuals which
         //is slightly inefficient.
         for(unsigned m=0;m<n_dof;m++) {newres[m] = 0.0;}
         AdvectionDiffusionReactionEquations<2,DIM>::
          fill_in_contribution_to_residuals(newres);

         //Now fill in the Advection-Diffusion-Reaction sub-block
         //of the jacobian
         for(unsigned m=0;m<n_node;m++)
          {
           for(unsigned r=0;r<2;r++)
            {
             //Local equation for temperature or solute
             local_eqn = this->nodal_local_eqn(m,C_nodal_adv_diff_react[r]);

             //If it's not a boundary condition
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


     //Calculate the contribution of the temperature to the Navier--Stokes
     //equations
     for(unsigned r=0;r<2;r++)
      {
       //Get the local equation number
       local_unknown = this->nodal_local_eqn(n,C_nodal_adv_diff_react[r]);
       
       //If it's not pinned
       if(local_unknown >= 0)
        {
         //Get a pointer to the concentration value
         double *value_pt = 
          this->node_pt(n)->value_pt(C_nodal_adv_diff_react[r]); 
         
         //Save the old value
         double old_var = *value_pt;
         
         //Increment the value (Really need access)
         *value_pt += fd_step;
         
         //Get the altered Navier--Stokes residuals
         //Do this using fill_in because get_residuals has never been 
         //overloaded, and will actually compute all residuals which
         //is slightly inefficient.
         for(unsigned m=0;m<n_dof;m++) {newres[m] = 0.0;}
         NavierStokesEquations<DIM>::fill_in_contribution_to_residuals(newres);
         
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
      } //End of loop over reagents
     
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
   AdvectionDiffusionReactionEquations<2,DIM>::
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
   
   AdvectionDiffusionReactionEquations<2,DIM>::
    fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);

   fill_in_off_diagonal_jacobian_blocks_by_fd(residuals,jacobian);
  }

 void integrated_C_and_M(double &int_C, double &int_M)
 {
  //Vector to store the integrals
  Vector<double> sum(2,0.0);

  //Find the number of nodes
  const unsigned n_node = this->nnode();
  //Storage for the shape functions
  Shape psi(n_node);

  unsigned C_index[2];
  for(unsigned r=0;r<2;++r)
   {C_index[r] = this->c_index_adv_diff_react(r);}
  
  //Loop over the integration points
  const unsigned n_intpt = this->integral_pt()->nweight();
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Get the integral weight
    double w = this->integral_pt()->weight(ipt);
    //Get the value of the Jacobian of the mapping to global coordinates
    double J = this->J_eulerian_at_knot(ipt);
    double W = w*J;
    //Get the shape function at the know
    this->shape_at_knot(ipt,psi);

    //Get the interpolated values
    Vector<double> interpolated_C(2,0.0);
    for(unsigned l=0;l<n_node;++l)
     {
      const double psi_ = psi(l);
      for(unsigned r=0;r<2;++r)
       {
        interpolated_C[r] += this->nodal_value(l,C_index[r])*psi_;
       }
     }
          
    for(unsigned r=0;r<2;++r)
     {
      sum[r] += interpolated_C[r]*W;
     }
   } //End of integration loop

  //Return the values
  int_C = sum[0]; int_M = sum[1];
 }

};

//=========================================================================
/// Set the default physical value to be one
//=========================================================================
template<>
double DoubleBuoyantQCrouzeixRaviartElement<2>::Default_Physical_Constant_Value=1.0;


//=======================================================================
/// Face geometry of the 2D DoubleBuoyant Crouzeix_Raviart elements
//=======================================================================
template<unsigned int DIM>
class FaceGeometry<DoubleBuoyantQCrouzeixRaviartElement<DIM> >: 
public virtual QElement<DIM-1,3>
{
  public:
 FaceGeometry() : QElement<DIM-1,3>() {}
};

//=======================================================================
/// Face geometry of the 2D DoubleBuoyant Crouzeix_Raviart elements
//=======================================================================
template<>
class FaceGeometry<FaceGeometry<DoubleBuoyantQCrouzeixRaviartElement<2> > >: 
public virtual PointElement
{
  public:
 FaceGeometry() : PointElement() {}
};

} //End of the oomph namespace

#endif
