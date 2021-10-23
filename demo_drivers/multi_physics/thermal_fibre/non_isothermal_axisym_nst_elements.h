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
#define USE_FD_JACOBIAN_FOR_NONISOTHERMAL_AXISYMMETRIC_NAVIER_STOKES

//======================class definition======================
///A class that solves the nonisothermal melt spinning problem using the 
///axisymmetric Navier--Stokes and energy equations by coupling two pre-
///existing classes. 
///The QSteadyAxisymAdvectionDiffusionElement<Nnode1d> with bi-quadratic interpolation 
///for the scalar variable (temperature) and QCrouzeixRaviartElement
///which solves the Navier--Stokes equations
///using bi-quadratic interpolation for the velocities and a discontinuous
///bi-linear interpolation for the pressure. Note that we are free to 
///choose the order in which we store the variables at the nodes. In this
///case we choose to store the variables in the order fluid velocities
///followed by temperature. We must, therefore, overload the function
///AxisymAdvectionDiffusionEquations::u_index_axisym_adv_diff() to indicate that
///the temperature is stored at the 3-th position not the 0-th. We do not
///need to overload the corresponding function in the 
///AxisymmetricNavierStokesEquations class because the velocities are stored
///first.
//================================================================
class NonIsothermalAxisymmetricQCrouzeixRaviartElement :
 public virtual QSteadyAxisymAdvectionDiffusionElement<3> ,
 public virtual AxisymmetricQCrouzeixRaviartElement
{

 public:


 /// Function pointer to the function that specifies the
 /// viscosity ratio as function of the temperature. 
 typedef void (*ViscosityRatioFctPt)(double& temperature,
                                     double& result);

 /// Constructor: call the QSteadyAxisymAdvectionDiffusionElement
 /// and AxisymmetricQCrouzeixRaviartElement
 NonIsothermalAxisymmetricQCrouzeixRaviartElement() :  
  QSteadyAxisymAdvectionDiffusionElement<3>(),
  AxisymmetricQCrouzeixRaviartElement(),
  Viscosity_ratio_fct_pt(0)
  {}


 ///The required number of values stored at the nodes is the sum of the
 ///required values of the two single-physics  elements. Note that this step is
 ///generic for any multi-physics element of this type.
 unsigned required_nvalue(const unsigned &n) const
  {
   return ( QSteadyAxisymAdvectionDiffusionElement<3>::required_nvalue(n) +
            AxisymmetricQCrouzeixRaviartElement::required_nvalue(n) );
  }

 ///Overload the standard output function with the broken default
 void output(ostream &outfile) 
  {
   FiniteElement::output(outfile);
  }

 /// Output function:  
 ///  Output r, z, u, w, v, p, theta at Nplot^2 plot points
 // Start of output function
 void output(ostream &outfile, const unsigned &nplot)
  {
   //vector of local coordinates
   Vector<double> s(2);
   
   // Tecplot header info
   outfile << this->tecplot_zone_string(nplot);
   
   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,nplot,s);
     
     // Output the position (r,z) of the plot point
     for(unsigned i=0;i<2;i++) 
      {
       outfile << this->interpolated_x(s,i) << " ";
      }
     
     // Output the fluid velocities (u,w,v) at the plot point
     for(unsigned i=0;i<3;i++) 
      {
       outfile << this->interpolated_u_axi_nst(s,i) << " ";
      }
     
     // Output the fluid pressure at the plot point
     outfile << this->interpolated_p_axi_nst(s)  << " ";
  
     // Output the temperature (the advected variable) at the plot point
     outfile << this->interpolated_u_adv_diff(s) << std::endl;   
    }
   outfile << std::endl;
   
   // Write tecplot footer (e.g. FE connectivity lists)
   this->write_tecplot_zone_footer(outfile,nplot);
  } //End of output function


 /// C-style output function: Broken default
 void output(FILE* file_pt)
  {
   FiniteElement::output(file_pt);
  }

 ///  C-style output function: Broken default
 void output(FILE* file_pt, const unsigned &n_plot)
  {
   FiniteElement::output(file_pt,n_plot);
  }

 /// Output function for an exact solution: Broken default
 void output_fct(ostream &outfile, 
                 const unsigned &Nplot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
   FiniteElement::output_fct(outfile,Nplot,exact_soln_pt);
  }


 ///Overload the index at which the temperature 
 ///variable is stored. We choose to store it after the fluid velocities.
 unsigned u_index_axisym_adv_diff() const 
  {
   return 3;
  }

 /// Validate against exact solution.
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points and compute L2 error
 /// and L2 norm of velocity solution over element
 /// Call the broken default
 void compute_error(ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm)
 {
  FiniteElement::compute_error(outfile,exact_soln_pt,error,norm);
 }


 /// Access function: Pointer to viscosity ratio function
 ViscosityRatioFctPt& viscosity_ratio_fct_pt() {return Viscosity_ratio_fct_pt;}

 /// Access function: Pointer to viscosity ratio function. Const version
 ViscosityRatioFctPt viscosity_ratio_fct_pt() const 
  {return Viscosity_ratio_fct_pt;}


 /// Overload the viscosity ratio in the Navier-Stokes equations
 /// This provides the coupling from axisymmetric advection-diffusion equations
 /// to the Axisymmetric Navier--Stokes equations.
 inline virtual void get_viscosity_ratio_axisym_nst(const unsigned& ipt,
                                                    const Vector<double> &s, 
                                                    const Vector<double> &x,
                                                    double &visco)
  {
   //If the function pointer is zero return zero
   if (Viscosity_ratio_fct_pt == 0)
    {
     visco = *Viscosity_Ratio_pt;
    }
   //Otherwise call the function
   else
    {
     //Get the temperature
     double interpolated_t = this->interpolated_u_adv_diff(s);

     // Valuate viscosity ratio
     (*Viscosity_ratio_fct_pt)(interpolated_t,visco);
    }

/*    //Get the temperature */
/*    const double interpolated_t = this->interpolated_u_adv_diff(s); */

   
/*    // Parameters that modeling the viscosity */
/*    double G = 1.0; */
/*    double eta = 0.0; */

/*    visco = G*exp(eta*(1.0-interpolated_t)); */

  } // end of get_visco_axisym_nst



 /// Overload the wind function in the advection-diffusion equations.
 /// This provides the coupling from the Navier--Stokes equations to the
 /// advection-diffusion equations because the wind is the fluid velocity.
 void get_wind_axisym_adv_diff(const unsigned& ipt,
                               const Vector<double> &s, 
                               const Vector<double>& x,
                               Vector<double>& wind) const
 {
   //The wind function is simply the velocity at the points
   this->interpolated_u_axi_nst(s,wind);
 } // end of get_wind_axisym_adv_diff

 /// Calculate the element's contribution to the residual vector.
 /// Recall that fill_in_* functions MUST NOT initialise the entries 
 /// in the vector to zero. This allows us to call the 
 /// fill_in_* functions of the constituent single-physics elements
 /// sequentially, without wiping out any previously computed entries.
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Fill in the residuals of the Navier-Stokes equations
   AxisymmetricNavierStokesEquations::fill_in_contribution_to_residuals(residuals);

   //Fill in the residuals of the advection-diffusion eqautions
   SteadyAxisymAdvectionDiffusionEquations::fill_in_contribution_to_residuals(residuals);
  }


//-----------Finite-difference the entire jacobian-----------------------
//-----------------------------------------------------------------------
#ifdef USE_FD_JACOBIAN_FOR_NONISOTHERMAL_AXISYMMETRIC_NAVIER_STOKES


 ///Compute the element's residual vector and the Jacobian matrix.
 /// Jacobian is computed by finite-differencing.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
  {
   // This function computes the Jacobian by finite-differencing
   FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
  }

#else

//--------Finite-difference off-diagonal-blocks in jacobian--------------
//-----------------------------------------------------------------------
#ifdef USE_OFF_DIAGONAL_FD_JACOBIAN_FOR_NONISOTHERMAL_AXISYMMETRIC_NAVIER_STOKES

 ///Helper function to get the off-diagonal blocks of the Jacobian
 ///matrix by finite differences
 void fill_in_off_diagonal_jacobian_blocks_by_fd(Vector<double> &residuals,
                                                 DenseMatrix<double> &jacobian)
  {
   //Local storage for the index in the nodes at which the
   //Navier-Stokes velocities are stored (we know that this should be 0,1,2)
   unsigned u_nodal_nst[2];
   for(unsigned i=0;i<2;i++) 
    {
     u_nodal_nst[i] = this->u_index_axi_nst(i);
    }

   //Local storage for the  index at which the temperature is stored
   unsigned u_nodal_adv_diff = this->u_index_axisym_adv_diff();

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
     //There are 2 values of the velocities
     for(unsigned i=0;i<2;i++)
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
         //which must be done using fill_in_contribution because
         //get_residuals will always return the full residuals
         //because the appropriate fill_in function is overloaded above
         for(unsigned m=0;m<n_dof;m++) {newres[m] = 0.0;}
         SteadyAxisymAdvectionDiffusionEquations::
          fill_in_contribution_to_residuals(newres);

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
        //which must be done using fill_in_contribution because
        //get_residuals will always return the full residuals
        //because the appropriate fill_in function is overloaded above
        for(unsigned m=0;m<n_dof;m++) {newres[m] = 0.0;}
        AxisymmetricNavierStokesEquations::
         fill_in_contribution_to_residuals(newres);
        
        //Now fill in the Navier-Stokes sub-block
        for(unsigned m=0;m<n_node;m++)
         {
          //Loop over the fluid velocities
          for(unsigned j=0;j<2;j++)
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


 ///Compute the element's residual Vector and the Jacobian matrix.
 /// Use finite-differencing only for the off-diagonal blocks.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
  {

   //Calculate the Navier-Stokes contributions (diagonal block and residuals)
   AxisymmetricNavierStokesEquations::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //Calculate the advection-diffusion contributions 
   //(diagonal block and residuals)
   SteadyAxisymAdvectionDiffusionEquations::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //We now fill in the off-diagonal (interaction) blocks by finite
   //differences.
   fill_in_off_diagonal_jacobian_blocks_by_fd(residuals,jacobian);
  } //End of jacobian calculation


 //--------------------Analytic jacobian---------------------------------
 //----------------------------------------------------------------------
#else

 ///Helper function to get the off-diagonal blocks of the Jacobian
 ///matrix analytically
 void fill_in_off_diagonal_jacobian_blocks_analytic(Vector<double> &residuals, 
                                                    DenseMatrix<double> &jacobian)
  {
   //We now fill in the off-diagonal (interaction) blocks analytically
   //This requires knowledge of exactly how the residuals are assembled
   //within the parent elements and involves yet another loop over
   //the integration points!
   
   //Local storage for the index in the nodes at which the
   //Navier-Stokes velocities are stored (we know that this should be 0,1,2)
   unsigned u_nodal_nst[2];
   for(unsigned i=0;i<2;i++) 
    {
     u_nodal_nst[i] = this->u_index_axi_nst(i);
    }
   
   //Local storage for the  index at which the temperature is stored
   const unsigned u_nodal_adv_diff = this->u_index_axisym_adv_diff();

   //Find out how many nodes there are
   const unsigned n_node = this->nnode();
   
   //Set up memory for the shape and test functions and their derivatives
   Shape psif(n_node), testf(n_node);
   DShape dpsifdx(n_node,2), dtestfdx(n_node,2);
   
   //Number of integration points
   const unsigned n_intpt = this->integral_pt()->nweight();
   
   //Get Physical Variables from Element
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
     double J = this->dshape_and_dtest_eulerian_at_knot_axi_nst(ipt,
                                                                psif,
                                                                dpsifdx,
                                                                testf,
                                                                dtestfdx);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;
     
     //Calculate local values of temperature derivatives
     //Allocate
     Vector<double> interpolated_du_adv_diff_dx(2,0.0);
     
     // Loop over nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       //Get the nodal value
       double u_value = this->raw_nodal_value(l,u_nodal_adv_diff);
       //Loop over the derivative directions
       for(unsigned j=0;j<2;j++)
        {
         interpolated_du_adv_diff_dx[j] += u_value*dpsifdx(l,j);
        }
      }
     
     //Assemble the jacobian terms
     
     //Loop over the test functions
     for(unsigned l=0;l<n_node;l++)
      {       
       //Assemble the contributions of the fluid velocities to the
       //advection-diffusion equation for the temperature
       
        local_eqn = this->nodal_local_eqn(l,u_nodal_adv_diff);
        //If it's not pinned
        if(local_eqn >= 0)
         {
          //Loop over the shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
            //Loop over the velocity degrees of freedom
            for(unsigned i2=0;i2<2;i2++)
             {
              //Get the local unknown
              local_unknown = this->nodal_local_eqn(l2,u_nodal_nst[i2]);
              //If it's not pinned
              if(local_unknown >= 0)
               {
                //Add the contribution to the jacobian matrix (It is necessary multiply by r)
                jacobian(local_eqn,local_unknown)
                 -= Pe*psif(l2)*interpolated_du_adv_diff_dx[i2]*testf(l)*W;
               }
             }
           }
         }
       
      }
    }
  }


 ///Compute the element's residual Vector and the Jacobian matrix.
 /// Use analytic expressions for the off-diagonal blocks
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
  {
   
   //Calculate the Navier-Stokes contributions (diagonal block and residuals)
   AxisymmetricNavierStokesEquations::fill_in_contribution_to_jacobian(residuals,jacobian);

   //Calculate the advection-diffusion contributions 
   //(diagonal block and residuals)
   SteadyAxisymAdvectionDiffusionEquations::fill_in_contribution_to_jacobian(residuals,jacobian);

   //Fill in the off diagonal blocks analytically
   fill_in_off_diagonal_jacobian_blocks_analytic(residuals,jacobian);
  }

 
#endif
#endif

  private:


 /// Function pointer to the function that specifies the viscosity ratio
 /// as a function of the temperature
 ViscosityRatioFctPt Viscosity_ratio_fct_pt;

};


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

namespace oomph
{
 
//=======================================================================
/// Face geometry for the 
/// NonIsothermalAxisymmetricQCrouzeixRaviartElementelements: 
/// The spatial dimension of the face elements is one lower than that 
/// of the bulk element but they have the same number of points along 
/// their 1D edges.
//=======================================================================
template<>
class FaceGeometry<NonIsothermalAxisymmetricQCrouzeixRaviartElement> : 
 public virtual QElement<1,3>
 {
  
 public:
  
  /// Constructor: Call the constructor for the
  /// appropriate lower-dimensional QElement
  FaceGeometry() : QElement<1,3>() {}
  
 };

template<>
class FaceGeometry<FaceGeometry<NonIsothermalAxisymmetricQCrouzeixRaviartElement> >: 
 public virtual PointElement
 {
  public:
  FaceGeometry() : PointElement() {}
 };
 
}
