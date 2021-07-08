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
// Non-inline functions for elements that solve the equations of linear
// elasticity in cartesian coordinates

#include "linear_elasticity_elements.h"


namespace oomph
{


/// Static default value for timescale ratio (1.0 -- for natural scaling) 
template <unsigned DIM>
double LinearElasticityEquationsBase<DIM>::Default_lambda_sq_value=1.0;


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

//======================================================================
/// Compute the strain tensor at local coordinate s
//======================================================================
 template<unsigned DIM>
 void LinearElasticityEquationsBase<DIM>::get_strain(
   const Vector<double> &s,
   DenseMatrix<double> &strain) const
 {
#ifdef PARANOID
  if ((strain.ncol()!=DIM)||(strain.nrow()!=DIM))
   {
    std::ostringstream error_message;
    error_message << "Strain matrix is " << strain.ncol() << " x " 
                  << strain.nrow() << ", but dimension of the equations is " 
                  << DIM << std::endl;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }

  //Find out how many position types there are
  unsigned n_position_type = this->nnodal_position_type();
  
  if(n_position_type != 1)
   {
    throw OomphLibError(
     "LinearElasticity is not yet implemented for more than one position type",
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
#endif
  

  //Find out how many nodes there are in the element
  unsigned n_node = nnode();
    
  //Find the indices at which the local velocities are stored
  unsigned u_nodal_index[DIM];
  for(unsigned i=0;i<DIM;i++) {u_nodal_index[i] = u_index_linear_elasticity(i);}
  
  //Set up memory for the shape and derivative functions
  Shape psi(n_node);
  DShape dpsidx(n_node,DIM);
  
  //Call the derivatives of the shape functions
  (void) dshape_eulerian(s,psi,dpsidx);
  
  //Calculate interpolated values of the derivative of global position
  DenseMatrix<double> interpolated_dudx(DIM,DIM,0.0);
  
  //Loop over nodes
  for(unsigned l=0;l<n_node;l++) 
   {
    //Loop over velocity components
    for(unsigned i=0;i<DIM;i++)
     {
      //Get the nodal value
      const double u_value = this->nodal_value(l,u_nodal_index[i]);

      //Loop over derivative directions
      for(unsigned j=0;j<DIM;j++)
       {
        interpolated_dudx(i,j) += u_value*dpsidx(l,j);
       }
    }
   }
  
  ///Now fill in the entries of the strain tensor
  for(unsigned i=0;i<DIM;i++)
   {
    //Do upper half of matrix
    //Note that j must be signed here for the comparison test to work
    //Also i must be cast to an int
    for(int j=(DIM-1);j>=static_cast<int>(i);j--)
     {
      //Off diagonal terms
      if(static_cast<int>(i)!=j) 
       {
        strain(i,j) = 
         0.5*(interpolated_dudx(i,j) + interpolated_dudx(j,i));
       }
      //Diagonal terms will including growth factor when it comes back in
      else
       {
        strain(i,i) = interpolated_dudx(i,i);
       }
     }
    //Matrix is symmetric so just copy lower half
    for(int j=(i-1);j>=0;j--)
     {
      strain(i,j) = strain(j,i);
     }
   }
 }
 
 



//======================================================================
/// Compute the Cauchy stress tensor at local coordinate s for 
/// displacement formulation.
//======================================================================
template<unsigned DIM>
void LinearElasticityEquations<DIM>::get_stress(const Vector<double> &s,
                                                DenseMatrix<double> &stress)
 const
{
#ifdef PARANOID
 if ((stress.ncol()!=DIM)||(stress.nrow()!=DIM))
  {
   std::ostringstream error_message;
   error_message << "Stress matrix is " << stress.ncol() << " x " 
                 << stress.nrow() << ", but dimension of the equations is " 
                 << DIM << std::endl;
   throw OomphLibError(error_message.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 // Get strain
 DenseMatrix<double> strain(DIM,DIM);
 this->get_strain(s,strain);
 
 // Now fill in the entries of the stress tensor without exploiting
 // symmetry -- sorry too lazy. This fct is only used for
 // postprocessing anyway...
 for(unsigned i=0;i<DIM;i++)
  {
   for(unsigned j=0;j<DIM;j++)
    {
     stress(i,j) = 0.0;
     for (unsigned k=0;k<DIM;k++)
      {
       for (unsigned l=0;l<DIM;l++)
        {
         stress(i,j)+=this->E(i,j,k,l)*strain(k,l);
        }
      }
    }
  }
}
 
 
//=======================================================================
/// Compute the residuals for the linear elasticity equations in 
/// cartesian coordinates. Flag indicates if we want Jacobian too.
//=======================================================================
 template <unsigned DIM>
 void LinearElasticityEquations<DIM>::
 fill_in_generic_contribution_to_residuals_linear_elasticity(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,unsigned flag)
 {
  //Find out how many nodes there are
  unsigned n_node = this->nnode();
  
#ifdef PARANOID
  //Find out how many positional dofs there are
  unsigned n_position_type = this->nnodal_position_type();
  
  if(n_position_type != 1)
   {
    throw OomphLibError(
     "LinearElasticity is not yet implemented for more than one position type",
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }

  //Throw and error if an elasticity tensor has not been set
  if(this->Elasticity_tensor_pt==0)
   {
    throw OomphLibError(
     "No elasticity tensor set.",
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  //Find the indices at which the local displacements are stored
  unsigned u_nodal_index[DIM];
  for(unsigned i=0;i<DIM;i++) 
   {u_nodal_index[i] = this->u_index_linear_elasticity(i);}
  
  // Timescale ratio (non-dim density)
  double Lambda_sq = this->lambda_sq();
  
  //Set up memory for the shape functions
  Shape psi(n_node);
  DShape dpsidx(n_node,DIM);
  
  //Set the value of Nintpt -- the number of integration points
  unsigned n_intpt = this->integral_pt()->nweight();
  
  //Set the vector to hold the local coordinates in the element
  Vector<double> s(DIM);
  
  //Integer to store the local equation number
  int local_eqn=0, local_unknown=0;
  
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
    Vector<double> accel(DIM,0.0);
    
    //Calculate displacements and derivatives and lagrangian coordinates
    for(unsigned l=0;l<n_node;l++)
     {
      //Loop over displacement components (deformed position)
      for(unsigned i=0;i<DIM;i++)
       {
        //Calculate the Lagrangian coordinates and the accelerations
        interpolated_x[i] += this->raw_nodal_position(l,i)*psi(l);
        
        // Only compute accelerations if inertia is switched on
        if(this->Unsteady)
         {
          accel[i] += this->d2u_dt2_linear_elasticity(l,i)*psi(l);
         }
        
        //Get the nodal displacements
        const double u_value = this->raw_nodal_value(l,u_nodal_index[i]);
        
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
    
//=====EQUATIONS OF LINEAR ELASTICITY ========
    
    //Loop over the test functions, nodes of the element
    for(unsigned l=0;l<n_node;l++)
     {
      //Loop over the displacement components
      for(unsigned a=0;a<DIM;a++)
       {
        //Get the equation number
        local_eqn = this->nodal_local_eqn(l,u_nodal_index[a]);
        /*IF it's not a boundary condition*/
        if(local_eqn >= 0)
         {
          // Acceleration and body force
          residuals[local_eqn] += 
           (Lambda_sq*accel[a]-b[a])*psi(l)*W;
          
          // Stress term
          for(unsigned b=0;b<DIM;b++)
           {
            for(unsigned c=0;c<DIM;c++)
             {
              for(unsigned d=0;d<DIM;d++)
               {
                //Add the stress terms to the residuals
                residuals[local_eqn] +=
                 this->E(a,b,c,d)*interpolated_dudx(c,d)*dpsidx(l,b)*W;
               }
             }
           }

          //Jacobian entries
          if(flag)
           {
            //Loop over the displacement basis functions again
            for(unsigned l2=0;l2<n_node;l2++)
             {
              //Loop over the displacement components again
              for(unsigned c=0;c<DIM;c++)
               {
                local_unknown = this->nodal_local_eqn(l2,u_nodal_index[c]);
                //If it's not pinned
                if(local_unknown >= 0)
                 {
                  // Inertial term
                  if (a==c)
                   {
                    jacobian(local_eqn,local_unknown) +=
                     Lambda_sq*
                     this->node_pt(l2)->time_stepper_pt()->weight(2,0)*
                     psi(l)*psi(l2)*W;
                   }

                  for(unsigned b=0;b<DIM;b++)
                   {
                    for(unsigned d=0;d<DIM;d++)
                     {
                      //Add the contribution to the Jacobian matrix
                      jacobian(local_eqn,local_unknown) +=
                       this->E(a,b,c,d)*dpsidx(l2,d)*dpsidx(l,b)*W;
                     }
                   }
                 } //End of if not boundary condition
                
               }
             }
           } //End of jacobian calculation
          
         } //End of if not boundary condition
       } //End of loop over coordinate directions
     } //End of loop over shape functions
   } //End of loop over integration points
  
 }





//=======================================================================
/// Output exact solution x,y,[z],u,v,[w]
//=======================================================================
 template <unsigned DIM>
 void LinearElasticityEquations<DIM>::output_fct(std::ostream &outfile, 
                                                 const unsigned &nplot, 
                  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Vector for coordintes
 Vector<double> x(DIM);
 
 // Tecplot header info
 outfile << this->tecplot_zone_string(nplot);
 
 // Exact solution Vector 
 Vector<double> exact_soln(DIM);
 
 // Loop over plot points
 unsigned num_plot_points=this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   
   // Get x position as Vector
   this->interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,u_exact,...
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << exact_soln[i] << " ";
    }
   outfile << std::endl;  
  }
 
 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(outfile,nplot);

}



//=======================================================================
/// Output exact solution x,y,[z],u,v,[w] (unsteady version)
//=======================================================================
 template <unsigned DIM>
 void LinearElasticityEquations<DIM>::output_fct(std::ostream &outfile,
                                                 const unsigned &nplot,
                                                 const double &time,
                  FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Vector for coordintes
 Vector<double> x(DIM);
 
 // Tecplot header info
 outfile << this->tecplot_zone_string(nplot);
 
 // Exact solution Vector 
 Vector<double> exact_soln(DIM);
 
 // Loop over plot points
 unsigned num_plot_points=this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   
   // Get x position as Vector
   this->interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(time,x,exact_soln);
   
   //Output x,y,...,u_exact,...
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << exact_soln[i] << " ";
    }
   outfile << std::endl;  
  }
 
 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(outfile,nplot);

}



//=======================================================================
/// Output: x,y,[z],u,v,[w]
//=======================================================================
 template <unsigned DIM>
 void LinearElasticityEquations<DIM>::output(std::ostream &outfile, 
                                             const unsigned &nplot)
 {
  //Set output Vector
  Vector<double> s(DIM);
  Vector<double> x(DIM);
  Vector<double> u(DIM);
  
  // Tecplot header info
  outfile << this->tecplot_zone_string(nplot);
  
  // Loop over plot points
  unsigned num_plot_points=this->nplot_points(nplot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
   {
    
    // Get local coordinates of plot point
    this->get_s_plot(iplot,nplot,s);
    
    // Get Eulerian coordinates and displacements
    this->interpolated_x(s,x);
    this->interpolated_u_linear_elasticity(s,u);
    
    //Output the x,y,..
    for(unsigned i=0;i<DIM;i++) 
     {outfile << x[i] << " ";}

    // Output u,v,..
    for(unsigned i=0;i<DIM;i++) 
     {outfile << u[i] << " ";} 
    
    outfile << std::endl;
   }
  
  // Write tecplot footer (e.g. FE connectivity lists)
  this->write_tecplot_zone_footer(outfile,nplot);
 }
 



//=======================================================================
/// C-style output: x,y,[z],u,v,[w]
//=======================================================================
template <unsigned DIM>
void LinearElasticityEquations<DIM>::output(FILE* file_pt, 
                                            const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Tecplot header info
 fprintf(file_pt,"%s",this->tecplot_zone_string(nplot).c_str());
 
 // Loop over plot points
 unsigned num_plot_points=this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {   
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   
   // Coordinates
   for(unsigned i=0;i<DIM;i++) 
    {
     fprintf(file_pt,"%g ",this->interpolated_x(s,i));
    }
   
   // Displacement
   for(unsigned i=0;i<DIM;i++) 
    {
     fprintf(file_pt,"%g ",this->interpolated_u_linear_elasticity(s,i));
    }
  }
 fprintf(file_pt,"\n");
 
 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(file_pt,nplot);
 
}


//======================================================================
/// Validate against exact velocity solution
/// Solution is provided via function pointer.
/// Plot at a given number of plot points and compute L2 error
/// and L2 norm of velocity solution over element.
//=======================================================================
template<unsigned DIM>
void LinearElasticityEquations<DIM>::compute_error(
 std::ostream &outfile,
 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
 double& error, double& norm)
{
 
 error=0.0;
 norm=0.0;

 //Vector of local coordinates
 Vector<double> s(DIM);

 // Vector for coordinates
 Vector<double> x(DIM);

 //Set the value of n_intpt
 unsigned n_intpt = this->integral_pt()->nweight();
   
 outfile << "ZONE" << std::endl;
 
 // Exact solution Vector (u,v,[w])
 Vector<double> exact_soln(DIM);
   
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

   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);

   // Displacement error
   for(unsigned i=0;i<DIM;i++)
    {
     norm+=exact_soln[i]*exact_soln[i]*W;
     error+=(exact_soln[i]-this->interpolated_u_linear_elasticity(s,i))*
      (exact_soln[i]-this->interpolated_u_linear_elasticity(s,i))*W;
    }

   //Output x,y,[z]
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }

   //Output u_error,v_error,[w_error]
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << exact_soln[i]-this->interpolated_u_linear_elasticity(s,i) 
             << " ";
    }
   outfile << std::endl;   
  }
}

//======================================================================
/// Validate against exact velocity solution
/// Solution is provided via function pointer.
/// Plot at a given number of plot points and compute L2 error
/// and L2 norm of velocity solution over element. Unsteady version
//=======================================================================
template<unsigned DIM>
void LinearElasticityEquations<DIM>::compute_error(
 std::ostream &outfile,
 FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
 const double& time, double& error, double& norm)
{
 
 error=0.0;
 norm=0.0;

 //Vector of local coordinates
 Vector<double> s(DIM);

 // Vector for coordinates
 Vector<double> x(DIM);

 //Set the value of n_intpt
 unsigned n_intpt = this->integral_pt()->nweight();
   
 outfile << "ZONE" << std::endl;
 
 // Exact solution Vector (u,v,[w])
 Vector<double> exact_soln(DIM);
   
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

   // Get exact solution at this point
   (*exact_soln_pt)(time,x,exact_soln);

   // Displacement error
   for(unsigned i=0;i<DIM;i++)
    {
     norm+=exact_soln[i]*exact_soln[i]*W;
     error+=(exact_soln[i]-this->interpolated_u_linear_elasticity(s,i))*
      (exact_soln[i]-this->interpolated_u_linear_elasticity(s,i))*W;
    }

   //Output x,y,[z]
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }

   //Output u_error,v_error,[w_error]
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << exact_soln[i]-this->interpolated_u_linear_elasticity(s,i) 
             << " ";
    }
   outfile << std::endl;   
  }
}


//Instantiate the required elements
template class LinearElasticityEquationsBase<2>;
template class LinearElasticityEquations<2>;

template class QLinearElasticityElement<3,3>;
template class LinearElasticityEquationsBase<3>;
template class LinearElasticityEquations<3>;



}
