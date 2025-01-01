//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
//Non-inline functions for elements that solve the principle of virtual
//equations of solid mechanics

#include "homo_lin_elasticity_elements.h"

namespace oomph
{

/// Static default value for timescale ratio (1.0 -- for natural scaling) 
double HomogenisedLinearElasticityEquationsBase::Default_lambda_sq_value=1.0;

/// ///////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////

//=======================================================================
/// Compute the residuals for the discretised principle of 
/// virtual displacements.
//=======================================================================
 template<unsigned DIM>
 void HomogenisedLinearElasticityEquations<DIM>::
 fill_in_generic_contribution_to_residuals_linear_elasticity(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,unsigned flag)
 {
  //Find out how many nodes there are
  unsigned n_node = this->nnode();
  
  //Find out how many positional dofs there are
  unsigned n_position_type = this->nnodal_position_type();
  
  if(n_position_type != 1)
   {
    throw OomphLibError(
     "HomogenisedLinearElasticity is not yet implemented for more than one position type",
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }

  //Create pointer for the elasticity tensor
  ElasticityTensor* E_pt;
  
  //Find the indices at which the local velocities are stored
  //Even in "2" dimensions, there are always three displacement components
  unsigned u_nodal_index[3];
  for(unsigned i=0;i<3;i++) 
   {u_nodal_index[i] = this->u_index_linear_elasticity(i);}
  
  //Set up memory for the shape functions
  Shape psi(n_node);
  DShape dpsidx(n_node,DIM);
  
  //Set the value of Nintpt -- the number of integration points
  unsigned n_intpt = this->integral_pt()->nweight();
  
  //Set the vector to hold the local coordinates in the element
  Vector<double> s(DIM);
  
  //Integer to store the local equation number
  int local_eqn=0, local_unknown=0;

  //Read out the values of p and m for the specific cell problem
  const unsigned p = this->get_p(); 
  const unsigned m = this->get_m();

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
    DenseMatrix<double> interpolated_dudx(3,DIM,0.0);
    
    //Calculate displacements and derivatives and Eulerian coordinates
    for(unsigned l=0;l<n_node;l++)
     {
      //Loop over coordinates
      for(unsigned i=0;i<DIM;i++)
       {
        //Calculate the Lagrangian coordinates and the accelerations
        interpolated_x[i] += this->raw_nodal_position(l,i)*psi(l);
       }
     
      //Loop over displacement components (always 3)
      for(unsigned i=0;i<3;i++)
       {
        //Get the nodal displacements
        const double u_value = this->raw_nodal_value(l,u_nodal_index[i]);
        
        //Loop over derivative directions
        for(unsigned j=0;j<DIM;j++)
         {
          interpolated_dudx(i,j) += u_value*dpsidx(l,j);
         }
      }
     }

    //Premultiply the weights and the Jacobian
    double W = w*J; 
    
    //Get a pointer to  the elasticity tensor
    this->get_E_pt(interpolated_x,E_pt);
    
//=====EQUATIONS OF ELASTICITY FROM PRINCIPLE OF VIRTUAL DISPLACEMENTS========
    
    //Loop over the test functions, nodes of the element
    for(unsigned l=0;l<n_node;l++)
     {
      //Loop over the three displacement components
      for(unsigned a=0;a<3;a++)
       {
        //Get the equation number
        local_eqn = this->nodal_local_eqn(l,u_nodal_index[a]);
        /*IF it's not a boundary condition*/
        if(local_eqn >= 0)
         {
          // Stress term
          for(unsigned b=0;b<DIM;b++)
           {
            for(unsigned c=0;c<3;c++)
             {
              for(unsigned d=0;d<DIM;d++)
               {
                //Add the stress terms to the residuals
                residuals[local_eqn] +=
                 (*E_pt)(a,b,c,d)*interpolated_dudx(c,d)*dpsidx(l,b)*W;
               }
             }
            //Add the additional term
            residuals[local_eqn] += (*E_pt)(a,b,p,m)*dpsidx(l,b)*W;
           }
          
          //Jacobian entries
          if(flag)
           {
            //Loop over the displacement basis functions again
            for(unsigned l2=0;l2<n_node;l2++)
             {
              //Loop over the displacement components again
              for(unsigned c2=0;c2<3;c2++)
               {
                local_unknown = this->nodal_local_eqn(l2,u_nodal_index[c2]);
                //If it's not pinned
                if(local_unknown >= 0)
                 {
                  for(unsigned b=0;b<DIM;b++)
                   {
                    for(unsigned d=0;d<DIM;d++)
                     {
                      //Add the contribution to the Jacobian matrix
                      jacobian(local_eqn,local_unknown) +=
                       (*E_pt)(a,b,c2,d)*dpsidx(l2,d)*dpsidx(l,b)*W;
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
/// Compute the current contribution to the effective elasticity tensor
//=======================================================================
 template<unsigned DIM>
 void HomogenisedLinearElasticityEquations<DIM>::
 calculate_effective_modulus(DenseMatrix<double> &H)
 {
  //Find out how many nodes there are
  unsigned n_node = this->nnode();
  
  //Find out how many positional dofs there are
  unsigned n_position_type = this->nnodal_position_type();
  
  if(n_position_type != 1)
   {
    throw OomphLibError(
     "HomogenisedLinearElasticity is not yet implemented for more than one position type",
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }

  //Create pointer for the elasticity tensor
  ElasticityTensor* E_pt;
  
  //Find the indices at which the local velocities are stored
  //Even in "2" dimensions, there are always three displacement components
  unsigned u_nodal_index[3];
  for(unsigned i=0;i<3;i++) 
   {u_nodal_index[i] = this->u_index_linear_elasticity(i);}
  
  //Set up memory for the shape functions
  Shape psi(n_node);
  DShape dpsidx(n_node,DIM);
  
  //Set the value of Nintpt -- the number of integration points
  unsigned n_intpt = this->integral_pt()->nweight();
  
  //Set the vector to hold the local coordinates in the element
  Vector<double> s(DIM);

  //Read out the values of p and m for the specific cell problem
  const unsigned p = this->get_p(); 
  const unsigned m = this->get_m();

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
    DenseMatrix<double> interpolated_dudx(3,DIM,0.0);
    
    //Calculate displacements and derivatives and Eulerian coordinates
    for(unsigned l=0;l<n_node;l++)
     {
      //Loop over coordinates
      for(unsigned i=0;i<DIM;i++)
       {
        //Calculate the Lagrangian coordinates and the accelerations
        interpolated_x[i] += this->raw_nodal_position(l,i)*psi(l);
       }
     
      //Loop over displacement components (always 3)
      for(unsigned i=0;i<3;i++)
       {
        //Get the nodal displacements
        const double u_value = this->raw_nodal_value(l,u_nodal_index[i]);
        
        //Loop over derivative directions
        for(unsigned j=0;j<DIM;j++)
         {
          interpolated_dudx(i,j) += u_value*dpsidx(l,j);
         }
      }
     }

    //Premultiply the weights and the Jacobian
    double W = w*J; 
    
    //Get a pointer to  the elasticity tensor
    this->get_E_pt(interpolated_x,E_pt);

    //Loop over the three displacement components
    for(unsigned a=0;a<3;a++)
     {
      // Stress term (including all 3 terms here)
      for(unsigned b=0;b<3;b++)
       {
        for(unsigned c=0;c<3;c++)
         {
          for(unsigned d=0;d<DIM;d++)
           {
            //Add the stress terms to the residuals
            H(a,b) +=
             (*E_pt)(a,b,c,d)*interpolated_dudx(c,d)*W;
           }
         }
        //Add the additional term
        H(a,b) += (*E_pt)(a,b,p,m)*W;
       }
     } //End of loop over coordinate directions
    
   } //End of loop over integration points
 }


//=======================================================================
/// Output: x,y,[z],xi0,xi1,[xi2],gamma
//=======================================================================
 template<unsigned DIM>
 void HomogenisedLinearElasticityEquations<DIM>::output(std::ostream &outfile, 
                                                        const unsigned &n_plot)
 {
  //Set output Vector
  Vector<double> s(DIM);
  Vector<double> x(DIM);
  Vector<double> u(3);
  
  //Tecplot header info 
  outfile << tecplot_zone_string(n_plot);
  
  // Loop over plot points
  unsigned num_plot_points=nplot_points(n_plot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
   {
    // Get local coordinates of plot point
    get_s_plot(iplot,n_plot,s);
    
    // Get Eulerian and Lagrangian coordinates
    this->interpolated_x(s,x);
    this->interpolated_u_linear_elasticity(s,u);
     
    //Output the x,y,..
    for(unsigned i=0;i<DIM;i++) 
     {outfile << x[i] << " ";}
    // Output xi0,xi1,..
    for(unsigned i=0;i<3;i++) 
     {outfile << u[i] << " ";} 
    outfile << std::endl;
   }
  
  //Write tecplot footer
  write_tecplot_zone_footer(outfile,n_plot);
 }



//=======================================================================
/// C-style output: x,y,[z],xi0,xi1,[xi2],gamma
//=======================================================================
 template<unsigned DIM>
 void HomogenisedLinearElasticityEquations<DIM>::output(FILE* file_pt, 
                                                        const unsigned &n_plot)
 {
  //Set output Vector
  Vector<double> s(DIM);
  Vector<double> x(DIM);
  Vector<double> u(3);
  
  //Tecplot header info 
  fprintf(file_pt,"%s",tecplot_zone_string(n_plot).c_str());
  
  // Loop over plot points
  unsigned num_plot_points=nplot_points(n_plot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
   {
    // Get local coordinates of plot point
    get_s_plot(iplot,n_plot,s);

    // Get Eulerian and Lagrangian coordinates
    this->interpolated_x(s,x);
    this->interpolated_u_linear_elasticity(s,u);
    
    //Output the x,y,..
    for(unsigned i=0;i<DIM;i++) 
     {
      fprintf(file_pt,"%g ",x[i]);
     }
    // Output xi0,xi1,..
    for(unsigned i=0;i<3;i++) 
     {
      //outfile << xi[i] << " ";
      fprintf(file_pt,"%g ",u[i]);
     } 
   }
  
  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(file_pt,n_plot);
 }


 //Force instantiation of the classes
 template class HomogenisedLinearElasticityEquations<2>;
 template class HomogenisedLinearElasticityEquations<3>;

}
