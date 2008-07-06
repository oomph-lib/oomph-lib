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
//Non-inline functions for Poisson elements
#include "poisson_elements.h"


namespace oomph
{


//======================================================================
/// Set the data for the number of Variables at each node, always one
/// in every case
//======================================================================
 template<unsigned DIM, unsigned NNODE_1D>
 const unsigned QPoissonElement<DIM,NNODE_1D>::Initial_Nvalue = 1;


//======================================================================
/// Compute element residual Vector and/or element Jacobian matrix 
/// 
/// flag=1: compute both
/// flag=0: compute only residual Vector
///
/// Pure version without hanging nodes
//======================================================================
template <unsigned DIM>
void  PoissonEquations<DIM>::
fill_in_generic_residual_contribution_poisson(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              unsigned flag) 
{
 //Find out how many nodes there are
 const unsigned n_node = nnode();

 //Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

 //Index at which the poisson unknown is stored
 const unsigned u_nodal_index = u_index_poisson();
 
 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();

 //Integers to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Call the derivatives of the shape and test functions
   double J = dshape_and_dtest_eulerian_at_knot_poisson(ipt,psi,dpsidx,
                                                        test,dtestdx);
       
   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Calculate local values of unknown
   //Allocate and initialise to zero
   double interpolated_u=0.0;
   Vector<double> interpolated_x(DIM,0.0);
   Vector<double> interpolated_dudx(DIM,0.0);
   
   //Calculate function value and derivatives:
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Get the nodal value of the poisson unknown
     double u_value = raw_nodal_value(l,u_nodal_index);
     interpolated_u += u_value*psi(l);
     // Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       interpolated_x[j] += raw_nodal_position(l,j)*psi(l);
       interpolated_dudx[j] += u_value*dpsidx(l,j);
      }
    }


   //Get source function
   //-------------------
   double source;
   get_source_poisson(interpolated_x,source);

   // Assemble residuals and Jacobian
   //--------------------------------
       
   // Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     //Get the local equation
     local_eqn = nodal_local_eqn(l,u_nodal_index);
     /*IF it's not a boundary condition*/
     if(local_eqn >= 0)
      {
       // Add body force/source term here 
       residuals[local_eqn] += source*test(l)*W;
             
       // The Poisson bit itself
       for(unsigned k=0;k<DIM;k++)
        {
         residuals[local_eqn] += interpolated_dudx[k]*dtestdx(l,k)*W;
        }

       // Calculate the jacobian
       //-----------------------
       if(flag)
        {
         //Loop over the velocity shape functions again
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           local_unknown = nodal_local_eqn(l2,u_nodal_index);
           //If at a non-zero degree of freedom add in the entry
           if(local_unknown >= 0)
            {
             //Add contribution to Elemental Matrix
             for(unsigned i=0;i<DIM;i++)
              {
               jacobian(local_eqn,local_unknown) 
                += dpsidx(l2,i)*dtestdx(l,i)*W;
              }
            }
          }
        }
      }
    }

  } // End of loop over integration points
}   


//======================================================================
/// Compute derivatives of elemental residual vector with respect
/// to nodal coordinates. 
/// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
/// Overloads the FD-based version in the FE base class.
//======================================================================
template <unsigned DIM>
void  PoissonEquations<DIM>::get_dresidual_dnodal_coordinates(
 RankThreeTensor<double>&
 dresidual_dnodal_coordinates)
{

 //Find out how many nodes there are
 const unsigned n_node = nnode();

 //Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);
 DShape dpsidx_pls(n_node,DIM), dtestdx_pls(n_node,DIM);

 // Deriatives of shape fct derivatives w.r.t. nodal coords
 RankFourTensor<double> d_dpsidx_dX(DIM,n_node,n_node,DIM);
 RankFourTensor<double> d_dtestdx_dX(DIM,n_node,n_node,DIM);

 // Derivative of Jacobian of mapping w.r.t. to nodal coords
 DenseMatrix<double> dJ_dX(DIM,n_node);

 // Derivatives of derivative of u w.r.t. nodal coords
 RankThreeTensor<double> d_dudx_dX(DIM,n_node,DIM);

 // Gradient of source fct
 Vector<double> d_source_dx(DIM);

 //Index at which the poisson unknown is stored
 const unsigned u_nodal_index = u_index_poisson();
 
 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();

 //Integers to store the local equation number
 int local_eqn=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Call the derivatives of the shape and test functions
   double J = dshape_and_dtest_eulerian_at_knot_poisson(ipt,psi,dpsidx,
                                                        test,dtestdx);
       
   //Calculate local values 
   //Allocate and initialise to zero
   Vector<double> interpolated_x(DIM,0.0);
   Vector<double> interpolated_dudx(DIM,0.0);
   
   //Calculate function value and derivatives:
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Get the nodal value of the Poisson unknown
     double u_value = raw_nodal_value(l,u_nodal_index);
     // Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       interpolated_x[j] += raw_nodal_position(l,j)*psi(l);
       interpolated_dudx[j] += u_value*dpsidx(l,j);
      }
    }

   //Get source function
   //-------------------
   double source;
   get_source_poisson(interpolated_x,source);

   // FD step 
   double eps_fd=GeneralisedElement::Default_fd_jacobian_step;
   
   // Do FD loop
   for (unsigned jj=0;jj<n_node;jj++)
    {
     // Get node
     Node* nod_pt=node_pt(jj);
     
     // Loop over coordinate directions
     for (unsigned ii=0;ii<DIM;ii++)
      {
       // Make backup
       double backup=nod_pt->x(ii);
       
       // Do FD step. No node update required as we're
       // attacking the coordinate directly...
       nod_pt->x(ii)+=eps_fd;
       
       //Call the derivatives of the shape and test functions
       //at advanced level
       double J_pls = 
        dshape_and_dtest_eulerian_at_knot_poisson(ipt,psi,dpsidx_pls,
                                                  test,dtestdx_pls);
       
       // Assign
       dJ_dX(ii,jj)=(J_pls-J)/eps_fd;
       for (unsigned i=0;i<DIM;i++)
        {
         for (unsigned j=0;j<n_node;j++)
          {
           d_dpsidx_dX(ii,jj,j,i)=(dpsidx_pls(j,i)-dpsidx(j,i))/eps_fd;
           d_dtestdx_dX(ii,jj,j,i)=(dtestdx_pls(j,i)-dtestdx(j,i))/eps_fd;
          }
        }

       // Shape deriv of du/dx_i
       for (unsigned i=0;i<DIM;i++)
        {
         double aux=0.0;
         for (unsigned j_nod=0;j_nod<n_node;j_nod++)
          {
           aux+=raw_nodal_value(j_nod,u_nodal_index)*
            d_dpsidx_dX(ii,jj,j_nod,i);
          }
         d_dudx_dX(ii,jj,i)=aux;
        }
  
       // Reset coordinate. No node update required as we're
       // attacking the coordinate directly...
       nod_pt->x(ii)=backup;
      }
    }

   // Get gradient of source function
   get_source_gradient_poisson(interpolated_x, d_source_dx);


   // Assemble shape derivatives
   //---------------------------
       
   // Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     //Get the local equation
     local_eqn = nodal_local_eqn(l,u_nodal_index);

     /*IF it's not a boundary condition*/
     if(local_eqn >= 0)
      {
       // Loop over coordinate directions
       for (unsigned ii=0;ii<DIM;ii++)
        {              
         // Loop over nodes
         for (unsigned jj=0;jj<n_node;jj++)
          {       
           double sum=source*psi(l)*dJ_dX(ii,jj)+
            d_source_dx[ii]*psi(l)*psi(jj)*J;

           for (unsigned k=0;k<DIM;k++)
            {
             sum+=interpolated_dudx[k]*(dtestdx(l,k)*dJ_dX(ii,jj)+
                                        d_dtestdx_dX(ii,jj,l,k)*J)
              + d_dudx_dX(ii,jj,k)*dtestdx(l,k)*J;
            }

           // Multiply through by integration weight
           dresidual_dnodal_coordinates(local_eqn,ii,jj)+=sum*w;
          }
        }
      }
    }

  } // End of loop over integration points
}   


//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
template <unsigned DIM>
unsigned  PoissonEquations<DIM>::self_test()
{

 bool passed=true;

 // Check lower-level stuff
 if (FiniteElement::self_test()!=0)
  {
   passed=false;
  }

 // Return verdict
 if (passed)
  {
   return 0;
  }
 else
  {
   return 1;
  }
   
}



//======================================================================
/// Output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void  PoissonEquations<DIM>::output(std::ostream &outfile, 
                                    const unsigned &nplot)
{

 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   for(unsigned i=0;i<DIM;i++) 
    {
     outfile << interpolated_x(s,i) << " ";
    }
   outfile << interpolated_u_poisson(s) << std::endl;   
   
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);

}


//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void  PoissonEquations<DIM>::output(FILE* file_pt,
                                    const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Tecplot header info
 fprintf(file_pt,"%s",tecplot_zone_string(nplot).c_str());

 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   for(unsigned i=0;i<DIM;i++) 
    {
     fprintf(file_pt,"%g ",interpolated_x(s,i));
    }
   fprintf(file_pt,"%g \n",interpolated_u_poisson(s));
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(file_pt,nplot);
}



//======================================================================
 /// Output exact solution
 /// 
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points.
 ///
 ///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
template <unsigned DIM>
void PoissonEquations<DIM>::output_fct(std::ostream &outfile, 
                                       const unsigned &nplot, 
                  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
  
  // Vector for coordintes
  Vector<double> x(DIM);
  
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(1);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,u_exact
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   outfile << exact_soln[0] << std::endl;  
  }
 
 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);
}




//======================================================================
 /// Validate against exact solution
 /// 
 /// Solution is provided via function pointer.
 /// Plot error at a given number of plot points.
 ///
//======================================================================
template <unsigned DIM>
void PoissonEquations<DIM>::compute_error(std::ostream &outfile, 
                                          FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                                          double& error, double& norm)
{ 
 
 // Initialise
 error=0.0;
 norm=0.0;
 
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Vector for coordintes
 Vector<double> x(DIM);
 
 //Find out how many nodes there are in the element
 unsigned n_node = nnode();
 
 Shape psi(n_node);
 
 //Set the value of n_intpt
 unsigned n_intpt = integral_pt()->nweight();
  
 // Tecplot 
 outfile << "ZONE" << std::endl;
 
 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(1);
 
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   
   //Assign values of s
   for(unsigned i=0;i<DIM;i++)
    {
     s[i] = integral_pt()->knot(ipt,i);
    }
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   // Get jacobian of mapping
   double J=J_eulerian(s);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   // Get FE function value
   double u_fe=interpolated_u_poisson(s);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,error
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   outfile << exact_soln[0] << " " << exact_soln[0]-u_fe << std::endl;  
   
   // Add to error and norm
   norm+=exact_soln[0]*exact_soln[0]*W;
   error+=(exact_soln[0]-u_fe)*(exact_soln[0]-u_fe)*W;
   
  }
}





//====================================================================
// Force build of templates
//====================================================================
template class QPoissonElement<1,2>;
template class QPoissonElement<1,3>;
template class QPoissonElement<1,4>;

template class QPoissonElement<2,2>;
template class QPoissonElement<2,3>;
template class QPoissonElement<2,4>;

template class QPoissonElement<3,2>;
template class QPoissonElement<3,3>;
template class QPoissonElement<3,4>;

}
