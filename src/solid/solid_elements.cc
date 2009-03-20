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
//Non-inline functions for elements that solve the principle of virtual
//equations of solid mechanics

#include "solid_elements.h"


namespace oomph
{


/// Static default value for timescale ratio (1.0 -- for natural scaling) 
template <unsigned DIM>
double PVDEquationsBase<DIM>::Default_lambda_sq_value=1.0;


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

//======================================================================
/// Compute the strain tensor at local coordinate s
//======================================================================
template<unsigned DIM>
void PVDEquationsBase<DIM>::get_strain(const Vector<double> &s,
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
                       "PVDEquationsBase<DIM>::get_strain()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Find out how many nodes there are in the element
 unsigned n_node = nnode();

 //Find out how many position types there are
 unsigned n_position_type = this->nnodal_position_type();
 
 //Set up memory for the shape and test functions
 Shape psi(n_node,n_position_type);
 DShape dpsidxi(n_node,n_position_type,DIM);
 
 //Call the derivatives of the shape functions
 (void) dshape_lagrangian(s,psi,dpsidxi);
 
 //Calculate interpolated values of the derivative of global position
 DenseMatrix<double> interpolated_G(DIM);

 //Initialise to zero
 for(unsigned i=0;i<DIM;i++)
  {
   for(unsigned j=0;j<DIM;j++) {interpolated_G(i,j) = 0.0;}
  }
 
 //Storage for Lagrangian coordinates (initialised to zero)
 Vector<double> interpolated_xi(DIM,0.0);
 
 //Loop over nodes
 for(unsigned l=0;l<n_node;l++) 
  {
   //Loop over the positional dofs
   for(unsigned k=0;k<n_position_type;k++)
    {
     //Loop over velocity components
     for(unsigned i=0;i<DIM;i++)
      {
       //Calculate the Lagrangian coordinates
       interpolated_xi[i] += this->lagrangian_position_gen(l,k,i)*psi(l,k);
     
       //Loop over derivative directions
       for(unsigned j=0;j<DIM;j++)
        {                               
         interpolated_G(j,i) += 
          this->nodal_position_gen(l,k,i)*dpsidxi(l,k,j);
        }
      }
    }
  }
 
 //Get isotropic growth factor
 double gamma=1.0;
 get_isotropic_growth(s,interpolated_xi,gamma);
 
 // We use Cartesian coordinates as the reference coordinate
 // system. In this case the undeformed metric tensor is always
 // the identity matrix -- stretched by the isotropic growth
 double diag_entry=pow(gamma,2.0/double(DIM)); 
 DenseMatrix<double> g(DIM);
 for(unsigned i=0;i<DIM;i++)
  {
   for(unsigned j=0;j<DIM;j++)
    {
     if(i==j) {g(i,j) = diag_entry;} else {g(i,j) = 0.0;}
    }
  }


 
 //Declare and calculate the deformed metric tensor
 DenseMatrix<double> G(DIM);
 
 //Assign values of G
 for(unsigned i=0;i<DIM;i++)
  {
   //Do upper half of matrix
   for(unsigned j=i;j<DIM;j++) 
    {
     //Initialise G(i,j) to zero
     G(i,j) = 0.0;
     //Now calculate the dot product
     for(unsigned k=0;k<DIM;k++)
      {
       G(i,j) += interpolated_G(i,k)*interpolated_G(j,k);
      }
    }
   //Matrix is symmetric so just copy lower half
   for(unsigned j=0;j<i;j++) 
    {
     G(i,j) = G(j,i);
    }
  }
 
 //Fill in the strain tensor
 for(unsigned i=0;i<DIM;i++)
  {
   for(unsigned j=0;j<DIM;j++) {strain(i,j)= 0.5*(G(i,j) - g(i,j));}
  }
 
}

//=======================================================================
/// Compute the residuals for the discretised principle of 
/// virtual displacements.
//=======================================================================
template <unsigned DIM>
void PVDEquations<DIM>::
fill_in_generic_contribution_to_residuals_pvd(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              const unsigned& flag)
{

#ifdef PARANOID
 // Check if the constitutive equation requires the explicit imposition of an
 // incompressibility constraint
 if (this->Constitutive_law_pt->requires_incompressibility_constraint())
  {
   throw OomphLibError(
    "PVDEquations cannot be used with incompressible constitutive laws.",
    "PVDEquations<DIM>::fill_in_generic_contribution_to_residuals_pvd()",
    OOMPH_EXCEPTION_LOCATION);   
  }
#endif

 // Simply set up initial condition?
 if (this->Solid_ic_pt!=0)
  {
   this->fill_in_residuals_for_solid_ic(residuals);
   return;
  }

 // Storage for symmetric part of the Jacobian
 unsigned n_dof=jacobian.nrow();
 DenseMatrix<double> sym_jacobian(n_dof,n_dof,0.0);
 
 // Storage of eqn numbers
 Vector<unsigned> temp_eqn_numbers;
 temp_eqn_numbers.reserve(n_dof);

 //Find out how many nodes there are
 unsigned n_node = this->nnode();

 //Find out how many positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();
 
 //Set up memory for the shape functions
 Shape psi(n_node,n_position_type);
 DShape dpsidxi(n_node,n_position_type,DIM);

 //Set the value of Nintpt -- the number of integration points
 unsigned n_intpt = this->integral_pt()->nweight();
   
 //Set the vector to hold the local coordinates in the element
 Vector<double> s(DIM);
 
 // Timescale ratio (non-dim density)
 double lambda_sq = this->lambda_sq();
 
 // Time factor
 double time_factor=0.0;
 if (lambda_sq>0)
  {
   time_factor=this->node_pt(0)->position_time_stepper_pt()->weight(2,0);
  }
 
 //Integer to store the local equation number
 int local_eqn=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Assign the values of s
   for(unsigned i=0;i<DIM;++i) {s[i] = this->integral_pt()->knot(ipt,i);}

   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);

   //Call the derivatives of the shape functions (and get Jacobian)
   double J = this->dshape_lagrangian_at_knot(ipt,psi,dpsidxi);

   //Calculate interpolated values of the derivative of global position
   //wrt lagrangian coordinates
   DenseMatrix<double> interpolated_G(DIM);

   // Setup memory for accelerations
   Vector<double> accel(DIM);
      
   //Initialise to zero
   for(unsigned i=0;i<DIM;i++)
    {
     // Initialise acclerations
     accel[i]=0.0;
     for(unsigned j=0;j<DIM;j++)
     {
      interpolated_G(i,j) = 0.0;
     }
    }

   //Storage for Lagrangian coordinates (initialised to zero)
   Vector<double> interpolated_xi(DIM,0.0);

   //Calculate displacements and derivatives and lagrangian coordinates
   for(unsigned l=0;l<n_node;l++)
    {
     //Loop over positional dofs
     for(unsigned k=0;k<n_position_type;k++)
      {
       //Loop over displacement components (deformed position)
       for(unsigned i=0;i<DIM;i++)
        {
         //Calculate the Lagrangian coordinates and the accelerations
         interpolated_xi[i] += this->lagrangian_position_gen(l,k,i)*psi(l,k);

         // Only compute accelerations if inertia is switched on
         if ((lambda_sq>0.0)&&(this->unsteady()))
          {
           accel[i] += this->dnodal_position_gen_dt(2,l,k,i)*psi(l,k);
          }

         //Loop over derivative directions
         for(unsigned j=0;j<DIM;j++)
          {
           interpolated_G(j,i) += 
            this->nodal_position_gen(l,k,i)*dpsidxi(l,k,j);
          }
        }
      }
    }
       
   //Get isotropic growth factor
   double gamma=1.0;
   this->get_isotropic_growth(s,interpolated_xi,gamma);


   //Get body force at current time
   Vector<double> b(DIM);
   this->body_force(interpolated_xi,b);

   // We use Cartesian coordinates as the reference coordinate
   // system. In this case the undeformed metric tensor is always
   // the identity matrix -- stretched by the isotropic growth
   double diag_entry=pow(gamma,2.0/double(DIM)); 
   DenseMatrix<double> g(DIM);
   for(unsigned i=0;i<DIM;i++)
    {
     for(unsigned j=0;j<DIM;j++)
      {
       if(i==j) {g(i,j) = diag_entry;}
       else {g(i,j) = 0.0;}
      }
    }

   //Premultiply the undeformed volume ratio (from the isotropic
   // growth), the weights and the Jacobian
   double W = gamma*w*J; 

   //Declare and calculate the deformed metric tensor
   DenseMatrix<double> G(DIM);

   //Assign values of G
   for(unsigned i=0;i<DIM;i++)
    {
     //Do upper half of matrix
     for(unsigned j=i;j<DIM;j++) 
      {
       //Initialise G(i,j) to zero
       G(i,j) = 0.0;
       //Now calculate the dot product
       for(unsigned k=0;k<DIM;k++)
        {
         G(i,j) += interpolated_G(i,k)*interpolated_G(j,k);
        }
      }
     //Matrix is symmetric so just copy lower half
     for(unsigned j=0;j<i;j++) 
      {
       G(i,j) = G(j,i);
      }
    }

   //Now calculate the stress tensor from the constitutive law
   DenseMatrix<double> sigma(DIM);
   get_stress(g,G,sigma);



   // Get stress derivative by FD only needed for Jacobian
   //-----------------------------------------------------

   // Stress derivative
   RankFourTensor<double> d_stress_dG(DIM,DIM,DIM,DIM,0.0);
   
   // Derivative of metric tensor w.r.t. to nodal coords
   RankFiveTensor<double> d_G_dX(n_node,n_position_type,DIM,DIM,DIM,0.0);
   
   // Get Jacobian too?
   if (flag==1)
    {     
     // Derivative of metric tensor w.r.t. to discrete positional dofs
     // NOTE: Since G is symmetric we only compute the upper triangle
     //       and DO NOT copy the entries across. Subsequent computations
     //       must (and, in fact, do) therefore only operate with upper
     //       triangular entries
     for (unsigned ll=0;ll<n_node;ll++)
      {
       for (unsigned kk=0;kk<n_position_type;kk++)
        {         
         for (unsigned ii=0;ii<DIM;ii++)
          {
           for (unsigned aa=0;aa<DIM;aa++)
            {
             for (unsigned bb=aa;bb<DIM;bb++) 
              {
               d_G_dX(ll,kk,ii,aa,bb)=
                interpolated_G(aa,ii)*dpsidxi(ll,kk,bb)+
                interpolated_G(bb,ii)*dpsidxi(ll,kk,aa);
              }
            }
          }
        }
      }

     // FD step 
     double eps_fd=GeneralisedElement::Default_fd_jacobian_step;
     
     //Advanced metric tensor
     DenseMatrix<double> G_pls(DIM,DIM);
     DenseMatrix<double> sigma_pls(DIM,DIM);
     
     // Copy across
     for (unsigned i=0;i<DIM;i++)
      {
       for (unsigned j=0;j<DIM;j++)
        {
         G_pls(i,j)=G(i,j);
        }
      }
     
     // Do FD -- only w.r.t. to upper indices, exploiting symmetry.
     // NOTE: We exploit the symmetry of the stress and metric tensors
     //       by incrementing G(i,j) and G(j,i) simultaenously and
     //       only fill in the "upper" triangles without copying things
     //       across the lower triangle. This is taken into account
     //       in the remaining code further below.
     for(unsigned i=0;i<DIM;i++)
      {
       for(unsigned j=i;j<DIM;j++) 
        {
         G_pls(i,j) += eps_fd;
         G_pls(j,i) = G_pls(i,j);
         
         // Get advanced stress
         get_stress(g,G_pls,sigma_pls);
         
         for (unsigned ii=0;ii<DIM;ii++)
          {
           for (unsigned jj=ii;jj<DIM;jj++) 
            {
             d_stress_dG(ii,jj,i,j)=(sigma_pls(ii,jj)-sigma(ii,jj))/eps_fd;
            }   
          }
         
         // Reset 
         G_pls(i,j) = G(i,j);
         G_pls(j,i) = G(j,i);

        }
      }

    }
   
//=====EQUATIONS OF ELASTICITY FROM PRINCIPLE OF VIRTUAL DISPLACEMENTS========
       
   //Loop over the test functions, nodes of the element
   for(unsigned l=0;l<n_node;l++)
    {
     //Loop of types of dofs
     for(unsigned k=0;k<n_position_type;k++)
      {
       // Offset for faster access
       const unsigned offset5=dpsidxi.offset(l ,k);

       //Loop over the displacement components
       for(unsigned i=0;i<DIM;i++)
        {
         //Get the equation number
         local_eqn = this->position_local_eqn(l,k,i);

         /*IF it's not a boundary condition*/
         if(local_eqn >= 0)
          {
           // Record eqn number
           if (ipt==0) temp_eqn_numbers.push_back(local_eqn);

           double sum=0.0;

           // Acceleration and body force
           sum+=(lambda_sq*accel[i]-b[i])*psi(l,k);

           // Stress term
           for(unsigned a=0;a<DIM;a++)
            {
             unsigned count=offset5;
             for(unsigned b=0;b<DIM;b++) 
              {
               //Add the stress terms to the residuals
               sum+=sigma(a,b)*interpolated_G(a,i)*
                dpsidxi.raw_direct_access(count);
               ++count;
              }
            }
           residuals[local_eqn] += W*sum;

           // Get Jacobian too?
           if (flag==1)
            {

             // Offset for faster access in general stress loop
             const unsigned offset1=d_G_dX.offset( l, k, i);

             //Loop over the nodes of the element again
             for(unsigned ll=0;ll<n_node;ll++)
              {
               //Loop of types of dofs again
               for(unsigned kk=0;kk<n_position_type;kk++)
                {
                 //Loop over the displacement components again
                 for(unsigned ii=0;ii<DIM;ii++)
                  {
                   //Get the number of the unknown
                   int local_unknown = this->position_local_eqn(ll,kk,ii);
                   
                   /*IF it's not a boundary condition*/
                   if(local_unknown >= 0)
                    {                               
                     // Offset for faster access in general stress loop
                     const unsigned offset2=d_G_dX.offset( ll, kk, ii);
                     const unsigned offset4=dpsidxi.offset(ll, kk);

                     // General stress term
                     //--------------------
                     double sum=0.0;
                     unsigned count1=offset1;
                     for(unsigned a=0;a<DIM;a++)
                      {
                       // Bump up direct access because we're only
                       // accessing upper triangle
                       count1+=a;
                       for(unsigned b=a;b<DIM;b++)
                        {
                         double factor=d_G_dX.raw_direct_access(count1);
                         if (a==b) factor*=0.5;

                         // Offset for faster access
                         unsigned offset3=d_stress_dG.offset(a,b);
                         unsigned count2=offset2;
                         unsigned count3=offset3;
                         
                         for(unsigned aa=0;aa<DIM;aa++)
                          {
                           // Bump up direct access because we're only
                           // accessing upper triangle
                           count2+=aa;
                           count3+=aa;

                           // Only upper half of derivatives w.r.t. symm tensor
                           for(unsigned bb=aa;bb<DIM;bb++)
                            {                             
                             sum+=factor*
                              d_stress_dG.raw_direct_access(count3)*
                              d_G_dX.raw_direct_access(count2); 
                             ++count2;
                             ++count3;
                            }
                          }
                         ++count1;
                        }

                      }

                     // Multiply by weight and add contribution
                     // (Add directly because this bit is nonsymmetric)
                     jacobian(local_eqn,local_unknown)+=sum*W;
                                        
                     // Only upper triangle (no separate test for bc as
                     // local_eqn is already nonnegative)
                     if (local_unknown >= local_eqn)
                      {
                       double sum=0.0;

                       // Add diagonal terms
                       if (i==ii)
                        {
                         // Inertia term
                         sum+=lambda_sq*time_factor*psi(ll,kk)*psi(l,k);
                         
                         // Stress term                       
                         unsigned count4=offset4;
                         for(unsigned a=0;a<DIM;a++)
                          {
                           //Cache term
                           const double factor=
                            dpsidxi.raw_direct_access(count4);// ll ,kk 
                           ++count4;
                           
                           unsigned count5=offset5;
                           for(unsigned b=0;b<DIM;b++) 
                            {
                             sum+=sigma(a,b)*factor*
                              dpsidxi.raw_direct_access(count5); // l  ,k
                             ++count5;
                            }
                          }
                        }
                       
                                              
                       // Multiply by weight and add contribution
                       sym_jacobian(local_eqn,local_unknown)+=sum*W;
                       
                      } // endif for upper half
    
                    } //End of if not boundary condition
                  }
                }
              }
            }
          
          } //End of if not boundary condition
             
        } //End of loop over coordinate directions
      } //End of loop over type of dof
    } //End of loop over shape functions
  } //End of loop over integration points
        
 // Fill in symmetric lower half
 if (flag==1)
  {
   unsigned n_eqn=temp_eqn_numbers.size();
   for (unsigned i=0;i<n_eqn;i++)
    {
     unsigned local_eqn=temp_eqn_numbers[i];
     for (unsigned j=0;j<n_eqn;j++)
      {
       unsigned local_unknown=temp_eqn_numbers[j];
       
       // Lower bit
       if (local_unknown < local_eqn)
        {
         jacobian(local_eqn,local_unknown)+=
          sym_jacobian(local_unknown,local_eqn);
        }
       // Upper bit and diagonal
       else
        {
         jacobian(local_eqn,local_unknown)+=
          sym_jacobian(local_eqn,local_unknown);
        }
      }
    }
  }
}


//=======================================================================
/// Output: x,y,[z],xi0,xi1,[xi2],gamma
//=======================================================================
template <unsigned DIM>
void PVDEquations<DIM>::output(std::ostream &outfile, const unsigned &n_plot)
{
 
 Vector<double> x(DIM);
 Vector<double> xi(DIM);
 Vector<double> s(DIM);
 
 // Tecplot header info
 outfile << this->tecplot_zone_string(n_plot);
 
 // Loop over plot points
 unsigned num_plot_points=this->nplot_points(n_plot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local coordinates of plot point
   this->get_s_plot(iplot,n_plot,s);
   
   // Get Eulerian and Lagrangian coordinates
   this->interpolated_x(s,x);
   this->interpolated_xi(s,xi);
   
   // Get isotropic growth
   double gamma;
   this->get_isotropic_growth(s,xi,gamma);
   
   //Output the x,y,..
   for(unsigned i=0;i<DIM;i++) 
    {outfile << x[i] << " ";}
   
   // Output xi0,xi1,..
   for(unsigned i=0;i<DIM;i++) 
    {outfile << xi[i] << " ";} 
   
   // Output growth
   outfile << gamma << " ";
   outfile << std::endl;
  }
 

 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(outfile,n_plot);
 outfile << std::endl;
}




//=======================================================================
/// C-style output: x,y,[z],xi0,xi1,[xi2],gamma
//=======================================================================
template <unsigned DIM>
void PVDEquations<DIM>::output(FILE* file_pt, const unsigned &n_plot)
  {
   //Set output Vector
   Vector<double> s(DIM);
   Vector<double> x(DIM);
   Vector<double> xi(DIM);

   switch(DIM)
    {

    case 2:

     //Tecplot header info 
     //outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
     fprintf(file_pt,"ZONE I=%i, J=%i\n",n_plot,n_plot);

     //Loop over element nodes
     for(unsigned l2=0;l2<n_plot;l2++)
      {
       s[1] = -1.0 + l2*2.0/(n_plot-1);
       for(unsigned l1=0;l1<n_plot;l1++)
        {
         s[0] = -1.0 + l1*2.0/(n_plot-1);
         
         // Get Eulerian and Lagrangian coordinates
         this->interpolated_x(s,x);
         this->interpolated_xi(s,xi);

         // Get isotropic growth
         double gamma;
         this->get_isotropic_growth(s,xi,gamma);

         //Output the x,y,..
         for(unsigned i=0;i<DIM;i++) 
          {
           //outfile << x[i] << " ";
           fprintf(file_pt,"%g ",x[i]);
          }
         // Output xi0,xi1,..
         for(unsigned i=0;i<DIM;i++) 
          {
           //outfile << xi[i] << " ";
           fprintf(file_pt,"%g ",xi[i]);
          } 
         // Output growth
         //outfile << gamma << " ";
         //outfile << std::endl;
         fprintf(file_pt,"%g \n",gamma);
        }
      }
     //outfile << std::endl;
     fprintf(file_pt,"\n");

     break;
     
    case 3:

     //Tecplot header info 
     //outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
     fprintf(file_pt,"ZONE I=%i, J=%i, K=%i \n",n_plot,n_plot,n_plot);
     
     //Loop over element nodes
     for(unsigned l3=0;l3<n_plot;l3++)
      {
       s[2] = -1.0 + l3*2.0/(n_plot-1);
       for(unsigned l2=0;l2<n_plot;l2++)
        {
         s[1] = -1.0 + l2*2.0/(n_plot-1);
         for(unsigned l1=0;l1<n_plot;l1++)
          {
           s[0] = -1.0 + l1*2.0/(n_plot-1);
           
           // Get Eulerian and Lagrangian coordinates
           this->interpolated_x(s,x);
           this->interpolated_xi(s,xi);
           
           // Get isotropic growth
           double gamma;
           this->get_isotropic_growth(s,xi,gamma);
           
           //Output the x,y,z
           for(unsigned i=0;i<DIM;i++) 
            {
             //outfile << x[i] << " ";
             fprintf(file_pt,"%g ",x[i]);
            }
           // Output xi0,xi1,xi2
           for(unsigned i=0;i<DIM;i++) 
            {
             //outfile << xi[i] << " ";
             fprintf(file_pt,"%g ",xi[i]);
            } 
           // Output growth
           //outfile << gamma << " ";
           //outfile << std::endl;
           fprintf(file_pt,"%g \n",gamma);
          }
        }
      }
     //outfile << std::endl;
     fprintf(file_pt,"\n");

     break;
     
    default:
     std::ostringstream error_message;
     error_message << "No output routine for PVDEquations<" << 
      DIM << "> elements --  write it yourself!" << std::endl;
     throw OomphLibError(error_message.str(),"PVDEquations<DIM>::output()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }


//=======================================================================
/// Compute the contravariant second Piola Kirchoff stress at a given local
/// coordinate. Note: this replicates a lot of code that is already
/// coontained in get_residuals() but without sacrificing efficiency
/// (re-computing the shape functions several times) or creating
/// helper functions with horrendous interfaces (to pass all the
/// functions which shouldn't be recomputed) about this is
/// unavoidable.
//=======================================================================
template <unsigned DIM>
void PVDEquations<DIM>::get_stress(const Vector<double> &s, 
                                   DenseMatrix<double> &sigma)
{
 //Find out how many nodes there are
 unsigned n_node = this->nnode();

 //Find out how many positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();
 
 //Set up memory for the shape functions
 Shape psi(n_node,n_position_type);
 DShape dpsidxi(n_node,n_position_type,DIM);
 
 //Call the derivatives of the shape functions (ignore Jacobian)
 (void) this->dshape_lagrangian(s,psi,dpsidxi);
 
 // Lagrangian coordinates
 Vector<double> xi(DIM);
 this->interpolated_xi(s,xi);
 
 //Get isotropic growth factor
 double gamma;
 this->get_isotropic_growth(s,xi,gamma);
 
 // We use Cartesian coordinates as the reference coordinate
 // system. In this case the undeformed metric tensor is always
 // the identity matrix -- stretched by the isotropic growth
 double diag_entry=pow(gamma,2.0/double(DIM));
 DenseMatrix<double> g(DIM);
 for(unsigned i=0;i<DIM;i++)
  {
   for(unsigned j=0;j<DIM;j++)
    {
     if(i==j) {g(i,j) = diag_entry;}
     else {g(i,j) = 0.0;}
    }
  }
 
 
 //Calculate interpolated values of the derivative of global position
 //wrt lagrangian coordinates
 DenseMatrix<double> interpolated_G(DIM);
 
 //Initialise to zero
 for(unsigned i=0;i<DIM;i++)
  {for(unsigned j=0;j<DIM;j++) {interpolated_G(i,j) = 0.0;}}
 
 //Calculate displacements and derivatives
 for(unsigned l=0;l<n_node;l++)
  {
   //Loop over positional dofs
   for(unsigned k=0;k<n_position_type;k++)
    {
     //Loop over displacement components (deformed position)
     for(unsigned i=0;i<DIM;i++)
      {
       //Loop over derivative directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_G(j,i) += 
          this->nodal_position_gen(l,k,i)*dpsidxi(l,k,j);
        }
      }
    }
  }
 
 //Declare and calculate the deformed metric tensor
 DenseMatrix<double> G(DIM);
 //Assign values of G
 for(unsigned i=0;i<DIM;i++)
  {
   //Do upper half of matrix
   //Note that j must be signed here for the comparison test to work
   //Also i must be cast to an int
   for(int j=(DIM-1);j>=static_cast<int>(i);j--)
    {
     //Initialise G(i,j) to zero
     G(i,j) = 0.0;
     //Now calculate the dot product
     for(unsigned k=0;k<DIM;k++)
      {
       G(i,j) += interpolated_G(i,k)*interpolated_G(j,k);
      }
    }
   //Matrix is symmetric so just copy lower half
   for(int j=(i-1);j>=0;j--)
    {
     G(i,j) = G(j,i);
    }
  }
 
 //Now calculate the stress tensor from the constitutive law
 get_stress(g,G,sigma);
 
}



//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////



//=======================================================================
/// Compute principal stress vectors and (scalar) principal stresses
/// at specified local coordinate: \c  principal_stress_vector(i,j)
/// is the j-th component of the i-th principal stress vector.
//=======================================================================
template <unsigned DIM>
void PVDEquationsBase<DIM>::get_principal_stress(
 const Vector<double> &s, DenseMatrix<double>& principal_stress_vector,
 Vector<double>& principal_stress)
{

 /// Compute contravariant ("upper") 2nd Piola Kirchhoff stress 
 DenseDoubleMatrix sigma(2,2);
 get_stress(s,sigma);

 // Get covariant base vectors in deformed configuration
 DenseMatrix<double> lower_deformed_basis(DIM);
 get_deformed_covariant_basis_vectors(s,lower_deformed_basis);

 // Work out covariant ("lower") metric tensor
 DenseDoubleMatrix lower_metric(DIM);
 for (unsigned i=0;i<DIM;i++)
  {
   for (unsigned j=0;j<DIM;j++)
    {
     lower_metric(i,j)=0.0;
     for (unsigned k=0;k<DIM;k++)
      {
       lower_metric(i,j)+=
        lower_deformed_basis(i,k)*lower_deformed_basis(j,k);
      }
    }
  }

 // Work out cartesian components of contravariant ("upper") basis vectors
 DenseMatrix<double> upper_deformed_basis(DIM);
 
 // Loop over RHSs
 Vector<double> rhs(DIM);
 Vector<double> aux(DIM);
 for (unsigned k=0;k<DIM;k++)
  {

   for (unsigned l=0;l<DIM;l++)
    {
     rhs[l] = lower_deformed_basis(l,k);
    }

   lower_metric.solve(rhs,aux);

   for (unsigned l=0;l<DIM;l++)
    {     
     upper_deformed_basis(l,k) = aux[l];
    }
  }
 
 // Eigenvalues (=principal stresses) and eigenvectors
 DenseMatrix<double> ev(DIM);
 
 // Get eigenvectors of contravariant 2nd Piola Kirchoff stress
 sigma.eigenvalues_by_jacobi(principal_stress,ev);
 
 // ev(j,i) is the i-th component of the j-th eigenvector
 // relative to the deformed "lower variance" basis! 
 // Work out cartesian components of eigenvectors by multiplying
 // the "lower variance components" by these "upper variance" basis 
 // vectors
 
 // Loop over cartesian compnents
 for(unsigned i=0;i<DIM;i++)
  {
   // Initialise the row
   for(unsigned j=0;j<DIM;j++) {principal_stress_vector(j,i)=0.0;}
   
   // Loop over basis vectors
   for(unsigned j=0;j<DIM;j++)
    {
     for(unsigned k=0;k<DIM;k++)
      {
       principal_stress_vector(j,i) += upper_deformed_basis(k,i)*ev(j,k);
      }
    }
  }
 
 // Scaling factor to turn these vectors into unit vectors
 Vector<double> norm(DIM);
 for (unsigned i=0;i<DIM;i++)
  {
   norm[i]=0.0;
   for (unsigned j=0;j<DIM;j++)
    {
     norm[i] += pow(principal_stress_vector(i,j),2);
    }
   norm[i] = sqrt(norm[i]);
  }
 
 
 // Scaling and then multiplying by eigenvalue gives the principal stress
 // vectors
 for(unsigned i=0;i<DIM;i++)
  {
   for (unsigned j=0;j<DIM;j++)
    {
     principal_stress_vector(j,i) = ev(j,i)/norm[j]*principal_stress[j];
    }
  }
 
}



//=======================================================================
/// Return the deformed covariant basis vectors
/// at specified local coordinate:  \c def_covariant_basis(i,j)
/// is the j-th component of the i-th basis vector.
//=======================================================================
template <unsigned DIM>
void PVDEquationsBase<DIM>::get_deformed_covariant_basis_vectors(
 const Vector<double> &s, DenseMatrix<double>& def_covariant_basis)
{

 //Find out how many nodes there are
 unsigned n_node = nnode();

 //Find out how many positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();
 
 //Set up memory for the shape functions
 Shape psi(n_node,n_position_type);
 DShape dpsidxi(n_node,n_position_type,DIM);


 //Call the derivatives of the shape functions (ignore Jacobian)
 (void) dshape_lagrangian(s,psi,dpsidxi);
 
 
 //Initialise to zero
 for(unsigned i=0;i<DIM;i++)
  {for(unsigned j=0;j<DIM;j++) {def_covariant_basis(i,j) = 0.0;}}
 
 //Calculate displacements and derivatives
 for(unsigned l=0;l<n_node;l++)
  {
   //Loop over positional dofs
   for(unsigned k=0;k<n_position_type;k++)
    {
     //Loop over displacement components (deformed position)
     for(unsigned i=0;i<DIM;i++)
      {
       //Loop over derivative directions (i.e. base vectors)
       for(unsigned j=0;j<DIM;j++)
        {
         def_covariant_basis(j,i) += 
          nodal_position_gen(l,k,i)*dpsidxi(l,k,j);
        }
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//=====================================================================
/// "Magic" number that indicates that the solid pressure is not stored
/// at a node. It is a negative number that cannot be -1 because that is
/// used to represent the positional hanging scheme in Hanging_pt objects
//======================================================================
template<unsigned DIM>
int PVDEquationsWithPressure<DIM>::Solid_pressure_not_stored_at_node = -100;






//=======================================================================
/// Fill in element's contribution to the elemental 
/// residual vector and/or Jacobian matrix.
/// flag=0: compute only residual vector
/// flag=1: compute both, fully analytically
/// flag=2: compute both, using FD for the derivatives w.r.t. to the
///         discrete displacment dofs.
//=======================================================================
template <unsigned DIM>
void PVDEquationsWithPressure<DIM>::
fill_in_generic_residual_contribution_pvd_with_pressure(
 Vector<double> &residuals,DenseMatrix<double> &jacobian, const unsigned& flag)
{

#ifdef PARANOID
 // Check if the constitutive equation requires the explicit imposition of an
 // incompressibility constraint
 if (this->Constitutive_law_pt->requires_incompressibility_constraint()&&
     (!Incompressible))
  {
   throw OomphLibError(
    "The constitutive law requires the use of the incompressible formulation by setting the element's member function incompressible()",
    "PVDEquationsWithPressure<DIM>::fill_in_generic_contribution_to_residuals_pvd_with_pressure()",
    OOMPH_EXCEPTION_LOCATION);   
  }
#endif


 // Simply set up initial condition?
 if (this->Solid_ic_pt!=0)
  {
   this->get_residuals_for_solid_ic(residuals);
   return;
  }
 
 // Storage for symmetric part of the Jacobian
 unsigned n_dof=jacobian.nrow();
 DenseMatrix<double> sym_jacobian(n_dof,n_dof,0.0);
 
 // Storage of eqn numbers
 Vector<unsigned> temp_eqn_numbers;
 temp_eqn_numbers.reserve(n_dof);
 
 //Find out how many nodes there are
 unsigned n_node = this->nnode();

 //Find out how many position types of dof there are
 unsigned n_position_type = this->nnodal_position_type();

 //Find out how many pressure dofs there are
 unsigned n_solid_pres = npres_solid();

 //Set up memory for the shape functions
 Shape psi(n_node,n_position_type);
 DShape dpsidxi(n_node,n_position_type,DIM);
 
 //Set up memory for the pressure shape functions
 Shape psisp(n_solid_pres);

 //Set the value of n_intpt
 unsigned n_intpt = this->integral_pt()->nweight();

 //Set the vector to hold the local coordinates in the element
 Vector<double> s(DIM);

 // Timescale ratio (non-dim density)
 double lambda_sq = this->lambda_sq();

 // Time factor
 double time_factor=0.0;
 if (lambda_sq>0)
  {
   time_factor=this->node_pt(0)->position_time_stepper_pt()->weight(2,0);
  }

 //Integers to hold the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Assign the values of s
   for(unsigned i=0;i<DIM;++i)
    {
     s[i] = this->integral_pt()->knot(ipt,i);
    }
   
   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);

   //Call the derivatives of the shape functions
   double J = this->dshape_lagrangian_at_knot(ipt,psi,dpsidxi);

   //Call the pressure shape functions
   solid_pshape_at_knot(ipt,psisp);

   //Storage for Lagrangian coordinates (initialised to zero)
   Vector<double> interpolated_xi(DIM,0.0);

   // Deformed tangent vectors
   DenseMatrix<double> interpolated_G(DIM);

   // Setup memory for accelerations
   Vector<double> accel(DIM);
      
   //Initialise to zero
   for(unsigned i=0;i<DIM;i++)
    {
     // Initialise acclerations
     accel[i]=0.0;
     for(unsigned j=0;j<DIM;j++)
      {
       interpolated_G(i,j) = 0.0;
      }
    }

   //Calculate displacements and derivatives and lagrangian coordinates
   for(unsigned l=0;l<n_node;l++)
    {
     //Loop over positional dofs
     for(unsigned k=0;k<n_position_type;k++)
      {
       //Loop over displacement components (deformed position)
       for(unsigned i=0;i<DIM;i++)
        {
         //Calculate the lagrangian coordinates and the accelerations
         interpolated_xi[i] += this->lagrangian_position_gen(l,k,i)*psi(l,k);

         // Only compute accelerations if inertia is switched on
         // otherwise the timestepper might not be able to 
         // work out dx_gen_dt(2,...)
         if ((lambda_sq>0.0)&&(this->unsteady()))
          {
           accel[i] += this->dnodal_position_gen_dt(2,l,k,i)*psi(l,k);
          }
         
         //Loop over derivative directions
         for(unsigned j=0;j<DIM;j++)
          {
           interpolated_G(j,i) += 
            this->nodal_position_gen(l,k,i)*dpsidxi(l,k,j);
          }
        }
      }
    }

   //Get isotropic growth factor
   double gamma=1.0;
   this->get_isotropic_growth(s,interpolated_xi,gamma);

   //Get body force at current time
   Vector<double> b(DIM);
   this->body_force(interpolated_xi,b);

   // We use Cartesian coordinates as the reference coordinate
   // system. In this case the undeformed metric tensor is always
   // the identity matrix -- stretched by the isotropic growth
   double diag_entry=pow(gamma,2.0/double(DIM)); 
   DenseMatrix<double> g(DIM);
   for(unsigned i=0;i<DIM;i++)
    {
     for(unsigned j=0;j<DIM;j++)
      {
       if(i==j) {g(i,j) = diag_entry;}
       else {g(i,j) = 0.0;}
      }
    }
     
   //Premultiply the undeformed volume ratio (from the isotropic
   // growth), the weights and the Jacobian
   double W = gamma*w*J; 
     
   //Calculate the interpolated solid pressure
   double interpolated_solid_p=0.0;
   for(unsigned l=0;l<n_solid_pres;l++)
    {
     interpolated_solid_p += solid_p(l)*psisp[l];
    }


   //Declare and calculate the deformed metric tensor
   DenseMatrix<double> G(DIM);

   //Assign values of G
   for(unsigned i=0;i<DIM;i++)
    {
     //Do upper half of matrix
     for(unsigned j=i;j<DIM;j++) 
      {
       //Initialise G(i,j) to zero
       G(i,j) = 0.0;
       //Now calculate the dot product
       for(unsigned k=0;k<DIM;k++)
        {
         G(i,j) += interpolated_G(i,k)*interpolated_G(j,k);
        }
      }
     //Matrix is symmetric so just copy lower half
     for(unsigned j=0;j<i;j++) 
      {
       G(i,j) = G(j,i);
      }
    }

   //Now calculate the deviatoric stress and all pressure-related
   //quantitites
   DenseMatrix<double> sigma(DIM,DIM), sigma_dev(DIM,DIM), Gup(DIM,DIM);
   double detG = 0.0;
   double gen_dil=0.0;
   double inv_kappa=0.0;

   // Get stress derivative by FD only needed for Jacobian
      
   // Stress etc derivatives
   RankFourTensor<double> d_stress_dG(DIM,DIM,DIM,DIM,0.0);
   RankFourTensor<double> d_Gup_dG(DIM,DIM,DIM,DIM,0.0);
   DenseMatrix<double> d_detG_dG(DIM,DIM,0.0);
   DenseMatrix<double> d_gen_dil_dG(DIM,DIM,0.0);

   // Derivative of metric tensor w.r.t. to nodal coords
   RankFiveTensor<double> d_G_dX(n_node,n_position_type,DIM,DIM,DIM,0.0);

   // Get Jacobian too?
   if (flag==1) 
    {     
     // Derivative of metric tensor w.r.t. to discrete positional dofs
     // NOTE: Since G is symmetric we only compute the upper triangle
     //       and DO NOT copy the entries across. Subsequent computations
     //       must (and, in fact, do) therefore only operate with upper
     //       triangular entries
     for (unsigned ll=0;ll<n_node;ll++)
      {
       for (unsigned kk=0;kk<n_position_type;kk++)
        {
         for (unsigned ii=0;ii<DIM;ii++)
          {
           for (unsigned aa=0;aa<DIM;aa++)
            {
             for (unsigned bb=aa;bb<DIM;bb++) 
              {
               d_G_dX(ll,kk,ii,aa,bb)=
                interpolated_G(aa,ii)*dpsidxi(ll,kk,bb)+
                interpolated_G(bb,ii)*dpsidxi(ll,kk,aa);
              }
            }
          }
        }
      }
    }


   // Incompressible: Compute the deviatoric part of the stress tensor, the
   // contravariant deformed metric tensor and the determinant
   // of the deformed covariant metric tensor.
   if(Incompressible)
    {
     get_stress(g,G,sigma_dev,Gup,detG);
     
     // Get full stress
     for (unsigned a=0;a<DIM;a++)
      {
       for (unsigned b=0;b<DIM;b++)
        {         
         sigma(a,b)=sigma_dev(a,b) - interpolated_solid_p*Gup(a,b);
        }
      }

     // Get Jacobian too?
     if (flag==1)
      {
       // FD step
       double eps_fd=GeneralisedElement::Default_fd_jacobian_step;
       
       //Advanced metric tensor etc.
       DenseMatrix<double> G_pls(DIM,DIM);
       DenseMatrix<double> sigma_dev_pls(DIM,DIM);
       DenseMatrix<double> sigma_pls(DIM,DIM);
       DenseMatrix<double> Gup_pls(DIM,DIM);
       double detG_pls;

       // Copy across
       for (unsigned i=0;i<DIM;i++)
        {
         for (unsigned j=0;j<DIM;j++)
          {
           G_pls(i,j)=G(i,j);
          }
        }
       

       // Do FD -- only w.r.t. to upper indices, exploiting symmetry.
       // NOTE: We exploit the symmetry of the stress and metric tensors
       //       by incrementing G(i,j) and G(j,i) simultaenously and
       //       only fill in the "upper" triangles without copying things
       //       across the lower triangle. This is taken into account
       //       in the remaining code further below.
       for(unsigned i=0;i<DIM;i++)
        {
         for (unsigned j=i;j<DIM;j++)
          {
           G_pls(i,j) += eps_fd;
           G_pls(j,i) = G_pls(i,j);

           // Get advanced stress
           get_stress(g,G_pls,sigma_dev_pls,Gup_pls,detG_pls);
           
           // Derivative of determinant of deformed metric tensor
           d_detG_dG(i,j)=(detG_pls-detG)/eps_fd;

           // Derivatives of deviatoric stress and "upper" deformed metric
           // tensor
           for (unsigned ii=0;ii<DIM;ii++)
            {
             for (unsigned jj=ii;jj<DIM;jj++)
              {
               d_stress_dG(ii,jj,i,j)=(
                sigma_dev_pls(ii,jj)-interpolated_solid_p*Gup_pls(ii,jj)-
                sigma(ii,jj))/eps_fd;

               d_Gup_dG(ii,jj,i,j)=(Gup_pls(ii,jj)-Gup(ii,jj))/eps_fd;
              }  
            }
           
           // Reset 
           G_pls(i,j) = G(i,j);
           G_pls(j,i) = G(j,i);

          }
         
        }
      }
    }
   // Nearly incompressible: Compute the deviatoric part of the 
   // stress tensor, the contravariant deformed metric tensor,
   // the generalised dilatation and the inverse bulk modulus.
   else
    {
     get_stress(g,G,sigma_dev,Gup,gen_dil,inv_kappa);

     // Get full stress
     for (unsigned a=0;a<DIM;a++)
      {
       for (unsigned b=0;b<DIM;b++)
        {         
         sigma(a,b)=sigma_dev(a,b) - interpolated_solid_p*Gup(a,b);
        }
      }

     // Get Jacobian too?
     if (flag==1) 
      {
       // FD step 
       double eps_fd=GeneralisedElement::Default_fd_jacobian_step;
       
       //Advanced metric tensor etc
       DenseMatrix<double> G_pls(DIM,DIM);
       DenseMatrix<double> sigma_dev_pls(DIM,DIM);
       DenseMatrix<double> Gup_pls(DIM,DIM);
       double gen_dil_pls;

       // Copy across
       for (unsigned i=0;i<DIM;i++)
        {
         for (unsigned j=0;j<DIM;j++)
          {
           G_pls(i,j)=G(i,j);
          }
        }

       // Do FD -- only w.r.t. to upper indices, exploiting symmetry.
       // NOTE: We exploit the symmetry of the stress and metric tensors
       //       by incrementing G(i,j) and G(j,i) simultaenously and
       //       only fill in the "upper" triangles without copying things
       //       across the lower triangle. This is taken into account
       //       in the remaining code further below.
       for(unsigned i=0;i<DIM;i++)
        {
         for (unsigned j=i;j<DIM;j++)
          {
           G_pls(i,j) += eps_fd;
           G_pls(j,i) = G_pls(i,j);

           // Get advanced stress
           get_stress(g,G_pls,sigma_dev_pls,Gup_pls,gen_dil_pls,inv_kappa);

           // Derivative of generalised dilatation
           d_gen_dil_dG(i,j)=(gen_dil_pls-gen_dil)/eps_fd;

           // Derivatives of deviatoric stress and "upper" deformed metric
           // tensor
           for (unsigned ii=0;ii<DIM;ii++)
            {
             for (unsigned jj=ii;jj<DIM;jj++)
              {
               d_stress_dG(ii,jj,i,j)=(
                sigma_dev_pls(ii,jj)-interpolated_solid_p*Gup_pls(ii,jj)-
                sigma(ii,jj))/eps_fd;

               d_Gup_dG(ii,jj,i,j)=(Gup_pls(ii,jj)-Gup(ii,jj))/eps_fd;
              }        
            }
           
           // Reset 
           G_pls(i,j) = G(i,j); 
           G_pls(j,i) = G(j,i); 
          }
        }
      }
    }


//=====EQUATIONS OF ELASTICITY FROM PRINCIPLE OF VIRTUAL DISPLACEMENTS========
       
   //Loop over the test functions, nodes of the element
   for(unsigned l=0;l<n_node;l++)
    {
     //Loop over the types of dof
     for(unsigned k=0;k<n_position_type;k++)
      {
       // Offset for faster access
       const unsigned offset5=dpsidxi.offset(l ,k);

       //Loop over the displacement components
       for(unsigned i=0;i<DIM;i++)
        {
         //Get the equation number
         local_eqn = this->position_local_eqn(l,k,i);

         /*IF it's not a boundary condition*/
         if(local_eqn >= 0)
          {
           // Record eqn number
           if (ipt==0) temp_eqn_numbers.push_back(local_eqn);
           
           double sum=0.0;
           
           // Acceleration and body force
           sum+=(lambda_sq*accel[i]-b[i])*psi(l,k);
           
           // Stress term
           for(unsigned a=0;a<DIM;a++)
            {
             unsigned count=offset5;
             for(unsigned b=0;b<DIM;b++) 
              {
               //Add the stress terms to the residuals
               sum+=sigma(a,b)*interpolated_G(a,i)*
                dpsidxi.raw_direct_access(count);
               ++count;
              }
            }
           residuals[local_eqn] += W*sum;
                      
           //Add in the jacobian terms
           if(flag==1)
            {
             // Offset for faster access in general stress loop
             const unsigned offset1=d_G_dX.offset( l, k, i);

             //Loop over the nodes of the element again
             for(unsigned ll=0;ll<n_node;ll++)
              {
               //Loop of types of dofs again
               for(unsigned kk=0;kk<n_position_type;kk++)
                {
                 //Loop over the displacement components again
                 for(unsigned ii=0;ii<DIM;ii++)
                  {
                   //Get the number of the unknown
                   int local_unknown = this->position_local_eqn(ll,kk,ii);
                   
                   /*IF it's not a boundary condition*/
                   if(local_unknown >= 0)
                    {
                    
                     // Offset for faster access in general stress loop
                     const unsigned offset2=d_G_dX.offset( ll, kk, ii);
                     const unsigned offset4=dpsidxi.offset(ll, kk);
                     
                     // General stress term
                     //--------------------
                     double sum=0.0;
                     unsigned count1=offset1;
                     for(unsigned a=0;a<DIM;a++)
                      {
                       // Bump up direct access because we're only
                       // accessing upper triangle
                       count1+=a;
                       for(unsigned b=a;b<DIM;b++)
                        {
                         double factor=d_G_dX.raw_direct_access(count1);
                         if (a==b) factor*=0.5;
                         
                         // Offset for faster access
                         unsigned offset3=d_stress_dG.offset(a,b);
                         unsigned count2=offset2;
                         unsigned count3=offset3;
                         
                         for(unsigned aa=0;aa<DIM;aa++)
                          {
                           // Bump up direct access because we're only
                           // accessing upper triangle
                           count2+=aa;
                           count3+=aa;
                           
                           // Only upper half of derivatives w.r.t. symm tensor
                           for(unsigned bb=aa;bb<DIM;bb++)
                            {                             
                             sum+=factor*
                              d_stress_dG.raw_direct_access(count3)*
                              d_G_dX.raw_direct_access(count2); 
                             ++count2;
                             ++count3;
                            }
                          }
                         ++count1;
                        }
                       
                      }
                     
                     // Multiply by weight and add contribution
                     // (Add directly because this bit is nonsymmetric)
                     jacobian(local_eqn,local_unknown)+=sum*W;
                                          
                     // Only upper triangle (no separate test for bc as
                     // local_eqn is already nonnegative)
                     if (local_unknown >= local_eqn)
                      {
                       double sum=0.0;
                       
                       // Add diagonal terms
                       if (i==ii)
                        {
                         // Inertia term
                         sum+=lambda_sq*time_factor*psi(ll,kk)*psi(l,k);
                         
                         // Stress term                       
                         unsigned count4=offset4;
                         for(unsigned a=0;a<DIM;a++)
                          {
                           //Cache term
                           const double factor=
                            dpsidxi.raw_direct_access(count4);// ll ,kk 
                           ++count4;
                           
                           unsigned count5=offset5;
                           for(unsigned b=0;b<DIM;b++) 
                            {
                             sum+=sigma(a,b)*factor*
                              dpsidxi.raw_direct_access(count5); // l  ,k
                             ++count5;
                            }
                          }
                        }
                       
                       // Multiply by weight and add contribution
                       sym_jacobian(local_eqn,local_unknown)+=sum*W;
                       
                      } // endif for upper half
    
                    } //End of if not boundary condition
                  }
                }
              }
            }

           // Derivatives w.r.t. pressure dofs
           if (flag>0) 
            {
             //Loop over the pressure dofs for unknowns
             for(unsigned l2=0;l2<n_solid_pres;l2++)
              {
               local_unknown = solid_p_local_eqn(l2);

               //If it's not a boundary condition
               if(local_unknown >= 0)
                {

                 //Add the pressure terms to the jacobian
                 for(unsigned a=0;a<DIM;a++)
                  {
                   for(unsigned b=0;b<DIM;b++)
                    {
                     jacobian(local_eqn,local_unknown) -=
                      psisp[l2]*Gup(a,b)*
                      interpolated_G(a,i)*dpsidxi(l,k,b)*W; 
                    }
                  }
                }
              }
            } //End of Jacobian terms
           
          } //End of if not boundary condition
        } //End of loop over coordinate directions
      } //End of loop over types of dof
    } //End of loop over shape functions

   //==============CONSTRAINT EQUATIONS FOR PRESSURE=====================
       
   //Now loop over the pressure degrees of freedom
   for(unsigned l=0;l<n_solid_pres;l++)
    {
     local_eqn = solid_p_local_eqn(l);

     // Pinned (unlikely, actually) or real dof?
     if(local_eqn >= 0)
      {
       // Record eqn number
       if (ipt==0) temp_eqn_numbers.push_back(local_eqn);

       //For true incompressibility we need to conserve volume
       //so the determinant of the deformed metric tensor
       //needs to be equal to that of the undeformed one, which
       //is equal to the volumetric growth factor
       if(Incompressible)
        {
         residuals[local_eqn] += (detG - gamma)*psisp[l]*W;
             

         // Get Jacobian too?
         if (flag==1) 
          {
           //Loop over the nodes of the element again
           for(unsigned ll=0;ll<n_node;ll++)
            {
             //Loop of types of dofs again
             for(unsigned kk=0;kk<n_position_type;kk++)
              {
               //Loop over the displacement components again
               for(unsigned ii=0;ii<DIM;ii++)
                {
                 //Get the number of the unknown
                 int local_unknown = this->position_local_eqn(ll,kk,ii);
                 
                 /*IF it's not a boundary condition*/
                 if(local_unknown >= 0)
                  {
          
                   // Offset for faster access
                   const unsigned offset=d_G_dX.offset( ll, kk, ii);

                   // General stress term
                   double sum=0.0;
                   unsigned count=offset;
                   for(unsigned aa=0;aa<DIM;aa++)
                    {
                     // Bump up direct access because we're only
                     // accessing upper triangle
                     count+=aa;

                     // Only upper half
                     for(unsigned bb=aa;bb<DIM;bb++)
                      {          
                       sum+=d_detG_dG(aa,bb)*
                        d_G_dX.raw_direct_access(count)*psisp(l);
                       ++count;
                      }
                    }
                   jacobian(local_eqn,local_unknown)+=sum*W;
                  }
                }
              }
            }
          
           //No Jacobian terms due to pressure since it does not feature
           //in the incompressibility constraint
          }
        }
       //Nearly incompressible: (Neg.) pressure given by product of
       //bulk modulus and generalised dilatation
       else
        {
         residuals[local_eqn] += 
          (inv_kappa*interpolated_solid_p + gen_dil)*psisp[l]*W;
             
         //Add in the jacobian terms
         if (flag==1)
          {
           //Loop over the nodes of the element again
           for(unsigned ll=0;ll<n_node;ll++)
            {
             //Loop of types of dofs again
             for(unsigned kk=0;kk<n_position_type;kk++)
              {
               //Loop over the displacement components again
               for(unsigned ii=0;ii<DIM;ii++)
                {
                 //Get the number of the unknown
                 int local_unknown = this->position_local_eqn(ll,kk,ii);
                 
                 /*IF it's not a boundary condition*/
                 if(local_unknown >= 0)
                  {
                   // Offset for faster access
                   const unsigned offset=d_G_dX.offset( ll, kk, ii);

                   // General stress term
                   double sum=0.0;
                   unsigned count=offset;
                   for(unsigned aa=0;aa<DIM;aa++)
                    {
                     // Bump up direct access because we're only
                     // accessing upper triangle
                     count+=aa;

                     // Only upper half
                     for(unsigned bb=aa;bb<DIM;bb++)
                      {
                       sum+=d_gen_dil_dG(aa,bb)*
                        d_G_dX.raw_direct_access(count)*psisp(l);          
                       ++count;
                      }
                    }
                   jacobian(local_eqn,local_unknown)+=sum*W;
                  }
                }
              }
            }
          }

         // Derivatives w.r.t. pressure dofs
         if(flag>0) 
          {
           //Loop over the pressure nodes again
           for(unsigned l2=0;l2<n_solid_pres;l2++)
            {
             local_unknown = solid_p_local_eqn(l2);
             //If not pinnned 
             if(local_unknown >= 0)
              {
               jacobian(local_eqn,local_unknown)
                += inv_kappa*psisp[l2]*psisp[l]*W;
              }
            }
          } //End of jacobian terms

        } //End of else
       
      } //End of if not boundary condition
    }
   
  } //End of loop over integration points


 // Fill in symmetric lower half
 if (flag==1)
  {
   unsigned n_eqn=temp_eqn_numbers.size();
   for (unsigned i=0;i<n_eqn;i++)
    {
     unsigned local_eqn=temp_eqn_numbers[i];
     for (unsigned j=0;j<n_eqn;j++)
      {
       unsigned local_unknown=temp_eqn_numbers[j];
       
       // Lower bit
       if (local_unknown < local_eqn)
        {
         jacobian(local_eqn,local_unknown)+=
          sym_jacobian(local_unknown,local_eqn);
        }
       // Upper bit and diagonal
       else
        {
         jacobian(local_eqn,local_unknown)+=
          sym_jacobian(local_eqn,local_unknown);
        }
      }
    }
  }
}



//=======================================================================
/// Output: x,y,[z],xi0,xi1,[xi2],p,gamma
//=======================================================================
template <unsigned DIM>
void PVDEquationsWithPressure<DIM>::output(std::ostream &outfile, 
                                           const unsigned &n_plot)
{
 //Set output Vector
 Vector<double> s(DIM);
 Vector<double> x(DIM);
 Vector<double> xi(DIM);

 switch(DIM)
  {
  case 2:
   //Tecplot header info 
   outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
   
   //Loop over element nodes
   for(unsigned l2=0;l2<n_plot;l2++)
    {
     s[1] = -1.0 + l2*2.0/(n_plot-1);
     for(unsigned l1=0;l1<n_plot;l1++)
      {
       s[0] = -1.0 + l1*2.0/(n_plot-1);
       
       // Get Eulerian and Lagrangian coordinates
       this->interpolated_x(s,x);
       this->interpolated_xi(s,xi);
       
       // Get isotropic growth
       double gamma;
       this->get_isotropic_growth(s,xi,gamma);
       
       //Output the x,y,..
       for(unsigned i=0;i<DIM;i++) 
        {outfile << x[i] << " ";}
       // Output xi0,xi1,..
       for(unsigned i=0;i<DIM;i++) 
        {outfile << xi[i] << " ";} 
       // Output growth
       outfile << gamma << " ";
       // Output pressure
       outfile << interpolated_solid_p(s) << " ";
       outfile << std::endl;
      }
    }

   break;
   
  case 3:
   //Tecplot header info 
   outfile << "ZONE I=" << n_plot 
           << ", J=" << n_plot 
           << ", K=" << n_plot << std::endl;
   
     //Loop over element nodes
     for(unsigned l3=0;l3<n_plot;l3++)
      {
       s[2] = -1.0 + l3*2.0/(n_plot-1);
       for(unsigned l2=0;l2<n_plot;l2++)
        {
         s[1] = -1.0 + l2*2.0/(n_plot-1);
         for(unsigned l1=0;l1<n_plot;l1++)
          {
           s[0] = -1.0 + l1*2.0/(n_plot-1);
           
           // Get Eulerian and Lagrangian coordinates
           this->interpolated_x(s,x);
           this->interpolated_xi(s,xi);
           
           // Get isotropic growth
           double gamma;
           this->get_isotropic_growth(s,xi,gamma);
           
           //Output the x,y,..
           for(unsigned i=0;i<DIM;i++) 
            {outfile << x[i] << " ";}
           // Output xi0,xi1,..
           for(unsigned i=0;i<DIM;i++) 
            {outfile << xi[i] << " ";} 
           // Output growth
           outfile << gamma << " ";
           // Output pressure
           outfile << interpolated_solid_p(s) << " ";
           outfile << std::endl;           
          }
        }
      }
     break;
     
    default:
     std::ostringstream error_message;
     error_message << "No output routine for PVDEquationsWithPressure<" << 
      DIM << "> elements. Write it yourself!" << std::endl;
     throw OomphLibError(error_message.str(),
                         "PVDEquationsWithPressure<DIM>::output()",
                         OOMPH_EXCEPTION_LOCATION);
  }
}




//=======================================================================
/// C-stsyle output: x,y,[z],xi0,xi1,[xi2],p,gamma
//=======================================================================
template <unsigned DIM>
void PVDEquationsWithPressure<DIM>::output(FILE* file_pt,
                                           const unsigned &n_plot)
{
 //Set output Vector
 Vector<double> s(DIM);
 Vector<double> x(DIM);
 Vector<double> xi(DIM);

 switch(DIM)
  {
  case 2:
   //Tecplot header info 
   //outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
   fprintf(file_pt,"ZONE I=%i, J=%i\n",n_plot,n_plot);

   //Loop over element nodes
   for(unsigned l2=0;l2<n_plot;l2++)
    {
     s[1] = -1.0 + l2*2.0/(n_plot-1);
     for(unsigned l1=0;l1<n_plot;l1++)
      {
       s[0] = -1.0 + l1*2.0/(n_plot-1);
       
       // Get Eulerian and Lagrangian coordinates
       this->interpolated_x(s,x);
       this->interpolated_xi(s,xi);
       
       // Get isotropic growth
       double gamma;
       this->get_isotropic_growth(s,xi,gamma);
       
       //Output the x,y,..
       for(unsigned i=0;i<DIM;i++) 
        {
         //outfile << x[i] << " ";
         fprintf(file_pt,"%g ",x[i]);
        }
       // Output xi0,xi1,..
       for(unsigned i=0;i<DIM;i++) 
        {
         //outfile << xi[i] << " ";
         fprintf(file_pt,"%g ",xi[i]);
        } 
       // Output growth
       //outfile << gamma << " ";
       fprintf(file_pt,"%g ",gamma);

       // Output pressure
       //outfile << interpolated_solid_p(s) << " ";
       //outfile << std::endl;
       fprintf(file_pt,"%g \n",interpolated_solid_p(s));

      }
    }

   break;
   
  case 3:
   //Tecplot header info 
   //outfile << "ZONE I=" << n_plot 
   //        << ", J=" << n_plot 
   //        << ", K=" << n_plot << std::endl;
   fprintf(file_pt,"ZONE I=%i, J=%i, K=%i \n",n_plot,n_plot,n_plot);

     //Loop over element nodes
     for(unsigned l3=0;l3<n_plot;l3++)
      {
       s[2] = -1.0 + l3*2.0/(n_plot-1);
       for(unsigned l2=0;l2<n_plot;l2++)
        {
         s[1] = -1.0 + l2*2.0/(n_plot-1);
         for(unsigned l1=0;l1<n_plot;l1++)
          {
           s[0] = -1.0 + l1*2.0/(n_plot-1);
           
           // Get Eulerian and Lagrangian coordinates
           this->interpolated_x(s,x);
           this->interpolated_xi(s,xi);
           
           // Get isotropic growth
           double gamma;
           this->get_isotropic_growth(s,xi,gamma);
           
           //Output the x,y,..
           for(unsigned i=0;i<DIM;i++) 
            {
             //outfile << x[i] << " ";
             fprintf(file_pt,"%g ",x[i]);
            }
           // Output xi0,xi1,..
           for(unsigned i=0;i<DIM;i++) 
            {
             //outfile << xi[i] << " ";
             fprintf(file_pt,"%g ",xi[i]);
            } 
           // Output growth
           //outfile << gamma << " ";
           fprintf(file_pt,"%g ",gamma);

           // Output pressure
           //outfile << interpolated_solid_p(s) << " ";
           //outfile << std::endl;           
           fprintf(file_pt,"%g \n",interpolated_solid_p(s));
          }
        }
      }
     break;
     
    default:
     std::ostringstream error_message;
     error_message << "No output routine for PVDEquationsWithPressure<" << 
      DIM << "> elements. Write it yourself!" << std::endl;
     throw OomphLibError(error_message.str(),
                         "PVDEquationsWithPressure<DIM>::output()",
                         OOMPH_EXCEPTION_LOCATION);
  }
}



//=======================================================================
/// Compute the contravariant second Piola Kirchoff stress at a given local
/// coordinate. Note: this replicates a lot of code that is already
/// coontained in get_residuals() but without sacrificing efficiency
/// (re-computing the shape functions several times) or creating
/// helper functions with horrendous interfaces (to pass all the
/// functions which shouldn't be recomputed) about this is
/// unavoidable.
//=======================================================================
template <unsigned DIM>
void PVDEquationsWithPressure<DIM>::get_stress(const Vector<double> &s, 
                                               DenseMatrix<double> &sigma)
{
 //Find out how many nodes there are
 unsigned n_node = this->nnode();

 //Find out how many positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();
 
 //Find out how many pressure dofs there are
 unsigned n_solid_pres = npres_solid();

 //Set up memory for the shape functions
 Shape psi(n_node,n_position_type);
 DShape dpsidxi(n_node,n_position_type,DIM);

 //Set up memory for the pressure shape functions
 Shape psisp(n_solid_pres);

 //Find values of shape function
 solid_pshape(s,psisp);

 //Call the derivatives of the shape functions (ignore Jacobian)
 (void) this->dshape_lagrangian(s,psi,dpsidxi);
 
 // Lagrangian coordinates
 Vector<double> xi(DIM);
 this->interpolated_xi(s,xi);
 
 //Get isotropic growth factor
 double gamma;
 this->get_isotropic_growth(s,xi,gamma);
 
 // We use Cartesian coordinates as the reference coordinate
 // system. In this case the undeformed metric tensor is always
 // the identity matrix -- stretched by the isotropic growth
 double diag_entry=pow(gamma,2.0/double(DIM));
 DenseMatrix<double> g(DIM);
 for(unsigned i=0;i<DIM;i++)
  {
   for(unsigned j=0;j<DIM;j++)
    {
     if(i==j) {g(i,j) = diag_entry;}
     else {g(i,j) = 0.0;}
    }
  }
 
 
 //Calculate interpolated values of the derivative of global position
 //wrt lagrangian coordinates
 DenseMatrix<double> interpolated_G(DIM);
 
 //Initialise to zero
 for(unsigned i=0;i<DIM;i++)
  {for(unsigned j=0;j<DIM;j++) {interpolated_G(i,j) = 0.0;}}
 
 //Calculate displacements and derivatives
 for(unsigned l=0;l<n_node;l++)
  {
   //Loop over positional dofs
   for(unsigned k=0;k<n_position_type;k++)
    {
     //Loop over displacement components (deformed position)
     for(unsigned i=0;i<DIM;i++)
      {
       //Loop over derivative directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_G(j,i) += this->nodal_position_gen(l,k,i)*dpsidxi(l,k,j);
        }
      }
    }
  }
 
 //Declare and calculate the deformed metric tensor
 DenseMatrix<double> G(DIM);

 //Assign values of G
 for(unsigned i=0;i<DIM;i++)
  {
   //Do upper half of matrix
   //Note that j must be signed here for the comparison test to work
   //Also i must be cast to an int
   for(int j=(DIM-1);j>=static_cast<int>(i);j--)
    {
     //Initialise G(i,j) to zero
     G(i,j) = 0.0;
     //Now calculate the dot product
     for(unsigned k=0;k<DIM;k++)
      {
       G(i,j) += interpolated_G(i,k)*interpolated_G(j,k);
      }
    }
   //Matrix is symmetric so just copy lower half
   for(int j=(i-1);j>=0;j--)
    {
     G(i,j) = G(j,i);
    }
  }
 

 //Calculate the interpolated solid pressure
 double interpolated_solid_p=0.0;
 for(unsigned l=0;l<n_solid_pres;l++)
  {
   interpolated_solid_p += solid_p(l)*psisp[l];
  }

 //Now calculate the deviatoric stress and all pressure-related
 //quantitites
 DenseMatrix<double> sigma_dev(DIM), Gup(DIM);
 double detG = 0.0;
 double gen_dil=0.0;
 double inv_kappa=0.0;
 
 // Incompressible: Compute the deviatoric part of the stress tensor, the
 // contravariant deformed metric tensor and the determinant
 // of the deformed covariant metric tensor.
 
 if(Incompressible)
  {
   get_stress(g,G,sigma_dev,Gup,detG);
  }
 // Nearly incompressible: Compute the deviatoric part of the 
 // stress tensor, the contravariant deformed metric tensor,
 // the generalised dilatation and the inverse bulk modulus.
 else
  {
   get_stress(g,G,sigma_dev,Gup,gen_dil,inv_kappa);
  }

 // Get complete stress
 for (unsigned i=0;i<DIM;i++)
  {
   for (unsigned j=0;j<DIM;j++)
    {
     sigma(i,j) = -interpolated_solid_p*Gup(i,j)+sigma_dev(i,j);
    }
  }

}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



//====================================================================
/// Data for the number of Variables at each node
//====================================================================
template<>
const unsigned QPVDElementWithContinuousPressure<2>::Initial_Nvalue[9]=
{1,0,1,0,0,0,1,0,1};

//==========================================================================
/// Conversion from pressure dof to Node number at which pressure is stored
//==========================================================================
template<>
const unsigned QPVDElementWithContinuousPressure<2>::Pconv[4] =
{0,2,6,8};

//====================================================================
/// Data for the number of Variables at each node
//====================================================================
template<>
const unsigned QPVDElementWithContinuousPressure<3>::Initial_Nvalue[27]=
{1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1};

//==========================================================================
/// Conversion from pressure dof to Node number at which pressure is stored
//==========================================================================
template<>
const unsigned QPVDElementWithContinuousPressure<3>::Pconv[8] =
{0,2,6,8,18,20,24,26};


//Instantiate the required elements
template class QPVDElementWithPressure<2>;
template class QPVDElementWithContinuousPressure<2>;
template class PVDEquationsBase<2>;
template class PVDEquations<2>;
template class PVDEquationsWithPressure<2>;

template class QPVDElementWithPressure<3>;
template class QPVDElementWithContinuousPressure<3>;
template class PVDEquationsBase<3>;
template class PVDEquations<3>;
template class PVDEquationsWithPressure<3>;


}
