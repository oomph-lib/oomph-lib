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

//============start_element_class============================================
///A RefineableElement class that solves the 
///Boussinesq approximation of the Navier--Stokes
///and energy equations by coupling two pre-existing classes. 
///The RefineableQAdvectionDiffusionElement 
///with bi-quadratic interpolation for the
///scalar variable (temperature) and
///RefineableQCrouzeixRaviartElement which solves the Navier--Stokes equations
///using bi-quadratic interpolation for the velocities and a discontinuous
///bi-linear interpolation for the pressure. Note that we are free to 
///choose the order in which we store the variables at the nodes. In this
///case we choose to store the variables in the order fluid velocities
///followed by temperature. We must, therefore, overload the function
///AdvectionDiffusionEquations<DIM>::u_index_adv_diff() to indicate that
///the temperature is stored at the DIM-th position not the 0-th. We do not
///need to overload the corresponding function in the 
///NavierStokesEquations<DIM> class because the velocities are stored
///first. Finally, we choose to use the flux-recovery calculation from the
///fluid velocities to provide the error used in the mesh adaptation.
//==========================================================================
template<unsigned DIM>
class RefineableBuoyantQCrouzeixRaviartElement: 
public virtual RefineableQAdvectionDiffusionElement<DIM,3>,
public virtual RefineableQCrouzeixRaviartElement<DIM>
{

private:

 ///Pointer to a new physical variable, the Rayleigh number
 double* Ra_pt;

 /// The static default value of the Rayleigh number
 static double Default_Physical_Constant_Value;

public:
 ///\short Constructor: call the underlying constructors and 
 ///initialise the pointer to the Rayleigh number to address the default
 ///value of 0.0
 RefineableBuoyantQCrouzeixRaviartElement() : 
  RefineableQAdvectionDiffusionElement<DIM,3>(),
  RefineableQCrouzeixRaviartElement<DIM>()
  {
   Ra_pt = &Default_Physical_Constant_Value;
  }

 ///\short The required number of values stored at the nodes is 
 ///the sum of the required values of the two single-physics elements. This
 ///step is generic for any composed element of this type.
 inline unsigned required_nvalue(const unsigned &n) const
  {return (RefineableQAdvectionDiffusionElement<DIM,3>::required_nvalue(n) +
           RefineableQCrouzeixRaviartElement<DIM>::required_nvalue(n));}

 ///Access function for the Rayleigh number (const version)
 const double &ra() const {return *Ra_pt;}

 ///Access function for the pointer to the Rayleigh number
 double* &ra_pt() {return Ra_pt;}


 /// Final override for disable ALE
 void disable_ALE() 
  {
   //Disable ALE in both sets of equations
   RefineableNavierStokesEquations<DIM>::disable_ALE();
   RefineableAdvectionDiffusionEquations<DIM>::disable_ALE();
  }

 /// Final override for enable ALE
 void enable_ALE() 
  {
   //Enable ALE in both sets of equations
   RefineableNavierStokesEquations<DIM>::enable_ALE();
   RefineableAdvectionDiffusionEquations<DIM>::enable_ALE();
  }


 ///  Overload the standard output function with the broken default
 void output(ostream &outfile)
  {FiniteElement::output(outfile);}

 /// \short Output function:  
 ///  x,y,u   or    x,y,z,u at Nplot^DIM plot points
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
  
     // Output the temperature at the plot point
     outfile << this->interpolated_u_adv_diff(s) << " " << std::endl;
    }
   outfile << std::endl;
   
   // Write tecplot footer (e.g. FE connectivity lists)
   this->write_tecplot_zone_footer(outfile,nplot);
   
  }

 /// \short C-style output function:  Broken default
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


 /// \short Output function for a time-dependent exact solution.
 /// Broken default
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
 unsigned u_index_adv_diff() const 
  {return DIM;}
 
 /// \short Number of vertex nodes in the element is obtained from the
 /// geometric element.
 unsigned nvertex_node() const
  {return QElement<DIM,3>::nvertex_node();}

 /// \short Pointer to the j-th vertex node in the element,
 /// Call the geometric element's function.
 Node* vertex_node_pt(const unsigned& j) const
  {return QElement<DIM,3>::vertex_node_pt(j);}

 /// \short The total number of continously interpolated values is
 /// DIM+1 (DIM fluid velocities and one temperature).
 unsigned ncont_interpolated_values() const 
  {return DIM+1;}


 /// \short Get the continuously interpolated values at the local coordinate s.
 /// We choose to put the fluid velocities first, followed by the
 /// temperature.
 void get_interpolated_values(const Vector<double>&s,  Vector<double>& values)
  {
   //Storage for the fluid velocities
   Vector<double> nst_values;

   //Get the fluid velocities from the fluid element
   RefineableQCrouzeixRaviartElement<DIM>::
    get_interpolated_values(s,nst_values);

   //Storage for the temperature
   Vector<double> advection_values;

   //Get the temperature from the advection-diffusion element
   RefineableQAdvectionDiffusionElement<DIM,3>::
    get_interpolated_values(s,advection_values);
 
   //Add the fluid velocities to the values vector
   for(unsigned i=0;i<DIM;i++) {values.push_back(nst_values[i]);}  

   //Add the concentration to the end
   values.push_back(advection_values[0]);
  }


 
 /// \short Get all continuously interpolated values at the local 
 /// coordinate s at time level t (t=0: present; t>0: previous).
 /// We choose to put the fluid velocities first, followed by the
 /// temperature
 void get_interpolated_values(const unsigned& t, const Vector<double>&s,
                              Vector<double>& values)
  {
   //Storage for the fluid velocities
   Vector<double> nst_values;

   //Get the fluid velocities from the fluid element
   RefineableQCrouzeixRaviartElement<DIM>::
    get_interpolated_values(t,s,nst_values);

   //Storage for the temperature
   Vector<double> advection_values;

   //Get the temperature from the advection-diffusion element
   RefineableQAdvectionDiffusionElement<DIM,3>::
    get_interpolated_values(s,advection_values);

   //Add the fluid velocities to the values vector
   for(unsigned i=0;i<DIM;i++) {values.push_back(nst_values[i]);}   

   //Add the concentration to the end
   values.push_back(advection_values[0]);

  } // end of get_interpolated_values


 
 /// \short The additional hanging node information must be set up
 /// for both single-physics elements.
 void further_setup_hanging_nodes()
  {
   RefineableQCrouzeixRaviartElement<DIM>::further_setup_hanging_nodes();
   RefineableQAdvectionDiffusionElement<DIM,3>::further_setup_hanging_nodes();
  }


 
 /// \short Call the rebuild_from_sons functions for each of the
 /// constituent multi-physics elements.
 void rebuild_from_sons(Mesh* &mesh_pt) 
  {
   RefineableQAdvectionDiffusionElement<DIM,3>::rebuild_from_sons(mesh_pt);
   RefineableQCrouzeixRaviartElement<DIM>::rebuild_from_sons(mesh_pt);
  }
  


 /// \short Call the underlying single-physics element's further_build()
 /// functions and make sure that the pointer to the Rayleigh number
 /// is passed to the sons
 void further_build()
  {
   RefineableQCrouzeixRaviartElement<DIM>::further_build();
   RefineableQAdvectionDiffusionElement<DIM,3>::further_build();

   //Cast the pointer to the father element to the specific
   //element type
   RefineableBuoyantQCrouzeixRaviartElement<DIM>* cast_father_element_pt
    = dynamic_cast<RefineableBuoyantQCrouzeixRaviartElement<DIM>*>(
     this->father_element_pt());

   //Set the pointer to the Rayleigh number to be the same as that in
   //the father
   this->Ra_pt = cast_father_element_pt->ra_pt();
  } //end of further build


 
 /// The recovery order is that of the NavierStokes elements.
 unsigned nrecovery_order() 
  {return RefineableQCrouzeixRaviartElement<DIM>::nrecovery_order();}


 /// \short The number of Z2 flux terms is the same as that in 
 /// the fluid element.
 unsigned num_Z2_flux_terms()
  {
   return RefineableQCrouzeixRaviartElement<DIM>::num_Z2_flux_terms();
  }



 /// Get the Z2 flux from the fluid element
 void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
  {
   RefineableQCrouzeixRaviartElement<DIM>::get_Z2_flux(s,flux);
  } //end of get_Z2_flux



 /// \short Validate against exact solution at given time
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points and compute L2 error
 /// and L2 norm of velocity solution over element
 /// Overload to broken default
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
 /// Overload to broken default.
void compute_error(ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm)
  {FiniteElement::
   compute_error(outfile,exact_soln_pt,error,norm);}
 
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
 
 
 /// \short Overload the body force in the navier-stokes equations
 /// This provides the coupling from the advection-diffusion equations
 /// to the Navier--Stokes equations, the body force is the
 /// temperature multiplied by the Rayleigh number acting in the
 /// direction opposite to gravity.
 void get_body_force_nst(const double& time, const unsigned& ipt,
                         const Vector<double> &s, const Vector<double> &x,
                         Vector<double> &result)
  {
   // Get vector that indicates the direction of gravity from
   // the Navier-Stokes equations
   Vector<double> gravity(NavierStokesEquations<DIM>::g());
   
   // Temperature-dependent body force:
   for (unsigned i=0;i<DIM;i++)
    {
     result[i] = -gravity[i]*this->interpolated_u_adv_diff(s)*ra();
    }
  }

 /// Fill in the constituent elements' contribution to the residual vector.
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the residuals of the Navier-Stokes equations
   RefineableNavierStokesEquations<DIM>::fill_in_contribution_to_residuals(
    residuals);

   //Call the residuals of the advection-diffusion equations
   RefineableAdvectionDiffusionEquations<DIM>::
    fill_in_contribution_to_residuals(residuals);
  }


 ///\short Compute the element's residual Vector and the jacobian matrix
 ///using full finite differences, the default implementation
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
#ifdef USE_FD_JACOBIAN_FOR_REFINEABLE_BUOYANT_Q_ELEMENT
   FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
#else
   //Calculate the Navier-Stokes contributions (diagonal block and residuals)
   RefineableNavierStokesEquations<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //Calculate the advection-diffusion contributions 
   //(diagonal block and residuals)
   RefineableAdvectionDiffusionEquations<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //We now fill in the off-diagonal (interaction) blocks analytically
   this->fill_in_off_diagonal_jacobian_blocks_analytic(residuals,jacobian);
#endif
  } //End of jacobian calculation

 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   //Call the (broken) version in the base class
   FiniteElement::fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);
  }

 /// \short Compute the contribution of the off-diagonal blocks
 /// analytically.
 void fill_in_off_diagonal_jacobian_blocks_analytic(
  Vector<double> &residuals, DenseMatrix<double> &jacobian)
  {
   //Perform another loop over the integration loops using the information
   //from the original elements' residual assembly loops to determine
   //the conributions to the jacobian
   
   // Local storage for pointers to hang_info objects
   HangInfo *hang_info_pt=0, *hang_info2_pt=0;   
   
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
       double u_value = this->nodal_value(l,u_nodal_adv_diff);
       //Loop over the derivative directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_du_adv_diff_dx[j] += u_value*dpsifdx(l,j);
        }
      }
     
     //Assemble the Jacobian terms
     //---------------------------
     
     //Loop over the test functions/eqns
     for(unsigned l=0;l<n_node;l++)
      {
       //Local variables to store the number of master nodes and 
       //the weight associated with the shape function if the node is hanging
       unsigned n_master=1; 
       double hang_weight=1.0;
       
       //Local bool (is the node hanging)
       bool is_node_hanging = this->node_pt(l)->is_hanging();
       
       //If the node is hanging, get the number of master nodes
       if(is_node_hanging)
        {
         hang_info_pt = this->node_pt(l)->hanging_pt();
         n_master = hang_info_pt->nmaster();
        }
       //Otherwise there is just one master node, the node itself
       else 
        {
         n_master = 1;
        }
       
       //Loop over the master nodes
       for(unsigned m=0;m<n_master;m++)
        {
         //If the node is hanging get weight from master node
         if(is_node_hanging)
          {
           //Get the hang weight from the master node
           hang_weight = hang_info_pt->master_weight(m);
          }
         else
          {
           // Node contributes with full weight
           hang_weight = 1.0;
          }
         
         
         //Assemble derivatives of Navier Stokes momentum w.r.t. temperature
         //-----------------------------------------------------------------
         
         // Loop over velocity components for equations
         for(unsigned i=0;i<DIM;i++)
          {
           
           //Get the equation number
           if(is_node_hanging)
            {
             //Get the equation number from the master node
             local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                              u_nodal_nst[i]);
            }
           else
            {
             // Local equation number
             local_eqn = this->nodal_local_eqn(l,u_nodal_nst[i]);
            }
           
           if(local_eqn >= 0)
            {
             //Local variables to store the number of master nodes
             //and the weights associated with each hanging node
             unsigned n_master2=1; 
             double hang_weight2=1.0;
             
             //Loop over the nodes for the unknowns
             for(unsigned l2=0;l2<n_node;l2++)
              { 
               //Local bool (is the node hanging)
               bool is_node2_hanging = this->node_pt(l2)->is_hanging();
               
               //If the node is hanging, get the number of master nodes
               if(is_node2_hanging)
                {
                 hang_info2_pt = this->node_pt(l2)->hanging_pt();
                 n_master2 = hang_info2_pt->nmaster();
                }
               //Otherwise there is one master node, the node itself
               else
                {
                 n_master2 = 1;
                }
               
               //Loop over the master nodes
               for(unsigned m2=0;m2<n_master2;m2++)
                {
                 if(is_node2_hanging)
                  {
                   //Read out the local unknown from the master node
                   local_unknown = 
                    this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                         u_nodal_adv_diff);
                   //Read out the hanging weight from the master node
                   hang_weight2 = hang_info2_pt->master_weight(m2);
                  }
                 else
                  {
                   //The local unknown number comes from the node
                   local_unknown = this->nodal_local_eqn(l2,u_nodal_adv_diff);
                   //The hang weight is one
                   hang_weight2 = 1.0;
                  }
                 
                 if(local_unknown >= 0)
                  {
                   //Add contribution to jacobian matrix
                   jacobian(local_eqn,local_unknown) 
                    += -gravity[i]*psif(l2)*Ra*testf(l)* 
                    W*hang_weight*hang_weight2;
                  }
                }
              }
            }
          }
         
         
         //Assemble derivative of adv diff eqn w.r.t. fluid veloc
         //------------------------------------------------------
         {
          //Get the equation number
          if(is_node_hanging)
           {
            //Get the equation number from the master node
            local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                             u_nodal_adv_diff);
           }
          else
           {
            // Local equation number
            local_eqn = this->nodal_local_eqn(l,u_nodal_adv_diff);
           }
          
          //If it's not pinned
          if(local_eqn >= 0)
           {
            //Local variables to store the number of master nodes
            //and the weights associated with each hanging node
            unsigned n_master2=1; 
            double hang_weight2=1.0;
            
            //Loop over the nodes for the unknowns
            for(unsigned l2=0;l2<n_node;l2++)
             { 
              //Local bool (is the node hanging)
              bool is_node2_hanging = this->node_pt(l2)->is_hanging();
              
              //If the node is hanging, get the number of master nodes
              if(is_node2_hanging)
               {
                hang_info2_pt = this->node_pt(l2)->hanging_pt();
                n_master2 = hang_info2_pt->nmaster();
               }
              //Otherwise there is one master node, the node itself
              else
               {
                n_master2 = 1;
               }
              
              //Loop over the master nodes
              for(unsigned m2=0;m2<n_master2;m2++)
               {
                //If the node is hanging
                if(is_node2_hanging)
                 {
                  //Read out the hanging weight from the master node
                  hang_weight2 = hang_info2_pt->master_weight(m2);
                 }
                //If the node is not hanging
                else
                 {
                  //The hang weight is one
                  hang_weight2 = 1.0;
                 }
                
                //Loop over the velocity degrees of freedom
                for(unsigned i2=0;i2<DIM;i2++)
                 {
                  //If the node is hanging
                  if(is_node2_hanging)
                   {
                    //Read out the local unknown from the master node
                    local_unknown = 
                     this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                          u_nodal_nst[i2]);
                   }
                  else
                   {
                    //The local unknown number comes from the node
                    local_unknown=this->nodal_local_eqn(l2, u_nodal_nst[i2]);
                   }
                  
                  //If it's not pinned
                  if(local_unknown >= 0)
                   {
                    //Add the contribution to the jacobian matrix
                    jacobian(local_eqn,local_unknown)
                     -= Pe*psif(l2)*interpolated_du_adv_diff_dx[i2]*testf(l)
                     *W*hang_weight*hang_weight2;
                   }
                 }
               }
             }
           }
         }
         
        }
      }
    }
  } //End of function


};


//===================================================================
///Set the default value of the Rayleigh number to be zero
//=================================================================== 
template<>
double RefineableBuoyantQCrouzeixRaviartElement<2>::
Default_Physical_Constant_Value=0.0;



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



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



//=======start_of_problem_class=======================================
/// 2D Convection  problem on rectangular domain, discretised 
/// with refineable elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class RefineableConvectionProblem : public Problem
{

public:

 ///Constructor
 RefineableConvectionProblem();

 /// Destructor. Empty
 ~RefineableConvectionProblem() {}

 /// \short Update the problem specs before solve:
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// \short Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  } //end of access function to specic mesh

 /// Actions before adapt:(empty)
 void actions_before_adapt() {}

 /// \short Actions after adaptation,
 /// Re-pin a single pressure degree of freedom
 void actions_after_adapt()
  {
   //Unpin all the pressures to avoid pinning two pressures
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());

   //Pin the zero-th pressure dof in the zero-th element and set
   // its value to zero
   fix_pressure(0,0,0.0);
  }

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

 
 /// \short Access function to boolean flag that determines if the
 /// boundary condition on the upper wall is perturbed slightly
 /// to force the solution into the symmetry broken state.
 bool& use_imperfection() {return Imperfect;}
 
 /// \short Doc the solution.
 void doc_solution();
 
private:
 
 /// DocInfo object
 DocInfo Doc_info;
 
 /// Is the boundary condition imperfect or not
 bool Imperfect;

}; // end of problem class


//=======start_of_constructor=============================================
/// Constructor for adaptive thermal convection problem
//========================================================================
template<class ELEMENT>
RefineableConvectionProblem<ELEMENT>::
RefineableConvectionProblem() : Imperfect(false)
{ 
 // Set output directory
 Doc_info.set_directory("RESLT");
 
 // # of elements in x-direction
 unsigned n_x=9;

 // # of elements in y-direction
 unsigned n_y=8;

 // Domain length in x-direction
 double l_x=3.0;

 // Domain length in y-direction
 double l_y=1.0;
 
 // Build the mesh
 RefineableRectangularQuadMesh<ELEMENT>* cast_mesh_pt =
  new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 //Set the problem's mesh pointer
 Problem::mesh_pt() = cast_mesh_pt;


 // Create/set error estimator
 cast_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Set error targets for adaptive refinement
 cast_mesh_pt->max_permitted_error()=1.0e-4; 
 cast_mesh_pt->min_permitted_error()=0.5e-4; 


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
 
 // Pin the zero-th pressure value in the zero-th element and
 // set its value to zero.
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

   // Set the Reynolds number (1/Pr in our non-dimensionalisation)
   el_pt->re_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set the Rayleigh number
   el_pt->ra_pt() = &Global_Physical_Variables::Rayleigh;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor


//===================start_actions_before_newton_solve===========================
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to include an imperfection (or not) depending on the control flag.
//========================================================================
template<class ELEMENT>
void RefineableConvectionProblem<ELEMENT>::actions_before_newton_solve()
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
     //If we are on the side walls we only pin the x-velocity.
     if((ibound==1) || (ibound==3)) {vel_max = 1;}
     //Set the pinned velocities to zero
     for(unsigned j=0;j<vel_max;j++) {nod_pt->set_value(j,0.0);}

     //If we are on the top boundary
     if(ibound==2) 
      {
       //Set the temperature to -0.5 (cooled)
       nod_pt->set_value(2,-0.5);
       //Add small velocity imperfection if desired
       if(Imperfect)
        {
         //Read out the x position
         double x = nod_pt->x(0);
         //Set a sinusoidal perturbation in the vertical velocity
         //This perturbation is mass conserving
         double value = sin(2.0*3.141592654*x/3.0);
         nod_pt->set_value(1,value);
        }
      }

     //If we are on the bottom boundary, set the temperature
     //to 0.5 (heated)
     if(ibound==0) {nod_pt->set_value(2,0.5);}
    }
  }

}  // end of actions before solve



//====================start_of_doc_solution===============================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableConvectionProblem<ELEMENT>::doc_solution()
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


//===============start_of_main========================================
/// Driver code for 2D Boussinesq convection problem with 
/// adaptivity.
//====================================================================
int main()
{

 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;

 // Create the problem with 2D nine-node refineable elements.
 RefineableConvectionProblem<
  RefineableBuoyantQCrouzeixRaviartElement<2> > problem;
 
 // Apply a perturbation to the upper boundary condition to
 // force the solution into the symmetry-broken state.
 problem.use_imperfection() = true;
 
 //Solve the problem with (up to) two levels of adaptation
 problem.newton_solve(2);
 
 //Document the solution
 problem.doc_solution();
 
 // Make the boundary conditions perfect and solve again. 
 // Now the slightly perturbed symmetry broken state computed
 // above is used as the initial condition and the Newton solver
 // converges to the symmetry broken solution, even without
 // the perturbation
 problem.use_imperfection() = false;
 problem.newton_solve(2);
 problem.doc_solution();

} // end of main









