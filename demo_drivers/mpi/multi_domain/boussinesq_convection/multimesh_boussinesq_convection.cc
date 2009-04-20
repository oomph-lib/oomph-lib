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
//Driver for a multi-physics problem that couples a Navier--Stokes
//mesh to an advection diffusion mesh, giving Boussinesq convection

//Oomph-lib headers, we require the generic, advection-diffusion,
//and navier-stokes elements.
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"

// Both meshes are the standard rectangular quadmesh
#include "meshes/rectangular_quadmesh.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;


//======================class definitions==============================
/// Build MyNavierStokesElement that inherits from ElementWithExternalElement
/// so that it can "communicate" with MyAdvectionDiffusionElement
//=====================================================================
template<unsigned DIM>
class MyNavierStokesElement : 
 public virtual QCrouzeixRaviartElement<DIM>,
 public virtual ElementWithExternalElement
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
 MyNavierStokesElement() : QCrouzeixRaviartElement<DIM>(),
                           ElementWithExternalElement()
  {
   Ra_pt = &Default_Physical_Constant_Value;

   // Setup the storage for the interaction between elements
   unsigned n_interaction=1;
   unsigned nint_pt=integral_pt()->nweight();
   // The dimension of the source element is the same as this element
   unsigned n_dim_source=ndim();

   initialise_external_element_storage(n_interaction,nint_pt,n_dim_source);
  } 

 ///\short The required number of values stored at the nodes is the number of
 ///required values of the CrouzeixRaviartElement.
 unsigned required_nvalue(const unsigned &n) const
  {return QCrouzeixRaviartElement<DIM>::required_nvalue(n);}

 ///Access function for the Rayleigh number (const version)
 const double &ra() const {return *Ra_pt;}

 ///Access function for the pointer to the Rayleigh number
 double* &ra_pt() {return Ra_pt;}

 ///  Overload the standard output function with the broken default
 void output(ostream &outfile) {FiniteElement::output(outfile);}

 /// \short Output function:  
 ///  Output x, y, u, v, p at Nplot^DIM plot points
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
     outfile << this->interpolated_p_nst(s)  << " " << std::endl;;
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

 /// \short global position vector at local s (return interpolated_x)
 void position(const Vector<double>& s, Vector<double>& r) const
  {
   // Get the position vector using interpolated_x
   interpolated_x(s,r);
  }

 // Overload get_body_force_nst to get the temperature "body force"
 // from the "source" AdvectionDiffusion element via current integration point
 void get_body_force_nst(const double& time, const unsigned& ipt, 
                         const Vector<double> &s, const Vector<double> &x, 
                         Vector<double> &result);

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

};

//======================class definitions==============================
/// Build MyAdvectionDiffusionElement that inherits from 
/// ElementWithExternalElement
/// so that it can "communicate" with the MyNavierStokesElement
//=====================================================================
template<unsigned DIM>
class MyAdvectionDiffusionElement : 
 public virtual QAdvectionDiffusionElement<DIM,3>,
 public virtual ElementWithExternalElement
{

public:

 /// \short Constructor: call the underlying constructors
 MyAdvectionDiffusionElement() : QAdvectionDiffusionElement<DIM,3>(),
                                 ElementWithExternalElement()
  { 
   // Setup the storage for the interaction between elements
   unsigned n_interaction=1;
   unsigned nint_pt=integral_pt()->nweight();
   // The dimension of the source element is the same as this element
   unsigned n_dim_source=ndim();

   initialise_external_element_storage(n_interaction,nint_pt,n_dim_source);
  } 

 ///\short The required number of values stored at the nodes is the number of
 ///required values of the AdvectionDiffusionElement.
 unsigned required_nvalue(const unsigned &n) const
  {return QAdvectionDiffusionElement<DIM,3>::required_nvalue(n);}

 ///  Overload the standard output function with the broken default
 void output(ostream &outfile) {FiniteElement::output(outfile);}

 /// \short Output function:  
 ///  Output x, y, theta at Nplot^DIM plot points
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

  /// \short global position vector at local s (return interpolated_x)
 void position(const Vector<double>& s, Vector<double>& r) const
  {
   // Get the position vector using interpolated_x
   interpolated_x(s,r);
  }

 /// \short Overload the wind function in the advection-diffusion equations.
 /// This provides the coupling from the Navier--Stokes equations to the
 /// advection-diffusion equations because the wind is the fluid velocity,
 /// obtained from the source element in the other mesh
 void get_wind_adv_diff(const unsigned& ipt, const Vector<double> &s, 
                        const Vector<double>& x, Vector<double>& wind) const;

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

};

//============================================================
// Overload get_body_force_nst to get the temperature "body force"
// from the "source" AdvectionDiffusion element via current integration point
//========================================================
template<unsigned DIM>
void MyNavierStokesElement<DIM>::get_body_force_nst(const double& time, 
                                                    const unsigned& ipt,
                                                    const Vector<double> &s,
                                                    const Vector<double> &x,
                                                    Vector<double> &result)
{
 // The interaction index is 0 in this case
 unsigned interaction=0;

 // Dynamic cast "other" element to correct type
 MyAdvectionDiffusionElement<DIM>* source_el_pt=
  dynamic_cast<MyAdvectionDiffusionElement<DIM>*>
  (external_element_pt(interaction,ipt));

 // Get vector that indicates the direction of gravity from
 // the Navier-Stokes equations
 Vector<double> gravity(NavierStokesEquations<DIM>::g());
   
 // Temperature-dependent body force:
 for (unsigned i=0;i<DIM;i++)
  {
   result[i] = -gravity[i]*source_el_pt->interpolated_u_adv_diff
    (external_element_local_coord(interaction,ipt))*ra();
  }
}

//==========================================================================
/// \short Overload the wind function in the advection-diffusion equations.
/// This provides the coupling from the Navier--Stokes equations to the
/// advection-diffusion equations because the wind is the fluid velocity,
/// obtained from the source elements in the other mesh
//==========================================================================
template<unsigned DIM>
void MyAdvectionDiffusionElement<DIM>::get_wind_adv_diff
(const unsigned& ipt,const Vector<double> &s,const Vector<double>& x, 
 Vector<double>& wind) const
{
 // The interaction index is 0 in this case
 unsigned interaction=0;

 // Dynamic cast "other" element to correct type
 MyNavierStokesElement<DIM>* source_el_pt=
  dynamic_cast<MyNavierStokesElement<DIM>*>
  (external_element_pt(interaction,ipt));

 //The wind function is simply the velocity at the points of the "other" el
 source_el_pt->interpolated_u_nst
  (external_element_local_coord(interaction,ipt),wind);
}  

//=========================================================================
/// Set the default physical value to be zero
//=========================================================================
template<>
double MyNavierStokesElement<2>::Default_Physical_Constant_Value=0.0;


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
/// 2D Convection  problem on two rectangular domains, discretised 
/// with Navier-Stokes and Advection-Diffusion elements. The specific type
/// of elements is specified via the template parameters.
//====================================================================
template<class NST_ELEMENT,class AD_ELEMENT> 
class ConvectionProblem : public Problem //public virtual HelmholtzProblem
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

 /// Actions after distribute: set sources
 void actions_after_distribute()
  {
   // Set binning parameters
//   Multi_domain_functions::Setup_bins=true;
//   Multi_domain_functions::N_bin_dim=5;

   // Set interaction indices (default now)
//   unsigned interaction_nst=0;
//   unsigned interaction_ad=0;

   // Set sources
   Multi_domain_functions::set_sources<NST_ELEMENT,AD_ELEMENT,2,2>
    (this,nst_mesh_pt(),adv_diff_mesh_pt());
  }

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
   //Cast to specific element and fix pressure (only for NST element...?)
   dynamic_cast<NST_ELEMENT*>(nst_mesh_pt()->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

 /// \short Doc the solution.
 void doc_solution();

 /// \short Set the boundary conditions
 void set_boundary_conditions(const double &time);

 /// \short Access function to the Navier-Stokes mesh
 RectangularQuadMesh<NST_ELEMENT>* nst_mesh_pt() 
  {
   return dynamic_cast<RectangularQuadMesh<NST_ELEMENT>*>(Nst_mesh_pt);
  }

 /// \short Access function to the Advection-Diffusion mesh
 RectangularQuadMesh<AD_ELEMENT>* adv_diff_mesh_pt() 
  {
   return dynamic_cast<RectangularQuadMesh<AD_ELEMENT>*>(Adv_diff_mesh_pt);
  }
 
private:
 
 /// DocInfo object
 DocInfo Doc_info;

protected:

 RectangularQuadMesh<NST_ELEMENT>* Nst_mesh_pt;
 RectangularQuadMesh<AD_ELEMENT>* Adv_diff_mesh_pt;

}; // end of problem class

//===========start_of_constructor=========================================
/// \short Constructor for convection problem
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
ConvectionProblem<NST_ELEMENT,AD_ELEMENT>::ConvectionProblem()
{
 // Suppress warnings about repeated external data
 GeneralisedElement::Suppress_warning_about_repeated_external_data=true;

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

 // Build two standard rectangular quadmesh
 Nst_mesh_pt = 
  new RectangularQuadMesh<NST_ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());
 Adv_diff_mesh_pt = 
  new RectangularQuadMesh<AD_ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here

 //Loop over the boundaries
 unsigned num_bound = nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the maximum index to be pinned (all values by default)
   unsigned val_max;//=3; (fine for combined element... !)

   //Loop over the number of nodes on the boundry
   unsigned num_nod= nst_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //If we are on the side-walls, the v-velocity and temperature
     //satisfy natural boundary conditions, so we only pin the
     //first value
     if ((ibound==1) || (ibound==3)) 
      {
       val_max=1;
      }
     else // pin all values
      {
       val_max=nst_mesh_pt()->boundary_node_pt(ibound,inod)->nvalue();
      }

     //Loop over the desired values stored at the nodes and pin
     for(unsigned j=0;j<val_max;j++)
      {
       nst_mesh_pt()->boundary_node_pt(ibound,inod)->pin(j);
      }
    }
  }

 //Pin the zero-th pressure dof in element 0 and set its value to
 //zero:
 fix_pressure(0,0,0.0);

 //Loop over the boundaries of the AD mesh
 num_bound = adv_diff_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the maximum index to be pinned (all values by default)
   unsigned val_max;//=3;

   //Loop over the number of nodes on the boundry
   unsigned num_nod= adv_diff_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //If we are on the side-walls, the v-velocity and temperature
     //satisfy natural boundary conditions, so we don't pin anything
     // in this mesh
     if ((ibound==1) || (ibound==3)) 
      {
       val_max=0;
      }
     else // pin all values
      {
       val_max=adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->nvalue();
       //Loop over the desired values stored at the nodes and pin
       for(unsigned j=0;j<val_max;j++)
        {
         adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->pin(j);
        }
      }
    }
  }


 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructors. 
 unsigned n_nst_element = nst_mesh_pt()->nelement();
 for(unsigned i=0;i<n_nst_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   NST_ELEMENT *el_pt = dynamic_cast<NST_ELEMENT*>
    (nst_mesh_pt()->element_pt(i));

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

 unsigned n_ad_element = adv_diff_mesh_pt()->nelement();
 for(unsigned i=0;i<n_ad_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   AD_ELEMENT *el_pt = dynamic_cast<AD_ELEMENT*>
    (adv_diff_mesh_pt()->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;

   // Set the Peclet number multiplied by the Strouhal number
   el_pt->pe_st_pt() =&Global_Physical_Variables::Peclet;

   //The mesh is fixed, so we can disable ALE
   el_pt->disable_ALE();

   // Set pointer to the continuous time
   el_pt->time_pt() = time_pt();
  }

 // combine the submeshes
 add_sub_mesh(Nst_mesh_pt);
 add_sub_mesh(Adv_diff_mesh_pt);
 build_global_mesh();

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor



//===========start_of_set_boundary_conditions================
/// Set the boundary conditions as a function of continuous 
/// time
//===========================================================
template<class NST_ELEMENT,class AD_ELEMENT>
void ConvectionProblem<NST_ELEMENT,AD_ELEMENT>::set_boundary_conditions(
 const double &time)
{
 // Loop over all the boundaries on the NST mesh
 unsigned num_bound=nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=nst_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=nst_mesh_pt()->boundary_node_pt(ibound,inod);

     //Set the number of velocity components
     unsigned vel_max=2;

     //If we are on the side walls we only set the x-velocity.
     if((ibound==1) || (ibound==3)) {vel_max = 1;}

     //Set the pinned velocities to zero on NST mesh
     for(unsigned j=0;j<vel_max;j++) {nod_pt->set_value(j,0.0);}

     //If we are on the top boundary
     if(ibound==2) 
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

    }
  }

 // Loop over all the boundaries on the AD mesh
 num_bound=adv_diff_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=adv_diff_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=adv_diff_mesh_pt()->boundary_node_pt(ibound,inod);

     //If we are on the top boundary, set the temperature 
     //to -0.5 (cooled)
     if(ibound==2) {nod_pt->set_value(0,-0.5);}

     //If we are on the bottom boundary, set the temperature
     //to 0.5 (heated)
     if(ibound==0) {nod_pt->set_value(0,0.5);}
    }
  }


} // end_of_set_boundary_conditions

//===============start_doc_solution=======================================
/// Doc the solution
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
void ConvectionProblem<NST_ELEMENT,AD_ELEMENT>::doc_solution()
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output whole solution (this will output elements from one mesh
 //----------------------  followed by the other mesh at the moment...?)
 sprintf(filename,"%s/soln%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),MPI_Helpers::My_rank);
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

 // Output solution for each mesh
 //------------------------------
 sprintf(filename,"%s/vel_soln%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),MPI_Helpers::My_rank);
 some_file.open(filename);
 nst_mesh_pt()->output(some_file,npts);
 some_file.close();

 sprintf(filename,"%s/temp_soln%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),MPI_Helpers::My_rank);
 some_file.open(filename);
 adv_diff_mesh_pt()->output(some_file,npts);
 some_file.close();

 Doc_info.number()++;
} // end of doc


//=======start_of_main================================================
/// Driver code for 2D Boussinesq convection problem
//====================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;

 //Construct our problem
 ConvectionProblem<MyNavierStokesElement<2>,MyAdvectionDiffusionElement<2> > 
  problem;

 // Apply the boundary condition at time zero
 problem.set_boundary_conditions(0.0);

 // Distribute the problem (and set the sources)
#ifdef OOMPH_HAS_MPI
 DocInfo mesh_doc_info;
 bool report_stats=true;
 problem.distribute(mesh_doc_info,report_stats);
#endif

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
// unsigned n_steps = 20;

 //If we have a command line argument, perform fewer steps 
 //(used for self-test runs)
 if(argc > 1) {n_steps = 5;}

 //Perform n_steps timesteps
 for(unsigned i=0;i<n_steps;++i)
  {
   problem.unsteady_newton_solve(dt);
   problem.doc_solution();
  }

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

} // end of main









