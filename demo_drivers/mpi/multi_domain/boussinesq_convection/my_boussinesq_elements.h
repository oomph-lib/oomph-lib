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
//Header for a multi-physics problem that couples a Navier--Stokes
//mesh to an advection diffusion mesh, giving Boussinesq convection

//Oomph-lib headers, we require the generic, advection-diffusion
//and navier-stokes elements.
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;

//======================class definitions==============================
/// Build MyRefineableNavierStokesElement that inherits from 
/// ElementWithExternalElement so that it can "communicate" with 
/// MyRefineableAdvectionDiffusionElement
//=====================================================================
template<unsigned DIM>
class MyRefineableNavierStokesElement : 
 public virtual RefineableQCrouzeixRaviartElement<DIM>,
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
 MyRefineableNavierStokesElement() : RefineableQCrouzeixRaviartElement<DIM>(),
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
  {return RefineableQCrouzeixRaviartElement<DIM>::required_nvalue(n);}

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

 /// \short Call the underlying single-physics element's further_build()
 /// functions and make sure that the pointer to the Rayleigh number
 /// is passed to the sons
 void further_build()
  {
   RefineableQCrouzeixRaviartElement<DIM>::further_build();

   //Cast the pointer to the father element to the specific
   //element type
   MyRefineableNavierStokesElement<DIM>* cast_father_element_pt
    = dynamic_cast<MyRefineableNavierStokesElement<DIM>*>(
     this->father_element_pt());

   //Set the pointer to the Rayleigh number to be the same as that in
   //the father
   this->Ra_pt = cast_father_element_pt->ra_pt();
  }

 /// Does something need to be done about the Z2 flux? In order to test
 /// against the single-mesh problem this just needs to come from this mesh,
 /// so don't overload this one - overload the other one to have the same
 /// values as this one, somehow...

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

 /// Fill in the constituent elements' contribution to the residual vector.
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the residuals of the Navier-Stokes equations
   RefineableNavierStokesEquations<DIM>::fill_in_contribution_to_residuals(
    residuals);
  }

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

 /// \short Overload assign_all_generic_local_equation_numbers to
 ///        strip out external data and add back unique data
 void assign_all_generic_local_eqn_numbers()
  {
   //External data may not be distinct from nodal data depending upon
   //the source element, so call helper to remove non-unique external data
   assign_unique_external_data_helper();
   //Now call the refineable element equation num. (int, ext, nodal)
   RefineableElement::assign_all_generic_local_eqn_numbers();
  }

};

//======================class definitions==============================
/// Build MyAdvectionDiffusionElement that inherits from 
/// ElementWithExternalElement so that it can "communicate" with the 
/// MyNavierStokesElement
//=====================================================================
template<unsigned DIM> // second templated argument?
class MyRefineableAdvectionDiffusionElement : 
 public virtual RefineableQAdvectionDiffusionElement<DIM,3>,
 public virtual ElementWithExternalElement
{

public:

 /// \short Constructor: call the underlying constructors
 MyRefineableAdvectionDiffusionElement() : 
  RefineableQAdvectionDiffusionElement<DIM,3>(), ElementWithExternalElement()
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
  {return RefineableQAdvectionDiffusionElement<DIM,3>::required_nvalue(n);}

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

 /// Testing firstly what happens when each mesh has the same refinement
 /// pattern as the single mesh with the combined element.  In that case
 /// the recovery order and flux terms were solely based on the fluid.

//  /// The recovery order is that of the NavierStokes elements.
//  unsigned nrecovery_order() 
//   {return RefineableQCrouzeixRaviartElement<DIM>::nrecovery_order();}

//  /// \short The number of Z2 flux terms is the same as that in 
//  /// the NavierStokes (fluid) element.
//  unsigned num_Z2_flux_terms()
//   {
//    return RefineableQCrouzeixRaviartElement<DIM>::num_Z2_flux_terms();
//   }

//  /// Get the Z2 flux from the fluid element
//  void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
//   {
//    RefineableQCrouzeixRaviartElement<DIM>::get_Z2_flux(s,flux);
//   } //end of get_Z2_flux

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

 /// Just call the fill_in_residuals for AdvDiff
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   RefineableAdvectionDiffusionEquations<DIM>::
    fill_in_contribution_to_residuals(residuals);
  }

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

 /// \short Overload assign_all_generic_local_equation_numbers to
 ///        strip out external data and add back unique data
 void assign_all_generic_local_eqn_numbers()
  {
   //External data may not be distinct from nodal data depending upon
   //the source element, so call helper to remove non-unique external data
   assign_unique_external_data_helper();
   //Now call the refineable element equation num. (int, ext, nodal)
   RefineableElement::assign_all_generic_local_eqn_numbers();
  }

};

//============================================================
// Overload get_body_force_nst to get the temperature "body force"
// from the "source" AdvectionDiffusion element via current integration point
//========================================================
template<unsigned DIM>
void MyRefineableNavierStokesElement<DIM>::get_body_force_nst
(const double& time,const unsigned& ipt,const Vector<double> &s,
 const Vector<double> &x,Vector<double> &result)
{
 // The interaction is stored at index 0 of the AD element
 unsigned interaction=0;

 // Dynamic cast the source element at this integration point to correct type
 MyRefineableAdvectionDiffusionElement<DIM>* source_el_pt=
  dynamic_cast<MyRefineableAdvectionDiffusionElement<DIM>*>
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
void MyRefineableAdvectionDiffusionElement<DIM>::get_wind_adv_diff
(const unsigned& ipt,const Vector<double> &s,const Vector<double>& x, 
 Vector<double>& wind) const
{
 // The interatction is stored at index 0 of the NST element
 unsigned interaction=0;

 // Dynamic cast "other" element to correct type
 MyRefineableNavierStokesElement<DIM>* source_el_pt=
  dynamic_cast<MyRefineableNavierStokesElement<DIM>*>
  (external_element_pt(interaction,ipt));

 //The wind function is simply the velocity at the points of the source element
 source_el_pt->interpolated_u_nst
  (external_element_local_coord(interaction,ipt),wind);
}  


//=========================================================================
/// Set the default physical value to be zero in 2D and 3D
//=========================================================================
template<>
double MyRefineableNavierStokesElement<2>::Default_Physical_Constant_Value=0.0;

template<>
double MyRefineableNavierStokesElement<3>::Default_Physical_Constant_Value=0.0;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


