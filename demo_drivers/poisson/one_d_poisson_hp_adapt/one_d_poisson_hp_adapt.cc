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
// Driver for a simple 1D Poisson problem with hp-adaptive mesh refinement

// Generic oomph-lib routines
#include "generic.h"

// Poisson elements/equations
#include "poisson.h"

// The mesh
#include "meshes/one_d_mesh.h"

using namespace std;

using namespace oomph;



//======start_of_namespace================================================
/// Namespace for exact solution for Poisson equation with "sharp step" 
//========================================================================
namespace ArcTanSolnForPoisson
{
 
 /// Parameter for steepness of "step"
 double Alpha=100.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = atan(Alpha*(x[0]-0.5));
 }

 /// Exact gradient as a Vector
 void get_exact_gradient(const Vector<double>& x, Vector<double>& dudx)
 {
  dudx[0] = Alpha/(1.0+(Alpha*(x[0]-0.5))*(Alpha*(x[0]-0.5)));
 }

 /// Source function required to make the solution above an exact solution 
 void get_source(const Vector<double>& x, double& source)
 {
  // Cache tan( atan(Alpha*(x[0]-0.5)) ) term
  double tan_term = tan( atan(Alpha*(x[0]-0.5)) );
  
  // Compute source function
  source = -(2.0 * Alpha*Alpha * tan_term) /
   ( (1.0 + tan_term*tan_term)*(1.0 + tan_term*tan_term) );
 }
 
} // End of namespace



//=============================================================
/// A class for all elements that solve the Poisson equations.
/// \f[ 
/// \frac{\partial^2 u}{\partial x_i^2} = f(x_j)
/// \f] 
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
template <unsigned DIM>
class ModalPoissonEquations : public virtual FiniteElement
{
 
public:
 
 /// Function pointer to source function fct(x,f(x)) -- 
 /// x is a Vector! 
 typedef void (*PoissonSourceFctPt)(const Vector<double>& x, double& f);
 
 /// Function pointer to a diffusion function fxt(x,f(x))
 typedef void (*PoissonDiffFctPt)(const Vector<double> &x, double &f);
 
 /// Constructor (must initialise the Source_fct_pt to null)
 ModalPoissonEquations() : Source_fct_pt(0), Diff_fct_pt(0)  {}
 
 /// Access function: Nodal function value at local node n
 /// Uses suitably interpolated value for hanging nodes.
 virtual double u(const unsigned& n) const = 0;
 
 /// Number of basis functions
 virtual unsigned nbasis() const=0;
 
 virtual void basis(const Vector<double> &s,
                        Shape &basis) const=0;
 
 /// Output with default number of plot points
 void output(ostream &outfile) 
  {
   unsigned nplot=5;
   output(outfile,nplot);
  }
 
 /// Output FE representation of soln: x,y,u or x,y,z,u at 
 /// Nplot^DIM plot points
 void output(ostream &outfile, const unsigned &nplot)
  {
   //Vector of local coordinates
   Vector<double> s(DIM);
   
   //Find out how many nodes there are
   unsigned n_node = nnode();
   //Find out the total number of basis and test functions there are
   unsigned n_basis = nbasis();
   
   //Set up memory for the shape and test functions
   Shape psi(n_node), basis(n_basis), test(n_basis);
   DShape dpsidx(n_node,DIM), dbasisdx(n_basis,DIM), dtestdx(n_basis,DIM);
   
   // Tecplot header info
   outfile << tecplot_zone_string(nplot);
   
   // Loop over plot points
   unsigned num_plot_points=nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     
     // Get local coordinates of plot point
     get_s_plot(iplot,nplot,s);

     //get_basis(s,basis);
     dshape_dbasis_and_dtest_eulerian(s,psi,dpsidx,basis,
                                      dbasisdx,test,dtestdx);

     for(unsigned i=0;i<DIM;i++) 
      {outfile << interpolated_x(s,i) << " ";}
     
     outfile << interpolated_u(s) << std::endl;   
    }

  }



 /// Output exact soln: x,y,u_exact or x,y,z,u_exact at nplot^DIM plot points
 void output_fct(ostream &outfile, const unsigned &nplot, 
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
  }


 /// Output exact soln: x,y,u_exact or x,y,z,u_exact at 
 /// nplot^DIM plot points (dummy time-dependent version to 
 /// keep intel compiler happy)
 virtual void output_fct(ostream &outfile, const unsigned &nplot,
                         const double& time, 
                         FiniteElement::UnsteadyExactSolutionFctPt 
                         exact_soln_pt)
  {

   oomph_info << "No time-dep. output_fct() for Poisson elements " << std::endl;
   exit(1);
  }


 /// Get error against and norm of exact solution
 void compute_error(ostream &outfile, 
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
   
   //Set the value of n_intpt
   unsigned n_intpt = integral_pt()->nweight();
   
   // Setup output structure: Conversion is fishy but it's only output...
   unsigned nplot;
   if (DIM==1)
    {
     nplot=n_intpt;
    }
   else 
    {
     nplot=unsigned(pow(n_intpt,1.0/double(DIM)));
    }
   
   // Tecplot header info
   outfile << tecplot_zone_string(nplot);
   
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
     double u_fe=interpolated_u(s);
     
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


 /*/// Get error against and norm of exact solution
 virtual void compute_energy_error(ostream &outfile, 
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
   
   //Set the value of n_intpt
   unsigned n_intpt = integral_pt()->nweight();
   
   // Setup output structure: Conversion is fishy but it's only output...
   unsigned nplot;
   if (DIM==1)
    {
     nplot=n_intpt;
    }
   else 
    {
     nplot=unsigned(pow(n_intpt,1.0/double(DIM)));
    }
   
   // Tecplot header info
   outfile << tecplot_zone_string(nplot);
   
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
     double u_fe=interpolated_u(s);
     
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
  }*/
 
 
 /// Dummy, time dependent error checker
 void compute_error(ostream &outfile, 
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time, double& error, double& norm)
  {
   oomph_info << "No time-dep. compute_error() for Poisson elements " << std::endl;
   exit(1);
  }
 
 
 
 /// Access function: Pointer to source function
 PoissonSourceFctPt& source_fct_pt() 
  {
   return Source_fct_pt;
  }
 
 
 /// Access function: Pointer to source function. Const version
 PoissonSourceFctPt source_fct_pt() const
  {
   return Source_fct_pt;
  }
 
 
 /// Get source term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the source function might be determined by
 /// another system of equations 
 inline virtual void get_source_poisson(const Vector<double>& x, double& source) const
  {
   //If no source function has been set, return zero
   if(Source_fct_pt==0) {source = 0.0;}
   else
    {
     // Get source strength
     (*Source_fct_pt)(x,source);
    }
  }
 
 /// Access function: Pointer to diffusivity function
 PoissonDiffFctPt& diff_fct_pt() 
  {
   return Diff_fct_pt;
  }
 
 
 /// Access function: Pointer to source function. Const version
 PoissonDiffFctPt diff_fct_pt() const
  {
   return Diff_fct_pt;
  }
 
 
 /// Get diffusivity term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the diffusivity function might be determined by
 /// another system of equations 
 inline virtual void get_diff(const Vector<double>& x, double& diff) const
  {
   //If no source function has been set, return one
   if(Diff_fct_pt==0) {diff = 1.0;}
   else
    {
     // Get source strength
     (*Diff_fct_pt)(x,diff);
    }
  }


 /// Get flux: flux[i] = du/dx_i
 // Needs to know the p_order of the element to use the basis functions
 void get_flux(const Vector<double>& s, Vector<double>& flux, unsigned p_order) const
  {
   //Find out how many nodes there are
   unsigned n_node = nnode();
   //Find out the total number of basis and test functions there are
   unsigned n_basis = nbasis();
   
   //Set up memory for the shape and test functions
   Shape psi(n_node), basis(n_basis), test(n_basis);
   DShape dpsidx(n_node,DIM), dbasisdx(n_basis,DIM), dtestdx(n_basis,DIM);
   
   //Call the derivatives of the shape and test functions
   //dshape_eulerian(s,psi,dpsidx);
   dshape_dbasis_and_dtest_eulerian(s, psi, dpsidx, basis, dbasisdx, test, dtestdx);
     
   //Initialise to zero
   for(unsigned j=0;j<DIM;j++)
    {
     flux[j] = 0.0;
    }
   
   // Loop over modes
   for(unsigned l=0;l<p_order;l++) 
    {
     //Loop over derivative directions
     for(unsigned j=0;j<DIM;j++)
      {
       //flux[j] += u(l)*dpsidx(l,j);
       flux[j] += u(l)*dbasisdx(l,j);
      }
    }
  }
 
 
 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Create a dummy matrix
   DenseMatrix<double> dummy(1);
 
   //Call the generic residuals function with flag set to 0
   add_generic_residual_contribution(residuals,dummy,0);
  }
 
 
 /// Add the element's contribution to its residual vector and 
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   add_generic_residual_contribution(residuals,jacobian,1);
  }
 

 /// Return FE representation of function value u(s) at local coordinate s
 inline double interpolated_u(const Vector<double> &s) const
  {
   //Find number of basis functions
   unsigned n_basis = nbasis();

   //Local basis functions
   Shape psi(n_basis);

   //Find values of basis function
   basis(s,psi);

   //Initialise value of u
   double interpolated_u = 0.0;

   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_basis;l++) 
    {
     interpolated_u+=u(l)*psi(l);
    }

   return(interpolated_u);
  }


 /// Self-test: Return 0 for OK
 unsigned self_test() {return 0;}

protected:

 /// Shape/test functions and derivs w.r.t. to global coords at 
 /// local coord. s; return  Jacobian of mapping
 virtual double dshape_dbasis_and_dtest_eulerian(const Vector<double> &s, 
                                                 Shape &psi, 
                                                 DShape &dpsidx, 
                                                 Shape &basis,
                                                 DShape &dbasisdx,
                                                 Shape &test, 
                                                 DShape &dtestdx) const=0;

 /// Shape/test functions and derivs w.r.t. to global coords at 
 /// integration point ipt; return  Jacobian of mapping
 virtual double dshape_dbasis_and_dtest_eulerian_at_knot(const unsigned &ipt, 
                                                 Shape &psi, 
                                                 DShape &dpsidx,
                                                 Shape &basis,
                                                 DShape &dbasisdx,
                                                 Shape &test, 
                                                 DShape &dtestdx) const=0;

 /// Compute element residual Vector only (if flag=and/or element 
 /// Jacobian matrix 
 virtual void add_generic_residual_contribution(Vector<double> &residuals, 
                                                 DenseMatrix<double> &jacobian, 
                                                 unsigned flag) 
  {
   //Find out how many nodes there are
   unsigned n_node = nnode();
   //Find out the total number of basis and test functions there are
   unsigned n_basis = nbasis();
   
   //Set up memory for the shape and test functions
   Shape psi(n_node), basis(n_basis), test(n_basis);
   DShape dpsidx(n_node,DIM), dbasisdx(n_basis,DIM), dtestdx(n_basis,DIM);

   //Set the value of n_intpt
   unsigned n_intpt = integral_pt()->nweight();
   
   //Set the Vector to hold local coordinates
   Vector<double> s(DIM);
   
   //Integers to store the local equation and unknown numbers
   int local_eqn=0, local_unknown=0;
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Assign values of s
     for(unsigned i=0;i<DIM;i++) s[i] = integral_pt()->knot(ipt,i);

     //Get the integral weight
     double w = integral_pt()->weight(ipt);
     
     //Call the derivatives of the shape and test functions
     double J = dshape_dbasis_and_dtest_eulerian_at_knot(ipt,
                                                         psi,dpsidx,
                                                         basis,dbasisdx,
                                                         test,dtestdx);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;

     //Calculate local values of the pressure and velocity components
     //Allocate and initialise to zero
     double interpolated_u=0.0;
     Vector<double> interpolated_x(DIM,0.0);
     Vector<double> interpolated_dudx(DIM,0.0);
   
     //Calculate function value and derivatives:
     //-----------------------------------------
     // Loop over nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       // Loop over directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_x[j] += nodal_position(l,j)*psi[l];
        }
      }


     for(unsigned l=0;l<n_basis;l++)
      {
       interpolated_u += u(l)*basis[l];
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_dudx[j] += u(l)*dbasisdx(l,j);
        }
      }

     //Get source function
     //-------------------
     double source;
     get_source_poisson(interpolated_x,source);

     //Get the diffusivity
     double diff;
     get_diff(interpolated_x,diff);
     
     // Assemble residuals and Jacobian
     //--------------------------------
       
     // Loop over the test functions
     for(unsigned l=0;l<n_basis;l++)
      {
       //Get the local equation
       local_eqn = Local_eqn[l];
       //IF it's not a boundary condition
       if(local_eqn >= 0)
        {
         // Add body force/source term here 
         residuals[local_eqn] += source*test[l]*W;
             
         // The Poisson bit itself
         for(unsigned k=0;k<DIM;k++)
          {
           residuals[local_eqn] += diff*interpolated_dudx[k]*dtestdx(l,k)*W;
          }

         // Calculate the jacobian
         //-----------------------
         if(flag)
          {
           //Loop over the velocity shape functions again
           for(unsigned l2=0;l2<n_basis;l2++)
            { 
             local_unknown = Local_eqn[l2];
             //If at a non-zero degree of freedom add in the entry
             if(local_unknown >= 0)
              {
               //Add contribution to Elemental Matrix
               for(unsigned i=0;i<DIM;i++)
                {
                 jacobian(local_eqn,local_unknown) 
                  += diff*dbasisdx(l2,i)*dtestdx(l,i)*W;
                }
              }
            }
          }
        }
      }
     
    } // End of loop over integration points
   
  }
 
 /// Pointer to source function:
 PoissonSourceFctPt Source_fct_pt;
 
 /// Pointer to the diffusivity function
 PoissonDiffFctPt Diff_fct_pt;

 /// Array to hold local eqn numbers: Local_eqn[n] (=-1 for BC) 
 Vector<int> Local_eqn;
 
  private:
 
};



//======================================================================
/// Refineable version of Poisson equations
//======================================================================
template <unsigned DIM>
class RefineableModalPoissonEquations :
 public virtual ModalPoissonEquations<DIM>,
 public virtual RefineableElement,
 public virtual ElementWithZ2ErrorEstimator
{
  public:

 /// Constructor, simply call other constructors
 RefineableModalPoissonEquations() : ModalPoissonEquations<DIM>(),
  RefineableElement(), ElementWithZ2ErrorEstimator() 
  { } 

 /// Broken copy constructor
 RefineableModalPoissonEquations(const RefineableModalPoissonEquations<DIM>&
                                       dummy) 
  { 
   BrokenCopy::broken_copy("RefineableModalPoissonEquations");
  } 
 
 virtual inline unsigned u_index_poisson() const {return 0;}
 
 /// Number of 'flux' terms for Z2 error estimation 
 unsigned num_Z2_flux_terms() {return DIM;}

 /// Get 'flux' for Z2 error recovery:  Standard flux.from Poisson equations
 void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
  {this->get_flux(s,flux,2);}

 /// Get error against and norm of exact flux
 //void compute_exact_Z2_error(
 // std::ostream &outfile,
 // FiniteElement::SteadyExactSolutionFctPt exact_flux_pt,
 // double& error, double& norm);
 
/// Get the function value u in Vector.
/// Note: Given the generality of the interface (this function
/// is usually called from black-box documentation or interpolation routines),
/// the values Vector sets its own size in here.
void get_interpolated_values(const Vector<double>&s,  Vector<double>& values)
 {
  // Set size of Vector: u
  values.resize(1);
  
  //Find number of nodes
  unsigned n_node = nnode();
  
  //Local shape function
  Shape psi(n_node);
  
  //Find values of shape function
  shape(s,psi);
  
  //Initialise value of u
  values[0] = 0.0;

  //Find the index at which the poisson unknown is stored
  unsigned u_nodal_index = this->u_index_poisson();
  
  //Loop over the local nodes and sum up the values
  for(unsigned l=0;l<n_node;l++)
   {
    values[0] += this->nodal_value(l,u_nodal_index)*psi[l];
   }
 }


 /// Get the function value u in Vector.
 /// Note: Given the generality of the interface (this function
 /// is usually called from black-box documentation or interpolation routines),
 /// the values Vector sets its own size in here.
 void get_interpolated_values(const unsigned& t, const Vector<double>&s, 
                              Vector<double>& values)
  {
   if (t!=0)
    {
     std::string error_message =
      "Time-dependent version of get_interpolated_values() ";
     error_message += "not implemented for this element \n";
     throw 
      OomphLibError(error_message,
                    "RefineablePoissonEquations::get_interpolated_values()",
                    OOMPH_EXCEPTION_LOCATION);
    }
   else
    {
     //Make sure that we call this particular object's steady 
     //get_interpolated_values (it could get overloaded lower down)
     this->get_interpolated_values(s,values);
    }
  }

 
 ///  Further build: Copy source function pointer from father element
 void further_build()
  {
   this->Source_fct_pt=dynamic_cast<RefineableModalPoissonEquations<DIM>*>
                                (this->father_element_pt())->source_fct_pt();
  }


  private:

 /// Add element's contribution to elemental residual vector and/or 
 /// Jacobian matrix 
 /// flag=1: compute both
 /// flag=0: compute only residual vector
 void fill_in_generic_residual_contribution(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  const unsigned& flag); 

 /// Compute derivatives of elemental residual vector with respect
 /// to nodal coordinates. Overwrites default implementation in 
 /// FiniteElement base class.
 /// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
 virtual void get_dresidual_dnodal_coordinates(RankThreeTensor<double>&
                                               dresidual_dnodal_coordinates);
  
 void get_source_gradient(unsigned int&, oomph::Vector<double>&, oomph::Vector<double>&)
  {
   std::string error_message = "get_source_gradient() is ";
   error_message += "not implemented for this element \n";
   throw 
    OomphLibError(error_message,
                  "RefineablePoissonEquations::get_source_gradient()",
                  OOMPH_EXCEPTION_LOCATION);
  }
};

//======================================================================
/// Compute derivatives of elemental residual vector with respect
/// to nodal coordinates. 
/// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
/// Overloads the FD-based version in the FE base class.
//======================================================================
template <unsigned DIM>
void RefineableModalPoissonEquations<DIM>::get_dresidual_dnodal_coordinates(
                       RankThreeTensor<double>&
                       dresidual_dnodal_coordinates)
{

 //Find out how many nodes there are
 unsigned n_node = nnode();
 //Find out the total number of basis and test functions there are
 unsigned n_basis = this->nbasis();
 
 //Set up memory for the shape and test functions
 Shape psi(n_node), basis(n_basis), test(n_basis);
 DShape dpsidx(n_node,DIM), dbasisdx(n_basis,DIM), dtestdx(n_basis,DIM);
 DShape dpsidx_pls(n_node,DIM), dbasisdx_pls(n_basis,DIM), dtestdx_pls(n_basis,DIM);

 // Get number of shape controlling nodes 
 unsigned n_shape_controlling_node=nshape_controlling_nodes();

 // Deriatives of shape fct derivatives w.r.t. nodal coords
 RankFourTensor<double> d_dpsidx_dX(DIM,n_shape_controlling_node,n_node,DIM);
 RankFourTensor<double> d_dtestdx_dX(DIM,n_shape_controlling_node,n_node,DIM);

 // Derivative of Jacobian of mapping w.r.t. to nodal coords
 DenseMatrix<double> dJ_dX(DIM,n_shape_controlling_node);

 // Derivatives of derivative of u w.r.t. nodal coords
 RankThreeTensor<double> d_dudx_dX(DIM,n_shape_controlling_node,DIM);

 // Gradient of source fct
 Vector<double> d_source_dx(DIM);

 //Index at which the poisson unknown is stored
 const unsigned u_nodal_index = this->u_index_poisson();
 
 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();

 //Integers to store the local equation number
 int local_eqn=0;

 // Local storage for pointers to hang_info object
 HangInfo *hang_info_pt=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Call the derivatives of the shape and test functions
   double J = this->dshape_dbasis_and_dtest_eulerian_at_knot(ipt, 
                                                             psi, 
                                                             dpsidx,
                                                             basis,
                                                             dbasisdx,
                                                             test, 
                                                             dtestdx);
       
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
     double u_value = nodal_value(l,u_nodal_index);
     // Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       interpolated_x[j] += nodal_position(l,j)*psi(l);
       interpolated_dudx[j] += u_value*dpsidx(l,j);
      }
    }

   //Get source function
   //-------------------
   double source;
   //this->get_source(ipt,interpolated_x,source);
   this->get_source_poisson(interpolated_x,source);

   // FD step 
   double eps_fd=GeneralisedElement::Default_fd_jacobian_step;
   
   std::map<Node*,unsigned> local_shape_controlling_node_lookup= 
    shape_controlling_node_lookup();

   // FD loop over shape-controlling nodes
   for (std::map<Node*,unsigned>::iterator it=
        local_shape_controlling_node_lookup.begin();
        it!=local_shape_controlling_node_lookup.end();
        it++)
    {  
     // Get node
     Node* nod_pt=it->first;
     
     // Get its number
     unsigned jj=it->second;
          
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
        this->dshape_dbasis_and_dtest_eulerian_at_knot(ipt, 
                                                       psi, 
                                                       dpsidx_pls,
                                                       basis,
                                                       dbasisdx_pls,
                                                       test, 
                                                       dtestdx_pls);
       
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
           aux+=nodal_value(j_nod,u_nodal_index)*
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
   this->get_source_gradient(ipt, interpolated_x, d_source_dx);

   
   // Assemble shape derivatives
   //---------------------------
       
   // Loop over the nodes for the test functions 
   for(unsigned l=0;l<n_node;l++)
    {
     //Local variables used to store the number of master nodes and the
     //weight associated with the shape function if the node is hanging
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
       //Get the local equation number and hang_weight
       //If the node is hanging
       if(is_node_hanging)
        {
         //Read out the local equation number from the m-th master node
         local_eqn =  this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                           u_nodal_index);
         
         //Read out the weight from the master node
         hang_weight = hang_info_pt->master_weight(m);
        }
       //If the node is not hanging
       else
        {
         //The local equation number comes from the node itself
         local_eqn = this->nodal_local_eqn(l,u_nodal_index);
         //The hang weight is one
         hang_weight = 1.0;
        }
       
       //If the nodal equation is not a boundary condition
       if(local_eqn >= 0)
        {
         // Loop over coordinate directions
         for (unsigned ii=0;ii<DIM;ii++)
          {              
           // Loop over shape controlling nodes
           for (unsigned jj=0;jj<n_shape_controlling_node;jj++)
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
             dresidual_dnodal_coordinates(local_eqn,ii,jj)+=sum*w*hang_weight;
            }
          }
        }
      }
    }
   
  } // End of loop over integration points
 
}




template<unsigned DIM>
class ModalPRefineableQElement : public PRefineableQElement<DIM,2>
{
public:
 ModalPRefineableQElement() : PRefineableQElement<DIM,2>() {}
 
 void initial_setup(Tree* const &adopted_father_pt=0, const unsigned &initial_p_order=0);
	void pre_build(Mesh*&, Vector<Node*>&);
 void further_build()
  {
   PRefineableQElement<DIM,2>::further_build(); //(Empty)
  }
 
 // p-refine the element
 void p_refine(const int &inc,
               Mesh* const &mesh_pt,
               GeneralisedElement* const &clone_pt);
 
 // Overload the shape and basis functions
 void shape(const Vector<double> &s, Shape &psi) const;
 void basis(const Vector<double> &s, Shape &basis) const;
 void dshape_local(const Vector<double> &s, Shape &psi, DShape &dpsi) const;
 void dbasis_local(const Vector<double> &s, Shape &basis, DShape &dbasis) const;
 void d2shape_local(const Vector<double> &s, Shape &psi, DShape &dpsids, DShape &d2psids) const;
 void d2basis_local(const Vector<double> &s, Shape &basis, DShape &dbasisds, DShape &d2basisds) const;
};

template<>
void ModalPRefineableQElement<1>::initial_setup(Tree* const &adopted_father_pt, const unsigned &initial_p_order)
{
 // Create storage for internal data
 if (this->ninternal_data()==0)
  {
   this->add_internal_data(new Data(0));
  }

 //Storage for pointer to my father (in binarytree impersonation)
 BinaryTree* father_pt;
 
 // Check if an adopted father has been specified
 if (adopted_father_pt!=0)
  {
   //Get pointer to my father (in binarytree impersonation)
   father_pt = dynamic_cast<BinaryTree*>(adopted_father_pt);
  }
 // Check if element is in a tree
 else if (Tree_pt!=0)
  {
   //Get pointer to my father (in binarytree impersonation)
   father_pt = dynamic_cast<BinaryTree*>(binary_tree_pt()->father_pt());
  }
 else
  {
   throw OomphLibError(
          "Element not in a tree, and no adopted father has been specified!",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
  }

 // Check if element has father
 if (father_pt!=0)
  {
   if (PRefineableQElement<1,2>* father_el_pt =
          dynamic_cast<PRefineableQElement<1,2>*>
          (father_pt->object_pt()))
    {
     unsigned father_p_order = father_el_pt->p_order();
     // Set the correct p-order of the element
     this->p_order() = father_p_order;
     
     // Now sort out the element...
     
     // Set integration scheme
     delete this->integral_pt();
     switch(this->p_order())
     {
     case 2:
      this->set_integration_scheme(new GaussLobattoLegendre<1,2>);
      break;
     case 3:
      this->set_integration_scheme(new GaussLobattoLegendre<1,3>);
      break;
     case 4:
      this->set_integration_scheme(new GaussLobattoLegendre<1,4>);
      break;
     case 5:
      this->set_integration_scheme(new GaussLobattoLegendre<1,5>);
      break;
     case 6:
      this->set_integration_scheme(new GaussLobattoLegendre<1,6>);
      break;
     case 7:
      this->set_integration_scheme(new GaussLobattoLegendre<1,7>);
      break;
     default:
      oomph_info << "\n ERROR: Exceeded maximum polynomial order for";
      oomph_info << "\n        integration scheme." << std::endl;
      break;
     }
     
     // Set the size of the internal data
     if (this->internal_data_pt(0)->nvalue() <= this->p_order()-this->nnode())
      {
       this->internal_data_pt(0)->resize(this->p_order()-this->nnode());
      }
     else
      {
       Data* new_data_pt = new Data(this->p_order()-this->nnode());
       delete internal_data_pt(0);
       internal_data_pt(0) = new_data_pt;
      }
     if(this->ninternal_data()==0) delete this->internal_data_pt(0);
     // Interpolate initial guess from father to modes
     /// / How to interpolate for the initial guess?
     // Set to zero for now, don't do projection problem
     for (unsigned i=0; i<this->p_order()-this->nnode_1d(); i++)
      {
       this->internal_data_pt(0)->set_value(i, 0.0);
      }
    }
   else
    {
        oomph_info << "Dynamic cast FAILURE :-(" << endl;
    }
  }
}

template<>
void ModalPRefineableQElement<1>::pre_build(Mesh*&, Vector<Node*>&) {}

template<>
void ModalPRefineableQElement<1>::p_refine(
      const int &inc,
      Mesh* const &mesh_pt,
      GeneralisedElement* const &clone_pt)
{
 // BENFLAG: In this case we do not need the pointer to a clone of the
 //          element, or to the mesh, but they are required in the
 //          interface because we are overloading this function from
 //          the class PRefineableElement in which they are present.
 
 // Create storage for modes if none exists
 if (this->ninternal_data()==0)
  {
   this->add_internal_data(new Data(0));
  }
 
 // Increment p-order of the element
 this->p_order() += inc;
 
 // Change integration scheme
 delete this->integral_pt();
 switch(this->p_order())
 {
 case 2:
  this->set_integration_scheme(new GaussLobattoLegendre<1,2>);
  break;
 case 3:
  this->set_integration_scheme(new GaussLobattoLegendre<1,3>);
  break;
 case 4:
  this->set_integration_scheme(new GaussLobattoLegendre<1,4>);
  break;
 case 5:
  this->set_integration_scheme(new GaussLobattoLegendre<1,5>);
  break;
 case 6:
  this->set_integration_scheme(new GaussLobattoLegendre<1,6>);
  break;
 case 7:
  this->set_integration_scheme(new GaussLobattoLegendre<1,7>);
  break;
 default:
  oomph_info << "\n ERROR: Exceeded maximum polynomial order for";
  oomph_info << "\n        integration scheme." << std::endl;
  break;
 }
 
 // Add/remove new modes
 if (this->internal_data_pt(0)->nvalue() <= this->p_order()-this->nnode())
  {
   this->internal_data_pt(0)->resize(this->p_order()-this->nnode());
  }
 else
  {
   Data* new_data_pt = new Data(this->p_order()-this->nnode());
   delete internal_data_pt(0);
   internal_data_pt(0) = new_data_pt;
  }
 if(this->ninternal_data()==0) delete this->internal_data_pt(0);
}

template<unsigned DIM>
void ModalPRefineableQElement<DIM>::shape(const Vector<double> &s, Shape &psi) const
{
 // Shape functions
 psi(0) = 0.5*(1.0 - s[0]);
 psi(1) = 0.5*(1.0 + s[0]);
}

template<unsigned DIM>
void ModalPRefineableQElement<DIM>::basis(const Vector<double> &s, Shape &basis) const
{
 // Get nnode-1d and p-order
 unsigned p_order = this->p_order();
 
 //Call the shape functions
 OneDimensionalModalShape psi1(p_order, s[0]);
 
 // Loop over shapes and copy across
 for(unsigned i=0;i<p_order;i++) 
  {
   basis(i) = psi1[i];
  }
}

template<unsigned DIM>
void ModalPRefineableQElement<DIM>::dshape_local(const Vector<double> &s, Shape &psi, DShape &dpsi) const
{
 // Shape functions
 psi(0) = 0.5*(1.0 - s[0]);
 psi(1) = 0.5*(1.0 + s[0]);
 
 // dShape/dlocal
 dpsi(0,0) = -0.5;
 dpsi(1,0) = 0.5;
}

template<unsigned DIM>
void ModalPRefineableQElement<DIM>::dbasis_local(const Vector<double> &s, Shape &basis, DShape &dbasis) const
{
 // Get nnode-1d and p-order
 unsigned p_order = this->p_order();
 
 //Call the shape functions and derivatives
 OneDimensionalModalShape psi1(p_order, s[0]);
 OneDimensionalModalDShape dpsi1ds(p_order, s[0]);
 
 // Loop over shapes and copy across
 for(unsigned i=0;i<p_order;i++) 
  {
   basis(i) = psi1[i];
   dbasis(i,0) = dpsi1ds[i];
  }
 
}

template<unsigned DIM>
void ModalPRefineableQElement<DIM>::d2shape_local(const Vector<double> &s, Shape &psi, 
                                         DShape &dpsids, DShape &d2psids) const
{
 std::ostringstream error_message;
 error_message <<"\nd2shape_local currently not implemented for this element\n";
 throw OomphLibError(error_message.str(),
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
}

template<unsigned DIM>
void ModalPRefineableQElement<DIM>::d2basis_local(const Vector<double> &s, Shape &basis, 
                                         DShape &dbasisds, DShape &d2basisds) const
{
 std::ostringstream error_message;
 error_message <<"\nd2basis_local currently not implemented for this element\n";
 throw OomphLibError(error_message.str(),
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
}





template<unsigned DIM>
class ModalPRefineableQPoissonElement :
 public QPoissonElement<DIM,2>,
 public virtual RefineableModalPoissonEquations<DIM>,
 public virtual ModalPRefineableQElement<DIM>
{
   public:

 /// Constructor, simply call the other constructors 
 ModalPRefineableQPoissonElement() : 
  RefineableElement(),
  RefineableModalPoissonEquations<DIM>(),
  ModalPRefineableQElement<DIM>(),
  QPoissonElement<DIM,2>()
   {
    // Set integration scheme
    // (To avoid memory leaks in pre-build and p-refine where new integration schemes are created)
    this->set_integration_scheme(new GaussLobattoLegendre<DIM,2>);
   }


 /// Broken copy constructor
 ModalPRefineableQPoissonElement(const ModalPRefineableQPoissonElement<DIM>& 
                           dummy) 
  { 
   BrokenCopy::broken_copy("ModalPRefineableQPoissonElement");
  } 
 
 // Overload nnode_1d() const
 unsigned nnode_1d() const {return 2;}
 
 virtual void further_build();
 
 /// Number of continuously interpolated values: 1
 unsigned ncont_interpolated_values() const {return 1;}

 /// Number of vertex nodes in the element
 unsigned nvertex_node() const
  {return QPoissonElement<DIM,2>::nvertex_node();}

 /// Pointer to the j-th vertex node in the element
 Node* vertex_node_pt(const unsigned& j) const
  {return QPoissonElement<DIM,2>::vertex_node_pt(j);}

 /// Order of recovery shape functions for Z2 error estimation:
 /// - Same order as shape functions.
 //unsigned nrecovery_order()
 // {
 //  if(this->nnode_1d() < 4) {return (this->nnode_1d()-1);}
 //  else {return 3;}
 // }
 /// - Constant recovery order, since recovery order of the first element
 ///   is used for the whole mesh.
 unsigned nrecovery_order() {return 3;}
 
 void basis(const Vector<double> &s, Shape &basis) const
  {ModalPRefineableQElement<DIM>::basis(s, basis);}
 
 void dbasis_local(const Vector<double> &s, Shape &basis, DShape &dbasis) const
  {ModalPRefineableQElement<DIM>::dbasis_local(s, basis, dbasis);}
 
 void d2basis_local(const Vector<double> &s, Shape &basis, DShape &dbasisds,
                               DShape &d2basisds) const
  {ModalPRefineableQElement<DIM>::d2basis_local(s, basis, dbasisds, d2basisds);}
 
 //Extra modal stuff:

 /// Function pointer to source function fct(x,f(x)) -- 
 /// x is a Vector! 
 typedef void (*PoissonSourceFctPt)(const Vector<double>& x, double& f);
 
 inline double u(const unsigned& n) const
  {
   if(n < this->nnode())
    {
     return this->nodal_value(n,0);
    }
   else
    {
     return this->internal_data_pt(0)->value(n-this->nnode());
    }
  }
 
 virtual unsigned nbasis() const
  {return this->p_order();}
 
 void output(ostream &outfile)
  {RefineableModalPoissonEquations<DIM>::output(outfile);}
 
 void output(ostream &outfile, const unsigned &nplot)
  {RefineableModalPoissonEquations<DIM>::output(outfile, nplot);}
 
 void output_fct(ostream &outfile, const unsigned &nplot, 
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {RefineableModalPoissonEquations<DIM>::output_fct(outfile, nplot,
                 exact_soln_pt);}
 
 void output_fct(ostream &outfile, const unsigned &nplot,
                         const double& time, 
                         FiniteElement::UnsteadyExactSolutionFctPt 
                         exact_soln_pt)
  {RefineableModalPoissonEquations<DIM>::output_fct(outfile, nplot, time,
                         exact_soln_pt);}
 
 void compute_error(ostream &outfile, 
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm)
  {RefineableModalPoissonEquations<DIM>::compute_error(outfile,
                    exact_soln_pt, error, norm);}
 
 /// Get error against and norm of exact solution
 void compute_energy_error(ostream &outfile, 
                    FiniteElement::SteadyExactSolutionFctPt exact_grad_pt,
                    double& error, double& norm)
  {

   // Initialise
   error=0.0;
   norm=0.0;

   //Vector of local coordinates
   Vector<double> s(DIM);

   // Vector for coordintes
   Vector<double> x(DIM);
   
   //Set the value of n_intpt
   unsigned n_intpt = this->integral_pt()->nweight();
   
   // Setup output structure: Conversion is fishy but it's only output...
   unsigned nplot;
   if (DIM==1)
    {
     nplot=n_intpt;
    }
   else 
    {
     nplot=unsigned(pow(n_intpt,1.0/double(DIM)));
    }
   
   // Tecplot header info
   outfile << this->tecplot_zone_string(nplot);
   
   // Exact solution Vector (here a scalar)
   Vector<double> exact_grad(1);
   
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
     
     // Get FE du/dx
     Vector<double> dudx_fe(1);
     ModalPoissonEquations<1>::get_flux(s,dudx_fe,this->p_order());
     
     // Get exact solution at this point
     (*exact_grad_pt)(x,exact_grad);
     
     //Output x,y,...,error
     for(unsigned i=0;i<DIM;i++)
      {
       outfile << x[i] << " ";
      }
     outfile << exact_grad[0] << " " << exact_grad[0]-dudx_fe[0] << std::endl;
     
     // Add to error and norm
     norm+=exact_grad[0]*exact_grad[0]*W;
     error+=(exact_grad[0]-dudx_fe[0])*(exact_grad[0]-dudx_fe[0])*W;
     
    }
  }
 
 void compute_error(ostream &outfile, 
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time, double& error, double& norm)
  {RefineableModalPoissonEquations<DIM>::compute_error(outfile,
                    exact_soln_pt, time, error, norm);}
 
 PoissonSourceFctPt& source_fct_pt() 
  {return RefineableModalPoissonEquations<DIM>::source_fct_pt();}
 
 PoissonSourceFctPt source_fct_pt() const
  {return RefineableModalPoissonEquations<DIM>::source_fct_pt();}
 
 unsigned self_test()
 {return RefineableModalPoissonEquations<DIM>::self_test();}
 
//=========================================================================
/// Compute the geometric basis functions and also
/// first derivatives w.r.t. global coordinates at local coordinate s;
/// Returns Jacobian of mapping from global to local coordinates.
/// Most general form of the function, but may be over-loaded, if desired
//=========================================================================
 virtual double dbasis_eulerian(const Vector<double> &s,
                                Shape &basis,
                                DShape &dbasis) const;
 
//========================================================================
/// Compute the geometric shape functions and also first
/// derivatives w.r.t. global coordinates at integration point ipt.
/// Most general form of function, but may be over-loaded if desired
//========================================================================
 virtual double dbasis_eulerian_at_knot(const unsigned &ipt,
                                              Shape &basis, 
                                              DShape &dbasis) const;

//=========================================================================
/// Return the shape function and its derivatives w.r.t. the local
/// coordinates at the ipt-th integration point.
//=========================================================================
 virtual void dbasis_local_at_knot(const unsigned &ipt, Shape &basis,
                                          DShape &dbasisds) const
 {
  //Find the dimension of the element
  const unsigned el_dim = this->dim();
  //Storage for the local coordinates of the integration point
  Vector<double> s(el_dim); 
  //Set the local coordinate
  for(unsigned i=0;i<el_dim;i++) {s[i] = this->integral_pt()->knot(ipt,i);}
  //Get the shape function and derivatives
  this->dbasis_local(s,basis,dbasisds);
 }
 
 virtual double dshape_dbasis_and_dtest_eulerian(const Vector<double> &s, 
                                                 Shape &psi, 
                                                 DShape &dpsidx, 
                                                 Shape &basis,
                                                 DShape &dbasisdx,
                                                 Shape &test, 
                                                 DShape &dtestdx) const
  {
   //Call the geometrical basis functions and derivatives
   double J = dbasis_eulerian(s,basis,dbasisdx);
   
   //Loop over the nodes and set the first nodal shapes equal to the basis functions
   for(unsigned i=0;i<this->nnode_1d();i++)
    {
     psi[i] = basis[i];
     dpsidx(i,0) = dbasisdx(i,0);
    }
   
   //Loop over the test functions and derivatives and set them equal to the
   //basis functions
   for(unsigned i=0;i<this->p_order();i++)
    {
     test[i] = basis[i]; 
     dtestdx(i,0) = dbasisdx(i,0);
    }
  
   //Return the jacobian
   return J;
  }
 
 virtual double dshape_dbasis_and_dtest_eulerian_at_knot(const unsigned &ipt, 
                                                         Shape &psi, 
                                                         DShape &dpsidx,
                                                         Shape &basis,
                                                         DShape &dbasisdx,
                                                         Shape &test, 
                                                         DShape &dtestdx) const
 {
  //Call the geometrical shape functions and derivatives  
  double J = dbasis_eulerian_at_knot(ipt,basis,dbasisdx);
  
  Vector<double> s(1);
  s[0] = this->integral_pt()->knot(ipt,0);
  
  //Loop over the nodes and set the first nodal shapes equal to the basis functions
  for(unsigned i=0;i<this->nnode_1d();i++)
   {
    psi[i] = basis[i];
    dpsidx(i,0) = dbasisdx(i,0);
   }
  
  //Loop over the test functions and derivatives and set them equal to the
  //basis functions
  for(unsigned i=0;i<this->p_order();i++)
   {
    test[i] = basis[i]; 
    dtestdx(i,0) = dbasisdx(i,0);
   }
  
  //Return the jacobian
  return J;
 }
 
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {RefineableModalPoissonEquations<DIM>::
      fill_in_contribution_to_residuals(residuals);}
 
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
  {RefineableModalPoissonEquations<DIM>::
      fill_in_contribution_to_jacobian(residuals, jacobian);}
 
 void fill_in_generic_residual_contribution(Vector<double> &residuals,
                                            DenseMatrix<double> &jacobian, 
                                            const unsigned& flag)
  {RefineableModalPoissonEquations<DIM>::
      fill_in_generic_residual_contribution(residuals, jacobian, flag);}
 
 virtual void get_dresidual_dnodal_coordinates(RankThreeTensor<double>&
                                               dresidual_dnodal_coordinates)
  {ModalPoissonEquations<DIM>::
      get_dresidual_dnodal_coordinates(dresidual_dnodal_coordinates);}
 
 void assign_additional_local_eqn_numbers();
};

//=========================================================================
/// Compute the geometric basis functions and also
/// first derivatives w.r.t. global coordinates at local coordinate s;
/// Returns Jacobian of mapping from global to local coordinates.
/// Most general form of the function, but may be over-loaded, if desired
//=========================================================================
template<unsigned DIM>
double ModalPRefineableQPoissonElement<DIM>::dbasis_eulerian(
                                      const Vector<double> &s, 
                                      Shape &basis,
                                      DShape &dbasis) const
{
 //Find the element dimension
 const unsigned el_dim = this->dim();

 //Get the values of the shape functions and their local derivatives
 //Temporarily stored in dpsi
 this->dbasis_local(s,basis,dbasis);
 
 //Allocate memory for the inverse jacobian
 DenseMatrix<double> inverse_jacobian(el_dim);
 //Now calculate the inverse jacobian
 const double det = this->local_to_eulerian_mapping(dbasis,inverse_jacobian);
 
 //Now set the values of the derivatives to be dpsidx
 this->transform_derivatives(inverse_jacobian,dbasis);
 //Return the determinant of the jacobian
 
 return det;
}

//========================================================================
/// Compute the geometric shape functions and also first
/// derivatives w.r.t. global coordinates at integration point ipt.
/// Most general form of function, but may be over-loaded if desired
//========================================================================
template<unsigned DIM>
double ModalPRefineableQPoissonElement<DIM>::dbasis_eulerian_at_knot(const unsigned &ipt,
                                              Shape &basis, 
                                              DShape &dbasis) const
{
 //Find the element dimension
 const unsigned el_dim = this->dim();

 //Get the values of the shape function and local derivatives
 //Temporarily store it in dpsi
 this->dbasis_local_at_knot(ipt,basis,dbasis);

 //Allocate memory for the inverse jacobian
 DenseMatrix<double> inverse_jacobian(el_dim);
 //Now calculate the inverse jacobian
 const double det = this->local_to_eulerian_mapping(dbasis,inverse_jacobian);

 //Now set the values of the derivatives to dpsidx
 this->transform_derivatives(inverse_jacobian,dbasis);
 //Return the determinant of the jacobian
 
 return det;
}

//======================================================================
/// Setup the local equation numbering schemes:
///
///                     /  local equation number.
///   Local_eqn[n] = |
///                     \  -1 if boundary condition.
///
/// Pure version without hanging nodes
//======================================================================
template<unsigned DIM>
void ModalPRefineableQPoissonElement<DIM>::assign_additional_local_eqn_numbers()
{
 //Resize the equation counters
 this->Local_eqn.resize(this->p_order());
 
 //Loop over the nodes
 for(unsigned i=0;i<this->nnode_1d();i++)
  {
   //Set the local equation number to be the first value stored at the node
   this->Local_eqn[i] = this->nodal_local_eqn(i,0);
  }
 //Continue over the internal data
 for(unsigned i=this->nnode_1d();i<this->p_order();i++)
  {
   //Set the local equation number to be the ith value of the internal data
   this->Local_eqn[i] = this->internal_local_eqn(0,i-this->nnode_1d());
  }
}

template<unsigned DIM>
void ModalPRefineableQPoissonElement<DIM>::further_build()
{
 // Do the ModalPRefineableQElement further build
 ModalPRefineableQElement<DIM>::further_build();
 // Needed to set the source function pointer
 RefineableModalPoissonEquations<DIM>::further_build();
}



//======start_of_problem_class============================================
/// 1D Poisson problem discretised with refineable 1D QPoisson elements.
/// The specific type of element is specified via the template parameter.
//========================================================================
template<class ELEMENT> 
class PRefineableOneDPoissonProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source function
 PRefineableOneDPoissonProblem(PoissonEquations<1>::PoissonSourceFctPt 
                               source_fct_pt);
 
 /// Destructor (empty)
 ~PRefineableOneDPoissonProblem()
  {
   delete mesh_pt()->spatial_error_estimator_pt();
   delete Problem::mesh_pt();
  }
 
 /// Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();
 
 /// Update the problem after solve (empty)
 void actions_after_newton_solve()
  {}
 
 /// Doc the solution. DocInfo object stores flags/labels for where
 /// the output gets written to.
 void doc_solution(DocInfo& doc_info);
 
 /// Overloaded version of the Problem's access function to the mesh.
 /// Recasts the pointer to the base Mesh object to the actual mesh type.
 RefineableOneDMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableOneDMesh<ELEMENT>*>(Problem::mesh_pt());
  }
 
private:
 
 /// Pointer to source function
 PoissonEquations<1>::PoissonSourceFctPt Source_fct_pt;
 
}; // End of problem class



//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
PRefineableOneDPoissonProblem<ELEMENT>::
PRefineableOneDPoissonProblem(PoissonEquations<1>::PoissonSourceFctPt 
                              source_fct_pt)
 : Source_fct_pt(source_fct_pt)
{ 
 
 // Set up mesh
 // -----------

 // Number of elements
 const unsigned n = 2;

 // Domain length
 const double length = 1.0;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableOneDMesh<ELEMENT>(n,length);

 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt() = new Z2ErrorEstimator;
  
 // Set the boundary conditions for this problem. All nodes are free by
 // default so only need to pin the ones that have Dirichlet conditions here.
 const unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
    }
  }

 // Complete build of all elements so they are fully functional

 // Loop over the elements to set up element-specific things that cannot
 // be handled by the (argument-free!) ELEMENT constructor: Pass pointer
 // to source function
 const unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;
  }

 // Set up equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

} // End of constructor




//=====start_of_actions_before_newton_solve===============================
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void PRefineableOneDPoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 // Determine the number of mesh boundaries
 const unsigned n_boundary = mesh_pt()->nboundary();
 
 // Loop over these boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine the number of nodes on this boundary b
   const unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
   
   // Loop over these nodes
   for (unsigned n=0;n<n_boundary_node;n++)
    {
     // Get pointer to node
     Node* nod_pt = mesh_pt()->boundary_node_pt(b,n);
     
     // Extract nodal coordinates from node:
     Vector<double> x(1);
     x[0] = nod_pt->x(0);

     // Compute the value of the exact solution at the nodal point
     Vector<double> u(1);
     ArcTanSolnForPoisson::get_exact_u(x,u);

     // Assign the value to the one (and only) nodal value at this node
     nod_pt->set_value(0,u[0]);
    }
  }
}  // End of actions_before_newton_solve




//=====start_of_doc_solution==============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void PRefineableOneDPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 
 // Declare output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 const unsigned npts = 5;

 // Output solution 
 // ---------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
 // Output exact solution 
 // ---------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,ArcTanSolnForPoisson::get_exact_u); 
 some_file.close();

 // Doc error and return of the square of the L2 error
 // --------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,ArcTanSolnForPoisson::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

} // End of doc_solution



//======start_of_main=====================================================
/// Driver code for 1D Poisson problem
//========================================================================
int main()
{
 // Create label for output
 // -----------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;

 // Choose a large value for the steepness of the "step"
 ArcTanSolnForPoisson::Alpha=100.0;
 
 
 cout << "\n====================================\nNodal elements:\n" << endl;
 
 //Set up the problem
 //------------------
 
 // Create the problem with 1D three-node refineable elements from the
 // PRefineableQPoissonElement<1> family. Pass pointer to source function. 
 PRefineableOneDPoissonProblem<PRefineableQPoissonElement<1> > 
  nodal_problem(&ArcTanSolnForPoisson::get_source);
 
 // Check if we're ready to go:
 // ---------------------------
 cout << "\n\n\nProblem self-test ";
 if (nodal_problem.self_test()==0) 
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("Self test failed",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Refine problem uniformly 2 times
 for(unsigned i=0;i<2;i++)
  {
   cout << "p_refining:" << endl;
   nodal_problem.p_refine_uniformly();
  }
 for(unsigned i=0;i<2;i++)
  {
   cout << "h_refining:" << endl;
   nodal_problem.refine_uniformly();
  }
 
 doc_info.number()=2;
 nodal_problem.newton_solve();
 nodal_problem.doc_solution(doc_info);

 nodal_problem.adapt();
 doc_info.number()=3;
 nodal_problem.newton_solve();
 nodal_problem.doc_solution(doc_info);
 
 nodal_problem.p_adapt();
 doc_info.number()=4;
 nodal_problem.newton_solve();
 nodal_problem.doc_solution(doc_info);
 
 nodal_problem.adapt();
 doc_info.number()=5;
 nodal_problem.newton_solve();
 nodal_problem.doc_solution(doc_info);
 
 nodal_problem.p_adapt();
 doc_info.number()=6;
 nodal_problem.newton_solve();
 nodal_problem.doc_solution(doc_info);
 
 nodal_problem.adapt();
 doc_info.number()=7;
 nodal_problem.newton_solve();
 nodal_problem.doc_solution(doc_info);
 
 nodal_problem.p_adapt();
 doc_info.number()=8;
 nodal_problem.newton_solve();
 nodal_problem.doc_solution(doc_info);
 
 doc_info.number()=0;
 nodal_problem.doc_solution(doc_info);
 
 cout << "\n====================================\nModal elements:\n" << endl;
 
 //Set up the problem
 //------------------
 
 // Create the problem with 1D three-node refineable elements from the
 // RefineableLinePoissonElement family. Pass pointer to source function. 
 PRefineableOneDPoissonProblem<ModalPRefineableQPoissonElement<1> > 
  modal_problem(&ArcTanSolnForPoisson::get_source);
 
 // Check if we're ready to go:
 // ---------------------------
 cout << "\n\n\nProblem self-test ";
 if (modal_problem.self_test()==0) 
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("Self test failed",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Refine problem uniformly 2 times
 for(unsigned i=0;i<2;i++)
  {
   cout << "p_refining:" << endl;
   modal_problem.p_refine_uniformly();
  }
 for(unsigned i=0;i<2;i++)
  {
   cout << "h_refining:" << endl;
   modal_problem.refine_uniformly();
  }
 
 doc_info.number()=2;
 modal_problem.newton_solve();
 modal_problem.doc_solution(doc_info);
 
 modal_problem.adapt();
 doc_info.number()=3;
 modal_problem.newton_solve();
 modal_problem.doc_solution(doc_info);
 
 modal_problem.p_adapt();
 doc_info.number()=4;
 modal_problem.newton_solve();
 modal_problem.doc_solution(doc_info);
 
 modal_problem.adapt();
 doc_info.number()=5;
 modal_problem.newton_solve();
 modal_problem.doc_solution(doc_info);
 
 modal_problem.p_adapt();
 doc_info.number()=6;
 modal_problem.newton_solve();
 modal_problem.doc_solution(doc_info);

 modal_problem.adapt();
 doc_info.number()=7;
 modal_problem.newton_solve();
 modal_problem.doc_solution(doc_info);
 
 modal_problem.p_adapt();
 doc_info.number()=8;
 modal_problem.newton_solve();
 modal_problem.doc_solution(doc_info);
 
 doc_info.number()=1;
 modal_problem.doc_solution(doc_info);
 
} // End of main
