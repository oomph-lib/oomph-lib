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
/// Diver to solve the orr-sommerfeld equation for Poiseuille flow
/// using the default eigensolver.

// Generic oomph-lib routines only
#include "generic.h"

// Include the 1D mesh
#include "meshes/one_d_mesh.h"

//Use standard namespace
using namespace std;

using namespace oomph;

//Global Parameters
namespace Param
{
 double Re=6100.0;
 double A_imag = 0.0;
 unsigned N_element=50;
};

//================================================================
/// A class for all elements that solve the simple orr-sommerfeld
/// equations for a given base flow. The formulation is in primitive
/// variables and so contains quite a large pressure NULL space. 
/// The equations are derived from the Navier--Stokes equations by
/// posing a perturbation to an axially uniform base-flow of the form:
/// \f[
///  u = U(z) + \hat{u}, \quad v = \hat{v}, \quad w = \hat{w}
/// \f]
/// where the perturbation  is of the form
/// \f[
///  \hat{u} = \mbox{e}^{\lambda t + \mbox{i}\alpha x}
/// \f]
/// 
/// The equations are split into real and imaginary parts and is
/// currently only implemented for the two-dimensional equations.
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
template <unsigned DIM>
class OrrSommerfeldEquations : public virtual FiniteElement
{
public:
 /// Constructor
 OrrSommerfeldEquations() {}
 
 /// Access function: Nodal function value at local node n
 /// Uses suitably interpolated value for hanging nodes.
 virtual double u(const unsigned& n, const unsigned &i) const = 0;
 
 virtual double p(const unsigned &n, const unsigned &i) const = 0;

 virtual unsigned npres() const=0;

 virtual void pshape(const Vector<double> &s, Shape &psi) const=0;

 void get_base_flow(const Vector<double> &s, double &U, double &Uz)
  {
   //Get the coordinate
   double z = interpolated_x(s,0);
   //Return the value of the base flow
   U = (1.0-z*z);
   Uz = -2.0*z;
  }

 //Access function for the Reynolds number
 double* &re_pt() {return Re_pt;}

 double* &a_real_pt() {return A_Real_pt;}
 
 double* &a_imag_pt() {return A_Imag_pt;}

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
     
     for(unsigned i=0;i<4;i++)
      {
       outfile << interpolated_u(s,i) << " ";
       }

     outfile << interpolated_p(s,0) << " " 
             << interpolated_p(s,1) << " ";
     outfile << std::endl;
    }

   // Write tecplot footer (e.g. FE connectivity lists)
   write_tecplot_zone_footer(outfile,nplot);
  }

 /// Assemble the contributions to the jacobian and mass matrices
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals,
  DenseMatrix<double> &jacobian, DenseMatrix<double> &mass_matrix)
  {
   //Find out how many nodes there are
   unsigned n_node = nnode();
   //Number of pressure degrees of freedom
   unsigned n_pres = npres();

   //Set up memory for the shape functions
   Shape psi(n_node);
   DShape dpsidx(n_node,DIM);
   Shape psip(n_pres);

   //Set the value of n_intpt
   unsigned n_intpt = integral_pt()->nweight();
   
   //Set the Vector to hold local coordinates
   Vector<double> s(DIM);
   
   //Integers to store the local equation and unknown numbers
   int local_eqn=0, local_unknown=0;

   //Storage for the base flow
   double Ubar=0.0, Ubardz=0.0;

   //Get the Reynolds number
   double Re = *Re_pt;
   double A_real = *A_Real_pt;
   double A_imag = *A_Imag_pt;

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Assign values of s
     for(unsigned i=0;i<DIM;i++) s[i] = integral_pt()->knot(ipt,i);

     //Get the integral weight
     double w = integral_pt()->weight(ipt);
     
     //Call the derivatives of the shape and test functions
     double J = dshape_eulerian_at_knot(ipt,
                                        psi,dpsidx);

     //Get the pressure shape functions
     pshape(s,psip);

     //Get the base velocity field
     get_base_flow(s,Ubar,Ubardz);
     
     //Calculate the real and imaginary parts of our variable B
     double B_real = (A_real*A_real - A_imag*A_imag)/Re - A_imag*Ubar;
     double B_imag = 2.0*A_real*A_imag/Re + A_real*Ubar;

     //Premultiply the weights and the Jacobian
     double W = w*J;

     //Assemble the contributions to the matrices
     
     // Loop over the test functions
     for(unsigned l=0;l<n_node;l++)
      {
       //REAL PART OF THE U EQUATION
       //Get the local equation
       local_eqn = U_local_eqn(l,0);
       /*IF it's not a boundary condition*/
       if(local_eqn >= 0)
        {
         //Loop over the shape functions
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           //Get the contribution from the real part of u
           local_unknown = U_local_eqn(l2,0);
           //If at a non-zero degree of freedom add in the entry
           if(local_unknown >= 0)
            {
             //Mass matrix
             mass_matrix(local_eqn,local_unknown) += psi(l)*psi(l2)*W;
             //Jacobian terms
             jacobian(local_eqn,local_unknown) -= 
              B_real*psi(l2)*psi(l)*W;
             jacobian(local_eqn,local_unknown) -= 
              (1.0/Re)*dpsidx(l2,0)*dpsidx(l,0)*W;
            }

           //Get the contribution from the imaginary part of u
           local_unknown = U_local_eqn(l2,1);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += B_imag*psi(l2)*psi(l)*W;
            }

           //Get the contribution from the real part of w
           local_unknown = U_local_eqn(l2,2);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -= Ubardz*psi(l2)*psi(l)*W;
            }

          }
         
         //Loop over the pressure shape functions
         for(unsigned l2=0;l2<n_pres;l2++)
          {
           //Get the real pressure contribution to the equation
           local_unknown = P_local_eqn(l2,0);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += 
              A_imag*psip(l2)*psi(l)*W;
            }
           //Get the imaginary pressure contribution to the equation
           local_unknown = P_local_eqn(l2,1);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += 
              A_real*psip(l2)*psi(l)*W;
            }
          }
        }

       //IMAGINARY PART OF THE U EQUATION
       //Get the local equation
       local_eqn = U_local_eqn(l,1);
       /*IF it's not a boundary condition*/
       if(local_eqn >= 0)
        {
         //Loop over the shape functions
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           //Get the contribution from the real part of u
           local_unknown = U_local_eqn(l2,0);
           //If at a non-zero degree of freedom add in the entry
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -= B_imag*psi(l2)*psi(l)*W;
            }

           //Get the contribution from the imaginary part of u
           local_unknown = U_local_eqn(l2,1);
           if(local_unknown >= 0)
            {
             //Mass matrix
             mass_matrix(local_eqn,local_unknown) += psi(l)*psi(l2)*W;
             //Jacobian terms
             jacobian(local_eqn,local_unknown) -= B_real*psi(l2)*psi(l)*W;
             jacobian(local_eqn,local_unknown) -=
              (1.0/Re)*dpsidx(l2,0)*dpsidx(l,0)*W;
            }

           //Get the contibution from the imaginary part of w
           local_unknown = U_local_eqn(l2,3);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -= Ubardz*psi(l2)*psi(l)*W;
            }
          }

         //Loop over the pressure shape functions
         for(unsigned l2=0;l2<n_pres;l2++)
          {
           //Get the real pressure contribution to the equation
           local_unknown = P_local_eqn(l2,0);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -= 
              A_real*psip(l2)*psi(l)*W;
            }
           //Get the imaginary pressure contribution to the equation
           local_unknown = P_local_eqn(l2,1);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              A_imag*psip(l2)*psi(l)*W;
            }
          }
        }

       
       //REAL PART OF THE W EQUATION
       //Get the local equation
       local_eqn = U_local_eqn(l,2);
       /*IF it's not a boundary condition*/
       if(local_eqn >= 0)
        {
         //Loop over the shape functions
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           //Get the contribution from the real part of w
           local_unknown = U_local_eqn(l2,2);
           //If at a non-zero degree of freedom add in the entry
           if(local_unknown >= 0)
            {
             //Mass matrix
             mass_matrix(local_eqn,local_unknown) += psi(l)*psi(l2)*W;
             //Jacobian terms
             jacobian(local_eqn,local_unknown) -= B_real*psi(l2)*psi(l)*W;
             jacobian(local_eqn,local_unknown) -= 
              (1.0/Re)*dpsidx(l2,0)*dpsidx(l,0)*W;
            }

           //Get the contribution from the imaginary part of w
           local_unknown = U_local_eqn(l2,3);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += B_imag*psi(l2)*psi(l)*W;
            }

          } //End of loop over the shape functions
         

         //Loop over the pressure shape functions
         for(unsigned l2=0;l2<n_pres;l2++)
          {
           //Get the real pressure contribution to the equation
           local_unknown = P_local_eqn(l2,0);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += 
              psip(l2)*dpsidx(l,0)*W;
            }
          }
        }

       //IMAGINARY PART OF THE W EQUATION
       //Get the local equation
       local_eqn = U_local_eqn(l,3);
       /*IF it's not a boundary condition*/
       if(local_eqn >= 0)
        {
         //Loop over the shape functions
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           //Get the contribution from the real part of w
           local_unknown = U_local_eqn(l2,2);
           //If at a non-zero degree of freedom add in the entry
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -= B_imag*psi(l2)*psi(l)*W;
            }

           //Get the contribution from the imaginary part of w
           local_unknown = U_local_eqn(l2,3);
           if(local_unknown >= 0)
            {
             //Mass matrix
             mass_matrix(local_eqn,local_unknown) += psi(l)*psi(l2)*W;
             //Jacobian terms
             jacobian(local_eqn,local_unknown) -= B_real*psi(l2)*psi(l)*W;
             jacobian(local_eqn,local_unknown) -= 
              (1.0/Re)*dpsidx(l2,0)*dpsidx(l,0)*W;
            }
          }

         //Loop over the pressure shape functions
         for(unsigned l2=0;l2<n_pres;l2++)
          {
           //Get the imaginary pressure contribution to the equation
           local_unknown = P_local_eqn(l2,1);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) +=
              psip(l2)*dpsidx(l,0)*W;
            }
          }
        }
      }

     //Loop over the pressure shape functions
     for(unsigned l=0;l<n_pres;l++)
      {
       //Real part of the continuity equation
       local_eqn = P_local_eqn(l,0);
       if(local_eqn >= 0)
        {
         //Loop over the shape functions
         for(unsigned l2=0;l2<n_node;l2++)
          {
           //REAL U conributions
           local_unknown = U_local_eqn(l2,0);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=  
              A_imag*psi(l2)*psip(l)*W;
            }

           //IMAGINARY U contributions
           local_unknown = U_local_eqn(l2,1);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=
              A_real*psi(l2)*psip(l)*W;
            }

           //REAL W contribution
           local_unknown = U_local_eqn(l2,2);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += dpsidx(l2,0)*psip(l)*W;
            }
          }
        }


       //Imaginary part of the continuity equation
       local_eqn = P_local_eqn(l,1);
       if(local_eqn >= 0)
        {
         //Loop over the shape functions
         for(unsigned l2=0;l2<n_node;l2++)
          {
           //REAL U conributions
           local_unknown = U_local_eqn(l2,0);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += 
              A_real*psi(l2)*psip(l)*W;
            }

           //IMAGINARY U contributions
           local_unknown = U_local_eqn(l2,1);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -=
              A_imag*psi(l2)*psip(l)*W;
            }

           //IMAGINARY W contribution
           local_unknown = U_local_eqn(l2,3);
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += dpsidx(l2,0)*psip(l)*W;
            }
          }
        }
      } //End of loop over pressure shape functions
    }
  }

 /// Return FE representation of the pressure 
 inline double interpolated_p(const Vector<double> &s, const unsigned &i) const
  {
   unsigned n_pres = npres();

   //Pressure shape functions
   Shape psip(n_pres);

   //Find values of shape function
   this->pshape(s,psip);

   //Initialise value of p
   double interpolated_p = 0.0;

   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_pres;l++) 
    {
     interpolated_p+=p(l,i)*psip[l];
    }

   return(interpolated_p);
  }



 /// Return FE representation of function value u(s) at local coordinate s
 inline double interpolated_u(const Vector<double> &s, const unsigned &i) const
  {
   unsigned n_node = nnode();

   //Local shape functions
   Shape psi(n_node);

   //Find values of basis function
   this->shape(s,psi);

   //Initialise value of u
   double interpolated_u = 0.0;

   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_node;l++) 
    {
     interpolated_u+=u(l,i)*psi[l];
    }

   return(interpolated_u);
  }

protected:

 /// Shape/test functions and derivs w.r.t. to global coords at 
 /// local coord. s; return  Jacobian of mapping
 virtual double dshape_eulerian(const Vector<double> &s, 
                                Shape &psi, 
                                DShape &dpsidx) const=0;

 /// Shape/test functions and derivs w.r.t. to global coords at 
 /// integration point ipt; return  Jacobian of mapping
 virtual double dshape_eulerian_at_knot(const unsigned &ipt, 
                                        Shape &psi, 
                                        DShape &dpsidx) const=0;
 
 /// Array to hold local eqn numbers: Local_eqn[n] (=-1 for BC) 
 DenseMatrix<int> U_local_eqn;
 
 //Array to hold the pressure equations
 DenseMatrix<int> P_local_eqn;

private:

 double *Re_pt;
 double *A_Real_pt;
 double *A_Imag_pt;
 
};



/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////



//======================================================================
/// QXPoissonElement elements are linear/quadrilateral/brick-shaped 
/// Poisson elements with isoparametric interpolation for the function.
///
/// Empty, just establishes the template parameters
///
///
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
class QOrrSommerfeldElement
{

private:

 // Note: this is just in here to stop doxygen from barking.
 // This is a policy class and we provide specialised versions
 // for all dimensions below, so this one never actually gets built

 /// Static array of ints to hold number of variables at 
 /// nodes: Initial_Nvalue[n]
 static const unsigned Initial_Nvalue[];
 
};

//======================================================================
/// QXPoissonElement<1,NNODE_1D> elements are 1D Poisson elements with 
/// NNODE_1D nodal points in each coordinate direction. Inherits
/// from QElement and PoissonEquations
///
//======================================================================
template <unsigned NNODE_1D>
class QOrrSommerfeldElement<1,NNODE_1D> : public virtual QElement<1,NNODE_1D>, 
 public OrrSommerfeldEquations<1>
{
 
  public:


 /// Function to assign the local equation numbers arrays.
 /// Can be overloaded for hangeable version of element
 void assign_additional_local_eqn_numbers();


 /// Constructor: Call constructors for QElement and 
 /// Poisson equations
 QOrrSommerfeldElement() : QElement<1,NNODE_1D>(), OrrSommerfeldEquations<1>()
  {
   //Allocate two data each of size two to store the complex
   //pressure perturbation field.
   for(unsigned i=0;i<2;i++) {add_internal_data(new Data(2));}
  }



 ///  Required  # of `values' (pinned or dofs) 
 /// at node n
 inline unsigned required_nvalue(const unsigned &n) const 
  {
   return Initial_Nvalue[n];
  }



 /// Access functions for the nodal function values: Nodal 
 /// value at local node n.
 /// Uses suitably interpolated value for hanging nodes.
 inline double u(const unsigned &n, const unsigned &i) const
  {return nodal_value(n,i);}

 /// Return value of p
 inline double p(const unsigned &n, const unsigned &i) const
  {return this->internal_data_pt(n)->value(i);}

 /// Output function:  
 ///  x,y,u   or    x,y,z,u
 void output(ostream &outfile)
  {
   OrrSommerfeldEquations<1>::output(outfile);
  }

  unsigned npres() const {return ninternal_data();}


 ///  Output function:  
 ///   x,y,u   or    x,y,z,u at Nplot^DIM plot points
 void output(ostream &outfile, const unsigned &Nplot)
  {
   OrrSommerfeldEquations<1>::output(outfile,Nplot);
  }


 //Set the pshape
 void pshape(const Vector<double> &s, Shape &psi) const
  {
   //psi[0] = 0.5*(1.0 -s[0]);
   //psi[1] = 0.5*(1.0 + s[0]);
   psi[0] = 1.0;
   psi[1] = s[0];
  }

protected:

/// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double dshape_eulerian(const Vector<double> &s, 
                               Shape &psi, 
                               DShape &dpsidx) const;
 

 /// Shape, test functions & derivs. w.r.t. to global coords. at
 /// integration point ipt. Return Jacobian.
 inline double dshape_eulerian_at_knot(const unsigned& ipt,
                                       Shape &psi, 
                                       DShape &dpsidx) const;
private:
 
 /// Static array of ints to hold number of variables
 /// at nodes: Initial_Nvalue[n]
 static const unsigned Initial_Nvalue[];
 

};




//Inline functions:


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned NNODE_1D>
double QOrrSommerfeldElement<1,NNODE_1D>::
dshape_eulerian(const Vector<double> &s,
                Shape &psi, 
                DShape &dpsidx) const
{
 //Call the geometrical shape functions and derivatives  
 double J = QElement<1,NNODE_1D>::dshape_eulerian(s,psi,dpsidx);

 //Return the jacobian
 return J;
}


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned NNODE_1D>
double QOrrSommerfeldElement<1,NNODE_1D>::dshape_eulerian_at_knot(
 const unsigned &ipt,
 Shape &psi, 
 DShape &dpsidx) const
{
 //Call the geometrical shape functions and derivatives  
 double J = QElement<1,NNODE_1D>::dshape_eulerian_at_knot(ipt,psi,dpsidx);

 //Return the jacobian
 return J;
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
template<unsigned NNODE_1D>
void QOrrSommerfeldElement<1,NNODE_1D>::assign_additional_local_eqn_numbers()
{
 //Resize the equation counters
 U_local_eqn.resize(NNODE_1D,4);
 
 //Loop over the nodes
 for(unsigned n=0;n<NNODE_1D;n++)
  {
   for(unsigned i=0;i<4;i++)
    {
     //Set the local equation number to be the first value stored at the node
     U_local_eqn(n,i) = nodal_local_eqn(n,i);
    }
  }

 unsigned n_pres = npres();
 //Resize the pressure equation counter
 P_local_eqn.resize(n_pres,2);
 for(unsigned n=0;n<n_pres;n++)
  {
   for(unsigned i=0;i<2;i++)
    {
     P_local_eqn(n,i) = internal_local_eqn(n,i);
    }
  }

}


//======================================================================
// Set the data for the number of Variables at each node
//======================================================================
template<>
const unsigned QOrrSommerfeldElement<1,4>::Initial_Nvalue[4]={4,4,4,4};
template<>
const unsigned QOrrSommerfeldElement<1,3>::Initial_Nvalue[3]={4,4,4};
template<>
const unsigned QOrrSommerfeldElement<1,2>::Initial_Nvalue[2]={4,4};


//==start_of_problem_class============================================
/// OrrSommerfeld problem in unit interval.
//====================================================================
template<class ELEMENT> 
class OrrSommerfeldProblem : public Problem
{
private:
 double *A_real_pt;

 EigenSolver *eigensolver_pt;

private:
 Gauss<1,4> integral;

public:

 /// Constructor: Pass number of elements and pointer to source function
 OrrSommerfeldProblem(const unsigned& n_element, double *const &A_Real_pt);

 /// Destructor (empty)
 ~OrrSommerfeldProblem(){}
 
 /// Access function for the Real part of the wavenumber
 double* &a_real_pt() {return A_real_pt;}

 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve() {}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution, pass the number of the case considered,
 /// so that output files can be distinguished.
 void doc_solution(const unsigned& label);

}; // end of problem class





//=====start_of_constructor===============================================
/// Constructor for 1D Poisson problem in unit interval.
/// Discretise the 1D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function. 
//========================================================================
template<class ELEMENT>
OrrSommerfeldProblem<ELEMENT>::OrrSommerfeldProblem(const unsigned& n_element,
                                                double* const &A_REAL_PT)
 : A_real_pt(A_REAL_PT)
{ 
 // hierher Andrew: switched this to anasazi because QZ dies (I assume the
 // one of the two matrices doesn't have the required properties.
 eigen_solver_pt()=new LAPACK_QZ;
  
 //Set the shift to be zero (the default)
 eigen_solver_pt()->set_shift(0.0);

 //Set the minimum arc-length to be quite large, so that the system bails
 //quickly if there are problems
 Minimum_ds = 100.0;

 // Build mesh and store pointer in Problem
 Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element,-1.0,1.0);
 
 // Set the boundary conditions for this problem: By default, all nodal
 // values are free -- we only need to pin the ones that have 
 // Dirichlet conditions. 

 //Loop over the real and imaginary parts of the velocity
 for(unsigned i=0;i<4;i++)
  {
   // Pin the nodal values at the single node on mesh 
   // boundary 0 (= the left domain boundary at x=0)
   mesh_pt()->boundary_node_pt(0,0)->pin(i);
   
   // Pin the nodal values at the single node on mesh 
   // boundary 1 (= the right domain boundary at x=1)
   mesh_pt()->boundary_node_pt(1,0)->pin(i);
  }

 //We don't need to pin the pressure, I'm not quite sure why
 //this is. I suspect that it's because by NOT integrating the 
 //pressure by parts in the x-momentum equation, I have lost that
 //freedom. Don't really understand though. The freedom has been lost
 //because of the travelling-wave-type expansion of p which means that
 //when differentiated with respect to x the absolute value of p
 //matters!
 
 //mesh_pt()->element_pt(0)->internal_data_pt(0)->pin(0);
 //mesh_pt()->element_pt(0)->internal_data_pt(0)->pin(1);

 //Set up the elements
 for(unsigned e=0;e<n_element;e++)
  {
   ELEMENT* temp_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   temp_pt->re_pt() = &Param::Re;
   temp_pt->a_real_pt() = A_real_pt;
   temp_pt->a_imag_pt() = &Param::A_imag;
   //temp_pt->set_integration_scheme(&integral);
  }

 // Setup equation numbering scheme
 cout << assign_eqn_numbers() << " equations in eigenproblem" << std::endl;

} // end of constructor



//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT>
void OrrSommerfeldProblem<ELEMENT>::doc_solution(const unsigned& label)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution with specified number of plot points per element
 sprintf(filename,"soln%i.dat",label);
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

} // end of doc

 

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

/// Mini Element to implement the eigenvalue constraint
class ZeroResidualGenElement : public GeneralisedElement
{
 //This Element actually stores a pointer to the eigenvalue problem
 OrrSommerfeldProblem<QOrrSommerfeldElement<1,3> > *eigenvalue_problem_pt;
 
public:
 /// Simple constructor
 ZeroResidualGenElement()
  {
   //One item of internal data will store the real part of the eigenvalue
   add_internal_data(new Data(1));
   
   //Construct the eigenvalue problem
   eigenvalue_problem_pt = new 
    OrrSommerfeldProblem<QOrrSommerfeldElement<1,3> >
    (Param::N_element,internal_data_pt(0)->value_pt(0));
  }

 /// Destructor
 virtual ~ZeroResidualGenElement() 
  {delete eigenvalue_problem_pt;}

 /// Access function for the the eigenvalue
 double& eval() {return *(internal_data_pt(0)->value_pt(0));}

 double* eval_pt() {return internal_data_pt(0)->value_pt(0);}

 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //There is only one value
   unsigned n_dof = ndof();
   //If we're pinned return
   if(n_dof==0) return;
   
   //Solve the eigenvalue problem
   Vector<std::complex<double>> alpha;
   Vector<double> beta;
   eigenvalue_problem_pt->solve_eigenproblem(4,alpha,beta);
   
   //We need to make sure that we pick the eigenvalue with the largest real 
   //part that is not infinite
   unsigned index=0;
   double e_real = 0.0;
   const unsigned n_eval = alpha.size();
   //Find the first non-infinite value
   for(unsigned i=0;i<n_eval;++i)
    {
     if(beta[i]!=0.0)
      {
       index = i; 
       e_real = alpha[i].real()/beta[i];
       break;
      }
    }
   

   unsigned start = index;
   //Loop over the eigenvalues and find the largest (non-infinite)  magnitude
   for(unsigned i=start;i<n_eval;i++)
    {
     //Check we are not infinite
     if(beta[i] != 0.0)
      {
       //Real part
       double new_e_real = alpha[i].real()/beta[i];
       if(new_e_real > e_real) {e_real = new_e_real; index=i;}
      }
    }
   //Now set the residuals to be the eigenvalue with largest real part
   residuals[0] = alpha[index].real()/beta[index];
  }

 //Get the residuals
 void get_residuals(Vector<double> &residuals)
  {
   fill_in_contribution_to_residuals(residuals);
  }

 void output(std::ostream &output, const unsigned &)
  {
   eigenvalue_problem_pt->mesh_pt()->output(output,5);
  }

};


/// Mini problem class for the continuation problem
class ContinuationProblem : public Problem
{

 ZeroResidualGenElement *element_pt;

public:

 ContinuationProblem(const unsigned &n_element)
  {
   Param::N_element = n_element;
   //Create our single element
   element_pt = new ZeroResidualGenElement;

   //Create a mesh
   mesh_pt() = new Mesh;

   //Now add a single element to our mesh
   mesh_pt()->add_element_pt(element_pt);
   //Set an initial value for the parameter
   element_pt->eval() = 1.1;
   //Assign the equations
   cout << assign_eqn_numbers() << std::endl;
  }

 void doc_solution()
  {
   cout << element_pt->eval() << std::endl;
  }
   void actions_before_newton_solve() {}
   void actions_after_newton_solve() {}
   

 void solve()
  {
   ofstream trace("neutral.dat");
   //ofstream output_file;
   //char filename[100];

   //Desired_newton_iterations_ds = 3;
   //Newton_solver_tolerance = 1.0e-6;
   double ds = -100.0;

   //Do four Continuation steps
   //More can be done, but there appears to be a difference between
   //results between 64- and 32-bit machines. The step size is
   //then different enough that the validation tests fail, even though
   //exactly the same results are computed
   for(unsigned i=0;i<3;i++)
    {
     //Param::Re += 100.0;
     //Lets take an arc-length step
     ds = arc_length_step_solve(&Param::Re,ds);
     //newton_solve();

     //Output to a trace file
     trace << Param::Re << " " << element_pt->eval() << std::endl;

     //sprintf(filename,"output.%g.dat",Param::Re);
     //output_file.open(filename);
     //element_pt->output(output_file,5);
     //output_file.close();
    }
  }
};



//======start_of_main==================================================
/// Driver for 1D Poisson problem
//=====================================================================
int main(int argc, char* argv[])
{
 
//Want to test Trilinos if we have it, so we must initialise MPI
//if we have compiled with it
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 
 ContinuationProblem problem(50);

 problem.solve();

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

} // end of main









