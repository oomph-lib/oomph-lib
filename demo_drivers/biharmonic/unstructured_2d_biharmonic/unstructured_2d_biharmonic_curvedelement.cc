//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
//Driver for an unstructuredmesh 2D-Biharmonic problem with the curved elements

// Generic oomph-lib routines
#include "generic.h"

// Include the mesh
#include "meshes/one_d_mesh.h"
#include "meshes/simple_rectangular_tri_mesh.h"
#include "meshes/simple_rectangular_quadmesh.h"
#include "generic/geom_objects.h"
#include "meshes/triangle_mesh.h"

using namespace std;

using namespace oomph;


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//Header file for C1Curved biharmonic elements
#ifndef OOMPH_BIHARMONIC_CURVED_ELEMENTS_HEADER
#define OOMPH_BIHARMONIC_CURVED_ELEMENTS_HEADER

namespace oomph
{

//=============================================================
/// A class for all subparametric elements that solve the 
/// Biharmonic equations.
/// \f[ 
/// \frac{\partial^4 u}{\partial x_i^4} = f(x_j)
/// \f] 
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
template <unsigned DIM, unsigned NNODE_1D>
class MyBiharmonicEquations : public virtual C1CurvedElement<DIM,NNODE_1D>
{

public:
 
 /// \short Function pointer to source function fct(x,f(x)) -- 
 /// x is a Vector! 
 typedef void (*SourceFctPt)(const Vector<double>& x, double& f);


 /// \short Function pointer to gradient of source function  fct(x,g(x)) -- 
 /// x is a Vector! 
 typedef void (*SourceFctGradientPt)(const Vector<double>& x, 
                                            Vector<double>& gradient);

 /// \short Function pointer to the exact solutions -- 
 /// x is a Vector! 
 typedef void (*ExactSolnPt)(const Vector<double>& x, Vector<double>& exact_soln);

 /// Constructor (must initialise the Source_fct_pt to null)
 MyBiharmonicEquations() : Source_fct_pt(0), Source_fct_gradient_pt(0)
  {}
 
 /// Broken copy constructor
 MyBiharmonicEquations(const MyBiharmonicEquations& dummy) 
  { 
   BrokenCopy::broken_copy("MyBiharmonicEquations");
  } 
 
 /// Broken assignment operator
 void operator=(const MyBiharmonicEquations&) 
  {
   BrokenCopy::broken_assign("MyBiharmonicEquations");
  }

 /// \short Return the index at which the unknown value
 /// is stored. 
 /// In derived multi-physics elements, this function should be overloaded
 /// to reflect the chosen storage scheme. Note that these equations require
 /// that the unknown is always stored at the same index at each node.
  virtual inline unsigned u_index() const {return this->required_nvalue(0);}
 
 /// Output with default number of plot points
 void output(std::ostream &outfile) 
  {
   const unsigned n_plot=5;
   output(outfile,n_plot);
  }

 /// \short Output FE representation of soln: x,y,u or x,y,z,u at 
 /// n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot);

 /// C_style output with default number of plot points
 void output(FILE* file_pt)
  {
   const unsigned n_plot=5;
   output(file_pt,n_plot);
  }

 /// \short C-style output FE representation of soln: x,y,u or x,y,z,u at 
 /// n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot);

 /// Output exact soln: x,y,u_exact or x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot, 
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

 /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at 
 /// n_plot^DIM plot points (dummy time-dependent version to 
 /// keep intel compiler happy)
 virtual void output_fct(std::ostream &outfile, const unsigned &n_plot,
                         const double& time, 
                         FiniteElement::UnsteadyExactSolutionFctPt 
                         exact_soln_pt)
  {
   throw OomphLibError(
    "There is no time-dependent output_fct() for these elements ",
    "MyShellEquations<DIM>::output_fct()",
    OOMPH_EXCEPTION_LOCATION);
  }


 /// Get error against and norm of exact solution
 void compute_error(std::ostream &outfile, 
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm);


 /// Dummy, time dependent error checker
 void compute_error(std::ostream &outfile, 
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time, double& error, double& norm)
  {
   throw OomphLibError(
    "There is no time-dependent compute_error() for these elements",
    "MyBiharmonicEquations<DIM>::compute_error()",
    OOMPH_EXCEPTION_LOCATION);
  }

 /// Access function: Pointer to source function
 ExactSolnPt& exact_pt() {return Exact_soln_pt;}

 /// Access function: Pointer to source function. Const version
 ExactSolnPt exact_pt() const {return Exact_soln_pt;}

 /// Access function: Pointer to source function
 SourceFctPt& source_fct_pt() {return Source_fct_pt;}

 /// Access function: Pointer to source function. Const version
 SourceFctPt source_fct_pt() const {return Source_fct_pt;}

 /// Access function: Pointer to gradient of source function
 SourceFctGradientPt& source_fct_gradient_pt() 
  {return Source_fct_gradient_pt;}

 /// Access function: Pointer to gradient source function. Const version
 SourceFctGradientPt source_fct_gradient_pt() const 
  {return Source_fct_gradient_pt;}

 /// Get exact solutions at (Eulerian) position x.
 inline virtual void get_exact_solution(const Vector<double>& x,
                                        Vector<double>& exact_soln) const
  {
   (*Exact_soln_pt)(x,exact_soln);
  }

 
 inline virtual void get_source_function(const unsigned& ipt,
                                         const Vector<double>& x,
                                        double& source) const
  {
   //If no source function has been set, return zero
   if(Source_fct_pt==0) 
    {
     source = 0.0;
    }
   else
    {
     // Get source strength
     (*Source_fct_pt)(x,source);
    }
  }

 /// Get flux: flux[i] = du/dx_i
 void get_flux(const Vector<double>& s, Vector<double>& flux) const
  {}


 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_biharmonic(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }

 
 /// Add the element's contribution to its residual vector and 
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_biharmonic(residuals,jacobian,1);
  }


 /// \short Return FE representation of function value 
 /// at local coordinate s
 inline Vector<double> interpolated_u(const Vector<double> &s) 
  {
   //Find number of nodes
   const unsigned n_node = this->nnode();
   
   //Find number of position dofs
   const unsigned n_position_type = this->nnodal_position_type();
   
   //Get the index at which the unknown is stored
   const unsigned u_nodal_index = u_index();
   
   //Local c1-shape funtion
   Shape psi1(3,n_position_type),psi(n_node);
   DShape dpsidxi1(3,n_position_type,DIM),dpsidxi(3,DIM);
   DShape d2psidxi1(3,n_position_type,3);
   Vector<double> xi(DIM,0.0);
   
   //Initialise value of u
   Vector<double> interpolated_u(u_nodal_index,0.0);

   Vector<double> u(1,0.0);
   DenseMatrix<double> interpolated_dudxi(1,2,0.0), interpolated_d2udxi(1,3,0.0);
   DenseMatrix<double> M(21,36,0.0);
   this->my_interpolated_u_normal1(s,u,interpolated_dudxi,interpolated_d2udxi,M);
   
   interpolated_u[0] = u[0];
   interpolated_u[1] = interpolated_dudxi(0,0);
   interpolated_u[2] = interpolated_dudxi(0,1);
   interpolated_u[3] = interpolated_d2udxi(0,0);
   interpolated_u[4] = interpolated_d2udxi(0,1);
   interpolated_u[5] = interpolated_d2udxi(0,2);

   return(interpolated_u);
  }

 /// \short Self-test: Return 0 for OK
 unsigned self_test();


protected:

 /// \short Shape/test functions and derivs w.r.t. to global coords at 
 /// local coord. s; return  Jacobian of mapping
 virtual double d2shape_and_d2test_eulerian_biharmonic(const Vector<double> &s, 
                                                  Shape &psi, 
                                                  DShape &dpsidx, DShape &d2psidx, Shape &test, 
                                                  DShape &dtestdx, DShape &d2testdx) const=0;
 virtual double dshape_and_dtest_eulerian_biharmonic(const Vector<double> &s,
                                                  Shape &psi,
                                                  DShape &dpsidx, Shape &test,
                                                  DShape &dtestdx) const=0;

 /// \short Shape/test functions and derivs w.r.t. to global coords at 
 /// integration point ipt; return  Jacobian of mapping
 virtual double d2shape_and_d2test_eulerian_at_knot_biharmonic(const unsigned &ipt, 
                                                          Shape &psi, 
                                                          DShape &dpsidx,
                                                          DShape &d2psidx,
                                                          Shape &test, 
                                                          DShape &dtestdx,
                                                          DShape &d2testdx) 
  const=0;

 virtual double dshape_and_dtest_eulerian_at_knot_biharmonic(const unsigned &ipt,
                                                          Shape &psi,
                                                          DShape &dpsidx,
                                                          Shape &test,
                                                          DShape &dtestdx)
  const=0;

 /// \short Compute element residual Vector only (if flag=and/or element 
 /// Jacobian matrix 
 virtual void fill_in_generic_residual_contribution_biharmonic(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  const unsigned& flag); 
 
 /// Pointer to source function:
 SourceFctPt Source_fct_pt;

 /// Pointer to the exact solution:
 ExactSolnPt Exact_soln_pt;

 /// Pointer to gradient of source function
 SourceFctGradientPt Source_fct_gradient_pt;

};


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

//======================================================================
/// BiharmonicCurvedElement elements are a subparametric scheme with 
/// a Lagrange interpolation for approximating the geometry and 
/// the C1-functions for approximating variables.
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
class BiharmonicCurvedElement : public virtual MyBiharmonicEquations<DIM,NNODE_1D>
{

private:

 /// \short Static int that holds the number of variables at 
 /// nodes: always the same
 static const unsigned Initial_Nvalue;
 
  public:


 ///\short  Constructor: Call constructors for C1CurvedElement and 
 /// Biharmonic equations
 BiharmonicCurvedElement() : MyBiharmonicEquations<DIM,NNODE_1D>()
  {}
 
 /// Broken copy constructor
 BiharmonicCurvedElement(const BiharmonicCurvedElement<DIM,NNODE_1D>& dummy) 
  { 
   BrokenCopy::broken_copy("BiharmonicCurvedElement");
  } 
 
 /// Broken assignment operator
 void operator=(const BiharmonicCurvedElement<DIM,NNODE_1D>&) 
  {
   BrokenCopy::broken_assign("BiharmonicCurvedElement");
  }


 /// \short  Required  # of `values' (pinned or dofs) 
 /// at node n
 inline unsigned required_nvalue(const unsigned &n) const 
  {return Initial_Nvalue;}

 /// \short Output function:  
 ///  x,y,u   or    x,y,z,u
 void output(std::ostream &outfile)
  {MyBiharmonicEquations<DIM,NNODE_1D>::output(outfile);}


 ///  \short Output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot)
  {MyBiharmonicEquations<DIM,NNODE_1D>::output(outfile,n_plot);}


 /// \short C-style output function:  
 ///  x,y,u   or    x,y,z,u
 void output(FILE* file_pt)
  {MyBiharmonicEquations<DIM,NNODE_1D>::output(file_pt);}


 ///  \short C-style output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot)
  {MyBiharmonicEquations<DIM,NNODE_1D>::output(file_pt,n_plot);}


 /// \short Output function for an exact solution:
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {MyBiharmonicEquations<DIM,NNODE_1D>::output_fct(outfile,n_plot,exact_soln_pt);}



 /// \short Output function for a time-dependent exact solution.
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 /// (Calls the steady version)
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {MyBiharmonicEquations<DIM,NNODE_1D>::output_fct(outfile,n_plot,time,exact_soln_pt);}

 
protected:

/// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double d2shape_and_d2test_eulerian_biharmonic(
  const Vector<double> &s, Shape &psi, DShape &dpsidx, DShape &d2psidx, 
  Shape &test, DShape &dtestdx, DShape &d2testdx) const;

 inline double dshape_and_dtest_eulerian_biharmonic(
  const Vector<double> &s, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const;


 /// \short Shape, test functions & derivs. w.r.t. to global coords. at
 /// integration point ipt. Return Jacobian.
 inline double d2shape_and_d2test_eulerian_at_knot_biharmonic(const unsigned& ipt,
                                                         Shape &psi, 
                                                         DShape &dpsidx, 
                                                         DShape &d2psidx,
                                                         Shape &test,
                                                         DShape &dtestdx,
                                                         DShape &d2testdx) 
  const;

 inline double dshape_and_dtest_eulerian_at_knot_biharmonic(const unsigned &ipt,
                                                         Shape &psi,
                                                         DShape &dpsidx,
                                                         Shape &test,
                                                         DShape &dtestdx)
  const;

};


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


//=======================================================================
/// Face geometry for the CurvedTElement elements: The spatial 
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<BiharmonicCurvedElement<DIM,NNODE_1D> >: 
 public virtual TElement<DIM-1,NNODE_1D>
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional TElement
 FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}

};



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


//Inline functions:


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
 double BiharmonicCurvedElement<DIM,NNODE_1D>::dshape_and_dtest_eulerian_biharmonic(
  const Vector<double> &s,
  Shape &psi,
  DShape &dpsidx,
  Shape &test,
  DShape &dtestdx) const
{
 const double J = this->dbasis_c0_eulerian(s,psi,dpsidx);
 test = psi;
 dtestdx = dpsidx;

 return J;
}

template<unsigned DIM, unsigned NNODE_1D>
 double BiharmonicCurvedElement<DIM,NNODE_1D>::d2shape_and_d2test_eulerian_biharmonic(
  const Vector<double> &s,
  Shape &psi, 
  DShape &dpsidx,
  DShape &d2psidx,
  Shape &test, 
  DShape &dtestdx,
  DShape &d2testdx) const
{
 //Call the geometrical shape functions and derivatives  
 const double J = this->d2basis_c0_eulerian(s,psi,dpsidx,d2psidx);

 //Set the test functions equal to the shape functions
 test = psi;
 dtestdx= dpsidx;
 d2testdx = d2psidx;
 
 //Return the jacobian
 return J;
}

//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
double BiharmonicCurvedElement<DIM,NNODE_1D>::
 dshape_and_dtest_eulerian_at_knot_biharmonic(
  const unsigned &ipt,
  Shape &psi,
  DShape &dpsidx,
  Shape &test,
  DShape &dtestdx) const
{
 const double J = this->dshape_and_dtest_eulerian_at_knot(ipt,psi,dpsidx,test,dtestdx);

 return J;
}

template<unsigned DIM, unsigned NNODE_1D>
double BiharmonicCurvedElement<DIM,NNODE_1D>::
 d2shape_and_d2test_eulerian_at_knot_biharmonic(
  const unsigned &ipt,
  Shape &psi, 
  DShape &dpsidx,
  DShape &d2psidx,
  Shape &test, 
  DShape &dtestdx,
  DShape &d2testdx) const
{
 //Call the geometrical shape functions and derivatives  
 const double J = this->d2shape_and_d2test_eulerian_at_knot(ipt,psi,dpsidx,d2psidx,test,dtestdx,d2testdx);
 
 //Return the jacobian
 return J;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


}

#endif




namespace oomph
{
//======================================================================
/// Set the data for the number of Variables at each node
//======================================================================
 template<unsigned DIM, unsigned NNODE_1D>
 const unsigned BiharmonicCurvedElement<DIM,NNODE_1D>::Initial_Nvalue = 6;

//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void  MyBiharmonicEquations<DIM,NNODE_1D>::
fill_in_generic_residual_contribution_biharmonic(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              const unsigned& flag) 
{
 //Find out how many nodes there are
 const unsigned n_node = this->nnode();
 //Find out how many nodes positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();
  
 //Set up memory for the shape and test functions
 Shape psi(n_node,n_position_type), test(n_node);
 DShape dpsidxi(n_node,n_position_type,DIM), dtestdxi(n_node,DIM);
 DShape d2psidxi(n_node,n_position_type,3), d2testdxi(n_node,3);
 
 //Set the value of n_intpt
 const unsigned n_intpt = this->integral_pt()->nweight();
 
 //Integers to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;
      
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);
   
   //Call the derivatives of the shape and test functions for the unknown
   double J = d2shape_and_d2test_eulerian_at_knot_biharmonic(ipt,psi,dpsidxi,d2psidxi,test,dtestdxi,d2testdxi);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   //Calculate values of unknown
   Vector<double> interpolated_u(1,0.0);
   DenseMatrix<double> interpolated_dudxi(1,DIM,0.0);
   DenseMatrix<double> interpolated_d2udxi(1,3,0.0);
   
   //Allocate and initialise to zero
   Vector<double> interpolated_x(DIM,0.0);
   
   Vector<double>  s(2);
   s[0] = this->integral_pt()->knot(ipt,0);
   s[1] = this->integral_pt()->knot(ipt,1);
   this->my_interpolated_x(s,interpolated_x);
  
   //Calculate function value and derivatives
   //-----------------------------------------
   Vector<double> tt(2), xi(2);
   Vector<double> exact_soln(6);
   DenseMatrix<double> M(21,36,0.0);
   this->my_interpolated_u_normal1(s,interpolated_u,interpolated_dudxi,interpolated_d2udxi,M);

   //Get source function
   //-------------------
   double source;
   get_source_function(ipt,interpolated_x,source);
   
   //Check whether this element is a boundary element or not
   unsigned bd_element = this->is_boundary_element();

   //This is the case that it is NOT a boundary element. We can use 
   //the existing bell shape functions as a test function.
   if(bd_element==0)
    {
     // Loop over the test functions
     for(unsigned l=0;l<3;l++)
      {
       for(unsigned k=0;k<n_position_type;k++)
        {
         //Get the local equation
         local_eqn = this->nodal_local_eqn(l,k);
         //IF it's not a boundary condition
         if(local_eqn >= 0)
          {
           // Add body force/source term here 
           residuals[local_eqn] -= (source*psi(l,k)*W);
           residuals[local_eqn] += (interpolated_d2udxi(0,0)*d2psidxi(l,k,0)*W + interpolated_d2udxi(0,1)*d2psidxi(l,k,1)*W + 2.0*interpolated_d2udxi(0,2)*d2psidxi(l,k,2)*W);
           
           // Calculate the jacobian
           if(flag)
            {
             //Loop over the test functions again
             for(unsigned l2=0;l2<3;l2++)
              { 
               // Loop over position dofs
               for(unsigned k2=0;k2<n_position_type;k2++)
                {
                 local_unknown = this->nodal_local_eqn(l2,k2);
                 //If at a non-zero degree of freedom add in the entry
                 if(local_unknown >= 0)
                  {
                   //Add contribution to Elemental Matrix
                   jacobian(local_eqn,local_unknown) += d2psidxi(l2,k2,0)*d2psidxi(l,k,0)*W + d2psidxi(l2,k2,1)*d2psidxi(l,k,1)*W + 2.0*d2psidxi(l2,k2,2)*d2psidxi(l,k,2)*W;
                   
                  } 
                } 
              } 
            } // End of flag
          }
        }
      }
    }

   //This is the case that IT IS A CURVED-BOUNDARY ELEMENT. 
   if(bd_element==1)
    {
     // Assign the nodal local equation numbers
     DenseMatrix<int> nodal_eqn(6,n_position_type,-1);
     for(unsigned k=0;k<n_position_type;k++)
      {
       nodal_eqn(0,k) = this->nodal_local_eqn(0,k);
       nodal_eqn(1,k) = this->nodal_local_eqn(1,k);
       nodal_eqn(2,k) = this->nodal_local_eqn(2,k);
       nodal_eqn(3,k) = this->nodal_local_eqn(3,k);
       nodal_eqn(4,k) = this->nodal_local_eqn(4,k);
       nodal_eqn(5,k) = this->nodal_local_eqn(5,k);
      }            
    
     //Assign storage for the local shape function of curved elements
     //The shape functions are complete polynomials of degree 7 with 36 dofs
     Shape psi_curve(36);
     DShape dpsi_curve(36,2), d2psi_curve(36,3);
     this->d2basis_eulerian_curve(s,psi_curve,dpsi_curve,d2psi_curve); 

     //The following loops is to find the residuals.
     // Loop over the test functions
     // 1) which associate with values at interior nodes ei
     for(unsigned p=0;p<3;p++)
      {
       unsigned l=0,row=0;
       if(p==0){l=3;row=18;}
       else if(p==1){l=4;row=19;}
       else if(p==2){l=5;row=20;}
       
       //Get the local equation
       local_eqn = nodal_eqn(l,0);
       //IF it's not a boundary condition
       if(local_eqn >= 0)
        {
         for(unsigned j=0;j<36;j++)
          {
           if(M(row,j)!=0)
            {
             // Add body force/source term here 
             residuals[local_eqn] -= (source*psi_curve[j]*W)*M(row,j);
             
             residuals[local_eqn] += (interpolated_d2udxi(0,0)*d2psi_curve(j,0) + interpolated_d2udxi(0,1)*d2psi_curve(j,1) + 2.0*interpolated_d2udxi(0,2)*d2psi_curve(j,2))*W*M(row,j);
             
             // Calculate the jacobian
             //-----------------------
             if(flag)
              {
               //Loop over test functions associated with internal nodes
               for(unsigned p2=0;p2<3;p2++)
                {
                 unsigned l2=0,row2=0;
                 if(p2==0){l2=3;row2=18;}
                 else if(p2==1){l2=4;row2=19;}
                 else if(p2==2){l2=5;row2=20;}
                 //Get the local equation
                 local_unknown = nodal_eqn(l2,0);
                 //IF it's not a boundary condition
                 if(local_unknown >= 0)
                  {
                   for(unsigned j2=0;j2<36;j2++)
                    {
                     if(M(row2,j2)!=0)
                      {
                       //Add contribution to Elemental Matrix
                       jacobian(local_eqn,local_unknown) += (d2psi_curve(j2,0)*d2psi_curve(j,0) + d2psi_curve(j2,1)*d2psi_curve(j,1) + 2.0*d2psi_curve(j2,2)*d2psi_curve(j,2))*W*M(row2,j2)*M(row,j);
                      
                      }
                    } 
                  } 
                } 
                 
               //Loop over test functions associated with vertices
               for(unsigned l2=0;l2<3;l2++)
                {
                 for(unsigned k2=0;k2<6;k2++)
                  {
                   unsigned row2;
                   if(k2==0){row2 = l2;}
                   else if((k2==1) || (k2==2)){row2 = 2*(l2+1)+k2;}
                   else if(k2==3){row2 = 3*(l2+2)+3;}
                   else if(k2==4){row2 = 3*(l2+2)+5;}
                   else if(k2==5){row2 = 3*(l2+2)+4;}
                   //Get the local unknown
                   local_unknown = nodal_eqn(l2,k2);
                   //IF it's not a boundary condition
                   if(local_unknown >= 0)
                    {
                     for(unsigned j2=0;j2<36;j2++)
                      {
                       if(M(row2,j2)!=0)
                        {
                         //Add contribution to Elemental Matrix
                         jacobian(local_eqn,local_unknown) += (d2psi_curve(j2,0)*d2psi_curve(j,0) + d2psi_curve(j2,1)*d2psi_curve(j,1) + 2.0*d2psi_curve(j2,2)*d2psi_curve(j,2))*W*M(row2,j2)*M(row,j);
                         
                        }
                      }
                    } 
                  } 
                } 
              } // End of flag
            }
          }
        }
      }
     // 2) which associate with vertices
     for(unsigned l=0;l<3;l++)
      {
       for(unsigned k=0;k<6;k++)
        {
         unsigned row;
         if(k==0){row = l;}
         else if((k==1) || (k==2)){row = 2*(l+1)+k;}
         else if(k==3){row = 3*(l+2)+3;}
         else if(k==4){row = 3*(l+2)+5;}
         else if(k==5){row = 3*(l+2)+4;}
         //Get the local equation
         local_eqn = nodal_eqn(l,k);
         //IF it's not a boundary condition
         if(local_eqn >= 0)
          {
           for(unsigned j=0;j<36;j++)
            {
             if(M(row,j)!=0)
              {
               // Add body force/source term here 
               residuals[local_eqn] -= (source*psi_curve[j]*W)*M(row,j);
               
               residuals[local_eqn] += (interpolated_d2udxi(0,0)*d2psi_curve(j,0) + interpolated_d2udxi(0,1)*d2psi_curve(j,1) + 2.0*interpolated_d2udxi(0,2)*d2psi_curve(j,2))*W*M(row,j);
              
               // Calculate the jacobian
               //-----------------------
               if(flag)
                {
                 //Loop over test functions associated with internal nodes
                 for(unsigned p2=0;p2<3;p2++)
                  {
                   unsigned l2=0,row2=0;
                   if(p2==0){l2=3;row2=18;}
                   else if(p2==1){l2=4;row2=19;}
                   else if(p2==2){l2=5;row2=20;}
                   //Get the local equation
                   local_unknown = nodal_eqn(l2,0);
                   //IF it's not a boundary condition
                   if(local_unknown >= 0)
                    {
                     for(unsigned j2=0;j2<36;j2++)
                      {
                       if(M(row2,j2)!=0)
                        {
                         //Add contribution to Elemental Matrix
                         jacobian(local_eqn,local_unknown) += (d2psi_curve(j2,0)*d2psi_curve(j,0) + d2psi_curve(j2,1)*d2psi_curve(j,1) + 2.0*d2psi_curve(j2,2)*d2psi_curve(j,2))*W*M(row2,j2)*M(row,j);
                         
                        }
                      } 
                    } 
                  } 
                 //Loop over test functions associated with vertices
                 for(unsigned l2=0;l2<3;l2++)
                  {
                   for(unsigned k2=0;k2<6;k2++)
                    {
                     unsigned row2;
                     if(k2==0){row2 = l2;}
                     else if((k2==1) || (k2==2)){row2 = 2*(l2+1)+k2;}
                     else if(k2==3){row2 = 3*(l2+2)+3;}
                     else if(k2==4){row2 = 3*(l2+2)+5;}
                     else if(k2==5){row2 = 3*(l2+2)+4;}
                     //Get the local unknown
                     local_unknown = nodal_eqn(l2,k2);
                     //IF it's not a boundary condition
                     if(local_unknown >= 0)
                      {
                       for(unsigned j2=0;j2<36;j2++)
                        {
                         if(M(row2,j2)!=0)
                          {
                           //Add contribution to Elemental Matrix
                           jacobian(local_eqn,local_unknown) += (d2psi_curve(j2,0)*d2psi_curve(j,0) + d2psi_curve(j2,1)*d2psi_curve(j,1) + 2.0*d2psi_curve(j2,2)*d2psi_curve(j,2))*W*M(row2,j2)*M(row,j);
                           
                          }
                        } 
                      }
                    } 
                  } 
                } // End of flag
              }
            }
          }
        }
      }
    }
  } // End of loop over integration points 
}   



//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
unsigned  MyBiharmonicEquations<DIM,NNODE_1D>::self_test()
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
template <unsigned DIM, unsigned NNODE_1D>
void  MyBiharmonicEquations<DIM,NNODE_1D>::output(std::ostream &outfile, 
                                    const unsigned &nplot)
{

 //Vector of local coordinates
 Vector<double> s(DIM),x(DIM);
 
 // Tecplot header info
 //outfile << tecplot_zone_string(nplot);
 
 // Loop over plot points
 Vector<double> u(this->required_nvalue(0),0.0);
 unsigned num_plot_points=this->nplot_points(nplot);
 Vector<double> r(3);
 Vector<double> interpolated_xi(DIM,0.0);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   u = interpolated_u(s);
   // Get x position as Vector
   this->my_interpolated_x(s,x);
   
   for(unsigned i=0;i<DIM;i++) 
    {
     outfile << x[i] << " ";
    }

   // Loop for variables
   for(unsigned j=0;j<this->required_nvalue(0);j++)
    {
     outfile << u[j] << " " ;  
    }
    
   outfile << std::endl;
  }
}

 
//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void  MyBiharmonicEquations<DIM,NNODE_1D>::output(FILE* file_pt,
                                    const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Tecplot header info
 fprintf(file_pt,"%s",this->tecplot_zone_string(nplot).c_str());

 // Loop over plot points
 Vector<double> u(this->required_nvalue(0),0.0);
 unsigned num_plot_points=this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   
   for(unsigned i=0;i<DIM;i++) 
    {
     fprintf(file_pt,"%g ",this->interpolated_x(s,i));
    }
   u = interpolated_u(s);
   fprintf(file_pt,"%g \n",u[0]);
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 //write_tecplot_zone_footer(file_pt,nplot);
}



//======================================================================
 /// Output exact solution
 /// 
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points.
 ///
 ///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void MyBiharmonicEquations<DIM,NNODE_1D>::output_fct(std::ostream &outfile, 
                                       const unsigned &nplot, 
                  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
  
  // Vector for coordintes
  Vector<double> x(DIM);
  
 // Tecplot header info
 //outfile << tecplot_zone_string(nplot);
 
 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(this->required_nvalue(0),0.0);
 
 // Loop over plot points
 unsigned num_plot_points=this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   
   // Get x position as Vector
   this->my_interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,u_exact
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   // Loop over variables
   for(unsigned j=0;j<this->required_nvalue(0);j++)
   {
    outfile << exact_soln[j] << " ";
   }
   outfile <<  std::endl;  
  }
 
 // Write tecplot footer (e.g. FE connectivity lists)
 //write_tecplot_zone_footer(outfile,nplot);
}


//======================================================================
 /// Validate against exact solution
 /// 
 /// Solution is provided via function pointer.
 /// Plot error at a given number of plot points.
 ///
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void MyBiharmonicEquations<DIM,NNODE_1D>::compute_error(std::ostream &outfile, 
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
 unsigned n_intpt = this->integral_pt()->nweight();
 // Tecplot 
 //outfile << "ZONE" << std::endl;
 
 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(this->required_nvalue(0),0.0);
 
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
   double J = this->J_eulerian1(s);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   // Get x position as Vector
   this->my_interpolated_x(s,x);
    
   // Get FE function value
   Vector<double> u_fe(this->required_nvalue(0),0.0);
   u_fe = interpolated_u(s);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,error
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   
   for(unsigned ii=0;ii<this->required_nvalue(0);ii++)
    {
     outfile << exact_soln[ii] << " " << exact_soln[ii]-u_fe[ii] << " ";  
    }
   outfile << std::endl;
   
   // Loop over variables
   double tmp1=0.0, tmp2=0.0;
   for(unsigned ii=0;ii<1;ii++)
    {
     // Add to error and norm
     tmp1 = (exact_soln[ii]*exact_soln[ii]*W);
     tmp2 = ((exact_soln[ii]-u_fe[ii])*(exact_soln[ii]-u_fe[ii])*W);
     norm += tmp1;
     error += tmp2;
    }
  } //End of loop over integration pts
}





//====================================================================
// Force build of templates
//====================================================================
template class BiharmonicCurvedElement<2,2>;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



//==start_of_namespace================================================
/// Namespace for the solution of 2D Biharmonic equation
//====================================================================
namespace Physical_Variables
{
 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = cos(x[0])*exp(x[1]);
  u[1] = -1.0*exp(x[1])*sin(x[0]);
  u[2] = exp(x[1])*cos(x[0]);
  u[3] = -1.0*exp(x[1])*cos(x[0]);
  u[4] = exp(x[1])*cos(x[0]);
  u[5] = -1.0*exp(x[1])*sin(x[0]); 
 }

 /// Source function required to make the solution above an exact solution 
 void source_function(const Vector<double>& x, double& source)
 {
  source = 0.0;
 }
} // end of namespace



//==start_of_problem_class============================================
/// 2D Biharmonic problem.
//====================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D> 
class MyBiharmonicProblem : public Problem
{
public:

 /// Constructor: Pass number of elements and pointer to source function
 MyBiharmonicProblem(typename MyBiharmonicEquations<DIM,NNODE_1D>::SourceFctPt source_fct_pt,
                          typename MyBiharmonicEquations<DIM,NNODE_1D>::ExactSolnPt exact_pt,
                          const string& node_file_name,
                          const string& element_file_name,
                          const string& poly_file_name);

 /// Destructor (empty)
 ~MyBiharmonicProblem()
  {
   delete mesh_pt();
  }

 
 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve();

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// \short Doc the solution, pass the number of the case considered,
 /// so that output files can be distinguished.
 void doc_solution(DocInfo& doc_info);
 
 /// Perform a nodal permutation in curved elements
 void nodal_permutation();

private:

 /// Pointer to source function
 typename MyBiharmonicEquations<DIM,NNODE_1D>::SourceFctPt Source_fct_pt;

 /// Pointer to the exact solution
 typename MyBiharmonicEquations<DIM,NNODE_1D>::ExactSolnPt Exact_soln_pt;

}; // end of problem class





//=====start_of_constructor===============================================
/// \short Constructor for 2D Biharmonic problem.
/// Discretise the 2D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function. 
//========================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D>
MyBiharmonicProblem<ELEMENT,DIM,NNODE_1D>::MyBiharmonicProblem
(typename MyBiharmonicEquations<DIM,NNODE_1D>::SourceFctPt source_fct_pt,
 typename MyBiharmonicEquations<DIM,NNODE_1D>::ExactSolnPt exact_pt,
 const string& node_file_name,
 const string& element_file_name,
 const string& poly_file_name) : 
 Source_fct_pt(source_fct_pt), Exact_soln_pt(exact_pt)
{ 
 Problem::mesh_pt() = new TriangleMesh<ELEMENT>(node_file_name,element_file_name,poly_file_name);

 unsigned n_node = this->mesh_pt()->nnode();
 cout  << "nnode in the mesh = " << n_node << endl;
 // find number of elements in the mesh
 unsigned n_element = mesh_pt()->nelement();
 cout  << "nelement in the mesh = " << n_element << endl;
 
 // start_of_boundary_conditions
 // assign conditions at boundary 0
 unsigned i=0;
 unsigned  n_node1 = mesh_pt()->nboundary_node(i);
 cout  << "nnode on b0 = " << n_node1 << endl;

 for (unsigned n=0;n<n_node1;n++)
  {
   mesh_pt()->boundary_node_pt(i,n)->pin(0); 
   mesh_pt()->boundary_node_pt(i,n)->pin(1); 
   mesh_pt()->boundary_node_pt(i,n)->pin(2);
   mesh_pt()->boundary_node_pt(i,n)->pin(3);
   mesh_pt()->boundary_node_pt(i,n)->pin(4);
   mesh_pt()->boundary_node_pt(i,n)->pin(5);
  }
 // assign conditions at boundary 1
 i=1;
 n_node1 = mesh_pt()->nboundary_node(i);
 cout  << "nnode on b1 = " << n_node1 << endl;
 for (unsigned n=0;n<n_node1;n++)
  {
   mesh_pt()->boundary_node_pt(i,n)->pin(0); 
   mesh_pt()->boundary_node_pt(i,n)->pin(1); 
   mesh_pt()->boundary_node_pt(i,n)->pin(2);
   mesh_pt()->boundary_node_pt(i,n)->pin(3);
   mesh_pt()->boundary_node_pt(i,n)->pin(4);
   mesh_pt()->boundary_node_pt(i,n)->pin(5);
  }
 // assign conditions at boundary 2
 i=2;
 n_node1 = mesh_pt()->nboundary_node(2);
 cout  << "nnode on b2 = " << n_node1 << endl;
 for (unsigned n=0;n<n_node1;n++)
  {
   mesh_pt()->boundary_node_pt(i,n)->pin(0);
   mesh_pt()->boundary_node_pt(i,n)->pin(1); 
   mesh_pt()->boundary_node_pt(i,n)->pin(2);
   mesh_pt()->boundary_node_pt(i,n)->pin(3);
   mesh_pt()->boundary_node_pt(i,n)->pin(4);
   mesh_pt()->boundary_node_pt(i,n)->pin(5);
  }
 
 // end of boundary conditions
 
 // pinned for the extra node in straight-sided elements
 for(unsigned n=0;n<n_element;n++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(n));
   unsigned nnode = elem_pt->nnode();
   
   //Loop over nodes in an element to check whether it is a (curved)boundary element or not
   unsigned count = 0;
   // check just the vertice nodes
   for(unsigned i=0;i<3;i++)
    { 
     bool is_bd_node = elem_pt->node_pt(i)->is_on_boundary(1);
     if(is_bd_node==1)
      {
       count += 1;
      }
    }
   // pinned the unnecessary nodes in the straight-sided elements
   for(unsigned i=0;i<nnode;i++)
    {
     // take care internal nodes
     if(i>2)
      {
       if(count == 2) // it is a boundary element
        {
         //elem_pt->node_pt(i)->pin(0); // supposed to be free for the curved boundary
         elem_pt->node_pt(i)->pin(1);
         elem_pt->node_pt(i)->pin(2);
         elem_pt->node_pt(i)->pin(3);
         elem_pt->node_pt(i)->pin(4);
         elem_pt->node_pt(i)->pin(5);
        }
       else // it is NOT a boundary element: no internal nodes
        {
         elem_pt->node_pt(i)->pin(0);
         elem_pt->node_pt(i)->pin(1);
         elem_pt->node_pt(i)->pin(2);
         elem_pt->node_pt(i)->pin(3);
         elem_pt->node_pt(i)->pin(4);
         elem_pt->node_pt(i)->pin(5); 
        }
      }
    }
  }
  
 // Loop over elements and set pointers to Physical parameters
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   //Set the source function pointer and all physical variables
   elem_pt->source_fct_pt() = Source_fct_pt;
   
   //Set the exact solution pointer
   elem_pt->exact_pt() = Exact_soln_pt;
  } // end of pointers set up

  // Perform a nodal permutation
 nodal_permutation();

 // Setup equation numbering scheme
 assign_eqn_numbers();
} // end of constructor


//===start_of_nodal_permutation===========================================
/// \short Update the nodal permutation in curved elements: (Re)assign 
/// coordinates so that the first and the second nodes are associated
/// with the lower and upper position on a boundary, respectively. 
/// The third node is off-boundary node.
//========================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D>
 void MyBiharmonicProblem<ELEMENT,DIM,NNODE_1D>::nodal_permutation()
  {
   DenseMatrix<double> position(3,2);
   DenseMatrix<double> bd_position(20,2);
   Vector<double> x(2);
   unsigned bd_element;
  
   // Get number of elements in the mesh
   unsigned n_element = mesh_pt()->nelement();

   // Start a classification of nodes
   // Loop over elements to get coordinates of nodes on a boundary
   for(unsigned n=0;n<n_element;n++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(n));
     // we need just three vertices
     unsigned n_node = 3;
   
     for(unsigned l=0;l<n_node;l++)
      {
       for(unsigned j=0;j<DIM;j++)
        {
         position(l,j) = elem_pt->node_pt(l)->x(j);
        }
      }
    
     int count=0;
     
     // check whether this element is on boundary or not
     for(unsigned i=0;i<n_node;i++)
      {
       // check whether this node is on 1st boundary (curve boundary)
       bool bd_node1 = elem_pt->node_pt(i)->is_on_boundary(1);
       if(bd_node1==1)
        {
         count += 1;
         for(unsigned j=0;j<DIM;j++)
          {
           bd_position(count-1,j) = elem_pt->node_pt(i)->x(j);
          }
        }
       else
        {
         for(unsigned j=0;j<DIM;j++)
          {
           x[j] = elem_pt->node_pt(i)->x(j);
          }
        }
      }
     
     if(count == 2) // this case is a boundary element
      {
       bd_element = 1;
      }          
     else
      {
       bd_element = 0;           
      } // end_of_classification
     
     /// Permutation just for the curved elements
     if(bd_element==1)
      {
       unsigned n_position_type = 6;
       DenseMatrix<double> nodal_value(3,n_position_type,0.0);
       double angle_1 = atan(bd_position(0,1)/bd_position(0,0));
       double angle_2 = atan(bd_position(1,1)/bd_position(1,0));
       
       if(angle_1 > angle_2)
        {
         if(position(2,0)==x[0] && position(2,1)==x[1])                                  
          {
           Node* first_pt = elem_pt->node_pt(0);
           elem_pt->node_pt(0) = elem_pt->node_pt(1);
           elem_pt->node_pt(1) = first_pt;
           
           Node* second_pt = elem_pt->node_pt(3);
           elem_pt->node_pt(3) = elem_pt->node_pt(4);
           elem_pt->node_pt(5) = second_pt;                         
          }
         else if(position(1,0)==bd_position(0,0) && position(1,1)==bd_position(0,1))                                  
          {
           Node* first_pt = elem_pt->node_pt(2);
           elem_pt->node_pt(0) = elem_pt->node_pt(2);
           elem_pt->node_pt(2) = first_pt;
           
           Node* second_pt = elem_pt->node_pt(3);
           elem_pt->node_pt(3) = elem_pt->node_pt(5);
           elem_pt->node_pt(5) = second_pt;                         
          }
         else 
          { 
           Node* first_pt = elem_pt->node_pt(0);
           elem_pt->node_pt(0) = elem_pt->node_pt(2);
           elem_pt->node_pt(2) = elem_pt->node_pt(1);
           elem_pt->node_pt(1) = first_pt;
           
           Node* second_pt = elem_pt->node_pt(3);
           elem_pt->node_pt(3) = elem_pt->node_pt(5);
           elem_pt->node_pt(5) = elem_pt->node_pt(4);
           elem_pt->node_pt(4) = second_pt;
          }
        }
       else if(angle_1 < angle_2)
        {
         if(position(0,0)==bd_position(0,0) && position(0,1)==bd_position(0,1))
          {
           if(position(2,0)==bd_position(1,0) && position(2,1)==bd_position(1,1))                                  
            {
             Node* first_pt = elem_pt->node_pt(1);
             elem_pt->node_pt(1) = elem_pt->node_pt(2);
             elem_pt->node_pt(2) = first_pt;
             
             Node* second_pt = elem_pt->node_pt(4);
             elem_pt->node_pt(4) = elem_pt->node_pt(5);
             elem_pt->node_pt(5) = second_pt;                         
            }
          }
         else 
          {
           
           Node* first_pt = elem_pt->node_pt(0);
           elem_pt->node_pt(0) = elem_pt->node_pt(1);
           elem_pt->node_pt(1) = elem_pt->node_pt(2);
           elem_pt->node_pt(2) = first_pt;
           
           Node* second_pt = elem_pt->node_pt(3);
           elem_pt->node_pt(3) = elem_pt->node_pt(4);
           elem_pt->node_pt(4) = elem_pt->node_pt(5);
           elem_pt->node_pt(5) = second_pt;
          }   
        }
      }
    }
  }  // end of permutation
 


//===start_of_actions_before_newton_solve========================================
/// \short Update the problem specs before solve: (Re)set boundary values
/// from the exact solution. 
//========================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D>
void MyBiharmonicProblem<ELEMENT,DIM,NNODE_1D>::actions_before_newton_solve()
{
 // How many boundaries are there?
 unsigned n_bound = mesh_pt()->nboundary();
 
 //Loop over the boundaries
 for(unsigned i=0;i<n_bound;i++)
  {
   // How many nodes are there on this boundary?
   unsigned n_node = mesh_pt()->nboundary_node(i);

   // Loop over the nodes on boundary
   for (unsigned n=0;n<n_node;n++)
    {
     // Get pointer to node
     Node* nod_pt=mesh_pt()->boundary_node_pt(i,n);

     // Extract nodal coordinates from node:
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     // Compute the value of the exact solution at the nodal point
     Vector<double> u(6);
     for(unsigned i=0;i<6;i++)
      {
       Physical_Variables::get_exact_u(x,u);
       
       // Assign the value to the one (and only) nodal value at this node
       nod_pt->set_value(i,u[i]);
      }
    }
  }
} // end of actions before solve


//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D>
void MyBiharmonicProblem<ELEMENT,DIM,NNODE_1D>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=10; 

 // Output solution with specified number of plot points per element
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 some_file.precision(20);
 mesh_pt()->output(some_file,npts);
 some_file.close();

 // Output exact solution at much higher resolution (so we can
 // see how well the solutions agree between nodal points)
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,Physical_Variables::get_exact_u); 
 some_file.close();

 // Doc pointwise error and compute norm of error and of the solution
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,Physical_Variables::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc error norm:
 cout << "\nNorm of error    : " << sqrt(error) << std::endl; 
 cout << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
 cout << std::endl;

} // end of doc


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//======start_of_main==================================================
/// Driver for 2D Biharmonic problem
//=====================================================================
int main(int argc, char* argv[])
{
 //---------------------------------
 // Check number of command line arguments: Need exactly two.
 if (argc!=4)
  {
   std::string error_message =
    "Wrong number of command line arguments.\n";
   error_message +=
    "Must specify the following file names  \n";
   error_message += 
    "filename.node then filename.ele then filename.poly\n";

   throw OomphLibError(error_message,
                       "main()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // Convert arguments to strings that specify th input file names
 string node_file_name(argv[1]);
 string element_file_name(argv[2]);
 string poly_file_name(argv[3]);

 // Set up the problem: 
 MyBiharmonicProblem<BiharmonicCurvedElement<2,2>,2,2> //Element type as template parameter
  problem(Physical_Variables::source_function,Physical_Variables::get_exact_u,
          node_file_name,element_file_name,poly_file_name);

 // Check whether the problem can be solved
 cout << "\n\n\nProblem self-test ";
 if (problem.self_test()==0)  
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("failed!","main()",OOMPH_EXCEPTION_LOCATION);
  }
 
 // Set up doc info
 DocInfo doc_info;
 doc_info.set_directory("RESLT_curved");
 doc_info.number()=0;

 // Solve the problem
 problem.newton_solve();

 //Output the solution
 problem.doc_solution(doc_info);

} // end of main




