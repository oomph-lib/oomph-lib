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


//Header file for 2d linear shell problem: circular plate bending
#ifndef OOMPH_LINEAR_SHELL_ELEMENTS_HEADER
#define OOMPH_LINEAR_SHELL_ELEMENTS_HEADER

namespace oomph
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Geometric object
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//=========================================================================
/// \short Circle in 3D space.
/// \f[ x = X_c   \f]
/// \f[ y = Y_c   \f]
/// \f[ z = 0.0
//=========================================================================
class CircularPlate : public GeomObject
{

public:

 /// Constructor: Specify radius
 CircularPlate(const double& a, const double& b) : 
  GeomObject(2,3), A(a), B(b){}
 
 /// Broken copy constructor
 CircularPlate(const CircularPlate& node) 
  { 
   BrokenCopy::broken_copy("CircularPlate");
  } 
 
 /// Broken assignment operator
 void operator=(const CircularPlate&) 
  {
   BrokenCopy::broken_assign("CircularPlate");
  }

 /// Access function to x-center
 double& a() {return A;}

 /// Access function to y-center
 double& b() {return B;}

 /// Position vector
 void position(const Vector<double>& zeta, Vector<double>& r)const
  {
   r[0]=A + zeta[0];
   r[1]=B + zeta[1];
   r[2]=0.0;
  }

 /// \short How many items of Data does the shape of the object depend on?
 virtual unsigned ngeom_data() const
  {
   return 0;
  }

 /// \short Position Vector and 1st and 2nd derivs w.r.t. zeta.
 void d2position(const Vector<double> &zeta,
                 RankThreeTensor<double> &ddrdzeta) const
  {
   ddrdzeta(0,0,0) = 0.0;
   ddrdzeta(0,0,1) = 0.0;
   ddrdzeta(0,0,2) = 0.0;
   
   ddrdzeta(1,1,0) = 0.0;
   ddrdzeta(1,1,1) = 0.0;
   ddrdzeta(1,1,2) = 0.0;
   
   //Mixed derivatives
   ddrdzeta(0,1,0) = ddrdzeta(1,0,0) = 0.0;
   ddrdzeta(0,1,1) = ddrdzeta(1,0,1) = 0.0;
   ddrdzeta(0,1,2) = ddrdzeta(1,0,2) = 0.0;
  }

 /// \short Position Vector and 1st and 2nd derivs w.r.t. zeta.
 void d2position(const Vector<double>& zeta, Vector<double>& r,
                 DenseMatrix<double> &drdzeta,
                 RankThreeTensor<double> &ddrdzeta) const
  {
   //Let's just do a simple tube
   r[0]=A + zeta[0];
   r[1]=B + zeta[1];
   r[2]=0.0;

   //Do the derivatives drdzeta0
   drdzeta(0,0) = 1.0;
   drdzeta(0,1) = 0.0;
   drdzeta(0,2) = 0.0;
   
   //Do the azimuthal derivatives drdzeta1
   drdzeta(1,0) = 0.0;
   drdzeta(1,1) = 1.0;
   drdzeta(1,2) = 0.0;
   
   //Now let's do the second derivatives
   ddrdzeta(0,0,0) = 0.0;
   ddrdzeta(0,0,1) = 0.0;
   ddrdzeta(0,0,2) = 0.0;
   
   ddrdzeta(1,1,0) = 0.0;
   ddrdzeta(1,1,1) = 0.0;
   ddrdzeta(1,1,2) = 0.0;
   
   //Mixed derivatives
   ddrdzeta(0,1,0) = ddrdzeta(1,0,0) = 0.0;
   ddrdzeta(0,1,1) = ddrdzeta(1,0,1) = 0.0;
   ddrdzeta(0,1,2) = ddrdzeta(1,0,2) = 0.0;
  }

private:

 /// x-center
 double A;

 /// y-center
 double B;

};


//=============================================================
/// A class for all subparametric elements that solve the 
/// linear shell equations.
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
template <unsigned DIM, unsigned NNODE_1D>
class MyShellEquations : public virtual C1CurvedElement<DIM,NNODE_1D>
{

public:
 
 /// \short Function pointer to source function fct(x,f(x)) -- 
 /// x is a Vector! 
 typedef void (*SourceFctPt)(const Vector<double>& x, const Vector<double>& unit_n, Vector<double>& f);


 /// \short Function pointer to gradient of source function  fct(x,g(x)) -- 
 /// x is a Vector! 
 typedef void (*SourceFctGradientPt)(const Vector<double>& x, 
                                            Vector<double>& gradient);

 /// \short Function pointer to the exact solutions -- 
 /// x is a Vector! 
 typedef void (*ExactSolnPt)(const Vector<double>& x, Vector<double>& exact_soln);

 /// Constructor (must initialise the Source_fct_pt to null)
 MyShellEquations() : Source_fct_pt(0), Source_fct_gradient_pt(0)
  {}
 
 /// Broken copy constructor
 MyShellEquations(const MyShellEquations& dummy) 
  { 
   BrokenCopy::broken_copy("MyShellEquations");
  } 
 
 /// Broken assignment operator
 void operator=(const MyShellEquations&) 
  {
   BrokenCopy::broken_assign("MyShellEquations");
  }

 /// \short Return the index at which the unknown value
 /// is stored. 
 /// In derived multi-physics elements, this function should be overloaded
 /// to reflect the chosen storage scheme. Note that these equations require
 /// that the unknown is always stored at the same index at each node.
  virtual inline unsigned u_index_shell() const {return this->required_nvalue(0);}

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
    "There is no time-dependent output_fct() for the elements ",
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
    "There is no time-dependent compute_error() for the elements",
    "MyShellEquations<DIM>::compute_error()",
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

 /// Access function: Undeformed shell 
 GeomObject*& undeformed_midplane_pt() {return Undeformed_midplane_pt;}
 

 /// Get exact solutions at (Eulerian) position x.
 inline virtual void get_exact_solution(const Vector<double>& x,
                                        Vector<double>& exact_soln) const
  {
   (*Exact_soln_pt)(x,exact_soln);
  }

 /// Get source term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the source function might be determined by
 /// another system of equations.
 inline virtual void get_source_shell(const unsigned& ipt,
                                        const Vector<double>& x,
                                        const Vector<double>& unit_n,
                                        Vector<double>& source) const
  {
   //If no source function has been set, return zero
   if(Source_fct_pt==0) 
    {
     for(unsigned i=0;i<(u_index_shell())/2;i++)
      {
       source[i] = 0.0;
      }
    }
   else
    {
     // Get source strength
     (*Source_fct_pt)(x,unit_n,source);
    }
  }

 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_shell(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }

 
 /// Add the element's contribution to its residual vector and 
 /// element Jacobian matrix (wrapper)
 /*void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                     DenseMatrix<double> &jacobian)
  {
   Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_shell(residuals,jacobian,1);
   }*/
 


 /// \short Return FE representation of function value u(s) 
 /// at local coordinate s
 inline Vector<double> interpolated_u_shell(const Vector<double> &s) const
  {
   //Get the index at which the unknown is stored
   const unsigned u_nodal_index = u_index_shell();
   //Initialise value of u
   Vector<double> interpolated_u(u_nodal_index,0.0);

   //Interpolated unknowns
   Vector<double> u(1,0.0);
   DenseMatrix<double> interpolated_dudxi(1,2,0.0), interpolated_d2udxi(1,3,0.0);
  
   //Interpolated unknown in the normal direction
   this->my_interpolated_u_normal(s,u,interpolated_dudxi,interpolated_d2udxi);
   interpolated_u[2] = u[0];
   interpolated_u[3] = interpolated_dudxi(0,0);
   interpolated_u[4] = interpolated_dudxi(0,1);
   interpolated_u[5] = interpolated_d2udxi(0,0);
   interpolated_u[6] = interpolated_d2udxi(0,1);
   interpolated_u[7] = interpolated_d2udxi(0,2);

   //Interpolated unknowns in the tangential directions
   Vector<double> ut(2,0.0);
   DenseMatrix<double> interpolated_dutdxi(2,2,0.0);
   this->my_interpolated_u_tangential(s,ut,interpolated_dutdxi);
   interpolated_u[0] = ut[0];
   interpolated_u[1] = ut[1];
 
   return(interpolated_u);
  }

 /// \short Self-test: Return 0 for OK
 unsigned self_test();


protected:

 /// \short Shape/test functions and derivs w.r.t. to global coords at 
 /// local coord. s; return  Jacobian of mapping
 virtual double d2shape_and_d2test_eulerian_shell(const Vector<double> &s, 
                                                  Shape &psi, 
                                                  DShape &dpsidx, DShape &d2psidx, Shape &test, 
                                                  DShape &dtestdx, DShape &d2testdx) const=0;
 virtual double dshape_and_dtest_eulerian_shell(const Vector<double> &s,
                                                  Shape &psi,
                                                  DShape &dpsidx, Shape &test,
                                                  DShape &dtestdx) const=0;

 /// \short Shape/test functions and derivs w.r.t. to global coords at 
 /// integration point ipt; return  Jacobian of mapping
 virtual double d2shape_and_d2test_eulerian_at_knot_shell(const unsigned &ipt, 
                                                          Shape &psi, 
                                                          DShape &dpsidx,
                                                          DShape &d2psidx,
                                                          Shape &test, 
                                                          DShape &dtestdx,
                                                          DShape &d2testdx) 
  const=0;

 virtual double dshape_and_dtest_eulerian_at_knot_shell(const unsigned &ipt,
                                                          Shape &psi,
                                                          DShape &dpsidx,
                                                          Shape &test,
                                                          DShape &dtestdx)
  const=0;

 /// \short Compute element residual Vector only (if flag=and/or element 
 /// Jacobian matrix 
 virtual void fill_in_generic_residual_contribution_shell(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  const unsigned& flag); 
 
 /// Pointer to source function:
 SourceFctPt Source_fct_pt;
 
 /// Pointer to the exact solution:
 ExactSolnPt Exact_soln_pt;

 /// Pointer to gradient of source function
 SourceFctGradientPt Source_fct_gradient_pt;

 /// Pointer to undeformed shell:
 GeomObject* Undeformed_midplane_pt;
 
 /// Pointer to axial prestress
 double* Sigma0_pt;
 
 /// Pointer to wall thickness
 double* H_pt;
 
 /// Pointer to Timescale ratio
 double* Lambda_sq_pt;
};

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//======================================================================
/// C1CurvedShellElement elements are with subparametric interpolation for the function.
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
class C1CurvedShellElement : public virtual MyShellEquations<DIM,NNODE_1D>
{

private:

 /// \short Static int that holds the number of variables at 
 /// nodes: always the same
 static const unsigned Initial_Nvalue;
 
  public:


 ///\short  Constructor: Call constructors for C1CurvedShellElement and 
 /// MyShell equations
 C1CurvedShellElement() : MyShellEquations<DIM,NNODE_1D>()
  {}
 
 /// Broken copy constructor
 C1CurvedShellElement(const C1CurvedShellElement<DIM,NNODE_1D>& dummy) 
  { 
   BrokenCopy::broken_copy("C1CurvedShellElement");
  } 
 
 /// Broken assignment operator
 void operator=(const C1CurvedShellElement<DIM,NNODE_1D>&) 
  {
   BrokenCopy::broken_assign("C1CurvedShellElement");
  }


 /// \short  Required  # of `values' (pinned or dofs) 
 /// at node n
 inline unsigned required_nvalue(const unsigned &n) const 
  {return Initial_Nvalue;}

 /// \short Output function:  
 ///  x,y,u   or    x,y,z,u
 void output(std::ostream &outfile)
  {MyShellEquations<DIM,NNODE_1D>::output(outfile);}


 ///  \short Output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot)
  {MyShellEquations<DIM,NNODE_1D>::output(outfile,n_plot);}


 /// \short C-style output function:  
 ///  x,y,u   or    x,y,z,u
 void output(FILE* file_pt)
  {MyShellEquations<DIM,NNODE_1D>::output(file_pt);}


 ///  \short C-style output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot)
  {MyShellEquations<DIM,NNODE_1D>::output(file_pt,n_plot);}


 /// \short Output function for an exact solution:
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {MyShellEquations<DIM,NNODE_1D>::output_fct(outfile,n_plot,exact_soln_pt);}



 /// \short Output function for a time-dependent exact solution.
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 /// (Calls the steady version)
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {MyShellEquations<DIM,NNODE_1D>::output_fct(outfile,n_plot,time,exact_soln_pt);}


protected:

/// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double d2shape_and_d2test_eulerian_shell(
  const Vector<double> &s, Shape &psi, DShape &dpsidx, DShape &d2psidx, 
  Shape &test, DShape &dtestdx, DShape &d2testdx) const;

 inline double dshape_and_dtest_eulerian_shell(
  const Vector<double> &s, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const;


 /// \short Shape, test functions & derivs. w.r.t. to global coords. at
 /// integration point ipt. Return Jacobian.
 inline double d2shape_and_d2test_eulerian_at_knot_shell(const unsigned& ipt,
                                                         Shape &psi, 
                                                         DShape &dpsidx, 
                                                         DShape &d2psidx,
                                                         Shape &test,
                                                         DShape &dtestdx,
                                                         DShape &d2testdx) 
  const;

 inline double dshape_and_dtest_eulerian_at_knot_shell(const unsigned &ipt,
                                                         Shape &psi,
                                                         DShape &dpsidx,
                                                         Shape &test,
                                                         DShape &dtestdx)
  const;

};




//Inline functions:

//=======================================================================
/// Face geometry for the C1CurvedShellElement elements: The spatial 
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<C1CurvedShellElement<DIM,NNODE_1D> >: 
 public virtual TElement<DIM-1,NNODE_1D>
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional TElement
 FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}

};


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
 double C1CurvedShellElement<DIM,NNODE_1D>::dshape_and_dtest_eulerian_shell(
  const Vector<double> &s,
  Shape &psi,
  DShape &dpsidx,
  Shape &test,
  DShape &dtestdx) const
{
 const double J = this->dshape_eulerian(s,psi,dpsidx);
 test = psi;
 dtestdx = dpsidx;
 return J;
}

template<unsigned DIM, unsigned NNODE_1D>
 double C1CurvedShellElement<DIM,NNODE_1D>::d2shape_and_d2test_eulerian_shell(
  const Vector<double> &s,
  Shape &psi, 
  DShape &dpsidx,
  DShape &d2psidx,
  Shape &test, 
  DShape &dtestdx,
  DShape &d2testdx) const
{
 //Call the geometrical shape functions and derivatives  
 const double J = this->d2shape_eulerian(s,psi,dpsidx,d2psidx);

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
double C1CurvedShellElement<DIM,NNODE_1D>::
 dshape_and_dtest_eulerian_at_knot_shell(
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
double C1CurvedShellElement<DIM,NNODE_1D>::
 d2shape_and_d2test_eulerian_at_knot_shell(
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
 const unsigned C1CurvedShellElement<DIM,NNODE_1D>::Initial_Nvalue = 8;

//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void  MyShellEquations<DIM,NNODE_1D>::
fill_in_generic_residual_contribution_shell(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              const unsigned& flag) 
{
 //Find out how many nodes there are
 const unsigned n_node = this->nnode()-3;
 
 //Find out how many nodes positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();
 
 //Set the dimension of the global coordinates
 unsigned n_dim = Undeformed_midplane_pt->ndim();
  
 //Set the number of lagrangian coordinates
 unsigned n_lagrangian =  Undeformed_midplane_pt->nlagrangian();
  
 //Set up memory for the shape and test functions
 Shape psi1(3,n_position_type), test(n_node), psi(n_node);
 DShape dpsidxi1(3,n_position_type,DIM), dtestdxi(n_node,DIM), dpsidxi(n_node,DIM);
 DShape d2psidxi1(3,n_position_type,3), d2testdxi(n_node,3), d2psidxi(n_node,3);
 
 
  //Set the value of n_intpt
 const unsigned n_intpt = this->integral_pt()->nweight();

 //Integers to store the local equation and unknown numbers
 int local_eqn=0;
      
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);
  
   //Call the derivatives of the c0-shape and test functions for tangential direcion--> Lagrange shape function
   double J = dshape_and_dtest_eulerian_at_knot_shell(ipt,psi,dpsidxi,test,dtestdxi);
   
   //Call the derivatives of the c1-shape and test functions for normal direction--> Bell shape function for straight triangular element
   double J1 = d2shape_and_d2test_eulerian_at_knot_shell(ipt,psi1,dpsidxi1,d2psidxi1,test,dtestdxi,d2testdxi);

   //Premultiply the weights and the Jacobian
   double W1 = w*J1;
   double W = w*J;
   
   double H = 0.01;
   
   //Calculate local values of unknowns
   Vector<double> u_value(n_position_type,0.0);
   Vector<double> interpolated_u(3,0.0);
   DenseMatrix<double> interpolated_dudxi(3,DIM,0.0);
   DenseMatrix<double> interpolated_d2udxi(3,3,0.0);

   //Calculate displacement value and derivatives:
   //-----------------------------------------
           
   // in normal direction
   Vector<double> u(1,0.0);
   DenseMatrix<double> dudxi(1,2,0.0);
   DenseMatrix<double> d2udxi(1,3,0.0);

   Vector<double> integration_point(2);
   integration_point[0] = this->integral_pt()->knot(ipt,0);
   integration_point[1] = this->integral_pt()->knot(ipt,1);

   this->my_interpolated_u_normal(integration_point,u,dudxi,d2udxi);
  
   interpolated_u[2] = u[0];
       
   // Loop over directions
   for(unsigned j=0;j<n_lagrangian;j++)
    {
     interpolated_dudxi(2,j) = dudxi(0,j);
    }
   
   // Loop over the second derivative directions
   for(unsigned j=0;j<3;j++)
    {
     interpolated_d2udxi(2,j) = d2udxi(0,j);        
    } 

   //in both tangential directions
   //Interpolated in tangential direction
   Vector<double> ut(2,0.0);
   DenseMatrix<double> interpolated_dutdxi(2,2,0.0);
   this->my_interpolated_u_tangential(integration_point,ut,interpolated_dutdxi);
    
   for(unsigned i=0;i<2;i++) 
    {
     interpolated_u[i] = ut[i];
     // Loop over the first derivative directions
     for(unsigned j=0;j<n_lagrangian;j++)
      {
       interpolated_dudxi(i,j) = interpolated_dutdxi(i,j);
      }        
    }
    
   //-------------------------------------------------------
   // setup position vector and derivatives of undeformed configuration
   Vector<double> r(n_dim);
   DenseMatrix<double> a(n_lagrangian,n_dim), a_tn(n_lagrangian,n_dim);
   RankThreeTensor<double> dadxi(n_lagrangian,n_lagrangian,n_dim);
   
   Vector<double> interpolated_xi(DIM,0.0);   
   this->my_interpolated_x(integration_point,interpolated_xi);
   
   // get the undeformed geometry
   Undeformed_midplane_pt->d2position(interpolated_xi,r,a,dadxi);
   
   // calculate the covariant metric tensors
   double tensor_a[2][2], aup[2][2];

   for(unsigned i=0;i<2;i++)
    {
     for(unsigned j=0;j<2;j++)
      {
       /// initialise tensors to zero
       tensor_a[i][j]=0.0;
       for(unsigned k=0;k<n_dim;k++)
        {
         tensor_a[i][j] += a(i,k)*a(j,k); 
        }
      }
    }
    
   // Calculate determinant
   double adet = tensor_a[0][0]*tensor_a[1][1] - tensor_a[0][1]*tensor_a[1][0];
   
   // calculate entries of the inverse
   aup[0][0] = tensor_a[1][1]/adet;
   aup[0][1] = -1*tensor_a[1][0]/adet;
   aup[1][0] = -1*tensor_a[0][1]/adet;
   aup[1][1] = tensor_a[0][0]/adet;

   // calculate the unit normal vectors
   Vector<double> unit_n(n_dim),n(n_dim);
   double magnitude_n;
   n[0] = (a(0,1)*a(1,2) - a(0,2)*a(1,1));
   n[1] = (a(0,2)*a(1,0) - a(0,0)*a(1,2));
   n[2] = (a(0,0)*a(1,1) - a(0,1)*a(1,0)); 
   magnitude_n = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
   unit_n[0] = n[0]/magnitude_n;
   unit_n[1] = n[1]/magnitude_n;
   unit_n[2] = n[2]/magnitude_n;
   
   // allocate all tangential and normal based vectors
   DenseMatrix<double> t(3,n_dim,0.0);
   //Calculate for tangent based vector
   for(unsigned i=0;i<3;i++)
    {
     for(unsigned j=0;j<n_dim;j++)
      {
       if(i==2)
        {
         t(2,j) = unit_n[j];
        }
       else
        {
         t(i,j) = a(i,j);     
        }
      }        
    }
   
   // compute the derivatives of normal vector
   DenseMatrix<double> dndxi(n_lagrangian,n_dim);
   Vector<double> vec(3);
   // loop over lagrangian coordinates
   for(unsigned i=0;i<n_lagrangian;i++)
    {
     vec[0] = (dadxi(0,i,1)*t(1,2) -  dadxi(0,i,2)*t(1,1)) - (dadxi(1,i,1)*t(0,2) -  dadxi(1,i,2)*t(0,1)); 
     vec[1] = (dadxi(0,i,2)*t(1,0) -  dadxi(0,i,0)*t(1,2)) - (dadxi(1,i,2)*t(0,0) -  dadxi(1,i,0)*t(0,2)); 
     vec[2] = (dadxi(0,i,0)*t(1,1) -  dadxi(0,i,1)*t(1,0)) - (dadxi(1,i,0)*t(0,1) -  dadxi(1,i,1)*t(0,0)); 
     
     dndxi(i,0) = vec[0]/magnitude_n ;
     dndxi(i,1) = vec[1]/magnitude_n ;
     dndxi(i,2) = vec[2]/magnitude_n ;
    }
   
   
   // compute for the derivative of the based vectors in the tangential and normal directions -> Christoffel symbols
   // \Gamma^{m}_{i,alpha} = dtdxi(i,alpha,m) = (dadxi(i,alpha,m)t(m,k))
   RankThreeTensor<double> dtdxi(n_dim,n_lagrangian,n_dim);
   
   // Initialize matrix
    for(unsigned i=0;i<n_dim;i++)
    {
     for(unsigned j=0;j<n_lagrangian;j++)
      {
       for(unsigned m=0;m<n_dim;m++)
        {
         dtdxi(i,j,m) = 0.0; 
        }            
      }            
    }
   
   for(unsigned i=0;i<n_lagrangian;i++)
    {
     for(unsigned j=0;j<n_lagrangian;j++)
      {
       for(unsigned m=0;m<n_dim;m++)
        {
         for(unsigned k=0;k<n_dim;k++)
          {
           dtdxi(i,j,m) += dadxi(i,j,k)*t(m,k); 
          }
        }            
      }            
    }
   
   for(unsigned i=0;i<n_lagrangian;i++)
    {
     for(unsigned m=0;m<n_dim;m++)
      {
       for(unsigned k=0;k<n_dim;k++)
        {
         dtdxi(2,i,m) += dndxi(i,k)*t(m,k); 
        }
      }
    }
   
   // compute the first derivative of a displacement 
   // in tangential and normal coordinate system
   DenseMatrix<double> interpolated_dudxi_tn(n_dim,n_lagrangian,0.0);
  
   for(unsigned i=0;i<n_lagrangian;i++)
    {
     for(unsigned k=0;k<n_dim;k++)
      {
       interpolated_dudxi_tn(k,i) = interpolated_dudxi(k,i);
       for(unsigned j=0;j<n_dim;j++)
        {
         interpolated_dudxi_tn(k,i) += interpolated_u[j]*dtdxi(j,i,k);
        }
      }
    }
   
   // Compute the linearised strain tensor
   double gamma[2][2];
   
   // Initialise the matrix
   for(unsigned al=0;al<2;al++)
    {
     for(unsigned be=0;be<2;be++)
      {
       gamma[al][be] = 0.0;
      }
    }
   
   for(unsigned al=0;al<2;al++)
    {
     for(unsigned be=0;be<2;be++)
      {
       gamma[al][be] = ( interpolated_dudxi_tn(al,be) + interpolated_dudxi_tn(be,al))/2.0; 
      
      }
    }
  
   // compute the second derivative of a displacement 
   // in tangential and normal coordinate system
   DenseMatrix<double> interpolated_d2udxi_tn(n_dim,3,0.0);
   
   unsigned al,be;
   // Loop over derivative
   for(unsigned c=0;c<3;c++)
   {
    if(c==0)
     { 
      al=0;
      be=0;
     }
    else if(c==1)
     {
      al=1;
      be=1;    
     }
    else if(c==2)
     {
      al=0;
      be=1;    
     }

    // Loop over components
    for(unsigned k=2;k<=2;k++)
     {
      interpolated_d2udxi_tn(k,c) += interpolated_d2udxi(k,c);

      for(unsigned j=0;j<n_dim;j++)
       {
        interpolated_d2udxi_tn(k,c) += interpolated_dudxi_tn(j,al)*dtdxi(j,be,k);
       }
      for(unsigned i=0;i<n_dim;i++)
       {
        interpolated_d2udxi_tn(k,c) += interpolated_dudxi(i,be)*dtdxi(i,al,k);  
       }  
     }
   }

   // Pre-compute for cross product terms
   Vector<double> L(n_dim,0.0),temp(n_dim,0.0);
   
   L[0] = tensor_a[0][1]*interpolated_dudxi_tn(2,1) - tensor_a[1][1]*interpolated_dudxi_tn(2,0);
   L[1] = -1.0*tensor_a[0][0]*interpolated_dudxi_tn(2,1) + tensor_a[1][0]*interpolated_dudxi_tn(2,0);
   L[2] = tensor_a[0][0]*interpolated_dudxi_tn(1,1) - tensor_a[0][1]*interpolated_dudxi_tn(0,1) 
          + tensor_a[1][1]*interpolated_dudxi_tn(0,0) - tensor_a[1][0]*interpolated_dudxi_tn(1,0);
   
   // calculate for the linearised bending tensor kappa 
   double kappa[2][2];
   
   // Initialise the bending matrix
   for(unsigned al=0;al<2;al++)
    {
     for(unsigned be=0;be<2;be++)
      {
       kappa[al][be] = 0.0;
      }
    }
   
   unsigned c;
   for(unsigned al=0;al<2;al++)
    {
     for(unsigned be=0;be<2;be++)
      {
       if(al==0 && be==0) {c=0;}
       else if(al==1 && be==1) {c=1;}
       else if(al!=be) {c=2;}
       
       for(unsigned k=0;k<n_dim;k++)
        {
         kappa[al][be] -= 1.0/adet*(dtdxi(al,be,k)*L[k]);
        }
       kappa[al][be] += 1.0/adet/adet*L[2]*dtdxi(al,be,2);
       kappa[al][be] -= interpolated_d2udxi_tn(2,c);
      }
    }
   
   // calculate the plane stress stiffness tensor
   double Et[2][2][2][2];
   double Nu = 0.3;
   // some constants
   double C1 = 0.5*(1.0-Nu);
   double C2 = Nu;

   //Loop over first index
   for(unsigned i=0;i<2;i++)
    {
     for(unsigned j=0;j<2;j++)
      {
       for(unsigned k=0;k<2;k++)
        {
         for(unsigned l=0;l<2;l++)
          {
           Et[i][j][k][l] = C1*(aup[i][k]*aup[j][l] + aup[i][l]*aup[j][k]) + C2*aup[i][j]*aup[k][l];
          }
        }
      }
    }

   //Get source function
   //-------------------
   Vector<double> source(n_dim),f(n_dim);
   get_source_shell(ipt,interpolated_xi,unit_n,source);

   // decompose in the tangential and normal coordinates
   f[0] = source[0]*t(0,0) + source[1]*t(0,1) + source[2]*t(0,2);
   f[1] = source[0]*t(1,0) + source[1]*t(1,1) + source[2]*t(1,2);
   f[2] = source[0]*t(2,0) + source[1]*t(2,1) + source[2]*t(2,2);
   source = f;
 
   //--------------------------------
   // Assemble residuals and Jacobian
   //--------------------------------
   //Check whether this element is a boundary element or not
   unsigned bd_element = this->is_boundary_element();

   //-----------------------------------------------------------------
   //This is the case that it is NOT a boundary element. We can use 
   //the existing bell shape functions as a test function without any 
   //modification.
   //-----------------------------------------------------------------
   // Loop over directions
   for(unsigned i=0;i<n_dim;i++)
    {
     // compute for variations
     // This is the case for normal direction which we have to concern
     // whether the element is curve or straight side
     if(i==2)
      {
       // striaght sided element
       if(bd_element==0) 
        {
         // Loop over the test functions
         for(unsigned l=0;l<3;l++)
          {
           for(unsigned p=0;p<n_position_type;p++)
            {
             // initialise 
             DenseMatrix<double> delta_kappa(2,2,0.0);
             DenseMatrix<double> delta_gamma(2,2,0.0);
             
             for(unsigned al=0;al<2;al++)
              {
               for(unsigned be=0;be<2;be++)
                {
                 // compute for variations of the strain tensor
                 delta_gamma(al,be) = (dtdxi(2,be,al) + dtdxi(2,al,be))/2.0*psi1(l,p);
                 
                 // compute for variations of the bending tensor
                 if((al==0) & (be==0)) {c=0;}
                 else if((al==1) & (be==1)) {c=1;}
                 else if(al!=be) {c=2;}
                 
                 delta_kappa(al,be) -= d2psidxi1(l,p,c);  
                
                 delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,0)*( dtdxi(2,1,2)*tensor_a[0][1]*psi1(l,p)-dtdxi(2,0,2)*tensor_a[1][1]*psi1(l,p) );
                 
                 delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,1)*( -1*dtdxi(2,1,2)*tensor_a[0][0]*psi1(l,p) + dtdxi(2,0,2)*tensor_a[1][0]*psi1(l,p) );
                 
                 delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,2)*( dtdxi(2,1,1)*tensor_a[0][0]*psi1(l,p)-dtdxi(2,1,0)*tensor_a[0][1]*psi1(l,p) 
                                                                  + dtdxi(2,0,0)*tensor_a[1][1]*psi1(l,p)-dtdxi(2,0,1)*tensor_a[1][0]*psi1(l,p) );
                 
                 delta_kappa(al,be) += dtdxi(al,be,2)/adet/adet*( dtdxi(2,1,1)*tensor_a[0][0]*psi1(l,p)-dtdxi(2,1,0)*tensor_a[0][1]*psi1(l,p) 
                                                                  + dtdxi(2,0,0)*tensor_a[1][1]*psi1(l,p)-dtdxi(2,0,1)*tensor_a[1][0]*psi1(l,p) );
                 
                 delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,0)*( tensor_a[0][1]*dpsidxi1(l,p,1) - tensor_a[1][1]*dpsidxi1(l,p,0));
                 
                 delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,1)*( -1*tensor_a[0][0]*dpsidxi1(l,p,1) + tensor_a[1][0]*dpsidxi1(l,p,0) );
                 
                 // for the zero derivative
                 for(unsigned k=0;k<n_dim;k++)
                  {
                   delta_kappa(al,be) -= (dtdxi(2,al,k)*dtdxi(k,be,2))*psi1(l,p);    
                  }
                 
                 // for the first derivative with respect to \al
                 delta_kappa(al,be) -= (dtdxi(2,be,2))*dpsidxi1(l,p,al);
                 
                 // for the first derivative with respect to \be
                 delta_kappa(al,be) -= (dtdxi(2,al,2))*dpsidxi1(l,p,be);
                }
              }// end of compute for variations
             
             local_eqn = this->nodal_local_eqn(l,2+p);
             if(local_eqn >= 0)
              {
               // add in external forcing     
               residuals[local_eqn] -= (1.0/H)*source[i]*psi1(l,p)*W1*sqrt(adet);
               // Loop over the Greek indicies
               for(unsigned al=0;al<2;al++)
                {
                 for(unsigned be=0;be<2;be++)
                  {
                   for(unsigned ga=0;ga<2;ga++)
                    {
                     for(unsigned de=0;de<2;de++)
                      {
                       residuals[local_eqn] += Et[al][be][ga][de]*(gamma[al][be]*delta_gamma(ga,de) + 1.0/12.0*H*H*kappa[al][be]*delta_kappa(ga,de))*W1*sqrt(adet);
                      }
                    }            
                  }
                } 
              } // end of IF
            }
          } // end of loop over n_shape
        }
       // curve sided element
       else if(bd_element==1)
        {
         //Assign storage for the local shape function of curved elements
         //The shape functions are complete polynomials of degree 7 with 36 dofs
         Shape psi_curve(36);
         DShape dpsi_curve(36,2), d2psi_curve(36,3);
         this->d2basis_eulerian_curve(integration_point,psi_curve,dpsi_curve,d2psi_curve); 
         
         //----------------
         // Rearrange nodes to be in the right position
         DenseMatrix<double> D(21,21), B(21,36), position(3,2), bd_position(20,2);
         Vector<double> x(2,0.0);
         this->get_value_transform_matrix(D,B,position,bd_element,bd_position,x);
         DenseMatrix<double> M(21,36,0.0);
         for(unsigned ii=0;ii<21;ii++)
          {
           for(unsigned j=0;j<36;j++)
            {
             for(unsigned k=0;k<21;k++)
              {
               M(ii,j) += D(ii,k)*B(k,j);
              }
            }
          }
         
         // Assign nodal equations for each nodes
         DenseMatrix<int> nodal_eqn(this->nnode(),n_position_type,-1);
         this->get_nodal_eqn_curve(nodal_eqn);
       
         //--------------------------------------------
         //The following loops is to find the residuals.
         // Loop over the test functions
         // values at vertices
         //--------------------------------------------
         // At interior nodes ei
         for(unsigned p=0;p<3;p++)
          {
           unsigned node,row;
           if(p==0){node=this->nnode()-3;row=20;}
           else if(p==1){node=this->nnode()-2;row=18;}
           else if(p==2){node=this->nnode()-1;row=19;}
           //Get the local equation
           local_eqn = nodal_eqn(node,0);
           if(local_eqn >= 0)
            {
             // Loop over the hermite test functions
             for(unsigned l=0;l<36;l++)
              {
               // initialise 
               DenseMatrix<double> delta_kappa(2,2,0.0);
               DenseMatrix<double> delta_gamma(2,2,0.0);
               for(unsigned al=0;al<2;al++)
                {
                 for(unsigned be=0;be<2;be++)
                  {
                   // compute for variations of the strain tensor
                   delta_gamma(al,be) = (dtdxi(2,be,al) + dtdxi(2,al,be))/2.0*psi_curve[l];
                   
                   // compute for variations of the bending tensor
                   if((al==0) & (be==0)) {c=0;}
                   else if((al==1) & (be==1)) {c=1;}
                   else if(al!=be) {c=2;}
                  
                   delta_kappa(al,be) -= d2psi_curve(l,c);  
                     
                   delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,0)*( dtdxi(2,1,2)*tensor_a[0][1]*psi_curve[l]-dtdxi(2,0,2)*tensor_a[1][1]*psi_curve[l] );
                   
                   delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,1)*( -1*dtdxi(2,1,2)*tensor_a[0][0]*psi_curve[l] + dtdxi(2,0,2)*tensor_a[1][0]*psi_curve[l] );
                   
                   delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,2)*( dtdxi(2,1,1)*tensor_a[0][0]*psi_curve[l]-dtdxi(2,1,0)*tensor_a[0][1]*psi_curve[l] 
                                                                    + dtdxi(2,0,0)*tensor_a[1][1]*psi_curve[l]-dtdxi(2,0,1)*tensor_a[1][0]*psi_curve[l] );
                   
                   delta_kappa(al,be) += dtdxi(al,be,2)/adet/adet*( dtdxi(2,1,1)*tensor_a[0][0]*psi_curve[l]-dtdxi(2,1,0)*tensor_a[0][1]*psi_curve[l] 
                                                                    + dtdxi(2,0,0)*tensor_a[1][1]*psi_curve[l]-dtdxi(2,0,1)*tensor_a[1][0]*psi_curve[l] );
                   
                   delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,0)*( tensor_a[0][1]*dpsi_curve(l,1) - tensor_a[1][1]*dpsi_curve(l,0));
                   
                   delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,1)*( -1*tensor_a[0][0]*dpsi_curve(l,1) + tensor_a[1][0]*dpsi_curve(l,0) );
                     
                   // for the zero derivative
                   for(unsigned k=0;k<n_dim;k++)
                    {
                     delta_kappa(al,be) -= (dtdxi(2,al,k)*dtdxi(k,be,2))*psi_curve[l];    
                    }
                   
                   // for the first derivative with respect to \al
                   delta_kappa(al,be) -= (dtdxi(2,be,2))*dpsi_curve(l,al);
                   
                   // for the first derivative with respect to \be
                   delta_kappa(al,be) -= (dtdxi(2,al,2))*dpsi_curve(l,be);
                  }
                }// end of compute for variations
               
               if(M(row,l)!=0)
                {
                 // add in external forcing     
                 residuals[local_eqn] -= (1.0/H)*source[i]*psi_curve[l]*W1*sqrt(adet)*M(row,l);
                 // Loop over the Greek indicies
                 for(unsigned al=0;al<2;al++)
                  {
                   for(unsigned be=0;be<2;be++)
                    {
                     for(unsigned ga=0;ga<2;ga++)
                      {
                       for(unsigned de=0;de<2;de++)
                        {
                         residuals[local_eqn] += Et[al][be][ga][de]*(gamma[al][be]*delta_gamma(ga,de) + 1.0/12.0*H*H*kappa[al][be]*delta_kappa(ga,de))*W1*sqrt(adet)*M(row,l); 
                        }
                      }            
                    }
                  } 
                } 
              }
            }
          } 
         //-------------------
         // At vertices ai
         for(unsigned m=0;m<3;m++)
          {
           for(unsigned type=0;type<6;type++)
            {
             unsigned row;
             if(type==0){row = m;}
             else if((type==1) || (type==2)){row = 2*(m+1)+type;}
             else if(type==3){row = 3*(m+2)+3;}
             else if(type==4){row = 3*(m+2)+5;}
             else if(type==5){row = 3*(m+2)+4;}
             //Get the local equation
             local_eqn = nodal_eqn(m,type);
             //IF it's not a boundary condition
             if(local_eqn >= 0)
              {
               for(unsigned l=0;l<36;l++)
                {

                 DenseMatrix<double> delta_kappa(2,2,0.0);
                 DenseMatrix<double> delta_gamma(2,2,0.0);
                 for(unsigned al=0;al<2;al++)
                  {
                   for(unsigned be=0;be<2;be++)
                    {
                     // compute for variations of the strain tensor
                     delta_gamma(al,be) = (dtdxi(2,be,al) + dtdxi(2,al,be))/2.0*psi_curve[l];
                     
                     // compute for variations of the bending tensor
                     if((al==0) & (be==0)) {c=0;}
                     else if((al==1) & (be==1)) {c=1;}
                     else if(al!=be) {c=2;}
                     
                     delta_kappa(al,be) -= d2psi_curve(l,c);  
                     
                     delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,0)*( dtdxi(2,1,2)*tensor_a[0][1]*psi_curve[l]-dtdxi(2,0,2)*tensor_a[1][1]*psi_curve[l] );
                     
                     delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,1)*( -1*dtdxi(2,1,2)*tensor_a[0][0]*psi_curve[l] + dtdxi(2,0,2)*tensor_a[1][0]*psi_curve[l] );
                     
                     delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,2)*( dtdxi(2,1,1)*tensor_a[0][0]*psi_curve[l]-dtdxi(2,1,0)*tensor_a[0][1]*psi_curve[l] 
                                                                      + dtdxi(2,0,0)*tensor_a[1][1]*psi_curve[l]-dtdxi(2,0,1)*tensor_a[1][0]*psi_curve[l] );
                     
                     delta_kappa(al,be) += dtdxi(al,be,2)/adet/adet*( dtdxi(2,1,1)*tensor_a[0][0]*psi_curve[l]-dtdxi(2,1,0)*tensor_a[0][1]*psi_curve[l] 
                                                                      + dtdxi(2,0,0)*tensor_a[1][1]*psi_curve[l]-dtdxi(2,0,1)*tensor_a[1][0]*psi_curve[l] );
                     
                     delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,0)*( tensor_a[0][1]*dpsi_curve(l,1) - tensor_a[1][1]*dpsi_curve(l,0));
                     
                     delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,1)*( -1*tensor_a[0][0]*dpsi_curve(l,1) + tensor_a[1][0]*dpsi_curve(l,0) );
                     
                     // for the zero derivative
                     for(unsigned k=0;k<n_dim;k++)
                      {
                       delta_kappa(al,be) -= (dtdxi(2,al,k)*dtdxi(k,be,2))*psi_curve[l];    
                      }
                     
                     // for the first derivative with respect to \al
                     delta_kappa(al,be) -= (dtdxi(2,be,2))*dpsi_curve(l,al);
                     
                     // for the first derivative with respect to \be
                     delta_kappa(al,be) -= (dtdxi(2,al,2))*dpsi_curve(l,be);
                    }
                  }// end of compute for variations
             
                 if(M(row,l)!=0)
                  {
                   // add in external forcing     
                   residuals[local_eqn] -= (1.0/H)*source[i]*psi_curve[l]*W1*sqrt(adet)*M(row,l);
                   // Loop over the Greek indicies
                   for(unsigned al=0;al<2;al++)
                    {
                     for(unsigned be=0;be<2;be++)
                      {
                       for(unsigned ga=0;ga<2;ga++)
                        {
                         for(unsigned de=0;de<2;de++)
                          {
                           residuals[local_eqn] += Et[al][be][ga][de]*(gamma[al][be]*delta_gamma(ga,de) + 1.0/12.0*H*H*kappa[al][be]*delta_kappa(ga,de))*W1*sqrt(adet)*M(row,l); 
                          }
                        }
                      }            
                    }
                  } 
                }
              } 
            }
          } // end of loop over n_shape
         //-------------------
        }
      }
     else /// find residuals for the first and second directions(tangents)
      {
       // Assign nodal equations for each nodes
       DenseMatrix<int> nodal_eqn(this->nnode(),n_position_type,-1);
       if(bd_element==1)
        {
         this->get_nodal_eqn_curve(nodal_eqn);
        }
       
       // Loop over the test functions
       for(unsigned l=0;l<this->nnode()-3;l++)
        {
         // initialise 
         DenseMatrix<double> delta_kappa(2,2,0.0);
         DenseMatrix<double> delta_gamma(2,2,0.0);
         
         // compute for variations of the strain tensor
         if(i==0)
          {
           delta_gamma(0,0) = dpsidxi(l,0) + dtdxi(i,0,0)*psi[l];
           delta_gamma(0,1) = dpsidxi(l,1)/2.0 + (dtdxi(i,1,0) + dtdxi(i,0,1))/2.0*psi[l];
           delta_gamma(1,0) = dpsidxi(l,1)/2.0 + (dtdxi(i,1,0) + dtdxi(i,0,1))/2.0*psi[l];
           delta_gamma(1,1) = dtdxi(i,1,1)*psi[l];
          }
         else if(i==1)
          {
           delta_gamma(0,0) = dtdxi(i,0,0)*psi[l];
           delta_gamma(0,1) = dpsidxi(l,0)/2.0 + (dtdxi(i,1,0) + dtdxi(i,0,1))/2.0*psi[l];
           delta_gamma(1,0) = dpsidxi(l,0)/2.0 + (dtdxi(i,1,0) + dtdxi(i,0,1))/2.0*psi[l];
           delta_gamma(1,1) = dpsidxi(l,1) + dtdxi(i,1,1)*psi[l]; 
          }
         
         // compute for variations of the bending tensor
         for(unsigned al=0;al<2;al++)
          {
           for(unsigned be=0;be<2;be++)
            {
             // compute for variations of the bending tensor
             if((al==0) & (be==0)) {c=0;}
             else if((al==1) & (be==1)) {c=1;}
             else if(al!=be) {c=2;}
             
             if(i==0)
              {
               delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,0)*( dtdxi(0,1,2)*tensor_a[0][1]*psi[l]-dtdxi(0,0,2)*tensor_a[1][1]*psi[l] );
               
               delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,1)*( -1*dtdxi(0,1,2)*tensor_a[0][0]*psi[l] + dtdxi(0,0,2)*tensor_a[1][0]*psi[l] );
               
               delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,2)*( dtdxi(0,1,1)*tensor_a[0][0]*psi[l]-dtdxi(0,1,0)*tensor_a[0][1]*psi[l] 
                                                                + dtdxi(0,0,0)*tensor_a[1][1]*psi[l]-dtdxi(0,0,1)*tensor_a[1][0]*psi[l] );
               
               delta_kappa(al,be) += dtdxi(al,be,2)/adet/adet*( dtdxi(0,1,1)*tensor_a[0][0]*psi[l]-dtdxi(0,1,0)*tensor_a[0][1]*psi[l] 
                                                                + dtdxi(0,0,0)*tensor_a[1][1]*psi[l]-dtdxi(0,0,1)*tensor_a[1][0]*psi[l] );
               
               delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,2)*( -1*tensor_a[0][1]*dpsidxi(l,1) + tensor_a[1][1]*dpsidxi(l,0));
               
               delta_kappa(al,be) += dtdxi(al,be,2)/adet/adet*( -1*tensor_a[0][1]*dpsidxi(l,1) + tensor_a[1][1]*dpsidxi(l,0) );
               
              }
             else if(i==1)
              {
               delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,0)*( dtdxi(1,1,2)*tensor_a[0][1]*psi[l]-dtdxi(1,0,2)*tensor_a[1][1]*psi[l] );
               
               delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,1)*( -1*dtdxi(1,1,2)*tensor_a[0][0]*psi[l] + dtdxi(1,0,2)*tensor_a[1][0]*psi[l] );
               
               delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,2)*( dtdxi(1,1,1)*tensor_a[0][0]*psi[l]-dtdxi(1,1,0)*tensor_a[0][1]*psi[l] 
                                                                + dtdxi(1,0,0)*tensor_a[1][1]*psi[l]-dtdxi(1,0,1)*tensor_a[1][0]*psi[l] );
               
               delta_kappa(al,be) += dtdxi(al,be,2)/adet/adet*( dtdxi(1,1,1)*tensor_a[0][0]*psi[l]-dtdxi(1,1,0)*tensor_a[0][1]*psi[l] 
                                                                + dtdxi(1,0,0)*tensor_a[1][1]*psi[l]-dtdxi(1,0,1)*tensor_a[1][0]*psi[l] );
               
               delta_kappa(al,be) += -1.0/adet*dtdxi(al,be,2)*( tensor_a[0][0]*dpsidxi(l,1) - tensor_a[1][0]*dpsidxi(l,0));
               
               delta_kappa(al,be) += dtdxi(al,be,2)/adet/adet*( tensor_a[0][0]*dpsidxi(l,1) - tensor_a[1][0]*dpsidxi(l,0) );
              }

             // for the zero derivative
             for(unsigned k=0;k<n_dim;k++)
              {
               delta_kappa(al,be) -= (dtdxi(i,al,k)*dtdxi(k,be,2))*psi[l];    
              }
             
             // for the first derivative with respect to \al
             delta_kappa(al,be) -= (dtdxi(i,be,2))*dpsidxi(l,al);
             
             // for the first derivative with respect to \be
             delta_kappa(al,be) -= (dtdxi(i,al,2))*dpsidxi(l,be);
            }
          }   // end of compute for variations
         
         //Get the local equation
         local_eqn = this->nodal_local_eqn(l,i);
         
         if(local_eqn >= 0)
          {
           // add in external forcing     
           residuals[local_eqn] -= (1.0/H)*source[i]*psi[l]*W*sqrt(adet);
           // Loop over the Greek indicies
           for(unsigned al=0;al<2;al++)
            {
             
             for(unsigned be=0;be<2;be++)
              {
               for(unsigned ga=0;ga<2;ga++)
                {
                 for(unsigned de=0;de<2;de++)
                  {
                   residuals[local_eqn] += Et[al][be][ga][de]*(gamma[al][be]*delta_gamma(ga,de) + 1.0/12.0*H*H*kappa[al][be]*delta_kappa(ga,de))*W1*sqrt(adet); 
                  }
                }
              }            
            }
          } // End of if  
        }
      }// end of else if (of direction)
    }
  } // End of loop over integration points 
}   



//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
unsigned  MyShellEquations<DIM,NNODE_1D>::self_test()
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
void  MyShellEquations<DIM,NNODE_1D>::output(std::ostream &outfile, 
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
   u = interpolated_u_shell(s);
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

 // Write tecplot footer (e.g. FE connectivity lists)
 //write_tecplot_zone_footer(outfile,nplot);
 
}

 
//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void  MyShellEquations<DIM,NNODE_1D>::output(FILE* file_pt,
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
   u = interpolated_u_shell(s);
   fprintf(file_pt,"%g \n",u[0]);//interpolated_u_poisson(s));
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
void MyShellEquations<DIM,NNODE_1D>::output_fct(std::ostream &outfile, 
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
    outfile << exact_soln[j] << "\t";
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
void MyShellEquations<DIM,NNODE_1D>::compute_error(std::ostream &outfile, 
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
   double J=this->J_eulerian(s);
  
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   // Get x position as Vector
   this->my_interpolated_x(s,x);
   
   // Get FE function value
   Vector<double> u_fe(this->required_nvalue(0),0.0);
   u_fe =interpolated_u_shell(s);
   
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
   double tmp1 = 0.0, tmp2 =0.0;
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

template class C1CurvedShellElement<2,2>;
template class C1CurvedShellElement<2,3>;
}






//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//==start_of_namespace================================================
/// Namespace for solutions of 2D Linear shell equation
//====================================================================
namespace Physical_Variables
{
 /// Pressure load
 double P_ext;
 
 double epsilon = 1.0e-7;
 
 /// Exact solution
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = 0.0;
  u[1] = 0.0;
  u[2] = 0.0;
  u[3] = 0.0;
  u[4] = 0.0;
  u[5] = 0.0;
  u[6] = 0.0;
  u[7] = 0.0;
 }

 /// Source function 
 void source_function(const Vector<double>& x, const Vector<double>& unit_n, Vector<double>& source)
 {
  for(unsigned i=0;i<3;i++)
   {
    source[i] = 1.0*epsilon*P_ext*unit_n[i];
   }
 }

} // end of namespace


//==start_of_problem_class============================================
/// 2D linearised shell problem.
//====================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D> 
class MyLinearisedShellProblem : public Problem
{

public:

 /// Constructor: Pass number of elements and pointer to source function
 MyLinearisedShellProblem(typename MyShellEquations<DIM,NNODE_1D>::SourceFctPt source_fct_pt,
                          typename MyShellEquations<DIM,NNODE_1D>::ExactSolnPt exact_pt,
                          const string& node_file_name,
                          const string& element_file_name,
                          const string& poly_file_name);
 /// Destructor (empty)
 ~MyLinearisedShellProblem()
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
 
 void parameter_study();

private:

 /// Pointer to source function
 typename MyShellEquations<DIM,NNODE_1D>::SourceFctPt Source_fct_pt;
 
  /// Pointer to the exact solution
 typename MyShellEquations<DIM,NNODE_1D>::ExactSolnPt Exact_soln_pt;

 /// Pointer to geometric object that represents the shell's undeformed shape
 GeomObject* Undef_midplane_pt;

}; // end of problem class

//=====start_of_constructor===============================================
/// \short Constructor for 2D Shell problem.
/// Discretise the 2D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function. 
//========================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D>
MyLinearisedShellProblem<ELEMENT,DIM,NNODE_1D>::MyLinearisedShellProblem
(typename MyShellEquations<DIM,NNODE_1D>::SourceFctPt source_fct_pt, 
 typename MyShellEquations<DIM,NNODE_1D>::ExactSolnPt exact_pt,
 const string& node_file_name,
 const string& element_file_name,
 const string& poly_file_name) : 
 Source_fct_pt(source_fct_pt), Exact_soln_pt(exact_pt)
{ 
 // set the undeformed shell
 Undef_midplane_pt = new CircularPlate(0.0,0.0);

 
 // Build mesh and store pointer in Problem
 Problem::mesh_pt() = new TriangleMesh<ELEMENT>(node_file_name,element_file_name,poly_file_name);

 unsigned n_node = this->mesh_pt()->nnode();
 cout << "\n\nNumber of nodes in the mesh = " << n_node << endl;
 
 // find number of elements in the mesh
 unsigned n_element = mesh_pt()->nelement();
 cout  << "Number of nelements in the mesh = " << n_element << endl;
 
 // pinned for the middle node in each element for normal direction
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
  // pinned the unnecessary nodes
  for(unsigned i=0;i<nnode;i++)
   {
    if(i>2 && i<(nnode-3))
     {
      elem_pt->node_pt(i)->pin(2);
      elem_pt->node_pt(i)->pin(3);
      elem_pt->node_pt(i)->pin(4);
      elem_pt->node_pt(i)->pin(5); 
      elem_pt->node_pt(i)->pin(6);
      elem_pt->node_pt(i)->pin(7);
     }
    else if(i>(nnode-4))
     {
      if(count == 2)
       {
        elem_pt->node_pt(i)->pin(0);
        elem_pt->node_pt(i)->pin(1);
        //elem_pt->node_pt(i)->pin(2); // supposed to be free
        elem_pt->node_pt(i)->pin(3);
        elem_pt->node_pt(i)->pin(4);
        elem_pt->node_pt(i)->pin(5);
        elem_pt->node_pt(i)->pin(6);
        elem_pt->node_pt(i)->pin(7);
       }
      else
       {
        elem_pt->node_pt(i)->pin(0);
        elem_pt->node_pt(i)->pin(1);
        elem_pt->node_pt(i)->pin(2);
        elem_pt->node_pt(i)->pin(3);
        elem_pt->node_pt(i)->pin(4);
        elem_pt->node_pt(i)->pin(5);
        elem_pt->node_pt(i)->pin(6);
        elem_pt->node_pt(i)->pin(7);
       }
     }
   }
 }

 // Set the boundary conditions for this problem: By default, all nodal
 // values are free -- we only need to pin the ones that have 
 // Dirichlet conditions. 
 unsigned n_side0 = mesh_pt()->nboundary_node(0);
 unsigned n_side1 = mesh_pt()->nboundary_node(1);
 unsigned n_side2 = mesh_pt()->nboundary_node(2);
   
 // Pin the single nodal value at the single node on mesh 
 // boundary 0 (= the left domain boundary at x=0)
 /// loop over the nodes on the boundary
 for(unsigned i=0;i<n_side0;i++)
  {
   //mesh_pt()->boundary_node_pt(0,i)->pin(0); //
   mesh_pt()->boundary_node_pt(0,i)->pin(1); 
   //mesh_pt()->boundary_node_pt(0,i)->pin(2);  //
   
   //mesh_pt()->boundary_node_pt(0,i)->pin(3); //
   mesh_pt()->boundary_node_pt(0,i)->pin(4); 
   //mesh_pt()->boundary_node_pt(0,i)->pin(5); //
   mesh_pt()->boundary_node_pt(0,i)->pin(6); 
   mesh_pt()->boundary_node_pt(0,i)->pin(7); //
  }

 // Pin the single nodal value at the single node on mesh 
 // boundary 1 (= the right domain boundary at x=1)
 for(unsigned i=0;i<n_side1;i++)
  {
   mesh_pt()->boundary_node_pt(1,i)->pin(0); 
   mesh_pt()->boundary_node_pt(1,i)->pin(1); 
   mesh_pt()->boundary_node_pt(1,i)->pin(2); 

   mesh_pt()->boundary_node_pt(1,i)->pin(3);
   mesh_pt()->boundary_node_pt(1,i)->pin(4); 
   mesh_pt()->boundary_node_pt(1,i)->pin(5); 
   mesh_pt()->boundary_node_pt(1,i)->pin(6);
   mesh_pt()->boundary_node_pt(1,i)->pin(7); 
  }

 for(unsigned i=0;i<n_side2;i++)
  {
   mesh_pt()->boundary_node_pt(2,i)->pin(0); 
   //mesh_pt()->boundary_node_pt(2,i)->pin(1); //
   //mesh_pt()->boundary_node_pt(2,i)->pin(2); //
   
   mesh_pt()->boundary_node_pt(2,i)->pin(3); 
   //mesh_pt()->boundary_node_pt(2,i)->pin(4); //
   mesh_pt()->boundary_node_pt(2,i)->pin(5); 
   //mesh_pt()->boundary_node_pt(2,i)->pin(6); //
   mesh_pt()->boundary_node_pt(2,i)->pin(7); //
  } 

 // Loop over elements and set pointers to Physical parameters
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   //Set the source function pointer and all physical variables
   elem_pt->source_fct_pt() = Source_fct_pt;
   
   elem_pt->exact_pt() = Exact_soln_pt;

   elem_pt->undeformed_midplane_pt() = Undef_midplane_pt;
  }

 //--------------------------------
 // Setup equation numbering scheme
 assign_eqn_numbers();
 
} // end of constructor




//===start_of_actions_before_newton_solve========================================
/// \short Update the problem specs before solve: (Re)set boundary values
/// from the exact solution. 
//========================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D>
void MyLinearisedShellProblem<ELEMENT,DIM,NNODE_1D>::actions_before_newton_solve()
{
 
 // Assign boundary values for this problem by reading them out
 // from the exact solution.
 for(unsigned n=0;n<mesh_pt()->nboundary();n++)
  {
   // find number of nodes in each boundary
   unsigned n_node = mesh_pt()->nboundary_node(n);
   /// loop over the nodes on the boundary
   for(unsigned j=0;j<n_node;j++)
    {
     // Left boundary is every nodes on the left boundary 
     Node* node_pt=mesh_pt()->boundary_node_pt(n,j);
   
     // Loop for variables u
     ELEMENT e;
     Vector<double> u((e.required_nvalue(0)));
     for(unsigned i=0;i<(e.required_nvalue(0));i++)
      {
       // Determine the position of the boundary node (the exact solution
       // requires the coordinate in a 1D vector!)
       Vector<double> x(2);
       x[0]=node_pt->x(0);
       x[1]=node_pt->x(1);
       
       // Boundary value (read in from exact solution which returns
       // the solution in a 1D vector)
       Physical_Variables::get_exact_u(x,u);
  
       // Assign the boundary condition to one (and only) nodal value
       node_pt->set_value(i,u[i]);
      }
    }
  }
} // end of actions before solve


//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D>
void MyLinearisedShellProblem<ELEMENT,DIM,NNODE_1D>::doc_solution(DocInfo& doc_info)
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
} // end of doc

//===start_of_parameter_study=========================================================
/// Solver loop to perform parameter study
//========================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D>
void MyLinearisedShellProblem<ELEMENT,DIM,NNODE_1D>::parameter_study()
{ 
 // over-ride the default maximum value for the residuals
 Problem::Max_residuals = 1.0e10;
 
 unsigned nstep = 3;
 
 // set the increments in control parameters
 double pext_increment = 1.0;
 // set initial values for control parameters
 Physical_Variables::P_ext = 0.0 - pext_increment;

 // Set up doc info
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 
 // Loop over parameter increments
 for(unsigned i=0;i<nstep;i++)
   {
    Physical_Variables::P_ext += pext_increment;
    
    // solve the system
    newton_solve();
    
    //Output the solution
    doc_info.number() = i;
    doc_solution(doc_info);    
   }
} //end of parameter study

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//======start_of_main==================================================
/// Driver for 2D shell problem
//=====================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

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

 // Set up the problem: Element type as template parameter
 MyLinearisedShellProblem<C1CurvedShellElement<2,3>,2,3> 
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
 
 problem.parameter_study();

} // end of main









