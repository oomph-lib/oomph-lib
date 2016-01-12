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


//Header file for 2d linear shell problem: square-plate bending
#ifndef OOMPH_LINEAR_SHELL_ELEMENTS_HEADER
#define OOMPH_LINEAR_SHELL_ELEMENTS_HEADER

namespace oomph
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Geometric object
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//===========================================================
/// \short Elliptical tube with half axes a and b.
/// 
/// \f[ {\bf r} = ( a \cos(\zeta_1), b \sin(zeta_1), \zeta_0)^T \f]
/// 
//===========================================================
class Plate : public GeomObject
{
public:

 /// Constructor: Specify radius
 Plate(const double& a, const double& b) : 
  GeomObject(2,3), A(a), B(b) {}
 
 /// Broken copy constructor
 Plate(const Plate& node) 
  { 
   BrokenCopy::broken_copy("Plate");
  } 
 
 /// Broken assignment operator
 void operator=(const Plate&) 
  {
   BrokenCopy::broken_assign("Plate");
  }

 /// Access function to x-half axis
 double& a() {return A;}

 /// Access function to y-half axis
 double& b() {return B;}

 /// Position vector
 void position(const Vector<double>& zeta, Vector<double>& r)const
  {
   r[0]=zeta[1];
   r[1]=0.0;
   r[2]=zeta[0];
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
   r[0]=zeta[1];
   r[1]=0.0;
   r[2]=zeta[0];

   //Do the azetaal derivatives
   drdzeta(0,0) = 0.0;
   drdzeta(0,1) = 0.0;
   drdzeta(0,2) = 1.0;
   
   //Do the azimuthal derivatives
   drdzeta(1,0) = 1.0;
   drdzeta(1,1) = 0.0;
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

 /// x-half axis
 double A;

 /// x-half axis
 double B;

};


//=============================================================
/// A class for all subparametric elements that solve the 
/// linear shell equations.
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
template <unsigned DIM, unsigned NNODE_1D>
class MyShellEquations : public virtual BellElement<DIM,NNODE_1D>
{

public:
 
 /// \short Function pointer to source function fct(x,f(x)) -- 
 /// x is a Vector! 
 typedef void (*SourceFctPt)(const Vector<double>& x, const Vector<double>& unit_n, Vector<double>& f);


 /// \short Function pointer to gradient of source function  fct(x,g(x)) -- 
 /// x is a Vector! 
 typedef void (*SourceFctGradientPt)(const Vector<double>& x, 
                                            Vector<double>& gradient);


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
    "There is no time-dependent output_fct() for linear shell elements ",
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
    "There is no time-dependent compute_error() for shell elements",
    "MyShellEquations<DIM>::compute_error()",
    OOMPH_EXCEPTION_LOCATION);
  }

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

 /// Get source term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the source function might be determined by
 /// another system of equations.
 inline virtual void get_source_function(const unsigned& ipt,
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
 //void fill_in_contribution_to_jacobian(Vector<double> &residuals,
 //                                  DenseMatrix<double> &jacobian)
 // {
   //Call the generic routine with the flag set to 1
  // fill_in_generic_residual_contribution_shell(residuals,jacobian,1);
  //}
 


 /// \short Return FE representation of function value u(s) 
 /// at local coordinate s
 inline Vector<double> interpolated_u_shell(const Vector<double> &s) const
  {
   //Find number of nodes
   const unsigned n_node = this->nnode();
   
   // Find number of dimension for eulerian coordinates
   const unsigned n_dim = Undeformed_midplane_pt->ndim();

   // Find number of dimension for eulerian coordinates
   const unsigned n_lagrangian = Undeformed_midplane_pt->nlagrangian();

   //Find number of position dofs
   const unsigned n_position_type = this->nnodal_position_type();
   
   //Get the index at which the unknown is stored
   const unsigned u_nodal_index = u_index_shell();
   
   //Local shape function
   Shape psi1(3,n_position_type),psi(n_node);
   DShape dpsidxi1(3,n_position_type,DIM);
   DShape d2psidxi1(3,n_position_type,3);

   //Initialise value of u
   Vector<double> interpolated_u(u_nodal_index,0.0);
   this->basis_c0(s,psi);
   
   //Interpolated in tangential direction
   for(unsigned l=0;l<n_node;l++) 
    {
     interpolated_u[0] += this->nodal_value(l,0)*psi[l];          
     interpolated_u[1] += this->nodal_value(l,1)*psi[l];  
    }

   //Find values of shape function
   this->d2basis_local(s,psi1,dpsidxi1,d2psidxi1);
   
   //Interpolated in normal direction
   for(unsigned l=0;l<3;l++) 
   {
    for(unsigned k=0;k<n_position_type;k++)
     {
      interpolated_u[2] += this->nodal_value(l,2+k)*psi1(l,k);
      interpolated_u[3] += this->nodal_value(l,2+k)*dpsidxi1(l,k,0);
      interpolated_u[4] += this->nodal_value(l,2+k)*dpsidxi1(l,k,1);
      interpolated_u[5] += this->nodal_value(l,2+k)*d2psidxi1(l,k,0);
      interpolated_u[6] += this->nodal_value(l,2+k)*d2psidxi1(l,k,1);
      interpolated_u[7] += this->nodal_value(l,2+k)*d2psidxi1(l,k,2);
     }
   } // End of varible loop
   
    
   // setup position vector and derivatives of undeformed configuration
   Vector<double> r(n_dim);
   DenseMatrix<double> a(n_lagrangian,n_dim);
   RankThreeTensor<double> dadxi(n_lagrangian,n_lagrangian,n_dim);
   
   Vector<double> interpolated_xi(DIM,0.0);
   this->my_interpolated_x(s,interpolated_xi);

   // get the undeformed geometry
   Undeformed_midplane_pt->d2position(interpolated_xi,r,a,dadxi);
   
   // calculate the covariant metric tensors
   double tensor_a[2][2];

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
   
   // calculate the contravariant metric tensor
   //Calculate determinant
   double adet = tensor_a[0][0]*tensor_a[1][1] - tensor_a[0][1]*tensor_a[1][0];
   
   // calculate the unit normal vectors
   Vector<double> unit_n(n_dim);
   unit_n[0] = 1.0/sqrt(adet)*(a(0,1)*a(1,2) - a(0,2)*a(1,1));
   unit_n[1] = 1.0/sqrt(adet)*(a(0,2)*a(1,0) - a(0,0)*a(1,2));
   unit_n[2] = 1.0/sqrt(adet)*(a(0,0)*a(1,1) - a(0,1)*a(1,0));

   // calculate the non-unit tangential vectors
   DenseMatrix<double> t(n_lagrangian,n_dim);
   for(unsigned i=0;i<n_lagrangian;i++)
    {
     for(unsigned j=0;j<n_dim;j++)
      {
       t(i,j) = a(i,j);            
      }            
    }
    
   // /// compute displacement u in cartesian coordinate
   // /// displacement component in y direction is identical with normal dirention
   // double x_dir;
   // double y_dir;
   // double z_dir;
   
   // x_dir = interpolated_u[0]*t(0,0) + interpolated_u[1]*t(1,0) + interpolated_u[2]*unit_n[0];
   // y_dir = interpolated_u[0]*t(0,1) + interpolated_u[1]*t(1,1) + interpolated_u[2]*unit_n[1];
   // z_dir = interpolated_u[0]*t(0,2) + interpolated_u[1]*t(1,2) + interpolated_u[2]*unit_n[2];
   
   /*interpolated_u[0] = x_dir;
   interpolated_u[1] = y_dir;
   interpolated_u[2] = z_dir;
   */
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
/// BellShellElement elements are with subparametric interpolation for the function.
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
class BellShellElement : public virtual MyShellEquations<DIM,NNODE_1D>
{

private:

 /// \short Static int that holds the number of variables at 
 /// nodes: always the same
 static const unsigned Initial_Nvalue;
 
  public:


 ///\short  Constructor: Call constructors for BellElement and 
 /// Shell equations
 BellShellElement() : MyShellEquations<DIM,NNODE_1D>()
  {}
 
 /// Broken copy constructor
 BellShellElement(const BellShellElement<DIM,NNODE_1D>& dummy) 
  { 
   BrokenCopy::broken_copy("BellShellElement");
  } 
 
 /// Broken assignment operator
 void operator=(const BellShellElement<DIM,NNODE_1D>&) 
  {
   BrokenCopy::broken_assign("BellShellElement");
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


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


//=======================================================================
/// Face geometry for the BellShellElement elements: The spatial 
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<BellShellElement<DIM,NNODE_1D> >: 
 public virtual TElement<DIM-1,NNODE_1D>
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional TElement
 FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}

};


//Inline functions:


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
 double BellShellElement<DIM,NNODE_1D>::dshape_and_dtest_eulerian_shell(
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
 double BellShellElement<DIM,NNODE_1D>::d2shape_and_d2test_eulerian_shell(
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
double BellShellElement<DIM,NNODE_1D>::
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
double BellShellElement<DIM,NNODE_1D>::
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
 const unsigned BellShellElement<DIM,NNODE_1D>::Initial_Nvalue = 8;

//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void  MyShellEquations<DIM,NNODE_1D>::
fill_in_generic_residual_contribution_shell(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              const unsigned& flag) 
{
 //Find out how many nodes there are
 const unsigned n_node = this->nnode();
 
 //Find out how many nodes positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();
 
 //Set the dimension of the global coordinates
 unsigned n_dim = Undeformed_midplane_pt->ndim();
  
 //Set the number of lagrangian coordinates
 unsigned n_lagrangian =  Undeformed_midplane_pt->nlagrangian();
  
 //Set up memory for the shape and test functions
 Shape psi1(3,n_position_type), test1(3,n_position_type), test(n_node), psi(n_node);
 DShape dpsidxi1(3,n_position_type,DIM), dtestdxi1(3,n_position_type,DIM), dtestdxi(n_node,DIM), dpsidxi(n_node,DIM);
 DShape d2psidxi1(3,n_position_type,3), d2testdxi1(3,n_position_type,3), d2testdxi(n_node,3), d2psidxi(n_node,3);
 
 //Set the value of n_intpt
 const unsigned n_intpt = this->integral_pt()->nweight();

 //Integers to store the local equation and unknown numbers
 int local_eqn=0;//, local_unknown=0;
      
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);
  

   //Call the derivatives of the c0-shape and test functions for tangential direcion
   double J = dshape_and_dtest_eulerian_at_knot_shell(ipt,psi,dpsidxi,test,dtestdxi);
   
   //Call the derivatives of the c1-shape and test functions for normal direction
   double J1 = d2shape_and_d2test_eulerian_at_knot_shell(ipt,psi1,dpsidxi1,d2psidxi1,test,dtestdxi,d2testdxi);

   //Premultiply the weights and the Jacobian
   double W1 = w*J1;
   double W = w*J;
   
   double H = 0.01;
   
   //Calculate local values of unknown
   Vector<double> u_value(n_position_type,0.0);
   Vector<double> interpolated_u(3,0.0);
   DenseMatrix<double> interpolated_dudxi(3,DIM,0.0);
   DenseMatrix<double> interpolated_d2udxi(3,3,0.0);
   Vector<double> interpolated_xi(DIM,0.0);
   
   Vector<double>  s(2);
   s[0] = this->integral_pt()->knot(ipt,0);
   s[1] = this->integral_pt()->knot(ipt,1);
   this->my_interpolated_x(s,interpolated_xi);

   
   //Calculate displacement value and derivatives:
   //-----------------------------------------
   //in both tangential directions
   for(unsigned i=0;i<2;i++) 
    {
     //Loop over lagrangian coordinate directions
     for(unsigned l=0;l<n_node;l++)
      {
        //Get the nodal value of the unknown
        u_value[i] = this->raw_nodal_value(l,i);
        interpolated_u[i] += u_value[i]*psi[l];
        
        // Loop over the first derivative directions
        for(unsigned j=0;j<n_lagrangian;j++)
         {
          interpolated_dudxi(i,j) += u_value[i]*dpsidxi(l,j);
         }
        
        // Loop over the second derivative directions
        for(unsigned j=0;j<3;j++)
         {
          interpolated_d2udxi(i,j) += u_value[i]*d2psidxi(l,j);            
         }
      }        
    }
           
   // in normal direction
   for(unsigned l=0;l<3;l++) 
    {
     for(unsigned k=0;k<n_position_type;k++)
      {
       //Get the nodal value of the unknown
       u_value[k] = this->raw_nodal_value(l,2+k);
       interpolated_u[2] += u_value[k]*psi1(l,k);
       
       // Loop over directions
       for(unsigned j=0;j<n_lagrangian;j++)
        {
         interpolated_dudxi(2,j) += u_value[k]*dpsidxi1(l,k,j);
        }
        
       // Loop over the second derivative directions
       for(unsigned j=0;j<3;j++)
        {
         interpolated_d2udxi(2,j) += u_value[k]*d2psidxi1(l,k,j);        
        }
      }
    } 
    
   //-------------------------------------------------------
   // setup position vector and derivatives of undeformed configuration
   Vector<double> r(n_dim);
   DenseMatrix<double> a(n_lagrangian,n_dim), a_tn(n_lagrangian,n_dim);
   RankThreeTensor<double> dadxi(n_lagrangian,n_lagrangian,n_dim);
   
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
   aup[0][1] = -1*tensor_a[0][1]/adet;
   aup[1][0] = -1*tensor_a[1][0]/adet;
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
   
   // compute for the derivative of the shell normal vector
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
   
   
   // compute for the tangential gradient in tangential and normal directions -> Christoffel symbol
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
    for(unsigned k=0;k<n_dim;k++)
     {
      interpolated_d2udxi_tn(k,c) = interpolated_d2udxi(k,c);

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
   
   unsigned c=0;
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
   get_source_function(ipt,interpolated_xi,unit_n,source);
   // decompose in the tangential and normal coordinates
   f[0] = source[0]*t(0,0) + source[1]*t(0,1) + source[2]*t(0,2);
   f[1] = source[0]*t(1,0) + source[1]*t(1,1) + source[2]*t(1,2);
   f[2] = source[0]*t(2,0) + source[1]*t(2,1) + source[2]*t(2,2);
   source = f;
   //--------------------------------
   // Assemble residuals and Jacobian
   //--------------------------------
   // Loop over directions
   for(unsigned i=0;i<n_dim;i++)
    {
     // compute for variations
     if(i==2)
      {
       // Loop over the c1-test functions
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
              } // End of if
            } // end of loop over position type
          }
        } // end of loop over n_shape
      }

     else /// find residuals for the first and second directions
      {
       // Loop over the test functions
       for(unsigned l=0;l<n_node;l++)
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
 //outfile << this->tecplot_zone_string(nplot);
 
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
 //this->write_tecplot_zone_footer(outfile,nplot);
 
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
 outfile << this->tecplot_zone_string(nplot);
 
 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(this->required_nvalue(0));
 
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
 this->write_tecplot_zone_footer(outfile,nplot);
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
 outfile << "ZONE" << std::endl;
 
 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(this->required_nvalue(0));
 
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
   for(unsigned ii=0;ii<1/*required_nvalue(0)*/;ii++)
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
template class BellShellElement<2,2>;
template class BellShellElement<2,3>;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//=======================start_of_namespace===========================
/// Namespace for the solution of 2D linear shell equation
//====================================================================
namespace Physical_Variables
{
 /// Pressure load
 double P_ext;
 
 double epsilon = 1.0e-8;
 
 /// Exact solution as a vector
 /// differentiate u with respect to global coordinates 
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

 /// Source function applied in the normal vector 
 void source_function(const Vector<double>& x, const Vector<double>& unit_n, Vector<double>& source)
 {
  for(unsigned i=0;i<3;i++)
   {
    source[i] = epsilon*P_ext*unit_n[i];
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
const string& node_file_name,const string& element_file_name,
const string& poly_file_name) : 
 Source_fct_pt(source_fct_pt)
{ 
 // Build mesh and store pointer in Problem
 Problem::mesh_pt() = new TriangleMesh<ELEMENT>(node_file_name,element_file_name,poly_file_name);

 // set the undeformed shell
 Undef_midplane_pt = new Plate(1.0,1.0);
 
 unsigned n_node = this->mesh_pt()->nnode();
 cout << "\n\nNumber of nodes in the mesh = " << n_node << endl;
 
 // find number of elements in the mesh
 unsigned n_element = mesh_pt()->nelement();
 cout  << "Number of nelements in the mesh = " << n_element << endl;
 
 // pinning the middle nodes in each element for normal direction
 for(unsigned n=0;n<n_element;n++)
 {
  // Upcast from GeneralisedElement to the present element
  ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(n));
   
  unsigned nnode = elem_pt->nnode();
  for(unsigned i=0;i<nnode;i++)
   {
    if((i==0) || (i==1) || (i==2))
     {
     }
    else
     {
      elem_pt->node_pt(i)->pin(2);
      elem_pt->node_pt(i)->pin(3);
      elem_pt->node_pt(i)->pin(4);
      elem_pt->node_pt(i)->pin(5); 
      elem_pt->node_pt(i)->pin(6);
      elem_pt->node_pt(i)->pin(7);
     }
   }
 } // end of the middle node pinning
 
 // start_of_boundary_conditions
 // Set the boundary conditions for this problem: By default, all nodal
 // values are free -- we only need to pin the ones that have 
 // Dirichlet conditions. 
// unsigned n_side0 = mesh_pt()->nboundary_node(0);
 unsigned n_side1 = mesh_pt()->nboundary_node(1);
 //unsigned n_side2 = mesh_pt()->nboundary_node(2);
 unsigned n_side3 = mesh_pt()->nboundary_node(3);
 
 //------ PLATE BENDING PROBLEM -------------
 // Pin the single nodal value at the single node on mesh 
 // boundary 1 (= the right domain boundary at x=l)
 /// loop over the nodes on the boundary
 for(unsigned i=0;i<n_side1;i++)
  {
   // loop over the degrees of freedom that need to be
   // taken care for the Dirichlet boundary conditions
   for(unsigned j=0;j<8;j++)
    {
     mesh_pt()->boundary_node_pt(1,i)->pin(j); 
    }
  } 

 // boundary 3 (= the left domain boundary at x=0)
 /// loop over the nodes on the boundary
 for(unsigned i=0;i<n_side3;i++)
  {
   // loop over the degrees of freedom that need to be
   // taken care for the Dirichlet boundary conditions
   for(unsigned j=0;j<8;j++)
    {
     mesh_pt()->boundary_node_pt(3,i)->pin(j); 
    }
  } // end of boundary conditions

 // Loop over elements and set pointers to Physical parameters
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   //Set the source function pointer and all physical variables
   elem_pt->source_fct_pt() = Source_fct_pt;
   
   //Set the pointer to the undeformed mid-plane geometry
   elem_pt->undeformed_midplane_pt() = Undef_midplane_pt;
  } // end of loop over elements

 // Setup equation numbering scheme
 assign_eqn_numbers();

} // end of constructor

//===start_of_actions_before_newton_solve=================================
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
   unsigned n_side = mesh_pt()->nboundary_node(n);
   /// loop over the nodes on the boundary
   for(unsigned j=0;j<n_side;j++)
    {
     // Left boundary is every nodes on the left boundary 
     Node* left_node_pt=mesh_pt()->boundary_node_pt(n,j);
   
     // Loop for variables u
     ELEMENT e;
     Vector<double> u((e.required_nvalue(0)));
     for(unsigned i=0;i<(e.required_nvalue(0));i++)
      {
       // Determine the position of the boundary node (the exact solution
       // requires the coordinate in a 1D vector!)
       Vector<double> x(2);
       x[0]=left_node_pt->x(0);
       x[1]=left_node_pt->x(1);
       
       // Boundary value (read in from exact solution)
       Physical_Variables::get_exact_u(x,u);
  
       // Assign the boundary condition to nodal values
       left_node_pt->set_value(i,u[i]);
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

//===start_of_parameter_study=============================================
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
 doc_info.set_directory("RESLT_unstructured_plate");
 
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
/// Driver for 2D linearised shell problem: square flat plate
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

 // Set up the problem: 
 MyLinearisedShellProblem<BellShellElement<2,3>,2,3> //Element type as template parameter
  problem(Physical_Variables::source_function,node_file_name,element_file_name,poly_file_name);
 
 
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
  
 // Solve the problem
 problem.parameter_study();
} //end of main










