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
//Driver for a simple 2D structured mesh biharmonic problem

// Generic oomph-lib routines
#include "generic.h"

// Include the mesh
#include "meshes/one_d_mesh.h"
#include "meshes/simple_rectangular_tri_mesh.h"
#include "meshes/simple_rectangular_quadmesh.h"
#include "generic/geom_objects.h"

using namespace std;

using namespace oomph;


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//Header file for (structured mesh) Biharmonic elements
#ifndef OOMPH_BIHARMONIC_ELEMENTS_HEADER
#define OOMPH_BIHARMONIC_ELEMENTS_HEADER

namespace oomph
{

//=============================================================
/// A class for all subparametric elements that solve the 2D-
/// Biharmonic equations.
/// \f[ 
/// \frac{\partial^4 u}{\partial x_i^4} = f(x_j)
/// \f] 
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
template <unsigned DIM, unsigned NNODE_1D>
class MyBiharmonicEquations : public virtual BellElement<DIM,NNODE_1D>
{

public:
 
 /// \short Function pointer to source function fct(x,f(x)) -- 
 /// x is a Vector! 
 typedef void (*SourceFctPt)(const Vector<double>& x, double& f);


 /// \short Function pointer to gradient of source function  fct(x,g(x)) -- 
 /// x is a Vector! 
 typedef void (*SourceFctGradientPt)(const Vector<double>& x, 
                                            Vector<double>& gradient);


 /// Constructor (must initialise the Source_fct_pt to null)
 MyBiharmonicEquations() : Source_fct_pt(0), Source_fct_gradient_pt(0)
  {}
 
 /// Broken copy constructor
 MyBiharmonicEquations(const MyBiharmonicEquations& dummy) 
  { 
   BrokenCopy::broken_copy("MyBiharmonicEquations");
  } 

 /// \short Return the index at which the unknown value
 /// is stored. The default value, 0, is appropriate for single-physics
 /// problems, when there is only one variable, the value that satisfies
 /// the biharmonic equation. 
 /// In derived multi-physics elements, this function should be overloaded
 /// to reflect the chosen storage scheme. Note that these equations require
 /// that the unknown is always stored at the same index at each node.
  virtual inline unsigned u_index_biharmonic() const {return this->required_nvalue(0);}

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
    "There is no time-dependent output_fct() for Bell elements ",
    "MyBiharmonicEquations<DIM>::output_fct()",
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
    "There is no time-dependent compute_error() for Bell elements",
    "MyBiharmonicEquations<DIM>::compute_error()",
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
 

 /// Get source term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the source function might be determined by
 /// another system of equations.
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

 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution(residuals,GeneralisedElement::Dummy_matrix,0);
  }

 
 /// Add the element's contribution to its residual vector and 
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution(residuals,jacobian,1);
  }
 
 /// \short Return FE representation of function value unknown(s) 
 /// at local coordinate s
 inline Vector<double> interpolated_u_biharmonic(const Vector<double> &s) const
  {
   //Find number of position dofs
   const unsigned n_position_type = this->nnodal_position_type();
   
   //Get the index at which the unknown is stored
   const unsigned u_nodal_index = u_index_biharmonic();
   
   //Local c1-shape funtion
   Shape psi(3,n_position_type);
   DShape dpsidxi(3,n_position_type,DIM);
   DShape d2psidxi(3,n_position_type,3);
   Vector<double> interpolated_xi(DIM,0.0);
   
   //Initialise value of u
   Vector<double> interpolated_u(u_nodal_index,0.0);

   //Find values of shape function
   this->d2basis_local(s,psi,dpsidxi,d2psidxi);
   
   //Interpolated unknown
   for(unsigned l=0;l<3;l++) 
   {
    for(unsigned k=0;k<n_position_type;k++)
     {
      interpolated_u[0] += this->nodal_value(l,k)*psi(l,k);
      interpolated_u[1] += this->nodal_value(l,k)*dpsidxi(l,k,0);
      interpolated_u[2] += this->nodal_value(l,k)*dpsidxi(l,k,1);
      interpolated_u[3] += this->nodal_value(l,k)*d2psidxi(l,k,0);
      interpolated_u[4] += this->nodal_value(l,k)*d2psidxi(l,k,1);
      interpolated_u[5] += this->nodal_value(l,k)*d2psidxi(l,k,2);
     }
   } // End of varible loop
  
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
 virtual void fill_in_generic_residual_contribution(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  const unsigned& flag); 
 
 /// Pointer to source function:
 SourceFctPt Source_fct_pt;

 /// Pointer to gradient of source function
 SourceFctGradientPt Source_fct_gradient_pt;
};






///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//======================================================================
/// BiharmonicBellElement elements are a subparametric scheme with 
/// linear Lagrange interpolation for approximating the geometry and 
/// the C1-functions for approximating variables.
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
class BiharmonicBellElement : public virtual MyBiharmonicEquations<DIM,NNODE_1D>
{

private:

 /// \short Static int that holds the number of variables at 
 /// nodes: always the same
 static const unsigned Initial_Nvalue;
 
  public:


 ///\short  Constructor: Call constructors for BellElement and 
 /// Biharmonic equations
 BiharmonicBellElement() : MyBiharmonicEquations<DIM,NNODE_1D>()
  {}
 
 /// Broken copy constructor
 BiharmonicBellElement(const BiharmonicBellElement<DIM,NNODE_1D>& dummy) 
  { 
   BrokenCopy::broken_copy("BiharmonicBellElement");
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




//Inline functions:


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
 double BiharmonicBellElement<DIM,NNODE_1D>::dshape_and_dtest_eulerian_biharmonic(
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
 double BiharmonicBellElement<DIM,NNODE_1D>::d2shape_and_d2test_eulerian_biharmonic(
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
double BiharmonicBellElement<DIM,NNODE_1D>::
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
double BiharmonicBellElement<DIM,NNODE_1D>::
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
 const unsigned BiharmonicBellElement<DIM,NNODE_1D>::Initial_Nvalue = 6;

//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
void  MyBiharmonicEquations<DIM,NNODE_1D>::
fill_in_generic_residual_contribution(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              const unsigned& flag) 
{
 //Find out how many nodes there are
 const unsigned n_node = this->nnode();
 //Find out how many nodes positional dofs there are
 unsigned n_position_type = this->nnodal_position_type();
 
 //Set up memory for the shape and test functions
 Shape psi(3,n_position_type), test(n_node);
 DShape dpsidxi(3,n_position_type,DIM), dtestdxi(n_node,DIM);
 DShape d2psidxi(3,n_position_type,3), d2testdxi(n_node,3);
 
 //Set the value of n_intpt
 const unsigned n_intpt = this->integral_pt()->nweight();
 
 //Integers to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;
      
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);
   
   //Call the derivatives of the shape and test functions for normal direction
   double J = d2shape_and_d2test_eulerian_at_knot_biharmonic(ipt,psi,dpsidxi,d2psidxi,test,dtestdxi,d2testdxi);
  
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   //Calculate values of unknown
   Vector<double> interpolated_u(1,0.0);
   DenseMatrix<double> interpolated_dudxi(1,DIM,0.0);
   DenseMatrix<double> interpolated_d2udxi(1,3,0.0);
   
   //Calculate local values of unknown
   //Allocate and initialise to zero
   Vector<double> interpolated_x(DIM,0.0);
   Vector<double>  s(2);
   s[0] = this->integral_pt()->knot(ipt,0);
   s[1] = this->integral_pt()->knot(ipt,1);
   this->my_interpolated_x(s,interpolated_x);
   
   //Calculate function value and derivatives
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<3;l++) 
    {
     for(unsigned k=0;k<n_position_type;k++)
      {
       //Get the nodal value of the unknown
       double u_value = this->raw_nodal_value(l,k);
       interpolated_u[0] += u_value*psi(l,k);
       // Loop over directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_dudxi(0,j) += u_value*dpsidxi(l,k,j);
        }
       for(unsigned j=0;j<3;j++)
        {
         interpolated_d2udxi(0,j) += u_value*d2psidxi(l,k,j);
        }
      }
    }

   //Get source function
   //-------------------
   double source;
   get_source_function(ipt,interpolated_x,source);
   // Assemble residuals and Jacobian
   //--------------------------------
   for(unsigned l=0;l<3;l++)
    {
     // Loop over the test functions
     for(unsigned k=0;k<n_position_type;k++)
      {
       //Get the local equation
       local_eqn = this->nodal_local_eqn(l,k);
       //IF it's not a boundary condition
       if(local_eqn >= 0)
        {
         // Add body force/source term here 
         residuals[local_eqn] -= source*psi(l,k)*W;
         
         residuals[local_eqn] += interpolated_d2udxi(0,0)*d2psidxi(l,k,0)*W + interpolated_d2udxi(0,1)*d2psidxi(l,k,1)*W + interpolated_d2udxi(0,1)*d2psidxi(l,k,0)*W + interpolated_d2udxi(0,0)*d2psidxi(l,k,1)*W;

         //-----------------------
         if(flag)
          {
           //Loop over the velocity shape functions again
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
                 jacobian(local_eqn,local_unknown) += d2psidxi(l2,k2,0)*d2psidxi(l,k,0)*W + d2psidxi(l2,k2,1)*d2psidxi(l,k,1)*W + 2.0*d2psidxi(l2,k2,1)*d2psidxi(l,k,0)*W;
                } 
              } 
            } 
          } // End of flag
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
 
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   u = interpolated_u_biharmonic(s);

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
   u = interpolated_u_biharmonic(s);
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
   double J=this->J_eulerian(s);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   // Get x position as Vector
   this->my_interpolated_x(s,x);
   
   // Get FE function value
   Vector<double> u_fe(this->required_nvalue(0),0.0);
   u_fe = interpolated_u_biharmonic(s);
   
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
template class BiharmonicBellElement<2,2>;
}



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////






//==start_of_namespace================================================
/// Namespace for the solution of 2D Biharmonic equation
//====================================================================
namespace Physical_Variables
{
 
 /// number of elements in the x direction
 unsigned n_x = 4;
 
 /// number of elements in the y direction
 unsigned n_y = 2;

 /// length in the x direction
 double l_x = 2.0;
 
 /// length in the y direction
 double l_y = 1.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = cos(x[0])*exp(x[1]);
  u[1] = -sin(x[0])*exp(x[1]);
  u[2] = cos(x[0])*exp(x[1]);
  u[3] = -cos(x[0])*exp(x[1]);
  u[4] = cos(x[0])*exp(x[1]);
  u[5] = -sin(x[0])*exp(x[1]);
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
class MyBiharmonicBellProblem : public Problem
{

public:

 /// Constructor: Pass number of elements and pointer to source function
 MyBiharmonicBellProblem(typename MyBiharmonicEquations<DIM,NNODE_1D>::SourceFctPt source_fct_pt);

 /// Destructor (empty)
 ~MyBiharmonicBellProblem()
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
 
private:

 /// Pointer to source function
 typename MyBiharmonicEquations<DIM,NNODE_1D>::SourceFctPt Source_fct_pt;
 
}; // end of problem class



//=====start_of_constructor===============================================
/// \short Constructor for 2D Biharmonic problem.
/// Discretise the 2D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function. 
//========================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D>
MyBiharmonicBellProblem<ELEMENT,DIM,NNODE_1D>::MyBiharmonicBellProblem
(typename MyBiharmonicEquations<DIM,NNODE_1D>::SourceFctPt source_fct_pt) : 
 Source_fct_pt(source_fct_pt)
{ 
 // Build mesh and store pointer in Problem
 Problem::mesh_pt() = new SimpleRectangularTriMesh<ELEMENT>
(Physical_Variables::n_x,Physical_Variables::n_y,Physical_Variables::l_x,Physical_Variables::l_y);

 unsigned n_node = this->mesh_pt()->nnode();
 cout  << "nnode in the mesh = " << n_node << endl;
 // find number of elements in the mesh
 unsigned n_element = mesh_pt()->nelement();
 cout  << "nelement in the mesh = " << n_element << endl;
 
 // start_of_boundary_conditions
 unsigned n_bound = mesh_pt()->nboundary();
 for(unsigned i=0;i<n_bound;i++)
  {
   unsigned n_node = mesh_pt()->nboundary_node(i);
   for (unsigned n=0;n<n_node;n++)
    {
     mesh_pt()->boundary_node_pt(i,n)->pin(0); 
     mesh_pt()->boundary_node_pt(i,n)->pin(1); 
     mesh_pt()->boundary_node_pt(i,n)->pin(2);

     mesh_pt()->boundary_node_pt(i,n)->pin(3);
     mesh_pt()->boundary_node_pt(i,n)->pin(4);
     mesh_pt()->boundary_node_pt(i,n)->pin(5);
    }
  }
 // end of boundary conditionss 

 // Loop over elements and set pointers to Physical parameters
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   //Set the source function pointer and all physical variables
   elem_pt->source_fct_pt() = Source_fct_pt;
  }
 // end of pointers set up

 // Setup equation numbering scheme
 assign_eqn_numbers();

} // end of constructor




//===start_of_actions_before_newton_solve========================================
/// \short Update the problem specs before solve: (Re)set boundary values
/// from the exact solution. 
//========================================================================
template<class ELEMENT, unsigned DIM, unsigned NNODE_1D>
void MyBiharmonicBellProblem<ELEMENT,DIM,NNODE_1D>::actions_before_newton_solve()
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
void MyBiharmonicBellProblem<ELEMENT,DIM,NNODE_1D>::doc_solution(DocInfo& doc_info)
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

 // Output exact solution 
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());some_file.open(filename);
 some_file.open(filename);
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
int main()
{

 // Set up the problem: 
 // Solve a 2D Biharmonic problem 
 MyBiharmonicBellProblem<BiharmonicBellElement<2,2>,2,2> //Element type as template parameter
  problem(Physical_Variables::source_function);
 
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
 doc_info.set_directory("RESLT");
 doc_info.number()=0;

 // Solve the problem
 problem.newton_solve();
 
 //Output the solution
 problem.doc_solution(doc_info);

} // end of main









