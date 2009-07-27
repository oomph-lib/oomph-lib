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
// include required libraries
#include "math.h"
#include "generic.h"
#include "biharmonic.h"

// use oomph namespace
using namespace oomph;





//=============================================================================
// two dimensional biharmonic plate problem 1 namespace - contains all the 
// problem functions
//=============================================================================
namespace BiharmonicTestFunctions1
{

 // DIRICHLET BOUNDARY CONDITIONS

 void u_NE(const double& s, double& u)
 {
  double x = (s+1) / 2;
  u = x*x*x;
 }
 void dudn_NE(const double& s, double& dudn)
 {
  double x = (s+1) / 2;
  dudn = 3*x*x*x;
 }
 void u_SW(const double& s, double& u)
 {
  u = 0;
 }
 void dudn_SW(const double& s, double& dudn)
 {
  dudn = 0;
 }  

 // SURFACE LOAD FUNCTION

 void surface_load(const Vector<double>& x, double& f)
 {
  f = 72*x[0]*x[1];
 }


 // NEUMANN BOUNDARY CONDITIONS

 void flux1_NE(const double& s, double& flux1)
 {
  double x = (s+1) / 2;
  flux1 = 6*x*x*x + 18*x;
 }
 void flux0_NE(const double& s, double& flux0)
 {
  double x = (s+1) / 2;
  flux0 = 6*(x*x*x + x);
 }
 

 // SOLUTION

 void solution(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = x[0]*x[0]*x[0]*x[1]*x[1]*x[1];
 }
}



//=============================================================================
// Two Dimensional Biharmonic Test Problem (square)
// All Edges Clamped
//=============================================================================
class BiharmonicTestProblem1 : public BiharmonicProblem<2>
{

private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem1(const unsigned n_element)
  {
   // force linear
   Problem::Problem_is_nonlinear = false;

   // use test functions 1 namespace
   using namespace BiharmonicTestFunctions1;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(1,1);

   // assemble mesh	
   this->build_bulk_mesh(n_element, n_element, Domain_pt);

   // set the bulk element source function
   set_source_function(surface_load);

   // clamp edge on all boundaries
   set_dirichlet_boundary_condition(0,u_SW,dudn_SW);
   set_dirichlet_boundary_condition(1,u_NE,dudn_NE);
   set_dirichlet_boundary_condition(2,u_NE,dudn_NE);
   set_dirichlet_boundary_condition(3,u_SW,dudn_SW);

   // assign equation numbers
   this->build_global_mesh_and_assign_eqn_numbers();
  }

 /// Destructor - just deletes domain pt
 virtual ~BiharmonicTestProblem1()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};


//=============================================================================
// TWO DIMENSIONAL BIHARMONIC TEST 2 (PLATE PROBLEM)
//  u = (r-2)^3-r
//   theta in [0,1/4*pi]
//   r in [pi,3pi]
//=============================================================================



// two dimensional biharmonic plate problem 2 namespace - contains all the 
// problem functions
namespace BiharmonicTestFunctions2
{


 // PARAMETERS

 double Pi = MathematicalConstants::Pi;
 double theta = Pi/4;
 double r_min = sqrt(1*Pi);
 double r_max = sqrt(3*Pi);


 // BOUNDARIES
 
 void boundary_N(const double& s, Vector<double>& r)
 {  
  r[0] = (r_min + (0.5*(s+1)*(r_max-r_min)))*cos(theta);  
  r[1] = (r_min + (0.5*(s+1)*(r_max-r_min)))*sin(theta);
 }
 void boundary_E(const double& s, Vector<double>& r)
 {    
  r[0] = r_max * cos((s+1)*theta/2);
  r[1] = r_max * sin((s+1)*theta/2);
 } 
 void boundary_S(const double& s, Vector<double>& r)
 {
  r[0] = r_min + (0.5*(s+1)*(r_max-r_min));  
  r[1] = 0.0;
 }
 void boundary_W(const double& s, Vector<double>& r)
 {
  r[0] = r_min * cos((s+1)*theta/2);
  r[1] = r_min * sin((s+1)*theta/2);
 }


 // NORMALS

 void normal_N(const double& s, Vector<double>& n)
 {
  n[0] = -sin(theta);
  n[1] = cos(theta);
 }
 void normal_E(const double& s, Vector<double>& n)
 {
  double t = (s+1)*theta/2; 
  n[0] = cos(t);
  n[1] = sin(t);
 }
 void normal_S(const double& s, Vector<double>& n)
 {
  n[0] = 0.0;
  n[1] = -1.0;
 }
 void normal_W(const double& s, Vector<double>& n)
 {
  double t = (s+1)*theta/2; 
  n[0] = -cos(t);
  n[1] = -sin(t);
 }


 // DIRICHLET BCs

 void u_N(const double& s, double& u)
 {  
  double r = r_min+0.5*(s+1)*(r_max-r_min);
  u = sin(r*r)*tan(theta);
 }
 void u_E(const double& s, double& u)
 {
  double t = (s+1)*theta/2;
  u = sin(r_max*r_max)*tan(t);
 } 
 void u_S(const double& s, double& u)
 {
  u = 0;
 }
 void u_W(const double& s, double& u)
 {  
  double t = (s+1)*theta/2;
  u = sin(r_min*r_min)*tan(t);
 }
 double dudx_0(const Vector<double> x)
 {
  return(2*cos(x[0]*x[0]+x[1]*x[1])*x[1]
         -sin(x[0]*x[0]+x[1]*x[1])*x[1]/(x[0]*x[0]));
 }
 double dudx_1(const Vector<double> x)
 {
  return(2*cos(x[0]*x[0]+x[1]*x[1])*x[1]*x[1]/x[0]
         +sin(x[0]*x[0]+x[1]*x[1])/x[0]);
 }
 void dudn_N(const double& s, double& dudn)
 {  
  Vector<double> x(2);
  boundary_N(s,x);
  Vector<double> n(2);
  normal_N(s,n);
  dudn = dudx_0(x)*n[0] + dudx_1(x)*n[1];
 }
 void dudn_E(const double& s, double& dudn)
 {  
  Vector<double> x(2);
  boundary_E(s,x);
  Vector<double> n(2);
  normal_E(s,n);
  dudn = dudx_0(x)*n[0] + dudx_1(x)*n[1];
 }
 void dudn_S(const double& s, double& dudn)
 {
  double x = r_min+0.5*(s+1)*(r_max-r_min);
  dudn = -sin(x*x)/x;
 }
 void dudn_W(const double& s, double& dudn)
 {  
  Vector<double> x(2);
  boundary_W(s,x);
  Vector<double> n(2);
  normal_W(s,n);
  dudn = dudx_0(x)*n[0] + dudx_1(x)*n[1];
 } 

 // SURFACE LOAD FUNCTION
 
 void surface_load(const Vector<double>& x, double& f)
 {
  double sinr2 = sin(x[0]*x[0]+x[1]*x[1]);
  double cosr2 = cos(x[0]*x[0]+x[1]*x[1]);
  f = (16*sinr2*x[0]*x[0]*x[0]*x[1]
       -64*cosr2*x[0]*x[1]
       -48*sinr2*x[1]/x[0]
       +24*sinr2*x[1]/(x[0]*x[0]*x[0]*x[0]*x[0])
       +16*sinr2*x[1]*x[1]*x[1]*x[1]*x[1]/x[0]
       -64*cosr2*x[1]*x[1]*x[1]/x[0]
       +32*sinr2*x[1]*x[1]*x[1]*x[0]
       -16*sinr2*x[1]*x[1]*x[1]/(x[0]*x[0]*x[0]));
 }

 // SOLUTION

 void solution(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = sin(x[0]*x[0] + x[1]*x[1])*x[1]/x[0];
 }
}



//=============================================================================
// Two Dimensional Biharmonic Test Problem
// All Edges Clamped
//=============================================================================
class BiharmonicTestProblem2 : public BiharmonicProblem<2>
{

private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem2(const unsigned n_element)
  {
   // force linear
   Problem::Problem_is_nonlinear = false;

   // use test functions 1 namespace
   using namespace BiharmonicTestFunctions2;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(boundary_N,boundary_E,
                                                  boundary_S,boundary_W);

   // assemble mesh	
   this->build_bulk_mesh(n_element, n_element, Domain_pt);

   // set the bulk element source function
   set_source_function(surface_load);

   // clamp edge on all boundaries
   set_dirichlet_boundary_condition(0,u_S,dudn_S);
   set_dirichlet_boundary_condition(1,u_E,dudn_E);
   set_dirichlet_boundary_condition(2,u_N,dudn_N);
   set_dirichlet_boundary_condition(3,u_W,dudn_W);

   // assign equation numbers
   this->build_global_mesh_and_assign_eqn_numbers();
  }

 /// Destructor - just deletes domain pt
 virtual ~BiharmonicTestProblem2()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};



//=============================================================================
/// main
//=============================================================================
int main(int argc, char *argv[])
{
  // number of element
  unsigned n_element = 20;
  
  // Set up doc info
  DocInfo doc_info;
  doc_info.set_directory("RESLT");
  doc_info.number()=0;

  // Biharmonic Problem 1 (square)
  // Exact Biharmonic Preconditioner
  {
    oomph_info 
      << "/////////////////////////////////////////////////////////////////////"
      << std::endl;
    oomph_info << "TESTING: Square 2D Biharmonic Problem w/ "
	       << "Exact Preconditioning"
	       << std::endl;
    oomph_info 
      << "/////////////////////////////////////////////////////////////////////"
      << std::endl;

      // create the problem
    BiharmonicTestProblem1 problem(n_element);

    // setup the preconditioner
    BiharmonicPreconditioner* prec_pt = new BiharmonicPreconditioner;
    prec_pt->bulk_element_mesh_pt() = problem.bulk_element_mesh_pt();
    prec_pt->preconditioner_type() = 0;

    // setup the solver
    IterativeLinearSolver* solver_pt = new CG<CRDoubleMatrix>;  
    solver_pt->preconditioner_pt() = prec_pt;
 
    // apply the solver
    problem.linear_solver_pt() = solver_pt;
    problem.newton_solve();

    // ouput the solution
    problem.doc_solution(doc_info,BiharmonicTestFunctions1::solution);
    doc_info.number()++;

    // clean up
    delete solver_pt;
    delete prec_pt;
  }

  // Biharmonic Problem 2 (section of annulus)
  // Inexact Biharmonic Preconditioner w/ SuperLU
  {
    oomph_info 
      << "/////////////////////////////////////////////////////////////////////"
      << std::endl;
    oomph_info << "TESTING: Curved 2D Biharmonic Problem w/ "
	       << "Inexact Preconditioning"
	       << std::endl;
    oomph_info 
      << "/////////////////////////////////////////////////////////////////////"
      << std::endl;

    // create the problem
    BiharmonicTestProblem2 problem(n_element);

    // setup the preconditioner
    BiharmonicPreconditioner* prec_pt = new BiharmonicPreconditioner;
    prec_pt->bulk_element_mesh_pt() = problem.bulk_element_mesh_pt();
    prec_pt->preconditioner_type() = 1;

    // setup the solver
    IterativeLinearSolver* solver_pt = new CG<CRDoubleMatrix>;  
    solver_pt->preconditioner_pt() = prec_pt;
 
    // apply the solver
    problem.linear_solver_pt() = solver_pt;
    problem.newton_solve();

    // ouput the solution
    problem.doc_solution(doc_info,BiharmonicTestFunctions2::solution);
    doc_info.number()++;

    // clean up
    delete solver_pt;
    delete prec_pt;
  }
}
